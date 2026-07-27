"""VALIDATION RUN #1 — the model's drug-partition prediction vs MEASURED human
glioma tissue:plasma exposure, region by region.

This is the first time tmesim's molecular/transport layer is scored against real
*measured drug concentrations in human brain tumour tissue*. It is a delivery /
PARTITION test only: it asks whether the fraction of plasma drug that shows up
in glioma tissue is right. It says nothing about response, kill, or potency.

Nothing here modifies core / physics / molecular / reconstruct / engine. The
module only *reads* them (``molecular.descriptors``, ``molecular.permeability``,
``molecular.library``) and *reads* two on-disk data tables.

THE TWO REAL DATA SOURCES
-------------------------
1. ``data_cache/hunt/autoradiography_microdialysis/parsed/
   regional_intratumoral_PK_anchors_literature.csv``
   27 rows hand-extracted from six peer-reviewed papers, every row carrying a
   PMID and a DOI. Human rows report ``tissue_to_plasma_auc_pct`` — a percentage
   AUC ratio between an intratumoral/peritumoral compartment and plasma — plus
   the region the probe or specimen came from. Rat rows (Grossman 2012/2013)
   report tissue AUC only, with no paired plasma, and are therefore EXCLUDED
   from the scorecard (they can only support a within-animal ordering).

2. ``data_cache/hunt/autoradiography_microdialysis/parsed/
   oh2026_neurooncol_TableS3_per_sample_tumor_drug_ratios.csv``
   Oh et al., Neuro-Oncology 2026 (``noaf204``, PMC12962645), Table S3:
   **125 surgical samples from 42 high-grade-glioma cases**. Each sample was
   taken from an MRI-defined region during resection and assayed for two probe
   drugs given intra-operatively:
       ``lev`` = **levetiracetam** (small, neutral, BBB-penetrant antiepileptic)
       ``cef`` = **cefazolin**     (anionic cephalosporin, BBB-excluded)
   Columns ``lev_tumor_to_plasma_ratio`` / ``cef_tumor_to_plasma_ratio`` are
   single-timepoint tumour-tissue-concentration : plasma-concentration ratios
   (dimensionless), not AUC ratios. ``region_contrast_enhancing`` (C) and
   ``region_flair`` (F) are ``+``/``-`` MRI flags; ``region_necrosis`` marks
   samples the pathologist called necrotic (those rows have C/F = N/A).
   ``mannitol_status`` is either "Not administered", "Unknown", or the number of
   minutes between mannitol dosing and sampling — mannitol is an osmotic agent
   that can itself open the barrier, so it is carried as a confounder.

HOW A ``bbb_penetrance`` BECOMES A PREDICTED TISSUE:PLASMA RATIO  (ASSUMPTION A1)
--------------------------------------------------------------------------------
The frozen engine's delivery term is (``sim/engine.py`` line 45)::

    source = vascular_density * bbb_penetrance * plasma_mean      [uM/s]

and ``core/config.py`` defines the field as "[0,1] fraction of vascular
concentration that crosses into" the tissue. So under the engine's own
semantics, with ``vascular_density`` normalised to 1 and the tissue held at
quasi-equilibrium with a sustained plasma level (the Stage-0 plasma-averaging
approximation the engine already makes), the tissue:plasma exposure ratio the
model implies is the identity::

    R_predicted  =  bbb_penetrance                                  (A1)

That is the mapping used throughout this module. Three consequences follow
directly and are reported rather than hidden:

  * **R_predicted is bounded by 1.** Any measured ratio above 1 is structurally
    unreachable by the model, no matter how it is parameterised.
  * **R_predicted has no region term.** ``bbb_penetrance`` is a function of the
    drug's 2D descriptors alone; it is identical in contrast-enhancing,
    non-enhancing and necrotic tissue. The model's predicted regional fold
    change is therefore *exactly 1.0*.
  * **R_predicted is a passive free-drug crossing fraction.** Microdialysis
    measures free extracellular drug and is the matched comparator. Tissue
    *homogenate* assays measure total (free + intracellular + lipid-bound) drug
    and can legitimately exceed 1. Comparing A1 to a homogenate ratio is a
    known apples-to-oranges step; every row in the scorecard is labelled with
    its measurement modality so this is visible.

Units: partition ratios are dimensionless fractions throughout (the literature
CSV's percentages are divided by 100 on load).
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from ..molecular.descriptors import from_smiles
from ..molecular.library import GBM_DRUGS
from ..molecular.permeability import (
    _radius_nm_from_descriptors,
    bbb_penetrance,
    membrane_permeability,
)
from ..reconstruct.barrier import (
    L_PARACELLULAR_LEAK,
    REGION_DISRUPTION,
    effective_penetrance,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
_PARSED = os.path.join(
    REPO_ROOT, "data_cache", "hunt", "autoradiography_microdialysis", "parsed"
)
LIT_CSV = os.path.join(_PARSED, "regional_intratumoral_PK_anchors_literature.csv")
PER_SAMPLE_CSV = os.path.join(
    _PARSED, "oh2026_neurooncol_TableS3_per_sample_tumor_drug_ratios.csv"
)

RESULTS_DIR = os.path.join(REPO_ROOT, "results")
REPORT_MD = os.path.join(RESULTS_DIR, "glioma_regional_pk.md")
REPORT_PNG = os.path.join(RESULTS_DIR, "glioma_regional_pk.png")

# ---------------------------------------------------------------------------
# Drugs we can map onto the molecular layer.
#
# temozolomide and methotrexate come straight out of `molecular.library`.
# The other four are measured in the data but absent from the library; rather
# than modify the library (frozen for this run) we carry their VERIFIED SMILES
# here. Each was pulled from PubChem PUG-REST by name and is checked at import
# time against the literature molecular formula + MW (`_verify_smiles`), so a
# silently-wrong structure cannot reach the scorecard.
# ---------------------------------------------------------------------------
_EXTRA_SMILES: Dict[str, Tuple[str, str, float]] = {
    # name: (SMILES, molecular formula, literature MW g/mol)
    "levetiracetam": (
        "CC[C@@H](C(=O)N)N1CCCC1=O",
        "C8H14N2O2",
        170.21,
    ),
    "cefazolin": (
        "CC1=NN=C(S1)SCC2=C(N3[C@@H]([C@@H](C3=O)NC(=O)CN4C=NN=N4)SC2)C(=O)O",
        "C14H14N8O4S3",
        454.51,
    ),
    "adavosertib": (
        "CC(C)(C1=NC(=CC=C1)N2C3=NC(=NC=C3C(=O)N2CC=C)NC4=CC=C(C=C4)N5CCN(CC5)C)O",
        "C27H32N8O2",
        500.60,
    ),
    # Not a glioma drug; carried only so the permeability-ceiling before/after
    # table in the post-fix section can be computed rather than quoted. Imatinib
    # is the molecule that surfaced the defect (results/gist_imatinib.md).
    "imatinib": (
        "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
        "C29H31N7O",
        493.60,
    ),
    "eribulin": (
        "C[C@@H]1C[C@@H]2CC[C@H]3C(=C)C[C@@H](O3)CC[C@]45C[C@@H]6[C@H](O4)[C@H]7"
        "[C@@H](O6)[C@@H](O5)[C@@H]8[C@@H](O7)CC[C@@H](O8)CC(=O)C[C@H]9[C@H]"
        "(C[C@H](C1=C)O2)O[C@@H]([C@@H]9OC)C[C@@H](CN)O",
        "C40H59NO11",
        729.91,
    ),
}

# Short pharmacological note per drug, for the report.
_DRUG_NOTE: Dict[str, str] = {
    "temozolomide": "first-line GBM alkylator; small, neutral, low logP",
    "methotrexate": "antifolate; large, very polar, 5 HBD, actively effluxed",
    "levetiracetam": "antiepileptic BBB-penetrant probe (Oh 2026 probe 1)",
    "cefazolin": "anionic cephalosporin, BBB-excluded probe (Oh 2026 probe 2)",
    "adavosertib": "WEE1 inhibitor (ABTC trial); lipophilic, P-gp substrate",
    "eribulin": "microtubule inhibitor; large macrocycle, strong P-gp substrate",
}

# Order-of-magnitude bands used for pass/fail.
BAND_TIGHT = 3.0
BAND_OOM = 10.0


# ---------------------------------------------------------------------------
# Molecular layer
# ---------------------------------------------------------------------------
def _verify_smiles(name: str) -> None:
    """Assert an `_EXTRA_SMILES` entry reproduces its literature MW/formula."""
    smiles, formula, mw_lit = _EXTRA_SMILES[name]
    mol = from_smiles(name, smiles)
    if mol.source != "rdkit":
        raise ValueError(f"{name}: SMILES did not parse in RDKit (source={mol.source})")
    mw = float(mol.get("mw"))
    if abs(mw - mw_lit) > 0.5:
        raise ValueError(f"{name}: RDKit MW {mw:.2f} != literature {mw_lit} ({formula})")


def _bbb_penetrance_arithmetic(descriptors: Dict[str, float]) -> float:
    """FROZEN reference implementation of the **pre-fix** ``bbb_penetrance``.

    Byte-for-byte the aggregation that shipped before validation run #1's fixes:
    the four CNS-MPO desirabilities combined by an ARITHMETIC mean and squashed
    through the same logistic. It is kept here — in the validation module, not in
    the molecular layer — for exactly one purpose: so the "before" column of the
    before/after scorecard is **computed live from the same descriptors**, rather
    than pasted in from a previous run where it could silently rot.

    Do NOT call this from model code. It is the defect, retained as evidence:
    an arithmetic mean cannot express a veto, so a molecule with two
    disqualifying properties is rescued by two favourable ones (eribulin, 66x).

    The per-descriptor desirability curves are IDENTICAL to the shipped
    :func:`tmesim.molecular.permeability.bbb_penetrance`; only the aggregation
    differs, which is the whole point of the comparison.
    """
    if descriptors.get("is_antibody"):
        return 0.0
    mw = float(descriptors.get("mw", 300.0))
    logp = float(descriptors.get("logp", 0.0))
    tpsa = float(descriptors.get("tpsa", 75.0))
    hbd = float(descriptors.get("hbd", 1.0))
    if mw >= 1000.0:
        return 0.0
    d_mw = 1.0 / (1.0 + math.exp((mw - 450.0) / 60.0))
    d_logp = math.exp(-((logp - 2.5) ** 2) / (2.0 * 1.7 ** 2))
    d_tpsa = 1.0 / (1.0 + math.exp((tpsa - 90.0) / 12.0))
    d_hbd = 1.0 / (1.0 + math.exp((hbd - 3.0) / 0.8))
    mean_desir = (d_mw + d_logp + d_tpsa + d_hbd) / 4.0
    score = 1.0 / (1.0 + math.exp(-8.0 * (mean_desir - 0.5)))
    return max(0.0, min(1.0, score))


#: ``aggregation`` keyword -> penetrance function. ``"arithmetic"`` is the
#: pre-fix defect; ``"geometric"`` is what the molecular layer ships now.
_AGGREGATIONS = {
    "arithmetic": _bbb_penetrance_arithmetic,
    "geometric": bbb_penetrance,
}


def _membrane_permeability_unbounded(descriptors: Dict[str, float],
                                     temperature_K: float = 310.0) -> float:
    """FROZEN reference implementation of the **pre-fix** ``membrane_permeability``.

    ``P = Kp * Dm / d`` with logP clamped only at 5 and a membrane viscosity of
    0.1 Pa*s, and NO ceiling of any kind. Returned 4.23e+07 um/s (42 m/s) for
    imatinib. Kept here purely so the before/after table is computed rather than
    quoted. Do NOT call this from model code.
    """
    if descriptors.get("is_antibody"):
        return 0.0
    logp = max(-3.0, min(5.0, float(descriptors.get("logp", 0.0))))
    kp = 10.0 ** logp
    r_m = max(_radius_nm_from_descriptors(descriptors), 0.1) * 1e-9
    dm = (1.380649e-23 * temperature_K) / (6.0 * math.pi * 0.10 * r_m)
    return max(0.0, (kp * dm / (4.0 * 1e-9)) * 1e6)


#: Drugs to show in the permeability-ceiling before/after table, with the SMILES
#: needed for any that are outside ``molecular.library``.
_PERM_CHECK = ["temozolomide", "cisplatin", "methotrexate", "doxorubicin",
               "lomustine", "imatinib", "bevacizumab"]


def molecular_profile(drug: str, aggregation: str = "geometric") -> Dict:
    """Descriptors + the two structure-only transport numbers for one drug.

    Returns mw / logp / tpsa / hbd, ``bbb_penetrance`` (the predicted partition
    under assumption A1) and passive ``membrane_permeability`` (um/s).

    ``aggregation`` selects how the four CNS-MPO desirabilities are combined:
    ``"geometric"`` (the shipped :func:`~tmesim.molecular.permeability.bbb_penetrance`,
    a veto) or ``"arithmetic"`` (:func:`_bbb_penetrance_arithmetic`, the pre-fix
    defect, retained so before/after is computed rather than quoted).
    """
    key = drug.strip().lower()
    if key in GBM_DRUGS:
        mol = from_smiles(key, GBM_DRUGS[key])
        origin = "molecular.library"
    elif key in _EXTRA_SMILES:
        _verify_smiles(key)
        mol = from_smiles(key, _EXTRA_SMILES[key][0])
        origin = "PubChem SMILES (MW-verified in this module)"
    else:
        raise KeyError(f"drug {drug!r} not mappable to the molecular layer")
    d = mol.descriptors
    try:
        penetrance_fn = _AGGREGATIONS[aggregation]
    except KeyError:
        raise ValueError(
            f"aggregation must be one of {sorted(_AGGREGATIONS)}; got {aggregation!r}"
        ) from None
    return {
        "drug": key,
        "smiles_origin": origin,
        "descriptor_source": mol.source,
        "aggregation": aggregation,
        "mw": float(d["mw"]),
        "logp": float(d["logp"]),
        "tpsa": float(d["tpsa"]),
        "hbd": int(d["hbd"]),
        "predicted_partition": float(penetrance_fn(d)),  # ASSUMPTION A1
        "membrane_permeability_um_s": float(membrane_permeability(d)),
        "note": _DRUG_NOTE.get(key, ""),
    }


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_literature_anchors() -> pd.DataFrame:
    """Load the 27-row PMID-anchored literature table.

    Adds ``tissue_to_plasma_ratio`` (the ``_pct`` column / 100) and, for the
    Eroglu eribulin row whose value is a *range string* ("0.13-1.99"), the
    geometric midpoint of the range. Rows with no ratio at all (the rat
    Grossman rows) keep NaN.
    """
    if not os.path.exists(LIT_CSV):
        raise FileNotFoundError(f"literature anchor CSV missing: {LIT_CSV}")
    df = pd.read_csv(LIT_CSV, dtype={"tissue_to_plasma_auc_pct": str})

    def _pct_to_ratio(v) -> float:
        if v is None or (isinstance(v, float) and math.isnan(v)):
            return float("nan")
        s = str(v).strip()
        if not s or s.lower() == "nan":
            return float("nan")
        if "-" in s and not s.startswith("-"):
            lo, hi = (float(x) for x in s.split("-", 1))
            return math.sqrt(max(lo, 1e-12) * max(hi, 1e-12)) / 100.0  # geometric mid
        return float(s) / 100.0

    df["tissue_to_plasma_ratio"] = df["tissue_to_plasma_auc_pct"].map(_pct_to_ratio)
    return df


def load_per_sample_ratios() -> pd.DataFrame:
    """Load Oh 2026 Table S3 and attach a single tidy ``region`` label.

    Region resolution (mutually exclusive, covers all 125 rows):
        ``necrosis``          -- ``region_necrosis == 'Necrosis'`` (C/F are N/A)
        ``contrast_enhancing``-- ``region_contrast_enhancing == '+'``
        ``non_enhancing``     -- ``region_contrast_enhancing == '-'``
    ``non_enhancing`` is further split by the FLAIR flag into
    ``non_enhancing_flair`` / ``non_enhancing_normal`` in ``region_fine``.

    Also adds ``mannitol_exposed`` (True where a minutes-since-mannitol number
    was recorded) and ``n_censored`` bookkeeping is left to the caller: a ratio
    of exactly 0.00 is a below-limit-of-quantification report, not a true zero.
    """
    if not os.path.exists(PER_SAMPLE_CSV):
        raise FileNotFoundError(f"per-sample CSV missing: {PER_SAMPLE_CSV}")
    df = pd.read_csv(PER_SAMPLE_CSV)
    df["case"] = df["case"].astype(str).str.replace(" ", "", regex=False)

    def _region(row) -> str:
        if str(row["region_necrosis"]).strip().lower().startswith("necro"):
            return "necrosis"
        ce = str(row["region_contrast_enhancing"]).strip()
        if ce == "+":
            return "contrast_enhancing"
        if ce == "-":
            return "non_enhancing"
        return "unclassified"

    def _region_fine(row) -> str:
        base = _region(row)
        if base != "non_enhancing":
            return base
        return (
            "non_enhancing_flair"
            if str(row["region_flair"]).strip() == "+"
            else "non_enhancing_normal"
        )

    df["region"] = df.apply(_region, axis=1)
    df["region_fine"] = df.apply(_region_fine, axis=1)

    def _mannitol(v) -> bool:
        s = str(v).strip()
        try:
            float(s)
            return True
        except ValueError:
            return False

    df["mannitol_exposed"] = df["mannitol_status"].map(_mannitol)
    return df


# ---------------------------------------------------------------------------
# Measured summaries
# ---------------------------------------------------------------------------
def _summarise(vals: Sequence[float]) -> Dict:
    a = np.asarray([v for v in vals if v is not None and not math.isnan(v)], float)
    if a.size == 0:
        return {"n": 0, "mean": float("nan"), "median": float("nan"),
                "min": float("nan"), "max": float("nan")}
    return {
        "n": int(a.size),
        "mean": float(a.mean()),
        "median": float(np.median(a)),
        "min": float(a.min()),
        "max": float(a.max()),
    }


def measured_partitions() -> List[Dict]:
    """Per drug x region measured tissue:plasma partition, with provenance.

    Human data only. Each returned row carries n, mean, median, range, method
    (microdialysis = free ECF vs homogenate = total drug), PMID and citation.
    """
    lit = load_literature_anchors()
    human = lit[(lit["species"] == "human") & lit["tissue_to_plasma_ratio"].notna()]

    rows: List[Dict] = []
    for (drug, region, method, pmid, study), g in human.groupby(
        ["drug", "region", "method", "source_pmid", "study"], dropna=False
    ):
        s = _summarise(g["tissue_to_plasma_ratio"].tolist())
        s.update(
            drug=str(drug), region=str(region), method=str(method),
            pmid=str(pmid), citation=str(study),
            unit="tissue:plasma AUC ratio",
            source="literature_anchors_csv",
        )
        rows.append(s)

    ps = load_per_sample_ratios()
    for drug, col in (("levetiracetam", "lev_tumor_to_plasma_ratio"),
                      ("cefazolin", "cef_tumor_to_plasma_ratio")):
        for region, g in ps.groupby("region"):
            s = _summarise(g[col].tolist())
            if s["n"] == 0:
                continue
            s.update(
                drug=drug, region=region, method="tissue_homogenate_LCMS",
                pmid="PMC12962645", citation="Oh 2026 Neuro-Oncology (noaf204) Table S3",
                unit="tissue:plasma concentration ratio (single timepoint)",
                source="oh2026_per_sample_csv",
            )
            rows.append(s)
    return rows


def _drug_level_measured() -> List[Dict]:
    """One pooled measured partition per drug (all human regions, equal weight
    per reported subject/cohort value). This is what the headline scorecard uses.
    """
    lit = load_literature_anchors()
    human = lit[(lit["species"] == "human") & lit["tissue_to_plasma_ratio"].notna()]
    out: List[Dict] = []
    for drug, g in human.groupby("drug"):
        s = _summarise(g["tissue_to_plasma_ratio"].tolist())
        s.update(
            drug=str(drug),
            method="|".join(sorted(set(g["method"]))),
            pmids="|".join(sorted(set(g["source_pmid"].astype(str)))),
            citations="; ".join(sorted(set(g["study"]))),
            regions="|".join(sorted(set(g["region"]))),
            source="literature_anchors_csv",
        )
        out.append(s)

    ps = load_per_sample_ratios()
    for drug, col in (("levetiracetam", "lev_tumor_to_plasma_ratio"),
                      ("cefazolin", "cef_tumor_to_plasma_ratio")):
        s = _summarise(ps[col].tolist())
        s.update(
            drug=drug, method="tissue_homogenate_LCMS", pmids="PMC12962645",
            citations="Oh 2026 Neuro-Oncology (noaf204) Table S3",
            regions="contrast_enhancing|necrosis|non_enhancing",
            source="oh2026_per_sample_csv",
        )
        out.append(s)
    return out


# ---------------------------------------------------------------------------
# TEST 1 -- the honest scorecard
# ---------------------------------------------------------------------------
#: MRI/pathology region labels that correspond to a DISRUPTED barrier. Necrosis
#: is included because it sits inside the enhancing rim and measures as leaky
#: (see ``reconstruct.barrier.REGION_DISRUPTION``); ``non_enhancing`` /
#: ``nonenhancing`` / ``noncontrast_enhancing`` / ``peritumoral_brain`` are the
#: intact-barrier reference.
_DISRUPTED_REGIONS = {"contrast_enhancing", "necrosis"}


def measured_sample_disruption() -> Dict[str, float]:
    """Mean ``bbb_disruption`` implied by each drug's MEASURED sample composition.

    Each drug's scorecard "measured mean" POOLS whatever regions that study
    sampled. Levetiracetam's 0.697, for instance, is the mean over 39
    contrast-enhancing + 17 necrotic + 69 non-enhancing specimens. A model with a
    regional barrier term must therefore be evaluated at the SAME regional mix,
    or the comparison silently penalises it for a difference in sampling rather
    than in physics.

    Returns drug -> fraction of that drug's measured samples drawn from
    disrupted-barrier regions, in [0, 1]. This is a composition weight, NOT a
    fitted quantity: it is read straight off the region column of the source
    tables.
    """
    out: Dict[str, float] = {}
    ps = load_per_sample_ratios()
    for drug, col in _DRUG_COL.items():
        sub = ps[ps[col].notna()]
        out[drug] = float(sub["region"].isin(_DISRUPTED_REGIONS).mean()) if len(sub) else 0.0
    lit = load_literature_anchors()
    lit = lit[(lit["species"] == "human") & lit["tissue_to_plasma_ratio"].notna()]
    for drug, g in lit.groupby("drug"):
        out[str(drug)] = float(g["region"].isin(_DISRUPTED_REGIONS).mean())
    return out


def run_scorecard(aggregation: str = "geometric", apply_leak: bool = False) -> Dict:
    """Predicted (A1) vs measured tissue:plasma partition, per drug.

    Reports the fold error ``predicted / measured_mean`` and pass/fail against a
    3x band and a 10x (order-of-magnitude) band, plus the aggregate systematic
    bias as the geometric mean of the fold errors.

    Parameters
    ----------
    aggregation : {"geometric", "arithmetic"}
        Desirability aggregation inside ``bbb_penetrance``. ``"arithmetic"`` is
        the pre-fix defect, retained for the before/after comparison.
    apply_leak : bool
        If True, add the fitted additive paracellular leak
        (``reconstruct.barrier.L_PARACELLULAR_LEAK``) weighted by each drug's
        measured regional composition (:func:`measured_sample_disruption`), i.e.
        score the model AS IT NOW EXISTS rather than its structure-only core.

    Two aggregate bias figures are reported and they answer different questions:

    * ``geometric_mean_fold_error`` — the geometric mean of SIGNED fold errors.
      Over- and under-predictions CANCEL, so a value near 1.0 can be a genuinely
      unbiased model or two large opposite errors annihilating.
    * ``geometric_mean_abs_fold_error`` — the geometric mean of
      ``max(fold, 1/fold)``. Nothing cancels; this is the honest accuracy figure
      and the one to read when the signed figure looks suspiciously good.
    """
    leak_weights = measured_sample_disruption() if apply_leak else {}
    rows: List[Dict] = []
    for m in _drug_level_measured():
        drug = m["drug"]
        try:
            prof = molecular_profile(drug, aggregation=aggregation)
        except KeyError:
            continue
        pred = prof["predicted_partition"]
        prof["structure_only_partition"] = pred
        disruption = float(leak_weights.get(drug, 0.0)) if apply_leak else 0.0
        prof["sample_mean_bbb_disruption"] = disruption
        if apply_leak:
            pred = float(effective_penetrance(pred, disruption))
            prof["predicted_partition"] = pred
        meas = m["mean"]
        fold = pred / meas if meas > 0 else float("inf")
        rows.append({
            **prof,
            "measured_mean": meas,
            "measured_median": m["median"],
            "measured_min": m["min"],
            "measured_max": m["max"],
            "n_measured": m["n"],
            "method": m["method"],
            "regions": m["regions"],
            "pmids": m["pmids"],
            "citations": m["citations"],
            "fold_error": fold,
            "log10_fold": math.log10(fold) if fold > 0 and math.isfinite(fold) else float("nan"),
            "pass_3x": (fold <= BAND_TIGHT) and (fold >= 1.0 / BAND_TIGHT),
            "pass_10x": (fold <= BAND_OOM) and (fold >= 1.0 / BAND_OOM),
            "free_drug_comparator": "microdialysis" in str(m["method"]),
        })
    rows.sort(key=lambda r: r["measured_mean"])

    logs = [r["log10_fold"] for r in rows if math.isfinite(r["log10_fold"])]
    gm = 10.0 ** (sum(logs) / len(logs)) if logs else float("nan")
    # bias excluding the single worst outlier, to show how much it drives things
    worst = max(rows, key=lambda r: abs(r["log10_fold"])) if rows else None
    logs_x = [r["log10_fold"] for r in rows
              if r is not worst and math.isfinite(r["log10_fold"])]
    gm_x = 10.0 ** (sum(logs_x) / len(logs_x)) if logs_x else float("nan")

    # free-drug (microdialysis) subset only -- the modality-matched comparison
    md = [r for r in rows if r["free_drug_comparator"] and math.isfinite(r["log10_fold"])]
    gm_md = 10.0 ** (sum(r["log10_fold"] for r in md) / len(md)) if md else float("nan")
    md_x = [r for r in md if r is not worst]
    gm_md_x = (10.0 ** (sum(r["log10_fold"] for r in md_x) / len(md_x))
               if md_x else float("nan"))

    # UNSIGNED accuracy: |log10 fold|, so over- and under-prediction cannot
    # cancel each other out into a flattering aggregate.
    alogs = [abs(r["log10_fold"]) for r in rows if math.isfinite(r["log10_fold"])]
    gma = 10.0 ** (sum(alogs) / len(alogs)) if alogs else float("nan")
    alogs_x = [abs(r["log10_fold"]) for r in rows
               if r is not worst and math.isfinite(r["log10_fold"])]
    gma_x = 10.0 ** (sum(alogs_x) / len(alogs_x)) if alogs_x else float("nan")
    worst_abs = max((abs(r["fold_error"]) if r["fold_error"] >= 1.0
                     else 1.0 / r["fold_error"]) for r in rows) if rows else float("nan")

    # Explicitly ex-ERIBULIN, which is NOT always the same as ex-worst: eribulin
    # is the worst drug pre-fix, but under the geometric mean alone methotrexate
    # overtakes it. The narrative is about eribulin, so it needs its own figure
    # rather than silently switching which drug got dropped.
    def _ex(drug: str, absolute: bool):
        v = [(abs(r["log10_fold"]) if absolute else r["log10_fold"]) for r in rows
             if r["drug"] != drug and math.isfinite(r["log10_fold"])]
        return 10.0 ** (sum(v) / len(v)) if v else float("nan")

    return {
        "rows": rows,
        "aggregation": aggregation,
        "leak_applied": bool(apply_leak),
        "leak_constant": L_PARACELLULAR_LEAK if apply_leak else 0.0,
        "n_drugs": len(rows),
        "n_pass_3x": sum(r["pass_3x"] for r in rows),
        "n_pass_10x": sum(r["pass_10x"] for r in rows),
        "geometric_mean_fold_error": gm,
        "geometric_mean_fold_error_ex_worst": gm_x,
        "geometric_mean_abs_fold_error": gma,
        "geometric_mean_abs_fold_error_ex_worst": gma_x,
        "geometric_mean_fold_error_ex_eribulin": _ex("eribulin", False),
        "geometric_mean_abs_fold_error_ex_eribulin": _ex("eribulin", True),
        "worst_abs_fold_error": worst_abs,
        "worst_drug": worst["drug"] if worst else None,
        "worst_fold": worst["fold_error"] if worst else float("nan"),
        "geometric_mean_fold_error_microdialysis_only": gm_md,
        "geometric_mean_fold_error_microdialysis_ex_worst": gm_md_x,
        "n_microdialysis": len(md),
    }


# ---------------------------------------------------------------------------
# TEST 2 -- the regional test (125 samples / 42 cases)
# ---------------------------------------------------------------------------
def _cliffs_delta(a: Sequence[float], b: Sequence[float]) -> float:
    """Cliff's delta in [-1, 1]: P(a>b) - P(a<b). Nonparametric, tolerates the
    many below-LLOQ zeros in the cefazolin column."""
    x, y = np.asarray(a, float), np.asarray(b, float)
    if x.size == 0 or y.size == 0:
        return float("nan")
    gt = (x[:, None] > y[None, :]).sum()
    lt = (x[:, None] < y[None, :]).sum()
    return float((gt - lt) / (x.size * y.size))


def _hedges_g(a: Sequence[float], b: Sequence[float]) -> float:
    x, y = np.asarray(a, float), np.asarray(b, float)
    if x.size < 2 or y.size < 2:
        return float("nan")
    n1, n2 = x.size, y.size
    sp = math.sqrt(((n1 - 1) * x.var(ddof=1) + (n2 - 1) * y.var(ddof=1)) / (n1 + n2 - 2))
    if sp == 0:
        return float("nan")
    d = (x.mean() - y.mean()) / sp
    j = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)  # small-sample correction
    return float(d * j)


def _perm_p_naive(a: Sequence[float], b: Sequence[float], n_perm: int = 20000,
                  seed: int = 0) -> float:
    """Two-sided permutation p on the difference of means, shuffling labels
    freely across ALL samples. This IGNORES the fact that samples are clustered
    within patients and is therefore ANTI-CONSERVATIVE; reported only alongside
    the cluster-aware test below."""
    x, y = np.asarray(a, float), np.asarray(b, float)
    if x.size == 0 or y.size == 0:
        return float("nan")
    obs = abs(x.mean() - y.mean())
    pool = np.concatenate([x, y])
    n1, n2 = x.size, y.size
    total = pool.sum()
    rng = np.random.default_rng(seed)
    # Vectorised: draw n_perm independent shuffles at once via argsort of noise.
    idx = np.argsort(rng.random((n_perm, pool.size)), axis=1)[:, :n1]
    sum_a = pool[idx].sum(axis=1)
    diff = sum_a / n1 - (total - sum_a) / n2
    hits = int((np.abs(diff) >= obs - 1e-12).sum())
    return (hits + 1) / (n_perm + 1)


def _perm_p_within_case(df: pd.DataFrame, value_col: str, region_col: str,
                        ra: str, rb: str, n_perm: int = 20000,
                        seed: int = 0) -> Dict:
    """Cluster-aware permutation: shuffle the REGION LABEL within each patient.

    Multiple samples come from the same tumour, so the free permutation above
    over-counts independent observations. Permuting labels only *within* a case
    preserves the patient as the unit of clustering and tests the paired,
    within-tumour regional contrast — the biologically correct null. Cases
    contributing only one of the two labels carry no information and drop out
    automatically (their permuted labels are unchanged).
    """
    sub = df[df[region_col].isin([ra, rb])][["case", region_col, value_col]].dropna()
    if sub.empty:
        return {"p": float("nan"), "n_informative_cases": 0, "observed_diff": float("nan")}
    groups = [g for _, g in sub.groupby("case")]
    informative = [g for g in groups if g[region_col].nunique() > 1]

    # Under within-case permutation the per-case COUNT of each label is fixed, so
    # the pooled group sizes n_a / n_b never change and only the pooled sum of the
    # A-labelled values varies. That lets the whole null be built vectorised:
    # per case, draw n_perm random k_g-subsets of its values and sum them.
    rng = np.random.default_rng(seed)
    n_a = n_b = 0
    obs_sum_a = 0.0
    total = 0.0
    sum_a_null = np.zeros(n_perm)
    for g in groups:
        v = g[value_col].to_numpy(dtype=float)
        is_a = (g[region_col].to_numpy() == ra)
        k = int(is_a.sum())
        n_a += k
        n_b += int(v.size - k)
        obs_sum_a += float(v[is_a].sum())
        total += float(v.sum())
        if k == 0:
            continue
        if k == v.size:
            sum_a_null += v.sum()
            continue
        idx = np.argsort(rng.random((n_perm, v.size)), axis=1)[:, :k]
        sum_a_null += v[idx].sum(axis=1)
    if n_a == 0 or n_b == 0:
        return {"p": float("nan"), "n_informative_cases": len(informative),
                "observed_diff": float("nan")}
    obs = obs_sum_a / n_a - (total - obs_sum_a) / n_b
    null = sum_a_null / n_a - (total - sum_a_null) / n_b
    hits = int((np.abs(null) >= abs(obs) - 1e-12).sum())
    return {
        "p": (hits + 1) / (n_perm + 1),
        "n_informative_cases": len(informative),
        "observed_diff": float(obs),
    }


_REGION_PAIRS = [
    ("contrast_enhancing", "non_enhancing"),
    ("contrast_enhancing", "necrosis"),
    ("necrosis", "non_enhancing"),
]

_DRUG_COL = {
    "levetiracetam": "lev_tumor_to_plasma_ratio",
    "cefazolin": "cef_tumor_to_plasma_ratio",
}


def run_regional_test(n_perm: int = 20000, seed: int = 0) -> Dict:
    """Do MEASURED tumour:plasma ratios differ by MRI region, and does the model
    reproduce it?

    For each probe drug and each region pair, reports group means/medians, the
    ratio of means, Cliff's delta, Hedges' g, a naive sample-level permutation p
    and a cluster-aware within-case permutation p. The model side is trivial and
    is reported as such: assumption A1 gives one drug-level number with no
    region dependence, so its predicted regional fold change is exactly 1.0.
    """
    ps = load_per_sample_ratios()
    out: Dict = {
        "n_samples": int(len(ps)),
        "n_cases": int(ps["case"].nunique()),
        "region_counts": ps["region"].value_counts().to_dict(),
        "region_fine_counts": ps["region_fine"].value_counts().to_dict(),
        "n_mannitol_exposed": int(ps["mannitol_exposed"].sum()),
        "drugs": {},
        "model_predicted_regional_fold_change": 1.0,
        "model_has_region_specific_bbb": False,
    }

    for drug, col in _DRUG_COL.items():
        d: Dict = {
            "column": col,
            "n_nonnull": int(ps[col].notna().sum()),
            "n_zero_below_lloq": int((ps[col] == 0).sum()),
            "by_region": {},
            "contrasts": [],
            "by_region_fine": {},
        }
        for region, g in ps.groupby("region"):
            d["by_region"][region] = _summarise(g[col].tolist())
        for region, g in ps.groupby("region_fine"):
            d["by_region_fine"][region] = _summarise(g[col].tolist())

        for ra, rb in _REGION_PAIRS:
            a = ps.loc[ps["region"] == ra, col].dropna().to_numpy()
            b = ps.loc[ps["region"] == rb, col].dropna().to_numpy()
            if a.size == 0 or b.size == 0:
                continue
            wc = _perm_p_within_case(ps, col, "region", ra, rb,
                                     n_perm=n_perm, seed=seed)
            d["contrasts"].append({
                "region_a": ra, "region_b": rb,
                "n_a": int(a.size), "n_b": int(b.size),
                "mean_a": float(a.mean()), "mean_b": float(b.mean()),
                "median_a": float(np.median(a)), "median_b": float(np.median(b)),
                "ratio_of_means": float(a.mean() / b.mean()) if b.mean() > 0 else float("inf"),
                "difference_of_means": float(a.mean() - b.mean()),
                "cliffs_delta": _cliffs_delta(a, b),
                "hedges_g": _hedges_g(a, b),
                "p_perm_naive": _perm_p_naive(a, b, n_perm=n_perm, seed=seed),
                "p_perm_within_case": wc["p"],
                "n_informative_cases": wc["n_informative_cases"],
                "model_ratio_of_means": 1.0,
            })

        # mannitol sensitivity: repeat the headline CE-vs-non-enhancing contrast
        # on the mannitol-naive subset only.
        naive = ps[~ps["mannitol_exposed"]]
        a = naive.loc[naive["region"] == "contrast_enhancing", col].dropna().to_numpy()
        b = naive.loc[naive["region"] == "non_enhancing", col].dropna().to_numpy()
        if a.size and b.size:
            wc = _perm_p_within_case(naive, col, "region", "contrast_enhancing",
                                     "non_enhancing", n_perm=n_perm, seed=seed)
            d["mannitol_naive_ce_vs_nonenhancing"] = {
                "n_a": int(a.size), "n_b": int(b.size),
                "mean_a": float(a.mean()), "mean_b": float(b.mean()),
                "ratio_of_means": float(a.mean() / b.mean()) if b.mean() > 0 else float("inf"),
                "cliffs_delta": _cliffs_delta(a, b),
                "p_perm_within_case": wc["p"],
                "n_informative_cases": wc["n_informative_cases"],
            }
        out["drugs"][drug] = d

    # Literature-side regional contrast (independent of Oh 2026): methotrexate
    # and adavosertib were each measured in CE vs non-enhancing tissue.
    lit = load_literature_anchors()
    human = lit[(lit["species"] == "human") & lit["tissue_to_plasma_ratio"].notna()]
    lit_rows = []
    for drug, g in human.groupby("drug"):
        ce = g[g["region"] == "contrast_enhancing"]["tissue_to_plasma_ratio"]
        ne = g[g["region"].isin(["nonenhancing", "noncontrast_enhancing"])][
            "tissue_to_plasma_ratio"]
        if len(ce) and len(ne):
            lit_rows.append({
                "drug": str(drug),
                "n_ce": int(len(ce)), "n_ne": int(len(ne)),
                "mean_ce": float(ce.mean()), "mean_ne": float(ne.mean()),
                "ratio_ce_over_ne": float(ce.mean() / ne.mean()) if ne.mean() > 0 else float("inf"),
                "difference_ce_minus_ne": float(ce.mean() - ne.mean()),
                "pmid": "|".join(sorted(set(g["source_pmid"].astype(str)))),
                "method": "|".join(sorted(set(g["method"]))),
            })
    out["literature_regional"] = lit_rows
    return out


# ---------------------------------------------------------------------------
# TEST 3 -- what a region-aware BBB term would have to look like
# ---------------------------------------------------------------------------
def required_bbb_disruption(regional: Dict) -> Dict:
    """Quantify the missing mechanism.

    The model has no region term. Two candidate one-parameter fixes are scored
    against the two Oh 2026 probes plus the Blakeley methotrexate microdialysis
    pair, using the CE-vs-non-enhancing contrast:

      MULTIPLICATIVE   R_region = R_base * f_region
                       -> f = mean_CE / mean_nonCE
      ADDITIVE LEAK    R_region = R_base + L_region
                       -> L = mean_CE - mean_nonCE

    A mechanism that is really a property of the *barrier* (not of the drug)
    should give a drug-INDEPENDENT parameter. The spread across drugs (max/min)
    is the diagnostic: the smaller spread wins.
    """
    obs: List[Dict] = []
    for drug in _DRUG_COL:
        for c in regional["drugs"][drug]["contrasts"]:
            if c["region_a"] == "contrast_enhancing" and c["region_b"] == "non_enhancing":
                obs.append({"drug": drug, "mean_ce": c["mean_a"], "mean_ne": c["mean_b"],
                            "source": "Oh 2026 (homogenate)"})
    for r in regional.get("literature_regional", []):
        if r["drug"] == "methotrexate":
            obs.append({"drug": r["drug"], "mean_ce": r["mean_ce"], "mean_ne": r["mean_ne"],
                        "source": "Blakeley 2008 (microdialysis)"})

    for o in obs:
        o["multiplicative_f"] = o["mean_ce"] / o["mean_ne"] if o["mean_ne"] > 0 else float("inf")
        o["additive_L"] = o["mean_ce"] - o["mean_ne"]

    def _spread(key: str) -> float:
        v = [o[key] for o in obs if math.isfinite(o[key]) and o[key] > 0]
        return (max(v) / min(v)) if v else float("nan")

    mult = [o["multiplicative_f"] for o in obs if math.isfinite(o["multiplicative_f"])]
    addi = [o["additive_L"] for o in obs]
    return {
        "observations": obs,
        "multiplicative_spread": _spread("multiplicative_f"),
        "additive_spread": _spread("additive_L"),
        "multiplicative_mean": float(np.mean(mult)) if mult else float("nan"),
        "additive_mean": float(np.mean(addi)) if addi else float("nan"),
        "additive_min": float(np.min(addi)) if addi else float("nan"),
        "additive_max": float(np.max(addi)) if addi else float("nan"),
        "winner": ("additive_leak" if _spread("additive_L") < _spread("multiplicative_f")
                   else "multiplicative"),
    }


# ---------------------------------------------------------------------------
# TEST 4 -- POST-FIX: re-run the regional test with the barrier term ACTIVE
# ---------------------------------------------------------------------------

#: MRI/pathology region label -> the ``bbb_disruption`` value the model assigns
#: it, via :data:`tmesim.reconstruct.barrier.REGION_DISRUPTION`. The MRI labels
#: are the measured side; the Ivy GAP zone names are the model side; the mapping
#: between them is the anatomic argument documented in ``reconstruct.barrier``.
_REGION_TO_DISRUPTION = {
    "contrast_enhancing": REGION_DISRUPTION["perivascular"],   # 1.00, the CE anchor
    "necrosis": REGION_DISRUPTION["necrotic"],                 # 0.48, measured
    "non_enhancing": REGION_DISRUPTION["margin"],              # 0.00, the reference
    # literature-table spellings of the non-enhancing reference
    "nonenhancing": REGION_DISRUPTION["margin"],
    "noncontrast_enhancing": REGION_DISRUPTION["margin"],
    "peritumoral_brain": REGION_DISRUPTION["margin"],
}


def run_post_fix_regional(regional: Dict, aggregation: str = "geometric") -> Dict:
    """Does the model reproduce the measured regional contrast NOW?

    Pre-fix the answer was trivially no: ``bbb_penetrance`` has no region
    argument, so every predicted regional fold change was exactly 1.000. With
    ``reconstruct.barrier`` active the model predicts, per region,

        R_region = clip(bbb_penetrance + L * D_region, 0, 1)

    with ``L = 0.238`` and ``D`` from :data:`_REGION_TO_DISRUPTION`. This
    function re-scores every contrast in ``regional`` against that prediction.

    **What is and is not a free fit here.** ``L`` WAS fitted to the
    CE-vs-non-enhancing *differences* of these very drugs, so the model
    reproducing a difference of ~0.24 is close to circular and is reported as a
    consistency check, not a prediction. What is NOT circular, and is the real
    test, is:

    * the **ratios**, which depend on the base ``bbb_penetrance`` the leak is
      added to — a quantity fitted to nothing;
    * the **absolute per-region levels**, likewise;
    * the **necrosis** contrasts under a single ``D = 0.48``;
    * whether the ceiling clip at 1.0 breaks the penetrant drug.

    Every one of those can fail with ``L`` exactly right, and the report says so.
    """
    out: Dict = {"leak_constant": L_PARACELLULAR_LEAK,
                 "region_disruption": dict(_REGION_TO_DISRUPTION),
                 "aggregation": aggregation, "drugs": {}, "literature": []}

    base: Dict[str, float] = {}
    for drug in list(_DRUG_COL) + ["methotrexate", "adavosertib", "temozolomide"]:
        try:
            base[drug] = molecular_profile(drug, aggregation=aggregation)[
                "predicted_partition"]
        except KeyError:
            continue

    def _pred(drug: str, region: str) -> float:
        return float(effective_penetrance(base[drug], _REGION_TO_DISRUPTION[region]))

    for drug in _DRUG_COL:
        if drug not in base:
            continue
        d = regional["drugs"][drug]
        entry: Dict = {"bbb_penetrance": base[drug], "by_region": {}, "contrasts": []}
        for region, s in d["by_region"].items():
            if region not in _REGION_TO_DISRUPTION:
                continue
            p = _pred(drug, region)
            entry["by_region"][region] = {
                "measured_mean": s["mean"], "predicted": p, "n": s["n"],
                "fold_error": p / s["mean"] if s["mean"] > 0 else float("inf"),
            }
        for c in d["contrasts"]:
            ra, rb = c["region_a"], c["region_b"]
            if ra not in _REGION_TO_DISRUPTION or rb not in _REGION_TO_DISRUPTION:
                continue
            pa, pb = _pred(drug, ra), _pred(drug, rb)
            entry["contrasts"].append({
                "region_a": ra, "region_b": rb,
                "measured_ratio": c["ratio_of_means"],
                "measured_difference": c["difference_of_means"],
                "model_ratio_prefix": 1.0,
                "model_ratio_postfix": pa / pb if pb > 0 else float("inf"),
                "model_difference_postfix": pa - pb,
                "p_perm_within_case": c["p_perm_within_case"],
                "ratio_error_prefix": (1.0 / c["ratio_of_means"]
                                       if c["ratio_of_means"] > 0 else float("nan")),
                "ratio_error_postfix": ((pa / pb) / c["ratio_of_means"]
                                        if pb > 0 and c["ratio_of_means"] > 0
                                        else float("nan")),
                "ceiling_clipped": bool(base[drug] + L_PARACELLULAR_LEAK
                                        * _REGION_TO_DISRUPTION[ra] > 1.0),
            })
        out["drugs"][drug] = entry

    for r in regional.get("literature_regional", []):
        drug = r["drug"]
        if drug not in base:
            continue
        pa = _pred(drug, "contrast_enhancing")
        pb = _pred(drug, "nonenhancing")
        out["literature"].append({
            "drug": drug, "measured_ratio": r["ratio_ce_over_ne"],
            "measured_difference": r["difference_ce_minus_ne"],
            "model_ratio_postfix": pa / pb if pb > 0 else float("inf"),
            "model_difference_postfix": pa - pb,
            "method": r["method"], "pmid": r["pmid"],
        })

    # Headline: how far off are the two anchor CE-vs-non-enhancing ratios now?
    head = []
    for drug in _DRUG_COL:
        for c in out["drugs"].get(drug, {}).get("contrasts", []):
            if c["region_a"] == "contrast_enhancing" and c["region_b"] == "non_enhancing":
                head.append({"drug": drug, **c})
    out["headline"] = head
    errs = [abs(math.log10(c["ratio_error_postfix"])) for c in head
            if c["ratio_error_postfix"] > 0 and math.isfinite(c["ratio_error_postfix"])]
    out["headline_geometric_mean_abs_ratio_error"] = (
        10.0 ** (sum(errs) / len(errs)) if errs else float("nan"))
    return out


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_figure(score: Dict, regional: Dict, disruption: Dict,
                path: str = REPORT_PNG) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    ps = load_per_sample_ratios()
    fig = plt.figure(figsize=(16.5, 5.4), dpi=150)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 0.95], wspace=0.30)

    # ---- Panel A: predicted vs measured, per drug x region -----------------
    ax = fig.add_subplot(gs[0, 0])
    lo, hi = 1e-3, 1e1
    ax.fill_between([lo, hi], [lo / BAND_OOM, hi / BAND_OOM],
                    [lo * BAND_OOM, hi * BAND_OOM], color="#c9d9ee", alpha=0.45,
                    lw=0, zorder=0)
    ax.fill_between([lo, hi], [lo / BAND_TIGHT, hi / BAND_TIGHT],
                    [lo * BAND_TIGHT, hi * BAND_TIGHT], color="#8fb2dd", alpha=0.40,
                    lw=0, zorder=0)
    ax.plot([lo, hi], [lo, hi], "k-", lw=1.3, zorder=1)

    pred_by_drug = {r["drug"]: r["predicted_partition"] for r in score["rows"]}
    marker = {"contrast_enhancing": "^", "non_enhancing": "v",
              "noncontrast_enhancing": "v", "nonenhancing": "v",
              "peritumoral_brain": "o", "necrosis": "s"}
    colors = plt.get_cmap("tab10")
    dcol = {d: colors(i) for i, d in enumerate(sorted(pred_by_drug))}

    for m in measured_partitions():
        d = m["drug"]
        if d not in pred_by_drug or not math.isfinite(m["mean"]) or m["mean"] <= 0:
            continue
        hollow = "microdialysis" not in m["method"]
        ax.scatter(m["mean"], pred_by_drug[d], s=95,
                   marker=marker.get(m["region"], "P"),
                   facecolor="none" if hollow else dcol[d],
                   edgecolor=dcol[d], linewidths=1.8, zorder=3)
    # Manual label offsets so the drug names do not sit on top of their own points.
    _lbl_off = {"eribulin": (9, 7), "cefazolin": (-6, -17), "methotrexate": (9, -13),
                "temozolomide": (-72, 14), "levetiracetam": (10, 7),
                "adavosertib": (-14, 13)}
    for d, p in pred_by_drug.items():
        row = next(r for r in score["rows"] if r["drug"] == d)
        ax.scatter(row["measured_mean"], p, s=185, marker="*", color=dcol[d],
                   edgecolor="k", linewidths=0.6, zorder=4)
        ax.annotate(d, (row["measured_mean"], p), textcoords="offset points",
                    xytext=_lbl_off.get(d, (7, 6)), fontsize=8.5, color=dcol[d],
                    fontweight="bold", zorder=5)

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lo, hi); ax.set_ylim(1e-2, 1.6)
    ax.axhline(1.0, color="0.45", ls=":", lw=1.0)
    ax.text(1.2e-3, 1.03, "model ceiling: bbb_penetrance <= 1", fontsize=7.2,
            color="0.35")
    ax.set_xlabel("MEASURED tissue:plasma ratio (human glioma)")
    ax.set_ylabel("PREDICTED partition = bbb_penetrance  (A1)")
    ax.set_title("A. Predicted vs measured drug partition\n"
                 f"identity line, {BAND_TIGHT:g}x and {BAND_OOM:g}x bands", fontsize=10)
    ax.legend(handles=[
        Line2D([], [], marker="*", ls="", color="0.3", ms=12, label="drug-level pooled"),
        Line2D([], [], marker="^", ls="", mfc="none", mec="0.3", color="0.3",
               label="contrast-enhancing"),
        Line2D([], [], marker="v", ls="", mfc="none", mec="0.3", color="0.3",
               label="non-enhancing"),
        Line2D([], [], marker="o", ls="", mfc="0.3", mec="0.3", color="0.3",
               label="filled = microdialysis (free)"),
        Line2D([], [], marker="o", ls="", mfc="none", mec="0.3", color="0.3",
               label="hollow = homogenate (total)"),
    ], fontsize=7.0, loc="lower right", framealpha=0.92)
    ax.grid(alpha=0.25, which="both", lw=0.4)

    # ---- Panel B: regional distributions -----------------------------------
    ax = fig.add_subplot(gs[0, 1])
    order = ["contrast_enhancing", "necrosis", "non_enhancing"]
    labels = ["contrast-\nenhancing", "necrosis", "non-\nenhancing"]
    rng = np.random.default_rng(7)
    pos = 0
    ticks, ticklabels = [], []
    for drug, col in _DRUG_COL.items():
        base_c = dcol.get(drug, colors(0))
        for k, reg in enumerate(order):
            v = ps.loc[ps["region"] == reg, col].dropna().to_numpy()
            if v.size == 0:
                continue
            pos += 1
            ax.boxplot([v], positions=[pos], widths=0.6, showfliers=False,
                       medianprops=dict(color="k", lw=1.6),
                       boxprops=dict(color=base_c, lw=1.2),
                       whiskerprops=dict(color=base_c, lw=1.0),
                       capprops=dict(color=base_c, lw=1.0))
            ax.scatter(pos + rng.uniform(-0.17, 0.17, v.size), v, s=9,
                       color=base_c, alpha=0.55, zorder=3, linewidths=0)
            ticks.append(pos); ticklabels.append(labels[k])
        # model prediction: one flat line across all regions of this drug
        p = pred_by_drug.get(drug)
        if p is not None:
            ax.hlines(p, pos - 2.6, pos + 0.45, color=base_c, ls="--", lw=1.8)
            ax.text(pos + 0.5, p, f" model\n {drug[:4]}. {p:.2f}", fontsize=7.2,
                    color=base_c, va="center")
        ax.text(pos - 1.0, 2.18, drug, ha="center", fontsize=9.5,
                color=base_c, fontweight="bold")
        pos += 1.2
    ax.set_xticks(ticks); ax.set_xticklabels(ticklabels, fontsize=7.4)
    ax.set_ylim(-0.06, 2.35)
    ax.set_ylabel("measured tumour:plasma ratio")
    ax.set_title("B. Regional partition, 125 samples / 42 GBM cases\n"
                 "(Oh 2026 Table S3; dashed = model, flat by construction)",
                 fontsize=10)
    ax.grid(alpha=0.25, axis="y", lw=0.4)

    # ---- Panel C: what the missing region term must be ---------------------
    ax = fig.add_subplot(gs[0, 2])
    obs = disruption["observations"]
    names = [o["drug"] for o in obs]
    x = np.arange(len(obs))
    ax.bar(x - 0.2, [o["multiplicative_f"] for o in obs], width=0.38,
           color="#d98d8d", edgecolor="k", lw=0.5, label="multiplicative f\n(CE/nonCE)")
    ax.bar(x + 0.2, [o["additive_L"] for o in obs], width=0.38,
           color="#7fa6c9", edgecolor="k", lw=0.5, label="additive leak L\n(CE-nonCE)")
    ax.axhline(disruption["additive_mean"], color="#2b5d85", ls="--", lw=1.4)
    ax.text(-0.45, disruption["additive_mean"] * 1.30,
            f"mean L = {disruption['additive_mean']:.2f}  (spread "
            f"{disruption['additive_spread']:.1f}x)", fontsize=7.6,
            color="#2b5d85", ha="left")
    ax.axhline(1.0, color="0.6", lw=0.8, ls=":")
    ax.set_yscale("log")
    ax.set_ylim(0.02, 12.0)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=7.6, rotation=18, ha="right")
    ax.set_ylabel("fitted region parameter")
    ax.set_title("C. The missing mechanism\n"
                 f"additive spread {disruption['additive_spread']:.1f}x vs "
                 f"multiplicative {disruption['multiplicative_spread']:.1f}x",
                 fontsize=10)
    ax.legend(fontsize=7.2, loc="upper left")
    ax.grid(alpha=0.25, axis="y", lw=0.4)

    fig.suptitle("tmesim validation run #1 — drug partition in human glioma "
                 "(delivery axis only; regional, not spatial; no response)",
                 fontsize=11.5, y=1.02)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _fmt(x: float, n: int = 3) -> str:
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "n/a"
    return f"{x:.{n}f}"


def build_post_fix(regional: Dict) -> Dict:
    """Assemble everything the POST-FIX report section needs.

    Three scorecards (pre-fix / FIX 1 only / FIX 1+2), the aggregation bake-off,
    the ordering sanity checks, and the re-run regional test.
    """
    score_pre = run_scorecard(aggregation="arithmetic", apply_leak=False)
    score_post_structure = run_scorecard(aggregation="geometric", apply_leak=False)
    score_post = run_scorecard(aggregation="geometric", apply_leak=True)

    bakeoff = [
        ("arithmetic mean (pre-fix defect)", score_pre),
        ("geometric mean, structure-only", score_post_structure),
        ("**geometric mean + leak (CHOSEN)**", score_post),
    ]

    # --- ordering sanity ---------------------------------------------------
    orderings: List[str] = []
    ab = {"mw": 149000.0, "logp": -1.0, "tpsa": 6000.0, "hbd": 200,
          "hba": 400, "is_antibody": True}
    orderings.append(
        f"Antibody (bevacizumab-like) penetrance = {bbb_penetrance(ab):.4f} — "
        "still a hard 0, both via the `is_antibody` gate and the MW >= 1000 gate.")
    sm = {r["drug"]: r["structure_only_partition"] for r in score_post_structure["rows"]}
    ranked = sorted(sm, key=sm.get)
    orderings.append(
        "Least- to most-penetrant among the six small molecules (structure-only): "
        + " < ".join(f"{d} {sm[d]:.3f}" for d in ranked) + ".")
    orderings.append(
        f"Methotrexate is the least penetrant small molecule: "
        f"**{ranked[0] == 'methotrexate'}** (it was also least pre-fix, so the veto "
        "preserves this).")
    # cross-drug separation, the strongest pre-fix positive result
    lev = sm.get("levetiracetam", float("nan"))
    cef = sm.get("cefazolin", float("nan"))
    orderings.append(
        f"BBB-penetrant vs BBB-excluded probe separation (Oh 2026, same tissue "
        f"samples, measured 4.66x pooled / 9.3x in intact-barrier tissue): model now "
        f"predicts {lev / cef:.2f}x (levetiracetam {lev:.3f} / cefazolin {cef:.3f}), "
        "up from 3.65x pre-fix — the veto widens the gap toward the measured "
        "intact-barrier value, which is the right direction.")

    # --- permeability ceiling, before/after -------------------------------
    perm_rows: List[Dict] = []
    for name in _PERM_CHECK:
        if name in GBM_DRUGS:
            mol = from_smiles(name, GBM_DRUGS[name])
        else:
            _verify_smiles(name)
            mol = from_smiles(name, _EXTRA_SMILES[name][0])
        d = mol.descriptors
        perm_rows.append({
            "drug": name,
            "mw": float(d["mw"]),
            "logp": float(d.get("logp", 0.0)),
            "before": _membrane_permeability_unbounded(d),
            "after": float(membrane_permeability(d)),
        })

    return {
        "score_pre": score_pre,
        "score_post_structure": score_post_structure,
        "score_post": score_post,
        "bakeoff": bakeoff,
        "orderings": orderings,
        "regional_post": run_post_fix_regional(regional),
        "permeability": perm_rows,
        "leak": L_PARACELLULAR_LEAK,
    }


def _write_post_fix(A, pf: Dict) -> None:
    """Append the clearly-marked POST-FIX section. ``A`` is a line-appender."""
    pre, post_s, post = pf["score_pre"], pf["score_post_structure"], pf["score_post"]
    reg = pf["regional_post"]
    by = {r["drug"]: r for r in pre["rows"]}
    by_s = {r["drug"]: r for r in post_s["rows"]}
    by_p = {r["drug"]: r for r in post["rows"]}

    A("")
    A("# POST-FIX — the same tests, re-run after both fixes")
    A("")
    A("Two fixes were implemented from the numbers above. **The pre-fix numbers in "
      "sections 1-6 are deliberately kept: the before/after IS the evidence.**")
    A("")
    A("- **FIX 1 — `bbb_penetrance` can now express a veto.** "
      "`tmesim/molecular/permeability.py` combines the four CNS-MPO desirabilities "
      "with a **geometric** mean instead of an arithmetic one. No new parameter; for "
      "a molecule with four equal desirabilities the two forms are identical, so the "
      "change bites only on HETEROGENEOUS descriptor profiles — which is exactly the "
      "veto that was missing.")
    A(f"- **FIX 2 — a barrier-integrity term now exists.** New module "
      f"`tmesim/reconstruct/barrier.py`: "
      f"`effective_penetrance = clip(bbb_penetrance + L * bbb_disruption, 0, 1)` with "
      f"`L = {pf['leak']:.3f}` (mean of the three fitted values 0.211 / 0.230 / 0.272 "
      f"from section 5). Wired through `treatment/drugs.py::delivery_source` behind a "
      f"default-off argument; `core`, `physics` and `sim/engine.py` are untouched.")
    A("")

    # ---- FIX 1: aggregation bake-off ---------------------------------------
    A("## P1. FIX 1 — which aggregation, and why")
    A("")
    A("Both candidates the audit named were evaluated against the same 6-drug "
      "measured scorecard. `pure min` and a `0.7*min + 0.3*mean` variant were "
      "evaluated too and were worse than both; the two finalists are below.")
    A("")
    A("| aggregation | eribulin | worst \\|fold\\| | within 10x | within 3x | GM signed bias | GM \\|fold error\\| |")
    A("|---|---|---|---|---|---|---|")
    for label, s in pf["bakeoff"]:
        er = {r["drug"]: r for r in s["rows"]}["eribulin"]
        A(f"| {label} | {er['fold_error']:.2f}x | {s['worst_abs_fold_error']:.2f}x | "
          f"{s['n_pass_10x']}/{s['n_drugs']} | {s['n_pass_3x']}/{s['n_drugs']} | "
          f"{s['geometric_mean_fold_error']:.3f}x | "
          f"{s['geometric_mean_abs_fold_error']:.3f}x |")
    A("")
    A("Read the last two columns together. **GM signed bias** lets over- and "
      "under-prediction cancel; **GM |fold error|** does not, and is the honest "
      "accuracy number. The pre-fix model's flattering *ex-eribulin* bias of 1.01x "
      "was precisely such a cancellation — temozolomide +2.86x and cefazolin +1.66x "
      "annihilating adavosertib 0.26x — while its unsigned error over the same five "
      "drugs was 2.06x.")
    A("")
    A("**Chosen: the geometric mean.** Against the full post-fix model it gives the "
      "lowest unsigned error (1.99x vs 2.28x for min-weighted), the lowest worst-case "
      "(6.9x vs 17.8x), and it is the only form that brings all 6 drugs inside the "
      "order-of-magnitude band. It also adds no tunable parameter, whereas every "
      "min-weighted form introduces a weight `w` that would have to be fitted to "
      "these same six drugs.")
    A("")

    # ---- FIX 1: before/after scorecard -------------------------------------
    A("## P2. The 6-drug scorecard, BEFORE and AFTER")
    A("")
    A("Three columns, because the two fixes must be separable. *pre-fix* is the "
      "arithmetic mean with no barrier term (what sections 1-6 scored). "
      "*FIX 1 only* is the geometric mean, still with no barrier term. *FIX 1 + 2* is "
      "the model as it now stands, with the additive leak weighted by each drug's "
      "**own measured regional composition** (`measured_sample_disruption`) — because "
      "each measured mean pools whatever regions that study sampled, and scoring a "
      "region-aware model at the wrong regional mix would penalise it for sampling, "
      "not physics.")
    A("")
    A("| drug | measured | modality | disrupted-region fraction of samples | pre-fix pred | pre-fix fold | FIX 1 pred | FIX 1 fold | FIX 1+2 pred | **FIX 1+2 fold** |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    order = ["eribulin", "cefazolin", "methotrexate", "temozolomide",
             "levetiracetam", "adavosertib"]
    for dn in order:
        a, b, c = by[dn], by_s[dn], by_p[dn]
        mod = "microdialysis (free)" if a["free_drug_comparator"] else "homogenate (total)"
        A(f"| {dn} | {a['measured_mean']:.4f} | {mod} | "
          f"{c['sample_mean_bbb_disruption']:.2f} | "
          f"{a['predicted_partition']:.3f} | {a['fold_error']:.2f}x | "
          f"{b['predicted_partition']:.3f} | {b['fold_error']:.2f}x | "
          f"{c['predicted_partition']:.3f} | **{c['fold_error']:.2f}x** |")
    A("")
    A("| aggregate | pre-fix | FIX 1 only | FIX 1 + 2 |")
    A("|---|---|---|---|")
    A(f"| within 10x | {pre['n_pass_10x']}/6 | {post_s['n_pass_10x']}/6 | "
      f"{post['n_pass_10x']}/6 |")
    A(f"| within 3x | {pre['n_pass_3x']}/6 | {post_s['n_pass_3x']}/6 | "
      f"{post['n_pass_3x']}/6 |")
    A(f"| worst \\|fold error\\| | {pre['worst_abs_fold_error']:.2f}x | "
      f"{post_s['worst_abs_fold_error']:.2f}x | {post['worst_abs_fold_error']:.2f}x |")
    A(f"| GM signed bias | {pre['geometric_mean_fold_error']:.3f}x | "
      f"{post_s['geometric_mean_fold_error']:.3f}x | "
      f"{post['geometric_mean_fold_error']:.3f}x |")
    A(f"| GM signed bias, ex-eribulin | "
      f"{pre['geometric_mean_fold_error_ex_eribulin']:.3f}x | "
      f"{post_s['geometric_mean_fold_error_ex_eribulin']:.3f}x | "
      f"{post['geometric_mean_fold_error_ex_eribulin']:.3f}x |")
    A(f"| **GM \\|fold error\\|** | **{pre['geometric_mean_abs_fold_error']:.3f}x** | "
      f"**{post_s['geometric_mean_abs_fold_error']:.3f}x** | "
      f"**{post['geometric_mean_abs_fold_error']:.3f}x** |")
    A(f"| GM \\|fold error\\|, ex-eribulin | "
      f"{pre['geometric_mean_abs_fold_error_ex_eribulin']:.3f}x | "
      f"{post_s['geometric_mean_abs_fold_error_ex_eribulin']:.3f}x | "
      f"{post['geometric_mean_abs_fold_error_ex_eribulin']:.3f}x |")
    A("")
    A("**Honest reading of the middle column.** FIX 1 *alone* does not simply "
      "improve things. It fixes the outlier — eribulin 65.6x -> 6.9x, worst case "
      "65.6x -> 7.9x, and 6/6 drugs inside 10x for the first time — but it converts "
      "the model from an over-predictor into an under-predictor: the ex-eribulin "
      "signed bias moves from 1.01x to 0.36x and the tight 3x band drops from 4/6 to "
      "2/6. That is a real cost and it is not hidden here.")
    A("")
    A("**Why the third column recovers it.** The veto pushes methotrexate and "
      "cefazolin down because both genuinely carry disqualifying descriptors "
      "(TPSA 210 / 156). But their measured means POOL contrast-enhancing tissue "
      "where the barrier is broken — about half the samples for methotrexate, 44% "
      "for cefazolin. The structure-only score cannot represent that tissue at all. "
      "Once the additive leak supplies it, cefazolin lands at 1.02x and methotrexate "
      "at 0.79x. The two fixes are not independent improvements; **FIX 1 is only "
      "correct in the presence of FIX 2**, because the veto is a statement about the "
      "*intact* barrier and most of the measurements are not.")
    A("")
    A("**Requirement check, stated plainly.** The brief required that aggregate "
      "fold-bias not worsen against the pre-fix 2.02x / 1.01x-ex-eribulin. On the "
      "signed measure the full model gives "
      f"{post['geometric_mean_fold_error']:.2f}x (better than 2.02x) but "
      f"{post['geometric_mean_fold_error_ex_eribulin']:.2f}x ex-eribulin, which is "
      "farther from 1.0 than the pre-fix 1.01x. **That requirement is not met on the "
      "signed ex-eribulin measure.** It is met on every non-cancelling measure: "
      f"unsigned error {pre['geometric_mean_abs_fold_error']:.2f}x -> "
      f"{post['geometric_mean_abs_fold_error']:.2f}x overall and "
      f"{pre['geometric_mean_abs_fold_error_ex_eribulin']:.2f}x -> "
      f"{post['geometric_mean_abs_fold_error_ex_eribulin']:.2f}x ex-eribulin, worst case "
      f"{pre['worst_abs_fold_error']:.1f}x -> {post['worst_abs_fold_error']:.1f}x, "
      f"10x band {pre['n_pass_10x']}/6 -> {post['n_pass_10x']}/6. The pre-fix 1.01x "
      "was two large errors cancelling, so it is the wrong thing to protect.")
    A("")

    # ---- orderings ---------------------------------------------------------
    A("### Orderings preserved")
    A("")
    for line in pf["orderings"]:
        A(f"- {line}")
    A("")

    # ---- temozolomide ------------------------------------------------------
    tz_pre, tz_s, tz_p = by["temozolomide"], by_s["temozolomide"], by_p["temozolomide"]
    A("### Temozolomide — reported, not targeted")
    A("")
    A(f"TMZ was **not** tuned for, and nothing in either fix mentions it. It moves "
      f"from {tz_pre['predicted_partition']:.3f} ({tz_pre['fold_error']:.2f}x, a 2.9x "
      f"OVER-prediction of the measured {tz_pre['measured_mean']:.3f}) to "
      f"{tz_s['predicted_partition']:.3f} ({tz_s['fold_error']:.2f}x) under FIX 1 "
      f"alone, and to {tz_p['predicted_partition']:.3f} "
      f"({tz_p['fold_error']:.2f}x) with the leak. **The reshaping improves it**, from "
      f"{max(tz_pre['fold_error'], 1 / tz_pre['fold_error']):.2f}x error to "
      f"{max(tz_p['fold_error'], 1 / tz_p['fold_error']):.2f}x. This is a "
      "consequence of the veto (TMZ's logP -2.08 gives it a desirability of 0.027, "
      "which the arithmetic mean was averaging away), not a TMZ-specific "
      "adjustment.")
    A("")

    # ---- FIX 2: regional ---------------------------------------------------
    A("## P3. FIX 2 — the regional test, re-run")
    A("")
    A(f"Model side is now `R_region = clip(bbb_penetrance + {pf['leak']:.3f} * "
      f"D_region, 0, 1)` with "
      f"`D = 1.00` for contrast-enhancing, `0.48` for necrosis, `0.00` for "
      "non-enhancing. Pre-fix every model ratio in this table was exactly **1.00x**.")
    A("")
    A("**What is circular and what is not.** `L` was fitted to the "
      "CE-vs-non-enhancing *differences* of these same drugs, so the model "
      "reproducing a difference near 0.24 is a consistency check, not a prediction. "
      "The **ratios** are not circular: they depend on the base `bbb_penetrance` the "
      "leak is added to, and that was fitted to nothing.")
    A("")
    A("| drug | contrast | measured ratio | model ratio (pre-fix) | **model ratio (post-fix)** | ratio error | measured diff | model diff | p (within-case) |")
    A("|---|---|---|---|---|---|---|---|---|")
    for drug in _DRUG_COL:
        for c in reg["drugs"].get(drug, {}).get("contrasts", []):
            A(f"| {drug} | {c['region_a']} vs {c['region_b']} | "
              f"{c['measured_ratio']:.2f}x | 1.00x | "
              f"**{c['model_ratio_postfix']:.2f}x** | "
              f"{c['ratio_error_postfix']:.2f}x | "
              f"{c['measured_difference']:+.3f} | "
              f"{c['model_difference_postfix']:+.3f} | "
              f"{c['p_perm_within_case']:.4f} |")
    A("")
    A("Absolute per-region levels (the harder test — these are not differences, so "
      "the fitted `L` cannot rescue them):")
    A("")
    A("| drug | region | n | measured mean | model | fold error |")
    A("|---|---|---|---|---|---|")
    for drug in _DRUG_COL:
        e = reg["drugs"].get(drug, {})
        for r in ("contrast_enhancing", "necrosis", "non_enhancing"):
            s = e.get("by_region", {}).get(r)
            if s:
                A(f"| {drug} | {r} | {s['n']} | {s['measured_mean']:.3f} | "
                  f"{s['predicted']:.3f} | {s['fold_error']:.2f}x |")
    A("")
    A("Independent literature replication (different patients, different assay, and "
      "`L` was fitted on only one of these three):")
    A("")
    A("| drug | measured CE/non-enh | model CE/non-enh (post-fix) | measured diff | model diff | method | PMID |")
    A("|---|---|---|---|---|---|---|")
    for r in reg["literature"]:
        A(f"| {r['drug']} | {r['measured_ratio']:.2f}x | "
          f"{r['model_ratio_postfix']:.2f}x | {r['measured_difference']:+.3f} | "
          f"{r['model_difference_postfix']:+.3f} | {r['method']} | {r['pmid']} |")
    A("")

    # ---- honest verdict ----------------------------------------------------
    head = {c["drug"]: c for c in reg["headline"]}
    lev, cef = head.get("levetiracetam"), head.get("cefazolin")
    A("### Honest verdict on FIX 2")
    A("")
    A("**The model can now express the answer, and it lands the right order of "
      "magnitude — but it does not hit either anchor.**")
    A("")
    if lev:
        A(f"- **levetiracetam CE/non-enhancing: measured 1.45x, model "
          f"{lev['model_ratio_postfix']:.2f}x — UNDERSHOOT "
          f"({lev['ratio_error_postfix']:.2f}x).** The cause is the ceiling, and it "
          f"is diagnosable rather than mysterious: levetiracetam's base penetrance is "
          f"{reg['drugs']['levetiracetam']['bbb_penetrance']:.3f}, so "
          f"{reg['drugs']['levetiracetam']['bbb_penetrance']:.3f} + 0.238 = "
          f"{reg['drugs']['levetiracetam']['bbb_penetrance'] + 0.238:.3f} is clipped "
          f"to 1.0 by the free-drug ceiling. The model delivers a difference of "
          f"{lev['model_difference_postfix']:+.3f} where the measurement says "
          f"{lev['measured_difference']:+.3f}. Note the measured CE mean is 0.875 and "
          "the max is 2.05 — a *homogenate* (total-drug) assay that legitimately "
          "exceeds 1, which assumption A1 structurally cannot. This is the same "
          "ceiling artefact that sinks adavosertib in section 6, not a new failure.")
    if cef:
        A(f"- **cefazolin CE/non-enhancing: measured 4.26x, model "
          f"{cef['model_ratio_postfix']:.2f}x — OVERSHOOT "
          f"({cef['ratio_error_postfix']:.2f}x).** Here the difference is right "
          f"({cef['model_difference_postfix']:+.3f} vs measured "
          f"{cef['measured_difference']:+.3f}, 1.13x) and the ratio is inflated "
          "because the geometric-mean veto puts cefazolin's intact-barrier base at "
          f"{reg['drugs']['cefazolin']['bbb_penetrance']:.3f} against a measured "
          "non-enhancing mean of 0.065 — i.e. the *denominator* is 25% low. The "
          "regional error is inherited from the partition error, not from `L`.")
    mtx = next((r for r in reg["literature"] if r["drug"] == "methotrexate"), None)
    if mtx:
        A(f"- **methotrexate (independent, Blakeley 2008): measured "
          f"{mtx['measured_ratio']:.2f}x, model {mtx['model_ratio_postfix']:.2f}x — "
          f"OVERSHOOT {mtx['model_ratio_postfix'] / mtx['measured_ratio']:.1f}x.** "
          "Same mechanism as cefazolin and worse, because the veto drives "
          f"methotrexate's base to {reg['drugs'].get('methotrexate', {}).get('bbb_penetrance', 0.0225):.4f} "
          "against a measured non-enhancing 0.063. When the denominator is 3x too "
          "small, a correct numerator still gives a 2-3x wrong ratio.")
    A(f"- **Necrosis, the one contrast nothing was fitted to as a ratio.** Measured "
      "necrosis/non-enhancing is 1.12x (levetiracetam) and 3.40x (cefazolin); the "
      "model gives "
      + ", ".join(
          f"{c['model_ratio_postfix']:.2f}x ({d})"
          for d in _DRUG_COL
          for c in reg["drugs"].get(d, {}).get("contrasts", [])
          if c["region_a"] == "necrosis" and c["region_b"] == "non_enhancing")
      + ". The single `D = 0.48` reproduces both, which is the strongest "
        "evidence in this section that the additive form generalises — it was "
        "derived from the necrosis *differences* and it lands the *ratios* for two "
        "drugs 3x apart in penetrance.")
    A("")
    A(f"Aggregate: geometric-mean absolute error on the two anchor CE/non-enhancing "
      f"ratios is **{reg['headline_geometric_mean_abs_ratio_error']:.2f}x**, against "
      f"a pre-fix error of "
      f"{10.0 ** ((abs(math.log10(1.0 / lev['measured_ratio'])) + abs(math.log10(1.0 / cef['measured_ratio']))) / 2):.2f}x "
      "(the pre-fix model predicted 1.00x for everything). The direction, the sign, "
      "the drug-dependence and the necrosis ordering are all now correct; the "
      "magnitudes are off by 1.2-1.4x on the anchors and ~2.5x on the independent "
      "methotrexate replication.")
    A("")
    # ---- FIX 3: the permeability ceiling -----------------------------------
    if pf.get("permeability"):
        pr = pf["permeability"]
        A("## P4. A third defect in the same file — `membrane_permeability` had no ceiling")
        A("")
        A("Surfaced by the GIST imatinib validation (`results/gist_imatinib.md` item 7), "
          "fixed in this pass because it is the **same root pathology** as the eribulin "
          "veto bug: unbounded or unshaped use of a single descriptor. "
          "`membrane_permeability` computed `P = Kp * Dm / d` with "
          "`Kp = 10**logP` clamped only at logP 5, so it returned "
          f"**{next(r['before'] for r in pr if r['drug'] == 'imatinib'):.3g} um/s for "
          "imatinib** — 42 metres per second.")
        A("")
        A("Three bounds, all physical:")
        A("")
        A("1. **logP driving `Kp` clamped at 3.0** (was 5.0). Bilayer/water partitioning "
          "plateaus above logP ~3 while octanol/water keeps climbing — the bilayer is a "
          "more ordered, more polar solvent than bulk octanol. Independent support: "
          "RDKit's Crippen logP for imatinib is 4.59 against an experimental ~3.0, i.e. "
          "the estimator runs ~1.6 log units lipophilic exactly where the clamp bites.")
        A("2. **Unstirred aqueous boundary layer in SERIES** — "
          "`1/P = 1/P_membrane + 1/P_aqueous` with `P_aqueous = D_water / h_UWL`, "
          "`h_UWL = 30 um`. A solute must cross the stagnant water film before it "
          "reaches the bilayer, so this is the *real* physical ceiling on passive "
          "permeability, and it is why measured Caco-2/PAMPA permeabilities saturate. "
          "30 um is the THIN end of the reported 30-1500 um range, i.e. the most "
          "permissive ceiling the literature allows.")
        A("3. **Membrane viscosity 0.1 -> 1.0 Pa*s**, the top of the 0.1-1 Pa*s range "
          "the module already cited. At 0.1 the prefactor `Dm/d` alone was ~2.5e3 um/s, "
          "so every drug with `Kp >= 1` already exceeded any achievable permeability and "
          "the boundary layer would have flattened the whole library to one value.")
        A("")
        A("| drug | logP | MW | BEFORE (um/s) | AFTER (um/s) | physically absurd before? |")
        A("|---|---|---|---|---|---|")
        for r in pr:
            A(f"| {r['drug']} | {r['logp']:.2f} | {r['mw']:.1f} | {r['before']:.4g} | "
              f"{r['after']:.4g} | {'**YES**' if r['before'] > 1e3 else 'no'} |")
        A("")
        mx = max(r["after"] for r in pr)
        A(f"Maximum over the set after the fix: **{mx:.1f} um/s**. Nothing is absurd; "
          "the antibody stays a hard 0.")
        A("")
        A("**Honest limitation — the ceiling binds, and the top of the ranking is now "
          "compressed.** Cisplatin (22.2), imatinib (21.0) and methotrexate (19.6 um/s) "
          "sit within 15% of each other because all three are boundary-layer limited "
          "under this parameterisation, whereas physically they differ by ORDERS of "
          "magnitude in measured apparent permeability — cisplatin is famously almost "
          "membrane-impermeant. Two residual defects cause this and neither is fixed "
          "here: the membrane term is still too permissive in absolute terms (real "
          "Caco-2 Papp for these drugs is 0.001-1 um/s, 1-4 orders below what this "
          "returns), and Crippen logP is simply wrong for some species (cisplatin +0.20 "
          "vs experimental ~-2.2). The output is a bounded, correctly-ordered-at-the-"
          "bottom RANKING, not a calibrated rate. That is a smaller error than 42 m/s, "
          "but it is not zero.")
        A("")
        A("**A gap left deliberately open.** Both functions in `permeability.py` are "
          "neutral-species estimates with no pKa term, so neither can express "
          "pH-partitioning or acidic-compartment ion-trapping. Two independent "
          "validations now converge on that same gap: the GIST run found measured "
          "imatinib tumour:normal liver ratios of 0.023-0.142 lying BELOW the transport "
          "model's structural floor of 0.208 — unreachable at ANY parameterisation, with "
          "hepatic uptake transporters and lysosomal ion-trapping of a weak base named "
          "as the missing mechanisms — and the Moncayo explant run is independently "
          "testing the same effect (weak bases niraparib/rucaparib sequester, neutral "
          "olaparib does not). It is NOT implemented: there is no measured tumour "
          "acidic-phospholipid content in usable units and no tumour-specific lysosomal "
          "parameters in any published PBPK model, so it must be designed rather than "
          "parameterised. Neither fix above forecloses it — a neutral-fraction factor "
          "multiplies `Kp` and a compartment term composes downstream of the series "
          "resistance. The limitation is recorded in the module docstring so it is "
          "visible where the code is read.")
        A("")

    A("**What FIX 2 does NOT fix.** It does not move the model's free-drug ceiling of "
      "1.0, so homogenate ratios above 1 (adavosertib 4.18, levetiracetam CE max "
      "2.05) remain unreachable — that needs the separate total-vs-free bookkeeping "
      "named in section 5 item 3. It does not model efflux, which is still the "
      "largest residual on eribulin. And the disruption FIELD is a documented proxy: "
      "`REGION_DISRUPTION` maps Ivy GAP transcriptional zones onto MRI enhancement by "
      "an anatomic argument (microvascular proliferation is what enhances; leading "
      "edge is what does not), and **no study has assayed Ivy GAP zone labels and "
      "drug partition in the same tissue**. The `core` zone weight of 0.50 is an "
      "unmeasured interpolation and is the weakest number in the module.")
    A("")


def write_report(score: Dict, regional: Dict, disruption: Dict,
                 png: str, path: str = REPORT_MD,
                 post_fix: Optional[Dict] = None) -> str:
    L: List[str] = []
    A = L.append

    A("# Validation run #1 — predicted vs MEASURED drug partition in human glioma")
    A("")
    A("**What this is.** The first time tmesim's molecular/transport layer has been "
      "scored against measured drug concentrations in human brain-tumour tissue. It "
      "tests ONE axis: the **partition / delivery** axis — what fraction of the drug "
      "in plasma actually reaches glioma tissue, and whether that fraction tracks "
      "MRI-defined regional blood-brain-barrier disruption.")
    A("")
    A(f"![predicted vs measured]({os.path.basename(png)})")
    A("")

    # --- provenance -------------------------------------------------------
    A("## 1. Data provenance")
    A("")
    lit = load_literature_anchors()
    ps = load_per_sample_ratios()
    A(f"**Source 1 — literature anchors.** `{os.path.relpath(LIT_CSV, REPO_ROOT)}`, "
      f"{len(lit)} rows, every row carrying a PMID and a DOI. "
      f"{int((lit['species'] == 'human').sum())} human rows / "
      f"{int((lit['species'] == 'rat').sum())} rat rows.")
    A("")
    A("| study | PMID | drug | species | n rows | region(s) | method |")
    A("|---|---|---|---|---|---|---|")
    for (study, pmid, drug, sp, meth), g in lit.groupby(
            ["study", "source_pmid", "drug", "species", "method"]):
        A(f"| {study} | {pmid} | {drug} | {sp} | {len(g)} | "
          f"{', '.join(sorted(set(g['region'])))} | {meth} |")
    A("")
    A("Column meanings: `plasma_cmax` / `tissue_cmax` are peak concentrations, "
      "`plasma_auc` / `tissue_auc` are areas under the concentration-time curve in "
      "the units given by `units_conc` / `units_auc` (they differ per study — "
      "ug/mL vs uM — which is exactly why only the **dimensionless ratio** "
      "`tissue_to_plasma_auc_pct` is used here). That column is divided by 100 on "
      "load to give a fraction. Two rows are not per-subject: Lee 2025 reports a "
      "cohort mean intratumoral:plasma ratio x100 (418 % and 74 %), and Eroglu 2024 "
      "reports a range string `0.13-1.99` (%), from which the geometric midpoint "
      "0.51 % is taken.")
    A("")
    A(f"**Rat rows are excluded from the scorecard.** The "
      f"{int((lit['species'] == 'rat').sum())} Grossman 2012/2013 rows "
      "report intratumoral / contralateral ECF AUC with **no paired plasma AUC**, so "
      "no tissue:plasma ratio exists. They can only support a within-animal ordering "
      "(tumour ECF 3.35 vs contralateral normal-brain ECF 3.22 ug*h/mL pre-bevacizumab "
      "— i.e. essentially NO tumour-vs-normal-brain enrichment for TMZ in that model), "
      "which is itself a useful negative but is not scoreable as a partition.")
    A("")
    A(f"**Source 2 — Oh et al. 2026 Neuro-Oncology (`noaf204`, PMC12962645) Table S3.** "
      f"`{os.path.relpath(PER_SAMPLE_CSV, REPO_ROOT)}`: **{len(ps)} surgical samples "
      f"from {ps['case'].nunique()} high-grade-glioma cases**, a window-of-opportunity "
      "study of BBB permeability heterogeneity. Two probe drugs were given "
      "intra-operatively and assayed in MRI-region-targeted resection specimens:")
    A("")
    A("- `lev_tumor_to_plasma_ratio` = **levetiracetam**, small/neutral, BBB-penetrant "
      f"({int(ps['lev_tumor_to_plasma_ratio'].notna().sum())} values, mean "
      f"{ps['lev_tumor_to_plasma_ratio'].mean():.3f}, range "
      f"{ps['lev_tumor_to_plasma_ratio'].min():.2f}-{ps['lev_tumor_to_plasma_ratio'].max():.2f})")
    A("- `cef_tumor_to_plasma_ratio` = **cefazolin**, anionic cephalosporin, "
      f"BBB-excluded ({int(ps['cef_tumor_to_plasma_ratio'].notna().sum())} values, mean "
      f"{ps['cef_tumor_to_plasma_ratio'].mean():.3f}, range "
      f"{ps['cef_tumor_to_plasma_ratio'].min():.2f}-{ps['cef_tumor_to_plasma_ratio'].max():.2f}) "
      f"— {int((ps['cef_tumor_to_plasma_ratio'] == 0).sum())} of these are reported as "
      "exactly 0.00, i.e. **left-censored at the assay LLOQ**, not true zeros; they are "
      "kept at 0 here, which biases the cefazolin mean and any parametric test "
      "downward and is one reason the nonparametric statistics below are the primary "
      "ones.")
    A("")
    A(f"These ratios are **single-timepoint tissue-concentration : plasma-concentration** "
      "ratios from tissue homogenate, NOT AUC ratios — sampling times (`lev_collection_"
      "time_min`, `cef_collection_time_min`) differ per sample, so each ratio is a "
      "snapshot on a distribution/elimination curve rather than an exposure ratio. "
      "Region flags: `region_contrast_enhancing` (T1-Gd) "
      f"{ps['region_contrast_enhancing'].eq('+').sum()}(+)/"
      f"{ps['region_contrast_enhancing'].eq('-').sum()}(-), `region_flair` "
      f"{ps['region_flair'].eq('+').sum()}(+)/{ps['region_flair'].eq('-').sum()}(-), "
      f"and {int((ps['region'] == 'necrosis').sum())} pathologist-called necrotic "
      "samples (whose C/F flags are N/A). "
      f"`mannitol_status` records minutes since intra-operative mannitol for "
      f"{int(ps['mannitol_exposed'].sum())} samples — mannitol osmotically opens the "
      "barrier, so it is a confounder and a mannitol-naive sensitivity analysis is run.")
    A("")

    # --- mapping ----------------------------------------------------------
    A("## 2. How the model produces a number (assumption A1)")
    A("")
    A("The frozen engine's delivery term is `source = vascular_density * "
      "bbb_penetrance * plasma_mean` (`tmesim/sim/engine.py` L45), and "
      "`tmesim/core/config.py` defines `bbb_penetrance` as the *\"[0,1] fraction of "
      "vascular concentration that crosses into\"* tissue. With `vascular_density` "
      "normalised to 1 and the tissue at quasi-equilibrium with a sustained plasma "
      "level — the Stage-0 plasma-averaging approximation the engine already makes — "
      "the tissue:plasma exposure ratio the model implies is the identity:")
    A("")
    A("> **A1:  R_predicted = `bbb_penetrance(descriptors)`**")
    A("")
    A("`bbb_penetrance` is a CNS-MPO-style structure-only score over MW / logP / TPSA "
      "/ HBD (`tmesim/molecular/permeability.py`). No fitting, no calibration to any "
      "of the data below. Three consequences of A1, stated up front:")
    A("")
    A("1. **R_predicted <= 1 by construction.** Any measured ratio above 1 is "
      "unreachable by the model at any parameterisation. Adavosertib's measured "
      "contrast-enhancing ratio is 4.18. That is a structural ceiling failure, not a "
      "tuning failure.")
    A("2. **R_predicted carries no region term.** `bbb_penetrance` is a function of "
      "the drug's descriptors only. The model's predicted CE-vs-non-enhancing fold "
      "change is **exactly 1.000** for every drug.")
    A("3. **R_predicted is a passive free-drug crossing fraction.** Microdialysis "
      "(Portnow, Blakeley, Jackson, Eroglu) measures *free* extracellular drug and is "
      "the modality-matched comparator. Tissue homogenate (Lee 2025, Oh 2026) measures "
      "*total* drug including intracellular and lipid-bound fractions and can exceed "
      "1. Every scorecard row is labelled with its modality.")
    A("")

    # --- scorecard --------------------------------------------------------
    A("## 3. The scorecard — predicted vs measured, per drug")
    A("")
    A("| drug | MW | logP | TPSA | HBD | **predicted** (A1) | **measured** mean | measured range | n | modality | fold error | 3x | 10x | PMID |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in score["rows"]:
        mod = "microdialysis (free)" if r["free_drug_comparator"] else "homogenate (total)"
        A(f"| {r['drug']} | {r['mw']:.1f} | {r['logp']:.2f} | {r['tpsa']:.1f} | "
          f"{r['hbd']} | **{r['predicted_partition']:.3f}** | "
          f"**{r['measured_mean']:.4f}** | {r['measured_min']:.4f}-{r['measured_max']:.3f} | "
          f"{r['n_measured']} | {mod} | **{r['fold_error']:.2f}x** | "
          f"{'PASS' if r['pass_3x'] else 'FAIL'} | "
          f"{'PASS' if r['pass_10x'] else 'FAIL'} | {r['pmids']} |")
    A("")
    A(f"- **Order-of-magnitude (10x) band: {score['n_pass_10x']}/{score['n_drugs']} pass.**")
    A(f"- Tight (3x) band: {score['n_pass_3x']}/{score['n_drugs']} pass.")
    A(f"- **Aggregate systematic bias: geometric-mean fold error = "
      f"{score['geometric_mean_fold_error']:.2f}x** "
      f"({'OVER' if score['geometric_mean_fold_error'] > 1 else 'UNDER'}-prediction of "
      "brain-tumour drug exposure).")
    A(f"- Excluding the single worst drug ({score['worst_drug']}, "
      f"{score['worst_fold']:.1f}x): geometric-mean fold error = "
      f"{score['geometric_mean_fold_error_ex_worst']:.2f}x — i.e. the aggregate bias is "
      "dominated by one drug, not spread evenly.")
    A(f"- **Modality-matched subset** (microdialysis / free drug only — the only "
      f"apples-to-apples comparison, n={score['n_microdialysis']}): geometric-mean fold "
      f"error = {score['geometric_mean_fold_error_microdialysis_only']:.2f}x, i.e. "
      "*worse* than the pooled figure, again because eribulin lives in this subset. "
      f"Excluding it: {score['geometric_mean_fold_error_microdialysis_ex_worst']:.2f}x "
      "(temozolomide and methotrexate only, n=2 — too few drugs to call this a "
      "calibration).")
    A("")
    A("Per-drug, honestly:")
    A("")
    for r in score["rows"]:
        A(f"- **{r['drug']}** ({r['note']}): predicted {r['predicted_partition']:.3f} vs "
          f"measured {r['measured_mean']:.4f} = **{r['fold_error']:.2f}x**. "
          f"{r['citations']}.")
    A("")

    # --- regional ---------------------------------------------------------
    A("## 4. The regional test")
    A("")
    A(f"n = {regional['n_samples']} samples from {regional['n_cases']} cases. Region "
      f"counts: " + ", ".join(f"`{k}` {v}" for k, v in
                              sorted(regional["region_counts"].items())) + ".")
    A("")
    A("Two permutation tests are reported for every contrast. The **naive** one shuffles "
      "region labels freely across all samples and treats each sample as independent — "
      "it is ANTI-CONSERVATIVE here, because up to 8 samples come from the same tumour. "
      "The **within-case** one shuffles the region label only *within* each patient, "
      "preserving the patient as the clustering unit and testing the paired "
      "within-tumour contrast. **The within-case p is the one to believe.** Both use "
      "20,000 permutations, seeded.")
    A("")
    for drug in _DRUG_COL:
        d = regional["drugs"][drug]
        A(f"### {drug}")
        A("")
        A("| region | n | mean | median | min | max |")
        A("|---|---|---|---|---|---|")
        for reg in ("contrast_enhancing", "necrosis", "non_enhancing"):
            s = d["by_region"].get(reg)
            if not s or s["n"] == 0:
                continue
            A(f"| {reg} | {s['n']} | {s['mean']:.3f} | {s['median']:.3f} | "
              f"{s['min']:.2f} | {s['max']:.2f} |")
        A("")
        A("| contrast | mean A | mean B | ratio A/B | diff A-B | Cliff's d | Hedges' g | p (naive) | p (within-case) | informative cases | MODEL ratio |")
        A("|---|---|---|---|---|---|---|---|---|---|---|")
        for c in d["contrasts"]:
            A(f"| {c['region_a']} vs {c['region_b']} | {c['mean_a']:.3f} | "
              f"{c['mean_b']:.3f} | **{c['ratio_of_means']:.2f}x** | "
              f"{c['difference_of_means']:+.3f} | {c['cliffs_delta']:+.3f} | "
              f"{_fmt(c['hedges_g'], 2)} | {c['p_perm_naive']:.4f} | "
              f"**{c['p_perm_within_case']:.4f}** | {c['n_informative_cases']} | "
              f"{c['model_ratio_of_means']:.2f}x |")
        A("")
        mn = d.get("mannitol_naive_ce_vs_nonenhancing")
        if mn:
            A(f"Mannitol-naive sensitivity (CE vs non-enhancing, n={mn['n_a']} vs "
              f"{mn['n_b']}): ratio {mn['ratio_of_means']:.2f}x, Cliff's d "
              f"{mn['cliffs_delta']:+.3f}, within-case p = {mn['p_perm_within_case']:.4f} "
              f"({mn['n_informative_cases']} informative cases). The regional effect "
              "does not depend on mannitol exposure.")
            A("")
    A("**Independent literature replication of the regional ordering** (different "
      "patients, different drugs, different assay):")
    A("")
    A("| drug | n CE | n non-enh | mean CE | mean non-enh | CE/non-enh | CE - non-enh | method | PMID |")
    A("|---|---|---|---|---|---|---|---|---|")
    for r in regional["literature_regional"]:
        A(f"| {r['drug']} | {r['n_ce']} | {r['n_ne']} | {r['mean_ce']:.3f} | "
          f"{r['mean_ne']:.3f} | **{r['ratio_ce_over_ne']:.2f}x** | "
          f"{r['difference_ce_minus_ne']:+.3f} | {r['method']} | {r['pmid']} |")
    A("")

    # --- missing mechanism -------------------------------------------------
    A("## 5. The model has NO region-specific BBB mechanism — and here is what one "
      "would have to be")
    A("")
    A("Stated plainly: **`bbb_penetrance` is a pure function of the drug's 2D "
      "descriptors. There is no term in the molecular layer, the transport PDE, or "
      "`DrugParams` that represents regional barrier integrity.** The engine's only "
      "spatial modulation of delivery is `vascular_density`, which is vessel *density*, "
      "not barrier *leakiness* — a region can have identical vessel density and a "
      "completely different permeability-surface-area product. So the model's predicted "
      "CE-vs-non-enhancing fold change is exactly 1.00 against measured fold changes of "
      + " and ".join(f"{c['ratio_of_means']:.2f}x ({drug})"
                     for drug in _DRUG_COL
                     for c in regional["drugs"][drug]["contrasts"]
                     if c["region_a"] == "contrast_enhancing"
                     and c["region_b"] == "non_enhancing") +
      ". **The model fails the regional test completely — it cannot even express the "
      "answer.**")
    A("")
    A("What would have to be added, and what its value would have to be. Fit two "
      "candidate one-parameter region terms to the CE-vs-non-enhancing contrast in "
      "three independent drug/assay combinations. A term that is genuinely a property "
      "of the *barrier* rather than of the *drug* must come out drug-INDEPENDENT; the "
      "spread across drugs is the diagnostic.")
    A("")
    A("| drug | source | mean CE | mean non-enh | multiplicative f = CE/nonCE | additive leak L = CE-nonCE |")
    A("|---|---|---|---|---|---|")
    for o in disruption["observations"]:
        A(f"| {o['drug']} | {o['source']} | {o['mean_ce']:.3f} | {o['mean_ne']:.3f} | "
          f"{o['multiplicative_f']:.2f} | {o['additive_L']:.3f} |")
    A("")
    A(f"- multiplicative `f` spread across drugs: **{disruption['multiplicative_spread']:.1f}x**")
    A(f"- additive leak `L` spread across drugs: **{disruption['additive_spread']:.1f}x** "
      f"(L = {disruption['additive_min']:.3f} to {disruption['additive_max']:.3f}, "
      f"mean **{disruption['additive_mean']:.3f}**)")
    A("")
    A(f"**The additive form wins** ({disruption['additive_spread']:.1f}x spread vs "
      f"{disruption['multiplicative_spread']:.1f}x). That is mechanistically "
      "interpretable and it is the concrete specification of the missing feature: "
      "contrast enhancement corresponds to a **drug-independent parallel leak "
      "conductance** that adds roughly a constant increment (~"
      f"{disruption['additive_mean']:.2f} in tissue:plasma units) to the partition "
      "regardless of the molecule's own physicochemistry — exactly what a "
      "size-nonselective paracellular opening should do. A multiplicative barrier "
      "factor cannot work: it would have to be ~"
      + f"{max(o['multiplicative_f'] for o in disruption['observations']):.1f}x"
      " for the barrier-excluded drugs and only ~"
      + f"{min(o['multiplicative_f'] for o in disruption['observations']):.1f}x"
      " for the already-penetrant one, because a drug that is already near-freely "
      "permeable has no headroom left to gain.")
    A("")
    A("Concretely, to pass this test the model would need:")
    A("")
    A("1. A **regional barrier field** — a per-voxel/per-region `bbb_disruption` in "
      "[0,1] derived from contrast enhancement (or a Ktrans-like DCE-MRI parameter), "
      "carried alongside `vascular_density` in the microenvironment fields.")
    A("2. A **delivery term of the form** `source = vascular_density * "
      "(bbb_penetrance + L * bbb_disruption) * plasma`, i.e. an additive leak in "
      "parallel with the drug's intrinsic transcellular route, with `L` calibrated to "
      f"~{disruption['additive_mean']:.2f} and capped so the partition cannot exceed "
      "the free-drug ceiling.")
    A("3. A **separate total-vs-free bookkeeping** so that homogenate-measured ratios "
      "above 1 (adavosertib 4.18) are represented as tissue binding rather than as "
      "impossible barrier crossing.")
    A("")

    # --- verdict ----------------------------------------------------------
    A("## 6. The honest verdict")
    A("")
    A(f"**Partition magnitude: partial pass.** {score['n_pass_10x']}/{score['n_drugs']} "
      f"drugs land inside an order of magnitude and {score['n_pass_3x']}/"
      f"{score['n_drugs']} inside 3x. The aggregate geometric-mean fold error is "
      f"{score['geometric_mean_fold_error']:.2f}x — the model **over-predicts** human "
      "glioma drug exposure — but that aggregate is dominated by one drug "
      f"({score['worst_drug']}, {score['worst_fold']:.0f}x); with it removed the bias "
      f"falls to {score['geometric_mean_fold_error_ex_worst']:.2f}x, i.e. close to "
      "unbiased on the remaining drugs.")
    A("")
    A("**The two failures are two DIFFERENT failures, and neither is a mis-tuned "
      "constant.** They should not be lumped together:")
    A("")
    A("- **Eribulin, 66x over-predicted — a missing efflux mechanism, exposed by an "
      "averaging rule that cannot express a veto.** Eribulin is a 730 Da macrocycle and "
      "a well-characterised P-glycoprotein substrate; its measured peritumoral "
      "brain:plasma ratio is 0.13-1.99 % (Eroglu 2024, PMID 39422741). "
      "`bbb_penetrance` still scores it **0.334**. Two things go wrong at once. "
      "(i) The function combines its four descriptor desirabilities by taking their "
      "*mean*, so eribulin's two disqualifying properties — MW 730 (desirability ~0.01) "
      "and TPSA 146 (~0.01) — are averaged against its two favourable ones — logP 3.44 "
      "(~0.86) and 2 HBD (~0.78) — and the molecule is rescued to a mid-range score "
      "instead of being vetoed. A mean cannot express \"this is disqualifying\"; a "
      "minimum, or a product, could. (ii) Even a correct passive score would be too "
      "high, because the real limit on eribulin is active efflux, which the function's "
      "own docstring concedes it \"cannot capture\". Efflux is the single largest source "
      "of error in the delivery axis, and this is the cheapest available fix to the "
      "second-largest.")
    A("- **Adavosertib, 3.8x under-predicted — a ceiling artefact, not a physics "
      "error.** Its measured contrast-enhancing ratio is 4.18, above the model's "
      "structural maximum of 1.0. But that 4.18 is a *tissue homogenate* number "
      "(total drug, including everything bound inside cells); a free-drug partition of "
      "4.18 would require active uptake INTO brain. The model cannot be right here "
      "because it is not modelling the same quantity. Its predicted 0.639 is in fact "
      "close to adavosertib's *non-enhancing* measured ratio (0.74, 0.86x).")
    A("")
    A("**A caveat in the model's own docstring turns out to be backwards here.** "
      "`permeability.bbb_penetrance` documents that it UNDER-predicts temozolomide, "
      "because TMZ's logP and TPSA are outside the CNS-favourable window even though "
      "TMZ is clinically a BBB-crossing drug. Against these measurements it does the "
      "opposite: predicted 0.558 vs a measured mean of 0.195 across 14 human subjects "
      "(Portnow 2009 + Jackson 2018), a 2.9x OVER-prediction. The documented caveat is "
      "about TMZ's *rank* among CNS-like molecules; as an *absolute partition* the score "
      "is too generous. Both statements can be true, but only the second one is the one "
      "that matters when the number is used as a delivery fraction.")
    A("")
    lev_ne = regional["drugs"]["levetiracetam"]["by_region"]["non_enhancing"]["mean"]
    cef_ne = regional["drugs"]["cefazolin"]["by_region"]["non_enhancing"]["mean"]
    lev_pred = molecular_profile("levetiracetam")["predicted_partition"]
    cef_pred = molecular_profile("cefazolin")["predicted_partition"]
    lev_pooled = next(r["measured_mean"] for r in score["rows"] if r["drug"] == "levetiracetam")
    cef_pooled = next(r["measured_mean"] for r in score["rows"] if r["drug"] == "cefazolin")
    A("**Ordering across drugs: pass — and this is the strongest positive result in the "
      "run.** Oh 2026 measured a BBB-penetrant and a BBB-excluded drug in the *same* "
      f"tissue samples, which removes every between-patient confounder. Measured "
      f"separation: {lev_pooled / cef_pooled:.2f}x pooled ({lev_pooled:.3f} vs "
      f"{cef_pooled:.3f}), rising to {lev_ne / cef_ne:.1f}x in non-enhancing tissue where "
      "the barrier is intact. The model, from 2D structure alone with no fitting, "
      f"predicts {lev_pred:.3f} vs {cef_pred:.3f} = {lev_pred / cef_pred:.2f}x — right "
      "direction, right order of magnitude for the separation. It does under-state the "
      "gap in intact-barrier tissue, which is the same efflux-blindness that sinks "
      "eribulin, just in a milder form.")
    A("")
    ce_bits = []
    for drug in _DRUG_COL:
        for c in regional["drugs"][drug]["contrasts"]:
            if (c["region_a"] == "contrast_enhancing"
                    and c["region_b"] == "non_enhancing"):
                ce_bits.append(f"{c['ratio_of_means']:.2f}x for {drug} "
                               f"(within-case p = {c['p_perm_within_case']:.4f}, "
                               f"Cliff's d {c['cliffs_delta']:+.2f})")
    A("**Regional ordering: outright fail.** The measured regional effect is real and "
      "strong — contrast-enhancing tissue carries " + " and ".join(ce_bits)
      + " the partition of non-enhancing tissue — and it replicates independently in "
      "Blakeley 2008's methotrexate "
      "microdialysis. The model predicts a fold change of exactly 1.00 because it has "
      "no region-specific BBB mechanism at all. Section 5 gives the concrete "
      "specification of what to add.")
    A("")
    A("**Necrosis does not behave as naively expected.** The prior hypothesis was "
      "CE > non-enhancing > necrosis. Measured, necrosis sits "
      "*between* CE and non-enhancing (and for cefazolin it is statistically "
      "indistinguishable from CE), which makes sense — necrotic cores sit inside the "
      "enhancing region and share its disrupted barrier, even though their perfusion is "
      "poor. Any region term added to the model should be driven by barrier integrity, "
      "not by a tissue-viability label.")
    A("")
    A("### What this DOES establish")
    A("")
    A("- The molecular layer's structure-only `bbb_penetrance`, used as a tissue:plasma "
      "partition, is within an order of magnitude of measured human glioma values for "
      f"{score['n_pass_10x']}/{score['n_drugs']} drugs, and correctly rank-orders a "
      "BBB-penetrant against a BBB-excluded drug measured in the *same* human tissue "
      "samples.")
    A("- The measured regional heterogeneity of drug delivery in human GBM is real, "
      "large, and reproducible across two independent studies and four drugs — it is a "
      "legitimate, quantified validation target.")
    A("- The specific missing mechanism has been identified and its magnitude "
      f"estimated (an additive, drug-independent leak of ~{disruption['additive_mean']:.2f} "
      "in contrast-enhancing tissue).")
    A("")
    A("### What this does NOT establish")
    A("")
    A("- **Nothing about response.** No cell kill, no survival, no potency. Delivery is "
      "necessary, not sufficient.")
    A("- **Nothing spatial.** These data are *regional* (an MRI-defined label per "
      "specimen); there are no coordinates, no distance-to-vessel, no gradients. The "
      "sim's perivascular penetration lengths are untouched by this run.")
    A("- **Nothing about time.** The Oh 2026 ratios are single-timepoint snapshots at "
      "varying collection times, not AUC ratios; assumption A1 equates a steady-state "
      "partition to both, which is only approximately right.")
    A("- **Nothing about the free/total distinction.** Two of the six drugs are scored "
      "against tissue homogenate (total drug), which is not the quantity A1 predicts. "
      "Their fold errors are contaminated by tissue binding.")
    A("- **Nothing about n.** Four of the six drugs rest on n <= 14 subjects, two of "
      "them on a single cohort summary. Only levetiracetam/cefazolin have real sample "
      "size, and those two are clustered within 42 patients.")
    A("")
    A("---")
    A("")
    A("Sections 1-6 above are validation run #1 **as it was run, before any fix**. "
      "The predicted column is recomputed live from the frozen pre-fix aggregation "
      "(`_bbb_penetrance_arithmetic`), not pasted, so it cannot drift. Everything "
      "below is the post-fix re-run.")

    if post_fix is not None:
        _write_post_fix(A, post_fix)

    A("")
    A("---")
    A("")
    A(f"Generated by `tmesim/validation/glioma_regional_pk.py`. `core`, `physics` and "
      "`sim.engine` are untouched. All statistics use 20,000 seeded permutations.")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    return path


# ---------------------------------------------------------------------------
def main(n_perm: int = 20000, seed: int = 0) -> Dict:
    for name in _EXTRA_SMILES:
        _verify_smiles(name)

    # Sections 1-6 are the PRE-FIX run, recomputed live from the frozen
    # arithmetic aggregation so the "before" side of the evidence cannot rot.
    score = run_scorecard(aggregation="arithmetic")
    regional = run_regional_test(n_perm=n_perm, seed=seed)
    disruption = required_bbb_disruption(regional)

    post_fix = build_post_fix(regional)
    png = make_figure(score, regional, disruption)
    md = write_report(score, regional, disruption, png, post_fix=post_fix)

    print("=" * 92)
    print("VALIDATION RUN #1 -- predicted vs MEASURED drug partition in human glioma")
    print("=" * 92)
    print(f"\n{'drug':15s}{'pred(A1)':>10s}{'measured':>11s}{'range':>18s}"
          f"{'n':>5s}{'fold':>9s}{'3x':>5s}{'10x':>5s}  modality")
    print("-" * 92)
    for r in score["rows"]:
        mod = "microdialysis" if r["free_drug_comparator"] else "homogenate"
        print(f"{r['drug']:15s}{r['predicted_partition']:>10.3f}"
              f"{r['measured_mean']:>11.4f}"
              f"{r['measured_min']:>9.4f}-{r['measured_max']:<8.3f}"
              f"{r['n_measured']:>5d}{r['fold_error']:>8.2f}x"
              f"{'  Y' if r['pass_3x'] else '  N':>5s}"
              f"{'  Y' if r['pass_10x'] else '  N':>5s}  {mod}")
    print("-" * 92)
    print(f"pass 10x: {score['n_pass_10x']}/{score['n_drugs']}   "
          f"pass 3x: {score['n_pass_3x']}/{score['n_drugs']}   "
          f"geo-mean fold error: {score['geometric_mean_fold_error']:.2f}x   "
          f"(ex-{score['worst_drug']}: "
          f"{score['geometric_mean_fold_error_ex_worst']:.2f}x)")

    print(f"\nREGIONAL TEST  n={regional['n_samples']} samples / "
          f"{regional['n_cases']} cases")
    print("-" * 92)
    for drug in _DRUG_COL:
        d = regional["drugs"][drug]
        print(f"\n{drug}:")
        for reg in ("contrast_enhancing", "necrosis", "non_enhancing"):
            s = d["by_region"].get(reg)
            if s and s["n"]:
                print(f"   {reg:20s} n={s['n']:3d}  mean={s['mean']:.3f}  "
                      f"median={s['median']:.3f}")
        for c in d["contrasts"]:
            print(f"   {c['region_a'][:14]:>15s} vs {c['region_b'][:14]:<15s} "
                  f"ratio={c['ratio_of_means']:5.2f}x  cliff_d={c['cliffs_delta']:+.3f}  "
                  f"p_naive={c['p_perm_naive']:.4f}  p_within_case="
                  f"{c['p_perm_within_case']:.4f}  [model ratio=1.00]")

    print(f"\nMISSING MECHANISM: additive leak spread {disruption['additive_spread']:.1f}x "
          f"vs multiplicative {disruption['multiplicative_spread']:.1f}x -> "
          f"{disruption['winner']} (mean L = {disruption['additive_mean']:.3f})")

    # ---------------- POST-FIX ------------------------------------------
    print("\n" + "=" * 92)
    print("POST-FIX -- geometric-mean veto (FIX 1) + additive barrier leak (FIX 2)")
    print("=" * 92)
    sp, ss, so = (post_fix["score_pre"], post_fix["score_post_structure"],
                  post_fix["score_post"])
    bp = {r["drug"]: r for r in sp["rows"]}
    bs = {r["drug"]: r for r in ss["rows"]}
    bo = {r["drug"]: r for r in so["rows"]}
    print(f"\n{'drug':15s}{'measured':>10s}{'pre-fix':>10s}{'fold':>8s}"
          f"{'FIX1':>9s}{'fold':>8s}{'FIX1+2':>9s}{'fold':>8s}")
    print("-" * 92)
    for dn in ("eribulin", "cefazolin", "methotrexate", "temozolomide",
               "levetiracetam", "adavosertib"):
        a, b, c = bp[dn], bs[dn], bo[dn]
        print(f"{dn:15s}{a['measured_mean']:>10.4f}"
              f"{a['predicted_partition']:>10.3f}{a['fold_error']:>7.2f}x"
              f"{b['predicted_partition']:>9.3f}{b['fold_error']:>7.2f}x"
              f"{c['predicted_partition']:>9.3f}{c['fold_error']:>7.2f}x")
    print("-" * 92)
    for label, s in (("pre-fix ", sp), ("FIX1    ", ss), ("FIX1+2  ", so)):
        print(f"{label} 10x:{s['n_pass_10x']}/6  3x:{s['n_pass_3x']}/6  "
              f"worst:{s['worst_abs_fold_error']:6.2f}x  "
              f"GM-signed:{s['geometric_mean_fold_error']:.3f}x "
              f"(ex-eri {s['geometric_mean_fold_error_ex_eribulin']:.3f}x)  "
              f"GM-|err|:{s['geometric_mean_abs_fold_error']:.3f}x "
              f"(ex-eri {s['geometric_mean_abs_fold_error_ex_eribulin']:.3f}x)")

    rp = post_fix["regional_post"]
    print(f"\nREGIONAL, POST-FIX (L = {rp['leak_constant']}):")
    print("-" * 92)
    for c in rp["headline"]:
        print(f"  {c['drug']:15s} CE/non-enh  measured {c['measured_ratio']:5.2f}x   "
              f"model pre-fix 1.00x -> post-fix {c['model_ratio_postfix']:5.2f}x   "
              f"({c['ratio_error_postfix']:.2f}x off"
              f"{', CEILING-CLIPPED' if c['ceiling_clipped'] else ''})")
    for r in rp["literature"]:
        print(f"  {r['drug']:15s} CE/non-enh  measured {r['measured_ratio']:5.2f}x   "
              f"model post-fix {r['model_ratio_postfix']:5.2f}x   [independent]")
    print(f"  anchor GM |ratio error| = "
          f"{rp['headline_geometric_mean_abs_ratio_error']:.2f}x")

    print(f"\nwrote {md}\nwrote {png}")
    return {"scorecard": score, "regional": regional, "disruption": disruption,
            "post_fix": post_fix, "md": md, "png": png}


if __name__ == "__main__":
    main()
