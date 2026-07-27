"""Load REAL, OPEN glioblastoma DRUG-RESPONSE outcome data cached under
``data_cache/drugresponse/``.

This is the *outcome / ground-truth* layer: measured drug sensitivity that a
simulation must eventually be scored against. Every source here is genuinely
open (direct download, no login / DUA). The full ranked landscape — including
the access-GATED sources that would close the gold-standard validation loop —
is in ``docs/real_data/drug_response.md`` and mirrored in :data:`SOURCES`.

Two open datasets are loaded programmatically:

1. **HGCC** — Johansson et al. 2020, *Cell Reports* 32(2):107897,
   "A Patient-Derived Cell Atlas Informs Precision Targeting of Glioblastoma"
   (doi:10.1016/j.celrep.2020.107897). The Human Glioma Cell Culture biobank:
   low-passage, **patient-derived** GBM cell cultures (NOT immortalised lines),
   spanning all transcriptional subtypes. We cache the portal's open JSON
   (http://hgcc.se/portal, no login):
     * ``hgcc_drugresponse.json`` — drug-response matrix, 262 drugs x 85 cultures.
       Value = normalised **survival / viability relative to vehicle control**
       (~1.0 = no effect; **lower = more sensitive**; e.g. the proteasome
       inhibitor Bortezomib averages ~0.42). This is an AUC-like sensitivity
       readout, not an LN_IC50.
     * ``hgcc_bootstrap.json`` — bundled per-culture clinical (molecular
       ``subtype``, sex, age, overall ``survival`` days, vital status) and drug
       metadata. This is the MOLECULAR/CLINICAL pairing for the same cultures.
   Matched omics (expression / CNA / mutation) for the SAME cultures are equally
   open at ``http://hgcc.se/portal/json/{expr,cna,mutation}.json`` and GEO
   (GSE72218 expression, GSE72209 CNA); their large matrices are not cached here
   by default (URLs documented) to keep the cache small.

2. **GDSC** — Sanger GDSC2 fitted dose-response (release 8.5), already cached by
   ``tmesim.data.acquire`` as ``data_cache/gdsc_glioma_dose_response.csv``.
   Glioma **cell lines** (GBM + LGG), quantitative ``LN_IC50`` and ``AUC`` per
   (cell line, drug). Fully open, but immortalised lines with NO microenvironment.

Honest scope: both are *dissociated cell* models — real measured drug response
paired with real molecular profiles, but with **no tissue microenvironment**
and no spatial context. Patient tissue paired with that patient's measured drug
response (organoid/PDX + spatial) is largely GATED — see :data:`SOURCES` and the
companion doc for the exact access steps. Only pandas/numpy/stdlib are used;
nothing here raises on a missing cache (returns an empty frame + a warning).
"""
from __future__ import annotations

import json
import os
from typing import Dict, List

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
_CACHE = os.path.join(_REPO_ROOT, "data_cache", "drugresponse")

HGCC_DRUGRESPONSE_JSON = os.path.join(_CACHE, "hgcc_drugresponse.json")
HGCC_BOOTSTRAP_JSON = os.path.join(_CACHE, "hgcc_bootstrap.json")
GDSC_CSV = os.path.join(_REPO_ROOT, "data_cache", "gdsc_glioma_dose_response.csv")


def _warn(msg: str) -> None:
    print(f"[drug_response] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Ranked source map (open vs gated). Mirrors docs/real_data/drug_response.md.
# ---------------------------------------------------------------------------
SOURCES: List[Dict] = [
    {
        "rank": 1,
        "name": "HGCC / Johansson 2020 (Cell Reports)",
        "kind": "patient-derived GBM cell cultures",
        "open": True,
        "readout": "normalised survival/viability vs control (AUC-like); lower=more sensitive",
        "n_models": 85,
        "n_drugs": 262,
        "molecular_paired": "OPEN: subtype+survival (bundled) + expression/CNA/mutation "
                            "(hgcc.se portal JSON; GEO GSE72218/GSE72209)",
        "microenvironment": False,
        "pmid": "32668248",
        "doi": "10.1016/j.celrep.2020.107897",
        "url": "http://hgcc.se/portal",
        "downloaded": ["hgcc_drugresponse.json", "hgcc_bootstrap.json"],
    },
    {
        "rank": 2,
        "name": "GDSC2 glioma (Sanger, release 8.5)",
        "kind": "immortalised glioma cell lines (GBM+LGG)",
        "open": True,
        "readout": "LN_IC50 + AUC (fitted dose-response)",
        "n_models": 51,
        "n_drugs": 286,
        "molecular_paired": "OPEN via Cell Model Passports / COSMIC (expression, "
                            "mutation, CNA) — but generic lines, not this cohort's tissue",
        "microenvironment": False,
        "pmid": "23180760",  # Yang et al. 2013, GDSC, Nucleic Acids Research
        "doi": "10.1093/nar/gks1111",
        "url": "https://www.cancerrxgene.org",
        "downloaded": ["../gdsc_glioma_dose_response.csv"],
    },
    {
        "rank": 3,
        "name": "Jacob 2020 GBO biobank (Cell)",
        "kind": "patient-derived glioblastoma organoids (GBOs)",
        "open": True,
        "readout": "QUALITATIVE — KI67% / GBO size in figures; genotype->drug "
                   "correlations. NO numeric per-organoid IC50/AUC table.",
        "n_models": 70,
        "n_drugs": 4,  # TMZ+RT, gefitinib, trametinib, everolimus
        "molecular_paired": "OPEN: mutation panel (Table S1), exome AF (Table S3), "
                            "RNA-seq GSE141947",
        "microenvironment": "partial (organoid preserves some cellular diversity; no spatial map)",
        "pmid": "31883794",
        "doi": "10.1016/j.cell.2019.11.036",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141947",
        "downloaded": [f"jacob2020_TableS{i}.xlsx" for i in range(1, 7)],
    },
    {
        "rank": 4,
        "name": "Golebiewska 2020 PDOX/organoid (Acta Neuropathol)",
        "kind": "patient-derived orthotopic xenografts + organoids",
        "open": "molecular OPEN; drug matrix figure-level; models on request",
        "readout": "IC50 + AUC (CellTiter-Glo) for targeted panels + one 41-drug screen",
        "n_models": "~40 PDOX (drug subset smaller)",
        "n_drugs": "~9 targeted + 41-drug array on one model",
        "molecular_paired": "OPEN: aCGH GSE137959, methylation GSE137845, expression "
                            "GSE134470, scRNA GSE128195, seq PRJNA627814",
        "microenvironment": "PDOX preserves stroma/vasculature in vivo",
        "pmid": "33009951",
        "doi": "10.1007/s00401-020-02226-7",
        "url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134470",
        "downloaded": [],
    },
    {
        "rank": 5,
        "name": "HCMI glioblastoma models (GDC)",
        "kind": "next-gen patient-derived models (organoid/cell)",
        "open": True,
        "readout": "NONE systematic — molecular characterisation only; no drug screen",
        "n_models": "dozens of GBM models",
        "n_drugs": 0,
        "molecular_paired": "OPEN: WGS/WXS/RNA-seq/methylation via GDC + open-tier clinical",
        "microenvironment": False,
        "pmid": "",
        "doi": "",
        "url": "https://hcmi-searchable-catalog.nci.nih.gov/",
        "downloaded": [],
    },
    {
        "rank": 6,
        "name": "NCI PDMR (patient-derived models repository)",
        "kind": "patient-derived xenografts / organoids",
        "open": "molecular open; drug response partially open; models via MTA",
        "readout": "in-vivo/ex-vivo drug response (subset public)",
        "n_models": "GBM PDX/organoid subset",
        "n_drugs": "varies",
        "molecular_paired": "OPEN: WES + RNA-seq in public DB",
        "microenvironment": "PDX preserves microenvironment in vivo",
        "pmid": "",
        "doi": "",
        "url": "https://pdmr.cancer.gov",
        "downloaded": [],
    },
    {
        "rank": 7,
        "name": "HTAN GBM spatial multi-omics",
        "kind": "patient tissue (CODEX/MIBI/Visium + single-cell)",
        "open": False,  # GATED — dbGaP + HTAN DUA for imaging/sequence levels
        "readout": "no per-sample measured drug kill; standard-of-care + outcome",
        "n_models": "cohort",
        "n_drugs": 0,
        "molecular_paired": "GATED: dbGaP + HTAN DUA",
        "microenvironment": "TRUE spatial tissue microenvironment",
        "pmid": "",
        "doi": "",
        "url": "https://data.humantumoratlas.org",
        "downloaded": [],
    },
]


# ---------------------------------------------------------------------------
# HGCC (best open patient-derived paired benchmark)
# ---------------------------------------------------------------------------
def load_hgcc_clinical(cache_dir: str = _CACHE) -> pd.DataFrame:
    """Per-culture clinical / molecular-subtype table for the HGCC cohort.

    Columns: ``culture, subtype, sex, age, survival_days, vital_status``.
    ``subtype`` is the transcriptional subtype (CL=classical, MS=mesenchymal,
    PN=proneural, NL=neural, NA=unassigned). Returns an empty frame if the cache
    is missing.
    """
    path = os.path.join(cache_dir, "hgcc_bootstrap.json")
    if not os.path.exists(path):
        _warn(f"missing {path}; returning empty clinical frame")
        return pd.DataFrame(
            columns=["culture", "subtype", "sex", "age", "survival_days", "vital_status"]
        )
    with open(path) as fh:
        boot = json.load(fh)
    rows = boot.get("clinical", [])
    df = pd.DataFrame(rows).rename(
        columns={"name": "culture", "survival": "survival_days", "vitalstat": "vital_status"}
    )
    keep = ["culture", "subtype", "sex", "age", "survival_days", "vital_status"]
    return df[[c for c in keep if c in df.columns]]


def load_hgcc(cache_dir: str = _CACHE, join_clinical: bool = True) -> pd.DataFrame:
    """Long-format HGCC drug-response table (the best OPEN paired benchmark).

    One row per (culture, drug). Columns::

        culture, drug, response, [subtype, sex, age, survival_days, vital_status]

    ``response`` = normalised survival/viability vs vehicle control from the
    HGCC portal (~1.0 = no effect; **lower = more drug-sensitive**). This is an
    AUC-like sensitivity metric, not an LN_IC50. When ``join_clinical`` is True,
    the per-culture molecular subtype + survival columns are joined on.

    Returns an empty frame (never raises) if the cache is missing.
    """
    path = os.path.join(cache_dir, "hgcc_drugresponse.json")
    if not os.path.exists(path):
        _warn(f"missing {path}; run the drug_response acquisition first")
        return pd.DataFrame(columns=["culture", "drug", "response"])
    with open(path) as fh:
        matrix = json.load(fh)  # {drug: {culture: value}}
    records = [
        {"culture": culture, "drug": drug, "response": float(val)}
        for drug, per_culture in matrix.items()
        for culture, val in per_culture.items()
        if val is not None
    ]
    df = pd.DataFrame.from_records(records)
    if join_clinical and not df.empty:
        clin = load_hgcc_clinical(cache_dir)
        if not clin.empty:
            df = df.merge(clin, on="culture", how="left")
    return df


# ---------------------------------------------------------------------------
# GDSC (open quantitative cell-line reference; already cached by acquire.py)
# ---------------------------------------------------------------------------
def load_gdsc(gdsc_csv: str = GDSC_CSV, gbm_only: bool = False) -> pd.DataFrame:
    """Cached GDSC2 glioma dose-response as a tidy frame.

    Columns::

        cell_line, sanger_model_id, tcga_desc, drug, target, pathway, ln_ic50, auc

    ``gbm_only=True`` restricts to ``TCGA_DESC == 'GBM'`` (drops LGG). Returns an
    empty frame if the cache is missing.
    """
    if not os.path.exists(gdsc_csv):
        _warn(f"missing {gdsc_csv}; run tmesim.data.acquire")
        return pd.DataFrame(
            columns=["cell_line", "sanger_model_id", "tcga_desc", "drug",
                     "target", "pathway", "ln_ic50", "auc"]
        )
    df = pd.read_csv(gdsc_csv)
    if gbm_only:
        df = df[df["TCGA_DESC"] == "GBM"]
    out = df.rename(
        columns={
            "CELL_LINE_NAME": "cell_line",
            "SANGER_MODEL_ID": "sanger_model_id",
            "TCGA_DESC": "tcga_desc",
            "DRUG_NAME": "drug",
            "PUTATIVE_TARGET": "target",
            "PATHWAY_NAME": "pathway",
            "LN_IC50": "ln_ic50",
            "AUC": "auc",
        }
    )
    keep = ["cell_line", "sanger_model_id", "tcga_desc", "drug",
            "target", "pathway", "ln_ic50", "auc"]
    return out[[c for c in keep if c in out.columns]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
def summarize() -> Dict:
    """Return + print a short summary of what actually loaded from the cache."""
    hgcc = load_hgcc()
    gdsc = load_gdsc()
    summary = {
        "hgcc": {
            "rows": int(len(hgcc)),
            "cultures": int(hgcc["culture"].nunique()) if not hgcc.empty else 0,
            "drugs": int(hgcc["drug"].nunique()) if not hgcc.empty else 0,
        },
        "gdsc": {
            "rows": int(len(gdsc)),
            "cell_lines": int(gdsc["cell_line"].nunique()) if not gdsc.empty else 0,
            "drugs": int(gdsc["drug"].nunique()) if not gdsc.empty else 0,
        },
    }
    print("=== real drug-response cache ===")
    print(f"HGCC (patient-derived, OPEN, best benchmark): "
          f"{summary['hgcc']['rows']} rows, "
          f"{summary['hgcc']['cultures']} cultures x {summary['hgcc']['drugs']} drugs")
    print(f"GDSC (cell lines, OPEN): "
          f"{summary['gdsc']['rows']} rows, "
          f"{summary['gdsc']['cell_lines']} lines x {summary['gdsc']['drugs']} drugs")
    print(f"sources catalogued: {len(SOURCES)} "
          f"({sum(1 for s in SOURCES if s['open'] is True)} fully open, "
          f"rest partial/gated — see docs/real_data/drug_response.md)")
    return summary


if __name__ == "__main__":
    summarize()
