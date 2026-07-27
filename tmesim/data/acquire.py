"""Best-effort LIVE download of OPEN glioblastoma data into ``./data_cache/``.

This module attempts, in order, to pull *real* open-access glioblastoma data
from public APIs and file stores. Everything is wrapped in ``try/except`` with
explicit ``requests`` timeouts and a single retry: a failed or partial download
must NEVER crash the pipeline — callers degrade gracefully and simply
back to synthetic data if nothing usable landed on disk.

Sources attempted (all verified live against their current URLs/APIs before this
module was written — see ``docs/data_acquisition_plan.md``):

1. ``tcga_gbm_clinical`` — TCGA-GBM clinical records via the GDC REST API
   (https://api.gdc.cancer.gov/cases, open tier). Age, vital status, survival,
   and per-patient treatment/therapeutic-agent lists. Used to build patient
   features and a *survival-based* response proxy (TCGA open tier has no
   per-drug RECIST response, so we approximate — documented in loaders).

2. ``tcga_gbm_expression`` — a modest sample of OPEN STAR gene-count expression
   files (``rna_seq.augmented_star_gene_counts.tsv``) via the GDC files + data
   endpoints, each mapped back to its patient ``submitter_id``. A handful of
   marker genes (MGMT, EGFR, PTEN, VEGFA, CD68, MKI67 ...) parameterise the
   coarse microenvironment.

3. ``ivy_gap`` — Ivy Glioblastoma Atlas Project (Allen Institute) open RNA-seq
   gene-expression matrix zip (anatomic-structure RNA-seq). Contextual
   expression atlas; cached for downstream feature work.

4. ``gdsc_glioma`` — GDSC2 fitted dose-response for glioma cell lines
   (Sanger, open CSV). LN_IC50 / AUC per (cell line, drug) — the cleanest
   drug-response GROUND TRUTH available in the open tier. Streamed and filtered
   to glioma rows (TCGA_DESC in {GBM, LGG}) so the cache stays small.

Run standalone::

    python -m tmesim.data.acquire
"""
from __future__ import annotations

import json
import os
import time
import zipfile
from typing import Dict, List, Optional

import requests

# ---------------------------------------------------------------------------
# Constants — verified live 2026-07-25 (GDC Data Release 45.0, GDSC release 8.5)
# ---------------------------------------------------------------------------
GDC_CASES = "https://api.gdc.cancer.gov/cases"
GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"

IVY_GAP_EXPR_ZIP = (
    "https://glioblastoma.alleninstitute.org/api/v2/"
    "well_known_file_download/305873915"
)

GDSC2_DOSE_RESPONSE_CSV = (
    "https://cog.sanger.ac.uk/cancerrxgene/GDSC_release8.5/"
    "GDSC2_fitted_dose_response_27Oct23.csv"
)

# How much to pull (kept modest on purpose).
N_CLINICAL = 80          # TCGA-GBM patients of clinical metadata
N_EXPRESSION_FILES = 8   # open expression TSVs to download
GLIOMA_TCGA_DESC = ("GBM", "LGG")  # glioma rows to keep from GDSC

# Marker genes extracted from each expression TSV (GBM-relevant biology).
MARKER_GENES = (
    "MGMT",    # methylation/temozolomide-resistance proxy (low expr ~ methylated)
    "EGFR",    # classical amplification
    "PTEN",    # tumour suppressor loss
    "NF1",     # mesenchymal
    "VEGFA",   # angiogenesis / vascularity
    "CD68",    # tumour-associated macrophages (immune infiltrate)
    "MKI67",   # proliferation
    "OLIG2",   # proneural
    "CHI3L1",  # mesenchymal (YKL-40)
    "CD3E",    # T cells
)


def _log(msg: str) -> None:
    print(f"[acquire] {msg}", flush=True)


def _get(url: str, timeout: int, retries: int = 1, **kwargs) -> Optional[requests.Response]:
    """GET with a timeout and a single retry. Returns None on failure (never raises)."""
    last = None
    for attempt in range(retries + 1):
        try:
            resp = requests.get(url, timeout=timeout, **kwargs)
            resp.raise_for_status()
            return resp
        except Exception as exc:  # noqa: BLE001 - degrade gracefully by design
            last = exc
            if attempt < retries:
                time.sleep(1.5)
    _log(f"GET failed {url}: {last}")
    return None


# ---------------------------------------------------------------------------
# Source 1: TCGA-GBM clinical (GDC API)
# ---------------------------------------------------------------------------
def _acquire_tcga_clinical(cache_dir: str, timeout: int) -> Dict:
    out_path = os.path.join(cache_dir, "tcga_gbm_clinical.json")
    fields = ",".join([
        "submitter_id",
        "demographic.gender",
        "demographic.vital_status",
        "demographic.days_to_death",
        "diagnoses.age_at_diagnosis",
        "diagnoses.days_to_last_follow_up",
        "diagnoses.treatments.treatment_type",
        "diagnoses.treatments.treatment_or_therapy",
        "diagnoses.treatments.therapeutic_agents",
    ])
    filt = {
        "op": "in",
        "content": {"field": "project.project_id", "value": ["TCGA-GBM"]},
    }
    params = {
        "filters": json.dumps(filt),
        "fields": fields,
        "format": "json",
        "size": str(N_CLINICAL),
    }
    resp = _get(GDC_CASES, timeout=timeout, params=params)
    if resp is None:
        return {"status": "failed", "path": None, "note": "GDC cases request failed", "n_records": 0}
    try:
        hits = resp.json()["data"]["hits"]
    except Exception as exc:  # noqa: BLE001
        return {"status": "failed", "path": None, "note": f"bad GDC JSON: {exc}", "n_records": 0}
    with open(out_path, "w") as fh:
        json.dump(hits, fh)
    return {
        "status": "ok",
        "path": out_path,
        "note": f"{len(hits)} TCGA-GBM clinical cases",
        "n_records": len(hits),
    }


# ---------------------------------------------------------------------------
# Source 2: TCGA-GBM open expression files (GDC API)
# ---------------------------------------------------------------------------
def _acquire_tcga_expression(cache_dir: str, timeout: int) -> Dict:
    expr_dir = os.path.join(cache_dir, "tcga_expression")
    os.makedirs(expr_dir, exist_ok=True)
    manifest_path = os.path.join(expr_dir, "manifest.json")

    content = [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-GBM"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        {"op": "in", "content": {"field": "access", "value": ["open"]}},
    ]
    # Restrict to patients we already have clinical for, so expression and
    # clinical describe the SAME patients (enables per-patient feature joins).
    clin_path = os.path.join(cache_dir, "tcga_gbm_clinical.json")
    if os.path.exists(clin_path):
        try:
            with open(clin_path) as fh:
                clin = json.load(fh)
            submitters = [c["submitter_id"] for c in clin if c.get("submitter_id")]
            if submitters:
                content.append({
                    "op": "in",
                    "content": {"field": "cases.submitter_id", "value": submitters},
                })
        except Exception as exc:  # noqa: BLE001
            _log(f"could not read clinical for expression join: {exc}")

    filt = {"op": "and", "content": content}
    params = {
        "filters": json.dumps(filt),
        "fields": "file_id,file_name,cases.submitter_id",
        "format": "json",
        "size": str(N_EXPRESSION_FILES),
    }
    resp = _get(GDC_FILES, timeout=timeout, params=params)
    if resp is None:
        return {"status": "failed", "path": None, "note": "GDC files listing failed", "n_records": 0}
    try:
        hits = resp.json()["data"]["hits"]
    except Exception as exc:  # noqa: BLE001
        return {"status": "failed", "path": None, "note": f"bad GDC files JSON: {exc}", "n_records": 0}

    manifest: Dict[str, str] = {}
    n_ok = 0
    for h in hits:
        fid = h.get("file_id")
        cases = h.get("cases") or [{}]
        submitter = (cases[0] or {}).get("submitter_id", f"unknown_{fid}")
        dest = os.path.join(expr_dir, f"{submitter}.star_counts.tsv")
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            manifest[submitter] = os.path.basename(dest)
            n_ok += 1
            continue
        r = _get(f"{GDC_DATA}/{fid}", timeout=timeout, stream=True)
        if r is None:
            continue
        try:
            with open(dest, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 16):
                    if chunk:
                        fh.write(chunk)
            manifest[submitter] = os.path.basename(dest)
            n_ok += 1
        except Exception as exc:  # noqa: BLE001
            _log(f"expression download failed for {submitter}: {exc}")
        finally:
            r.close()

    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh)

    if n_ok == 0:
        return {"status": "failed", "path": None, "note": "no expression files downloaded", "n_records": 0}
    return {
        "status": "ok",
        "path": expr_dir,
        "note": f"{n_ok} open STAR-count expression files",
        "n_records": n_ok,
    }


# ---------------------------------------------------------------------------
# Source 3: Ivy GAP expression matrix zip (Allen Institute)
# ---------------------------------------------------------------------------
def _acquire_ivy_gap(cache_dir: str, timeout: int) -> Dict:
    ivy_dir = os.path.join(cache_dir, "ivy_gap")
    os.makedirs(ivy_dir, exist_ok=True)
    zip_path = os.path.join(ivy_dir, "gene_expression_matrix.zip")

    if not (os.path.exists(zip_path) and os.path.getsize(zip_path) > 1_000_000):
        # Longer timeout for the ~37 MB zip.
        r = _get(IVY_GAP_EXPR_ZIP, timeout=max(timeout, 120), stream=True)
        if r is None:
            return {"status": "failed", "path": None, "note": "Ivy GAP zip download failed", "n_records": 0}
        try:
            with open(zip_path, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 16):
                    if chunk:
                        fh.write(chunk)
        except Exception as exc:  # noqa: BLE001
            r.close()
            return {"status": "failed", "path": None, "note": f"Ivy GAP write failed: {exc}", "n_records": 0}
        finally:
            r.close()

    # Peek inside the zip to count samples (columns) without full extraction.
    n_samples = 0
    try:
        with zipfile.ZipFile(zip_path) as zf:
            names = zf.namelist()
            cols_name = next((n for n in names if "columns" in n.lower()), None)
            if cols_name:
                with zf.open(cols_name) as ch:
                    first = ch.read().decode("utf-8", "ignore")
                # columns-samples.csv: one row per RNA-seq sample.
                n_samples = max(0, len(first.splitlines()) - 1)
    except Exception as exc:  # noqa: BLE001
        _log(f"Ivy GAP zip peek failed (kept file anyway): {exc}")

    return {
        "status": "ok",
        "path": zip_path,
        "note": f"Ivy GAP gene-expression matrix zip ({n_samples} RNA-seq samples)",
        "n_records": n_samples,
    }


# ---------------------------------------------------------------------------
# Source 4: GDSC2 glioma cell-line drug sensitivity (Sanger, open CSV)
# ---------------------------------------------------------------------------
def _acquire_gdsc_glioma(cache_dir: str, timeout: int) -> Dict:
    out_path = os.path.join(cache_dir, "gdsc_glioma_dose_response.csv")

    r = _get(GDSC2_DOSE_RESPONSE_CSV, timeout=max(timeout, 120), stream=True)
    if r is None:
        return {"status": "failed", "path": None, "note": "GDSC CSV download failed", "n_records": 0}

    n_rows = 0
    header = None
    tcga_idx = None
    try:
        with open(out_path, "w") as fh:
            for i, raw in enumerate(r.iter_lines(decode_unicode=True)):
                if raw is None:
                    continue
                if i == 0:
                    header = raw.split(",")
                    try:
                        tcga_idx = header.index("TCGA_DESC")
                    except ValueError:
                        tcga_idx = 6  # known position, fallback
                    fh.write(raw + "\n")
                    continue
                if not raw:
                    continue
                fields = raw.split(",")
                if len(fields) > tcga_idx and fields[tcga_idx] in GLIOMA_TCGA_DESC:
                    fh.write(raw + "\n")
                    n_rows += 1
    except Exception as exc:  # noqa: BLE001
        r.close()
        if n_rows == 0:
            return {"status": "failed", "path": None, "note": f"GDSC stream failed: {exc}", "n_records": 0}
    finally:
        r.close()

    if n_rows == 0:
        return {"status": "failed", "path": None, "note": "no glioma rows in GDSC stream", "n_records": 0}
    return {
        "status": "ok",
        "path": out_path,
        "note": f"{n_rows} glioma (GBM/LGG) cell-line dose-response rows",
        "n_records": n_rows,
    }


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------
_SOURCES = [
    ("tcga_gbm_clinical", _acquire_tcga_clinical),
    ("tcga_gbm_expression", _acquire_tcga_expression),
    ("ivy_gap", _acquire_ivy_gap),
    ("gdsc_glioma", _acquire_gdsc_glioma),
]


def acquire_all(cache_dir: str = "data_cache", timeout: int = 60) -> Dict[str, Dict]:
    """Best-effort download of every open glioblastoma source into ``cache_dir``.

    Returns a status dict::

        {source_name: {"status": "ok"|"failed"|"skipped",
                       "path": <str|None>, "note": <str>, "n_records": <int>}}

    Never raises: any per-source exception is caught and recorded as ``failed``.
    """
    os.makedirs(cache_dir, exist_ok=True)
    status: Dict[str, Dict] = {}
    for name, fn in _SOURCES:
        _log(f"acquiring {name} ...")
        try:
            status[name] = fn(cache_dir, timeout)
        except Exception as exc:  # noqa: BLE001 - final safety net
            status[name] = {
                "status": "failed",
                "path": None,
                "note": f"unhandled error: {exc}",
                "n_records": 0,
            }
        _log(f"  -> {status[name]['status']}: {status[name]['note']}")
    return status


if __name__ == "__main__":
    st = acquire_all()
    print("\n=== acquire_all status ===")
    print(json.dumps(st, indent=2))
