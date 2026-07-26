"""Parse whatever ``acquire.py`` cached into ``list[PatientRecord]``.

Real open glioblastoma data has NO spatial imaging attached (TCGA-GBM open tier
is bulk expression + clinical; GDSC is bulk cell-line drug screens). So we cannot
reconstruct a true tissue microenvironment. Instead we build a **coarse,
documented-approximate** :class:`Microenvironment` per record whose spatial
layout is *parameterised by the real molecular/clinical features*:

  * tumour-core radius, immune fraction and vascularity are driven by real
    marker-gene expression (VEGFA -> vascular density, CD68/CD3E -> immune
    infiltrate, MKI67 -> tumour mass) and clinical fields;
  * per-cell ``sensitivity`` is anchored to the real drug-response signal
    (GDSC LN_IC50 for cell lines; MGMT-methylation + survival proxy for TCGA);
  * ``observed_response`` is set from the real response/sensitivity data.

This is explicitly an APPROXIMATION: the *biology-to-response* labels are real,
but the 2-D spatial arrangement is synthetic. It exists so the Stage-0 pipeline
can be exercised on genuine molecular cohorts; it is NOT a substitute for the
controlled-access spatial atlases documented in ``docs/data_acquisition_plan.md``.

If nothing usable was cached, :func:`load_real_cohort` returns ``[]`` and the
caller falls back to the fully synthetic cohort. Nothing here ever raises.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from typing import Dict, List, Optional

import numpy as np

from ..core import (
    CellPopulation,
    CellState,
    CellType,
    Domain,
    Microenvironment,
    PatientRecord,
    TissueFields,
)
from .acquire import MARKER_GENES

# Grid for the coarse reconstructed microenvironment.
_GRID = (60, 60)
_DX = 15.0
_N_CELLS = 1400


def _log(msg: str) -> None:
    print(f"[loaders] {msg}", flush=True)


def _seed_from(name: str) -> int:
    """Deterministic 32-bit seed from a patient/cell-line id."""
    return int(hashlib.md5(name.encode()).hexdigest()[:8], 16)


# ---------------------------------------------------------------------------
# Coarse microenvironment builder (documented approximation)
# ---------------------------------------------------------------------------
def _build_coarse_micro(
    patient_id: str,
    sensitivity_mean: float,
    vascularity: float,
    immune_frac: float,
    tumor_burden: float,
    metadata: Dict,
) -> Microenvironment:
    """Build a synthetic 2-D TME whose parameters come from real features.

    Parameters
    ----------
    sensitivity_mean : float in [0,1]
        Anchor for per-cell intrinsic drug sensitivity (from real response data).
    vascularity : float in [0,1]
        Overall vessel density (driven by VEGFA / real angiogenesis signal).
    immune_frac : float in [0,1]
        Fraction of non-tumour cells that are immune (CD68/CD3E signal).
    tumor_burden : float in [0,1]
        Relative size of the tumour core (proliferation signal).
    """
    rng = np.random.default_rng(_seed_from(patient_id))
    ny, nx = _GRID
    domain = Domain(nx=nx, ny=ny, dx=_DX)
    X, Y = domain.coords_um()
    cx, cy = X.mean(), Y.mean()
    R = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    Rmax = R.max()
    core_radius = Rmax * float(np.clip(0.25 + 0.45 * tumor_burden, 0.2, 0.85))

    # --- Continuum fields -------------------------------------------------
    # Vascular density: baseline vascularity, SPARSER in the hypoxic core.
    core_deficit = np.exp(-(R / (0.6 * core_radius + 1e-9)) ** 2)
    vascular = vascularity * (1.0 - 0.7 * core_deficit)
    vascular += 0.05 * rng.standard_normal(vascular.shape)
    vascular = np.clip(vascular, 0.0, 1.0)

    # ECM density: denser in the tumour core (raises tortuosity).
    ecm = 0.3 + 0.5 * np.exp(-(R / (0.9 * core_radius + 1e-9)) ** 2)
    ecm += 0.03 * rng.standard_normal(ecm.shape)
    ecm = np.clip(ecm, 0.0, 1.0)

    fields = TissueFields(ecm_density=ecm, vascular_density=vascular)

    # --- Cell population --------------------------------------------------
    n = _N_CELLS
    # Place ~65% cells inside/around the tumour core, rest scattered as stroma/immune.
    n_tumor = int(n * float(np.clip(0.4 + 0.4 * tumor_burden, 0.35, 0.8)))
    n_other = n - n_tumor

    # Tumour cells: gaussian cloud around the core (infiltrative tail).
    ang = rng.uniform(0, 2 * np.pi, n_tumor)
    rad = np.abs(rng.normal(0, 0.55 * core_radius, n_tumor))
    tx = np.clip(cx + rad * np.cos(ang), 0, domain.extent_um[0] - 1)
    ty = np.clip(cy + rad * np.sin(ang), 0, domain.extent_um[1] - 1)

    # Other cells: uniform over the domain.
    ox = rng.uniform(0, domain.extent_um[0] - 1, n_other)
    oy = rng.uniform(0, domain.extent_um[1] - 1, n_other)

    xs = np.concatenate([tx, ox])
    ys = np.concatenate([ty, oy])

    ctype = np.empty(n, dtype=np.int8)
    ctype[:n_tumor] = CellType.TUMOR
    # Split "other" into stroma / immune / endothelial (probs normalised to 1).
    p_immune = float(np.clip(immune_frac, 0.05, 0.7))
    p_endo = 0.15
    p_stroma = max(0.05, 1.0 - p_immune - p_endo)
    probs = np.array([p_stroma, p_immune, p_endo], dtype=float)
    probs /= probs.sum()
    other_types = rng.choice(
        [CellType.STROMAL, CellType.IMMUNE, CellType.ENDOTHELIAL],
        size=n_other,
        p=probs,
    ).astype(np.int8)
    ctype[n_tumor:] = other_types

    state = np.full(n, int(CellState.VIABLE), dtype=np.int8)

    # Per-cell sensitivity anchored to the real response, with spatial + clonal
    # variation; hypoxic/core cells slightly less sensitive.
    col, row = domain.to_index(xs, ys)
    local_vasc = vascular[row, col]
    sens = (
        sensitivity_mean
        + 0.10 * rng.standard_normal(n)
        + 0.10 * (local_vasc - vascular.mean())  # better perfused -> more sensitive
    )
    sens = np.clip(sens, 0.02, 0.98)

    cells = CellPopulation(
        x=xs, y=ys, cell_type=ctype, state=state, sensitivity=sens,
    )

    return Microenvironment(
        domain=domain, fields=fields, cells=cells,
        patient_id=patient_id, metadata=metadata,
    )


# ---------------------------------------------------------------------------
# Expression TSV parsing (marker genes only)
# ---------------------------------------------------------------------------
def _parse_expression_tsv(path: str) -> Dict[str, float]:
    """Return {gene_name: fpkm_unstranded} for MARKER_GENES only. Robust."""
    wanted = set(MARKER_GENES)
    out: Dict[str, float] = {}
    try:
        with open(path, "r") as fh:
            header_cols = None
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if header_cols is None:
                    header_cols = parts
                    continue
                if len(parts) < 2 or parts[0].startswith("N_"):
                    continue
                gene = parts[1]
                if gene in wanted:
                    try:
                        # fpkm_unstranded is column index 7 in the STAR-count TSV.
                        out[gene] = float(parts[7])
                    except (IndexError, ValueError):
                        pass
                if len(out) == len(wanted):
                    break
    except Exception as exc:  # noqa: BLE001
        _log(f"expression parse failed {path}: {exc}")
    return out


# ---------------------------------------------------------------------------
# Loader: TCGA-GBM clinical + expression -> PatientRecords
# ---------------------------------------------------------------------------
def _load_tcga(cache_dir: str) -> List[PatientRecord]:
    clin_path = os.path.join(cache_dir, "tcga_gbm_clinical.json")
    if not os.path.exists(clin_path):
        return []
    try:
        with open(clin_path) as fh:
            cases = json.load(fh)
    except Exception as exc:  # noqa: BLE001
        _log(f"clinical json unreadable: {exc}")
        return []

    # Load expression manifest + parse marker genes per patient (if available).
    expr_dir = os.path.join(cache_dir, "tcga_expression")
    expr_by_patient: Dict[str, Dict[str, float]] = {}
    manifest_path = os.path.join(expr_dir, "manifest.json")
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path) as fh:
                manifest = json.load(fh)
            for submitter, fname in manifest.items():
                genes = _parse_expression_tsv(os.path.join(expr_dir, fname))
                if genes:
                    # match on the 12-char patient barcode prefix
                    expr_by_patient[submitter[:12]] = genes
        except Exception as exc:  # noqa: BLE001
            _log(f"expression manifest parse failed: {exc}")

    # Cohort-level MGMT median (for methylation proxy) once we have expression.
    mgmt_vals = [g["MGMT"] for g in expr_by_patient.values() if "MGMT" in g]
    mgmt_median = float(np.median(mgmt_vals)) if mgmt_vals else None

    # First pass: gather survival for the median-split response proxy.
    survivals: List[float] = []
    parsed: List[Dict] = []
    for c in cases:
        submitter = c.get("submitter_id", "")
        demo = c.get("demographic") or {}
        diags = c.get("diagnoses") or [{}]
        d0 = diags[0] if diags else {}
        vital = demo.get("vital_status")
        days_death = demo.get("days_to_death")
        days_fu = None
        for d in diags:
            if d.get("days_to_last_follow_up") is not None:
                days_fu = d.get("days_to_last_follow_up")
                break
        surv = None
        if days_death is not None:
            surv = float(days_death)
        elif days_fu is not None:
            surv = float(days_fu)
        age_days = None
        for d in diags:
            if d.get("age_at_diagnosis") is not None:
                age_days = d.get("age_at_diagnosis")
                break
        # Treatment agents (flatten across diagnoses).
        agents = []
        for d in diags:
            for t in (d.get("treatments") or []):
                ag = t.get("therapeutic_agents")
                if ag and ag not in ("Not Reported", "null"):
                    agents.append(str(ag))
        parsed.append({
            "submitter": submitter,
            "vital": vital,
            "surv": surv,
            "age_days": age_days,
            "gender": demo.get("gender"),
            "agents": agents,
        })
        if surv is not None:
            survivals.append(surv)

    surv_median = float(np.median(survivals)) if survivals else None

    records: List[PatientRecord] = []
    for p in parsed:
        submitter = p["submitter"]
        pid12 = submitter[:12]
        genes = expr_by_patient.get(pid12, {})

        age_years = (p["age_days"] / 365.25) if p["age_days"] else 60.0

        # --- response proxy (documented approximation) -------------------
        # TCGA open tier lacks per-drug RECIST response. Proxy: a patient is a
        # "responder" if their overall survival exceeds the cohort median
        # (longer survival ~ better therapy response, imperfectly).
        if p["surv"] is not None and surv_median is not None:
            observed = 1.0 if p["surv"] >= surv_median else 0.0
        else:
            observed = 0.0

        # MGMT methylation proxy: low MGMT expression ~ methylated ~ TMZ-sensitive.
        mgmt_expr = genes.get("MGMT")
        if mgmt_expr is not None and mgmt_median is not None:
            mgmt_methylated = bool(mgmt_expr < mgmt_median)
        else:
            mgmt_methylated = False

        # Feature-driven microenvironment parameters (normalised heuristics).
        vegfa = genes.get("VEGFA", 8.0)
        cd68 = genes.get("CD68", 8.0)
        cd3e = genes.get("CD3E", 2.0)
        mki67 = genes.get("MKI67", 5.0)
        vascularity = float(np.clip(0.25 + 0.02 * vegfa, 0.15, 0.8))
        immune_frac = float(np.clip(0.1 + 0.01 * (cd68 + cd3e), 0.1, 0.7))
        tumor_burden = float(np.clip(0.3 + 0.05 * mki67, 0.3, 0.9))
        # sensitivity anchor: methylated MGMT -> more sensitive to TMZ-like drug.
        sensitivity_mean = 0.55 if mgmt_methylated else 0.35

        features = {
            "source_tcga": 1.0,
            "age_years": float(age_years),
            "is_male": 1.0 if p["gender"] == "male" else 0.0,
            "is_dead": 1.0 if p["vital"] == "Dead" else 0.0,
            "survival_days": float(p["surv"]) if p["surv"] is not None else float("nan"),
            "n_therapeutic_agents": float(len(set(p["agents"]))),
            "received_temozolomide": 1.0 if any("temozolomide" in a.lower() for a in p["agents"]) else 0.0,
            "mgmt_methylated": 1.0 if mgmt_methylated else 0.0,
        }
        for g in MARKER_GENES:
            features[f"expr_{g}"] = float(genes.get(g, float("nan")))

        metadata = {
            "source": "TCGA-GBM (GDC open)",
            "subtype": "mesenchymal" if genes.get("CHI3L1", 0) > genes.get("OLIG2", 0) else "proneural",
            "mgmt_methylated": mgmt_methylated,
            "expression": genes,
            "therapeutic_agents": sorted(set(p["agents"])),
            "response_proxy": "overall-survival median split (approximation)",
        }
        micro = _build_coarse_micro(
            patient_id=f"TCGA::{submitter}",
            sensitivity_mean=sensitivity_mean,
            vascularity=vascularity,
            immune_frac=immune_frac,
            tumor_burden=tumor_burden,
            metadata=metadata,
        )
        records.append(PatientRecord(
            patient_id=f"TCGA::{submitter}",
            micro=micro,
            observed_response=observed,
            features=features,
        ))
    _log(f"built {len(records)} TCGA-GBM patient records")
    return records


# ---------------------------------------------------------------------------
# Loader: GDSC glioma cell lines -> PatientRecords (drug-response ground truth)
# ---------------------------------------------------------------------------
_REFERENCE_DRUG = "Temozolomide"


def _load_gdsc(cache_dir: str) -> List[PatientRecord]:
    path = os.path.join(cache_dir, "gdsc_glioma_dose_response.csv")
    if not os.path.exists(path):
        return []

    # Aggregate per cell line: reference-drug LN_IC50 + mean AUC across all drugs.
    per_line: Dict[str, Dict] = {}
    try:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                line = row.get("CELL_LINE_NAME")
                if not line:
                    continue
                d = per_line.setdefault(line, {
                    "aucs": [], "tcga": row.get("TCGA_DESC"),
                    "ref_ln_ic50": None, "ref_auc": None, "n_drugs": 0,
                })
                try:
                    d["aucs"].append(float(row["AUC"]))
                except (KeyError, ValueError):
                    pass
                d["n_drugs"] += 1
                if row.get("DRUG_NAME") == _REFERENCE_DRUG:
                    try:
                        d["ref_ln_ic50"] = float(row["LN_IC50"])
                        d["ref_auc"] = float(row["AUC"])
                    except (KeyError, ValueError):
                        pass
    except Exception as exc:  # noqa: BLE001
        _log(f"GDSC csv parse failed: {exc}")
        return []

    if not per_line:
        return []

    # Response label from the reference drug LN_IC50 where present, else mean AUC.
    ref_ic50s = [d["ref_ln_ic50"] for d in per_line.values() if d["ref_ln_ic50"] is not None]
    ic50_median = float(np.median(ref_ic50s)) if ref_ic50s else None

    records: List[PatientRecord] = []
    for line, d in sorted(per_line.items()):
        mean_auc = float(np.mean(d["aucs"])) if d["aucs"] else float("nan")
        ref_ic50 = d["ref_ln_ic50"]
        ref_auc = d["ref_auc"]

        # observed_response in [0,1]: 1 = sensitive/responder.
        # Prefer the reference-drug signal; lower LN_IC50 == more sensitive.
        if ref_ic50 is not None and ic50_median is not None:
            observed = 1.0 if ref_ic50 < ic50_median else 0.0
            sensitivity_mean = float(np.clip(0.5 - 0.15 * (ref_ic50 - ic50_median), 0.1, 0.9))
        elif not np.isnan(mean_auc):
            # AUC in [0,1]; lower AUC (more killing) == more sensitive.
            observed = 1.0 if mean_auc < 0.85 else 0.0
            sensitivity_mean = float(np.clip(1.0 - mean_auc, 0.1, 0.9))
        else:
            continue

        # Cell-line screens have no tissue architecture; use neutral-but-varied
        # spatial parameters seeded by the line name (documented approximation).
        rng = np.random.default_rng(_seed_from(line))
        features = {
            "source_gdsc": 1.0,
            "ref_drug_ln_ic50": float(ref_ic50) if ref_ic50 is not None else float("nan"),
            "ref_drug_auc": float(ref_auc) if ref_auc is not None else float("nan"),
            "mean_auc_all_drugs": mean_auc,
            "n_drugs_screened": float(d["n_drugs"]),
            "is_gbm": 1.0 if d["tcga"] == "GBM" else 0.0,
        }
        metadata = {
            "source": "GDSC2 (Sanger open) glioma cell line",
            "cell_line": line,
            "tcga_desc": d["tcga"],
            "reference_drug": _REFERENCE_DRUG,
            "response_ground_truth": "cell-line LN_IC50 median split (real drug screen)",
        }
        micro = _build_coarse_micro(
            patient_id=f"GDSC::{line}",
            sensitivity_mean=sensitivity_mean,
            vascularity=float(rng.uniform(0.3, 0.6)),
            immune_frac=0.1,  # monoculture cell line: minimal immune infiltrate
            tumor_burden=float(rng.uniform(0.6, 0.85)),
            metadata=metadata,
        )
        records.append(PatientRecord(
            patient_id=f"GDSC::{line}",
            micro=micro,
            observed_response=observed,
            features=features,
        ))
    _log(f"built {len(records)} GDSC glioma cell-line records")
    return records


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def load_real_cohort(cache_dir: str = "data_cache") -> List[PatientRecord]:
    """Parse cached real data into PatientRecords. Returns [] if nothing usable.

    Combines two real open cohorts, each carrying a genuine response label and a
    coarse (synthetic-spatial, real-feature-parameterised) microenvironment:

      * TCGA-GBM patients — real clinical/expression features, survival-based
        response proxy.
      * GDSC glioma cell lines — real drug-screen LN_IC50 response ground truth.

    Never raises.
    """
    records: List[PatientRecord] = []
    try:
        records.extend(_load_gdsc(cache_dir))
    except Exception as exc:  # noqa: BLE001
        _log(f"GDSC loader errored (skipped): {exc}")
    try:
        records.extend(_load_tcga(cache_dir))
    except Exception as exc:  # noqa: BLE001
        _log(f"TCGA loader errored (skipped): {exc}")
    _log(f"total real cohort: {len(records)} records")
    return records


if __name__ == "__main__":
    from .acquire import acquire_all

    status = acquire_all()
    print("\n=== acquire_all status ===")
    print(json.dumps(status, indent=2))

    cohort = load_real_cohort()
    print(f"\n=== load_real_cohort -> {len(cohort)} PatientRecords ===")
    for r in cohort[:5]:
        print(f"  {r.patient_id:28s} y={r.observed_response:.1f} "
              f"n_cells={r.micro.cells.n} "
              f"src={r.micro.metadata.get('source')}")
    if cohort:
        import collections
        srcs = collections.Counter(
            "GDSC" if r.patient_id.startswith("GDSC") else "TCGA" for r in cohort
        )
        ys = [r.observed_response for r in cohort]
        print(f"\n  sources: {dict(srcs)}")
        print(f"  responders: {int(sum(ys))}/{len(ys)}")
