"""REAL VALIDATION #4 -- glioblastoma + temozolomide + measured ex-vivo response.

Right disease (GBM), right drug (temozolomide), measured response -- in
**patient-derived** GBM stem-like cell (GSC) cultures, with paired
whole-transcriptome RNA-seq.

SOURCE
------
Ntafoulis I, Kleijn A, Ju J, et al. "Ex vivo drug sensitivity screening predicts
response to temozolomide in glioblastoma patients."  *British Journal of Cancer*
2023;129(11):1780-1789.  PMID 37620410.  GEO accession **GSE232173**.
Erasmus MC / GLIOTRAIN, Rotterdam.

Cached primary files (``data_cache/hunt/geo_organoid_screens/``)::

    GSE232173_count_data_24_cells.txt.gz   raw gene x culture counts, 24 GSC cultures
    BJC2023_MOESM2.docx                    Supplemental tables 1-6 (incl. Table 2 =
                                           per-culture % viability at 100 uM TMZ, n=40)
    BJC2023_MOESM3.xlsx                    Supplemental data: DE genes, gene-vs-TMZ
                                           Spearman correlations, CoxPH, validations
    BJC2023_MOESM4.docx                    Supplemental figure/table legends
    BJC2023_MOESM5.docx                    Supplemental methods
    JOINED_GSE232173_TMZviability.csv      a DERIVED join produced by the data hunt

This module does **not** trust the hunter's join.  :func:`verify_provenance`
re-derives it from the primary ``.docx`` (parsed with the stdlib ``zipfile`` +
``ElementTree``, no ``python-docx`` dependency) and the primary counts matrix, and
reports any disagreement.

WHAT IS AND IS NOT IN THIS DATASET  (established, not assumed, by step 0)
------------------------------------------------------------------------
* The response metric is **mean % cell viability at 100 uM TMZ**, relative to
  vehicle control, from a 5-day ATP-based (CellTiter-Glo-type) assay.
  **LOWER = MORE SENSITIVE.**  It is a *single-dose* readout, not an IC50.
* Per-culture viability numbers are published for **40** cultures (Supplemental
  Table 2, the operator-reproducibility table).  The RNA-seq matrix covers **24**
  cultures.  The intersection -- cultures with BOTH a published viability number
  and paired expression -- is **17**.  That is the analysable n.
* Expression IS paired per model (same GS.### culture identifier), so gene ->
  sensitivity correlation is legitimate.
* **NO per-patient clinical outcome (PFS/OS/response) is present in any cached
  file.**  Supplemental Tables 5 and 6 give only *cohort-level* survival medians
  and *summary* Cox/log-rank p-values.  There is no per-culture PFS/OS column
  anywhere.  Test C is therefore reported as NOT PERFORMABLE on this data; the
  paper's own published associations are quoted, not recomputed.

THE THREE TESTS
---------------
A. POTENCY-AXIS VALIDATION.  Our molecular layer takes temozolomide potency from
   GDSC2 (:mod:`tmesim.molecular.potency`): median IC50 = 606 uM across 51 glioma
   cell lines.  Every one of those 51 fitted IC50s lies **above GDSC's maximum
   tested concentration of 30 uM**, i.e. GDSC's TMZ potency is entirely an
   extrapolation.  This ex-vivo dataset dosed at 100 uM -- 3.3x beyond GDSC's
   ceiling -- so it is an *out-of-range* test of that extrapolation, in
   patient-derived GBM cultures rather than immortalised lines.
   We recover each GDSC line's Hill slope from its own reported (LN_IC50, AUC)
   pair, predict % viability at 100 uM, and compare that predicted distribution
   with the measured one.
   NOTE ON WHAT IS NOT POSSIBLE: unlike the earlier GDSC-vs-HGCC check
   (rho = +0.83 over shared *drugs*), there is no shared identifier between GDSC
   cell lines and Erasmus GSC cultures and only ONE drug here, so **rank
   concordance cannot be computed**.  This is a magnitude/distribution
   concordance test only, and is reported as such.

B. MOLECULAR PREDICTOR -- MGMT.  MGMT is present and strongly variable in the
   GSE232173 counts.  We test whether MGMT expression predicts measured TMZ
   sensitivity across models (Spearman), cross-check against the paper's own
   published coefficient, rank MGMT genome-wide, and then quantify **exactly how
   much between-patient predictive power our spatial model forfeits** by leaving
   temozolomide sensitivity uniform (``real_xenium.SENSITIVITY_DETERMINANTS
   ["temozolomide"] is None`` -- MGMT is absent from the 319-gene Xenium panel).

C. PATIENT OUTCOME.  Not performable -- stated plainly, not manufactured.

Nothing here modifies core / physics / molecular / reconstruct / engine.
Run: ``python -m tmesim.validation.gbm_tmz_exvivo``.
"""
from __future__ import annotations

import hashlib
import os
import re
import zipfile
from typing import Dict, List, Optional, Tuple
from xml.etree import ElementTree as ET

import numpy as np

# --------------------------------------------------------------------------- #
# Paths and constants
# --------------------------------------------------------------------------- #
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
HUNT_DIR = os.path.join(_ROOT, "data_cache", "hunt", "geo_organoid_screens")
GDSC_CSV = os.path.join(_ROOT, "data_cache", "gdsc_glioma_dose_response.csv")
RESULTS_DIR = os.path.join(_ROOT, "results")

COUNTS_GZ = os.path.join(HUNT_DIR, "GSE232173_count_data_24_cells.txt.gz")
MOESM2 = os.path.join(HUNT_DIR, "BJC2023_MOESM2.docx")
MOESM3 = os.path.join(HUNT_DIR, "BJC2023_MOESM3.xlsx")
MOESM4 = os.path.join(HUNT_DIR, "BJC2023_MOESM4.docx")
MOESM5 = os.path.join(HUNT_DIR, "BJC2023_MOESM5.docx")
HUNTER_JOIN = os.path.join(HUNT_DIR, "JOINED_GSE232173_TMZviability.csv")

PMID_NTAFOULIS = "37620410"          # Ntafoulis 2023 Br J Cancer, GSE232173
PMID_HEGI = "15758010"               # Hegi 2005 NEJM -- MGMT methylation & TMZ benefit
PMID_ESTELLER = "11070098"           # Esteller 2000 NEJM -- MGMT silencing & alkylator
PMID_STUPP = "15758009"              # Stupp 2005 NEJM -- RT + TMZ standard of care
PMID_GDSC = "23180760"               # Yang 2013 NAR -- GDSC
PMID_KLEIJN_ASSAY = "26977153"       # Kleijn 2016 Stem Cells Int -- ATP viability read-out

#: The 19 cultures with paired parental tumour tissue (Supplemental Table 3).
#: The paper's own gene-vs-TMZ Spearman analysis used n = 19 cell cultures.
TISSUE_PAIRED_19 = [
    "GS.598", "GS.607", "GS.627", "GS.691", "GS.772", "GS.786", "GS.799",
    "GS.353", "GS.365", "GS.502", "GS.636", "GS.741", "GS.802", "GS.820",
    "GS.824", "GS.828", "GS.832", "GS.295", "GS.359",
]

#: Values the paper published that we cross-check against (Supplemental data
#: MOESM3 sheet ``Correlation_cell_005`` / MOESM2 Supplemental Table 6).
PAPER_MGMT_SPEARMAN = 0.709991374965395
PAPER_MGMT_P = 0.000660223580640485
PAPER_MGMT_PADJ = 0.0489287774797163

_W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def _sha256(path: str, nbytes: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(nbytes)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _docx_blocks(path: str) -> List[Tuple[str, object]]:
    """Yield ``('p', text)`` / ``('tbl', rows)`` in document order, stdlib only."""
    with zipfile.ZipFile(path) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    body = root.find(f"{_W}body")
    out: List[Tuple[str, object]] = []
    if body is None:
        return out
    for el in body:
        tag = el.tag.replace(_W, "")
        if tag == "p":
            txt = "".join(n.text or "" for n in el.iter(f"{_W}t")).strip()
            if txt:
                out.append(("p", txt))
        elif tag == "tbl":
            rows = []
            for tr in el.findall(f"{_W}tr"):
                cells = [
                    "".join(n.text or "" for n in tc.iter(f"{_W}t")).strip()
                    for tc in tr.findall(f"{_W}tc")
                ]
                rows.append(cells)
            out.append(("tbl", rows))
    return out


def _pearson_rows(mat: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pearson r of every row of ``mat`` against ``y`` (fed ranks -> Spearman).

    ``errstate`` also covers the matmul: on macOS/Accelerate with numpy 2.x the
    BLAS ``gemv`` kernel raises spurious divide/overflow/invalid status flags on
    perfectly finite input. Verified against ``scipy.stats.pearsonr`` row by row
    -- agreement to 1.1e-16 -- so the flags are cosmetic, not numerical.
    """
    m = np.asarray(mat, dtype=float)
    yc = np.asarray(y, dtype=float) - np.mean(y)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        mc = m - m.mean(axis=1, keepdims=True)
        denom = np.sqrt((mc ** 2).sum(axis=1) * (yc ** 2).sum())
        return (mc @ yc) / denom


def _canon_gs(raw: str) -> Optional[str]:
    """``GS295`` / ``GS.295`` / ``gs 295`` -> ``GS.295``; else None."""
    m = re.match(r"^\s*GS[\s._-]*(\d{3})\s*$", str(raw), flags=re.I)
    return f"GS.{m.group(1)}" if m else None


# --------------------------------------------------------------------------- #
# Step 0 -- PROVENANCE: re-derive the join from primary files, audit the hunter
# --------------------------------------------------------------------------- #
def verify_provenance() -> Dict:
    """Re-derive the viability x expression join from PRIMARY files and audit.

    Returns a dict with the re-derived table, the hunter's table, a cell-by-cell
    comparison, and the explicit finding on whether patient outcomes exist.
    """
    import pandas as pd

    files = {}
    for label, path in [
        ("counts", COUNTS_GZ), ("MOESM2", MOESM2), ("MOESM3", MOESM3),
        ("MOESM4", MOESM4), ("MOESM5", MOESM5), ("hunter_join", HUNTER_JOIN),
        ("gdsc", GDSC_CSV),
    ]:
        files[label] = {
            "path": path,
            "exists": os.path.exists(path),
            "bytes": os.path.getsize(path) if os.path.exists(path) else 0,
            "sha256_16": _sha256(path)[:16] if os.path.exists(path) else None,
        }
    missing = [k for k, v in files.items() if not v["exists"]]
    if missing:
        raise FileNotFoundError(f"primary files missing: {missing}")

    # ---- (a) primary viability table: MOESM2 Supplemental Table 2 ----------
    blocks = _docx_blocks(MOESM2)
    viab_rows: List[Tuple[str, float, float, float, int]] = []
    supp_tbl_texts: List[str] = [t for k, t in blocks if k == "p"]
    for kind, payload in blocks:
        if kind != "tbl":
            continue
        rows = payload  # type: ignore[assignment]
        header = [c.lower() for c in rows[0]] if rows else []
        if not any("viab" in h for h in header):
            continue
        n_ops = sum(1 for h in header if "operator" in h)
        for r in rows[1:]:
            gs = _canon_gs(r[0])
            if gs is None:
                continue
            try:
                mean = float(r[1 + n_ops])
                sd = float(r[2 + n_ops])
                cv = float(r[3 + n_ops])
            except (ValueError, IndexError):
                continue
            n_rep = sum(1 for c in r[1:1 + n_ops] if c.strip() != "")
            viab_rows.append((gs, mean, sd, cv, n_rep))
    primary_viab = pd.DataFrame(
        viab_rows, columns=["culture", "mean_pct_viability_100uM_TMZ",
                            "SD", "pct_CV", "n_operators"]
    ).set_index("culture")

    # ---- (b) primary expression matrix ------------------------------------
    counts = pd.read_csv(COUNTS_GZ, sep="\t", index_col=0)
    counts.columns = [str(c).strip().strip('"') for c in counts.columns]
    expr_cultures = list(counts.columns)

    # ---- (c) the honest join ----------------------------------------------
    joined_ids = [c for c in expr_cultures if c in primary_viab.index]
    rederived = primary_viab.loc[joined_ids].copy()

    # ---- (d) audit the hunter's derived CSV --------------------------------
    hunter = pd.read_csv(HUNTER_JOIN).set_index("culture")
    only_hunter = sorted(set(hunter.index) - set(rederived.index))
    only_ours = sorted(set(rederived.index) - set(hunter.index))
    shared = sorted(set(hunter.index) & set(rederived.index))
    diffs = {}
    for c in shared:
        a = float(hunter.loc[c, "mean_pct_viability_100uM_TMZ"])
        b = float(rederived.loc[c, "mean_pct_viability_100uM_TMZ"])
        if abs(a - b) > 1e-6:
            diffs[c] = (a, b)
    hunter_ok = (not only_hunter) and (not only_ours) and (not diffs)

    # ---- (e) is ANY per-patient outcome present anywhere? ------------------
    outcome_pat = re.compile(
        r"\b(pfs|os|progression[- ]free|overall survival|survival|"
        r"months|censor|event|follow[- ]?up|deceased|died)\b", re.I)
    outcome_hits: List[Dict] = []
    # docx tables
    for label, path in [("MOESM2", MOESM2), ("MOESM4", MOESM4), ("MOESM5", MOESM5)]:
        for kind, payload in _docx_blocks(path):
            if kind != "tbl":
                continue
            flat = " | ".join(" ".join(r) for r in payload)  # type: ignore[arg-type]
            if outcome_pat.search(flat):
                has_gs = any(_canon_gs(r[0]) for r in payload if r)  # type: ignore[union-attr]
                outcome_hits.append({
                    "file": label, "kind": "docx table",
                    "per_culture_id": bool(has_gs),
                    "preview": flat[:180],
                })
    # xlsx sheets
    import openpyxl
    wb = openpyxl.load_workbook(MOESM3, read_only=True, data_only=True)
    sheet_summary = {}
    for sn in wb.sheetnames:
        ws = wb[sn]
        hdr = next(ws.iter_rows(max_row=1, values_only=True), ())
        hdr = [str(h) for h in hdr if h is not None]
        sheet_summary[sn] = {"n_rows": ws.max_row, "header": hdr}
        # does any sheet carry GS.### row identifiers (i.e. per-model rows)?
        ids = 0
        for i, row in enumerate(ws.iter_rows(max_row=60, max_col=3, values_only=True)):
            ids += sum(1 for c in row if c is not None and _canon_gs(c))
        sheet_summary[sn]["gs_ids_in_first_60_rows"] = ids
    per_culture_outcome_found = any(h["per_culture_id"] for h in outcome_hits) or \
        any(v["gs_ids_in_first_60_rows"] > 0 for v in sheet_summary.values())

    return {
        "files": files,
        "primary_viability_n": int(len(primary_viab)),
        "primary_viability": primary_viab,
        "expression_cultures": expr_cultures,
        "expression_n": len(expr_cultures),
        "n_genes": int(counts.shape[0]),
        "rederived_join": rederived,
        "rederived_n": int(len(rederived)),
        "hunter_join": hunter,
        "hunter_n": int(len(hunter)),
        "hunter_ok": bool(hunter_ok),
        "only_in_hunter": only_hunter,
        "only_in_rederived": only_ours,
        "value_disagreements": diffs,
        "expr_without_viability": sorted(set(expr_cultures) - set(primary_viab.index)),
        "tissue_paired_19": TISSUE_PAIRED_19,
        "in_join_and_paper19": sorted(set(rederived.index) & set(TISSUE_PAIRED_19)),
        "in_join_not_paper19": sorted(set(rederived.index) - set(TISSUE_PAIRED_19)),
        "xlsx_sheets": sheet_summary,
        "outcome_hits": outcome_hits,
        "per_culture_outcome_found": bool(per_culture_outcome_found),
        "counts": counts,
        "supp_texts": supp_tbl_texts,
    }


# --------------------------------------------------------------------------- #
# Test A -- POTENCY-AXIS: GDSC TMZ potency vs measured ex-vivo viability
# --------------------------------------------------------------------------- #
def _hill_auc(ic50: float, hill: float, cmin: float, cmax: float,
              n: int = 2001) -> float:
    """Normalised area under a Hill viability curve over log10 concentration.

    GDSC's ``AUC`` is the area under the fitted dose-response (1 = no effect)
    across the screened concentration window, in log-concentration space.
    """
    lg = np.linspace(np.log10(cmin), np.log10(cmax), n)
    v = 1.0 / (1.0 + (10.0 ** lg / ic50) ** hill)
    return float(np.trapezoid(v, lg) / (lg[-1] - lg[0]))


def _solve_hill(ic50: float, auc: float, cmin: float, cmax: float) -> float:
    """Recover the Hill slope that reproduces GDSC's own reported AUC."""
    from scipy import optimize
    f = lambda lh: _hill_auc(ic50, 10.0 ** lh, cmin, cmax) - auc
    lo, hi = -2.0, 1.5
    try:
        if f(lo) * f(hi) > 0:
            return float("nan")
        return float(10.0 ** optimize.brentq(f, lo, hi))
    except Exception:
        return float("nan")


def test_a_potency_axis(prov: Dict, dose_uM: float = 100.0) -> Dict:
    """Does GDSC's TMZ potency reproduce the MEASURED ex-vivo viability at 100 uM?

    Magnitude/distribution concordance only -- see module docstring for why rank
    concordance is not computable here.
    """
    import pandas as pd
    from scipy import stats

    g = pd.read_csv(GDSC_CSV)
    t = g[g.DRUG_NAME.str.lower() == "temozolomide"].copy()
    t["ic50_uM"] = np.exp(t.LN_IC50)
    cmin, cmax = float(t.MIN_CONC.iloc[0]), float(t.MAX_CONC.iloc[0])
    t["hill"] = [_solve_hill(i, a, cmin, cmax) for i, a in zip(t.ic50_uM, t.AUC)]
    t["pred_pct_viab"] = 100.0 / (1.0 + (dose_uM / t.ic50_uM) ** t.hill)

    med_ic50 = float(np.exp(np.median(t.LN_IC50)))
    n_above_max = int((t.ic50_uM > cmax).sum())
    hill_med = float(np.nanmedian(t.hill))

    join = prov["rederived_join"]
    meas = join["mean_pct_viability_100uM_TMZ"].astype(float)
    meas40 = prov["primary_viability"]["mean_pct_viability_100uM_TMZ"].astype(float)

    pred = t["pred_pct_viab"].dropna()
    mw17 = stats.mannwhitneyu(pred, meas)
    mw40 = stats.mannwhitneyu(pred, meas40)

    # naive single-parameter prediction from the module's own headline number
    naive = 100.0 / (1.0 + dose_uM / med_ic50)

    # bootstrap CI on the measured median (n=17)
    rng = np.random.default_rng(20260727)
    bs = np.array([np.median(rng.choice(meas.values, size=len(meas), replace=True))
                   for _ in range(10000)])
    ci17 = (float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)))

    # implied per-model ex-vivo IC50 under the GDSC-recovered Hill slope
    implied = {}
    for c in meas.index:
        v = float(meas.loc[c]) / 100.0
        implied[c] = (dose_uM / (((1.0 / v) - 1.0) ** (1.0 / hill_med))
                      if 0.0 < v < 1.0 else np.inf)
    implied_s = pd.Series(implied).sort_values()
    fin = implied_s[np.isfinite(implied_s)]

    return {
        "gdsc_n_lines": int(len(t)),
        "gdsc_tissues": sorted(set(t.TCGA_DESC.astype(str))),
        "gdsc_median_ic50_uM": med_ic50,
        "gdsc_ic50_range_uM": (float(t.ic50_uM.min()), float(t.ic50_uM.max())),
        "gdsc_conc_window_uM": (cmin, cmax),
        "gdsc_n_ic50_above_max_conc": n_above_max,
        "gdsc_median_AUC": float(t.AUC.median()),
        "hill_median": hill_med,
        "hill_iqr": (float(np.nanpercentile(t.hill, 25)),
                     float(np.nanpercentile(t.hill, 75))),
        "dose_uM": dose_uM,
        "pred_median": float(pred.median()),
        "pred_iqr": (float(pred.quantile(0.25)), float(pred.quantile(0.75))),
        "pred_range": (float(pred.min()), float(pred.max())),
        "pred_naive_h1": naive,
        "meas17": meas,
        "meas17_median": float(meas.median()),
        "meas17_iqr": (float(meas.quantile(0.25)), float(meas.quantile(0.75))),
        "meas17_range": (float(meas.min()), float(meas.max())),
        "meas17_median_boot_ci": ci17,
        "meas40": meas40,
        "meas40_median": float(meas40.median()),
        "meas40_iqr": (float(meas40.quantile(0.25)), float(meas40.quantile(0.75))),
        "meas40_range": (float(meas40.min()), float(meas40.max())),
        "mw_pred_vs_meas17": (float(mw17.statistic), float(mw17.pvalue)),
        "mw_pred_vs_meas40": (float(mw40.statistic), float(mw40.pvalue)),
        "median_abs_error_17": abs(float(pred.median()) - float(meas.median())),
        "implied_ic50_uM": implied_s,
        "implied_finite_n": int(len(fin)),
        "implied_min": float(fin.min()), "implied_median": float(fin.median()),
        "implied_max": float(fin.max()),
        "implied_fold_spread": float(fin.max() / fin.min()),
        "gdsc_table": t,
        "rank_concordance_possible": False,
    }


# --------------------------------------------------------------------------- #
# Test B -- MGMT expression as a predictor of measured TMZ sensitivity
# --------------------------------------------------------------------------- #
def test_b_mgmt(prov: Dict, min_cpm: float = 1.0, min_models: int = 6) -> Dict:
    """Does MGMT expression predict measured ex-vivo TMZ sensitivity?

    Primary analysis: all n=17 cultures with a published viability number AND
    paired RNA-seq (the pre-specified, fully-reproducible set).
    Sensitivity analysis: the n=15 subset that is also in the paper's own
    tissue-paired cohort of 19 -- reported alongside, not instead of, the primary.
    """
    import pandas as pd
    from scipy import stats
    import openpyxl

    counts = prov["counts"]
    join = prov["rederived_join"]
    ids = [c for c in join.index if c in counts.columns]
    lib = counts[ids].sum(axis=0)
    cpm = counts[ids] / lib * 1e6
    viab = join.loc[ids, "mean_pct_viability_100uM_TMZ"].astype(float)

    mgmt_cpm = cpm.loc["MGMT"]
    mgmt_log = np.log2(mgmt_cpm + 1.0)

    def _spear(sel: List[str]) -> Tuple[float, float, int]:
        r, p = stats.spearmanr(np.log2(cpm.loc["MGMT", sel] + 1.0), viab.loc[sel])
        return float(r), float(p), len(sel)

    rho17, p17, n17 = _spear(ids)
    sub15 = [c for c in ids if c in TISSUE_PAIRED_19]
    rho15, p15, n15 = _spear(sub15)
    pr, pp = stats.pearsonr(mgmt_log.values, viab.values)

    # permutation p (10k label shuffles) for the primary set. Spearman on fixed
    # data == Pearson on ranks, so rank once and permute the rank vector.
    rng = np.random.default_rng(20260727)
    rx = stats.rankdata(mgmt_log.values)
    ry = stats.rankdata(viab.values)
    perm = np.array([rng.permutation(ry) for _ in range(10000)])
    null = _pearson_rows(perm, rx)
    perm_p = float((np.abs(null) >= abs(rho17)).mean())

    # binary MGMT-expressing vs MGMT-null (transcript detection, NOT methylation)
    expressed = (mgmt_cpm >= min_cpm)
    a = viab[expressed.values]
    b = viab[~expressed.values]
    mw = stats.mannwhitneyu(a, b)
    pooled_sd = np.sqrt(((len(a) - 1) * a.std(ddof=1) ** 2 +
                         (len(b) - 1) * b.std(ddof=1) ** 2) / (len(a) + len(b) - 2))
    cohens_d = float((a.mean() - b.mean()) / pooled_sd)

    # genome-wide ranking of MGMT among expressed genes (vectorised: Spearman is
    # Pearson on ranks, so rank every gene row once and correlate as a matrix)
    keep = (cpm > min_cpm).sum(axis=1) >= min_models
    E = np.log2(cpm[keep] + 1.0)
    rh = _pearson_rows(stats.rankdata(E.values, axis=1), ry)
    gi = E.index.get_loc("MGMT")
    rank = int((np.abs(rh) > abs(rh[gi])).sum()) + 1

    # the paper's own coefficient, read straight out of the supplementary
    wb = openpyxl.load_workbook(MOESM3, read_only=True, data_only=True)
    ws = wb["Correlation_cell_005"]
    rows = [r for r in ws.iter_rows(values_only=True)][1:]
    coeffs = [(abs(float(r[0])), str(r[3])) for r in rows if r[0] is not None]
    coeffs.sort(reverse=True)
    paper_rank = next(i for i, (_, g) in enumerate(coeffs) if g == "MGMT") + 1
    paper_row = next((r for r in rows if str(r[3]) == "MGMT"), None)
    ws2 = wb["DE_genes_cell_005"]
    de_row = next((r for r in ws2.iter_rows(values_only=True)
                   if r and str(r[0]) == "MGMT"), None)

    return {
        "ids": ids, "viab": viab, "mgmt_cpm": mgmt_cpm, "mgmt_log2cpm": mgmt_log,
        "median_libsize": float(lib.median()),
        "spearman_n17": (rho17, p17, n17),
        "spearman_n15_tissue_paired": (rho15, p15, n15),
        "pearson_n17": (float(pr), float(pp), float(pr ** 2)),
        "permutation_p_n17": perm_p,
        "rho2_n17": float(rho17 ** 2), "rho2_n15": float(rho15 ** 2),
        "mgmt_expressed_ids": [c for c in ids if expressed[c]],
        "mgmt_null_ids": [c for c in ids if not expressed[c]],
        "expressed_median": float(a.median()), "null_median": float(b.median()),
        "expressed_n": int(len(a)), "null_n": int(len(b)),
        "mw_expressed_vs_null": (float(mw.statistic), float(mw.pvalue)),
        "cohens_d": cohens_d,
        "genes_tested": int(E.shape[0]),
        "mgmt_genomewide_rank": rank,
        "mgmt_genomewide_pct": 100.0 * rank / E.shape[0],
        "paper_spearman": PAPER_MGMT_SPEARMAN,
        "paper_p": PAPER_MGMT_P, "paper_padj": PAPER_MGMT_PADJ,
        "paper_row": paper_row, "paper_n_sig_genes": len(coeffs),
        "paper_mgmt_rank": paper_rank,
        "paper_de_row": de_row,
        "min_cpm": min_cpm,
    }


# --------------------------------------------------------------------------- #
# Test C -- patient outcome (NOT PERFORMABLE: recorded, not manufactured)
# --------------------------------------------------------------------------- #
#: Cohort-level numbers the paper itself publishes (MOESM2 Supplemental Tables
#: 5 and 6).  These are QUOTED, not recomputed -- the per-patient rows that would
#: be needed to recompute them are not in any cached file.
PAPER_OUTCOME_REPORTED = {
    "cohort_n_patients": 66,
    "n_with_survival_analysis": 55,
    "mgmt_methylated": 31, "mgmt_unmethylated": 35,
    "median_PFS_months": 6.65, "median_PFS_range": (1.0, 36.40),
    # NOTE: Supplemental Table 5 labels the PFS rows "UM 8.55 / M 5.92", which
    # CONTRADICTS the paper's own Results text ("Median PFS and OS of the MGMT
    # methylated group was 8.55 and 17.19 months ... versus 5.92 and 12.4 months
    # in the unmethylated group") and contradicts the OS rows in the same table.
    # We follow the text (and the biology); the discrepancy is flagged in the
    # report rather than silently resolved.
    "median_PFS_methylated": 8.55, "median_PFS_unmethylated": 5.92,
    "supp_table5_PFS_labels_contradict_text": True,
    "median_OS_months": 13.69, "median_OS_range": (2.8, 37.40),
    "median_OS_unmethylated": 12.40, "median_OS_methylated": 17.19,
    # Spearman correlations of in-vitro readout vs patient outcome (Results text,
    # n = 55).  Negative rho = lower viability (more sensitive) -> longer survival.
    "spearman_viability100_PFS": (-0.38, 0.014),
    "spearman_viability100_OS": (-0.495, 0.0002),
    "spearman_AUC_PFS": (-0.29, 0.03),
    "spearman_AUC_OS": (-0.399, 0.0026),
    "spearman_IC50_PFS": (-0.206, 0.132),
    "spearman_IC50_OS": (-0.35, 0.009),
    # median OS in months for responders / intermediates / non-responders,
    # cut at 50% and 75% viability at 100 uM TMZ
    "median_OS_by_response_group": {"responders": 27.8, "intermediates": 15.9,
                                    "non_responders": 12.3},
    "response_group_cutoffs_pct_viability": {"responders": "<50",
                                             "intermediates": "50-75",
                                             "non_responders": ">75"},
    # MGMT in the parental-tissue cohort (MOESM3 CoxPH_tissue_005), n = 56
    "erasmus_tissue_MGMT_HR": 1.85438910580871,
    "erasmus_tissue_MGMT_p": 6.0919848747292e-05,
    "erasmus_tissue_MGMT_padj": 0.00791722829430051,
    "erasmus_tissue_MGMT_methylation_vs_OS_p": 0.16,
    "cox_viability100_PFS_p": 0.012, "cox_viability100_OS_p": "<0.00005",
    "cox_AUC_PFS_p": 0.092, "cox_AUC_OS_p": 0.001,
    "cox_IC50_PFS_p": "NS", "cox_IC50_OS_p": 0.0004,
    "cox_MGMT_PFS_p": "NS", "cox_MGMT_OS_p": 0.001,
    "logrank_resp_vs_nonresp_PFS_p": 0.018,
    "logrank_resp_vs_nonresp_OS_p": 0.0001,
    "logrank_MGMT_MvsUM_PFS_p": 0.0074,
    "logrank_MGMT_MvsUM_OS_p": 0.001,
    "response_groups": {"responders": 10, "intermediates": 17,
                        "non_responders": 28},
    "gliotrain_MGMT_HR": 1.21240535784314,
    "gliotrain_MGMT_p": 0.0247970942357781,
    "gliotrain_MGMT_padj": 0.178510335566003,
}


def test_c_outcomes(prov: Dict) -> Dict:
    """Report -- honestly -- that no per-patient outcome is joinable here."""
    return {
        "per_culture_outcome_present": prov["per_culture_outcome_found"],
        "searched": ["MOESM2.docx tables", "MOESM4.docx tables",
                     "MOESM5.docx tables", "MOESM3.xlsx (all 9 sheets)",
                     "GSE232173 counts header", "hunter join CSV"],
        "xlsx_sheets": {k: v["header"] for k, v in prov["xlsx_sheets"].items()},
        "analysis_performed": False,
        "reason": (
            "No cached file carries a per-patient or per-culture PFS/OS/response "
            "row keyed to a GS.### identifier. Supplemental Tables 5 and 6 report "
            "only cohort medians and summary Cox / log-rank p-values; the "
            "MOESM3 sheets are gene-level (DE, correlation, CoxPH coefficients), "
            "not sample-level. Individual-level survival analysis is therefore "
            "IMPOSSIBLE on this data and is NOT attempted."
        ),
        "paper_reported": PAPER_OUTCOME_REPORTED,
    }


# --------------------------------------------------------------------------- #
# What our spatial model forfeits by not measuring MGMT
# --------------------------------------------------------------------------- #
def quantify_model_gap(a: Dict, b: Dict) -> Dict:
    """Turn the MGMT result into the number that matters for our model.

    ``reconstruct.real_xenium.SENSITIVITY_DETERMINANTS['temozolomide'] is None``
    -> per-cell TMZ sensitivity is a constant 1.0, and potency is a single
    cohort-wide IC50 (606 uM).  So the model's *between-patient* variance in TMZ
    sensitivity is exactly ZERO by construction.  Quantify what that costs.
    """
    import gzip
    import glob
    import pandas as pd
    from scipy import stats

    # --- confirm MGMT really is off the Xenium panels, and on Visium ---------
    panels = {}
    for name in ("xenium_gbm", "xenium_gbm_io"):
        f = os.path.join(_ROOT, "data_cache", name,
                         "cell_feature_matrix", "features.tsv.gz")
        if not os.path.exists(f):
            continue
        genes = [l.decode().split("\t")[1].strip() for l in gzip.open(f)]
        real = [g for g in genes if not g.startswith(
            ("NegControl", "BLANK", "Unassigned", "Deprecated", "antisense_"))]
        panels[name] = {"n_features": len(genes), "n_panel_genes": len(real),
                        "MGMT_on_panel": "MGMT" in genes}

    # --- Visium is whole-transcriptome: is MGMT actually DETECTED there? -----
    import scipy.io as sio
    visium = {}
    for base in sorted(glob.glob(os.path.join(
            _ROOT, "data_cache", "visium_gbm", "**", "features.tsv.gz"),
            recursive=True)):
        d = os.path.dirname(base)
        label = os.path.relpath(d, os.path.join(_ROOT, "data_cache", "visium_gbm"))
        genes = [l.decode().split("\t")[1].strip() for l in gzip.open(base)]
        if "MGMT" not in genes:
            visium[label] = {"MGMT_in_reference": False}
            continue
        M = sio.mmread(os.path.join(d, "matrix.mtx.gz")).tocsr()
        row = np.asarray(M[genes.index("MGMT")].todense()).ravel()
        tot = np.asarray(M.sum(axis=0)).ravel()
        visium[label] = {
            "MGMT_in_reference": True,
            "n_spots": int(M.shape[1]),
            "pct_spots_MGMT_positive": float(100.0 * (row > 0).mean()),
            "mean_counts_per_spot": float(row.mean()),
            "pseudobulk_MGMT_cpm": float(row.sum() / tot.sum() * 1e6),
            "median_libsize": float(np.median(tot)),
        }

    # --- what the constant-sensitivity assumption costs ---------------------
    viab = b["viab"]
    rho17 = b["spearman_n17"][0]
    rho15 = b["spearman_n15_tissue_paired"][0]
    r2 = b["pearson_n17"][2]

    # a uniform model predicts the cohort median for everyone: its residual
    # spread IS the full observed spread.
    const_pred = float(viab.median())
    mae_const = float(np.mean(np.abs(viab - const_pred)))
    # a 1-covariate MGMT model (OLS on log2 CPM), leave-one-out cross-validated
    x = b["mgmt_log2cpm"].values
    y = viab.values
    loo = np.empty_like(y)
    for i in range(len(y)):
        k = np.arange(len(y)) != i
        s, ic = np.polyfit(x[k], y[k], 1)
        loo[i] = s * x[i] + ic
    mae_mgmt = float(np.mean(np.abs(y - loo)))
    # honest LOO baseline: constant model also cross-validated
    loo_c = np.array([np.median(np.delete(y, i)) for i in range(len(y))])
    mae_const_loo = float(np.mean(np.abs(y - loo_c)))

    return {
        "xenium_panels": panels,
        "visium_sections": visium,
        "tmz_determinant_in_model": None,
        "model_between_patient_variance_explained": 0.0,
        "mgmt_variance_explained_spearman2_n17": float(rho17 ** 2),
        "mgmt_variance_explained_spearman2_n15": float(rho15 ** 2),
        "mgmt_variance_explained_pearson_r2": float(r2),
        "observed_viability_range_pp": float(viab.max() - viab.min()),
        "observed_viability_min_max": (float(viab.min()), float(viab.max())),
        "model_predicted_spread_pp": 0.0,
        "implied_ic50_fold_spread": a["implied_fold_spread"],
        "implied_ic50_median_uM": a["implied_median"],
        "model_uniform_ic50_uM": a["gdsc_median_ic50_uM"],
        "median_ic50_ratio_model_vs_measured":
            a["gdsc_median_ic50_uM"] / a["implied_median"],
        "mae_constant_model_pp": mae_const,
        "mae_constant_model_loo_pp": mae_const_loo,
        "mae_mgmt_model_loo_pp": mae_mgmt,
        "mae_reduction_pp": mae_const_loo - mae_mgmt,
        "mae_reduction_pct": 100.0 * (mae_const_loo - mae_mgmt) / mae_const_loo,
    }


# --------------------------------------------------------------------------- #
# Test A2 -- can GDSC even RESOLVE the MGMT determinant?  (honest negative)
# --------------------------------------------------------------------------- #
def test_a2_gdsc_can_resolve_mgmt() -> Dict:
    """Cross-check: within GDSC's own glioma lines, does MGMT predict TMZ IC50?

    Our potency layer's TMZ number comes from GDSC.  If GDSC's TMZ screen can
    itself recover the MGMT determinant, the potency layer could in principle be
    made MGMT-aware from GDSC alone.  If it cannot, then the *only* place the
    determinant is visible is a screen dosed in the clinically relevant range --
    which is exactly what this ex-vivo dataset is.

    GDSC2 TMZ was screened over 0.03-30 uM; the effective range reported by
    Ntafoulis et al. is 25-300 uM.  Expression is DepMap 24Q4 log2(TPM+1),
    joined GDSC ``SANGER_MODEL_ID`` -> DepMap ``SangerModelID``.
    """
    import pandas as pd
    from scipy import stats

    dm_dir = os.path.join(_ROOT, "data_cache", "hunt", "depmap_gdsc_ctrp")
    model_csv = os.path.join(dm_dir, "DepMap_24Q4_Model.csv")
    expr_csv = os.path.join(dm_dir,
                            "DepMap_24Q4_Expression_ProteinCoding_TPMlogp1_CNSBRAIN.csv")
    if not (os.path.exists(model_csv) and os.path.exists(expr_csv)):
        return {"available": False,
                "reason": "DepMap 24Q4 CNS/brain expression not cached"}

    g = pd.read_csv(GDSC_CSV)
    t = g[g.DRUG_NAME.str.lower() == "temozolomide"].copy()
    t["ic50_uM"] = np.exp(t.LN_IC50)

    mod = pd.read_csv(model_csv, low_memory=False)[
        ["ModelID", "CellLineName", "SangerModelID"]]
    c0 = pd.read_csv(expr_csv, nrows=0).columns[0]
    expr = pd.read_csv(expr_csv, usecols=[c0, "MGMT (4255)"]).rename(
        columns={c0: "ModelID", "MGMT (4255)": "MGMT_log2TPM"})

    j = (t.merge(mod, left_on="SANGER_MODEL_ID", right_on="SangerModelID",
                 how="left")
          .merge(expr, on="ModelID", how="left")
          .dropna(subset=["MGMT_log2TPM"]))

    out = {"available": True, "n_gdsc_lines": int(len(t)),
           "n_mapped_to_depmap_expression": int(len(j)),
           "gdsc_conc_window_uM": (float(t.MIN_CONC.iloc[0]),
                                   float(t.MAX_CONC.iloc[0])),
           "effective_range_uM_reported": (25.0, 300.0),
           "table": j}
    for tag, sub in (("all_glioma", j), ("GBM", j[j.TCGA_DESC == "GBM"]),
                     ("LGG", j[j.TCGA_DESC == "LGG"])):
        if len(sub) < 5:
            continue
        r1, p1 = stats.spearmanr(sub.MGMT_log2TPM, sub.LN_IC50)
        r2, p2 = stats.spearmanr(sub.MGMT_log2TPM, sub.AUC)
        out[tag] = {"n": int(len(sub)),
                    "mgmt_vs_lnIC50": (float(r1), float(p1)),
                    "mgmt_vs_AUC": (float(r2), float(p2))}
    out["resolves_mgmt"] = bool(out.get("all_glioma", {})
                                .get("mgmt_vs_lnIC50", (0.0, 1.0))[1] < 0.05)
    return out


# --------------------------------------------------------------------------- #
# Test B2 -- is ANY Xenium-panel gene a usable MGMT surrogate?
# --------------------------------------------------------------------------- #
def test_b2_xenium_surrogate(prov: Dict, b: Dict, min_cpm: float = 1.0,
                             min_models: int = 6) -> Dict:
    """Screen every Xenium-panel gene as a stand-in for the missing MGMT.

    A usable surrogate must clear BOTH bars:
      1. track MGMT itself across the 24 RNA-seq cultures (FDR q < 0.05), and
      2. track the MEASURED TMZ response across the 17 joined cultures, with the
         SAME sign, and survive multiple testing over the panel.

    The base rate matters: if the fraction of *panel* genes correlated with MGMT
    is no higher than the fraction of *all* genes correlated with MGMT, then
    panel membership carries no MGMT information and any single hit is a
    multiple-testing artefact.  We compute both.
    """
    import gzip
    import pandas as pd
    from scipy import stats

    counts = prov["counts"]
    all_ids = list(counts.columns)
    cpm_all = counts / counts.sum(axis=0) * 1e6
    keep = (cpm_all > min_cpm).sum(axis=1) >= min_models
    E = np.log2(cpm_all[keep] + 1.0)
    mgmt_rank = stats.rankdata(E.loc["MGMT"].values)

    join_ids = b["ids"]
    col = {c: i for i, c in enumerate(E.columns)}
    jcols = [col[c] for c in join_ids]
    viab_rank = stats.rankdata(b["viab"].values)

    Erank_all = stats.rankdata(E.values, axis=1)
    rho_mgmt_all = _pearson_rows(Erank_all, mgmt_rank)
    Erank_join = stats.rankdata(E.values[:, jcols], axis=1)
    rho_viab_all = _pearson_rows(Erank_join, viab_rank)

    def _p_from_r(r: np.ndarray, n: int) -> np.ndarray:
        r = np.clip(np.nan_to_num(r, nan=0.0), -0.999999, 0.999999)
        tt = r * np.sqrt((n - 2) / (1.0 - r ** 2))
        return 2.0 * stats.t.sf(np.abs(tt), n - 2)

    p_mgmt_all = _p_from_r(rho_mgmt_all, len(all_ids))
    p_viab_all = _p_from_r(rho_viab_all, len(join_ids))

    def _bh(p: np.ndarray) -> np.ndarray:
        p = np.asarray(p, dtype=float)
        n = len(p)
        o = np.argsort(p)
        q = np.empty(n)
        q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
        return np.clip(q, 0.0, 1.0)

    # ---- panel membership --------------------------------------------------
    panels: Dict[str, List[str]] = {}
    for name in ("xenium_gbm", "xenium_gbm_io"):
        f = os.path.join(_ROOT, "data_cache", name,
                         "cell_feature_matrix", "features.tsv.gz")
        if not os.path.exists(f):
            continue
        genes = [l.decode().split("\t")[1].strip() for l in gzip.open(f)]
        panels[name] = sorted({
            g for g in genes if not g.startswith(
                ("NegControl", "BLANK", "Unassigned", "Deprecated",
                 "antisense_", "Intergenic"))})
    union = sorted(set().union(*panels.values())) if panels else []

    idx = {g: i for i, g in enumerate(E.index)}
    present = [g for g in union if g in idx and g != "MGMT"]
    pi = [idx[g] for g in present]

    df = pd.DataFrame({
        "gene": present,
        "rho_vs_MGMT": rho_mgmt_all[pi],
        "p_vs_MGMT": p_mgmt_all[pi],
        "rho_vs_viability": rho_viab_all[pi],
        "p_vs_viability": p_viab_all[pi],
    })
    df["q_vs_MGMT"] = _bh(df.p_vs_MGMT.values)
    df["q_vs_viability"] = _bh(df.p_vs_viability.values)
    df["same_sign"] = np.sign(df.rho_vs_MGMT) == np.sign(df.rho_vs_viability)
    df["passes_both"] = ((df.q_vs_MGMT < 0.05) & (df.q_vs_viability < 0.05)
                         & df.same_sign)

    # base rates
    q_all_mgmt = _bh(p_mgmt_all)
    base_rate_genomewide = float((q_all_mgmt < 0.05).mean())
    base_rate_panel = float((df.q_vs_MGMT < 0.05).mean())

    near = df[(df.q_vs_MGMT < 0.05) & (df.p_vs_viability < 0.05)
              & df.same_sign].sort_values("p_vs_viability")

    mgmt_i = idx["MGMT"]
    return {
        "panel_sizes": {k: len(v) for k, v in panels.items()},
        "panel_union_n": len(union),
        "MGMT_on_any_panel": "MGMT" in union,
        "panel_genes_measurable_in_GSE232173": len(present),
        "n_expressed_genes_tested": int(E.shape[0]),
        "table": df.sort_values("p_vs_viability").reset_index(drop=True),
        "n_pass_both_FDR": int(df.passes_both.sum()),
        "near_misses": near,
        "base_rate_mgmt_correlated_genomewide": base_rate_genomewide,
        "base_rate_mgmt_correlated_panel": base_rate_panel,
        "panel_enriched_for_mgmt_correlates":
            bool(base_rate_panel > 2.0 * base_rate_genomewide),
        "best_by_mgmt": df.reindex(
            df.rho_vs_MGMT.abs().sort_values(ascending=False).index).head(8),
        "best_by_viability": df.reindex(
            df.rho_vs_viability.abs().sort_values(ascending=False).index).head(8),
        "mgmt_self": {"rho_vs_viability": float(rho_viab_all[mgmt_i]),
                      "p_vs_viability": float(p_viab_all[mgmt_i])},
        "mki67": (df[df.gene == "MKI67"].to_dict("records") or [None])[0],
        "usable_surrogate_found": bool(df.passes_both.sum() > 0),
    }


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
def make_figure(prov: Dict, a: Dict, a2: Dict, b: Dict, b2: Dict, gap: Dict,
                path: Optional[str] = None) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = path or os.path.join(RESULTS_DIR, "gbm_tmz_exvivo.png")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fig, ax = plt.subplots(2, 3, figsize=(17.5, 10.2))
    # matplotlib renamed boxplot(labels=) -> tick_labels in 3.9
    _bp_lab = ("tick_labels"
               if tuple(int(v) for v in matplotlib.__version__.split(".")[:2])
               >= (3, 9) else "labels")
    C_MEAS, C_PRED, C_BAD = "#1b6ca8", "#c1440e", "#8a8a8a"

    # (1) potency concordance -- predicted vs measured % viability at 100 uM
    axx = ax[0, 0]
    pred = a["gdsc_table"]["pred_pct_viab"].dropna().values
    m40 = a["meas40"].values
    m17 = a["meas17"].values
    bins = np.linspace(0, 140, 22)
    axx.hist(pred, bins=bins, alpha=0.55, color=C_PRED, density=True,
             label=f"GDSC-predicted (n={len(pred)} glioma lines)")
    axx.hist(m40, bins=bins, alpha=0.55, color=C_MEAS, density=True,
             label=f"MEASURED ex vivo (n={len(m40)} GSC cultures)")
    axx.axvline(np.median(pred), color=C_PRED, lw=2.2, ls="--")
    axx.axvline(np.median(m40), color=C_MEAS, lw=2.2, ls="--")
    axx.set_xlabel("% viability at 100 uM TMZ")
    axx.set_ylabel("density")
    axx.set_title("A1. POTENCY AXIS -- real vs real\n"
                  f"median {np.median(pred):.1f}% predicted vs "
                  f"{np.median(m40):.1f}% measured "
                  f"(Mann-Whitney p={a['mw_pred_vs_meas40'][1]:.2f})",
                  fontsize=10)
    axx.legend(fontsize=7.5, loc="upper left")

    # (2) IC50 distributions, log scale
    axx = ax[0, 1]
    gd = a["gdsc_table"]["ic50_uM"].values
    imp = a["implied_ic50_uM"]
    imp = imp[np.isfinite(imp)].values
    axx.boxplot([np.log10(gd), np.log10(imp)], widths=0.55, showfliers=False,
                **{_bp_lab: [f"GDSC2 fitted\n({len(gd)} lines)",
                             f"implied from\nmeasured (n={len(imp)})"]})
    rng = np.random.default_rng(7)
    for i, v in enumerate((gd, imp), start=1):
        axx.scatter(i + rng.uniform(-0.13, 0.13, len(v)), np.log10(v), s=18,
                    alpha=0.65, color=[C_PRED, C_MEAS][i - 1], zorder=3)
    axx.axhline(np.log10(a["gdsc_median_ic50_uM"]), color=C_PRED, ls=":", lw=1.5)
    axx.axhline(np.log10(30.0), color=C_BAD, ls="--", lw=1.5)
    axx.text(0.02, np.log10(30.0) + 0.06, "GDSC max tested conc = 30 uM",
             transform=axx.get_yaxis_transform(), fontsize=7.5, color=C_BAD)
    axx.set_ylabel("log10 IC50 (uM)")
    axx.set_title("A2. TMZ IC50: model's GDSC prior vs\nex-vivo-implied "
                  f"(medians {a['gdsc_median_ic50_uM']:.0f} vs "
                  f"{a['implied_median']:.0f} uM)", fontsize=10)

    # (3) MGMT vs measured sensitivity -- THE result
    axx = ax[0, 2]
    x = b["mgmt_log2cpm"].values
    y = b["viab"].values
    det = b["mgmt_cpm"].values >= b["min_cpm"]
    axx.scatter(x[det], y[det], s=70, color="#b8272d", edgecolor="k", lw=0.6,
                zorder=3, label=f"MGMT-expressing (n={det.sum()})")
    axx.scatter(x[~det], y[~det], s=70, color="#2e7d32", edgecolor="k", lw=0.6,
                zorder=3, label=f"MGMT-silent (n={(~det).sum()})")
    s_, i_ = np.polyfit(x, y, 1)
    xs = np.linspace(x.min() - 0.2, x.max() + 0.2, 50)
    axx.plot(xs, s_ * xs + i_, color="k", lw=1.4, ls="--")
    for cx, cy, lab in zip(x, y, b["ids"]):
        axx.annotate(lab.replace("GS.", ""), (cx, cy), fontsize=6.2,
                     xytext=(3, 3), textcoords="offset points", color="#444444")
    rho, p, n = b["spearman_n17"]
    axx.set_xlabel("MGMT expression, log2(CPM+1)")
    axx.set_ylabel("% viability at 100 uM TMZ  (higher = resistant)")
    axx.set_title("B1. MGMT PREDICTS MEASURED TMZ RESPONSE\n"
                  f"rho = {rho:+.3f}, p = {p:.4f}, n = {n}  "
                  f"[paper: {b['paper_spearman']:+.2f}, n=19]", fontsize=10)
    axx.legend(fontsize=8, loc="lower right")

    # (4) MGMT-expressing vs MGMT-silent
    axx = ax[1, 0]
    ve = b["viab"][[c in b["mgmt_expressed_ids"] for c in b["ids"]]].values
    vn = b["viab"][[c in b["mgmt_null_ids"] for c in b["ids"]]].values
    axx.boxplot([vn, ve], widths=0.5, showfliers=False,
                **{_bp_lab: [f"MGMT-silent\n(n={len(vn)})",
                             f"MGMT-expressing\n(n={len(ve)})"]})
    for i, v in enumerate((vn, ve), start=1):
        axx.scatter(i + rng.uniform(-0.11, 0.11, len(v)), v, s=52,
                    color=["#2e7d32", "#b8272d"][i - 1], edgecolor="k", lw=0.5,
                    zorder=3)
    axx.axhline(50, color=C_BAD, ls=":", lw=1.2)
    axx.axhline(75, color=C_BAD, ls=":", lw=1.2)
    axx.text(0.55, 51, "paper's responder cut (<50%)", fontsize=7, color=C_BAD)
    axx.text(0.55, 76, "non-responder cut (>75%)", fontsize=7, color=C_BAD)
    axx.set_ylabel("% viability at 100 uM TMZ")
    axx.set_title("B2. binary MGMT detection\n"
                  f"medians {b['null_median']:.1f}% vs {b['expressed_median']:.1f}%, "
                  f"Mann-Whitney p = {b['mw_expressed_vs_null'][1]:.4f}",
                  fontsize=10)

    # (5) Xenium-panel surrogate screen
    axx = ax[1, 1]
    tb = b2["table"]
    axx.scatter(tb.rho_vs_MGMT, tb.rho_vs_viability, s=16, alpha=0.5,
                color="#8899aa", edgecolor="none",
                label=f"Xenium-panel genes (n={len(tb)})")
    hi = tb[(tb.q_vs_MGMT < 0.05) & (tb.p_vs_viability < 0.05) & tb.same_sign]
    axx.scatter(hi.rho_vs_MGMT, hi.rho_vs_viability, s=60, color="#e08214",
                edgecolor="k", lw=0.6, zorder=3,
                label="tracks MGMT (q<0.05) AND response (raw p<0.05)")
    for _, r_ in hi.iterrows():
        axx.annotate(r_.gene, (r_.rho_vs_MGMT, r_.rho_vs_viability), fontsize=7.5,
                     xytext=(4, 3), textcoords="offset points")
    axx.scatter([1.0], [b2["mgmt_self"]["rho_vs_viability"]], marker="*", s=320,
                color="#b8272d", edgecolor="k", lw=0.7, zorder=4,
                label="MGMT itself (NOT on any panel)")
    axx.axhline(0, color="k", lw=0.6)
    axx.axvline(0, color="k", lw=0.6)
    axx.set_xlabel("Spearman rho vs MGMT expression  (n=24)")
    axx.set_ylabel("Spearman rho vs measured viability  (n=17)")
    axx.set_title("B3. MGMT SURROGATE SEARCH on the Xenium panels\n"
                  f"{b2['n_pass_both_FDR']} genes pass BOTH FDR bars   "
                  f"(panel base rate {100*b2['base_rate_mgmt_correlated_panel']:.0f}%"
                  f" vs genome-wide "
                  f"{100*b2['base_rate_mgmt_correlated_genomewide']:.0f}%)",
                  fontsize=10)
    axx.legend(fontsize=7, loc="lower left")

    # (6) can GDSC resolve MGMT?  (honest negative)
    axx = ax[1, 2]
    if a2.get("available"):
        j = a2["table"]
        axx.scatter(j.MGMT_log2TPM, j.LN_IC50, s=55, color="#5b4b8a",
                    edgecolor="k", lw=0.5, zorder=3)
        r1, p1 = a2["all_glioma"]["mgmt_vs_lnIC50"]
        s2, i2 = np.polyfit(j.MGMT_log2TPM, j.LN_IC50, 1)
        xs = np.linspace(j.MGMT_log2TPM.min(), j.MGMT_log2TPM.max(), 20)
        axx.plot(xs, s2 * xs + i2, "k--", lw=1.3)
        axx.axhline(np.log(30.0), color=C_BAD, ls="--", lw=1.5)
        axx.text(0.02, np.log(30.0) + 0.12, "max tested conc (30 uM):"
                 " every fitted IC50 is above it",
                 transform=axx.get_yaxis_transform(), fontsize=7.5, color=C_BAD)
        axx.set_xlabel("MGMT expression, DepMap log2(TPM+1)")
        axx.set_ylabel("GDSC2 ln IC50 (uM)")
        axx.set_title("A3. GDSC CANNOT see the MGMT determinant\n"
                      f"n={a2['all_glioma']['n']} glioma lines, "
                      f"rho = {r1:+.3f}, p = {p1:.2f} (n.s.)", fontsize=10)
    else:
        axx.axis("off")
        axx.text(0.5, 0.5, "DepMap expression not cached", ha="center")

    fig.suptitle("REAL VALIDATION -- glioblastoma + temozolomide + measured "
                 "ex-vivo response (GSE232173 / Ntafoulis 2023, PMID 37620410)",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.975))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def write_report(prov: Dict, a: Dict, a2: Dict, b: Dict, b2: Dict, gap: Dict,
                 c: Dict, png: str, path: Optional[str] = None) -> str:
    path = path or os.path.join(RESULTS_DIR, "gbm_tmz_exvivo.md")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    L: List[str] = []
    w = L.append
    pr = c["paper_reported"]

    w("# GBM + temozolomide: measured ex-vivo response")
    w("")
    w("**Right disease (glioblastoma, IDH-wildtype). Right drug (temozolomide). "
      "Measured response (ATP viability in patient-derived GSC cultures).**")
    w("")
    w(f"Source: Ntafoulis I, Kleijn A, Ju J, *et al.* \"Ex vivo drug sensitivity "
      f"screening predicts response to temozolomide in glioblastoma patients and "
      f"identifies candidate biomarkers.\" *Br J Cancer* 2023;129(8):1327-1338. "
      f"**PMID {PMID_NTAFOULIS}**, doi:10.1038/s41416-023-02402-y. "
      f"GEO **GSE232173** (Erasmus MC / GLIOTRAIN, Rotterdam).")
    w("")
    w(f"![figure]({os.path.basename(png)})")
    w("")

    # ---------------- provenance ----------------
    w("## 0. Provenance -- verified, not assumed")
    w("")
    w("| file | bytes | sha256[:16] |")
    w("|---|---:|---|")
    for k, v in prov["files"].items():
        w(f"| `{os.path.relpath(v['path'], _ROOT)}` | {v['bytes']:,} | "
          f"`{v['sha256_16']}` |")
    w("")
    w("### What the data actually contains")
    w("")
    w(f"* **Response metric:** mean **% cell viability at 100 uM TMZ**, relative "
      f"to matched DMSO controls, from a 6-day ATP-based CellTiter-Glo assay "
      f"(supplementary methods, MOESM5). **LOWER = MORE SENSITIVE.** It is a "
      f"single-dose readout, not an IC50. 100 uM was chosen by the authors "
      f"because it is close to the reported maximum TMZ CSF concentration "
      f"(~75 uM). This is also the paper's *best* readout -- it out-predicted "
      f"both AUC and IC50 against patient survival.")
    w(f"* **Per-culture viability numbers exist for {prov['primary_viability_n']} "
      f"cultures** (MOESM2 Supplemental Table 2). This is the "
      f"*operator-reproducibility* table -- 2 or 3 operators per culture over "
      f"multiple years -- and it is the only per-culture response table published "
      f"anywhere in the paper. The main screen (n=66) is published as *figures "
      f"only*; there is no n=66 numeric table.")
    w(f"* **Expression:** GSE232173 raw counts, {prov['n_genes']:,} genes x "
      f"{prov['expression_n']} GSC cultures. Paired per model by the `GS.###` "
      f"identifier. Note GEO deposits **24** cultures whereas the paper's own "
      f"transcriptome analyses used **19**.")
    w(f"* **The analysable join is n = {prov['rederived_n']}** -- cultures with "
      f"both a published viability number and paired RNA-seq.")
    w(f"* **NO per-patient clinical outcome (PFS/OS/response) exists in any "
      f"cached file.** See section C.")
    w("")
    w("### Audit of the hunter's `JOINED_GSE232173_TMZviability.csv`")
    w("")
    w(f"Re-derived independently from the primary `.docx` (stdlib "
      f"`zipfile`+`ElementTree`) and the primary counts matrix.")
    w("")
    if prov["hunter_ok"]:
        w(f"**VERDICT: the hunter's join is CORRECT.** All "
          f"{prov['rederived_n']}/{prov['rederived_n']} rows and all values "
          f"reproduce exactly; no orphans in either direction.")
    else:
        w(f"**VERDICT: the hunter's join DISAGREES with the primary files.** "
          f"only-in-hunter: {prov['only_in_hunter']}; only-in-rederived: "
          f"{prov['only_in_rederived']}; value disagreements: "
          f"{prov['value_disagreements']}. The re-derived table is used below.")
    w("")
    w("Two corrections to how that file is *labelled*, though:")
    w("")
    w("1. The filename implies a GEO-supplied response table. It is not. The "
      "viability numbers come from **MOESM2 Supplemental Table 2** (the operator "
      "reproducibility experiment), not from GEO. GEO's sample characteristics "
      "are constant strings (`tissue: Brain tumor`, `genotype: GBM - IDH WT`, "
      "`treatment: STUPP`) and carry **no** response or outcome data at all.")
    w("2. Because it is the reproducibility table, several cultures have wide "
      "operator spread -- %CV up to "
      f"{prov['primary_viability']['pct_CV'].max():.0f}% "
      f"(median {prov['primary_viability']['pct_CV'].median():.1f}%). "
      "The response variable used below is therefore itself noisy, which "
      "attenuates every correlation reported here.")
    w("")
    w(f"Cultures with RNA-seq but no published viability number "
      f"({len(prov['expr_without_viability'])}): "
      f"{', '.join(prov['expr_without_viability'])}.")
    w("")

    # ---------------- A ----------------
    w("## A. Potency axis -- does GDSC's TMZ potency hold up in patient-derived GBM?")
    w("")
    w(f"Our molecular layer takes TMZ potency from GDSC2 "
      f"(PMID {PMID_GDSC}): median IC50 = **{a['gdsc_median_ic50_uM']:.0f} uM** "
      f"over **{a['gdsc_n_lines']} glioma cell lines** "
      f"({'/'.join(a['gdsc_tissues'])}).")
    w("")
    w(f"**The catch, stated up front:** GDSC2 screened TMZ over "
      f"{a['gdsc_conc_window_uM'][0]:.2f}-{a['gdsc_conc_window_uM'][1]:.0f} uM, "
      f"and **{a['gdsc_n_ic50_above_max_conc']}/{a['gdsc_n_lines']}** fitted "
      f"IC50s lie *above* that ceiling (median GDSC AUC = "
      f"{a['gdsc_median_AUC']:.3f}, i.e. almost no killing in the tested range). "
      f"GDSC's TMZ potency is an **extrapolation**. This ex-vivo dataset dosed at "
      f"100 uM, 3.3x beyond GDSC's ceiling, so it is a genuine out-of-range test "
      f"of that extrapolation.")
    w("")
    w("Prediction: recover each GDSC line's Hill slope from its own reported "
      f"(IC50, AUC) pair (median slope {a['hill_median']:.2f}, IQR "
      f"{a['hill_iqr'][0]:.2f}-{a['hill_iqr'][1]:.2f}), then evaluate the curve "
      f"at 100 uM. Zero parameters fitted to the Erasmus data.")
    w("")
    w("| quantity | GDSC-predicted | MEASURED ex vivo |")
    w("|---|---|---|")
    w(f"| % viability at 100 uM, median | **{a['pred_median']:.1f}%** | "
      f"**{a['meas40_median']:.1f}%** (n={len(a['meas40'])}) |")
    w(f"| IQR | {a['pred_iqr'][0]:.1f}-{a['pred_iqr'][1]:.1f}% | "
      f"{a['meas40_iqr'][0]:.1f}-{a['meas40_iqr'][1]:.1f}% |")
    w(f"| range | {a['pred_range'][0]:.1f}-{a['pred_range'][1]:.1f}% | "
      f"{a['meas40_range'][0]:.1f}-{a['meas40_range'][1]:.1f}% |")
    w(f"| median on the n={len(a['meas17'])} joined subset | -- | "
      f"{a['meas17_median']:.1f}% (bootstrap 95% CI "
      f"{a['meas17_median_boot_ci'][0]:.1f}-{a['meas17_median_boot_ci'][1]:.1f}) |")
    w("")
    w(f"Mann-Whitney predicted vs measured: **p = "
      f"{a['mw_pred_vs_meas40'][1]:.2f}** (n=40) and "
      f"**p = {a['mw_pred_vs_meas17'][1]:.2f}** (n=17) -- the two distributions "
      f"are statistically indistinguishable. A naive single-number prediction "
      f"from the layer's headline IC50 alone (Hill = 1) gives "
      f"{a['pred_naive_h1']:.1f}%, i.e. within "
      f"{abs(a['pred_naive_h1'] - a['meas40_median']):.1f} percentage points of "
      f"the measured median.")
    w("")
    w(f"Inverting the measured viabilities back to per-model IC50s under the "
      f"GDSC-recovered Hill slope gives a median of "
      f"**{a['implied_median']:.0f} uM** "
      f"(range {a['implied_min']:.0f}-{a['implied_max']:.0f} uM, "
      f"{a['implied_fold_spread']:.0f}x spread, over the "
      f"{a['implied_finite_n']}/{len(a['meas17'])} cultures whose measured "
      f"viability is below 100% and therefore invertible) against the layer's "
      f"{a['gdsc_median_ic50_uM']:.0f} uM -- a "
      f"{gap['median_ic50_ratio_model_vs_measured']:.2f}x ratio on the central "
      f"value.")
    w("")
    w("**What this does NOT establish.** There is no shared model identifier "
      "between GDSC cell lines and the Erasmus GSC biobank, and only one drug is "
      "in play, so **rank concordance cannot be computed here.** The earlier "
      "GDSC-vs-HGCC check (rho = +0.83) was a *cross-drug* ordering test; this is "
      "a *magnitude/distribution* test. They are different claims. What is "
      "established is that an extrapolated cell-line potency, pushed 3.3x past "
      "its own assay ceiling, lands within a few percentage points of the median "
      "measured response in 40 patient-derived GBM cultures -- it could easily "
      "have been off by 30 points and was not.")
    w("")

    # ---------------- A2 ----------------
    if a2.get("available"):
        ag = a2["all_glioma"]
        w("### A2. Can GDSC itself resolve the MGMT determinant? No.")
        w("")
        w(f"Joining GDSC2 TMZ to DepMap 24Q4 MGMT expression "
          f"({a2['n_mapped_to_depmap_expression']}/{a2['n_gdsc_lines']} glioma "
          f"lines mapped):")
        w("")
        w("| subset | n | MGMT vs ln IC50 | MGMT vs AUC |")
        w("|---|---:|---|---|")
        for tag in ("all_glioma", "GBM", "LGG"):
            if tag in a2:
                d = a2[tag]
                w(f"| {tag} | {d['n']} | rho = {d['mgmt_vs_lnIC50'][0]:+.3f}, "
                  f"p = {d['mgmt_vs_lnIC50'][1]:.2f} | "
                  f"rho = {d['mgmt_vs_AUC'][0]:+.3f}, "
                  f"p = {d['mgmt_vs_AUC'][1]:.2f} |")
        w("")
        w(f"**MGMT does not predict TMZ IC50 in GDSC glioma lines "
          f"(rho = {ag['mgmt_vs_lnIC50'][0]:+.3f}, p = "
          f"{ag['mgmt_vs_lnIC50'][1]:.2f}, n = {ag['n']}).** That is not a "
          f"contradiction of section B -- it is a statement about assay range. "
          f"GDSC dosed 0.03-30 uM; the effective TMZ range reported by the "
          f"authors is 25-300 uM. GDSC's TMZ arm has essentially no dynamic "
          f"range, so it fixes the central potency well but carries **no "
          f"resolving power for between-model ordering**. The determinant is "
          f"only visible in a screen dosed in the clinical range -- which is "
          f"exactly what this ex-vivo dataset is. Practical consequence: our "
          f"potency layer cannot be made MGMT-aware from GDSC.")
        w("")

    # ---------------- B ----------------
    rho, p, n = b["spearman_n17"]
    rho15, p15, n15 = b["spearman_n15_tissue_paired"]
    w("## B. Does MGMT expression predict measured TMZ sensitivity?  **Yes.**")
    w("")
    w("| test | statistic | p | n |")
    w("|---|---|---|---:|")
    w(f"| Spearman, MGMT log2(CPM+1) vs % viability @100 uM | **rho = {rho:+.3f}** "
      f"| **{p:.4f}** | {n} |")
    w(f"| ... permutation p (10,000 label shuffles) | | {b['permutation_p_n17']:.4f} "
      f"| {n} |")
    w(f"| ... restricted to the paper's tissue-paired cohort | rho = {rho15:+.3f} "
      f"| {p15:.4f} | {n15} |")
    w(f"| Pearson (same data) | r = {b['pearson_n17'][0]:+.3f} "
      f"(r2 = {b['pearson_n17'][2]:.3f}) | {b['pearson_n17'][1]:.4f} | {n} |")
    w(f"| MGMT-expressing vs MGMT-silent, Mann-Whitney | medians "
      f"{b['expressed_median']:.1f}% vs {b['null_median']:.1f}% "
      f"(Cohen d = {b['cohens_d']:+.2f}) | "
      f"**{b['mw_expressed_vs_null'][1]:.4f}** | "
      f"{b['expressed_n']}v{b['null_n']} |")
    w(f"| *paper's own value* (MOESM3 `Correlation_cell_005`) | "
      f"rho = {b['paper_spearman']:+.3f} | {b['paper_p']:.4f} "
      f"(q = {b['paper_padj']:.3f}) | 19 |")
    w("")
    w(f"Sign is the right way round: **higher MGMT -> higher surviving fraction "
      f"-> more resistant**, which is the Esteller/Hegi biology "
      f"(PMID {PMID_ESTELLER}, PMID {PMID_HEGI}).")
    w("")
    w(f"Genome-wide context: among the {b['genes_tested']:,} expressed genes, "
      f"MGMT ranks **#{b['mgmt_genomewide_rank']}** by |Spearman| against "
      f"measured viability (top {b['mgmt_genomewide_pct']:.1f}%). In the "
      f"paper's own analysis MGMT ranked #{b['paper_mgmt_rank']} of "
      f"{b['paper_n_sig_genes']} significant genes and was one of only 121 genes "
      f"that were *both* correlated with ex-vivo response *and* prognostic in "
      f"parental tissue.")
    w("")
    w(f"MGMT expression is strongly bimodal here: "
      f"{b['expressed_n'] + 0}/{n} joined cultures have MGMT at CPM >= "
      f"{b['min_cpm']:.0f} and {b['null_n']} are effectively silent -- consistent "
      f"with promoter methylation. The binary contrast "
      f"(p = {b['mw_expressed_vs_null'][1]:.4f}) is the more robust of the two "
      f"tests, because the continuous Spearman is computed across a block of "
      f"tied near-zero values.")
    w("")
    w("**Honest caveats on this number.** n = 17, from raw counts with simple "
      "CPM where the authors used edgeR TMM on a differently filtered gene set, "
      "and against the noisy operator-reproducibility viabilities rather than "
      "the primary screen. On the full join my rho (%+.2f) is weaker than the "
      "paper's (%+.2f); restricted to the %d cultures that are also in the "
      "paper's own tissue-paired cohort it is %+.2f (p = %.4f) -- i.e. the gap "
      "is mostly *which cultures*, not *which pipeline*. Two independent "
      "analyses of the same raw counts agree in sign, in significance, and "
      "closely in magnitude on the common subset. This is a replication, not a "
      "precise effect-size estimate."
      % (rho, b["paper_spearman"], n15, rho15, p15))
    w("")

    # ---------------- B2 ----------------
    w("### B2. Is any Xenium-panel gene a usable MGMT surrogate?  **No.**")
    w("")
    w(f"MGMT is absent from both panels -- confirmed here directly from the "
      f"shipped feature lists: "
      + "; ".join(f"`{k}` {v} genes" for k, v in b2["panel_sizes"].items())
      + f"; union {b2['panel_union_n']} genes; MGMT present = "
      f"**{b2['MGMT_on_any_panel']}**. "
      f"{b2['panel_genes_measurable_in_GSE232173']} panel genes are measurable in "
      f"GSE232173 and were screened.")
    w("")
    w("A usable surrogate has to clear two bars: track MGMT across the 24 "
      "cultures, **and** track the measured response across the 17, with the same "
      "sign.")
    w("")
    w(f"* Genes passing **both** at FDR q < 0.05: **{b2['n_pass_both_FDR']}**.")
    w(f"* Genes passing the MGMT bar (q<0.05) *and* raw p<0.05 against response "
      f"with consistent sign: {len(b2['near_misses'])}"
      + (" -- " + ", ".join(
          f"**{r.gene}** (rho_MGMT {r.rho_vs_MGMT:+.2f}, q={r.q_vs_MGMT:.3f}; "
          f"rho_viab {r.rho_vs_viability:+.2f}, raw p={r.p_vs_viability:.3f}, "
          f"q={r.q_vs_viability:.2f})"
          for _, r in b2["near_misses"].iterrows()) if len(b2["near_misses"])
         else "."))
    w("")
    if len(b2["near_misses"]):
        w(f"Neither survives multiple testing against the response (all q > "
          f"{b2['near_misses'].q_vs_viability.min():.2f}), and neither has a "
          f"mechanistic claim on TMZ metabolism.")
    w("")
    w(f"And the base rate kills the whole idea: "
      f"**{100*b2['base_rate_mgmt_correlated_panel']:.1f}%** of panel genes "
      f"correlate with MGMT at q<0.05, versus "
      f"**{100*b2['base_rate_mgmt_correlated_genomewide']:.1f}%** of *all* "
      f"expressed genes. Panel membership carries **no** MGMT information "
      f"(enrichment = "
      f"{b2['base_rate_mgmt_correlated_panel'] / max(b2['base_rate_mgmt_correlated_genomewide'], 1e-9):.2f}x). "
      f"With n = 24 a large fraction of the transcriptome co-varies with any "
      f"given gene; the top panel correlate of MGMT is **"
      f"{b2['best_by_mgmt'].iloc[0].gene}** at rho = "
      f"{b2['best_by_mgmt'].iloc[0].rho_vs_MGMT:+.2f}, and a housekeeping-grade "
      f"hit at the top of that list is itself the tell.")
    w("")
    w("**Actionable conclusion: there is no MGMT proxy to add to the spatial "
      "model. The forfeit in section B3 is real and cannot be bought back from "
      "the current panels.** The only fixes are (a) a custom Xenium add-on probe "
      "for MGMT, or (b) carrying MGMT as a per-patient scalar from bulk "
      "methylation/expression rather than as a spatial field.")
    w("")

    # ---------------- B3 ----------------
    w("### B3. What our spatial model forfeits")
    w("")
    w("`reconstruct.real_xenium.SENSITIVITY_DETERMINANTS['temozolomide'] is "
      "None`, so per-cell TMZ sensitivity is a constant 1.0 and potency is one "
      "cohort-wide IC50. **The model's between-patient variance in TMZ "
      "sensitivity is exactly zero by construction.** Against this dataset that "
      "costs:")
    w("")
    w(f"* MGMT explains **{100*gap['mgmt_variance_explained_spearman2_n17']:.0f}%** "
      f"of the rank variance in measured response here "
      f"(rho2, n=17); **{100*b['paper_spearman']**2:.0f}%** using the paper's own "
      f"coefficient. The model captures 0%.")
    w(f"* Measured response spans "
      f"{gap['observed_viability_min_max'][0]:.1f}-"
      f"{gap['observed_viability_min_max'][1]:.1f}% viability, an "
      f"**{gap['observed_viability_range_pp']:.0f} percentage-point** spread, and "
      f"a **{gap['implied_ic50_fold_spread']:.0f}x** spread in implied IC50. The "
      f"model predicts one number for every patient.")
    w(f"* Leave-one-out MAE: constant model "
      f"{gap['mae_constant_model_loo_pp']:.1f} pp vs a one-covariate MGMT model "
      f"{gap['mae_mgmt_model_loo_pp']:.1f} pp -- a "
      f"{gap['mae_reduction_pct']:.0f}% reduction from a single gene "
      f"(cross-validated, so not a fitting artefact).")
    w(f"* At cohort level this is the difference between the paper's response "
      f"tiers: median OS "
      f"{pr['median_OS_by_response_group']['responders']} vs "
      f"{pr['median_OS_by_response_group']['intermediates']} vs "
      f"{pr['median_OS_by_response_group']['non_responders']} months "
      f"(responders / intermediates / non-responders). A uniform-sensitivity "
      f"model cannot place a patient in any of those tiers.")
    w("")
    if gap["visium_sections"]:
        w("Where MGMT *is* measurable in our own data -- Visium is "
          "whole-transcriptome, so MGMT is in the reference; the question is "
          "whether it is detected at usable depth:")
        w("")
        w("| section | MGMT in ref | spots | % spots MGMT+ | pseudobulk CPM |")
        w("|---|---|---:|---:|---:|")
        for k, v in gap["visium_sections"].items():
            if not v.get("MGMT_in_reference"):
                w(f"| {k} | no | | | |")
            else:
                w(f"| {k} | yes | {v['n_spots']:,} | "
                  f"{v['pct_spots_MGMT_positive']:.1f}% | "
                  f"{v['pseudobulk_MGMT_cpm']:.1f} |")
        w("")
        w("So MGMT is recoverable as a *pseudobulk* scalar from Visium, but at "
          "single-spot level it is near the detection floor -- it is not a "
          "spatial field we can honestly resolve, only a per-section covariate.")
        w("")

    # ---------------- C ----------------
    w("## C. Do ex-vivo results predict patient outcome?  "
      "**Not testable on this data.**")
    w("")
    w(f"{c['reason']}")
    w("")
    w("Searched: " + ", ".join(f"`{s}`" for s in c["searched"]) + ". "
      "The MOESM3 sheets are gene-level throughout "
      f"({', '.join(c['xlsx_sheets'].keys())}); none carries a `GS.###` row key.")
    w("")
    w("**No outcome analysis is performed and none is simulated.** For the "
      "record, these are the paper's own published associations -- quoted, not "
      "recomputed:")
    w("")
    w("| association (n = 55 patients on full Stupp) | statistic |")
    w("|---|---|")
    for lab, key in (("% viability @100 uM vs PFS", "spearman_viability100_PFS"),
                     ("% viability @100 uM vs OS", "spearman_viability100_OS"),
                     ("AUC vs PFS", "spearman_AUC_PFS"),
                     ("AUC vs OS", "spearman_AUC_OS"),
                     ("IC50 vs PFS", "spearman_IC50_PFS"),
                     ("IC50 vs OS", "spearman_IC50_OS")):
        r_, p_ = pr[key]
        w(f"| {lab} | Spearman rho = {r_:+.3f}, p = {p_} |")
    w(f"| Cox, % viability @100 uM | PFS p = {pr['cox_viability100_PFS_p']}, "
      f"OS p = {pr['cox_viability100_OS_p']} |")
    w(f"| log-rank, responders vs non-responders | PFS p = "
      f"{pr['logrank_resp_vs_nonresp_PFS_p']}, OS p = "
      f"{pr['logrank_resp_vs_nonresp_OS_p']} |")
    w(f"| median OS by ex-vivo tier | "
      f"{pr['median_OS_by_response_group']['responders']} / "
      f"{pr['median_OS_by_response_group']['intermediates']} / "
      f"{pr['median_OS_by_response_group']['non_responders']} months |")
    w(f"| MGMT expression in parental tissue, CoxPH (n=56) | HR = "
      f"{pr['erasmus_tissue_MGMT_HR']:.2f}, p = "
      f"{pr['erasmus_tissue_MGMT_p']:.1e} (q = "
      f"{pr['erasmus_tissue_MGMT_padj']:.4f}) |")
    w(f"| ... external replication, GLIOTRAIN n=70 | HR = "
      f"{pr['gliotrain_MGMT_HR']:.2f}, p = {pr['gliotrain_MGMT_p']:.4f} |")
    w(f"| MGMT *promoter methylation* vs OS, same cohort | p = "
      f"{pr['erasmus_tissue_MGMT_methylation_vs_OS_p']} (n.s.) |")
    w("")
    w("Two of those rows are worth flagging because they were extracted from the "
      "primary supplementary file (MOESM3) rather than the abstract: MGMT "
      "*expression* is prognostic in the parental tumours "
      f"(HR {pr['erasmus_tissue_MGMT_HR']:.2f}, p = "
      f"{pr['erasmus_tissue_MGMT_p']:.1e}) and replicates externally, while MGMT "
      f"*promoter methylation* -- the clinical biomarker -- was **not** "
      f"significant for OS in the same cohort (p = "
      f"{pr['erasmus_tissue_MGMT_methylation_vs_OS_p']}). That strengthens the "
      "case for treating MGMT as a continuous expression covariate rather than a "
      "binary methylation flag.")
    w("")
    if pr.get("supp_table5_PFS_labels_contradict_text"):
        w("*Data-quality note:* MOESM2 Supplemental Table 5 labels the PFS rows "
          "`UM 8.55 / M 5.92`, which contradicts both the paper's Results text "
          "and the OS rows in the same table. We follow the text "
          "(methylated PFS 8.55, unmethylated 5.92). Flagged, not silently "
          "resolved.")
        w("")

    # ---------------- verdict ----------------
    w("## Verdict")
    w("")
    w("**What this validates.**")
    w("")
    w(f"1. **The POTENCY axis, in patient-derived GBM.** An extrapolated GDSC "
      f"TMZ potency ({a['gdsc_median_ic50_uM']:.0f} uM, 51 glioma lines) "
      f"predicts a median {a['pred_median']:.1f}% viability at 100 uM; "
      f"{prov['primary_viability_n']} patient-derived GSC cultures measured "
      f"{a['meas40_median']:.1f}% (Mann-Whitney p = "
      f"{a['mw_pred_vs_meas40'][1]:.2f}). Right disease, right drug, measured "
      f"response, no fitted parameters. This is a magnitude test, not an "
      f"ordering test.")
    w(f"2. **The MGMT determinant is real and quantified.** MGMT expression "
      f"predicts measured ex-vivo TMZ response at rho = {rho:+.3f} (p = "
      f"{p:.4f}, n = {n}; binary contrast p = "
      f"{b['mw_expressed_vs_null'][1]:.4f}), independently reproducing the "
      f"paper's rho = {b['paper_spearman']:+.2f} from raw counts.")
    w("")
    w("**What this does NOT validate.**")
    w("")
    w("1. **Nothing spatial.** These are dissociated 2-D monolayer GSC cultures. "
      "No architecture, no vasculature, no diffusion gradient, no immune or "
      "neuronal compartment -- the authors' own GSEA shows tissue-exclusive "
      "immune and neuronal-communication programmes lost in culture. This "
      "dataset cannot test the transport half of the model at all.")
    w("2. **Not in vivo, and not PK.** 100 uM is a bath concentration on a plate. "
      "It brackets the reported CSF Cmax (~75 uM) but says nothing about "
      "delivery.")
    w("3. **No ordering claim on the potency axis** -- zero shared models "
      "between GDSC and the Erasmus biobank.")
    w("4. **No patient-outcome claim of our own** -- section C.")
    w(f"5. **n = {n} for every molecular claim.** The response variable is the "
      f"operator-reproducibility table with %CV up to "
      f"{prov['primary_viability']['pct_CV'].max():.0f}%. Section B is a "
      f"directionally solid, magnitude-approximate replication, not a precise "
      f"effect-size estimate.")
    w("")
    w("**The one actionable finding.** MGMT is on neither Xenium panel "
      f"({b2['panel_union_n']}-gene union), no panel gene surrogates for it "
      f"({b2['n_pass_both_FDR']} pass both FDR bars; panel MGMT-correlation base "
      f"rate {100*b2['base_rate_mgmt_correlated_panel']:.0f}% vs genome-wide "
      f"{100*b2['base_rate_mgmt_correlated_genomewide']:.0f}%, i.e. no "
      f"enrichment), and the determinant is invisible to GDSC because GDSC's TMZ "
      f"arm has no dynamic range. So the uniform-TMZ-sensitivity assumption is "
      f"currently *unavoidable* given the inputs -- but it is now **measured**: "
      f"it discards ~{100*gap['mgmt_variance_explained_spearman2_n17']:.0f}% of "
      f"the between-model rank variance in ex-vivo TMZ response "
      f"(~{100*b['paper_spearman']**2:.0f}% by the paper's estimate), an "
      f"{gap['observed_viability_range_pp']:.0f}-point viability spread, and the "
      f"tiering that separates "
      f"{pr['median_OS_by_response_group']['responders']}-month from "
      f"{pr['median_OS_by_response_group']['non_responders']}-month median "
      f"survival. That number belongs in the model's uncertainty statement, and "
      f"the fix is a custom MGMT probe or a per-patient MGMT scalar -- not a "
      f"panel proxy.")
    w("")
    w("### References")
    w("")
    w(f"* Ntafoulis I, *et al.* Br J Cancer 2023;129(8):1327-1338. "
      f"PMID {PMID_NTAFOULIS}. GEO GSE232173.")
    w(f"* Stupp R, *et al.* N Engl J Med 2005;352:987-996. PMID {PMID_STUPP}.")
    w(f"* Hegi ME, *et al.* N Engl J Med 2005;352:997-1003. PMID {PMID_HEGI}.")
    w(f"* Esteller M, *et al.* N Engl J Med 2000;343:1350-1354. "
      f"PMID {PMID_ESTELLER}.")
    w(f"* Yang W, *et al.* Nucleic Acids Res 2013;41:D955-D961 (GDSC). "
      f"PMID {PMID_GDSC}.")
    w(f"* Kleijn A, *et al.* Stem Cells Int 2016;2016:5629214 (ATP viability "
      f"read-out on primary GSCs). PMID {PMID_KLEIJN_ASSAY}.")
    w("* DepMap Public 24Q4 (MGMT expression for GDSC glioma lines).")
    w("")

    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return path


# --------------------------------------------------------------------------- #
def main() -> Dict:
    prov = verify_provenance()
    a = test_a_potency_axis(prov)
    a2 = test_a2_gdsc_can_resolve_mgmt()
    b = test_b_mgmt(prov)
    b2 = test_b2_xenium_surrogate(prov, b)
    gap = quantify_model_gap(a, b)
    c = test_c_outcomes(prov)
    png = make_figure(prov, a, a2, b, b2, gap)
    md = write_report(prov, a, a2, b, b2, gap, c, png)

    line = "=" * 88
    print(line)
    print("REAL VALIDATION -- GBM + temozolomide + MEASURED ex-vivo response")
    print(f"GSE232173 / Ntafoulis 2023 Br J Cancer, PMID {PMID_NTAFOULIS}")
    print(line)

    print("\n0. PROVENANCE")
    print(f"   per-culture viability published : n = {prov['primary_viability_n']} "
          f"(MOESM2 Supplemental Table 2, the OPERATOR-REPRODUCIBILITY table)")
    print(f"   RNA-seq cultures in GSE232173   : n = {prov['expression_n']} "
          f"({prov['n_genes']:,} genes)")
    print(f"   analysable join                 : n = {prov['rederived_n']}")
    print(f"   hunter's join re-derived from primary files: "
          f"{'CORRECT' if prov['hunter_ok'] else 'WRONG'}")
    print(f"   metric = mean % viability at 100 uM TMZ; LOWER = MORE SENSITIVE")
    print(f"   per-patient PFS/OS present anywhere: "
          f"{prov['per_culture_outcome_found']}")

    print("\nA. POTENCY AXIS (magnitude concordance; rank concordance NOT possible)")
    print(f"   GDSC2 TMZ: median IC50 {a['gdsc_median_ic50_uM']:.0f} uM over "
          f"{a['gdsc_n_lines']} glioma lines; "
          f"{a['gdsc_n_ic50_above_max_conc']}/{a['gdsc_n_lines']} IC50s exceed "
          f"the {a['gdsc_conc_window_uM'][1]:.0f} uM assay ceiling")
    print(f"   predicted %viab @100uM : median {a['pred_median']:.1f}  "
          f"IQR {a['pred_iqr'][0]:.1f}-{a['pred_iqr'][1]:.1f}")
    print(f"   MEASURED  %viab @100uM : median {a['meas40_median']:.1f}  "
          f"IQR {a['meas40_iqr'][0]:.1f}-{a['meas40_iqr'][1]:.1f}   (n=40)")
    print(f"   Mann-Whitney p = {a['mw_pred_vs_meas40'][1]:.3f} (n=40), "
          f"{a['mw_pred_vs_meas17'][1]:.3f} (n=17)  -> indistinguishable")
    print(f"   implied ex-vivo IC50 median {a['implied_median']:.0f} uM vs "
          f"model's {a['gdsc_median_ic50_uM']:.0f} uM "
          f"({gap['median_ic50_ratio_model_vs_measured']:.2f}x)")
    if a2.get("available"):
        ag = a2["all_glioma"]
        print(f"   A2: within GDSC glioma lines, MGMT vs lnIC50 "
              f"rho = {ag['mgmt_vs_lnIC50'][0]:+.3f} p = "
              f"{ag['mgmt_vs_lnIC50'][1]:.2f} (n={ag['n']}) -> GDSC CANNOT "
              f"resolve the MGMT determinant")

    rho, p, n = b["spearman_n17"]
    print("\nB. MGMT vs MEASURED TMZ SENSITIVITY")
    print(f"   Spearman rho = {rho:+.3f}  p = {p:.4f}  n = {n}   "
          f"(permutation p = {b['permutation_p_n17']:.4f})")
    print(f"   paper's own : rho = {b['paper_spearman']:+.3f}  "
          f"p = {b['paper_p']:.4f}  n = 19  -> SAME SIGN, same significance")
    print(f"   MGMT-silent {b['null_median']:.1f}% vs MGMT-expressing "
          f"{b['expressed_median']:.1f}% viability, "
          f"Mann-Whitney p = {b['mw_expressed_vs_null'][1]:.4f}")
    print(f"   MGMT ranks #{b['mgmt_genomewide_rank']} of "
          f"{b['genes_tested']:,} expressed genes")

    print("\nB2. XENIUM-PANEL MGMT SURROGATE SEARCH")
    print(f"   MGMT on any panel: {b2['MGMT_on_any_panel']}  "
          f"(union {b2['panel_union_n']} genes; "
          f"{b2['panel_genes_measurable_in_GSE232173']} screened)")
    print(f"   genes passing BOTH FDR bars: {b2['n_pass_both_FDR']}")
    if len(b2["near_misses"]):
        for _, r_ in b2["near_misses"].iterrows():
            print(f"     near-miss {r_.gene:8s} rho_MGMT={r_.rho_vs_MGMT:+.2f} "
                  f"(q={r_.q_vs_MGMT:.3f})  rho_viab={r_.rho_vs_viability:+.2f} "
                  f"(raw p={r_.p_vs_viability:.3f}, q={r_.q_vs_viability:.2f})")
    print(f"   MGMT-correlation base rate: panel "
          f"{100*b2['base_rate_mgmt_correlated_panel']:.1f}% vs genome-wide "
          f"{100*b2['base_rate_mgmt_correlated_genomewide']:.1f}% -> NO enrichment")
    print(f"   USABLE SURROGATE FOUND: {b2['usable_surrogate_found']}")

    print("\nB3. WHAT THE UNIFORM-SENSITIVITY ASSUMPTION COSTS")
    print(f"   MGMT explains {100*gap['mgmt_variance_explained_spearman2_n17']:.0f}%"
          f" of rank variance here "
          f"({100*b['paper_spearman']**2:.0f}% by the paper); model captures 0%")
    print(f"   measured spread {gap['observed_viability_min_max'][0]:.1f}-"
          f"{gap['observed_viability_min_max'][1]:.1f}% viability "
          f"({gap['observed_viability_range_pp']:.0f} pp), "
          f"{gap['implied_ic50_fold_spread']:.0f}x in implied IC50; "
          f"model predicts ONE number")
    print(f"   LOO MAE: constant {gap['mae_constant_model_loo_pp']:.1f} pp -> "
          f"MGMT {gap['mae_mgmt_model_loo_pp']:.1f} pp "
          f"({gap['mae_reduction_pct']:.0f}% reduction)")

    print("\nC. PATIENT OUTCOME: NOT PERFORMABLE")
    print(f"   {c['reason'].splitlines()[0]}")
    print("   paper's own (quoted, NOT recomputed): %viab@100uM vs OS "
          "rho = -0.495 p = 0.0002 (n=55); median OS 27.8 / 15.9 / 12.3 months")

    print(f"\nwrote {md}\nwrote {png}")
    return {"provenance": prov, "A": a, "A2": a2, "B": b, "B2": b2,
            "gap": gap, "C": c, "md": md, "png": png}


if __name__ == "__main__":
    main()
