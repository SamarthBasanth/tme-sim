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
    """Pearson r of every row of ``mat`` against ``y`` (fed ranks -> Spearman)."""
    m = np.asarray(mat, dtype=float)
    yc = np.asarray(y, dtype=float) - np.mean(y)
    mc = m - m.mean(axis=1, keepdims=True)
    denom = np.sqrt((mc ** 2).sum(axis=1) * (yc ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
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
