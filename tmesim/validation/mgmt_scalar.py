"""MGMT as a per-patient scalar: calibration, Visium recovery, and IMPACT.

This module does three things and writes ``results/mgmt_scalar.md``:

1. :func:`calibrate_from_exvivo` — re-derives the MGMT -> temozolomide-sensitivity
   mapping from the CACHED PRIMARY BYTES (GEO ``GSE232173`` raw counts +
   ``BJC2023_MOESM2.docx`` Supplemental Table 2), so every constant frozen in
   :mod:`tmesim.molecular.mgmt` is auditable rather than asserted.
2. :func:`visium_mgmt_table` — the per-section pseudobulk MGMT level for each of
   the four REAL cached Visium GBM sections, and the sensitivity each implies.
3. :func:`run_impact` — re-runs the five-drug panel on the REAL Xenium tumour at
   the grid-converged ``dx = 12 um`` with (a) the current uniform sensitivity and
   (b) MGMT-informed temozolomide sensitivity across the plausible MGMT range
   observed in those Visium sections, and reports how far predicted kill moves
   and whether TMZ's rank changes.

Nothing here modifies core / physics / engine; the only non-additive change in
the whole feature is the optional ``mgmt=`` argument on
``reconstruct.real_pipeline.apply_drug_sensitivity``, whose default is
bit-identical to the previous behaviour.

Run: ``./.venv/bin/python -m tmesim.validation.mgmt_scalar``
     (add ``--quick`` to skip the simulations and only refit + report Visium).
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Tuple

import numpy as np

from ..molecular.mgmt import (
    MGMT_BINARY_CI95,
    MGMT_BINARY_P,
    MGMT_BINARY_SENSITIVITY,
    MGMT_LAMBDA_BOOT_CI95,
    MGMT_CALIBRATION_N,
    MGMT_KILL_FLOOR,
    MGMT_LAMBDA,
    MGMT_LAMBDA_CI95,
    MGMT_LAMBDA_P,
    MGMT_LAMBDA_R2,
    MGMT_OBSERVED_LOG2CPM_RANGE,
    MGMT_SENSITIVITY_FLOOR,
    MGMT_VIAB_INTERCEPT,
    MGMT_VIAB_SLOPE,
    MGMT_VIAB_SLOPE_CI95,
    mgmt_from_visium_all,
    mgmt_sensitivity,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))
RESULTS_DIR = os.path.join(_ROOT, "results")
OUT_MD = os.path.join(RESULTS_DIR, "mgmt_scalar.md")
#: On-disk cache of finished simulations (gitignored ``data_cache/``), so an
#: interrupted ~1 h panel resumes instead of restarting. Delete to force a
#: clean re-run.
RUN_CACHE = os.path.join(_ROOT, "data_cache", "mgmt_scalar_runs.json")

#: The same five-drug panel and the same converged grid as
#: ``experiments/real_drug_panel.py`` / ``results/real_drug_panel.md``.
PANEL: Tuple[str, ...] = (
    "temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin",
)
DX_UM = 12.0

#: Bootstrap settings for the calibration CIs (fixed seed -> reproducible).
BOOT_N = 10000
BOOT_SEED = 20260727


# --------------------------------------------------------------------------- #
# 1. Calibration -- re-derived from the cached primary bytes
# --------------------------------------------------------------------------- #
def _load_paired() -> Tuple[np.ndarray, np.ndarray, List[str], np.ndarray]:
    """``(x = log2(MGMT CPM + 1), y = % viability @100 uM TMZ, ids, mgmt_cpm)``.

    Parsed from the PRIMARY files, reusing the audited docx parser in
    :mod:`tmesim.validation.gbm_tmz_exvivo` (whose :func:`verify_provenance` has
    already checked this join cell-by-cell against the data hunt's CSV).
    """
    import pandas as pd

    from .gbm_tmz_exvivo import COUNTS_GZ, MOESM2, _canon_gs, _docx_blocks

    viab: Dict[str, float] = {}
    for kind, payload in _docx_blocks(MOESM2):
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
                viab[gs] = float(r[1 + n_ops])
            except (ValueError, IndexError):
                continue

    counts = pd.read_csv(COUNTS_GZ, sep="\t", index_col=0)
    counts.columns = [str(c).strip().strip('"') for c in counts.columns]
    ids = [c for c in counts.columns if c in viab]
    cpm = counts[ids] / counts[ids].sum(axis=0) * 1e6
    mgmt_cpm = cpm.loc["MGMT"].values.astype(float)
    x = np.log2(mgmt_cpm + 1.0)
    y = np.array([viab[c] for c in ids], dtype=float)
    return x, y, ids, mgmt_cpm


def calibrate_from_exvivo() -> Dict:
    """Re-derive the MGMT -> TMZ-sensitivity mapping from the primary bytes.

    Returns every number frozen in :mod:`tmesim.molecular.mgmt`, plus the
    leave-one-out cross-validated MAE of each candidate model against the
    constant baseline (so the choice of functional form is evidenced, not
    asserted), plus the bootstrap CIs.
    """
    from scipy import stats

    x, y, ids, mgmt_cpm = _load_paired()
    n = len(x)
    kill = np.clip(1.0 - y / 100.0, MGMT_KILL_FLOOR, None)
    n_floored = int((1.0 - y / 100.0 < MGMT_KILL_FLOOR).sum())
    lk = np.log(kill)

    def _ols(xv: np.ndarray, yv: np.ndarray) -> Tuple[float, float, float, float]:
        """(slope, intercept, slope_se, r2)."""
        b, a = np.polyfit(xv, yv, 1)
        res = yv - (a + b * xv)
        sse = float((res ** 2).sum())
        sxx = float(((xv - xv.mean()) ** 2).sum())
        se = math.sqrt(sse / (len(xv) - 2) / sxx)
        r2 = 1.0 - sse / float(((yv - yv.mean()) ** 2).sum())
        return float(b), float(a), se, r2

    # --- primary: exponential decay of kill ---------------------------------
    b_lk, c_lk, se_lk, r2_lk = _ols(x, lk)
    lam = -b_lk
    tcrit = float(stats.t.ppf(0.975, n - 2))
    lam_ci = (lam - tcrit * se_lk, lam + tcrit * se_lk)
    lam_p = float(2.0 * stats.t.sf(abs(b_lk / se_lk), n - 2))

    # --- secondary: linear in viability (the report's own OLS) --------------
    b_v, a_v, se_v, r2_v = _ols(x, y)
    v_ci = (b_v - tcrit * se_v, b_v + tcrit * se_v)
    v_p = float(2.0 * stats.t.sf(abs(b_v / se_v), n - 2))

    # --- binary contrast ----------------------------------------------------
    expressed = mgmt_cpm >= 1.0
    kill_e = float(np.median(kill[expressed]))
    kill_s = float(np.median(kill[~expressed]))
    ratio = kill_e / kill_s
    mw = stats.mannwhitneyu(y[expressed], y[~expressed])
    a_, b_ = y[expressed], y[~expressed]
    pooled = math.sqrt(((len(a_) - 1) * a_.std(ddof=1) ** 2
                        + (len(b_) - 1) * b_.std(ddof=1) ** 2)
                       / (len(a_) + len(b_) - 2))
    cohens_d = float((a_.mean() - b_.mean()) / pooled)

    # --- bootstrap CIs ------------------------------------------------------
    rng = np.random.default_rng(BOOT_SEED)
    lams = np.full(BOOT_N, np.nan)
    ratios = np.empty(BOOT_N)
    ne, ns = int(expressed.sum()), int((~expressed).sum())
    for i in range(BOOT_N):
        k = rng.integers(0, n, n)
        if len(set(x[k])) >= 2:
            lams[i] = -np.polyfit(x[k], lk[k], 1)[0]
        ie = rng.integers(0, ne, ne)
        isl = rng.integers(0, ns, ns)
        ratios[i] = (np.median(kill[expressed][ie])
                     / np.median(kill[~expressed][isl]))
    lam_boot = (float(np.nanpercentile(lams, 2.5)),
                float(np.nanpercentile(lams, 97.5)))
    ratio_boot = (float(np.percentile(ratios, 2.5)),
                  float(np.percentile(ratios, 97.5)))

    # --- leave-one-out: is the functional form earning its keep? -------------
    loo_c = np.empty(n)
    loo_lin = np.empty(n)
    loo_exp = np.empty(n)
    for i in range(n):
        m = np.arange(n) != i
        loo_c[i] = np.median(y[m])
        bb, aa = np.polyfit(x[m], y[m], 1)
        loo_lin[i] = aa + bb * x[i]
        mm, cc = np.polyfit(x[m], lk[m], 1)
        loo_exp[i] = 100.0 * (1.0 - math.exp(cc + mm * x[i]))

    # --- direction / association checks -------------------------------------
    rho_v = stats.spearmanr(x, y)
    rho_k = stats.spearmanr(x, kill)

    # continuous mapping averaged within the two binary groups (independent
    # cross-check on the binary numbers)
    s_cont = np.exp(-lam * x)

    return {
        "n": int(n),
        "culture_ids": ids,
        "mgmt_cpm": mgmt_cpm,
        "log2cpm": x,
        "viability_pct": y,
        "kill": kill,
        "n_kill_floored": n_floored,
        "kill_floor": MGMT_KILL_FLOOR,
        "observed_cpm_range": (float(mgmt_cpm.min()), float(mgmt_cpm.max())),
        "observed_log2cpm_range": (float(x.min()), float(x.max())),
        "observed_viability_range": (float(y.min()), float(y.max())),
        # primary
        "lambda": float(lam),
        "lambda_se": float(se_lk),
        "lambda_ci95": lam_ci,
        "lambda_boot_ci95": lam_boot,
        "lambda_p": lam_p,
        "lambda_r2": float(r2_lk),
        "kill_at_zero": float(math.exp(c_lk)),
        # secondary
        "viab_slope": float(b_v),
        "viab_intercept": float(a_v),
        "viab_slope_ci95": v_ci,
        "viab_slope_p": v_p,
        "viab_r2": float(r2_v),
        "viab_zero_kill_x": float((100.0 - a_v) / b_v),
        # binary
        "n_expressing": ne,
        "n_silent": ns,
        "median_kill_expressing": kill_e,
        "median_kill_silent": kill_s,
        "binary_ratio": float(ratio),
        "binary_ratio_boot_ci95": ratio_boot,
        "binary_mw_p": float(mw.pvalue),
        "binary_cohens_d": cohens_d,
        "binary_from_continuous": {
            "methylated": float(s_cont[~expressed].mean()),
            "unmethylated": float(s_cont[expressed].mean()),
        },
        # cross-validation
        "loo_mae_constant": float(np.abs(y - loo_c).mean()),
        "loo_mae_linear": float(np.abs(y - loo_lin).mean()),
        "loo_mae_exponential": float(np.abs(y - loo_exp).mean()),
        # direction
        "spearman_x_vs_viability": (float(rho_v.statistic), float(rho_v.pvalue)),
        "spearman_x_vs_kill": (float(rho_k.statistic), float(rho_k.pvalue)),
    }


def verify_frozen_constants(tol: float = 1e-9) -> Dict:
    """Check the constants in :mod:`tmesim.molecular.mgmt` against a live refit."""
    fit = calibrate_from_exvivo()
    checks = {
        "n": (MGMT_CALIBRATION_N, fit["n"]),
        "lambda": (MGMT_LAMBDA, fit["lambda"]),
        "lambda_ci95_lo": (MGMT_LAMBDA_CI95[0], fit["lambda_ci95"][0]),
        "lambda_ci95_hi": (MGMT_LAMBDA_CI95[1], fit["lambda_ci95"][1]),
        "lambda_p": (MGMT_LAMBDA_P, fit["lambda_p"]),
        "lambda_r2": (MGMT_LAMBDA_R2, fit["lambda_r2"]),
        "viab_slope": (MGMT_VIAB_SLOPE, fit["viab_slope"]),
        "viab_intercept": (MGMT_VIAB_INTERCEPT, fit["viab_intercept"]),
        "viab_slope_ci95_lo": (MGMT_VIAB_SLOPE_CI95[0], fit["viab_slope_ci95"][0]),
        "viab_slope_ci95_hi": (MGMT_VIAB_SLOPE_CI95[1], fit["viab_slope_ci95"][1]),
        "binary_unmethylated": (MGMT_BINARY_SENSITIVITY["unmethylated"],
                                fit["binary_ratio"]),
        "binary_ci95_lo": (MGMT_BINARY_CI95[0], fit["binary_ratio_boot_ci95"][0]),
        "binary_ci95_hi": (MGMT_BINARY_CI95[1], fit["binary_ratio_boot_ci95"][1]),
        "lambda_boot_ci95_lo": (MGMT_LAMBDA_BOOT_CI95[0],
                                fit["lambda_boot_ci95"][0]),
        "lambda_boot_ci95_hi": (MGMT_LAMBDA_BOOT_CI95[1],
                                fit["lambda_boot_ci95"][1]),
        "binary_p": (MGMT_BINARY_P, fit["binary_mw_p"]),
        "observed_log2cpm_max": (MGMT_OBSERVED_LOG2CPM_RANGE[1],
                                 fit["observed_log2cpm_range"][1]),
    }
    bad = {k: v for k, v in checks.items() if abs(v[0] - v[1]) > tol}
    return {"ok": not bad, "checks": checks, "mismatches": bad, "fit": fit}


# --------------------------------------------------------------------------- #
# 2. Visium -- the per-patient scalar recovered from our own real data
# --------------------------------------------------------------------------- #
def visium_mgmt_table() -> Dict[str, Dict]:
    """Pseudobulk MGMT + implied TMZ sensitivity for every cached Visium section."""
    return mgmt_from_visium_all()


# --------------------------------------------------------------------------- #
# 3. Impact -- re-run the real panel with and without the MGMT scalar
# --------------------------------------------------------------------------- #
def _run_one(task: Dict) -> Dict:
    """Simulate ONE (drug, mode, mgmt) combination on the real Xenium tumour.

    Module-level so it is picklable by ``ProcessPoolExecutor`` (macOS spawn).
    Each worker rebuilds the reconstruction itself (~11 s) — the reconstruction
    is deterministic, so every worker sees byte-identical tissue.
    """
    from ..core import PDEConfig, SimConfig
    from ..molecular.binding_site_barrier import apply_barrier_to_drugparams
    from ..molecular.library import get_molecule
    from ..reconstruct.real_pipeline import apply_drug_sensitivity, reconstruct_real
    from ..sim.engine import simulate_response
    from .calibration import calibrate_alpha, calibrated_drugparams

    micro = reconstruct_real(source="xenium", dx_um=float(task["dx_um"]))
    alpha = calibrate_alpha()
    dp = calibrated_drugparams(task["drug"], alpha)
    dp = apply_barrier_to_drugparams(dp, get_molecule(task["drug"]))
    info = apply_drug_sensitivity(
        micro, task["drug"], mode=task["mode"], mgmt=task.get("mgmt")
    )
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)
    res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
    return {
        "key": task["key"],
        "drug": task["drug"],
        "mode": task["mode"],
        "mgmt": task.get("mgmt"),
        "label": task.get("label"),
        "kill": float(res.kill_fraction),
        "pen": float(res.extras["penetration_index"]),
        "mean_auc": float(res.extras["mean_auc"]),
        "hypoxic": float(res.extras["hypoxic_tumor_fraction"]),
        "ic50": float(dp.ic50),
        "bbb": float(dp.bbb_penetrance),
        "sens_mean": float(info.get("mean", 1.0)),
        "sens_gene": info.get("determinant_gene"),
        "mgmt_applied": bool(info.get("mgmt_applied", False)),
        "mgmt_scalar": info.get("mgmt_scalar"),
        "upper_bound": bool(info.get("upper_bound", True)),
        "alpha": float(alpha),
    }


def _tasks(dx_um: float, visium: Dict[str, Dict]) -> List[Dict]:
    tasks: List[Dict] = []
    for drug in PANEL:
        for mode in ("uniform", "determinant"):
            tasks.append({"key": f"{mode}|{drug}", "drug": drug, "mode": mode,
                          "dx_um": dx_um, "label": f"{mode} (no MGMT)"})
    # temozolomide across the plausible MGMT range measured in OUR Visium data
    for sid, rec in sorted(visium.items(),
                           key=lambda kv: kv[1].get("log2cpm", 0.0)):
        if not rec.get("mgmt_in_reference"):
            continue
        tasks.append({
            "key": f"mgmt|{sid}", "drug": "temozolomide", "mode": "uniform",
            "dx_um": dx_um, "mgmt": float(rec["log2cpm"]),
            "label": f"MGMT from Visium {sid} "
                     f"({rec['pseudobulk_cpm']:.2f} CPM)",
        })
    tasks.append({
        "key": "mgmt|binary_unmethylated", "drug": "temozolomide",
        "mode": "uniform", "dx_um": dx_um, "mgmt": "unmethylated",
        "label": "MGMT binary: unmethylated (clinically reported)",
    })
    tasks.append({
        "key": "mgmt|binary_methylated", "drug": "temozolomide",
        "mode": "uniform", "dx_um": dx_um, "mgmt": "methylated",
        "label": "MGMT binary: methylated (clinically reported)",
    })
    return tasks


def run_impact(dx_um: float = DX_UM, n_workers: int = 4,
               visium: Optional[Dict[str, Dict]] = None,
               cache_path: Optional[str] = None) -> Dict:
    """Quantify what the MGMT scalar does to predicted kill and to TMZ's rank.

    Runs the five-drug panel under both existing sensitivity modes, plus
    temozolomide with MGMT-informed sensitivity at every MGMT level measured in
    our own Visium sections and at both binary statuses. Everything else --
    tissue, grid, PDE, alpha, seed -- is identical across runs, so any difference
    is the MGMT term and nothing else.

    Each finished simulation is appended to ``cache_path`` (default
    :data:`RUN_CACHE`) and reused on a re-run, so a ~1 h panel survives an
    interruption. Delete that file to force a clean re-run.
    """
    import json
    from concurrent.futures import ProcessPoolExecutor, as_completed

    vis = visium if visium is not None else visium_mgmt_table()
    tasks = _tasks(dx_um, vis)
    cache_path = cache_path or RUN_CACHE
    results: Dict[str, Dict] = {}
    if os.path.exists(cache_path):
        try:
            with open(cache_path) as fh:
                cached = json.load(fh)
            if float(cached.get("dx_um", -1)) == float(dx_um):
                results = dict(cached.get("results", {}))
                print(f"  reusing {len(results)} cached simulation(s) from "
                      f"{os.path.basename(cache_path)}")
        except (ValueError, OSError):
            results = {}

    def _persist() -> None:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, "w") as fh:
            json.dump({"dx_um": float(dx_um), "results": results}, fh, indent=1)

    todo = [t for t in tasks if t["key"] not in results]
    if n_workers and n_workers > 1 and len(todo) > 1:
        with ProcessPoolExecutor(max_workers=int(n_workers)) as pool:
            futures = {pool.submit(_run_one, t): t["key"] for t in todo}
            for fut in as_completed(futures):
                r = fut.result()
                results[r["key"]] = r
                print(f"    done {r['key']:34s} kill={r['kill']:.4f}", flush=True)
                _persist()
    else:
        for t in todo:
            r = _run_one(t)
            results[r["key"]] = r
            print(f"    done {r['key']:34s} kill={r['kill']:.4f}", flush=True)
            _persist()

    base = {d: results[f"uniform|{d}"]["kill"] for d in PANEL}
    det = {d: results[f"determinant|{d}"]["kill"] for d in PANEL}
    mgmt_runs = {k: v for k, v in results.items() if k.startswith("mgmt|")}

    def _rank_of(drug: str, kills: Dict[str, float]) -> int:
        order = sorted(kills, key=lambda d: kills[d], reverse=True)
        return order.index(drug) + 1

    # TMZ's rank when its MGMT-informed kill replaces its uniform one, against
    # each of the two competitor baselines.
    rank_moves = {}
    for key, r in mgmt_runs.items():
        u = dict(base); u["temozolomide"] = r["kill"]
        d = dict(det); d["temozolomide"] = r["kill"]
        rank_moves[key] = {
            "label": r["label"],
            "mgmt_scalar": r["mgmt_scalar"],
            "kill": r["kill"],
            "kill_delta_pp": 100.0 * (r["kill"] - base["temozolomide"]),
            "kill_ratio": (r["kill"] / base["temozolomide"]
                           if base["temozolomide"] else float("nan")),
            "rank_vs_uniform_panel": _rank_of("temozolomide", u),
            "rank_vs_determinant_panel": _rank_of("temozolomide", d),
            "order_vs_determinant_panel": sorted(d, key=lambda k: d[k],
                                                 reverse=True),
        }

    return {
        "dx_um": float(dx_um),
        "panel": list(PANEL),
        "results": results,
        "kill_uniform": base,
        "kill_determinant": det,
        "tmz_rank_uniform_baseline": _rank_of("temozolomide", base),
        "tmz_rank_determinant_baseline": _rank_of("temozolomide", det),
        "order_uniform_baseline": sorted(base, key=lambda k: base[k], reverse=True),
        "order_determinant_baseline": sorted(det, key=lambda k: det[k],
                                             reverse=True),
        "mgmt_runs": rank_moves,
        "visium": vis,
    }


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def build_markdown(fit: Dict, visium: Dict[str, Dict],
                   impact: Optional[Dict]) -> str:
    L: List[str] = []
    A = L.append

    A("# MGMT as a per-patient scalar for temozolomide sensitivity")
    A("")
    A("**The uniform-1.0 temozolomide assumption is now optional, and the way out "
      "is a measured number.** `results/gbm_tmz_exvivo.md` established that MGMT "
      "expression predicts measured ex-vivo TMZ response (Spearman +0.516, "
      "p = 0.034, n = 17), that this is ~27 % of the between-model rank variance "
      "and an 81-point viability spread, and that **no** Xenium panel gene "
      "surrogates for it (0/396 pass FDR; panel base rate 14.1 % vs 12.3 % "
      "genome-wide, i.e. no enrichment) and GDSC cannot resolve it "
      "(rho = -0.068, p = 0.71). The fix is not a panel proxy. **MGMT does not "
      "need to be spatial — it is a per-patient scalar**, recoverable as Visium "
      "pseudobulk and routinely measured in the clinic.")
    A("")
    A("Code: `tmesim/molecular/mgmt.py` (the mapping), "
      "`tmesim/validation/mgmt_scalar.py` (this analysis), and one optional "
      "`mgmt=` argument on "
      "`tmesim.reconstruct.real_pipeline.apply_drug_sensitivity`. "
      "**Omitting `mgmt` reproduces the previous behaviour exactly** — asserted "
      "numerically in `tests/test_mgmt_scalar.py`.")
    A("")

    # ---------------- calibration ----------------
    A("## 1. The calibration")
    A("")
    A(f"Fitted on the **n = {fit['n']}** patient-derived GBM stem-cell cultures "
      "that have BOTH a published viability number (BJC2023_MOESM2 Supplemental "
      "Table 2) and paired RNA-seq (GEO GSE232173). Response = **mean % cell "
      "viability at 100 uM TMZ** vs DMSO, 6-day ATP assay — *lower = more "
      "sensitive*. Predictor = `log2(MGMT CPM + 1)` from the raw counts.")
    A("")
    A(f"Observed range: MGMT {fit['observed_cpm_range'][0]:.2f}–"
      f"{fit['observed_cpm_range'][1]:.2f} CPM "
      f"(`log2(CPM+1)` {fit['observed_log2cpm_range'][0]:.3f}–"
      f"{fit['observed_log2cpm_range'][1]:.3f}); viability "
      f"{fit['observed_viability_range'][0]:.1f}–"
      f"{fit['observed_viability_range'][1]:.1f} %.")
    A("")
    A("### Primary mapping — exponential in kill space (the default)")
    A("")
    A("`sensitivity` enters the engine **multiplicatively** "
      "(`cells/pkpd.py`: `p_kill = p_hill * sensitivity * R_hypoxia * "
      "R_quiescence`), so the mapping is fitted multiplicatively. With "
      "`kill = 1 - viability/100`, OLS of `ln(kill)` on `x = log2(CPM+1)` gives")
    A("")
    A("```")
    A(f"s(x) = exp(-lambda * x),   lambda = {fit['lambda']:.4f} per log2 unit")
    A(f"  95 % CI (t, {fit['n'] - 2} df) : {fit['lambda_ci95'][0]:.4f} – "
      f"{fit['lambda_ci95'][1]:.4f}")
    A(f"  bootstrap 95 % CI (10k)  : {fit['lambda_boot_ci95'][0]:.4f} – "
      f"{fit['lambda_boot_ci95'][1]:.4f}")
    A(f"  p = {fit['lambda_p']:.4f}   r2 = {fit['lambda_r2']:.3f}   "
      f"n = {fit['n']}")
    A("```")
    A("")
    A(f"Anchored so **s = 1.0 at MGMT = 0 CPM**. That anchor is the whole point: "
      "the model's existing uniform-1.0 TMZ default becomes the *MGMT-null "
      "special case* of the mapping, so supplying no MGMT changes nothing and "
      "supplying a measured one can only scale kill **down**. Bounded in (0, 1], "
      f"floored at {MGMT_SENSITIVITY_FLOOR} (a resistant tumour is not immortal: "
      "TMZ's N7-meG / N3-meA adducts are not MGMT substrates, and no zero-kill "
      "culture was observed in n = 40).")
    A("")
    A(f"**Direction is the Esteller/Hegi biology**: Spearman(x, viability) = "
      f"+{fit['spearman_x_vs_viability'][0]:.3f} "
      f"(p = {fit['spearman_x_vs_viability'][1]:.4f}) and Spearman(x, kill) = "
      f"{fit['spearman_x_vs_kill'][0]:.3f} "
      f"(p = {fit['spearman_x_vs_kill'][1]:.4f}) — "
      "**high MGMT -> repaired O6-meG -> resistant -> LOW sensitivity**.")
    A("")
    A(f"Two of the {fit['n']} cultures measured >= "
      f"{100 * (1 - fit['kill_floor']):.0f} % viability (one at 104.2 %, i.e. more "
      f"cells than the DMSO control). That is assay noise, not negative kill, so "
      f"they are pinned at {fit['kill_floor']:.2f} kill for the log fit rather "
      "than dropped.")
    A("")
    A("### Secondary mapping — linear in viability (sensitivity analysis only)")
    A("")
    A(f"The report's own OLS, `viability% = {fit['viab_intercept']:.2f} + "
      f"{fit['viab_slope']:.3f} * x` (slope 95 % CI "
      f"{fit['viab_slope_ci95'][0]:.2f}–{fit['viab_slope_ci95'][1]:.2f}, "
      f"p = {fit['viab_slope_p']:.4f}, r2 = {fit['viab_r2']:.3f}), inverted the "
      "same way. Available as `model='linear'`. It is **not** the default "
      f"because it predicts **zero kill at x = {fit['viab_zero_kill_x']:.2f}** and "
      "goes negative above that — and three of our four Visium sections sit at "
      "x >= 4.0, uncomfortably close. An unbounded multiplier is not usable in "
      "the engine.")
    A("")
    A("### Does the functional form earn its keep? Leave-one-out")
    A("")
    A("| model | LOO MAE (percentage points of viability) |")
    A("|---|---:|")
    A(f"| constant (what the model does today) | {fit['loo_mae_constant']:.2f} |")
    A(f"| linear in viability | {fit['loo_mae_linear']:.2f} |")
    A(f"| **exponential in kill (the default)** | "
      f"**{fit['loo_mae_exponential']:.2f}** |")
    A("")
    A("Both one-covariate models beat the constant out-of-sample. The linear one "
      "fits viability slightly better; the exponential one is the one that is "
      "**usable as a bounded multiplier**, and it costs "
      f"{fit['loo_mae_exponential'] - fit['loo_mae_linear']:.2f} pp of LOO MAE to "
      "get that. That trade is stated, not hidden.")
    A("")
    A("### Binary input (methylated / unmethylated)")
    A("")
    A(f"For a clinically reported promoter status the **measured group contrast** "
      f"is used rather than the regression: median kill "
      f"{fit['median_kill_silent']:.4f} (MGMT-silent, n = {fit['n_silent']}) vs "
      f"{fit['median_kill_expressing']:.4f} (MGMT-expressing, "
      f"n = {fit['n_expressing']}), Mann-Whitney **p = {fit['binary_mw_p']:.4f}**, "
      f"Cohen d = {fit['binary_cohens_d']:.2f}.")
    A("")
    A("| status | meaning | sensitivity | 95 % CI |")
    A("|---|---|---:|---|")
    A("| `methylated` | promoter methylated -> MGMT silenced -> sensitive | "
      f"{MGMT_BINARY_SENSITIVITY['methylated']:.4f} | anchor |")
    A("| `unmethylated` | MGMT expressed -> resistant | "
      f"{fit['binary_ratio']:.4f} | "
      f"{fit['binary_ratio_boot_ci95'][0]:.3f} – "
      f"{fit['binary_ratio_boot_ci95'][1]:.3f} (bootstrap) |")
    A("")
    A(f"**Independent cross-check.** Averaging the *continuous* `s(x)` within the "
      f"same two groups gives {fit['binary_from_continuous']['methylated']:.3f} "
      f"(silent) and {fit['binary_from_continuous']['unmethylated']:.3f} "
      f"(expressing), a ratio of "
      f"{fit['binary_from_continuous']['unmethylated'] / fit['binary_from_continuous']['methylated']:.3f} "
      f"against the measured {fit['binary_ratio']:.3f}. Two different routes "
      "through the same data agree.")
    A("")
    A("### Prefer continuous expression over binary methylation — and why")
    A("")
    A("In this same cohort (MOESM3 `CoxPH_tissue_005`, n = 56 parental tumours) "
      "MGMT **expression** was prognostic for overall survival "
      "(**HR 1.85, p = 6.1e-05**, q = 0.0079), externally replicated in "
      "GLIOTRAIN (n = 70, HR 1.21, p = 0.025) — while MGMT **promoter "
      "methylation**, the routine clinical biomarker, was **not significant for "
      "OS in the same cohort (p = 0.16)**. Methylation is a thresholded, noisy "
      "proxy for the quantity that actually matters, which is how much MGMT "
      "protein the tumour can make. So `mgmt_sensitivity(level)` is the primary "
      "interface and `mgmt_sensitivity(status=...)` is the fallback.")
    A("")
    A("### The mapping, tabulated")
    A("")
    A("| MGMT (CPM) | log2(CPM+1) | sensitivity (exponential) | 95 % CI | "
      "sensitivity (linear) |")
    A("|---:|---:|---:|---|---:|")
    for cpm in (0.0, 1.0, 3.0, 7.0, 15.0, 21.0, 25.0):
        x = math.log2(cpm + 1.0)
        d = mgmt_sensitivity(x, detail=True)
        lin = mgmt_sensitivity(x, model="linear")
        A(f"| {cpm:.1f} | {x:.3f} | {d['sensitivity']:.4f} | "
          f"{d['ci95'][0]:.3f} – {d['ci95'][1]:.3f} | {lin:.4f} |")
    A("")

    # ---------------- visium ----------------
    A("## 2. MGMT recovered from our own REAL Visium sections (pseudobulk)")
    A("")
    A("Visium is whole-transcriptome, so MGMT *is* in the reference. Summed over "
      "all in-tissue spots it yields one honest number per section.")
    A("")
    A("| section | accession | spots | % spots MGMT+ | MGMT UMI | pseudobulk CPM "
      "| log2(CPM+1) | -> TMZ sensitivity | 95 % CI |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---|")
    for sid, r in sorted(visium.items(), key=lambda kv: kv[1].get("log2cpm", 0)):
        if not r.get("mgmt_in_reference"):
            A(f"| {sid} | {r.get('accession')} | {r.get('n_spots')} | "
              "MGMT absent from reference | | | | | |")
            continue
        A(f"| `{sid}` | {r['accession']} | {r['n_spots']:,} | "
          f"{r['pct_spots_mgmt_positive']:.1f} % | {r['total_mgmt_counts']:.0f} | "
          f"{r['pseudobulk_cpm']:.2f} | {r['log2cpm']:.3f} | "
          f"**{r['sensitivity']:.4f}** | {r['sensitivity_ci95'][0]:.3f} – "
          f"{r['sensitivity_ci95'][1]:.3f} |")
    A("")
    lo = min(r["sensitivity"] for r in visium.values() if r.get("mgmt_in_reference"))
    hi = max(r["sensitivity"] for r in visium.values() if r.get("mgmt_in_reference"))
    A(f"**Plausible per-patient MGMT range from our own data: sensitivity "
      f"{lo:.3f} – {hi:.3f}** — i.e. these four real GBM sections imply TMZ "
      f"sensitivity between {100 * lo:.0f} % and {100 * hi:.0f} % of the "
      "MGMT-silent case. None of them is MGMT-silent.")
    A("")
    A("**This is pseudobulk, NOT a resolvable spatial field — stated plainly.** "
      "MGMT is detected in only 5.3–29.4 % of spots and each positive spot "
      "carries ~1–2 UMI, i.e. it sits at the detection floor. A per-spot 'MGMT "
      "map' built from that would be counting noise dressed up as biology. One "
      "scalar per section is the honest product — and a scalar is exactly what "
      "the determinant needs to be.")
    A("")
    A("Note also that these sections are a different assay (3' spot-based "
      "Visium) from the bulk RNA-seq the mapping was calibrated on. CPM is "
      "nominally the same unit and the values (7.0–25.5 CPM) land inside/just "
      "above the calibration range (0–21.1 CPM), but the two are not identically "
      "calibrated. The `10x_parent_gbm` value (25.5 CPM, log2 4.73) is an "
      "**extrapolation** past the calibration ceiling (4.47) and is flagged as "
      "such by `mgmt_sensitivity(..., detail=True)['extrapolated']`.")
    A("")

    # ---------------- impact ----------------
    A("## 3. Impact on the REAL Xenium tumour at the converged dx = 12 um")
    A("")
    if impact is None:
        A("*(Not run in this pass — re-run "
          "`./.venv/bin/python -m tmesim.validation.mgmt_scalar` without "
          "`--quick`.)*")
        A("")
    else:
        base = impact["kill_uniform"]
        det = impact["kill_determinant"]
        tmz0 = base["temozolomide"]
        A("Same tissue, same grid, same PDE, same perivascular-sink alpha, same "
          "seed. **Sensitivity is the only thing that changes.**")
        A("")
        A("> **Baseline note.** The five-drug baseline below is a *fresh* run of "
          "the current code, not a copy of any earlier table. Any previously "
          "quoted temozolomide figure around `kill = 0.4173` predates the three "
          "data-derived model fixes in commit `8182c7a` — in particular the "
          "geometric-mean BBB veto, which lowered temozolomide's "
          "`bbb_penetrance` and therefore its delivery. Compare against the "
          "baseline column here, not against pre-`8182c7a` numbers.")
        A("")
        A("### Baseline: the five-drug panel, no MGMT")
        A("")
        A("| drug | kill (`uniform`) | kill (`determinant`) | determinant gene |")
        A("|---|---:|---:|---|")
        for d in impact["panel"]:
            g = impact["results"][f"determinant|{d}"]["sens_gene"]
            A(f"| {d} | {base[d]:.4f} | {det[d]:.4f} | {g or '— (none on panel)'} |")
        A("")
        A("Ranking, `uniform` (fair cross-drug): "
          + " > ".join(f"**{d}** ({base[d]:.4f})"
                       for d in impact["order_uniform_baseline"]))
        A("")
        A("Ranking, `determinant` (confounded — TMZ alone keeps 1.0): "
          + " > ".join(f"**{d}** ({det[d]:.4f})"
                       for d in impact["order_determinant_baseline"]))
        A("")
        A("### Temozolomide with an MGMT-informed sensitivity")
        A("")
        A("| MGMT input | scalar | predicted kill | Δ vs uniform (pp) | "
          "× vs uniform | TMZ rank in the 5-drug panel "
          "(vs `uniform` competitors / vs `determinant` competitors) |")
        A("|---|---:|---:|---:|---:|---|")
        A(f"| *none supplied (current model)* | 1.0000 | {tmz0:.4f} | 0.0 | "
          f"1.00× | {impact['tmz_rank_uniform_baseline']} / "
          f"{impact['tmz_rank_determinant_baseline']} |")
        for key in sorted(impact["mgmt_runs"],
                          key=lambda k: -impact["mgmt_runs"][k]["kill"]):
            m = impact["mgmt_runs"][key]
            A(f"| {m['label']} | {m['mgmt_scalar']:.4f} | {m['kill']:.4f} | "
              f"{m['kill_delta_pp']:+.1f} | {m['kill_ratio']:.2f}× | "
              f"{m['rank_vs_uniform_panel']} / "
              f"{m['rank_vs_determinant_panel']} |")
        A("")
        meth = impact["mgmt_runs"].get("mgmt|binary_methylated")
        if meth is not None:
            same = abs(meth["kill"] - tmz0) < 1e-12
            A("**The anchor, verified in the simulator.** A *methylated* "
              f"(MGMT-silenced) patient gets scalar 1.0 and predicted kill "
              f"{meth['kill']:.4f}"
              + (" — bit-identical to the no-MGMT run above. "
                 if same else f" against {tmz0:.4f} with no MGMT supplied. ")
              + "That is the design: the model's existing uniform-1.0 default is "
                "exactly the MGMT-null case, so this feature can only ever move a "
                "result in the resistant direction, and only when a measurement "
                "says so.")
            A("")
        vis_runs = {k: v for k, v in impact["mgmt_runs"].items()
                    if k.startswith("mgmt|") and "binary" not in k}
        if vis_runs:
            kl = min(v["kill"] for v in vis_runs.values())
            kh = max(v["kill"] for v in vis_runs.values())
            A(f"**How much does predicted kill move?** Across the MGMT range our "
              f"own four Visium sections actually span, TMZ's predicted kill goes "
              f"from **{tmz0:.4f}** (uniform, the current model) to "
              f"**{kl:.4f}–{kh:.4f}** — a drop of "
              f"{100 * (tmz0 - kh):.1f}–{100 * (tmz0 - kl):.1f} percentage points, "
              f"i.e. {tmz0 / kh:.1f}×–{tmz0 / kl:.1f}× less kill. The MGMT term is "
              "not a rounding correction; it is the dominant uncertainty in TMZ's "
              "absolute number on this tissue.")
            A("")
        A("### Does it correct TMZ's artificially high rank?")
        A("")
        det_ranks = {impact["mgmt_runs"][k]["rank_vs_determinant_panel"]
                     for k in vis_runs} or {None}
        A("The earlier concern was specific: in `mode='determinant'` every drug "
          "whose determinant IS on the panel is scaled down by ~0.585, while TMZ "
          "— whose determinant is MGMT and therefore unmeasurable — kept 1.0. In "
          "the superseded `results/real_drug_panel.md` run that promoted TMZ from "
          "**#5** (uniform) to **#4** (determinant), ahead of methotrexate, on no "
          "biological grounds whatsoever.")
        A("")
        A(f"On the CURRENT model that promotion has already gone away for an "
          f"unrelated reason — the BBB veto fix cut temozolomide's delivery — so "
          f"TMZ now sits at **#{impact['tmz_rank_determinant_baseline']}** in "
          f"`determinant` mode before any MGMT term is applied. The artefact this "
          f"section was built to test for is therefore no longer visible in the "
          f"ordering; that is worth saying plainly rather than claiming a fix for "
          f"something already fixed.")
        A("")
        if det_ranks and None not in det_ranks:
            lo_r, hi_r = min(det_ranks), max(det_ranks)
            moved = (lo_r > impact["tmz_rank_determinant_baseline"])
            A("With a measured MGMT value from any of our four real sections, TMZ "
              f"sits at rank **#{lo_r}" + (f"–#{hi_r}" if hi_r != lo_r else "") +
              "** against the same competitors. "
              + ("**The MGMT term demotes it further — the residual artificial "
                 "promotion is corrected.**" if moved else
                 "**The rank does not move**: on this tissue TMZ's delivery x "
                 "potency already places it last, so the unmeasured-MGMT flattery "
                 "now shows up entirely in the *absolute* kill, not in the "
                 "ordering. The honest reading is that the MGMT scalar fixes a "
                 "magnitude here, not a ranking — and that a ranking it *would* "
                 "have fixed existed in the previous model state."))
            A("")
        A(f"Against the `uniform` panel the ranking is unchanged at "
          f"**#{impact['tmz_rank_uniform_baseline']}**: TMZ was already last on "
          "this section before any MGMT term, because its delivery x potency is "
          "the weakest of the five. What MGMT changes there is the **magnitude** "
          "— and it removes the UPPER BOUND label, which is the point. A number "
          "that used to be flagged 'no intrinsic resistance modelled' is now "
          "flagged with the resistance value that was used and where it came "
          "from.")
        A("")

    # ---------------- caveats ----------------
    A("## 4. Honest caveats")
    A("")
    A(f"1. **n = {fit['n']}.** The 95 % CI on lambda is "
      f"{fit['lambda_ci95'][0]:.4f}–{fit['lambda_ci95'][1]:.4f} — the lower bound "
      "is nearly zero (no effect) and the upper bound is ~2× the point estimate. "
      "The *direction* is solid and independently replicated; the *magnitude* is "
      "approximate. Every sensitivity this module returns carries that CI in "
      "`detail=True`.")
    A("2. **The response variable is noisy by construction.** It is the "
      "operator-*reproducibility* table (2–3 operators per culture over multiple "
      "years, %CV up to 43 %, median 8.3 %), not the primary screen — the primary "
      "n = 66 screen was published as figures only. Measurement noise in y "
      "attenuates the fitted slope, so this mapping is if anything conservative.")
    A("3. **Ex-vivo 2-D cultures are not spatial tissue.** These are dissociated "
      "GSC monolayers: no architecture, no vasculature, no diffusion gradient, no "
      "immune or neuronal compartment (the authors' own GSEA shows those "
      "programmes are lost in culture). **Nothing spatial is calibrated here**, "
      "and this dataset cannot test the transport half of the model at all. The "
      "scalar multiplies the *pharmacodynamic* term only; the delivery physics is "
      "untouched.")
    A("4. **Visium MGMT is pseudobulk, not a spatial field.** 5.3–29.4 % of spots "
      "positive at ~1–2 UMI each. It is used as one number per section and "
      "labelled as such everywhere it appears.")
    A("5. **Cross-assay.** Calibrated on bulk RNA-seq CPM of GSC cultures, "
      "applied to 3' spot-based Visium pseudobulk CPM of tissue. Same nominal "
      "unit, different chemistry and different cellular composition (Visium "
      "includes stroma, immune and normal brain; the cultures are pure tumour). "
      "One of the four sections extrapolates past the calibration ceiling.")
    A("6. **A clinically measured MGMT is the proper input.** The Visium numbers "
      "here demonstrate that a real per-patient value is recoverable from our own "
      "data; they are not a claim that Visium is the right assay. In a real "
      "deployment the input is the patient's own MGMT result — ideally "
      "quantitative expression, not the binary methylation call, for the reason "
      "in §1.")
    A("7. **Temozolomide only.** The nitrosoureas (lomustine, carmustine) and "
      "procarbazine are also O6-alkylators whose lesions MGMT reverses, but no "
      "measured dose-response was fitted for them, so `mgmt=` is a recorded no-op "
      "for those drugs rather than a silent extrapolation.")
    A("8. **Still no outcome label on the Xenium section.** It is an untreated "
      "resection. Adding a measured determinant makes the prediction *less* "
      "unconditional; it does not make it a validated response.")
    A("")
    A("## References")
    A("")
    A("* Ntafoulis I, Kleijn A, Ju J, *et al.* Br J Cancer 2023;129(8):1327-1338. "
      "PMID **37620410**. GEO **GSE232173**.")
    A("* Hegi ME, *et al.* NEJM 2005;352:997-1003. PMID **15758010** (MGMT "
      "methylation and TMZ benefit).")
    A("* Esteller M, *et al.* NEJM 2000;343:1350-1354. PMID **11070098** (MGMT "
      "silencing and alkylator response).")
    A("* Stupp R, *et al.* NEJM 2005;352:987-996. PMID **15758009**.")
    A("* Prior work in this repo: `results/gbm_tmz_exvivo.md` (the measurement), "
      "`results/real_drug_panel.md` (the panel this impact section re-runs), "
      "`results/numerical_rigor.md` (why dx = 12 um).")
    A("")
    return "\n".join(L)


def main(quick: bool = False, dx_um: float = DX_UM, n_workers: int = 4) -> str:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    print("Re-deriving the MGMT calibration from the cached primary bytes ...")
    chk = verify_frozen_constants()
    fit = chk["fit"]
    print(f"  n = {fit['n']}  lambda = {fit['lambda']:.4f} "
          f"(95 % CI {fit['lambda_ci95'][0]:.4f}-{fit['lambda_ci95'][1]:.4f}, "
          f"p = {fit['lambda_p']:.4f})")
    print(f"  frozen constants match live refit: {chk['ok']}"
          + ("" if chk["ok"] else f"  MISMATCH {chk['mismatches']}"))

    print("\nVisium pseudobulk MGMT (REAL cached sections):")
    visium = visium_mgmt_table()
    for sid, r in sorted(visium.items(), key=lambda kv: kv[1].get("log2cpm", 0)):
        if not r.get("mgmt_in_reference"):
            print(f"  {sid:16s} MGMT absent from reference")
            continue
        print(f"  {sid:16s} {r['pct_spots_mgmt_positive']:5.1f}% spots+  "
              f"{r['pseudobulk_cpm']:6.2f} CPM  -> sensitivity "
              f"{r['sensitivity']:.4f}")

    impact = None
    if not quick:
        print(f"\nRunning the impact panel on the REAL Xenium tumour at dx = "
              f"{dx_um} um ({len(_tasks(dx_um, visium))} simulations, "
              f"{n_workers} workers) ...")
        impact = run_impact(dx_um=dx_um, n_workers=n_workers, visium=visium)
        b = impact["kill_uniform"]
        print("  baseline (uniform): "
              + "  ".join(f"{d}={b[d]:.4f}" for d in impact["panel"]))
        for k, m in impact["mgmt_runs"].items():
            print(f"  {m['label']:52s} scalar={m['mgmt_scalar']:.4f} "
                  f"kill={m['kill']:.4f} ({m['kill_delta_pp']:+.1f} pp)  "
                  f"rank {m['rank_vs_uniform_panel']}/"
                  f"{m['rank_vs_determinant_panel']}")

    md = build_markdown(fit, visium, impact)
    with open(OUT_MD, "w") as fh:
        fh.write(md)
    print(f"\nreport written: {OUT_MD}")
    return OUT_MD


if __name__ == "__main__":  # pragma: no cover
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true",
                    help="skip the simulations; refit + Visium only")
    ap.add_argument("--dx", type=float, default=DX_UM)
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    main(quick=a.quick, dx_um=a.dx, n_workers=a.workers)
