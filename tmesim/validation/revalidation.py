"""RE-VALIDATION HARNESS — do the project's external claims still hold?

Every external validation in this repo was scored BEFORE four substantive
changes landed in the engine and reconstruction layers:

  (1) ``4ad221f`` **classifier arity fix** (``matched_null_z``) —
      ``reconstruct/marker_panels.py``, ``reconstruct/real_xenium.py``.
      Endothelial fractions moved 9.74 -> 10.50, 4.18 -> 7.07, 12.19 -> 8.75
      across the three real sections (``results/classifier_fix.md``).
  (2) ``01e6eb3`` **vessel geometry fix** — per-vessel rather than per-cell
      footprints, data-derived radius (``reconstruct/vessel_geometry.py``,
      ``reconstruct/vessel_network.py``, ``results/vessel_geometry_fix.md``).
  (3) ``ceb96c0`` **discrete vascular field promoted to default** —
      ``reconstruct/real_fields.py``; changes the UNIT of absolute exposure by
      3.36x (``results/vessel_promotion.md``).
  (4) ``0eca2fc`` **time-resolved PK** — ``core/config.py``,
      ``physics/diffusion.py``, ``sim/engine.py``, ``reconstruct/real_pipeline.py``.
      ``PDEConfig.time_resolved`` defaults OFF but the engine path changed
      (``results/dose_schedule.md``).

A validation that has not been re-run after the code changed is not a
validation. This module re-runs the five external validations under the
**current defaults** and prints a before/after table.

DESIGN CONSTRAINT — THIS MODULE ONLY READS
------------------------------------------
It imports the five validation modules and calls their *compute* entry points
only. It never calls their ``write_report`` / ``make_figure`` writers, so
re-validating does not clobber the "before" evidence on disk
(``results/glioma_regional_pk.md`` and friends stay exactly as they were when
the claim was made). Nothing in ``core``, ``physics``, ``molecular``,
``reconstruct``, ``sim`` or the five validation modules is modified.

The BEFORE numbers below are transcribed from the committed reports and are
labelled with the file and line they came from, so the comparison cannot drift
silently.

Usage::

    python -m tmesim.validation.revalidation            # all five
    python -m tmesim.validation.revalidation glioma     # one at a time
    python -m tmesim.validation.revalidation --write-md # also write results/revalidation.md
"""
from __future__ import annotations

import json
import math
import os
import sys
import time
import traceback
from typing import Callable, Dict, List, Optional

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(_REPO_ROOT, "results")
REPORT_MD = os.path.join(RESULTS_DIR, "revalidation.md")
CACHE_JSON = os.path.join(RESULTS_DIR, ".revalidation_cache.json")


# ---------------------------------------------------------------------------
# The BEFORE side of the evidence, transcribed with provenance.
# ---------------------------------------------------------------------------
#: name -> (before value, tolerance as a RELATIVE fraction, provenance)
#: A tolerance of 0.0 means "must be byte-identical to the reported precision".
BEFORE = {
    # ---- (a) glioma_regional_pk ------------------------------------------
    "glioma.gm_abs_fold_error_postfix": (1.995, "results/glioma_regional_pk.md:227"),
    "glioma.n_pass_10x_postfix": (6, "results/glioma_regional_pk.md:201"),
    "glioma.n_pass_3x_postfix": (4, "results/glioma_regional_pk.md:201"),
    "glioma.worst_abs_fold_error_postfix": (6.94, "results/glioma_regional_pk.md:201"),
    "glioma.lev_model_ratio": (1.19, "results/glioma_regional_pk.md:255"),
    "glioma.cef_model_ratio": (5.87, "results/glioma_regional_pk.md:258"),
    "glioma.lev_measured_ratio": (1.45, "results/glioma_regional_pk.md:86"),
    "glioma.cef_measured_ratio": (4.26, "results/glioma_regional_pk.md:102"),
    "glioma.anchor_gm_ratio_error": (1.30, "results/glioma_regional_pk.md:290"),
    # ---- (b) pk_multidrug --------------------------------------------------
    "pk.alpha": (148.6, "results/pk_multidrug.md:11"),
    "pk.mitoxantrone_L_um": (48.2, "results/pk_multidrug.md:16"),
    "pk.topotecan_L_um": (83.6, "results/pk_multidrug.md:17"),
    "pk.n_oos_pass": (2, "results/pk_multidrug.md:50"),
    # ---- (c) gist_imatinib -------------------------------------------------
    "gist.m1_geomean_fold_error": (7.92, "results/gist_imatinib.md:168"),
    "gist.n": (6, "results/gist_imatinib.md:168"),
    "gist.null_geomean_fold_error": (12.96, "results/gist_imatinib.md:174"),
    "gist.m1_spearman_rho": (0.655, "results/gist_imatinib.md:168"),
    # ---- (d) reconstruction_fidelity --------------------------------------
    "recon.pearson_4type": (0.997, "results/reconstruction_fidelity.md:16"),
    "recon.pearson_nontumor": (0.945, "results/reconstruction_fidelity.md:17"),
    "recon.xenium_endothelial_pct": (9.74, "results/reconstruction_fidelity.md:13"),
    "recon.mean_abs_diff_pct": (4.44, "results/reconstruction_fidelity.md:15"),
    # Not published as a number, but it IS the report's stated verdict: "two
    # independent real assays + estimators agree on the composition ORDERING".
    # Recomputed from the published before-table (73.80/7.92/8.54/9.74 vs
    # 82.66/1.04/6.53/9.77), where the two rankings are identical -> rho = 1.0.
    "recon.ordering_spearman": (1.0, "results/reconstruction_fidelity.md:11-14 (derived)"),
    "recon.PECAM1_d": (1.29, "results/reconstruction_fidelity.md:31"),
    "recon.VEGFA_d": (1.35, "results/reconstruction_fidelity.md:33"),
    "recon.CA9_d": (0.67, "results/reconstruction_fidelity.md:32"),
    "recon.PECAM1_p": (1.0e-4, "results/reconstruction_fidelity.md:31"),
    "recon.VEGFA_p": (1.0e-4, "results/reconstruction_fidelity.md:33"),
    "recon.CA9_p": (1.0e-4, "results/reconstruction_fidelity.md:32"),
    # ---- (e) pharmacology_real --------------------------------------------
    "pharm.spearman": (0.827, "results/pharmacology_real.md:55"),
    "pharm.p_value": (2.35e-05, "results/pharmacology_real.md:55"),
    "pharm.n_overlap": (18, "results/pharmacology_real.md:53"),
}

#: Which of the four changes each validation's dependency cone can even see.
#: Derived from the import graph, not asserted. Used to attribute movement.
DEPENDENCY_NOTE = {
    "glioma": (
        "imports molecular.descriptors / molecular.library / molecular.permeability "
        "/ reconstruct.barrier ONLY. None of those files appear in 4ad221f, "
        "01e6eb3, ceb96c0 or 0eca2fc. Structurally blind to all four changes."
    ),
    "pk": (
        "imports molecular.{descriptors,library,compile_drug}, core.DrugParams, "
        "validation.{pk_anchors,calibration,pk_reality}. calibration.calibrate_alpha "
        "solves L = sqrt(D_eff/k_eff) analytically on a homogeneous slab -- it does "
        "not build a reconstruction, so it never touches the vascular field, the "
        "classifier or the engine's time loop."
    ),
    "gist": (
        "imports molecular.{descriptors,compile_drug}, validation.{calibration,"
        "pk_reality}, plus its own on-disk imatinib MALDI-MSI tables. Same "
        "analytic-slab path as pk_multidrug; no reconstruction, no engine."
    ),
    "recon": (
        "CHECK 1 calls reconstruct.real_pipeline.reconstruct_real('xenium'), which "
        "runs the full Xenium reconstruction -- marker_panels + real_xenium (the "
        "classifier, change 1), vessel_geometry/vessel_network (change 2) and "
        "real_fields (change 3). This is the ONLY one of the five inside the blast "
        "radius. CHECK 2 (region markers) uses reconstruct.region_label on Visium "
        "expression and is independent of all four."
    ),
    "pharm": (
        "imports data.real.drug_response, molecular.descriptors, molecular.library. "
        "Pure table join over two wet-lab screens; no model in the loop for the "
        "concordance number."
    ),
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _rel(before: float, after: float) -> float:
    """Relative change, guarding a zero baseline."""
    if before == 0:
        return float("inf") if after != 0 else 0.0
    return (after - before) / abs(before)


def _verdict_for(key: str, before, after) -> str:
    """UNCHANGED / MOVED, at the precision the claim was published to."""
    if isinstance(before, int) and isinstance(after, (int, float)):
        return "UNCHANGED" if int(round(after)) == before else "MOVED"
    if not isinstance(after, (int, float)) or not math.isfinite(after):
        return "ERROR"
    # Compare at the number of significant digits the report published.
    txt = f"{before:g}"
    if "." in txt:
        dp = len(txt.split(".")[1].rstrip("0"))
        dp = max(dp, 0)
    else:
        dp = 0
    if abs(before) < 1e-3:  # p-values published as 1e-4 / 2.35e-05
        return "UNCHANGED" if abs(_rel(before, after)) < 0.05 else "MOVED"
    return "UNCHANGED" if round(after, dp) == round(before, dp) else "MOVED"


def _retry_once(fn: Callable, *a, **kw):
    """Run ``fn``; on a network/IO-shaped failure, retry exactly once.

    The reconstruction and drug-response loaders hit the data cache and, on a
    cold cache, the network. A single retry covers a transient fetch; anything
    that fails twice is a real defect and is re-raised.
    """
    try:
        return fn(*a, **kw)
    except Exception as exc:  # noqa: BLE001 - deliberately broad, then re-raised
        transient = any(
            t in f"{type(exc).__name__}: {exc}".lower()
            for t in ("timeout", "connection", "network", "temporarily",
                      "http", "urlopen", "ssl", "resolve")
        )
        if not transient:
            raise
        print(f"  [retry] transient failure ({type(exc).__name__}: {exc}); retrying once")
        time.sleep(2.0)
        return fn(*a, **kw)


# ---------------------------------------------------------------------------
# (a) glioma_regional_pk — the 6-drug measured human glioma scorecard
# ---------------------------------------------------------------------------
def revalidate_glioma(n_perm: int = 20000, seed: int = 0) -> Dict:
    """Re-run the scorecard + the regional CE-vs-non-enhancing test.

    Calls ``run_scorecard`` / ``run_regional_test`` / ``build_post_fix`` -- the
    same three functions ``glioma_regional_pk.main`` calls -- but NOT
    ``write_report`` or ``make_figure``.
    """
    from . import glioma_regional_pk as G

    regional = _retry_once(G.run_regional_test, n_perm=n_perm, seed=seed)
    post_fix = _retry_once(G.build_post_fix, regional)

    sp = post_fix["score_post"]
    rp = post_fix["regional_post"]
    head = {c["drug"]: c for c in rp["headline"]}

    out = {
        "glioma.gm_abs_fold_error_postfix": float(sp["geometric_mean_abs_fold_error"]),
        "glioma.n_pass_10x_postfix": int(sp["n_pass_10x"]),
        "glioma.n_pass_3x_postfix": int(sp["n_pass_3x"]),
        "glioma.worst_abs_fold_error_postfix": float(sp["worst_abs_fold_error"]),
        "glioma.lev_model_ratio": float(head["levetiracetam"]["model_ratio_postfix"]),
        "glioma.cef_model_ratio": float(head["cefazolin"]["model_ratio_postfix"]),
        "glioma.lev_measured_ratio": float(head["levetiracetam"]["measured_ratio"]),
        "glioma.cef_measured_ratio": float(head["cefazolin"]["measured_ratio"]),
        "glioma.anchor_gm_ratio_error": float(
            rp["headline_geometric_mean_abs_ratio_error"]),
    }
    out["_detail"] = {
        "per_drug": [
            {"drug": r["drug"], "measured": r["measured_mean"],
             "predicted": r["predicted_partition"], "fold": r["fold_error"]}
            for r in sp["rows"]
        ],
        "leak_constant": rp.get("leak_constant"),
    }
    return out


# ---------------------------------------------------------------------------
# (b) pk_multidrug — the out-of-sample penetration anchors
# ---------------------------------------------------------------------------
def revalidate_pk() -> Dict:
    from . import pk_multidrug as P

    pen = _retry_once(P.run_penetration)
    rows = {r["drug"]: r for r in pen["rows"]}
    out = {
        "pk.alpha": float(pen["alpha"]),
        "pk.mitoxantrone_L_um": float(rows["mitoxantrone"]["predicted_L_um"]),
        "pk.topotecan_L_um": float(rows["topotecan"]["predicted_L_um"]),
        "pk.n_oos_pass": int(pen["n_oos_pass"]),
    }
    out["_detail"] = {
        "rows": [
            {"drug": d, "L_um": r["predicted_L_um"],
             "lequiv_band_um": list(r["lequiv_band_um"]),
             "verdict": r["verdict"], "in_sample": r["in_sample"]}
            for d, r in rows.items()
        ],
        "n_oos": pen["n_oos"],
        "dox_target_um": pen["dox_target_um"],
    }
    return out


# ---------------------------------------------------------------------------
# (c) gist_imatinib — the tumour:normal partition
# ---------------------------------------------------------------------------
def revalidate_gist() -> Dict:
    from . import gist_imatinib as GI

    char = _retry_once(GI.characterize)
    targets = GI.measured_targets(char)
    null = GI.plasma_null(char)
    compiled = GI.compile_imatinib()
    pred = GI.predict(compiled)
    sc = GI.score(targets, compiled, pred, null)

    out = {
        "gist.m1_geomean_fold_error": float(sc["M1"]["geomean_fold_error"]),
        "gist.n": int(sc["n"]),
        "gist.null_geomean_fold_error": float(sc["null_unity"]["geomean_fold_error"]),
        "gist.m1_spearman_rho": float(sc["M1"]["spearman_rho"]),
    }
    out["_detail"] = {
        "max_fold": sc["M1"]["max_fold_error"],
        "min_fold": sc["M1"]["min_fold_error"],
        "n_within_3x": sc["M1"]["n_within_3x"],
        "beats_plasma_null": sc["beats_plasma_null"],
    }
    return out


# ---------------------------------------------------------------------------
# (d) reconstruction_fidelity — cross-platform + region markers
# ---------------------------------------------------------------------------
def revalidate_recon(n_perm: int = 10000) -> Dict:
    from . import reconstruction_fidelity as RF

    from scipy.stats import spearmanr

    cp = _retry_once(RF.check_cross_platform_composition)
    rm = _retry_once(RF.check_region_markers, n_perm=n_perm)
    t = rm["tests"]

    # The report's stated verdict is about ORDERING, not the Pearson magnitude:
    # "two independent real assays + estimators agree on the composition
    # ordering (tumour-dominant, minority immune/endothelial/stromal)". Score it
    # directly as the rank correlation of the two 4-class rankings.
    order_types = list(RF._SHARED_TYPES)
    xo = [cp["xenium_fractions"][k] for k in order_types]
    vo = [cp["visium_fractions"][k] for k in order_types]
    ordering_rho = float(spearmanr(xo, vo).correlation)

    out = {
        "recon.ordering_spearman": ordering_rho,
        "recon.pearson_4type": float(cp["pearson_4type"]),
        "recon.pearson_nontumor": float(cp["pearson_nontumor"]),
        "recon.xenium_endothelial_pct": 100.0 * float(cp["xenium_fractions"]["ENDOTHELIAL"]),
        "recon.mean_abs_diff_pct": 100.0 * float(cp["mean_abs_diff"]),
        "recon.PECAM1_d": float(t["PECAM1@perivascular"]["cohens_d"]),
        "recon.VEGFA_d": float(t["VEGFA@necrotic"]["cohens_d"]),
        "recon.CA9_d": float(t["CA9@necrotic"]["cohens_d"]),
        "recon.PECAM1_p": float(t["PECAM1@perivascular"]["perm_p_value"]),
        "recon.VEGFA_p": float(t["VEGFA@necrotic"]["perm_p_value"]),
        "recon.CA9_p": float(t["CA9@necrotic"]["perm_p_value"]),
    }
    out["_detail"] = {
        "xenium_n_cells": cp["xenium_n_cells"],
        "visium_n_spots": cp["visium_n_spots"],
        "xenium_fractions_pct": {k: 100.0 * v for k, v in cp["xenium_fractions"].items()},
        "visium_fractions_pct": {k: 100.0 * v for k, v in cp["visium_fractions"].items()},
        "region_distribution": rm["region_distribution"],
    }
    return out


# ---------------------------------------------------------------------------
# (e) pharmacology_real — GDSC vs HGCC concordance
# ---------------------------------------------------------------------------
def revalidate_pharm() -> Dict:
    from . import pharmacology_real as PR

    hgcc = _retry_once(PR.hgcc_mean_sensitivity)
    conc = _retry_once(PR.gdsc_hgcc_concordance, hgcc)
    out = {
        "pharm.spearman": float(conc["spearman"]),
        "pharm.p_value": float(conc["p_value"]),
        "pharm.n_overlap": int(conc["n_overlap"]),
    }
    out["_detail"] = {
        "resolve_ic50_measured_library_drugs":
            conc.get("resolve_ic50_measured_library_drugs", []),
    }
    return out


# ---------------------------------------------------------------------------
CHECKS = {
    "glioma": ("(a) glioma_regional_pk -- 6-drug measured human glioma scorecard "
               "+ regional CE test", revalidate_glioma),
    "pk": ("(b) pk_multidrug -- out-of-sample penetration anchors", revalidate_pk),
    "gist": ("(c) gist_imatinib -- tumour:normal partition", revalidate_gist),
    "recon": ("(d) reconstruction_fidelity -- cross-platform + region markers",
              revalidate_recon),
    "pharm": ("(e) pharmacology_real -- GDSC vs HGCC concordance", revalidate_pharm),
}


def run_one(name: str) -> Dict:
    title, fn = CHECKS[name]
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")
    t0 = time.time()
    try:
        res = fn()
        res["_status"] = "ok"
    except Exception as exc:  # noqa: BLE001
        res = {"_status": "CRASH",
               "_error": f"{type(exc).__name__}: {exc}",
               "_traceback": traceback.format_exc()}
        print(f"  !! CRASH  {res['_error']}")
    res["_seconds"] = round(time.time() - t0, 1)
    print(f"  ({res['_seconds']}s)")
    for k, v in res.items():
        if k.startswith("_"):
            continue
        before, prov = BEFORE.get(k, (None, ""))
        if before is None:
            continue
        vd = _verdict_for(k, before, v)
        bt = f"{before:g}" if isinstance(before, float) else str(before)
        at = f"{v:.6g}" if isinstance(v, float) else str(v)
        print(f"    {k:42s} before={bt:>10s}  after={at:>12s}  {vd}")
    return res


def _load_cache() -> Dict:
    if os.path.exists(CACHE_JSON):
        try:
            with open(CACHE_JSON) as f:
                return json.load(f)
        except Exception:  # noqa: BLE001
            return {}
    return {}


def _save_cache(d: Dict) -> None:
    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(CACHE_JSON, "w") as f:
        json.dump(d, f, indent=1, default=str)


def main(names: Optional[List[str]] = None, use_cache: bool = True) -> Dict:
    names = names or list(CHECKS)
    cache = _load_cache() if use_cache else {}
    for n in names:
        cache[n] = run_one(n)
        _save_cache(cache)
    return cache


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    main(args or None, use_cache="--no-cache" not in sys.argv)
