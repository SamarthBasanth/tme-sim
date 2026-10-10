"""Fit and HONESTLY score surrogates of the transport solver.

Reads ``sweep_data.json`` (produced by ``sweep.py``, every row a real
``execute_run``) and answers one question per target: can a model predict the
SOLVER's number closely enough to be worth showing, on tissue it has never
seen?

THE SPLIT IS THE WHOLE EXPERIMENT
---------------------------------
Three splits are scored for every target, deliberately, because the gap between
them is the result:

``random``
    5-fold KFold over rows, ignoring tissue. This is the WRONG split and it is
    computed only to show how much it flatters: rows from the same section are
    in train and test, so the model is scored on a tissue it has memorised.
``section``
    Leave-one-section-out. Still leaky on GBM-Space, which ships up to five
    sections per patient.
``donor``
    Leave-one-DONOR-out. The honest one. Every section of a patient is held out
    together, so the fold is a tissue the model has never seen in any form.
    This is the number to quote.

A fourth split, ``compound``, leaves one compound out entirely, because
bring-your-own-SMILES is a shipped product surface and a surrogate that only
works on the six library compounds it was trained on cannot serve it.

BASELINES THE MODEL MUST BEAT
-----------------------------
``mean``
    The training set's central value for the target (geometric, for the targets
    fitted in log space). No features, no model. A surrogate that does not beat
    this is worthless, and reporting R-squared instead of this comparison is how
    a worthless surrogate gets shipped.
``analytic``
    For the penetration length ONLY: ``validation.delivery_metrics.
    penetration_length_um(dp, ecm_density)``, the closed form the product
    ALREADY computes for free on every run. For that target this, not the mean,
    is the thing to beat -- it costs nothing and it is already shipped.

Writes ``experiments/surrogate_feasibility/fit_results.json``. Prints a table.
"""
from __future__ import annotations

import json
import math
import os
import time
from typing import Dict, List, Optional, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_IN = os.path.join(_HERE, "sweep_data.json")
_OUT = os.path.join(_HERE, "fit_results.json")

# --------------------------------------------------------------------------- #
# Targets. ``log`` means the model is fitted on log10(y) -- the right space for
# a strictly positive ratio or length -- and every error is reported AFTER
# inverting the transform, in the unit a user would read.
# --------------------------------------------------------------------------- #
TARGETS = {
    "E_p90_over_p10": {
        "unit": "x (dimensionless exposure spread)",
        "log": True,
        "desc": "ratio of the p90 to the p10 of per-tumour-cell exposure",
    },
    "penetration_length_um": {
        "unit": "um",
        "log": True,
        "desc": "fitted exp+floor length of AUC vs distance-from-vessel",
    },
    "median_auc_over_peak": {
        "unit": "s (median tumour AUC per uM of plasma peak)",
        "log": True,
        "desc": "dose-invariant normalised exposure magnitude",
    },
}

# Features the Studio has WITHOUT the transport solve: the cached
# reconstruction's summary statistics, the compound's descriptors, and the
# RunSpec's PK numbers. Nothing here is a solve output.
TISSUE_FEATURES = [
    "n_cells", "n_tumor", "frac_TUMOR", "frac_STROMAL", "frac_IMMUNE",
    "frac_ENDOTHELIAL", "frac_NECROTIC", "mean_vascular_density",
    "mean_ecm_density", "mean_oxygen", "tissue_coverage_fraction",
    "max_darcy_speed_um_s",
]
COMPOUND_FEATURES = [
    "drug_diffusivity_um2_s", "drug_binding_1_per_s",
    "drug_clearance_1_per_s", "drug_bbb_penetrance", "drug_mw", "drug_logp",
]
PK_FEATURES = ["half_life_s", "n_doses", "total_time_s", "plasma_peak"]


def _f(v, default=np.nan) -> float:
    try:
        x = float(v)
        return x if math.isfinite(x) else default
    except (TypeError, ValueError):
        return default


def build(rows: List[Dict]) -> Tuple[np.ndarray, List[str], Dict, Dict]:
    """Feature matrix, feature names, targets, and the row metadata."""
    feats, names = [], []
    for r in rows:
        v, n = [], []
        for k in TISSUE_FEATURES:
            v.append(_f(r.get(k))); n.append(k)
        grid = r.get("grid") or [np.nan, np.nan]
        ext = r.get("extent_um") or [np.nan, np.nan]
        v += [_f(grid[0]), _f(grid[1]),
              _f(ext[0]), _f(ext[1]),
              _f(ext[0]) * _f(ext[1]) / 1e6,
              math.log10(max(_f(r.get("n_cells"), 1.0), 1.0))]
        n += ["grid_rows", "grid_cols", "extent_x_um", "extent_y_um",
              "area_mm2", "log10_n_cells"]
        for k in COMPOUND_FEATURES:
            v.append(_f(r.get(k))); n.append(k)
        D = _f(r.get("drug_diffusivity_um2_s"))
        clr = _f(r.get("drug_clearance_1_per_s"))
        bnd = _f(r.get("drug_binding_1_per_s"))
        hl = _f(r.get("half_life_s"))
        tt = _f(r.get("total_time_s"))
        # Physics-motivated derived features. A linear model on raw D and
        # clearance cannot express sqrt(D/k); handing it the group explicitly
        # is the difference between testing the features and testing the
        # learner.
        v += [math.log10(D) if D > 0 else np.nan,
              math.log10(clr) if clr > 0 else -12.0,
              math.log10(bnd) if bnd > 0 else -12.0,
              math.sqrt(D / clr) if (D > 0 and clr > 0) else np.nan,
              math.sqrt(D * tt) if (D > 0 and tt > 0) else np.nan,
              math.log10(hl) if hl > 0 else np.nan,
              math.log10(tt) if tt > 0 else np.nan,
              tt / hl if hl > 0 else np.nan,
              _f(r.get("penetration_length_analytic_um"))]
        n += ["log10_D", "log10_clearance", "log10_binding", "sqrt_D_over_k",
              "sqrt_D_t", "log10_half_life", "log10_total_time",
              "t_over_half_life", "L_analytic_um"]
        for k in PK_FEATURES:
            v.append(_f(r.get(k))); n.append(k)
        feats.append(v)
        if not names:
            names = n
    X = np.asarray(feats, dtype=float)

    y: Dict[str, np.ndarray] = {}
    e10 = np.array([_f(r.get("E_p10_uM")) for r in rows])
    e90 = np.array([_f(r.get("E_p90_uM")) for r in rows])
    with np.errstate(divide="ignore", invalid="ignore"):
        y["E_p90_over_p10"] = e90 / e10
    y["penetration_length_um"] = np.array(
        [_f(r.get("penetration_length_um")) for r in rows])
    med = np.array([_f(r.get("median_tumor_auc")) for r in rows])
    pk = np.array([_f(r.get("plasma_peak_uM")) for r in rows])
    with np.errstate(divide="ignore", invalid="ignore"):
        y["median_auc_over_peak"] = med / pk

    meta = {
        "donor": np.array([r.get("donor") for r in rows]),
        "section": np.array([r.get("section") for r in rows]),
        "drug": np.array([r.get("drug") for r in rows]),
        "valid": np.array([bool(r.get("spatial_spread_valid")) for r in rows]),
        "source": np.array([r.get("source") for r in rows]),
        "L_analytic": np.array(
            [_f(r.get("penetration_length_analytic_um")) for r in rows]),
    }
    return X, names, y, meta


# --------------------------------------------------------------------------- #
# Models
# --------------------------------------------------------------------------- #
def _models():
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.impute import SimpleImputer
    from sklearn.linear_model import RidgeCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    return {
        "ridge": lambda: make_pipeline(
            SimpleImputer(strategy="median"), StandardScaler(),
            RidgeCV(alphas=np.logspace(-3, 4, 40))),
        # Hyperparameters are DECLARED, not tuned against any held-out score,
        # following the standing practice in tmesim/ml/panel_imputation.py:
        # tuning against the numbers reported here would make them optimistic
        # by construction. early_stopping is off so a fold's score does not
        # depend on an internal split of that fold.
        "hgb": lambda: HistGradientBoostingRegressor(
            max_iter=200, learning_rate=0.08, max_depth=None,
            min_samples_leaf=8, l2_regularization=1.0,
            early_stopping=False, random_state=0),
    }


def _folds(split: str, meta: Dict, n: int):
    """(train_idx, test_idx, fold_label) for one split strategy."""
    if split == "random":
        from sklearn.model_selection import KFold
        for i, (tr, te) in enumerate(
                KFold(n_splits=5, shuffle=True, random_state=0).split(
                    np.arange(n))):
            yield tr, te, f"kfold{i}"
        return
    key = {"section": "section", "donor": "donor", "compound": "drug"}[split]
    groups = meta[key]
    for g in sorted(set(groups.tolist())):
        te = np.where(groups == g)[0]
        tr = np.where(groups != g)[0]
        if len(te) and len(tr):
            yield tr, te, str(g)


def evaluate(X, y, meta, target: str, split: str, subset=None) -> Dict:
    """Out-of-fold predictions and errors, in the target's own unit."""
    spec = TARGETS[target]
    use_log = bool(spec["log"])
    yy = y[target]
    ok = np.isfinite(yy)
    if use_log:
        ok &= yy > 0
    if subset is not None:
        ok &= subset
    idx = np.where(ok)[0]
    # +-inf would break both learners; NaN is handled (SimpleImputer for ridge,
    # natively by HistGradientBoosting). Rows are never dropped for a missing
    # FEATURE -- only for a missing TARGET, which is a run the solver refused.
    Xs = np.where(np.isfinite(X[idx]), X[idx], np.nan)
    ys = yy[idx]
    sub_meta = {k: v[idx] for k, v in meta.items()}
    t = np.log10(ys) if use_log else ys

    preds = {m: np.full(len(idx), np.nan) for m in
             (list(_models()) + ["mean"] +
              (["analytic"] if target == "penetration_length_um" else []))}
    fold_of = np.empty(len(idx), dtype=object)
    builders = _models()
    train_s = {m: 0.0 for m in builders}
    infer_s = {m: 0.0 for m in builders}
    n_folds = 0
    for tr, te, lab in _folds(split, sub_meta, len(idx)):
        n_folds += 1
        fold_of[te] = lab
        preds["mean"][te] = t[tr].mean()
        if target == "penetration_length_um":
            la = sub_meta["L_analytic"][te]
            with np.errstate(divide="ignore", invalid="ignore"):
                preds["analytic"][te] = np.log10(np.where(la > 0, la, np.nan))
        for name, build_m in builders.items():
            m = build_m()
            t0 = time.perf_counter()
            m.fit(Xs[tr], t[tr])
            train_s[name] += time.perf_counter() - t0
            t0 = time.perf_counter()
            preds[name][te] = m.predict(Xs[te])
            infer_s[name] += time.perf_counter() - t0

    out = {"target": target, "split": split, "n": int(len(idx)),
           "n_folds": n_folds, "unit": spec["unit"], "log_fitted": use_log,
           "y_min": float(np.min(ys)), "y_max": float(np.max(ys)),
           "y_median": float(np.median(ys)),
           "y_spread_p90_over_p10": float(
               np.percentile(ys, 90) / np.percentile(ys, 10))
           if np.percentile(ys, 10) > 0 else None,
           "models": {}}
    # Across-GROUP spread of the target: the number an error has to be small
    # against. Computed on per-group medians so within-group config variation
    # does not inflate it.
    gk = {"random": "section", "section": "section", "donor": "donor",
          "compound": "drug"}[split]
    gmed = np.array([np.median(ys[sub_meta[gk] == g])
                     for g in sorted(set(sub_meta[gk].tolist()))])
    out["across_group_median_min"] = float(np.min(gmed))
    out["across_group_median_max"] = float(np.max(gmed))
    out["across_group_fold_spread"] = (float(np.max(gmed) / np.min(gmed))
                                       if np.min(gmed) > 0 else None)
    out["across_group_abs_range"] = float(np.max(gmed) - np.min(gmed))

    for name, p in preds.items():
        good = np.isfinite(p)
        if good.sum() == 0:
            out["models"][name] = {"n_scored": 0}
            continue
        yhat = (10.0 ** p[good]) if use_log else p[good]
        ytrue = ys[good]
        ae = np.abs(yhat - ytrue)
        with np.errstate(divide="ignore", invalid="ignore"):
            fold_err = np.maximum(yhat / ytrue, ytrue / yhat)
        rec = {
            "n_scored": int(good.sum()),
            "mae": float(np.mean(ae)),
            "median_ae": float(np.median(ae)),
            "p90_ae": float(np.percentile(ae, 90)),
            "max_ae": float(np.max(ae)),
            "median_fold_error": float(np.median(fold_err[np.isfinite(fold_err)])),
            "p90_fold_error": float(np.percentile(
                fold_err[np.isfinite(fold_err)], 90)),
            "r2": float(1.0 - np.sum((yhat - ytrue) ** 2)
                        / max(np.sum((ytrue - ytrue.mean()) ** 2), 1e-30)),
        }
        if name in train_s:
            rec["train_s_total"] = round(train_s[name], 3)
            rec["infer_s_total"] = round(infer_s[name], 4)
            rec["infer_s_per_row"] = round(
                infer_s[name] / max(int(good.sum()), 1), 6)
        out["models"][name] = rec

    base = out["models"].get("mean", {}).get("mae")
    for name, rec in out["models"].items():
        if base and rec.get("mae") is not None:
            rec["skill_vs_mean"] = round(1.0 - rec["mae"] / base, 4)
    an = out["models"].get("analytic", {}).get("mae")
    if an:
        for name, rec in out["models"].items():
            if rec.get("mae") is not None:
                rec["skill_vs_analytic"] = round(1.0 - rec["mae"] / an, 4)
    return out


#: Which feature-name prefixes belong to which group, for the ablation.
GROUPS = {
    "tissue": set(TISSUE_FEATURES) | {
        "grid_rows", "grid_cols", "extent_x_um", "extent_y_um", "area_mm2",
        "log10_n_cells"},
    "compound": set(COMPOUND_FEATURES) | {
        "log10_D", "log10_clearance", "log10_binding", "sqrt_D_over_k",
        "L_analytic_um"},
    "pk": set(PK_FEATURES) | {
        "log10_half_life", "log10_total_time", "t_over_half_life", "sqrt_D_t"},
}


def ablate(X, y, meta, names: List[str], split: str = "donor") -> Dict:
    """Does the model learn the TISSUE, or only the compound?

    A surrogate whose whole skill comes from the compound and PK columns is a
    per-compound lookup table. It would reproduce the solver's answer for a
    compound averaged over tissues and be blind to the one thing a user is
    asking the section for -- and the project's central finding is precisely
    that its answers do not travel across donors. So each feature group is
    scored alone, and all of them together, on the honest split.
    """
    out: Dict = {}
    idx_of = {n: i for i, n in enumerate(names)}
    combos = {
        "compound_only": GROUPS["compound"] | GROUPS["pk"],
        "tissue_only": GROUPS["tissue"] | GROUPS["pk"],
        "pk_only": GROUPS["pk"],
        "all": set(names),
    }
    for target in TARGETS:
        rec: Dict = {}
        for cname, keep in combos.items():
            cols = [idx_of[n] for n in names if n in keep]
            try:
                e = evaluate(X[:, cols], y, meta, target, split)
            except Exception as exc:  # noqa: BLE001
                rec[cname] = {"error": f"{type(exc).__name__}: {exc}"}
                continue
            rec[cname] = {
                "n_features": len(cols),
                **{m: {"mae": e["models"][m].get("mae"),
                       "skill_vs_mean": e["models"][m].get("skill_vs_mean")}
                   for m in ("mean", "ridge", "hgb") if m in e["models"]},
            }
        out[target] = rec
    return out


def cost_of_final_models(X, y, meta) -> Dict:
    """Train on EVERYTHING and measure what shipping one would cost.

    Training time, serialised size, and SINGLE-ROW inference latency -- a
    preview predicts one run, not a batch, so a per-row figure taken from a
    batch predict would understate it.
    """
    import pickle

    out: Dict = {}
    for target, spec in TARGETS.items():
        yy = y[target]
        ok = np.isfinite(yy) & (yy > 0 if spec["log"] else True)
        idx = np.where(ok)[0]
        t = np.log10(yy[idx]) if spec["log"] else yy[idx]
        rec: Dict = {"n_train": int(len(idx))}
        for name, build_m in _models().items():
            m = build_m()
            t0 = time.perf_counter()
            m.fit(X[idx], t)
            fit_s = time.perf_counter() - t0
            blob = pickle.dumps(m)
            one = X[idx][:1]
            m.predict(one)  # warm any lazy allocation
            lat = []
            for _ in range(200):
                t0 = time.perf_counter()
                m.predict(one)
                lat.append(time.perf_counter() - t0)
            rec[name] = {
                "train_s": round(fit_s, 3),
                "pickle_bytes": len(blob),
                "single_row_predict_ms_median": round(
                    float(np.median(lat)) * 1e3, 4),
                "single_row_predict_ms_p90": round(
                    float(np.percentile(lat, 90)) * 1e3, 4),
            }
        out[target] = rec
    return out


def main() -> int:
    with open(_IN) as fh:
        blob = json.load(fh)
    rows = [r for r in blob["rows"] if r.get("status") == "ok"]
    print(f"[fit] {len(rows)} ok rows of {len(blob['rows'])}; "
          f"alpha={blob.get('alpha')}")
    X, names, y, meta = build(rows)
    print(f"[fit] X = {X.shape}, {len(set(meta['donor'].tolist()))} donors, "
          f"{len(set(meta['section'].tolist()))} sections, "
          f"{len(set(meta['drug'].tolist()))} compounds")
    print(f"[fit] spatial_spread_valid rows: {int(meta['valid'].sum())} "
          f"({len(set(meta['section'][meta['valid']].tolist()))} sections)")

    out = {"alpha": blob.get("alpha"), "dx_um": blob.get("dx_um"),
           "n_rows_ok": len(rows), "feature_names": names,
           "n_donors": len(set(meta["donor"].tolist())),
           "n_sections": len(set(meta["section"].tolist())),
           "n_compounds": len(set(meta["drug"].tolist())),
           "evals": []}

    subsets = {"all": None, "valid_only": meta["valid"]}
    for target in TARGETS:
        for split in ("random", "section", "donor", "compound"):
            for sname, sub in subsets.items():
                if sname == "valid_only" and split in ("section", "donor"):
                    nv = len(set(meta["donor"][meta["valid"]].tolist()))
                    if nv < 3:
                        continue
                try:
                    e = evaluate(X, y, meta, target, split, sub)
                except Exception as exc:  # noqa: BLE001
                    print(f"  !! {target}/{split}/{sname}: "
                          f"{type(exc).__name__}: {exc}")
                    continue
                e["subset"] = sname
                out["evals"].append(e)
                mm = e["models"]
                fs = e["across_group_fold_spread"]
                print(f"\n{target} [{split}/{sname}] n={e['n']} "
                      f"folds={e['n_folds']}  "
                      f"y in [{e['y_min']:.4g}, {e['y_max']:.4g}] "
                      f"across-group spread "
                      f"{('%.3g' % fs) if fs else 'n/a'}x "
                      f"(range {e['across_group_abs_range']:.4g} "
                      f"{e['unit'].split(' ')[0]})")
                for mname in ("mean", "analytic", "ridge", "hgb"):
                    r = mm.get(mname)
                    if not r or not r.get("n_scored"):
                        continue
                    print(f"    {mname:9s} MAE={r['mae']:10.4g} "
                          f"medAE={r['median_ae']:10.4g} "
                          f"p90AE={r['p90_ae']:10.4g} "
                          f"foldErr(med)={r['median_fold_error']:7.3f} "
                          f"R2={r['r2']:8.3f} "
                          f"skill_vs_mean="
                          f"{r.get('skill_vs_mean', float('nan')):7.3f}")

    out["ablation_donor_split"] = ablate(X, y, meta, names, "donor")
    print("\n[ablation] leave-one-DONOR-out, skill vs the mean baseline:")
    for target, rec in out["ablation_donor_split"].items():
        print(f"  {target}")
        for cname, r in rec.items():
            if "error" in r:
                print(f"      {cname:14s} {r['error'][:70]}")
                continue
            print(f"      {cname:14s} nfeat={r['n_features']:3d} "
                  f"ridge_skill="
                  f"{(r.get('ridge') or {}).get('skill_vs_mean')} "
                  f"hgb_skill={(r.get('hgb') or {}).get('skill_vs_mean')}")

    out["final_model_cost"] = cost_of_final_models(X, y, meta)
    print("\n[cost] final models (trained on every row):")
    for target, rec in out["final_model_cost"].items():
        for name in ("ridge", "hgb"):
            c = rec.get(name) or {}
            print(f"    {target:24s} {name:6s} n={rec['n_train']:5d} "
                  f"train={c.get('train_s')}s "
                  f"size={c.get('pickle_bytes')}B "
                  f"predict1={c.get('single_row_predict_ms_median')}ms")

    with open(_OUT, "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print(f"\n[fit] -> {_OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
