"""Shared, honest evaluation for the drug-selection sweep.

An agent supplies only a ground-truth kill matrix gt[P,D] (its independent
biological model of how well each drug works in each patient). This module scores
the fixed policies against it identically for everyone:

  - sim         : mechanistic PDE-sim prediction (precomputed, no training)
  - potency     : delivery-agnostic potency shortcut (precomputed, no training)
  - feature_ml  : cheap cross-validated classifier on OBSERVABLE patient features
                  (the real competitor -- if it matches the sim, physics adds nothing)
  - random/oracle references

Metrics: top-1 accuracy (picked the true-best drug) and mean regret (true kill
lost by the pick). Plus diagnostics on how delivery-limited the ground truth is,
so the synthesis can test whether the sim's advantage is real or circular.
"""
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

_CACHE = None


def load_cache(path="results/selection_cache.npz"):
    global _CACHE
    if _CACHE is None:
        z = np.load(path, allow_pickle=True)
        _CACHE = {k: z[k] for k in z.files}
    return _CACHE


def _feature_ml_picks(feat, true_best, seed=0):
    """Out-of-fold multiclass predictions of the best drug from patient features."""
    P = len(true_best)
    picks = np.array(true_best, copy=True)
    classes, counts = np.unique(true_best, return_counts=True)
    if len(classes) < 2 or counts.min() < 2:
        # degenerate: not enough per-class samples to CV -> majority class
        return np.full(P, classes[np.argmax(counts)])
    n_splits = int(min(5, counts.min()))
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for tr, te in skf.split(feat, true_best):
        sc = StandardScaler().fit(feat[tr])
        clf = LogisticRegression(max_iter=500, C=1.0)
        try:
            clf.fit(sc.transform(feat[tr]), true_best[tr])
            picks[te] = clf.predict(sc.transform(feat[te]))
        except Exception:
            picks[te] = np.bincount(true_best[tr]).argmax()
    return picks


def score_policies(gt_kill, seed=0, cache_path="results/selection_cache.npz"):
    """gt_kill: (P,D) ground-truth kill for each patient x drug. Returns metrics."""
    c = load_cache(cache_path)
    K_sim, potency, feat = c["K_sim"], c["potency"], c["feat"]
    drug_ic50 = c["drug_ic50"]
    P, D = gt_kill.shape
    true_best = gt_kill.argmax(1)

    def acc_reg(picks):
        picks = np.asarray(picks)
        acc = float(np.mean(picks == true_best))
        regret = float(np.mean(gt_kill[np.arange(P), true_best] - gt_kill[np.arange(P), picks]))
        return acc, regret

    sim_pick = K_sim.argmax(1)
    pot_pick = potency.argmax(1)
    ml_pick = _feature_ml_picks(feat, true_best, seed=seed)
    most_potent = int(np.argmin(drug_ic50))
    rand_regret = float(np.mean(gt_kill[np.arange(P), true_best] - gt_kill.mean(1)))

    a_sim, r_sim = acc_reg(sim_pick)
    a_pot, r_pot = acc_reg(pot_pick)
    a_ml, r_ml = acc_reg(ml_pick)

    # diagnostics: is this ground truth delivery-limited (does the best drug differ
    # from the globally-most-potent drug)? and does sim's ranking track the GT?
    frac_not_most_potent = float(np.mean(true_best != most_potent))
    # per-patient Spearman-ish: mean rank correlation of sim scores vs gt
    from scipy.stats import spearmanr
    rho = np.nanmean([spearmanr(K_sim[p], gt_kill[p]).correlation for p in range(P)])

    return {
        "acc": {"sim": a_sim, "potency": a_pot, "feature_ml": a_ml,
                "random": 1.0 / D, "oracle": 1.0},
        "regret": {"sim": r_sim, "potency": r_pot, "feature_ml": r_ml,
                   "random": rand_regret, "oracle": 0.0},
        "diag": {"frac_best_not_most_potent": frac_not_most_potent,
                 "sim_vs_gt_rank_corr": float(rho),
                 "best_distribution": [int((true_best == d).sum()) for d in range(D)],
                 "P": P, "D": D},
        "sim_beats_ml": bool(a_sim > a_ml and r_sim < r_ml),
        "sim_beats_potency": bool(a_sim > a_pot and r_sim < r_pot),
    }


def print_json(result):
    import json
    print("RESULT_JSON " + json.dumps(result, default=float))
