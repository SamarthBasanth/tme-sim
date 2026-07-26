"""Scoring metrics for the Stage 0 kill-test validation harness.

These functions score a set of predicted responses (``y_score``, typically in
[0, 1]) against observed outcomes (``y_true``, either binary responder labels or
a continuous sensitivity such as ``-log(IC50)``). The harness deliberately has to
handle *both* label kinds and a lot of degenerate cohorts (all responders, all
non-responders, a baseline that emits a constant), so every function here is
written to be **nan-safe and never crash** — a broken metric must not be able to
flatter or sink the verdict by throwing.

We prefer scikit-learn / scipy implementations but guard their edge cases
explicitly, because their default behaviour on degenerate input is to raise
(e.g. ``roc_auc_score`` with a single class) rather than return a neutral value.

Conventions
-----------
- Higher ``y_score`` means "more likely to respond / more sensitive".
- ``roc_auc`` and ``brier`` need a *binary* ``y_true``. If a continuous
  ``y_true`` is passed, it is thresholded at its own median (values strictly
  greater than the median become the positive class). This is documented on each
  function and applied consistently.
- ``spearman`` is the metric for genuinely continuous outcomes and never
  binarizes.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import brier_score_loss, roc_auc_score

__all__ = [
    "roc_auc",
    "brier",
    "calibration_curve",
    "spearman",
    "score_all",
]


def _as_1d_float(a) -> np.ndarray:
    """Coerce input to a flat float64 array (copy), for defensive arithmetic."""
    return np.asarray(a, dtype=float).ravel()


def _binarize(y_true: np.ndarray) -> np.ndarray:
    """Return ``y_true`` as a {0, 1} array.

    If ``y_true`` already contains only 0/1 values it is returned unchanged
    (cast to int). Otherwise it is treated as continuous and thresholded at its
    median: ``value > median -> 1`` (positive), else 0. Thresholding strictly
    above the median keeps the split stable and avoids assigning the whole
    median-tie group to the positive class.
    """
    y = _as_1d_float(y_true)
    uniq = np.unique(y[np.isfinite(y)])
    if uniq.size <= 2 and np.all(np.isin(uniq, (0.0, 1.0))):
        return (y > 0.5).astype(int)
    med = np.nanmedian(y)
    return (y > med).astype(int)


def roc_auc(y_true, y_score) -> float:
    """Area under the ROC curve for a binary outcome.

    ``y_true`` may be binary (0/1) or continuous; a continuous ``y_true`` is
    thresholded at its **median** (value > median => positive class) before
    scoring — see :func:`_binarize`.

    Degenerate cases return the chance value **0.5** rather than raising:
    - ``y_true`` has only one class after binarizing (AUC undefined),
    - ``y_score`` is constant,
    - fewer than 2 finite paired samples.

    Parameters
    ----------
    y_true : array-like
        Observed outcomes (binary or continuous).
    y_score : array-like
        Predicted responses / scores (higher = more likely positive).

    Returns
    -------
    float
        ROC AUC in [0, 1], or 0.5 for degenerate input.
    """
    yt = _as_1d_float(y_true)
    ys = _as_1d_float(y_score)
    mask = np.isfinite(yt) & np.isfinite(ys)
    if mask.sum() < 2:
        return 0.5
    yb = _binarize(yt[mask])
    ss = ys[mask]
    if np.unique(yb).size < 2:
        return 0.5
    if np.ptp(ss) == 0:  # constant score -> chance
        return 0.5
    try:
        return float(roc_auc_score(yb, ss))
    except ValueError:
        return 0.5


def brier(y_true, y_score) -> float:
    """Brier score (mean squared error of probabilistic predictions).

    Lower is better; 0 is perfect. ``y_score`` is clipped to [0, 1] and treated
    as a probability. ``y_true`` may be binary or continuous; a continuous
    ``y_true`` is thresholded at its median (see :func:`_binarize`).

    Returns ``nan`` only if there are no finite paired samples; otherwise always
    returns a finite value and never raises.

    Parameters
    ----------
    y_true : array-like
        Observed outcomes (binary or continuous).
    y_score : array-like
        Predicted probabilities (clipped to [0, 1]).

    Returns
    -------
    float
        Brier score in [0, 1], or ``nan`` if no usable samples.
    """
    yt = _as_1d_float(y_true)
    ys = _as_1d_float(y_score)
    mask = np.isfinite(yt) & np.isfinite(ys)
    if mask.sum() == 0:
        return float("nan")
    yb = _binarize(yt[mask])
    ss = np.clip(ys[mask], 0.0, 1.0)
    try:
        # Pin the label set so a single-class y_true doesn't confuse sklearn.
        return float(brier_score_loss(yb, ss, pos_label=1))
    except (ValueError, TypeError):
        # Manual fallback: MSE against the binary labels.
        return float(np.mean((ss - yb) ** 2))


def calibration_curve(y_true, y_score, bins: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """Reliability-diagram data: mean predicted vs observed frequency per bin.

    Predictions are bucketed into ``bins`` equal-width bins over [0, 1]. For each
    non-empty bin we report the mean predicted score and the observed positive
    frequency. Empty bins are dropped, so the two returned arrays are aligned and
    have the same (<= ``bins``) length.

    ``y_true`` may be binary or continuous (thresholded at its median, see
    :func:`_binarize`). Never raises; returns two empty arrays if there are no
    usable samples.

    Parameters
    ----------
    y_true : array-like
        Observed outcomes.
    y_score : array-like
        Predicted probabilities (clipped to [0, 1] for binning).
    bins : int, default 10
        Number of equal-width probability bins.

    Returns
    -------
    (mean_predicted, observed_frequency) : tuple of np.ndarray
        Per non-empty bin, aligned. Empty if no usable samples.
    """
    yt = _as_1d_float(y_true)
    ys = _as_1d_float(y_score)
    mask = np.isfinite(yt) & np.isfinite(ys)
    if mask.sum() == 0 or bins < 1:
        return np.array([]), np.array([])
    yb = _binarize(yt[mask])
    ss = np.clip(ys[mask], 0.0, 1.0)

    edges = np.linspace(0.0, 1.0, bins + 1)
    # np.digitize with the interior edges yields bin indices 0..bins-1;
    # clip guards the value exactly at 1.0 into the last bin.
    idx = np.clip(np.digitize(ss, edges[1:-1], right=False), 0, bins - 1)

    mean_pred = []
    obs_freq = []
    for b in range(bins):
        sel = idx == b
        if not np.any(sel):
            continue
        mean_pred.append(ss[sel].mean())
        obs_freq.append(yb[sel].mean())
    return np.asarray(mean_pred), np.asarray(obs_freq)


def spearman(y_true, y_score) -> float:
    """Spearman rank correlation between observed and predicted, for continuous
    outcomes.

    This is the metric for genuinely continuous ``y_true`` (e.g. ``-log(IC50)``)
    and never binarizes anything. Higher = better monotone agreement.

    Degenerate cases return **0.0** (no detectable monotone relationship) rather
    than raising or returning nan:
    - fewer than 2 finite paired samples,
    - either input constant (rank correlation undefined),
    - scipy returns nan.

    Parameters
    ----------
    y_true : array-like
        Observed outcomes.
    y_score : array-like
        Predicted responses / scores.

    Returns
    -------
    float
        Spearman rho in [-1, 1], or 0.0 for degenerate input.
    """
    yt = _as_1d_float(y_true)
    ys = _as_1d_float(y_score)
    mask = np.isfinite(yt) & np.isfinite(ys)
    if mask.sum() < 2:
        return 0.0
    a = yt[mask]
    b = ys[mask]
    if np.ptp(a) == 0 or np.ptp(b) == 0:  # a constant vector has no rank order
        return 0.0
    rho = spearmanr(a, b).correlation
    if rho is None or not np.isfinite(rho):
        return 0.0
    return float(rho)


def score_all(y_true, y_score) -> dict:
    """Score predictions with all headline metrics, robustly.

    Returns a dict ``{'roc_auc', 'brier', 'spearman'}``. Every degenerate case
    (single-class ``y_true``, constant ``y_score``, empty overlap, nans) is
    absorbed by the individual metrics, so this never raises:
    - ``roc_auc`` -> 0.5 on degenerate input,
    - ``brier``  -> nan only when there are no usable samples,
    - ``spearman`` -> 0.0 on degenerate input.

    Parameters
    ----------
    y_true : array-like
        Observed outcomes (binary or continuous).
    y_score : array-like
        Predicted responses / scores.

    Returns
    -------
    dict
        {'roc_auc': float, 'brier': float, 'spearman': float}
    """
    return {
        "roc_auc": roc_auc(y_true, y_score),
        "brier": brier(y_true, y_score),
        "spearman": spearman(y_true, y_score),
    }


if __name__ == "__main__":
    rng = np.random.default_rng(0)

    # --- Known-answer checks -------------------------------------------------
    # 1. Perfect separation -> AUC 1.0, spearman 1.0, low brier.
    y_true_bin = np.array([0, 0, 0, 1, 1, 1])
    y_perfect = np.array([0.01, 0.02, 0.10, 0.90, 0.95, 0.99])
    auc_perfect = roc_auc(y_true_bin, y_perfect)
    assert abs(auc_perfect - 1.0) < 1e-9, auc_perfect
    assert brier(y_true_bin, y_perfect) < 0.02
    print(f"[perfect]   auc={auc_perfect:.3f}  brier={brier(y_true_bin, y_perfect):.4f}")

    # 2. Random scores on a balanced label -> AUC ~ 0.5.
    yt = (rng.random(2000) > 0.5).astype(int)
    ys = rng.random(2000)
    auc_rand = roc_auc(yt, ys)
    assert 0.4 < auc_rand < 0.6, auc_rand
    print(f"[random]    auc={auc_rand:.3f}  (expect ~0.5)")

    # 3. Constant score -> exactly 0.5 (degenerate guard).
    auc_const = roc_auc(yt, np.full_like(ys, 0.7))
    assert auc_const == 0.5, auc_const
    print(f"[constant]  auc={auc_const:.3f}  (expect 0.5)")

    # 4. Monotone continuous relationship -> spearman near 1.
    x = np.linspace(0, 10, 100)
    y_mono = 2.0 * x + 3.0 + rng.normal(0, 0.05, x.size)
    rho = spearman(x, y_mono)
    assert rho > 0.99, rho
    print(f"[monotone]  spearman={rho:.3f}  (expect ~1.0)")

    # 5. Perfectly anti-correlated -> spearman -1.
    rho_neg = spearman(x, -x)
    assert abs(rho_neg + 1.0) < 1e-9, rho_neg
    print(f"[anti]      spearman={rho_neg:.3f}  (expect -1.0)")

    # 6. Continuous y_true thresholded at median for AUC (perfect ranking).
    y_cont = np.array([-2.0, -1.0, 0.0, 1.0, 2.0, 3.0])  # median 0.5
    auc_cont = roc_auc(y_cont, y_cont)  # score == truth -> perfect ranking
    assert abs(auc_cont - 1.0) < 1e-9, auc_cont
    print(f"[cont->bin] auc={auc_cont:.3f}  (median-thresholded, expect 1.0)")

    # 7. Calibration curve: well-calibrated data lands near the diagonal.
    p = rng.random(20000)
    labels = (rng.random(20000) < p).astype(int)
    mean_pred, obs_freq = calibration_curve(labels, p, bins=10)
    max_gap = float(np.max(np.abs(mean_pred - obs_freq)))
    assert max_gap < 0.05, max_gap
    print(f"[calib]     bins={len(mean_pred)}  max|pred-obs|={max_gap:.4f}  (expect <0.05)")

    # 8. Degenerate cases must never crash.
    single = score_all(np.ones(10), rng.random(10))          # single-class y_true
    empty = score_all(np.array([]), np.array([]))            # no samples
    withnan = score_all([0, 1, np.nan, 1], [0.2, 0.8, 0.5, np.nan])
    assert single["roc_auc"] == 0.5 and single["spearman"] == 0.0
    assert empty["roc_auc"] == 0.5 and np.isnan(empty["brier"])
    assert np.isfinite(withnan["roc_auc"])
    print(f"[degenerate] single={single}")
    print(f"[degenerate] empty ={empty}")
    print(f"[degenerate] w/nan ={withnan}")

    # 9. score_all on the perfect case.
    print(f"[score_all] {score_all(y_true_bin, y_perfect)}")

    print("\nAll known-answer checks passed.")
