"""Stage 0 kill-test harness — the honest, adversarial experiment.

This module runs the single question the whole scaffold exists to answer
(``ARCHITECTURE_CONTRACT.md``):

    *Does a mechanistic microenvironment + drug-transport simulation predict
    observed drug response BETTER than a cheap non-spatial baseline?*

It is deliberately written to be able to say **no**. If the simulation does not
clear the best baseline on the primary metric by a real margin, the ``verdict``
says so in plain language and ``beat_baseline`` is ``False``. A loss is not spun
as a win — a failed kill-test is a valid (and useful) Stage 0 outcome.

Fair-comparison design
----------------------
Every model is scored on the SAME pooled per-patient predictions against the
SAME observed labels:

* The **simulation** and the **fit-free baselines** (``mean``,
  ``intrinsic_sensitivity``) need no training, so they predict every patient
  directly.
* The one **fit-requiring baseline** (``expression_logreg``) is evaluated with
  K-fold cross-validation (``K=5``, ``sklearn.KFold(random_state=0)``) producing
  **out-of-fold** predictions, so it is never scored on data it trained on. This
  is the strongest honest version of the ML competitor — handicapping it would
  make a sim "win" meaningless.

The pass rule (``beat_baseline``) is intentionally conservative: the sim must
*exceed* the best baseline's ROC-AUC by a margin (``>= 0.02``), not merely tie.
"""
from __future__ import annotations

import json
import os
from typing import Dict, List, Optional, Sequence

import numpy as np
from sklearn.model_selection import KFold

from tmesim.core import ExperimentResult, SimConfig
from tmesim.validation.baselines import all_baselines
from tmesim.validation.metrics import score_all

__all__ = ["run_stage0_experiment", "PASS_MARGIN", "N_SPLITS"]

# Minimum ROC-AUC advantage the sim must have over the BEST baseline to "pass".
PASS_MARGIN: float = 0.02
# Cross-validation folds for the fit-requiring baselines.
N_SPLITS: int = 5


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _observed(records: Sequence) -> np.ndarray:
    """Ground-truth labels aligned to ``records`` as a float array."""
    return np.asarray([float(r.observed_response) for r in records], dtype=np.float64)


def _cross_val_predict(fn, records: List) -> np.ndarray:
    """Out-of-fold predictions for a fit-requiring baseline.

    Runs ``KFold(n_splits=min(N_SPLITS, n), shuffle=True, random_state=0)`` over
    the cohort; for each fold the baseline is fit on the train records and used
    to predict the held-out test records, whose predictions are written back into
    their original positions. No patient is ever scored on a model trained on it.

    Falls back gracefully for cohorts too small to split (``n < 2``): the
    baseline is fit and predicted on the whole cohort (degenerate, but never
    crashes), which the honest metrics will treat as leaky/optimistic — such a
    cohort is too small to conclude anything anyway.

    Parameters
    ----------
    fn : Callable[[list, list], np.ndarray]
        A ``needs_split`` baseline: ``fn(train, test) -> predictions``.
    records : list[PatientRecord]

    Returns
    -------
    np.ndarray
        Out-of-fold predictions aligned to ``records``.
    """
    n = len(records)
    oof = np.full(n, np.nan, dtype=np.float64)
    if n == 0:
        return oof
    if n < 2:
        # Cannot cross-validate a single sample; degenerate fit-on-all.
        oof[:] = np.asarray(fn(records, records), dtype=np.float64)
        return oof

    k = min(N_SPLITS, n)
    kf = KFold(n_splits=k, shuffle=True, random_state=0)
    idx = np.arange(n)
    for train_idx, test_idx in kf.split(idx):
        train = [records[i] for i in train_idx]
        test = [records[i] for i in test_idx]
        pred = np.asarray(fn(train, test), dtype=np.float64)
        oof[test_idx] = pred
    return oof


def _simulate_cohort(records: List, cfg: SimConfig) -> np.ndarray:
    """Run the full mechanistic simulation for every patient.

    ``pred = simulate_response(record.micro, cfg).predicted_response`` for each
    record, aligned to ``records``. Imported lazily so this module still imports
    even while ``tmesim.sim.engine`` is landing concurrently.

    Parameters
    ----------
    records : list[PatientRecord]
    cfg : SimConfig
        The simulation configuration (its ``seed`` makes the run deterministic).

    Returns
    -------
    np.ndarray
        Predicted responses in ``[0, 1]`` aligned to ``records``.
    """
    from tmesim.sim.engine import simulate_response  # lazy: engine may be WIP

    preds = np.empty(len(records), dtype=np.float64)
    for i, r in enumerate(records):
        result = simulate_response(r.micro, cfg)
        preds[i] = float(result.predicted_response)
    return preds


def _to_jsonable(obj):
    """Recursively convert numpy scalars/arrays to JSON-serializable Python."""
    if isinstance(obj, dict):
        return {str(k): _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _to_jsonable(obj.tolist())
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, float):
        return obj if np.isfinite(obj) else None
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


# --------------------------------------------------------------------------- #
# The experiment
# --------------------------------------------------------------------------- #
def run_stage0_experiment(
    records: List, cfg: SimConfig, out_dir: str = "results"
) -> ExperimentResult:
    """Run the Stage 0 kill-test across a cohort and render an honest verdict.

    For every patient the mechanistic simulation predicts a response, which is
    scored against the observed outcome alongside each baseline. The
    fit-requiring baseline is cross-validated (out-of-fold); the sim and fit-free
    baselines predict directly. Everything is scored on the same pooled
    predictions with :func:`~tmesim.validation.metrics.score_all`.

    The sim "passes" (``beat_baseline=True``) only if it EXCEEDS the best
    baseline's ROC-AUC by at least :data:`PASS_MARGIN`. The ``verdict`` string
    states the outcome plainly whether it is a win or a loss.

    Parameters
    ----------
    records : list[PatientRecord]
        The cohort. Each must carry ``micro``, ``observed_response`` and
        ``features``.
    cfg : SimConfig
        Simulation configuration (deterministic via ``cfg.seed``).
    out_dir : str, default "results"
        Directory to write ``stage0_result.json`` into (created if missing).

    Returns
    -------
    ExperimentResult
        Fully populated: ``model_scores`` (sim), ``baseline_scores`` (per
        baseline), ``per_patient`` rows, ``primary_metric='roc_auc'``,
        ``verdict`` and ``beat_baseline``.

    Raises
    ------
    ValueError
        If ``records`` is empty (nothing to test).
    """
    if not records:
        raise ValueError("run_stage0_experiment: empty cohort — nothing to test.")

    primary_metric = "roc_auc"
    y_true = _observed(records)

    # --- 1. simulation predictions ---------------------------------------- #
    y_sim = _simulate_cohort(records, cfg)

    # --- 2. baseline predictions (fit-free direct; fit-needing cross-val) -- #
    baseline_preds: Dict[str, np.ndarray] = {}
    for name, spec in all_baselines().items():
        if spec.needs_split:
            baseline_preds[name] = _cross_val_predict(spec.fn, records)
        else:
            baseline_preds[name] = np.asarray(spec.fn(records), dtype=np.float64)

    # --- 3. score everything on the same labels --------------------------- #
    model_scores = score_all(y_true, y_sim)
    baseline_scores: Dict[str, Dict[str, float]] = {
        name: score_all(y_true, pred) for name, pred in baseline_preds.items()
    }

    # --- 4. verdict: does the sim beat the BEST baseline by the margin? ---- #
    sim_auc = float(model_scores[primary_metric])
    best_name = max(baseline_scores, key=lambda k: baseline_scores[k][primary_metric])
    best_auc = float(baseline_scores[best_name][primary_metric])
    delta = sim_auc - best_auc
    beat_baseline = bool(delta >= PASS_MARGIN)

    if beat_baseline:
        verdict = (
            f"Sim AUC {sim_auc:.2f} vs best baseline ({best_name}) {best_auc:.2f} "
            f"-> BEATS baseline by {delta:+.2f} (>= {PASS_MARGIN:.2f} margin); "
            f"mechanistic transport modeling adds measurable value on this cohort "
            f"(n={len(records)})."
        )
    else:
        verdict = (
            f"Sim AUC {sim_auc:.2f} vs best baseline ({best_name}) {best_auc:.2f} "
            f"-> DOES NOT beat baseline (delta {delta:+.2f} < {PASS_MARGIN:.2f} "
            f"margin); transport modeling adds no measurable value beyond the "
            f"baseline on this cohort (n={len(records)})."
        )

    # --- 5. per-patient rows (also carry y_true & sim score for the ROC) --- #
    per_patient: List[Dict] = []
    for i, r in enumerate(records):
        row: Dict[str, object] = {
            "patient_id": r.patient_id,
            "observed": float(y_true[i]),
            "sim_pred": float(y_sim[i]),
        }
        for name, pred in baseline_preds.items():
            row[name] = float(pred[i]) if np.isfinite(pred[i]) else None
        per_patient.append(row)

    result = ExperimentResult(
        model_scores=model_scores,
        baseline_scores=baseline_scores,
        per_patient=per_patient,
        primary_metric=primary_metric,
        verdict=verdict,
        beat_baseline=beat_baseline,
    )

    # --- 6. persist an honest, JSON-serializable summary ------------------ #
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "stage0_result.json")
    summary = {
        "primary_metric": primary_metric,
        "beat_baseline": beat_baseline,
        "verdict": verdict,
        "pass_margin": PASS_MARGIN,
        "n_patients": len(records),
        "cv_folds": min(N_SPLITS, len(records)),
        "best_baseline": best_name,
        "model_scores": model_scores,
        "baseline_scores": baseline_scores,
        "per_patient": per_patient,
    }
    with open(out_path, "w") as fh:
        json.dump(_to_jsonable(summary), fh, indent=2)

    return result


# --------------------------------------------------------------------------- #
# Self-test / smoke run
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    import copy

    from tmesim.data.synthetic import generate_synthetic_cohort

    print("Stage 0 harness self-test: building 24-patient synthetic GBM cohort...")
    cohort = generate_synthetic_cohort(n=24, seed=0)

    # Small/fast config: shrink the PDE window so 24 patients finish quickly.
    # (The transport solver is implicit/unconditionally stable, so a short
    # window is a valid smoke test — it just samples less exposure.)
    cfg = SimConfig()
    cfg = copy.deepcopy(cfg)
    cfg.pde.total_time = 600.0   # seconds (default is 12 h) -> fast
    cfg.pde.save_every = 200

    out_dir = "results"
    exp = run_stage0_experiment(cohort, cfg, out_dir=out_dir)

    print("\n=== STAGE 0 VERDICT ===")
    print(exp.verdict)
    print(f"beat_baseline = {exp.beat_baseline}")

    print("\n--- Simulation scores ---")
    for k, v in exp.model_scores.items():
        print(f"  sim.{k:9s} = {v:.4f}")
    print("--- Baseline scores ---")
    for name, sc in exp.baseline_scores.items():
        print(f"  {name}:")
        for k, v in sc.items():
            print(f"      {k:9s} = {v:.4f}")

    # --- assert the result is well-formed -------------------------------- #
    assert exp.verdict.strip(), "verdict must be non-empty"
    assert exp.primary_metric == "roc_auc"
    for sc in [exp.model_scores, *exp.baseline_scores.values()]:
        assert np.isfinite(sc["roc_auc"]), "roc_auc must be finite"
        assert np.isfinite(sc["spearman"]), "spearman must be finite"
    assert len(exp.per_patient) == len(cohort)
    for row in exp.per_patient:
        assert np.isfinite(row["observed"]) and np.isfinite(row["sim_pred"])

    json_path = os.path.join(out_dir, "stage0_result.json")
    assert os.path.exists(json_path), "stage0_result.json must be written"
    with open(json_path) as fh:
        loaded = json.load(fh)
    assert loaded["verdict"] == exp.verdict
    print(f"\n[OK] well-formed; json written to {json_path}")
