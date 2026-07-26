"""Cheap, non-spatial baselines the mechanistic simulation must beat.

These are the *null / competitor* models for the Stage 0 kill-test
(``ARCHITECTURE_CONTRACT.md``). The whole point of the mechanistic
microenvironment + drug-transport simulation is to predict observed drug
response *better than something cheap*. If it cannot, the idea is not ready.
So these baselines are built to be **rigorous and fair** — each is the strongest
honest version of its idea, given only non-spatial information:

- ``baseline_mean`` — a proper null model: predict the cohort mean response for
  everyone. Any model that cannot beat this has learned nothing.
- ``baseline_intrinsic_sensitivity`` — the KEY control. It uses ONLY the mean
  per-cell intrinsic tumor ``sensitivity`` (MGMT status, resistance clones, etc.
  already baked into the biology) with **no transport / spatial physics at all**.
  If the full simulation cannot beat this, then the expensive spatial
  drug-transport modelling is adding no value beyond intrinsic sensitivity.
- ``baseline_expression_logreg`` — a standard ML baseline: a (regularized)
  logistic / linear regression on the flat non-spatial ``features`` dict, with
  proper standardization and train/test separation.

All predictions are returned aligned (same order, same length) to the records
passed in, in ``[0, 1]``, and are fully deterministic (``random_state=0``).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Sequence

import numpy as np
from sklearn.linear_model import LinearRegression, LogisticRegression, Ridge
from sklearn.preprocessing import StandardScaler

# PatientRecord is a plain dataclass in tmesim.core.results; we only touch
# attributes (patient_id, micro, observed_response, features) so we avoid a hard
# import to keep this module light and testable in isolation.

__all__ = [
    "baseline_mean",
    "baseline_intrinsic_sensitivity",
    "baseline_expression_logreg",
    "all_baselines",
    "BaselineSpec",
]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _observed(records: Sequence) -> np.ndarray:
    """Ground-truth labels as a float array aligned to ``records``."""
    return np.asarray([float(r.observed_response) for r in records], dtype=np.float64)


def _is_binary(y: np.ndarray) -> bool:
    """True if every label is exactly 0 or 1 (a classification outcome).

    Continuous outcomes (e.g. ``-log(IC50)``) fail this and are handled by
    regression instead.
    """
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return False
    return bool(np.all(np.isin(finite, (0.0, 1.0))))


def _mean_tumor_sensitivity(record) -> float:
    """Mean intrinsic per-cell ``sensitivity`` over this patient's TUMOR cells.

    No transport, no spatial physics — just the biology encoded per cell.
    Falls back to all cells (then to 0.5) if a record has no tumour cells, so
    the baseline never emits NaN.
    """
    cells = record.micro.cells
    sens = np.asarray(cells.sensitivity, dtype=np.float64)
    if sens.size == 0:
        return 0.5
    tum = np.asarray(cells.tumor_mask())
    if tum.any():
        return float(np.nanmean(sens[tum]))
    return float(np.nanmean(sens))


def _numeric_feature_keys(records: Sequence) -> List[str]:
    """Sorted union of feature keys that carry a numeric (non-None) value.

    Sorting makes the feature ordering deterministic and stable across records,
    which is what lets us align train/test matrices column-for-column.
    """
    keys = set()
    for r in records:
        for k, v in (r.features or {}).items():
            if v is None:
                continue
            if isinstance(v, (bool, int, float, np.integer, np.floating)):
                keys.add(k)
    return sorted(keys)


def _feature_matrix(
    records: Sequence, keys: Sequence[str], fill: Dict[str, float]
) -> np.ndarray:
    """Build an ``(n_records, len(keys))`` matrix, imputing missing keys.

    Robust to records that lack some keys entirely: a missing/None/non-numeric
    entry is filled with ``fill[key]`` (the train-column mean), so a feature that
    is present for some patients but not others still contributes.
    """
    n = len(records)
    X = np.empty((n, len(keys)), dtype=np.float64)
    for i, r in enumerate(records):
        feats = r.features or {}
        for j, k in enumerate(keys):
            v = feats.get(k, None)
            if isinstance(v, (bool, int, float, np.integer, np.floating)) and v is not None:
                X[i, j] = float(v)
            else:
                X[i, j] = fill.get(k, 0.0)
    return X


# --------------------------------------------------------------------------- #
# Baselines
# --------------------------------------------------------------------------- #
def baseline_mean(records: List) -> np.ndarray:
    """Null model: predict the cohort mean observed response for every patient.

    This is the honest floor. A model with no skill scores exactly here on
    ranking metrics (constant score) and this is the best possible *constant*
    Brier prediction. Anything the simulation does must clear this.

    Parameters
    ----------
    records : list[PatientRecord]
        The cohort to predict.

    Returns
    -------
    np.ndarray
        Shape ``(len(records),)``, a constant equal to the mean observed
        response, in ``[0, 1]`` (observed responses are assumed in ``[0, 1]``).
    """
    n = len(records)
    if n == 0:
        return np.zeros(0, dtype=np.float64)
    mean = float(np.nanmean(_observed(records)))
    if not np.isfinite(mean):
        mean = 0.5
    return np.full(n, mean, dtype=np.float64)


def baseline_intrinsic_sensitivity(records: List) -> np.ndarray:
    """KEY baseline: mean per-cell tumor ``sensitivity`` only — no transport.

    For each patient we average the intrinsic drug ``sensitivity`` over its
    tumour cells (the biology a pure transport model cannot see: MGMT status,
    resistance clones, etc.), then min-max normalize the cohort into ``[0, 1]``
    so the baseline uses its full dynamic range. This isolates the question the
    whole project turns on: **does spatial drug-TRANSPORT modelling add anything
    beyond intrinsic sensitivity?** If the mechanistic sim cannot beat this, the
    answer is no.

    The normalization is monotone (rank-preserving), so it does not change
    ranking metrics (ROC-AUC / Spearman) versus the raw mean sensitivity — it
    only rescales for probability-style metrics. That keeps the baseline fair,
    not handicapped.

    Parameters
    ----------
    records : list[PatientRecord]

    Returns
    -------
    np.ndarray
        Shape ``(len(records),)`` in ``[0, 1]``.
    """
    n = len(records)
    if n == 0:
        return np.zeros(0, dtype=np.float64)
    raw = np.asarray([_mean_tumor_sensitivity(r) for r in records], dtype=np.float64)
    lo, hi = float(np.nanmin(raw)), float(np.nanmax(raw))
    if not np.isfinite(lo) or not np.isfinite(hi) or (hi - lo) < 1e-12:
        # Degenerate spread (or single patient): mean sensitivity is already in
        # [0,1], so just clip rather than divide by ~0.
        return np.clip(np.nan_to_num(raw, nan=0.5), 0.0, 1.0)
    return (raw - lo) / (hi - lo)


def baseline_expression_logreg(train: List, test: List) -> np.ndarray:
    """Standardized (regularized) regression on the flat ``features`` dict.

    A conventional non-spatial ML baseline. Fits on ``train`` and predicts on
    ``test``:

    - Features are the numeric entries of ``PatientRecord.features``, aligned by
      key across records; keys missing from a record are imputed with the train
      column mean, so heterogeneous feature dicts still work.
    - Features are standardized (fit on train only — no leakage from test).
    - Outcome type is detected from the *train* labels:
        * binary (all labels in {0,1})  -> ``LogisticRegression`` (L2),
          returns ``predict_proba`` of the positive class.
        * continuous                    -> ``Ridge`` (L2) with a
          ``LinearRegression`` fallback, prediction clipped to ``[0, 1]``.

    Deterministic (``random_state=0``). Robust to degenerate cases (no features,
    single-class train set) by falling back to the train mean response.

    Parameters
    ----------
    train, test : list[PatientRecord]

    Returns
    -------
    np.ndarray
        Shape ``(len(test),)`` in ``[0, 1]``, aligned to ``test``.
    """
    n_test = len(test)
    if n_test == 0:
        return np.zeros(0, dtype=np.float64)

    y_train = _observed(train)
    # Fallback prediction if we cannot fit anything meaningful.
    fallback_val = float(np.nanmean(y_train)) if len(train) else 0.5
    if not np.isfinite(fallback_val):
        fallback_val = 0.5
    fallback = np.full(n_test, np.clip(fallback_val, 0.0, 1.0), dtype=np.float64)

    if len(train) == 0:
        return fallback

    keys = _numeric_feature_keys(train)
    if not keys:
        # No usable features -> degenerate to the constant-mean prediction.
        return fallback

    # Train-column means for imputation (computed only over records that have
    # the key with a numeric value; default 0.0 if none do).
    fill: Dict[str, float] = {}
    for k in keys:
        vals = [
            float(r.features[k])
            for r in train
            if r.features
            and isinstance(r.features.get(k), (bool, int, float, np.integer, np.floating))
        ]
        fill[k] = float(np.mean(vals)) if vals else 0.0

    X_train = _feature_matrix(train, keys, fill)
    X_test = _feature_matrix(test, keys, fill)

    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    if _is_binary(y_train):
        classes = np.unique(y_train[np.isfinite(y_train)])
        if classes.size < 2:
            # Only one class present -> the honest prediction is that class.
            return np.full(n_test, float(classes[0]) if classes.size else fallback_val,
                           dtype=np.float64)
        model = LogisticRegression(random_state=0, max_iter=1000)
        model.fit(X_train_s, y_train.astype(int))
        pos_idx = int(np.where(model.classes_ == 1)[0][0])
        proba = model.predict_proba(X_test_s)[:, pos_idx]
        return np.clip(proba.astype(np.float64), 0.0, 1.0)

    # Continuous outcome -> regularized linear regression.
    try:
        model = Ridge(alpha=1.0, random_state=0)
        model.fit(X_train_s, y_train)
    except TypeError:
        # Older sklearn Ridge without random_state, or solver issue -> plain OLS.
        model = LinearRegression()
        model.fit(X_train_s, y_train)
    pred = model.predict(X_test_s).astype(np.float64)
    return np.clip(pred, 0.0, 1.0)


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class BaselineSpec:
    """Description of a baseline for the experiment harness to iterate over.

    Attributes
    ----------
    name : str
        Registry key / display name.
    fn : Callable
        Either ``fn(records) -> np.ndarray`` (when ``needs_split`` is False) or
        ``fn(train, test) -> np.ndarray`` (when ``needs_split`` is True).
    needs_split : bool
        Whether this baseline must be trained on a disjoint train set and
        evaluated on a held-out test set (i.e. cross-validated). Fit-free
        baselines (mean, intrinsic sensitivity) can score the whole cohort.
    description : str
        One-line rationale.
    """

    name: str
    fn: Callable
    needs_split: bool
    description: str


def all_baselines() -> Dict[str, BaselineSpec]:
    """Registry of baselines for the Stage 0 harness.

    Returns
    -------
    dict[str, BaselineSpec]
        Keyed by baseline name. The harness iterates this, checks
        ``spec.needs_split`` to decide whether to cross-validate, and calls
        ``spec.fn`` accordingly.
    """
    return {
        "mean": BaselineSpec(
            name="mean",
            fn=baseline_mean,
            needs_split=False,
            description="Constant cohort-mean response (null model).",
        ),
        "intrinsic_sensitivity": BaselineSpec(
            name="intrinsic_sensitivity",
            fn=baseline_intrinsic_sensitivity,
            needs_split=False,
            description="Mean per-cell tumor sensitivity only; no transport/spatial physics.",
        ),
        "expression_logreg": BaselineSpec(
            name="expression_logreg",
            fn=baseline_expression_logreg,
            needs_split=True,
            description="Standardized (logistic/linear) regression on non-spatial features.",
        ),
    }


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    # Build a tiny synthetic cohort of PatientRecords inline (no data.synthetic).
    from tmesim.core import (
        CellPopulation,
        CellType,
        CellState,
        Domain,
        Microenvironment,
        PatientRecord,
        TissueFields,
    )

    rng = np.random.default_rng(0)

    def _make_record(pid: str, base_sens: float, resp: float) -> PatientRecord:
        """Fabricate a minimal but valid PatientRecord."""
        dom = Domain(nx=8, ny=8, dx=15.0)
        n_cells = 40
        x = rng.uniform(0, dom.extent_um[0], n_cells)
        y = rng.uniform(0, dom.extent_um[1], n_cells)
        # Half tumor, half stromal.
        ctype = np.where(np.arange(n_cells) < n_cells // 2,
                         int(CellType.TUMOR), int(CellType.STROMAL))
        state = np.full(n_cells, int(CellState.VIABLE))
        sens = np.clip(base_sens + 0.05 * rng.standard_normal(n_cells), 0.0, 1.0)
        cells = CellPopulation(x=x, y=y, cell_type=ctype, state=state, sensitivity=sens)
        fields = TissueFields(ecm_density=dom.ones() * 0.5,
                              vascular_density=dom.ones() * 0.2)
        micro = Microenvironment(domain=dom, fields=fields, cells=cells, patient_id=pid)
        feats = {
            "mgmt_methylated": float(base_sens > 0.5),
            "vascular_density": 0.2 + 0.1 * base_sens,
            "expr_gene_a": base_sens + 0.02 * rng.standard_normal(),
            "expr_gene_b": 1.0 - base_sens,
        }
        # Deliberately drop a key from one record to exercise robustness.
        if pid == "P3":
            feats.pop("expr_gene_b")
        return PatientRecord(patient_id=pid, micro=micro,
                             observed_response=resp, features=feats)

    # Binary-outcome cohort.
    records_bin = [
        _make_record("P1", 0.2, 0.0),
        _make_record("P2", 0.8, 1.0),
        _make_record("P3", 0.3, 0.0),
        _make_record("P4", 0.9, 1.0),
        _make_record("P5", 0.6, 1.0),
        _make_record("P6", 0.4, 0.0),
    ]
    # Continuous-outcome cohort (e.g. -log IC50 scaled into [0,1]).
    records_cont = [
        _make_record("C1", 0.2, 0.15),
        _make_record("C2", 0.8, 0.82),
        _make_record("C3", 0.5, 0.44),
        _make_record("C4", 0.9, 0.91),
        _make_record("C5", 0.35, 0.30),
        _make_record("C6", 0.65, 0.70),
    ]

    def _check(name: str, pred: np.ndarray, n_expected: int) -> None:
        pred = np.asarray(pred)
        ok_len = pred.shape == (n_expected,)
        ok_finite = bool(np.all(np.isfinite(pred)))
        ok_range = bool(np.all((pred >= -1e-9) & (pred <= 1.0 + 1e-9)))
        status = "PASS" if (ok_len and ok_finite and ok_range) else "FAIL"
        rng_str = (f"[{pred.min():.3f},{pred.max():.3f}]" if pred.size else "[]")
        print(f"  [{status}] {name:32s} len={pred.shape} "
              f"finite={ok_finite} in[0,1]={ok_range} range={rng_str}")
        assert ok_len and ok_finite and ok_range, f"{name} failed self-test"

    print("Self-test: fit-free baselines (binary cohort, n=6)")
    _check("baseline_mean", baseline_mean(records_bin), 6)
    _check("baseline_intrinsic_sensitivity",
           baseline_intrinsic_sensitivity(records_bin), 6)

    print("Self-test: expression_logreg (binary train/test split)")
    train_b, test_b = records_bin[:4], records_bin[4:]
    _check("expression_logreg[binary]",
           baseline_expression_logreg(train_b, test_b), len(test_b))

    print("Self-test: expression_logreg (continuous train/test split)")
    train_c, test_c = records_cont[:4], records_cont[4:]
    _check("expression_logreg[continuous]",
           baseline_expression_logreg(train_c, test_c), len(test_c))

    print("Self-test: determinism check (logreg run twice)")
    p1 = baseline_expression_logreg(train_b, test_b)
    p2 = baseline_expression_logreg(train_b, test_b)
    assert np.array_equal(p1, p2), "logreg not deterministic"
    print(f"  [PASS] identical across runs: {p1}")

    print("Self-test: edge cases")
    _check("empty-mean", baseline_mean([]), 0)
    _check("single-record-intrinsic",
           baseline_intrinsic_sensitivity([records_bin[0]]), 1)
    # Single-class train set -> should not crash.
    single_class = [records_bin[0], records_bin[2], records_bin[5]]  # all resp=0
    _check("single-class-logreg",
           baseline_expression_logreg(single_class, test_b), len(test_b))

    print("Self-test: registry")
    reg = all_baselines()
    for k, spec in reg.items():
        print(f"  {k:24s} needs_split={spec.needs_split}  -- {spec.description}")
    assert set(reg) == {"mean", "intrinsic_sensitivity", "expression_logreg"}
    assert reg["expression_logreg"].needs_split is True
    assert reg["mean"].needs_split is False

    print("\nAll baseline self-tests passed.")
