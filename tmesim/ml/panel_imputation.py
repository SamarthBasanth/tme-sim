"""ML imputation of genes/states a targeted panel does not carry — SWEPT, not single-target.

WHY THIS MODULE EXISTS
-----------------------
Three real GBM cell/spot-resolved sections carry a targeted transcript panel
missing a clinically or mechanistically load-bearing gene: the 319-gene Xenium
Human Brain panel (``data_cache/xenium_gbm``) carries NONE of the hypoxia genes
VEGFA/CA9/ADM/NDRG1/SLC2A1/HIF1A and no MGMT; the 480-gene Immuno-Oncology panel
(``data_cache/xenium_gbm_io``) carries VEGFA/HIF1A but not CA9/ADM/NDRG1/SLC2A1,
not MGMT, not TOP2A/PCNA, not most of ``real_pipeline.ECM_MARKERS``, and not
PECAM1. This module asks: given the OTHER genes a panel actually measures, can a
model predict the value of a gene it does not carry — and if so, how well,
HONESTLY validated on sections the model never saw during training.

THE ONE NON-NEGOTIABLE RULE
----------------------------
Every validation number in this module is computed on a SECTION the model was
not trained on — never a train/test split of cells or spots WITHIN one section.
A within-section split would let the model exploit patient- and batch-specific
structure it will never see again outside that patient, and would look far
better than the number that matters: does this generalise to a NEW section.
:func:`leave_one_section_out_cv` enforces this by construction (each fold trains
on two of the three real GBM-Space sections and tests on the third, entire,
untouched section) and :func:`cross_panel_check` goes further, testing on a
completely different ASSAY (Xenium single-cell, not Visium spots) that took no
part in training at all.

GROUND-TRUTH RESOURCE: GBM-Space (Visium, whole transcriptome, 3 real sections)
---------------------------------------------------------------------------
``data_cache/gbmspace_visium/{AT10-BRA-5-FO-2_1,AT10-BRA-5-FO-3_1,AT10-BRA-5-FO-4_2}``
(De Jong et al., bioRxiv 2025, via ``tmesim.data.adapters.gbmspace``) is whole
TRANSCRIPTOME (36,591 genes) Visium from three real, distinct GBM tissue blocks
— the only real, multi-SECTION resource this project has that measures every
gene this sweep cares about. It is spot-resolution (~1-10 cells pooled per
55 um spot), NOT single-cell, which is the honest limit stated everywhere below:
primary validation here is spot-resolution.

SECONDARY, HARDER CHECK: cross-panel / cross-assay generalisation
-------------------------------------------------------------------
Where a real single-cell Xenium section (the 319-gene brain panel, the 480-gene
IO panel, or the 366-gene GSE324268 section) actually measures the swept gene,
:func:`cross_panel_check` refits the SAME model family on ALL THREE GBM-Space
sections pooled, restricted to whichever feature genes that Xenium panel also
carries, and scores it against that panel's own REAL measured value — a Visium
spot -> Xenium single-cell, cross-assay transfer. This is the number closest to
"if I imputed this gene on my actual targeted panel, would it be any good", and
it is expected to be, and often is, worse than the within-Visium number: assay
chemistry and resolution both differ. Both numbers are reported; neither is
hidden because it looked worse.

MODEL: RandomForestRegressor (scikit-learn 1.6.1, already a project dependency)
---------------------------------------------------------------------------
Chosen over gradient boosting only by a coin's width (both are defensible here);
over anything from a deep-learning stack for the reason stated everywhere else in
this project: three independent tissue sections is not enough data to fit or
honestly validate a network, and PyTorch is not an installed dependency of this
project. Hyperparameters (``n_estimators=200``, ``max_depth=10``,
``min_samples_leaf=4``) are DECLARED, not tuned on any held-out score — tuning
against the numbers reported here would make them optimistic by construction.

THE CHEAP BASELINE
-------------------
Every leave-one-section-out fold also fits the best SINGLE gene (chosen by
absolute Pearson correlation on the TRAINING sections only, then evaluated,
unselected, on the held-out section) as a trivial one-feature linear
comparison. A forest that cannot beat its own best single input gene is not
buying anything with the other ~250-450 genes it was given.

OPT-IN ONLY
-----------
Nothing in this module is called by any existing reconstruction path.
:func:`impute_gene_score` must be called explicitly; :mod:`reconstruct.real_pipeline`
only reaches it through the explicit ``state_mode='ml_imputed'`` route (hypoxia
only — the one target this module wires all the way into the CellState schema).
Every other swept target is reachable only through :func:`impute_gene_score`
directly, tagged with provenance ``kind='ml_imputed'`` — a tier that is neither
``measured`` nor ``computed`` nor ``assumed`` (see ``app/web/js/ui/provenance.js``
for the existing taxonomy this deliberately does NOT reuse) — and it REFUSES
(raises ``ValueError``) rather than silently degrading when the calling panel
does not carry enough of the required feature genes.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import sparse as _sp

__all__ = [
    "TargetSpec",
    "TARGET_REGISTRY",
    "TARGETS_BY_NAME",
    "panel_genes",
    "composite_score",
    "leave_one_section_out_cv",
    "cross_panel_check",
    "run_sweep",
    "impute_gene_score",
    "MIN_FEATURE_FLOOR",
    "GBMSPACE_SECTIONS",
    "XENIUM_PANELS",
]

# ---------------------------------------------------------------------------
# Real data locations (all pre-existing caches; nothing is downloaded here).
# ---------------------------------------------------------------------------
GBMSPACE_SECTIONS: Tuple[str, ...] = (
    "AT10-BRA-5-FO-2_1", "AT10-BRA-5-FO-3_1", "AT10-BRA-5-FO-4_2",
)
_GBMSPACE_ROOT = os.path.join("data_cache", "gbmspace_visium")

#: The three real Xenium single-cell sections this project has cached, keyed by
#: the panel name used throughout ``reconstruct.marker_panels``.
XENIUM_PANELS: Dict[str, str] = {
    "brain_319": os.path.join("data_cache", "xenium_gbm"),
    "io_480": os.path.join("data_cache", "xenium_gbm_io"),
    "gse_366": os.path.join("data_cache", "gse324268_xenium"),
}

#: A cross-panel check below this many overlapping feature genes is refused
#: (skipped, with a reason recorded) rather than run on too little to mean
#: anything.
MIN_FEATURE_FLOOR = 20

#: Below this many present feature genes on a CALLING panel, ``impute_gene_score``
#: refuses outright (raises) rather than predicting off a near-empty input.
_IMPUTE_MIN_FEATURES = 20


# ---------------------------------------------------------------------------
# Data loading (cached; all reads go through the existing, unmodified loaders).
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _load_xenium(panel_key: str) -> Dict[str, object]:
    from ..data.real.singlecell_spatial import load
    return load(path=XENIUM_PANELS[panel_key], dense=False)


@lru_cache(maxsize=None)
def _load_gbmspace(section_id: str) -> Dict[str, object]:
    from ..data.real.visium import load_from_dirs
    base = os.path.join(_GBMSPACE_ROOT, section_id)
    return load_from_dirs(
        os.path.join(base, "spatial"), os.path.join(base, "filtered_feature_bc_matrix"),
        sample_id=section_id, marker_report=False,
    )


def panel_genes(panel_key: str) -> List[str]:
    """The full gene list of a real cached Xenium panel (319 / 480 / 366)."""
    return list(np.asarray(_load_xenium(panel_key)["genes"], dtype=str))


def clear_cache() -> None:
    """Drop cached loads (tests use this so different cache states are isolated)."""
    _load_xenium.cache_clear()
    _load_gbmspace.cache_clear()
    _fit_pooled_gbmspace_model.cache_clear()


# ---------------------------------------------------------------------------
# Feature / target extraction — small, local, self-contained (does not import
# private helpers from reconstruct.real_xenium / cell_state, to keep this
# package's only coupling to the rest of tmesim at the PUBLIC loader level).
# ---------------------------------------------------------------------------
def _extract_columns(
    expression, genes: Sequence[str], wanted: Sequence[str]
) -> Tuple[List[str], np.ndarray]:
    """Dense ``(n_cells, k)`` float64 sub-matrix of ``wanted`` genes actually on
    the panel, plus the list of which of ``wanted`` were present (panel order)."""
    genes_arr = np.asarray(genes, dtype=str)
    idx = {g: j for j, g in enumerate(genes_arr)}
    present = [g for g in wanted if g in idx]
    if not present:
        return present, np.zeros((int(expression.shape[0]), 0), dtype=np.float64)
    cols = [idx[g] for g in present]
    if _sp.issparse(expression):
        sub = np.asarray(expression[:, cols].todense(), dtype=np.float64)
    else:
        sub = np.asarray(expression, dtype=np.float64)[:, cols]
    return present, sub


def _log1p_zscore(dense: np.ndarray) -> np.ndarray:
    """log1p then per-column z-score across the given population. Local copy of
    the same convention ``reconstruct.real_xenium._zscore_columns`` uses,
    kept self-contained rather than importing a private helper."""
    x = np.log1p(np.asarray(dense, dtype=np.float64))
    if x.shape[1] == 0:
        return x
    mean = x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    std[std == 0.0] = 1.0
    return (x - mean) / std


def composite_score(
    expression, genes: Sequence[str], target_genes: Sequence[str]
) -> Tuple[Optional[np.ndarray], List[str], Dict[str, float]]:
    """Mean log1p z-score over whichever of ``target_genes`` the panel measures.

    Returns ``(score, present_genes, detection_rates)``. ``score`` is ``None``
    (not assigned to any cell) if NONE of ``target_genes`` is on the panel — the
    same refusal cell_state.py's ``marker_detection`` applies, generalised from a
    binary positivity call to a continuous composite. When only a SUBSET of
    ``target_genes`` is present, the composite is honestly computed over that
    subset alone and ``present_genes`` says so; the caller must read it.
    """
    present, sub = _extract_columns(expression, genes, target_genes)
    if not present:
        return None, present, {}
    z = _log1p_zscore(sub)
    score = z.mean(axis=1)
    detection = {g: float((sub[:, i] > 0).mean()) for i, g in enumerate(present)}
    return score, present, detection


# ---------------------------------------------------------------------------
# The model.
# ---------------------------------------------------------------------------
def _default_regressor(random_state: int = 0):
    from sklearn.ensemble import RandomForestRegressor
    return RandomForestRegressor(
        n_estimators=200, max_depth=10, min_samples_leaf=4,
        n_jobs=-1, random_state=random_state,
    )


def _r2_pearson(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[float, float]:
    from scipy.stats import pearsonr
    from sklearn.metrics import r2_score
    r2 = float(r2_score(y_true, y_pred))
    if np.std(y_true) > 0 and np.std(y_pred) > 0:
        r, _ = pearsonr(y_true, y_pred)
        r = float(r)
    else:
        r = float("nan")
    return r2, r


def _best_single_feature_baseline(
    Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, yte: np.ndarray,
    feat_names: Sequence[str],
) -> Tuple[float, Optional[str]]:
    """The cheap baseline: pick the single feature most correlated with the
    target ON THE TRAINING SECTIONS ONLY, fit a 1-D OLS line to it, and score
    that (unselected-on-test) line on the held-out section. A forest that does
    not beat this is not using its other features for anything."""
    from sklearn.metrics import r2_score
    if Xtr.shape[1] == 0:
        return float("nan"), None
    corrs = np.zeros(Xtr.shape[1])
    for j in range(Xtr.shape[1]):
        xj = Xtr[:, j]
        if xj.std() == 0 or ytr.std() == 0:
            continue
        c = np.corrcoef(xj, ytr)[0, 1]
        corrs[j] = 0.0 if np.isnan(c) else c
    j_best = int(np.argmax(np.abs(corrs)))
    xj_tr, xj_te = Xtr[:, j_best], Xte[:, j_best]
    if xj_tr.std() == 0:
        return float("nan"), feat_names[j_best]
    b, a = np.polyfit(xj_tr, ytr, 1)
    pred = a + b * xj_te
    return float(r2_score(yte, pred)), feat_names[j_best]


# ---------------------------------------------------------------------------
# Target registry.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TargetSpec:
    """One swept prediction target.

    ``gap_panels``
        Real cached panels that do NOT measure this target — the genuine
        product gap. The FIRST entry supplies the feature-gene set (that
        panel's full gene list, minus the target genes) used for the primary
        GBM-Space leave-one-section-out CV, so the CV answers "could THIS
        panel's other genes predict what it is missing".
    ``ground_truth_panels``
        Real cached panels that DO measure this target — used by
        :func:`cross_panel_check` as an independent, single-cell, cross-assay
        test of the model trained on GBM-Space.
    ``is_real_gap``
        False for the handful of targets (e.g. MKI67) already measured on
        every real cached panel — these are included as a PIPELINE SANITY
        CHECK / positive control, not because any current panel needs them
        imputed.
    """
    name: str
    category: str
    target_genes: Tuple[str, ...]
    gap_panels: Tuple[str, ...]
    ground_truth_panels: Tuple[str, ...] = field(default_factory=tuple)
    is_real_gap: bool = True
    notes: str = ""


TARGET_REGISTRY: Tuple[TargetSpec, ...] = (
    # --- hypoxia -------------------------------------------------------
    TargetSpec(
        "hypoxia_composite", "hypoxia",
        ("VEGFA", "CA9", "ADM", "NDRG1", "SLC2A1"),  # == cell_state.HYPOXIA_MARKERS
        gap_panels=("brain_319",), ground_truth_panels=("gse_366", "io_480"),
        notes="io_480 ground truth is computed from VEGFA alone (the only one "
              "of the 5 it carries) — a 1-gene composite, reported as such.",
    ),
    TargetSpec("hypoxia_VEGFA", "hypoxia", ("VEGFA",),
               gap_panels=("brain_319",), ground_truth_panels=("io_480", "gse_366")),
    TargetSpec("hypoxia_CA9", "hypoxia", ("CA9",),
               gap_panels=("brain_319", "io_480"), ground_truth_panels=("gse_366",)),
    TargetSpec("hypoxia_HIF1A", "hypoxia", ("HIF1A",),
               gap_panels=("brain_319",), ground_truth_panels=("io_480", "gse_366")),
    TargetSpec("hypoxia_ADM", "hypoxia", ("ADM",),
               gap_panels=("brain_319", "io_480"), ground_truth_panels=("gse_366",)),
    TargetSpec("hypoxia_NDRG1", "hypoxia", ("NDRG1",),
               gap_panels=("brain_319", "io_480"), ground_truth_panels=(),
               notes="Not measured by ANY cached single-cell Xenium panel — "
                     "GBM-Space leave-one-section-out is the only check available."),
    TargetSpec("hypoxia_SLC2A1", "hypoxia", ("SLC2A1",),
               gap_panels=("brain_319", "io_480"), ground_truth_panels=("gse_366",)),

    # --- resistance ------------------------------------------------------
    TargetSpec(
        "mgmt", "resistance", ("MGMT",),
        gap_panels=("brain_319", "io_480"), ground_truth_panels=("gse_366",),
        notes="TMZ resistance determinant; absent from BOTH incumbent flagship "
              "panels. Distinct from molecular.mgmt, which takes a per-patient "
              "MGMT SCALAR as an explicit user input rather than predicting one "
              "from spatial co-expression.",
    ),

    # --- proliferation -----------------------------------------------------
    TargetSpec("proliferation_composite", "proliferation", ("MKI67", "TOP2A"),
               gap_panels=("io_480",), ground_truth_panels=("brain_319", "gse_366")),
    TargetSpec("proliferation_MKI67", "proliferation", ("MKI67",),
               gap_panels=("brain_319",), ground_truth_panels=("io_480", "gse_366"),
               is_real_gap=False,
               notes="Measured on every cached panel already — sanity/positive "
                     "control, not a real gap."),
    TargetSpec("proliferation_TOP2A", "proliferation", ("TOP2A",),
               gap_panels=("io_480",), ground_truth_panels=("brain_319", "gse_366")),

    # --- ECM / matrix --------------------------------------------------------
    TargetSpec("ecm_composite", "ecm",
               ("COL1A2", "COL12A1", "COL25A1", "POSTN", "VCAN", "DCN", "FBLN1",
                "LAMA2", "THBS1", "SPON1", "TGFBI"),  # == real_pipeline.ECM_MARKERS
               gap_panels=("io_480",), ground_truth_panels=("brain_319", "gse_366")),
    TargetSpec("ecm_COL1A2", "ecm", ("COL1A2",),
               gap_panels=("io_480", "gse_366"), ground_truth_panels=("brain_319",)),
    TargetSpec("ecm_POSTN", "ecm", ("POSTN",),
               gap_panels=("io_480",), ground_truth_panels=("brain_319", "gse_366")),

    # --- efflux transporters ------------------------------------------------
    TargetSpec("efflux_ABCB1", "efflux", ("ABCB1",),
               gap_panels=("brain_319", "io_480", "gse_366"), ground_truth_panels=(),
               notes="On no cached targeted panel at all — GBM-Space internal "
                     "CV is the only check this data snapshot supports."),
    TargetSpec("efflux_ABCG2", "efflux", ("ABCG2",),
               gap_panels=("brain_319", "io_480", "gse_366"), ground_truth_panels=()),

    # --- vascular / mural ----------------------------------------------------
    TargetSpec("vascular_PECAM1", "vascular", ("PECAM1",),
               gap_panels=("io_480",), ground_truth_panels=("brain_319", "gse_366")),
    TargetSpec("vascular_composite", "vascular",
               ("CD34", "CLDN5", "VWF", "PLVAP", "KDR", "FLT1", "ACKR1", "DLL4"),
               gap_panels=("brain_319",), ground_truth_panels=("gse_366", "io_480"),
               notes="io_480 carries only VWF/PLVAP/FLT1/ACKR1 of the 8 — a "
                     "4-gene partial composite there, reported as such."),
    TargetSpec("mural_composite", "vascular", ("RGS5", "PDGFRB", "NOTCH3", "ACTA2"),
               gap_panels=("brain_319",), ground_truth_panels=("io_480", "gse_366"),
               notes="gse_366 carries RGS5/PDGFRB/ACTA2 but not NOTCH3 — a "
                     "3-gene partial composite there, reported as such."),
)

TARGETS_BY_NAME: Dict[str, TargetSpec] = {t.name: t for t in TARGET_REGISTRY}


# ---------------------------------------------------------------------------
# Primary validation: leave-one-GBM-Space-section-out.
# ---------------------------------------------------------------------------
def _gbmspace_feature_target(
    section_id: str, feature_genes: Sequence[str], target_genes: Sequence[str]
):
    d = _load_gbmspace(section_id)
    genes, expr = d["genes"], d["expression"]
    present_feat, Xsub = _extract_columns(expr, genes, feature_genes)
    Xz = _log1p_zscore(Xsub)
    y, present_target, det = composite_score(expr, genes, target_genes)
    return present_feat, Xz, y, present_target, det, int(expr.shape[0])


def leave_one_section_out_cv(
    target: TargetSpec,
    feature_genes: Sequence[str],
    model_factory=_default_regressor,
    random_state: int = 0,
) -> Dict[str, object]:
    """Honest primary validation: train on 2 of the 3 real GBM-Space sections,
    test on the ENTIRE third, for all three choices of held-out section.

    Never splits cells/spots within one section — the non-negotiable rule this
    module exists to keep. Returns per-fold and pooled-out-of-fold R^2 /
    Pearson r, plus the cheap single-gene baseline for each fold.
    """
    per_section = {
        s: _gbmspace_feature_target(s, feature_genes, target.target_genes)
        for s in GBMSPACE_SECTIONS
    }
    feat_names = per_section[GBMSPACE_SECTIONS[0]][0]
    if any(per_section[s][0] != feat_names for s in GBMSPACE_SECTIONS[1:]):
        # Whole-transcriptome sections on the same platform; should never
        # differ, but fail loudly rather than silently misalign columns.
        raise ValueError(
            f"{target.name}: GBM-Space sections disagree on which feature genes "
            "are present — refusing to run a misaligned CV"
        )
    if not feat_names:
        return {
            "target": target.name, "skipped": True,
            "reason": "no feature genes present in GBM-Space", "n_features": 0,
        }
    target_present_any = per_section[GBMSPACE_SECTIONS[0]][3]
    if per_section[GBMSPACE_SECTIONS[0]][2] is None:
        return {
            "target": target.name, "skipped": True,
            "reason": "target gene(s) not present in GBM-Space", "n_features": len(feat_names),
        }

    folds = []
    oof_true, oof_pred = [], []
    for held in GBMSPACE_SECTIONS:
        train_secs = [s for s in GBMSPACE_SECTIONS if s != held]
        Xtr = np.concatenate([per_section[s][1] for s in train_secs], axis=0)
        ytr = np.concatenate([per_section[s][2] for s in train_secs], axis=0)
        Xte, yte = per_section[held][1], per_section[held][2]

        model = model_factory(random_state=random_state)
        model.fit(Xtr, ytr)
        pred = model.predict(Xte)
        r2, r = _r2_pearson(yte, pred)
        base_r2, base_gene = _best_single_feature_baseline(Xtr, ytr, Xte, yte, feat_names)

        folds.append({
            "held_out_section": held,
            "n_train": int(Xtr.shape[0]), "n_test": int(Xte.shape[0]),
            "r2": r2, "pearson_r": r,
            "baseline_r2_single_gene": base_r2, "baseline_gene": base_gene,
            "beats_baseline": bool(np.isfinite(r2) and np.isfinite(base_r2) and r2 > base_r2),
        })
        oof_true.append(yte)
        oof_pred.append(pred)

    oof_true_arr = np.concatenate(oof_true)
    oof_pred_arr = np.concatenate(oof_pred)
    pooled_r2, pooled_r = _r2_pearson(oof_true_arr, oof_pred_arr)

    return {
        "target": target.name,
        "skipped": False,
        "n_features": len(feat_names),
        "feature_genes": feat_names,
        "target_genes_present": target_present_any,
        "folds": folds,
        "pooled_r2": pooled_r2,
        "pooled_pearson_r": pooled_r,
        "n_pooled": int(oof_true_arr.size),
    }


# ---------------------------------------------------------------------------
# Secondary validation: cross-panel / cross-assay generalisation.
# ---------------------------------------------------------------------------
def cross_panel_check(
    target: TargetSpec,
    feature_genes: Sequence[str],
    ground_truth_panel: str,
    model_factory=_default_regressor,
    random_state: int = 0,
) -> Dict[str, object]:
    """Train on ALL THREE GBM-Space sections pooled (Visium, spot-resolution),
    test on a REAL single-cell Xenium panel section that took no part in
    training at all — a harder, cross-assay, cross-resolution check of whether
    the imputation would actually help the panel it is aimed at.
    """
    gt = _load_xenium(ground_truth_panel)
    gt_genes_set = set(np.asarray(gt["genes"], dtype=str).tolist())
    feat = [g for g in feature_genes if g in gt_genes_set]
    if len(feat) < MIN_FEATURE_FLOOR:
        return {
            "target": target.name, "ground_truth_panel": ground_truth_panel,
            "skipped": True,
            "reason": (f"only {len(feat)} feature genes overlap with "
                       f"{ground_truth_panel} (floor {MIN_FEATURE_FLOOR})"),
        }

    Xtr_list, ytr_list = [], []
    for s in GBMSPACE_SECTIONS:
        present_feat, Xz, y, _present_t, _det, _n = _gbmspace_feature_target(
            s, feat, target.target_genes
        )
        if present_feat != feat:
            # GBM-Space is whole-transcriptome; this would mean a requested
            # feature gene is missing there, which every check in this module
            # has verified is not the case for the registry above. Fail loud
            # rather than silently drop columns and misalign across sections.
            raise ValueError(
                f"{target.name}: GBM-Space section {s!r} is missing feature "
                f"genes {sorted(set(feat) - set(present_feat))}"
            )
        Xtr_list.append(Xz)
        ytr_list.append(y)
    Xtr = np.concatenate(Xtr_list, axis=0)
    ytr = np.concatenate(ytr_list, axis=0)
    if ytr is None or np.std(ytr) == 0:
        return {
            "target": target.name, "ground_truth_panel": ground_truth_panel,
            "skipped": True, "reason": "degenerate training target (zero variance)",
        }

    model = model_factory(random_state=random_state)
    model.fit(Xtr, ytr)

    _present_feat_gt, Xte_sub = _extract_columns(gt["expression"], gt["genes"], feat)
    Xte = _log1p_zscore(Xte_sub)
    yte, present_target_gt, det_gt = composite_score(gt["expression"], gt["genes"], target.target_genes)
    if yte is None:
        return {
            "target": target.name, "ground_truth_panel": ground_truth_panel,
            "skipped": True, "reason": "target genes not present in ground-truth panel",
        }

    pred = model.predict(Xte)
    r2, r = _r2_pearson(yte, pred)
    return {
        "target": target.name,
        "ground_truth_panel": ground_truth_panel,
        "skipped": False,
        "n_features_used": len(feat),
        "n_test_cells": int(yte.size),
        "target_genes_present_in_ground_truth": present_target_gt,
        "target_gene_detection_rates": det_gt,
        "r2": r2,
        "pearson_r": r,
    }


# ---------------------------------------------------------------------------
# The sweep.
# ---------------------------------------------------------------------------
def run_sweep(
    targets: Optional[Sequence[TargetSpec]] = None,
    random_state: int = 0,
) -> Dict[str, Dict[str, object]]:
    """Run every registered target's primary CV and (where available) its
    cross-panel checks. Returns ``{target_name: {...}}``; nothing is averaged
    across targets — each is its own honest number, per the design brief."""
    targets = list(targets) if targets is not None else list(TARGET_REGISTRY)
    out: Dict[str, Dict[str, object]] = {}
    for spec in targets:
        gap_panel = spec.gap_panels[0]
        feature_genes = [g for g in panel_genes(gap_panel) if g not in spec.target_genes]
        loso = leave_one_section_out_cv(spec, feature_genes, random_state=random_state)
        cross = {
            gt: cross_panel_check(spec, feature_genes, gt, random_state=random_state)
            for gt in spec.ground_truth_panels
        }
        out[spec.name] = {
            "spec": spec,
            "gap_panel_used_for_features": gap_panel,
            "n_feature_genes_requested": len(feature_genes),
            "loso_cv": loso,
            "cross_panel": cross,
        }
    return out


# ---------------------------------------------------------------------------
# Opt-in inference API.
# ---------------------------------------------------------------------------
@lru_cache(maxsize=None)
def _fit_pooled_gbmspace_model(target_name: str, feature_genes: Tuple[str, ...],
                                random_state: int = 0):
    """Model trained on ALL THREE GBM-Space sections, for the given exact
    feature-gene tuple. Cached so repeated ``impute_gene_score`` calls for the
    same target/panel do not refit. NOT itself a validation function — the
    honest score for this exact configuration comes from :func:`run_sweep`."""
    spec = TARGETS_BY_NAME[target_name]
    Xtr_list, ytr_list = [], []
    for s in GBMSPACE_SECTIONS:
        present_feat, Xz, y, _pt, _det, _n = _gbmspace_feature_target(
            s, feature_genes, spec.target_genes
        )
        if present_feat != list(feature_genes):
            raise ValueError(
                f"{target_name}: GBM-Space is missing requested feature genes "
                f"{sorted(set(feature_genes) - set(present_feat))}"
            )
        Xtr_list.append(Xz)
        ytr_list.append(y)
    Xtr = np.concatenate(Xtr_list, axis=0)
    ytr = np.concatenate(ytr_list, axis=0)
    model = _default_regressor(random_state=random_state)
    model.fit(Xtr, ytr)
    return model


def impute_gene_score(
    expression,
    genes: Sequence[str],
    target: str,
    gap_panel: str = "brain_319",
    validation: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    """Predict ``target``'s composite score from THIS panel's other genes. OPT-IN.

    Never called automatically by any reconstruction path. Raises rather than
    degrading silently if the calling panel does not carry enough of the
    required feature genes (:data:`_IMPUTE_MIN_FEATURES` floor).

    Parameters
    ----------
    expression, genes : the CALLING panel's own real counts and gene symbols —
        NOT GBM-Space; this is the panel asking for an estimate of a gene it
        does not measure.
    target : one of :data:`TARGET_REGISTRY`'s names (e.g. ``'hypoxia_composite'``,
        ``'mgmt'``).
    gap_panel : which registered real panel's gene list to use as the feature
        template (default ``'brain_319'``, the panel with the fewest measured
        hypoxia/MGMT/ECM genes of the three). The features actually used are
        the intersection of that template with what ``genes`` truly carries.
    validation : the honest score to attach to ``provenance['validation']``
        (typically ``run_sweep()[target]``, or a cached summary of it — see
        ``results/panel_gap_imputation_sweep.md``). If omitted, the returned
        provenance says explicitly that no validation number was attached
        rather than inventing one.

    Returns
    -------
    dict with ``score`` (n_cells,), ``target``, ``target_genes``,
    ``n_features_used``, ``present_feature_genes``, and ``provenance`` — a
    block carrying ``kind='ml_imputed'`` (a tier distinct from measured /
    literature / computed / assumed), what the model is, what it was trained
    on, and (if supplied) the honest validation numbers.

    Raises
    ------
    ValueError
        If ``target`` is unknown, or if the calling panel carries fewer than
        :data:`_IMPUTE_MIN_FEATURES` of the required feature genes.
    """
    if target not in TARGETS_BY_NAME:
        raise ValueError(f"unknown target {target!r}; known: {sorted(TARGETS_BY_NAME)}")
    spec = TARGETS_BY_NAME[target]
    feature_genes = tuple(g for g in panel_genes(gap_panel) if g not in spec.target_genes)

    present_feat, Xsub = _extract_columns(expression, genes, feature_genes)
    if len(present_feat) < _IMPUTE_MIN_FEATURES:
        raise ValueError(
            f"impute_gene_score refuses: this panel carries only "
            f"{len(present_feat)} of the {len(feature_genes)} feature genes "
            f"{target!r} needs (floor {_IMPUTE_MIN_FEATURES}). This is not a "
            "supported use of ml_imputed on this panel."
        )
    Xz = _log1p_zscore(Xsub)

    model = _fit_pooled_gbmspace_model(target, tuple(present_feat))
    pred = model.predict(Xz)

    val_note = (
        "NO VALIDATION NUMBER WAS ATTACHED to this call — pass validation= "
        "(e.g. from run_sweep()) so a caller can see how much to trust it."
        if validation is None else None
    )

    return {
        "score": pred,
        "target": target,
        "target_genes": spec.target_genes,
        "n_features_used": len(present_feat),
        "present_feature_genes": present_feat,
        "provenance": {
            "kind": "ml_imputed",
            "model": (
                "RandomForestRegressor (scikit-learn), trained on GBM-Space "
                "Visium (3 real sections, whole transcriptome, spot-resolution)"
            ),
            "trained_on": "GBM-Space Visium — a DIFFERENT assay and resolution "
                           "than the panel this score is being applied to",
            "validation": validation,
            "validation_note": val_note,
            "caveat": (
                f"MODEL-IMPUTED, NOT MEASURED. This is a machine-learning "
                f"ESTIMATE of what {list(spec.target_genes)} would read on this "
                "specimen, inferred from correlated genes this panel DOES "
                "measure. It is not a transcript count. Read the attached "
                "validation number (if any) before trusting its magnitude — "
                "see the module docstring for why cross-assay generalisation "
                "is the harder, more relevant number."
            ),
        },
    }
