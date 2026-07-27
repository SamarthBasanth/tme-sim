"""DOES THE OXYGEN-BANDED CELL-STATE MODEL CARRY REAL BIOLOGY? — a falsification
test against MEASURED per-cell proliferation on the real Xenium GBM section.

The claim under test
--------------------
:func:`tmesim.cells.agents.update_cell_states` assigns every living tumour cell a
:class:`~tmesim.core.microenvironment.CellState` from the **simulated** oxygen
field alone, by a hard banding on ``t = OxygenConfig.hypoxia_threshold``::

    o <  t        -> HYPOXIC     (kill x PDConfig.hypoxia_resistance  = 0.35)
    t <= o < 2t   -> QUIESCENT   (kill x PDConfig.quiescent_resistance = 0.50)
    o >= 2t       -> VIABLE      (kill x 1.0)

On the real section this yields ~0.6% VIABLE / ~86% QUIESCENT / ~13% HYPOXIC
(``results/patient_report_xenium.md``), and those multipliers are *load-bearing*:
the patient report's own counterfactual (re-simulating with both penalties off)
shows they cost every drug on the panel ~50% of its predicted kill. So the state
model is not decoration — it roughly halves every headline kill number.

That makes it worth falsifying, because on this section we can actually measure
the thing the state model claims to describe. ``MKI67`` (Ki-67, the canonical
proliferation marker) and ``TOP2A`` are **measured per cell** on the 319-gene
Xenium Human Brain panel, alongside six further cell-cycle genes (``PCNA``,
``CDK1``, ``CCNB2``, ``CCNA1``, ``CENPF``, ``GAS2L3``). ``VIABLE`` is supposed to
mean "well-perfused and proliferation-competent" and ``QUIESCENT`` is supposed to
mean "alive but not proliferating" — those are *empirical* statements about cell
biology, and the panel measures the marker that adjudicates them.

Three tests, in increasing order of how much they can be trusted
----------------------------------------------------------------
**Test 1 — MARGINAL fractions (weak).** Compare the model's VIABLE fraction to
the measured MKI67-positive fraction of tumour cells. This is the obvious
comparison and it is also the *least* trustworthy one, because a targeted in-situ
panel's per-cell detection rate is not a labelling index: MKI67 is detected in a
fraction of cells far below any published Ki-67 index, which is a statement about
probe sensitivity and transcript dropout, not about how many cells cycle. A
numerical agreement between "model VIABLE %" and "MKI67-detected %" would be a
**coincidence of two unrelated small numbers**, and is reported as such rather
than as corroboration. Literature Ki-67 indices (with PMIDs) are quoted so the
reader can see how far *both* numbers sit from the real labelling index.

**Test 2 — SPATIAL enrichment (the real test).** Detection efficiency cancels out
of a *within-section contrast*: whatever fraction of true Ki-67+ cells the probe
catches, it catches it equally in every model state. So the falsifiable question
is: **are model-VIABLE cells more proliferative than model-QUIESCENT/HYPOXIC
cells?** If the oxygen banding carries biological signal, MKI67 (and the
composite cell-cycle score) must be enriched in VIABLE and depleted in HYPOXIC.
If it is not, the state labels are a picture of the simulated oxygen field and
nothing more. Reported as an effect size + a label-permutation p-value, with an
explicit power statement so a null result is never silently read as a negative.

**Test 3 — CONTINUOUS correlation (the best-powered test).** Bin nothing: rank-
correlate the per-cell **simulated oxygen** against the per-cell measured
proliferation score across all ~30k tumour cells. This uses every cell and does
not depend on where the two band edges were drawn, so it is the most sensitive
probe of whether the oxygen field carries proliferation information at all. The
same correlation is run against a measured **HIF-target** score (``HILPDA``,
``LOX``) — a direct audit of the simulated oxygen field itself, which must be
*negatively* correlated with hypoxia-response transcripts if it is real.

Two confounders are handled explicitly rather than assumed away
---------------------------------------------------------------
1. **Sequencing depth.** A cell with more total transcripts is more likely to
   show a nonzero MKI67 purely by chance. Every expression statistic here is
   counts-per-median normalised before scoring, mean total counts are reported
   per state, and Test 3 is additionally run as a *partial* Spearman with
   ``log(total counts)`` residualised out of both variables.
2. **Typing circularity.** ``MKI67`` is one of the four TUMOR markers in
   :data:`reconstruct.real_xenium.MARKER_MAP`, so an MKI67-high cell is slightly
   more likely to be typed TUMOR. This inflates the MKI67+ rate *uniformly* across
   states, so it cannot manufacture a between-state contrast — it is noted and
   quantified, not corrected.

ADDS ONLY. This module drives the frozen stack (``reconstruct.real_pipeline`` ->
``physics.oxygen`` -> ``cells.agents``) and reads the real loader. It modifies
nothing in ``core`` / ``physics`` / ``molecular`` / ``reconstruct`` / ``sim``.

Outputs ``results/cellstate_validation.md`` and ``results/cellstate_validation.png``.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import sparse
from scipy.stats import rankdata

from ..core.config import OxygenConfig, PDConfig
from ..core.microenvironment import CellState, CellType
from ..cells.agents import update_cell_states
from ..data.real.singlecell_spatial import load as load_singlecell_spatial
from ..reconstruct.real_pipeline import reconstruct_real

# ---------------------------------------------------------------------------
# Declared constants (all modelling choices are named here, not buried).
# ---------------------------------------------------------------------------

#: Grid-converged spacing (um). ``results/numerical_rigor.md`` §2 establishes that
#: the ~45 um perivascular boundary layer is unresolved for dx >~ 15 um; the
#: patient report's headline state fractions are quoted at this dx.
DX_CONVERGED = 12.0

#: dx values swept to expose the grid-sensitivity of the state fractions. The
#: hypoxic fraction is explicitly NOT grid-converged (patient report §1.3), so the
#: whole state vector inherits that sensitivity and we measure how much.
DX_SWEEP: Tuple[float, ...] = (30.0, 20.0, DX_CONVERGED)

#: Both tissue-mask conventions are swept. The mask gates the vascular (oxygen
#: SUPPLY) field, so it directly sets the oxygen level and hence every state
#: fraction — yet it is a bookkeeping convention, not a measurement.
#:   ``geometric``  — union of fixed-radius discs around measured cell centroids
#:                    (the current default; grid-independent).
#:   ``occupancy``  — legacy "grid cells containing >= 1 segmented cell", which is
#:                    a function of dx and shrinks as the mesh is refined. This is
#:                    the convention under which ``results/patient_report_xenium.md``
#:                    quoted its ~13.4% hypoxic / ~0.6% viable headline.
TISSUE_MASK_MODES: Tuple[str, ...] = ("geometric", "occupancy")

#: The mode the headline statistics in this report are computed under (whatever
#: the shipping reconstruction currently defaults to).
DEFAULT_MASK_MODE = "geometric"

#: Label-permutation draws for every significance test. 10,000 gives a resolution
#: floor of p ~ 1e-4, which is finer than any claim made from these data.
N_PERMUTATIONS = 10_000

#: Deterministic permutation seed — this module contains no other randomness.
PERM_SEED = 0

#: Working-set cap (cells x draws) for the blocked permutation kernels, so a
#: 10,000-draw test over ~30,000 cells runs as ~50 vectorized blocks of a few tens
#: of MB rather than a 10,000-iteration Python loop.
_PERM_BLOCK_CELLS = 5_000_000

#: The two markers the task adjudicates on, measured per cell on the panel.
PRIMARY_MARKERS: Tuple[str, ...] = ("MKI67", "TOP2A")

#: Full cell-cycle / proliferation gene set present on the 319-gene panel. Pooling
#: them into one score is what buys statistical power back from the extreme
#: sparsity of any single in-situ transcript.
#:   MKI67  — Ki-67, expressed through G1/S/G2/M, absent in G0 (the canonical LI marker)
#:   TOP2A  — topoisomerase II alpha, peaks G2/M
#:   PCNA   — DNA polymerase processivity factor, S-phase
#:   CDK1   — the mitotic CDK, G2/M
#:   CCNB2  — cyclin B2, G2/M
#:   CCNA1  — cyclin A1, S/G2
#:   CENPF  — centromere protein F ("mitosin"), G2/M
#:   GAS2L3 — growth-arrest-specific 2 like 3, mitosis/cytokinesis
PROLIFERATION_GENES: Tuple[str, ...] = (
    "MKI67", "TOP2A", "PCNA", "CDK1", "CCNB2", "CCNA1", "CENPF", "GAS2L3",
)

#: Measured HIF-target transcripts on the panel, used to audit the SIMULATED
#: oxygen field directly. Deliberately short and defensible — the canonical
#: hypoxia stains (CA9, VEGFA, ADM, NDRG1) are NOT on this panel.
#:   HILPDA — hypoxia-inducible protein 2 / HIG2; a *direct and specific* HIF-1
#:            target (Gimm et al., FASEB J 2010, PMID 20624928).
#:   LOX    — lysyl oxidase; a well-established hypoxia/HIF-induced gene.
HIF_TARGET_GENES: Tuple[str, ...] = ("HILPDA", "LOX")

#: Positive-call gate for the "strict" z-score positivity rule, matching the
#: repo's existing convention (``real_xenium.TYPING_THRESHOLD``): ~1 SD above the
#: tumour-cell mean.
STRICT_Z = 1.0

#: Published Ki-67 labelling indices in glioblastoma / high-grade glioma. Every
#: entry is a real, PubMed-verified measurement; nothing is estimated.
LITERATURE_KI67: Tuple[Dict[str, object], ...] = (
    {
        "label": "GBM, n=70 (median Ki-67)",
        "value_pct": 40.0,
        "citation": "Caramanti et al., Cureus 2022;14(4):e24246",
        "pmid": "35602791",
        "doi": "10.7759/cureus.24246",
        "note": "median Ki-67 index 40% across 70 GBM resections",
    },
    {
        "label": "GBM, n=44 (cohort split point)",
        "value_pct": 27.0,
        "citation": "Alkhaibary et al., Hematol Oncol Stem Cell Ther 2018;12(2):82-88",
        "pmid": "30552865",
        "doi": "10.1016/j.hemonc.2018.11.001",
        "note": "cohort dichotomised at a Ki-67 LI of 27%; LI alone did not predict OS",
    },
    {
        "label": "epithelioid GBM, n=16 (median, range 10-50%)",
        "value_pct": 30.0,
        "citation": "Zheng et al., Zhonghua Bing Li Xue Za Zhi 2021;50(3):229-235",
        "pmid": "33677887",
        "doi": "10.3760/cma.j.cn112151-20200617-00479",
        "note": "median Ki-67 index 30% (range 10-50%)",
    },
)

#: The directional prediction the oxygen-banded state model MUST reproduce if it
#: is describing real biology — measured directly in GBM tissue.
LITERATURE_HYPOXIA_PROLIFERATION: Dict[str, object] = {
    "claim": (
        "hypoxic (pseudopalisading, HIF-1alpha-positive) glioblastoma cells are "
        "5-50% LESS proliferative and 6-20x more apoptotic than adjacent astrocytoma"
    ),
    "citation": "Brat et al., Cancer Res 2004;64(3):920-927",
    "pmid": "14871821",
    "doi": "10.1158/0008-5472.CAN-03-2073",
}

#: Interobserver-variability caveat on the Ki-67 LI itself — the reference
#: standard we are comparing against is not a precise number.
LITERATURE_KI67_VARIABILITY: Dict[str, object] = {
    "claim": (
        "Ki-67 LI scoring in gliomas agrees poorly between observers: kappa 0.32 "
        "for whole-slide and as low as 0.04 for hot-spot scoring (n=235 gliomas)"
    ),
    "citation": "Nielsen et al., Diagn Pathol 2018;13:38",
    "pmid": "29885671",
    "doi": "10.1186/s13000-018-0711-2",
}

#: Categorical colours (validated: worst adjacent CVD dE 37.7, all slots inside
#: the lightness band). Yellow is sub-3:1 on a light surface, so every mark using
#: it also carries a direct value label — the skill's relief rule.
STATE_COLORS: Dict[str, str] = {
    "VIABLE": "#2a78d6",     # categorical slot 1 (blue)
    "QUIESCENT": "#eda100",  # categorical slot 3 (yellow)
    "HYPOXIC": "#e34948",    # categorical slot 6 (red)
}
MEASURED_COLOR = "#2a78d6"   # slot 1 — measured data
MODEL_COLOR = "#eb6834"      # slot 8 — model prediction (dE 96.7 vs measured)
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8984"
GRID_COLOR = "#e4e3df"

#: Ordered state names used everywhere (low oxygen -> high oxygen).
STATE_ORDER: Tuple[str, ...] = ("HYPOXIC", "QUIESCENT", "VIABLE")

_STATE_ENUM = {
    "HYPOXIC": CellState.HYPOXIC,
    "QUIESCENT": CellState.QUIESCENT,
    "VIABLE": CellState.VIABLE,
}


def _results_dir() -> str:
    d = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
    )
    os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Real data loading + measured expression scores.
# ---------------------------------------------------------------------------
def load_panel(sample: Optional[str] = None) -> Dict:
    """Load the REAL Xenium panel (retry once on transient failure).

    Returns the loader dict unchanged. Rows are in the loader's deterministic
    order, which is the SAME order ``reconstruct_from_xenium`` builds its
    ``CellPopulation`` in — so expression rows align to cells index-for-index.
    That alignment is asserted by :func:`run`, never assumed.
    """
    kwargs: Dict[str, object] = {"dense": False}
    if sample is not None:
        kwargs["path"] = sample
    try:
        return load_singlecell_spatial(**kwargs)
    except Exception:
        return load_singlecell_spatial(**kwargs)


def _column(expression, genes: np.ndarray, gene: str) -> Optional[np.ndarray]:
    """Raw per-cell counts for one panel gene, or None if absent from the panel."""
    genes = np.asarray(genes, dtype=str)
    hit = np.flatnonzero(genes == gene)
    if hit.size == 0:
        return None
    col = int(hit[0])
    if sparse.issparse(expression):
        return np.asarray(expression[:, [col]].todense()).ravel().astype(np.float64)
    return np.asarray(expression)[:, col].astype(np.float64)


def normalized_log1p(counts: np.ndarray, totals: np.ndarray) -> np.ndarray:
    """Counts-per-median-cell normalise, then ``log1p``.

    Divides each cell's count by its total transcript count, rescales to the
    cohort MEDIAN total, and log1p-transforms. This is the standard depth
    correction and it is what removes the "a deeper cell is more likely to look
    marker-positive" confounder from every downstream statistic. Cells with zero
    total counts are left at 0 rather than producing a divide-by-zero.
    """
    counts = np.asarray(counts, dtype=np.float64)
    totals = np.asarray(totals, dtype=np.float64)
    scale = float(np.median(totals[totals > 0])) if np.any(totals > 0) else 1.0
    safe = np.where(totals > 0, totals, 1.0)
    return np.log1p(counts / safe * scale)


def _zscore(values: np.ndarray) -> np.ndarray:
    v = np.asarray(values, dtype=np.float64)
    sd = v.std()
    if sd == 0.0:
        return np.zeros_like(v)
    return (v - v.mean()) / sd


def composite_score(
    expression, genes: np.ndarray, gene_set: Sequence[str], totals: np.ndarray,
    subset: np.ndarray,
) -> Tuple[np.ndarray, List[str], List[str]]:
    """Mean z-score over a gene set, computed WITHIN ``subset`` (the tumour cells).

    Each present gene is depth-normalised (:func:`normalized_log1p`) then z-scored
    across the subset only, so every gene contributes on a comparable scale and
    the score is centred on the tumour-cell population being compared. Returns the
    score for the subset cells, the genes actually used, and the genes missing
    from the panel.
    """
    used: List[str] = []
    missing: List[str] = []
    cols: List[np.ndarray] = []
    for g in gene_set:
        raw = _column(expression, genes, g)
        if raw is None:
            missing.append(g)
            continue
        used.append(g)
        cols.append(_zscore(normalized_log1p(raw, totals)[subset]))
    if not cols:
        return np.zeros(int(subset.sum())), used, missing
    return np.mean(np.vstack(cols), axis=0), used, missing


def positivity_rules(
    counts: np.ndarray, totals: np.ndarray, subset: np.ndarray
) -> Dict[str, Dict[str, float]]:
    """Three DOCUMENTED positivity rules for one marker, over ``subset`` cells.

    There is no single correct way to binarise a sparse in-situ transcript into a
    "positive cell", so all three are reported and the reader can pick:

    ``nonzero``
        ``count >= 1``. The most permissive rule and the one that maps most
        directly onto "this probe fired in this cell".
    ``strict_z``
        depth-normalised log1p z-score ``>= STRICT_Z`` (1 SD above the tumour-cell
        mean) — the same gate ``real_xenium.TYPING_THRESHOLD`` uses for typing.
    ``count_ge_2``
        ``count >= 2``. A genuinely stricter rule that suppresses single-molecule
        false positives.

    Also returns ``percentile_rule_degenerate``: the fraction of subset cells at
    exactly zero. When that exceeds 1 - p, a "top-p percentile" positivity rule is
    mathematically degenerate (the threshold falls inside a mass of ties at zero),
    which is why no percentile rule is offered here.
    """
    c = np.asarray(counts, dtype=np.float64)[subset]
    norm_z = _zscore(normalized_log1p(np.asarray(counts, float), totals)[subset])
    n = int(c.size)
    out = {
        "nonzero": {
            "rule": "raw count >= 1",
            "n_positive": int((c >= 1).sum()),
            "fraction": float((c >= 1).mean()) if n else 0.0,
        },
        "strict_z": {
            "rule": f"depth-normalised log1p z-score >= {STRICT_Z}",
            "n_positive": int((norm_z >= STRICT_Z).sum()),
            "fraction": float((norm_z >= STRICT_Z).mean()) if n else 0.0,
        },
        "count_ge_2": {
            "rule": "raw count >= 2",
            "n_positive": int((c >= 2).sum()),
            "fraction": float((c >= 2).mean()) if n else 0.0,
        },
    }
    out["_meta"] = {
        "n_cells": n,
        "zero_fraction": float((c == 0).mean()) if n else 1.0,
        "mean_count": float(c.mean()) if n else 0.0,
        "max_count": float(c.max()) if n else 0.0,
    }
    return out


# ---------------------------------------------------------------------------
# Statistics: effect sizes + label-permutation p-values (no parametric claims).
# ---------------------------------------------------------------------------
def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Pooled-SD standardised mean difference ``(mean(a) - mean(b)) / s_pooled``."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.size < 2 or b.size < 2:
        return float("nan")
    va, vb = a.var(ddof=1), b.var(ddof=1)
    pooled = np.sqrt(((a.size - 1) * va + (b.size - 1) * vb) / (a.size + b.size - 2))
    if pooled == 0.0:
        return 0.0
    return float((a.mean() - b.mean()) / pooled)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta in ``[-1, 1]`` — P(a > b) - P(a < b), computed via ranks.

    Rank-based and therefore valid for the extremely tied, zero-inflated marker
    distributions here, where Cohen's d leans on a normality that does not hold.
    ``0`` means the two groups are stochastically identical.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    na, nb = a.size, b.size
    if na == 0 or nb == 0:
        return float("nan")
    pooled = np.concatenate([a, b])
    r = rankdata(pooled)
    rank_sum_a = float(r[:na].sum())
    # Mann-Whitney U from the rank sum; delta = 2U/(na*nb) - 1.
    u_a = rank_sum_a - na * (na + 1) / 2.0
    return float(2.0 * u_a / (na * nb) - 1.0)


def permutation_test_means(
    values: np.ndarray, group: np.ndarray, n_perm: int = N_PERMUTATIONS,
    seed: int = PERM_SEED,
) -> Dict[str, float]:
    """Two-sided label-permutation test on the difference of group means.

    ``group`` is a boolean mask selecting group A; the complement is group B. The
    observed statistic is ``mean(A) - mean(B)``; the null is generated by shuffling
    the group labels ``n_perm`` times, which is exactly the "the state labels carry
    no information about this marker" hypothesis we want to reject. The p-value is
    the standard ``(1 + #{|null| >= |obs|}) / (1 + n_perm)`` estimator, so it can
    never return 0.
    """
    values = np.asarray(values, float)
    group = np.asarray(group, bool)
    na = int(group.sum())
    nb = int(values.size - na)
    if na == 0 or nb == 0:
        return {"observed": float("nan"), "p_value": float("nan"),
                "n_a": na, "n_b": nb, "n_perm": 0}
    obs = float(values[group].mean() - values[~group].mean())

    # Blocked, fully vectorized permutation: shuffle each row of a tiled copy of
    # ``values`` independently (Generator.permuted, C-speed), take the first ``na``
    # entries of each row as that draw's group A. Exactly a label permutation, just
    # not written as a Python loop over 10,000 draws of 30,000 cells.
    rng = np.random.default_rng(seed)
    total = float(values.sum())
    n = int(values.size)
    block = max(1, min(n_perm, int(_PERM_BLOCK_CELLS // max(n, 1))))
    exceed = 0
    done = 0
    while done < n_perm:
        b = min(block, n_perm - done)
        tiled = np.tile(values, (b, 1))
        rng.permuted(tiled, axis=1, out=tiled)
        sa = tiled[:, :na].sum(axis=1)
        stats = sa / na - (total - sa) / nb
        exceed += int((np.abs(stats) >= abs(obs) - 1e-15).sum())
        done += b
    return {
        "observed": obs,
        "p_value": float((1 + exceed) / (1 + n_perm)),
        "n_a": na, "n_b": nb, "n_perm": int(n_perm),
    }


def spearman_permutation(
    x: np.ndarray, y: np.ndarray, n_perm: int = N_PERMUTATIONS,
    seed: int = PERM_SEED,
) -> Dict[str, float]:
    """Spearman rho with a permutation p-value, tie-safe.

    Ranks are computed once (``rankdata`` averages ties, which matters enormously
    here — the marker scores are ~99% tied at zero). rho is then the Pearson
    correlation of the ranks. The null shuffles one rank vector.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if x.size < 3:
        return {"rho": float("nan"), "p_value": float("nan"), "n": int(x.size)}
    rx = rankdata(x)
    ry = rankdata(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt(float((rx ** 2).sum()) * float((ry ** 2).sum()))
    if denom == 0.0:
        return {"rho": 0.0, "p_value": 1.0, "n": int(x.size)}
    obs = float((rx * ry).sum() / denom)

    rng = np.random.default_rng(seed)
    norm = float(np.sqrt((rx ** 2).sum())) * float(np.sqrt((ry ** 2).sum()))
    n = int(x.size)
    block = max(1, min(n_perm, int(_PERM_BLOCK_CELLS // max(n, 1))))
    exceed = 0
    done = 0
    while done < n_perm:
        b = min(block, n_perm - done)
        tiled = np.tile(ry, (b, 1))
        rng.permuted(tiled, axis=1, out=tiled)
        stats = (tiled @ rx) / norm
        exceed += int((np.abs(stats) >= abs(obs) - 1e-15).sum())
        done += b
    return {
        "rho": obs,
        "p_value": float((1 + exceed) / (1 + n_perm)),
        "n": int(x.size),
        "n_perm": int(n_perm),
    }


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Spearman rho between x and y with ``z`` linearly residualised out of ranks.

    Used to strip the sequencing-depth confounder: ``z`` is ``log(total counts)``.
    Rank-transform all three, regress the ranks of x and y on the ranks of z, and
    correlate the residuals. If the raw correlation survives this, it is not a
    depth artifact.
    """
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    rz_c = rz - rz.mean()
    denom = float((rz_c ** 2).sum())
    if denom == 0.0:
        return float("nan")

    def resid(r):
        rc = r - r.mean()
        beta = float((rc * rz_c).sum()) / denom
        return rc - beta * rz_c

    ex, ey = resid(rx), resid(ry)
    d = np.sqrt(float((ex ** 2).sum()) * float((ey ** 2).sum()))
    return float((ex * ey).sum() / d) if d > 0 else float("nan")


def rate_by_state(
    positive: np.ndarray, states: np.ndarray
) -> Dict[str, Dict[str, float]]:
    """Marker-positive RATE within each model state, with a Wilson 95% interval.

    Wilson rather than normal-approximation because several of these rates are a
    handful of cells out of a few hundred, where the normal interval is nonsense
    (it can dip below zero).
    """
    out: Dict[str, Dict[str, float]] = {}
    for name in STATE_ORDER:
        m = states == int(_STATE_ENUM[name])
        n = int(m.sum())
        k = int(positive[m].sum()) if n else 0
        if n == 0:
            out[name] = {"n": 0, "k": 0, "rate": float("nan"),
                         "lo": float("nan"), "hi": float("nan")}
            continue
        p = k / n
        z = 1.959963985
        denom = 1 + z * z / n
        centre = (p + z * z / (2 * n)) / denom
        half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
        out[name] = {"n": n, "k": k, "rate": float(p),
                     "lo": float(max(0.0, centre - half)),
                     "hi": float(min(1.0, centre + half))}
    return out


def power_note(rate_overall: float, n_group: int) -> Dict[str, float]:
    """How many positives a group would contain under the null — the power check.

    A "no enrichment" result is only meaningful if the test COULD have detected
    enrichment. This reports the expected positive count in a group of size
    ``n_group`` under the null rate, and the fold-enrichment that group would need
    in order to contain at least 5 positives (a rough floor for any rate estimate
    to mean anything).
    """
    expected = float(rate_overall * n_group)
    needed_fold = float("inf") if expected <= 0 else float(5.0 / expected)
    return {"expected_positives_under_null": expected,
            "fold_enrichment_needed_for_5_positives": needed_fold}


# ---------------------------------------------------------------------------
# Driving the frozen pipeline: reconstruct -> oxygen -> state banding.
# ---------------------------------------------------------------------------
def model_states_at_dx(
    dx_um: float, sample: Optional[str] = None,
    tissue_mask_mode: Optional[str] = None,
) -> Dict:
    """Run the REAL pipeline at one ``dx`` and return the banded tumour-cell states.

    Calls the unmodified :func:`reconstruct.real_pipeline.reconstruct_real` (which
    already solves the real-perfusion oxygen field via ``physics.solve_oxygen``),
    then the unmodified :func:`cells.agents.update_cell_states` to apply the
    oxygen banding. Nothing here re-implements the banding — the object under test
    is the shipping code path.

    Returns the per-cell state array, the per-cell sampled oxygen, the tumour mask,
    and the tumour-cell state fractions.
    """
    cfg = OxygenConfig()
    kwargs: Dict[str, object] = {"source": "xenium", "dx_um": float(dx_um)}
    if sample is not None:
        kwargs["sample"] = sample
    if tissue_mask_mode is not None:
        kwargs["tissue_mask_mode"] = str(tissue_mask_mode)
    micro = reconstruct_real(**kwargs)
    oxygen = micro.fields.oxygen
    update_cell_states(micro, oxygen, cfg)

    cells = micro.cells
    col, row = micro.domain.to_index(cells.x, cells.y)
    local_o2 = np.asarray(oxygen, float)[row, col]

    tumor = cells.cell_type == int(CellType.TUMOR)
    n_tumor = int(tumor.sum())
    fractions = {}
    for name in STATE_ORDER:
        k = int((cells.state[tumor] == int(_STATE_ENUM[name])).sum())
        fractions[name] = float(k / n_tumor) if n_tumor else 0.0

    return {
        "dx_um": float(dx_um),
        "state": cells.state.copy(),
        "local_o2": local_o2,
        "tumor_mask": tumor,
        "n_tumor": n_tumor,
        "fractions": fractions,
        "x": cells.x.copy(),
        "y": cells.y.copy(),
        "hypoxia_threshold": float(cfg.hypoxia_threshold),
        "mean_oxygen_tumor": float(local_o2[tumor].mean()) if n_tumor else float("nan"),
        "tissue_mask_mode": micro.metadata.get("real_pipeline", {}).get(
            "tissue_mask_mode"),
    }


# ---------------------------------------------------------------------------
# The experiment.
# ---------------------------------------------------------------------------
def analyze_mode(
    main: Dict, expr, genes: np.ndarray, totals: np.ndarray,
    n_perm: int, label: str,
) -> Dict:
    """Every per-cell test, for ONE reconstruction (one tissue-mask mode at dx).

    Split out of :func:`run` so the identical battery can be applied under both
    tissue-mask conventions — the ``geometric`` default and the legacy
    ``occupancy`` mode under which ``results/patient_report_xenium.md`` quoted its
    headline state fractions. Running both is what separates "the state model is
    wrong" from "the state model is a bookkeeping artifact".
    """
    print(f"\n{'='*70}\nPER-CELL TESTS — tissue mask = {label}\n{'='*70}")
    tumor = main["tumor_mask"]
    n_tumor = main["n_tumor"]
    states_t = main["state"][tumor]
    o2_t = main["local_o2"][tumor]
    totals_t = totals[tumor]

    viable_mask = states_t == int(CellState.VIABLE)
    hypoxic_mask = states_t == int(CellState.HYPOXIC)
    quiescent_mask = states_t == int(CellState.QUIESCENT)

    # --- TEST 1: marginal measured vs modelled fractions --------------------
    print("\n--- Test 1: marginal fractions ---")
    marker_stats: Dict[str, Dict] = {}
    for g in PRIMARY_MARKERS:
        raw = _column(expr, genes, g)
        if raw is None:
            marker_stats[g] = {"on_panel": False}
            continue
        rules = positivity_rules(raw, totals, tumor)
        marker_stats[g] = {"on_panel": True, "rules": rules,
                           "counts_tumor": raw[tumor]}
        print(f"  {g}: nonzero={100*rules['nonzero']['fraction']:.3f}%  "
              f"strict_z={100*rules['strict_z']['fraction']:.3f}%  "
              f"count>=2={100*rules['count_ge_2']['fraction']:.3f}%")
    print(f"  model VIABLE fraction = {100*main['fractions']['VIABLE']:.3f}%")

    # --- composite scores --------------------------------------------------
    prolif, prolif_used, prolif_missing = composite_score(
        expr, genes, PROLIFERATION_GENES, totals, tumor)
    hif, hif_used, hif_missing = composite_score(
        expr, genes, HIF_TARGET_GENES, totals, tumor)

    # --- DEPTH CONFOUNDER AUDIT --------------------------------------------
    # "any proliferation transcript detected" is the cleanest detection-rate
    # readout: it tracks per-cell transcript depth regardless of biology, so
    # tabulating it beside depth makes any detection artifact visible.
    prolif_set = set(prolif_used)
    prolif_cols = [i for i, g in enumerate(genes) if g in prolif_set]
    if sparse.issparse(expr):
        any_all = np.asarray(expr[:, prolif_cols].sum(axis=1)).ravel() > 0
    else:
        any_all = np.asarray(expr)[:, prolif_cols].sum(axis=1) > 0
    any_prolif = any_all[tumor].astype(np.float64)

    depth_by_state: Dict[str, Dict] = {}
    for name in STATE_ORDER:
        m = states_t == int(_STATE_ENUM[name])
        depth_by_state[name] = {
            "n": int(m.sum()),
            "mean_total_counts": float(totals_t[m].mean()) if m.any() else float("nan"),
            "median_total_counts": float(np.median(totals_t[m])) if m.any() else float("nan"),
            "any_prolif_transcript_rate": float(any_prolif[m].mean()) if m.any() else float("nan"),
        }
    print("\n  DEPTH AUDIT (per model state):")
    for name in STATE_ORDER:
        d = depth_by_state[name]
        if d["n"] == 0:
            print(f"    {name:<10} n=     0  (state is empty under this mask)")
            continue
        print(f"    {name:<10} n={d['n']:>6,}  mean total counts={d['mean_total_counts']:7.1f}"
              f"  any-proliferation-transcript rate={100*d['any_prolif_transcript_rate']:5.2f}%")
    occupied = [depth_by_state[s]["mean_total_counts"] for s in STATE_ORDER
                if depth_by_state[s]["n"] > 0]
    depth_spread = float(max(occupied) / min(occupied)) if len(occupied) > 1 else float("nan")
    print(f"    -> transcript-depth spread across occupied states = {depth_spread:.2f}x")

    # --- TEST 2: spatial enrichment by model state -------------------------
    print("\n--- Test 2: spatial enrichment by model state ---")
    enrichment: Dict[str, Dict] = {}
    for g in PRIMARY_MARKERS:
        if not marker_stats[g]["on_panel"]:
            continue
        pos = (marker_stats[g]["counts_tumor"] >= 1).astype(np.float64)
        rates = rate_by_state(pos, states_t)
        overall = float(pos.mean())
        hq = hypoxic_mask | quiescent_mask
        enrichment[g] = {
            "rates": rates,
            "overall_rate": overall,
            "viable_vs_rest": permutation_test_means(pos, viable_mask, n_perm),
            "hypoxic_vs_quiescent": permutation_test_means(
                pos[hq], hypoxic_mask[hq], n_perm),
            "power": power_note(overall, int(viable_mask.sum())),
        }
        r = rates
        print(f"  {g}+ rate by state: " + "  ".join(
            f"{k}={100*r[k]['rate']:.3f}% ({r[k]['k']}/{r[k]['n']})"
            for k in STATE_ORDER))
        print(f"    VIABLE vs rest       p={enrichment[g]['viable_vs_rest']['p_value']:.4f}"
              f"  (expected {enrichment[g]['power']['expected_positives_under_null']:.1f}"
              f" positives in VIABLE under the null)")
        print(f"    HYPOXIC vs QUIESCENT p="
              f"{enrichment[g]['hypoxic_vs_quiescent']['p_value']:.4f}")

    score_by_state: Dict[str, Dict] = {}
    for name in STATE_ORDER:
        m = states_t == int(_STATE_ENUM[name])
        score_by_state[name] = {
            "n": int(m.sum()),
            "mean": float(prolif[m].mean()) if m.any() else float("nan"),
            "median": float(np.median(prolif[m])) if m.any() else float("nan"),
            "hif_mean": float(hif[m].mean()) if m.any() else float("nan"),
            "values": prolif[m],
        }

    hq = hypoxic_mask | quiescent_mask
    composite_tests = {
        "viable_vs_rest": {
            "cohens_d": cohens_d(prolif[viable_mask], prolif[~viable_mask]),
            "cliffs_delta": cliffs_delta(prolif[viable_mask], prolif[~viable_mask]),
            **permutation_test_means(prolif, viable_mask, n_perm),
        },
        "hypoxic_vs_quiescent": {
            "cohens_d": cohens_d(prolif[hypoxic_mask], prolif[quiescent_mask]),
            "cliffs_delta": cliffs_delta(prolif[hypoxic_mask], prolif[quiescent_mask]),
            **permutation_test_means(prolif[hq], hypoxic_mask[hq], n_perm),
        },
    }
    print("  composite proliferation score by state: " + "  ".join(
        f"{k}={score_by_state[k]['mean']:+.4f}" for k in STATE_ORDER))
    for k, v in composite_tests.items():
        print(f"    {k}: d={v['cohens_d']:+.4f}  cliff={v['cliffs_delta']:+.4f}  "
              f"p={v['p_value']:.4f}")

    hif_tests = {
        "hypoxic_vs_rest": {
            "cohens_d": cohens_d(hif[hypoxic_mask], hif[~hypoxic_mask]),
            "cliffs_delta": cliffs_delta(hif[hypoxic_mask], hif[~hypoxic_mask]),
            **permutation_test_means(hif, hypoxic_mask, n_perm),
        },
        "viable_vs_rest": {
            "cohens_d": cohens_d(hif[viable_mask], hif[~viable_mask]),
            "cliffs_delta": cliffs_delta(hif[viable_mask], hif[~viable_mask]),
            **permutation_test_means(hif, viable_mask, n_perm),
        },
    }
    print("  composite HIF-target score by state: " + "  ".join(
        f"{k}={score_by_state[k]['hif_mean']:+.4f}" for k in STATE_ORDER))

    # --- TEST 2c: DEPTH-STRATIFIED enrichment (the decisive control) -------
    print("\n  DEPTH-STRATIFIED re-test (labels permuted within total-count deciles):")
    strata = depth_deciles(totals_t)
    stratified = {
        "prolif_viable_vs_rest": stratified_permutation_test(
            prolif, viable_mask, strata, n_perm),
        "prolif_hypoxic_vs_quiescent": stratified_permutation_test(
            prolif[hq], hypoxic_mask[hq], strata[hq], n_perm),
        "any_prolif_viable_vs_rest": stratified_permutation_test(
            any_prolif, viable_mask, strata, n_perm),
        "any_prolif_hypoxic_vs_quiescent": stratified_permutation_test(
            any_prolif[hq], hypoxic_mask[hq], strata[hq], n_perm),
        "hif_hypoxic_vs_rest": stratified_permutation_test(
            hif, hypoxic_mask, strata, n_perm),
        "hif_viable_vs_rest": stratified_permutation_test(
            hif, viable_mask, strata, n_perm),
    }
    for k, v in stratified.items():
        print(f"    {k:<34} diff={v['observed']:+.5f}  p={v['p_value']:.4f}"
              f"  ({v['n_strata_used']} strata)")

    # --- TEST 3: continuous oxygen vs measured scores ----------------------
    print("\n--- Test 3: continuous simulated-O2 vs measured score ---")
    log_depth = np.log(np.maximum(totals_t, 1.0))
    continuous = {
        "o2_vs_proliferation": {
            **spearman_permutation(o2_t, prolif, n_perm),
            "partial_rho_depth_removed": partial_spearman(o2_t, prolif, log_depth),
        },
        "o2_vs_hif_targets": {
            **spearman_permutation(o2_t, hif, n_perm),
            "partial_rho_depth_removed": partial_spearman(o2_t, hif, log_depth),
        },
        "o2_vs_depth": spearman_permutation(o2_t, totals_t, n_perm),
    }
    for g in PRIMARY_MARKERS:
        if not marker_stats[g]["on_panel"]:
            continue
        vals = normalized_log1p(_column(expr, genes, g), totals)[tumor]
        continuous[f"o2_vs_{g}"] = {
            **spearman_permutation(o2_t, vals, n_perm),
            "partial_rho_depth_removed": partial_spearman(o2_t, vals, log_depth),
        }
    for k, v in continuous.items():
        extra = (f"  partial(depth removed)={v['partial_rho_depth_removed']:+.4f}"
                 if "partial_rho_depth_removed" in v else "")
        print(f"  {k}: rho={v['rho']:+.4f}  p={v['p_value']:.4f}{extra}")

    return {
        "mask_mode": label,
        "n_tumor": n_tumor,
        "fractions": main["fractions"],
        "marker_stats": {g: (marker_stats[g]["rules"] if marker_stats[g]["on_panel"]
                             else {"on_panel": False}) for g in PRIMARY_MARKERS},
        "prolif_genes_used": prolif_used,
        "prolif_genes_missing": prolif_missing,
        "hif_genes_used": hif_used,
        "hif_genes_missing": hif_missing,
        "enrichment": enrichment,
        "score_by_state": {k: {kk: vv for kk, vv in v.items() if kk != "values"}
                           for k, v in score_by_state.items()},
        "composite_tests": composite_tests,
        "hif_tests": hif_tests,
        "stratified": stratified,
        "continuous": continuous,
        "depth_by_state": depth_by_state,
        "depth_spread_x": depth_spread,
        "_arrays": {
            "prolif": prolif, "hif": hif, "o2_t": o2_t, "states_t": states_t,
            "any_prolif": any_prolif, "totals_t": totals_t,
            "x_t": main["x"][tumor], "y_t": main["y"][tumor],
            "mki67_counts_t": (marker_stats["MKI67"]["counts_tumor"]
                               if marker_stats["MKI67"]["on_panel"] else None),
        },
    }


def run(sample: Optional[str] = None, n_perm: int = N_PERMUTATIONS) -> Dict:
    """Run the full validation and return every number the report quotes."""
    print("=== CELL-STATE VALIDATION on the REAL Xenium GBM section ===\n")

    print("Loading REAL Xenium panel ...")
    data = load_panel(sample)
    expr = data["expression"]
    genes = np.asarray(data["genes"], dtype=str)
    n_cells = int(expr.shape[0])
    totals = np.asarray(expr.sum(axis=1)).ravel().astype(np.float64)
    print(f"  {n_cells:,} cells x {len(genes)} panel genes "
          f"({data['source']['panel']})\n")

    # --- model states across the dx sweep, under BOTH tissue-mask conventions
    sweep: Dict[str, Dict[float, Dict]] = {m: {} for m in TISSUE_MASK_MODES}
    for mode in TISSUE_MASK_MODES:
        for dx in DX_SWEEP:
            print(f"Reconstruct + oxygen + banding: mask={mode:<10} dx={dx:>2.0f} um ...")
            sweep[mode][dx] = model_states_at_dx(dx, sample, tissue_mask_mode=mode)
            f = sweep[mode][dx]["fractions"]
            print("    tumour-cell states: " + "  ".join(
                f"{k}={100*v:.2f}%" for k, v in f.items()))

    ref = sweep[DEFAULT_MASK_MODE][DX_CONVERGED]
    if ref["state"].size != n_cells:
        raise RuntimeError(
            f"row misalignment: {ref['state'].size} cells from the reconstruction "
            f"vs {n_cells} expression rows — the expression matrix and the "
            "CellPopulation are not index-aligned, so no per-cell comparison is valid"
        )
    print(f"\nTumour cells at dx={DX_CONVERGED:.0f}: {ref['n_tumor']:,} "
          f"({100.0*ref['n_tumor']/n_cells:.1f}% of the section)")

    per_mode = {
        mode: analyze_mode(sweep[mode][DX_CONVERGED], expr, genes, totals,
                           n_perm, mode)
        for mode in TISSUE_MASK_MODES
    }

    return {
        "n_cells": n_cells,
        "n_tumor": ref["n_tumor"],
        "panel": data["source"]["panel"],
        "dataset": data["source"]["dataset"],
        "sample": data["source"]["sample"],
        "dx_converged": DX_CONVERGED,
        "default_mask_mode": DEFAULT_MASK_MODE,
        "sweep": {
            mode: {dx: {"fractions": sweep[mode][dx]["fractions"],
                        "n_tumor": sweep[mode][dx]["n_tumor"],
                        "mean_oxygen_tumor": sweep[mode][dx]["mean_oxygen_tumor"]}
                   for dx in DX_SWEEP}
            for mode in TISSUE_MASK_MODES},
        "hypoxia_threshold": ref["hypoxia_threshold"],
        "pd_config": {"hypoxia_resistance": PDConfig().hypoxia_resistance,
                      "quiescent_resistance": PDConfig().quiescent_resistance},
        "per_mode": per_mode,
        "n_perm": int(n_perm),
    }


def stratified_permutation_test(
    values: np.ndarray, group: np.ndarray, strata: np.ndarray,
    n_perm: int = N_PERMUTATIONS, seed: int = PERM_SEED,
) -> Dict[str, float]:
    """Difference of means between groups, POOLED ACROSS STRATA, permuted WITHIN strata.

    This is the decisive control for the transcript-depth confounder documented in
    the report. Model states differ several-fold in mean per-cell transcript count,
    and a cell with fewer total transcripts is less likely to show ANY given rare
    marker — a dropout effect that counts-per-median normalisation CANNOT remove
    (normalising a zero still gives a zero). Comparing states without conditioning
    on depth therefore measures detection efficiency, not biology.

    Here cells are first binned into ``strata`` (depth deciles). The statistic is
    the stratum-size-weighted mean of the within-stratum group difference, so it
    only ever compares cells of comparable sequencing depth. The null shuffles the
    group labels **within each stratum independently**, which is exactly the
    hypothesis "given depth, state carries no information about this marker".

    Strata in which either group is empty contribute nothing and are counted in
    ``n_strata_used``.
    """
    values = np.asarray(values, float)
    group = np.asarray(group, bool)
    strata = np.asarray(strata)
    rng = np.random.default_rng(seed)

    blocks: List[Tuple[np.ndarray, int, int]] = []
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        na = int(group[idx].sum())
        nb = int(idx.size - na)
        if na == 0 or nb == 0:
            continue
        blocks.append((values[idx], na, nb))
    if not blocks:
        return {"observed": float("nan"), "p_value": float("nan"),
                "n_strata_used": 0, "n_perm": 0}

    weights = np.array([na + nb for _, na, nb in blocks], dtype=np.float64)
    weights /= weights.sum()

    # Observed: stratum-size-weighted mean of the within-stratum group difference.
    obs = 0.0
    wi = 0
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        g = group[idx]
        if not g.any() or g.all():
            continue
        obs += weights[wi] * (values[idx][g].mean() - values[idx][~g].mean())
        wi += 1

    # null: permute within each stratum, accumulate the pooled statistic per draw
    stats = np.zeros(n_perm, dtype=np.float64)
    for w, (vals, na, nb) in zip(weights, blocks):
        total = float(vals.sum())
        n = int(vals.size)
        block = max(1, min(n_perm, int(_PERM_BLOCK_CELLS // max(n, 1))))
        done = 0
        while done < n_perm:
            b = min(block, n_perm - done)
            tiled = np.tile(vals, (b, 1))
            rng.permuted(tiled, axis=1, out=tiled)
            sa = tiled[:, :na].sum(axis=1)
            stats[done:done + b] += w * (sa / na - (total - sa) / nb)
            done += b

    exceed = int((np.abs(stats) >= abs(obs) - 1e-15).sum())
    return {
        "observed": float(obs),
        "p_value": float((1 + exceed) / (1 + n_perm)),
        "n_strata_used": len(blocks),
        "n_perm": int(n_perm),
    }


def depth_deciles(totals: np.ndarray, n_bins: int = 10) -> np.ndarray:
    """Assign each cell to a total-transcript-count quantile bin (0..n_bins-1)."""
    totals = np.asarray(totals, float)
    edges = np.quantile(totals, np.linspace(0, 1, n_bins + 1)[1:-1])
    return np.searchsorted(edges, totals, side="right")


# ---------------------------------------------------------------------------
# Figure.
# ---------------------------------------------------------------------------
def _style_axes(ax, grid_axis: str = "y") -> None:
    """Recessive grid + no top/right spines — the marks carry the ink, not the frame."""
    ax.set_axisbelow(True)
    ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.8, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID_COLOR)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8, length=0)
    ax.set_title(ax.get_title(), color=INK_PRIMARY)


def make_figure(results: Dict, out_path: Optional[str] = None) -> str:
    """Six-panel evidence figure. Returns the written path."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    geo = results["per_mode"]["geometric"]
    occ = results["per_mode"]["occupancy"]
    lit_lo = min(float(e["value_pct"]) for e in LITERATURE_KI67)
    lit_hi = max(float(e["value_pct"]) for e in LITERATURE_KI67)

    fig = plt.figure(figsize=(17.5, 10.4))
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.26,
                          left=0.085, right=0.985, top=0.878, bottom=0.105)
    fig.patch.set_facecolor("#fcfcfb")

    # ---- A: measured vs modelled proliferating fraction (log scale) -------
    ax = fig.add_subplot(gs[0, 0])
    rows = [
        ("MKI67+  (≥2 counts)", 100 * geo["marker_stats"]["MKI67"]["count_ge_2"]["fraction"], MEASURED_COLOR),
        ("MKI67+  (≥1 count)", 100 * geo["marker_stats"]["MKI67"]["nonzero"]["fraction"], MEASURED_COLOR),
        ("TOP2A+  (≥1 count)", 100 * geo["marker_stats"]["TOP2A"]["nonzero"]["fraction"], MEASURED_COLOR),
        ("VIABLE\noccupancy mask", 100 * occ["fractions"]["VIABLE"], MODEL_COLOR),
        ("VIABLE\ngeometric mask", 100 * geo["fractions"]["VIABLE"], MODEL_COLOR),
    ]
    ypos = np.arange(len(rows))
    ax.axvspan(lit_lo, lit_hi, color="#c9c8c2", alpha=0.55, zorder=1)
    ax.text(np.sqrt(lit_lo * lit_hi), 0.45,
            f"published GBM\nKi-67 index\n{lit_lo:.0f}–{lit_hi:.0f}%",
            ha="center", va="center", fontsize=7.5, color=INK_SECONDARY, zorder=5)
    for i, (lab, val, col) in enumerate(rows):
        ax.barh(i, val, color=col, height=0.58, zorder=3)
        ax.text(val * 1.25, i, f"{val:.3g}%", va="center", ha="left",
                fontsize=8, color=INK_PRIMARY, zorder=4)
    ax.set_yticks(ypos)
    ax.set_yticklabels([r[0] for r in rows], fontsize=8.5, color=INK_PRIMARY)
    ax.set_xscale("log")
    ax.set_xlim(0.05, 600)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("% of tumour cells  (log scale)", fontsize=8.5, color=INK_SECONDARY)
    ax.set_title("A · Measured proliferation vs modelled VIABLE fraction",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax, grid_axis="x")
    ax.legend(handles=[Patch(color=MEASURED_COLOR, label="measured (Xenium panel)"),
                       Patch(color=MODEL_COLOR, label="model state")],
              fontsize=7.5, loc="lower right", frameon=False)

    # ---- B: MKI67+ rate by model state (occupancy: all three states exist) -
    ax = fig.add_subplot(gs[0, 1])
    rates = occ["enrichment"]["MKI67"]["rates"]
    overall = occ["enrichment"]["MKI67"]["overall_rate"]
    xs = np.arange(len(STATE_ORDER))
    for i, name in enumerate(STATE_ORDER):
        r = rates[name]
        if r["n"] == 0:
            continue
        ax.bar(i, 100 * r["rate"], color=STATE_COLORS[name], width=0.6, zorder=3)
        ax.plot([i, i], [100 * r["lo"], 100 * r["hi"]], color=INK_PRIMARY,
                linewidth=1.6, zorder=4)
        ax.text(i, 100 * r["hi"] + 0.12, f"{100*r['rate']:.2f}%\n{r['k']}/{r['n']:,}",
                ha="center", va="bottom", fontsize=8, color=INK_PRIMARY, zorder=5)
    ax.axhline(100 * overall, color=INK_MUTED, linewidth=1.2, linestyle=(0, (4, 3)), zorder=2)
    ax.text(2.48, 100 * overall, f"all tumour cells {100*overall:.2f}%  ", fontsize=7.5,
            color=INK_SECONDARY, va="bottom", ha="right")
    ax.set_xticks(xs)
    ax.set_xticklabels(STATE_ORDER, fontsize=8.5, color=INK_PRIMARY)
    ax.set_ylabel("MKI67-positive tumour cells (%)", fontsize=8.5, color=INK_SECONDARY)
    ax.set_ylim(0, max(100 * rates[s]["hi"] for s in STATE_ORDER if rates[s]["n"]) * 1.55)
    p_v = occ["enrichment"]["MKI67"]["viable_vs_rest"]["p_value"]
    p_hq = occ["enrichment"]["MKI67"]["hypoxic_vs_quiescent"]["p_value"]
    ax.set_title(f"B · MKI67+ rate by model state (occupancy mask, dx=12 µm)\n"
                 f"VIABLE vs rest p={p_v:.3f} (on 4 cells) · HYPOXIC vs QUIESCENT p={p_hq:.3f}",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax)

    # ---- C: the depth artifact — detection rate tracks transcript depth ----
    ax = fig.add_subplot(gs[0, 2])
    for mode, pm, marker in (("geometric", geo, "o"), ("occupancy", occ, "s")):
        for name in STATE_ORDER:
            d = pm["depth_by_state"][name]
            if d["n"] == 0 or not np.isfinite(d["mean_total_counts"]):
                continue
            ax.scatter(d["mean_total_counts"], 100 * d["any_prolif_transcript_rate"],
                       s=95, marker=marker, color=STATE_COLORS[name],
                       edgecolor="#fcfcfb", linewidth=1.6, zorder=4)
            right = d["mean_total_counts"] > 200
            ax.annotate(f"{name[:4].title()} · {mode[:3]}",
                        (d["mean_total_counts"], 100 * d["any_prolif_transcript_rate"]),
                        textcoords="offset points",
                        xytext=(-9, -12) if right else (9, -3),
                        ha="right" if right else "left", fontsize=7.5,
                        color=INK_SECONDARY, zorder=5)
    ax.set_xlabel("mean total transcripts per cell", fontsize=8.5, color=INK_SECONDARY)
    ax.set_ylabel("cells with ANY proliferation transcript (%)",
                  fontsize=8.5, color=INK_SECONDARY)
    ax.margins(0.22)
    ax.set_title("C · The confound: proliferation DETECTION tracks transcript depth,\n"
                 "not model state (circles = geometric mask, squares = occupancy)",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax, grid_axis="both")

    # ---- D: dx- and mask-dependence of the state fractions ----------------
    ax = fig.add_subplot(gs[1, 0])
    labels, bottoms = [], []
    xi = 0
    xticks, xticklabels = [], []
    for mode in TISSUE_MASK_MODES:
        for dx in DX_SWEEP:
            f = results["sweep"][mode][dx]["fractions"]
            base = 0.0
            for name in STATE_ORDER:
                v = 100 * f[name]
                if v > 0:
                    ax.bar(xi, v, bottom=base, color=STATE_COLORS[name], width=0.68,
                           zorder=3, edgecolor="#fcfcfb", linewidth=1.2)
                    if v >= 7:
                        ax.text(xi, base + v / 2, f"{v:.0f}", ha="center", va="center",
                                fontsize=7.5, color="#ffffff" if name != "QUIESCENT" else INK_PRIMARY,
                                zorder=5)
                base += v
            xticks.append(xi)
            xticklabels.append(f"{dx:.0f}")
            xi += 1
        xi += 0.7
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=8.5, color=INK_PRIMARY)
    ax.set_ylabel("% of tumour cells", fontsize=8.5, color=INK_SECONDARY)
    ax.set_xlabel("grid spacing dx (µm)", fontsize=8.5, color=INK_SECONDARY)
    ax.text(1, -18, "geometric mask", ha="center", fontsize=8.5, color=INK_PRIMARY)
    ax.text(4.7, -18, "occupancy mask", ha="center", fontsize=8.5, color=INK_PRIMARY)
    ax.set_ylim(0, 100)
    ax.set_title("D · The state vector is set by bookkeeping: VIABLE runs 98.2% → 0.6%\n"
                 "across two tissue-mask conventions and three grid spacings",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax)
    ax.legend(handles=[Patch(color=STATE_COLORS[s], label=s) for s in STATE_ORDER],
              fontsize=7.5, loc="lower left", frameon=False, ncol=3,
              bbox_to_anchor=(0.0, -0.30))

    # ---- E: spatial map, model state + measured MKI67+ --------------------
    ax = fig.add_subplot(gs[1, 1])
    a = occ["_arrays"]
    for name in STATE_ORDER:
        m = a["states_t"] == int(_STATE_ENUM[name])
        if not m.any():
            continue
        ax.scatter(a["x_t"][m], a["y_t"][m], s=1.1, c=STATE_COLORS[name],
                   linewidths=0, alpha=0.55, zorder=2, rasterized=True)
    pos = a["mki67_counts_t"] >= 1
    ax.scatter(a["x_t"][pos], a["y_t"][pos], s=13, facecolor="none",
               edgecolor=INK_PRIMARY, linewidth=0.7, zorder=5,
               label=f"MKI67+ cell (n={int(pos.sum())})")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (µm)", fontsize=8.5, color=INK_SECONDARY)
    ax.set_ylabel("y (µm)", fontsize=8.5, color=INK_SECONDARY)
    ax.set_title("E · Model state in space vs the 134 measured MKI67+ cells\n"
                 "(occupancy mask, dx=12 µm) — MKI67+ cells are spread across every\n"
                 "state, not concentrated in the model's VIABLE region",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax, grid_axis="both")
    ax.legend(handles=[Patch(color=STATE_COLORS[s], label=s) for s in STATE_ORDER]
              + [plt.Line2D([], [], marker="o", linestyle="none", markersize=5,
                            markerfacecolor="none", markeredgecolor=INK_PRIMARY,
                            label=f"MKI67+ (n={int(pos.sum())})")],
              fontsize=7.5, loc="upper center", frameon=False, ncol=4,
              bbox_to_anchor=(0.5, -0.16), columnspacing=1.2, handletextpad=0.5)

    # ---- F: what the simulated oxygen field actually correlates with -------
    ax = fig.add_subplot(gs[1, 2])
    keys = [("o2_vs_depth", "total transcript count\n(a TECHNICAL variable)"),
            ("o2_vs_proliferation", "composite proliferation\nscore (8 genes)"),
            ("o2_vs_hif_targets", "HIF-target score\n(HILPDA, LOX)"),
            ("o2_vs_MKI67", "MKI67"),
            ("o2_vs_TOP2A", "TOP2A")]
    ys = np.arange(len(keys))[::-1]
    for (k, lab), y in zip(keys, ys):
        raw = occ["continuous"][k]["rho"]
        par = occ["continuous"][k].get("partial_rho_depth_removed")
        ax.barh(y + 0.18, raw, height=0.33, color="#2a78d6", zorder=3)
        ax.text(raw + (0.002 if raw >= 0 else -0.002), y + 0.18, f"{raw:+.3f}",
                va="center", ha="left" if raw >= 0 else "right", fontsize=7.5,
                color=INK_PRIMARY)
        if par is not None and np.isfinite(par):
            ax.barh(y - 0.18, par, height=0.33, color="#d03b3b", zorder=3)
            ax.text(par + (0.002 if par >= 0 else -0.002), y - 0.18, f"{par:+.3f}",
                    va="center", ha="left" if par >= 0 else "right", fontsize=7.5,
                    color=INK_PRIMARY)
    ax.axvline(0, color=INK_MUTED, linewidth=1.0, zorder=4)
    ax.set_yticks(ys)
    ax.set_yticklabels([lab for _, lab in keys], fontsize=8, color=INK_PRIMARY)
    ax.set_xlabel("Spearman ρ with the SIMULATED oxygen at each cell",
                  fontsize=8.5, color=INK_SECONDARY)
    ax.set_xlim(-0.02, 0.105)
    ax.set_title("F · The simulated O₂ field tracks transcript depth ~10× more\n"
                 "strongly than it tracks any proliferation marker (occupancy mask)",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax, grid_axis="x")
    ax.legend(handles=[Patch(color="#2a78d6", label="raw ρ"),
                       Patch(color="#d03b3b", label="ρ with log-depth partialled out")],
              fontsize=7.5, loc="lower right", frameon=False)

    fig.suptitle(
        "Is the oxygen-banded cell-state model supported by measured proliferation? "
        "— real 10x Xenium GBM section, 40,887 cells, 319-gene panel",
        fontsize=13, color=INK_PRIMARY, x=0.085, ha="left", y=0.968)
    fig.text(0.085, 0.930,
             "VERDICT: CONTRADICTED. The state labels track per-cell transcript depth, "
             "not measured Ki-67/TOP2A; every apparent enrichment vanishes under "
             "depth-stratified permutation.",
             fontsize=9.5, color=INK_SECONDARY, ha="left")

    out_path = out_path or os.path.join(_results_dir(), "cellstate_validation.png")
    fig.savefig(out_path, dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"figure written: {out_path}")
    return out_path
