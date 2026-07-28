"""Reconstruct a REAL cell-resolved microenvironment from the Xenium GBM section.

This is the flagship *real* reconstruction (see ``docs/REAL_BUILD_PLAN.md``): it
turns the real 40,887-cell 10x Xenium single-cell-resolution glioblastoma
section into a real, cell-resolved :class:`tmesim.core.Microenvironment`. There
is NO synthetic data anywhere in this path — every cell position and every
expression count is a measured byte on disk
(``tmesim.data.real.singlecell_spatial.load``).

What is real vs. provisional
----------------------------
REAL (measured):
  * per-cell (x, y) centroids in micrometers;
  * per-cell targeted-panel counts for the 319-gene Xenium Human Brain panel;
  * the marker-gene phenotype of each cell (derived only from that real signal).
PROVISIONAL (a documented modelling choice, NOT a real outcome):
  * ``ecm_density`` / ``vascular_density`` — Gaussian-smoothed spatial densities
    of the stromal / endothelial cell phenotypes (a proxy field, not a stain).

Per-cell drug sensitivity — now EXPLICIT, DRUG-SPECIFIC and deterministic
--------------------------------------------------------------------------
The reconstruction's ``CellPopulation.sensitivity`` is **uniform 1.0**: the
geometry is drug-agnostic, and "1.0" states plainly that no per-cell resistance
information is being asserted. It replaces an earlier placeholder that scaled
sensitivity by the mean TUMOR-marker (EGFR/SOX2/OLIG2/MKI67) z-score — a lineage
signal with no mechanistic link to chemosensitivity, which nonetheless silently
scaled every ABSOLUTE ``kill_fraction`` the model has ever reported.

Drug-specific sensitivity is applied at SIMULATION time by
:func:`per_cell_sensitivity`, which uses a measured determinant only where the
319-gene panel actually carries one (see :data:`SENSITIVITY_DETERMINANTS`), and
returns uniform 1.0 otherwise. Nothing is drawn from an RNG anywhere.

READ THIS BEFORE QUOTING AN ABSOLUTE KILL FRACTION: every absolute kill number is
CONDITIONAL on this assumption. With uniform sensitivity = 1.0 the predicted kill
is pure delivery + potency physics with no intrinsic-resistance heterogeneity —
i.e. an upper bound in which every tumour cell is fully drug-competent. With a
marker-derived sensitivity it is conditional on that marker being a valid
mechanistic proxy AND on the logistic mapping constants, neither of which is
calibrated against a measured drug response. This section is an untreated
resection with no outcome label; none of this is validated.

Phenotyping (marker -> CellType) — MATCHED-NULL GATED
-----------------------------------------------------
Markers are log1p-transformed then z-scored across all cells (so each marker
contributes on a comparable scale) and a type's raw score is the *mean* z-score
over its present markers. That raw mean-z is then converted, per cell, into a z
against a **matched random gene set of the same size and the same per-marker
detection rates** (:func:`_matched_null_moments`). A cell is assigned a type only
when that matched-null z both (a) clears :data:`TYPING_THRESHOLD` (>= 1.0) AND
(b) is the maximum among the gated types. Cells that clear no gate default to
``TUMOR`` (the malignant background of a GBM section).

WHY THE MATCHED NULL — the ARITY BUG it fixes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The previous rule gated the RAW mean-z at 1.0. Every marker z-score has unit
variance by construction, so a type's mean-z has standard deviation

    sd = sqrt( (1 + (n - 1) * rho) / n )

for ``n`` markers of mean pairwise correlation ``rho`` — i.e. it SHRINKS with the
number of markers. Measured on the three real sections (empirical sd matches that
formula to 3 decimal places):

    section              type          n   sd     % of cells clearing z >= 1
    brain 319   (s1)     ENDOTHELIAL   1   1.000  11.76
                         TUMOR         4   0.607   7.29
    IO 480      (s2)     ENDOTHELIAL   3   0.821   4.57
    GSE324268   (s3)     ENDOTHELIAL   1   1.000  14.03

So a FIXED absolute gate on a score whose spread depends on the marker count is
mechanically easier to clear for a one-marker type. ENDOTHELIAL is defined by one
marker (PECAM1) on the brain panel and by three (VWF/FLT1/PLVAP) on the IO panel,
so the bias differed per section and the reconstructed endothelial fraction swung
**9.74 / 4.18 / 12.19 %** (a 2.91x spread) across three sections whose true
vascularity is unlikely to differ by 3x. That was a CLASSIFIER DEFECT, not
biology, and not merely a panel artefact (s1 and s3 share the PECAM1 gate and
still differ 1.25x).

The matched null removes it at the source: the observed marker-set score is
compared against what a *random* gene set of the same size and the same detection
rates produces IN THE SAME CELL, so arity, per-marker sparsity and per-cell
sequencing depth are all controlled. The gate then means "this cell's marker-set
signal is >= 1 SD above a size- and abundance-matched random signature", which is
a statement of the same strength for a 1-marker and a 4-marker type. It is the
same control used successfully by the MGMT spatial work
(``validation/mgmt_spatial.py``, abundance-matched reference genes).

Measured effect (see ``results/classifier_fix.md`` for the full comparison):

    ENDOTHELIAL %        s1      s2      s3     max/min
    mean_z (old)        9.74    4.18   12.19     2.91x
    matched_null_z     10.50    7.07    8.75     1.49x

with the independent Visium NNLS deconvolution (``reconstruct/deconvolve.py``,
which uses no marker gate at all) putting ENDOTHELIAL at 9.77% on the parent
Visium sample. See :data:`TYPING_AGGREGATION` for the alternatives that were
evaluated and rejected.

HONEST RESIDUAL, opposite in sign. The matched null draws its ``n`` reference
genes INDEPENDENTLY, while a real marker set is a co-expressed programme (mean
pairwise marker correlation 0.16-0.48 on section 1). Correlated markers make the
observed score more variable than an independent null predicts, so the null SD is
slightly too small for a coherent multi-marker programme and its z is nudged UP —
which is part of why IMMUNE and STROMAL grow under the new rule. It is
second-order next to the ``sqrt(n)`` effect it replaces and it is NOT corrected
here: doing so needs a null reproducing the marker set's covariance, which
detection-rate matching alone cannot build. Treat a multi-marker class fraction
from this classifier as a mild over-estimate.

The map (only markers actually present in the panel are used):

  TUMOR       <- EGFR, SOX2, OLIG2, MKI67        (malignant glioma)
  IMMUNE      <- PTPRC, CD68, AIF1, CD163        (TAM / lymphoid)
  ENDOTHELIAL <- PECAM1                          (vessel-lining)
  STROMAL     <- MOG (oligodendrocyte),          (normal / stroma niche)
                 AQP4, GJA1 (astrocyte)

NECROTIC is a valid :class:`CellType` but is NOT assigned here: the panel
carries no necrosis marker, so no cell is phenotyped necrotic from expression
(reported honestly rather than invented).

Tissue mask — GEOMETRIC by default (a dx-independent physical domain)
--------------------------------------------------------------------
``Domain.tissue_mask`` used to be "grid cells containing at least one segmented
cell". That definition is a function of the GRID, not of the tissue: as ``dx``
shrinks, each cell occupies a smaller tile and the covered fraction collapses
(50.6% / 32.7% / 13.9% of the bounding box at dx = 30 / 20 / 12 um).
``real_fields._kde_on_grid`` multiplies the vascular field by that mask, so the
total O2 SUPPLY — and hence the mean oxygen and the hypoxic fraction — shrank as
the mesh was refined. That is a discretisation artifact masquerading as biology.

The default is now :func:`tissue_mask_from_cells` with ``mode='geometric'``: a
grid cell is tissue when its CENTRE lies within :data:`TISSUE_DISC_RADIUS_UM`
of at least one measured cell centroid (the union of fixed-radius physical
discs). The domain is then a property of the measured cell positions alone, so
refining the mesh only resolves the same physical region more finely. Measured
coverage is 67.2% / 67.4% / 67.7% at dx = 30 / 20 / 12 um (a 1.007x spread).
The legacy behaviour is still reachable via ``mode='occupancy'`` for comparison.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional

import numpy as np
from scipy import sparse
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from ..core.domain import Domain
from ..core.microenvironment import (
    CellPopulation,
    CellState,
    CellType,
    Microenvironment,
    TissueFields,
)
from ..data.real.singlecell_spatial import load as load_singlecell_spatial
from .marker_panels import BRAIN_PANEL_MARKER_MAP, select_marker_map

# ---------------------------------------------------------------------------
# Marker -> CellType map (documented; only present markers are used at runtime).
#
# UNCHANGED IN CONTENT: this is an ALIAS for the same dict, moved verbatim (same
# entries, same key order) into ``reconstruct.marker_panels`` so a SECOND panel
# can be registered alongside it. ``real_xenium.MARKER_MAP`` remains the 319-gene
# Xenium Human Brain map. Under the LEGACY ``aggregation='mean_z'`` rule the
# 319-panel reconstruction is still bit-for-bit its historical self (TUMOR 73.80 /
# STROMAL 7.92 / IMMUNE 8.54 / ENDOTHELIAL 9.74 %); the default aggregation has
# since changed to ``'matched_null_z'`` (TUMOR 58.24 / STROMAL 11.78 / IMMUNE
# 19.48 / ENDOTHELIAL 10.50 %) to remove the marker-count bias — see
# :data:`TYPING_AGGREGATION` and ``results/classifier_fix.md``.
#
# The map used at runtime is chosen from the PANEL that was actually loaded (see
# :func:`marker_panels.select_marker_map`), because this map is meaningless on
# the 480-gene Immuno-Oncology section, which carries none of PECAM1 / SOX2 /
# OLIG2 / AQP4 / GJA1 / MOG.
# ---------------------------------------------------------------------------
MARKER_MAP: Dict[CellType, List[str]] = BRAIN_PANEL_MARKER_MAP

# Bandwidth (um) for the Gaussian density fields (ECM / vascular). ~one cell
# diameter scale of spatial smoothing so a field varies over tissue niches.
_DENSITY_BANDWIDTH_UM = 50.0

# Positive marker-signal gate a type must clear to be assigned: the type's
# aggregated score must be >= 1.0. Under the default ``'matched_null_z'``
# aggregation that reads "the marker-set score is at least ~1 SD above what a
# size- and detection-rate-matched RANDOM gene set gives in this same cell".
# UNCHANGED at 1.0 across every aggregation evaluated in
# ``results/classifier_fix.md``, deliberately, so the comparison isolates the
# aggregation rule and nothing is tuned to the 8.0% volume-EM anchor.
TYPING_THRESHOLD = 1.0

#: How a type's per-cell score is aggregated from its marker z-scores. Four rules
#: were measured on all three real sections (``results/classifier_fix.md``);
#: ENDOTHELIAL % on brain-319 / IO-480 / GSE324268 and the resulting spread:
#:
#:   rule               s1     s2     s3    max/min  verdict
#:   'mean_z'          9.74   4.18  12.19    2.91x   LEGACY — the arity bug
#:   'max_z'           6.70   8.62   7.31    1.29x   rejected: arity-biased the
#:                                                   OTHER way (E[max of n] grows
#:                                                   with n), so its tightness
#:                                                   here is cancellation, not
#:                                                   control; also types every
#:                                                   cell off ONE extreme count
#:   'marker_vote'    10.37   8.37  13.77    1.65x   rejected: "fraction of
#:                                                   markers positive" saturates
#:                                                   at 1.0 for a 1-marker type,
#:                                                   which is arity bias again
#:   'mean_z_norm'     8.65   5.34  10.73    2.01x   partial: fixes arity but not
#:                                                   per-marker sparsity / depth
#:   'matched_null_z' 10.50   7.07   8.75    1.49x   CHOSEN
#:
#: ``'matched_null_z'`` is option (c): the per-cell z of the marker-set mean
#: against a same-size, detection-rate-matched random gene set, which controls
#: arity AND per-marker sparsity AND per-cell depth simultaneously. The other
#: rules remain reachable so the comparison stays reproducible and so the legacy
#: numbers can be regenerated exactly.
TYPING_AGGREGATION = "matched_null_z"

#: Number of detection-rate-matched panel genes forming the null pool for EACH
#: marker slot under ``'matched_null_z'``. The null moments are computed
#: ANALYTICALLY over the whole pool (uniform sampling with replacement), so there
#: is no RNG and no Monte-Carlo noise. 25 is a declared constant, not fitted:
#: ENDOTHELIAL % on the three sections moves 10.50/7.07/8.75 (K=25) ->
#: 10.26/7.65/8.43 (K=10) -> 10.40/6.88/8.67 (K=50), i.e. the choice is not
#: load-bearing (the 3-section spread is 1.49x / 1.34x / 1.51x respectively).
MATCHED_NULL_POOL_SIZE = 25

#: Floor applied to the matched-null standard deviation, as a fraction of its
#: median over cells. Prevents a degenerate (near-constant) matched pool from
#: manufacturing an unbounded z for a single cell.
_MATCHED_NULL_SIGMA_FLOOR_FRAC = 0.1

#: Every aggregation rule ``phenotype_cells`` accepts.
TYPING_AGGREGATIONS = ("matched_null_z", "mean_z", "mean_z_norm", "max_z",
                       "marker_vote")

# ---------------------------------------------------------------------------
# Tissue mask: a PHYSICAL domain, independent of the numerical grid.
# ---------------------------------------------------------------------------
#: Radius (um) of the disc placed around every measured cell centroid; their
#: union is the simulated tissue domain. A DECLARED MODELLING CONSTANT with a
#: physical, grid-free justification measured on this section:
#:
#:   * nearest-neighbour spacing between the 40,887 measured Xenium centroids is
#:     median **16.7 um** (p10 10.4, p90 24.2, p99 34.2 um);
#:   * 30 um is ~1.8x the median spacing and just above the p90 spacing, so the
#:     union of discs BRIDGES more than 90% of real inter-cell gaps and yields a
#:     connected tissue domain instead of a dotted one — which is what the
#:     diffusion operator needs, since an isolated tissue pixel cannot exchange
#:     drug or oxygen with anything;
#:   * it is also only ~1 cell diameter (a glioma soma is ~15-20 um) beyond the
#:     outermost measured centroid, so the tissue is not inflated into regions
#:     where nothing was actually segmented.
#:
#: It is NOT fitted to any outcome. Coverage is insensitive to the exact value
#: at the ~10 um level (0.571 / 0.677 / 0.710 of the bounding box at r = 20 / 30
#: / 40 um) but, crucially, is dx-INDEPENDENT at every radius.
TISSUE_DISC_RADIUS_UM = 30.0

#: Default tissue-mask construction. ``'geometric'`` = union of
#: :data:`TISSUE_DISC_RADIUS_UM` discs around measured cells (grid-independent).
#: ``'occupancy'`` = the legacy "grid cells containing >= 1 cell" rule, kept only
#: so the grid-dependence can be reproduced and compared against.
TISSUE_MASK_MODE = "geometric"


def tissue_mask_from_cells(
    x_um: np.ndarray,
    y_um: np.ndarray,
    nx: int,
    ny: int,
    dx_um: float,
    mode: str = TISSUE_MASK_MODE,
    radius_um: float = TISSUE_DISC_RADIUS_UM,
) -> np.ndarray:
    """Boolean (ny, nx) tissue mask from measured cell positions.

    Parameters
    ----------
    x_um, y_um : (n,) float
        Measured cell centroids in micrometers, already shifted so the minimum
        coordinate is 0 (the Domain frame).
    nx, ny, dx_um : grid definition.
    mode : {'geometric', 'occupancy'}
        ``'geometric'`` (default, GRID-INDEPENDENT) — a grid cell is tissue when
        its CENTRE ``((i+0.5)dx, (j+0.5)dx)`` lies within ``radius_um`` of at
        least one measured centroid, i.e. the mask is the union of fixed physical
        discs. The domain is then a property of the measured tissue, not of the
        mesh, so the covered PHYSICAL AREA is invariant to ``dx`` up to the
        O(dx) sampling of the disc boundary.

        ``'occupancy'`` (LEGACY, GRID-DEPENDENT) — a grid cell is tissue when it
        contains at least one measured centroid. Covered area collapses as ``dx``
        shrinks (each cell can only ever mark one tile), which made total O2
        supply — and therefore mean oxygen and the hypoxic fraction — a function
        of the mesh. Retained only for comparison; do not use for production.
    radius_um : float
        Disc radius for ``'geometric'`` (default :data:`TISSUE_DISC_RADIUS_UM`).

    Returns
    -------
    (ny, nx) bool ndarray.
    """
    x_um = np.asarray(x_um, dtype=np.float64)
    y_um = np.asarray(y_um, dtype=np.float64)
    nx, ny = int(nx), int(ny)
    dx_um = float(dx_um)
    m = str(mode).strip().lower()

    if m == "occupancy":
        mask = np.zeros((ny, nx), dtype=bool)
        if x_um.size:
            col = np.clip((x_um / dx_um).astype(int), 0, nx - 1)
            row = np.clip((y_um / dx_um).astype(int), 0, ny - 1)
            mask[row, col] = True
        return mask

    if m != "geometric":
        raise ValueError(
            f"unknown tissue-mask mode {mode!r}; expected 'geometric' or 'occupancy'"
        )

    if x_um.size == 0:
        return np.zeros((ny, nx), dtype=bool)

    # Distance from every grid-cell CENTRE to the nearest measured centroid.
    gx = (np.arange(nx, dtype=np.float64) + 0.5) * dx_um
    gy = (np.arange(ny, dtype=np.float64) + 0.5) * dx_um
    GX, GY = np.meshgrid(gx, gy)
    tree = cKDTree(np.column_stack([x_um, y_um]))
    dist, _ = tree.query(np.column_stack([GX.ravel(), GY.ravel()]), k=1)
    return (dist.reshape(ny, nx) <= float(radius_um))


def _zscore_columns(dense: np.ndarray) -> np.ndarray:
    """log1p then per-column (per-gene) z-score across cells. Returns float64."""
    x = np.log1p(dense.astype(np.float64))
    mean = x.mean(axis=0, keepdims=True)
    std = x.std(axis=0, keepdims=True)
    std[std == 0.0] = 1.0  # constant genes contribute 0, not NaN
    return (x - mean) / std


def _column_reader(expression):
    """Return ``get(j) -> (n_cells,) float64 counts`` for one gene column.

    Converts a sparse matrix to CSC ONCE so the matched-null pool (which touches
    :data:`MATCHED_NULL_POOL_SIZE` x n_markers columns) does not pay CSR
    column-slicing costs per read.
    """
    if sparse.issparse(expression):
        csc = expression.tocsc()

        def get(j: int) -> np.ndarray:
            return np.asarray(csc[:, int(j)].todense()).ravel().astype(np.float64)
    else:
        arr = np.asarray(expression)

        def get(j: int) -> np.ndarray:
            return arr[:, int(j)].astype(np.float64)

    return get


def _zscore_one(col: np.ndarray) -> np.ndarray:
    """log1p + z-score a single gene column (the 1-D form of _zscore_columns)."""
    return _zscore_columns(np.asarray(col, dtype=np.float64)[:, None]).ravel()


def detection_rates(expression) -> np.ndarray:
    """Per-gene fraction of cells with a non-zero count.

    The abundance statistic the matched null is matched ON: it is the dominant
    driver of how a targeted-panel marker's z-score is distributed (a gene seen in
    0.2% of cells has a near-degenerate, heavy-tailed z; one seen in 99% is close
    to continuous), so a null that ignores it is not a fair reference.
    """
    if sparse.issparse(expression):
        b = expression.tocsc(copy=True)
        b.data = np.ones_like(b.data)
        return np.asarray(b.sum(axis=0)).ravel() / float(expression.shape[0])
    arr = np.asarray(expression)
    return (arr > 0).mean(axis=0).astype(np.float64)


def _matched_pool(
    det: np.ndarray, target: float, exclude: set, pool_size: int
) -> List[int]:
    """Panel-gene columns whose detection rate is closest to ``target``.

    Distance is on the LOGIT of the detection rate, so "4.2% vs 4.4%" is treated
    as close while "90% vs 96%" is treated as the large difference it is in
    z-score shape terms. Genes that are markers of ANY type are excluded (a null
    must not contain the signal it is a null for), as are all-zero genes.
    """
    eps = 1e-6
    det = np.asarray(det, dtype=np.float64)
    d = np.clip(det, eps, 1.0 - eps)
    t = float(np.clip(target, eps, 1.0 - eps))
    dist = np.abs(np.log(d / (1.0 - d)) - np.log(t / (1.0 - t)))
    if exclude:
        dist[sorted(exclude)] = np.inf
    dist[det <= 0.0] = np.inf
    order = np.argsort(dist, kind="stable")
    return [int(j) for j in order[: int(pool_size)] if np.isfinite(dist[j])]


def _matched_null_moments(
    get,
    det: np.ndarray,
    marker_cols: List[int],
    exclude: set,
    pool_size: int,
    n_cells: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell mean and SD of the mean-z score of a MATCHED random gene set.

    The reference set has the same SIZE as ``marker_cols`` and, slot by slot, the
    same per-marker DETECTION RATE, so it is the right null for both the arity
    bug and per-marker sparsity. The moments are computed analytically over the
    whole pool rather than by Monte-Carlo draws — for independent uniform draws
    from each slot's pool,

        E[score]   = (1/n)   * sum_s  mean_{g in pool_s} z_g(cell)
        Var[score] = (1/n^2) * sum_s  var_{g in pool_s}  z_g(cell)

    — so the result is exact, deterministic and RNG-free.
    """
    n = len(marker_cols)
    mu = np.zeros(n_cells, dtype=np.float64)
    var = np.zeros(n_cells, dtype=np.float64)
    for c in marker_cols:
        pool = _matched_pool(det, float(det[c]), exclude, pool_size)
        if not pool:
            continue
        s1 = np.zeros(n_cells, dtype=np.float64)
        s2 = np.zeros(n_cells, dtype=np.float64)
        for j in pool:
            zz = _zscore_one(get(j))
            s1 += zz
            s2 += zz * zz
        k = float(len(pool))
        m1 = s1 / k
        mu += m1
        var += np.maximum(s2 / k - m1 * m1, 0.0)
    return mu / n, np.sqrt(var) / n


def phenotype_cells(
    expression,
    genes: np.ndarray,
    threshold: float = TYPING_THRESHOLD,
    marker_map: Optional[Dict[CellType, List[str]]] = None,
    aggregation: str = TYPING_AGGREGATION,
    pool_size: int = MATCHED_NULL_POOL_SIZE,
) -> tuple[np.ndarray, Dict[CellType, List[str]], List[str]]:
    """Assign a CellType to each cell from its marker z-scores, ARITY-CONTROLLED.

    Markers are log1p-transformed and z-scored across all cells; a type's raw
    score is the MEAN z over its present markers. That raw score is then turned
    into the gated score by ``aggregation``. A cell is assigned the type with the
    highest aggregated score AMONG the types clearing ``threshold``; cells that
    clear no gate default to ``TUMOR`` (the malignant GBM background).

    THE ARITY BUG THIS FIXES. A mean of ``n`` unit-variance z-scores has standard
    deviation ``sqrt((1 + (n-1)*rho)/n) <= 1``, so a FIXED gate on the raw mean-z
    is systematically easier to clear the FEWER markers a type has. ENDOTHELIAL
    has one marker on the brain panel (PECAM1) and three on the IO panel
    (VWF/FLT1/PLVAP), which is why the reconstructed endothelial fraction used to
    swing 9.74 / 4.18 / 12.19 % over three real sections. See the module docstring
    and ``results/classifier_fix.md``.

    Parameters
    ----------
    expression : scipy.sparse or ndarray (n_cells, n_genes)
        Real per-cell targeted-panel counts. ALL columns are used by
        ``'matched_null_z'`` (it needs the non-marker genes to build its null);
        the other rules touch only the marker columns.
    genes : ndarray of str
        Column labels for ``expression``.
    threshold : float
        Gate a type's aggregated score must clear (default
        :data:`TYPING_THRESHOLD` = 1.0). Held FIXED across aggregations so the
        comparison in ``results/classifier_fix.md`` isolates the aggregation rule
        and nothing is tuned to the 8.0% volume-EM anchor.
    aggregation : str
        One of :data:`TYPING_AGGREGATIONS`; default :data:`TYPING_AGGREGATION`.

        * ``'matched_null_z'`` (**chosen**) — per-cell z of the marker-set mean
          against a same-size, detection-rate-matched RANDOM gene set. Controls
          arity, per-marker sparsity and per-cell depth simultaneously.
        * ``'mean_z'`` — LEGACY raw mean-z. Reproduces the historical numbers
          exactly; carries the arity bug. Kept only for reproducibility.
        * ``'mean_z_norm'`` — mean-z divided by its own SD across cells, i.e.
          mean-z with the ``1/sqrt(n_markers)`` variance correction. Fixes arity
          but not per-marker sparsity or per-cell depth.
        * ``'max_z'`` — largest single-marker z. Arity-biased the OTHER way
          (``E[max of n]`` grows with ``n``) and types each cell off one count.
        * ``'marker_vote'`` — fraction of a type's markers individually
          z >= ``threshold``, tie-broken by max-z, requiring >= 1 positive
          marker. Also arity-biased (a 1-marker type saturates the score at 1.0).
    pool_size : int
        Matched-null pool size per marker slot (default
        :data:`MATCHED_NULL_POOL_SIZE`). Ignored by the other aggregations.
    marker_map : dict CellType -> list of str, optional
        Marker map to type with. ``None`` (the default) resolves it from the
        PANEL via :func:`marker_panels.select_marker_map`, which returns
        :data:`MARKER_MAP` itself for the 319-gene brain panel and the
        Immuno-Oncology map for the 480-gene section, whose panel carries none of
        PECAM1 / SOX2 / OLIG2 / AQP4 / GJA1 / MOG and on which the brain map
        would collapse the whole reconstruction to TUMOR + IMMUNE with no
        vasculature (hence no oxygen source).

    Returns
    -------
    cell_type : ndarray (n_cells,) int8 of CellType values.
    used_map : dict CellType -> list of markers actually found in the panel.
    missing : list of markers requested but absent from the panel.
    """
    scores, types, used_map, missing = marker_scores(
        expression, genes, threshold=threshold, marker_map=marker_map,
        aggregation=aggregation, pool_size=pool_size,
    )

    # Threshold gate: a type is only assignable to a cell when its aggregated
    # score clears the positive threshold. Non-clearing scores are set to -inf so
    # they can never win the arg-max.
    gated = np.where(scores >= float(threshold), scores, -np.inf)
    winner = gated.argmax(axis=1)
    cell_type = np.array([int(types[i]) for i in winner], dtype=np.int8)

    # Cells that clear NO gate (all types below threshold) default to the
    # malignant background: TUMOR.
    ungated = ~np.isfinite(gated.max(axis=1))
    cell_type[ungated] = int(CellType.TUMOR)
    return cell_type, used_map, missing


def marker_scores(
    expression,
    genes: np.ndarray,
    threshold: float = TYPING_THRESHOLD,
    marker_map: Optional[Dict[CellType, List[str]]] = None,
    aggregation: str = TYPING_AGGREGATION,
    pool_size: int = MATCHED_NULL_POOL_SIZE,
) -> tuple[np.ndarray, List[CellType], Dict[CellType, List[str]], List[str]]:
    """The per-cell, per-type SCORES that :func:`phenotype_cells` gates.

    Split out so the arity audit (``validation/classifier_fix.py``) and the tests
    can inspect the score distributions directly instead of re-deriving them, and
    so the two can never drift apart. Arguments are exactly those of
    :func:`phenotype_cells`; ``threshold`` is only used by the ``'marker_vote'``
    rule, which needs it to make its per-marker positive calls.

    Returns
    -------
    scores : (n_cells, n_types) float64. ``-inf`` for a type with no measured
        marker on this panel.
    types : list of CellType, aligned to the score columns (map key order, which
        also sets the arg-max tie-break).
    used_map, missing : as :func:`phenotype_cells`.
    """
    agg = str(aggregation).strip().lower()
    if agg not in TYPING_AGGREGATIONS:
        raise ValueError(
            f"unknown aggregation {aggregation!r}; expected one of "
            f"{TYPING_AGGREGATIONS}"
        )

    genes = np.asarray(genes, dtype=str)
    gene_index = {g: j for j, g in enumerate(genes)}

    if marker_map is None:
        marker_map, _panel_name, _cov = select_marker_map(genes)

    used_map: Dict[CellType, List[str]] = {}
    missing: List[str] = []
    needed_cols: List[int] = []
    for ct, markers in marker_map.items():
        present = [m for m in markers if m in gene_index]
        absent = [m for m in markers if m not in gene_index]
        missing.extend(absent)
        used_map[ct] = present
        needed_cols.extend(gene_index[m] for m in present)

    # Densify only the marker columns (small; keeps memory light for 40k cells).
    needed_cols = sorted(set(needed_cols))
    if sparse.issparse(expression):
        sub = expression[:, needed_cols].toarray()
    else:
        sub = np.asarray(expression)[:, needed_cols]
    z = _zscore_columns(sub)  # (n_cells, len(needed_cols))
    col_pos = {c: k for k, c in enumerate(needed_cols)}

    n_cells = z.shape[0]
    types = list(used_map.keys())
    n_types = len(types)
    type_cols = [[col_pos[gene_index[m]] for m in used_map[ct]] for ct in types]

    mean_z = np.full((n_cells, n_types), -np.inf, dtype=np.float64)
    for t_idx, cols in enumerate(type_cols):
        if cols:
            mean_z[:, t_idx] = z[:, cols].mean(axis=1)

    # --- aggregate the raw mean-z into the score the gate is applied to ------
    if agg == "mean_z":
        scores = mean_z

    elif agg == "max_z":
        scores = np.full((n_cells, n_types), -np.inf, dtype=np.float64)
        for t_idx, cols in enumerate(type_cols):
            if cols:
                scores[:, t_idx] = z[:, cols].max(axis=1)

    elif agg == "mean_z_norm":
        # Every marker z has mean 0, so mean_z has mean 0 too; dividing by its own
        # SD puts a type with ANY number of markers on the same scale, which is
        # exactly the n_markers variance correction.
        scores = np.array(mean_z, copy=True)
        for t_idx, cols in enumerate(type_cols):
            if cols:
                s = scores[:, t_idx].std()
                if s > 0:
                    scores[:, t_idx] = scores[:, t_idx] / s

    elif agg == "marker_vote":
        # Fraction of a type's markers individually positive; max-z breaks ties;
        # a type needs at least one positive marker to be assignable at all.
        frac_pos = np.zeros((n_cells, n_types), dtype=np.float64)
        tie = np.zeros((n_cells, n_types), dtype=np.float64)
        for t_idx, cols in enumerate(type_cols):
            if not cols:
                continue
            frac_pos[:, t_idx] = (z[:, cols] >= float(threshold)).mean(axis=1)
            mx = z[:, cols].max(axis=1)
            rng = float(np.ptp(mx))
            tie[:, t_idx] = (mx - mx.min()) / (rng if rng > 0 else 1.0)
        # Encoded so the shared gate below fires exactly when a type has >= 1
        # positive marker, and the arg-max then prefers the higher vote fraction.
        scores = np.where(
            frac_pos > 0.0, float(threshold) + frac_pos + 1e-6 * tie, -np.inf
        )

    else:  # 'matched_null_z' — the default
        det = detection_rates(expression)
        get = _column_reader(expression)
        exclude = set(needed_cols)
        scores = np.full((n_cells, n_types), -np.inf, dtype=np.float64)
        for t_idx, ct in enumerate(types):
            cols = [gene_index[m] for m in used_map[ct]]
            if not cols:
                continue
            mu, sigma = _matched_null_moments(
                get, det, cols, exclude, int(pool_size), n_cells
            )
            pos = sigma[sigma > 0]
            floor = (
                _MATCHED_NULL_SIGMA_FLOOR_FRAC * float(np.median(pos))
                if pos.size else 1e-6
            )
            scores[:, t_idx] = (mean_z[:, t_idx] - mu) / np.maximum(sigma, floor)

    return scores, types, used_map, sorted(set(missing))


# ---------------------------------------------------------------------------
# Per-cell drug sensitivity: measured determinant where one exists, else 1.0.
# ---------------------------------------------------------------------------
#: Sensitivity floor. A marker-derived sensitivity is mapped into
#: ``[SENSITIVITY_FLOOR, 1.0]`` so a marker-negative cell is *less* sensitive but
#: never zero: "no transcript detected" is a dropout-prone statement about a
#: targeted panel, not a statement that the cell is immortal. A hard 0 would let
#: technical dropout fabricate perfectly resistant cells.
SENSITIVITY_FLOOR = 0.2

#: Logistic steepness applied to the marker z-score. 1.0 means a cell one SD above
#: the cohort mean sits at sigmoid(1) = 0.73 of the way up the band. Declared
#: modelling constant, not fitted.
SENSITIVITY_SLOPE = 1.0

#: Drug -> the measured panel gene that determines its per-cell sensitivity.
#:
#: A drug appears here ONLY when the 319-gene Xenium Human Brain panel carries a
#: gene that is a genuine mechanistic determinant of that drug's cell kill. A drug
#: mapped to ``None`` (or absent from this table) gets uniform 1.0 — no invented
#: heterogeneity.
#:
#: Panel audit (what is and is not measurable here):
#:   PRESENT: MKI67, TOP2A, EGFR, TP53, AQP4, PCNA, CDK1
#:   ABSENT : MGMT, ABCB1, ABCG2, CDKN2A, ATM, MSH2, MLH1, PARP1
#:
#: Per-drug rationale:
#:   * doxorubicin / etoposide -> TOP2A. Topoisomerase II alpha is the literal
#:     molecular TARGET of these drugs: they are TOP2 poisons that trap the
#:     enzyme-DNA cleavage complex, so cytotoxicity scales with how much TOP2A
#:     protein a cell actually has. TOP2A-low cells are correspondingly resistant.
#:     (Nitiss, Nat Rev Cancer 2009;9(5):338-50, PMID 19377506; Burgess et al.,
#:     PNAS 2008;105(26):9053-8, PMID 18574145 -- TOP2A level determines
#:     sensitivity to TOP2 poisons.) This is real target biology, not a proxy.
#:   * cisplatin / lomustine / carmustine / methotrexate / vincristine -> MKI67.
#:     These are cell-cycle-dependent cytotoxics (DNA cross-linkers, an antifolate
#:     that blocks S-phase dTMP synthesis, and a mitotic spindle poison); cycling
#:     cells are markedly more sensitive than quiescent ones. Ki-67 is the standard
#:     marker of a cell being in cycle. This is a MECHANISTIC PROXY (proliferation-
#:     gated cytotoxicity), NOT a measured sensitivity.
#:   * temozolomide -> None. TMZ's canonical resistance determinant is MGMT
#:     (O6-methylguanine-DNA methyltransferase) status, which repairs the O6-meG
#:     lesion TMZ depends on (Hegi et al., NEJM 2005;352:997-1003, PMID 15758010).
#:     MGMT is NOT on this panel and methylation is not measurable from
#:     transcript counts anyway, so TMZ per-cell sensitivity is UNMODELLED: uniform
#:     1.0, stated explicitly rather than faked from an unrelated marker.
#:   * bevacizumab -> None. An anti-VEGF antibody; it does not kill tumour cells
#:     directly, so a per-tumour-cell cytotoxic sensitivity is not meaningful.
SENSITIVITY_DETERMINANTS: Dict[str, Optional[Dict[str, str]]] = {
    "doxorubicin": {
        "gene": "TOP2A",
        "kind": "molecular target",
        "rationale": "topoisomerase II alpha is doxorubicin's direct target; "
                     "TOP2A level sets sensitivity to TOP2 poisons",
    },
    "etoposide": {
        "gene": "TOP2A",
        "kind": "molecular target",
        "rationale": "etoposide is a TOP2 poison; cytotoxicity scales with TOP2A",
    },
    "mitoxantrone": {
        "gene": "TOP2A",
        "kind": "molecular target",
        "rationale": "anthracenedione TOP2 poison; same target as doxorubicin",
    },
    "cisplatin": {
        "gene": "MKI67",
        "kind": "mechanistic proxy (proliferation-gated cytotoxicity)",
        "rationale": "DNA cross-linker; cycling cells are far more sensitive",
    },
    "lomustine": {
        "gene": "MKI67",
        "kind": "mechanistic proxy (proliferation-gated cytotoxicity)",
        "rationale": "nitrosourea DNA cross-linker; cell-cycle dependent",
    },
    "carmustine": {
        "gene": "MKI67",
        "kind": "mechanistic proxy (proliferation-gated cytotoxicity)",
        "rationale": "nitrosourea DNA cross-linker; cell-cycle dependent",
    },
    "methotrexate": {
        "gene": "MKI67",
        "kind": "mechanistic proxy (proliferation-gated cytotoxicity)",
        "rationale": "antifolate blocking S-phase dTMP synthesis; S-phase specific",
    },
    "vincristine": {
        "gene": "MKI67",
        "kind": "mechanistic proxy (proliferation-gated cytotoxicity)",
        "rationale": "mitotic spindle poison; only kills cells entering mitosis",
    },
    "temozolomide": None,   # determinant is MGMT -- NOT on this panel
    "procarbazine": None,   # determinant is MGMT/mismatch repair -- NOT on panel
    "bevacizumab": None,    # anti-VEGF antibody; no direct tumour-cell cytotoxicity
}


def _marker_zscore(expression, genes: np.ndarray, marker: str) -> Optional[np.ndarray]:
    """Per-cell log1p z-score of one panel gene, or None if it is not on the panel."""
    genes = np.asarray(genes, dtype=str)
    gene_index = {g: j for j, g in enumerate(genes)}
    if marker not in gene_index:
        return None
    col = gene_index[marker]
    if sparse.issparse(expression):
        sub = expression[:, [col]].toarray()
    else:
        sub = np.asarray(expression)[:, [col]]
    return _zscore_columns(sub).ravel()


def per_cell_sensitivity(
    expression, genes: np.ndarray, drug_name: str
) -> tuple[np.ndarray, Dict]:
    """DRUG-SPECIFIC per-cell sensitivity in ``[SENSITIVITY_FLOOR, 1.0]``.

    Deterministic (no RNG). Looks the drug up in :data:`SENSITIVITY_DETERMINANTS`:

    * **A determinant gene on the panel** -> sensitivity is the per-cell log1p
      z-score of that gene squashed through a logistic into the bounded band::

          s_c = FLOOR + (1 - FLOOR) * sigmoid(SLOPE * z_c)

      Bounded so the marker can never fabricate a perfectly resistant (0) cell out
      of transcript dropout, and monotone so a marker-high cell is always the more
      sensitive one. For ``TOP2A`` (doxorubicin/etoposide/mitoxantrone) this is the
      drug's actual molecular target; for ``MKI67`` (cisplatin, nitrosoureas,
      methotrexate, vincristine) it is a proliferation-gated *mechanistic proxy*,
      not a measured sensitivity.

    * **No determinant on the panel** (temozolomide -> MGMT; bevacizumab; anything
      unknown) -> uniform **1.0**. No per-cell information, so none is invented,
      and the prediction is pure delivery + potency physics.

    HONEST LIMITATION -- read this with the returned ``info``. On this targeted
    319-gene panel the determinant transcripts are extremely sparse: TOP2A is
    detected in ~1.1% of cells and MKI67 in ~0.35%. The z-score of a mostly-zero
    count vector is therefore near-binary: ~99% of cells sit at the single value
    ``sigmoid(z_zero)`` and a small positive tail sits near 1.0. So the marker-derived
    field behaves mostly as a near-uniform SCALE on kill (set by the arbitrary
    ``sigmoid(0) = 0.5`` midpoint of the mapping) plus a thin genuinely-high tail —
    it is NOT a rich per-cell heterogeneity map, and the absolute kill it produces
    is conditional on the mapping constants. ``info`` reports the detection rate so
    this is auditable at every call site.

    Parameters
    ----------
    expression : scipy.sparse or ndarray (n_cells, n_genes) -- real panel counts.
    genes : ndarray of str -- column labels for ``expression``.
    drug_name : str -- case-insensitive drug name.

    Returns
    -------
    sensitivity : (n_cells,) float64 in ``[SENSITIVITY_FLOOR, 1.0]``.
    info : dict describing exactly what was used (determinant gene, kind,
        rationale, detection rate, resulting distribution) for reporting.
    """
    n_cells = int(expression.shape[0])
    name = str(drug_name).strip().lower()
    spec = SENSITIVITY_DETERMINANTS.get(name)

    def _uniform(reason: str) -> tuple[np.ndarray, Dict]:
        sens = np.ones(n_cells, dtype=np.float64)
        return sens, {
            "drug": name,
            "determinant_gene": None,
            "kind": "none (uniform 1.0)",
            "rationale": reason,
            "marker_detection_rate": None,
            "mean": 1.0, "std": 0.0, "min": 1.0, "max": 1.0,
            "is_uniform": True,
        }

    if spec is None:
        if name in SENSITIVITY_DETERMINANTS:
            return _uniform(
                "no mechanistic determinant for this drug is measurable on the "
                "319-gene panel (e.g. temozolomide depends on MGMT, which is absent); "
                "sensitivity is UNMODELLED and set to uniform 1.0"
            )
        return _uniform(
            f"{name!r} has no entry in SENSITIVITY_DETERMINANTS; no per-cell "
            "sensitivity information is asserted"
        )

    z = _marker_zscore(expression, genes, spec["gene"])
    if z is None:
        return _uniform(
            f"determinant {spec['gene']} is not on this panel; sensitivity is "
            "UNMODELLED and set to uniform 1.0"
        )

    sens = SENSITIVITY_FLOOR + (1.0 - SENSITIVITY_FLOOR) / (
        1.0 + np.exp(-SENSITIVITY_SLOPE * z)
    )
    sens = np.clip(sens, SENSITIVITY_FLOOR, 1.0)

    # Detection rate of the determinant transcript (the honesty statistic).
    genes_arr = np.asarray(genes, dtype=str)
    col = int(np.flatnonzero(genes_arr == spec["gene"])[0])
    if sparse.issparse(expression):
        counts = np.asarray(expression[:, [col]].todense()).ravel()
    else:
        counts = np.asarray(expression)[:, col]
    detection = float((counts > 0).mean())

    return sens, {
        "drug": name,
        "determinant_gene": spec["gene"],
        "kind": spec["kind"],
        "rationale": spec["rationale"],
        "marker_detection_rate": detection,
        "mean": float(sens.mean()),
        "std": float(sens.std()),
        "min": float(sens.min()),
        "max": float(sens.max()),
        "is_uniform": False,
    }


def _density_field(
    x_um: np.ndarray,
    y_um: np.ndarray,
    domain: Domain,
    bandwidth_um: float = _DENSITY_BANDWIDTH_UM,
) -> np.ndarray:
    """Gaussian kernel-density estimate of point (x,y) on the grid, in [0,1].

    Points are binned onto the (ny, nx) grid then convolved with a Gaussian of
    the given bandwidth (a discretised Gaussian KDE), and rescaled by its max.
    """
    grid = domain.zeros()
    if len(x_um) > 0:
        col, row = domain.to_index(x_um, y_um)
        np.add.at(grid, (row, col), 1.0)
    sigma = max(bandwidth_um / domain.dx, 1e-6)
    dens = gaussian_filter(grid, sigma=sigma, mode="constant")
    m = dens.max()
    if m > 0:
        dens = dens / m
    return dens


def reconstruct_from_xenium(
    sample_path: str | None = None,
    dx_um: float = 30.0,
    tissue_mask_mode: str = TISSUE_MASK_MODE,
    tissue_disc_radius_um: float = TISSUE_DISC_RADIUS_UM,
    typing_aggregation: str = TYPING_AGGREGATION,
) -> Microenvironment:
    """Reconstruct the real cell-resolved GBM microenvironment from Xenium.

    Parameters
    ----------
    sample_path : str, optional
        Directory holding the cached Xenium files. ``None`` uses the loader
        default (``data_cache/xenium_gbm``).
    dx_um : float
        Grid spacing in micrometers for the continuum Domain / fields.
    tissue_mask_mode : {'geometric', 'occupancy'}
        How ``Domain.tissue_mask`` is built (see :func:`tissue_mask_from_cells`).
        Default ``'geometric'`` — the union of fixed-radius physical discs around
        the measured centroids, which makes the simulated domain independent of
        ``dx_um``. ``'occupancy'`` reproduces the legacy grid-dependent rule.
    tissue_disc_radius_um : float
        Disc radius for the geometric mask (default
        :data:`TISSUE_DISC_RADIUS_UM` = 30 um; see that constant for the
        measured justification).
    typing_aggregation : str
        Marker-score aggregation passed to :func:`phenotype_cells` (default
        :data:`TYPING_AGGREGATION` = ``'matched_null_z'``). Pass ``'mean_z'`` to
        regenerate the legacy, arity-biased composition exactly.

    Returns
    -------
    Microenvironment
        Real cells at real (x, y) positions with marker-derived CellTypes, a
        Domain covering the real tissue extent, and minimal continuum fields
        (vascular from endothelial density, ECM from stromal density; oxygen
        left None for the physics module to fill).
    """
    # --- load real data (retry once on transient failure) -------------------
    kwargs = {"dense": False}
    if sample_path is not None:
        kwargs["path"] = sample_path
    try:
        data = load_singlecell_spatial(**kwargs)
    except Exception:
        data = load_singlecell_spatial(**kwargs)

    xy = np.asarray(data["cell_xy_um"], dtype=np.float64)
    expression = data["expression"]
    genes = np.asarray(data["genes"], dtype=str)

    # --- phenotype every cell from real marker signal -----------------------
    # The marker map follows the PANEL that was actually loaded: the 319-gene
    # brain panel resolves to the unchanged MARKER_MAP, the 480-gene
    # Immuno-Oncology section to its own map (see reconstruct.marker_panels).
    panel_map, panel_name, panel_coverage = select_marker_map(genes)
    cell_type, used_map, missing = phenotype_cells(
        expression, genes, marker_map=panel_map, aggregation=typing_aggregation
    )

    # --- per-cell sensitivity: UNIFORM 1.0 (drug-agnostic geometry) ----------
    # Deliberately NOT random and NOT marker-scaled. 1.0 asserts "no per-cell
    # intrinsic-resistance information", so a kill fraction computed on this
    # reconstruction is pure delivery + potency physics. Drug-specific
    # sensitivity is layered on at simulation time via per_cell_sensitivity()
    # (see real_pipeline.apply_drug_sensitivity). ABSOLUTE kill fractions are
    # CONDITIONAL on whichever of the two is in force.
    sensitivity = np.ones(len(xy), dtype=np.float64)

    # --- Domain covering the real extent (shift so min coord -> 0) ----------
    x = xy[:, 0] - xy[:, 0].min()
    y = xy[:, 1] - xy[:, 1].min()
    nx = int(np.ceil(x.max() / dx_um)) + 1
    ny = int(np.ceil(y.max() / dx_um)) + 1

    # tissue_mask: GRID-INDEPENDENT by default (union of fixed-radius physical
    # discs around the measured centroids). See tissue_mask_from_cells.
    tissue_mask = tissue_mask_from_cells(
        x, y, nx, ny, dx_um,
        mode=tissue_mask_mode, radius_um=tissue_disc_radius_um,
    )

    domain = Domain(nx=nx, ny=ny, dx=dx_um, tissue_mask=tissue_mask)

    # --- CellPopulation: real positions + types, all VIABLE -----------------
    state = np.full(len(x), int(CellState.VIABLE), dtype=np.int8)
    cells = CellPopulation(
        x=x, y=y, cell_type=cell_type, state=state, sensitivity=sensitivity
    )

    # --- minimal continuum fields from real cell densities ------------------
    endo = cell_type == int(CellType.ENDOTHELIAL)
    strom = cell_type == int(CellType.STROMAL)
    vascular_density = _density_field(x[endo], y[endo], domain)
    ecm_density = _density_field(x[strom], y[strom], domain)
    fields = TissueFields(
        ecm_density=ecm_density,
        vascular_density=vascular_density,
        oxygen=None,  # filled later by the physics module
    )

    # --- type distribution for metadata / reporting -------------------------
    dist = {ct.name: int((cell_type == int(ct)).sum()) for ct in CellType}
    metadata = {
        "assay": data["source"]["dataset"],
        "panel": data["source"]["panel"],
        "n_panel_genes": int(len(genes)),
        "marker_map": {ct.name: used_map.get(ct, []) for ct in panel_map},
        "marker_panel_used": panel_name,
        "marker_panel_coverage": panel_coverage,
        "markers_missing_from_panel": missing,
        "typing_aggregation": str(typing_aggregation),
        "typing_threshold": float(TYPING_THRESHOLD),
        "typing_note": (
            "Cell types come from phenotype_cells with aggregation="
            f"{typing_aggregation!r}. The default 'matched_null_z' scores each "
            "type by the per-cell z of its marker-set mean against a same-size, "
            "detection-rate-matched RANDOM gene set, which removes the "
            "marker-COUNT (arity) bias of the legacy raw mean-z rule: a mean of n "
            "unit-variance z-scores has SD sqrt((1+(n-1)rho)/n), so a fixed gate "
            "was mechanically easier to clear for a 1-marker type (ENDOTHELIAL on "
            "the brain panel) than for a 3- or 4-marker one. Endothelial fraction "
            "over three real sections went 9.74/4.18/12.19 % (2.91x spread) -> "
            "10.50/7.07/8.75 % (1.49x). See results/classifier_fix.md."
        ),
        "celltype_distribution": dist,
        "tissue_mask_mode": str(tissue_mask_mode),
        "tissue_disc_radius_um": (
            float(tissue_disc_radius_um)
            if str(tissue_mask_mode).lower() == "geometric" else None
        ),
        "tissue_coverage_fraction": float(tissue_mask.mean()),
        "tissue_mask_note": (
            "GEOMETRIC: tissue = union of "
            f"{tissue_disc_radius_um:.0f} um discs around the measured cell "
            "centroids, so the simulated domain is a property of the measured "
            "tissue and NOT of the grid spacing (coverage 67.2/67.4/67.7% at "
            "dx = 30/20/12 um). Replaces the legacy 'grid cells containing >=1 "
            "cell' rule, whose coverage collapsed 50.6/32.7/13.9% over the same "
            "grids and made total O2 supply grid-dependent."
            if str(tissue_mask_mode).lower() == "geometric" else
            "LEGACY OCCUPANCY mask (grid cells containing >=1 measured cell). "
            "GRID-DEPENDENT — retained only for comparison."
        ),
        "sensitivity_model": "uniform_1.0",
        "sensitivity_note": (
            "UNIFORM 1.0 — no per-cell intrinsic-resistance information is asserted, "
            "so kill is pure delivery + potency physics. (Replaces an earlier "
            "placeholder scaled from the TUMOR-marker z-score, which had no "
            "mechanistic link to chemosensitivity yet scaled every absolute kill "
            "fraction.) Drug-specific sensitivity is applied at simulation time by "
            "reconstruct.real_xenium.per_cell_sensitivity; see "
            "SENSITIVITY_DETERMINANTS. ABSOLUTE kill fractions are CONDITIONAL on "
            "this assumption. This section is an untreated resection with no "
            "measured drug response."
        ),
        "fields_note": (
            "vascular_density = Gaussian KDE of ENDOTHELIAL cells; "
            "ecm_density = Gaussian KDE of STROMAL cells; oxygen left None."
        ),
        "has_outcome_label": bool(data["source"].get("has_outcome_label", False)),
    }

    return Microenvironment(
        domain=domain,
        fields=fields,
        cells=cells,
        patient_id="xenium-gbm-real",
        metadata=metadata,
    )


if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    micro = reconstruct_from_xenium()

    print("=== REAL Xenium GBM microenvironment ===")
    print("summary:")
    for k, v in micro.summary().items():
        print(f"  {k}: {v}")

    print(f"\ngrid dims (ny, nx): {micro.domain.shape}  dx={micro.domain.dx} um")
    print(f"extent (um):        {micro.domain.extent_um}")

    print("\nmarker -> type map (present markers):")
    for ct_name, markers in micro.metadata["marker_map"].items():
        print(f"  {ct_name:<12}: {markers}")
    if micro.metadata["markers_missing_from_panel"]:
        print("markers requested but ABSENT from panel:",
              micro.metadata["markers_missing_from_panel"])

    ct = micro.cells.cell_type
    n = micro.cells.n
    print(f"\ncell-type distribution among {n:,} real cells:")
    for c in CellType:
        k = int((ct == int(c)).sum())
        print(f"  {c.name:<12}: {k:>7,}  ({100.0 * k / n:5.2f}%)")

    # --- scatter colored by CellType --------------------------------------
    colors = {
        CellType.TUMOR: "#d62728",
        CellType.STROMAL: "#2ca02c",
        CellType.IMMUNE: "#1f77b4",
        CellType.ENDOTHELIAL: "#ff7f0e",
        CellType.NECROTIC: "#7f7f7f",
    }
    fig, ax = plt.subplots(figsize=(10, 10))
    for c in CellType:
        m = ct == int(c)
        if m.any():
            ax.scatter(
                micro.cells.x[m], micro.cells.y[m], s=2, c=colors[c],
                label=f"{c.name} ({int(m.sum()):,})", linewidths=0,
            )
    ax.set_aspect("equal")
    ax.invert_yaxis()  # image convention
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title("Real Xenium GBM microenvironment — cells by phenotype")
    ax.legend(markerscale=4, loc="upper right", fontsize=8)

    out_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "results"
    )
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "real_xenium_microenvironment.png")
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    print(f"\nfigure written: {out_path}")
