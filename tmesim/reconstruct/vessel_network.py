"""DISCRETE vessel objects from REAL endothelial cells — the structural fix.

Why this module exists
----------------------
Every drug in this project is delivered by ``source = vascular_density *
bbb_penetrance * plasma`` (``treatment.drugs.delivery_source``). Until now
``vascular_density`` has been a **40 µm-bandwidth Gaussian KDE of measured
endothelial-cell positions, divided by its own maximum**
(``reconstruct.real_fields.vascular_density_from_cells``). Four independent
analyses in ``results/`` have converged on the same defect:

* ``results/delivery_metrics.md`` — the fitted penetration length is 26.5–28.2 µm
  for **every** drug in the five-drug panel against a **source-null** of 23.7 µm.
  The number being measured is the KDE bandwidth, not drug transport.
* ``results/dose_schedule.md`` — the exposure-vs-distance profile is a steep
  perivascular drop sitting on a **flat floor at ~12% of the perivascular value**
  that persists to the far edge of the section, because a continuous source field
  keeps feeding tissue a millimetre from any vessel. Hence "delivery is never the
  binding constraint".
* ``results/niche_response.md`` — "perivascular tissue gets more drug" is
  circular: the region label and the source field are both functions of the same
  endothelial marker signal.
* ``results/vascular_validation.md`` §7.4 — "a Gaussian KDE with a 40 µm bandwidth
  applied to points that are ~44 µm apart … the smoothing kernel is the same size
  as the structure it is meant to resolve, so it necessarily erases it."

Real tumour vasculature is **discrete**: distinct vessel objects with genuinely
avascular tissue in between. This module builds that representation.

What it does
------------
1. :func:`segment_vessels` — single-linkage / DBSCAN clustering of the measured
   ``CellType.ENDOTHELIAL`` positions at a documented physical linking distance,
   yielding discrete putative vessel objects with a size distribution and an
   inter-vessel distance distribution that can be compared to published human GBM
   microvascular morphometry.
2. :func:`discrete_vascular_field` — a source field that is nonzero **only** within
   a documented vessel radius of a segmented vessel and **exactly zero** elsewhere.
   No Gaussian tail, no floor.

Both operate on the CORE interface only (``Domain`` + ``CellPopulation``), exactly
like ``reconstruct.real_fields``. Nothing in ``core`` / ``physics`` /
``sim.engine`` / ``molecular`` is touched.

Parameter provenance (both constants are load-bearing; both are cited)
----------------------------------------------------------------------
``LINK_DISTANCE_UM = 20.0`` — the distance below which two endothelial cells are
taken to belong to the same vessel object.

* ``results/vascular_validation.md`` §3.2 measured, on this very section, a median
  nearest-other-endothelial-cell distance of **26.4 µm**, and §3.3 measured a
  radial same-type enrichment of **2.8× within 10 µm and 2.8× at 10–20 µm**,
  decaying to the all-cell null by ~60 µm — the short-range signature of cells
  sharing a vessel wall. 20 µm sits inside that enriched shell and outside the
  regime where the enrichment has decayed away. That report also adopts 20 µm as
  its reference linking distance, so this module stays consistent with it.
* Physically, endothelial cells are elongated along the vessel axis with a
  cell length of order 20–30 µm, so consecutive nuclei in one vessel wall are
  ~20 µm apart.
* It is a free parameter and it is **swept**, not asserted: pass
  ``link_distances`` to :func:`link_distance_sweep`.

``VESSEL_RADIUS_UM = 7.5`` — the radius around a vessel's endothelial cells inside
which the source is nonzero. Interpreted as the **outer wall radius of the
perfused vessel profile**, i.e. the annulus of tissue in immediate contact with
leaking endothelium, not the lumen alone.

* In-vivo human vessel-size imaging gives a mean vessel **caliber (diameter)** of
  37.4 ± 12.9 µm in cortical grey matter, 20.7 ± 8.8 µm in thalamus and
  15.0 ± 5.1 µm in white matter, with **all brain tumours larger than
  normal-appearing tissue** — Hsu YY et al., J Magn Reson Imaging 2009;30(5):
  1078-84, PMID 19856441, doi:10.1002/jmri.21960. That is a radius of ~7.5–19 µm,
  and 7.5 µm is the **low** end of it.
* In an orthotopic U87 glioma, MRI vessel size validated against
  immunohistochemistry was 18.1 ± 2.4 µm diameter untreated (12.8 ± 2.0 µm after
  bevacizumab) — Kording F et al., J Magn Reson Imaging 2014;40(6):1310-8,
  PMID 24390982, doi:10.1002/jmri.24513. Radius ~9 µm. (Mouse xenograft.)
* The normal-brain lower bound is a mean capillary radius of 2.39 ± 0.08 µm —
  Ji X et al., Neuron 2021;109(7):1168-1187, PMID 33657412,
  doi:10.1016/j.neuron.2021.02.006. (Mouse, light-sheet whole-brain.)
* **SUPERSEDED — see** :mod:`reconstruct.vessel_geometry` **and**
  ``results/vessel_geometry_fix.md``. The paragraph that used to stand here said
  that at ``dx = 12 µm`` a grid cell has an equivalent disc radius of 6.8 µm so
  "no vessel smaller than that is representable", and that 7.5 µm was therefore
  "the smallest radius that reliably paints ≥1 grid cell per vessel". Two
  corrections:

  1. That reasoning does not make the *area fraction* grid-dependent. Sweeping
     ``dx = 30/20/12/8/6 µm`` on both production sections moves it by ≤1.03×,
     because point-sampling a fixed physical region on a lattice is an unbiased
     area estimator at any pitch. What the coarse grid destroys is the field's
     *geometry* — at ``dx = 30 µm`` only 0.26 pixels are painted per vessel, so
     most vessels vanish from the source entirely.
  2. A radius justified by the grid still cannot be used to produce a tissue
     measurement, and — the defect that actually moves the number — the disc is
     painted around every endothelial **cell**, so the area fraction tracks
     endothelial-cell density rather than vessel density. The two production
     sections have near-identical vessel densities (111 vs 104 /mm²) but 2.7×
     different endothelial densities, and their reported area fractions differed
     by exactly that ratio.

  :data:`VESSEL_RADIUS_UM` and the binary field below are therefore kept
  **bit-identical** (so no published result silently moves) but are no longer the
  route to the area fraction. :func:`segmentation_report` now reports the sub-grid,
  data-derived area fraction from :func:`reconstruct.vessel_geometry.
  vascular_area_fraction`, still scored against the measured GBM whole-section
  CD34 area fraction of 2.1 ± 0.3% (Angelucci C et al., Cells 2022;12(1):11,
  PMID 36611806, doi:10.3390/cells12010011).

Run:  ./.venv/bin/python -m tmesim.reconstruct.vessel_network
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from ..core import CellPopulation, CellType, Domain

__all__ = [
    "LINK_DISTANCE_UM",
    "VESSEL_RADIUS_UM",
    "PUBLISHED",
    "segment_vessels",
    "discrete_vascular_field",
    "link_distance_sweep",
    "segmentation_report",
]

#: Single-linkage distance (µm) below which two endothelial cells are taken to lie
#: on the same vessel. See the module docstring for the full justification.
LINK_DISTANCE_UM = 20.0

#: Radius (µm) of the nonzero source annulus around a segmented vessel's
#: endothelial cells. See the module docstring for the full justification and for
#: the dx = 12 µm resolution caveat.
VESSEL_RADIUS_UM = 7.5

#: Published human morphometry this segmentation is compared against. Every entry
#: is a measured value from the literature with its PMID; nothing here is estimated
#: by this project. The MVD anchor and the Poisson conversion are the same ones
#: ``results/vascular_validation.md`` uses, so the two reports cannot drift apart.
#:
#: **These are POPULATION REFERENCES, NOT TARGETS.** Every value here is a cohort
#: mean over tens of patients. This project reconstructs *one section of one
#: tumour*, and glioblastoma vasculature is abnormal and heterogeneous by
#: definition — that abnormality is diagnostic of the disease, not an error in the
#: reconstruction. A patient-specific reconstruction is therefore *expected* to
#: depart from these numbers, and pulling it towards them would erase exactly the
#: patient-specific signal the model exists to carry. Read the ``verdict`` field as
#: "inside / above / below the published cohort range", which is a description of
#: where this tumour sits, not a pass or a failure. The quantities that must
#: actually hold are (a) grid-independence — a number that claims to describe
#: tissue must not move when the mesh does — and (b) fidelity to the sample, i.e.
#: the field must reproduce the measured endothelial distribution including its
#: unevenness (:func:`reconstruct.vessel_geometry.vascular_heterogeneity`).
PUBLISHED: Dict[str, Dict] = {
    "mvd_per_mm2": {
        "value": 65.0, "sd": 30.0, "range": (35.0, 95.0),
        "what": "microvessel density, human GBM (n=20), CD31 whole-slide, "
                "unbiased systematic-uniform-random sampling",
        "source": "Marien KM et al., PLoS One 2016;11(9):e0161496",
        "pmid": "27583442", "doi": "10.1371/journal.pone.0161496",
        "measures": (
            "count of CD31-stained vessel OBJECTS per mm^2 of tissue, under an "
            "explicit connectivity rule quoted verbatim from the source's "
            "Methods: 'every CD31-positive object, no matter how small, should "
            "be counted, except suspected CD31-positive monocytes, macrophages "
            "and tumor cells. Furthermore, if CD31-positive objects were "
            "connected, they were considered a single object, while absence of "
            "staining defined two or more separate objects.'"),
        "commensurable_with_model": False,
        "caveat": (
            "The measurement's connectivity is CONTINUITY OF STAIN along the "
            "vessel wall; this module's is PROXIMITY OF NUCLEI at "
            "LINK_DISTANCE_UM = 20 um. Those cannot agree, and the direction of "
            "the disagreement is fixed by this package's own measured constant: "
            "vessel_geometry.EC_AXIAL_PITCH_UM = 50 um is the spacing of "
            "consecutive endothelial nuclei along one vessel (PMID 38723085). "
            "Linking at 20 um is provably below it, so consecutive nuclei on a "
            "SINGLE stained vessel are split into separate objects and MVD is "
            "over-counted BY CONSTRUCTION, in the same direction on every "
            "section. The source also excludes CD31+ monocytes, macrophages and "
            "tumour cells, which the marker-gated typing does not do, and notes "
            "independently that 'one vessel can appear as separate shapes in "
            "the two-dimensional section'. See "
            "results/vascular_segmentation_audit.md."),
    },
    "intervessel_spacing_um": {
        "value": 124.0, "range": (102.6, 169.0),
        "what": "1/sqrt(MVD) spacing implied by the MVD anchor above",
        "source": "derived from PMID 27583442; see vascular_validation.md §2",
        "pmid": "27583442", "doi": "10.1371/journal.pone.0161496",
    },
    "median_tumor_to_vessel_um": {
        "value": 58.0, "range": (48.2, 79.4),
        "what": "median tumour-cell-to-nearest-vessel distance implied by the MVD "
                "anchor under a random (Poisson) vessel arrangement, "
                "median = sqrt(ln2 / (pi q))",
        "source": "derived from PMID 27583442; see vascular_validation.md §2",
        "pmid": "27583442", "doi": "10.1371/journal.pone.0161496",
    },
    "vascular_area_fraction_pct": {
        "value": 2.1, "sd": 0.3, "range": (1.8, 2.4),
        "what": "CD34+ STAIN-POSITIVITY percentage, human IDH-wt GBM (n=30), "
                "Aperio Positive Pixel Count over 10 selected 40x fields in the "
                "TUMOUR CENTRE of 4 um sections",
        "source": "Tamma R, Ingravallo G, Annese T, d'Amati A, Lorusso L, "
                  "Ribatti D, Cells 2022;12(1):11",
        "pmid": "36611806", "doi": "10.3390/cells12010011",
        # --- provenance verified against the source's own Methods -----------
        "measures": (
            "fraction of image PIXELS scoring positive for CD34 immunostain. "
            "The source defines it verbatim as 'the number of positively "
            "stained pixels on the total pixels in the image'."),
        "commensurable_with_model": False,
        "caveat": (
            "NOT the same physical quantity as this module's vascular area "
            "fraction. The model paints the geometric cross-section of a vessel "
            "PROFILE (lumen + wall); CD34 stains only the endothelium, so the "
            "lumen contributes no positive pixels. The measurement is therefore "
            "a lower bound on vessel-profile area fraction, by an unquantified "
            "factor that grows with vessel caliber. Two further mismatches: it "
            "is 10 operator-selected 40x fields in the tumour centre, not a "
            "whole section or an unbiased sample; and the reported +/-0.3 is a "
            "STANDARD ERROR (the source's Statistical Analysis states 'means "
            "+/- SE', n=30), so the 1.8-2.4 band below is a +/-1 SE band on the "
            "cohort MEAN, not a cohort range. The implied cohort SD is "
            "0.3*sqrt(30) = 1.64, i.e. a +/-1 SD range of roughly 0.5-3.7. The "
            "band is left UNCHANGED here on purpose: widening it would convert "
            "the development section's 1.25% from a miss into a pass, which is "
            "a scoring decision for a human reviewer, not a silent edit. See "
            "results/vascular_segmentation_audit.md. Use "
            "'vessel_areal_fraction_pct' below for the commensurable comparison."),
        "dispersion_is": "standard error of the mean (n=30), NOT a cohort SD",
    },
    # The COMMENSURABLE areal-fraction anchor, from the SAME paper this module
    # already cites for MVD. Its stereological point-count rule counts a grid
    # point as vascular when it falls on "an endothelial cell OR a vessel lumen"
    # -- i.e. the whole vessel profile, which is exactly what
    # vessel_geometry.vascular_area_fraction paints. Reported alongside the CD34
    # stain-positivity entry above rather than replacing it, so that no existing
    # verdict moves and the two can be compared.
    "vessel_areal_fraction_pct": {
        "value": 4.70, "sd": 2.38, "range": (2.32, 7.08),
        "what": "microvessel AREAL FRACTION A_A, human GBM (n=20), CD31 "
                "whole-slide, systematic-uniform-random sampling, stereological "
                "point counting (grid point counts when it falls on an "
                "endothelial cell OR a vessel lumen)",
        "source": "Marien KM et al., PLoS One 2016;11(9):e0161496",
        "pmid": "27583442", "doi": "10.1371/journal.pone.0161496",
        "measures": (
            "areal fraction of VESSEL PROFILES (lumen included) per unit tissue "
            "area, by unbiased point counting. Same specimens, same sampling "
            "design and same denominator convention as the MVD anchor above."),
        "commensurable_with_model": True,
        "caveat": (
            "Cohort mean +/- 1 SD over 20 GBM patients; population context, not "
            "a target. The source itself reports this parameter as the LESS "
            "robust of its two: inter-ROI ICC for A_A in GBM was 0.6 against "
            "1.0 for the profile count, because point counting samples only a "
            "small area. Range shown is mean +/- 1 SD."),
    },
    "vessel_diameter_um": {
        "value": 15.0, "range": (15.0, 37.4),
        "what": "in-vivo vessel-size-imaging caliber in human brain "
                "(WM 15.0+/-5.1, thalamus 20.7+/-8.8, cortical GM 37.4+/-12.9 um); "
                "all brain tumours measured LARGER than normal-appearing tissue",
        "source": "Hsu YY et al., J Magn Reson Imaging 2009;30(5):1078-84",
        "pmid": "19856441", "doi": "10.1002/jmri.21960",
    },
    "glioma_vessel_diameter_um": {
        "value": 18.1, "sd": 2.4, "range": (12.8, 18.1),
        "what": "MRI vessel size in orthotopic U87 glioma, validated against IHC "
                "(18.1+/-2.4 um untreated; 12.8+/-2.0 um after bevacizumab). MOUSE.",
        "source": "Kording F et al., J Magn Reson Imaging 2014;40(6):1310-8",
        "pmid": "24390982", "doi": "10.1002/jmri.24513",
    },
    "diffusion_limit_um": {
        "value": 150.0, "range": (100.0, 180.0),
        "what": "distance from a vessel at which necrosis appears; viable tumour "
                "cord rim never exceeded ~180 um (human bronchial carcinoma)",
        "source": "Thomlinson RH & Gray LH, Br J Cancer 1955;9(4):539-549",
        "pmid": "13304213", "doi": "10.1038/bjc.1955.55",
    },
}


# --------------------------------------------------------------------------- #
# 1. Segmentation
# --------------------------------------------------------------------------- #
def _endothelial_xy(cells: CellPopulation) -> np.ndarray:
    endo = np.asarray(cells.cell_type) == int(CellType.ENDOTHELIAL)
    return np.column_stack([
        np.asarray(cells.x, dtype=np.float64)[endo],
        np.asarray(cells.y, dtype=np.float64)[endo],
    ])


def _single_linkage_labels(xy: np.ndarray, eps: float) -> np.ndarray:
    """Connected components of the ``eps``-neighbourhood graph.

    Identical to ``sklearn.cluster.DBSCAN(eps=eps, min_samples=1)`` — every point
    is a core point, so the clustering is exactly single-linkage at threshold
    ``eps`` with no points labelled noise. Implemented on scipy so the package
    picks up no new dependency.
    """
    if xy.shape[0] == 0:
        return np.zeros(0, dtype=int)
    tree = cKDTree(xy)
    graph = tree.sparse_distance_matrix(tree, float(eps), output_type="coo_matrix")
    _, labels = connected_components(graph, directed=False)
    return np.asarray(labels, dtype=int)


def _median_finite_or_inf(a: np.ndarray) -> float:
    """Median over the finite entries; ``+inf`` when every entry is infinite.

    Used for the cluster aspect ratio, where a perfectly collinear object gives
    ``sv[1] == 0`` and therefore an infinite ratio. Returning ``nan`` there would
    be wrong twice over: it is not missing data, and the case it describes
    (perfectly straight vessels) is the most elongated one possible.
    """
    a = np.asarray(a, dtype=np.float64)
    finite = a[np.isfinite(a)]
    if finite.size:
        return float(np.median(finite))
    return float("inf") if a.size else float("nan")


def segment_vessels(
    cells: CellPopulation,
    domain: Domain,
    link_distance_um: float = LINK_DISTANCE_UM,
    tissue_area_mm2: Optional[float] = None,
) -> Dict:
    """Cluster the REAL endothelial cells into DISCRETE putative vessel objects.

    Parameters
    ----------
    cells : CellPopulation
        The reconstructed real population. Only ``CellType.ENDOTHELIAL`` members
        are used; their measured ``x``/``y`` (µm) are the input.
    domain : Domain
        Supplies ``dx`` and ``tissue_mask``, used only to convert counts into
        areal densities. The clustering itself is computed from cell coordinates
        and is therefore **exactly invariant to** ``dx``.
    link_distance_um : float
        Single-linkage threshold (µm). Default :data:`LINK_DISTANCE_UM` = 20 µm.
    tissue_area_mm2 : float, optional
        Tissue area to divide counts by. Defaults to the rasterised
        ``domain.tissue_mask`` area. (At the default geometric mask this is
        dx-stable to well under a percent, but pass it explicitly if you want a
        dx-free denominator.)

    Returns
    -------
    dict
        ``n_vessels``, ``n_endothelial_cells``, per-vessel arrays
        (``cells_per_vessel``, ``centroids``, ``extent_um``, ``aspect_ratio``),
        the vessel-to-nearest-vessel distance distribution
        (``intervessel_nn_um`` + percentiles), the tumour-cell-to-nearest-vessel
        distance distribution (``tumor_to_vessel_um`` + percentiles), ``mvd_per_mm2``,
        and ``labels`` (the per-endothelial-cell cluster id, aligned with ``xy``).

    Notes
    -----
    Distances are measured to the nearest **endothelial cell of another vessel**,
    not to the cluster centroid, because a vessel is a spatially extended object
    and its centroid is not where the drug leaks from.
    """
    xy = _endothelial_xy(cells)
    n_endo = int(xy.shape[0])
    labels = _single_linkage_labels(xy, float(link_distance_um))
    n_vessels = int(labels.max() + 1) if n_endo else 0

    if tissue_area_mm2 is None:
        tissue = np.asarray(domain.tissue_mask, dtype=bool)
        tissue_area_mm2 = float(tissue.sum()) * float(domain.dx) ** 2 / 1e6

    counts = (np.bincount(labels, minlength=max(n_vessels, 1)).astype(int)
              if n_endo else np.zeros(0, dtype=int))
    centroids = np.zeros((n_vessels, 2), dtype=np.float64)
    extent = np.zeros(n_vessels, dtype=np.float64)
    aspect = np.ones(n_vessels, dtype=np.float64)
    for k in range(n_vessels):
        pts = xy[labels == k]
        centroids[k] = pts.mean(axis=0)
        if pts.shape[0] >= 2:
            # Max pairwise separation = the object's length along its own axis.
            d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1)
            extent[k] = float(d.max())
            # Aspect ratio from the principal axes of the point cloud.
            c = pts - pts.mean(axis=0)
            sv = np.linalg.svd(c, compute_uv=False)
            aspect[k] = float(sv[0] / sv[1]) if sv.size > 1 and sv[1] > 1e-9 else float("inf")

    # --- vessel-to-nearest-other-vessel distance (cell-to-cell, cross-label) --- #
    inter_nn = np.full(n_vessels, np.nan)
    if n_vessels > 1:
        tree = cKDTree(xy)
        # Query enough neighbours to be sure of finding one in a different vessel.
        kmax = int(min(max(counts.max() + 1, 2), n_endo))
        dists, idx = tree.query(xy, k=kmax)
        dists = np.atleast_2d(dists)
        idx = np.atleast_2d(idx)
        per_cell = np.full(n_endo, np.inf)
        for j in range(dists.shape[1]):
            other = labels[idx[:, j]] != labels
            upd = other & (dists[:, j] < per_cell)
            per_cell[upd] = dists[upd, j]
        for k in range(n_vessels):
            v = per_cell[labels == k]
            v = v[np.isfinite(v)]
            if v.size:
                inter_nn[k] = float(v.min())

    # --- tumour-cell-to-nearest-vessel distance (the delivery-relevant one) ---- #
    tumor = np.asarray(cells.cell_type) == int(CellType.TUMOR)
    tumor_xy = np.column_stack([np.asarray(cells.x)[tumor], np.asarray(cells.y)[tumor]])
    tumor_d = np.zeros(0)
    if tumor_xy.shape[0] and n_endo:
        tumor_d = cKDTree(xy).query(tumor_xy, k=1)[0].astype(np.float64)

    def _pct(a: np.ndarray) -> Dict[str, float]:
        a = np.asarray(a, dtype=np.float64)
        a = a[np.isfinite(a)]
        if a.size == 0:
            return {k: float("nan") for k in
                    ("median", "mean", "p10", "p25", "p75", "p90", "p95", "p99", "max")}
        return {
            "median": float(np.median(a)), "mean": float(a.mean()),
            "p10": float(np.percentile(a, 10)), "p25": float(np.percentile(a, 25)),
            "p75": float(np.percentile(a, 75)), "p90": float(np.percentile(a, 90)),
            "p95": float(np.percentile(a, 95)), "p99": float(np.percentile(a, 99)),
            "max": float(a.max()),
        }

    mvd = float(n_vessels) / tissue_area_mm2 if tissue_area_mm2 > 0 else float("nan")
    multi = counts >= 2
    return {
        "link_distance_um": float(link_distance_um),
        "n_endothelial_cells": n_endo,
        "n_vessels": n_vessels,
        "labels": labels,
        "xy": xy,
        "cells_per_vessel": counts,
        "centroids": centroids,
        "extent_um": extent,
        "aspect_ratio": aspect,
        "tissue_area_mm2": float(tissue_area_mm2),
        "mvd_per_mm2": mvd,
        "spacing_from_mvd_um": (1e3 / np.sqrt(mvd)) if mvd > 0 else float("nan"),
        "singleton_fraction": float((counts == 1).mean()) if n_vessels else float("nan"),
        "mean_cells_per_vessel": float(counts.mean()) if n_vessels else float("nan"),
        "n_vessels_ge4_cells": int((counts >= 4).sum()),
        "median_extent_multi_um": (float(np.median(extent[multi])) if multi.any()
                                   else float("nan")),
        # Aspect ratio is taken over the FINITE values only, because any 2-point
        # object is exactly collinear (sv[1] == 0 -> inf) and carries no shape
        # information; this keeps the statistic comparable to the "median aspect
        # ratio of multi-cell objects" in results/vascular_validation.md §3.3.
        # If EVERY multi-cell object is perfectly collinear the finite subset is
        # empty, and the honest answer is +inf (maximally elongated), not nan.
        "median_aspect_multi": _median_finite_or_inf(aspect[multi]) if multi.any()
                               else float("nan"),
        "intervessel_nn_um": inter_nn,
        "intervessel_nn": _pct(inter_nn),
        "tumor_to_vessel_um": tumor_d,
        "tumor_to_vessel": _pct(tumor_d),
        "tumor_frac_beyond_100um": (float((tumor_d > 100).mean()) if tumor_d.size
                                    else float("nan")),
        "tumor_frac_beyond_150um": (float((tumor_d > 150).mean()) if tumor_d.size
                                    else float("nan")),
        "n_tumor": int(tumor.sum()),
    }


def link_distance_sweep(
    cells: CellPopulation,
    domain: Domain,
    link_distances: Sequence[float] = (10.0, 15.0, 20.0, 25.0, 30.0, 40.0),
) -> List[Dict]:
    """Re-segment at several linking distances. The linking distance is a free
    parameter, so the honest presentation is the sweep, not the single value."""
    out = []
    for d in link_distances:
        v = segment_vessels(cells, domain, link_distance_um=float(d))
        out.append({k: v[k] for k in (
            "link_distance_um", "n_vessels", "mvd_per_mm2", "mean_cells_per_vessel",
            "singleton_fraction", "spacing_from_mvd_um", "n_vessels_ge4_cells",
            "median_extent_multi_um")}
            | {"intervessel_nn_median_um": v["intervessel_nn"]["median"],
               "tumor_to_vessel_median_um": v["tumor_to_vessel"]["median"],
               "tumor_frac_beyond_100um": v["tumor_frac_beyond_100um"]})
    return out


# --------------------------------------------------------------------------- #
# 2. The discrete source field
# --------------------------------------------------------------------------- #
def discrete_vascular_field(
    domain: Domain,
    vessels: Dict,
    radius_um: float = VESSEL_RADIUS_UM,
    amplitude: str = "unit",
    match_total_to: Optional[np.ndarray] = None,
    geometry: str = "binary",
    radius_mode: str = "data",
) -> np.ndarray:
    """Source field (ny, nx) that is nonzero ONLY within ``radius_um`` of a
    segmented vessel, and EXACTLY zero everywhere else.

    A grid cell is "vessel" when its **centre** lies within ``radius_um`` of any
    endothelial cell belonging to a segmented vessel object. The result is then
    masked to ``domain.tissue_mask``. There is no Gaussian tail and therefore no
    far-field floor: ``field[d > radius_um] == 0.0`` exactly.

    ``geometry='subgrid'`` — the FIXED geometry
    -------------------------------------------
    **ADDITIVE AND OPT-IN; the default is bit-identical to the historical field.**
    With ``geometry='subgrid'`` the field is instead the fractional area coverage
    of :func:`reconstruct.vessel_geometry.subgrid_coverage_field`: one footprint
    per segmented **vessel object** (not per endothelial cell), at a radius derived
    from that object's own endothelial content, painted as a partial-area weight so
    a vessel smaller than a grid cell is still represented. Values are then in
    ``[0, 1]`` continuously rather than binary. ``radius_um`` is ignored unless
    ``radius_mode='constant'``. See ``results/vessel_geometry_fix.md``.

    The default is left at ``'binary'`` because ``vascular_density_from_cells``
    (and therefore every reconstruction and every result already published against
    it) routes through this function; flipping it is a deliberate re-run, not a
    silent one.

    Normalisation — what the numbers MEAN
    -------------------------------------
    ``amplitude='unit'`` (default) sets the field to **1.0 inside a vessel and 0.0
    outside**. This is a **relative** field, in exactly the same sense the existing
    max-normalised KDE is relative, but with a cleaner referent: 1.0 means *"a
    perfused vessel wall leaking at the reference rate"*, whereas 1.0 in the KDE
    means *"the densest endothelial-KDE spot anywhere in this particular section"*.
    Neither is an absolute perfused fraction, and no result built on either should
    be read as an absolute delivery prediction. The engine's delivery term is
    ``vascular_density * bbb_penetrance * plasma``, so with ``amplitude='unit'``
    the physical reading is: *drug enters the interstitium at the plasma-driven
    leak rate at the vessel wall and nowhere else*, which is what a leaky-vessel
    boundary condition actually says.

    ``match_total_to=<array>`` additionally rescales the field so its
    tissue-integrated total equals that of the supplied field (normally the KDE).
    This is the **controlled** comparison: it holds the total drug entering the
    tissue fixed and changes only *where* it enters, so any difference in outcome
    is attributable to geometry rather than to dose. Use it as a control, and say
    which one a number came from.

    Parameters
    ----------
    domain : Domain
        Grid + tissue mask.
    vessels : dict
        Output of :func:`segment_vessels`.
    radius_um : float
        Vessel radius (µm). Default :data:`VESSEL_RADIUS_UM` = 7.5 µm — see the
        module docstring for the citations and for the dx = 12 µm resolution
        caveat (a single grid cell already has an equivalent radius of 6.8 µm).
    amplitude : {'unit', 'per_vessel_cells'}
        ``'unit'`` — every vessel pixel is 1.0. ``'per_vessel_cells'`` — a vessel
        pixel carries the number of endothelial cells in its vessel object
        divided by the section's max, i.e. bigger vessels leak proportionally
        more. ``'unit'`` is the default because vessel *caliber* is not measured
        by Xenium: the endothelial cell count per object is dominated by the
        section plane and the PECAM1 gate, so weighting by it would import an
        unvalidated assumption.
    match_total_to : ndarray, optional
        Rescale so ``field.sum()`` over tissue equals this field's sum over
        tissue.

    Returns
    -------
    np.ndarray
        (ny, nx) non-negative field, zero outside the vessel radius and outside
        tissue.
    """
    geom_mode = str(geometry).strip().lower()
    if geom_mode not in {"binary", "subgrid"}:
        raise ValueError(f"unknown geometry {geometry!r}; "
                         "expected 'binary' or 'subgrid'")

    field = domain.zeros()
    xy = np.asarray(vessels["xy"], dtype=np.float64)
    if xy.shape[0] == 0:
        return field

    if geom_mode == "subgrid":
        from .vessel_geometry import subgrid_coverage_field, vessel_geometry
        if amplitude not in {"unit", "per_vessel_cells"}:
            raise ValueError(f"unknown amplitude {amplitude!r}; "
                             "expected 'unit' or 'per_vessel_cells'")
        geom = vessel_geometry(
            vessels, radius_mode=radius_mode,
            constant_radius_um=(radius_um if radius_mode == "constant" else None))
        field = subgrid_coverage_field(domain, geom)
        if amplitude == "per_vessel_cells":
            raise ValueError(
                "amplitude='per_vessel_cells' is not defined for "
                "geometry='subgrid': the sub-grid footprint already carries the "
                "object's endothelial content through its radius, so scaling by "
                "the cell count again would double-count it")
        tissue = np.asarray(domain.tissue_mask, dtype=bool)
        if match_total_to is not None:
            ref = np.asarray(match_total_to, dtype=np.float64)
            ref_total = float(ref[tissue].sum())
            cur_total = float(field[tissue].sum())
            if cur_total > 0 and ref_total > 0:
                field = field * (ref_total / cur_total)
        return field

    ny, nx = domain.shape
    dx = float(domain.dx)
    # Grid-cell CENTRES, matching Domain.to_index's convention.
    gx = (np.arange(nx, dtype=np.float64) + 0.5) * dx
    gy = (np.arange(ny, dtype=np.float64) + 0.5) * dx
    GX, GY = np.meshgrid(gx, gy)
    grid_pts = np.column_stack([GX.ravel(), GY.ravel()])

    tree = cKDTree(xy)
    if amplitude == "per_vessel_cells":
        counts = np.asarray(vessels["cells_per_vessel"], dtype=np.float64)
        labels = np.asarray(vessels["labels"], dtype=int)
        d, j = tree.query(grid_pts, k=1, distance_upper_bound=float(radius_um))
        inside = np.isfinite(d)
        vals = np.zeros(grid_pts.shape[0], dtype=np.float64)
        vals[inside] = counts[labels[j[inside]]]
        if vals.max() > 0:
            vals /= vals.max()
        field = vals.reshape(ny, nx)
    elif amplitude == "unit":
        hit = tree.query_ball_point(grid_pts, r=float(radius_um),
                                    return_length=True) > 0
        field = hit.reshape(ny, nx).astype(np.float64)
    else:
        raise ValueError(f"unknown amplitude {amplitude!r}; "
                         "expected 'unit' or 'per_vessel_cells'")

    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    field = field * tissue

    if match_total_to is not None:
        ref = np.asarray(match_total_to, dtype=np.float64)
        ref_total = float(ref[tissue].sum())
        cur_total = float(field[tissue].sum())
        if cur_total > 0 and ref_total > 0:
            field = field * (ref_total / cur_total)
    return field


# --------------------------------------------------------------------------- #
# 3. Scoring the segmentation against published morphometry
# --------------------------------------------------------------------------- #
def segmentation_report(vessels: Dict, domain: Domain,
                        radius_um: float = VESSEL_RADIUS_UM,
                        area_fraction_geometry: str = "subgrid",
                        radius_mode: str = "data") -> Dict:
    """Score a :func:`segment_vessels` result against :data:`PUBLISHED`.

    Returns one row per comparable quantity with the reconstructed value, the
    published value + range, the ratio, and a PASS/HIGH/LOW verdict (PASS = inside
    the published range).

    Parameters
    ----------
    radius_um : float
        Radius for the LEGACY binary field, which is still computed and reported
        as ``legacy_binary_area_fraction_pct`` so the two can be compared.
    area_fraction_geometry : {'subgrid', 'binary'}
        Which geometry supplies the headline ``vessel_area_fraction_pct`` that is
        scored against the published 2.1 ± 0.3 %. Default ``'subgrid'`` — the
        sub-grid, data-derived footprint of
        :func:`reconstruct.vessel_geometry.vascular_area_fraction`, which is
        grid-independent by construction and uses no hardcoded radius.
        ``'binary'`` reproduces the historical (per-endothelial-cell, fixed
        7.5 µm, centre-in-disc) number.
    radius_mode : {'data', 'stereological', 'constant'}
        Passed through to :func:`reconstruct.vessel_geometry.vessel_geometry`.
    """
    from .vessel_geometry import vascular_area_fraction

    field = discrete_vascular_field(domain, vessels, radius_um=radius_um)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    legacy_area_frac_pct = (
        100.0 * float((field > 0)[tissue].sum()) / float(tissue.sum())
        if tissue.any() else float("nan"))

    geom_mode = str(area_fraction_geometry).strip().lower()
    if geom_mode not in {"subgrid", "binary"}:
        raise ValueError(f"unknown area_fraction_geometry {area_fraction_geometry!r};"
                         " expected 'subgrid' or 'binary'")
    af = vascular_area_fraction(
        domain, vessels, radius_mode=radius_mode,
        constant_radius_um=(radius_um if radius_mode == "constant" else None))
    af.pop("coverage", None)
    area_frac_pct = (af["area_fraction_pct"] if geom_mode == "subgrid"
                     else legacy_area_frac_pct)

    rows = []

    def _row(key: str, recon: float, units: str):
        pub = PUBLISHED[key]
        lo, hi = pub["range"]
        if not np.isfinite(recon):
            verdict = "n/a"
        elif recon < lo:
            verdict = "LOW"
        elif recon > hi:
            verdict = "HIGH"
        else:
            verdict = "PASS"
        rows.append({
            "quantity": pub["what"], "key": key, "units": units,
            "reconstructed": float(recon),
            "published": pub["value"], "published_range": (float(lo), float(hi)),
            "ratio": float(recon / pub["value"]) if pub["value"] else float("nan"),
            # NB: 'HIGH'/'LOW' mean "outside the published COHORT range", not
            # "wrong". See the note on PUBLISHED: these are population references
            # for context, not targets for a single-patient reconstruction.
            "verdict": verdict,
            "reference_kind": "population cohort mean — context, not a target",
            "pmid": pub["pmid"], "doi": pub["doi"],
            "source": pub["source"],
            # Is the model computing the SAME physical quantity the source
            # measured? Where this is False the ratio above is not an error
            # measurement and must not be read as one. ``None`` means the
            # question has not been adjudicated for that row.
            "measures": pub.get("measures"),
            "commensurable": pub.get("commensurable_with_model"),
            "caveat": pub.get("caveat"),
        })

    _row("mvd_per_mm2", vessels["mvd_per_mm2"], "vessels/mm^2")
    _row("intervessel_spacing_um", vessels["spacing_from_mvd_um"], "um")
    _row("median_tumor_to_vessel_um", vessels["tumor_to_vessel"]["median"], "um")
    _row("vascular_area_fraction_pct", area_frac_pct, "% of tissue area")
    # The SAME reconstructed number scored against the commensurable anchor.
    # Both rows are emitted so the size of the modality correction is always
    # visible next to the historical comparison; neither replaces the other.
    _row("vessel_areal_fraction_pct", area_frac_pct, "% of tissue area")

    # Poisson (complete-spatial-randomness) expectation AT THE RECONSTRUCTION'S OWN
    # density — the test that ``vascular_validation.md`` shows the KDE fails. Real
    # GBM is clustered, so a real section must have a HEAVIER far tail than this.
    q = vessels["mvd_per_mm2"] / 1e6  # vessels per um^2
    csr_median = float(np.sqrt(np.log(2.0) / (np.pi * q))) if q > 0 else float("nan")
    csr_beyond_100 = float(np.exp(-np.pi * (100.0 ** 2) * q)) if q > 0 else float("nan")
    obs_median = vessels["tumor_to_vessel"]["median"]
    obs_beyond_100 = vessels["tumor_frac_beyond_100um"]

    return {
        "rows": rows,
        "vessel_area_fraction_pct": area_frac_pct,
        "area_fraction_geometry": geom_mode,
        # The historical per-endothelial-cell, fixed-7.5 µm, centre-in-disc number,
        # kept so the size of the correction is always visible next to the fix.
        "legacy_binary_area_fraction_pct": legacy_area_frac_pct,
        "area_fraction_detail": af,
        "radius_um": float(radius_um),
        # Only meaningful for the legacy binary field; the sub-grid painter has no
        # minimum representable radius, which is the point of it.
        "grid_min_representable_radius_um": float(np.sqrt(domain.dx ** 2 / np.pi)),
        "csr_median_um": csr_median,
        "csr_frac_beyond_100um": csr_beyond_100,
        "clustering_index_median": (float(obs_median / csr_median)
                                    if np.isfinite(csr_median) and csr_median > 0
                                    else float("nan")),
        "clustering_index_tail": (float(obs_beyond_100 / csr_beyond_100)
                                  if np.isfinite(csr_beyond_100) and csr_beyond_100 > 0
                                  else float("nan")),
    }


# --------------------------------------------------------------------------- #
# Self-test on REAL Xenium GBM data
# --------------------------------------------------------------------------- #
if __name__ == "__main__":  # pragma: no cover
    import time

    from .real_pipeline import reconstruct_real

    t0 = time.time()
    try:
        micro = reconstruct_real("xenium", dx_um=12.0)
    except Exception as exc:
        print(f"load failed ({exc!r}); retrying once ...")
        time.sleep(3)
        micro = reconstruct_real("xenium", dx_um=12.0)
    print(f"REAL Xenium GBM reconstructed in {time.time()-t0:.0f}s  "
          f"grid={micro.domain.shape}  dx={micro.domain.dx:g} um")

    v = segment_vessels(micro.cells, micro.domain)
    print(f"\n{v['n_endothelial_cells']:,} endothelial cells -> "
          f"{v['n_vessels']:,} discrete vessels at "
          f"{v['link_distance_um']:.0f} um linking")
    print(f"  MVD {v['mvd_per_mm2']:.1f}/mm^2   "
          f"cells/vessel {v['mean_cells_per_vessel']:.2f}   "
          f"singletons {v['singleton_fraction']*100:.0f}%")
    print(f"  intervessel NN median {v['intervessel_nn']['median']:.1f} um   "
          f"tumour->vessel median {v['tumor_to_vessel']['median']:.1f} um   "
          f">100 um {v['tumor_frac_beyond_100um']*100:.1f}%")

    print("\nlink-distance sweep:")
    for r in link_distance_sweep(micro.cells, micro.domain):
        print(f"  {r['link_distance_um']:5.0f} um -> {r['n_vessels']:6,} vessels  "
              f"MVD {r['mvd_per_mm2']:6.1f}/mm^2  "
              f"cells/vessel {r['mean_cells_per_vessel']:.2f}  "
              f"interNN {r['intervessel_nn_median_um']:6.1f} um")

    f = discrete_vascular_field(micro.domain, v)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    kde = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    print(f"\ndiscrete field: {(f>0)[tissue].sum():,} vessel px "
          f"({100*(f>0)[tissue].mean():.2f}% of tissue area), "
          f"mean over tissue {f[tissue].mean():.4f}")
    print(f"KDE field     : {100*(kde>0)[tissue].mean():.2f}% of tissue nonzero, "
          f"mean over tissue {kde[tissue].mean():.4f}")

    # The whole point: the discrete field is EXACTLY zero away from vessels.
    assert f.min() == 0.0, "discrete field must reach exactly zero"
    assert (f[tissue] == 0).mean() > 0.5, (
        "discrete field should leave most of the tissue with NO source")
    assert (kde[tissue] == 0).mean() < 0.05, (
        "sanity: the KDE field should be nonzero almost everywhere (the defect)")
    assert f[~tissue].sum() == 0.0, "no source outside tissue"

    rep = segmentation_report(v, micro.domain)
    print("\nvs published human GBM morphometry:")
    for r in rep["rows"]:
        print(f"  {r['key']:32s} {r['reconstructed']:9.2f} {r['units']:16s} vs "
              f"{r['published']:.2f} [{r['published_range'][0]:.1f}-"
              f"{r['published_range'][1]:.1f}]  {r['verdict']}  (PMID {r['pmid']})")
    print(f"\nclustering index (obs/CSR): median {rep['clustering_index_median']:.2f}x, "
          f">100 um tail {rep['clustering_index_tail']:.2f}x")
    print("\nSELF-TEST PASSED (real Xenium GBM data)")
