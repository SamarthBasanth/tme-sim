"""Oxygen SURROGATES that are not a solved oxygen tension — and are better at it.

WHY THIS MODULE EXISTS
----------------------
``results/state_validation.md`` §1.2 falsified the reaction-diffusion oxygen solve
**as value-add** against measured hypoxia markers in 816,769 real cells (10x Xenium
Human Brain Cancer GBM FFPE, IO panel). Sign-aligned to oxygen, per-cell Spearman
rho against the measured VEGFA/HIF1A hypoxia score:

    distance to nearest ENDOTHELIAL cell   -0.1503   <- one cKDTree query
    the PDE's own UNDIFFUSED source (KDE)  -0.1270   <- the solver's input
    the SOLVED oxygen field (the PDE)      -0.1132   <- the sparse linear solve
    local cell density                     -0.0192   <- the sanity floor

Two facts, not one. (1) A distance transform beats the solve by 33%. (2) The solve
is WORSE than its own input: diffusing the source destroys predictive information
rather than adding it. The nested partials say the same thing — conditioning on
distance leaves the PDE only -0.066, conditioning on the PDE leaves distance
-0.119, so the PDE is the redundant member of the pair. And it is grid-converged
(-0.1132 at dx = 30 / 20 / 12 um), so refinement cannot rescue it.

This module therefore supplies the two HONEST alternatives:

``oxygen_from_vessel_distance``
    A monotone map from "distance to the nearest perfused vessel" to a
    dimensionless [0, 1] **oxygen surrogate**. It is NOT a solved oxygen tension
    and this module never claims it is; it is a geometric proxy whose only defence
    is that it out-predicts the solve on the one section where we can measure.

``hypoxia_from_markers``
    Hypoxia read DIRECTLY off measured hypoxia-response transcripts, with no
    physics at all. This is the only route on this list that is a MEASUREMENT.
    Where the panel does not carry a hypoxia-response gene it returns ``None``
    with an explicit reason instead of silently manufacturing a field — the
    319-gene brain panel carries no VEGFA, CA9, ADM, NDRG1 or SLC2A1, so on that
    section hypoxia is UNMEASURED and must be reported as such.

WHAT IS AND IS NOT CLAIMED
--------------------------
* The surrogate is **relative and dimensionless**, exactly like the PDE's output.
  Neither is a pO2 in mmHg. ``OxygenConfig.consumption`` has never been fitted to a
  measured oxygen tension (``physics.oxygen`` docstring), so the absolute level of
  the solved field is uncalibrated too. Demoting the PDE does not calibrate
  anything; it only stops paying for a solve that loses to a kd-tree.
* The surrogate is **monotone in distance by construction**. That has a sharp
  consequence for the benchmark: *every* strictly-monotone map of distance has the
  SAME Spearman rho. The mapping's shape and its length scale therefore cannot be
  tuned to improve the rank correlation — they only set where the HYPOXIC /
  QUIESCENT / VIABLE band edges fall. This is why the calibration below is
  anchored to published morphometry rather than fitted to the marker data: fitting
  it would be unfalsifiable decoration.
* VEGFA is a hypoxia-**RESPONSE** gene, not an oxygen tension. It is the
  transcriptional output of a pathway with its own kinetics and its own
  non-hypoxic inducers. ``hypoxia_from_markers`` returns a measured
  hypoxia-response field, and the docstring says exactly that.

THE MAPPING, AND HOW IT IS CALIBRATED
-------------------------------------
Steady 1-D diffusion away from a line source into consuming tissue decays with a
single length scale, so the minimal dimensionally-sane monotone map is

    O_surrogate(d) = exp(-d / L)                      d = um to nearest vessel

with ONE free constant, the oxygen diffusion length ``L`` (um). It is fixed by a
published anchor rather than fitted:

* Thomlinson & Gray, Br J Cancer 1955;9(4):539-549 (PMID 13304213,
  doi:10.1038/bjc.1955.55) measured, in human lung carcinoma, that tumour cords
  are viable out to ~100-180 um from the stroma and necrotic beyond it — the
  original oxygen-diffusion-distance measurement. :data:`HYPOXIA_ONSET_UM` = 150
  um is the centre of that band.
* The band edge must line up with the banding the rest of the codebase already
  uses: ``cells.agents.update_cell_states`` calls a cell HYPOXIC below
  ``OxygenConfig.hypoxia_threshold`` = 0.2. Requiring
  ``exp(-HYPOXIA_ONSET_UM / L) == 0.2`` gives

      L = 150 / ln(5) = 93.2 um

  so a cell 150 um from the nearest vessel lands exactly on the existing hypoxia
  threshold, 2x that threshold (the QUIESCENT/VIABLE edge) falls at 85.9 um, and
  the surrogate needs no new tuning knob anywhere downstream.
* Diagnostic, not calibration — and it does NOT agree: the PDE's own effective
  diffusion length is ``sqrt(D / (consumption * rho_tum))``
  (:func:`pde_diffusion_length_um`). At the default ``OxygenConfig`` and the real
  sections' measured tumour density that is ~310 um on the IO section and ~800 um
  on the brain section, i.e. **2-5x the published oxygen-diffusion distance**.
  The uncalibrated ``OxygenConfig.consumption``
  therefore makes the solved field far smoother than tissue oxygen is measured to
  be, which is a mechanistic account of §1.2's finding that the diffusion step
  DESTROYS predictive information: it low-pass-filters the source at a length scale
  the tissue does not have. The surrogate does not inherit that error, but neither
  does it fix the underlying problem — ``consumption`` is still unfitted.
* Independent support that the vessel-distance axis is the physically right one:
  Helmlinger et al., Nat Med 1997;3(2):177-182 (PMID 9018236,
  doi:10.1038/nm0297-177) measured interstitial pO2 profiles BETWEEN ADJACENT
  VESSELS at <=10 um resolution in vivo. That paper also documents the honest
  limit of the proxy: profile shapes are heterogeneous and some well-perfused
  vessels are themselves hypoxic, so distance-to-a-vessel is an upper bound on
  what geometry alone can know.

``L`` is a keyword argument everywhere, and :func:`vessel_distance_um` is exposed
separately so a caller can apply any other monotone map without touching the rank
statistics.

COST
----
The PDE is a sparse direct solve on ny*nx unknowns (905 x 956 = 865,180 at the
production dx = 12 um). The surrogate is one ``cKDTree`` build over the
endothelial cells plus one nearest-neighbour query per grid point. See
``results/oxygen_demotion.md`` for the measured wall-clock comparison.

WHAT IS NOT TOUCHED
-------------------
``core``, ``sim.engine`` and ``physics.oxygen`` are unmodified. This module is
purely additive; the shipping default oxygen route is still the PDE
(``reconstruct.real_pipeline.reconstruct_real(oxygen_mode='pde')``).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy import sparse as _sp
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from ..core import CellPopulation, CellType, Domain, Microenvironment, OxygenConfig

__all__ = [
    "HYPOXIA_ONSET_UM",
    "HYPOXIA_MARKER_GENES",
    "OXYGEN_MODES",
    "VESSEL_SOURCES",
    "diffusion_length_um",
    "hypoxia_from_markers",
    "oxygen_from_vessel_distance",
    "pde_diffusion_length_um",
    "vessel_distance_um",
]

#: Distance (um) from a perfused vessel at which tissue is taken to have reached
#: the HYPOXIC band. Thomlinson & Gray, Br J Cancer 1955;9(4):539-549 (PMID
#: 13304213, doi:10.1038/bjc.1955.55): tumour cords in human lung carcinoma are
#: viable to ~100-180 um from the stroma and necrotic beyond. 150 um is the centre
#: of that measured band. This is the ONLY calibration constant in the mapping,
#: and it is a published morphometric distance, not a fit to the marker data.
HYPOXIA_ONSET_UM = 150.0

#: Canonical HIF-target hypoxia-RESPONSE transcripts, in the order preferred when
#: a panel carries several. NOT an oxygen tension: these are the output of a
#: transcriptional programme with its own kinetics and its own non-hypoxic
#: inducers (see the module docstring).
#:
#: HIF1A is deliberately ABSENT. Its mRNA ranked 249/480 at rho = +0.0085 against
#: the simulated field (``results/state_validation.md`` §1.1) for a known
#: mechanistic reason — oxygen sensing acts on the PROTEIN via O2-dependent prolyl
#: hydroxylation and VHL-mediated proteolysis (Ivan et al., Science
#: 2001;292(5516):464-468, PMID 11292862, doi:10.1126/science.1059817; Jaakkola et
#: al., Science 2001;292(5516):468-472, PMID 11292861, doi:10.1126/science.1059796)
#: — so its transcript has no obligation to track oxygen, and including it dilutes
#: the composite from rho = -0.148 (VEGFA alone) to -0.113. Pass it explicitly via
#: ``gene_set=`` if you want the diluted historical score back.
HYPOXIA_MARKER_GENES: Tuple[str, ...] = (
    "VEGFA",     # Fischer et al., Brain Pathol 2005 (PMID 16389942): induced in
                 # perinecrotic pseudopalisading glioma cells specifically.
    "CA9",
    "ADM",
    "NDRG1",
    "SLC2A1",
    "LDHA",
    "PGK1",
    "ANGPTL4",
)

#: The three oxygen routes this project can take. ``'pde'`` is and remains the
#: shipping default (see :func:`reconstruct.real_pipeline.reconstruct_real`).
OXYGEN_MODES: Tuple[str, ...] = ("pde", "vessel_distance", "markers")

#: Where the "vessel" points come from for the distance transform.
#:
#: ``'endothelial'``
#:     Every measured ``CellType.ENDOTHELIAL`` cell is a vessel point. This is the
#:     exact predictor that scored -0.1503 in ``results/state_validation.md``, so
#:     it is the DEFAULT: it is the one that was actually benchmarked.
#: ``'vessel_network'``
#:     Endothelial cells first clustered into DISCRETE vessel objects by
#:     :func:`reconstruct.vessel_network.segment_vessels`, keeping only objects
#:     with at least ``min_cells_per_vessel`` members. More physical — a lone
#:     PECAM1+ cell is more likely a segmentation artefact or a non-perfused
#:     sprout than a perfused vessel — but it discards real cells on the strength
#:     of a free linking-distance parameter, so it is offered, not defaulted.
VESSEL_SOURCES: Tuple[str, ...] = ("endothelial", "vessel_network")

#: Minimum endothelial cells in a segmented object for it to count as a perfused
#: vessel under ``source='vessel_network'``. 2 only excludes true singletons; the
#: segmentation's own report treats >=4 cells as the "well-resolved" population.
MIN_CELLS_PER_VESSEL = 2


# ---------------------------------------------------------------------------
# Length scales
# ---------------------------------------------------------------------------
def diffusion_length_um(
    hypoxia_onset_um: float = HYPOXIA_ONSET_UM,
    hypoxia_threshold: float = OxygenConfig().hypoxia_threshold,
) -> float:
    """The surrogate's length scale ``L``, derived (not fitted) from two anchors.

    Solves ``exp(-hypoxia_onset_um / L) == hypoxia_threshold`` for ``L``::

        L = hypoxia_onset_um / ln(1 / hypoxia_threshold)

    so that tissue at the published oxygen-diffusion distance lands exactly on the
    hypoxia threshold the rest of the codebase already bands at. At the defaults
    (150 um, 0.2) this is 93.2 um.

    Raises ``ValueError`` if the threshold is not strictly inside ``(0, 1)`` — the
    equation has no positive solution there.
    """
    t = float(hypoxia_threshold)
    d = float(hypoxia_onset_um)
    if not (0.0 < t < 1.0):
        raise ValueError(
            f"hypoxia_threshold must be in (0, 1) to derive a length scale; got {t}"
        )
    if not d > 0.0:
        raise ValueError(f"hypoxia_onset_um must be > 0; got {d}")
    return d / float(np.log(1.0 / t))


def pde_diffusion_length_um(
    micro: Microenvironment, cfg: OxygenConfig = OxygenConfig()
) -> Dict[str, float]:
    """The PDE's OWN effective diffusion length, for an order-of-magnitude check.

    The solved balance is ``D*lap(O) + supply*vasc*(1-O) - consumption*rho_tum*O = 0``,
    whose homogeneous decay length away from a source is ``sqrt(D / k)`` with
    ``k = consumption * rho_tum`` (the areal tumour density of
    :func:`physics.oxygen._tumor_density`, in cells per 144 um^2 reference tile).

    Reported in ``results/oxygen_demotion.md`` as a DIAGNOSTIC, not a calibration.
    On the real sections it comes out at ~310-800 um against the ~150 um
    Thomlinson & Gray oxygen-diffusion distance, so the solved field is smoothed
    at a length scale 2-5x longer than tissue oxygen is measured to vary over.
    That is the mechanism behind "the diffusion step destroys information": it
    low-pass-filters its own source. Neither number is anchored to a measured
    oxygen tension — ``OxygenConfig.consumption`` has never been fitted — so this
    is a statement about the model's internal scale, not about the tumour.

    Returns ``{'L_mean_um', 'L_median_um', 'mean_rho_tum', 'median_rho_tum'}``,
    evaluated over the tissue mask. Entries are ``inf`` where the tumour density
    is zero (no sink -> no decay).
    """
    from .oxygen import _tumor_density  # local: avoid a circular import at module load

    rho = np.asarray(_tumor_density(micro), dtype=np.float64)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    vals = rho[tissue] if tissue.any() else rho.ravel()
    D = float(cfg.diffusivity)
    k = float(cfg.consumption)

    def _L(r: float) -> float:
        return float(np.sqrt(D / (k * r))) if (k > 0.0 and r > 0.0) else float("inf")

    mean_rho = float(vals.mean()) if vals.size else 0.0
    med_rho = float(np.median(vals)) if vals.size else 0.0
    return {
        "L_mean_um": _L(mean_rho),
        "L_median_um": _L(med_rho),
        "mean_rho_tum": mean_rho,
        "median_rho_tum": med_rho,
        "note": (  # type: ignore[dict-item]
            "sqrt(D / (consumption * rho_tum)); consumption is UNCALIBRATED "
            "against any measured oxygen tension, so this is an order-of-magnitude "
            "consistency check only"
        ),
    }


# ---------------------------------------------------------------------------
# The distance transform
# ---------------------------------------------------------------------------
def _vessel_points(
    cells: CellPopulation,
    domain: Domain,
    source: str,
    link_distance_um: Optional[float],
    min_cells_per_vessel: int,
) -> Tuple[np.ndarray, Dict[str, object]]:
    """Vessel point cloud (m, 2) in um plus a provenance block. See :data:`VESSEL_SOURCES`."""
    s = str(source).strip().lower()
    endo = np.asarray(cells.cell_type) == int(CellType.ENDOTHELIAL)
    xy_endo = np.column_stack([
        np.asarray(cells.x, dtype=np.float64)[endo],
        np.asarray(cells.y, dtype=np.float64)[endo],
    ])

    if s == "endothelial":
        return xy_endo, {
            "source": "endothelial",
            "what": "every measured CellType.ENDOTHELIAL cell centroid",
            "n_points": int(xy_endo.shape[0]),
            "n_endothelial_cells": int(xy_endo.shape[0]),
            "n_vessels": None,
            "link_distance_um": None,
            "min_cells_per_vessel": None,
        }

    if s == "vessel_network":
        # Lazy import: keeps physics free of a hard dependency on reconstruct.
        from ..reconstruct.vessel_network import LINK_DISTANCE_UM, segment_vessels

        link = float(LINK_DISTANCE_UM if link_distance_um is None else link_distance_um)
        vessels = segment_vessels(cells, domain, link_distance_um=link)
        xy = np.asarray(vessels["xy"], dtype=np.float64)
        labels = np.asarray(vessels["labels"], dtype=int)
        counts = np.asarray(vessels["cells_per_vessel"], dtype=int)
        keep_vessel = counts >= int(min_cells_per_vessel)
        keep = keep_vessel[labels] if labels.size else np.zeros(0, dtype=bool)
        return xy[keep], {
            "source": "vessel_network",
            "what": (
                "endothelial cells clustered into DISCRETE vessel objects by "
                "reconstruct.vessel_network.segment_vessels, keeping objects with "
                f">= {int(min_cells_per_vessel)} cells"
            ),
            "n_points": int(keep.sum()),
            "n_endothelial_cells": int(xy.shape[0]),
            "n_vessels": int(keep_vessel.sum()),
            "n_vessels_total": int(counts.size),
            "link_distance_um": link,
            "min_cells_per_vessel": int(min_cells_per_vessel),
        }

    raise ValueError(f"unknown vessel source {source!r}; expected one of {VESSEL_SOURCES}")


def vessel_distance_um(
    micro: Microenvironment,
    source: str = "endothelial",
    link_distance_um: Optional[float] = None,
    min_cells_per_vessel: int = MIN_CELLS_PER_VESSEL,
    at: str = "grid",
    detail: bool = False,
):
    """Euclidean distance (um) to the nearest perfused vessel — the raw geometry.

    Exposed separately from :func:`oxygen_from_vessel_distance` because the
    distance is the object that carries the predictive signal; the [0, 1] map on
    top of it is a monotone re-labelling that cannot change any rank statistic.
    A caller who wants a different mapping should take this array and apply it.

    Parameters
    ----------
    micro : Microenvironment
        Supplies ``domain`` (grid + dx) and ``cells``.
    source : {'endothelial', 'vessel_network'}
        See :data:`VESSEL_SOURCES`. Default ``'endothelial'`` is the predictor
        actually benchmarked in ``results/state_validation.md``.
    link_distance_um : float, optional
        Single-linkage threshold forwarded to
        :func:`reconstruct.vessel_network.segment_vessels`; ignored unless
        ``source='vessel_network'`` (default: that module's 20 um).
    min_cells_per_vessel : int
        Discard segmented objects smaller than this; ignored unless
        ``source='vessel_network'``.
    at : {'grid', 'cells'}
        ``'grid'`` -> distance evaluated at every grid-cell CENTRE, shape
        ``domain.shape``, matching ``Domain.to_index``'s ``(col + 0.5) * dx``
        convention. ``'cells'`` -> distance at each cell's exact measured
        coordinate, shape ``(n_cells,)``; use this for a per-cell benchmark, where
        binning to the mesh would only add discretisation noise.
    detail : bool
        Also return the provenance block.

    Returns
    -------
    np.ndarray, or (np.ndarray, dict) when ``detail``
        Distances in um. If the section contains NO vessel points at all, every
        entry is ``+inf`` (there is no vessel to be near), which
        :func:`oxygen_from_vessel_distance` maps to an oxygen surrogate of 0.
    """
    domain = micro.domain
    xy, info = _vessel_points(
        micro.cells, domain, source, link_distance_um, int(min_cells_per_vessel)
    )

    mode = str(at).strip().lower()
    if mode == "grid":
        ny, nx = domain.shape
        dx = float(domain.dx)
        gx = (np.arange(nx, dtype=np.float64) + 0.5) * dx
        gy = (np.arange(ny, dtype=np.float64) + 0.5) * dx
        GX, GY = np.meshgrid(gx, gy)
        query = np.column_stack([GX.ravel(), GY.ravel()])
        out_shape: Tuple[int, ...] = (ny, nx)
    elif mode == "cells":
        query = np.column_stack([
            np.asarray(micro.cells.x, dtype=np.float64),
            np.asarray(micro.cells.y, dtype=np.float64),
        ])
        out_shape = (int(micro.cells.n),)
    else:
        raise ValueError(f"unknown 'at' {at!r}; expected 'grid' or 'cells'")

    if xy.shape[0] == 0:
        dist = np.full(out_shape, np.inf, dtype=np.float64)
    else:
        d, _ = cKDTree(xy).query(query, k=1)
        dist = np.asarray(d, dtype=np.float64).reshape(out_shape)

    info = dict(info)
    info["at"] = mode
    info["median_distance_um"] = (
        float(np.median(dist[np.isfinite(dist)])) if np.isfinite(dist).any()
        else float("inf")
    )
    return (dist, info) if detail else dist


def oxygen_from_vessel_distance(
    micro: Microenvironment,
    cfg: OxygenConfig = OxygenConfig(),
    source: str = "endothelial",
    diffusion_length_um_: Optional[float] = None,
    hypoxia_onset_um: float = HYPOXIA_ONSET_UM,
    link_distance_um: Optional[float] = None,
    min_cells_per_vessel: int = MIN_CELLS_PER_VESSEL,
    detail: bool = False,
):
    """Oxygen SURROGATE field (ny, nx) in [0, 1] from distance to the nearest vessel.

    **This is not a solved oxygen tension and must never be labelled one.** It is a
    dimensionless, monotone, published-distance-calibrated proxy that out-predicts
    the reaction-diffusion solve against measured hypoxia markers by 33%
    (``results/state_validation.md`` §1.2; ``results/oxygen_demotion.md``).

    The map is ``O(d) = exp(-d / L)``, clipped to [0, 1], where ``d`` is the
    distance from each grid-cell centre to the nearest vessel point and ``L`` comes
    from :func:`diffusion_length_um` — i.e. it is chosen so that tissue at the
    published oxygen-diffusion distance ``hypoxia_onset_um`` sits exactly on
    ``cfg.hypoxia_threshold``. Consequences, all deliberate:

    * ``O = 1`` exactly AT a vessel point and decays smoothly away from it;
    * ``O = cfg.hypoxia_threshold`` at ``hypoxia_onset_um`` (150 um by default), so
      ``cells.agents.update_cell_states`` bands HYPOXIC beyond that distance and
      VIABLE inside ``L * ln(1 / (2 * threshold))`` (85.9 um at the defaults) with
      no new tuning knob anywhere downstream;
    * the field is NOT masked to tissue. The mask is a bookkeeping convention that
      gates the PDE's supply and therefore its absolute level
      (``results/state_validation.md`` §2.2); the surrogate has no supply term to
      gate, and every measured cell lies inside the mask anyway, so masking here
      would only introduce a discontinuity with no readout attached to it.

    Parameters
    ----------
    micro : Microenvironment
        Grid + measured cells. Not modified; the caller stores the result.
    cfg : OxygenConfig
        Only ``hypoxia_threshold`` is read, and only to derive ``L``. No
        diffusivity / supply / consumption is used — there is no solve here.
    source : {'endothelial', 'vessel_network'}
        See :data:`VESSEL_SOURCES`.
    diffusion_length_um_ : float, optional
        Override ``L`` (um) directly. Trailing underscore because
        :func:`diffusion_length_um` is a function in this namespace. Overriding
        ``L`` moves the band edges; it CANNOT change any rank correlation, because
        the map stays monotone.
    hypoxia_onset_um : float
        The published anchor distance used to derive ``L``. Ignored when
        ``diffusion_length_um_`` is given.
    link_distance_um, min_cells_per_vessel
        Forwarded to :func:`vessel_distance_um`; only used by
        ``source='vessel_network'``.
    detail : bool
        Also return a provenance dict (point source, counts, ``L``, band-edge
        distances, resulting band fractions over tissue).

    Returns
    -------
    np.ndarray, or (np.ndarray, dict) when ``detail``
        (ny, nx) surrogate in [0, 1]. All-zero if the section has no vessel points
        (no supply anywhere -> maximally hypoxic everywhere), which is the
        physically honest answer rather than a silent fallback.
    """
    dist, info = vessel_distance_um(
        micro, source=source, link_distance_um=link_distance_um,
        min_cells_per_vessel=int(min_cells_per_vessel), at="grid", detail=True,
    )

    t = float(cfg.hypoxia_threshold)
    L = (float(diffusion_length_um_) if diffusion_length_um_ is not None
         else diffusion_length_um(hypoxia_onset_um, t))
    if not L > 0.0:
        raise ValueError(f"diffusion length must be > 0 um; got {L}")

    with np.errstate(over="ignore"):
        field = np.exp(-dist / L)
    field = np.clip(np.nan_to_num(field, nan=0.0, posinf=0.0, neginf=0.0), 0.0, 1.0)

    if not detail:
        return field

    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    sel = field[tissue] if tissue.any() else field.ravel()
    hypoxic_edge = L * float(np.log(1.0 / t)) if 0.0 < t < 1.0 else float("nan")
    viable_edge = (L * float(np.log(1.0 / (2.0 * t)))
                   if 0.0 < 2.0 * t < 1.0 else float("nan"))
    info = dict(info)
    info.update({
        "kind": "SURROGATE (geometric proxy), NOT a solved oxygen tension",
        "mapping": "O(d) = exp(-d / L), clipped to [0, 1]",
        "diffusion_length_um": float(L),
        "hypoxia_onset_um": float(hypoxia_onset_um),
        "hypoxia_threshold": t,
        "hypoxic_beyond_um": hypoxic_edge,
        "viable_within_um": viable_edge,
        "calibration": (
            "L = hypoxia_onset_um / ln(1/hypoxia_threshold); hypoxia_onset_um is "
            "the Thomlinson & Gray 1955 oxygen-diffusion distance (PMID 13304213), "
            "NOT a fit to the marker data. Any strictly monotone map of distance "
            "gives an identical Spearman rho, so this constant sets the band edges "
            "and nothing else."
        ),
        "mean_over_tissue": float(sel.mean()) if sel.size else 0.0,
        "frac_below_threshold_over_tissue": (
            float((sel < t).mean()) if sel.size else 0.0
        ),
    })
    return field, info


# ---------------------------------------------------------------------------
# Hypoxia read directly off measured transcripts
# ---------------------------------------------------------------------------
def _counts_column(expression, genes: np.ndarray, gene: str) -> Optional[np.ndarray]:
    """Raw per-cell counts for one panel gene, or None if the panel lacks it."""
    genes = np.asarray(genes, dtype=str)
    hit = np.flatnonzero(genes == gene)
    if hit.size == 0:
        return None
    j = int(hit[0])
    if _sp.issparse(expression):
        return np.asarray(expression[:, [j]].todense()).ravel().astype(np.float64)
    return np.asarray(expression)[:, j].astype(np.float64)


def _depth_normalised_log1p(counts: np.ndarray, totals: np.ndarray) -> np.ndarray:
    """Counts-per-median-cell normalise then log1p — the repo's standard depth
    correction (identical estimator to
    ``validation.cellstate_validation.normalized_log1p``, reproduced here only so
    ``physics`` takes no dependency on ``validation``)."""
    counts = np.asarray(counts, dtype=np.float64)
    totals = np.asarray(totals, dtype=np.float64)
    scale = float(np.median(totals[totals > 0])) if np.any(totals > 0) else 1.0
    safe = np.where(totals > 0, totals, 1.0)
    return np.log1p(counts / safe * scale)


def _zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    s = float(v.std())
    return (v - float(v.mean())) / s if s > 0 else np.zeros_like(v)


def _nadaraya_watson(
    domain: Domain, x: np.ndarray, y: np.ndarray, values: np.ndarray,
    bandwidth_um: float,
) -> np.ndarray:
    """Smoothed per-cell MEAN onto the grid (not a density).

    ``gaussian_filter(sum of values) / gaussian_filter(count)`` — the Nadaraya-
    Watson kernel-regression estimator. This is deliberately NOT the KDE used for
    the vascular field: a density would confound "these cells express a lot of
    VEGFA" with "there are a lot of cells here", and local cell density is itself
    one of the trivial baselines the PDE was benchmarked against. Grid cells with
    no smoothed support fall back to the section mean, which is the least
    informative value available and therefore the honest filler.
    """
    num = domain.zeros()
    den = domain.zeros()
    col, row = domain.to_index(np.asarray(x, float), np.asarray(y, float))
    np.add.at(num, (row, col), np.asarray(values, dtype=np.float64))
    np.add.at(den, (row, col), 1.0)
    sigma = max(float(bandwidth_um) / float(domain.dx), 1e-6)
    num = gaussian_filter(num, sigma=sigma, mode="nearest")
    den = gaussian_filter(den, sigma=sigma, mode="nearest")
    fallback = float(np.mean(values)) if np.size(values) else 0.0
    out = np.full(domain.shape, fallback, dtype=np.float64)
    ok = den > 1e-12
    out[ok] = num[ok] / den[ok]
    return out


#: Kernel bandwidth (um) for smoothing the per-cell marker score onto the grid.
#: Matches ``reconstruct.real_fields._DEFAULT_BANDWIDTH_UM`` so the marker field
#: and the vascular field are smoothed at the same physical scale and any
#: difference between them is signal rather than bandwidth.
MARKER_BANDWIDTH_UM = 40.0


def hypoxia_from_markers(
    micro: Microenvironment,
    expression,
    genes: Sequence[str],
    gene_set: Sequence[str] = HYPOXIA_MARKER_GENES,
    bandwidth_um: float = MARKER_BANDWIDTH_UM,
    as_oxygen: bool = True,
    detail: bool = False,
) -> Union[None, np.ndarray, Dict[str, object]]:
    """Hypoxia read DIRECTLY off measured hypoxia-response transcripts, or ``None``.

    The only route in this module that is a MEASUREMENT rather than a model. Each
    panel-present gene of ``gene_set`` is depth-normalised (counts-per-median-cell,
    ``log1p``) and z-scored; the per-cell score is their mean; the score is
    smoothed onto the grid with :func:`_nadaraya_watson` (a kernel-regressed MEAN,
    not a density, so it does not re-encode local cell crowding).

    With ``as_oxygen=True`` (default) the returned field is sign-aligned to OXYGEN
    — larger means better perfused — by mapping the smoothed score through its own
    empirical rank to [0, 1] and inverting: ``O = 1 - rank_pct(score)``. The rank
    map is used rather than a z-score-to-[0,1] squash because the marker score has
    no natural scale and a rank map is the only monotone normalisation that cannot
    be gamed by an outlier. It is a UNIFORM field by construction, so its ABSOLUTE
    level is meaningless and it must not be compared to the PDE's mean oxygen; only
    its spatial ORDERING is measured. With ``as_oxygen=False`` the raw smoothed
    hypoxia-response score is returned instead (larger = MORE hypoxic).

    **HONEST CAVEATS.**

    1. VEGFA and its co-markers are hypoxia-**RESPONSE** genes, i.e. the output of
       a transcriptional programme with its own kinetics and its own non-hypoxic
       inducers — not an oxygen tension. A field built from them is a measured
       *hypoxia response*, and calling it "oxygen" is a convenience of sign, not a
       claim about pO2.
    2. Targeted in-situ panels drop out. A zero is "no probe fired", not "no
       transcript" (VEGFA fires in 63.4% of cells on the IO panel).
    3. Feeding this back into a state model that is then validated against the SAME
       markers would be circular. Use it as the hypoxia covariate OR as the
       validation target, never both in one claim.

    Parameters
    ----------
    micro : Microenvironment
        Supplies the grid and the cell coordinates. ``expression`` rows must be
        index-aligned to ``micro.cells`` — asserted, never assumed.
    expression : scipy.sparse matrix or ndarray, shape (n_cells, n_genes)
        Raw per-cell counts, as returned by ``data.real.singlecell_spatial.load``.
    genes : sequence of str, length n_genes
        Panel gene symbols, aligned to ``expression``'s columns.
    gene_set : sequence of str
        Hypoxia-response genes to look for, in preference order. Default
        :data:`HYPOXIA_MARKER_GENES` (HIF1A deliberately excluded — see its
        docstring).
    bandwidth_um : float
        Gaussian bandwidth for the grid smoother.
    as_oxygen : bool
        Sign convention of the returned field (see above).
    detail : bool
        Return the full dict instead of the bare field / ``None``.

    Returns
    -------
    None
        When the panel carries NO gene of ``gene_set``. This is the expected result
        on the 319-gene Xenium Human Brain panel, which has no VEGFA, CA9, ADM,
        NDRG1 or SLC2A1: on that section hypoxia is **UNMEASURED**, and returning
        ``None`` says so instead of manufacturing a field. With ``detail=True`` a
        dict is returned whose ``'field'`` is ``None`` and whose ``'reason'`` names
        every gene that was looked for and not found.
    np.ndarray
        Otherwise, an (ny, nx) field. With ``detail=True``, a dict with ``field``,
        ``per_cell`` (the un-smoothed per-cell score), ``genes_used``,
        ``genes_missing``, ``detection`` (per-gene nonzero fraction), ``reason``
        (``None`` on success) and ``kind``.

    Raises
    ------
    ValueError
        If ``expression`` has a different number of rows than ``micro.cells`` — a
        misalignment would silently attach one cell's transcripts to another's
        coordinates, so it is a hard error rather than a warning.
    """
    genes_arr = np.asarray(genes, dtype=str)
    n_cells = int(micro.cells.n)
    n_rows = int(expression.shape[0])
    if n_rows != n_cells:
        raise ValueError(
            f"expression rows ({n_rows}) != cells ({n_cells}); the reconstruction "
            "and the expression matrix are not aligned"
        )

    if _sp.issparse(expression):
        totals = np.asarray(expression.sum(axis=1)).ravel().astype(np.float64)
    else:
        totals = np.asarray(expression, dtype=np.float64).sum(axis=1)

    used: List[str] = []
    missing: List[str] = []
    per_gene: Dict[str, np.ndarray] = {}
    detection: Dict[str, float] = {}
    for g in gene_set:
        raw = _counts_column(expression, genes_arr, g)
        if raw is None:
            missing.append(str(g))
            continue
        used.append(str(g))
        detection[str(g)] = float((raw > 0).mean())
        per_gene[str(g)] = _zscore(_depth_normalised_log1p(raw, totals))

    if not used:
        reason = (
            "no hypoxia-response gene of "
            f"{list(gene_set)} is on this {genes_arr.size}-gene panel "
            f"(all absent: {missing}). Hypoxia is UNMEASURED on this section — it "
            "has not been simulated away and must not be reported as if it had. "
            "Use oxygen_mode='vessel_distance' for a geometric surrogate, or state "
            "that hypoxia is unmeasured."
        )
        if detail:
            return {
                "field": None, "per_cell": None, "genes_used": [],
                "genes_missing": missing, "detection": {}, "reason": reason,
                "kind": "UNAVAILABLE",
            }
        return None

    score = np.mean(np.vstack([per_gene[g] for g in used]), axis=0)
    smoothed = _nadaraya_watson(
        micro.domain, micro.cells.x, micro.cells.y, score, float(bandwidth_um)
    )

    if as_oxygen:
        flat = smoothed.ravel()
        order = np.argsort(np.argsort(flat, kind="stable"), kind="stable")
        pct = order.astype(np.float64) / max(flat.size - 1, 1)
        field = (1.0 - pct).reshape(smoothed.shape)
    else:
        field = smoothed

    if not detail:
        return field
    return {
        "field": field,
        "per_cell": score,
        "genes_used": used,
        "genes_missing": missing,
        "detection": detection,
        "reason": None,
        "as_oxygen": bool(as_oxygen),
        "bandwidth_um": float(bandwidth_um),
        "kind": (
            "MEASURED hypoxia-RESPONSE (not an oxygen tension). Sign-flipped and "
            "rank-normalised to an oxygen-like [0, 1] field; its ABSOLUTE level is "
            "meaningless by construction (uniform by rank) — only its spatial "
            "ordering is measured."
            if as_oxygen else
            "MEASURED hypoxia-RESPONSE score (larger = more hypoxic); not an "
            "oxygen tension."
        ),
    }


if __name__ == "__main__":  # pragma: no cover - exercised on the real cached section
    import time

    from ..reconstruct.real_pipeline import reconstruct_real

    print("=== oxygen_proxy self-test on the REAL Xenium GBM section ===")
    L = diffusion_length_um()
    print(f"derived diffusion length L = {L:.2f} um "
          f"(onset {HYPOXIA_ONSET_UM:.0f} um, threshold "
          f"{OxygenConfig().hypoxia_threshold})")
    assert abs(np.exp(-HYPOXIA_ONSET_UM / L) - OxygenConfig().hypoxia_threshold) < 1e-12

    t0 = time.time()
    micro = reconstruct_real("xenium", dx_um=30.0)
    print(f"reconstructed {micro.cells.n:,} cells, grid {micro.domain.shape} "
          f"[{time.time() - t0:.1f}s]")

    t0 = time.time()
    field, info = oxygen_from_vessel_distance(micro, detail=True)
    dt_proxy = time.time() - t0
    print(f"\nvessel-distance surrogate [{dt_proxy:.3f}s]")
    for k in ("n_points", "diffusion_length_um", "hypoxic_beyond_um",
              "viable_within_um", "median_distance_um", "mean_over_tissue",
              "frac_below_threshold_over_tissue"):
        print(f"  {k:36s} {info[k]}")
    assert field.shape == micro.domain.shape
    assert float(field.min()) >= 0.0 and float(field.max()) <= 1.0

    # Monotonicity: the surrogate must decrease with distance, exactly.
    d = vessel_distance_um(micro)
    o = np.exp(-d / info["diffusion_length_um"])
    assert np.allclose(o, field), "mapping is not exp(-d/L)"
    order = np.argsort(d.ravel())
    assert np.all(np.diff(field.ravel()[order]) <= 1e-12), "surrogate not monotone"

    # The discrete-vessel variant runs and keeps fewer points.
    f2, i2 = oxygen_from_vessel_distance(micro, source="vessel_network", detail=True)
    print(f"\nvessel_network source: {i2['n_points']:,} points from "
          f"{i2['n_vessels']:,}/{i2['n_vessels_total']:,} objects "
          f"(link {i2['link_distance_um']:.0f} um)")
    assert i2["n_points"] <= info["n_points"]

    # The PDE's own length scale, for the order-of-magnitude check.
    print("\nPDE own diffusion length:", pde_diffusion_length_um(micro))

    print("\noxygen_proxy self-test PASSED")
