"""SUB-GRID vessel footprints with a DATA-DERIVED radius — the geometry fix.

Why this module exists
----------------------
``reconstruct.vessel_network.discrete_vascular_field`` paints, around **every
endothelial cell**, a binary disc of a single hardcoded radius
(``VESSEL_RADIUS_UM = 7.5 µm``) and calls a grid cell "vessel" when its *centre*
falls inside that disc. ``segmentation_report`` then divides the painted pixel
count by the tissue pixel count and compares the result to the measured GBM CD34
vascular area fraction of 2.1 ± 0.3 % (Angelucci/Tamma C et al., Cells 2022;12(1):11,
PMID 36611806). Two things are wrong with that quantity, and one thing that was
*believed* wrong turns out not to be:

**NOT wrong (measured, see results/vessel_geometry_fix.md §1).** The binary
centre-in-disc rule was suspected of making the area fraction an artefact of the
mesh, on the grounds that at ``dx = 12 µm`` a single grid cell already has an
equivalent-disc radius of 6.8 µm, so a 7.5 µm vessel "cannot be represented".
That reasoning is wrong for an *aggregate* quantity. Point-sampling a fixed
physical region on a regular lattice is an **unbiased** estimator of that
region's area at any pitch: the expected number of lattice centres inside a
region of area ``A`` is ``A/dx²`` whatever ``dx`` is. Sweeping the production
sections over ``dx = 30/20/12/8/6 µm`` moves the area fraction by at most 1.03×,
i.e. it was already inside the project's 1.05× grid-independence convention.
What *is* grid-dependent under the binary rule is **which vessels survive**: at
``dx = 30 µm`` the mean painted-pixel count is 0.26 per vessel object, so ~three
quarters of the segmented vessels contribute no source at all and the field's
*geometry* — not its integral — is a function of the mesh.

**Wrong #1 — the radius is a constant, and the constant was chosen partly from
the mesh.** ``VESSEL_RADIUS_UM``'s own docstring says it is "the smallest radius
that reliably paints ≥1 grid cell per vessel while staying inside the published
range". A parameter justified by the grid cannot then be used to produce a
tissue measurement. :func:`vessel_geometry` replaces it with a per-object radius
derived from the section's own endothelial content through two **measured human**
conversions (below).

**Wrong #2 — the footprint is painted per endothelial CELL, not per vessel.**
This is the defect that actually moves the number. The painted area is
essentially ``n_EC × π r²``, so the reported area fraction tracks the
**endothelial-cell density**, which is a property of the marker panel and the
typing threshold, not of the vasculature. On the two production sections the
vessel densities are nearly identical (111 vs 104 objects/mm²) while the
endothelial-cell densities differ by 2.7× (143 vs 382 EC/mm²) — and the reported
area fractions differed by exactly that ratio (2.50 % vs 5.79 %). A vessel
profile lined by ``k`` endothelial cells is **one** vessel, not ``k`` of them.

Scope — "Wrong #2" is an argument about an AREAL FRACTION, and it does NOT
transfer to the delivery source
---------------------------------------------------------------------------
Measured in ``results/delivery_source_geometry.md``; recorded here because the
gap between the corrected geometry and the field that actually delivers drug
looks like an oversight and is not one.

The engine's delivery term is ``vascular_density * bbb_penetrance * plasma``, a
**transvascular flux**, whose geometric weight is the vessel **WALL** area per
unit tissue (``S/V``). This module paints a vessel **CROSS-SECTION** area. Those
are different powers of the radius, and with this module's own ``r = r1 * k``
they are different powers of ``k``::

    wall area          ~  2 pi r L      ~  k        (what a flux BC wants)
    cross-section area ~  pi r^2 + 2rL  ~  k^2      (what this module paints)

Endothelial-cell COUNT is a direct wall-area proxy in both regimes, from the two
measured human constants below: axially, one cell lines
:data:`EC_AXIAL_PITCH_UM` = 50 um of wall (PMID 38723085); circumferentially,
``2 pi r = k w`` makes the perimeter proportional to ``k``. So a source that
deposits a CONSTANT per endothelial cell implements the wall-area law, and the
legacy binary field in ``vessel_network`` does exactly that.

Measured consequence of routing the delivery source through :func:`vessel_geometry`
instead: the source deposited per endothelial cell stops being constant and
spans **8.06x / 8.40x / 8.27x / 8.43x** across ``k`` on the four cached real
sections (30.4 um^2/EC at ``k = 1``, ~250 um^2/EC at ``k >= 5``). On the module's
own fidelity score, :func:`vascular_heterogeneity`, the legacy field scores
0.94-1.01 against a target of 1.0 and the sub-grid field 1.45-2.26 — i.e. the
``r ∝ k`` over-dispersion documented below reappears as sharpened architecture.
The 7.5 um radius of the legacy field is arbitrary, but it is a **gauge**
constant that ``PS`` absorbs exactly (``results/transvascular_exchange_revalidation.md``
section 5); the ``k^2`` scaling is not, because it re-weights vessels relative to
one another inside one section.

**So the two representations answer two different questions.** A stereological
areal fraction counts each vessel profile once (Marien's rule: a grid point is
vascular when it falls on "an endothelial cell or a vessel lumen"), which is what
:func:`vascular_area_fraction` computes and what "Wrong #2" correctly fixes. A
leak source is a wall-area quantity, and for it the per-endothelial-cell painting
is the RIGHT law with an arbitrary constant. The default of
``discrete_vascular_field`` is therefore left at ``geometry='binary'``
deliberately, not by omission.

What this module does
---------------------
1. :func:`vessel_geometry` — a per-vessel-object radius and footprint from the
   DATA (endothelial cells per object + the object's own spatial extent) via two
   measured conversions, with no constant chosen to hit any target.
2. :func:`subgrid_coverage_field` — paints each footprint as a **fractional area
   coverage** on a fixed *physical* sub-lattice, so a vessel far smaller than a
   grid cell still contributes exactly its own area. The painted total is
   computed with no reference to ``dx`` at all and is therefore grid-independent
   by construction, not by luck.
3. :func:`vascular_area_fraction` — the resulting area fraction, which is now a
   property of the tissue.

Nothing in ``core`` / ``physics`` / ``sim.engine`` / ``molecular`` is touched, and
the legacy binary field in ``vessel_network`` is left bit-identical so that no
already-published result silently moves. The corrected geometry is opt-in at the
field level (``discrete_vascular_field(..., geometry='subgrid')``) and is the
default for the *reported* area fraction, which is the quantity being fixed.

Parameter provenance — every constant is measured, none is fitted
-----------------------------------------------------------------
``EC_PER_MM_VESSEL = 20.0`` — endothelial cells per millimetre of vasculature.
  Measured directly in **human** cortex: the H01 petavoxel reconstruction of
  1.05 mm³ of human temporal cortex contains ~230 mm of blood vessel lined by
  **4,604 endothelial cells, "~20 per mm of vasculature"** — Shapson-Coe A et al.,
  Science 2024;384(6696):eadk4858, PMID 38723085, doi:10.1126/science.adk4858.
  Its reciprocal, :data:`EC_AXIAL_PITCH_UM` = 50 µm, is the length of vessel that
  one endothelial cell lines when the wall is one cell thick.

``VESSEL_LENGTH_DENSITY_MM_PER_MM3 = 230/1.05 = 219`` and
``CORTEX_VASCULAR_VOLUME_FRACTION = 0.015`` — from the *same* human dataset
  (PMID 38723085): 230 mm of vessel in 1.05 mm³, and blood vessels occupying
  ~1.5 % of neuropil volume. Using one dataset for both keeps them mutually
  consistent. Their ratio fixes the **mean vessel radius** with no free
  parameter, because a volume fraction is a length density times a
  cross-sectional area::

      VV = L_V · π r̄²   ⟹   r̄ = sqrt(VV / (π L_V)) = 4.67 µm

  (The H01 vessels were fixed in a collapsed state, so 4.67 µm is a **lower**
  bound on the perfused radius; see the sensitivity in
  ``results/vessel_geometry_fix.md``.)

``MEAN_EC_PER_CROSS_SECTION = 1.5`` — endothelial cells per vessel cross-section,
  averaged over the microvasculature. Measured by volume EM: brain capillaries
  "are composed of either one or two endothelial cells, with roughly equal
  proportions" (hence 1.5), while transitional vessels near the penetrating
  arteriole and ascending venule carry two to six — Sargent SM et al., Vasc Biol
  2023;5(1):e230010, PMID 37582180, doi:10.1530/VB-23-0010. (Mouse cortex,
  MICrONS volume EM.)

``SINGLE_EC_VESSEL_RADIUS_UM`` — **derived, not chosen**: the ``k`` endothelial
  cells of a cross-section tile its circumference, so ``2πr = k·w`` and
  ``r ∝ k``. Dividing the measured mean radius by the measured mean
  cells-per-cross-section gives the radius of a vessel lined by a single
  endothelial cell::

      r₁ = r̄ / k̄ = 4.669 / 1.5 = 3.11 µm

  Cross-check, independent of the arithmetic above: Sargent 2023 measures a
  one-endothelial-cell capillary **lumen** diameter of ~2–5 µm, i.e. a lumen
  radius of 1.0–2.5 µm plus wall — the same scale. Second cross-check, from the
  sections themselves: a two-endothelial-cell profile is two cells interlocked
  around one lumen (Sargent 2023), so their nuclear centroids sit across a
  diameter; the measured median separation of two-cell objects is 11.3 µm
  (xenium_io) / 12.7 µm (GSE324268), implying ``r(k=2) ≈ 5.7–6.4 µm`` against
  this model's ``2 r₁ = 6.2 µm``.

  **Documented limitation.** Sargent 2023 also reports that, *within* the
  capillary zone, two-endothelial-cell capillaries are only ~0.4 µm wider than
  one-cell capillaries — much less than the factor 2 that strict circumferential
  tiling predicts — and states explicitly that endothelial number "alone is
  insufficient to explain the full range of capillary diameters". The ``r ∝ k``
  model therefore over-disperses radii at the capillary end. ``radius_mode=
  'stereological'`` gives the opposite extreme (every object at the measured mean
  r̄ = 4.67 µm, no k dependence); both are reported in the results file so the
  reader can see how much of the answer rides on the scaling.

Run:  ./.venv/bin/python -m tmesim.reconstruct.vessel_geometry
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from ..core import Domain

__all__ = [
    "EC_PER_MM_VESSEL",
    "EC_AXIAL_PITCH_UM",
    "VESSEL_LENGTH_DENSITY_MM_PER_MM3",
    "CORTEX_VASCULAR_VOLUME_FRACTION",
    "MEAN_EC_PER_CROSS_SECTION",
    "MEAN_VESSEL_RADIUS_UM",
    "SINGLE_EC_VESSEL_RADIUS_UM",
    "RADIUS_MODES",
    "vessel_geometry",
    "subgrid_coverage_field",
    "vascular_area_fraction",
    "vascular_heterogeneity",
]

# --------------------------------------------------------------------------- #
# Measured constants (see the module docstring for the full provenance)
# --------------------------------------------------------------------------- #
#: Endothelial cells per mm of vasculature — HUMAN cortex, PMID 38723085.
EC_PER_MM_VESSEL = 20.0

#: Length of vessel lined by one endothelial cell (µm) = 1000 / EC_PER_MM_VESSEL.
EC_AXIAL_PITCH_UM = 1000.0 / EC_PER_MM_VESSEL

#: Vessel length density of human cortex (mm of vessel per mm³) — 230 mm in
#: 1.05 mm³, PMID 38723085.
VESSEL_LENGTH_DENSITY_MM_PER_MM3 = 230.0 / 1.05

#: Vascular volume fraction of human cortex — PMID 38723085 ("collapsed blood
#: vessels ~1.5 %" of neuropil volume; a LOWER bound on the perfused fraction).
CORTEX_VASCULAR_VOLUME_FRACTION = 0.015

#: Endothelial cells per vessel cross-section, microvascular mean — PMID 37582180.
MEAN_EC_PER_CROSS_SECTION = 1.5

#: Mean vessel radius (µm), DERIVED from the two human-cortex measurements above
#: via VV = L_V · π r̄². Nothing here is fitted.
MEAN_VESSEL_RADIUS_UM = 1e3 * float(np.sqrt(
    CORTEX_VASCULAR_VOLUME_FRACTION
    / (np.pi * VESSEL_LENGTH_DENSITY_MM_PER_MM3)
))  # 4.669 µm

#: Radius (µm) of a vessel whose cross-section is lined by ONE endothelial cell.
#: DERIVED: circumferential tiling gives r ∝ k, so r₁ = r̄ / k̄.
SINGLE_EC_VESSEL_RADIUS_UM = MEAN_VESSEL_RADIUS_UM / MEAN_EC_PER_CROSS_SECTION

#: How :func:`vessel_geometry` assigns a radius to each vessel object.
#:
#: ``'data'``      — per object, r = r₁ · k where k is that object's endothelial
#:                   cells per cross-section, measured from its own cell count and
#:                   spatial extent. THE FIX.
#: ``'stereological'`` — every object at the measured human mean r̄ = 4.67 µm. No
#:                   k dependence; the conservative end of the Sargent caveat.
#: ``'constant'``  — a caller-supplied constant, for reproducing the legacy
#:                   7.5 µm behaviour under the new sub-grid painting.
RADIUS_MODES = ("data", "stereological", "constant")


# --------------------------------------------------------------------------- #
# 1. Per-object radius and footprint FROM THE DATA
# --------------------------------------------------------------------------- #
def vessel_geometry(
    vessels: Dict,
    radius_mode: str = "data",
    constant_radius_um: Optional[float] = None,
    max_ec_per_cross_section: Optional[float] = None,
) -> Dict:
    """Per-vessel-object radius, axis and footprint area, derived from the data.

    The model, in one paragraph
    ---------------------------
    A segmented object is a piece of one vessel seen in a thin section. Its
    ``N`` endothelial cells are distributed two ways: **around** the wall
    (``k`` per cross-section) and **along** the axis (one every
    :data:`EC_AXIAL_PITCH_UM` = 50 µm, measured in human cortex, PMID 38723085).
    The object's own in-plane extent ``E`` says how much of the split is axial,
    so the number of cells per cross-section is::

        k = N / max(1, E / 50 µm)

    and the radius follows from circumferential tiling, ``r = r₁ · k`` with
    ``r₁`` = :data:`SINGLE_EC_VESSEL_RADIUS_UM` = 3.11 µm derived from human EM.
    The footprint is then everything within ``r`` of the vessel **axis**, which
    is a segment of length ``max(0, E − 2r)`` through the object's centroid along
    its principal direction. The two limits are the two things a section can
    show, and the formula reproduces both without a switch:

    * a **transverse cut** has ``E ≈ 2r`` (the cells sit across one lumen), the
      axis segment collapses to a point, and the footprint is a single disc of
      radius ``r`` — *one* vessel profile, not ``N`` of them. This is the case the
      per-endothelial-cell painting got wrong.
    * an **in-plane run** has ``E ≫ 2r``, ``k`` falls towards 1, and the footprint
      is a stadium of half-width ``r`` along the run.

    Parameters
    ----------
    vessels : dict
        Output of :func:`~reconstruct.vessel_network.segment_vessels`. Uses
        ``xy``, ``labels``, ``cells_per_vessel`` and ``extent_um``.
    radius_mode : {'data', 'stereological', 'constant'}
        See :data:`RADIUS_MODES`.
    constant_radius_um : float, optional
        Required for ``radius_mode='constant'``.
    max_ec_per_cross_section : float, optional
        Cap on ``k``. ``None`` (default) applies no cap; the honest reading of an
        object with a very large ``k`` is that the linking distance merged several
        vessels, and capping would hide that rather than fix it. Sargent 2023
        (PMID 37582180) observed at most 5–6 cells per cross-section in mouse
        cortex, so pass ``6.0`` for the literature-bounded variant.

    Returns
    -------
    dict
        ``radius_um``, ``ec_per_cross_section``, ``axis_length_um``,
        ``axis_p0`` / ``axis_p1`` (n, 2) segment endpoints, ``footprint_area_um2``
        per object, and the section totals.
    """
    mode = str(radius_mode).strip().lower()
    if mode not in RADIUS_MODES:
        raise ValueError(f"unknown radius_mode {radius_mode!r}; expected one of "
                         f"{RADIUS_MODES}")
    if mode == "constant" and constant_radius_um is None:
        raise ValueError("radius_mode='constant' requires constant_radius_um")

    xy = np.asarray(vessels["xy"], dtype=np.float64)
    labels = np.asarray(vessels["labels"], dtype=int)
    counts = np.asarray(vessels["cells_per_vessel"], dtype=np.float64)
    extent = np.asarray(vessels["extent_um"], dtype=np.float64)
    n_v = int(counts.size)

    if n_v == 0:
        z2 = np.zeros((0, 2), dtype=np.float64)
        z1 = np.zeros(0, dtype=np.float64)
        return {
            "radius_mode": mode, "n_vessels": 0,
            "radius_um": z1, "ec_per_cross_section": z1, "axis_length_um": z1,
            "axis_p0": z2, "axis_p1": z2, "footprint_area_um2": z1,
            "total_footprint_area_um2": 0.0,
            "single_ec_radius_um": float(SINGLE_EC_VESSEL_RADIUS_UM),
            "ec_axial_pitch_um": float(EC_AXIAL_PITCH_UM),
        }

    # --- endothelial cells per cross-section, from the object's own geometry --- #
    axial_units = np.maximum(1.0, extent / EC_AXIAL_PITCH_UM)
    k = counts / axial_units
    k = np.maximum(k, 1.0)          # a vessel wall has at least one cell in it
    if max_ec_per_cross_section is not None:
        k = np.minimum(k, float(max_ec_per_cross_section))

    if mode == "data":
        radius = SINGLE_EC_VESSEL_RADIUS_UM * k
    elif mode == "stereological":
        radius = np.full(n_v, MEAN_VESSEL_RADIUS_UM, dtype=np.float64)
    else:
        radius = np.full(n_v, float(constant_radius_um), dtype=np.float64)

    # --- the vessel AXIS: a segment through the centroid along the principal
    #     direction, of length max(0, extent - 2r) so a transverse cut collapses
    #     to a point and an in-plane run keeps its length.
    axis_len = np.maximum(0.0, extent - 2.0 * radius)
    p0 = np.zeros((n_v, 2), dtype=np.float64)
    p1 = np.zeros((n_v, 2), dtype=np.float64)
    order = np.argsort(labels, kind="stable")
    starts = np.searchsorted(labels[order], np.arange(n_v), side="left")
    ends = np.searchsorted(labels[order], np.arange(n_v), side="right")
    for kk in range(n_v):
        pts = xy[order[starts[kk]:ends[kk]]]
        c = pts.mean(axis=0)
        if axis_len[kk] > 0.0 and pts.shape[0] >= 2:
            d = pts - c
            _, _, vt = np.linalg.svd(d, full_matrices=False)
            u = vt[0]
        else:
            u = np.array([1.0, 0.0])
        half = 0.5 * axis_len[kk] * u
        p0[kk] = c - half
        p1[kk] = c + half

    # stadium area: a rectangle of length axis_len and width 2r, plus a full disc
    area = np.pi * radius ** 2 + 2.0 * radius * axis_len
    return {
        "radius_mode": mode,
        "n_vessels": n_v,
        "radius_um": radius,
        "ec_per_cross_section": k,
        "axis_length_um": axis_len,
        "axis_p0": p0,
        "axis_p1": p1,
        "footprint_area_um2": area,
        "total_footprint_area_um2": float(area.sum()),
        "single_ec_radius_um": float(SINGLE_EC_VESSEL_RADIUS_UM),
        "ec_axial_pitch_um": float(EC_AXIAL_PITCH_UM),
        "max_ec_per_cross_section": (None if max_ec_per_cross_section is None
                                     else float(max_ec_per_cross_section)),
    }


# --------------------------------------------------------------------------- #
# 2. Sub-grid fractional-area painting
# --------------------------------------------------------------------------- #
def _stadium_disc_centres(geom: Dict) -> "tuple[np.ndarray, np.ndarray]":
    """Each stadium as a chain of overlapping discs (centres, radii).

    A stadium of half-width ``r`` around a segment is exactly the union of discs
    of radius ``r`` centred on that segment, so sampling the segment at a pitch
    ``<= r`` and taking the union of those discs reproduces it to within the
    sagitta ``r(1 - cos(theta))``, which is under 3 % of ``r`` at that pitch.
    Reducing everything to "union of discs" is what lets the painter below be a
    single vectorised stencil operation.
    """
    p0 = np.asarray(geom["axis_p0"], dtype=np.float64)
    p1 = np.asarray(geom["axis_p1"], dtype=np.float64)
    r = np.asarray(geom["radius_um"], dtype=np.float64)
    L = np.asarray(geom["axis_length_um"], dtype=np.float64)
    if r.size == 0:
        return np.zeros((0, 2)), np.zeros(0)

    n_disc = np.maximum(1, np.ceil(L / np.maximum(r, 1e-9)).astype(int) + 1)
    total = int(n_disc.sum())
    centres = np.empty((total, 2), dtype=np.float64)
    radii = np.empty(total, dtype=np.float64)
    off = 0
    for i in range(r.size):
        m = int(n_disc[i])
        t = np.linspace(0.0, 1.0, m)[:, None] if m > 1 else np.zeros((1, 1))
        centres[off:off + m] = p0[i] + t * (p1[i] - p0[i])
        radii[off:off + m] = r[i]
        off += m
    return centres, radii


def _subgrid_coverage(
    domain: Domain,
    geom: Dict,
    subpitch_um: Optional[float] = None,
    mask_to_tissue: bool = True,
) -> "tuple[np.ndarray, np.ndarray]":
    """``(coverage, cell_area_um2)`` — see :func:`subgrid_coverage_field`.

    ``cell_area_um2`` is the physical area each grid cell actually owns on the
    sub-lattice. It equals ``dx²`` only when ``dx/subpitch`` is an integer, and
    callers that need an *area* rather than a per-cell fraction must weight by it
    (:func:`vascular_area_fraction` does).
    """
    ny, nx = domain.shape
    dx = float(domain.dx)
    out = np.zeros((ny, nx), dtype=np.float64)
    flat_area = np.full((ny, nx), dx * dx, dtype=np.float64)

    centres, radii = _stadium_disc_centres(geom)
    if centres.shape[0] == 0:
        return out, flat_area

    if subpitch_um is None:
        subpitch_um = float(np.clip(radii.min() / 4.0, 0.25, 1.0))
    h = float(subpitch_um)

    W, H = nx * dx, ny * dx
    nfx, nfy = int(np.ceil(W / h)), int(np.ceil(H / h))
    hit = np.zeros(nfy * nfx, dtype=bool)

    # Each disc is rasterised from a precomputed stencil. Two things decide which
    # stencil, and both exist to keep the painted area unbiased:
    #
    #  * RADIUS is grouped LOGARITHMICALLY at 1 %, a *relative* tolerance, so the
    #    stencil radius stays within 1 % (2 % in area) of the true radius over the
    #    whole range. An absolute quantisation to whole sub-pixels would cost up to
    #    half a sub-pixel and would be worst exactly where it matters most — on the
    #    smallest vessels, which are the majority.
    #  * SUB-PIXEL PHASE is kept (to 1/4 of a sub-pixel) for small discs instead of
    #    snapping every centre to a lattice point. Snapping is not merely
    #    imprecise, it is *biased*: every disc of a given radius then gets the
    #    identical lattice count, so the O(1) point-count error of a small disc
    #    stops averaging out over the section and becomes a systematic offset.
    #    Above PHASE_FREE_RHO the count error is already below a percent and the
    #    phase is not worth the stencils.
    n_phase, phase_free_rho = 4, 16.0
    rho = np.maximum(radii / h, 1.0)
    step = np.log1p(0.01)
    bin_of = np.rint(np.log(rho) / step).astype(np.int64)

    fx = centres[:, 0] / h - 0.5
    fy = centres[:, 1] / h - 0.5
    ci, cj = np.floor(fx).astype(np.int64), np.floor(fy).astype(np.int64)
    small = rho <= phase_free_rho
    pi = np.where(small, np.minimum((( fx - ci) * n_phase).astype(np.int64),
                                    n_phase - 1), 0)
    pj = np.where(small, np.minimum(((fy - cj) * n_phase).astype(np.int64),
                                    n_phase - 1), 0)
    # centres without a phase are snapped instead, which is what rint() does
    ci = np.where(small, ci, np.rint(fx).astype(np.int64))
    cj = np.where(small, cj, np.rint(fy).astype(np.int64))

    key = (bin_of * n_phase + pi) * n_phase + pj
    order = np.argsort(key, kind="stable")
    bounds = np.flatnonzero(np.r_[True, np.diff(key[order]) != 0, True])
    for s0, s1 in zip(bounds[:-1], bounds[1:]):
        idx = order[s0:s1]
        b, p_i, p_j = bin_of[idx[0]], pi[idx[0]], pj[idx[0]]
        rb = float(np.exp(b * step))
        ox = (p_i + 0.5) / n_phase if rb <= phase_free_rho else 0.0
        oy = (p_j + 0.5) / n_phase if rb <= phase_free_rho else 0.0
        rr = int(np.ceil(rb + 1.0))
        d = np.arange(-rr, rr + 1)
        DI, DJ = np.meshgrid(d, d)
        inside = ((DI - ox) ** 2 + (DJ - oy) ** 2) <= rb ** 2
        di, dj = DI[inside], DJ[inside]
        if di.size == 0:
            continue
        # chunk so the (n_centres x n_stencil) broadcast stays bounded
        chunk = max(1, int(4_000_000 // di.size))
        for s in range(0, idx.size, chunk):
            c = idx[s:s + chunk]
            ii = ci[c][:, None] + di[None, :]
            jj = cj[c][:, None] + dj[None, :]
            ok = (ii >= 0) & (ii < nfx) & (jj >= 0) & (jj < nfy)
            hit[(jj[ok] * nfx + ii[ok])] = True

    # sub-pixel centre -> grid cell
    col_of = np.minimum(
        (((np.arange(nfx) + 0.5) * h) / dx).astype(np.int64), nx - 1)
    row_of = np.minimum(
        (((np.arange(nfy) + 0.5) * h) / dx).astype(np.int64), ny - 1)
    # The sub-pixels each cell ACTUALLY owns. This equals (dx/h)^2 only when dx/h
    # is an integer; when it is not, cells own one more or one fewer sub-pixel per
    # axis, and treating them all as dx^2 is what leaves a residual dx trend in the
    # area fraction (and forces a clip at 1.0). Carrying the true owned area makes
    # both the per-cell fraction and the summed area exact.
    owned = (np.bincount(col_of, minlength=nx).astype(np.float64)[None, :]
             * np.bincount(row_of, minlength=ny).astype(np.float64)[:, None])
    cell_area = owned * h * h

    lin = np.flatnonzero(hit)
    if lin.size == 0:
        return out, cell_area
    fj, fi = np.divmod(lin, nfx)
    covered = np.bincount(row_of[fj] * nx + col_of[fi],
                          minlength=ny * nx).astype(np.float64)
    out = covered.reshape(ny, nx) / np.maximum(owned, 1.0)
    if mask_to_tissue:
        out *= np.asarray(domain.tissue_mask, dtype=bool)
    return out, cell_area


def subgrid_coverage_field(
    domain: Domain,
    geom: Dict,
    subpitch_um: Optional[float] = None,
    mask_to_tissue: bool = True,
) -> np.ndarray:
    """Fractional vessel-area coverage per grid cell, painted SUB-GRID.

    Every vessel footprint is rasterised onto a lattice whose pitch is a fixed
    **physical** length (``subpitch_um``, chosen from the smallest vessel radius),
    and the resulting sub-pixels are then binned into the ``domain`` grid. A cell
    returns the *fraction of its area* that the vasculature covers, so:

    * a vessel far smaller than a grid cell still deposits exactly its own area —
      the "no vessel below 6.8 µm is representable" limitation of the binary
      centre-in-disc rule is gone;
    * the total painted area is computed **without reference to** ``domain.dx``,
      so it is grid-independent by construction. Only the assignment of sub-pixels
      to cells, and the tissue mask, involve ``dx`` at all.

    Parameters
    ----------
    domain : Domain
        Grid + tissue mask.
    geom : dict
        Output of :func:`vessel_geometry`.
    subpitch_um : float, optional
        Sub-lattice pitch. Default ``clip(min_radius/4, 0.25, 1.0)`` µm — set by
        the vessels, never by ``dx``.
    mask_to_tissue : bool
        Zero the coverage outside ``domain.tissue_mask`` (default True).

    Returns
    -------
    (ny, nx) float64 in [0, 1].
    """
    return _subgrid_coverage(domain, geom, subpitch_um=subpitch_um,
                             mask_to_tissue=mask_to_tissue)[0]


# --------------------------------------------------------------------------- #
# 3. The corrected area fraction
# --------------------------------------------------------------------------- #
def vascular_area_fraction(
    domain: Domain,
    vessels: Dict,
    radius_mode: str = "data",
    constant_radius_um: Optional[float] = None,
    max_ec_per_cross_section: Optional[float] = None,
    subpitch_um: Optional[float] = None,
) -> Dict:
    """Vascular area fraction (%) of the tissue, sub-grid and data-derived.

    Returns the coverage-weighted fraction together with everything needed to
    audit it: the per-object radius distribution, the analytic footprint total
    (which involves no grid at all and is the reference the rasterised number is
    checked against), and the constants used.
    """
    geom = vessel_geometry(vessels, radius_mode=radius_mode,
                           constant_radius_um=constant_radius_um,
                           max_ec_per_cross_section=max_ec_per_cross_section)
    cov, cell_area = _subgrid_coverage(domain, geom, subpitch_um=subpitch_um)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    # AREA over AREA, not a mean of per-cell fractions: cells own slightly
    # different physical areas whenever dx is not a whole number of sub-pitches.
    tissue_area_um2 = float(cell_area[tissue].sum())
    frac = (float((cov * cell_area)[tissue].sum()) / tissue_area_um2
            if tissue_area_um2 > 0 else float("nan"))
    r = np.asarray(geom["radius_um"], dtype=np.float64)
    return {
        "area_fraction_pct": 100.0 * frac,
        "radius_mode": geom["radius_mode"],
        "n_vessels": geom["n_vessels"],
        "median_radius_um": float(np.median(r)) if r.size else float("nan"),
        "mean_radius_um": float(r.mean()) if r.size else float("nan"),
        "p90_radius_um": float(np.percentile(r, 90)) if r.size else float("nan"),
        "median_ec_per_cross_section": (
            float(np.median(geom["ec_per_cross_section"])) if r.size else float("nan")),
        # Analytic (grid-free) footprint total. Differs from the rasterised value
        # only by the overlap between neighbouring footprints, which the raster
        # de-duplicates and this sum does not.
        "analytic_area_fraction_pct": (
            100.0 * geom["total_footprint_area_um2"] / tissue_area_um2
            if tissue_area_um2 > 0 else float("nan")),
        "tissue_area_mm2": tissue_area_um2 / 1e6,
        "single_ec_radius_um": float(SINGLE_EC_VESSEL_RADIUS_UM),
        "mean_vessel_radius_um": float(MEAN_VESSEL_RADIUS_UM),
        "coverage": cov,
    }


# --------------------------------------------------------------------------- #
# 4. Fidelity to the sample: is the measured heterogeneity preserved?
# --------------------------------------------------------------------------- #
def vascular_heterogeneity(
    domain: Domain,
    field: np.ndarray,
    endothelial_xy: np.ndarray,
    bin_um: float = 100.0,
    min_tissue_frac: float = 0.5,
) -> Dict:
    """Does a vascular field reproduce the MEASURED spatial heterogeneity?

    The point of a patient-specific reconstruction is to carry this tumour's own
    vascular architecture, including how uneven it is — in glioblastoma the
    unevenness is diagnostic, not noise. A field can have a perfectly correct area
    fraction and still be useless if it has smeared that architecture flat. This
    scores the field against the raw endothelial cell positions it came from.

    Both the field and the measured endothelial cells are aggregated into bins of
    a fixed **physical** size (``bin_um``, default 100 µm — the scale of the
    published 124 µm inter-vessel spacing), so the comparison is between two
    quantities on the same footing and does not itself depend on ``dx``.

    Returns
    -------
    dict
        ``cv_measured`` — coefficient of variation of the measured endothelial
        count per bin. This is the target.
        ``cv_field`` — the same for the field's mass per bin.
        ``heterogeneity_ratio`` — ``cv_field / cv_measured``. 1.0 means the field
        carries exactly the measured unevenness; well below 1.0 means it has been
        smoothed away; above 1.0 means the painting has sharpened it.
        ``pearson_r`` / ``spearman_r`` — bin-wise agreement between field mass and
        measured endothelial count: does the field put its vasculature *where the
        cells actually are*?
        ``empty_frac_measured`` / ``empty_frac_field`` — fraction of tissue bins
        with no vasculature at all. A real section has genuinely avascular tissue;
        a smooth kernel does not.
    """
    field = np.asarray(field, dtype=np.float64)
    xy = np.asarray(endothelial_xy, dtype=np.float64).reshape(-1, 2)
    dx = float(domain.dx)
    ny, nx = domain.shape
    tissue = np.asarray(domain.tissue_mask, dtype=bool)

    nbx = max(1, int(np.ceil(nx * dx / bin_um)))
    nby = max(1, int(np.ceil(ny * dx / bin_um)))

    gx = ((np.arange(nx) + 0.5) * dx / bin_um).astype(np.int64).clip(0, nbx - 1)
    gy = ((np.arange(ny) + 0.5) * dx / bin_um).astype(np.int64).clip(0, nby - 1)
    flat = (gy[:, None] * nbx + gx[None, :]).ravel()

    n_cells = np.bincount(flat, minlength=nbx * nby).astype(np.float64)
    n_tissue = np.bincount(flat, weights=tissue.ravel().astype(np.float64),
                           minlength=nbx * nby)
    f_mass = np.bincount(flat, weights=(field * tissue).ravel(),
                         minlength=nbx * nby)

    if xy.shape[0]:
        bx = (xy[:, 0] / bin_um).astype(np.int64).clip(0, nbx - 1)
        by = (xy[:, 1] / bin_um).astype(np.int64).clip(0, nby - 1)
        ec = np.bincount(by * nbx + bx, minlength=nbx * nby).astype(np.float64)
    else:
        ec = np.zeros(nbx * nby, dtype=np.float64)

    keep = (n_cells > 0) & (n_tissue / np.maximum(n_cells, 1.0) >= min_tissue_frac)
    ec, f_mass = ec[keep], f_mass[keep]

    def _cv(a):
        m = a.mean()
        return float(a.std() / m) if m > 0 else float("nan")

    def _rank(a):
        o = np.argsort(a, kind="stable")
        r = np.empty(a.size, dtype=np.float64)
        r[o] = np.arange(a.size, dtype=np.float64)
        return r

    def _corr(a, b):
        if a.size < 3 or a.std() == 0 or b.std() == 0:
            return float("nan")
        return float(np.corrcoef(a, b)[0, 1])

    cv_m, cv_f = _cv(ec), _cv(f_mass)
    return {
        "bin_um": float(bin_um),
        "n_bins": int(ec.size),
        "cv_measured": cv_m,
        "cv_field": cv_f,
        "heterogeneity_ratio": (cv_f / cv_m if np.isfinite(cv_m) and cv_m > 0
                                else float("nan")),
        "pearson_r": _corr(ec, f_mass),
        "spearman_r": _corr(_rank(ec), _rank(f_mass)),
        "empty_frac_measured": float((ec == 0).mean()) if ec.size else float("nan"),
        "empty_frac_field": (float((f_mass == 0).mean()) if f_mass.size
                             else float("nan")),
    }


# --------------------------------------------------------------------------- #
# Self-test on a synthetic section with an ANALYTICALLY KNOWN area
# --------------------------------------------------------------------------- #
if __name__ == "__main__":  # pragma: no cover
    from ..core import CellPopulation, CellType
    from .vessel_network import segment_vessels

    print(f"r_bar  (human EM, derived)      = {MEAN_VESSEL_RADIUS_UM:.3f} um")
    print(f"r_1    (one EC per x-section)   = {SINGLE_EC_VESSEL_RADIUS_UM:.3f} um")
    print(f"EC axial pitch                  = {EC_AXIAL_PITCH_UM:.1f} um")

    # Isolated single endothelial cells on a 200 um lattice: N discs of radius r_1,
    # no overlap, so the area fraction is known in closed form.
    g = np.arange(6) * 200.0 + 100.0
    X, Y = np.meshgrid(g, g)
    x, y = X.ravel(), Y.ravel()
    cells = CellPopulation(
        x=x, y=y,
        cell_type=np.full(x.size, int(CellType.ENDOTHELIAL), dtype=np.int8),
        state=np.zeros(x.size, dtype=np.int8), sensitivity=np.ones(x.size))
    exact = x.size * np.pi * SINGLE_EC_VESSEL_RADIUS_UM ** 2 / (1200.0 ** 2)
    print(f"\nexact area fraction of {x.size} isolated vessels = {100*exact:.4f}%")
    for dx in (30.0, 20.0, 12.0, 8.0, 6.0):
        n = int(np.ceil(1200.0 / dx))
        dom = Domain(nx=n, ny=n, dx=dx)
        v = segment_vessels(cells, dom)
        af = vascular_area_fraction(dom, v)
        print(f"  dx={dx:5.1f}  subgrid {af['area_fraction_pct']:7.4f}%  "
              f"(ratio to exact {af['area_fraction_pct']/(100*exact):.4f})")
    print("\nSELF-TEST DONE")
