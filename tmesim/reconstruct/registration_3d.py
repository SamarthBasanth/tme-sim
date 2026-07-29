"""Serial-section 3D registration: align two independently-reconstructed 2D
tissue sections into one shared (x, y, z) frame (`reconstruct` stage, additive).

Why this module exists
-----------------------
``ARCHITECTURE_CONTRACT.md``'s ``Domain`` is (ny, nx) - there is no z-axis
anywhere in ``tmesim.core``, and this module does not add one there. tme-sim's
stated 3D strategy is instead to BUILD the third dimension itself, by taking
multiple 2D sections cut at a known z-spacing through the same tissue block and
registering their (x, y) tissue coordinate frames onto each other - the way
PASTE2 and STalign do - rather than reading a native 3D coordinate from a
vendor file (there is no such vendor file; see
``results/landscape_spatial_platforms.md`` and
``results/landscape_serial_sections.md``). This module is that registration
step. It is a NEW module; it does not touch ``reconstruct``'s existing 2D
path (``assemble.py``, ``spatial.py``, ``vasculature.py``, ...), ``core``,
``physics``, ``sim`` or ``molecular``.

What "registration" means here, precisely - read before trusting the output
-----------------------------------------------------------------------------
Two adjacent 5 um Xenium sections of the same FFPE block do NOT contain the
same cells: a cell present in Rep1 is, in general, physically absent from
Rep2 (it was sectioned through, or sits one slice away). So there is no such
thing as ground-truth cell-to-cell correspondence between the two point
clouds, and nothing in this module claims to find one. "Registration" here
means finding the rigid-plus-uniform-scale (similarity) transform of the
(x, y) TISSUE COORDINATE FRAME of the moving section that best overlaps its
TISSUE GEOMETRY (cell density / outline) with the reference section's - i.e.
the same class of alignment PASTE2 calls the "spatial" step and STalign calls
diffeomorphic image registration, done here with the simpler and more
honestly-scoped tool: a MULTI-START search (identity, centroid-only
translation, and several coarse rotation/scale grid-search seeds) each
refined by trimmed point-cloud ICP (Besl & McKay 1992) with a Umeyama (1991)
closed-form similarity-transform fit at each iteration, with the final answer
picked by a fine boundary-outline metric rather than the coarse search score.
The multi-start design is not cosmetic - a single-seed version of this
algorithm was tried first, converged to a plausible-looking but WRONG 20-
degree rotation on the real breast-cancer data this module targets (the
coarse correlation grid aliased on repeated lobular tissue structure), and
the multi-start fix is documented in detail, with the actual before/after
numbers, in ``results/registration_3d_breast_poc.md``. No optimal transport,
no gene-expression term, no non-rigid warp - see that same file for what a
fuller PASTE2-style method would add and why it was not attempted this round.

The z-axis this module produces is NOT measured. It is an ASSUMED constant
offset (``z_spacing_um``) between the reference section (z=0) and the moving
section (z=z_spacing_um), because Janesick et al. 2023 documents that the two
Xenium replicates are serial 5 um sections cut from the same block but does
NOT publish the exact number of microns separating Rep1 from Rep2 specifically
(``results/landscape_serial_sections.md``, "Sectioning / z-spacing" - flagged
UNVERIFIED there). Every caller of :func:`register_serial_sections` must pass
an explicit, labelled ``z_spacing_um`` for this reason; there is no silent
default presented as measured.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

__all__ = [
    "SimilarityTransform",
    "RegistrationQuality",
    "SerialSectionRegistration",
    "rasterize_density",
    "tissue_boundary_points",
    "density_correlation",
    "chamfer_distance",
    "hausdorff_percentile",
    "coarse_rotation_search",
    "icp_similarity",
    "register_serial_sections",
    "Z_SPACING_UM_ASSUMED",
    "Z_SPACING_CAVEAT",
]

#: Placeholder z-spacing between the Janesick et al. 2023 Xenium Rep1/Rep2
#: sections, in micrometres. NOT read from the paper - the paper states each
#: section is 5 um thick and that the two replicates are serial, but not how
#: many microns (how many intervening sections/curls) separate them
#: specifically. 10 um is the midpoint of the plausible single-digit-to-tens-
#: of-microns range one serial 5 um section apart would imply if the two
#: replicates were cut back-to-back; treat it as a labelled placeholder, not a
#: measurement. See Z_SPACING_CAVEAT and results/registration_3d_breast_poc.md.
Z_SPACING_UM_ASSUMED: float = 10.0

Z_SPACING_CAVEAT: str = (
    "z-spacing between Xenium Rep1 and Rep2 is UNVERIFIED: Janesick et al. "
    "2023 documents 5 um section thickness and that both replicates are "
    "serial sections of the same FFPE block, but does not publish the exact "
    "micron offset between Rep1 and Rep2 specifically (results/"
    "landscape_serial_sections.md). The value used here is an assumed "
    "placeholder in the plausible 5-15 um range, not a measured quantity."
)


# --------------------------------------------------------------------------- #
# Small linear-algebra pieces
# --------------------------------------------------------------------------- #
def _rotation_matrix(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]], dtype=np.float64)


@dataclass
class SimilarityTransform:
    """A 2D rotation + uniform scale + translation: ``y = s * R @ x + t``."""

    rotation: np.ndarray  # (2, 2)
    scale: float
    translation: np.ndarray  # (2,)

    def apply(self, xy: np.ndarray) -> np.ndarray:
        xy = np.asarray(xy, dtype=np.float64)
        # errstate: numpy's Accelerate-backed matmul on Apple Silicon emits
        # spurious divide/overflow/invalid RuntimeWarnings on some (n, 2)
        # shapes even though every output value is finite (checked directly -
        # see results/registration_3d_breast_poc.md); silenced here rather
        # than left to alarm a reader about a numerical fault that isn't one.
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            return self.scale * (xy @ self.rotation.T) + self.translation

    @property
    def rotation_deg(self) -> float:
        return float(np.degrees(np.arctan2(self.rotation[1, 0], self.rotation[0, 0])))

    def to_dict(self) -> Dict:
        return {
            "rotation_deg": self.rotation_deg,
            "scale": float(self.scale),
            "translation_um": [float(self.translation[0]), float(self.translation[1])],
        }


def _umeyama(source: np.ndarray, target: np.ndarray, allow_scale: bool = True
             ) -> SimilarityTransform:
    """Closed-form similarity transform mapping ``source`` onto ``target``.

    Umeyama, 1991, "Least-squares estimation of transformation parameters
    between two point patterns" - the standard SVD solution for the rigid (or
    rigid+uniform-scale) transform minimizing sum-of-squared point-to-point
    distances GIVEN a correspondence (here, the correspondence is whatever the
    caller matched - a coarse-search grid or an ICP nearest-neighbour pairing,
    never a claim of biological identity).
    """
    source = np.asarray(source, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    mu_s = source.mean(axis=0)
    mu_t = target.mean(axis=0)
    sc = source - mu_s
    tc = target - mu_t
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        # See the matching comment in SimilarityTransform.apply - benign
        # Accelerate/numpy matmul warning on some shapes, not a real fault.
        cov = (tc.T @ sc) / source.shape[0]
    u, d, vt = np.linalg.svd(cov)
    sgn = np.eye(2)
    if np.linalg.det(u) * np.linalg.det(vt) < 0:
        sgn[1, 1] = -1.0
    r = u @ sgn @ vt
    if allow_scale:
        var_s = (sc ** 2).sum() / source.shape[0]
        s = float(np.trace(np.diag(d) @ sgn) / var_s) if var_s > 0 else 1.0
    else:
        s = 1.0
    t = mu_t - s * (r @ mu_s)
    return SimilarityTransform(rotation=r, scale=s, translation=t)


# --------------------------------------------------------------------------- #
# Rasterization + quality metrics (used both for the search and for reporting)
# --------------------------------------------------------------------------- #
def rasterize_density(
    xy: np.ndarray, bin_um: float, extent: Optional[Tuple[float, float, float, float]] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """2D histogram of point counts, in um bins. Returns ``(grid, x_edges, y_edges)``.

    ``grid`` has shape ``(n_y_bins, n_x_bins)`` (row-major, matching the
    ``(ny, nx)`` convention the rest of the codebase uses for fields), even
    though this module deliberately does not build a ``core.Domain`` - two
    whole sections at cell-cloud resolution is not the simulation grid.
    """
    xy = np.asarray(xy, dtype=np.float64)
    if extent is None:
        pad = 2.0 * bin_um
        x_lo, x_hi = xy[:, 0].min() - pad, xy[:, 0].max() + pad
        y_lo, y_hi = xy[:, 1].min() - pad, xy[:, 1].max() + pad
    else:
        x_lo, x_hi, y_lo, y_hi = extent
    nx = max(int(np.ceil((x_hi - x_lo) / bin_um)), 1)
    ny = max(int(np.ceil((y_hi - y_lo) / bin_um)), 1)
    x_edges = x_lo + bin_um * np.arange(nx + 1)
    y_edges = y_lo + bin_um * np.arange(ny + 1)
    grid, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=(y_edges, x_edges))
    return grid, x_edges, y_edges


def density_correlation(xy_a: np.ndarray, xy_b: np.ndarray, bin_um: float = 60.0) -> float:
    """Pearson correlation of two point clouds' density fields on a SHARED grid.

    The grid covers the union bounding box of both point sets (with margin),
    at ``bin_um`` resolution, so ``xy_a`` and ``xy_b`` are rasterized onto
    IDENTICAL bins and are directly comparable cell-by-cell. This is the
    "does the tissue mass line up" metric: 1.0 means the two density fields
    are a perfect linear match, 0.0 no linear relationship, and it is
    computed identically before and after registration so the before/after
    delta is meaningful.
    """
    xy_a = np.asarray(xy_a, dtype=np.float64)
    xy_b = np.asarray(xy_b, dtype=np.float64)
    pad = 2.0 * bin_um
    x_lo = min(xy_a[:, 0].min(), xy_b[:, 0].min()) - pad
    x_hi = max(xy_a[:, 0].max(), xy_b[:, 0].max()) + pad
    y_lo = min(xy_a[:, 1].min(), xy_b[:, 1].min()) - pad
    y_hi = max(xy_a[:, 1].max(), xy_b[:, 1].max()) + pad
    extent = (x_lo, x_hi, y_lo, y_hi)
    grid_a, _, _ = rasterize_density(xy_a, bin_um, extent)
    grid_b, _, _ = rasterize_density(xy_b, bin_um, extent)
    fa, fb = grid_a.ravel(), grid_b.ravel()
    if fa.std() <= 0 or fb.std() <= 0:
        return 0.0
    r = float(np.corrcoef(fa, fb)[0, 1])
    return 0.0 if np.isnan(r) else r


def tissue_boundary_points(
    xy: np.ndarray, bin_um: float = 30.0, close_iterations: int = 1
) -> np.ndarray:
    """Approximate tissue OUTLINE as the boundary pixels of an occupancy mask.

    Rasterizes ``xy`` at ``bin_um`` resolution, morphologically closes small
    gaps between cells (``close_iterations`` rounds of binary closing) so the
    mask is one solid tissue silhouette rather than a Swiss-cheese of
    individual-cell pixels, then returns the (x, y) um centres of pixels on
    the mask's boundary (mask AND NOT eroded-mask). This is a cheap, honest
    stand-in for a proper alpha-shape / concave hull - fine for a Hausdorff/
    Chamfer distance where what matters is "does the silhouette overlap",
    not a millimetre-precise contour.
    """
    grid, x_edges, y_edges = rasterize_density(xy, bin_um)
    mask = grid > 0
    if close_iterations > 0:
        mask = ndimage.binary_closing(mask, iterations=close_iterations)
    eroded = ndimage.binary_erosion(mask)
    boundary = mask & ~eroded
    rows, cols = np.nonzero(boundary)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
    return np.column_stack([x_centers[cols], y_centers[rows]])


def chamfer_distance(pts_a: np.ndarray, pts_b: np.ndarray) -> float:
    """Symmetric mean nearest-neighbour distance between two point sets (um).

    The standard Chamfer distance: mean over ``a`` of its distance to the
    nearest ``b`` point, plus the mirror, averaged. Reported in preference to
    the raw (max) Hausdorff distance as the HEADLINE boundary metric because
    a raw max is dominated by a single worst-case outlier pixel on a noisy,
    real-tissue boundary; :func:`hausdorff_percentile` is reported alongside
    it for the outlier-sensitive view.
    """
    pts_a = np.asarray(pts_a, dtype=np.float64)
    pts_b = np.asarray(pts_b, dtype=np.float64)
    tree_a = cKDTree(pts_a)
    tree_b = cKDTree(pts_b)
    d_ab, _ = tree_b.query(pts_a, k=1)
    d_ba, _ = tree_a.query(pts_b, k=1)
    return float(0.5 * (d_ab.mean() + d_ba.mean()))


def hausdorff_percentile(pts_a: np.ndarray, pts_b: np.ndarray, pct: float = 95.0) -> float:
    """The ``pct``-th percentile of the symmetric nearest-neighbour distances (um).

    A robust stand-in for the classical (max) Hausdorff distance, which a
    single stray boundary pixel can dominate. ``pct=100`` reproduces the
    classical directed-then-symmetrized Hausdorff distance exactly.
    """
    pts_a = np.asarray(pts_a, dtype=np.float64)
    pts_b = np.asarray(pts_b, dtype=np.float64)
    tree_a = cKDTree(pts_a)
    tree_b = cKDTree(pts_b)
    d_ab, _ = tree_b.query(pts_a, k=1)
    d_ba, _ = tree_a.query(pts_b, k=1)
    return float(np.percentile(np.concatenate([d_ab, d_ba]), pct))


# --------------------------------------------------------------------------- #
# Coarse global search + ICP refinement
# --------------------------------------------------------------------------- #
def coarse_rotation_search(
    xy_ref: np.ndarray,
    xy_mov: np.ndarray,
    bin_um: float = 100.0,
    n_angles: int = 72,
    scales: Tuple[float, ...] = (1.0,),
    top_k: int = 1,
) -> List[SimilarityTransform]:
    """Grid search over rotation (and optionally scale) to seed ICP.

    Cell point clouds from a tumour section are often roughly blob-shaped
    with no single obvious landmark, so a naive ICP started from identity can
    converge to a plausible-looking but wrong local optimum. This picks
    candidate starting rotations (and, if ``scales`` has more than one entry,
    scales) by centring both point clouds on their centroids and ranking
    ``(angle, scale)`` pairs by :func:`density_correlation` on a COARSE grid
    (fast: ``n_angles`` x ``len(scales)`` rasterizations at ``bin_um``
    resolution, not a per-point search).

    IMPORTANT CAVEAT, found the hard way on the real breast-cancer data this
    module was built against: breast tissue is full of similarly-sized,
    roughly circular glandular/lobular structures repeated across the
    section. A coarse density grid can ALIAS on that repetition - i.e. a
    wrong rotation can score a higher coarse correlation than the true
    alignment simply because it phase-matches one lobule to a different,
    similar-looking lobule. Returning only the single best-scoring candidate
    was observed to do exactly that (see
    ``results/registration_3d_breast_poc.md``): it seeded ICP into a
    20-degree rotation that looked BETTER on the coarse grid but produced a
    visually and metrically WORSE alignment than doing nothing. The fix is
    architectural, not a threshold tweak: return the top ``top_k``
    candidates (plus the caller always separately tries identity/
    centroid-only - see :func:`register_serial_sections`) and let the FINAL
    boundary-based quality metric, computed after ICP on each candidate,
    decide the winner - never trust the coarse correlation ranking alone.
    """
    xy_ref = np.asarray(xy_ref, dtype=np.float64)
    xy_mov = np.asarray(xy_mov, dtype=np.float64)
    ref_centroid = xy_ref.mean(axis=0)
    mov_centroid = xy_mov.mean(axis=0)
    ref_c = xy_ref - ref_centroid
    mov_c = xy_mov - mov_centroid

    pad = 2.0 * bin_um
    half = max(np.abs(ref_c).max(), np.abs(mov_c).max()) + pad
    x_edges = np.arange(-half, half + bin_um, bin_um)
    y_edges = np.arange(-half, half + bin_um, bin_um)
    extent = (x_edges[0], x_edges[-1], y_edges[0], y_edges[-1])
    ref_grid, _, _ = rasterize_density(ref_c, bin_um, extent)
    ref_flat = ref_grid.ravel()
    ref_ok = ref_flat.std() > 0

    scored: List[Tuple[float, float, float]] = []  # (corr, angle, scale)
    for theta in np.linspace(0.0, 2.0 * np.pi, n_angles, endpoint=False):
        rot = _rotation_matrix(theta)
        with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
            # See the matching comment in SimilarityTransform.apply.
            rotated = mov_c @ rot.T
        for s in scales:
            mov_grid, _, _ = rasterize_density(rotated * s, bin_um, extent)
            mov_flat = mov_grid.ravel()
            if not ref_ok or mov_flat.std() <= 0:
                corr = 0.0
            else:
                corr = float(np.corrcoef(ref_flat, mov_flat)[0, 1])
                if np.isnan(corr):
                    corr = 0.0
            scored.append((corr, theta, s))

    scored.sort(key=lambda c: c[0], reverse=True)
    out: List[SimilarityTransform] = []
    for corr, theta, s in scored[:max(top_k, 1)]:
        r0 = _rotation_matrix(theta)
        t0 = ref_centroid - s * (r0 @ mov_centroid)
        out.append(SimilarityTransform(rotation=r0, scale=s, translation=t0))
    return out


def icp_similarity(
    xy_ref: np.ndarray,
    xy_mov: np.ndarray,
    init: Optional[SimilarityTransform] = None,
    max_iterations: int = 40,
    tol_um: float = 0.05,
    sample_size: int = 6000,
    trim_fraction: float = 0.8,
    allow_scale: bool = True,
    seed: int = 0,
) -> Tuple[SimilarityTransform, List[Dict]]:
    """Trimmed point-to-point ICP (Besl & McKay 1992) with a Umeyama refit.

    NOT cell matching. Each iteration finds, for a FIXED random subsample of
    ``xy_mov`` (fixed by ``seed``, so the whole procedure is deterministic),
    the current nearest neighbour in ``xy_ref`` under the running transform;
    the closest ``trim_fraction`` of those pairs (by distance) are kept and
    refit with :func:`_umeyama` to produce the next transform. Trimming keeps
    the fit from being dragged around by the fraction of each section that
    the other section's tissue does not cover (real serial sections overlap
    partially, not perfectly - PASTE2's "partial" in its name is the same
    problem; this is the cheaper, non-optimal-transport answer to it).

    Returns ``(transform, history)`` where ``history`` is one dict per
    iteration (``{"rmse_um", "n_pairs", "rotation_deg", "scale"}``) so
    convergence can be inspected/plotted.
    """
    xy_ref = np.asarray(xy_ref, dtype=np.float64)
    xy_mov = np.asarray(xy_mov, dtype=np.float64)
    rng = np.random.default_rng(seed)
    n = xy_mov.shape[0]
    if sample_size < n:
        idx = rng.choice(n, size=sample_size, replace=False)
        idx.sort()  # deterministic, reproducible ordering
    else:
        idx = np.arange(n)
    mov_sample = xy_mov[idx]

    ref_tree_size = min(xy_ref.shape[0], 40000)
    if ref_tree_size < xy_ref.shape[0]:
        ref_idx = rng.choice(xy_ref.shape[0], size=ref_tree_size, replace=False)
        ref_idx.sort()
        ref_for_tree = xy_ref[ref_idx]
    else:
        ref_for_tree = xy_ref
    tree = cKDTree(ref_for_tree)

    transform = init if init is not None else SimilarityTransform(
        rotation=np.eye(2), scale=1.0, translation=np.zeros(2))

    history: List[Dict] = []
    prev_rmse = np.inf
    for it in range(max_iterations):
        transformed = transform.apply(mov_sample)
        dist, nn_idx = tree.query(transformed, k=1)
        keep_n = max(int(np.ceil(trim_fraction * dist.size)), 3)
        keep = np.argsort(dist)[:keep_n]
        matched_ref = ref_for_tree[nn_idx[keep]]
        source_pts = mov_sample[keep]
        transform = _umeyama(source_pts, matched_ref, allow_scale=allow_scale)
        rmse = float(np.sqrt(np.mean(dist[keep] ** 2)))
        history.append({
            "iteration": it, "rmse_um": rmse, "n_pairs": int(keep_n),
            "rotation_deg": transform.rotation_deg, "scale": float(transform.scale),
        })
        if abs(prev_rmse - rmse) < tol_um:
            break
        prev_rmse = rmse

    return transform, history


# --------------------------------------------------------------------------- #
# Top-level result + entry point
# --------------------------------------------------------------------------- #
@dataclass
class RegistrationQuality:
    """Before/after quality metrics. Every field is measured, none assumed."""

    density_correlation_before: float
    density_correlation_after: float
    chamfer_um_before: float
    chamfer_um_after: float
    hausdorff95_um_before: float
    hausdorff95_um_after: float

    def to_dict(self) -> Dict:
        return {
            "density_correlation_before": self.density_correlation_before,
            "density_correlation_after": self.density_correlation_after,
            "chamfer_um_before": self.chamfer_um_before,
            "chamfer_um_after": self.chamfer_um_after,
            "hausdorff95_um_before": self.hausdorff95_um_before,
            "hausdorff95_um_after": self.hausdorff95_um_after,
        }


@dataclass
class SerialSectionRegistration:
    """Two serial sections aligned into one (x, y, z) point cloud.

    ``points_xyz`` is ``(n_ref + n_mov, 3)``: the reference section's own
    (untouched) (x, y) at z=0, stacked with the moving section's REGISTERED
    (x, y) at z=``z_spacing_um``. ``section_id`` is 0 for reference rows, 1
    for moving rows, same row order. This is a 2.5D reconstruction, not a
    volumetric one - see the module docstring and
    ``results/registration_3d_breast_poc.md`` for exactly what that does and
    does not establish.
    """

    ref_sample_id: str
    mov_sample_id: str
    transform: SimilarityTransform
    z_spacing_um: float
    z_spacing_caveat: str
    points_xyz: np.ndarray
    section_id: np.ndarray
    cell_ids: np.ndarray
    quality: RegistrationQuality
    icp_history: List[Dict]
    candidates: List[Dict] = field(default_factory=list)
    method: str = field(default=(
        "multi-start: identity, centroid-only-translation, and the top coarse "
        "rotation/scale grid-search seeds (density-field correlation) are each "
        "refined by trimmed point-to-point ICP (Besl & McKay 1992, closed-form "
        "Umeyama 1991 similarity-transform refit per iteration); the FINAL "
        "candidate is selected by post-ICP boundary Chamfer distance (not by "
        "the coarse correlation score, which was found to alias on repeated "
        "tissue structure - see coarse_rotation_search's docstring and "
        "results/registration_3d_breast_poc.md). Aligns TISSUE GEOMETRY (cell "
        "density / outline), not individual cells - no cell-to-cell "
        "correspondence is claimed or used."
    ))

    def summary(self) -> Dict:
        return {
            "ref_sample_id": self.ref_sample_id,
            "mov_sample_id": self.mov_sample_id,
            "transform": self.transform.to_dict(),
            "z_spacing_um": self.z_spacing_um,
            "z_spacing_caveat": self.z_spacing_caveat,
            "n_points": int(self.points_xyz.shape[0]),
            "n_ref": int((self.section_id == 0).sum()),
            "n_mov": int((self.section_id == 1).sum()),
            "quality": self.quality.to_dict(),
            "candidates": self.candidates,
            "method": self.method,
        }


def register_serial_sections(
    xy_ref: np.ndarray,
    xy_mov: np.ndarray,
    z_spacing_um: float,
    z_spacing_caveat: str = Z_SPACING_CAVEAT,
    ref_sample_id: str = "ref",
    mov_sample_id: str = "mov",
    ref_cell_ids: Optional[np.ndarray] = None,
    mov_cell_ids: Optional[np.ndarray] = None,
    coarse_bin_um: float = 100.0,
    coarse_n_angles: int = 72,
    coarse_scales: Tuple[float, ...] = (0.95, 1.0, 1.05),
    coarse_top_k: int = 5,
    icp_max_iterations: int = 40,
    icp_sample_size: int = 6000,
    icp_trim_fraction: float = 0.8,
    allow_scale: bool = True,
    quality_bin_um: float = 60.0,
    boundary_bin_um: float = 30.0,
    seed: int = 0,
) -> SerialSectionRegistration:
    """Register the moving section's (x, y) frame onto the reference section's.

    Deterministic given ``seed`` (the ICP subsampling is the only randomness
    in the pipeline; the coarse search is an exhaustive grid, not sampled).

    Parameters
    ----------
    xy_ref, xy_mov : (n, 2) arrays, micrometres
        Per-cell centroids of the two sections, e.g.
        ``ingest.SampleRecord.cell_xy_um`` for each replicate. ``xy_ref`` is
        left untouched (it defines the shared frame); ``xy_mov`` is
        transformed to align onto it.
    z_spacing_um : float
        REQUIRED, no default: the assumed z-offset between the two sections
        in micrometres. Pass :data:`Z_SPACING_UM_ASSUMED` explicitly if using
        the module's placeholder value, so the caller cannot receive an
        unlabelled number - see :data:`Z_SPACING_CAVEAT`.
    z_spacing_caveat : str
        The prose caveat attached to the returned
        :class:`SerialSectionRegistration`, defaulting to
        :data:`Z_SPACING_CAVEAT` (the Janesick et al. 2023 breast Xenium
        Rep1/Rep2 caveat this module was first built against). A DIFFERENT
        dataset's z-spacing provenance is NOT the same claim - e.g. the
        Greenwald et al. 2024 GBM Visium/CODEX pair
        (``results/gbm_3d_registration.md``) is "near-adjacent sections,
        exact micron offset not published" for a different reason (a
        different paper, different platforms) than the breast pair is. This
        parameter exists so a caller registering a NEW dataset must pass its
        OWN caveat rather than silently inheriting one written about a
        different dataset - this is a documentation-string parameter, not a
        change to the registration algorithm itself.
    allow_scale : bool
        If True (default) the fit is a similarity transform (rotation + one
        uniform scale + translation); if False it is rigid (scale pinned to
        1.0). Two independent Xenium runs can differ by a small calibration
        scale factor, so allowing it is the more defensible default for
        cross-run registration - but it is reported either way in
        ``transform.scale`` so a rigid-only reader can check how far from 1.0
        it moved.

    Returns
    -------
    SerialSectionRegistration
    """
    xy_ref = np.asarray(xy_ref, dtype=np.float64)
    xy_mov = np.asarray(xy_mov, dtype=np.float64)
    if xy_ref.ndim != 2 or xy_ref.shape[1] != 2:
        raise ValueError(f"xy_ref must be (n, 2); got {xy_ref.shape}")
    if xy_mov.ndim != 2 or xy_mov.shape[1] != 2:
        raise ValueError(f"xy_mov must be (n, 2); got {xy_mov.shape}")

    # --- quality BEFORE any transform ------------------------------------- #
    corr_before = density_correlation(xy_ref, xy_mov, bin_um=quality_bin_um)
    bpts_ref = tissue_boundary_points(xy_ref, bin_um=boundary_bin_um)
    bpts_mov_before = tissue_boundary_points(xy_mov, bin_um=boundary_bin_um)
    chamfer_before = chamfer_distance(bpts_ref, bpts_mov_before)
    hd95_before = hausdorff_percentile(bpts_ref, bpts_mov_before, pct=95.0)

    # --- multi-start: identity, centroid-only, and top coarse-search seeds -- #
    # A single "best coarse correlation" seed was found (on this dataset) to
    # alias onto a WORSE alignment than doing nothing, because breast tissue
    # has repeated similarly-shaped lobular structure - see
    # coarse_rotation_search's docstring. So every plausible seed is refined
    # by ICP and the FINAL candidate is chosen by the boundary Chamfer metric,
    # never by the coarse correlation score alone.
    seeds: List[SimilarityTransform] = [
        SimilarityTransform(rotation=np.eye(2), scale=1.0, translation=np.zeros(2)),
        SimilarityTransform(rotation=np.eye(2), scale=1.0,
                            translation=xy_ref.mean(axis=0) - xy_mov.mean(axis=0)),
    ]
    seeds.extend(coarse_rotation_search(
        xy_ref, xy_mov, bin_um=coarse_bin_um, n_angles=coarse_n_angles,
        scales=coarse_scales if allow_scale else (1.0,), top_k=coarse_top_k,
    ))

    candidates: List[Dict] = []
    best_transform: Optional[SimilarityTransform] = None
    best_history: List[Dict] = []
    best_chamfer = np.inf
    for seed_idx, seed_transform in enumerate(seeds):
        cand_transform, cand_history = icp_similarity(
            xy_ref, xy_mov, init=seed_transform, max_iterations=icp_max_iterations,
            sample_size=icp_sample_size, trim_fraction=icp_trim_fraction,
            allow_scale=allow_scale, seed=seed,
        )
        cand_xy = cand_transform.apply(xy_mov)
        cand_bpts = tissue_boundary_points(cand_xy, bin_um=boundary_bin_um)
        cand_chamfer = chamfer_distance(bpts_ref, cand_bpts)
        cand_corr = density_correlation(xy_ref, cand_xy, bin_um=quality_bin_um)
        candidates.append({
            "seed_index": seed_idx,
            "seed_rotation_deg": seed_transform.rotation_deg,
            "final_rotation_deg": cand_transform.rotation_deg,
            "final_scale": float(cand_transform.scale),
            "final_chamfer_um": cand_chamfer,
            "final_density_correlation": cand_corr,
        })
        if cand_chamfer < best_chamfer:
            best_chamfer = cand_chamfer
            best_transform = cand_transform
            best_history = cand_history

    transform = best_transform
    history = best_history
    xy_mov_registered = transform.apply(xy_mov)

    # --- quality AFTER, same metrics, same definitions ---------------------- #
    corr_after = density_correlation(xy_ref, xy_mov_registered, bin_um=quality_bin_um)
    bpts_mov_after = tissue_boundary_points(xy_mov_registered, bin_um=boundary_bin_um)
    chamfer_after = chamfer_distance(bpts_ref, bpts_mov_after)
    hd95_after = hausdorff_percentile(bpts_ref, bpts_mov_after, pct=95.0)

    quality = RegistrationQuality(
        density_correlation_before=corr_before,
        density_correlation_after=corr_after,
        chamfer_um_before=chamfer_before,
        chamfer_um_after=chamfer_after,
        hausdorff95_um_before=hd95_before,
        hausdorff95_um_after=hd95_after,
    )

    n_ref, n_mov = xy_ref.shape[0], xy_mov.shape[0]
    z = np.concatenate([np.zeros(n_ref), np.full(n_mov, float(z_spacing_um))])
    points_xyz = np.column_stack([
        np.concatenate([xy_ref[:, 0], xy_mov_registered[:, 0]]),
        np.concatenate([xy_ref[:, 1], xy_mov_registered[:, 1]]),
        z,
    ])
    section_id = np.concatenate([np.zeros(n_ref, dtype=int), np.ones(n_mov, dtype=int)])

    if ref_cell_ids is None:
        ref_cell_ids = np.array([f"{ref_sample_id}:{i}" for i in range(n_ref)], dtype=str)
    if mov_cell_ids is None:
        mov_cell_ids = np.array([f"{mov_sample_id}:{i}" for i in range(n_mov)], dtype=str)
    cell_ids = np.concatenate([np.asarray(ref_cell_ids, dtype=str),
                                np.asarray(mov_cell_ids, dtype=str)])

    return SerialSectionRegistration(
        ref_sample_id=ref_sample_id,
        mov_sample_id=mov_sample_id,
        transform=transform,
        z_spacing_um=float(z_spacing_um),
        z_spacing_caveat=z_spacing_caveat,
        points_xyz=points_xyz,
        section_id=section_id,
        cell_ids=cell_ids,
        quality=quality,
        icp_history=history,
        candidates=candidates,
    )
