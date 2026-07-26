"""Nucleus segmentation and phenotyping from tissue images.

This is the entry point of the ``reconstruct`` pipeline: it turns a (possibly
multi-channel) tissue image into a :class:`~tmesim.core.CellPopulation` — a set
of cell-center positions (in micrometers), phenotype labels, dynamic states and
intrinsic drug sensitivities — that downstream reconstruction, physics and
simulation modules consume.

Pipeline (classic marker-controlled watershed):

    1. build a single-channel "nuclear signal" (max projection over channels,
       or the grayscale image itself),
    2. Gaussian smooth to suppress speckle,
    3. Otsu threshold -> binary foreground mask,
    4. Euclidean distance transform of the mask,
    5. peak_local_max on the distance map -> one seed per nucleus,
    6. watershed to split touching nuclei,
    7. regionprops for per-nucleus centroids.

Each nucleus is phenotyped by sampling the marker channels at its centroid and
taking the strongest marker (``markers`` maps channel index -> phenotype name).
Plain grayscale images have no markers, so a deterministic intensity-banded mix
is assigned instead. Everything is deterministic and vectorized where possible.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from scipy import ndimage as ndi
from skimage.feature import peak_local_max
from skimage.filters import gaussian, threshold_otsu
from skimage.measure import regionprops
from skimage.segmentation import watershed

from tmesim.core import CellPopulation, CellState, CellType

# Default channel -> phenotype mapping used when the caller does not supply one.
# Values are channel indices into a (H, W, C) image.
DEFAULT_MARKERS: Dict[str, int] = {
    "tumor": 0,
    "stroma": 1,
    "immune": 2,
    "endothelial": 3,
}

# Marker-name -> CellType. Names are matched case-insensitively and by prefix so
# common aliases ("tumour", "endo", "macrophage") still resolve.
_MARKER_TO_CELLTYPE: Dict[str, CellType] = {
    "tumor": CellType.TUMOR,
    "tumour": CellType.TUMOR,
    "malignant": CellType.TUMOR,
    "stroma": CellType.STROMAL,
    "stromal": CellType.STROMAL,
    "immune": CellType.IMMUNE,
    "macrophage": CellType.IMMUNE,
    "tcell": CellType.IMMUNE,
    "endothelial": CellType.ENDOTHELIAL,
    "endo": CellType.ENDOTHELIAL,
    "vessel": CellType.ENDOTHELIAL,
    "necrotic": CellType.NECROTIC,
}


def _resolve_marker_type(name: str) -> CellType:
    """Map a marker name to a :class:`CellType`, tolerant of aliases/case."""
    key = name.strip().lower()
    if key in _MARKER_TO_CELLTYPE:
        return _MARKER_TO_CELLTYPE[key]
    for prefix, ctype in _MARKER_TO_CELLTYPE.items():
        if key.startswith(prefix) or prefix.startswith(key):
            return ctype
    # Unknown marker -> treat as tumor (the phenotype the pipeline cares about).
    return CellType.TUMOR


def _empty_population() -> CellPopulation:
    """A valid, zero-cell CellPopulation (for empty/degenerate images)."""
    z = np.zeros(0, dtype=np.float64)
    return CellPopulation(
        x=z.copy(),
        y=z.copy(),
        cell_type=np.zeros(0, dtype=np.int8),
        state=np.zeros(0, dtype=np.int8),
        sensitivity=z.copy(),
    )


def _nuclear_signal(image: np.ndarray) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Reduce an image to a 2D nuclear signal and return (signal, channels).

    For a (H, W, C) stack the per-pixel max projection across channels is used
    as the segmentation signal (a nucleus bright in *any* marker is detected);
    ``channels`` is the (H, W, C) float stack for later phenotyping. For a 2D
    grayscale image the image itself is the signal and ``channels`` is None.
    """
    arr = np.asarray(image, dtype=np.float64)
    if arr.ndim == 2:
        return arr, None
    if arr.ndim == 3:
        # Treat a trailing singleton channel as grayscale.
        if arr.shape[2] == 1:
            return arr[..., 0], None
        signal = arr.max(axis=2)
        return signal, arr
    raise ValueError(
        f"image must be 2D grayscale or 3D (H,W,C); got shape {arr.shape}"
    )


def _normalize(a: np.ndarray) -> np.ndarray:
    """Scale an array to [0, 1]; a flat array maps to all-zeros."""
    a = np.asarray(a, dtype=np.float64)
    lo = float(a.min())
    hi = float(a.max())
    if hi - lo <= 0.0:
        return np.zeros_like(a)
    return (a - lo) / (hi - lo)


def _segment_nuclei(
    signal: np.ndarray, dx: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Segment nuclei from a 2D signal.

    Returns ``(centroids_rc, labels)`` where ``centroids_rc`` is an (N, 2) array
    of (row, col) centroids in pixels and ``labels`` is the (H, W) label image.
    Returns an empty (0, 2) array when nothing is found.
    """
    empty = np.zeros((0, 2), dtype=np.float64)
    sig = _normalize(signal)
    if not np.any(sig > 0):
        return empty, np.zeros(signal.shape, dtype=np.int32)

    # 1-2. Gaussian smooth (sigma ~ half a typical 10 um nucleus, in pixels).
    sigma = max(1.0, 5.0 / max(dx, 1e-6))
    smooth = gaussian(sig, sigma=sigma, preserve_range=True)

    # 3. Otsu threshold. Guard against a degenerate (single-value) histogram.
    try:
        thr = threshold_otsu(smooth)
    except ValueError:
        return empty, np.zeros(signal.shape, dtype=np.int32)
    mask = smooth > thr
    # Drop tiny specks below a plausible nucleus area (~ (0.5 nucleus)^2 px).
    min_area = max(1, int(round((3.0 / max(dx, 1e-6)) ** 2)))
    mask = ndi.binary_opening(mask, iterations=1)
    if not mask.any():
        return empty, np.zeros(signal.shape, dtype=np.int32)

    # 4. Distance transform interior -> ridges at nucleus centers.
    distance = ndi.distance_transform_edt(mask)

    # 5. One seed per nucleus. min_distance sets the expected nucleus spacing.
    min_distance = max(1, int(round(6.0 / max(dx, 1e-6))))
    peaks = peak_local_max(
        distance,
        min_distance=min_distance,
        labels=mask,
        exclude_border=False,
    )
    markers_img = np.zeros(distance.shape, dtype=np.int32)
    if peaks.shape[0] == 0:
        # No interior peaks (e.g. thin structures): fall back to CC labeling.
        labels, n = ndi.label(mask)
    else:
        markers_img[tuple(peaks.T)] = np.arange(1, peaks.shape[0] + 1)
        # 6. Watershed on the inverted distance splits touching nuclei.
        labels = watershed(-distance, markers_img, mask=mask)
        n = int(labels.max())

    if n == 0:
        return empty, labels.astype(np.int32)

    # 7. Centroids via regionprops, filtering sub-nucleus fragments.
    props = regionprops(labels)
    centroids = [p.centroid for p in props if p.area >= min_area]
    if not centroids:
        # Everything got filtered; keep the largest region so we return >=1 cell.
        largest = max(props, key=lambda p: p.area)
        centroids = [largest.centroid]
    return np.asarray(centroids, dtype=np.float64), labels.astype(np.int32)


def _phenotype(
    centroids_rc: np.ndarray,
    channels: Optional[np.ndarray],
    signal: np.ndarray,
    markers: Optional[Dict[str, int]],
) -> Tuple[np.ndarray, np.ndarray]:
    """Assign a CellType and a drug-sensitivity to each nucleus.

    Returns ``(cell_type, sensitivity)`` arrays. When marker channels are
    available each nucleus takes the phenotype of its strongest marker channel;
    otherwise a deterministic intensity-banded mix is assigned. Sensitivity is a
    tumor-marker-correlated heuristic in [0, 1] (more tumor-like -> more
    sensitive, non-tumor cells resistant with mild jitter).
    """
    rows = np.clip(np.round(centroids_rc[:, 0]).astype(int), 0, signal.shape[0] - 1)
    cols = np.clip(np.round(centroids_rc[:, 1]).astype(int), 0, signal.shape[1] - 1)
    n = rows.shape[0]

    if channels is not None and markers:
        # Build parallel arrays of (channel_index, CellType) for valid channels.
        n_chan = channels.shape[2]
        names, chans, ctypes = [], [], []
        for name, ch in markers.items():
            if 0 <= int(ch) < n_chan:
                names.append(name)
                chans.append(int(ch))
                ctypes.append(int(_resolve_marker_type(name)))
        if chans:
            chan_idx = np.asarray(chans, dtype=int)
            ctype_arr = np.asarray(ctypes, dtype=np.int8)
            # (n, k) marker intensities at each centroid, per-channel normalized
            # so channels with different dynamic ranges compete fairly.
            per_chan = channels[rows[:, None], cols[:, None], chan_idx[None, :]]
            norm = np.empty_like(per_chan)
            for j in range(chan_idx.shape[0]):
                norm[:, j] = _normalize(channels[..., chan_idx[j]])[rows, cols]
            winner = np.argmax(norm, axis=1)
            cell_type = ctype_arr[winner]

            # Tumor-marker-correlated sensitivity.
            tumor_cols = [
                i for i, name in enumerate(names)
                if _resolve_marker_type(name) == CellType.TUMOR
            ]
            if tumor_cols:
                tumor_sig = norm[:, tumor_cols].max(axis=1)
            else:
                tumor_sig = _normalize(signal)[rows, cols]
            sensitivity = _sensitivity_from_tumor(tumor_sig, cell_type)
            return cell_type.astype(np.int8), sensitivity

    # Grayscale / no usable markers: deterministic intensity-banded mix.
    intensity = _normalize(signal)[rows, cols]
    # Rank-based banding keeps a mix regardless of the intensity distribution.
    order = np.argsort(np.argsort(intensity))  # 0..n-1 rank of each cell
    frac = order / max(n - 1, 1)
    cell_type = np.empty(n, dtype=np.int8)
    # Brightest band -> tumor, then immune, endothelial, dimmest -> stroma.
    cell_type[frac >= 0.55] = int(CellType.TUMOR)
    cell_type[(frac >= 0.35) & (frac < 0.55)] = int(CellType.IMMUNE)
    cell_type[(frac >= 0.20) & (frac < 0.35)] = int(CellType.ENDOTHELIAL)
    cell_type[frac < 0.20] = int(CellType.STROMAL)
    sensitivity = _sensitivity_from_tumor(intensity, cell_type)
    return cell_type, sensitivity


def _sensitivity_from_tumor(
    tumor_sig: np.ndarray, cell_type: np.ndarray
) -> np.ndarray:
    """Map a tumor-marker signal in [0,1] to drug sensitivity in [0,1].

    Tumor cells inherit sensitivity ~0.3-0.9 rising with tumor-marker intensity
    (a proxy for proliferative, drug-targetable cells); non-tumor cells are set
    resistant (~0.1-0.3) since cytotoxic response is largely a tumor readout.
    Deterministic — no RNG.
    """
    tumor_sig = np.clip(np.asarray(tumor_sig, dtype=np.float64), 0.0, 1.0)
    sens = 0.3 + 0.6 * tumor_sig
    is_tumor = cell_type == int(CellType.TUMOR)
    # Non-tumor cells: low, deterministic spread derived from their own signal.
    sens = np.where(is_tumor, sens, 0.1 + 0.2 * tumor_sig)
    return np.clip(sens, 0.0, 1.0)


def segment_and_phenotype(
    image: np.ndarray, dx: float = 15.0, markers: Optional[Dict[str, int]] = None
) -> CellPopulation:
    """Segment nuclei from a tissue image and assign phenotypes.

    Parameters
    ----------
    image:
        Either a 2D grayscale array ``(H, W)`` or a multi-channel stack
        ``(H, W, C)``. Multi-channel images are segmented on the per-pixel max
        projection so a nucleus visible in any channel is detected.
    dx:
        Physical pixel size in micrometers. Centroid pixel coordinates are
        converted to micrometers with this spacing, and it also sets the
        expected nucleus size/spacing used during segmentation.
    markers:
        Mapping of phenotype name -> channel index, e.g.
        ``{'tumor': 1, 'immune': 2, 'endothelial': 3, 'stroma': 0}``. Ignored
        for grayscale images. If None and the image is multi-channel,
        :data:`DEFAULT_MARKERS` is used (clipped to the available channels).

    Returns
    -------
    CellPopulation
        Cells with positions in micrometers, ``CellType`` phenotypes,
        ``state = VIABLE`` and tumor-correlated ``sensitivity`` in [0, 1].
        An empty (but valid) population is returned for empty/degenerate images.
    """
    if image is None:
        return _empty_population()
    arr = np.asarray(image)
    if arr.size == 0 or arr.ndim not in (2, 3):
        return _empty_population()

    signal, channels = _nuclear_signal(arr)
    if channels is not None and markers is None:
        # Default markers, clipped to channels that actually exist.
        markers = {
            k: v for k, v in DEFAULT_MARKERS.items() if v < channels.shape[2]
        }

    centroids_rc, _labels = _segment_nuclei(signal, dx)
    if centroids_rc.shape[0] == 0:
        return _empty_population()

    cell_type, sensitivity = _phenotype(centroids_rc, channels, signal, markers)

    # Pixel (row, col) -> micrometer (x, y). Column is x, row is y.
    x_um = (centroids_rc[:, 1] + 0.5) * dx
    y_um = (centroids_rc[:, 0] + 0.5) * dx
    state = np.full(cell_type.shape[0], int(CellState.VIABLE), dtype=np.int8)

    return CellPopulation(
        x=x_um,
        y=y_um,
        cell_type=cell_type.astype(np.int8),
        state=state,
        sensitivity=sensitivity.astype(np.float64),
    )


def _synth_test_image(
    seed: int = 0,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """Synthesize a small 4-channel test image with blobs per channel.

    Returns ``(image, markers)`` where image is ``(H, W, 4)`` with tumor,
    stroma, immune and endothelial channels each carrying a few Gaussian blobs.
    """
    rng = np.random.default_rng(seed)
    h, w, c = 160, 160, 4
    img = np.zeros((h, w, c), dtype=np.float64)
    yy, xx = np.mgrid[0:h, 0:w]
    markers = {"tumor": 0, "stroma": 1, "immune": 2, "endothelial": 3}
    # Blobs per channel: tumor gets the most (dominant phenotype).
    counts = {"tumor": 10, "stroma": 5, "immune": 4, "endothelial": 3}
    radius = 5.0
    for name, ch in markers.items():
        for _ in range(counts[name]):
            cy = rng.uniform(15, h - 15)
            cx = rng.uniform(15, w - 15)
            amp = rng.uniform(0.7, 1.0)
            img[..., ch] += amp * np.exp(
                -((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * radius ** 2)
            )
    img += rng.normal(0, 0.01, size=img.shape)
    return np.clip(img, 0, None), markers


if __name__ == "__main__":
    from collections import Counter

    # --- Multi-channel image -------------------------------------------------
    image, markers = _synth_test_image(seed=0)
    total_blobs = 10 + 5 + 4 + 3
    pop = segment_and_phenotype(image, dx=15.0, markers=markers)
    dist = Counter(int(t) for t in pop.cell_type)
    dist_named = {CellType(k).name: v for k, v in sorted(dist.items())}
    print(f"[multi-channel] n cells found : {pop.n} (planted {total_blobs} blobs)")
    print(f"[multi-channel] type distribution: {dist_named}")
    print(
        f"[multi-channel] sensitivity range: "
        f"[{pop.sensitivity.min():.3f}, {pop.sensitivity.max():.3f}]"
    )

    # Recovered a sensible cell count (watershed may merge overlapping blobs).
    assert pop.n >= total_blobs // 2, f"too few cells: {pop.n}"
    assert pop.n <= total_blobs * 3, f"implausibly many cells: {pop.n}"
    # A mix of types (not all one phenotype).
    assert len(dist) >= 2, f"expected a mix of cell types, got {dist_named}"
    # Tumor should be present and the states/sensitivities are valid.
    assert dist.get(int(CellType.TUMOR), 0) > 0, "no tumor cells assigned"
    assert np.all(pop.state == int(CellState.VIABLE))
    assert np.all((pop.sensitivity >= 0) & (pop.sensitivity <= 1))
    # Positions in micrometers within the field of view.
    assert np.all(pop.x >= 0) and np.all(pop.x <= 160 * 15.0)
    assert np.all(pop.y >= 0) and np.all(pop.y <= 160 * 15.0)

    # --- Grayscale image -----------------------------------------------------
    gray = image.max(axis=2)
    gpop = segment_and_phenotype(gray, dx=15.0)
    gdist = Counter(int(t) for t in gpop.cell_type)
    gdist_named = {CellType(k).name: v for k, v in sorted(gdist.items())}
    print(f"[grayscale]    n cells found : {gpop.n}")
    print(f"[grayscale]    type distribution: {gdist_named}")
    assert gpop.n > 0, "grayscale segmentation found no cells"
    assert len(gdist) >= 2, f"grayscale should still yield a mix, got {gdist_named}"

    # --- Degenerate inputs ---------------------------------------------------
    assert segment_and_phenotype(np.zeros((50, 50))).n == 0
    assert segment_and_phenotype(np.zeros((50, 50, 4)), markers=markers).n == 0
    assert segment_and_phenotype(np.array([])).n == 0
    assert segment_and_phenotype(np.full((40, 40), 3.0)).n == 0  # flat image

    print("\nAll segmentation self-tests passed.")
