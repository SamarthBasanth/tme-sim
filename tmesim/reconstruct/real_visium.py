"""Reconstruct a real spot-level ``Microenvironment`` from a REAL Visium sample.

This is the L2 reconstruction step for spot-resolution data: it consumes the
per-spot cell-type fractions from :mod:`tmesim.reconstruct.deconvolve` (which are
themselves derived only from real Visium counts + the real Darmanis reference)
and assembles them into the core :class:`Microenvironment` object.

Pipeline:

1. Deconvolve the sample -> ``(coords_um, fractions[S, K], type_names)``.
2. Build a :class:`Domain` at the Visium spot pitch (``dx_um``, default 100 um).
3. Place a small bag of discrete cells at each spot, apportioned by that spot's
   fractions (dominant + minor types), jittered inside the spot footprint.
4. Build continuum fields directly from the real fractions via a KDE onto the
   grid: ``vascular_density`` from the ENDOTHELIAL fraction, ``ecm_density`` from
   the STROMAL fraction.
5. Assemble + validate into a :class:`Microenvironment`.

Cell ``sensitivity`` is set to a neutral 1.0 here: reconstruction captures real
*composition*, not resistance biology (which the sim stage parameterises).
No synthetic data is introduced at any step.
"""
from __future__ import annotations

from typing import Dict, List

import numpy as np
from scipy.ndimage import gaussian_filter

from ..core.domain import Domain
from ..core.microenvironment import CellPopulation, CellState, CellType, TissueFields
from .assemble import assemble_microenvironment
from .deconvolve import deconvolve_sample

# Cells instantiated per spot (apportioned across types by the spot fractions).
_CELLS_PER_SPOT: int = 12
# KDE bandwidth (um) for scattering per-spot fractions onto the continuum grid.
_KDE_BANDWIDTH_UM: float = 100.0


def _normalize01(arr: np.ndarray) -> np.ndarray:
    lo = float(arr.min())
    hi = float(arr.max())
    if hi - lo <= 0.0:
        return np.zeros_like(arr)
    return (arr - lo) / (hi - lo)


def _field_from_fraction(
    domain: Domain, coords_um: np.ndarray, values: np.ndarray
) -> np.ndarray:
    """KDE a per-spot scalar (e.g. a cell-type fraction) onto the grid, [0, 1].

    Scatters each spot's value at its grid cell, Gaussian-smooths, then scales to
    [0, 1] and confines to the tissue mask.
    """
    grid = domain.zeros()
    if len(coords_um):
        col, row = domain.to_index(coords_um[:, 0], coords_um[:, 1])
        np.add.at(grid, (row, col), np.asarray(values, dtype=np.float64))
    sigma = _KDE_BANDWIDTH_UM / domain.dx
    grid = gaussian_filter(grid, sigma=sigma)
    grid = _normalize01(grid)
    return np.where(domain.tissue_mask, grid, 0.0)


def _place_cells(
    coords_um: np.ndarray,
    fractions: np.ndarray,
    type_names: List[str],
    extent: tuple,
    jitter_um: float,
    seed: int = 0,
) -> CellPopulation:
    """Instantiate discrete cells per spot, apportioned by fractions.

    Each spot contributes ``round(frac * _CELLS_PER_SPOT)`` cells of each type,
    placed at the spot centroid with uniform jitter inside a spot-sized disc and
    clamped to the domain extent.
    """
    rng = np.random.default_rng(seed)
    type_codes = [int(CellType[name]) for name in type_names]
    width, height = extent

    xs: List[float] = []
    ys: List[float] = []
    codes: List[int] = []
    for s in range(coords_um.shape[0]):
        cx, cy = coords_um[s]
        counts = np.round(fractions[s] * _CELLS_PER_SPOT).astype(int)
        for k, n_k in enumerate(counts):
            if n_k <= 0:
                continue
            ang = rng.uniform(0, 2 * np.pi, size=n_k)
            rad = jitter_um * np.sqrt(rng.uniform(0, 1, size=n_k))
            xs.extend(np.clip(cx + rad * np.cos(ang), 0.0, width))
            ys.extend(np.clip(cy + rad * np.sin(ang), 0.0, height))
            codes.extend([type_codes[k]] * n_k)

    n = len(xs)
    return CellPopulation(
        x=np.asarray(xs, dtype=np.float64),
        y=np.asarray(ys, dtype=np.float64),
        cell_type=np.asarray(codes, dtype=np.int8),
        state=np.full(n, int(CellState.VIABLE), dtype=np.int8),
        sensitivity=np.ones(n, dtype=np.float64),
    )


def reconstruct_from_visium(sample_id: str, dx_um: float = 100.0):
    """Build a real spot-level :class:`Microenvironment` from a Visium sample.

    Parameters
    ----------
    sample_id : one of ``tmesim.data.real.visium.list_samples()``.
    dx_um : grid spacing in micrometers; defaults to the Visium v1 spot pitch
        (100 um), so one grid cell corresponds to roughly one spot.

    Returns
    -------
    Microenvironment
        With discrete cells placed by real per-spot fractions, ``vascular_density``
        from the real ENDOTHELIAL fraction and ``ecm_density`` from the real
        STROMAL fraction. ``metadata`` records the type names, mean fractions and
        the shared-gene count used.
    """
    coords_um, fractions, type_names = deconvolve_sample(sample_id)

    # Domain sized to cover all spots (+1 cell margin), tissue mask = spot support.
    xmax = float(coords_um[:, 0].max())
    ymax = float(coords_um[:, 1].max())
    nx = int(np.floor(xmax / dx_um)) + 1
    ny = int(np.floor(ymax / dx_um)) + 1
    domain = Domain(nx=nx, ny=ny, dx=dx_um)

    name_to_col = {n: i for i, n in enumerate(type_names)}

    def frac_of(name: str) -> np.ndarray:
        if name in name_to_col:
            return fractions[:, name_to_col[name]]
        return np.zeros(fractions.shape[0])

    vascular_density = _field_from_fraction(
        domain, coords_um, frac_of(CellType.ENDOTHELIAL.name)
    )
    ecm_density = _field_from_fraction(
        domain, coords_um, frac_of(CellType.STROMAL.name)
    )
    fields = TissueFields(
        ecm_density=ecm_density, vascular_density=vascular_density
    )

    cells = _place_cells(
        coords_um,
        fractions,
        type_names,
        extent=domain.extent_um,
        jitter_um=min(27.5, 0.45 * dx_um),
    )

    metadata: Dict = {
        "source": "real_visium",
        "sample_id": sample_id,
        "type_names": type_names,
        "mean_fractions": {
            n: float(m) for n, m in zip(type_names, fractions.mean(axis=0))
        },
        "n_spots": int(coords_um.shape[0]),
    }
    return assemble_microenvironment(
        domain, cells, fields, patient_id=f"visium:{sample_id}", metadata=metadata
    )


if __name__ == "__main__":
    import os

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sample_id = "10x_parent_gbm"

    # --- 1. Deconvolution sanity (real data) --------------------------------
    coords_um, fractions, type_names = deconvolve_sample(sample_id)
    n_spots = fractions.shape[0]
    print(f"=== REAL Visium deconvolution: {sample_id} ===")
    print(f"spots                 : {n_spots}")
    print(f"cell types            : {type_names}")
    print("mean cell-type fractions across spots:")
    for name, m in zip(type_names, fractions.mean(axis=0)):
        print(f"  {name:12s} {m:.3f}")
    print(f"per-spot fraction sums in [{fractions.sum(1).min():.3f}, "
          f"{fractions.sum(1).max():.3f}] (expect ~1.0)")

    # --- 2. Reconstruct the real Microenvironment ----------------------------
    micro = reconstruct_from_visium(sample_id)
    print("\n=== Reconstructed Microenvironment summary ===")
    for k, v in micro.summary().items():
        print(f"  {k}: {v}")

    # --- 3. Save fraction-field heatmaps ------------------------------------
    results_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "results"
    )
    results_dir = os.path.abspath(results_dir)
    os.makedirs(results_dir, exist_ok=True)
    out_png = os.path.join(results_dir, "real_visium_fractions.png")

    name_to_col = {n: i for i, n in enumerate(type_names)}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    # (a) spot scatter coloured by TUMOR fraction
    tum = fractions[:, name_to_col.get("TUMOR", 0)]
    sc = axes[0].scatter(coords_um[:, 0], coords_um[:, 1], c=tum, s=6,
                         cmap="magma", vmin=0, vmax=1)
    axes[0].set_title("spot TUMOR fraction")
    axes[0].set_aspect("equal"); axes[0].invert_yaxis()
    fig.colorbar(sc, ax=axes[0], fraction=0.046)
    # (b) vascular_density field (from ENDOTHELIAL fraction)
    im1 = axes[1].imshow(micro.fields.vascular_density, origin="lower",
                        cmap="viridis")
    axes[1].set_title("vascular_density (ENDOTHELIAL)")
    fig.colorbar(im1, ax=axes[1], fraction=0.046)
    # (c) ecm_density field (from STROMAL fraction)
    im2 = axes[2].imshow(micro.fields.ecm_density, origin="lower", cmap="cividis")
    axes[2].set_title("ecm_density (STROMAL)")
    fig.colorbar(im2, ax=axes[2], fraction=0.046)
    fig.suptitle(f"Real Visium reconstruction: {sample_id}")
    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    print(f"\nsaved: {out_png}")
