"""Synthetic vascular-network generation for the tumor microenvironment.

The vasculature is the drug's delivery highway and the tissue's oxygen supply.
For the Stage 0 kill-test we need a *spatially structured* vascular-density field
rather than a smooth blob, because the whole point of a mechanistic model is that
delivery is heterogeneous: tortuous vessels seed the periphery and thin out toward
a poorly-perfused tumor core.

This module grows branching, tortuous vessel trees as random walks anchored at the
domain periphery, rasterizes them with thickness onto the grid, and gaussian-blurs
the result into a perfused-vessel *fraction* per cell, normalized to ``[0, 1]``.

Biological anchor: solid tumor cores are abnormally and sparsely vascularized
(chaotic angiogenesis, vessel collapse under interstitial pressure). We therefore
*suppress* vessel growth and density inside ``tumor_mask`` so the downstream oxygen
and transport solvers see the hypoxic, drug-starved core that real tumors present.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter

from ..core.domain import Domain


def build_vascular_network(
    domain: Domain,
    seed: int = 0,
    tumor_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Build a tortuous vascular-density field on ``domain``.

    Seeds vessel roots at the domain periphery and grows branching random walks
    inward, rasterizes them with vessel thickness, and gaussian-blurs the mask
    into a smooth perfused-vessel fraction in ``[0, 1]``. Growth is actively
    suppressed inside ``tumor_mask`` (and its density down-weighted afterward) so
    the tumor core is sparsely vascularized, driving downstream hypoxia and poor
    drug delivery.

    Parameters
    ----------
    domain : Domain
        The (ny, nx) grid to build the field on.
    seed : int, default 0
        Seed for ``numpy.random.default_rng`` — fully determines the output.
    tumor_mask : np.ndarray of bool, shape (ny, nx), optional
        True inside the tumor. Vessels are grown sparser here. If None, a roughly
        uniform tortuous network is produced across the whole tissue.

    Returns
    -------
    np.ndarray
        A ``(ny, nx)`` float64 vascular-density field, values in ``[0, 1]``,
        where the value is the perfused-vessel fraction of that grid cell.
    """
    rng = np.random.default_rng(seed)
    ny, nx = domain.shape

    # ---- tumor suppression map ------------------------------------------------
    # A soft [0,1] "growth permissiveness" field: ~1 in healthy tissue, small in
    # the tumor core, with a graded margin so vessels thin out rather than stop
    # at a hard edge. distance_transform_edt gives distance-into-tumor in cells.
    if tumor_mask is not None:
        tumor_mask = np.asarray(tumor_mask, dtype=bool)
        if tumor_mask.shape != (ny, nx):
            raise ValueError(
                f"tumor_mask shape {tumor_mask.shape} != {(ny, nx)}"
            )
        depth = distance_transform_edt(tumor_mask)  # cells inside the tumor
        margin = max(3.0, 0.10 * min(ny, nx))  # graded margin width in cells
        # 1.0 outside tumor -> decays toward core_perfusion deep inside.
        core_perfusion = 0.12
        permissive = core_perfusion + (1.0 - core_perfusion) * np.exp(-depth / margin)
    else:
        permissive = np.ones((ny, nx), dtype=np.float64)

    # ---- grow branching vessel trees -----------------------------------------
    thickness = domain.zeros()  # accumulates vessel "mass" before blurring
    n_roots = max(6, int(round(0.9 * np.sqrt(nx * ny) / 4)))
    roots = _seed_periphery_roots(n_roots, nx, ny, rng)

    max_len = int(1.5 * (nx + ny))  # generous per-vessel step budget
    for (rx, ry, dirx, diry) in roots:
        _grow_vessel(
            thickness,
            permissive,
            x=rx,
            y=ry,
            dx=dirx,
            dy=diry,
            steps_left=max_len,
            radius=1.6,
            depth=0,
            rng=rng,
        )

    # Chaotic, abnormal intratumoral micro-vessels: real tumor cores are not
    # perfectly avascular but carry sparse, disorganized vessels. Seed a few short
    # tortuous fragments *inside* the tumor so the core is sparse-but-present.
    if tumor_mask is not None and tumor_mask.any():
        tumor_rows, tumor_cols = np.nonzero(tumor_mask)
        n_abnormal = max(4, int(round(0.02 * tumor_mask.sum() ** 0.5 * 6)))
        pick = rng.integers(0, tumor_rows.size, size=n_abnormal)
        for p in pick:
            ang = rng.uniform(0.0, 2.0 * np.pi)
            _grow_vessel(
                thickness,
                permissive,
                x=float(tumor_cols[p]),
                y=float(tumor_rows[p]),
                dx=float(np.cos(ang)),
                dy=float(np.sin(ang)),
                steps_left=int(0.15 * (nx + ny)),
                radius=1.0,
                depth=3,  # near max_depth: minimal branching, stays sparse
                rng=rng,
            )

    # ---- rasterized mass -> smooth density ------------------------------------
    density = gaussian_filter(thickness, sigma=max(1.0, 0.012 * min(ny, nx)))

    # Re-apply tumor suppression on the density itself (blurring bleeds healthy
    # vessels into the core; this keeps the core genuinely sparse).
    density *= permissive

    # Respect the tissue mask (no vessels outside simulated tissue).
    density *= domain.tissue_mask.astype(np.float64)

    # ---- normalize to [0, 1] --------------------------------------------------
    dmax = float(density.max())
    if dmax > 0.0:
        density = density / dmax
    return np.clip(density, 0.0, 1.0)


def _seed_periphery_roots(
    n_roots: int, nx: int, ny: int, rng: np.random.Generator
) -> List[Tuple[float, float, float, float]]:
    """Place ``n_roots`` vessel roots on the domain border, aimed inward.

    Returns a list of ``(x, y, dx, dy)`` tuples with a unit-ish inward direction.
    """
    roots: List[Tuple[float, float, float, float]] = []
    sides = rng.integers(0, 4, size=n_roots)
    for side in sides:
        if side == 0:  # left edge, aim right
            x, y, dx, dy = 0.5, rng.uniform(0, ny), 1.0, rng.uniform(-0.4, 0.4)
        elif side == 1:  # right edge, aim left
            x, y, dx, dy = nx - 1.5, rng.uniform(0, ny), -1.0, rng.uniform(-0.4, 0.4)
        elif side == 2:  # top edge, aim down
            x, y, dx, dy = rng.uniform(0, nx), 0.5, rng.uniform(-0.4, 0.4), 1.0
        else:  # bottom edge, aim up
            x, y, dx, dy = rng.uniform(0, nx), ny - 1.5, rng.uniform(-0.4, 0.4), -1.0
        norm = np.hypot(dx, dy) + 1e-12
        roots.append((float(x), float(y), float(dx / norm), float(dy / norm)))
    return roots


def _grow_vessel(
    thickness: np.ndarray,
    permissive: np.ndarray,
    x: float,
    y: float,
    dx: float,
    dy: float,
    steps_left: int,
    radius: float,
    depth: int,
    rng: np.random.Generator,
) -> None:
    """Grow one tortuous vessel segment via a persistent random walk, recursing at
    branch points. Deposits thickness into ``thickness`` (in place). Vessels
    terminate early with higher probability where ``permissive`` is low (tumor).
    """
    ny, nx = thickness.shape
    heading = np.arctan2(dy, dx)
    max_depth = 5

    for _ in range(steps_left):
        ix, iy = int(round(x)), int(round(y))
        if not (0 <= ix < nx and 0 <= iy < ny):
            return  # walked off the grid

        perm = permissive[iy, ix]

        # In poorly-perfused (tumor) tissue vessels are likely to prune/terminate.
        if rng.random() < 0.05 + 0.30 * (1.0 - perm):
            return

        # Deposit a thickness disk; thinner (lower mass) in low-perfusion tissue.
        _stamp_disk(thickness, x, y, radius, weight=0.5 + 0.5 * perm)

        # Persistent random walk: tortuosity via correlated heading jitter.
        heading += rng.normal(0.0, 0.35)
        x += np.cos(heading)
        y += np.sin(heading)

        # Taper the vessel as it travels (arteriole -> capillary).
        radius = max(0.7, radius * 0.998)

        # Branch: spawn a child vessel that diverges from the parent heading.
        if depth < max_depth and rng.random() < 0.03 * perm:
            branch_heading = heading + rng.choice([-1.0, 1.0]) * rng.uniform(0.5, 1.1)
            _grow_vessel(
                thickness,
                permissive,
                x=x,
                y=y,
                dx=float(np.cos(branch_heading)),
                dy=float(np.sin(branch_heading)),
                steps_left=int(steps_left * 0.5),
                radius=max(0.7, radius * 0.8),
                depth=depth + 1,
                rng=rng,
            )


def _stamp_disk(
    field: np.ndarray, cx: float, cy: float, radius: float, weight: float
) -> None:
    """Add a soft disk of the given ``radius`` and ``weight`` to ``field`` at
    ``(cx, cy)`` (in place), vectorized over the local bounding box."""
    ny, nx = field.shape
    r = int(np.ceil(radius))
    x0, x1 = max(0, int(cx) - r), min(nx, int(cx) + r + 1)
    y0, y1 = max(0, int(cy) - r), min(ny, int(cy) + r + 1)
    if x0 >= x1 or y0 >= y1:
        return
    ys, xs = np.mgrid[y0:y1, x0:x1]
    dist2 = (xs - cx) ** 2 + (ys - cy) ** 2
    field[y0:y1, x0:x1] += weight * np.exp(-dist2 / (2.0 * radius ** 2))


if __name__ == "__main__":
    # Self-test: vessels must be SPARSER inside the tumor core than in the margin.
    dom = Domain(nx=120, ny=120, dx=15.0)
    X, Y = dom.coords_um()
    cx_um, cy_um = dom.extent_um[0] / 2.0, dom.extent_um[1] / 2.0
    core_radius = 0.22 * min(dom.extent_um)
    tumor_mask = (X - cx_um) ** 2 + (Y - cy_um) ** 2 <= core_radius ** 2

    field = build_vascular_network(dom, seed=0, tumor_mask=tumor_mask)

    # Sanity on the field itself.
    assert field.shape == dom.shape, "field shape must match domain"
    assert field.min() >= 0.0 and field.max() <= 1.0, "field must be in [0,1]"
    assert np.isclose(field.max(), 1.0), "field should be normalized to reach 1.0"

    # Compare core vs. a margin annulus just outside the tumor.
    dist_um = np.sqrt((X - cx_um) ** 2 + (Y - cy_um) ** 2)
    margin_mask = (dist_um > core_radius) & (dist_um <= 1.8 * core_radius)
    core_mean = float(field[tumor_mask].mean())
    margin_mean = float(field[margin_mask].mean())

    print(f"grid                 : {dom.shape}, dx={dom.dx} um")
    print(f"core mean density    : {core_mean:.4f}")
    print(f"margin mean density  : {margin_mean:.4f}")
    print(f"whole-field mean      : {field.mean():.4f}")
    print(f"core / margin ratio  : {core_mean / margin_mean:.3f}")

    assert core_mean < margin_mean, (
        f"tumor core ({core_mean:.4f}) must be sparser than margin ({margin_mean:.4f})"
    )
    print("SELF-TEST PASSED: tumor core is more sparsely vascularized than the margin.")
