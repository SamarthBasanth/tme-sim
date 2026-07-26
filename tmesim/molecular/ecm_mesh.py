"""Mesh-scale ECM hindrance: SIZE-dependent gating of solute diffusion by the
fibre matrix of the extracellular space.

Why this module exists
-----------------------
`brownian.py` models obstacles as ~300 nm disks (cells / large ECM aggregates).
Against those, a 0.4 nm small molecule and a 5.5 nm antibody are *both*
effectively points, so that model barely differentiates drug SIZE. But the real
size gate in tissue is the ECM FIBRE MESH: collagen / proteoglycan fibres a few
nm thick, spaced tens of nm apart. A solute's radius relative to that mesh
determines both how much of the space it can enter (equilibrium partition) and
how much drag the nearby fibres add (hydrodynamic hindrance). This module adds
that missing physics with two complementary, independent methods:

  1. `ogston_renkin_hindrance` -- fast analytic fibre-matrix theory
     (Ogston partition x Renkin/Faxen hydrodynamic drag).
  2. `brownian_mesh_hindrance` -- genuine nm-scale Brownian dynamics of a
     FINITE-radius tracer through fibre cross-sections (excluded-volume MSD).

Both take the ECM as a random matrix of cylindrical fibres of radius
`fiber_radius_nm` (~2 nm, collagen-scale) at volume fraction
`fiber_volume_fraction` (~0.05, loose interstitial matrix; higher = denser mesh).

Units: length nm (mesh scale), time s, diffusivity um^2/s (tissue convention;
1 um^2 = 1e6 nm^2). Deterministic given `seed`. No external dependencies beyond
numpy. Does NOT import or modify any existing module.

References
----------
* Ogston, A.G. (1958) "The spaces in a uniform random suspension of fibres",
  Trans. Faraday Soc. 54, 1754-1757.  -> partition / available-volume law.
* Renkin, E.M. (1954) "Filtration, diffusion, and molecular sieving through
  porous cellulose membranes", J. Gen. Physiol. 38, 225-243. -> hindered
  diffusion factor for a sphere in a pore (steric x hydrodynamic).
* Faxen drag polynomial as compiled by Deen, W.M. (1987) AIChE J. 33, 1409-1425
  ("Hindered transport of large molecules in liquid-filled pores").
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

import numpy as np

# nm^2 per um^2, for converting diffusivities between the mesh (nm) working scale
# and the tissue (um^2/s) reporting convention.
_NM2_PER_UM2 = 1.0e6


# ---------------------------------------------------------------------------
# Geometry helper: characteristic mesh pore radius of a random fibre matrix.
# ---------------------------------------------------------------------------
def mesh_pore_radius_nm(fiber_radius_nm: float, fiber_volume_fraction: float) -> float:
    """Characteristic pore radius of a random fibre matrix, in nm.

    Defined self-consistently from the Ogston partition law (below): the pore
    radius r_pore is the SOLUTE radius at which the Ogston available-volume
    fraction falls to 1/e. Setting

        exp(-phi * (1 + r_pore/r_f)^2) = 1/e
        => phi * (1 + r_pore/r_f)^2 = 1
        => r_pore = r_f * (1/sqrt(phi) - 1).

    This ties the effective "pore" directly to the same fibre geometry the
    partition uses, and gives physically sane tens-of-nm pores for loose matrix
    (e.g. r_f=2 nm, phi=0.05 -> r_pore ~ 6.9 nm; mesh diameter ~14 nm).

    Clamped to a small positive floor so lambda = r_s/r_pore stays finite.
    """
    if fiber_radius_nm <= 0:
        raise ValueError(f"fiber_radius_nm must be positive, got {fiber_radius_nm!r}")
    if not (0.0 < fiber_volume_fraction < 1.0):
        raise ValueError(
            f"fiber_volume_fraction must be in (0,1), got {fiber_volume_fraction!r}"
        )
    r_pore = fiber_radius_nm * (1.0 / math.sqrt(fiber_volume_fraction) - 1.0)
    return max(r_pore, 1e-6)


def _faxen_drag(lam: float) -> float:
    """Renkin/Faxen HYDRODYNAMIC drag factor Kd(lambda) for a sphere in a pore.

    lambda = solute_radius / pore_radius. This is the drag (enhanced-friction)
    part of Renkin's (1954) hindered-diffusion expression, as compiled by
    Deen (1987):

        Kd(lambda) = 1 - 2.104 lambda + 2.09 lambda^3 - 0.95 lambda^5,   0<=lambda<1

    It multiplies the *equilibrium partition* (supplied here by Ogston, not by
    Renkin's own steric (1-lambda)^2 term, to avoid double-counting sterics).
    For lambda >= 1 the solute cannot pass and Kd is taken as 0. Clamped to
    [0,1]. Monotonically decreasing in lambda on [0,1).
    """
    if lam >= 1.0:
        return 0.0
    if lam <= 0.0:
        return 1.0
    kd = 1.0 - 2.104 * lam + 2.09 * lam ** 3 - 0.95 * lam ** 5
    return float(min(max(kd, 0.0), 1.0))


def ogston_partition(
    solute_radius_nm: float,
    fiber_radius_nm: float = 2.0,
    fiber_volume_fraction: float = 0.05,
) -> float:
    """Ogston (1958) available-volume fraction (equilibrium partition) in [0,1].

        phi_avail = exp( -phi_f * ((r_s + r_f)/r_f)^2 )
                  = exp( -phi_f * (1 + r_s/r_f)^2 )

    Physical picture: fraction of the fibre-matrix void space whose centre a
    sphere of radius r_s can occupy without its surface overlapping a fibre of
    radius r_f. Monotonically DECREASING in both r_s and phi_f; -> 1 as r_s->0
    and phi_f->0; -> 0 for a large solute in a dense mesh (strong size sieving).
    """
    if solute_radius_nm < 0:
        raise ValueError(f"solute_radius_nm must be >= 0, got {solute_radius_nm!r}")
    ratio = (solute_radius_nm + fiber_radius_nm) / fiber_radius_nm
    return float(math.exp(-fiber_volume_fraction * ratio ** 2))


def ogston_renkin_hindrance(
    solute_radius_nm: float,
    fiber_radius_nm: float = 2.0,
    fiber_volume_fraction: float = 0.05,
) -> Dict[str, float]:
    """Analytic SIZE-dependent ECM hindrance: Ogston partition x Renkin drag.

    Combines two independent, established relations for a sphere in a random
    fibre matrix:

    * partition (Ogston 1958): equilibrium available-volume fraction
          phi_avail = exp(-phi_f (1 + r_s/r_f)^2).
    * hindrance_diffusion (Renkin 1954 / Faxen drag, Deen 1987): hydrodynamic
          slowdown from fibre walls, Kd(lambda) with lambda = r_s/r_pore, where
          r_pore = mesh_pore_radius_nm(r_f, phi_f).

    effective_hindrance = partition * hindrance_diffusion  (in [0,1]): the total
    factor by which the free diffusivity is reduced by the mesh (both fewer
    accessible sites AND more drag). Larger solutes -> both factors collapse, so
    effective hindrance drops steeply with size.

    Returns
    -------
    dict: {partition, hindrance_diffusion, effective_hindrance, mesh_pore_nm,
           lambda, solute_radius_nm}.
    """
    r_pore = mesh_pore_radius_nm(fiber_radius_nm, fiber_volume_fraction)
    lam = solute_radius_nm / r_pore
    partition = ogston_partition(solute_radius_nm, fiber_radius_nm, fiber_volume_fraction)
    kd = _faxen_drag(lam)
    return {
        "solute_radius_nm": float(solute_radius_nm),
        "partition": partition,
        "hindrance_diffusion": kd,
        "effective_hindrance": partition * kd,
        "mesh_pore_nm": r_pore,
        "lambda": float(lam),
    }


# ---------------------------------------------------------------------------
# Method 2: genuine nm-scale Brownian dynamics through fibre cross-sections.
# ---------------------------------------------------------------------------
def _linfit(x: np.ndarray, y: np.ndarray):
    """OLS line fit -> (slope, intercept, R^2)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = slope * x + intercept
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(slope), float(intercept), float(r2)


def brownian_mesh_hindrance(
    solute_radius_nm: float,
    fiber_radius_nm: float = 2.0,
    fiber_volume_fraction: float = 0.05,
    free_D_um2s: float = 300.0,
    box_nm: float = 80.0,
    n_particles: int = 400,
    steps: int = 4000,
    dt_s: Optional[float] = None,
    seed: int = 0,
) -> Dict:
    """Overdamped Brownian dynamics of a FINITE-radius tracer through an ECM
    fibre mesh, at the NM scale -- the size-resolving counterpart of brownian.py.

    Geometry
    --------
    The ECM fibres are cylinders of radius `fiber_radius_nm`. We simulate the 2D
    plane PERPENDICULAR to the fibres, in which each fibre appears as a disk of
    radius r_f. For parallel cylinders the disk area fraction in this plane
    equals the fibre VOLUME fraction, so the number of fibre cross-sections in a
    periodic box of side L is

        n_fibres = round( phi_f * L^2 / (pi r_f^2) ).

    The tracer is a sphere of radius r_s. Its centre is excluded from within
    (r_f + r_s) of any fibre centre -- this FINITE tracer radius is what makes the
    hindrance size-dependent: a bigger solute sees fatter obstacles and a smaller
    accessible space. Fibre centres are Poisson-placed (i.i.d. uniform; overlaps
    allowed), matching the random-matrix assumption behind the Ogston law.

    Dynamics (Euler-Maruyama, overdamped Langevin)
    ----------------------------------------------
    Per-axis step std = sqrt(2 free_D dt). A proposed move whose (periodic-image)
    centre lands within (r_f + r_s) of any fibre is REJECTED (tracer stays put) --
    the standard Monte-Carlo scheme for obstructed diffusion (Saxton 1994).
    Effective_D is read from the long-lag slope of the time-origin-averaged MSD of
    UNWRAPPED positions: MSD(t) = 4 D_eff t in 2D. hindrance = D_eff / free_D.

    Working units are nm and s (free_D converted from um^2/s via 1 um^2 = 1e6 nm^2);
    the returned diffusivities are converted back to um^2/s.

    Returns
    -------
    dict: {effective_D_um2s, hindrance, free_D_um2s, msd_nm2, times_s,
           mesh_pore_nm, excluded_radius_nm, n_fibres, accessible_fraction,
           fit_r2, dt_s, solute_radius_nm}.
    """
    if solute_radius_nm < 0:
        raise ValueError(f"solute_radius_nm must be >= 0, got {solute_radius_nm!r}")
    rng = np.random.default_rng(seed)

    free_D_nm2s = free_D_um2s * _NM2_PER_UM2  # um^2/s -> nm^2/s
    excl = fiber_radius_nm + solute_radius_nm  # excluded centre-distance
    area = box_nm ** 2

    # number of fibre cross-sections from the (in-plane) area = volume fraction
    n_fib = int(round(fiber_volume_fraction * area / (math.pi * fiber_radius_nm ** 2)))
    n_fib = max(n_fib, 1)
    fibers = rng.uniform(0.0, box_nm, size=(n_fib, 2))

    # time step: keep rms per-axis step small vs the fibre radius so obstacles are
    # well resolved and cannot be tunnelled through.
    if dt_s is None:
        step_std_target = 0.04 * fiber_radius_nm  # nm
        dt = step_std_target ** 2 / (2.0 * free_D_nm2s)
    else:
        dt = float(dt_s)
    step_std = math.sqrt(2.0 * free_D_nm2s * dt)  # nm

    half = box_nm / 2.0

    def _inside(points: np.ndarray) -> np.ndarray:
        """Boolean mask: is each (wrapped) point within excl of any fibre?
        Minimum-image periodic distance. points: (N,2) -> (N,) bool."""
        dx = points[:, None, 0] - fibers[None, :, 0]
        dy = points[:, None, 1] - fibers[None, :, 1]
        dx -= box_nm * np.round(dx / box_nm)
        dy -= box_nm * np.round(dy / box_nm)
        d2 = dx * dx + dy * dy
        return np.any(d2 < excl * excl, axis=1)

    # Estimate accessible-area fraction (Monte-Carlo) and place tracers there.
    probe = rng.uniform(0.0, box_nm, size=(4000, 2))
    accessible_fraction = float(np.mean(~_inside(probe)))

    pos = rng.uniform(0.0, box_nm, size=(n_particles, 2))
    for _ in range(400):
        bad = _inside(pos)
        if not bad.any():
            break
        pos[bad] = rng.uniform(0.0, box_nm, size=(int(bad.sum()), 2))
    unwrapped = pos.copy()

    record_every = max(1, steps // 200)
    frames = [unwrapped.copy()]
    for s in range(1, steps + 1):
        delta = rng.normal(0.0, step_std, size=(n_particles, 2))
        prop = unwrapped + delta
        pw = prop - box_nm * np.floor(prop / box_nm)  # wrap into [0, box)
        accept = ~_inside(pw)
        unwrapped[accept] = prop[accept]
        if s % record_every == 0:
            frames.append(unwrapped.copy())

    traj = np.stack(frames)  # (n_frames, n_particles, 2), unwrapped, nm
    n_frames = traj.shape[0]
    frame_dt = record_every * dt
    lags = np.arange(1, n_frames)
    msd = np.empty(lags.size, float)
    for i, lag in enumerate(lags):
        diff = traj[lag:] - traj[:-lag]
        msd[i] = float(np.mean(np.sum(diff ** 2, axis=2)))
    times = lags * frame_dt

    lo = len(times) // 4
    hi = (3 * len(times)) // 4
    slope, _b, r2 = _linfit(times[lo:hi], msd[lo:hi])  # nm^2/s
    eff_D_nm2s = max(slope / 4.0, 0.0)
    eff_D_um2s = eff_D_nm2s / _NM2_PER_UM2

    return {
        "solute_radius_nm": float(solute_radius_nm),
        "effective_D_um2s": eff_D_um2s,
        "hindrance": eff_D_um2s / free_D_um2s if free_D_um2s > 0 else 0.0,
        "free_D_um2s": float(free_D_um2s),
        "msd_nm2": msd,
        "times_s": times,
        "mesh_pore_nm": mesh_pore_radius_nm(fiber_radius_nm, fiber_volume_fraction),
        "excluded_radius_nm": float(excl),
        "n_fibres": n_fib,
        "accessible_fraction": accessible_fraction,
        "fit_r2": r2,
        "dt_s": dt,
    }


def hindrance_by_drug_size(
    solute_radii_nm: List[float],
    fiber_radius_nm: float = 2.0,
    fiber_volume_fraction: float = 0.05,
    free_D_um2s: float = 300.0,
    brownian_kw: Optional[dict] = None,
    seed: int = 0,
) -> Dict:
    """Map several solute radii to effective hindrance by BOTH methods.

    Returns
    -------
    dict: {
      'fiber_radius_nm', 'fiber_volume_fraction', 'mesh_pore_nm',
      'rows': [ {solute_radius_nm, ogston_partition, ogston_hindrance_diffusion,
                 ogston_hindrance, brownian_hindrance, brownian_eff_D_um2s},
                ... one per radius ]
    }
    ogston_hindrance is the effective (partition x drag) analytic factor;
    brownian_hindrance is the BD effective_D/free_D. Both are in [0,1] and drop
    with solute size.
    """
    bkw = dict(brownian_kw or {})
    rows = []
    for i, r_s in enumerate(solute_radii_nm):
        an = ogston_renkin_hindrance(r_s, fiber_radius_nm, fiber_volume_fraction)
        bd = brownian_mesh_hindrance(
            r_s,
            fiber_radius_nm=fiber_radius_nm,
            fiber_volume_fraction=fiber_volume_fraction,
            free_D_um2s=free_D_um2s,
            seed=seed + i,
            **bkw,
        )
        rows.append({
            "solute_radius_nm": float(r_s),
            "ogston_partition": an["partition"],
            "ogston_hindrance_diffusion": an["hindrance_diffusion"],
            "ogston_hindrance": an["effective_hindrance"],
            "brownian_hindrance": bd["hindrance"],
            "brownian_eff_D_um2s": bd["effective_D_um2s"],
            "brownian_accessible_fraction": bd["accessible_fraction"],
        })
    return {
        "fiber_radius_nm": fiber_radius_nm,
        "fiber_volume_fraction": fiber_volume_fraction,
        "mesh_pore_nm": mesh_pore_radius_nm(fiber_radius_nm, fiber_volume_fraction),
        "rows": rows,
    }


if __name__ == "__main__":
    import time as _time

    t0 = _time.time()
    R_F = 2.0
    PHI = 0.05
    FREE_D = 300.0

    print("=" * 74)
    print("ECM mesh-scale hindrance self-test")
    print(f"  fibre radius = {R_F} nm, fibre volume fraction = {PHI}, "
          f"free_D = {FREE_D} um^2/s")
    print(f"  mesh pore radius = {mesh_pore_radius_nm(R_F, PHI):.2f} nm "
          f"(diameter ~ {2*mesh_pore_radius_nm(R_F, PHI):.1f} nm)")
    print("=" * 74)

    SMALL = 0.4   # nm, typical small-molecule drug
    ANTIBODY = 5.5  # nm, IgG hydrodynamic radius

    # ---- [1] Ogston partition properties: in [0,1], monotone decreasing ------
    radii = [0.2, 0.4, 0.8, 1.5, 3.0, 5.5, 8.0]
    parts = [ogston_partition(r, R_F, PHI) for r in radii]
    print("\n[1] Ogston partition vs solute radius:")
    for r, p in zip(radii, parts):
        print(f"    r_s = {r:5.2f} nm  ->  partition = {p:.4f}")
    assert all(0.0 <= p <= 1.0 for p in parts), "partition out of [0,1]"
    assert all(parts[i] > parts[i + 1] for i in range(len(parts) - 1)), \
        "partition must decrease with solute radius"
    # monotone decreasing in fibre volume fraction (fix r_s = ANTIBODY)
    phis = [0.02, 0.05, 0.10, 0.20, 0.30]
    parts_phi = [ogston_partition(ANTIBODY, R_F, f) for f in phis]
    print("\n    Ogston partition of antibody (5.5 nm) vs fibre volume fraction:")
    for f, p in zip(phis, parts_phi):
        print(f"    phi_f = {f:4.2f}  ->  partition = {p:.4f}")
    assert all(parts_phi[i] > parts_phi[i + 1] for i in range(len(parts_phi) - 1)), \
        "partition must decrease with fibre volume fraction"

    # ---- [2] Antibody strongly excluded from a tight mesh -------------------
    part_ab_tight = ogston_partition(ANTIBODY, R_F, 0.30)
    print(f"\n[2] Antibody (5.5 nm) partition in tight mesh (phi_f=0.30): "
          f"{part_ab_tight:.4f}")
    assert part_ab_tight < 0.05, \
        f"antibody not strongly excluded from tight mesh: {part_ab_tight}"

    # ---- [3] Size table: BOTH methods, loose interstitium (phi_f=0.05) ------
    table_radii = [SMALL, 1.0, 2.0, 3.5, ANTIBODY]
    res = hindrance_by_drug_size(
        table_radii, fiber_radius_nm=R_F, fiber_volume_fraction=PHI,
        free_D_um2s=FREE_D,
        brownian_kw=dict(n_particles=400, steps=4000, box_nm=80.0),
        seed=7,
    )
    print(f"\n[3] Hindrance vs solute size, loose matrix (phi_f={PHI}):")
    print(f"    {'r_s(nm)':>8} {'ogston_part':>12} {'ogston_hind':>12} "
          f"{'brown_hind':>11} {'brown_D':>9} {'brown_acc':>10}")
    for row in res["rows"]:
        print(f"    {row['solute_radius_nm']:8.2f} {row['ogston_partition']:12.4f} "
              f"{row['ogston_hindrance']:12.4f} {row['brownian_hindrance']:11.4f} "
              f"{row['brownian_eff_D_um2s']:9.2f} "
              f"{row['brownian_accessible_fraction']:10.3f}")

    og = {row["solute_radius_nm"]: row["ogston_hindrance"] for row in res["rows"]}
    bd = {row["solute_radius_nm"]: row["brownian_hindrance"] for row in res["rows"]}

    og_ratio = og[SMALL] / max(og[ANTIBODY], 1e-12)
    bd_drop = bd[SMALL] / max(bd[ANTIBODY], 1e-12)
    print(f"\n    small/antibody ratio  ->  Ogston: {og_ratio:.1f}x   "
          f"Brownian(tortuosity): {bd_drop:.2f}x")

    # both methods AGREE IN TREND: hindrance monotonically decreasing with size
    og_seq = [og[r] for r in table_radii]
    bd_seq = [bd[r] for r in table_radii]
    assert all(og_seq[i] >= og_seq[i + 1] - 1e-9 for i in range(len(og_seq) - 1)), \
        f"Ogston hindrance not decreasing with size: {og_seq}"
    assert all(bd_seq[i] >= bd_seq[i + 1] - 5e-3 for i in range(len(bd_seq) - 1)), \
        f"Brownian hindrance not decreasing with size: {bd_seq}"
    # STRONG analytic size gate (partition x near-wall drag): antibody << small
    assert og_ratio >= 5.0, f"Ogston size dependence too weak: {og_ratio:.1f}x"
    # BD (pure excluded-volume tortuosity) is milder in a LOOSE mesh but still drops
    assert bd[ANTIBODY] < bd[SMALL], "Brownian hindrance must drop for the antibody"

    # ---- [4] Dense tumour-ECM mesh: the geometric BD size gate turns ON ------
    # Pure excluded-volume tortuosity only gates size strongly once the mesh is
    # dense enough that the large solute's accessible space approaches percolation.
    PHI_DENSE = 0.24
    small_bd = brownian_mesh_hindrance(
        SMALL, fiber_radius_nm=R_F, fiber_volume_fraction=PHI_DENSE,
        free_D_um2s=FREE_D, n_particles=400, steps=4000, box_nm=80.0, seed=21)
    ab_bd = brownian_mesh_hindrance(
        ANTIBODY, fiber_radius_nm=R_F, fiber_volume_fraction=PHI_DENSE,
        free_D_um2s=FREE_D, n_particles=400, steps=4000, box_nm=80.0, seed=22)
    bd_dense_ratio = small_bd["hindrance"] / max(ab_bd["hindrance"], 1e-12)
    print(f"\n[4] Dense tumour-ECM mesh (phi_f={PHI_DENSE}), Brownian dynamics:")
    print(f"    small (0.4 nm) : hindrance={small_bd['hindrance']:.3f}  "
          f"eff_D={small_bd['effective_D_um2s']:.1f} um^2/s  "
          f"accessible={small_bd['accessible_fraction']:.3f}")
    print(f"    antibody(5.5nm): hindrance={ab_bd['hindrance']:.3f}  "
          f"eff_D={ab_bd['effective_D_um2s']:.1f} um^2/s  "
          f"accessible={ab_bd['accessible_fraction']:.3f}")
    print(f"    small/antibody Brownian ratio = {bd_dense_ratio:.1f}x  "
          f"(fit R^2 small={small_bd['fit_r2']:.3f}, ab={ab_bd['fit_r2']:.3f})")
    assert bd_dense_ratio >= 2.0, \
        f"Dense-mesh Brownian size gate too weak: {bd_dense_ratio:.1f}x"
    assert ab_bd["accessible_fraction"] < 0.10, \
        f"antibody not near-excluded from dense mesh: {ab_bd['accessible_fraction']}"

    print("\nSummary: SIZE strongly gates ECM diffusion. Analytic Ogston-Renkin gives")
    print(f"  a {og_ratio:.0f}x small-molecule/antibody hindrance gap even in loose matrix;")
    print("  genuine nm-scale excluded-volume BD confirms the same trend and sharpens")
    print(f"  to {bd_dense_ratio:.1f}x with near-total antibody exclusion in a dense mesh.")
    print(f"\nAll asserts passed. Runtime {_time.time() - t0:.2f} s")
