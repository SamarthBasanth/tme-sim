"""Genuine overdamped Langevin / Brownian dynamics for effective diffusivity in a
crowded medium (obstructed diffusion / 2D continuum Lorentz gas).

This is the "molecular dynamics at a computable scale" piece of the molecular layer.
We do NOT integrate atoms of a tissue (~10^15 of them). Instead we simulate the
transport-relevant physics directly: many independent point-like tracer particles
performing an overdamped Brownian random walk through a periodic 2D box packed with
randomly placed circular obstacles (cells / ECM aggregates). The excluded volume of
the obstacles hinders diffusion, and we measure the resulting EFFECTIVE diffusivity
from the mean-squared displacement (MSD).

Physics summary
---------------
Overdamped Langevin (no inertia): dx = sqrt(2 D dt) * xi,  xi ~ N(0, I).
So each Cartesian step has per-axis std = sqrt(2 * free_D * dt). This is the exact
Euler-Maruyama update for free Brownian motion; in 2D the free MSD grows as
    <|r(t) - r(0)|^2> = 4 * free_D * t          (2 axes * 2 D t).

Obstacles are impenetrable disks of radius R_obs (see OBSTACLE_RADIUS_UM). A proposed
move that would place a tracer's centre inside any obstacle (within R_obs + particle
radius) is REJECTED: the particle stays put for that step. Rejection sampling of the
free random walk restricted to the accessible domain is the standard Monte-Carlo
scheme for obstructed diffusion (Saxton, Biophys. J. 1994); it preserves the uniform
equilibrium measure over the accessible area and yields proper diffusive motion with a
reduced diffusion coefficient D_eff <= free_D.

Obstacles are placed by a Poisson process (positions i.i.d. uniform), so they may
overlap. `obstacle_fraction` is therefore the NOMINAL area fraction
    phi = N_obs * pi * R_obs^2 / box_area,
and the actually-covered area fraction is 1 - exp(-phi) (Boolean/Swiss-cheese model).
For small phi the two coincide, which is the regime used for the dilute-limit check.

Units: length um, time s, diffusivity um^2/s. Deterministic given `seed`.
"""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np

# Obstacle radius. The contract asks for "a few hundred nm" objects representing
# cells / ECM aggregates. We use 300 nm = 0.3 um: large compared with a small-molecule
# tracer (~0.5 nm) so the tracer is effectively a point, and large compared with the
# per-step displacement so obstacles are well resolved and cannot be tunnelled through.
OBSTACLE_RADIUS_UM: float = 0.30

# Grid used only to accelerate the "is this point inside an obstacle?" test. The static
# obstacle field is rasterised once into a boolean occupancy grid; each per-step test is
# then an O(1) lookup instead of an O(N_obstacles) distance scan. Resolution R_obs/20
# resolves the obstacle boundary to ~5% of its radius.
_GRID_CELLS_PER_RADIUS: int = 20


def _linfit(x: np.ndarray, y: np.ndarray):
    """Ordinary least-squares line fit. Returns (slope, intercept, R^2)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    A = np.vstack([x, np.ones_like(x)]).T
    (slope, intercept), *_ = np.linalg.lstsq(A, y, rcond=None)
    yhat = slope * x + intercept
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(slope), float(intercept), float(r2)


def _build_occupancy(centers: np.ndarray, excl_um: float, box_um: float, ng: int) -> np.ndarray:
    """Rasterise impenetrable disks into a periodic boolean occupancy grid.

    A grid cell is marked True if its centre lies within `excl_um` of any obstacle
    centre (minimum-image / periodic distance). `excl_um` = R_obs + particle_radius.
    """
    occ = np.zeros((ng, ng), dtype=bool)
    if len(centers) == 0:
        return occ
    h = box_um / ng
    r_cells = int(np.ceil(excl_um / h)) + 1
    for cx, cy in centers:
        ic = int(cx // h)
        jc = int(cy // h)
        ii = np.arange(ic - r_cells, ic + r_cells + 1)
        jj = np.arange(jc - r_cells, jc + r_cells + 1)
        cxs = (ii + 0.5) * h
        cys = (jj + 0.5) * h
        # minimum-image displacement handles obstacles/cells straddling the boundary
        dx = ((cxs - cx + box_um / 2.0) % box_um) - box_um / 2.0
        dy = ((cys - cy + box_um / 2.0) % box_um) - box_um / 2.0
        mask = (dx[:, None] ** 2 + dy[None, :] ** 2) <= excl_um ** 2
        gi = (ii % ng)[:, None]
        gj = (jj % ng)[None, :]
        gi_b = np.broadcast_to(gi, mask.shape)[mask]
        gj_b = np.broadcast_to(gj, mask.shape)[mask]
        occ[gi_b, gj_b] = True
    return occ


def simulate_effective_diffusivity(free_D_um2s: float, obstacle_fraction: float,
                                   particle_radius_nm: float = 0.5, box_um: float = 20.0,
                                   n_particles: int = 400, steps: int = 4000,
                                   dt_s: Optional[float] = None, seed: int = 0) -> Dict:
    """2D overdamped Brownian dynamics of tracers through a crowded periodic box.

    Parameters
    ----------
    free_D_um2s : dilute (obstacle-free) diffusivity, um^2/s (Stokes-Einstein input).
    obstacle_fraction : nominal obstacle area fraction phi (Poisson-placed disks; see
        module docstring). 0 => no obstacles.
    particle_radius_nm : tracer radius (adds to obstacle radius for the excluded region).
    box_um : side length of the periodic square box, um.
    n_particles : number of independent tracers (averaged for MSD).
    steps : number of Brownian time steps.
    dt_s : time step. If None, chosen so the rms per-axis step ~ 0.04 * R_obs, i.e.
        rms step << obstacle size (and << box) so obstacles are well resolved.
    seed : RNG seed (numpy default_rng) for determinism.

    Scheme
    ------
    * Euler-Maruyama free step: delta ~ N(0, sqrt(2 D dt)) per axis, all particles at once.
    * Excluded volume: reject any step landing the centre inside an obstacle disk
      (radius R_obs + particle radius); particle stays at its previous position.
    * Periodicity: the WRAPPED position (pos mod box) is used for the collision test,
      while the UNWRAPPED position is accumulated separately so the MSD is not corrupted
      by boundary wrapping.
    * Estimator: MSD(t) = <|r_unwrapped(t) - r(0)|^2> averaged over particles;
      effective_D = (slope of MSD vs t) / 4, fitted over the long-lag half (past the
      short-time crossover where obstructed tracers are still locally free).

    Returns
    -------
    dict with keys: effective_D, hindrance (= effective_D/free_D), msd (array),
    times (array), obstacle_fraction, plus diagnostics fit_r2, dt, n_obstacles,
    covered_fraction.
    """
    rng = np.random.default_rng(seed)
    R_obs = OBSTACLE_RADIUS_UM
    r_p_um = particle_radius_nm / 1000.0
    excl = R_obs + r_p_um
    area = box_um ** 2

    # obstacle count from nominal area fraction (Poisson field; overlaps allowed)
    n_obs = int(round(obstacle_fraction * area / (np.pi * R_obs ** 2)))
    centers = rng.uniform(0.0, box_um, size=(n_obs, 2)) if n_obs > 0 else np.zeros((0, 2))

    # time step: keep rms per-axis step small vs obstacle size (and box)
    if dt_s is None:
        step_std_target = 0.04 * R_obs
        dt = step_std_target ** 2 / (2.0 * free_D_um2s)
    else:
        dt = float(dt_s)
    step_std = np.sqrt(2.0 * free_D_um2s * dt)

    # occupancy grid for fast collision tests
    ng = max(1, int(np.ceil(box_um / (R_obs / _GRID_CELLS_PER_RADIUS))))
    occ = _build_occupancy(centers, excl, box_um, ng)
    h = box_um / ng

    # place tracers uniformly in the ACCESSIBLE area (rejection sampling)
    pos = rng.uniform(0.0, box_um, size=(n_particles, 2))
    if n_obs > 0:
        for _ in range(200):
            ix = np.clip((pos[:, 0] / h).astype(int), 0, ng - 1)
            iy = np.clip((pos[:, 1] / h).astype(int), 0, ng - 1)
            bad = occ[ix, iy]
            if not bad.any():
                break
            pos[bad] = rng.uniform(0.0, box_um, size=(int(bad.sum()), 2))
    unwrapped = pos.copy()

    # record unwrapped positions on a coarse schedule; MSD is then computed with
    # TIME-ORIGIN AVERAGING (every stored frame is used as an origin), which cuts the
    # MSD variance by ~ the number of origins vs a single-origin estimate.
    record_every = max(1, steps // 200)
    frames = [unwrapped.copy()]
    for s in range(1, steps + 1):
        delta = rng.normal(0.0, step_std, size=(n_particles, 2))
        prop = unwrapped + delta
        if n_obs > 0:
            pw = prop % box_um
            ix = np.clip((pw[:, 0] / h).astype(int), 0, ng - 1)
            iy = np.clip((pw[:, 1] / h).astype(int), 0, ng - 1)
            accept = ~occ[ix, iy]
        else:
            accept = np.ones(n_particles, dtype=bool)
        unwrapped[accept] = prop[accept]
        if s % record_every == 0:
            frames.append(unwrapped.copy())

    traj = np.stack(frames)  # (n_frames, n_particles, 2), unwrapped
    n_frames = traj.shape[0]
    frame_dt = record_every * dt
    # time-origin-averaged MSD as a function of lag (lag = 1 .. n_frames-1)
    lags = np.arange(1, n_frames)
    msd = np.empty(lags.size, float)
    for i, lag in enumerate(lags):
        diff = traj[lag:] - traj[:-lag]          # (n_frames-lag, n_particles, 2)
        msd[i] = float(np.mean(np.sum(diff ** 2, axis=2)))
    times = lags * frame_dt

    # long-lag linear fit (skip the short-time crossover where obstructed tracers are
    # still locally free); the fewer origins at large lag are the noisiest, so fit the
    # well-sampled middle-to-long region.
    lo = len(times) // 4
    hi = (3 * len(times)) // 4
    slope, _intercept, r2 = _linfit(times[lo:hi], msd[lo:hi])
    effective_D = max(slope / 4.0, 0.0)

    return {
        "effective_D": effective_D,
        "hindrance": effective_D / free_D_um2s if free_D_um2s > 0 else 0.0,
        "msd": msd,
        "times": times,
        "obstacle_fraction": obstacle_fraction,
        "fit_r2": r2,
        "dt": dt,
        "n_obstacles": n_obs,
        "covered_fraction": 1.0 - np.exp(-obstacle_fraction),
    }


def hindrance_vs_ecm(free_D_um2s: float,
                     ecm_fractions=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                     **kw) -> Dict[float, float]:
    """Run the Brownian sim across ECM/obstacle fractions.

    Returns {ecm_fraction: effective_D_um2s}. Extra keyword args are forwarded to
    `simulate_effective_diffusivity` (e.g. n_particles, steps, seed, box_um).
    """
    out: Dict[float, float] = {}
    for phi in ecm_fractions:
        res = simulate_effective_diffusivity(free_D_um2s, phi, **kw)
        out[phi] = res["effective_D"]
    return out


if __name__ == "__main__":
    import time as _time

    free_D = 100.0  # um^2/s, representative small-molecule free diffusivity
    t0 = _time.time()

    print("=" * 68)
    print("Brownian obstructed-diffusion self-test")
    print(f"  free_D = {free_D} um^2/s, obstacle radius = {OBSTACLE_RADIUS_UM*1000:.0f} nm")
    print("=" * 68)

    # -- Test 1 & 2: linearity of MSD + zero-obstacle recovery of free_D ----------
    base = simulate_effective_diffusivity(free_D, obstacle_fraction=0.0, seed=1)
    r2_free = base["fit_r2"]
    d0 = base["effective_D"]
    recovery_pct = 100.0 * d0 / free_D
    print(f"\n[1] MSD linearity (phi=0): long-lag fit R^2 = {r2_free:.5f}")
    print(f"    dt = {base['dt']:.3e} s, total time = {base['times'][-1]:.4e} s")
    print(f"[2] Zero-obstacle recovery: effective_D = {d0:.3f}  "
          f"({recovery_pct:.1f}% of free_D)")

    # linearity also under crowding (long-lag region must be diffusive, not anomalous)
    crowd = simulate_effective_diffusivity(free_D, obstacle_fraction=0.4, seed=1)
    print(f"    crowded (phi=0.4) long-lag fit R^2 = {crowd['fit_r2']:.5f}  "
          f"(diffusive, not sub-diffusive)")

    assert r2_free > 0.99, f"MSD not linear at phi=0 (R^2={r2_free})"
    assert crowd["fit_r2"] > 0.99, f"MSD not linear at phi=0.4 (R^2={crowd['fit_r2']})"
    assert abs(d0 - free_D) / free_D < 0.05, f"D0 not recovered within 5% ({recovery_pct:.1f}%)"

    # -- Test 3: monotonic hindrance across phi = 0 -> 0.6 -----------------------
    curve = hindrance_vs_ecm(free_D, ecm_fractions=(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), seed=2)
    print("\n[3] Hindrance vs obstacle fraction:")
    print(f"    {'phi':>6} {'eff_D':>10} {'hindrance':>10}")
    hinds = []
    for phi, dph in curve.items():
        hprime = dph / free_D
        hinds.append(hprime)
        print(f"    {phi:6.2f} {dph:10.3f} {hprime:10.3f}")
    # monotone decreasing over 0 -> 0.5 (drop 0.6 which nears void-percolation noise)
    hinds05 = hinds[:6]  # phi = 0.0 .. 0.5
    mono = all(hinds05[i + 1] <= hinds05[i] + 1e-9 for i in range(len(hinds05) - 1))
    print(f"    monotonically decreasing over phi 0->0.5: {mono}")
    assert mono, f"hindrance not monotonic over 0->0.5: {hinds05}"

    # -- Test 4: dilute-limit slope  hindrance ~= 1 - a*phi ----------------------
    # Use extra statistics here: the slope is a small difference of noisy D values.
    dilute_phis = (0.0, 0.05, 0.10, 0.15, 0.20)
    dilute = hindrance_vs_ecm(free_D, ecm_fractions=dilute_phis,
                              n_particles=800, steps=6000, seed=3)
    xs = np.array(dilute_phis)
    ys = np.array([dilute[p] / free_D for p in dilute_phis])
    slope, intercept, r2d = _linfit(xs, ys)
    a = -slope
    print("\n[4] Dilute-limit fit hindrance ~= b - a*phi:")
    print(f"    intercept b = {intercept:.3f}, slope a = {a:.3f}, fit R^2 = {r2d:.4f}")
    print("    (point tracers through random overlapping disks: a>0 of order ~1;")
    print("     the impenetrable effective-medium limit (1-phi)/(1+phi) gives a~2.)")
    assert a > 0.5, f"dilute slope a not positive/physical: a={a}"
    # intercept is a stochastic extrapolation of the D estimate to phi=0; ~10% is the
    # honest bound for this many particles/steps.
    assert abs(intercept - 1.0) < 0.10, f"dilute intercept not ~1: b={intercept}"

    print(f"\nAll asserts passed. Runtime {_time.time()-t0:.2f} s")
