"""Core numerics: convection-diffusion-reaction drug transport on the TME grid.

This is the numerical heart of tmesim. It solves, for a drug concentration
field ``C(x, y, t)`` (uM) on the ``(ny, nx)`` microenvironment grid,

    dC/dt = div(D_eff grad C)  -  div(v C)  +  S  -  (clearance + binding) C

subject to **no-flux (Neumann) boundaries** at the domain edges *and* at every
tissue/non-tissue interface (``micro.domain.tissue_mask == False``).

Numerical scheme (and why it is stable)
---------------------------------------
We use **operator splitting with an implicit (backward-Euler) diffusion+reaction
step and an explicit first-order upwind advection step** (an IMEX scheme):

    (I  -  dt * L  +  dt * K) C^{n+1}  =  C^n  +  dt * (S  -  Adv(C^n))

* ``L`` is the discrete variable-coefficient diffusion operator, built in
  **finite-volume form** so that the flux leaving cell P through a face equals
  minus the flux entering the neighbour through that same face. ``L`` is
  symmetric, its rows sum to zero, and face coefficients are set to **zero** at
  domain edges and tissue-mask boundaries -> this *is* the no-flux condition and
  makes the scheme **exactly mass-conserving** (see below).
* Because diffusion and the linear decay ``K = clearance + binding`` are treated
  **implicitly**, the diffusion update is **unconditionally stable**: the system
  matrix ``A = I - dt*L + dt*K`` is symmetric, diagonally dominant and positive
  definite for any ``dt > 0`` (there is *no* ``dt <= dx^2/(4D)`` restriction).
  ``A`` is time-independent, so we LU-factorize it **once** and reuse the factors
  for every step -> a few thousand steps on a 120x120 grid run in seconds.
* Advection uses a stable **upwind** discretisation and is the only explicit
  piece, so it carries the mild CFL condition ``dt*(|vx|+|vy|)/dx <= 1``. Because
  interstitial velocities are tiny (~<1 um/s), this is satisfied with room to
  spare; we nonetheless *check it and raise* a clear error if it is violated.

Why mass is conserved exactly (clearance=binding=0, no source): ``L`` has zero
row/column sums, so ``1^T L = 0`` and hence ``1^T A = 1^T``. From ``A C^{n+1} =
C^n`` this gives ``1^T C^{n+1} = 1^T A^{-1} C^n = 1^T C^n``. Upwind advection uses
antisymmetric face fluxes that vanish on the boundary, so it too conserves the
tissue-total mass. Conservation therefore holds to the linear-solver tolerance.

Units: length um, time s, concentration uM, diffusivity um^2/s.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import factorized

from ..core import (
    CellPopulation,
    DrugParams,
    Microenvironment,
    PDEConfig,
    TransportResult,
)

# ECM tortuosity coupling: D_eff = D * (1 - ALPHA_TORTUOSITY * ecm_density).
ALPHA_TORTUOSITY = 0.6


def _effective_diffusivity(drug: DrugParams, ecm_density: np.ndarray) -> np.ndarray:
    """Spatially varying effective diffusivity reduced by ECM tortuosity.

    ``D_eff = drug.diffusivity * (1 - ALPHA_TORTUOSITY * ecm_density)``, floored
    at a small positive value so the diffusion operator stays well-posed even in
    dense-ECM regions where the raw expression could go to (or below) zero.
    """
    ecm = np.clip(np.asarray(ecm_density, dtype=np.float64), 0.0, 1.0)
    d = drug.diffusivity * (1.0 - ALPHA_TORTUOSITY * ecm)
    floor = max(1e-6, 1e-3 * drug.diffusivity)
    return np.maximum(d, floor)


def _build_diffusion_operator(
    d_eff: np.ndarray, tissue: np.ndarray, dx: float
) -> csc_matrix:
    """Assemble the symmetric finite-volume diffusion operator ``L`` (N x N).

    ``(L C)[P] = (1/dx^2) * sum_faces D_face * (C_neighbour - C_P)`` with the
    harmonic mean of the two adjacent cell diffusivities at each face. Faces to
    the domain exterior or to non-tissue cells get coefficient 0 (no-flux). Rows
    of non-tissue cells are all-zero (they are decoupled and held at 0).
    """
    ny, nx = d_eff.shape
    d2 = dx * dx
    idx = np.arange(ny * nx).reshape(ny, nx)
    diag = np.zeros((ny, nx), dtype=np.float64)

    rows_list = []
    cols_list = []
    data_list = []

    def _harmonic(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return 2.0 * a * b / (a + b)

    # --- x-direction faces: between (i, j) and (i, j+1) ---
    if nx > 1:
        df = _harmonic(d_eff[:, :-1], d_eff[:, 1:]) / d2  # (ny, nx-1)
        active = tissue[:, :-1] & tissue[:, 1:]
        df = np.where(active, df, 0.0)
        p = idx[:, :-1].ravel()
        e = idx[:, 1:].ravel()
        w = df.ravel()
        rows_list += [p, e]
        cols_list += [e, p]
        data_list += [w, w]
        np.subtract.at(diag.ravel(), p, w)
        np.subtract.at(diag.ravel(), e, w)

    # --- y-direction faces: between (i, j) and (i+1, j) ---
    if ny > 1:
        df = _harmonic(d_eff[:-1, :], d_eff[1:, :]) / d2  # (ny-1, nx)
        active = tissue[:-1, :] & tissue[1:, :]
        df = np.where(active, df, 0.0)
        p = idx[:-1, :].ravel()
        s = idx[1:, :].ravel()
        w = df.ravel()
        rows_list += [p, s]
        cols_list += [s, p]
        data_list += [w, w]
        np.subtract.at(diag.ravel(), p, w)
        np.subtract.at(diag.ravel(), s, w)

    rows_list.append(idx.ravel())
    cols_list.append(idx.ravel())
    data_list.append(diag.ravel())

    rows = np.concatenate(rows_list)
    cols = np.concatenate(cols_list)
    data = np.concatenate(data_list)
    n = ny * nx
    return csc_matrix((data, (rows, cols)), shape=(n, n))


def _advective_divergence(
    c: np.ndarray, vx: np.ndarray, vy: np.ndarray, tissue: np.ndarray, dx: float
) -> np.ndarray:
    """First-order upwind approximation of ``div(v C)`` on the grid (ny, nx).

    Face velocities are the average of the two adjacent cell velocities; the
    upwind (donor-cell) concentration is used for the flux. Fluxes across
    domain-edge and tissue/non-tissue faces are zeroed (no-flux), so the total
    over tissue cells is conserved.
    """
    ny, nx = c.shape
    div = np.zeros((ny, nx), dtype=np.float64)

    if nx > 1:
        vf = 0.5 * (vx[:, :-1] + vx[:, 1:])              # (ny, nx-1)
        c_up = np.where(vf > 0.0, c[:, :-1], c[:, 1:])
        flux = vf * c_up
        active = tissue[:, :-1] & tissue[:, 1:]
        flux = np.where(active, flux, 0.0) / dx
        div[:, :-1] += flux   # east-face flux out of cell j
        div[:, 1:] -= flux    # same face is the west face of cell j+1

    if ny > 1:
        vf = 0.5 * (vy[:-1, :] + vy[1:, :])              # (ny-1, nx)
        c_up = np.where(vf > 0.0, c[:-1, :], c[1:, :])
        flux = vf * c_up
        active = tissue[:-1, :] & tissue[1:, :]
        flux = np.where(active, flux, 0.0) / dx
        div[:-1, :] += flux
        div[1:, :] -= flux

    return div


def simulate_drug_transport(
    micro: Microenvironment,
    drug: DrugParams,
    source: np.ndarray,
    pde: PDEConfig,
    velocity: Optional[np.ndarray] = None,
) -> TransportResult:
    """Solve convection-diffusion-reaction for drug concentration ``C(x, y, t)``.

    Governing equation (per :mod:`tmesim` architecture contract)::

        dC/dt = div(D_eff grad C) - div(v C) + S - (clearance + binding) C

    Parameters
    ----------
    micro : Microenvironment
        Provides the ``domain`` (grid + ``tissue_mask``) and ``fields``
        (``ecm_density`` sets the tortuosity that reduces diffusivity).
    drug : DrugParams
        ``diffusivity`` (um^2/s), ``clearance`` and ``binding`` (1/s) first-order
        loss rates. The effective, spatially varying diffusivity is
        ``D_eff = diffusivity * (1 - 0.6 * ecm_density)`` floored positive.
    source : np.ndarray, shape (ny, nx)
        Base spatial delivery-rate field (uM/s), e.g. from vasculature *
        bbb_penetrance * plasma. **Time profile:** we treat the plasma level as
        constant over the window (``plasma(t) = 1``), i.e. ``S(x,y,t) =
        source(x,y)`` applied every step. A caller wanting a time-varying
        infusion should pre-scale ``source`` by the (constant) plasma level for
        this window; the engine drives the full multi-dose profile by calling
        this solver per dose window. Source is zeroed outside the tissue mask.
    pde : PDEConfig
        ``dt``, ``total_time``, ``convection``, ``save_every``.
    velocity : np.ndarray, optional, shape (2, ny, nx)
        Interstitial velocity ``(vx, vy)`` in um/s. Advection is included only
        when ``pde.convection`` is True and ``velocity`` is not None (otherwise
        it falls back to ``micro.fields.velocity`` if present).

    Returns
    -------
    TransportResult
        With ``auc`` (integral of C dt, uM*s), ``final_concentration``,
        ``times``, ``mean_concentration_history`` and ``penetration_index``.

    Raises
    ------
    ValueError
        On shape mismatches, non-positive ``dt``/``total_time``, or if an
        explicit advection step would violate the upwind CFL condition.
    """
    domain = micro.domain
    ny, nx = domain.shape
    dx = float(domain.dx)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)

    source = np.asarray(source, dtype=np.float64)
    if source.shape != (ny, nx):
        raise ValueError(f"source shape {source.shape} != grid {(ny, nx)}")
    if pde.dt <= 0.0:
        raise ValueError(f"pde.dt must be > 0, got {pde.dt}")
    if pde.total_time <= 0.0:
        raise ValueError(f"pde.total_time must be > 0, got {pde.total_time}")

    dt = float(pde.dt)

    # --- effective diffusivity + reaction rate (both zeroed outside tissue) ---
    d_eff = _effective_diffusivity(drug, micro.fields.ecm_density)
    k = np.full((ny, nx), float(drug.clearance) + float(drug.binding))
    k[~tissue] = 0.0
    src = np.where(tissue, source, 0.0)

    # --- advection setup ---
    use_convection = bool(pde.convection)
    if velocity is None:
        velocity = getattr(micro.fields, "velocity", None)
    if velocity is None:
        use_convection = False
    if use_convection:
        velocity = np.asarray(velocity, dtype=np.float64)
        if velocity.shape != (2, ny, nx):
            raise ValueError(
                f"velocity shape {velocity.shape} != {(2, ny, nx)}"
            )
        vx, vy = velocity[0], velocity[1]
        # Explicit upwind CFL check (diffusion is implicit -> no diffusion CFL).
        max_speed = float(np.max(np.abs(vx) + np.abs(vy))) if velocity.size else 0.0
        cfl = max_speed * dt / dx
        if cfl > 1.0:
            raise ValueError(
                f"Advection CFL violated: dt*(|vx|+|vy|)/dx = {cfl:.3f} > 1. "
                f"Reduce pde.dt below {dx / max_speed:.4g} s or disable convection."
            )
    else:
        vx = vy = None

    # --- assemble and factorize the implicit operator A = I - dt*L + dt*K ---
    n = ny * nx
    lap = _build_diffusion_operator(d_eff, tissue, dx)
    a = csc_matrix(
        (np.ones(n), (np.arange(n), np.arange(n))), shape=(n, n)
    )  # identity
    a = a - dt * lap
    a = a + csc_matrix(
        (dt * k.ravel(), (np.arange(n), np.arange(n))), shape=(n, n)
    )
    # Non-tissue rows: L rows already zero and k=0, so A row = identity -> held 0.
    a = a.tocsc()
    solve = factorized(a)

    # --- time stepping ---
    n_steps = int(round(pde.total_time / dt))
    save_every = max(1, int(pde.save_every))

    c = domain.zeros()  # initial condition: no drug in tissue
    auc = domain.zeros()

    times = [0.0]
    mean_hist = [float(c[tissue].mean()) if tissue.any() else 0.0]

    src_flat = src.ravel()
    for step in range(1, n_steps + 1):
        rhs = c.copy()
        rhs += dt * src
        if use_convection:
            rhs -= dt * _advective_divergence(c, vx, vy, tissue, dx)
        c_new = solve(rhs.ravel()).reshape(ny, nx)
        # Guard tiny negative undershoot from upwind/round-off; keep non-negative.
        np.maximum(c_new, 0.0, out=c_new)
        c_new[~tissue] = 0.0

        # Trapezoidal accumulation of the exposure integral (AUC).
        auc += 0.5 * (c + c_new) * dt
        c = c_new

        if step % save_every == 0 or step == n_steps:
            times.append(step * dt)
            mean_hist.append(float(c[tissue].mean()) if tissue.any() else 0.0)

    final_concentration = c

    # --- penetration index in [0,1]: mean AUC relative to a robust peak ---
    # A well-mixed (perfectly penetrating) field has mean == peak -> index 1.
    # A field concentrated near vessels has mean << peak -> index near 0. Higher
    # ECM tortuosity lowers D_eff, steepens gradients, and lowers this index.
    if tissue.any():
        auc_t = auc[tissue]
        peak = float(np.percentile(auc_t, 99.0))
        if peak > 0.0:
            penetration_index = float(np.clip(auc_t.mean() / peak, 0.0, 1.0))
        else:
            penetration_index = 0.0
    else:
        penetration_index = 0.0

    return TransportResult(
        auc=auc,
        final_concentration=final_concentration,
        times=np.asarray(times, dtype=np.float64),
        mean_concentration_history=np.asarray(mean_hist, dtype=np.float64),
        penetration_index=penetration_index,
    )


# ---------------------------------------------------------------------------
# Self-tests
# ---------------------------------------------------------------------------
def _make_micro(ny, nx, dx, ecm_value, tissue_mask=None):
    """Minimal hand-built Microenvironment for self-tests (no data generator)."""
    from ..core import Domain, TissueFields

    domain = Domain(nx=nx, ny=ny, dx=dx, tissue_mask=tissue_mask)
    ecm = np.full((ny, nx), float(ecm_value))
    vasc = np.zeros((ny, nx))
    fields = TissueFields(ecm_density=ecm, vascular_density=vasc)
    cells = CellPopulation(
        x=np.array([nx * dx / 2]), y=np.array([ny * dx / 2]),
        cell_type=np.array([0]), state=np.array([0]),
        sensitivity=np.array([0.5]),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells)


def _test_mass_conservation():
    ny = nx = 60
    dx = 15.0
    micro = _make_micro(ny, nx, dx, ecm_value=0.3)
    drug = DrugParams(name="inert", diffusivity=300.0, clearance=0.0, binding=0.0)

    # Positive blob, no source, no reaction, no convection.
    c0 = micro.domain.zeros()
    yy, xx = np.mgrid[0:ny, 0:nx]
    c0 += 10.0 * np.exp(-(((xx - nx / 2) ** 2 + (yy - ny / 2) ** 2) / (2 * 4.0 ** 2)))
    mass0 = float(c0.sum())

    # Inject the blob as the initial condition by running with a tailored solver:
    # simulate_drug_transport starts from zero, so we test the operator directly
    # by seeding C via a one-off backward-Euler stepping loop mirroring the code.
    pde = PDEConfig(dt=10.0, total_time=3600.0, convection=False, save_every=50)

    # Re-run the internal machinery with a non-zero initial condition.
    domain = micro.domain
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    d_eff = _effective_diffusivity(drug, micro.fields.ecm_density)
    n = ny * nx
    lap = _build_diffusion_operator(d_eff, tissue, dx)
    from scipy.sparse import identity as _identity
    a = (_identity(n, format="csc") - pde.dt * lap).tocsc()
    solve = factorized(a)
    c = c0.copy()
    n_steps = int(round(pde.total_time / pde.dt))
    min_val = c.min()
    for _ in range(n_steps):
        c = solve(c.ravel()).reshape(ny, nx)
        min_val = min(min_val, float(c.min()))
    mass1 = float(c.sum())

    rel_err = abs(mass1 - mass0) / mass0
    spread = float(c.max()) < float(c0.max())  # blob should have spread/flattened
    print(f"[test1a] blob diffusion: mass0={mass0:.6f} mass1={mass1:.6f} "
          f"rel_err={rel_err:.2e} min_seen={min_val:.2e} spread={spread}")
    assert rel_err < 1e-9, f"mass not conserved: rel_err={rel_err:.2e}"
    assert min_val >= -1e-9, f"concentration went negative: {min_val:.2e}"
    assert spread, "blob did not diffuse/spread"

    # (b) End-to-end through the PUBLIC solver: with clearance=binding=0 and a
    # localized source, no-flux boundaries mean nothing escapes, so the final
    # total mass must equal exactly the injected mass = n_steps*dt*sum(source).
    source = micro.domain.zeros()
    source[ny // 2, nx // 2] = 2.0  # uM/s at one vessel
    res = simulate_drug_transport(micro, drug, source, pde)
    injected = int(round(pde.total_time / pde.dt)) * pde.dt * float(source.sum())
    final_mass = float(res.final_concentration.sum())
    api_err = abs(final_mass - injected) / injected
    neg = float(res.final_concentration.min())
    print(f"[test1b] public solver: injected={injected:.4f} final={final_mass:.4f} "
          f"rel_err={api_err:.2e} min_C={neg:.2e}")
    assert api_err < 1e-9, f"public-solver mass not conserved: {api_err:.2e}"
    assert neg >= -1e-9, f"public solver produced negative C: {neg:.2e}"
    return rel_err


def _test_penetration_monotonicity():
    ny = nx = 80
    dx = 15.0
    drug = DrugParams(name="tmz", diffusivity=300.0, clearance=3e-4, binding=0.0)
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 2, convection=False, save_every=100)

    # A single central vessel as the drug source.
    source = np.zeros((ny, nx))
    source[ny // 2, nx // 2] = 1.0

    results = {}
    for ecm in (0.1, 0.9):
        micro = _make_micro(ny, nx, dx, ecm_value=ecm)
        res = simulate_drug_transport(micro, drug, source, pde)
        results[ecm] = res
        near = res.final_concentration[ny // 2, nx // 2]
        print(f"[test2] ecm={ecm:.1f} penetration={res.penetration_index:.4f} "
              f"near_vessel_C={near:.4f} mean_auc={res.auc.mean():.4g}")

    pen_low = results[0.1].penetration_index
    pen_high = results[0.9].penetration_index
    # Source raises concentration near the vessel.
    assert results[0.1].final_concentration[ny // 2, nx // 2] > 0.0
    # Higher ECM tortuosity -> lower penetration.
    assert pen_high < pen_low, (
        f"penetration not monotone: high-ECM {pen_high:.4f} !< low-ECM {pen_low:.4f}"
    )
    return pen_low, pen_high


def _test_runtime_120():
    import time
    ny = nx = 120
    dx = 15.0
    micro = _make_micro(ny, nx, dx, ecm_value=0.4)
    drug = DRUG_LIBRARY_get()
    source = np.zeros((ny, nx))
    source[::20, ::20] = 1.0  # scattered vessels
    # velocity field to exercise the convection path.
    vel = np.zeros((2, ny, nx))
    vel[0] = 0.05  # 0.05 um/s outward-ish uniform drift
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 3, convection=True, save_every=200)
    t0 = time.time()
    res = simulate_drug_transport(micro, drug, source, pde, velocity=vel)
    elapsed = time.time() - t0
    n_steps = int(round(pde.total_time / pde.dt))
    print(f"[test3] 120x120, {n_steps} steps (convection ON), "
          f"elapsed={elapsed:.2f}s, penetration={res.penetration_index:.4f}")
    assert elapsed < 10.0, f"too slow: {elapsed:.2f}s"
    assert res.auc.min() >= -1e-9
    return elapsed


def DRUG_LIBRARY_get():
    from ..core import DRUG_LIBRARY
    return DRUG_LIBRARY["temozolomide-like"]


if __name__ == "__main__":
    print("Running diffusion.py self-tests...\n")
    e1 = _test_mass_conservation()
    print()
    p_low, p_high = _test_penetration_monotonicity()
    print()
    rt = _test_runtime_120()
    print("\nAll self-tests PASSED.")
    print(f"  mass-conservation rel err : {e1:.2e}")
    print(f"  penetration low/high ECM  : {p_low:.4f} / {p_high:.4f} "
          f"(monotone decreasing: {p_high < p_low})")
    print(f"  120x120 runtime           : {rt:.2f}s")
