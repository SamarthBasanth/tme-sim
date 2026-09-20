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

Time-resolved sources
---------------------
``S`` may be either **static** (a single ``(ny, nx)`` field applied unchanged at
every step — the historical behaviour and still the default) or **time-resolved**
via the optional ``source_profile`` argument, a per-step scalar multiplier so
that ``S(x, y, t_n) = source_profile[n] * source(x, y)``. This is what lets the
engine drive the solver with the real multi-dose plasma curve ``plasma(t)``
instead of its time-average, which is the only way concentration can genuinely
rise and fall and therefore the only way *when* a drug is given can change the
answer.

It costs essentially nothing: the operator ``A = I - dt*L + dt*K`` is
time-independent, so the LU factorization is still computed once and reused; a
time-varying source only changes the RHS by one extra scalar multiply per step.

Why mass is conserved exactly (clearance=binding=0, no source): ``L`` has zero
row/column sums, so ``1^T L = 0`` and hence ``1^T A = 1^T``. From ``A C^{n+1} =
C^n`` this gives ``1^T C^{n+1} = 1^T A^{-1} C^n = 1^T C^n``. Upwind advection uses
antisymmetric face fluxes that vanish on the boundary, so it too conserves the
tissue-total mass. Conservation therefore holds to the linear-solver tolerance.

Units: length um, time s, concentration uM, diffusivity um^2/s.
"""
from __future__ import annotations

from typing import Dict, Optional

import warnings

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
    source_profile: Optional[np.ndarray] = None,
    capture_frames: Optional[int] = None,
    exchange_k: Optional[np.ndarray] = None,
) -> TransportResult:
    """Solve convection-diffusion-reaction for drug concentration ``C(x, y, t)``.

    Governing equation (per :mod:`tmesim` architecture contract)::

        dC/dt = div(D_eff grad C) - div(v C) + S - (clearance + binding) C

    With ``exchange_k`` supplied, the reaction rate becomes
    ``(clearance + binding + exchange_k)`` and the equation is the genuine
    transvascular-exchange form ``S = PS*(Cp - C)`` rather than a one-way leak;
    see the comment at the ``exchange_k`` handling below for why those are the
    same thing and what it fixes.

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
        bbb_penetrance * plasma. **Time profile:** with ``source_profile=None``
        (the default) the plasma level is held constant over the window
        (``plasma(t) = 1``), i.e. ``S(x,y,t) = source(x,y)`` applied every step,
        and the caller is expected to have pre-scaled ``source`` by a single
        representative plasma level. Pass ``source_profile`` to drive a genuine
        time-varying infusion (see below). Source is zeroed outside the tissue
        mask.
    pde : PDEConfig
        ``dt``, ``total_time``, ``convection``, ``save_every``. Note this solver
        does NOT read ``pde.time_resolved`` — time resolution is requested here by
        passing ``source_profile``. The flag lives one level up, in
        :func:`tmesim.sim.engine.simulate_response`, which is what turns
        ``cfg.pk`` into a profile and passes it down.
    velocity : np.ndarray, optional, shape (2, ny, nx)
        Interstitial velocity ``(vx, vy)`` in um/s. Advection is included only
        when ``pde.convection`` is True and ``velocity`` is not None (otherwise
        it falls back to ``micro.fields.velocity`` if present).
    source_profile : np.ndarray, optional, shape (n_steps,) or longer
        **Time-resolved delivery.** A dimensionless per-step scale factor applied
        to ``source``, so the PDE actually integrated becomes

            S(x, y, t_n) = source_profile[n] * source(x, y)

        with ``t_n = n*dt``, ``n = 0 .. n_steps-1``. Typically
        ``source`` is the *unit* delivery field (``delivery_source(..., 1.0)``)
        and ``source_profile`` is ``plasma_profile(pk, t)`` in uM, so the product
        is the true instantaneous delivery rate. Only the first ``n_steps``
        entries are used; a longer array is accepted (and truncated) so callers
        may pass a profile sampled on a finer/longer grid. Must be finite and
        non-negative.

        The source enters the RHS, so it is evaluated at the OLD time level
        ``t_{n}`` for the step ``t_n -> t_{n+1}`` — the same explicit convention
        the upwind advection term already uses. Consequence: the total drug
        injected over the window is ``dt * sum(source_profile[:n_steps]) *
        sum(source)``, which is *exactly* the static-source total when
        ``source_profile`` has mean 1 on those same sample points. That identity
        is what makes "matched cumulative exposure" comparisons meaningful.

        **Backward compatibility is exact.** ``source_profile=None`` runs the
        original statement ``rhs += dt * src`` untouched, so the default path is
        bit-for-bit what it was. ``source_profile`` of all ones is *also*
        bit-identical (``dt * 1.0 == dt`` in IEEE-754); both are covered by
        tests.

        Note the operator ``A = I - dt*L + dt*K`` does not depend on the source,
        so the single LU factorization is still built once and reused for every
        step: a time-varying source costs one extra scalar multiply per step and
        nothing else.
    capture_frames : int, optional
        **Time-lapse.** Number of concentration snapshots to keep, evenly spaced
        over the window and always including ``t = 0`` (the zero initial
        condition) and ``t = total_time``. The requested count is clipped to
        ``n_steps + 1``, so a short solve returns every step rather than
        duplicating frames. ``None`` (the default) keeps nothing and is the
        historical path; capture is a pure ``copy`` at selected steps, so it
        changes no arithmetic and costs ``n_frames * ny * nx * 8`` bytes and
        nothing else. The snapshots land in ``TransportResult.frames`` /
        ``.frame_times``.

        The snapshots are the concentration ``C``, not the exposure ``AUC``:
        AUC only ever rises, so it cannot show drug leaving. Watching drug
        actually *clear* additionally needs a time-resolved plasma driver
        (``source_profile``), because a static source holds the delivery rate up
        for the whole window by construction.

    Returns
    -------
    TransportResult
        With ``auc`` (integral of C dt, uM*s), ``final_concentration``,
        ``times``, ``mean_concentration_history`` and ``penetration_index``.

    Raises
    ------
    ValueError
        On shape mismatches, non-positive ``dt``/``total_time``, a
        ``source_profile`` that is not 1-D or is shorter than ``n_steps``, or if
        an explicit advection step would violate the upwind CFL condition.
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
    # Optional transvascular BACK-FLUX. The historical source is one-way:
    # S = PS*Cp with no return path, so the steady state is set by
    # (influx)/(loss) and has no reason to sit below the plasma level -- it can
    # exceed it without bound (measured: results/breast_kp_validation.md, where
    # 5-FU predicts 7.8x plasma against 0.61x measured). A genuine exchange
    # term S = PS*(Cp - C) is algebraically just that same source MINUS a
    # first-order loss PS*C, so it needs no new solver machinery: the caller
    # passes the SAME PS field here that it used to build ``source``, and it is
    # added to the existing reaction rate. The steady state then becomes
    # C_ss = PS*Cp/(PS + clearance + binding) <= Cp, i.e. bounded by the plasma
    # level BY CONSTRUCTION rather than by luck.
    #
    # ``None`` (the default) short-circuits entirely and reproduces the
    # one-way behaviour bit-for-bit, so every stored reference result is
    # unaffected unless a caller opts in.
    if exchange_k is not None:
        ex = np.asarray(exchange_k, dtype=np.float64)
        if ex.shape != (ny, nx):
            raise ValueError(f"exchange_k shape {ex.shape} != grid {(ny, nx)}")
        if not np.all(np.isfinite(ex)):
            raise ValueError("exchange_k contains non-finite values")
        if np.any(ex < 0.0):
            raise ValueError(
                f"exchange_k must be non-negative (min {float(ex.min()):.4g}); "
                "a negative exchange rate would pump drug INTO the tissue "
                "against its own gradient."
            )
        k = k + ex
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
        #
        # TWO QUANTITIES, AND THEY ARE NOT THE SAME NUMBER.
        #
        # `max_speed` is a max over CELLS of the POINTWISE SUM. The condition
        # that actually guarantees positivity takes the two component maxima
        # INDEPENDENTLY and allows half as much:
        #
        #     dt * (max|vx| + max|vy|) / dx  <=  1/2
        #
        # They coincide only when the two maxima sit in the same cell. A
        # divergent cross -- four neighbours each flowing away, no single cell
        # carrying both components -- separates them by exactly 2x: every face
        # of the centre cell is outflow, the step removes 2*v*dt/dx of its
        # content, and the guard below sees only v. MEASURED 2026-09-20: at
        # CFL = 1, which this guard ACCEPTS, the undershoot is -100% of the
        # peak, and the clamp at the end of the step then CREATES mass equal to
        # everything injected. The measured undershoot follows -(2*CFL - 1) to
        # 4.0e-07.
        #
        # THE REFUSAL IS DELIBERATELY NOT TIGHTENED TO THE STRICT BOUND. On all
        # five real sections, at dx = 80/30/12 and dt = 0.5/10, the measured
        # raw undershoot before any clamp is 0.000e+00: real Darcy velocity
        # fields are smooth, so their component maxima are near-colocated and
        # the two quantities nearly agree. Refusing at 1/2 would reject this
        # product's own shipped default spec, which has never produced an
        # undershoot -- a guard that refuses the thing it is meant to protect
        # is not a stricter guard, it is a broken one.
        #
        # What IS reported is the margin, because the margin is gone: at the
        # shipped dx = 12 um and dt = 10 s the strict quantity is 0.6506,
        # already past 1/2, and refining the grid erodes it further. Production
        # is currently safe by the smoothness of the flow field and not by this
        # check, and that distinction belongs in the record rather than in
        # somebody's memory.
        max_speed = float(np.max(np.abs(vx) + np.abs(vy))) if velocity.size else 0.0
        cfl = max_speed * dt / dx
        if cfl > 1.0:
            raise ValueError(
                f"Advection CFL violated: dt*(|vx|+|vy|)/dx = {cfl:.3f} > 1. "
                f"Reduce pde.dt below {dx / max_speed:.4g} s or disable convection."
            )
        strict_speed = (float(np.max(np.abs(vx))) + float(np.max(np.abs(vy)))
                        if velocity.size else 0.0)
        strict_cfl = strict_speed * dt / dx
        if strict_cfl > 0.5:
            warnings.warn(
                "advection is past the STRICT positivity bound: "
                f"dt*(max|vx| + max|vy|)/dx = {strict_cfl:.4f} > 0.5. The "
                f"necessary CFL check above passed at {cfl:.4f} <= 1 and the "
                "solve proceeds, because on real Darcy fields the two maxima "
                "are near-colocated and no undershoot has been measured. But "
                "the guarantee is gone: positivity now rests on the smoothness "
                "of this velocity field, not on the check. If the "
                "post-step max(C, 0) clamp ever fires, this is why.",
                RuntimeWarning,
                stacklevel=2,
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

    # --- optional time-resolved source scaling (None => original static path) ---
    if source_profile is None:
        src_scale = None
    else:
        src_scale = np.asarray(source_profile, dtype=np.float64)
        if src_scale.ndim != 1:
            raise ValueError(
                f"source_profile must be 1-D, got shape {src_scale.shape}"
            )
        if src_scale.size < n_steps:
            raise ValueError(
                f"source_profile has {src_scale.size} samples but the solve takes "
                f"{n_steps} steps; supply at least one scale factor per step."
            )
        used = src_scale[:n_steps]
        if not np.all(np.isfinite(used)):
            raise ValueError("source_profile contains non-finite values")
        # A negative plasma level is meaningless, and because the step clamps C at
        # zero it would silently break mass balance rather than fail loudly.
        if float(used.min()) < 0.0:
            raise ValueError(
                f"source_profile must be non-negative (min {float(used.min()):.4g}): "
                f"a negative plasma level would silently break mass conservation, "
                f"because the solver clamps concentration at zero each step."
            )

    c = domain.zeros()  # initial condition: no drug in tissue
    auc = domain.zeros()

    times = [0.0]
    mean_hist = [float(c[tissue].mean()) if tissue.any() else 0.0]

    # --- optional time-lapse capture (None => nothing allocated, nothing kept) --
    frames = frame_times = None
    frame_interval_s = None
    capture_at: Dict[int, int] = {}
    if capture_frames is not None:
        n_frames = int(capture_frames)
        if n_frames < 2:
            raise ValueError(
                f"capture_frames must be >= 2 (a time-lapse needs at least a "
                f"first and a last frame); got {capture_frames!r}"
            )
        n_frames = min(n_frames, n_steps + 1)
        steps_wanted = np.unique(
            np.rint(np.linspace(0.0, float(n_steps), n_frames)).astype(int))
        capture_at = {int(s): i for i, s in enumerate(steps_wanted)}
        n_frames = steps_wanted.size
        frames = np.zeros((n_frames, ny, nx), dtype=np.float64)
        frame_times = steps_wanted.astype(np.float64) * dt
        frame_interval_s = (
            float(pde.total_time) / (n_frames - 1) if n_frames > 1 else None)
        frames[0] = c  # t = 0, the zero initial condition

    src_flat = src.ravel()
    for step in range(1, n_steps + 1):
        rhs = c.copy()
        if src_scale is None:
            rhs += dt * src
        else:
            # Source at the OLD time level, matching the explicit advection term.
            rhs += (dt * src_scale[step - 1]) * src
        if use_convection:
            rhs -= dt * _advective_divergence(c, vx, vy, tissue, dx)
        c_new = solve(rhs.ravel()).reshape(ny, nx)
        # Guard tiny negative undershoot from upwind/round-off; keep non-negative.
        np.maximum(c_new, 0.0, out=c_new)
        c_new[~tissue] = 0.0

        # Trapezoidal accumulation of the exposure integral (AUC).
        auc += 0.5 * (c + c_new) * dt
        c = c_new

        if capture_at:
            slot = capture_at.get(step)
            if slot is not None:
                frames[slot] = c

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
        frames=frames,
        frame_times=frame_times,
        frame_interval_s=frame_interval_s,
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


def _test_time_resolved_source():
    """The ``source_profile`` path: inert by default, exact when it acts.

    (a) an all-ones profile is BIT-identical to ``source_profile=None``;
    (b) mass balance still holds exactly with a strongly time-varying source;
    (c) two profiles with the same window MEAN inject the same total drug, but
        produce a different exposure field — i.e. the schedule reaches the answer;
    (d) an early (bolus-like) profile realises MORE exposure inside the window
        than a late one at identical total dose, because less drug is left in
        flight at T. This is the direction the terminal-residual argument
        predicts (see the module docstring of ``tmesim.sim.engine``).
    """
    ny = nx = 40
    dx = 15.0
    micro = _make_micro(ny, nx, dx, ecm_value=0.3)
    inert = DrugParams(name="inert", diffusivity=300.0, clearance=0.0, binding=0.0)
    reactive = DrugParams(name="rx", diffusivity=300.0, clearance=2e-3, binding=0.0)
    pde = PDEConfig(dt=10.0, total_time=3600.0, convection=False, save_every=50)
    n_steps = int(round(pde.total_time / pde.dt))
    source = micro.domain.zeros()
    source[ny // 2, nx // 2] = 2.0

    # (a) inert default path
    base = simulate_drug_transport(micro, reactive, source, pde)
    ones = simulate_drug_transport(micro, reactive, source, pde,
                                   source_profile=np.ones(n_steps))
    identical = (np.array_equal(base.auc, ones.auc)
                 and np.array_equal(base.final_concentration,
                                    ones.final_concentration))
    print(f"[test4a] ones-profile bit-identical to None: {identical}")
    assert identical, "an all-ones source_profile must reproduce the static path exactly"

    # A ramp and its mirror: same window mean, opposite ordering in time.
    early = np.linspace(2.0, 0.0, n_steps)
    late = early[::-1].copy()
    assert abs(early.mean() - late.mean()) < 1e-12

    # (b) + (c) mass balance and non-degeneracy, on the conservative tracer
    r_early = simulate_drug_transport(micro, inert, source, pde,
                                      source_profile=early)
    r_late = simulate_drug_transport(micro, inert, source, pde,
                                     source_profile=late)
    injected = pde.dt * float(early.sum()) * float(source.sum())
    m_e = float(r_early.final_concentration.sum())
    m_l = float(r_late.final_concentration.sum())
    err = max(abs(m_e - injected), abs(m_l - injected)) / injected
    print(f"[test4b] time-varying source mass balance: injected={injected:.6f} "
          f"early={m_e:.6f} late={m_l:.6f} rel_err={err:.2e}")
    assert err < 1e-9, f"mass not conserved with a time-varying source: {err:.2e}"
    assert not np.allclose(r_early.auc, r_late.auc), (
        "same total dose, different timing -> the exposure field must differ"
    )

    # (d) direction, on the reactive drug (a loss term is what makes exposure
    #     depend on WHEN, via the residual left un-eliminated at T).
    a_early = simulate_drug_transport(micro, reactive, source, pde,
                                      source_profile=early).auc.sum()
    a_late = simulate_drug_transport(micro, reactive, source, pde,
                                     source_profile=late).auc.sum()
    print(f"[test4d] exposure: early-weighted={a_early:.6g} "
          f"late-weighted={a_late:.6g}  ratio={a_early/a_late:.6f}")
    assert a_early > a_late, (
        "front-loaded delivery must realise more exposure inside the window than "
        "back-loaded delivery at identical total dose"
    )
    return err, float(a_early / a_late)


if __name__ == "__main__":
    print("Running diffusion.py self-tests...\n")
    e1 = _test_mass_conservation()
    print()
    p_low, p_high = _test_penetration_monotonicity()
    print()
    rt = _test_runtime_120()
    print()
    tr_err, tr_ratio = _test_time_resolved_source()
    print("\nAll self-tests PASSED.")
    print(f"  mass-conservation rel err : {e1:.2e}")
    print(f"  penetration low/high ECM  : {p_low:.4f} / {p_high:.4f} "
          f"(monotone decreasing: {p_high < p_low})")
    print(f"  120x120 runtime           : {rt:.2f}s")
    print(f"  time-resolved mass err    : {tr_err:.2e}  "
          f"early/late exposure ratio: {tr_ratio:.4f}")
