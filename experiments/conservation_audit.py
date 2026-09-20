"""CONSERVATION AUDIT — the physical invariants the transport solve must satisfy,
measured as PROPERTIES over classes of input rather than as fixtures over one run.

WHAT THIS IS
------------
``tests/test_physics.py`` already pins mass conservation and non-negativity — on
one 40x40 grid, with one drug, at one ``dt``. That is a fixture: it catches a
regression at that point and says nothing about the input class. This driver asks
the property form of the same questions:

    for ANY drug, ANY grid, ANY dt, ANY source profile, ANY Peclet number in the
    ladder, ANY of the 12 compiled GBM compounds on ANY of the 5 available real
    sections -- does the invariant hold, and if not, by how much and on what
    minimal input?

Nothing in ``tmesim/`` or ``app/`` is modified, tuned or fitted by this file. It
only drives the shipped code and reports. ``calibrate_alpha()`` is asserted at
22.223477711172563 at the start and the end of every run.

THE SOLVER MIRROR, AND WHY IT IS NEEDED
---------------------------------------
Two of the five invariants cannot be measured from ``TransportResult`` alone:

* the **discrete mass budget** needs the per-step loss ``dt * sum(k * c^{n+1})``,
  and ``TransportResult`` carries only the trapezoidal ``auc`` (which, as sec. 1
  shows, is the WRONG quadrature for the budget and misses by O(dt));
* the **non-negativity** question is about the field BEFORE
  ``np.maximum(c_new, 0.0, out=c_new)`` clamps it. The shipped result can never
  be negative; the scheme's overshoot is what "does this scheme overshoot?" means.

:func:`march_mirror` re-runs the loop of
:func:`tmesim.physics.diffusion.simulate_drug_transport` using **that module's own
private operators** (``_effective_diffusivity``, ``_build_diffusion_operator``,
``_advective_divergence``) so no discretisation is re-implemented, with a
``clamp`` switch and per-step accounting added. :func:`assert_mirror_is_production`
checks the mirror reproduces the public solver **bit-for-bit** (``np.array_equal``
on both ``final_concentration`` and ``auc``) on the static-source, time-resolved,
convection-off and convection-on paths, and is called at the head of every stage
that uses the mirror. If the mirror ever stops being the solver, every number it
produced stops being about the solver, and the run aborts.

THE PRE-REGISTRATION
--------------------
Every expectation and threshold below was fixed in :data:`PREREG` before the
definitive run and is asserted at the start of the run. Where an expectation is
a DERIVATION rather than a guess, the derivation is given, because that is what
makes a failed threshold informative:

A. MASS BALANCE. Summing the shipped step
   ``(I - dt*L + dt*K) c^{n+1} = c^n + dt*P_n*S - dt*Adv(c^n)``
   over every cell: ``L`` has zero row sums by construction (finite-volume, with
   face coefficients zeroed at the domain edge and at every tissue/non-tissue
   face), and ``Adv`` is assembled from antisymmetric face fluxes that are zeroed
   on exactly those same faces, so ``1^T L = 0`` and ``1^T Adv(c) = 0``. What
   survives is

       sum(c^{n+1}) + dt * sum(k * c^{n+1}) = sum(c^n) + dt * P_n * sum(S)

   and, cumulatively from ``c^0 = 0``,

       M_final = INJECTED - REMOVED
       INJECTED = dt * sum_n P_n * sum(S)
       REMOVED  = dt * sum_n sum(k * c^{n+1})          <- RIGHT-ENDPOINT, not trapezoid

   REMOVED splits into "cleared" and "bound" in the exact ratio
   ``clearance : binding`` (and, with ``transvascular_exchange`` on, a third
   ``exchange`` share), because ``k`` is the sum of those rates on the same cells.
   The loss quadrature is fixed by the scheme: the sink is evaluated IMPLICITLY,
   at the new time level, so the budget's loss term is a right-endpoint sum. The
   ``auc`` the solver returns is a TRAPEZOIDAL integral, so ``k * sum(auc)`` is a
   different quadrature of the same integral and must miss by O(dt). That is
   prediction A2, and it is the reason this budget is derived from the code
   rather than from the continuum equation.

   Tolerance 1e-10 relative. Reason: the operator is LU-factorised ONCE by
   SuperLU and applied by triangular solves, so each step's residual is
   ~eps*cond(A); ``A = I - dt*L + dt*K`` is a symmetric diagonally-dominant
   M-matrix whose condition number is ~1 + 8*dt*D_eff/dx^2 (a few hundred at the
   ladder's worst point), and N <= ~4000 steps accumulate at worst ~N*eps*cond
   ~ 1e-10. Anything at 1e-13 is round-off; anything at 1e-6 is a real leak.

B. NON-NEGATIVITY. ``A = I - dt*L + dt*K`` has strictly positive diagonal and
   non-positive off-diagonal entries (``-dt`` times ``L``'s positive face
   coefficients) and is strictly diagonally dominant, so it is an M-matrix and
   ``A^{-1} >= 0`` elementwise. Therefore ``c^{n+1} = A^{-1} rhs`` can only be
   negative if ``rhs`` is, and ``rhs = c^n + dt*P_n*S - dt*Adv(c^n)`` with the
   first two terms non-negative. **The implicit half cannot produce an
   undershoot at any Peclet number or any dt; only the explicit upwind advection
   can.** Prediction B1: the Peclet ladder with a smooth flow field shows
   undershoot exactly 0, not "small".

   Prediction B2 is the sharp one, and it was derived before it was run. The
   outflow fraction removed from cell P in one explicit step is

       (dt/dx) * [ max(vf_e,0) + max(-vf_w,0) + max(vf_n,0) + max(-vf_s,0) ]

   with each face velocity the mean of its two cells. The x pair is bounded by
   ``(vx_E - vx_W)/2 <= max|vx|`` and the y pair by ``max|vy|``, so positivity
   needs ``dt*(max|vx| + max|vy|)/dx <= 1`` with the two maxima taken
   INDEPENDENTLY over the grid. The shipped guard is
   ``dt*max(|vx| + |vy|)/dx <= 1``, a max over cells of the POINTWISE SUM, which
   is up to **2x smaller** whenever the two maxima are not co-located. So the
   guard admits ``CFL`` up to 1 where positivity needs 1/2, and on the minimal
   input that realises the factor exactly the undershoot must be

       min(c_raw) / max(c) = -(2*CFL - 1)      for CFL > 1/2

   i.e. **-100% of the peak at the guard's own accepted limit**. A scheme that
   produces -1e-9 is fine; this is the other kind. Threshold: any measured
   relative undershoot below -1e-3 is reported VIOLATED.

   B3 asks whether it is REACHED. A measured violation on a constructed field is
   a fact about the guard; whether production crosses the derived 1/2 bound is a
   separate measurement and is made on all five real sections at the shipped
   ``dx`` and ``dt``.

C. MONOTONICITY. C1 (more dose never less exposure) and C2 (longer window never
   less AUC) are expected EXACT: the PDE is linear in ``C`` and the source enters
   linearly, and ``auc`` is a cumulative integral of a non-negative field whose
   trajectory does not depend on where the window ends. C3 is pre-registered as
   **EXPECTED FALSE AS STATED**, from a derivation: at the fixed point,
   ``1^T (K - L + Adv) c* = 1^T S`` reduces (zero row sums again) to
   ``sum(k*c*) = sum(S)``, so with uniform ``k`` the TOTAL steady exposure is
   ``sum(S)/k`` -- **exactly independent of the diffusivity**. If raising D
   raised the exposure everywhere, the total would have to rise. It cannot. So
   the curves must cross, and the only property that can hold is the repaired
   one: monotone in D beyond a crossover radius ``r*``, which is measured.

D. DETERMINISM. Expected exact, and the interesting part is the search: which
   objects in the pipeline are ordered by a ``set`` or an unseeded generator. The
   run re-executes an identical spec in fresh subprocesses under four different
   ``PYTHONHASHSEED`` values (``0``, ``1``, ``12345``, ``random``) and compares
   SHA-256 of the raw bytes of every array. D4 is pre-registered as EXPECTED
   NON-DETERMINISTIC: ``cells.pkpd.apply_drug_kill`` defaults ``rng`` to
   ``np.random.default_rng()`` with no seed, so a direct caller who omits it gets
   a different kill every time. Whether that is reachable through
   ``sim.engine.simulate_response`` is a separate check and is made.

E. SCALE INVARIANCE. Re-derived here rather than read off
   ``results/reach_invariances_measure.md`` or ``results/steady_state_reach.md``.
   E1/E2 (exact proportionality to any scalar source factor -- which is what
   makes the reach map barrier-blind and dose-free -- and superposition) follow
   from linearity. E3 is the capacity identity, and it is testable WITHOUT a
   capacity term: capacity ``R`` multiplies only the time derivative, so the
   capacity-R problem at ``(T, dt)`` is the ordinary problem with
   ``(D, k, v, S)`` all divided by ``R`` and the window and step likewise -- i.e.
   scaling all four by ``s`` and ``(T, dt)`` by ``1/s`` must scale ``auc`` by
   ``1/s``. E4 is its steady-state corollary: ``(K - L + Adv) c* = S`` is
   homogeneous of degree one in that same common factor, so ``c*`` is EXACTLY
   invariant -- ``R`` is absent from the steady problem, not small in it. All
   four thresholds 1e-12 relative.

Units: length um, time s, concentration uM, diffusivity um^2/s, rates 1/s.

Run:  .venv/bin/python experiments/conservation_audit.py
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import os
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import ndimage
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import factorized

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from tmesim.core import (  # noqa: E402
    CellPopulation,
    CellState,
    CellType,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    PKConfig,
    SimConfig,
    TissueFields,
)
from tmesim.physics.diffusion import (  # noqa: E402
    _advective_divergence,
    _build_diffusion_operator,
    _effective_diffusivity,
    simulate_drug_transport,
)

# --------------------------------------------------------------------------- #
# The one constant this repo pins everywhere.
# --------------------------------------------------------------------------- #
ALPHA_PINNED = 22.223477711172563


def assert_alpha(where: str) -> float:
    """Assert ``calibrate_alpha()`` is unmoved. Called at the start AND the end."""
    from tmesim.validation.calibration import calibrate_alpha

    a = calibrate_alpha()
    assert a == ALPHA_PINNED, f"{where}: calibrate_alpha() = {a!r} != {ALPHA_PINNED!r}"
    print(f"[alpha {where}] calibrate_alpha() = {a!r}  OK")
    return a


# --------------------------------------------------------------------------- #
# PRE-REGISTRATION. Fixed before the definitive run; asserted at the start of it.
# A threshold that fails is reported FAILED. None of these is swapped.
# --------------------------------------------------------------------------- #
PREREG: Dict[str, Dict] = {
    "A1_mass_budget_backward_euler": {
        "claim": "M_final == INJECTED - REMOVED with REMOVED the right-endpoint "
                 "sum dt*sum_n sum(k*c^{n+1}); holds for any drug/grid/dt/profile/"
                 "convection/exchange",
        "metric": "|M_final - (INJECTED - REMOVED)| / INJECTED",
        "threshold": 1e-10,
        "expect": "PASS",
    },
    "A2_trapezoid_budget_does_not_close": {
        "claim": "the intuitive budget k*sum(auc) uses the WRONG quadrature and "
                 "misses by O(dt)",
        "metric": "|M_final - (INJECTED - k*sum(auc))| / INJECTED",
        "threshold": 1e-5,          # must be ABOVE this, i.e. must visibly fail
        "expect": "EXCEEDS (a pass here would falsify the A1 derivation)",
    },
    "A3_zero_clamp_injects_nothing": {
        "claim": "the solver's np.maximum(c,0) clamp is never reached on a "
                 "physical input, so it is never a silent mass source",
        "metric": "sum_n sum(max(-c_raw, 0)) / INJECTED",
        "threshold": 0.0,
        "expect": "EXACTLY 0",
    },
    "B1_nonneg_smooth_flow": {
        "claim": "A is an M-matrix so the implicit half cannot undershoot; a "
                 "smooth flow field gives no undershoot at any Peclet number",
        "metric": "min_n min(c_raw) / max(c)",
        "threshold": -1e-9,
        "expect": "EXACTLY 0",
    },
    "B2_cfl_guard_sufficient_for_positivity": {
        "claim": "the shipped guard dt*max(|vx|+|vy|)/dx <= 1 is sufficient for "
                 "positivity",
        "metric": "min(c_raw)/max(c) on the divergent-cross field",
        "threshold": -1e-3,
        "expect": "VIOLATED; derived law min(c_raw)/max(c) = -(2*CFL - 1) for "
                  "CFL > 1/2, so -1.0 at the guard's own accepted limit",
    },
    "B3_strict_bound_in_production": {
        "claim": "the derived strict bound dt*(max|vx| + max|vy|)/dx <= 1/2 "
                 "holds at the shipped dx and dt on every real section",
        "metric": "dt*(max|vx| + max|vy|)/dx",
        "threshold": 0.5,
        "expect": "measure; report per section and per dt",
    },
    "C1_dose_monotone_and_linear": {
        "claim": "more dose never less exposure anywhere, and auc is exactly "
                 "proportional to the dose scale, for every compound x section",
        "metric": "max relative deviation from an exact rescale; and "
                  "all(auc_hi >= auc_lo) elementwise",
        "threshold": 1e-12,
        "expect": "PASS",
    },
    "C2_window_monotone": {
        "claim": "a longer window never gives less AUC anywhere, for every "
                 "compound x section",
        "metric": "min(auc(T2) - auc(T1)) over the grid, T2 > T1",
        "threshold": 0.0,
        "expect": "PASS (>= 0)",
    },
    "C3_diffusivity_monotone_everywhere": {
        "claim": "higher diffusivity never makes any cell see less",
        "metric": "fraction of tissue cells where auc(D_hi) < auc(D_lo)",
        "threshold": 0.0,
        "expect": "FALSE BY DERIVATION: sum(k*c*) = sum(S) makes the total "
                  "steady exposure exactly D-independent, so the curves must "
                  "cross. Repaired property: monotone beyond a crossover r*.",
    },
    "C3R_diffusivity_monotone_beyond_crossover": {
        "claim": "beyond the measured crossover radius r*, higher diffusivity "
                 "never gives less",
        "metric": "fraction of cells at distance > r* where auc(D_hi) < auc(D_lo)",
        "threshold": 0.0,
        "vacuity_guard": 0.01,
        "expect": "PASS -- BUT this is true by the definition of r*, so it is "
                  "only worth anything if the region beyond r* is a real part "
                  "of the tissue. Any row with < 1% of tissue beyond r* is "
                  "reported VACUOUS, not PASS, and the distance-binned profile "
                  "is reported instead because it cannot be vacuous.",
    },
    "C3a_single_vessel_crossover_is_clean": {
        "claim": "in single-vessel geometry the D-crossover is a SINGLE radius "
                 "and the region beyond it is a real part of the domain",
        "metric": "sign crossings along a ray; fraction of the domain beyond r*",
        "threshold": "crossings == 1 and fraction beyond r* > 0.5",
        "expect": "PASS -- this is the geometry in which the repaired property "
                  "is worth stating",
    },
    "C3b_real_tissue_crossover": {
        "claim": "on real tissue there is a distance beyond which higher "
                 "diffusivity never gives less",
        "metric": "fraction of cells in the OUTERMOST distance band that lose "
                  "exposure when D rises; measured at two grid spacings",
        "threshold": 0.0,
        "expect": "measure. Many vessels at many spacings mean distance-to-"
                  "nearest-vessel need not order cells by near/far field, so "
                  "this may well fail; if it does, the repaired property does "
                  "not transfer off single-vessel geometry and must not be "
                  "claimed for a real section.",
    },
    "C4_kill_monotone_in_dose": {
        "claim": "at a fixed seed, kill_fraction is non-decreasing in dose",
        "metric": "min over the ladder of (kill[i+1] - kill[i])",
        "threshold": 0.0,
        "expect": "PASS (>= 0)",
    },
    "D1_same_process": {
        "claim": "same spec, same seed, bit-identical within one process",
        "metric": "SHA-256 of every array",
        "threshold": 0.0,
        "expect": "IDENTICAL",
    },
    "D2_across_restarts_and_hash_seeds": {
        "claim": "bit-identical across process restarts and independent of "
                 "PYTHONHASHSEED (i.e. of set/dict iteration order)",
        "metric": "SHA-256 of every array, 4 hash seeds, fresh subprocesses",
        "threshold": 0.0,
        "expect": "IDENTICAL",
    },
    "D4_unseeded_rng_exists": {
        "claim": "cells.pkpd.apply_drug_kill defaults rng to an UNSEEDED "
                 "np.random.default_rng(); a caller that omits rng is "
                 "non-deterministic",
        "metric": "spread of kill_fraction over repeated identical calls",
        "threshold": 0.0,
        "expect": "NON-ZERO SPREAD (the defect is reached), and NOT reachable "
                  "through sim.engine.simulate_response",
    },
    "E1_source_scale_exact": {
        "claim": "auc is exactly proportional to any scalar factor on the "
                 "source (plasma level, bbb_penetrance) -- this is what makes "
                 "the reach ratio dose-free and barrier-blind",
        "metric": "max relative deviation from an exact rescale",
        "threshold": 1e-12,
        "expect": "PASS",
    },
    "E2_superposition": {
        "claim": "auc(S1 + S2) == auc(S1) + auc(S2)",
        "metric": "max relative deviation",
        "threshold": 1e-12,
        "expect": "PASS",
    },
    "E3_capacity_time_rescale": {
        "claim": "auc(D,k,v,S; T,dt) == s * auc(sD,sk,sv,sS; T/s, dt/s) -- the "
                 "capacity identity, re-derived without a capacity term",
        "metric": "max relative deviation",
        "threshold": 1e-12,
        "expect": "PASS",
    },
    "E4_steady_state_capacity_invariant": {
        "claim": "C_ss is EXACTLY invariant to a common scaling of (D,k,v,S), "
                 "i.e. capacity R is absent from the steady problem",
        "metric": "max relative deviation",
        "threshold": 1e-12,
        "expect": "PASS",
    },
    "E5_march_reaches_steady_state": {
        "claim": "the march's fixed point is the steady solve",
        "metric": "max relative deviation of final C from C_ss at T >= 24 h",
        "threshold": 1e-12,
        "expect": "PASS",
    },
}


def assert_prereg() -> None:
    """The pre-registration is data, and is asserted rather than described."""
    assert PREREG["A1_mass_budget_backward_euler"]["threshold"] == 1e-10
    assert PREREG["A2_trapezoid_budget_does_not_close"]["threshold"] == 1e-5
    assert PREREG["A3_zero_clamp_injects_nothing"]["threshold"] == 0.0
    assert PREREG["B1_nonneg_smooth_flow"]["threshold"] == -1e-9
    assert PREREG["B2_cfl_guard_sufficient_for_positivity"]["threshold"] == -1e-3
    assert PREREG["B3_strict_bound_in_production"]["threshold"] == 0.5
    assert PREREG["C1_dose_monotone_and_linear"]["threshold"] == 1e-12
    assert PREREG["C3R_diffusivity_monotone_beyond_crossover"]["vacuity_guard"] == 0.01
    for key in ("E1_source_scale_exact", "E2_superposition",
                "E3_capacity_time_rescale", "E4_steady_state_capacity_invariant",
                "E5_march_reaches_steady_state"):
        assert PREREG[key]["threshold"] == 1e-12, key
    print(f"[prereg] {len(PREREG)} expectations fixed and asserted")


# --------------------------------------------------------------------------- #
# The solver mirror.
# --------------------------------------------------------------------------- #
def march_mirror(
    micro: Microenvironment,
    drug: DrugParams,
    source: np.ndarray,
    pde: PDEConfig,
    velocity: Optional[np.ndarray] = None,
    source_profile: Optional[np.ndarray] = None,
    exchange_k: Optional[np.ndarray] = None,
    clamp: bool = True,
) -> Dict:
    """Re-run the shipped march with per-step accounting and an optional clamp.

    Uses ``diffusion``'s OWN ``_effective_diffusivity``, ``_build_diffusion_operator``
    and ``_advective_divergence``, and reproduces the loop of
    :func:`~tmesim.physics.diffusion.simulate_drug_transport` statement for
    statement. With ``clamp=True`` it is bit-identical to that function (asserted
    by :func:`assert_mirror_is_production`); with ``clamp=False`` it is the same
    scheme with ``np.maximum(c_new, 0.0)`` removed, which is the only way to see
    what the scheme actually produced before the shipped code hid it.

    Returns a dict with ``c`` (final field), ``auc``, ``injected``, ``removed``
    (the RIGHT-ENDPOINT loss sum, see module docstring sec. A), ``clamp_mass``
    (total mass the clamp created), ``min_raw`` (most negative value the scheme
    produced at any step) and ``max_c``.
    """
    domain = micro.domain
    ny, nx = domain.shape
    dx = float(domain.dx)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)

    source = np.asarray(source, dtype=np.float64)
    dt = float(pde.dt)

    d_eff = _effective_diffusivity(drug, micro.fields.ecm_density)
    k = np.full((ny, nx), float(drug.clearance) + float(drug.binding))
    if exchange_k is not None:
        k = k + np.asarray(exchange_k, dtype=np.float64)
    k[~tissue] = 0.0
    src = np.where(tissue, source, 0.0)

    use_convection = bool(pde.convection)
    if velocity is None:
        velocity = getattr(micro.fields, "velocity", None)
    if velocity is None:
        use_convection = False
    if use_convection:
        velocity = np.asarray(velocity, dtype=np.float64)
        vx, vy = velocity[0], velocity[1]
    else:
        vx = vy = None

    n = ny * nx
    lap = _build_diffusion_operator(d_eff, tissue, dx)
    a = csc_matrix((np.ones(n), (np.arange(n), np.arange(n))), shape=(n, n))
    a = a - dt * lap
    a = a + csc_matrix((dt * k.ravel(), (np.arange(n), np.arange(n))), shape=(n, n))
    solve = factorized(a.tocsc())

    n_steps = int(round(pde.total_time / dt))
    src_scale = (None if source_profile is None
                 else np.asarray(source_profile, dtype=np.float64))

    c = domain.zeros()
    auc = domain.zeros()
    injected = 0.0
    removed = 0.0
    clamp_mass = 0.0
    min_raw = 0.0
    max_c = 0.0
    for step in range(1, n_steps + 1):
        rhs = c.copy()
        if src_scale is None:
            rhs += dt * src
            injected += dt * float(src.sum())
        else:
            rhs += (dt * src_scale[step - 1]) * src
            injected += float((dt * src_scale[step - 1]) * src.sum())
        if use_convection:
            rhs -= dt * _advective_divergence(c, vx, vy, tissue, dx)
        c_new = solve(rhs.ravel()).reshape(ny, nx)

        raw_min = float(c_new.min())
        if raw_min < min_raw:
            min_raw = raw_min
        clamp_mass += -float(np.minimum(c_new, 0.0).sum())

        if clamp:
            np.maximum(c_new, 0.0, out=c_new)
        c_new[~tissue] = 0.0

        removed += dt * float((k * c_new).sum())
        auc += 0.5 * (c + c_new) * dt
        c = c_new
        max_c = max(max_c, float(c.max()))

    return {
        "c": c, "auc": auc, "injected": injected, "removed": removed,
        "clamp_mass": clamp_mass, "min_raw": min_raw, "max_c": max_c,
        "n_steps": n_steps, "k_uniform": float(drug.clearance) + float(drug.binding),
    }


def assert_mirror_is_production(verbose: bool = True) -> None:
    """Abort unless the mirror reproduces the shipped solver BIT-FOR-BIT.

    Four paths, because they are four different statements in the solver:
    static source / time-resolved profile, crossed with convection off / on.
    """
    rng = np.random.default_rng(20260920)
    ny = nx = 34
    dx = 15.0
    micro = flat_micro(ny, nx, dx, ecm=0.37)
    drug = DrugParams(name="mirror-probe", diffusivity=280.0,
                      clearance=2.3e-4, binding=6.1e-4)
    src = micro.domain.zeros()
    src[::7, ::8] = 1.3
    checked = 0
    for convection in (False, True):
        vel = None
        if convection:
            vel = np.zeros((2, ny, nx))
            vel[0] = 0.03 * rng.standard_normal((ny, nx))
            vel[1] = 0.03 * rng.standard_normal((ny, nx))
        pde = PDEConfig(dt=9.0, total_time=900.0, convection=convection,
                        save_every=17)
        n_steps = int(round(pde.total_time / pde.dt))
        for profile in (None, np.abs(rng.standard_normal(n_steps)) + 0.05):
            prod = simulate_drug_transport(micro, drug, src, pde, velocity=vel,
                                           source_profile=profile)
            mir = march_mirror(micro, drug, src, pde, velocity=vel,
                               source_profile=profile, clamp=True)
            assert np.array_equal(prod.final_concentration, mir["c"]), (
                f"mirror diverged from the solver (final C), convection={convection}, "
                f"profile={'yes' if profile is not None else 'no'}")
            assert np.array_equal(prod.auc, mir["auc"]), (
                f"mirror diverged from the solver (auc), convection={convection}, "
                f"profile={'yes' if profile is not None else 'no'}")
            checked += 1
    if verbose:
        print(f"[mirror] bit-identical to simulate_drug_transport on "
              f"{checked}/4 paths (final C and auc)")


# --------------------------------------------------------------------------- #
# Controlled geometry (no synthetic patient generator, no invented ground truth).
# --------------------------------------------------------------------------- #
def flat_micro(ny: int, nx: int, dx: float, ecm: float = 0.3,
               tissue_mask: Optional[np.ndarray] = None) -> Microenvironment:
    """A homogeneous slab used only to exercise the numerics."""
    domain = Domain(nx=nx, ny=ny, dx=dx, tissue_mask=tissue_mask)
    fields = TissueFields(ecm_density=np.full((ny, nx), float(ecm)),
                          vascular_density=np.zeros((ny, nx)))
    cells = CellPopulation(
        x=np.array([nx * dx / 2]), y=np.array([ny * dx / 2]),
        cell_type=np.array([int(CellType.TUMOR)]),
        state=np.array([int(CellState.VIABLE)]),
        sensitivity=np.array([0.5]),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells)


def divergent_cross_velocity(ny: int, nx: int, v: float,
                             i: int, j: int) -> np.ndarray:
    """THE MINIMAL REPRODUCING INPUT for the positivity violation (sec. B2).

    Four cells around ``(i, j)`` carry velocity ``v``, and nothing else moves::

        vx[i, j-1] = -v      vx[i, j+1] = +v
        vy[i-1, j] = -v      vy[i+1, j] = +v

    Every face of ``(i, j)`` is therefore an OUTFLOW face -- the east face
    velocity is ``(0 + v)/2 = +v/2`` and the west face is ``(-v + 0)/2 = -v/2``,
    both pointing away from ``(i, j)``, and likewise in y -- so the explicit step
    removes ``dt * 2v/dx * c`` from that one cell. But at no single cell does
    ``|vx| + |vy|`` exceed ``v``, so the shipped CFL guard sees ``CFL = dt*v/dx``
    and accepts anything up to 1. The factor of 2 between what the guard measures
    and what the stencil removes is the whole defect, and this is the smallest
    field that realises it exactly.
    """
    vel = np.zeros((2, ny, nx), dtype=np.float64)
    vel[0][i, j - 1] = -float(v)
    vel[0][i, j + 1] = +float(v)
    vel[1][i - 1, j] = -float(v)
    vel[1][i + 1, j] = +float(v)
    return vel


def cross_undershoot_law(cfl: float) -> float:
    """Derived (sec. B2): ``min(c_raw)/max(c) = -(2*CFL - 1)`` for ``CFL > 1/2``."""
    return -(2.0 * float(cfl) - 1.0) if cfl > 0.5 else 0.0


def strict_positivity_cfl(velocity: np.ndarray, dt: float, dx: float) -> float:
    """``dt*(max|vx| + max|vy|)/dx`` -- the two maxima taken INDEPENDENTLY.

    This is the bound the upwind stencil actually needs (<= 1/2 for positivity in
    the worst case; <= 1 for stability). The shipped guard uses
    ``max(|vx| + |vy|)`` -- a max over cells of the pointwise sum -- which is the
    same number only when both maxima sit on the same cell.
    """
    vel = np.asarray(velocity, dtype=np.float64)
    return float(dt) * (float(np.abs(vel[0]).max()) + float(np.abs(vel[1]).max())) / float(dx)


def shipped_cfl(velocity: np.ndarray, dt: float, dx: float) -> float:
    """The guard's own number, reproduced from ``diffusion.simulate_drug_transport``."""
    vel = np.asarray(velocity, dtype=np.float64)
    return float(np.max(np.abs(vel[0]) + np.abs(vel[1]))) * float(dt) / float(dx)


# --------------------------------------------------------------------------- #
# Real sections and real compounds.
# --------------------------------------------------------------------------- #
#: The five sections the repo can currently load, with the reconstruction call
#: each needs. dx is coarse so the whole audit runs in minutes; every invariant
#: here is a statement about the SCHEME, not about resolution.
SECTIONS: List[Tuple[str, Dict]] = [
    ("xenium_319", dict(source="xenium", dx_um=80.0)),
    ("xenium_io", dict(source="xenium", sample="data_cache/xenium_gbm_io", dx_um=80.0)),
    ("gse324268", dict(source="xenium", sample="data_cache/gse324268_xenium", dx_um=80.0)),
    ("visium_10x_parent_gbm", dict(source="visium", dx_um=100.0)),
    ("visium_gbm1", dict(source="visium", sample="GBM1", dx_um=100.0)),
]

_SECTION_CACHE: Dict[str, Tuple[Microenvironment, np.ndarray]] = {}
_DRUG_CACHE: Dict[str, DrugParams] = {}


def load_section(name: str, kwargs: Dict) -> Tuple[Microenvironment, np.ndarray]:
    """Reconstruct a section once and attach its IFP/velocity fields."""
    if name in _SECTION_CACHE:
        return _SECTION_CACHE[name]
    from tmesim.physics.pressure import solve_ifp
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    micro = reconstruct_real(**kwargs)
    ifp, velocity = solve_ifp(micro)
    micro.fields.ifp = ifp
    micro.fields.velocity = velocity
    _SECTION_CACHE[name] = (micro, velocity)
    return _SECTION_CACHE[name]


def compiled_drugs() -> Dict[str, DrugParams]:
    """All 12 GBM compounds through the flagship compile path the panel uses."""
    if _DRUG_CACHE:
        return _DRUG_CACHE
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.validation.calibration import calibrated_drugparams

    for name in sorted(library.GBM_DRUGS):
        dp = calibrated_drugparams(name)
        _DRUG_CACHE[name] = apply_barrier_to_drugparams(dp, library.get_molecule(name))
    return _DRUG_CACHE


def vessel_distance_field(micro: Microenvironment, source: np.ndarray,
                          quantile: float = 0.90) -> np.ndarray:
    """Distance (um) to the nearest high-delivery cell.

    "Vessel" is defined as the top decile of the delivery source over tissue,
    which is a definition internal to this audit and stated so it can be argued
    with. It is used ONLY to order cells by distance for the crossover
    measurement (sec. C3); no reported invariant depends on where the cut sits.
    """
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    src = np.asarray(source, dtype=np.float64)
    vals = src[tissue]
    if vals.size == 0 or float(vals.max()) <= 0.0:
        return np.full(src.shape, np.nan)
    thr = float(np.quantile(vals, quantile))
    vessel = (src >= thr) & tissue
    if not vessel.any():
        return np.full(src.shape, np.nan)
    return ndimage.distance_transform_edt(~vessel, sampling=float(micro.domain.dx))


def _relmax(a: np.ndarray, b: np.ndarray) -> float:
    """max |a - b| / |b| over the entries where ``b`` is non-zero. 0 if all zero."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    den = np.abs(b)
    mask = den > 0.0
    if not mask.any():
        return 0.0 if np.array_equal(a, b) else float("inf")
    return float((np.abs(a - b)[mask] / den[mask]).max())


# --------------------------------------------------------------------------- #
# Stage A -- mass balance.
# --------------------------------------------------------------------------- #
def stage_mass_balance() -> Dict:
    """A1/A2/A3 over a ladder of drugs, grids, dt, profiles and flags, then on
    every real section."""
    print("\n=== STAGE A -- MASS BALANCE ===")
    assert_mirror_is_production()
    rows: List[Dict] = []

    rng = np.random.default_rng(1)
    ny = nx = 36
    dx = 15.0
    holes = np.ones((ny, nx), dtype=bool)
    holes[8:14, 9:17] = False          # an internal non-tissue hole: no-flux inside
    ladder = [
        # (label, D, clearance, binding, dt, T, convection, profile?, mask, ecm)
        ("inert, no sink", 300.0, 0.0, 0.0, 10.0, 1800.0, False, False, None, 0.30),
        ("reactive", 300.0, 2e-4, 5e-4, 10.0, 1800.0, False, False, None, 0.30),
        ("fine dt", 300.0, 2e-4, 5e-4, 0.5, 600.0, False, False, None, 0.30),
        ("coarse dt", 300.0, 2e-4, 5e-4, 60.0, 7200.0, False, False, None, 0.30),
        ("fast sink", 50.0, 1e-2, 1e-2, 10.0, 1800.0, False, False, None, 0.55),
        ("slow diffusion", 3.0, 3e-2, 1e-3, 5.0, 1200.0, False, False, None, 0.80),
        ("convection", 300.0, 2e-4, 5e-4, 10.0, 1800.0, True, False, None, 0.30),
        ("time-resolved", 300.0, 2e-4, 5e-4, 10.0, 1800.0, True, True, None, 0.30),
        ("internal hole", 300.0, 2e-4, 5e-4, 10.0, 1800.0, True, True, holes, 0.30),
        ("zero ecm", 300.0, 2e-4, 5e-4, 10.0, 1800.0, False, False, None, 0.0),
    ]
    for (label, D, cl, bd, dt, T, conv, prof, mask, ecm) in ladder:
        micro = flat_micro(ny, nx, dx, ecm=ecm, tissue_mask=mask)
        drug = DrugParams(name=label, diffusivity=D, clearance=cl, binding=bd)
        src = micro.domain.zeros()
        src[::7, ::9] = 1.5
        vel = None
        if conv:
            vel = np.zeros((2, ny, nx))
            vel[0] = 0.04 * rng.standard_normal((ny, nx))
            vel[1] = 0.04 * rng.standard_normal((ny, nx))
        n_steps = int(round(T / dt))
        profile = (np.abs(rng.standard_normal(n_steps)) + 0.02) if prof else None
        pde = PDEConfig(dt=dt, total_time=T, convection=conv, save_every=10 ** 9)
        m = march_mirror(micro, drug, src, pde, velocity=vel,
                         source_profile=profile, clamp=True)
        rows.append(_budget_row(label, m, cl, bd))

    # exchange term on (adds to k; the budget must absorb it identically)
    micro = flat_micro(ny, nx, dx, ecm=0.3)
    drug = DrugParams(name="exchange", diffusivity=300.0, clearance=2e-4, binding=5e-4)
    src = micro.domain.zeros()
    src[::7, ::9] = 1.5
    ex = np.full((ny, nx), 4e-3)
    pde = PDEConfig(dt=10.0, total_time=1800.0, convection=False, save_every=10 ** 9)
    m = march_mirror(micro, drug, src, pde, exchange_k=ex, clamp=True)
    rows.append(_budget_row("transvascular exchange ON", m,
                            drug.clearance, drug.binding, exchange=4e-3))

    # real sections x a few compounds
    print("  -- real sections --")
    from tmesim.treatment.drugs import delivery_source
    drugs = compiled_drugs()
    for sec_name, kwargs in SECTIONS:
        micro, vel = load_section(sec_name, kwargs)
        for dname in ("temozolomide", "doxorubicin", "lomustine"):
            drug = drugs[dname]
            src = delivery_source(micro, drug, 1.0)
            pde = PDEConfig(dt=10.0, total_time=1800.0, convection=True,
                            save_every=10 ** 9)
            m = march_mirror(micro, drug, src, pde, velocity=vel, clamp=True)
            rows.append(_budget_row(f"{sec_name}/{dname}", m,
                                    drug.clearance, drug.binding))

    worst_be = max(r["rel_resid_backward_euler"] for r in rows)
    worst_trap = max(r["rel_resid_trapezoid"] for r in rows)
    best_trap = min(r["rel_resid_trapezoid"] for r in rows if r["k"] > 0)
    worst_clamp = max(r["clamp_mass_over_injected"] for r in rows)
    verdict = {
        "A1": "PASS" if worst_be < PREREG["A1_mass_budget_backward_euler"]["threshold"] else "FAILED",
        "A2": ("CONFIRMED (does not close)"
               if best_trap > PREREG["A2_trapezoid_budget_does_not_close"]["threshold"]
               else "FAILED (it closed, so the A1 derivation is wrong)"),
        "A3": "PASS" if worst_clamp == 0.0 else "VIOLATED",
    }
    for r in rows:
        print(f"  {r['label']:34s} relBE={r['rel_resid_backward_euler']:.3e} "
              f"relTRAP={r['rel_resid_trapezoid']:.3e} clamp/inj={r['clamp_mass_over_injected']:.3e}")
    print(f"  worst backward-Euler residual  = {worst_be:.3e}  "
          f"(threshold {PREREG['A1_mass_budget_backward_euler']['threshold']:.0e}) -> {verdict['A1']}")
    print(f"  trapezoid residual range over reactive rows = "
          f"[{best_trap:.3e}, {worst_trap:.3e}] -> {verdict['A2']}")
    print(f"  worst clamp mass / injected    = {worst_clamp:.3e} -> {verdict['A3']}")
    return {"rows": rows, "worst_backward_euler": worst_be,
            "trapezoid_min": best_trap, "trapezoid_max": worst_trap,
            "worst_clamp_over_injected": worst_clamp, "verdict": verdict}


def _budget_row(label: str, m: Dict, clearance: float, binding: float,
                exchange: float = 0.0) -> Dict:
    mass = float(m["c"].sum())
    inj = m["injected"]
    k = m["k_uniform"] + exchange
    resid_be = mass - (inj - m["removed"])
    resid_trap = mass - (inj - k * float(m["auc"].sum()))
    denom = inj if inj > 0 else 1.0
    share = {}
    if k > 0:
        share = {
            "cleared": m["removed"] * clearance / k,
            "bound": m["removed"] * binding / k,
            "exchanged": m["removed"] * exchange / k,
        }
    return {
        "label": label, "k": k, "n_steps": m["n_steps"],
        "injected": inj, "removed_backward_euler": m["removed"],
        "final_mass": mass,
        "rel_resid_backward_euler": abs(resid_be) / denom,
        "rel_resid_trapezoid": abs(resid_trap) / denom,
        "clamp_mass_over_injected": m["clamp_mass"] / denom,
        "min_raw": m["min_raw"],
        "removed_split": share,
    }


# --------------------------------------------------------------------------- #
# Stage B -- non-negativity.
# --------------------------------------------------------------------------- #
PECLET_LADDER = (0.0, 0.1, 1.0, 10.0, 100.0, 1000.0)
CFL_LADDER = (0.4, 0.5, 0.6, 0.75, 0.9, 1.0)


def stage_nonnegativity() -> Dict:
    """B1 (smooth-flow Peclet ladder), B2 (the constructed counterexample) and
    B3 (is the derived bound crossed in production?)."""
    print("\n=== STAGE B -- NON-NEGATIVITY ===")
    assert_mirror_is_production(verbose=False)

    ny = nx = 40
    dx = 15.0
    ecm = 0.35
    micro = flat_micro(ny, nx, dx, ecm=ecm)
    drug = DrugParams(name="ladder", diffusivity=300.0, clearance=2e-4, binding=5e-4)
    d_eff = float(_effective_diffusivity(drug, micro.fields.ecm_density).max())
    src = micro.domain.zeros()
    src[20, 20] = 2.0

    smooth: List[Dict] = []
    print(f"  -- B1: Peclet ladder, uniform flow, CFL held at 0.9 (D_eff={d_eff:g}) --")
    for pe in PECLET_LADDER:
        v = pe * d_eff / dx
        vel = None
        dt = 1.0
        if v > 0:
            vel = np.zeros((2, ny, nx))
            vel[0] = v
            dt = 0.9 * dx / v
        pde = PDEConfig(dt=dt, total_time=200 * dt, convection=v > 0,
                        save_every=10 ** 9)
        m = march_mirror(micro, drug, src, pde, velocity=vel, clamp=False)
        rel = m["min_raw"] / m["max_c"] if m["max_c"] > 0 else 0.0
        smooth.append({"kind": "uniform", "peclet": pe, "v_um_s": v, "dt": dt,
                       "min_raw": m["min_raw"], "max_c": m["max_c"],
                       "rel_undershoot": rel})
        print(f"    Pe={pe:8.1f} v={v:11.3f} dt={dt:.4g} min_raw={m['min_raw']:.3e} "
              f"rel={rel:.3e}")

    print("  -- B1b: Peclet ladder, smooth diverging radial flow, CFL 0.9 and 0.99 --")
    X, Y = micro.domain.coords_um()
    cx, cy = X.mean(), Y.mean()
    radius = max(float(np.abs(X - cx).max()), float(np.abs(Y - cy).max()))
    for pe in PECLET_LADDER[1:]:
        vmax = pe * d_eff / dx
        vel = np.zeros((2, ny, nx))
        vel[0] = vmax * (X - cx) / radius
        vel[1] = vmax * (Y - cy) / radius
        speed = float(np.max(np.abs(vel[0]) + np.abs(vel[1])))
        for target in (0.9, 0.99):
            dt = target * dx / speed
            pde = PDEConfig(dt=dt, total_time=200 * dt, convection=True,
                            save_every=10 ** 9)
            m = march_mirror(micro, drug, src, pde, velocity=vel, clamp=False)
            rel = m["min_raw"] / m["max_c"] if m["max_c"] > 0 else 0.0
            smooth.append({"kind": "radial-diverging", "peclet": pe,
                           "cfl": target, "dt": dt, "min_raw": m["min_raw"],
                           "max_c": m["max_c"], "rel_undershoot": rel})
            print(f"    Pe={pe:8.1f} CFL={target} dt={dt:.4g} "
                  f"min_raw={m['min_raw']:.3e} rel={rel:.3e}")

    worst_smooth = min(r["rel_undershoot"] for r in smooth)
    b1 = ("PASS" if worst_smooth >= PREREG["B1_nonneg_smooth_flow"]["threshold"]
          else "VIOLATED")
    print(f"  B1 worst relative undershoot over the smooth ladder = {worst_smooth:.3e}"
          f" -> {b1}")

    # ---------------- B2: the constructed counterexample -------------------- #
    print("  -- B2: divergent-cross field, the guard's own accepted CFL range --")
    cross: List[Dict] = []
    for (gy, gx, gdx, v, D) in [(9, 9, 15.0, 30.0, 1e-6),
                                (11, 13, 7.5, 2.0, 1e-6),
                                (21, 21, 25.0, 100.0, 1e-6),
                                (9, 9, 15.0, 30.0, 1e-3),
                                (9, 9, 15.0, 30.0, 1e-1)]:
        cm = flat_micro(gy, gx, gdx, ecm=0.0)
        cd = DrugParams(name="tracer", diffusivity=D, clearance=0.0, binding=0.0)
        i, j = gy // 2, gx // 2
        vel = divergent_cross_velocity(gy, gx, v, i, j)
        csrc = cm.domain.zeros()
        csrc[i, j] = 1.0
        speed = float(np.max(np.abs(vel[0]) + np.abs(vel[1])))
        for cfl in CFL_LADDER:
            dt = cfl * gdx / speed
            n_steps = 8
            pde = PDEConfig(dt=dt, total_time=n_steps * dt, convection=True,
                            save_every=10 ** 9)
            profile = np.zeros(n_steps)
            profile[0] = 1.0            # one bolus, then the source is off
            accepted = True
            try:
                simulate_drug_transport(cm, cd, csrc, pde, velocity=vel,
                                        source_profile=profile)
            except ValueError:
                accepted = False
            raw = march_mirror(cm, cd, csrc, pde, velocity=vel,
                               source_profile=profile, clamp=False)
            clamped = march_mirror(cm, cd, csrc, pde, velocity=vel,
                                   source_profile=profile, clamp=True)
            rel = raw["min_raw"] / raw["max_c"] if raw["max_c"] > 0 else 0.0
            cross.append({
                "grid": f"{gy}x{gx}", "dx": gdx, "v": v, "D": D, "cfl": cfl,
                "dt": dt, "solver_accepts": accepted,
                "min_raw": raw["min_raw"], "max_c": raw["max_c"],
                "rel_undershoot": rel, "law": cross_undershoot_law(cfl),
                "clamp_mass_over_injected": (clamped["clamp_mass"]
                                             / clamped["injected"]),
                "shipped_cfl": shipped_cfl(vel, dt, gdx),
                "strict_cfl": strict_positivity_cfl(vel, dt, gdx),
            })
    for r in cross:
        if r["D"] == 1e-6 and r["grid"] == "9x9":
            print(f"    CFL={r['cfl']:4.2f} accepted={r['solver_accepts']} "
                  f"rel={r['rel_undershoot']:+.6f} law={r['law']:+.6f} "
                  f"clamp/inj={r['clamp_mass_over_injected']:.4f}")
    tracer = [r for r in cross if r["D"] == 1e-6]
    law_err = max(abs(r["rel_undershoot"] - r["law"]) for r in tracer)
    worst_cross = min(r["rel_undershoot"] for r in cross)
    accepted_and_negative = [r for r in cross
                             if r["solver_accepts"] and r["rel_undershoot"] < -1e-3]
    b2 = ("VIOLATED" if accepted_and_negative else "PASS")
    print(f"    derived law reproduced to {law_err:.3e} (D->0 rows); worst "
          f"relative undershoot {worst_cross:.6f}")
    print(f"    {len(accepted_and_negative)}/{len(cross)} rows are ACCEPTED by the "
          f"shipped CFL guard and undershoot past {PREREG['B2_cfl_guard_sufficient_for_positivity']['threshold']:.0e}"
          f" -> B2 {b2}")

    # ---------------- B3: is it reached in production? ---------------------- #
    # The velocity magnitude rises as the grid is refined, and the shipped
    # PDEConfig() default dt is 0.5 s while app/server/jobs.py runs PDE_DT_S =
    # 10.0 s at whatever dx the run spec asks for (12 um is the grid the
    # steady-state and resolution work calls the shipping one). So B3 is a
    # statement about (section, dx, dt) and is measured on all three axes, not
    # only at the coarse dx this audit uses elsewhere.
    print("  -- B3: the strict bound across (section, dx, dt) --")
    prod: List[Dict] = []
    arms: List[Tuple[str, Dict]] = list(SECTIONS) + [
        ("xenium_319 @ dx=30", dict(source="xenium", dx_um=30.0)),
        ("xenium_319 @ dx=12 (shipping)", dict(source="xenium", dx_um=12.0)),
    ]
    for sec_name, kwargs in arms:
        micro_r, vel_r = load_section(sec_name, kwargs)
        dxr = float(micro_r.domain.dx)
        for dt in (0.5, 10.0):
            prod.append({
                "section": sec_name, "dx_um": dxr, "dt_s": dt,
                "shipped_cfl": shipped_cfl(vel_r, dt, dxr),
                "strict_cfl": strict_positivity_cfl(vel_r, dt, dxr),
                "max_speed_um_s": float(np.max(np.abs(vel_r[0]) + np.abs(vel_r[1]))),
            })
            print(f"    {sec_name:30s} dx={dxr:5.1f} dt={dt:5.1f} "
                  f"shipped_CFL={prod[-1]['shipped_cfl']:.4f} "
                  f"strict_CFL={prod[-1]['strict_cfl']:.4f}")

    # and: does the field on real tissue actually undershoot at the shipped dt?
    print("  -- B3b: raw (pre-clamp) undershoot on real tissue at dt = 0.5 and 10 --")
    from tmesim.treatment.drugs import delivery_source
    drugs = compiled_drugs()
    real_rows: List[Dict] = []
    for sec_name, kwargs in arms:
        micro_r, vel_r = load_section(sec_name, kwargs)
        for dname in ("temozolomide", "lomustine"):
            drug_r = drugs[dname]
            src_r = delivery_source(micro_r, drug_r, 1.0)
            for dt in (0.5, 10.0):
                n_steps = 60 if float(micro_r.domain.dx) >= 50.0 else 25
                pde = PDEConfig(dt=dt, total_time=n_steps * dt, convection=True,
                                save_every=10 ** 9)
                profile = np.zeros(n_steps)
                profile[0] = 50.0            # bolus then nothing: the worst case
                m = march_mirror(micro_r, drug_r, src_r, pde, velocity=vel_r,
                                 source_profile=profile, clamp=False)
                rel = m["min_raw"] / m["max_c"] if m["max_c"] > 0 else 0.0
                real_rows.append({"section": sec_name, "drug": dname, "dt": dt,
                                  "min_raw": m["min_raw"], "rel_undershoot": rel})
    worst_real = min(r["rel_undershoot"] for r in real_rows)
    print(f"    worst relative undershoot on real tissue = {worst_real:.3e}")

    crossed = [r for r in prod if r["strict_cfl"] > PREREG["B3_strict_bound_in_production"]["threshold"]]
    b3 = ("PASS (margin intact)" if not crossed
          else f"MARGIN LOST on {len(crossed)}/{len(prod)} (section, dt) pairs")
    print(f"  B3 -> {b3}")
    return {"smooth": smooth, "cross": cross, "production": prod,
            "real_undershoot": real_rows,
            "worst_smooth_rel_undershoot": worst_smooth,
            "worst_cross_rel_undershoot": worst_cross,
            "law_max_error": law_err, "worst_real_rel_undershoot": worst_real,
            "verdict": {"B1": b1, "B2": b2, "B3": b3}}


# --------------------------------------------------------------------------- #
# Stage C -- monotonicity.
# --------------------------------------------------------------------------- #
DOSE_LADDER = (1e-3, 0.1, 1.0, 10.0, 1e3)
WINDOW_LADDER = (600.0, 1800.0, 3600.0)
D_LADDER = (0.25, 0.5, 1.0, 2.0, 4.0)


def stage_monotonicity() -> Dict:
    """C1/C2/C3 over every compound x every section; C4 on the engine."""
    print("\n=== STAGE C -- MONOTONICITY (12 compounds x 5 sections) ===")
    from tmesim.treatment.drugs import delivery_source

    drugs = compiled_drugs()
    dose_rows: List[Dict] = []
    window_rows: List[Dict] = []
    diff_rows: List[Dict] = []
    degenerate: List[str] = []

    for sec_name, kwargs in SECTIONS:
        micro, vel = load_section(sec_name, kwargs)
        tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
        t0 = time.time()
        for dname in sorted(drugs):
            drug = drugs[dname]
            src = delivery_source(micro, drug, 1.0)
            if float(np.abs(src).max()) == 0.0:
                # bbb_penetrance == 0 -> the source is identically zero, so every
                # field is exactly zero and the monotonicity questions are
                # degenerate rather than passing. Recorded, not silently skipped.
                degenerate.append(f"{sec_name}/{dname}")
                continue

            # --- C1: dose ---
            pde = PDEConfig(dt=10.0, total_time=1800.0, convection=True,
                            save_every=10 ** 9)
            base = simulate_drug_transport(micro, drug, src, pde, velocity=vel).auc
            prev = None
            worst_lin = 0.0
            ordering_ok = True
            for lam in DOSE_LADDER:
                a = simulate_drug_transport(micro, drug, lam * src, pde,
                                            velocity=vel).auc
                worst_lin = max(worst_lin, _relmax(a, lam * base))
                if prev is not None:
                    ordering_ok = ordering_ok and bool(np.all(a >= prev))
                prev = a
            dose_rows.append({"section": sec_name, "drug": dname,
                              "max_rel_dev_from_rescale": worst_lin,
                              "pointwise_ordering_ok": ordering_ok})

            # --- C2: window ---
            aucs = []
            for T in WINDOW_LADDER:
                p = PDEConfig(dt=10.0, total_time=T, convection=True,
                              save_every=10 ** 9)
                aucs.append(simulate_drug_transport(micro, drug, src, p,
                                                    velocity=vel).auc)
            min_gap = min(float((aucs[i + 1] - aucs[i]).min())
                          for i in range(len(aucs) - 1))
            # The grid-wide minimum is exactly 0 on every non-tissue cell (both
            # windows hold those at 0), which would make "never decreases" true
            # for a vacuous reason. The tissue-restricted minimum is the one
            # that can actually go negative, so both are recorded.
            min_gap_tissue = min(float((aucs[i + 1] - aucs[i])[tissue].min())
                                 for i in range(len(aucs) - 1))
            window_rows.append({"section": sec_name, "drug": dname,
                                "min_increment": min_gap,
                                "min_increment_tissue": min_gap_tissue,
                                "ok": min_gap >= 0.0 and min_gap_tissue >= 0.0})

            # --- C3: diffusivity ---
            dist = vessel_distance_field(micro, src)
            d_aucs = []
            for f in D_LADDER:
                dr = dataclasses.replace(drug, diffusivity=drug.diffusivity * f)
                p = PDEConfig(dt=10.0, total_time=1800.0, convection=True,
                              save_every=10 ** 9)
                d_aucs.append(simulate_drug_transport(micro, dr, src, p,
                                                      velocity=vel).auc)
            totals = [float(a[tissue].sum()) for a in d_aucs]
            total_spread = ((max(totals) - min(totals)) / max(totals)
                            if max(totals) > 0 else 0.0)
            worst_frac_down = 0.0
            crossovers: List[float] = []
            for i in range(len(d_aucs) - 1):
                lo, hi = d_aucs[i], d_aucs[i + 1]
                down = (hi < lo) & tissue
                worst_frac_down = max(worst_frac_down,
                                      float(down.sum()) / float(tissue.sum()))
                # r* = the OUTERMOST radius at which raising D still takes drug
                # away. Beyond it, nothing decreases, by construction of r*.
                d_down = dist[down]
                crossovers.append(float(d_down.max()) if d_down.size else 0.0)
            r_star = max(crossovers) if crossovers else float("nan")
            beyond = tissue & (dist > r_star)
            n_beyond = int(beyond.sum())
            # The repaired property "monotone beyond r*" is TRUE BY THE
            # DEFINITION of r* and is therefore worthless unless the region it
            # speaks about is a real part of the tissue. The share of tissue
            # beyond r* is recorded next to it so the claim can be read at the
            # strength it has; where that share is ~0 the row is VACUOUS and is
            # reported as vacuous, not as a pass.
            frac_tissue_beyond = n_beyond / float(tissue.sum())
            frac_down_beyond = 0.0
            if n_beyond:
                for i in range(len(d_aucs) - 1):
                    down = (d_aucs[i + 1] < d_aucs[i]) & beyond
                    frac_down_beyond = max(frac_down_beyond,
                                           float(down.sum()) / float(n_beyond))
            # The same question asked so it cannot be vacuous: the fraction of
            # cells that DECREASE, binned by distance, out to the 95th
            # percentile of the distance field.
            edges = np.quantile(dist[tissue], np.linspace(0.0, 1.0, 11))
            profile = []
            for b in range(10):
                band = tissue & (dist >= edges[b]) & (dist <= edges[b + 1])
                if not band.any():
                    profile.append(None)
                    continue
                worst = 0.0
                for i in range(len(d_aucs) - 1):
                    down = (d_aucs[i + 1] < d_aucs[i]) & band
                    worst = max(worst, float(down.sum()) / float(band.sum()))
                profile.append({"d_lo_um": float(edges[b]),
                                "d_hi_um": float(edges[b + 1]),
                                "n": int(band.sum()),
                                "fraction_decreasing": worst})
            diff_rows.append({
                "section": sec_name, "drug": dname,
                "L_um": float(np.sqrt(drug.diffusivity
                                      / (drug.clearance + drug.binding))),
                "total_auc_spread_over_D": total_spread,
                "worst_fraction_decreasing": worst_frac_down,
                "crossover_r_star_um": r_star,
                "n_cells_beyond_r_star": n_beyond,
                "fraction_of_tissue_beyond_r_star": frac_tissue_beyond,
                "fraction_decreasing_beyond_r_star": frac_down_beyond,
                "distance_profile": profile,
            })
        print(f"  {sec_name}: {len(drugs)} compounds in {time.time() - t0:.1f}s")

    worst_lin = max(r["max_rel_dev_from_rescale"] for r in dose_rows)
    all_ordered = all(r["pointwise_ordering_ok"] for r in dose_rows)
    worst_window = min(r["min_increment"] for r in window_rows)
    worst_window_tissue = min(r["min_increment_tissue"] for r in window_rows)
    worst_down = max(r["worst_fraction_decreasing"] for r in diff_rows)
    worst_down_beyond = max(r["fraction_decreasing_beyond_r_star"] for r in diff_rows)
    worst_total_spread = max(r["total_auc_spread_over_D"] for r in diff_rows)
    r_stars = [r["crossover_r_star_um"] for r in diff_rows
               if np.isfinite(r["crossover_r_star_um"])]
    beyond_shares = [r["fraction_of_tissue_beyond_r_star"] for r in diff_rows]
    vacuous = [f"{r['section']}/{r['drug']}" for r in diff_rows
               if r["fraction_of_tissue_beyond_r_star"] < 0.01]

    print(f"  C1 worst rel dev from exact rescale = {worst_lin:.3e}; pointwise "
          f"ordering holds everywhere: {all_ordered}")
    print(f"  C2 worst AUC increment over a longer window = {worst_window:.3e} "
          f"(grid-wide), {worst_window_tissue:.3e} (tissue only)")
    print(f"  C3 worst fraction of tissue where higher D gives LESS = {worst_down:.4f}")
    print(f"     total AUC spread across the D ladder (should be ~0) = "
          f"{worst_total_spread:.3e}")
    print(f"     crossover r* range = [{min(r_stars):.1f}, {max(r_stars):.1f}] um; "
          f"fraction decreasing beyond r* = {worst_down_beyond:.4f}")
    print(f"     share of tissue beyond r* = [{min(beyond_shares):.4f}, "
          f"{max(beyond_shares):.4f}] (median "
          f"{float(np.median(beyond_shares)):.4f}); "
          f"{len(vacuous)}/{len(diff_rows)} rows are VACUOUS (< 1% of tissue "
          f"beyond r*)")
    # a real curve, not a summary: how the decrease fraction falls with distance
    rep = next(r for r in diff_rows
               if r["section"] == "xenium_319" and r["drug"] == "temozolomide")
    print("     distance profile (xenium_319 / temozolomide), fraction of cells "
          "where higher D gives LESS:")
    for band in rep["distance_profile"]:
        if band:
            print(f"       {band['d_lo_um']:7.1f}-{band['d_hi_um']:7.1f} um  "
                  f"n={band['n']:6d}  down={band['fraction_decreasing']:.4f}")
    if degenerate:
        print(f"  degenerate (source identically zero, bbb_penetrance = 0): "
              f"{sorted(set(d.split('/')[1] for d in degenerate))}")

    # --- C4: kill monotone in dose, on the engine, every compound ----------- #
    print("  -- C4: kill_fraction vs dose at a fixed seed (engine) --")
    print("     Only the SIGN of each rung is printed. The invariant is the "
          "ordering WITHIN one compound across its own dose ladder; it needs no "
          "absolute kill fraction, and absolute kill fraction is not a "
          "defensible output of this project. Compounds are listed "
          "alphabetically and nothing here compares one to another.")
    from tmesim.sim.engine import simulate_response
    micro, _ = load_section(*SECTIONS[0])
    kill_rows: List[Dict] = []
    for dname in sorted(drugs):
        ladder = []
        for peak in (1.0, 10.0, 50.0, 250.0, 1000.0):
            cfg = SimConfig(drug=drugs[dname], pk=PKConfig(plasma_peak=peak),
                            pde=PDEConfig(dt=10.0, total_time=1800.0,
                                          convection=True, save_every=180),
                            seed=3)
            ladder.append(float(simulate_response(micro, cfg).kill_fraction))
        gaps = [ladder[i + 1] - ladder[i] for i in range(len(ladder) - 1)]
        kill_rows.append({"drug": dname, "n_rungs": len(ladder),
                          "min_gap": min(gaps), "ok": min(gaps) >= 0.0})
        print(f"    {dname:14s} {len(gaps)}/{len(gaps)} rungs non-decreasing; "
              f"smallest rung {min(gaps):+.3e}; monotone={min(gaps) >= 0.0}")
    worst_kill_gap = min(r["min_gap"] for r in kill_rows)

    crossover = stage_crossover()

    verdict = {
        "C1": ("PASS" if worst_lin < PREREG["C1_dose_monotone_and_linear"]["threshold"]
               and all_ordered else "FAILED"),
        "C2": ("PASS" if worst_window >= 0.0 and worst_window_tissue >= 0.0
               else "VIOLATED"),
        "C3": ("FALSE AS STATED (as pre-registered)" if worst_down > 0.0
               else "PASS (which would falsify the derivation)"),
        "C3R": (("VACUOUS on %d/%d rows" % (len(vacuous), len(diff_rows)))
                if vacuous else
                ("PASS" if worst_down_beyond == 0.0 else "VIOLATED")),
        "C3a": crossover["verdict_single_vessel"],
        "C3b": crossover["verdict_real_tissue"],
        "C4": "PASS" if worst_kill_gap >= 0.0 else "VIOLATED",
    }
    return {"dose": dose_rows, "window": window_rows, "diffusivity": diff_rows,
            "kill": kill_rows, "degenerate": degenerate,
            "worst_linearity": worst_lin, "all_ordered": all_ordered,
            "worst_window_increment": worst_window,
            "worst_window_increment_tissue": worst_window_tissue,
            "worst_fraction_decreasing": worst_down,
            "worst_fraction_decreasing_beyond_r_star": worst_down_beyond,
            "worst_total_auc_spread_over_D": worst_total_spread,
            "r_star_min_um": min(r_stars), "r_star_max_um": max(r_stars),
            "beyond_share_min": min(beyond_shares),
            "beyond_share_max": max(beyond_shares),
            "beyond_share_median": float(np.median(beyond_shares)),
            "vacuous_rows": vacuous,
            "representative_profile": rep["distance_profile"],
            "crossover": crossover,
            "worst_kill_gap": worst_kill_gap, "verdict": verdict}


def stage_crossover() -> Dict:
    """C3a / C3b -- the crossover, asked so it CANNOT come back vacuous.

    The ``r*`` in C3 is defined as the outermost radius at which raising D still
    takes drug away, so "monotone beyond r*" is true by definition and is worth
    nothing unless the region beyond r* is a real part of the tissue. On the real
    sections it is not (median 0.12% of tissue). So the question is re-asked two
    ways that can each fail:

    C3a, single-vessel controlled geometry, at the STEADY STATE (the exact fixed
    point of the same stepper -- E5 -- so no finite-window confound): is the
    crossing SINGLE, and is the region beyond it a real fraction of the domain?

    C3b, real tissue: does the fraction of cells that LOSE exposure when D rises
    fall to zero in the outermost distance band, and does the answer survive
    refining the grid?
    """
    print("  -- C3a: single-vessel geometry, steady state: is the crossing clean? --")
    from app.server.reach import steady_state_concentration

    ny = nx = 121
    dx = 8.0
    micro = flat_micro(ny, nx, dx, ecm=0.0)
    X, Y = micro.domain.coords_um()
    r = np.hypot(X - X[60, 60], Y - Y[60, 60])
    src = micro.domain.zeros()
    src[60, 60] = 1.0

    single: List[Dict] = []
    for k in (7e-4, 5e-3):
        fields = {}
        for D in (25.0, 50.0, 100.0, 200.0, 400.0):
            drug = DrugParams(name="one-vessel", diffusivity=D,
                              clearance=0.2 * k, binding=0.8 * k)
            fields[D] = steady_state_concentration(micro, drug, src)
        ds = sorted(fields)
        # the conservation identity, measured: sum(k*C_ss) == sum(S) for every D
        identity = [float(fields[D].sum()) * k for D in ds]
        for i in range(len(ds) - 1):
            lo, hi = fields[ds[i]], fields[ds[i + 1]]
            down = hi < lo
            r_star = float(r[down].max()) if down.any() else 0.0
            frac_beyond = float((r > r_star).sum()) / r.size
            ray = hi[60, 60:] - lo[60, 60:]
            sign = np.sign(ray[np.abs(ray) > 0])
            crossings = int((np.diff(sign) != 0).sum())
            single.append({
                "k": k, "D_lo": ds[i], "D_hi": ds[i + 1],
                "L_lo_um": float(np.sqrt(ds[i] / k)),
                "r_star_um": r_star,
                "r_star_over_L_lo": r_star / float(np.sqrt(ds[i] / k)),
                "fraction_of_domain_beyond_r_star": frac_beyond,
                "sign_crossings_along_ray": crossings,
                "identity_sum_k_Css": identity[i],
            })
            print(f"    k={k:g} D {ds[i]:>5.0f}->{ds[i + 1]:<5.0f} r*={r_star:7.1f} um "
                  f"({single[-1]['r_star_over_L_lo']:.2f} x L_lo)  beyond r*: "
                  f"{frac_beyond:.4f} of the domain  crossings={crossings}  "
                  f"sum(k*C_ss)={identity[i]:.10f}")
    max_crossings = max(row["sign_crossings_along_ray"] for row in single)
    min_beyond = min(row["fraction_of_domain_beyond_r_star"] for row in single)
    identity_err = max(abs(row["identity_sum_k_Css"] - 1.0) for row in single)
    v_single = ("PASS (single crossing, non-vacuous region beyond it)"
                if max_crossings == 1 and min_beyond > 0.5
                else f"FAILED (crossings={max_crossings}, beyond={min_beyond:.4f})")
    print(f"    conservation identity sum(k*C_ss) == sum(S) to {identity_err:.3e}"
          f" -> C3a {v_single}")

    print("  -- C3b: real tissue, far band, and whether dx changes the answer --")
    from tmesim.treatment.drugs import delivery_source
    drug = compiled_drugs()["temozolomide"]
    real: List[Dict] = []
    for label, kwargs in [("dx=80", dict(source="xenium", dx_um=80.0)),
                          ("dx=30", dict(source="xenium", dx_um=30.0))]:
        m, vel = load_section(f"xenium_319 @ {label}", kwargs)
        tissue = np.asarray(m.domain.tissue_mask, dtype=bool)
        s = delivery_source(m, drug, 1.0)
        dist = vessel_distance_field(m, s)
        css = {}
        for f in (1.0, 4.0):
            dr = dataclasses.replace(drug, diffusivity=drug.diffusivity * f)
            css[f] = steady_state_concentration(m, dr, s, velocity=vel)
        down = (css[4.0] < css[1.0]) & tissue
        edges = np.quantile(dist[tissue], np.linspace(0.0, 1.0, 11))
        bands = []
        for b in range(10):
            band = tissue & (dist >= edges[b]) & (dist <= edges[b + 1])
            if band.any():
                bands.append({"d_lo_um": float(edges[b]), "d_hi_um": float(edges[b + 1]),
                              "n": int(band.sum()),
                              "fraction_decreasing": float((down & band).sum()) / float(band.sum())})
        real.append({"grid": label, "dx_um": float(m.domain.dx),
                     "overall_fraction_decreasing": float(down.sum()) / float(tissue.sum()),
                     "far_band_fraction_decreasing": bands[-1]["fraction_decreasing"],
                     "n_distinct_distances": int(np.unique(dist[tissue]).size),
                     "bands": bands})
        print(f"    {label}: overall down={real[-1]['overall_fraction_decreasing']:.4f}, "
              f"far band down={real[-1]['far_band_fraction_decreasing']:.4f}, "
              f"{real[-1]['n_distinct_distances']} distinct distances")
        for band in bands:
            print(f"       {band['d_lo_um']:7.1f}-{band['d_hi_um']:7.1f} um  "
                  f"n={band['n']:6d}  down={band['fraction_decreasing']:.4f}")
    far_nonzero = all(row["far_band_fraction_decreasing"] > 0.0 for row in real)
    dx_stable = abs(real[0]["overall_fraction_decreasing"]
                    - real[1]["overall_fraction_decreasing"]) < 0.01
    v_real = ("NO CLEAN CROSSOVER (the far band still loses exposure; "
              + ("stable under refinement)" if dx_stable else "NOT stable under refinement)")
              ) if far_nonzero else "clean crossover on real tissue"
    print(f"    C3b {v_real}")
    return {"single_vessel": single, "real_tissue": real,
            "identity_max_error": identity_err,
            "verdict_single_vessel": v_single, "verdict_real_tissue": v_real}


# --------------------------------------------------------------------------- #
# Stage D -- determinism.
# --------------------------------------------------------------------------- #
_DET_WORKER = r'''
import hashlib, json, os, sys
sys.path.insert(0, {repo!r})
import numpy as np
from tmesim.reconstruct.real_pipeline import reconstruct_real
from tmesim.sim.engine import simulate_response
from tmesim.core import SimConfig, PDEConfig
from tmesim.molecular import library
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

def h(a):
    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(a, dtype=np.float64)).tobytes()
    ).hexdigest()

micro = reconstruct_real(source="xenium", dx_um=80.0)
drug = apply_barrier_to_drugparams(calibrated_drugparams("temozolomide"),
                                   library.get_molecule("temozolomide"))
cfg = SimConfig(drug=drug,
                pde=PDEConfig(dt=10.0, total_time=1800.0, convection=True,
                              save_every=60),
                seed=7)
res = simulate_response(micro, cfg)
print("RESULT " + json.dumps({{
    "auc": h(res.transport.auc),
    "final_concentration": h(res.transport.final_concentration),
    "oxygen": h(micro.fields.oxygen),
    "vascular_density": h(micro.fields.vascular_density),
    "ecm_density": h(micro.fields.ecm_density),
    "cells_x": h(micro.cells.x),
    "cells_sensitivity": h(micro.cells.sensitivity),
    "viable_history": h(res.viable_history),
    "kill_fraction": repr(res.kill_fraction),
    "alpha": repr(calibrate_alpha()),
    "hashseed": os.environ.get("PYTHONHASHSEED"),
}}))
'''

HASH_SEEDS = ("0", "1", "12345", "random")


def stage_determinism() -> Dict:
    """D1/D2 (bit-identity within and across processes, under four hash seeds)
    and D4 (the unseeded generator: is it there, does it bite, is it reachable
    through the engine?)."""
    print("\n=== STAGE D -- DETERMINISM ===")
    script = _DET_WORKER.format(repo=_REPO)
    # The worker goes to the OS temp dir, never into the repo: this audit owns
    # exactly three files and a scratch script is not one of them.
    import tempfile
    fd, worker = tempfile.mkstemp(suffix="_det_worker.py")
    with os.fdopen(fd, "w") as fh:
        fh.write(script)

    runs: List[Dict] = []
    try:
        for seed in HASH_SEEDS:
            env = dict(os.environ)
            env["PYTHONHASHSEED"] = seed
            out = subprocess.run([sys.executable, worker], env=env,
                                 capture_output=True, text=True, timeout=900)
            line = [ln for ln in out.stdout.splitlines() if ln.startswith("RESULT ")]
            assert line, (f"determinism worker produced no RESULT line "
                          f"(PYTHONHASHSEED={seed}); stderr tail:\n"
                          + "\n".join(out.stderr.splitlines()[-8:]))
            runs.append(json.loads(line[-1][len("RESULT "):]))
            print(f"  PYTHONHASHSEED={seed:7s} auc={runs[-1]['auc'][:16]}... "
                  f"kill={runs[-1]['kill_fraction']}")
    finally:
        try:
            os.unlink(worker)
        except OSError:
            pass

    keys = [k for k in runs[0] if k != "hashseed"]
    differing = [k for k in keys if len({r[k] for r in runs}) != 1]
    d2 = "PASS" if not differing else f"VIOLATED on {differing}"
    print(f"  D2 across 4 process restarts / hash seeds: "
          f"{len(keys) - len(differing)}/{len(keys)} objects identical -> {d2}")
    assert runs[0]["alpha"] == repr(ALPHA_PINNED), runs[0]["alpha"]

    # --- D1: same process ---
    from tmesim.sim.engine import simulate_response
    micro, _ = load_section(*SECTIONS[0])
    drugs = compiled_drugs()
    cfg = SimConfig(drug=drugs["temozolomide"],
                    pde=PDEConfig(dt=10.0, total_time=1800.0, convection=True,
                                  save_every=60), seed=7)
    a = simulate_response(micro, cfg)
    b = simulate_response(micro, cfg)
    d1_ok = (np.array_equal(a.transport.auc, b.transport.auc)
             and np.array_equal(a.transport.final_concentration,
                                b.transport.final_concentration)
             and a.kill_fraction == b.kill_fraction
             and np.array_equal(a.viable_history, b.viable_history))
    print(f"  D1 same process, two calls, bit-identical: {d1_ok}")

    # --- D4: the unseeded generator ---
    from tmesim.cells.agents import update_cell_states
    from tmesim.cells.pkpd import apply_drug_kill
    unseeded: List[float] = []
    for _ in range(6):
        work = dataclasses.replace(micro, cells=micro.cells.copy())
        update_cell_states(work, a.fields["oxygen"], cfg.oxygen)
        apply_drug_kill(work, a.transport, cfg.pd, drugs["temozolomide"])  # no rng
        unseeded.append(1.0 - float(work.cells.viable_tumor_fraction()))
    spread = max(unseeded) - min(unseeded)
    engine_repeats = [float(simulate_response(micro, cfg).kill_fraction)
                      for _ in range(3)]
    reachable_via_engine = len(set(engine_repeats)) != 1
    print(f"  D4 apply_drug_kill(rng=None) over 6 identical calls: "
          f"{min(unseeded):.6f} .. {max(unseeded):.6f} (spread {spread:.3e}, "
          f"{100 * spread / max(unseeded):.2f}% of the value)")
    print(f"     reachable through simulate_response: {reachable_via_engine}")

    # --- ordering search ---
    order_note = ("searched tmesim/ and app/ for iteration over an unordered "
                  "set and for an unseeded generator; the only order-sensitive "
                  "sites found are membership tests (`x in set(...)`) inside "
                  "comprehensions that iterate an ORDERED list, which cannot "
                  "reorder anything, and every set that feeds an output is "
                  "wrapped in sorted().")
    return {"runs": runs, "differing": differing, "d1_ok": d1_ok,
            "unseeded_kills": unseeded, "unseeded_spread": spread,
            "engine_repeats": engine_repeats,
            "unseeded_reachable_via_engine": reachable_via_engine,
            "ordering_note": order_note,
            "verdict": {"D1": "PASS" if d1_ok else "VIOLATED", "D2": d2,
                        "D4": ("CONFIRMED (defect present and reached, not "
                               "reachable through simulate_response)"
                               if spread > 0.0 and not reachable_via_engine
                               else "unexpected")}}


# --------------------------------------------------------------------------- #
# Stage E -- scale invariance.
# --------------------------------------------------------------------------- #
SCALE_LADDER = (1e-3, 0.1, 0.5, 2.0, 10.0, 1e3)
CAPACITY_LADDER = (0.25, 0.5, 2.0, 4.0, 10.0)


def stage_scale_invariance() -> Dict:
    """E1-E5, re-derived here rather than read off the stored records."""
    print("\n=== STAGE E -- SCALE INVARIANCE (re-derived) ===")
    from app.server.reach import steady_state_concentration

    ny = nx = 41
    dx = 12.0
    micro = flat_micro(ny, nx, dx, ecm=0.25)
    drug = DrugParams(name="scale", diffusivity=120.0, clearance=1e-4, binding=6e-4)
    src = micro.domain.zeros()
    src[::7, ::7] = 0.8
    rng = np.random.default_rng(11)
    vel = np.zeros((2, ny, nx))
    vel[0] = 0.02 * rng.standard_normal((ny, nx))
    vel[1] = 0.02 * rng.standard_normal((ny, nx))

    pde = PDEConfig(dt=4.0, total_time=3600.0, convection=True, save_every=10 ** 9)
    base = simulate_drug_transport(micro, drug, src, pde, velocity=vel).auc

    e1: List[Dict] = []
    for lam in SCALE_LADDER:
        a = simulate_drug_transport(micro, drug, lam * src, pde, velocity=vel).auc
        e1.append({"lambda": lam, "max_rel_dev": _relmax(a, lam * base)})
        print(f"  E1 lambda={lam:<8g} max rel dev = {e1[-1]['max_rel_dev']:.3e}")
    worst_e1 = max(r["max_rel_dev"] for r in e1)

    s1 = micro.domain.zeros()
    s1[::7, ::7] = 0.8
    s2 = micro.domain.zeros()
    s2[3::9, 2::9] = 1.7
    p = PDEConfig(dt=4.0, total_time=1800.0, convection=True, save_every=10 ** 9)
    a1 = simulate_drug_transport(micro, drug, s1, p, velocity=vel).auc
    a2 = simulate_drug_transport(micro, drug, s2, p, velocity=vel).auc
    a12 = simulate_drug_transport(micro, drug, s1 + s2, p, velocity=vel).auc
    worst_e2 = _relmax(a1 + a2, a12)
    print(f"  E2 superposition max rel dev = {worst_e2:.3e}")

    e3: List[Dict] = []
    for s in CAPACITY_LADDER:
        d2 = dataclasses.replace(drug, diffusivity=drug.diffusivity * s,
                                 clearance=drug.clearance * s,
                                 binding=drug.binding * s)
        scaled = simulate_drug_transport(
            micro, d2, s * src,
            PDEConfig(dt=4.0 / s, total_time=3600.0 / s, convection=True,
                      save_every=10 ** 9),
            velocity=s * vel).auc
        e3.append({"s": s, "R_equivalent": 1.0 / s,
                   "max_rel_dev": _relmax(s * scaled, base)})
        print(f"  E3 s={s:<6g} (capacity R={1 / s:<8.4g}) max rel dev = "
              f"{e3[-1]['max_rel_dev']:.3e}")
    worst_e3 = max(r["max_rel_dev"] for r in e3)

    c0 = steady_state_concentration(micro, drug, src, velocity=vel)
    e4: List[Dict] = []
    for s in SCALE_LADDER:
        d2 = dataclasses.replace(drug, diffusivity=drug.diffusivity * s,
                                 clearance=drug.clearance * s,
                                 binding=drug.binding * s)
        c2 = steady_state_concentration(micro, d2, s * src, velocity=s * vel)
        e4.append({"s": s, "max_rel_dev": _relmax(c2, c0)})
        print(f"  E4 s={s:<8g} steady-state max rel dev = "
              f"{e4[-1]['max_rel_dev']:.3e}")
    worst_e4 = max(r["max_rel_dev"] for r in e4)

    e5: List[Dict] = []
    for T in (3600.0, 86400.0, 4 * 86400.0):
        a = simulate_drug_transport(
            micro, drug, src,
            PDEConfig(dt=10.0, total_time=T, convection=True, save_every=10 ** 9),
            velocity=vel)
        e5.append({"T_s": T, "max_rel_dev": _relmax(a.final_concentration, c0)})
        print(f"  E5 T={T / 3600:>5.1f} h  max rel dev of final C from C_ss = "
              f"{e5[-1]['max_rel_dev']:.3e}")
    worst_e5 = min(r["max_rel_dev"] for r in e5)

    tol = 1e-12
    verdict = {
        "E1": "PASS" if worst_e1 < tol else "FAILED",
        "E2": "PASS" if worst_e2 < tol else "FAILED",
        "E3": "PASS" if worst_e3 < tol else "FAILED",
        "E4": "PASS" if worst_e4 < tol else "FAILED",
        "E5": "PASS" if worst_e5 < tol else "FAILED",
    }
    return {"e1": e1, "e2": worst_e2, "e3": e3, "e4": e4, "e5": e5,
            "worst": {"E1": worst_e1, "E2": worst_e2, "E3": worst_e3,
                      "E4": worst_e4, "E5": worst_e5},
            "verdict": verdict}


# --------------------------------------------------------------------------- #
# Stage P -- the plants. A green check that cannot go red is not evidence.
# --------------------------------------------------------------------------- #
def _leaky_march(micro, drug, src, pde, eps: float) -> float:
    """The A1 march with an UNBUDGETED per-step loss of ``eps`` planted in it.

    The plant is here, in the audit's own copy of the loop -- never in
    ``tmesim/physics/diffusion.py``, which this work does not own. Returns the
    relative budget residual the A1 guard would see.
    """
    domain = micro.domain
    ny, nx = domain.shape
    dx = float(domain.dx)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    dt = float(pde.dt)
    d_eff = _effective_diffusivity(drug, micro.fields.ecm_density)
    k = np.full((ny, nx), float(drug.clearance) + float(drug.binding))
    k[~tissue] = 0.0
    s = np.where(tissue, src, 0.0)
    n = ny * nx
    lap = _build_diffusion_operator(d_eff, tissue, dx)
    a = csc_matrix((np.ones(n), (np.arange(n), np.arange(n))), shape=(n, n))
    a = a - dt * lap
    a = a + csc_matrix((dt * k.ravel(), (np.arange(n), np.arange(n))), shape=(n, n))
    solve = factorized(a.tocsc())
    c = domain.zeros()
    injected = 0.0
    removed = 0.0
    for _ in range(int(round(pde.total_time / dt))):
        rhs = c + dt * s
        injected += dt * float(s.sum())
        c_new = solve(rhs.ravel()).reshape(ny, nx)
        c_new *= (1.0 - eps)                       # <-- THE PLANT
        np.maximum(c_new, 0.0, out=c_new)
        c_new[~tissue] = 0.0
        removed += dt * float((k * c_new).sum())
        c = c_new
    return abs(float(c.sum()) - (injected - removed)) / injected


def stage_plants() -> Dict:
    """Break each invariant on purpose and confirm the guard FIRES.

    Every plant is in this file or in the INPUT handed to the shipped code.
    Nothing in ``tmesim/`` or ``app/`` is touched -- this audit has no plant
    licence there, and a guard that could only be proved by editing someone
    else's live file would be reported as unproved instead.

    Each row carries a CONTROL (the unplanted input, which must stay green) as
    well as the plant, because a guard that is red on everything is as useless
    as one that is green on everything.
    """
    print("\n=== STAGE P -- THE PLANTS (every guard watched go red) ===")
    from app.server.reach import steady_state_concentration

    rows: List[Dict] = []

    def record(name: str, fired: bool, detail: str, expect_red: bool) -> None:
        ok = (fired == expect_red)
        rows.append({"guard": name, "fired": fired, "detail": detail,
                     "expected_red": expect_red, "as_expected": ok})
        print(f"  {'RED  ' if fired else 'GREEN'} {'' if ok else '<UNEXPECTED> '}"
              f"{name}: {detail}")

    micro = flat_micro(30, 30, 15.0, ecm=0.3)
    drug = DrugParams(name="plant", diffusivity=300.0, clearance=2e-4, binding=5e-4)
    src = micro.domain.zeros()
    src[::6, ::7] = 1.2
    pde = PDEConfig(dt=10.0, total_time=600.0, convection=False, save_every=10 ** 9)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    a1_tol = PREREG["A1_mass_budget_backward_euler"]["threshold"]

    # --- A1, and its detection floor -------------------------------------- #
    floor: List[Dict] = []
    for eps in (0.0, 1e-13, 3e-12, 1e-11, 1e-10, 1e-9):
        rel = _leaky_march(micro, drug, src, pde, eps)
        floor.append({"eps_per_step": eps, "rel_residual": rel,
                      "fires": rel >= a1_tol})
        record(f"A1 mass budget, planted per-step leak eps={eps:g}",
               rel >= a1_tol, f"rel residual {rel:.4e} vs {a1_tol:.0e}",
               expect_red=(eps >= 1e-11))
    detected = [r["eps_per_step"] for r in floor if r["fires"]]
    missed = [r["eps_per_step"] for r in floor if not r["fires"] and r["eps_per_step"] > 0]
    print(f"    -> A1 detection floor at 60 steps: fires at eps >= "
          f"{min(detected):g}, blind at eps <= {max(missed):g}")

    # --- A3 and B1, on the divergent cross --------------------------------- #
    cm = flat_micro(9, 9, 15.0, ecm=0.0)
    cd = DrugParams(name="tracer", diffusivity=1e-6, clearance=0.0, binding=0.0)
    vel = divergent_cross_velocity(9, 9, 30.0, 4, 4)
    csrc = cm.domain.zeros()
    csrc[4, 4] = 1.0
    speed = float(np.max(np.abs(vel[0]) + np.abs(vel[1])))
    for cfl in (0.4, 0.6, 1.0):
        dt = cfl * 15.0 / speed
        n_steps = 8
        p = PDEConfig(dt=dt, total_time=n_steps * dt, convection=True,
                      save_every=10 ** 9)
        prof = np.zeros(n_steps)
        prof[0] = 1.0
        cl = march_mirror(cm, cd, csrc, p, velocity=vel, source_profile=prof,
                          clamp=True)
        raw = march_mirror(cm, cd, csrc, p, velocity=vel, source_profile=prof,
                           clamp=False)
        created = cl["clamp_mass"] / cl["injected"]
        rel = raw["min_raw"] / raw["max_c"] if raw["max_c"] > 0 else 0.0
        record(f"A3 zero-clamp, CFL={cfl}", created != 0.0,
               f"clamp created {created:.4f} of injected", expect_red=(cfl > 0.5))
        record(f"B1 non-negativity, CFL={cfl}",
               rel < PREREG["B1_nonneg_smooth_flow"]["threshold"],
               f"relative undershoot {rel:+.6f}", expect_red=(cfl > 0.5))

    # --- C1 / C2 / C3 ------------------------------------------------------ #
    base = simulate_drug_transport(micro, drug, src, pde).auc
    up = simulate_drug_transport(micro, drug, 2.0 * src, pde).auc
    down = simulate_drug_transport(micro, drug, 0.5 * src, pde).auc
    record("C1 dose ordering, control (2x source)", not bool(np.all(up >= base)),
           "ordering holds", expect_red=False)
    record("C1 dose ordering, PLANTED decrease (0.5x source)",
           not bool(np.all(down >= base)),
           f"{int((down < base).sum())} grid points fell", expect_red=True)

    short = simulate_drug_transport(
        micro, drug, src,
        PDEConfig(dt=10.0, total_time=600.0, convection=False,
                  save_every=10 ** 9)).auc
    long = simulate_drug_transport(
        micro, drug, src,
        PDEConfig(dt=10.0, total_time=1200.0, convection=False,
                  save_every=10 ** 9)).auc
    record("C2 window, control (longer window)",
           float((long - short)[tissue].min()) < 0.0,
           f"min increment {float((long - short)[tissue].min()):+.4e}",
           expect_red=False)
    record("C2 window, PLANTED shorter window",
           float((short - long)[tissue].min()) < 0.0,
           f"min increment {float((short - long)[tissue].min()):+.4e}",
           expect_red=True)

    long_pde = PDEConfig(dt=20.0, total_time=20000.0, convection=False,
                         save_every=10 ** 9)
    for label, move_k, expect in (("control (D alone)", False, False),
                                  ("PLANTED k moved with D", True, True)):
        totals = []
        for f in (0.25, 1.0, 4.0):
            dr = dataclasses.replace(
                drug, diffusivity=drug.diffusivity * f,
                binding=(drug.binding * f if move_k else drug.binding))
            totals.append(float(
                simulate_drug_transport(micro, dr, src, long_pde).auc[tissue].sum()))
        spread = (max(totals) - min(totals)) / max(totals)
        record(f"C3 total-exposure invariance, {label}", spread >= 1e-6,
               f"spread {spread:.4e}", expect_red=expect)

    # --- D --------------------------------------------------------------- #
    def sha(arr) -> str:
        return hashlib.sha256(
            np.ascontiguousarray(np.asarray(arr, dtype=np.float64)).tobytes()
        ).hexdigest()

    same = simulate_drug_transport(micro, drug, src, pde).auc
    record("D determinism, control (two identical solves)",
           sha(base) != sha(same), "hashes equal", expect_red=False)
    # A GENUINE one-ulp change. The obvious plant, src*(1 + 1e-16), is a NO-OP:
    # for values near 1.2 it rounds back to src exactly, so the guard staying
    # green would prove nothing about the guard. That is a plant APPLIED but not
    # REACHED, and it is recorded here rather than quietly replaced.
    noop = src * (1.0 + 1e-16)
    record("D determinism, PLANTED src*(1+1e-16) -- a NO-OP plant",
           sha(base) != sha(simulate_drug_transport(micro, drug, noop, pde).auc),
           f"the 'perturbed' source is bit-identical to the original: "
           f"{np.array_equal(src, noop)}", expect_red=False)
    ulp = np.nextafter(src, np.inf)
    record("D determinism, PLANTED src at +1 ulp (a real change)",
           sha(base) != sha(simulate_drug_transport(micro, drug, ulp, pde).auc),
           f"source arrays differ: {not np.array_equal(src, ulp)}",
           expect_red=True)

    # --- E ----------------------------------------------------------------- #
    tol = 1e-12
    lam = 10.0
    scaled = simulate_drug_transport(micro, drug, lam * src, pde).auc
    record("E1 source scale, control", _relmax(scaled, lam * base) >= tol,
           f"rel {_relmax(scaled, lam * base):.4e}", expect_red=False)
    record("E1 source scale, PLANTED 1e-9 relative dose error",
           _relmax(scaled, lam * (1 + 1e-9) * base) >= tol,
           f"rel {_relmax(scaled, lam * (1 + 1e-9) * base):.4e}", expect_red=True)

    s = 2.0
    d2 = dataclasses.replace(drug, diffusivity=drug.diffusivity * s,
                             clearance=drug.clearance * s,
                             binding=drug.binding * s)
    good = simulate_drug_transport(
        micro, d2, s * src,
        PDEConfig(dt=10.0 / s, total_time=600.0 / s, convection=False,
                  save_every=10 ** 9)).auc
    record("E3 capacity identity, control", _relmax(s * good, base) >= tol,
           f"rel {_relmax(s * good, base):.4e}", expect_red=False)
    bad = simulate_drug_transport(
        micro, d2, s * src,
        PDEConfig(dt=10.0, total_time=600.0 / s, convection=False,
                  save_every=10 ** 9)).auc
    record("E3 capacity identity, PLANTED dt not rescaled",
           _relmax(s * bad, base) >= tol, f"rel {_relmax(s * bad, base):.4e}",
           expect_red=True)

    c0 = steady_state_concentration(micro, drug, src)
    record("E4 steady invariance, control",
           _relmax(steady_state_concentration(micro, d2, s * src), c0) >= tol,
           f"rel {_relmax(steady_state_concentration(micro, d2, s * src), c0):.4e}",
           expect_red=False)
    d3 = dataclasses.replace(drug, diffusivity=drug.diffusivity * s)  # k left alone
    record("E4 steady invariance, PLANTED k not scaled",
           _relmax(steady_state_concentration(micro, d3, s * src), c0) >= tol,
           f"rel {_relmax(steady_state_concentration(micro, d3, s * src), c0):.4e}",
           expect_red=True)

    # --- the instrument's own guard ---------------------------------------- #
    prod = simulate_drug_transport(micro, drug, src, pde)
    record("mirror bit-identity, control",
           not np.array_equal(prod.auc, march_mirror(micro, drug, src, pde,
                                                     clamp=True)["auc"]),
           "mirror reproduces the solver", expect_red=False)
    off = PDEConfig(dt=float(np.nextafter(10.0, 11.0)), total_time=600.0,
                    convection=False, save_every=10 ** 9)
    record("mirror bit-identity, PLANTED dt + 1 ulp",
           not np.array_equal(prod.auc, march_mirror(micro, drug, src, off,
                                                     clamp=True)["auc"]),
           "mirror run at dt + 1 ulp", expect_red=True)

    surprises = [r["guard"] for r in rows if not r["as_expected"]]
    planted = [r for r in rows if r["expected_red"]]
    reached = [r for r in planted if r["fired"]]
    verdict = (f"{len(reached)}/{len(planted)} plants REACHED their guard"
               + ("" if not surprises else f"; UNEXPECTED: {surprises}"))
    print(f"  {verdict}")
    return {"rows": rows, "a1_floor": floor,
            "n_planted": len(planted), "n_reached": len(reached),
            "unexpected": surprises, "verdict": {"P": verdict}}


# --------------------------------------------------------------------------- #
def main() -> Dict:
    t0 = time.time()
    assert_alpha("START")
    assert_prereg()
    out = {
        "A": stage_mass_balance(),
        "B": stage_nonnegativity(),
        "C": stage_monotonicity(),
        "D": stage_determinism(),
        "E": stage_scale_invariance(),
        "P": stage_plants(),
    }
    print("\n=== VERDICT TABLE ===")
    for stage in "ABCDEP":
        for key, value in out[stage]["verdict"].items():
            print(f"  {key:5s} {value}")
    print(f"\ntotal {time.time() - t0:.1f}s")
    assert_alpha("END")
    return out


if __name__ == "__main__":
    main()
