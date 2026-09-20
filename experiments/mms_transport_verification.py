"""METHOD OF MANUFACTURED SOLUTIONS on the transport solver.

Why this file exists
--------------------
This repository has grid-convergence studies (``experiments/grid_convergence_*.py``,
``tmesim/validation/numerical_rigor.py``) but they are all Richardson
extrapolation on *derived quantities of interest* — penetration index, kill
fraction — against no known answer. A Richardson slope on a QoI tells you the
answer is still moving with the mesh; it cannot tell you whether the PDE solver
solves the PDE it claims to. Nothing in the tree verified that. Searched for
before writing this: ``manufactured``, ``MMS``, ``exact solution``,
``analytic_solution``, ``observed order`` — the only hits are Richardson-on-QoI
and one comment in ``experiments/boundary_layer_resolution.py``.

MMS is the right tool here precisely because it needs no tissue, no patient, no
anchor and no calibration argument: pick an exact solution, substitute it into
the governing equation to derive the source that makes it exact, run the solver
with that source, and watch the error norm fall as dx does. The observed order
either matches the scheme's formal order or the implementation has a bug. That
verdict is independent of every accuracy dispute in ``results/``.

The equation as IMPLEMENTED (read from ``tmesim/physics/diffusion.py``)
-----------------------------------------------------------------------
    dC/dt = div(D_eff grad C) - div(v C) + S - k C
    D_eff = drug.diffusivity * (1 - 0.6 * ecm), floored at max(1e-6, 1e-3*D)
    k     = clearance + binding (+ exchange_k when supplied)
    BC    : zero TOTAL flux across every domain-edge face and every
            tissue/non-tissue face (the face coefficient is set to zero)
    IC    : C == 0, not settable through the public API
    Scheme: IMEX -- backward Euler on diffusion+reaction (implicit, single LU,
            reused), first-order donor-cell upwind on advection (explicit),
            source evaluated at the OLD time level; post-step clamp
            C <- max(C, 0) and C <- 0 outside tissue.
    Grid  : cell-centred, x_j = (j + 1/2) dx, extent nx*dx (Domain.coords_um),
            which is exactly the convention the flux-zeroing BC implies.

Formal orders that scheme should have: 2 in space for diffusion (5-point finite
volume with harmonic face conductance and the no-flux condition imposed AT the
face), 1 in space for advection (donor cell), 1 in time (backward Euler with a
lagged source).

Design
------
Four of the five cases are STEADY manufactured solutions. That is deliberate:
the fixed point of the solver's own iteration is

    C* = (K - L)^{-1} S

which is the *discrete steady problem* and is exactly independent of dt, so a
converged steady run isolates the SPATIAL discretisation with no temporal
contamination at all. The fifth case refines dt instead and uses a spatially
uniform exact solution (annihilated exactly by L, whose rows sum to zero), which
isolates the TIME integrator and the source convention.

The manufactured solutions are chosen with du/dn = 0 on every no-flux face, so
they are compatible with the only boundary condition the solver can impose, and
strictly positive, so the non-negativity clamp provably never binds (see V1).

Usage
-----
    .venv/bin/python experiments/mms_transport_verification.py           # all
    .venv/bin/python experiments/mms_transport_verification.py --case E3
    .venv/bin/python experiments/mms_transport_verification.py --case E1 --quick
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from tmesim.core import (
    CellPopulation,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    TissueFields,
)
from tmesim.physics.diffusion import ALPHA_TORTUOSITY, simulate_drug_transport

# --------------------------------------------------------------------------
# Pre-registered thresholds. Fixed before any error was measured; see
# results/mms_transport_verification.md for the frozen pre-registration.
# --------------------------------------------------------------------------
ORDER2_BAND = (1.8, 2.2)     # second-order cases E1, E2, E5
ORDER1_BAND = (0.8, 1.2)     # first-order cases E3 (upwind), E4 (backward Euler)
RESIDUAL_FACTOR = 1e-3       # V2: steady residual must be this much below the error
ALGEBRA_TOL = 1e-6           # V3: closed-form S vs nested finite differences

L_DOMAIN = 200.0             # um, fixed physical extent for every steady case
DX_LADDER = [4.0, 2.0, 1.0, 0.5]
DT_LADDER = [4.0, 2.0, 1.0, 0.5]
QUICK_DX_LADDER = [4.0, 2.0, 1.0]
QUICK_DT_LADDER = [4.0, 2.0, 1.0]


# --------------------------------------------------------------------------
# Manufactured solution, coefficients and the derived source (closed form)
# --------------------------------------------------------------------------
@dataclass
class Manufactured:
    """u, D_eff, v and the derived source S for one steady MMS case.

    ``x0, y0`` is the lower-left no-flux corner and ``lx, ly`` the extent of the
    no-flux box (the whole domain, or the tissue sub-block for E5).
    """

    amp_a: float = 2.0            # A in u = A + B cos cos  (A - B > 0 keeps u > 0)
    amp_b: float = 1.0            # B
    modes: int = 2                # half-waves per box side; du/dn = 0 at both ends
    x0: float = 0.0
    y0: float = 0.0
    lx: float = L_DOMAIN
    ly: float = L_DOMAIN

    d0: float = 300.0             # drug.diffusivity, um^2/s
    ecm_mean: float = 0.0         # ecm = ecm_mean + ecm_amp*sin(2pi X/lx)sin(2pi Y/ly)
    ecm_amp: float = 0.0
    reaction: float = 1.0         # k = clearance + binding, 1/s
    vel: float = 0.0              # V in v = (V sin(pi X/lx), V sin(pi Y/ly)), um/s

    # ---- exact solution -------------------------------------------------
    def kx(self) -> float:
        return math.pi * self.modes / self.lx

    def ky(self) -> float:
        return math.pi * self.modes / self.ly

    def u(self, x, y):
        return self.amp_a + self.amp_b * np.cos(self.kx() * (x - self.x0)) * np.cos(
            self.ky() * (y - self.y0)
        )

    def _grad_u(self, x, y):
        kx, ky = self.kx(), self.ky()
        sx = np.sin(kx * (x - self.x0))
        cx = np.cos(kx * (x - self.x0))
        sy = np.sin(ky * (y - self.y0))
        cy = np.cos(ky * (y - self.y0))
        return (-self.amp_b * kx * sx * cy, -self.amp_b * ky * cx * sy)

    def _lap_u(self, x, y):
        kx, ky = self.kx(), self.ky()
        return -self.amp_b * (kx * kx + ky * ky) * np.cos(kx * (x - self.x0)) * np.cos(
            ky * (y - self.y0)
        )

    # ---- coefficients ---------------------------------------------------
    def ecm(self, x, y):
        if self.ecm_amp == 0.0:
            return np.full_like(np.asarray(x, dtype=np.float64), self.ecm_mean)
        a = 2.0 * math.pi / self.lx
        b = 2.0 * math.pi / self.ly
        return self.ecm_mean + self.ecm_amp * np.sin(a * (x - self.x0)) * np.sin(
            b * (y - self.y0)
        )

    def d_eff(self, x, y):
        """Must reproduce ``_effective_diffusivity`` exactly (it does: the ecm
        field stays inside [0,1] and well above the floor, so neither clip nor
        floor is active — asserted in ``_check_coefficients``)."""
        return self.d0 * (1.0 - ALPHA_TORTUOSITY * self.ecm(x, y))

    def _grad_d(self, x, y):
        if self.ecm_amp == 0.0:
            z = np.zeros_like(np.asarray(x, dtype=np.float64))
            return (z, z)
        a = 2.0 * math.pi / self.lx
        b = 2.0 * math.pi / self.ly
        c = -self.d0 * ALPHA_TORTUOSITY * self.ecm_amp
        return (
            c * a * np.cos(a * (x - self.x0)) * np.sin(b * (y - self.y0)),
            c * b * np.sin(a * (x - self.x0)) * np.cos(b * (y - self.y0)),
        )

    def velocity(self, x, y):
        """v with v.n = 0 on every face of the no-flux box, so the solver's
        zeroing of the boundary advective flux is the CORRECT condition for this
        manufactured solution rather than an inconsistency."""
        if self.vel == 0.0:
            z = np.zeros_like(np.asarray(x, dtype=np.float64))
            return (z, z)
        return (
            self.vel * np.sin(math.pi * (x - self.x0) / self.lx),
            self.vel * np.sin(math.pi * (y - self.y0) / self.ly),
        )

    def _div_v(self, x, y):
        if self.vel == 0.0:
            z = np.zeros_like(np.asarray(x, dtype=np.float64))
            return z
        return self.vel * (math.pi / self.lx) * np.cos(
            math.pi * (x - self.x0) / self.lx
        ) + self.vel * (math.pi / self.ly) * np.cos(math.pi * (y - self.y0) / self.ly)

    # ---- derived source: S = k u - div(D grad u) + div(v u) --------------
    def div_d_grad_u(self, x, y):
        ux, uy = self._grad_u(x, y)
        dx_, dy_ = self._grad_d(x, y)
        return self.d_eff(x, y) * self._lap_u(x, y) + dx_ * ux + dy_ * uy

    def div_v_u(self, x, y):
        if self.vel == 0.0:
            return np.zeros_like(np.asarray(x, dtype=np.float64))
        ux, uy = self._grad_u(x, y)
        vx, vy = self.velocity(x, y)
        return vx * ux + vy * uy + self.u(x, y) * self._div_v(x, y)

    def source(self, x, y):
        return self.reaction * self.u(x, y) - self.div_d_grad_u(x, y) + self.div_v_u(x, y)


# --------------------------------------------------------------------------
# V3 -- guard the hand algebra with nested finite differences
# --------------------------------------------------------------------------
def check_algebra(m: Manufactured, n_probe: int = 7, h: float = 0.05) -> Tuple[float, float]:
    """Compare the closed-form ``div(D grad u)`` and ``div(v u)`` against nested
    finite differences of the SAME functions. Returns the two relative errors.

    This does not test the solver; it tests that the manufactured source really
    is the source that makes ``u`` an exact solution. If it were not, ``u`` would
    not be a solution at all and every order below would be meaningless.
    """
    rng = np.random.default_rng(0)
    xs = m.x0 + m.lx * (0.15 + 0.7 * rng.random(n_probe))
    ys = m.y0 + m.ly * (0.15 + 0.7 * rng.random(n_probe))

    def dgu(x, y):
        fe = m.d_eff(x + h / 2, y) * (m.u(x + h, y) - m.u(x, y)) / h
        fw = m.d_eff(x - h / 2, y) * (m.u(x, y) - m.u(x - h, y)) / h
        fn = m.d_eff(x, y + h / 2) * (m.u(x, y + h) - m.u(x, y)) / h
        fs = m.d_eff(x, y - h / 2) * (m.u(x, y) - m.u(x, y - h)) / h
        return (fe - fw) / h + (fn - fs) / h

    def dvu(x, y):
        vxp, _ = m.velocity(x + h, y)
        vxm, _ = m.velocity(x - h, y)
        _, vyp = m.velocity(x, y + h)
        _, vym = m.velocity(x, y - h)
        return (vxp * m.u(x + h, y) - vxm * m.u(x - h, y)) / (2 * h) + (
            vyp * m.u(x, y + h) - vym * m.u(x, y - h)
        ) / (2 * h)

    a_exact = m.div_d_grad_u(xs, ys)
    a_fd = dgu(xs, ys)
    err_a = float(np.max(np.abs(a_exact - a_fd)) / max(np.max(np.abs(a_exact)), 1e-300))

    if m.vel == 0.0:
        err_b = 0.0
    else:
        b_exact = m.div_v_u(xs, ys)
        b_fd = dvu(xs, ys)
        err_b = float(np.max(np.abs(b_exact - b_fd)) / max(np.max(np.abs(b_exact)), 1e-300))
    return err_a, err_b


# --------------------------------------------------------------------------
# Grid plumbing
# --------------------------------------------------------------------------
def build_micro(
    dx: float, m: Manufactured, masked: bool
) -> Tuple[Microenvironment, np.ndarray, np.ndarray]:
    n = int(round(L_DOMAIN / dx))
    assert abs(n * dx - L_DOMAIN) < 1e-12, (dx, n)

    xs = (np.arange(n) + 0.5) * dx
    ys = (np.arange(n) + 0.5) * dx
    X, Y = np.meshgrid(xs, ys)

    if masked:
        tissue = (
            (X > m.x0) & (X < m.x0 + m.lx) & (Y > m.y0) & (Y < m.y0 + m.ly)
        )
        # The sub-block must land exactly on cell faces on every rung of the
        # ladder, otherwise the "no-flux box" moves with dx and the manufactured
        # solution stops satisfying the boundary condition.
        assert abs(m.x0 / dx - round(m.x0 / dx)) < 1e-12
        assert abs((m.x0 + m.lx) / dx - round((m.x0 + m.lx) / dx)) < 1e-12
    else:
        tissue = np.ones((n, n), dtype=bool)

    domain = Domain(nx=n, ny=n, dx=dx, tissue_mask=tissue)
    fields = TissueFields(
        ecm_density=np.clip(m.ecm(X, Y), 0.0, 1.0),
        vascular_density=np.zeros((n, n)),
    )
    cells = CellPopulation(
        x=np.array([L_DOMAIN / 2]),
        y=np.array([L_DOMAIN / 2]),
        cell_type=np.array([0]),
        state=np.array([0]),
        sensitivity=np.array([0.5]),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells), X, Y


def _check_coefficients(m: Manufactured, micro: Microenvironment, X, Y) -> None:
    """The solver's own D_eff must equal the analytic D_eff to machine precision
    (i.e. neither the [0,1] ecm clip nor the positive floor is active)."""
    from tmesim.physics.diffusion import _effective_diffusivity

    drug = DrugParams(name="mms", diffusivity=m.d0, clearance=m.reaction, binding=0.0)
    d_solver = _effective_diffusivity(drug, micro.fields.ecm_density)
    d_exact = m.d_eff(X, Y)
    rel = float(np.max(np.abs(d_solver - d_exact)) / np.max(np.abs(d_exact)))
    assert rel < 1e-14, f"solver D_eff != analytic D_eff (rel {rel:.3e})"


def boundary_ring(tissue: np.ndarray) -> np.ndarray:
    """Tissue cells with at least one neighbour that is off-domain or non-tissue.

    These are exactly the cells whose stencil has a zeroed face, i.e. the cells
    the no-flux condition actually touches. Everything else is interior.
    """
    t = tissue
    pad = np.zeros((t.shape[0] + 2, t.shape[1] + 2), dtype=bool)
    pad[1:-1, 1:-1] = t
    nb = (
        pad[:-2, 1:-1] & pad[2:, 1:-1] & pad[1:-1, :-2] & pad[1:-1, 2:]
    )  # all four neighbours are tissue AND on-domain
    return t & ~nb


# --------------------------------------------------------------------------
# One rung of a steady ladder
# --------------------------------------------------------------------------
@dataclass
class Level:
    dx: float
    n: int
    l2_all: float
    linf_all: float
    l2_bnd: float
    linf_bnd: float
    l2_int: float
    linf_int: float
    residual: float
    smin: float
    cmin: float
    cfl: float
    seconds: float


def run_steady_level(
    dx: float, m: Manufactured, dt: float, total_time: float, masked: bool
) -> Level:
    t0 = time.time()
    micro, X, Y = build_micro(dx, m, masked)
    _check_coefficients(m, micro, X, Y)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)

    drug = DrugParams(name="mms", diffusivity=m.d0, clearance=m.reaction, binding=0.0)
    src = m.source(X, Y)
    smin = float(src[tissue].min())

    use_conv = m.vel != 0.0
    if use_conv:
        vx, vy = m.velocity(X, Y)
        velocity = np.stack([vx, vy])
        cfl = float(np.max(np.abs(vx) + np.abs(vy)) * dt / dx)
    else:
        velocity = None
        cfl = 0.0

    pde = PDEConfig(
        dt=dt, total_time=total_time, convection=use_conv, save_every=10 ** 9
    )
    res = simulate_drug_transport(
        micro, drug, src, pde, velocity=velocity, capture_frames=3
    )
    c = res.final_concentration
    # frames are at t = 0, T/2, T -> V2 steady-state residual over the last half
    residual = float(np.max(np.abs(res.frames[2] - res.frames[1])[tissue]))

    err = np.abs(c - m.u(X, Y))
    bnd = boundary_ring(tissue)
    inte = tissue & ~bnd

    def norms(mask):
        e = err[mask]
        return float(np.sqrt(np.mean(e ** 2))), float(np.max(e))

    l2a, lia = norms(tissue)
    l2b, lib = norms(bnd)
    l2i, lii = norms(inte)
    return Level(
        dx=dx,
        n=micro.domain.nx,
        l2_all=l2a,
        linf_all=lia,
        l2_bnd=l2b,
        linf_bnd=lib,
        l2_int=l2i,
        linf_int=lii,
        residual=residual,
        smin=smin,
        cmin=float(c[tissue].min()),
        cfl=cfl,
        seconds=time.time() - t0,
    )


# --------------------------------------------------------------------------
# Temporal ladder (E4): spatially uniform exact solution
# --------------------------------------------------------------------------
@dataclass
class TLevel:
    dt: float
    err: float
    smin: float
    seconds: float


def run_temporal_level(dt: float, k: float, cmax: float, tau: float, total_time: float) -> TLevel:
    t0 = time.time()
    n = 16
    domain = Domain(nx=n, ny=n, dx=10.0)
    fields = TissueFields(
        ecm_density=np.zeros((n, n)), vascular_density=np.zeros((n, n))
    )
    cells = CellPopulation(
        x=np.array([80.0]), y=np.array([80.0]), cell_type=np.array([0]),
        state=np.array([0]), sensitivity=np.array([0.5]),
    )
    micro = Microenvironment(domain=domain, fields=fields, cells=cells)
    drug = DrugParams(name="mms-t", diffusivity=300.0, clearance=k, binding=0.0)

    n_steps = int(round(total_time / dt))
    assert abs(n_steps * dt - total_time) < 1e-12
    t_old = np.arange(n_steps) * dt            # source is read at the OLD level
    # u(t) = cmax (1 - exp(-t/tau));  S = u' + k u
    prof = (cmax / tau) * np.exp(-t_old / tau) + k * cmax * (1.0 - np.exp(-t_old / tau))

    pde = PDEConfig(dt=dt, total_time=total_time, convection=False, save_every=10 ** 9)
    res = simulate_drug_transport(
        micro, drug, np.ones((n, n)), pde, source_profile=prof
    )
    exact = cmax * (1.0 - math.exp(-total_time / tau))
    err = float(np.max(np.abs(res.final_concentration - exact)))
    return TLevel(dt=dt, err=err, smin=float(prof.min()), seconds=time.time() - t0)


# --------------------------------------------------------------------------
# Order bookkeeping / reporting
# --------------------------------------------------------------------------
def orders(errs: List[float], ratio: float = 2.0) -> List[float]:
    out = []
    for a, b in zip(errs[:-1], errs[1:]):
        if b <= 0.0 or a <= 0.0:
            out.append(float("nan"))
        else:
            out.append(math.log(a / b) / math.log(ratio))
    return out


def verdict(ps: List[float], band: Tuple[float, float]) -> bool:
    return all(band[0] <= p <= band[1] for p in ps if not math.isnan(p)) and not any(
        math.isnan(p) for p in ps
    )


def report_steady(name: str, desc: str, levels: List[Level], band, algebra, check_bnd=True):
    lo, hi = band
    print(f"\n{'=' * 78}\n{name}: {desc}\n{'=' * 78}")
    print(f"  V3 algebra check (closed-form S vs nested FD): "
          f"div(D grad u) rel {algebra[0]:.3e}, div(vu) rel {algebra[1]:.3e}  "
          f"[tol {ALGEBRA_TOL:.0e}] -> {'OK' if max(algebra) < ALGEBRA_TOL else 'VOID'}")
    print(f"\n  {'dx':>6} {'n':>5} {'L2(all)':>12} {'Linf(all)':>12} "
          f"{'Linf(bnd)':>12} {'Linf(int)':>12} {'resid/err':>10} {'min S':>9} "
          f"{'min C':>9} {'CFL':>6} {'s':>6}")
    for lv in levels:
        rr = lv.residual / lv.linf_all if lv.linf_all > 0 else float("inf")
        print(f"  {lv.dx:>6.3g} {lv.n:>5d} {lv.l2_all:>12.5e} {lv.linf_all:>12.5e} "
              f"{lv.linf_bnd:>12.5e} {lv.linf_int:>12.5e} {rr:>10.2e} "
              f"{lv.smin:>9.4g} {lv.cmin:>9.4g} {lv.cfl:>6.3f} {lv.seconds:>6.1f}")

    series = {
        "L2   (global)": [lv.l2_all for lv in levels],
        "Linf (global)": [lv.linf_all for lv in levels],
    }
    if check_bnd:
        series["Linf (boundary)"] = [lv.linf_bnd for lv in levels]
        series["Linf (interior)"] = [lv.linf_int for lv in levels]
        series["L2   (boundary)"] = [lv.l2_bnd for lv in levels]

    print(f"\n  observed order per refinement step (r = 2), pre-registered band "
          f"[{lo}, {hi}]:")
    allpass = True
    for label, errs in series.items():
        ps = orders(errs)
        ok = verdict(ps, band)
        # L2(boundary) is reported for information only; the pre-registration
        # names L2(global), Linf(global), Linf(boundary), Linf(interior).
        graded = label != "L2   (boundary)"
        if graded:
            allpass = allpass and ok
        tag = ("PASS" if ok else "FAIL") if graded else ("ok" if ok else "--")
        steps = "  ".join(
            f"{lv0.dx:g}->{lv1.dx:g}: {p:6.3f}"
            for lv0, lv1, p in zip(levels[:-1], levels[1:], ps)
        )
        print(f"    {label:<16} {steps}   [{tag}]")

    resid_ok = all(
        lv.residual < RESIDUAL_FACTOR * lv.linf_all for lv in levels
    )
    s_ok = all(lv.smin > 0.0 for lv in levels)
    c_ok = all(lv.cmin > 0.0 for lv in levels)
    cfl_ok = all(lv.cfl <= 1.0 for lv in levels)
    alg_ok = max(algebra) < ALGEBRA_TOL
    print(f"\n  guards: V1 min(S)>0 {s_ok} (min C>0 {c_ok})   "
          f"V2 residual<{RESIDUAL_FACTOR:g}*Linf {resid_ok}   "
          f"V3 algebra {alg_ok}   V4 CFL<=1 {cfl_ok}")
    void = not (resid_ok and s_ok and alg_ok and cfl_ok)
    print(f"  {name} VERDICT: {'VOID (guard failed)' if void else ('PASS' if allpass else 'FAILED')}")
    return (not void) and allpass, void


# --------------------------------------------------------------------------
# Cases
# --------------------------------------------------------------------------
def case_E1(ladder):
    m = Manufactured(d0=300.0, ecm_mean=0.0, ecm_amp=0.0, reaction=1.0, vel=0.0)
    alg = check_algebra(m)
    levels = [run_steady_level(dx, m, dt=1.0, total_time=100.0, masked=False)
              for dx in ladder]
    return report_steady(
        "E1", "pure diffusion, CONSTANT D, domain-edge no-flux (expect order 2)",
        levels, ORDER2_BAND, alg)


def case_E2(ladder):
    m = Manufactured(d0=300.0, ecm_mean=0.4, ecm_amp=0.3, reaction=1.0, vel=0.0)
    alg = check_algebra(m)
    levels = [run_steady_level(dx, m, dt=1.0, total_time=100.0, masked=False)
              for dx in ladder]
    return report_steady(
        "E2", "pure diffusion, VARIABLE D (harmonic face averaging) (expect order 2)",
        levels, ORDER2_BAND, alg)


def case_E3(ladder):
    m = Manufactured(d0=2.0, ecm_mean=0.0, ecm_amp=0.0, reaction=0.1, vel=0.2)
    alg = check_algebra(m)
    levels = [run_steady_level(dx, m, dt=1.0, total_time=600.0, masked=False)
              for dx in ladder]
    return report_steady(
        "E3", "CONVECTION + diffusion, donor-cell upwind (expect order 1)",
        levels, ORDER1_BAND, alg, check_bnd=False)


def case_E5(ladder):
    m = Manufactured(
        d0=300.0, ecm_mean=0.0, ecm_amp=0.0, reaction=2.0, vel=0.0,
        x0=40.0, y0=40.0, lx=120.0, ly=120.0,
    )
    alg = check_algebra(m)
    levels = [run_steady_level(dx, m, dt=1.0, total_time=100.0, masked=True)
              for dx in ladder]
    return report_steady(
        "E5", "TISSUE-MASK internal no-flux boundary, constant D (expect order 2)",
        levels, ORDER2_BAND, alg)


def case_E4(ladder):
    k, cmax, tau, T = 0.01, 10.0, 50.0, 320.0
    lvls = [run_temporal_level(dt, k, cmax, tau, T) for dt in ladder]
    lo, hi = ORDER1_BAND
    print(f"\n{'=' * 78}\nE4: TEMPORAL order, spatially uniform exact solution "
          f"(expect order 1)\n{'=' * 78}")
    print(f"  u(t) = {cmax} (1 - exp(-t/{tau:g})), S(t) = u'(t) + k u(t), "
          f"k = {k}, T = {T:g} s")
    print(f"\n  {'dt':>6} {'|C(T)-u(T)|':>14} {'min S':>10} {'s':>6}")
    for lv in lvls:
        print(f"  {lv.dt:>6.3g} {lv.err:>14.6e} {lv.smin:>10.4g} {lv.seconds:>6.1f}")
    ps = orders([lv.err for lv in lvls])
    ok = verdict(ps, ORDER1_BAND)
    steps = "  ".join(
        f"{a.dt:g}->{b.dt:g}: {p:6.3f}" for a, b, p in zip(lvls[:-1], lvls[1:], ps)
    )
    print(f"\n  observed order per step (r = 2), pre-registered band [{lo}, {hi}]:")
    print(f"    |C(T)-u(T)|      {steps}   [{'PASS' if ok else 'FAIL'}]")
    s_ok = all(lv.smin >= 0.0 for lv in lvls)
    print(f"\n  guards: V1 min(S)>=0 {s_ok}")
    print(f"  E4 VERDICT: {'VOID (guard failed)' if not s_ok else ('PASS' if ok else 'FAILED')}")
    return ok and s_ok, not s_ok


def case_E4b(ladder):
    """POST-HOC diagnosis of the E4 failure. Not part of the pre-registration.

    The pre-registered E4 used k = 0.01 with tau = 50 s. For this scheme the
    leading local truncation error is (dt/2)(S' + k u'), and with
    u = Cmax(1-exp(-t/tau)), S = u' + k u it equals
        (dt/2) (Cmax/tau) exp(-t/tau) (2k - 1/tau),
    which vanishes identically when 2k = 1/tau -- exactly the pre-registered
    pair. So E4's order-2 result is a degeneracy in the MANUFACTURED SOLUTION,
    not evidence that backward Euler is second order. Predictions registered
    before this ran: order in [0.8, 1.2] whenever 2k != 1/tau, order ~2 at
    k = 0.01, and the dt = 4 error roughly proportional to |2k - 1/tau|.
    """
    cmax, tau, T = 10.0, 50.0, 320.0
    print(f"\n{'=' * 78}\nE4b (POST-HOC, not pre-registered): temporal order vs the "
          f"degeneracy 2k = 1/tau\n{'=' * 78}")
    print(f"  leading LTE coefficient ~ (Cmax/(2 tau)) exp(-t/tau) (2k - 1/tau); "
          f"1/tau = {1.0 / tau:g}")
    print(f"\n  {'k':>8} {'2k-1/tau':>10} " +
          "".join(f"{('err dt=' + format(d, 'g')):>14}" for d in ladder) +
          "   observed order per step")
    rows = []
    for k in (0.005, 0.01, 0.02, 0.04):
        lvls = [run_temporal_level(dt, k, cmax, tau, T) for dt in ladder]
        errs = [lv.err for lv in lvls]
        ps = orders(errs)
        rows.append((k, errs, ps))
        steps = "  ".join(f"{p:5.3f}" for p in ps)
        print(f"  {k:>8.4g} {2 * k - 1 / tau:>10.4g} " +
              "".join(f"{e:>14.6e}" for e in errs) + f"   {steps}")
    ok = all(
        verdict(ps, ORDER1_BAND) for k, _, ps in rows if abs(2 * k - 1 / tau) > 1e-12
    )
    degen_ok = all(
        max(ps) > 1.8 for k, _, ps in rows if abs(2 * k - 1 / tau) <= 1e-12
    )

    # Prediction 2 as registered: err(dt=4) proportional to |2k - 1/tau|.
    ref = [r for r in rows if abs(2 * r[0] - 1 / tau - (-0.01)) < 1e-12][0]
    hit = [r for r in rows if abs(2 * r[0] - 1 / tau - 0.02) < 1e-12][0]
    ratio_meas = hit[1][0] / ref[1][0]
    print(f"\n  prediction 1 (order in [0.8,1.2] whenever 2k != 1/tau): {ok}")
    print(f"  prediction 3 (order ~2 at the degenerate k = 1/(2 tau)):  {degen_ok}")
    print(f"  prediction 2 (err(dt=4) ~ |2k - 1/tau|): predicted ratio 2.00 for "
          f"k=0.02 vs k=0.005, MEASURED {ratio_meas:.3f} -> "
          f"{'holds' if 1.6 < ratio_meas < 2.5 else 'FAILED'}")

    # Post-hoc: the registered prediction 2 ignored how a per-step error is
    # PROPAGATED to t = T by the reaction, exp(-k(T-s)). Including it,
    #   E(T) = (dt/2)(Cmax/tau)(2k - 1/tau) exp(-kT) (exp((k-1/tau)T)-1)/(k-1/tau)
    # This closed form is fitted to nothing; it is stated and then compared.
    print(f"\n  POST-HOC closed form including the exp(-k(T-s)) propagation "
          f"(no fitted constant):")
    print(f"    {'k':>8} {'predicted':>14} {'measured':>14} {'ratio':>8}   (dt = "
          f"{ladder[0]:g} s)")
    worst = 0.0
    for k, errs, _ps in rows:
        d = k - 1.0 / tau
        integ = T if abs(d) < 1e-15 else (math.exp(d * T) - 1.0) / d
        pred = abs(0.5 * ladder[0] * (cmax / tau) * (2 * k - 1 / tau)
                   * math.exp(-k * T) * integ)
        r = pred / errs[0] if errs[0] > 0 else float("nan")
        worst = max(worst, abs(math.log(r))) if r > 0 else worst
        print(f"    {k:>8.4g} {pred:>14.6e} {errs[0]:>14.6e} {r:>8.3f}")
    print(f"    worst |log ratio| = {worst:.3f} "
          f"({100 * (math.exp(worst) - 1):.1f}% off)")
    return ok and degen_ok, False


CASES: Dict[str, Callable] = {
    "E1": case_E1, "E2": case_E2, "E3": case_E3, "E4": case_E4,
    "E4b": case_E4b, "E5": case_E5,
}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", action="append", choices=sorted(CASES),
                    help="run only this case (repeatable); default: all")
    ap.add_argument("--quick", action="store_true",
                    help="drop the finest rung of the ladder (plant tests)")
    args = ap.parse_args(argv)

    from tmesim.validation.calibration import calibrate_alpha
    alpha = calibrate_alpha()
    print(f"calibrate_alpha() = {alpha!r}")
    assert alpha == 22.223477711172563, alpha

    dxs = QUICK_DX_LADDER if args.quick else DX_LADDER
    dts = QUICK_DT_LADDER if args.quick else DT_LADDER
    names = args.case or sorted(CASES)

    results = {}
    for nm in names:
        ladder = dts if nm.startswith("E4") else dxs
        results[nm] = CASES[nm](ladder)

    print(f"\n{'=' * 78}\nSUMMARY\n{'=' * 78}")
    for nm in names:
        ok, void = results[nm]
        print(f"  {nm}: {'VOID' if void else ('PASS' if ok else 'FAILED')}")
    overall = all(ok for ok, void in results.values())
    print(f"\n  OVERALL: {'PASS' if overall else 'FAILED'}")
    print(f"calibrate_alpha() = {calibrate_alpha()!r}")
    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
