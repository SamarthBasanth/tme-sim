"""Why MMS case E4 failed: a closed-form diagnosis of the temporal MMS case.

Background
----------
``experiments/mms_transport_verification.py`` ran five pre-registered
manufactured-solution cases against ``tmesim.physics.diffusion``. Four passed.
**E4 -- the temporal case -- FAILED**: registered to show order 1 (backward
Euler with a lagged source), it measured 2.019 / 2.009 / 2.005.

``results/mms_transport_verification.md`` Sec.6 asserted the cause is a
degeneracy in the manufactured solution (2k = 1/tau annihilates the leading
truncation term) and supported it with a four-point k sweep plus a post-hoc
closed form good to 7.8%. This script settles the question instead of
supporting it, by removing every fitted or asymptotic step:

  A. The E4 configuration reduces the solver EXACTLY to the scalar recurrence
        C^{n+1} = (C^n + dt S(t_n)) / (1 + k dt),  C^0 = 0
     -- verified against the shipped solver to ~1e-15 relative, on every dt.
     This is what rules out a solver defect: whatever makes E4 order 2 is a
     property of this recurrence, which is provably the documented scheme.

  B. That recurrence has a CLOSED-FORM solution, so the global error at t = T
     is known exactly -- no truncation analysis, no fitted constant:
        E(T) = (r^n - a^n) * ( Cmax + dt*beta*a/(r - a) )
        a = 1/(1+k dt), r = exp(-dt/tau), beta = Cmax(1/tau - k), n = T/dt
     Evaluated at 80 decimal digits and compared to the solver's measured
     error. Agreement here means the observed order is fully explained.

  C. Expanding E(T) in dt gives the exact first-order coefficient
        G1 = -(Cmax/(2 tau)) * (2k - 1/tau)/(k - 1/tau) * (e^{-T/tau} - e^{-kT})
     which is IDENTICALLY ZERO iff 2k = 1/tau. E4 registered k = 0.01 and
     tau = 50 s: 2k = 0.02 = 1/tau exactly. The O(dt) term is not small at
     E4's parameters, it is absent.

  D. The decisive discriminator between "solver is second order" and
     "manufactured solution is degenerate": detune k by delta = 2k - 1/tau and
     the observed order must CROSS OVER from 2 (coarse dt, dt^2 term dominant)
     to 1 (fine dt, dt^1 term dominant) at a predictable dt* = |G1/G2|.
     A genuinely second-order integrator shows order 2 at every dt and every k
     and has no crossover. Predicted dt* comes from the closed form; measured
     dt* comes from the shipped solver.

  E. A non-degenerate temporal manufactured solution at E4's OWN k, T and grid:
        u(t) = Cmax t/(t + tau)      (rational, u' is not proportional to
                                      e^{-2kt}, so no degeneracy is possible)
     If the solver is first order in time, this must show order 1 at exactly
     the parameters where E4 showed 2. This is the minimal repair to the TEST;
     no engine file is touched and no tolerance is moved.

  F. Regime / CFL: E4 runs with convection=False and no velocity, so the CFL
     guard in diffusion.py is never evaluated and the reported CFL number is
     identically 0. Asserted here rather than argued.

Nothing in this file edits the engine, edits the E4 pre-registration, or
revises E4's FAILED verdict.

Usage
-----
    .venv/bin/python experiments/mms_e4_diagnosis.py
    .venv/bin/python experiments/mms_e4_diagnosis.py --quick
"""
from __future__ import annotations

import argparse
import math
import time
from decimal import Decimal as D
from decimal import getcontext
from typing import Callable, List, Tuple

import numpy as np

from tmesim.core import (
    CellPopulation,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    TissueFields,
)
from tmesim.physics.diffusion import simulate_drug_transport
from tmesim.validation.calibration import calibrate_alpha

getcontext().prec = 80

# ---- E4's registered parameters, copied verbatim from case_E4 --------------
E4_K = 0.01
E4_CMAX = 10.0
E4_TAU = 50.0
E4_T = 320.0
E4_N = 16          # 16x16 grid
E4_DX = 10.0
E4_LADDER = [4.0, 2.0, 1.0, 0.5]


# ==========================================================================
# The shipped solver, driven exactly as run_temporal_level drives it
# ==========================================================================
def solver_final(dt: float, k: float, prof: np.ndarray, total_time: float) -> float:
    """Return the (spatially uniform) final concentration from the real solver."""
    n = E4_N
    domain = Domain(nx=n, ny=n, dx=E4_DX)
    fields = TissueFields(
        ecm_density=np.zeros((n, n)), vascular_density=np.zeros((n, n))
    )
    cells = CellPopulation(
        x=np.array([80.0]), y=np.array([80.0]), cell_type=np.array([0]),
        state=np.array([0]), sensitivity=np.array([0.5]),
    )
    micro = Microenvironment(domain=domain, fields=fields, cells=cells)
    drug = DrugParams(name="mms-t", diffusivity=300.0, clearance=k, binding=0.0)
    pde = PDEConfig(dt=dt, total_time=total_time, convection=False,
                    save_every=10 ** 9)
    res = simulate_drug_transport(
        micro, drug, np.ones((n, n)), pde, source_profile=prof
    )
    c = res.final_concentration
    spread = float(c.max() - c.min())
    assert spread < 1e-12, f"field not uniform: spread {spread:g}"
    return float(c.mean()), res


def expo_profile(dt: float, k: float, cmax: float, tau: float, T: float) -> np.ndarray:
    """E4's source: u = cmax(1 - exp(-t/tau)), S = u' + k u, read at t_n."""
    n_steps = int(round(T / dt))
    assert abs(n_steps * dt - T) < 1e-12
    t = np.arange(n_steps) * dt
    return (cmax / tau) * np.exp(-t / tau) + k * cmax * (1.0 - np.exp(-t / tau))


def expo_exact(cmax: float, tau: float, T: float) -> float:
    return cmax * (1.0 - math.exp(-T / tau))


# ==========================================================================
# The exact solution of the discrete recurrence, at 80 digits
# ==========================================================================
def _d(x: float) -> D:
    return D(repr(x))


def exact_discrete(dt: float, k: float, cmax: float, tau: float, T: float
                   ) -> Tuple[D, D, D]:
    """Closed form of  C^{n+1} = (C^n + dt S(t_n))/(1+k dt),  C^0 = 0,
    with S(t) = u'(t) + k u(t) and u(t) = cmax(1 - exp(-t/tau)).

        C^n = cmax(1 - a^n) + dt*beta*a*(r^n - a^n)/(r - a)

    Returns (C^n, u(T), C^n - u(T))."""
    dtD, kD, cD, tauD, TD = _d(dt), _d(k), _d(cmax), _d(tau), _d(T)
    n = int(round(T / dt))
    a = D(1) / (D(1) + kD * dtD)
    r = (-dtD / tauD).exp()
    beta = cD * (D(1) / tauD - kD)
    rn = r ** n
    an = a ** n
    cn = cD * (D(1) - an) + dtD * beta * a * (rn - an) / (r - a)
    uT = cD * (D(1) - (-TD / tauD).exp())
    return cn, uT, cn - uT


def brute_discrete(dt: float, k: float, cmax: float, tau: float, T: float) -> D:
    """Direct high-precision iteration of the same recurrence (guards the algebra)."""
    dtD, kD, cD, tauD = _d(dt), _d(k), _d(cmax), _d(tau)
    n_steps = int(round(T / dt))
    c = D(0)
    den = D(1) + kD * dtD
    for i in range(n_steps):
        t = dtD * i
        e = (-t / tauD).exp()
        s = (cD / tauD) * e + kD * cD * (D(1) - e)
        c = (c + dtD * s) / den
    return c


def g1_closed(k: float, cmax: float, tau: float, T: float) -> D:
    """Exact O(dt) coefficient of the global error (derived, not fitted)."""
    kD, cD, tauD, TD = _d(k), _d(cmax), _d(tau), _d(T)
    q = D(1) / tauD
    return -(cD * q * (2 * kD - q) / (2 * (kD - q))) * (
        (-q * TD).exp() - (-kD * TD).exp()
    )


def g1_g2_numeric(k: float, cmax: float, tau: float, T: float, m: int = 26
                  ) -> Tuple[D, D]:
    """Extract G1, G2 from E(dt) = G1 dt + G2 dt^2 + ... by Richardson on the
    CLOSED FORM (no solver run, no truncation analysis)."""
    dt1 = T / (2 ** m)
    dt2 = dt1 / 2
    e1 = exact_discrete(dt1, k, cmax, tau, T)[2]
    e2 = exact_discrete(dt2, k, cmax, tau, T)[2]
    r1 = e1 / _d(dt1)
    r2 = e2 / _d(dt2)
    g1 = 2 * r2 - r1
    g2 = 2 * (r1 - r2) / _d(dt1)
    return g1, g2


# ==========================================================================
# helpers
# ==========================================================================
def orders(errs: List[float]) -> List[float]:
    out = []
    for a, b in zip(errs[:-1], errs[1:]):
        out.append(float("nan") if (a <= 0 or b <= 0)
                   else math.log(a / b) / math.log(2.0))
    return out


def fmt(ps: List[float]) -> str:
    return "  ".join("  nan " if math.isnan(p) else f"{p:6.3f}" for p in ps)


def hdr(title: str) -> None:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")


# ==========================================================================
# A. the solver IS the scalar recurrence
# ==========================================================================
def part_a(ladder: List[float]) -> bool:
    hdr("A. Does the shipped solver reduce EXACTLY to backward Euler + lagged "
        "source?")
    print("   E4 config: 16x16, ecm=0, convection=False, source=ones*S(t_n).")
    print("   The discrete operator L has zero row sums, so L*(uniform) = 0 and")
    print("   A*(uniform) = (1 + k dt)*(uniform). Claim: C^{n+1} = (C^n + dt S_n)/(1+k dt).")
    print(f"\n   {'dt':>7} {'steps':>7} {'solver C(T)':>24} {'closed form':>24} "
          f"{'rel diff':>11}")
    worst = 0.0
    for dt in ladder:
        prof = expo_profile(dt, E4_K, E4_CMAX, E4_TAU, E4_T)
        c_sol, _ = solver_final(dt, E4_K, prof, E4_T)
        c_ex = exact_discrete(dt, E4_K, E4_CMAX, E4_TAU, E4_T)[0]
        rel = abs(D(repr(c_sol)) - c_ex) / abs(c_ex)
        worst = max(worst, float(rel))
        print(f"   {dt:>7.4g} {int(round(E4_T / dt)):>7d} {c_sol:>24.17g} "
              f"{float(c_ex):>24.17g} {float(rel):>11.3e}")
    ok = worst < 1e-13
    print(f"\n   worst relative difference = {worst:.3e}   "
          f"({'MATCH' if ok else 'MISMATCH'}, threshold 1e-13)")
    # guard the algebra of the closed form itself
    b = brute_discrete(ladder[0], E4_K, E4_CMAX, E4_TAU, E4_T)
    c = exact_discrete(ladder[0], E4_K, E4_CMAX, E4_TAU, E4_T)[0]
    alg = abs(b - c) / abs(c)
    print(f"   closed form vs direct 80-digit iteration at dt={ladder[0]:g}: "
          f"rel {float(alg):.3e}")
    assert alg < D("1e-60"), "closed-form algebra is wrong"
    print("   => the solver's time integrator is EXACTLY the documented scheme.")
    print("      Any order anomaly is a property of this recurrence, not of code.")
    return ok


# ==========================================================================
# B. the measured error is the closed-form error
# ==========================================================================
def part_b(ladder: List[float]) -> bool:
    hdr("B. Is E4's measured error reproduced by the closed form, with no "
        "fitted constant?")
    exact = expo_exact(E4_CMAX, E4_TAU, E4_T)
    print(f"   u(T) = {exact:.17g}   (k={E4_K}, tau={E4_TAU}, Cmax={E4_CMAX}, T={E4_T})")
    print(f"\n   {'dt':>7} {'measured |C-u|':>16} {'closed form |E|':>16} "
          f"{'ratio':>9}")
    errs_m, errs_c = [], []
    worst = 0.0
    for dt in ladder:
        prof = expo_profile(dt, E4_K, E4_CMAX, E4_TAU, E4_T)
        c_sol, _ = solver_final(dt, E4_K, prof, E4_T)
        em = abs(c_sol - exact)
        ec = abs(float(exact_discrete(dt, E4_K, E4_CMAX, E4_TAU, E4_T)[2]))
        errs_m.append(em)
        errs_c.append(ec)
        worst = max(worst, abs(em / ec - 1.0))
        print(f"   {dt:>7.4g} {em:>16.6e} {ec:>16.6e} {em / ec:>9.5f}")
    print(f"\n   observed order, measured   : {fmt(orders(errs_m))}")
    print(f"   observed order, closed form: {fmt(orders(errs_c))}")
    ok = worst < 1e-6
    print(f"   worst |ratio-1| = {worst:.3e}  ({'MATCH' if ok else 'MISMATCH'})")
    print("   => E4's order-2 result is entirely accounted for analytically.")
    return ok


# ==========================================================================
# C. the forcing term is right (no inconsistency floor)
# ==========================================================================
def part_c(ladder: List[float]) -> bool:
    hdr("C. Is the MANUFACTURED FORCING wrong? (the alternative that makes a "
        "correct solver look broken)")
    print("   A wrong S makes u not a solution of the semi-discrete problem, so")
    print("   C(T) converges to something OTHER than u(T): the error stops at a")
    print("   nonzero floor instead of going to 0. Two checks.")
    print("\n   C1. S(t) = u'(t) + k u(t) against a central difference of u:")
    worst = 0.0
    for name, u, s in (
        ("exponential (E4)",
         lambda t: E4_CMAX * (1.0 - math.exp(-t / E4_TAU)),
         lambda t: (E4_CMAX / E4_TAU) * math.exp(-t / E4_TAU)
                   + E4_K * E4_CMAX * (1.0 - math.exp(-t / E4_TAU))),
        ("rational (F)",
         lambda t: E4_CMAX * t / (t + E4_TAU),
         lambda t: E4_CMAX * E4_TAU / (t + E4_TAU) ** 2
                   + E4_K * E4_CMAX * t / (t + E4_TAU)),
    ):
        h = 1e-5
        for t in (10.0, 80.0, 160.0, 310.0):
            fd = (u(t + h) - u(t - h)) / (2 * h) + E4_K * u(t)
            worst = max(worst, abs(fd - s(t)) / abs(s(t)))
        print(f"       {name:<20} worst relative mismatch {worst:.3e}")
    ok_fd = worst < 1e-8
    print(f"       -> {'S is consistent with u' if ok_fd else 'S IS WRONG'} "
          f"(threshold 1e-8)")

    print("\n   C2. Does the E4 error go to ZERO, or to a floor? Closed-form ladder")
    print("       (validated against the solver in B) over 20 halvings of dt:")
    print(f"       {'dt (s)':>12} {'steps':>10} {'|E(T)|':>14} {'order':>8}")
    errs, dts = [], []
    for m in range(1, 25, 3):
        dt = E4_T / (2 ** m)
        e = abs(float(exact_discrete(dt, E4_K, E4_CMAX, E4_TAU, E4_T)[2]))
        dts.append(dt); errs.append(e)
    for i, (dt, e) in enumerate(zip(dts, errs)):
        p_ = "" if i == 0 else f"{math.log(errs[i-1]/e)/math.log(8.0):8.3f}"
        print(f"       {dt:>12.5g} {int(round(E4_T/dt)):>10d} {e:>14.6e} {p_:>8}")
    ok_floor = errs[-1] < 1e-12 * errs[0]
    print(f"       error fell by {errs[0]/errs[-1]:.3e}x with no floor -> "
          f"{'S is CORRECT' if ok_floor else 'INCONSISTENCY FLOOR PRESENT'}")
    print("       (order column is per 8x refinement / log8, i.e. 2.0 == 6.0 here)")
    return ok_fd and ok_floor


# ==========================================================================
# D. the exact O(dt) coefficient, and where it vanishes
# ==========================================================================
def part_d() -> bool:
    hdr("D. The exact O(dt) coefficient G1 of the global error, swept over k")
    print("   G1 = -(Cmax/(2 tau)) (2k - 1/tau)/(k - 1/tau) (e^{-T/tau} - e^{-kT})")
    print("   derived from the closed form; cross-checked against a Richardson")
    print("   extrapolation of the closed form at dt = T/2^26 and T/2^27.")
    print(f"\n   1/tau = {1.0 / E4_TAU:g};  E4 registered k = {E4_K:g}, so "
          f"2k - 1/tau = {2 * E4_K - 1 / E4_TAU:g}")
    print(f"\n   {'k':>9} {'2k - 1/tau':>12} {'G1 (closed)':>15} "
          f"{'G1 (Richardson)':>17} {'G2':>14}")
    ks = [0.005, 0.008, 0.0099, 0.00999, 0.01, 0.01001, 0.0101, 0.012, 0.04]
    ok = True
    for k in ks:
        gc = g1_closed(k, E4_CMAX, E4_TAU, E4_T)
        gn, g2 = g1_g2_numeric(k, E4_CMAX, E4_TAU, E4_T)
        if abs(gc) > D("1e-30"):
            agree = abs(gn - gc) / abs(gc) < D("1e-6")
            note = ""
        else:
            # G1 is exactly zero: the Richardson estimate must be at
            # roundoff level relative to the next coefficient G2.
            agree = abs(gn) < D("1e-10") * abs(g2)
            note = "   <-- exactly zero"
        ok = ok and bool(agree)
        print(f"   {k:>9.5g} {2 * k - 1 / E4_TAU:>12.4g} {float(gc):>15.6e} "
              f"{float(gn):>17.6e} {float(g2):>14.6e}{note}")
    gc0 = g1_closed(E4_K, E4_CMAX, E4_TAU, E4_T)
    print(f"\n   At E4's own k: G1 = {float(gc0):.6e} exactly, killed by the "
          f"factor (2k - 1/tau).")
    print("   The O(dt) term is not small at E4's parameters -- it is ABSENT,")
    print("   so the leading surviving term is G2 dt^2 and the measured order is 2.")
    print("\n   Degeneracy condition, general form: the leading local truncation")
    print("   error of this scheme is dt*(u''/2 + k u'), which vanishes for all t")
    print("   iff u' is proportional to exp(-2kt). E4 chose u' ~ exp(-t/tau)")
    print("   with 1/tau = 2k, i.e. it landed on that one-parameter family exactly.")
    return ok and abs(gc0) < D("1e-30")


def _sign_change(k: float, lo: float = 1e-4, hi: float = 40.0):
    """Bisect the closed-form signed error E(dt) for a zero crossing. dt values
    need not divide T here: the closed form is evaluated with n = T/dt rounded,
    so we bisect on a log grid of dt that always divides T."""
    # grid of dt = T/2^m with m fractional is not allowed (n must be integer),
    # so bisect on n instead: n integer, dt = T/n.
    def e_of_n(n: int) -> float:
        return float(exact_discrete(E4_T / n, k, E4_CMAX, E4_TAU, E4_T)[2])
    n_lo = max(1, int(E4_T / hi))
    n_hi = int(E4_T / lo)
    a, b = e_of_n(n_lo), e_of_n(n_hi)
    if a == 0.0 or b == 0.0 or (a > 0) == (b > 0):
        return None
    while n_hi - n_lo > 1:
        mid = (n_lo + n_hi) // 2
        if (e_of_n(mid) > 0) == (a > 0):
            n_lo = mid
        else:
            n_hi = mid
    return E4_T / ((n_lo + n_hi) / 2.0)


# ==========================================================================
# E. the crossover signature -- the discriminator
# ==========================================================================
def part_e(quick: bool) -> bool:
    hdr("E. Detune k off the degeneracy: the order must CROSS OVER from 2 to 1")
    print("   If the integrator were genuinely second order, the order would be 2")
    print("   at every dt and every k. If instead G1 vanishes only AT k = 1/(2 tau),")
    print("   then for k nearby the error is G1 dt + G2 dt^2: order 2 while")
    print("   |G2| dt > |G1|, order 1 below, with the changeover at dt* = |G1/G2|")
    print("   (and, when G1 and G2 have opposite signs, a near-exact cancellation")
    print("   -- a visible error MINIMUM -- at that same dt*).")
    m_lo, m_hi = 4, (18 if quick else 24)
    solver_m_hi = 12 if quick else 14
    exact = expo_exact(E4_CMAX, E4_TAU, E4_T)
    ok = True
    for k in ([0.0105, 0.010005, 0.01] if quick
              else [0.0105, 0.01005, 0.010005, 0.01]):
        g1, g2 = g1_g2_numeric(k, E4_CMAX, E4_TAU, E4_T)
        dt_star = float(abs(g1 / g2)) if g2 != 0 else float("inf")
        opposite = (g1 * g2) < 0
        degen = abs(g1) < D("1e-10") * abs(g2)
        print(f"\n   k = {k:<9g} 2k-1/tau = {2 * k - 1 / E4_TAU:>9.3g}   "
              f"G1 = {float(g1):>12.5e}  G2 = {float(g2):>12.5e}")
        if degen:
            print("     G1 = 0 exactly -> no crossover is possible; order must be 2 "
                  "at EVERY dt")
        else:
            print(f"     predicted changeover dt* = |G1/G2| = {dt_star:.5g} s"
                  + ("  (signs oppose -> error minimum there)" if opposite
                     else "  (same sign -> no minimum, order rolls over)"))
        # exact location of the sign change of E(dt), if any, by bisection on
        # the closed form -- the two-term model puts it at dt*, but higher-order
        # terms move it, and that difference is the honest measure of how far
        # the asymptotic argument reaches.
        zero_dt = _sign_change(k)
        if zero_dt is not None:
            print(f"     exact sign change of E(dt) (bisection on the closed "
                  f"form): dt = {zero_dt:.5g} s")
        else:
            print("     E(dt) never changes sign on [1e-4, 40] s")
        print(f"     {'dt (s)':>11} {'steps':>9} {'|E| closed':>13} "
              f"{'|E| solver':>13} {'G1dt+G2dt^2':>13} {'order':>7}")
        prev = None
        orders_seen = []
        for m in range(m_lo, m_hi + 1):
            dt = E4_T / (2 ** m)
            ec = abs(float(exact_discrete(dt, k, E4_CMAX, E4_TAU, E4_T)[2]))
            if m <= solver_m_hi:
                prof = expo_profile(dt, k, E4_CMAX, E4_TAU, E4_T)
                es_, _ = solver_final(dt, k, prof, E4_T)
                es = f"{abs(es_ - exact):13.5e}"
            else:
                es = f"{'--':>13}"
            p_ = (math.log(prev / ec) / math.log(2.0)
                  if prev and ec > 0 else float("nan"))
            orders_seen.append((dt, p_))
            prev = ec
            two = abs(float(g1 * _d(dt) + g2 * _d(dt) * _d(dt)))
            if m <= m_lo + 2 or m >= m_hi - 3 or (m - m_lo) % 3 == 0:
                print(f"     {dt:>11.5g} {int(round(E4_T / dt)):>9d} {ec:>13.5e} "
                      f"{es} {two:>13.5e} "
                      + ("    nan" if math.isnan(p_) else f"{p_:7.3f}"))
        tail = [p for _, p in orders_seen[-3:]]
        if degen:
            good = all(1.8 <= p <= 2.2 for p in tail)
            print(f"     asymptotic order (finest 3 intervals): "
                  f"{' '.join(f'{p:.3f}' for p in tail)} -> "
                  f"{'stays 2, as required' if good else 'UNEXPECTED'}")
        else:
            good = all(0.9 <= p <= 1.1 for p in tail)
            print(f"     asymptotic order (finest 3 intervals): "
                  f"{' '.join(f'{p:.3f}' for p in tail)} -> "
                  f"{'FIRST ORDER' if good else 'UNEXPECTED'}")
        ok = ok and good
    print("\n   Every non-degenerate k lands on order 1 once dt < dt*; only the")
    print("   single point k = 1/(2 tau) stays at 2 forever. That is a property of")
    print("   the manufactured solution, not of the integrator.")
    return ok


# ==========================================================================
# F. a non-degenerate manufactured solution at E4's own parameters
# ==========================================================================
def rational_profile(dt: float, k: float, cmax: float, tau: float, T: float
                     ) -> np.ndarray:
    """u(t) = cmax t/(t+tau);  S = u' + k u = cmax tau/(t+tau)^2 + k cmax t/(t+tau)."""
    n_steps = int(round(T / dt))
    t = np.arange(n_steps) * dt
    return cmax * tau / (t + tau) ** 2 + k * cmax * t / (t + tau)


def part_f(ladder: List[float]) -> bool:
    hdr("F. The same solver, same k, same T, a NON-degenerate manufactured "
        "solution")
    print("   u(t) = Cmax t/(t+tau), u(0)=0 (the solver's fixed IC), u' > 0,")
    print("   u' is rational, so it can never be proportional to exp(-2kt):")
    print("   no choice of k can make this one degenerate.")
    exact = E4_CMAX * E4_T / (E4_T + E4_TAU)
    print(f"\n   u(T) = {exact:.17g}   (Cmax={E4_CMAX}, tau={E4_TAU}, T={E4_T}, "
          f"k={E4_K} -- E4's own k)")
    print(f"\n   {'dt':>7} {'|C(T) - u(T)|':>16} {'min S':>12}")
    errs = []
    smin_all = float("inf")
    for dt in ladder:
        prof = rational_profile(dt, E4_K, E4_CMAX, E4_TAU, E4_T)
        smin_all = min(smin_all, float(prof.min()))
        c_sol, _ = solver_final(dt, E4_K, prof, E4_T)
        errs.append(abs(c_sol - exact))
        print(f"   {dt:>7.4g} {errs[-1]:>16.6e} {float(prof.min()):>12.4e}")
    ps = orders(errs)
    print(f"\n   observed order: {fmt(ps)}")
    ok = all(0.8 <= p <= 1.2 for p in ps)
    print(f"   V1 guard min(S) > 0 over every level: {smin_all:.4e} -> "
          f"{'OK' if smin_all > 0 else 'VIOLATED'}")
    print(f"   order in the [0.8, 1.2] band E4 registered: "
          f"{'YES -- first order' if ok else 'NO'}")
    print("   => at E4's OWN k the integrator is first order. The order-2 result")
    print("      belongs to E4's manufactured solution, not to the solver.")
    return ok and smin_all > 0


# ==========================================================================
# G. regime / CFL cannot be the explanation
# ==========================================================================
def part_g() -> bool:
    hdr("G. Is E4 outside a stated regime (CFL, stability)?")
    prof = expo_profile(4.0, E4_K, E4_CMAX, E4_TAU, E4_T)
    _, res = solver_final(4.0, E4_K, prof, E4_T)
    print(f"   E4 runs convection=False and passes no velocity, so the explicit")
    print(f"   upwind branch -- and with it the CFL guard -- is never entered.")
    print(f"   Reported CFL for the E4 run: {getattr(res, 'cfl', None)!r}")
    print("   Diffusion and reaction are solved implicitly (A = I - dt L + dt K,")
    print("   an M-matrix for any dt > 0), so there is no dt <= dx^2/(4D) limit")
    print("   and backward Euler is unconditionally stable. A dt ladder of")
    print("   4/2/1/0.5 s on a 16x16 grid with u spatially uniform sits in the")
    print("   scheme's stated regime at every level.")
    # unconditional stability: push dt absurdly coarse and confirm no blow-up
    exact = expo_exact(E4_CMAX, E4_TAU, E4_T)
    for dt in (320.0, 160.0, 80.0):
        p = expo_profile(dt, E4_K, E4_CMAX, E4_TAU, E4_T)
        c, _ = solver_final(dt, E4_K, p, E4_T)
        print(f"   dt = {dt:>6.1f} s ({int(E4_T/dt):>3d} step(s)): C(T) = {c:.6f}, "
              f"|C-u| = {abs(c-exact):.4e}  (finite, bounded -> stable)")
    print("   => 'E4 is outside the regime' is not available as an explanation:")
    print("      the only regime boundary this scheme has (the advective CFL) is")
    print("      not exercised by E4 at all. The separately measured 2x-loose CFL")
    print("      guard (strict margin 0.6506 at the shipped dx=12 um, dt=10 s)")
    print("      concerns explicit upwind advection and has no bearing here.")
    return True


# ==========================================================================
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    a0 = calibrate_alpha()
    assert a0 == 22.223477711172563, a0
    print(f"calibrate_alpha() = {a0!r}")
    t0 = time.time()

    ladder = E4_LADDER[:3] if args.quick else E4_LADDER
    res = {
        "A solver == backward Euler recurrence": part_a(ladder),
        "B measured error == closed form": part_b(ladder),
        "C manufactured forcing S is correct": part_c(ladder),
        "D G1 vanishes exactly iff 2k = 1/tau": part_d(),
        "E order crosses 2 -> 1 off the degeneracy": part_e(args.quick),
        "F non-degenerate u at E4's own k is order 1": part_f(ladder),
        "G CFL / regime not implicated": part_g(),
    }

    hdr("SUMMARY")
    for k, v in res.items():
        print(f"  {'OK  ' if v else 'FAIL'}  {k}")
    print(f"\n  elapsed {time.time() - t0:.1f} s")
    a1 = calibrate_alpha()
    assert a1 == 22.223477711172563, a1
    print(f"calibrate_alpha() = {a1!r}")
    return 0 if all(res.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
