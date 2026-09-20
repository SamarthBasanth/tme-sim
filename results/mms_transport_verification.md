# Method of manufactured solutions on the transport solver

**Headline: 4 of 5 pre-registered cases PASS, and E4 FAILED.**
The spatial discretisation is verified at its formal order — second order for
diffusion with constant D, with variable D, at the domain-edge no-flux boundary
and at the tissue-mask no-flux boundary; first order for the donor-cell upwind
advection, which is the order that scheme is *supposed* to have. The temporal
case **E4 FAILED its pre-registered band**: it was registered to show order 1
and measured **2.019 / 2.009 / 2.005**. §6 shows the cause is a degeneracy in
*my manufactured solution* (`2k = 1/tau` annihilates the leading truncation
term), not a solver defect — but the pre-registered verdict stands as FAILED and
no substitute metric has been put in its place.

Driver: `experiments/mms_transport_verification.py`.
Raw run: reproduced verbatim in §5 and §7 below.
`calibrate_alpha() = 22.223477711172563` asserted at the start and the end of
every run quoted in this file.

> **Why this file exists.** The repository already had grid-convergence work —
> `experiments/grid_convergence_gbm.py`, `experiments/grid_convergence_breast.py`,
> `tmesim/validation/numerical_rigor.py` — but every one of those is Richardson
> extrapolation on a *derived quantity of interest* (penetration index, kill
> fraction) against **no known answer**. That tells you whether an answer is
> still moving with the mesh. It cannot tell you whether the PDE solver solves
> the PDE it claims to. A search of the tree for `manufactured`, `MMS`,
> `exact solution`, `analytic_solution` and `observed order` before this work
> returned only Richardson-on-QoI plus one passing comment in
> `experiments/boundary_layer_resolution.py`. **There was no solver
> verification in this repository.** There is now.
>
> MMS is the right instrument here because it needs no tissue, no patient, no
> anchor and no calibration argument. Its verdict is independent of every
> accuracy dispute in `results/`.

---

## 1. Pre-registration

Written and frozen **before any error was measured**
(`PREREG_mms.md`, 2026-09-20 03:08:03; first measurement 03:09). Reproduced
verbatim; §6 and §7 carry two later, separately-dated pre-registrations, each
labelled where it appears.

> ### Equation as implemented (read from source, not from docs)
> ```
>     dC/dt = div(D_eff grad C) - div(v C) + S - k C
>     D_eff = D * (1 - 0.6*ecm), floored at max(1e-6, 1e-3*D)
>     k     = clearance + binding (+ exchange_k when supplied)
>     BC    : zero TOTAL flux on every domain-edge face and every
>             tissue/non-tissue face (face coefficients set to 0)
>     IC    : C == 0
>     Scheme: IMEX. Backward Euler on diffusion+reaction (implicit, one LU reuse),
>             first-order donor-cell upwind on advection (explicit),
>             source evaluated at the OLD time level,
>             post-step clamp C <- max(C,0) and C <- 0 outside tissue.
>     Grid  : cell-centred, x_j=(j+1/2)dx, extent nx*dx (Domain.coords_um).
> ```
>
> ### Formal orders claimed by that scheme
> ```
>   spatial diffusion : 2   (5-point FV, harmonic face conductance, flux-zero BC)
>   spatial advection : 1   (donor-cell upwind)
>   temporal          : 1   (backward Euler, source lagged)
> ```
>
> ### Manufactured solutions
> Steady cases use `u(x,y) = A + B cos(kx x) cos(ky y)`, A=2, B=1, so u in
> [1,3]>0 and du/dn = 0 EXACTLY on every no-flux face. Source is derived, never
> fitted:
> ```
>       S = k u - div(D_eff grad u) + div(v u)
> ```
> Temporal case uses the spatially UNIFORM `u(t) = Cmax(1 - exp(-t/tau))`,
> Cmax=10, tau=50 s, u(0)=0 (the solver's fixed IC), `S(t) = u'(t) + k u(t)`
> passed as `source=ones * source_profile=S(t_n)`.
> ```
>   E1 diffusion, constant D:  L=200, kx=ky=2pi/L, D=300, k=1.0, dt=1, T=100
>   E2 diffusion, variable D:  as E1 but ecm = 0.4 + 0.3 sin(2pi x/L) sin(2pi y/L)
>   E3 convection+diffusion:   L=200, D=2, k=0.1, v=(V sin(pi x/L), V sin(pi y/L)),
>                              V=0.2 um/s (v.n = 0 on every edge), dt=1, T=600
>   E4 temporal:               16x16 grid, k=0.01, T=320 s
>   E5 tissue-mask no-flux:    L=200, tissue = [40,160]^2 (face-aligned on every
>                              dx of the ladder), kx=ky=2pi/120, D=300, k=2.0
> ```
>
> ### Ladders (refinement ratio r = 2)
> ```
>   E1,E2,E3,E5: dx = 4, 2, 1, 0.5 um  (nx = 50, 100, 200, 400) -> 3 observed orders
>   E4:          dt = 4, 2, 1, 0.5 s   (T = 320 s)               -> 3 observed orders
> ```
>
> ### Norms (reported per level AND as per-step observed order, not only a slope)
> `L2 = sqrt(mean over the norm set of (C-u)^2)` (RMS);
> `Linf = max over the norm set of |C-u|`.
> Norm sets: GLOBAL (all tissue cells); BOUNDARY (tissue cells with at least one
> off-domain or non-tissue neighbour); INTERIOR (the rest).
> Observed order between consecutive levels: `p = log(E_coarse/E_fine)/log(2)`.
>
> ### PASS / FAIL thresholds — fixed now
> * **E1 PASS iff EVERY consecutive p in [1.8, 2.2]** for ALL FOUR of: global L2,
>   global Linf, boundary Linf, interior Linf.
> * **E2 PASS**: identical thresholds to E1.
> * **E5 PASS**: identical thresholds to E1.
> * **E3 PASS iff EVERY consecutive p in [0.8, 1.2]** for global L2 and global
>   Linf. An order near 2 is ALSO a FAIL: it would mean the advective flux is not
>   the donor-cell upwind the module documents.
> * **E4 PASS iff EVERY consecutive p in [0.8, 1.2]** on `|C(T) - u(T)|`.
>
> Anything outside these intervals is reported as FAILED in the headline. No
> substitute metric will be introduced.
>
> ### Validity guards — a case is VOID (neither PASS nor FAIL) if one fails
> * **V1** `min(S) > 0` on every grid. With S >= 0, `A = I - dt L + dt K` is an
>   M-matrix and upwind is positivity-preserving at CFL <= 1, so C >= 0 for every
>   step and the solver's `max(C,0)` clamp provably never binds. If it could
>   bind, the solver would not be solving the linear PDE being verified.
> * **V2** steady-state residual `max|C(T) - C(T/2)| < 1e-3 * Linf` error at that
>   level, so the reported error is spatial, not leftover transient. (E1,E2,E3,E5.)
> * **V3** the closed-form S agrees with a nested finite-difference evaluation of
>   `div(D grad u)` and `div(v u)` to < 1e-6 relative. Guards MY algebra: a wrong
>   S makes u not a solution and the order meaningless.
> * **V4** reported CFL <= 1 on every E3 level.
>

**Correction to V1's justification, added 2026-09-20 AFTER the fact — the
pre-registration above is deliberately left unedited.**

V1 argues that "upwind is positivity-preserving at CFL <= 1, so the solver's
`max(C,0)` clamp provably never binds". The M-matrix half is right. The upwind
half is **false in general**: the solver's CFL number is a maximum over CELLS of
the pointwise sum `|vx| + |vy|`, while positivity requires the two component
maxima taken INDEPENDENTLY and allows only half as much. A divergent cross —
four neighbours each flowing away, no single cell carrying both components —
separates the two quantities by exactly 2x, and at CFL = 1 the undershoot is
-100% of the peak. Measured and reproduced to 4.0e-07 by
`experiments/conservation_audit.py` on 2026-09-20.

**Nothing in this document's results changes.** Its E3 velocity field is a sine
whose two component maxima are co-located, giving a measured ratio of 1.0000
against the counterexample's 2.0000, and its own V1 check recorded `min(S) > 0`
with no undershoot at CFL up to 1. The guard did its job on the fields it ran
on; its stated reason was broader than the thing it establishes.

This is appended rather than corrected in place because a pre-registration that
is edited after the data are in is not a pre-registration.

> ### Plant tests (the check must go red when its subject is broken)
> ```
>   P1 -> E1/E5 boundary: no-flux mirrored about the boundary cell CENTRE instead
>         of the cell FACE (edge face coefficient doubled into the diagonal).
>   P2 -> E2 averaging: harmonic face diffusivity replaced by the one-sided
>         (donor-cell) value. E1 is expected to stay GREEN under P2 (constant D
>         makes the two identical) — that is the point: E1 does not test averaging.
>   P3 -> E3 advection: _advective_divergence returns zeros.
>   P4 -> E4 source/time: the source contribution scaled by 0.99.
> ```
> Each plant is first inserted as a bare raise to prove the line is REACHED,
> then as the defect, then reverted.
>
> ### What a PASS would NOT show
> MMS verifies the SOLVER against the equation the module says it solves. It says
> nothing about whether that equation describes tumour tissue, nothing about any
> parameter value, and nothing about any published result in `results/`.

---

## 2. The equation, as implemented

Read out of `tmesim/physics/diffusion.py` rather than out of the docs. The
module docstring and the code agree; this is recorded so the reader can check
what was verified against what.

```
dC/dt = div(D_eff grad C) - div(v C) + S - k C          (uM/s)
```

| Piece | Implementation | Formal order |
|---|---|---|
| `div(D_eff grad C)` | `_build_diffusion_operator`: 5-point finite volume, face conductance = **harmonic mean** of the two adjacent cell `D_eff`, divided by `dx^2`; symmetric, rows sum to zero | **2** in space |
| `D_eff` | `_effective_diffusivity`: `D*(1 - 0.6*ecm)`, ecm clipped to [0,1], floored at `max(1e-6, 1e-3 D)` | — |
| `div(v C)` | `_advective_divergence`: face velocity = arithmetic mean of the two cell velocities, **donor-cell (upwind)** concentration, explicit | **1** in space |
| `S` | added to the RHS at the **old** time level, optionally scaled per step by `source_profile[n]` | **1** in time |
| `k C` | `clearance + binding (+ exchange_k)`, implicit | exact (linear) |
| Boundaries | **zero total flux**: the face coefficient is set to 0 on every domain-edge face and every tissue/non-tissue face. On a cell-centred grid this places the no-flux condition exactly **at the face**, i.e. at `x = 0` and `x = nx*dx` | **2** |
| Time | backward Euler on diffusion+reaction, explicit upwind advection (IMEX); `A = I - dt L + dt K` LU-factorised once | **1** in time |
| Initial condition | `C == 0`, not settable through the public API | — |
| Post-step | `C <- max(C, 0)`; `C <- 0` outside the tissue mask | nonlinear — see V1 |

One thing worth recording because it could have gone either way and did not:
`Domain.coords_um()` uses `x_j = (j + 1/2) dx` with extent `nx*dx`. That is
**exactly** the convention the flux-zeroing boundary implies. Had the repo used
`x_j = j*dx`, the manufactured solution's `du/dn = 0` would have been imposed
half a cell away from where the scheme imposes it, and every case below would
have degraded to first order. The two conventions agree here.

---

## 3. Method

Five cases. Four are **steady** manufactured solutions, which is deliberate: the
fixed point of the solver's own iteration is

```
C* = (K - L)^{-1} S
```

— the discrete **steady** problem, **exactly independent of dt**. So a converged
steady run isolates the spatial discretisation with the temporal error *removed*
rather than merely bounded. Convergence to that fixed point is not assumed; it
is measured (guard V2, the `resid/err` column).

| Case | Subject | Boundary exercised | Expected order |
|---|---|---|---|
| **E1** | diffusion, **constant** D | domain edge | 2 |
| **E2** | diffusion, **variable** D (smooth ECM field) — exercises the harmonic face averaging | domain edge | 2 |
| **E3** | **convection** + diffusion, `v.n = 0` on every edge | domain edge | 1 |
| **E4** | **time** integration + source convention (spatially uniform u, so `L u = 0` exactly — the spatial operator is annihilated by construction) | none | 1 |
| **E5** | diffusion, constant D, on a **tissue-mask sub-block** `[40,160]^2` | tissue/non-tissue mask | 2 |

The manufactured solutions satisfy `du/dn = 0` on every no-flux face, which is
the only boundary condition the solver can impose, and are strictly positive
(`u in [1,3]`). Positivity matters: with `S >= 0` everywhere (guard V1), `A` is
an M-matrix and upwind is positivity-preserving under CFL <= 1, so `C >= 0` at
every step and the solver's `max(C,0)` clamp **provably never binds**. The
measured `min C` column confirms it: it never drops below **1.000** in any of
the five clean cases.

For E3 the velocity is `v = (V sin(pi x/L), V sin(pi y/L))`, whose normal
component vanishes on all four edges — so the solver zeroing the boundary
advective flux is the *correct* condition for this solution rather than an
inconsistency introduced by the test. E3 also runs at `D = 2 um^2/s` so that the
O(dx) advective truncation dominates the O(dx^2) diffusive one by ~10^2 across
the ladder; otherwise E3 would read as order 2 and would be testing diffusion
while claiming to test advection.

`S` is derived in closed form and checked against a nested finite-difference
evaluation of the same continuous operators (guard V3) — worst relative
disagreement across all cases **5.7e-07**, against a 1e-06 tolerance. Nothing in
this study is fitted to anything.

---

## 4. Verdict per case

| Case | Subject | Pre-registered band | Measured (3 consecutive steps) | Verdict |
|---|---|---|---|---|
| E1 | diffusion, constant D | [1.8, 2.2] | Linf 1.996 / 1.999 / 2.000 | **PASS** |
| E2 | diffusion, variable D | [1.8, 2.2] | Linf 1.997 / 1.999 / 2.000 | **PASS** |
| E3 | upwind advection | [0.8, 1.2] | Linf 0.986 / 0.994 / 0.997 | **PASS** |
| E4 | time + source | [0.8, 1.2] | 2.019 / 2.009 / 2.005 | **FAILED** |
| E5 | tissue-mask boundary | [1.8, 2.2] | Linf 1.988 / 1.997 / 1.999 | **PASS** |

**Boundary vs interior, answered separately as required.** In E1, E2 and E5 the
boundary-restricted Linf order is second order and is never the worse of the
two (E1: 1.996/1.999/2.000 in both; E2: boundary 1.997/1.999/2.000, interior
1.981/1.995/1.999, i.e. the boundary is marginally BETTER; E5: identical in
both, 1.988/1.997/1.999). There is **no order loss at
either no-flux boundary** — neither the domain edge nor the tissue mask. §7
shows this split is not vacuous: plant P1 drops the boundary Linf order from
1.996 to **0.964** while the check is looking at it.

---

## 5. Raw table

Verbatim from a single run of the unmodified tree,
`.venv/bin/python experiments/mms_transport_verification.py`. `resid/err` is
guard V2 (steady-state residual over the last half window, divided by the Linf
error at that level); `min S` is guard V1; `min C` shows the clamp never bound;
`CFL` is guard V4.

```
calibrate_alpha() = 22.223477711172563

==============================================================================
E1: pure diffusion, CONSTANT D, domain-edge no-flux (expect order 2)
==============================================================================
  V3 algebra check (closed-form S vs nested FD): div(D grad u) rel 2.055e-07, div(vu) rel 0.000e+00  [tol 1e-06] -> OK

      dx     n      L2(all)    Linf(all)    Linf(bnd)    Linf(int)  resid/err     min S     min C    CFL      s
       4    50  2.44710e-04  4.87491e-04  4.87491e-04  4.87491e-04   6.38e-12    0.4141     1.003  0.000    0.0
       2   100  6.11793e-05  1.22238e-04  1.22238e-04  1.22238e-04   2.91e-11    0.4094     1.001  0.000    0.1
       1   200  1.52949e-05  3.05823e-05  3.05823e-05  3.05823e-05   1.31e-10    0.4082         1  0.000    0.5
     0.5   400  3.82374e-06  7.64700e-06  7.64700e-06  7.64700e-06   6.97e-10    0.4079         1  0.000    2.7

  observed order per refinement step (r = 2), pre-registered band [1.8, 2.2]:
    L2   (global)    4->2:  2.000  2->1:  2.000  1->0.5:  2.000   [PASS]
    Linf (global)    4->2:  1.996  2->1:  1.999  1->0.5:  2.000   [PASS]
    Linf (boundary)  4->2:  1.996  2->1:  1.999  1->0.5:  2.000   [PASS]
    Linf (interior)  4->2:  1.996  2->1:  1.999  1->0.5:  2.000   [PASS]
    L2   (boundary)  4->2:  1.990  2->1:  1.996  1->0.5:  1.998   [ok]

  guards: V1 min(S)>0 True (min C>0 True)   V2 residual<0.001*Linf True   V3 algebra True   V4 CFL<=1 True
  E1 VERDICT: PASS

==============================================================================
E2: pure diffusion, VARIABLE D (harmonic face averaging) (expect order 2)
==============================================================================
  V3 algebra check (closed-form S vs nested FD): div(D grad u) rel 2.275e-07, div(vu) rel 0.000e+00  [tol 1e-06] -> OK

      dx     n      L2(all)    Linf(all)    Linf(bnd)    Linf(int)  resid/err     min S     min C    CFL      s
       4    50  2.12242e-04  4.78644e-04  4.78644e-04  4.71313e-04   7.42e-12    0.5548     1.004  0.000    0.0
       2   100  5.30582e-05  1.19871e-04  1.19871e-04  1.19394e-04   2.96e-11    0.5512     1.001  0.000    0.1
       1   200  1.32644e-05  2.99814e-05  2.99814e-05  2.99510e-05   1.33e-10    0.5503         1  0.000    0.5
     0.5   400  3.31609e-06  7.49623e-06  7.49623e-06  7.49431e-06   6.52e-10      0.55         1  0.000    2.8

  observed order per refinement step (r = 2), pre-registered band [1.8, 2.2]:
    L2   (global)    4->2:  2.000  2->1:  2.000  1->0.5:  2.000   [PASS]
    Linf (global)    4->2:  1.997  2->1:  1.999  1->0.5:  2.000   [PASS]
    Linf (boundary)  4->2:  1.997  2->1:  1.999  1->0.5:  2.000   [PASS]
    Linf (interior)  4->2:  1.981  2->1:  1.995  1->0.5:  1.999   [PASS]
    L2   (boundary)  4->2:  1.991  2->1:  1.996  1->0.5:  1.998   [ok]

  guards: V1 min(S)>0 True (min C>0 True)   V2 residual<0.001*Linf True   V3 algebra True   V4 CFL<=1 True
  E2 VERDICT: PASS

==============================================================================
E3: CONVECTION + diffusion, donor-cell upwind (expect order 1)
==============================================================================
  V3 algebra check (closed-form S vs nested FD): div(D grad u) rel 2.055e-07, div(vu) rel 5.166e-07  [tol 1e-06] -> OK

      dx     n      L2(all)    Linf(all)    Linf(bnd)    Linf(int)  resid/err     min S     min C    CFL      s
       4    50  2.72226e-03  7.32459e-03  4.44152e-03  7.32459e-03   5.63e-10     0.093     1.008  0.100    0.1
       2   100  1.37230e-03  3.69881e-03  2.20543e-03  3.69881e-03   1.13e-09   0.09286     1.003  0.200    0.6
       1   200  6.88855e-04  1.85733e-03  1.09957e-03  1.85733e-03   2.27e-09   0.09281     1.001  0.400    3.3
     0.5   400  3.45093e-04  9.30547e-04  5.49204e-04  9.30547e-04   4.54e-09   0.09281     1.001  0.800   13.9

  observed order per refinement step (r = 2), pre-registered band [0.8, 1.2]:
    L2   (global)    4->2:  0.988  2->1:  0.994  1->0.5:  0.997   [PASS]
    Linf (global)    4->2:  0.986  2->1:  0.994  1->0.5:  0.997   [PASS]

  guards: V1 min(S)>0 True (min C>0 True)   V2 residual<0.001*Linf True   V3 algebra True   V4 CFL<=1 True
  E3 VERDICT: PASS

==============================================================================
E4: TEMPORAL order, spatially uniform exact solution (expect order 1)
==============================================================================
  u(t) = 10.0 (1 - exp(-t/50)), S(t) = u'(t) + k u(t), k = 0.01, T = 320 s

      dt    |C(T)-u(T)|      min S      s
       4   4.280997e-04     0.1002    0.0
       2   1.056540e-04     0.1002    0.0
       1   2.624076e-05     0.1002    0.0
     0.5   6.538511e-06     0.1002    0.0

  observed order per step (r = 2), pre-registered band [0.8, 1.2]:
    |C(T)-u(T)|      4->2:  2.019  2->1:  2.009  1->0.5:  2.005   [FAIL]

  guards: V1 min(S)>=0 True
  E4 VERDICT: FAILED

==============================================================================
E5: TISSUE-MASK internal no-flux boundary, constant D (expect order 2)
==============================================================================
  V3 algebra check (closed-form S vs nested FD): div(D grad u) rel 5.710e-07, div(vu) rel 0.000e+00  [tol 1e-06] -> OK

      dx     n      L2(all)    Linf(all)    Linf(bnd)    Linf(int)  resid/err     min S     min C    CFL      s
       4    50  8.24985e-04  1.63194e-03  1.63194e-03  1.63194e-03   0.00e+00    0.3949     1.009  0.000    0.0
       2   100  2.06217e-04  4.11305e-04  4.11305e-04  4.11305e-04   0.00e+00     0.365     1.002  0.000    0.0
       1   200  5.15525e-05  1.03034e-04  1.03034e-04  1.03034e-04   0.00e+00    0.3576     1.001  0.000    0.2
     0.5   400  1.28880e-05  2.57716e-05  2.57716e-05  2.57716e-05   0.00e+00    0.3557         1  0.000    1.0

  observed order per refinement step (r = 2), pre-registered band [1.8, 2.2]:
    L2   (global)    4->2:  2.000  2->1:  2.000  1->0.5:  2.000   [PASS]
    Linf (global)    4->2:  1.988  2->1:  1.997  1->0.5:  1.999   [PASS]
    Linf (boundary)  4->2:  1.988  2->1:  1.997  1->0.5:  1.999   [PASS]
    Linf (interior)  4->2:  1.988  2->1:  1.997  1->0.5:  1.999   [PASS]
    L2   (boundary)  4->2:  1.982  2->1:  1.992  1->0.5:  1.997   [ok]

  guards: V1 min(S)>0 True (min C>0 True)   V2 residual<0.001*Linf True   V3 algebra True   V4 CFL<=1 True
  E5 VERDICT: PASS

==============================================================================
SUMMARY
==============================================================================
  E1: PASS
  E2: PASS
  E3: PASS
  E4: FAILED
  E5: PASS

  OVERALL: FAILED
calibrate_alpha() = 22.223477711172563
```

(E5's `resid/err` is exactly `0.00e+00`: with `k = 2.0` and `dt = 1` the slowest
mode is damped by `3^-50 ~ 1e-24` per half-window, so `C(T)` and `C(T/2)` are
bit-identical. `0 < 1e-3 * Linf` satisfies V2.)

---

## 6. The E4 failure, and its diagnosis

**E4 FAILED.** Registered band [0.8, 1.2]; measured 2.019 / 2.009 / 2.005. That
verdict is not revised below and no other metric is substituted for it.

A second pre-registration (`PREREG_E4b.md`, 2026-09-20 03:11:54, written after
seeing the failure and **before** running any diagnostic) states the hypothesis
and what would falsify it. For the scheme

```
(C^{n+1} - C^n)/dt = S(t_n) - k C^{n+1}        exact:  u' = S(t) - k u
```

the leading local truncation error is `(dt/2)u'' + dt k u' = (dt/2)(S' + k u')`.
With `u = Cmax(1 - exp(-t/tau))` and `S = u' + k u`,

```
S' + k u' = (Cmax/tau) exp(-t/tau) * (2k - 1/tau).
```

The pre-registered E4 used `k = 0.01`, `tau = 50`, i.e. `2k = 0.02 = 1/tau`
**exactly**. The leading term is identically zero for all t, so the scheme is
accidentally second order *on that one problem*. The registered predictions
were: order in [0.8,1.2] whenever `2k != 1/tau`; order ~2 at `k = 1/(2 tau)`;
and `err(dt=4)` roughly proportional to `|2k - 1/tau|`.

```
==============================================================================
E4b (POST-HOC, not pre-registered): temporal order vs the degeneracy 2k = 1/tau
==============================================================================
  leading LTE coefficient ~ (Cmax/(2 tau)) exp(-t/tau) (2k - 1/tau); 1/tau = 0.02

         k   2k-1/tau       err dt=4      err dt=2      err dt=1    err dt=0.5   observed order per step
     0.005      -0.01   5.497619e-02  2.709220e-02  1.344744e-02  6.699096e-03   1.021  1.011  1.005
      0.01          0   4.280997e-04  1.056540e-04  2.624076e-05  6.538511e-06   2.019  2.009  2.005
      0.02       0.02   4.572895e-03  2.206425e-03  1.083277e-03  5.366647e-04   1.051  1.026  1.013
      0.04       0.06   2.146076e-03  1.033079e-03  5.069571e-04  2.511327e-04   1.055  1.027  1.013

  prediction 1 (order in [0.8,1.2] whenever 2k != 1/tau): True
  prediction 3 (order ~2 at the degenerate k = 1/(2 tau)):  True
  prediction 2 (err(dt=4) ~ |2k - 1/tau|): predicted ratio 2.00 for k=0.02 vs k=0.005, MEASURED 0.083 -> FAILED

  POST-HOC closed form including the exp(-k(T-s)) propagation (no fitted constant):
           k      predicted       measured    ratio   (dt = 4 s)
       0.005   5.339599e-02   5.497619e-02    0.971
        0.01   0.000000e+00   4.280997e-04    0.000
        0.02   4.253587e-03   4.572895e-03    0.930
        0.04   1.990556e-03   2.146076e-03    0.928
    worst |log ratio| = 0.075 (7.8% off)
```

Predictions 1 and 3 — the ones about **order**, which is what the case exists to
measure — hold: the order is 1.021/1.011/1.005 at `k = 0.005`, 1.051/1.026/1.013
at `k = 0.02`, 1.055/1.027/1.013 at `k = 0.04`, and ~2.0 only at the degenerate
point. **Prediction 2 FAILED**, by a factor of 24 in the wrong direction. It was
wrong because it ignored how a per-step error is *propagated* to `t = T` by the
reaction, `exp(-k(T-s))`: larger k damps the accumulated error far faster than it
raises the per-step term. The corrected closed form,

```
E(T) = (dt/2)(Cmax/tau)(2k - 1/tau) exp(-kT) (exp((k - 1/tau)T) - 1)/(k - 1/tau)
```

reproduces the three non-degenerate measured errors to within **7.8 %** with no
fitted constant, and predicts exactly 0 at `k = 0.01` (where the measured
4.28e-04 is the next-order term). That formula is **post-hoc** and is labelled
so; the failure of prediction 2 as written is not retro-fitted away.

**Conclusion of E4/E4b:** the time integrator is first order, as documented.
The pre-registered E4 nonetheless failed, because a manufactured solution can be
accidentally degenerate with respect to the scheme under test, and a
single-parameter temporal MMS with no sensitivity sweep will not notice. That is
a defect in how the test was designed, and the record of it stays.

---

## 7. Plant tests

Each plant goes into `tmesim/physics/diffusion.py`, **first as a bare `raise` to
prove the line is reached**, then as the defect, then reverted. After every
revert `shasum -a 256 tmesim/physics/diffusion.py` returned
`233b9664272953c4f5b94594c4842ad1c9f296f5e42e83cfb1ff4a4f3e3f50ea`, the
pre-plant hash, and `git status` showed the file unmodified.

### P1 — no-flux mirrored about the cell CENTRE instead of the cell FACE

A ghost cell `C_{-1} = C_{+1}` instead of `C_{-1} = C_0`. This is the classic
"second order in the interior, first order at the boundary" defect the boundary
norm exists to catch.

Reach:
```
  File ".../tmesim/physics/diffusion.py", line 126, in _build_diffusion_operator
    raise RuntimeError("PLANT P1 REACHED")
RuntimeError: PLANT P1 REACHED
```

Red (E1), with E5 run alongside:
```
      dx     n      L2(all)    Linf(all)    Linf(bnd)    Linf(int)  resid/err     min S     min C    CFL      s
       4    50  1.07203e-02  5.75538e-02  5.75538e-02  4.35681e-02   6.17e-14    0.4141     1.033  0.000    0.0
       2   100  5.14623e-03  2.94983e-02  2.94983e-02  2.57647e-02   1.35e-13    0.4094     1.016  0.000    0.1
       1   200  2.51449e-03  1.48946e-02  1.48946e-02  1.39336e-02   2.68e-13    0.4082     1.008  0.000    0.7

  observed order per refinement step (r = 2), pre-registered band [1.8, 2.2]:
    L2   (global)    4->2:  1.059  2->1:  1.033   [FAIL]
    Linf (global)    4->2:  0.964  2->1:  0.986   [FAIL]
    Linf (boundary)  4->2:  0.964  2->1:  0.986   [FAIL]
    Linf (interior)  4->2:  0.758  2->1:  0.887   [FAIL]
  E1 VERDICT: FAILED
  ...
  E5 VERDICT: PASS
```

Order 2 -> 1 and the error up by two orders of magnitude. **E5 stayed PASS with
numbers identical to the clean run in every printed digit** (8.24985e-04 / 2.06217e-04 /
5.15525e-05) — P1 touches only the *domain-edge* rows, and E5's no-flux boundary
is the *tissue mask*, which is interior to the domain. E1 and E5 are therefore
verifying two genuinely different boundary implementations, and the pair of runs
demonstrates it rather than asserting it.

### P2 — harmonic face diffusivity replaced by the one-sided (donor) value

Reach:
```
  File ".../tmesim/physics/diffusion.py", line 111, in _harmonic
    raise RuntimeError("PLANT P2 REACHED")
RuntimeError: PLANT P2 REACHED
```

Red (E2), with E1 run alongside:
```
E1 ...
       4    50  2.44710e-04  4.87491e-04  4.87491e-04  4.87491e-04   6.38e-12
    Linf (global)    4->2:  1.996  2->1:  1.999   [PASS]
  E1 VERDICT: PASS

E2 ...
       4    50  5.10639e-04  1.19273e-03  1.19175e-03  1.19273e-03   2.23e-12
       2   100  2.38378e-04  5.15886e-04  5.14269e-04  5.15886e-04   6.89e-12
       1   200  1.16976e-04  2.38143e-04  2.37131e-04  2.38143e-04   1.86e-11
    L2   (global)    4->2:  1.099  2->1:  1.027   [FAIL]
    Linf (global)    4->2:  1.209  2->1:  1.115   [FAIL]
    Linf (boundary)  4->2:  1.212  2->1:  1.117   [FAIL]
    Linf (interior)  4->2:  1.209  2->1:  1.115   [FAIL]
  E2 VERDICT: FAILED
```

**This is the point of running both.** E1's numbers match the clean run in every
printed digit and E1 stays green, because with constant D the harmonic mean and the
one-sided value are the same number. A constant-D MMS is *green while not
testing its subject*. E2 is what tests the averaging, and it goes red.

### P3 — the advective divergence dropped entirely

Reach:
```
  File ".../tmesim/physics/diffusion.py", line 164, in _advective_divergence
    raise RuntimeError("PLANT P3 REACHED")
RuntimeError: PLANT P3 REACHED
```

Red (E3):
```
       4    50  7.07635e-02  1.82462e-01  1.82462e-01  1.78112e-01   5.11e-12     0.093    0.9674  0.100
       2   100  7.07618e-02  1.82830e-01  1.82830e-01  1.81727e-01   5.10e-12   0.09286     0.966  0.200
       1   200  7.07614e-02  1.82923e-01  1.82923e-01  1.82646e-01   5.12e-12   0.09281    0.9656  0.400

    L2   (global)    4->2:  0.000  2->1:  0.000   [FAIL]
    Linf (global)    4->2: -0.003  2->1: -0.001   [FAIL]
  E3 VERDICT: FAILED
```

The error stops falling entirely (order 0.000): the solver is now converging to
a different PDE. E3 genuinely exercises the advection term.

### P4 — the source contribution scaled by 0.99

Reach:
```
  File ".../tmesim/physics/diffusion.py", line 450, in simulate_drug_transport
    raise RuntimeError("PLANT P4 REACHED")
RuntimeError: PLANT P4 REACHED
```

Red (E4):
```
      dt    |C(T)-u(T)|      min S      s
       4   9.941003e-02     0.1002    0.0
       2   9.972925e-02     0.1002    0.0
       1   9.980787e-02     0.1002    0.0
     0.5   9.982737e-02     0.1002    0.0

    |C(T)-u(T)|      4->2: -0.005  2->1: -0.001  1->0.5: -0.000   [FAIL]
  E4 VERDICT: FAILED
```

The error stagnates at 9.94e-02 — exactly 1 % of the solution value ~9.98, the
planted deficit — and the order collapses to 0. Note that E4 catches this
*despite* having failed its own pre-registered band on the clean tree: a
degenerate manufactured solution weakened the case, it did not empty it.

### P5 — a plant predicted to stay GREEN, and it did

Pre-registered in `PREREG_E4b.md` at 2026-09-20 03:17:37, **before** running:
flipping the donor cell to the acceptor cell (`c_up = where(vf>0, c[:,1:],
c[:,:-1])` — i.e. **downwind**, the flux taken from the wrong side) is still a
first-order-accurate approximation of the face value, so E3 was predicted to
stay green. It did:

```
       4    50  2.80582e-03  7.56523e-03  4.40406e-03  7.56523e-03   5.70e-10
       2   100  1.39322e-03  3.75568e-03  2.19392e-03  3.75568e-03   1.14e-09
       1   200  6.94089e-04  1.87156e-03  1.09667e-03  1.87156e-03   2.27e-09
     0.5   400  3.46402e-04  9.34181e-04  5.48531e-04  9.34181e-04   4.54e-09

    L2   (global)    4->2:  1.010  2->1:  1.005  1->0.5:  1.003   [PASS]
    Linf (global)    4->2:  1.010  2->1:  1.005  1->0.5:  1.002   [PASS]
  E3 VERDICT: PASS
```

Order 1.010/1.005/1.003 and the error moved by 3 %. **Order verification cannot
see the upwind direction.** This is recorded as a measured limit of the method,
not a hypothetical one. (The direction is what makes the scheme *stable* at high
cell Peclet, and stability is not what an order study measures.)

### Regression tests

Two fast MMS guards were added to `tests/test_physics.py`
(`test_mms_diffusion_is_second_order`, `test_mms_upwind_advection_is_first_order`,
0.6 s combined). They were plant-tested on the same principle:

```
# under P1 (x edges only):
E  AssertionError: orders [0.92988237955428, 0.9736524353741294] from errors
   [0.055534928276474105, 0.02915034684858564, 0.014843801505180831]
1 failed, 1 passed

# under P3:
E  AssertionError: orders [-0.011621221462213348, -0.0029076601982810325] from errors
   [0.18099836690494442, 0.18246223598334055, 0.18283034785630026]
1 failed, 1 passed
```

Each plant reddens exactly the test whose subject it breaks and leaves the other
green. On the reverted tree, `tests/test_physics.py` is 21 passed.

---

## 8. What this does NOT show

**MMS verifies the SOLVER, not the model.** It answers exactly one question:
does the code solve the equation the module says it solves? It says nothing
whatever about whether that equation describes drug transport in tumour tissue.
Spelled out:

* **It is not validation.** Nothing here was compared to any measurement of any
  tissue. There is no anchor, no patient, no compound in this study. A solver
  can converge at its formal order to the solution of an equation that is the
  wrong equation for the physics.
* **It does not touch any parameter value.** `D`, the 0.6 tortuosity
  coefficient, the clearance and binding rates, the permeability, the
  BBB-penetrance term — none of them is tested, and the values used here
  (`D = 2` and `D = 300 um^2/s`, `k` up to 2 /s) were chosen to make the
  *numerics* clean, not to be physiological. MMS is indifferent to whether a
  coefficient is right.
* **It does not revive, support or bear on any published `results/` claim.** In
  particular it says nothing about the withdrawn transport-capacity anchor, the
  GMFE floor, cross-compound ranking, kill fraction, or any accuracy number in
  this repository. A verified solver and a wrong answer are entirely compatible.
* **It does not verify order for the QoIs.** The convergence measured here is of
  the concentration field `C` against a known `u`. `penetration_index`, `auc`
  and every downstream kill statistic are *functionals* of `C` computed by other
  code (percentiles, clipping, thresholds), and several are non-smooth. Their
  observed orders in `results/grid_convergence_gbm.md` (0.96, −3.33, both flagged
  degenerate) are not contradicted by, nor explained by, anything here.
* **It does not test the upwind direction** — plant P5 above, measured.
* **It does not test stability limits, the CFL guard, the mass-conservation
  identity, the non-negativity clamp, `exchange_k`, `capture_frames`, or the
  3-D case** (there is no 3-D case; the solver is 2-D by design).
* **It does not test the AUC time integral.** Only `final_concentration` was
  compared. The trapezoidal AUC accumulation is untested by this study.
* **It tests an unrealistically smooth problem.** The manufactured solutions are
  `C^infinity` cosines and the coefficients are smooth. Real runs have a binary
  tissue mask with a ragged perimeter, an ECM field from image data, and point-
  like vascular sources. Second order on a smooth problem is a *necessary*
  condition for a correct implementation, not a claim about the order achieved
  on a real section — where a non-smooth `D` or a stair-stepped mask will
  generally degrade it.
* **A single-mode manufactured solution is a weak boundary probe on its own.**
  In E1 and E5 the error field is essentially one cosine mode, so its Linf is
  attained at the corners *and* at interior extrema, and boundary and interior
  Linf coincide to three decimals. What makes the boundary split informative
  here is plant P1 moving it (1.996 -> 0.964), not the coincidence of the two
  numbers in the clean run.
* **One of the five cases failed on the clean tree** (E4), for a reason internal
  to the test design. The four spatial PASSes were not re-run, re-tuned or
  re-registered after that; they are as first measured.

---

## 9. Reproduce

```
cd ~/tme-sim
.venv/bin/python experiments/mms_transport_verification.py            # all cases, ~30 s
.venv/bin/python experiments/mms_transport_verification.py --case E3  # one case
.venv/bin/python experiments/mms_transport_verification.py --quick    # drop the finest rung
.venv/bin/python -m pytest tests/test_physics.py -q -k mms            # 0.6 s regression guards
```

Exit code is 0 only if every requested case passes its pre-registered band. The
script asserts `calibrate_alpha() == 22.223477711172563` before it measures
anything and prints it again at the end.
