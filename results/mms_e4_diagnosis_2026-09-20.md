# Why MMS case E4 failed — a closed-form diagnosis

**Date: 2026-09-20.** Driver: `experiments/mms_e4_diagnosis.py`. Raw run
reproduced verbatim in §6. `calibrate_alpha() = 22.223477711172563` asserted at
the start and the end of the run.

**This file does not modify `results/mms_transport_verification.md` and does not
revise E4's verdict. E4 remains FAILED.** No engine file was edited, no
tolerance was moved, no parameter was tuned.

---

## 1. Headline

**E4 failed because the manufactured solution it registered is exactly
degenerate for the scheme under test. The shipped time integrator has no
defect: it is first order, as documented.**

The evidence is not an asymptotic argument and contains no fitted constant. At
E4's configuration the shipped solver reduces *exactly* to a scalar recurrence
whose solution is known in closed form. That closed form reproduces E4's
measured errors to 2.3e-08 relative — every digit of the 2.019 / 2.009 / 2.005
that failed the test — and its O(dt) coefficient carries the factor
`(2k - 1/tau)`, which is **zero** at E4's registered `k = 0.01`, `tau = 50 s`.

The three candidate explanations were separated as follows.

| Candidate | Status | Evidence that settles it |
|---|---|---|
| The solver has a real defect (the time integrator is not the documented scheme) | **RULED OUT** | §3A: solver output matches the backward-Euler-with-lagged-source recurrence to 1.5e-14 relative at every dt. §3F: at E4's *own* k, a non-degenerate manufactured solution measures order 0.940 / 0.971 / 0.986. |
| The manufactured solution or its forcing is wrong | **HALF TRUE, and precisely which half is now known** | The forcing `S = u' + k u` is CORRECT (§3C: matches a central difference to 5.9e-10; the error falls 7.6e+12x over 20 halvings of dt with no inconsistency floor). The *choice of (u, k)* is the defect: §3D. |
| The case sits outside the regime the scheme claims | **RULED OUT** | §3G: E4 runs `convection=False` with no velocity, so the CFL branch is never entered and the reported CFL is `None`. Diffusion+reaction is implicit and unconditionally stable; E4 is stable and bounded even at dt = 320 s (one step). There is no regime boundary for E4 to be outside of. |

---

## 2. What E4 actually is

From `case_E4` in `experiments/mms_transport_verification.py`:

```
grid      16 x 16, dx = 10 um, ecm = 0, tissue = all, convection = False
drug      diffusivity 300 um^2/s, clearance k = 0.01 1/s, binding = 0
exact     u(t) = Cmax (1 - exp(-t/tau)),  Cmax = 10, tau = 50 s,  u(0) = 0
source    S(t) = u'(t) + k u(t), passed as source = ones * source_profile[n],
          read at the OLD time level t_n = n dt
window    T = 320 s;   ladder dt = 4, 2, 1, 0.5 s
metric    |C(T) - u(T)| at the single (uniform) value
registered  PASS iff every consecutive observed order lies in [0.8, 1.2]
measured    2.019 / 2.009 / 2.005   ->  FAILED
```

`u` is spatially uniform, so the discrete diffusion operator annihilates it
(`L` has zero row sums) and the case isolates the time integrator and the
source convention. That part of the design is sound — the isolation works.

---

## 3. The diagnosis

### A. The E4 configuration reduces the solver EXACTLY to a scalar recurrence

With `L u = 0`, `A = I - dt L + dt K` acting on a uniform field is just
multiplication by `(1 + k dt)`, so one solver step is

```
C^{n+1} = (C^n + dt S(t_n)) / (1 + k dt),        C^0 = 0
```

Measured against the shipped `simulate_drug_transport`:

```
     dt   steps              solver C(T)              closed form    rel diff
      4      80       9.9838125269792251       9.9838125269790865   1.390e-14
      2     160       9.9834900812510163        9.983490081250892   1.248e-14
      1     320       9.9834106680271333       9.9834106680269983   1.358e-14
    0.5     640       9.9833909657794919       9.9833909657793409   1.519e-14
```

This is the step that rules out a solver defect, and it has to be done *before*
any truncation argument. A defect that turned the scheme into Crank–Nicolson,
or that evaluated the source at the new time level, would also have produced
order 2 at E4's parameters. Only the exact match above distinguishes
"degenerate manufactured solution" from "the integrator is not what the module
says it is". It matches. The integrator is what the module says it is.

### B. The failing numbers are reproduced analytically

That recurrence has a closed-form solution (geometric sums; verified against a
direct 80-digit iteration to 8.0e-80 relative):

```
C^n = Cmax (1 - a^n) + dt * beta * a * (r^n - a^n)/(r - a)
a = 1/(1 + k dt),  r = exp(-dt/tau),  beta = Cmax (1/tau - k),  n = T/dt
E(T) = C^n - u(T) = (r^n - a^n) * ( Cmax + dt * beta * a / (r - a) )
```

Evaluated at 80 decimal digits against the solver's measured error:

```
     dt   measured |C-u|  closed form |E|     ratio
      4     4.280997e-04     4.280997e-04   1.00000
      2     1.056540e-04     1.056540e-04   1.00000
      1     2.624076e-05     2.624076e-05   1.00000
    0.5     6.538511e-06     6.538511e-06   1.00000

 observed order, measured   :  2.019   2.009   2.005
 observed order, closed form:  2.019   2.009   2.005
 worst |ratio-1| = 2.316e-08
```

The failing result is accounted for to eight significant figures with no fitted
constant. (`results/mms_transport_verification.md` §6 reached the same
conclusion with a post-hoc formula good to 7.8 %; the closed form above is
exact and needs no error budget.)

### C. The forcing term is correct

Two checks, because "the manufactured forcing is wrong" is the alternative that
makes a correct solver look broken:

* `S(t) = u'(t) + k u(t)` against a central difference of `u`: worst relative
  mismatch **4.2e-10** (E4's exponential) and **5.9e-10** (the rational
  solution of §3F).
* A wrong `S` makes `u` not a solution, so `C(T)` converges to the wrong limit
  and the error stops at a nonzero floor. It does not: over 20 halvings of dt
  (2 steps to 4 194 304 steps) the error falls from 1.155e+00 to 1.517e-13, a
  factor **7.6e+12**, at a clean order 2 throughout. There is no floor, so `u`
  *is* the exact solution of the semi-discrete problem and `S` is right.

### D. The O(dt) term is not small at E4's parameters — it is absent

Expanding `E(T)` in dt gives the exact first-order coefficient

```
G1 = -(Cmax/(2 tau)) * (2k - 1/tau)/(k - 1/tau) * (exp(-T/tau) - exp(-kT))
```

(cross-checked against Richardson extrapolation of the closed form at
dt = T/2^26 and T/2^27 — agreement to <1e-6 relative at every non-degenerate k).

```
        k   2k - 1/tau     G1 (closed)   G1 (Richardson)             G2
    0.005        -0.01    1.334900e-02      1.334900e-02   9.833573e-05
    0.008       -0.004    2.521439e-03      2.521439e-03   5.160108e-05
   0.0099      -0.0002    8.005175e-05      8.005175e-05   2.720619e-05
  0.00999       -2e-05    7.838421e-06      7.838421e-06   2.618031e-05
     0.01            0   -0.000000e+00     -1.979582e-18   2.606710e-05  <-- E4
  0.01001        2e-05   -7.801885e-06     -7.801885e-06   2.595405e-05
   0.0101       0.0002   -7.639779e-05     -7.639779e-05   2.494362e-05
    0.012        0.004   -9.916022e-04     -9.916022e-04   6.640109e-06
     0.04         0.06   -4.976390e-04     -4.976390e-04  -9.187820e-06
```

E4 registered `k = 0.01` with `tau = 50 s`, so `2k = 0.02 = 1/tau` **exactly**,
`G1 = 0`, and the leading surviving term is `G2 dt^2`. The measured order is 2
for the same reason a quadrature rule is exact on a polynomial it happens to
integrate exactly: not because the method is better than claimed, but because
this particular problem cannot see the method's leading error.

**The general condition.** The leading local truncation error of backward Euler
with a lagged source is `dt * (u''/2 + k u')`. It vanishes for all `t` iff
`u' ∝ exp(-2kt)`, i.e. iff `u ∈ span{1, exp(-2kt)}`. E4 chose
`u' ∝ exp(-t/tau)` and then set `1/tau = 2k`. It landed on that one-parameter
blind family exactly.

### E. The discriminator: order must cross over 2 → 1 as soon as k leaves the point

A genuinely second-order integrator shows order 2 at every dt and every k. A
degenerate manufactured solution shows `E = G1 dt + G2 dt^2`: order 2 while
`|G2| dt > |G1|`, order 1 below, changing over at `dt* = |G1/G2|` — a dt* that
shrinks in proportion to the detuning `2k - 1/tau`. Measured on the shipped
solver, with the closed form continuing the ladder past where the solver is
affordable:

| k | 2k − 1/τ | predicted dt\* = \|G1/G2\| | exact sign change of E(dt) | asymptotic order |
|---|---|---|---|---|
| 0.0105   | 1e-03 | 16.894 s | 14.884 s | **1.000 1.000 1.000** |
| 0.01005  | 1e-04 | 1.5154 s | 1.4988 s | **1.000 1.000 1.000** |
| 0.010005 | 1e-05 | 0.15015 s | 0.14999 s | **0.999 1.000 1.000** |
| 0.01 (E4) | 0 | — (G1 = 0) | never changes sign | **2.000 2.000 2.000** |

The predicted changeover and the exact sign change of `E(dt)` agree to 12 %,
1.1 % and 0.1 % respectively — the agreement improves as dt\* moves into the
asymptotic range, which is what a two-term model should do and a coincidence
should not. The order tables in §6 show the full rollover, and the
`G1dt+G2dt^2` column marks where the two-term model is and is not valid.

Only the single point `k = 1/(2 tau)` stays at order 2 all the way down to
dt = 1.9e-05 s (16 777 216 steps). The set of k for which this integrator looks
second order has measure zero, and E4 sits on it.

### F. Same solver, same k, a non-degenerate manufactured solution → order 1

The cleanest single demonstration. Keep E4's grid, its `k = 0.01`, its `T`,
its source convention. Change only the manufactured solution to one whose
derivative is rational and therefore can never be proportional to `exp(-2kt)`:

```
u(t) = Cmax t/(t + tau)      u(0) = 0 (the solver's fixed IC), u' > 0
S(t) = Cmax tau/(t+tau)^2 + k Cmax t/(t+tau)     min S = 8.75e-02 > 0 (V1 holds)

     dt    |C(T) - u(T)|
      4     9.649807e-03
      2     5.031077e-03
      1     2.565954e-03
    0.5     1.295438e-03
 observed order:  0.940   0.971   0.986      -> inside E4's registered [0.8, 1.2]
```

At exactly the parameters where E4 measured 2, a non-degenerate manufactured
solution measures 1. The order-2 result is a property of E4's `u`, not of the
code. This is the minimal repair to the *test*; it is recorded here and is not
retro-fitted into the pre-registration.

### G. CFL and regime are not implicated

E4 passes `convection=False` and no velocity, so `simulate_drug_transport` never
enters the explicit-upwind branch, never evaluates the CFL guard, and reports
`cfl = None`. Diffusion and reaction are implicit (`A = I - dt L + dt K`, an
M-matrix for any dt > 0), so there is no `dt <= dx^2/(4D)` restriction and
backward Euler is unconditionally stable — confirmed by running E4 at dt = 320,
160 and 80 s (1, 2 and 4 steps): the answer stays finite and bounded
(|C−u| = 5.25e+00, 1.16e+00, 2.40e-01, falling monotonically).

**Connection to the CFL finding of the same session.** The measured
2x-looseness of the CFL guard in `tmesim/physics/diffusion.py` (strict margin
0.6506 at the shipped dx = 12 um, dt = 10 s, deliberately not tightened because
tightening would refuse the product's own default configuration) concerns the
**explicit upwind advection term**. E4 does not exercise it. The two findings do
not interact numerically, and nothing here contradicts that decision.

They do rhyme, though, and the rhyme is the transferable lesson: in both cases
a check's *stated justification* was broader than the thing the check actually
establishes. The V1 guard's erratum in `results/mms_transport_verification.md`
says the same. E4's registration said "this measures the temporal order"; what
it measured was the temporal order *of one problem that cannot see the leading
temporal error*.

---

## 4. Does this imply a defect in the shipped solver?

**No.** The transport solver's time integrator is backward Euler on
diffusion+reaction with the source lagged to the old time level, verified here
to 1.5e-14 against its own closed form, and first order in time, verified at
E4's own parameters by a non-degenerate manufactured solution (0.940 / 0.971 /
0.986). Nothing in the engine was changed and nothing needs to be.

The defect is in the **test**, and it is a real one worth recording: a
single-parameter temporal MMS case with no sensitivity sweep cannot tell a
first-order integrator from a degenerate problem. Three consequences:

1. **E4's FAILED verdict is correct and stays.** The test failed its registered
   band. What it failed to do was measure what it was registered to measure.
2. **E4 retains its plant-test power but not its verification power.** The P4
   plant (source scaled by 0.99) still collapses E4's order to 0, as recorded in
   `results/mms_transport_verification.md` §7. What E4 at `k = 1/(2 tau)` cannot
   do is *confirm* the integrator's order — it reports 2 whether the scheme is
   backward Euler or something genuinely second order.
3. **The four passing cases are not at risk from the same failure mode.** A
   degenerate manufactured solution shows up as an order *higher* than claimed —
   that is the entire E4 signature. E1, E2 and E5 measured 2 where 2 was
   claimed, E3 measured 1 where 1 was claimed. None of them shows the
   signature, so none of them is degenerate. (For the record, the 5-point
   Laplacian's leading error `(dx^2/12)(u_xxxx + u_yyyy)` applied to
   `u = A + B cos(kx x) cos(ky y)` gives `B(kx^4 + ky^4) cos cos`, which is not
   identically zero for any admissible A, B, kx, ky — the steady cases have no
   degenerate family to land on.)

**A temporal MMS case should sweep at least one coefficient.** If E4 had run
its ladder at two values of k rather than one, the degeneracy would have been
visible on the first run instead of needing a closed-form solution to settle.

---

## 5. What this does NOT establish

* Nothing about tissue, any parameter value, or any published number in
  `results/`. This is verification, not validation.
* Nothing about whether the shipped default dt is *accurate enough* for any
  simulation. It establishes the *order*, which says how error scales with dt,
  not how large it is at production settings. That is a separate question and is
  not answered here.
* Nothing about the spatial discretisation beyond what
  `results/mms_transport_verification.md` already recorded.

---

## 6. Raw run, verbatim

`.venv/bin/python experiments/mms_e4_diagnosis.py`

```
calibrate_alpha() = 22.223477711172563

==============================================================================
A. Does the shipped solver reduce EXACTLY to backward Euler + lagged source?
==============================================================================
   E4 config: 16x16, ecm=0, convection=False, source=ones*S(t_n).
   The discrete operator L has zero row sums, so L*(uniform) = 0 and
   A*(uniform) = (1 + k dt)*(uniform). Claim: C^{n+1} = (C^n + dt S_n)/(1+k dt).

        dt   steps              solver C(T)              closed form    rel diff
         4      80       9.9838125269792251       9.9838125269790865   1.390e-14
         2     160       9.9834900812510163        9.983490081250892   1.248e-14
         1     320       9.9834106680271333       9.9834106680269983   1.358e-14
       0.5     640       9.9833909657794919       9.9833909657793409   1.519e-14

   worst relative difference = 1.519e-14   (MATCH, threshold 1e-13)
   closed form vs direct 80-digit iteration at dt=4: rel 8.013e-80
   => the solver's time integrator is EXACTLY the documented scheme.
      Any order anomaly is a property of this recurrence, not of code.

==============================================================================
B. Is E4's measured error reproduced by the closed form, with no fitted constant?
==============================================================================
   u(T) = 9.9833844272682608   (k=0.01, tau=50.0, Cmax=10.0, T=320.0)

        dt   measured |C-u|  closed form |E|     ratio
         4     4.280997e-04     4.280997e-04   1.00000
         2     1.056540e-04     1.056540e-04   1.00000
         1     2.624076e-05     2.624076e-05   1.00000
       0.5     6.538511e-06     6.538511e-06   1.00000

   observed order, measured   :  2.019   2.009   2.005
   observed order, closed form:  2.019   2.009   2.005
   worst |ratio-1| = 2.316e-08  (MATCH)
   => E4's order-2 result is entirely accounted for analytically.

==============================================================================
C. Is the MANUFACTURED FORCING wrong? (the alternative that makes a correct solver look broken)
==============================================================================
   A wrong S makes u not a solution of the semi-discrete problem, so
   C(T) converges to something OTHER than u(T): the error stops at a
   nonzero floor instead of going to 0. Two checks.

   C1. S(t) = u'(t) + k u(t) against a central difference of u:
       exponential (E4)     worst relative mismatch 4.178e-10
       rational (F)         worst relative mismatch 5.864e-10
       -> S is consistent with u (threshold 1e-8)

   C2. Does the E4 error go to ZERO, or to a floor? Closed-form ladder
       (validated against the solver in B) over 20 halvings of dt:
             dt (s)      steps         |E(T)|    order
                160          2   1.155034e+00         
                 20         16   1.174955e-02    2.206
                2.5        128   1.656220e-04    2.050
             0.3125       1024   2.550924e-06    2.007
           0.039062       8192   3.978561e-08    2.001
          0.0048828      65536   6.215083e-10    2.000
         0.00061035     524288   9.710790e-12    2.000
         7.6294e-05    4194304   1.517306e-13    2.000
       error fell by 7.612e+12x with no floor -> S is CORRECT
       (order column is per 8x refinement / log8, i.e. 2.0 == 6.0 here)

==============================================================================
D. The exact O(dt) coefficient G1 of the global error, swept over k
==============================================================================
   G1 = -(Cmax/(2 tau)) (2k - 1/tau)/(k - 1/tau) (e^{-T/tau} - e^{-kT})
   derived from the closed form; cross-checked against a Richardson
   extrapolation of the closed form at dt = T/2^26 and T/2^27.

   1/tau = 0.02;  E4 registered k = 0.01, so 2k - 1/tau = 0

           k   2k - 1/tau     G1 (closed)   G1 (Richardson)             G2
       0.005        -0.01    1.334900e-02      1.334900e-02   9.833573e-05
       0.008       -0.004    2.521439e-03      2.521439e-03   5.160108e-05
      0.0099      -0.0002    8.005175e-05      8.005175e-05   2.720619e-05
     0.00999       -2e-05    7.838421e-06      7.838421e-06   2.618031e-05
        0.01            0   -0.000000e+00     -1.979582e-18   2.606710e-05   <-- exactly zero
     0.01001        2e-05   -7.801885e-06     -7.801885e-06   2.595405e-05
      0.0101       0.0002   -7.639779e-05     -7.639779e-05   2.494362e-05
       0.012        0.004   -9.916022e-04     -9.916022e-04   6.640109e-06
        0.04         0.06   -4.976390e-04     -4.976390e-04  -9.187820e-06

   At E4's own k: G1 = -0.000000e+00 exactly, killed by the factor (2k - 1/tau).
   The O(dt) term is not small at E4's parameters -- it is ABSENT,
   so the leading surviving term is G2 dt^2 and the measured order is 2.

   Degeneracy condition, general form: the leading local truncation
   error of this scheme is dt*(u''/2 + k u'), which vanishes for all t
   iff u' is proportional to exp(-2kt). E4 chose u' ~ exp(-t/tau)
   with 1/tau = 2k, i.e. it landed on that one-parameter family exactly.

==============================================================================
E. Detune k off the degeneracy: the order must CROSS OVER from 2 to 1
==============================================================================
   If the integrator were genuinely second order, the order would be 2
   at every dt and every k. If instead G1 vanishes only AT k = 1/(2 tau),
   then for k nearby the error is G1 dt + G2 dt^2: order 2 while
   |G2| dt > |G1|, order 1 below, with the changeover at dt* = |G1/G2|
   (and, when G1 and G2 have opposite signs, a near-exact cancellation
   -- a visible error MINIMUM -- at that same dt*).

   k = 0.0105    2k-1/tau =     0.001   G1 = -3.48144e-04  G2 =  2.06077e-05
     predicted changeover dt* = |G1/G2| = 16.894 s  (signs oppose -> error minimum there)
     exact sign change of E(dt) (bisection on the closed form): dt = 14.884 s
          dt (s)     steps    |E| closed    |E| solver   G1dt+G2dt^2   order
              20        16   2.63024e-03   2.63024e-03   1.28021e-03     nan
              10        32   1.24759e-03   1.24759e-03   1.42067e-03   1.076
               5        64   1.20361e-03   1.20361e-03   1.22553e-03   0.052
             2.5       128   7.38805e-04   7.38805e-04   7.41562e-04   0.704
          0.3125      1024   1.06777e-04   1.06777e-04   1.06783e-04   0.972
        0.039062      8192   1.35679e-05   1.35679e-05   1.35679e-05   0.997
       0.0048828     65536   1.69943e-06            --   1.69943e-06   1.000
      0.00061035    524288   2.12483e-07            --   2.12483e-07   1.000
      0.00015259   2097152   5.31221e-08            --   5.31221e-08   1.000
      7.6294e-05   4194304   2.65612e-08            --   2.65612e-08   1.000
      3.8147e-05   8388608   1.32806e-08            --   1.32806e-08   1.000
      1.9073e-05  16777216   6.64032e-09            --   6.64032e-09   1.000
     asymptotic order (finest 3 intervals): 1.000 1.000 1.000 -> FIRST ORDER

   k = 0.01005   2k-1/tau =    0.0001   G1 = -3.86469e-05  G2 =  2.55034e-05
     predicted changeover dt* = |G1/G2| = 1.5154 s  (signs oppose -> error minimum there)
     exact sign change of E(dt) (bisection on the closed form): dt = 1.4988 s
          dt (s)     steps    |E| closed    |E| solver   G1dt+G2dt^2   order
              20        16   1.07544e-02   1.07544e-02   9.42842e-03     nan
              10        32   2.33390e-03   2.33390e-03   2.16387e-03   2.204
               5        64   4.65883e-04   4.65883e-04   4.44351e-04   2.325
             2.5       128   6.54884e-05   6.54884e-05   6.27791e-05   2.831
          0.3125      1024   9.58126e-06   9.58126e-06   9.58658e-06   0.562
        0.039062      8192   1.47072e-06   1.47072e-06   1.47073e-06   0.961
       0.0048828     65536   1.88097e-07            --   1.88097e-07   0.995
      0.00061035    524288   2.35787e-08            --   2.35787e-08   0.999
      0.00015259   2097152   5.89645e-09            --   5.89645e-09   1.000
      7.6294e-05   4194304   2.94837e-09            --   2.94837e-09   1.000
      3.8147e-05   8388608   1.47422e-09            --   1.47422e-09   1.000
      1.9073e-05  16777216   7.37121e-10            --   7.37121e-10   1.000
     asymptotic order (finest 3 intervals): 1.000 1.000 1.000 -> FIRST ORDER

   k = 0.010005  2k-1/tau =     1e-05   G1 = -3.90550e-06  G2 =  2.60106e-05
     predicted changeover dt* = |G1/G2| = 0.15015 s  (signs oppose -> error minimum there)
     exact sign change of E(dt) (bisection on the closed form): dt = 0.14999 s
          dt (s)     steps    |E| closed    |E| solver   G1dt+G2dt^2   order
              20        16   1.16492e-02   1.16492e-02   1.03261e-02     nan
              10        32   2.73165e-03   2.73165e-03   2.56200e-03   2.092
               5        64   6.52221e-04   6.52221e-04   6.30736e-04   2.066
             2.5       128   1.55506e-04   1.55506e-04   1.52802e-04   2.068
          0.3125      1024   1.32493e-06   1.32493e-06   1.31962e-06   2.550
        0.039062      8192   1.12859e-07   1.12860e-07   1.12870e-07   0.374
       0.0048828     65536   1.84497e-08            --   1.84497e-08   0.951
      0.00061035    524288   2.37404e-09            --   2.37404e-09   0.994
      0.00015259   2097152   5.95327e-10            --   5.95327e-10   0.999
      7.6294e-05   4194304   2.97815e-10            --   2.97815e-10   0.999
      3.8147e-05   8388608   1.48945e-10            --   1.48945e-10   1.000
      1.9073e-05  16777216   7.44821e-11            --   7.44821e-11   1.000
     asymptotic order (finest 3 intervals): 0.999 1.000 1.000 -> FIRST ORDER

   k = 0.01      2k-1/tau =         0   G1 = -1.97958e-18  G2 =  2.60671e-05
     G1 = 0 exactly -> no crossover is possible; order must be 2 at EVERY dt
     E(dt) never changes sign on [1e-4, 40] s
          dt (s)     steps    |E| closed    |E| solver   G1dt+G2dt^2   order
              20        16   1.17496e-02   1.17496e-02   1.04268e-02     nan
              10        32   2.77631e-03   2.77631e-03   2.60671e-03   2.081
               5        64   6.73156e-04   6.73156e-04   6.51677e-04   2.044
             2.5       128   1.65622e-04   1.65622e-04   1.62919e-04   2.023
          0.3125      1024   2.55092e-06   2.55092e-06   2.54562e-06   2.003
        0.039062      8192   3.97856e-08   3.97832e-08   3.97752e-08   2.000
       0.0048828     65536   6.21508e-10            --   6.21488e-10   2.000
      0.00061035    524288   9.71079e-12            --   9.71075e-12   2.000
      0.00015259   2097152   6.06923e-13            --   6.06922e-13   2.000
      7.6294e-05   4194304   1.51731e-13            --   1.51730e-13   2.000
      3.8147e-05   8388608   3.79326e-14            --   3.79326e-14   2.000
      1.9073e-05  16777216   9.48316e-15            --   9.48316e-15   2.000
     asymptotic order (finest 3 intervals): 2.000 2.000 2.000 -> stays 2, as required

   Every non-degenerate k lands on order 1 once dt < dt*; only the
   single point k = 1/(2 tau) stays at 2 forever. That is a property of
   the manufactured solution, not of the integrator.

==============================================================================
F. The same solver, same k, same T, a NON-degenerate manufactured solution
==============================================================================
   u(t) = Cmax t/(t+tau), u(0)=0 (the solver's fixed IC), u' > 0,
   u' is rational, so it can never be proportional to exp(-2kt):
   no choice of k can make this one degenerate.

   u(T) = 8.6486486486486491   (Cmax=10.0, tau=50.0, T=320.0, k=0.01 -- E4's own k)

        dt    |C(T) - u(T)|        min S
         4     9.649807e-03   8.7501e-02
         2     5.031077e-03   8.7500e-02
         1     2.565954e-03   8.7500e-02
       0.5     1.295438e-03   8.7500e-02

   observed order:  0.940   0.971   0.986
   V1 guard min(S) > 0 over every level: 8.7500e-02 -> OK
   order in the [0.8, 1.2] band E4 registered: YES -- first order
   => at E4's OWN k the integrator is first order. The order-2 result
      belongs to E4's manufactured solution, not to the solver.

==============================================================================
G. Is E4 outside a stated regime (CFL, stability)?
==============================================================================
   E4 runs convection=False and passes no velocity, so the explicit
   upwind branch -- and with it the CFL guard -- is never entered.
   Reported CFL for the E4 run: None
   Diffusion and reaction are solved implicitly (A = I - dt L + dt K,
   an M-matrix for any dt > 0), so there is no dt <= dx^2/(4D) limit
   and backward Euler is unconditionally stable. A dt ladder of
   4/2/1/0.5 s on a 16x16 grid with u spatially uniform sits in the
   scheme's stated regime at every level.
   dt =  320.0 s (  1 step(s)): C(T) = 15.238095, |C-u| = 5.2547e+00  (finite, bounded -> stable)
   dt =  160.0 s (  2 step(s)): C(T) = 11.138418, |C-u| = 1.1550e+00  (finite, bounded -> stable)
   dt =   80.0 s (  4 step(s)): C(T) = 10.223654, |C-u| = 2.4027e-01  (finite, bounded -> stable)
   => 'E4 is outside the regime' is not available as an explanation:
      the only regime boundary this scheme has (the advective CFL) is
      not exercised by E4 at all. The separately measured 2x-loose CFL
      guard (strict margin 0.6506 at the shipped dx=12 um, dt=10 s)
      concerns explicit upwind advection and has no bearing here.

==============================================================================
SUMMARY
==============================================================================
  OK    A solver == backward Euler recurrence
  OK    B measured error == closed form
  OK    C manufactured forcing S is correct
  OK    D G1 vanishes exactly iff 2k = 1/tau
  OK    E order crosses 2 -> 1 off the degeneracy
  OK    F non-degenerate u at E4's own k is order 1
  OK    G CFL / regime not implicated

  elapsed 2.0 s
calibrate_alpha() = 22.223477711172563
```
