"""THE SCALE-FREE REACH MAP — how much of the tumour the drug never reaches.

Why this module exists
----------------------
Vivum has stopped ranking compounds. The cross-compound ordering was measured
and does not hold (``results/pk_ranking_robustness.md``,
``results/capacity_ranking_robustness.md``), so the product is the SINGLE-compound
map: one real tissue, one compound, where the drug reaches and how much of the
tumour it never reaches.

The replacement headline therefore has to be **SCALE-FREE**, and that constraint
is the entire design.
``validation.delivery_metrics.therapeutic_coverage`` cannot be it. Its own
docstring lists "(d) **absolute** — scaling the field scales the answer" as a
property, and the two quantities it reads — the IC50 and the plasma
concentration — are exactly the two the robustness studies measured as uncertain
(per-compound plasma exposure across the retained GBM compounds spans
**113.1x**, ``results/pk_ranking_robustness.md`` sec. 0). Shipping an absolute
number as the headline would rebuild the same failure one compound at a time.

WHICH OBJECT: the map is a STEADY-STATE object (2026-08-22)
-----------------------------------------------------------
The claim this module makes changed from *"the reach this dose achieved in 6
hours"* to **"the reach this tissue's geometry permits"** — vessel spacing and
the compound's transport, with no dosing schedule entangled. Three measurements
pointed here and a fourth, made while implementing, changed one of their
reasons. All four are in ``results/steady_state_reach.md``:

1. The **numerical route is a single sparse solve, not a longer window.** The
   shipped stepper advances ``(I - dt*L + dt*K) c_new = c + dt*S - dt*Adv(c)``,
   so its fixed point is the LINEAR SYSTEM ``(K - L + Adv) c = S`` with every
   ``dt`` cancelled. :func:`steady_state_concentration` solves it once, using
   the engine's own ``L``. There is no window and no ``dt`` in the answer, so
   there is no convergence criterion to defend.
2. Obtaining the field that way is **cheaper than marching to it**, which
   inverts the obvious prior: one factorisation and one triangular solve
   against one factorisation and 2,160 triangular solves — 0.90 s against
   43.70 s at the shipping grid. **This is not a performance win for the
   product and must not be quoted as one.** The run still marches
   (``simulate_response`` is required for ``kill_fraction``, for
   ``delivery_metrics`` and for the oxygen field this map's own viable band
   reads), so the steady solve is computed IN ADDITION and *adds* 1.854 s,
   **+4.2 %**, to a run at dx = 12. What the cost measurement decided was
   whether the route was shippable at all; it is, with room to spare
   (``results/steady_state_reach.md`` sec. 2.2).
3. A capacity term ``R`` multiplies the TIME DERIVATIVE and nothing else, so it
   is **absent from the steady problem by inspection**. The finite-window
   residual it removes is large and measured.
4. **The previous disclosure's evidence was about the estimator, not the
   field.** It reported the shipping window as unconverged, citing ``281/T``
   with no plateau. ``AUC(x,T) = C_ss(x)*T + J(x,T)`` with ``J`` bounded, so an
   AUC-based ratio ALWAYS has an ``O(1/T)`` tail however converged the field is.
   The field at 6 h is converged to ``1.112e-05`` of its own maximum. The
   capacity residual is real anyway, because ``R`` shrinks the EFFECTIVE window
   (``AUC_R(x,T) = R*AUC_1(x, T/R)``): at ``R = 39`` the 6 h window becomes
   554 s. Capacity manufactures an unconverged window rather than revealing one.

The finite-window object is still reachable (omit ``drug``) because the
production caller cannot yet supply the compiled ``DrugParams``. Every payload
says which one it is in ``window_mode``, and carries the disclosure that matches.

What IS protected is a **RATIO BETWEEN TWO POINTS OF ONE MAP**, and only to the
extent stated here. Two perturbations were swept and measured:

* **Plasma / dose scale — EXACT, for the finite-window AUC field.** The delivery
  source is one multiplication by the plasma level and the transport PDE is
  linear in ``C``, so ``AUC`` is exactly proportional to it. Measured at the
  shipping configuration (real Xenium GBM section, dx = 12 um, dt = 10 s,
  convection ON, 6 h window) over plasma factors 0.005 to 26 — a 5,200x span —
  the worst relative deviation of the AUC field from an exact rescale anywhere
  was **1.41e-14** (``results/pk_ranking_robustness.md`` sec. 1.2). A factor that
  multiplies every point of the map divides out of a ratio between points.
* **Tissue capacity — EXACT for a CONVERGED gradient shape, NOT for a finite
  window.** ``results/capacity_ranking_robustness.md`` sec. 2b states it at that
  strength and no more: *"steady-state and equilibrium-referenced outputs —
  penetration length, interstitial:plasma ratio, converged gradient shape — are
  EXACTLY invariant to tissue capacity (verified to 1e-9). The finite-window kill
  ranking is not."* A capacity term ``R`` enters as ``dC/dt = (1/R)[...]``, i.e.
  a pure rescale of TIME, so ``C_R(x, t) = C_1(x, t/R)`` and over a finite window
  ``AUC_R(x, T) = R * AUC_1(x, T/R)``. The spatial shape is then the shape at a
  DIFFERENT effective window: the between-point ratio is exactly
  capacity-invariant only where the field has converged, and carries a residual
  where it has not. This product's window is finite (6 h by default), so the
  reach ratio sits in the "converged gradient shape" family only to the extent
  the field has converged, and that is not asserted here.

**Neither invariance protects against error in the transport law itself**, whose
descriptor path is unscored (``results/descriptor_falsification.md``). And the
two bullets above are a LIST OF THE AXES THAT WERE SWEPT, not a coverage claim:
silence about an unswept axis here is silence, not a pass.

The quantities (exact definitions)
----------------------------------
Let ``AUC(x)`` be the transport solve's exposure field (uM*s).

* **Perivascular reference exposure** ``E_ref`` (uM*s) — the mean AUC over tissue
  in the FIRST usable bin of the distance-from-nearest-endothelial-cell profile.
  Bins are ``vessel_profile.BIN_UM`` = 10 um wide and a bin is usable only with at
  least ``vessel_profile.MIN_BIN_COUNT`` = 30 grid points in it. This is ALREADY
  computed as ``delivery_metrics.penetration_length_fit(...)["drug"]
  ["perivascular_mean"]`` and that value is REUSED here verbatim — never
  recomputed — so the reach number and the published penetration length can
  never drift apart.
* **Reach ratio** ``r_i = AUC_i / E_ref`` for tumour cell ``i``. Dimensionless.
  Numerator and denominator are both linear in any plasma scale factor, so it is
  invariant to dose; both are read off the same field at the same window, so a
  capacity term cancels from it exactly wherever that field has converged (see
  above for where it does not).
* **Under-reached fraction** ``F(x)`` — the fraction of VIABLE tumour cells with
  ``r_i < x``. "Viable" is defined EXACTLY as
  ``delivery_metrics.hypoxic_vs_normoxic_exposure_ratio`` defines it — local
  oxygen ``o >= 2 * OxygenConfig.hypoxia_threshold`` at the cell's own grid
  point, the engine's own banding (``tmesim.cells.agents.update_cell_states``) —
  so the two cannot disagree about which cells they are talking about.
* **The curve** — ``F(x)`` over the fixed log ladder :data:`REACH_LADDER`.
* **The headline** — ``F(0.1)``, picked off that ladder, not recomputed.
* **The scale-free falloff** — the same distance bins that produced ``E_ref``,
  each bin's mean AUC divided by ``E_ref``: exposure vs distance from the nearest
  vessel, normalised, with the bin distances and per-bin counts.

**The 0.1 threshold is a REPORTING CONVENTION, not a physical constant.** It is
returned in the payload as :data:`CONVENTION_NOTE` so a reader meets the sentence
where they meet the number, and the whole ladder is returned beside the headline
so the choice stays visible. That is deliberate: ``1x IC50`` acquired an
authority it had not earned by being quoted alone, and ``1/10`` must not repeat
it.

REFUSAL, not defaulting
-----------------------
Whenever a quantity is undefined — no tumour cells, no viable tumour cells, no
usable distance bin, ``E_ref`` non-finite or non-positive — this returns
``null_reason`` and **no number**. Nothing is substituted, extrapolated or
imputed. This is the repo's standing convention (see
``vessel_profile.D_HALF_NEVER`` for the same rule applied to a half-distance
that does not exist).

Nothing in ``tmesim/`` is modified by this module; it only consumes a
``Microenvironment``, a ``SimulationResult`` and the already-computed
``penetration_length_fit`` block.

Units: distance um, reach ratio and F(x) dimensionless in both modes. ``E_ref``
and the field it is taken from are **uM*s in finite-window mode and uM in
steady-state mode** — an exposure against a concentration. The wire key stays
``E_ref_uM_s`` because ``app/web`` and ``app/server/report.py`` read that name,
so the payload carries the true unit in ``E_ref_units`` and in ``units['E_ref']``
and a consumer must read it there rather than from the key.
"""
from __future__ import annotations

from typing import Dict, Optional, Sequence

import numpy as np
import scipy.sparse as _sp
from scipy.sparse.linalg import factorized as _factorized

from tmesim.core import DrugParams, Microenvironment, OxygenConfig, SimulationResult
# Imported, never re-derived: the steady-state solve below must be the fixed
# point of THE SHIPPED time-stepper, so it assembles its operator from the same
# two functions that time-stepper assembles its own from.
from tmesim.physics.diffusion import (
    _build_diffusion_operator,
    _effective_diffusivity,
)
from tmesim.validation.delivery_metrics import (
    penetration_length_fit,
    tumor_cell_exposure,
)
# Imported deliberately rather than re-implemented: this is the same
# nearest-grid-cell sampling map (``domain.to_index``) that
# ``hypoxic_vs_normoxic_exposure_ratio`` uses to read a tumour cell's local
# oxygen, and the "viable" band here has to be the identical set of cells.
from tmesim.validation.delivery_metrics import _sample_at_cells
from tmesim.validation.vessel_profile import (
    BIN_UM,
    MIN_BIN_COUNT,
    distance_profile,
    usable_bins,
)

__all__ = [
    "REACH_LADDER",
    "HEADLINE_X",
    "CONVENTION_NOTE",
    "SCALE_INVARIANCE_NOTE",
    "STEADY_SCALE_INVARIANCE_NOTE",
    "STEADY_STATE_NOTE",
    "FINITE_WINDOW_NOTE",
    "WINDOW_MODES",
    "compute_reach",
    "null_reach",
    "steady_state_concentration",
]

#: The fixed log ladder of reach ratios at which ``F(x)`` is reported. Returned
#: in full, always, so the headline threshold is visibly one point on a curve
#: rather than a constant of nature.
REACH_LADDER: Sequence[float] = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)

#: The headline point of the ladder. It is a member of :data:`REACH_LADDER` and
#: the headline is read OFF the ladder, never recomputed.
HEADLINE_X: float = 0.1

#: Shipped with the number, every time. The threshold is a choice.
CONVENTION_NOTE: str = (
    "x = 0.1 is a REPORTING CONVENTION, not a physical constant. Nothing in the "
    "transport physics, the tissue or the compound marks it out, and no "
    "measurement in this repo selects it. It is NO LONGER HEADLINED ANYWHERE, "
    "and the reason is worth stating because it is a result rather than a "
    "preference: measured on the flagship Xenium GBM section at dx = 12 um, "
    "F(0.1) is 0.0005 — the compound reaches essentially the whole tumour, and "
    "a headline that reads 0.05% on the specimen this product ships with is not "
    "a headline. Retuning the threshold to a value that reads better was "
    "declined: any threshold chosen AFTER seeing what the field does is "
    "post-hoc selection wearing a convention's clothes, which is exactly how "
    "'1x IC50' came to be defended. The headline is instead the delivery spread "
    "(P90/P10 of tumour-cell exposure), which needs no threshold at all and "
    "carries signal even where the map is flat. The full F(x) ladder is still "
    "returned, and is still the evidence — read the curve, not any one point "
    "on it."
)

#: What the ratio IS protected against, at the strength each was measured, and
#: what it is not protected against. Returned in the payload.
SCALE_INVARIANCE_NOTE: str = (
    "The reach ratio is a ratio between two points of ONE map, so a factor that "
    "multiplies every point of the map divides out of it. Two axes were swept. "
    "(1) PLASMA/DOSE SCALE: exact for the finite-window AUC field — worst "
    "relative deviation from an exact rescale 1.41e-14 over plasma factors "
    "0.005 to 26 (a 5,200x span), measured at the shipping configuration "
    "(results/pk_ranking_robustness.md sec. 1.2). "
    "(2) TISSUE CAPACITY: the physics of the cited invariance is CORRECT and the "
    "REGIME CLAIM WAS NOT — read that distinction before reading the source. "
    "capacity_ranking_robustness.md sec. 2b says capacity is exact for a "
    "CONVERGED gradient; that is true. What was wrong, in an earlier version of "
    "this very note, was the assumption that this product's configuration IS a "
    "converged gradient. It is not, and the next reader of that document will "
    "make the same inference unless the gap between its claim and our window is "
    "written down. So what follows is NOT an invariance — it is a BOUNDED "
    "RESIDUAL, "
    "and the bound differs by an order of magnitude between the two readings "
    "below. A capacity term acts as a rescale of time "
    "(AUC_R(x,T) = R * AUC_1(x,T/R)), so it is exact only for a CONVERGED "
    "gradient shape (1e-9, results/capacity_ranking_robustness.md sec. 2b) and "
    "THIS FIELD, AS THIS MODE READS IT, IS NOT AT ITS CONVERGED SHAPE — but NOT "
    "for the reason an earlier version of this note gave, and the correction "
    "matters because the wrong reason is the more quotable one. That version "
    "said the shipping WINDOW is unconverged, citing a gap of 0.01302 to the "
    "infinite-window limit decaying as 281/T with no plateau. The arithmetic "
    "was right and the SUBJECT was wrong: AUC(x,T) = C_ss(x)*T + J(x,T) with J "
    "bounded, so ANY AUC-based ratio carries an O(1/T) tail that can never "
    "plateau no matter how converged the field is — 281/T described the "
    "ESTIMATOR. The field itself, at that same 6 h, is converged to 1.112e-05 "
    "of its own maximum and converges exponentially at the reaction rate "
    "k = 5.00761e-04 1/s (tau = 1997.0 s, so 6 h = 10.82 tau); measured, "
    "results/steady_state_reach.md sec. 1. The capacity residual below is "
    "REAL anyway, and here is the actual mechanism: capacity acts as "
    "AUC_R(x,T) = R * AUC_1(x, T/R), so R SHRINKS THE EFFECTIVE WINDOW — at "
    "R = 39 the shipping 6 h becomes 554 s = 0.28 tau, which genuinely is far "
    "from steady. Capacity does not expose an unconverged window; it "
    "manufactures one. Consequences, measured — and stated PER RUNG, because a "
    "single number over the ladder reports its worst rung while the headline is "
    "its best (a previous version of this note said 'at most 0.0011 over R = 1 "
    "to 39, flat to R = 60, up to 0.80'; re-measured 2026-08-22, no arm and no "
    "reading reproduces any of those three, and read as 'the whole ladder' the "
    "first understates by 365x — results/reach_invariances_measure.md): the "
    "HEADLINE RUNG F(0.1) and every rung below it are BIT-IDENTICAL over "
    "R = 1 to 39 (0.000000) and flat out to R = 80, while the WHOLE 7-RUNG "
    "LADDER moves 0.4021 over that same range because the top rung F(1.0) "
    "carries nearly all of it (0.4056 to 0.8077). Over the full measured "
    "capacity range (3.2 to 4184 mL/g, n = 33) the headline rung moves 0.7409 "
    "and the whole ladder 0.9427, both WORSE than the compound-ordering failure "
    "this product withdrew. The "
    "NORMALISED FALLOFF moves 0.3772 absolute at R = 39 on a curve whose whole "
    "range is 1.0 to 0.51, with the fitted penetration length collapsing 279.7 "
    "to 76.7 um: the SHAPE claim inherits the same disease as the ranking and "
    "MUST NOT be presented as capacity-invariant. ALL OF THAT GOES AWAY IN "
    "STEADY-STATE MODE, which is why the product moved: supply the compiled "
    "drug and the same reductions are taken off the converged field, where a "
    "capacity term is absent from the equations rather than merely small. "
    "(3) TRANSVASCULAR PENETRANCE: EXACT IN EXACT ARITHMETIC, because it is a "
    "pure global multiplier on the source — but the realised arithmetic is "
    "not uniform across the four published quantities and an earlier wording "
    "of this clause said 'bit-identical' about all of them, which is false "
    "for two. Measured at bbb_penetrance 0.1 against 1.0 on the flagship at "
    "three grids: the UNDER-REACHED LADDER and the headline F(0.1) are "
    "bit-identical (max abs difference 0.000e+00, 0 of 7 rungs differing); "
    "E_ref rescales to 1.4e-16 relative; the NORMALISED FALLOFF differs in "
    "10 of 15 bins by <= 3.331e-16, i.e. 1-2 ulp; and the FITTED "
    "PENETRATION LENGTH differs by 3.8e-10 relative at dx = 15 and 2.9e-05 "
    "at dx = 30, because the exponential fit amplifies that round-off. So "
    "anyone diffing two runs at a coarse grid WILL see the fitted length "
    "move, and that is round-off, not physics. That whole invariance is a "
    "LIMIT AS MUCH AS A PROPERTY - this map cannot see a "
    "permeability barrier at all, and depends on the compound only through "
    "L = sqrt(D_eff/k). Do not read it as a barrier or brain-penetration "
    "finding. "
    "Neither invariance protects against error in the transport law itself, "
    "whose descriptor path is unscored (results/descriptor_falsification.md). "
    "This is the list of axes that were SWEPT, not a coverage claim: an axis "
    "not named here was not tested."
)

#: What a STEADY-STATE payload is a statement about. Returned as
#: ``window_note`` so the reframing meets the reader where the number is.
STEADY_STATE_NOTE: str = (
    "STEADY STATE. This map is NOT 'the reach this dose achieved in 6 hours'. "
    "It is 'the REACH THIS TISSUE'S GEOMETRY PERMITS' — a statement about "
    "vessel spacing and this compound's transport, with no dose, no schedule "
    "and no exposure window entangled in it. The field is the converged "
    "concentration C_ss, obtained as the EXACT fixed point of the shipped "
    "time-stepper: physics.diffusion.simulate_drug_transport advances "
    "(I - dt*L + dt*K) c_new = c + dt*S - dt*Adv(c), so setting c_new = c "
    "cancels every dt and leaves the linear system (K - L + Adv) c = S, which "
    "is solved once. It carries no dt and no window, so there is nothing left "
    "to converge and no window convention to defend. The diffusion operator L "
    "is the engine's own (diffusion._build_diffusion_operator, imported); the "
    "upwind advection operator is the only piece re-expressed as a matrix here "
    "and it is checked against diffusion._advective_divergence itself to a "
    "worst relative disagreement of 2.140e-16 over 6 random grids, masks and "
    "velocity fields. WHAT THIS COSTS YOU IS REAL AND IS NOT HIDDEN: the "
    "steady map answers a DIFFERENT QUESTION from the finite-window one, and "
    "it is the wrong question if what you actually want to know is what a "
    "specific schedule delivered in a specific number of hours. A compound "
    "that never reaches steady state in a patient because it is cleared from "
    "plasma first will still be described here by the steady field it would "
    "reach if the plasma level were sustained. Read it as a property of the "
    "TISSUE AND THE COMPOUND'S TRANSPORT, never as a property of a regimen."
)

#: What a FINITE-WINDOW payload is a statement about.
FINITE_WINDOW_NOTE: str = (
    "FINITE WINDOW. This map is the exposure (AUC) field this run accumulated "
    "over its own window, so it IS entangled with the dose schedule and the "
    "window length. It is retained because the production caller cannot yet "
    "supply the compiled DrugParams the steady-state route needs, NOT because "
    "it is the preferred object. Its tissue-capacity behaviour is a bounded "
    "residual rather than an invariance and the bound is large; see "
    "scale_invariance."
)

#: The STEADY-STATE reading of the same three axes. This is the THIRD version of
#: this disclosure in two days; the second is preserved verbatim as
#: :data:`SCALE_INVARIANCE_NOTE` for the finite-window mode, and the correction
#: between them is stated inside both rather than quietly applied. The repo rule
#: adopted this week — PIN THE MAGNITUDE, NOT THE DIRECTION, AND PIN IT AT THE
#: SITE THAT QUOTES THE MEASUREMENT — is why every clause below carries a number.
STEADY_SCALE_INVARIANCE_NOTE: str = (
    "WHAT THIS MAP DOES NOT BUY — FIRST, because two previous versions of this "
    "note were confidently wrong and the ordering should match that history. "
    "(a) It does not make the transport law right. The descriptor path is "
    "unscored (results/descriptor_falsification.md), and an exactly "
    "capacity-invariant answer computed from a wrong D or a wrong k is exactly "
    "as wrong — invariance is a statement about what CANNOT move the answer, "
    "never about the answer. "
    "(b) It does not buy the ABSOLUTE claim, and the absolute claim is "
    "independently believed to be bad in a direction that FLATTERS the "
    "compound: results/kyle2004_penetration.md sec. 6 measured this engine's "
    "perivascular gradient as ~14x too shallow for the one compound where it "
    "is a genuine out-of-sample prediction, and ~3.7x too compressed between "
    "compounds. Steady state changes none of that. A drug that looks like it "
    "reaches the whole tumour here may not. "
    "(c) It does not buy validation against measurement. No measurement of "
    "these compounds' delivery heterogeneity on this tissue exists to score "
    "this map against; robustness is not accuracy "
    "(results/delivery_spread_noise.md sec. 3). "
    "(d) It does not touch the vessel-map inheritance limit, and it cannot see "
    "a permeability barrier at all — see (3). "
    "(e) It answers a DIFFERENT QUESTION from the finite-window map. Read "
    "window_note before comparing the two. "
    "NOW WHAT SURVIVES. The reach ratio is a ratio between two points of ONE "
    "map, so a factor that multiplies every point of the map divides out of "
    "it. Three axes were swept and each is stated at the strength it was "
    "MEASURED at, not at the strength the physics suggests. "
    "(1) PLASMA/DOSE SCALE: exact — worst relative deviation from an exact "
    "rescale 1.41e-14 over plasma factors 0.005 to 26 (a 5,200x span), measured "
    "at the shipping configuration (results/pk_ranking_robustness.md sec. 1.2). "
    "In this mode it is stronger than an invariance: the steady field is solved "
    "with a UNIT source, so no dose and no schedule enters the computation at "
    "all. "
    "(2) TISSUE CAPACITY — A STATEMENT ABOUT THE MAP: THE MAP IS EXACTLY "
    "CAPACITY-INVARIANT. A capacity term R enters as "
    "R*dC/dt = div(D grad C) - div(v C) + S - k*C: it multiplies the TIME "
    "DERIVATIVE AND NOTHING ELSE, and the steady problem has no time "
    "derivative, so R is absent from it by inspection — not small, not "
    "bounded, ABSENT. Measured on the flagship Xenium GBM section over "
    "R = 1 to 4184, on the two quantities this payload publishes: the "
    "UNDER-REACHED LADDER moves 0.000e+00 — exactly zero at every R, because it "
    "counts cells across a threshold and no cell crossed one — and the "
    "NORMALISED FALLOFF moves at most 4.452e-14, which is float64 round-off in "
    "forming D/R, k/R and S/R, not a residual. The same sweep on the "
    "FINITE-WINDOW field moves that falloff 0.3731 at R = 39 alone "
    "(results/steady_state_reach.md sec. 4; "
    "results/capacity_ranking_robustness.md reports 0.3772 for the same "
    "quantity at its own dx and dt, and sec. 2b verifies the steady-state "
    "invariance to 1e-9 on a homogeneous probe). "
    "(2b) TISSUE CAPACITY — A STATEMENT ABOUT AN INSTRUMENT, NOT A WALK-BACK "
    "OF (2): THE PUBLISHED CAPACITY SWEEP IS NOT A FAITHFUL CAPACITY CHANGE "
    "WHEN CONVECTION IS ON. That sweep emulates R by dividing D, binding, "
    "clearance and bbb_penetrance, and it CANNOT divide the velocity field "
    "(results/capacity_ranking_robustness.md sec. 2(d): it is built inside "
    "simulate_response and is unreachable from outside). Multiply the emulated "
    "steady equation through by R and the advection term keeps its R: "
    "div(D grad C) - R*div(v C) + S - k*C = 0. So with convection ON that "
    "sweep applies a capacity change AND an R-fold boost to the Peclet number, "
    "and the map moves — measured on the flagship, the normalised falloff "
    "moves 1.226 at R = 39 and 4.380e+03 at R = 4184, against 2.548e-14 and "
    "1.221e-14 for the same sweep with convection off "
    "(results/steady_state_reach.md sec. 5). Convection OFF hides "
    "this, and THE SHIPPING CONFIGURATION IS CONVECTION ON "
    "(app.server.jobs.PDE_CONVECTION), so anyone re-running that sweep against "
    "this map at shipped settings will see a large residual and should not "
    "read it as a failure of (2). Dividing the velocity by R as well — the "
    "faithful emulation — restores the invariance to 2.787e-14 and "
    "1.177e-14 at those same two R, which is the measurement that "
    "identifies the residual as the instrument's. Both sentences are true "
    "at once because "
    "they have different subjects: (2) is about the model, (2b) is about the "
    "sweep. "
    "(2c) SUPERSESSION: this note replaces a 2026-08-14 version which said the "
    "6 h window was not converged, citing a gap of 0.01302 decaying as 281/T "
    "with no plateau. That was measuring the ESTIMATOR, not the field — the "
    "derivation, the window ladder and the field's actual convergence are in "
    "results/steady_state_reach.md sec. 1. "
    "(3) TRANSVASCULAR PENETRANCE: EXACT IN EXACT ARITHMETIC, because it is a "
    "pure global multiplier on the source — but the realised arithmetic is not "
    "uniform across the four published quantities, and an earlier wording said "
    "'bit-identical' about all of them, which is false for two. Measured at "
    "bbb_penetrance 0.1 against 1.0 on the flagship at three grids: the "
    "UNDER-REACHED LADDER and the headline F(0.1) are bit-identical (max abs "
    "difference 0.000e+00, 0 of 7 rungs differing); E_ref rescales to 1.4e-16 "
    "relative; the NORMALISED FALLOFF differs in 10 of 15 bins by <= 3.331e-16, "
    "i.e. 1-2 ulp; and the FITTED PENETRATION LENGTH differs by 3.8e-10 "
    "relative at dx = 15 and 2.9e-05 at dx = 30, because the exponential fit "
    "amplifies that round-off. Anyone diffing two runs at a coarse grid WILL "
    "see the fitted length move, and that is round-off, not physics. The "
    "invariance is a LIMIT AS MUCH AS A PROPERTY - this map cannot see a "
    "permeability barrier at all, and depends on the compound only through "
    "L = sqrt(D_eff/k). Do not read it as a barrier or brain-penetration "
    "finding. "
    "This is the list of axes that were SWEPT, not a coverage claim: an axis "
    "not named here was not tested."
)

#: How ``E_ref`` is obtained. Recorded in the payload so a consumer can see that
#: the reach number and the penetration length read the same bin.
E_REF_SOURCE: str = (
    "delivery_metrics.penetration_length_fit(...)['drug']['perivascular_mean'] "
    f"— mean AUC over tissue in the first distance bin ({BIN_UM:g} um wide, "
    f">= {MIN_BIN_COUNT} grid points) of the distance-from-nearest-"
    "endothelial-cell profile, REUSED verbatim, never recomputed"
)

#: How ``E_ref`` is obtained in STEADY-STATE mode — and what is LOST by it.
E_REF_SOURCE_STEADY: str = (
    "penetration_length_fit(C_ss, ...)['drug']['perivascular_mean'] — mean "
    f"STEADY CONCENTRATION over tissue in the first distance bin ({BIN_UM:g} um "
    f"wide, >= {MIN_BIN_COUNT} grid points) of the distance-from-nearest-"
    "endothelial-cell profile. NOTE THE LOSS: in finite-window mode this value "
    "was REUSED verbatim from the caller's own delivery_metrics fit, which made "
    "it impossible for the reach number and the published penetration length to "
    "drift apart. That guarantee DOES NOT HOLD HERE. delivery_metrics fits the "
    "finite-window AUC profile and this map is normalised by the STEADY "
    "profile; they are different numbers off different fields, a caller-passed "
    "fit is DISCARDED rather than trusted (see penetration_reused), and the two "
    "penetration lengths are measured against each other in "
    "results/steady_state_reach.md sec. 6 instead of being guaranteed equal."
)

#: How "viable" is defined, in the payload, in the words of the metric it copies.
VIABLE_DEFINITION: str = (
    "a TUMOR cell whose local oxygen at its own grid point satisfies "
    "o >= 2 * OxygenConfig.hypoxia_threshold — the engine's own banding "
    "(tmesim.cells.agents.update_cell_states), identical to the 'viable' band of "
    "delivery_metrics.hypoxic_vs_normoxic_exposure_ratio, so the two cannot "
    "disagree about which cells they describe"
)


#: The two objects this module can produce, and the key that says which one a
#: payload is. A consumer must never have to infer it.
WINDOW_MODES = ("steady_state", "finite_window")


class _Unset:
    """Sentinel for ``convection``: 'the caller did not say'.

    ``convection`` must match the run's ``PDEConfig.convection`` and CANNOT be
    read off ``sim_result``, so a default value would be a guess dressed as a
    contract. A default of ``True`` happens to agree with
    ``app.server.jobs.PDE_CONVECTION`` today, and an edit to either constant
    would silently end that agreement and produce a reach map describing a
    regime the run did not use — the same shape as the capacity/Peclet trap in
    ``results/steady_state_reach.md`` sec. 5, where an unstated convection
    setting made an emulation stop being the thing it claimed to be.

    So in steady-state mode this is REQUIRED, and omitting it REFUSES rather
    than defaults. Finite-window mode never reads it, so it stays optional
    there.
    """

    def __repr__(self) -> str:                        # pragma: no cover
        return "<convection not specified>"


UNSET = _Unset()


def _advection_matrix(vx: np.ndarray, vy: np.ndarray, tissue: np.ndarray,
                      dx: float) -> "_sp.csc_matrix":
    """The linear operator ``diffusion._advective_divergence`` APPLIES, as a matrix.

    ``_advective_divergence(c, ...)`` is linear in ``c`` for a frozen velocity
    field — the upwind switch depends on the sign of the FACE VELOCITY, not on
    ``c`` — so it has an exact matrix representation, and the steady-state solve
    needs it as a matrix rather than as an action. This is the only piece of the
    engine's discretisation that is re-expressed here rather than imported, and
    it is therefore the only piece that can silently disagree with the engine.
    That risk is not argued away: ``tests/test_reach.py::TestSteadyStateOperator``
    applies both to random fields and requires them equal to machine precision.

    Face by face, exactly as the engine does it: ``vf`` is the mean of the two
    adjacent cell velocities, the donor cell supplies the concentration, the
    flux is zeroed on domain-edge and tissue/non-tissue faces, and the same flux
    leaves ``P`` and enters ``N``.
    """
    ny, nx = tissue.shape
    idx = np.arange(ny * nx).reshape(ny, nx)
    rows: list = []
    cols: list = []
    data: list = []

    def _face(p_i: np.ndarray, n_i: np.ndarray, f: np.ndarray) -> None:
        # flux = f * c_upwind;  div[P] += flux;  div[N] -= flux
        from_p = np.where(f > 0.0, f, 0.0)     # donor is P
        from_n = np.where(f > 0.0, 0.0, f)     # donor is N
        rows.extend([p_i, n_i, p_i, n_i])
        cols.extend([p_i, p_i, n_i, n_i])
        data.extend([from_p, -from_p, from_n, -from_n])

    if nx > 1:
        vf = 0.5 * (vx[:, :-1] + vx[:, 1:])
        active = tissue[:, :-1] & tissue[:, 1:]
        f = np.where(active, vf, 0.0) / dx
        _face(idx[:, :-1].ravel(), idx[:, 1:].ravel(), f.ravel())
    if ny > 1:
        vf = 0.5 * (vy[:-1, :] + vy[1:, :])
        active = tissue[:-1, :] & tissue[1:, :]
        f = np.where(active, vf, 0.0) / dx
        _face(idx[:-1, :].ravel(), idx[1:, :].ravel(), f.ravel())

    n = ny * nx
    if not rows:
        return _sp.csc_matrix((n, n))
    return _sp.csc_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n))


def steady_state_concentration(
    micro: Microenvironment,
    drug: DrugParams,
    source: np.ndarray,
    velocity: Optional[np.ndarray] = None,
    exchange_k: Optional[np.ndarray] = None,
) -> np.ndarray:
    """The tissue's converged drug concentration field ``C_ss`` (uM). ONE solve.

    **This is not an approximation of the engine — it is the engine's own limit,
    written down.** ``diffusion.simulate_drug_transport`` advances

    ``(I - dt*L + dt*K) c_new  =  c + dt*S - dt*Adv(c)``

    (backward Euler on diffusion and the first-order sink, explicit upwind
    advection, source at the old time level). Set ``c_new = c = c*`` and every
    ``dt`` cancels identically, leaving the LINEAR SYSTEM

    ``(K - L + Adv) c*  =  S``

    with the SAME ``L`` (:func:`diffusion._build_diffusion_operator`, imported),
    the same ``K = clearance + binding`` and the same upwind ``Adv``. So ``c*``
    is the exact fixed point of the shipped stepper, and — this is the point —
    it carries **no ``dt`` and no window**. There is nothing left to converge.

    Why this is the cheap route as well as the exact one
    ---------------------------------------------------
    The time march pays one sparse factorisation AND ``total_time/dt`` triangular
    solves (2,160 of them at the shipping 6 h / 10 s). This pays the SAME ONE
    factorisation and a SINGLE triangular solve — 0.90 s against 43.70 s at the
    shipping grid, and more accurate as well (the march's agreement with its own
    fixed point floors at 1.99e-13 and stops improving past ~24 h of simulated
    time, because that is backward-Euler error at ``dt`` and more time does not
    remove it).

    **Beware the obvious misreading.** That ratio is about THE FIELD. A run is
    not 48x faster; a run still marches, because ``kill_fraction``,
    ``delivery_metrics`` and the oxygen field all need it. This is computed in
    addition, and adds +4.2 % to a run at dx = 12. See
    ``results/steady_state_reach.md`` sec. 2 and 2.2 — 2.2 exists because the
    first draft of that file quoted the ratio in a way that invited exactly this
    error.

    Non-tissue cells
    ----------------
    ``L``'s non-tissue rows are already all-zero, ``K`` is zeroed there and
    ``Adv`` only has entries on tissue-tissue faces, so those rows would make the
    matrix singular. A unit diagonal is ADDED on them (their off-diagonal entries
    are already zero) and their right-hand side is zero, which holds them at
    ``C = 0`` — the same thing the stepper's identity rows do.

    Parameters
    ----------
    micro : Microenvironment
        Supplies the grid, the tissue mask and ``fields.ecm_density``.
    drug : DrugParams
        Supplies ``diffusivity`` (through the SAME tortuosity law the engine
        uses, :func:`diffusion._effective_diffusivity`), ``clearance`` and
        ``binding``.
    source : np.ndarray
        The delivery field ``S`` (uM/s), i.e. what
        :func:`treatment.drugs.delivery_source` returns. Any SCALAR factor in it
        — the plasma level, ``bbb_penetrance`` — divides out of every ratio this
        module reports, so a caller that only wants the map may pass the unit
        source. That is a measured invariance, not an assumption; see
        :data:`SCALE_INVARIANCE_NOTE`.
    velocity : np.ndarray, optional
        ``(2, ny, nx)`` interstitial velocity from
        :func:`physics.pressure.solve_ifp`. ``None`` solves without convection.
    exchange_k : np.ndarray, optional
        The transvascular back-flux rate, added to ``K`` exactly as the engine
        adds it. ``None`` (the shipped default,
        ``PDEConfig.transvascular_exchange = False``) reproduces the one-way
        source.

    Returns
    -------
    np.ndarray
        ``(ny, nx)`` steady concentration in uM, zero outside the tissue mask.
    """
    domain = micro.domain
    ny, nx = domain.ny, domain.nx
    dx = float(domain.dx)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    n = ny * nx

    d_eff = _effective_diffusivity(drug, micro.fields.ecm_density)
    k = np.full((ny, nx), float(drug.clearance) + float(drug.binding))
    if exchange_k is not None:
        k = k + np.asarray(exchange_k, dtype=np.float64)
    k[~tissue] = 0.0

    src = np.asarray(source, dtype=np.float64)
    if src.shape != (ny, nx):
        raise ValueError(f"source shape {src.shape} != grid {(ny, nx)}")
    src = np.where(tissue, src, 0.0)

    diag = k.ravel() + (~tissue).ravel().astype(np.float64)
    a = (-_build_diffusion_operator(d_eff, tissue, dx)
         + _sp.csc_matrix((diag, (np.arange(n), np.arange(n))), shape=(n, n)))
    if velocity is not None:
        vel = np.asarray(velocity, dtype=np.float64)
        if vel.shape != (2, ny, nx):
            raise ValueError(f"velocity shape {vel.shape} != {(2, ny, nx)}")
        a = a + _advection_matrix(vel[0], vel[1], tissue, dx)

    c = _factorized(a.tocsc())(src.ravel()).reshape(ny, nx)
    c[~tissue] = 0.0
    return c


def null_reach(reason: str, steady: bool = False, **known) -> Dict:
    """A refusal. Every quantity that could not be computed is absent, not
    defaulted; anything that genuinely IS known (cell counts, which object was
    being built, the convention text) is still reported.

    ``steady`` selects the same mode-dependent disclosures the success payload
    carries, so a consumer rendering a refusal shows the SAME sentences it
    would have shown beside the number — including which object failed to be
    computed."""
    out: Dict = {
        "null_reason": reason,
        "headline_under_reached_fraction": None,
        "headline_x": float(HEADLINE_X),
        "curve": None,
        "falloff": None,
        "E_ref_uM_s": None,
        "E_ref_source": (E_REF_SOURCE_STEADY if steady else E_REF_SOURCE),
        "E_ref_units": ("uM" if steady else "uM*s"),
        "window_mode": ("steady_state" if steady else "finite_window"),
        "window_s": None,
        "window_note": (STEADY_STATE_NOTE if steady else FINITE_WINDOW_NOTE),
        "penetration_reused": False,
        "n_viable": None,
        "n_tumor": None,
        "viable_definition": VIABLE_DEFINITION,
        "hypoxia_threshold": None,
        "convention": CONVENTION_NOTE,
        "scale_invariance": (STEADY_SCALE_INVARIANCE_NOTE if steady
                             else SCALE_INVARIANCE_NOTE),
        "penetration_fit_error": None,
        "units": (_UNITS_STEADY if steady else _UNITS),
    }
    out.update(known)
    return out


_UNITS = {
    "E_ref": "uM*s",
    "reach_ratio": "dimensionless (AUC_i / E_ref)",
    "under_reached_fraction": "fraction of VIABLE tumour cells, in [0, 1]",
    "distance_um": "um",
    "normalised_mean_auc": "dimensionless (bin mean AUC / E_ref)",
}

#: The steady-state reading of the same reductions. ``E_ref`` is a
#: CONCENTRATION, not an exposure, so its unit genuinely changes — the wire key
#: is still ``E_ref_uM_s`` because ``app/web`` and ``app/server/report.py`` read
#: that name, and breaking a wire contract to fix a unit label would be the
#: worse trade. The payload therefore carries the true unit HERE and in
#: ``E_ref_units``, and a consumer that hard-codes "uM.s" instead of reading it
#: is wrong in steady-state mode. Two such sites exist and are named in
#: ``results/steady_state_reach.md`` sec. 7.
_UNITS_STEADY = {
    "E_ref": "uM",
    "reach_ratio": "dimensionless (C_ss,i / E_ref)",
    "under_reached_fraction": "fraction of VIABLE tumour cells, in [0, 1]",
    "distance_um": "um",
    "normalised_mean_auc": "dimensionless (bin mean C_ss / E_ref)",
}


def compute_reach(
    micro: Microenvironment,
    sim_result: SimulationResult,
    distance_um: Optional[np.ndarray] = None,
    penetration: Optional[Dict] = None,
    ladder: Sequence[float] = REACH_LADDER,
    hypoxia_threshold: float = OxygenConfig().hypoxia_threshold,
    drug: Optional[DrugParams] = None,
    velocity: Optional[np.ndarray] = None,
    exchange_k: Optional[np.ndarray] = None,
    convection=UNSET,
) -> Dict:
    """The scale-free reach map for one (tissue, simulation) pair.

    WHICH OBJECT THIS RETURNS — the payload always says, in ``window_mode``
    ---------------------------------------------------------------------
    * ``drug`` GIVEN -> ``window_mode = "steady_state"``. The map is built on
      the converged concentration field :func:`steady_state_concentration`,
      and the claim it supports is **"the reach this tissue's geometry
      permits"** — a statement about vessel spacing and this compound's
      transport, with no dosing schedule and no window in it. Tissue capacity
      is then EXACTLY invariant rather than bounded, because a capacity term
      multiplies the time derivative and the steady problem has none.
    * ``drug`` OMITTED -> ``window_mode = "finite_window"``. The historical
      object: the map of ``sim_result.fields['auc']``, i.e. the reach this dose
      achieved in this run's window. Every quantity is the same reduction of a
      different field. It is retained because the only production caller
      (``app/server/jobs.py``) cannot yet supply ``drug``; it is NOT the
      preferred object and its capacity residual is large and measured.

    The two are not interchangeable and are never mixed in one payload.

    Parameters
    ----------
    micro : Microenvironment
        The reconstructed real microenvironment the simulation ran on.
    sim_result : SimulationResult
        Output of :func:`tmesim.sim.engine.simulate_response`. ``fields
        ['oxygen']`` is required in BOTH modes (it defines the viable band);
        ``fields['auc']`` is required only in finite-window mode.
    drug : DrugParams, optional
        Supplying it selects steady-state mode. It must be the SAME compiled
        ``DrugParams`` the run used — the one ``simulate_response`` resolved —
        or the map describes a different compound from the rest of the payload.
    velocity : np.ndarray, optional
        ``(2, ny, nx)`` interstitial velocity, steady-state mode only. When
        omitted it is obtained from :func:`physics.pressure.solve_ifp` on
        ``micro``, which is BIT-IDENTICAL to the array the engine's own
        transport call received: ``solve_ifp`` reads only the grid,
        ``vascular_density``, ``ecm_density`` and the cells' ``cell_type``, and
        the engine changes none of those before calling it. Asserted in
        ``tests/test_reach.py::TestSteadyStateField``.
    convection : bool
        **REQUIRED in steady-state mode; there is no default.** It must match
        the run's ``PDEConfig.convection`` and it cannot be read off
        ``sim_result``, so omitting it REFUSES (``null_reason``) rather than
        assuming. A default of ``True`` would agree with the shipped
        ``app.server.jobs.PDE_CONVECTION`` today and would silently stop
        agreeing the moment either constant moved, which is exactly how the
        capacity emulation in ``results/steady_state_reach.md`` sec. 5 stopped
        being a capacity emulation. ``False`` drops the advection term entirely
        and is what the capacity studies use, for the reason
        ``results/capacity_ranking_robustness.md`` sec. 2(d) gives. Ignored in
        finite-window mode, which never builds a field.
    exchange_k : np.ndarray, optional
        Transvascular back-flux rate, steady-state mode only. Must match the
        run: pass :func:`treatment.drugs.exchange_coefficient` when the run set
        ``PDEConfig.transvascular_exchange``, and ``None`` otherwise — which is
        the shipped default (``transvascular_exchange = False``).
    distance_um : np.ndarray, optional
        Pre-computed distance-from-nearest-vessel map (um). The distance
        transform is the expensive part of this computation and it is pure
        geometry, so a caller that has already built one (for
        :func:`delivery_metrics.compute_all`, say) should pass the SAME array
        here rather than paying for it twice.
    penetration : dict, optional
        An already-computed :func:`delivery_metrics.penetration_length_fit`
        result. **Finite-window mode only, and this is a real loss under steady
        state.** Pass it when the caller has one — ``E_ref`` is read out of
        ``penetration['drug']['perivascular_mean']`` either way, so passing it
        only avoids re-running the fit. When omitted the fit is run here with
        ``include_source_null=False`` (the source null plays no part in this
        metric and would double the binning work); ``['drug']
        ['perivascular_mean']`` is the identical number either way.

        In STEADY-STATE mode a passed ``penetration`` is **ignored**, and the
        payload says so in ``penetration_reused``. It has to be: the caller's
        fit is a fit to the finite-window AUC profile, and the steady map's
        ``E_ref`` is the first bin of the STEADY profile. Those are different
        numbers off different fields. The consequence is that the guarantee
        "the reach number and the published penetration length can never drift
        apart" — which held by construction in finite-window mode — DOES NOT
        HOLD in steady-state mode, because ``delivery_metrics`` publishes an
        AUC-fitted penetration length and this map no longer reads it. That is
        a genuine cost of the change and it is disclosed rather than papered
        over; the two lengths are measured against each other in
        ``results/steady_state_reach.md`` sec. 6.
    ladder : sequence of float
        Reach ratios at which ``F`` is evaluated. Defaults to
        :data:`REACH_LADDER`; the headline is picked off it.
    hypoxia_threshold : float
        Passed through to the "viable" band, defaulting to
        ``OxygenConfig.hypoxia_threshold`` — the same default
        :func:`delivery_metrics.hypoxic_vs_normoxic_exposure_ratio` uses.

    Returns
    -------
    dict
        ``headline_under_reached_fraction`` (``F(0.1)``), ``headline_x``,
        ``curve`` (``{'x': [...], 'under_reached_fraction': [...]}``),
        ``falloff`` (``{'distance_um', 'normalised_mean_auc', 'count',
        'bin_um', 'min_bin_count'}``), ``E_ref_uM_s`` + ``E_ref_source``,
        ``n_viable``, ``n_tumor``, ``viable_definition``, ``convention``,
        ``scale_invariance``, ``penetration_fit_error``, ``units`` and
        ``null_reason``.

        ``null_reason`` is ``None`` on success and a sentence otherwise, and
        when it is set every derived number is ``None``. See the module
        docstring: this refuses rather than defaults.

    Note on ``penetration_fit_error``
    ---------------------------------
    A failure of the exponential-plus-floor fit is REPORTED (as
    ``penetration_fit_error``) but does not null this result, because no
    quantity here reads a fitted parameter: ``E_ref`` is a bin mean, and the
    ladder and the falloff are ratios of measured bin/cell values. Nulling on it
    would be discarding a number that is defined, which is the opposite of the
    refusal rule.
    """
    x = np.asarray(list(ladder), dtype=np.float64)
    # Which object is being built. Fixed HERE, before anything can fail, so
    # every refusal below carries the right disclosures too.
    steady = drug is not None

    # --- 1. the tumour cells and their exposure ---------------------------- #
    exp = tumor_cell_exposure(micro, sim_result)
    n_tumor = int(exp["n_tumor"])
    auc_cells = np.asarray(exp["auc_uM_s"], dtype=np.float64)
    if n_tumor == 0:
        return null_reach(
            "no TUMOR cells in this reconstruction, so there is no tumour for "
            "the drug to reach or miss", n_tumor=0, steady=steady)

    # --- 2. the viable band, defined exactly as metric 5 defines it --------- #
    oxygen = sim_result.fields.get("oxygen")
    if oxygen is None:
        return null_reach(
            "sim_result.fields carries no 'oxygen' field, so the VIABLE band "
            "(o >= 2 * hypoxia_threshold) that this fraction is taken over "
            "cannot be identified; no substitute population is used",
            n_tumor=n_tumor, steady=steady)
    o2 = _sample_at_cells(oxygen, micro, exp["tumor_mask"])
    viable = np.asarray(o2, dtype=np.float64) >= 2.0 * float(hypoxia_threshold)
    n_viable = int(viable.sum())
    if n_viable == 0:
        return null_reach(
            "no VIABLE tumour cells in this reconstruction (no tumour cell has "
            f"local oxygen >= 2 * hypoxia_threshold = "
            f"{2.0 * float(hypoxia_threshold):g}), so a fraction OF the viable "
            "population is undefined; it is not fallen back to all tumour cells",
            n_tumor=n_tumor, n_viable=0, steady=steady,
            hypoxia_threshold=float(hypoxia_threshold))

    # --- 2b. WHICH FIELD this is a map of ----------------------------------- #
    # The two modes differ HERE and only here: every reduction below is the
    # same reduction of whichever field this block selects.
    if steady:
        if isinstance(convection, _Unset):
            return null_reach(
                "steady-state mode was selected (a drug was supplied) but "
                "`convection` was not, and it cannot be recovered from "
                "sim_result: PDEConfig.convection is not carried on a "
                "SimulationResult. The steady field would then be built for a "
                "transport regime that may not be the one the run used, so no "
                "map is returned. Pass convection=<the run's "
                "PDEConfig.convection>.",
                steady=True, n_tumor=n_tumor, n_viable=n_viable,
                hypoxia_threshold=float(hypoxia_threshold))
        if not convection:
            velocity = None
        elif velocity is None:
            from tmesim.physics.pressure import solve_ifp
            _ifp, velocity = solve_ifp(micro)
        from tmesim.treatment.drugs import delivery_source
        # Unit plasma: every scalar factor in the source — the plasma level and
        # bbb_penetrance both are — multiplies every point of the field and so
        # divides out of every ratio reported below. That is measured, not
        # assumed (SCALE_INVARIANCE_NOTE items 1 and 3), and it is exactly why
        # the steady map needs no dose and no schedule to be well defined.
        field = steady_state_concentration(
            micro, drug, delivery_source(micro, drug, 1.0),
            velocity=velocity, exchange_k=exchange_k)
    else:
        field = np.asarray(sim_result.fields["auc"], dtype=np.float64)
    values_at_cells = (
        np.asarray(_sample_at_cells(field, micro, exp["tumor_mask"]),
                   dtype=np.float64) if steady else auc_cells)

    # --- 3. the perivascular reference -------------------------------------- #
    # The distance transform is the expensive step and it is pure geometry, so
    # it is resolved ONCE here and handed to both the fit and the falloff.
    if distance_um is None:
        distance_um = _distance_for(micro)
    penetration_reused = bool(penetration is not None and not steady)
    if steady or penetration is None:
        # In steady-state mode a caller-supplied fit is DISCARDED rather than
        # trusted: it is a fit to the finite-window AUC profile and this map's
        # reference is the first bin of the STEADY profile. Reusing it would
        # normalise the steady field by a number taken off a different field.
        penetration = penetration_length_fit(
            field, micro, distance_um=distance_um, include_source_null=False)
    fit_error = penetration["drug"].get("error")
    e_ref = float(penetration["drug"]["perivascular_mean"])
    if not np.isfinite(e_ref) or e_ref <= 0.0:
        return null_reach(
            f"the perivascular reference exposure E_ref is {e_ref!r}: the "
            "distance-from-vessel profile has no usable first bin (a bin needs "
            f">= {MIN_BIN_COUNT} grid points) or its mean AUC is not positive, "
            "so every reach ratio would be a division by it",
            n_tumor=n_tumor, n_viable=n_viable, steady=steady,
            hypoxia_threshold=float(hypoxia_threshold),
            penetration_fit_error=fit_error)

    # --- 4. the curve ------------------------------------------------------- #
    # F(x) = fraction of VIABLE tumour cells with r_i < x. Evaluated by binary
    # search on the sorted ratios, which makes monotonicity in x exact rather
    # than incidental, and keeps the strict "<" of the definition.
    r = values_at_cells[viable] / e_ref
    r_sorted = np.sort(r)
    counts = np.searchsorted(r_sorted, x, side="left")
    F = counts.astype(np.float64) / float(n_viable)

    on_ladder = np.flatnonzero(x == float(HEADLINE_X))
    if on_ladder.size == 0:
        return null_reach(
            f"the headline threshold x = {HEADLINE_X:g} is not on the ladder "
            f"{[float(v) for v in x]} this call was given, and it is read OFF "
            "the ladder rather than recomputed, so no headline is reported",
            n_tumor=n_tumor, n_viable=n_viable, E_ref_uM_s=e_ref,
            steady=steady, hypoxia_threshold=float(hypoxia_threshold),
            penetration_fit_error=fit_error)
    headline_idx = int(on_ladder[0])

    # --- 5. the scale-free falloff ----------------------------------------- #
    # The SAME bins that produced E_ref: same distance map, same binning, same
    # usable-bin rule, so falloff[0] is 1.0 by construction.
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    prof = distance_profile(field, distance_um, tissue)
    d_bins, m_bins, n_bins = usable_bins(prof)
    if m_bins.size == 0 or float(m_bins[0]) != e_ref:
        # The falloff must be the same profile E_ref came from. If it is not,
        # something upstream has changed the binning between the two calls and
        # the honest answer is to say so, not to publish two profiles that
        # disagree.
        return null_reach(
            "the distance profile rebuilt here does not reproduce the "
            "perivascular mean that delivery_metrics.penetration_length_fit "
            f"reported (E_ref = {e_ref!r}); the falloff and the reference "
            "would come from different binnings, so neither is reported",
            n_tumor=n_tumor, n_viable=n_viable, E_ref_uM_s=e_ref,
            steady=steady, hypoxia_threshold=float(hypoxia_threshold),
            penetration_fit_error=fit_error)

    falloff = m_bins / e_ref

    return {
        "null_reason": None,
        "headline_under_reached_fraction": float(F[headline_idx]),
        "headline_x": float(HEADLINE_X),
        "curve": {
            "x": [float(v) for v in x],
            "under_reached_fraction": [float(v) for v in F],
        },
        "falloff": {
            "distance_um": [float(v) for v in d_bins],
            "normalised_mean_auc": [float(v) for v in falloff],
            "count": [int(v) for v in n_bins],
            "bin_um": float(prof["bin_um"]),
            "min_bin_count": int(MIN_BIN_COUNT),
        },
        "E_ref_uM_s": e_ref,
        "E_ref_source": (E_REF_SOURCE_STEADY if steady else E_REF_SOURCE),
        "E_ref_units": ("uM" if steady else "uM*s"),
        "window_mode": ("steady_state" if steady else "finite_window"),
        "window_s": (None if steady else _window_of(sim_result)),
        "window_note": (STEADY_STATE_NOTE if steady else FINITE_WINDOW_NOTE),
        "penetration_reused": penetration_reused,
        "n_viable": n_viable,
        "n_tumor": n_tumor,
        "viable_definition": VIABLE_DEFINITION,
        "hypoxia_threshold": float(hypoxia_threshold),
        "convention": CONVENTION_NOTE,
        "scale_invariance": (STEADY_SCALE_INVARIANCE_NOTE if steady
                             else SCALE_INVARIANCE_NOTE),
        "penetration_fit_error": fit_error,
        "units": (_UNITS_STEADY if steady else _UNITS),
    }


def _window_of(sim_result: SimulationResult) -> Optional[float]:
    """The exposure window the finite-window field was accumulated over (s),
    read off the transport result's own last time stamp rather than from a
    config the caller might have changed after the solve."""
    t = getattr(getattr(sim_result, "transport", None), "times", None)
    if t is None or len(t) == 0:
        return None
    return float(np.asarray(t)[-1])


def _distance_for(micro: Microenvironment) -> np.ndarray:
    """The PRIMARY distance-from-nearest-vessel map for this tissue.

    Only reached when a caller supplied no ``distance_um``; it is the same
    direct nearest-endothelial-cell distance
    (:func:`vessel_profile.direct_vessel_distance_um`) that
    :func:`delivery_metrics.penetration_length_fit` defaults to, so the two
    profiles are the same profile.
    """
    from tmesim.validation.vessel_profile import (
        direct_vessel_distance_um,
        endothelial_xy,
    )

    return direct_vessel_distance_um(micro.domain, endothelial_xy(micro.cells))
