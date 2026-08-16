"""Time-resolved pharmacokinetics: the drug layer's schedule capability.

Before this, ``sim.engine.simulate_response`` collapsed the multi-dose plasma
profile to a single time-averaged level and drove the transport solver with a
STATIC source. Two schedules with equal cumulative plasma exposure were then
bit-identical inputs, so fractionation and infusion duration produced a kill
spread of exactly 0.0 — a structural inability to answer a scheduling question,
not a finding about biology (see ``results/dose_schedule.md``).

``PDEConfig.time_resolved=True`` (and, at the solver level, the ``source_profile``
argument of :func:`tmesim.physics.diffusion.simulate_drug_transport`) drives the
source with the real ``plasma(t)`` instead. This module pins:

**Backward compatibility (must be BIT-identical, not merely close)**
  1. the default static path reproduces a stored pre-change reference exactly,
  2. an all-ones ``source_profile`` is bit-identical to ``source_profile=None``,
  3. ``PDEConfig.time_resolved`` defaults to False.

**The physics (checks that fail if the implementation is wrong)**
  4. *Convergence.* As the plasma level is held constant, the time-resolved path
     converges to the static-source result.
  5. *Mass balance.* Conservation still holds exactly with a time-varying source,
     and matched-cumulative-exposure schedules inject exactly the same total drug.
  6. *Depth / direction.* At matched total exposure over a fixed window the
     schedule genuinely changes the exposure field, and it does so in a
     DERIVABLE direction — see :func:`test_earlier_delivery_reaches_deeper`.
  7. *Capability.* Fractionation at matched cumulative exposure no longer gives a
     spread of exactly zero.

Why check 6 goes the opposite way to the naive expectation
----------------------------------------------------------
The usual transport-limited intuition says a long slow infusion should reach
DEEPER than a short bolus. For this (linear) transport model that is provably
false, and the proof also says exactly how big the schedule effect can be.
Integrate ``dC/dt = div(D grad C) - div(vC) + S - kC`` over the window ``[0, T]``
and write ``A = AUC = int_0^T C dt``::

    C(T) - C(0) = div(D grad A) - div(v A) - k A + int_0^T S dt

With ``C(0) = 0`` and ``M := k - div(D grad .) + div(v .)``::

    A = M^{-1} ( S_total - C(T) ),      S_total = int_0^T S dt

``S_total`` is IDENTICAL for any two schedules at matched cumulative plasma
exposure. So the entire schedule effect flows through ``C(T)``, the drug still
in flight at the end of the window. A bolus given early has washed out by ``T``
(``C(T) ~ 0``) and therefore realises the full steady-state exposure field; a
sustained infusion still running at ``T`` leaves ``C(T) > 0`` un-counted, so it
loses exposure. Two consequences, both tested here:

* Earlier delivery yields MORE exposure and, because the un-counted term
  ``M^{-1}C(T)`` is itself a doubly-smoothed (flatter) field, also a relatively
  FLATTER profile — i.e. the bolus reaches deeper, not the infusion.
* The effect is bounded by the terminal residual and therefore VANISHES as the
  window grows past washout. Schedule matters exactly insofar as the window
  truncates the response.

The 2026-08-12 re-pin of ``tests/data/static_source_reference.json``
-------------------------------------------------------------------
The stored reference exists so that nobody can change static-source behaviour
without it showing up. That only works if a re-pin is loud. This one was, and is
recorded here as well as in the reference file's own ``_README._repin_history``.

**Cause.** ``tmesim/physics/pressure.py`` smoothed the tumour-occupancy field
``phi`` with ``_TUMOR_SIGMA = 5.0``, commented "grid cells" but handed straight
to ``scipy.ndimage.gaussian_filter``, whose ``sigma`` is in PIXELS. ``phi`` sets
the vessel-leak gain and the lymphatic suppression, so it is a MODEL parameter:
the physical smoothing length was therefore ``5 * dx`` — 60 um at the production
grid but 75 um at ``dx = 15`` and 400 um at ``dx = 80``. Refining the mesh
silently redefined the model. It is now ``_TUMOR_SIGMA_UM = 60.0``, converted to
pixels at the point of use.

**What moved, and why that is the predicted result rather than a regression.**
``phi`` reaches a transport result by exactly one route — the Darcy velocity
``v = -K grad p`` — so it can only bite when convection is enabled. That is
precisely the observed split: the three ``dx = 15`` convection-ON cases and the
``dx = 80`` real-tissue case moved (by 1-3 ULP in the transport sums; more in
``penetration_index``, which is a shape statistic); ``weak_10dose_noconv``
(convection OFF) and ``solver_flat`` (raw solver, no IFP at all) are byte-for-byte
unchanged and were not regenerated. ``plasma_mean``, a pure-PK quantity with no
spatial dependence, is unchanged in every case. The superseded values are kept
verbatim in the reference file under
``_superseded_2026-08-12_grid_cell_tumor_sigma``.

**The guard was strengthened, not weakened.** ``60.0 / 12.0 == 5.0`` exactly in
IEEE-754, so the production grid is the fixed point of the change and every
shipped ``dx = 12`` result is bit-identical — but until now ``dx = 12`` was the
one grid with *no* bitwise pin, while the grids that did have pins were the ones
whose values were wrong. :func:`test_production_dx12_static_path_matches_stored_reference_bitwise`
closes that hole. Its values were generated BEFORE the fix and verified
identical AFTER it, so it pins the invariance claim itself. Full write-up:
``results/tumor_sigma_units.md``.

The 2026-08-13 re-pin — the two-route transport model became the default
-----------------------------------------------------------------------
**Cause.** ``compile_drug.build_profile`` now applies a two-route cell term
(paracellular between cells plus transcellular through membranes), so
temozolomide's compiled diffusivity went ``277.985 -> 42.4245 um^2/s`` and the
single perivascular-sink constant ``alpha``, fitted in closed form from
``L = sqrt(D_eff/k_eff)`` against the same measured 64.92 um anchor, went
``148.60815765445545 -> 22.223477711172563``. Both cases below use a COMPILED,
alpha-calibrated temozolomide, so both had to move.

**Why the new values are right rather than merely different.** The superseded
values encode a diffusivity law whose only solute input was molecular weight,
which was scored blind against 71 measured compounds and REFUTED at R^2 = 0.022
(``results/diffusivity_falsification.md``). The law now in force scores
R^2 = 0.718 on the same 71 compounds, blind and unfitted, and brings the engine's
dimensionless multicellular-layer hindrance factor inside its measured band for
the first time (1.77-fold -> 10.33-fold against a measured 12-fold, PMID
9806526). These cases now run a transport coefficient that has been tested
externally and largely agrees; they did not before.

**What moved, and why that is the predicted result.** Only a case whose drug
comes from the molecular compile can move, and only those did:
``real_xenium_dx80_tmz_3dose`` and ``production_xenium_dx12_tmz_3dose`` were
regenerated; ``default_1dose``, ``weak_1dose``, ``weak_3dose``,
``weak_10dose_noconv`` and ``solver_flat`` all use hand-built ``DrugParams``
literals with no compile in the path and are byte-for-byte unchanged; and
:func:`test_production_dx12_ifp_solve_matches_stored_reference_bitwise`, which
has no drug in it at all, passed untouched. ``plasma_mean`` — pure PK, no spatial
dependence — is BIT-IDENTICAL in the cases that moved, which is the positive
control.

**The magnitudes are a finding, not a footnote.** Unlike the 2026-08-12 ULP-level
re-pin, ``kill_fraction`` moves 2.76x at ``dx = 12`` and 5.29x at ``dx = 80``.
Diffusivity FELL 6.6x and exposure ROSE, because the Primeau anchor pins the
ratio ``D/k`` and the fitted loss rate fell with ``D``. The p99-normalised ratio
the breast panel scores is immune to that common-mode change — 0 of 12 of those
verdicts moved — but threshold statistics such as kill fraction are not. Full
write-up: ``results/two_route_promotion.md``; superseded values under
``_superseded_2026-08-13_mw_only_diffusivity``.

The 2026-08-14 re-pin — the area-conserving vessel painting became the default
------------------------------------------------------------------------------
**Cause.** ``reconstruct.vessel_network.discrete_vascular_field``'s ``geometry``
default moved ``'binary' -> 'area'``
(``results/area_conservation_promotion.md``). The binary painting asked each grid
cell a YES/NO question — *"is your centre within 7.5 µm of an endothelial
cell?"* — and deposited a whole ``dx²`` cell of vessel per yes. That is a COUNT
dressed as a field: it carries a ``dx²`` factor, so the total painted vessel area
was a function of the mesh, spanning **4.7312 %** non-monotonically over
``dx = 30…5 µm``. It is the same defect ``physics.oxygen._tumor_density`` already
fixed on the oxygen SINK by dividing a count by the grid-cell area.
``geometry='area'`` paints the same discs at the same radius around the same
endothelial cells and gives each cell the covered-area FRACTION, conserving the
total exactly — 742,713.0 µm² at every rung, span 6.3e-14 %.

**Why the new values are right rather than merely different — and this one is
checkable arithmetic rather than a judgement.** The transport integrals are
LINEAR in the source total, so they must move by exactly the source total's
change. They do:

===================  ==========================  ====================
grid                 tissue-integrated source    ``mean_auc`` / ``auc_sum`` / ``final_c_sum``
===================  ==========================  ====================
``dx = 12 µm``       5199.0000 → 5157.7292       all three **−0.794 %**
                     (**−0.7938 %**)
``dx = 80 µm``       104.0000 → 113.9567         all three **+9.574 %**
                     (**+9.5738 %**)
===================  ==========================  ====================

Three independent integrals reproducing the source change to three decimals, on
two grids, with opposite signs, is not something a regression does. The signs
differ because the binary total is a *different number at every rung*: at
``dx = 80 µm`` a 7.5 µm disc is 1/36 of a grid cell so the centre-vote misses
most vessels, whereas ``dx = 12 µm`` happens to sit near the TOP of the
quantisation lottery (748,656 µm² against a ladder mean of 741,882).

**What moved, and why that is the predicted result.** The painting reaches a
result only through a reconstruction built by
``real_fields.vascular_density_from_cells``, so only REAL-TISSUE cases can move.
Exactly those did: ``real_xenium_dx80_tmz_3dose`` and
``production_xenium_dx12_tmz_3dose`` were regenerated, including this time the
IFP pin — which was untouched by the 2026-08-13 re-pin because it has no drug in
it, and moves now because it reads ``fields.vascular_density`` directly.
``default_1dose``, ``weak_1dose``, ``weak_3dose``, ``weak_10dose_noconv``
(hand-built synthetic microenvironments) and ``solver_flat`` (raw solver, no
reconstruction at all) are byte-for-byte unchanged. ``plasma_mean`` is
BIT-IDENTICAL in both moved cases — the positive control again.

**What moved non-proportionally, named so it is not mistaken for drift.**
``kill_fraction`` moves −3.361 % at ``dx = 12`` and +19.614 % at ``dx = 80``
because kill is a saturating Hill function of exposure, not a linear one.
``penetration_index`` rises +1.639 % at ``dx = 12`` because the area painting
spreads the same source over **3.1x** the tissue support (0.0268 → 0.0838), so
the exposure field is less peaked. The IFP solve moves in two directions at once
and both are expected: peak and mean pressure RISE (+4.110 %, +4.079 %) while the
peak Darcy speed FALLS (−6.891 %), because a smoother source builds a slightly
higher, flatter plateau with gentler gradients and ``v = -K grad p`` reads the
gradient. That speed drop is what re-pinned
``test_app_server.py::test_cfl_is_flagged_once_a_speed_has_been_measured``.

**No verdict moved.** Re-solved live on the post-two-route engine, both arms in
one process: zero of the nine published 3-D:2-D verdicts flip, worst ratio move
−2.603 %. ``OxygenConfig.supply`` was NOT re-tuned to absorb the −0.794 %, and
``calibrate_alpha()`` still returns 22.223477711172563 bit-exactly. Full
write-up: ``results/area_conservation_promotion.md``; superseded values under
``_superseded_2026-08-14_binary_vessel_painting``.
"""
from __future__ import annotations

import json
import math
import os
import platform
import sys

import numpy as np
import pytest

from tmesim.core import (
    CellPopulation,
    Domain,
    DrugParams,
    Microenvironment,
    PDEConfig,
    PKConfig,
    SimConfig,
    TissueFields,
)
from tmesim.physics.diffusion import simulate_drug_transport
from tmesim.sim.engine import _build_selftest_micro, simulate_response
from tmesim.treatment.drugs import delivery_source, plasma_profile

_REF_PATH = os.path.join(os.path.dirname(__file__), "data",
                         "static_source_reference.json")

# A deliberately weak drug so the kill fraction sits mid-curve and is genuinely
# sensitive to the exposure field (the library default saturates on this micro).
WEAK_DRUG = DrugParams(name="weak", diffusivity=200.0, ic50=4000.0, hill=1.2,
                       bbb_penetrance=0.15, clearance=3e-3, binding=0.0)

#: The exact cases captured in the stored reference, BEFORE the time-resolved
#: source existed. Any change to these definitions invalidates the reference.
_REFERENCE_CASES = [
    ("default_1dose", None,
     PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=1,
              dose_interval=3600.0 * 12),
     PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=200)),
    ("weak_1dose", WEAK_DRUG,
     PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=1,
              dose_interval=3600.0 * 12),
     PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=200)),
    ("weak_3dose", WEAK_DRUG,
     PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=3,
              dose_interval=3600.0 * 2),
     PDEConfig(dt=20.0, total_time=3600.0 * 6, convection=True, save_every=50)),
    ("weak_10dose_noconv", WEAK_DRUG,
     PKConfig(plasma_peak=20.0, half_life=1800.0, n_doses=10,
              dose_interval=1200.0),
     PDEConfig(dt=10.0, total_time=3600.0 * 4, convection=False, save_every=100)),
]

#: The three pinned cases whose PDE configuration is not in ``_REFERENCE_CASES``.
#: Hoisted here, rather than left inline in their tests, so that EVERY pinned
#: configuration is enumerable in one place — the cross-platform bound below is
#: derived from the longest time integration among them, and a case pinned with
#: a longer horizon or a finer dt must raise that bound rather than silently
#: leave it too strict.
_SOLVER_FLAT_PDE = PDEConfig(dt=10.0, total_time=1800.0, convection=False,
                             save_every=30)
_REAL_TISSUE_PDE = PDEConfig(dt=10.0, total_time=3600.0 * 12, convection=True,
                             save_every=100)
#: The exact production configuration: the same reconstruction call the landing
#: asset and every shipped result use. ``dt`` must stay at 10 s — at dx = 12 um
#: the advection CFL limit on this section is ~18.9 s, so dt = 20 raises.
_DX12_PDE = PDEConfig(dt=10.0, total_time=1800.0, convection=True, save_every=30)

_PINNED_PDE_CONFIGS = ([pde for _, _, _, pde in _REFERENCE_CASES]
                       + [_SOLVER_FLAT_PDE, _REAL_TISSUE_PDE, _DX12_PDE])


def _integration_steps(pde: PDEConfig) -> int:
    """Number of sequential time steps a case accumulates over.

    This is the ``n`` the cross-platform bound is derived from, so it is
    computed from the configuration rather than written down.
    """
    return int(round(float(pde.total_time) / float(pde.dt)))


@pytest.fixture(scope="module")
def reference():
    with open(_REF_PATH) as f:
        return json.load(f)


def _hex(x: str) -> float:
    """Reference values are stored as IEEE-754 hex so the comparison is exact."""
    return float.fromhex(x)


# --------------------------------------------------------------------------- #
# HOW A STORED REFERENCE IS COMPARED — two contracts, and which platform gets
# which. READ THIS BEFORE CONCLUDING THAT THE TOLERANCE BELOW HIDES A BUG.
# --------------------------------------------------------------------------- #
#: True on the machine every value in ``tests/data/static_source_reference.json``
#: was generated on: macOS / arm64, the reference environment named in the header
#: of ``constraints.txt``. Nothing else is inspected here — version drift is what
#: ``constraints.txt`` exists to remove, and a wrong-version run is meant to be a
#: loud failure under either contract, not a silently tolerated one.
_ON_REFERENCE_PLATFORM = (sys.platform == "darwin"
                          and platform.machine() == "arm64")

# WHY THERE ARE TWO CONTRACTS AT ALL.
#
# Bitwise equality is the RIGHT contract on the machine the reference was made on:
# it pins a change of one ULP, and the module docstring's 2026-08-12 entry is a
# case where a real physics fix moved these transport sums by only 1-3 ULP. Give
# that up and you give up the only guard that catches that class of change.
#
# Bitwise equality is an UNSATISFIABLE contract anywhere else, and not for a
# reason anyone can fix. Floating-point addition is not associative, so a
# different vectorised reduction width, a different BLAS thread count or an FMA
# contraction the other compiler did not emit changes the last few ULPs of every
# reduction. There is no bug to find and no version to pin: the numbers are
# physically the same result computed in a different order. So off the reference
# platform this module asserts a TIGHT RELATIVE BOUND instead, and the rest of
# this comment is how that bound was chosen.
#
# X — THE OBSERVED CROSS-PLATFORM NOISE. Two independent non-reference platforms,
# both running the constraints.txt versions (numpy 2.0.2 / scipy 1.13.1):
#
#   * The linux/amd64 GitHub runner, run 31932705658. Worst relative deviation
#     **8.662e-15** (48 ULP), over only FIVE field observations -- the old
#     assert-chain short-circuited at the first mismatching field, so that run
#     measured ``mean_auc`` four times and ``auc_sum`` once and NOTHING else.
#     This is why :class:`_ReferenceCheck` no longer short-circuits.
#   * A qemu-emulated linux/amd64 container (python:3.11-slim, the same install
#     line the workflow uses), which reaches every field: 25 of 59 differ, worst
#     **1.941e-14** (87 ULP), on ``solver_flat.mean_hist[6]``.
#
# Those two platforms do not agree with each OTHER bitwise either -- the emulated
# run reproduces none of the runner's five values exactly. So the emulation is a
# third platform, not a stand-in for the runner, and the worst observation over
# all 64 field-observations is 1.941e-14. Everything below the bound is
# measurement; the bound itself is not fitted to it.
#
# Y — THE BOUND, and why it is not "1e-12 because that looks safe". These
# quantities are time integrals: a sequential accumulation over ``n`` steps of
# strictly non-negative terms (concentrations), finished by a pairwise spatial
# reduction whose contribution is only ``O(log2 N)`` ULP. For a non-negative sum
# the classical worst-case round-off is ``(n-1)u`` with ``u = 2**-53``, so two
# different orderings of the same sum can differ by at most ``2(n-1)u``. The
# bound is that expression at the longest ``n`` this module pins -- today the
# 12 h window at ``dt = 10 s`` of
# ``test_real_tissue_static_path_matches_stored_reference_bitwise``, n = 4320:
#
#     _XPLATFORM_RTOL = 2 * 4320 * DBL_EPSILON = 1.918e-12
#
# with nothing rounded and nothing added. It is the worst case ARITHMETIC allows
# for the largest sum here, not a multiple of what was measured; it happens to
# sit 99x above the worst deviation observed on either platform, which is the
# cross-check rather than the derivation.
#
# ``n`` IS DERIVED, NOT WRITTEN DOWN. It is ``max`` over ``_PINNED_PDE_CONFIGS``,
# which enumerates every configuration this module pins a reference for, so a
# case added with a longer horizon or a finer dt raises the bound with it. Left
# as a literal, the bound would silently become too STRICT the day someone pins a
# longer integration -- it would then fail on a platform difference the theory
# permits, and the next person would widen it by hand without knowing why. And
# because a case could still be compared without going through that registry,
# :meth:`_ReferenceCheck.assert_matches` refuses outright to judge a case that
# integrates for longer than the bound was derived for.
#
# It is deliberately not a PER-CASE bound, which would be tighter (``solver_flat``
# integrates 180 steps -> 8.0e-14). Two reasons: the reordering model covers the
# final summation, not per-step FMA contraction inside the stencil or a library
# taking a different internal path, and the tightest per-case bound would leave
# ~4x margin over an observation we could only make under emulation. A single
# bound at the longest ``n`` is the conservative reading of the same theory.
#
# The argument depends on the solver not AMPLIFYING round-off, which it does not:
# explicit diffusion plus first-order decay at a CFL-stable ``dt`` is a
# contractive operator, so a perturbation introduced at step k decays by step n
# instead of growing. On an amplifying (chaotic, or unstably stepped) system this
# derivation would be invalid and no fixed relative bound would be defensible.
#
# Z — WHAT A REAL CHANGE LOOKS LIKE, and the honest limit of this bound. The two
# most recent deliberate changes to these very integrals, both recorded in the
# module docstring with their write-ups, moved them by **-0.794 %** (2026-08-14,
# dx = 12) and **+9.574 %** (2026-08-14, dx = 80); the 2026-08-13 two-route
# promotion moved ``kill_fraction`` by 2.76x and 5.29x. The smallest of those,
# 7.94e-3, is 4.1e9 times the bound. Anything of that kind fails here loudly.
#
# But the 2026-08-12 tumour-sigma fix moved the same sums by 1-3 ULP, i.e. by
# ~2e-16 to 7e-16 relative -- BELOW the cross-platform noise floor measured above.
# No tolerance can separate a change that size from a change of reduction order,
# so off the reference platform that class of change is undetectable by anyone, at
# any bound. That is a real hole and it is why the bitwise contract stays exactly
# as it was on macOS/arm64: the tolerant contract is a second, weaker instrument
# for a platform the strict one cannot run on, not a replacement for it.
#
# A NOTE ON THE FIELDS THAT NEVER MOVE. ``kill_fraction`` and ``viable_history``
# are counts over a fixed cell population (multiples of 1/n_cells) and
# ``plasma_mean`` is a pure-PK scalar with no spatial reduction; all three were
# bit-identical on both non-reference platforms. Because they are quantised, if
# one ever does move it moves by a whole quantum -- >= 1/300 = 3.3e-3 for the
# selftest micro, nine orders of magnitude above the bound -- so a genuine
# threshold flip fails this check rather than hiding inside it. That is the
# desired behaviour: a cell crossing the kill threshold because of a last-ULP
# difference is benign, but it is worth a human look, not a silent pass.

#: The longest sequential time integration among the pinned cases. Derived from
#: the configurations themselves (today: the real-tissue case, 12 h at
#: ``dt = 10 s``, n = 4320) so that pinning a longer one raises the bound too.
_LONGEST_TIME_INTEGRATION_STEPS = max(_integration_steps(p)
                                      for p in _PINNED_PDE_CONFIGS)

#: 2 * n * u -- see the derivation above. == 1.9185e-12 at n = 4320.
_XPLATFORM_RTOL = (2.0 * _LONGEST_TIME_INTEGRATION_STEPS
                   * sys.float_info.epsilon)


def _relative_deviation(got: float, want: float) -> float:
    """``|got - want| / |want|``, with the degenerate cases spelled out.

    0.0 for an exact match (including ``0.0 == -0.0``, which the bitwise
    contract accepts too), and ``inf`` for a mismatch against a zero reference
    or for a non-finite result — both of which must fail under either contract
    rather than divide their way to a small-looking number.
    """
    if got == want:
        return 0.0
    if want == 0.0 or not (math.isfinite(got) and math.isfinite(want)):
        return math.inf
    return abs(got - want) / abs(want)


def _record_deviations(case: str, rows) -> None:
    """Append this case's per-field deviations to ``$CODON_REFERENCE_DEVIATIONS``.

    The linux/amd64 job in ``.github/workflows/docker.yml`` is a MEASUREMENT of
    cross-platform float reproducibility, so it should produce the measurement
    when it PASSES, not only when it fails — otherwise the only way to learn the
    size of the divergence is to keep the job red, which is how it stayed red.

    Best-effort by construction: no failure in here may change whether a test
    passes. Unset variable, unwritable path, anything — the measurement is lost,
    the assertions are not.
    """
    path = os.environ.get("CODON_REFERENCE_DEVIATIONS")
    if not path:
        return
    try:
        with open(path, "a") as fh:
            for name, got, want, rel in rows:
                fh.write(json.dumps({
                    "case": case, "field": name, "rel": rel,
                    "got": float(got).hex(), "want": float(want).hex(),
                }) + "\n")
    except Exception:  # pragma: no cover - diagnostics must never fail a test
        pass


# DOES MY NEW BITWISE ASSERTION NEED THIS? Usually not. There are two kinds of
# bitwise assertion in this suite and only one of them has a platform problem:
#
#   1. Bitwise against a value computed HERE AND NOW -- two code paths, one
#      process, one machine. `assert a.tobytes() == b.tobytes()` where both
#      sides were just computed. This is platform-INDEPENDENT: whatever
#      reduction order the machine picks, it picks the same one twice, so the
#      assertion is satisfiable everywhere and must stay bitwise. Widening one
#      of these to a tolerance would destroy it for nothing.
#      Examples that correctly do NOT use `_ReferenceCheck`:
#      ``test_immune_sink.py`` (sink OFF must reproduce the engine byte for
#      byte) and ``test_mgmt_scalar.py::test_omitting_mgmt_is_bit_identical_to_
#      previous_behaviour`` (mgmt=None must not perturb one bit).
#
#   2. Bitwise against a value STORED FROM ANOTHER MACHINE -- a literal in
#      ``tests/data/``, generated once on macOS/arm64. This is satisfiable ONLY
#      on the machine that generated it, for the reasons above. It needs
#      `_ReferenceCheck`.
#
# The whole of category 2 in the linux/amd64 job's selected module set is the
# five tests in this file -- checked, not assumed. If you add a sixth, it goes
# through `_ReferenceCheck` and its PDEConfig goes in `_PINNED_PDE_CONFIGS`.
class _ReferenceCheck:
    """Every field of one stored-reference case, compared and reported together.

    WHY A COLLECTOR RATHER THAN A CHAIN OF BARE ``assert``\\ s. A chain
    short-circuits at the first mismatching field, so a run on a non-reference
    platform reports ONE number and hides every other field's deviation. That is
    not a hypothetical: run 31932705658 of the linux/amd64 job reported
    ``mean_auc`` for four cases and ``auc_sum`` for the solver, and told us
    nothing about ``penetration_index``, ``final_c_sum`` or ``viable_history``
    on that platform, because the tests died before reaching them — which is
    exactly the information needed to set the bound below. This compares every
    field first and names all the offenders in one failure.
    """

    def __init__(self, case: str, steps: int) -> None:
        #: ``steps`` is how many sequential time steps this case accumulates
        #: over — the ``n`` the bound is derived from. Required, not optional,
        #: so that pinning a new case forces its author to state it and the
        #: staleness check in :meth:`assert_matches` can do its job. 0 for a
        #: steady-state solve with no time integration at all.
        self.case = case
        self.steps = int(steps)
        self.rows = []  # (field, got, want, relative deviation)

    def field(self, name: str, got, want_hex: str) -> None:
        got = float(got)
        want = _hex(want_hex)
        self.rows.append((name, got, want, _relative_deviation(got, want)))

    def sequence(self, name: str, got_seq, want_hexes) -> None:
        """A whole history, element by element. Length is part of the contract:
        the list equality this replaces would also have failed on a short one."""
        got_seq = [float(v) for v in got_seq]
        assert len(got_seq) == len(want_hexes), (
            f"{self.case}: {name} has {len(got_seq)} entries, the stored "
            f"reference has {len(want_hexes)} — that is a shape change, not a "
            f"float difference, and no tolerance applies to it"
        )
        for i, (g, w) in enumerate(zip(got_seq, want_hexes)):
            self.field(f"{name}[{i}]", g, w)

    def assert_matches(self) -> None:
        # The bound is 2*n*u for the LONGEST accumulation in
        # ``_PINNED_PDE_CONFIGS``. A case that integrates for longer than that
        # is one the bound was not derived for, and applying it anyway would
        # silently be too strict. Refuse rather than judge — on every platform,
        # so it cannot lie in wait until someone runs this on linux.
        if self.steps > _LONGEST_TIME_INTEGRATION_STEPS:
            raise AssertionError(
                f"{self.case}: integrates {self.steps} time steps, but "
                f"_XPLATFORM_RTOL was derived for at most "
                f"{_LONGEST_TIME_INTEGRATION_STEPS} — so the bound in force is "
                "too STRICT for this case and would fail on a difference the "
                "theory permits. Add this case's PDEConfig to "
                "_PINNED_PDE_CONFIGS, which re-derives the bound. Do not raise "
                "the bound by hand."
            )

        worst = max((r[3] for r in self.rows), default=0.0)
        _record_deviations(self.case, self.rows)

        if _ON_REFERENCE_PLATFORM:
            bad = [r for r in self.rows if r[1] != r[2]]
            contract = ("BITWISE equality — this IS the platform "
                        "tests/data/static_source_reference.json was generated "
                        "on, so there is no float noise to allow for")
        else:
            bad = [r for r in self.rows if not r[3] <= _XPLATFORM_RTOL]
            contract = (f"relative deviation <= {_XPLATFORM_RTOL:g} — this is "
                        "NOT the platform the reference was generated on; see "
                        "the derivation of the bound in this module")
        if not bad:
            return

        lines = [
            f"{self.case}: the stored reference was not reproduced.",
            f"  platform : {sys.platform}/{platform.machine()}",
            f"  contract : {contract}",
            f"  {len(bad)} of {len(self.rows)} compared fields are outside it "
            f"(worst relative deviation across ALL fields: {worst:.6e}):",
        ]
        for name, got, want, rel in bad:
            lines.append(f"    {name:22s} got {got!r}  want {want!r}  "
                         f"rel {rel:.6e}")
        lines.append(
            "  A deviation far above the bound is a physics or wiring change, "
            "not float noise: find what moved and re-pin deliberately (the "
            "module docstring records how the last three re-pins were "
            "justified). Do NOT widen the bound to make this green."
        )
        raise AssertionError("\n".join(lines))


def _slab_micro(nx: int = 160, ny: int = 3, dx: float = 10.0):
    """A 1-D transport slab: a single vessel column on the left edge.

    Deliberately simple geometry so the perivascular-to-deep profile is a clean
    1-D decay and a "depth" statistic means what it says.
    """
    domain = Domain(nx=nx, ny=ny, dx=dx)
    vasc = np.zeros((ny, nx))
    vasc[:, 0] = 1.0
    fields = TissueFields(ecm_density=np.zeros((ny, nx)), vascular_density=vasc)
    cells = CellPopulation(
        x=np.array([nx * dx / 2]), y=np.array([ny * dx / 2]),
        cell_type=np.array([0]), state=np.array([0]),
        sensitivity=np.array([0.5]),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells)


def _relative_depth_um(auc_row: np.ndarray, dx: float, frac: float = 0.1) -> float:
    """Distance at which AUC has fallen to ``frac`` of its perivascular value.

    Normalised by the perivascular value, so this measures the SHAPE of the
    profile, not its amplitude — raising the dose cannot move it.
    """
    target = auc_row[0] * frac
    below = np.nonzero(auc_row < target)[0]
    assert below.size, "profile never falls below the depth threshold"
    i = int(below[0])
    a, b = float(auc_row[i - 1]), float(auc_row[i])
    return (i - 1 + (a - target) / (a - b)) * dx


def _matched_peak(target_mean: float, half_life: float, total_time: float,
                  dt: float) -> float:
    """Per-dose peak that makes the window-mean plasma equal ``target_mean``."""
    n = int(round(total_time / dt))
    unit = plasma_profile(
        PKConfig(plasma_peak=1.0, half_life=half_life, n_doses=1,
                 dose_interval=total_time),
        np.arange(n, dtype=float) * dt,
    ).mean()
    return float(target_mean) / float(unit)


# --------------------------------------------------------------------------- #
# (0) The comparison machinery itself
# --------------------------------------------------------------------------- #
def test_the_cross_platform_bound_is_the_derived_number():
    """The bound is a derivation, not a literal — so assert the derivation.

    Widening it now costs an argument, not a digit: this test fails, and the
    only honest way to make it pass is to change the reasoning in the comment
    above it. The two inequalities are the measurements that reasoning rests on:
    the worst cross-platform deviation ever observed (1.941e-14, qemu-emulated
    linux/amd64, ``solver_flat.mean_hist[6]``) and the smallest relative change
    a deliberate physics promotion has made to these same integrals (7.94e-3,
    the 2026-08-14 area-conserving vessel painting at dx = 12).
    """
    assert _LONGEST_TIME_INTEGRATION_STEPS == 4320
    assert _XPLATFORM_RTOL == 2.0 * 4320 * 2.0 ** -52
    assert 1.941e-14 < _XPLATFORM_RTOL < 7.94e-3
    # ...and 4320 is the real-tissue case's 12 h at dt = 10 s, reached through
    # the registry rather than written down, so pinning a longer case moves it.
    assert _LONGEST_TIME_INTEGRATION_STEPS == _integration_steps(_REAL_TISSUE_PDE)
    assert all(_integration_steps(p) <= _LONGEST_TIME_INTEGRATION_STEPS
               for p in _PINNED_PDE_CONFIGS)


def test_the_deviation_log_records_every_field_even_when_the_run_passes(
        tmp_path, monkeypatch):
    """``.github/workflows/docker.yml`` reads this file to build the summary
    table for its run, so its shape is a contract rather than a debug aid — and
    the point of it is that a PASSING run still says how close to the bound the
    platform sat."""
    log = tmp_path / "deviations.jsonl"
    monkeypatch.setenv("CODON_REFERENCE_DEVIATIONS", str(log))
    monkeypatch.setattr(sys.modules[__name__], "_ON_REFERENCE_PLATFORM", False)

    check = _ReferenceCheck("synthetic", steps=0)
    check.field("exact", 1.0, float(1.0).hex())
    check.field("one_ulp", math.nextafter(1.0, 2.0), float(1.0).hex())
    check.assert_matches()  # passes: one ULP is inside the bound

    rows = [json.loads(l) for l in log.read_text().splitlines() if l.strip()]
    assert [r["field"] for r in rows] == ["exact", "one_ulp"]
    assert set(rows[0]) == {"case", "field", "rel", "got", "want"}
    assert rows[0]["rel"] == 0.0
    assert 0.0 < rows[1]["rel"] < _XPLATFORM_RTOL


def test_a_case_longer_than_the_bound_was_derived_for_is_refused():
    """The staleness guard: a case that integrates longer than the bound's ``n``
    must not be judged by that bound, on any platform.

    Without this, adding a longer-horizon pin would leave the bound silently too
    STRICT — the failure mode where a guard drifts out from under its own
    justification and the next person widens it by hand.
    """
    check = _ReferenceCheck("synthetic_longer_than_the_bound",
                            steps=_LONGEST_TIME_INTEGRATION_STEPS + 1)
    check.field("exact", 1.0, float(1.0).hex())
    with pytest.raises(AssertionError, match="derived for at most"):
        check.assert_matches()


def test_both_reference_contracts_are_enforced_and_name_every_bad_field(
        monkeypatch):
    """One ULP must fail on the reference platform and pass off it; a deviation
    far above the bound must fail on both; and the failure must name EVERY bad
    field, which is the property the old assert-chain did not have."""
    mod = sys.modules[__name__]
    one = float(1.0).hex()
    # The instrument must not measure ITSELF: these deliberately-wrong values
    # would otherwise be appended to the deviation log and would then dominate
    # the "worst deviation this run" line in CI's summary, which is exactly the
    # number a reader trusts.
    monkeypatch.delenv("CODON_REFERENCE_DEVIATIONS", raising=False)

    def build():
        check = _ReferenceCheck("synthetic", steps=0)
        check.field("exact", 1.0, one)
        check.field("one_ulp", math.nextafter(1.0, 2.0), one)
        check.field("ten_x_bound", 1.0 + 10.0 * _XPLATFORM_RTOL, one)
        check.field("way_off", 2.0, one)
        return check

    monkeypatch.setattr(mod, "_ON_REFERENCE_PLATFORM", True)
    with pytest.raises(AssertionError) as bitwise:
        build().assert_matches()
    for field in ("one_ulp", "ten_x_bound", "way_off"):
        assert field in str(bitwise.value)
    assert "3 of 4" in str(bitwise.value)

    monkeypatch.setattr(mod, "_ON_REFERENCE_PLATFORM", False)
    with pytest.raises(AssertionError) as tolerant:
        build().assert_matches()
    # one ULP is exactly what the tolerant contract exists to absorb
    assert "one_ulp" not in str(tolerant.value)
    assert "ten_x_bound" in str(tolerant.value)
    assert "way_off" in str(tolerant.value)
    assert "2 of 4" in str(tolerant.value)

    # A shape change is not a float difference and no tolerance applies to it.
    for on_reference in (True, False):
        monkeypatch.setattr(mod, "_ON_REFERENCE_PLATFORM", on_reference)
        with pytest.raises(AssertionError, match="shape change"):
            _ReferenceCheck("synthetic", steps=0).sequence("h", [1.0, 2.0], [one])


# --------------------------------------------------------------------------- #
# (1-3) Backward compatibility — BIT-identical, not approximately equal
# --------------------------------------------------------------------------- #
def test_time_resolved_defaults_to_off():
    """The new capability must be opt-in, or every stored result in the repo
    silently changes meaning."""
    assert PDEConfig().time_resolved is False
    assert SimConfig().pde.time_resolved is False


@pytest.mark.parametrize("label,drug,pk,pde", _REFERENCE_CASES,
                         ids=[c[0] for c in _REFERENCE_CASES])
def test_static_path_matches_stored_reference_bitwise(reference, label, drug,
                                                      pk, pde):
    """The default (``time_resolved=False``) engine path reproduces values
    captured BEFORE the time-resolved source was added — bit for bit on the
    reference platform, within ``_XPLATFORM_RTOL`` anywhere else."""
    micro = _build_selftest_micro()
    kw = dict(pk=pk, pde=pde, seed=0)
    if drug is not None:
        kw["drug"] = drug
    res = simulate_response(micro, SimConfig(**kw))
    ref = reference[label]

    check = _ReferenceCheck(label, _integration_steps(pde))
    check.field("kill_fraction", res.kill_fraction, ref["kill_fraction"])
    check.field("mean_auc", res.extras["mean_auc"], ref["mean_auc"])
    check.field("plasma_mean", res.extras["plasma_mean"], ref["plasma_mean"])
    check.field("penetration_index", res.extras["penetration_index"],
                ref["penetration_index"])
    check.field("auc_sum", np.asarray(res.transport.auc).sum(), ref["auc_sum"])
    check.field("final_c_sum",
                np.asarray(res.transport.final_concentration).sum(),
                ref["final_c_sum"])
    check.sequence("viable_history", res.viable_history, ref["viable_history"])
    check.assert_matches()


def test_static_solver_matches_stored_reference_bitwise(reference,
                                                        flat_micro_factory):
    """Same pin one level down, on the raw solver."""
    micro = flat_micro_factory(ny=40, nx=40, dx=15.0, ecm_value=0.3)
    drug = DrugParams(name="tmz", diffusivity=300.0, clearance=3e-4, binding=1e-4)
    pde = _SOLVER_FLAT_PDE
    source = micro.domain.zeros()
    source[::12, ::12] = 1.0

    res = simulate_drug_transport(micro, drug, source, pde)
    ref = reference["solver_flat"]

    check = _ReferenceCheck("solver_flat",
                            _integration_steps(_SOLVER_FLAT_PDE))
    check.field("auc_sum", res.auc.sum(), ref["auc_sum"])
    check.field("auc_max", res.auc.max(), ref["auc_max"])
    check.field("final_c_sum", res.final_concentration.sum(), ref["final_c_sum"])
    check.field("penetration_index", res.penetration_index,
                ref["penetration_index"])
    check.sequence("mean_hist", res.mean_concentration_history,
                   ref["mean_hist"])
    check.assert_matches()


def test_real_tissue_static_path_matches_stored_reference_bitwise(reference,
                                                                  real_micro,
                                                                  real_drug):
    """And on the REAL Xenium section, through the whole engine."""
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=3,
                  dose_interval=3600.0 * 4)
    pde = _REAL_TISSUE_PDE
    res = simulate_response(real_micro, SimConfig(drug=real_drug, pk=pk, pde=pde,
                                                  seed=0))
    ref = reference["real_xenium_dx80_tmz_3dose"]

    check = _ReferenceCheck("real_xenium_dx80_tmz_3dose",
                            _integration_steps(_REAL_TISSUE_PDE))
    check.field("kill_fraction", res.kill_fraction, ref["kill_fraction"])
    check.field("mean_auc", res.extras["mean_auc"], ref["mean_auc"])
    check.field("auc_sum", np.asarray(res.transport.auc).sum(), ref["auc_sum"])
    check.field("final_c_sum",
                np.asarray(res.transport.final_concentration).sum(),
                ref["final_c_sum"])
    check.assert_matches()


# --------------------------------------------------------------------------- #
# (3b) The PRODUCTION grid — added 2026-08-12, see the module docstring
# --------------------------------------------------------------------------- #
#: ``_DX12_PDE`` — the exact production configuration — now lives with the other
#: pinned PDE configurations at the top of this module, because the
#: cross-platform bound is derived from all of them at once.
#:
#: The cached section the dx = 12 pins read. Named so the fixture's skip guard
#: and the fixture's own reconstruct call can never disagree about the path.
_DX12_SAMPLE = "data_cache/xenium_gbm"


@pytest.fixture(scope="module")
def production_micro_dx12():
    """The real Xenium GBM section on the PRODUCTION grid (dx = 12 um).

    Identical arguments to ``experiments/export_landing_asset.py``, so this pins
    the configuration the shipped IFP figures are actually computed on. Costs a
    couple of seconds to assemble; module-scoped so both dx = 12 pins share it.

    SKIPS, rather than errors, when the section's bytes are not cached. A
    missing dataset is not a numerics signal: without this guard both dx = 12
    pins ERRORED with ``IngestError`` on the cache-less linux/amd64 runner
    (``.github/workflows/docker.yml``), a job whose entire purpose is
    cross-platform float reproducibility. The check asks the real loader
    (:func:`tmesim.data.ingest.resolve_xenium_files`) about the exact path this
    fixture is about to read, so it skips only when the bytes are genuinely
    absent -- an ingest failure on a cache that IS present still errors.
    """
    from tmesim.data.ingest import resolve_xenium_files
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    if resolve_xenium_files(_DX12_SAMPLE).get("cells") is None:
        pytest.skip(f"no cached Xenium GBM section at {_DX12_SAMPLE!r} "
                    "(see `codon sections`) -- the dx = 12 production pins need "
                    "the real bytes")

    return reconstruct_real(source="xenium", sample=_DX12_SAMPLE,
                            dx_um=12.0, oxygen_mode="pde", state_mode="oxygen")


def test_production_dx12_ifp_solve_matches_stored_reference_bitwise(
        reference, production_micro_dx12):
    """The IFP solve itself, on the production grid, pinned bit for bit.

    This is the level the tumour-sigma units fix acts at, so it is the level the
    production guard belongs at: if ``_TUMOR_SIGMA_UM`` is ever changed, or the
    um -> pixel conversion is dropped and the constant goes back to meaning grid
    cells, these three numbers move immediately — no transport solve required to
    see it.
    """
    from tmesim.physics.pressure import solve_ifp

    p, velocity = solve_ifp(production_micro_dx12)
    mask = np.asarray(production_micro_dx12.domain.tissue_mask, dtype=bool)
    ref = reference["production_xenium_dx12_tmz_3dose"]

    check = _ReferenceCheck("production_xenium_dx12_tmz_3dose/ifp", steps=0)
    check.field("ifp_max_in_tissue_mmhg", p[mask].max(),
                ref["ifp_max_in_tissue_mmhg"])
    check.field("ifp_mean_in_tissue_mmhg", p[mask].mean(),
                ref["ifp_mean_in_tissue_mmhg"])
    check.field("ifp_speed_max_um_s",
                np.hypot(velocity[0], velocity[1])[mask].max(),
                ref["ifp_speed_max_um_s"])
    check.assert_matches()


def test_production_dx12_static_path_matches_stored_reference_bitwise(
        reference, production_micro_dx12):
    """End-to-end on the PRODUCTION grid, convection ON.

    The reference file used to pin only dx = 15 and dx = 80 — neither of which
    anything ships on. This closes that hole: the grid every published result in
    ``results/`` is computed on now has a bitwise pin too, running the same
    Darcy-convection path that the tumour-sigma fix touches.
    """
    from tmesim.molecular import library
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.validation.calibration import calibrated_drugparams

    drug = apply_barrier_to_drugparams(calibrated_drugparams("temozolomide"),
                                       library.get_molecule("temozolomide"))
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=3,
                  dose_interval=3600.0 * 4)
    res = simulate_response(
        production_micro_dx12,
        SimConfig(drug=drug, pk=pk, pde=_DX12_PDE, seed=0),
    )
    ref = reference["production_xenium_dx12_tmz_3dose"]

    check = _ReferenceCheck("production_xenium_dx12_tmz_3dose",
                            _integration_steps(_DX12_PDE))
    check.field("kill_fraction", res.kill_fraction, ref["kill_fraction"])
    check.field("mean_auc", res.extras["mean_auc"], ref["mean_auc"])
    check.field("plasma_mean", res.extras["plasma_mean"], ref["plasma_mean"])
    check.field("penetration_index", res.extras["penetration_index"],
                ref["penetration_index"])
    check.field("auc_sum", np.asarray(res.transport.auc).sum(), ref["auc_sum"])
    check.field("final_c_sum",
                np.asarray(res.transport.final_concentration).sum(),
                ref["final_c_sum"])
    check.sequence("viable_history", res.viable_history, ref["viable_history"])
    check.assert_matches()


def test_tumor_occupancy_smoothing_is_a_physical_length():
    """The units fix, asserted directly so it cannot silently regress.

    ``_TUMOR_SIGMA_UM`` is a MODEL parameter (it sets the vessel-leak gain and
    the lymphatic suppression), so its value must be a physical length and the
    pixel sigma handed to ``gaussian_filter`` must scale as ``1/dx``. Asserting
    the conversion rather than the constant means a deliberate recalibration of
    the 60 um is still allowed; reverting to grid-cell units is not.

    Also pins the fixed point that makes the whole change safe: at the
    production grid the conversion is EXACT in IEEE-754, so ``gaussian_filter``
    receives a bit-identical argument and no shipped result moves.
    """
    from tmesim.physics import pressure as _pressure

    assert not hasattr(_pressure, "_TUMOR_SIGMA"), (
        "the grid-cell tumour-sigma constant is back; the smoothing length is a "
        "physical length, not a pixel count"
    )
    sigma_um = _pressure._TUMOR_SIGMA_UM
    assert sigma_um > 0.0
    # dx = 12 um is the exact fixed point -> production is bit-identical.
    assert sigma_um / 12.0 == 5.0
    assert float.hex(sigma_um / 12.0) == float.hex(5.0)
    # The physical smoothing length is the SAME on every mesh.
    for dx in (8.0, 12.0, 15.0, 18.0, 24.0, 80.0):
        assert (sigma_um / dx) * dx == pytest.approx(sigma_um, rel=1e-12)


def test_unit_source_profile_is_bit_identical_to_none(flat_micro_factory):
    """``source_profile = ones`` must reproduce the static path EXACTLY (not just
    to a tolerance): ``dt * 1.0 == dt`` in IEEE-754, so the arithmetic is the
    same arithmetic. This is what makes the new branch provably inert."""
    micro = flat_micro_factory(ny=40, nx=40, dx=15.0, ecm_value=0.3)
    drug = DrugParams(name="tmz", diffusivity=300.0, clearance=3e-4, binding=1e-4)
    pde = _SOLVER_FLAT_PDE
    source = micro.domain.zeros()
    source[::12, ::12] = 1.0
    n_steps = int(round(pde.total_time / pde.dt))

    a = simulate_drug_transport(micro, drug, source, pde)
    b = simulate_drug_transport(micro, drug, source, pde,
                                source_profile=np.ones(n_steps))
    assert np.array_equal(a.auc, b.auc)
    assert np.array_equal(a.final_concentration, b.final_concentration)
    assert np.array_equal(a.mean_concentration_history,
                          b.mean_concentration_history)
    assert a.penetration_index == b.penetration_index


def test_source_profile_validation(flat_micro_factory):
    micro = flat_micro_factory(ny=16, nx=16, dx=15.0, ecm_value=0.3)
    drug = DrugParams(name="tmz", diffusivity=300.0, clearance=3e-4)
    pde = PDEConfig(dt=10.0, total_time=200.0, convection=False, save_every=5)
    source = micro.domain.ones()
    n_steps = int(round(pde.total_time / pde.dt))

    with pytest.raises(ValueError, match="at least one scale factor per step"):
        simulate_drug_transport(micro, drug, source, pde,
                                source_profile=np.ones(n_steps - 1))
    with pytest.raises(ValueError, match="1-D"):
        simulate_drug_transport(micro, drug, source, pde,
                                source_profile=np.ones((n_steps, 2)))
    with pytest.raises(ValueError, match="non-finite"):
        bad = np.ones(n_steps)
        bad[3] = np.nan
        simulate_drug_transport(micro, drug, source, pde, source_profile=bad)
    with pytest.raises(ValueError, match="non-negative"):
        neg = np.ones(n_steps)
        neg[3] = -0.5
        simulate_drug_transport(micro, drug, source, pde, source_profile=neg)
    # A LONGER profile is fine — only the first n_steps entries are consumed.
    long_p = np.concatenate([np.ones(n_steps), np.full(50, 1e6)])
    ok = simulate_drug_transport(micro, drug, source, pde, source_profile=long_p)
    base = simulate_drug_transport(micro, drug, source, pde)
    assert np.array_equal(ok.auc, base.auc)


# --------------------------------------------------------------------------- #
# (4) Convergence: constant plasma -> the static-source answer
# --------------------------------------------------------------------------- #
def test_time_resolved_converges_to_static_for_constant_plasma():
    """PHYSICS CHECK (b). Hold the plasma level constant and the time-resolved
    path must reproduce the static-source result. Verified two ways: with a
    literally constant profile (exact to round-off) and with a plasma half-life
    made long compared with the window, which is the physical limit the schedule
    sweep actually approaches."""
    micro = _slab_micro()
    drug = DrugParams(name="tl", diffusivity=100.0, clearance=1e-3, binding=0.0)
    dt, total_time = 5.0, 3600.0 * 6
    pde = PDEConfig(dt=dt, total_time=total_time, convection=False,
                    save_every=200)
    n_steps = int(round(total_time / dt))
    unit = delivery_source(micro, drug, 1.0)
    level = 10.0

    static = simulate_drug_transport(micro, drug, unit * level, pde)

    # (i) literally constant profile
    flat = simulate_drug_transport(micro, drug, unit, pde,
                                   source_profile=np.full(n_steps, level))
    assert np.allclose(flat.auc, static.auc, rtol=1e-12, atol=0.0)

    # (ii) the physical limit: half-life >> window, peak rescaled to matched
    #      cumulative exposure. Convergence must be MONOTONE in half-life.
    errs = []
    for half_life_h in (4.0, 48.0, 480.0, 4800.0):
        hl = half_life_h * 3600.0
        pk = PKConfig(plasma_peak=_matched_peak(level, hl, total_time, dt),
                      half_life=hl, n_doses=1, dose_interval=total_time)
        prof = plasma_profile(pk, np.arange(n_steps, dtype=float) * dt)
        assert abs(float(prof.mean()) - level) < 1e-9 * level  # matched exposure
        res = simulate_drug_transport(micro, drug, unit, pde, source_profile=prof)
        errs.append(abs(float(res.auc.sum()) / float(static.auc.sum()) - 1.0))

    assert errs == sorted(errs, reverse=True), (
        f"convergence to the static limit is not monotone in half-life: {errs}"
    )
    assert errs[-1] < 1e-4, (
        f"a 4800 h half-life over a 6 h window is effectively a constant plasma "
        f"level but still differs from the static source by {errs[-1]:.2e}"
    )


# --------------------------------------------------------------------------- #
# (5) Mass balance with a time-varying source
# --------------------------------------------------------------------------- #
def test_mass_balance_with_time_varying_source(inert_drug, flat_micro_factory):
    """PHYSICS CHECK (c). With clearance = binding = 0 and no-flux boundaries,
    nothing escapes, so the final tissue mass must equal exactly the mass
    injected by the (now time-varying) source: ``dt * sum_n p[n] * sum(S)``."""
    micro = flat_micro_factory(ny=40, nx=40, dx=15.0, ecm_value=0.3)
    pde = PDEConfig(dt=10.0, total_time=3600.0, convection=False, save_every=20)
    source = micro.domain.zeros()
    source[20, 20] = 2.0
    n_steps = int(round(pde.total_time / pde.dt))

    pk = PKConfig(plasma_peak=50.0, half_life=600.0, n_doses=4,
                  dose_interval=600.0)
    prof = plasma_profile(pk, np.arange(n_steps, dtype=float) * pde.dt)
    assert prof.max() > 1.5 * prof.min(), "profile must actually vary"

    res = simulate_drug_transport(micro, inert_drug, source, pde,
                                  source_profile=prof)
    injected = pde.dt * float(prof.sum()) * float(source.sum())
    final_mass = float(res.final_concentration.sum())

    assert abs(final_mass - injected) / injected < 1e-9
    assert res.final_concentration.min() >= -1e-9
    assert res.auc.min() >= -1e-9


def test_matched_cumulative_exposure_injects_identical_total(inert_drug,
                                                             flat_micro_factory):
    """The identity the whole "matched cumulative exposure" comparison rests on:
    a time-resolved schedule whose window-mean plasma equals ``plasma_mean``
    delivers exactly as much drug as the static source at ``plasma_mean``. If
    this ever drifted, every schedule comparison would be confounded by dose."""
    micro = flat_micro_factory(ny=30, nx=30, dx=15.0, ecm_value=0.2)
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 2, convection=False,
                    save_every=50)
    source = micro.domain.zeros()
    source[15, 15] = 1.0
    n_steps = int(round(pde.total_time / pde.dt))
    times = np.arange(n_steps, dtype=float) * pde.dt

    level = 7.5
    static = simulate_drug_transport(micro, inert_drug, source * level, pde)
    static_mass = float(static.final_concentration.sum())

    for n_doses, interval, half_life in ((1, pde.total_time, 3600.0),
                                         (3, 2400.0, 1800.0),
                                         (10, 720.0, 600.0)):
        peak = level / plasma_profile(
            PKConfig(plasma_peak=1.0, half_life=half_life, n_doses=n_doses,
                     dose_interval=interval), times).mean()
        prof = plasma_profile(
            PKConfig(plasma_peak=peak, half_life=half_life, n_doses=n_doses,
                     dose_interval=interval), times)
        res = simulate_drug_transport(micro, inert_drug, source, pde,
                                      source_profile=prof)
        mass = float(res.final_concentration.sum())
        assert abs(mass - static_mass) / static_mass < 1e-9, (
            f"n_doses={n_doses} delivered {mass:.6e} vs static {static_mass:.6e}"
        )


# --------------------------------------------------------------------------- #
# (6) Depth: the schedule changes the exposure field, in a derivable direction
# --------------------------------------------------------------------------- #
def test_earlier_delivery_reaches_deeper():
    """PHYSICS CHECK (a), with the direction ASSERTED.

    Naive expectation: at matched total exposure a long slow infusion reaches
    DEEPER than a short bolus. For this linear transport model that is provably
    the wrong way round (see the module docstring): the entire schedule effect is
    the terminal residual ``C(T)``, and the earlier the drug is given the smaller
    that residual is. So the BOLUS wins — in absolute exposure everywhere, and
    in relative depth too, because the un-counted term is a doubly-smoothed
    (flatter) field whose removal steepens what is left.

    This test asserts that direction on a transport-limited slab, and asserts the
    accompanying bound: the effect must SHRINK as the observation window grows
    past washout, which is what identifies it as a residual effect rather than a
    numerical artefact.
    """
    micro = _slab_micro()
    dx = 10.0
    # Transport-limited: L = sqrt(D/k) = 316 um against a 1600 um slab, and a
    # loss time-scale 1/k = 1000 s comparable to the short window below.
    drug = DrugParams(name="tl", diffusivity=100.0, clearance=1e-3, binding=0.0)
    dt = 5.0
    unit = delivery_source(micro, drug, 1.0)
    level = 10.0
    mid = micro.domain.shape[0] // 2

    def profile_result(half_life_h, total_time):
        n = int(round(total_time / dt))
        hl = half_life_h * 3600.0
        pk = PKConfig(plasma_peak=_matched_peak(level, hl, total_time, dt),
                      half_life=hl, n_doses=1, dose_interval=total_time)
        prof = plasma_profile(pk, np.arange(n, dtype=float) * dt)
        assert abs(float(prof.mean()) - level) < 1e-9 * level
        pde = PDEConfig(dt=dt, total_time=total_time, convection=False,
                        save_every=200)
        return simulate_drug_transport(micro, drug, unit, pde,
                                       source_profile=prof)

    # --- short window (2 h): comparable to washout, so schedule bites ---------
    short_T = 3600.0 * 2
    bolus = profile_result(0.05, short_T)      # 3 min half-life  -> all up front
    sustained = profile_result(48.0, short_T)  # 48 h half-life   -> ~flat infusion

    d_bolus = _relative_depth_um(bolus.auc[mid], dx)
    d_sustained = _relative_depth_um(sustained.auc[mid], dx)

    # 1. The schedule changes the answer AT ALL (this is the whole point).
    assert not np.allclose(bolus.auc, sustained.auc), (
        "matched-exposure schedules still produce an identical exposure field — "
        "the time-resolved source is not actually driving the solve"
    )
    # 2. Direction, on total exposure.
    assert bolus.auc.sum() > sustained.auc.sum() * 1.01
    # 3. Direction, on SHAPE — the depth at which exposure falls to 10% of its
    #    perivascular value, which is invariant to amplitude.
    assert d_bolus > d_sustained + 1.0, (
        f"bolus depth {d_bolus:.1f} um is not deeper than sustained "
        f"{d_sustained:.1f} um"
    )
    # 4. Same conclusion from the solver's own shape statistic.
    assert bolus.penetration_index > sustained.penetration_index
    # 5. The mechanism: the sustained infusion leaves drug in flight at T.
    assert (sustained.final_concentration.sum()
            > 10.0 * bolus.final_concentration.sum())

    # --- long window (96 h): past washout, so the effect must collapse --------
    long_T = 3600.0 * 96
    bolus_l = profile_result(0.05, long_T)
    sustained_l = profile_result(48.0, long_T)
    gap_short = bolus.auc.sum() / sustained.auc.sum() - 1.0
    gap_long = bolus_l.auc.sum() / sustained_l.auc.sum() - 1.0
    assert gap_long > 0.0
    assert gap_long < 0.1 * gap_short, (
        f"the schedule effect ({gap_long:.4f}) did not collapse relative to the "
        f"short window ({gap_short:.4f}) — it is not a terminal-residual effect"
    )


# --------------------------------------------------------------------------- #
# (7) Capability: the previously-zero spread is no longer zero
# --------------------------------------------------------------------------- #
def test_fractionation_at_matched_exposure_is_no_longer_degenerate(real_micro,
                                                                   real_drug):
    """The regression this whole change exists to fix, on the REAL section.

    With the static source, 1 / 3 / 10 fractions at matched cumulative plasma
    exposure were bit-identical inputs and the kill spread was exactly
    ``0.00e+00``. Assert (i) that the static path STILL behaves that way — the
    old behaviour is preserved bit-for-bit when the flag is off — and (ii) that
    the time-resolved path does not, and moves in the derivable direction: more
    fractions means later average delivery means a bigger terminal residual means
    less exposure realised inside the window.

    Run on ``real_micro`` rather than the 300-cell hand-built micro because kill
    fraction is quantised at 1/n_tumor and the schedule effect here is ~1%.
    """
    total_time, dt = 3600.0 * 12, 20.0
    half_life = 3600.0 * 4
    times = np.arange(int(round(total_time / dt)), dtype=float) * dt
    target_mean = plasma_profile(
        PKConfig(plasma_peak=300.0, half_life=half_life, n_doses=1,
                 dose_interval=total_time), times).mean()

    schedules = [(1, total_time), (3, 3600.0 * 4), (10, 3600.0 * 1.2)]
    for flag in (False, True):
        pde = PDEConfig(dt=dt, total_time=total_time, convection=True,
                        save_every=200, time_resolved=flag)
        kills, aucs = [], []
        for n, interval in schedules:
            unit = plasma_profile(
                PKConfig(plasma_peak=1.0, half_life=half_life, n_doses=n,
                         dose_interval=interval), times).mean()
            pk = PKConfig(plasma_peak=float(target_mean / unit),
                          half_life=half_life, n_doses=n, dose_interval=interval)
            res = simulate_response(real_micro,
                                    SimConfig(drug=real_drug, pk=pk, pde=pde,
                                              seed=0))
            kills.append(float(res.kill_fraction))
            aucs.append(float(np.asarray(res.transport.auc).sum()))

        if not flag:
            # Static source: matched cumulative exposure -> identical inputs.
            assert max(aucs) - min(aucs) == 0.0
            assert max(kills) - min(kills) == 0.0
        else:
            # Time-resolved: fractionation is now visible, and monotone in n.
            assert aucs == sorted(aucs, reverse=True), aucs
            assert max(aucs) - min(aucs) > 0.0
            assert kills == sorted(kills, reverse=True), kills
            assert max(kills) - min(kills) > 0.0, (
                "fractionation still produces a zero kill spread with "
                "time_resolved=True"
            )


def test_time_resolved_engine_reports_matched_total_exposure():
    """The engine must deliver the same total drug on both paths, so that a
    static-vs-resolved comparison is a comparison of SCHEDULE and not of dose."""
    micro = _build_selftest_micro()
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0 * 2, n_doses=4,
                  dose_interval=3600.0 * 2)
    kw = dict(drug=WEAK_DRUG, pk=pk, seed=0)
    base = PDEConfig(dt=20.0, total_time=3600.0 * 8, convection=True,
                     save_every=100)
    res_static = simulate_response(micro, SimConfig(pde=base, **kw))
    res_tr = simulate_response(
        micro,
        SimConfig(pde=PDEConfig(dt=base.dt, total_time=base.total_time,
                                convection=True, save_every=base.save_every,
                                time_resolved=True), **kw),
    )
    assert res_static.extras["time_resolved"] == 0.0
    assert res_tr.extras["time_resolved"] == 1.0
    # Same window-mean plasma reported on both paths -> same total delivered.
    assert res_tr.extras["plasma_delivered_mean"] == pytest.approx(
        res_static.extras["plasma_mean"], rel=1e-12
    )
    # ... but a genuinely different exposure field.
    assert res_tr.extras["mean_auc"] != res_static.extras["mean_auc"]


def test_time_resolved_is_deterministic():
    micro = _build_selftest_micro()
    pde = PDEConfig(dt=20.0, total_time=3600.0 * 4, convection=True,
                    save_every=100, time_resolved=True)
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0, n_doses=3,
                  dose_interval=3600.0)
    cfg = SimConfig(drug=WEAK_DRUG, pk=pk, pde=pde, seed=0)
    a = simulate_response(micro, cfg)
    b = simulate_response(micro, cfg)
    assert a.kill_fraction == b.kill_fraction
    assert np.array_equal(a.transport.auc, b.transport.auc)
