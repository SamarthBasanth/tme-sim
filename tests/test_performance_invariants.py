"""ALGORITHMIC performance invariants — the things that must not silently
regress between here and the next profile.

WHY THERE IS NO WALL-CLOCK ASSERTION IN THIS FILE
-------------------------------------------------
A wall-clock threshold is a measurement of the machine, not of the code. The
profile these tests were cut from (``results/performance_profile_2026-09-20.md``)
recorded the SAME spec, on the SAME commit, taking 40.7 s and 98.1 s in
``simulate`` inside one three-minute window, because the box was at load 14 on
10 CPUs. A ``assert elapsed < 60`` would have been green for the first and red
for the second, and would therefore have been reporting the neighbouring agent's
build rather than this repo's solver. That test is worse than no test: it is a
check that passes without testing its subject, which is the failure mode this
repo has already been bitten by seven times in one week.

So every assertion here is one of three kinds, and all three are EXACT and
load-independent:

1. **Counts.** How many times the sparse operator is factorised, how many
   triangular solves a step takes, how many times the vessel-distance transform
   is computed in one run. A loaded machine cannot change an integer.

2. **Dtypes.** ``float64`` where ``float64`` is load-bearing. A ``float32``
   downcast halves the memory and is the single most tempting "optimisation"
   on the transport path — and it is not an optimisation, it is a model change,
   because it changes every number the run reports. It must fail here rather
   than pass quietly and show up as a drifted ``kill_fraction`` six weeks later.

3. **Allocation scaling.** Peak traced bytes at grid size ``N`` versus ``4N``.
   Linear gives ~4x. Quadratic gives ~16x. The gap is wide enough that no
   amount of machine noise closes it, because ``tracemalloc`` counts bytes, not
   seconds. This is how "there is no O(n^2) in a path that should be O(n)" is
   pinned without a stopwatch.

MEASURED BASELINES (this file's numbers came from a run, not from a guess)
-------------------------------------------------------------------------
Taken on 2026-09-20, numpy 2.0.2, scipy per ``constraints.txt``:

* ``_build_diffusion_operator`` peak traced bytes / N: 276.2 (64x64),
  271.7 (128x128), 272.3 (256x256); 4x-grid ratios 3.93 and 4.01.
* ``_advective_divergence`` peak traced bytes / N: 49.3 (64x64), 48.7 (128x128).
* ``simulate_drug_transport`` peak traced bytes at 50 / 200 / 800 steps:
  682,268 / 681,180 / 681,180 — flat in ``n_steps`` to 3 significant figures.
* ``direct_vessel_distance_um`` peak traced bytes: 221,466 at (N=4096, M=500)
  and 220,013 at (N=4096, M=2000) — flat in M, which is the cKDTree doing its
  job. A pairwise ``(N, M)`` distance matrix would have gone up 4x.
* Operator nnz on an all-tissue grid: exactly ``5*ny*nx - 2*(ny+nx)``
  (460 / 2,900 / 7,840 at 10x10 / 20x30 / 40x40).
* One ``simulate_drug_transport`` call: 1 factorisation, ``n_steps`` triangular
  solves, ``n_steps`` advection calls with convection on and 0 with it off.

The thresholds below sit WIDE of these baselines on purpose. They are not
regression bounds on a constant factor — a 20% allocation increase is a code
review's business, not a test's. They are bounds on the COMPLEXITY CLASS, and
they are placed where only a change of class can cross them.
"""
from __future__ import annotations

import gc
import tracemalloc
from typing import Callable, Dict, Tuple

import numpy as np
import pytest

from tmesim.core import DrugParams, Domain, PDEConfig
from tmesim.physics import diffusion as diffusion_mod
from tmesim.physics.diffusion import (
    _advective_divergence,
    _build_diffusion_operator,
    _effective_diffusivity,
    simulate_drug_transport,
)
from tmesim.validation.vessel_profile import direct_vessel_distance_um

from .conftest import make_flat_micro

# The numpy data-allocation domain. numpy registers array buffers with
# tracemalloc under this domain; without it, array DATA is invisible to
# tracemalloc (it is malloc'd, not PyMem'd) and every scaling test below would
# be measuring Python object headers instead of the arrays — which is exactly
# the kind of check that passes without testing its subject.
#
# This is deliberately NOT wrapped in a skip. If a numpy upgrade removes the
# attribute, this module must go RED and be rewritten, not go green having
# measured nothing.
_NUMPY_TRACEMALLOC_DOMAIN = np.lib.tracemalloc_domain

INERT = DrugParams(name="inert", diffusivity=300.0, clearance=0.0, binding=0.0)


# --------------------------------------------------------------------------- #
# Measurement helpers (shared by the tests AND by the red-plant script, so a
# planted defect is scored by the same code the green run is)
# --------------------------------------------------------------------------- #
def peak_traced_bytes(fn: Callable, *args, **kwargs) -> Tuple[int, object]:
    """Peak bytes traced by ``tracemalloc`` across one call to ``fn``.

    ``gc.collect()`` first so a previous test's garbage is not charged to this
    call, and ``reset_peak()`` so the peak is this call's and not the
    session's. Returns ``(peak_bytes, return_value)``.
    """
    if tracemalloc.is_tracing():  # pragma: no cover - defensive
        tracemalloc.stop()
    gc.collect()
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        out = fn(*args, **kwargs)
        _cur, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return int(peak), out


def scaling_ratio(fn: Callable, small_kwargs: Dict, big_kwargs: Dict,
                  make_args: Callable) -> Tuple[float, int, int]:
    """Peak-bytes ratio between a 4x-larger problem and a baseline one.

    ``make_args(**kwargs)`` builds the call arguments so the construction cost
    is NOT charged to the measured call.
    """
    a_small = make_args(**small_kwargs)
    a_big = make_args(**big_kwargs)
    # One warm-up each: the first call through a scipy/numpy path can allocate
    # module-level lookup tables that belong to no particular problem size.
    fn(*a_small)
    fn(*a_big)
    peak_small, _ = peak_traced_bytes(fn, *a_small)
    peak_big, _ = peak_traced_bytes(fn, *a_big)
    return peak_big / peak_small, peak_small, peak_big


class TransportCallCounter:
    """Counts sparse factorisations, triangular solves and advection calls
    inside :func:`simulate_drug_transport`, by wrapping the two names the
    module resolves at call time.

    Adds timing-free wrappers and changes no arithmetic: the counted run
    returns the same arrays as an uncounted one, which
    ``test_counting_does_not_change_the_answer`` asserts rather than assumes.
    """

    def __init__(self) -> None:
        self.factorisations = 0
        self.solves = 0
        self.advections = 0
        self._real_factorized = None
        self._real_advect = None

    def __enter__(self) -> "TransportCallCounter":
        self._real_factorized = diffusion_mod.factorized
        self._real_advect = diffusion_mod._advective_divergence
        real_f = self._real_factorized
        real_a = self._real_advect

        def counted_factorized(a):
            self.factorisations += 1
            inner = real_f(a)

            def solve(b):
                self.solves += 1
                return inner(b)

            return solve

        def counted_advect(*a, **kw):
            self.advections += 1
            return real_a(*a, **kw)

        diffusion_mod.factorized = counted_factorized
        diffusion_mod._advective_divergence = counted_advect
        return self

    def __exit__(self, *exc) -> None:
        diffusion_mod.factorized = self._real_factorized
        diffusion_mod._advective_divergence = self._real_advect
        return None


def _transport_micro(ny: int = 24, nx: int = 24):
    micro = make_flat_micro(ny=ny, nx=nx, dx=15.0, ecm_value=0.3)
    source = micro.domain.zeros()
    source[ny // 2, nx // 2] = 2.0
    return micro, source


# --------------------------------------------------------------------------- #
# 1. COUNTS — one factorisation per solve, one solve per step
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("total_time,save_every", [(500.0, 10), (2000.0, 500)])
def test_exactly_one_factorisation_per_transport_solve(total_time, save_every):
    """The sparse LU is built ONCE and reused for every step.

    This is the single largest algorithmic property of the transport solve: the
    operator ``A = I - dt*L + dt*K`` does not depend on ``t``, so it is
    factorised once and back-substituted ``n_steps`` times. The measured cost
    on the flagship 286,200-voxel grid is 1.3 s for the factorisation against
    68.4 s for the 2,160 back-substitutions
    (results/performance_profile_2026-09-20.md §3). Moving ``factorized(a)``
    inside the loop — which is what any edit that makes ``k`` or ``d_eff``
    time-dependent would do — would turn 1.3 s into ~2,900 s and not change a
    single output, so nothing else in this suite would notice.
    """
    micro, source = _transport_micro()
    pde = PDEConfig(dt=10.0, total_time=total_time, convection=False,
                    save_every=save_every)
    n_steps = int(round(total_time / 10.0))

    with TransportCallCounter() as counter:
        simulate_drug_transport(micro, INERT, source, pde)

    assert counter.factorisations == 1, (
        f"{counter.factorisations} sparse factorisations for {n_steps} steps; "
        "the implicit operator is time-invariant and must be factorised once")
    assert counter.solves == n_steps, (
        f"{counter.solves} triangular solves for {n_steps} steps: the implicit "
        "step is ONE back-substitution, not an iteration")


def test_advection_cost_is_one_pass_per_step_and_zero_when_off():
    """Upwind advection is the second-largest term in the solve (15.8 % of
    ``simulate`` at dx = 12, results/performance_profile_2026-09-20.md §3), and
    it is EXPLICIT — one pass per step, and none at all when convection is off.

    The "and none at all" half is the one that matters: ``PDE_CONVECTION`` is a
    shipped constant, and a refactor that computes the divergence and then
    multiplies it by zero would cost 15.8 % of every run forever while every
    number stayed bit-identical.
    """
    ny = nx = 24
    micro, source = _transport_micro(ny, nx)
    velocity = np.zeros((2, ny, nx), dtype=np.float64)
    n_steps = 50
    pde_on = PDEConfig(dt=10.0, total_time=n_steps * 10.0, convection=True,
                       save_every=10)
    pde_off = PDEConfig(dt=10.0, total_time=n_steps * 10.0, convection=False,
                        save_every=10)

    with TransportCallCounter() as on:
        simulate_drug_transport(micro, INERT, source, pde_on, velocity=velocity)
    with TransportCallCounter() as off:
        simulate_drug_transport(micro, INERT, source, pde_off, velocity=None)

    assert on.advections == n_steps, (
        f"convection on: {on.advections} divergence passes for {n_steps} steps")
    assert off.advections == 0, (
        f"convection OFF still evaluated the advective divergence "
        f"{off.advections} times — that is pure cost with no effect")


def test_counting_does_not_change_the_answer():
    """The counter is instrumentation, not a model change.

    Planted deliberately: if ``TransportCallCounter`` ever perturbed the solve,
    every count above would be a count of a different program. The AUC field is
    compared bit-for-bit, not approximately.
    """
    micro, source = _transport_micro()
    pde = PDEConfig(dt=10.0, total_time=300.0, convection=False, save_every=10)

    plain = simulate_drug_transport(micro, INERT, source, pde)
    with TransportCallCounter():
        counted = simulate_drug_transport(micro, INERT, source, pde)

    assert np.array_equal(plain.auc, counted.auc), (
        "the call counter changed the AUC field; it is not instrumentation")
    assert plain.penetration_index == counted.penetration_index


def test_saved_history_is_o_steps_over_save_every_not_o_steps():
    """``save_every`` is what keeps the history from growing with the window.

    ``PDE_SAVE_EVERY`` is 500 and the shipped window is 2,160 steps, so the
    history is 6 samples and not 2,160 full-grid means. This pins that the
    knob is still connected: a regression that appends every step would be
    invisible in every output (the history is diagnostic) and would add an
    O(n_steps) allocation to a solve that currently has none.

    The exact count is ``1 + n_steps//save_every`` plus one more when the last
    step is not itself a save step, because the loop saves on
    ``step % save_every == 0 OR step == n_steps``. Spelling it the obvious way
    (``n_steps // save_every + 1``) is WRONG whenever the window does not
    divide evenly, which the shipped 2,160/500 does not — this test was written
    with the obvious spelling and went red on its first run, which is the only
    reason the formula here is the real one.
    """
    micro, source = _transport_micro()
    for total_time, save_every in ((2000.0, 500), (2000.0, 50), (8000.0, 500),
                                   (21600.0, 500)):
        pde = PDEConfig(dt=10.0, total_time=total_time, convection=False,
                        save_every=save_every)
        res = simulate_drug_transport(micro, INERT, source, pde)
        n_steps = int(round(total_time / 10.0))
        expected = (1 + n_steps // save_every
                    + (1 if n_steps % save_every else 0))
        assert res.mean_concentration_history.size == expected, (
            f"history has {res.mean_concentration_history.size} samples for "
            f"n_steps={n_steps}, save_every={save_every}; expected {expected}")
        assert res.times.size == res.mean_concentration_history.size
        # The point of the knob: the history must be O(n_steps/save_every), so
        # a 2,160-step window is single digits and not thousands.
        assert res.times.size <= n_steps // save_every + 2


def test_no_timelapse_allocates_no_frames():
    """``capture_frames=None`` must allocate NOTHING, not allocate and discard.

    24 frames of the flagship grid is 55 MB of float64 — more than a tenth of a
    512 MB instance — and the default server path asks for them. A run that
    does not ask must not pay.
    """
    micro, source = _transport_micro()
    pde = PDEConfig(dt=10.0, total_time=500.0, convection=False, save_every=10)
    res = simulate_drug_transport(micro, INERT, source, pde)
    assert res.frames is None
    assert res.frame_times is None
    assert res.frame_interval_s is None


# --------------------------------------------------------------------------- #
# 2. DTYPES — float64 is load-bearing, not a default
# --------------------------------------------------------------------------- #
def test_transport_outputs_are_float64():
    """Every array the transport solve hands downstream is float64.

    ``auc`` feeds ``delivery_metrics.compute_all``, ``compute_reach`` and the
    kill; ``final_concentration`` is a payload. float32 would halve the peak
    RSS of the solve — the most tempting single change available on the memory
    cliff — and would change ``kill_fraction`` in the 7th significant figure
    and the penetration length in the 3rd. That is a MODEL CHANGE and it must
    not arrive disguised as a memory fix.
    """
    micro, source = _transport_micro()
    pde = PDEConfig(dt=10.0, total_time=300.0, convection=False, save_every=10)
    res = simulate_drug_transport(micro, INERT, source, pde)

    for name in ("auc", "final_concentration", "times",
                 "mean_concentration_history"):
        arr = getattr(res, name)
        assert np.asarray(arr).dtype == np.float64, (
            f"TransportResult.{name} is {np.asarray(arr).dtype}, not float64")
    assert isinstance(res.penetration_index, float)


def test_operator_and_diffusivity_are_float64():
    """The operator's data and the effective-diffusivity field stay float64.

    A float32 operator would be factorised in single precision by SuperLU, and
    the implicit solve's conditioning is set by ``dt*D/dx^2`` — 2.1e2 on the
    shipped grid. There is no accuracy budget here to spend.
    """
    micro = make_flat_micro(ny=32, nx=32, dx=15.0, ecm_value=0.3)
    d_eff = _effective_diffusivity(INERT, micro.fields.ecm_density)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    operator = _build_diffusion_operator(d_eff, tissue, 15.0)

    assert d_eff.dtype == np.float64, f"d_eff is {d_eff.dtype}"
    assert operator.dtype == np.float64, f"operator is {operator.dtype}"
    assert operator.indices.dtype.kind == "i"


def test_vessel_distance_map_is_float64():
    """The distance transform is computed ONCE and shared by the penetration
    length and the reach map (``app/server/jobs.py`` §5). Both fit exponentials
    to it; float32 distances quantise the near-vessel bins the fit is most
    sensitive to."""
    domain = Domain(nx=32, ny=32, dx=15.0)
    rng = np.random.default_rng(0)
    endo = rng.uniform(0.0, 32 * 15.0, size=(64, 2))
    dist = direct_vessel_distance_um(domain, endo)
    assert dist.dtype == np.float64, f"distance map is {dist.dtype}"
    assert dist.shape == (32, 32)


# --------------------------------------------------------------------------- #
# 3. ALLOCATION SCALING — the complexity class, measured in bytes
# --------------------------------------------------------------------------- #
#: A 4x grid. Linear allocation gives ~4.0; quadratic gives ~16. The bound sits
#: at 8.0 — twice linear, half quadratic — so it cannot be crossed by a
#: constant-factor change and cannot be missed by a class change.
LINEAR_4X_CEILING = 8.0


def _operator_args(ny: int, nx: int):
    micro = make_flat_micro(ny=ny, nx=nx, dx=15.0, ecm_value=0.3)
    d_eff = _effective_diffusivity(INERT, micro.fields.ecm_density)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    return (d_eff, tissue, 15.0)


def test_diffusion_operator_allocation_is_linear_in_grid_size():
    """Assembling the 5-point operator is O(N) in memory, and must stay so.

    Measured at 276.2 / 271.7 / 272.3 bytes per voxel at 64x64 / 128x128 /
    256x256. A dense or banded-dense assembly — the natural thing to reach for
    if someone "simplifies" the COO triples — would be O(N^2) and would put the
    flagship grid's operator at 655 TB.
    """
    ratio, small, big = scaling_ratio(
        _build_diffusion_operator, {"ny": 64, "nx": 64}, {"ny": 128, "nx": 128},
        _operator_args)
    assert ratio < LINEAR_4X_CEILING, (
        f"operator assembly allocated {big:,} bytes at 4x the grid vs "
        f"{small:,} at 1x (ratio {ratio:.2f}); linear is ~4, quadratic is ~16")


def test_diffusion_operator_nnz_is_the_five_point_stencil_exactly():
    """On an all-tissue grid the operator has exactly ``5*N - 2*(ny+nx)``
    stored entries — N diagonals plus two per interior face.

    An exact identity, not a bound, because it is cheap to state exactly and
    every way of getting it wrong is interesting: a wider stencil changes the
    physics, a denser one changes the factorisation cost by orders of
    magnitude, and explicitly-stored structural zeros would inflate both.
    """
    for ny, nx in ((10, 10), (20, 30), (40, 40)):
        d_eff, tissue, dx = _operator_args(ny, nx)
        assert tissue.all(), "this identity is for an all-tissue grid"
        operator = _build_diffusion_operator(d_eff, tissue, dx)
        expected = 5 * ny * nx - 2 * (ny + nx)
        assert operator.nnz == expected, (
            f"{ny}x{nx}: nnz={operator.nnz}, five-point stencil predicts "
            f"{expected}")
        assert operator.nnz <= 5 * ny * nx


def _advect_args(ny: int, nx: int):
    return (np.zeros((ny, nx), dtype=np.float64),
            np.zeros((ny, nx), dtype=np.float64),
            np.zeros((ny, nx), dtype=np.float64),
            np.ones((ny, nx), dtype=bool),
            15.0)


def test_advective_divergence_allocation_is_linear_in_grid_size():
    """One advection pass allocates a FIXED NUMBER of grid-sized temporaries.

    It runs once per step — 2,160 times on the shipped window — so its
    allocation profile is multiplied by the step count in a way nothing else in
    the solve is. Measured at 49.3 / 48.7 bytes per voxel at 64x64 / 128x128.
    """
    ratio, small, big = scaling_ratio(
        _advective_divergence, {"ny": 64, "nx": 64}, {"ny": 128, "nx": 128},
        _advect_args)
    assert ratio < LINEAR_4X_CEILING, (
        f"advective divergence allocated {big:,} bytes at 4x the grid vs "
        f"{small:,} at 1x (ratio {ratio:.2f})")


def test_transport_peak_memory_does_not_grow_with_the_window():
    """The solve's memory is O(grid), NOT O(grid x steps).

    A 6 h window is 2,160 steps and a 24 h window (the protocol ceiling,
    ``TMESIM_PROTOCOL_MAX_WINDOW_S``) is 8,640. If peak memory scaled with the
    window, the 512 MB instance would be OOM-killed by the CLOCK rather than by
    the mesh, and the ``TMESIM_MAX_GRID_CELLS`` ceiling — which only knows
    about voxels — could not see it coming.

    Measured flat: 682,268 / 681,180 / 681,180 bytes at 50 / 200 / 800 steps.
    The ceiling is 1.5x across a 16x step increase, which no constant-factor
    change reaches and no per-step accumulation survives.
    """
    micro, source = _transport_micro(ny=48, nx=48)

    def run(total_time: float):
        pde = PDEConfig(dt=10.0, total_time=total_time, convection=False,
                        save_every=500)
        return simulate_drug_transport(micro, INERT, source, pde)

    run(500.0)  # warm-up
    peak_short, _ = peak_traced_bytes(run, 500.0)     # 50 steps
    peak_long, _ = peak_traced_bytes(run, 8000.0)     # 800 steps, 16x

    ratio = peak_long / peak_short
    assert ratio < 1.5, (
        f"peak traced bytes rose {ratio:.2f}x when the window grew 16x "
        f"({peak_short:,} -> {peak_long:,}); the solve is accumulating "
        "per-step state")


def _distance_args(ny: int, nx: int, n_endo: int):
    domain = Domain(nx=nx, ny=ny, dx=15.0)
    rng = np.random.default_rng(0)
    endo = rng.uniform(0.0, max(ny, nx) * 15.0, size=(n_endo, 2))
    return (domain, endo)


def test_vessel_distance_is_not_quadratic_in_cell_count():
    """``direct_vessel_distance_um`` must stay a cKDTree query, not a pairwise
    distance matrix.

    The flagship section has 286,200 grid centres and ~1,900 endothelial cells;
    the pairwise matrix is 4.4 GB and the tree query is 0.15 s. This is the
    single most inviting place in the repo to write an O(N x M) broadcast,
    because the broadcast spelling is one line and is CORRECT — it would pass
    every other test in this suite and only fail on a real section, on a real
    instance, as an OOM.

    Quadrupling M with N held fixed must not move the allocation: measured
    221,466 bytes at M=500 and 220,013 at M=2,000.
    """
    ratio, small, big = scaling_ratio(
        direct_vessel_distance_um,
        {"ny": 64, "nx": 64, "n_endo": 500},
        {"ny": 64, "nx": 64, "n_endo": 2000},
        _distance_args)
    assert ratio < 2.0, (
        f"quadrupling the endothelial-cell count moved the allocation "
        f"{small:,} -> {big:,} bytes (ratio {ratio:.2f}); an O(N*M) pairwise "
        "distance would have gone up ~4x")


def test_vessel_distance_is_linear_in_grid_size():
    """And linear in the grid it is evaluated on: 4x the voxels, ~4x the
    bytes (measured 221,466 -> 797,837)."""
    ratio, small, big = scaling_ratio(
        direct_vessel_distance_um,
        {"ny": 64, "nx": 64, "n_endo": 500},
        {"ny": 128, "nx": 128, "n_endo": 500},
        _distance_args)
    assert ratio < LINEAR_4X_CEILING, (
        f"distance map allocated {big:,} bytes at 4x the grid vs {small:,} "
        f"at 1x (ratio {ratio:.2f})")


# --------------------------------------------------------------------------- #
# 4. THE WHOLE-RUN INVARIANT — the one the server's own comment promises
# --------------------------------------------------------------------------- #
class _SilentCtx:
    """The three methods ``execute_run`` calls on its ctx, and nothing else."""

    def log(self, msg: str) -> None:
        return None

    def stage(self, name: str) -> None:
        return None

    def check_cancelled(self) -> None:
        return None


def test_one_run_factorises_once_and_maps_vessel_distance_once():
    """``app/server/jobs.py`` §5 states, in a comment, that the
    distance-from-nearest-endothelial-cell transform is "computed ONCE and
    handed to both, so adding the reach map does not double the cost of a run".

    That is a performance CLAIM and it had no test. It is also fragile in a
    specific way: ``penetration_length_fit`` computes the same map internally
    when it is handed nothing, so a future edit that stops threading
    ``distance_um=dist`` through ``compute_all`` or ``compute_reach`` would
    silently restore the double cost and change no number at all.

    Run at a coarse mesh and a short window: the invariant is a COUNT, so it
    does not need the shipped resolution to be true, and a test that needed
    47 s of solve to assert an integer would not be run.
    """
    pytest.importorskip("scipy")
    from app.server.cache import ReconstructionCache, make_key
    from app.server import jobs
    from tmesim.validation import vessel_profile

    cache = ReconstructionCache(maxsize=1)
    spec = {
        "section": "xenium_319",
        "dx_um": 80.0,
        "drug": "temozolomide",
        "sensitivity_mode": "uniform",
        "vascular_mode": "discrete",
        "oxygen_mode": "pde",
        "state_mode": "oxygen",
        "mgmt": None,
        "pk": {"plasma_peak": 50.0, "half_life_s": 14400.0,
               "n_doses": 1, "dose_interval_s": 86400.0},
        "total_time_s": 600.0,
        "timelapse_frames": 0,
    }
    try:
        cache.get_or_build(make_key(spec["section"], spec["dx_um"],
                                    spec["vascular_mode"], spec["oxygen_mode"],
                                    spec["state_mode"]))
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"section {spec['section']!r} unavailable here: {exc}")

    # THREE namespaces, because the function is reached through three of them
    # and patching only the definition site would count only some of the calls
    # — a counter that cannot see the second computation is exactly the
    # vacuous check this file exists to avoid. ``jobs.execute_run`` and
    # ``reach._distance_map`` import it inside the function (so the patch on
    # ``vessel_profile`` reaches them), while ``delivery_metrics`` binds it at
    # module scope and needs its own.
    from tmesim.validation import delivery_metrics

    calls = {"distance": 0}
    real_distance = vessel_profile.direct_vessel_distance_um

    def counted(*a, **kw):
        calls["distance"] += 1
        return real_distance(*a, **kw)

    patched = [(vessel_profile, "direct_vessel_distance_um"),
               (delivery_metrics, "direct_vessel_distance_um")]
    originals = [(mod, name, getattr(mod, name)) for mod, name in patched]
    for mod, name in patched:
        setattr(mod, name, counted)
    try:
        with TransportCallCounter() as counter:
            body = jobs.execute_run(spec, _SilentCtx(), cache=cache)
    finally:
        for mod, name, original in originals:
            setattr(mod, name, original)

    assert body["metrics"]["kill_fraction"] is not None
    assert calls["distance"] == 1, (
        f"one run computed the vessel-distance transform {calls['distance']} "
        "times; jobs.py promises exactly once, shared by the penetration "
        "length and the reach map")
    assert counter.factorisations == 1, (
        f"one run performed {counter.factorisations} sparse factorisations; "
        "one run is one transport solve is one factorisation")
