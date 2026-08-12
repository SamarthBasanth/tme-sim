"""Keystone integrator: full per-patient microenvironment drug-response simulation.

This is the wave-2 integrator that wires the (independently built and self-tested)
physics, cell, and treatment modules into ONE deterministic pipeline. Given a
reconstructed :class:`~tmesim.core.microenvironment.Microenvironment` and a
:class:`~tmesim.core.config.SimConfig`, :func:`simulate_response` predicts how
much of the tumor a drug kills.

Pipeline / wiring order
-----------------------
1. **Isolation.** Work on a private copy of the microenvironment so the caller's
   object is never mutated (the cell states and tissue fields are written to
   during the solve). See :func:`_working_copy`.
2. **Oxygen.** Steady-state O2 field (:func:`~tmesim.physics.oxygen.solve_oxygen`)
   stored into ``work.fields.oxygen``; then
   :func:`~tmesim.cells.agents.update_cell_states` bands every living cell into
   VIABLE / QUIESCENT / HYPOXIC (drives downstream drug resistance).
3. **Pressure.** Interstitial fluid pressure + Darcy velocity
   (:func:`~tmesim.physics.pressure.solve_ifp`) stored into ``work.fields.ifp`` /
   ``work.fields.velocity`` (the Jain-effect convective field).
4. **Delivery source.** The systemic PK profile
   (:func:`~tmesim.treatment.drugs.plasma_profile`) over the exposure window is
   turned into a spatial delivery source
   (:func:`~tmesim.treatment.drugs.delivery_source`) — either **time-averaged**
   into a static source (default) or driven **step by step** by ``plasma(t)``
   when ``cfg.pde.time_resolved`` is set. See "Plasma driver" below.
5. **Transport.** Convection-diffusion-reaction PDE
   (:func:`~tmesim.physics.diffusion.simulate_drug_transport`) integrates the
   drug field and returns the exposure ``auc`` (uM*s) at every grid point.
6. **Kill.** :func:`~tmesim.cells.pkpd.apply_drug_kill` maps each cell's local AUC
   through a Hill dose-response (modulated by sensitivity + hypoxia/quiescence
   resistance) and marks killed cells DEAD.
7. **Readouts.** Final viable tumor fraction, kill fraction, and the predicted
   response (see the response map below), plus a short viable-fraction trajectory.

Plasma driver: ``cfg.pde.time_resolved``
----------------------------------------
**Default (``time_resolved=False``) — Stage 0 plasma-averaging.** The multi-dose
plasma profile over the exposure window is collapsed to a single representative
**sustained** level, its time-average, and the solver is driven by the resulting
*static* source::

    plasma_mean = mean_t plasma_profile(pk, times)
    PS          = PS_SCALE_S_INV * vascular_density * bbb_penetrance   (1/s)
    source      = PS * plasma_mean                                     (uM/s)

This preserves the *total* delivered drug (mean level x window duration) and the
correct relative ranking of drugs/patients. Its cost is that the schedule is
erased before the solver ever sees it: any two schedules with the same cumulative
plasma exposure become **bit-identical inputs**, so fractionation and infusion
duration are unresolvable *by construction*, not as a finding.

**Where the ``1/s`` comes from.** The two spatial factors carry no units at all:
``vascular_density`` is a volume fraction and ``bbb_penetrance`` is a
dimensionless brain:plasma PARTITION coefficient, not a permeability and not a
rate (see ``core/config.py``). The whole ``1/s`` is carried by
:data:`tmesim.treatment.drugs.PS_SCALE_S_INV`, the transvascular rate scale, which
turns that dimensionless product into the first-order rate ``PS`` the transport
PDE needs. It defaults to **1.0** — exactly the implicit factor this engine used
before the constant was named, so every stored result is bit-identical — and its
calibrated value is **0.16122 1/s**, fit by :mod:`tmesim.validation.ps_scale`
against measured rodent brain PS. See
``results/transvascular_exchange_revalidation.md`` section 6.2.

**``time_resolved=True`` — the real plasma curve drives the solve.** The unit
delivery field is built once and rescaled by ``plasma(t_n)`` at every step::

    PS            = PS_SCALE_S_INV * vascular_density * bbb_penetrance   (1/s)
    S(x, y, t_n)  = plasma_profile(pk, t_n) * PS(x, y)                   (uM/s)

``PS`` is the same field as above — a rate, i.e. uM/s per uM of plasma — and it is
exactly what the unit delivery field holds, so ``PS_SCALE_S_INV`` multiplies it
once and switching plasma drivers cannot change the units of the source.

Tissue concentration can then genuinely rise and fall, so *when* drug is given
changes the answer. It costs essentially nothing: the implicit operator
``A = I - dt*L + dt*K`` does not involve the source, so the LU factorization is
still built once and reused (see :mod:`tmesim.physics.diffusion`).

Both paths deliver the SAME total drug over the window (the profile is sampled on
the solver's own step grid, whose mean is exactly ``plasma_mean``), which is what
makes "matched cumulative exposure" comparisons between them meaningful. What
differs is only the residual drug still in flight at the end of the window —
which is precisely the channel through which schedule acts in a linear transport
model. In the limit of a constant plasma level the two paths coincide.

``time_resolved`` defaults to **False** so every stored result, calibration and
regression in the repo reproduces bit-for-bit; a regression test pins that
against a stored reference (``tests/data/static_source_reference.json``).

predicted_response map
----------------------
``predicted_response`` is defined explicitly as the **fraction of tumor cells
killed**::

    viable_tumor_fraction = alive_tumor / original_tumor        (in [0, 1])
    kill_fraction         = 1 - viable_tumor_fraction           (in [0, 1])
    predicted_response    = kill_fraction                       (in [0, 1])

So ``predicted_response = 1`` means a complete kill (maximal response) and ``0``
means no tumor cell was killed. This is the number the validation harness scores
against each patient's observed response.

Determinism
-----------
Every stochastic step is driven by ``np.random.default_rng(cfg.seed)`` (or an
explicitly supplied generator). Given the same ``cfg.seed`` the result is
bit-for-bit reproducible.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from ..cells.agents import update_cell_states
from ..cells.pkpd import apply_drug_kill
from ..core import (
    DrugParams,
    Microenvironment,
    SimConfig,
    SimulationResult,
    TissueFields,
    TransportResult,
)
from ..physics.diffusion import simulate_drug_transport
from ..physics.oxygen import solve_oxygen
from ..physics.pressure import solve_ifp
from ..treatment.drugs import delivery_source, exchange_coefficient, get_drug

# Number of points in the diagnostic viable-fraction trajectory.
_N_HISTORY = 6


def _working_copy(micro: Microenvironment) -> Microenvironment:
    """Return a private copy of ``micro`` safe to mutate during the solve.

    The domain (grid + tissue mask) is immutable input and is shared by
    reference. The cell population is deep-copied (it is relabeled and killed in
    place) and the tissue fields are placed in a fresh :class:`TissueFields`
    container so the oxygen / IFP / velocity fields we fill do not leak back into
    the caller's object. The ECM and vascular-density arrays are read-only inputs
    and may be shared by reference.
    """
    fields = micro.fields
    work_fields = TissueFields(
        ecm_density=fields.ecm_density,
        vascular_density=fields.vascular_density,
        oxygen=None,
        ifp=None,
        velocity=None,
    )
    return Microenvironment(
        domain=micro.domain,
        fields=work_fields,
        cells=micro.cells.copy(),
        patient_id=micro.patient_id,
        metadata=micro.metadata,
    )


def _resolve_drug(cfg: SimConfig) -> DrugParams:
    """Resolve the drug to simulate from ``cfg.drug`` (a :class:`DrugParams` or a
    library name string)."""
    drug = cfg.drug
    if isinstance(drug, str):
        return get_drug(drug)
    return drug


def _viable_history(
    labeled_cells,
    micro: Microenvironment,
    transport: TransportResult,
    cfg: SimConfig,
    drug: DrugParams,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a deterministic viable-tumor-fraction trajectory over the exposure.

    We reconstruct the cumulative kill at a handful of exposure fractions
    ``f in linspace(0, 1, K)`` by scaling the transport exposure ``auc * f`` and
    re-running the exact Hill kill mapping on a FRESH copy of the (oxygen-labeled,
    pre-kill) cells each time — so there is no double-killing.

    To make the trajectory a *coherent* cumulative-kill curve (monotonically
    non-increasing) rather than independent noisy snapshots, every fraction is
    thresholded against a SINGLE shared per-cell uniform draw ``u``: a cell is
    considered dead by fraction ``f`` iff ``u < p_kill(f)``. Since the Hill kill
    probability increases with exposure, the dead set only grows as ``f`` rises.
    The kill probabilities themselves come straight from
    :func:`~tmesim.cells.pkpd.apply_drug_kill` (its return value), so the pkpd
    mapping is not re-implemented here. ``f`` maps linearly to wall-clock time.

    Returns ``(times, viable_history)``, both length ``K``, with
    ``viable_history[0] == 1.0`` (no exposure) and ``viable_history[-1]`` the
    full-exposure viable fraction.
    """
    window = float(transport.times[-1]) if transport.times.size else 0.0
    fractions = np.linspace(0.0, 1.0, _N_HISTORY)
    times = window * fractions
    viable_history = np.empty(_N_HISTORY, dtype=np.float64)

    tumor = labeled_cells.tumor_mask()
    n_tumor = int(tumor.sum())
    # A single shared per-cell threshold makes the curve monotone by construction.
    draws = np.random.default_rng(int(cfg.seed) + 1_000).random(labeled_cells.n)

    base_auc = np.asarray(transport.auc, dtype=np.float64)
    for i, f in enumerate(fractions):
        if n_tumor == 0:
            viable_history[i] = 0.0
            continue
        scaled = TransportResult(
            auc=base_auc * float(f),
            final_concentration=transport.final_concentration,
            times=transport.times,
            mean_concentration_history=transport.mean_concentration_history,
            penetration_index=transport.penetration_index,
        )
        snapshot = Microenvironment(
            domain=micro.domain,
            fields=micro.fields,
            cells=labeled_cells.copy(),
            patient_id=micro.patient_id,
            metadata=micro.metadata,
        )
        # Reuse pkpd's exact Hill mapping; ignore its internal draw and instead
        # threshold against the shared per-cell draw for a coherent curve.
        p_kill = apply_drug_kill(snapshot, scaled, cfg.pd, drug, rng=rng_dummy())
        killed_tumor = int((tumor & (draws < p_kill)).sum())
        viable_history[i] = (n_tumor - killed_tumor) / n_tumor

    return times, viable_history


def rng_dummy() -> np.random.Generator:
    """A throwaway generator for :func:`apply_drug_kill` calls whose stochastic
    outcome we discard (we only consume the returned per-cell probabilities)."""
    return np.random.default_rng(0)


def simulate_response(
    micro: Microenvironment,
    cfg: SimConfig,
    rng: Optional[np.random.Generator] = None,
    attribute: bool = False,
    capture_frames: Optional[int] = None,
    thresholds: Optional["MechanismThresholds"] = None,
) -> SimulationResult:
    """Run the full microenvironment drug-response simulation for one patient.

    Executes the oxygen -> pressure -> delivery -> transport -> kill pipeline
    (see the module docstring for the wiring order, the plasma-averaging
    simplification, and the ``predicted_response`` map) on a private copy of
    ``micro`` — the caller's microenvironment is never mutated.

    Parameters
    ----------
    micro : Microenvironment
        The reconstructed patient microenvironment (grid, tissue fields, cells).
        Not modified.
    cfg : SimConfig
        Bundled per-stage configuration: ``drug`` (a :class:`DrugParams` or a
        drug-library name), ``pk``, ``pde``, ``oxygen``, ``pd``, and ``seed``.
    rng : np.random.Generator, optional
        RNG for the stochastic kill draws. Defaults to
        ``np.random.default_rng(cfg.seed)`` for reproducibility.
    attribute : bool, default False
        Compute the per-cell **mechanism attribution** (see
        :mod:`tmesim.sim.mechanism`) and return it as ``result.mechanism``. It is
        built from the arrays :func:`~tmesim.cells.pkpd.apply_drug_kill` actually
        decided on, consumes no extra RNG draw, and changes no arithmetic — the
        run is bit-identical with it on or off, which is why it defaults to off
        rather than being made unconditional. Raises
        :class:`~tmesim.sim.mechanism.MechanismReconciliationError` if the labels
        do not reproduce ``kill_fraction`` bit-for-bit.
    capture_frames : int, optional
        Number of drug-concentration snapshots to keep through the solve (the
        **time-lapse**), passed straight to
        :func:`~tmesim.physics.diffusion.simulate_drug_transport`. ``None``
        keeps none. To see the drug actually clear rather than plateau, combine
        this with ``cfg.pde.time_resolved = True``.
    thresholds : MechanismThresholds, optional
        Bands for the attribution. ``None`` uses the defaults, whose two
        otherwise-arbitrary values are derived from ``cfg.pd``.

    Returns
    -------
    SimulationResult
        ``predicted_response`` (= kill fraction, in [0, 1]),
        ``viable_tumor_fraction``, ``kill_fraction``, a short ``viable_history``
        trajectory with matching ``times``, the driving ``transport`` result,
        ``fields`` for plotting (oxygen / ifp / auc / vascular_density),
        diagnostic ``extras``, and — when ``attribute=True`` — ``mechanism``.
    """
    if rng is None:
        rng = np.random.default_rng(cfg.seed)

    # 1. Isolation: never mutate the caller's microenvironment.
    work = _working_copy(micro)

    # 2. Oxygen field + hypoxia/quiescence banding.
    oxygen = solve_oxygen(work, cfg.oxygen)
    work.fields.oxygen = oxygen
    update_cell_states(work, oxygen, cfg.oxygen)

    # Snapshot the oxygen-labeled, pre-kill cells for the history reconstruction.
    labeled_cells = work.cells.copy()

    # 3. Interstitial fluid pressure + convective velocity (Jain effect).
    ifp, velocity = solve_ifp(work)
    work.fields.ifp = ifp
    work.fields.velocity = velocity

    # 4. Drug + delivery source.
    #    cfg.pde.time_resolved is False  -> Stage 0 plasma-averaging (default,
    #                                       bit-identical to the historical path)
    #    cfg.pde.time_resolved is True   -> the real plasma(t) drives the solver
    drug = _resolve_drug(cfg)
    times = np.arange(0.0, cfg.pde.total_time, cfg.pde.dt)
    plasma = plasma_profile_safe(cfg, times)
    plasma_mean = float(plasma.mean()) if plasma.size else 0.0

    time_resolved = bool(getattr(cfg.pde, "time_resolved", False))
    if time_resolved:
        # Unit delivery field; the plasma level is reapplied step by step.
        source = delivery_source(work, drug, 1.0)
        source_profile = _step_plasma_profile(cfg, plasma)
        plasma_delivered_mean = (
            float(source_profile.mean()) if source_profile.size else 0.0
        )
    else:
        source = delivery_source(work, drug, plasma_mean)
        source_profile = None
        plasma_delivered_mean = plasma_mean

    # 4b. Optional transvascular BACK-FLUX. The source above is one-way; with
    #     this on, the same permeability field is also subtracted as a
    #     first-order loss, which is algebraically S = PS*(Cp - C) and bounds
    #     tissue concentration by the plasma level. Same PS the source used --
    #     no new parameter. Default off keeps every stored result bit-exact.
    #     See PDEConfig.transvascular_exchange and results/breast_kp_validation.md.
    exchange_k = (exchange_coefficient(work, drug)
                  if bool(getattr(cfg.pde, "transvascular_exchange", False))
                  else None)

    # 5. Transport PDE -> exposure (AUC) field.
    transport = simulate_drug_transport(
        work,
        drug,
        source,
        cfg.pde,
        velocity=(velocity if cfg.pde.convection else None),
        source_profile=source_profile,
        capture_frames=capture_frames,
        exchange_k=exchange_k,
    )

    # 6. Pharmacodynamic kill over the treatment course.
    #    When an attribution is asked for, hand apply_drug_kill a dict to fill
    #    with the arrays it decides on, so the attribution reads the kill's own
    #    numbers instead of re-deriving (and possibly drifting from) them.
    kill_diag: Optional[dict] = {} if attribute else None
    apply_drug_kill(work, transport, cfg.pd, drug, rng=rng, diagnostics=kill_diag)

    # 7. Readouts.
    viable = work.cells.viable_tumor_fraction()
    kill_fraction = 1.0 - viable
    predicted_response = kill_fraction  # documented map: response == fraction killed

    # 8. Diagnostic viable-fraction trajectory.
    hist_times, viable_history = _viable_history(
        labeled_cells, work, transport, cfg, drug
    )

    # 9. Package.
    tissue = np.asarray(work.domain.tissue_mask, dtype=bool)
    mean_oxygen = float(oxygen[tissue].mean()) if tissue.any() else 0.0
    tumor = work.cells.tumor_mask()
    n_tumor = int(tumor.sum())
    hypoxic_tumor_frac = (
        float((tumor & (labeled_cells.state == 1)).sum()) / n_tumor
        if n_tumor
        else 0.0
    )

    # 10. Mechanism attribution (opt-in; asserts its own faithfulness).
    mechanism = None
    if kill_diag is not None:
        from .mechanism import attribute as _attribute

        mechanism = _attribute(
            work.cells.cell_type,
            kill_diag,
            kill_fraction,
            pd=cfg.pd,
            thresholds=thresholds,
        )

    return SimulationResult(
        predicted_response=float(predicted_response),
        viable_tumor_fraction=float(viable),
        kill_fraction=float(kill_fraction),
        viable_history=viable_history,
        times=hist_times,
        transport=transport,
        fields={
            "oxygen": oxygen,
            "ifp": ifp,
            "auc": transport.auc,
            "vascular_density": work.fields.vascular_density,
        },
        extras={
            "penetration_index": float(transport.penetration_index),
            "plasma_mean": plasma_mean,
            "plasma_delivered_mean": float(plasma_delivered_mean),
            "plasma_peak": float(plasma.max()) if plasma.size else 0.0,
            "time_resolved": float(time_resolved),
            "mean_auc": float(transport.auc[tissue].mean()) if tissue.any() else 0.0,
            "mean_oxygen": mean_oxygen,
            "hypoxic_tumor_fraction": hypoxic_tumor_frac,
            "bbb_penetrance": float(drug.bbb_penetrance),
            "n_tumor": float(n_tumor),
        },
        mechanism=mechanism,
    )


def _step_plasma_profile(cfg: SimConfig, plasma: np.ndarray) -> np.ndarray:
    """Plasma concentration (uM) sampled at the solver's own step times.

    The transport solver takes ``n_steps = round(total_time / dt)`` steps and
    applies the source at the OLD time level of each step, i.e. at
    ``t_n = n*dt`` for ``n = 0 .. n_steps-1``. This helper returns
    ``plasma_profile(cfg.pk, t_n)`` on exactly that grid.

    ``plasma`` is the profile the engine has already evaluated on
    ``np.arange(0, total_time, dt)``, which *is* that grid whenever
    ``total_time`` is an exact multiple of ``dt`` (the normal case) — it is then
    reused unchanged, so the delivered window-mean equals the ``plasma_mean``
    the static path uses and the two paths deliver the SAME total drug. Only if
    floating-point ``arange`` produced a different count do we re-evaluate on the
    exact step grid.
    """
    n_steps = int(round(float(cfg.pde.total_time) / float(cfg.pde.dt)))
    if plasma.size == n_steps:
        return plasma
    step_times = np.arange(n_steps, dtype=np.float64) * float(cfg.pde.dt)
    return plasma_profile_safe(cfg, step_times)


def plasma_profile_safe(cfg: SimConfig, times: np.ndarray) -> np.ndarray:
    """Evaluate the systemic plasma profile, importing lazily to avoid a hard
    module-load dependency cycle. Thin wrapper over
    :func:`tmesim.treatment.drugs.plasma_profile`."""
    from ..treatment.drugs import plasma_profile

    return plasma_profile(cfg.pk, times)


def _build_selftest_micro() -> Microenvironment:
    """A tiny hand-built Microenvironment (no data generator) for the self-test.

    A central dense tumor mass on a poorly-perfused core, ringed by endothelial
    cells feeding a peripheral vascular field, with denser ECM in the core — just
    enough real structure to exercise oxygen -> pressure -> transport -> kill.
    """
    from ..core import (
        CellPopulation,
        CellState,
        CellType,
        Domain,
        TissueFields,
    )

    rng = np.random.default_rng(0)
    ny = nx = 48
    dx = 15.0
    domain = Domain(nx=nx, ny=ny, dx=dx)
    X, Y = domain.coords_um()
    cx, cy = X.mean(), Y.mean()
    R = np.hypot(X - cx, Y - cy)
    core_r = 0.35 * R.max()

    # Vascular density: sparse in the core, richer at the periphery.
    vascular = np.clip(0.7 * (1.0 - np.exp(-(R / (0.9 * core_r)) ** 2)) + 0.1, 0.0, 1.0)
    # ECM: denser in the core.
    ecm = np.clip(0.3 + 0.5 * np.exp(-(R / (1.0 * core_r)) ** 2), 0.0, 1.0)
    fields = TissueFields(ecm_density=ecm, vascular_density=vascular)

    # Cells: tumor mass in the core, endothelial cells at the periphery.
    n_tumor, n_endo = 300, 80
    ang = rng.uniform(0, 2 * np.pi, n_tumor)
    rad = np.abs(rng.normal(0, 0.5 * core_r, n_tumor))
    tx = np.clip(cx + rad * np.cos(ang), 0, domain.extent_um[0] - 1e-6)
    ty = np.clip(cy + rad * np.sin(ang), 0, domain.extent_um[1] - 1e-6)
    ex = rng.uniform(0, domain.extent_um[0] - 1e-6, n_endo)
    ey = rng.uniform(0, domain.extent_um[1] - 1e-6, n_endo)

    x = np.concatenate([tx, ex])
    y = np.concatenate([ty, ey])
    cell_type = np.concatenate([
        np.full(n_tumor, int(CellType.TUMOR), dtype=np.int8),
        np.full(n_endo, int(CellType.ENDOTHELIAL), dtype=np.int8),
    ])
    state = np.full(x.size, int(CellState.VIABLE), dtype=np.int8)
    sensitivity = np.clip(0.5 + 0.1 * rng.standard_normal(x.size), 0.0, 1.0)
    cells = CellPopulation(x=x, y=y, cell_type=cell_type, state=state,
                           sensitivity=sensitivity)
    return Microenvironment(domain=domain, fields=fields, cells=cells,
                            patient_id="ENGINE-SELFTEST")


if __name__ == "__main__":
    # -----------------------------------------------------------------------
    # Self-test: run one tiny hand-built microenvironment end to end and assert
    # the integrator is well-formed, deterministic, and that stronger BBB
    # penetrance (better delivery) yields at least as much kill as weak
    # penetrance — the sanity check that transport is actually coupled into the
    # kill. No synthetic data generator is used.
    # -----------------------------------------------------------------------
    import time

    from ..core import PDEConfig

    micro = _build_selftest_micro()
    print("patient summary:", micro.summary())

    # Keep the transport solve to a few seconds: coarser dt than the fine default,
    # over a 6h window. (The plasma-averaging simplification makes the coarse dt
    # harmless for the Stage 0 readouts.)
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=200)

    # Snapshot caller state to prove non-mutation.
    state_before = micro.cells.state.copy()
    oxygen_before = micro.fields.oxygen

    cfg = SimConfig(pde=pde, seed=0)
    t0 = time.time()
    result = simulate_response(micro, cfg)
    elapsed = time.time() - t0

    print(f"\n=== engine self-test ({micro.domain.shape}, {elapsed:.2f}s) ===")
    print(f"predicted_response  = {result.predicted_response:.4f}")
    print(f"kill_fraction       = {result.kill_fraction:.4f}")
    print(f"viable_tumor_frac   = {result.viable_tumor_fraction:.4f}")
    print(f"penetration_index   = {result.extras['penetration_index']:.4f}")
    print(f"plasma_mean/peak    = {result.extras['plasma_mean']:.3f} / "
          f"{result.extras['plasma_peak']:.3f} uM")
    print(f"mean_auc            = {result.extras['mean_auc']:.4g} uM*s")
    print(f"mean_oxygen         = {result.extras['mean_oxygen']:.4f}")
    print(f"hypoxic_tumor_frac  = {result.extras['hypoxic_tumor_fraction']:.4f}")
    print(f"viable_history      = {np.round(result.viable_history, 4).tolist()}")
    print(f"history times (s)   = {np.round(result.times, 1).tolist()}")

    # --- structural asserts --------------------------------------------------
    assert isinstance(result, SimulationResult)
    assert 0.0 <= result.predicted_response <= 1.0, result.predicted_response
    assert 0.0 <= result.viable_tumor_fraction <= 1.0
    assert abs(result.kill_fraction - (1.0 - result.viable_tumor_fraction)) < 1e-12
    assert abs(result.predicted_response - result.kill_fraction) < 1e-12
    for key in ("oxygen", "ifp", "auc", "vascular_density"):
        assert key in result.fields, f"missing field {key}"
        assert result.fields[key].shape == micro.domain.shape, key
    assert result.viable_history.shape == result.times.shape
    assert result.viable_history[0] == 1.0, "zero-exposure viable frac must be 1.0"
    # Monotone non-increasing viable fraction as exposure accumulates.
    assert np.all(np.diff(result.viable_history) <= 1e-9), result.viable_history
    assert result.transport is not None

    # --- non-mutation of the caller -----------------------------------------
    assert np.array_equal(micro.cells.state, state_before), "caller cells mutated!"
    assert micro.fields.oxygen is oxygen_before, "caller fields.oxygen mutated!"

    # --- determinism ---------------------------------------------------------
    result_b = simulate_response(micro, SimConfig(pde=pde, seed=0))
    assert result.predicted_response == result_b.predicted_response, "non-deterministic"

    # --- transport coupling: high BBB penetrance kills >= low penetrance -----
    # Two drugs identical except bbb_penetrance, so any kill difference is due to
    # the delivery -> transport -> exposure coupling alone.
    drug_lo = DrugParams(name="lowbbb", diffusivity=300.0, ic50=10.0, hill=1.3,
                         bbb_penetrance=0.05, clearance=3e-4)
    drug_hi = DrugParams(name="highbbb", diffusivity=300.0, ic50=10.0, hill=1.3,
                         bbb_penetrance=0.5, clearance=3e-4)
    kill_lo = simulate_response(
        micro, SimConfig(drug=drug_lo, pde=pde, seed=0)
    ).kill_fraction
    kill_hi = simulate_response(
        micro, SimConfig(drug=drug_hi, pde=pde, seed=0)
    ).kill_fraction
    print(f"\nBBB coupling: kill(bbb=0.05)={kill_lo:.4f}  "
          f"kill(bbb=0.50)={kill_hi:.4f}")
    assert kill_hi >= kill_lo, (
        f"higher BBB penetrance ({kill_hi:.4f}) should kill >= lower "
        f"({kill_lo:.4f}) — transport is not coupled into the kill"
    )

    # --- time-resolved plasma driver -----------------------------------------
    # (a) it delivers exactly the same TOTAL drug as the static path, so a
    #     difference between them is a schedule effect and not a dose effect;
    # (b) it actually changes the exposure field (the static path cannot);
    # (c) with a constant plasma level the two paths coincide.
    from ..core import PKConfig

    pk_multi = PKConfig(plasma_peak=50.0, half_life=3600.0 * 2, n_doses=4,
                        dose_interval=3600.0 * 1.5)
    pde_tr = PDEConfig(dt=pde.dt, total_time=pde.total_time, convection=True,
                       save_every=pde.save_every, time_resolved=True)
    res_static = simulate_response(micro, SimConfig(pk=pk_multi, pde=pde, seed=0))
    res_tr = simulate_response(micro, SimConfig(pk=pk_multi, pde=pde_tr, seed=0))
    print(f"\ntime-resolved PK: mean_auc static={res_static.extras['mean_auc']:.6g}  "
          f"time-resolved={res_tr.extras['mean_auc']:.6g}  "
          f"(delivered plasma mean {res_static.extras['plasma_mean']:.6f} vs "
          f"{res_tr.extras['plasma_delivered_mean']:.6f} uM)")
    assert res_static.extras["time_resolved"] == 0.0
    assert res_tr.extras["time_resolved"] == 1.0
    assert abs(res_tr.extras["plasma_delivered_mean"]
               - res_static.extras["plasma_mean"]) < 1e-9, "total dose not matched"
    assert res_tr.extras["mean_auc"] != res_static.extras["mean_auc"], (
        "time_resolved=True produced an identical exposure field — the plasma "
        "profile is not reaching the solver"
    )
    # Constant plasma level (single dose, half-life astronomically longer than the
    # window, so plasma(t) is flat to ~1e-14). The two paths are then analytically
    # identical and must agree to round-off.
    pk_flat = PKConfig(plasma_peak=50.0, half_life=1e18, n_doses=1,
                       dose_interval=pde.total_time)
    auc_flat_static = simulate_response(
        micro, SimConfig(pk=pk_flat, pde=pde, seed=0)).extras["mean_auc"]
    auc_flat_tr = simulate_response(
        micro, SimConfig(pk=pk_flat, pde=pde_tr, seed=0)).extras["mean_auc"]
    rel = abs(auc_flat_tr - auc_flat_static) / auc_flat_static
    print(f"constant-plasma limit: static={auc_flat_static:.12g} "
          f"time-resolved={auc_flat_tr:.12g}  rel diff={rel:.2e}")
    assert rel < 1e-12, f"paths do not coincide at constant plasma: {rel:.2e}"

    print("\nengine self-test PASSED")
