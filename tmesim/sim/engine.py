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
   reduced to a single **representative sustained plasma level** — its
   time-average — and turned into a static spatial delivery source
   (:func:`~tmesim.treatment.drugs.delivery_source`). See the "Stage 0
   plasma-averaging simplification" note below.
5. **Transport.** Convection-diffusion-reaction PDE
   (:func:`~tmesim.physics.diffusion.simulate_drug_transport`) integrates the
   drug field and returns the exposure ``auc`` (uM*s) at every grid point.
6. **Kill.** :func:`~tmesim.cells.pkpd.apply_drug_kill` maps each cell's local AUC
   through a Hill dose-response (modulated by sensitivity + hypoxia/quiescence
   resistance) and marks killed cells DEAD.
7. **Readouts.** Final viable tumor fraction, kill fraction, and the predicted
   response (see the response map below), plus a short viable-fraction trajectory.

Stage 0 plasma-averaging simplification
---------------------------------------
The transport solver applies a *static* spatial source every timestep (its plasma
time-profile is held at 1). A fully time-resolved infusion would re-scale the
source by ``plasma(t)`` each step. For Stage 0 we instead drive the solver with a
single representative **sustained** plasma level equal to the time-average of the
multi-dose plasma profile over the exposure window::

    plasma_mean = mean_t plasma_profile(pk, times)
    source      = vascular_density * bbb_penetrance * plasma_mean   (uM/s)

This preserves the *total* delivered drug (mean level x window duration) and the
correct relative ranking of drugs/patients — which is all the Stage 0 kill-test
scores — while keeping the solve to a single factorized PDE integration. A
time-resolved driver is a Stage-1 refinement.

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
from ..treatment.drugs import delivery_source, get_drug

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

    Returns
    -------
    SimulationResult
        ``predicted_response`` (= kill fraction, in [0, 1]),
        ``viable_tumor_fraction``, ``kill_fraction``, a short ``viable_history``
        trajectory with matching ``times``, the driving ``transport`` result,
        ``fields`` for plotting (oxygen / ifp / auc / vascular_density), and
        diagnostic ``extras``.
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

    # 4. Drug + delivery source (Stage 0 plasma-averaging simplification).
    drug = _resolve_drug(cfg)
    times = np.arange(0.0, cfg.pde.total_time, cfg.pde.dt)
    plasma = plasma_profile_safe(cfg, times)
    plasma_mean = float(plasma.mean()) if plasma.size else 0.0
    source = delivery_source(work, drug, plasma_mean)

    # 5. Transport PDE -> exposure (AUC) field.
    transport = simulate_drug_transport(
        work,
        drug,
        source,
        cfg.pde,
        velocity=(velocity if cfg.pde.convection else None),
    )

    # 6. Pharmacodynamic kill over the treatment course.
    apply_drug_kill(work, transport, cfg.pd, drug, rng=rng)

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
            "plasma_peak": float(plasma.max()) if plasma.size else 0.0,
            "mean_auc": float(transport.auc[tissue].mean()) if tissue.any() else 0.0,
            "mean_oxygen": mean_oxygen,
            "hypoxic_tumor_fraction": hypoxic_tumor_frac,
            "bbb_penetrance": float(drug.bbb_penetrance),
            "n_tumor": float(n_tumor),
        },
    )


def plasma_profile_safe(cfg: SimConfig, times: np.ndarray) -> np.ndarray:
    """Evaluate the systemic plasma profile, importing lazily to avoid a hard
    module-load dependency cycle. Thin wrapper over
    :func:`tmesim.treatment.drugs.plasma_profile`."""
    from ..treatment.drugs import plasma_profile

    return plasma_profile(cfg.pk, times)


if __name__ == "__main__":
    # -----------------------------------------------------------------------
    # Self-test: run one synthetic 120x120 patient end to end and assert the
    # integrator is well-formed, deterministic, and that stronger BBB penetrance
    # (better delivery) yields at least as much kill as weak penetrance — the
    # sanity check that transport is actually coupled into the kill.
    # -----------------------------------------------------------------------
    import time

    from ..core import PDEConfig
    from ..data.synthetic import generate_synthetic_patient

    rec = generate_synthetic_patient("ENGINE-SELFTEST", seed=0)
    micro = rec.micro
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

    print("\nengine self-test PASSED")
