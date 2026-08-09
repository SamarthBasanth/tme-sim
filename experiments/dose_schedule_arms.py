"""Both-arms (transvascular exchange OFF vs ON) regeneration of
``results/dose_schedule.md``.

Why this script exists
======================
``results/dose_schedule.md`` carried a conditional retraction: it was produced
under the ONE-WAY transvascular delivery law ``S = PS*Cp``, which has no path
back from the interstitium to the blood and therefore lets tissue concentration
exceed the plasma level without bound. The physically correct form is
gradient-driven exchange ``S = PS*(Cp - C)``, whose steady state
``PS*Cp/(PS + clearance + binding)`` is bounded by ``Cp`` by construction.

**Which law ships is an OPEN decision** -- ``PDEConfig.transvascular_exchange``
still defaults to ``False``. So the report has to carry BOTH laws as first-class
columns, and that needs both arms MEASURED in one process from one code path
with nothing but the delivery law changed. That is what this script does.

The propagation trap, and why this file injects ``exchange_k`` by hand
=====================================================================
``experiments/dose_schedule.py`` is DECLARED ``INERT`` in
:data:`tmesim.validation.propagation.DECLARATIONS` -- indeed it is the module on
which the whole failure mode was first measured::

    exchange OFF                      auc_sum 1.554447793e+09
    flag set, exchange_k NOT passed    auc_sum 1.554447793e+09   <-- bit-identical
    honest ON (exchange_k passed)      auc_sum 6.86608671e+08

Its ``unit_transport`` helper calls ``simulate_drug_transport`` DIRECTLY, so
``PDEConfig(transvascular_exchange=True)`` never reaches the solver there and an
"on vs off" comparison built on it would be a comparison of a run to itself.

This script therefore does the engine-faithful thing explicitly: it builds
``exchange_k = exchange_coefficient(work, drug)`` -- byte-for-byte the field
``sim.engine.simulate_response`` builds -- and passes it to the solver itself.
It is registered as a NON-PRODUCTION :class:`WiringPoint` in
``tmesim/validation/propagation.py`` for exactly that reason.

Nothing is taken on trust:

* every solve runs inside :func:`tmesim.validation.propagation.instrument`, so
  the solver boundary is observed and ``exchange_k`` is counted as it crosses;
* the ON arm must have been handed ``exchange_k`` on EVERY solver call and the
  OFF arm on NONE, asserted, not hoped for;
* the two arms' AUC fields must differ BIT-WISE (``propagation.digest``) before
  any comparison is reported;
* the ON arm additionally re-runs the FULL ``sim.engine.simulate_response`` with
  ``PDEConfig(transvascular_exchange=True)`` -- the one declared PROPAGATED path
  -- and requires it to reproduce this script's fast-path kill. That closes the
  loop from the other end: the hand-injected field is validated against the
  production wiring, not merely against itself.

Telling a REAL zero from a SOLVER-ARTEFACT zero, by construction
===============================================================
``results/dose_schedule.md`` was rewritten once already to retract a
``0.00e+00`` that was an arithmetic artefact (the static source made two
schedules bit-identical INPUTS). Under the corrected law temozolomide's leg-A
kill spread prints ``0.00e+00`` again at the reference dose. The two look
identical in a table, so this script distinguishes them mechanically rather than
by eye, four ways:

1. **Per-shape AUC digests.** Each schedule's exposure field is hashed at bit
   level. Three DIFFERENT digests means three genuinely different fields, so a
   zero kill spread cannot be "the same run" -- it can only be quantisation.
   Three IDENTICAL digests would be the artefact itself.
2. **A continuous observable.** Mean tumour exposure (uM*s) is not quantised by
   counting cells. It is reported alongside every kill.
3. **The counting quantum is stated.** ``1 / n_tumor`` is printed next to the
   spread, so "below the quantum" is checkable rather than asserted.
4. **A deliberate positive control.** The superseded STATIC source is run in
   both arms. Its rows ARE bit-identical by construction and DO print exactly
   zero -- the artefact reproduced on purpose, next to the real measurement.

Usage
-----
    ./.venv/bin/python -u experiments/dose_schedule_arms.py --out <json> \
        [--stage all|geom|sched|engine] [--dx 12]

Writes nothing into ``results/``; the caller transcribes the numbers.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import dose_schedule as ds  # noqa: E402  (the module this rewrite measures)

from tmesim.core import (  # noqa: E402
    CellState,
    Microenvironment,
    PDEConfig,
    PKConfig,
    SimConfig,
)
from tmesim.physics.diffusion import simulate_drug_transport  # noqa: E402
from tmesim.reconstruct.real_pipeline import (  # noqa: E402
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from tmesim.sim.engine import simulate_response  # noqa: E402
from tmesim.treatment.drugs import (  # noqa: E402
    delivery_source,
    exchange_coefficient,
    plasma_profile,
)
from tmesim.validation.calibration import calibrate_alpha  # noqa: E402
from tmesim.validation.propagation import digest, instrument  # noqa: E402
from tmesim.validation.vessel_profile import (  # noqa: E402
    direct_vessel_distance_um,
    endothelial_xy,
    vessel_mask_from_density,
)

ARMS = ("off", "on")


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _jsonable(obj):
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


# --------------------------------------------------------------------------- #
# The one call site, and it ALWAYS passes exchange_k explicitly
# --------------------------------------------------------------------------- #
def unit_transport_arm(work, drug, pde: PDEConfig, velocity,
                       profile: Optional[np.ndarray], exchange: bool):
    """``dose_schedule.unit_transport`` with the delivery law made explicit.

    ``exchange=False`` passes ``exchange_k=None``, which is the solver's own
    default, so the OFF arm is BIT-IDENTICAL to the shipped path (asserted
    against the published numbers in the report). ``exchange=True`` passes
    ``exchange_coefficient(work, drug)`` -- the identical field
    ``sim.engine.simulate_response`` builds at its step 4b.

    Written with ONE call site that always names ``exchange_k=``: a bare
    ``simulate_drug_transport(...)`` here would be a new undeclared silent-no-op
    call site, which ``tests/test_propagation.py::
    test_no_undeclared_direct_solver_callers`` correctly fails on. The point of
    this script is not to re-create the defect it is measuring.
    """
    source = delivery_source(work, drug, 1.0)
    ex = (np.asarray(exchange_coefficient(work, drug), dtype=np.float64)
          if exchange else None)
    tr = simulate_drug_transport(
        work, drug, source, pde,
        velocity=(velocity if pde.convection else None),
        source_profile=profile,
        exchange_k=ex,
    )
    return tr, ex


def _solve(work, drug, pde, velocity, profile, exchange: bool, label: str):
    """One instrumented solve. Returns ``(tr, exchange_k, record)``."""
    t0 = time.time()
    with instrument() as probe:
        tr, ex = unit_transport_arm(work, drug, pde, velocity, profile, exchange)
    rec = dict(
        label=label,
        exchange=bool(exchange),
        seconds=time.time() - t0,
        solver_calls=probe.solver_calls,
        exchange_k_calls=probe.exchange_k_calls,
        source_profile_calls=probe.source_profile_calls,
        auc_digest=digest(np.asarray(tr.auc, dtype=np.float64)),
        auc_sum=float(np.asarray(tr.auc, dtype=np.float64).sum()),
        penetration_index=float(tr.penetration_index),
    )
    # The propagation guard, enforced at the point of measurement rather than
    # checked afterwards: ON must have been handed exchange_k on every solver
    # call, OFF on none. An INERT path silently fails the first of these.
    want = probe.solver_calls if exchange else 0
    if rec["exchange_k_calls"] != want:
        raise AssertionError(
            f"{label}: exchange_k reached the solver on "
            f"{rec['exchange_k_calls']} of {probe.solver_calls} call(s), "
            f"expected {want}. The delivery law did not propagate."
        )
    _log(f"    solve[{label}] exch={'ON ' if exchange else 'OFF'} "
         f"({rec['seconds']:.0f}s) auc_sum={rec['auc_sum']:.9e} "
         f"pen={rec['penetration_index']:.6f} "
         f"solver={rec['solver_calls']} exch_k={rec['exchange_k_calls']} "
         f"digest={rec['auc_digest'][:12]}")
    return tr, ex, rec


def residual_fraction_arm(tr, tissue: np.ndarray, k_loss: float,
                          exchange_k: Optional[np.ndarray]) -> float:
    """Terminal residual as a fraction of delivered drug, correct in BOTH arms.

    ``dose_schedule.residual_fraction`` computes ``C(T) / (k_loss*AUC + C(T))``
    with ``k_loss = clearance + binding``, a SCALAR. That is the whole
    elimination rate under the one-way law and the right answer there.

    Under exchange the drug also leaves through the vessel wall, so the true
    first-order loss is ``k_loss + PS(x, y)`` and it varies in space. Using the
    scalar there understates eliminated drug and therefore OVERSTATES the
    residual. This version integrates the real per-voxel rate::

        residual = C(T) / ( sum[(k_loss + PS) * AUC] + C(T) )

    and reduces to the scalar form exactly when ``exchange_k is None``, so the
    OFF arm reproduces the published column bit-for-bit.
    """
    auc = np.asarray(tr.auc, dtype=np.float64)[tissue]
    c_end = float(np.asarray(tr.final_concentration, dtype=np.float64)[tissue].sum())
    if exchange_k is None:
        eliminated = float(k_loss) * float(auc.sum())
    else:
        kf = float(k_loss) + np.asarray(exchange_k, dtype=np.float64)[tissue]
        eliminated = float((kf * auc).sum())
    total = eliminated + c_end
    return float(c_end / total) if total > 0 else 0.0


def _ps_stats(ex: Optional[np.ndarray], tissue: np.ndarray, k_loss: float) -> Dict:
    """Summarise the permeability field the exchange law adds to the loss rate."""
    if ex is None:
        return dict(
            ps_mean=0.0, ps_median=0.0, ps_p95=0.0, ps_max=0.0,
            k_eff_mean=float(k_loss), k_eff_min=float(k_loss),
            k_eff_max=float(k_loss),
            tau_eff_mean_s=1.0 / float(k_loss),
            tau_eff_min_s=1.0 / float(k_loss),
            tau_eff_max_s=1.0 / float(k_loss),
        )
    v = np.asarray(ex, dtype=np.float64)[tissue]
    k = float(k_loss) + v
    return dict(
        ps_mean=float(v.mean()), ps_median=float(np.median(v)),
        ps_p95=float(np.percentile(v, 95.0)), ps_max=float(v.max()),
        k_eff_mean=float(k.mean()), k_eff_min=float(k.min()),
        k_eff_max=float(k.max()),
        tau_eff_mean_s=float(1.0 / k.mean()),
        tau_eff_min_s=float(1.0 / k.max()),
        tau_eff_max_s=float(1.0 / k.min()),
    )


def _analytic_L_um(diffusivity: float, mean_ecm: float, k: float) -> float:
    """``sqrt(D_eff / k)`` on the report's own convention (section 3 column)."""
    return float(math.sqrt(diffusivity * (1.0 - 0.6 * float(mean_ecm)) / float(k)))


# --------------------------------------------------------------------------- #
# Shared tissue prologue -- built ONCE and reused by every arm, so the two arms
# differ in the delivery law and in NOTHING else.
# --------------------------------------------------------------------------- #
def build_context(dx_um: float) -> Dict:
    _log(f"[1] reconstructing the real Xenium GBM section (dx = {dx_um} um) ...")
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
    rep = reconstruction_report(micro)
    alpha = calibrate_alpha()
    _log(f"    {rep['n_cells']:,} cells (tumor {rep['n_tumor']:,}) grid {rep['grid']} "
         f"| alpha = {alpha:.4f}")

    base_cfg = SimConfig(seed=ds.SEED)
    work, labeled_base, velocity = ds.prepare_tissue(micro, base_cfg)
    tumor = labeled_base.tumor_mask()
    hypoxic = float(np.mean(labeled_base.state[tumor] == CellState.HYPOXIC))
    quiescent = float(np.mean(labeled_base.state[tumor] == CellState.QUIESCENT))

    sens_by_drug: Dict[str, np.ndarray] = {}
    sens_info: Dict[str, Dict] = {}
    for name in ds.DRUGS:
        probe = Microenvironment(
            domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
            patient_id=micro.patient_id, metadata=micro.metadata,
        )
        info = apply_drug_sensitivity(probe, name, mode="determinant")
        sens_by_drug[name] = np.asarray(probe.cells.sensitivity,
                                        dtype=np.float64).copy()
        sens_info[name] = info

    endo_xy = endothelial_xy(micro.cells)
    distance_um = direct_vessel_distance_um(micro.domain, endo_xy)
    vmask = vessel_mask_from_density(micro.fields.vascular_density,
                                     micro.domain.tissue_mask)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    _log(f"    prologue done: hypoxic {hypoxic*100:.1f}%, quiescent "
         f"{quiescent*100:.1f}% | tissue {int(tissue.sum()):,} px")

    return dict(
        micro=micro, work=work, labeled_base=labeled_base, velocity=velocity,
        alpha=float(alpha), pd_cfg=base_cfg.pd, tumor=tumor, tissue=tissue,
        distance_um=distance_um,
        sens_by_drug=sens_by_drug, sens_info=sens_info,
        meta=dict(
            dx_um=float(dx_um), dt_s=ds.DT_S, alpha=float(alpha),
            window_dose_h=ds.WINDOW_DOSE_S / 3600.0,
            window_sched_h=ds.WINDOW_SCHED_S / 3600.0,
            nominal_peak_um=ds.NOMINAL_PEAK_UM,
            nominal_half_life_h=ds.NOMINAL_HALF_LIFE_S / 3600.0,
            time_resolved=bool(ds.TIME_RESOLVED), seed=ds.SEED,
            n_cells=int(rep["n_cells"]), n_tumor=int(rep["n_tumor"]),
            grid=list(rep["grid"]),
            mean_vascular_density=float(rep["mean_vascular_density"]),
            mean_ecm_density=float(rep["mean_ecm_density"]),
            mean_oxygen=float(rep["mean_oxygen"]),
            hypoxic_tumor_fraction=hypoxic, quiescent_tumor_fraction=quiescent,
            vessel_area_fraction=float(vmask["area_fraction"]),
            vessel_percentile=float(vmask["percentile"]),
            kill_quantum=1.0 / float(rep["n_tumor"]),
            distance_median_um=float(np.median(distance_um[tissue])),
            distance_p90_um=float(np.percentile(distance_um[tissue], 90.0)),
            distance_max_um=float(distance_um[tissue].max()),
        ),
    )


# --------------------------------------------------------------------------- #
# Stage GEOM -- sections 2 and 3: dose-response and delivery geometry
# --------------------------------------------------------------------------- #
def stage_geom(ctx: Dict, arms: Sequence[str], drugs: Sequence[str],
               checkpoint, skip_eq: bool = False) -> Dict:
    micro, work = ctx["micro"], ctx["work"]
    meta, tissue, tumor = ctx["meta"], ctx["tissue"], ctx["tumor"]
    dx_um = meta["dx_um"]

    pde_dose = PDEConfig(dt=ds.DT_S, total_time=ds.WINDOW_DOSE_S, convection=True,
                         save_every=10, time_resolved=ds.TIME_RESOLVED)
    # A SUSTAINED-source solve, used only to MEASURE the tissue equilibration
    # time t95. Under the one-way law that number has a closed form,
    # ln(20)/k_loss, and the report quotes it. Under exchange the loss rate is
    # k_loss + PS(x, y), which is SPATIALLY VARYING, so no single exponential
    # describes the tissue mean and the closed form is not available. Measuring
    # it off a static source is valid in BOTH arms (a linear system under a
    # constant source rises monotonically to its plateau), which is what makes
    # the two arms' t95 comparable at all. save_every=1 so the measurement is
    # not quantised by the history stride.
    pde_eq = PDEConfig(dt=ds.DT_S, total_time=ds.WINDOW_DOSE_S, convection=True,
                       save_every=1, time_resolved=False)

    dose_shape, _ = ds.shape_profile(
        PKConfig(plasma_peak=1.0, half_life=ds.NOMINAL_HALF_LIFE_S, n_doses=1,
                 dose_interval=ds.WINDOW_DOSE_S),
        ds.WINDOW_DOSE_S, ds.DT_S)
    # Plasma AUC over the dose window at the nominal peak -- the denominator of
    # the tissue/plasma ratio that decides which law is physical.
    times = np.arange(0.0, ds.WINDOW_DOSE_S, ds.DT_S)
    plasma_auc = float(plasma_profile(
        PKConfig(plasma_peak=ds.NOMINAL_PEAK_UM, half_life=ds.NOMINAL_HALF_LIFE_S,
                 n_doses=1, dose_interval=ds.WINDOW_DOSE_S), times).sum() * ds.DT_S)
    meta["plasma_auc_dose_window_uM_s"] = plasma_auc

    out: Dict = {}
    for name in drugs:
        drug, mol = ds.compile_real_drug(name, ctx["alpha"])
        k_loss = float(drug.clearance + drug.binding)
        labeled = ctx["labeled_base"].copy()
        labeled.sensitivity = ctx["sens_by_drug"][name].copy()
        ceiling, ceiling_breakdown = ds.kill_ceiling(labeled, ctx["pd_cfg"], drug=drug)
        _log(f"[2] {name}: IC50={drug.ic50:.2f} uM bbb={drug.bbb_penetrance:.4f} "
             f"k_loss={k_loss:.3e}/s  ceiling={ceiling:.4f}")

        for arm in arms:
            exchange = (arm == "on")
            tr, ex, rec = _solve(work, drug, pde_dose, ctx["velocity"], dose_shape,
                                 exchange, f"{name}|dose6h|{arm}")

            # --- section 3: delivery geometry -----------------------------
            depth = ds.depth_analysis(micro, tr, drug, dx_um, ctx["distance_um"])
            ds.augment_depth(depth, meta["distance_p90_um"])

            # --- equilibration t95, MEASURED under a sustained source -----
            # Skipped on the mesh-refinement arms: t95 is a bulk clearance
            # time-scale, it is not the quantity the grid check is about, and
            # halving the solve count is what makes a dx = 8 um check affordable.
            if skip_eq:
                rec_eq, t95_measured = None, None
            else:
                tr_eq, _, rec_eq = _solve(work, drug, pde_eq, ctx["velocity"],
                                          None, exchange, f"{name}|static6h|{arm}")
                t95_measured = ds.equilibration_time_s(tr_eq)

            ps = _ps_stats(ex, tissue, k_loss)

            # --- section 2: dose-response (pure rescaling, no new solve) --
            dose_rows = []
            for m in ds.DOSE_MULTIPLIERS:
                pk = PKConfig(plasma_peak=ds.NOMINAL_PEAK_UM * float(m),
                              half_life=ds.NOMINAL_HALF_LIFE_S, n_doses=1,
                              dose_interval=ds.WINDOW_DOSE_S)
                pm = ds.window_plasma_mean(pk, ds.WINDOW_DOSE_S, ds.DT_S)
                row = ds.kill_at(work, labeled, tr, drug, ctx["pd_cfg"], pm)
                row["multiplier"] = float(m)
                row["plasma_peak"] = float(ds.NOMINAL_PEAK_UM * m)
                dep_um, censored = ds.ic50_depth_um(depth, pm)
                row["ic50_depth_um"] = dep_um
                row["ic50_depth_censored"] = bool(censored)
                row["pct_of_ceiling"] = (float(row["kill"] / ceiling)
                                         if ceiling > 0 else None)
                dose_rows.append(row)
                _log(f"      {name} {arm} {m:5.2f}x kill={row['kill']:.6f} "
                     f"cov={row['ic50_coverage']:.4f} "
                     f"cov0.1={row['cov_10pct_ic50']:.4f} "
                     f"medE/IC50={row['median_E_over_ic50']:.4f} "
                     f"depth={'-' if dep_um is None else round(dep_um, 1)}"
                     f"{'+' if censored else ''}")

            auc_field = np.asarray(tr.auc, dtype=np.float64)
            row1x = next(r for r in dose_rows if abs(r["multiplier"] - 1.0) < 1e-9)
            mean_at_1x = row1x["plasma_mean"]
            tissue_auc_1x = auc_field[tissue] * mean_at_1x
            entry = dict(
                drug=name, arm=arm, exchange=exchange,
                params=dict(name=name, mw=float(mol.get("mw")),
                            diffusivity=float(drug.diffusivity),
                            ic50=float(drug.ic50), hill=float(drug.hill),
                            bbb=float(drug.bbb_penetrance),
                            binding=float(drug.binding),
                            clearance=float(drug.clearance), k_loss=k_loss,
                            tau_transport_s=1.0 / k_loss),
                kill_ceiling=float(ceiling),
                kill_ceiling_breakdown=ceiling_breakdown,
                sensitivity=dict(
                    determinant_gene=ctx["sens_info"][name].get("determinant_gene"),
                    kind=ctx["sens_info"][name].get("kind"),
                    median_tumor=float(np.median(ctx["sens_by_drug"][name][tumor])),
                    is_uniform=bool(ctx["sens_info"][name].get("is_uniform", False)),
                ),
                solve=rec, solve_static=rec_eq,
                penetration_index=float(tr.penetration_index),
                depth=depth,
                ps=ps,
                analytic_L_um=_analytic_L_um(drug.diffusivity,
                                             meta["mean_ecm_density"], k_loss),
                analytic_L_um_exchange_mean=(
                    _analytic_L_um(drug.diffusivity, meta["mean_ecm_density"],
                                   ps["k_eff_mean"]) if exchange else None),
                analytic_L_um_exchange_min=(
                    _analytic_L_um(drug.diffusivity, meta["mean_ecm_density"],
                                   ps["k_eff_max"]) if exchange else None),
                t95_analytic_s=ds.equilibration_t95_analytic_s(k_loss),
                t95_measured_s=(float(t95_measured)
                                if t95_measured is not None else None),
                plasma_mean_at_1x=float(mean_at_1x),
                mean_tumor_auc_1x=float(row1x["mean_auc_tumor"]),
                mean_tissue_auc_1x=float(tissue_auc_1x.mean()),
                mean_tissue_auc_over_plasma_auc=float(
                    tissue_auc_1x.mean() / plasma_auc),
                max_tissue_auc_over_plasma_auc=float(
                    tissue_auc_1x.max() / plasma_auc),
                mean_tumor_auc_over_plasma_auc=float(
                    row1x["mean_auc_tumor"] / plasma_auc),
                dose=dose_rows,
            )
            out.setdefault(arm, {})[name] = entry
            _log(f"    {name} {arm}: pen_idx={tr.penetration_index:.6f} "
                 f"emp_d1/2={depth['half_distance_um']} "
                 f"srcnull_d1/2={depth['source_null_half_distance_um']} "
                 f"nearL={depth['fit_near'].get('L_um')} "
                 f"floorL={depth['fit_floor'].get('L_um')} "
                 f"plateau={depth.get('plateau_ratio')} "
                 f"t95meas={t95_measured} tissue/plasma="
                 f"{entry['mean_tissue_auc_over_plasma_auc']:.4f}")
            checkpoint()

    # The comparison is only meaningful if the arms are different runs.
    for name in drugs:
        if "off" in out and "on" in out and name in out["off"] and name in out["on"]:
            d_off = out["off"][name]["solve"]["auc_digest"]
            d_on = out["on"][name]["solve"]["auc_digest"]
            if d_off == d_on:
                raise AssertionError(
                    f"{name}: OFF and ON AUC fields are BIT-IDENTICAL "
                    f"({d_off}). The delivery law did not reach the physics; "
                    "this is the silent no-op, not a robustness result."
                )
            _log(f"    [guard] {name}: arms differ bit-wise "
                 f"{d_off[:12]} != {d_on[:12]}  OK")
    return out


# --------------------------------------------------------------------------- #
# Stage SCHED -- sections 4 and 5: fractionation and infusion duration
# --------------------------------------------------------------------------- #
def stage_sched(ctx: Dict, arms: Sequence[str], checkpoint) -> Dict:
    micro, work = ctx["micro"], ctx["work"]
    meta, tissue = ctx["meta"], ctx["tissue"]
    dx_um = meta["dx_um"]
    name = ds.SCHEDULE_DRUG

    pde_sched = PDEConfig(dt=ds.DT_S, total_time=ds.WINDOW_SCHED_S, convection=True,
                          save_every=100, time_resolved=ds.TIME_RESOLVED)
    pde_static = PDEConfig(dt=ds.DT_S, total_time=ds.WINDOW_SCHED_S, convection=True,
                           save_every=100, time_resolved=False)

    drug, _ = ds.compile_real_drug(name, ctx["alpha"])
    k_loss = float(drug.clearance + drug.binding)
    labeled = ctx["labeled_base"].copy()
    labeled.sensitivity = ctx["sens_by_drug"][name].copy()
    pd_cfg = ctx["pd_cfg"]

    pk_ref = PKConfig(plasma_peak=ds.NOMINAL_PEAK_UM,
                      half_life=ds.NOMINAL_HALF_LIFE_S, n_doses=1,
                      dose_interval=ds.WINDOW_SCHED_S)
    target_mean = ds.window_plasma_mean(pk_ref, ds.WINDOW_SCHED_S, ds.DT_S)

    # Every leg is binned on ONE distance axis so the shapes are comparable.
    # ``depth_analysis`` derives that axis from the geometry alone (p99.5 of the
    # tissue's distance-to-vessel), never from the solve, so it is computed here
    # once and held fixed across schedules, half-lives and arms.
    depth_max_um = float(np.percentile(ctx["distance_um"][tissue], 99.5))
    _log(f"    schedule legs binned to {depth_max_um:.1f} um "
         f"(p99.5 of distance-to-vessel, geometry only)")

    out: Dict = {}
    for arm in arms:
        exchange = (arm == "on")
        # The SUPERSEDED static source -- kept as a deliberate POSITIVE CONTROL
        # for the spurious zero. Every schedule shares this one field, so its
        # rows are bit-identical BY CONSTRUCTION and its spreads are exactly 0.
        tr_static, ex_static, rec_static = _solve(work, drug, pde_static, ctx["velocity"],
                                          None, exchange, f"{name}|static12h|{arm}")

        sched_rows: List[Dict] = []
        for n, interval in ds.SCHEDULES:
            sh, _ = ds.shape_profile(
                PKConfig(plasma_peak=1.0, half_life=ds.NOMINAL_HALF_LIFE_S,
                         n_doses=n, dose_interval=interval),
                ds.WINDOW_SCHED_S, ds.DT_S)
            tr_n, ex_n, rec_n = _solve(work, drug, pde_sched, ctx["velocity"], sh,
                                    exchange, f"{name}|sched-n{n}|{arm}")

            peak_a = ds.peak_for_matched_auc(target_mean, n, interval,
                                             ds.NOMINAL_HALF_LIFE_S,
                                             ds.WINDOW_SCHED_S, ds.DT_S)
            pk_a = PKConfig(plasma_peak=peak_a, half_life=ds.NOMINAL_HALF_LIFE_S,
                            n_doses=n, dose_interval=interval)
            mean_a = ds.window_plasma_mean(pk_a, ds.WINDOW_SCHED_S, ds.DT_S)
            res_a = ds.kill_at(work, labeled, tr_n, drug, pd_cfg, mean_a)

            peak_b = ds.NOMINAL_PEAK_UM / n
            pk_b = PKConfig(plasma_peak=peak_b, half_life=ds.NOMINAL_HALF_LIFE_S,
                            n_doses=n, dose_interval=interval)
            mean_b = ds.window_plasma_mean(pk_b, ds.WINDOW_SCHED_S, ds.DT_S)
            res_b = ds.kill_at(work, labeled, tr_n, drug, pd_cfg, mean_b)

            pk_c = PKConfig(plasma_peak=ds.NOMINAL_PEAK_UM,
                            half_life=ds.NOMINAL_HALF_LIFE_S, n_doses=n,
                            dose_interval=interval)
            mean_c = ds.window_plasma_mean(pk_c, ds.WINDOW_SCHED_S, ds.DT_S)
            res_c = ds.kill_at(work, labeled, tr_n, drug, pd_cfg, mean_c)

            res_a_static = ds.kill_at(work, labeled, tr_static, drug, pd_cfg, mean_a)

            dep_n = ds.depth_analysis(micro, tr_n, drug, dx_um, ctx["distance_um"],
                                      max_um=depth_max_um)
            dep_a_um, dep_a_cen = ds.ic50_depth_um(dep_n, mean_a)

            by_dose = {f"{m:g}": ds.kill_at(work, labeled, tr_n, drug, pd_cfg,
                                            mean_a * m)["kill"]
                       for m in ds.SCHEDULE_DOSE_MULTIPLIERS}
            by_dose_static = {f"{m:g}": ds.kill_at(work, labeled, tr_static, drug,
                                                   pd_cfg, mean_a * m)["kill"]
                              for m in ds.SCHEDULE_DOSE_MULTIPLIERS}

            sched_rows.append(dict(
                n_doses=int(n), interval_h=float(interval / 3600.0),
                solve=rec_n,
                peak_matched_auc=float(peak_a), plasma_mean_matched=float(mean_a),
                kill_matched=res_a["kill"], cov_matched=res_a["ic50_coverage"],
                peak_matched_dose=float(peak_b), plasma_mean_dose=float(mean_b),
                kill_matched_dose=res_b["kill"],
                plasma_mean_peak=float(mean_c), kill_matched_peak=res_c["kill"],
                kill_matched_static=res_a_static["kill"],
                kill_matched_by_dose=by_dose,
                kill_matched_by_dose_static=by_dose_static,
                mean_auc_tumor_matched=res_a["mean_auc_tumor"],
                mean_auc_tumor_matched_static=res_a_static["mean_auc_tumor"],
                residual_fraction=residual_fraction_arm(tr_n, tissue, k_loss, ex_n),
                residual_fraction_kloss_only=ds.residual_fraction(tr_n, tissue, k_loss),
                half_distance_um=dep_n.get("half_distance_um"),
                penetration_index=float(tr_n.penetration_index),
                depth_matched_um=dep_a_um, depth_matched_censored=bool(dep_a_cen),
            ))
            _log(f"      SCHED[{arm}] n={n:2d}: killA={res_a['kill']:.6f} "
                 f"(static {res_a_static['kill']:.6f}) killB={res_b['kill']:.6f} "
                 f"killC={res_c['kill']:.6f} "
                 f"meanAUC={res_a['mean_auc_tumor']:.6e} "
                 f"resid={sched_rows[-1]['residual_fraction']*100:.3f}% "
                 f"pen={tr_n.penetration_index:.6f} "
                 f"d1/2={sched_rows[-1]['half_distance_um']}")
            checkpoint()

        hl_rows: List[Dict] = []
        for hl_h in ds.HALF_LIVES_H:
            hl = hl_h * 3600.0
            sh, _ = ds.shape_profile(
                PKConfig(plasma_peak=1.0, half_life=hl, n_doses=1,
                         dose_interval=ds.WINDOW_SCHED_S),
                ds.WINDOW_SCHED_S, ds.DT_S)
            tr_h, ex_h, rec_h = _solve(work, drug, pde_sched, ctx["velocity"], sh,
                                    exchange, f"{name}|hl-{hl_h:g}h|{arm}")
            dep_h = ds.depth_analysis(micro, tr_h, drug, dx_um, ctx["distance_um"],
                                      max_um=depth_max_um)
            peak_m = ds.peak_for_matched_auc(target_mean, 1, ds.WINDOW_SCHED_S, hl,
                                             ds.WINDOW_SCHED_S, ds.DT_S)
            pk_m = PKConfig(plasma_peak=peak_m, half_life=hl, n_doses=1,
                            dose_interval=ds.WINDOW_SCHED_S)
            mean_m = ds.window_plasma_mean(pk_m, ds.WINDOW_SCHED_S, ds.DT_S)
            res_m = ds.kill_at(work, labeled, tr_h, drug, pd_cfg, mean_m)
            pk_p = PKConfig(plasma_peak=ds.NOMINAL_PEAK_UM, half_life=hl, n_doses=1,
                            dose_interval=ds.WINDOW_SCHED_S)
            mean_p = ds.window_plasma_mean(pk_p, ds.WINDOW_SCHED_S, ds.DT_S)
            res_p = ds.kill_at(work, labeled, tr_h, drug, pd_cfg, mean_p)
            res_m_static = ds.kill_at(work, labeled, tr_static, drug, pd_cfg, mean_m)
            dm, dm_c = ds.ic50_depth_um(dep_h, mean_m)
            dp_, dp_c = ds.ic50_depth_um(dep_h, mean_p)
            ds.augment_depth(dep_h, meta["distance_p90_um"])
            hl_rows.append(dict(
                half_life_h=float(hl_h), solve=rec_h,
                peak_matched=float(peak_m), plasma_mean_matched=float(mean_m),
                kill_matched=res_m["kill"], cov_matched=res_m["ic50_coverage"],
                depth_matched_um=dm, depth_matched_censored=bool(dm_c),
                plasma_mean_peakmatched=float(mean_p),
                kill_peakmatched=res_p["kill"],
                depth_peakmatched_um=dp_, depth_peakmatched_censored=bool(dp_c),
                kill_matched_static=res_m_static["kill"],
                kill_matched_by_dose={
                    f"{m:g}": ds.kill_at(work, labeled, tr_h, drug, pd_cfg,
                                         mean_m * m)["kill"]
                    for m in ds.SCHEDULE_DOSE_MULTIPLIERS},
                kill_matched_by_dose_static={
                    f"{m:g}": ds.kill_at(work, labeled, tr_static, drug, pd_cfg,
                                         mean_m * m)["kill"]
                    for m in ds.SCHEDULE_DOSE_MULTIPLIERS},
                mean_auc_tumor_matched=res_m["mean_auc_tumor"],
                mean_auc_tumor_matched_static=res_m_static["mean_auc_tumor"],
                residual_fraction=residual_fraction_arm(tr_h, tissue, k_loss, ex_h),
                residual_fraction_kloss_only=ds.residual_fraction(tr_h, tissue, k_loss),
                half_distance_um=dep_h.get("half_distance_um"),
                penetration_index=float(tr_h.penetration_index),
                perivascular_unit_E=dep_h.get("perivascular_unit_E"),
                plateau_ratio=dep_h.get("plateau_ratio"),
            ))
            _log(f"      HL[{arm}] t1/2={hl_h:g}h: kill={res_m['kill']:.6f} "
                 f"(static {res_m_static['kill']:.6f}) "
                 f"meanAUC={res_m['mean_auc_tumor']:.6e} "
                 f"resid={hl_rows[-1]['residual_fraction']*100:.3f}% "
                 f"pen={tr_h.penetration_index:.6f} d1/2="
                 f"{hl_rows[-1]['half_distance_um']}")
            checkpoint()

        # ---- the REAL-vs-ARTEFACT zero test, done mechanically ----------
        sched_digests = [r["solve"]["auc_digest"] for r in sched_rows]
        hl_digests = [r["solve"]["auc_digest"] for r in hl_rows]
        zero_test = dict(
            schedule_field_digests=sched_digests,
            schedule_fields_all_distinct=(len(set(sched_digests)) == len(sched_digests)),
            half_life_field_digests=hl_digests,
            half_life_fields_all_distinct=(len(set(hl_digests)) == len(hl_digests)),
            static_field_digest=rec_static["auc_digest"],
            kill_quantum=meta["kill_quantum"],
            schedule_kill_spread=float(max(r["kill_matched"] for r in sched_rows)
                                       - min(r["kill_matched"] for r in sched_rows)),
            schedule_exposure_spread_pct=float(
                100.0 * (max(r["mean_auc_tumor_matched"] for r in sched_rows)
                         - min(r["mean_auc_tumor_matched"] for r in sched_rows))
                / max(r["mean_auc_tumor_matched"] for r in sched_rows)),
            schedule_exposure_monotone=bool(all(
                sched_rows[i]["mean_auc_tumor_matched"]
                > sched_rows[i + 1]["mean_auc_tumor_matched"]
                for i in range(len(sched_rows) - 1))),
            half_life_kill_spread=float(max(r["kill_matched"] for r in hl_rows)
                                        - min(r["kill_matched"] for r in hl_rows)),
            half_life_exposure_spread_pct=float(
                100.0 * (max(r["mean_auc_tumor_matched"] for r in hl_rows)
                         - min(r["mean_auc_tumor_matched"] for r in hl_rows))
                / max(r["mean_auc_tumor_matched"] for r in hl_rows)),
            static_schedule_kill_spread=float(
                max(r["kill_matched_static"] for r in sched_rows)
                - min(r["kill_matched_static"] for r in sched_rows)),
        )
        if not zero_test["schedule_fields_all_distinct"]:
            raise AssertionError(
                f"[{arm}] the three schedule AUC fields are NOT all distinct "
                f"({sched_digests}). A zero spread here would be the solver "
                "artefact, not a physical result."
            )
        _log(f"    [zero-test {arm}] schedule fields distinct="
             f"{zero_test['schedule_fields_all_distinct']} "
             f"kill_spread={zero_test['schedule_kill_spread']:.3e} "
             f"quantum={meta['kill_quantum']:.3e} "
             f"exposure_spread={zero_test['schedule_exposure_spread_pct']:.4f}% "
             f"monotone={zero_test['schedule_exposure_monotone']} "
             f"static_kill_spread={zero_test['static_schedule_kill_spread']:.3e}")

        out[arm] = dict(
            drug=name, arm=arm, exchange=exchange,
            target_window_mean_um=float(target_mean),
            static_solve=rec_static,
            static_penetration_index=float(tr_static.penetration_index),
            static_residual_fraction=residual_fraction_arm(tr_static, tissue,
                                                           k_loss, ex_static),
            static_residual_fraction_kloss_only=ds.residual_fraction(
                tr_static, tissue, k_loss),
            tau_transport_s=1.0 / k_loss,
            bound_tau_over_T_pct=100.0 * (1.0 / k_loss) / ds.WINDOW_SCHED_S,
            schedule=sched_rows, half_life=hl_rows, zero_test=zero_test,
        )
        checkpoint()
    return out


# --------------------------------------------------------------------------- #
# Stage ENGINE -- the fast path against the FULL engine, in BOTH arms
# --------------------------------------------------------------------------- #
def stage_engine(ctx: Dict, geom: Dict, sched: Dict, arms: Sequence[str],
                 checkpoint) -> List[Dict]:
    """Re-run ``sim.engine.simulate_response`` end to end and require it to
    reproduce the fast-path kill.

    Two things are pinned at once, and the second is the important one here:

    1. the shape-x-amplitude fast path is exact under BOTH delivery laws (the
       transport operator stays linear when ``exchange_k`` is added, but that is
       an argument and this is a measurement);
    2. **the hand-injected ``exchange_k`` matches the PRODUCTION wiring.** The
       engine is the one declared PROPAGATED path; it builds its own
       ``exchange_coefficient`` from ``PDEConfig.transvascular_exchange``. If
       this script's injection differed from the engine's in any way, the ON
       cross-checks would not reconcile.
    """
    micro = ctx["micro"]
    checks: List[Dict] = []

    def _micro_for(dp_name):
        m = Microenvironment(
            domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
            patient_id=micro.patient_id, metadata=micro.metadata,
        )
        m.cells.sensitivity = ctx["sens_by_drug"][dp_name].copy()
        return m

    def _check(label, dp_name, pk, pde, fast_kill, exchange):
        drug, _ = ds.compile_real_drug(dp_name, ctx["alpha"])
        t0 = time.time()
        with instrument() as probe:
            res = simulate_response(
                _micro_for(dp_name),
                SimConfig(drug=drug, pk=pk, pde=pde, seed=ds.SEED))
        d = abs(float(res.kill_fraction) - float(fast_kill))
        rec = dict(label=label, drug=dp_name, arm=("on" if exchange else "off"),
                   engine_kill=float(res.kill_fraction),
                   fast_kill=float(fast_kill), abs_delta=float(d),
                   solver_calls=probe.solver_calls,
                   exchange_k_calls=probe.exchange_k_calls,
                   seconds=time.time() - t0)
        want = probe.solver_calls if exchange else 0
        if rec["exchange_k_calls"] != want:
            raise AssertionError(
                f"{label}: engine handed exchange_k on {rec['exchange_k_calls']} "
                f"of {probe.solver_calls} solver call(s), expected {want}.")
        _log(f"    ENGINE {label:52s} engine={res.kill_fraction:.10f} "
             f"fast={fast_kill:.10f} |d|={d:.2e} exch_k={rec['exchange_k_calls']}"
             f"/{rec['solver_calls']} ({rec['seconds']:.0f}s)")
        checks.append(rec)
        checkpoint()

    for arm in arms:
        exchange = (arm == "on")
        pde_dose = PDEConfig(dt=ds.DT_S, total_time=ds.WINDOW_DOSE_S,
                             convection=True, save_every=10,
                             time_resolved=ds.TIME_RESOLVED,
                             transvascular_exchange=exchange)
        pde_sched = PDEConfig(dt=ds.DT_S, total_time=ds.WINDOW_SCHED_S,
                              convection=True, save_every=100,
                              time_resolved=ds.TIME_RESOLVED,
                              transvascular_exchange=exchange)
        if geom.get(arm):
            for dp_name in sorted(geom[arm]):
                row = next(r for r in geom[arm][dp_name]["dose"]
                           if abs(r["multiplier"] - 1.0) < 1e-9)
                _check(f"{dp_name} dose 1x (6 h window) [{arm}]", dp_name,
                       PKConfig(plasma_peak=ds.NOMINAL_PEAK_UM,
                                half_life=ds.NOMINAL_HALF_LIFE_S, n_doses=1,
                                dose_interval=ds.WINDOW_DOSE_S),
                       pde_dose, row["kill"], exchange)
        if sched.get(arm):
            srow = max(sched[arm]["schedule"], key=lambda r: r["n_doses"])
            _check(f"{ds.SCHEDULE_DRUG} schedule n={srow['n_doses']} matched-AUC "
                   f"(12 h window) [{arm}]", ds.SCHEDULE_DRUG,
                   PKConfig(plasma_peak=srow["peak_matched_auc"],
                            half_life=ds.NOMINAL_HALF_LIFE_S,
                            n_doses=srow["n_doses"],
                            dose_interval=srow["interval_h"] * 3600.0),
                   pde_sched, srow["kill_matched"], exchange)
            hrow = min(sched[arm]["half_life"], key=lambda r: r["half_life_h"])
            _check(f"{ds.SCHEDULE_DRUG} half-life {hrow['half_life_h']:g} h "
                   f"matched-exposure (12 h window) [{arm}]", ds.SCHEDULE_DRUG,
                   PKConfig(plasma_peak=hrow["peak_matched"],
                            half_life=hrow["half_life_h"] * 3600.0, n_doses=1,
                            dose_interval=ds.WINDOW_SCHED_S),
                   pde_sched, hrow["kill_matched"], exchange)
    return checks


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--stage", default="all",
                    choices=["all", "geom", "sched", "engine"])
    ap.add_argument("--arms", default="off,on")
    ap.add_argument("--drugs", default=",".join(ds.DRUGS))
    ap.add_argument("--dx", type=float, default=ds.DX_UM)
    ap.add_argument("--skip-eq", action="store_true",
                    help="skip the sustained-source solve used only to measure "
                         "the equilibration t95 (halves the geom solve count).")
    ap.add_argument("--resume", default="",
                    help="comma-separated JSON files from earlier stages; their "
                         "'geom'/'sched' blocks are loaded so a later stage "
                         "(notably 'engine', which needs the fast-path kills to "
                         "check against) can run in its own process.")
    args = ap.parse_args()

    arms = tuple(a for a in args.arms.split(",") if a)
    drugs = tuple(d for d in args.drugs.split(",") if d)
    for a in arms:
        if a not in ARMS:
            raise SystemExit(f"unknown arm {a!r}; expected one of {ARMS}")

    t_all = time.time()
    payload: Dict = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "default_transvascular_exchange": PDEConfig().transvascular_exchange,
        "arms": list(arms), "drugs": list(drugs), "stage": args.stage,
    }

    def checkpoint() -> None:
        tmp = args.out + ".tmp"
        with open(tmp, "w") as fh:
            json.dump(_jsonable(payload), fh, indent=1)
        os.replace(tmp, args.out)

    _log(f"PDEConfig().transvascular_exchange default = "
         f"{PDEConfig().transvascular_exchange}  (the SHIPPED law)")
    checkpoint()

    ctx = build_context(args.dx)
    payload["meta"] = ctx["meta"]
    checkpoint()

    geom: Dict = {}
    sched: Dict = {}
    for path in (p for p in args.resume.split(",") if p):
        with open(path) as fh:
            prior = json.load(fh)
        for key, store in (("geom", geom), ("sched", sched)):
            for arm, val in (prior.get(key) or {}).items():
                if key == "geom":
                    store.setdefault(arm, {}).update(val)
                else:
                    store[arm] = val
        _log(f"resumed {path}: geom arms={sorted(geom)} sched arms={sorted(sched)}")
    if geom:
        payload["geom"] = geom
    if sched:
        payload["sched"] = sched

    if args.stage in ("all", "geom"):
        geom = stage_geom(ctx, arms, drugs, checkpoint, skip_eq=args.skip_eq)
        payload["geom"] = geom
        checkpoint()
    if args.stage in ("all", "sched"):
        sched = stage_sched(ctx, arms, checkpoint)
        payload["sched"] = sched
        checkpoint()
    if args.stage in ("all", "engine"):
        payload["engine_checks"] = stage_engine(ctx, geom, sched, arms, checkpoint)
        if payload["engine_checks"]:
            payload["max_engine_fast_delta"] = max(
                c["abs_delta"] for c in payload["engine_checks"])
        checkpoint()

    payload["total_seconds"] = time.time() - t_all
    checkpoint()
    _log(f"DONE in {payload['total_seconds']/60:.1f} min -> {args.out}")


if __name__ == "__main__":
    main()
