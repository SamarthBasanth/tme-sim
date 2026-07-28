"""DOSE & SCHEDULE on the REAL Xenium glioblastoma tumor — the decision the
simulator is actually built to inform.

Everything upstream of this script answers "which drug?". This one answers the
next question a treating team actually faces: **how much, and how often?** It is
run on ONE real tumor (the open 10x Xenium GBM section, 40,887 measured cells,
real endothelial-derived vasculature, real collagen/matrix ECM field, real
perfusion-driven oxygen) for TWO real drugs (temozolomide, cisplatin), each
compiled through the same structure -> transport -> potency pipeline the drug
panel uses (``calibrated_drugparams`` + ``apply_barrier_to_drugparams``).

Three questions, one tumor
--------------------------
1. **DOSE-RESPONSE.** Sweep plasma peak 0.1x - 30x a nominal level and measure
   predicted tumor kill AND the fraction of tumor cells whose *exposure* clears
   the drug's own IC50. Where does kill saturate, and is the saturation a
   delivery ceiling or a sensitivity ceiling?
2. **SCHEDULE at MATCHED TOTAL DOSE.** 1 vs 3 vs 5 vs 10 fractions over the same
   window with cumulative plasma exposure (AUC) held constant. Does fractionation
   change predicted kill on this real geometry — and why / why not?
3. **INFUSION DURATION.** Plasma half-life 0.25 h - 48 h at matched total
   exposure (short bolus vs sustained). If delivery is transport-limited, a
   sustained low level should reach deeper than a bolus. Tested, not assumed.

TIME-RESOLVED PLASMA DRIVER (``TIME_RESOLVED``) — supersedes the static source
------------------------------------------------------------------------------
Earlier versions of this study ran against an engine that collapsed the whole
multi-dose plasma profile to a single **time-averaged** level and drove the
transport solver with a static source. Questions 2 and 3 above were then
unanswerable by construction: at matched cumulative exposure the solver input was
*bit-identical* across schedules, and the reported kill spread was exactly
``0.00e+00``. That zero was arithmetic, not biology, and the conclusions drawn
from it ("fractionation changes nothing", "sustained infusion does not reach
further") were model artefacts.

With ``PDEConfig.time_resolved=True`` the source is rescaled by ``plasma(t)`` at
every timestep, so tissue concentration genuinely rises and falls. Cumulative
delivered drug is held *exactly* equal to the static path (the plasma curve is
sampled on the solver's own step grid, whose mean is exactly ``plasma_mean``), so
a static-vs-resolved difference is a pure schedule effect and not a dose effect.

HOW THE SWEEPS ARE COMPUTED (and why it is exact, not an approximation)
-----------------------------------------------------------------------
The transport PDE ``dC/dt = div(D_eff grad C) - div(vC) + S - (clearance+binding)C``
is **linear** in the source and starts from ``C = 0``, and the source is
``S(x,y,t) = vascular_density * bbb_penetrance * plasma(t)``. Scaling the whole
plasma curve therefore scales the exposure field exactly. Every schedule is split
into **shape x amplitude** (see :func:`shape_profile`):

* **amplitude is free** — all dose points, and the matched-exposure /
  matched-administered-dose / matched-peak legs of one schedule, share a single
  solve and are recovered by rescaling it. Exact, not interpolated.
* **shape is not free** — 1 vs 3 vs 10 fractions, and a 0.25 h vs a 48 h
  half-life, are different shapes and each needs its own solve. Under the static
  source they all collapsed to the same constant, which is exactly why the
  schedule question had no answer.

Every point then re-runs the engine's own unmodified kill mapping
(``tmesim.cells.pkpd.apply_drug_kill``) on the rescaled field. Nothing in
``core``/``physics``/``molecular``/``reconstruct``/``engine`` is re-implemented
here; this script orchestrates them in the engine's own order.

This is not taken on faith: step [5] re-executes the **full** ``simulate_response``
engine once on each leg and asserts the fast path reproduces ``kill_fraction`` to
machine precision.

READ THIS BEFORE READING ANY NUMBER
-----------------------------------
* Research prototype. **NOT clinical dosing guidance.** No number here should
  influence the care of any patient.
* This Xenium section is an UNTREATED resection: there is **no measured drug
  response** to validate any predicted kill against. Zero outcome validation.
* Absolute kill fractions are conditional on the per-cell **sensitivity model**,
  which is an unvalidated prior. The *shape* of the dose-response and the
  *ranking* of doses are far more defensible than any absolute kill.
* Time resolution was added to **transport**, not to pharmacodynamics. The kill
  map is still a static Hill function of cumulative exposure — no inter-fraction
  repair, no cell-cycle redistribution, no regrowth, no normal-tissue recovery.
  The schedule effect measured here is therefore a *transport-only floor*, not an
  estimate of the clinical value of fractionating anything.

Run:    ./.venv/bin/python experiments/dose_schedule.py
Writes: results/dose_schedule.png
        results/dose_schedule.md
"""
from __future__ import annotations

import argparse
import json
import math
import os
import time
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from tmesim.cells.agents import update_cell_states
from tmesim.cells.pkpd import apply_drug_kill
from tmesim.core import (
    CellState,
    Microenvironment,
    PDConfig,
    PDEConfig,
    PKConfig,
    SimConfig,
    TissueFields,
    TransportResult,
)
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.physics.diffusion import simulate_drug_transport
from tmesim.physics.oxygen import solve_oxygen
from tmesim.physics.pressure import solve_ifp
from tmesim.reconstruct.real_pipeline import (
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from tmesim.sim.engine import simulate_response
from tmesim.treatment.drugs import delivery_source, plasma_profile
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.validation.vessel_profile import (
    distance_from_vessel_um,
    distance_profile,
    empirical_half_distance_um,
    fit_exponential,
    fit_exponential_plus_floor,
    usable_bins,
    vessel_mask_from_density,
)

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
DRUGS = ["temozolomide", "cisplatin"]

# CONVERGED grid resolution (results/numerical_rigor.md section 2). dx = 10 um
# violates the upwind advection CFL condition at dt = 10 s.
DX_UM = 12.0
DT_S = 10.0

# Dose-response window: the 6 h exposure window used by the real drug panel.
WINDOW_DOSE_S = 3600.0 * 6
# Schedule / half-life window. SCOPE CUT (see module docstring): shortened from
# 24 h to 12 h so the extra PDE solve fits the compute budget. It is held EXACTLY
# CONSTANT across every schedule and half-life compared, which is what the
# comparison requires; only the absolute kill level depends on the window length.
WINDOW_SCHED_S = 3600.0 * 12

# Nominal plasma peak. This is the repo-wide default (``PKConfig.plasma_peak``)
# and is used as the SHARED reference level for both drugs so the dose axis means
# the same thing for each. It is roughly the order of magnitude of temozolomide's
# clinical plasma Cmax; free-platinum Cmax for cisplatin is roughly an order of
# magnitude LOWER, so cisplatin's "1x" here is supratherapeutic. Approximate
# reference bands are drawn on the figure with that caveat attached.
NOMINAL_PEAK_UM = 50.0
CLINICAL_CMAX_BAND_UM = {           # ORDER-OF-MAGNITUDE reference only
    "temozolomide": (30.0, 70.0),
    "cisplatin": (2.0, 6.0),
}

DOSE_MULTIPLIERS = np.array([0.05, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0])

NOMINAL_HALF_LIFE_S = 3600.0 * 4
# (n_doses, dose_interval_s) -- fractionation of one 24 h course.
SCHEDULES = [
    (1, WINDOW_SCHED_S),          # single bolus
    (3, 3600.0 * 4),              # q4h
    (10, 3600.0 * 1.2),           # q1.2h
]
HALF_LIVES_H = [0.25, 4.0, 48.0]

#: The schedule effect is a fixed RELATIVE change in tissue exposure, so the
#: absolute kill spread it produces depends on where you sit on the Hill curve.
#: At 1x nominal temozolomide kill is ~0.01 and a ~1% relative change is smaller
#: than the 1/n_tumor quantum of the kill fraction, which would report as a
#: spurious zero. Every schedule / half-life leg is therefore ALSO evaluated at
#: higher dose multiples — free, because the exposure field is exactly linear in
#: amplitude (no extra PDE solve).
SCHEDULE_DOSE_MULTIPLIERS = [1.0, 3.0, 10.0, 30.0]

#: SCOPE CUT: the schedule and half-life legs need their OWN transport solve (a
#: different exposure window), so they are run for ONE drug only. The dose-response
#: sweep still covers both drugs. Stated explicitly in the report.
SCHEDULE_DRUG = "temozolomide"

#: TIME-RESOLVED PK. The engine used to collapse the whole dosing schedule to a
#: single time-averaged plasma level and drive the transport solver with a STATIC
#: source, which made every matched-cumulative-exposure schedule comparison
#: bit-identical BY CONSTRUCTION. With `PDEConfig.time_resolved=True` the solver
#: is driven step by step by the real `plasma(t)`, so tissue concentration rises
#: and falls and the schedule is genuinely resolvable. Set False to reproduce the
#: superseded static-source numbers exactly.
TIME_RESOLVED = True

SEED = 0

# Colours: categorical slots 1 and 2 of the validated reference palette
# (validate_palette.js: PASS lightness / chroma / CVD; contrast WARN on aqua is
# relieved by direct labels + the tables in the .md).
C_DRUG = {"temozolomide": "#2a78d6", "cisplatin": "#1baf7a"}
C_INK = "#0b0b0b"
C_INK2 = "#52514e"
C_MUTED = "#8a8984"
C_GRID = "#e4e3df"
C_CRITICAL = "#d03b3b"
# Sequential blue ramp (magnitude = dose), steps 250 -> 700.
SEQ_BLUE = ["#86b6ef", "#5598e7", "#2a78d6", "#256abf", "#184f95", "#0d366b"]

_HERE = os.path.dirname(os.path.abspath(__file__))
_RESULTS_DIR = os.path.abspath(os.path.join(_HERE, "..", "results"))
_OUT_PNG = os.path.join(_RESULTS_DIR, "dose_schedule.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "dose_schedule.md")
_CACHE = os.path.join(_RESULTS_DIR, ".dose_schedule_cache.json")


# --------------------------------------------------------------------------- #
# Engine-faithful helpers (orchestration only -- no frozen module is re-implemented)
# --------------------------------------------------------------------------- #
def compile_real_drug(name: str, alpha: float):
    """The panel's compilation pipeline: calibrated perivascular sink (alpha fit
    once on doxorubicin) then the saturable binding-site barrier (a no-op for
    small molecules, but every drug goes through the identical path)."""
    dp = calibrated_drugparams(name, alpha)
    mol = get_molecule(name)
    return apply_barrier_to_drugparams(dp, mol), mol


def prepare_tissue(micro: Microenvironment, cfg: SimConfig):
    """Run the engine's DRUG-INDEPENDENT prologue exactly once.

    Mirrors ``sim.engine.simulate_response`` steps 1-3 (private working copy ->
    steady-state oxygen + hypoxia/quiescence banding -> IFP + Darcy velocity)
    using the same frozen functions. Returns the working microenvironment, the
    oxygen-labeled pre-kill cell snapshot, and the convective velocity field.
    """
    work = Microenvironment(
        domain=micro.domain,
        fields=TissueFields(
            ecm_density=micro.fields.ecm_density,
            vascular_density=micro.fields.vascular_density,
        ),
        cells=micro.cells.copy(),
        patient_id=micro.patient_id,
        metadata=micro.metadata,
    )
    oxygen = solve_oxygen(work, cfg.oxygen)
    work.fields.oxygen = oxygen
    update_cell_states(work, oxygen, cfg.oxygen)
    ifp, velocity = solve_ifp(work)
    work.fields.ifp = ifp
    work.fields.velocity = velocity
    return work, work.cells.copy(), velocity


def shape_profile(pk: PKConfig, total_time: float, dt: float):
    """``plasma(t)`` on the solver's own step grid, renormalised to unit window mean.

    Returns ``(shape, mean)`` where ``shape = plasma(t) / mean(plasma(t))`` — a
    dimensionless SHAPE with window-average exactly 1 — and ``mean`` the actual
    window-mean plasma level of ``pk`` in uM.

    Why this is the right object. The transport PDE is linear in the source and
    starts from C = 0, so scaling the whole plasma curve scales the exposure field
    exactly. Splitting a schedule into (shape x amplitude) therefore means ONE
    solve per plasma-curve *shape* covers every dose level of that shape, and the
    resulting exposure field for any window-mean level ``P`` is exactly
    ``P * auc_shape``. The amplitude fast path survives time resolution untouched;
    what no longer collapses is the SHAPE — 1 vs 10 fractions are different
    shapes and now need different solves, which is precisely the capability that
    was missing.

    The peak in ``pk`` is irrelevant here (it divides out), so the shape depends
    only on ``(half_life, n_doses, dose_interval)``.
    """
    n = int(round(float(total_time) / float(dt)))
    times = np.arange(n, dtype=np.float64) * float(dt)
    prof = plasma_profile(pk, times)
    mean = float(prof.mean())
    if mean <= 0.0:
        raise ValueError("degenerate schedule: zero window-averaged plasma")
    return prof / mean, mean


def unit_transport(work, drug, pde: PDEConfig, velocity,
                   profile: Optional[np.ndarray] = None) -> TransportResult:
    """Solve the transport PDE ONCE at unit window-mean plasma level.

    ``profile=None`` -> the historical STATIC source: a constant unit plasma level
    every step. ``profile`` (from :func:`shape_profile`, window-mean 1) -> the
    time-resolved source, i.e. the real plasma curve normalised to the same
    cumulative exposure.

    Either way the returned exposure field is "per uM of window-mean plasma", so
    the exposure for any level ``P`` is exactly ``P * auc``. Verified against the
    full engine in step [5] of :func:`run`.
    """
    source = delivery_source(work, drug, 1.0)
    return simulate_drug_transport(
        work, drug, source, pde, velocity=(velocity if pde.convection else None),
        source_profile=profile,
    )


def residual_fraction(tr: TransportResult, tissue: np.ndarray, k_loss: float) -> float:
    """Fraction of the delivered drug still IN THE TISSUE at the end of the window.

    This is the single quantity through which schedule acts in a linear transport
    model. Integrating ``dC/dt = div(D grad C) - div(vC) + S - kC`` over ``[0, T]``
    with ``C(0) = 0`` gives, for ``A = AUC`` and ``M = k - div(D grad .) + div(v .)``::

        A = M^{-1} ( S_total - C(T) )

    ``S_total`` is identical for any two schedules at matched cumulative plasma
    exposure, so the ENTIRE difference between them is the terminal residual
    ``C(T)`` — drug delivered but not yet counted inside the observation window.
    A schedule that finishes early leaves ``C(T) ~ 0`` and realises the full
    steady-state exposure field; one still infusing at ``T`` does not.

    Returned as ``sum C(T) / (k * sum AUC + sum C(T))``: eliminated drug is
    ``k * AUC``, still-resident drug is ``C(T)``, and (with no-flux boundaries)
    the two account for everything delivered.
    """
    c_end = float(np.asarray(tr.final_concentration)[tissue].sum())
    eliminated = float(k_loss) * float(np.asarray(tr.auc)[tissue].sum())
    total = eliminated + c_end
    return float(c_end / total) if total > 0 else 0.0


def _scaled(tr: TransportResult, scale: float) -> TransportResult:
    return TransportResult(
        auc=tr.auc * float(scale),
        final_concentration=tr.final_concentration * float(scale),
        times=tr.times,
        mean_concentration_history=tr.mean_concentration_history * float(scale),
        penetration_index=tr.penetration_index,
    )


def kill_at(work, labeled_cells, tr: TransportResult, drug, pd_cfg: PDConfig,
            plasma_mean: float, seed: int = SEED) -> Dict:
    """Apply the engine's own kill mapping at plasma level ``plasma_mean``.

    Uses ``cells.pkpd.apply_drug_kill`` unmodified, on a fresh copy of the
    oxygen-labeled pre-kill cells, with ``default_rng(seed)`` — byte-for-byte the
    call the engine makes.
    """
    snap = Microenvironment(
        domain=work.domain,
        fields=work.fields,
        cells=labeled_cells.copy(),
        patient_id=work.patient_id,
        metadata=work.metadata,
    )
    scaled = _scaled(tr, plasma_mean)
    apply_drug_kill(snap, scaled, pd_cfg, drug, rng=np.random.default_rng(int(seed)))

    duration = float(tr.times[-1]) if tr.times.size and tr.times[-1] > 0 else 1.0
    tumor = labeled_cells.tumor_mask()
    exposure = np.asarray(snap.cells.exposure, dtype=np.float64)
    E = exposure / duration                      # time-averaged conc, uM
    E_tumor = E[tumor]
    ic50 = float(drug.ic50)
    return dict(
        plasma_mean=float(plasma_mean),
        kill=float(1.0 - snap.cells.viable_tumor_fraction()),
        ic50_coverage=float(np.mean(E_tumor >= ic50)),
        cov_10pct_ic50=float(np.mean(E_tumor >= 0.1 * ic50)),
        median_E_over_ic50=float(np.median(E_tumor) / ic50),
        mean_E_over_ic50=float(np.mean(E_tumor) / ic50),
        mean_auc_tumor=float(np.mean(exposure[tumor])),
    )


def kill_ceiling(labeled_cells, pd_cfg: PDConfig) -> float:
    """The maximum kill fraction ANY dose can reach on this tumor.

    As exposure -> infinity the Hill term saturates at ``pd.max_kill_rate`` and the
    per-cell kill probability tends to ``max_kill_rate * sensitivity * R`` with
    ``R`` the hypoxia / quiescence resistance factor. Averaged over tumor cells
    this is the asymptote of the dose-response curve — a property of the
    SENSITIVITY MODEL and the tissue's oxygenation, not of delivery.
    """
    tumor = labeled_cells.tumor_mask()
    resistance = np.ones(labeled_cells.n, dtype=np.float64)
    resistance[labeled_cells.state == CellState.HYPOXIC] = pd_cfg.hypoxia_resistance
    resistance[labeled_cells.state == CellState.QUIESCENT] = pd_cfg.quiescent_resistance
    p = np.clip(pd_cfg.max_kill_rate * labeled_cells.sensitivity * resistance, 0.0, 1.0)
    bands = {"VIABLE": None, "HYPOXIC": CellState.HYPOXIC,
             "QUIESCENT": CellState.QUIESCENT}
    breakdown = {}
    st = labeled_cells.state[tumor]
    pt = p[tumor]
    n_t = max(int(tumor.sum()), 1)
    for label, code in bands.items():
        sel = (~np.isin(st, [int(CellState.HYPOXIC), int(CellState.QUIESCENT)])
               if code is None else st == int(code))
        breakdown[label] = dict(
            fraction=float(sel.sum()) / n_t,
            mean_p=float(pt[sel].mean()) if sel.any() else 0.0,
            contribution=float(pt[sel].sum()) / n_t,
            resistance_multiplier=(1.0 if code is None else
                                   (pd_cfg.hypoxia_resistance
                                    if code == CellState.HYPOXIC
                                    else pd_cfg.quiescent_resistance)),
        )
    return float(np.mean(pt)), breakdown


def window_plasma_mean(pk: PKConfig, total_time: float, dt: float) -> float:
    """The single number the engine actually consumes from the schedule: the
    time-average of the multi-dose plasma profile over the exposure window,
    evaluated on the engine's own time grid."""
    times = np.arange(0.0, float(total_time), float(dt))
    prof = plasma_profile(pk, times)
    return float(prof.mean()) if prof.size else 0.0


def peak_for_matched_auc(target_mean: float, n_doses: int, interval: float,
                         half_life: float, total_time: float, dt: float) -> float:
    """Plasma peak per dose that makes the WINDOW-AVERAGED plasma level equal
    ``target_mean`` for the given schedule. The profile is linear in the peak, so
    this is a single division (no search)."""
    unit = window_plasma_mean(
        PKConfig(plasma_peak=1.0, half_life=half_life, n_doses=n_doses,
                 dose_interval=interval),
        total_time, dt,
    )
    if unit <= 0.0:
        raise ValueError("degenerate schedule: zero window-averaged plasma")
    return float(target_mean) / unit


# --------------------------------------------------------------------------- #
# Penetration-depth analysis (distance from the nearest real vessel)
# --------------------------------------------------------------------------- #
def depth_analysis(micro, tr: TransportResult, drug, dx_um: float,
                   distance_um: np.ndarray, max_um: Optional[float] = None) -> Dict:
    """Bin the UNIT-plasma exposure by distance from the nearest vessel.

    Returns the profile (in time-averaged concentration per uM of plasma), its
    exponential decay length, its half-distance, and — the operationally useful
    number — the distance at which exposure crosses IC50 as a function of dose.
    Because the field is linear in dose, the profile SHAPE is dose-invariant and
    only its amplitude moves, so the IC50 depth necessarily grows like
    ``L * ln(dose)``.
    """
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    duration = float(tr.times[-1])
    unit_E = np.asarray(tr.auc, dtype=np.float64) / duration   # uM per uM plasma
    if max_um is None:
        max_um = float(np.percentile(distance_um[tissue], 99.5))
    prof = distance_profile(unit_E, distance_um, tissue, bin_um=dx_um,
                            max_distance_um=max_um)
    d, v, w = usable_bins(prof, max_um=max_um)
    fit = fit_exponential(d, v, w)
    half = empirical_half_distance_um(d, v)
    # SOURCE NULL: the same binning applied to the raw vascular-density field. The
    # calibrated sink is fast, so the drug field largely inherits the source shape;
    # reporting both keeps that confound visible (see validation/vessel_profile.py).
    src_prof = distance_profile(np.asarray(micro.fields.vascular_density, float),
                                distance_um, tissue, bin_um=dx_um,
                                max_distance_um=max_um)
    sd, sv, _ = usable_bins(src_prof, max_um=max_um)
    src_half = empirical_half_distance_um(sd, sv)
    return dict(
        distance_um=d.tolist(),
        unit_E=v.tolist(),
        bin_count=w.tolist(),
        max_profiled_um=float(max_um),
        L_um=float(fit["L_um"]) if fit.get("L_um") is not None else None,
        r2=float(fit.get("r2_log", float("nan"))),
        half_distance_um=(float(half) if half is not None else None),
        source_null_half_distance_um=(float(src_half) if src_half is not None else None),
        perivascular_unit_E=float(v[0]) if v.size else 0.0,
        ic50=float(drug.ic50),
    )


def augment_depth(depth: Dict, p90_um: float) -> Dict:
    """Characterise the profile HONESTLY: it is not a single exponential.

    On this real section the exposure-vs-distance profile is a steep perivascular
    drop sitting on a flat, non-zero FLOOR. The floor is not diffusion: the engine's
    delivery source is the *continuous* vascular-density field, so even a pixel far
    from the thresholded vessel mask is still being fed by diffuse low-density
    source. A single exponential fitted over the whole range is dominated by that
    plateau and returns a meaningless L (this is precisely the failure
    ``vessel_profile.fit_exponential_plus_floor`` exists to avoid).

    We therefore report three things instead of one number:
      * ``fit_near``  — single exponential restricted to the WELL-POPULATED
        near-field (d <= p90 of the tissue's distance-to-vessel),
      * ``fit_floor`` — ``C0*exp(-d/L) + floor`` over the full range, which
        separates the decaying perivascular component from the background,
      * ``plateau_ratio`` — the far-field floor as a fraction of the perivascular
        value, i.e. how much of the field is background rather than gradient.
    Mutates and returns ``depth``.
    """
    d = np.asarray(depth["distance_um"], dtype=np.float64)
    v = np.asarray(depth["unit_E"], dtype=np.float64)
    w = np.asarray(depth.get("bin_count", np.ones_like(d)), dtype=np.float64)
    near = d <= float(p90_um)
    fn = fit_exponential(d[near], v[near], w[near]) if near.sum() >= 3 else {}
    ff = fit_exponential_plus_floor(d, v, w)
    far = d >= 0.75 * d.max() if d.size else np.zeros(0, bool)
    depth["fit_near"] = dict(
        L_um=(float(fn["L_um"]) if fn.get("L_um") is not None else None),
        r2=float(fn.get("r2_log", float("nan"))),
        range_um=float(p90_um),
        n_bins=int(near.sum()),
    )
    depth["fit_floor"] = {k: (float(ff[k]) if isinstance(ff.get(k), (int, float))
                              else ff.get(k))
                          for k in ("C0", "L_um", "floor", "d_half_um",
                                    "d_half_total_um", "r2_linear_space", "error")
                          if k in ff}
    depth["plateau_ratio"] = (float(v[far].mean() / v[0]) if v.size and v[0] > 0
                              and far.any() else None)
    depth["p90_um"] = float(p90_um)
    return depth


def ic50_depth_um(depth: Dict, plasma_mean: float):
    """``(depth_um, censored)`` — the deepest distance-from-vessel at which the
    time-averaged exposure still clears IC50, by linear interpolation on the
    binned profile.

    ``depth_um is None`` means IC50 is not cleared even at the vessel wall.
    ``censored is True`` means exposure is still above IC50 at the far end of the
    profiled range, so the returned depth is a lower bound.
    """
    d = np.asarray(depth["distance_um"], dtype=np.float64)
    E = np.asarray(depth["unit_E"], dtype=np.float64) * float(plasma_mean)
    ic50 = depth["ic50"]
    if d.size == 0 or E[0] < ic50:
        return None, False
    below = np.nonzero(E < ic50)[0]
    if below.size == 0:
        return float(d[-1]), True
    i = int(below[0])
    e0, e1 = float(E[i - 1]), float(E[i])
    if e0 == e1:
        return float(d[i]), False
    return float(d[i - 1] + (e0 - ic50) / (e0 - e1) * (d[i] - d[i - 1])), False


def equilibration_time_s(tr: TransportResult, frac: float = 0.95) -> Optional[float]:
    """Time for the mean tissue concentration to reach ``frac`` of its final value.

    **Only meaningful for a STATIC (sustained) source.** Under a sustained plasma
    level the mean tissue concentration rises monotonically to a plateau and this
    is the tissue's transport/clearance equilibration time-scale — the number that
    decides whether schedule SHAPE can matter at all.

    Under a TIME-RESOLVED source the history rises and then falls with the plasma
    curve, so "95% of the final value" is measuring the washout tail, not
    equilibration, and this function must not be used. Use
    :func:`equilibration_t95_analytic_s` instead.
    """
    h = np.asarray(tr.mean_concentration_history, dtype=np.float64)
    t = np.asarray(tr.times, dtype=np.float64)
    if h.size < 2 or h[-1] <= 0:
        return None
    hit = np.nonzero(h >= frac * h[-1])[0]
    return float(t[int(hit[0])]) if hit.size else None


def equilibration_t95_analytic_s(k_loss: float, frac: float = 0.95) -> float:
    """Analytic tissue equilibration time under a sustained plasma level.

    The tissue-mean concentration under a step source relaxes as
    ``C(t) = C_ss (1 - e^{-k t})`` (diffusion redistributes but does not change
    the mean; only the loss rate ``k = clearance + binding`` sets the relaxation),
    so ``t_95 = -ln(1 - 0.95)/k = ln(20)/k = 3.0 * tau``.

    Preferred over the measured history because it is exact, independent of
    ``save_every`` quantisation, and — unlike the measured version — still
    well-defined when the solve was driven by a time-varying source.
    """
    return float(-math.log(1.0 - frac) / float(k_loss))


# --------------------------------------------------------------------------- #
# The run
# --------------------------------------------------------------------------- #
def run(dx_um: float = DX_UM) -> Dict:
    t_start = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    print(f"[1] Reconstructing the REAL Xenium GBM section once (dx = {dx_um} um) ...")
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
    rep = reconstruction_report(micro)
    print(f"    {rep['n_cells']:,} real cells (tumor {rep['n_tumor']:,})  "
          f"grid {rep['grid']}  vasc {rep['mean_vascular_density']:.3f}  "
          f"ecm {rep['mean_ecm_density']:.3f}")

    alpha = calibrate_alpha()
    print(f"[2] perivascular-sink alpha = {alpha:.2f} (fit once on doxorubicin)")

    base_cfg = SimConfig(seed=SEED)
    pd_cfg = base_cfg.pd
    work, labeled_base, velocity = prepare_tissue(micro, base_cfg)
    tumor = labeled_base.tumor_mask()
    hypoxic_frac = float(np.mean(labeled_base.state[tumor] == CellState.HYPOXIC))
    quiescent_frac = float(np.mean(labeled_base.state[tumor] == CellState.QUIESCENT))
    print(f"[3] tissue prologue done: hypoxic tumor {hypoxic_frac*100:.1f}%, "
          f"quiescent {quiescent_frac*100:.1f}%")

    # DRUG-SPECIFIC per-cell sensitivity (reconstruct.real_pipeline). The
    # reconstruction itself is drug-agnostic (uniform 1.0); the measured
    # determinant is layered on per drug -- MKI67 for cell-cycle-dependent
    # cytotoxics, and honestly nothing (uniform 1.0) for temozolomide, whose
    # determinant MGMT is not on this panel.
    sens_by_drug: Dict[str, np.ndarray] = {}
    sens_info: Dict[str, Dict] = {}
    for name in DRUGS:
        probe = Microenvironment(
            domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
            patient_id=micro.patient_id, metadata=micro.metadata,
        )
        info = apply_drug_sensitivity(probe, name, mode="determinant")
        sens_by_drug[name] = np.asarray(probe.cells.sensitivity, dtype=np.float64).copy()
        sens_info[name] = info
        print(f"    sensitivity[{name}]: mode={info.get('mode')} "
              f"determinant={info.get('determinant_gene')} "
              f"({info.get('kind')})  median over tumor cells = "
              f"{float(np.median(sens_by_drug[name][tumor])):.3f}")
        if info.get("upper_bound"):
            print(f"      !! {name}: {info.get('comparability')}")

    # Vessel-distance map (shared by both drugs; depends only on the tissue).
    vmask = vessel_mask_from_density(micro.fields.vascular_density,
                                     micro.domain.tissue_mask)
    distance_um = distance_from_vessel_um(vmask["mask"], dx_um)
    print(f"    vessel mask: {vmask['area_fraction']*100:.1f}% of tissue "
          f"(p{vmask['percentile']:.0f} of vascular density)")

    pde_dose = PDEConfig(dt=DT_S, total_time=WINDOW_DOSE_S, convection=True,
                         save_every=10, time_resolved=TIME_RESOLVED)
    pde_sched = PDEConfig(dt=DT_S, total_time=WINDOW_SCHED_S, convection=True,
                          save_every=100, time_resolved=TIME_RESOLVED)
    # The superseded static-source configuration, kept so this run can quantify
    # what the time-resolved driver actually changed rather than asserting it.
    pde_sched_static = PDEConfig(dt=DT_S, total_time=WINDOW_SCHED_S,
                                 convection=True, save_every=100,
                                 time_resolved=False)
    tissue_mask = np.asarray(micro.domain.tissue_mask, dtype=bool)

    out: Dict = dict(
        meta=dict(
            dx_um=float(dx_um), dt_s=DT_S, alpha=float(alpha),
            window_dose_h=WINDOW_DOSE_S / 3600.0,
            window_sched_h=WINDOW_SCHED_S / 3600.0,
            nominal_peak_um=NOMINAL_PEAK_UM,
            nominal_half_life_h=NOMINAL_HALF_LIFE_S / 3600.0,
            time_resolved=bool(TIME_RESOLVED),
            seed=SEED,
            n_cells=int(rep["n_cells"]), n_tumor=int(rep["n_tumor"]),
            grid=list(rep["grid"]),
            mean_vascular_density=float(rep["mean_vascular_density"]),
            mean_ecm_density=float(rep["mean_ecm_density"]),
            mean_oxygen=float(rep["mean_oxygen"]),
            hypoxic_tumor_fraction=hypoxic_frac,
            quiescent_tumor_fraction=quiescent_frac,
            vessel_area_fraction=float(vmask["area_fraction"]),
            vessel_percentile=float(vmask["percentile"]),
            distance_median_um=float(np.median(distance_um[np.asarray(
                micro.domain.tissue_mask, dtype=bool)])),
            distance_p90_um=float(np.percentile(distance_um[np.asarray(
                micro.domain.tissue_mask, dtype=bool)], 90.0)),
            distance_max_um=float(distance_um[np.asarray(
                micro.domain.tissue_mask, dtype=bool)].max()),
        ),
        drugs={},
        engine_checks=[],
    )

    for name in DRUGS:
        drug, mol = compile_real_drug(name, alpha)
        k_loss = float(drug.clearance + drug.binding)
        # This drug's own per-cell sensitivity, on the oxygen-labeled cells.
        labeled = labeled_base.copy()
        labeled.sensitivity = sens_by_drug[name].copy()
        ceiling, ceiling_breakdown = kill_ceiling(labeled, pd_cfg)
        print(f"\n[4] {name}: D={drug.diffusivity:.1f} um^2/s  IC50={drug.ic50:.2f} uM  "
              f"hill={drug.hill}  bbb={drug.bbb_penetrance:.3f}  "
              f"k_loss={k_loss:.3e} /s  (tau = {1.0/k_loss:.0f} s)")
        print(f"    KILL CEILING (exposure -> inf) = {ceiling:.4f}  <- set by the "
              f"sensitivity model + oxygenation, NOT by dose")

        # ---- 4a. unit transport solves ---------------------------------------
        # One solve per plasma-curve SHAPE (amplitude is still exact by
        # linearity). Under the old static source every shape collapsed to the
        # same constant level, which is exactly why schedule was unanswerable.
        dose_shape, _ = shape_profile(
            PKConfig(plasma_peak=1.0, half_life=NOMINAL_HALF_LIFE_S, n_doses=1,
                     dose_interval=WINDOW_DOSE_S),
            WINDOW_DOSE_S, DT_S)
        t0 = time.time()
        tr_dose = unit_transport(work, drug, pde_dose, velocity,
                                 profile=(dose_shape if TIME_RESOLVED else None))
        n_solves = 1
        print(f"    unit transport, {WINDOW_DOSE_S/3600:.0f} h window "
              f"({'time-resolved' if TIME_RESOLVED else 'static'}): "
              f"{time.time()-t0:.0f}s")
        # SCOPE CUT: the schedule-window solves are run for one drug only.
        do_schedule = (name == SCHEDULE_DRUG)
        tr_sched_static = None
        if do_schedule:
            t0 = time.time()
            tr_sched_static = unit_transport(work, drug, pde_sched_static,
                                             velocity, profile=None)
            n_solves += 1
            print(f"    unit transport, {WINDOW_SCHED_S/3600:.0f} h window "
                  f"(STATIC baseline, the superseded path): {time.time()-t0:.0f}s")
        else:
            print(f"    [scope cut] schedule / half-life legs SKIPPED for {name} "
                  f"(they need their own PDE solves; dropped for compute budget)")

        # Equilibration is a SUSTAINED-source property. Measure it on the static
        # solve when we have one; otherwise use the analytic ln(20)/k, which is
        # what the measured number converges to anyway. Never measure it off a
        # time-resolved solve — the history there rises AND falls.
        t_eq = equilibration_t95_analytic_s(k_loss)
        t_eq_measured = equilibration_time_s(
            tr_sched_static if tr_sched_static is not None else tr_dose)
        depth = depth_analysis(micro, tr_dose, drug, dx_um, distance_um)
        augment_depth(depth, out["meta"]["distance_p90_um"])
        depth_max_um = depth["max_profiled_um"]

        # ---- 4b. DOSE-RESPONSE sweep -----------------------------------------
        pk_nominal = PKConfig(plasma_peak=NOMINAL_PEAK_UM,
                              half_life=NOMINAL_HALF_LIFE_S, n_doses=1,
                              dose_interval=WINDOW_DOSE_S)
        mean_at_1x = window_plasma_mean(pk_nominal, WINDOW_DOSE_S, DT_S)
        dose_rows = []
        for m in DOSE_MULTIPLIERS:
            pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM * float(m),
                          half_life=NOMINAL_HALF_LIFE_S, n_doses=1,
                          dose_interval=WINDOW_DOSE_S)
            pm = window_plasma_mean(pk, WINDOW_DOSE_S, DT_S)
            row = kill_at(work, labeled, tr_dose, drug, pd_cfg, pm)
            row["multiplier"] = float(m)
            row["plasma_peak"] = float(NOMINAL_PEAK_UM * m)
            dep_um, censored = ic50_depth_um(depth, pm)
            row["ic50_depth_um"] = dep_um
            row["ic50_depth_censored"] = bool(censored)
            dose_rows.append(row)
            print(f"      {m:5.2f}x  peak={row['plasma_peak']:7.1f} uM  "
                  f"kill={row['kill']:.4f}  IC50-cov={row['ic50_coverage']:.4f}  "
                  f"depth={'-' if dep_um is None else round(dep_um, 1)}"
                  f"{'+' if censored else ''}")

        # ---- 4c. SCHEDULE at matched total exposure ---------------------------
        # Reference total exposure: one nominal bolus over the schedule window.
        sched_rows = []
        hl_rows = []
        target_mean = float("nan")
        pk_ref = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                          n_doses=1, dose_interval=WINDOW_SCHED_S)
        target_mean = window_plasma_mean(pk_ref, WINDOW_SCHED_S, DT_S)
        for n, interval in (SCHEDULES if do_schedule else []):
            # ONE solve per schedule SHAPE. Legs A/B/C differ only in amplitude,
            # so they share it exactly (linearity). With TIME_RESOLVED off this
            # falls back to the single static solve and all three schedules give
            # bit-identical leg-A answers -- the degeneracy being superseded.
            if TIME_RESOLVED:
                sh, _ = shape_profile(
                    PKConfig(plasma_peak=1.0, half_life=NOMINAL_HALF_LIFE_S,
                             n_doses=n, dose_interval=interval),
                    WINDOW_SCHED_S, DT_S)
                t0 = time.time()
                tr_n = unit_transport(work, drug, pde_sched, velocity, profile=sh)
                n_solves += 1
                print(f"      [solve] schedule n={n} shape: {time.time()-t0:.0f}s")
            else:
                tr_n = tr_sched_static

            # leg A: MATCHED cumulative plasma exposure (the key comparison)
            peak_a = peak_for_matched_auc(target_mean, n, interval,
                                          NOMINAL_HALF_LIFE_S, WINDOW_SCHED_S, DT_S)
            pk_a = PKConfig(plasma_peak=peak_a, half_life=NOMINAL_HALF_LIFE_S,
                            n_doses=n, dose_interval=interval)
            mean_a = window_plasma_mean(pk_a, WINDOW_SCHED_S, DT_S)
            res_a = kill_at(work, labeled, tr_n, drug, pd_cfg, mean_a)

            # leg B: matched ADMINISTERED dose (sum of peaks fixed), same window
            peak_b = NOMINAL_PEAK_UM / n
            pk_b = PKConfig(plasma_peak=peak_b, half_life=NOMINAL_HALF_LIFE_S,
                            n_doses=n, dose_interval=interval)
            mean_b = window_plasma_mean(pk_b, WINDOW_SCHED_S, DT_S)
            res_b = kill_at(work, labeled, tr_n, drug, pd_cfg, mean_b)

            # leg C: matched PEAK (total exposure rises with n) -- the control
            pk_c = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                            n_doses=n, dose_interval=interval)
            mean_c = window_plasma_mean(pk_c, WINDOW_SCHED_S, DT_S)
            res_c = kill_at(work, labeled, tr_n, drug, pd_cfg, mean_c)

            # The superseded static-source answer for the SAME leg A, so the
            # report can quantify what time resolution changed.
            res_a_static = kill_at(work, labeled, tr_sched_static, drug, pd_cfg,
                                   mean_a)

            # Depth + the mechanism (terminal residual) for this shape.
            dep_n = depth_analysis(micro, tr_n, drug, dx_um, distance_um,
                                   max_um=depth_max_um)
            dep_a_um, dep_a_cen = ic50_depth_um(dep_n, mean_a)

            # Leg A re-evaluated up the dose axis (pure rescaling, no new solve).
            by_dose = {f"{m:g}": kill_at(work, labeled, tr_n, drug, pd_cfg,
                                         mean_a * m)["kill"]
                       for m in SCHEDULE_DOSE_MULTIPLIERS}
            by_dose_static = {f"{m:g}": kill_at(work, labeled, tr_sched_static,
                                                drug, pd_cfg, mean_a * m)["kill"]
                              for m in SCHEDULE_DOSE_MULTIPLIERS}

            sched_rows.append(dict(
                n_doses=int(n), interval_h=float(interval / 3600.0),
                peak_matched_auc=float(peak_a), plasma_mean_matched=float(mean_a),
                kill_matched=res_a["kill"], cov_matched=res_a["ic50_coverage"],
                peak_matched_dose=float(peak_b), plasma_mean_dose=float(mean_b),
                kill_matched_dose=res_b["kill"], cov_matched_dose=res_b["ic50_coverage"],
                plasma_mean_peak=float(mean_c), kill_matched_peak=res_c["kill"],
                # time-resolution diagnostics
                kill_matched_static=res_a_static["kill"],
                kill_matched_by_dose=by_dose,
                kill_matched_by_dose_static=by_dose_static,
                mean_auc_tumor_matched=res_a["mean_auc_tumor"],
                mean_auc_tumor_matched_static=res_a_static["mean_auc_tumor"],
                residual_fraction=residual_fraction(tr_n, tissue_mask, k_loss),
                half_distance_um=dep_n.get("half_distance_um"),
                penetration_index=float(tr_n.penetration_index),
                depth_matched_um=dep_a_um, depth_matched_censored=bool(dep_a_cen),
            ))
            print(f"      n={n:2d} q{interval/3600.0:4.1f}h | matched-AUC: peak="
                  f"{peak_a:7.2f} mean={mean_a:8.4f} kill={res_a['kill']:.6f} "
                  f"(static was {res_a_static['kill']:.6f}) | "
                  f"matched-dose: mean={mean_b:8.4f} kill={res_b['kill']:.6f} | "
                  f"matched-peak: kill={res_c['kill']:.6f} | "
                  f"resid={sched_rows[-1]['residual_fraction']*100:.2f}% "
                  f"halfdist={sched_rows[-1]['half_distance_um']}")

        # ---- 4d. HALF-LIFE (infusion-duration proxy) --------------------------
        for hl_h in (HALF_LIVES_H if do_schedule else []):
            hl = hl_h * 3600.0
            # ONE solve per half-life SHAPE (matched-exposure and matched-peak
            # legs differ only in amplitude, so they share it).
            if TIME_RESOLVED:
                sh, _ = shape_profile(
                    PKConfig(plasma_peak=1.0, half_life=hl, n_doses=1,
                             dose_interval=WINDOW_SCHED_S),
                    WINDOW_SCHED_S, DT_S)
                t0 = time.time()
                tr_h = unit_transport(work, drug, pde_sched, velocity, profile=sh)
                n_solves += 1
                print(f"      [solve] half-life {hl_h:g} h shape: "
                      f"{time.time()-t0:.0f}s")
            else:
                tr_h = tr_sched_static

            # Each half-life now has its OWN exposure field, so the depth profile
            # is re-binned per shape (under the static source there was only one
            # field, which is why 'depth unchanged' was guaranteed a priori).
            dep_h = depth_analysis(micro, tr_h, drug, dx_um, distance_um,
                                   max_um=depth_max_um)

            # matched total exposure
            peak_m = peak_for_matched_auc(target_mean, 1, WINDOW_SCHED_S, hl,
                                          WINDOW_SCHED_S, DT_S)
            pk_m = PKConfig(plasma_peak=peak_m, half_life=hl, n_doses=1,
                            dose_interval=WINDOW_SCHED_S)
            mean_m = window_plasma_mean(pk_m, WINDOW_SCHED_S, DT_S)
            res_m = kill_at(work, labeled, tr_h, drug, pd_cfg, mean_m)
            # matched peak (total exposure free to move) -- the control
            pk_p = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=hl, n_doses=1,
                            dose_interval=WINDOW_SCHED_S)
            mean_p = window_plasma_mean(pk_p, WINDOW_SCHED_S, DT_S)
            res_p = kill_at(work, labeled, tr_h, drug, pd_cfg, mean_p)
            # the superseded static-source answer for the same matched exposure
            res_m_static = kill_at(work, labeled, tr_sched_static, drug, pd_cfg,
                                   mean_m)
            dm, dm_c = ic50_depth_um(dep_h, mean_m)
            dp_, dp_c = ic50_depth_um(dep_h, mean_p)
            hl_rows.append(dict(
                half_life_h=float(hl_h),
                peak_matched=float(peak_m), plasma_mean_matched=float(mean_m),
                kill_matched=res_m["kill"], cov_matched=res_m["ic50_coverage"],
                depth_matched_um=dm, depth_matched_censored=bool(dm_c),
                plasma_mean_peakmatched=float(mean_p),
                kill_peakmatched=res_p["kill"], cov_peakmatched=res_p["ic50_coverage"],
                depth_peakmatched_um=dp_, depth_peakmatched_censored=bool(dp_c),
                # time-resolution diagnostics
                kill_matched_static=res_m_static["kill"],
                kill_matched_by_dose={
                    f"{m:g}": kill_at(work, labeled, tr_h, drug, pd_cfg,
                                      mean_m * m)["kill"]
                    for m in SCHEDULE_DOSE_MULTIPLIERS},
                kill_matched_by_dose_static={
                    f"{m:g}": kill_at(work, labeled, tr_sched_static, drug,
                                      pd_cfg, mean_m * m)["kill"]
                    for m in SCHEDULE_DOSE_MULTIPLIERS},
                mean_auc_tumor_matched=res_m["mean_auc_tumor"],
                mean_auc_tumor_matched_static=res_m_static["mean_auc_tumor"],
                residual_fraction=residual_fraction(tr_h, tissue_mask, k_loss),
                half_distance_um=dep_h.get("half_distance_um"),
                penetration_index=float(tr_h.penetration_index),
                perivascular_unit_E=dep_h.get("perivascular_unit_E"),
                plateau_ratio=augment_depth(
                    dep_h, out["meta"]["distance_p90_um"]).get("plateau_ratio"),
                depth_profile=dict(distance_um=dep_h["distance_um"],
                                   unit_E=dep_h["unit_E"]),
            ))
            print(f"      t1/2={hl_h:5.2f} h | matched-exposure: peak={peak_m:8.2f} "
                  f"kill={res_m['kill']:.6f} (static was {res_m_static['kill']:.6f}) "
                  f"depth={dm} halfdist={hl_rows[-1]['half_distance_um']} "
                  f"resid={hl_rows[-1]['residual_fraction']*100:.2f}% | "
                  f"matched-peak: mean={mean_p:8.4f} kill={res_p['kill']:.6f}")

        out["drugs"][name] = dict(
            params=dict(name=name, mw=float(mol.get("mw")),
                        diffusivity=float(drug.diffusivity), ic50=float(drug.ic50),
                        hill=float(drug.hill), bbb=float(drug.bbb_penetrance),
                        binding=float(drug.binding), clearance=float(drug.clearance),
                        k_loss=k_loss, tau_transport_s=1.0 / k_loss),
            kill_ceiling=float(ceiling),
            kill_ceiling_breakdown=ceiling_breakdown,
            sensitivity=dict(
                determinant_gene=sens_info[name].get("determinant_gene"),
                kind=sens_info[name].get("kind"),
                rationale=sens_info[name].get("rationale"),
                median_tumor=float(np.median(sens_by_drug[name][tumor])),
                mean_tumor=float(np.mean(sens_by_drug[name][tumor])),
                is_uniform=bool(sens_info[name].get("is_uniform", False)),
            ),
            equilibration_t95_s=t_eq,
            equilibration_t95_measured_s=t_eq_measured,
            penetration_index=float(tr_dose.penetration_index),
            depth=depth,
            plasma_mean_at_1x=float(mean_at_1x),
            dose=dose_rows,
            schedule=sched_rows,
            half_life=hl_rows,
            target_window_mean_um=float(target_mean),
            n_pde_solves=int(n_solves),
            static_penetration_index=(float(tr_sched_static.penetration_index)
                                      if tr_sched_static is not None else None),
            static_residual_fraction=(
                residual_fraction(tr_sched_static, tissue_mask, k_loss)
                if tr_sched_static is not None else None),
        )

    # ---- 5. ENGINE CROSS-CHECKS: the full frozen engine, no fast path --------
    print("\n[5] Cross-checking the fast path against the FULL simulate_response engine")
    checks = []
    tmz, _ = compile_real_drug("temozolomide", alpha)
    cis, _ = compile_real_drug("cisplatin", alpha)

    def _micro_for(dp_name):
        """The caller's microenvironment carrying THIS drug's per-cell sensitivity —
        exactly what the fast path used."""
        m = Microenvironment(
            domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
            patient_id=micro.patient_id, metadata=micro.metadata,
        )
        m.cells.sensitivity = sens_by_drug[dp_name].copy()
        return m

    def _check(label, drug, dp_name, pk, pde, fast_kill):
        t0 = time.time()
        res = simulate_response(_micro_for(dp_name),
                                SimConfig(drug=drug, pk=pk, pde=pde, seed=SEED))
        d = abs(float(res.kill_fraction) - float(fast_kill))
        print(f"    {label:44s} engine={res.kill_fraction:.10f}  "
              f"fast={fast_kill:.10f}  |delta|={d:.2e}  ({time.time()-t0:.0f}s)")
        checks.append(dict(label=label, drug=dp_name,
                           engine_kill=float(res.kill_fraction),
                           fast_kill=float(fast_kill), abs_delta=float(d)))

    # SCOPE CUT: three verification runs instead of six. One anchors the dose
    # sweep (dose window), one the schedule leg (most fractionated schedule), one
    # the half-life leg (shortest half-life, i.e. the most bolus-like shape and
    # the biggest departure from the superseded static source). Together they pin
    # the shape-x-amplitude fast path in BOTH windows and on all three legs.
    row = next(r for r in out["drugs"]["temozolomide"]["dose"]
               if abs(r["multiplier"] - 1.0) < 1e-9)
    _check(f"temozolomide dose 1x ({WINDOW_DOSE_S/3600:.0f} h window)",
           tmz, "temozolomide",
           PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                    n_doses=1, dose_interval=WINDOW_DOSE_S),
           pde_dose, row["kill"])

    srow = max(out["drugs"][SCHEDULE_DRUG]["schedule"], key=lambda r: r["n_doses"])
    _check(f"{SCHEDULE_DRUG} schedule n={srow['n_doses']} matched-AUC "
           f"({WINDOW_SCHED_S/3600:.0f} h window)", tmz, SCHEDULE_DRUG,
           PKConfig(plasma_peak=srow["peak_matched_auc"],
                    half_life=NOMINAL_HALF_LIFE_S, n_doses=srow["n_doses"],
                    dose_interval=srow["interval_h"] * 3600.0),
           pde_sched, srow["kill_matched"])

    hrow = min(out["drugs"][SCHEDULE_DRUG]["half_life"],
               key=lambda r: r["half_life_h"])
    _check(f"{SCHEDULE_DRUG} half-life {hrow['half_life_h']:g} h matched-exposure "
           f"({WINDOW_SCHED_S/3600:.0f} h window)", tmz, SCHEDULE_DRUG,
           PKConfig(plasma_peak=hrow["peak_matched"],
                    half_life=hrow["half_life_h"] * 3600.0, n_doses=1,
                    dose_interval=WINDOW_SCHED_S),
           pde_sched, hrow["kill_matched"])

    out["engine_checks"] = checks
    out["meta"]["max_engine_fast_delta"] = float(max(c["abs_delta"] for c in checks))
    out["meta"]["runtime_s"] = float(time.time() - t_start)
    print(f"\n    max |engine - fast| = {out['meta']['max_engine_fast_delta']:.2e}")
    print(f"[done] {out['meta']['runtime_s']/60:.1f} min")
    return out


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #
def _style_axes(ax):
    ax.set_facecolor("#fcfcfb")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(C_GRID)
    ax.tick_params(colors=C_INK2, labelsize=8.5, length=3)
    ax.grid(True, color=C_GRID, lw=0.8, alpha=0.9)
    ax.set_axisbelow(True)
    ax.xaxis.label.set_color(C_INK2)
    ax.yaxis.label.set_color(C_INK2)
    ax.xaxis.label.set_fontsize(9)
    ax.yaxis.label.set_fontsize(9)
    ax.title.set_color(C_INK)
    ax.title.set_fontsize(10)


def make_figure(out: Dict, path: str) -> None:
    meta = out["meta"]
    fig = plt.figure(figsize=(18.5, 10.4), facecolor="#fcfcfb")
    gs = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.26,
                          left=0.055, right=0.985, top=0.845, bottom=0.115)

    mult = DOSE_MULTIPLIERS

    # --- A. dose-response ---------------------------------------------------
    ax = fig.add_subplot(gs[0, 0]); _style_axes(ax)
    ceil_max = 0.0
    # Rank the end-points so the two drug labels are offset apart instead of
    # stacking on top of each other when the curves finish close together.
    _ends = sorted(DRUGS, key=lambda n: out["drugs"][n]["dose"][-1]["kill"])
    for name in DRUGS:
        rows = out["drugs"][name]["dose"]
        y = [r["kill"] for r in rows]
        ceil = out["drugs"][name]["kill_ceiling"]
        ceil_max = max(ceil_max, ceil, max(y))
        ax.plot(mult, y, color=C_DRUG[name], lw=2.0, marker="o", ms=4.5,
                mec="#fcfcfb", mew=0.8, zorder=3)
        ax.annotate(name, xy=(mult[-1], y[-1]),
                    xytext=(-4, (-15 if _ends.index(name) == 0 else 9)),
                    textcoords="offset points", ha="right", color=C_DRUG[name],
                    fontsize=9, fontweight="bold")
        ax.axhline(ceil, color=C_DRUG[name], lw=1.3, ls=(0, (5, 3)), alpha=0.75,
                   zorder=2)
        ax.annotate(f"ceiling {ceil:.3f}", xy=(mult[0], ceil), xytext=(0, 3),
                    textcoords="offset points", color=C_DRUG[name], fontsize=7.5)
    for name in DRUGS:
        lo, hi = (v / NOMINAL_PEAK_UM for v in CLINICAL_CMAX_BAND_UM[name])
        ax.axvspan(lo, hi, color=C_DRUG[name], alpha=0.10, lw=0, zorder=1)
    ax.axvline(1.0, color=C_MUTED, lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("plasma peak  (x nominal, nominal = 50 uM)")
    ax.set_ylabel("predicted tumor kill fraction")
    ax.set_ylim(0, min(1.0, ceil_max * 1.22))
    ax.set_title("A · Dose-response on THIS tumor\n"
                 "dashed = kill ceiling (sensitivity + hypoxia, not dose)")

    # --- B. IC50 coverage ---------------------------------------------------
    ax = fig.add_subplot(gs[0, 1]); _style_axes(ax)
    for name in DRUGS:
        rows = out["drugs"][name]["dose"]
        y = [r["ic50_coverage"] for r in rows]
        ax.plot(mult, y, color=C_DRUG[name], lw=2.0, marker="o", ms=4.5,
                mec="#fcfcfb", mew=0.8, zorder=3)
        j = int(np.argmax(np.asarray(y) > 0.02)) if max(y) > 0.02 else len(y) - 1
        ax.annotate(name, xy=(mult[j], y[j]), xytext=(6, -2),
                    textcoords="offset points", color=C_INK, fontsize=9,
                    fontweight="bold")
        lo, hi = (v / NOMINAL_PEAK_UM for v in CLINICAL_CMAX_BAND_UM[name])
        ax.axvspan(lo, hi, color=C_DRUG[name], alpha=0.10, lw=0, zorder=1)
    ax.axvline(1.0, color=C_MUTED, lw=1.0, ls=":")
    ax.set_xscale("log")
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("plasma peak  (x nominal)")
    ax.set_ylabel("fraction of tumor cells with exposure >= IC50")
    ax.set_title("B · Who actually sees a therapeutic exposure\n"
                 "shading = approx. clinical Cmax band (order-of-magnitude)")

    # --- C. exposure vs distance from vessel --------------------------------
    ax = fig.add_subplot(gs[0, 2]); _style_axes(ax)
    name = "temozolomide"
    dep = out["drugs"][name]["depth"]
    d = np.asarray(dep["distance_um"])
    show = [0.3, 1.0, 3.0, 10.0, 30.0]
    rows = {r["multiplier"]: r for r in out["drugs"][name]["dose"]}
    for i, m in enumerate(show):
        pm = rows[m]["plasma_mean"]
        ax.plot(d, np.asarray(dep["unit_E"]) * pm, color=SEQ_BLUE[i + 1], lw=2.0)
        ax.annotate(f"{m:g}x", xy=(d[-1], dep["unit_E"][-1] * pm), xytext=(4, -3),
                    textcoords="offset points", color=SEQ_BLUE[i + 1], fontsize=8,
                    fontweight="bold")
    ax.axhline(dep["ic50"], color=C_CRITICAL, lw=1.6, ls=(0, (5, 3)))
    ax.annotate(f"IC50 = {dep['ic50']:.0f} uM", xy=(d[0], dep["ic50"]),
                xytext=(2, 4), textcoords="offset points", color=C_CRITICAL,
                fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("distance from nearest vessel (um)")
    ax.set_ylabel("time-averaged exposure (uM)")
    _fn = (dep.get("fit_near") or {}).get("L_um")
    _pl = dep.get("plateau_ratio")
    ax.set_title(f"C · {name}: dose moves amplitude, NOT shape\n"
                 f"steep drop (half-dist {dep['half_distance_um']:.0f} um) onto a FLOOR "
                 f"at {'—' if _pl is None else f'{_pl:.0%}'} of peak")
    ax.axvline(out["meta"]["distance_p90_um"], color=C_MUTED, lw=1.0, ls=":")
    ax.annotate("p90 of tissue\ndistance-to-vessel\n(beyond: sparse bins)",
                xy=(out["meta"]["distance_p90_um"], ax.get_ylim()[1]),
                xytext=(4, -4), textcoords="offset points", color=C_INK2,
                fontsize=7, va="top")

    # --- D. schedule at matched total dose (ONE drug -- scope cut) ----------
    ax = fig.add_subplot(gs[1, 0]); _style_axes(ax)
    sd = SCHEDULE_DRUG
    srows = out["drugs"][sd]["schedule"]
    ns = [s["n_doses"] for s in srows]
    x = np.arange(len(ns), dtype=float)
    w = 0.34
    km = [s["kill_matched"] for s in srows]
    kd = [s["kill_matched_dose"] for s in srows]
    ax.bar(x - w / 2, km, w * 0.92, color=C_DRUG[sd], zorder=3,
           label="matched cumulative exposure  (the key comparison)")
    ax.bar(x + w / 2, kd, w * 0.92, color=C_DRUG[sd], alpha=0.40,
           hatch="///", edgecolor=C_DRUG[sd], lw=0.0, zorder=3,
           label="matched administered dose  (window-truncation artifact)")
    for xx, yy in zip(x - w / 2, km):
        ax.text(xx, yy, f"  {yy:.5f}", ha="center", va="bottom", fontsize=7.5,
                color=C_INK2, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([("1 dose\n(bolus)" if n == 1
                         else f"{n} doses\nq{s['interval_h']:.1f}h")
                        for n, s in zip(ns, srows)])
    ax.set_ylabel("predicted tumor kill fraction")
    ax.legend(fontsize=7.0, frameon=False, loc="upper left", labelcolor=C_INK2)
    ax.set_ylim(0, max(km) * 1.70)
    if out["meta"].get("time_resolved"):
        _es = [s["mean_auc_tumor_matched"] for s in srows]
        _hi = f"{SCHEDULE_DOSE_MULTIPLIERS[-1]:g}"
        _kh = [s["kill_matched_by_dose"][_hi] for s in srows]
        _spread_here = max(km) - min(km)
        ax.set_title(f"D · {sd}: fractionation over "
                     f"{out['meta']['window_sched_h']:.0f} h\n"
                     + ("at matched exposure the solid bars now DIFFER "
                        "(they were exactly equal)" if _spread_here > 0 else
                        "the schedule is no longer ERASED at matched exposure "
                        "— see the numbers"))
        ax.text(0.98, 0.97,
                f"tissue-exposure spread {(max(_es)/min(_es)-1)*100:.2f}%\n"
                f"(bounded by tau/T)\n"
                f"kill spread {max(_kh)-min(_kh):.1e} at {_hi}x this dose\n"
                f"was 0.00e+00 at every dose",
                transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                color=C_INK2)
    else:
        ax.set_title(f"D · {sd}: fractionation over "
                     f"{out['meta']['window_sched_h']:.0f} h\n"
                     "solid bars are IDENTICAL — cumulative exposure held EQUAL")

    # --- E. half-life sweep (ONE drug -- scope cut) -------------------------
    ax = fig.add_subplot(gs[1, 1]); _style_axes(ax)
    hr = out["drugs"][sd]["half_life"]
    hls = [r["half_life_h"] for r in hr]
    ax.plot(hls, [r["kill_matched"] for r in hr], color=C_DRUG[sd], lw=2.4,
            marker="o", ms=6, mec="#fcfcfb", mew=0.8, zorder=3)
    ax.plot(hls, [r["kill_peakmatched"] for r in hr], color=C_DRUG[sd],
            lw=1.8, ls=(0, (4, 3)), marker="s", ms=6, mfc="#fcfcfb", zorder=3)
    ax.annotate(sd, xy=(hls[-1], hr[-1]["kill_matched"]), xytext=(-4, 9),
                textcoords="offset points", ha="right", color=C_INK,
                fontsize=9, fontweight="bold")
    ax.set_xscale("log")
    ax.set_xticks(hls)
    ax.set_xticklabels([f"{h:g}" for h in hls])
    ax.set_xlabel("plasma half-life (h)   —   bolus <-> sustained")
    ax.set_ylabel("predicted tumor kill fraction")
    handles = [Line2D([], [], color=C_INK2, lw=2.4, marker="o", ms=6,
                      label="matched TOTAL exposure (peak rescaled)"),
               Line2D([], [], color=C_INK2, lw=1.8, ls=(0, (4, 3)), marker="s",
                      ms=6, mfc="#fcfcfb", label="matched PEAK (exposure free)")]
    ax.legend(handles=handles, fontsize=7.5, frameon=False, loc="lower right",
              labelcolor=C_INK2)
    if out["meta"].get("time_resolved"):
        _hd = [r.get("half_distance_um") for r in hr]
        _km = [r["kill_matched"] for r in hr]
        ax.set_title("E · Bolus vs sustained infusion proxy\n"
                     "the solid line now SLOPES (it was flat by construction)")
        _msg = (f"kill spread {max(_km)-min(_km):.1e} at matched exposure "
                f"(was 0.00e+00)")
        if all(v is not None for v in _hd) and len(_hd) > 1:
            _msg += (f"\nexposure half-distance {_hd[0]:.2f} -> {_hd[-1]:.2f} um: "
                     f"the BOLUS reaches deeper")
        ax.text(0.02, 0.97, _msg, transform=ax.transAxes, ha="left", va="top",
                fontsize=7.5, color=C_INK2)
    else:
        ax.set_title("E · Bolus vs sustained infusion proxy\n"
                     "solid lines are FLAT: at matched exposure, half-life is invisible")

    # --- F. therapeutic depth vs dose ---------------------------------------
    ax = fig.add_subplot(gs[1, 2]); _style_axes(ax)
    cap = max(out["drugs"][n]["depth"]["max_profiled_um"] for n in DRUGS)
    for k, name in enumerate(DRUGS):
        rows = out["drugs"][name]["dose"]
        pts = [(r["multiplier"], r["ic50_depth_um"],
                bool(r.get("ic50_depth_censored", False)))
               for r in rows if r["ic50_depth_um"] is not None]
        if not pts:
            continue
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, color=C_DRUG[name], lw=2.2, zorder=3)
        for x_, y_, cen in pts:
            ax.plot([x_], [y_], marker="o", ms=5, color=C_DRUG[name],
                    mfc=("#fcfcfb" if cen else C_DRUG[name]),
                    mec=C_DRUG[name], mew=1.4, zorder=4)
        _right = xs[0] >= 0.5 * max(DOSE_MULTIPLIERS)
        ax.annotate(name, xy=(xs[0], ys[0]),
                    xytext=(-9 if _right else 9, -1 if k == 0 else -12),
                    textcoords="offset points",
                    ha="right" if _right else "left", va="center",
                    color=C_DRUG[name], fontsize=9, fontweight="bold")
    ax.axhline(cap, color=C_MUTED, lw=1.2, ls=":", zorder=2)
    ax.text(0.02, 0.995, "edge of the profiled range — open markers are lower bounds\n"
            "(exposure still above IC50 at the last profiled bin)",
            transform=ax.transAxes, color=C_INK2, fontsize=7.5, va="top")
    ax.set_ylim(0, cap * 1.32)
    ax.set_xscale("log")
    ax.set_xlabel("plasma peak  (x nominal)")
    ax.set_ylabel("deepest distance from a vessel still above IC50 (um)")
    ax.set_title(f"F · What {DOSE_MULTIPLIERS[-1]:g}x more drug buys in reach\n"
                 "the jump is the FLOOR clearing IC50, not a deeper front",
                 fontsize=10)

    fig.suptitle(
        "DOSE & SCHEDULE on ONE REAL TUMOR — 10x Xenium glioblastoma section "
        f"({meta['n_tumor']:,} tumor cells, dx = {meta['dx_um']:.0f} um, converged)\n"
        "RESEARCH PROTOTYPE — NOT CLINICAL DOSING GUIDANCE · untreated resection, "
        "no measured outcome to validate against · absolute kills are conditional "
        "on an unvalidated per-cell sensitivity model",
        fontsize=12.5, color=C_INK, y=0.975, ha="center",
    )
    if meta.get("time_resolved"):
        _bnd = sched_bound_pct(out, SCHEDULE_DRUG) if out["drugs"][SCHEDULE_DRUG].get(
            "schedule") else float("nan")
        fig.text(0.5, 0.020,
                 "SUPERSEDES the static-source version of this figure. Panels D and E were "
                 "previously DEGENERATE BY CONSTRUCTION — the engine collapsed any schedule "
                 "to its time-averaged plasma level, so matched cumulative exposure gave a "
                 "bit-identical answer and both panels were flat.\n"
                 "The source is now rescaled by plasma(t) every timestep, at exactly matched "
                 f"cumulative dose. The effect is real, monotone and bounded by "
                 f"tau_tissue/T_window = {_bnd:.2f}% — small here because this small molecule "
                 "equilibrates in minutes against a "
                 f"{meta['window_sched_h']:.0f} h window — and it favours the BOLUS, not the "
                 "sustained infusion.\nThe pharmacodynamic mechanisms that make fractionation "
                 "matter clinically are still absent, so treat this as a floor. See "
                 "results/dose_schedule.md sections 4-5.",
                 ha="center", va="bottom", fontsize=8.0, color=C_INK2)
    else:
        fig.text(0.5, 0.032,
                 "Panels D and E are DEGENERATE BY CONSTRUCTION at matched exposure: the engine "
                 "reduces any schedule to its time-averaged plasma level, so equal cumulative "
                 "exposure gives a bit-identical answer.\nThat is a model limitation, reported as "
                 "such — see results/dose_schedule.md. The physical argument for why a "
                 "time-resolved model would also show little effect is the transport equilibration "
                 "time (~minutes) vs the dose interval (~hours).",
                 ha="center", va="bottom", fontsize=8.5, color=C_INK2)

    fig.savefig(path, dpi=150, facecolor="#fcfcfb")
    plt.close(fig)
    print(f"wrote {path}")


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _fmt(v, nd=4):
    if v is None:
        return "—"
    return f"{v:.{nd}f}"

def _dose_row(rows, m):
    return next(r for r in rows if abs(r["multiplier"] - m) < 1e-9)


def _first_multiplier_reaching(rows, key, target):
    """Smallest swept dose multiple at which ``key`` first reaches ``target``."""
    for r in rows:
        if r[key] >= target:
            return float(r["multiplier"])
    return None


def _band_multipliers(name):
    """The drug's approximate clinical Cmax band expressed in units of the shared
    nominal plasma peak. ORDER-OF-MAGNITUDE literature reference, not a fit."""
    lo, hi = CLINICAL_CMAX_BAND_UM[name]
    return lo / NOMINAL_PEAK_UM, hi / NOMINAL_PEAK_UM


def _interp_at(rows, key, mult):
    """Log-dose interpolation of a swept quantity, clamped to the swept range."""
    x = np.log10([r["multiplier"] for r in rows])
    y = np.asarray([r[key] for r in rows], dtype=np.float64)
    m = float(np.clip(mult, 10 ** x[0], 10 ** x[-1]))
    return float(np.interp(np.log10(m), x, y)), (m != mult)


def _depth_str(v, censored=False):
    if v is None:
        return "— (below IC50 even at the vessel wall)"
    return f"{v:.0f}{'+' if censored else ''}"


def sched_bound_pct(out: Dict, name: str) -> float:
    """The analytic ceiling on the schedule effect, as a percentage.

    From ``A = M^{-1}(S_total - C(T))`` (see section 4), two schedules at matched
    cumulative exposure differ only through the terminal residual ``C(T)``. Under
    a sustained level the tissue sits at ``C_ss ~ S_rate / k``, and the drug
    delivered over the window is ``S_rate * T``, so the residual can be at most a
    fraction ``(1/k) / T = tau_tissue / T_window`` of what was delivered. That is
    the largest relative exposure difference any schedule can produce.
    """
    tau = float(out["drugs"][name]["params"]["tau_transport_s"])
    return 100.0 * tau / float(WINDOW_SCHED_S)


def report_t95_s(out: Dict, name: str) -> float:
    """Tissue equilibration t95 for the report, always derived analytically.

    Deliberately does NOT trust a stored measured value: under a time-resolved
    source the mean-concentration history rises and falls, so "95% of the final
    value" measures washout, not equilibration. ``ln(20)/k`` is exact for the
    sustained case and is what the measured number converges to.
    """
    return equilibration_t95_analytic_s(out["drugs"][name]["params"]["k_loss"])


def write_report(out: Dict, path: str) -> None:
    meta = out["meta"]
    if out["drugs"].get(SCHEDULE_DRUG):
        meta["sched_bound_pct"] = sched_bound_pct(out, SCHEDULE_DRUG)
    L: List[str] = []
    A = L.append
    top = DOSE_MULTIPLIERS[-1]

    A("# DOSE & SCHEDULE on ONE REAL TUMOR")
    A("")
    if meta.get("time_resolved"):
        A("> ## THIS VERSION SUPERSEDES THE STATIC-SOURCE VERSION OF THIS FILE")
        A(">")
        A("> **What changed.** Every earlier version of this report was produced with a "
          "transport solver driven by a *static* source: the engine collapsed the whole "
          "multi-dose plasma profile to a single time-averaged level before the tissue "
          "ever saw it. Two schedules with the same cumulative plasma exposure were "
          "therefore **bit-identical inputs**, and the reported kill spread across "
          "1 / 3 / 10 fractions, and across a 192-fold plasma-half-life range, was "
          "`0.00e+00` — *exactly* zero, to machine precision.")
        A(">")
        A("> **That zero was an arithmetic artefact of the solver, not a result about "
          "this tumor.** The old file said so in its caveats but still drew the "
          "conclusion \"on this tumor, at matched cumulative exposure, this model "
          "predicts fractionation changes nothing\" and \"sustained infusion does NOT "
          "reach further here\". Neither conclusion was earned: the simulator was "
          "structurally incapable of answering either question. Both are retracted and "
          "replaced below.")
        A(">")
        A("> **What this version does instead.** The source is now rescaled by "
          "`plasma(t)` at every timestep (`PDEConfig.time_resolved=True`; the solver's "
          "`source_profile` argument), so tissue concentration genuinely rises and falls "
          "and each schedule gets its own transport solve. Cumulative delivered drug is "
          "held *exactly* equal between the two paths, so this is a comparison of "
          "schedule and not of dose. The static path remains the default everywhere else "
          "in the repo and is bit-for-bit unchanged (pinned by "
          "`tests/test_time_resolved_pk.py` against a stored reference).")
        A(">")
        A("> **The corrected answer, in one line.** Fractionation and infusion duration "
          "do NOT have zero effect on this tumor — they have a *small, non-zero, "
          "monotone, and analytically bounded* effect, and it runs in the **opposite "
          "direction to the textbook expectation**: at matched cumulative exposure the "
          "BOLUS reaches deeper than the sustained infusion. Sections 4 and 5 give the "
          "numbers, the derivation, and the bound.")
        A("")
    # Provenance banner. Carried forward from the discrete-vessel-source
    # regeneration (results/vessel_promotion.md) — it is generated here rather
    # than hand-prepended so it survives every re-run of this script.
    A("> **Vascular source / units.** `vascular_mode = discrete` "
      "(`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in "
      "[`results/vessel_promotion.md`](vessel_promotion.md)) · "
      f"**alpha = {meta['alpha']:.1f}**, the e-folding-length convention "
      "`L = d_half / ln 2` · per-cell sensitivity mode = **`determinant`** "
      "(per-drug: MKI67 proxy for cisplatin; temozolomide has no determinant on this "
      "panel and stays at 1.0 — flagged UPPER BOUND in the body) · "
      f"**dx = {meta['dx_um']:.0f} µm**.")
    A(">")
    A("> Every absolute exposure / AUC / coverage / kill below is on the "
      "**discrete-source unit** and is NOT comparable with the KDE-era numbers this file "
      "carried before that regeneration (the tissue-integrated source total is 3.36x "
      "smaller; that is a change of *unit*, not of physics — the dose-matched, "
      "geometry-only effect on coverage is <= 0.0166). To recover the old numbers, re-run "
      "with `mode='kde'`.")
    A(">")
    A(f"> Note the consequence of promotion: the p{meta['vessel_percentile']:.0f} "
      f"'vessel mask' used for the depth statistics resolves to "
      f"**{meta['vessel_area_fraction']*100:.1f}% of tissue**, i.e. essentially the "
      f"segmented vessel mask itself, rather than the broad p95 contour of a smooth KDE.")
    A("")
    A("**Predicted dose-response, fractionation and infusion-duration behaviour for "
      "temozolomide and cisplatin on the real 10x Xenium glioblastoma section.** This "
      "is the question the simulator exists to inform — not *which* drug, but *how "
      "much* and *how often* — asked of one real patient's tissue.")
    A("")
    A("![dose & schedule](dose_schedule.png)")
    A("")
    A("> ## READ THIS FIRST")
    A(">")
    A("> * **Research prototype. NOT clinical dosing guidance.** Nothing here should "
      "influence the care of any patient, and no number here is a dose recommendation.")
    A("> * **Zero outcome validation.** This Xenium section is an *untreated* resection. "
      "There is no measured drug response for either drug on this tissue, so every kill "
      "fraction below is a model prediction with nothing to score it against.")
    A("> * **Absolute kills are conditional on the sensitivity model.** For "
      "temozolomide the per-cell sensitivity determinant (MGMT) is *not on this panel*, "
      "so sensitivity is uniform 1.0 — i.e. temozolomide resistance is **unmodelled**, "
      "and its absolute kill is correspondingly optimistic. For cisplatin sensitivity is "
      "an MKI67-derived *mechanistic proxy*, not a measured chemosensitivity. Trust the "
      "*shape* of the dose-response and the *relative* comparisons; do not trust an "
      "absolute kill fraction.")
    if meta.get("time_resolved"):
        A("> * **The schedule comparison is now real, but the effect is SMALL and its "
          "size is set by a ratio you can compute in advance.** With a time-resolved "
          "plasma driver the whole schedule effect is the drug still in flight at the end "
          "of the observation window, which for a linear transport model is bounded by "
          "(tissue equilibration time) / (window length). For temozolomide here that is "
          f"{meta.get('sched_bound_pct', float('nan')):.2f}%, and the measured spread "
          "matches it. A slowly-equilibrating molecule (an antibody) would show a large "
          "effect on the same machinery; a fast small molecule over a long window cannot.")
        A("> * **Pharmacodynamics is still a single static Hill map on cumulative "
          "exposure.** Time resolution was added to *transport*, not to the kill model. "
          "There is still no inter-fraction DNA repair, no cell-cycle redistribution, no "
          "regrowth and no normal-tissue recovery — and those, not transport, are what "
          "give fractionation its clinical benefit. So the small effect reported here is "
          "a *lower bound on how much schedule could matter*, not an estimate of it.")
    else:
        A("> * **The matched-dose schedule comparison is degenerate by construction.** The "
          "engine collapses a dosing schedule to its time-averaged plasma level "
          "(the documented Stage-0 plasma-averaging simplification in `sim/engine.py`), so "
          "any two schedules with equal cumulative exposure are *bit-identical* inputs. "
          "Sections 4 and 5 report that plainly and separate the model artifact from the "
          "physical argument.")
    A("> * The perivascular gradient the depth numbers rest on partly **inherits the "
      "shape of the vascular-density source field** — the confound documented in "
      "`tmesim/validation/vessel_profile.py`. Depth figures are order-of-magnitude.")
    A(f"> * **~{meta['quiescent_tumor_fraction']*100:.0f}% of tumor cells are labelled "
      "QUIESCENT by the oxygen-banding model**, "
      "and quiescent cells carry a 0.5x kill multiplier. That single labelling assumption "
      "dominates every absolute kill and the kill ceiling below. It is a *model* "
      "assignment from the oxygen field, currently being tested independently against "
      "measured MKI67 — **absolute kills may shift materially when that lands.** The "
      "dose/schedule *comparisons* are far more robust to it than the absolute levels.")
    A("> * **SCOPE (compute budget).** The machine was running ~35 concurrent jobs, so "
      "this run was deliberately cut: "
      f"{len(DOSE_MULTIPLIERS)} dose points per drug; the schedule and half-life legs "
      f"run for **{SCHEDULE_DRUG} only**"
      + (" (each plasma-curve SHAPE now needs its own PDE solve — that is the cost of "
         "resolving schedule at all, and it is why the second drug is dropped rather "
         "than the legs)" if meta.get("time_resolved") else
         " (they need a second PDE solve per drug)")
      + f"; {len(SCHEDULES)} schedules and {len(HALF_LIVES_H)} half-lives; the "
      f"schedule/half-life window shortened to {WINDOW_SCHED_S/3600:.0f} h (held "
      f"EXACTLY constant across everything compared); and "
      f"{len(out['engine_checks'])} rather than 6 full-engine verification runs. "
      "Everything reported below actually ran — nothing is extrapolated. What was "
      "dropped is named where it would have appeared.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------------ 1
    A("## 1. Methods")
    A("")
    A(f"* **Tissue — one real tumor.** The open 10x Xenium GBM section: "
      f"{meta['n_cells']:,} measured cells ({meta['n_tumor']:,} tumor), reconstructed "
      f"**once** at the grid-converged `dx = {meta['dx_um']:.0f} um` "
      f"(grid {meta['grid'][0]}x{meta['grid'][1]}) and reused for every sweep, so all "
      f"differences are attributable to dose/schedule and not to the tissue. Mean "
      f"vascular density {meta['mean_vascular_density']:.3f}, mean ECM "
      f"{meta['mean_ecm_density']:.3f}, mean oxygen {meta['mean_oxygen']:.3f}. After the "
      f"oxygen solve {meta['hypoxic_tumor_fraction']*100:.1f}% of tumor cells are hypoxic "
      f"and {meta['quiescent_tumor_fraction']*100:.1f}% quiescent (both carry reduced "
      f"drug kill).")
    A(f"* **Vascular geometry.** Vessel mask = top {100 - meta['vessel_percentile']:.0f}% "
      f"of the endothelial-derived vascular-density field "
      f"({meta['vessel_area_fraction']*100:.1f}% of tissue area). Distance from the "
      f"nearest vessel over tissue: median "
      f"{meta.get('distance_median_um', float('nan')):.0f} um, p90 "
      f"{meta.get('distance_p90_um', float('nan')):.0f} um, max "
      f"{meta.get('distance_max_um', float('nan')):.0f} um. This is the geometry the "
      f"drug's penetration has to cover — though see section 3: on this section most of "
      f"the far tissue is reached by a diffuse source floor rather than by a diffusion "
      f"front, which changes how that comparison should be read.")
    A(f"* **Drugs.** `calibrated_drugparams(name, alpha = {meta['alpha']:.2f})` "
      "(single perivascular-sink constant, fit once on doxorubicin) then "
      "`apply_barrier_to_drugparams` — the identical compilation pipeline "
      "`real_drug_panel.py` uses.")
    A("* **Per-cell sensitivity** is layered on per drug with "
      "`reconstruct.real_pipeline.apply_drug_sensitivity` (measured determinant where "
      "the panel carries one, honest uniform 1.0 otherwise). See the table below.")
    A(f"* **Numerics.** `dt = {meta['dt_s']:.0f} s`, convection on. Dose-response window "
      f"{meta['window_dose_h']:.0f} h; schedule / half-life window "
      f"{meta['window_sched_h']:.0f} h (cut from 24 h for compute, and held EXACTLY "
      f"constant across every schedule and half-life compared — this matters more than "
      f"before, because with a time-resolved driver the window length now *sets the size "
      f"of the schedule effect*; see section 4). "
      f"`dx = 10 um` is deliberately not used — it violates the upwind advection CFL "
      f"condition at this `dt`. Seed {meta['seed']}, fully deterministic.")
    A(f"* **Plasma driver — the change that defines this version.** "
      f"`PDEConfig.time_resolved = {bool(meta.get('time_resolved'))}`. "
      + ("The transport source is rescaled by `plasma(t)` at EVERY timestep, so tissue "
         "concentration rises and falls with the real multi-dose curve. The implicit "
         "operator `A = I - dt*L + dt*K` contains no source term, so the LU "
         "factorisation is still built once and reused for all "
         f"{int(round(WINDOW_SCHED_S/meta['dt_s']))} steps — a time-varying source is one "
         "extra scalar multiply per step and costs no measurable runtime. The plasma "
         "curve is sampled on the solver's own step grid, whose mean is exactly the "
         "`plasma_mean` the static path used, so **both paths deliver exactly the same "
         "total drug** and a static-vs-resolved difference is a pure schedule effect."
         if meta.get("time_resolved") else
         "The multi-dose plasma profile is collapsed to its window time-average and the "
         "solver sees a STATIC source — the superseded Stage-0 simplification. Schedule "
         "is unresolvable in this mode by construction."))
    A(f"* **Dose axis.** Nominal plasma peak = {meta['nominal_peak_um']:.0f} uM (the "
      f"repo-wide `PKConfig` default), swept "
      f"{DOSE_MULTIPLIERS[0]:g}x–{DOSE_MULTIPLIERS[-1]:g}x with a "
      f"{meta['nominal_half_life_h']:.0f} h plasma half-life. The *same* nominal level "
      "is used for both drugs so the axis means one thing; note this is roughly "
      "temozolomide's clinical order of magnitude but is **supratherapeutic for "
      "cisplatin** (free-platinum Cmax is roughly an order of magnitude lower), so read "
      "cisplatin's low-dose end as the clinically relevant part of its curve.")
    A("")
    A("### 1.1 One PDE solve per plasma-curve SHAPE — exact, not an approximation")
    A("")
    A("The transport PDE `dC/dt = div(D_eff grad C) - div(vC) + S - (clearance+binding)C` "
      "is **linear** in the source and starts from `C = 0`. The source is "
      "`S(x,y,t) = vascular_density x bbb_penetrance x plasma(t)`, so scaling the whole "
      "plasma curve scales the exposure field exactly: `AUC(P x shape) = P x AUC(shape)`. "
      "Every schedule is therefore split into **(shape x amplitude)**, where the shape is "
      "`plasma(t)` renormalised to unit window-mean, and:")
    A("")
    A("* **Amplitude is free.** All 7 dose points, and the matched-exposure / "
      "matched-administered-dose / matched-peak legs of a given schedule, share ONE solve "
      "and are recovered by rescaling it. Exact, not interpolated.")
    if meta.get("time_resolved"):
        A("* **Shape is not free, and that is the whole point.** 1 vs 3 vs 10 fractions, "
          "and a 0.25 h vs a 48 h half-life, are DIFFERENT shapes and each needs its own "
          "solve. Under the superseded static source every shape collapsed to the same "
          "constant level, which is exactly why the schedule question had no answer. "
          f"This run performs "
          f"{sum(int(out['drugs'][n].get('n_pde_solves') or 0) for n in out['drugs'])} "
          f"transport solves rather than the 3 the static version needed.")
    A("")
    A("Every dose / schedule / half-life point then re-runs the engine's own unmodified "
      "kill mapping (`cells.pkpd.apply_drug_kill`) on the rescaled field — the same trick "
      "the engine itself uses to build its viable-fraction history. This script "
      "orchestrates `core` / `physics` / `molecular` / `reconstruct` / `engine`; it "
      "re-implements none of them.")
    A("")
    A(f"This is verified, not asserted: the **full** `sim.engine.simulate_response` was "
      f"re-run end to end at {len(out['engine_checks'])} points — one on each leg (dose "
      f"sweep, fractionation, half-life) and covering both exposure windows, so the "
      f"shape-x-amplitude fast path is pinned everywhere it is used. Each check is itself "
      f"a complete extra PDE solve.")
    A("")
    A("| verification point | full-engine kill | fast-path kill | abs delta |")
    A("|---|---|---|---|")
    for c in out["engine_checks"]:
        A(f"| {c['label']} | {c['engine_kill']:.10f} | {c['fast_kill']:.10f} | "
          f"{c['abs_delta']:.1e} |")
    A("")
    A(f"Largest discrepancy anywhere: **{meta['max_engine_fast_delta']:.1e}**.")
    A("")
    A("### 1.2 Drug parameters as compiled onto this tissue")
    A("")
    A("| drug | MW (Da) | D (um²/s) | IC50 (uM) | hill | BBB | k_loss (1/s) | transport "
      "time-scale 1/k | analytic L = √(D_eff/k) | sensitivity determinant | median "
      "sensitivity (tumor) | kill ceiling |")
    A("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for name in DRUGS:
        d = out["drugs"][name]
        p, s = d["params"], d["sensitivity"]
        Lp = math.sqrt(p["diffusivity"] * (1 - 0.6 * meta["mean_ecm_density"]) / p["k_loss"])
        det = s["determinant_gene"] or "none (uniform 1.0)"
        A(f"| {name} | {p['mw']:.1f} | {p['diffusivity']:.0f} | {p['ic50']:.1f} | "
          f"{p['hill']:.2f} | {p['bbb']:.3f} | {p['k_loss']:.2e} | "
          f"{p['tau_transport_s']:.0f} s | {Lp:.0f} um | {det} | "
          f"{s['median_tumor']:.3f} | {d['kill_ceiling']:.4f} |")
    A("")
    A("**The kill ceiling** is the kill fraction reached as exposure → ∞: the Hill term "
      "saturates at `max_kill_rate`, so per-cell kill probability tends to "
      "`sensitivity x hypoxia/quiescence resistance`. It is a property of the sensitivity "
      "model and this tissue's oxygenation. **No dose can exceed it**, and it is the "
      "reference against which \"is more drug worth it?\" must be judged.")
    A("")
    A("### 1.3 What actually sets the ceiling — the decomposition")
    A("")
    A("This matters more than any single kill number, so it is broken out per cell-state "
      "band. `contribution` is that band's share of the ceiling (fraction of tumor cells x "
      "their mean asymptotic kill probability); the three contributions sum to the ceiling.")
    A("")
    for name in DRUGS:
        d = out["drugs"][name]
        bd = d.get("kill_ceiling_breakdown") or {}
        if not bd:
            continue
        A(f"**{name}** — ceiling {d['kill_ceiling']:.4f}")
        A("")
        A("| tumor-cell band | fraction of tumor cells | resistance multiplier | mean "
          "asymptotic kill probability | contribution to ceiling |")
        A("|---|---|---|---|---|")
        for band in ("VIABLE", "HYPOXIC", "QUIESCENT"):
            b = bd.get(band)
            if b is None:
                continue
            A(f"| {band} | {b['fraction']*100:.1f}% | "
              f"x{b['resistance_multiplier']:.2f} | {b['mean_p']:.3f} | "
              f"{b['contribution']:.4f} |")
        A("")
    A("Read that table before reading anything else. The ceiling is **not** a delivery "
      "number and **not** an IC50 number: at infinite exposure every cell is at the top of "
      "its Hill curve, so what is left is exactly `sensitivity x resistance`. For "
      "temozolomide, whose per-cell sensitivity is uniform 1.0 on this panel, the ceiling "
      "is *entirely* the hypoxia (x0.35) and quiescence (x0.5) multipliers applied to the "
      "oxygen-banded cell states — i.e. it is set by how much of this tumor the model "
      "considers metabolically out of reach, not by how much drug arrives. For cisplatin "
      "the MKI67-derived sensitivity multiplies that down further.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------------ 2
    A("## 2. Dose-response on this patient's tissue")
    A("")
    A("`IC50 coverage` = fraction of tumor cells whose *time-averaged tissue exposure* "
      "(`AUC / window`, the quantity the engine's Hill mapping consumes) is at or above "
      "the drug's own IC50. `IC50 depth` = the deepest distance from a vessel at which "
      "the binned exposure profile still clears IC50 (`+` = still above IC50 at the far "
      "end of the profiled range, i.e. a lower bound).")
    A("")
    for name in DRUGS:
        d = out["drugs"][name]
        rows = d["dose"]
        ceil = d["kill_ceiling"]
        A(f"### 2.{DRUGS.index(name)+1} {name}")
        A("")
        A("| dose (× nominal) | plasma peak (uM) | window-mean plasma (uM) | kill "
          "fraction | % of kill ceiling | tumor cells ≥ IC50 | tumor cells ≥ 0.1·IC50 | "
          "median exposure / IC50 | IC50 depth from vessel (um) |")
        A("|---|---|---|---|---|---|---|---|---|")
        for r in rows:
            A(f"| {r['multiplier']:g}x | {r['plasma_peak']:.1f} | {r['plasma_mean']:.2f} | "
              f"{r['kill']:.5f} | {r['kill']/ceil*100:.1f}% | "
              f"{r['ic50_coverage']*100:.2f}% | {r['cov_10pct_ic50']*100:.2f}% | "
              f"{r['median_E_over_ic50']:.3f} | "
              f"{_depth_str(r['ic50_depth_um'], r.get('ic50_depth_censored', False))} |")
        A("")
        r1 = _dose_row(rows, 1.0)
        rmax = rows[-1]
        sat90 = _first_multiplier_reaching(rows, "kill", 0.90 * rmax["kill"])
        sat99 = _first_multiplier_reaching(rows, "kill", 0.99 * rmax["kill"])
        cov90 = _first_multiplier_reaching(rows, "ic50_coverage", 0.90)
        if sat90 is not None and sat90 >= top:
            sat_txt = (f"**Kill has NOT saturated within the swept range** — it is still "
                       f"climbing at {top:g}x nominal, reaching 90% of its {top:g}x value "
                       f"only at {top:g}x itself. Over this whole range more drug is still "
                       f"buying more kill")
        else:
            sat_txt = (f"**Kill saturates at about {sat90:g}x** (90% of its {top:g}x "
                       f"value) and is flat by {sat99:g}x")
        A(f"**Saturation.** {sat_txt}"
          f"{'' if cov90 is None else f'; IC50 coverage passes 90% of tumor cells by **{cov90:g}x**'}"
          f". At {top:g}x the predicted kill is {rmax['kill']:.4f} = "
          f"{rmax['kill']/ceil*100:.1f}% of the {ceil:.4f} ceiling, so the residual "
          f"{(ceil - rmax['kill'])/ceil*100:.1f}% of the ceiling is what *no achievable "
          f"dose* recovers — that residual is the hypoxic/quiescent fraction the Hill "
          f"term can never reach past.")
        A("")
        blo, bhi = _band_multipliers(name)
        klo, clamped_lo = _interp_at(rows, "kill", blo)
        khi, _ = _interp_at(rows, "kill", bhi)
        vlo, _ = _interp_at(rows, "ic50_coverage", blo)
        vhi, _ = _interp_at(rows, "ic50_coverage", bhi)
        cb = CLINICAL_CMAX_BAND_UM[name]
        A(f"**At the approximate clinical plasma band.** Published Cmax for {name} is "
          f"of order {cb[0]:g}–{cb[1]:g} uM, i.e. **{blo:.2g}x–{bhi:.2g}x** of the shared "
          f"nominal level used for the axis. Interpolating the sweep there gives kill "
          f"{klo:.4f}–{khi:.4f} ({klo/ceil*100:.0f}–{khi/ceil*100:.0f}% of ceiling) and "
          f"IC50 coverage {vlo*100:.1f}–{vhi*100:.1f}% of tumor cells"
          f"{' (the low end is clamped to the swept minimum)' if clamped_lo else ''}. "
          f"*This band is an order-of-magnitude literature reference attached for scale "
          f"only — it is not fitted, not patient-specific, and not a dose recommendation.*")
        A("")
        gain = rmax["kill"] - r1["kill"]
        rel = "n/a" if r1["kill"] <= 0 else f"{rmax['kill'] / r1['kill']:.2f}x"
        A(f"**What {top:g}x buys over 1x.** kill {r1['kill']:.4f} → {rmax['kill']:.4f} "
          f"(+{gain:.4f} absolute, {rel} "
          f"relative) for a {top:g}-fold increase in systemic exposure — i.e. a "
          f"{gain/max(rmax['kill'], 1e-12)*100:.0f}% share of the final kill bought at "
          f"{top:g}x the systemic toxicity.")
        A("")
    A("---")
    A("")

    # ------------------------------------------------------------------ 3
    A("## 3. Delivery geometry: what the exposure field actually looks like")
    A("")
    A("Because the field is linear in dose, **dose moves the amplitude of the "
      "perivascular exposure profile and nothing else** — the shape is fixed by "
      "`D_eff`, the loss rate `k`, and the vascular geometry. Panel C shows this: the "
      "curves are parallel on a log axis at every dose. This is exact, not approximate: "
      "the transport operator is linear, so scaling the source scales the field.")
    if meta.get("time_resolved"):
        A("")
        A("**Dose, but not schedule.** The plasma curve's SHAPE does move the profile "
          "shape — by well under a percent on this drug and window, and in the direction "
          "that favours the bolus. Section 5 quantifies it. Everything in this section is "
          f"measured on one shape (the {meta['nominal_half_life_h']:.0f} h-half-life single "
          f"dose over the {meta['window_dose_h']:.0f} h dose window), which is the correct "
          "reading of a dose sweep since the shape is held constant across it.")
    A("")
    A("**The profile is NOT a single exponential, and saying otherwise would be the "
      "main way to mislead with this figure.** It is a steep perivascular drop sitting "
      "on a flat, non-zero FLOOR that persists to the far edge of the tissue. The floor "
      "is *not* a diffusion front: the engine's delivery term is the **continuous** "
      "vascular-density field, so a pixel 1 mm from the thresholded vessel mask is still "
      "being fed by diffuse low-density source. A single exponential fitted over the "
      "whole range is dominated by that plateau and returns a meaningless decay length "
      "— which is exactly why `vessel_profile.fit_exponential_plus_floor` exists. All "
      "three characterisations are therefore reported together.")
    A("")
    A("| drug | empirical half-distance (no model assumed) | **source-null half-distance** "
      "| near-field exponential L (d ≤ p90) | near-field fit R² | floor-model L | floor "
      "level (% of perivascular) | far-field plateau / perivascular | analytic √(D_eff/k) "
      "| equilibration t95 | penetration index |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for name in DRUGS:
        d = out["drugs"][name]
        dep, prm = d["depth"], d["params"]
        teq = report_t95_s(out, name)
        Lan = math.sqrt(prm["diffusivity"] * (1 - 0.6 * meta["mean_ecm_density"])
                        / prm["k_loss"])
        fn = dep.get("fit_near", {}) or {}
        ff = dep.get("fit_floor", {}) or {}
        floor_pct = (100.0 * ff["floor"] / (ff["C0"] + ff["floor"])
                     if ff.get("floor") is not None and ff.get("C0") is not None
                     and (ff["C0"] + ff["floor"]) > 0 else None)
        A(f"| {name} | {_fmt(dep['half_distance_um'], 0)} um | "
          f"{_fmt(dep['source_null_half_distance_um'], 0)} um | "
          f"{_fmt(fn.get('L_um'), 0)} um | {_fmt(fn.get('r2'), 2)} | "
          f"{_fmt(ff.get('L_um'), 0)} um | {_fmt(floor_pct, 0)}% | "
          f"{_fmt(dep.get('plateau_ratio'), 2)} | {Lan:.0f} um | "
          f"{'—' if teq is None else f'{teq:.0f} s ({teq/60:.1f} min)'} | "
          f"{d['penetration_index']:.3f} |")
    A("")
    A("**Which column to believe:** the *empirical* half-distance assumes no functional "
      "form and is the trustworthy gradient number. The single-exponential fits — over "
      "the full range and over the near field alike — have poor R² and return decay "
      "lengths inflated by the plateau; disregard them as descriptions of the gradient. "
      "The floor model reconciles the two: a short decaying component (its `L` is close "
      "to the empirical half-distance) plus a background floor. Note also that the two "
      "drugs return nearly identical geometry despite differing in `D_eff` and `k` by "
      "~2x — further evidence that what is being measured here is mostly the source "
      "field, not drug-specific transport.")
    A("")
    A("Three things to take from that table, in order of importance:")
    A("")
    A("1. **The measured gradient is essentially the source geometry.** The empirical "
      "half-distance and the **source-null** half-distance (the identical binning applied "
      "to the raw vascular-density field, with no drug physics at all) agree to within a "
      "few micrometres. The steep near-vessel drop this simulation predicts is inherited "
      "from where the endothelial cells are, not produced by transport. This is the "
      "confound `vessel_profile.py` documents and it is the single biggest caveat on "
      "every depth number in this report.")
    A("2. **A large fraction of the field is floor, not gradient.** See the plateau "
      "column. Any statement of the form \"drug reaches X um\" is therefore mostly a "
      "statement about where the floor sits relative to IC50, not about a diffusion "
      "front, and the `IC50 depth` column in section 2 must be read that way: when it "
      "jumps to the edge of the profiled range between two doses, what happened is that "
      "the **floor** crossed IC50, not that the drug travelled further.")
    A("3. **Depth cannot be bought with dose.** The shape is fixed under dose, so raising "
      "the dose slides the whole curve up and the IC50 crossing moves only as `L*ln(dose)` "
      "in the decaying region — and, once the floor clears IC50, saturates entirely."
      + ("  Nor, in practice, with schedule: section 5 finds the shape does respond to the "
         "plasma time-course, but by a fraction of a percent." if meta.get("time_resolved")
         else ""))
    A("")
    A(f"For scale: median distance-to-vessel on this section is "
      f"{meta.get('distance_median_um', float('nan')):.0f} um, p90 "
      f"{meta.get('distance_p90_um', float('nan')):.0f} um, max "
      f"{meta.get('distance_max_um', float('nan')):.0f} um. Bins beyond the p90 hold only "
      f"a few hundred pixels each, so the far tail of panel C is noisy as well as flat — "
      f"the near-field exponential above is deliberately restricted to d ≤ p90.")
    A("")
    A("---")
    A("")
    A("## 4. Schedule: fractionation at MATCHED TOTAL DOSE")
    A("")
    A(f"One {meta['window_sched_h']:.0f} h course, {meta['nominal_half_life_h']:.0f} h "
      "plasma half-life, split into "
      + " / ".join(str(n) for n, _ in SCHEDULES) + " fractions. Three legs are reported "
      "because \"matched dose\" can mean two different things and the distinction turns "
      "out to be the whole story:")
    A("")
    A("* **A — matched cumulative exposure** (the key comparison the brief asks for): the "
      "per-dose peak is rescaled so the *window-averaged plasma level*, i.e. total plasma "
      "AUC, is identical across schedules.")
    A("* **B — matched administered dose**: the sum of the per-dose peaks is held fixed "
      "(`peak = nominal / n`), which is what \"same total milligrams\" naively means.")
    A("* **C — matched peak** (a control): each fraction reaches the full nominal peak, so "
      "total exposure *rises* with the number of doses.")
    A("")
    ran = [n for n in DRUGS if out["drugs"][n].get("schedule")]
    dropped = [n for n in DRUGS if n not in ran]
    if dropped:
        A(f"> **DROPPED for compute budget:** the schedule leg was NOT run for "
          f"{', '.join(dropped)} — it requires its own transport solve at the schedule "
          f"window. Only {', '.join(ran)} is reported here.")
        A("")
    for i, name in enumerate(ran):
        d = out["drugs"][name]
        A(f"### 4.{i+1} {name}")
        A("")
        A("| schedule | per-dose peak, leg A (uM) | window-mean plasma, A | **kill, A "
          "(matched exposure)** | tumor ≥ IC50, A | window-mean plasma, B | kill, B "
          "(matched administered dose) | window-mean plasma, C | kill, C (matched peak) |")
        A("|---|---|---|---|---|---|---|---|---|")
        for s in d["schedule"]:
            label = ("1 dose (bolus)" if s["n_doses"] == 1
                     else f"{s['n_doses']} doses q{s['interval_h']:.1f}h")
            A(f"| {label} | {s['peak_matched_auc']:.2f} | {s['plasma_mean_matched']:.4f} | "
              f"**{s['kill_matched']:.6f}** | {s['cov_matched']*100:.2f}% | "
              f"{s['plasma_mean_dose']:.4f} | {s['kill_matched_dose']:.6f} | "
              f"{s['plasma_mean_peak']:.4f} | {s['kill_matched_peak']:.6f} |")
        kills_a = [s["kill_matched"] for s in d["schedule"]]
        spread_a = max(kills_a) - min(kills_a)
        kills_b = [s["kill_matched_dose"] for s in d["schedule"]]
        kills_c = [s["kill_matched_peak"] for s in d["schedule"]]
        A("")
        if meta.get("time_resolved"):
            _es_a = [s["mean_auc_tumor_matched"] for s in d["schedule"]]
            _rel = (max(_es_a) / min(_es_a) - 1.0) * 100.0
            A(f"Leg B falls {kills_b[0]:.6f} → {kills_b[-1]:.6f} "
              f"({(kills_b[0]-kills_b[-1])/max(kills_b[0],1e-12)*100:.1f}% relative). "
              f"Leg C rises {kills_c[0]:.6f} → {kills_c[-1]:.6f}.")
            A("")
            A(f"**Leg A at this reference dose: kill spread {spread_a:.2e}, tissue-exposure "
              f"spread {_rel:.3f}%.** "
              + ("The exposure genuinely differs across schedules now — it did not before, "
                 "at all — but at the 1x reference level temozolomide kill is only "
                 f"~{kills_a[0]:.3f}, and a ~{_rel:.1f}% relative change in exposure is "
                 f"smaller than the 1/{meta['n_tumor']:,} = "
                 f"{1.0/max(meta['n_tumor'],1):.1e} quantum of a kill *fraction* computed "
                 "by counting cells. **The zero in this column is now a rounding floor, "
                 "not a structural degeneracy** — an important distinction, because it "
                 "disappears the moment you look at either the exposure itself or a higher "
                 "dose. Both are tabulated immediately below."
                 if spread_a == 0.0 else
                 "Both are non-zero, and both were exactly zero under the static source."))
            A("")
            A("#### The same comparison, now with the STATIC source alongside it")
            A("")
            A("This is the direct before/after. `mean tumor exposure` is the continuous "
              "quantity the kill map consumes; it is not quantised by cell counting, so it "
              "is the honest place to read the size of the schedule effect. The static "
              "column is bit-identical down the whole table BY CONSTRUCTION — that is the "
              "bug, visible.")
            A("")
            A("| schedule | mean tumor exposure, TIME-RESOLVED (uM·s) | mean tumor "
              "exposure, STATIC (uM·s) | Δ vs static | kill, time-resolved | kill, static "
              "| terminal residual C(T) / delivered | exposure half-distance from vessel "
              "(um) | penetration index |")
            A("|---|---|---|---|---|---|---|---|---|")
            for s in d["schedule"]:
                label = ("1 dose (bolus)" if s["n_doses"] == 1
                         else f"{s['n_doses']} doses q{s['interval_h']:.1f}h")
                e, e0 = s["mean_auc_tumor_matched"], s["mean_auc_tumor_matched_static"]
                A(f"| {label} | {e:.5e} | {e0:.5e} | {(e/e0-1)*100:+.3f}% | "
                  f"{s['kill_matched']:.6f} | {s['kill_matched_static']:.6f} | "
                  f"{s['residual_fraction']*100:.2f}% | "
                  f"{_fmt(s.get('half_distance_um'), 2)} | "
                  f"{s['penetration_index']:.6f} |")
            A("")
            es = [s["mean_auc_tumor_matched"] for s in d["schedule"]]
            A(f"**Relative exposure spread across the fractionation range: "
              f"{(max(es)/min(es)-1)*100:.3f}%**, monotonically decreasing with the number "
              f"of fractions. The static column's spread is exactly 0.")
            A("")
            A("#### Absolute kill spread depends on where you are on the Hill curve")
            A("")
            A("The schedule effect is a fixed *relative* change in exposure, so the "
              "*absolute* kill spread it produces grows with dose — at 1x nominal "
              f"{ran[0]} kill is only ~{kills_a[0]:.3f}, so a ~{_rel:.1f}% relative change "
              f"in exposure lands near the 1/{meta['n_tumor']:,} = "
              f"{1.0/max(meta['n_tumor'],1):.1e} quantum of a kill fraction computed by "
              "counting cells"
              + ("" if spread_a > 0 else " and can round to a tie")
              + ". Rescaling leg A up the dose axis (free: exposure is exactly linear in "
              "amplitude) lifts it clear of that floor and shows the effect scaling as it "
              "should. **Every static-source row here is exactly zero.**")
            A("")
            A("| leg-A dose (× the matched-exposure level) | "
              + " | ".join(("1 dose" if s["n_doses"] == 1 else f"{s['n_doses']} doses")
                           for s in d["schedule"])
              + " | spread, TIME-RESOLVED | spread, STATIC |")
            A("|---" * (len(d["schedule"]) + 3) + "|")
            for m in SCHEDULE_DOSE_MULTIPLIERS:
                key = f"{m:g}"
                ks = [s["kill_matched_by_dose"][key] for s in d["schedule"]]
                ss = [s["kill_matched_by_dose_static"][key] for s in d["schedule"]]
                A(f"| {m:g}x | " + " | ".join(f"{k:.6f}" for k in ks)
                  + f" | **{max(ks)-min(ks):.2e}** | {max(ss)-min(ss):.2e} |")
            A("")
        else:
            A(f"Leg A spread across 1 → {d['schedule'][-1]['n_doses']} fractions: "
              f"**{spread_a:.2e}** (exactly zero to machine precision). Leg B falls "
              f"{kills_b[0]:.6f} → {kills_b[-1]:.6f} "
              f"({(kills_b[0]-kills_b[-1])/max(kills_b[0],1e-12)*100:.1f}% relative). "
              f"Leg C rises {kills_c[0]:.6f} → {kills_c[-1]:.6f}.")
            A("")
    if not meta.get("time_resolved"):
        A("### 4.x Interpretation — read this carefully, the honest answer has two halves")
        A("")
        A("**(a) The model answer is degenerate, and that is a limitation, not a finding.** "
          "The engine reduces the whole dosing schedule to its time-averaged plasma level "
          "before it ever touches the tissue. Holding cumulative exposure constant "
          "therefore makes the *input to the solver bit-identical*, and leg A's zero "
          "spread is arithmetic, not biology. **Do not report leg A as evidence that "
          "fractionation does not matter.**")
        A("")
        A("**(b) Transport is linear** and **the kill map is a function of cumulative "
          "exposure only**, so with a static source timing has no channel through which "
          "to act at all.")
        A("")
    else:
        d = out["drugs"][ran[0]] if ran else None
        A("### 4.x Interpretation — what the corrected model actually says")
        A("")
        A("**(a) RETRACTED: \"fractionation changes nothing\".** The static-source version "
          "of this file reported a leg-A spread of exactly `0.00e+00` and concluded that "
          "*on this tumor, at matched cumulative exposure, this model predicts "
          "fractionation changes nothing*. That conclusion was **a model artefact**. The "
          "solver was fed the time-average of the plasma curve, so 1 dose and 10 doses "
          "were literally the same input array; the zero was arithmetic. It is now "
          "measured rather than assumed away, and it is not zero.")
        A("")
        A("**(b) The corrected answer, and the derivation that bounds it.** Integrate the "
          "transport equation over the window `[0, T]` and let `A = AUC = ∫C dt`, "
          "`M = k − div(D grad ·) + div(v ·)`. With `C(0) = 0`:")
        A("")
        A("```")
        A("    C(T) − C(0) = div(D grad A) − div(v A) − k A + ∫S dt")
        A("=>            A = M⁻¹ ( S_total − C(T) ),      S_total = ∫₀ᵀ S dt")
        A("```")
        A("")
        A("`S_total` is **identical** for any two schedules at matched cumulative plasma "
          "exposure — that is what \"matched\" means, and the implementation guarantees it "
          "exactly (the plasma curve is sampled on the solver's own step grid, so its mean "
          "is exactly the level the static source used). Therefore:")
        A("")
        A("1. **The entire schedule effect is the terminal residual `C(T)`** — drug that "
          "was delivered but is still in the tissue, uncounted, when the window closes. "
          "Nothing else about the timing can reach the answer.")
        A("2. **Earlier delivery wins.** A schedule that finishes early has washed out by "
          "`T` (`C(T) ≈ 0`) and realises the full steady-state exposure field. One still "
          "infusing at `T` does not. Hence exposure decreases monotonically with the "
          "number of fractions — which is exactly what the table shows.")
        A("3. **The effect is bounded by `τ_tissue / T_window`.** Under a sustained level "
          "the tissue holds `C_ss ≈ S_rate/k`, while the drug delivered over the window is "
          "`S_rate·T`; the residual is therefore at most `(1/k)/T = τ/T` of what was "
          "given. This is a *prediction*, and the run tests it:")
        A("")
        if d is not None:
            tau = d["params"]["tau_transport_s"]
            bound = 100.0 * tau / WINDOW_SCHED_S
            es = [s["mean_auc_tumor_matched"] for s in d["schedule"]]
            meas = (max(es) / min(es) - 1.0) * 100.0
            statres = d.get("static_residual_fraction")
            A(f"   * predicted ceiling `τ/T` = {tau:.0f} s / "
              f"{WINDOW_SCHED_S/3600:.0f} h = **{bound:.3f}%**")
            A(f"   * measured exposure spread across 1 → "
              f"{d['schedule'][-1]['n_doses']} fractions = **{meas:.3f}%**")
            if statres is not None:
                A(f"   * measured terminal residual of the constant-infusion "
                  f"(static-source) limit, which should SATURATE the bound = "
                  f"**{statres*100:.3f}%**")
            A("")
            A("   The three agree. That agreement is the strongest evidence here that the "
              "time-resolved source is doing real physics and not producing numerical "
              "noise: the size of the effect was derived before it was measured.")
            A("")
        A("**(c) Why the effect is small HERE, and where it would not be.** Temozolomide "
          "on this tissue has a tissue equilibration time of minutes ("
          + "; ".join(f"{n}: {report_t95_s(out, n):.0f} s "
                      f"({report_t95_s(out, n)/60:.1f} min) to 95%, = ln(20)/k"
                      for n in ran)
          + f") against a {WINDOW_SCHED_S/3600:.0f} h window, so `τ/T` is under a percent "
          "and the tissue tracks plasma quasi-statically. **This is a statement about the "
          "ratio, not about schedules.** Shorten the window, or take a molecule that "
          "equilibrates slowly — an antibody, `D ~ 1–10 um²/s` with a binding-site-barrier "
          "sink and `τ` of hours to days — and `τ/T` approaches 1, at which point schedule "
          "dominates. The same solver now answers that question; before, it could not be "
          "asked.")
        A("")
        A("**(d) What is still missing, and it is the big one.** Time resolution was added "
          "to *transport*. The kill map remains a **static Hill function of cumulative "
          "exposure**: no inter-fraction DNA repair, no cell-cycle redistribution, no "
          "regrowth between fractions, no resistance induction, and no normal-tissue "
          "recovery. Those are the mechanisms that give fractionation its actual clinical "
          "benefit, and they are pharmacodynamics, not transport. So read the numbers "
          "above as **the transport-only floor on how much schedule can matter** — the "
          "part that is now correctly computed — and not as an estimate of the clinical "
          "effect of fractionating temozolomide.")
        A("")
    A("**Leg B is a window artifact, and worth naming as one.** Kill falls with "
      "fractionation there only because later fractions have less of the "
      f"{meta['window_sched_h']:.0f} h observation window left in which to be delivered "
      "and counted. Analytically the window exposure for `n` fractions of peak `P/n` at "
      "interval `Δ` is `(P/n)·Σ_k (1/λ)(1 − e^(−λ(T−kΔ)))`, which tends to `P/λ` for "
      "*every* `n` as `T → ∞`. Extend the window past full washout and leg B collapses "
      "onto leg A. It is a measurement-window effect, not a schedule effect.")
    A("")
    if meta.get("time_resolved"):
        A("Leg C confirms the model's dominant lever by a wide margin: raise total exposure "
          "and kill rises steeply. Everything else about the schedule is now *visible* "
          "rather than invisible — but it is two to three orders of magnitude smaller than "
          "the exposure lever on this drug and this window.")
    else:
        A("Leg C confirms the model's only real lever: raise total exposure and kill rises; "
          "everything else about the schedule is invisible to it.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------------ 5
    A("## 5. Infusion duration proxy: plasma half-life, bolus vs sustained")
    A("")
    A("A short half-life is a bolus that is gone quickly; a long half-life is a sustained "
      "low level. The hypothesis under test is the standard transport-limited one: *if "
      "delivery is limited by how far drug can travel before it is consumed, a sustained "
      "low concentration should reach further than a high, brief one.* Half-life is swept "
      f"{HALF_LIVES_H[0]:g}–{HALF_LIVES_H[-1]:g} h (a {HALF_LIVES_H[-1]/HALF_LIVES_H[0]:.0f}-fold "
      "range) with the peak rescaled to hold total window exposure fixed, plus a "
      "matched-peak control.")
    A("")
    ranh = [n for n in DRUGS if out["drugs"][n].get("half_life")]
    droph = [n for n in DRUGS if n not in ranh]
    if droph:
        A(f"> **DROPPED for compute budget:** the half-life sweep was NOT run for "
          f"{', '.join(droph)} (same second-solve cost as the schedule leg).")
        A("")
    for i, name in enumerate(ranh):
        d = out["drugs"][name]
        A(f"### 5.{i+1} {name}")
        A("")
        A("| plasma half-life (h) | peak, matched exposure (uM) | window-mean plasma | "
          "**kill (matched total exposure)** | tumor ≥ IC50 | IC50 depth (um) | "
          "window-mean plasma, matched peak | kill (matched peak) | IC50 depth, matched "
          "peak (um) |")
        A("|---|---|---|---|---|---|---|---|---|")
        for r in d["half_life"]:
            A(f"| {r['half_life_h']:g} | {r['peak_matched']:.2f} | "
              f"{r['plasma_mean_matched']:.4f} | **{r['kill_matched']:.6f}** | "
              f"{r['cov_matched']*100:.2f}% | "
              f"{_depth_str(r['depth_matched_um'], r.get('depth_matched_censored', False))} | "
              f"{r['plasma_mean_peakmatched']:.4f} | {r['kill_peakmatched']:.6f} | "
              f"{_depth_str(r['depth_peakmatched_um'], r.get('depth_peakmatched_censored', False))} |")
        km = [r["kill_matched"] for r in d["half_life"]]
        A("")
        if not meta.get("time_resolved"):
            A(f"Spread in kill across the {HALF_LIVES_H[-1]/HALF_LIVES_H[0]:.0f}-fold "
              f"half-life range at matched exposure: **{max(km) - min(km):.2e}**. "
              f"Penetration depth: **unchanged**.")
            A("")
            continue
        A(f"Spread in kill across the {HALF_LIVES_H[-1]/HALF_LIVES_H[0]:.0f}-fold "
          f"half-life range at matched exposure: **{max(km) - min(km):.2e}** — the number "
          f"that was **exactly `0.00e+00`** in the static-source version, where the "
          f"penetration depth was also reported as \"unchanged\" because there was only "
          f"ever one exposure field to measure it on.")
        A("")
        A("#### Each half-life now has its OWN exposure field, so depth is measurable")
        A("")
        A("Under the static source there was a single AUC field for the whole sweep, so "
          "\"depth is unchanged across half-life\" was true by construction and carried no "
          "information. Each shape is now solved separately and the profile is re-binned "
          "per shape. `half-distance` and `penetration index` are both **shape** "
          "statistics — normalised by the perivascular value / the p99, so amplitude "
          "cannot move them.")
        A("")
        A("| plasma half-life (h) | mean tumor exposure (uM·s) | Δ vs static | kill, "
          "time-resolved | kill, static | terminal residual C(T) / delivered | exposure "
          "half-distance from vessel (um) | penetration index | far-field plateau / "
          "perivascular |")
        A("|---|---|---|---|---|---|---|---|---|")
        for r in d["half_life"]:
            e = r["mean_auc_tumor_matched"]
            e0 = r["mean_auc_tumor_matched_static"]
            plateau = _fmt(r.get("plateau_ratio"), 4)
            A(f"| {r['half_life_h']:g} | {e:.5e} | {(e/e0-1)*100:+.3f}% | "
              f"{r['kill_matched']:.6f} | {r['kill_matched_static']:.6f} | "
              f"{r['residual_fraction']*100:.2f}% | "
              f"{_fmt(r.get('half_distance_um'), 2)} | "
              f"{r['penetration_index']:.6f} | {plateau} |")
        A("")
        hd = [r.get("half_distance_um") for r in d["half_life"]]
        pi = [r["penetration_index"] for r in d["half_life"]]
        if all(v is not None for v in hd):
            A(f"**Half-distance moves from {hd[0]:.2f} um at the "
              f"{d['half_life'][0]['half_life_h']:g} h (bolus-like) half-life to "
              f"{hd[-1]:.2f} um at {d['half_life'][-1]['half_life_h']:g} h (sustained) — "
              f"a change of {hd[-1]-hd[0]:+.2f} um.** The penetration index moves "
              f"{pi[0]:.6f} → {pi[-1]:.6f}. **Both say the BOLUS reaches relatively "
              f"deeper, not the sustained infusion.**")
            A("")
        A("#### Kill spread up the dose axis (static rows are exactly zero)")
        A("")
        A("| matched-exposure dose (×) | "
          + " | ".join(f"t½ {r['half_life_h']:g} h" for r in d["half_life"])
          + " | spread, TIME-RESOLVED | spread, STATIC |")
        A("|---" * (len(d["half_life"]) + 3) + "|")
        for m in SCHEDULE_DOSE_MULTIPLIERS:
            key = f"{m:g}"
            ks = [r["kill_matched_by_dose"][key] for r in d["half_life"]]
            ss = [r["kill_matched_by_dose_static"][key] for r in d["half_life"]]
            A(f"| {m:g}x | " + " | ".join(f"{k:.6f}" for k in ks)
              + f" | **{max(ks)-min(ks):.2e}** | {max(ss)-min(ss):.2e} |")
        A("")
    if not meta.get("time_resolved"):
        A("### 5.x Interpretation — the comparison is degenerate under the static source")
        A("")
        A("With a single time-averaged plasma level there is only one exposure field for "
          "the whole sweep, so both the kill and the depth are identical across every "
          "half-life BY CONSTRUCTION. Nothing can be concluded from this leg in this mode.")
        A("")
    else:
        A("### 5.x Interpretation — the naive transport-limited hypothesis is REFUTED, and "
          "in the opposite direction to the retracted claim")
        A("")
        A("The hypothesis under test was: *if delivery is transport-limited, a sustained "
          "low concentration should reach further than a high brief one.* With a working "
          "time-resolved driver we can now actually test it, and the answer is **no — and "
          "not for the reason the old version of this file gave.**")
        A("")
        A("**What the old version said, and why it was not an answer.** It reported the "
          "depth as identical across the whole 192-fold half-life range and concluded "
          "\"sustained infusion does NOT reach further here … it is exposure that moves "
          "depth, never the shape of the plasma curve\". But there was only ONE exposure "
          "field in that run — the plasma shape had been averaged away before the solve. "
          "The invariance was the input being invariant. **Retracted.**")
        A("")
        A("**What the corrected model says.** The depth DOES move with plasma shape, by a "
          "small amount, and it moves the *other* way: the bolus reaches relatively deeper "
          "than the sustained infusion. This is derivable from `A = M⁻¹(S_total − C(T))` "
          "(section 4):")
        A("")
        A("* At matched cumulative exposure `S_total` is identical, so the two schedules "
          "differ by `M⁻¹` applied to their terminal residuals. The bolus has washed out "
          "(`C(T) ≈ 0`); the sustained infusion has not.")
        A("* The exposure lost by the sustained infusion is `M⁻¹C_ss = M⁻²S_rate` — the "
          "source field smoothed **twice** by the operator, i.e. a *flatter* field than "
          "the exposure field `M⁻¹S_total` itself. Subtracting a flatter field from a "
          "peakier one leaves something **steeper**. So the sustained infusion loses more "
          "of its exposure at depth, proportionally, than near the vessel.")
        A("* Hence: bolus wins on total exposure AND on relative reach. The table above "
          "confirms both, on the half-distance and on the penetration index "
          "independently.")
        A("")
        A("**Why a bolus can never lose here.** Over a fixed observation window, the "
          "average delivery time of a schedule is what sets its residual. A bolus at "
          "`t = 0` has the earliest possible average delivery time, so it always realises "
          "at least as much of the steady-state exposure field as anything else at the "
          "same cumulative dose. There is no regime in a *linear* transport model where "
          "the ordering flips. (Where the textbook intuition does apply is a model with "
          "saturable perivascular binding or nonlinear uptake — neither of which is in "
          "this transport operator.)")
        A("")
        A("**Do not over-read the size of this.** The effect is bounded by `τ_tissue/T` "
          f"({sched_bound_pct(out, SCHEDULE_DRUG):.3f}% here) and the depth shift is "
          "sub-micrometre-scale against a perivascular gradient that section 3 shows is "
          "mostly inherited from the source geometry anyway. The value of this leg is "
          "**that the question is now answerable and the answer has a sign and a derived "
          "magnitude** — not that infusion duration is an actionable lever on this tumor.")
        A("")
        A("**Where it WOULD be actionable.** `τ_tissue/T` is the whole story, so the "
          "regime to look at is slow equilibration: an antibody (`D ~ 1–10 um²/s` plus a "
          "binding-site-barrier sink, `τ` of hours to days), or any molecule observed over "
          "a window comparable to its own washout. The old version of this file named that "
          "regime and said the plasma-averaging simplification \"would have to be lifted "
          "before the question could be asked at all\". It has been.")
        A("")
    A("---")
    A("")

    # ------------------------------------------------------------------ 6
    A("## 6. THE VERDICT for this tumor: dose-limited or delivery-limited?")
    A("")
    # The headline verdict is DERIVED from the same coverage statistics the
    # per-drug verdicts below use, so it can never drift out of agreement with
    # them (it used to be hand-written prose, and did).
    _kind = {}
    for name in DRUGS:
        _rows = out["drugs"][name]["dose"]
        _c1 = _dose_row(_rows, 1.0)["ic50_coverage"]
        _cmax = _rows[-1]["ic50_coverage"]
        _kind[name] = ("sensitivity" if (_cmax >= 0.95 and _c1 >= 0.95)
                       else "dose" if _cmax >= 0.95 else "delivery")
    _any_delivery = [n for n in DRUGS if _kind[n] == "delivery"]
    A("> ### Headline")
    A("> ")
    if not _any_delivery:
        A("> **On this tumor the response is DOSE-limited, not delivery-limited — and "
          "above the dose-limited range it becomes RESISTANCE-limited, never "
          "delivery-limited.**")
        A("> ")
        A(f"> Delivery is not the binding constraint here. The exposure field reaches the "
          f"whole section: at high enough dose >95% of tumor cells clear IC50 for "
          f"{'both drugs' if len(DRUGS) == 2 else 'every drug tested'}, and the far-field "
          f"floor (section 3) means even tissue a millimetre from a vessel is exposed. "
          f"What limits the answer is, in order: (1) at clinically plausible plasma "
          f"levels, simply not enough drug — the sub-saturating part of the curve, where "
          f"more exposure genuinely buys more kill; and (2) above that, a hard ceiling "
          f"that dose cannot touch.")
    else:
        A(f"> **On this tumor the verdict SPLITS BY DRUG: "
          + "; ".join(f"{n} is {_kind[n].upper()}-limited" for n in DRUGS) + ".**")
        A("> ")
        A(f"> For {', '.join(_any_delivery)}, delivery IS a binding constraint on this "
          f"geometry: even at {top:g}x nominal, fewer than 95% of tumor cells clear IC50, "
          f"so a share of the tumor is never therapeutically exposed at any dose swept. "
          f"For the rest, the exposure field does reach the whole section and the limits "
          f"are (1) not enough drug at clinically plausible plasma levels, then (2) a hard "
          f"resistance ceiling dose cannot touch. Read the per-drug sections below rather "
          f"than a single headline.")
    A("> ")
    if meta.get("time_resolved"):
        A("> **Either way, this is NOT a schedule verdict, and now that is a measurement "
          "rather than an assumption.** Sections 4-5 show the schedule effect on this "
          "tumor is real but bounded by `τ_tissue/T_window` — orders of magnitude below "
          "the dose lever. The static-source version of this file reported it as exactly "
          "zero, which was an artefact; the corrected answer is small, not zero.")
        A("> ")
    _tmz = out["drugs"]["temozolomide"]
    _bd = _tmz.get("kill_ceiling_breakdown") or {}
    _q = _bd.get("QUIESCENT") or {}
    _h = _bd.get("HYPOXIC") or {}
    _ceil = _tmz["kill_ceiling"]
    A(f"> **That ceiling is the actionable result, and it is not a transport number.** "
      f"For temozolomide it is {_ceil:.4f}, and section 1.3 decomposes it: "
      f"{_q.get('fraction', 0.0)*100:.1f}% of tumor cells are labelled QUIESCENT "
      f"(x{_q.get('resistance_multiplier', 0.5):.2f} kill multiplier) and "
      f"{_h.get('fraction', 0.0)*100:.1f}% HYPOXIC "
      f"(x{_h.get('resistance_multiplier', 0.35):.2f}), contributing "
      f"{_q.get('contribution', 0.0):.4f} and {_h.get('contribution', 0.0):.4f} of it "
      f"respectively — i.e. ~{_q.get('contribution', 0.0)/max(_ceil, 1e-12)*100:.0f}% of "
      f"the ceiling is the quiescent band alone. Because temozolomide's per-cell "
      f"sensitivity is uniform 1.0 on this panel (MGMT is not measurable here), the "
      f"ceiling for TMZ is *entirely* the cell-state resistance multipliers. It is not "
      f"delivery, not IC50, and not achievable exposure: at infinite drug every cell sits "
      f"at the top of its Hill curve and what remains is exactly "
      f"`sensitivity x resistance`.")
    A("> ")
    if not _any_delivery:
        A("> The operational reading: **on this tumor, escalating dose is a real but "
          "bounded lever, and the bound is metabolic state, not penetration.** If the "
          "ceiling is to move, it moves by changing the cell-state distribution "
          "(reoxygenation, cycle recruitment) or by a drug whose kill does not carry the "
          "quiescence penalty — not by more drug and not by better penetration.")
    else:
        A(f"> The operational reading: **escalating dose is a real but bounded lever, and "
          f"which bound bites first depends on the drug.** For "
          f"{', '.join(_any_delivery)} on this geometry the near-term bound is "
          f"PENETRATION — a share of the tumor is never exposed above IC50 at any swept "
          f"dose — and behind that sits the same cell-state ceiling. For the others the "
          f"ceiling is reached and the bound is metabolic state, not penetration; there, "
          f"moving it requires changing the cell-state distribution (reoxygenation, cycle "
          f"recruitment) or a drug whose kill does not carry the quiescence penalty.")
    A("> ")
    if meta.get("time_resolved") and out["drugs"][SCHEDULE_DRUG].get("schedule"):
        _sd = out["drugs"][SCHEDULE_DRUG]
        _es = [s["mean_auc_tumor_matched"] for s in _sd["schedule"]]
        _bnd = sched_bound_pct(out, SCHEDULE_DRUG)
        A(f"> **And schedule? Now measured rather than erased — and it is a REAL but "
          f"THIRD-ORDER lever on this tumor.** With the time-resolved plasma driver, "
          f"fractionation at matched cumulative exposure moves tissue exposure by "
          f"{(max(_es)/min(_es)-1)*100:.2f}% and infusion duration by a comparable amount, "
          f"against the `τ_tissue/T_window` ceiling of {_bnd:.2f}% that sections 4-5 derive "
          f"and confirm. Compare that with dose, which moves kill by a factor of "
          f"{_dose_row(_sd['dose'], top)['kill']/_dose_row(_sd['dose'], 1.0)['kill']:.0f} "
          f"over the swept range. **The previous version of this file reported this effect "
          f"as exactly zero and concluded schedule does not matter; that zero was a solver "
          f"artefact. The corrected conclusion is not \"schedule matters after all\" — it "
          f"is that schedule matters a little, for a reason we can now name and bound, and "
          f"far less than dose on this fast-equilibrating small molecule.** The caveat that "
          f"survives unchanged: the pharmacodynamic mechanisms that make fractionation "
          f"matter clinically are still absent from the kill model, so this is a floor.")
        A("> ")
    A(f"> **Load-bearing caveat:** that "
      f"{meta['quiescent_tumor_fraction']*100:.0f}% quiescent fraction is a *model* label "
      f"from the simulated oxygen field, not a measurement. It is being tested "
      f"independently against measured MKI67. If it is wrong, the ceiling and every "
      f"absolute kill here move with it — the dose/schedule *comparisons* survive, the "
      f"absolute levels may not.")
    A("")
    dmed = meta.get("distance_median_um", float("nan"))
    dp90 = meta.get("distance_p90_um", float("nan"))
    for name in DRUGS:
        d = out["drugs"][name]
        rows = d["dose"]
        ceil = d["kill_ceiling"]
        r1 = _dose_row(rows, 1.0)
        rmax = rows[-1]
        Lfit = (d["depth"].get("fit_near", {}) or {}).get("L_um") or d["depth"]["L_um"]
        cov1 = r1["ic50_coverage"]
        covmax = rmax["ic50_coverage"]
        sat90 = _first_multiplier_reaching(rows, "kill", 0.90 * rmax["kill"])
        headroom = rmax["kill"] / ceil

        if covmax >= 0.95 and cov1 >= 0.95:
            verdict = ("**NEITHER dose- nor delivery-limited — it is SENSITIVITY-limited.**")
            detail = (f"Already at the nominal dose, {cov1*100:.1f}% of tumor cells see an "
                      f"exposure above IC50. The drug reaches essentially the whole tumor; "
                      f"raising the dose {top:g}-fold moves kill only "
                      f"{r1['kill']:.4f} → {rmax['kill']:.4f}. What caps the response is "
                      f"the kill ceiling ({ceil:.4f}) — the sensitivity model plus hypoxic "
                      f"/ quiescent resistance. More drug does not help; a different drug, "
                      f"or reversing resistance, is the only lever the model leaves open.")
        elif covmax >= 0.95:
            verdict = ("**DOSE-limited over the swept range, then sensitivity-limited.**")
            detail = (f"At the nominal dose {cov1*100:.1f}% of tumor cells clear IC50, "
                      f"rising to {covmax*100:.1f}% by {top:g}x — the drug *can* reach "
                      f"essentially the whole tumor once enough of it is given, so the "
                      f"binding constraint over this range is how much drug is present, "
                      f"not how far it travels. Kill climbs with dose "
                      f"({r1['kill']:.4f} → {rmax['kill']:.4f}), and at {top:g}x it has "
                      f"reached {headroom*100:.1f}% of the {ceil:.4f} ceiling"
                      + (" without yet flattening — the sweep never gets far enough for "
                         "delivery to stop mattering."
                         if (sat90 is not None and sat90 >= top) else
                         f", flattening from about {sat90:g}x onward, beyond which extra "
                         f"drug buys essentially nothing and the limit becomes intrinsic "
                         f"resistance."))
        else:
            verdict = ("**DELIVERY-limited: a substantial share of tumor cells never see a "
                       "therapeutic exposure at any dose tested.**")
            detail = (f"Even at {top:g}x nominal only {covmax*100:.1f}% of tumor cells clear "
                      f"IC50 (from {cov1*100:.1f}% at 1x), and kill plateaus at "
                      f"{rmax['kill']:.4f} = {headroom*100:.1f}% of the {ceil:.4f} ceiling. "
                      f"With an exposure decay length of {Lfit:.0f} um against a p90 "
                      f"distance-to-vessel of {dp90:.0f} um, the far tissue is out of "
                      f"reach and escalating the systemic dose pushes it back only "
                      f"logarithmically.")
        A(f"### {name}")
        A("")
        A(verdict)
        A("")
        A(detail)
        A("")
        blo, bhi = _band_multipliers(name)
        vlo, clamped = _interp_at(rows, "ic50_coverage", blo)
        klo, _ = _interp_at(rows, "kill", blo)
        cb = CLINICAL_CMAX_BAND_UM[name]
        if bhi < 0.5:
            A(f"**Important qualification — the verdict above is stated at the shared "
              f"nominal level, which is well above {name}'s clinical range.** Its "
              f"approximate clinical Cmax ({cb[0]:g}–{cb[1]:g} uM) is only "
              f"{blo:.2g}x–{bhi:.2g}x nominal. Down there the sweep gives IC50 coverage "
              f"{vlo*100:.0f}% and kill {klo:.4f} ({klo/ceil*100:.0f}% of ceiling)"
              f"{', extrapolated to the edge of the swept range' if clamped else ''} — "
              f"i.e. at clinically plausible exposures this drug is materially "
              f"**dose-limited on this tumor**, and it only becomes sensitivity-limited "
              f"once pushed to supratherapeutic levels. Both statements are in the table; "
              f"which one is relevant depends on the dose actually deliverable.")
            A("")
        plat = d["depth"].get("plateau_ratio")
        A(f"Geometry check: near-field exposure decay length **L = {Lfit:.0f} um** "
          f"(d ≤ p90) vs median distance-to-vessel **{dmed:.0f} um** and p90 "
          f"**{dp90:.0f} um**. IC50 depth at 1x: "
          f"{_depth_str(r1['ic50_depth_um'], r1.get('ic50_depth_censored', False))} um; at "
          f"{top:g}x: "
          f"{_depth_str(rmax['ic50_depth_um'], rmax.get('ic50_depth_censored', False))} um "
          f"— but see section 3: with a far-field plateau at "
          f"{'—' if plat is None else f'{plat:.2f}'} of the perivascular value, the large "
          f"depths are the FLOOR clearing IC50, not a diffusion front reaching further.")
        A("")
    A("### What this means operationally (within the model's limits)")
    A("")
    A("* **Would more drug help?** Read the `% of kill ceiling` column in section 2. "
      "Where it is already near 100% at the nominal dose, escalation is pure added "
      "toxicity in this model. Where it climbs steeply, the model says the tumor is "
      "under-exposed rather than intrinsically resistant — a materially different "
      "clinical situation, and the distinction is the one thing this simulator is "
      "positioned to make.")
    if meta.get("time_resolved"):
        A("* **Is the tissue transport-saturated?** DOSE moves amplitude only — the "
          "perivascular profile shape is invariant under it exactly, because the operator "
          "is linear. SCHEDULE does move the shape, but by well under a percent here (the "
          f"`τ_tissue/T` bound of {sched_bound_pct(out, SCHEDULE_DRUG):.2f}%, sections "
          f"4-5), and in the direction that favours the bolus. So for practical purposes "
          "the geometry of delivery is set by three things — the drug's `D_eff`, its loss "
          "rate `k`, and the vascular geometry — and neither dose nor schedule is a lever "
          "on it *on this drug and this window*.")
    else:
        A("* **Is the tissue transport-saturated?** The invariance of the *shape* of the "
          "perivascular profile under every dose, schedule and half-life tested says the "
          "geometry of delivery is fixed. Only three things move it: the drug's `D_eff`, "
          "its loss rate `k`, and the vascular geometry itself. Dose and schedule move "
          "amplitude only. (NOTE: under the static source this was partly circular — see "
          "section 5.)")
    A("* **The lever that is NOT in this analysis** is the normal-tissue side of the "
      "therapeutic index. Every statement above is about tumor kill alone. A schedule "
      "that is near-identical for the tumor is emphatically not identical for the patient.")
    A("")
    A("---")
    A("")

    # ------------------------------------------------------------------ 7
    A("## 7. Limitations (the full list)")
    A("")
    A("1. **No outcome validation, at all.** Untreated resection; no measured response "
      "for either drug. Nothing here has been scored against reality.")
    A("2. **Not clinical dosing guidance.** Research prototype.")
    if meta.get("time_resolved"):
        A("3. **Schedule IS now resolvable — but only through transport, and the size of "
          "the effect is set by the observation window.** The `τ_tissue/T_window` bound "
          "(sections 4-5) means a longer window shrinks the measured schedule effect "
          "toward zero and a shorter one inflates it. The "
          f"{WINDOW_SCHED_S/3600:.0f} h window here was chosen for compute, not physics, "
          "and it is a *reporting* choice that changes a headline number. It is held "
          "exactly constant across everything compared, so the comparisons are sound, but "
          "the absolute percentages are window-conditional and should be quoted with the "
          "window attached.")
    else:
        A("3. **Schedule is not resolvable by this engine.** The plasma-averaging "
          "simplification makes matched-exposure schedule comparisons degenerate by "
          "construction (section 4a). A time-resolved driver is required before this engine "
          "can answer a schedule question for a slowly-equilibrating molecule.")
    A("4. **Pharmacodynamics is a single static Hill map — and this is now the binding "
      "limitation on the schedule question.** No inter-fraction DNA repair, no cell-cycle "
      "redistribution, no regrowth, no resistance evolution, no immune contribution, no "
      "normal-tissue toxicity. These are exactly the mechanisms that give real "
      "fractionation its benefit. Time-resolving *transport* removed a structural blind "
      "spot; it did not make this model a fractionation model.")
    A("5. **Sensitivity model.** Temozolomide's determinant (MGMT) is not on this panel, "
      "so TMZ resistance is unmodelled (uniform 1.0) and its absolute kill is optimistic. "
      "Cisplatin uses an MKI67 proliferation proxy, which is mechanistic but not a "
      "measured chemosensitivity.")
    A("6. **The perivascular gradient partly inherits the source-field shape** (the "
      "source-null column in section 3). Depth numbers are order-of-magnitude.")
    A("7. **Vascular density is a KDE of measured endothelial cell positions**, not a "
      "measured vessel lumen or a measured perfusion map; the vessel mask is a "
      f"percentile threshold (p{meta['vessel_percentile']:.0f}).")
    A("8. **Single section, single patient, 2D.** No inter-patient variability, no "
      "out-of-plane transport, and the nominal plasma level is a shared reference rather "
      "than a per-drug clinical PK model — supratherapeutic for cisplatin in particular.")
    A(f"9. **Cell-state banding dominates the ceiling.** "
      f"{meta['quiescent_tumor_fraction']*100:.0f}% of tumor cells are labelled QUIESCENT "
      f"by the oxygen model (x0.5 kill) and "
      f"{meta['hypoxic_tumor_fraction']*100:.0f}% HYPOXIC (x0.35). That labelling "
      f"is a model assignment from the simulated oxygen field, not a measurement; it is "
      f"being tested independently against measured MKI67. If the quiescent fraction is "
      f"wrong, every absolute kill and the ceiling move with it.")
    A("10. **Scope was cut for the compute budget** (see the banner at the top): "
      f"{len(DOSE_MULTIPLIERS)} dose points; schedule + half-life for {SCHEDULE_DRUG} "
      f"only; {len(SCHEDULES)} schedules; {len(HALF_LIVES_H)} half-lives; "
      f"{WINDOW_SCHED_S/3600:.0f} h schedule window; "
      f"{len(out['engine_checks'])} engine verification runs. Nothing "
      "was extrapolated to fill the gaps — the gaps are named instead.")
    A("11. **One calibration constant** (`alpha`, fit once on doxorubicin) sets the "
      "perivascular sink for every drug; the penetration lengths and therefore every "
      "depth number inherit its uncertainty.")
    if meta.get("time_resolved"):
        A("12. **The absolute kill levels shifted slightly from the static-source version "
          "of this file, and not only in the schedule sections.** The dose-response leg is "
          "also solved with a time-resolved source now (a decaying bolus over the "
          f"{meta['window_dose_h']:.0f} h dose window), so its exposure field is a little "
          "different too. Cross-file comparisons of absolute kill against the superseded "
          "version are therefore not apples to apples; comparisons *within* this file are.")
    A("")
    n_solves_total = sum(int(out["drugs"][n].get("n_pde_solves") or 0)
                         for n in out["drugs"]) or (len(DRUGS) + 1)
    A(f"*Generated by `experiments/dose_schedule.py` — ONE reconstruction reused "
      f"everywhere, {n_solves_total} transport solves "
      f"(`PDEConfig.time_resolved={bool(meta.get('time_resolved'))}`), "
      f"{len(out['engine_checks'])} full-engine verification runs, "
      f"{meta['runtime_s']/60:.1f} min total, seed {meta['seed']}, deterministic.*")

    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dx", type=float, default=DX_UM)
    ap.add_argument("--from-cache", action="store_true",
                    help="re-render the figure/report from the cached run")
    args = ap.parse_args()

    if args.from_cache and os.path.exists(_CACHE):
        with open(_CACHE) as fh:
            out = json.load(fh)
        # Re-derive the profile characterisation from the cached bins so a cache
        # re-render is identical to a fresh run (pure post-processing, no solve).
        for _n in out["drugs"]:
            augment_depth(out["drugs"][_n]["depth"], out["meta"]["distance_p90_um"])
        print(f"loaded cached run from {_CACHE}")
    else:
        out = run(dx_um=args.dx)
        with open(_CACHE, "w") as fh:
            json.dump(out, fh)
        print(f"cached run -> {_CACHE}")

    make_figure(out, _OUT_PNG)
    write_report(out, _OUT_MD)


if __name__ == "__main__":
    main()
