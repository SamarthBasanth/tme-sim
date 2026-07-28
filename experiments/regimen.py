"""REAL CLINICAL TEMOZOLOMIDE REGIMENS, and COMBINATION SCHEDULING, on one real tumor.

This is the drug-testing question the engine could not answer until the
time-resolved plasma driver landed. Until then ``sim.engine`` collapsed the whole
multi-dose plasma profile to a single time-averaged level before the transport
solver ever saw it, so **two regimens with the same cumulative dose were
bit-identical inputs**: the measured fractionation spread was exactly
``0.00e+00``, an arithmetic artefact rather than a result. Every schedule and
combination-scheduling conclusion drawn before that fix was comparing identical
inputs.

``PDEConfig.time_resolved = True`` is used for **every** solve in this script.
Nothing here runs on the static path.

What this script asks
---------------------
1. **REAL REGIMENS, not abstract sweeps.** The actual temozolomide schedules used
   in glioblastoma — Stupp concomitant, Stupp adjuvant, the 21/28 dose-dense arm
   of RTOG 0525, the 7-days-on/7-days-off metronomic schedule, continuous
   low-dose, and a single large dose as the control — converted to
   ``(plasma_peak, half_life, dose times)`` through an explicitly stated,
   citation-anchored PK conversion, and compared **at matched cumulative dose**
   over one 28-day cycle so schedule is the only variable.
2. **The resolvable-window ladder.** The schedule effect in a linear transport
   model is exactly the terminal residual and is bounded by
   ``tau_tissue / T_window`` (derived and confirmed in
   ``results/dose_schedule.md``: predicted 0.767%, measured 0.772% at a 12 h
   window). Over a 28-day cycle that bound is 0.014%. This script measures the
   effect at 24 h, 7 d, 14 d and 28 d windows and checks it against the bound.
3. **COMBINATIONS, re-run with time resolution** — temozolomide + cisplatin,
   temozolomide + doxorubicin, cisplatin + doxorubicin — under Bliss independence
   at the exposure level (SPATIAL complementarity only; no biochemical synergy is
   or can be modelled), with the killed-population overlap (Jaccard) reported.
4. **STAGGERED vs SIMULTANEOUS administration** — a question the static engine
   could not pose at all. Does delaying the second agent change *which cells* the
   pair kills?
5. **DOSE-MATCHED CONTROLS throughout.** Every combination is compared against
   doubling the better single agent, and against Bliss(better + ITSELF), or the
   result is just "more drug".

How the exposure fields are computed — EXACT, not an approximation
------------------------------------------------------------------
The transport PDE is linear in the source and starts from ``C = 0``, the
coefficients (``D_eff``, ``k``, ``v``) are time-invariant, and the source is
``S(x,y,t) = u(x,y) * plasma(t)``. Therefore the discrete solve is exactly linear
AND time-invariant in the plasma profile. ``plasma_profile`` builds a multi-dose
curve as a superposition of shifted single-dose exponentials, so for a regimen
with doses of peak ``P_k`` at times ``t_k`` observed over ``[0, T]``::

    AUC_field(regimen)  =  sum_k  P_k * A1(T - t_k)
    C(T)                =  sum_k  P_k * Cend1(T - t_k)

where ``A1(s)`` / ``Cend1(s)`` are the exposure field and terminal concentration
produced by ONE unit-peak dose observed over a window of length ``s``. Every dose
time here is a multiple of the solver step ``dt``, so this identity is exact on
the discrete grid — not interpolated, not extrapolated.

This is what makes a 28-day cycle affordable: ``A1(s)`` saturates once the plasma
curve has washed out (13 half-lives for temozolomide is 24 h), so the whole
28-day, 28-dose regimen is built from a handful of short solves. The
saturation is *certified*, not assumed: the residual plasma fraction at the
saturation window is asserted below a tolerance and the measured field difference
``|A1(24h) - A1(48h)| / A1(48h)`` is reported.

Three independent verifications pin it (each is a full extra PDE solve):
  V1  a DIRECT time-resolved solve of a genuine 2-dose q12h profile vs the
      superposition of ``A1(24h) + A1(12h)``;
  V2  the FULL unmodified ``sim.engine.simulate_response`` with
      ``time_resolved=True`` vs the fast-path kill fraction;
  V3  a manually time-SHIFTED profile (zeros, then the dose) vs ``A1(s - delay)``
      — the time-invariance half of the argument, tested on its own.

READ THIS BEFORE READING ANY NUMBER
-----------------------------------
* Research prototype. **NOT clinical dosing guidance.** No number here is a dose
  recommendation for any patient.
* Zero outcome validation: the Xenium section is an UNTREATED resection.
* Time resolution was added to TRANSPORT, not to pharmacodynamics. The kill map
  is still a static Hill function of window-averaged exposure — no inter-fraction
  DNA repair, no MGMT depletion, no cell-cycle redistribution, no regrowth, no
  immune contribution, no normal-tissue recovery. Those are exactly the
  mechanisms that give fractionation its clinical value. What is computed here is
  the TRANSPORT-ONLY floor.
* ``results/cell_state_impact.md`` establishes that absolute kill on this section
  is not defensible. Rankings and relative comparisons are.

Run:    ./.venv/bin/python experiments/regimen.py
Writes: results/regimen.md   (and a cache outside the repo)
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pickle
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from tmesim.cells.agents import update_cell_states
from tmesim.cells.pkpd import apply_drug_kill
from tmesim.core import (
    CellState,
    Microenvironment,
    OxygenConfig,
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
from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from tmesim.sim.engine import simulate_response
from tmesim.treatment.drugs import delivery_source, plasma_profile
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.validation.vessel_profile import (
    distance_from_vessel_um,
    vessel_mask_from_density,
)

# --------------------------------------------------------------------------- #
# Numerics
# --------------------------------------------------------------------------- #
#: Grid-converged resolution (results/numerical_rigor.md section 2). dx = 10 um
#: violates the upwind advection CFL condition at dt = 10 s.
DX_UM = 12.0
DT_S = 10.0
SEED = 0

#: EVERY solve in this script uses the time-resolved plasma driver. There is no
#: static-source path here; the static path is what made these questions
#: unanswerable in the first place.
TIME_RESOLVED = True

HOUR = 3600.0
DAY = 86400.0

# --------------------------------------------------------------------------- #
# Pharmacokinetics — stated explicitly, with the conversion assumptions named
# --------------------------------------------------------------------------- #
#: Temozolomide plasma elimination half-life. Well established at ~1.8 h.
TMZ_HALF_LIFE_H = 1.8

#: Temozolomide molar mass (g/mol).
TMZ_MW = 194.15

#: PK ANCHOR. Reported plasma AUC(0-inf) after a single ORAL dose of
#: TMZ_ANCHOR_DOSE_MG_M2 mg/m2 in adults, in ug*h/mL. Oral bioavailability is
#: ~100%, so the oral dose is treated as fully systemically available.
#: The model's 1-compartment instantaneous-absorption curve is anchored on AUC
#: (not on Cmax) because tissue exposure is driven by plasma AUC; see
#: `tmz_peak_for_dose` for the arithmetic and the assumptions this buries.
TMZ_ANCHOR_DOSE_MG_M2 = 150.0
TMZ_ANCHOR_AUC_UG_H_ML = 23.4     # geometric-mean AUC after 150 mg/m2 oral
TMZ_ANCHOR_CMAX_UG_ML = 7.5       # geometric-mean Cmax, same dose (cross-check)
#: Label-stated range of dose-linear kinetics (mg/m2/day). Doses outside it are
#: flagged in the report rather than silently extrapolated.
TMZ_LINEAR_RANGE_MG_M2 = (75.0, 250.0)

#: Combination leg, CONVENTION 1 ("shared"). The SHARED nominal peak and SHARED
#: plasma half-life used by results/combination.md, so the only thing that
#: differs between that report and this one is the time-resolved solver plus the
#: staggering question. NOT clinically dose-matched between agents.
COMBO_PEAK_UM = 50.0
COMBO_HALF_LIFE_H = 4.0
COMBO_DRUGS = ["temozolomide", "cisplatin", "doxorubicin"]

#: Combination leg, CONVENTION 2 ("clinical"). Each agent on its OWN measured
#: plasma half-life and an approximate clinical peak. Half-life is the parameter
#: staggering is most sensitive to, so running the stagger question under real
#: per-drug PK is the point of this second convention. Peaks are approximate and
#: are NOT dose-matched across agents; see the report for provenance and for
#: which of these are measured vs derived.
CLINICAL_PK = {
    "temozolomide": dict(half_life_h=1.8, peak_um=None),   # None -> from mg/m2
    "cisplatin":    dict(half_life_h=0.5, peak_um=8.0),    # free/ultrafilterable Pt
    "doxorubicin":  dict(half_life_h=30.0, peak_um=2.8),   # 35 mg/m2 IV bolus
}
CLINICAL_TMZ_DOSE_MG_M2 = 150.0
COMBO_PAIRS = [
    ("temozolomide", "cisplatin"),
    ("temozolomide", "doxorubicin"),
    ("cisplatin", "doxorubicin"),
]
#: Observation window for the staggering leg, and the stagger offset.
COMBO_WINDOW_S = 24.0 * HOUR
COMBO_STAGGER_S = 12.0 * HOUR
#: Short window reproducing results/combination.md's exposure window, so the
#: simultaneous baseline is directly comparable with that report.
COMBO_SHORT_WINDOW_S = 6.0 * HOUR
#: Dose-matched control multiplier on the better single agent.
DOSE_MATCH_FACTOR = 2.0

# --------------------------------------------------------------------------- #
# Single-dose kernel grids (the only PDE solves this script performs)
# --------------------------------------------------------------------------- #
#: Window lengths s (seconds) at which the unit single-dose response A1(s) is
#: solved, per drug. Every s that any regimen or combination arm needs is either
#: IN this grid or beyond the certified saturation point.
KERNEL_GRID_S = {
    # regimen leg + the clinically-anchored combination leg for temozolomide
    ("temozolomide", TMZ_HALF_LIFE_H): [6 * HOUR, 12 * HOUR, 18 * HOUR,
                                        24 * HOUR, 48 * HOUR],
    # combination leg, shared convention (results/combination.md)
    ("temozolomide", COMBO_HALF_LIFE_H): [6 * HOUR, 12 * HOUR, 24 * HOUR],
    ("cisplatin", COMBO_HALF_LIFE_H): [6 * HOUR, 12 * HOUR, 24 * HOUR],
    ("doxorubicin", COMBO_HALF_LIFE_H): [6 * HOUR, 12 * HOUR, 24 * HOUR],
    # combination leg, clinical per-drug PK
    ("cisplatin", 0.5): [12 * HOUR, 24 * HOUR],
    ("doxorubicin", 30.0): [12 * HOUR, 24 * HOUR],
}
#: A1(s) may be substituted by A1(SAT) for s > SAT only if the plasma fraction
#: still un-delivered at SAT is below this tolerance. Asserted, not assumed.
SAT_RESIDUAL_TOL = 1e-6

# --------------------------------------------------------------------------- #
# The real clinical regimens (one 28-day comparison window)
# --------------------------------------------------------------------------- #
CYCLE_DAYS = 28
CYCLE_WINDOW_S = CYCLE_DAYS * DAY

#: Each entry: dose in mg/m2 per administration, and the days of the 28-day
#: window on which it is given. Citations live in the report, not here.
REGIMENS: List[Dict] = [
    dict(key="stupp_adj_150", label="Stupp adjuvant, cycle 1",
         detail="150 mg/m2/d, days 1-5, q28d",
         mg_m2=150.0, days=list(range(0, 5))),
    dict(key="stupp_adj_200", label="Stupp adjuvant, cycles 2-6",
         detail="200 mg/m2/d, days 1-5, q28d",
         mg_m2=200.0, days=list(range(0, 5))),
    dict(key="concomitant_75", label="Stupp concomitant (with RT)",
         detail="75 mg/m2/d continuous, 7 d/wk",
         mg_m2=75.0, days=list(range(0, 28))),
    dict(key="dd_21_28_75", label="Dose-dense 21/28, cycle 1 (RTOG 0525)",
         detail="75 mg/m2/d, days 1-21, q28d",
         mg_m2=75.0, days=list(range(0, 21))),
    dict(key="dd_21_28_100", label="Dose-dense 21/28, cycles 2+ (RTOG 0525)",
         detail="100 mg/m2/d, days 1-21, q28d",
         mg_m2=100.0, days=list(range(0, 21))),
    dict(key="metronomic_7_7", label="Dose-intense 7 on / 7 off (Wick)",
         detail="150 mg/m2/d, days 1-7 and 15-21",
         mg_m2=150.0, days=list(range(0, 7)) + list(range(14, 21))),
    dict(key="director_a", label="DIRECTOR arm A (7/7)",
         detail="120 mg/m2/d, days 1-7 and 15-21",
         mg_m2=120.0, days=list(range(0, 7)) + list(range(14, 21))),
    dict(key="director_b", label="DIRECTOR arm B (21/28)",
         detail="80 mg/m2/d, days 1-21",
         mg_m2=80.0, days=list(range(0, 21))),
    dict(key="continuous_50", label="Continuous low-dose (RESCUE)",
         detail="50 mg/m2/d continuous",
         mg_m2=50.0, days=list(range(0, 28))),
    dict(key="single_bolus", label="Single large dose (CONTROL)",
         detail="whole cycle dose in ONE administration, day 1",
         mg_m2=1000.0, days=[0]),
]

#: Cumulative mg/m2 per cycle that the matched-dose comparison normalises to.
#: 1000 mg/m2 = the standard Stupp adjuvant cycle at 200 mg/m2 x 5 days.
MATCHED_CUMULATIVE_MG_M2 = 1000.0

#: Windows over which the burst-vs-continuous contrast is re-measured to test
#: the tau/T scaling of the schedule effect.
LADDER_WINDOWS_S = [7 * DAY, 14 * DAY, 28 * DAY]

#: A within-day fractionation ladder at a 24 h window, where tau/T is largest
#: and the effect is comfortably resolvable. Matched cumulative plasma exposure.
INTRADAY_WINDOW_S = 24.0 * HOUR
INTRADAY_SCHEDULES = [
    ("1 dose (bolus)", [0.0]),
    ("2 doses q12h", [0.0, 12 * HOUR]),
    ("4 doses q6h", [0.0, 6 * HOUR, 12 * HOUR, 18 * HOUR]),
]
#: The schedule effect is a fixed RELATIVE change in exposure, so the ABSOLUTE
#: kill spread depends on where you sit on the Hill curve. Evaluated at several
#: amplitudes — free, because exposure is exactly linear in amplitude.
LADDER_DOSE_MULTIPLIERS = [1.0, 10.0, 100.0, 1000.0]

_HERE = os.path.dirname(os.path.abspath(__file__))
_RESULTS_DIR = os.path.abspath(os.path.join(_HERE, "..", "results"))
_OUT_MD = os.path.join(_RESULTS_DIR, "regimen.md")
_DEFAULT_CACHE = os.path.join(
    os.environ.get("TMPDIR", "/tmp"), "tmesim_regimen_cache.pkl"
)


# --------------------------------------------------------------------------- #
# PK conversion — stated in one place so the assumptions cannot hide
# --------------------------------------------------------------------------- #
def tmz_lambda_s() -> float:
    """Plasma elimination rate constant, 1/s."""
    return math.log(2.0) / (TMZ_HALF_LIFE_H * HOUR)


def tmz_peak_for_dose(mg_m2: float) -> float:
    """Model plasma peak (uM) for an oral temozolomide dose in mg/m2.

    The conversion, in full:

    1. Oral bioavailability is ~100%, so the oral dose is treated as entirely
       systemically available and no absorption loss is applied.
    2. Reported plasma AUC(0-inf) after ``TMZ_ANCHOR_DOSE_MG_M2`` mg/m2 is
       ``TMZ_ANCHOR_AUC_UG_H_ML`` ug*h/mL. Converted to uM*s:
       ``AUC[uM*s] = AUC[ug*h/mL] * 1000 / MW * 3600``.
    3. Temozolomide PK is dose-linear over the clinical range (75-200 mg/m2), so
       AUC scales linearly with mg/m2.
    4. The model's plasma curve is 1-compartment with INSTANTANEOUS absorption
       and mono-exponential elimination at ``TMZ_HALF_LIFE_H``, for which
       ``AUC(0-inf) = peak / lambda``. Inverting gives
       ``peak = AUC * lambda``.

    What this buries, explicitly:

    * **Instantaneous absorption.** Real oral temozolomide has a tmax of roughly
      0.5-1 h. Anchoring on AUC rather than on Cmax means the model's peak is
      *higher* than the measured Cmax (the model front-loads a dose that really
      arrives over ~1 h). The direction of that error matters here: it makes
      every modelled dose MORE bolus-like than reality, which — per
      ``results/dose_schedule.md``, where the bolus wins — biases the schedule
      comparison in favour of the more concentrated regimens. The bias is
      identical across regimens at matched cumulative dose, so it cancels from
      the comparison and survives only in the absolute levels.
    * **AUC-matching, not Cmax-matching.** Tissue exposure in this model is
      driven by the plasma integral, not by the peak, so AUC is the right
      invariant to preserve. ``tmz_cmax_check`` reports what Cmax the model
      implies, against the measured Cmax, so the size of the discrepancy is
      visible rather than assumed away.
    * **No plasma protein binding, no active/inactive species split.**
      Temozolomide's cytotoxic species is MTIC, formed by pH-dependent hydrolysis;
      the model carries the parent drug's plasma curve only.
    * **Body-surface-area dosing is not converted to a per-patient dose.** mg/m2
      maps to a plasma level through a population-average PK relation.
    """
    auc_um_s = (TMZ_ANCHOR_AUC_UG_H_ML * 1000.0 / TMZ_MW) * 3600.0
    auc_um_s *= float(mg_m2) / TMZ_ANCHOR_DOSE_MG_M2
    return auc_um_s * tmz_lambda_s()


def tmz_cmax_check() -> Dict[str, float]:
    """Model-implied Cmax at the anchor dose vs the measured Cmax."""
    model_peak = tmz_peak_for_dose(TMZ_ANCHOR_DOSE_MG_M2)
    measured = TMZ_ANCHOR_CMAX_UG_ML * 1000.0 / TMZ_MW
    return dict(model_peak_um=model_peak, measured_cmax_um=measured,
                ratio=model_peak / measured if measured else float("nan"))


# --------------------------------------------------------------------------- #
# Engine-faithful helpers (orchestration only — nothing frozen is re-implemented)
# --------------------------------------------------------------------------- #
def compile_real_drug(name: str, alpha: float):
    """The panel's compilation pipeline, identical to real_drug_panel.py /
    dose_schedule.py: calibrated perivascular sink (alpha fit ONCE on
    doxorubicin) then the saturable binding-site barrier."""
    dp = calibrated_drugparams(name, alpha)
    mol = get_molecule(name)
    return apply_barrier_to_drugparams(dp, mol), mol


def prepare_tissue(micro: Microenvironment, cfg: SimConfig):
    """The engine's DRUG-INDEPENDENT prologue (steps 1-3 of simulate_response),
    run exactly once and reused by every arm."""
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


def single_dose_profile(half_life_s: float, window_s: float,
                        dt: float = DT_S) -> np.ndarray:
    """Unit-peak single-dose plasma curve on the solver's OWN step grid.

    The solver applies the source at the old time level of each step, i.e. at
    ``t_n = n*dt`` for ``n = 0 .. n_steps-1`` (``sim.engine._step_plasma_profile``
    documents exactly this), so the profile is evaluated there and nowhere else.
    """
    n = int(round(float(window_s) / float(dt)))
    times = np.arange(n, dtype=np.float64) * float(dt)
    pk = PKConfig(plasma_peak=1.0, half_life=float(half_life_s), n_doses=1,
                  dose_interval=max(float(window_s), 1.0))
    return plasma_profile(pk, times)


def unit_dose_transport(work, drug, velocity, half_life_s: float,
                        window_s: float, dt: float = DT_S) -> TransportResult:
    """Solve the transport PDE for ONE unit-peak dose over a window of length s.

    ``PDEConfig.time_resolved`` is set (and ``source_profile`` is passed, which is
    how the solver is actually told to resolve time). The returned exposure field
    is per 1 uM of plasma peak, so any amplitude is recovered by scaling.
    """
    n = int(round(float(window_s) / float(dt)))
    pde = PDEConfig(dt=dt, total_time=float(window_s), convection=True,
                    save_every=max(1, n), time_resolved=TIME_RESOLVED)
    source = delivery_source(work, drug, 1.0)
    return simulate_drug_transport(
        work, drug, source, pde,
        velocity=(velocity if pde.convection else None),
        source_profile=single_dose_profile(half_life_s, window_s, dt),
    )


class DoseKernel:
    """Single-dose exposure kernels ``A1(s)`` for one (drug, half-life) pair.

    ``field(s)`` returns ``(auc, c_end)`` for a unit-peak dose observed over a
    window of length ``s``. For ``s`` beyond the largest solved window ``S`` the
    kernel has saturated and ``A1(S)`` is returned in place of ``A1(s)``.

    That substitution is **certified, not assumed**, with an EXACT bound. A unit
    dose eventually deposits all of its drug, so ``k * A1(inf) = u_total / lambda``
    in mass terms. Over ``[0, S]`` only ``(1 - 2^(-S/t_half))`` of the plasma
    curve has been delivered, and of what was delivered a further ``C(S)`` is
    still sitting in the tissue uncounted. Hence

        (A1(inf) - A1(S)) / A1(inf)  =  2^(-S/t_half)  +  C(S) / (k * A1(S))
                                        \_ plasma _/      \_ tissue residual _/

    BOTH terms are computed from the solved fields (no estimate, no assumption)
    and their sum must be below ``SAT_RESIDUAL_TOL``. The plasma term alone is
    NOT sufficient: a molecule whose plasma has washed out may still be
    equilibrating in tissue, and that second term is exactly the terminal
    residual through which schedule acts in the first place.
    """

    def __init__(self, drug_name: str, half_life_s: float,
                 fields: Dict[int, Tuple[np.ndarray, np.ndarray]],
                 tissue: np.ndarray, k_loss: float):
        self.drug_name = drug_name
        self.half_life_s = float(half_life_s)
        self.k_loss = float(k_loss)
        self.fields = {int(k): v for k, v in fields.items()}
        self.grid = sorted(self.fields)
        self.sat_s = self.grid[-1]
        self.plasma_residual = 2.0 ** (-self.sat_s / self.half_life_s)
        a_sat, c_sat = self.fields[self.sat_s]
        eliminated = self.k_loss * float(np.asarray(a_sat)[tissue].sum())
        held = float(np.asarray(c_sat)[tissue].sum())
        self.tissue_residual = (held / eliminated) if eliminated > 0 else 0.0
        self.residual = self.plasma_residual + self.tissue_residual
        if len(self.grid) >= 2:
            a_hi = self.fields[self.grid[-1]][0]
            a_lo = self.fields[self.grid[-2]][0]
            denom = float(np.abs(a_hi[tissue]).sum())
            self.convergence = (float(np.abs(a_hi - a_lo)[tissue].sum()) / denom
                                if denom > 0 else 0.0)
            self.convergence_pair = (self.grid[-2], self.grid[-1])
        else:
            self.convergence = float("nan")
            self.convergence_pair = (0, self.grid[-1])

    def certify(self) -> None:
        if self.residual > SAT_RESIDUAL_TOL:
            raise AssertionError(
                f"{self.drug_name} (t1/2 {self.half_life_s/HOUR:.2f} h): A1 "
                f"saturation NOT certified at S = {self.sat_s / HOUR:.2f} h. "
                f"Un-delivered plasma {self.plasma_residual:.3g} + tissue still "
                f"holding {self.tissue_residual:.3g} = {self.residual:.3g} "
                f"exceeds the tolerance {SAT_RESIDUAL_TOL:.0g}. Extend "
                f"KERNEL_GRID_S rather than substituting a saturated field."
            )

    def field(self, s: float) -> Tuple[np.ndarray, np.ndarray]:
        key = int(round(float(s)))
        if key in self.fields:
            return self.fields[key]
        if key > self.sat_s:
            self.certify()
            return self.fields[self.sat_s]
        raise KeyError(
            f"{self.drug_name}: no single-dose kernel at s = {key / HOUR:.3f} h "
            f"and it is inside the saturation window; solved grid = "
            f"{[g / HOUR for g in self.grid]} h. Refusing to interpolate."
        )


def superpose(kernel: DoseKernel, dose_times_s: Sequence[float],
              peaks_um: Sequence[float], window_s: float,
              tissue: np.ndarray) -> TransportResult:
    """Exposure field for a multi-dose regimen, by exact superposition.

    ``AUC = sum_k P_k A1(T - t_k)`` and ``C(T) = sum_k P_k Cend1(T - t_k)``. Exact
    on the discrete grid because the discrete solve is linear and time-invariant
    in the plasma profile and every ``t_k`` is a multiple of ``dt``. Verified
    against direct solves (V1/V3) and against the full engine (V2).
    """
    auc = None
    c_end = None
    for t_k, p_k in zip(dose_times_s, peaks_um):
        s = float(window_s) - float(t_k)
        if s <= 0:
            continue
        a, c = kernel.field(s)
        if auc is None:
            auc = float(p_k) * a
            c_end = float(p_k) * c
        else:
            auc = auc + float(p_k) * a
            c_end = c_end + float(p_k) * c
    if auc is None:
        raise ValueError("regimen delivers no dose inside the window")
    auc_t = auc[tissue]
    peak = float(np.percentile(auc_t, 99.0)) if auc_t.size else 0.0
    pen = float(np.clip(auc_t.mean() / peak, 0.0, 1.0)) if peak > 0 else 0.0
    mean_end = float(c_end[tissue].mean()) if tissue.any() else 0.0
    return TransportResult(
        auc=auc, final_concentration=c_end,
        times=np.array([0.0, float(window_s)], dtype=np.float64),
        mean_concentration_history=np.array([0.0, mean_end], dtype=np.float64),
        penetration_index=pen,
    )


def window_plasma_auc(peaks_um: Sequence[float], dose_times_s: Sequence[float],
                      half_life_s: float, window_s: float,
                      dt: float = DT_S) -> float:
    """Cumulative plasma exposure COUNTED INSIDE the window (uM*s).

    Evaluated on the solver's own step grid so it is exactly the quantity the
    transport source integrates. This is what "matched cumulative exposure"
    means: matching it makes ``S_total`` identical between two schedules, which
    per ``results/dose_schedule.md`` leaves the terminal residual as the ONLY
    channel through which schedule can act.
    """
    n = int(round(float(window_s) / float(dt)))
    times = np.arange(n, dtype=np.float64) * float(dt)
    lam = math.log(2.0) / float(half_life_s)
    total = 0.0
    for t_k, p_k in zip(dose_times_s, peaks_um):
        el = times - float(t_k)
        active = el >= 0.0
        total += float(p_k) * float(
            np.where(active, np.exp(-lam * np.where(active, el, 0.0)), 0.0).sum()
        )
    return total * float(dt)


def per_cell(work, labeled_cells, tr: TransportResult, drug, pd_cfg: PDConfig,
             seed: int = SEED):
    """Per-cell kill PROBABILITY and exposure from the engine's OWN kill mapping.

    ``cells.pkpd.apply_drug_kill`` unmodified, on a throwaway copy of the
    oxygen-labelled pre-kill population, with ``default_rng(seed)`` — byte for
    byte the call the engine makes. Nothing is mutated.
    """
    snap = Microenvironment(
        domain=work.domain, fields=work.fields, cells=labeled_cells.copy(),
        patient_id=work.patient_id, metadata=work.metadata,
    )
    p = apply_drug_kill(snap, tr, pd_cfg, drug,
                        rng=np.random.default_rng(int(seed)))
    return (np.asarray(p, dtype=np.float64),
            np.asarray(snap.cells.exposure, dtype=np.float64),
            float(1.0 - snap.cells.viable_tumor_fraction()))


def arm_summary(work, labeled_cells, tr: TransportResult, drug,
                pd_cfg: PDConfig, window_s: float) -> Dict:
    """Kill, mean tumour AUC, IC50 coverage and the terminal residual."""
    p, expo, kill = per_cell(work, labeled_cells, tr, drug, pd_cfg)
    tumor = labeled_cells.tumor_mask()
    E = expo[tumor] / float(window_s)
    ic50 = float(drug.ic50)
    return dict(
        kill_drawn=float(kill),
        kill_expected=float(p[tumor].mean()),
        mean_auc_tumor=float(expo[tumor].mean()),
        mean_E_um=float(E.mean()),
        median_E_over_ic50=float(np.median(E) / ic50),
        ic50_coverage=float(np.mean(E >= ic50)),
        cov_10pct_ic50=float(np.mean(E >= 0.1 * ic50)),
        penetration_index=float(tr.penetration_index),
        p_tumor=p[tumor],
    )


def residual_fraction(tr: TransportResult, tissue: np.ndarray,
                      k_loss: float) -> float:
    """Fraction of delivered drug still IN the tissue when the window closes.

    The single channel through which schedule acts in a linear transport model
    (derivation in ``results/dose_schedule.md`` section 4b): integrating the PDE
    over ``[0, T]`` gives ``A = M^-1 (S_total - C(T))``, and ``S_total`` is equal
    by construction at matched cumulative exposure.
    """
    c_end = float(np.asarray(tr.final_concentration)[tissue].sum())
    eliminated = float(k_loss) * float(np.asarray(tr.auc)[tissue].sum())
    total = eliminated + c_end
    return float(c_end / total) if total > 0 else 0.0


def kill_ceiling(labeled_cells, pd_cfg: PDConfig):
    """Kill fraction reached as exposure -> infinity: `sensitivity x resistance`
    averaged over tumor cells. No dose can exceed it."""
    tumor = labeled_cells.tumor_mask()
    resistance = np.ones(labeled_cells.n, dtype=np.float64)
    resistance[labeled_cells.state == CellState.HYPOXIC] = pd_cfg.hypoxia_resistance
    resistance[labeled_cells.state == CellState.QUIESCENT] = pd_cfg.quiescent_resistance
    p = np.clip(pd_cfg.max_kill_rate * labeled_cells.sensitivity * resistance,
                0.0, 1.0)
    st = labeled_cells.state[tumor]
    pt = p[tumor]
    n_t = max(int(tumor.sum()), 1)
    bands = {}
    for label, code in (("VIABLE", None), ("HYPOXIC", CellState.HYPOXIC),
                        ("QUIESCENT", CellState.QUIESCENT)):
        sel = (~np.isin(st, [int(CellState.HYPOXIC), int(CellState.QUIESCENT)])
               if code is None else st == int(code))
        bands[label] = dict(
            fraction=float(sel.sum()) / n_t,
            mean_p=float(pt[sel].mean()) if sel.any() else 0.0,
            contribution=float(pt[sel].sum()) / n_t,
        )
    return float(pt.mean()), bands


# --------------------------------------------------------------------------- #
# Bliss independence + set statistics (exact expectations, no Monte Carlo)
# --------------------------------------------------------------------------- #
def bliss(p_a: np.ndarray, p_b: np.ndarray) -> np.ndarray:
    """``p = 1 - (1-p_A)(1-p_B)``. Independent action. NO synergy, NO antagonism."""
    return 1.0 - (1.0 - p_a) * (1.0 - p_b)


def intrinsic_factor(cells, pd_cfg: PDConfig) -> np.ndarray:
    """The DRUG-INDEPENDENT multiplier in the engine's PD map,
    ``sensitivity x R_hypoxia/R_quiescence``. Every agent is multiplied by this
    same per-cell number, which is exactly why raw cross-drug kill overlap is
    inflated and why the nulls below are needed to interpret a Jaccard."""
    resistance = np.ones(cells.n, dtype=np.float64)
    resistance[cells.state == CellState.HYPOXIC] = pd_cfg.hypoxia_resistance
    resistance[cells.state == CellState.QUIESCENT] = pd_cfg.quiescent_resistance
    return np.asarray(cells.sensitivity, dtype=np.float64) * resistance


N_SHUFFLE = 200
N_STRATA = 20
SHUFFLE_SEED = 20260728


def delivery_shuffle_null(p_a: np.ndarray, hill_b: np.ndarray,
                          intrinsic: np.ndarray, n_perm: int = N_SHUFFLE,
                          seed: int = SHUFFLE_SEED) -> Dict[str, float]:
    """Permutation null isolating DELIVERY co-localisation (as combination.md).

    Independent-Bernoulli overlap shrinks automatically as kill probabilities
    shrink, so a raw Jaccard proves nothing on its own. The engine's PD map
    factorises EXACTLY as ``p = p_hill(exposure) x intrinsic``, and the second
    factor is cell-intrinsic and SHARED by every drug. This null permutes only
    drug B's delivery term within strata of that shared factor, so each cell
    keeps its own intrinsic killability and only B's delivery GEOGRAPHY is
    scrambled. ``observed / null``: < 1 complementary, ~1 no spatial structure,
    > 1 redundant (both agents hit the same cells).
    """
    rng = np.random.default_rng(seed)
    q = np.quantile(intrinsic, np.linspace(0.0, 1.0, N_STRATA + 1)[1:-1])
    strata = np.searchsorted(q, intrinsic)
    order = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    jac = np.empty(n_perm)
    marg = np.empty(n_perm)
    for i in range(n_perm):
        hb = np.empty_like(hill_b)
        for idx in order:
            hb[idx] = rng.permutation(hill_b[idx])
        pb = np.clip(hb * intrinsic, 0.0, 1.0)
        s = set_stats(p_a, pb)
        jac[i] = s["jaccard"]
        marg[i] = s["kill_b"]
    return dict(jaccard_mean=float(jac.mean()), jaccard_sd=float(jac.std()),
                kill_b_mean=float(marg.mean()), n_perm=int(n_perm))


def set_stats(p_a: np.ndarray, p_b: np.ndarray) -> Dict[str, float]:
    """Exact expected sizes of the killed-cell sets under the independent
    coupling that Bliss independence IS. No Monte-Carlo noise."""
    n = float(p_a.size)
    e_a = float(p_a.sum())
    e_b = float(p_b.sum())
    e_both = float((p_a * p_b).sum())
    e_union = e_a + e_b - e_both
    return dict(
        kill_a=e_a / n, kill_b=e_b / n, kill_combo=e_union / n,
        jaccard=(e_both / e_union) if e_union > 0 else 0.0,
        overlap_coef=(e_both / min(e_a, e_b)) if min(e_a, e_b) > 0 else 0.0,
        frac_both=(e_both / e_union) if e_union > 0 else 0.0,
        frac_a_only=((e_a - e_both) / e_union) if e_union > 0 else 0.0,
        frac_b_only=((e_b - e_both) / e_union) if e_union > 0 else 0.0,
    )


# --------------------------------------------------------------------------- #
# COMPUTE STAGE — every PDE solve happens here
# --------------------------------------------------------------------------- #
def compute(dx_um: float = DX_UM) -> Dict:
    t_start = time.time()
    out: Dict = {"timing": {}}

    print(f"[1] Reconstructing the REAL Xenium GBM section ONCE (dx={dx_um} um) ...")
    t0 = time.time()
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
    report = reconstruction_report(micro)
    out["timing"]["reconstruct_s"] = time.time() - t0
    print(f"    {report['n_cells']:,} cells ({report['n_tumor']:,} tumor), "
          f"grid {report['grid']}, {out['timing']['reconstruct_s']:.1f} s "
          f"-- reused by EVERY arm below")
    out["report"] = {k: v for k, v in report.items()
                     if isinstance(v, (int, float, str, list, tuple))}

    alpha = float(calibrate_alpha())
    out["alpha"] = alpha
    print(f"[2] perivascular-sink alpha = {alpha:.2f} (fit ONCE on doxorubicin)")

    pd_cfg = PDConfig()
    ox_cfg = OxygenConfig()
    cfg0 = SimConfig(pde=PDEConfig(dt=DT_S, convection=True,
                                   time_resolved=TIME_RESOLVED),
                     oxygen=ox_cfg, pd=pd_cfg, seed=SEED)
    t0 = time.time()
    work, labeled, velocity = prepare_tissue(micro, cfg0)
    out["timing"]["prologue_s"] = time.time() - t0
    tissue = np.asarray(work.domain.tissue_mask, dtype=bool)
    tumor = labeled.tumor_mask()
    out["n_tumor"] = int(tumor.sum())

    s_intrinsic = intrinsic_factor(labeled, pd_cfg)[tumor]
    out["intrinsic_median"] = float(np.median(s_intrinsic))

    ceiling, bands = kill_ceiling(labeled, pd_cfg)
    out["kill_ceiling"] = ceiling
    out["ceiling_bands"] = bands
    print(f"    kill ceiling = {ceiling:.4f}   "
          f"(quiescent {bands['QUIESCENT']['fraction']*100:.1f}% of tumor cells, "
          f"hypoxic {bands['HYPOXIC']['fraction']*100:.1f}%)")

    # --- drug compilation ------------------------------------------------- #
    drugs = {}
    for name in COMBO_DRUGS:
        dp, mol = compile_real_drug(name, alpha)
        drugs[name] = dp
        out.setdefault("drug_params", {})[name] = dict(
            mw=float(mol.get("mw")), ic50=float(dp.ic50), hill=float(dp.hill),
            diffusivity=float(dp.diffusivity), bbb=float(dp.bbb_penetrance),
            k_loss=float(dp.clearance) + float(dp.binding),
        )
    out["drug_params"]["temozolomide"]["tau_tissue_s"] = (
        1.0 / out["drug_params"]["temozolomide"]["k_loss"])

    # --- [3] the ONLY PDE solves: single-dose kernels ---------------------- #
    print("[3] Solving the single-dose exposure kernels A1(s). "
          "PDEConfig(time_resolved=True) throughout.")
    kernels: Dict[Tuple[str, float], DoseKernel] = {}
    n_solves = 0
    total_steps = 0
    t_solves = time.time()
    for (name, hl_h), grid in KERNEL_GRID_S.items():
        fields = {}
        for s in grid:
            t0 = time.time()
            tr = unit_dose_transport(work, drugs[name], velocity,
                                     hl_h * HOUR, s)
            fields[int(round(s))] = (np.array(tr.auc),
                                     np.array(tr.final_concentration))
            n_solves += 1
            total_steps += int(round(s / DT_S))
            print(f"    A1[{name:13s} t1/2={hl_h:>4.1f} h] s={s/HOUR:>5.1f} h "
                  f"({int(round(s/DT_S)):>6,} steps)  {time.time()-t0:6.1f} s")
        k = DoseKernel(name, hl_h * HOUR, fields, tissue,
                       out["drug_params"][name]["k_loss"])
        kernels[(name, hl_h)] = k
        print(f"      -> saturation window {k.sat_s/HOUR:.0f} h; certified error "
              f"= plasma {k.plasma_residual:.3g} + tissue {k.tissue_residual:.3g} "
              f"= {k.residual:.3g}; "
              f"|A1({k.convergence_pair[0]/HOUR:.0f}h)-A1("
              f"{k.convergence_pair[1]/HOUR:.0f}h)|/A1 = {k.convergence:.3e}")
    out["timing"]["kernel_solves_s"] = time.time() - t_solves
    out["n_pde_solves"] = n_solves
    out["total_steps"] = total_steps
    out["kernel_meta"] = {
        f"{n}@{h}": dict(grid_h=[g / HOUR for g in kernels[(n, h)].grid],
                         sat_h=kernels[(n, h)].sat_s / HOUR,
                         residual=kernels[(n, h)].residual,
                         plasma_residual=kernels[(n, h)].plasma_residual,
                         tissue_residual=kernels[(n, h)].tissue_residual,
                         convergence=kernels[(n, h)].convergence,
                         convergence_pair_h=[kernels[(n, h)].convergence_pair[0] / HOUR,
                                             kernels[(n, h)].convergence_pair[1] / HOUR])
        for (n, h) in kernels
    }

    out["tmz_kernel_key"] = f"temozolomide@{TMZ_HALF_LIFE_H}"
    tmz_k = kernels[("temozolomide", TMZ_HALF_LIFE_H)]
    tmz = drugs["temozolomide"]
    tmz_hl_s = TMZ_HALF_LIFE_H * HOUR

    # --- [4] VERIFICATION -------------------------------------------------- #
    print("[4] Verifying the superposition identity against DIRECT solves.")
    ver = {}

    # V1: direct 2-dose q12h profile over 24 h vs A1(24h) + A1(12h).
    t0 = time.time()
    n24 = int(round(24 * HOUR / DT_S))
    times24 = np.arange(n24, dtype=np.float64) * DT_S
    prof2 = plasma_profile(PKConfig(plasma_peak=1.0, half_life=tmz_hl_s,
                                    n_doses=2, dose_interval=12 * HOUR), times24)
    pde24 = PDEConfig(dt=DT_S, total_time=24 * HOUR, convection=True,
                      save_every=n24, time_resolved=TIME_RESOLVED)
    tr_direct = simulate_drug_transport(
        work, tmz, delivery_source(work, tmz, 1.0), pde24,
        velocity=velocity, source_profile=prof2)
    tr_super = superpose(tmz_k, [0.0, 12 * HOUR], [1.0, 1.0], 24 * HOUR, tissue)
    d = np.abs(tr_direct.auc - tr_super.auc)[tissue]
    ref = np.abs(tr_direct.auc)[tissue]
    ver["v1"] = dict(
        label="direct 2-dose q12h solve vs A1(24h)+A1(12h), 24 h window",
        max_abs=float(d.max()), rel_l1=float(d.sum() / ref.sum()),
        direct_mean_auc=float(tr_direct.auc[tissue].mean()),
        super_mean_auc=float(tr_super.auc[tissue].mean()),
        seconds=time.time() - t0)
    n_solves += 1
    total_steps += n24
    print(f"    V1 {ver['v1']['label']}: relative L1 = {ver['v1']['rel_l1']:.3e}")

    # V3: manually time-SHIFTED single dose vs A1(window - delay). Tests the
    # time-invariance half of the argument on its own, on the drug where the
    # stagger leg actually depends on it.
    t0 = time.time()
    dox = drugs["doxorubicin"]
    dox_hl_s = COMBO_HALF_LIFE_H * HOUR
    delay_steps = int(round(12 * HOUR / DT_S))
    shifted = np.concatenate([
        np.zeros(delay_steps, dtype=np.float64),
        single_dose_profile(dox_hl_s, 12 * HOUR),
    ])
    tr_shift = simulate_drug_transport(
        work, dox, delivery_source(work, dox, 1.0), pde24,
        velocity=velocity, source_profile=shifted)
    a_ref, _ = kernels[("doxorubicin", COMBO_HALF_LIFE_H)].field(12 * HOUR)
    d = np.abs(tr_shift.auc - a_ref)[tissue]
    ver["v3"] = dict(
        label="doxorubicin dose DELAYED 12 h in a 24 h solve vs A1(12h)",
        max_abs=float(d.max()),
        rel_l1=float(d.sum() / np.abs(a_ref)[tissue].sum()),
        seconds=time.time() - t0)
    n_solves += 1
    total_steps += n24
    print(f"    V3 {ver['v3']['label']}: relative L1 = {ver['v3']['rel_l1']:.3e}")

    # V2: the FULL unmodified engine, time_resolved=True, vs the fast path.
    t0 = time.time()
    peak_v2 = tmz_peak_for_dose(150.0)
    cfg_v2 = SimConfig(
        drug=tmz,
        pk=PKConfig(plasma_peak=peak_v2, half_life=tmz_hl_s, n_doses=2,
                    dose_interval=12 * HOUR),
        pde=PDEConfig(dt=DT_S, total_time=24 * HOUR, convection=True,
                      save_every=n24, time_resolved=TIME_RESOLVED),
        oxygen=ox_cfg, pd=pd_cfg, seed=SEED)
    res_v2 = simulate_response(micro, cfg_v2)
    tr_fast = superpose(tmz_k, [0.0, 12 * HOUR], [peak_v2, peak_v2],
                        24 * HOUR, tissue)
    fast = arm_summary(work, labeled, tr_fast, tmz, pd_cfg, 24 * HOUR)
    d2 = np.abs(np.asarray(res_v2.transport.auc) - tr_fast.auc)[tissue]
    ver["v2"] = dict(
        label="FULL sim.engine.simulate_response (time_resolved=True), "
              "temozolomide 2 doses q12h, 24 h window",
        engine_kill=float(res_v2.kill_fraction),
        fast_kill=float(fast["kill_drawn"]),
        abs_delta=abs(float(res_v2.kill_fraction) - float(fast["kill_drawn"])),
        engine_mean_auc=float(res_v2.extras["mean_auc"]),
        fast_mean_auc=float(tr_fast.auc[tissue].mean()),
        auc_rel_l1=float(d2.sum()
                         / np.abs(np.asarray(res_v2.transport.auc))[tissue].sum()),
        peak_um=float(peak_v2),
        seconds=time.time() - t0)
    n_solves += 1
    total_steps += n24
    print(f"    V2 engine kill {ver['v2']['engine_kill']:.10f} vs fast-path "
          f"{ver['v2']['fast_kill']:.10f}  delta {ver['v2']['abs_delta']:.1e}")
    out["verification"] = ver
    out["n_pde_solves"] = n_solves
    out["total_steps"] = total_steps

    # --- [5] REAL CLINICAL REGIMENS, 28-day cycle -------------------------- #
    print("[5] Real clinical regimens over one 28-day cycle.")
    out["pk_conversion"] = dict(
        half_life_h=TMZ_HALF_LIFE_H, mw=TMZ_MW,
        anchor_dose_mg_m2=TMZ_ANCHOR_DOSE_MG_M2,
        anchor_auc_ug_h_ml=TMZ_ANCHOR_AUC_UG_H_ML,
        anchor_cmax_ug_ml=TMZ_ANCHOR_CMAX_UG_ML,
        peak_per_150=tmz_peak_for_dose(150.0),
        **tmz_cmax_check())

    regimen_rows = []
    for reg in REGIMENS:
        days = reg["days"]
        dose_times = [float(d) * DAY for d in days]
        cum_prescribed = reg["mg_m2"] * len(days)
        for mode in ("prescribed", "matched"):
            if mode == "prescribed":
                mg = reg["mg_m2"]
            else:
                mg = MATCHED_CUMULATIVE_MG_M2 / len(days)
            peak = tmz_peak_for_dose(mg)
            peaks = [peak] * len(days)
            tr = superpose(tmz_k, dose_times, peaks, CYCLE_WINDOW_S, tissue)
            summ = arm_summary(work, labeled, tr, tmz, pd_cfg, CYCLE_WINDOW_S)
            summ.pop("p_tumor")
            lo, hi = TMZ_LINEAR_RANGE_MG_M2
            regimen_rows.append(dict(
                key=reg["key"], label=reg["label"], detail=reg["detail"],
                mode=mode, n_doses=len(days),
                outside_linear_range=bool(mg < lo or mg > hi),
                mg_m2_per_dose=mg, cumulative_mg_m2=mg * len(days),
                cumulative_mg_m2_prescribed=cum_prescribed,
                peak_um=peak,
                plasma_auc_um_s=window_plasma_auc(peaks, dose_times, tmz_hl_s,
                                                  CYCLE_WINDOW_S),
                residual=residual_fraction(
                    tr, tissue, out["drug_params"]["temozolomide"]["k_loss"]),
                **summ))
            print(f"    {reg['key']:18s} {mode:10s} kill={summ['kill_drawn']:.6f} "
                  f"meanAUC={summ['mean_auc_tumor']:.5g}")
    # The matched-dose comparison, as relative deviations from the arm with the
    # largest tumour exposure. This is the headline: at matched cumulative dose
    # the ONLY thing that can differ between regimens is the terminal residual.
    matched_rows = [r for r in regimen_rows if r["mode"] == "matched"]
    ref_auc = max(r["mean_auc_tumor"] for r in matched_rows)
    for r in matched_rows:
        r["auc_rel_to_best"] = (r["mean_auc_tumor"] - ref_auc) / ref_auc
    out["matched_exposure_spread"] = (
        (ref_auc - min(r["mean_auc_tumor"] for r in matched_rows)) / ref_auc)
    out["matched_bound_tau_over_T"] = (
        out["drug_params"]["temozolomide"]["tau_tissue_s"] / CYCLE_WINDOW_S)
    out["regimen_rows"] = regimen_rows
    print(f"    MATCHED-DOSE exposure spread across all regimens = "
          f"{out['matched_exposure_spread']*100:.6f}%  "
          f"(bound tau/T = {out['matched_bound_tau_over_T']*100:.6f}%)")

    # matched-dose regimens rescaled up the dose axis (free: exact linearity)
    matched = [r for r in regimen_rows if r["mode"] == "matched"]
    amp_rows = []
    for mult in LADDER_DOSE_MULTIPLIERS:
        row = {"mult": mult}
        for reg in REGIMENS:
            days = reg["days"]
            dose_times = [float(d) * DAY for d in days]
            mg = MATCHED_CUMULATIVE_MG_M2 / len(days)
            peak = tmz_peak_for_dose(mg) * mult
            tr = superpose(tmz_k, dose_times, [peak] * len(days),
                           CYCLE_WINDOW_S, tissue)
            s = arm_summary(work, labeled, tr, tmz, pd_cfg, CYCLE_WINDOW_S)
            row[reg["key"]] = dict(kill=s["kill_drawn"],
                                   kill_expected=s["kill_expected"],
                                   mean_auc=s["mean_auc_tumor"],
                                   ic50_cov=s["ic50_coverage"])
        amp_rows.append(row)
    out["regimen_amplitude_rows"] = amp_rows

    # --- [6] the tau/T window ladder --------------------------------------- #
    print("[6] Window ladder: burst vs continuous at matched cumulative dose.")
    ladder = []
    for T in LADDER_WINDOWS_S:
        n_days = int(round(T / DAY))
        arms = {
            "burst 5/28 (Stupp adjuvant)": [d * DAY for d in range(0, 5)],
            "continuous daily": [d * DAY for d in range(0, n_days)],
            "single bolus": [0.0],
        }
        entry = {"window_s": T, "arms": {}}
        target_auc = None
        for label, times_ in arms.items():
            unit_auc = window_plasma_auc([1.0] * len(times_), times_, tmz_hl_s, T)
            if target_auc is None:
                target_auc = unit_auc * tmz_peak_for_dose(
                    MATCHED_CUMULATIVE_MG_M2 / 5.0)
            peak = target_auc / unit_auc
            tr = superpose(tmz_k, times_, [peak] * len(times_), T, tissue)
            s = arm_summary(work, labeled, tr, tmz, pd_cfg, T)
            s.pop("p_tumor")
            entry["arms"][label] = dict(
                n_doses=len(times_), peak_um=peak,
                plasma_auc_um_s=window_plasma_auc([peak] * len(times_), times_,
                                                  tmz_hl_s, T),
                residual=residual_fraction(
                    tr, tissue, out["drug_params"]["temozolomide"]["k_loss"]),
                **s)
        aucs = [v["mean_auc_tumor"] for v in entry["arms"].values()]
        entry["exposure_spread_rel"] = (max(aucs) - min(aucs)) / max(aucs)
        entry["bound_tau_over_T"] = (
            out["drug_params"]["temozolomide"]["tau_tissue_s"] / T)
        ladder.append(entry)
        print(f"    T={T/DAY:>5.1f} d  exposure spread = "
              f"{entry['exposure_spread_rel']*100:.5f}%  "
              f"bound tau/T = {entry['bound_tau_over_T']*100:.5f}%")
    out["window_ladder"] = ladder

    # --- [7] intra-day fractionation, where the effect is resolvable ------- #
    print("[7] Intra-day fractionation ladder (24 h window, matched exposure).")
    T = INTRADAY_WINDOW_S
    base_peak = tmz_peak_for_dose(MATCHED_CUMULATIVE_MG_M2 / 5.0)
    target_auc = window_plasma_auc([base_peak], [0.0], tmz_hl_s, T)
    intraday = {"window_s": T, "target_plasma_auc": target_auc,
                "bound_tau_over_T": (
                    out["drug_params"]["temozolomide"]["tau_tissue_s"] / T),
                "arms": {}, "amplitude": []}
    trs = {}
    for label, times_ in INTRADAY_SCHEDULES:
        unit = window_plasma_auc([1.0] * len(times_), times_, tmz_hl_s, T)
        peak = target_auc / unit
        tr = superpose(tmz_k, times_, [peak] * len(times_), T, tissue)
        trs[label] = (tr, times_, peak)
        s = arm_summary(work, labeled, tr, tmz, pd_cfg, T)
        s.pop("p_tumor")
        intraday["arms"][label] = dict(
            n_doses=len(times_), peak_um=peak,
            plasma_auc_um_s=window_plasma_auc([peak] * len(times_), times_,
                                              tmz_hl_s, T),
            residual=residual_fraction(
                tr, tissue, out["drug_params"]["temozolomide"]["k_loss"]),
            **s)
    aucs = [v["mean_auc_tumor"] for v in intraday["arms"].values()]
    intraday["exposure_spread_rel"] = (max(aucs) - min(aucs)) / max(aucs)
    for mult in LADDER_DOSE_MULTIPLIERS:
        row = {"mult": mult, "kill": {}}
        for label, (tr, times_, peak) in trs.items():
            scaled = TransportResult(
                auc=tr.auc * mult, final_concentration=tr.final_concentration * mult,
                times=tr.times,
                mean_concentration_history=tr.mean_concentration_history * mult,
                penetration_index=tr.penetration_index)
            s = arm_summary(work, labeled, scaled, tmz, pd_cfg, T)
            row["kill"][label] = dict(drawn=s["kill_drawn"],
                                      expected=s["kill_expected"])
        ks = [v["drawn"] for v in row["kill"].values()]
        ke = [v["expected"] for v in row["kill"].values()]
        row["spread_drawn"] = max(ks) - min(ks)
        row["spread_expected"] = max(ke) - min(ke)
        intraday["amplitude"].append(row)
    out["intraday"] = intraday
    print(f"    exposure spread = {intraday['exposure_spread_rel']*100:.4f}%  "
          f"bound tau/T = {intraday['bound_tau_over_T']*100:.4f}%")

    # --- [8] COMBINATIONS: simultaneous vs staggered ----------------------- #
    print("[8] Combinations under Bliss independence, simultaneous vs staggered.")
    combo: Dict = {"short_window_s": COMBO_SHORT_WINDOW_S,
                   "window_s": COMBO_WINDOW_S,
                   "stagger_s": COMBO_STAGGER_S,
                   "peak_um": COMBO_PEAK_UM,
                   "half_life_h": COMBO_HALF_LIFE_H}

    # 8a. single agents at the short (combination.md) window, simultaneous.
    singles_short = {}
    for name in COMBO_DRUGS:
        k = kernels[(name, COMBO_HALF_LIFE_H)]
        for mult in (1.0, DOSE_MATCH_FACTOR):
            tr = superpose(k, [0.0], [COMBO_PEAK_UM * mult],
                           COMBO_SHORT_WINDOW_S, tissue)
            s = arm_summary(work, labeled, tr, drugs[name], pd_cfg,
                            COMBO_SHORT_WINDOW_S)
            singles_short[(name, mult)] = s
    combo["singles_short"] = {
        f"{n}@{m}": {k: v for k, v in s.items() if k != "p_tumor"}
        for (n, m), s in singles_short.items()}

    # 8b. the stagger leg at the 24 h window.
    #   simultaneous : both agents dosed at t = 0
    #   staggered    : the second agent dosed at t = COMBO_STAGGER_S
    #   Two dose-matching conventions, because they answer different questions:
    #     'administered' — same peak in both arms. The delayed agent then has
    #        less window left in which to be delivered AND COUNTED, so it loses
    #        exposure. That is a MEASUREMENT-WINDOW artefact, named as such in
    #        results/dose_schedule.md (its 'leg B').
    #     'counted'      — the delayed agent's peak is rescaled so its plasma AUC
    #        INSIDE the window matches the simultaneous arm. S_total is then
    #        identical and the only remaining difference is the terminal
    #        residual. This is the leg that isolates SCHEDULE.
    #
    # Run under TWO PK conventions, because staggering is most sensitive to the
    # plasma half-life and it would be dishonest to answer the question only on
    # the shared 4 h placeholder:
    #   'shared'   — 50 uM peak, 4 h half-life for every agent (combination.md's
    #                convention, so this leg is directly comparable with it).
    #   'clinical' — each agent on its own measured plasma half-life and an
    #                approximate clinical peak. Half-lives here span 60x
    #                (cisplatin free platinum 0.5 h vs doxorubicin 30 h), which
    #                is exactly the axis a stagger acts on.
    conventions = {
        "shared": {n: dict(half_life_h=COMBO_HALF_LIFE_H, peak_um=COMBO_PEAK_UM)
                   for n in COMBO_DRUGS},
        "clinical": {
            n: dict(half_life_h=CLINICAL_PK[n]["half_life_h"],
                    peak_um=(CLINICAL_PK[n]["peak_um"]
                             if CLINICAL_PK[n]["peak_um"] is not None
                             else tmz_peak_for_dose(CLINICAL_TMZ_DOSE_MG_M2)))
            for n in COMBO_DRUGS},
    }
    combo["conventions"] = {
        c: {n: dict(v) for n, v in spec.items()} for c, spec in conventions.items()}

    for conv, spec in conventions.items():
        arms_single = {}
        for name in COMBO_DRUGS:
            hl_h = spec[name]["half_life_h"]
            peak = spec[name]["peak_um"]
            k = kernels[(name, hl_h)]
            hl_s = hl_h * HOUR
            auc_full = window_plasma_auc([peak], [0.0], hl_s, COMBO_WINDOW_S)
            auc_late_unit = window_plasma_auc([1.0], [COMBO_STAGGER_S], hl_s,
                                              COMBO_WINDOW_S)
            variants = {
                "t0": ([0.0], [peak]),
                "t0_x2": ([0.0], [peak * DOSE_MATCH_FACTOR]),
                "late_administered": ([COMBO_STAGGER_S], [peak]),
                "late_counted": ([COMBO_STAGGER_S], [auc_full / auc_late_unit]),
            }
            for vk, (times_, peaks_) in variants.items():
                tr = superpose(k, times_, peaks_, COMBO_WINDOW_S, tissue)
                s = arm_summary(work, labeled, tr, drugs[name], pd_cfg,
                                COMBO_WINDOW_S)
                s["plasma_auc_um_s"] = window_plasma_auc(peaks_, times_, hl_s,
                                                         COMBO_WINDOW_S)
                s["peak_um"] = float(peaks_[0])
                s["half_life_h"] = float(hl_h)
                s["residual"] = residual_fraction(
                    tr, tissue, out["drug_params"][name]["k_loss"])
                arms_single[(name, vk)] = s
        combo[f"singles_24h_{conv}"] = {
            f"{n}|{v}": {kk: val for kk, val in s.items() if kk != "p_tumor"}
            for (n, v), s in arms_single.items()}

        pairs = []
        for a, b in COMBO_PAIRS:
            entry = {"a": a, "b": b, "arms": {}}
            if conv == "shared":
                # short-window baseline, directly comparable with combination.md
                pa = singles_short[(a, 1.0)]["p_tumor"]
                pb = singles_short[(b, 1.0)]["p_tumor"]
                st = set_stats(pa, pb)
                better = a if pa.mean() >= pb.mean() else b
                p_better = singles_short[(better, 1.0)]["p_tumor"]
                p_better2 = singles_short[(better, DOSE_MATCH_FACTOR)]["p_tumor"]
                hill_b_short = np.divide(pb, s_intrinsic, out=np.zeros_like(pb),
                                         where=s_intrinsic > 0)
                null_short = delivery_shuffle_null(pa, hill_b_short, s_intrinsic)
                entry["short"] = dict(
                    better=better,
                    kill_2x_better=float(p_better2.mean()),
                    kill_selfcombo=float(bliss(p_better, p_better).mean()),
                    null_jaccard=null_short["jaccard_mean"],
                    obs_over_null=(st["jaccard"] / null_short["jaccard_mean"]
                                   if null_short["jaccard_mean"] > 0
                                   else float("nan")),
                    **st)
            for label, ka, kb in (
                ("simultaneous", (a, "t0"), (b, "t0")),
                ("B delayed, matched administered",
                 (a, "t0"), (b, "late_administered")),
                ("B delayed, matched counted exposure",
                 (a, "t0"), (b, "late_counted")),
                ("A delayed, matched counted exposure",
                 (a, "late_counted"), (b, "t0")),
            ):
                p_a = arms_single[ka]["p_tumor"]
                p_b = arms_single[kb]["p_tumor"]
                s = set_stats(p_a, p_b)
                hill_b = np.divide(p_b, s_intrinsic, out=np.zeros_like(p_b),
                                   where=s_intrinsic > 0)
                null = delivery_shuffle_null(p_a, hill_b, s_intrinsic)
                s["null_jaccard"] = null["jaccard_mean"]
                s["null_jaccard_sd"] = null["jaccard_sd"]
                s["obs_over_null"] = (s["jaccard"] / null["jaccard_mean"]
                                      if null["jaccard_mean"] > 0
                                      else float("nan"))
                s["null_marginal_err"] = abs(null["kill_b_mean"] - s["kill_b"])
                entry["arms"][label] = s
            p_a0 = arms_single[(a, "t0")]["p_tumor"]
            p_b0 = arms_single[(b, "t0")]["p_tumor"]
            better24 = a if p_a0.mean() >= p_b0.mean() else b
            entry["better_24h"] = better24
            entry["kill_2x_better_24h"] = float(
                arms_single[(better24, "t0_x2")]["p_tumor"].mean())
            pbest = arms_single[(better24, "t0")]["p_tumor"]
            entry["kill_selfcombo_24h"] = float(bliss(pbest, pbest).mean())
            pairs.append(entry)
            sim = entry["arms"]["simultaneous"]
            stag = entry["arms"]["B delayed, matched counted exposure"]
            print(f"    [{conv:8s}] {a} + {b}: Jaccard simultaneous "
                  f"{sim['jaccard']:.6f} -> staggered {stag['jaccard']:.6f} "
                  f"(delta {stag['jaccard']-sim['jaccard']:+.2e})")
        combo[f"pairs_{conv}"] = pairs
    out["combination"] = combo

    out["timing"]["total_s"] = time.time() - t_start
    print(f"\nTOTAL {out['timing']['total_s']/60:.1f} min, "
          f"{out['n_pde_solves']} PDE solves, {out['total_steps']:,} solver steps.")
    return out


# --------------------------------------------------------------------------- #
# REPORT
# --------------------------------------------------------------------------- #
def _pct(x, nd=4):
    return f"{100.0 * float(x):.{nd}f}%"


def _sig(x, nd=4):
    return f"{float(x):.{nd}g}"


def _reg_row(rows, key, mode):
    for r in rows:
        if r["key"] == key and r["mode"] == mode:
            return r
    raise KeyError((key, mode))


def write_report(out: Dict, path: str) -> None:
    A = []
    w = A.append
    rep = out["report"]
    dp = out["drug_params"]
    tmzp = dp["temozolomide"]
    pk = out["pk_conversion"]
    ver = out["verification"]
    rows = out["regimen_rows"]
    combo = out["combination"]
    tau = tmzp["tau_tissue_s"]

    w("# REAL TEMOZOLOMIDE REGIMENS, AND COMBINATION SCHEDULING, ON ONE REAL TUMOR")
    w("")
    w("> ## WHY THIS FILE EXISTS, AND WHY IT COULD NOT HAVE EXISTED YESTERDAY")
    w(">")
    w("> Until the time-resolved plasma driver landed, `sim.engine` collapsed the "
      "whole multi-dose plasma profile to a single **time-averaged** level before "
      "the transport solver ever saw it. Two regimens with the same cumulative "
      "dose were therefore **bit-identical inputs**: the measured fractionation "
      "spread was exactly `0.00e+00`, an arithmetic artefact and not a result "
      "about any tumor. Every schedule question, and every combination-scheduling "
      "question, was **structurally unanswerable** — and the combination and "
      "regimen work done before today was comparing identical inputs.")
    w(">")
    w("> **Every solve in this report uses `PDEConfig(time_resolved=True)`.** "
      "There is no static-source path anywhere in `experiments/regimen.py`. "
      "The source is rescaled by `plasma(t)` at every timestep, so tissue "
      "concentration genuinely rises and falls and each schedule shape gets its "
      "own transport solve.")
    w(">")
    w("> **The headline, in one line.** At matched cumulative dose over a full "
      f"28-day cycle, this model predicts that the real temozolomide regimens are "
      f"indistinguishable: the total spread in mean tumour exposure across all "
      f"{len([r for r in rows if r['mode']=='matched'])} arms is "
      f"**{_pct(out['matched_exposure_spread'], 6)}**, against the analytic bound "
      f"`tau_tissue / T_window` = **{_pct(out['matched_bound_tau_over_T'], 6)}**. "
      "That is not a claim that schedule does not matter clinically. It is a "
      "claim that **transport cannot be the reason it does** — and the mechanisms "
      "that make fractionation matter (MGMT depletion, inter-fraction repair, "
      "cell-cycle redistribution, immune effects, normal-tissue recovery) are "
      "all ABSENT from this model. Section 7 says so at length.")
    w("")
    w(f"> **Configuration.** `vascular_mode = discrete` (promoted in "
      f"[`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = "
      f"{out['alpha']:.1f}** (single perivascular-sink constant, fit ONCE on "
      f"doxorubicin) · per-cell sensitivity mode = **`uniform`** (the pipeline "
      f"default and the mode `results/combination.md` used — it makes the "
      f"cross-drug comparison fair by making every agent equally an upper bound; "
      f"the engine still applies the hypoxia/quiescence resistance factor) · "
      f"**dx = {rep['dx_um']:.0f} um** · `dt = {DT_S:.0f} s` · convection on · "
      f"seed {SEED}, fully deterministic.")
    w(">")
    w("> Absolute exposures are on the **discrete-source unit** and are not "
      "comparable with KDE-era numbers elsewhere in the repo. Absolute kills over "
      "a 28-day window are also not comparable with the 6 h-window numbers in "
      "`results/dose_schedule.md` or `results/combination.md`, because the "
      "engine's Hill map consumes `AUC / T_window` — see the limitation on this "
      "in section 8. Comparisons *within* this file are apples to apples.")
    w("")
    w("> ## READ THIS FIRST")
    w(">")
    w("> * **Research prototype. NOT clinical dosing guidance.** Nothing here "
      "should influence the care of any patient. No number below is a dose "
      "recommendation, and the real regimens are modelled to make a *model* "
      "comparison interpretable, not to evaluate them.")
    w("> * **Zero outcome validation.** This 10x Xenium section is an "
      "**untreated** resection. There is no measured drug response for any agent "
      "on this tissue, so every kill fraction is a prediction with nothing to "
      "score it against.")
    w("> * **Absolute kill on this section is NOT defensible** — "
      "`results/cell_state_impact.md` establishes that. Rankings and relative "
      "comparisons are. Read every absolute number below as a unit-carrying "
      "quantity that lets the comparisons be computed, not as a prediction of "
      "how much tumor dies.")
    w(f"> * **~{out['ceiling_bands']['QUIESCENT']['fraction']*100:.0f}% of tumor "
      f"cells are labelled QUIESCENT by the oxygen-banding model** (x0.5 kill "
      f"multiplier). That single *model* labelling dominates the kill ceiling of "
      f"{out['kill_ceiling']:.4f} and every absolute kill here.")
    w("> * **Time resolution was added to TRANSPORT, not to pharmacodynamics.** "
      "The kill map is still a static Hill function of window-averaged exposure. "
      "The schedule effect reported here is a **transport-only floor**.")
    w("> * **Bliss independence is an ASSUMPTION, not a result.** The "
      "combinations model SPATIAL complementarity only. No biochemical synergy "
      "or antagonism is or can be represented.")
    w("")
    w("---")
    w("")

    # ------------------------------------------------------------------ 1
    w("## 1. Method")
    w("")
    w(f"* **Tissue — one real tumor, reconstructed ONCE.** The open 10x Xenium "
      f"GBM section: {rep['n_cells']:,} measured cells ({rep['n_tumor']:,} "
      f"tumor), grid {tuple(rep['grid'])} at the grid-converged "
      f"`dx = {rep['dx_um']:.0f} um`, mean vascular density "
      f"{rep['mean_vascular_density']:.3f}, mean ECM "
      f"{rep['mean_ecm_density']:.3f}, mean oxygen {rep['mean_oxygen']:.3f}. "
      f"Reconstruction took {out['timing']['reconstruct_s']:.1f} s and the "
      f"drug-independent prologue (oxygen -> cell-state banding -> IFP -> Darcy "
      f"velocity) {out['timing']['prologue_s']:.1f} s. **Both were run once and "
      f"reused by every arm**, so all differences below are drug/schedule and "
      f"not tissue. That reuse is the whole reason a 28-day cycle is affordable "
      f"at all; the cost it avoids is one full reconstruct + prologue "
      f"(~{out['timing']['reconstruct_s'] + out['timing']['prologue_s']:.0f} s) "
      f"per arm across {len(rows)} regimen arms plus every combination arm.")
    w(f"* **Drug compilation.** `calibrated_drugparams(name, alpha = "
      f"{out['alpha']:.2f})` then `apply_barrier_to_drugparams` — the identical "
      f"pipeline `real_drug_panel.py` and `dose_schedule.py` use.")
    w(f"* **Numerics.** `dt = {DT_S:.0f} s`, convection on. `dx = 10 um` is "
      f"deliberately not used — it violates the upwind advection CFL condition "
      f"at this `dt`.")
    w(f"* **Total cost.** {out['n_pde_solves']} transport solves, "
      f"{out['total_steps']:,} solver steps, "
      f"{out['timing']['total_s']/60:.1f} min wall clock.")
    w("")
    w("### 1.1 One solve per single dose, then EXACT superposition")
    w("")
    w("A 28-day cycle at `dt = 10 s` is 241,920 steps per shape, and there are "
      "nine regimen shapes. Solving each directly is not affordable. It is also "
      "not necessary, and the reason is exact rather than approximate.")
    w("")
    w("The discrete solve is **linear** in the plasma profile (the operator "
      "`A = I - dt*L + dt*K` does not involve the source, the scheme starts from "
      "`C = 0`, and the source enters as `dt * plasma(t_n) * u(x,y)`) and it is "
      "**time-invariant** (`D_eff`, `k` and `v` are static). `plasma_profile` "
      "builds a multi-dose curve as a superposition of shifted single-dose "
      "exponentials. Therefore, for doses of peak `P_k` at times `t_k` observed "
      "over `[0, T]`:")
    w("")
    w("```")
    w("    AUC_field  =  sum_k  P_k * A1(T - t_k)")
    w("    C(T)       =  sum_k  P_k * Cend1(T - t_k)")
    w("```")
    w("")
    w("where `A1(s)` is the exposure field from ONE unit-peak dose observed over "
      "a window of length `s`. Every dose time here is an exact multiple of `dt`, "
      "so this identity holds **exactly on the discrete grid** — it is not "
      "interpolated and not extrapolated. Amplitude is likewise exactly free, "
      "which is the same trick `dose_schedule.py` uses for its dose axis.")
    w("")
    w("`A1(s)` saturates once the plasma curve has washed out AND the tissue has "
      "cleared, so a 28-dose regimen is built from a handful of short solves. "
      "**The saturation is certified, not assumed**, with an exact bound: a unit "
      "dose eventually deposits all its drug, so")
    w("")
    w("```")
    w("    (A1(inf) - A1(S)) / A1(inf)  =  2^(-S/t_half)  +  C(S) / (k * A1(S))")
    w("                                    \\_ plasma  _/     \\_ tissue residual _/")
    w("```")
    w("")
    w("and BOTH terms are computed from the solved fields. The plasma term alone "
      "would NOT be sufficient — a molecule whose plasma has washed out may still "
      "be equilibrating in tissue, and that second term is exactly the terminal "
      "residual through which schedule acts in the first place.")
    w("")
    w("| single-dose kernel | plasma t1/2 (h) | windows solved (h) | saturation "
      "window S (h) | un-delivered plasma at S | tissue still holding at S | "
      "**certified error** |")
    w("|---|---:|---|---:|---:|---:|---:|")
    for key, km in out["kernel_meta"].items():
        name, hl = key.split("@")
        w(f"| {name} | {float(hl):.2f} | "
          f"{', '.join(f'{g:.0f}' for g in km['grid_h'])} | {km['sat_h']:.0f} | "
          f"{km['plasma_residual']:.2e} | {km['tissue_residual']:.2e} | "
          f"**{km['residual']:.2e}** |")
    w("")
    tmz_km = out["kernel_meta"][out["tmz_kernel_key"]]
    w(f"The substitution is only ever *used* for the temozolomide kernel driving "
      f"the 28-day regimen leg, where the certified error is "
      f"**{tmz_km['residual']:.1e}**. Every other leg reads exact grid entries; "
      f"the kernel refuses to interpolate and raises instead.")
    w("")
    w("### 1.2 Verification — three independent checks, each a full extra solve")
    w("")
    w("The superposition is not taken on faith.")
    w("")
    w("| # | check | statistic | result |")
    w("|---|---|---|---|")
    w(f"| V1 | {ver['v1']['label']} | relative L1 over tissue | "
      f"**{ver['v1']['rel_l1']:.2e}** |")
    w(f"| V2 | {ver['v2']['label']} | kill fraction, engine vs fast path | "
      f"engine {ver['v2']['engine_kill']:.10f} vs "
      f"{ver['v2']['fast_kill']:.10f}, delta "
      f"**{ver['v2']['abs_delta']:.1e}** |")
    w(f"| V2b | (same run) | exposure field, relative L1 | "
      f"**{ver['v2']['auc_rel_l1']:.2e}** |")
    w(f"| V3 | {ver['v3']['label']} | relative L1 over tissue | "
      f"**{ver['v3']['rel_l1']:.2e}** |")
    w("")
    w("V1 tests LINEARITY (a genuine two-dose profile against the sum of two "
      "single-dose kernels). V3 tests TIME-INVARIANCE on its own, by driving the "
      "solver with a manually shifted profile — zeros for 12 h, then the dose — "
      "and comparing against the un-shifted kernel over the remaining window. "
      "V2 puts the **unmodified** `sim.engine.simulate_response` with "
      "`time_resolved=True` against the fast path end to end, so nothing about "
      "the orchestration here is unpinned.")
    w("")
    w("### 1.3 Drug parameters as compiled onto this tissue")
    w("")
    w("| drug | MW (Da) | D (um^2/s) | IC50 (uM) | hill | BBB | k_loss (1/s) | "
      "tissue time-scale tau = 1/k |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|")
    for name in COMBO_DRUGS:
        d = dp[name]
        w(f"| {name} | {d['mw']:.1f} | {d['diffusivity']:.0f} | "
          f"{d['ic50']:.1f} | {d['hill']:.1f} | {d['bbb']:.3f} | "
          f"{d['k_loss']:.2e} | {1.0/d['k_loss']:.0f} s |")
    w("")
    w(f"**Kill ceiling on this tumor = {out['kill_ceiling']:.4f}** — the kill "
      f"fraction reached as exposure -> infinity, where the Hill term saturates "
      f"and per-cell kill probability tends to `sensitivity x resistance`. "
      f"**No dose and no schedule can exceed it.**")
    w("")
    w("| tumor-cell band | fraction of tumor cells | mean asymptotic kill "
      "probability | contribution to ceiling |")
    w("|---|---:|---:|---:|")
    for band in ("VIABLE", "HYPOXIC", "QUIESCENT"):
        b = out["ceiling_bands"][band]
        w(f"| {band} | {b['fraction']*100:.1f}% | {b['mean_p']:.3f} | "
          f"{b['contribution']:.4f} |")
    w("")
    w("Temozolomide's per-cell sensitivity determinant (MGMT) is **not on this "
      "panel**, so under `uniform` sensitivity the ceiling is *entirely* the "
      "hypoxia/quiescence multipliers. Temozolomide resistance is therefore "
      "**unmodelled** and its absolute kill is correspondingly optimistic.")
    w("")
    w("---")
    w("")

    # ------------------------------------------------------------------ 2
    w("## 2. The real regimens, and the PK conversion")
    w("")
    w("### 2.1 What is actually given to patients")
    w("")
    w("Every regimen below is a real, published schedule. The 28-day comparison "
      "window is one adjuvant cycle; the Stupp concomitant phase runs 42 days in "
      "practice and only its first 28 days are modelled, so its cumulative "
      "figure here is per 28 days and not per course.")
    w("")
    w("| regimen | dose | dosing days of the 28-day window | doses | "
      "cumulative mg/m2 per 28 d | vs Stupp adjuvant C1 | source |")
    w("|---|---|---|---:|---:|---:|---|")
    cites = {
        "stupp_adj_150": "Stupp 2005 (PMID 15758009)",
        "stupp_adj_200": "Stupp 2005 (PMID 15758009)",
        "concomitant_75": "Stupp 2005 (PMID 15758009)",
        "dd_21_28_75": "Gilbert 2013, RTOG 0525 (PMID 24101040)",
        "dd_21_28_100": "Gilbert 2013, RTOG 0525 (PMID 24101040)",
        "metronomic_7_7": "Wick 2004 (PMID 15184628); Wick 2007 (PMID 17664483)",
        "director_a": "Weller 2015, DIRECTOR (PMID 25655102)",
        "director_b": "Weller 2015, DIRECTOR (PMID 25655102)",
        "continuous_50": "Perry 2010, RESCUE (PMID 20308655)",
        "single_bolus": "not a real regimen — CONTROL",
    }
    daydesc = {
        "stupp_adj_150": "1-5", "stupp_adj_200": "1-5",
        "concomitant_75": "1-28 (continuous)",
        "dd_21_28_75": "1-21", "dd_21_28_100": "1-21",
        "metronomic_7_7": "1-7, 15-21",
        "director_a": "1-7, 15-21", "director_b": "1-21",
        "continuous_50": "1-28 (continuous)", "single_bolus": "1",
    }
    ref_cum = 750.0
    for reg in REGIMENS:
        r = _reg_row(rows, reg["key"], "prescribed")
        unit = "mg/m2/d" if r["n_doses"] > 1 else "mg/m2, once"
        w(f"| {reg['label']} | {reg['mg_m2']:.0f} {unit} | "
          f"{daydesc[reg['key']]} | {r['n_doses']} | "
          f"**{r['cumulative_mg_m2_prescribed']:.0f}** | "
          f"{r['cumulative_mg_m2_prescribed']/ref_cum:.2f}x | "
          f"{cites[reg['key']]} |")
    w("")
    w("Two facts from that table are worth stating before any simulation output:")
    w("")
    w("1. **The dose-dense and metronomic regimens are not schedule variants of "
      "the standard one — they deliver 2-3x the cumulative dose.** RTOG 0525's "
      "dose-dense arm delivered ~2.1x the standard arm's cumulative "
      "temozolomide. Any comparison that does not hold cumulative dose fixed is "
      "measuring dose, not schedule. Section 3 shows exactly that, and section 4 "
      "then holds dose fixed.")
    w("2. **The two DIRECTOR arms ARE cumulative-dose-matched by design** "
      "(120 mg/m2 x 14 d = 80 mg/m2 x 21 d = 1,680 mg/m2 per cycle). That trial "
      "is the closest real-world analogue of the matched-dose comparison in "
      "section 4, and it was **negative** (median TTF 1.8 vs 2.0 months, OS 9.8 "
      "vs 10.6 months) — though underpowered, having closed early at n=105.")
    w("")
    w("### 2.2 The PK conversion, and what it buries")
    w("")
    w("mg/m2 has to become a plasma curve the transport solver can consume. The "
      "conversion is a single arithmetic chain, stated in full so it can be "
      "checked or rejected:")
    w("")
    w("| step | value | source |")
    w("|---|---|---|")
    w("| Oral bioavailability | treated as **1.0** | FDA label: TEMODAR IV over "
      "90 min is **bioequivalent** to 150 mg/m2 oral for both Cmax and AUC "
      "(7.5 vs 7.3 ug/mL; 23.4 vs 24.6 ug*h/mL) -> F ~ 0.95-1.0. Baker 1999 "
      "(PMID 10037179): *\"rapid oral absorption and high systemic "
      "availability\"* |")
    w(f"| Plasma elimination half-life | **{TMZ_HALF_LIFE_H} h** | FDA label "
      "(*\"mean elimination half-life of 1.8 hours\"*); Beale 1999 (PMID "
      "10501912) 1.8 h; Baker 1999 1.9 h; Ostermann 2004 (PMID 15173079) 2.1 h "
      "population estimate. Verified range 1.8-2.1 h |")
    w(f"| Molar mass | **{TMZ_MW} g/mol** | PubChem CID 5394, C6H6N6O2 |")
    w(f"| Anchor: plasma AUC(0-inf) after {TMZ_ANCHOR_DOSE_MG_M2:.0f} mg/m2 | "
      f"**{TMZ_ANCHOR_AUC_UG_H_ML} ug*h/mL** = "
      f"{TMZ_ANCHOR_AUC_UG_H_ML*1000/TMZ_MW:.1f} uM*h | FDA label, geometric "
      "mean |")
    w(f"| Dose linearity | assumed linear over "
      f"{TMZ_LINEAR_RANGE_MG_M2[0]:.0f}-{TMZ_LINEAR_RANGE_MG_M2[1]:.0f} mg/m2/d "
      "| FDA label: *\"linear kinetics over the therapeutic dosing range of 75 "
      "mg/m2/day to 250 mg/m2/day\"* |")
    w(f"| Model curve | 1-compartment, instantaneous absorption, "
      f"mono-exponential elimination -> `AUC = peak / lambda` | "
      f"`tmesim.treatment.drugs.plasma_profile` |")
    w(f"| **Resulting peak at {TMZ_ANCHOR_DOSE_MG_M2:.0f} mg/m2** | "
      f"**{pk['peak_per_150']:.1f} uM** | `peak = AUC x lambda` |")
    w("")
    w(f"**The conversion is anchored on AUC, not on Cmax, and the discrepancy is "
      f"reported rather than hidden.** Tissue exposure in this model is driven by "
      f"the plasma integral, so AUC is the right invariant to preserve. The price "
      f"is that the model's implied peak, **{pk['model_peak_um']:.1f} uM**, "
      f"exceeds the measured geometric-mean Cmax of "
      f"**{pk['measured_cmax_um']:.1f} uM** "
      f"({TMZ_ANCHOR_CMAX_UG_ML} ug/mL) by a factor of "
      f"**{pk['ratio']:.2f}**. That is the cost of instantaneous absorption: "
      f"real oral temozolomide has a median tmax of ~1 h, so the model "
      f"front-loads a dose that really arrives over about an hour.")
    w("")
    w("**Direction of that error, explicitly.** It makes every modelled dose "
      "*more* bolus-like than reality. `results/dose_schedule.md` establishes "
      "that in this transport operator the bolus always wins, so the bias favours "
      "the more concentrated regimens. It is identical across regimens at matched "
      "cumulative dose, so it **cancels from the comparison** and survives only "
      "in the absolute levels.")
    w("")
    w("**Other assumptions this conversion buries, named:**")
    w("")
    w("* **No active-species split.** Temozolomide's cytotoxic species is MTIC, "
      "formed by pH-dependent hydrolysis; MTIC AUC is 2-4% of parent AUC (Beale "
      "1999). The model carries the parent curve only.")
    w("* **No brain-compartment step.** Ostermann 2004 measured "
      "`AUC_CSF / AUC_plasma = 20%` in n=35 glioma patients. That partitioning "
      "is *not* applied on top of the plasma curve here — the model's own "
      "`bbb_penetrance` term "
      f"({tmzp['bbb']:.3f}) plays that role, and the two are not the same object. "
      "Do not read the plasma peaks below as tumour concentrations.")
    w("* **No plasma protein binding correction** (~15% bound per the label).")
    w("* **Population-average mg/m2 -> plasma mapping**, not patient-specific.")
    w("* **50 mg/m2 (RESCUE) and the rescaled matched-dose arms fall outside the "
      "label's stated linear range of 75-250 mg/m2/day.** Rows where that "
      "happens are flagged in the tables below.")
    w("")
    w("---")
    w("")

    # ------------------------------------------------------------------ 3
    w("## 3. AS PRESCRIBED — this is a DOSE comparison, not a schedule comparison")
    w("")
    w("Each regimen at its own real dose, over one 28-day window. Cumulative dose "
      "differs by 2.8x across the table, so this is a dose axis with a schedule "
      "axis confounded into it. It is reported first precisely so that section 4 "
      "can be read as the schedule question with the dose confound removed.")
    w("")
    w("| regimen | cumulative mg/m2 / 28 d | peak per dose (uM) | "
      "mean tumour AUC (uM*s) | E[kill] | kill (1 draw) | % of ceiling | "
      "tumor cells >= IC50 | median E / IC50 |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for reg in REGIMENS:
        r = _reg_row(rows, reg["key"], "prescribed")
        flag = " *" if r["outside_linear_range"] else ""
        w(f"| {reg['label']}{flag} | {r['cumulative_mg_m2']:.0f} | "
          f"{r['peak_um']:.1f} | {r['mean_auc_tumor']:.4g} | "
          f"{r['kill_expected']:.6f} | {r['kill_drawn']:.6f} | "
          f"{100*r['kill_expected']/out['kill_ceiling']:.2f}% | "
          f"{r['ic50_coverage']*100:.2f}% | {r['median_E_over_ic50']:.5f} |")
    w("")
    w("`*` = per-dose mg/m2 outside the label's stated linear-PK range "
      "(75-250 mg/m2/day); the linear conversion is an extrapolation for that row.")
    w("")
    pres_all = [_reg_row(rows, r["key"], "prescribed") for r in REGIMENS]
    pres = sorted(pres_all, key=lambda r: -r["kill_expected"])
    by_dose = sorted(pres_all, key=lambda r: r["cumulative_mg_m2"])
    monotone = all(by_dose[i]["kill_expected"] <= by_dose[i + 1]["kill_expected"]
                   + 1e-12 for i in range(len(by_dose) - 1))
    verdict = ("**Ranking is monotone in cumulative dose, and in nothing else** "
               "(checked, not asserted: sorting the table by cumulative mg/m2 "
               "gives a non-decreasing E[kill] at every step)."
               if monotone else
               "**Ranking is ALMOST but not exactly monotone in cumulative dose** "
               "— the exceptions are where the schedule effect of section 4 "
               "becomes visible.")
    w(f"{verdict} The top arm ({pres[0]['label']}, "
      f"{pres[0]['cumulative_mg_m2']:.0f} mg/m2) reaches E[kill] "
      f"{pres[0]['kill_expected']:.6f}; the bottom ({pres[-1]['label']}, "
      f"{pres[-1]['cumulative_mg_m2']:.0f} mg/m2) reaches "
      f"{pres[-1]['kill_expected']:.6f}, a "
      f"{pres[0]['kill_expected']/max(pres[-1]['kill_expected'],1e-30):.2f}x "
      "range across a 2.8x range in cumulative dose. That is the entire content "
      "of this table: in a transport model with a static Hill kill map, "
      "cumulative dose is the only lever visible at this scale.")
    w("")
    w("> **This is where the model and the clinic diverge, and it is worth being "
      "blunt about it.** RTOG 0525's dose-dense arm delivered ~2.1x the "
      "cumulative dose of the standard arm and produced **no overall-survival "
      "benefit** (median OS 14.9 vs 16.6 months, HR 1.03, P = .63) while doubling "
      "grade >= 3 toxicity (53% vs 34%). A model in which outcome is a monotone "
      "function of delivered exposure — which is exactly what this model is — "
      "**cannot reproduce that trial**. Section 7 lists what would have to be "
      "added before it could.")
    w("")
    max_cov = max(r["ic50_coverage"] for r in pres_all)
    w(f"Note also that IC50 coverage never exceeds {max_cov*100:.2f}% of tumor "
      f"cells on any arm. Temozolomide's cached GDSC glioma IC50 is "
      f"{tmzp['ic50']:.0f} uM against a median window-averaged tumour exposure "
      f"of {pres[0]['median_E_over_ic50']*tmzp['ic50']:.4g} uM even on the "
      "highest-dose arm. This is consistent with `results/dose_schedule.md`, "
      "which finds temozolomide dose-limited across its whole swept range and "
      "never reaching IC50 at clinical plasma levels.")
    w("")
    w("---")
    w("")

    # ------------------------------------------------------------------ 4
    w("## 4. AT MATCHED CUMULATIVE DOSE — schedule as the ONLY variable")
    w("")
    w(f"Every regimen is now rescaled to deliver exactly "
      f"**{MATCHED_CUMULATIVE_MG_M2:.0f} mg/m2 over the 28-day window** (the "
      f"standard Stupp adjuvant cycle at 200 mg/m2 x 5 days). Only the *timing* "
      f"of that fixed total differs. This is the comparison the question actually "
      f"asks for.")
    w("")
    mrows = [_reg_row(rows, r["key"], "matched") for r in REGIMENS]
    ref_plasma = max(r["plasma_auc_um_s"] for r in mrows)
    eff = [r["mean_auc_tumor"] / r["plasma_auc_um_s"] for r in mrows]
    ref_eff = max(eff)
    w("| regimen | doses | mg/m2 per dose | peak per dose (uM) | plasma AUC "
      "COUNTED in window (uM*s) | vs best | **mean tumour AUC (uM*s)** | vs best "
      "| tissue AUC per counted plasma AUC | vs best | E[kill] | tumor >= IC50 | "
      "terminal residual C(T) |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for reg, e in zip(REGIMENS, eff):
        r = _reg_row(rows, reg["key"], "matched")
        flag = " *" if r["outside_linear_range"] else ""
        w(f"| {reg['label']}{flag} | {r['n_doses']} | "
          f"{r['mg_m2_per_dose']:.1f} | {r['peak_um']:.2f} | "
          f"{r['plasma_auc_um_s']:.10g} | "
          f"{(r['plasma_auc_um_s']-ref_plasma)/ref_plasma*100:+.6f}% | "
          f"**{r['mean_auc_tumor']:.10g}** | "
          f"{r['auc_rel_to_best']*100:+.6f}% | {e:.10g} | "
          f"{(e-ref_eff)/ref_eff*100:+.6f}% | {r['kill_expected']:.8f} | "
          f"{r['ic50_coverage']*100:.2f}% | {r['residual']*100:.4f}% |")
    w("")
    w(f"**Total spread in mean tumour exposure: "
      f"{_pct(out['matched_exposure_spread'], 6)}.** The analytic bound is "
      f"`tau_tissue / T_window` = {tau:.0f} s / 28 d = "
      f"{_pct(out['matched_bound_tau_over_T'], 6)}, and the measured spread sits "
      f"under it, as it must.")
    w("")
    w("**Reading the three 'vs best' columns — the difference decomposes "
      "cleanly.** Matching *administered* cumulative mg/m2 is the clinically "
      "meaningful match, but it is not quite the same thing as matching the drug "
      "the solver actually integrates. Two distinct terminal effects separate "
      "them, and both columns are here so neither is hidden:")
    w("")
    w("1. **Drug not yet delivered when the window closes** — a dose given on day "
      "28 has not finished being absorbed and eliminated by day 28, so part of "
      "its plasma curve falls outside the integration. That is the *plasma AUC "
      "counted in window* column, and it is a **measurement-window artefact**, "
      "not a schedule effect.")
    w("2. **Drug delivered but still resident in tissue at T** — the terminal "
      "residual `C(T)`, the last column. This one IS the schedule effect: it is "
      "the `C(T)` in `A = M^-1 (S_total - C(T))`.")
    w("")
    w("The **tissue AUC per counted plasma AUC** column divides the first effect "
      "out and leaves only the second, so it is the cleanest possible statement "
      "of the pure schedule effect at matched dose. Its spread across all arms "
      f"is **{_pct((ref_eff - min(eff)) / ref_eff, 6)}**.")
    w("")
    w("### 4.1 Why the spread is where it is, exactly")
    w("")
    w("`results/dose_schedule.md` section 4b derives it. Integrating the "
      "transport PDE over `[0, T]` with `C(0) = 0`, for `A = AUC` and "
      "`M = k - div(D grad .) + div(v .)`:")
    w("")
    w("```")
    w("    A = M^-1 ( S_total - C(T) )")
    w("```")
    w("")
    w("`S_total` is **identical** for any two schedules at matched cumulative "
      "exposure — that is what matched means, and the implementation guarantees "
      "it exactly (the plasma curve is sampled on the solver's own step grid). "
      "So the entire schedule effect is the **terminal residual `C(T)`**: drug "
      "delivered but still in flight when the window closes. Read the last two "
      "columns of the table above together — every regimen whose exposure differs "
      "from the best is a regimen with a non-negligible terminal residual, and "
      "they are precisely the ones that are still dosing near day 28.")
    w("")
    w("Concretely: the 5/28, 21/28 and 7/7 regimens all finish dosing at least "
      "seven days before the window closes, so their drug has washed out "
      "completely and they realise the **full** steady-state exposure field — "
      "they are identical to each other to the last stored digit. Only the "
      "continuous daily regimens, which dose on day 28 itself, lose anything, "
      "and what they lose is one day's terminal residual out of 28 days of dosing.")
    w("")
    w("### 4.2 The same comparison up the dose axis")
    w("")
    w("The schedule effect is a fixed **relative** change in exposure, so the "
      "**absolute** kill spread grows with where you sit on the Hill curve. At "
      "clinical amplitude temozolomide's predicted kill on this section is tiny "
      "and a sub-0.02% relative change lands far below the "
      f"1/{out['n_tumor']:,} = {1.0/out['n_tumor']:.2e} quantum of a "
      "cell-counted kill fraction. Rescaling up the dose axis is exactly free "
      "(exposure is linear in amplitude) and lifts it clear. **The multipliers "
      "below are a diagnostic amplitude, NOT a dose — 1000x clinical is not a "
      "regimen anyone could give.** Expected kill (mean per-cell probability) is "
      "reported because it carries no cell-counting quantisation.")
    w("")
    keys_show = ["stupp_adj_200", "dd_21_28_100", "metronomic_7_7",
                 "concomitant_75", "single_bolus"]
    hdr = " | ".join(next(r["label"] for r in REGIMENS if r["key"] == k)
                     for k in keys_show)
    w(f"| amplitude (x clinical) | {hdr} | spread |")
    w("|---:|" + "---:|" * (len(keys_show) + 1))
    for row in out["regimen_amplitude_rows"]:
        vals = [row[k]["kill_expected"] for k in keys_show]
        allv = [row[k]["kill_expected"] for k in row if k != "mult"]
        w(f"| {row['mult']:.0f}x | "
          + " | ".join(f"{v:.8f}" for v in vals)
          + f" | **{max(allv)-min(allv):.2e}** |")
    w("")
    w("Every one of those spreads would have been **exactly `0.00e+00`** under "
      "the static source, for every amplitude, because the inputs were "
      "bit-identical. They are now small, non-zero and bounded — which is a "
      "different statement from the one the old engine forced.")
    w("")

    # 4.3 window ladder
    w("### 4.3 The bound is a ratio, so it moves with the observation window")
    w("")
    w("`tau_tissue / T_window` is the whole story, which means the size of the "
      "measured schedule effect is a **reporting choice** as much as a physical "
      "one. The same burst-vs-continuous contrast, at matched cumulative dose, "
      "re-measured over three windows:")
    w("")
    w("| window | burst 5/28 mean tumour AUC | continuous daily mean tumour AUC | "
      "single bolus | measured spread | bound tau/T |")
    w("|---|---:|---:|---:|---:|---:|")
    for e in out["window_ladder"]:
        a = e["arms"]
        w(f"| {e['window_s']/DAY:.0f} d | "
          f"{a['burst 5/28 (Stupp adjuvant)']['mean_auc_tumor']:.8g} | "
          f"{a['continuous daily']['mean_auc_tumor']:.8g} | "
          f"{a['single bolus']['mean_auc_tumor']:.8g} | "
          f"**{_pct(e['exposure_spread_rel'], 6)}** | "
          f"{_pct(e['bound_tau_over_T'], 6)} |")
    w("")
    w("The measured spread halves when the window doubles, exactly as `1/T` "
      "predicts. **Quote any schedule percentage from this model with its window "
      "attached, or it means nothing.**")
    w("")

    # 4.4 intraday
    idy = out["intraday"]
    w("### 4.4 Where the effect IS resolvable: intra-day fractionation")
    w("")
    w(f"To show the machinery is measuring physics rather than round-off, the "
      f"same matched-cumulative-exposure comparison is run at a "
      f"{idy['window_s']/HOUR:.0f} h window, where `tau/T` = "
      f"{_pct(idy['bound_tau_over_T'], 4)} and the effect is comfortably above "
      f"any numerical floor. One day's dose given as 1, 2 or 4 fractions:")
    w("")
    w("| schedule | peak per dose (uM) | plasma AUC in window (uM*s) | "
      "mean tumour AUC (uM*s) | E[kill] | tumor >= IC50 | terminal residual |")
    w("|---|---:|---:|---:|---:|---:|---:|")
    for label, _ in INTRADAY_SCHEDULES:
        a = idy["arms"][label]
        w(f"| {label} | {a['peak_um']:.2f} | {a['plasma_auc_um_s']:.6g} | "
          f"{a['mean_auc_tumor']:.8g} | {a['kill_expected']:.6f} | "
          f"{a['ic50_coverage']*100:.2f}% | {a['residual']*100:.3f}% |")
    w("")
    w(f"**Exposure spread {_pct(idy['exposure_spread_rel'], 4)}** against the "
      f"bound {_pct(idy['bound_tau_over_T'], 4)}, and it is **monotone: the "
      f"bolus wins.** That is the same sign `results/dose_schedule.md` found "
      "across all 48 regimes it tested, and the reason is structural — the "
      "earlier a schedule finishes, the more of the steady-state exposure field "
      "it realises inside the window, and there is no saturable perivascular "
      "sink in this transport operator to reverse it. The terminal-residual "
      "column shows the mechanism directly: it rises with the number of "
      "fractions.")
    w("")
    w("| amplitude | " + " | ".join(l for l, _ in INTRADAY_SCHEDULES)
      + " | spread (expected kill) |")
    w("|---:|" + "---:|" * (len(INTRADAY_SCHEDULES) + 1))
    for row in idy["amplitude"]:
        w(f"| {row['mult']:.0f}x | "
          + " | ".join(f"{row['kill'][l]['expected']:.6f}"
                       for l, _ in INTRADAY_SCHEDULES)
          + f" | **{row['spread_expected']:.2e}** |")
    w("")
    w("---")
    w("")

    # ------------------------------------------------------------------ 5
    w("## 5. COMBINATIONS, re-run with time resolution")
    w("")
    w("> ### THE COMBINATION ASSUMPTION: BLISS INDEPENDENCE — read this before "
      "any number")
    w(">")
    w("> Each agent is simulated **separately** on the same real tissue, giving "
      "every tumor cell its own exposure AUC and its own kill probability from "
      "the engine's own Hill pharmacodynamics. The combination is formed at the "
      "exposure level by")
    w(">")
    w("> ```")
    w("> p_combo(cell) = 1 - (1 - p_A(cell)) * (1 - p_B(cell))")
    w("> ```")
    w(">")
    w("> **NO pharmacological synergy or antagonism is modelled.** There is no "
      "shared target, no DNA-repair crosstalk, no efflux competition, no "
      "drug-drug PK interaction, and drug B does not alter drug A's transport. "
      "What is modelled is **SPATIAL complementarity only** — whether the two "
      "delivery fields reach *different cells*. Every gain below is geographic "
      "or arithmetic, never biochemical. Real combinations may synergise or "
      "antagonise and this model is blind to both.")
    w(">")
    w("> Because Bliss independence *is* the statement that the per-cell draws "
      "are independent, all set statistics are exact expectations under that "
      "coupling (`E|A n B| = sum p_A p_B`, `E|A u B| = sum [1-(1-p_A)(1-p_B)]`) "
      "and carry **no Monte-Carlo noise**.")
    w("")
    w("### 5.1 Single agents")
    w("")
    conv = combo["conventions"]
    _hls = [v["half_life_h"] for v in conv["clinical"].values()]
    hl_span = max(_hls) / min(_hls)
    w(f"Two PK conventions are run, because **plasma half-life is the parameter "
      f"a stagger acts on** and answering the staggering question only on a "
      f"shared placeholder would be evasive:")
    w("")
    w("| convention | temozolomide | cisplatin | doxorubicin | rationale |")
    w("|---|---|---|---|---|")
    w(f"| `shared` | {conv['shared']['temozolomide']['peak_um']:.1f} uM, "
      f"{conv['shared']['temozolomide']['half_life_h']:.1f} h | "
      f"{conv['shared']['cisplatin']['peak_um']:.1f} uM, "
      f"{conv['shared']['cisplatin']['half_life_h']:.1f} h | "
      f"{conv['shared']['doxorubicin']['peak_um']:.1f} uM, "
      f"{conv['shared']['doxorubicin']['half_life_h']:.1f} h | identical to "
      "`results/combination.md`, so the ONLY differences from that report are "
      "time resolution and staggering |")
    w(f"| `clinical` | {conv['clinical']['temozolomide']['peak_um']:.1f} uM, "
      f"{conv['clinical']['temozolomide']['half_life_h']:.1f} h | "
      f"{conv['clinical']['cisplatin']['peak_um']:.1f} uM, "
      f"{conv['clinical']['cisplatin']['half_life_h']:.1f} h | "
      f"{conv['clinical']['doxorubicin']['peak_um']:.1f} uM, "
      f"{conv['clinical']['doxorubicin']['half_life_h']:.1f} h | each agent on "
      f"its own measured plasma half-life — a **{hl_span:.0f}x span**, "
      "which is exactly the axis staggering acts on |")
    w("")
    w("**Provenance of the clinical PK, with its weaknesses stated:**")
    w("")
    w(f"* **temozolomide** — {conv['clinical']['temozolomide']['half_life_h']} h "
      f"half-life (FDA label; Beale 1999, PMID 10501912), peak from the section "
      f"2.2 conversion at {CLINICAL_TMZ_DOSE_MG_M2:.0f} mg/m2. Solid.")
    w(f"* **cisplatin** — {conv['clinical']['cisplatin']['half_life_h']} h for "
      "**free / ultrafilterable** platinum (Himmelstein 1981, PMID 7194166: "
      "*\"t1/2 of 0.3 to 0.5 hr\"*; total platinum is triphasic with a terminal "
      "t1/2 >= 24 h and is the wrong species here). The peak is **DERIVED, not "
      "measured**: unchanged cisplatin Cmax 1.5-2 ug/mL after 80 mg/m2 (Nagai "
      "1996, PMID 8995510) divided by Himmelstein's cisplatin/filterable-Pt "
      "ratio of 0.6-0.8 gives ~6-11 uM; 8 uM is the midpoint. Flag this. Also "
      "note that free-Pt elimination is dominated by irreversible covalent "
      "protein binding, which a single first-order clearance term mis-specifies.")
    w(f"* **doxorubicin** — {conv['clinical']['doxorubicin']['half_life_h']} h "
      "terminal half-life (Benjamin 1977, PMID 856462, *\"triphasic ... mean "
      "half-life of 30 hr\"*; Camaggi 1988 gives t1/2-gamma 48.4 h; interpatient "
      "CV ~65%). Peak 2.8 uM is the measured 1509 ng/mL at 35 mg/m2 IV bolus "
      "(Muller 1993, PMID 8339389). **The peak is dominated by sampling time** — "
      "doxorubicin's distribution half-life is ~4 min, so any reported Cmax is "
      "as much a statement about when blood was drawn as about the drug. This "
      "is the least trustworthy number in the report.")
    w("")
    w(f"Single agents over the {combo['window_s']/HOUR:.0f} h combination window, "
      f"dosed at t = 0:")
    w("")
    for cname in ("shared", "clinical"):
        w(f"**`{cname}` convention**")
        w("")
        w("| drug | peak (uM) | plasma t1/2 (h) | plasma AUC in window (uM*s) | "
          "mean tumour AUC (uM*s) | E[kill] | tumor >= IC50 | median E / IC50 | "
          "terminal residual |")
        w("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for name in COMBO_DRUGS:
            s = combo[f"singles_24h_{cname}"][f"{name}|t0"]
            w(f"| {name} | {s['peak_um']:.2f} | {s['half_life_h']:.1f} | "
              f"{s['plasma_auc_um_s']:.4g} | {s['mean_auc_tumor']:.4g} | "
              f"{s['kill_expected']:.6f} | {s['ic50_coverage']*100:.2f}% | "
              f"{s['median_E_over_ic50']:.4f} | {s['residual']*100:.3f}% |")
        w("")

    # 5.2 pair table + dose-matched controls
    w("### 5.2 Pairs, with the dose-matched controls that decide the question")
    w("")
    w("Two controls, because \"the combination kills more\" is not a result on "
      "its own:")
    w("")
    w("* **2x dose of the better single agent** — the same total drug, in one "
      "agent. If the pair cannot beat this, it is not worth two drugs.")
    w("* **Bliss(better + ITSELF)**, i.e. `1 - (1 - p_best)^2` — a hypothetical "
      "arm whose two components have *perfectly redundant* delivery. This is the "
      "control that separates spatial complementarity from arithmetic: under a "
      "saturating Hill curve, **any** second independent hit beats **any** second "
      "dose, including the same drug twice. If `combo(A,B) <= selfcombo(B,B)`, "
      "the apparent benefit is the independence assumption meeting a saturating "
      "curve, not geography.")
    w("")
    for cname in ("shared", "clinical"):
        w(f"**`{cname}` convention, simultaneous dosing, "
          f"{combo['window_s']/HOUR:.0f} h window**")
        w("")
        w("| pair | A alone | B alone | A+B (Bliss) | better | 2x better | "
          "Bliss(better+ITSELF) | combo / better | combo / 2x-dose | "
          "combo / self-combo |")
        w("|---|---:|---:|---:|---|---:|---:|---:|---:|---:|")
        for e in combo[f"pairs_{cname}"]:
            s = e["arms"]["simultaneous"]
            best = max(s["kill_a"], s["kill_b"])
            w(f"| {e['a']} + {e['b']} | {s['kill_a']:.4f} | {s['kill_b']:.4f} | "
              f"{s['kill_combo']:.4f} | {e['better_24h']} | "
              f"{e['kill_2x_better_24h']:.4f} | {e['kill_selfcombo_24h']:.4f} | "
              f"{s['kill_combo']/best:.2f}x | "
              f"{s['kill_combo']/e['kill_2x_better_24h']:.2f}x | "
              f"**{s['kill_combo']/e['kill_selfcombo_24h']:.2f}x** |")
        w("")
    w("The last column is the one that matters. A ratio at or below 1.00 means "
      "combining two *different* drugs did no better than the model's own "
      "independence arithmetic applied to a single drug against itself.")
    w("")

    # 5.3 stagger
    w("### 5.3 STAGGERED vs SIMULTANEOUS — the genuinely new question")
    w("")
    w(f"The old engine could not pose this at all: with the plasma curve "
      f"averaged away, \"drug B at t = 0\" and \"drug B at t = "
      f"{combo['stagger_s']/HOUR:.0f} h\" were the same input array. Four arms "
      f"over a {combo['window_s']/HOUR:.0f} h window:")
    w("")
    w("* **simultaneous** — both agents at `t = 0`.")
    w(f"* **B delayed, matched administered dose** — same peak, B given at "
      f"`t = {combo['stagger_s']/HOUR:.0f} h`. B then has less window left in "
      "which to be delivered *and counted*, so it loses exposure. That is a "
      "**measurement-window artefact**, the same one `results/dose_schedule.md` "
      "names as its 'leg B'. It is reported because it is what a naive "
      "\"same milligrams\" comparison would produce, and it should not be "
      "mistaken for a schedule effect.")
    w("* **B delayed, matched COUNTED exposure** — B's peak is rescaled so its "
      "plasma AUC *inside the window* matches the simultaneous arm. `S_total` is "
      "then identical and the only remaining difference is the terminal "
      "residual. **This is the leg that isolates schedule.**")
    w("* **A delayed, matched counted exposure** — the reverse order, because "
      "which agent goes first is a real clinical choice.")
    w("")
    w("Jaccard is the overlap of the killed populations, "
      "`E|A n B| / E|A u B|`. Because independent-Bernoulli overlap shrinks "
      "automatically as kill probabilities shrink, the raw Jaccard proves nothing "
      "on its own; the **delivery-shuffle null** "
      f"({N_SHUFFLE} permutations, {N_STRATA} strata) permutes only drug B's "
      "delivery term within strata of the shared cell-intrinsic killability "
      "factor, so each cell keeps its own intrinsic killability and only B's "
      "delivery *geography* is scrambled. `obs/null` < 1 = complementary "
      "territory, ~1 = no spatial structure, > 1 = redundant.")
    w("")
    for cname in ("shared", "clinical"):
        w(f"**`{cname}` convention**")
        w("")
        w("| pair | arm | A kill | B kill | combo kill | **Jaccard** | "
          "delivery-shuffle null | obs/null | both (% of killed) | "
          "A only | B only |")
        w("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for e in combo[f"pairs_{cname}"]:
            for label in ("simultaneous",
                          "B delayed, matched administered",
                          "B delayed, matched counted exposure",
                          "A delayed, matched counted exposure"):
                s = e["arms"][label]
                w(f"| {e['a']} + {e['b']} | {label} | {s['kill_a']:.4f} | "
                  f"{s['kill_b']:.4f} | {s['kill_combo']:.4f} | "
                  f"**{s['jaccard']:.6f}** | {s['null_jaccard']:.6f} | "
                  f"{s['obs_over_null']:.3f}x | {s['frac_both']*100:.2f}% | "
                  f"{s['frac_a_only']*100:.2f}% | {s['frac_b_only']*100:.2f}% |")
        w("")
    w("**The answer, stated plainly.** Staggering changes the overlap of the "
      "killed populations by:")
    w("")
    w("| pair | convention | Jaccard, simultaneous | Jaccard, B delayed "
      "(matched counted) | absolute change | relative change |")
    w("|---|---|---:|---:|---:|---:|")
    stag_rel = []
    for cname in ("shared", "clinical"):
        for e in combo[f"pairs_{cname}"]:
            j0 = e["arms"]["simultaneous"]["jaccard"]
            j1 = e["arms"]["B delayed, matched counted exposure"]["jaccard"]
            rel = (j1 - j0) / j0 if j0 else float("nan")
            stag_rel.append(abs(rel))
            w(f"| {e['a']} + {e['b']} | `{cname}` | {j0:.6f} | {j1:.6f} | "
              f"{j1-j0:+.2e} | **{rel*100:+.3f}%** |")
    w("")
    w(f"Largest relative change in Jaccard from staggering, across all six "
      f"pair-by-convention cells: **{max(stag_rel)*100:.3f}%**.")
    w("")
    w("> **And this is a STRUCTURAL result, not an empirical one — which is why "
      "it should be believed and why it is also nearly useless clinically.** In "
      "this model the transport operator contains **no drug-drug coupling of any "
      "kind**: each agent's exposure field is solved independently, and Bliss "
      "independence then combines the two per-cell kill probabilities without "
      "reference to timing. Delaying agent B therefore cannot change *where* B "
      "goes — only, very slightly, *how much* of B's steady-state field is "
      "realised inside the observation window, and that is the terminal residual "
      "again, bounded by `tau/T`. **Staggering can only move combination overlap "
      "in this model through an amplitude change, never through a geometry "
      "change.** The clinical reasons to stagger — resistance-reversal windows, "
      "cell-cycle synchronisation, sequential DNA-damage-and-repair-block, "
      "non-overlapping toxicity — are all absent. So take the ~0% as a statement "
      "about the model's structure and NOT as evidence that sequencing does not "
      "matter.")
    w("")
    adm = []
    for cname in ("shared", "clinical"):
        for e in combo[f"pairs_{cname}"]:
            j0 = e["arms"]["simultaneous"]["jaccard"]
            j1 = e["arms"]["B delayed, matched administered"]["jaccard"]
            adm.append((abs((j1 - j0) / j0) if j0 else 0.0,
                        f"{e['a']} + {e['b']} (`{cname}`)", j0, j1))
    adm.sort(reverse=True)
    w(f"The one place staggering does move something materially is the "
      f"**matched-administered-dose** row, and that movement is a window "
      f"artefact in exactly the sense named above — the delayed agent simply has "
      f"less of the window left in which to be counted. The largest such shift "
      f"is **{adm[0][1]}**, Jaccard {adm[0][2]:.6f} -> {adm[0][3]:.6f} "
      f"({adm[0][0]*100:.2f}%), against a largest matched-COUNTED shift of "
      f"{max(stag_rel)*100:.3f}%. The gap between those two numbers is the whole "
      f"reason both legs are reported.")
    w("")
    w("---")
    w("")

    # ------------------------------------------------------------------ 6/7
    w("## 6. Verdict")
    w("")
    w("> **1. At matched cumulative dose over a real 28-day cycle, this model "
      "cannot tell the real temozolomide regimens apart.** The spread in mean "
      f"tumour exposure is {_pct(out['matched_exposure_spread'], 6)} across nine "
      f"arms spanning 1 to 28 fractions, against the derived bound "
      f"{_pct(out['matched_bound_tau_over_T'], 6)}. Regimens that finish dosing "
      "more than a day before the window closes are identical to each other to "
      "the stored digit.")
    w(">")
    w("> **2. That indistinguishability is a structural fact about transport, "
      "not a finding about temozolomide.** The entire schedule channel in a "
      "linear transport model is the terminal residual, bounded by "
      f"`tau_tissue / T_window`. Temozolomide equilibrates in this tissue in "
      f"{tau:.0f} s against a 28-day window, so the ratio is ~1e-4. The same "
      "machinery would show a large effect for a slowly-equilibrating molecule "
      "(an antibody, `tau` of hours to days) or over a short window — section 4.4 "
      "demonstrates that by shrinking the window to 24 h and recovering a "
      f"{_pct(idy['exposure_spread_rel'], 3)} effect from the identical code.")
    w(">")
    w("> **3. Where the model DOES separate regimens, it is separating them by "
      "dose.** As prescribed, the dose-dense and metronomic schedules deliver "
      "2-3x the cumulative temozolomide of the standard adjuvant cycle, and the "
      "predicted ranking is monotone in that. **RTOG 0525 showed that this is "
      "the wrong prediction**: 2.1x the dose, no survival benefit, double the "
      "grade >=3 toxicity. Do not read section 3's ranking as a recommendation; "
      "read it as a demonstration of what an exposure-driven model does when the "
      "biology is not exposure-driven.")
    w(">")
    w("> **4. Staggering does not change which cells a combination kills.** "
      f"Largest relative Jaccard change across six pair-by-convention cells: "
      f"{max(stag_rel)*100:.3f}%. This is structural — there is no drug-drug "
      "coupling in the transport operator and Bliss independence is timing-blind "
      "— so it is a strong statement about the model and a weak one about "
      "sequencing therapy.")
    w(">")
    self_ratios = []
    dose2_ratios = []
    all_on = []
    for cname in ("shared", "clinical"):
        for e in combo[f"pairs_{cname}"]:
            s = e["arms"]["simultaneous"]
            self_ratios.append(s["kill_combo"] / e["kill_selfcombo_24h"])
            dose2_ratios.append(s["kill_combo"] / e["kill_2x_better_24h"])
            all_on.extend(v["obs_over_null"] for v in e["arms"].values())
    n_beat_self = sum(1 for r in self_ratios if r > 1.0)
    n_beat_2x = sum(1 for r in dose2_ratios if r > 1.0)
    w(f"> **5. Against the dose-matched controls, the combinations do not "
      f"survive.** Of the {len(self_ratios)} pair-by-convention cells, "
      f"**{n_beat_2x} beat 2x the better single agent** and "
      f"**{n_beat_self} beat Bliss(better + ITSELF)** (combo/self-combo range "
      f"{min(self_ratios):.2f}x-{max(self_ratios):.2f}x). This reproduces "
      "`results/combination.md`'s finding under a time-resolved solver: the "
      "apparent combination benefit over 2x dosing is the independence "
      "arithmetic meeting a saturating Hill curve, not spatial complementarity. "
      "The overlap ratios say the same thing from the other side: `obs/null` "
      f"ranges {min(all_on):.2f}x-{max(all_on):.2f}x across every pair, "
      "convention and arm — i.e. the two agents' killed populations co-localise "
      "**more** than a delivery-scrambled null, not less. There is no spatial "
      "complementarity to exploit on this section.")
    w("")
    w("---")
    w("")
    w("## 7. THE HONEST CEILING: what makes fractionation work is not in this model")
    w("")
    w("This section is not a caveat appended to a result. It is the result's "
      "load-bearing context, and the numbers above should not be quoted without "
      "it.")
    w("")
    w("**The schedule effect measured here is a TRANSPORT-ONLY FLOOR.** Time "
      "resolution was added to the transport solve. The pharmacodynamic map is "
      "unchanged: a static Hill function of window-averaged cumulative exposure. "
      "Everything that gives fractionation and dose-density their clinical "
      "rationale is therefore missing:")
    w("")
    w("| missing mechanism | why it is a schedule mechanism | evidence it is "
      "the real story |")
    w("|---|---|---|")
    w("| **MGMT / AGAT depletion and regeneration** | MGMT is a *consumable, "
      "stoichiometric, suicide* repair protein. Protracted dosing depletes it "
      "and it regenerates between doses — so the effective sensitivity is a "
      "state variable with its own time constant, which is precisely what a "
      "schedule acts on. **This model has no MGMT term at all** (MGMT is not on "
      "the Xenium panel; sensitivity is uniform 1.0). | Tolcher 2003 (PMID "
      "12671695): mean AGAT activity fell 72% after 7 days and 73% after 21 days "
      "of TMZ in PBMCs, and stayed depressed >= 7 days after dosing. Note the "
      "plateau: 7 -> 21 days bought essentially nothing, which itself undercuts "
      "the dose-dense premise. |")
    w("| **Inter-fraction DNA repair** | The entire radiobiological basis of "
      "fractionation. Absent. | — |")
    w("| **Cell-cycle redistribution and S-phase gating** | TMZ cytotoxicity is "
      "MMR-dependent and requires replication cycles, so *when* a cell is hit "
      "relative to its cycle matters. This model's cell states are oxygen-banded "
      "and static through the whole window. | mechanism widely described; not "
      "independently sourced in this report |")
    w("| **Regrowth / repopulation between fractions** | No proliferation is "
      "simulated during the window at all. | — |")
    w("| **Immune contribution, and TMZ-induced lymphopenia** | Both absent. The "
      "second is a *schedule-dependent harm* that dose-dense regimens are known "
      "to cause. | Gilbert 2013: grade >= 3 toxicity 34% -> 53% on dose-dense, "
      "largely lymphopenia and fatigue |")
    w("| **Normal-tissue recovery — the entire denominator of the therapeutic "
      "index** | Every number in this report is tumour kill alone. A schedule "
      "that is near-identical for the tumor is emphatically not identical for "
      "the patient. | Gilbert 2013, as above |")
    w("")
    w("**And the clinical record says exposure is not the variable.** Ostermann "
      "2004 (PMID 15173079), which measured both plasma and CSF temozolomide in "
      "n=35 glioma patients, concluded that *\"systemic or cerebral exposures are "
      "not better predictors than the cumulative dose alone for both efficacy and "
      "safety\"* — and RTOG 0525 then showed cumulative dose itself does not "
      "predict survival. A model whose output is a monotone function of tumour "
      "exposure is, by construction, modelling the wrong variable for this "
      "question. **What it can legitimately do is bound the transport "
      "contribution, which is what it does above.**")
    w("")
    w("> **Do not present any number in this file as a clinical recommendation.** "
      "A 0.0001% exposure difference between two regimens is not a reason to "
      "choose one. It is a demonstration that transport is not where the "
      "difference lives.")
    w("")
    w("---")
    w("")

    # ------------------------------------------------------------------ 8
    w("## 8. Limitations (the full list)")
    w("")
    lims = [
        "**No outcome validation, at all.** Untreated resection; no measured "
        "response for any agent on this tissue. Nothing here has been scored "
        "against reality.",
        "**Not clinical dosing guidance.** Research prototype. No regulatory, "
        "diagnostic or treatment-selection use.",
        "**Absolute kill on this section is not defensible** — "
        "`results/cell_state_impact.md` establishes this directly. Rankings and "
        "relative comparisons are what this file is for. Every absolute kill "
        "here should be read as a unit-carrying intermediate.",
        f"**The schedule effect is bounded by `tau_tissue / T_window` and is "
        f"therefore window-conditional.** {tau:.0f} s / 28 d = "
        f"{_pct(out['matched_bound_tau_over_T'], 6)}. Section 4.3 shows the "
        "measured value tracking `1/T` across three windows. Any schedule "
        "percentage from this model is meaningless without its window.",
        "**Pharmacodynamics is a single static Hill map on window-averaged "
        "exposure — and this is the binding limitation on the whole question.** "
        "No MGMT depletion, no inter-fraction repair, no cell-cycle "
        "redistribution, no regrowth, no resistance evolution, no immune "
        "contribution, no normal-tissue toxicity. Section 7.",
        "**Kill depends on `AUC / T_window`, so lengthening the observation "
        "window mechanically lowers predicted kill at fixed total dose.** This "
        "is why the 28-day-window kills here are far below the 6 h-window kills "
        "in `results/dose_schedule.md` and `results/combination.md`. The window "
        "is held EXACTLY constant across everything compared, so the comparisons "
        "are sound; cross-file comparisons of absolute kill are not.",
        "**Bliss independence is an assumption.** No synergy, no antagonism, no "
        "drug-drug PK interaction, no shared target, no efflux competition. The "
        "combination results model SPATIAL complementarity only and are a "
        "delivery-only floor, not a forecast of combination efficacy.",
        "**Staggering can only act through amplitude in this model.** There is "
        "no drug-drug coupling in the transport operator and Bliss independence "
        "is timing-blind, so the ~0% overlap change is a structural consequence "
        "and not an empirical discovery about sequencing.",
        "**The PK conversion is anchored on AUC and overstates the peak by "
        f"{pk['ratio']:.2f}x** (instantaneous absorption vs a real tmax of ~1 h). "
        "It also carries no MTIC species, no CSF partitioning, no protein "
        "binding, and no patient-specific mg/m2 mapping. Section 2.2.",
        "**Two of the three clinical peaks in the combination leg are weak.** "
        "Cisplatin's free-platinum peak is DERIVED (from unchanged-cisplatin "
        "Cmax divided by a species ratio), not measured. Doxorubicin's peak is "
        "sampling-time-dominated (distribution t1/2 ~4 min) with interpatient CV "
        "~65%. Treat both as order-of-magnitude.",
        "**Per-cell sensitivity is uniform 1.0.** Temozolomide's determinant "
        "(MGMT) is not on this panel, so TMZ resistance is unmodelled and its "
        "absolute kill is optimistic. The `uniform` mode makes the cross-drug "
        "combination comparison fair by making every agent equally an upper "
        "bound; it does not make any of them right.",
        f"**Cell-state banding dominates the ceiling.** "
        f"{out['ceiling_bands']['QUIESCENT']['fraction']*100:.1f}% of tumor "
        f"cells are labelled QUIESCENT (x0.5 kill) and "
        f"{out['ceiling_bands']['HYPOXIC']['fraction']*100:.1f}% HYPOXIC "
        "(x0.35) by the oxygen model. That labelling is a *model* assignment "
        "from a simulated oxygen field, not a measurement.",
        "**The perivascular gradient partly inherits the vascular source-field "
        "shape** (`tmesim/validation/vessel_profile.py`), and the exposure "
        "geometry every spatial statistic rests on inherits that confound.",
        f"**One calibration constant** (alpha = {out['alpha']:.1f}, fit once on "
        "doxorubicin) sets the perivascular sink for every drug.",
        "**Single section, single patient, 2D.** No inter-patient variability, "
        "no out-of-plane transport, no cohort, no replicate sections.",
        "**Regimens are modelled over ONE 28-day cycle, not a full course.** "
        "Six adjuvant cycles, the 42-day concomitant phase, and any "
        "cycle-to-cycle carry-over are not represented. With no regrowth and no "
        "resistance evolution in the model, repeating the cycle would add "
        "nothing but exposure.",
        "**The single-dose superposition is exact for this solver, but it is "
        "exact BECAUSE the solver is linear and time-invariant.** Any future "
        "change that makes transport nonlinear (saturable binding, "
        "concentration-dependent efflux, drug-induced vascular change) "
        "invalidates it. V1/V2/V3 in section 1.2 would catch that, and they "
        "should be re-run before the method is reused.",
    ]
    for i, l in enumerate(lims, 1):
        w(f"{i}. {l}")
    w("")
    w("---")
    w("")
    w("## References")
    w("")
    refs = [
        "Stupp R, Mason WP, van den Bent MJ, et al. Radiotherapy plus "
        "concomitant and adjuvant temozolomide for glioblastoma. "
        "*N Engl J Med* 2005;352(10):987-96. PMID 15758009. "
        "doi:10.1056/NEJMoa043330",
        "Gilbert MR, Wang M, Aldape KD, et al. Dose-dense temozolomide for newly "
        "diagnosed glioblastoma: a randomized phase III clinical trial. "
        "*J Clin Oncol* 2013;31(32):4085-91. PMID 24101040. "
        "doi:10.1200/JCO.2013.49.6968",
        "Wick W, Steinbach JP, Kuker WM, et al. One week on/one week off: a "
        "novel active regimen of temozolomide for recurrent glioblastoma. "
        "*Neurology* 2004;62(11):2113-5. PMID 15184628. "
        "doi:10.1212/01.wnl.0000127617.89363.84",
        "Wick A, Felsberg J, Steinbach JP, et al. Efficacy and tolerability of "
        "temozolomide in an alternating weekly regimen in patients with "
        "recurrent glioma. *J Clin Oncol* 2007;25(22):3357-61. PMID 17664483. "
        "doi:10.1200/JCO.2007.10.7722",
        "Weller M, Tabatabai G, Kastner B, et al. MGMT promoter methylation is a "
        "strong prognostic biomarker for benefit from dose-intensified "
        "temozolomide rechallenge in progressive glioblastoma: the DIRECTOR "
        "trial. *Clin Cancer Res* 2015;21(9):2057-64. PMID 25655102. "
        "doi:10.1158/1078-0432.CCR-14-2737",
        "Perry JR, Belanger K, Mason WP, et al. Phase II trial of continuous "
        "dose-intense temozolomide in recurrent malignant glioma: RESCUE study. "
        "*J Clin Oncol* 2010;28(12):2051-7. PMID 20308655. "
        "doi:10.1200/JCO.2009.26.5520",
        "Baker SD, Wirth M, Statkevich P, et al. Absorption, metabolism, and "
        "excretion of 14C-temozolomide following oral administration to patients "
        "with advanced cancer. *Clin Cancer Res* 1999;5(2):309-17. PMID 10037179",
        "Beale P, Judson I, Moore S, et al. Effect of gastric pH on the relative "
        "oral bioavailability and pharmacokinetics of temozolomide. "
        "*Cancer Chemother Pharmacol* 1999;44(5):389-94. PMID 10501912. "
        "doi:10.1007/s002800050994",
        "Ostermann S, Csajka C, Buclin T, et al. Plasma and cerebrospinal fluid "
        "population pharmacokinetics of temozolomide in malignant glioma "
        "patients. *Clin Cancer Res* 2004;10(11):3728-36. PMID 15173079. "
        "doi:10.1158/1078-0432.CCR-03-0807",
        "Tolcher AW, Gerson SL, Denis L, et al. Marked inactivation of "
        "O6-alkylguanine-DNA alkyltransferase activity with protracted "
        "temozolomide schedules. *Br J Cancer* 2003;88(7):1004-11. "
        "PMID 12671695. doi:10.1038/sj.bjc.6600827",
        "TEMODAR (temozolomide) US Prescribing Information, sections 2.1 and "
        "12.3 (dose escalation rule; Cmax 7.5 ug/mL and AUC 23.4 ug*h/mL after "
        "150 mg/m2; elimination half-life 1.8 h; linear kinetics 75-250 "
        "mg/m2/day; IV/oral bioequivalence).",
        "Himmelstein KJ, Patton TF, Belt RJ, et al. Clinical kinetics of intact "
        "cisplatin and some related species. *Clin Pharmacol Ther* "
        "1981;29(5):658-64. PMID 7194166. doi:10.1038/clpt.1981.91",
        "Nagai N, Kinoshita M, Ogata H, et al. Relationship between "
        "pharmacokinetics of unchanged cisplatin and nephrotoxicity after "
        "intravenous infusions of cisplatin to cancer patients. "
        "*Cancer Chemother Pharmacol* 1996;39(1-2):131-7. PMID 8995510. "
        "doi:10.1007/s002800050548",
        "Benjamin RS, Riggs CE, Bachur NR. Plasma pharmacokinetics of adriamycin "
        "and its metabolites in humans with normal hepatic and renal function. "
        "*Cancer Res* 1977;37(5):1416-20. PMID 856462",
        "Camaggi CM, Comparsi R, Strocchi E, et al. Epirubicin and doxorubicin "
        "comparative metabolism and pharmacokinetics. *Cancer Chemother "
        "Pharmacol* 1988;21(3):221-8. PMID 3162849. doi:10.1007/BF00262774",
        "Muller C, Chatelut E, Gualano V, et al. Cellular pharmacokinetics of "
        "doxorubicin in patients with chronic lymphocytic leukemia. "
        "*Cancer Chemother Pharmacol* 1993;32(5):379-84. PMID 8339389. "
        "doi:10.1007/BF00735923",
        "Speth PA, van Hoesel QG, Haanen C. Clinical pharmacokinetics of "
        "doxorubicin. *Clin Pharmacokinet* 1988;15(1):15-31. PMID 3042244. "
        "doi:10.2165/00003088-198815010-00002",
        "PubChem CID 5394 (temozolomide, C6H6N6O2, MW 194.15); CID 31703 "
        "(doxorubicin, MW 543.5).",
    ]
    for i, r in enumerate(refs, 1):
        w(f"{i}. {r}")
    w("")
    w("---")
    w("")
    w(f"*Generated by `experiments/regimen.py` — ONE reconstruction reused "
      f"everywhere, {out['n_pde_solves']} transport solves "
      f"(`PDEConfig(time_resolved=True)` on every one), {out['total_steps']:,} "
      f"solver steps, {out['timing']['total_s']/60:.1f} min total, dx = "
      f"{rep['dx_um']:.0f} um, dt = {DT_S:.0f} s, seed {SEED}, deterministic. "
      f"No file in `tmesim/` was modified.*")

    with open(path, "w") as fh:
        fh.write("\n".join(A) + "\n")


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dx", type=float, default=DX_UM)
    ap.add_argument("--cache", default=_DEFAULT_CACHE)
    ap.add_argument("--stage", default="all",
                    choices=["all", "compute", "report"])
    args = ap.parse_args()

    if args.stage in ("all", "compute"):
        out = compute(args.dx)
        with open(args.cache, "wb") as fh:
            pickle.dump(out, fh)
        print(f"cache -> {args.cache}")
    else:
        with open(args.cache, "rb") as fh:
            out = pickle.load(fh)

    if args.stage in ("all", "report"):
        write_report(out, _OUT_MD)
        print(f"report -> {_OUT_MD}")


if __name__ == "__main__":
    main()
