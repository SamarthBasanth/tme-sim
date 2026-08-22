"""Localise the "~14x too shallow" perivascular-gradient gap of
``results/kyle2004_penetration.md`` sec. 6 / C4. DECOMPOSE ONLY — CLOSE NOTHING.

WHAT THIS SCRIPT IS FOR
-----------------------
``results/kyle2004_penetration.md`` sec. 6.1 reports daunorubicin's modelled
near:far exposure ratio across a 140 um slab as **3.6457** against Kyle 2004's
measured **50 +- 15**, and calls that "~14x too shallow". It has been quoted all
week and never decomposed. This script decomposes it, and before decomposing it
verifies that the 14x is what it claims to be.

WHAT IS CHANGED
---------------
NOTHING. No file outside ``experiments/`` and ``results/accuracy_gap_*`` is
written. No default, constant or shipped parameter is modified.
``calibrate_alpha()`` is asserted live at ``22.223477711172563`` in every stage.

ON THE DIAGNOSTIC SWEEPS
------------------------
Several stages scale ``k_loss`` or ``D_eff`` on a LOCAL COPY of a ``DrugParams``
over a grid declared in this file's constants. That is a SENSITIVITY
MEASUREMENT of the readout, not a fit: the grid is fixed in advance, the whole
curve is reported, and NO value is adopted, written back, or recommended. The
project rule is "never tune a parameter to hit the anchor" and nothing here does
— the point of tracing R(lambda) is to show how badly the READOUT amplifies a
lambda error, which is the finding.

Run:
    ./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage verify
    ./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage anchors
    ./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage offslab
    ./.venv/bin/python experiments/accuracy_gap_decomposition.py --stage all
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import replace
from typing import Dict, List, Sequence

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import matplotlib
matplotlib.use("Agg")

import isomer_penetration_map as IPM  # noqa: E402  (experiments/ is on sys.path)

from tmesim.core import PDEConfig, SimConfig  # noqa: E402
from tmesim.physics.diffusion import ALPHA_TORTUOSITY  # noqa: E402
from tmesim.sim.engine import simulate_response  # noqa: E402
from tmesim.validation.calibration import calibrate_alpha  # noqa: E402

# --------------------------------------------------------------------------- #
# Constants — every one of these is stated in the same sentence as its numbers
# in results/accuracy_gap_decomposition.md.
# --------------------------------------------------------------------------- #
ALPHA_EXPECTED = 22.223477711172563

#: Kyle 2004's apparatus, as results/kyle2004_penetration.md sec. 6 reproduces it.
T_UM = 140.0            # measured MCC thickness
W_UM = 30.0             # "the first 30 um of tissue" on each side
ECM_ARM_A = 0.08747329648865268   # the GBM Xenium section's mean ECM — see sec. 3
DX_ARM_A = 2.0
EXPOSURE_S = 3600.0     # 1 h, matching the measurement

#: Kyle 2004 Table 1, drug-derived-fluorescence column (near:far, dimensionless).
KYLE_FLUOR = {"doxorubicin": (11.0, 5.0), "epirubicin": (12.0, 3.0),
              "daunorubicin": (50.0, 15.0)}

#: DIAGNOSTIC sensitivity grid. Declared here, before any run. Multiplies k_loss
#: on a LOCAL DrugParams copy. NOTHING on this grid is adopted.
K_SCALE_GRID: Sequence[float] = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0)
#: DIAGNOSTIC grid for the off-slab (real-section) transfer function. Kept short
#: because each entry is a full Xenium GBM solve.
K_SCALE_GRID_OFFSLAB: Sequence[float] = (1.0, 4.0, 16.0)
#: DIAGNOSTIC grid for the retardation/transient reading of Kyle's 1 h snapshot.
RETARDATION_GRID: Sequence[float] = (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 90.0,
                                     100.0, 200.0)
#: DIAGNOSTIC grid for the preparation inputs neither of which was measured on an
#: HCT-116 MCC.
ECM_GRID: Sequence[float] = (0.0, ECM_ARM_A, 0.2, 0.4)
PHI_GRID: Sequence[float] = (0.70, 0.7892, 0.85, 0.90)

_RESULTS = os.path.join(_ROOT, "results")
_DATA_PATH = os.path.join(_RESULTS, "accuracy_gap_decomposition_data.json")
_KYLE_DATA = os.path.join(_RESULTS, "kyle2004_penetration_data.json")


def _assert_alpha() -> float:
    a = float(calibrate_alpha())
    if abs(a - ALPHA_EXPECTED) > 1e-9:
        raise AssertionError(
            f"calibrate_alpha() = {a!r}, expected {ALPHA_EXPECTED!r}. Something "
            f"upstream moved; every number in this file is void.")
    return a


# --------------------------------------------------------------------------- #
# The analytic transfer function: lambda -> Kyle's near:far ratio
# --------------------------------------------------------------------------- #
def R_of_lambda(lam: float, thickness_um: float = T_UM,
                window_um: float = W_UM) -> float:
    """Steady-state near:far ratio for Kyle's geometry.

    One-dimensional slab of thickness ``T``, FIXED CONCENTRATION on the exposed
    face (Kyle's stirred reservoir, a Dirichlet condition) and NO FLUX on the
    sealed far face, with a first-order sink. The steady state is
    ``c(x) ~ cosh((T-x)/lam)``, and Kyle's readout is the mean over the first
    ``W`` um divided by the mean over the last ``W`` um:

        R(lam) = [sinh(T/lam) - sinh((T-W)/lam)] / sinh(W/lam)

    This is EXPONENTIAL in ``T/lam``. That is the whole point of sec. 1 of the
    report: a modest error in the length ``lam`` presents as a large error in
    ``R``.
    """
    T, W = float(thickness_um), float(window_um)
    try:
        return ((math.sinh(T / lam) - math.sinh((T - W) / lam))
                / math.sinh(W / lam))
    except OverflowError:
        # lam << W: sinh(z) -> e^z/2 and the expression -> e^((T-W)/lam), which
        # is astronomically large. Return inf rather than invent a number; the
        # bisection in lambda_of_R only needs the sign.
        return float("inf")


def lambda_of_R(R: float, thickness_um: float = T_UM,
                window_um: float = W_UM) -> float:
    """Invert :func:`R_of_lambda` by bisection. Monotone decreasing in lam."""
    lo, hi = 1e-3, 1e6
    for _ in range(400):
        mid = math.sqrt(lo * hi)
        if R_of_lambda(mid, thickness_um, window_um) > R:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def lambda_engine_um(dp, ecm: float) -> float:
    """The engine's own transport length for a slab at ECM ``ecm``.

    ``sqrt(D_eff * (1 - ALPHA_TORTUOSITY*ecm) / (binding + clearance))`` — the
    same ``D_eff`` the slab solver is handed by
    ``isomer_penetration_map.mcc_ratio``.
    """
    d = float(dp.diffusivity) * max(1.0 - ALPHA_TORTUOSITY * float(ecm), 1e-6)
    return math.sqrt(d / float(dp.binding + dp.clearance))


# --------------------------------------------------------------------------- #
# Slab readouts, driven through the frozen solver
# --------------------------------------------------------------------------- #
def slab_readouts(dp, ecm=None, dx_um=None, thickness_um=None, exposure_s=None,
                  window_um=None, source_density: float = 1.0) -> Dict:
    """Both estimators on one solve: the AUC ratio and the FINAL-SNAPSHOT ratio.

    ``isomer_penetration_map.mcc_ratio`` reports the AUC (time-integrated)
    ratio. Kyle froze the cultures immediately after exposure and read
    fluorescence, which is a SNAPSHOT. Reporting both separates the estimator
    term from the physics term.
    """
    # Resolved at CALL time, not bound as default arguments at import time.
    # Binding them as defaults made a planted defect (mutate the module
    # constant, re-run, expect the reproduction check to fire) silently miss —
    # the plant could not reach the code. See sec. 8 of the report.
    ecm = ECM_ARM_A if ecm is None else float(ecm)
    dx_um = DX_ARM_A if dx_um is None else float(dx_um)
    thickness_um = T_UM if thickness_um is None else float(thickness_um)
    exposure_s = EXPOSURE_S if exposure_s is None else float(exposure_s)
    window_um = W_UM if window_um is None else float(window_um)
    micro, nx = IPM._slab(ecm, thickness_um, dx_um)
    if source_density != 1.0:
        np.asarray(micro.fields.vascular_density)[:, 0] = float(source_density)
    d_eff = float(dp.diffusivity) * max(1.0 - ALPHA_TORTUOSITY * ecm, 1e-6)
    dt = min(2.0, IPM.MCC_CFL_SAFETY * dx_um * dx_um / max(d_eff, 1e-9))
    pde = PDEConfig(convection=False, dt=dt, total_time=float(exposure_s),
                    save_every=100000)
    res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
    auc = np.asarray(res.fields["auc"], dtype=np.float64).mean(axis=0)
    conc = np.asarray(res.transport.final_concentration,
                      dtype=np.float64).mean(axis=0)
    x = np.arange(nx) * float(dx_um)
    near = x <= window_um
    far = x >= (thickness_um - window_um)
    out = {"dx_um": float(dx_um), "ecm": float(ecm),
           "thickness_um": float(thickness_um), "exposure_s": float(exposure_s),
           "window_um": float(window_um), "source_density": float(source_density),
           "lambda_um": lambda_engine_um(dp, ecm)}
    for lab, f in (("auc", auc), ("conc", conc)):
        out[f"R_{lab}"] = float(f[near].mean() / f[far].mean())
    out["R_analytic_at_lambda"] = R_of_lambda(out["lambda_um"], thickness_um,
                                              window_um)
    return out


def _scaled(dp, k_scale: float = 1.0, d_scale: float = 1.0):
    """A LOCAL copy of ``dp`` with k and/or D multiplied. Diagnostic only."""
    return replace(dp,
                   diffusivity=float(dp.diffusivity) * float(d_scale),
                   binding=float(dp.binding) * float(k_scale),
                   clearance=float(dp.clearance) * float(k_scale))


# --------------------------------------------------------------------------- #
# The transient / store reading of Kyle's 1 h snapshot
# --------------------------------------------------------------------------- #
def transient_near_far(d_app_um2_s: float, t_s: float = EXPOSURE_S,
                       thickness_um: float = T_UM, window_um: float = W_UM,
                       n_terms: int = 600, nx: int = 2001) -> float:
    """near:far at time ``t`` for PURE diffusion, NO loss, in Kyle's geometry.

    Dirichlet ``c=1`` at ``x=0``, no flux at ``x=T``, ``c(x>0, t=0) = 0``:

        c(x,t)/c0 = 1 - sum_n (4/((2n+1)pi)) sin(l_n x) exp(-l_n^2 D t),
        l_n = (2n+1) pi / (2T)

    The steady state of this problem is UNIFORM (ratio 1). Any ratio it produces
    at finite ``t`` is a TIME-SCALE effect and says nothing about a length scale.
    """
    x = np.linspace(0.0, thickness_um, nx)
    s = np.zeros_like(x)
    for n in range(n_terms):
        lam = (2 * n + 1) * math.pi / (2.0 * thickness_um)
        s += (4.0 / ((2 * n + 1) * math.pi)) * np.sin(lam * x) * math.exp(
            -lam * lam * float(d_app_um2_s) * float(t_s))
    c = 1.0 - s
    return float(c[x <= window_um].mean()
                 / c[x >= (thickness_um - window_um)].mean())


# --------------------------------------------------------------------------- #
# STAGE verify
# --------------------------------------------------------------------------- #
def stage_verify(verbose: bool = True) -> Dict:
    alpha = _assert_alpha()
    compiled = IPM.compile_all(alpha, two_route=True)
    dox = compiled["doxorubicin"]["drug_params"]
    dau = compiled["daunorubicin"]["drug_params"]
    out: Dict = {"alpha": alpha, "config": {
        "arm": "A (uniform slab, Kyle's geometry) unless stated",
        "dx_um": DX_ARM_A, "ecm": ECM_ARM_A, "thickness_um": T_UM,
        "window_um": W_UM, "exposure_s": EXPOSURE_S, "convection": False,
        "solver": "simulate_response, frozen"}}

    # --- V1. reproduce the published arm-A numbers before arguing on them ----
    rep = {n: slab_readouts(compiled[n]["drug_params"])
           for n in ("doxorubicin", "daunorubicin", "mitoxantrone")}
    published = {"doxorubicin": 2.997645248639645,
                 "daunorubicin": 3.6456683362306537}
    out["V1_reproduction"] = {
        n: {"R_auc_live": rep[n]["R_auc"],
            "R_auc_published": published.get(n),
            "rel_diff": (abs(rep[n]["R_auc"] - published[n]) / published[n]
                         if n in published else None),
            "R_conc_live": rep[n]["R_conc"],
            "lambda_um": rep[n]["lambda_um"]}
        for n in rep}

    # --- V2. what the 14x is a ratio of -------------------------------------
    R_dau = rep["daunorubicin"]["R_auc"]
    R_dox = rep["doxorubicin"]["R_auc"]
    m_dau, sd_dau = KYLE_FLUOR["daunorubicin"]
    m_dox, sd_dox = KYLE_FLUOR["doxorubicin"]
    out["V2_what_the_14x_is"] = {
        "numerator": "Kyle 2004 Table 1 fluorescence near:far for daunorubicin "
                     "(dimensionless concentration ratio)",
        "denominator": "engine arm-A AUC near:far for daunorubicin "
                       "(dimensionless exposure ratio)",
        "value": m_dau / R_dau,
        "band_from_measured_sd": [(m_dau - sd_dau) / R_dau,
                                  (m_dau + sd_dau) / R_dau],
        "note": "BOTH quantities are dimensionless. This is a ratio of ratios, "
                "NOT a ratio of two lengths.",
        "the_other_number_3p7x": {
            "R_sim": R_dau / R_dox, "R_measured": m_dau / m_dox,
            "value": (m_dau / m_dox) / (R_dau / R_dox),
            "note": "sec. 6's C3 '~3.7x compressed' is the FOLD comparison. It "
                    "is contained in the 14x: 14x = 3.7x (compound "
                    "under-separation) x the doxorubicin-row residual."},
    }

    # --- V3. the transfer function, and its amplification --------------------
    lam_lo, lam_hi = 5.0, 200.0
    curve = []
    lam = lam_lo
    while lam <= lam_hi + 1e-9:
        curve.append({"lambda_um": lam, "R": R_of_lambda(lam)})
        lam *= 1.15
    out["V3_transfer_function"] = {
        "form": "R(lam) = [sinh(T/lam) - sinh((T-W)/lam)] / sinh(W/lam), "
                f"T={T_UM} um, W={W_UM} um, steady state, Dirichlet near face, "
                "sealed far face",
        "curve": curve,
        "implied_lambda_um": {
            "engine_daunorubicin_R_auc": lambda_of_R(R_dau),
            "engine_doxorubicin_R_auc": lambda_of_R(R_dox),
            "kyle_daunorubicin_50": lambda_of_R(m_dau),
            "kyle_daunorubicin_35_minus1sd": lambda_of_R(m_dau - sd_dau),
            "kyle_daunorubicin_65_plus1sd": lambda_of_R(m_dau + sd_dau),
            "kyle_doxorubicin_11": lambda_of_R(m_dox),
            "kyle_doxorubicin_6_minus1sd": lambda_of_R(m_dox - sd_dox),
            "kyle_doxorubicin_16_plus1sd": lambda_of_R(m_dox + sd_dox),
        },
    }
    lam_e = rep["daunorubicin"]["lambda_um"]
    lam_k = lambda_of_R(m_dau)
    out["V3_transfer_function"]["headline"] = {
        "engine_lambda_um_at_arm_A_ecm": lam_e,
        "kyle_implied_lambda_um": lam_k,
        "lambda_error_fold": lam_e / lam_k,
        "lambda_error_band": [lam_e / lambda_of_R(m_dau - sd_dau),
                              lam_e / lambda_of_R(m_dau + sd_dau)],
        "D_over_k_error_fold": (lam_e / lam_k) ** 2,
        "mean_amplification_dlnR_dlnlambda": (
            math.log(m_dau / R_dau) / math.log(lam_e / lam_k)),
        "local_amplification_at_engine_lambda": (
            (math.log(R_of_lambda(lam_e * 1.01))
             - math.log(R_of_lambda(lam_e / 1.01))) / math.log(1.01 ** 2)),
        "local_amplification_at_kyle_lambda": (
            (math.log(R_of_lambda(lam_k * 1.01))
             - math.log(R_of_lambda(lam_k / 1.01))) / math.log(1.01 ** 2)),
    }

    # --- V4. does the ENGINE follow the analytic transfer function? ----------
    sweep = []
    for s in K_SCALE_GRID:
        r = slab_readouts(_scaled(dau, k_scale=s))
        sweep.append({"k_scale": s, "lambda_um": r["lambda_um"],
                      "R_auc_engine": r["R_auc"], "R_conc_engine": r["R_conc"],
                      "R_analytic": r["R_analytic_at_lambda"],
                      "engine_auc_over_analytic": (r["R_auc"]
                                                   / r["R_analytic_at_lambda"]),
                      "engine_conc_over_analytic": (
                          r["R_conc"] / r["R_analytic_at_lambda"])})
    out["V4_engine_vs_analytic"] = {
        "drug": "daunorubicin", "k_scaled_on_a_local_copy_only": True,
        "sweep": sweep}

    # --- V5. estimator terms -------------------------------------------------
    out["V5_estimator"] = {
        "auc_over_snapshot": {n: rep[n]["R_auc"] / rep[n]["R_conc"]
                              for n in rep},
        "note": "the engine reports an AUC ratio; Kyle read a SNAPSHOT of "
                "frozen tissue. At the engine's own settings these differ by "
                "the factor above.",
        "dx_sensitivity": [
            {"dx_um": d,
             "R_auc": slab_readouts(dau, dx_um=d)["R_auc"]}
            for d in (1.0, 2.0, 4.0)],
        "window_sensitivity": [
            {"window_um": w, "R_auc": slab_readouts(dau, window_um=w)["R_auc"],
             "R_analytic": R_of_lambda(lam_e, window_um=w)}
            for w in (15.0, 30.0, 45.0)],
        "exposure_sensitivity": [
            {"exposure_s": t, "R_auc": slab_readouts(dau, exposure_s=t)["R_auc"]}
            for t in (600.0, 3600.0, 3600.0 * 6)],
    }

    # --- V6. source geometry / source strength ------------------------------
    src = [{"source_density": s,
            "R_auc": slab_readouts(dau, source_density=s)["R_auc"]}
           for s in (0.25, 1.0)]
    out["V6_source"] = {
        "sweep": src,
        "max_rel_spread": (max(r["R_auc"] for r in src)
                           / min(r["R_auc"] for r in src) - 1.0),
        "note": "arm A has NO vessel map and NO delivery geometry — it is a "
                "uniform slab with a planar source in column 0. The source "
                "candidate is ruled out for the 14x BY CONSTRUCTION; this "
                "sweep confirms the readout is invariant to source LEVEL too."}

    # --- V7. preparation inputs neither measured on an MCC -------------------
    out["V7_preparation_inputs"] = {
        "ecm_sweep": [{"ecm": e, "lambda_um": lambda_engine_um(dau, e),
                       "R_analytic": R_of_lambda(lambda_engine_um(dau, e))}
                      for e in ECM_GRID],
        "phi_sweep": [],
        "note": "arm A's slab uses ecm = 0.0875, the mean of the 10x Xenium "
                "human GBM SECTION, and a FIXED cell volume fraction phi = "
                "0.7892 (two_route.DEFAULT_CELL_VOLUME_FRACTION). Neither was "
                "measured on an HCT-116 multicell layer, which is the "
                "preparation being scored."}
    from tmesim.molecular import ecm_mesh
    base_h = float(ecm_mesh.cell_exclusion(0.7892)["hindrance"])
    for phi in PHI_GRID:
        h = float(ecm_mesh.cell_exclusion(phi)["hindrance"])
        lam_phi = lambda_engine_um(_scaled(dau, d_scale=h / base_h), ECM_ARM_A)
        out["V7_preparation_inputs"]["phi_sweep"].append(
            {"phi": phi, "hindrance": h, "hindrance_rel_to_default": h / base_h,
             "lambda_um": lam_phi, "R_analytic": R_of_lambda(lam_phi)})

    # --- V8. the transient / store reading -----------------------------------
    d_eff_dau = (float(dau.diffusivity)
                 * max(1.0 - ALPHA_TORTUOSITY * ECM_ARM_A, 1e-6))
    out["V8_transient_store"] = {
        "D_eff_um2_s_daunorubicin_at_arm_A_ecm": d_eff_dau,
        "engine_is_equilibrated": {
            "one_over_k_s": 1.0 / float(dau.binding + dau.clearance),
            "T2_over_D_s": T_UM ** 2 / d_eff_dau,
            "exposure_s": EXPOSURE_S},
        "retardation_sweep": [
            {"R_f": rf, "D_apparent_um2_s": d_eff_dau / rf,
             "tau_dimensionless": d_eff_dau * EXPOSURE_S / (rf * T_UM ** 2),
             "near_far_at_1h_no_loss": transient_near_far(d_eff_dau / rf)}
            for rf in RETARDATION_GRID],
        "note": "PURE retarded diffusion with ZERO loss. Its steady state is "
                "uniform (ratio 1). Any ratio here is a TIME-SCALE effect."}

    # --- V9. the log-additive split of the lambda error ----------------------
    lam_e_dox = rep["doxorubicin"]["lambda_um"]
    lam_k_dox = lambda_of_R(m_dox)
    tot = math.log(lam_e / lam_k)
    dox_term = math.log(lam_e_dox / lam_k_dox)
    cmp_term = tot - dox_term
    # The engine's arm-A slab runs at ecm 0.0875 while calibrate_alpha() pins
    # doxorubicin's L at ecm_reference 0.4. Quantify that inconsistency.
    lam_dau_at_cal_ecm = lambda_engine_um(dau, 0.4)
    out["V9_split"] = {
        "space": "log of the transport length lambda; terms are additive",
        "total_ln": tot, "total_fold": math.exp(tot),
        "doxorubicin_inherited": {
            "ln": dox_term, "fold": math.exp(dox_term),
            "share": abs(dox_term) / (abs(dox_term) + abs(cmp_term)),
            "what": "the residual on the compound alpha is FIT to. Not a "
                    "prediction: calibrate_alpha() pins doxorubicin's L to the "
                    "Primeau 2005 40-50 um half-distance anchor (PMID "
                    "16361566). This term IS the Primeau-vs-Kyle anchor "
                    "conflict."},
        "daunorubicin_specific": {
            "ln": cmp_term, "fold": math.exp(cmp_term),
            "share": abs(cmp_term) / (abs(dox_term) + abs(cmp_term)),
            "what": "the genuine out-of-sample prediction error — the "
                    "compound-to-compound under-separation sec. 6 reports as "
                    "'~3.7x compressed', expressed as a length"},
        "primeau_band_L_um": [40.0 / math.log(2.0), 50.0 / math.log(2.0)],
        "kyle_dox_band_lambda_um": [lambda_of_R(m_dox + sd_dox),
                                    lambda_of_R(m_dox - sd_dox)],
        "bands_overlap": not (lambda_of_R(m_dox - sd_dox)
                              < 40.0 / math.log(2.0)),
        "arm_A_ecm_inconsistency": {
            "lambda_at_arm_A_ecm_0p0875": lam_e,
            "lambda_at_calibration_ecm_0p4": lam_dau_at_cal_ecm,
            "R_at_arm_A_ecm": R_of_lambda(lam_e),
            "R_at_calibration_ecm": R_of_lambda(lam_dau_at_cal_ecm),
            "fold_on_R": R_of_lambda(lam_dau_at_cal_ecm) / R_of_lambda(lam_e),
            "what": "results/kyle2004_penetration.md sec. 4's L column is "
                    "computed at ecm 0.4; arm A's slab is solved at ecm "
                    "0.0875. The file prints an L next to a run that does not "
                    "have that L."},
    }

    if verbose:
        _print_verify(out)
        v9 = out["V9_split"]
        print(f"\nV9 log-additive split of the {v9['total_fold']:.3f}x lambda "
              f"error (daunorubicin):")
        print(f"   doxorubicin-inherited (alpha/Primeau anchor): "
              f"{v9['doxorubicin_inherited']['fold']:.3f}x  "
              f"({v9['doxorubicin_inherited']['share']:.1%})")
        print(f"   daunorubicin-specific (out-of-sample):        "
              f"{v9['daunorubicin_specific']['fold']:.3f}x  "
              f"({v9['daunorubicin_specific']['share']:.1%})")
        print(f"   Primeau L band {v9['primeau_band_L_um'][0]:.1f}-"
              f"{v9['primeau_band_L_um'][1]:.1f} um vs Kyle-dox lambda band "
              f"{v9['kyle_dox_band_lambda_um'][0]:.1f}-"
              f"{v9['kyle_dox_band_lambda_um'][1]:.1f} um  -> overlap: "
              f"{v9['bands_overlap']}")
        print(f"   arm-A ECM inconsistency worth "
              f"{v9['arm_A_ecm_inconsistency']['fold_on_R']:.3f}x on R")
    return out


def _print_verify(o: Dict) -> None:
    print("\n" + "=" * 84)
    print("STAGE verify — is the 14x what it claims to be?")
    print("=" * 84)
    print(f"\nalpha = {o['alpha']!r}  (asserted)")
    print(f"config: {o['config']}")
    print("\nV1 reproduction of results/kyle2004_penetration.md sec. 6.1:")
    for n, v in o["V1_reproduction"].items():
        rd = "" if v["rel_diff"] is None else f"  rel-diff {v['rel_diff']:.2e}"
        print(f"  {n:<15} R_auc {v['R_auc_live']:8.4f} vs published "
              f"{v['R_auc_published']}{rd}   R_conc {v['R_conc_live']:8.4f}"
              f"   lambda {v['lambda_um']:7.2f} um")
    v2 = o["V2_what_the_14x_is"]
    print(f"\nV2 the 14x = {v2['value']:.3f}  (band {v2['band_from_measured_sd'][0]:.2f}"
          f"-{v2['band_from_measured_sd'][1]:.2f}) — a ratio of two DIMENSIONLESS ratios")
    print(f"   the 3.7x   = {v2['the_other_number_3p7x']['value']:.3f}")
    h = o["V3_transfer_function"]["headline"]
    print(f"\nV3 engine lambda {h['engine_lambda_um_at_arm_A_ecm']:.2f} um vs "
          f"Kyle-implied {h['kyle_implied_lambda_um']:.2f} um")
    print(f"   lambda error {h['lambda_error_fold']:.3f}x "
          f"(band {h['lambda_error_band'][0]:.2f}-{h['lambda_error_band'][1]:.2f}x)"
          f"   D/k error {h['D_over_k_error_fold']:.2f}x")
    print(f"   amplification dlnR/dlnlambda: mean over the interval "
          f"{h['mean_amplification_dlnR_dlnlambda']:.2f}, local at engine "
          f"{h['local_amplification_at_engine_lambda']:.2f}, local at Kyle "
          f"{h['local_amplification_at_kyle_lambda']:.2f}")
    print("\nV4 engine slab vs analytic transfer function (daunorubicin, k scaled "
          "on a LOCAL copy):")
    print(f"   {'k_scale':>8}{'lambda':>10}{'R_auc':>10}{'R_conc':>10}"
          f"{'R_analytic':>12}{'auc/an':>9}{'conc/an':>9}")
    for r in o["V4_engine_vs_analytic"]["sweep"]:
        print(f"   {r['k_scale']:>8.2f}{r['lambda_um']:>10.2f}"
              f"{r['R_auc_engine']:>10.3f}{r['R_conc_engine']:>10.3f}"
              f"{r['R_analytic']:>12.3f}{r['engine_auc_over_analytic']:>9.4f}"
              f"{r['engine_conc_over_analytic']:>9.4f}")
    print("\nV5 estimator terms:")
    print(f"   AUC/snapshot: {o['V5_estimator']['auc_over_snapshot']}")
    print(f"   dx:        {o['V5_estimator']['dx_sensitivity']}")
    print(f"   window:    {o['V5_estimator']['window_sensitivity']}")
    print(f"   exposure:  {o['V5_estimator']['exposure_sensitivity']}")
    print(f"\nV6 source level: spread {o['V6_source']['max_rel_spread']:.2e} "
          f"over {[r['source_density'] for r in o['V6_source']['sweep']]}")
    print("\nV7 preparation inputs:")
    for r in o["V7_preparation_inputs"]["ecm_sweep"]:
        print(f"   ecm {r['ecm']:.4f}  lambda {r['lambda_um']:7.2f} um  "
              f"R {r['R_analytic']:7.3f}")
    for r in o["V7_preparation_inputs"]["phi_sweep"]:
        print(f"   phi {r['phi']:.4f}  hindrance {r['hindrance']:.4f} "
              f"({r['hindrance_rel_to_default']:.3f}x default)  lambda "
              f"{r['lambda_um']:7.2f} um  R {r['R_analytic']:7.3f}")
    v8 = o["V8_transient_store"]
    print(f"\nV8 transient/store. engine 1/k = "
          f"{v8['engine_is_equilibrated']['one_over_k_s']:.0f} s, T^2/D = "
          f"{v8['engine_is_equilibrated']['T2_over_D_s']:.0f} s, exposure "
          f"{v8['engine_is_equilibrated']['exposure_s']:.0f} s -> equilibrated")
    for r in v8["retardation_sweep"]:
        print(f"   R_f {r['R_f']:7.1f}  D_app {r['D_apparent_um2_s']:9.4f}  "
              f"tau {r['tau_dimensionless']:8.4f}  near:far@1h "
              f"{r['near_far_at_1h_no_loss']:12.3f}")


# --------------------------------------------------------------------------- #
# STAGE anchors — the only thing in the repo that separates D_eff from k
# --------------------------------------------------------------------------- #
def stage_anchors(verbose: bool = True) -> Dict:
    """lambda for tirapazamine, from two SEPARATELY measured external anchors.

    ``lambda = sqrt(D/k)`` cannot separate ``D`` from ``k``. Tirapazamine can,
    because its ``D_eff`` and its ``k`` were each measured, and both are already
    in ``tmesim/validation/scorecard.py`` as separate rows. Nobody has combined
    them into a length. This does — LIVE, not from the prose.
    """
    _assert_alpha()
    from tmesim.molecular import compile_drug, descriptors
    from tmesim.validation.scorecard import _import_experiment
    AP = _import_experiment("anchor_prototype")

    mol = descriptors.from_smiles("tirapazamine", AP.TPZ_SMILES)
    prof = compile_drug.build_profile(mol, ecm_reference=AP.ECM_REFERENCE)
    alpha = float(calibrate_alpha())
    d_eng = float(prof.effective_diffusivity)
    k_eng = float(prof.clearance + alpha * prof.tissue_binding_rate)

    d_meas = {ln: float(v[0]) for ln, v in AP.D_MCL_MEASURED_UM2_S.items()}
    k_meas = {ln: float(v[0]) for ln, v in AP.K_CONSUMPTION_MEASURED_PER_S.items()}
    d_lo, d_hi = min(d_meas.values()), max(d_meas.values())
    k_lo, k_hi = min(k_meas.values()), max(k_meas.values())

    # The MCL-APPROPRIATE arm. results/cellular_exclusion.md sec. 4.1 and
    # scorecard.py both say the arm to quote against an MCL is the
    # IMPERMEABLE-CELL TERM ALONE, because an MCL is packed cells with little
    # matrix. Computed live rather than quoted, the same way scorecard.py's
    # d_armB is.
    from tmesim.molecular import ecm_mesh
    from tmesim.validation.scorecard import PHI_CELL_BREAST_SECTION
    cell_scale = float(ecm_mesh.cell_exclusion(PHI_CELL_BREAST_SECTION)
                       ["hindrance"])
    d_armB = float(prof.free_diffusivity) * cell_scale

    lam_eng = math.sqrt(d_eng / k_eng)
    lam_armB = math.sqrt(d_armB / k_eng)
    lam_meas_lo = math.sqrt(d_lo / k_hi)   # shortest measured length
    lam_meas_hi = math.sqrt(d_hi / k_lo)   # longest measured length

    out = {
        "compound": "tirapazamine",
        "ecm_reference": float(AP.ECM_REFERENCE),
        "engine": {"D_eff_um2_s": d_eng, "k_per_s": k_eng,
                   "lambda_um": lam_eng, "alpha": alpha},
        "measured": {"D_eff_um2_s": d_meas, "k_per_s": k_meas,
                     "D_band": [d_lo, d_hi], "k_band": [k_lo, k_hi],
                     "lambda_band_um": [lam_meas_lo, lam_meas_hi]},
        "D_error_fold_engine_over_measured": [d_eng / d_hi, d_eng / d_lo],
        "k_error_fold_engine_over_measured": [k_eng / k_hi, k_eng / k_lo],
        "lambda_error_fold_engine_over_measured": [lam_eng / lam_meas_hi,
                                                   lam_eng / lam_meas_lo],
        "mcl_appropriate_arm": {
            "D_eff_um2_s": d_armB, "k_per_s": k_eng, "lambda_um": lam_armB,
            "D_inside_measured_band": bool(d_lo <= d_armB <= d_hi),
            "lambda_error_fold": [lam_armB / lam_meas_hi,
                                  lam_armB / lam_meas_lo],
            "what": "the IMPERMEABLE-cell term alone, which "
                    "results/cellular_exclusion.md sec. 4.1 and scorecard.py "
                    "both name as the MCL-appropriate arm. If D_eff lands "
                    "INSIDE its measured band on this arm, then on the "
                    "preparation Kyle actually used the D term contributes "
                    "NOTHING to the lambda error and k carries all of it."},
        "log_split": {},
        "caveat": ("scorecard.py declares the k row NOT_TESTABLE because the "
                   "measured constant is anoxic bioreductive metabolism of a "
                   "hypoxia-activated prodrug and the engine's k is a generic "
                   "binding+clearance sink. That objection is about MECHANISM, "
                   "not about commensurability of the number: both sit in the "
                   "-kC slot of the same diffusion-reaction equation and both "
                   "were fitted from the same MCL slab flux, so the LENGTH "
                   "sqrt(D/k) they jointly determine is the quantity that "
                   "experiment actually measures well. It is reported here as "
                   "a DERIVED row and must be quoted with this caveat."),
    }
    # Log-additive split of the lambda error into its D and k halves.
    for lab, dm, km in (("vs_low_D_high_k", d_lo, k_hi),
                        ("vs_high_D_low_k", d_hi, k_lo)):
        t_k = 0.5 * math.log(km / k_eng)
        t_d = -0.5 * math.log(dm / d_eng)
        out["log_split"][lab] = {
            "k_term_ln": t_k, "D_term_ln": t_d, "total_ln": t_k + t_d,
            "lambda_fold": math.exp(t_k + t_d),
            "k_share_of_magnitude": abs(t_k) / (abs(t_k) + abs(t_d)),
            "D_share_of_magnitude": abs(t_d) / (abs(t_k) + abs(t_d)),
            "D_term_sign": "same direction" if t_d > 0 else "OPPOSITE — "
                           "D_eff being too SMALL partially masks the k error"}
    if verbose:
        print("\n" + "=" * 84)
        print("STAGE anchors — separating D_eff from k with tirapazamine")
        print("=" * 84)
        print(f"  engine  D_eff {d_eng:8.2f} um^2/s   k {k_eng:.4e} /s   "
              f"lambda {lam_eng:7.2f} um   (ecm_reference "
              f"{AP.ECM_REFERENCE}, alpha {alpha:.6f})")
        print(f"  measured D_eff {d_lo:.0f}-{d_hi:.0f} um^2/s   k "
              f"{k_lo:.4e}-{k_hi:.4e} /s   lambda {lam_meas_lo:.2f}-"
              f"{lam_meas_hi:.2f} um")
        print(f"  D error {out['D_error_fold_engine_over_measured'][0]:.3f}-"
              f"{out['D_error_fold_engine_over_measured'][1]:.3f}x   "
              f"k error {out['k_error_fold_engine_over_measured'][0]:.4f}-"
              f"{out['k_error_fold_engine_over_measured'][1]:.4f}x   "
              f"lambda error "
              f"{out['lambda_error_fold_engine_over_measured'][0]:.3f}-"
              f"{out['lambda_error_fold_engine_over_measured'][1]:.3f}x")
        mb = out["mcl_appropriate_arm"]
        print(f"  MCL-APPROPRIATE arm (cell term alone): D_eff "
              f"{mb['D_eff_um2_s']:8.2f} um^2/s  INSIDE the measured "
              f"{d_lo:.0f}-{d_hi:.0f} band? {mb['D_inside_measured_band']}  "
              f"lambda {mb['lambda_um']:7.2f} um  error "
              f"{mb['lambda_error_fold'][0]:.2f}-"
              f"{mb['lambda_error_fold'][1]:.2f}x")
        for lab, v in out["log_split"].items():
            print(f"  {lab}: lambda {v['lambda_fold']:.3f}x, k share "
                  f"{v['k_share_of_magnitude']:.1%}, D share "
                  f"{v['D_share_of_magnitude']:.1%} ({v['D_term_sign']})")
    return out


# --------------------------------------------------------------------------- #
# STAGE prep — the preparation-composition term, promoted into the ledger
# --------------------------------------------------------------------------- #
#: DECLARED BEFORE THE RUN, with sources, exactly as the project requires.
#:
#: ECM for an HCT-116 multicell layer -> 0.0. NOT a free choice and NOT fitted.
#: results/cellular_exclusion.md sec. 4.1: "an MCL is packed cells with little
#: matrix and the repo's ECM legs describe matrix the MCL does not have", and
#: scorecard.py's own arm-B note: "the IMPERMEABLE-cell term alone, the
#: MCL-appropriate arm, since an MCL is packed cells with little matrix". This
#: is the repository's declared position about this preparation, applied.
ECM_MCL_DECLARED = 0.0
ECM_MCL_SOURCE = ("results/cellular_exclusion.md sec. 4.1 and "
                  "tmesim/validation/scorecard.py's arm-B note — an MCL is "
                  "packed cells with little matrix, so the ECM tortuosity leg "
                  "describes matrix this preparation does not have")
#: PHI for an HCT-116 multicell layer -> THERE IS NO SOURCED VALUE, and one is
#: NOT invented here. results/cellular_exclusion.md sec. 4.1, verbatim: "A
#: literature search for a measured cell (or extracellular) volume fraction of a
#: Wilson/Hicks MCL was run and returned nothing usable; that search did not
#: complete cleanly and is recorded as an open item rather than as a null."
#: phi therefore stays at its shipped value in the corrected ledger and appears
#: ONLY as a signed open-ended sensitivity. The DIRECTION is sourced from the
#: same paragraph: "If an MCL is denser than this section — which is the likely
#: direction, since an MCL has no vasculature, no necrosis and little stroma".
PHI_MCL_DECLARED = None
PHI_MCL_SOURCE = ("NO SOURCED VALUE EXISTS. results/cellular_exclusion.md "
                  "sec. 4.1 records a literature search that returned nothing "
                  "usable and logs it as an open item. Direction only is "
                  "sourced: an MCL is likely DENSER than the section phi was "
                  "measured on.")
#: REGISTERED PREDICTION, written before stage_prep was run: correcting ecm from
#: 0.0875 to 0.0 RAISES D_eff by 1/(1-0.6*0.0875) = 1.0554, RAISES lambda by
#: 1.0273, and therefore LOWERS R. The gap gets WORSE by about 1.05x. If the run
#: disagrees with this, the run is what is reported.
PREP_REGISTERED_PREDICTION = (
    "correcting ecm 0.0875 -> 0.0 lowers R by ~1.05x, i.e. the gap gets WORSE")


def stage_prep(verbose: bool = True) -> Dict:
    """Is the 14x partly an artefact of simulating the wrong preparation?"""
    alpha = _assert_alpha()
    compiled = IPM.compile_all(alpha, two_route=True)
    dau = compiled["daunorubicin"]["drug_params"]
    dox = compiled["doxorubicin"]["drug_params"]

    shipped = slab_readouts(dau)
    corrected = slab_readouts(dau, ecm=ECM_MCL_DECLARED)
    shipped_dox = slab_readouts(dox)
    corrected_dox = slab_readouts(dox, ecm=ECM_MCL_DECLARED)

    m_dau = KYLE_FLUOR["daunorubicin"][0]
    m_dox = KYLE_FLUOR["doxorubicin"][0]
    out = {
        "declared_before_the_run": {
            "ecm_mcl": ECM_MCL_DECLARED, "ecm_source": ECM_MCL_SOURCE,
            "phi_mcl": PHI_MCL_DECLARED, "phi_source": PHI_MCL_SOURCE,
            "registered_prediction": PREP_REGISTERED_PREDICTION},
        "config": "arm A, dx = 2 um, 140 um, 30 um windows, 1 h, convection OFF",
        "daunorubicin": {
            "shipped_ecm_0p0875": {"lambda_um": shipped["lambda_um"],
                                   "R_auc": shipped["R_auc"],
                                   "gap_vs_measured": m_dau / shipped["R_auc"]},
            "corrected_ecm_0p0": {"lambda_um": corrected["lambda_um"],
                                  "R_auc": corrected["R_auc"],
                                  "gap_vs_measured": m_dau / corrected["R_auc"]},
            "fold_on_R": corrected["R_auc"] / shipped["R_auc"],
            "fold_on_lambda": corrected["lambda_um"] / shipped["lambda_um"]},
        "doxorubicin": {
            "shipped_ecm_0p0875": {"lambda_um": shipped_dox["lambda_um"],
                                   "R_auc": shipped_dox["R_auc"],
                                   "gap_vs_measured": m_dox / shipped_dox["R_auc"]},
            "corrected_ecm_0p0": {"lambda_um": corrected_dox["lambda_um"],
                                  "R_auc": corrected_dox["R_auc"],
                                  "gap_vs_measured": m_dox / corrected_dox["R_auc"]},
            "fold_on_R": corrected_dox["R_auc"] / shipped_dox["R_auc"]},
        "phi_is_NOT_corrected_because": PHI_MCL_SOURCE,
        "phi_open_ended_sensitivity": [],
    }
    from tmesim.molecular import ecm_mesh
    base_h = float(ecm_mesh.cell_exclusion(0.7892)["hindrance"])
    for phi in PHI_GRID:
        h = float(ecm_mesh.cell_exclusion(phi)["hindrance"])
        lam = lambda_engine_um(_scaled(dau, d_scale=h / base_h),
                               ECM_MCL_DECLARED)
        out["phi_open_ended_sensitivity"].append(
            {"phi": phi, "hindrance": h, "lambda_um": lam,
             "R_analytic": R_of_lambda(lam),
             "residual_gap": m_dau / R_of_lambda(lam),
             "phi_at_which_gap_closes": None})
    # The phi that would close the gap entirely — reported ONLY to show that
    # phi can absorb an ARBITRARY amount of it, which is the argument for NOT
    # letting it absorb any until it is measured. NOT adopted, NOT recommended.
    lo, hi = 0.7892, 0.9999
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        h = float(ecm_mesh.cell_exclusion(mid)["hindrance"])
        lam = lambda_engine_um(_scaled(dau, d_scale=h / base_h),
                               ECM_MCL_DECLARED)
        if R_of_lambda(lam) < m_dau:
            lo = mid
        else:
            hi = mid
    out["phi_that_would_close_the_gap_DO_NOT_ADOPT"] = {
        "phi": 0.5 * (lo + hi),
        "why_reported": "to demonstrate that an UNMEASURED parameter can "
                        "absorb the entire residual. That is the argument for "
                        "measuring it, and the argument against using it. This "
                        "value is not adopted, not written back, and must "
                        "never be quoted as a result."}
    if verbose:
        print("\n" + "=" * 84)
        print("STAGE prep — was the wrong preparation simulated?")
        print("=" * 84)
        print(f"  DECLARED BEFORE THE RUN: ecm_MCL = {ECM_MCL_DECLARED} "
              f"({ECM_MCL_SOURCE[:60]}...)")
        print(f"  DECLARED BEFORE THE RUN: phi_MCL = {PHI_MCL_DECLARED} "
              f"-- no sourced value exists")
        print(f"  REGISTERED PREDICTION:   {PREP_REGISTERED_PREDICTION}")
        for n, v in (("daunorubicin", out["daunorubicin"]),
                     ("doxorubicin", out["doxorubicin"])):
            s, c = v["shipped_ecm_0p0875"], v["corrected_ecm_0p0"]
            print(f"  {n:<14} shipped ecm 0.0875: lambda {s['lambda_um']:6.2f} "
                  f"R {s['R_auc']:6.3f} gap {s['gap_vs_measured']:6.2f}x   ->  "
                  f"CORRECTED ecm 0.0: lambda {c['lambda_um']:6.2f} R "
                  f"{c['R_auc']:6.3f} gap {c['gap_vs_measured']:6.2f}x   "
                  f"({v['fold_on_R']:.4f}x on R)")
        print("  phi open-ended sensitivity (NOT a correction — no source):")
        for r in out["phi_open_ended_sensitivity"]:
            print(f"    phi {r['phi']:.4f}  lambda {r['lambda_um']:6.2f} um  "
                  f"R {r['R_analytic']:8.3f}  residual gap "
                  f"{r['residual_gap']:7.2f}x")
        print(f"    phi that would close it entirely: "
              f"{out['phi_that_would_close_the_gap_DO_NOT_ADOPT']['phi']:.4f} "
              f"— DO NOT ADOPT, see the note in the JSON")
    return out


# --------------------------------------------------------------------------- #
# STAGE offslab — does the lambda argument survive the real vessel field?
# --------------------------------------------------------------------------- #
def stage_offslab(verbose: bool = True) -> Dict:
    """Trace d_half(lambda) on the reconstructed section, not on a slab.

    Two independent measurements of the same transfer function:
      (a) FREE — the four compounds already saved in
          ``results/kyle2004_penetration_data.json`` span lambda 54-72 um at one
          fixed geometry (the tissue is built once and reused), so their
          d_half values trace the curve at no compute cost.
      (b) PAID — a diagnostic k-scaling on the real section, which extends the
          curve down towards Kyle's implied lambda. Each point is a full solve.
    """
    _assert_alpha()
    out: Dict = {"config": {
        "arm": "B (10x Xenium human GBM, hBrain 319-gene panel, dx12_repaired)",
        "dx_um": 12.0, "convection": True, "total_time_s": 3600.0 * 6,
        "dt_s": 10.0, "bin_um": float(IPM.VP.BIN_UM),
        "readout": "model-free empirical half-distance of the AUC field binned "
                   "by direct_vessel_distance_um"}}

    # (a) free — from the saved data
    with open(_KYLE_DATA) as fh:
        saved = json.load(fh)
    ecm_b = float(saved["map"]["ecm_density"])
    free = []
    for n in ("doxorubicin", "daunorubicin", "mitoxantrone"):
        sh = saved["compounds"]["shipping"][n]
        lam = math.sqrt(sh["D_eff"] * max(1.0 - ALPHA_TORTUOSITY * ecm_b, 1e-6)
                        / sh["k_loss"])
        free.append({"drug": n, "lambda_um": lam,
                     "d_half_um": saved["map"]["drugs"][n]["auc_primary"][
                         "empirical_d_half_um"],
                     "fit_L_um": saved["map"]["drugs"][n]["auc_primary"][
                         "fit_L_um"]})
    out["free_points"] = {"ecm_section_mean": ecm_b, "points": free}
    if len(free) >= 2:
        xs = np.log([p["lambda_um"] for p in free])
        ys = np.log([p["d_half_um"] for p in free])
        slope, intercept = np.polyfit(xs, ys, 1)
        out["free_points"]["log_log_slope_dlnd_half_dlnlambda"] = float(slope)
        out["free_points"]["d_half_over_lambda"] = [
            p["d_half_um"] / p["lambda_um"] for p in free]
        out["free_points"]["note"] = (
            "d_half ~ lambda^%.2f over the narrow range lambda %.1f-%.1f um. "
            "For a 1-D exponential the slope would be exactly 1 and "
            "d_half/lambda = ln2 = 0.693; on the real vessel field it is "
            "%.3f-%.3f." % (slope, min(p['lambda_um'] for p in free),
                            max(p['lambda_um'] for p in free),
                            min(out['free_points']['d_half_over_lambda']),
                            max(out['free_points']['d_half_over_lambda'])))

    # (b) paid — diagnostic k scaling on the real section
    micro, _cn, _cr = IPM.VP.build_arm_microenvironment(dx_um=12.0, repair=True)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    endo = IPM.VP.endothelial_xy(micro.cells)
    dist = IPM.VP.direct_vessel_distance_um(micro.domain, endo)
    ecm_live = float(np.asarray(micro.fields.ecm_density)[tissue].mean())
    alpha = float(calibrate_alpha())
    dau = IPM.compile_all(alpha, two_route=True)["daunorubicin"]["drug_params"]
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True,
                    save_every=500)
    paid = []
    for s in K_SCALE_GRID_OFFSLAB:
        dp = _scaled(dau, k_scale=s)
        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
        auc = np.asarray(res.fields["auc"], dtype=np.float64)
        prof = IPM._profile_all(auc, dist, tissue, IPM.VP.BIN_UM)
        lam = lambda_engine_um(dp, ecm_live)
        paid.append({"k_scale": s, "lambda_um": lam,
                     "d_half_um": prof["empirical_d_half_um"],
                     "fit_L_um": prof["fit_L_um"],
                     "d_half_over_lambda": prof["empirical_d_half_um"] / lam})
        if verbose:
            print(f"  k x{s:<5.1f}  lambda {lam:7.2f} um  d_half "
                  f"{prof['empirical_d_half_um']:7.2f} um  fit_L "
                  f"{prof['fit_L_um']:7.2f} um  d_half/lambda "
                  f"{paid[-1]['d_half_over_lambda']:.3f}")
    out["paid_points"] = {"ecm_section_mean": ecm_live,
                          "drug": "daunorubicin", "points": paid}
    if len(paid) >= 2:
        xs = np.log([p["lambda_um"] for p in paid])
        ys = np.log([p["d_half_um"] for p in paid])
        slope, _ = np.polyfit(xs, ys, 1)
        out["paid_points"]["log_log_slope_dlnd_half_dlnlambda"] = float(slope)
    return out


# --------------------------------------------------------------------------- #
def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items() if not k.startswith("_")}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def main() -> Dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="verify",
                    choices=("verify", "anchors", "prep", "offslab",
                             "all"))
    args = ap.parse_args()
    res: Dict = {}
    if args.stage in ("verify", "all"):
        res["verify"] = stage_verify()
    if args.stage in ("anchors", "all"):
        res["anchors"] = stage_anchors()
    if args.stage in ("prep", "all"):
        res["prep"] = stage_prep()
    if args.stage in ("offslab", "all"):
        res["offslab"] = stage_offslab()
    if args.stage == "all":
        with open(_DATA_PATH, "w") as fh:
            json.dump(_jsonable(res), fh, indent=1)
        print(f"\nwrote {_DATA_PATH}")
    return res


if __name__ == "__main__":
    main()
