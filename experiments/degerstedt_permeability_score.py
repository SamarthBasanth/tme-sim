"""Score the engine's doxorubicin MEMBRANE permeability against Degerstedt 2024.

PRE-REGISTERED. Stage `prereg` freezes the criteria to
`experiments/degerstedt_prereg_2026-10-10.json` and REFUSES to overwrite it.
Stage `score` refuses to run unless that file exists, records its SHA-256 in the
output, and never reads anything from it except to echo it.

NOTHING IN `tmesim/` IS TOUCHED. No default moves. `alpha` is never recomputed
with a changed input; it is asserted live against 22.223477711172563 at the start
and end of every stage. The measured permeability is substituted into a LOCAL
call to `two_route.two_route_diffusivity` as a probe; the shipped compile path is
re-run unmodified beside it so the two can be differenced.

The measurement (verified verbatim from the CC BY 4.0 full text, PMC10927899):

  Degerstedt O, O'Callaghan P, Clavero AL, Grasjo J, Eriksson O, Sjogren E,
  Hansson P, Heindryckx F, Kreuger J, Lennernas H. "Quantitative imaging of
  doxorubicin diffusion and cellular uptake in biomimetic gels with human liver
  tumor cells." Drug Deliv Transl Res 2024;14(4):970-983.
  PMID 37824040, PMCID PMC10927899, doi:10.1007/s13346-023-01445-1.

  Eq. 2, verbatim:  P_cell = (k_cell x V_cell) / (A_cell x 60)
  "where V_cell and A_cell are the volume and area of the cells respectively.
   Approximating the cells are spherical, the ratio V_cell / A_cell can be
   simplified to r/3 where r is the cell radius."
  k_cell (HD HepG2) = 0.0163 +/- 0.001 min^-1  ->  P_cell = 9.00 +/- 0.74e-4 um/s
  k_cell (LD HepG2) = 0.0120 +/- 0.002 min^-1  ->  P_cell = 6.49 +/- 1.25e-4 um/s

Run:
    .venv/bin/python experiments/degerstedt_permeability_score.py --stage prereg
    .venv/bin/python experiments/degerstedt_permeability_score.py --stage score
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from typing import Dict

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, _REPO)

from tmesim.molecular import compile_drug, library, permeability, two_route  # noqa: E402
from tmesim.validation import calibration as C                              # noqa: E402

ALPHA_FROZEN = 22.223477711172563
PREREG_PATH = os.path.join(_REPO, "experiments", "degerstedt_prereg_2026-10-10.json")
SCORE_PATH = os.path.join(_REPO, "experiments", "degerstedt_score_2026-10-10.json")

# --- the measurement, as verified from the primary source -------------------
P_CELL_HD_UM_S = 9.00e-4          # HepG2, 2e6 cells/mL, +/- 0.74e-4
P_CELL_LD_UM_S = 6.49e-4          # HepG2, 1e6 cells/mL, +/- 1.25e-4
K_CELL_HD_PER_MIN = 0.0163        # +/- 0.001
# Secondary, quoted by Degerstedt as their ref [65]; NOT verified at its own
# source. DOX across an MDCK monolayer lining a cylindrical microvessel.
P_MDCK_MICROVESSEL_UM_S = 32e-4   # +/- 23e-4


def _guard() -> float:
    a = C.calibrate_alpha()
    assert a == ALPHA_FROZEN, f"alpha drifted: {a!r}"
    return a


def _sha256(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


# --------------------------------------------------------------------------- #
# Stage 1: freeze the pre-registration.
# --------------------------------------------------------------------------- #
def stage_prereg() -> None:
    alpha = _guard()
    if os.path.exists(PREREG_PATH):
        raise SystemExit(
            f"REFUSING to overwrite an existing pre-registration: {PREREG_PATH}\n"
            f"sha256={_sha256(PREREG_PATH)}"
        )

    prereg = {
        "frozen_utc_date": "2026-10-10",
        "alpha_at_freeze": alpha,
        "subject": (
            "The engine's doxorubicin CELL-MEMBRANE permeability, scored against "
            "the apparent cell permeability measured by Degerstedt 2024 "
            "(PMID 37824040, PMC10927899, CC BY 4.0)."
        ),
        "disclosure_this_prereg_is_not_blind": {
            "note": (
                "Full disclosure: the engine values and the measured value were "
                "BOTH in hand when these criteria were frozen, and a "
                "back-of-envelope conversion of P_cell into a tissue sink rate "
                "had already been done by hand. What this pre-registration buys "
                "is therefore NOT blindness; it is that the DECISION RULES and "
                "the predicted directions/bands below cannot be moved after the "
                "numbers land. A reader who wants a blind test should treat "
                "TEST 1 and TEST 2 as bookkeeping and TEST 3 and TEST 4 as the "
                "registered predictions, since neither of those had been "
                "computed with the real shipped inputs at freeze time."
            ),
            "in_hand_at_freeze": {
                "engine_two_route_P_mem_um_s": 9.63424115672293e-07,
                "engine_permeability_py_P_um_s": 17.07427760724692,
                "phi_cell": 0.7892,
                "cell_radius_um": 7.8517,
                "alpha": ALPHA_FROZEN,
                "measured_P_cell_HD_um_s": P_CELL_HD_UM_S,
            },
        },
        "calibration_guard": {
            "expression": "tmesim.validation.calibration.calibrate_alpha()",
            "must_equal": ALPHA_FROZEN,
        },
        "hard_constraints_acknowledged": [
            "No parameter is tuned to hit this anchor.",
            "alpha is NOT changed; it is asserted live.",
            "Nothing under tmesim/ is edited.",
            "A failure is a recorded failure.",
            "No absolute kill fraction is produced; no cross-compound ranking "
            "claim is made (WITHDRAWN 2026-08-21).",
        ],
        "TEST_0_commensurability_gate": {
            "decided_by": "reading the primary source BEFORE any computation",
            "criteria_all_must_hold": [
                "both quantities are in um/s",
                "both are defined as a flux per unit concentration difference "
                "across a CELL membrane",
                "the length scale converting a first-order rate into a "
                "permeability is a cell volume/area ratio in both",
            ],
            "verdict_recorded_at_freeze": "PASS",
            "evidence": {
                "paper_eq2": "P_cell = (k_cell x V_cell)/(A_cell x 60), V/A = r/3",
                "paper_k_cell_units": "min^-1 (stated verbatim)",
                "engine_entry_point": (
                    "two_route.two_route_factor: biot = sigma_c/(cell_radius_um "
                    "* P), i.e. 1/P is an interfacial resistance in series with "
                    "the inclusion interior -> P is flux per unit "
                    "concentration difference across the cell membrane, um/s"
                ),
            },
            "if_this_had_failed": (
                "every test below is VOID and the answer is reported as 'the "
                "comparison is not commensurable'"
            ),
        },
        "TEST_1_the_ratio": {
            "quantity": "R1 = P_cell_measured_HD / P_mem_engine_two_route",
            "predicted_direction": "engine LOWER, i.e. R1 > 1",
            "predicted_band": [400.0, 2000.0],
            "prior": "the scouting pass reported 935x; this re-derives it",
            "pass_if": "R1 > 1 and 400 <= R1 <= 2000",
        },
        "TEST_2_which_engine_term": {
            "quantity": "R2 = P_permeability_py_engine / P_cell_measured_HD",
            "predicted_direction": "engine HIGHER, i.e. R2 > 1e4",
            "status": "value was in hand at freeze; registered for the record, "
                      "not as a blind prediction",
            "registered_consequence_if_R1_gt_1_and_R2_gt_1": (
                "the engine carries TWO shipped doxorubicin 'membrane "
                "permeability' numbers that STRADDLE the measurement, so no "
                "one-sided statement of the form 'the engine's permeability is "
                "X-fold off' is available, and the headline must be the "
                "INTERNAL INCONSISTENCY rather than a single gap"
            ),
        },
        "TEST_3_consequence_for_the_shipped_transport_number": {
            "procedure": (
                "re-call two_route.two_route_diffusivity for doxorubicin with "
                "membrane_permeability_um_s overridden to the measured 9.00e-4 "
                "um/s, every other input at its shipped value, and difference "
                "the resulting effective_D against the shipped compile. Also "
                "compute the P -> infinity ceiling."
            ),
            "predicted_direction": "shipped D_eff moves UP, negligibly",
            "predicted_band": "|delta D_eff| / D_eff < 0.02",
            "reason_stated_in_advance": (
                "the dimensionless membrane resistance sigma_c/(a*P) is 2.18e7 "
                "at the shipped P, so a ~934x rise in P leaves it near 2.3e4, "
                "still >> 1, and the inclusion phase stays effectively "
                "impermeable"
            ),
            "pass_if": "relative change < 0.02",
            "registered_consequence_of_a_pass": (
                "the 10^3 permeability gap has NO consequence for any shipped "
                "doxorubicin transport number, which forbids calling it 'the "
                "most consequential number' on D_eff grounds and forces the "
                "alpha argument to stand entirely on its own"
            ),
        },
        "TEST_4_the_alpha_test": {
            "why": (
                "calibration.py attributes a ~300x sink to 'rapid cellular "
                "uptake plus high-capacity intracellular / DNA binding' and "
                "explicitly declines to use permeability as the proxy. The "
                "scouting pass claims the measurement REFRAMES alpha as a patch "
                "for a measurably wrong permeability. This test decides that."
            ),
            "construction": (
                "membrane-limited uptake cannot remove extracellular solute "
                "faster than P * (S/V). At the ENGINE's own geometry, "
                "S/V_tissue = 3*phi/a. Report BOTH conventions: "
                "k_ceiling_tissue = P * 3*phi/a and "
                "k_ceiling_interstitial = k_ceiling_tissue/(1-phi). The "
                "interstitial one is the PERMISSIVE end and the verdict uses it."
            ),
            "compared_against": {
                "k_struct": "clearance + binding, UNCALIBRATED compile",
                "k_eff_alpha": "clearance + alpha*uptake_proxy, alpha = 22.2235",
            },
            "decision_rule_fixed_now": {
                "A_reframing_SUPPORTED": "k_ceiling_interstitial >= k_eff_alpha",
                "B_partial_overlap": "k_struct <= k_ceiling_interstitial < k_eff_alpha",
                "C_reframing_REFUTED": "k_ceiling_interstitial < k_struct",
            },
            "predicted_outcome": "B",
            "predicted_band_k_ceiling_interstitial_per_s": [8e-4, 1.5e-3],
            "predicted_shortfall_vs_k_eff_alpha": "3x to 7x BELOW",
            "registered_reading_of_outcome_B": (
                "Under B, alpha asserts a perivascular sink ABOVE the ceiling "
                "that a MEASURED doxorubicin cell permeability can sustain. A "
                "binding CAPACITY cannot rescue that: capacity sets how long a "
                "flux can be absorbed, permeability sets the flux. So B is NOT "
                "'the same hole counted twice'. It is a NEW, independent, "
                "measured UPPER BOUND that alpha violates, and it is a "
                "different finding from the one the scouting pass asserted."
            ),
            "sensitivity_arm": (
                "repeat the ceiling with P_MDCK_MICROVESSEL = 32e-4 um/s, the "
                "highest measured DOX cell permeability Degerstedt cite (their "
                "ref 65, NOT verified at its own source), to establish whether "
                "the verdict's DIRECTION survives the top of the measured range"
            ),
        },
        "TEST_5_the_transvascular_term": {
            "status": "definitional correction, established by reading code, not a test",
            "finding_at_freeze": (
                "tmesim/treatment/transvascular.py states in its own docstring "
                "that the repository has NEITHER factor of PS = P*(S/V): the "
                "shipped entry term is PS = vascular_density * bbb_penetrance * "
                "PS_SCALE_S_INV with PS_SCALE_S_INV = 1 s^-1 undocumented, and "
                "no measured vessel radius exists. There is therefore NO "
                "transvascular wall permeability in the engine to score."
            ),
            "and": (
                "Degerstedt's P_cell is a CELL permeability, not a vessel-wall "
                "permeability, so it would be the wrong anchor for one even if "
                "the term existed. The task's 'transvascular permeability term' "
                "framing is corrected to 'cell-membrane permeability term'."
            ),
        },
        "what_this_cannot_establish_whatever_the_numbers": [
            "HepG2 hepatoma in a collagen/fibrinogen gel is not GBM in brain; "
            "nothing here transfers a doxorubicin number to another compound or "
            "another tissue.",
            "P_cell is APPARENT: it lumps active transport, binding and "
            "partition, and doxorubicin is a P-gp substrate and a DNA "
            "intercalator. The engine's P_mem is a PASSIVE permeability. These "
            "are not strictly the same object.",
            "The authors state the fluorescence self-quenches above ~25 uM and "
            "that they are 'more likely to underestimate rather than "
            "overestimate' P_cell, so the measurement is a LOWER bound and every "
            "ceiling derived from it is a lower bound on a ceiling.",
            "No absolute output of the simulator is validated or invalidated by "
            "this; the subject is one shipped parameter and one derived rate.",
        ],
    }

    with open(PREREG_PATH, "w", encoding="utf-8") as fh:
        json.dump(prereg, fh, indent=2, sort_keys=False)
        fh.write("\n")
    _guard()
    print(f"pre-registration FROZEN -> {PREREG_PATH}")
    print(f"sha256 = {_sha256(PREREG_PATH)}")


# --------------------------------------------------------------------------- #
# Stage 2: score.
# --------------------------------------------------------------------------- #
def stage_score() -> Dict:
    alpha = _guard()
    if not os.path.exists(PREREG_PATH):
        raise SystemExit(
            "no pre-registration found. Run --stage prereg first. Post-hoc is "
            "no pass at all here."
        )
    prereg_sha = _sha256(PREREG_PATH)
    print(f"pre-registration sha256 = {prereg_sha}\n")

    mol = library.get_molecule("doxorubicin")
    desc = mol.descriptors
    radius_nm = float(desc["radius_hydrodynamic_nm"])
    logp = float(desc["logp"])
    hbd = float(desc["hbd"])

    # --- the engine, as shipped ------------------------------------------- #
    prof_on = compile_drug.build_profile(mol, two_route=True)
    prof_off = compile_drug.build_profile(mol, two_route=False)
    interstitial_D = prof_off.effective_diffusivity  # the two_route input

    shipped = two_route.two_route_diffusivity(interstitial_D, radius_nm, logp, hbd)
    p_engine_two_route = float(shipped["membrane_permeability_um_s"])
    p_engine_perm_py = float(prof_on.membrane_permeability)
    phi = float(shipped["cell_volume_fraction"])
    a_um = two_route.DEFAULT_CELL_RADIUS_UM

    # --- TEST 1 ------------------------------------------------------------ #
    r1 = P_CELL_HD_UM_S / p_engine_two_route
    t1_pass = (r1 > 1.0) and (400.0 <= r1 <= 2000.0)

    # --- TEST 2 ------------------------------------------------------------ #
    r2 = p_engine_perm_py / P_CELL_HD_UM_S
    t2_pass = r2 > 1e4
    straddle = (r1 > 1.0) and (r2 > 1.0)
    internal_spread = p_engine_perm_py / p_engine_two_route

    # --- TEST 3: probe the measured P through the shipped law -------------- #
    probe = two_route.two_route_factor(
        cell_volume_fraction=phi,
        interstitial_D_um2s=interstitial_D,
        cytoplasm_D_um2s=float(shipped["cytoplasm_D_um2s"]),
        membrane_permeability_um_s=P_CELL_HD_UM_S,
        cell_radius_um=a_um,
    )
    ceiling = two_route.two_route_factor(
        cell_volume_fraction=phi,
        interstitial_D_um2s=interstitial_D,
        cytoplasm_D_um2s=float(shipped["cytoplasm_D_um2s"]),
        membrane_permeability_um_s=float("inf"),
        cell_radius_um=a_um,
    )
    d_shipped = float(shipped["effective_D_um2s"])
    d_probe = float(probe["effective_D_um2s"])
    d_ceiling = float(ceiling["effective_D_um2s"])
    rel_change = (d_probe - d_shipped) / d_shipped
    t3_pass = abs(rel_change) < 0.02

    # --- TEST 4: the alpha test ------------------------------------------- #
    dp = C._compile("doxorubicin")
    k_struct = float(dp.clearance) + float(dp.binding)
    k_eff_alpha = C.k_eff(dp, alpha)
    s_over_v = 3.0 * phi / a_um  # 1/um

    def ceil_rates(p_um_s: float) -> Dict[str, float]:
        kt = p_um_s * s_over_v
        return {
            "P_um_s": p_um_s,
            "k_ceiling_tissue_per_s": kt,
            "k_ceiling_interstitial_per_s": kt / (1.0 - phi),
        }

    ceil_measured = ceil_rates(P_CELL_HD_UM_S)
    ceil_measured_ld = ceil_rates(P_CELL_LD_UM_S)
    ceil_mdck = ceil_rates(P_MDCK_MICROVESSEL_UM_S)
    ceil_engine = ceil_rates(p_engine_two_route)

    k_ceil = ceil_measured["k_ceiling_interstitial_per_s"]
    if k_ceil >= k_eff_alpha:
        t4_outcome = "A_reframing_SUPPORTED"
    elif k_ceil >= k_struct:
        t4_outcome = "B_partial_overlap"
    else:
        t4_outcome = "C_reframing_REFUTED"
    t4_predicted = "B"
    t4_as_predicted = t4_outcome.startswith(t4_predicted)
    t4_in_band = 8e-4 <= k_ceil <= 1.5e-3

    # independent cross-check that Eq.2 reproduces the published P_cell
    # with a HepG2-scale radius: P = (k_cell/60) * r/3
    implied_r_um = P_CELL_HD_UM_S * 60.0 * 3.0 / K_CELL_HD_PER_MIN

    out = {
        "prereg_sha256": prereg_sha,
        "alpha": alpha,
        "doxorubicin_descriptors": {
            "mw": desc["mw"], "logp": logp, "hbd": hbd,
            "tpsa": desc["tpsa"], "radius_hydrodynamic_nm": radius_nm,
        },
        "engine_geometry": {"phi_cell": phi, "cell_radius_um": a_um,
                            "S_over_V_per_um": s_over_v},
        "engine_permeability_terms_um_s": {
            "two_route.membrane_permeability_um_s": p_engine_two_route,
            "two_route_log_k_membrane": float(shipped["log_k_membrane"]),
            "permeability.membrane_permeability": p_engine_perm_py,
            "internal_spread_fold": internal_spread,
        },
        "measurement_um_s": {
            "P_cell_HD_HepG2": P_CELL_HD_UM_S,
            "P_cell_LD_HepG2": P_CELL_LD_UM_S,
            "k_cell_HD_per_min": K_CELL_HD_PER_MIN,
            "radius_um_implied_by_eq2": implied_r_um,
            "P_MDCK_microvessel_unverified_at_source": P_MDCK_MICROVESSEL_UM_S,
        },
        "TEST_0": {"verdict": "PASS", "decided": "before computation, by reading"},
        "TEST_1": {"R1_measured_over_engine": r1, "pass": t1_pass,
                   "predicted_band": [400.0, 2000.0]},
        "TEST_2": {"R2_engine_over_measured": r2, "pass": t2_pass,
                   "straddle": straddle},
        "TEST_3": {
            "D_eff_shipped_um2s": d_shipped,
            "D_eff_at_measured_P_um2s": d_probe,
            "D_eff_at_infinite_P_um2s": d_ceiling,
            "relative_change": rel_change,
            "transcellular_gain_shipped": float(shipped["transcellular_gain"]),
            "transcellular_gain_at_measured_P": float(probe["transcellular_gain"]),
            "membrane_resistance_ratio_shipped": float(shipped["membrane_resistance_ratio"]),
            "membrane_resistance_ratio_at_measured_P": float(probe["membrane_resistance_ratio"]),
            "pass": t3_pass,
        },
        "TEST_4": {
            "k_struct_per_s": k_struct,
            "k_eff_alpha_per_s": k_eff_alpha,
            "alpha_amplification_on_k": k_eff_alpha / k_struct,
            "ceiling_from_engine_own_P": ceil_engine,
            "ceiling_from_measured_P_HD": ceil_measured,
            "ceiling_from_measured_P_LD": ceil_measured_ld,
            "ceiling_from_MDCK_unverified": ceil_mdck,
            "outcome": t4_outcome,
            "predicted_outcome": t4_predicted,
            "as_predicted": t4_as_predicted,
            "predicted_band_hit": t4_in_band,
            "k_eff_alpha_over_ceiling_interstitial": k_eff_alpha / k_ceil,
            "k_eff_alpha_over_ceiling_tissue":
                k_eff_alpha / ceil_measured["k_ceiling_tissue_per_s"],
            "k_eff_alpha_over_ceiling_MDCK_interstitial":
                k_eff_alpha / ceil_mdck["k_ceiling_interstitial_per_s"],
            "engine_own_binding_over_its_own_P_ceiling_interstitial":
                float(dp.binding) / ceil_engine["k_ceiling_interstitial_per_s"],
        },
    }

    print(json.dumps(out, indent=2))
    with open(SCORE_PATH, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
        fh.write("\n")
    print(f"\nwrote {SCORE_PATH}")
    _guard()
    print(f"alpha OK {alpha}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["prereg", "score"], required=True)
    args = ap.parse_args()
    if args.stage == "prereg":
        stage_prereg()
    else:
        stage_score()


if __name__ == "__main__":
    main()
