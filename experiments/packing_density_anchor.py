#!/usr/bin/env python3
"""Score the cell-volume-fraction term `phi` against a measured packing-density
perturbation — Grantab & Tannock 2012, BMC Cancer 12:214 (PMID 22672469,
PMC3407510, open access).

ROUTE R4 of results/verification_routes_2026-09-20.md.

WHAT IS UNDER TEST
------------------
`two_route.DEFAULT_CELL_VOLUME_FRACTION = 0.7892` enters EVERY compiled compound's
effective diffusivity through the Maxwell/Hasselman-Johnson hindrance

    hindrance = (1 + 2*beta*phi) / (1 - beta*phi)

and, in the impermeable-cell limit that all three of this paper's compounds occupy,

    hindrance = (1 - phi) / (1 + phi/2).

0.7892 was obtained by Delesse's principle on THIS REPOSITORY'S OWN Xenium breast
section. That makes it self-consistent, not verified. It has never been scored
against an external measurement.

WHAT THE PAPER SUPPLIES
-----------------------
Two perturbations of packing density (EGTA, bortezomib), measured as % nuclear
area before and after, with the matched change in drug penetration for the SAME
drug at the SAME dose through the SAME multicellular culture. The reported fold is

    C2(6 h)_treated / C2(6 h)_control

which is dimensionless and within-drug, so dose asymmetry, detectability and
absolute calibration cancel identically. That is what keeps this route clear of the
cross-compound ranking withdrawal of 2026-08-21.

PRE-REGISTRATION
----------------
results/packing_density_anchor.md section 0 was written to disk BEFORE this script
was run. It fixes the mapping, the scored observable, the criteria P1/P2/P3 and the
sensitivities S1-S7, and discloses exactly what had already been computed when the
thresholds were set. Nothing in this file may move a threshold.

NO PARAMETER IS TUNED HERE. The only engine value used is the shipped 0.7892, and
the only thing done to it is to scale it by a ratio the paper measured.

Usage:  .venv/bin/python experiments/packing_density_anchor.py
"""
from __future__ import annotations

import json
import math
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.linalg import solve_banded

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from tmesim.molecular import compile_drug, stokes_einstein, two_route  # noqa: E402
from tmesim.molecular.types import Molecule  # noqa: E402
from tmesim.validation.calibration import calibrate_alpha  # noqa: E402

ALPHA_FROZEN = 22.223477711172563
PHI_SHIPPED = 0.7892

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_JSON = os.path.join(REPO, "results", "packing_density_anchor_data.json")

PMC_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC3407510/"


def _hr(ch: str = "-") -> None:
    print(ch * 78)


# --------------------------------------------------------------------------- #
# THE MEASUREMENT. Every number below was read by me out of the PMC full text at
# PMC_URL on 2026-09-20 (HTTP 200, 195,531 bytes). The sentence each came from is
# quoted in results/packing_density_anchor.md section 0.3 and is reproduced in the
# `sentence` field so the code and the document cannot drift apart.
# --------------------------------------------------------------------------- #

PACKING = {
    # perturbation, cell line -> (NA_control_pct, NA_treated_pct, source sentence)
    ("EGTA", "HCT-8Ea"): (
        50.0, 40.0,
        "at this dose mean packing density (+/-SD) was reduced from 50 +/- 7.3% to "
        "40 +/- 1.6% and from 40 +/- 6.8% to 31 +/- 4.5%, in MCC derived from "
        "HCT-8Ea and HCT-8E11 cells respectively",
    ),
    ("EGTA", "HCT-8E11"): (
        40.0, 31.0,
        "at this dose mean packing density (+/-SD) was reduced from 50 +/- 7.3% to "
        "40 +/- 1.6% and from 40 +/- 6.8% to 31 +/- 4.5%, in MCC derived from "
        "HCT-8Ea and HCT-8E11 cells respectively",
    ),
    ("bortezomib", "HCT-8Ea"): (
        53.2, 44.6,
        "Table 1, 'MCC Packing Density (% nuclear area)', HCT-8Ea: 0 uM = 53.2 +/- 5.5, "
        "1.0 uM = 44.6 +/- 5.7*",
    ),
    ("bortezomib", "HCT-8E11"): (
        50.8, 42.5,
        "Table 1, 'MCC Packing Density (% nuclear area)', HCT-8E11: 0 uM = 50.8 +/- 3.8, "
        "1.0 uM = 42.5 +/- 3.8*",
    ),
}

# Standard deviations on packing density, same sentences/table. Used only for the
# propagated-uncertainty report, never for scoring.
PACKING_SD = {
    ("EGTA", "HCT-8Ea"): (7.3, 1.6),
    ("EGTA", "HCT-8E11"): (6.8, 4.5),
    ("bortezomib", "HCT-8Ea"): (5.5, 5.7),
    ("bortezomib", "HCT-8E11"): (3.8, 3.8),
}

#: The eight distinct reported fold-changes. `lo`/`hi` bracket the printed range;
#: a single printed value has lo == hi. `primary` marks the n=8 scored set.
MEASURED: List[dict] = [
    dict(row=1, pert="EGTA", line="HCT-8Ea", drug="5-fluorouracil",
         lo=3.3, hi=4.4, primary=True, lumped="printed once for 5-FU and gemcitabine",
         sentence="EGTA pre-treatment improved penetration of 5-fluorouracil and "
                  "gemcitabine (measured at 6 h) by ~3.3-4.4 fold in MCC derived from "
                  "HCT-8Ea cells"),
    dict(row=2, pert="EGTA", line="HCT-8Ea", drug="gemcitabine",
         lo=3.3, hi=4.4, primary=False, lumped="same printed value as row 1",
         sentence="EGTA pre-treatment improved penetration of 5-fluorouracil and "
                  "gemcitabine (measured at 6 h) by ~3.3-4.4 fold in MCC derived from "
                  "HCT-8Ea cells"),
    dict(row=3, pert="EGTA", line="HCT-8E11", drug="5-fluorouracil",
         lo=1.5, hi=1.9, primary=True, lumped="printed once for 5-FU and gemcitabine",
         sentence="and by ~1.5-1.9 fold in MCC derived from HCT-8E11 cells"),
    dict(row=4, pert="EGTA", line="HCT-8E11", drug="gemcitabine",
         lo=1.5, hi=1.9, primary=False, lumped="same printed value as row 3",
         sentence="and by ~1.5-1.9 fold in MCC derived from HCT-8E11 cells"),
    dict(row=5, pert="EGTA", line="HCT-8Ea", drug="doxorubicin",
         lo=2.0, hi=2.0, primary=True, lumped="printed once for both cell lines",
         note="data not shown; no sucrose internal standard in doxorubicin runs",
         sentence="EGTA also led to ~2-fold increase in penetration of doxorubicin "
                  "(p = 0.02 for HCT-8E11, p = 0.002 for HCT-8Ea; data not shown)"),
    dict(row=6, pert="EGTA", line="HCT-8E11", drug="doxorubicin",
         lo=2.0, hi=2.0, primary=True, lumped="printed once for both cell lines",
         note="data not shown; no sucrose internal standard in doxorubicin runs",
         sentence="EGTA also led to ~2-fold increase in penetration of doxorubicin "
                  "(p = 0.02 for HCT-8E11, p = 0.002 for HCT-8Ea; data not shown)"),
    dict(row=7, pert="bortezomib", line="HCT-8Ea", drug="5-fluorouracil",
         lo=1.7, hi=1.7, primary=True, lumped=None,
         sentence="Bortezomib pre-treatment improved penetration of 5-fluorouracil by "
                  "~1.7-fold and gemcitabine by ~3-fold in MCC derived from HCT-8Ea cells"),
    dict(row=8, pert="bortezomib", line="HCT-8Ea", drug="gemcitabine",
         lo=3.0, hi=3.0, primary=True, lumped=None,
         sentence="Bortezomib pre-treatment improved penetration of 5-fluorouracil by "
                  "~1.7-fold and gemcitabine by ~3-fold in MCC derived from HCT-8Ea cells"),
    dict(row=9, pert="bortezomib", line="HCT-8E11", drug="5-fluorouracil",
         lo=2.0, hi=2.0, primary=True, lumped="printed once for both drugs",
         sentence="and increased penetration of both drugs by ~2-fold in MCC derived "
                  "from HCT-8E11 cells at 6 h"),
    dict(row=10, pert="bortezomib", line="HCT-8E11", drug="gemcitabine",
         lo=2.0, hi=2.0, primary=True, lumped="printed once for both drugs",
         sentence="and increased penetration of both drugs by ~2-fold in MCC derived "
                  "from HCT-8E11 cells at 6 h"),
    dict(row=11, pert="bortezomib", line="HCT-8Ea", drug="doxorubicin",
         lo=1.8, hi=1.8, primary=False, lumped="printed once for both cell lines",
         note="data not shown; no sucrose internal standard in doxorubicin runs",
         sentence="Penetration of doxorubicin was also increased by about 1.8-fold at "
                  "6 h (p = 0.02 for HCT-8Ea, p = 0.005 for HCT-8E11, data not shown)"),
    dict(row=12, pert="bortezomib", line="HCT-8E11", drug="doxorubicin",
         lo=1.8, hi=1.8, primary=False, lumped="printed once for both cell lines",
         note="data not shown; no sucrose internal standard in doxorubicin runs",
         sentence="Penetration of doxorubicin was also increased by about 1.8-fold at "
                  "6 h (p = 0.02 for HCT-8Ea, p = 0.005 for HCT-8E11, data not shown)"),
]

#: Descriptors. logP for 5-FU and doxorubicin are the engine's OWN experimental
#: entries (lipophilicity.EXPERIMENTAL_LOGP, Hansch 1995 via PubChem/HSDB).
#: Gemcitabine is NOT in the engine's library and has no experimental entry; the
#: value below is a placeholder whose irrelevance is DEMONSTRATED by sensitivity S6
#: (all three compounds sit at the impermeable limit, where logP does not enter).
DRUGS = {
    "5-fluorouracil": dict(mw=130.077, logp=-0.89, hbd=2,
                           logp_src="lipophilicity.EXPERIMENTAL_LOGP, Hansch 1995 "
                                    "via PubChem/HSDB 3228"),
    "doxorubicin": dict(mw=543.52, logp=1.27, hbd=6,
                        logp_src="lipophilicity.EXPERIMENTAL_LOGP, Hansch 1995"),
    "gemcitabine": dict(mw=263.20, logp=-1.4, hbd=4,
                        logp_src="NOT an engine value and NOT scored: gemcitabine has "
                                 "no entry in lipophilicity.EXPERIMENTAL_LOGP. "
                                 "Irrelevance demonstrated by sensitivity S6."),
}

MEASUREMENT_NOISE_FLOOR = math.sqrt(2.0) * 0.20  # from the paper's own +/-20% sucrose gate


# --------------------------------------------------------------------------- #
# ENGINE
# --------------------------------------------------------------------------- #

_D_CACHE: Dict[float, float] = {}


def engine_interstitial_D(mw: float) -> float:
    """The shipped interstitial (matrix) diffusivity for a molecular weight.

    `two_route=False` is REQUIRED: since results/two_route_promotion.md
    `build_profile` applies the two-route term itself, so the default would return a
    D that already carries it and the law would be applied twice. Same call as
    experiments/two_route_transport.py::_engine_interstitial_D.

    Memoised on MW only. The call is deterministic and depends on nothing else
    (phi does NOT enter it -- the two-route term is applied afterwards), so the
    cache is exact, not an approximation. Same caching pattern as
    experiments/two_route_transport.py::_build_rows.
    """
    mw = float(mw)
    if mw not in _D_CACHE:
        _D_CACHE[mw] = compile_drug.build_profile(
            Molecule(name=f"mw{mw}", smiles="", descriptors={"mw": mw}),
            two_route=False,
        ).effective_diffusivity
    return _D_CACHE[mw]


def hindrance(phi: float, drug: str, logp_override: Optional[float] = None) -> dict:
    """The shipped two-route hindrance at a given phi, through the shipped path."""
    d = DRUGS[drug]
    logp = d["logp"] if logp_override is None else logp_override
    D = engine_interstitial_D(d["mw"])
    r = stokes_einstein.hydrodynamic_radius_nm(d["mw"])
    return two_route.two_route_diffusivity(
        D, r, logp, d["hbd"], cell_volume_fraction=phi
    )


def phi_treated(phi_control: float, na_ctrl: float, na_trt: float) -> float:
    """PRE-REGISTERED mapping (results/packing_density_anchor.md 0.4).

    Nuclear volume fraction nu = f_N * phi. If neither nuclear nor cell size changes,
    f_N cancels from the ratio and phi_trt/phi_ctrl = NA_trt/NA_ctrl. The paper
    explicitly did not verify constant nuclear size ("without taking into account
    differences in nuclear size"), and that is carried as a limit, not buried.
    """
    return phi_control * (na_trt / na_ctrl)


def predicted_ratio(drug: str, pert: str, line: str, phi_control: float = PHI_SHIPPED,
                    logp_override: Optional[float] = None) -> float:
    na_c, na_t, _ = PACKING[(pert, line)]
    phi_t = phi_treated(phi_control, na_c, na_t)
    h_c = hindrance(phi_control, drug, logp_override)["hindrance"]
    h_t = hindrance(phi_t, drug, logp_override)["hindrance"]
    return h_t / h_c


# --------------------------------------------------------------------------- #
# SCORING
# --------------------------------------------------------------------------- #

def gmfe(pred: List[float], meas: List[float]) -> float:
    return math.exp(sum(abs(math.log(p / m)) for p, m in zip(pred, meas)) / len(pred))


def geo_mid(lo: float, hi: float) -> float:
    return math.sqrt(lo * hi)


def spread(xs: List[float]) -> float:
    return max(xs) / min(xs)


# --------------------------------------------------------------------------- #
# S3 -- the UNSCORED transient arm.
#
# A 1-D slab between a donor and a stirred receiver. The paper reports NEITHER the
# slab thickness NOR the membrane area, so this arm has free parameters and is
# therefore reported as a band and NEVER scored. See section 0.6 of the document:
# using it as the scored arm would be indistinguishable from tuning to the anchor.
# --------------------------------------------------------------------------- #

def two_compartment_receiver(D_um2_s: float, L_um: float, area_cm2: float,
                             v1_ml: float, v2_ml: float, t_end_s: float,
                             n_nodes: int = 120, dt_s: float = 5.0,
                             retardation: float = 1.0) -> float:
    """C2(t_end)/C_inf for a slab between a well-mixed donor and a stirred receiver.

    Backward Euler on the slab, explicit coupling of the two well-mixed compartments
    through the end-face fluxes. Verified in `verify_transient_solver()` against the
    quasi-steady analytic limit and against global mass conservation.

    Units: um, s, um^2/s. Volumes are converted from mL, area from cm^2.
    """
    n = int(n_nodes)
    dx = L_um / (n - 1)
    a = D_um2_s / retardation
    lam = a * dt_s / dx ** 2

    A_um2 = area_cm2 * 1.0e8          # 1 cm^2 = 1e8 um^2
    V1 = v1_ml * 1.0e12               # 1 mL = 1 cm^3 = 1e12 um^3
    V2 = v2_ml * 1.0e12

    c = np.zeros(n)                   # slab, starts drug-free
    c1, c2 = 1.0, 0.0                 # donor at unit concentration, receiver empty

    ab = np.zeros((3, n))
    ab[0, 1:] = -lam
    ab[1, :] = 1.0 + 2.0 * lam
    ab[2, :-1] = -lam
    # Dirichlet rows at both faces
    ab[1, 0] = 1.0
    ab[0, 1] = 0.0
    ab[1, -1] = 1.0
    ab[2, -2] = 0.0

    steps = int(round(t_end_s / dt_s))
    for _ in range(steps):
        rhs = c.copy()
        rhs[0] = c1
        rhs[-1] = c2
        c = solve_banded((1, 1), ab, rhs)
        # one-sided end-face gradients -> fluxes -> compartment updates
        j_in = -D_um2_s * (c[1] - c[0]) / dx      # into the slab at x=0
        j_out = -D_um2_s * (c[-1] - c[-2]) / dx   # out of the slab at x=L
        c1 -= A_um2 * j_in * dt_s / V1
        c2 += A_um2 * j_out * dt_s / V2

    c_inf = v1_ml / (v1_ml + v2_ml)   # equilibrium between the two compartments
    return c2 / c_inf


def verify_transient_solver() -> dict:
    """Verification of the transient harness. Not a result about the engine.

    (i) Thin-slab quasi-steady limit: when the slab's own capacity is negligible the
        approach to equilibrium is a single exponential with rate
        k = (A*D/L) * (1/V1 + 1/V2), so C2/C_inf = 1 - exp(-k t) exactly.
    (ii) Grid independence: halving dx and dt must not move the answer materially.
    """
    D, L, A, V1, V2 = 50.0, 20.0, 1.0, 0.5, 18.0
    t = 3600.0
    A_um2, V1u, V2u = A * 1e8, V1 * 1e12, V2 * 1e12
    k = (A_um2 * D / L) * (1.0 / V1u + 1.0 / V2u)
    analytic = 1.0 - math.exp(-k * t)
    numeric = two_compartment_receiver(D, L, A, V1, V2, t, n_nodes=200, dt_s=1.0)
    err = abs(numeric - analytic) / analytic

    coarse = two_compartment_receiver(50.0, 150.0, 1.0, 0.5, 18.0, 21600.0,
                                      n_nodes=80, dt_s=10.0)
    fine = two_compartment_receiver(50.0, 150.0, 1.0, 0.5, 18.0, 21600.0,
                                    n_nodes=240, dt_s=2.5)
    grid = abs(fine - coarse) / fine
    return {"quasi_steady_analytic": analytic, "quasi_steady_numeric": numeric,
            "quasi_steady_rel_err": err, "grid_coarse": coarse, "grid_fine": fine,
            "grid_rel_diff": grid}


# --------------------------------------------------------------------------- #

def main() -> int:
    alpha_start = calibrate_alpha()
    assert alpha_start == ALPHA_FROZEN, f"calibrate_alpha START {alpha_start!r}"
    assert two_route.DEFAULT_CELL_VOLUME_FRACTION == PHI_SHIPPED, "phi moved"
    print(f"calibrate_alpha START = {alpha_start!r}  (asserted == {ALPHA_FROZEN!r})")
    print(f"two_route.DEFAULT_CELL_VOLUME_FRACTION = {PHI_SHIPPED}")
    print(f"paper: Grantab & Tannock 2012, BMC Cancer 12:214, PMC3407510")
    print(f"       {PMC_URL}")

    out: Dict[str, object] = {
        "route": "R4 of results/verification_routes_2026-09-20.md",
        "date": "2026-09-20",
        "paper": {"pmid": 22672469, "pmc": "PMC3407510",
                  "doi": "10.1186/1471-2407-12-214", "url": PMC_URL},
        "calibrate_alpha_start": alpha_start,
        "phi_shipped": PHI_SHIPPED,
        "measurement_noise_floor_rel": MEASUREMENT_NOISE_FLOOR,
    }

    # ---------------- compound compile: establish the impermeable limit -------
    _hr("=")
    print("STAGE 1 -- COMPOUND COMPILE THROUGH THE SHIPPED PATH (phi = 0.7892)")
    _hr("=")
    print(f"  {'compound':<16}{'MW':>8}{'logP':>7}{'D_int':>9}{'P(um/s)':>11}"
          f"{'memb.resist':>13}{'hindrance':>12}{'imperm':>12}{'gain':>11}")
    compile_rows = []
    for name in ("5-fluorouracil", "doxorubicin", "gemcitabine"):
        h = hindrance(PHI_SHIPPED, name)
        row = dict(drug=name, mw=DRUGS[name]["mw"], logp=DRUGS[name]["logp"],
                   logp_source=DRUGS[name]["logp_src"],
                   interstitial_D_um2s=h["sigma_matrix_um2s"],
                   membrane_permeability_um_s=h["membrane_permeability_um_s"],
                   membrane_resistance_ratio=h["membrane_resistance_ratio"],
                   hindrance=h["hindrance"],
                   impermeable_hindrance=h["impermeable_hindrance"],
                   transcellular_gain=h["transcellular_gain"])
        compile_rows.append(row)
        print(f"  {name:<16}{row['mw']:>8.2f}{row['logp']:>7.2f}"
              f"{row['interstitial_D_um2s']:>9.1f}"
              f"{row['membrane_permeability_um_s']:>11.3g}"
              f"{row['membrane_resistance_ratio']:>13.3g}"
              f"{row['hindrance']:>12.8f}{row['impermeable_hindrance']:>12.8f}"
              f"{row['transcellular_gain']:>11.7f}")
    out["compile"] = compile_rows
    max_gain = max(r["transcellular_gain"] for r in compile_rows)
    print(f"\n  All three sit at the impermeable-cell limit: the largest transcellular")
    print(f"  gain is {max_gain:.7f}, i.e. the transcellular route contributes at most")
    print(f"  {100*(max_gain-1):.4f}% of the hindrance. The predicted RATIO is therefore")
    print(f"  a function of phi alone and is the SAME for all three compounds.")

    # ---------------- the phi mapping ----------------------------------------
    _hr("=")
    print("STAGE 2 -- THE PRE-REGISTERED phi MAPPING")
    _hr("=")
    print(f"  {'perturbation':<13}{'line':<11}{'NA_ctrl':>9}{'NA_trt':>9}"
          f"{'NA ratio':>10}{'phi_ctrl':>10}{'phi_trt':>10}")
    map_rows = []
    for (pert, line), (na_c, na_t, sent) in PACKING.items():
        sd_c, sd_t = PACKING_SD[(pert, line)]
        pt = phi_treated(PHI_SHIPPED, na_c, na_t)
        mr = dict(pert=pert, line=line, na_control_pct=na_c, na_treated_pct=na_t,
                  na_control_sd=sd_c, na_treated_sd=sd_t,
                  na_ratio=na_t / na_c, phi_control=PHI_SHIPPED, phi_treated=pt,
                  sentence=sent)
        map_rows.append(mr)
        print(f"  {pert:<13}{line:<11}{na_c:>9.1f}{na_t:>9.1f}{na_t/na_c:>10.5f}"
              f"{PHI_SHIPPED:>10.4f}{pt:>10.6f}")
    out["mapping"] = map_rows
    na_ratios = [m["na_ratio"] for m in map_rows]
    print(f"\n  The four NA ratios span {min(na_ratios):.4f}-{max(na_ratios):.4f}, "
          f"a range of {spread(na_ratios):.3f}x.")
    print(f"  Every phi_treated is inside [0,1): admissible.")

    # ---------------- the scored comparison ----------------------------------
    _hr("=")
    print("STAGE 3 -- PREDICTED vs MEASURED (scored arm: the diffusivity ratio r_D)")
    _hr("=")
    print(f"  {'#':>3} {'pert':<11}{'line':<11}{'drug':<16}{'meas':>10}"
          f"{'pred':>9}{'pred/meas':>11}  {'primary':<8}")
    table = []
    for m in MEASURED:
        pred = predicted_ratio(m["drug"], m["pert"], m["line"])
        mid = geo_mid(m["lo"], m["hi"])
        meas_s = f"{m['lo']:.2f}" if m["lo"] == m["hi"] else f"{m['lo']:.1f}-{m['hi']:.1f}"
        r = dict(row=m["row"], pert=m["pert"], line=m["line"], drug=m["drug"],
                 measured_lo=m["lo"], measured_hi=m["hi"], measured_geomid=mid,
                 predicted=pred, fold_error=pred / mid, primary=m["primary"],
                 lumped=m.get("lumped"), note=m.get("note"), sentence=m["sentence"])
        table.append(r)
        print(f"  {m['row']:>3} {m['pert']:<11}{m['line']:<11}{m['drug']:<16}"
              f"{meas_s:>10}{pred:>9.4f}{pred/mid:>11.3f}  "
              f"{'YES' if m['primary'] else 'no':<8}")
    out["table"] = table

    prim = [r for r in table if r["primary"]]
    pred_p = [r["predicted"] for r in prim]
    meas_p = [r["measured_geomid"] for r in prim]
    n = len(prim)
    print(f"\n  n(primary) = {n}")

    # ---------------- criteria ------------------------------------------------
    _hr("=")
    print("STAGE 4 -- THE PRE-REGISTERED CRITERIA")
    _hr("=")

    g_eng = gmfe(pred_p, meas_p)
    g_null = gmfe([1.0] * n, meas_p)
    p1 = g_eng <= 2.0
    p2 = g_eng <= 0.8 * g_null
    s_meas, s_pred = spread(meas_p), spread(pred_p)
    explained = math.log(s_pred) / math.log(s_meas)
    p3 = explained >= 0.5
    p4 = all(p > 1.0 for p in pred_p)

    print(f"  P1  accuracy        GMFE(engine) = {g_eng:.4f}  bar <= 2.0        "
          f"-> {'PASS' if p1 else 'FAIL'}   [DEMOTED, see 0.1]")
    print(f"  P2  anti-vacuity    GMFE(engine) = {g_eng:.4f}  vs 0.8 x GMFE(null)")
    print(f"                      GMFE(null, phi-blind, ratio==1) = {g_null:.4f}; "
          f"bar <= {0.8*g_null:.4f}  -> {'PASS' if p2 else 'FAIL'}")
    print(f"  P3  dispersion      measured spread {s_meas:.4f}x, predicted spread "
          f"{s_pred:.4f}x")
    print(f"                      ln(S_pred)/ln(S_meas) = {explained:.4f}  bar >= 0.50"
          f"  -> {'PASS' if p3 else 'FAIL'}")
    print(f"  P4  direction       all predicted > 1: {p4}  "
          f"-> VACUOUS BY CONSTRUCTION, NOT COUNTED")

    verdict = "PASS" if (p1 and p2 and p3) else "FAIL"
    out["criteria"] = dict(
        n_primary=n, gmfe_engine=g_eng, gmfe_null_phi_blind=g_null,
        P1_accuracy=dict(value=g_eng, bar=2.0, pass_=p1, demoted=True),
        P2_anti_vacuity=dict(value=g_eng, bar=0.8 * g_null, gmfe_null=g_null, pass_=p2),
        P3_dispersion=dict(measured_spread=s_meas, predicted_spread=s_pred,
                           log_fraction_explained=explained, bar=0.5, pass_=p3),
        P4_direction=dict(value=p4, counted=False,
                          why="vacuous by construction: phi falls in every row and "
                              "hindrance is monotone decreasing in phi"),
        verdict=verdict,
    )
    _hr()
    print(f"  VERDICT (P1 and P2 and P3) = {verdict}")
    _hr()

    # ---------------- S1: implied phi_control --------------------------------
    _hr("=")
    print("STAGE 5 -- S1: WHAT phi_control WOULD REPRODUCE EACH MEASUREMENT?")
    _hr("=")
    print("  Invert the law per row. Admissible band is (NA_control, 1): a cell volume")
    print("  fraction cannot be below the nuclear volume fraction it contains.")
    print(f"  {'#':>3} {'pert':<11}{'line':<11}{'drug':<16}{'meas':>8}"
          f"{'phi_implied':>13}{'band':>16}{'admissible':>12}")
    s1 = []
    for r in table:
        na_c, na_t, _ = PACKING[(r["pert"], r["line"])]
        lo_band = na_c / 100.0
        target = r["measured_geomid"]
        drug = r["drug"]

        def f(phi: float) -> float:
            hc = hindrance(phi, drug)["hindrance"]
            ht = hindrance(phi_treated(phi, na_c, na_t), drug)["hindrance"]
            return ht / hc - target

        # r_D(phi_control) rises monotonically from 1 (phi -> 0) to unbounded
        # (phi -> 1), so a sign change brackets the unique root.
        lo, hi = 1e-6, 0.999999
        f_lo, f_hi = f(lo), f(hi)
        implied: Optional[float] = None
        if f_lo * f_hi <= 0:
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                f_mid = f(mid)
                if f_lo * f_mid <= 0:
                    hi, f_hi = mid, f_mid
                else:
                    lo, f_lo = mid, f_mid
            implied = 0.5 * (lo + hi)
        adm = implied is not None and lo_band < implied < 1.0
        s1.append(dict(row=r["row"], pert=r["pert"], line=r["line"], drug=drug,
                       measured=target, phi_implied=implied,
                       admissible_band=[lo_band, 1.0], admissible=bool(adm),
                       primary=r["primary"]))
        ps = f"{implied:.6f}" if implied is not None else "none in (0,1)"
        print(f"  {r['row']:>3} {r['pert']:<11}{r['line']:<11}{drug:<16}"
              f"{target:>8.3f}{ps:>13}"
              f"{f'({lo_band:.3f}, 1)':>16}{('yes' if adm else 'NO'):>12}")
    out["S1_implied_phi"] = s1
    got = [x["phi_implied"] for x in s1 if x["primary"] and x["phi_implied"] is not None]
    n_inadm = sum(1 for x in s1 if x["primary"] and not x["admissible"])
    if got:
        print(f"\n  primary rows with an implied phi in (0,1): {len(got)}/{n}; "
              f"range {min(got):.4f}-{max(got):.4f}")
    print(f"  primary rows with NO admissible phi_control at all: {n_inadm}/{n}")
    out["criteria"]["P5_admissibility"] = dict(
        n_primary=n, n_inadmissible=n_inadm,
        implied_phi_min=min(got) if got else None,
        implied_phi_max=max(got) if got else None)

    # ---------------- S2: the inadmissible identity mapping ------------------
    _hr("=")
    print("STAGE 6 -- S2: IDENTITY MAPPING phi_control = NA_control")
    _hr("=")
    print("  PHYSICALLY INADMISSIBLE (it asserts the cells are entirely nucleus).")
    print("  Computed only as the extreme low-sensitivity bound.")
    s2_pred, s2_meas = [], []
    print(f"  {'#':>3} {'pert':<11}{'line':<11}{'drug':<16}{'phi_ctrl':>10}"
          f"{'pred':>9}{'meas':>8}")
    s2_rows = []
    for r in table:
        na_c, _, _ = PACKING[(r["pert"], r["line"])]
        pc = na_c / 100.0
        p = predicted_ratio(r["drug"], r["pert"], r["line"], phi_control=pc)
        s2_rows.append(dict(row=r["row"], phi_control=pc, predicted=p,
                            measured=r["measured_geomid"], primary=r["primary"]))
        if r["primary"]:
            s2_pred.append(p)
            s2_meas.append(r["measured_geomid"])
        print(f"  {r['row']:>3} {r['pert']:<11}{r['line']:<11}{r['drug']:<16}"
              f"{pc:>10.4f}{p:>9.4f}{r['measured_geomid']:>8.3f}")
    g_s2 = gmfe(s2_pred, s2_meas)
    print(f"\n  GMFE under the identity mapping = {g_s2:.4f} "
          f"(shipped mapping {g_eng:.4f}, phi-blind null {g_null:.4f})")
    out["S2_identity_mapping"] = dict(rows=s2_rows, gmfe=g_s2,
                                      beats_null=bool(g_s2 <= 0.8 * g_null))

    # ---------------- S3: the unscored transient band ------------------------
    _hr("=")
    print("STAGE 7 -- S3: UNSCORED TRANSIENT BAND (geometry is NOT reported)")
    _hr("=")
    ver = verify_transient_solver()
    print(f"  harness verification:")
    print(f"    quasi-steady analytic {ver['quasi_steady_analytic']:.8f} vs numeric "
          f"{ver['quasi_steady_numeric']:.8f}  rel err "
          f"{ver['quasi_steady_rel_err']:.3e}")
    print(f"    grid halving: {ver['grid_coarse']:.8f} -> {ver['grid_fine']:.8f}  "
          f"rel diff {ver['grid_rel_diff']:.3e}")
    out["S3_solver_verification"] = ver

    print(f"\n  Sweep: slab thickness 40-400 um, area 1 cm^2, V1 = 0.5 mL, V2 = 18 mL,")
    print(f"  t = 6 h, D_control = the engine's 5-FU tissue D_eff at phi = 0.7892.")
    h5 = hindrance(PHI_SHIPPED, "5-fluorouracil")
    D_ctrl = h5["effective_D_um2s"]
    print(f"  D_control = {D_ctrl:.4f} um^2/s")
    print(f"  {'L (um)':>8}{'C2/Cinf ctrl':>15}{'C2/Cinf EGTA-Ea':>18}{'ratio':>9}"
          f"{'r_D':>9}{'amplification':>15}")
    r_D_egta_ea = predicted_ratio("5-fluorouracil", "EGTA", "HCT-8Ea")
    D_trt = D_ctrl * r_D_egta_ea
    s3 = []
    for L in (40.0, 60.0, 100.0, 150.0, 200.0, 300.0, 400.0, 600.0, 800.0,
              1200.0, 1600.0, 2400.0, 3200.0):
        cc = two_compartment_receiver(D_ctrl, L, 1.0, 0.5, 18.0, 6 * 3600.0,
                                      n_nodes=160, dt_s=5.0)
        ct = two_compartment_receiver(D_trt, L, 1.0, 0.5, 18.0, 6 * 3600.0,
                                      n_nodes=160, dt_s=5.0)
        ratio = ct / cc
        s3.append(dict(thickness_um=L, c2_over_cinf_control=cc,
                       c2_over_cinf_treated=ct, transient_ratio=ratio,
                       r_D=r_D_egta_ea, amplification=ratio / r_D_egta_ea))
        print(f"  {L:>8.0f}{cc:>15.6f}{ct:>18.6f}{ratio:>9.4f}"
              f"{r_D_egta_ea:>9.4f}{ratio/r_D_egta_ea:>15.4f}")
    out["S3_transient_band"] = dict(
        rows=s3, D_control_um2s=D_ctrl, r_D=r_D_egta_ea,
        note="UNSCORED. The paper reports neither slab thickness nor membrane area; "
             "this arm therefore has free parameters and could be made to produce "
             "almost any ratio. Reported as a band only.")
    amps = [x["amplification"] for x in s3]
    print(f"\n  amplification over the sweep: {min(amps):.4f}-{max(amps):.4f}x")
    reach = [x for x in s3 if x["transient_ratio"] >= 3.811]
    if reach:
        print(f"  thinnest slab at which the transient arm reaches the measured 3.811x:"
              f" {min(x['thickness_um'] for x in reach):.0f} um")
    else:
        print(f"  the transient arm NEVER reaches the measured 3.811x anywhere in the")
        print(f"  sweep (max {max(x['transient_ratio'] for x in s3):.4f}x at "
              f"{max(x['thickness_um'] for x in s3):.0f} um). The unreported geometry")
        print(f"  is NOT an escape route for the largest measured row.")
    out["S3_transient_band"]["reaches_largest_measured"] = bool(reach)

    # ---------------- POST-HOC D1: can this experiment see the LAW's form? ---
    _hr("=")
    print("STAGE 7b -- POST-HOC DIAGNOSTIC D1 (not a pre-registered criterion)")
    _hr("=")
    print("  Does this measurement discriminate the Maxwell/Hasselman-Johnson form")
    print("  from cruder laws with the same monotonicity? Three laws, same mapping,")
    print("  same shipped phi_control = 0.7892:")
    print("    HJ        (1-phi)/(1+phi/2)      -- the shipped law, impermeable limit")
    print("    exclusion (1-phi)                -- naive volume exclusion, no tortuosity")
    print("    MG2D      (1-phi)/(1+phi)        -- 2-D Maxwell-Garnett")
    laws = {
        "hasselman_johnson_shipped": lambda p: (1.0 - p) / (1.0 + p / 2.0),
        "naive_volume_exclusion": lambda p: (1.0 - p),
        "maxwell_garnett_2d": lambda p: (1.0 - p) / (1.0 + p),
    }
    d1 = {}
    print(f"  {'law':<28}{'GMFE':>9}{'spread':>9}{'ln-frac':>10}")
    for lname, fn in laws.items():
        pv = []
        for r in prim:
            na_c, na_t, _ = PACKING[(r["pert"], r["line"])]
            pt = phi_treated(PHI_SHIPPED, na_c, na_t)
            pv.append(fn(pt) / fn(PHI_SHIPPED))
        g = gmfe(pv, meas_p)
        sp = spread(pv)
        ex = math.log(sp) / math.log(s_meas)
        d1[lname] = dict(gmfe=g, predicted_spread=sp, log_fraction_explained=ex,
                         predicted=pv)
        print(f"  {lname:<28}{g:>9.4f}{sp:>9.4f}{ex:>10.4f}")
    out["D1_law_discrimination"] = dict(
        laws=d1,
        note="POST-HOC diagnostic, not a pre-registered criterion. It asks whether "
             "this measurement could distinguish the shipped functional form from "
             "cruder ones, i.e. what the test is capable of seeing.")

    # ---------------- POST-HOC D2: where the measured spread lives -----------
    _hr("=")
    print("STAGE 7c -- POST-HOC DIAGNOSTIC D2 (not a pre-registered criterion)")
    _hr("=")
    print("  The engine is a function of phi alone, so it must predict ONE ratio per")
    print("  (perturbation, cell line) regardless of compound. Where does the measured")
    print("  dispersion actually sit?")
    print(f"  {'perturbation':<13}{'line':<11}{'NA ratio':>10}{'pred':>9}"
          f"{'measured rows':>32}{'meas spread':>13}")
    d2 = []
    for (pert, line), (na_c, na_t, _) in PACKING.items():
        sub = [r for r in prim if r["pert"] == pert and r["line"] == line]
        if not sub:
            continue
        mv = [r["measured_geomid"] for r in sub]
        pv = [r["predicted"] for r in sub]
        lbl = ", ".join(f"{r['drug'][:4]} {r['measured_geomid']:.2f}" for r in sub)
        d2.append(dict(pert=pert, line=line, na_ratio=na_t / na_c,
                       predicted_spread=spread(pv), measured_spread=spread(mv),
                       predicted_mean=sum(pv) / len(pv), measured=mv))
        print(f"  {pert:<13}{line:<11}{na_t/na_c:>10.5f}{sum(pv)/len(pv):>9.4f}"
              f"{lbl:>32}{spread(mv):>13.4f}")
    within = max(x["measured_spread"] for x in d2)
    print(f"\n  Largest measured spread WITHIN a single (perturbation, cell line) cell,")
    print(f"  where the engine's predicted spread is {max(x['predicted_spread'] for x in d2):.4f}x:"
          f" {within:.4f}x")
    print(f"  GUARDRAIL: these are within-drug RESPONSE RATIOS (same drug, same dose,")
    print(f"  before vs after). Nothing here orders compounds by penetration or by")
    print(f"  deliverability, and no two compounds' concentrations are compared.")
    out["D2_dispersion_location"] = dict(
        cells=d2, max_within_cell_measured_spread=within,
        max_within_cell_predicted_spread=max(x["predicted_spread"] for x in d2),
        note="POST-HOC diagnostic. Within-drug response ratios only; this is a "
             "statement about the engine's dispersion, NOT a cross-compound ordering.")

    # ---------------- S4/S5/S6/S7 --------------------------------------------
    _hr("=")
    print("STAGE 8 -- S4/S5/S6/S7")
    _hr("=")

    # S4: drop doxorubicin (no sucrose internal standard, data not shown)
    nd = [r for r in prim if r["drug"] != "doxorubicin"]
    g_s4 = gmfe([r["predicted"] for r in nd], [r["measured_geomid"] for r in nd])
    gn_s4 = gmfe([1.0] * len(nd), [r["measured_geomid"] for r in nd])
    sp_s4 = spread([r["predicted"] for r in nd])
    sm_s4 = spread([r["measured_geomid"] for r in nd])
    ex_s4 = math.log(sp_s4) / math.log(sm_s4)
    print(f"  S4  doxorubicin rows dropped (n={len(nd)}): GMFE {g_s4:.4f} vs null "
          f"{gn_s4:.4f}; P2 {'PASS' if g_s4 <= 0.8*gn_s4 else 'FAIL'}; "
          f"ln-spread {ex_s4:.4f}, P3 {'PASS' if ex_s4 >= 0.5 else 'FAIL'}")
    out["S4_no_doxorubicin"] = dict(n=len(nd), gmfe=g_s4, gmfe_null=gn_s4,
                                    P2_pass=bool(g_s4 <= 0.8 * gn_s4),
                                    log_fraction_explained=ex_s4,
                                    P3_pass=bool(ex_s4 >= 0.5))

    # S5: range endpoints
    s5 = {}
    for end in ("lo", "hi"):
        mv = [r[f"measured_{end}"] for r in prim]
        g = gmfe(pred_p, mv)
        gn = gmfe([1.0] * n, mv)
        ex = math.log(spread(pred_p)) / math.log(spread(mv))
        s5[end] = dict(gmfe=g, gmfe_null=gn, P1_pass=bool(g <= 2.0),
                       P2_pass=bool(g <= 0.8 * gn), log_fraction_explained=ex,
                       P3_pass=bool(ex >= 0.5))
        print(f"  S5  measured at range {end:>2}: GMFE {g:.4f} (null {gn:.4f}), "
              f"P1 {'PASS' if g<=2.0 else 'FAIL'}, "
              f"P2 {'PASS' if g<=0.8*gn else 'FAIL'}, "
              f"ln-spread {ex:.4f} P3 {'PASS' if ex>=0.5 else 'FAIL'}")
    out["S5_range_endpoints"] = s5

    # S6: gemcitabine logP
    base = predicted_ratio("gemcitabine", "EGTA", "HCT-8Ea")
    alt = predicted_ratio("gemcitabine", "EGTA", "HCT-8Ea", logp_override=0.5)
    print(f"  S6  gemcitabine logP -1.4 -> +0.5: predicted ratio {base:.8f} -> "
          f"{alt:.8f} (rel change {abs(alt-base)/base:.3e}). "
          f"The descriptor gap is not load-bearing.")
    out["S6_gemcitabine_logp"] = dict(base_logp=-1.4, alt_logp=0.5,
                                      predicted_base=base, predicted_alt=alt,
                                      rel_change=abs(alt - base) / base)

    # S7: EGTA and bortezomib separately
    s7 = {}
    for pert in ("EGTA", "bortezomib"):
        sub = [r for r in prim if r["pert"] == pert]
        pv = [r["predicted"] for r in sub]
        mv = [r["measured_geomid"] for r in sub]
        g = gmfe(pv, mv)
        gn = gmfe([1.0] * len(sub), mv)
        ex = (math.log(spread(pv)) / math.log(spread(mv))) if spread(mv) > 1.0 else None
        s7[pert] = dict(n=len(sub), gmfe=g, gmfe_null=gn,
                        P1_pass=bool(g <= 2.0), P2_pass=bool(g <= 0.8 * gn),
                        measured_spread=spread(mv), predicted_spread=spread(pv),
                        log_fraction_explained=ex)
        exs = f"{ex:.4f}" if ex is not None else "n/a"
        print(f"  S7  {pert:<11} n={len(sub)}: GMFE {g:.4f} (null {gn:.4f}), "
              f"P1 {'PASS' if g<=2.0 else 'FAIL'}, "
              f"P2 {'PASS' if g<=0.8*gn else 'FAIL'}, ln-spread {exs}")
    out["S7_by_perturbation"] = s7

    # propagated packing-density uncertainty on the predicted ratio
    _hr("=")
    print("STAGE 9 -- PROPAGATED PACKING-DENSITY UNCERTAINTY ON THE PREDICTION")
    _hr("=")
    print("  The paper gives SD on every packing density. Propagate the extremes of")
    print("  NA_ctrl +/- SD and NA_trt +/- SD through the mapping (5-FU; all three")
    print("  compounds give the same ratio at the impermeable limit).")
    print(f"  {'pert':<13}{'line':<11}{'r_D nominal':>13}{'r_D low':>10}{'r_D high':>11}")
    unc = []
    for (pert, line), (na_c, na_t, _) in PACKING.items():
        sd_c, sd_t = PACKING_SD[(pert, line)]
        nom = predicted_ratio("5-fluorouracil", pert, line)
        vals = []
        for cc in (na_c - sd_c, na_c + sd_c):
            for tt in (na_t - sd_t, na_t + sd_t):
                pt = phi_treated(PHI_SHIPPED, cc, tt)
                if not (0.0 <= pt < 1.0):
                    continue
                hc = hindrance(PHI_SHIPPED, "5-fluorouracil")["hindrance"]
                ht = hindrance(pt, "5-fluorouracil")["hindrance"]
                vals.append(ht / hc)
        unc.append(dict(pert=pert, line=line, nominal=nom,
                        low=min(vals), high=max(vals)))
        print(f"  {pert:<13}{line:<11}{nom:>13.4f}{min(vals):>10.4f}{max(vals):>11.4f}")
    out["packing_sd_propagation"] = unc

    alpha_end = calibrate_alpha()
    assert alpha_end == ALPHA_FROZEN, f"calibrate_alpha END {alpha_end!r}"
    assert two_route.DEFAULT_CELL_VOLUME_FRACTION == PHI_SHIPPED, "phi moved"
    out["calibrate_alpha_end"] = alpha_end
    _hr("=")
    print(f"calibrate_alpha END   = {alpha_end!r}  (asserted == {ALPHA_FROZEN!r})")
    print(f"two_route.DEFAULT_CELL_VOLUME_FRACTION = "
          f"{two_route.DEFAULT_CELL_VOLUME_FRACTION} (unchanged)")
    _hr("=")

    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=False)
    print(f"wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
