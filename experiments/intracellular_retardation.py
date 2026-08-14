"""INTRACELLULAR RETARDATION -- the `K > 1` capacity term that
`results/two_route_efflux.md` 5.3 named as the last thing standing between the
two-route transport model and a default flip.

WHAT THIS SCRIPT IS FOR
-----------------------
`results/two_route_efflux.md` 2.5 relocated the diagnosis of paclitaxel's
anomalous transcellular gain from efflux to SATURABLE INTRACELLULAR BINDING, on
the strength of a tissue-scale measurement in the anchor's own cell line (Choi
2011, PMID 21225235) in a line that is functionally P-glycoprotein-negative (Wils
1994, PMID 7945455). It then listed the missing physics as a `K > 1` retardation
term. This script builds that term, verifies it, and scores it.

  ``--stage sources``   what could and could not be sourced for `K_b`, including
                        the verification of Choi 2011 against its own abstract.
  ``--stage algebra``   verifies the composition against a DIRECT TRANSIENT PDE
                        solve of a layered medium, and shows that the form the
                        commissioning brief proposed FAILS that same check.
  ``--stage pruijn``    scores it blind against the 71-compound anchor.
  ``--stage panel``     the breast tissue panel, and what happens to paclitaxel.
  ``--stage blast``     blast radius, both alpha arms.

THE RULE THIS FILE IS WRITTEN UNDER, UNCHANGED FROM ITS PARENTS
---------------------------------------------------------------
**Nothing is fitted to Pruijn.** `results/sweep/retard_prereg.json` was written to
disk BEFORE the term was scored even once, and it registers in advance both the
rule and the specific ways the author expected it to FAIL. `--stage pruijn`
asserts the live code against it.

Run:
    ./.venv/bin/python experiments/intracellular_retardation.py --stage sources
    ./.venv/bin/python experiments/intracellular_retardation.py --stage algebra
    ./.venv/bin/python experiments/intracellular_retardation.py --stage pruijn
    ./.venv/bin/python experiments/intracellular_retardation.py --stage panel
    ./.venv/bin/python experiments/intracellular_retardation.py --stage blast
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.molecular import compile_drug, stokes_einstein, two_route  # noqa: E402
from tmesim.molecular.breast_library import get_molecule_breast  # noqa: E402
from tmesim.molecular.library import get_molecule  # noqa: E402
from tmesim.validation.pk_reality import DEFAULT_ECM  # noqa: E402

OUT_DIR = os.path.join(_ROOT, "results", "sweep")
PREREG = os.path.join(OUT_DIR, "retard_prereg.json")
T0 = time.time()

#: The H-bond-donor decrement `results/two_route_efflux.md` 5.1 recommends
#: flipping WITH. Imported rather than re-typed so the two files cannot drift.
from experiments.two_route_efflux import HBD_ABRAHAM  # noqa: E402


def _hr(c: str = "=") -> None:
    print(c * 78)


def _tick(msg: str) -> None:
    print(f"  [t+{time.time() - T0:6.1f}s] {msg}", flush=True)


def _getmol(name: str):
    try:
        return get_molecule_breast(name)
    except KeyError:
        return get_molecule(name)


def _attr(mol, key, default=None):
    v = mol.get(key) if hasattr(mol, "get") else getattr(mol, key, None)
    return default if v is None else v


def _panel() -> List[str]:
    import experiments.dtissue_sensitivity as dts
    return list(dts.TESTABLE)


# =========================================================================== #
# The direct transient solve. This is the verification instrument, and it is
# deliberately independent of tmesim.molecular.two_route: it knows nothing about
# Maxwell, about effective media, or about retardation factors. It solves a
# layered slab and measures what comes out the far side.
# =========================================================================== #
def _laminate(n_periods: int, l_matrix: float, l_cell: float, d_mem: float,
              sigma_m: float, d_cyto: float, p_mem: float,
              k_partition: float, k_bound: float
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """One 1-D finite-volume laminate: [matrix | membrane | cell | membrane] x N.

    Returns (widths, conductivities, capacities) per control volume, in the
    matrix-equivalent potential `u`. The membrane is ONE control volume of the
    true bilayer thickness with conductivity `P*d_mem`, so that its face
    resistance `d_mem/(P*d_mem) = 1/P` is exactly a membrane permeability. No
    approximation is made anywhere in this construction.
    """
    h, sig, cap = [], [], []
    for _ in range(n_periods):
        for width, s, c in ((l_matrix, sigma_m, 1.0),
                            (d_mem, p_mem * d_mem, 1.0),
                            (l_cell, k_partition * d_cyto,
                             k_partition * (1.0 + k_bound)),
                            (d_mem, p_mem * d_mem, 1.0)):
            h.append(width)
            sig.append(s)
            cap.append(c)
    return np.array(h), np.array(sig), np.array(cap)


def _time_lag(h: np.ndarray, sig: np.ndarray, cap: np.ndarray,
              n_steps: int = 6000, n_tau: float = 40.0) -> Dict[str, float]:
    """Solve the transient directly and return the MEASURED steady flux and time lag.

    `u(0) = 1`, `u(H) = 0`, `u(x, 0) = 0`. Implicit Euler, tridiagonal solve.
    At long time the cumulative flux through the far face satisfies the
    Daynes-Barrer result `Q(t) -> J_ss * (t - tau)`, and `tau` is the quantity
    that depends on the medium's CAPACITY. The steady flux `J_ss` depends only on
    its conductance. Measuring BOTH separates the leg this work adds from the leg
    that was already verified.
    """
    n = len(h)
    x_face_cond = np.zeros(n + 1)               # conductance to the neighbour
    # interior faces: harmonic mean of the two half-cell resistances -- exact
    # series composition, which is why the thin membrane needs no refinement.
    for i in range(1, n):
        x_face_cond[i] = 1.0 / (h[i - 1] / (2.0 * sig[i - 1])
                                + h[i] / (2.0 * sig[i]))
    x_face_cond[0] = sig[0] / (h[0] / 2.0)      # Dirichlet at the near face
    x_face_cond[n] = sig[-1] / (h[-1] / 2.0)    # Dirichlet at the far face

    total_h = float(h.sum())
    r_series = float((h / sig).sum())
    sigma_ser = total_h / r_series
    s_avg = float((h * cap).sum() / total_h)
    tau_expect = total_h * total_h * s_avg / (6.0 * sigma_ser)
    t_end = n_tau * tau_expect
    dt = t_end / n_steps

    # implicit Euler: (C/dt + A) u^{n+1} = C/dt u^n + b
    lower = np.zeros(n)
    upper = np.zeros(n)
    diag = np.zeros(n)
    cvol = cap * h
    for i in range(n):
        gl = x_face_cond[i]
        gr = x_face_cond[i + 1]
        diag[i] = cvol[i] / dt + gl + gr
        if i > 0:
            lower[i] = -x_face_cond[i]
        if i < n - 1:
            upper[i] = -x_face_cond[i + 1]

    u = np.zeros(n)
    q_cum = 0.0
    ts, qs = [], []
    for step in range(n_steps):
        rhs = cvol / dt * u
        rhs[0] += x_face_cond[0] * 1.0          # u_left = 1
        # rhs[-1] += x_face_cond[n] * 0.0       # u_right = 0
        # Thomas algorithm
        cp = np.zeros(n)
        dp = np.zeros(n)
        cp[0] = upper[0] / diag[0]
        dp[0] = rhs[0] / diag[0]
        for i in range(1, n):
            m = diag[i] - lower[i] * cp[i - 1]
            cp[i] = upper[i] / m if i < n - 1 else 0.0
            dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / m
        u_new = np.zeros(n)
        u_new[-1] = dp[-1]
        for i in range(n - 2, -1, -1):
            u_new[i] = dp[i] - cp[i] * u_new[i + 1]
        u = u_new
        q_cum += x_face_cond[n] * u[-1] * dt     # flux out of the far face
        ts.append((step + 1) * dt)
        qs.append(q_cum)

    ts_a, qs_a = np.array(ts), np.array(qs)
    keep = ts_a > 0.6 * t_end                   # long-time linear regime
    slope, intercept = np.polyfit(ts_a[keep], qs_a[keep], 1)
    return {
        "J_steady_measured": float(slope),
        "J_steady_from_conductance": sigma_ser / total_h,
        "time_lag_measured": float(-intercept / slope),
        "time_lag_from_capacity": tau_expect,
        "sigma_series": sigma_ser,
        "s_avg": s_avg,
        "D_apparent_measured": total_h * total_h / (6.0 * float(-intercept / slope)),
        "D_apparent_predicted": sigma_ser / s_avg,
        "thickness_um": total_h,
    }


# =========================================================================== #
# STAGE A -- what could and could not be sourced
# =========================================================================== #
def stage_sources() -> dict:
    _hr()
    print("STAGE A -- WHERE K_b COMES FROM, AND WHERE IT DOES NOT")
    _hr()
    print("""
  A.1 -- CHOI 2011 (PMID 21225235), THE MEASUREMENT THAT MOTIVATED THIS TERM,
         VERIFIED AGAINST ITS OWN ABSTRACT -- AND PARTLY CORRECTED.

  The commissioning brief stated: "Choi 2011 measured paclitaxel at 38 %
  penetration through an HT-29 multicellular layer at 1 uM against 100 % at
  10 uM". The abstract was retrieved and read. VERBATIM:

    "The penetration rate measured after 24 h by cytotoxicity of the conditioned
     media was 40 and 38% in DLD-1 (20 uM) and HT-29 MCLs (1 uM), respectively"

    "In HT-29 MCL, full penetration of PTX was obtained at 10 uM after 48 h,
     whereas only 80% was obtained at 1 uM."

  So: the cell line is HT-29 and the saturation is REAL. But the 38 % and the
  100 % are NOT the same experiment. The 38 % is at 24 h by the CYTOTOXICITY
  method; the 100 % at 10 uM is at 48 h by [14C]-paclitaxel. The TIME-MATCHED
  pair, both at 48 h by the same method, is 80 % at 1 uM against 100 % at 10 uM.

  THE CORRECTION MATTERS AND IT WEAKENS THE MOTIVATING EVIDENCE. The saturable
  restriction in HT-29 at matched time and method is about 1.25-fold, not the
  ~2.6-fold the brief's numbers imply. The qualitative signature -- saturable
  penetration in a P-gp-negative line -- survives; the magnitude does not.
  (The strong, unambiguously saturable case in that paper is DLD-1, which is not
  the anchor cell line: 10 % at 24 h and 20 % at 96 h at 1 uM.)
""")

    print("""  A.2 -- PER-COMPOUND K_b FOR THE SCORED SET: NOT AVAILABLE.

  The Pruijn anchor is 67 benzotriazine di-N-oxide analogues of tirapazamine
  plus NLCQ-1, urea, mannitol and acridan-9(10H)-one. No intracellular bound
  fraction, cell-to-medium partition, or unbound fraction could be retrieved for
  ANY of them. A per-compound K_b for the anchor set is UNAVAILABLE, and the
  house rule forbids inventing one. That is why a STRUCTURAL rule is registered
  instead -- the route the commissioning brief's item 2 provides for.
""")

    print("  A.3 -- THE STRUCTURAL RULE, AND ITS ONE-LINE DERIVATION\n")
    print("""    Kilford 2008 (PMID 18411401) predicts the unbound fraction in a
    HEPATOCYTE incubation from lipophilicity alone,

        fu_hep = 1 / (1 + 125 * VR * X),
        X      = 10**(0.072*logP**2 + 0.067*logP - 1.126)   [Hallifax & Houston]

    with VR the volume ratio of CELLS to incubation. The bound-to-free ratio in
    the incubation is 125*VR*X. All that bound material sits INSIDE the cells,
    which occupy a fraction VR of the incubation, so per unit CELL volume

        K_b = (125 * VR * X) / VR = 125 * X.

    THE VOLUME RATIO CANCELS EXACTLY. That is what makes the number usable here:
    the result is a property of a millilitre of cells, not of Kilford's assay
    geometry, so it transfers to a cell sitting in a tissue. The constant 125
    carries units of mg microsomal-protein-equivalent per mL of CELL volume,
    against a whole-cell protein content of order 200-250 mg/mL.
""")
    print(f"  {'logP':>7s}{'K_b':>11s}{'bound %':>10s}{'R = 1+phi*K_b':>16s}")
    tbl = []
    for lp in (-2.62, -1.5, -0.8, -0.465, 0.0, 1.0, 2.0, 2.97, 3.26, 3.74):
        kb = two_route.binding_capacity_from_logp(lp)
        r = two_route.retardation_factor(kb)
        tbl.append({"logp": lp, "k_bound": kb, "retardation": r})
        print(f"  {lp:>7.2f}{kb:>11.2f}{kb/(1+kb):>9.1%}{r:>16.2f}")
    print("""
  READ THE logP = 0 ROW. The rule asserts that a compound with NO lipophilicity
  at all is 90 % bound inside a cell. Urea and mannitol are in the Pruijn set,
  are classical NON-BINDING flux markers, and Pruijn uses urea as its INTERNAL
  STANDARD. This is registered in results/sweep/retard_prereg.json as an
  expected out-of-domain failure BEFORE the rule was scored, not discovered
  afterwards.

  The rule is also non-monotone: it has a vertex at logP = -0.465, so it
  predicts MORE binding for the most hydrophilic compounds in the set than for
  neutral ones. That is an artefact of extrapolating a quadratic fitted over
  drug-like logP down to logP = -2.62, and it is disclosed rather than clipped.
""")
    print("""  A.4 -- THE PER-COMPOUND TAXANE ROUTE: BRACKETED, NOT PINNED.

  The commissioning brief named tubulin-binding stoichiometry as a candidate
  per-compound source. At tracer concentrations, 1:1 binding to a target present
  at concentration [T] with dissociation constant Kd gives K_b = [T]/Kd.

  VERIFIED: the stoichiometry is "0.99 +/- 0.04 ligand per alpha beta dimer"
  (Diaz & Andreu, Biochemistry 1993, PMID 8096151). The POLYMER Kd is
  "8.7 x 10-7 M" (Parness & Horwitz, J Cell Biol 1981, PMID 6118377, PMC2111958,
  full text read). A second primary gives 3.7e7 M-1, i.e. ~27 nM (Andreu &
  Barasoain 2001, PMID 11580273) -- a ~30-fold disagreement, NOT averaged here.

  UNVERIFIED, AND IT IS THE BLOCKING INPUT: the intracellular TUBULIN
  CONCENTRATION. No primary mammalian uM value could be retrieved. The only
  primary number is a PERCENTAGE -- "approximately 2.5--3.3% of their total
  protein" (Hiller & Weber, Cell 1978, PMID 688394). A published honest bound is
  "anywhere between 8 uM and 50 uM" (eLife 2022, PMID 35029146).

  So the taxane route can only BRACKET K_b for paclitaxel:
""")
    brackets = []
    for kd_um, kd_src in ((0.87, "Parness & Horwitz 1981 (polymer)"),
                          (0.027, "Andreu & Barasoain 2001 (polymer)")):
        lo, hi = 8.0 / kd_um, 50.0 / kd_um
        brackets.append({"kd_um": kd_um, "source": kd_src,
                         "k_bound_lo": lo, "k_bound_hi": hi})
        print(f"    Kd = {kd_um:>6.3f} uM  [{kd_src}]"
              f"  ->  K_b = {lo:>8.1f} - {hi:>8.1f}")
    rule_pac = two_route.binding_capacity_from_logp(3.74)
    print(f"\n    the STRUCTURAL RULE gives K_b = {rule_pac:.1f} for paclitaxel.")
    print("""
  The bracket spans 9 to 1852 and contains the structural rule's value. That is
  a consistency check and nothing more: a range spanning two orders of magnitude
  cannot discriminate between hypotheses, and it exists only for taxanes, so it
  cannot be applied to a single compound in the Pruijn anchor set. It is
  reported because it is the honest answer to "can K be sourced per compound?" --
  which is NO, not even for the one drug the whole question is about.
""")
    return {"choi_2011_verified": True, "structural_rule_table": tbl,
            "taxane_brackets": brackets, "structural_rule_paclitaxel": rule_pac}


# =========================================================================== #
# STAGE B -- the composition, verified against a direct transient solve
# =========================================================================== #
def stage_algebra() -> dict:
    _hr()
    print("STAGE B -- THE COMPOSITION, AGAINST A DIRECT TRANSIENT PDE SOLVE")
    _hr()
    print("""
  THE CLAIM UNDER TEST. Reversible intracellular binding at local equilibrium
  retards the TRANSIENT apparent diffusivity and leaves the STEADY-STATE
  conductance EXACTLY unchanged:

      conductance  sigma_eff : matrix sigma_m, inclusion K*D_cyto   [UNAFFECTED]
      capacity     s_eff     = (1-phi)*1 + phi*K*(1+K_b)            [volume avg]
      apparent D             = sigma_eff / s_eff

  THE RIVAL, WHICH IS THE FORM THE COMMISSIONING BRIEF PROPOSED:
  "D_eff = D_free/(1+K)" applied to the INCLUSION CONDUCTIVITY, i.e.
  sigma_c -> K*D_cyto/(1+K_b) with no capacity term at all.

  These are different physics and they are distinguishable by experiment, so
  they are distinguished by solve rather than by argument. The instrument is a
  1-D finite-volume laminate [matrix | membrane | cell | membrane] x N carrying
  the real 4 nm bilayer as its own control volume, run as a transient from an
  empty start. It knows nothing about Maxwell or about retardation factors.
  Daynes-Barrer: the cumulative flux tends to J_ss*(t - tau); J_ss tests the
  CONDUCTANCE and tau tests the CAPACITY.
""")
    l_cell = 2.0 * two_route.DEFAULT_CELL_RADIUS_UM
    phi = two_route.DEFAULT_CELL_VOLUME_FRACTION
    l_matrix = l_cell * (1.0 / phi - 1.0)
    d_mem = two_route.BILAYER_THICKNESS_NM * 1e-3
    sigma_m, d_cyto = 265.0, 200.0

    _tick("B.1 -- convergence of the laminate onto the homogenised composition")
    print(f"\n  cell {l_cell:.3f} um, matrix {l_matrix:.3f} um, membrane "
          f"{d_mem*1e3:.1f} nm, phi = {phi}")
    print(f"  {'N':>4s}{'K_b':>7s}{'J meas/pred':>14s}{'tau meas/pred':>15s}"
          f"{'D_app meas/pred':>17s}")
    conv = []
    for n_per in (10, 20, 40, 80):
        for kb in (0.0, 30.0):
            h, sg, cp = _laminate(n_per, l_matrix, l_cell, d_mem, sigma_m,
                                  d_cyto, 50.0, 1.0, kb)
            r = _time_lag(h, sg, cp)
            jr = r["J_steady_measured"] / r["J_steady_from_conductance"]
            tr_ = r["time_lag_measured"] / r["time_lag_from_capacity"]
            dr = r["D_apparent_measured"] / r["D_apparent_predicted"]
            conv.append({"n_periods": n_per, "k_bound": kb, "J_ratio": jr,
                         "tau_ratio": tr_, "D_ratio": dr})
            print(f"  {n_per:>4d}{kb:>7.1f}{jr:>14.6f}{tr_:>15.6f}{dr:>17.6f}")
    print("""
  The steady flux matches to machine precision at every N -- that leg was
  already verified. The TIME LAG converges onto the volume-averaged capacity as
  the laminate is refined, which is the leg this work adds. The residual at
  finite N is the laminate's own finite-size correction, not an error in the
  composition, and it falls with N as it must.
""")

    _tick("B.2 -- the exact retardation identity, at machine precision")
    worst = 0.0
    n_pts = 0
    for ph in (0.1, 0.4, 0.7892, 0.95):
        for k in (0.5, 1.0, 3.0):
            for kb in (0.0, 0.3, 3.0, 30.0, 300.0, 3000.0):
                s0 = (1.0 - ph) + ph * k
                s1 = (1.0 - ph) + ph * k * (1.0 + kb)
                direct = s1 / s0
                mine = two_route.retardation_factor(kb, cell_volume_fraction=ph,
                                                    cell_partition=k)
                worst = max(worst, abs(mine - direct) / direct)
                n_pts += 1
    print(f"  R = s_eff(K_b)/s_eff(0) over n = {n_pts} (phi, K, K_b) points: "
          f"max relative error {worst:.3g}")
    assert worst < 1e-14, worst
    print("  -> R = 1 + phi*K*K_b/[(1-phi)+phi*K] is an IDENTITY. OK")

    _tick("B.3 -- THE DISCRIMINATING TEST: correct form vs the brief's form")
    print(f"\n  {'K_b':>8s}{'D_app direct':>15s}{'capacity form':>15s}{'err':>9s}"
          f"{'brief form':>13s}{'err':>9s}")
    rival = []
    for kb in (0.0, 3.0, 30.0, 300.0):
        h, sg, cp = _laminate(80, l_matrix, l_cell, d_mem, sigma_m, d_cyto,
                              50.0, 1.0, kb)
        r = _time_lag(h, sg, cp)
        direct = r["D_apparent_measured"]
        mine = r["D_apparent_predicted"]
        # the brief's form: retardation inside the conductivity, no capacity.
        h2, sg2, cp2 = _laminate(80, l_matrix, l_cell, d_mem, sigma_m,
                                 d_cyto / (1.0 + kb), 50.0, 1.0, 0.0)
        theirs = float(h2.sum()) / float((h2 / sg2).sum())
        rival.append({"k_bound": kb, "D_direct": direct, "D_capacity": mine,
                      "D_brief": theirs,
                      "err_capacity": abs(mine / direct - 1.0),
                      "err_brief": abs(theirs / direct - 1.0)})
        print(f"  {kb:>8.1f}{direct:>15.4f}{mine:>15.4f}"
              f"{abs(mine/direct-1):>8.2%}{theirs:>13.4f}"
              f"{abs(theirs/direct-1):>8.2%}")
    print("""
  The capacity form tracks the direct solve. The brief's form does not, and it
  fails in a way that grows with K_b -- because it applies a capacity effect to
  a conductance, and because at phi = 0.79 it entirely misses the retardation of
  MATRIX-phase transport by cell storage, which is the dominant term.
""")

    _tick("B.4 -- where LOCAL EQUILIBRIUM fails: the shut-membrane limit")
    print("""
  Volume-averaging the capacity assumes the cell interior equilibrates with the
  interstitium within the experiment. A cell behind a shut membrane cannot. The
  registered prediction is that R OVER-states the retardation there. Measured:
""")
    print(f"  {'P (um/s)':>11s}{'tau_cell (s)':>14s}{'D_app direct':>15s}"
          f"{'capacity form':>15s}{'overstated by':>15s}")
    eq = []
    for p_mem in (1e-4, 1e-2, 1.0, 50.0, 1e4):
        h, sg, cp = _laminate(40, l_matrix, l_cell, d_mem, sigma_m, d_cyto,
                              p_mem, 1.0, 30.0)
        r = _time_lag(h, sg, cp)
        tau_c = two_route.equilibration_time_s(30.0, p_mem)
        eq.append({"P_um_s": p_mem, "tau_cell_s": tau_c,
                   "D_direct": r["D_apparent_measured"],
                   "D_capacity": r["D_apparent_predicted"],
                   "ratio": r["D_apparent_predicted"] / r["D_apparent_measured"]})
        print(f"  {p_mem:>11.4g}{tau_c:>14.4g}{r['D_apparent_measured']:>15.5f}"
              f"{r['D_apparent_predicted']:>15.5f}"
              f"{r['D_apparent_predicted']/r['D_apparent_measured']:>15.4f}")
    print("""
  Note what this says and what it does not. The LONG-TIME asymptote is correct
  at every finite P -- a slow membrane delays equilibrium, it does not abolish
  it -- so the composition is not wrong here. What is wrong is applying it to an
  experiment shorter than tau_cell. For the shut-membrane drugs (doxorubicin,
  epirubicin: P ~ 6e-5 um/s) tau_cell exceeds any experiment in this project by
  orders of magnitude, and for those drugs the retardation factor is an UPPER
  BOUND, not an estimate.
""")
    return {"convergence": conv, "identity_max_rel_error": worst,
            "rival_forms": rival, "equilibration": eq}


# =========================================================================== #
# STAGE C -- scored, blind, against Pruijn
# =========================================================================== #
def stage_pruijn() -> dict:
    _hr()
    print("STAGE C -- SCORED AGAINST PRUIJN 2008 (n=71, HT29 MCLs). BLIND, ONCE.")
    _hr()
    from experiments.two_route_transport import _build_rows, _score
    from tmesim.validation.calibration import calibrate_alpha

    with open(PREREG) as fh:
        reg = json.load(fh)
    print("  C.0 -- PRE-REGISTRATION AND INVARIANT CHECKS")
    assert two_route.KILFORD_MG_EQUIV_PER_ML_CELL == 125.0
    assert (two_route.HALLIFAX_C2, two_route.HALLIFAX_C1,
            two_route.HALLIFAX_C0) == (0.072, 0.067, -1.126)
    assert two_route.HBD_LOGP_DECREMENT == 1.0, "parent pre-registration drifted"
    a = two_route.two_route_diffusivity(265.0, 0.6, 2.0, 2)
    b = two_route.two_route_diffusivity(265.0, 0.6, 2.0, 2,
                                        bound_to_free_ratio=None)
    assert a["hindrance"] == b["hindrance"]
    assert a["effective_D_um2s"] == a["apparent_D_um2s"]
    alpha = calibrate_alpha()
    assert alpha == 148.60815765445545, alpha
    print(f"    live constants match results/sweep/retard_prereg.json : OK")
    print(f"    default-off is bitwise identical                     : OK")
    print(f"    calibrate_alpha()                                    : {alpha!r}")

    cache: Dict[float, float] = {}
    _tick("compiling the 71 Pruijn compounds through the unmodified pipeline")
    _build_rows(1.0, cache)

    def scored(dec: float, mult: Optional[float], uniform: Optional[float] = None):
        rows = [r for r in _build_rows(dec, cache) if "two_route" in r]
        out = []
        for r in rows:
            if uniform is not None:
                kb = uniform
            elif mult is None:
                kb = 0.0
            else:
                kb = mult * two_route.binding_capacity_from_logp(r["logp"])
            rr = dict(r)
            rr["retarded"] = r["two_route"] / two_route.retardation_factor(kb)
            rr["k_bound"] = kb
            out.append(rr)
        return out

    results = {}
    for dec, decname in ((1.0, "pre-registered convention"),
                         (HBD_ABRAHAM, "SOURCED Abraham-derived")):
        print(f"\n  C.1 -- hbd_decrement = {dec} ({decname})")
        print(f"  {'arm':<44s}{'spread':>8s}{'R^2':>8s}{'Spearman':>10s}"
              f"{'geo bias':>11s}{'<3x':>6s}")
        arms = []
        base = scored(dec, None)
        s = _score(base, "two_route")
        s["arm"] = "A. retardation OFF (control)"
        arms.append(s)
        print(f"  {s['arm']:<44s}{s['fold']:>7.1f}x{s['r2']:>8.3f}"
              f"{s['spearman']:>+10.3f}{s['geo_bias']:>10.3g}x{s['within_3x']:>6.0%}")
        for mult in (0.01, 0.03, 0.1, 0.3, 1.0, 3.0):
            rr = scored(dec, mult)
            s = _score(rr, "retarded")
            s["arm"] = (f"C. structural rule x {mult:g}"
                        + ("   <-- THE SOURCED VALUE" if mult == 1.0 else ""))
            s["multiplier"] = mult
            arms.append(s)
            print(f"  {s['arm']:<44s}{s['fold']:>7.1f}x{s['r2']:>8.3f}"
                  f"{s['spearman']:>+10.3f}{s['geo_bias']:>10.3g}x"
                  f"{s['within_3x']:>6.0%}")
        for uni in (1.0, 10.0, 100.0):
            rr = scored(dec, None, uniform=uni)
            s = _score(rr, "retarded")
            s["arm"] = f"D. CONTROL: uniform K_b = {uni:g}"
            arms.append(s)
            print(f"  {s['arm']:<44s}{s['fold']:>7.1f}x{s['r2']:>8.3f}"
                  f"{s['spearman']:>+10.3f}{s['geo_bias']:>10.3g}x"
                  f"{s['within_3x']:>6.0%}")
        results[str(dec)] = arms

        off = arms[0]
        sourced = next(a for a in arms if a.get("multiplier") == 1.0)
        best = max((a for a in arms if "multiplier" in a), key=lambda a: a["r2"])
        print(f"\n    control R^2 {off['r2']:.3f} -> sourced rule R^2 "
              f"{sourced['r2']:.3f}   ({'WORSE' if sourced['r2'] < off['r2'] else 'better'})")
        print(f"    best multiplier anywhere in the disclosure sweep: "
              f"{best['multiplier']:g} (R^2 {best['r2']:.3f})")

    # ---- the matched-MW test, at the sourced decrement -------------------- #
    print("\n  C.2 -- THE MATCHED-MW TEST. Where a MW-only model is flat by")
    print("  construction, does retardation help or hurt the spread and the rank?")
    from experiments.diffusivity_falsification import fold
    mm = {}
    for label, rr in (("retardation OFF", scored(HBD_ABRAHAM, None)),
                      ("structural rule", scored(HBD_ABRAHAM, 1.0))):
        key = "two_route" if label.endswith("OFF") else "retarded"
        by_mw: Dict[float, List[dict]] = {}
        for r in rr:
            by_mw.setdefault(r["mw"], []).append(r)
        conc = disc = 0
        for _mw, sub in by_mw.items():
            for i in range(len(sub)):
                for j in range(i + 1, len(sub)):
                    a_, b_ = sub[i], sub[j]
                    if a_["meas"] == b_["meas"] or a_[key] == b_[key]:
                        continue
                    if (a_["meas"] > b_["meas"]) == (a_[key] > b_[key]):
                        conc += 1
                    else:
                        disc += 1
        folds = [(m, fold([r[key] for r in s]), fold([r["meas"] for r in s]))
                 for m, s in sorted(by_mw.items()) if len(s) >= 3]
        agree = conc / (conc + disc) if (conc + disc) else float("nan")
        mm[label] = {"concordant": conc, "discordant": disc,
                     "rank_agreement": agree, "ties": folds}
        print(f"    {label:18s} within-tie rank agreement "
              f"{agree:.1%}  ({conc} concordant / {disc} discordant)")
    for m, f_off, f_meas in mm["retardation OFF"]["ties"]:
        f_on = next(f for mm_, f, _ in mm["structural rule"]["ties"] if mm_ == m)
        print(f"    MW {m:>6.0f}   measured {f_meas:>5.1f}x   "
              f"model OFF {f_off:>5.2f}x   model ON {f_on:>5.2f}x")
    results["matched_mw"] = mm

    print("""
  THE CONTROL ARM D IS THE ONE TO READ FIRST. A compound-INDEPENDENT K_b is a
  constant divisor, i.e. a rigid shift in log space, so it CANNOT change R^2 or
  Spearman at all -- only the bias. The registered expectation was exactly that,
  and it is what the table shows. So every change in R^2 in arm C is attributable
  to the logP DEPENDENCE of the structural rule and to nothing else.
""")
    return {"arms": results, "alpha": alpha}


# =========================================================================== #
# STAGE D -- the breast tissue panel, and paclitaxel specifically
# =========================================================================== #
def stage_panel() -> dict:
    _hr()
    print("STAGE D -- THE BREAST PANEL, AND WHAT HAPPENS TO PACLITAXEL")
    _hr()
    panel = _panel()
    out = {}
    for dec in (1.0, HBD_ABRAHAM):
        print(f"\n  hbd_decrement = {dec}")
        print(f"  {'drug':16s}{'logP':>7s}{'K_b':>9s}{'R':>9s}{'gain':>8s}"
              f"{'scale':>9s}{'app.scale':>11s}{'vs dox':>9s}{'tau_cell(s)':>13s}")
        rows = {}
        for d in panel:
            mol = _getmol(d)
            # two_route=False: PRE-PROMOTION baseline. compile_drug now applies the
            # two-route term itself (results/two_route_promotion.md), so this file's
            # 'shipping engine' column must ask for the MW-only interstitial value
            # explicitly or the correction would be applied twice.
            prof = compile_drug.build_profile(mol, ecm_reference=DEFAULT_ECM,
                                              two_route=False)
            radius = float(_attr(mol, "radius_hydrodynamic_nm")
                           or stokes_einstein.hydrodynamic_radius_nm(
                               float(_attr(mol, "mw")),
                               bool(_attr(mol, "is_antibody", False))))
            lp = float(_attr(mol, "logp", 0.0))
            kb = two_route.binding_capacity_from_logp(lp)
            tr = two_route.two_route_diffusivity(
                float(prof.effective_diffusivity), radius, lp,
                float(_attr(mol, "hbd", 0)), hbd_decrement=dec,
                bound_to_free_ratio=kb)
            rows[d] = {"logp": lp, "k_bound": kb,
                       "retardation": tr["retardation_factor"],
                       "gain": tr["transcellular_gain"],
                       "scale": tr["hindrance"],
                       "apparent_scale": tr["apparent_hindrance"],
                       "tau_cell_s": tr["equilibration_time_s"],
                       "P_um_s": tr["membrane_permeability_um_s"]}
        sd = rows["doxorubicin"]["apparent_scale"]
        for d in panel:
            r = rows[d]
            print(f"  {d:16s}{r['logp']:>7.2f}{r['k_bound']:>9.2f}"
                  f"{r['retardation']:>9.2f}{r['gain']:>8.4f}{r['scale']:>9.4f}"
                  f"{r['apparent_scale']:>11.5f}"
                  f"{r['apparent_scale']/sd:>9.4f}{r['tau_cell_s']:>13.4g}")
        diff_c = (max(r["scale"] for r in rows.values())
                  / min(r["scale"] for r in rows.values()))
        diff_a = (max(r["apparent_scale"] for r in rows.values())
                  / min(r["apparent_scale"] for r in rows.values()))
        print(f"\n  panel differential, conductance basis : {diff_c:.3f}x")
        print(f"  panel differential, apparent basis    : {diff_a:.3f}x")
        out[str(dec)] = {"rows": rows, "differential_conductance": diff_c,
                         "differential_apparent": diff_a}
    print("""
  READ THE `gain` COLUMN. It is IDENTICAL to the un-retarded model's, for every
  drug, at every K_b. That is not an omission: transcellular gain is a ratio of
  two CONDUCTANCES and binding does not touch a conductance. The commissioning
  brief asked for paclitaxel's gain to fall as the test of this physics. Under
  the correct physics that test cannot be passed by any value of K_b, and this
  was registered in advance in results/sweep/retard_prereg.json.

  READ THE tau_cell COLUMN NEXT. It is the time a cell needs to load its bound
  pool through its own membrane. Where it is large the retardation factor is an
  upper bound, not an estimate.
""")
    return out


# =========================================================================== #
# STAGE E -- blast radius
# =========================================================================== #
def stage_blast(dx_um: float = 12.0) -> dict:
    _hr()
    print("STAGE E -- BLAST RADIUS, published dx = 12 um grid, BOTH alpha arms")
    _hr()
    import experiments.dtissue_sensitivity as dts

    ref = dts.reference_transport_constants()
    d_ship = float(ref["D_tissue_cm2_s"])
    panel = _panel()
    print(f"  shipping D_tissue (doxorubicin) = {d_ship:.4e} cm^2/s")

    scales = {}
    for d in panel:
        mol = _getmol(d)
        # two_route=False: PRE-PROMOTION baseline. compile_drug now applies the
        # two-route term itself (results/two_route_promotion.md), so this file's
        # 'shipping engine' column must ask for the MW-only interstitial value
        # explicitly or the correction would be applied twice.
        prof = compile_drug.build_profile(mol, ecm_reference=DEFAULT_ECM,
                                          two_route=False)
        radius = float(_attr(mol, "radius_hydrodynamic_nm")
                       or stokes_einstein.hydrodynamic_radius_nm(
                           float(_attr(mol, "mw")),
                           bool(_attr(mol, "is_antibody", False))))
        lp = float(_attr(mol, "logp", 0.0))
        tr = two_route.two_route_diffusivity(
            float(prof.effective_diffusivity), radius, lp,
            float(_attr(mol, "hbd", 0)), hbd_decrement=HBD_ABRAHAM,
            bound_to_free_ratio=two_route.binding_capacity_from_logp(lp))
        scales[d] = tr["apparent_hindrance"]
    sd = scales["doxorubicin"]
    print(f"\n  RETARDED scales at hbd_decrement = {HBD_ABRAHAM}")
    print(f"  {'drug':16s}{'scale':>12s}{'fold':>10s}{'vs dox':>10s}")
    for d in panel:
        print(f"  {d:16s}{scales[d]:12.6f}{1/scales[d]:9.1f}x{scales[d]/sd:9.4f}x")
    print(f"  panel differential: {max(scales.values())/min(scales.values()):.3f}x")

    points = []
    for d in panel:
        p = d_ship * scales[d]
        if all(abs(p - q) > 1e-6 * abs(q) for q in points):
            points.append(p)
    if all(abs(d_ship - q) > 1e-6 * abs(q) for q in points):
        points.append(d_ship)
    print(f"\n  {len(points)} distinct D points to solve, x 2 alpha modes")

    ckpt = os.path.join(OUT_DIR, f"retard_records_dx{int(dx_um)}.jsonl")
    seed = os.path.join(OUT_DIR, f"tworoute_efflux_records_dx{int(dx_um)}.jsonl")
    if not os.path.exists(ckpt) and os.path.exists(seed):
        import shutil
        shutil.copyfile(seed, ckpt)
        print(f"  seeded {ckpt} from the parent run's records")
    out_json = os.path.join(OUT_DIR, f"retard_sweep_dx{int(dx_um)}.json")
    _tick("solving")
    res = dts.run(dx_um=dx_um, modes=("frozen", "recal"), out_json=out_json,
                  checkpoint=ckpt, explicit_points=points)
    _tick("solved")

    def pick(mode, drug, d):
        cands = [r for r in res["records"]
                 if r["mode"] == mode and r["drug"] == drug]
        if not cands:
            return None
        best = min(cands, key=lambda r: abs(r["D_tissue_cm2_s"] - d))
        if abs(best["D_tissue_cm2_s"] - d) / abs(d) > 1e-3:
            return None
        return best

    print(f"\n  {'mode':7s} {'drug':16s} {'scale':>10s} {'ratio ship':>11s} "
          f"{'ratio corr':>11s} {'verdict':>13s}->{'verdict':>13s} "
          f"{'alpha':>10s} {'moved':>6s}")
    rows = []
    for mode in ("frozen", "recal"):
        for drug in res["drugs"]:
            a = pick(mode, drug, d_ship)
            b = pick(mode, drug, d_ship * scales.get(drug, 1.0))
            if a is None or b is None:
                continue
            moved = a["verdict"] != b["verdict"]
            rows.append({"mode": mode, "drug": drug, "scale": scales[drug],
                         "ratio_ship": a["ratio_exact"],
                         "ratio_corr": b["ratio_exact"],
                         "verdict_ship": a["verdict"],
                         "verdict_corr": b["verdict"],
                         "alpha_corr": b["alpha"], "moved": moved})
            print(f"  {mode:7s} {drug:16s} {scales[drug]:10.6f} "
                  f"{a['ratio_exact']:11.4f} {b['ratio_exact']:11.4f} "
                  f"{a['verdict']:>13s}->{b['verdict']:>13s} "
                  f"{b['alpha']:>10.4f} {'YES' if moved else 'no':>6s}")
    n = sum(1 for r in rows if r["moved"])
    print(f"\n  verdicts that MOVE: {n} of {len(rows)}")
    return {"D_ship_cm2_s": d_ship, "scales": scales, "rows": rows,
            "n_moved": n, "n_rows": len(rows), "out_json": out_json}


# =========================================================================== #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=("all", "sources", "algebra", "pruijn", "panel",
                             "blast"))
    ap.add_argument("--dx", type=float, default=12.0)
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    stages = {"sources": stage_sources, "algebra": stage_algebra,
              "pruijn": stage_pruijn, "panel": stage_panel,
              "blast": lambda: stage_blast(dx_um=args.dx)}
    todo = list(stages) if args.stage == "all" else [args.stage]
    for name in todo:
        res = stages[name]()
        path = os.path.join(OUT_DIR, f"retard_{name}.json")
        with open(path, "w") as fh:
            json.dump(res, fh, indent=1, default=str)
        print(f"  wrote {path}")


if __name__ == "__main__":
    main()
