"""logP sensitivity sweep: how far does a wrong logP move a simulated output?

WHY THIS IS CHEAP (the structural result that makes a brute-force sweep unnecessary)
------------------------------------------------------------------------------------
logP never enters the transport coefficient. Tracing
``molecular.compile_drug.build_profile``:

    DrugParams.diffusivity  <- profile.effective_diffusivity
                            <- brownian sim x ecm_mesh, driven ONLY by
                               radius_hydrodynamic_nm <- f(MW)          [NO logP]
    DrugParams.binding      <- binding.tissue_binding_rate(logp, hbd)   [logP]
    DrugParams.clearance    <- binding.clearance_rate(logp, mw)         [logP]
    DrugParams.bbb_penetrance <- permeability.bbb_penetrance(mw, logp,
                                 tpsa, hbd)                             [logP]

So a logP perturbation cannot move D_eff. It moves only the first-order loss
rate k = clearance + binding and the delivery-source scalar bbb_penetrance.
That means we can compile each drug ONCE (paying the Brownian sim exactly once)
and then re-evaluate every logP scenario in closed form, because the only
quantities that change are pure functions of the descriptor dict.

The headline metric is the analytic reaction-diffusion penetration length

    L = sqrt(D_eff / k_eff),   k_eff = k_passive + alpha * uptake_proxy

which is the quantity ``validation.calibration`` is calibrated against
(doxorubicin, Primeau et al. 2005, PMID 16361566, 40-50 um half-distance
converted to a ~64.9 um e-folding length).

THE ALPHA TRAP (why two columns)
--------------------------------
``calibration.calibrate_alpha`` fits the single constant alpha on ONE anchor
drug: DOXORUBICIN. alpha = (k_target - k_passive_dox) / uptake_proxy_dox.

  * Perturbing a NON-anchor drug's logP leaves alpha untouched, so the
    "alpha fixed" and "alpha refit" columns are IDENTICAL by construction.
  * Perturbing DOXORUBICIN's logP changes uptake_proxy_dox, so alpha moves
    inversely. Doxorubicin's own L snaps back to the fit target (that is the
    fit, not a prediction), while EVERY OTHER DRUG's k_eff is rescaled by the
    new alpha. A naive per-drug sweep would misattribute this global re-export
    as a doxorubicin-local effect.

Both columns are reported so the distinction is visible rather than assumed.

PROVENANCE OF EVERY ALTERNATIVE logP IN THIS FILE
-------------------------------------------------
No number below is invented. Each scenario is tagged with its tier:
  measured        - a real octanol/water shake-flask log Kow
  measured_logD   - a real OECD 107 logD at a stated pH (NOT a logP)
  refuted         - a value that could not be found in the document cited for it
  propagation     - another molecule's measured value, quoted for this one
  prediction      - a calculator output (RDKit Crippen, PHYSPROP/OPERA, etc.)
  unsourced       - a plausible-error band, used only where no number exists

Run:  ./.venv/bin/python experiments/logp_sensitivity.py
"""
from __future__ import annotations

import math
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

from tmesim.molecular import binding as _binding
from tmesim.molecular import breast_library, library, permeability
from tmesim.molecular.compile_drug import molecule_to_drugparams
from tmesim.physics.diffusion import ALPHA_TORTUOSITY
from tmesim.validation.calibration import DOX_TARGET_EFOLD_UM
from tmesim.validation.pk_reality import DEFAULT_ECM

ANCHOR = "doxorubicin"


# ---------------------------------------------------------------------------
# Scenario table. Every entry: (logp_value, tier, label).
# ---------------------------------------------------------------------------
# Baselines are filled in at runtime from the repo's own RDKit Crippen value so
# this file cannot drift from what tmesim actually computes.
SCENARIOS: Dict[str, List[Tuple[float, str, str]]] = {
    "cisplatin": [
        (-2.53, "measured", "Screnci 2000 Br J Cancer 82(4):966-72, PMID 10732773, "
                            "Table 3, shake-flask n=5 (octanol/0.9% SALINE, total Pt "
                            "by ICP-MS -- not textbook octanol/water)"),
        (-2.25, "measured", "Screnci 2000 upper 1-SD bound (-2.53 + 0.28)"),
        (-2.81, "measured", "Screnci 2000 lower 1-SD bound (-2.53 - 0.28)"),
    ],
    "palbociclib": [
        (0.99, "refuted", "DrugBank value attributed to the IBRANCE monograph; "
                          "absent from Health Canada monograph rev. 2024-09-04 s.13 "
                          "and from FDA label 207103s008"),
        (0.228, "measured_logD", "EMA Ibrance EPAR EMA/652627/2016 p.28, OECD 107, pH 5"),
        (1.11, "measured_logD", "EMA Ibrance EPAR EMA/652627/2016 p.28, OECD 107, pH 7"),
        (2.26, "measured_logD", "EMA Ibrance EPAR EMA/652627/2016 p.28, OECD 107, pH 9"),
    ],
    "epirubicin": [
        (1.27, "propagation", "doxorubicin's HSDB log Kow (Hansch p.187) quoted in "
                              "epirubicin narrative text; epirubicin is the 4'-epimer "
                              "and has NO Hansch-backed field of its own"),
    ],
    "docetaxel": [
        (4.26, "prediction", "EMA EPAR figure; the two EPARs proved to be the same "
                             "ASMF boilerplate quoted twice, no stated method -- "
                             "downgraded mid-audit"),
        (3.2, "prediction", "Sanofi Taxotere SDS, conflicts with 4.26, no method"),
    ],
    "tamoxifen": [
        (7.1, "prediction", "PHYSPROP/OPERA 'experimental' field, verified from the "
                            "OPERA/PHYSPROP files to be a PREDICTION laundered into "
                            "an experimental slot"),
        (7.87, "prediction", "El-Kattan 2001 PMID 11250108; real text but a calculator "
                             "output -- the same table reports nicardipine "
                             "HYDROCHLORIDE logP -0.99, impossible for a salt"),
    ],
    "doxorubicin": [
        (1.27, "measured", "HSDB canonical field 'log Kow = 1.27' citing Hansch p.187 "
                           "(PubChem CID 31703) -- THE CALIBRATION ANCHOR"),
    ],
    # No usable retrieved number for these four: paclitaxel has no Hansch-backed
    # field at all; vinorelbine, irinotecan and ribociclib have real OECD 107 logD
    # but the specific pH-resolved values were not retrieved in this audit, so
    # inventing them here would be exactly the failure this audit exists to stop.
    "paclitaxel": [],
    "vinorelbine": [],
    "irinotecan": [],
    "ribociclib": [],
}

#: Applied to every drug so a uniform "what if this is off by 1 log unit"
#: comparison exists even where no alternative literature number was retrieved.
GENERIC_BAND = (-1.0, +1.0)

AUDITED = [
    "cisplatin", "paclitaxel", "docetaxel", "epirubicin", "vinorelbine",
    "tamoxifen", "irinotecan", "palbociclib", "ribociclib",
]


# ---------------------------------------------------------------------------
# Compile once per drug (this is the only expensive step).
# ---------------------------------------------------------------------------
def _molecule(name: str):
    """Resolve a drug name across the two disjoint libraries."""
    if name in library.GBM_DRUGS:
        return library.get_molecule(name)
    if name in breast_library.BREAST_DRUGS:
        return breast_library.get_molecule_breast(name)
    raise KeyError(f"{name!r} is in neither GBM_DRUGS nor BREAST_DRUGS")


def compile_once(name: str) -> Dict:
    """Compile a drug and cache the logP-INDEPENDENT parts.

    ``diffusivity`` is the only expensive quantity and it does not depend on
    logP, so it is computed once and reused for every scenario.
    """
    mol = _molecule(name)
    dp = molecule_to_drugparams(mol, ecm_reference=DEFAULT_ECM)
    d_eff = float(dp.diffusivity) * (
        1.0 - ALPHA_TORTUOSITY * float(np.clip(DEFAULT_ECM, 0.0, 1.0))
    )
    return {
        "name": name,
        "descriptors": dict(mol.descriptors),
        "source": mol.source,
        "baseline_logp": float(mol.descriptors["logp"]),
        "diffusivity": float(dp.diffusivity),
        "d_eff": max(d_eff, 1e-12),
    }


def evaluate(rec: Dict, logp: float, alpha: float) -> Dict:
    """Recompute every logP-dependent quantity in closed form."""
    d = dict(rec["descriptors"])
    d["logp"] = float(logp)
    bind = _binding.tissue_binding_rate(d)
    clr = _binding.clearance_rate(d)
    bbb = permeability.bbb_penetrance(d)
    perm = permeability.membrane_permeability(d)
    ff = _binding.plasma_free_fraction(d)
    k = clr + alpha * bind
    L = math.sqrt(rec["d_eff"] / k) if k > 0 else float("inf")
    return {
        "logp": float(logp), "binding": bind, "clearance": clr, "k_eff": k,
        "L_um": L, "bbb": bbb, "perm": perm, "free_fraction": ff,
    }


def alpha_for(dox_rec: Dict, dox_logp: float) -> float:
    """Refit the single calibration constant on doxorubicin at a given logP.

    Mirrors ``calibration.calibrate_alpha`` exactly:
        k_target = D_eff / L*^2 ;  alpha = (k_target - clearance) / binding
    """
    d = dict(dox_rec["descriptors"])
    d["logp"] = float(dox_logp)
    proxy = _binding.tissue_binding_rate(d)
    kp = _binding.clearance_rate(d)
    k_target = dox_rec["d_eff"] / (DOX_TARGET_EFOLD_UM ** 2)
    return (k_target - kp) / proxy


# ---------------------------------------------------------------------------
def main() -> int:
    print("=" * 100)
    print("logP SENSITIVITY SWEEP")
    print("=" * 100)
    print(f"anchor drug              : {ANCHOR} (the ONLY drug alpha is fit on)")
    print(f"fit target (e-fold L*)   : {DOX_TARGET_EFOLD_UM:.4f} um")
    print(f"ECM density              : {DEFAULT_ECM}")
    print()

    names = sorted(set(AUDITED) | {ANCHOR})
    print("compiling (one Brownian sim per drug; logP does not affect this step)...")
    recs: Dict[str, Dict] = {}
    for n in names:
        recs[n] = compile_once(n)
        print(f"  {n:16s} src={recs[n]['source']:14s} "
              f"logP(Crippen)={recs[n]['baseline_logp']:+7.3f}  "
              f"D={recs[n]['diffusivity']:7.2f} um^2/s")
    print()

    alpha0 = alpha_for(recs[ANCHOR], recs[ANCHOR]["baseline_logp"])
    print(f"alpha at baseline        : {alpha0:.6f}")
    print(f"(repo DRUG_LIBRARY froze : 148.60815765445545)")
    print()

    # --- SANITY: does the analytic route reproduce the frozen library rows? ---
    print("-" * 100)
    print("SANITY CHECK vs the frozen tmesim.core.config.DRUG_LIBRARY snapshots")
    print("-" * 100)
    from tmesim.core import DRUG_LIBRARY
    ok = True
    for lib_name, drug in [("palbociclib-like", "palbociclib"),
                           ("ribociclib-like", "ribociclib"),
                           ("irinotecan-like", "irinotecan")]:
        if lib_name not in DRUG_LIBRARY or drug not in recs:
            continue
        frozen = DRUG_LIBRARY[lib_name]
        got = evaluate(recs[drug], recs[drug]["baseline_logp"], 148.60815765445545)
        # The frozen row stores the CALIBRATED sink, i.e. alpha * tissue_binding_rate.
        bind_cal = 148.60815765445545 * got["binding"]
        # Some library rows were written to 4 significant figures (irinotecan),
        # others to full precision (palbociclib/ribociclib) -- use a relative
        # tolerance that accepts 4-sig-fig rounding.
        def _rel(a, b):
            return abs(a - b) / max(abs(b), 1e-30)
        db = _rel(bind_cal, frozen.binding)
        dc = _rel(got["clearance"], frozen.clearance)
        dbbb = _rel(got["bbb"], frozen.bbb_penetrance)
        good = db < 1e-3 and dc < 1e-3 and dbbb < 1e-3
        ok &= good
        print(f"  {lib_name:18s} alpha*binding {bind_cal:.12f} vs {frozen.binding:.12f} "
              f"| clearance {got['clearance']:.4e} vs {frozen.clearance:.4e} "
              f"| bbb {got['bbb']:.6f} vs {frozen.bbb_penetrance:.6f}  "
              f"-> {'MATCH' if good else 'MISMATCH'}")
    print(f"  => the frozen library rows ARE alpha * tissue_binding_rate(Crippen logP): "
          f"{'CONFIRMED' if ok else 'NOT CONFIRMED'}")
    print()

    # --- Per-drug scenario sweep -------------------------------------------
    print("-" * 100)
    print("PER-DRUG SCENARIOS  (L = penetration length, um; bbb = delivery scalar)")
    print("  'alpha fixed' = alpha frozen at the baseline fit")
    print("  'alpha refit' = calibrate_alpha() re-run under the SAME perturbation")
    print("-" * 100)
    hdr = (f"{'drug':13s} {'logP':>7s} {'tier':>14s} {'k_eff':>10s} "
           f"{'L_fix':>8s} {'dL_fix':>8s} {'L_ref':>8s} {'dL_ref':>8s} "
           f"{'bbb':>7s} {'dbbb':>8s}")
    rows_out = []
    for name in AUDITED + [ANCHOR]:
        rec = recs[name]
        base_lp = rec["baseline_logp"]
        base = evaluate(rec, base_lp, alpha0)
        print()
        print(f"### {name}   (baseline = RDKit Crippen {base_lp:+.3f}, "
              f"L={base['L_um']:.2f} um, bbb={base['bbb']:.4f})")
        print(hdr)
        scen = list(SCENARIOS.get(name, []))
        scen += [(base_lp + GENERIC_BAND[0], "unsourced", "generic -1 log unit"),
                 (base_lp + GENERIC_BAND[1], "unsourced", "generic +1 log unit")]
        for lp, tier, note in scen:
            a_ref = alpha_for(recs[ANCHOR], lp) if name == ANCHOR else alpha0
            fix = evaluate(rec, lp, alpha0)
            ref = evaluate(rec, lp, a_ref)
            print(f"{'':13s} {lp:+7.3f} {tier:>14s} {fix['k_eff']:10.3e} "
                  f"{fix['L_um']:8.2f} {fix['L_um']/base['L_um']:7.3f}x "
                  f"{ref['L_um']:8.2f} {ref['L_um']/base['L_um']:7.3f}x "
                  f"{fix['bbb']:7.4f} {fix['bbb']-base['bbb']:+8.4f}")
            print(f"{'':13s}   ^ {note}")
            rows_out.append((name, lp, tier, fix["L_um"] / base["L_um"],
                             ref["L_um"] / base["L_um"], fix["bbb"] - base["bbb"]))

    # --- The anchor's global re-export -------------------------------------
    print()
    print("-" * 100)
    print("ANCHOR RE-EXPORT: correcting DOXORUBICIN's logP to the measured 1.27")
    print("  (this is the one perturbation that moves EVERY other drug)")
    print("-" * 100)
    a_new = alpha_for(recs[ANCHOR], 1.27)
    print(f"  alpha {alpha0:.4f} -> {a_new:.4f}   ({a_new/alpha0:.4f}x)")
    dox_base = evaluate(recs[ANCHOR], recs[ANCHOR]["baseline_logp"], alpha0)
    dox_new = evaluate(recs[ANCHOR], 1.27, a_new)
    print(f"  doxorubicin's OWN L: {dox_base['L_um']:.3f} -> {dox_new['L_um']:.3f} um "
          f"(pinned to the fit target by construction -- NOT a prediction)")
    print()
    print(f"  {'drug':16s} {'L_before':>10s} {'L_after':>10s} {'ratio':>8s}")
    for name in AUDITED:
        rec = recs[name]
        b = evaluate(rec, rec["baseline_logp"], alpha0)
        a = evaluate(rec, rec["baseline_logp"], a_new)
        print(f"  {name:16s} {b['L_um']:10.2f} {a['L_um']:10.2f} "
              f"{a['L_um']/b['L_um']:7.3f}x")

    # --- Composite: delivery x penetration ----------------------------------
    # bbb_penetrance is a LINEAR multiplier on the delivery source
    # (treatment.drugs.delivery_source: S = vascular_density * bbb_penetrance
    # * plasma), and L sets the exponential decay length. The drug delivered per
    # unit vessel length at steady state therefore scales as
    #     integral_0^inf  C(0) * exp(-r/L) dr   ~   bbb * L
    # so the two routes MULTIPLY. This is the honest blast-radius number.
    print()
    print("-" * 100)
    print("COMPOSITE BLAST RADIUS: (bbb x L) relative to the Crippen baseline")
    print("  bbb is a linear multiplier on delivery; L is the decay length; they multiply.")
    print("-" * 100)
    print(f"  {'drug':13s} {'logP':>7s} {'tier':>14s} {'bbb_x':>8s} {'L_x':>8s} {'bbb*L_x':>9s}")
    for name in AUDITED:
        rec = recs[name]
        base = evaluate(rec, rec["baseline_logp"], alpha0)
        for lp, tier, _ in SCENARIOS.get(name, []):
            v = evaluate(rec, lp, alpha0)
            br = v["bbb"] / base["bbb"] if base["bbb"] > 0 else float("nan")
            lr = v["L_um"] / base["L_um"]
            print(f"  {name:13s} {lp:+7.3f} {tier:>14s} {br:7.3f}x {lr:7.3f}x "
                  f"{br * lr:8.3f}x")

    # --- Headline summary ---------------------------------------------------
    print()
    print("-" * 100)
    print("HEADLINE")
    print("-" * 100)
    worst = max(rows_out, key=lambda r: abs(math.log10(r[3])) if r[3] > 0 else 0)
    print(f"  largest |L| move from a NON-anchor logP change: {worst[0]} "
          f"logP={worst[1]:+.3f} ({worst[2]}) -> {worst[3]:.3f}x")
    print(f"  alpha-fixed and alpha-refit columns are IDENTICAL for every "
          f"non-anchor drug (alpha depends only on {ANCHOR}).")
    print(f"  logP CANNOT move diffusivity: D_eff is a function of MW/radius only.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
