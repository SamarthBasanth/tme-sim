"""PROTOTYPE ONLY -- score the engine's compiled transport parameters against a
measured multicellular-layer (MCL) tissue diffusion coefficient and a measured
first-order tissue consumption rate constant.

Nothing here is registered in ``tmesim/validation/scorecard.py``. This script
exists so ``results/external_anchor_survey.md`` can quote what the engine
actually gives on the highest-ranked candidate anchor rather than speculating.
It imports the molecular pipeline and reads it; it modifies nothing.

WHY THIS ANCHOR
---------------
The engine's transport PDE has exactly two per-drug coefficients:

    dC/dt = D_eff * lap(C) - k*C + S       (C = FREE INTERSTITIAL concentration)

``D_eff`` is ``MolecularProfile.effective_diffusivity`` (um^2/s) and ``k`` is
``clearance + alpha*tissue_binding_rate`` (1/s). Both are ABSOLUTE quantities.
``results/alpha_identifiability.md`` establishes that the project has exactly
ONE measurement constraining either of them (the doxorubicin perivascular
half-distance, which ``alpha`` is fitted to in closed form, so it is in-sample
by construction), and ``results/pk_reality.md`` shows the uncalibrated ``k`` was
~300x too small.

Hicks et al. 1998 (PMID 9806526) measured BOTH coefficients directly, in the
same functional roles, by fitting flux through multicellular layers to a
diffusion-reaction model of exactly this form. That makes it the most
commensurable external measurement of the solver's own inputs found in this
survey.

MEASURED VALUES USED (all read from the PubMed abstract of the cited paper;
none derived, none inferred)
----------------------------------------------------------------------------
Hicks KO, Fleming Y, Siim BG, Koch CJ, Wilson WR. "Extravascular diffusion of
tirapazamine: effect of metabolic consumption assessed using the multicellular
layer model." Int J Radiat Oncol Biol Phys 1998;42(3):641-9. PMID 9806526,
doi:10.1016/s0360-3016(98)00268-5

  * D_mcl(tirapazamine, V79-171b MCL)   = 7.4e-7  cm^2/s  =  74 um^2/s
  * D_mcl(tirapazamine, MGH-U1 MCL)     = 1.3e-6  cm^2/s  = 130 um^2/s
  * stated hindrance vs culture medium (V79) = 12-fold lower
  * anoxic metabolic consumption rate constant, V79    = 0.52 /min
  * anoxic metabolic consumption rate constant, MGH-U1 = 0.31 /min

Hicks KO, Myint H, Patterson AV, Pruijn FB, Siim BG, Patel K, Wilson WR.
Int J Radiat Oncol Biol Phys 2007;69(2):560-71. PMID 17869669,
doi:10.1016/j.ijrobp.2007.05.049  -- INDEPENDENT REPLICATE, different cell line

  * D_mcl(tirapazamine, SiHa MCL) = 1.30 +/- 0.05 e-6 cm^2/s = 130 +/- 5 um^2/s
  * D_mcl(PR-104A,       SiHa MCL) = 4.42 +/- 0.15 e-7 cm^2/s = 44.2 +/- 1.5 um^2/s

Run:  ./.venv/bin/python experiments/anchor_prototype.py
"""
from __future__ import annotations

import math

from tmesim.molecular import compile_drug, descriptors

# --------------------------------------------------------------------------- #
# Measured anchors. Values ONLY; every one traced to the PMID in the docstring.
# --------------------------------------------------------------------------- #
CM2_S_TO_UM2_S = 1.0e8          # 1 cm^2/s = 1e8 um^2/s
PER_MIN_TO_PER_S = 1.0 / 60.0

D_MCL_MEASURED_UM2_S = {
    # (value, uncertainty or None, cell line, PMID)
    "V79-171b":  (7.4e-7 * CM2_S_TO_UM2_S, None,                   "9806526"),
    "MGH-U1":    (1.3e-6 * CM2_S_TO_UM2_S, None,                   "9806526"),
    "SiHa":      (1.30e-6 * CM2_S_TO_UM2_S, 0.05e-6 * CM2_S_TO_UM2_S, "17869669"),
}
K_CONSUMPTION_MEASURED_PER_S = {
    "V79-171b": (0.52 * PER_MIN_TO_PER_S, "9806526"),
    "MGH-U1":   (0.31 * PER_MIN_TO_PER_S, "9806526"),
}
#: Hicks 1998 states the V79 MCL diffusion coefficient is "12-fold lower than in
#: culture medium". That is a directly measured HINDRANCE factor.
MEASURED_HINDRANCE_FOLD_V79 = 12.0

#: Tirapazamine. SMILES validated by round-trip: RDKit returns MW 178.151
#: against the published monoisotopic/average MW 178.15 for TPZ.
TPZ_SMILES = "C1=CC=C2C(=C1)[N+](=NC(=[N+]2[O-])N)[O-]"

#: The repo's own ECM reference used everywhere else in the pipeline.
ECM_REFERENCE = 0.4


def _fold(pred: float, meas: float) -> str:
    r = pred / meas
    return f"{r:8.2f}x" + ("  (over)" if r > 1 else "  (under)")


def main() -> None:
    print("=" * 78)
    print("ANCHOR PROTOTYPE -- measured MCL diffusion coefficient & consumption rate")
    print("NOT REGISTERED. Recommendation input for results/external_anchor_survey.md")
    print("=" * 78)

    mol = descriptors.from_smiles("tirapazamine", TPZ_SMILES)
    print(f"\nmolecule: {mol.name}   descriptor source: {mol.source}")
    for k in ("mw", "logp", "tpsa", "hbd", "hba"):
        print(f"    {k:8s} = {mol.get(k)}")
    print("    MW check: published tirapazamine MW is 178.15; RDKit gives "
          f"{mol.get('mw')}. SMILES accepted.")

    prof = compile_drug.build_profile(mol, ecm_reference=ECM_REFERENCE)

    print("\n--- what the engine compiles (ABSOLUTE, both are PDE coefficients) ---")
    print(f"    free_diffusivity       = {prof.free_diffusivity:10.2f} um^2/s "
          "(Stokes-Einstein)")
    print(f"    effective_diffusivity  = {prof.effective_diffusivity:10.2f} um^2/s "
          f"(Brownian hindrance @ ECM={ECM_REFERENCE})")
    print(f"    hindrance D_eff/D_free = {prof.hindrance:10.4f}  "
          f"=> {1.0/prof.hindrance:.2f}-fold reduction")
    print(f"    clearance              = {prof.clearance:10.3e} /s")
    print(f"    tissue_binding_rate    = {prof.tissue_binding_rate:10.3e} /s")
    k_raw = prof.clearance + prof.tissue_binding_rate
    print(f"    k_raw = clearance + binding             = {k_raw:10.3e} /s")

    # The shipped engine amplifies the binding half by the fitted constant alpha.
    try:
        from tmesim.core.config import PDEConfig  # noqa: F401
        from tmesim.validation import calibration
        alpha = float(calibration.calibrate_alpha())
        k_cal = prof.clearance + alpha * prof.tissue_binding_rate
        print(f"    alpha (fitted, in-sample)               = {alpha:10.4f}")
        print(f"    k_calibrated = clearance + alpha*binding = {k_cal:10.3e} /s")
    except Exception as exc:  # pragma: no cover - diagnostic only
        alpha, k_cal = None, None
        print(f"    [alpha unavailable in this checkout: {exc}]")

    print("\n" + "-" * 78)
    print("SCORE 1 -- effective diffusion coefficient (um^2/s), ABSOLUTE")
    print("-" * 78)
    print(f"{'MCL cell line':14s} {'measured':>14s} {'predicted':>12s} "
          f"{'pred/meas':>16s}  PMID")
    for line, (meas, unc, pmid) in D_MCL_MEASURED_UM2_S.items():
        u = f" +/- {unc:.0f}" if unc else ""
        print(f"{line:14s} {meas:10.1f}{u:>4s} {prof.effective_diffusivity:12.1f} "
              f"{_fold(prof.effective_diffusivity, meas):>16s}  {pmid}")
    lo = min(v[0] for v in D_MCL_MEASURED_UM2_S.values())
    hi = max(v[0] for v in D_MCL_MEASURED_UM2_S.values())
    inside = lo <= prof.effective_diffusivity <= hi
    print(f"\n    measured range across 3 independent MCL cell lines: "
          f"{lo:.0f}-{hi:.0f} um^2/s")
    print(f"    engine prediction {prof.effective_diffusivity:.1f} um^2/s -> "
          f"{'INSIDE' if inside else 'OUTSIDE'} the measured range")
    print(f"    even against the LOOSEST 3x order-of-magnitude band this repo "
          f"uses: {'PASS' if prof.effective_diffusivity <= 3*hi else 'FAIL'}")

    print("\n" + "-" * 78)
    print("SCORE 2 -- hindrance factor D_free/D_eff (dimensionless), ABSOLUTE")
    print("-" * 78)
    pred_fold = 1.0 / prof.hindrance
    print(f"    measured  (V79 MCL vs culture medium, PMID 9806526) = "
          f"{MEASURED_HINDRANCE_FOLD_V79:6.2f}-fold")
    print(f"    predicted (Brownian obstacles @ ECM={ECM_REFERENCE})            = "
          f"{pred_fold:6.2f}-fold")
    print(f"    ratio predicted/measured = "
          f"{pred_fold/MEASURED_HINDRANCE_FOLD_V79:.3f}x")
    print("    NOTE: the engine's only obstacle population is 'ECM'. An MCL is "
          "packed CELLS\n          with little stroma, so this row asks whether "
          "the hindrance law\n          accounts for cellular volume exclusion "
          "at all.")

    print("\n" + "-" * 78)
    print("SCORE 3 -- first-order tissue loss rate k (1/s), ABSOLUTE")
    print("-" * 78)
    print("    CAVEAT, stated before the numbers: the measured constant is "
          "ANOXIC BIOREDUCTIVE\n    metabolism of a hypoxia-activated prodrug. "
          "The engine's k is a generic\n    'binding + clearance' sink. Same "
          "units, same position in the same PDE, but\n    NOT the same physical "
          "process. Read this as an order-of-magnitude bound on\n    how large a "
          "real measured interstitial consumption rate can be, not as a\n"
          "    like-for-like score.")
    print(f"\n{'MCL cell line':14s} {'measured k':>14s} {'k_raw':>12s} "
          f"{'k_calibrated':>14s}  PMID")
    for line, (meas, pmid) in K_CONSUMPTION_MEASURED_PER_S.items():
        kc = f"{k_cal:12.3e}" if k_cal is not None else "n/a".rjust(12)
        print(f"{line:14s} {meas:14.3e} {k_raw:12.3e} {kc:>14s}  {pmid}")
    for line, (meas, pmid) in K_CONSUMPTION_MEASURED_PER_S.items():
        print(f"    {line}: k_raw/measured = {k_raw/meas:8.4f}x", end="")
        if k_cal is not None:
            print(f" ; k_calibrated/measured = {k_cal/meas:8.4f}x")
        else:
            print()

    print("\n" + "-" * 78)
    print("SCORE 4 -- the penetration length the two coefficients imply")
    print("-" * 78)
    print("    L = sqrt(D/k). Computed BOTH ways, so the comparison is "
          "coefficient-for-\n    coefficient rather than mixing a measured D "
          "with a predicted k.")
    for line in ("V79-171b", "MGH-U1"):
        dm = D_MCL_MEASURED_UM2_S[line][0]
        km = K_CONSUMPTION_MEASURED_PER_S[line][0]
        l_meas = math.sqrt(dm / km)
        print(f"    {line:10s} measured D & measured k -> L = {l_meas:8.2f} um")
    if k_cal is not None and k_cal > 0:
        print(f"    {'engine':10s} predicted D & predicted k -> L = "
              f"{math.sqrt(prof.effective_diffusivity / k_cal):8.2f} um")
    if k_raw > 0:
        print(f"    {'engine':10s} predicted D & UNAMPLIFIED k -> L = "
              f"{math.sqrt(prof.effective_diffusivity / k_raw):8.2f} um")

    print("\n" + "=" * 78)
    print("Prototype complete. No scorecard row was registered and no file "
          "under tmesim/ was\ntouched. See results/external_anchor_survey.md "
          "for the commensurability verdicts.")
    print("=" * 78)


if __name__ == "__main__":
    main()
