"""Passive membrane permeability and blood-brain-barrier (BBB) penetrance estimators.

Both functions are structure-only, physicochemical estimates. They take a
`descriptors` dict (the keys produced by `molecular.descriptors`: mw, logp, tpsa,
hbd, hba, ..., is_antibody) and return single physical numbers the tissue engine
can consume.

Honesty line (see MOLECULAR_CONTRACT.md): permeability and BBB crossing are NOT
reliably predictable from 2D descriptors alone. Active transport, efflux pumps
(P-gp), carrier-mediated uptake, and metabolism all move real drugs off the
physicochemical trend line. These estimators reproduce the *passive* solubility-
diffusion physics and the well-known CNS-drug "sweet spot" in descriptor space,
and they are documented where they are known to mis-rank real drugs (temozolomide
in particular — see `bbb_penetrance`).

KNOWN MISSING MECHANISM — ionisation and sequestration (pH partitioning)
------------------------------------------------------------------------
Both functions here are **neutral-species** estimates. Neither takes a pKa, so
neither can express pH-dependent partitioning, and nothing downstream of them
represents the acidic-compartment ion-trapping that concentrates weak bases in
tissue. Two independent validations have now converged on this same gap:

* ``results/gist_imatinib.md`` — measured tumour:normal ratios for imatinib in
  human liver metastases are 0.023-0.142, BELOW the transport model's structural
  floor of 0.208. No parameter choice reaches them. Hepatic uptake transporters
  and **lysosomal ion-trapping of a weak base** (imatinib, basic pKa 8.07) in
  normal liver are the named missing mechanisms.
* ``results/moncayo_explant.md`` — independently tests the same effect in HGSOC
  explants, where the weak bases niraparib/rucaparib sequester and the neutral
  olaparib does not.

This is **deliberately not implemented yet**: a literature sweep established
there is no measured tumour acidic-phospholipid content in usable units and no
tumour-specific lysosomal parameters in any published PBPK model, so a
pH-partitioning term would have to be designed rather than parameterised. The
note is here so the limitation is visible where the code is read. The structure
above does not foreclose it — a neutral-fraction factor multiplies ``Kp`` and a
compartment term composes downstream of the series resistance in
``membrane_permeability``.
"""
from __future__ import annotations

import math
from typing import Dict

# Physical constants (SI unless noted)
_KB = 1.380649e-23  # Boltzmann constant, J/K
_ETA_MEMBRANE_PA_S = 1.0  # effective viscosity inside a lipid bilayer core, Pa*s.
# The bilayer interior is far more viscous than water (~7e-4 Pa*s at 310 K);
# membrane-core microviscosities of ~0.1-1 Pa*s are commonly reported. We take the
# TOP of that reported range. It used to be 0.1 (the bottom), which put the
# prefactor Dm/d at ~2.5e3 um/s, i.e. every drug with Kp >= 1 already exceeded any
# physically achievable permeability before its lipophilicity was even applied.
_BILAYER_THICKNESS_NM = 4.0  # d, hydrophobic bilayer thickness, nm (~4 nm typical)

# --- Ceilings on the passive estimate (see `membrane_permeability`) ----------
# Upper clamp on the logP used to build the membrane/water partition coefficient.
# NOT a numerical guard: octanol/water and BILAYER/water partitioning decouple
# above logP ~ 3 (the bilayer is a more ordered, more polar solvent than bulk
# octanol, so liposome/water partitioning plateaus while octanol/water keeps
# rising). Extrapolating 10**logP past that point is the single biggest source of
# absurd permeabilities. Independent support from our own runs: RDKit's Crippen
# logP for imatinib is 4.59 against an experimental octanol/water ~3.0
# (`results/gist_imatinib.md`), i.e. the estimator itself runs ~1.6 log units
# lipophilic exactly where this clamp bites.
_LOGP_KP_MAX = 3.0
_LOGP_KP_MIN = -3.0

# Unstirred aqueous boundary layer (UWL) thickness, um. A solute must cross a
# stagnant water film before it ever reaches the bilayer, so the two resistances
# add in series and the UWL sets the HARD physical ceiling on passive
# permeability, no matter how lipophilic the molecule. Reported UWL thicknesses
# span ~30 um (well-stirred / in vivo) to >1000 um (unstirred Caco-2 wells); we
# take 30 um, the THIN end, so the ceiling is as permissive as the literature
# allows and only removes values that are absurd rather than merely high.
_UWL_THICKNESS_UM = 30.0
_ETA_WATER_PA_S = 6.9e-4  # water viscosity at 310 K, Pa*s

# Numerical floor on a single descriptor desirability inside the geometric mean in
# `bbb_penetrance`. Only guards log(0) when a sigmoid underflows; it is orders of
# magnitude below anything the output logistic can resolve.
_DESIRABILITY_FLOOR = 1e-12


def _radius_nm_from_descriptors(descriptors: Dict[str, float]) -> float:
    """Solute radius in nm.

    Prefer a precomputed hydrodynamic radius; otherwise estimate from MW assuming
    a spherical solute at ~1.0 g/cm^3 packing density:
        V = mw / (rho * N_A)  ->  r = (3V / 4pi)^(1/3)
    Small organics come out near 0.3-0.5 nm, which is the right ballpark.
    """
    r = descriptors.get("radius_hydrodynamic_nm")
    if r:
        return float(r)
    mw = float(descriptors.get("mw", 300.0))
    rho = 1.0e3  # kg/m^3 (1 g/cm^3)
    n_a = 6.02214076e23
    mass_kg = (mw * 1e-3) / n_a  # g/mol -> kg per molecule
    vol_m3 = mass_kg / rho
    r_m = (3.0 * vol_m3 / (4.0 * math.pi)) ** (1.0 / 3.0)
    return r_m * 1e9  # m -> nm


def membrane_permeability(descriptors: Dict[str, float], temperature_K: float = 310.0) -> float:
    """Passive permeability across a lipid bilayer, in **um/s**.

    Solubility-diffusion model (Overton / Meyer picture) in series with the
    unstirred aqueous boundary layer:

        P_membrane = Kp * Dm / d
        P_aqueous  = D_water / h_UWL
        1 / P      = 1 / P_membrane + 1 / P_aqueous

    where
      * Kp  : membrane/water partition coefficient (dimensionless), from
              ``Kp = 10**logP`` with logP clamped into
              ``[_LOGP_KP_MIN, _LOGP_KP_MAX]`` = [-3, 3]. The UPPER clamp is
              physical, not numerical — see :data:`_LOGP_KP_MAX`.
      * Dm  : diffusion coefficient of the solute *inside* the bilayer core,
              m^2/s, Stokes-Einstein with a membrane-core viscosity:
                  Dm = kB*T / (6*pi*eta_membrane*r)
      * d   : bilayer hydrophobic thickness, ~4 nm.
      * D_water : Stokes-Einstein in water at ``temperature_K``.
      * h_UWL   : unstirred water layer thickness, :data:`_UWL_THICKNESS_UM`.

    Units bookkeeping:
        Dm [m^2/s] / d [m]  ->  [m/s];  * Kp (dimensionless)  ->  [m/s].
        Convert m/s -> um/s by * 1e6.

    Antibodies / macromolecules (`is_antibody`) do not cross by passive diffusion:
    return ~0.

    WHY THERE IS A CEILING AT ALL (a measured defect, not a style preference)
    ------------------------------------------------------------------------
    The previous form was ``P = Kp * Dm / d`` with logP clamped only at 5 and a
    membrane viscosity of 0.1 Pa*s. It returned **4.23e+07 um/s for imatinib** —
    42 metres per second — because it exponentiates a Crippen logP of 4.59 with
    no bound (`results/gist_imatinib.md` item 7). Nothing in the solubility-
    diffusion equation stops that: it is a bare exponential in a descriptor.

    This is the SAME root pathology as the eribulin failure in
    :func:`bbb_penetrance` — unbounded or unshaped use of a single descriptor —
    and it is fixed the same way, by writing down the physics that actually
    limits the quantity:

    1. **Kp cannot follow octanol/water forever.** Bilayer/water partitioning
       plateaus above logP ~ 3 while octanol/water keeps climbing, so the logP
       driving Kp is clamped at 3 (:data:`_LOGP_KP_MAX`).
    2. **A solute must cross the stagnant water film first.** The unstirred
       aqueous boundary layer is in SERIES with the membrane, so it caps P at
       ``D_water / h_UWL`` (~20-30 um/s for a small organic) no matter how
       lipophilic the molecule is. This is the real physical ceiling on passive
       permeability and it is why measured Caco-2/PAMPA permeabilities saturate.
    3. **The membrane prefactor was at the permissive end.** ``eta_membrane``
       moved from 0.1 to 1.0 Pa*s, the top of the 0.1-1 Pa*s range this module
       already cited. At 0.1 the prefactor ``Dm/d`` alone was ~2.5e3 um/s, so
       every drug with Kp >= 1 was already above any achievable permeability and
       the boundary-layer ceiling would have flattened the entire library to a
       single value, destroying the cross-drug ordering that is the only thing
       this estimator is actually used for.

    HONEST LIMITATION — the ceiling binds, and the top of the ranking is
    compressed. Cisplatin (22.2 um/s), imatinib (21.0) and methotrexate (19.6)
    now sit within 15% of each other, because all three are boundary-layer
    limited under this parameterisation. Physically they differ by ORDERS of
    magnitude in measured apparent permeability (cisplatin is famously almost
    membrane-impermeant). The ordering is no longer absurd, but it is no longer
    informative at the top either. Two separate residual defects cause this and
    neither is fixed here: the membrane term is still too permissive in absolute
    terms (real Caco-2 Papp for these drugs is 0.001-1 um/s, i.e. 1-4 orders
    below what this returns), and Crippen logP is simply wrong for some species
    (cisplatin's +0.20 vs an experimental ~-2.2). Treat the output as a bounded,
    correctly-ordered-at-the-bottom RANKING, not as a calibrated rate.

    HONEST LIMITATION — no ionisation, no sequestration. This is a neutral-species
    passive estimate. It has no pKa term, so it cannot express pH-partitioning,
    and no lysosomal/acidic-compartment term, so it cannot express the
    ion-trapping that concentrates weak bases in normal tissue. Two independent
    validations now point at that same gap: the GIST imatinib run
    (`results/gist_imatinib.md`) found measured tumour:normal liver ratios of
    0.023-0.142 lying BELOW the transport model's structural floor of 0.208 —
    unreachable at ANY parameterisation, with hepatic uptake transporters and
    lysosomal ion-trapping of a weak base named as the missing mechanisms — and
    the Moncayo explant run is independently testing the same effect (weak bases
    niraparib/rucaparib sequester, neutral olaparib does not). Nothing in this
    function's structure forecloses adding a pH-partitioning factor later: it
    would multiply ``Kp`` by the neutral fraction at the relevant pH and add a
    compartment term downstream, both of which compose with the series
    resistance above.

    Returns a non-negative float in um/s.
    """
    if descriptors.get("is_antibody"):
        return 0.0

    # --- Kp from logP, clamped (upper clamp is PHYSICAL, see _LOGP_KP_MAX) ---
    logp = float(descriptors.get("logp", 0.0))
    logp_clamped = max(_LOGP_KP_MIN, min(_LOGP_KP_MAX, logp))
    kp = 10.0 ** logp_clamped  # dimensionless, in [1e-3, 1e3]

    # --- Dm from size (Stokes-Einstein in the membrane core) ---
    r_nm = _radius_nm_from_descriptors(descriptors)
    r_m = max(r_nm, 0.1) * 1e-9  # guard against tiny/zero radius
    dm_m2s = (_KB * temperature_K) / (6.0 * math.pi * _ETA_MEMBRANE_PA_S * r_m)

    # --- Membrane leg: P = Kp * Dm / d ---
    d_m = _BILAYER_THICKNESS_NM * 1e-9
    p_membrane_um_s = (kp * dm_m2s / d_m) * 1e6  # m/s -> um/s
    if p_membrane_um_s <= 0.0:
        return 0.0

    # --- Aqueous leg: the unstirred boundary layer, in SERIES ---
    d_water_m2s = (_KB * temperature_K) / (6.0 * math.pi * _ETA_WATER_PA_S * r_m)
    p_aqueous_um_s = (d_water_m2s * 1e12) / _UWL_THICKNESS_UM  # um^2/s / um -> um/s

    p_um_s = 1.0 / (1.0 / p_membrane_um_s + 1.0 / p_aqueous_um_s)
    return max(0.0, p_um_s)


def bbb_penetrance(descriptors: Dict[str, float]) -> float:
    """Physicochemical estimate of blood-brain-barrier penetrance, in **[0, 1]**.

    Structure-only score built from the classic CNS multiparameter-optimization
    (CNS-MPO) style desirability of four descriptors, combined and squashed through
    a logistic so the output is a smooth fraction, not a hard pass/fail:

        favor  mw   < 450   (smaller crosses more easily)
        favor  1 < logp < 4 (moderately lipophilic; too polar or too greasy both hurt)
        favor  tpsa < 90    (low polar surface area crosses)
        favor  hbd  < 3     (few H-bond donors)

    Each descriptor maps to a smooth [0,1] desirability; the four desirabilities are
    combined with a **GEOMETRIC mean** and passed through a logistic centered so a
    "CNS-like" molecule scores high and a large polar molecule scores low.
    Antibodies -> ~0 (hard gate).

    WHY A GEOMETRIC MEAN AND NOT AN ARITHMETIC ONE (a veto, not an average)
    ----------------------------------------------------------------------
    Crossing the barrier is a CONJUNCTION of requirements: a molecule must be small
    enough AND lipophilic enough AND not too polar AND not too H-bonding. An
    ARITHMETIC mean cannot express that — it lets two favourable properties buy back
    two disqualifying ones. A geometric mean (equivalently, a sum of log
    desirabilities) is the standard aggregation for an AND of independent
    requirements: any single desirability near 0 drives the product near 0, i.e. a
    disqualifying property can VETO the molecule.

    This was not a stylistic choice — it is the fix for a measured failure. Against
    the 6-drug human-glioma tissue:plasma scorecard in
    ``results/glioma_regional_pk.md`` (Oh 2026 Neuro-Oncology PMC12962645; Eroglu
    2024 PMID 39422741; Blakeley 2008 PMID 18787762; Jackson 2018 PMID 29332604;
    Portnow 2009 PMID 19861433; Lee 2025 PMID 39820282), the arithmetic mean scored
    eribulin — a 730 Da macrocycle (MW desirability 0.009) with TPSA 146 (desirability
    0.009) but a favourable logP 3.44 (0.859) and 2 HBD (0.777) — at **0.334** against
    a MEASURED peritumoral brain:plasma ratio of **0.0051**, a 66x over-prediction and
    the single worst miss in the run. Two disqualifying properties were averaged away
    by two good ones. Under the geometric mean eribulin scores **0.035** (6.9x), and
    the worst-case error over the 6 drugs falls from 65.6x to 6.9x; all 6 drugs move
    inside the order-of-magnitude band (was 5/6). See the "post-fix" section of
    ``results/glioma_regional_pk.md`` for the full before/after table, including the
    drugs this makes *worse*.

    The change is deliberately minimal and adds NO tunable parameter: for a molecule
    whose four desirabilities are equal, geometric mean == arithmetic mean, so the
    score is unchanged. The two forms differ only for HETEROGENEOUS descriptor
    profiles — which is exactly the veto behaviour that was missing.

    HONEST LIMITATION — temozolomide:
        Temozolomide (mw ~194, logP ~ -2 to -1, TPSA ~106, hbd ~1) is clinically a
        first-line GBM drug BECAUSE it crosses the BBB well. But its logP is very low
        and its TPSA is above the 90 cutoff, so a PURE physicochemical rule set
        mis-scores it. We do NOT special-case TMZ here — the function stays honest —
        but among the small molecules it still yields a SANE ORDERING by
        physicochemical favorability, and it correctly puts antibodies and large polar
        drugs at ~0. Empirically the geometric mean moves TMZ from 0.558 (2.9x OVER
        the measured 0.195) to 0.125 (0.64x), i.e. closer to the measurement — but
        that is a *consequence* of the veto reshaping, not a TMZ-specific tuning, and
        TMZ's rank relative to cisplatin is unchanged.

    HONEST LIMITATION — the veto is descriptor-blind, and two drugs get worse:
        Methotrexate (TPSA 210, 5 HBD) and cefazolin (TPSA 156) also carry
        disqualifying descriptors, so they are vetoed too: methotrexate goes from
        0.65x to 0.13x of its measured pooled tumour value and cefazolin from 1.66x
        to 0.33x. Both of those measured means POOL contrast-enhancing tissue, where
        the barrier is disrupted; against the intact-barrier (non-enhancing) subset
        the vetoed scores are much closer (see ``reconstruct.barrier``, which supplies
        the additive regional leak that accounts for the difference). Still, a
        structure-only score cannot see active efflux or carrier-mediated uptake, and
        the residual error here is real.
    """
    if descriptors.get("is_antibody"):
        return 0.0

    mw = float(descriptors.get("mw", 300.0))
    logp = float(descriptors.get("logp", 0.0))
    tpsa = float(descriptors.get("tpsa", 75.0))
    hbd = float(descriptors.get("hbd", 1.0))

    # Hard gate: anything genuinely macromolecular can't cross passively.
    if mw >= 1000.0:
        return 0.0

    # --- Per-descriptor smooth desirabilities in [0,1] ---

    # MW: 1 below ~360, ramping down to 0 by ~700 (sigmoid centered at 450).
    d_mw = 1.0 / (1.0 + math.exp((mw - 450.0) / 60.0))

    # logP: bell centered on the CNS sweet spot (~2.5), penalizing both extremes.
    # Gaussian with width ~1.7 so 1<logp<4 stays high, and very polar (logp<0) or
    # very greasy (logp>5) falls off.
    d_logp = math.exp(-((logp - 2.5) ** 2) / (2.0 * 1.7 ** 2))

    # TPSA: 1 below ~60, ramping down to 0 past ~120 (sigmoid centered at 90).
    d_tpsa = 1.0 / (1.0 + math.exp((tpsa - 90.0) / 12.0))

    # HBD: 1 for 0-2 donors, falling off after 3 (sigmoid centered at 3).
    d_hbd = 1.0 / (1.0 + math.exp((hbd - 3.0) / 0.8))

    # --- Combine (GEOMETRIC mean desirability -> a veto) and squash ---
    # Product^(1/4), computed in log space for numerical stability. The floor keeps
    # a desirability that has underflowed to exactly 0 from producing log(0); it is
    # far below any score the logistic can distinguish, so it does not soften the
    # veto (a molecule at the floor still lands at ~0).
    desirabilities = (d_mw, d_logp, d_tpsa, d_hbd)
    geo_desir = math.exp(
        sum(math.log(max(d, _DESIRABILITY_FLOOR)) for d in desirabilities)
        / len(desirabilities)
    )
    # Logistic on the aggregated desirability: center at 0.5 so a balanced-favorable
    # molecule lands near the steep part of the curve; slope k=8 gives a smooth
    # spread across [0,1] without saturating too hard. Unchanged from the arithmetic
    # version, so a molecule with uniform desirabilities scores exactly as before.
    score = 1.0 / (1.0 + math.exp(-8.0 * (geo_desir - 0.5)))

    return max(0.0, min(1.0, score))


if __name__ == "__main__":
    # Contrasting molecules via inline descriptor dicts (no library import, to avoid
    # coupling). Values are reasonable literature-ballpark descriptors.
    molecules = {
        # Small, lipophilic, CNS-like (e.g. a lomustine/CCNU-ish nitrosourea profile).
        "small_lipophilic_CNS": dict(
            mw=233.0, logp=2.8, tpsa=62.0, hbd=1, hba=3, is_antibody=False,
        ),
        # Temozolomide-like: small but low logP and TPSA above the cutoff.
        "temozolomide_like": dict(
            mw=194.0, logp=-1.5, tpsa=106.0, hbd=1, hba=5, is_antibody=False,
        ),
        # Doxorubicin-like: large, polar, many donors -> poor passive/BBB.
        "doxorubicin_like": dict(
            mw=543.0, logp=1.3, tpsa=206.0, hbd=6, hba=12, is_antibody=False,
        ),
        # Antibody (bevacizumab-like): macromolecule, no passive crossing.
        "antibody": dict(
            mw=149000.0, logp=-1.0, tpsa=6000.0, hbd=200, hba=400, is_antibody=True,
        ),
    }

    print(f"{'molecule':<22} {'perm(um/s)':>14} {'bbb_penetrance':>16}")
    print("-" * 54)
    perm = {}
    bbb = {}
    for name, d in molecules.items():
        perm[name] = membrane_permeability(d)
        bbb[name] = bbb_penetrance(d)
        print(f"{name:<22} {perm[name]:>14.4e} {bbb[name]:>16.4f}")

    # --- Ordering assertions (the key validity checks) ---

    # 1) Antibody permeability and BBB are ~0, and << any small molecule.
    assert perm["antibody"] == 0.0, "antibody permeability must be ~0"
    assert bbb["antibody"] == 0.0, "antibody BBB penetrance must be ~0"
    for sm in ("small_lipophilic_CNS", "temozolomide_like", "doxorubicin_like"):
        assert bbb[sm] > bbb["antibody"], f"{sm} BBB must exceed antibody"
        assert perm[sm] > perm["antibody"], f"{sm} perm must exceed antibody"

    # 2) A lipophilic small CNS-like molecule scores HIGHER (BBB) than a large polar one.
    assert bbb["small_lipophilic_CNS"] > bbb["doxorubicin_like"], (
        "small lipophilic CNS molecule should out-score large polar molecule on BBB"
    )
    assert perm["small_lipophilic_CNS"] > perm["doxorubicin_like"], (
        "small lipophilic CNS molecule should out-permeate large polar molecule"
    )

    # 3) Sane ordering among the small molecules by physicochemical favorability:
    #    the CNS-like molecule is the most favorable of the small set.
    assert bbb["small_lipophilic_CNS"] > bbb["temozolomide_like"], (
        "physicochemical rules rank CNS-like > TMZ-like (documented TMZ under-prediction)"
    )

    # 4) The permeability CEILING: no molecule, however lipophilic, may exceed the
    #    unstirred-boundary-layer limit. Imatinib (Crippen logP 4.59) used to
    #    return 4.23e+07 um/s = 42 m/s.
    imatinib_like = dict(mw=493.6, logp=4.59, tpsa=86.3, hbd=2, hba=7,
                         is_antibody=False)
    p_imatinib = membrane_permeability(imatinib_like)
    assert p_imatinib < 100.0, f"imatinib-like permeability absurd: {p_imatinib:.4g} um/s"
    # A 10-log-unit molecule must not go any higher than a 4-log-unit one.
    absurd = dict(imatinib_like, logp=10.0)
    assert membrane_permeability(absurd) <= p_imatinib * 1.001, (
        "permeability must saturate, not keep climbing with logP"
    )
    for name, d in molecules.items():
        assert perm[name] < 100.0, f"{name} permeability {perm[name]:.4g} um/s exceeds ceiling"

    print("-" * 54)
    print(f"imatinib-like (logP 4.59): {p_imatinib:.3f} um/s "
          f"(was 4.23e+07 um/s = 42 m/s before the UWL ceiling)")
    print(f"logP 10 molecule:          {membrane_permeability(absurd):.3f} um/s (saturated)")
    print("All ordering assertions passed.")
