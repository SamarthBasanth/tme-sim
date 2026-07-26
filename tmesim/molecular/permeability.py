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
"""
from __future__ import annotations

import math
from typing import Dict

# Physical constants (SI unless noted)
_KB = 1.380649e-23  # Boltzmann constant, J/K
_ETA_MEMBRANE_PA_S = 0.10  # effective viscosity inside a lipid bilayer core, Pa*s.
# The bilayer interior is far more viscous than water (~1e-3 Pa*s); membrane-core
# microviscosities of ~0.1-1 Pa*s are commonly reported. We use 0.1 Pa*s so that
# Dm lands in a physically sane ~1e-11 m^2/s (1e1 um^2/s) range for a small solute.
_BILAYER_THICKNESS_NM = 4.0  # d, hydrophobic bilayer thickness, nm (~4 nm typical)


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

    Solubility-diffusion model (Overton / Meyer picture):

        P = Kp * Dm / d

    where
      * Kp  : membrane/water partition coefficient (dimensionless). Driven by
              lipophilicity. From logP:  Kp ~= 10**logP, clamped to a sane window
              so that a very lipophilic molecule does not blow the estimate up and a
              very polar one still gets a tiny finite value. We clamp logP into
              [-3, 5] before exponentiating (Kp in [1e-3, 1e5]).
      * Dm  : lateral/normal diffusion coefficient of the solute *inside* the
              bilayer core, m^2/s. Stokes-Einstein with a membrane-core viscosity:
                  Dm = kB*T / (6*pi*eta_membrane*r)
              Smaller solute -> larger Dm (faster). r comes from the hydrodynamic
              radius or an MW->radius estimate.
      * d   : bilayer hydrophobic thickness, ~4 nm.

    Units bookkeeping:
        Dm [m^2/s] / d [m]  ->  [m/s];  * Kp (dimensionless)  ->  [m/s].
        Convert m/s -> um/s by * 1e6.

    Antibodies / macromolecules (`is_antibody`) do not cross by passive diffusion:
    return ~0.

    Returns a non-negative float in um/s.
    """
    if descriptors.get("is_antibody"):
        return 0.0

    # --- Kp from logP, clamped ---
    logp = float(descriptors.get("logp", 0.0))
    logp_clamped = max(-3.0, min(5.0, logp))
    kp = 10.0 ** logp_clamped  # dimensionless, in [1e-3, 1e5]

    # --- Dm from size (Stokes-Einstein in the membrane core) ---
    r_nm = _radius_nm_from_descriptors(descriptors)
    r_m = max(r_nm, 0.1) * 1e-9  # guard against tiny/zero radius
    dm_m2s = (_KB * temperature_K) / (6.0 * math.pi * _ETA_MEMBRANE_PA_S * r_m)

    # --- Assemble P = Kp * Dm / d ---
    d_m = _BILAYER_THICKNESS_NM * 1e-9
    p_m_s = kp * dm_m2s / d_m  # m/s
    p_um_s = p_m_s * 1e6  # um/s

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

    Each descriptor maps to a smooth [0,1] desirability; their mean is rescaled and
    passed through a logistic centered so a "CNS-like" molecule scores high and a
    large polar molecule scores low. Antibodies -> ~0 (hard gate).

    HONEST LIMITATION — temozolomide:
        Temozolomide (mw ~194, logP ~ -2 to -1, TPSA ~106, hbd ~1) is clinically a
        first-line GBM drug BECAUSE it crosses the BBB well. But its logP is very low
        and its TPSA is above the 90 cutoff, so a PURE physicochemical rule set
        UNDER-predicts its penetrance. This is a real, well-documented limitation of
        structure-only BBB estimation (TMZ's small size and specific transport/
        stability behavior aren't captured by 2D descriptors). We do NOT special-case
        TMZ here — the function stays honest — but among the small molecules it still
        yields a SANE ORDERING by physicochemical favorability, and it correctly puts
        antibodies and large polar drugs at ~0.
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

    # --- Combine (mean desirability) and squash ---
    mean_desir = (d_mw + d_logp + d_tpsa + d_hbd) / 4.0
    # Logistic on the mean desirability: center at 0.5 so a balanced-favorable
    # molecule lands near the steep part of the curve; slope k=8 gives a smooth
    # spread across [0,1] without saturating too hard.
    score = 1.0 / (1.0 + math.exp(-8.0 * (mean_desir - 0.5)))

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

    print("-" * 54)
    print("All ordering assertions passed.")
