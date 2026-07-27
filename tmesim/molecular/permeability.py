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

    print("-" * 54)
    print("All ordering assertions passed.")
