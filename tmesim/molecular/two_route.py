"""TWO-ROUTE tissue diffusivity: a PARACELLULAR route through the interstitium
between cells, plus a TRANSCELLULAR route THROUGH cells across their membranes.

Why this module exists
----------------------
`ecm_mesh.cell_exclusion` treats cells as IMPERMEABLE inclusions. That is the
right zeroth-order correction -- and `results/cellular_exclusion.md` shows it
moves two independent measured signatures a long way in the right direction --
but it is size- and chemistry-INDEPENDENT by construction: Maxwell's factor
`(1-phi)/(1+phi/2)` is identical for cisplatin and for docetaxel. It therefore
cannot produce ANY compound-to-compound spread, which is precisely what
`results/diffusivity_falsification.md` measured as the model's decisive failure:
across 71 compounds measured by one method in one cell line (Pruijn 2008,
doi:10.1071/CH08240) the measurement spans 29.7-fold and the engine 2.4-fold, at
R^2 = 0.028, and among four compounds of IDENTICAL molecular weight the
measurement still moves 3.7-fold while the engine is flat.

The authors of that dataset state the mechanism themselves: hydrophilic
compounds are confined to the PARACELLULAR route between cells, while lipophilic
ones additionally access the TRANSCELLULAR route through cell membranes. The
engine has no membrane, so it has no transcellular route at all -- and whether a
molecule can cross a bilayer is exactly what `logP` and hydrogen-bond donor count
predict. Both are computed by `descriptors.compute_descriptors` today and both
are discarded before the diffusivity calculation.

This module adds the missing route. It is ADDITIVE and DEFAULT-OFF: nothing in
`compile_drug.py` calls it, so `DrugParams.diffusivity` and every frozen
`DRUG_LIBRARY` snapshot are bit-identical with this file present.

THE PHYSICS
-----------
The tissue is a two-phase composite: an interstitial MATRIX of conductivity
`sigma_m` (the ECM-hindered interstitial diffusivity the engine already computes)
containing spherical CELL inclusions of radius `a`, volume fraction `phi`, and
conductivity `sigma_c = K * D_cyto`, separated from the matrix by a membrane of
finite permeability `P`.

Diffusion in a two-phase medium solves the same Laplace problem as conduction, so
the conduction result transfers unchanged (Hashin & Shtrikman 1962 state this
explicitly for diffusivity; see `ecm_mesh`). Maxwell's (1873) result for a
dispersion of spheres of FINITE conductivity is

    sigma_eff / sigma_m = (1 + 2*beta*phi) / (1 - beta*phi),
    beta = (sigma_c - sigma_m) / (sigma_c + 2*sigma_m).

A membrane of permeability `P` is an interfacial resistance in series with the
inclusion interior. Its effect is exactly a renormalisation of the inclusion
conductivity,

    sigma_c*  =  sigma_c / (1 + sigma_c/(a*P)),

i.e. the interior resistance `a/sigma_c` and the interfacial resistance `1/P` add.
Substituting `sigma_c*` for `sigma_c` in Maxwell's expression reproduces, term for
term, the composite conductivity of

  Hasselman DPH, Johnson LF. "Effective Thermal Conductivity of Composites with
  Interfacial Thermal Barrier Resistance." J Composite Materials 1987;21(6):
  508-515. doi:10.1177/002199838702100602

(their Eqn 8; their interfacial conductance `h` is this module's `P`, and their
`k_p/(a*h)` is this module's dimensionless membrane resistance `sigma_c/(a*P)`).
The algebraic identity between the two forms is asserted in the self-test below
and in `tests/test_two_route.py`, not asserted in prose.

THE THREE LIMITS, ALL CHECKED IN CODE
-------------------------------------
1. `P -> 0` (impermeable membrane): `sigma_c* -> 0`, `beta -> -1/2`, and the
   factor collapses to `(1-phi)/(1+phi/2)` -- EXACTLY
   `ecm_mesh.cell_exclusion(phi)["hindrance"]`, the law already shipped. So this
   module is a strict generalisation of what is already there and cannot
   contradict it.
2. `phi -> 0` (no cells): factor -> 1, the pure-interstitium answer.
3. `sigma_c* -> sigma_m` (cells indistinguishable from interstitium):
   `beta -> 0`, factor -> 1. A tissue with no contrast is transparent.

WHAT SETS THE CEILING, AND WHY IT IS NOT A FREE PARAMETER
---------------------------------------------------------
As `P -> infinity` the factor saturates at the Maxwell value for permeable
inclusions of conductivity `K*D_cyto`. So the MOST the transcellular route can
ever buy, at fixed phi, is set by the ratio of cytoplasmic to interstitial
diffusivity -- a measured quantity (Kao 1993, below), not a knob. The model
therefore predicts a BOUNDED spread, and that bound is a falsifiable statement.

Units: length um (nm where noted), time s, diffusivity um^2/s, permeability um/s.
No fitting, no dependence on `ecm_mesh` (import direction is ecm_mesh -> here).
"""
from __future__ import annotations

import math
from typing import Dict, Optional

# --------------------------------------------------------------------------- #
# Constants. Every one is either derived from the repo's own data, adopted
# unchanged from a module already in the package, or taken from a primary source
# quoted verbatim. The single exception is HBD_LOGP_DECREMENT, which is labelled
# a CONVENTION and has no source for its magnitude. See results/sweep/
# tworoute_prereg.json, which freezes all of these before any scoring, and
# results/two_route_transport.md.
# --------------------------------------------------------------------------- #

_KB = 1.380649e-23  # Boltzmann constant, J/K (CODATA exact)

#: Cell volume fraction. REUSED from results/cellular_exclusion.md 3.2 rather
#: than re-derived: Delesse's principle on the repo's own Xenium breast section
#: (167,780 segmented cells, sum(cell_area) = 3.2495e7 um^2, bounding-box tissue
#: area 4.118e7 um^2). An UPPER estimate -- Xenium's expansion segmentation tiles
#: the section. Literature tumour band implies 0.40-0.72.
DEFAULT_CELL_VOLUME_FRACTION = 0.7892

#: Cell radius, um. Equivalent-circle radius of the MEAN segmented cell area on
#: the same section: sqrt(193.675 um^2 / pi). Derived from the repo's own data.
#: It enters only through the dimensionless group sigma_c/(a*P), and the same
#: segmentation bias that inflates phi inflates `a`, which INCREASES the membrane
#: resistance and so UNDER-states the transcellular route. The bias is
#: conservative for this model.
DEFAULT_CELL_RADIUS_UM = 7.8517

#: Cytoplasmic diffusion coefficient of a small solute, relative to WATER.
#: Kao HP, Abney JR, Verkman AS. "Determinants of the translational mobility of a
#: small solute in cell cytoplasm." J Cell Biol 1993;120(1):175-184.
#: PMID 8416987, doi:10.1083/jcb.120.1.175. From the abstract, verbatim: "the
#: relative diffusion coefficients for BCECF translation and rotation in
#: cytoplasm were 0.27 +/- 0.01 (SEM, n = 24, 23 degrees C) and 0.78 +/- 0.03".
#: The paper decomposes it into fluid-phase viscosity (0.78), transient binding
#: (0.81) and collisions with cell solids (0.40), whose product 0.25 reproduces
#: the measurement -- i.e. it is a mechanistic number, not a black-box ratio.
CYTOPLASM_D_OVER_WATER = 0.27

#: Water viscosity at 310 K, Pa.s. Same value as permeability.py::_ETA_WATER_PA_S.
#: Used because Kao's ratio is against WATER; the engine's free_diffusivity uses
#: an interstitial-fluid viscosity of 1.2e-3 Pa.s, which is a claim about the
#: interstitium and must not be reused as a claim about cytoplasm.
ETA_WATER_PA_S = 6.9e-4

#: Bilayer-core viscosity, Pa.s. ADOPTED UNCHANGED from permeability.py::
#: _ETA_MEMBRANE_PA_S, which documents it as the TOP of the commonly reported
#: 0.1-1 Pa.s range and records the measured failure that made the bottom of the
#: range untenable. Adopting the codebase's already-published value -- fixed long
#: before the Pruijn dataset was retrieved -- is what keeps the membrane leg out
#: of reach of any tuning.
ETA_MEMBRANE_PA_S = 1.0

#: Hydrophobic bilayer thickness, nm. Unchanged from permeability.py.
BILAYER_THICKNESS_NM = 4.0

#: Upper clamp on log10 of the membrane/water partition coefficient. Unchanged
#: from permeability.py::_LOGP_KP_MAX, documented there as PHYSICAL: bilayer/water
#: partitioning plateaus above logP ~3 while octanol/water keeps climbing.
#: NOTE the lower clamp at -3 in permeability.py is deliberately NOT adopted --
#: it is a numerical guard, and flooring K at 1e-3 would force the transcellular
#: route OPEN for exactly the hydrophilic compounds whose exclusion from it is
#: the mechanism this module exists to represent.
LOGK_MEMBRANE_MAX = 3.0

#: Decrement applied to octanol logP per hydrogen-bond DONOR before it is used as
#: a bilayer/water partition coefficient, in log10 units.
#:
#: *** THIS IS A DECLARED CONVENTION, NOT A LITERATURE VALUE. ***
#:
#: The STRUCTURE is forced by two verified primary results:
#:   * Walter A, Gutknecht J. "Permeability of small nonelectrolytes through
#:     lipid bilayer membranes." J Membr Biol 1986;90(3):207-217, PMID 3735402,
#:     doi:10.1007/BF01870127 -- across 22 solutes, bilayer permeability
#:     correlates with HEXADECANE/water partition at r = 0.95 but with
#:     OCTANOL/water at only r = 0.75. Octanol is the wrong solvent for a bilayer
#:     because, unlike an alkane, wet octanol accepts hydrogen bonds; the
#:     correction from octanol to alkane is therefore dominated by the solute's
#:     H-bond DONOR capacity.
#:   * Conradi RA, Hilgers AR, Ho NF, Burton PS. Pharm Res 1991;8(12):1453-1460
#:     (PMID 1808606) and Pharm Res 1992;9(3):435-439 (PMID 1614980) -- over
#:     peptides spanning logP -2.2 to +2.8, "no correlation was found between the
#:     observed flux and the apparent lipophilicity", while "a strong correlation
#:     was found for the flux ... and the total number of hydrogen bonds the
#:     peptide could potentially make with water"; sequential N-methylation
#:     raised permeability with "only slight difference in the octanol-water
#:     partition coefficient".
#:
#: The MAGNITUDE has no primary source this work could verify. 1.0 log unit per
#: donor is a round number fixed before the model was ever scored (see
#: results/sweep/tworoute_prereg.json). results/two_route_transport.md sweeps it
#: and discloses the value that would maximise agreement with the measurement, so
#: a reader can confirm the adopted value is not that one.
HBD_LOGP_DECREMENT = 1.0


# --------------------------------------------------------------------------- #
# Leg 1: the membrane. logP + HBD -> a passive bilayer permeability.
# --------------------------------------------------------------------------- #
def membrane_permeability_um_s(
    logp: float,
    hbd: float,
    solute_radius_nm: float,
    temperature_K: float = 310.0,
    hbd_decrement: float = HBD_LOGP_DECREMENT,
    logk_max: float = LOGK_MEMBRANE_MAX,
) -> Dict[str, float]:
    """Passive permeability of ONE cell membrane, um/s, from logP and HBD only.

    Solubility-diffusion (Overton/Meyer; Walter & Gutknecht 1986 confirm this
    model form against 22 measured bilayer permeabilities):

        log10 K_mem = min(logP - hbd_decrement * n_HBD,  logk_max)
        D_mem       = kB*T / (6*pi*eta_mem*r)            [Stokes-Einstein in the
                                                          bilayer core]
        P           = K_mem * D_mem / d

    Why this is NOT `permeability.membrane_permeability`
    ---------------------------------------------------
    That function is the right model for a solute arriving at a cell from a
    STIRRED BATH, and it therefore places a 30 um unstirred aqueous boundary
    layer in series with the bilayer. Inside tissue there is no such layer: the
    aqueous resistance around a cell IS the interstitium, and the interstitium is
    already the matrix phase of the effective-medium problem this permeability
    feeds. Reusing that function here would double-count the aqueous leg and, far
    worse, would cap every lipophilic compound at the same boundary-layer-limited
    ~22 um/s -- destroying the very compound-to-compound discrimination this
    module exists to produce. That function's own docstring records the resulting
    compression ("the top of the ranking is compressed ... within 15% of each
    other") as a known defect.

    Every CONSTANT is nevertheless taken unchanged from that module (bilayer
    viscosity, thickness, upper partition clamp), so no new tunable number is
    introduced; only the series topology differs, and it differs for a stated
    physical reason.

    Returns a dict so the caller can see the partition and the diffusion leg
    separately: {log_k_membrane, k_membrane, d_membrane_um2s, permeability_um_s,
    clamped}.
    """
    if solute_radius_nm <= 0:
        raise ValueError(f"solute_radius_nm must be positive, got {solute_radius_nm!r}")
    if temperature_K <= 0:
        raise ValueError(f"temperature_K must be positive, got {temperature_K!r}")
    if hbd < 0:
        raise ValueError(f"hbd must be >= 0, got {hbd!r}")

    raw = float(logp) - float(hbd_decrement) * float(hbd)
    log_k = min(raw, float(logk_max))
    k_mem = 10.0 ** log_k

    r_m = solute_radius_nm * 1e-9
    d_mem_m2s = (_KB * temperature_K) / (6.0 * math.pi * ETA_MEMBRANE_PA_S * r_m)
    d_mem_um2s = d_mem_m2s * 1e12                      # m^2/s -> um^2/s
    d_um = BILAYER_THICKNESS_NM * 1e-3                 # nm -> um

    return {
        "log_k_membrane": log_k,
        "log_k_membrane_unclamped": raw,
        "k_membrane": k_mem,
        "d_membrane_um2s": d_mem_um2s,
        "permeability_um_s": k_mem * d_mem_um2s / d_um,
        "clamped": raw > logk_max,
    }


def cytoplasm_diffusivity_um2s(
    solute_radius_nm: float,
    temperature_K: float = 310.0,
    d_over_water: float = CYTOPLASM_D_OVER_WATER,
) -> float:
    """Cytoplasmic diffusivity, um^2/s = ``d_over_water`` x Stokes-Einstein in water.

    The reference fluid is WATER (eta = 6.9e-4 Pa.s at 310 K), not the engine's
    1.2e-3 Pa.s "interstitial fluid", because Kao 1993's measured ratio is against
    water. Using the interstitial viscosity here would silently apply an
    interstitial claim to the cytoplasm.
    """
    if solute_radius_nm <= 0:
        raise ValueError(f"solute_radius_nm must be positive, got {solute_radius_nm!r}")
    r_m = solute_radius_nm * 1e-9
    d_water_um2s = (_KB * temperature_K) / (
        6.0 * math.pi * ETA_WATER_PA_S * r_m
    ) * 1e12
    return float(d_over_water) * d_water_um2s


# --------------------------------------------------------------------------- #
# Leg 2: the effective-medium law with a permeable inclusion phase.
# --------------------------------------------------------------------------- #
def two_route_factor(
    cell_volume_fraction: float,
    interstitial_D_um2s: float,
    cytoplasm_D_um2s: float,
    membrane_permeability_um_s: float,
    cell_radius_um: float = DEFAULT_CELL_RADIUS_UM,
    cell_partition: float = 1.0,
) -> Dict[str, float]:
    """Maxwell effective medium with a MEMBRANE-limited permeable inclusion phase.

    Returns the factor by which the INTERSTITIAL diffusivity must be multiplied
    to give the tissue-scale effective diffusivity, together with every
    intermediate so the algebra can be audited.

    ``hindrance`` is on the same FLUX/conductivity basis as
    ``ecm_mesh.cell_exclusion(phi)["hindrance"]``, which it reduces to exactly at
    zero membrane permeability -- see ``impermeable_hindrance`` and
    ``transcellular_gain`` in the returned dict.

    Parameters
    ----------
    cell_volume_fraction : float
        phi, in [0, 1).
    interstitial_D_um2s : float
        sigma_m -- the matrix conductivity. This is the diffusivity of the
        interstitial space AFTER the ECM legs the engine already applies
        (Brownian micro-tortuosity x Ogston-Renkin mesh gate), i.e. exactly
        today's ``DrugParams.diffusivity``. Scale separation licenses the
        composition: the fibre mesh has a ~7 nm pore and cells are ~16 um across.
    cytoplasm_D_um2s : float
        D_cyto -- diffusivity INSIDE a cell. Not ECM-hindered: cytoplasmic
        crowding is a different obstacle population and is already inside Kao's
        measured 0.27.
    membrane_permeability_um_s : float
        P. Zero gives the impermeable-cell law exactly.
    cell_radius_um : float
        a. Enters only through the dimensionless resistance ratio.
    cell_partition : float
        K, the cytoplasm/interstitium equilibrium partition of the solute. The
        inclusion conductivity is K*D_cyto because the driving potential across a
        partitioning interface is the concentration normalised by K.
    """
    phi = float(cell_volume_fraction)
    if not (0.0 <= phi < 1.0):
        raise ValueError(f"cell_volume_fraction must be in [0,1), got {phi!r}")
    if interstitial_D_um2s <= 0:
        raise ValueError(
            f"interstitial_D_um2s must be positive, got {interstitial_D_um2s!r}"
        )
    if cytoplasm_D_um2s < 0:
        raise ValueError(f"cytoplasm_D_um2s must be >= 0, got {cytoplasm_D_um2s!r}")
    if membrane_permeability_um_s < 0:
        raise ValueError(
            f"membrane_permeability_um_s must be >= 0, got "
            f"{membrane_permeability_um_s!r}"
        )
    if cell_radius_um <= 0:
        raise ValueError(f"cell_radius_um must be positive, got {cell_radius_um!r}")

    sigma_m = float(interstitial_D_um2s)
    sigma_c = float(cell_partition) * float(cytoplasm_D_um2s)
    p = float(membrane_permeability_um_s)

    # Dimensionless membrane resistance: interior resistance a/sigma_c against
    # interfacial resistance 1/P. Large -> membrane-limited -> impermeable cells.
    if sigma_c <= 0.0:
        biot = float("inf")
        sigma_c_eff = 0.0
    elif p <= 0.0:
        biot = float("inf")
        sigma_c_eff = 0.0
    else:
        biot = sigma_c / (cell_radius_um * p)
        sigma_c_eff = sigma_c / (1.0 + biot)

    beta = (sigma_c_eff - sigma_m) / (sigma_c_eff + 2.0 * sigma_m)
    hind = (1.0 + 2.0 * beta * phi) / (1.0 - beta * phi)

    impermeable = (1.0 - phi) / (1.0 + phi / 2.0)
    return {
        "cell_volume_fraction": phi,
        "sigma_matrix_um2s": sigma_m,
        "sigma_cell_um2s": sigma_c,
        "sigma_cell_effective_um2s": sigma_c_eff,
        "membrane_resistance_ratio": biot,
        "beta": beta,
        "hindrance": hind,
        "fold_reduction": (1.0 / hind) if hind > 0 else float("inf"),
        "impermeable_hindrance": impermeable,
        "transcellular_gain": hind / impermeable if impermeable > 0 else float("nan"),
        "effective_D_um2s": sigma_m * hind,
        "paracellular_D_um2s": sigma_m * impermeable,
    }


def two_route_ceiling(
    cell_volume_fraction: float,
    interstitial_D_um2s: float,
    cytoplasm_D_um2s: float,
    cell_partition: float = 1.0,
) -> Dict[str, float]:
    """The P -> infinity limit: the most the transcellular route can ever buy.

    A registered, falsifiable prediction of the model: no compound, however
    lipophilic, can exceed this. Depends only on phi and the cytoplasm/interstitium
    conductivity contrast -- both measured, neither tunable.
    """
    return two_route_factor(
        cell_volume_fraction,
        interstitial_D_um2s,
        cytoplasm_D_um2s,
        membrane_permeability_um_s=float("inf"),
        cell_radius_um=DEFAULT_CELL_RADIUS_UM,
        cell_partition=cell_partition,
    )


def hasselman_johnson_hindrance(
    cell_volume_fraction: float,
    k_matrix: float,
    k_inclusion: float,
    inclusion_radius: float,
    interfacial_conductance: float,
) -> float:
    """Hasselman & Johnson (1987) Eqn 8, written out literally, for cross-checking.

    Their published form, in their own symbols (k_m matrix, k_p particle, a
    radius, h interfacial conductance, phi particle volume fraction):

        k_eff/k_m = [ 2*(k_p/k_m - k_p/(a*h) - 1)*phi + k_p/k_m + 2*k_p/(a*h) + 2 ]
                  / [   (1 - k_p/k_m + k_p/(a*h))*phi + k_p/k_m + 2*k_p/(a*h) + 2 ]

    `two_route_factor` computes the same quantity by the equivalent and more
    transparent route of substituting a series-resistance inclusion conductivity
    into Maxwell's expression. The two agree to machine precision; that identity
    is the arithmetic check that the substitution is legitimate, and it is
    asserted in the self-test and in tests/test_two_route.py rather than claimed.
    """
    phi = float(cell_volume_fraction)
    kp_km = k_inclusion / k_matrix
    x = k_inclusion / (inclusion_radius * interfacial_conductance)
    num = 2.0 * (kp_km - x - 1.0) * phi + kp_km + 2.0 * x + 2.0
    den = (1.0 - kp_km + x) * phi + kp_km + 2.0 * x + 2.0
    return num / den


def two_route_diffusivity(
    interstitial_D_um2s: float,
    solute_radius_nm: float,
    logp: float,
    hbd: float,
    cell_volume_fraction: float = DEFAULT_CELL_VOLUME_FRACTION,
    cell_radius_um: float = DEFAULT_CELL_RADIUS_UM,
    temperature_K: float = 310.0,
    cell_partition: float = 1.0,
    hbd_decrement: float = HBD_LOGP_DECREMENT,
    d_over_water: float = CYTOPLASM_D_OVER_WATER,
) -> Dict[str, float]:
    """End-to-end: engine interstitial D + (logP, HBD, radius) -> tissue D_eff.

    This is the one call an experiment or a wiring point needs. It is NOT called
    from `compile_drug.py`; the two-route model is default-off.
    """
    mem = membrane_permeability_um_s(
        logp, hbd, solute_radius_nm, temperature_K=temperature_K,
        hbd_decrement=hbd_decrement,
    )
    d_cyto = cytoplasm_diffusivity_um2s(
        solute_radius_nm, temperature_K=temperature_K, d_over_water=d_over_water
    )
    res = two_route_factor(
        cell_volume_fraction,
        interstitial_D_um2s,
        d_cyto,
        mem["permeability_um_s"],
        cell_radius_um=cell_radius_um,
        cell_partition=cell_partition,
    )
    res.update({
        "logp": float(logp),
        "hbd": float(hbd),
        "solute_radius_nm": float(solute_radius_nm),
        "membrane_permeability_um_s": mem["permeability_um_s"],
        "log_k_membrane": mem["log_k_membrane"],
        "membrane_clamped": mem["clamped"],
        "cytoplasm_D_um2s": d_cyto,
    })
    return res


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    import time as _time

    t0 = _time.time()
    print("=" * 78)
    print("two_route self-test: paracellular + transcellular effective diffusivity")
    print("=" * 78)

    PHI = DEFAULT_CELL_VOLUME_FRACTION
    SIGMA_M = 265.0   # um^2/s, an ECM-hindered interstitial D for a small molecule
    D_CYTO = 200.0    # um^2/s

    # ---- [1] LIMIT: P -> 0 reproduces ecm_mesh.cell_exclusion EXACTLY ---------
    impermeable = (1.0 - PHI) / (1.0 + PHI / 2.0)
    r0 = two_route_factor(PHI, SIGMA_M, D_CYTO, 0.0)
    print(f"\n[1] P -> 0 limit at phi = {PHI}")
    print(f"    two_route hindrance        = {r0['hindrance']:.17g}")
    print(f"    (1-phi)/(1+phi/2)          = {impermeable:.17g}")
    assert abs(r0["hindrance"] - impermeable) < 1e-15, r0["hindrance"]
    tiny = two_route_factor(PHI, SIGMA_M, D_CYTO, 1e-12)
    assert abs(tiny["hindrance"] - impermeable) < 1e-9, tiny["hindrance"]
    print("    -> reduces to the impermeable-cell law already shipped. OK")

    # ---- [2] LIMIT: phi -> 0, and zero-contrast --------------------------------
    assert abs(two_route_factor(0.0, SIGMA_M, D_CYTO, 1.0)["hindrance"] - 1.0) < 1e-15
    same = two_route_factor(PHI, SIGMA_M, SIGMA_M, float("inf"))
    print(f"[2] phi -> 0 gives 1.0; zero-contrast tissue gives "
          f"{same['hindrance']:.17g}")
    assert abs(same["hindrance"] - 1.0) < 1e-15

    # ---- [3] IDENTITY with Hasselman & Johnson (1987) Eqn 8 --------------------
    print("[3] Maxwell-with-series-resistance vs Hasselman & Johnson Eqn 8:")
    worst = 0.0
    for p in (1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e4):
        mine = two_route_factor(PHI, SIGMA_M, D_CYTO, p)["hindrance"]
        theirs = hasselman_johnson_hindrance(PHI, SIGMA_M, D_CYTO,
                                             DEFAULT_CELL_RADIUS_UM, p)
        worst = max(worst, abs(mine - theirs))
        print(f"    P = {p:>9.4g} um/s   mine = {mine:.12f}   HJ87 = {theirs:.12f}")
    assert worst < 1e-14, f"HJ87 identity fails at {worst:.3g}"
    print(f"    max |difference| = {worst:.3g}  -> identical. OK")

    # ---- [4] MONOTONE and BOUNDED ---------------------------------------------
    ps = [0.0, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1e3, 1e6]
    hs = [two_route_factor(PHI, SIGMA_M, D_CYTO, p)["hindrance"] for p in ps]
    assert all(hs[i] <= hs[i + 1] + 1e-15 for i in range(len(hs) - 1)), hs
    ceiling = two_route_ceiling(PHI, SIGMA_M, D_CYTO)["hindrance"]
    assert all(h <= ceiling + 1e-12 for h in hs)
    print(f"[4] hindrance is monotone in P, from {hs[0]:.4f} to a ceiling of "
          f"{ceiling:.4f} ({ceiling/impermeable:.2f}x the impermeable value)")

    # ---- [5] the membrane leg discriminates on logP and HBD -------------------
    print("\n[5] membrane permeability from (logP, HBD) at r = 0.45 nm:")
    print(f"    {'logP':>6} {'HBD':>4} {'log K_mem':>10} {'P (um/s)':>12} "
          f"{'hindrance':>10} {'gain':>6}")
    rows = [(-2.01, 4), (-2.62, 2), (-0.38, 2), (0.0, 1), (1.5, 1), (2.96, 1)]
    prev_h = -1.0
    for lp, hd in rows:
        r = two_route_diffusivity(SIGMA_M, 0.45, lp, hd, cell_volume_fraction=PHI)
        print(f"    {lp:>6.2f} {hd:>4} {r['log_k_membrane']:>10.2f} "
              f"{r['membrane_permeability_um_s']:>12.4g} {r['hindrance']:>10.4f} "
              f"{r['transcellular_gain']:>6.2f}")
    spread = max(
        two_route_diffusivity(SIGMA_M, 0.45, lp, hd)["hindrance"] for lp, hd in rows
    ) / min(
        two_route_diffusivity(SIGMA_M, 0.45, lp, hd)["hindrance"] for lp, hd in rows
    )
    print(f"    spread across these six at FIXED size = {spread:.2f}x "
          f"(the impermeable law gives exactly 1.00x)")
    assert spread > 1.5, f"two-route model produces no spread at fixed size: {spread}"

    # more donors must never help; more logP must never hurt
    a = two_route_diffusivity(SIGMA_M, 0.45, 0.0, 1)["hindrance"]
    b = two_route_diffusivity(SIGMA_M, 0.45, 0.0, 3)["hindrance"]
    c = two_route_diffusivity(SIGMA_M, 0.45, 2.0, 1)["hindrance"]
    assert b < a <= c, (a, b, c)
    print("    monotone: more HBD -> slower, higher logP -> faster. OK")

    print(f"\nAll asserts passed. Runtime {_time.time() - t0:.2f} s")
