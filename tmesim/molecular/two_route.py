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


def efflux_reduced_partition(cell_partition: float, efflux_ratio: float) -> float:
    """Cytoplasm/interstitium partition ``K`` in the presence of active efflux.

    ``K_eff = K / ER``. This is NOT an assumption; it is the exact steady-state
    solution of a cell with efflux, and it is DERIVED below rather than adopted.

    THE DERIVATION
    --------------
    Take one cell as a slab crossed by the transport direction: an entry
    membrane, an interior of conductivity ``sigma_c`` over a path length ``L``,
    and an exit membrane. Each membrane has passive permeability ``P`` and an
    ADDITIONAL outward active clearance ``E`` (the pump). With outside
    concentrations ``C1`` and ``C2`` and just-inside concentrations ``Ca``, ``Cb``,
    the three steady-state flux statements are

        J = P*(C1 - Ca) - E*Ca          [entry face: passive in, pump out]
        J = (sigma_c/L)*(Ca - Cb)       [the interior]
        J = P*(Cb - C2) + E*Cb          [exit face]

    Eliminating ``Ca`` and ``Cb`` gives, in closed form,

        (C1 - C2)/J  =  2/P  +  (L/sigma_c) * ER,        ER = 1 + E/P

    i.e. **efflux multiplies the INTERIOR resistance by the efflux ratio and
    leaves the membrane resistance untouched.** In this module's
    series-resistance parameterisation

        1/sigma_c*  =  1/sigma_c + 1/(a*P)      ->      ER/sigma_c + 1/(a*P)

    which is exactly ``sigma_c -> sigma_c/ER``, and since ``sigma_c = K*D_cyto``
    it is exactly ``K -> K/ER`` with ``P`` unchanged. So efflux enters this model
    at a parameter that ALREADY EXISTS and introduces **no new degree of
    freedom**. The identity is asserted numerically in the self-test below and in
    ``tests/test_two_route.py``, not claimed in prose.

    WHY THAT IS THE SAME OBJECT ``efflux.FORM_KPUU`` ALREADY MODELS
    ---------------------------------------------------------------
    ``efflux.py``'s default form is ``f = Kp,uu = P_pass/(P_pass + CL_eff) =
    1/ER`` -- an UNBOUND PARTITION RATIO across one membrane at steady state.
    ``K`` in this module is also an unbound partition ratio across one membrane
    (cytoplasm over interstitium). They are the same quantity, so composing them
    is multiplication of a partition by a partition, not a coincidence of form.
    The derivation above is compartment-agnostic and does not depend on which
    membrane is meant.

    THE CONSEQUENCE THAT MATTERS, AND IT IS FALSIFIABLE
    ---------------------------------------------------
    The interior resistance is the term efflux multiplies. When the membrane is
    SHUT (``1/(a*P) >> 1/sigma_c``, i.e. a large ``membrane_resistance_ratio``)
    the interior contributes nothing to the series and efflux therefore changes
    the tissue diffusivity by **exactly nothing**. Doxorubicin sits there
    (``membrane_resistance_ratio`` ~ 2e5), so efflux cannot move it however large
    its efflux ratio is. Under the PREVIOUS transport model, which had no
    membrane and no cell inclusions at all, ``K`` never entered the calculation,
    so this term was identically zero **for every drug**. That is the precise
    sense in which efflux has nothing to act on until a transcellular route
    exists -- and it is a prediction of the composition, not a rationalisation
    of one.

    Parameters
    ----------
    cell_partition : float
        ``K``, the passive cytoplasm/interstitium equilibrium partition.
    efflux_ratio : float
        ``ER = Papp(B->A)/Papp(A->B) = 1 + CL_efflux/P_passive >= 1``. A value
        below 1 means net absorptive transport, which this module does not
        represent; it is clamped to 1 (never a boost), matching
        ``efflux._apply_form``.
    """
    if cell_partition < 0:
        raise ValueError(f"cell_partition must be >= 0, got {cell_partition!r}")
    er = float(efflux_ratio)
    if er <= 0:
        raise ValueError(f"efflux_ratio must be positive, got {efflux_ratio!r}")
    return float(cell_partition) / max(1.0, er)


# --------------------------------------------------------------------------- #
# Leg 3: INTRACELLULAR RETARDATION. Reversible binding inside the cell.
#
# See results/intracellular_retardation.md. Everything below is ADDITIVE and
# DEFAULT-OFF: `bound_to_free_ratio=None` is bitwise identical to omitting it.
# --------------------------------------------------------------------------- #

#: Intracellular "microsomal-protein-equivalent" binding capacity, mg/mL of CELL
#: volume, in Kilford's hepatocyte extension of the Hallifax & Houston microsomal
#: binding equation. See :func:`binding_capacity_from_logp` for the derivation and
#: for why the incubation volume ratio cancels exactly.
#:
#: Kilford PJ, Gertz M, Houston JB, Galetin A. "Hepatocellular binding of drugs:
#: correction for unbound fraction in hepatocyte incubations using microsomal
#: binding or drug lipophilicity data." Drug Metab Dispos 2008;36(7):1194-1197.
#: PMID 18411401, doi:10.1124/dmd.108.020834.
KILFORD_MG_EQUIV_PER_ML_CELL = 125.0

#: Coefficients of the Hallifax & Houston (2006) log-quadratic binding function
#: `10**(c2*logP**2 + c1*logP + c0)`, which is the `fu,mic` term Kilford substitutes.
#:
#: Hallifax D, Houston JB. "Binding of drugs to hepatic microsomes: comment and
#: assessment of current prediction methodology with recommendation for
#: improvement." Drug Metab Dispos 2006;34(4):724-726. PMID 16552024,
#: doi:10.1124/dmd.105.007658.
HALLIFAX_C2 = 0.072
HALLIFAX_C1 = 0.067
HALLIFAX_C0 = -1.126


def binding_capacity_from_logp(
    logp: float,
    mg_equiv_per_ml_cell: float = KILFORD_MG_EQUIV_PER_ML_CELL,
) -> float:
    """Intracellular bound-to-free ratio ``K_b`` from logP alone. STRUCTURAL RULE.

    ``K_b`` is the ratio of BOUND to FREE solute inside a cell at equilibrium. It
    is the quantity that sets the retardation factor (:func:`retardation_factor`),
    and it is the one input of the retardation term.

    WHY THIS RULE, AND WHAT IT IS DERIVED FROM
    ------------------------------------------
    No per-compound intracellular bound fraction exists for the 71-compound
    benzotriazine di-N-oxide set this model is scored against, so a per-compound
    ``K_b`` cannot be sourced. What CAN be sourced is a STRUCTURAL rule from a
    descriptor the engine already computes. Kilford 2008 predicts the fraction
    unbound in a hepatocyte incubation from lipophilicity alone:

        fu_hep = 1 / (1 + 125 * VR * X),
        X      = 10**(0.072*logP**2 + 0.067*logP - 1.126)     [Hallifax & Houston]

    where ``VR`` is the volume ratio of CELLS to incubation. The bound-to-free
    ratio in the incubation is therefore ``125 * VR * X``. That bound material all
    sits INSIDE the cells, which occupy a fraction ``VR`` of the incubation, so per
    unit CELL volume the bound-to-free ratio is

        K_b = (125 * VR * X) / VR = 125 * X.

    **The volume ratio cancels exactly.** That is what makes this usable here: the
    result is a property of a millilitre of cells, not of Kilford's assay geometry,
    so it transfers to a cell sitting in a tissue. ``125`` carries units of mg
    microsomal-protein-equivalent per mL of cell volume, against a whole-cell
    protein content of order 200-250 mg/mL -- i.e. it is a plausible magnitude and
    not a free parameter.

    WHAT THIS RULE IS NOT, STATED HERE RATHER THAN BURIED
    -----------------------------------------------------
    * It is an **empirical QSAR fitted to drug-like compounds binding to liver
      microsomes**, not a measurement on any compound in this model's compound set,
      and not a measurement on HT-29.
    * Its intercept is large: at ``logP = 0`` it returns ``K_b = 125*10**-1.126 =
      9.35``, i.e. it asserts that even a compound with no lipophilicity at all is
      90 % bound inside a cell. Applied to a genuinely non-binding flux marker such
      as urea or mannitol -- both of which are in the Pruijn set, and one of which
      that paper uses as its INTERNAL STANDARD -- that is almost certainly out of
      domain. This is registered as an expected failure mode BEFORE scoring, in
      results/sweep/retard_prereg.json, not discovered afterwards.
    * Kilford's own abstract records that prediction accuracy is "high for
      relatively hydrophilic drugs (logP/D <= 2.5)" and that "less accurate
      prediction of the fu(hep) is observed for lipophilic drugs (logP > 3)".
    * As published the rule uses **logP for BASES and logD7.4 for ACIDS AND
      NEUTRALS** -- not, as is often misremembered, logP for neutrals. This
      function is given whatever lipophilicity its caller holds and does not know
      the compound's acid/base class, so it cannot apply that distinction. For a
      genuinely neutral compound logP and logD7.4 coincide and the point is moot;
      for the weakly basic members of the anchor set it is a real, uncorrected
      approximation.

    PROVENANCE OF THE COEFFICIENTS, STATED EXACTLY
    ----------------------------------------------
    **Both primary papers are paywalled and neither was read.** Kilford 2008 and
    Hallifax & Houston 2006 are both closed at ASPET (HTTP 403), neither has a
    PMC record, and Unpaywall reports ``is_oa: false`` for both DOIs. Kilford's
    ABSTRACT is verified and confirms the substitution architecture in the
    authors' own words, but carries no coefficients; Hallifax & Houston is a
    Letter and PubMed records "[Abstract not available]".

    The equation and all four constants were therefore read VERBATIM from two
    INDEPENDENT open-access secondary sources, which agree exactly:
      * Eur J Pharm Sci 2018, PMID 30098391, PMC6365006,
        doi:10.1016/j.ejps.2018.08.008 -- gives both equations, defines ``C`` as
        "the microsomal protein concentration in mg/mL", and defines ``VR`` as
        "the ratio of cell volume to incubation volume ... 0.005 at a cell
        concentration of 10^6 cells/mL".
      * Yim DS, Bae SH, Choi S. Transl Clin Pharmacol 2021;29(2), PMID 34235120,
        PMC8255549, doi:10.12793/tcp.2021.29.e12 -- same equation, same
        constants, same ``VR`` definition and value.
    This is a weaker provenance than a primary read and it is labelled as such
    rather than dressed up.

    Finally: the reading of ``125`` as mg protein-equivalent per mL of CELL volume
    follows from the algebra of the substitution (``125*VR`` occupies exactly the
    slot ``C`` in mg/mL occupies, and ``VR`` is dimensionless), and nothing
    retrieved contradicts it -- but **no source states that the authors intended
    it as a physical protein concentration.** Kilford's own abstract calls it "A
    new nonlinear empirical equation". It is an empirical fit that is
    dimensionally interpretable this way, not a measured intracellular protein
    content, and that distinction is recorded rather than glossed.

    The function is pure and total: it never returns a negative ratio.
    """
    lp = float(logp)
    x = 10.0 ** (HALLIFAX_C2 * lp * lp + HALLIFAX_C1 * lp + HALLIFAX_C0)
    return float(mg_equiv_per_ml_cell) * x


def retardation_factor(
    bound_to_free_ratio: float,
    cell_volume_fraction: float = DEFAULT_CELL_VOLUME_FRACTION,
    cell_partition: float = 1.0,
) -> float:
    """Retardation factor ``R`` by which reversible intracellular binding divides
    the tissue-scale APPARENT diffusivity. ``R = 1 + phi*K*K_b``.

    THE PHYSICS, AND WHY IT IS NOT WHAT IT LOOKS LIKE
    -------------------------------------------------
    The tempting form is "binding slows the molecule down, so put
    ``D_cyto -> D_cyto/(1+K_b)`` into the inclusion conductivity". **That is wrong,
    and the error matters.** Reversible binding at local equilibrium does not
    change a STEADY-STATE conductance at all: at steady state the bound pool is
    time-invariant, carries no flux, and drops out of the flux balance entirely.
    Only the FREE profile drives transport, and the free profile obeys the
    unretarded conductivity. What binding changes is CAPACITY, and capacity only
    shows up in a TRANSIENT.

    Write the homogenised problem in matrix-equivalent potential ``u`` (the
    concentration the interstitium would hold in equilibrium with the local state).
    Two coefficients appear, and they homogenise differently:

        conductance  sigma_eff : matrix sigma_m, inclusion K*D_cyto
                                 -> Maxwell / Hasselman & Johnson, ALREADY in this
                                    module, and UNAFFECTED by binding.
        capacity     s_eff     : matrix 1, inclusion K*(1 + K_b)
                                 -> an exact VOLUME AVERAGE, because the amount
                                    stored per unit volume at equilibrium is a
                                    volume average and nothing else.

        s_eff = (1 - phi)*1 + phi*K*(1 + K_b)

    and the apparent (transient) diffusivity is ``D_app = sigma_eff / s_eff``. That
    is the standard retardation-factor result, and `results/two_route_efflux.md`
    2.2 (iii) already named its absence as the model's missing "capacity term".

    THE MINIMAL ADDITIVE FORM, WHICH IS WHAT THIS FUNCTION RETURNS
    --------------------------------------------------------------
    This module's published output is a CONDUCTANCE, computed against a convention
    in which the no-binding capacity is 1. Relitigating that convention would
    change what the model already ships. So the retardation term is introduced as
    the RATIO of the capacity with binding to the capacity without it, which
    isolates exactly the new physics and nothing else:

        R = s_eff(K_b) / s_eff(K_b = 0)
          = [(1-phi) + phi*K*(1+K_b)] / [(1-phi) + phi*K]
          = 1 + phi*K*K_b / [(1-phi) + phi*K]

    and at the pre-registered ``K = 1`` the denominator is exactly 1, giving the
    quoted ``R = 1 + phi*K_b``. ``K_b = 0`` gives ``R = 1`` exactly, so the term is
    bitwise off by default.

    THE CONSEQUENCE THAT MUST NOT BE GLOSSED
    -----------------------------------------
    Because binding does not touch conductance, ``transcellular_gain`` -- which is
    a ratio of two conductances -- is **identically unchanged** by this term. Any
    claim that retardation "reduces paclitaxel's transcellular gain" is false as
    stated. What retardation reduces is the APPARENT diffusivity, and it reduces
    it for every compound that partitions into cells, lipophilic ones most.

    THE LIMIT OF VALIDITY, STATED HERE BECAUSE IT IS LOAD-BEARING
    -------------------------------------------------------------
    Volume-averaging the capacity assumes the cell interior is in LOCAL EQUILIBRIUM
    with the interstitium on the timescale of the measurement. A cell behind a shut
    membrane cannot store solute any faster than the membrane admits it, so for
    large ``membrane_resistance_ratio`` the full capacity is not accessed within a
    finite experiment and this factor OVER-states the retardation.
    :func:`equilibration_time_s` returns the timescale that has to be short
    compared with the experiment for this term to be applicable, and
    ``experiments/intracellular_retardation.py --stage algebra`` measures the
    breakdown directly rather than asserting it.
    """
    kb = float(bound_to_free_ratio)
    if kb < 0:
        raise ValueError(f"bound_to_free_ratio must be >= 0, got {kb!r}")
    phi = float(cell_volume_fraction)
    if not (0.0 <= phi < 1.0):
        raise ValueError(f"cell_volume_fraction must be in [0,1), got {phi!r}")
    k = float(cell_partition)
    if k < 0:
        raise ValueError(f"cell_partition must be >= 0, got {k!r}")
    s0 = (1.0 - phi) + phi * k
    if s0 <= 0.0:
        return 1.0
    return 1.0 + phi * k * kb / s0


def equilibration_time_s(
    bound_to_free_ratio: float,
    membrane_permeability_um_s: float,
    cell_radius_um: float = DEFAULT_CELL_RADIUS_UM,
    cell_partition: float = 1.0,
) -> float:
    """Time for one cell to load its bound pool through its own membrane, seconds.

    A sphere of radius ``a`` with capacity ``K*(1+K_b)`` per unit volume, fed
    through a membrane of permeability ``P``, has volume-to-area ratio ``a/3`` and
    so relaxes exponentially with time constant

        tau = (a/3) * K*(1 + K_b) / P.

    :func:`retardation_factor` assumes local equilibrium, i.e. ``tau`` SHORT
    compared with the experiment being modelled. Where ``tau`` is long the
    retardation factor is an upper bound rather than an estimate, and the honest
    reading is that the cell's capacity is only partly accessible. Returns
    ``inf`` at zero permeability, which is the correct statement that a sealed
    cell never equilibrates.
    """
    p = float(membrane_permeability_um_s)
    if p <= 0.0:
        return float("inf")
    kb = float(bound_to_free_ratio)
    return (float(cell_radius_um) / 3.0) * float(cell_partition) * (1.0 + kb) / p


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
    efflux_ratio: Optional[float] = None,
    bound_to_free_ratio: Optional[float] = None,
) -> Dict[str, float]:
    """End-to-end: engine interstitial D + (logP, HBD, radius) -> tissue D_eff.

    This is the one call an experiment or a wiring point needs. It is NOT called
    from `compile_drug.py`; the two-route model is default-off.

    ``efflux_ratio`` is OPT-IN and defaults to ``None``, which is bit-identical
    to omitting it. When supplied it is applied through
    :func:`efflux_reduced_partition` (``K -> K/ER``), which that function derives
    in closed form rather than assumes. Pass the value
    ``efflux.efflux_ratio_for(drug)`` returns; ``None`` there means *no published
    ratio*, which must stay ``None`` here rather than becoming 1.0, so that "no
    efflux measured" and "measured to have no efflux" remain distinguishable in
    the caller's reporting.

    ``bound_to_free_ratio`` (``K_b``) is OPT-IN and defaults to ``None``, which is
    bitwise identical to omitting it. When supplied it adds the intracellular
    retardation term of :func:`retardation_factor`, which acts on CAPACITY and
    therefore appears ONLY in the new ``apparent_*`` outputs. ``hindrance``,
    ``effective_D_um2s`` and ``transcellular_gain`` are conductances and are
    **unchanged by it, by construction** -- see :func:`retardation_factor` for why
    that is the physics and not an omission.
    """
    mem = membrane_permeability_um_s(
        logp, hbd, solute_radius_nm, temperature_K=temperature_K,
        hbd_decrement=hbd_decrement,
    )
    d_cyto = cytoplasm_diffusivity_um2s(
        solute_radius_nm, temperature_K=temperature_K, d_over_water=d_over_water
    )
    k_eff = float(cell_partition)
    if efflux_ratio is not None:
        k_eff = efflux_reduced_partition(k_eff, efflux_ratio)
    res = two_route_factor(
        cell_volume_fraction,
        interstitial_D_um2s,
        d_cyto,
        mem["permeability_um_s"],
        cell_radius_um=cell_radius_um,
        cell_partition=k_eff,
    )
    res.update({
        "logp": float(logp),
        "hbd": float(hbd),
        "solute_radius_nm": float(solute_radius_nm),
        "membrane_permeability_um_s": mem["permeability_um_s"],
        "log_k_membrane": mem["log_k_membrane"],
        "membrane_clamped": mem["clamped"],
        "cytoplasm_D_um2s": d_cyto,
        "cell_partition_passive": float(cell_partition),
        "cell_partition_effective": k_eff,
        "efflux_ratio": None if efflux_ratio is None else float(efflux_ratio),
    })

    # ---- intracellular retardation: CAPACITY, not conductance. Default-off. --
    if bound_to_free_ratio is None:
        r_factor = 1.0
    else:
        r_factor = retardation_factor(
            bound_to_free_ratio, cell_volume_fraction=cell_volume_fraction,
            cell_partition=k_eff)
    res.update({
        "bound_to_free_ratio": (None if bound_to_free_ratio is None
                                else float(bound_to_free_ratio)),
        "retardation_factor": r_factor,
        "apparent_hindrance": res["hindrance"] / r_factor,
        "apparent_D_um2s": res["effective_D_um2s"] / r_factor,
        "equilibration_time_s": (
            float("nan") if bound_to_free_ratio is None else
            equilibration_time_s(bound_to_free_ratio,
                                 mem["permeability_um_s"],
                                 cell_radius_um=cell_radius_um,
                                 cell_partition=k_eff)),
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

    # ---- [6] EFFLUX: K -> K/ER is DERIVED, not assumed -------------------------
    # Solve the three-compartment slab directly and check the closed form
    #     (C1-C2)/J = 2/P + (L/sigma_c)*ER,   ER = 1 + E/P
    # which is what makes `efflux_reduced_partition` an identity rather than a
    # modelling choice. Pure Python, no numpy.
    print("\n[6] efflux enters as K -> K/ER: direct slab solve vs the closed form")
    worst_e = 0.0
    for p_mem in (1e-3, 0.1, 1.0, 10.0, 1e3):
        for sig_c in (0.5, 20.0, 500.0):
            for length in (0.5, 5.0, 20.0):
                for e_pump in (0.0, 0.5, 10.0, 500.0):
                    c1, c2 = 1.0, 0.0
                    # J = P(C1-Ca) - E*Ca ; J = (sig/L)(Ca-Cb) ; J = P(Cb-C2)+E*Cb
                    j = p_mem * (c1 - c2) / (2.0 + (p_mem + e_pump) * length / sig_c)
                    ca = (p_mem * c1 - j) / (p_mem + e_pump)
                    cb = (j + p_mem * c2) / (p_mem + e_pump)
                    # the three statements must all hold at this J
                    assert abs(j - (p_mem * (c1 - ca) - e_pump * ca)) < 1e-9 * max(1.0, j)
                    assert abs(j - (sig_c / length) * (ca - cb)) < 1e-9 * max(1.0, j)
                    assert abs(j - (p_mem * (cb - c2) + e_pump * cb)) < 1e-9 * max(1.0, j)
                    er = 1.0 + e_pump / p_mem
                    r_direct = (c1 - c2) / j
                    r_form = 2.0 / p_mem + (length / sig_c) * er
                    worst_e = max(worst_e, abs(r_direct - r_form) / r_form)
    assert worst_e < 1e-12, worst_e
    print(f"    max relative error over 180 (P, sigma_c, L, E) points = {worst_e:.3g}")
    print("    -> efflux multiplies the INTERIOR resistance by ER and leaves the")
    print("       membrane resistance alone, i.e. K -> K/ER exactly. OK")

    # opt-in and bit-identical when off
    base = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
    assert two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=1.0)[
        "hindrance"] == base["hindrance"], "ER = 1 must be a no-op"
    assert two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=0.2)[
        "hindrance"] == base["hindrance"], "ER < 1 must never boost"

    # THE PREDICTION: efflux is invisible where the membrane is shut.
    shut = two_route_diffusivity(SIGMA_M, 0.539, 0.0, 6)          # doxorubicin-like
    shut_e = two_route_diffusivity(SIGMA_M, 0.539, 0.0, 6, efflux_ratio=1e4)
    assert abs(shut_e["hindrance"] / shut["hindrance"] - 1.0) < 1e-6, (
        "a shut membrane must be immune to efflux")
    open_ = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)        # paclitaxel-like
    open_e = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4, efflux_ratio=39.0)
    print(f"    shut membrane  (HBD 6, logP 0.0): gain "
          f"{shut['transcellular_gain']:.4f} -> {shut_e['transcellular_gain']:.4f} "
          f"at ER = 1e4")
    print(f"    open membrane  (HBD 4, logP 3.74): gain "
          f"{open_['transcellular_gain']:.4f} -> {open_e['transcellular_gain']:.4f} "
          f"at ER = 39")
    assert open_e["transcellular_gain"] < open_["transcellular_gain"] * 0.5

    # ---- [7] RETARDATION: capacity, not conductance ---------------------------
    # The whole content of this term is that binding changes the amount a tissue
    # can HOLD, not the rate at which it conducts. The two are separate
    # coefficients of the homogenised problem and they compose differently: the
    # conductance by Maxwell, the capacity by an exact volume average. The check
    # below is that the retardation factor really is the ratio of two volume
    # averages, and that it is arithmetically incapable of moving a conductance.
    print("\n[7] intracellular retardation: R = s_eff(K_b)/s_eff(0)")
    worst_r = 0.0
    for ph in (0.1, 0.4, PHI, 0.95):
        for kk in (0.5, 1.0, 3.0):
            for kb in (0.0, 0.3, 3.0, 30.0, 300.0, 3000.0):
                s0 = (1.0 - ph) + ph * kk
                s1 = (1.0 - ph) + ph * kk * (1.0 + kb)
                mine = retardation_factor(kb, cell_volume_fraction=ph,
                                          cell_partition=kk)
                worst_r = max(worst_r, abs(mine - s1 / s0) / (s1 / s0))
    assert worst_r < 1e-14, worst_r
    print(f"    max relative error over 72 (phi, K, K_b) points = {worst_r:.3g}")

    # default-off must be BITWISE off
    _b = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4)
    assert _b["retardation_factor"] == 1.0
    assert _b["apparent_D_um2s"] == _b["effective_D_um2s"]
    assert two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4,
                                 bound_to_free_ratio=None)["hindrance"] == \
        _b["hindrance"]

    # THE PREDICTION: no value of K_b may move a conductance.
    for _kb in (0.0, 1.0, 1e2, 1e4, 1e8):
        _r = two_route_diffusivity(SIGMA_M, 0.626, 3.74, 4,
                                   bound_to_free_ratio=_kb)
        assert _r["transcellular_gain"] == _b["transcellular_gain"]
        assert _r["hindrance"] == _b["hindrance"]
    print("    transcellular_gain is INVARIANT under K_b, as it must be:")
    print(f"      paclitaxel-like gain = {_b['transcellular_gain']:.4f} at K_b = 0")
    print(f"      paclitaxel-like gain = {_r['transcellular_gain']:.4f} at K_b = 1e8")
    print("    -> 'retardation reduces the transcellular gain' is FALSE as stated.")
    print("       What it reduces is the APPARENT diffusivity. See")
    print("       results/intracellular_retardation.md.")

    print(f"\n    the structural rule, and the row that condemns it:")
    print(f"      {'logP':>7}{'K_b':>10}{'bound %':>10}{'R':>10}")
    for _lp in (-2.62, 0.0, 2.0, 3.74):
        _kb2 = binding_capacity_from_logp(_lp)
        print(f"      {_lp:>7.2f}{_kb2:>10.2f}{_kb2/(1+_kb2):>9.1%}"
              f"{retardation_factor(_kb2):>10.2f}")
    print("      logP 0 -> 90% bound. Urea and mannitol are in the anchor set.")

    print(f"\nAll asserts passed. Runtime {_time.time() - t0:.2f} s")
