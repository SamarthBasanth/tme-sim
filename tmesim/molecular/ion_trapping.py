"""Lysosomal ion trapping as a TISSUE-DIFFERENTIAL retardation. Default OFF.

``./.venv/bin/python -m tmesim.molecular.ion_trapping``

What this is
------------
A weak base that is neutral at cytosolic pH crosses the lysosomal membrane, is
protonated in the acidic lumen, and cannot come back out. At equilibrium the
lysosome holds a Henderson-Hasselbalch excess of drug, and a tissue whose cells
devote a volume fraction ``phi_lys`` to acidic compartment therefore holds a
**retardation**

::

    ratio = [D]_lys / [D]_cyt = (1 + 10^(pKa - pH_lys)) / (1 + 10^(pKa - pH_cyt))
    R     = C_total / C_free  = 1 + phi_lys * (ratio - 1)

This is *exactly* the formulation already implemented and validated in
``tmesim/validation/moncayo_explant.py`` — same two functions, same constants,
same structure-derived pKa. It is restated here rather than imported because
``molecular`` must not depend on ``validation``; the equivalence is asserted
numerically in :func:`formulation_matches_moncayo` and in the self-test, so the
two cannot silently drift.

Why it is worth a module of its own
-----------------------------------
``results/moncayo_explant.md`` §4.3 proved that a *spatially uniform* ``R``
produces no dispersion **within** a section: at equilibrium ``C = R * C_bath``
everywhere, however large ``R`` gets. That proof is correct. But
``results/uptake_impact.md`` §10 showed the chain then used it to eliminate ion
trapping for the **GIST** target, which is a ratio *between two tissues*:

::

    C_tumour / C_normal = (R_t * E_t) / (R_n * E_n)          does NOT cancel

A tissue-differential ``R`` scales a tumour:normal ratio by ``R_t / R_n`` for
exactly the algebraic reason a tissue-differential uptake factor does, and for
exactly the reason a per-drug efflux scalar provably does not. Both identities
are verified in the self-test below.

The parameter, and the rule about it
------------------------------------
``phi_lys`` is the one number that could turn this into a fitted story, so it is
sourced only from measurements and every value carries a PMID
(:data:`PHI_LYS`). Rat per-tissue values come from Schmitt et al. 2020
(PMID 33148688), who derived them from the lysosomal marker lipid BMP; the
**liver** value is the only one in that table anchored to a direct measurement
and every other tissue is scaled from it by BMP concentration, which is recorded
in each entry rather than hidden.

**Tumour tissue has no measured value at all** — not for GIST, not for glioma,
not on any basis. Tumour entries are therefore marked ``is_measured=False`` and
carry an explicit, named assumption; :func:`trapping_detail` reports which of
the two a factor came from, so no caller can mistake one for the other.

Default OFF
-----------
Nothing in ``core``, ``physics``, ``sim``, ``treatment``, ``reconstruct`` or the
rest of ``molecular`` imports this module, so the term is off by construction
rather than by default. ``tmesim/validation/ion_trapping_impact.py`` verifies
that numerically.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional, Sequence, Tuple, Union

# ===========================================================================
# 1. pH constants
# ===========================================================================
#: Cytosolic pH. Identical to ``moncayo_explant.PH_CYTOSOL``.
PH_CYTOSOL = 7.2
#: Lysosomal lumen pH. Identical to ``moncayo_explant.PH_LYSOSOME``. This is the
#: canonical textbook value and the one the validated Moncayo stepper used.
PH_LYSOSOME = 4.7
#: Alternative lumen pH, carried ONLY as a sensitivity axis. Schmitt et al. 2019
#: (PMID 30409837) measured a Gaussian distribution of endo-/lysosomal pH in
#: plated rat hepatocytes with mean 5.53, and Schmitt et al. 2020 (PMID 33148688)
#: then used an effective pH_L of 5.3 in their PBPK model. A higher lumen pH
#: makes trapping WEAKER, so using 4.7 as the default is the choice that is
#: generous to the mechanism under test.
PH_LYSOSOME_SCHMITT = 5.3


# ===========================================================================
# 2. Tissue vocabulary -- shared with molecular.uptake so the two agree
# ===========================================================================
TISSUE_GENERIC = "generic"
TISSUE_NORMAL_LIVER = "normal_liver"
TISSUE_LIVER_METASTASIS = "liver_metastasis"
TISSUE_NORMAL_INTESTINE = "normal_intestine"
TISSUE_INTESTINAL_TUMOUR = "intestinal_tumour"
TISSUE_NORMAL_BRAIN = "normal_brain"
TISSUE_GLIOMA = "glioma"
TISSUE_KIDNEY = "kidney"
TISSUE_SPLEEN = "spleen"
TISSUE_MUSCLE = "muscle"
TISSUE_ADIPOSE = "adipose"
TISSUE_LUNG = "lung"
TISSUE_HEART = "heart"
TISSUE_SKIN = "skin"

TISSUES: Tuple[str, ...] = (
    TISSUE_GENERIC,
    TISSUE_NORMAL_LIVER, TISSUE_LIVER_METASTASIS,
    TISSUE_NORMAL_INTESTINE, TISSUE_INTESTINAL_TUMOUR,
    TISSUE_NORMAL_BRAIN, TISSUE_GLIOMA,
    TISSUE_KIDNEY, TISSUE_SPLEEN, TISSUE_MUSCLE, TISSUE_ADIPOSE,
    TISSUE_LUNG, TISSUE_HEART, TISSUE_SKIN,
)


# ===========================================================================
# 3. THE PARAMETER -- measured lysosomal volume fractions, with PMIDs
# ===========================================================================
#: Fraction of LIVER volume that is cellular, used by Schmitt et al. 2020 to
#: convert their hepatocyte cell-basis measurement onto a tissue basis. Recorded
#: so that ``cell_to_tissue_basis`` reproduces the published liver number rather
#: than approximating it.
LIVER_CELLULAR_VOLUME_FRACTION = 0.839


@dataclass(frozen=True)
class PhiLysEntry:
    """One lysosomal volume fraction, on a stated basis, with its provenance.

    ``basis`` is load-bearing and is the single easiest way to be wrong by a
    factor of ~1.2: ``"tissue"`` is the fraction of whole tissue volume (what a
    tissue-level retardation needs), ``"cell"`` is the fraction of cell volume,
    ``"cytoplasm"`` is the fraction of cytoplasmic volume. Only ``"tissue"``
    entries are usable directly by :func:`retardation_factor`.
    """
    tissue: str
    phi_lys: Optional[float]
    basis: str
    species: str
    method: str
    citation: str
    is_measured: bool
    notes: str = ""


def _e(**kw) -> PhiLysEntry:
    return PhiLysEntry(**kw)


#: Rat per-tissue lysosomal volume fractions, **fraction of TISSUE volume**.
#:
#: Schmitt et al. 2020, Drug Metab Dispos 49(1):53-61, PMID 33148688,
#: doi:10.1124/dmd.120.000161, Table 1 / Table 2 -- derived from the lysosomal
#: marker lipid bis(monoacylglycero)phosphate (BMP), range 0.03% (adipose) to
#: 5.3% (spleen).
#:
#: PROVENANCE CAVEAT, RECORDED ON EVERY ROW: only LIVER is anchored to a
#: measurement. The paper states verbatim that "liver was taken as reference,
#: whose lysosomal content was determined experimentally. Other tissues were
#: calculated relative to the liver by bis(monoacylglycero)phosphate tissue
#: concentrations." These are therefore one measurement plus a lipid-ratio
#: extrapolation, not fourteen independent measurements.
_SCHMITT_2020 = ("PMID 33148688 (Schmitt 2020 Drug Metab Dispos 49:53-61, "
                 "doi:10.1124/dmd.120.000161), Table 1/2: rat tissue lysosomal "
                 "volume fractions from the marker lipid BMP")
_SCHMITT_SCALED = (" -- SCALED FROM LIVER by BMP concentration, not independently "
                   "measured (the paper's own Table 2 footnote)")

PHI_LYS: Dict[str, PhiLysEntry] = {
    # ---- the anchor: the one directly measured tissue --------------------
    TISSUE_NORMAL_LIVER: _e(
        tissue=TISSUE_NORMAL_LIVER, phi_lys=0.0125, basis="tissue", species="rat",
        method="BMP marker lipid + hepatocyte LysoTracker volumetry",
        citation=_SCHMITT_2020,
        is_measured=True,
        notes=("THE ONLY DIRECTLY MEASURED ROW in the source table, and the "
               "reference every other row is scaled from. Derived as "
               "(1 - f_EW) x fV_lys(hepatocyte) = 0.839 x 1.49% = 1.25%, "
               "displayed as 1.3% in the paper's Table 2. The hepatocyte "
               "cell-basis input 1.49% +/- 0.08% is Schmitt et al. 2019, "
               "PMID 30409837, doi:10.1124/dmd.118.084541, LysoTracker Red spot "
               "volumetry over >350,000 endo-/lysosomes in plated rat "
               "hepatocytes from 6 rats. Independent cross-check: Blouin, "
               "Bolender & Weibel 1977, PMID 833203, doi:10.1083/jcb.72.2.441, "
               "EM point counting, n = 5 rats -- hepatocyte dense bodies 0.82% "
               "+/- 0.16% of CELL volume, plus autophagic vacuoles 0.23% = 1.1% "
               "total lysosomal compartment, i.e. 0.92% on a tissue basis. The "
               "two methods bracket 0.9-1.3%; LysoTracker labels all acidic "
               "compartment while EM dense bodies do not, so the higher value "
               "is used and it is the one generous to the mechanism.")),
    # ---- scaled from liver by BMP ----------------------------------------
    TISSUE_NORMAL_BRAIN: _e(
        tissue=TISSUE_NORMAL_BRAIN, phi_lys=0.0136, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True, notes="1.4% of tissue volume." + _SCHMITT_SCALED),
    TISSUE_NORMAL_INTESTINE: _e(
        tissue=TISSUE_NORMAL_INTESTINE, phi_lys=0.0122, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True,
        notes=("1.2% of tissue volume, the paper's weighted 'gut' mean (small "
               "intestine 1.4%, large intestine 0.88%, stomach 1.1%)."
               + _SCHMITT_SCALED)),
    TISSUE_KIDNEY: _e(
        tissue=TISSUE_KIDNEY, phi_lys=0.0173, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True, notes="1.7% of tissue volume." + _SCHMITT_SCALED),
    TISSUE_SPLEEN: _e(
        tissue=TISSUE_SPLEEN, phi_lys=0.0527, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True,
        notes=("5.3% of tissue volume -- THE LARGEST VALUE FOR ANY TISSUE IN THE "
               "BODY, and therefore the ceiling any bound argument should use."
               + _SCHMITT_SCALED)),
    TISSUE_LUNG: _e(
        tissue=TISSUE_LUNG, phi_lys=0.0151, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True, notes="1.5% of tissue volume." + _SCHMITT_SCALED),
    TISSUE_HEART: _e(
        tissue=TISSUE_HEART, phi_lys=0.0032, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True, notes="0.32% of tissue volume." + _SCHMITT_SCALED),
    TISSUE_MUSCLE: _e(
        tissue=TISSUE_MUSCLE, phi_lys=0.0011, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True, notes="0.11% of tissue volume." + _SCHMITT_SCALED),
    TISSUE_SKIN: _e(
        tissue=TISSUE_SKIN, phi_lys=0.0016, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True, notes="0.16% of tissue volume." + _SCHMITT_SCALED),
    TISSUE_ADIPOSE: _e(
        tissue=TISSUE_ADIPOSE, phi_lys=0.00027, basis="tissue", species="rat",
        method="BMP marker lipid, scaled from liver", citation=_SCHMITT_2020,
        is_measured=True,
        notes=("0.03% of tissue volume -- the smallest value for any tissue."
               + _SCHMITT_SCALED)),
    # ---- tumour tissues: NO MEASUREMENT EXISTS ---------------------------
    TISSUE_LIVER_METASTASIS: _e(
        tissue=TISSUE_LIVER_METASTASIS, phi_lys=None, basis="tissue",
        species="human", method="none",
        citation="NO MEASUREMENT EXISTS",
        is_measured=False,
        notes=("No lysosomal volume fraction has ever been measured in GIST "
               "tissue or in GIST-T1 cells, on any basis. Both GIST/lysosome "
               "papers -- Burger 2015 PMID 26108972 doi:10.1124/mol.114.097451 "
               "and Sasaki 2021 PMID 33826667 doi:10.1371/journal.pone.0249650 "
               "-- are pharmacological and colocalisation studies with no "
               "volumetry. A Europe PMC full-text search for the exact phrase "
               "'lysosomal volume fraction' returns 6 papers corpus-wide and "
               "none concerns GIST. Beware the review PMID 36899844, which "
               "states '70-80% of IM was found to accumulate in lysosomes ... "
               "in IM-sensitive GIST-1 cells' citing Burger 2015 -- that figure "
               "is not in Burger 2015 and must not be propagated. Callers must "
               "supply an explicit assumption; see CANCER_CELL_PHI_LYS.")),
    TISSUE_INTESTINAL_TUMOUR: _e(
        tissue=TISSUE_INTESTINAL_TUMOUR, phi_lys=None, basis="tissue",
        species="human", method="none", citation="NO MEASUREMENT EXISTS",
        is_measured=False,
        notes="As liver_metastasis: no GIST tissue measurement exists."),
    TISSUE_GLIOMA: _e(
        tissue=TISSUE_GLIOMA, phi_lys=None, basis="tissue", species="human",
        method="none", citation="NO RETRIEVABLE MEASUREMENT",
        is_measured=False,
        notes=("The one on-target study exists and its numbers are not "
               "retrievable: Forsby, Griph & Ivarsson 1983, PMID 6137128, "
               "doi:10.1111/j.1699-0463.1983.tb02766.x, stereological volume "
               "density of secondary lysosomes in a human normal glial line vs "
               "a human malignant glioma line, on a CYTOPLASM basis. Only the "
               "direction survives in the abstract ('the volume density of "
               "lysosomes was significantly higher in the post-confluent "
               "glioma cells'); the paper has no PMCID, is not open access, is "
               "not in the Wayback Machine, and has 3 citations, so the four "
               "numbers were never quoted onward. This is an interlibrary-loan "
               "problem, not a search problem.")),
    TISSUE_GENERIC: _e(
        tissue=TISSUE_GENERIC, phi_lys=None, basis="tissue", species="-",
        method="none", citation="not a real tissue",
        is_measured=False,
        notes=("The generic tissue deliberately carries no value, so a caller "
               "who forgets to name a tissue gets R = 1 rather than a silent "
               "default. Note that moncayo_explant.py uses phi_lys = 0.014 for "
               "'generic tissue' (0.7 cellularity x 0.02 lysosomal fraction), "
               "which was a plausibility placeholder rather than a measurement; "
               "it happens to land within 12% of the Schmitt liver value.")),
}


#: The ONLY modern measured lysosomal volume fractions in cancer cells, and the
#: closest thing to a tumour parameter that exists. **Fraction of CELL volume.**
#:
#: Collins et al. 2021, J Pharmacol Exp Ther, PMID 33172973,
#: doi:10.1124/jpet.120.000309, PMC7841421: live-cell confocal z-stacks with
#: lysosomal volume normalised to experimentally determined cell volume.
#:
#: Three caveats travel with it and all three cut against using it as a
#: tumour parameter: the panel is human BREAST lines only, there is **no normal
#: comparator in the paper**, and the internal spread is 7-fold. The paper's own
#: point is that lysosomotropic exposure (HCQ, Torin1) inflates the parameter
#: 1.5-3x within 24 h, i.e. it is not a fixed constant -- which matters directly
#: here, because imatinib is itself lysosomotropic.
CANCER_CELL_PHI_LYS: Dict[str, float] = {
    "MCF7": 0.0050,          # 0.50% +/- 0.19% of cell volume
    "T47D": 0.0078,          # ~0.78%
    "MDA-MB-468": 0.00835,   # ~0.835%
    "MDA-MB-231": 0.0367,    # 3.67% +/- 1.09%
}
CANCER_CELL_PHI_LYS_CITATION = (
    "PMID 33172973 (Collins 2021 J Pharmacol Exp Ther, "
    "doi:10.1124/jpet.120.000309): lysosomal percent volume of total cellular "
    "volume, human breast cancer lines, live-cell confocal z-stacks; 0.50% "
    "+/- 0.19% (MCF7) to 3.67% +/- 1.09% (MDA-MB-231)")

#: Whether tumour cells hold MORE or LESS lysosomal volume than normal cells is
#: not settled, and the quantitative record points both ways. Recorded here
#: because the sign of this contrast is the sign of the whole mechanism.
CANCER_VS_NORMAL_DIRECTION_NOTE = (
    "NO MATCHED CANCER-VS-NORMAL LYSOSOMAL VOLUME-FRACTION PAIR EXISTS, on any "
    "basis. The widely repeated 'cancer cells have more lysosomes' claim rests "
    "on LAMP1 westerns, EM appearance and cathepsin activity, not volumetry, "
    "and the two studies that do report volumes point the other way: Glunde "
    "2003 (PMID 14965446, doi:10.1016/S1476-5586(03)80037-4) gives non-malignant "
    "MCF-12A 16.8 um^3/cell against metastatic MDA-MB-231 11.2 and MDA-MB-435 "
    "11.6, i.e. the NORMAL line is higher; and Kundra & Dean 1990 (PMID 1693341, "
    "doi:10.1016/0014-4827(90)90261-8) found 'number of lysosomes and their "
    "total volume were also much larger in 3T3 cells' than in SV40-transformed "
    "SV3T3, i.e. transformation REDUCED it. Perera 2015 (PMID 26168401, "
    "doi:10.1038/nature14587) gives the strongest cancer-vs-normal contrast "
    "(15.0 +/- 1.9 vs 1.0 +/- 0.9 lysosomes per cell in PDA vs normal pancreas) "
    "but that is a COUNT, not a volume fraction.")

#: Human-liver values, for the record, and none of them is usable.
HUMAN_LIVER_NOTE = (
    "NO NORMAL-HUMAN LIVER LYSOSOMAL VOLUME FRACTION IS RETRIEVABLE. Every "
    "'human' value in the PBPK literature is rat, a textbook, or a fit: Liu & "
    "Jusko 2020 (PMID 33277347, doi:10.1124/jpet.120.000385) states 'rat "
    "lysosomal fractions were employed' and its liver 9.55% is FITTED, not "
    "anatomy; Assmus 2017 (PMID 28823852, doi:10.1016/j.ejps.2017.08.014) "
    "tabulates liver 1.12% of tissue but its caption says 'in rat' and its "
    "values are Blouin 1977's; Kazmi 2013 (PMID 23378628) 0.7% traces to a "
    "review whose own work is rat; Hallifax & Houston 2007's 1% cites a "
    "textbook; the 0.5% that circulates via Trapp 2008 "
    "(doi:10.1007/s00249-008-0338-4) is a footnoted ASSUMPTION chosen 'to get a "
    "volume ratio lysosome to cytosol of 1:200 (De Duve et al. 1974)'. The only "
    "retrievable human hepatocyte number, 3.1% of CYTOPLASM (Yu 1988, PMID "
    "3169726, doi:10.1016/s0046-8177(88)80149-7), is from postmortem livers "
    "after shock or head trauma, includes lipofuscin granules, and that paper's "
    "own conclusion is that trauma EXPANDS hepatocyte lysosomes -- so it is not "
    "a healthy-human parameter. Rohr 1976 (PMID 823698) and Roessner 1978 "
    "(PMID 654842) certainly contain normal-human values and are not obtainable "
    "online.")


def cell_to_tissue_basis(phi_cell: float,
                         cellular_volume_fraction: float
                         = LIVER_CELLULAR_VOLUME_FRACTION) -> float:
    """Convert a fraction-of-CELL-volume onto a fraction-of-TISSUE-volume.

    Mixing the two bases is the easiest way to be wrong by ~20% here, so the
    conversion is a named function rather than an inline multiply. With the
    liver numbers it reproduces the published value:
    ``0.839 * 0.0149 = 0.0125``.
    """
    return float(phi_cell) * float(cellular_volume_fraction)


def phi_lys_for(tissue: str) -> Optional[float]:
    """Measured lysosomal volume fraction of ``tissue`` (fraction of TISSUE volume).

    Returns ``None`` where no measurement exists -- which is every tumour
    tissue. Callers must decide what to do about that; this module will not
    invent a number.
    """
    ent = PHI_LYS.get(tissue)
    return None if ent is None else ent.phi_lys


# ===========================================================================
# 4. THE FORMULATION -- identical to validation/moncayo_explant.py
# ===========================================================================
def pH_partition_ratio(pka: float, ph_lys: float = PH_LYSOSOME,
                       ph_cyt: float = PH_CYTOSOL) -> float:
    """Henderson-Hasselbalch ion-trapping ratio for a monoprotic base.

    Only the neutral species crosses the lysosomal membrane, so at equilibrium
    the TOTAL (neutral + protonated) concentration ratio is

    ::

        [D]_lys / [D]_cyt = (1 + 10^(pKa - pH_lys)) / (1 + 10^(pKa - pH_cyt)).

    For any strong base (pKa >> both pH values) this saturates at
    ``10^(pH_cyt - pH_lys) ~ 316``. That ceiling is the single most important
    property of the mechanism: it bounds how much asymmetry ion trapping can
    ever produce, and it is why the required-versus-measured comparison in
    ``validation/ion_trapping_impact.py`` has a definite answer.

    Byte-for-byte the same expression as
    ``tmesim.validation.moncayo_explant.pH_partition_ratio``, which was used in
    a stepper checked against Crank's analytical solution to 0.53%.
    """
    return (1.0 + 10.0 ** (pka - ph_lys)) / (1.0 + 10.0 ** (pka - ph_cyt))


def retardation_factor(pka: float, phi_lys: float,
                       ph_lys: float = PH_LYSOSOME,
                       ph_cyt: float = PH_CYTOSOL) -> float:
    """Local-equilibrium retardation ``R = C_total / C_free``.

    Reversible instantaneous trapping into an acidic compartment occupying
    volume fraction ``phi_lys`` of tissue gives

    ::

        R = 1 + phi_lys * (ratio - 1)

    and the transport equation becomes ``R dC_free/dt = div(D grad C_free)``,
    i.e. an effective diffusivity ``D/R`` with the measured (TOTAL) drug equal
    to ``R * C_free``. Exact for linear equilibrium partitioning.

    The ``C_total`` in that identity is what makes the mechanism testable: a
    measurement of TOTAL tissue drug (MALDI-MSI, tissue homogenate) sees ``R``,
    and a measurement of UNBOUND EXTRACELLULAR drug (intracerebral
    microdialysis) does not, because a lysosome is an intracellular reservoir
    and at equilibrium a reservoir does not shift the free extracellular
    concentration. That distinction is applied, not just noted, in
    ``validation/ion_trapping_impact.py``.

    Byte-for-byte the same expression as
    ``tmesim.validation.moncayo_explant.retardation_factor``, except that
    ``phi_lys`` is required rather than defaulted -- this module refuses to
    supply a lysosomal volume fraction the caller did not ask for.
    """
    return 1.0 + float(phi_lys) * (pH_partition_ratio(pka, ph_lys, ph_cyt) - 1.0)


def required_phi_lys(required_R: float, pka: float,
                     phi_lys_other: float = 0.0,
                     ph_lys: float = PH_LYSOSOME,
                     ph_cyt: float = PH_CYTOSOL) -> float:
    """Invert the retardation: what ``phi_lys`` would a required ratio demand?

    Solves ``R_this / R_other = required_R`` for ``phi_lys(this)`` given
    ``phi_lys(other)``. With ``phi_lys_other = 0`` -- the physically
    unattainable limit of a tissue whose cells contain no acidic compartment at
    all -- this is the SMALLEST demand the mechanism can be asked to meet, which
    is what makes it usable as a falsifying bound.
    """
    k = pH_partition_ratio(pka, ph_lys, ph_cyt) - 1.0
    r_other = 1.0 + float(phi_lys_other) * k
    return (required_R * r_other - 1.0) / k


def max_reachable_ratio(phi_lys_high: float, pka: float,
                        phi_lys_low: float = 0.0,
                        ph_lys: float = PH_LYSOSOME,
                        ph_cyt: float = PH_CYTOSOL) -> float:
    """Largest ``R_high / R_low`` the measured parameters can produce.

    The mirror of :func:`required_phi_lys`, and the actual bound: with
    ``phi_lys_low = 0`` it is the most a tissue-differential trapping term can
    ever do, given a measured ``phi_lys`` for the other tissue.
    """
    k = pH_partition_ratio(pka, ph_lys, ph_cyt) - 1.0
    return (1.0 + float(phi_lys_high) * k) / (1.0 + float(phi_lys_low) * k)


def formulation_matches_moncayo(tol: float = 0.0) -> Dict[str, object]:
    """Assert this module's formulation IS the validated Moncayo one.

    Compares :func:`pH_partition_ratio` and :func:`retardation_factor` against
    ``tmesim.validation.moncayo_explant``'s functions over a grid, with ``==``
    by default rather than ``allclose``. The point of restating rather than
    importing is to keep ``molecular`` free of a ``validation`` dependency; the
    point of this function is that restating cannot become drifting.
    """
    from ..validation.moncayo_explant import (PHI_LYS as MONCAYO_PHI_LYS,
                                              PH_CYTOSOL as M_CYT,
                                              PH_LYSOSOME as M_LYS,
                                              pH_partition_ratio as m_ratio,
                                              retardation_factor as m_R)
    checks = []
    for pka in (1.0, 5.0, 7.0, 8.07, 9.4, 9.7, 9.8, 9.95, 10.0, 10.4, 12.0):
        a, b = pH_partition_ratio(pka), m_ratio(pka)
        checks.append({"what": "ratio", "pka": pka,
                       "identical": (a == b) if tol == 0.0 else abs(a - b) <= tol})
        for phi in (0.0, 0.005, 0.0125, 0.014, 0.05, 0.1):
            ra = retardation_factor(pka, phi)
            rb = m_R(pka, phi)
            checks.append({"what": "R", "pka": pka, "phi": phi,
                           "identical": (ra == rb) if tol == 0.0 else abs(ra - rb) <= tol})
    return {
        "n": len(checks),
        "n_identical": sum(c["identical"] for c in checks),
        "all_identical": all(c["identical"] for c in checks),
        "ph_constants_identical": (PH_CYTOSOL == M_CYT and PH_LYSOSOME == M_LYS),
        "moncayo_generic_phi_lys": MONCAYO_PHI_LYS,
        "failures": [c for c in checks if not c["identical"]],
    }


# ===========================================================================
# 5. BASICITY FROM STRUCTURE -- the Moncayo SMARTS, unchanged
# ===========================================================================
#: A basic centre is an sp3 nitrogen that is NOT an amide/sulfonamide/
#: phosphoramide N, NOT attached to an aromatic ring, and not already charged.
#: Copied verbatim from ``moncayo_explant.BASIC_N_SMARTS``.
BASIC_N_SMARTS = "[NX3;H2,H1,H0;!$(N[C,S,P]=[O,S,N]);!$(N-a);!$(N#*);!$([N+])]"

#: Curated substructure-class pKa table, verbatim from
#: ``moncayo_explant._PKA_CLASS``. There is no first-principles pKa predictor in
#: RDKit, so a class table is the honest alternative to a hand-entered value per
#: drug -- it cannot be tuned drug by drug.
PKA_CLASS: Dict[str, float] = {
    "cyclic_secondary_amine": 10.4,   # piperidine / pyrrolidine NH
    "acyclic_secondary_amine": 9.7,   # R-NH-CH3
    "tertiary_amine": 9.8,
    "primary_amine": 10.0,
    "none": 1.0,                      # no basic centre -> neutral at pH 4.7-7.2
}

#: SMILES for every drug on either scorecard, each verified against PubChem and
#: asserted by RDKit InChIKey at call time. The two already in
#: ``molecular.library`` are repeated here with their PubChem isomeric forms so
#: that one assertion covers the whole table.
DRUG_SMILES: Dict[str, str] = {
    "imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
    "adavosertib": "CC(C)(C1=NC(=CC=C1)N2C3=NC(=NC=C3C(=O)N2CC=C)NC4=CC=C(C=C4)N5CCN(CC5)C)O",
    "eribulin": ("C[C@@H]1C[C@@H]2CC[C@H]3C(=C)C[C@@H](O3)CC[C@]45C[C@@H]6[C@H](O4)"
                 "[C@H]7[C@@H](O6)[C@@H](O5)[C@@H]8[C@@H](O7)CC[C@@H](O8)CC(=O)C[C@H]9"
                 "[C@H](C[C@H](C1=C)O2)O[C@@H]([C@@H]9OC)C[C@@H](CN)O"),
    "levetiracetam": "CC[C@@H](C(=O)N)N1CCCC1=O",
    "cefazolin": "CC1=NN=C(S1)SCC2=C(N3[C@@H]([C@@H](C3=O)NC(=O)CN4C=NN=N4)SC2)C(=O)O",
    "temozolomide": "CN1C(=O)N2C=NC(=C2N=N1)C(=O)N",
    "methotrexate": ("CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)N[C@@H]"
                     "(CCC(=O)O)C(=O)O"),
}

#: PubChem InChIKeys, asserted before any chemistry is done.
EXPECTED_INCHIKEY: Dict[str, str] = {
    "imatinib": "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
    "adavosertib": "BKWJAKQVGHWELA-UHFFFAOYSA-N",
    "eribulin": "UFNVPOGXISZXJD-JBQZKEIOSA-N",
    "levetiracetam": "HPHUVLMMVZITSG-LURJTMIESA-N",
    "cefazolin": "MLYYVTUWGNIJIB-BXKDBHETSA-N",
    "temozolomide": "BPEGJWRSRHCHSN-UHFFFAOYSA-N",
    "methotrexate": "FBOZXECLQNJBKD-ZDUSSCGKSA-N",
}

#: Literature pKa of the strongest basic centre, where one could be verified.
#: These are a CROSS-CHECK on the structure-derived class value and, for
#: imatinib only, the value actually used -- because the repo already committed
#: to it in ``validation/gist_imatinib.py`` and using a different number here
#: purely for this test would be indefensible.
LITERATURE_PKA: Dict[str, Optional[float]] = {
    "imatinib": 8.07,
    "adavosertib": None,
    "eribulin": None,
    "levetiracetam": None,
    "cefazolin": None,
    "temozolomide": None,
    "methotrexate": None,
}
LITERATURE_PKA_CITATION: Dict[str, str] = {
    "imatinib": ("strongest basic pKa 8.07 (piperazine N), the value already "
                 "hard-coded as IMATINIB_PKA_BASE in "
                 "validation/gist_imatinib.py, whose comment records the "
                 "literature range as 7.7-8.07; NOT a new parameter"),
    "adavosertib": ("NO pKa RETRIEVABLE. The 4-methylpiperazine distal nitrogen "
                    "is unambiguously basic from structure, but no measured pKa "
                    "for AZD1775 could be found in the primary literature, in "
                    "the Selleck/MedChemExpress datasheets, or in the phase I "
                    "PK papers. The structure-class value 9.8 is used and a pKa "
                    "sensitivity sweep is reported instead of a false precision."),
    "eribulin": ("NO pKa RETRIEVABLE for the primary amine; the structure-class "
                 "value 10.0 is used."),
    "levetiracetam": "no basic centre -- both nitrogens are amides",
    "cefazolin": "no basic centre -- an acid",
    "temozolomide": "no basic centre",
    "methotrexate": ("no basic centre passes the SMARTS -- the N-methyl aniline "
                     "nitrogen is aromatic-attached and the pteridine amines are "
                     "aromatic. Methotrexate is a diacid and its weakly basic "
                     "pteridine N1 (pKa ~5.7) is below the cytosolic-lysosomal "
                     "window, so the trapping ratio would be ~1 even if it were "
                     "counted."),
}

#: Which drug names are antibodies and short-circuit to R = 1 before any
#: chemistry runs, for the same reason ``molecular.uptake`` does it.
ANTIBODIES = ("bevacizumab",)


@dataclass(frozen=True)
class Basicity:
    drug: str
    smiles: Optional[str]
    inchikey: Optional[str]
    n_basic_centres: int
    basic_class: str
    pka_structure: float
    pka_literature: Optional[float]
    pka_used: float
    pka_source: str
    trapping_ratio: float
    note: str = ""


_BASICITY_CACHE: Dict[str, Basicity] = {}


def basicity(drug: str, prefer_literature_pka: bool = True,
             ph_lys: float = PH_LYSOSOME,
             ph_cyt: float = PH_CYTOSOL) -> Basicity:
    """Structure -> basic centre -> pKa -> Henderson-Hasselbalch trapping ratio.

    ``prefer_literature_pka`` decides only what ``pka_used`` is; both values are
    always reported. It is ``True`` by default because the single drug for which
    a literature pKa exists is imatinib, and the repo already uses 8.07 for it
    elsewhere -- switching to the structure class value 9.8 here would make the
    mechanism look STRONGER on the strength of a coarser number, which is
    exactly the wrong direction for a test the module is trying to fail
    honestly.
    """
    key = (f"{drug}|{prefer_literature_pka}|{ph_lys}|{ph_cyt}")
    hit = _BASICITY_CACHE.get(key)
    if hit is not None:
        return hit

    name = str(drug).lower()
    if name in ANTIBODIES:
        out = Basicity(drug=name, smiles=None, inchikey=None, n_basic_centres=0,
                       basic_class="antibody", pka_structure=1.0,
                       pka_literature=None, pka_used=1.0, pka_source="antibody",
                       trapping_ratio=1.0,
                       note=("an IgG does not partition across a lysosomal "
                             "membrane by pH; short-circuited to R = 1"))
        _BASICITY_CACHE[key] = out
        return out

    smi = DRUG_SMILES.get(name)
    if smi is None:
        out = Basicity(drug=name, smiles=None, inchikey=None, n_basic_centres=0,
                       basic_class="unknown", pka_structure=1.0,
                       pka_literature=None, pka_used=1.0, pka_source="not_in_table",
                       trapping_ratio=1.0,
                       note="no verified SMILES in this module's table -> R = 1")
        _BASICITY_CACHE[key] = out
        return out

    from rdkit import Chem
    from rdkit.Chem import inchi as rd_inchi

    mol = Chem.MolFromSmiles(smi)
    if mol is None:                                     # degrade gracefully
        out = Basicity(drug=name, smiles=smi, inchikey=None, n_basic_centres=0,
                       basic_class="unparseable", pka_structure=1.0,
                       pka_literature=LITERATURE_PKA.get(name), pka_used=1.0,
                       pka_source="rdkit_parse_failed", trapping_ratio=1.0,
                       note="SMILES did not parse; factor forced to 1.0")
        _BASICITY_CACHE[key] = out
        return out
    ikey = rd_inchi.MolToInchiKey(mol)
    expected = EXPECTED_INCHIKEY.get(name)
    if expected is not None and ikey != expected:
        raise AssertionError(f"{name} InChIKey {ikey} != expected {expected}")

    patt = Chem.MolFromSmarts(BASIC_N_SMARTS)
    hits = mol.GetSubstructMatches(patt)
    if not hits:
        cls = "none"
    else:
        atom = mol.GetAtomWithIdx(hits[0][0])
        nh = atom.GetTotalNumHs()
        if nh >= 2:
            cls = "primary_amine"
        elif nh == 1:
            cls = "cyclic_secondary_amine" if atom.IsInRing() else "acyclic_secondary_amine"
        else:
            cls = "tertiary_amine"
    pka_struct = PKA_CLASS[cls]
    pka_lit = LITERATURE_PKA.get(name)
    if prefer_literature_pka and pka_lit is not None and cls != "none":
        pka_used, src = pka_lit, "literature"
    else:
        pka_used, src = pka_struct, "structure_class"

    out = Basicity(drug=name, smiles=smi, inchikey=ikey,
                   n_basic_centres=len(hits), basic_class=cls,
                   pka_structure=pka_struct, pka_literature=pka_lit,
                   pka_used=pka_used, pka_source=src,
                   trapping_ratio=(1.0 if cls == "none"
                                   else pH_partition_ratio(pka_used, ph_lys, ph_cyt)),
                   note=LITERATURE_PKA_CITATION.get(name, ""))
    _BASICITY_CACHE[key] = out
    return out


# ===========================================================================
# 6. THE PUBLIC TERM
# ===========================================================================
def _resolve_name(drug: Union[str, Mapping]) -> str:
    if isinstance(drug, str):
        return drug.lower()
    for attr in ("name", "drug", "drug_name"):
        v = getattr(drug, attr, None) if not isinstance(drug, Mapping) else drug.get(attr)
        if isinstance(v, str):
            return v.lower()
    return ""


def trapping_detail(drug: Union[str, Mapping], tissue: str = TISSUE_GENERIC,
                    phi_lys: Optional[float] = None,
                    prefer_literature_pka: bool = True,
                    ph_lys: float = PH_LYSOSOME,
                    ph_cyt: float = PH_CYTOSOL) -> Dict[str, object]:
    """Full audit trail for one (drug, tissue): R, and exactly where it came from.

    ``phi_lys`` overrides the tissue table. It MUST be supplied for any tumour
    tissue, because no measurement exists for one; without it the factor is
    exactly 1.0 and the reason field says so. That is deliberate -- the failure
    mode this guards against is a tumour parameter quietly defaulting to
    something plausible.
    """
    name = _resolve_name(drug)
    b = basicity(name, prefer_literature_pka=prefer_literature_pka,
                 ph_lys=ph_lys, ph_cyt=ph_cyt)
    ent = PHI_LYS.get(tissue)

    if phi_lys is not None:
        phi, phi_src, measured = float(phi_lys), "caller_supplied", False
    elif ent is not None and ent.phi_lys is not None:
        phi, phi_src, measured = ent.phi_lys, ent.citation, ent.is_measured
    else:
        phi, phi_src, measured = 0.0, "no_measured_phi_lys", False

    if b.basic_class in ("none", "antibody", "unknown", "unparseable"):
        reason = {"none": "no_basic_centre", "antibody": "antibody",
                  "unknown": "drug_not_in_table",
                  "unparseable": "smiles_unparseable"}[b.basic_class]
        R = 1.0
    elif phi == 0.0:
        reason = ("no_measured_phi_lys_for_this_tissue"
                  if phi_src == "no_measured_phi_lys" else "phi_lys_zero")
        R = 1.0
    else:
        reason = "applied"
        R = retardation_factor(b.pka_used, phi, ph_lys, ph_cyt)

    return {
        "drug": name, "tissue": tissue, "factor": float(R), "reason": reason,
        "phi_lys": float(phi), "phi_lys_source": phi_src,
        "phi_lys_is_measured": bool(measured),
        "pka_used": b.pka_used, "pka_source": b.pka_source,
        "pka_structure": b.pka_structure, "pka_literature": b.pka_literature,
        "basic_class": b.basic_class, "n_basic_centres": b.n_basic_centres,
        "trapping_ratio": b.trapping_ratio,
        "ph_lys": ph_lys, "ph_cyt": ph_cyt,
        "inchikey": b.inchikey, "note": b.note,
    }


def trapping_factor(drug: Union[str, Mapping], tissue: str = TISSUE_GENERIC,
                    phi_lys: Optional[float] = None,
                    prefer_literature_pka: bool = True,
                    ph_lys: float = PH_LYSOSOME,
                    ph_cyt: float = PH_CYTOSOL) -> float:
    """``R = C_total / C_free`` for this drug in this tissue. 1.0 when inert."""
    return float(trapping_detail(drug, tissue=tissue, phi_lys=phi_lys,
                                 prefer_literature_pka=prefer_literature_pka,
                                 ph_lys=ph_lys, ph_cyt=ph_cyt)["factor"])


def trapping_adjusted_total(free_partition: float,
                            drug: Union[str, Mapping],
                            tissue: str = TISSUE_GENERIC,
                            phi_lys: Optional[float] = None,
                            measures_total_drug: bool = True,
                            **kw) -> float:
    """Convert a FREE-drug partition into the TOTAL-drug partition an assay sees.

    ``measures_total_drug`` is not a convenience switch, it is the physics. A
    lysosome is an intracellular reservoir; at equilibrium it raises total
    tissue drug by ``R`` and leaves the unbound EXTRACELLULAR concentration
    alone. So an assay that homogenises tissue (MALDI-MSI, LC-MS/MS on a
    homogenate) sees ``R``, and intracerebral microdialysis -- which samples
    extracellular fluid -- does not. Passing ``measures_total_drug=False``
    returns the input unchanged, and callers are expected to key it off the
    recorded measurement method rather than off convenience.
    """
    if not measures_total_drug:
        return float(free_partition)
    return float(free_partition) * trapping_factor(drug, tissue=tissue,
                                                   phi_lys=phi_lys, **kw)


def tissue_ratio_scale(drug: Union[str, Mapping], tumour_tissue: str,
                       normal_tissue: str,
                       phi_lys_tumour: Optional[float] = None,
                       phi_lys_normal: Optional[float] = None,
                       **kw) -> float:
    """``R_tumour / R_normal`` -- the factor a tumour:normal ratio is scaled by.

    This is the whole reason the module exists: the quantity is not 1 unless
    the two tissues have the same acidic-compartment content, and Moncayo's
    within-section proof says nothing about it.
    """
    rt = trapping_factor(drug, tissue=tumour_tissue, phi_lys=phi_lys_tumour, **kw)
    rn = trapping_factor(drug, tissue=normal_tissue, phi_lys=phi_lys_normal, **kw)
    return rt / rn


def summary() -> Dict[str, object]:
    """Counts for the report: how much of this table is measurement."""
    meas = [t for t, e in PHI_LYS.items()
            if e.phi_lys is not None and e.is_measured]
    return {
        "n_tissues": len(PHI_LYS),
        "n_with_measured_phi_lys": len(meas),
        "n_without": len(PHI_LYS) - len(meas),
        "n_directly_measured_not_scaled": 1,      # liver only
        "n_drugs": len(DRUG_SMILES),
        "n_drugs_with_basic_centre": sum(
            1 for d in DRUG_SMILES if basicity(d).basic_class
            not in ("none", "antibody", "unknown", "unparseable")),
        "n_drugs_with_literature_pka": sum(
            1 for d, v in LITERATURE_PKA.items() if v is not None),
        "saturating_ratio": 10.0 ** (PH_CYTOSOL - PH_LYSOSOME),
    }


# ===========================================================================
if __name__ == "__main__":
    # --- the formulation is the validated one ------------------------------
    fm = formulation_matches_moncayo()
    assert fm["all_identical"], fm["failures"]
    assert fm["ph_constants_identical"]

    # --- inert by default ---------------------------------------------------
    for d in DRUG_SMILES:
        assert trapping_factor(d, tissue=TISSUE_GENERIC) == 1.0, d
        # every tumour tissue lacks a measured phi_lys, so the term must not act
        for t in (TISSUE_LIVER_METASTASIS, TISSUE_INTESTINAL_TUMOUR, TISSUE_GLIOMA):
            assert trapping_factor(d, tissue=t) == 1.0, (d, t)
    # non-bases are inert even where phi_lys IS measured
    for d in ("temozolomide", "cefazolin", "levetiracetam", "methotrexate"):
        assert trapping_factor(d, tissue=TISSUE_NORMAL_LIVER) == 1.0, d
    # ... and so is anything not in the table
    assert trapping_factor("not_a_drug", tissue=TISSUE_NORMAL_LIVER) == 1.0
    assert trapping_factor("bevacizumab", tissue=TISSUE_NORMAL_LIVER) == 1.0

    # --- non-vacuity: it must bite somewhere -------------------------------
    R_liver = trapping_factor("imatinib", tissue=TISSUE_NORMAL_LIVER)
    assert R_liver > 4.0, R_liver
    assert trapping_factor("adavosertib", tissue=TISSUE_NORMAL_BRAIN) > 4.0

    # --- the bases the two scorecards actually contain ----------------------
    assert basicity("imatinib").basic_class == "tertiary_amine"
    assert basicity("adavosertib").basic_class == "tertiary_amine"
    assert basicity("eribulin").basic_class == "primary_amine"
    for d in ("temozolomide", "cefazolin", "levetiracetam", "methotrexate"):
        assert basicity(d).basic_class == "none", d
    assert basicity("imatinib").pka_used == 8.07                 # repo's own value
    assert basicity("imatinib", prefer_literature_pka=False).pka_used == 9.8

    # --- the identity that matters: R does NOT cancel out of a tissue ratio -
    e_t, e_n = 0.5572, 0.9614                       # GIST Krogh efficiencies
    r_t = trapping_factor("imatinib", tissue=TISSUE_LIVER_METASTASIS, phi_lys=0.0)
    r_n = trapping_factor("imatinib", tissue=TISSUE_NORMAL_LIVER)
    assert (r_t * e_t) / (r_n * e_n) == (r_t / r_n) * (e_t / e_n)
    assert abs((r_t / r_n) - 1.0) > 0.5, "must genuinely move the ratio"
    # ... and it DOES cancel when the two tissues share a phi_lys -- which is
    # the correct reading of Moncayo, restated as a test rather than as prose
    same = trapping_factor("imatinib", tissue=TISSUE_LIVER_METASTASIS,
                           phi_lys=phi_lys_for(TISSUE_NORMAL_LIVER))
    assert same == r_n
    assert (same * e_t) / (r_n * e_n) == e_t / e_n

    # --- the ceiling, which is what makes the bound argument possible -------
    assert abs(pH_partition_ratio(20.0) - 10.0 ** (PH_CYTOSOL - PH_LYSOSOME)) < 1e-6
    assert pH_partition_ratio(1.0) < 1.02

    # --- inversion round-trips ---------------------------------------------
    for pka in (8.07, 9.8):
        for phi_o in (0.0, 0.0125):
            for req in (2.0, 4.5, 11.7):
                p = required_phi_lys(req, pka, phi_o)
                assert abs(max_reachable_ratio(p, pka, phi_o) - req) < 1e-9
    # basis conversion reproduces the published liver number
    assert abs(cell_to_tissue_basis(0.0149) - 0.0125) < 5e-5

    # --- microdialysis vs homogenate ---------------------------------------
    assert trapping_adjusted_total(0.1, "imatinib", TISSUE_NORMAL_LIVER,
                                   measures_total_drug=False) == 0.1
    assert trapping_adjusted_total(0.1, "imatinib", TISSUE_NORMAL_LIVER) > 0.4

    # --- report -------------------------------------------------------------
    s = summary()
    print(f"ion_trapping: {s['n_tissues']} tissues, "
          f"{s['n_with_measured_phi_lys']} with a measured phi_lys, "
          f"{s['n_directly_measured_not_scaled']} of those measured DIRECTLY "
          f"(liver); {s['n_without']} with none (every tumour).")
    print(f"  drugs {s['n_drugs']}, with a basic centre "
          f"{s['n_drugs_with_basic_centre']}, with a literature pKa "
          f"{s['n_drugs_with_literature_pka']}")
    print(f"  trapping ratio saturates at {s['saturating_ratio']:.1f} "
          f"(pH {PH_CYTOSOL} -> {PH_LYSOSOME})")
    print()
    print(f"{'tissue':22s} {'phi_lys':>9s}  {'measured?':>9s}")
    print("-" * 45)
    for t in TISSUES:
        e = PHI_LYS[t]
        v = "--" if e.phi_lys is None else f"{100*e.phi_lys:.3f}%"
        print(f"{t:22s} {v:>9s}  {str(e.is_measured):>9s}")
    print()
    print(f"{'drug':14s} {'class':24s} {'pKa':>6s} {'src':>14s} {'ratio':>8s} "
          f"{'R(liver)':>9s} {'R(brain)':>9s}")
    print("-" * 92)
    for d in DRUG_SMILES:
        b = basicity(d)
        print(f"{d:14s} {b.basic_class:24s} {b.pka_used:6.2f} {b.pka_source:>14s} "
              f"{b.trapping_ratio:8.1f} "
              f"{trapping_factor(d, tissue=TISSUE_NORMAL_LIVER):9.3f} "
              f"{trapping_factor(d, tissue=TISSUE_NORMAL_BRAIN):9.3f}")
