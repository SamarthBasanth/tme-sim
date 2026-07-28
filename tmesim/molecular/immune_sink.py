"""Immune-cell (TAM) drug sequestration as a SPATIALLY-VARYING first-order sink.

WHAT THIS IS
------------
The reconstruction places real IMMUNE cells — 11.4-19.5% of the segmented cells
across the three real sections under the shipping ``matched_null_z`` classifier
(``results/classifier_fix.md`` §3), gated on PTPRC / CD68 / AIF1 / CD163 and
dominated on the IO panel by CD163, whose detection rate there is 18.2%
(``results/classifier_fix.md`` §1) — but until now they were pharmacologically
INERT: they occupied space and consumed nothing. This module
adds the one thing a microenvironment model can say about them that a
compartmental model cannot: **drug is lost where the macrophages are**, so the
loss is a FIELD rather than a number.

    dC/dt = div(D_eff grad C) - div(vC) + S - (clearance + binding) C
                                            - k_sink(x, y) C          <-- NEW

with

    k_sink(x, y) = k_max * phi_imm(x, y)                        [1/s]

* ``phi_imm`` is the local IMMUNE-cell AREA FRACTION in [0, 1], estimated from
  the MEASURED immune-cell centroids (:func:`immune_area_fraction`). It is
  dimensionless: cells/um^2 (a Gaussian kernel density of the measured points)
  times the per-cell cross-sectional area ``pi * r^2``.
* ``k_max`` [1/s] is the first-order loss constant that a voxel packed ENTIRELY
  with macrophages would impose. It is the only new number, and it is NOT
  fitted: it is swept (see :data:`SEQUESTRATION_SWEEP`) and reported as a
  sensitivity, because the literature supports the DIRECTION but not the
  MAGNITUDE for free drug (see "WHAT THE LITERATURE ACTUALLY SAYS" below).

Units check: ``[k_max] = 1/s``, ``[phi_imm] = 1`` (dimensionless), so
``[k_sink] = 1/s`` and ``k_sink * C`` is uM/s — the same shape as the engine's
existing ``(clearance + binding) * C`` term, which is exactly the point: this is
that term made SPATIAL.

DEFAULT: OFF, AND OFF BY CONSTRUCTION
-------------------------------------
Following the escalation set by ``molecular/uptake.py`` and
``molecular/ion_trapping.py``: **no shipped code path imports this module.**
Nothing under ``core/``, ``physics/``, ``sim/``, ``treatment/`` or
``reconstruct/`` was modified or even touched, so the term is off by
construction rather than off by default. Additionally
:func:`transport_with_immune_sink` with ``sink_rate=None`` (or an all-zero
field) is verified BIT-IDENTICAL to ``physics.diffusion.simulate_drug_transport``
in ``validation/immune_sink_impact.py`` §1 and ``tests/test_immune_sink.py``.

WHAT THE LITERATURE ACTUALLY SAYS — read this before quoting a magnitude
-----------------------------------------------------------------------
This was searched before it was coded, and the honest summary is that the
well-quantified TAM drug-handling literature is about NANOPARTICLE and LIPOSOMAL
formulations, and that where it is quantified the sign is a **DEPOT** (TAMs
help) at least as often as a **SINK** (TAMs hurt). For FREE small-molecule
cytotoxics the evidence is directional at best and in one direct measurement is
NEGATIVE. The table:

  1. **Miller et al., Nat Commun 2015;6:8692 (PMID 26503691, doi
     10.1038/ncomms9692) — the best-quantified case, and it is a DEPOT.**
     Single-cell imaging of a PLGA-b-PEG nanoparticle carrying a fluorescent
     Pt(IV) pro-drug. Measured, verbatim: TAMs were **4% of cells in the tumour
     mass but accumulated 30% of the total intratumoral nanoparticle** (and
     >40% relative to tumour cells); payload concentration was **>3-fold
     elevated within roughly one cell-length (15 um) of a TAM**, with gradients
     resolved over **<25 um**; and **macrophage depletion LOWERED intratumoral
     Pt by 24% (5.0 -> 3.8 uM) and made the tumours grow FASTER**. So in that
     system the macrophages were a slow-release reservoir that fed neighbouring
     tumour cells, not a sink that starved them. Two numbers here are used
     below: the 4%/30% partition (:data:`MILLER2015_TAM_CELL_FRACTION` /
     :data:`MILLER2015_TAM_DRUG_FRACTION`) as the ONLY quantitative anchor for
     how large a diverted-drug fraction is physically plausible, and the <25 um
     gradient length as the anchor for :data:`IMMUNE_KERNEL_BANDWIDTH_UM`.

  2. **Islam et al., Nanotheranostics 2022;6(4):451-464 (PMID 36105861, doi
     10.7150/ntno.75045) — the direct free-doxorubicin measurement, and it is
     NEGATIVE.** In a WEHI-164 fibrosarcoma heavily infiltrated by TAMs, both
     liposomal doxorubicin and a liposomal alendronate/doxorubicin co-formulation
     "shifted cellular drug uptake to TAM and to monocytic MDSC", while **"free
     doxorubicin cellular drug uptake was below detectable levels"**. i.e. the
     TAM-uptake phenomenon in that model is a property of the CARRIER, not of
     doxorubicin.

  3. **Mayer et al., J Pharmacol Exp Ther 1997;280(3):1406-1414 (PMID 9067330)
     — the tumour-level test, and it too fails to find TAM uptake.** Two murine
     fibrosarcoma lines differing 4-fold in TAM content. TAM-rich tumours took
     up 1.5-2x MORE liposomal doxorubicin, but fluorescence microscopy showed a
     **poor correlation between CD11b+ cells and doxorubicin signal**, and the
     authors concluded that "uptake of liposomal drugs by TAMs does not account
     for the enhanced accumulation", attributing it instead to TAM-mediated
     vascular permeability.

  4. **What DOES support the direction, mechanistically, for free doxorubicin.**
     Doxorubicin is a hydrophobic weak base (amino-sugar pKa ~8.3) and is a
     documented substrate for passive ion-trapping into acidic organelles:
     Zhitomirsky & Assaraf, Oncotarget 2015;6(2):1143-1156 (PMID 25544758, doi
     10.18632/oncotarget.2732) show nanomolar doxorubicin drives TFEB-mediated
     lysosomal biogenesis and that lysosome-sequestered drug is held away from
     its nuclear target; reviewed in Zhitomirsky & Assaraf, Drug Resist Updat
     2016;24:23-33 (PMID 26830313, doi 10.1016/j.drup.2015.11.004). Macrophages
     are the most lysosome-rich cell in the niche, so an acidic-compartment sink
     should be LARGER in them. Consistent with this, Aghai & Tokes, Cancer
     Chemother Pharmacol 1990;25(4):295-298 (PMID 2295115, doi
     10.1007/BF00684889) imaged doxorubicin compartmentalisation in normal human
     marrow and found the nuclei of **macrophages and activated phagocytes were
     the LEAST accessible to doxorubicin** of any lineage — i.e. in macrophages
     the drug is held in the cytoplasmic/vesicular compartment rather than
     reaching chromatin, which is what a sequestering sink looks like at the
     single-cell level. NOTE that this repo already models the intracellular
     half of that mechanism separately and generically in
     ``molecular/ion_trapping.py``; this module is the EXTRACELLULAR,
     SPATIALLY-RESOLVED complement (loss from the interstitium where the
     macrophages are), not a duplicate of it.

  5. **The direction is also asserted at the systems level.** Miller's own
     discussion notes that "TAMs also provide chemoprotective cues to tumour
     cells, and TAM infiltration following chemotherapy correlates with poor
     treatment response" (PMID 26503691). That is a correlation and a mechanism
     statement, not an uptake ratio.

**Therefore: no drug in :data:`SUBSTRATE_TABLE` carries a measured free-drug
TAM uptake ratio, and :func:`immune_sink_strength_for` returns exactly 0.0 for
every drug.** The magnitude is a SWEPT parameter, never a fitted one, and every
number this module produces must be quoted with the strength it was produced at.

WHAT WOULD MAKE THE MAGNITUDE MEASURED
--------------------------------------
A single experiment closes it: free doxorubicin (or its fluorescent analogue)
dosed into a TAM-containing tumour, then flow-sorted CD163+/TREM2+ macrophages
versus tumour cells assayed for drug content per cell, giving a per-cell
enrichment ratio; or single-cell intravital imaging of the free drug rather than
the carrier, as Miller 2015 did for the carrier. Until then this term has a
literature-supported SIGN and an unmeasured SIZE.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import factorized

from ..core import (
    CellType,
    DrugParams,
    Microenvironment,
    PDEConfig,
    TransportResult,
)
from ..physics.diffusion import (
    _advective_divergence,
    _build_diffusion_operator,
    _effective_diffusivity,
)

__all__ = [
    "IMMUNE_CELL_RADIUS_UM",
    "IMMUNE_KERNEL_BANDWIDTH_UM",
    "MILLER2015_TAM_CELL_FRACTION",
    "MILLER2015_TAM_DRUG_FRACTION",
    "MILLER2015_DEPLETION_DROP",
    "MILLER2015_GRADIENT_LENGTH_UM",
    "SEQUESTRATION_SWEEP",
    "SUBSTRATE_TABLE",
    "immune_area_fraction",
    "immune_sink_rate",
    "immune_sink_strength_for",
    "immune_sink_detail",
    "strength_from_diverted_fraction",
    "diverted_fraction_from_strength",
    "k_max_from_strength",
    "transport_with_immune_sink",
    "substrate_summary",
    "field_summary",
]


# --------------------------------------------------------------------------- #
# Declared geometric constants (NOT fitted to anything)
# --------------------------------------------------------------------------- #
#: Cross-sectional radius (um) attributed to one immune cell when converting a
#: measured centroid DENSITY (cells/um^2) into an AREA FRACTION. A human
#: monocyte-derived macrophage / microglial soma is ~15 um across, the same
#: order as the ~15-20 um glioma soma this repo already assumes for
#: ``TISSUE_DISC_RADIUS_UM``'s justification (``reconstruct/real_xenium.py``).
#: It is a DECLARED MODELLING CONSTANT: it sets the ABSOLUTE scale of
#: ``phi_imm`` and therefore trades off exactly against ``k_max``, which is why
#: ``k_max`` is swept rather than asserted. Only the PRODUCT ``k_max * phi_imm``
#: is physically meaningful, and the sweep spans it.
IMMUNE_CELL_RADIUS_UM = 7.5

#: Gaussian bandwidth (um) of the immune-density kernel — i.e. the length scale
#: over which one macrophage is taken to deplete drug around itself.
#: LITERATURE-ANCHORED: Miller et al. 2015 (PMID 26503691) measured >3-fold
#: elevated payload within "roughly one cell-length (15 um)" of a TAM and
#: resolved the payload concentration gradients around TAMs over length scales
#: of **<25 um**. 25 um is the outer edge of that measured gradient, i.e. the
#: most generous (largest) spatial footprint the measurement supports.
IMMUNE_KERNEL_BANDWIDTH_UM = 25.0

# --------------------------------------------------------------------------- #
# The ONE quantitative literature anchor (a NANOPARTICLE, and a DEPOT)
# --------------------------------------------------------------------------- #
#: Miller 2015 (PMID 26503691): TAMs were 4% of cells in the tumour mass.
MILLER2015_TAM_CELL_FRACTION = 0.04
#: Miller 2015 (PMID 26503691): those TAMs accumulated 30% of total intratumoral
#: therapeutic nanoparticle. This is the largest defensible "fraction of the
#: delivered drug that ends up inside macrophages" number in the literature, and
#: it is for a PLGA-b-PEG nanoparticle, not for a free small molecule.
MILLER2015_TAM_DRUG_FRACTION = 0.30
#: Miller 2015 (PMID 26503691): macrophage depletion lowered intratumoral Pt
#: concentration by 24% (5.0 -> 3.8 uM), an independent estimate of the same
#: diverted fraction.
MILLER2015_DEPLETION_DROP = 0.24
#: Miller 2015 (PMID 26503691): payload concentration gradients around TAMs were
#: resolved over <25 um; >3-fold elevation within ~15 um.
MILLER2015_GRADIENT_LENGTH_UM = 25.0

#: Dimensionless sequestration strengths to sweep. ``S = k_max / k_base`` where
#: ``k_base = drug.clearance + drug.binding`` is the engine's EXISTING scalar
#: first-order loss, so ``S`` reads "how many times the drug's own tissue sink
#: would a fully-macrophage voxel add". ``S = 0`` is OFF and is required to be
#: bit-identical to the unmodified solver. The range is chosen to BRACKET the
#: Miller anchor once it is converted through each section's own measured immune
#: area fraction (:func:`strength_from_diverted_fraction`), which lands at
#: S ~ 5-15 on these sections — see ``results/immune_sink_impact.md`` §3.
SEQUESTRATION_SWEEP: Tuple[float, ...] = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0)

#: Per-drug evidence for TAM sequestration of the FREE (non-encapsulated) drug.
#: ``ratio`` is a measured macrophage:tumour-cell or intracellular:extracellular
#: uptake ratio; **it is None for every drug**, which is the finding. ``sign``
#: records the direction the literature supports, which is not the same thing.
SUBSTRATE_TABLE: Dict[str, Dict] = {
    "doxorubicin": {
        "ratio": None,
        "sign": "sink (mechanistic, unquantified for free drug)",
        "mechanism": (
            "hydrophobic weak base (amino-sugar pKa ~8.3); passive ion-trapping "
            "into acidic organelles, of which macrophages carry the largest "
            "complement. In normal human marrow the nuclei of macrophages and "
            "activated phagocytes were the LEAST doxorubicin-accessible of any "
            "lineage, i.e. the drug is held in the cytoplasmic/vesicular "
            "compartment"
        ),
        "pmids": ["25544758", "26830313", "2295115"],
        "against": (
            "PMID 36105861 measured FREE doxorubicin cellular uptake in TAMs of "
            "a TAM-rich fibrosarcoma as BELOW DETECTABLE LEVELS (only the "
            "liposomal formulations shifted uptake to TAM/MDSC); PMID 9067330 "
            "found a poor correlation between CD11b+ cells and doxorubicin "
            "fluorescence and concluded TAM uptake does not explain tumour drug "
            "accumulation"
        ),
    },
    "cisplatin": {
        "ratio": None,
        "sign": "depot, not sink, and only when nano-encapsulated",
        "mechanism": (
            "a small neutral/aquated platinum species with no basic centre, so "
            "no ion-trapping route. The one quantified TAM result for a "
            "platinum agent is Miller 2015's Pt(IV) PRO-DRUG delivered in a "
            "PLGA-b-PEG nanoparticle, where TAMs acted as a slow-release DEPOT "
            "that raised neighbouring tumour-cell payload; unencapsulated "
            "cisplatin supernatant from the same macrophages at 1-300 uM was "
            "NOT cytotoxic to tumour cells"
        ),
        "pmids": ["26503691"],
        "against": (
            "the effect in PMID 26503691 is a property of the nanoparticle "
            "vehicle and its sign is OPPOSITE to a sink"
        ),
    },
    "temozolomide": {
        "ratio": None,
        "sign": "no evidence in either direction",
        "mechanism": (
            "small (194 Da), neutral, non-basic, freely diffusing, and it is a "
            "PRO-DRUG that decomposes spontaneously at physiological pH rather "
            "than being taken up and held. There is no proposed transporter, no "
            "basic centre for ion-trapping, and no phagocytic route. A TAM sink "
            "for temozolomide would be an invention"
        ),
        "pmids": [],
        "against": "nothing to argue against; there is no claim to test",
    },
}


# --------------------------------------------------------------------------- #
# 1. The FIELD: measured immune cells -> a dimensionless area fraction
# --------------------------------------------------------------------------- #
def immune_area_fraction(
    micro: Microenvironment,
    bandwidth_um: float = IMMUNE_KERNEL_BANDWIDTH_UM,
    cell_radius_um: float = IMMUNE_CELL_RADIUS_UM,
    cell_type: int = int(CellType.IMMUNE),
    positions: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> np.ndarray:
    """Local IMMUNE-cell area fraction in [0, 1] on the microenvironment grid.

    The estimator, stated so it can be checked:

    1. bin the MEASURED immune-cell centroids onto the (ny, nx) grid (counts);
    2. convolve with an isotropic Gaussian of standard deviation
       ``bandwidth_um`` (``gaussian_filter`` conserves the total count, so this
       is a kernel density estimate of the point pattern, not a rescaling);
    3. divide by the voxel area ``dx^2`` to get a NUMBER DENSITY in cells/um^2;
    4. multiply by the per-cell cross-sectional area ``pi * cell_radius_um^2``
       to get a dimensionless AREA FRACTION;
    5. clip to [0, 1] and zero outside the tissue mask.

    Steps 3-4 are the only unit conversion, and it is exact:
    ``(cells/um^2) * (um^2/cell) = 1``.

    Parameters
    ----------
    micro : Microenvironment
        Supplies the grid and (unless ``positions`` is given) the measured cells.
    bandwidth_um : float
        Gaussian bandwidth (um); default :data:`IMMUNE_KERNEL_BANDWIDTH_UM`.
    cell_radius_um : float
        Per-cell cross-sectional radius (um); default
        :data:`IMMUNE_CELL_RADIUS_UM`.
    cell_type : int
        Which ``CellType`` counts as immune. Exposed only so a permutation null
        can be built; the default is the only meaningful value.
    positions : (x, y) arrays, optional
        Override the immune positions entirely — this is how the spatial
        permutation null in ``validation/immune_sink_impact.py`` randomises the
        immune compartment while holding its COUNT fixed.

    Returns
    -------
    (ny, nx) float64 array in [0, 1], zero outside tissue.
    """
    domain = micro.domain
    ny, nx = domain.shape
    dx = float(domain.dx)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)

    if positions is None:
        cells = micro.cells
        sel = np.asarray(cells.cell_type) == int(cell_type)
        x = np.asarray(cells.x, dtype=np.float64)[sel]
        y = np.asarray(cells.y, dtype=np.float64)[sel]
    else:
        x = np.asarray(positions[0], dtype=np.float64)
        y = np.asarray(positions[1], dtype=np.float64)

    counts = domain.zeros()
    if x.size:
        col, row = domain.to_index(x, y)
        np.add.at(counts, (row, col), 1.0)

    sigma = max(float(bandwidth_um) / dx, 1e-6)
    smoothed = gaussian_filter(counts, sigma=sigma, mode="constant")

    number_density = smoothed / (dx * dx)                    # cells / um^2
    area_per_cell = np.pi * float(cell_radius_um) ** 2       # um^2 / cell
    phi = number_density * area_per_cell                     # dimensionless

    phi = np.clip(phi, 0.0, 1.0)
    phi[~tissue] = 0.0
    return phi


def immune_sink_rate(
    micro: Microenvironment,
    k_max_per_s: float,
    bandwidth_um: float = IMMUNE_KERNEL_BANDWIDTH_UM,
    cell_radius_um: float = IMMUNE_CELL_RADIUS_UM,
    phi: Optional[np.ndarray] = None,
    **phi_kwargs,
) -> np.ndarray:
    """The per-voxel first-order loss rate field ``k_sink = k_max * phi_imm``.

    Units: ``k_max_per_s`` is 1/s, ``phi_imm`` is dimensionless, so the result is
    1/s and is directly addable to the engine's scalar
    ``drug.clearance + drug.binding``.

    ``k_max_per_s = 0`` returns an exact all-zero array, which is what makes the
    OFF path bit-identical rather than merely close.
    """
    if phi is None:
        phi = immune_area_fraction(
            micro, bandwidth_um=bandwidth_um, cell_radius_um=cell_radius_um,
            **phi_kwargs,
        )
    return float(k_max_per_s) * np.asarray(phi, dtype=np.float64)


# --------------------------------------------------------------------------- #
# 2. The MAGNITUDE: swept, never fitted
# --------------------------------------------------------------------------- #
def k_max_from_strength(drug: DrugParams, strength: float) -> float:
    """``k_max = S * (drug.clearance + drug.binding)`` in 1/s.

    ``S`` is dimensionless and reads "a voxel packed entirely with macrophages
    removes drug ``S`` times as fast as the drug's own tissue sink already
    does". Expressing the magnitude relative to a term the model ALREADY has
    means the sweep is comparable across drugs whose absolute clearances differ
    by orders of magnitude.
    """
    return float(strength) * (float(drug.clearance) + float(drug.binding))


def strength_from_diverted_fraction(
    diverted_fraction: float, mean_phi: float
) -> float:
    """Convert a measured "fraction of tissue drug held by macrophages" into ``S``.

    At quasi-steady state two first-order losses in parallel split the drug in
    proportion to their rate constants, so if a fraction ``f`` of the disposition
    goes through the immune sink,

        f = <k_sink> / (<k_sink> + k_base),   <k_sink> = k_max * <phi_imm>

    =>  S = k_max / k_base = f / ((1 - f) * <phi_imm>).

    This is how the Miller 2015 anchor (f = 0.30 from the 4%-of-cells/30%-of-drug
    partition, or f = 0.24 from the depletion experiment) is mapped onto THIS
    section's own measured mean immune area fraction, rather than being imported
    as a number. It is a bracketing calculation, not a calibration: the anchor is
    a nanoparticle in a murine fibrosarcoma and the drugs here are free small
    molecules in human glioma.
    """
    f = float(diverted_fraction)
    phi = float(mean_phi)
    if not (0.0 <= f < 1.0):
        raise ValueError(f"diverted_fraction must be in [0,1), got {f}")
    if phi <= 0.0:
        return 0.0
    return f / ((1.0 - f) * phi)


def diverted_fraction_from_strength(strength: float, mean_phi: float) -> float:
    """Inverse of :func:`strength_from_diverted_fraction`: ``f = S phi/(1 + S phi)``."""
    s_phi = float(strength) * float(mean_phi)
    return s_phi / (1.0 + s_phi) if s_phi > 0.0 else 0.0


def immune_sink_strength_for(drug: str) -> float:
    """The MEASURED sequestration strength for a drug — **exactly 0.0, always**.

    This is the module's identity, in the sense
    ``molecular/uptake.py::uptake_factor`` returns exactly 1.0 for everything
    unproven. No entry in :data:`SUBSTRATE_TABLE` carries a measured free-drug
    macrophage:tumour uptake ratio, so no drug gets a non-zero sink from the
    literature and the whole module is inert unless a strength is passed in
    explicitly by a sweep. Call :func:`immune_sink_detail` to see WHY for a
    given drug — the three reasons are not equivalent.
    """
    return 0.0


def immune_sink_detail(drug: str) -> Dict:
    """Provenance for :func:`immune_sink_strength_for` — which reason applies."""
    name = str(drug).strip().lower()
    row = SUBSTRATE_TABLE.get(name)
    if row is None:
        return {
            "drug": name,
            "strength": 0.0,
            "applied": False,
            "ratio": None,
            "reason": "not in SUBSTRATE_TABLE; no evidence was assessed",
            "sign": "unknown",
            "source": None,
            "pmids": [],
        }
    return {
        "drug": name,
        "strength": 0.0,
        "applied": False,
        "ratio": row["ratio"],
        "reason": (
            "no measured FREE-drug macrophage:tumour-cell uptake ratio exists; "
            "the literature supports a direction only, so the magnitude is a "
            "SWEPT parameter (SEQUESTRATION_SWEEP) and never a fitted one"
        ),
        "sign": row["sign"],
        "mechanism": row["mechanism"],
        "against": row["against"],
        "source": "literature review, see module docstring",
        "pmids": list(row["pmids"]),
    }


# --------------------------------------------------------------------------- #
# 3. The SOLVER: physics.diffusion.simulate_drug_transport with a per-voxel k
# --------------------------------------------------------------------------- #
def transport_with_immune_sink(
    micro: Microenvironment,
    drug: DrugParams,
    source: np.ndarray,
    pde: PDEConfig,
    velocity: Optional[np.ndarray] = None,
    source_profile: Optional[np.ndarray] = None,
    sink_rate: Optional[np.ndarray] = None,
) -> TransportResult:
    """``simulate_drug_transport`` with an EXTRA per-voxel first-order loss.

    WHY A RE-SOLVE AND NOT AN AUC RESCALE. The engine's solver builds its loss
    term as ``k = np.full((ny, nx), drug.clearance + drug.binding)`` — a SCALAR
    broadcast to the grid — and exposes no per-voxel hook. A spatially varying
    sink cannot be recovered by post-multiplying the AUC field, because drug
    diffuses INTO the depleted region from its neighbours: the depletion is
    non-local and that non-locality is precisely the microenvironment effect
    being measured. So the transport is re-solved.

    HOW IT AVOIDS DRIFT. ``physics/diffusion.py`` is NOT modified and NOT copied:
    its three numerical primitives (``_effective_diffusivity``,
    ``_build_diffusion_operator``, ``_advective_divergence``) are imported and
    used unchanged, exactly as ``validation/moncayo_explant.py::bath_diffusion``
    already does. The single functional difference from
    ``simulate_drug_transport`` is one line::

        k = np.full((ny, nx), clearance + binding) + sink_rate

    and since ``x + 0.0 == x`` exactly in IEEE-754 for every finite ``x``, a
    ``sink_rate`` of ``None`` or of exact zeros reproduces the shipped solver
    BIT-FOR-BIT. That is asserted numerically, not assumed
    (``tests/test_immune_sink.py``, ``validation/immune_sink_impact.py`` §1).

    Parameters
    ----------
    micro, drug, source, pde, velocity, source_profile
        Exactly as ``physics.diffusion.simulate_drug_transport``.
    sink_rate : (ny, nx) float array, optional
        Extra first-order loss in 1/s, typically from :func:`immune_sink_rate`.
        ``None`` (default) means OFF.

    Returns
    -------
    TransportResult, with the same fields the shipped solver returns.
    """
    domain = micro.domain
    ny, nx = domain.shape
    dx = float(domain.dx)
    tissue = np.asarray(domain.tissue_mask, dtype=bool)

    source = np.asarray(source, dtype=np.float64)
    if source.shape != (ny, nx):
        raise ValueError(f"source shape {source.shape} != grid {(ny, nx)}")
    if pde.dt <= 0.0:
        raise ValueError(f"pde.dt must be > 0, got {pde.dt}")
    if pde.total_time <= 0.0:
        raise ValueError(f"pde.total_time must be > 0, got {pde.total_time}")

    dt = float(pde.dt)

    d_eff = _effective_diffusivity(drug, micro.fields.ecm_density)
    k = np.full((ny, nx), float(drug.clearance) + float(drug.binding))
    if sink_rate is not None:
        sink = np.asarray(sink_rate, dtype=np.float64)
        if sink.shape != (ny, nx):
            raise ValueError(f"sink_rate shape {sink.shape} != grid {(ny, nx)}")
        if np.any(sink < 0.0):
            raise ValueError("sink_rate must be non-negative (it is a LOSS rate)")
        if not np.all(np.isfinite(sink)):
            raise ValueError("sink_rate contains non-finite values")
        k = k + sink
    k[~tissue] = 0.0
    src = np.where(tissue, source, 0.0)

    use_convection = bool(pde.convection)
    if velocity is None:
        velocity = getattr(micro.fields, "velocity", None)
    if velocity is None:
        use_convection = False
    if use_convection:
        velocity = np.asarray(velocity, dtype=np.float64)
        if velocity.shape != (2, ny, nx):
            raise ValueError(f"velocity shape {velocity.shape} != {(2, ny, nx)}")
        vx, vy = velocity[0], velocity[1]
        max_speed = float(np.max(np.abs(vx) + np.abs(vy))) if velocity.size else 0.0
        cfl = max_speed * dt / dx
        if cfl > 1.0:
            raise ValueError(
                f"Advection CFL violated: dt*(|vx|+|vy|)/dx = {cfl:.3f} > 1. "
                f"Reduce pde.dt below {dx / max_speed:.4g} s or disable convection."
            )
    else:
        vx = vy = None

    n = ny * nx
    lap = _build_diffusion_operator(d_eff, tissue, dx)
    a = csc_matrix((np.ones(n), (np.arange(n), np.arange(n))), shape=(n, n))
    a = a - dt * lap
    a = a + csc_matrix(
        (dt * k.ravel(), (np.arange(n), np.arange(n))), shape=(n, n)
    )
    a = a.tocsc()
    solve = factorized(a)

    n_steps = int(round(pde.total_time / dt))
    save_every = max(1, int(pde.save_every))

    if source_profile is None:
        src_scale = None
    else:
        src_scale = np.asarray(source_profile, dtype=np.float64)
        if src_scale.ndim != 1:
            raise ValueError(
                f"source_profile must be 1-D, got shape {src_scale.shape}"
            )
        if src_scale.size < n_steps:
            raise ValueError(
                f"source_profile has {src_scale.size} samples but the solve takes "
                f"{n_steps} steps; supply at least one scale factor per step."
            )
        used = src_scale[:n_steps]
        if not np.all(np.isfinite(used)):
            raise ValueError("source_profile contains non-finite values")
        if float(used.min()) < 0.0:
            raise ValueError(
                f"source_profile must be non-negative (min {float(used.min()):.4g})"
            )

    c = domain.zeros()
    auc = domain.zeros()

    times = [0.0]
    mean_hist = [float(c[tissue].mean()) if tissue.any() else 0.0]

    for step in range(1, n_steps + 1):
        rhs = c.copy()
        if src_scale is None:
            rhs += dt * src
        else:
            rhs += (dt * src_scale[step - 1]) * src
        if use_convection:
            rhs -= dt * _advective_divergence(c, vx, vy, tissue, dx)
        c_new = solve(rhs.ravel()).reshape(ny, nx)
        np.maximum(c_new, 0.0, out=c_new)
        c_new[~tissue] = 0.0

        auc += 0.5 * (c + c_new) * dt
        c = c_new

        if step % save_every == 0 or step == n_steps:
            times.append(step * dt)
            mean_hist.append(float(c[tissue].mean()) if tissue.any() else 0.0)

    final_concentration = c

    if tissue.any():
        auc_t = auc[tissue]
        peak = float(np.percentile(auc_t, 99.0))
        penetration_index = (
            float(np.clip(auc_t.mean() / peak, 0.0, 1.0)) if peak > 0.0 else 0.0
        )
    else:
        penetration_index = 0.0

    return TransportResult(
        auc=auc,
        final_concentration=final_concentration,
        times=np.asarray(times, dtype=np.float64),
        mean_concentration_history=np.asarray(mean_hist, dtype=np.float64),
        penetration_index=penetration_index,
    )


# --------------------------------------------------------------------------- #
# 4. Reporting helpers
# --------------------------------------------------------------------------- #
def substrate_summary() -> Dict[str, int]:
    """How many drugs carry a MEASURED free-drug TAM uptake ratio. (Zero.)"""
    n = len(SUBSTRATE_TABLE)
    measured = sum(1 for r in SUBSTRATE_TABLE.values() if r["ratio"] is not None)
    return {
        "n_drugs_assessed": n,
        "n_with_measured_ratio": measured,
        "n_direction_only": sum(
            1 for r in SUBSTRATE_TABLE.values()
            if r["ratio"] is None and "no evidence" not in r["sign"]
        ),
        "n_no_evidence": sum(
            1 for r in SUBSTRATE_TABLE.values() if "no evidence" in r["sign"]
        ),
    }


def field_summary(phi: np.ndarray, tissue: np.ndarray) -> Dict[str, float]:
    """Descriptive statistics of an immune area-fraction field over tissue."""
    v = np.asarray(phi, dtype=np.float64)[np.asarray(tissue, dtype=bool)]
    if v.size == 0:
        return {"mean": 0.0, "median": 0.0, "p90": 0.0, "p99": 0.0, "max": 0.0}
    return {
        "mean": float(v.mean()),
        "median": float(np.median(v)),
        "p90": float(np.percentile(v, 90.0)),
        "p99": float(np.percentile(v, 99.0)),
        "max": float(v.max()),
    }


if __name__ == "__main__":  # pragma: no cover - self-test
    print("=== tmesim.molecular.immune_sink self-test ===")
    print("substrate summary:", substrate_summary())
    for d in SUBSTRATE_TABLE:
        det = immune_sink_detail(d)
        print(f"  {d:<14} strength={det['strength']:.1f}  sign={det['sign']}")
    for f in (MILLER2015_DEPLETION_DROP, MILLER2015_TAM_DRUG_FRACTION):
        for phi in (0.02, 0.05, 0.10):
            s = strength_from_diverted_fraction(f, phi)
            back = diverted_fraction_from_strength(s, phi)
            print(f"  f={f:.2f} phi={phi:.2f} -> S={s:7.2f} (round-trip f={back:.4f})")
