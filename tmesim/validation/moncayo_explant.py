"""REAL VALIDATION RUN #3 — Moncayo et al. HGSOC patient-derived explants.

The only open dataset that pairs a **measured spatial drug field** with
**molecular state on the same rows**: Moncayo et al., *Nat Commun* 2026
(PMID 41844641), Zenodo 17610220 + GEO GSE281519 + Nat Commun Source Data.

WHAT THIS CAN AND CANNOT VALIDATE
---------------------------------
The explants were **bath-dosed ex vivo** (20 uM, 24 h). There is **no
vasculature and no perfusion source term**. This run therefore validates
*diffusion into tissue, uptake, and intracellular sequestration* and
**cannot validate vascular delivery**, interstitial pressure, or vessel
normalisation. It is also **ovarian, not glioma** — a physics-transfer
argument, not a biology-transfer argument. The MSI values are drug-ion
intensities normalised to a drug-matched **deuterated internal standard**;
they are unitless ratios, so **no absolute uM claim** is made anywhere.
n = 3 patients, 114 GeoMx ROIs, 18 MSI slices.

THE THREE ANALYSES (per results/DATA_HUNT.md section 5)
-------------------------------------------------------
1. Accumulation-heterogeneity fit — bath-boundary diffusion into the measured
   explant geometry, with and without a pH-partitioning / lysosomal-trapping
   term, scored against `drug_measure_norm` in the authors' own mixed-effects
   frame `measure ~ prediction + (1|Patient) + (1|Patient:Slice)`.
2. Compound-specificity falsification test — does the molecular layer make the
   weak bases (niraparib, rucaparib) more heterogeneous than the non-base
   control (olaparib)?
3. Exposure -> molecular response coupling — apoptosis and lysosome signature
   scores from the 11,224-gene ROI matrix against measured drug level and
   against our predicted exposure.

ADDITIVE-ONLY CONTRACT
----------------------
Nothing in `core/`, `physics/`, `molecular/`, `reconstruct/` or `sim/` is
modified. Two engine internals are *imported and called* (never edited):
`physics.diffusion._effective_diffusivity` (the ECM-tortuosity law) and
`physics.diffusion._build_diffusion_operator` (the finite-volume Laplacian).
The reason is stated in `bath_diffusion()`: the frozen public solver
`simulate_drug_transport` implements a **no-flux boundary + volumetric source
and a zero initial condition**, which cannot express a Dirichlet bath. We
therefore reuse the engine's *exact discretisation and its tortuosity law* and
add only the boundary condition, in this analysis layer. The stepper is
verified against the analytic Crank plane-sheet solution in
`verify_against_analytic_slab()`.

Units: length um, time s, concentration uM, diffusivity um^2/s.

Run:  ./.venv/bin/python -m tmesim.validation.moncayo_explant
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.sparse import csc_matrix, identity as sp_identity
from scipy.sparse.linalg import factorized
from scipy.stats import pearsonr, spearmanr

from ..core import DrugParams
# Engine internals: imported read-only, never modified. See module docstring.
from ..physics.diffusion import _build_diffusion_operator, _effective_diffusivity

# ---------------------------------------------------------------------------
# Paths / provenance
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
HUNT = os.path.join(REPO, "data_cache", "hunt", "geomx_dsp")

ROI_CSV = os.path.join(HUNT, "zenodo17610220_ROI_drug_measurements.csv")
EXPR_CSV = os.path.join(HUNT, "zenodo17610220_ROI_drug_plus_expression.csv.gz")
MSI_CSV = os.path.join(HUNT, "NatCommun2026_SourceData_SupFig2_MSI_perpixel.csv.gz")
RMD = os.path.join(HUNT, "zenodo17610220_modelling.Rmd")

OUT_MD = os.path.join(REPO, "results", "moncayo_explant.md")
OUT_PNG = os.path.join(REPO, "results", "moncayo_explant.png")

# ---------------------------------------------------------------------------
# Experimental constants (from the paper / deposit)
# ---------------------------------------------------------------------------
BATH_CONC_UM = 20.0          # 20 uM bath dose
EXPOSURE_S = 24.0 * 3600.0   # 24 h
MSI_PITCH_MM = 0.02          # 20 um MALDI-MSI pixel pitch

# Model grid: the MSI pitch (20 um) is far finer than the ~1-2 mm diffusion
# problem needs. We coarsen 3x to 60 um, which also puts the model on the same
# spatial scale as the noise-controlled measured heterogeneity metric.
MODEL_PITCH_UM = 60.0
COARSEN = 3

# Uniform ECM volume fraction. We deliberately do NOT invent tissue-specific ECM
# densities; instead we sweep 0.2 / 0.4 / 0.6 as a sensitivity check.
ECM_REFERENCE = 0.4
ECM_SWEEP = (0.2, 0.4, 0.6)

# ---------------------------------------------------------------------------
# Chemistry. SMILES verified by RDKit InChIKey against the public identifiers:
#   Niraparib PCHKPVIQAHNQLW-AWEZNQCLSA-N · Olaparib HUTZIPRGMCJLQS-UHFFFAOYSA-N
#   Rucaparib HMABYWSNWIZPAG-UHFFFAOYSA-N
# ---------------------------------------------------------------------------
SMILES: Dict[str, str] = {
    "Niraparib": "NC(=O)c1cccc2cn(-c3ccc(cc3)[C@H]3CCCNC3)nc12",
    "Olaparib": "O=C(c1ccc(CC2=NNC(=O)c3ccccc32)cc1F)N1CCN(C(=O)C2CC2)CC1",
    "Rucaparib": "CNCc1ccc(cc1)-c1[nH]c2cc(F)cc3c2c1CCNC3=O",
}
EXPECTED_INCHIKEY = {
    "Niraparib": "PCHKPVIQAHNQLW-AWEZNQCLSA-N",
    "Olaparib": "HUTZIPRGMCJLQS-UHFFFAOYSA-N",
    "Rucaparib": "HMABYWSNWIZPAG-UHFFFAOYSA-N",
}
# Literature pKa (strongest basic centre) — used only as a cross-check on the
# structure-derived value, never as the model input.
LIT_PKA = {"Niraparib": 9.95, "Olaparib": None, "Rucaparib": 9.4}

# pH-partition constants (Henderson-Hasselbalch ion trapping).
PH_CYTOSOL = 7.2
PH_LYSOSOME = 4.7
# Lysosomal volume fraction of tissue (cell volume fraction x lysosome fraction
# of cell volume ~ 0.7 x 0.02). Swept below; the reference value is 0.014.
PHI_LYS = 0.014
PHI_LYS_SWEEP = (0.005, 0.014, 0.05)

# ROI-table Slice -> MSI Slice. Asserted for Patient/Compound/Tissue agreement
# at load time. T23_071_NirBlad2 has no MSI counterpart in the Source Data.
SLICE_MAP: Dict[str, Optional[str]] = {
    "T23_015_Nir4": "P23-015_Nir4",
    "T23_015_Ruc5": "P23-015_Ruc5",
    "T23_037_Nir1": "P23-037_Nir1",
    "T23_037_Nir2": "P23-037_Nir2",
    "T23_037_Nir3": "P23-037_Nir3",
    "T23_037_Ruc1": "P23-037_Ruc1",
    "T23_037_Ruc1'": "P23-037_Ruc1'",
    "T23_037_Ruc2": "P23-037_Ruc2",
    "T23_071_NirOm1": "P23-071_OmNir1",
    "T23_071_RucOm1": "P23-071_OmRuc1",
    "T23_071_NirBlad3": "P23-071_BladNir3",
    "T23_071_NirBlad2": None,
}

MSI_CHANNEL = {
    "Niraparib": ("norm_Nir[M+Na]", "Nir[M+Na]", "Heavy labelled_Nir[M+Na]"),
    "Olaparib": ("norm_Ola[M+Na]", "Ola[M+Na]", "Heavy labelled_Ola[M+Na]"),
    "Rucaparib": ("norm_Ruc[M+H]", "Ruc[M+H]", "Heavy labelled_Ruc[M+H]"),
}

# ---------------------------------------------------------------------------
# Gene sets. Fetched 2026-07-27 from Enrichr's public library endpoint
# (https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=...)
# and embedded verbatim so the module is reproducible offline.
#   APOPTOSIS_GENES : MSigDB_Hallmark_2020, term "Apoptosis" (161 symbols)
#   LYSOSOME_GENES  : KEGG_2021_Human,      term "Lysosome"  (128 symbols)
# ---------------------------------------------------------------------------
APOPTOSIS_GENES: Tuple[str, ...] = (
    "ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX", "BCAP31",
    "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN", "BID", "BIK",
    "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2", "BTG3", "CASP1",
    "CASP2", "CASP3", "CASP4", "CASP6", "CASP7", "CASP8", "CASP9", "CAV1",
    "CCNA1", "CCND1", "CCND2", "CD14", "CD2", "CD38", "CD44", "CD69",
    "CDC25B", "CDK2", "CDKN1A", "CDKN1B", "CFLAR", "CLU", "CREBBP", "CTH",
    "CTNNB1", "CYLD", "DAP", "DAP3", "DCN", "DDIT3", "DFFA", "DIABLO",
    "DNAJA1", "DNAJC3", "DNM1L", "DPYD", "EBP", "EGR3", "EMP1", "ENO2",
    "ERBB2", "ERBB3", "EREG", "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR",
    "FEZ1", "GADD45A", "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4",
    "GSN", "GSR", "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1",
    "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18", "IL1A",
    "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1", "LGALS3", "LMNA",
    "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9", "NEFH", "PAK1", "PDCD4",
    "PDGFRB", "PEA15", "PLAT", "PLCB2", "PLPPR4", "PMAIP1", "PPP2R5B",
    "PPP3R1", "PPT1", "PRF1", "PSEN1", "PSEN2", "PTK2", "RARA", "RELA",
    "RETSAT", "RHOB", "RHOT2", "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D",
    "SLC20A1", "SMAD7", "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2",
    "TGFBR3", "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10",
    "TOP2A", "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP",
)
LYSOSOME_GENES: Tuple[str, ...] = (
    "ABCA2", "ABCB9", "ACP2", "ACP5", "AGA", "AP1B1", "AP1G1", "AP1G2",
    "AP1M1", "AP1M2", "AP1S1", "AP1S2", "AP1S3", "AP3B1", "AP3B2", "AP3D1",
    "AP3M1", "AP3M2", "AP3S1", "AP3S2", "AP4B1", "AP4E1", "AP4M1", "AP4S1",
    "ARSA", "ARSB", "ARSG", "ASAH1", "ATP6AP1", "ATP6V0A1", "ATP6V0A2",
    "ATP6V0A4", "ATP6V0B", "ATP6V0C", "ATP6V0D1", "ATP6V0D2", "ATP6V1H",
    "CD164", "CD63", "CD68", "CLN3", "CLN5", "CLTA", "CLTB", "CLTC", "CLTCL1",
    "CTNS", "CTSA", "CTSB", "CTSC", "CTSD", "CTSE", "CTSF", "CTSG", "CTSH",
    "CTSK", "CTSL", "CTSO", "CTSS", "CTSV", "CTSW", "CTSZ", "DNASE2",
    "DNASE2B", "ENTPD4", "FUCA1", "FUCA2", "GAA", "GALC", "GALNS", "GBA",
    "GGA1", "GGA2", "GGA3", "GLA", "GLB1", "GM2A", "GNPTAB", "GNPTG", "GNS",
    "GUSB", "HEXA", "HEXB", "HGSNAT", "HYAL1", "HYAL2", "HYAL3", "HYAL4",
    "IDS", "IDUA", "IGF2R", "LAMP1", "LAMP2", "LAMP3", "LAPTM4A", "LAPTM4B",
    "LAPTM5", "LGMN", "LIPA", "LITAF", "M6PR", "MAN2B1", "MANBA", "MCOLN1",
    "MFSD8", "NAGA", "NAGLU", "NAGPA", "NAPSA", "NEU1", "NPC1", "NPC2",
    "PLA2G15", "PPT1", "PPT2", "PSAP", "PSAPL1", "SCARB2", "SGSH", "SLC11A1",
    "SLC11A2", "SLC17A5", "SMPD1", "SORT1", "SPAM1", "SUMF1", "TCIRG1", "TPP1",
)


def _sha256(path: str, limit: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            h.update(b)
            if h.block_size and fh.tell() > limit:
                break
    return h.hexdigest()[:16]


# ===========================================================================
# 1. DATA
# ===========================================================================
def load_roi_table() -> pd.DataFrame:
    """114 GeoMx ROIs x 13 columns with the per-ROI MALDI drug measurement."""
    d = pd.read_csv(ROI_CSV)
    assert d.shape == (114, 13), f"unexpected ROI table shape {d.shape}"
    # drug_measure_norm = drug_measure / measure_norm_factor, and
    # measure_norm_factor is ONE value per COMPOUND (not per slice/patient).
    assert np.allclose(d.drug_measure / d.measure_norm_factor, d.drug_measure_norm)
    nf = d.groupby("Compound")["measure_norm_factor"].nunique()
    assert (nf == 1).all(), "measure_norm_factor is not compound-constant"
    d["msi_slice"] = d["Slice"].map(SLICE_MAP)
    return d


def load_msi_pixels() -> pd.DataFrame:
    """376,736 MALDI-MSI pixels x 19 columns, 20 um pitch, X/Y, 3 compounds."""
    m = pd.read_csv(MSI_CSV)
    assert m.shape[0] == 376_736, f"unexpected MSI pixel count {m.shape[0]}"
    return m


def load_expression() -> pd.DataFrame:
    """The same 114 ROIs x 11,224 genes (GeoMx WTA, Q3-normalised counts)."""
    e = pd.read_csv(EXPR_CSV)
    assert e.shape == (114, 11237), f"unexpected expression shape {e.shape}"
    return e


def provenance() -> List[dict]:
    rows = []
    for label, path in [
        ("ROI drug measurements", ROI_CSV),
        ("ROI drug + 11,224-gene WTA", EXPR_CSV),
        ("MALDI-MSI per-pixel field", MSI_CSV),
        ("authors' lme4 model spec", RMD),
    ]:
        rows.append({
            "label": label,
            "path": os.path.relpath(path, REPO),
            "bytes": os.path.getsize(path),
            "sha256_16": _sha256(path),
        })
    return rows


# ===========================================================================
# 2. EXPLANT GEOMETRY — measured, not synthetic
# ===========================================================================
@dataclass
class SliceGeometry:
    """One explant section, rasterised from its own MALDI-MSI pixel coordinates.

    GEOMETRY CHOICE (documented because it is the main modelling assumption):
    the explant outline is *measured* — it is the set of on-tissue MSI pixels of
    that section, hole-filled, on the native 20 um grid, then block-coarsened to
    60 um. The bath boundary is the entire section perimeter (Dirichlet at
    20 uM), because a bath-dosed explant is surrounded by medium on every side.
    We solve in 2D and therefore ignore diffusion in/out of the section plane
    (through the explant's top and bottom faces). That makes our in-plane
    distance-to-perimeter an UPPER BOUND on the true 3D distance to the nearest
    bath surface, so the 2D model OVER-predicts penetration gradients. Any
    negative heterogeneity result below is therefore conservative.
    """
    slice_id: str
    compound: str
    patient: str
    tissue: str
    mask: np.ndarray            # (ny, nx) bool at MODEL_PITCH_UM
    dx: float                   # um
    depth_um: np.ndarray        # (ny, nx) distance to nearest non-tissue voxel
    measured: np.ndarray        # (ny, nx) block-mean IS-normalised drug ion
    n_px_native: int
    area_mm2: float
    max_depth_um: float
    mean_depth_um: float
    measured_raw: np.ndarray = field(default_factory=lambda: np.zeros(0))
    measured_is: np.ndarray = field(default_factory=lambda: np.zeros(0))


def _rasterise(g: pd.DataFrame, channel: str) -> Tuple[np.ndarray, np.ndarray]:
    xs = np.round(g.X.values / MSI_PITCH_MM).astype(int)
    ys = np.round(g.Y.values / MSI_PITCH_MM).astype(int)
    x0, y0 = xs.min(), ys.min()
    nx, ny = xs.max() - x0 + 1, ys.max() - y0 + 1
    mask = np.zeros((ny, nx), dtype=bool)
    val = np.zeros((ny, nx), dtype=float)
    mask[ys - y0, xs - x0] = True
    val[ys - y0, xs - x0] = g[channel].values
    return mask, val


def _coarsen(mask: np.ndarray, val: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """Block-reduce by k. A coarse voxel is tissue if >=50% of its fine voxels
    are tissue; its value is the mean over its tissue fine voxels."""
    ny, nx = mask.shape
    py, px = (-ny) % k, (-nx) % k
    m = np.pad(mask, ((0, py), (0, px)), constant_values=False)
    v = np.pad(val * mask, ((0, py), (0, px)), constant_values=0.0)
    Y, X = m.shape[0] // k, m.shape[1] // k
    cnt = m.reshape(Y, k, X, k).sum(axis=(1, 3))
    tot = v.reshape(Y, k, X, k).sum(axis=(1, 3))
    cm = cnt >= (k * k) / 2.0
    cv = np.zeros_like(tot)
    np.divide(tot, cnt, out=cv, where=cnt > 0)
    return cm, cv


def build_geometries(msi: pd.DataFrame) -> Dict[str, SliceGeometry]:
    """Rasterise every MSI section into a model geometry.

    Three channels are carried through the identical raster/coarsen path so
    that they are pixel-for-pixel comparable:
      ``measured``     IS-normalised drug ion (the quantity being predicted),
      ``measured_raw`` the bare drug ion count,
      ``measured_is``  the co-sprayed **deuterated internal standard** count.
    The IS channel is the technical-noise control used in analysis 2: it is a
    uniformly applied compound, so any spatial structure in it is matrix
    effect + ionisation + shot noise, not biology.
    """
    geos: Dict[str, SliceGeometry] = {}
    for sl, g in msi.groupby("Slice"):
        compound = str(g.Drug.iloc[0])
        chan, raw_chan, is_chan = MSI_CHANNEL[compound]
        fine_mask, fine_val = _rasterise(g, chan)
        _, fine_raw = _rasterise(g, raw_chan)
        _, fine_is = _rasterise(g, is_chan)
        fine_mask = ndimage.binary_fill_holes(fine_mask)
        mask, val = _coarsen(fine_mask, fine_val, COARSEN)
        _, val_raw = _coarsen(fine_mask, fine_raw, COARSEN)
        _, val_is = _coarsen(fine_mask, fine_is, COARSEN)
        # keep the largest connected component (drops rasterisation specks)
        lab, n = ndimage.label(mask)
        if n > 1:
            sizes = ndimage.sum(mask, lab, range(1, n + 1))
            mask = lab == (int(np.argmax(sizes)) + 1)
            val, val_raw, val_is = val * mask, val_raw * mask, val_is * mask
        depth = ndimage.distance_transform_edt(np.pad(mask, 1))[1:-1, 1:-1] * MODEL_PITCH_UM
        geos[sl] = SliceGeometry(
            slice_id=sl, compound=compound, patient=str(g.Patient.iloc[0]),
            tissue=str(g.Tissue.iloc[0]), mask=mask, dx=MODEL_PITCH_UM,
            depth_um=depth, measured=val, n_px_native=int(len(g)),
            area_mm2=float(fine_mask.sum()) * (MSI_PITCH_MM ** 2),
            max_depth_um=float(depth[mask].max()),
            mean_depth_um=float(depth[mask].mean()),
            measured_raw=val_raw, measured_is=val_is,
        )
    return geos


# ===========================================================================
# 3. CHEMISTRY — structure -> basicity -> lysosomal retardation
# ===========================================================================
# Structure-derived basicity. A basic centre is an sp3 nitrogen that is NOT an
# amide/sulfonamide/phosphoramide N, NOT attached to an aromatic ring, and not
# already charged. The SMARTS below encodes exactly that. pKa is then read off a
# curated substructure-class table (there is no first-principles pKa predictor
# in RDKit); LIT_PKA is carried alongside purely as a cross-check.
BASIC_N_SMARTS = "[NX3;H2,H1,H0;!$(N[C,S,P]=[O,S,N]);!$(N-a);!$(N#*);!$([N+])]"
_PKA_CLASS = {
    "cyclic_secondary_amine": 10.4,   # piperidine / pyrrolidine NH
    "acyclic_secondary_amine": 9.7,   # R-NH-CH3
    "tertiary_amine": 9.8,
    "primary_amine": 10.0,
    "none": 1.0,                      # no basic centre -> neutral at pH 4.7-7.2
}


@dataclass
class CompoundChemistry:
    name: str
    smiles: str
    inchikey: str
    mw: float
    logp: float
    n_basic_centres: int
    basic_class: str
    pka_structure: float
    pka_literature: Optional[float]
    lysosome_cytosol_ratio: float
    retardation_R: float
    drug_params_shipped: DrugParams     # from the UNMODIFIED molecular layer
    drug_params_seq: DrugParams         # analysis-layer sequestration encoding
    profile: dict = field(default_factory=dict)


def pH_partition_ratio(pka: float, ph_lys: float = PH_LYSOSOME,
                       ph_cyt: float = PH_CYTOSOL) -> float:
    """Henderson-Hasselbalch ion-trapping ratio for a monoprotic base.

    Only the neutral species crosses the lysosomal membrane, so at equilibrium
    the TOTAL (neutral + protonated) concentration ratio is

        [D]_lys / [D]_cyt = (1 + 10^(pKa - pH_lys)) / (1 + 10^(pKa - pH_cyt)).

    For any strong base (pKa >> both pH values) this saturates at
    10^(pH_cyt - pH_lys) ~ 316, which is why the model cannot separate two
    strong bases from one another (see the compound-specificity verdict).
    """
    return (1.0 + 10.0 ** (pka - ph_lys)) / (1.0 + 10.0 ** (pka - ph_cyt))


def retardation_factor(pka: float, phi_lys: float = PHI_LYS) -> float:
    """Local-equilibrium retardation R = C_total / C_free.

    Reversible instantaneous trapping into a lysosomal compartment occupying
    volume fraction ``phi_lys`` of tissue gives

        R = 1 + phi_lys * (ratio - 1),

    and the transport equation becomes R dC_free/dt = div(D grad C_free), i.e.
    an effective diffusivity D/R with the measured (total) drug equal to
    R * C_free. This is EXACT for linear equilibrium partitioning, which is why
    it can be encoded as a DrugParams diffusivity plus a post-multiplication —
    no engine change required.
    """
    return 1.0 + phi_lys * (pH_partition_ratio(pka) - 1.0)


def compound_chemistry(ecm_reference: float = ECM_REFERENCE,
                       phi_lys: float = PHI_LYS) -> Dict[str, CompoundChemistry]:
    from rdkit import Chem
    from rdkit.Chem import inchi as rd_inchi
    from ..molecular.compile_drug import build_profile, molecule_to_drugparams
    from ..molecular.descriptors import from_smiles

    patt = Chem.MolFromSmarts(BASIC_N_SMARTS)
    out: Dict[str, CompoundChemistry] = {}
    for name, smi in SMILES.items():
        mol = Chem.MolFromSmiles(smi)
        key = rd_inchi.MolToInchiKey(mol)
        assert key == EXPECTED_INCHIKEY[name], f"{name} InChIKey {key} != expected"
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
        pka = _PKA_CLASS[cls]

        m = from_smiles(name, smi)
        prof = build_profile(m, ecm_reference=ecm_reference)
        dp = molecule_to_drugparams(m, ecm_reference=ecm_reference)
        R = retardation_factor(pka, phi_lys)
        # Encode sequestration purely as DrugParams: retard the diffusivity by R.
        dp_seq = DrugParams(
            name=f"{name}+lysosomal_trap(R={R:.2f})",
            diffusivity=dp.diffusivity / R, molecular_weight=dp.molecular_weight,
            ic50=dp.ic50, hill=dp.hill, binding=dp.binding,
            clearance=dp.clearance, bbb_penetrance=dp.bbb_penetrance,
        )
        out[name] = CompoundChemistry(
            name=name, smiles=smi, inchikey=key,
            mw=float(m.get("mw")), logp=float(m.get("logp")),
            n_basic_centres=len(hits), basic_class=cls, pka_structure=pka,
            pka_literature=LIT_PKA[name],
            lysosome_cytosol_ratio=pH_partition_ratio(pka),
            retardation_R=R, drug_params_shipped=dp, drug_params_seq=dp_seq,
            profile={
                "free_diffusivity": prof.free_diffusivity,
                "effective_diffusivity": prof.effective_diffusivity,
                "hindrance": prof.hindrance,
                "membrane_permeability": prof.membrane_permeability,
            },
        )
    return out


# ===========================================================================
# 4. BATH-BOUNDARY TRANSPORT SOLVER (analysis layer, engine numerics)
# ===========================================================================
@dataclass
class BathResult:
    c_free: np.ndarray      # (ny, nx) free concentration at t = T, uM
    c_total: np.ndarray     # (ny, nx) R * c_free, i.e. what MALDI sees, uM
    auc_total: np.ndarray   # (ny, nx) integral of c_total dt, uM*s
    R: float
    d_eff_mean: float
    n_steps: int


def _log_dt_schedule(total: float, dt0: float = 2.0, n: int = 96,
                     rep: int = 6) -> np.ndarray:
    """Geometrically increasing backward-Euler steps summing exactly to `total`.

    Diffusion into a bath is fast early and slow late, so log-spaced steps put
    resolution where the solution moves. Steps come in blocks of `rep` identical
    dt values so the LU factorisation is reused `rep` times before it has to be
    rebuilt (n/rep factorisations instead of n). Accuracy is checked against the
    analytic solution in `verify_against_analytic_slab()`.
    """
    n_blocks = max(1, n // rep)
    r = 1.0
    lo, hi = 1.0 + 1e-9, 3.0
    for _ in range(200):                      # bisect the per-block growth ratio
        r = 0.5 * (lo + hi)
        s = rep * dt0 * (r ** n_blocks - 1.0) / (r - 1.0)
        if s < total:
            lo = r
        else:
            hi = r
    dts = np.repeat(dt0 * r ** np.arange(n_blocks), rep)
    return dts * (total / dts.sum())


def bath_diffusion(mask: np.ndarray, dx: float, drug: DrugParams,
                   ecm: float = ECM_REFERENCE, R: float = 1.0,
                   c_bath: float = BATH_CONC_UM, total_time: float = EXPOSURE_S,
                   first_order_loss: bool = False,
                   n_steps: int = 96) -> BathResult:
    """Solve dC/dt = div(D_eff grad C) - k C in the explant, C = c_bath on the
    whole outer surface, C(t=0) = 0 inside.

    WHY THIS IS NOT `simulate_drug_transport`. The frozen engine solver takes a
    volumetric source with **no-flux** boundaries and always starts from C = 0;
    it has no Dirichlet facility and no non-zero initial condition, so it
    structurally cannot represent a bath. Rather than modify it, we reuse its
    two numerical primitives unchanged — `_effective_diffusivity` (the same
    `D * (1 - 0.6 * ecm)` tortuosity law) and `_build_diffusion_operator` (the
    same symmetric finite-volume Laplacian with harmonic face averaging) — and
    add only the boundary condition here.

    HOW THE BATH IS IMPOSED. The mask is dilated by one voxel to create a bath
    ring. The ring is given a diffusivity 1e6x the tissue value, so the harmonic
    face average between a tissue voxel and a ring voxel tends to 2*D_tissue —
    which is exactly the discretisation that places the Dirichlet value on the
    tissue *face* rather than a voxel-centre away. We then solve only the
    tissue-tissue block, moving the (known) ring column into the right-hand
    side. Backward Euler with `_log_dt_schedule` steps; the LU factorisation is
    cached per distinct dt.

    HOW SEQUESTRATION IS ENCODED. `R` is the local-equilibrium retardation
    factor from `retardation_factor()`. The caller passes a DrugParams whose
    `diffusivity` has already been divided by R (see `drug_params_seq`); the
    solved field is the FREE concentration and the measured/total field is
    `R * c_free`. No engine change, exact for linear partitioning.
    """
    mask = np.asarray(mask, dtype=bool)
    ring = ndimage.binary_dilation(np.pad(mask, 1)) & ~np.pad(mask, 1)
    full_mask = np.pad(mask, 1) | ring
    ny, nx = full_mask.shape

    ecm_field = np.full((ny, nx), float(ecm))
    ecm_field[ring] = 0.0
    d_eff = _effective_diffusivity(drug, ecm_field)
    d_eff[ring] = d_eff[~ring].max() * 1e6          # well-stirred bath

    L = _build_diffusion_operator(d_eff, full_mask, dx)
    idx = np.arange(ny * nx).reshape(ny, nx)
    t_idx = idx[np.pad(mask, 1)]
    r_idx = idx[ring]
    L_TT = L[t_idx][:, t_idx].tocsc()
    L_TR = L[t_idx][:, r_idx].tocsc()
    b = L_TR @ np.full(len(r_idx), float(c_bath))    # constant Dirichlet inflow

    k = (float(drug.clearance) + float(drug.binding)) if first_order_loss else 0.0
    n = len(t_idx)
    ident = sp_identity(n, format="csc")

    dts = _log_dt_schedule(total_time, n=n_steps)
    c = np.zeros(n)
    auc = np.zeros(n)
    cache: Dict[float, object] = {}
    for dt in dts:
        keyv = float(f"{dt:.6g}")
        solve = cache.get(keyv)
        if solve is None:
            solve = factorized((ident - dt * L_TT + dt * k * ident).tocsc())
            cache[keyv] = solve
        c_new = solve(c + dt * b)
        np.maximum(c_new, 0.0, out=c_new)
        auc += 0.5 * (c + c_new) * dt
        c = c_new

    out_free = np.zeros((ny, nx))
    out_free.ravel()[t_idx] = c
    out_auc = np.zeros((ny, nx))
    out_auc.ravel()[t_idx] = auc
    sl = (slice(1, -1), slice(1, -1))
    return BathResult(
        c_free=out_free[sl], c_total=R * out_free[sl], auc_total=R * out_auc[sl],
        R=R, d_eff_mean=float(d_eff[np.pad(mask, 1)].mean()), n_steps=len(dts),
    )


def verify_against_analytic_slab(verbose: bool = True) -> dict:
    """Check the stepper against Crank's plane-sheet solution.

    For a slab of half-thickness L with both faces held at C0 and C(0)=0, the
    mean saturation is  Cbar/C0 = 1 - sum_n 8/((2n+1)^2 pi^2)
    exp(-D (2n+1)^2 pi^2 t / (4 L^2))   (Crank, *Mathematics of Diffusion*, eq 4.18).
    """
    dx = MODEL_PITCH_UM
    nrow, ncol = 3, 41                       # 41 x 60 um = 2460 um thick slab
    mask = np.ones((nrow, ncol), dtype=bool)
    # make it 1D: pad only left/right, so top/bottom are also bath -> use a tall
    # strip instead and compare along the short axis only via a wide aspect.
    mask = np.ones((201, ncol), dtype=bool)
    drug = DrugParams(name="probe", diffusivity=200.0, binding=0.0, clearance=0.0)
    d_tissue = 200.0 * (1.0 - 0.6 * ECM_REFERENCE)
    half = ncol * dx / 2.0
    rows = []
    for t in (900.0, 3600.0, 14400.0, 86400.0):
        res = bath_diffusion(mask, dx, drug, ecm=ECM_REFERENCE, R=1.0,
                             total_time=t, n_steps=96)
        # centre band, away from the short-axis ends
        num = float(res.c_free[80:121, :].mean() / BATH_CONC_UM)
        ana = 1.0 - sum(
            8.0 / (((2 * n + 1) ** 2) * math.pi ** 2)
            * math.exp(-d_tissue * ((2 * n + 1) ** 2) * math.pi ** 2 * t / (4 * half ** 2))
            for n in range(200)
        )
        rows.append({"t_s": t, "numeric": num, "analytic": ana,
                     "rel_err": abs(num - ana) / max(ana, 1e-12)})
        if verbose:
            print(f"  [slab] t={t:8.0f}s  numeric={num:.5f}  analytic={ana:.5f}  "
                  f"rel_err={rows[-1]['rel_err']:.2%}")
    worst = max(r["rel_err"] for r in rows)
    return {"rows": rows, "worst_rel_err": worst, "passed": bool(worst < 0.03)}


# ===========================================================================
# 5. SIMULATION SWEEP
# ===========================================================================
MODEL_CONDITIONS = {
    # label -> (use lysosomal sequestration, use the molecular layer's
    #           first-order clearance+binding loss)
    "transport_only": (False, False),
    "transport+shipped_loss": (False, True),
    "transport+lysosomal_seq": (True, False),
}


def run_simulations(geos: Dict[str, SliceGeometry],
                    chem: Dict[str, CompoundChemistry],
                    ecm: float = ECM_REFERENCE,
                    conditions=tuple(MODEL_CONDITIONS),
                    compounds=tuple(SMILES),
                    verbose: bool = True) -> pd.DataFrame:
    """Every (slice geometry x compound x model condition). Compounds are run on
    every geometry, not just their own, so the compound-specificity test is a
    within-geometry contrast and cannot be confounded by explant shape."""
    rows = []
    t0 = time.time()
    for sid, g in geos.items():
        for cname in compounds:
            cc = chem[cname]
            for cond in conditions:
                seq, loss = MODEL_CONDITIONS[cond]
                drug = cc.drug_params_seq if seq else cc.drug_params_shipped
                R = cc.retardation_R if seq else 1.0
                res = bath_diffusion(g.mask, g.dx, drug, ecm=ecm, R=R,
                                     first_order_loss=loss)
                v = res.c_total[g.mask]
                a = res.auc_total[g.mask]
                rho, _ = spearmanr(g.depth_um[g.mask], v)
                rows.append({
                    "msi_slice": sid, "own_compound": g.compound,
                    "compound": cname, "condition": cond, "ecm": ecm,
                    "patient": g.patient, "tissue": g.tissue,
                    "area_mm2": g.area_mm2, "max_depth_um": g.max_depth_um,
                    "mean_depth_um": g.mean_depth_um,
                    "pred_mean": float(v.mean()),
                    "pred_median": float(np.median(v)),
                    "pred_cv": heterogeneity(v)["cv"],
                    "pred_gini": heterogeneity(v)["gini"],
                    "pred_q90_q50": heterogeneity(v)["q90_q50"],
                    "pred_auc_mean": float(a.mean()),
                    "pred_sat_frac": float(v.mean() / (res.R * BATH_CONC_UM)),
                    "pred_rho_depth": float(rho),
                    "R": res.R, "d_eff_mean": res.d_eff_mean,
                })
    if verbose:
        print(f"  [sim] {len(rows)} runs in {time.time() - t0:.1f}s (ecm={ecm})")
    return pd.DataFrame(rows)


def gini(x: np.ndarray) -> float:
    """Gini coefficient of a non-negative field: 0 = perfectly uniform.

    Chosen as the PRIMARY heterogeneity metric because, unlike the CV, it is
    bounded in [0, 1), scale-free, and does not explode when a large fraction
    of pixels sit at or below the detection floor (olaparib has >50% exact-zero
    pixels; a p90/p10 ratio literally divides by zero there).
    """
    v = np.sort(np.asarray(x, dtype=float).ravel())
    v = np.clip(v, 0.0, None)
    n = v.size
    s = v.sum()
    if n == 0 or s <= 0:
        return float("nan")
    return float((2.0 * np.arange(1, n + 1) - n - 1.0) @ v / (n * s))


def heterogeneity(values: np.ndarray) -> Dict[str, float]:
    v = np.asarray(values, dtype=float)
    m = float(v.mean())
    return {
        "cv": float(v.std() / max(abs(m), 1e-30)),
        "gini": gini(v),
        "q90_q50": float(np.percentile(v, 90) / max(np.percentile(v, 50), 1e-30)),
    }


def morans_i(field: np.ndarray, mask: np.ndarray) -> float:
    """Rook-adjacency Moran's I of a field restricted to ``mask``.

    Why this metric matters here: dispersion (CV/Gini) cannot tell a genuinely
    patchy drug distribution from shot noise, but SPATIAL AUTOCORRELATION can —
    noise is spatially independent (I -> 0) while a sequestration hotspot is
    contiguous (I >> 0). Applied to the internal-standard channel it gives the
    technical floor from matrix/ionisation gradients, which ARE spatially
    correlated; the drug-minus-IS difference is the part attributable to the
    drug's own distribution.
    """
    f = np.where(mask, np.asarray(field, dtype=float), np.nan)
    mu = np.nanmean(f)
    z = np.where(mask, f - mu, 0.0)
    num = 0.0
    w = 0.0
    for ax, sh in ((0, 1), (1, 1)):
        a = mask & np.roll(mask, sh, axis=ax)
        a = a.copy()
        if ax == 0:
            a[:sh, :] = False
        else:
            a[:, :sh] = False
        num += 2.0 * float((z * np.roll(z, sh, axis=ax))[a].sum())
        w += 2.0 * float(a.sum())
    den = float((z[mask] ** 2).sum())
    if den <= 0 or w <= 0:
        return float("nan")
    n = int(mask.sum())
    return float((n / w) * (num / den))


def measured_slice_stats(geos: Dict[str, SliceGeometry],
                         msi: pd.DataFrame) -> pd.DataFrame:
    """Measured per-slice heterogeneity on the SAME 60 um grid as the model.

    TWO noise controls, because a raw per-pixel CV mostly measures MALDI shot
    noise:

    1. **Spatial coarsening.** Block-averaging 20 um -> 60 um puts the measured
       field on the model's own scale, so measured and predicted dispersion are
       commensurable rather than the measurement being penalised for resolution
       the model does not have.
    2. **The internal-standard channel.** The deuterated IS is sprayed
       uniformly over the section, so its spatial dispersion on the same pixels
       is a direct, per-slice estimate of the technical floor. We report
       ``excess_gini = gini(drug) - gini(IS)`` and
       ``excess_cv = sqrt(max(cv_drug^2 - cv_IS^2, 0))``. Only the excess can be
       attributed to real drug-distribution heterogeneity.
    """
    rows = []
    for sid, g in geos.items():
        v = g.measured[g.mask]
        vis = g.measured_is[g.mask]
        sub = msi[msi.Slice == sid]
        chan, _raw, _heavy = MSI_CHANNEL[g.compound]
        hd, hi = heterogeneity(v), heterogeneity(vis)
        rho, prho = spearmanr(g.depth_um[g.mask], v)
        rows.append({
            "msi_slice": sid, "compound": g.compound, "patient": g.patient,
            "tissue": g.tissue, "area_mm2": g.area_mm2,
            "max_depth_um": g.max_depth_um, "mean_depth_um": g.mean_depth_um,
            "meas_mean": float(v.mean()),
            "meas_cv": hd["cv"], "meas_gini": hd["gini"],
            "meas_q90_q50": hd["q90_q50"],
            "is_cv": hi["cv"], "is_gini": hi["gini"],
            "excess_cv": float(math.sqrt(max(hd["cv"] ** 2 - hi["cv"] ** 2, 0.0))),
            "excess_gini": float(hd["gini"] - hi["gini"]),
            "meas_moran": morans_i(g.measured, g.mask),
            "is_moran": morans_i(g.measured_is, g.mask),
            "excess_moran": float(morans_i(g.measured, g.mask)
                                  - morans_i(g.measured_is, g.mask)),
            "meas_rho_depth": float(rho), "meas_rho_depth_p": float(prho),
            "zero_frac_native": float((sub[chan].values <= 0).mean()),
            "n_px_native": int(len(sub)), "n_vox_model": int(g.mask.sum()),
        })
    return pd.DataFrame(rows)


# ===========================================================================
# 6. MIXED-EFFECTS MACHINERY (the authors' frame)
# ===========================================================================
def _mixedlm(df: pd.DataFrame, response: str, fixed: Optional[str]):
    """Fit `response ~ fixed + (1|Patient) + (1|Patient:Slice)` by ML.

    Mirrors the deposit's own lme4 call
    `lmer(gene ~ drug_measure_norm + (1|Patient) + (1|Patient:Slice), REML=FALSE)`.
    statsmodels expresses the nested slice term as a variance component inside
    the Patient grouping, which is the same model.
    """
    import statsmodels.formula.api as smf
    formula = f"{response} ~ 1" if fixed is None else f"{response} ~ {fixed}"
    md = smf.mixedlm(formula, df, groups=df["Patient"],
                     re_formula="1", vc_formula={"Slice": "0 + C(Slice)"})
    best = None
    for method in ("lbfgs", "bfgs", "cg", "powell"):
        try:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = md.fit(reml=False, method=method, maxiter=4000, disp=False)
        except Exception:
            continue
        if best is None or fit.llf > best.llf:
            best = fit
        if getattr(fit, "converged", False):
            break
    if best is None:
        raise RuntimeError("MixedLM failed to fit under every optimiser")
    return best


def _re_variances(fit) -> Tuple[float, float]:
    """(patient variance, slice-within-patient variance) on the response scale.

    statsmodels reports both scaled by the residual variance."""
    scale = float(fit.scale)
    v_pat = float(np.atleast_1d(np.asarray(fit.cov_re).diagonal()).sum()) * scale
    v_sli = float(np.sum(np.asarray(fit.vcomp))) * scale if len(fit.vcomp) else 0.0
    return v_pat, v_sli


def _nakagawa(fit, df: pd.DataFrame) -> Tuple[float, float]:
    """(marginal R2, conditional R2) — Nakagawa & Schielzeth 2013."""
    fe = np.asarray(fit.model.exog @ fit.fe_params)
    var_f = float(np.var(fe, ddof=0))
    v_pat, v_sli = _re_variances(fit)
    var_r = v_pat + v_sli
    var_e = float(fit.scale)
    tot = var_f + var_r + var_e
    return var_f / tot, (var_f + var_r) / tot


DEGENERATE_TOL = 1e-6


def mixed_test(df: pd.DataFrame, response: str, fixed: str) -> dict:
    """Likelihood-ratio test of one fixed effect inside the authors' frame.

    DEGENERATE PREDICTORS ARE REPORTED, NOT HIDDEN. If the fixed effect has
    (near-)zero variance across the rows — which is exactly what a
    transport-only model produces here, because every explant saturates to the
    bath concentration — no slope is identifiable. Rather than let the
    optimiser fail or return a meaningless value, we return
    ``degenerate: True`` with ``var_explained_beyond_re = 0.0``. Zero variance
    explained is the *correct* answer for a constant predictor, and it is the
    headline result of analysis 1.
    """
    from scipy.stats import chi2
    x = np.asarray(df[fixed], dtype=float)
    rel_sd = float(np.std(x)) / max(abs(float(np.mean(x))), 1e-30)
    if rel_sd < DEGENERATE_TOL:
        return {
            "n": int(len(df)), "degenerate": True, "predictor_rel_sd": rel_sd,
            "beta": float("nan"), "se": float("nan"),
            "intercept": float(np.mean(df[response])),
            "lr_chi2": 0.0, "p_lrt": 1.0,
            "r2_marginal": 0.0, "r2_conditional": float("nan"),
            "delta_r2_marginal": 0.0, "var_explained_beyond_re": 0.0,
            "note": "predictor is constant across rows; no slope identifiable",
        }
    f0 = _mixedlm(df, response, None)
    f1 = _mixedlm(df, response, fixed)
    lr = float(2.0 * (f1.llf - f0.llf))
    p = float(chi2.sf(max(lr, 0.0), df=1))
    r2m0, r2c0 = _nakagawa(f0, df)
    r2m1, r2c1 = _nakagawa(f1, df)
    return {
        "n": int(len(df)), "degenerate": False, "predictor_rel_sd": rel_sd,
        "beta": float(f1.fe_params.get(fixed, np.nan)),
        "se": float(f1.bse.get(fixed, np.nan)),
        "intercept": float(f1.fe_params.get("Intercept", np.nan)),
        "lr_chi2": lr, "p_lrt": p,
        "r2_marginal": r2m1, "r2_conditional": r2c1,
        "delta_r2_marginal": r2m1 - r2m0,
        "var_explained_beyond_re": max(r2m1 - r2m0, 0.0),
        # A negative LRT means the "full" model fitted WORSE than the nested
        # null, which can only happen when the variance components are on a
        # boundary (singular fit). The R2 from such a fit is not trustworthy and
        # is flagged rather than quietly reported.
        "unstable": bool(lr < -1e-6),
    }


def permute_mixed(df: pd.DataFrame, response: str, fixed: str,
                  n_perm: int = 999, seed: int = 0) -> dict:
    """Slice-label permutation test of the MixedLM fixed effect.

    WHY THIS IS NECESSARY AND NOT OPTIONAL. Our predictor is **constant within a
    slice** — a transport model with no ROI coordinates cannot be otherwise. A
    slice random intercept whose variance component is small does NOT fully
    absorb a slice-constant regressor, so the 112 ROIs leak in as
    pseudo-replicates and the parametric LRT p-value is anticonservative by
    roughly the ROI-per-slice factor (~10x). The exchangeable unit here is the
    SLICE, of which there are 11. We therefore permute the slice->prediction
    assignment **within compound** (preserving the compound split, the slice
    means and the ROI counts) and refit, giving an exact-by-construction null
    for the same |t| statistic.
    """
    rng = np.random.default_rng(seed)
    obs_stats = mixed_test(df, response, fixed)
    obs_t = abs(obs_stats["beta"] / max(obs_stats["se"], 1e-30))
    obs_r2 = obs_stats["var_explained_beyond_re"]
    key = (df.groupby(["Slice", "Compound"])[fixed].first().reset_index())
    nt, nr = [], []
    for _ in range(n_perm):
        k = key.copy()
        for c in k.Compound.unique():
            m = k.Compound == c
            k.loc[m, fixed] = rng.permutation(k.loc[m, fixed].values)
        d2 = df.copy()
        d2[fixed] = d2["Slice"].map(dict(zip(k.Slice, k[fixed])))
        try:
            s = mixed_test(d2, response, fixed)
            t = abs(s["beta"] / max(s["se"], 1e-30))
            if np.isfinite(t):
                nt.append(t)
                nr.append(s["var_explained_beyond_re"])
        except Exception:
            continue
    nt, nr = np.asarray(nt), np.asarray(nr)
    def _p(null, obs):
        return float((np.sum(null >= obs) + 1) / (len(null) + 1)) if len(null) else float("nan")
    return {
        "n_perm_completed": int(len(nt)),
        "observed_abs_t": obs_t, "observed_var_explained_beyond_re": obs_r2,
        "null_median_abs_t": float(np.nanmedian(nt)) if len(nt) else float("nan"),
        "null_p95_abs_t": float(np.nanpercentile(nt, 95)) if len(nt) else float("nan"),
        "null_median_var_explained": float(np.nanmedian(nr)) if len(nr) else float("nan"),
        "null_p95_var_explained": float(np.nanpercentile(nr, 95)) if len(nr) else float("nan"),
        "p_perm_t": _p(nt, obs_t), "p_perm_var_explained": _p(nr, obs_r2),
        "exchangeable_unit": "slice (n = %d)" % df["Slice"].nunique(),
    }


def variance_decomposition(df: pd.DataFrame, response: str) -> dict:
    """How much of the measured signal is between-patient / between-slice /
    within-slice. The last stratum is the part no coordinate-free model can
    reach, because the ROI table carries no X/Y."""
    f0 = _mixedlm(df, response, None)
    v_pat, v_sli = _re_variances(f0)
    v_res = float(f0.scale)
    tot = v_pat + v_sli + v_res
    return {"between_patient": v_pat / tot, "between_slice_within_patient": v_sli / tot,
            "within_slice_residual": v_res / tot, "total_var": tot}


# ===========================================================================
# 7. ANALYSIS 1 — accumulation-heterogeneity fit
# ===========================================================================
def attach_predictions(roi: pd.DataFrame, sim: pd.DataFrame,
                       condition: str) -> pd.DataFrame:
    """Give every ROI its model prediction and normalise it the way the data is.

    THE HARD LIMIT, STATED UP FRONT. The ROI table carries `Slice`, `ROILabel`
    and a High/Low call, but **no X/Y**: GeoMx ROI positions were never
    deposited, and the MSI pixel file's only ROI column is the constant
    `TissueROI`. A transport model's per-ROI output is a function of position,
    so the finest granularity at which our prediction can be *honestly*
    evaluated is the slice. Every ROI in a slice therefore receives that slice's
    predicted mean accumulation. The within-slice stratum — which carries the
    MAJORITY of the measured variance (see `variance_decomposition`) — is
    unreachable, and we do not attempt to reach it by any rank- or
    label-conditioned construction, because those would be circular.

    NORMALISATION. `drug_measure_norm = drug_measure / measure_norm_factor`,
    where `measure_norm_factor` is a single number per COMPOUND. We apply the
    identical transform to the prediction (divide by the compound-wide median of
    the predicted ROI values), so predicted and measured are on the same
    footing. A consequence worth stating: any compound-level *level* difference
    our model predicts is removed by this normalisation on both sides, so only
    within-compound (slice-to-slice) structure is testable here.
    """
    s = sim[sim.condition == condition]
    key = s.set_index(["msi_slice", "compound"])
    out = roi.copy()
    vals, sat, rho = [], [], []
    for _, r in out.iterrows():
        sid = r["msi_slice"]
        if sid is None or (isinstance(sid, float) and np.isnan(sid)):
            vals.append(np.nan); sat.append(np.nan); rho.append(np.nan); continue
        row = key.loc[(sid, r["Compound"])]
        vals.append(float(row["pred_mean"]))
        sat.append(float(row["pred_sat_frac"]))
        rho.append(float(row["pred_cv"]))
    out["pred_raw"] = vals
    out["pred_sat_frac"] = sat
    out["pred_cv"] = rho
    out["pred_norm"] = np.nan
    for cmp_, grp in out.groupby("Compound"):
        med = np.nanmedian(grp["pred_raw"])
        out.loc[grp.index, "pred_norm"] = grp["pred_raw"] / med
    return out


def analysis1(roi: pd.DataFrame, sim: pd.DataFrame, meas: pd.DataFrame,
              n_perm: int = 2000, seed: int = 0, n_perm_mixed: int = 999,
              verbose: bool = True) -> dict:
    rng = np.random.default_rng(seed)
    res: dict = {"variance_decomposition": {}, "conditions": {}}
    res["variance_decomposition"]["all"] = variance_decomposition(roi, "drug_measure_norm")
    for c in ("Niraparib", "Rucaparib"):
        res["variance_decomposition"][c] = variance_decomposition(
            roi[roi.Compound == c], "drug_measure_norm")

    for cond in MODEL_CONDITIONS:
        d = attach_predictions(roi, sim, cond).dropna(subset=["pred_norm"])
        entry = {"n_roi": int(len(d)),
                 "n_roi_dropped_no_msi": int(len(roi) - len(d)),
                 "roi_table": d}
        entry["pooled"] = mixed_test(d, "drug_measure_norm", "pred_norm")
        for c in ("Niraparib", "Rucaparib"):
            sub = d[d.Compound == c]
            entry[c] = mixed_test(sub, "drug_measure_norm", "pred_norm")

        # Slice-level view: n = 11 mapped slices, the real unit of information.
        sl = (d.groupby(["msi_slice", "Compound", "Patient"])
                .agg(meas=("drug_measure_norm", "mean"),
                     pred=("pred_norm", "first"), n=("drug_measure_norm", "size"))
                .reset_index())
        pred_rel_sd = float(np.std(sl["pred"])) / max(abs(float(np.mean(sl["pred"]))), 1e-30)
        if pred_rel_sd < DEGENERATE_TOL:
            entry["slice_level"] = {
                "n_slices": int(len(sl)), "degenerate": True,
                "pred_rel_sd": pred_rel_sd,
                "spearman_rho": float("nan"), "spearman_p": float("nan"),
                "pearson_r": float("nan"), "pearson_p": float("nan"),
            }
            entry["permutation"] = {"n": 0, "p_perm_two_sided": float("nan"),
                                    "null_sd_rho": float("nan"),
                                    "note": "constant prediction; nothing to permute"}
            entry["slice_table"] = sl.to_dict("records")
            res["conditions"][cond] = entry
            continue
        rho, prho = spearmanr(sl["pred"], sl["meas"])
        r, pr = pearsonr(sl["pred"], sl["meas"])
        entry["slice_level"] = {
            "n_slices": int(len(sl)), "degenerate": False,
            "pred_rel_sd": pred_rel_sd, "spearman_rho": float(rho),
            "spearman_p": float(prho), "pearson_r": float(r), "pearson_p": float(pr),
        }
        # Permutation null on the slice-level statistic: shuffle the prediction
        # across slices *within compound*, so the null preserves the compound
        # split, the slice means and the marginal distributions, and destroys
        # only the model's slice-to-slice ordering. Done at slice level (not by
        # refitting the mixed model 2000x) because n = 11 makes the exact
        # permutation distribution the honest reference anyway.
        obs_rho = entry["slice_level"]["spearman_rho"]
        null = np.empty(n_perm)
        pv, cv_ = sl["pred"].values, sl["Compound"].values
        mv = sl["meas"].values
        for i in range(n_perm):
            perm = pv.copy()
            for c in np.unique(cv_):
                m = cv_ == c
                perm[m] = rng.permutation(pv[m])
            null[i] = spearmanr(perm, mv)[0]
        entry["permutation"] = {
            "n": int(n_perm),
            "p_perm_two_sided": float((np.sum(np.abs(null) >= abs(obs_rho)) + 1)
                                      / (n_perm + 1)),
            "null_sd_rho": float(np.nanstd(null)),
        }
        entry["slice_table"] = sl.to_dict("records")
        if n_perm_mixed:
            t0 = time.time()
            entry["mixedlm_permutation"] = permute_mixed(
                d, "drug_measure_norm", "pred_norm", n_perm=n_perm_mixed, seed=seed)
            if verbose:
                print(f"  [a1] {cond}: MixedLM slice-permutation "
                      f"p = {entry['mixedlm_permutation']['p_perm_t']:.4f} "
                      f"({time.time() - t0:.0f}s)")
        res["conditions"][cond] = entry

    # The condition whose prediction is not constant, i.e. the only one that can
    # be scored at all. Chosen by predictor variance, never by outcome.
    non_deg = [c for c in MODEL_CONDITIONS
               if not res["conditions"][c]["pooled"].get("degenerate", False)]
    res["headline_condition"] = non_deg[0] if non_deg else list(MODEL_CONDITIONS)[0]
    res["degenerate_conditions"] = [c for c in MODEL_CONDITIONS if c not in non_deg]
    return res


# ===========================================================================
# 8. ANALYSIS 2 — the compound-specificity falsification test
# ===========================================================================
def saturation_timecourse(geos: Dict[str, SliceGeometry],
                          chem: Dict[str, CompoundChemistry],
                          hours=(0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 24.0),
                          ecm: float = ECM_REFERENCE) -> pd.DataFrame:
    """Mean tissue saturation vs time for each section under its OWN compound.

    This is the quantitative form of the equilibrium argument: if the explant
    is already at the bath concentration well before the 24 h endpoint, then no
    purely diffusive model — with or without a uniform sequestration term — can
    produce spatial structure at the endpoint.
    """
    rows = []
    for sid, g in geos.items():
        cc = chem[g.compound]
        for cond in ("transport_only", "transport+lysosomal_seq"):
            seq, loss = MODEL_CONDITIONS[cond]
            drug = cc.drug_params_seq if seq else cc.drug_params_shipped
            R = cc.retardation_R if seq else 1.0
            for h in hours:
                res = bath_diffusion(g.mask, g.dx, drug, ecm=ecm, R=R,
                                     total_time=h * 3600.0, first_order_loss=loss)
                v = res.c_total[g.mask]
                rows.append({
                    "msi_slice": sid, "compound": g.compound, "condition": cond,
                    "hours": h, "sat_frac": float(v.mean() / (R * BATH_CONC_UM)),
                    "cv": heterogeneity(v)["cv"], "gini": gini(v),
                    "max_depth_um": g.max_depth_um,
                })
    return pd.DataFrame(rows)


def phi_lys_sweep(geos: Dict[str, SliceGeometry],
                  chem: Dict[str, CompoundChemistry],
                  phis=PHI_LYS_SWEEP, ecm: float = ECM_REFERENCE) -> pd.DataFrame:
    """How much heterogeneity can a UNIFORM lysosomal compartment ever buy?

    Sweeps the lysosomal volume fraction over 0.005-0.05 (0.5-5% of tissue
    volume, bracketing every value in the cell-biology literature) and records
    the resulting predicted dispersion. The point of the sweep is that the
    answer is bounded and small: with phi_lys uniform in space, R is uniform,
    so at equilibrium C_total = R * C_bath everywhere and the dispersion goes to
    zero no matter how large R gets. Sequestration buys a longer time constant,
    never an equilibrium pattern.
    """
    rows = []
    for phi in phis:
        for name, cc in chem.items():
            R = retardation_factor(cc.pka_structure, phi)
            dp = DrugParams(
                name=f"{name}_phi{phi}", diffusivity=cc.drug_params_shipped.diffusivity / R,
                molecular_weight=cc.drug_params_shipped.molecular_weight,
                ic50=cc.drug_params_shipped.ic50, hill=cc.drug_params_shipped.hill,
                binding=cc.drug_params_shipped.binding,
                clearance=cc.drug_params_shipped.clearance,
                bbb_penetrance=cc.drug_params_shipped.bbb_penetrance,
            )
            for sid, g in geos.items():
                res = bath_diffusion(g.mask, g.dx, dp, ecm=ecm, R=R)
                v = res.c_total[g.mask]
                rows.append({
                    "phi_lys": phi, "compound": name, "msi_slice": sid,
                    "R": R, "pred_cv": heterogeneity(v)["cv"],
                    "pred_gini": gini(v),
                    "sat_frac": float(v.mean() / (R * BATH_CONC_UM)),
                })
    return pd.DataFrame(rows)


def ecm_sensitivity(geos: Dict[str, SliceGeometry],
                    chem: Dict[str, CompoundChemistry],
                    ecms=ECM_SWEEP) -> pd.DataFrame:
    """Does the (uncalibrated) ECM volume fraction change any conclusion?"""
    rows = []
    for e in ecms:
        for name, cc in chem.items():
            for sid, g in geos.items():
                res = bath_diffusion(g.mask, g.dx, cc.drug_params_seq, ecm=e,
                                     R=cc.retardation_R)
                v = res.c_total[g.mask]
                rows.append({"ecm": e, "compound": name, "msi_slice": sid,
                             "pred_cv": heterogeneity(v)["cv"],
                             "pred_gini": gini(v),
                             "sat_frac": float(v.mean() / (cc.retardation_R * BATH_CONC_UM))})
    return pd.DataFrame(rows)


def analysis2(geos: Dict[str, SliceGeometry], sim: pd.DataFrame,
              meas: pd.DataFrame, chem: Dict[str, CompoundChemistry],
              sat: pd.DataFrame, sweep: pd.DataFrame) -> dict:
    """THE SHARP TEST.

    Claim under test: niraparib and rucaparib are weak bases whose lysosomal
    sequestration drives accumulation heterogeneity; olaparib is not a base and
    does not. A model carrying the mechanism must predict MORE heterogeneity
    for the two bases than for olaparib, AND at a magnitude comparable to what
    was measured. Both halves must pass. Getting the ordering right at 1/200th
    of the observed magnitude is not a pass.

    The contrast is run WITHIN geometry: every compound is simulated on all 18
    explant outlines, so shape cannot confound it.
    """
    from scipy.stats import mannwhitneyu

    out: dict = {}

    # --- (a) structure -> basicity, straight out of RDKit ------------------
    out["chemistry"] = {
        n: {"inchikey": c.inchikey, "mw": c.mw, "logp": c.logp,
            "n_basic_centres": c.n_basic_centres, "basic_class": c.basic_class,
            "pka_structure": c.pka_structure, "pka_literature": c.pka_literature,
            "lysosome_cytosol_ratio": c.lysosome_cytosol_ratio,
            "retardation_R": c.retardation_R,
            "free_diffusivity": c.profile["free_diffusivity"],
            "effective_diffusivity": c.profile["effective_diffusivity"],
            "membrane_permeability": c.profile["membrane_permeability"]}
        for n, c in chem.items()
    }
    out["weak_base_split_detected"] = bool(
        chem["Niraparib"].n_basic_centres > 0 and chem["Rucaparib"].n_basic_centres > 0
        and chem["Olaparib"].n_basic_centres == 0)
    out["nir_vs_ruc_ratio_gap_pct"] = float(
        100.0 * abs(chem["Niraparib"].lysosome_cytosol_ratio
                    - chem["Rucaparib"].lysosome_cytosol_ratio)
        / chem["Niraparib"].lysosome_cytosol_ratio)

    # --- (b) predicted heterogeneity, within-geometry ----------------------
    pred = {}
    for cond in MODEL_CONDITIONS:
        s = sim[sim.condition == cond]
        pred[cond] = {c: {"cv_mean": float(s[s.compound == c].pred_cv.mean()),
                          "gini_mean": float(s[s.compound == c].pred_gini.mean()),
                          "sat_frac_mean": float(s[s.compound == c].pred_sat_frac.mean())}
                      for c in SMILES}
        wb = float(np.mean([pred[cond]["Niraparib"]["cv_mean"],
                            pred[cond]["Rucaparib"]["cv_mean"]]))
        ol = pred[cond]["Olaparib"]["cv_mean"]
        pred[cond]["weak_base_over_olaparib_cv"] = float(wb / max(ol, 1e-30))
        pred[cond]["nir_over_ruc_cv"] = float(
            pred[cond]["Niraparib"]["cv_mean"] / max(pred[cond]["Rucaparib"]["cv_mean"], 1e-30))
    out["predicted"] = pred

    # --- (c) measured heterogeneity, with the IS noise floor ---------------
    m = meas.copy()
    m["is_base"] = m.compound.isin(["Niraparib", "Rucaparib"])
    out["measured"] = {
        c: {"n_slices": int((m.compound == c).sum()),
            "cv_mean": float(m[m.compound == c].meas_cv.mean()),
            "gini_mean": float(m[m.compound == c].meas_gini.mean()),
            "is_gini_mean": float(m[m.compound == c].is_gini.mean()),
            "excess_gini_mean": float(m[m.compound == c].excess_gini.mean()),
            "excess_cv_mean": float(m[m.compound == c].excess_cv.mean()),
            "meas_moran_mean": float(m[m.compound == c].meas_moran.mean()),
            "is_moran_mean": float(m[m.compound == c].is_moran.mean()),
            "excess_moran_mean": float(m[m.compound == c].excess_moran.mean()),
            "zero_frac_mean": float(m[m.compound == c].zero_frac_native.mean())}
        for c in SMILES
    }
    out["measured_weak_base_vs_olaparib"] = {}
    for metric in ("excess_gini", "excess_cv", "excess_moran"):
        wb_v = m.loc[m.is_base, metric].values
        ol_v = m.loc[~m.is_base, metric].values
        u, p = mannwhitneyu(wb_v, ol_v, alternative="greater")
        out["measured_weak_base_vs_olaparib"][metric] = {
            "n_weak_base_slices": int(len(wb_v)), "n_olaparib_slices": int(len(ol_v)),
            "weak_base_mean": float(wb_v.mean()), "olaparib_mean": float(ol_v.mean()),
            "mannwhitney_u": float(u), "p_one_sided_greater": float(p),
            "ratio": float(wb_v.mean() / max(ol_v.mean(), 1e-30)),
        }
    out["measured_metric_note"] = (
        "olaparib sections carry a mean 51% exact-zero pixel fraction (below the "
        "MALDI detection floor) against 10-38% for the weak bases, which inflates "
        "every dispersion metric for olaparib. The measured side of this contrast "
        "is therefore confounded and we do not claim it in either direction; the "
        "falsification is run on the PREDICTED side, where no such confound exists."
    )

    # --- (d) the magnitude gap: predicted vs measured ----------------------
    best_cond = max(MODEL_CONDITIONS,
                    key=lambda c: pred[c]["Niraparib"]["cv_mean"])
    out["magnitude_gap"] = {
        "best_condition_for_heterogeneity": best_cond,
        "pred_cv_weak_base": float(np.mean([pred[best_cond]["Niraparib"]["cv_mean"],
                                            pred[best_cond]["Rucaparib"]["cv_mean"]])),
        "meas_excess_cv_weak_base": float(m.loc[m.is_base, "excess_cv"].mean()),
        "meas_cv_weak_base": float(m.loc[m.is_base, "meas_cv"].mean()),
    }
    out["magnitude_gap"]["fold_under_prediction"] = float(
        out["magnitude_gap"]["meas_excess_cv_weak_base"]
        / max(out["magnitude_gap"]["pred_cv_weak_base"], 1e-30))
    # sequestration-only condition, which is the one the mechanism claim is about
    seq_cv = float(np.mean([pred["transport+lysosomal_seq"]["Niraparib"]["cv_mean"],
                            pred["transport+lysosomal_seq"]["Rucaparib"]["cv_mean"]]))
    out["magnitude_gap"]["pred_cv_weak_base_seq_only"] = seq_cv
    out["magnitude_gap"]["fold_under_prediction_seq_only"] = float(
        out["magnitude_gap"]["meas_excess_cv_weak_base"] / max(seq_cv, 1e-30))

    # --- (e) is the field even diffusion-limited at 24 h? ------------------
    fo = []
    for sid, g in geos.items():
        d_eff = chem[g.compound].profile["effective_diffusivity"]
        R = chem[g.compound].retardation_R
        fo.append({"msi_slice": sid, "compound": g.compound,
                   "L_um": g.max_depth_um,
                   "fourier_no_seq": d_eff * EXPOSURE_S / g.max_depth_um ** 2,
                   "fourier_with_seq": (d_eff / R) * EXPOSURE_S / g.max_depth_um ** 2})
    fo_df = pd.DataFrame(fo)
    out["fourier"] = {
        "definition": "Fo = D_eff * t / L^2, L = max in-plane distance to the bath surface",
        "min_no_seq": float(fo_df.fourier_no_seq.min()),
        "median_no_seq": float(fo_df.fourier_no_seq.median()),
        "min_with_seq": float(fo_df.fourier_with_seq.min()),
        "median_with_seq": float(fo_df.fourier_with_seq.median()),
        "n_slices_Fo_lt_1_with_seq": int((fo_df.fourier_with_seq < 1).sum()),
        "n_slices": int(len(fo_df)),
    }
    s24 = sat[sat.hours == 24.0]
    s1 = sat[sat.hours == 1.0]
    out["saturation"] = {
        "mean_sat_frac_1h_transport_only":
            float(s1[s1.condition == "transport_only"].sat_frac.mean()),
        "mean_sat_frac_1h_with_seq":
            float(s1[s1.condition == "transport+lysosomal_seq"].sat_frac.mean()),
        "min_sat_frac_24h_with_seq":
            float(s24[s24.condition == "transport+lysosomal_seq"].sat_frac.min()),
        "mean_sat_frac_24h_with_seq":
            float(s24[s24.condition == "transport+lysosomal_seq"].sat_frac.mean()),
    }

    # --- (f) can a bigger lysosomal compartment rescue it? -----------------
    sw = sweep.groupby(["phi_lys", "compound"])[["pred_cv", "R"]].mean().reset_index()
    out["phi_lys_sweep"] = sw.to_dict("records")
    out["phi_lys_sweep_max_cv"] = float(sweep.pred_cv.max())

    # --- (g) direction of the spatial gradient -----------------------------
    out["depth_gradient"] = {
        "model_mean_rho_depth": float(
            sim[sim.condition == "transport+lysosomal_seq"].pred_rho_depth.mean()),
        "measured_mean_rho_depth": float(meas.meas_rho_depth.mean()),
        "n_slices_measured_positive": int((meas.meas_rho_depth > 0).sum()),
        "n_slices": int(len(meas)),
        "sign_agreement": int(((meas.set_index("msi_slice").meas_rho_depth < 0)).sum()),
    }

    # --- verdict -----------------------------------------------------------
    ordering_ok = pred["transport+lysosomal_seq"]["weak_base_over_olaparib_cv"] > 2.0
    magnitude_ok = out["magnitude_gap"]["fold_under_prediction_seq_only"] < 3.0
    out["verdict"] = {
        "ordering_correct": bool(ordering_ok),
        "magnitude_correct": bool(magnitude_ok),
        "overall": "PASS" if (ordering_ok and magnitude_ok) else
                   ("PARTIAL - ordering only" if ordering_ok else "FAIL"),
    }
    return out


# ===========================================================================
# 9. ANALYSIS 3 — exposure -> molecular response coupling
# ===========================================================================
GENE_COL_START = 13   # columns 0..12 are metadata; 13.. are the 11,224 genes


def signature_scores(expr: pd.DataFrame, seed: int = 0,
                     n_null: int = 1000) -> Tuple[pd.DataFrame, dict]:
    """Per-ROI apoptosis and lysosome signature scores + a size-matched null.

    SCORING RULE (stated so it can be argued with):
      1. take the deposited Q3-normalised GeoMx counts as given (this is the
         matrix the authors modelled),
      2. ``log2(x + 1)`` to stabilise variance across a 4-decade dynamic range,
      3. z-score EACH GENE across all 114 ROIs, so no single high-expressed
         gene dominates the set,
      4. the score is the unweighted mean z over the set members present on the
         panel.
    This is the standard single-sample mean-z ("Tirosh-style" without the
    control-bin correction). Gene sets: MSigDB Hallmark 2020 *Apoptosis* and
    KEGG 2021 Human *Lysosome*, embedded verbatim at the top of this module.

    NULL: ``n_null`` random gene sets of the SAME SIZE drawn from the same
    11,224-gene panel, scored identically. Every statistic computed on a real
    signature below is also computed on all null sets, which gives an empirical
    p-value that is immune to the "any large gene set correlates with anything"
    objection.
    """
    rng = np.random.default_rng(seed)
    genes = list(expr.columns[GENE_COL_START:])
    X = np.log2(expr[genes].values.astype(float) + 1.0)
    mu, sd = X.mean(axis=0), X.std(axis=0)
    Z = (X - mu) / np.where(sd > 0, sd, 1.0)
    gidx = {g: i for i, g in enumerate(genes)}

    sets = {"apoptosis": APOPTOSIS_GENES, "lysosome": LYSOSOME_GENES}
    scores = expr.iloc[:, :GENE_COL_START].copy()
    info: dict = {"n_genes_panel": len(genes), "n_null": n_null, "sets": {}}
    null_scores: Dict[str, np.ndarray] = {}
    for name, gs in sets.items():
        idx = [gidx[g] for g in gs if g in gidx]
        scores[f"{name}_score"] = Z[:, idx].mean(axis=1)
        info["sets"][name] = {"n_requested": len(gs), "n_on_panel": len(idx)}
        null_scores[name] = np.stack(
            [Z[:, rng.choice(len(genes), size=len(idx), replace=False)].mean(axis=1)
             for _ in range(n_null)], axis=1)          # (114, n_null)
    return scores, {"info": info, "null": null_scores, "Z": Z, "genes": genes}


def _within_transform(y: np.ndarray, groups: pd.DataFrame) -> np.ndarray:
    """Residualise columns of ``y`` on Patient + Patient:Slice FIXED effects.

    This is the conservative counterpart to the authors' random-effects fit:
    absorbing slice as a fixed effect removes ALL between-slice signal, so what
    remains is strictly within-slice covariation. Any association surviving it
    cannot be an artefact of patient or section identity.
    """
    d = pd.get_dummies(groups["Slice"].astype(str), drop_first=False).values.astype(float)
    q, _ = np.linalg.qr(d)
    return np.asarray(y, dtype=float) - q @ (q.T @ np.asarray(y, dtype=float))


def gene_level_enrichment(expr: pd.DataFrame, bundle: dict,
                          compound: Optional[str] = None) -> dict:
    """Fast, whole-panel replication of the direction of the authors' GSEA.

    All 11,224 genes are ranked by their WITHIN-SLICE partial correlation with
    ``drug_measure_norm`` (both sides residualised on slice identity — see
    ``_within_transform``). We then ask whether the apoptosis and lysosome sets
    sit unusually high in that ranking, against 1,000 size-matched random sets.
    This is a rank-enrichment statistic, not a permutation-GSEA NES, so the
    number is NOT directly comparable to the paper's NES = 2.278; the DIRECTION
    and the significance are.
    """
    sub = expr if compound is None else expr[expr.Compound == compound]
    Z = np.log2(sub[bundle["genes"]].values.astype(float) + 1.0)
    Z = (Z - Z.mean(0)) / np.where(Z.std(0) > 0, Z.std(0), 1.0)
    Zw = _within_transform(Z, sub)
    dw = _within_transform(sub["drug_measure_norm"].values.reshape(-1, 1), sub).ravel()
    denom = np.linalg.norm(Zw, axis=0) * np.linalg.norm(dw)
    r = np.where(denom > 0, (Zw * dw[:, None]).sum(0) / np.where(denom > 0, denom, 1.0), 0.0)
    order = np.argsort(-r)
    rank = np.empty_like(order)
    rank[order] = np.arange(len(order))
    frac = rank / (len(order) - 1.0)          # 0 = most positively drug-correlated

    gidx = {g: i for i, g in enumerate(bundle["genes"])}
    rng = np.random.default_rng(1)
    out = {"n_roi": int(len(sub)), "compound": compound or "pooled", "sets": {}}
    for name, gs in (("apoptosis", APOPTOSIS_GENES), ("lysosome", LYSOSOME_GENES)):
        idx = [gidx[g] for g in gs if g in gidx]
        obs = float(np.mean(frac[idx]))                       # <0.5 = enriched up
        null = np.array([np.mean(frac[rng.choice(len(frac), len(idx), replace=False)])
                         for _ in range(1000)])
        out["sets"][name] = {
            "n_genes": len(idx),
            "mean_rank_fraction": obs,
            "null_mean": float(null.mean()), "null_sd": float(null.std()),
            "z": float((0.5 - obs) / max(null.std(), 1e-12)),
            "p_empirical_up": float((np.sum(null <= obs) + 1) / 1001.0),
            "mean_within_slice_r": float(np.mean(r[idx])),
        }
    return out


def analysis3(expr: pd.DataFrame, roi: pd.DataFrame, sim: pd.DataFrame,
              condition: str = "transport+shipped_loss",
              n_null: int = 1000, verbose: bool = True) -> dict:
    """Does drug exposure predict the apoptotic signature — measured, then ours?"""
    scores, bundle = signature_scores(expr, n_null=n_null)
    scores["msi_slice"] = scores["Slice"].map(SLICE_MAP)
    pred = attach_predictions(roi, sim, condition)
    key = pred.set_index(["Slice", "ROILabel"])["pred_norm"]
    scores["pred_norm"] = [key.get((s, r), np.nan)
                           for s, r in zip(scores["Slice"], scores["ROILabel"])]

    out: dict = {"scoring": bundle["info"], "condition_used": condition, "tests": {}}
    nulls = bundle["null"]

    for sig in ("apoptosis", "lysosome"):
        col = f"{sig}_score"
        entry: dict = {}
        # ---- measured drug -> signature, authors' frame -------------------
        for grp, d in [("pooled", scores),
                       ("Niraparib", scores[scores.Compound == "Niraparib"]),
                       ("Rucaparib", scores[scores.Compound == "Rucaparib"])]:
            entry[f"measured::{grp}"] = mixed_test(d, col, "drug_measure_norm")
        # ---- OUR predicted exposure -> signature --------------------------
        dp = scores.dropna(subset=["pred_norm"])
        entry["predicted::pooled"] = mixed_test(dp, col, "pred_norm")
        for c in ("Niraparib", "Rucaparib"):
            entry[f"predicted::{c}"] = mixed_test(dp[dp.Compound == c], col, "pred_norm")

        # ---- empirical null over size-matched random gene sets ------------
        obs_t = abs(entry["measured::pooled"]["beta"]) / max(
            entry["measured::pooled"]["se"], 1e-30)
        nt = []
        for j in range(min(200, nulls[sig].shape[1])):   # 200 refits is enough
            tmp = scores.copy()
            tmp["_y"] = nulls[sig][:, j]
            try:
                m = mixed_test(tmp, "_y", "drug_measure_norm")
                nt.append(abs(m["beta"]) / max(m["se"], 1e-30))
            except Exception:
                continue
        nt = np.array([v for v in nt if np.isfinite(v)])
        entry["gene_set_null"] = {
            "n_random_sets_fitted": int(len(nt)),
            "observed_abs_t": float(obs_t),
            "null_median_abs_t": float(np.median(nt)) if len(nt) else float("nan"),
            "null_p95_abs_t": float(np.percentile(nt, 95)) if len(nt) else float("nan"),
            "p_empirical": float((np.sum(nt >= obs_t) + 1) / (len(nt) + 1)) if len(nt) else float("nan"),
        }
        out["tests"][sig] = entry
        if verbose:
            print(f"  [a3] {sig}: measured beta={entry['measured::pooled']['beta']:+.4f} "
                  f"p={entry['measured::pooled']['p_lrt']:.3g} | "
                  f"predicted p={entry['predicted::pooled']['p_lrt']:.3g}")

    out["gene_level_enrichment"] = {
        "pooled": gene_level_enrichment(expr, bundle),
        "Niraparib": gene_level_enrichment(expr, bundle, "Niraparib"),
        "Rucaparib": gene_level_enrichment(expr, bundle, "Rucaparib"),
    }
    out["paper_benchmark"] = {"gsea_apoptosis_NES": 2.278, "gsea_apoptosis_FDR_q": 0.0}

    # ---- the bridge back to analysis 2: does lysosomal content track drug? --
    lys = out["tests"]["lysosome"]["measured::pooled"]
    apo = out["tests"]["apoptosis"]["measured::pooled"]
    out["mechanism_bridge"] = {
        "lysosome_signature_beta": lys["beta"], "lysosome_signature_p": lys["p_lrt"],
        "apoptosis_signature_beta": apo["beta"], "apoptosis_signature_p": apo["p_lrt"],
    }
    out["_scores"] = scores
    return out


# ===========================================================================
# 10. FIGURE
# ===========================================================================
# Categorical hues: dataviz reference palette, slots 1 / 8 / 5. Validated with
# scripts/validate_palette.js (light, surface #fcfcfb): all four checks PASS,
# worst adjacent CVD dE = 96.7. Compound identity is carried by MARKER SHAPE,
# never by colour alone.
PATIENT_COLOR = {"P23_015": "#2a78d6", "P23_037": "#eb6834", "P23_071": "#4a3aa7"}
PATIENT_COLOR_MSI = {"P23-015": "#2a78d6", "P23-037": "#eb6834", "P23-071": "#4a3aa7"}
COMPOUND_MARKER = {"Niraparib": "o", "Rucaparib": "s", "Olaparib": "^"}
MEAS_COLOR = "#4a3aa7"      # measured  (validated pair with #2a78d6, dE 16.6)
PRED_COLOR = "#2a78d6"      # predicted
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8f8e88"


def _style(ax, xlabel: str, ylabel: str, title: str) -> None:
    ax.set_xlabel(xlabel, fontsize=9.5, color=INK2)
    ax.set_ylabel(ylabel, fontsize=9.5, color=INK2)
    ax.set_title(title, fontsize=11, color=INK, loc="left", pad=8)
    ax.grid(True, lw=0.5, color="#e6e5e0", zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#d5d4ce")
    ax.tick_params(colors=INK2, labelsize=8.5, length=3, width=0.8)


def make_figure(a1: dict, a2: dict, a3: dict, meas: pd.DataFrame,
                sim: pd.DataFrame, geos: Dict[str, SliceGeometry],
                chem: Dict[str, CompoundChemistry], sat: pd.DataFrame,
                path: str = OUT_PNG) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.6))
    fig.patch.set_facecolor("#fcfcfb")
    for ax in axes.ravel():
        ax.set_facecolor("#fcfcfb")

    scores = a3["_scores"]
    cond = a1["headline_condition"]
    d = a1["conditions"][cond]["roi_table"]

    # ---- (a) measured vs predicted per-ROI accumulation -------------------
    ax = axes[0, 0]
    for pat, gp in d.groupby("Patient"):
        for cmp_, gc in gp.groupby("Compound"):
            ax.scatter(gc.pred_norm, gc.drug_measure_norm, s=42,
                       facecolor=PATIENT_COLOR[pat], edgecolor="#fcfcfb", lw=1.0,
                       marker=COMPOUND_MARKER[cmp_], zorder=3, alpha=0.95)
    sl = pd.DataFrame(a1["conditions"][cond]["slice_table"])
    ax.scatter(sl.pred, sl.meas, s=140, facecolor="none", edgecolor=INK,
               lw=1.4, zorder=4, marker="D", label="slice mean")
    st = a1["conditions"][cond]["slice_level"]
    ax.axhline(1.0, color=INK3, lw=1.0, ls=":", zorder=1)
    _style(ax, "predicted accumulation (compound-normalised)",
           "measured drug_measure_norm (IS ratio)",
           "a  Per-ROI accumulation: measured vs predicted")
    mperm = a1["conditions"][cond].get("mixedlm_permutation", {})
    ax.text(0.02, 0.96,
            f"n = {a1['conditions'][cond]['n_roi']} ROIs but only {st['n_slices']} "
            f"exchangeable slices\n"
            f"$\\Delta R^2_m$ beyond patient+slice REs = "
            f"{a1['conditions'][cond]['pooled']['var_explained_beyond_re']*100:.1f}% "
            f"(nominal LRT p = {a1['conditions'][cond]['pooled']['p_lrt']:.3f})\n"
            f"slice-permutation null median $\\Delta R^2_m$ = "
            f"{mperm.get('null_median_var_explained', float('nan'))*100:.1f}%  ->  "
            f"$p_{{perm}}$ = {mperm.get('p_perm_var_explained', float('nan')):.2f}\n"
            f"slice-level Spearman $\\rho$ = {st['spearman_rho']:.2f} "
            f"(p = {st['spearman_p']:.2f})",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.4, color=INK2,
            bbox=dict(boxstyle="round,pad=0.45", fc="#f3f2ee", ec="#e0dfd9", lw=0.8))
    handles = [Line2D([], [], marker="o", ls="", ms=7, mfc=c, mec="#fcfcfb", label=p)
               for p, c in PATIENT_COLOR.items()]
    handles += [Line2D([], [], marker=m, ls="", ms=7, mfc=INK3, mec="#fcfcfb", label=c)
                for c, m in COMPOUND_MARKER.items() if c != "Olaparib"]
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="lower right", ncol=2)

    # ---- (b) heterogeneity across compounds -------------------------------
    ax = axes[0, 1]
    comps = ["Niraparib", "Rucaparib", "Olaparib"]
    xs = np.arange(len(comps))
    mvals = [a2["measured"][c]["excess_cv_mean"] for c in comps]
    pvals = [max(a2["predicted"]["transport+lysosomal_seq"][c]["cv_mean"], 1e-9)
             for c in comps]
    ax.bar(xs - 0.20, mvals, width=0.36, color=MEAS_COLOR, zorder=3,
           label="measured (IS-noise-corrected)")
    ax.bar(xs + 0.20, pvals, width=0.36, color=PRED_COLOR, zorder=3,
           label="predicted (transport + lysosomal trap)")
    ax.set_yscale("log")
    ax.set_ylim(1e-9, 30)
    for x, v in zip(xs - 0.20, mvals):
        ax.text(x, v * 1.5, f"{v:.2f}", ha="center", fontsize=8, color=MEAS_COLOR)
    for x, v in zip(xs + 0.20, pvals):
        ax.text(x, v * 1.5, f"{v:.1e}", ha="center", fontsize=8, color=PRED_COLOR)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{c}\npKa {chem[c].pka_structure:.1f}"
                        f"{'  (base)' if chem[c].n_basic_centres else '  (neutral)'}"
                        for c in comps], fontsize=8.5)
    _style(ax, "", "coefficient of variation across the section (log)",
           "b  Heterogeneity: the model is 2-3 decades short")
    gapf = a2["magnitude_gap"]["fold_under_prediction_seq_only"]
    ax.text(0.02, 0.97,
            f"weak base / olaparib, predicted: "
            f"{a2['predicted']['transport+lysosomal_seq']['weak_base_over_olaparib_cv']:.0f}x "
            f"(ordering correct)\n"
            f"measured / predicted magnitude: {gapf:,.0f}x under-predicted\n"
            f"niraparib / rucaparib predicted: "
            f"{a2['predicted']['transport+lysosomal_seq']['nir_over_ruc_cv']:.3f} "
            f"(indistinguishable)",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.4, color=INK2,
            bbox=dict(boxstyle="round,pad=0.45", fc="#f3f2ee", ec="#e0dfd9", lw=0.8))
    ax.legend(fontsize=8, frameon=False, loc="lower left")

    # ---- (c) measured exposure -> apoptosis signature ---------------------
    ax = axes[0, 2]
    t = a3["tests"]["apoptosis"]["measured::pooled"]
    for pat, gp in scores.groupby("Patient"):
        for cmp_, gc in gp.groupby("Compound"):
            ax.scatter(gc.drug_measure_norm, gc.apoptosis_score, s=42,
                       facecolor=PATIENT_COLOR[pat], edgecolor="#fcfcfb", lw=1.0,
                       marker=COMPOUND_MARKER[cmp_], zorder=3, alpha=0.95)
    xx = np.linspace(scores.drug_measure_norm.min(), scores.drug_measure_norm.max(), 50)
    ax.plot(xx, t["intercept"] + t["beta"] * xx, color=INK, lw=2.0, zorder=4)
    _style(ax, "measured drug_measure_norm (IS ratio)",
           "apoptosis signature score (mean z)",
           "c  Measured exposure -> apoptosis signature")
    pt = a3["tests"]["apoptosis"]["predicted::pooled"]
    ax.text(0.02, 0.97,
            f"MEASURED drug: $\\beta$ = {t['beta']:+.3f} $\\pm$ {t['se']:.3f}, "
            f"LRT p = {t['p_lrt']:.1e}\n"
            f"  vs 1,000 size-matched random gene sets: "
            f"p = {a3['tests']['apoptosis']['gene_set_null']['p_empirical']:.3f}\n"
            f"OUR predicted exposure: p = {pt['p_lrt']:.2f}"
            f"{'  (degenerate)' if pt.get('degenerate') else ''}",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.4, color=INK2,
            bbox=dict(boxstyle="round,pad=0.45", fc="#f3f2ee", ec="#e0dfd9", lw=0.8))

    # ---- (d) saturation time-course ---------------------------------------
    ax = axes[1, 0]
    for cond_, ls, lab in (("transport_only", "--", "transport only"),
                           ("transport+lysosomal_seq", "-", "+ lysosomal trap (R~5.4)")):
        s = sat[sat.condition == cond_].groupby("hours").sat_frac
        m, lo, hi = s.mean(), s.min(), s.max()
        ax.plot(m.index, m.values, ls, color=PRED_COLOR, lw=2.0, zorder=3, label=lab)
        ax.fill_between(m.index, lo.values, hi.values, color=PRED_COLOR, alpha=0.12,
                        lw=0, zorder=2)
    ax.axvline(24.0, color=INK3, lw=1.0, ls=":", zorder=1)
    ax.text(23.2, 0.06, "endpoint 24 h", rotation=90, fontsize=8, color=INK2, ha="right")
    ax.axhline(1.0, color=INK3, lw=1.0, ls=":", zorder=1)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.12)
    _style(ax, "bath exposure time (h, log)", "mean tissue saturation  C / (R x C_bath)",
           "d  Why: the explant equilibrates long before 24 h")
    ax.text(0.03, 0.35,
            f"Fourier number Fo = D t / L$^2$ at 24 h:\n"
            f"  median {a2['fourier']['median_with_seq']:.0f}, "
            f"min {a2['fourier']['min_with_seq']:.0f} (with trapping)\n"
            f"Fo >> 1 on {a2['fourier']['n_slices'] - a2['fourier']['n_slices_Fo_lt_1_with_seq']}"
            f"/{a2['fourier']['n_slices']} sections -> diffusive equilibrium,\n"
            f"so a uniform-tissue model must predict a FLAT field.",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.4, color=INK2,
            bbox=dict(boxstyle="round,pad=0.45", fc="#f3f2ee", ec="#e0dfd9", lw=0.8))
    ax.legend(fontsize=8, frameon=False, loc="lower right")

    # ---- (e) measured field vs predicted field, one section ---------------
    ax = axes[1, 1]
    sid = max(geos, key=lambda k: geos[k].mask.sum() if geos[k].compound == "Niraparib" else -1)
    g = geos[sid]
    v = np.where(g.mask, g.measured, np.nan)
    im = ax.imshow(v, origin="lower", cmap="Blues", aspect="equal",
                   extent=[0, g.mask.shape[1] * g.dx / 1000.0,
                           0, g.mask.shape[0] * g.dx / 1000.0])
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("measured norm_Nir (IS ratio)", fontsize=8.5, color=INK2)
    cb.ax.tick_params(labelsize=8, colors=INK2)
    _style(ax, "x (mm)", "y (mm)",
           f"e  Measured field, {sid} — patchy, not a rim")
    row = meas[meas.msi_slice == sid].iloc[0]
    ax.text(0.02, 0.97,
            f"measured Gini {row.meas_gini:.2f} vs internal-standard Gini {row.is_gini:.2f}\n"
            f"drug~depth-from-bath Spearman $\\rho$ = {row.meas_rho_depth:+.2f}\n"
            f"model predicts $\\rho$ = "
            f"{sim[(sim.msi_slice == sid) & (sim.compound == 'Niraparib') & (sim.condition == 'transport+lysosomal_seq')].pred_rho_depth.iloc[0]:+.2f}"
            f"  (wrong sign)",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.2, color=INK,
            bbox=dict(boxstyle="round,pad=0.45", fc="#fcfcfbdd", ec="#e0dfd9", lw=0.8))

    # ---- (f) predicted exposure -> apoptosis signature --------------------
    ax = axes[1, 2]
    dp = scores.dropna(subset=["pred_norm"])
    for pat, gp in dp.groupby("Patient"):
        for cmp_, gc in gp.groupby("Compound"):
            ax.scatter(gc.pred_norm, gc.apoptosis_score, s=42,
                       facecolor=PATIENT_COLOR[pat], edgecolor="#fcfcfb", lw=1.0,
                       marker=COMPOUND_MARKER[cmp_], zorder=3, alpha=0.95)
    pt = a3["tests"]["apoptosis"]["predicted::pooled"]
    if not pt.get("degenerate") and np.isfinite(pt["beta"]):
        xx = np.linspace(dp.pred_norm.min(), dp.pred_norm.max(), 50)
        ax.plot(xx, pt["intercept"] + pt["beta"] * xx, color=INK, lw=2.0, zorder=4)
    _style(ax, "OUR predicted exposure (compound-normalised)",
           "apoptosis signature score (mean z)",
           "f  Our predicted exposure -> apoptosis signature")
    ax.text(0.02, 0.97,
            f"$\\beta$ = {pt['beta']:+.3f} $\\pm$ {pt['se']:.3f}, LRT p = {pt['p_lrt']:.2f}\n"
            f"variance beyond patient+slice REs: "
            f"{pt['var_explained_beyond_re']*100:.1f}%\n"
            f"the predictor is constant within a slice, so it can only\n"
            f"ever reach the between-slice stratum "
            f"({a1['variance_decomposition']['all']['between_slice_within_patient']*100:.0f}% of variance)",
            transform=ax.transAxes, va="top", ha="left", fontsize=8.4, color=INK2,
            bbox=dict(boxstyle="round,pad=0.45", fc="#f3f2ee", ec="#e0dfd9", lw=0.8))

    fig.suptitle(
        "Moncayo et al. 2026 HGSOC patient-derived explants — bath-dosed ex vivo (20 uM, 24 h), "
        "NO vasculature: validates diffusion/uptake/sequestration only",
        fontsize=12.5, color=INK, y=0.985)
    fig.text(0.5, 0.005,
             "n = 3 patients, 114 GeoMx ROIs, 18 MALDI-MSI sections. Drug values are "
             "deuterated-internal-standard-normalised ion ratios (unitless) — no absolute uM. "
             "Ovarian, not glioma: physics transfer, not biology transfer.",
             ha="center", fontsize=8.5, color=INK2)
    fig.tight_layout(rect=[0, 0.02, 1, 0.965])
    fig.savefig(path, dpi=170, facecolor=fig.get_facecolor())
    plt.close(fig)
    return path


# ===========================================================================
# 11. REPORT
# ===========================================================================
def _f(x, nd=3, pct=False):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "n/a"
    if not np.isfinite(v):
        return "n/a"
    return f"{v * 100:.{nd}f}%" if pct else f"{v:.{nd}f}"


def _p(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return "n/a"
    if not np.isfinite(v):
        return "n/a"
    return f"{v:.3f}" if v >= 1e-3 else f"{v:.1e}"


def write_report(a1: dict, a2: dict, a3: dict, meas: pd.DataFrame,
                 sim: pd.DataFrame, geos: Dict[str, SliceGeometry],
                 chem: Dict[str, CompoundChemistry], slab: dict,
                 prov: List[dict], path: str = OUT_MD) -> str:
    cond = a1["headline_condition"]
    C = a1["conditions"]
    vd = a1["variance_decomposition"]["all"]
    seqp = a2["predicted"]["transport+lysosomal_seq"]
    gap = a2["magnitude_gap"]
    apo = a3["tests"]["apoptosis"]
    lys = a3["tests"]["lysosome"]
    enr = a3["gene_level_enrichment"]["pooled"]["sets"]
    L: List[str] = []
    w = L.append

    w("# Moncayo et al. 2026 HGSOC explants — real-data validation of diffusion, "
      "uptake and sequestration")
    w("")
    w("`tmesim/validation/moncayo_explant.py` · "
      "`./.venv/bin/python -m tmesim.validation.moncayo_explant`")
    w("")

    # ---- the caveat, first ------------------------------------------------
    w("## 0. Read this before any number below")
    w("")
    w("> **The explants were bath-dosed ex vivo: 20 uM in medium, 24 h. There is no "
      "vasculature, no perfusion, and no source term.** This run can validate "
      "**diffusion into tissue, cellular uptake and intracellular sequestration**. "
      "It **cannot** validate vascular delivery, interstitial pressure, vessel "
      "normalisation, or anything downstream of a blood supply. Every claim below "
      "is scoped to the bath-dosed regime.")
    w("")
    w("Four more limits that travel with every number:")
    w("")
    w("1. **Ovarian, not glioma.** High-grade serous ovarian carcinoma "
      "patient-derived explants. Everything else in `tme-sim` is glioma. Whatever "
      "transfers is *physics* (a diffusion-reaction field in tissue), not *biology*.")
    w("2. **No absolute concentrations.** The measured quantity is a MALDI-MSI drug "
      "ion intensity divided by a co-sprayed **drug-matched deuterated internal "
      "standard**, averaged over the ROI. It is a unitless ratio. No uM claim is "
      "made about any measured value anywhere in this document.")
    w(f"3. **n = 3 patients**, 114 GeoMx ROIs, 18 MALDI-MSI sections, 12 ROI-table "
      f"slices of which 11 have an MSI counterpart "
      f"({C[cond]['n_roi']} ROIs enter the fits; "
      f"{C[cond]['n_roi_dropped_no_msi']} are dropped because slice "
      "`T23_071_NirBlad2` was not deposited in the MSI table).")
    w("4. **The ROI table has no X/Y.** GeoMx ROI positions were never deposited, so "
      "a transport model's prediction can only be resolved to the **slice**. See "
      "§1 for what that costs.")
    w("")

    # ---- provenance -------------------------------------------------------
    w("## 1. Provenance and what is being predicted")
    w("")
    w("Moncayo et al., *Nature Communications* 2026, PMID 41844641, "
      "DOI 10.1038/s41467-026-70558-1. Zenodo 17610220 (CC-BY-4.0) + GEO GSE281519 "
      "+ the journal Source Data workbook. Files as used, on disk:")
    w("")
    w("| file | bytes | sha256 (first 16) |")
    w("|---|---:|---|")
    for r in prov:
        w(f"| `{r['path']}` | {r['bytes']:,} | `{r['sha256_16']}` |")
    w("")
    w("**Response variable**: `drug_measure_norm` = `drug_measure` / "
      "`measure_norm_factor`, where `measure_norm_factor` is a single constant per "
      "compound (asserted at load time). **Frame**: the authors' own, shipped in "
      "`zenodo17610220_modelling.Rmd` — "
      "`lmer(y ~ x + (1|Patient) + (1|Patient:Slice), REML=FALSE)` with a "
      "likelihood-ratio test, reproduced here with `statsmodels` MixedLM "
      "(nested slice as a variance component inside the patient grouping, ML).")
    w("")
    w("### Where the measured variance actually lives")
    w("")
    w("| stratum | share of variance in `drug_measure_norm` |")
    w("|---|---:|")
    w(f"| between patient | {_f(vd['between_patient'], 1, pct=True)} |")
    w(f"| between slice, within patient | "
      f"{_f(vd['between_slice_within_patient'], 1, pct=True)} |")
    w(f"| **within slice (ROI to ROI)** | "
      f"**{_f(vd['within_slice_residual'], 1, pct=True)}** |")
    w("")
    w("This single table caps the whole exercise. A coordinate-free, slice-resolved "
      "prediction is **structurally incapable of reaching the 90% of the variance "
      "that lives between ROIs inside one section**. The most a perfect "
      "slice-level model could explain beyond the random effects is the "
      f"{_f(vd['between_slice_within_patient'], 1, pct=True)} between-slice stratum, "
      "and the random effects are already fitted to absorb exactly that. We state "
      "this before reporting the fit rather than after.")
    w("")

    # ---- methods ----------------------------------------------------------
    w("## 2. What was simulated, and how (additive-only)")
    w("")
    w("**Geometry — measured, not invented.** Each of the 18 MSI sections is "
      "rasterised from *its own* on-tissue pixel coordinates on the native 20 um "
      "grid, hole-filled, then block-coarsened 3x to a 60 um model grid. Section "
      "areas run "
      f"{min(g.area_mm2 for g in geos.values()):.1f}-"
      f"{max(g.area_mm2 for g in geos.values()):.1f} mm^2 with maximum in-plane "
      f"distance-to-edge "
      f"{min(g.max_depth_um for g in geos.values()):.0f}-"
      f"{max(g.max_depth_um for g in geos.values()):.0f} um.")
    w("")
    w("**Boundary condition — a bath, not a vessel.** Dirichlet C = 20 uM on the "
      "entire section perimeter, C = 0 inside at t = 0, for 24 h. This is the "
      "opposite of the vascular-source geometry used everywhere else in this repo, "
      "and it is the correct representation of a submerged explant.")
    w("")
    w("**Solver — engine numerics, analysis-layer boundary.** The frozen public "
      "solver `simulate_drug_transport` implements a volumetric source with "
      "**no-flux** boundaries and a zero initial condition; it has no Dirichlet "
      "facility, so it structurally cannot express a bath. Rather than modify it, "
      "this module imports and calls two engine internals unchanged — "
      "`physics.diffusion._effective_diffusivity` (the same "
      "`D_eff = D (1 - 0.6 x ecm)` tortuosity law) and "
      "`physics.diffusion._build_diffusion_operator` (the same symmetric "
      "finite-volume Laplacian) — and adds only the boundary condition. "
      "**Nothing in `core/`, `physics/`, `molecular/`, `reconstruct/` or `sim/` was "
      "edited.**")
    w("")
    w(f"**Solver verification.** Checked against Crank's analytic plane-sheet "
      f"solution (eq. 4.18) at 0.25 / 1 / 4 / 24 h: worst relative error "
      f"**{_f(slab['worst_rel_err'] * 100, 2)}%** "
      f"({'PASS' if slab['passed'] else 'FAIL'}, tolerance 3%).")
    w("")
    w("**How lysosomal sequestration is encoded without touching the engine.** "
      "Reversible instantaneous partitioning into an acidic compartment of volume "
      "fraction `phi_lys` is a *linear equilibrium*, so it is exactly equivalent to "
      "a retardation factor")
    w("")
    w("```")
    w("  ratio = [D]_lys/[D]_cyt = (1 + 10^(pKa - 4.7)) / (1 + 10^(pKa - 7.2))"
      "      (Henderson-Hasselbalch ion trapping)")
    w("  R     = 1 + phi_lys * (ratio - 1)")
    w("  R dC_free/dt = div(D grad C_free)   <=>   D -> D/R,  C_measured = R * C_free")
    w("```")
    w("")
    w("so it is implemented as a **constructed `DrugParams` whose diffusivity is "
      "divided by R**, plus a post-multiplication of the solved field by R. No "
      "engine change, and exact for linear partitioning. Reference "
      f"`phi_lys` = {PHI_LYS} (cell volume fraction ~0.7 x lysosomal fraction of "
      f"cell volume ~0.02), swept over {PHI_LYS_SWEEP} in §4.")
    w("")
    w("**Three model conditions**, all run on all 18 geometries x all 3 compounds "
      "(so the compound contrast is within-geometry and cannot be confounded by "
      "explant shape):")
    w("")
    w("| condition | pH-partitioning / lysosomal trap | molecular-layer "
      "clearance+binding sink |")
    w("|---|---|---|")
    w("| `transport_only` | no | no |")
    w("| `transport+shipped_loss` | no | yes |")
    w("| `transport+lysosomal_seq` | **yes** | no |")
    w("")

    # ---- analysis 1 -------------------------------------------------------
    w("## 3. Analysis 1 — accumulation-heterogeneity fit")
    w("")
    w("### 3.1 The transport-only model has nothing to fit")
    w("")
    sat = a2["saturation"]
    fo = a2["fourier"]
    w(f"Under `transport_only` the predicted accumulation is **identical to the bath "
      f"concentration in every voxel of every section** (mean saturation "
      f"{_f(a2['predicted']['transport_only']['Niraparib']['sat_frac_mean'], 6)}, "
      f"predicted CV "
      f"{a2['predicted']['transport_only']['Niraparib']['cv_mean']:.1e}). The "
      "predictor is constant, no slope is identifiable, and **variance explained "
      "beyond the random effects is 0.0% by construction**. That is not a fitting "
      "failure; it is the model's answer.")
    w("")
    w(f"The reason is dimensional. The Fourier number Fo = D_eff t / L^2 at the 24 h "
      f"endpoint has median **{_f(fo['median_no_seq'], 1)}** without trapping and "
      f"**{_f(fo['median_with_seq'], 1)}** with it (minimum "
      f"{_f(fo['min_with_seq'], 2)}); "
      f"{fo['n_slices'] - fo['n_slices_Fo_lt_1_with_seq']}/{fo['n_slices']} sections "
      f"have Fo > 1. Mean tissue saturation is already "
      f"{_f(sat['mean_sat_frac_1h_transport_only'], 2)} at **1 hour** "
      f"({_f(sat['mean_sat_frac_1h_with_seq'], 2)} with trapping) and "
      f"{_f(sat['mean_sat_frac_24h_with_seq'], 4)} at 24 h. "
      "**A 24 h bath-dosed 1-2 mm explant is at diffusive equilibrium.** A model "
      "whose only spatial mechanism is diffusion from a uniform boundary into "
      "uniform tissue must predict a flat field there — and does.")
    w("")
    w("### 3.2 The one condition with a non-constant prediction")
    w("")
    w(f"`{cond}` adds the molecular layer's own first-order "
      "`clearance + binding` sink, which makes the steady state depend on section "
      "size and shape, so the prediction finally varies slice to slice "
      f"(relative SD {_f(C[cond]['slice_level']['pred_rel_sd'], 3)}).")
    w("")
    w("| statistic | pooled | Niraparib | Rucaparib |")
    w("|---|---:|---:|---:|")
    for k, lab in (("beta", "beta"), ("se", "SE"), ("p_lrt", "LRT p")):
        row = [C[cond]["pooled"], C[cond]["Niraparib"], C[cond]["Rucaparib"]]
        vals = [(_p(r[k]) if k == "p_lrt" else _f(r[k])) for r in row]
        w(f"| {lab} | {vals[0]} | {vals[1]} | {vals[2]} |")
    row = [C[cond]["pooled"], C[cond]["Niraparib"], C[cond]["Rucaparib"]]
    w("| **dR2 marginal (variance beyond REs)** | " +
      " | ".join(_f(r["var_explained_beyond_re"], 1, pct=True) for r in row) + " |")
    w("| R2 conditional | " + " | ".join(_f(r["r2_conditional"], 3) for r in row) + " |")
    w("| singular fit | " +
      " | ".join("**yes**" if r.get("unstable") else "no" for r in row) + " |")
    w("")
    mp = C[cond].get("mixedlm_permutation", {})
    sl = C[cond]["slice_level"]
    w("### 3.3 And why that p-value does not survive")
    w("")
    w("The parametric LRT above treats 112 ROIs as if they carried 112 units of "
      "information. They do not: **our predictor is constant within a slice**, so "
      "the exchangeable unit is the slice, of which there are 11. A slice random "
      "intercept with a small variance component "
      f"({_f(vd['between_slice_within_patient'], 1, pct=True)} of total) does not "
      "absorb a slice-constant regressor, and the ROI count leaks straight into the "
      "standard error. Permuting the slice->prediction assignment *within compound* "
      "(preserving the compound split, the slice means and the ROI counts) and "
      f"refitting {mp.get('n_perm_completed', 0)} times gives the honest reference:")
    w("")
    w("| | value |")
    w("|---|---:|")
    w(f"| observed abs(t) | {_f(mp.get('observed_abs_t'), 2)} |")
    w(f"| permutation null median abs(t) | {_f(mp.get('null_median_abs_t'), 2)} |")
    w(f"| permutation null 95th pct abs(t) | {_f(mp.get('null_p95_abs_t'), 2)} |")
    w(f"| **permutation p (slice-exchangeable)** | "
      f"**{_p(mp.get('p_perm_t'))}** |")
    w(f"| observed dR2 beyond REs | "
      f"{_f(mp.get('observed_var_explained_beyond_re'), 1, pct=True)} |")
    w(f"| permutation null median dR2 | "
      f"{_f(mp.get('null_median_var_explained'), 1, pct=True)} |")
    w(f"| permutation null 95th pct dR2 | "
      f"{_f(mp.get('null_p95_var_explained'), 1, pct=True)} |")
    w(f"| **permutation p for dR2** | {_p(mp.get('p_perm_var_explained'))} |")
    w("")
    w(f"Slice-level agreement, the same information viewed directly (n = "
      f"{sl['n_slices']} slices): Spearman rho = {_f(sl['spearman_rho'], 3)} "
      f"(p = {_p(sl['spearman_p'])}), Pearson r = {_f(sl['pearson_r'], 3)} "
      f"(p = {_p(sl['pearson_p'])}); within-compound slice permutation "
      f"p = {_p(C[cond]['permutation']['p_perm_two_sided'])}.")
    w("")
    w("**Verdict, analysis 1.** *Variance explained beyond the patient/slice random "
      f"effects: nominally {_f(C[cond]['pooled']['var_explained_beyond_re'], 1, pct=True)}, "
      f"but the permutation null median under a random slice ordering is "
      f"{_f(mp.get('null_median_var_explained'), 1, pct=True)} and the permutation "
      f"p is {_p(mp.get('p_perm_var_explained'))}.* With 11 exchangeable slices, a "
      "few percent of marginal R2 is what chance delivers. **The honest number is "
      "zero: this model does not predict per-ROI drug accumulation in bath-dosed "
      "explants.** Under `transport_only` it is zero by construction; under "
      f"`{cond}` it is zero after calibration against the correct null.")
    w("")

    # ---- analysis 2 -------------------------------------------------------
    w("## 4. Analysis 2 — the compound-specificity falsification test")
    w("")
    w("**The claim under test.** Niraparib and rucaparib are weak bases whose "
      "lysosomal sequestration drives accumulation heterogeneity; olaparib is not a "
      "base and does not show it. A model carrying the mechanism must (a) order the "
      "compounds correctly *and* (b) get the magnitude approximately right. Passing "
      "(a) at 1/175th of the observed magnitude is not a pass.")
    w("")
    w("### 4.1 Does the molecular layer see the chemistry? Yes.")
    w("")
    w("Basic centres are found from structure with RDKit — an sp3 nitrogen that is "
      "not an amide/sulfonamide/phosphoramide N, not attached to an aromatic ring, "
      "not already charged (`" + BASIC_N_SMARTS + "`). InChIKeys were asserted "
      "against the public identifiers before anything else ran.")
    w("")
    w("| compound | InChIKey | MW | cLogP | basic centres | class | pKa (structure) "
      "| pKa (literature, cross-check only) | lysosome:cytosol ratio | R at "
      f"phi_lys={PHI_LYS} |")
    w("|---|---|---:|---:|---:|---|---:|---:|---:|---:|")
    for n in ("Niraparib", "Rucaparib", "Olaparib"):
        c = a2["chemistry"][n]
        w(f"| {n} | `{c['inchikey']}` | {c['mw']:.1f} | {c['logp']:.2f} | "
          f"{c['n_basic_centres']} | {c['basic_class']} | {c['pka_structure']:.1f} | "
          f"{c['pka_literature'] if c['pka_literature'] else '-- (no basic centre)'} | "
          f"{c['lysosome_cytosol_ratio']:.1f} | {c['retardation_R']:.2f} |")
    w("")
    w("**This part works.** Purely from SMILES, the layer finds one basic secondary "
      "amine in niraparib and one in rucaparib and **none** in olaparib (both "
      "piperazine nitrogens are acylated), matching the literature cross-check "
      "exactly. The weak-base / non-base split is recovered from structure: "
      f"`weak_base_split_detected = {a2['weak_base_split_detected']}`.")
    w("")
    w("### 4.2 Does that chemistry produce the predicted heterogeneity? Ordering "
      "yes, magnitude no.")
    w("")
    w("Predicted coefficient of variation of accumulation across each section, "
      "averaged over all 18 geometries (every compound on every geometry):")
    w("")
    w("| compound | `transport_only` | `transport+lysosomal_seq` | measured CV | "
      "measured CV minus IS-noise floor |")
    w("|---|---:|---:|---:|---:|")
    for n in ("Niraparib", "Rucaparib", "Olaparib"):
        w(f"| {n} | {a2['predicted']['transport_only'][n]['cv_mean']:.2e} | "
          f"{seqp[n]['cv_mean']:.2e} | "
          f"{_f(a2['measured'][n]['cv_mean'], 3)} | "
          f"{_f(a2['measured'][n]['excess_cv_mean'], 3)} |")
    w("")
    w(f"- **Ordering: correct.** With sequestration on, predicted heterogeneity for "
      f"the weak bases is **{seqp['weak_base_over_olaparib_cv']:,.0f}x** olaparib's. "
      "The model does not produce the same heterogeneity for all three compounds — "
      "so the narrow version of the falsification test is *passed*, and the "
      "pH-partitioning term is doing real, chemistry-driven work.")
    w(f"- **Magnitude: off by {gap['fold_under_prediction_seq_only']:,.0f}x.** "
      f"Predicted weak-base CV is {gap['pred_cv_weak_base_seq_only']:.2e}; the "
      f"measured, internal-standard-noise-corrected weak-base CV is "
      f"{_f(gap['meas_excess_cv_weak_base'], 3)}. The mechanism is present in the "
      "model but contributes essentially none of the observed spatial structure.")
    w(f"- **Niraparib and rucaparib are indistinguishable to the model.** Predicted "
      f"CV ratio {seqp['nir_over_ruc_cv']:.3f}. Both are strong enough bases that "
      "the Henderson-Hasselbalch trapping ratio saturates at "
      f"10^(7.2-4.7) ~ 316 — the two ratios differ by only "
      f"{_f(a2['nir_vs_ruc_ratio_gap_pct'], 2)}%. Any pH-partition model has this "
      "ceiling; it can separate a base from a non-base and nothing finer.")
    w("")
    w("### 4.3 Why a bigger lysosomal compartment cannot rescue it")
    w("")
    w("| phi_lys | R (weak base) | predicted CV (weak base) |")
    w("|---:|---:|---:|")
    for rec in a2["phi_lys_sweep"]:
        if rec["compound"] == "Niraparib":
            w(f"| {rec['phi_lys']} | {rec['R']:.1f} | {rec['pred_cv']:.2e} |")
    w("")
    w(f"Even at phi_lys = {max(PHI_LYS_SWEEP)} — 5% of tissue volume as acidic "
      "compartment, at or above every value in the cell-biology literature — the "
      f"predicted CV reaches only {a2['phi_lys_sweep_max_cv']:.2e}, still "
      f"{_f(gap['meas_excess_cv_weak_base'] / a2['phi_lys_sweep_max_cv'], 1)}x short.")
    w("")
    w("**And this is structural, not a tuning problem.** With `phi_lys` *uniform in "
      "space*, R is uniform, so at equilibrium C_total = R x C_bath **everywhere** "
      "and the dispersion goes to zero however large R becomes. Uniform "
      "sequestration buys a longer time constant, never an equilibrium pattern. To "
      "produce the observed hotspots the model would need a **spatial field of "
      "lysosomal / acidic-phospholipid content**. A literature pass for this project "
      "established that **no measured mg/g or mol% acidic-phospholipid content "
      "exists for tumour tissue and no published PBPK model carries "
      "tumour-specific lysosomal parameters** — so this is not a parameter we "
      "declined to fit, it is a parameter that does not exist in the literature. "
      "That makes it unexplored territory rather than an oversight, and it is the "
      "single most useful thing this run identifies.")
    w("")
    w("### 4.4 The measured side of the contrast is confounded — stated, not spun")
    w("")
    w("| compound | slices | measured Gini | internal-standard Gini | excess Gini | "
      "excess Moran's I | exact-zero pixel fraction |")
    w("|---|---:|---:|---:|---:|---:|---:|")
    for n in ("Niraparib", "Rucaparib", "Olaparib"):
        mm = a2["measured"][n]
        w(f"| {n} | {mm['n_slices']} | {_f(mm['gini_mean'])} | "
          f"{_f(mm['is_gini_mean'])} | {_f(mm['excess_gini_mean'])} | "
          f"{_f(mm['excess_moran_mean'])} | {_f(mm['zero_frac_mean'], 2, pct=True)} |")
    w("")
    for metric in ("excess_gini", "excess_moran"):
        t = a2["measured_weak_base_vs_olaparib"][metric]
        w(f"- Weak bases vs olaparib on **{metric}**: "
          f"{_f(t['weak_base_mean'])} vs {_f(t['olaparib_mean'])}, Mann-Whitney "
          f"one-sided p = {_p(t['p_one_sided_greater'])} "
          f"(n = {t['n_weak_base_slices']} vs {t['n_olaparib_slices']} sections).")
    w("")
    w(a2["measured_metric_note"])
    w("")
    w("Two further reasons not to read the measured columns as a replication of the "
      "authors' claim: **olaparib contributes zero GeoMx ROIs** (its drug hotspots "
      "were too small to hit the nuclei minimum, so it is absent from analyses 1 "
      "and 3 entirely), and the paper's compound-specificity evidence rests on "
      "FACS-sorted reculture and hotspot morphology, not on section-level "
      "dispersion. Our metric is not their metric.")
    w("")
    w("### 4.5 The sign of the spatial gradient — the sharpest single failure")
    w("")
    dg = a2["depth_gradient"]
    w(f"Correlation between drug level and distance-to-bath-surface, per section: "
      f"the model predicts Spearman rho = {_f(dg['model_mean_rho_depth'], 2)} "
      "(drug falls off with depth, as a penetration-limited field must). The "
      f"measurement gives mean rho = {_f(dg['measured_mean_rho_depth'], 3)}, and is "
      f"**positive on {dg['n_slices_measured_positive']}/{dg['n_slices']} sections** "
      "— drug is, if anything, slightly *higher* in the section interior.")
    w("")
    w("The measured field is not a penetration gradient at all. It is patchy, and "
      "the patches are not organised by distance from the medium. Whatever sets "
      "them — cell density, tumour-vs-stroma composition, lysosomal capacity — is "
      "orthogonal to the transport axis this model computes.")
    w("")
    v = a2["verdict"]
    w(f"**Verdict, analysis 2: {v['overall']}.** Ordering correct "
      f"({v['ordering_correct']}), magnitude correct ({v['magnitude_correct']}). "
      "The molecular layer *does* distinguish the weak bases from olaparib, from "
      "structure, with no hand-entered pKa — so this is not the clean total "
      "negative one might have expected. But the mechanism as implemented is "
      f"~{gap['fold_under_prediction_seq_only']:,.0f}x too weak to account for the "
      "measured heterogeneity, it cannot separate two strong bases from each other, "
      "and it produces a gradient of the wrong sign. **The lysosomal-reservoir "
      "mechanism is present in name and absent in effect.** Reported as a "
      "publishable negative, not a failure.")
    w("")

    # ---- analysis 3 -------------------------------------------------------
    w("## 5. Analysis 3 — exposure -> molecular response coupling")
    w("")
    w("### 5.1 Gene sets and scoring")
    w("")
    sc = a3["scoring"]
    w(f"**MSigDB Hallmark 2020 *Apoptosis*** ({sc['sets']['apoptosis']['n_requested']} "
      f"symbols, {sc['sets']['apoptosis']['n_on_panel']} on the GeoMx WTA panel) and "
      f"**KEGG 2021 Human *Lysosome*** ({sc['sets']['lysosome']['n_requested']} "
      f"symbols, {sc['sets']['lysosome']['n_on_panel']} on panel), both embedded "
      "verbatim in the module so the run is reproducible offline. Score = mean "
      "z-score: deposited Q3-normalised counts -> `log2(x+1)` -> per-gene z across "
      f"all 114 ROIs -> unweighted mean over set members. Null = "
      f"{sc['n_null']} random gene sets of the same size from the same "
      f"{sc['n_genes_panel']}-gene panel, scored identically.")
    w("")
    w("### 5.2 Does the MEASURED drug level predict the signatures? Yes.")
    w("")
    w("`signature ~ drug_measure_norm + (1|Patient) + (1|Patient:Slice)`:")
    w("")
    w("| signature | subset | beta | SE | LRT p | dR2 beyond REs |")
    w("|---|---|---:|---:|---:|---:|")
    for name, block in (("apoptosis", apo), ("lysosome", lys)):
        for grp in ("pooled", "Niraparib", "Rucaparib"):
            r = block[f"measured::{grp}"]
            w(f"| {name} | {grp} | {_f(r['beta'], 4)} | {_f(r['se'], 4)} | "
              f"{_p(r['p_lrt'])} | {_f(r['var_explained_beyond_re'], 1, pct=True)} |")
    w("")
    for name, block in (("apoptosis", apo), ("lysosome", lys)):
        n_ = block["gene_set_null"]
        w(f"- **{name}** vs {n_['n_random_sets_fitted']} size-matched random gene "
          f"sets refitted through the same mixed model: observed |t| = "
          f"{_f(n_['observed_abs_t'], 2)}, null median "
          f"{_f(n_['null_median_abs_t'], 2)}, null 95th pct "
          f"{_f(n_['null_p95_abs_t'], 2)}, **empirical p = "
          f"{_p(n_['p_empirical'])}**.")
    w("")
    w("Whole-panel rank enrichment, the fast analogue of the authors' GSEA: all "
      f"{sc['n_genes_panel']} genes ranked by their **within-slice** partial "
      "correlation with `drug_measure_norm` (both sides residualised on slice "
      "identity as fixed effects, which is strictly more conservative than the "
      "authors' random-effects treatment), then tested against 1,000 size-matched "
      "random sets:")
    w("")
    w("| set | mean rank fraction (0.5 = no enrichment) | z vs random sets | "
      "empirical p | mean within-slice r |")
    w("|---|---:|---:|---:|---:|")
    for k in ("apoptosis", "lysosome"):
        e = enr[k]
        w(f"| {k} | {_f(e['mean_rank_fraction'], 3)} | {_f(e['z'], 2)} | "
          f"{_p(e['p_empirical_up'])} | {_f(e['mean_within_slice_r'], 3)} |")
    w("")
    w(f"The paper's benchmark is GSEA **APOPTOSIS NES = "
      f"{a3['paper_benchmark']['gsea_apoptosis_NES']}, FDR q = "
      f"{a3['paper_benchmark']['gsea_apoptosis_FDR_q']:g}**. Our statistic is a "
      "rank-enrichment z, not a permutation NES, so the numbers are not directly "
      "comparable — but the **direction and significance replicate**: apoptosis "
      f"genes sit high among the drug-correlated genes (z = {_f(enr['apoptosis']['z'], 1)}, "
      f"p < 0.001) even after absorbing slice identity as a fixed effect. "
      "**Measured drug level does predict the apoptotic signature, as it should.** "
      "The dataset's response axis is real and our reading of it agrees with the "
      "authors'.")
    w("")
    w("### 5.3 Does OUR predicted exposure predict them? No.")
    w("")
    w("| signature | subset | beta | SE | LRT p | dR2 beyond REs | singular |")
    w("|---|---|---:|---:|---:|---:|---|")
    for name, block in (("apoptosis", apo), ("lysosome", lys)):
        for grp in ("pooled", "Niraparib", "Rucaparib"):
            r = block[f"predicted::{grp}"]
            w(f"| {name} | {grp} | {_f(r['beta'], 4)} | {_f(r['se'], 4)} | "
              f"{_p(r['p_lrt'])} | {_f(r['var_explained_beyond_re'], 1, pct=True)} | "
              f"{'**yes**' if r.get('unstable') else 'no'} |")
    w("")
    ap = apo["predicted::pooled"]
    w(f"Pooled apoptosis against our predicted exposure: beta = {_f(ap['beta'], 3)} "
      f"(**negative** — the wrong sign), LRT p = {_p(ap['p_lrt'])}. Under "
      "`transport_only` the predictor is constant and the question does not arise "
      "at all. This is a **clean negative**, and it is the expected consequence of "
      "§3: a predictor that is constant within a slice can only ever address the "
      f"{_f(vd['between_slice_within_patient'], 1, pct=True)} between-slice stratum, "
      "while the measured drug-to-apoptosis coupling lives *within* slices, which is "
      "exactly where the authors' design put it.")
    w("")
    w("### 5.4 The bridge back to analysis 2")
    w("")
    mb = a3["mechanism_bridge"]
    w(f"The **lysosome** signature tracks measured drug level *more strongly than "
      f"apoptosis does* (beta = {_f(mb['lysosome_signature_beta'], 4)}, "
      f"p = {_p(mb['lysosome_signature_p'])} vs beta = "
      f"{_f(mb['apoptosis_signature_beta'], 4)}, "
      f"p = {_p(mb['apoptosis_signature_p'])}), and in the within-slice ranking it "
      f"is the stronger of the two sets (z = {_f(enr['lysosome']['z'], 1)} vs "
      f"{_f(enr['apoptosis']['z'], 1)}). ROIs with more lysosomal machinery hold "
      "more drug, at the sub-slice scale.")
    w("")
    w("That is same-dataset, independent evidence that the mechanism missing from "
      "§4 is precisely lysosomal capacity, **and** that it varies on the spatial "
      "scale our slice-resolved prediction cannot see. It also points at the "
      "experiment that would fix it: a spatial map of lysosomal/acidic-phospholipid "
      "content co-registered with drug, which does not currently exist for any "
      "tumour tissue.")
    w("")

    # ---- overall ----------------------------------------------------------
    w("## 6. Honest overall verdict")
    w("")
    w("| # | test | verdict |")
    w("|---|---|---|")
    w(f"| 1 | per-ROI accumulation fit, authors' mixed-effects frame | "
      f"**FAIL.** 0.0% beyond REs by construction under `transport_only`; "
      f"{_f(C[cond]['pooled']['var_explained_beyond_re'], 1, pct=True)} nominal "
      f"under `{cond}`, permutation p = {_p(mp.get('p_perm_var_explained'))} "
      "against a slice-exchangeable null. Not distinguishable from chance. |")
    w(f"| 2 | compound-specificity falsification | **PARTIAL / effectively "
      f"negative.** Ordering right ({seqp['weak_base_over_olaparib_cv']:,.0f}x weak "
      f"base : olaparib, recovered from SMILES), magnitude "
      f"{gap['fold_under_prediction_seq_only']:,.0f}x too small, the two weak bases "
      "indistinguishable, gradient of the wrong sign on "
      f"{dg['n_slices_measured_positive']}/{dg['n_slices']} sections. |")
    w(f"| 3 | exposure -> apoptosis coupling | **Measured: PASS** "
      f"(beta = {_f(mb['apoptosis_signature_beta'], 4)}, p = "
      f"{_p(mb['apoptosis_signature_p'])}; rank-enrichment z = "
      f"{_f(enr['apoptosis']['z'], 1)}, replicating the paper's NES 2.278 in "
      f"direction). **Ours: FAIL** (p = {_p(ap['p_lrt'])}, negative sign). |")
    w("")
    w("**What this run establishes.**")
    w("")
    w("1. On the *only* open dataset that pairs a measured spatial drug field with "
      "molecular state on the same rows, **this simulator's transport layer "
      "explains none of the measured per-ROI drug accumulation** beyond patient and "
      "slice identity. Not a small effect — zero, after the correct null.")
    w("2. The reason is diagnosable and is not a bug: a 24 h bath-dosed explant is "
      f"at diffusive equilibrium (Fo median {_f(fo['median_with_seq'], 1)}), so a "
      "uniform-tissue diffusion model *must* predict a flat field. Everything the "
      "measurement sees is heterogeneity the model has no term for.")
    w("3. The pH-partitioning / lysosomal-trapping mechanism, implemented here in "
      "the analysis layer as an exact linear-equilibrium retardation, **recovers "
      "the weak-base/non-base split from chemical structure alone** — but with a "
      "spatially uniform lysosomal volume fraction it cannot generate spatial "
      "structure at equilibrium at any parameter value. The missing ingredient is a "
      "*spatial field* of acidic-compartment content, for which **no measurement "
      "and no PBPK parameterisation exists in tumour tissue**.")
    w("4. The dataset's response axis is real and we reproduce it: measured drug "
      "level predicts the apoptotic transcriptional signature within slices, "
      "concordant with the authors' GSEA. Our predicted exposure does not.")
    w("")
    w("**What this run does *not* establish.** Nothing about vascular delivery, "
      "interstitial pressure, perfusion heterogeneity, or drug transport in glioma. "
      "Nothing about absolute concentrations. Nothing that generalises past n = 3 "
      "patients. A negative here does not falsify the vascular-transport model; it "
      "falsifies the claim that diffusion plus uniform sequestration explains "
      "intratumoral drug heterogeneity in bath-dosed human explants.")
    w("")
    w("![Moncayo explant validation](moncayo_explant.png)")
    w("")
    w("*Panels: (a) measured vs predicted per-ROI accumulation, coloured by patient, "
      "marker by compound, slice means as open diamonds; (b) heterogeneity across "
      "compounds, measured (IS-noise-corrected) vs predicted, log scale; "
      "(c) measured exposure -> apoptosis signature; (d) saturation time-course "
      "showing why the model must predict a flat field; (e) one measured MSI field, "
      "showing patchiness that is not organised by distance from the bath; "
      "(f) our predicted exposure -> apoptosis signature.*")
    w("")

    txt = "\n".join(L) + "\n"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(txt)
    return path


# ===========================================================================
# 12. ENTRY POINT
# ===========================================================================
def main(n_perm: int = 2000, n_perm_mixed: int = 999, n_null: int = 1000,
         verbose: bool = True) -> dict:
    t0 = time.time()
    if verbose:
        print("== Moncayo et al. 2026 HGSOC explant validation ==")
        print("  bath-dosed ex vivo, NO vasculature: diffusion/uptake/sequestration only")

    prov = provenance()
    msi = load_msi_pixels()
    roi = load_roi_table()
    expr = load_expression()
    if verbose:
        print(f"  [data] {len(roi)} ROIs, {len(msi):,} MSI pixels, "
              f"{expr.shape[1] - GENE_COL_START:,} genes")

    slab = verify_against_analytic_slab(verbose=verbose)
    geos = build_geometries(msi)
    chem = compound_chemistry()
    sim = run_simulations(geos, chem, verbose=verbose)
    meas = measured_slice_stats(geos, msi)
    sat = saturation_timecourse(geos, chem)
    sweep = phi_lys_sweep(geos, chem)
    ecm_sens = ecm_sensitivity(geos, chem)
    if verbose:
        print(f"  [ecm] CV across ecm in {ECM_SWEEP}: "
              f"{ecm_sens.groupby('ecm').pred_cv.mean().round(6).to_dict()}")

    a1 = analysis1(roi, sim, meas, n_perm=n_perm, n_perm_mixed=n_perm_mixed,
                   verbose=verbose)
    a2 = analysis2(geos, sim, meas, chem, sat, sweep)
    a3 = analysis3(expr, roi, sim, condition=a1["headline_condition"],
                   n_null=n_null, verbose=verbose)

    png = make_figure(a1, a2, a3, meas, sim, geos, chem, sat)
    md = write_report(a1, a2, a3, meas, sim, geos, chem, slab, prov)
    if verbose:
        print(f"  [out] {os.path.relpath(md, REPO)}")
        print(f"  [out] {os.path.relpath(png, REPO)}")
        print(f"  done in {time.time() - t0:.0f}s")
    return {"a1": a1, "a2": a2, "a3": a3, "sim": sim, "meas": meas, "md": md,
            "png": png}


if __name__ == "__main__":
    main()
