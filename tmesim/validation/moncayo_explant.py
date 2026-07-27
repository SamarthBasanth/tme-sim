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
    geos: Dict[str, SliceGeometry] = {}
    for sl, g in msi.groupby("Slice"):
        compound = str(g.Drug.iloc[0])
        chan = MSI_CHANNEL[compound][0]
        fine_mask, fine_val = _rasterise(g, chan)
        fine_mask = ndimage.binary_fill_holes(fine_mask)
        mask, val = _coarsen(fine_mask, fine_val, COARSEN)
        # keep the largest connected component (drops rasterisation specks)
        lab, n = ndimage.label(mask)
        if n > 1:
            sizes = ndimage.sum(mask, lab, range(1, n + 1))
            mask = lab == (int(np.argmax(sizes)) + 1)
            val = val * mask
        depth = ndimage.distance_transform_edt(np.pad(mask, 1))[1:-1, 1:-1] * MODEL_PITCH_UM
        geos[sl] = SliceGeometry(
            slice_id=sl, compound=compound, patient=str(g.Patient.iloc[0]),
            tissue=str(g.Tissue.iloc[0]), mask=mask, dx=MODEL_PITCH_UM,
            depth_um=depth, measured=val, n_px_native=int(len(g)),
            area_mm2=float(fine_mask.sum()) * (MSI_PITCH_MM ** 2),
            max_depth_um=float(depth[mask].max()),
            mean_depth_um=float(depth[mask].mean()),
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
                    "pred_cv": float(v.std() / max(v.mean(), 1e-30)),
                    "pred_p10_p90": float(np.percentile(v, 90) / max(np.percentile(v, 10), 1e-30)),
                    "pred_auc_mean": float(a.mean()),
                    "pred_sat_frac": float(v.mean() / (res.R * BATH_CONC_UM)),
                    "pred_rho_depth": float(rho),
                    "R": res.R, "d_eff_mean": res.d_eff_mean,
                })
    if verbose:
        print(f"  [sim] {len(rows)} runs in {time.time() - t0:.1f}s (ecm={ecm})")
    return pd.DataFrame(rows)


def measured_slice_stats(geos: Dict[str, SliceGeometry],
                         msi: pd.DataFrame) -> pd.DataFrame:
    """Measured per-slice heterogeneity on the SAME 60 um grid as the model.

    Coarsening from 20 um to 60 um is itself the noise control: MALDI per-pixel
    counts are shot-noise dominated (olaparib has >50% exact-zero pixels), and
    raw per-pixel CV rewards noise. Block-averaging to the model's own scale
    makes the measured and predicted CVs commensurable.
    """
    rows = []
    for sid, g in geos.items():
        v = g.measured[g.mask]
        sub = msi[msi.Slice == sid]
        chan, raw, heavy = MSI_CHANNEL[g.compound]
        rho, _ = spearmanr(g.depth_um[g.mask], v)
        rows.append({
            "msi_slice": sid, "compound": g.compound, "patient": g.patient,
            "tissue": g.tissue, "area_mm2": g.area_mm2,
            "max_depth_um": g.max_depth_um, "mean_depth_um": g.mean_depth_um,
            "meas_mean": float(v.mean()),
            "meas_cv": float(v.std() / max(v.mean(), 1e-30)),
            "meas_p10_p90": float(np.percentile(v, 90) / max(np.percentile(v, 10), 1e-30)),
            "meas_rho_depth": float(rho),
            "meas_rho_depth_p": float(spearmanr(g.depth_um[g.mask], v)[1]),
            "zero_frac_native": float((sub[chan].values <= 0).mean()),
            "n_px_native": int(len(sub)),
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


def mixed_test(df: pd.DataFrame, response: str, fixed: str) -> dict:
    """Likelihood-ratio test of one fixed effect inside the authors' frame."""
    from scipy.stats import chi2
    f0 = _mixedlm(df, response, None)
    f1 = _mixedlm(df, response, fixed)
    lr = float(2.0 * (f1.llf - f0.llf))
    p = float(chi2.sf(max(lr, 0.0), df=1))
    r2m0, r2c0 = _nakagawa(f0, df)
    r2m1, r2c1 = _nakagawa(f1, df)
    return {
        "n": int(len(df)), "beta": float(f1.fe_params.get(fixed, np.nan)),
        "se": float(f1.bse.get(fixed, np.nan)),
        "lr_chi2": lr, "p_lrt": p,
        "r2_marginal": r2m1, "r2_conditional": r2c1,
        "delta_r2_marginal": r2m1 - r2m0,
        "var_explained_beyond_re": max(r2m1 - r2m0, 0.0),
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
              n_perm: int = 2000, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    res: dict = {"variance_decomposition": {}, "conditions": {}}
    res["variance_decomposition"]["all"] = variance_decomposition(roi, "drug_measure_norm")
    for c in ("Niraparib", "Rucaparib"):
        res["variance_decomposition"][c] = variance_decomposition(
            roi[roi.Compound == c], "drug_measure_norm")

    for cond in MODEL_CONDITIONS:
        d = attach_predictions(roi, sim, cond).dropna(subset=["pred_norm"])
        entry = {"n_roi": int(len(d)),
                 "n_roi_dropped_no_msi": int(len(roi) - len(d))}
        entry["pooled"] = mixed_test(d, "drug_measure_norm", "pred_norm")
        for c in ("Niraparib", "Rucaparib"):
            sub = d[d.Compound == c]
            entry[c] = mixed_test(sub, "drug_measure_norm", "pred_norm")

        # Slice-level view: n = 11 mapped slices, the real unit of information.
        sl = (d.groupby(["msi_slice", "Compound", "Patient"])
                .agg(meas=("drug_measure_norm", "mean"),
                     pred=("pred_norm", "first"), n=("drug_measure_norm", "size"))
                .reset_index())
        rho, prho = spearmanr(sl["pred"], sl["meas"])
        r, pr = pearsonr(sl["pred"], sl["meas"])
        entry["slice_level"] = {
            "n_slices": int(len(sl)), "spearman_rho": float(rho),
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
        res["conditions"][cond] = entry
    return res
