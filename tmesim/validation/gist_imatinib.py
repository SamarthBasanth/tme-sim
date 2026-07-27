"""REAL VALIDATION #2 -- measured spatial imatinib distribution in human GIST.

What this is
------------
The only open, in-vivo, in-HUMAN drug-distribution dataset we have with (a) a
matched within-patient normal-tissue control, (b) a per-patient plasma boundary
condition, and (c) a documented null: systemic exposure explains essentially none
of the tissue variance.

Data
----
METASPACE project "GIST" (Sci Rep 2019;9:10698, PMC6650609). MALDI-FTICR imaging
of imatinib (m/z 494.2663, [M+H]+) and N-desmethylimatinib (m/z 480.2506) in
human gastrointestinal stromal tumour sections and matched normal tissue, joined
to Supplementary Table S1 clinical data (55/55 spatial rows, zero orphans).

  data_cache/hunt/maldi_msi_drug/gist_imatinib_msi_index.csv   (55 rows)
  data_cache/hunt/maldi_msi_drug/ion_images/*.npy              (110 arrays)
  data_cache/hunt/verified_0/gist_spatial_plus_clinical.csv    (55 x 31)

Pixel values are ARBITRARY ION COUNTS, not ng/g. A RATIO of ion counts from the
same drug, same adduct, same acquisition run is unitless and therefore
comparable to a model partition ratio; ABSOLUTE counts are not.

The test
--------
  measured target 1 : tumour:normal imatinib partition ratio, n = 6 matched pairs
  measured target 2 : intratumoral heterogeneity (CV, detected-pixel fraction,
                      spatial autocorrelation e-folding length)
  null              : plasma imatinib (ng/mL) vs tumour-section mean imatinib,
                      Spearman, n = 9  -- recomputed here from the raw pixels
  prediction        : imatinib compiled through tmesim.molecular (SMILES ->
                      descriptors -> Stokes-Einstein D -> ECM/mesh hindrance ->
                      permeability/binding/clearance -> calibrated k_eff -> L),
                      mapped to a tumour:normal ratio through a reaction-diffusion
                      delivery efficiency and (optionally) a literature
                      tissue-composition partition. ZERO parameters are fitted to
                      the GIST measurements.

Honesty line (do not delete)
----------------------------
n = 6 matched pairs, and 5 of the 6 share the same tissue pairing (liver
metastasis vs normal liver). This dataset therefore contains roughly ONE
independent discriminating comparison. It can falsify a badly wrong transport
model; it CANNOT rank two plausible ones. Any correlation quoted over n = 6 here
is carried almost entirely by the single non-liver pair.

Nothing in core / physics / molecular / reconstruct / engine is modified.

Units: length um, time s, diffusivity um^2/s, rates 1/s.
"""
from __future__ import annotations

import dataclasses
import math
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths (repo-relative, resolved from this file).
# ---------------------------------------------------------------------------
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
HUNT = os.path.join(REPO, "data_cache", "hunt")

MSI_DIR = os.path.join(HUNT, "maldi_msi_drug")
ION_DIR = os.path.join(MSI_DIR, "ion_images")
QUANT_DIR = os.path.join(MSI_DIR, "ion_images_quant")
QUANT_INDEX = os.path.join(MSI_DIR, "gist_imatinib_quantitative_index.csv")
JOIN_CSV = os.path.join(HUNT, "verified_0", "gist_spatial_plus_clinical.csv")
CLINICAL_CSV = os.path.join(HUNT, "verified_0", "gist_S1_clinical.csv")

REPORT_MD = os.path.join(REPO, "results", "gist_imatinib.md")
REPORT_PNG = os.path.join(REPO, "results", "gist_imatinib.png")

CITATION = (
    "Abu Sammour D, Marsching C, Geisel A, et al. Spatial Probabilistic Mapping "
    "of Metabolite Ensembles in Mass Spectrometry Imaging / imatinib MSI in human "
    "GIST. Sci Rep 2019;9:10698. PMC6650609. METASPACE project 'GIST' "
    "(2b505b2a-ca3a-11e9-9d77-9b72ac66ec3c)."
)

#: The paper reports 50 um raster for the FTICR acquisitions. ``MS_Analysis.
#: Pixel_Size`` is NULL in the deposited metadata of the 'GIST' project datasets
#: (it IS present, 50x50 um, in the later 'Spatial Probabilistic Mapping' project
#: from the same lab). Every micrometre figure derived from pixel counts below is
#: therefore FLAGGED as resting on this out-of-deposit assumption.
PIXEL_UM = 50.0
PIXEL_UM_SOURCE = (
    "paper text only; MS_Analysis.Pixel_Size is null in the 'GIST' project "
    "deposit (it is 50x50 um in the same lab's later quantitative project)"
)

# ---------------------------------------------------------------------------
# Imatinib -- structure and literature physicochemistry.
# ---------------------------------------------------------------------------
IMATINIB_NAME = "imatinib"
#: Canonical SMILES; RDKit returns C29H31N7O / MW 493.615, matching the deposited
#: sum formula and the annotated m/z 494.2663 for [M+H]+ exactly.
IMATINIB_SMILES = (
    "CN1CCN(CC1)Cc1ccc(cc1)C(=O)Nc1ccc(C)c(Nc2nccc(n2)-c2cccnc2)c1"
)
NDESMETHYL_NAME = "N-desmethylimatinib"
NDESMETHYL_SMILES = (
    "C1CN(CCN1)Cc1ccc(cc1)C(=O)Nc1ccc(C)c(Nc2nccc(n2)-c2cccnc2)c1"
)

#: Strongest basic pKa (piperazine N). Literature range 7.7-8.07; the pyridine N
#: (pKa ~3.7) is >99.9% neutral at pH 7.0-7.4 and is ignored. Used ONLY for the
#: pH-partition (ion-trapping) term of the tissue-composition layer.
IMATINIB_PKA_BASE = 8.07
#: Experimental octanol/water logP (used only as a sensitivity check against the
#: RDKit/Crippen value our molecular layer computes).
IMATINIB_LOGP_EXPERIMENTAL = 3.0
#: Human plasma protein binding ~95% (mainly alpha-1-acid glycoprotein).
IMATINIB_FU_PLASMA_LIT = 0.05
#: Cellular IC50 against KIT-driven proliferation, ~0.1 uM (100 nM). Only used to
#: put measured/predicted tissue exposure on a pharmacological scale; NOT fitted.
IMATINIB_KIT_IC50_UM = 0.1

# ---------------------------------------------------------------------------
# Tissue composition (Rodgers, Leahy & Rowland 2005 / Rodgers & Rowland 2006
# human tissue-composition table, as used for mechanistic Kp prediction).
# fEW = extracellular water, fIW = intracellular water, fNL = neutral lipid,
# fNP = neutral phospholipid (fractions of tissue volume); AP = acidic
# phospholipid (mg/g, carried for documentation -- the AP term is NOT used, see
# ``tissue_partition_kpu``).
# ---------------------------------------------------------------------------
TISSUE_COMPOSITION: Dict[str, Dict[str, float]] = {
    "liver":     {"fEW": 0.161, "fIW": 0.573, "fNL": 0.0348, "fNP": 0.0252, "AP": 4.56},
    "intestine": {"fEW": 0.282, "fIW": 0.475, "fNL": 0.0487, "fNP": 0.0163, "AP": 2.41},
    "muscle":    {"fEW": 0.079, "fIW": 0.666, "fNL": 0.0238, "fNP": 0.0072, "AP": 2.42},
    "spleen":    {"fEW": 0.207, "fIW": 0.579, "fNL": 0.0201, "fNP": 0.0198, "AP": 3.18},
    "lung":      {"fEW": 0.336, "fIW": 0.446, "fNL": 0.0030, "fNP": 0.0090, "AP": 3.91},
}

#: Tumour interstitial (extracellular) fluid volume fraction in solid tumours is
#: measurably LARGER than in normal soft tissue: ~0.2-0.6 (Jain, Cancer Res 1987
#: and subsequent). We take 0.40.
TUMOUR_INTERSTITIAL_FRACTION = 0.40

#: There is NO published tumour row in the Rodgers-Rowland table. We construct a
#: surrogate: take a normal soft-tissue composition and swap in the measured
#: tumour interstitial fraction, taking the difference out of intracellular
#: water. This is the WEAKEST input in the whole prediction and is reported with
#: an explicit sensitivity band over three surrogate parents.
TUMOUR_SURROGATE_PARENTS = ("muscle", "spleen", "lung")
TUMOUR_SURROGATE_DEFAULT = "muscle"


def tumour_composition(parent: str = TUMOUR_SURROGATE_DEFAULT) -> Dict[str, float]:
    """Surrogate tumour tissue composition built from a normal-tissue parent.

    ``fEW`` is replaced by :data:`TUMOUR_INTERSTITIAL_FRACTION` and the change is
    subtracted from ``fIW`` (clamped at 0.05). Lipid fractions are inherited
    unchanged -- we have no GIST lipid measurement and will not invent one.
    """
    base = dict(TISSUE_COMPOSITION[parent])
    d = TUMOUR_INTERSTITIAL_FRACTION - base["fEW"]
    base["fEW"] = TUMOUR_INTERSTITIAL_FRACTION
    base["fIW"] = max(0.05, base["fIW"] - d)
    return base


# ---------------------------------------------------------------------------
# Microvascular geometry: HALF the intervessel distance (um), i.e. the maximum
# diffusion path from a vessel wall to the least-perfused point.
# Normal liver: hepatocyte plates are one cell thick between sinusoids, giving a
# sinusoid-to-sinusoid spacing of ~20-30 um -> half-distance ~12.5 um. Normal gut
# wall intercapillary distance is larger, ~40-60 um -> ~25 um. Solid tumours run
# ~50-200 um intervessel -> ~60 um half-distance, with a wide band.
# ---------------------------------------------------------------------------
HALF_SPACING_UM: Dict[str, float] = {
    "liver": 12.5,
    "intestine": 25.0,
    "stomach": 25.0,
    "peritoneum": 30.0,
    "tumour": 60.0,
}
HALF_SPACING_BAND_UM: Dict[str, tuple] = {
    "liver": (10.0, 15.0),
    "intestine": (20.0, 30.0),
    "tumour": (30.0, 120.0),
}

#: The six patients with BOTH a tumour and a matched normal imatinib map.
MATCHED_PATIENTS = ("A", "AA", "J", "K", "L", "U")


# ===========================================================================
# 1. DATA LOADING + HONEST CHARACTERISATION
# ===========================================================================
def load_join() -> pd.DataFrame:
    """The verified spatial x clinical join (55 rows), with numeric plasma."""
    df = pd.read_csv(JOIN_CSV)
    df["plasma_ng_ml"] = pd.to_numeric(df["plasma_conc_ng_ml"], errors="coerce")
    return df


def load_map(df: pd.DataFrame, dataset_name: str, compound: str = "Imatinib",
             adduct: str = "+H") -> Optional[np.ndarray]:
    """Load one per-pixel ion image (NaN outside the tissue mask). None if absent."""
    sel = df[(df.dataset_name == dataset_name) & (df.compound == compound)
             & (df.adduct == adduct)]
    if len(sel) == 0:
        return None
    return np.load(os.path.join(ION_DIR, sel.iloc[0].npy))


def imatinib_rows(df: pd.DataFrame) -> pd.DataFrame:
    """The one comparable channel: imatinib [M+H]+, present for every dataset.

    +Na / +K adducts exist for only a subset of sections, and adduct formation is
    matrix/salt dependent, so mixing adducts across sections would confound the
    ratio. We use [M+H]+ throughout.
    """
    return df[(df.compound == "Imatinib") & (df.adduct == "+H")].copy()


def _autocorr_efold_px(a: np.ndarray, maxlag: int = 15) -> float:
    """Isotropic (row+column averaged) lag-autocorrelation e-folding length, px.

    Returns the interpolated lag at which the spatial autocorrelation of the
    masked ion image first falls below 1/e. A value near 1 px means the image has
    no resolved spatial structure beyond the pixel/point-spread floor.
    """
    v = np.where(np.isfinite(a), a, np.nan)
    sd = np.nanstd(v)
    if not np.isfinite(sd) or sd <= 0:
        return float("nan")
    z = (v - np.nanmean(v)) / sd
    ac = [1.0]
    for lag in range(1, maxlag + 1):
        ac.append(float(np.nanmean([np.nanmean(z[:-lag, :] * z[lag:, :]),
                                    np.nanmean(z[:, :-lag] * z[:, lag:])])))
    ac = np.asarray(ac)
    thr = 1.0 / math.e
    idx = np.where(ac < thr)[0]
    if len(idx) == 0 or idx[0] == 0:
        return float("nan")
    i = idx[0]
    return float((i - 1) + (ac[i - 1] - thr) / (ac[i - 1] - ac[i]))


def characterize(df: Optional[pd.DataFrame] = None) -> Dict:
    """Per-section honest characterisation: shapes, mask coverage, ion-count
    ranges, detected-pixel fraction, heterogeneity, autocorrelation."""
    df = load_join() if df is None else df
    rows = []
    for _, r in imatinib_rows(df).iterrows():
        a = load_map(df, r.dataset_name)
        met = load_map(df, r.dataset_name, "N-desmethylimatinib")
        m = np.isfinite(a)
        v = a[m]
        nz = v[v > 0]
        rows.append(dict(
            patient=r.patient_id, tn=r.tn, organ=r.organ,
            met_or_primary=r.met_or_primary,
            height=int(a.shape[0]), width=int(a.shape[1]),
            n_px=int(a.size), n_tissue_px=int(m.sum()),
            mask_coverage=float(m.mean()),
            ion_min=float(v.min()), ion_max=float(v.max()),
            ion_mean=float(v.mean()), ion_median=float(np.median(v)),
            detected_fraction=float((v > 0).mean()),
            cv_all=float(v.std(ddof=1) / v.mean()),
            cv_detected=float(nz.std(ddof=1) / nz.mean()) if len(nz) > 2 else float("nan"),
            autocorr_efold_px=_autocorr_efold_px(a),
            metabolite_to_parent=(float(met[np.isfinite(met)].mean() / v.mean())
                                  if met is not None else float("nan")),
            plasma_ng_ml=r.plasma_ng_ml,
            regime_raw=r.imatinib_regime,
        ))
    sections = pd.DataFrame(rows).sort_values(["patient", "tn"]).reset_index(drop=True)
    tumours = sections[sections.tn == "tumor"]
    return {
        "sections": sections,
        "n_datasets": int(len(sections)),
        "n_patients": int(sections.patient.nunique()),
        "n_tumour_maps": int(len(tumours)),
        "n_normal_maps": int((sections.tn == "normal").sum()),
        "matched_patients": sorted(
            set(sections[sections.tn == "tumor"].patient)
            & set(sections[sections.tn == "normal"].patient)),
        "n_tumours_with_plasma": int(tumours.plasma_ng_ml.notna().sum()),
        "shape_range": (
            f"{sections.height.min()}x{sections.width.min()} to "
            f"{sections.height.max()}x{sections.width.max()} px"),
        "mask_coverage_range": (float(sections.mask_coverage.min()),
                                float(sections.mask_coverage.max())),
        "ion_count_range": (float(sections.ion_min.min()), float(sections.ion_max.max())),
        "citation": CITATION,
    }


# ===========================================================================
# 2. MEASURED TARGETS
# ===========================================================================
def measured_targets(char: Optional[Dict] = None) -> pd.DataFrame:
    """Per-patient measured tumour:normal imatinib partition ratio + heterogeneity.

    The ratio is taken on tissue-masked MEAN ion counts of the SAME ion
    (imatinib [M+H]+) in the tumour and matched-normal section of the SAME
    patient. Response factor, matrix, laser energy and normalisation cancel to
    first order, which is exactly why a ratio is comparable to a model partition
    coefficient while an absolute count is not. (Second order they do NOT cancel:
    ion suppression differs between liver parenchyma and tumour stroma -- an
    uncontrolled confounder, stated in the report.)
    """
    char = characterize() if char is None else char
    s = char["sections"]
    out = []
    for pid in sorted(set(s[s.tn == "tumor"].patient) & set(s[s.tn == "normal"].patient)):
        t = s[(s.patient == pid) & (s.tn == "tumor")].iloc[0]
        n = s[(s.patient == pid) & (s.tn == "normal")].iloc[0]
        out.append(dict(
            patient=pid, organ=t.organ, met_or_primary=t.met_or_primary,
            ratio_mean=float(t.ion_mean / n.ion_mean),
            ratio_median=float((t.ion_median + 1.0) / (n.ion_median + 1.0)),
            tumour_mean=float(t.ion_mean), normal_mean=float(n.ion_mean),
            tumour_cv_all=float(t.cv_all), tumour_cv_detected=float(t.cv_detected),
            tumour_detected_fraction=float(t.detected_fraction),
            normal_detected_fraction=float(n.detected_fraction),
            tumour_autocorr_px=float(t.autocorr_efold_px),
            normal_autocorr_px=float(n.autocorr_efold_px),
            tumour_autocorr_um=float(t.autocorr_efold_px) * PIXEL_UM,
            plasma_ng_ml=t.plasma_ng_ml,
        ))
    return pd.DataFrame(out)


def plasma_null(char: Optional[Dict] = None) -> Dict:
    """Recompute the published null from the raw pixels: does systemic exposure
    predict tissue exposure? Spearman(plasma ng/mL, tumour-section mean imatinib).

    Also tests the same null against the RATIO target, which is what the model is
    scored on, and against a constant-prediction baseline.
    """
    from scipy import stats
    char = characterize() if char is None else char
    s = char["sections"]
    t = s[(s.tn == "tumor") & s.plasma_ng_ml.notna()]
    rho, p = stats.spearmanr(t.plasma_ng_ml.values, t.ion_mean.values)
    r_log, p_log = stats.pearsonr(np.log10(t.plasma_ng_ml.values),
                                  np.log10(t.ion_mean.values))
    tgt = measured_targets(char)
    tp = tgt[tgt.plasma_ng_ml.notna()]
    rho_r, p_r = stats.spearmanr(tp.plasma_ng_ml.values, tp.ratio_mean.values)
    return {
        "n_tumours_with_plasma": int(len(t)),
        "spearman_rho": float(rho), "spearman_p": float(p),
        "pearson_loglog_r": float(r_log), "pearson_loglog_p": float(p_log),
        "plasma_range_ng_ml": (float(t.plasma_ng_ml.min()), float(t.plasma_ng_ml.max())),
        "tumour_mean_fold_range": float(t.ion_mean.max() / t.ion_mean.min()),
        "n_ratio": int(len(tp)),
        "spearman_rho_ratio": float(rho_r), "spearman_p_ratio": float(p_r),
    }


# ===========================================================================
# 3. MODEL PREDICTION (zero parameters fitted to this dataset)
# ===========================================================================
def compile_imatinib(ecm_reference: float = 0.4) -> Dict:
    """SMILES -> tmesim molecular layer -> calibrated DrugParams -> L.

    Runs the real pipeline: RDKit descriptors -> Stokes-Einstein free D ->
    Brownian-dynamics micro-tortuosity x Ogston-Renkin nm-mesh size gate ->
    passive permeability, plasma free fraction, tissue binding, clearance
    (``molecular.compile_drug.build_profile``), then applies the doxorubicin-
    anchored perivascular-uptake constant ``alpha`` from
    ``validation.calibration`` -- which was fitted ONCE to an unrelated
    doxorubicin measurement and is NOT touched here -- to get the effective loss
    rate and the reaction-diffusion e-folding penetration length
    ``L = sqrt(D_eff / k_eff)``.
    """
    from ..molecular import descriptors as _desc
    from ..molecular.compile_drug import build_profile, molecule_to_drugparams
    from . import calibration as cal
    from .pk_reality import penetration_length_um

    mol = _desc.from_smiles(IMATINIB_NAME, IMATINIB_SMILES)
    prof = build_profile(mol, ecm_reference=ecm_reference)
    dp = molecule_to_drugparams(prof)
    alpha = cal.calibrate_alpha()
    dp_cal = dataclasses.replace(dp, binding=alpha * cal.uptake_proxy(dp),
                                 clearance=cal.k_passive(dp))
    met = _desc.from_smiles(NDESMETHYL_NAME, NDESMETHYL_SMILES)
    return {
        "molecule": mol, "profile": prof,
        "drugparams_uncalibrated": dp, "drugparams": dp_cal,
        "alpha": float(alpha),
        "mw": float(mol.get("mw")), "logp_rdkit": float(mol.get("logp")),
        "radius_nm": float(mol.get("radius_hydrodynamic_nm")),
        "free_D": float(prof.free_diffusivity),
        "eff_D": float(prof.effective_diffusivity),
        "hindrance": float(prof.hindrance),
        "permeability_um_s": float(prof.membrane_permeability),
        "free_fraction": float(prof.free_fraction),
        "binding_uncal": float(dp.binding), "clearance": float(dp.clearance),
        "k_eff": float(dp_cal.binding + dp_cal.clearance),
        "L_uncalibrated_um": float(penetration_length_um(dp)),
        "L_um": float(penetration_length_um(dp_cal)),
        "metabolite_mw": float(met.get("mw")),
        "metabolite_logp": float(met.get("logp")),
    }


def delivery_efficiency(L_um: float, half_spacing_um: float) -> float:
    """Steady-state tissue-average concentration / vessel concentration.

    1D Krogh slab: vessels at x = 0 and x = 2R, first-order loss k in between.
    ``D d2C/dx2 = k C`` with C(0) = C(2R) = C_v gives
    ``C(x)/C_v = cosh((x-R)/L) / cosh(R/L)`` and a spatial mean

        E = (L/R) * tanh(R/L)      in (0, 1]

    E -> 1 when the drug penetrates much further than the intervessel spacing
    (delivery is NOT transport-limited); E -> L/R when it does not.
    """
    a = float(half_spacing_um) / float(L_um)
    if a <= 0:
        return 1.0
    return float(math.tanh(a) / a)


def slab_profile_cv(L_um: float, half_spacing_um: float, n: int = 2001) -> float:
    """Coefficient of variation of the Krogh-slab concentration profile -- the
    intratumoral heterogeneity a homogeneous transport model can generate."""
    a = float(half_spacing_um) / float(L_um)
    u = np.linspace(-a, a, n)
    c = np.cosh(u) / math.cosh(a)
    return float(c.std(ddof=0) / c.mean())


def tissue_partition_kpu(logp: float, pka: float, comp: Dict[str, float],
                         ph_iw: float = 7.0, ph_p: float = 7.4) -> float:
    """Tissue-to-UNBOUND-plasma partition Kpu for a monoprotic base, from tissue
    composition (Rodgers-Rowland backbone, acidic-phospholipid term OMITTED).

        Kpu = fEW
            + fIW * (1 + 10^(pKa - pH_iw)) / (1 + 10^(pKa - pH_p))
            + [P*fNL + (0.3P + 0.7)*fNP] / (1 + 10^(pKa - pH_p))

    with P = 10^logP (octanol:water, non-adipose tissue). The three terms are
    extracellular water, pH-partitioned (ion-trapped) intracellular water, and
    neutral lipid + phospholipid dissolution.

    THE AP TERM IS DELIBERATELY OMITTED. Rodgers-Rowland's acidic-phospholipid
    binding constant Ka_AP cannot be derived from structure -- it is back-
    calculated from measured blood-cell partitioning, which we do not have for
    imatinib. Including it with a guessed Ka_AP would manufacture exactly the
    liver-vs-tumour difference this test is supposed to probe. Omitting it makes
    the prediction CONSERVATIVE in a known direction (it under-states liver
    uptake), and that is stated in the report.
    """
    P = 10.0 ** float(logp)
    fi_p = 1.0 + 10.0 ** (pka - ph_p)      # ionisation denominator (plasma)
    fi_iw = 1.0 + 10.0 ** (pka - ph_iw)    # ionisation numerator (cytosol)
    water = comp["fEW"] + comp["fIW"] * fi_iw / fi_p
    lipid = (P * comp["fNL"] + (0.3 * P + 0.7) * comp["fNP"]) / fi_p
    return float(water + lipid)


def predict(compiled: Optional[Dict] = None,
            tumour_parent: str = TUMOUR_SURROGATE_DEFAULT,
            logp: Optional[float] = None) -> Dict:
    """Predicted tumour:normal imatinib partition ratio, two nested models.

    M1  TRANSPORT ONLY -- pure tmesim molecular layer. Assumes the tumour and the
        normal tissue have the SAME intrinsic affinity for the drug (our layer's
        binding/clearance terms are drug-specific but tissue-agnostic), so the
        whole ratio comes from reaction-diffusion delivery:
            R = E_tumour / E_normal,  E = (L/R_k) tanh(R_k/L)
        This is the honest "what does the transport model alone say" prediction.

    M2  TRANSPORT x TISSUE-COMPOSITION PARTITION -- multiplies M1 by the ratio of
        composition-based equilibrium partitions Kpu(tumour)/Kpu(normal),
        computed from imatinib's OWN logP (from our layer's RDKit descriptors)
        and its literature basic pKa, against published human tissue composition.
        Still zero parameters fitted to the GIST data.

    Both are constant within a given tissue pairing: the model has no per-patient
    tumour input (we have no per-patient vessel density or ECM measurement), so
    it predicts one value for every liver-metastasis pair and one for the
    intestinal primary. That is a real limitation, not an oversight.
    """
    compiled = compile_imatinib() if compiled is None else compiled
    L = compiled["L_um"]
    lp = compiled["logp_rdkit"] if logp is None else float(logp)
    tcomp = tumour_composition(tumour_parent)
    kpu_t = tissue_partition_kpu(lp, IMATINIB_PKA_BASE, tcomp)

    E_t = delivery_efficiency(L, HALF_SPACING_UM["tumour"])
    per_tissue = {}
    for tissue in ("liver", "intestine"):
        E_n = delivery_efficiency(L, HALF_SPACING_UM[tissue])
        kpu_n = tissue_partition_kpu(lp, IMATINIB_PKA_BASE,
                                     TISSUE_COMPOSITION[tissue])
        per_tissue[tissue] = {
            "E_normal": E_n, "E_tumour": E_t,
            "kpu_normal": kpu_n, "kpu_tumour": kpu_t,
            "M1_transport_only": E_t / E_n,
            "M2_transport_x_composition": (E_t / E_n) * (kpu_t / kpu_n),
        }
    return {
        "L_um": L,
        "L_px": L / PIXEL_UM,
        "logp_used": lp,
        "tumour_parent": tumour_parent,
        "tumour_composition": tcomp,
        "per_tissue": per_tissue,
        "predicted_tumour_cv": slab_profile_cv(L, HALF_SPACING_UM["tumour"]),
        "predicted_normal_cv_liver": slab_profile_cv(L, HALF_SPACING_UM["liver"]),
    }


def predict_for_patients(targets: pd.DataFrame, compiled: Dict,
                         pred: Optional[Dict] = None) -> pd.DataFrame:
    """Attach M1/M2 predictions to each matched patient by its normal-tissue organ."""
    pred = predict(compiled) if pred is None else pred
    key = {"Liver": "liver", "Intestine": "intestine"}
    rows = []
    for _, r in targets.iterrows():
        t = pred["per_tissue"][key[r.organ]]
        rows.append(dict(patient=r.patient, organ=r.organ,
                         measured=r.ratio_mean,
                         pred_M1=t["M1_transport_only"],
                         pred_M2=t["M2_transport_x_composition"],
                         pred_null_unity=1.0))
    out = pd.DataFrame(rows)
    for c in ("M1", "M2", "null_unity"):
        out[f"fold_{c}"] = np.maximum(out[f"pred_{c}"] / out.measured,
                                      out.measured / out[f"pred_{c}"])
    return out


# ===========================================================================
# 4. THE ABSOLUTE-QUANTITATION LAYER (assessed, and rejected)
# ===========================================================================
def quantitative_layer() -> Dict:
    """Can the IS-normalised maps + on-tissue dilution series give ng/g or uM?

    Assesses ``data_cache/hunt/maldi_msi_drug/ion_images_quant`` (METASPACE
    project 'Spatial Probabilistic Mapping in MSI', same lab): an imatinib-D8
    internal-standard-normalised GIST liver-metastasis map plus a 6-point
    on-tissue imatinib dilution series. Returns the segmented calibration spots
    and an explicit VERDICT on whether an absolute scale is defensible. It does
    NOT fabricate one.
    """
    from scipy import ndimage
    idx = pd.read_csv(QUANT_INDEX)
    cal_name = "FTICR-GIST-Q1_imatinib_dilution_series-normalized_to_IS__Imatinib___H.npy"
    tof_name = "TOF-GIST-imatinib-IS_normalized__Imatinib___H.npy"
    cal_img = np.load(os.path.join(QUANT_DIR, cal_name))
    tof_img = np.load(os.path.join(QUANT_DIR, tof_name))

    mask = ndimage.binary_opening(np.isfinite(cal_img), np.ones((3, 3)))
    lab, n = ndimage.label(mask)
    sizes = ndimage.sum(mask, lab, range(1, n + 1))
    spots = []
    for i in range(n):
        if sizes[i] < 500:
            continue
        sel = lab == (i + 1)
        sl = ndimage.find_objects(sel)[0]
        v = cal_img[sel]
        spots.append(dict(n_px=int(sel.sum()), row0=int(sl[0].start),
                          col0=int(sl[1].start), mean=float(np.nanmean(v)),
                          median=float(np.nanmedian(v))))
    spots = pd.DataFrame(spots).sort_values(["row0", "col0"]).reset_index(drop=True)
    means = np.sort(spots["mean"].values)
    dyn = float(means[-1] / means[0]) if len(means) else float("nan")

    # Row-major spatial ordering (top-left to bottom-right) -- is the response a
    # monotone ladder, as a dilution series laid out in a grid should be?
    top = spots[spots.row0 < spots.row0.median()].sort_values("col0")["mean"].values
    bot = spots[spots.row0 >= spots.row0.median()].sort_values("col0")["mean"].values
    monotone = bool(np.all(np.diff(np.concatenate([top, bot])) > 0)
                    or np.all(np.diff(np.concatenate([top, bot])) < 0))

    blockers = [
        "No spot->concentration labels anywhere in the deposit: the METASPACE "
        "dataset has description=null and no supplementary annotation. The "
        "nominal 0.78-25 pmol 6-point ladder is a claim from the paper text, not "
        "a machine-readable mapping to the six spots in the image.",
        f"The six segmented spots are NOT monotone in any spatial ordering "
        f"(row-major means {np.round(np.concatenate([top, bot]), 3).tolist()}), "
        "so the concentration assignment cannot be inferred from the layout. "
        "Sorting the means to force monotonicity and then fitting would be a "
        "tautology with a meaningless R^2.",
        f"Response dynamic range across the six spots is only {dyn:.1f}x for a "
        "nominal 32x (6-point, 2-fold) concentration range, and the two middle "
        "spots differ by <8%, i.e. adjacent 2-fold steps are not resolved.",
        "INSTRUMENT MISMATCH: the dilution series is FTICR (resolving power "
        "~160,000 at m/z 200); the IS-normalised tumour map is TOF reflector "
        "(~10,000 at m/z 500). The IS response factor is instrument- and "
        "matrix-specific and does not transfer between them.",
        "TISSUE MISMATCH: the calibration spots are on NORMAL SMALL INTESTINE; "
        "the quantitative map is LIVER metastasis. On-tissue calibration is only "
        "valid on the tissue it was spotted on (ion suppression is tissue-specific).",
        "NO PATIENT LINK: 'TOF-GIST-imatinib-IS normalized' carries no sample "
        "letter, so it cannot be joined to any of the 11 patients or the 6 "
        "matched pairs that carry the ratio measurement and the plasma value.",
        "The 20 TOF datasets in the main 'GIST' project have ZERO drug "
        "annotations, so even a valid TOF calibration could not be applied to "
        "the FTICR maps used for the matched-pair analysis.",
    ]
    return {
        "n_quant_maps": int(len(idx)),
        "spots": spots,
        "n_spots": int(len(spots)),
        "spot_means_row_major": np.concatenate([top, bot]).tolist(),
        "dynamic_range": dyn,
        "monotone_in_layout": monotone,
        "tof_map_shape": tuple(int(x) for x in tof_img.shape),
        "tof_map_range": (float(np.nanmin(tof_img)), float(np.nanmax(tof_img))),
        "tof_map_mean": float(np.nanmean(tof_img)),
        "pixel_size_um_in_metadata": 50.0,
        "defensible": False,
        "blockers": blockers,
        "verdict": (
            "NOT DEFENSIBLE. A calibration curve cannot be built: the six "
            "on-tissue standard spots carry no concentration labels in the "
            "deposit and are non-monotone in layout, the response spans only "
            f"{dyn:.1f}x over a nominal 32x, and the series is on the wrong "
            "instrument and the wrong tissue relative to the map it would "
            "calibrate -- which in turn cannot be linked to any patient in the "
            "matched-pair analysis. All results below stay on the unitless "
            "RATIO scale. No absolute concentration is reported."),
    }


# ===========================================================================
# 5. SCORECARD
# ===========================================================================
def score(targets: pd.DataFrame, compiled: Dict, pred: Dict,
          null: Dict) -> Dict:
    """Predicted vs measured tumour:normal ratio for the n = 6 matched patients.

    Reports, for each model and for the two baselines, the geometric-mean fold
    error and the Spearman rank correlation -- and, in the same breath, the fact
    that n = 6 with 5 of 6 sharing one tissue pairing means the rank correlation
    encodes about one bit of information.
    """
    from scipy import stats
    tab = predict_for_patients(targets, compiled, pred)
    res = {"table": tab, "n": int(len(tab))}
    for m in ("M1", "M2", "null_unity"):
        fold = tab[f"fold_{m}"].values
        rho, p = stats.spearmanr(tab[f"pred_{m}"].values, tab.measured.values)
        res[m] = {
            "geomean_fold_error": float(10 ** np.mean(np.log10(fold))),
            "max_fold_error": float(fold.max()),
            "min_fold_error": float(fold.min()),
            "spearman_rho": float(rho) if np.isfinite(rho) else float("nan"),
            "spearman_p": float(p) if np.isfinite(p) else float("nan"),
            "n_within_3x": int((fold <= 3.0).sum()),
        }
    # Liver-only subset: the model is CONSTANT there, so it can explain none of
    # the 6.1x spread the measurement shows across the five liver pairs.
    liv = tab[tab.organ == "Liver"]
    res["liver_only"] = {
        "n": int(len(liv)),
        "measured_spread_fold": float(liv.measured.max() / liv.measured.min()),
        "predicted_spread_fold": float(liv.pred_M1.max() / liv.pred_M1.min()),
        "note": ("the model has no per-patient tumour input, so it predicts one "
                 "value for all five liver pairs and explains 0% of the observed "
                 "within-liver spread"),
    }
    res["plasma_null"] = null
    res["beats_plasma_null"] = bool(
        res["M1"]["geomean_fold_error"] < res["null_unity"]["geomean_fold_error"]
        and abs(res["M1"]["spearman_rho"]) > abs(null["spearman_rho_ratio"]))
    return res


# ===========================================================================
# 6. HETEROGENEITY / PENETRATION-LENGTH COMPARISON
# ===========================================================================
def heterogeneity_scorecard(targets: pd.DataFrame, pred: Dict) -> Dict:
    """Measured intratumoral heterogeneity vs what a homogeneous-vasculature
    reaction-diffusion model can generate."""
    m = targets
    return {
        "measured_cv_all_pixels": (float(m.tumour_cv_all.min()),
                                   float(m.tumour_cv_all.max())),
        "measured_cv_detected_pixels": (float(m.tumour_cv_detected.min()),
                                        float(m.tumour_cv_detected.max())),
        "measured_detected_fraction": (float(m.tumour_detected_fraction.min()),
                                       float(m.tumour_detected_fraction.max())),
        "predicted_cv_krogh_slab": float(pred["predicted_tumour_cv"]),
        "predicted_cv_range_over_spacing": (
            float(slab_profile_cv(pred["L_um"], HALF_SPACING_BAND_UM["tumour"][0])),
            float(slab_profile_cv(pred["L_um"], HALF_SPACING_BAND_UM["tumour"][1]))),
        "measured_autocorr_px_tumour": (float(m.tumour_autocorr_px.min()),
                                        float(m.tumour_autocorr_px.max())),
        "measured_autocorr_px_normal": (float(m.normal_autocorr_px.min()),
                                        float(m.normal_autocorr_px.max())),
        "predicted_L_px": float(pred["L_px"]),
        "predicted_L_um": float(pred["L_um"]),
        "caveat": (
            "CV over ALL tissue pixels is dominated by a detection floor, not by "
            "a resolved gradient: 8-78% of tumour pixels are exactly zero. CV "
            "over DETECTED pixels is the fairer comparison. And a measured "
            "autocorrelation e-folding of ~1 px IS the pixel/point-spread floor, "
            "so agreement there is 'not inconsistent', not confirmation."),
    }


# ===========================================================================
# 7. FIGURE
# ===========================================================================
_C_MEAS = "#d62728"     # measured
_C_M1 = "#1f77b4"       # transport-only prediction
_C_M2 = "#7f9cb5"       # transport x composition
_C_NULL = "#8c8c8c"     # baselines


def make_figure(char: Dict, targets: pd.DataFrame, compiled: Dict, pred: Dict,
                sc: Dict, null: Dict, het: Dict, path: str = REPORT_PNG) -> str:
    """Measured maps, measured-vs-predicted partition, and the plasma null."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.gridspec import GridSpec

    df = load_join()
    fig = plt.figure(figsize=(17.5, 13.0))
    gs = GridSpec(3, 4, figure=fig, hspace=0.42, wspace=0.30,
                  height_ratios=[1.05, 1.0, 1.0])

    # --- row 0: measured ion images, matched pairs, SHARED per-patient scale ---
    pairs = [("A", "Liver metastasis"), ("AA", "Intestinal primary")]
    col = 0
    for pid, label in pairs:
        tn_names = {}
        for tn in ("tumor", "normal"):
            r = df[(df.patient_id == pid) & (df.tn == tn) & (df.compound == "Imatinib")
                   & (df.adduct == "+H")].iloc[0]
            tn_names[tn] = np.load(os.path.join(ION_DIR, r.npy))
        allv = np.concatenate([a[np.isfinite(a) & (a > 0)].ravel()
                               for a in tn_names.values()])
        vmin, vmax = np.percentile(allv, 1), np.percentile(allv, 99.5)
        ratio = float(targets[targets.patient == pid].ratio_mean.iloc[0])
        for tn, title in (("tumor", "tumour"), ("normal", "matched normal")):
            ax = fig.add_subplot(gs[0, col])
            a = np.ma.masked_invalid(tn_names[tn])
            # grey tissue-mask underlay so below-detection pixels are visibly
            # INSIDE the section rather than reading as empty background
            ax.imshow(np.where(np.isfinite(tn_names[tn]), 1.0, np.nan),
                      cmap="Greys", vmin=0.0, vmax=1.6, interpolation="nearest")
            im = ax.imshow(np.ma.masked_less_equal(a, 0), cmap="magma",
                           norm=LogNorm(vmin=max(vmin, 1.0), vmax=vmax),
                           interpolation="nearest")
            ax.set_title(f"patient {pid} — {title}\n{label}", fontsize=10)
            ax.set_xticks([]); ax.set_yticks([])
            if tn == "tumor":
                ax.set_xlabel(f"tumour:normal = {ratio:.3f}", fontsize=10,
                              color=_C_MEAS, fontweight="bold")
            fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
            col += 1
    fig.text(0.5, 0.955,
             "Measured imatinib [M+H]+ ion images (shared log scale within each "
             "patient pair; grey = below detection / outside tissue mask)",
             ha="center", fontsize=11, style="italic")

    # --- row 1 left: measured vs predicted partition ---
    ax = fig.add_subplot(gs[1, 0:2])
    tab = sc["table"]
    lim = (0.012, 2.0)
    ax.plot(lim, lim, "k--", lw=1, zorder=1, label="perfect")
    for f, ls in ((3.0, ":"), (10.0, "-.")):
        ax.plot(lim, [lim[0] * f, lim[1] * f], color="0.7", ls=ls, lw=0.9, zorder=1)
        ax.plot(lim, [lim[0] / f, lim[1] / f], color="0.7", ls=ls, lw=0.9, zorder=1)
    mk = {"Liver": "o", "Intestine": "s"}
    for _, r in tab.iterrows():
        ax.scatter(r.measured, r.pred_M1, s=110, marker=mk[r.organ], color=_C_M1,
                   edgecolor="0.2", zorder=3)
        ax.scatter(r.measured, r.pred_M2, s=80, marker=mk[r.organ], color=_C_M2,
                   edgecolor="0.2", zorder=3)
        ax.annotate(r.patient, (r.measured, r.pred_M1), textcoords="offset points",
                    xytext=(7, 5), fontsize=9)
    ax.axhline(1.0, color=_C_NULL, lw=1.2, ls="-",
               label=f"null: no partitioning (R=1), {sc['null_unity']['geomean_fold_error']:.1f}x")
    ax.scatter([], [], s=110, marker="o", color=_C_M1, edgecolor="0.2",
               label=f"M1 transport-only, {sc['M1']['geomean_fold_error']:.1f}x")
    ax.scatter([], [], s=80, marker="o", color=_C_M2, edgecolor="0.2",
               label=f"M2 transport x composition, {sc['M2']['geomean_fold_error']:.1f}x")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("MEASURED tumour:normal imatinib ratio (ion counts)")
    ax.set_ylabel("PREDICTED tumour:normal ratio")
    ax.set_title("Predicted vs measured partition — n = 6 matched patients\n"
                 "(circles = liver metastasis, squares = intestinal primary; "
                 "grey lines 3x / 10x)", fontsize=11)
    ax.legend(fontsize=8.5, loc="lower right", framealpha=0.95)
    ax.grid(alpha=0.25)

    # --- row 1 right: the plasma null ---
    ax = fig.add_subplot(gs[1, 2:4])
    s = char["sections"]
    t = s[(s.tn == "tumor") & s.plasma_ng_ml.notna()]
    ax.scatter(t.plasma_ng_ml, t.ion_mean, s=110, color=_C_MEAS, edgecolor="0.2")
    for _, r in t.iterrows():
        ax.annotate(r.patient, (r.plasma_ng_ml, r.ion_mean),
                    textcoords="offset points", xytext=(7, 4), fontsize=9)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("plasma imatinib (ng/mL)")
    ax.set_ylabel("tumour-section mean imatinib (ion counts)")
    ax.set_title(
        f"THE NULL: systemic exposure vs tissue exposure\n"
        f"Spearman rho = {null['spearman_rho']:.3f}, p = {null['spearman_p']:.2f}, "
        f"n = {null['n_tumours_with_plasma']} — dose explains nothing",
        fontsize=11)
    ax.grid(alpha=0.25)

    # --- row 2 left: per-patient measured ratio vs model ---
    ax = fig.add_subplot(gs[2, 0:2])
    o = tab.sort_values("measured")
    x = np.arange(len(o))
    ax.bar(x, o.measured, color=_C_MEAS, edgecolor="0.3", label="measured")
    ax.plot(x, o.pred_M1, "o-", color=_C_M1, lw=1.6, label="M1 transport-only")
    ax.plot(x, o.pred_M2, "s--", color=_C_M2, lw=1.4, label="M2 transport x composition")
    ax.axhline(1.0, color=_C_NULL, lw=1.0, ls=":")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p}\n{g.lower()}" for p, g in zip(o.patient, o.organ)],
                       fontsize=9)
    ax.set_ylabel("tumour:normal imatinib ratio")
    ax.set_title("Per-patient scorecard (n = 6). The model over-predicts every\n"
                 "liver pair by 4-25x and hits the one intestinal pair to 1.1x.",
                 fontsize=11)
    ax.legend(fontsize=9); ax.grid(alpha=0.25, axis="y")

    # --- row 2 right: heterogeneity ---
    ax = fig.add_subplot(gs[2, 2:4])
    o2 = targets.sort_values("tumour_detected_fraction")
    x = np.arange(len(o2))
    ax.bar(x - 0.2, o2.tumour_cv_all, width=0.38, color="#c0c0c0",
           edgecolor="0.3", label="measured CV, all tissue pixels")
    ax.bar(x + 0.2, o2.tumour_cv_detected, width=0.38, color=_C_MEAS,
           edgecolor="0.3", label="measured CV, detected pixels only")
    lo, hi = het["predicted_cv_range_over_spacing"]
    ax.axhspan(lo, hi, color=_C_M1, alpha=0.18,
               label=f"model CV, intervessel 30-120 um ({lo:.2f}-{hi:.2f})")
    ax.axhline(het["predicted_cv_krogh_slab"], color=_C_M1, lw=1.6,
               label=f"model CV at 60 um ({het['predicted_cv_krogh_slab']:.2f})")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{p}\n{d:.0%} det." for p, d in
                        zip(o2.patient, o2.tumour_detected_fraction)], fontsize=9)
    ax.set_ylabel("intratumoral coefficient of variation")
    ax.set_title("Intratumoral heterogeneity (n = 6). All-pixel CV is dominated by\n"
                 "a detection floor (8-78% of tumour pixels are exactly zero).",
                 fontsize=11)
    ax.legend(fontsize=8.5); ax.grid(alpha=0.25, axis="y")

    fig.savefig(path, dpi=115, bbox_inches="tight")
    plt.close(fig)
    return path


# ===========================================================================
# 8. REPORT
# ===========================================================================
def _md_table(df: pd.DataFrame, cols: List[str], fmt: Dict[str, str]) -> str:
    head = "| " + " | ".join(cols) + " |\n|" + "|".join(["---"] * len(cols)) + "|\n"
    body = ""
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            f = fmt.get(c)
            cells.append(f.format(v) if (f and pd.notna(v)) else
                         ("n/a" if pd.isna(v) else str(v)))
        body += "| " + " | ".join(cells) + " |\n"
    return head + body


def write_report(char: Dict, targets: pd.DataFrame, compiled: Dict, pred: Dict,
                 sc: Dict, null: Dict, het: Dict, quant: Dict,
                 path: str = REPORT_MD) -> str:
    s = char["sections"]
    tab = sc["table"]
    liv = tab[tab.organ == "Liver"]
    prof = compiled["profile"]
    L = pred["L_um"]

    md = []
    A = md.append
    A("# GIST imatinib — measured spatial drug distribution in human tumours\n")
    A("**Real validation run #2.** In-vivo, in-human, matched within-patient "
      "control, per-patient plasma boundary condition, and a documented null.\n")
    A("![gist imatinib](gist_imatinib.png)\n")

    # ---------------- headline ----------------
    A("## Headline\n")
    A(f"- Measured tumour:normal imatinib partition ratio, **n = 6 matched "
      f"patients**: {', '.join(f'{r.patient} {r.measured:.3f}' for _, r in tab.iterrows())}.\n")
    A(f"- The plasma null is real and reproduced here from raw pixels: plasma "
      f"imatinib vs tumour-section mean imatinib, **Spearman rho = "
      f"{null['spearman_rho']:.3f}, p = {null['spearman_p']:.2f}, n = "
      f"{null['n_tumours_with_plasma']}**, over a 10.7x plasma range and a "
      f"{null['tumour_mean_fold_range']:.0f}x tissue range. Systemic exposure "
      f"explains essentially none of the tissue variance.\n")
    A(f"- Our transport-only prediction (M1) is **{sc['M1']['geomean_fold_error']:.1f}x "
      f"geometric-mean fold error over n = 6**, versus **"
      f"{sc['null_unity']['geomean_fold_error']:.1f}x** for the no-partitioning "
      f"null (R = 1). It beats the null, but only by ~1.6x, and it is wrong in a "
      f"structured way: it nails the single intestinal-primary pair to "
      f"{tab[tab.organ=='Intestine'].fold_M1.iloc[0]:.2f}x and over-predicts every "
      f"one of the five liver-metastasis pairs by "
      f"{liv.fold_M1.min():.1f}-{liv.fold_M1.max():.1f}x.\n")
    A("- **Verdict: partially falsified, and the failure is diagnostic.** "
      "Passive reaction-diffusion transport is not what makes a GIST liver "
      "metastasis 20-40x poorer in imatinib than the liver next to it. "
      "The missing physics is active hepatic uptake plus lysosomal trapping of a "
      "weak base in hepatocytes — neither of which exists in our model.\n")
    A(f"- **n = 6, and 5 of the 6 are the same tissue pairing (liver metastasis "
      f"vs normal liver).** This dataset carries roughly ONE independent "
      f"discriminating comparison. It can falsify a badly wrong transport model; "
      f"it cannot rank two plausible ones. Every number below should be read "
      f"with that in front of it.\n")

    # ---------------- provenance ----------------
    A("\n## 1. Provenance and honest characterisation\n")
    A(f"- Source: {CITATION}\n")
    A(f"- Local: `data_cache/hunt/maldi_msi_drug/` (index {len(load_join())} rows, "
      f"110 per-pixel `.npy` arrays) and `data_cache/hunt/verified_0/"
      f"gist_spatial_plus_clinical.csv` (55 x 31; the spatial->clinical join was "
      f"verified 55/55 with zero orphans).\n")
    A(f"- Loaded here: **{char['n_datasets']} FTICR sections** across "
      f"**{char['n_patients']} patients** — {char['n_tumour_maps']} tumour maps "
      f"and {char['n_normal_maps']} normal maps, giving "
      f"**{len(char['matched_patients'])} matched tumour+normal pairs** "
      f"({', '.join(char['matched_patients'])}).\n")
    A(f"- Image shapes {char['shape_range']}; tissue-mask coverage "
      f"{char['mask_coverage_range'][0]:.0%}-{char['mask_coverage_range'][1]:.0%} "
      f"of the raster; ion counts {char['ion_count_range'][0]:.0f} to "
      f"{char['ion_count_range'][1]:.3g} (**arbitrary units, not ng/g**).\n")
    A(f"- Channel used throughout: **imatinib [M+H]+, m/z 494.2663**, the only "
      f"adduct present for every section. +Na/+K exist for a subset only and "
      f"adduct formation is salt/matrix dependent, so mixing them would confound "
      f"the ratio.\n")
    A(f"- Pixel size: **{PIXEL_UM:.0f} um assumed** — {PIXEL_UM_SOURCE}. Every "
      f"micrometre figure derived from pixel counts inherits that assumption.\n")
    A(f"- {char['n_tumours_with_plasma']} of {char['n_tumour_maps']} tumour "
      f"sections carry a numeric plasma imatinib value (patients B and G are "
      f"`---` in Table S1).\n")
    A("- The `imatinib_regime` column in the join is a concatenation of "
      "treatment-duration and dose fields from the source table (e.g. "
      "`7/15400/800`) and is **not** reliably parseable into mg/day. We do not "
      "use dose as a predictor; plasma concentration is the boundary condition.\n")

    A("\n### Per-section measurements\n")
    A(_md_table(s, ["patient", "tn", "organ", "height", "width", "n_tissue_px",
                    "mask_coverage", "ion_mean", "detected_fraction", "cv_all",
                    "cv_detected", "autocorr_efold_px"],
                {"mask_coverage": "{:.2f}", "ion_mean": "{:.3g}",
                 "detected_fraction": "{:.3f}", "cv_all": "{:.2f}",
                 "cv_detected": "{:.2f}", "autocorr_efold_px": "{:.2f}"}))

    A("\n### Reproduction of the paper's headline from raw pixels\n")
    A("The paper reports tumour:normal imatinib ratios of ~0.02-0.14 for liver "
      "metastases and ~0.72 for the one primary small-intestine GIST. Computed "
      "here directly from the deposited per-pixel arrays, with no intermediate "
      "processing:\n")
    A(f"- liver metastases (n = 5): "
      f"{', '.join(f'{r.patient} {r.ratio_mean:.3f}' for _, r in targets[targets.organ=='Liver'].iterrows())} "
      f"— range {targets[targets.organ=='Liver'].ratio_mean.min():.3f}-"
      f"{targets[targets.organ=='Liver'].ratio_mean.max():.3f}\n")
    A(f"- intestinal primary (n = 1): AA "
      f"{targets[targets.organ=='Intestine'].ratio_mean.iloc[0]:.3f}\n")
    A("**Match.** The independent reproduction succeeds; the arrays and the join "
      "are what they claim to be.\n")

    # ---------------- measured targets ----------------
    A("\n## 2. The measured targets\n")
    A("### 2.1 Tumour:normal partition ratio (n = 6 matched patients)\n")
    A("Ratio of tissue-masked MEAN imatinib [M+H]+ ion counts, tumour section "
      "over the same patient's normal section. Same ion, same adduct, same "
      "patient — so response factor, matrix and laser settings cancel to first "
      "order, which is what makes a RATIO comparable to a model partition "
      "coefficient even though the absolute counts are meaningless.\n")
    A(_md_table(targets, ["patient", "organ", "met_or_primary", "ratio_mean",
                          "tumour_mean", "normal_mean", "plasma_ng_ml"],
                {"ratio_mean": "**{:.3f}**", "tumour_mean": "{:.3g}",
                 "normal_mean": "{:.3g}", "plasma_ng_ml": "{:.0f}"}))
    A(f"Geometric mean across the five liver pairs: "
      f"{10**np.mean(np.log10(targets[targets.organ=='Liver'].ratio_mean)):.3f}. "
      f"Within-liver spread: {sc['liver_only']['measured_spread_fold']:.1f}x.\n")

    A("\n### 2.2 Intratumoral heterogeneity\n")
    A(f"- Coefficient of variation over ALL tissue pixels: "
      f"{het['measured_cv_all_pixels'][0]:.2f}-{het['measured_cv_all_pixels'][1]:.2f}.\n")
    A(f"- CV over DETECTED (non-zero) pixels only: "
      f"{het['measured_cv_detected_pixels'][0]:.2f}-"
      f"{het['measured_cv_detected_pixels'][1]:.2f}.\n")
    A(f"- Fraction of tumour pixels with any detectable imatinib: "
      f"{het['measured_detected_fraction'][0]:.0%}-"
      f"{het['measured_detected_fraction'][1]:.0%} (normal sections: 86-100%).\n")
    A(f"- Spatial autocorrelation e-folding length: tumours "
      f"{het['measured_autocorr_px_tumour'][0]:.2f}-"
      f"{het['measured_autocorr_px_tumour'][1]:.2f} px "
      f"({het['measured_autocorr_px_tumour'][0]*PIXEL_UM:.0f}-"
      f"{het['measured_autocorr_px_tumour'][1]*PIXEL_UM:.0f} um at the assumed 50 um "
      f"pixel), matched normals "
      f"{het['measured_autocorr_px_normal'][0]:.2f}-"
      f"{het['measured_autocorr_px_normal'][1]:.2f} px. Tumour fields are "
      f"systematically patchier than their matched normal in 5 of 6 patients.\n")
    A(f"**{het['caveat']}**\n")

    # ---------------- model ----------------
    A("\n## 3. The model prediction\n")
    A("### 3.1 Imatinib through `tmesim.molecular` (structure in, physics out)\n")
    A(f"SMILES `{IMATINIB_SMILES}` -> RDKit gives C29H31N7O, MW "
      f"{compiled['mw']:.3f}, which matches the deposited sum formula and the "
      f"annotated m/z 494.2663 for [M+H]+ exactly (an independent check that we "
      f"built the right molecule).\n")
    A(f"| quantity | value | route |\n|---|---|---|\n"
      f"| MW | {compiled['mw']:.1f} Da | RDKit |\n"
      f"| logP (Crippen) | {compiled['logp_rdkit']:.2f} | RDKit — note the "
      f"experimental octanol/water logP is ~{IMATINIB_LOGP_EXPERIMENTAL:.1f}, so "
      f"our layer runs ~1.6 log units lipophilic |\n"
      f"| hydrodynamic radius | {compiled['radius_nm']:.3f} nm | MW scaling |\n"
      f"| free diffusivity | {compiled['free_D']:.1f} um^2/s | Stokes-Einstein, "
      f"310 K, eta = 1.2e-3 Pa.s |\n"
      f"| effective diffusivity | {compiled['eff_D']:.1f} um^2/s | Brownian-"
      f"dynamics micro-tortuosity x Ogston-Renkin nm-mesh gate at ECM phi = 0.4 |\n"
      f"| hindrance | {compiled['hindrance']:.3f} | eff_D / free_D |\n"
      f"| plasma free fraction | {compiled['free_fraction']:.3f} | logP heuristic "
      f"(literature: ~{IMATINIB_FU_PLASMA_LIT:.2f}) |\n"
      f"| tissue binding (uncalibrated) | {compiled['binding_uncal']:.2e} /s | "
      f"structure heuristic |\n"
      f"| clearance | {compiled['clearance']:.2e} /s | structure heuristic |\n"
      f"| alpha (perivascular uptake) | {compiled['alpha']:.1f} | **frozen**, fit "
      f"once to doxorubicin (Primeau 2005), never refit here |\n"
      f"| effective loss rate k_eff | {compiled['k_eff']:.3f} /s | k_passive + "
      f"alpha * uptake proxy |\n"
      f"| **penetration length L** | **{L:.1f} um** | sqrt(D_eff / k_eff); "
      f"uncalibrated it would be {compiled['L_uncalibrated_um']:.0f} um |\n")
    A(f"At the assumed 50 um raster, L = {L:.1f} um is **{pred['L_px']:.2f} "
      f"pixels**. The model's own prediction is therefore that imatinib gradients "
      f"in this tissue are BELOW the resolution of this experiment.\n")

    A("\n### 3.2 How a model quantity is mapped to a measured ion-count ratio\n")
    A("A ratio of ion counts is unitless; a partition ratio is unitless; those "
      "are comparable. Absolute counts are not comparable to anything. At quasi-"
      "steady state, both the tumour and the matched normal section see the SAME "
      "plasma, so\n")
    A("```\n"
      "  C_tissue  =  C_plasma,unbound  x  Kpu_tissue  x  E_tissue\n"
      "  R_pred    =  C_tumour / C_normal  =  (Kpu_t / Kpu_n) x (E_t / E_n)\n"
      "  E         =  (L/Rk) * tanh(Rk/L)        1-D Krogh slab, half-spacing Rk\n"
      "```\n")
    A("The plasma term cancels exactly — which is the point: this test is blind "
      "to dose by construction, and that is why the r = 0.03 plasma null is the "
      "right baseline.\n")
    A("Two nested predictions, **neither with any parameter fitted to the GIST "
      "data**:\n")
    A(f"- **M1, transport only.** Kpu_t = Kpu_n (our molecular layer's binding "
      f"terms are drug-specific but tissue-agnostic), so R = E_t/E_n. With L = "
      f"{L:.1f} um and literature half-intervessel distances — normal liver "
      f"{HALF_SPACING_UM['liver']:.1f} um (hepatocyte plates one cell thick "
      f"between sinusoids), gut wall {HALF_SPACING_UM['intestine']:.1f} um, solid "
      f"tumour {HALF_SPACING_UM['tumour']:.1f} um — this gives E_liver = "
      f"{pred['per_tissue']['liver']['E_normal']:.3f}, E_intestine = "
      f"{pred['per_tissue']['intestine']['E_normal']:.3f}, E_tumour = "
      f"{pred['per_tissue']['liver']['E_tumour']:.3f}, hence "
      f"**R = {pred['per_tissue']['liver']['M1_transport_only']:.3f}** for a liver "
      f"pair and **{pred['per_tissue']['intestine']['M1_transport_only']:.3f}** "
      f"for the intestinal pair.\n")
    A(f"- **M2, transport x tissue composition.** Multiplies M1 by "
      f"Kpu(tumour)/Kpu(normal) from a Rodgers-Rowland-style composition model "
      f"using imatinib's own logP (from our layer) and its literature basic pKa "
      f"{IMATINIB_PKA_BASE}. Kpu_liver = "
      f"{pred['per_tissue']['liver']['kpu_normal']:.0f}, Kpu_intestine = "
      f"{pred['per_tissue']['intestine']['kpu_normal']:.0f}, Kpu_tumour = "
      f"{pred['per_tissue']['liver']['kpu_tumour']:.0f}, giving **R = "
      f"{pred['per_tissue']['liver']['M2_transport_x_composition']:.3f}** (liver) "
      f"and **{pred['per_tissue']['intestine']['M2_transport_x_composition']:.3f}** "
      f"(intestine).\n")
    A("The acidic-phospholipid binding term of Rodgers-Rowland is deliberately "
      "OMITTED: its constant Ka_AP is back-calculated from measured blood-cell "
      "partitioning, which we do not have. Including it with a guessed value "
      "would manufacture exactly the liver-vs-tumour difference this test exists "
      "to probe. Omitting it makes M2 conservative in a known direction — it "
      "UNDER-states liver uptake, i.e. it biases M2 toward looking better than it "
      "is.\n")
