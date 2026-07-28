"""THE PERIVASCULAR DRUG GRADIENT on REAL tissue — the core physical prediction,
measured the way the real experiment measures it.

Why this module exists
----------------------
The canonical experimental readout of tumor drug penetration is not a scalar, it is
a **profile**: quantitative immunofluorescence of drug intensity against **distance
from the nearest CD31+ blood vessel**. Primeau et al. (Clin Cancer Res 2005,
PMID 16361566) found doxorubicin intensity falling to half its perivascular value
within ~40-50 um; Patel/Tannock et al. (Cancer Chemother Pharmacol 2013,
PMID 23680920) replicated it (~25-75 um at 10 min, reported jointly for
doxorubicin/mitoxantrone/topotecan).

tmesim predicts a full spatial drug field on **real** tissue, so we can compute the
*same quantity* and compare the **functional form**, not just one number:

1. Reconstruct the REAL 10x Xenium GBM section at the grid-**converged** resolution
   ``dx = 12 um`` (``results/numerical_rigor.md`` section 2: the ~45 um perivascular
   boundary layer is unresolved for dx >~ 15 um, so a coarse grid flattens exactly
   the gradient measured here).
2. Threshold the REAL endothelial-derived ``vascular_density`` field into a vessel
   mask and run ``scipy.ndimage.distance_transform_edt`` (``sampling=dx``) to get
   **distance-from-nearest-vessel in um** at every grid point — the sim analogue of
   the CD31 distance map.
3. Bin the simulated drug **AUC** (and final concentration) by that distance ->
   mean +- inter-quartile spread vs distance.
4. Fit ``C(d) = C0*exp(-d/L)`` and report ``L`` and the half-distance
   ``d_half = L*ln2``; also report the **empirical** half-distance read straight off
   the profile with no functional form assumed.

THE CONFOUND THIS MODULE MUST NOT HIDE
--------------------------------------
The engine's delivery term is ``vascular_density * bbb_penetrance * plasma`` — the
source is **distributed over the whole continuous density field**, not confined to
the thresholded mask. Because the calibrated perivascular sink is fast, the tissue
field is close to locally reaction-dominated (``C ~ S/k``), so the drug field on
real tissue largely **inherits the shape of the source field**. We therefore run
two controls that are as prominent as the headline number:

* **SOURCE NULL** — the identical distance-binning applied to the raw
  ``vascular_density`` field itself. Whatever gradient survives *that* is vessel/KDE
  geometry, not drug transport. The drug/source ratio is reported bin by bin.
* **ISOLATED SINGLE-VESSEL geometry** — the same frozen transport solver, the same
  compiled drug, but one planar (line) vessel in a uniform medium at the real
  tissue's mean ECM. This is the clean perivascular decay the physics predicts with
  the diffuse-source confound removed, and it is where the drugs actually separate.

HONESTY (read before quoting any number)
----------------------------------------
* The vessel mask is a **threshold on an endothelial-cell KDE** (40 um Gaussian
  bandwidth), **not** a measured vessel-lumen segmentation. The bandwidth is the
  same order as the gradient being measured.
* 2D section; straight-line (not geodesic) distances; an untreated resection with no
  measured drug field of its own.
* **Doxorubicin agreement is partly by construction** — the single perivascular-sink
  constant alpha (``validation.calibration``) was fit to the *same* Primeau anchor.
  **Temozolomide's profile is the genuinely independent output.**

Nothing in ``core`` / ``physics`` / ``molecular`` / ``reconstruct`` is modified; this
module only drives them.

Run:    ./.venv/bin/python -m tmesim.validation.vessel_profile
Writes: results/vessel_profile.png
        results/vessel_profile.md

Units: distance um, AUC uM*s, concentration uM.
"""
from __future__ import annotations

import math
import os
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.ndimage import (
    binary_closing,
    binary_fill_holes,
    distance_transform_edt,
    generate_binary_structure,
    label as _label,
    sum as _ndsum,
)
from scipy.optimize import curve_fit

from ..core import (
    CellPopulation,
    Domain,
    Microenvironment,
    PDEConfig,
    SimConfig,
    TissueFields,
)
from ..molecular.binding_site_barrier import apply_barrier_to_drugparams
from ..molecular.library import get_molecule
from ..physics.diffusion import simulate_drug_transport
from ..reconstruct.real_fields import build_real_fields
from ..reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from ..sim.engine import simulate_response
from .calibration import calibrate_alpha, calibrated_drugparams
from .pk_anchors import (
    LN2,
    Q_PERIVASCULAR_HALF,
    anchors_for,
    efold_length_equivalent,
)

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
#: Grid-CONVERGED resolution (results/numerical_rigor.md section 2).
DX_UM = 12.0

#: The two drugs: the calibrated anchor drug, and the opposite-regime contrast.
DRUGS: List[str] = ["doxorubicin", "temozolomide"]

#: Default vessel-mask threshold as a PERCENTILE of the tissue vascular density.
#: p95 leaves ~5% of tissue area as "vessel" — in the range of measured tumor
#: vascular area fractions (~1-6%).
VESSEL_PERCENTILE = 95.0

#: Threshold sensitivity sweep.
VESSEL_PERCENTILE_SWEEP: Sequence[float] = (90.0, 95.0, 98.0, 99.0)

#: Distance-bin width (um); comfortably above dx = 12 um.
BIN_UM = 10.0

#: A distance bin enters a fit only if it holds at least this many grid points.
MIN_BIN_COUNT = 30

#: Near-field window (um) for the pure-exponential fit. The published
#: immunofluorescence profiles are read over roughly this range.
NEAR_WINDOW_UM = 150.0

#: Measured doxorubicin half-distance band (Primeau 2005, PMID 16361566) and the
#: wider replication band (Patel/Tannock 2013, PMID 23680920).
MEASURED_DOX_HALF_UM = (40.0, 50.0)
MEASURED_WIDE_HALF_UM = (25.0, 75.0)

#: Isolated single-vessel geometry: a uniform square domain with one planar
#: (line) vessel down the middle.
ISO_N = 301
ISO_MAX_UM = 900.0

#: The three arms. The cell-derived tissue mask FRAGMENTS as dx shrinks (see
#: :func:`mask_connectivity`), so the nominally-converged dx = 12 um grid solves the
#: PDE on an almost totally disconnected domain and suppresses lateral transport.
#: We therefore run three arms and report all of them.
ARMS: List[Dict] = [
    {"key": "dx12_repaired", "dx_um": 12.0, "repair": 2, "primary": True,
     "label": "dx = 12 µm, morphologically REPAIRED (connected) mask",
     "short": "dx12 repaired"},
    {"key": "dx12_native", "dx_um": 12.0, "repair": 0, "primary": False,
     "label": "dx = 12 µm, NATIVE cell-derived mask (fragmented)",
     "short": "dx12 native"},
    {"key": "dx30_native", "dx_um": 30.0, "repair": 0, "primary": False,
     "label": "dx = 30 µm, NATIVE mask (connected, but boundary layer coarse)",
     "short": "dx30 native"},
]

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "vessel_profile.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "vessel_profile.md")

_PDE = dict(dt=10.0, total_time=3600.0 * 6, save_every=500)


# --------------------------------------------------------------------------- #
# 1. Vessel mask + distance transform (the sim analogue of the CD31 distance map)
# --------------------------------------------------------------------------- #
def vessel_mask_from_density(
    vascular_density: np.ndarray,
    tissue_mask: np.ndarray,
    percentile: float = VESSEL_PERCENTILE,
) -> Dict:
    """Threshold the REAL endothelial-derived vascular-density field into a binary
    vessel mask.

    The threshold is the ``percentile``-th percentile of the density *over tissue*,
    so it adapts to the field's own normalisation rather than to an arbitrary
    absolute level. Returns the boolean ``mask``, the absolute ``threshold``, and the
    resulting ``area_fraction`` of tissue.

    NOT a measured vessel lumen when the field is the smooth KDE: there
    ``vascular_density`` is a 40 um-bandwidth Gaussian KDE of measured
    ENDOTHELIAL cell positions (``reconstruct.real_fields``, ``mode='kde'``).

    SPARSE-SOURCE PATH (the default field since ``results/vessel_promotion.md``).
    The discrete source is exactly zero over ~97.5% of the tissue, so its
    ``percentile``-th percentile is 0.0 and ``vasc >= 0`` would select the WHOLE
    tissue as "vessel" — silently destroying the distance transform. When the
    threshold comes out non-positive we therefore fall back to ``vasc > 0``,
    which for a segmented field IS the vessel set (this is
    ``vessel_network_impact.md`` §8.4 item 3: the mask becomes the segmented mask
    rather than a percentile on the field). The KDE path is untouched — its p95
    is strictly positive — so this changes no existing number.
    """
    vasc = np.asarray(vascular_density, dtype=np.float64)
    tissue = np.asarray(tissue_mask, dtype=bool)
    if not tissue.any():
        raise ValueError("empty tissue mask")
    thr = float(np.percentile(vasc[tissue], float(percentile)))
    if thr > 0.0:
        mask = (vasc >= thr) & tissue
        rule = f"p{float(percentile):g} threshold on the field"
    else:
        mask = (vasc > 0.0) & tissue
        rule = "field > 0 (sparse/segmented source; p-threshold degenerate at 0)"
    return {
        "mask": mask,
        "threshold": thr,
        "percentile": float(percentile),
        "rule": rule,
        "area_fraction": float(mask.sum()) / float(tissue.sum()),
        "n_vessel_px": int(mask.sum()),
    }


def distance_from_vessel_um(vessel_mask: np.ndarray, dx_um: float) -> np.ndarray:
    """Euclidean distance (um) from every grid point to the NEAREST vessel pixel.

    ``distance_transform_edt`` on the *complement* of the mask with ``sampling =
    dx_um``, so the result is already in micrometers; points inside the mask get 0
    (the "perivascular" bin). Straight-line, not geodesic.
    """
    mask = np.asarray(vessel_mask, dtype=bool)
    if not mask.any():
        raise ValueError("vessel mask is empty; lower the threshold")
    return distance_transform_edt(~mask, sampling=float(dx_um)).astype(np.float64)


# --------------------------------------------------------------------------- #
# 2. Binning: signal vs distance-from-vessel (the measured quantity)
# --------------------------------------------------------------------------- #
def distance_profile(
    field: np.ndarray,
    distance_um: np.ndarray,
    tissue_mask: np.ndarray,
    bin_um: float = BIN_UM,
    max_distance_um: Optional[float] = None,
) -> Dict:
    """Bin ``field`` by distance-from-nearest-vessel over the tissue mask.

    Returns per-bin centres with mean / median / std / 25th / 75th percentile and
    point count — the mean +- spread profile the immunofluorescence experiment
    reports.
    """
    f = np.asarray(field, dtype=np.float64)
    d = np.asarray(distance_um, dtype=np.float64)
    tissue = np.asarray(tissue_mask, dtype=bool)
    fv, dv = f[tissue], d[tissue]

    if max_distance_um is None:
        max_distance_um = float(np.percentile(dv, 99.5))
    n_bins = max(int(math.ceil(max_distance_um / float(bin_um))), 1)
    edges = np.arange(n_bins + 1, dtype=np.float64) * float(bin_um)
    idx = np.digitize(dv, edges) - 1

    keys = ("distance_um", "mean", "median", "std", "q25", "q75")
    acc: Dict[str, List[float]] = {k: [] for k in keys}
    count: List[int] = []
    for b in range(n_bins):
        sel = idx == b
        n = int(sel.sum())
        if n == 0:
            continue
        vals = fv[sel]
        acc["distance_um"].append(0.5 * (edges[b] + edges[b + 1]))
        acc["mean"].append(float(vals.mean()))
        acc["median"].append(float(np.median(vals)))
        acc["std"].append(float(vals.std()))
        acc["q25"].append(float(np.percentile(vals, 25.0)))
        acc["q75"].append(float(np.percentile(vals, 75.0)))
        count.append(n)

    out = {k: np.asarray(v, dtype=np.float64) for k, v in acc.items()}
    out["count"] = np.asarray(count, dtype=int)
    out["bin_um"] = float(bin_um)
    return out


def usable_bins(profile: Dict, max_um: Optional[float] = None) -> Tuple[np.ndarray, ...]:
    """``(distance, mean, count)`` for bins with enough points (and, optionally,
    within ``max_um``)."""
    use = profile["count"] >= MIN_BIN_COUNT
    if max_um is not None:
        use &= profile["distance_um"] <= float(max_um)
    return profile["distance_um"][use], profile["mean"][use], profile["count"][use]


# --------------------------------------------------------------------------- #
# 3. Fits and shape tests
# --------------------------------------------------------------------------- #
def fit_exponential(
    distance_um: np.ndarray,
    value: np.ndarray,
    weight: Optional[np.ndarray] = None,
) -> Dict:
    """Weighted least squares for ``log(value) = log(C0) - d/L``.

    Weights default to per-bin point counts. Returns ``C0``, ``L_um``, the derived
    ``d_half_um = L*ln2``, the weighted R^2 in **log space** (the exponentiality
    score) and in linear space.
    """
    d, v, w = _clean(distance_um, value, weight)
    if d.size < 3:
        return {"error": "fewer than 3 usable bins"}
    y = np.log(v)
    slope, intercept = np.polyfit(d, y, 1, w=np.sqrt(w))
    if slope >= 0:
        return {"error": f"fitted slope {slope:.3g} >= 0 (profile not decaying)"}
    L = -1.0 / float(slope)
    C0 = float(np.exp(intercept))
    return {
        "C0": C0,
        "L_um": float(L),
        "d_half_um": float(L * LN2),
        "r2_log": _r2(y, intercept + slope * d, w),
        "r2_linear_space": _r2(v, C0 * np.exp(-d / L), w),
        "n_bins": int(d.size),
        "fit_range_um": (float(d.min()), float(d.max())),
    }


def fit_exponential_plus_floor(
    distance_um: np.ndarray,
    value: np.ndarray,
    weight: Optional[np.ndarray] = None,
) -> Dict:
    """Fit ``C(d) = C0*exp(-d/L) + C_floor`` in LINEAR space.

    On real tissue the distributed vascular source puts a non-zero background under
    the perivascular decay; a pure exponential over the full range is then dominated
    by that plateau and reports a meaningless L. This three-parameter form separates
    the decaying perivascular component (``L``, ``d_half``) from the floor, and also
    reports the half-distance of the **total** profile
    ``d_half_total`` (where ``C0*exp(-d/L)+floor`` reaches half of ``C0+floor``),
    which is what an experimentalist reading a raw curve would quote.
    """
    d, v, w = _clean(distance_um, value, weight)
    if d.size < 4:
        return {"error": "fewer than 4 usable bins"}

    def model(x, c0, L, floor):
        return c0 * np.exp(-x / L) + floor

    floor0 = float(np.min(v))
    p0 = [max(float(v[0]) - floor0, 1e-12), 40.0, floor0]
    try:
        popt, _ = curve_fit(
            model, d, v, p0=p0, sigma=1.0 / np.sqrt(w), absolute_sigma=False,
            bounds=([0.0, 1.0, 0.0], [np.inf, 1e5, np.inf]), maxfev=40000,
        )
    except (RuntimeError, ValueError) as exc:
        return {"error": f"curve_fit failed: {exc}"}
    c0, L, floor = (float(x) for x in popt)
    total0 = c0 + floor
    # Half of the TOTAL perivascular signal, if the floor allows it to be reached.
    half_total = 0.5 * total0
    d_half_total = float(-L * math.log((half_total - floor) / c0)) if (
        c0 > 0 and half_total > floor
    ) else None
    return {
        "C0": c0,
        "L_um": L,
        "floor": floor,
        "floor_fraction": float(floor / total0) if total0 > 0 else float("nan"),
        "d_half_um": float(L * LN2),
        "d_half_total_um": d_half_total,
        "r2_linear_space": _r2(v, model(d, *popt), w),
        "n_bins": int(d.size),
        "fit_range_um": (float(d.min()), float(d.max())),
    }


def fit_linear_shape(
    distance_um: np.ndarray, value: np.ndarray, weight: Optional[np.ndarray] = None
) -> Dict:
    """Straight-line fit ``C(d) = a + b*d`` — the null model the exponential has to
    beat before any claim about "exponential shape" means anything."""
    d, v, w = _clean(distance_um, value, weight)
    b, a = np.polyfit(d, v, 1, w=np.sqrt(w))
    return {"a": float(a), "b": float(b),
            "r2_linear_space": _r2(v, a + b * d, w)}


def empirical_half_distance_um(
    distance_um: np.ndarray, value: np.ndarray
) -> Optional[float]:
    """The measured quantity with NO functional form assumed: the distance at which
    the profile first falls to HALF of its perivascular (first-bin) value, by linear
    interpolation between bracketing bins. ``None`` if it never halves."""
    d = np.asarray(distance_um, dtype=np.float64)
    v = np.asarray(value, dtype=np.float64)
    if d.size < 2 or v[0] <= 0:
        return None
    target = 0.5 * float(v[0])
    below = np.nonzero(v <= target)[0]
    if below.size == 0:
        return None
    i = int(below[0])
    if i == 0:
        return float(d[0])
    v0, v1 = float(v[i - 1]), float(v[i])
    if v0 == v1:
        return float(d[i])
    return float(d[i - 1] + (v0 - target) / (v0 - v1) * (d[i] - d[i - 1]))


def _clean(distance_um, value, weight):
    d = np.asarray(distance_um, dtype=np.float64)
    v = np.asarray(value, dtype=np.float64)
    ok = np.isfinite(d) & np.isfinite(v) & (v > 0.0)
    w = (np.ones_like(d) if weight is None
         else np.asarray(weight, dtype=np.float64))
    return d[ok], v[ok], w[ok]


def _r2(y, yhat, w) -> float:
    ybar = float(np.average(y, weights=w))
    ss_tot = float(np.sum(w * (y - ybar) ** 2))
    if ss_tot <= 0:
        return float("nan")
    return float(1.0 - float(np.sum(w * (y - yhat) ** 2)) / ss_tot)


def profile_fits(profile: Dict) -> Dict:
    """The standard bundle of fits/readouts applied to one distance profile."""
    d_all, m_all, w_all = usable_bins(profile)
    d_near, m_near, w_near = usable_bins(profile, NEAR_WINDOW_UM)
    return {
        "profile": profile,
        "near_exp": fit_exponential(d_near, m_near, w_near),
        "full_exp": fit_exponential(d_all, m_all, w_all),
        "exp_floor": fit_exponential_plus_floor(d_all, m_all, w_all),
        "linear_null": fit_linear_shape(d_near, m_near, w_near),
        "empirical_d_half_um": empirical_half_distance_um(d_all, m_all),
        "perivascular_mean": float(m_all[0]) if m_all.size else float("nan"),
        "plateau_fraction": (float(np.median(m_all[d_all > 200.0]) / m_all[0])
                             if (m_all.size and (d_all > 200.0).any()) else float("nan")),
    }


# --------------------------------------------------------------------------- #
# 4. Drive the real simulation
# --------------------------------------------------------------------------- #
def compile_real_drug(name: str, alpha: float):
    """The same compilation pipeline the flagship real runs use: alpha-calibrated
    perivascular sink + saturable binding-site barrier."""
    dp = calibrated_drugparams(name, alpha)
    mol = get_molecule(name)
    return apply_barrier_to_drugparams(dp, mol), mol


def analyze_drug(micro, name: str, alpha: float, distance_um: np.ndarray,
                 pde: PDEConfig) -> Dict:
    """Simulate ONE drug on the real tissue and reduce its AUC / final-concentration
    fields to distance-from-vessel profiles with fits."""
    dp, mol = compile_real_drug(name, alpha)
    res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)

    auc = np.asarray(res.fields["auc"], dtype=np.float64)
    conc = np.asarray(res.transport.final_concentration, dtype=np.float64)
    k = float(dp.binding) + float(dp.clearance)
    return {
        "drug": name,
        "mw": float(mol.get("mw")),
        "diffusivity": float(dp.diffusivity),
        "binding": float(dp.binding),
        "clearance": float(dp.clearance),
        "k_loss": k,
        "L_analytic_um": float(math.sqrt(dp.diffusivity / k)) if k > 0 else float("inf"),
        "bbb": float(dp.bbb_penetrance),
        "kill_fraction": float(res.kill_fraction),
        "penetration_index": float(res.extras["penetration_index"]),
        "mean_auc": float(res.extras["mean_auc"]),
        "auc_field": auc,
        "auc": profile_fits(distance_profile(auc, distance_um, tissue)),
        "conc": profile_fits(distance_profile(conc, distance_um, tissue)),
    }


# --------------------------------------------------------------------------- #
# 5. CONTROL: the isolated single-vessel geometry (diffuse source removed)
# --------------------------------------------------------------------------- #
def _uniform_micro(ecm_density: float, n: int = ISO_N, dx_um: float = DX_UM):
    """A uniform square Microenvironment (no cells to kill) at a fixed ECM density —
    used only as a clean substrate for the frozen transport solver."""
    dom = Domain(nx=n, ny=n, dx=float(dx_um))
    fields = TissueFields(
        ecm_density=np.full((n, n), float(ecm_density)),
        vascular_density=np.zeros((n, n)),
    )
    cells = CellPopulation(
        x=np.array([0.0]), y=np.array([0.0]), cell_type=np.array([0]),
        state=np.array([0]), sensitivity=np.array([1.0]),
    )
    return Microenvironment(domain=dom, fields=fields, cells=cells,
                            patient_id="isolated_vessel")


def isolated_vessel_profile(name: str, alpha: float, ecm_density: float,
                            dx_um: float = DX_UM) -> Dict:
    """The SAME frozen transport solver and the SAME compiled drug, but with a single
    planar (line) vessel in a uniform medium at the real tissue's mean ECM.

    This removes the distributed-source confound: the only thing shaping the profile
    is transport + the calibrated perivascular sink. For a planar source the steady
    solution is exactly ``C(d) = C0*exp(-d/L)`` with ``L = sqrt(D_eff/k)``, so this is
    the cleanest possible statement of "the perivascular gradient this drug's physics
    predicts", and the geometry the calibration constant was defined in.
    """
    dp, _ = compile_real_drug(name, alpha)
    micro = _uniform_micro(ecm_density, dx_um=dx_um)
    n = micro.domain.nx
    src = np.zeros((n, n), dtype=np.float64)
    src[:, n // 2] = 1.0
    pde = PDEConfig(convection=False, **_PDE)
    tr = simulate_drug_transport(micro, dp, src, pde)

    auc = np.asarray(tr.auc, dtype=np.float64)
    dist = np.abs(np.arange(n) - n // 2) * float(dx_um)
    keep = dist <= ISO_MAX_UM
    # Average the two sides; the geometry is symmetric about the vessel.
    curve = auc.mean(axis=0)[keep]
    d = dist[keep]
    order = np.argsort(d)
    d_s, c_s = d[order], curve[order]
    # Collapse the +/- pair at each distance.
    uniq = np.unique(d_s)
    mean = np.array([c_s[d_s == u].mean() for u in uniq])
    fit = fit_exponential(uniq[uniq > 0], mean[uniq > 0])
    return {
        "drug": name,
        "distance_um": uniq,
        "mean": mean,
        "fit": fit,
        "empirical_d_half_um": empirical_half_distance_um(uniq, mean),
        "L_analytic_um": float(math.sqrt(
            dp.diffusivity * max(1.0 - 0.6 * ecm_density, 1e-6)
            / (dp.binding + dp.clearance))),
    }


# --------------------------------------------------------------------------- #
# 6. Threshold sensitivity
# --------------------------------------------------------------------------- #
def threshold_sensitivity(micro, auc_field: np.ndarray,
                          percentiles: Sequence[float] = VESSEL_PERCENTILE_SWEEP
                          ) -> List[Dict]:
    """How far does the answer move if the (arbitrary) vessel threshold moves?
    Re-derives the mask, the distance map and the fits at each percentile."""
    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    rows = []
    for p in percentiles:
        vm = vessel_mask_from_density(vasc, tissue, percentile=p)
        dist = distance_from_vessel_um(vm["mask"], micro.domain.dx)
        pf = profile_fits(distance_profile(auc_field, dist, tissue))
        vpf = profile_fits(distance_profile(vasc, dist, tissue))
        rows.append({
            "percentile": float(p),
            "threshold": vm["threshold"],
            "area_fraction": vm["area_fraction"],
            "max_distance_um": float(dist[tissue].max()),
            "near_L_um": pf["near_exp"].get("L_um", float("nan")),
            "near_d_half_um": pf["near_exp"].get("d_half_um", float("nan")),
            "empirical_d_half_um": pf["empirical_d_half_um"],
            "source_empirical_d_half_um": vpf["empirical_d_half_um"],
            "r2_log_near": pf["near_exp"].get("r2_log", float("nan")),
        })
    return rows


# --------------------------------------------------------------------------- #
# 6b. Tissue-mask connectivity — the artefact that forced the three-arm design
# --------------------------------------------------------------------------- #
def mask_connectivity(tissue_mask: np.ndarray) -> Dict:
    """Connected-component statistics of a tissue mask.

    The reconstruction rasterises measured cell positions, so as ``dx`` shrinks the
    mask stops being a filled sheet and becomes a scatter of isolated pixels. A
    disconnected mask is fatal for a transport study: the no-flux boundary is applied
    at every mask edge, so an isolated pixel has **no neighbours to diffuse into**
    and its concentration collapses to the purely local ``C = S/k``. Any
    distance-from-vessel profile then measures the *source* field, not transport.
    """
    t = np.asarray(tissue_mask, dtype=bool)
    lab, n = _label(t)
    if n == 0:
        return {"n_components": 0, "largest_fraction": 0.0, "tissue_px": 0,
                "box_fraction": 0.0}
    sizes = _ndsum(t, lab, range(1, n + 1))
    return {
        "n_components": int(n),
        "largest_fraction": float(np.max(sizes) / t.sum()),
        "tissue_px": int(t.sum()),
        "box_fraction": float(t.sum() / t.size),
    }


def repair_tissue_mask(tissue_mask: np.ndarray, iterations: int = 2) -> np.ndarray:
    """Morphologically close (then hole-fill) a fragmented cell-derived tissue mask so
    the transport domain is connected.

    Pure image morphology on a *copy* — ``reconstruct`` is untouched. The added pixels
    are the interstitial gaps *between* measured cells, which are physically tissue
    the drug can diffuse through; the native mask excludes them only because the
    rasteriser had no cell centroid there.
    """
    t = np.asarray(tissue_mask, dtype=bool)
    if iterations <= 0:
        return t
    st = generate_binary_structure(2, 2)
    return binary_fill_holes(binary_closing(t, structure=st, iterations=int(iterations)))


def build_arm_microenvironment(dx_um: float, repair: int) -> Tuple:
    """Reconstruct the real section at ``dx_um`` and, if ``repair > 0``, rebuild it on
    a morphologically repaired tissue mask.

    When the mask is repaired the ECM / vascular fields are **re-derived on the new
    mask** with the same public ``reconstruct.real_fields.build_real_fields`` (the KDE
    is multiplied by the tissue mask at build time, so reusing the old fields would
    leave the added pixels with a spurious zero source).
    """
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
    native = np.asarray(micro.domain.tissue_mask, dtype=bool)
    conn_native = mask_connectivity(native)
    if repair <= 0:
        return micro, conn_native, conn_native
    repaired = repair_tissue_mask(native, repair)
    dom = Domain(nx=micro.domain.nx, ny=micro.domain.ny, dx=micro.domain.dx,
                 tissue_mask=repaired)
    fields = build_real_fields(dom, micro.cells)
    micro2 = Microenvironment(domain=dom, fields=fields, cells=micro.cells,
                              patient_id=micro.patient_id, metadata=micro.metadata)
    return micro2, conn_native, mask_connectivity(repaired)


# --------------------------------------------------------------------------- #
# 6c. KDE-bandwidth sweep — is the fitted decay OUR OWN smoothing?
# --------------------------------------------------------------------------- #
#: Bandwidths (um) for the endothelial-KDE sensitivity sweep. 40 um is the default
#: in ``reconstruct.real_fields``.
BANDWIDTH_SWEEP_UM: Sequence[float] = (25.0, 40.0, 60.0)


def bandwidth_sensitivity(
    dx_um: float,
    repair: int,
    alpha: float,
    drugs: Sequence[str],
    bandwidths: Sequence[float] = BANDWIDTH_SWEEP_UM,
) -> List[Dict]:
    """Rebuild the REAL fields at several endothelial-KDE bandwidths and re-measure the
    gradient.

    This is the decisive test of whether the fitted decay length is drug transport or
    **our own smoothing kernel**. If the fitted `L` tracks the bandwidth roughly 1:1
    and the drug/source-null ratio stays ~1, the spatial fit on reconstructed tissue is
    measuring the KDE, not the pharmacology — and must not be quoted as a validation.
    """
    base = reconstruct_real(source="xenium", dx_um=float(dx_um))
    native = np.asarray(base.domain.tissue_mask, dtype=bool)
    mask = repair_tissue_mask(native, repair) if repair > 0 else native
    dom = Domain(nx=base.domain.nx, ny=base.domain.ny, dx=base.domain.dx,
                 tissue_mask=mask)
    pde = PDEConfig(convection=True, **_PDE)

    rows = []
    for bw in bandwidths:
        fields = build_real_fields(dom, base.cells, bandwidth_um=float(bw))
        micro = Microenvironment(domain=dom, fields=fields, cells=base.cells,
                                 patient_id=base.patient_id, metadata=base.metadata)
        vasc = np.asarray(fields.vascular_density, dtype=np.float64)
        vm = vessel_mask_from_density(vasc, mask, VESSEL_PERCENTILE)
        dist = distance_from_vessel_um(vm["mask"], dom.dx)
        null = profile_fits(distance_profile(vasc, dist, mask))
        null_L = null["exp_floor"].get("L_um", float("nan"))
        entry = {"bandwidth_um": float(bw),
                 "source_null_L_um": null_L,
                 "source_null_d_half_um": null["empirical_d_half_um"],
                 "drugs": {}}
        for name in drugs:
            r = analyze_drug(micro, name, alpha, dist, pde)
            Lf = r["auc"]["exp_floor"].get("L_um", float("nan"))
            entry["drugs"][name] = {
                "L_um": Lf,
                "d_half_um": r["auc"]["empirical_d_half_um"],
                "ratio_to_null": (Lf / null_L) if null_L else float("nan"),
                "corr_with_source": float(np.corrcoef(
                    r["auc_field"][mask], vasc[mask])[0, 1]),
            }
        rows.append(entry)
        parts = "  ".join(f"{n}: L={entry['drugs'][n]['L_um']:.1f} "
                          f"({entry['drugs'][n]['ratio_to_null']:.2f}x null)"
                          for n in drugs)
        print(f"  KDE bandwidth {bw:4.0f} um -> source-null L = {null_L:6.1f} um   {parts}")
    return rows


# --------------------------------------------------------------------------- #
# 7. Verdict
# --------------------------------------------------------------------------- #
def verdict_vs_measurement(d_half_um: Optional[float]) -> Dict:
    """Compare a predicted half-distance to the measured 40-50 um band (Primeau 2005)
    and the wider 25-75 um replication band (Patel/Tannock 2013)."""
    lo, hi = MEASURED_DOX_HALF_UM
    wide_lo, wide_hi = MEASURED_WIDE_HALF_UM
    if d_half_um is None or not np.isfinite(d_half_um):
        return {"d_half_um": None, "in_band": False, "in_wide_band": False,
                "ratio_to_mid": float("nan"),
                "verdict": "NO DECAY (no half-distance could be measured)"}
    mid = 0.5 * (lo + hi)
    in_band = bool(lo <= d_half_um <= hi)
    in_wide = bool(wide_lo <= d_half_um <= wide_hi)
    ratio = float(d_half_um) / mid
    if in_band:
        v = "PASS (inside the measured 40-50 um band)"
    elif in_wide:
        v = "PASS(wide) (outside 40-50 um, inside the 25-75 um replication band)"
    elif 1 / 3 <= ratio <= 3:
        v = "PASS(3x) (outside both bands, same order of magnitude)"
    else:
        v = "FAIL (more than 3x off the measured band)"
    return {"d_half_um": float(d_half_um), "in_band": in_band,
            "in_wide_band": in_wide, "ratio_to_mid": ratio, "verdict": v}


def verdict_L_vs_efold_band(L_um: Optional[float]) -> Dict:
    """Compare a predicted **e-folding length** L to the *L-equivalent* of the
    measured half-distance band.

    ``pk_anchors`` is explicit that the published 40–50 µm is a HALF-distance, so the
    quantity directly comparable to an e-folding length is ``40–50 / ln2 ≈ 58–72 µm``.
    This is the stricter (and more correct) of the two comparisons; the existing
    ``pk_reality`` check compares L to the raw 40–50 µm, which is deliberately
    conservative.
    """
    prim = next((a for a in anchors_for("doxorubicin")
                 if a["quantity"] == Q_PERIVASCULAR_HALF and a["pmid"] == "16361566"),
                None)
    band = efold_length_equivalent(prim) if prim else None
    if band is None or L_um is None or not np.isfinite(L_um):
        return {"L_um": None, "band_um": band, "in_band": False,
                "ratio_to_mid": float("nan"), "verdict": "n/a"}
    lo, hi = band
    mid = 0.5 * (lo + hi)
    ratio = float(L_um) / mid
    in_band = bool(lo <= L_um <= hi)
    if in_band:
        v = f"PASS (inside the L-equivalent {lo:.0f}-{hi:.0f} um band)"
    elif 1 / 3 <= ratio <= 3:
        v = (f"PASS(3x) ({'steeper' if ratio < 1 else 'shallower'} than measured by "
             f"{abs(1 - ratio) * 100:.0f}%, same order of magnitude)")
    else:
        v = "FAIL (more than 3x off the L-equivalent band)"
    return {"L_um": float(L_um), "band_um": (float(lo), float(hi)),
            "in_band": in_band, "ratio_to_mid": ratio, "verdict": v}


# --------------------------------------------------------------------------- #
# 8. Run
# --------------------------------------------------------------------------- #
def run_arm(arm: Dict, alpha: float, drugs: Sequence[str]) -> Dict:
    """Run ONE arm: reconstruct (optionally repairing the mask), build the vessel mask
    and distance map, simulate every drug, and reduce everything to profiles + fits."""
    print(f"\n=== ARM {arm['key']}: {arm['label']} ===")
    micro, conn_native, conn_used = build_arm_microenvironment(arm["dx_um"], arm["repair"])
    report = reconstruction_report(micro)
    print(f"  grid {report['grid']}  dx={report['dx_um']:.0f} um  tissue "
          f"{conn_used['tissue_px']:,} px ({conn_used['box_fraction']*100:.0f}% of box)  "
          f"components={conn_used['n_components']:,}  largest="
          f"{conn_used['largest_fraction']*100:.1f}% of tissue")

    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    vm = vessel_mask_from_density(vasc, tissue, VESSEL_PERCENTILE)
    dist = distance_from_vessel_um(vm["mask"], micro.domain.dx)
    print(f"  vessel mask: density >= {vm['threshold']:.3f} (p{VESSEL_PERCENTILE:.0f}) "
          f"-> {vm['area_fraction']*100:.1f}% of tissue area; distance median "
          f"{np.median(dist[tissue]):.0f} um, max {dist[tissue].max():.0f} um")

    source_null = profile_fits(distance_profile(vasc, dist, tissue))
    print(f"  SOURCE NULL (vascular_density itself): empirical d_half = "
          f"{source_null['empirical_d_half_um']:.1f} um")

    pde = PDEConfig(convection=True, **_PDE)
    results = []
    for name in drugs:
        print(f"  simulating {name:14s} ...", end="", flush=True)
        r = analyze_drug(micro, name, alpha, dist, pde)
        r["corr_with_source"] = float(np.corrcoef(r["auc_field"][tissue], vasc[tissue])[0, 1])
        ef = r["auc"]["exp_floor"]
        print(f" d_half(emp)={r['auc']['empirical_d_half_um']:6.1f} um  "
              f"L(exp+floor)={ef.get('L_um', float('nan')):6.1f} um  "
              f"R2={ef.get('r2_linear_space', float('nan')):.4f}  "
              f"corr(AUC,source)={r['corr_with_source']:.3f}")
        results.append(r)

    by_name = {r["drug"]: r for r in results}
    dox = by_name.get("doxorubicin")
    return {
        **arm,
        "micro": micro,
        "report": report,
        "connectivity_native": conn_native,
        "connectivity_used": conn_used,
        "vessel_mask": vm,
        "distance_um": dist,
        "source_null": source_null,
        "results": results,
        "by_name": by_name,
        "verdict_empirical": verdict_vs_measurement(
            dox["auc"]["empirical_d_half_um"] if dox else None),
        "verdict_L": verdict_L_vs_efold_band(
            dox["auc"]["exp_floor"].get("L_um") if dox else None),
    }


def run(drugs: Sequence[str] = tuple(DRUGS), arms: Sequence[Dict] = tuple(ARMS)) -> Dict:
    """Full pipeline over all three arms + the isolated single-vessel control."""
    t0 = time.time()
    alpha = calibrate_alpha()
    print(f"perivascular-sink alpha = {alpha:.1f} (fit ONCE on the doxorubicin anchor)")

    # Connectivity-vs-dx diagnostic: WHY there are three arms.
    print("\nTissue-mask connectivity vs grid resolution (the cell-derived mask "
          "fragments as dx shrinks):")
    conn_sweep = []
    for dx in (60.0, 40.0, 30.0, 20.0, 15.0, 12.0):
        m = reconstruct_real(source="xenium", dx_um=dx)
        c = mask_connectivity(np.asarray(m.domain.tissue_mask, dtype=bool))
        c["dx_um"] = dx
        conn_sweep.append(c)
        print(f"  dx={dx:5.1f} um  tissue {c['tissue_px']:7,} px "
              f"({c['box_fraction']*100:4.1f}% of box)  components {c['n_components']:6,}  "
              f"largest {c['largest_fraction']*100:5.1f}% of tissue")

    arm_results = [run_arm(a, alpha, drugs) for a in arms]
    by_key = {a["key"]: a for a in arm_results}
    primary = next(a for a in arm_results if a.get("primary"))

    # Isolated single-vessel control (uniform medium; the drug physics alone).
    mean_ecm = float(primary["report"]["mean_ecm_density"])
    print(f"\nISOLATED single-vessel control (uniform medium, ECM={mean_ecm:.3f}):")
    isolated = []
    for name in drugs:
        iso = isolated_vessel_profile(name, alpha, mean_ecm, dx_um=DX_UM)
        f = iso["fit"]
        print(f"  {name:14s} L={f.get('L_um', float('nan')):7.1f} um  "
              f"d_half={f.get('d_half_um', float('nan')):7.1f} um  "
              f"R2(log)={f.get('r2_log', float('nan')):.5f}")
        isolated.append(iso)
    iso_by_name = {i["drug"]: i for i in isolated}

    print("\nThreshold sensitivity on the primary arm ...")
    sens = threshold_sensitivity(primary["micro"],
                                 primary["by_name"]["doxorubicin"]["auc_field"])

    print("\nKDE-BANDWIDTH sweep — is the fitted decay our own smoothing kernel?")
    print("  (primary arm geometry; the vascular/ECM fields are rebuilt at each "
          "bandwidth)")
    bw_primary = bandwidth_sensitivity(primary["dx_um"], primary["repair"], alpha, drugs)
    print("  ... and on the FRAGMENTED native dx = 12 um mask, for contrast:")
    bw_native = bandwidth_sensitivity(12.0, 0, alpha, drugs)

    return {
        "alpha": alpha,
        "connectivity_sweep": conn_sweep,
        "arms": arm_results,
        "by_key": by_key,
        "primary": primary,
        "isolated": isolated,
        "iso_by_name": iso_by_name,
        "verdict_isolated": verdict_vs_measurement(
            iso_by_name["doxorubicin"]["fit"].get("d_half_um")
            if "doxorubicin" in iso_by_name else None),
        "verdict_isolated_L": verdict_L_vs_efold_band(
            iso_by_name["doxorubicin"]["fit"].get("L_um")
            if "doxorubicin" in iso_by_name else None),
        "threshold_sensitivity": sens,
        "bandwidth_primary": bw_primary,
        "bandwidth_native": bw_native,
        "runtime_s": time.time() - t0,
    }


# --------------------------------------------------------------------------- #
# 9. Figure
# --------------------------------------------------------------------------- #
_COLORS = {"doxorubicin": "#d62728", "temozolomide": "#1f77b4"}
_SRC_COLOR = "#6a6a6a"


def make_figure(res: Dict, path: str = _OUT_PNG) -> str:
    """Nine panels. Row 1: the real tissue + vessel mask and the two simulated AUC
    maps. Row 2: the measured quantity (drug vs distance) with the source null and the
    40-50 um band, the same log-linear, and the isolated single-vessel control.
    Row 3: the mask-fragmentation artefact across arms, the summary of every fitted
    half-distance against the measured bands, and the connectivity-vs-dx diagnostic.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.gridspec import GridSpec

    prim = res["primary"]
    micro = prim["micro"]
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    vmask = prim["vessel_mask"]["mask"]
    W, H = micro.domain.extent_um
    ext = [0, W, H, 0]
    rows = prim["results"]
    lo, hi = MEASURED_DOX_HALF_UM
    wlo, whi = MEASURED_WIDE_HALF_UM

    fig = plt.figure(figsize=(20.5, 17.0))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1.0, 0.85, 0.85],
                  hspace=0.34, wspace=0.24)

    # --- (a) real tissue + vessel mask ------------------------------------- #
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(np.where(tissue, 1.0, np.nan), cmap="Greys", extent=ext, aspect="equal",
              vmin=0, vmax=6.0, interpolation="nearest")
    ax.imshow(np.where(tissue & ~vmask, vasc, np.nan), cmap="Blues", extent=ext,
              aspect="equal", vmin=0.0, vmax=0.45, interpolation="nearest")
    overlay = np.zeros(vmask.shape + (4,), dtype=float)
    overlay[vmask] = (0.85, 0.05, 0.05, 1.0)
    ax.imshow(overlay, extent=ext, aspect="equal", interpolation="nearest")
    cu = prim["connectivity_used"]
    ax.set_title(
        "(a) REAL tissue (grey) + endothelial-derived vascular density (blue)\n"
        f"+ VESSEL MASK (red) = density ≥ {prim['vessel_mask']['threshold']:.2f} "
        f"(p{prim['vessel_mask']['percentile']:.0f}) → "
        f"{prim['vessel_mask']['area_fraction']*100:.1f}% of tissue area\n"
        f"[{prim['short']}: largest connected component "
        f"{cu['largest_fraction']*100:.0f}% of tissue]", fontsize=9.5)
    ax.set_xticks([]); ax.set_yticks([])

    # --- (b), (c) simulated AUC maps ---------------------------------------- #
    for j, r in enumerate(rows[:2]):
        axm = fig.add_subplot(gs[0, 1 + j])
        f = np.where(tissue, r["auc_field"], np.nan)
        vmax = float(np.nanpercentile(f, 99.8)) or 1.0
        vmin = max(float(np.nanpercentile(f, 2.0)), vmax * 1e-4, 1e-30)
        im = axm.imshow(f, cmap="magma", extent=ext, aspect="equal",
                        norm=LogNorm(vmin=vmin, vmax=vmax), interpolation="nearest")
        axm.contour(vmask.astype(float), levels=[0.5], colors="#00e5ff",
                    linewidths=0.35, extent=ext, origin="upper")
        ef = r["auc"]["exp_floor"]
        axm.set_title(
            f"({'bc'[j]}) simulated {r['drug']} exposure AUC (µM·s, log scale)\n"
            f"empirical d½ = {r['auc']['empirical_d_half_um']:.0f} µm  ·  "
            f"L = {ef.get('L_um', float('nan')):.0f} µm  ·  far-field plateau "
            f"{r['auc']['plateau_fraction']*100:.0f}%\n(cyan = vessel mask)",
            fontsize=9.5)
        axm.set_xlim(0, W); axm.set_ylim(H, 0)
        axm.set_xticks([]); axm.set_yticks([])
        fig.colorbar(im, ax=axm, shrink=0.78, pad=0.02)

    # --- (d) THE measured quantity, primary arm ----------------------------- #
    axd = fig.add_subplot(gs[1, 0])
    axd.axvspan(lo, hi, color="#ffb000", alpha=0.35, zorder=0,
                label="measured doxorubicin d½ 40–50 µm\n(Primeau 2005, PMID 16361566)")
    axd.axvspan(wlo, whi, color="#ffb000", alpha=0.13, zorder=0,
                label="replication band 25–75 µm\n(Patel/Tannock 2013, PMID 23680920)")
    axd.axhline(0.5, color="0.35", lw=0.9, ls=":", zorder=1)
    for r in rows:
        p = r["auc"]["profile"]
        d, m, _ = usable_bins(p)
        use = p["count"] >= MIN_BIN_COUNT
        c0, col = float(m[0]), _COLORS.get(r["drug"], "0.3")
        axd.fill_between(d, p["q25"][use] / c0, p["q75"][use] / c0, color=col,
                         alpha=0.15, lw=0)
        axd.plot(d, m / c0, "o-", color=col, ms=3.0, lw=1.5,
                 label=f"{r['drug']} — sim mean ± IQR "
                       f"(d½={r['auc']['empirical_d_half_um']:.0f} µm)")
        ef = r["auc"]["exp_floor"]
        if "L_um" in ef:
            dd = np.linspace(0, 500, 300)
            axd.plot(dd, (ef["C0"] * np.exp(-dd / ef["L_um"]) + ef["floor"]) / c0,
                     "--", color=col, lw=1.2,
                     label=f"    fit C₀e^(−d/L)+floor: L={ef['L_um']:.0f} µm, "
                           f"R²={ef['r2_linear_space']:.3f}")
    sn = prim["source_null"]
    ds, ms, _ = usable_bins(sn["profile"])
    axd.plot(ds, ms / float(ms[0]), "-", color=_SRC_COLOR, lw=2.2, alpha=0.85, zorder=1,
             label=f"SOURCE NULL: vascular_density itself "
                   f"(d½={sn['empirical_d_half_um']:.0f} µm)")
    axd.set_xlim(0, 500); axd.set_ylim(0, 1.05)
    axd.set_xlabel("distance from nearest blood vessel (µm)")
    axd.set_ylabel("drug AUC / perivascular AUC")
    axd.set_title(f"(d) THE MEASURED QUANTITY — {prim['short']}\n"
                  "drug exposure vs distance from the nearest real vessel",
                  fontsize=10.5)
    axd.legend(fontsize=7.0, loc="upper right", framealpha=0.93)
    axd.grid(alpha=0.25)

    # --- (e) log-linear ------------------------------------------------------ #
    axl = fig.add_subplot(gs[1, 1])
    axl.axvspan(lo, hi, color="#ffb000", alpha=0.35, zorder=0)
    for r in rows:
        d, m, _ = usable_bins(r["auc"]["profile"])
        col = _COLORS.get(r["drug"], "0.3")
        axl.semilogy(d, m / float(m[0]), "o", color=col, ms=3.0, label=r["drug"])
        ef = r["auc"]["exp_floor"]
        if "L_um" in ef:
            dd = np.linspace(0, 500, 300)
            axl.semilogy(dd, (ef["C0"] * np.exp(-dd / ef["L_um"]) + ef["floor"])
                         / float(m[0]), "--", color=col, lw=1.2,
                         label=f"  exp+floor (R²={ef['r2_linear_space']:.3f})")
    axl.semilogy(ds, ms / float(ms[0]), "-", color=_SRC_COLOR, lw=2.0, alpha=0.85,
                 label="SOURCE NULL")
    axl.set_xlim(0, 500)
    axl.set_xlabel("distance from nearest vessel (µm)")
    axl.set_ylabel("normalised signal (log)")
    axl.set_title("(e) log-linear — a straight line is a true exponential.\n"
                  "Real tissue: exponential decay onto a distributed-source floor.",
                  fontsize=10.5)
    axl.legend(fontsize=7.5, loc="lower left"); axl.grid(alpha=0.25, which="both")

    # --- (f) isolated single-vessel control ---------------------------------- #
    axi = fig.add_subplot(gs[1, 2])
    axi.axvspan(lo, hi, color="#ffb000", alpha=0.35, zorder=0,
                label="measured doxorubicin d½ 40–50 µm")
    for iso in res["isolated"]:
        d, m = iso["distance_um"], iso["mean"]
        col = _COLORS.get(iso["drug"], "0.3")
        axi.semilogy(d, m / float(m[0]), "o", color=col, ms=3.0,
                     label=f"{iso['drug']} (isolated vessel)")
        fit = iso["fit"]
        if "L_um" in fit:
            dd = np.linspace(0, float(d.max()), 300)
            axi.semilogy(dd, fit["C0"] * np.exp(-dd / fit["L_um"]) / float(m[0]), "--",
                         color=col, lw=1.2,
                         label=f"  L={fit['L_um']:.0f} µm, d½={fit['d_half_um']:.0f} µm, "
                               f"R²={fit['r2_log']:.4f}")
    axi.set_xlim(0, ISO_MAX_UM); axi.set_ylim(1e-6, 1.5)
    axi.set_xlabel("distance from the single vessel (µm)")
    axi.set_ylabel("AUC / perivascular AUC (log)")
    axi.set_title("(f) CONTROL — ONE vessel, uniform medium, same solver.\n"
                  "Diffuse-source confound removed: a clean exponential,\n"
                  "and the drugs separate.", fontsize=10.5)
    axi.legend(fontsize=7.5, loc="upper right"); axi.grid(alpha=0.25, which="both")

    # --- (g) the mask-fragmentation artefact, across arms -------------------- #
    axg = fig.add_subplot(gs[2, 0])
    axg.axvspan(lo, hi, color="#ffb000", alpha=0.35, zorder=0)
    styles = {"dx12_repaired": "-", "dx12_native": "--", "dx30_native": ":"}
    for arm in res["arms"]:
        ls = styles.get(arm["key"], "-")
        for r in arm["results"]:
            d, m, _ = usable_bins(r["auc"]["profile"])
            axg.semilogy(d, m / float(m[0]), ls, color=_COLORS.get(r["drug"], "0.3"),
                         lw=1.6, alpha=0.9,
                         label=f"{r['drug'][:3]} · {arm['short']} "
                               f"(d½={r['auc']['empirical_d_half_um']:.0f} µm)")
        ds2, ms2, _ = usable_bins(arm["source_null"]["profile"])
        axg.semilogy(ds2, ms2 / float(ms2[0]), ls, color=_SRC_COLOR, lw=1.2, alpha=0.6,
                     label=f"source null · {arm['short']}")
    axg.set_xlim(0, 400); axg.set_ylim(5e-2, 1.3)
    axg.set_xlabel("distance from nearest vessel (µm)")
    axg.set_ylabel("normalised signal (log)")
    axg.set_title("(g) THE NUMERICAL ARTEFACT. At dx = 12 µm the cell-derived\n"
                  "tissue mask is 99.3% fragmented → no lateral diffusion →\n"
                  "both drugs collapse onto the source null.", fontsize=10)
    axg.legend(fontsize=6.2, loc="lower left", ncol=2, framealpha=0.9)
    axg.grid(alpha=0.25, which="both")

    # --- (h) every half-distance vs the measured bands ----------------------- #
    axh = fig.add_subplot(gs[2, 1])
    axh.axvspan(wlo, whi, color="#ffb000", alpha=0.15,
                label="25–75 µm (PMID 23680920)")
    axh.axvspan(lo, hi, color="#ffb000", alpha=0.40, label="40–50 µm (PMID 16361566)")
    labels, vals, cols = [], [], []
    for arm in res["arms"]:
        for r in arm["results"]:
            labels.append(f"{r['drug'][:12]} · {arm['short']}")
            vals.append(r["auc"]["empirical_d_half_um"])
            cols.append(_COLORS.get(r["drug"], "0.3"))
        labels.append(f"source null · {arm['short']}")
        vals.append(arm["source_null"]["empirical_d_half_um"])
        cols.append(_SRC_COLOR)
    for iso in res["isolated"]:
        labels.append(f"{iso['drug'][:12]} · ISOLATED vessel")
        vals.append(iso["fit"].get("d_half_um", float("nan")))
        cols.append(_COLORS.get(iso["drug"], "0.3"))
    for b in res.get("bandwidth_primary", []):
        for name, dd in b["drugs"].items():
            labels.append(f"{name[:12]} · KDE bw {b['bandwidth_um']:.0f} µm")
            vals.append(dd["d_half_um"])
            cols.append(_COLORS.get(name, "0.3"))
        labels.append(f"source null · KDE bw {b['bandwidth_um']:.0f} µm")
        vals.append(b["source_null_d_half_um"])
        cols.append(_SRC_COLOR)
    y = np.arange(len(labels))
    axh.barh(y, vals, color=cols, alpha=0.85, height=0.68)
    for yi, v in zip(y, vals):
        axh.text(v * 1.06, yi, f"{v:.0f}", va="center", fontsize=7.0)
    axh.set_yticks(y); axh.set_yticklabels(labels, fontsize=7.0)
    axh.invert_yaxis(); axh.set_xscale("log"); axh.set_xlim(5, 900)
    axh.set_xlabel("half-distance d½ (µm, log)")
    axh.set_title("(h) every predicted half-distance vs the measured bands\n"
                  "(doxorubicin is the CALIBRATION anchor — its agreement is\n"
                  "partly by construction; temozolomide is the free prediction)",
                  fontsize=10)
    axh.legend(fontsize=7.0, loc="lower right"); axh.grid(alpha=0.25, axis="x", which="both")

    # --- (i) connectivity vs dx ---------------------------------------------- #
    axc = fig.add_subplot(gs[2, 2])
    cs = res["connectivity_sweep"]
    dxs = [c["dx_um"] for c in cs]
    axc.plot(dxs, [c["largest_fraction"] * 100 for c in cs], "o-", color="#2ca02c",
             lw=1.8, label="largest connected component (% of tissue)")
    axc.plot(dxs, [c["box_fraction"] * 100 for c in cs], "s--", color="#8c564b",
             lw=1.4, label="tissue as % of bounding box")
    for arm in res["arms"]:
        if arm["repair"]:
            axc.plot([arm["dx_um"]], [arm["connectivity_used"]["largest_fraction"] * 100],
                     "*", color="#d62728", ms=16, zorder=5,
                     label=f"after morphological repair ({arm['short']})")
    axc.axvspan(10, 15, color="#d62728", alpha=0.10)
    axc.set_xlabel("grid spacing dx (µm)"); axc.set_ylabel("percent")
    axc.set_ylim(0, 105); axc.invert_xaxis()
    axc.set_title("(i) WHY three arms: the cell-derived tissue mask\n"
                  "fragments as dx shrinks. results/numerical_rigor.md asks for\n"
                  "dx ≲ 12 µm — exactly where the domain falls apart.", fontsize=10)
    axc.legend(fontsize=7.2, loc="lower left"); axc.grid(alpha=0.25)

    vp_, vl_ = prim["verdict_empirical"], prim["verdict_L"]
    vi = res["verdict_isolated"]
    fig.suptitle(
        "PERIVASCULAR DRUG GRADIENT on the REAL 10x Xenium glioblastoma section — "
        "the same quantity the immunofluorescence experiments measure\n"
        f"PRIMARY ARM: {prim['label']}  ·  doxorubicin d½ = "
        f"{vp_['d_half_um']:.0f} µm → {vp_['verdict']}  ·  isolated single vessel d½ = "
        f"{vi['d_half_um']:.0f} µm\n"
        "doxorubicin is the CALIBRATION anchor (agreement partly by construction); "
        "TEMOZOLOMIDE is the independent, unvalidated prediction",
        fontsize=12)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=115, bbox_inches="tight")
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# 10. Report
# --------------------------------------------------------------------------- #
def _f(x, nd=1):
    try:
        v = float(x)
        return "n/a" if not np.isfinite(v) else f"{v:.{nd}f}"
    except (TypeError, ValueError):
        return "n/a"


def write_report(res: Dict, path: str = _OUT_MD) -> str:
    """Write results/vessel_profile.md."""
    prim = res["primary"]
    rep, vm, sn = prim["report"], prim["vessel_mask"], prim["source_null"]
    rows = prim["results"]
    dox, tmz = prim["by_name"].get("doxorubicin"), prim["by_name"].get("temozolomide")
    iso_dox = res["iso_by_name"].get("doxorubicin")
    iso_tmz = res["iso_by_name"].get("temozolomide")
    vt, vl = prim["verdict_empirical"], prim["verdict_L"]
    vi, vil = res["verdict_isolated"], res["verdict_isolated_L"]
    native = res["by_key"]["dx12_native"]
    dx30 = res["by_key"]["dx30_native"]
    dist_t = prim["distance_um"][np.asarray(prim["micro"].domain.tissue_mask, dtype=bool)]

    prim_a = next((a for a in anchors_for("doxorubicin")
                   if a["quantity"] == Q_PERIVASCULAR_HALF and a["pmid"] == "16361566"), None)
    patel = next((a for a in anchors_for("doxorubicin")
                  if a["quantity"] == Q_PERIVASCULAR_HALF and a["pmid"] == "23680920"), None)

    def dhalf(arm, drug):
        r = arm["by_name"].get(drug)
        return r["auc"]["empirical_d_half_um"] if r else float("nan")

    def Lof(arm, drug):
        r = arm["by_name"].get(drug)
        return r["auc"]["exp_floor"].get("L_um", float("nan")) if r else float("nan")

    L: List[str] = [
        "# Perivascular drug gradient on the REAL tissue — predicted vs measured",
        "",
        "**Module:** `tmesim/validation/vessel_profile.py` "
        "(run: `./.venv/bin/python -m tmesim.validation.vessel_profile`)",
        "",
        "**Figure:** `results/vessel_profile.png`",
        "",
        "The canonical experimental readout of tumor drug penetration is a **profile**, "
        "not a scalar: quantitative immunofluorescence of drug intensity against "
        "**distance from the nearest CD31⁺ blood vessel**. tmesim predicts a full "
        "spatial drug field on **real** tissue, so we can compute *the same quantity* "
        "and compare the **functional form**.",
        "",
        "## Headline",
        "",
        f"| | doxorubicin | temozolomide |",
        "|---|---:|---:|",
        f"| **fitted L** (µm), real tissue, primary arm | **{_f(Lof(prim,'doxorubicin'))}** "
        f"| **{_f(Lof(prim,'temozolomide'))}** |",
        f"| **half-distance d½** (µm), real tissue, primary arm | "
        f"**{_f(dhalf(prim,'doxorubicin'))}** | **{_f(dhalf(prim,'temozolomide'))}** |",
        f"| fitted L (µm), isolated single vessel | "
        f"{_f(iso_dox['fit'].get('L_um')) if iso_dox else 'n/a'} | "
        f"{_f(iso_tmz['fit'].get('L_um')) if iso_tmz else 'n/a'} |",
        f"| d½ = L·ln2 (µm), isolated single vessel | "
        f"{_f(iso_dox['fit'].get('d_half_um')) if iso_dox else 'n/a'} | "
        f"{_f(iso_tmz['fit'].get('d_half_um')) if iso_tmz else 'n/a'} |",
        f"| measured | **40–50 µm** (PMID 16361566), 25–75 µm (PMID 23680920) | "
        "*no published perivascular-distance dataset* |",
        "",
        "### Read this before quoting any number",
        "",
        "Two different quantities are easy to conflate, and only one of them is clean:",
        "",
        "1. **The analytic / single-vessel penetration length** "
        "`L = √(D_eff / (clearance + binding))`, which is genuinely drug-specific "
        f"(doxorubicin {_f(iso_dox['L_analytic_um']) if iso_dox else 'n/a'} µm vs "
        f"temozolomide {_f(iso_tmz['L_analytic_um']) if iso_tmz else 'n/a'} µm) and is "
        "what the isolated single-vessel control reproduces with the real solver. "
        "**This is a legitimate model output and is the thing to compare to the "
        "literature.**",
        "2. **The spatial fit to drug-vs-distance on our reconstructed tissue.** This one "
        "is **confounded**: the vessel mask and the drug source are built from the *same* "
        "endothelial-cell KDE, so part of the fitted decay is our own 40 µm smoothing "
        f"bandwidth. The source null (d½ = {_f(sn['empirical_d_half_um'])} µm with **no "
        "drug physics at all**) and the bandwidth sweep below quantify how much. On the "
        "fragmented native dx = 12 µm mask the confound is *total* — every drug returns "
        "the same ~26–28 µm decay, within 1.1–1.2× of the null, reproducing the "
        "independent finding in `results/delivery_metrics.md`. **A match to the measured "
        "40–50 µm band read off the spatial fit is not a validation.** Quote the "
        "drug/source-null ratio instead.",
        "",
        f"**Shape:** exponential. On the isolated single vessel the predicted profile is "
        f"a pure exponential to five decimals "
        f"(log-space R² = {_f(iso_dox['fit'].get('r2_log'), 5) if iso_dox else 'n/a'}). "
        "On real tissue it is **exponential decaying onto a floor**, "
        f"`C(d) = C₀·exp(−d/L) + floor` (R² = "
        f"{_f(dox['auc']['exp_floor'].get('r2_linear_space'), 4) if dox else 'n/a'}), "
        "the floor being the drug delivered by the *distributed* vascular source. "
        "Published profiles also flatten into a background, so the qualitative form "
        "matches.",
        "",
        f"**Verdict for doxorubicin:** real tissue d½ = {_f(vt['d_half_um'])} µm → "
        f"*{vt['verdict']}*; the fitted `L = {_f(vl.get('L_um'))} µm` against the "
        "L-equivalent band `40–50/ln2 ≈ 58–72 µm` → "
        f"*{vl['verdict']}*. **This is not an independent validation** — α was fit to "
        "this same anchor. **Temozolomide "
        f"(L = {_f(Lof(prim,'temozolomide'))} µm, "
        f"{Lof(prim,'temozolomide')/Lof(prim,'doxorubicin'):.1f}× doxorubicin's) is the "
        "free prediction**, and it has no measurement to be scored against.",
        "",
        "---",
        "",
        "## Method",
        "",
        f"1. **Real geometry.** The open 10x Xenium GBM section "
        f"({rep['n_cells']:,} measured cells; tumor {rep['n_tumor']:,}) reconstructed by "
        f"`reconstruct.real_pipeline.reconstruct_real('xenium', dx_um=…)` → grid "
        f"{rep['grid']}, extent {rep['extent_um'][0]:.0f}×{rep['extent_um'][1]:.0f} µm, "
        f"mean ECM {rep['mean_ecm_density']:.3f}, mean vascular density "
        f"{rep['mean_vascular_density']:.3f}. **dx = 12 µm** is the grid-converged "
        "resolution `results/numerical_rigor.md` §2 asks for (the ~45 µm perivascular "
        "boundary layer is unresolved for dx ≳ 15 µm).",
        f"2. **Vessel mask.** The REAL endothelial-derived `vascular_density` field "
        f"thresholded at its {vm['percentile']:.0f}th tissue percentile "
        f"(`≥ {vm['threshold']:.3f}`) → **{vm['n_vessel_px']:,} vessel pixels = "
        f"{vm['area_fraction']*100:.1f}% of tissue area**, in the range of measured "
        "tumor vascular area fractions.",
        f"3. **Distance map.** `scipy.ndimage.distance_transform_edt` on the complement "
        f"of the mask with `sampling = dx` → distance-from-nearest-vessel in µm at "
        f"every grid point (over tissue: median {np.median(dist_t):.0f} µm, p95 "
        f"{np.percentile(dist_t, 95):.0f} µm, max {dist_t.max():.0f} µm).",
        f"4. **Simulation.** Each drug compiled through the standard pipeline "
        f"(`calibrate_alpha()` → α = {res['alpha']:.1f}, fit **once** on doxorubicin → "
        "`calibrated_drugparams` → `apply_barrier_to_drugparams`), run on the same "
        "tissue with `PDEConfig(dt=10 s, total_time=6 h, convection=True)`.",
        f"5. **Profile + fits.** AUC (and final concentration) binned in {BIN_UM:.0f} µm "
        f"distance bins (≥ {MIN_BIN_COUNT} grid points per bin), mean ± IQR. Readouts: "
        "the **empirical** half-distance (no functional form assumed), a "
        "`C₀·exp(−d/L) + floor` fit over the full range (the correct model on real "
        f"tissue), a pure exponential over the near field (d ≤ {NEAR_WINDOW_UM:.0f} µm), "
        "and a straight-line null.",
        "6. **Two controls:** the **source null** (identical binning applied to the raw "
        "`vascular_density` field) and the **isolated single-vessel** geometry.",
        "",
        "---",
        "",
        "## A numerical artefact we had to deal with first (please read)",
        "",
        "The reconstruction rasterises measured **cell positions**, so as `dx` shrinks "
        "the tissue mask stops being a filled sheet and becomes a scatter of isolated "
        "pixels. The transport solver applies a no-flux boundary at every mask edge, so "
        "an isolated pixel has **no neighbour to diffuse into** and its concentration "
        "collapses to the purely local `C = S/k`. Measured directly:",
        "",
        "| dx (µm) | tissue px | tissue as % of bounding box | connected components | "
        "largest component (% of tissue) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for c in res["connectivity_sweep"]:
        L.append(f"| {c['dx_um']:.0f} | {c['tissue_px']:,} | {c['box_fraction']*100:.1f}% "
                 f"| {c['n_components']:,} | {c['largest_fraction']*100:.1f}% |")
    L += [
        "",
        "**At dx = 12 µm — exactly the resolution `numerical_rigor.md` recommends — the "
        "domain is 26,058 disconnected fragments whose largest component holds 0.7% of "
        "the tissue.** Lateral diffusion is effectively switched off. That is why, on "
        "the native dx = 12 µm mask, doxorubicin and temozolomide produce *the same* "
        f"profile ({_f(dhalf(native,'doxorubicin'))} vs "
        f"{_f(dhalf(native,'temozolomide'))} µm) despite ~4× different transport "
        "lengths, and both sit on top of the source null "
        f"({_f(native['source_null']['empirical_d_half_um'])} µm) — correlation between "
        f"the AUC field and the source field {native['by_name']['doxorubicin']['corr_with_source']:.3f} "
        f"and {native['by_name']['temozolomide']['corr_with_source']:.3f}. This is a "
        "grid artefact, not physics.",
        "",
        "We therefore report **three arms** and take as primary the one that resolves "
        "the boundary layer *and* has a connected domain:",
        "",
        "| arm | dx (µm) | largest connected component | doxorubicin d½ (µm) | "
        "doxorubicin L (µm) | temozolomide d½ (µm) | temozolomide L (µm) | "
        "source-null d½ (µm) | corr(AUC, source): dox / TMZ |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for arm in res["arms"]:
        star = " **(PRIMARY)**" if arm.get("primary") else ""
        L.append(
            f"| {arm['label']}{star} | {arm['dx_um']:.0f} | "
            f"{arm['connectivity_used']['largest_fraction']*100:.1f}% | "
            f"**{_f(dhalf(arm,'doxorubicin'))}** | {_f(Lof(arm,'doxorubicin'))} | "
            f"**{_f(dhalf(arm,'temozolomide'))}** | {_f(Lof(arm,'temozolomide'))} | "
            f"*{_f(arm['source_null']['empirical_d_half_um'])}* | "
            f"{arm['by_name']['doxorubicin']['corr_with_source']:.3f} / "
            f"{arm['by_name']['temozolomide']['corr_with_source']:.3f} |")
    L += [
        "",
        "The **repair** is pure image morphology applied to a *copy* of the mask "
        "(`scipy.ndimage.binary_closing`, 2 iterations, then `binary_fill_holes`), which "
        "restores the largest component to "
        f"{prim['connectivity_used']['largest_fraction']*100:.1f}% of tissue; the ECM and "
        "vascular fields are then **re-derived on the repaired mask** with the public "
        "`reconstruct.real_fields.build_real_fields`, because the KDE is multiplied by "
        "the tissue mask at build time and reusing the old fields would leave the added "
        "pixels with a spurious zero source. **Nothing in `reconstruct` is modified.** "
        "The added pixels are the interstitial gaps *between* measured cells — real "
        "tissue the drug diffuses through, excluded by the native mask only because the "
        "rasteriser found no cell centroid there.",
        "",
        f"The dx = 30 µm arm (native mask, {dx30['connectivity_used']['largest_fraction']*100:.0f}% "
        "connected) is an independent cross-check with no morphology applied at all, and "
        "it agrees with the repaired arm qualitatively: doxorubicin "
        f"{_f(dhalf(dx30,'doxorubicin'))} µm vs temozolomide "
        f"{_f(dhalf(dx30,'temozolomide'))} µm. **Honest range across the two "
        "physically-sane arms: doxorubicin d½ ≈ "
        f"{min(dhalf(prim,'doxorubicin'), dhalf(dx30,'doxorubicin')):.0f}–"
        f"{max(dhalf(prim,'doxorubicin'), dhalf(dx30,'doxorubicin')):.0f} µm, "
        f"temozolomide d½ ≈ {min(dhalf(prim,'temozolomide'), dhalf(dx30,'temozolomide')):.0f}–"
        f"{max(dhalf(prim,'temozolomide'), dhalf(dx30,'temozolomide')):.0f} µm.**",
        "",
        "---",
        "",
        "## Result 1 — the gradient on real tissue (primary arm)",
        "",
        "**Model-free readout** (what the experiment reports — the distance at which the "
        "binned mean falls to half its perivascular value, no functional form assumed):",
        "",
        "| signal | perivascular → far-field plateau | **empirical d½** (µm) | "
        "corr(AUC, source field) |",
        "|---|---:|---:|---:|",
    ]
    for r in rows:
        L.append(f"| **{r['drug']}** AUC | 1.00 → {r['auc']['plateau_fraction']:.2f} | "
                 f"**{_f(r['auc']['empirical_d_half_um'])}** | "
                 f"{r['corr_with_source']:.3f} |")
    L.append(f"| *SOURCE NULL: `vascular_density` itself* | 1.00 → "
             f"{sn['plateau_fraction']:.2f} | *{_f(sn['empirical_d_half_um'])}* | 1.000 |")
    L += [
        "",
        "**Fits.** A *pure* exponential is the wrong model on real tissue: the "
        "distributed source puts a plateau under the perivascular decay and drags a "
        "single-exponential fit towards it. The three-parameter "
        "`C(d) = C₀·exp(−d/L) + floor` fits far better and its `L` is the real "
        "perivascular decay length. All R² are reported in **linear space** so the "
        "models are directly comparable:",
        "",
        "| signal | model | **L** (µm) | **d½ = L·ln2** (µm) | floor (% of perivascular) "
        "| **R² (linear space)** | R² (log space) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for r in rows + [{"drug": "SOURCE NULL", "auc": sn}]:
        a, nm = r["auc"], r["drug"]
        ef, ne, ln = a["exp_floor"], a["near_exp"], a["linear_null"]
        L += [
            f"| **{nm}** | `C₀·exp(−d/L) + floor` (full range) | "
            f"**{_f(ef.get('L_um'))}** | **{_f((ef.get('L_um') or float('nan')) * LN2)}** | "
            f"{_f((ef.get('floor_fraction', float('nan')) or 0) * 100)}% | "
            f"**{_f(ef.get('r2_linear_space'), 4)}** | — |",
            f"| {nm} | pure `C₀·exp(−d/L)`, near field d ≤ {NEAR_WINDOW_UM:.0f} µm | "
            f"{_f(ne.get('L_um'))} | {_f(ne.get('d_half_um'))} | — | "
            f"{_f(ne.get('r2_linear_space'), 4)} | {_f(ne.get('r2_log'), 4)} |",
            f"| {nm} | straight-line null `a + b·d` (near field) | — | — | — | "
            f"{_f(ln['r2_linear_space'], 4)} | — |",
        ]
    L += [
        "",
        "The **final-concentration** field gives essentially the same profile "
        f"(doxorubicin empirical d½ {_f(dox['conc']['empirical_d_half_um']) if dox else 'n/a'} µm, "
        f"temozolomide {_f(tmz['conc']['empirical_d_half_um']) if tmz else 'n/a'} µm) — "
        "expected, because with a constant-source window the field is quasi-steady, so "
        "`AUC ≈ C_steady · Δt` and the two share a spatial shape. It is *not* an "
        "independent check, but it does rule out the profile being an artefact of "
        "time-integration.",
        "",
        "### How much of the gradient is the drug, and how much is the vasculature?",
        "",
        "The engine's delivery term is `vascular_density · bbb_penetrance · plasma`: the "
        "source is spread over the **whole continuous density field**, not confined to "
        "the thresholded mask (mean vascular density over tissue = "
        f"{rep['mean_vascular_density']:.3f}; tissue 200 µm from the mask still carries "
        f"{sn['plateau_fraction']*100:.0f}% of the peak density). The SOURCE NULL row "
        "measures what that alone buys: an apparent half-distance of "
        f"**{_f(sn['empirical_d_half_um'])} µm** with no drug physics at all. On the "
        f"primary arm the drugs sit clearly above it "
        f"({_f(dhalf(prim,'doxorubicin'))} and {_f(dhalf(prim,'temozolomide'))} µm, "
        f"correlations {dox['corr_with_source']:.2f} and {tmz['corr_with_source']:.2f}), "
        "so transport is doing real work — but a meaningful part of the profile is still "
        "vessel/KDE geometry, and that share must be subtracted mentally before reading "
        "the agreement with the literature as evidence.",
        "",
        "---",
        "",
        "## Result 2 — CONTROL: the isolated single-vessel gradient",
        "",
        "To see the gradient the *physics* predicts, free of the distributed-source "
        "confound, we run the **same frozen transport solver** and the **same compiled "
        "drug** on one planar (line) vessel in a uniform medium at the real tissue's mean "
        f"ECM ({rep['mean_ecm_density']:.3f}). For a planar source the steady solution is "
        "exactly `C(d) = C₀·exp(−d/L)` with `L = √(D_eff/k)` — the geometry the "
        "calibration constant is defined in.",
        "",
        "| drug | MW (Da) | D (µm²/s) | k = clearance+binding (1/s) | **L** (µm) | "
        "**d½ = L·ln2** (µm) | R² (log) |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        fit = res["iso_by_name"].get(r["drug"], {}).get("fit", {})
        L.append(f"| **{r['drug']}** | {r['mw']:.0f} | {r['diffusivity']:.0f} | "
                 f"{r['k_loss']:.3g} | **{_f(fit.get('L_um'))}** | "
                 f"**{_f(fit.get('d_half_um'))}** | {_f(fit.get('r2_log'), 5)} |")
    if iso_dox and iso_tmz and iso_dox["fit"].get("L_um"):
        ratio = iso_tmz["fit"]["L_um"] / iso_dox["fit"]["L_um"]
        L += [
            "",
            f"**Temozolomide's gradient is {ratio:.1f}× longer than doxorubicin's** in "
            "identical geometry — the physically expected ordering: temozolomide is "
            f"small (MW {tmz['mw']:.0f} vs {dox['mw']:.0f} Da), diffuses faster "
            f"({tmz['diffusivity']:.0f} vs {dox['diffusivity']:.0f} µm²/s) and is not "
            f"consumed by a high-capacity intracellular/DNA sink (k = {tmz['k_loss']:.2g} "
            f"vs {dox['k_loss']:.2g} /s), so it escapes the perivascular cuff and reaches "
            "tissue the anthracycline never sees.",
            "",
            "This control is also where the exponential-shape claim is unambiguous: "
            "log-space R² ≈ 1.0 for both drugs, i.e. the predicted profile is a **true "
            "exponential**, the functional form the measurements report.",
        ]

    L += [
        "",
        "---",
        "",
        "## Comparison to the real measurements",
        "",
        "**The anchors (both PMID-verified, carried in `tmesim/validation/pk_anchors.py`):**",
        "",
    ]
    if prim_a:
        L += [f"> {prim_a['citation']} **PMID {prim_a['pmid']}** — doxorubicin intensity "
              f"falls to half its perivascular value within "
              f"**{prim_a['value'][0]:.0f}–{prim_a['value'][1]:.0f} µm** of CD31⁺ "
              f"vessels ({prim_a['method']}; {prim_a['model']}).", ""]
    if patel:
        L += [f"> {patel['citation']} **PMID {patel['pmid']}** — independent "
              f"replication; half-distance **{patel['value'][0]:.0f}–"
              f"{patel['value'][1]:.0f} µm** at 10 min, reported jointly for "
              "doxorubicin / mitoxantrone / topotecan.", ""]
    L += [
        "| doxorubicin readout | predicted | measured | verdict |",
        "|---|---:|---:|:--|",
        f"| real tissue (primary arm), **empirical d½** | **{_f(vt['d_half_um'])} µm** | "
        f"40–50 µm (25–75 µm) | {vt['verdict']} |",
        f"| real tissue (primary arm), fitted **L** vs L-equivalent band | "
        f"**{_f(vl.get('L_um'))} µm** | 58–72 µm | {vl['verdict']} |",
        f"| real tissue (dx = 30 µm cross-check), empirical d½ | "
        f"{_f(dhalf(dx30,'doxorubicin'))} µm | 40–50 µm (25–75 µm) | "
        f"{verdict_vs_measurement(dhalf(dx30,'doxorubicin'))['verdict']} |",
        f"| **isolated single vessel**, d½ | **{_f(vi['d_half_um'])} µm** | "
        f"40–50 µm (25–75 µm) | {vi['verdict']} |",
        f"| isolated single vessel, L vs L-equivalent band | {_f(vil.get('L_um'))} µm | "
        f"58–72 µm | {vil['verdict']} |",
        "",
        "**Why two bands.** `pk_anchors.py` is explicit that the published 40–50 µm is a "
        "*half*-distance, so the quantity comparable to an e-folding length is "
        "`40–50 / ln2 ≈ 58–72 µm`. Comparing an `L` to the raw 40–50 µm (as "
        "`pk_reality.py` does) is deliberately conservative. Both comparisons are given "
        "above so no reader has to guess which convention a number is in.",
        "",
        "### Is the shape exponential, as measured?",
        "",
    ]
    for r in rows:
        a = r["auc"]
        ne, ln, ef = a["near_exp"], a["linear_null"], a["exp_floor"]
        fit = res["iso_by_name"].get(r["drug"], {}).get("fit", {})
        L.append(
            f"- **{r['drug']}** — *isolated vessel*: log-space R² = "
            f"**{_f(fit.get('r2_log'), 5)}**, a clean exponential to five decimals. "
            f"*Real tissue*: **yes, exponential-plus-floor** — R²(linear) = "
            f"**{_f(ef.get('r2_linear_space'), 4)}** for `C₀·exp(−d/L)+floor`, versus "
            f"{_f(ne.get('r2_linear_space'), 4)} for a pure exponential and "
            f"{_f(ln['r2_linear_space'], 4)} for a straight line over the same near "
            f"field. The decaying component has L = {_f(ef.get('L_um'))} µm; beyond "
            f"~150 µm the profile plateaus at {a['plateau_fraction']*100:.0f}% of the "
            "perivascular value instead of continuing to decay.")
    L += [
        "",
        "So the predicted functional form **is** the measured one: exponential where the "
        "geometry is clean, exponential-onto-a-background on real tissue. The published "
        "profiles also flatten into a background at large distance — but ours is a "
        "modelling consequence of the diffuse source, not a measured autofluorescence "
        "floor, and the two should not be equated.",
        "",
        "---",
        "",
        "## THE DECISIVE CONTROL — is the fitted decay our own KDE bandwidth?",
        "",
        "`vascular_density` is a Gaussian KDE of endothelial-cell positions with a "
        "**40 µm bandwidth**, and the vessel mask is a threshold on that *same* smoothed "
        "field. If the fitted decay length simply tracks the bandwidth, the spatial fit "
        "is measuring our own smoothing kernel and **cannot** be quoted as agreement with "
        "the literature. We therefore rebuild the real fields at 25 / 40 / 60 µm "
        "bandwidth and re-measure everything.",
        "",
        "**Primary arm (dx = 12 µm, repaired connected mask):**",
        "",
        "| KDE bandwidth (µm) | source-null L (µm) | doxorubicin L (µm) | dox / null | "
        "temozolomide L (µm) | TMZ / null | corr(AUC, source): dox / TMZ |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for b in res["bandwidth_primary"]:
        dd, tt = b["drugs"].get("doxorubicin", {}), b["drugs"].get("temozolomide", {})
        L.append(f"| {b['bandwidth_um']:.0f} | {_f(b['source_null_L_um'])} | "
                 f"**{_f(dd.get('L_um'))}** | {_f(dd.get('ratio_to_null'), 2)}× | "
                 f"**{_f(tt.get('L_um'))}** | {_f(tt.get('ratio_to_null'), 2)}× | "
                 f"{_f(dd.get('corr_with_source'), 3)} / {_f(tt.get('corr_with_source'), 3)} |")
    L += [
        "",
        "**Fragmented native dx = 12 µm mask (the arm with no lateral transport), for "
        "contrast:**",
        "",
        "| KDE bandwidth (µm) | source-null L (µm) | doxorubicin L (µm) | dox / null | "
        "temozolomide L (µm) | TMZ / null |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for b in res["bandwidth_native"]:
        dd, tt = b["drugs"].get("doxorubicin", {}), b["drugs"].get("temozolomide", {})
        L.append(f"| {b['bandwidth_um']:.0f} | {_f(b['source_null_L_um'])} | "
                 f"{_f(dd.get('L_um'))} | {_f(dd.get('ratio_to_null'), 2)}× | "
                 f"{_f(tt.get('L_um'))} | {_f(tt.get('ratio_to_null'), 2)}× |")
    L += [
        "",
        "**How to read this.** On the **fragmented** mask the drug L sits within "
        "~1.1–1.2× of the source null at every bandwidth and the two drugs are "
        "indistinguishable — there, the fitted decay *is* the KDE, exactly as the "
        "independent `results/delivery_metrics.md` analysis found across all five drugs. "
        "On the **repaired (connected)** mask the drugs separate from the null and from "
        "each other, but the null itself still moves with the bandwidth, so **part of "
        "the fitted decay on real tissue is our smoothing kernel in every case**. Quote "
        "the drug/null ratio, never the absolute fitted L, as the evidence.",
        "",
        "---",
        "",
        "## Sensitivity to the (arbitrary) vessel threshold — primary arm",
        "",
        "| percentile | threshold | vessel area | max distance (µm) | near-field L (µm) | "
        "near-field d½ (µm) | **empirical d½ (µm)** | *source-null d½ (µm)* | R² (log) |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for sv in res["threshold_sensitivity"]:
        L.append(
            f"| p{sv['percentile']:.0f} | {sv['threshold']:.3f} | "
            f"{sv['area_fraction']*100:.1f}% | {sv['max_distance_um']:.0f} | "
            f"{_f(sv['near_L_um'])} | {_f(sv['near_d_half_um'])} | "
            f"**{_f(sv['empirical_d_half_um'])}** | "
            f"*{_f(sv['source_empirical_d_half_um'])}* | {_f(sv['r2_log_near'], 4)} |")

    L += [
        "",
        "---",
        "",
        "## Honest caveats",
        "",
        "- **The nominally-converged grid breaks the transport domain.** At dx = 12 µm "
        "the cell-derived tissue mask is 26,058 disconnected fragments (largest = 0.7% "
        "of tissue), so lateral diffusion is switched off and *both* drugs collapse onto "
        "the vascular source field. `results/numerical_rigor.md` converged "
        "`kill_fraction` and `penetration_index` at dx ≲ 12 µm but did not test mask "
        "connectivity; this module does, and the primary arm therefore repairs the mask "
        "morphologically. **Any perivascular-gradient number quoted from a native "
        "dx ≤ 20 µm run is a source-geometry measurement, not a transport measurement.**",
        "- **The repair is a modelling choice.** Closing the mask adds interstitial "
        "pixels the rasteriser missed; it is defensible (real tissue is continuous) but "
        "it is not measured, and the closing radius is a free parameter. The dx = 30 µm "
        "native-mask arm is included precisely because it needs no morphology at all, "
        "and the two agree only qualitatively — hence the honest ranges quoted above "
        "rather than a single number.",
        "- **Part of the gradient is still vessel geometry, not drug transport.** The "
        f"source null alone gives d½ = {_f(sn['empirical_d_half_um'])} µm. The engine "
        "delivers drug proportionally to the *continuous* density everywhere, so tissue "
        "scored as 'far from a vessel' still has a local source, which flattens the "
        "predicted profile relative to a true line-source geometry.",
        "- **The vessel mask is not a measured vessel lumen.** `vascular_density` is a "
        "40 µm-bandwidth Gaussian KDE of the positions of cells *typed* as endothelial "
        "from the Xenium marker panel (`reconstruct.real_fields."
        "vascular_density_from_cells`); the real experiments segment CD31-stained vessel "
        "**lumens**. Our distances are distances from an endothelial-density ridge, and "
        "the 40 µm bandwidth is the same order as the gradient being measured — it puts "
        "a floor on how sharp the near-vessel profile can be.",
        "- **The threshold is a free parameter**, fixed *a priori* at the 95th percentile "
        "(≈5% vessel area fraction) rather than chosen to match the literature; the "
        "sensitivity table shows what moves with it.",
        "- **Doxorubicin agreement is PARTLY BY CONSTRUCTION.** The single α in "
        "`validation.calibration` was fit to the *same* Primeau anchor, through the 1-D "
        "analytic relation `L = √(D_eff/k_eff)`. The isolated single-vessel control "
        "therefore recovers ~45 µm essentially by construction and is a **consistency "
        "check that the calibration survives the full numerical solver on real geometry, "
        "not an independent validation**. What is genuinely non-trivial: the shape comes "
        "out exponential, and the real heterogeneous 2-D geometry did not move the "
        "number by an order of magnitude.",
        "- **The α convention matters and has changed.** This run used "
        f"**α = {res['alpha']:.1f}**, calibrated against the *e-folding* target "
        "`45/ln2 ≈ 64.9 µm` (`validation/calibration.py`, `DOX_TARGET_EFOLD_UM`). An "
        "earlier convention fitted α to the raw 45 µm half-distance as if it were an "
        "e-folding length, giving α = 309.6 and a ~1.4× steeper gradient. Every absolute "
        "penetration number in this report inherits that single convention choice; "
        "check `calibrate_alpha()` before comparing to older results files.",
        "- **Temozolomide is the independent output** — α was never fit to it — but "
        "there is **no published temozolomide perivascular-distance dataset**, so that "
        "arm is a *prediction*, checkable only for direction and order of magnitude. "
        "(The measured temozolomide anchor in `pk_anchors.py`, PMID 19861433, is a "
        "brain-ECF:plasma **partition ratio** — a different physical quantity that must "
        "not be compared to a penetration length.)",
        "- **2D, straight-line distances.** The section is a 2D plane through a 3D "
        "vascular network: a vessel just above or below the plane is invisible, so grid "
        "points scored 'far' may really be near an out-of-plane vessel. Distances are "
        "Euclidean, not geodesic around the irregular tissue mask.",
        "- **Not time-matched.** Primeau's images are ~10 min post-injection snapshots; "
        "the primary profile integrates a 6 h simulated exposure (AUC). The "
        "final-concentration readout is the closer snapshot analogue and agrees closely.",
        "- **Untreated section, no measured drug field.** The Xenium GBM sample is an "
        "untreated resection; the comparison is to published measurements in other "
        "(murine / xenograft) tumors, not to drug imaged on *this* tissue.",
        "",
    ]
    text = "\n".join(L)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


# --------------------------------------------------------------------------- #
def _print_summary(res: Dict) -> None:
    prim = res["primary"]
    print("\n" + "=" * 100)
    print("PERIVASCULAR DRUG GRADIENT on the REAL Xenium GBM section")
    print("=" * 100)
    hdr = (f"  {'arm':34s} {'conn':>6} {'signal':14s} {'emp d_half':>11} "
           f"{'L(exp+floor)':>13} {'R2':>7} {'corr(src)':>10}")
    print(hdr); print("-" * 100)
    for arm in res["arms"]:
        tag = arm["short"] + (" *PRIMARY*" if arm.get("primary") else "")
        conn = f"{arm['connectivity_used']['largest_fraction']*100:.0f}%"
        for r in arm["results"]:
            ef = r["auc"]["exp_floor"]
            print(f"  {tag:34s} {conn:>6} {r['drug']:14s} "
                  f"{_f(r['auc']['empirical_d_half_um']):>11} "
                  f"{_f(ef.get('L_um')):>13} {_f(ef.get('r2_linear_space'), 3):>7} "
                  f"{r['corr_with_source']:10.3f}")
            tag, conn = "", ""
        sn = arm["source_null"]
        print(f"  {'':34s} {'':>6} {'SOURCE NULL':14s} "
              f"{_f(sn['empirical_d_half_um']):>11} {_f(sn['exp_floor'].get('L_um')):>13}")
        print("-" * 100)
    print("ISOLATED SINGLE VESSEL (uniform medium, same solver -- the drug physics alone):")
    for iso in res["isolated"]:
        f = iso["fit"]
        print(f"  {iso['drug']:16s} L = {_f(f.get('L_um')):>7} um   d_half = "
              f"{_f(f.get('d_half_um')):>7} um   R2(log) = {_f(f.get('r2_log'), 5)}")
    print("=" * 100)
    print("vs MEASURED doxorubicin half-distance 40-50 um (Primeau 2005, PMID 16361566;")
    print("    replicated 25-75 um, Patel/Tannock 2013, PMID 23680920)")
    print(f"  primary arm (empirical d_half) : {prim['verdict_empirical']['verdict']}")
    print(f"  primary arm (L vs 58-72 um)    : {prim['verdict_L']['verdict']}")
    print(f"  isolated vessel (d_half)       : {res['verdict_isolated']['verdict']}")
    print(f"  isolated vessel (L vs 58-72)   : {res['verdict_isolated_L']['verdict']}")
    print("=" * 100)
    print("KDE-BANDWIDTH sweep (the decisive control -- fitted L vs our own smoothing):")
    for tag, key in (("repaired/connected", "bandwidth_primary"),
                     ("native/fragmented ", "bandwidth_native")):
        for b in res.get(key, []):
            parts = "  ".join(
                f"{n[:3]} L={b['drugs'][n]['L_um']:6.1f} ({b['drugs'][n]['ratio_to_null']:.2f}x null)"
                for n in b["drugs"])
            print(f"  {tag}  bw={b['bandwidth_um']:4.0f} um  null L="
                  f"{b['source_null_L_um']:6.1f}   {parts}")
    print("=" * 100)
    print("HONESTY:\n"
          "  * dx=12 um NATIVE mask is 99.3% fragmented -> no lateral diffusion -> both\n"
          "    drugs collapse onto the vascular SOURCE field. That arm is an artefact.\n"
          "  * doxorubicin is the drug alpha was CALIBRATED on (same Primeau anchor), so\n"
          "    its agreement is partly by construction -- a solver/geometry consistency\n"
          "    check, not an independent validation.\n"
          "  * the vessel mask is an endothelial-KDE threshold (40 um bandwidth), not a\n"
          "    segmented vessel lumen; 2D section; untreated resection, no measured drug.\n"
          "  * TEMOZOLOMIDE's profile is the free prediction -- but no published\n"
          "    perivascular-distance dataset exists for it, so it cannot be scored.\n"
          "  * the fitted decay on reconstructed tissue partly tracks OUR OWN 40 um KDE\n"
          "    bandwidth (see the sweep above) -- quote the drug/source-null RATIO, not\n"
          "    the absolute L. The analytic L = sqrt(D_eff/k), which does differ per\n"
          "    drug, is the legitimate model output to compare to literature.")
    print("=" * 100)


def main() -> Dict:
    res = run()
    png = make_figure(res)
    md = write_report(res)
    _print_summary(res)
    print(f"\nfigure written: {png}")
    print(f"report written: {md}")
    print(f"[runtime {res['runtime_s']:.1f} s]")
    return res


if __name__ == "__main__":
    main()
