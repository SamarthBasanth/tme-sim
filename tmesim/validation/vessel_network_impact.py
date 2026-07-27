"""DOES THE DISCRETE VESSEL NETWORK CHANGE THE ANSWER? — the A/B quantification.

``results/SESSION_AUDIT.md`` §3 names one structural confound as the highest-value
outstanding fix: the drug source is a **40 um-bandwidth Gaussian KDE of measured
endothelial-cell positions**, max-normalised, so it is nonzero over ~99.7% of the
tissue. Three independent analyses converged on its consequences:

* ``results/delivery_metrics.md`` — the fitted penetration length is 26-28 um for
  EVERY drug against a **source-null of 23.7 um**: the number measures the KDE
  bandwidth, not drug transport.
* ``results/dose_schedule.md`` — the exposure-vs-distance profile sits on a **flat
  floor at ~12% of peak**, so no tissue is truly drug-starved and *delivery is
  never the binding constraint*.
* ``results/niche_response.md`` — "perivascular tissue gets more drug" is circular:
  the region label and the source share the endothelial signal.

:mod:`reconstruct.vessel_network` builds the alternative — measured endothelial
cells clustered into DISCRETE vessel objects, painted as a source that is exactly
zero away from a vessel. **This module measures whether that changes anything.**
Writing the alternative is cheap; the load-bearing question is whether the model's
conclusions move, and by how much.

The experiment
--------------
One reconstruction (REAL 10x Xenium GBM, dx = 12 um, the production resolution),
three drugs spanning the delivery axis (temozolomide / cisplatin / doxorubicin),
three source fields:

===  =========================  ===================================================
arm  field                      what it isolates
===  =========================  ===================================================
A    KDE (the current default)  the status quo
B    discrete, unit amplitude   the geometry fix, at its own natural amplitude
C    discrete, total-matched    the DOSE CONTROL — same total drug entering the
                                tissue as A, delivered through B's geometry, so any
                                A-vs-C difference is geometry alone
===  =========================  ===================================================

C exists because B and A do not deliver the same amount of drug (the discrete
field covers 2.5% of tissue at 1.0; the KDE averages 0.084 over 100% of it). Any
coverage change from A to B is therefore *part geometry, part dose*, and reporting
it without C would overstate the structural effect. Both are reported.

What is measured
----------------
1. **The floor.** Exposure vs distance-from-vessel on a COMMON distance axis (the
   segmented discrete vessels — the only axis that means "distance from a vessel"
   for both arms), plus the ``C0*exp(-d/L) + floor`` fit's ``floor_fraction``, and
   a fit-free empirical floor (mean over bins beyond 150 um / perivascular bin).
2. **Decay-length separation.** Per drug, ``L`` from the same fit, each against
   **its own arm's source-null** — the identical fit applied to that arm's own
   ``vascular_density``. Under A the null is 23.7 um and every drug lands on it.
3. **Therapeutic coverage @ 1x IC50**, from :mod:`validation.delivery_metrics`,
   unmodified.
4. **Essentially-zero-drug fraction** of tumour cells, at two transparent
   thresholds (E < 0.01x IC50; E < 1% of that arm's own p99 tumour exposure).
5. **Is delivery binding?** Coverage in the far-from-vessel tumour-cell quartile
   vs the near quartile. Under a flat floor these are the same number.

What is NOT touched
-------------------
``core`` / ``physics`` / ``sim.engine`` / ``molecular`` are unmodified. The only
change outside this file and ``reconstruct/vessel_network.py`` is an additive,
default-unchanged ``mode=`` keyword on
:func:`reconstruct.real_fields.vascular_density_from_cells`.

Run:    ./.venv/bin/python -m tmesim.validation.vessel_network_impact
Writes: results/vessel_network_impact.md
        results/vessel_network_impact.png

Units: AUC uM*s, exposure/IC50 uM, length um, coverage/fractions dimensionless.
"""
from __future__ import annotations

import os
import pickle
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..core import (
    CellType,
    Microenvironment,
    OxygenConfig,
    PDEConfig,
    SimConfig,
    TissueFields,
)
from ..reconstruct.real_fields import add_oxygen
from ..reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from ..reconstruct.vessel_network import (
    LINK_DISTANCE_UM,
    PUBLISHED,
    VESSEL_RADIUS_UM,
    discrete_vascular_field,
    link_distance_sweep,
    segment_vessels,
    segmentation_report,
)
from ..sim.engine import simulate_response
from .calibration import calibrate_alpha
from .delivery_metrics import (
    _PDE,
    compile_real_drug,
    compute_all,
    penetration_length_fit,
)
from .vessel_profile import (
    BIN_UM,
    MIN_BIN_COUNT,
    NEAR_WINDOW_UM,
    VESSEL_PERCENTILE,
    distance_from_vessel_um,
    distance_profile,
    empirical_half_distance_um,
    fit_exponential_plus_floor,
    usable_bins,
    vessel_mask_from_density,
)

__all__ = [
    "DRUGS",
    "DX_UM",
    "ARMS",
    "build_arm_fields",
    "arm_microenvironment",
    "profile_stats",
    "zero_drug_fractions",
    "coverage_by_distance",
    "run",
    "make_figure",
    "write_report",
    "main",
]


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
#: Three drugs spanning the delivery axis of ``results/real_drug_panel.md``:
#: temozolomide (small, BBB-penetrant, weak), cisplatin (small, poorly penetrant,
#: mid potency), doxorubicin (large, excluded, very potent). If the discrete field
#: is going to separate anything, it has to separate these.
DRUGS: List[str] = ["temozolomide", "cisplatin", "doxorubicin"]

#: PRODUCTION resolution (``results/numerical_rigor.md``): the finest grid that
#: satisfies the advection CFL at dt = 10 s. Single dx on purpose — this study is
#: an A/B on the source field, not a grid study.
DX_UM = 12.0

#: The three source-field arms. ``key`` is used in every table and figure.
ARMS: List[Dict] = [
    {"key": "A_kde", "label": "A · KDE (current default)",
     "short": "A (KDE)", "mode": "kde"},
    {"key": "B_discrete", "label": "B · discrete vessels (unit amplitude)",
     "short": "B (discrete)", "mode": "discrete"},
    {"key": "C_discrete_matched",
     "label": "C · discrete vessels, total matched to A (dose control)",
     "short": "C (discrete, dose-matched)", "mode": "discrete_matched"},
]

#: Distance beyond which the profile is called "far field" for the fit-free floor.
FAR_FIELD_UM = 150.0

#: "Essentially zero drug" thresholds. Both are stated, neither is tuned.
ZERO_IC50_MULT = 0.01     # E < 1% of the drug's own IC50
ZERO_PEAK_FRAC = 0.01     # E < 1% of that arm's own p99 tumour-cell exposure

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "vessel_network_impact.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "vessel_network_impact.md")

#: Reused verbatim from ``validation.delivery_metrics`` so the two reports colour
#: the same drug the same way (validated categorical palette, dataviz skill
#: reference instance; worst adjacent CVD dE 24.2).
_DRUG_COLORS = {
    "temozolomide": "#2a78d6",   # blue
    "lomustine": "#1baf7a",      # aqua
    "cisplatin": "#eda100",      # yellow
    "methotrexate": "#008300",   # green
    "doxorubicin": "#4a3aa7",    # violet
}
_ARM_COLORS = {"A_kde": "#e34948", "B_discrete": "#2a78d6",
               "C_discrete_matched": "#1baf7a"}
_INK = "#0b0b0b"
_INK_2 = "#52514e"
_GRID = "#dcdbd6"

#: Published reference values quoted in ``results/delivery_metrics.md`` /
#: ``results/dose_schedule.md`` for the STATUS QUO, so the A arm can be checked
#: against what was reported before rather than silently re-derived.
PRIOR_KDE = {
    "L_drug_range_um": (26.5, 28.2),
    "L_source_null_um": 23.7,
    "floor_fraction": 0.12,
}


# --------------------------------------------------------------------------- #
# 1. Building the three arms
# --------------------------------------------------------------------------- #
def build_arm_fields(micro: Microenvironment,
                     link_distance_um: float = LINK_DISTANCE_UM,
                     radius_um: float = VESSEL_RADIUS_UM) -> Tuple[Dict, Dict]:
    """Segment the real endothelium and return ``(vessels, fields)``.

    ``fields`` maps each arm key to its (ny, nx) source field:

    * ``A_kde`` — the reconstruction's own ``micro.fields.vascular_density``,
      untouched (40 um Gaussian KDE of endothelial positions, max-normalised).
    * ``B_discrete`` — 1.0 inside ``radius_um`` of a segmented vessel, exactly 0.0
      outside. Relative, like A, but with a cleaner referent: 1.0 = "a perfused
      vessel wall leaking at the reference rate" rather than "the densest
      endothelial-KDE spot in this particular section".
    * ``C_discrete_matched`` — B rescaled so its tissue-integrated total equals
      A's. Holds the total drug entering the tissue fixed and changes only WHERE
      it enters.
    """
    domain = micro.domain
    kde = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    vessels = segment_vessels(micro.cells, domain, link_distance_um=link_distance_um)
    disc = discrete_vascular_field(domain, vessels, radius_um=radius_um,
                                   amplitude="unit")
    disc_matched = discrete_vascular_field(domain, vessels, radius_um=radius_um,
                                           amplitude="unit", match_total_to=kde)
    return vessels, {"A_kde": kde, "B_discrete": disc,
                     "C_discrete_matched": disc_matched}


def arm_microenvironment(micro: Microenvironment, field: np.ndarray,
                         tag: str) -> Microenvironment:
    """A copy of ``micro`` whose ``vascular_density`` is ``field``, with oxygen
    re-solved by the UNMODIFIED physics solver on that field.

    Cells are copied so the three arms cannot contaminate each other, and the
    caller's ``micro`` is left untouched. Note the honest consequence: the
    vascular field is the source for oxygen and (through ``physics.pressure``) for
    IFP as well as for drug, so an arm is a *self-consistent* alternative tissue,
    not a drug-only swap. ``sim.engine`` re-solves both internally regardless.
    """
    fields = TissueFields(
        ecm_density=np.array(micro.fields.ecm_density, dtype=np.float64, copy=True),
        vascular_density=np.asarray(field, dtype=np.float64),
        oxygen=None, ifp=None,
    )
    out = Microenvironment(domain=micro.domain, fields=fields,
                           cells=micro.cells.copy(),
                           patient_id=f"{micro.patient_id}::{tag}",
                           metadata=dict(micro.metadata))
    add_oxygen(out, OxygenConfig())
    return out


def field_stats(field: np.ndarray, tissue: np.ndarray, dx: float) -> Dict:
    """Descriptive statistics of a source field over the tissue mask."""
    f = np.asarray(field, dtype=np.float64)
    t = np.asarray(tissue, dtype=bool)
    v = f[t]
    return {
        "mean": float(v.mean()),
        "max": float(v.max()),
        "nonzero_area_fraction": float((v > 0).mean()),
        "zero_area_fraction": float((v == 0).mean()),
        "total_over_tissue": float(v.sum()),
        "integral_um2": float(v.sum()) * float(dx) ** 2,
        "n_tissue_px": int(t.sum()),
    }


# --------------------------------------------------------------------------- #
# 2. The floor: exposure vs distance from a segmented vessel
# --------------------------------------------------------------------------- #
def profile_stats(field: np.ndarray, distance_um: np.ndarray,
                  tissue: np.ndarray, far_um: float = FAR_FIELD_UM) -> Dict:
    """Bin ``field`` by distance-from-vessel and quantify the FLOOR two ways.

    * ``floor_fraction`` — from the three-parameter ``C0*exp(-d/L) + floor`` fit
      (:func:`vessel_profile.fit_exponential_plus_floor`), i.e.
      ``floor / (C0 + floor)``. This is the number ``results/dose_schedule.md``
      reports as ~0.12 for the KDE arm.
    * ``empirical_floor`` — FIT-FREE: the count-weighted mean of the profile over
      bins beyond ``far_um``, divided by the perivascular (first) bin. No
      functional form is assumed, so it cannot be an artefact of the fit.

    Also returns the raw profile (for plotting), the fitted ``L_um``, and the
    fit-free ``empirical_half_distance_um``.
    """
    prof = distance_profile(np.asarray(field, dtype=np.float64),
                            np.asarray(distance_um, dtype=np.float64),
                            np.asarray(tissue, dtype=bool), bin_um=BIN_UM)
    d, m, w = usable_bins(prof)
    fit = fit_exponential_plus_floor(d, m, w)

    peri = float(m[0]) if m.size else float("nan")
    far = d >= float(far_um)
    if far.any() and np.isfinite(peri) and peri > 0:
        emp_floor = float(np.average(m[far], weights=w[far]) / peri)
    else:
        emp_floor = float("nan")

    return {
        "distance_um": d, "mean": m, "count": w,
        "perivascular_mean": peri,
        "far_field_mean": (float(np.average(m[far], weights=w[far]))
                           if far.any() else float("nan")),
        "empirical_floor": emp_floor,
        "far_um": float(far_um),
        "L_um": float(fit.get("L_um", np.nan)),
        "floor_fraction": float(fit.get("floor_fraction", np.nan)),
        "C0": float(fit.get("C0", np.nan)),
        "r2_linear": float(fit.get("r2_linear_space", np.nan)),
        "d_half_empirical_um": empirical_half_distance_um(d, m),
        "fit_error": fit.get("error"),
        "n_bins": int(d.size),
    }


# --------------------------------------------------------------------------- #
# 3. Who gets essentially nothing, and is the deficit spatially organised?
# --------------------------------------------------------------------------- #
def zero_drug_fractions(exposure_uM: np.ndarray, ic50_uM: float) -> Dict:
    """Fraction of tumour cells receiving *essentially zero* drug.

    Two thresholds, both stated up front rather than chosen after seeing the
    answer, because "zero" is not a property a diffusion solve ever produces
    exactly:

    * ``frac_below_0.01x_ic50`` — DRUG-RELEVANT. ``E < 0.01 * IC50``: two orders of
      magnitude below a half-killing concentration is pharmacologically nothing.
    * ``frac_below_1pct_peak`` — DRUG-AGNOSTIC / shape only. ``E < 0.01 * p99(E)``
      within the same arm, so it measures the *dynamic range* of the delivery
      field and is invariant to how the field was normalised.

    ``frac_exactly_zero`` is reported too, and is expected to be ~0 in every arm:
    6 h of diffusion fills in any hole, so a source that is exactly zero does NOT
    imply an exposure that is exactly zero. That distinction is the point.
    """
    E = np.asarray(exposure_uM, dtype=np.float64).ravel()
    E = E[np.isfinite(E)]
    if E.size == 0:
        return {k: float("nan") for k in
                ("frac_below_0.01x_ic50", "frac_below_1pct_peak",
                 "frac_exactly_zero", "p99_uM", "min_uM")}
    p99 = float(np.percentile(E, 99.0))
    return {
        "frac_below_0.01x_ic50": float((E < ZERO_IC50_MULT * float(ic50_uM)).mean()),
        "frac_below_1pct_peak": float((E < ZERO_PEAK_FRAC * p99).mean()),
        "frac_exactly_zero": float((E <= 0.0).mean()),
        "p99_uM": p99,
        "min_uM": float(E.min()),
        "p10_over_p99": (float(np.percentile(E, 10.0) / p99) if p99 > 0
                         else float("nan")),
        "ic50_uM": float(ic50_uM),
        "zero_ic50_mult": float(ZERO_IC50_MULT),
        "zero_peak_frac": float(ZERO_PEAK_FRAC),
    }


def coverage_by_distance(exposure_uM: np.ndarray, tumor_to_vessel_um: np.ndarray,
                         ic50_uM: float) -> Dict:
    """Coverage @ 1x IC50 split by the tumour cell's own distance to a vessel.

    THE test of "is delivery the binding constraint?". Distance is measured from
    each tumour cell's measured coordinates to the nearest measured ENDOTHELIAL
    cell (``segment_vessels``' ``tumor_to_vessel_um``) — cell-to-cell, so it is
    exactly invariant to ``dx`` and independent of any field threshold.

    If the far quartile is covered as well as the near quartile, delivery is not
    what is limiting the drug: potency is. If the far quartile collapses, geometry
    has become a real constraint.
    """
    E = np.asarray(exposure_uM, dtype=np.float64).ravel()
    d = np.asarray(tumor_to_vessel_um, dtype=np.float64).ravel()
    if E.size != d.size or E.size == 0:
        return {"error": f"length mismatch E={E.size} d={d.size}"}
    thr = float(ic50_uM)
    q = np.percentile(d, [25.0, 50.0, 75.0])
    bins = [d <= q[0], (d > q[0]) & (d <= q[1]),
            (d > q[1]) & (d <= q[2]), d > q[2]]
    cov = [float((E[b] > thr).mean()) if b.any() else float("nan") for b in bins]
    med = [float(np.median(d[b])) if b.any() else float("nan") for b in bins]
    mean_E = [float(E[b].mean()) if b.any() else float("nan") for b in bins]
    near, far = cov[0], cov[3]
    return {
        "quartile_edges_um": [float(x) for x in q],
        "quartile_median_distance_um": med,
        "coverage_by_quartile": cov,
        "mean_E_by_quartile": mean_E,
        "coverage_near_q1": near,
        "coverage_far_q4": far,
        "near_minus_far": float(near - far),
        "far_over_near": (float(far / near) if near > 0 else float("nan")),
        "mean_E_far_over_near": (float(mean_E[3] / mean_E[0]) if mean_E[0] > 0
                                 else float("nan")),
    }


# --------------------------------------------------------------------------- #
# 4. The experiment
# --------------------------------------------------------------------------- #
def run(drugs: Sequence[str] = tuple(DRUGS),
        dx_um: float = DX_UM,
        link_distance_um: float = LINK_DISTANCE_UM,
        radius_um: float = VESSEL_RADIUS_UM,
        cache_path: Optional[str] = None) -> Dict:
    """Reconstruct once, build the three source fields, solve every (arm, drug).

    ``len(ARMS) * len(drugs)`` transport solves — 9 by default. No parameter
    sweeps: the linking distance and the vessel radius are each reported with
    their sensitivity (``link_distance_sweep``, and the field's area fraction
    against the published CD34 value) but the impact study is run at one setting.
    """
    if cache_path and os.path.exists(cache_path):
        with open(cache_path, "rb") as fh:
            return pickle.load(fh)

    t0 = time.time()
    alpha = calibrate_alpha()
    print(f"perivascular-sink calibration alpha = {alpha:.2f} (fit ONCE on doxorubicin)")

    try:
        micro = reconstruct_real("xenium", dx_um=float(dx_um))
    except Exception as exc:  # pragma: no cover - network/IO robustness
        print(f"reconstruction failed ({exc!r}); retrying once ...")
        time.sleep(3)
        micro = reconstruct_real("xenium", dx_um=float(dx_um))
    rep = reconstruction_report(micro)
    domain = micro.domain
    tissue = np.asarray(domain.tissue_mask, dtype=bool)
    print(f"REAL Xenium GBM  grid {rep['grid']}  dx {domain.dx:g} um  "
          f"{rep['n_cells']:,} cells (tumor {rep['n_tumor']:,})")

    # --- segmentation ------------------------------------------------------- #
    vessels, fields = build_arm_fields(micro, link_distance_um, radius_um)
    seg_rep = segmentation_report(vessels, domain, radius_um=radius_um)
    sweep = link_distance_sweep(micro.cells, domain)
    print(f"{vessels['n_endothelial_cells']:,} endothelial cells -> "
          f"{vessels['n_vessels']:,} vessels at {link_distance_um:g} um linking; "
          f"MVD {vessels['mvd_per_mm2']:.1f}/mm^2")

    fstats = {k: field_stats(v, tissue, domain.dx) for k, v in fields.items()}
    for k, s in fstats.items():
        print(f"  {k:20s} nonzero {s['nonzero_area_fraction']*100:6.2f}% of tissue, "
              f"mean {s['mean']:.4f}, max {s['max']:.3f}")

    # --- distance axes ------------------------------------------------------ #
    # COMMON axis: distance from a segmented DISCRETE vessel. This is the only
    # axis on which "distance from a vessel" means the same thing in both arms.
    vessel_mask = fields["B_discrete"] > 0
    dist_common = distance_from_vessel_um(vessel_mask, domain.dx)
    # LEGACY axis: the 95th-percentile threshold on the KDE that
    # delivery_metrics / vessel_profile use. Kept so the A arm can be checked
    # against the 26-28 um vs 23.7 um numbers already in results/.
    vm_legacy = vessel_mask_from_density(fields["A_kde"], tissue,
                                         percentile=VESSEL_PERCENTILE)
    dist_legacy = distance_from_vessel_um(vm_legacy["mask"], domain.dx)
    print(f"common distance axis: {int(vessel_mask[tissue].sum()):,} vessel px "
          f"({100*vessel_mask[tissue].mean():.2f}% of tissue); "
          f"legacy p{VESSEL_PERCENTILE:.0f} axis: {vm_legacy['n_vessel_px']:,} px")

    # --- source-null profiles (identical pipeline applied to the SOURCE) ----- #
    source_null = {}
    for k, f in fields.items():
        source_null[k] = {
            "common": profile_stats(f, dist_common, tissue),
            "legacy": profile_stats(f, dist_legacy, tissue),
        }

    # --- solve --------------------------------------------------------------- #
    per_arm: Dict[str, Dict] = {}
    for arm in ARMS:
        key = arm["key"]
        m_arm = arm_microenvironment(micro, fields[key], key)
        print(f"\n--- {arm['label']} ---")
        rows: Dict[str, Dict] = {}
        for name in drugs:
            dp, mol = compile_real_drug(name, alpha)
            t1 = time.time()
            res = simulate_response(
                m_arm, SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0))
            met = compute_all(m_arm, res, dp, distance_um=dist_common)
            auc = np.asarray(res.fields["auc"], dtype=np.float64)
            met["mw"] = float(mol.get("mw"))
            met["profile_common"] = profile_stats(auc, dist_common, tissue)
            met["profile_legacy"] = profile_stats(auc, dist_legacy, tissue)
            met["penetration_length_legacy"] = penetration_length_fit(
                auc, m_arm, distance_um=dist_legacy)
            met["zero_drug"] = zero_drug_fractions(met["_exposure_uM"],
                                                   met["ic50_uM"])
            met["coverage_by_distance"] = coverage_by_distance(
                met["_exposure_uM"], vessels["tumor_to_vessel_um"], met["ic50_uM"])
            met["mean_oxygen"] = float(res.extras.get("mean_oxygen", np.nan))
            met["hypoxic_fraction"] = float(res.extras.get("hypoxic_fraction", np.nan))
            rows[name] = met
            p = met["profile_common"]
            z = met["zero_drug"]
            print(f"  {name:14s} cov@IC50={met['therapeutic_coverage']['coverage_1x']:.4f}  "
                  f"L={p['L_um']:6.1f}um  floor={p['floor_fraction']:.3f}  "
                  f"empFloor={p['empirical_floor']:.4f}  "
                  f"zero<0.01xIC50={z['frac_below_0.01x_ic50']:.3f}  "
                  f"[{time.time()-t1:.0f}s]")
        per_arm[key] = rows

    out = {
        "alpha": float(alpha),
        "dx_um": float(dx_um),
        "drugs": list(drugs),
        "arms": ARMS,
        "link_distance_um": float(link_distance_um),
        "radius_um": float(radius_um),
        "reconstruction": rep,
        "vessels": {k: v for k, v in vessels.items()
                    if k not in ("labels", "xy", "centroids", "extent_um",
                                 "aspect_ratio", "intervessel_nn_um")},
        "vessel_xy": vessels["xy"],
        "vessel_labels": vessels["labels"],
        "intervessel_nn_um": vessels["intervessel_nn_um"],
        "tumor_to_vessel_um": vessels["tumor_to_vessel_um"],
        "segmentation_report": seg_rep,
        "link_sweep": sweep,
        "fields": fields,
        "field_stats": fstats,
        "source_null": source_null,
        "legacy_vessel_mask": {k: v for k, v in vm_legacy.items() if k != "mask"},
        "per_arm": per_arm,
        "grid": tuple(domain.shape),
        "tissue_area_mm2": float(tissue.sum()) * domain.dx ** 2 / 1e6,
        "runtime_s": time.time() - t0,
        "_domain": domain,
    }
    if cache_path:
        os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
        with open(cache_path, "wb") as fh:
            pickle.dump(out, fh)
    return out


# --------------------------------------------------------------------------- #
# 5. Derived verdicts
# --------------------------------------------------------------------------- #
def _spread(vals: Sequence[float]) -> float:
    """max/min of a set of positive values — the cross-drug separation factor."""
    v = np.asarray([x for x in vals if np.isfinite(x) and x > 0], dtype=np.float64)
    return float(v.max() / v.min()) if v.size >= 2 else float("nan")


def verdicts(res: Dict) -> Dict:
    """Turn the raw run into the four yes/no answers the study was run to get."""
    drugs = res["drugs"]
    A, B, C = (res["per_arm"][k] for k in ("A_kde", "B_discrete",
                                           "C_discrete_matched"))
    out: Dict = {}

    # --- 1. does the flat floor disappear? ---------------------------------- #
    out["floor"] = {}
    for key in ("A_kde", "B_discrete", "C_discrete_matched"):
        rows = res["per_arm"][key]
        out["floor"][key] = {
            "floor_fraction": [rows[d]["profile_common"]["floor_fraction"] for d in drugs],
            "empirical_floor": [rows[d]["profile_common"]["empirical_floor"] for d in drugs],
            "source_floor_fraction":
                res["source_null"][key]["common"]["floor_fraction"],
            "source_empirical_floor":
                res["source_null"][key]["common"]["empirical_floor"],
        }
    a_emp = np.nanmean(out["floor"]["A_kde"]["empirical_floor"])
    b_emp = np.nanmean(out["floor"]["B_discrete"]["empirical_floor"])
    out["floor"]["A_mean_empirical"] = float(a_emp)
    out["floor"]["B_mean_empirical"] = float(b_emp)
    out["floor"]["fold_reduction"] = float(a_emp / b_emp) if b_emp > 0 else float("inf")
    out["floor"]["disappears"] = bool(b_emp < 0.5 * a_emp)

    # --- 2. do drugs separate on decay length? ------------------------------ #
    out["decay"] = {}
    for key in ("A_kde", "B_discrete", "C_discrete_matched"):
        rows = res["per_arm"][key]
        L = [rows[d]["profile_common"]["L_um"] for d in drugs]
        Lnull = res["source_null"][key]["common"]["L_um"]
        out["decay"][key] = {
            "L_um": L,
            "L_source_null_um": Lnull,
            "spread": _spread(L),
            "L_over_null": [(float(x / Lnull) if np.isfinite(Lnull) and Lnull > 0
                             else float("nan")) for x in L],
            "d_half_empirical_um": [rows[d]["profile_common"]["d_half_empirical_um"]
                                    for d in drugs],
        }
    # legacy axis, arm A only — the reproduction of results/delivery_metrics.md
    out["decay"]["A_legacy_axis"] = {
        "L_um": [A[d]["penetration_length_legacy"]["drug"]["L_um"] for d in drugs],
        "L_source_null_um":
            A[drugs[0]]["penetration_length_legacy"]["source_null"]["L_um"],
        "prior_reported": PRIOR_KDE,
    }
    sa = out["decay"]["A_kde"]["spread"]
    sb = out["decay"]["B_discrete"]["spread"]
    out["decay"]["separate"] = bool(np.isfinite(sb) and np.isfinite(sa) and sb > 1.5 * sa)

    # --- 3. coverage ---------------------------------------------------------- #
    out["coverage"] = {
        key: {d: res["per_arm"][key][d]["therapeutic_coverage"]["coverage_1x"]
              for d in drugs}
        for key in ("A_kde", "B_discrete", "C_discrete_matched")
    }
    out["coverage"]["delta_B_minus_A"] = {
        d: float(B[d]["therapeutic_coverage"]["coverage_1x"]
                 - A[d]["therapeutic_coverage"]["coverage_1x"]) for d in drugs}
    out["coverage"]["delta_C_minus_A"] = {
        d: float(C[d]["therapeutic_coverage"]["coverage_1x"]
                 - A[d]["therapeutic_coverage"]["coverage_1x"]) for d in drugs}

    # --- 4. is delivery the binding constraint? ------------------------------ #
    out["binding"] = {}
    for key in ("A_kde", "B_discrete", "C_discrete_matched"):
        rows = res["per_arm"][key]
        out["binding"][key] = {
            "coverage_near_q1": [rows[d]["coverage_by_distance"]["coverage_near_q1"]
                                 for d in drugs],
            "coverage_far_q4": [rows[d]["coverage_by_distance"]["coverage_far_q4"]
                                for d in drugs],
            "near_minus_far": [rows[d]["coverage_by_distance"]["near_minus_far"]
                               for d in drugs],
            "mean_E_far_over_near": [
                rows[d]["coverage_by_distance"]["mean_E_far_over_near"] for d in drugs],
            "frac_below_0.01x_ic50": [rows[d]["zero_drug"]["frac_below_0.01x_ic50"]
                                      for d in drugs],
            "frac_below_1pct_peak": [rows[d]["zero_drug"]["frac_below_1pct_peak"]
                                     for d in drugs],
            "p10_over_p99": [rows[d]["zero_drug"]["p10_over_p99"] for d in drugs],
        }
    return out


# --------------------------------------------------------------------------- #
# 6. Figure
# --------------------------------------------------------------------------- #
def _crop_window(res: Dict, size_um: float = 700.0) -> Tuple[float, float, float, float]:
    """A representative square window for the side-by-side source zoom.

    Chosen deterministically as the ``size_um`` box (on a grid aligned to
    ``size_um/2``) containing the MEDIAN number of endothelial cells among boxes
    that contain any — not the densest box, which would flatter both fields.
    """
    xy = np.asarray(res["vessel_xy"], dtype=np.float64)
    if xy.shape[0] == 0:
        return (0.0, size_um, 0.0, size_um)
    step = size_um / 2.0
    gx = np.floor(xy[:, 0] / step).astype(int)
    gy = np.floor(xy[:, 1] / step).astype(int)
    keys, counts = np.unique(np.column_stack([gx, gy]), axis=0, return_counts=True)
    order = np.argsort(counts)
    # Median-density box among the top half (so it is populated but not a hotspot).
    pick = keys[order[int(0.80 * order.size)]]
    x0, y0 = float(pick[0]) * step, float(pick[1]) * step
    return (x0, x0 + size_um, y0, y0 + size_um)


def _imcrop(field: np.ndarray, domain, win) -> Tuple[np.ndarray, List[float]]:
    x0, x1, y0, y1 = win
    dx = float(domain.dx)
    ny, nx = field.shape
    c0, c1 = int(max(x0 / dx, 0)), int(min(np.ceil(x1 / dx), nx))
    r0, r1 = int(max(y0 / dx, 0)), int(min(np.ceil(y1 / dx), ny))
    return field[r0:r1, c0:c1], [c0 * dx, c1 * dx, r1 * dx, r0 * dx]


def make_figure(res: Dict, path: str = _OUT_PNG) -> str:
    """Nine panels: the two source fields side by side, the source-level floor,
    the exposure-vs-distance profiles under each field, the floor per drug, the
    per-drug decay separation, the coverage change, and coverage vs distance."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    drugs = res["drugs"]
    dom = res["_domain"]
    vd = res.get("_verdicts") or verdicts(res)
    plt.rcParams.update({
        "font.size": 9.5, "axes.edgecolor": _GRID, "axes.labelcolor": _INK_2,
        "xtick.color": _INK_2, "ytick.color": _INK_2, "text.color": _INK,
        "axes.titlecolor": _INK, "figure.facecolor": "white",
        "axes.facecolor": "white",
    })
    fig = plt.figure(figsize=(18.5, 15.5))
    gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.26,
                  left=0.055, right=0.985, top=0.905, bottom=0.055)

    def _clean_ax(ax, ygrid=True):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if ygrid:
            ax.grid(axis="y", color=_GRID, lw=0.8, zorder=0)
            ax.set_axisbelow(True)

    # ---- (a,b) the two source fields, same window --------------------------- #
    win = _crop_window(res)
    xy = np.asarray(res["vessel_xy"], dtype=np.float64)
    sel = ((xy[:, 0] >= win[0]) & (xy[:, 0] <= win[1]) &
           (xy[:, 1] >= win[2]) & (xy[:, 1] <= win[3]))
    for i, (key, title) in enumerate((
            ("A_kde", "a · A — KDE source (40 µm bandwidth, max-normalised)"),
            ("B_discrete", "b · B — discrete vessel source (r = "
                           f"{res['radius_um']:g} µm, exactly 0 elsewhere)"))):
        ax = fig.add_subplot(gs[0, i])
        crop, ext = _imcrop(res["fields"][key], dom, win)
        im = ax.imshow(crop, cmap="inferno", extent=ext, vmin=0,
                       vmax=float(np.percentile(res["fields"]["A_kde"], 99.9)),
                       interpolation="nearest")
        ax.scatter(xy[sel, 0], xy[sel, 1], s=7, facecolors="none",
                   edgecolors="#7fd4ff", linewidths=0.7, zorder=3)
        fs = res["field_stats"][key]
        ax.set_title(f"{title}\n{fs['nonzero_area_fraction']*100:.1f}% of tissue "
                     f"nonzero · mean {fs['mean']:.4f} · cyan rings = measured "
                     f"endothelial cells", fontsize=9.5, loc="left")
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("y (µm)")
        fig.colorbar(im, ax=ax, shrink=0.82)

    # ---- (c) the floor is already in the SOURCE ----------------------------- #
    ax = fig.add_subplot(gs[0, 2])
    for key, lab in (("A_kde", "A · KDE source"),
                     ("B_discrete", "B · discrete source")):
        p = res["source_null"][key]["common"]
        y = p["mean"] / p["perivascular_mean"]
        ax.plot(p["distance_um"], np.maximum(y, 1e-6), lw=2.2,
                color=_ARM_COLORS[key], label=lab, zorder=3)
    ax.axhline(PRIOR_KDE["floor_fraction"], color=_INK_2, lw=1.0, ls=":", zorder=2)
    ax.text(ax.get_xlim()[1], PRIOR_KDE["floor_fraction"], "  12% floor reported\n"
            "  in dose_schedule.md", va="center", fontsize=8, color=_INK_2)
    ax.set_yscale("log")
    ax.set_ylim(1e-6, 2)
    ax.set_xlabel("distance from a segmented vessel (µm)")
    ax.set_ylabel("source / perivascular source")
    ax.legend(frameon=False, fontsize=9)
    _clean_ax(ax)
    ax.set_title("c · THE SOURCE ITSELF — the KDE never reaches zero;\n"
                 "the discrete field does, by construction",
                 fontsize=10, loc="left")

    # ---- (d,e) exposure vs distance under each field ------------------------ #
    for i, key in enumerate(("A_kde", "B_discrete")):
        ax = fig.add_subplot(gs[1, i])
        for d in drugs:
            p = res["per_arm"][key][d]["profile_common"]
            y = p["mean"] / p["perivascular_mean"]
            ax.plot(p["distance_um"], np.maximum(y, 1e-6), lw=2.2,
                    color=_DRUG_COLORS.get(d, _INK), label=d, zorder=3)
        p0 = res["source_null"][key]["common"]
        ax.plot(p0["distance_um"],
                np.maximum(p0["mean"] / p0["perivascular_mean"], 1e-6),
                lw=1.4, ls="--", color=_INK_2, label="source null", zorder=2)
        floors = vd["floor"][key]["empirical_floor"]
        ax.axhline(float(np.nanmean(floors)), color=_ARM_COLORS[key], lw=1.2,
                   ls=":", zorder=2)
        ax.set_yscale("log")
        ax.set_ylim(1e-4, 2)
        ax.set_xlabel("distance from a segmented vessel (µm)")
        ax.set_ylabel("exposure / perivascular exposure")
        ax.legend(frameon=False, fontsize=8.5, ncol=2)
        _clean_ax(ax)
        lbl = "d" if i == 0 else "e"
        ax.set_title(f"{lbl} · exposure profile under "
                     f"{'A — KDE' if i == 0 else 'B — DISCRETE'}\n"
                     f"mean far-field floor (>{FAR_FIELD_UM:g} µm) = "
                     f"{np.nanmean(floors)*100:.2f}% of perivascular",
                     fontsize=10, loc="left")

    # ---- (f) the floor per drug, all three arms ----------------------------- #
    ax = fig.add_subplot(gs[1, 2])
    x = np.arange(len(drugs))
    w = 0.26
    for j, key in enumerate(("A_kde", "B_discrete", "C_discrete_matched")):
        vals = vd["floor"][key]["empirical_floor"]
        ax.bar(x + (j - 1) * w, vals, w, color=_ARM_COLORS[key], zorder=3,
               label=next(a["short"] for a in ARMS if a["key"] == key))
        for xi, v in zip(x + (j - 1) * w, vals):
            ax.text(xi, v * 1.15, f"{v*100:.2f}%", ha="center", va="bottom",
                    fontsize=7.5, color=_INK)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(drugs, fontsize=9)
    ax.set_ylabel(f"far-field (> {FAR_FIELD_UM:g} µm) / perivascular")
    ax.legend(frameon=False, fontsize=8.5)
    _clean_ax(ax)
    ax.set_title("f · THE FLAT FLOOR, fit-free\n"
                 f"A→B changes it by {vd['floor']['fold_reduction']:.2f}×",
                 fontsize=10, loc="left")

    # ---- (g) decay length per drug, each vs its own source null ------------- #
    ax = fig.add_subplot(gs[2, 0])
    for j, key in enumerate(("A_kde", "B_discrete")):
        L = vd["decay"][key]["L_um"]
        ax.bar(x + (j - 0.5) * 0.34, L, 0.34, color=_ARM_COLORS[key], zorder=3,
               label=next(a["short"] for a in ARMS if a["key"] == key))
        for xi, v in zip(x + (j - 0.5) * 0.34, L):
            ax.text(xi, v, f"{v:.0f}", ha="center", va="bottom", fontsize=8,
                    color=_INK)
        Ln = vd["decay"][key]["L_source_null_um"]
        if np.isfinite(Ln):
            ax.axhline(Ln, color=_ARM_COLORS[key], lw=1.3, ls="--", zorder=2)
            ax.text(len(drugs) - 0.45, Ln, f" {key.split('_')[0]} source null "
                    f"{Ln:.1f} µm", fontsize=8, color=_ARM_COLORS[key], va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(drugs, fontsize=9)
    ax.set_ylabel("fitted decay length L (µm)")
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    _clean_ax(ax)
    sa, sb = vd["decay"]["A_kde"]["spread"], vd["decay"]["B_discrete"]["spread"]
    ax.set_title("g · DO DRUGS SEPARATE?  cross-drug spread (max/min)\n"
                 f"A {sa:.2f}×   →   B {sb:.2f}×", fontsize=10, loc="left")

    # ---- (h) therapeutic coverage ------------------------------------------- #
    ax = fig.add_subplot(gs[2, 1])
    for j, key in enumerate(("A_kde", "B_discrete", "C_discrete_matched")):
        vals = [vd["coverage"][key][d] for d in drugs]
        ax.bar(x + (j - 1) * w, vals, w, color=_ARM_COLORS[key], zorder=3,
               label=next(a["short"] for a in ARMS if a["key"] == key))
        for xi, v in zip(x + (j - 1) * w, vals):
            ax.text(xi, v + 0.02, f"{v:.3f}", ha="center", va="bottom",
                    fontsize=7.5, color=_INK)
    ax.set_xticks(x)
    ax.set_xticklabels(drugs, fontsize=9)
    ax.set_ylabel("fraction of tumour cells with E > IC50")
    ax.set_ylim(0, 1.14)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    _clean_ax(ax)
    ax.set_title("h · THERAPEUTIC COVERAGE @ 1× IC50\n"
                 "C isolates geometry: same total drug as A", fontsize=10, loc="left")

    # ---- (i) is delivery the binding constraint? ---------------------------- #
    ax = fig.add_subplot(gs[2, 2])
    for d in drugs:
        for key, ls in (("A_kde", "-"), ("B_discrete", "--")):
            cbd = res["per_arm"][key][d]["coverage_by_distance"]
            ax.plot(cbd["quartile_median_distance_um"], cbd["coverage_by_quartile"],
                    ls=ls, marker="o", ms=5, lw=2.0,
                    color=_DRUG_COLORS.get(d, _INK), zorder=3,
                    label=f"{d} · {'A' if key == 'A_kde' else 'B'}")
    ax.set_xlabel("median tumour-cell→vessel distance of the quartile (µm)")
    ax.set_ylabel("coverage @ 1× IC50 within quartile")
    ax.set_ylim(-0.03, 1.08)
    ax.legend(frameon=False, fontsize=7.5, ncol=2)
    _clean_ax(ax)
    ax.set_title("i · IS DELIVERY BINDING?  coverage by distance quartile\n"
                 "flat = delivery is not the constraint (solid A, dashed B)",
                 fontsize=10, loc="left")

    v = res["vessels"]
    fig.suptitle(
        "Does a DISCRETE vessel network change the model's conclusions?\n"
        f"REAL 10x Xenium GBM · dx = {res['dx_um']:g} µm · "
        f"{v['n_endothelial_cells']:,} measured endothelial cells → "
        f"{v['n_vessels']:,} discrete vessels at {res['link_distance_um']:g} µm "
        f"linking (MVD {v['mvd_per_mm2']:.0f}/mm², published GBM 65 ± 30, "
        f"PMID 27583442) · {len(ARMS)} source fields × {len(drugs)} drugs",
        fontsize=13.5, x=0.055, ha="left", y=0.975)

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# 7. Report
# --------------------------------------------------------------------------- #
def _f(x, nd=3, na="—"):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return na
    if not np.isfinite(v):
        return na
    return f"{v:,.{nd}f}"


def _pct(x, nd=2, na="—"):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return na
    return f"{v*100:.{nd}f}%" if np.isfinite(v) else na


def write_report(res: Dict, path: str = _OUT_MD) -> str:
    """Write ``results/vessel_network_impact.md``."""
    vd = res.get("_verdicts") or verdicts(res)
    drugs = res["drugs"]
    v = res["vessels"]
    seg = res["segmentation_report"]
    fs = res["field_stats"]
    L = []
    A = res["per_arm"]["A_kde"]
    B = res["per_arm"]["B_discrete"]
    C = res["per_arm"]["C_discrete_matched"]

    floor_A = vd["floor"]["A_mean_empirical"]
    floor_B = vd["floor"]["B_mean_empirical"]
    spread_A = vd["decay"]["A_kde"]["spread"]
    spread_B = vd["decay"]["B_discrete"]["spread"]
    covA = [vd["coverage"]["A_kde"][d] for d in drugs]
    covB = [vd["coverage"]["B_discrete"][d] for d in drugs]
    covC = [vd["coverage"]["C_discrete_matched"][d] for d in drugs]

    L.append("# Does a DISCRETE vessel network change the model's conclusions?")
    L.append("")
    L.append("**The single highest-value structural fix in `results/SESSION_AUDIT.md` §3, "
             "built and then measured.** The drug source in this model is a 40 µm-bandwidth "
             "Gaussian KDE of measured endothelial-cell positions. This report replaces it "
             "with discrete segmented vessels — a source that is **exactly zero** away from "
             "a vessel — and quantifies what moves.")
    L.append("")
    L.append(f"Generated by `tmesim.validation.vessel_network_impact` in "
             f"{res['runtime_s']:.0f} s on the REAL 10x Xenium GBM section "
             f"(`reconstruct_real('xenium', dx_um={res['dx_um']:g})`, grid "
             f"{res['grid'][0]}×{res['grid'][1]}, {res['reconstruction']['n_cells']:,} "
             f"measured cells, {res['reconstruction']['n_tumor']:,} tumour). "
             f"{len(ARMS)} source fields × {len(drugs)} drugs = "
             f"{len(ARMS)*len(drugs)} transport solves. "
             f"Figure: `vessel_network_impact.png`. "
             f"Modules: `reconstruct/vessel_network.py`, "
             f"`validation/vessel_network_impact.py`.")
    L.append("")
    L.append("## Verdict, up front")
    L.append("")
    return "\n".join(L)  # placeholder replaced below
