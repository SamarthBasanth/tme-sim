"""WHY IS THE RECONSTRUCTED-SECTION SPATIAL READOUT NOISY? — a stage-by-stage
decomposition of the arm-B failure in ``results/kyle2004_penetration.md`` §8.

The finding this file investigates
----------------------------------
``results/kyle2004_penetration.md`` C2 fails on the reconstructed section:
doxorubicin spans 39.4-63.1 um and daunorubicin 32.7-53.9 um across a registered
nuisance grid (tissue arm x field x distance metric x bin width), so the compound
difference sits inside its own methodological noise. The same file's arm A -- the
clean single-vessel slab -- is DISJOINT. That localises the problem to the
spatial reduction, not the transport physics.

The reduction has four stages, and a single moved number cannot say which one
moved:

    1. SOLVE      the PDE field on a tissue mask  (grid dx, mask repair)
    2. DISTANCE   assign every tissue pixel a distance-to-vessel  (metric)
    3. BIN        reduce (d_i, C_i) pixels to a binned profile     (bin width,
                                                                   truncation)
    4. ESTIMATE   reduce the profile to one number  (empirical d-half / exp+floor)

This module attributes the observed spread to each stage separately, and it does
so by CACHING the solved fields once and then replaying every downstream choice
on the identical fields. Stages 2-4 therefore cost no extra solve, and the
attribution is exact rather than inferred from marginal means.

NOTHING IS TUNED. No shipped default, constant or parameter is changed; the
alpha calibration is asserted live at 22.223477711172563 in every stage. The
estimator variants defined here are ADDITIONS that live next to the shipped one,
and the shipped one is called unmodified so that every published number
reproduces bit-for-bit inside this run.

Stages
------
``--stage bias``       closed-form + synthetic demonstration of the estimator's
                       bin-width bias. No PDE solve; seconds.
``--stage fields``     6 real-section solves (3 tissue arms x 2 compounds),
                       cached to disk. ~15 min.
``--stage decompose``  replay every reduction choice on the cached fields and
                       score every estimator. Seconds once ``fields`` has run.
``--stage converge``   fine-to-coarse bin sweep: which labelling agrees with the
                       ``bin_um -> 0`` limit? Uses the cache.
``--stage repair``     split the mask-repair effect into "the field changed" vs
                       "the pixel set changed". Uses the cache.
``--stage all``        bias -> fields -> decompose -> converge -> repair.

Writes ``results/spatial_noise_decomposition.json``. The field cache is large and
goes to a scratch directory, not to ``results/``.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.core import Domain, Microenvironment, PDEConfig, SimConfig  # noqa: E402
from tmesim.molecular import descriptors as _descriptors  # noqa: E402
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams  # noqa: E402
from tmesim.molecular.compile_drug import build_profile, molecule_to_drugparams  # noqa: E402
from tmesim.reconstruct.real_fields import build_real_fields  # noqa: E402
from tmesim.sim.engine import simulate_response  # noqa: E402
from tmesim.validation import vessel_profile as VP  # noqa: E402
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams  # noqa: E402

# --------------------------------------------------------------------------- #
# 0. Frozen constants — asserted, never redefined
# --------------------------------------------------------------------------- #
ALPHA_EXPECTED = 22.223477711172563

DOX = "doxorubicin"
DAU = "daunorubicin"
PAIR: Tuple[str, str] = (DOX, DAU)

#: Reused VERBATIM from experiments/isomer_penetration_map.py so this file cannot
#: introduce a second, silently diverging definition of the same molecule.
SMILES = {
    DOX: ("C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5"
          "cccc(c5C4=O)OC)O)(C(=O)CO)O)N)O"),
    DAU: ("C[C@H]1[C@H]([C@H](C[C@@H](O1)O[C@H]2C[C@@](Cc3c2c(c4c(c3O)C(=O)c5"
          "cccc(c5C4=O)OC)O)(C(=O)C)O)N)O"),
}
SHARED_IC50_UM = 1.0

#: The registered arm-B nuisance grid (kyle2004_penetration.md §8), reproduced
#: exactly, plus one axis that file left IMPLICIT: the profile's distance
#: truncation. ``distance_profile`` truncates at the 99.5th percentile of the
#: tissue distance distribution by default and never sweeps it.
NUISANCE_ARMS: Sequence[str] = ("dx12_repaired", "dx12_native", "dx30_native")
NUISANCE_FIELDS: Sequence[str] = ("auc", "conc")
NUISANCE_METRICS: Sequence[str] = ("primary", "legacy")
NUISANCE_BINS_UM: Sequence[float] = (5.0, 10.0, 20.0)
NUISANCE_TRUNC_PCT: Sequence[float] = (99.0, 99.5, 100.0)

_RESULTS = os.path.join(_ROOT, "results")
_OUT_JSON = os.path.join(_RESULTS, "spatial_noise_decomposition.json")

_SCRATCH = os.environ.get(
    "TMESIM_SCRATCH",
    "/private/tmp/claude-501/-Users-samarthbasanth/"
    "05d22d67-3289-4ed1-ba30-fe372c1e0d5e/scratchpad",
)
_CACHE_DIR = os.path.join(_SCRATCH, "spatial_noise_fields")

_PDE = dict(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)


def _assert_alpha() -> float:
    a = calibrate_alpha()
    if abs(a - ALPHA_EXPECTED) > 1e-12:
        raise AssertionError(
            f"calibrate_alpha() = {a!r}, expected {ALPHA_EXPECTED!r}. "
            "The calibration constant moved; every number in this file is void."
        )
    return float(a)


def compile_pair(alpha: float) -> Dict:
    """The identical compilation path ``isomer_penetration_map.compile_all`` uses."""
    out = {}
    for name in PAIR:
        mol = _descriptors.from_smiles(name, SMILES[name])
        prof = build_profile(mol, ic50=SHARED_IC50_UM, two_route=True)
        dp = molecule_to_drugparams(prof, ic50=SHARED_IC50_UM)
        dp = calibrated_drugparams(dp, alpha)
        out[name] = apply_barrier_to_drugparams(dp, mol)
    return out


# --------------------------------------------------------------------------- #
# 1. STAGE 4 — the estimators
# --------------------------------------------------------------------------- #
# The shipped estimator is VP.empirical_half_distance_um, called unmodified.
# Everything below is an addition placed beside it.

#: The anchor-corrected half-distance lives in ``validation/vessel_profile.py``
#: next to the shipped estimator it corrects, so the two cannot drift apart and
#: this file cannot carry a second, silently diverging definition of it. Its
#: derivation is in that function's docstring; ``--stage bias`` verifies the
#: closed form numerically.
anchor_corrected_half_distance_um = VP.anchor_corrected_half_distance_um


def exposure_weighted_depth(
    distance_um: np.ndarray, value: np.ndarray, weights: Optional[np.ndarray] = None
) -> Dict:
    """BINNING-FREE depth statistics computed straight from the pixel cloud.

    No bins, no bin width, no truncation, no anchor bin, no fit, no refusal
    branch. Every one of the reduction choices that the decomposition attributes
    spread to is structurally absent, in exactly the sense that
    ``results/delivery_spread_noise.md`` §1 documents for ``P90/P10``.

    Returns, for the distance distribution re-weighted by delivered signal:

    ``depth_mean_um``
        ``sum(C_i d_i) / sum(C_i)`` — the centroid of delivered exposure in
        distance space. For a planar vessel in a uniform 1-D medium with
        ``C = C0 exp(-d/L)`` this is exactly ``L``.
    ``depth_median_um``
        the weighted median. For the same geometry it is exactly ``L*ln2``, i.e.
        the half-distance itself.
    ``tissue_mean_um`` / ``tissue_median_um``
        the same two statistics with ``C`` replaced by 1 — the tissue's own
        distance geometry, which is what a drug that reached everywhere equally
        would score.
    ``depth_ratio`` / ``depth_median_ratio``
        the drug statistic over the tissue statistic. Dimensionless, and the
        tissue geometry — which is what the grid and the mask change — divides
        out. Values below 1 mean the delivered exposure sits closer to vessels
        than the tissue does.

    All of these are homogeneous of degree 0 in ``C``, so the AUC-vs-final-
    concentration axis is degenerate for them by the same algebra
    ``delivery_spread_noise.md`` §1 measures for ``P90/P10`` (``exposure = auc/T``
    with ``T`` one scalar per run).
    """
    d = np.asarray(distance_um, dtype=np.float64).ravel()
    c = np.asarray(value, dtype=np.float64).ravel()
    ok = np.isfinite(d) & np.isfinite(c) & (c > 0.0)
    d, c = d[ok], c[ok]
    if weights is not None:
        w = np.asarray(weights, dtype=np.float64).ravel()[ok]
    else:
        w = np.ones_like(d)
    if d.size == 0:
        return {}
    cw = c * w
    order = np.argsort(d)
    d_s, cw_s, w_s = d[order], cw[order], w[order]

    def _wmedian(x: np.ndarray, ww: np.ndarray) -> float:
        cs = np.cumsum(ww)
        tot = cs[-1]
        return float(np.interp(0.5 * tot, cs, x))

    depth_mean = float(np.sum(cw_s * d_s) / np.sum(cw_s))
    depth_med = _wmedian(d_s, cw_s)
    tiss_mean = float(np.sum(w_s * d_s) / np.sum(w_s))
    tiss_med = _wmedian(d_s, w_s)
    return {
        "depth_mean_um": depth_mean,
        "depth_median_um": depth_med,
        "tissue_mean_um": tiss_mean,
        "tissue_median_um": tiss_med,
        "depth_ratio": depth_mean / tiss_mean if tiss_mean else float("nan"),
        "depth_median_ratio": depth_med / tiss_med if tiss_med else float("nan"),
        "n_px": int(d.size),
    }


def pixel_log_slope_L_um(distance_um: np.ndarray, value: np.ndarray,
                         max_um: Optional[float] = None) -> Optional[float]:
    """``L`` from an UNBINNED least-squares fit of ``log C_i`` on ``d_i``.

    The control that separates "binning is the problem" from "the exponential
    shape is the problem": it is the same log-linear estimator the binned
    ``fit_exponential`` uses, with the binning stage deleted and nothing else
    changed. If the binned fit is unstable and this one is not, the instability
    is the binning; if both are unstable, it is the model shape.
    """
    d = np.asarray(distance_um, dtype=np.float64).ravel()
    c = np.asarray(value, dtype=np.float64).ravel()
    ok = np.isfinite(d) & np.isfinite(c) & (c > 0.0)
    if max_um is not None:
        ok &= d <= float(max_um)
    d, c = d[ok], c[ok]
    if d.size < 10:
        return None
    slope, _ = np.polyfit(d, np.log(c), 1)
    return float(-1.0 / slope) if slope < 0 else None


def p90_p10(value: np.ndarray) -> float:
    """The delivery-spread statistic of ``results/delivery_spread_noise.md``,
    computed over the same pixel set — a positive control for "a readout of this
    field that involves no distance reduction at all"."""
    v = np.asarray(value, dtype=np.float64).ravel()
    v = v[np.isfinite(v) & (v > 0)]
    if v.size == 0:
        return float("nan")
    lo = float(np.percentile(v, 10.0))
    return float(np.percentile(v, 90.0)) / lo if lo > 0 else float("nan")


# --------------------------------------------------------------------------- #
# 2. STAGE 3 — binning, with the bin's own mean distance recorded
# --------------------------------------------------------------------------- #
def binned_profile(distance_px: np.ndarray, value_px: np.ndarray,
                   bin_um: float, max_um: float) -> Dict:
    """Bin the pixel cloud exactly as ``VP.distance_profile`` does, and ALSO
    record each bin's count-weighted MEAN distance.

    ``VP.distance_profile`` labels every bin by its geometric centre. The pixels
    inside a bin are not uniformly distributed in distance (near a vessel the
    density of tissue at distance ``d`` rises with ``d``), so the centre is not
    where the bin's mean value actually sits. Recording ``distance_mean`` costs
    one extra reduction and lets the anchor correction use the true anchor
    location instead of ``b/2``. The shipped ``distance_um`` key is returned
    unchanged alongside it, so the shipped estimator can be evaluated on the
    identical bins.
    """
    d = np.asarray(distance_px, dtype=np.float64)
    v = np.asarray(value_px, dtype=np.float64)
    n_bins = max(int(math.ceil(float(max_um) / float(bin_um))), 1)
    edges = np.arange(n_bins + 1, dtype=np.float64) * float(bin_um)
    idx = np.digitize(d, edges) - 1
    centre, dmean, mean, count = [], [], [], []
    for b in range(n_bins):
        sel = idx == b
        n = int(sel.sum())
        if n == 0:
            continue
        centre.append(0.5 * (edges[b] + edges[b + 1]))
        dmean.append(float(d[sel].mean()))
        mean.append(float(v[sel].mean()))
        count.append(n)
    cnt = np.asarray(count, dtype=int)
    keep = cnt >= VP.MIN_BIN_COUNT
    return {
        "distance_um": np.asarray(centre)[keep],
        "distance_mean_um": np.asarray(dmean)[keep],
        "mean": np.asarray(mean)[keep],
        "count": cnt[keep],
        # The SAME bins with the MIN_BIN_COUNT filter NOT applied. MIN_BIN_COUNT
        # is a fixed absolute count applied to bins whose population scales with
        # the bin width, so it silently shortens the analysed distance range at
        # small bin widths. Keeping both lets that coupling be measured instead
        # of argued.
        "all_distance_um": np.asarray(centre),
        "all_mean": np.asarray(mean),
        "all_count": cnt,
        "n_dropped_by_min_count": int((~keep).sum()),
        # ALIASING DIAGNOSTIC. The distance map is sampled on a lattice of pitch
        # ``dx``, so when ``bin_um < dx`` the pixel-distance histogram is a comb:
        # some bins straddle a populated shell and some fall between shells and
        # are empty or near-empty. ``distance_profile`` drops empty bins SILENTLY
        # (``if n == 0: continue``), so an aliased profile is indistinguishable
        # from a well-sampled one in the shipped output. ``BIN_UM = 10`` is
        # documented as "comfortably above dx = 12 um" — which stops being true
        # the moment the dx = 30 um arm is run at the same bin width.
        "n_empty_bins": int(n_bins - cnt.size),
        "count_cv": (float(cnt.std() / cnt.mean()) if cnt.size and cnt.mean()
                     else float("nan")),
        "bin_um": float(bin_um),
    }


# --------------------------------------------------------------------------- #
# 3. STAGE `bias` — the closed form, verified against re-binned synthetics
# --------------------------------------------------------------------------- #
def stage_bias(verbose: bool = True) -> Dict:
    """Demonstrate the bin-width bias MECHANICALLY, on profiles whose true answer
    is known by construction.

    Two synthetic profiles are built from the same pixel cloud geometry and
    re-binned at every nuisance bin width:

    * ``pure`` — ``C = exp(-d/L)``, whose true half-distance is ``L*ln2`` exactly.
    * ``floor`` — ``C = (1-f)*exp(-d/L) + f``, the shape the real profile actually
      has (the shipped fit reports ``floor_fraction ~= 0.12`` on the flagship
      section), whose true half-distance is also known in closed form.

    No fitting, no tissue, no PDE. If the shipped estimator's error tracks the
    closed form on these, the bin-width sensitivity measured on the real section
    is the estimator's, not the tissue's.
    """
    L = 50.0
    out: Dict = {"L_true_um": L, "cases": []}
    # A 1-D pixel cloud on a fine lattice, uniform in distance — the geometry
    # under which the closed form is exact.
    d_px = np.arange(0.0, 400.0, 0.05) + 0.025
    for floor_frac in (0.0, 0.12, 0.30):
        c_px = (1.0 - floor_frac) * np.exp(-d_px / L) + floor_frac
        # True half-distance of THIS profile, from its own d=0 value.
        c0 = 1.0
        target = 0.5 * c0
        if target > floor_frac:
            d_true = float(-L * math.log((target - floor_frac) / (1.0 - floor_frac)))
        else:
            d_true = float("nan")
        case = {"floor_fraction": floor_frac, "d_half_true_um": d_true, "bins": []}
        for b in NUISANCE_BINS_UM:
            prof = binned_profile(d_px, c_px, b, max_um=400.0)
            shipped = VP.empirical_half_distance_um(prof["distance_um"], prof["mean"])
            corr_centre = anchor_corrected_half_distance_um(
                prof["distance_um"], prof["mean"])
            corr_mean = anchor_corrected_half_distance_um(
                prof["distance_um"], prof["mean"],
                anchor_distance_um=float(prof["distance_mean_um"][0]))
            free = exposure_weighted_depth(d_px, c_px)
            case["bins"].append({
                "bin_um": float(b),
                "shipped_d_half_um": shipped,
                "shipped_minus_true_um": (shipped - d_true) if shipped else None,
                "closed_form_bias_um": 0.5 * float(b),
                "anchor_corrected_centre_um": corr_centre,
                "anchor_corrected_mean_um": corr_mean,
                "binfree_depth_median_um": free["depth_median_um"],
            })
        out["cases"].append(case)
        if verbose:
            print(f"\n  floor fraction {floor_frac:.2f} — true d½ = {d_true:.3f} µm")
            print(f"    {'bin':>5} {'shipped':>9} {'err':>8} {'b/2':>6} "
                  f"{'corr(centre)':>13} {'corr(mean)':>11} {'bin-free med':>13}")
            for r in case["bins"]:
                print(f"    {r['bin_um']:5.0f} {r['shipped_d_half_um']:9.3f} "
                      f"{r['shipped_minus_true_um']:+8.3f} "
                      f"{r['closed_form_bias_um']:6.1f} "
                      f"{r['anchor_corrected_centre_um']:13.3f} "
                      f"{r['anchor_corrected_mean_um']:11.3f} "
                      f"{r['binfree_depth_median_um']:13.3f}")
    return out


# --------------------------------------------------------------------------- #
# 4. STAGE `fields` — solve once, cache, never solve again
# --------------------------------------------------------------------------- #
def _arm_paths(arm_key: str) -> str:
    return os.path.join(_CACHE_DIR, f"{arm_key}.npz")


def stage_fields(verbose: bool = True) -> Dict:
    """Solve the pair on each tissue arm ONCE and cache every array the
    downstream stages need. Six solves; everything after this is free."""
    alpha = _assert_alpha()
    compiled = compile_pair(alpha)
    os.makedirs(_CACHE_DIR, exist_ok=True)
    pde = PDEConfig(**_PDE)
    meta: Dict = {"alpha": alpha, "arms": {}}
    for arm_key in NUISANCE_ARMS:
        arm = next(a for a in VP.ARMS if a["key"] == arm_key)
        t0 = time.time()
        micro, conn_native, conn_used = VP.build_arm_microenvironment(
            dx_um=arm["dx_um"], repair=arm["repair"])
        tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
        endo = VP.endothelial_xy(micro.cells)
        dist_p = VP.direct_vessel_distance_um(micro.domain, endo)
        legacy = VP.vessel_mask_from_density(micro.fields.vascular_density, tissue)
        dist_l = VP.distance_from_vessel_um(legacy["mask"], float(arm["dx_um"]))
        vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
        store = {"tissue": tissue, "dist_primary": dist_p, "dist_legacy": dist_l,
                 "vascular_density": vasc}
        # The NATIVE (unrepaired) mask, on the SAME grid, so `--stage repair` can
        # separate "the field changed" from "the pixel set changed".
        if arm["repair"] > 0:
            base = VP.reconstruct_real(source="xenium", dx_um=float(arm["dx_um"]))
            store["tissue_native"] = np.asarray(base.domain.tissue_mask, dtype=bool)
        for name in PAIR:
            t1 = time.time()
            res = simulate_response(micro, SimConfig(drug=compiled[name], pde=pde,
                                                     seed=0))
            store[f"auc__{name}"] = np.asarray(res.fields["auc"], dtype=np.float64)
            store[f"conc__{name}"] = np.asarray(res.transport.final_concentration,
                                                dtype=np.float64)
            if verbose:
                print(f"  {arm_key:<15} {name:<14} solved in "
                      f"{time.time() - t1:5.1f}s")
        np.savez_compressed(_arm_paths(arm_key), **store)
        meta["arms"][arm_key] = {
            "dx_um": float(arm["dx_um"]), "repair": int(arm["repair"]),
            "n_tissue_px": int(tissue.sum()), "n_endothelial": int(endo.shape[0]),
            "conn_native": conn_native, "conn_used": conn_used,
            "legacy_area_fraction": legacy["area_fraction"],
            "legacy_rule": legacy["rule"],
            "elapsed_s": time.time() - t0,
        }
        if verbose:
            print(f"  {arm_key:<15} cached ({time.time() - t0:.0f}s total)")
    return meta


def _load_arm(arm_key: str) -> Dict:
    p = _arm_paths(arm_key)
    if not os.path.exists(p):
        raise SystemExit(
            f"missing field cache {p}; run `--stage fields` first")
    z = np.load(p)
    return {k: z[k] for k in z.files}


# --------------------------------------------------------------------------- #
# 5. STAGE `decompose` — replay every reduction choice on the cached fields
# --------------------------------------------------------------------------- #
def stage_decompose(verbose: bool = True) -> Dict:
    """Every (arm x field x metric x bin x truncation) cell, every estimator.

    The four shipped-grid axes are reproduced exactly; ``truncation`` is the axis
    ``kyle2004_penetration.md`` §8 left implicit inside ``distance_profile``.
    """
    rows: List[Dict] = []
    for arm_key in NUISANCE_ARMS:
        z = _load_arm(arm_key)
        tissue = z["tissue"].astype(bool)
        for name in PAIR:
            for fk in NUISANCE_FIELDS:
                field = z[f"{fk}__{name}"]
                fv = field[tissue]
                for mk in NUISANCE_METRICS:
                    dv = z[f"dist_{mk}"][tissue]
                    for trunc in NUISANCE_TRUNC_PCT:
                        max_um = float(np.percentile(dv, trunc))
                        sel = dv <= max_um
                        d_in, v_in = dv[sel], fv[sel]
                        free = exposure_weighted_depth(d_in, v_in)
                        base = {
                            "drug": name, "arm": arm_key, "field": fk,
                            "metric": mk, "trunc_pct": float(trunc),
                            "max_um": max_um,
                            "binfree_depth_mean_um": free.get("depth_mean_um"),
                            "binfree_depth_median_um": free.get("depth_median_um"),
                            "binfree_depth_ratio": free.get("depth_ratio"),
                            "binfree_depth_median_ratio":
                                free.get("depth_median_ratio"),
                            "binfree_tissue_mean_um": free.get("tissue_mean_um"),
                            "pixel_log_L_um": pixel_log_slope_L_um(d_in, v_in),
                            "p90_p10": p90_p10(v_in),
                        }
                        for b in NUISANCE_BINS_UM:
                            prof = binned_profile(d_in, v_in, b, max_um)
                            d, y, w = (prof["distance_um"], prof["mean"],
                                       prof["count"].astype(float))
                            shipped = VP.empirical_half_distance_um(d, y)
                            fit = (VP.fit_exponential_plus_floor(d, y, w)
                                   if d.size >= 4 else {"error": "too few bins"})
                            # The SAME fit, with each bin plotted at the distance
                            # its own mean value actually corresponds to instead
                            # of at the bin's geometric centre. The pixels inside
                            # a bin are not uniform in distance — on a 2-D
                            # section the tissue area at distance d rises with d
                            # — so the centre sits systematically NEARER the
                            # vessel than the value it is labelling, by an amount
                            # that grows with the bin width. Correctness test is
                            # convergence, not score: the b -> 0 limit is the
                            # definition, and the labelling that agrees with it
                            # at coarse b is the right one.
                            fit_dbar = (
                                VP.fit_exponential_plus_floor(
                                    prof["distance_mean_um"], y, w)
                                if d.size >= 4 else {"error": "too few bins"})
                            # The MATCHED comparator for pixel_log_L_um: the same
                            # pure-exponential log-space estimator, binned. The
                            # exp+floor fit is a DIFFERENT MODEL (three
                            # parameters, linear space), so comparing it against
                            # the unbinned log fit would confound "binning" with
                            # "floor term" — a pure exponential over a profile
                            # with a plateau reads systematically longer, which
                            # is exactly what vessel_profile documents for
                            # near_exp vs exp_floor.
                            fit_log = VP.fit_exponential(d, y, w)
                            da, ya, wa = (prof["all_distance_um"],
                                          prof["all_mean"],
                                          prof["all_count"].astype(float))
                            fit_nomin = (VP.fit_exponential_plus_floor(da, ya, wa)
                                         if da.size >= 4
                                         else {"error": "too few bins"})
                            rows.append({
                                **base, "bin_um": float(b),
                                "n_bins": int(d.size),
                                "shipped_d_half_um": shipped,
                                "anchor_corrected_centre_um":
                                    anchor_corrected_half_distance_um(d, y),
                                "anchor_corrected_mean_um":
                                    anchor_corrected_half_distance_um(
                                        d, y,
                                        anchor_distance_um=float(
                                            prof["distance_mean_um"][0])
                                        if d.size else None),
                                "fit_L_um": fit.get("L_um"),
                                "fit_floor_fraction": fit.get("floor_fraction"),
                                "fit_r2": fit.get("r2_linear_space"),
                                "fit_error": fit.get("error"),
                                # The window MIN_BIN_COUNT actually left for the
                                # fit. It is a FIXED ABSOLUTE count applied to
                                # bins whose population scales with bin width,
                                # so the analysed distance range is itself a
                                # function of the bin width — a coupling the
                                # shipped code never sweeps.
                                "fit_range_hi_um": (float(fit["fit_range_um"][1])
                                                    if "fit_range_um" in fit
                                                    else None),
                                "profile_hi_um": (float(d[-1]) if d.size
                                                  else None),
                                "usable_fraction_of_max": (
                                    float(d[-1] / max_um) if d.size and max_um
                                    else None),
                                "n_dropped_by_min_count":
                                    prof["n_dropped_by_min_count"],
                                "n_empty_bins": prof["n_empty_bins"],
                                "count_cv": prof["count_cv"],
                                "first_bin_count": (int(prof["count"][0])
                                                    if prof["count"].size else 0),
                                "first_bin_mean_d_um": (
                                    float(prof["distance_mean_um"][0])
                                    if prof["distance_mean_um"].size else None),
                                "fit_L_no_min_count_um": fit_nomin.get("L_um"),
                                "fit_L_binmean_um": fit_dbar.get("L_um"),
                                "fit_L_log_binned_um": fit_log.get("L_um"),
                                "fit_binmean_floor_fraction":
                                    fit_dbar.get("floor_fraction"),
                                "fit_binmean_r2": fit_dbar.get("r2_linear_space"),
                                "shipped_d_half_no_min_count_um":
                                    VP.empirical_half_distance_um(da, ya),
                            })
        if verbose:
            print(f"  {arm_key:<15} replayed")
    return {"rows": rows,
            "axes": {"arm": list(NUISANCE_ARMS), "field": list(NUISANCE_FIELDS),
                     "metric": list(NUISANCE_METRICS),
                     "bin_um": list(NUISANCE_BINS_UM),
                     "trunc_pct": list(NUISANCE_TRUNC_PCT)}}


# --------------------------------------------------------------------------- #
# 5b. STAGE `converge` — which labelling converges to the fine-bin limit?
# --------------------------------------------------------------------------- #
#: A fine-to-coarse bin sweep at the flagship configuration. This is the
#: CORRECTNESS test, and it is deliberately not a score: as ``bin_um -> 0`` the
#: binning stops doing anything and the profile becomes the pixel cloud itself,
#: so the fine-bin value IS the definition of the quantity. An estimator is right
#: if it still agrees with that limit at coarse bins, and wrong if it drifts away
#: from it. No estimator is preferred here because it produces a number anyone
#: wanted.
CONVERGENCE_BINS_UM: Sequence[float] = (2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 15.0,
                                        20.0, 30.0, 40.0)


def stage_converge(arm_key: str = "dx12_repaired", metric: str = "primary",
                   verbose: bool = True) -> Dict:
    z = _load_arm(arm_key)
    tissue = z["tissue"].astype(bool)
    out: Dict = {"arm": arm_key, "metric": metric, "drugs": {}}
    for name in PAIR:
        dv = z[f"dist_{metric}"][tissue]
        fv = z[f"auc__{name}"][tissue]
        max_um = float(np.percentile(dv, 99.5))
        sel = dv <= max_um
        d_in, v_in = dv[sel], fv[sel]
        rows = []
        for b in CONVERGENCE_BINS_UM:
            prof = binned_profile(d_in, v_in, b, max_um)
            d, y, w = prof["distance_um"], prof["mean"], prof["count"].astype(float)
            if d.size < 4:
                continue
            rows.append({
                "bin_um": float(b), "n_bins": int(d.size),
                "n_empty_bins": prof["n_empty_bins"],
                "first_bin_mean_d_um": float(prof["distance_mean_um"][0]),
                "shipped_d_half_um": VP.empirical_half_distance_um(d, y),
                "anchor_corrected_mean_um": anchor_corrected_half_distance_um(
                    d, y, anchor_distance_um=float(prof["distance_mean_um"][0])),
                "fit_L_centre_um": VP.fit_exponential_plus_floor(
                    d, y, w).get("L_um"),
                "fit_L_binmean_um": VP.fit_exponential_plus_floor(
                    prof["distance_mean_um"], y, w).get("L_um"),
            })
        out["drugs"][name] = rows
        if verbose:
            print(f"\n  {name} — {arm_key}, {metric} metric, AUC")
            print(f"    {'bin':>5} {'empty':>6} {'d̄₀':>6} {'d½ shipped':>11} "
                  f"{'d½ anchored':>12} {'L centre':>9} {'L binmean':>10}")
            for r in rows:
                print(f"    {r['bin_um']:5.1f} {r['n_empty_bins']:6d} "
                      f"{r['first_bin_mean_d_um']:6.2f} "
                      f"{r['shipped_d_half_um']:11.2f} "
                      f"{r['anchor_corrected_mean_um']:12.2f} "
                      f"{r['fit_L_centre_um']:9.2f} {r['fit_L_binmean_um']:10.2f}")
    return out


# --------------------------------------------------------------------------- #
# 6. STAGE `repair` — is mask repair a FIELD change or a PIXEL-SET change?
# --------------------------------------------------------------------------- #
def stage_repair(verbose: bool = True) -> Dict:
    """Split the mask-repair nuisance into its two mechanisms.

    Repairing the mask does two separable things:

    A. it changes the SOLVE — the added interstitial pixels connect the domain,
       so lateral transport is no longer suppressed by the no-flux boundary at
       every fragment edge, and ``build_real_fields`` re-derives the source on
       the new mask;
    B. it changes the REDUCTION SET — the added pixels enter the distance binning
       and the population statistics, and they are not drawn from the same
       distance distribution as the original pixels.

    Reducing the REPAIRED solve over the NATIVE pixel set isolates A; comparing
    the two reductions of the same repaired field isolates B.
    """
    rep = _load_arm("dx12_repaired")
    nat = _load_arm("dx12_native")
    t_rep = rep["tissue"].astype(bool)
    t_nat_on_rep = rep["tissue_native"].astype(bool)
    t_nat = nat["tissue"].astype(bool)
    if t_nat_on_rep.shape != t_nat.shape or not np.array_equal(t_nat_on_rep, t_nat):
        raise AssertionError(
            "the native mask cached inside the repaired arm is not the same array "
            "as the native arm's own mask; the repair decomposition is not valid")
    added = t_rep & ~t_nat
    out: Dict = {
        "n_native_px": int(t_nat.sum()), "n_repaired_px": int(t_rep.sum()),
        "n_added_px": int(added.sum()),
        "added_fraction": float(added.sum()) / float(t_rep.sum()),
        "drugs": {},
    }
    for mk in NUISANCE_METRICS:
        d_rep, d_nat = rep[f"dist_{mk}"], nat[f"dist_{mk}"]
        out.setdefault("distance", {})[mk] = {
            "mean_native_px_um": float(d_rep[t_nat].mean()),
            "mean_added_px_um": float(d_rep[added].mean())
            if added.any() else float("nan"),
            "median_native_px_um": float(np.median(d_rep[t_nat])),
            "median_added_px_um": float(np.median(d_rep[added]))
            if added.any() else float("nan"),
            "identical_maps": bool(np.allclose(d_rep, d_nat)),
        }
    for name in PAIR:
        f_rep = rep[f"auc__{name}"]
        f_nat = nat[f"auc__{name}"]
        rec: Dict = {}
        for mk in NUISANCE_METRICS:
            d_rep = rep[f"dist_{mk}"]
            d_nat = nat[f"dist_{mk}"]

            def _read(fld, dist, mask, label):
                dv, fv = dist[mask], fld[mask]
                max_um = float(np.percentile(dv, 99.5))
                sel = dv <= max_um
                prof = binned_profile(dv[sel], fv[sel], VP.BIN_UM, max_um)
                free = exposure_weighted_depth(dv[sel], fv[sel])
                return {
                    "label": label,
                    "shipped_d_half_um": VP.empirical_half_distance_um(
                        prof["distance_um"], prof["mean"]),
                    "anchor_corrected_mean_um":
                        anchor_corrected_half_distance_um(
                            prof["distance_um"], prof["mean"],
                            anchor_distance_um=float(prof["distance_mean_um"][0])),
                    "binfree_depth_median_um": free["depth_median_um"],
                    "binfree_depth_median_ratio": free["depth_median_ratio"],
                    "p90_p10": p90_p10(fv[sel]),
                    "mean_field": float(fv.mean()),
                }
            rec[mk] = {
                "native_solve_native_px": _read(f_nat, d_nat, t_nat,
                                                "native solve, native pixels"),
                "repaired_solve_native_px": _read(f_rep, d_rep, t_nat,
                                                  "repaired solve, native pixels"),
                "repaired_solve_repaired_px": _read(f_rep, d_rep, t_rep,
                                                    "repaired solve, repaired px"),
            }
        out["drugs"][name] = rec
        if verbose:
            r = rec["primary"]
            print(f"  {name}")
            for k in ("native_solve_native_px", "repaired_solve_native_px",
                      "repaired_solve_repaired_px"):
                v = r[k]
                print(f"    {v['label']:<32} d½={v['shipped_d_half_um']:6.2f}  "
                      f"binfree_med={v['binfree_depth_median_um']:6.2f}  "
                      f"P90/P10={v['p90_p10']:6.3f}")
    return out


# --------------------------------------------------------------------------- #
# 7. Attribution + separation scoring
# --------------------------------------------------------------------------- #
ESTIMATORS: Sequence[Tuple[str, str]] = (
    ("shipped_d_half_um", "empirical half-distance (SHIPPED)"),
    ("anchor_corrected_centre_um", "half-distance, anchor-corrected (bin centre)"),
    ("anchor_corrected_mean_um", "half-distance, anchor-corrected (bin mean d)"),
    ("fit_L_um", "L from C0*exp(-d/L)+floor (SHIPPED fit)"),
    ("fit_L_no_min_count_um", "same fit, MIN_BIN_COUNT filter removed"),
    ("fit_L_binmean_um", "same fit, bins plotted at their MEAN distance"),
    ("fit_L_log_binned_um", "L from BINNED log-linear fit (no floor term)"),
    ("pixel_log_L_um", "L from UNBINNED pixel log-linear fit (same model)"),
    ("binfree_depth_mean_um", "bin-free exposure-weighted mean depth"),
    ("binfree_depth_median_um", "bin-free exposure-weighted median depth"),
    ("binfree_depth_ratio", "bin-free mean depth / tissue mean depth"),
    ("binfree_depth_median_ratio", "bin-free median depth / tissue median depth"),
    ("p90_p10", "P90/P10 over the same pixels (no distance at all)"),
)

#: Axes a bin-free estimator is structurally blind to — reported as such rather
#: than as a zero measurement.
_BLIND: Dict[str, Sequence[str]] = {
    "binfree_depth_mean_um": ("bin_um",),
    "binfree_depth_median_um": ("bin_um",),
    "binfree_depth_ratio": ("bin_um",),
    "binfree_depth_median_ratio": ("bin_um",),
    "pixel_log_L_um": ("bin_um",),
    "p90_p10": ("bin_um",),
}


def _vals(rows, drug, key):
    v = [r.get(key) for r in rows if r["drug"] == drug]
    return np.array([x for x in v if x is not None and np.isfinite(x)],
                    dtype=np.float64)


def attribute(dec: Dict) -> Dict:
    """For every estimator: its span over the full grid, its PAIRED sensitivity to
    each axis (the median and max range across that axis with every other axis
    held fixed), and whether the two compounds stay disjoint."""
    rows = dec["rows"]
    axes = ["arm", "field", "metric", "bin_um", "trunc_pct"]
    out: Dict = {"estimators": {}}
    for key, label in ESTIMATORS:
        rec: Dict = {"label": label, "drugs": {}, "axes": {}}
        for drug in PAIR:
            R = [r for r in rows if r["drug"] == drug]
            v = _vals(rows, drug, key)
            n_missing = sum(1 for r in R
                            if r.get(key) is None or not np.isfinite(
                                r.get(key) or float("nan")))
            rec["drugs"][drug] = {
                "n": int(v.size), "n_missing": int(n_missing),
                "min": float(v.min()) if v.size else None,
                "max": float(v.max()) if v.size else None,
                "span": float(v.max() - v.min()) if v.size else None,
                "median": float(np.median(v)) if v.size else None,
            }
        a, b = (rec["drugs"][DOX], rec["drugs"][DAU])
        if a["min"] is not None and b["min"] is not None:
            disjoint = (a["min"] > b["max"]) or (b["min"] > a["max"])
            gap = (max(a["min"], b["min"]) - min(a["max"], b["max"]))
            sep = abs(a["median"] - b["median"])
            worst = max(a["span"], b["span"])
            rec["separation"] = {
                "disjoint": bool(disjoint),
                "gap": float(gap),
                "median_separation": float(sep),
                "worst_span": float(worst),
                "separation_over_span": float(sep / worst) if worst > 0
                else float("inf"),
                "central_dox": a["median"], "central_dau": b["median"],
            }
        # PAIRED analysis. The registered C2 bar in kyle2004_penetration.md §8 is
        # that the two compounds' MARGINAL intervals over the whole grid be
        # disjoint — which compares doxorubicin at one nuisance setting against
        # daunorubicin at a DIFFERENT one. If a nuisance moves both compounds the
        # same way (common mode), it inflates both marginals without touching the
        # comparison any single run actually makes, where the setting is shared.
        # Reporting the paired difference alongside the marginal one is not a
        # softer bar, it is a different question, and both are printed.
        pair_rows: Dict = {}
        for r in rows:
            k = tuple(r[a] for a in axes)
            pair_rows.setdefault(k, {})[r["drug"]] = r.get(key)
        diffs, ratios, n_inv = [], [], 0
        for k, dd in pair_rows.items():
            a_, b_ = dd.get(DOX), dd.get(DAU)
            if a_ is None or b_ is None:
                continue
            if not (np.isfinite(a_) and np.isfinite(b_)):
                continue
            diffs.append(float(a_ - b_))
            if b_ != 0:
                ratios.append(float(a_ / b_))
            if a_ <= b_:
                n_inv += 1
        if diffs:
            rec["paired"] = {
                "n_cells": len(diffs),
                "diff_min": float(np.min(diffs)),
                "diff_max": float(np.max(diffs)),
                "diff_median": float(np.median(diffs)),
                "ratio_min": float(np.min(ratios)) if ratios else None,
                "ratio_max": float(np.max(ratios)) if ratios else None,
                "ratio_median": float(np.median(ratios)) if ratios else None,
                "n_inversions": int(n_inv),
                "sign_stable": bool(n_inv == 0),
            }
        for ax in axes:
            if ax in _BLIND.get(key, ()):
                rec["axes"][ax] = {"blind": True}
                continue
            others = [o for o in axes if o != ax]
            grp: Dict = {}
            for r in rows:
                if r["drug"] != DOX:
                    continue
                val = r.get(key)
                if val is None or not np.isfinite(val):
                    continue
                grp.setdefault(tuple(r[o] for o in others), []).append(val)
            rngs = [max(g) - min(g) for g in grp.values() if len(g) > 1]
            rec["axes"][ax] = {
                "blind": False,
                "median_range": float(np.median(rngs)) if rngs else 0.0,
                "max_range": float(np.max(rngs)) if rngs else 0.0,
            }
        out["estimators"][key] = rec
    return out


def print_attribution(att: Dict) -> None:
    print("\n" + "=" * 100)
    print("STAGE ATTRIBUTION — paired range of the doxorubicin readout across each "
          "axis, others held fixed")
    print("=" * 100)
    hdr = (f"{'estimator':<44} {'arm':>9} {'field':>8} {'metric':>9} "
           f"{'bin':>9} {'trunc':>9}")
    print(hdr)
    print("-" * 100)
    for key, label in ESTIMATORS:
        rec = att["estimators"][key]
        cells = []
        for ax in ("arm", "field", "metric", "bin_um", "trunc_pct"):
            a = rec["axes"][ax]
            cells.append("  blind" if a.get("blind")
                         else f"{a['median_range']:9.3f}")
        print(f"{label[:43]:<44} " + " ".join(f"{c:>9}" for c in cells))
    print("\n" + "=" * 100)
    print("SEPARATION — does the estimator keep doxorubicin and daunorubicin apart "
          "over the WHOLE grid?")
    print("=" * 100)
    print(f"{'estimator':<44} {'DOX central':>12} {'DAU central':>12} "
          f"{'worst span':>11} {'sep/span':>9} {'disjoint':>9} {'missing':>8}")
    print("-" * 100)
    for key, label in ESTIMATORS:
        rec = att["estimators"][key]
        s = rec.get("separation")
        miss = (rec["drugs"][DOX]["n_missing"] + rec["drugs"][DAU]["n_missing"])
        if not s:
            print(f"{label[:43]:<44} {'n/a':>12}")
            continue
        print(f"{label[:43]:<44} {s['central_dox']:12.3f} {s['central_dau']:12.3f} "
              f"{s['worst_span']:11.3f} {s['separation_over_span']:9.3f} "
              f"{('YES' if s['disjoint'] else 'no'):>9} {miss:8d}")

    print("\n" + "=" * 100)
    print("PAIRED SEPARATION — doxorubicin MINUS daunorubicin computed WITHIN each "
          "nuisance cell")
    print("=" * 100)
    print(f"{'estimator':<44} {'diff min':>10} {'diff med':>10} {'diff max':>10} "
          f"{'ratio min':>10} {'ratio max':>10} {'inversions':>11}")
    print("-" * 100)
    for key, label in ESTIMATORS:
        p = att["estimators"][key].get("paired")
        if not p:
            continue
        print(f"{label[:43]:<44} {p['diff_min']:10.3f} {p['diff_median']:10.3f} "
              f"{p['diff_max']:10.3f} "
              f"{(p['ratio_min'] if p['ratio_min'] is not None else float('nan')):10.4f} "
              f"{(p['ratio_max'] if p['ratio_max'] is not None else float('nan')):10.4f} "
              f"{p['n_inversions']:>4d}/{p['n_cells']:<6d}")


# --------------------------------------------------------------------------- #
# 8. Entry point
# --------------------------------------------------------------------------- #
def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", default="all",
                    choices=("bias", "fields", "decompose", "converge",
                             "repair", "all"))
    args = ap.parse_args(argv)
    t0 = time.time()
    payload: Dict = {"alpha": _assert_alpha()}

    if args.stage in ("bias", "all"):
        print("\n" + "=" * 100)
        print("STAGE 1 — the estimator's bin-width bias, on profiles with a KNOWN "
              "answer")
        print("=" * 100)
        payload["bias"] = stage_bias()

    if args.stage in ("fields", "all"):
        print("\n" + "=" * 100)
        print("STAGE 2 — solving the pair on every tissue arm (cached)")
        print("=" * 100)
        payload["fields"] = stage_fields()

    if args.stage in ("decompose", "all"):
        print("\n" + "=" * 100)
        print("STAGE 3 — replaying every reduction choice on the cached fields")
        print("=" * 100)
        dec = stage_decompose()
        att = attribute(dec)
        print_attribution(att)
        payload["decompose"] = dec
        payload["attribution"] = att

    if args.stage in ("converge", "all"):
        print("\n" + "=" * 100)
        print("STAGE 3b — bin-width CONVERGENCE: which labelling agrees with the "
              "fine-bin limit?")
        print("=" * 100)
        payload["converge"] = stage_converge()

    if args.stage in ("repair", "all"):
        print("\n" + "=" * 100)
        print("STAGE 4 — mask repair: field change or pixel-set change?")
        print("=" * 100)
        payload["repair"] = stage_repair()

    payload["elapsed_s"] = time.time() - t0
    if args.stage != "fields":
        prev = {}
        if os.path.exists(_OUT_JSON):
            with open(_OUT_JSON) as fh:
                prev = json.load(fh)
        prev.update(payload)
        with open(_OUT_JSON, "w") as fh:
            json.dump(prev, fh, indent=1, default=float)
        print(f"\nwrote {_OUT_JSON}  ({time.time() - t0:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
