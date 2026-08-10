"""Score the SOLVED interstitial fluid pressure field against the one published
human intracranial-tumour IFP measurement.

The anchor
----------
Boucher Y, Salehi H, Witwer B, Harsh GR, Jain RK. "Interstitial fluid pressure
in intracranial tumours in patients and in rodents." *Br J Cancer*
1997;75(6):829-836. **PMID 9062403**, PMC2063404, doi:10.1038/bjc.1997.148.
Retrieved from PubMed/PMC and read in full; every number below is quoted from
the paper, not from a secondary summary.

Results, "Measurements in patients", verbatim:

    "The IFP was measured in 11 patients with an astrocytoma, glioblastoma,
     meningioma or ganglioglioma (Table 3). The IFP in 10 of 11 tumours varied
     between 0.0 and 3.0 mmHg. In one astrocytoma the IFP was 9.5 mmHg. The
     mean IFP was 2.0 +/- 2.5 mmHg."

Statistical analysis, verbatim: "The data are presented as the mean +/-
standard deviation." So **2.5 is an SD, not an SEM**, and the +/- interval is a
spread over 11 patients rather than a confidence interval on the mean. That is
recorded in :data:`ANCHOR` and is the single most important qualifier on any
verdict this file produces: a +/- 1 SD band on a non-negative quantity with mean
2.0 admits everything from 0 to 4.5 mmHg, so passing it is easy.

Table 3, transcribed in full (tumour IFP, mmHg, and the brain surrounding the
tumour where it was also taken):

    M/30  R-Glioblastoma      X       3.0   3.5
    M/39  R-Glioblastoma      X       1.5   2.5
    F/51  Glioblastoma        -       0.5   0.5
    F/48  Astrocytoma         -       9.5   1.5
    F/70  Meningioma          -       1.5    -
    M/48  R-Glioblastoma      S,X     0.0    -
    F/36  R-Glioblastoma      S,X     1.0   3.0
    M/49  Oligo/Astrocytoma   -       2.5    -
    F/30  R-Astrocytoma       S,C,X   2.5    -
    M/36  Ganglioglioma       -       0.5    -
    F/62  R-Astrocytoma       X,C,S   2.0    -

Technique: **wick-in-needle** (Boucher et al 1991), taken intraoperatively
"following the opening of the intracranial cavity and removal of the dura".

Confounds the paper states about itself:
  * every patient had dexamethasone 10 mg + mannitol 0.5 mg/kg + furosemide
    20 mg i.v. before surgery, plus hyperventilation to PaCO2 25-32 mmHg;
  * the authors' own conclusion is that "antioedema therapy in patients with
    brain tumours is responsible in part for the low tumour IFP";
  * the same paper measures a 17-50% fall in rodent tumour IFP from hypocapnia
    alone (27 +/- 15 -> 18 +/- 10 mmHg).
So 2.0 +/- 2.5 mmHg is a **treated, decompressed, open-cranium** number, not an
untreated intracranial-tumour baseline. The engine models none of that.

What this file computes
-----------------------
1. The commensurable statistic. A wick-in-needle probe is a 23-gauge needle
   whose side hole equilibrates with the mobile interstitial fluid over a
   sub-millimetre to millimetre contact region. It is therefore a **local
   volume mean at one placement inside the tumour mass** -- not a domain mean,
   not a field maximum, and emphatically not a rendering contour level. The
   model-side statistic is the solved field averaged over a wick-sized Gaussian
   window (a ladder of 125 / 250 / 500 um is swept, because the exact sampling
   volume is not derivable from the paper), restricted to tumour-bearing
   tissue, and reported as a distribution over placements.
2. A numerical-resolution diagnostic. The screened-Poisson operator has a
   screening length L = sqrt(K / (a + d)); if L < dx the mesh cannot resolve
   the operator and each node relaxes to its own local algebraic equilibrium
   a*p_v/(a+d) with the Laplacian playing almost no role. That ratio is
   measured, not assumed.
3. A grid ladder, because neither flagship panel is grid-converged
   (results/grid_convergence_breast.md, results/grid_convergence_gbm.md).

NOTHING IS TUNED. No coefficient in tmesim/physics/pressure.py is read, written
or fitted here; the module is imported and called exactly as it ships.

Run:
    ./.venv/bin/python experiments/ifp_validation.py
    ./.venv/bin/python experiments/ifp_validation.py --dx 24 12
    ./.venv/bin/python experiments/ifp_validation.py --report-only
Writes:
    results/sweep/ifp_validation.json     (raw records, restart-safe)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy import ndimage

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.physics import pressure as P                       # noqa: E402
from tmesim.physics.pressure import solve_ifp                  # noqa: E402
from tmesim.reconstruct.real_pipeline import reconstruct_real  # noqa: E402

# --------------------------------------------------------------------------- #
# The anchor, as read from the primary source.
# --------------------------------------------------------------------------- #
#: Per-patient tumour IFP in mmHg, Table 3, in table order.
BOUCHER_TUMOUR_IFP_MMHG = (3.0, 1.5, 0.5, 9.5, 1.5, 0.0, 1.0, 2.5, 2.5, 0.5, 2.0)
#: The same rows' tumour types, so the glioma-only subset can be taken.
BOUCHER_TUMOUR_TYPE = ("R-Glioblastoma", "R-Glioblastoma", "Glioblastoma",
                       "Astrocytoma", "Meningioma", "R-Glioblastoma",
                       "R-Glioblastoma", "Oligo/Astrocytoma", "R-Astrocytoma",
                       "Ganglioglioma", "R-Astrocytoma")

ANCHOR = {
    "pmid": "9062403",
    "pmc": "PMC2063404",
    "doi": "10.1038/bjc.1997.148",
    "citation": ("Boucher Y, Salehi H, Witwer B, Harsh GR, Jain RK. "
                 "Interstitial fluid pressure in intracranial tumours in "
                 "patients and in rodents. Br J Cancer 1997;75(6):829-836."),
    "quantity": "human intracranial-tumour interstitial fluid pressure",
    "mean_mmhg": 2.0,
    "spread_mmhg": 2.5,
    #: Quoted from the paper's Statistical analysis section, not inferred.
    "spread_is": "standard deviation",
    "spread_evidence": ('Statistical analysis, verbatim: "The data are '
                        'presented as the mean +/- standard deviation." An SEM '
                        'reading is additionally impossible: SEM 2.5 at n=11 '
                        'would imply SD 8.3, which cannot be produced by 11 '
                        'values of which 10 lie in 0.0-3.0.'),
    "n": 11,
    "technique": ("wick-in-needle (Boucher et al 1991), intraoperative, after "
                  "opening of the intracranial cavity and removal of the dura"),
    "observed_range_10_of_11_mmhg": (0.0, 3.0),
    "single_outlier_mmhg": 9.5,
    "confound": ("all 11 patients received dexamethasone 10 mg + mannitol "
                 "0.5 mg/kg + furosemide 20 mg i.v. pre-operatively and were "
                 "hyperventilated to PaCO2 25-32 mmHg; the authors conclude "
                 "antioedema therapy is responsible in part for the low IFP"),
}

#: The BREAST anchors, same technique, verified against PubMed the same way.
#:
#: These matter for two reasons. (1) The breast section is scored here too, and
#: scoring it against an INTRACRANIAL anchor would be exactly the modality
#: conflation this repo has been bitten by twice; the breast section must be
#: scored against breast measurements. (2) Less 1992 also measured NORMAL
#: tissue with the same wick-in-needle technique and an overlapping author list
#: (Boucher and Jain are on both papers), and got 2 +/- 4 mmHg, n = 11 --
#: statistically indistinguishable from Boucher 1997's intracranial TUMOUR
#: value of 2.0 +/- 2.5 mmHg, n = 11. That coincidence is the single most
#: important fact in this file: an anchor that cannot separate tumour from
#: normal tissue cannot discriminate a tumour-transport model from one with no
#: tumour in it.
BREAST_ANCHOR = {
    "pmid": "1423283",
    "citation": ("Less JR, Posner MC, Boucher Y, Borochovitz D, Wolmark N, "
                 "Jain RK. Interstitial hypertension in human breast and "
                 "colorectal tumors. Cancer Res 1992;52(22):6371-4."),
    "quantity": "human PRIMARY BREAST carcinoma interstitial fluid pressure",
    "mean_mmhg": 15.0, "spread_mmhg": 9.0, "spread_is": "standard deviation",
    "technique": "wick-in-needle",
    "all_neoplasms_mmhg": "20 +/- 13 (SD), range 4-50, n=41 measurements / 21 tumours",
    "normal_tissue_mmhg": "2 +/- 4 (SD), n=11",
}
#: Corroborating, but its uncertainty is an SE and must not be read as an SD.
BREAST_ANCHOR_SECONDARY = {
    "pmid": "7850532", "doi": "10.1007/BF03187139",
    "citation": ("Nathanson SD, Nelson L. Interstitial fluid pressure in breast "
                 "cancer, benign breast conditions, and breast parenchyma. "
                 "Ann Surg Oncol 1994;1(4):333-8."),
    "mean_mmhg": 29.0, "spread_mmhg": 3.0,
    "spread_is": "STANDARD ERROR (the paper writes '29 +/- 3 (SE)'), NOT an SD",
    "n": 25, "quantity": "invasive ductal carcinoma IFP",
    "normal_parenchyma_mmhg": "-0.3 +/- 0.1",
}

#: Sections. GBM is the section the landing page renders.
SECTIONS = {
    "gbm": {
        "label": "10x Xenium human GBM (landing-page section)",
        "sample": os.path.join(_ROOT, "data_cache", "xenium_gbm"),
        "tissue": "GBM",
    },
    "breast": {
        "label": "10x Xenium human breast, Janesick rep1",
        "sample": os.path.join(_ROOT, "data_cache", "xenium_breast_janesick",
                               "rep1"),
        "tissue": "breast",
    },
}

#: Grid ladder. 12 um is the production/landing spacing.
DX_LADDER = (24.0, 18.0, 12.0, 8.0)

#: Wick contact scales swept, um (Gaussian sigma). A 23-gauge needle is 0.64 mm
#: OD with a side hole of order 1 mm, so the sampled region is sub-mm to ~1 mm.
WICK_SIGMA_UM = (125.0, 250.0, 500.0)

#: A node counts as "in the tumour mass" -- i.e. somewhere a surgeon would put
#: the needle -- when the smoothed tumour-occupancy field exceeds this. 0.5 of
#: the section's own maximum occupancy; swept for sensitivity.
TUMOUR_PHI_THRESHOLDS = (0.25, 0.5, 0.75)

OUT_REL = os.path.join("results", "sweep", "ifp_validation.json")


# --------------------------------------------------------------------------- #
def coarse_grain(field, mask, sigma_px):
    """Local tissue-only mean, identical to experiments/export_landing_asset.py.

    Normalised convolution: smoothing the masked field and the mask itself and
    dividing, so p = 0 exterior nodes are not averaged in and no spurious
    falloff is invented at the section edge.
    """
    num = ndimage.gaussian_filter(np.where(mask, field, 0.0), sigma_px,
                                  mode="nearest")
    den = ndimage.gaussian_filter(mask.astype(np.float64), sigma_px,
                                  mode="nearest")
    return np.where(den > 1e-3, num / np.maximum(den, 1e-9), 0.0)


def coefficients(micro):
    """Re-derive solve_ifp's own spatial coefficients, by calling its own code.

    ``_tumor_field`` and the module constants are imported rather than copied,
    so this diagnostic cannot drift from the solver it is diagnosing.
    """
    ny, nx = micro.domain.shape
    vasc = np.clip(np.asarray(micro.fields.vascular_density, float), 0.0, None)
    ecm = np.asarray(micro.fields.ecm_density, float)
    ecm = np.clip(ecm, 0.0, 1.0) if ecm.shape == (ny, nx) else np.zeros((ny, nx))
    phi = P._tumor_field(micro)
    K = P._K0 * (1.0 - P._K_ECM_ATTEN * ecm)
    a = vasc * P._LP_NORMAL * (1.0 + P._LEAK_GAIN * phi)
    d = P._LAMBDA_LYMPH * (1.0 - phi)
    react = a + d
    with np.errstate(divide="ignore", invalid="ignore"):
        local = np.where(react > 0, a * P._PV / np.maximum(react, 1e-300), 0.0)
        screen_um = np.sqrt(K / np.maximum(react, 1e-300))
    return {"vasc": vasc, "ecm": ecm, "phi": phi, "K": K, "a": a, "d": d,
            "react": react, "local_equilibrium": local, "screen_um": screen_um}


def _pct(x, qs=(5, 25, 50, 75, 95, 99)):
    if x.size == 0:
        return {f"p{q}": None for q in qs}
    v = np.percentile(x, qs)
    return {f"p{q}": float(vv) for q, vv in zip(qs, v)}


def measure_section(key, dx):
    """Solve IFP on one section at one grid spacing and reduce it to statistics."""
    spec = SECTIONS[key]
    t0 = time.time()
    micro = reconstruct_real(source="xenium", sample=spec["sample"], dx_um=dx)
    t_recon = time.time() - t0
    t0 = time.time()
    ifp, velocity = solve_ifp(micro)
    t_solve = time.time() - t0

    dom = micro.domain
    mask = np.asarray(dom.tissue_mask, dtype=bool)
    co = coefficients(micro)
    phi = co["phi"]

    rec = {
        "section": key, "label": spec["label"], "tissue": spec["tissue"],
        "sample": os.path.relpath(spec["sample"], _ROOT),
        "dx_um": float(dx), "grid": [int(dom.shape[0]), int(dom.shape[1])],
        "extent_um": [float(v) for v in dom.extent_um],
        "tissue_fraction": float(mask.mean()),
        "n_cells": int(micro.cells.x.size),
        "t_reconstruct_s": round(t_recon, 2), "t_solve_s": round(t_solve, 2),
    }

    # --- raw field -------------------------------------------------------- #
    inside = ifp[mask]
    rec["raw"] = {
        "max_mmhg": float(ifp.max()),
        "max_in_tissue_mmhg": float(inside.max()),
        "mean_in_tissue_mmhg": float(inside.mean()),
        **{f"{k}_mmhg": v for k, v in _pct(inside).items()},
        "frac_tissue_above_0p5_mmhg": float((inside > 0.5).mean()),
        "frac_tissue_above_2_mmhg": float((inside > 2.0).mean()),
    }

    # --- resolution diagnostic -------------------------------------------- #
    L = co["screen_um"][mask]
    finite = L[np.isfinite(L)]
    rec["resolution"] = {
        "screening_length_um": _pct(finite),
        "frac_tissue_screening_below_dx": float((finite < dx).mean()),
        "frac_tissue_screening_below_2dx": float((finite < 2 * dx).mean()),
        "median_L_over_dx": float(np.median(finite) / dx),
        # Pointwise algebraic limit the discrete operator collapses to when the
        # reaction term dominates the diagonal (react*dx^2 >> K).
        "corr_solved_vs_local_equilibrium": float(
            np.corrcoef(ifp[mask], co["local_equilibrium"][mask])[0, 1]),
        "local_equilibrium_mean_in_tissue_mmhg": float(
            co["local_equilibrium"][mask].mean()),
        "local_equilibrium_max_mmhg": float(co["local_equilibrium"][mask].max()),
        "reaction_dominance_median": float(
            np.median((co["react"][mask] * dx * dx) / np.maximum(co["K"][mask],
                                                                 1e-30))),
    }

    # --- the commensurable statistic: wick-scale local means --------------- #
    rec["wick"] = {}
    for sig_um in WICK_SIGMA_UM:
        cg = coarse_grain(ifp, mask, sig_um / dx)
        entry = {"all_tissue": {**{f"{k}_mmhg": v for k, v in
                                   _pct(cg[mask]).items()},
                                "mean_mmhg": float(cg[mask].mean()),
                                "max_mmhg": float(cg[mask].max())}}
        for thr in TUMOUR_PHI_THRESHOLDS:
            sel = mask & (phi >= thr)
            entry[f"tumour_phi_ge_{thr}"] = {
                "n_nodes": int(sel.sum()),
                "area_frac_of_tissue": float(sel.sum() / max(mask.sum(), 1)),
                **{f"{k}_mmhg": v for k, v in _pct(cg[sel]).items()},
                "mean_mmhg": float(cg[sel].mean()) if sel.any() else None,
                "max_mmhg": float(cg[sel].max()) if sel.any() else None,
            }
        rec["wick"][f"sigma_{sig_um:g}um"] = entry

    rec["tumour_field"] = {
        "phi_max": float(phi.max()),
        "phi_median_in_tissue": float(np.median(phi[mask])),
        **{f"area_frac_phi_ge_{t}": float((mask & (phi >= t)).sum()
                                          / max(mask.sum(), 1))
           for t in TUMOUR_PHI_THRESHOLDS},
    }

    speed = np.hypot(velocity[0], velocity[1])
    rec["velocity_um_s"] = {"max": float(speed.max()),
                            "mean_in_tissue": float(speed[mask].mean())}
    return rec


# --------------------------------------------------------------------------- #
def score(records):
    """Score the production (dx = 12 um) run of each section against the anchor.

    Two bars are reported, deliberately:
      * ``sd_band``  -- mean +/- 1 SD = -0.5 to 4.5 mmHg. This is the band the
        anchor as usually quoted implies. Because IFP is non-negative and the
        solver clips at 0, its effective lower edge is 0, so ANY value in
        0-4.5 mmHg passes. Reported as a weak bar, by construction.
      * ``observed_range`` -- 0.0 to 3.0 mmHg, the interval the paper reports
        10 of its 11 tumours actually lay in. Tighter, and a real test.
    """
    vals = np.array(BOUCHER_TUMOUR_IFP_MMHG, float)
    glioma = np.array([v for v, t in zip(BOUCHER_TUMOUR_IFP_MMHG,
                                         BOUCHER_TUMOUR_TYPE)
                       if "lioblastoma" in t], float)
    out = {
        "anchor": dict(ANCHOR),
        "anchor_recomputed": {
            "table3_mean_mmhg": float(vals.mean()),
            "table3_sd_mmhg": float(vals.std(ddof=1)),
            "table3_median_mmhg": float(np.median(vals)),
            "glioblastoma_only_n": int(glioma.size),
            "glioblastoma_only_mean_mmhg": float(glioma.mean()),
            "glioblastoma_only_range_mmhg": [float(glioma.min()),
                                             float(glioma.max())],
        },
        "bars": {
            "sd_band_mmhg": [ANCHOR["mean_mmhg"] - ANCHOR["spread_mmhg"],
                             ANCHOR["mean_mmhg"] + ANCHOR["spread_mmhg"]],
            "sd_band_effective_mmhg": [0.0, ANCHOR["mean_mmhg"]
                                       + ANCHOR["spread_mmhg"]],
            "observed_range_mmhg": list(ANCHOR["observed_range_10_of_11_mmhg"]),
        },
        "sections": {},
    }
    lo_sd, hi_sd = 0.0, ANCHOR["mean_mmhg"] + ANCHOR["spread_mmhg"]
    lo_obs, hi_obs = ANCHOR["observed_range_10_of_11_mmhg"]
    for key in SECTIONS:
        prod = [r for r in records if r["section"] == key and r["dx_um"] == 12.0]
        if not prod:
            continue
        r = prod[0]
        blk = {}
        for sig in WICK_SIGMA_UM:
            w = r["wick"][f"sigma_{sig:g}um"]["tumour_phi_ge_0.5"]
            med = w["p50_mmhg"]
            blk[f"sigma_{sig:g}um"] = {
                "median_mmhg": med,
                "iqr_mmhg": [w["p25_mmhg"], w["p75_mmhg"]],
                "p95_mmhg": w["p95_mmhg"],
                "max_mmhg": w["max_mmhg"],
                "verdict_sd_band": ("PASS" if med is not None
                                    and lo_sd <= med <= hi_sd else "FAIL"),
                "verdict_observed_range": ("PASS" if med is not None
                                           and lo_obs <= med <= hi_obs
                                           else "FAIL"),
            }
        # grid movement of the headline statistic
        ladder = sorted([x for x in records if x["section"] == key],
                        key=lambda x: -x["dx_um"])
        blk["grid"] = [
            {"dx_um": x["dx_um"],
             "wick250_tumour_median_mmhg":
                 x["wick"]["sigma_250um"]["tumour_phi_ge_0.5"]["p50_mmhg"],
             "wick250_tumour_p95_mmhg":
                 x["wick"]["sigma_250um"]["tumour_phi_ge_0.5"]["p95_mmhg"],
             "raw_mean_in_tissue_mmhg": x["raw"]["mean_in_tissue_mmhg"],
             "raw_max_mmhg": x["raw"]["max_mmhg"],
             "median_L_over_dx": x["resolution"]["median_L_over_dx"]}
            for x in ladder]
        meds = [g["wick250_tumour_median_mmhg"] for g in blk["grid"]
                if g["wick250_tumour_median_mmhg"] is not None]
        if len(meds) > 1 and max(meds) > 0:
            blk["grid_spread_fold"] = float(max(meds) / max(min(meds), 1e-12))
        out["sections"][key] = blk

    # --- the breast section against the BREAST anchor ---------------------- #
    # The discriminating comparison. Unlike the intracranial value, the breast
    # measurement is far from zero and far from normal tissue, so a wrong model
    # cannot pass it by producing nothing.
    if "breast" in out["sections"]:
        med = out["sections"]["breast"]["sigma_250um"]["median_mmhg"]
        lo, hi = (BREAST_ANCHOR["mean_mmhg"] - BREAST_ANCHOR["spread_mmhg"],
                  BREAST_ANCHOR["mean_mmhg"] + BREAST_ANCHOR["spread_mmhg"])
        out["breast_vs_breast_anchor"] = {
            "anchor": dict(BREAST_ANCHOR),
            "corroborating": dict(BREAST_ANCHOR_SECONDARY),
            "bar_mmhg": [lo, hi],
            "predicted_mmhg": med,
            "fold_low_vs_mean": (BREAST_ANCHOR["mean_mmhg"] / med
                                 if med else None),
            "fold_low_vs_band_floor": lo / med if med else None,
            "verdict": "PASS" if (med is not None and lo <= med <= hi) else "FAIL",
        }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dx", type=float, nargs="+", default=list(DX_LADDER))
    ap.add_argument("--sections", nargs="+", default=list(SECTIONS))
    ap.add_argument("--out", default=os.path.join(_ROOT, OUT_REL))
    ap.add_argument("--report-only", action="store_true",
                    help="re-score the existing JSON without re-running solves")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    records = []
    if os.path.exists(args.out):
        try:
            records = json.load(open(args.out)).get("records", [])
        except Exception:
            records = []

    if not args.report_only:
        have = {(r["section"], r["dx_um"]) for r in records}
        for key in args.sections:
            if not os.path.isdir(SECTIONS[key]["sample"]):
                print(f"SKIP {key}: {SECTIONS[key]['sample']} absent")
                continue
            for dx in args.dx:
                if (key, float(dx)) in have:
                    print(f"  cached {key} dx={dx}")
                    continue
                t0 = time.time()
                print(f"  solving {key} dx={dx} ...", flush=True)
                rec = measure_section(key, float(dx))
                records.append(rec)
                # checkpoint after every rung
                json.dump({"records": records}, open(args.out, "w"), indent=1)
                w = rec["wick"]["sigma_250um"]["tumour_phi_ge_0.5"]["p50_mmhg"]
                print(f"    done {time.time() - t0:.1f}s  raw max "
                      f"{rec['raw']['max_mmhg']:.3f}  wick250 tumour median "
                      f"{w:.4f} mmHg", flush=True)

    doc = {"anchor": ANCHOR, "score": score(records), "records": records}
    json.dump(doc, open(args.out, "w"), indent=1)
    print(f"\nwrote {args.out}")
    print(json.dumps(doc["score"]["sections"], indent=1)[:4000])


if __name__ == "__main__":
    main()
