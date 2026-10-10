"""Evans 2008 (PMID 19138031) DFEOV replication on the engine's own oxygen field.

Pre-registration: ``results/evans2008_perivascular_sign_2026-10-10.md`` sections 1-4,
written and committed BEFORE this script produced a number.

What this does
--------------
Puts the engine's steady-state oxygen field through Evans' *distance from each observed
vessel* algorithm and extracts the one statistic Evans reports that transfers to this
engine without an EF5->pO2 conversion: **the per-vessel sign of the oxygen slope against
radial distance, and its split across vessels.**

It TUNES NOTHING. Every ``OxygenConfig`` value is either the shipped default or the
documented opt-in ``sink_form``. Nothing under ``tmesim/`` is imported for modification.

Run::

    ./.venv/bin/python experiments/evans2008/perivascular_sign.py \
        --sections xenium_319 --dx 8 --sink linear

Writes one JSON record per (section, dx, sink) into
``experiments/evans2008/out/perivascular_sign.json`` (append-safe: a killed run keeps
its finished records).
"""
from __future__ import annotations

import argparse
import json
import os
import time
import traceback
from typing import Dict, List, Optional

import numpy as np

ALPHA_EXPECTED = 22.223477711172563

#: section id -> loader cache dir, from ``tmesim/data/section_catalog.py``.
#: Xenium only: the four Visium sections are excluded by
#: ``app.server.api.SPOT_SPATIAL_REFUSAL`` (synthesised endothelial positions), see
#: pre-registration section 4.2.
XENIUM_SECTIONS: Dict[str, str] = {
    "xenium_319": os.path.join("data_cache", "xenium_gbm"),
    "xenium_io": os.path.join("data_cache", "xenium_gbm_io"),
    "gse324268": os.path.join("data_cache", "gse324268_xenium"),
}

BAND_UM = 10.0
N_BANDS = 5
BAND_CENTRES = np.array([5.0, 15.0, 25.0, 35.0, 45.0])

#: Evans Table 1 pO2 class lower bounds in mmHg, highest first.
EVANS_CLASSES = (
    ("oxygenated", 20.0),
    ("mild_hypoxia", 4.0),
    ("moderate_hypoxia", 0.75),
    ("severe_hypoxia", 0.0),
)

OUT_DIR = os.path.join("experiments", "evans2008", "out")
OUT_JSON = os.path.join(OUT_DIR, "perivascular_sign.json")


def assert_alpha(where: str) -> float:
    from tmesim.validation.calibration import calibrate_alpha

    a = float(calibrate_alpha())
    if repr(a) != repr(ALPHA_EXPECTED):
        raise SystemExit(
            f"[{where}] calibrate_alpha() = {a!r}, expected {ALPHA_EXPECTED!r}. "
            "The engine is not the one this experiment claims to measure; "
            "refusing to produce numbers."
        )
    print(f"[alpha] {where}: calibrate_alpha() = {a!r}  OK", flush=True)
    return a


def evans_class(po2_mmhg: float) -> str:
    for name, lo in EVANS_CLASSES:
        if po2_mmhg >= lo:
            return name
    return "anoxia"


def ols_slope(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Per-row OLS slope of y (n_rows, n_x) against x (n_x,)."""
    xc = x - x.mean()
    return (y * xc).sum(axis=1) / (xc * xc).sum()


def clopper_pearson(k: int, n: int) -> List[float]:
    from scipy.stats import beta

    if n == 0:
        return [float("nan"), float("nan")]
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return [lo, hi]


def band_table(
    oxygen: np.ndarray,
    micro,
    vessels: Dict,
    vessel_radius_um: float,
    min_count: int,
) -> Dict:
    """Evans' DFEOV bands on an arbitrary scalar field, per segmented vessel.

    Returns per-vessel band means, counts, slopes (bands 1-5 and 2-5) and the
    Evans Table 1 class of band 1, for every vessel whose five bands each hold at
    least ``min_count`` tissue grid-cell centres.
    """
    from scipy.spatial import cKDTree

    from tmesim.physics.oxygen import REFERENCE_PO2_MMHG

    dom = micro.domain
    xy = np.asarray(vessels["xy"], dtype=np.float64)
    labels = np.asarray(vessels["labels"], dtype=np.int64)
    n_vessels = int(vessels["n_vessels"])
    if n_vessels == 0 or xy.shape[0] == 0:
        return {"n_vessels": 0, "reason": "no segmented vessels"}

    X, Y = dom.coords_um()
    mask = np.asarray(dom.tissue_mask, dtype=bool)
    pts = np.column_stack([X[mask], Y[mask]])
    vals = np.asarray(oxygen, dtype=np.float64)[mask]

    # Nearest endothelial cell -> owning vessel. This implements Evans' band
    # truncation at neighbouring vessels' contours: every pixel is assessed once,
    # against the vessel it is nearest to.
    d, idx = cKDTree(xy).query(pts, k=1)
    owner = labels[idx]
    r = d - float(vessel_radius_um)          # distance from the painted periphery

    keep = (r >= 0.0) & (r < BAND_UM * N_BANDS)
    owner = owner[keep]
    band = np.minimum((r[keep] / BAND_UM).astype(np.int64), N_BANDS - 1)
    v = vals[keep]

    flat = owner * N_BANDS + band
    size = n_vessels * N_BANDS
    cnt = np.bincount(flat, minlength=size).reshape(n_vessels, N_BANDS)
    tot = np.bincount(flat, weights=v, minlength=size).reshape(n_vessels, N_BANDS)

    ok = (cnt >= int(min_count)).all(axis=1)
    vid = np.flatnonzero(ok)
    if vid.size == 0:
        return {"n_vessels": n_vessels, "n_analysable": 0,
                "reason": f"no vessel has >={min_count} grid centres in all 5 bands"}

    means = tot[vid] / cnt[vid]                       # (n_ok, 5) oxygen saturation
    mmhg = means * float(REFERENCE_PO2_MMHG)

    s15 = ols_slope(BAND_CENTRES, means)
    s25 = ols_slope(BAND_CENTRES[1:], means[:, 1:])

    counts_per_vessel = np.asarray(vessels["cells_per_vessel"], dtype=np.int64)[vid]
    extent = np.asarray(vessels["extent_um"], dtype=np.float64)[vid]
    classes = np.array([evans_class(float(m)) for m in mmhg[:, 0]])

    return {
        "n_vessels": n_vessels,
        "n_analysable": int(vid.size),
        "min_count": int(min_count),
        "vessel_ids": vid,
        "band_means_o2": means,
        "band_means_mmhg": mmhg,
        "band_counts": cnt[vid],
        "slope_15_o2_per_um": s15,
        "slope_25_o2_per_um": s25,
        "slope_15_mmhg_per_um": s15 * float(REFERENCE_PO2_MMHG),
        "classes": classes,
        "cells_per_vessel": counts_per_vessel,
        "extent_um": extent,
    }


def split_summary(slope: np.ndarray, classes: np.ndarray) -> Dict:
    """Evans' Table 4 statistic: counts by slope sign, overall and per class.

    ``expected`` = oxygen FALLS away from the vessel (slope < 0) == Evans'
    POSITIVE EF5 slope. ``anti`` = oxygen RISES away (slope > 0) == Evans'
    NEGATIVE EF5 slope. Exact zeros are counted separately and excluded from the
    proportion, never silently assigned to a side.
    """
    exp_m = slope < 0.0
    anti_m = slope > 0.0
    zero = int((slope == 0.0).sum())
    n_exp, n_anti = int(exp_m.sum()), int(anti_m.sum())
    n = n_exp + n_anti
    out = {
        "n_expected_sign": n_exp,
        "n_anti_gradient": n_anti,
        "n_exact_zero": zero,
        "n_signed": n,
        "p_expected": (n_exp / n) if n else float("nan"),
        "p_expected_cp95": clopper_pearson(n_exp, n) if n else [float("nan")] * 2,
        "per_class": {},
    }
    for name, _ in EVANS_CLASSES:
        sel = classes == name
        if not sel.any():
            continue
        e = int((slope[sel] < 0).sum())
        a = int((slope[sel] > 0).sum())
        out["per_class"][name] = {
            "n_expected_sign": e, "n_anti_gradient": a,
            "p_expected": (e / (e + a)) if (e + a) else float("nan"),
        }
    anox = classes == "anoxia"
    if anox.any():
        out["per_class"]["anoxia"] = {
            "n_expected_sign": int((slope[anox] < 0).sum()),
            "n_anti_gradient": int((slope[anox] > 0).sum()),
        }
    return out


def run_one(section: str, dx_um: float, sink: str, min_count: int) -> Dict:
    from tmesim.core import OxygenConfig
    from tmesim.physics.oxygen import REFERENCE_PO2_MMHG, solve_oxygen, to_mmhg
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.reconstruct.vessel_network import (LINK_DISTANCE_UM,
                                                   VESSEL_RADIUS_UM,
                                                   segment_vessels)

    key = f"{section}@dx{dx_um:g}@{sink}"
    rec: Dict = {
        "section": section, "dx_um": float(dx_um), "sink_form": sink,
        "min_count": int(min_count),
        "vessel_radius_um": float(VESSEL_RADIUS_UM),
        "link_distance_um": float(LINK_DISTANCE_UM),
        "reference_po2_mmhg": float(REFERENCE_PO2_MMHG),
        "band_um": BAND_UM, "n_bands": N_BANDS,
        "band_centres_um": BAND_CENTRES.tolist(),
        "status": "pending",
    }
    t0 = time.time()
    try:
        micro = reconstruct_real(source="xenium", sample=XENIUM_SECTIONS[section],
                                 dx_um=float(dx_um))
        rec["grid"] = list(micro.domain.shape)
        rec["extent_um"] = list(micro.domain.extent_um)
        rec["n_unknowns"] = int(micro.domain.nx * micro.domain.ny)
        rec["tissue_coverage_fraction"] = float(
            np.asarray(micro.domain.tissue_mask, dtype=bool).mean())

        # The field measured must be the ENGINE's own field (prereg 4.1.2).
        o_default = np.asarray(micro.fields.oxygen, dtype=np.float64)
        o_linear = solve_oxygen(micro, OxygenConfig())
        rec["pipeline_field_is_bit_identical_to_solve_oxygen"] = bool(
            np.array_equal(o_default, o_linear))
        if not rec["pipeline_field_is_bit_identical_to_solve_oxygen"]:
            rec["status"] = "refused"
            rec["refusal"] = ("fields.oxygen is not bit-identical to "
                              "solve_oxygen(micro, OxygenConfig()); the field this "
                              "pass would measure is not the engine's own. "
                              "Pre-registered branch R5.")
            return rec

        if sink == "linear":
            O = o_linear
        elif sink == "michaelis_menten":
            O = solve_oxygen(micro, OxygenConfig(sink_form="michaelis_menten"))
        else:
            raise ValueError(f"unknown sink {sink!r}")

        rec["field_o2_min"] = float(O.min())
        rec["field_o2_max"] = float(O.max())
        rec["field_mmhg_min"] = float(to_mmhg(O).min())
        rec["field_mmhg_max"] = float(to_mmhg(O).max())
        # Structural check: can ANY vessel act as an oxygen sink? The source term
        # is supply*vasc*(1-O), so only O > 1 anywhere could make it negative.
        rec["field_o2_exceeds_one_anywhere"] = bool((O > 1.0).any())

        vessels = segment_vessels(micro.cells, micro.domain)
        rec["vessel_segmentation"] = {
            "n_vessels": int(vessels["n_vessels"]),
            "n_endothelial_cells": int(vessels["n_endothelial_cells"]),
            "mvd_per_mm2": float(vessels["mvd_per_mm2"]),
            "singleton_fraction": float(vessels["singleton_fraction"]),
            "n_vessels_ge4_cells": int(vessels["n_vessels_ge4_cells"]),
            "intervessel_nn_median_um": float(vessels["intervessel_nn"]["median"]),
            "tissue_area_mm2": float(vessels["tissue_area_mm2"]),
        }

        bt = band_table(O, micro, vessels, VESSEL_RADIUS_UM, min_count)
        if bt.get("n_analysable", 0) == 0:
            rec["status"] = "void"
            rec["void_reason"] = bt.get("reason", "no analysable vessels")
            rec["n_analysable"] = 0
            return rec

        s15 = bt["slope_15_o2_per_um"]
        s25 = bt["slope_25_o2_per_um"]
        cls = bt["classes"]
        rec["n_analysable"] = int(bt["n_analysable"])
        rec["bands_1_5"] = split_summary(s15, cls)
        rec["bands_2_5"] = split_summary(s25, cls)

        ge4 = bt["cells_per_vessel"] >= 4
        rec["bands_1_5_vessels_ge4_cells"] = (
            split_summary(s15[ge4], cls[ge4]) if ge4.any()
            else {"reason": "no analysable vessel has >=4 endothelial cells"})
        rec["n_analysable_ge4_cells"] = int(ge4.sum())

        mm = bt["band_means_mmhg"]
        rec["perivascular_mmhg"] = {
            "band1_min": float(mm[:, 0].min()), "band1_max": float(mm[:, 0].max()),
            "band1_median": float(np.median(mm[:, 0])),
            "band5_min": float(mm[:, 4].min()), "band5_max": float(mm[:, 4].max()),
            "band5_median": float(np.median(mm[:, 4])),
            "band1_minus_band5_median": float(np.median(mm[:, 0] - mm[:, 4])),
            "band1_minus_band5_max_abs": float(np.abs(mm[:, 0] - mm[:, 4]).max()),
        }
        sl = bt["slope_15_mmhg_per_um"]
        rec["slope_mmhg_per_um"] = {
            "median": float(np.median(sl)), "mean": float(sl.mean()),
            "min": float(sl.min()), "max": float(sl.max()),
            "median_abs": float(np.median(np.abs(sl))),
        }
        rec["class_counts"] = {c: int((cls == c).sum())
                               for c in np.unique(cls).tolist()}
        rec["n_distinct_evans_classes"] = int(len(set(cls.tolist())))
        rec["analysable_vessel_extent_um"] = {
            "median": float(np.median(bt["extent_um"])),
            "max": float(bt["extent_um"].max()),
        }
        rec["analysable_cells_per_vessel"] = {
            "median": float(np.median(bt["cells_per_vessel"])),
            "max": int(bt["cells_per_vessel"].max()),
        }
        rec["status"] = "ok"
    except MemoryError:
        rec["status"] = "infeasible"
        rec["reason"] = "MemoryError during reconstruct/solve at this dx"
    except Exception as exc:                                  # noqa: BLE001
        rec["status"] = "error"
        rec["reason"] = f"{type(exc).__name__}: {exc}"
        rec["traceback"] = traceback.format_exc()[-2000:]
    rec["wall_clock_s"] = round(time.time() - t0, 2)
    print(f"[{key}] {rec['status']}  "
          f"n_analysable={rec.get('n_analysable')}  "
          f"p_expected={rec.get('bands_1_5', {}).get('p_expected')}  "
          f"{rec['wall_clock_s']}s", flush=True)
    return rec


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--sections", default="xenium_319")
    p.add_argument("--dx", default="8")
    p.add_argument("--sink", default="linear")
    p.add_argument("--min-count", type=int, default=3)
    p.add_argument("--out", default=OUT_JSON)
    a = p.parse_args(argv)

    assert_alpha("start")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    store: Dict = {}
    if os.path.exists(a.out):
        with open(a.out) as fh:
            store = json.load(fh)
    store.setdefault("alpha", ALPHA_EXPECTED)
    store.setdefault("runs", {})

    for sec in [s.strip() for s in a.sections.split(",") if s.strip()]:
        for dxs in [s.strip() for s in str(a.dx).split(",") if s.strip()]:
            for sink in [s.strip() for s in a.sink.split(",") if s.strip()]:
                dx = float(dxs)
                key = f"{sec}@dx{dx:g}@{sink}@mc{a.min_count}"
                rec = run_one(sec, dx, sink, a.min_count)
                store["runs"][key] = rec
                with open(a.out, "w") as fh:
                    json.dump(store, fh, indent=1, default=float)

    assert_alpha("end")
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
