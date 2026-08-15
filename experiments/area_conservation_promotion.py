"""Promote ``geometry='area'`` to the default -- and RE-VERIFY the evidence on
the POST-two-route engine before anything is changed.

WHY THIS EXISTS, AND WHY IT IS NOT A RE-RUN OF ITS TWO PREDECESSORS
-------------------------------------------------------------------
``results/vessel_source_conservation.md`` measured the area-conserving vessel
painting and landed it OPT-IN. ``results/area_conservation_regeneration.md``
measured the blast radius of promoting it and found zero verdict moves. **Both
were measured on the PRE-two-route engine**, at ``alpha =
148.60815765445545``. Commit ``3de97a1`` then flipped two-route transport ON and
moved ``alpha`` to ``22.223477711172563`` -- every drug's compiled diffusivity
changed, so every transport solve behind those nine ratios is a different solve.

**Evidence measured on a superseded engine is not evidence about this one.**
This script re-measures it here, on the current tree, BEFORE the default moves,
so that "zero verdicts move" is a live result rather than an inherited claim.
If a verdict now moves that did not move before, that is the finding and it
stops the promotion.

THE ARMS NAME THEIR GEOMETRY EXPLICITLY, WHICH IS THE WHOLE TRICK
-----------------------------------------------------------------
``experiments/delivery_source_geometry._source_field`` and
``experiments/vessel_source_conservation.source_field`` both build their
``binary`` control by calling ``discrete_vascular_field(domain, vessels)`` with
NO geometry keyword -- i.e. by the default. That is correct while the default is
``'binary'`` and silently WRONG the moment it is not: the control arm would
become the treatment arm and every measurement would read zero change. Every arm
in this file passes ``geometry=`` explicitly, so the same script measures the
same thing before and after the flip. :func:`stage_hygiene` asserts this
directly (``binary_arm_survives_flip``).

STAGES
------
``--stage hygiene``  Every constant this work must not move, checked; the
                     current ``alpha``; and whether the shipped production field
                     is the binary or the area painting. Run BEFORE and AFTER.
``--stage area``     The painted-vessel-area ladder, binary vs area vs subgrid.
                     Pure geometry: nothing here touches transport, so it must
                     reproduce ``vessel_source_conservation.md`` §1-2 exactly.
``--stage fstate``   The ``F_state`` ladder under the three source arms. Pure
                     oxygen: also independent of the two-route flip, so it too
                     must reproduce.
``--stage beta``     The sub-grid closure's ``beta`` on the field the promotion
                     ships, at production ``dx``. This is what makes
                     ``physics/subgrid.py`` staying default-off load-bearing.
``--stage ratios``   THE ONE THAT CAN ACTUALLY MOVE. The nine published rows,
                     re-solved live in both arms, on the post-flip engine.

Run::

    ./.venv/bin/python -u -m experiments.area_conservation_promotion --stage hygiene
    ./.venv/bin/python -u -m experiments.area_conservation_promotion --stage area
    ./.venv/bin/python -u -m experiments.area_conservation_promotion --stage fstate
    ./.venv/bin/python -u -m experiments.area_conservation_promotion --stage beta
    ./.venv/bin/python -u -m experiments.area_conservation_promotion --stage ratios
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

SWEEP = os.path.join(_ROOT, "results", "sweep")
_T0 = time.time()

#: What ``calibrate_alpha()`` returned BEFORE the two-route flip, and therefore
#: the alpha every number in ``results/vessel_source_conservation.md`` and
#: ``results/area_conservation_regeneration.md`` was measured at.
ALPHA_SUPERSEDED = 148.60815765445545

#: What it returns now (``results/two_route_promotion.md``, commit 3de97a1).
#: ``calibrate_alpha`` is a closed-form solve on the compiled anchor drug and the
#: Primeau half-distance; NO vascular field enters it, so this promotion must
#: leave it bit-identical. That is checked, not assumed.
ALPHA_CURRENT = 22.223477711172563

#: The ladder ``experiments/fstate_mechanism.py`` and
#: ``experiments/vessel_source_conservation.py`` use.
DX_LADDER = (30.0, 24.0, 18.0, 16.0, 12.0, 10.0, 8.0, 7.0, 6.0, 5.0)
PRODUCTION_DX_UM = 12.0

ARMS = ("binary", "area", "area_matched")

#: THE LIVE CONTROL, and the reason this constant has to exist.
#:
#: ``delivery_source_geometry.PUBLISHED_RATIO`` holds the ratios as published
#: BEFORE the two-route transport flip (commit ``3de97a1``), i.e. at
#: ``alpha = 148.60815765445545``. That flip changed every drug's compiled
#: diffusivity, so a correct run on the CURRENT engine *must* disagree with it.
#: The harness's ``binary_control_reproduces_published`` flag therefore comes out
#: **False**, and that False is EXPECTED AND BENIGN -- it measures 3de97a1's
#: documentation debt, not this promotion.
#:
#: A red flag that needs a footnote to be understood is a trap, so the real
#: control is named here instead: the six breast rows as re-published at the
#: CURRENT alpha in ``results/two_route_promotion.md`` §5.2 Arm B, which is the
#: shipping record. The ``binary`` arm must reproduce these, and does, to the
#: four decimals they are published at.
#:
#: The three GBM rows are deliberately absent: ``two_route_promotion.md`` scored
#: only the breast panel, so there is no post-flip published GBM value, and
#: inventing one here would manufacture a control rather than apply one.
SHIPPING_RATIO_POST_TWO_ROUTE = {
    ("breast", "doxorubicin"): 0.2330,
    ("breast", "paclitaxel"): 0.1747,
    ("breast", "docetaxel"): 0.1520,
    ("breast", "epirubicin"): 0.2330,
    ("breast", "5-fluorouracil"): 0.3355,
    ("breast", "cisplatin"): 0.3105,
}

#: Tolerance for the above: those values are quoted to four decimals, so this is
#: the tightest check their precision supports.
SHIPPING_RATIO_TOL = 5e-5


def log(msg: str) -> None:
    print(f"[acprom +{time.time() - _T0:7.1f}s] {msg}", flush=True)


def _dump(obj, name: str) -> str:
    os.makedirs(SWEEP, exist_ok=True)
    path = os.path.join(SWEEP, f"acprom_{name}")
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=float)
    log(f"wrote {os.path.relpath(path, _ROOT)}")
    return path


def _span_pct(a: Sequence[float]) -> float:
    x = np.asarray(list(a), dtype=np.float64)
    return float(100.0 * (x.max() - x.min()) / x.mean()) if x.size else float("nan")


# --------------------------------------------------------------------------- #
# The arms -- EVERY ONE NAMES ITS GEOMETRY
# --------------------------------------------------------------------------- #
def source_field(domain, vessels: Dict, arm: str) -> np.ndarray:
    """Paint ``arm``'s source field, naming the geometry explicitly.

    ``experiments/vessel_source_conservation.source_field`` builds its ``binary``
    arm from the DEFAULT. This one does not, which is what lets the identical
    measurement run on both sides of the default flip.
    """
    from tmesim.reconstruct.vessel_network import discrete_vascular_field
    binary = discrete_vascular_field(domain, vessels, geometry="binary")
    if arm == "binary":
        return binary
    if arm == "area":
        return discrete_vascular_field(domain, vessels, geometry="area")
    if arm == "area_matched":
        return discrete_vascular_field(domain, vessels, geometry="area",
                                       match_total_to=binary)
    if arm == "subgrid":
        return discrete_vascular_field(domain, vessels, geometry="subgrid")
    raise ValueError(f"unknown arm {arm!r}")


# --------------------------------------------------------------------------- #
# Stage 0 -- hygiene
# --------------------------------------------------------------------------- #
def stage_hygiene() -> Dict:
    """Every constant this work is forbidden to move, CHECKED rather than claimed."""
    import inspect

    from tmesim.core.config import OxygenConfig
    from tmesim.physics import subgrid as subgrid_mod
    from tmesim.reconstruct import real_fields
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.reconstruct.vessel_network import (
        LINK_DISTANCE_UM, SOURCE_SUBPIXEL_PITCH_UM, VESSEL_RADIUS_UM,
        discrete_vascular_field, segment_vessels)
    from tmesim.validation.calibration import calibrate_alpha
    from tmesim.validation.heldout_patient import DX_UM as _HELDOUT_DX_UM
    from tmesim.validation.mgmt_scalar import DX_UM as _MGMT_DX_UM

    ox = OxygenConfig()
    sig = inspect.signature(discrete_vascular_field)
    # The sub-grid closure is default-off through ONE keyword default, and that
    # default becomes load-bearing under this promotion (it was previously safe
    # twice over: default-off AND a provable no-op on the binary field).
    sg_default = inspect.signature(
        subgrid_mod.solve_oxygen_subgrid).parameters["enabled"].default
    out: Dict = {
        "alpha": float(calibrate_alpha()),
        "alpha_current_expected": ALPHA_CURRENT,
        "alpha_superseded_pre_two_route": ALPHA_SUPERSEDED,
        "oxygen_supply": float(ox.supply),
        "oxygen_consumption": float(ox.consumption),
        "hypoxia_threshold": float(ox.hypoxia_threshold),
        "default_vascular_mode": real_fields._DEFAULT_VASCULAR_MODE,
        "default_geometry": str(sig.parameters["geometry"].default),
        "vessel_radius_um": float(VESSEL_RADIUS_UM),
        "link_distance_um": float(LINK_DISTANCE_UM),
        "subpixel_pitch_um": float(SOURCE_SUBPIXEL_PITCH_UM),
        "dx_um_heldout_patient": float(_HELDOUT_DX_UM),
        "dx_um_mgmt_scalar": float(_MGMT_DX_UM),
        "subgrid_solve_enabled_default": bool(sg_default),
    }
    out["subgrid_default_off_ok"] = bool(sg_default is False)
    out["alpha_ok"] = bool(out["alpha"] == ALPHA_CURRENT)
    out["alpha_moved_by_this_work"] = bool(out["alpha"] != ALPHA_CURRENT)
    out["supply_ok"] = bool(out["oxygen_supply"] == 0.05)
    out["consumption_ok"] = bool(out["oxygen_consumption"] == 0.02)
    out["hypoxia_ok"] = bool(out["hypoxia_threshold"] == 0.2)
    out["dx_um_ok"] = bool(_HELDOUT_DX_UM == 12.0 and _MGMT_DX_UM == 12.0)
    out["mode_ok"] = bool(out["default_vascular_mode"] == "discrete")

    # What does PRODUCTION actually paint right now? The whole promotion is this
    # one bit, so it is measured on the real section rather than read off a
    # signature.
    micro = reconstruct_real(source="xenium", dx_um=PRODUCTION_DX_UM)
    v = segment_vessels(micro.cells, micro.domain)
    shipped = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    f_bin = source_field(micro.domain, v, "binary")
    f_area = source_field(micro.domain, v, "area")
    dx2 = PRODUCTION_DX_UM ** 2
    out.update({
        "shipped_equals_binary": bool(np.array_equal(shipped, f_bin)),
        "shipped_equals_area": bool(np.array_equal(shipped, f_area)),
        "shipped_is_binary_valued": bool(
            np.isin(shipped, (0.0, 1.0)).all()),
        "production_binary_area_um2": float(f_bin.sum() * dx2),
        "production_area_area_um2": float(f_area.sum() * dx2),
        "production_shipped_area_um2": float(shipped.sum() * dx2),
    })
    out["production_total_change_pct"] = float(
        100.0 * (out["production_area_area_um2"]
                 / out["production_binary_area_um2"] - 1.0))
    # The control arm must still be the control arm after the flip.
    out["binary_arm_survives_flip"] = bool(
        out["shipped_equals_binary"] != out["shipped_equals_area"]
        and not np.array_equal(f_bin, f_area))
    for k, val in out.items():
        log(f"  {k} = {val!r}")
    _dump(out, "hygiene.json")
    return out


# --------------------------------------------------------------------------- #
# Stage 1 -- painted area (pure geometry; must reproduce)
# --------------------------------------------------------------------------- #
def stage_area(dx_ladder: Sequence[float] = DX_LADDER) -> Dict:
    """Painted vessel area across the ladder, under each painting.

    Nothing solves a PDE here, so nothing the two-route flip touched can reach
    it. It is re-run anyway because ``vessel_source_conservation.md`` §1-2 is the
    load-bearing claim of the promotion and 'it cannot have changed' is a
    prediction, not a measurement.
    """
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.reconstruct.vessel_network import (
        VESSEL_RADIUS_UM, segment_vessels)

    rows: List[Dict] = []
    for dx in dx_ladder:
        dx = float(dx)
        t0 = time.time()
        micro = reconstruct_real(source="xenium", dx_um=dx)
        dom = micro.domain
        v = segment_vessels(micro.cells, dom)
        tissue = np.asarray(dom.tissue_mask, dtype=bool)
        f_bin = source_field(dom, v, "binary")
        f_area = source_field(dom, v, "area")
        f_sub = source_field(dom, v, "subgrid")
        cell_area = dx * dx
        rows.append({
            "dx_um": dx,
            "grid": list(dom.shape),
            "n_endothelial_cells": int(v["n_endothelial_cells"]),
            "n_vessels": int(v["n_vessels"]),
            "cell_area_um2": cell_area,
            "disc_over_cell_area": float(
                np.pi * VESSEL_RADIUS_UM ** 2 / cell_area),
            "binary_area_um2": float(f_bin.sum() * cell_area),
            "area_area_um2": float(f_area.sum() * cell_area),
            "subgrid_area_um2": float(f_sub.sum() * cell_area),
            "binary_support_frac": float((f_bin > 0)[tissue].mean()),
            "area_support_frac": float((f_area > 0)[tissue].mean()),
            "area_saturated_frac": float((f_area >= 1.0)[tissue].mean()),
            "secs": time.time() - t0,
        })
        log(f"dx={dx:g}: binary {rows[-1]['binary_area_um2']:.1f} | "
            f"area {rows[-1]['area_area_um2']:.1f} | "
            f"subgrid {rows[-1]['subgrid_area_um2']:.0f} um2 "
            f"({rows[-1]['secs']:.0f}s)")
        _dump({"rows": rows, "partial": True}, "area_ladder.json")

    binary = [r["binary_area_um2"] for r in rows]
    area = [r["area_area_um2"] for r in rows]
    prod = next(r for r in rows if r["dx_um"] == PRODUCTION_DX_UM)
    out = {
        "rows": rows,
        "binary_span_pct": _span_pct(binary),
        "area_span_pct": _span_pct(area),
        "subgrid_span_pct": _span_pct([r["subgrid_area_um2"] for r in rows]),
        "production_binary_area_um2": prod["binary_area_um2"],
        "production_area_area_um2": prod["area_area_um2"],
        "production_total_change_pct": float(
            100.0 * (prod["area_area_um2"] / prod["binary_area_um2"] - 1.0)),
        "binary_ladder_mean_um2": float(np.mean(binary)),
        "area_conserved_value_um2": float(np.mean(area)),
        "published_binary_span_pct": 4.7312,
        "published_area_span_pct": 6.3e-14,
        "published_production_total_change_pct": -0.794,
    }
    log(f"SPAN binary {out['binary_span_pct']:.4f} % (published 4.7312 %) | "
        f"area {out['area_span_pct']:.3e} % (published 6.3e-14 %)")
    log(f"production total change {out['production_total_change_pct']:+.3f} % "
        f"(published -0.794 %)")
    _dump(out, "area_ladder.json")
    return out


# --------------------------------------------------------------------------- #
# Stage 2 -- F_state ladder (pure oxygen; must reproduce)
# --------------------------------------------------------------------------- #
def stage_fstate(dx_ladder: Sequence[float] = DX_LADDER,
                 arms: Sequence[str] = ARMS) -> Dict:
    """``F_state`` across the ladder under each source arm.

    Reuses ``experiments.fstate_mechanism``'s own solver / sampler / banding
    helpers verbatim, exactly as ``experiments/vessel_source_conservation.py``
    does, so the ladder is comparable to ``results/fstate_mechanism.md`` rather
    than merely similar to it.
    """
    from experiments.fstate_mechanism import (
        f_state, sample_floor, solve_from_fields)
    from tmesim.core.config import OxygenConfig, PDConfig
    from tmesim.physics.oxygen import _tumor_density
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.reconstruct.vessel_network import segment_vessels

    cfg, pd_cfg = OxygenConfig(), PDConfig()
    rows: List[Dict] = []
    for dx in dx_ladder:
        dx = float(dx)
        t0 = time.time()
        micro = reconstruct_real(source="xenium", dx_um=dx)
        dom = micro.domain
        v = segment_vessels(micro.cells, dom)
        rho = _tumor_density(micro)
        tum = micro.cells.tumor_mask()
        xt, yt = micro.cells.x[tum], micro.cells.y[tum]
        rec: Dict = {"dx_um": dx, "grid": list(dom.shape)}
        for arm in arms:
            f = source_field(dom, v, arm)
            o2 = solve_from_fields(f, rho, float(dom.dx), cfg)
            o2c = sample_floor(o2, float(dom.dx), xt, yt)
            band = f_state(o2c, cfg, pd_cfg)
            rec[arm] = {
                "f_state": float(band["f_state"]),
                "quiescent_frac": float(band["quiescent_frac"]),
                "hypoxic_frac": float(band["hypoxic_frac"]),
                "viable_frac": float(band["viable_frac"]),
                "mean_o2": float(np.mean(o2c)),
                "source_sum": float(f[np.asarray(dom.tissue_mask, bool)].sum()),
                "source_area_um2": float(f.sum() * dx * dx),
            }
        rows.append(rec)
        log(f"dx={dx:g}: " + " | ".join(
            f"{a} F={rec[a]['f_state']:.4f} O2={rec[a]['mean_o2']:.5f}"
            for a in arms) + f"  ({time.time() - t0:.0f}s)")
        _dump({"rows": rows, "partial": True}, "fstate_ladder.json")

    summary: Dict = {}
    for arm in arms:
        fs = [r[arm]["f_state"] for r in rows]
        o2 = [r[arm]["mean_o2"] for r in rows]
        ar = [r[arm]["source_area_um2"] for r in rows]
        summary[arm] = {
            "f_state_span": float(max(fs) - min(fs)),
            "mean_o2_span": float(max(o2) - min(o2)),
            "source_area_span_pct": _span_pct(ar),
            "f_state_at_production": float(
                next(r for r in rows if r["dx_um"] == PRODUCTION_DX_UM)
                [arm]["f_state"]),
        }
    base = summary["binary"]["f_state_span"]
    for arm in arms:
        summary[arm]["drift_removed_pct"] = float(
            100.0 * (1.0 - summary[arm]["f_state_span"] / base))
    out = {"rows": rows, "summary": summary,
           "published_binary_f_state_span": 0.04609,
           "published_area_f_state_span": 0.000546,
           "published_area_drift_removed_pct": 98.82}
    for arm in arms:
        s = summary[arm]
        log(f"{arm:13s} F_state span {s['f_state_span']:.6f} | mean-O2 span "
            f"{s['mean_o2_span']:.6f} | drift removed "
            f"{s['drift_removed_pct']:.2f} %")
    _dump(out, "fstate_ladder.json")
    return out


# --------------------------------------------------------------------------- #
# Stage 3 -- what the closure would do to the field the promotion ships
# --------------------------------------------------------------------------- #
def stage_beta(dx_um: float = PRODUCTION_DX_UM) -> Dict:
    """``beta`` on each painting, at production ``dx``.

    ``results/subgrid_closure.md`` §5.0's control -- "on the shipped field the
    closure is an exact no-op, because ``beta == 1`` identically when the field
    is exactly 0 or 1" -- is a property of the BINARY painting. Promoting the
    area painting voids it, which is precisely why ``physics/subgrid.py`` staying
    default-off stops being a convenience and becomes load-bearing. This stage is
    the measurement that justifies the comment the promotion adds there.
    """
    from experiments.area_conservation_regeneration import _beta_diagnostics
    from tmesim.core.config import OxygenConfig
    from tmesim.physics.oxygen import _tumor_density
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.reconstruct.vessel_network import segment_vessels

    cfg = OxygenConfig()
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
    dom = micro.domain
    v = segment_vessels(micro.cells, dom)
    rho = _tumor_density(micro)
    out: Dict = {"dx_um": float(dx_um), "arms": {}}
    for arm in ("binary", "area"):
        f = source_field(dom, v, arm)
        d = _beta_diagnostics(f, float(dom.dx), float(cfg.supply),
                              float(cfg.consumption) * rho)
        src = f > 0
        d["n_source_cells"] = int(src.sum())
        d["n_partial_cells"] = int(((f > 0) & (f < 1.0)).sum())
        d["partial_frac_of_source"] = (
            float(d["n_partial_cells"] / d["n_source_cells"])
            if d["n_source_cells"] else 0.0)
        out["arms"][arm] = d
        log(f"{arm}: partial cells {d['n_partial_cells']} of "
            f"{d['n_source_cells']} source cells "
            f"({100.0 * d['partial_frac_of_source']:.2f} %); "
            f"beta max {d.get('beta_max')!r}")
    _dump(out, "beta.json")
    return out


# --------------------------------------------------------------------------- #
# Stage 4 -- THE nine published rows, re-solved on the post-flip engine
# --------------------------------------------------------------------------- #
def stage_ratios(dx_um: float = PRODUCTION_DX_UM,
                 tissues: Sequence[str] = ("gbm", "breast"),
                 arm_list: Sequence[str] = ("binary", "area")) -> Dict:
    """The nine published rows, RE-SOLVED under both paintings, at THIS alpha.

    Both arms are solved in the same process off the same reconstruction with
    the same alpha, so the binary/area comparison is internal and exact. The
    comparison against ``delivery_source_geometry.PUBLISHED_RATIO`` is reported
    but is NOT a control here: those ratios were published at ``alpha =
    148.60815765445545``, before commit ``3de97a1``. What must be true is that
    the two arms of THIS run agree on every verdict.

    ``experiments/delivery_source_geometry.py`` is NOT edited (it is owned
    elsewhere); its ``_source_field`` is patched in this process only, and the
    patch names every geometry explicitly so it cannot be fooled by the default.
    """
    import experiments.delivery_source_geometry as dsg
    import tmesim.validation.kill_calibration as kc
    from experiments.area_conservation_regeneration import (
        _norm_family, _ratio_table)
    from experiments.p99_convention import (
        HOMOGENEOUS_TAGS, SHIPPED_TAG, decompose_general)
    from tmesim.reconstruct.vessel_network import LINK_DISTANCE_UM
    from tmesim.validation.calibration import calibrate_alpha

    alpha = float(calibrate_alpha())
    if alpha != ALPHA_CURRENT:
        raise AssertionError(
            f"alpha is {alpha!r}, expected the post-two-route "
            f"{ALPHA_CURRENT!r}. This promotion must not move alpha.")
    log(f"alpha = {alpha!r} (guard OK; superseded value was "
        f"{ALPHA_SUPERSEDED!r})")

    original_source = dsg._source_field
    original_decompose = kc.decompose
    captured: List[Dict] = []

    def patched_source(micro, arm: str) -> np.ndarray:
        """Every arm names its geometry, INCLUDING the binary control."""
        v = dsg._segment_safe(micro.cells, micro.domain, float(LINK_DISTANCE_UM))
        if arm in ("binary", "area", "area_matched", "subgrid"):
            return source_field(micro.domain, v, arm)
        return original_source(micro, arm)

    def capturing_decompose(e_unit, sens, resist, c_bath, ic50, hill,
                            max_kill=1.0, **kw):
        """Record the scoring inputs, then delegate to the SHIPPED function."""
        captured.append({
            "e_unit": np.asarray(e_unit, dtype=np.float64).copy(),
            "sens": np.asarray(sens, dtype=np.float64).copy(),
            "resist": np.asarray(resist, dtype=np.float64).copy(),
            "c_bath": float(c_bath), "ic50": float(ic50), "hill": float(hill),
            "max_kill": float(max_kill),
        })
        return original_decompose(e_unit, sens, resist, c_bath, ic50, hill,
                                  max_kill, **kw)

    ck = os.path.join(SWEEP, f"acprom_ratios_dx{dx_um:g}.jsonl")
    done: Dict[str, Dict] = {}
    if os.path.exists(ck):
        with open(ck) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                p = json.loads(line)
                if (float(p.get("dx_um", -1)) == float(dx_um)
                        and float(p.get("alpha", -1)) == alpha):
                    done[f"{p['tissue']}|{p['arm']}"] = p
        log(f"checkpoint: {len(done)} panel(s) already complete at this alpha")

    dsg._source_field = patched_source
    kc.decompose = capturing_decompose
    panels: List[Dict] = []
    try:
        for tissue in tissues:
            drugs = (dsg.BREAST_DRUGS if tissue == "breast" else dsg.GBM_DRUGS)
            for arm in arm_list:
                key = f"{tissue}|{arm}"
                if key in done:
                    log(f"{key}: cached (normaliser family NOT available)")
                    panels.append(done[key])
                    continue
                log(f"=== panel {key} at dx={dx_um:g} (alpha={alpha:.6f}) ===")
                t0 = time.time()
                captured.clear()
                p = dsg._run_panel(tissue, arm, alpha, drugs, float(dx_um))
                if len(captured) == len(p["rows"]):
                    for row, cap in zip(p["rows"], captured):
                        row["normalisers"] = _norm_family(
                            cap, decompose_general, HOMOGENEOUS_TAGS,
                            SHIPPED_TAG, shipped_ratio=row["ratio_exact"])
                else:  # pragma: no cover - defensive
                    log(f"WARNING: {len(captured)} captures for "
                        f"{len(p['rows'])} rows; normaliser family skipped")
                panels.append(p)
                with open(ck, "a") as fh:
                    fh.write(json.dumps(p, default=float) + "\n")
                log(f"=== panel {key} done in {time.time() - t0:.0f}s ===")
    finally:
        dsg._source_field = original_source
        kc.decompose = original_decompose

    for p in panels:
        n = dsg.rescore_panel(p)
        log(f"{p['tissue']}|{p['arm']}: {n} verdict(s) re-derived through the "
            f"shared band decision")

    out = _ratio_table(panels, dsg, alpha, dx_um, arm_list)
    out["alpha_superseded_pre_two_route"] = ALPHA_SUPERSEDED
    out = _annotate_control(out, alpha)
    _dump(out, f"ratios_dx{dx_um:g}.json")
    return out


def _annotate_control(out: Dict, alpha: float) -> Dict:
    """Make the stale-``PUBLISHED_RATIO`` red flag explain ITSELF.

    ``_ratio_table`` scores the ``binary`` arm against
    ``delivery_source_geometry.PUBLISHED_RATIO``, which is PRE-two-route, so
    ``binary_control_reproduces_published`` is False on a correct run. Left
    alone, that is a red control that a checker has to reconstruct a whole
    conversation to dismiss -- exactly the trap this project keeps falling into.

    So: the misleading key is KEPT (deleting it would hide the staleness rather
    than explain it), renamed in place to say what it actually measures, and the
    genuine current-engine control is computed alongside it and asserted.
    """
    stale_flag = bool(out.get("binary_control_reproduces_published", False))
    out["binary_control_reproduces_PRE_TWO_ROUTE_published"] = stale_flag
    out["binary_control_reproduces_published_EXPECTED"] = False
    out["binary_control_vs_pre_two_route_note"] = (
        "EXPECTED FALSE, AND NOT A CONTROL FAILURE. "
        "delivery_source_geometry.PUBLISHED_RATIO holds the ratios as published "
        "before the two-route transport flip (commit 3de97a1), at alpha = "
        f"{ALPHA_SUPERSEDED!r}. This run is at alpha = {alpha!r}, where every "
        "drug's compiled diffusivity differs, so a CORRECT run must disagree "
        "with it. The disagreement measures 3de97a1's documentation debt, not "
        "the vessel-source promotion. The live control for THIS engine is "
        "binary_control_reproduces_shipping below. PUBLISHED_RATIO was not "
        "updated here: it is duplicated in experiments/delivery_source_"
        "geometry.py, grid_convergence_gbm.py and grid_convergence_breast.py, "
        "all owned elsewhere, and several reports compare against it precisely "
        "BECAUSE it is the historical record."
    )

    rows = []
    worst = 0.0
    for rec in out["table"]:
        key = (rec["tissue"], rec["drug"])
        ship = SHIPPING_RATIO_POST_TWO_ROUTE.get(key)
        if ship is None or "binary" not in rec:
            continue
        got = float(rec["binary"]["ratio"])
        dev = abs(got - float(ship))
        worst = max(worst, dev)
        rows.append({"row": f"{key[0]}|{key[1]}", "shipping_published": ship,
                     "binary_arm": got, "abs_dev": dev})
    reproduces = bool(rows) and worst < SHIPPING_RATIO_TOL
    out["shipping_control_rows"] = rows
    out["binary_control_reproduces_shipping"] = reproduces
    out["worst_shipping_control_deviation"] = float(worst)
    out["binary_control_reproduces_shipping_note"] = (
        "THE LIVE CONTROL. The binary arm, re-solved in this process at this "
        "alpha, against the six breast rows as re-published at the CURRENT "
        "alpha in results/two_route_promotion.md sec 5.2 Arm B. These are "
        "quoted to four decimals, hence the 5e-5 tolerance. The three GBM rows "
        "have no post-two-route published value and are deliberately excluded "
        "rather than invented."
    )

    log("")
    log("=== CONTROL SUMMARY =======================================")
    log(f"vs PRE-two-route PUBLISHED_RATIO : {stale_flag} "
        f"(worst dev {out.get('worst_binary_control_deviation', float('nan')):.2e}) "
        f"-- EXPECTED False, see binary_control_vs_pre_two_route_note")
    log(f"vs CURRENT shipping (two_route_promotion.md 5.2 Arm B) : "
        f"{reproduces} (worst dev {worst:.2e} over {len(rows)} rows)")
    for r in rows:
        log(f"    {r['row']:24s} shipping {r['shipping_published']:.4f} | "
            f"binary {r['binary_arm']:.6f} | dev {r['abs_dev']:.2e}")
    log("===========================================================")
    if not reproduces:
        raise AssertionError(
            "THE LIVE CONTROL FAILED: the binary arm does not reproduce the "
            "current shipping ratios of results/two_route_promotion.md 5.2 Arm "
            f"B (worst deviation {worst:.3e} > {SHIPPING_RATIO_TOL:g}). This is "
            "NOT the benign pre-two-route staleness; something in the source "
            "painting or the transport path has genuinely moved. STOP.")
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", required=True,
                    choices=("hygiene", "area", "fstate", "beta", "ratios"))
    ap.add_argument("--dx", type=float, default=PRODUCTION_DX_UM)
    a = ap.parse_args(argv)
    if a.stage == "hygiene":
        stage_hygiene()
    elif a.stage == "area":
        stage_area()
    elif a.stage == "fstate":
        stage_fstate()
    elif a.stage == "beta":
        stage_beta(dx_um=float(a.dx))
    elif a.stage == "ratios":
        stage_ratios(dx_um=float(a.dx))
    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
