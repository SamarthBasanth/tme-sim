"""Resolve the CFL discrepancy between two independent speed measurements.

The question
------------
Two records state the maximum interstitial speed of the SAME nominal
configuration and they straddle the advection stability limit at dx = 8 um:

* ``results/vessel_promotion_provenance.md`` (commit ``88e11d8``, recording the
  state of commit ``ceb96c0``): discrete field, 0.730 / 0.779 / 0.843 / 0.868
  um/s at dx = 15 / 12 / 8 / 6, i.e. ``dt_max = dx/speed`` of
  20.5 / 15.4 / **9.5** / 6.9 s -- so ``dt = 10 s`` VIOLATES the CFL at dx <= 8.
* ``results/grid_convergence_gbm.md`` (commits ``88c7882`` / ``823fcbe``): ran
  its whole ladder at ``dt = 10 s`` and measured 0.7513 um/s at dx = 8
  (CFL 0.939) and 0.6684 at dx = 7 (CFL 0.955).

At dx = 8 those are CFL 0.939 (admissible) versus 1.054 (not). Both cannot
describe the same field. It matters because the GBM ladder's conclusion -- that
the panel is NOT in an asymptotic convergence regime -- rests on its fine rungs.

What this script does
---------------------
``probe``  Re-measures ``max(|vx| + |vy|)`` -- the EXACT quantity
           ``physics/diffusion.py`` builds its CFL number from (line 361) -- on
           the shipped configuration at dx = 15 / 12 / 8 / 7 / 6, and on four
           deliberate variant arms that are the candidate causes of the
           disagreement: the KDE vascular field, the legacy ``'occupancy'``
           tissue mask, ``amplitude='per_vessel_cells'``, and a wider vessel
           link distance. Prologue only -- no transport solve.

``rerun``  If the fine rungs turn out inadmissible at dt = 10 s, re-runs
           dx = 8 and dx = 7 for the three tested drugs at a CFL-satisfying
           dt and re-scores them through ``grid_convergence_gbm.score`` -- the
           published scoring path, imported unmodified -- so the ladder's
           non-monotonicity can be re-read at a stable dt.

Nothing here writes to any file another study owns. Its own outputs are
``results/sweep/cfl_*`` and ``results/cfl_resolution.md``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from tmesim.core import (
    CellState,
    CellType,
    Microenvironment,
    PDConfig,
    PDEConfig,
    PKConfig,
    SimConfig,
)
from tmesim.reconstruct.real_fields import (
    _DEFAULT_BANDWIDTH_UM,
    _DEFAULT_VASCULAR_MODE,
    vascular_density_from_cells,
)
from tmesim.reconstruct.real_pipeline import reconstruct_real
from tmesim.reconstruct.vessel_network import (
    LINK_DISTANCE_UM,
    VESSEL_RADIUS_UM,
    segment_vessels,
)
from tmesim.validation.calibration import calibrate_alpha
from tmesim.validation.kill_calibration import (
    DT_S,
    NOMINAL_HALF_LIFE_S,
    NOMINAL_PEAK_UM,
    SAVE_EVERY,
    SEED,
    WINDOW_S,
    compile_drug,
    prepare_tissue,
    resistance_vector,
    unit_transport,
    window_plasma_mean,
)

import grid_convergence_gbm as G  # published scoring path, imported unmodified

OUT_DIR = os.path.join(_ROOT, "results", "sweep")
PROBE_JSONL = os.path.join(OUT_DIR, "cfl_speed_probe.jsonl")
PROBE_JSON = os.path.join(OUT_DIR, "cfl_speed_probe.json")
RERUN_JSONL = os.path.join(OUT_DIR, "cfl_rerun_records.jsonl")
OUT_MD = os.path.join(_ROOT, "results", "cfl_resolution.md")

DX_PROBE: Tuple[float, ...] = (15.0, 12.0, 8.0, 7.0, 6.0)
DT_CANDIDATES: Tuple[float, ...] = (10.0, 5.0)

_T0 = time.time()


def log(msg: str) -> None:
    print(f"[cflres {time.time()-_T0:7.1f}s] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# Arms -- each is a named, fully-stated configuration
# --------------------------------------------------------------------------- #
#: ``arm -> (tissue_mask_mode, vascular override kwargs or None)``. ``None``
#: means "leave ``reconstruct_real``'s own field untouched", i.e. the shipped
#: default; anything else rebuilds ``fields.vascular_density`` through the same
#: public constructor ``reconstruct_real`` itself calls (line 692), so the only
#: thing that differs is the named keyword.
ARMS: Dict[str, Dict] = {
    "shipped": {"mask": "geometric", "vasc": None,
                "what": "the shipped default: discrete vessels, geometric mask, "
                        "amplitude='unit', link_distance_um=20, radius 7.5 um"},
    "kde": {"mask": "geometric", "vasc": {"mode": "kde"},
            "what": "historical 40 um endothelial KDE vascular field"},
    "occupancy_mask": {"mask": "occupancy", "vasc": None,
                       "what": "shipped discrete field on the LEGACY grid-dependent "
                               "tissue mask"},
    "per_vessel_cells": {"mask": "geometric",
                         "vasc": {"mode": "discrete",
                                  "amplitude": "per_vessel_cells"},
                         "what": "discrete field weighted by endothelial cells per "
                                 "vessel object"},
    "link40": {"mask": "geometric",
               "vasc": {"mode": "discrete", "link_distance_um": 40.0},
               "what": "discrete field with the vessel link distance doubled to "
                       "40 um (fewer, larger vessel objects)"},
    # The state of the reconstruction as it was when the vessel promotion was
    # committed (ceb96c0, 2026-07-28 01:03) -- i.e. BEFORE the classifier arity
    # fix (4ad221f, 2026-07-28 02:42) replaced the raw mean-z marker gate with
    # matched_null_z. ENDOTHELIAL was 9.74 % of cells under the legacy rule and
    # 10.50 % after, so this arm changes WHICH cells are vessels and therefore
    # the IFP source field. Selected through the public
    # ``reconstruct_from_xenium(typing_aggregation=...)`` keyword; nothing in
    # ``tmesim/`` is modified.
    "legacy_typing": {"mask": "geometric", "vasc": None, "typing": "mean_z",
                      "what": "shipped discrete field, but cells typed with the "
                              "LEGACY 'mean_z' marker aggregation (the rule in "
                              "force at the vessel-promotion commit)"},
    "legacy_typing_kde": {"mask": "geometric", "vasc": {"mode": "kde"},
                          "typing": "mean_z",
                          "what": "KDE vascular field under the LEGACY 'mean_z' "
                                  "typing -- the other row of the "
                                  "vessel_promotion.md speed table"},
}


def _cells_census(cells) -> Dict:
    ct = np.asarray(cells.cell_type)
    return {
        "n_cells": int(ct.size),
        "n_tumor": int(np.count_nonzero(cells.tumor_mask())),
        "n_endothelial": int(np.count_nonzero(ct == int(CellType.ENDOTHELIAL))),
    }


def probe_one(dx: float, arm: str) -> Dict:
    """Reconstruct + run the drug-independent prologue; measure the speed field.

    The measured statistic is ``max(|vx| + |vy|)`` over the WHOLE grid, which is
    byte-for-byte the quantity ``physics/diffusion.py`` line 361 computes and
    compares against ``dx/dt``. The Euclidean max speed is reported alongside it
    because the two differ by up to sqrt(2) and confusing them is one candidate
    explanation for the disagreement.
    """
    spec = ARMS[arm]
    t0 = time.time()
    if spec.get("typing"):
        # ``reconstruct_real`` does not expose the typing rule, but the function
        # it delegates to does. Swap the bound name for the duration of this one
        # reconstruction and put it straight back.
        import tmesim.reconstruct.real_pipeline as RP
        _orig = RP.reconstruct_from_xenium
        RP.reconstruct_from_xenium = (
            lambda *a, **k: _orig(*a, **{**k,
                                         "typing_aggregation": spec["typing"]}))
        try:
            micro = reconstruct_real(source="xenium", dx_um=float(dx),
                                     tissue_mask_mode=spec["mask"])
        finally:
            RP.reconstruct_from_xenium = _orig
    else:
        micro = reconstruct_real(source="xenium", dx_um=float(dx),
                                 tissue_mask_mode=spec["mask"])
    vasc_shipped = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    rebuilt_identical = None
    if spec["vasc"] is not None:
        kw = dict(spec["vasc"])
        mode = kw.pop("mode")
        v = vascular_density_from_cells(micro.domain, micro.cells, mode=mode, **kw)
        micro.fields.vascular_density = v
    else:
        # Control: rebuilding through the public constructor with the shipped
        # defaults must reproduce reconstruct_real's own field exactly, or the
        # variant arms below are not comparable with this one.
        v = vascular_density_from_cells(micro.domain, micro.cells)
        rebuilt_identical = bool(np.array_equal(v, vasc_shipped))

    # Segmented for the record only (n_vessel_objects). Uses the same link
    # distance the arm's field was painted with, so the metadata describes the
    # field that was actually solved on rather than the default one.
    _seg_kw = {}
    if (spec["vasc"] or {}).get("link_distance_um") is not None:
        _seg_kw["link_distance_um"] = float(spec["vasc"]["link_distance_um"])
    vessels = segment_vessels(micro.cells, micro.domain, **_seg_kw)
    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True,
                    save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM, half_life=NOMINAL_HALF_LIFE_S,
                  n_doses=1, dose_interval=3600.0 * 24)
    work, labeled, velocity = prepare_tissue(
        micro, SimConfig(pde=pde, pk=pk, pd=PDConfig(), seed=SEED))
    vx = np.asarray(velocity[0], dtype=np.float64)
    vy = np.asarray(velocity[1], dtype=np.float64)
    l1 = float(np.max(np.abs(vx) + np.abs(vy)))
    l2 = float(np.max(np.hypot(vx, vy)))
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    vfield = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    state = np.asarray(labeled.state)
    tumor = labeled.tumor_mask()
    rec = {
        "arm": arm,
        "arm_what": spec["what"],
        "dx_um": float(dx),
        "tissue_mask_mode": spec["mask"],
        "vascular_mode": ("kde" if (spec["vasc"] or {}).get("mode") == "kde"
                          else "discrete"),
        "vascular_amplitude": (spec["vasc"] or {}).get("amplitude", "unit"),
        "link_distance_um": (spec["vasc"] or {}).get("link_distance_um",
                                                     float(LINK_DISTANCE_UM)),
        "vessel_radius_um": float(VESSEL_RADIUS_UM),
        "kde_bandwidth_um": float(_DEFAULT_BANDWIDTH_UM),
        "typing_aggregation": spec.get("typing", "matched_null_z"),
        "grid": [int(s) for s in micro.domain.shape],
        "tissue_px": int(np.count_nonzero(tissue)),
        "rebuilt_field_bit_identical": rebuilt_identical,
        "max_speed_l1_um_s": l1,
        "max_speed_l2_um_s": l2,
        "dt_max_s": (float(dx) / l1) if l1 > 0 else float("inf"),
        "cfl_dt10": l1 * 10.0 / float(dx),
        "cfl_dt5": l1 * 5.0 / float(dx),
        "admissible_dt10": bool(l1 * 10.0 / float(dx) <= 1.0),
        "admissible_dt5": bool(l1 * 5.0 / float(dx) <= 1.0),
        "vasc_mean_tissue": float(vfield[tissue].mean()),
        "vasc_frac_nonzero": float(np.count_nonzero(vfield[tissue]) / max(1, tissue.sum())),
        "n_vessel_objects": int(vessels["n_vessels"]) if "n_vessels" in vessels
                            else int(np.unique(np.asarray(vessels["labels"])).size),
        "ifp_max_mmhg": float(np.max(np.asarray(work.fields.ifp))),
        # CellState.HYPOXIC == 1 (VIABLE 0, HYPOXIC 1, QUIESCENT 2, DEAD 3).
        # Recorded only as a sanity check that the prologue labelled cells at
        # all; nothing in this study is scored on it.
        "hypoxic_frac": float(np.mean(state[tumor] == int(CellState.HYPOXIC)))
                        if tumor.any() else None,
        "prologue_s": time.time() - t0,
    }
    rec.update(_cells_census(labeled))
    log(f"{arm:17s} dx={dx:g} grid={rec['grid']} "
        f"L1max={l1:.4f} um/s  CFL(dt=10)={rec['cfl_dt10']:.3f} "
        f"{'OK' if rec['admissible_dt10'] else 'REFUSES'}  "
        f"CFL(dt=5)={rec['cfl_dt5']:.3f}  endo={rec['n_endothelial']}  "
        f"({rec['prologue_s']:.0f}s)")
    return rec


def _load_jsonl(path: str) -> List[Dict]:
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    out.append(json.loads(line))
                except json.JSONDecodeError:
                    pass
    return out


def probe(dx_list: Sequence[float] = DX_PROBE,
          arms: Sequence[str] = ("shipped",),
          checkpoint: str = PROBE_JSONL) -> List[Dict]:
    os.makedirs(OUT_DIR, exist_ok=True)
    done = {(r["arm"], float(r["dx_um"])) for r in _load_jsonl(checkpoint)}
    if done:
        log(f"checkpoint: {len(done)} probe point(s) already measured")
    with open(checkpoint, "a") as ck:
        for arm in arms:
            for dx in dx_list:
                if (arm, float(dx)) in done:
                    log(f"{arm} dx={dx:g}: from checkpoint")
                    continue
                rec = probe_one(float(dx), arm)
                ck.write(json.dumps(rec, default=float) + "\n")
                ck.flush()
    recs = _load_jsonl(checkpoint)
    with open(PROBE_JSON, "w") as fh:
        json.dump(recs, fh, indent=1, default=float)
    log(f"wrote {PROBE_JSON} ({len(recs)} records)")
    return recs


# --------------------------------------------------------------------------- #
# Re-run: the fine rungs at a second dt, scored through the published path
# --------------------------------------------------------------------------- #
def rerun(dx_list: Sequence[float] = (8.0, 7.0),
          dt_list: Sequence[float] = (5.0,),
          drugs: Sequence[str] = G.TESTABLE,
          checkpoint: str = RERUN_JSONL) -> List[Dict]:
    """Re-run fine rungs at a different dt and score them the published way.

    Uses ``grid_convergence_gbm``'s own ``score()``, ``build_bands()`` and
    ``sensitivity_uniform()`` -- so a ratio produced here is directly comparable
    with the ratios in ``results/grid_convergence_gbm.md`` §3, and any movement
    is attributable to dt alone.
    """
    os.makedirs(OUT_DIR, exist_ok=True)
    done = {(float(r["dx_um"]), float(r["dt_s"]), r["drug"])
            for r in _load_jsonl(checkpoint)}
    alpha = float(calibrate_alpha())
    log(f"alpha = {alpha:.10f} (published, unchanged)")
    compiled = {n: compile_drug(n, alpha) for n in drugs}
    bands = G.build_bands({n: (float(d.ic50), float(d.hill))
                           for n, d in compiled.items()})
    pd_cfg = PDConfig()
    with open(checkpoint, "a") as ck:
        for dx in dx_list:
            dx = float(dx)
            for dt in dt_list:
                dt = float(dt)
                pending = [n for n in drugs if (dx, dt, n) not in done]
                if not pending:
                    log(f"dx={dx:g} dt={dt:g}: all from checkpoint")
                    continue
                t0 = time.time()
                micro = reconstruct_real(source="xenium", dx_um=dx)
                pde = PDEConfig(dt=dt, total_time=WINDOW_S, convection=True,
                                save_every=SAVE_EVERY)
                pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM,
                              half_life=NOMINAL_HALF_LIFE_S, n_doses=1,
                              dose_interval=3600.0 * 24)
                work, labeled, velocity = prepare_tissue(
                    micro, SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
                tumor = labeled.tumor_mask()
                state_t = np.asarray(labeled.state)[tumor]
                plasma_mean = window_plasma_mean(pk, WINDOW_S, dt)
                vx, vy = velocity
                l1 = float(np.max(np.abs(np.asarray(vx)) + np.abs(np.asarray(vy))))
                ginfo = {
                    "dx_um": dx, "dt_s": dt,
                    "grid": [int(s) for s in micro.domain.shape],
                    "n_cells": int(labeled.x.size), "n_tumor": int(tumor.sum()),
                    "max_speed_l1_um_s": l1, "cfl": l1 * dt / dx,
                    "n_timesteps": int(round(WINDOW_S / dt)),
                    "prologue_s": time.time() - t0,
                }
                log(f"dx={dx:g} dt={dt:g}: grid={ginfo['grid']} "
                    f"CFL={ginfo['cfl']:.3f} steps={ginfo['n_timesteps']}")
                for name in pending:
                    drug = compiled[name]
                    resist = resistance_vector(state_t, pd_cfg, drug=drug)
                    sens = G.sensitivity_uniform(work, labeled, name, tumor)
                    ts = time.time()
                    tr = unit_transport(work, drug, pde, velocity)
                    solve_s = time.time() - ts
                    col, row = work.domain.to_index(labeled.x, labeled.y)
                    auc = np.asarray(tr.auc, dtype=np.float64)[row, col]
                    duration = (float(tr.times[-1])
                                if tr.times.size and tr.times[-1] > 0 else 1.0)
                    e_unit = (auc / duration)[tumor]
                    rec = G.score(name, e_unit, sens, resist, drug, pd_cfg,
                                  bands, plasma_mean)
                    rec.update(ginfo)
                    rec.update({"drug": name, "solve_s": solve_s,
                                "_e_unit_sum": float(e_unit.sum()),
                                "_e_unit_first": float(e_unit[0]),
                                "alpha": alpha})
                    ck.write(json.dumps(rec, default=float) + "\n")
                    ck.flush()
                    log(f"  dx={dx:g} dt={dt:g} {name}: "
                        f"ratio={rec['ratio_exact']:.6f} "
                        f"F_pen={rec['f_penetration']:.4f} "
                        f"F_state={rec['f_state']:.4f} {rec['verdict']} "
                        f"({solve_s:.0f}s)")
    return _load_jsonl(checkpoint)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["probe", "rerun"])
    ap.add_argument("--arms", default="shipped")
    ap.add_argument("--dx", default=",".join(f"{d:g}" for d in DX_PROBE))
    ap.add_argument("--dt", default="5")
    ap.add_argument("--drugs", default=",".join(G.TESTABLE),
                    help="comma-separated; one at a time keeps the peak "
                         "resident set to a single sparse factorisation")
    a = ap.parse_args(argv)
    dxs = [float(x) for x in a.dx.split(",") if x.strip()]
    if a.stage == "probe":
        probe(dxs, [s.strip() for s in a.arms.split(",") if s.strip()])
    else:
        rerun(dxs, [float(x) for x in a.dt.split(",") if x.strip()],
              drugs=[s.strip() for s in a.drugs.split(",") if s.strip()])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
