"""IS THE JOINTLY CALIBRATED OXYGEN CONFIGURATION GRID-STABLE?

`results/oxygen_supply_calibration.md` (commit 1c2bf3d) derived, and this file
does NOT re-derive, two things:

  * `supply` is an intravascular mass-transfer coefficient in s^-1 whose
    literature value on this GBM section is ~194 s^-1 (Secomb 2020, PMID
    31423930), against a shipped 0.05;
  * fixing it ALONE makes the oxygen field WORSE (mean tumour pO2 34.7 mmHg vs a
    measured human-GBM median of 5.1 mmHg, Lally 2006 PMID 17207315), because
    the shipped value is a COMPENSATING error -- the sink is also ~140x too
    small (Watabe 2017, PMID 28210996).

Pinning both gives median tumour pO2 4.5-8.2 mmHg with nothing fitted to the
measurement. That report's §6.3 states its own biggest gap verbatim:

    "The jointly calibrated configuration, which sits at F_state 0.42-0.52 with
    a genuine spread of cells across bands, was NOT run across the dx ladder and
    its grid behaviour is UNKNOWN. That is the single most important untested
    question this report leaves open."

This file answers exactly that question and nothing else.

THE SHARP FORM OF THE QUESTION. `F_state`'s lability was CREATED by the
discrete-vessel promotion because mean tumour O2 landed at 0.349 -- a hair
inside the QUIESCENT band edge at 0.4 (`results/oxygen_level_shift.md`). The
jointly calibrated field sits at O = 0.16-0.26, which is next to the OTHER edge
(HYPOXIC/QUIESCENT at 0.2). So the lability may simply MOVE rather than
disappear. The mean is not the right predictor of that; the predictor is how
much of the cell POPULATION sits within grid-noise distance of an edge, which
§3 measures directly.

NOTHING IS CHANGED. `OxygenConfig.supply` is still 0.05, `consumption` still
0.02, the band edges still 0.2/0.4, `REFERENCE_PO2_MMHG` still 35. The only
intervention is `oxygen_supply_calibration._force_cfg`, imported unmodified from
the work being extended: a harness-level wrapper that hands `solve_oxygen`
`dataclasses.replace(cfg, supply=..., consumption=...)`. The solve, the banding
and the whole scoring path are library code.

CFL. `dt = 10 s` is the published value and is held at dx = 30/18/12. Another
workstream is live on whether `dt = 10 s` is admissible for the DISCRETE oxygen
field at dx <= 8, so **dx = 8 is run at dt = 5 s** here (stated in the report)
rather than inheriting that open question. `grid_convergence_gbm.md` §7.3
measured dt = 10 -> 5 at dx = 8 to move the scored ratio by 0.0001 %, so the
rung stays comparable with its published counterpart.

STAGES (each checkpoints to results/sweep/oxjoint_<stage>.json):

  repro     -- dx = 12, all four configurations, oxygen only. ASSERTS the
               jointly calibrated headline numbers of
               `oxygen_supply_calibration.md` §5.1 / §6.1 are reproduced before
               anything else runs. If they are not, that is the finding.
  ladder    -- dx = 30/18/12/8 x {shipped, derived-supply-only, joint lo, joint
               hi}, oxygen only (the O2 solve is steady state, so no dt enters).
               Per-tumour-cell O2 vectors are persisted -- the cell population is
               dx-independent on this section, so they are directly comparable
               rung to rung, which is what makes §3 possible.
  bandedge  -- THE ANSWER. Where the field sits relative to 0.2 and 0.4; the
               per-cell grid-noise scale MEASURED across the ladder; the
               population fraction within that distance of an edge; and the
               count of cells that actually DO change band between rungs.
  verdict   -- reported LAST and separately: the scored 3-D:2-D ratio on the
               ladder under the jointly calibrated arms, through
               `grid_convergence_gbm.py`'s (i.e. `kill_calibration`'s) own
               unmodified scoring path.

Run: ./.venv/bin/python experiments/oxygen_joint_grid_test.py repro ladder bandedge verdict
"""
from __future__ import annotations

import json
import math
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

# Imported unmodified from the work this file extends. `_force_cfg` is the ONLY
# intervention; `measure` is the same reduction `oxygen_supply_calibration.md`
# and `oxygen_level_shift.md` both report, so the numbers are comparable.
from oxygen_supply_calibration import (  # noqa: E402
    BRAIN,
    DERIVED_SUPPLY,
    _force_cfg,
    measure,
)

_T0 = time.time()
SWEEP = os.path.join(_ROOT, "results", "sweep")
O2_STORE = os.path.join(SWEEP, "oxjoint_o2_percell.npz")
VERDICT_CK = os.path.join(SWEEP, "oxjoint_verdict_records.jsonl")

#: Shipped defaults, spelled out so every arm is explicit rather than implicit.
SHIPPED_SUPPLY, SHIPPED_CONSUMPTION = 0.05, 0.02

#: The two jointly calibrated sinks. NOT chosen here: these are the two rungs of
#: `oxygen_supply_calibration.md` §5.1 whose implied CMRO2 (2.227 and 2.748
#: cm3O2/100cm3/min) falls inside Watabe 2017's PET-measured C6-glioma TMRO2
#: band of 1.9-2.8. Nothing was fitted to the pO2 measurement they reproduce.
JOINT_CONSUMPTION_LO, JOINT_CONSUMPTION_HI = 10.0, 20.0

#: (label, supply, consumption). `shipped` reruns the published defaults through
#: the identical harness so any harness artefact would show up as a discrepancy
#: against the published ladder rather than hiding inside the new arms.
ARMS: Tuple[Tuple[str, float, float], ...] = (
    ("shipped", SHIPPED_SUPPLY, SHIPPED_CONSUMPTION),
    ("derived_supply_only", DERIVED_SUPPLY, SHIPPED_CONSUMPTION),
    ("joint_cmro2_lo", DERIVED_SUPPLY, JOINT_CONSUMPTION_LO),
    ("joint_cmro2_hi", DERIVED_SUPPLY, JOINT_CONSUMPTION_HI),
)
JOINT_ARMS = ("joint_cmro2_lo", "joint_cmro2_hi")

#: `grid_convergence_gbm.md`'s rungs, so the two ladders are comparable.
#: 18/12/8 is its constant-ratio (r = 1.5) Richardson triple.
DX_RUNGS: Tuple[float, ...] = (30.0, 18.0, 12.0, 8.0)

#: Band edges, READ from the shipped config rather than written to it.
#: `cells.agents.update_cell_states`: o < t HYPOXIC, [t, 2t) QUIESCENT, >= 2t
#: VIABLE. On the 35 mmHg anchor these are 7 and 14 mmHg.
def band_edges() -> Tuple[float, float]:
    from tmesim.core import OxygenConfig
    t = float(OxygenConfig().hypoxia_threshold)
    return t, 2.0 * t


#: Published dx = 12 headline numbers this file must reproduce before it is
#: allowed to say anything. Source: results/sweep/oxsupply_sinkmatch.json (the
#: raw record behind `oxygen_supply_calibration.md` §5.1) and §6.1's F_state
#: table. Tolerance is tight because the run should be deterministic and
#: bit-comparable -- nothing stochastic sits between the config and these.
REPRO_TARGETS = {
    "joint_cmro2_lo": {
        "mean_o2_tumour": 0.26147588888307727,
        "median_mmhg": 8.193206988360274,
        "f_state": 0.5248530152864102,
        "hypoxic_frac_tumour": 0.4022341676465648,
    },
    "joint_cmro2_hi": {
        "mean_o2_tumour": 0.16135059979584923,
        "median_mmhg": 4.525962348874388,
        "f_state": 0.41674365865949936,
        "hypoxic_frac_tumour": 0.7076264068536872,
    },
    "derived_supply_only": {
        "mean_o2_tumour": 0.9925, "f_state": 1.0000,
    },
    "shipped": {
        "mean_o2_tumour": 0.3486, "f_state": 0.5502,
    },
}
REPRO_RTOL = 1e-6      # applied to the four-decimal targets as 1e-3 instead
REPRO_RTOL_COARSE = 2e-3


def log(msg: str) -> None:
    print(f"[oxjoint {time.time() - _T0:7.1f}s] {msg}", flush=True)


def save(stage: str, payload: dict) -> str:
    os.makedirs(SWEEP, exist_ok=True)
    path = os.path.join(SWEEP, f"oxjoint_{stage}.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=float)
    log(f"checkpoint -> {os.path.relpath(path, _ROOT)}")
    return path


# --------------------------------------------------------------------------- #
# Per-cell O2 store. The whole of §3 depends on these being COMPARABLE across
# rungs, which they are for one specific reason: this Xenium section carries
# measured single-cell coordinates and `reconstruct_real`'s tissue mask is the
# dx-independent 'geometric' one, so refining the mesh changes the continuum
# field and NOT which cells exist (`grid_convergence_gbm.md` §2.1 asserts the
# identical 40,887 / 23,812 counts at every rung). Cell i at dx = 30 is
# therefore the same cell as cell i at dx = 8, and the difference between its
# two O2 values is grid error and nothing else.
# --------------------------------------------------------------------------- #
def load_o2(path: str = O2_STORE) -> Dict[str, np.ndarray]:
    if not os.path.exists(path):
        return {}
    with np.load(path) as z:
        return {k: z[k] for k in z.files}


def store_o2(key: str, o2: np.ndarray, path: str = O2_STORE) -> None:
    data = load_o2(path)
    data[key] = np.asarray(o2, dtype=np.float64)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **data)
    os.replace(tmp, path)


def o2_key(arm: str, dx: float) -> str:
    return f"{arm}|dx{dx:g}"


def _percell_o2(micro) -> np.ndarray:
    """Per-tumour-cell oxygen, sampled exactly as `measure` samples it."""
    from tmesim.core import CellType
    oxy = np.asarray(micro.fields.oxygen, dtype=np.float64)
    tum = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
    col, row = micro.domain.to_index(np.asarray(micro.cells.x)[tum],
                                     np.asarray(micro.cells.y)[tum])
    return oxy[row, col]


def _solve_arm(arm: str, supply: float, consumption: float, dx: float):
    """One oxygen solve. Returns (measure dict, per-cell O2)."""
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    with _force_cfg(supply=supply, consumption=consumption):
        micro = reconstruct_real("xenium", sample=BRAIN, dx_um=dx)
    rec = measure(micro)
    rec.update({"arm": arm, "supply": float(supply),
                "consumption": float(consumption), "dx_um": float(dx)})
    return rec, _percell_o2(micro)


# --------------------------------------------------------------------------- #
# Stage: repro -- the precondition
# --------------------------------------------------------------------------- #
def stage_repro() -> dict:
    """Reproduce `oxygen_supply_calibration.md`'s dx = 12 numbers, or stop.

    This is a precondition, not a result. If the jointly calibrated
    configuration cannot be re-created from that report and its generator, then
    every number downstream is measuring something else and the honest output of
    this file is that failure.
    """
    lo, hi = band_edges()
    out = {"dx_um": 12.0, "band_edges": [lo, hi],
           "reference_po2_mmhg": None, "arms": {}, "checks": []}
    from tmesim.physics.oxygen import REFERENCE_PO2_MMHG
    out["reference_po2_mmhg"] = float(REFERENCE_PO2_MMHG)

    for arm, sup, con in ARMS:
        log(f"repro: {arm} (supply={sup:g}, consumption={con:g}) at dx=12")
        rec, o2 = _solve_arm(arm, sup, con, 12.0)
        rec["median_mmhg"] = rec["tumour_o2_mmhg_percentiles"]["50"]
        out["arms"][arm] = rec
        store_o2(o2_key(arm, 12.0), o2)
        log(f"  meanO2={rec['mean_o2_tumour']:.6f} "
            f"({rec['mean_o2_tumour_mmhg']:.2f} mmHg)  "
            f"median={rec['median_mmhg']:.3f} mmHg  "
            f"H/Q/V {rec['band_fractions']['HYPOXIC']:.4f}/"
            f"{rec['band_fractions']['QUIESCENT']:.4f}/"
            f"{rec['band_fractions']['VIABLE']:.4f}  "
            f"F_state={rec['f_state']:.6f}")

        for field, target in REPRO_TARGETS[arm].items():
            got = float(rec[field])
            tol = REPRO_RTOL if arm in JOINT_ARMS else REPRO_RTOL_COARSE
            rel = abs(got - target) / max(abs(target), 1e-12)
            ok = rel <= tol
            out["checks"].append({"arm": arm, "field": field, "target": target,
                                  "got": got, "rel_error": rel,
                                  "tolerance": tol, "ok": bool(ok)})
            if not ok:
                raise AssertionError(
                    f"REPRODUCTION FAILED: {arm}.{field} = {got!r}, "
                    f"oxygen_supply_calibration.md reports {target!r} "
                    f"(relative error {rel:.3e} > {tol:.0e}). Everything "
                    f"downstream would be measuring a different configuration.")
    log(f"repro: all {len(out['checks'])} checks PASS -- the jointly calibrated "
        f"configuration is reproduced")
    return out


# --------------------------------------------------------------------------- #
# Stage: ladder -- the dx ladder on the jointly calibrated field
# --------------------------------------------------------------------------- #
def stage_ladder() -> dict:
    lo, hi = band_edges()
    out = {"sample": BRAIN, "dx_rungs": list(DX_RUNGS), "band_edges": [lo, hi],
           "arms": {a: {} for a, _, _ in ARMS}, "spans": {}}
    for arm, sup, con in ARMS:
        for dx in DX_RUNGS:
            log(f"ladder: {arm} dx={dx:g}")
            rec, o2 = _solve_arm(arm, sup, con, dx)
            rec["median_mmhg"] = rec["tumour_o2_mmhg_percentiles"]["50"]
            # Distance of the MEAN to each edge, and to the nearer one.
            m = rec["mean_o2_tumour"]
            rec["dist_mean_to_edge_lo"] = float(m - lo)
            rec["dist_mean_to_edge_hi"] = float(m - hi)
            rec["nearest_edge"] = lo if abs(m - lo) <= abs(m - hi) else hi
            rec["dist_mean_to_nearest_edge"] = float(
                min(abs(m - lo), abs(m - hi)))
            out["arms"][arm][f"dx{dx:g}"] = rec
            store_o2(o2_key(arm, dx), o2)
            log(f"  meanO2={rec['mean_o2_tumour']:.6f} F_state={rec['f_state']:.6f} "
                f"H/Q/V {rec['band_fractions']['HYPOXIC']:.4f}/"
                f"{rec['band_fractions']['QUIESCENT']:.4f}/"
                f"{rec['band_fractions']['VIABLE']:.4f}")
            save("ladder_partial", out)

    for arm, _, _ in ARMS:
        vals = [out["arms"][arm][f"dx{dx:g}"]["f_state"] for dx in DX_RUNGS]
        means = [out["arms"][arm][f"dx{dx:g}"]["mean_o2_tumour"] for dx in DX_RUNGS]
        out["spans"][arm] = {
            "f_state_values_coarse_to_fine": vals,
            "f_state_span": float(max(vals) - min(vals)),
            "f_state_span_relative": float(
                (max(vals) - min(vals)) / max(min(vals), 1e-12)),
            "f_state_monotone": _monotone(vals),
            "mean_o2_values_coarse_to_fine": means,
            "mean_o2_span": float(max(means) - min(means)),
            "mean_o2_monotone": _monotone(means),
        }
        log(f"ladder/{arm}: F_state span {out['spans'][arm]['f_state_span']:.4f}, "
            f"monotone={out['spans'][arm]['f_state_monotone']}")
    return out


def _monotone(vals: Sequence[float], tol: float = 0.0) -> bool:
    d = np.diff(np.asarray(vals, dtype=float))
    return bool(np.all(d >= -tol) or np.all(d <= tol))


# --------------------------------------------------------------------------- #
# Stage: bandedge -- the quantity that actually predicts lability
# --------------------------------------------------------------------------- #
def stage_bandedge() -> dict:
    """How much of the POPULATION is within grid-noise distance of a band edge?

    The mean tells you where the distribution sits; it does not tell you how many
    cells a mesh change can relabel. `oxygen_level_shift.md` §6's mechanism is a
    statement about the population near an edge, so that is what is measured.

    THE GRID-NOISE SCALE IS MEASURED, NOT ASSUMED. For each cell i, `sigma_i` is
    the sample standard deviation of its own O2 across the four rungs -- the
    spread the mesh alone puts on that cell's oxygen, with the cell held fixed
    (see the store's docstring for why cell i is the same cell at every rung).
    Two summaries are used: the population median `sigma` (a typical cell) and
    the population mean. A cell is 'edge-labile' if an edge lies within `k *
    sigma_i` of its own O2.

    AND THE DIRECT MEASUREMENT, which needs no noise model at all: the number of
    cells that actually DO change band between one rung and another. That is the
    ground truth `F_state` lability reduces to.
    """
    lo, hi = band_edges()
    store = load_o2()
    out = {"band_edges": [lo, hi], "dx_rungs": list(DX_RUNGS), "arms": {}}
    for arm, _, _ in ARMS:
        keys = [o2_key(arm, dx) for dx in DX_RUNGS]
        missing = [k for k in keys if k not in store]
        if missing:
            log(f"bandedge: SKIP {arm} -- missing {missing}; run `ladder` first")
            continue
        M = np.vstack([store[k] for k in keys])            # (n_rungs, n_cells)
        n_rung, n_cell = M.shape
        sigma = M.std(axis=0, ddof=1)                      # per-cell grid noise
        rec = {
            "n_cells": int(n_cell), "n_rungs": int(n_rung),
            "per_cell_grid_sigma": {
                "median": float(np.median(sigma)),
                "mean": float(sigma.mean()),
                "p90": float(np.percentile(sigma, 90)),
                "max": float(sigma.max()),
                "median_mmhg": float(np.median(sigma) * 35.0),
            },
            "mean_o2_per_rung": [float(r.mean()) for r in M],
        }
        # --- population within k*sigma of an edge, at the published rung ----- #
        ref = store[o2_key(arm, 12.0)]
        for k in (1.0, 2.0, 3.0):
            near_lo = np.abs(ref - lo) <= k * sigma
            near_hi = np.abs(ref - hi) <= k * sigma
            rec[f"frac_within_{k:g}sigma_of_an_edge"] = float(
                (near_lo | near_hi).mean())
            rec[f"frac_within_{k:g}sigma_of_edge_lo"] = float(near_lo.mean())
            rec[f"frac_within_{k:g}sigma_of_edge_hi"] = float(near_hi.mean())
        # Fixed-width alternative, so the answer does not depend on the noise
        # model: population inside a fixed +/- window around either edge.
        for w in (0.005, 0.01, 0.02, 0.05):
            rec[f"frac_within_{w:g}_absolute_of_an_edge"] = float(
                ((np.abs(ref - lo) <= w) | (np.abs(ref - hi) <= w)).mean())
        # --- the direct measurement: cells that actually change band --------- #
        def bands_of(v):
            return np.where(v < lo, 0, np.where(v < hi, 1, 2))
        B = np.vstack([bands_of(store[k]) for k in keys])
        pair = {}
        for i in range(n_rung):
            for j in range(i + 1, n_rung):
                ch = B[i] != B[j]
                pair[f"dx{DX_RUNGS[i]:g}_vs_dx{DX_RUNGS[j]:g}"] = {
                    "n_changed_band": int(ch.sum()),
                    "frac_changed_band": float(ch.mean()),
                    "delta_f_state": float(
                        _f_state_of(B[j]) - _f_state_of(B[i])),
                }
        rec["band_changes_between_rungs"] = pair
        rec["frac_changing_band_somewhere_on_ladder"] = float(
            (B.max(axis=0) != B.min(axis=0)).mean())
        rec["n_changing_band_somewhere_on_ladder"] = int(
            (B.max(axis=0) != B.min(axis=0)).sum())
        # Adjacent-rung worst case -- the number a refinement step can move.
        adj = [pair[f"dx{DX_RUNGS[i]:g}_vs_dx{DX_RUNGS[i+1]:g}"]
               for i in range(n_rung - 1)]
        rec["worst_adjacent_rung_frac_changed"] = float(
            max(a["frac_changed_band"] for a in adj))
        # Where the distribution actually is, relative to the edges.
        rec["quantiles_of_o2"] = {str(p): float(np.percentile(ref, p))
                                  for p in (1, 5, 10, 25, 50, 75, 90, 99)}
        rec["frac_below_edge_lo"] = float((ref < lo).mean())
        rec["frac_between_edges"] = float(((ref >= lo) & (ref < hi)).mean())
        rec["frac_above_edge_hi"] = float((ref >= hi).mean())
        # Local density of cells at each edge: dN/dO per unit O, estimated in a
        # +/-0.01 window. This is the derivative that converts a field shift into
        # a relabelled population, i.e. the mechanism itself.
        for name, e in (("lo", lo), ("hi", hi)):
            w = 0.01
            rec[f"density_at_edge_{name}_per_unit_O"] = float(
                ((np.abs(ref - e) <= w).sum() / (2 * w)) / n_cell)
        out["arms"][arm] = rec
        log(f"bandedge/{arm}: median per-cell grid sigma "
            f"{rec['per_cell_grid_sigma']['median']:.5f} "
            f"({rec['per_cell_grid_sigma']['median_mmhg']:.3f} mmHg); "
            f"{100*rec['frac_within_1sigma_of_an_edge']:.2f} % within 1 sigma "
            f"of an edge; {100*rec['frac_changing_band_somewhere_on_ladder']:.2f} % "
            f"actually change band on the ladder")
    return out


def _f_state_of(bands: np.ndarray) -> float:
    """F_state from a band-index vector (0 HYPOXIC, 1 QUIESCENT, 2 VIABLE).

    Uses the shipped `PDConfig` multipliers, the same way `measure` does, so this
    is the same scalar and not a re-definition of it.
    """
    from tmesim.core import PDConfig
    p = PDConfig()
    mult = np.array([float(p.hypoxia_resistance),
                     float(p.quiescent_resistance), 1.0])
    return float(mult[bands].mean())


# --------------------------------------------------------------------------- #
# Stage: verdict -- reported last, separately, and through the published path
# --------------------------------------------------------------------------- #
#: dt per rung. 10 s is the published value; dx = 8 is run at 5 s because
#: whether dt = 10 s is admissible for the discrete oxygen field at dx <= 8 is
#: an OPEN question owned by another workstream right now, and this ladder must
#: not inherit it. `grid_convergence_gbm.md` §7.3 measured that halving dt at
#: dx = 8 moves the scored ratio by 0.0001 %, so the rung stays comparable.
def dt_for(dx: float) -> float:
    return 5.0 if dx <= 8.0 else 10.0


def stage_verdict() -> dict:
    import grid_convergence_gbm as G
    from tmesim.core import (CellState, Microenvironment, PDConfig, PDEConfig,
                             PKConfig, SimConfig)
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.validation.calibration import _d_eff_tortuous, calibrate_alpha
    from tmesim.validation.kill_calibration import (NOMINAL_HALF_LIFE_S,
                                                    NOMINAL_PEAK_UM, SAVE_EVERY,
                                                    SEED, WINDOW_S, compile_drug,
                                                    prepare_tissue,
                                                    resistance_vector,
                                                    unit_transport,
                                                    window_plasma_mean)
    from tmesim.validation.pk_reality import DEFAULT_ECM

    alpha = float(calibrate_alpha())
    drugs = G.TESTABLE
    ic50_hill, compiled, static = {}, {}, {}
    for name in drugs:
        dp = compile_drug(name, alpha)
        compiled[name] = dp
        ic50_hill[name] = (float(dp.ic50), float(dp.hill))
        static[name] = {
            "ic50_uM": float(dp.ic50), "hill": float(dp.hill),
            "L_um": math.sqrt(_d_eff_tortuous(dp, DEFAULT_ECM)
                              / (float(dp.binding) + float(dp.clearance))),
        }
    bands = G.build_bands(ic50_hill)
    log(f"verdict: alpha = {alpha:.4f} (published, unchanged); bands frozen: "
        + ", ".join(f"{n} [{bands[n]['band']['ratio_lo']:.4f},"
                    f"{bands[n]['band']['ratio_hi']:.4f}]" for n in drugs))

    done: Dict[str, dict] = {}
    if os.path.exists(VERDICT_CK):
        with open(VERDICT_CK) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if {"arm", "dx_um", "drug", "ratio_exact"} <= set(r):
                    done[f"{r['arm']}|{r['dx_um']:g}|{r['drug']}"] = r
    if done:
        log(f"verdict: {len(done)} runs from checkpoint")
    os.makedirs(SWEEP, exist_ok=True)
    ck = open(VERDICT_CK, "a")

    records: List[dict] = []
    plan = [(a, dx) for a in JOINT_ARMS for dx in DX_RUNGS]
    total = len(plan) * len(drugs)
    n = 0
    for arm, dx in plan:
        sup = DERIVED_SUPPLY
        con = (JOINT_CONSUMPTION_LO if arm == "joint_cmro2_lo"
               else JOINT_CONSUMPTION_HI)
        dt = dt_for(dx)
        pending = [d for d in drugs if f"{arm}|{dx:g}|{d}" not in done]
        if not pending:
            for d in drugs:
                records.append(done[f"{arm}|{dx:g}|{d}"])
                n += 1
            log(f"verdict: {arm} dx={dx:g}: all from checkpoint ({n}/{total})")
            continue
        tg = time.time()
        with _force_cfg(supply=sup, consumption=con):
            micro = reconstruct_real(source="xenium", dx_um=dx)
            pde = PDEConfig(dt=dt, total_time=WINDOW_S, convection=True,
                            save_every=SAVE_EVERY)
            pk = PKConfig(plasma_peak=NOMINAL_PEAK_UM,
                          half_life=NOMINAL_HALF_LIFE_S, n_doses=1,
                          dose_interval=3600.0 * 24)
            pd_cfg = PDConfig()
            work, labeled, velocity = prepare_tissue(
                micro, SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
            tumor = labeled.tumor_mask()
            state_t = np.asarray(labeled.state)[tumor]
            plasma_mean = window_plasma_mean(pk, WINDOW_S, dt)
            vx, vy = velocity
            max_speed = float(np.max(np.abs(np.asarray(vx))
                                     + np.abs(np.asarray(vy))))
            cfl = max_speed * dt / dx
            ginfo = {
                "arm": arm, "supply": float(sup), "consumption": float(con),
                "dx_um": float(dx), "dt_s": float(dt),
                "grid": list(micro.domain.shape),
                "n_cells": int(labeled.n), "n_tumor": int(tumor.sum()),
                "state_hypoxic_frac": float(np.mean(
                    state_t == int(CellState.HYPOXIC))),
                "state_quiescent_frac": float(np.mean(
                    state_t == int(CellState.QUIESCENT))),
                "state_viable_frac": float(np.mean(~np.isin(
                    state_t, [int(CellState.HYPOXIC),
                              int(CellState.QUIESCENT)]))),
                "max_speed_um_s": max_speed, "cfl": cfl,
                "prologue_s": time.time() - tg,
            }
            log(f"verdict: {arm} dx={dx:g} dt={dt:g}: grid {micro.domain.shape}, "
                f"H/Q/V {ginfo['state_hypoxic_frac']:.4f}/"
                f"{ginfo['state_quiescent_frac']:.4f}/"
                f"{ginfo['state_viable_frac']:.4f}, CFL={cfl:.3f} "
                f"(prologue {ginfo['prologue_s']:.0f}s)")
            for name in drugs:
                key = f"{arm}|{dx:g}|{name}"
                if key in done:
                    records.append(done[key])
                    n += 1
                    continue
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
                rec = G.score(name, e_unit, sens, resist, drug, pd_cfg, bands,
                              plasma_mean)
                rec.update(ginfo)
                rec.update(static[name])
                rec.update({"drug": name, "solve_s": solve_s,
                            "n_timesteps": int(round(WINDOW_S / dt)),
                            "_e_unit_sum": float(e_unit.sum())})
                records.append(rec)
                ck.write(json.dumps(rec, default=float) + "\n")
                ck.flush()
                n += 1
                log(f"  {n}/{total} {arm} dx={dx:g} {name}: "
                    f"ratio={rec['ratio_exact']:.4f} "
                    f"F_pen={rec['f_penetration']:.4f} "
                    f"F_state={rec['f_state']:.4f} {rec['verdict']} "
                    f"({solve_s:.0f}s)")
        _ = Microenvironment
    ck.close()

    # The mesh really moved: exposure fingerprints must be distinct per arm.
    fp = {}
    for arm in JOINT_ARMS:
        for name in drugs:
            sums = [r["_e_unit_sum"] for r in records
                    if r["arm"] == arm and r["drug"] == name]
            fp[f"{arm}|{name}"] = {"n": len(sums), "n_distinct": len(set(sums))}
            assert len(set(sums)) == len(sums), (
                f"{arm}/{name}: exposure fields NOT distinct across the dx "
                f"ladder -- dx is a no-op")
    out = {"alpha": alpha, "drugs": list(drugs), "dx_rungs": list(DX_RUNGS),
           "dt_by_dx": {f"{dx:g}": dt_for(dx) for dx in DX_RUNGS},
           "bands": {n: bands[n]["band"] for n in drugs},
           "band_is_drug_matched": {n: bands[n]["band_is_drug_matched"]
                                    for n in drugs},
           "exposure_fingerprints": fp, "records": records,
           "summary": {}}
    for arm in JOINT_ARMS:
        for name in drugs:
            rows = sorted([r for r in records
                           if r["arm"] == arm and r["drug"] == name],
                          key=lambda r: -float(r["dx_um"]))
            rr = [float(r["ratio_exact"]) for r in rows]
            lo_b = float(bands[name]["band"]["ratio_lo"])
            hi_b = float(bands[name]["band"]["ratio_hi"])
            out["summary"][f"{arm}|{name}"] = {
                "dx_coarse_to_fine": [float(r["dx_um"]) for r in rows],
                "ratio_coarse_to_fine": rr,
                "f_state_coarse_to_fine": [float(r["f_state"]) for r in rows],
                "f_pen_coarse_to_fine": [float(r["f_penetration"]) for r in rows],
                "verdicts": [r["verdict"] for r in rows],
                "n_pass": sum(1 for r in rows if r["verdict"] == "PASS"),
                "n_rungs": len(rows),
                "span": float(max(rr) - min(rr)) if rr else None,
                "span_relative": (float((max(rr) - min(rr)) / min(rr))
                                  if rr and min(rr) > 0 else None),
                "monotone": _monotone(rr) if rr else None,
                "band": [lo_b, hi_b],
                "margin_to_floor_min": (float(min(rr) - lo_b) if rr else None),
                "margin_to_ceiling_min": (float(hi_b - max(rr)) if rr else None),
                "pct_above_floor_min": (
                    float(100.0 * (min(rr) - lo_b) / lo_b) if rr and lo_b > 0
                    else None),
                "movement_over_margin": (
                    float((max(rr) - min(rr))
                          / min(abs(min(rr) - lo_b), abs(hi_b - max(rr))))
                    if rr else None),
            }
    return out


# --------------------------------------------------------------------------- #
# Stage: fineprobe -- does the jointly calibrated F_state converge AT ALL?
# --------------------------------------------------------------------------- #
#: Extra rungs BELOW the comparable ladder. The oxygen solve is steady state, so
#: no CFL condition and no dt applies to it -- these rungs cost only a linear
#: solve and are not transport runs. They exist because the 30/18/12/8 ladder
#: leaves `F_state` still moving at its finest rung, and "still moving at dx = 8"
#: and "converging to a limit somewhere below dx = 8" are different findings.
DX_FINE: Tuple[float, ...] = (6.0, 5.0, 4.0)


def stage_fineprobe() -> dict:
    out = {"dx_fine": list(DX_FINE), "note": (
        "oxygen-only; the steady-state O2 solve carries no CFL condition and no "
        "dt. These rungs are NOT transport runs and no verdict is taken on "
        "them."), "arms": {}}
    for arm, sup, con in ARMS:
        if arm == "derived_supply_only":
            continue                       # already exactly 1.0000 everywhere
        out["arms"][arm] = {}
        for dx in DX_FINE:
            log(f"fineprobe: {arm} dx={dx:g}")
            t0 = time.time()
            rec, o2 = _solve_arm(arm, sup, con, dx)
            rec["median_mmhg"] = rec["tumour_o2_mmhg_percentiles"]["50"]
            rec["solve_s"] = time.time() - t0
            out["arms"][arm][f"dx{dx:g}"] = rec
            store_o2(o2_key(arm, dx), o2)
            log(f"  meanO2={rec['mean_o2_tumour']:.6f} "
                f"F_state={rec['f_state']:.6f} "
                f"H/Q/V {rec['band_fractions']['HYPOXIC']:.4f}/"
                f"{rec['band_fractions']['QUIESCENT']:.4f}/"
                f"{rec['band_fractions']['VIABLE']:.4f} "
                f"({rec['solve_s']:.0f}s)")
            save("fineprobe_partial", out)
    return out


# --------------------------------------------------------------------------- #
# Stage: mechanism -- WHY the jointly calibrated field is grid-sensitive
# --------------------------------------------------------------------------- #
def stage_mechanism() -> dict:
    """The oxygen boundary layer each configuration has to resolve.

    Same structural argument `grid_convergence_gbm.md` §2 makes for the DRUG
    field, applied to the oxygen field: the mesh has to resolve
    ``L = sqrt(D / k_sink)``, where `k_sink = consumption * rho_tum` is the local
    first-order sink rate. That length is a property of the CONFIGURATION, not of
    the mesh, and it is the thing the shipped and the jointly calibrated
    configurations differ in by two orders of magnitude.
    """
    from tmesim.core import OxygenConfig
    from tmesim.physics.oxygen import _tumor_density
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.reconstruct.vessel_network import VESSEL_RADIUS_UM

    cfg = OxygenConfig()
    D = float(cfg.diffusivity)
    micro = reconstruct_real("xenium", sample=BRAIN, dx_um=12.0)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    rho = np.asarray(_tumor_density(micro), dtype=np.float64)[tissue]
    out = {"diffusivity_um2_s": D, "vessel_radius_um": float(VESSEL_RADIUS_UM),
           "dx_rungs": list(DX_RUNGS),
           "rho_tum_stats": {"mean": float(rho.mean()),
                             "median": float(np.median(rho)),
                             "p90": float(np.percentile(rho, 90)),
                             "p99": float(np.percentile(rho, 99)),
                             "max": float(rho.max())},
           "arms": {}}
    for arm, sup, con in ARMS:
        r = {"supply": float(sup), "consumption": float(con)}
        # Source equilibration length: how far the source acts. sqrt(D/supply).
        r["source_equilibration_length_um"] = float(math.sqrt(D / sup))
        # Sink diffusion length at the tissue-mean and at a dense (p90) tile.
        for tag, rv in (("mean_rho", float(rho.mean())),
                        ("p90_rho", float(np.percentile(rho, 90))),
                        ("p99_rho", float(np.percentile(rho, 99)))):
            k = con * rv
            r[f"sink_rate_s^-1_at_{tag}"] = float(k)
            r[f"sink_diffusion_length_um_at_{tag}"] = (
                float(math.sqrt(D / k)) if k > 0 else None)
        L = math.sqrt(D / (con * float(np.percentile(rho, 90))))
        r["L_at_p90_rho_um"] = float(L)
        r["cells_per_L"] = {f"dx{dx:g}": float(L / dx) for dx in DX_RUNGS}
        out["arms"][arm] = r
        log(f"mechanism/{arm}: source length {r['source_equilibration_length_um']:.1f} um; "
            f"sink length at p90 rho {L:.1f} um = "
            f"{L/12.0:.2f} mesh cells at dx = 12")
    return out


STAGES = {"repro": stage_repro, "ladder": stage_ladder,
          "bandedge": stage_bandedge, "fineprobe": stage_fineprobe,
          "mechanism": stage_mechanism, "verdict": stage_verdict}


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv == ["all"]:
        argv = ["repro", "ladder", "bandedge", "verdict"]
    for stage in argv:
        if stage not in STAGES:
            print(f"unknown stage {stage!r}; expected {sorted(STAGES)}")
            return 2
        log(f"=== STAGE {stage} ===")
        save(stage, STAGES[stage]())
        log(f"=== STAGE {stage} done ===")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
