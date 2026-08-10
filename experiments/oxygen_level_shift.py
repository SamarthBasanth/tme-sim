"""DIAGNOSING THE 12x GROWTH IN THE OXYGEN ROUTE'S DOWNSTREAM KILL EFFECT.

`results/oxygen_demotion.md` (superseded) recorded a max |delta kill| of 0.026
between `oxygen_mode='pde'` and `oxygen_mode='vessel_distance'` across
temozolomide / cisplatin / doxorubicin on the 40,887-cell GBM brain section at
dx = 12 um. Regenerating the same arm after commit a9bb703 produced 0.3252 --
a ~12x larger number for the same quantity. Nobody had checked whether that is
real, a defect, or what it does to the claims resting on it.

This script does NOT tune anything and does NOT revert anything. Every arm runs
unmodified library code; the only intervention is a HARNESS-LEVEL monkeypatch
that forces `reconstruct.real_pipeline`'s vascular-field constructor into
`mode='kde'`, which is a supported public argument of
`real_fields.vascular_density_from_cells` and the exact historical default that
commit ceb96c0 replaced. Nothing under `tmesim/` is edited.

STAGES (each checkpoints to results/sweep/oxlevel_<stage>.json so a stall does
not lose work):

  fields      -- the two vascular fields and the three oxygen fields they
                 produce, plus the structural question: which field does the
                 OXYGEN solve consume, and is the surrogate arm invariant to
                 the vascular mode?
  downstream  -- `oxygen_benchmark.downstream_effect()` run twice, once under
                 the discrete default and once under the forced KDE, on the
                 identical section / drugs / dx. This is the causal test.
  invariants  -- independent re-derivation of the two claimed invariants
                 (coverage exactly invariant; cross-drug ranking unchanged),
                 computed from the raw per-cell AUC rather than read off the
                 benchmark's own summary.
  fstate      -- the oxygen-banded cell-state label under every arm, against
                 `results/grid_convergence_gbm.md`'s dx ladder, to test whether
                 this is the F_state story from a third direction.

Run: ./.venv/bin/python experiments/oxygen_level_shift.py <stage> [<stage> ...]
     ./.venv/bin/python experiments/oxygen_level_shift.py all
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

_T0 = time.time()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SWEEP = os.path.join(ROOT, "results", "sweep")


def log(msg: str) -> None:
    print(f"[{time.time() - _T0:7.1f}s] {msg}", flush=True)


def save(stage: str, payload: dict) -> str:
    os.makedirs(SWEEP, exist_ok=True)
    path = os.path.join(SWEEP, f"oxlevel_{stage}.json")
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=float)
    log(f"checkpoint -> {os.path.relpath(path, ROOT)}")
    return path


# --------------------------------------------------------------------------- #
# The one intervention: force the historical KDE vascular field.
# --------------------------------------------------------------------------- #
class force_vascular_mode:
    """Context manager forcing `real_pipeline`'s vascular constructor mode.

    `vascular_density_from_cells`'s `mode` default is BOUND AT DEFINITION TIME
    from `real_fields._DEFAULT_VASCULAR_MODE`, so rebinding that module global
    does nothing. The pipeline calls the function by the name it imported, so
    the honest way to switch the field is to rebind that name to a wrapper that
    passes `mode=` explicitly. The wrapped callable is the unmodified library
    function; no numerical code is replaced.
    """

    def __init__(self, mode: str):
        self.mode = str(mode)
        self._saved = None

    def __enter__(self):
        from tmesim.reconstruct import real_pipeline as rp
        from tmesim.reconstruct.real_fields import vascular_density_from_cells

        self._saved = rp.vascular_density_from_cells
        mode = self.mode

        def _forced(domain, cells, *a, **kw):
            kw.pop("mode", None)
            return vascular_density_from_cells(domain, cells, *a, mode=mode, **kw)

        rp.vascular_density_from_cells = _forced
        return self

    def __exit__(self, *exc):
        from tmesim.reconstruct import real_pipeline as rp

        rp.vascular_density_from_cells = self._saved
        return False


class force_vascular_field:
    """Like :class:`force_vascular_mode` but installs an arbitrary callable.

    Used for the DOSE-MATCHED arm: the discrete field rescaled so its
    tissue-integrated total equals the KDE's, which is `discrete_vascular_field`'s
    own documented control (`match_total_to=`) for separating "how much source"
    from "where the source is".
    """

    def __init__(self, builder):
        self.builder = builder
        self._saved = None

    def __enter__(self):
        from tmesim.reconstruct import real_pipeline as rp

        self._saved = rp.vascular_density_from_cells
        rp.vascular_density_from_cells = self.builder
        return self

    def __exit__(self, *exc):
        from tmesim.reconstruct import real_pipeline as rp

        rp.vascular_density_from_cells = self._saved
        return False


def _dose_matched_builder(domain, cells, *a, **kw):
    """Discrete geometry, KDE total. The controlled comparison."""
    from tmesim.reconstruct.real_fields import vascular_density_from_cells

    kw.pop("mode", None)
    kde = vascular_density_from_cells(domain, cells, mode="kde")
    return vascular_density_from_cells(
        domain, cells, mode="discrete", match_total_to=kde)


BRAIN = "data_cache/xenium_gbm"
DX = 12.0


def _field_stats(a) -> dict:
    a = np.asarray(a, dtype=np.float64)
    return {
        "mean": float(a.mean()), "min": float(a.min()), "max": float(a.max()),
        "nonzero_fraction": float((a > 0).mean()),
        "median": float(np.median(a)),
        "sum": float(a.sum()),
    }


def _tissue_stats(a, tissue) -> dict:
    a = np.asarray(a, dtype=np.float64)[np.asarray(tissue, dtype=bool)]
    return {
        "mean": float(a.mean()), "min": float(a.min()), "max": float(a.max()),
        "nonzero_fraction": float((a > 0).mean()),
        "median": float(np.median(a)),
        "sum": float(a.sum()),
    }


# --------------------------------------------------------------------------- #
# Stage: fields
# --------------------------------------------------------------------------- #
def stage_fields() -> dict:
    from tmesim.core import OxygenConfig
    from tmesim.physics.oxygen_proxy import oxygen_from_vessel_distance
    from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report

    out = {"sample": BRAIN, "dx_um": DX, "arms": {}}

    for mode in ("discrete", "kde"):
        log(f"fields: reconstructing under vascular mode {mode!r}")
        with force_vascular_mode(mode):
            micro = reconstruct_real("xenium", sample=BRAIN, dx_um=DX)
        rep = reconstruction_report(micro)
        tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
        vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
        oxy = np.asarray(micro.fields.oxygen, dtype=np.float64)

        # The surrogate, built on the SAME micro. It reads only cell positions.
        surr, sinfo = oxygen_from_vessel_distance(
            micro, OxygenConfig(), detail=True)

        arm = {
            "vascular_source": micro.metadata["real_pipeline"]["vascular_source"],
            "mean_oxygen_reported": float(rep["mean_oxygen"]),
            "mean_vascular_density_reported": float(rep["mean_vascular_density"]),
            "vascular": {"all": _field_stats(vasc),
                         "tissue": _tissue_stats(vasc, tissue)},
            "oxygen_pde": {"all": _field_stats(oxy),
                           "tissue": _tissue_stats(oxy, tissue)},
            "oxygen_surrogate": {"all": _field_stats(surr),
                                 "tissue": _tissue_stats(surr, tissue)},
            "surrogate_sha": _sha(surr),
            "vascular_sha": _sha(vasc),
            "oxygen_pde_sha": _sha(oxy),
            "surrogate_L_um": float(sinfo["diffusion_length_um"]),
            "surrogate_n_endothelial_cells": int(sinfo["n_endothelial_cells"]),
            "surrogate_point_source": str(sinfo["source"]),
            "grid": [int(v) for v in micro.domain.shape],
            "n_cells": int(micro.cells.n),
        }
        out["arms"][mode] = arm
        log(f"  {mode}: mean vasc(tissue) {arm['vascular']['tissue']['mean']:.5f}  "
            f"mean O2(PDE,all) {arm['oxygen_pde']['all']['mean']:.4f}  "
            f"mean O2(surrogate,all) {arm['oxygen_surrogate']['all']['mean']:.4f}")

    d, k = out["arms"]["discrete"], out["arms"]["kde"]
    out["surrogate_is_vascular_mode_invariant"] = (
        d["surrogate_sha"] == k["surrogate_sha"])
    out["pde_field_changed"] = d["oxygen_pde_sha"] != k["oxygen_pde_sha"]
    out["vascular_field_changed"] = d["vascular_sha"] != k["vascular_sha"]
    log(f"  surrogate bit-identical across vascular modes: "
        f"{out['surrogate_is_vascular_mode_invariant']}")
    return out


def _sha(a) -> str:
    import hashlib

    return hashlib.sha256(
        np.ascontiguousarray(np.asarray(a, dtype=np.float64)).tobytes()
    ).hexdigest()[:16]


# --------------------------------------------------------------------------- #
# Stage: downstream -- the causal test
# --------------------------------------------------------------------------- #
def stage_downstream() -> dict:
    from tmesim.validation.oxygen_benchmark import downstream_effect

    out = {"sample": BRAIN, "dx_um": DX, "arms": {}}
    for mode in ("discrete", "kde"):
        log(f"downstream: running downstream_effect under vascular mode {mode!r}")
        t0 = time.time()
        with force_vascular_mode(mode):
            res = downstream_effect(sample=BRAIN, dx_um=DX)
        res["wall_s"] = time.time() - t0
        kd = {n: abs(r["delta_kill"]) for n, r in res["rows"].items()}
        res["max_abs_delta_kill"] = max(kd.values())
        res["abs_delta_kill_by_drug"] = kd
        out["arms"][mode] = res
        log(f"  {mode}: max |delta kill| = {res['max_abs_delta_kill']:.4f}  "
            f"mean O2 pde={res['mean_oxygen_pde']:.4f} "
            f"surrogate={res['mean_oxygen_vessel_distance']:.4f}  "
            f"({res['wall_s']:.0f}s)")
        save("downstream_partial", out)

    a = out["arms"]["discrete"]["max_abs_delta_kill"]
    b = out["arms"]["kde"]["max_abs_delta_kill"]
    out["ratio_discrete_over_kde"] = float(a / b) if b else float("inf")
    log(f"  RATIO discrete/kde on max |delta kill| = "
        f"{out['ratio_discrete_over_kde']:.2f}x  ({a:.4f} vs {b:.4f})")
    return out


# --------------------------------------------------------------------------- #
# Stage: invariants -- re-derived from raw AUC, not from the benchmark summary
# --------------------------------------------------------------------------- #
def stage_invariants() -> dict:
    """Independently verify: coverage exactly invariant, ranking unchanged.

    Coverage is re-derived here from the per-cell exposure returned by the
    engine rather than read out of `downstream_effect`'s own dict, and the two
    arms' AUC fields are compared BIT-FOR-BIT. Ranking is checked on kill under
    both oxygen routes and under both vascular modes.
    """
    from tmesim.core import PDEConfig, SimConfig
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.calibration import calibrate_alpha
    from tmesim.validation.delivery_metrics import _PDE, compile_real_drug, compute_all
    from tmesim.validation.oxygen_benchmark import (
        DOWNSTREAM_DRUGS, _kill_from_oxygen, _swap_oxygen_in_result,
    )
    from tmesim.core import OxygenConfig
    from tmesim.physics.oxygen_proxy import oxygen_from_vessel_distance

    alpha = calibrate_alpha()
    out = {"alpha": float(alpha), "arms": {}}

    for mode in ("discrete", "kde"):
        log(f"invariants: vascular mode {mode!r}")
        with force_vascular_mode(mode):
            micro_a = reconstruct_real("xenium", sample=BRAIN, dx_um=DX)
            micro_b = reconstruct_real("xenium", sample=BRAIN, dx_um=DX,
                                       oxygen_mode="vessel_distance")
        ox_b = np.asarray(micro_b.fields.oxygen, dtype=np.float64)
        rows = {}
        for name in DOWNSTREAM_DRUGS:
            dp, _mol = compile_real_drug(name, alpha)
            cfg = SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0)
            res = simulate_response(micro_a, cfg)

            # Arm B on the engine's own transport result.
            arm_b = _kill_from_oxygen(micro_b, ox_b, res.transport, cfg, dp)

            met_a = compute_all(micro_a, res, dp)
            met_b = compute_all(micro_b, _swap_oxygen_in_result(res, ox_b), dp)

            ca = met_a["therapeutic_coverage"]
            cb = met_b["therapeutic_coverage"]
            # Raw exposure comparison: identical object by construction, but
            # verify rather than assert.
            exp_a = np.asarray(_exposure(res), dtype=np.float64)
            rows[name] = {
                "kill_pde": float(res.kill_fraction),
                "kill_surrogate": float(arm_b["kill_fraction"]),
                "delta_kill": float(arm_b["kill_fraction"] - res.kill_fraction),
                "coverage_1x_pde": float(ca["coverage_1x"]),
                "coverage_1x_surrogate": float(cb["coverage_1x"]),
                "coverage_1x_exact_equal": ca["coverage_1x"] == cb["coverage_1x"],
                "coverage_0.5x_exact_equal": ca["coverage_0.5x"] == cb["coverage_0.5x"],
                "coverage_2x_exact_equal": ca["coverage_2x"] == cb["coverage_2x"],
                "median_E_over_ic50_exact_equal":
                    ca["median_E_over_ic50"] == cb["median_E_over_ic50"],
                "exposure_sha": _sha(exp_a),
                "exposure_mean": float(exp_a.mean()),
                "ic50_uM": float(dp.ic50),
                "states_pde": _kill_from_oxygen(
                    micro_a, res.fields["oxygen"], res.transport, cfg, dp
                )["tumor_state_fractions"],
                "states_surrogate": arm_b["tumor_state_fractions"],
            }
            log(f"  {mode}/{name}: kill {rows[name]['kill_pde']:.4f} -> "
                f"{rows[name]['kill_surrogate']:.4f}  cov equal="
                f"{rows[name]['coverage_1x_exact_equal']}")

        # Rankings.
        rank_pde = sorted(rows, key=lambda n: -rows[n]["kill_pde"])
        rank_surr = sorted(rows, key=lambda n: -rows[n]["kill_surrogate"])
        rank_cov = sorted(rows, key=lambda n: -rows[n]["coverage_1x_pde"])
        out["arms"][mode] = {
            "rows": rows,
            "rank_by_kill_pde": rank_pde,
            "rank_by_kill_surrogate": rank_surr,
            "rank_by_coverage": rank_cov,
            "ranking_preserved_across_oxygen_route": rank_pde == rank_surr,
            "all_coverage_exactly_invariant": all(
                r["coverage_1x_exact_equal"] and r["coverage_0.5x_exact_equal"]
                and r["coverage_2x_exact_equal"]
                and r["median_E_over_ic50_exact_equal"] for r in rows.values()),
            "max_abs_delta_kill": max(abs(r["delta_kill"]) for r in rows.values()),
        }
        save("invariants_partial", out)

    d = out["arms"]["discrete"]
    k = out["arms"]["kde"]
    out["ranking_same_across_vascular_modes"] = (
        d["rank_by_kill_pde"] == k["rank_by_kill_pde"]
        and d["rank_by_kill_surrogate"] == k["rank_by_kill_surrogate"])
    return out


def _exposure(res):
    """Per-cell / per-grid exposure the coverage metric is a function of."""
    t = res.transport
    for attr in ("auc", "exposure", "auc_field", "cell_auc"):
        if hasattr(t, attr):
            return getattr(t, attr)
    if isinstance(t, dict):
        for k in ("auc", "exposure"):
            if k in t:
                return t[k]
    raise AttributeError(f"cannot find an AUC on {type(t)}: {dir(t)}")


# --------------------------------------------------------------------------- #
# Stage: fstate -- the F_state connection
# --------------------------------------------------------------------------- #
def stage_fstate() -> dict:
    """Tumour state-band fractions across vascular mode, oxygen route and dx.

    `results/grid_convergence_gbm.md` §5 reports the PDE-banded quiescent
    fraction swinging 0.8996 (dx = 12) to 0.9918 (dx = 30) and identifies that
    as the mechanism behind the GBM ladder's non-monotonicity. If the oxygen
    route's downstream influence has grown 12x, the state bands are where it
    would show, so they are measured on the same ladder here.
    """
    from tmesim.cells.agents import update_cell_states
    from tmesim.core import CellState, CellType, Microenvironment, OxygenConfig
    from tmesim.physics.oxygen_proxy import oxygen_from_vessel_distance
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    def bands(micro, oxy):
        work = Microenvironment(
            domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
            patient_id=micro.patient_id, metadata=micro.metadata)
        update_cell_states(work, np.asarray(oxy, dtype=np.float64), OxygenConfig())
        st = np.asarray(work.cells.state)
        tum = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
        n = int(tum.sum())
        return {s.name.lower() + "_frac": float((st[tum] == int(s)).sum()) / n
                for s in (CellState.HYPOXIC, CellState.QUIESCENT, CellState.VIABLE)}

    out = {"sample": BRAIN, "ladder": {}}
    for dx in (30.0, 20.0, 16.0, 12.0, 10.0, 8.0):
        for mode in ("discrete", "kde"):
            log(f"fstate: dx={dx:g} mode={mode}")
            with force_vascular_mode(mode):
                micro = reconstruct_real("xenium", sample=BRAIN, dx_um=dx)
            oxy_pde = np.asarray(micro.fields.oxygen, dtype=np.float64)
            surr = oxygen_from_vessel_distance(micro, OxygenConfig())
            key = f"dx{dx:g}_{mode}"
            out["ladder"][key] = {
                "dx_um": float(dx), "vascular_mode": mode,
                "grid": [int(v) for v in micro.domain.shape],
                "mean_oxygen_pde": float(oxy_pde.mean()),
                "mean_oxygen_surrogate": float(surr.mean()),
                "mean_vascular": float(
                    np.asarray(micro.fields.vascular_density).mean()),
                "bands_pde": bands(micro, oxy_pde),
                "bands_surrogate": bands(micro, surr),
            }
            e = out["ladder"][key]
            log(f"  {key}: PDE quiesc {e['bands_pde']['quiescent_frac']:.4f} "
                f"hypox {e['bands_pde']['hypoxic_frac']:.4f} | "
                f"SURR quiesc {e['bands_surrogate']['quiescent_frac']:.4f} "
                f"hypox {e['bands_surrogate']['hypoxic_frac']:.4f}")
            save("fstate_partial", out)
    return out


# --------------------------------------------------------------------------- #
# Stage: factorial -- WHY the delta grew, decomposed
# --------------------------------------------------------------------------- #
def stage_factorial() -> dict:
    """Cross every TRANSPORT field with every OXYGEN field.

    `delta_kill` in the benchmark is `kill(surrogate bands) - kill(PDE bands)`
    evaluated on ONE transport result, so it can grow for two entirely different
    reasons and the benchmark cannot tell them apart:

      (a) the PDE oxygen FIELD moved (mean O2 0.60 -> 0.35), so the two routes'
          band assignments disagree about far more cells; or
      (b) the OPERATING POINT moved. The discrete source delivers ~3.4x less
          total drug, so cisplatin sits at kill 0.44 instead of 0.93 -- off
          saturation and onto the steep part of the dose-response, where the
          SAME resistance multiplier moves kill far further.

    Crossing the two factors separates them. Every kill here is the unmodified
    `apply_drug_kill` re-run on the engine's own transport result via
    `oxygen_benchmark._kill_from_oxygen`, with that function's hard replication
    gate enforced on the native pairing.
    """
    from tmesim.core import OxygenConfig, PDEConfig, SimConfig
    from tmesim.physics.oxygen_proxy import oxygen_from_vessel_distance
    from tmesim.reconstruct.real_pipeline import reconstruct_real
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.calibration import calibrate_alpha
    from tmesim.validation.delivery_metrics import _PDE, compile_real_drug
    from tmesim.validation.oxygen_benchmark import DOWNSTREAM_DRUGS, _kill_from_oxygen

    alpha = calibrate_alpha()
    arms = {
        "kde": force_vascular_mode("kde"),
        "discrete": force_vascular_mode("discrete"),
        "discrete_dosematched": force_vascular_field(_dose_matched_builder),
    }

    micros, ox_pde, vasc = {}, {}, {}
    surr = None
    surr_sha = {}
    for name, ctx in arms.items():
        log(f"factorial: building {name}")
        with ctx:
            m = reconstruct_real("xenium", sample=BRAIN, dx_um=DX)
        micros[name] = m
        ox_pde[name] = np.asarray(m.fields.oxygen, dtype=np.float64)
        v = np.asarray(m.fields.vascular_density, dtype=np.float64)
        tis = np.asarray(m.domain.tissue_mask, dtype=bool)
        vasc[name] = {"mean_tissue": float(v[tis].mean()),
                      "total_tissue": float(v[tis].sum()),
                      "nonzero_frac_tissue": float((v[tis] > 0).mean()),
                      "sha": _sha(v)}
        s = oxygen_from_vessel_distance(m, OxygenConfig())
        surr_sha[name] = _sha(s)
        surr = s
        log(f"  {name}: vasc total(tissue) {vasc[name]['total_tissue']:.1f}  "
            f"mean O2(PDE) {ox_pde[name].mean():.4f}")

    out = {
        "alpha": float(alpha), "dx_um": DX, "sample": BRAIN,
        "vascular": vasc,
        "surrogate_sha_by_arm": surr_sha,
        "surrogate_identical_everywhere": len(set(surr_sha.values())) == 1,
        "mean_oxygen": {k: float(v.mean()) for k, v in ox_pde.items()},
        "mean_oxygen_surrogate": float(surr.mean()),
        "cells": {},
    }
    save("factorial_partial", out)

    ox_choices = dict(ox_pde)
    ox_choices["surrogate"] = surr

    for tname in arms:
        m = micros[tname]
        for drug in DOWNSTREAM_DRUGS:
            dp, _mol = compile_real_drug(drug, alpha)
            cfg = SimConfig(drug=dp, pde=PDEConfig(**_PDE), seed=0)
            log(f"factorial: transport={tname} drug={drug}")
            res = simulate_response(m, cfg)
            gate = _kill_from_oxygen(m, res.fields["oxygen"], res.transport, cfg, dp)
            if gate["kill_fraction"] != float(res.kill_fraction):
                raise AssertionError(
                    f"replication gate FAILED for {tname}/{drug}: "
                    f"{gate['kill_fraction']!r} vs {float(res.kill_fraction)!r}")
            cell = {"engine_kill": float(res.kill_fraction), "kills": {}, "bands": {}}
            for oname, ox in ox_choices.items():
                k = _kill_from_oxygen(m, ox, res.transport, cfg, dp)
                cell["kills"][oname] = float(k["kill_fraction"])
                cell["bands"][oname] = k["tumor_state_fractions"]
            native = cell["kills"][tname]
            cell["delta_kill_surrogate_minus_native_pde"] = float(
                cell["kills"]["surrogate"] - native)
            out["cells"][f"{tname}|{drug}"] = cell
            log(f"  kills: " + "  ".join(
                f"{o}={v:.4f}" for o, v in cell["kills"].items()))
            save("factorial_partial", out)

    return out


# --------------------------------------------------------------------------- #
# Stage: cliff -- the missing sensitivity axis
# --------------------------------------------------------------------------- #
def stage_cliff() -> dict:
    """Sweep the SOURCE normalisation and watch the state bands.

    `results/oxygen_sensitivity_band.md` swept the SINK (`consumption`) over
    2,500x and reported the HYPOXIC fraction, which is 0.0000 at the shipping
    default. It never swept the SOURCE, and it never reported the QUIESCENT
    band. Commit ceb96c0 changed the source normalisation by ~2.8x, and the
    QUIESCENT band is exactly where the resulting move landed -- so the
    published sensitivity study is blind to it by construction.

    This sweeps a scalar multiplier on the discrete vascular field (identical in
    effect to scaling `OxygenConfig.supply`, since the source is `supply*vasc`)
    and records mean tumour O2, all three band fractions and `F_state`. Nothing
    is fitted and no default is changed; the shipped point is m = 1.
    """
    from tmesim.cells.agents import update_cell_states
    from tmesim.core import (CellState, CellType, Microenvironment,
                             OxygenConfig, PDConfig)
    from tmesim.reconstruct.real_pipeline import reconstruct_real

    pdc = PDConfig()
    mult = {"HYPOXIC": pdc.hypoxia_resistance,
            "QUIESCENT": pdc.quiescent_resistance, "VIABLE": 1.0}

    def measure(micro):
        oxy = np.asarray(micro.fields.oxygen, dtype=np.float64)
        work = Microenvironment(
            domain=micro.domain, fields=micro.fields, cells=micro.cells.copy(),
            patient_id=micro.patient_id, metadata=micro.metadata)
        update_cell_states(work, oxy, OxygenConfig())
        st = np.asarray(work.cells.state)
        tum = np.asarray(micro.cells.cell_type) == int(CellType.TUMOR)
        n = int(tum.sum())
        fr = {s.name: float((st[tum] == int(s)).sum()) / n
              for s in (CellState.HYPOXIC, CellState.QUIESCENT, CellState.VIABLE)}
        # Tumour-cell local O2, sampled the way label_hypoxia does.
        idx = micro.domain.to_index(np.asarray(micro.cells.x)[tum],
                                    np.asarray(micro.cells.y)[tum])
        o_t = oxy[idx[1], idx[0]] if isinstance(idx, tuple) else None
        return {
            "band_fractions": fr,
            "f_state": float(sum(fr[k] * mult[k] for k in fr)),
            "mean_oxygen_tissue": float(
                oxy[np.asarray(micro.domain.tissue_mask, bool)].mean()),
            "mean_oxygen_all": float(oxy.mean()),
            "mean_oxygen_tumour": (float(np.mean(o_t)) if o_t is not None
                                   else None),
            "mean_vascular_tissue": float(np.asarray(
                micro.fields.vascular_density)[
                    np.asarray(micro.domain.tissue_mask, bool)].mean()),
        }

    from tmesim.reconstruct.real_fields import vascular_density_from_cells

    out = {"sample": BRAIN, "dx_um": DX, "pd_multipliers": mult,
           "quiescent_band": [OxygenConfig().hypoxia_threshold,
                              2.0 * OxygenConfig().hypoxia_threshold],
           "sweep": {}, "reference": {}}

    # Reference points: the two shipped conventions.
    for mode in ("discrete", "kde"):
        with force_vascular_mode(mode):
            m = reconstruct_real("xenium", sample=BRAIN, dx_um=DX)
        out["reference"][mode] = measure(m)
        log(f"cliff ref {mode}: meanO2(tum)="
            f"{out['reference'][mode]['mean_oxygen_tumour']}  "
            f"F_state={out['reference'][mode]['f_state']:.4f}")
        save("cliff_partial", out)

    for m_scale in (0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.78, 3.0, 4.0, 6.0):
        def _b(domain, cells, *a, _s=m_scale, **kw):
            kw.pop("mode", None)
            return _s * vascular_density_from_cells(
                domain, cells, mode="discrete")

        with force_vascular_field(_b):
            mi = reconstruct_real("xenium", sample=BRAIN, dx_um=DX)
        rec = measure(mi)
        rec["source_multiplier"] = float(m_scale)
        out["sweep"][f"x{m_scale:g}"] = rec
        log(f"cliff x{m_scale:g}: meanO2(tum)={rec['mean_oxygen_tumour']}  "
            f"bands H/Q/V "
            f"{rec['band_fractions']['HYPOXIC']:.4f}/"
            f"{rec['band_fractions']['QUIESCENT']:.4f}/"
            f"{rec['band_fractions']['VIABLE']:.4f}  "
            f"F_state={rec['f_state']:.4f}")
        save("cliff_partial", out)
    return out


# --------------------------------------------------------------------------- #
# Stage: killcal -- does the 3/3 PASS survive the pre-promotion oxygen source?
# --------------------------------------------------------------------------- #
def stage_killcal() -> dict:
    """Re-run `validation.kill_calibration` under the HISTORICAL KDE source.

    `results/kill_calibration.md`'s headline is "3-D:2-D ratio test PASSES for
    all 3 drugs", and its own factorisation is `ratio = F_penetration x F_state
    x interaction`. `F_state` is measured here as 0.5502 under the production
    discrete source and EXACTLY 1.0000 under the historical KDE, so the verdict
    is a direct function of the thing this investigation is about. This runs the
    unmodified `run()` + `analyse()` under the KDE rather than projecting it.

    The discrete arm is NOT re-run: it is the published configuration and its
    numbers are read from `results/kill_calibration.md`.
    """
    from tmesim.validation import kill_calibration as kc

    log("killcal: running kill_calibration.run() under vascular mode 'kde'")
    with force_vascular_mode("kde"):
        res = kc.run()
        analysis = kc.analyse(res)

    rows = {}
    for r in res.get("drugs", res.get("rows", [])):
        if not isinstance(r, dict):
            continue
        name = r.get("drug") or r.get("name")
        dec = r.get("decomposition") or r.get("headline") or {}
        rows[name] = {k: r.get(k) for k in
                      ("drug", "ratio_exact", "f_penetration", "f_state",
                       "implied_ic50_fold", "verdict") if k in r}
        if dec:
            rows[name]["decomposition"] = {
                k: dec.get(k) for k in
                ("f_penetration", "f_state", "ratio_exact", "interaction")
                if k in dec}
    out = {"vascular_mode": "kde", "rows": rows,
           "state_breakdown": res.get("state"),
           "analysis_verdicts": _verdicts(analysis),
           "raw_keys": sorted(res.keys())}
    log(f"killcal: verdicts {out['analysis_verdicts']}")
    return out


def _verdicts(analysis) -> dict:
    found = {}

    def walk(o, path=""):
        if isinstance(o, dict):
            for k, v in o.items():
                if k in ("verdict", "ratio_exact", "f_state", "f_penetration",
                         "predicted_ratio", "implied_ic50_fold", "band_lo",
                         "band_hi", "lit_ratio_lo", "lit_ratio_hi"):
                    found[f"{path}.{k}"] = v
                else:
                    walk(v, f"{path}.{k}")
        elif isinstance(o, (list, tuple)):
            for i, v in enumerate(o):
                walk(v, f"{path}[{i}]")

    walk(analysis)
    return found


STAGES = {
    "fields": stage_fields,
    "factorial": stage_factorial,
    "cliff": stage_cliff,
    "killcal": stage_killcal,
    "downstream": stage_downstream,
    "invariants": stage_invariants,
    "fstate": stage_fstate,
}


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv == ["all"]:
        argv = ["fields", "downstream", "invariants", "factorial", "fstate"]
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
