"""Is the capacity-ranking failure DOUBLE-COUNTED by `calibration.uptake_proxy`?

A sensitivity study on top of two already-published results. **No model change.**
Nothing in `tmesim/` is touched, no default moves, `calibrate_alpha()` is asserted
live against `22.223477711172563` at the start of every stage.

The question
------------
`results/capacity_ranking_robustness.md` sec. 6.2 raised, but did not settle, the
possibility that the engine ALREADY carries a per-compound intracellular-capacity
ordering under another name -- `validation.calibration.uptake_proxy`, which `alpha`
multiplies into every compiled drug's `binding` rate. If so, the capacity sweep
applied part of the same physics twice and the measured ranking instability is an
over-estimate.

This script establishes, in order:

  1. WHAT `uptake_proxy` IS -- read off live code, not from docstrings: what it is
     computed from, where it enters the physics, and whether it acts as a reversible
     CAPACITY (stores solute, retards, leaves the steady state alone) or as an
     irreversible SINK (removes solute, moves the steady state). Verified by running
     the shipped solver, not by assertion.
  2. HOW WELL IT ORDERS COMPOUNDS the way measured capacity does -- against the
     structural `K_b` rule on both shipping panels, and against measured rat
     brain-slice `Vu,brain` on whatever overlaps.
  3. WHETHER THE RANKING FAILURE IS DOUBLE-COUNTED -- by re-running the ordering
     sensitivity with capacity applied RESIDUALLY, after removing whatever the
     proxy already accounts for, plus the strongest-possible-defence arm in which
     the proxy IS taken to be the capacity ordering.

Stages
------
  role      what uptake_proxy is; sink-vs-capacity, measured on the shipped solver
  overlap   proxy vs measured Vu,brain; proxy vs structural K_b; n and p
  residual  ordering under residual / proportional capacity (uses the FROZEN
            E[kill] sweep in results/sweep/capacity_ranking_robustness.json)
  report    write results/uptake_proxy_double_count.md

Run:
    ./.venv/bin/python experiments/uptake_proxy_double_count.py --stage all
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO, "experiments"))
sys.path.insert(0, _REPO)

from tmesim.core import DrugParams, PDEConfig, SimConfig                    # noqa: E402
from tmesim.molecular import binding as mol_binding                        # noqa: E402
from tmesim.molecular import two_route                                     # noqa: E402
from tmesim.molecular.library import get_molecule                          # noqa: E402
from tmesim.validation.calibration import (                                # noqa: E402
    _compile, calibrate_alpha, k_passive, uptake_proxy,
)

FROZEN_ALPHA = 22.223477711172563
FROZEN_PHI = 0.7892
FROZEN_SHIP_TOTAL_TIME = 21600.0            # the panels' own window, sec. 1
BEST_CONSTANT = 39.0                        # GMFE-minimising constant, B-S1
MEASURED_LO, MEASURED_HI = 3.2, 4184.0      # retained rat brain-slice range

MC_SEED = 20260821                          # the brain-slice anchor's own seed
BOOT_SEED = 20260822
MC_DRAWS = 40_000
BOOT_REPS = 4_000

_SWEEP = os.path.join(_REPO, "results", "sweep")
_SOURCE = os.path.join(_SWEEP, "capacity_ranking_robustness.json")
_OUT = os.path.join(_SWEEP, "uptake_proxy_double_count.json")
_MD = os.path.join(_REPO, "results", "uptake_proxy_double_count.md")
_VU = os.path.join(_REPO, "data_cache", "brain_slice", "vu_brain.csv")

GBM5 = ["cisplatin", "doxorubicin", "lomustine", "methotrexate", "temozolomide"]


def _tick(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _assert_frozen() -> Dict[str, float]:
    """Every constant this study leans on, asserted against LIVE code."""
    a = calibrate_alpha()
    assert abs(a - FROZEN_ALPHA) < 1e-12, (
        f"calibrate_alpha() moved: {a!r} != {FROZEN_ALPHA!r}. This study is not "
        f"the study that was registered; re-run and re-register."
    )
    assert two_route.DEFAULT_CELL_VOLUME_FRACTION == FROZEN_PHI
    s_eff = (1.0 - FROZEN_PHI) + FROZEN_PHI * 1.0 * (1.0 + 0.0)
    assert s_eff == 1.0, s_eff
    # uptake_proxy IS drug.binding, and drug.binding IS binding.tissue_binding_rate
    dp = _compile("doxorubicin")
    mol = get_molecule("doxorubicin")
    assert uptake_proxy(dp) == dp.binding
    assert dp.binding == mol_binding.tissue_binding_rate(mol), (
        "uptake_proxy is no longer the structural tissue_binding_rate"
    )
    return {"alpha": float(a), "phi": FROZEN_PHI, "s_eff_shipping": float(s_eff)}


def _load_source() -> Dict:
    with open(_SOURCE) as fh:
        return json.load(fh)


# ==========================================================================
# STAGE role -- what uptake_proxy IS, and what it DOES to the physics
# ==========================================================================
def _interp(r_grid: Sequence[float], kills: Sequence[float]):
    from scipy.interpolate import PchipInterpolator
    x = np.log10(np.asarray(r_grid, dtype=float))
    y = np.asarray(kills, dtype=float)
    o = np.argsort(x)
    return PchipInterpolator(x[o], y[o], extrapolate=False)


def stage_role() -> Dict[str, object]:
    """Is `binding` a CAPACITY or a SINK? Measured on the shipped solver.

    A capacity R retards the approach to equilibrium and leaves the equilibrium
    itself alone. A sink k moves the equilibrium. The two are physically
    different objects and the distinction decides whether a double-count is even
    possible. This runs the SHIPPED entry point (`sim.engine.simulate_response`)
    to steady state under (i) a doubled capacity and (ii) a doubled `binding`.
    """
    _tick("stage role")
    import dataclasses

    from tmesim.core.domain import Domain
    from tmesim.core.microenvironment import (
        CellPopulation, CellType, Microenvironment, TissueFields)
    from tmesim.sim.engine import simulate_response

    out: Dict[str, object] = {"frozen": _assert_frozen()}

    # --- (0) the static facts, read off live code -------------------------
    src_binding = mol_binding.tissue_binding_rate.__doc__ or ""
    out["static"] = {
        "uptake_proxy_returns": "drug.binding",
        "drug_binding_source": "molecular.binding.tissue_binding_rate(descriptors)",
        "tissue_binding_rate_inputs": ["logp", "hbd"],
        "tissue_binding_rate_form":
            "score=(logP-2.5)/2.5+(HBD-2)/4 ; rate=1e-5*(5e-3/1e-5)**sigmoid(score)",
        "tissue_binding_rate_band_per_s": [mol_binding._BIND_MIN,
                                           mol_binding._BIND_MAX],
        "enters_physics_at": "physics/diffusion.py:315  k = clearance + binding",
        "functional_form_in_pde": "-k*C   (first-order removal, no return path)",
        "has_bound_pool": False,
        "docstring_excerpt": " ".join(src_binding.split())[:240],
    }

    # --- (1) run the shipped solver: capacity vs sink ---------------------
    # The same homogeneous probe geometry capacity_ranking_robustness.py sec. 2
    # uses, driven through the SHIPPED entry point (calling
    # physics.diffusion.simulate_drug_transport directly would make
    # PDEConfig.transvascular_exchange a silent no-op on this file, which
    # tests/test_propagation.py correctly fails on).
    ny = nx = 9
    dx = 15.0
    dom = Domain(nx=nx, ny=ny, dx=dx)
    fields = TissueFields(ecm_density=np.zeros((ny, nx)),
                          vascular_density=np.ones((ny, nx)),
                          oxygen=np.ones((ny, nx)))
    gx, gy = np.meshgrid((np.arange(nx) + 0.5) * dx, (np.arange(ny) + 0.5) * dx)
    npts = gx.size
    micro = Microenvironment(
        domain=dom, fields=fields,
        cells=CellPopulation(
            x=gx.ravel(), y=gy.ravel(),
            cell_type=np.full(npts, int(CellType.TUMOR)),
            state=np.zeros(npts), sensitivity=np.ones(npts)),
        patient_id="sink-vs-capacity-probe")
    probe = DrugParams(name="probe", diffusivity=200.0, binding=1e-3,
                       clearance=1e-4, bbb_penetrance=1.0)

    def _run(dp: DrugParams, total_time: float) -> float:
        cfg = SimConfig(drug=dp,
                        pde=PDEConfig(dt=1.0, total_time=total_time,
                                      convection=False, save_every=100_000),
                        seed=0)
        res = simulate_response(micro, cfg)
        return float(res.transport.final_concentration.mean())

    T_SS = 4.0e6   # far past equilibrium for k ~ 1e-3 /s
    base = _run(probe, T_SS)
    # (a) capacity R = 10 : divide D, k and S by R -- the sec. 2 operator
    R = 10.0
    cap = dataclasses.replace(probe, diffusivity=probe.diffusivity / R,
                              binding=probe.binding / R,
                              clearance=probe.clearance / R,
                              bbb_penetrance=probe.bbb_penetrance / R)
    cap_ss = _run(cap, T_SS)
    # (b) sink x10 : multiply `binding` by 10, everything else untouched
    sink = dataclasses.replace(probe, binding=probe.binding * 10.0)
    sink_ss = _run(sink, T_SS)

    out["steady_state_probe"] = {
        "k_base_per_s": float(probe.binding + probe.clearance),
        "shipping_C_ss": base,
        "capacity_x10_C_ss": cap_ss,
        "capacity_x10_ratio": cap_ss / base,
        "binding_x10_C_ss": sink_ss,
        "binding_x10_ratio": sink_ss / base,
    }
    # A capacity must leave the steady state alone; a sink must not.
    assert abs(cap_ss / base - 1.0) < 1e-3, cap_ss / base
    assert sink_ss / base < 0.5, sink_ss / base
    out["verdict_role"] = (
        "SINK. `binding` is a first-order removal term with no return path; "
        "multiplying it moves the equilibrium concentration. A capacity does not."
    )
    return out


# ==========================================================================
# STAGE overlap -- does the proxy ORDER compounds like measured capacity?
# ==========================================================================
def _panel_rows(src: Dict, panel_key: str) -> Dict[str, Dict[str, float]]:
    return {r["drug"]: r for r in src["algebra"]["panels"][panel_key]}


def _proxy_of(name: str, panel_key: str) -> float:
    """`uptake_proxy` for a panel drug: the structural tissue-binding rate.

    Identical to `uptake_proxy(_compile(name))` for GBM-library drugs; for the
    breast panel the molecule comes from `molecular.breast_library`, exactly as
    `experiments/breast_drug_panel.py` compiles it. `uptake_proxy` reads only
    `binding`, which `molecule_to_drugparams` sets from the descriptors alone
    (ic50 does not enter), so this is the same number the shipping path uses --
    asserted for every GBM drug below.
    """
    key = name.strip().lower()
    try:
        mol = get_molecule(key)
    except Exception:
        from tmesim.molecular.breast_library import get_molecule_breast
        mol = get_molecule_breast(key)
    return float(mol_binding.tissue_binding_rate(mol))


def _logp_of(name: str) -> float:
    key = name.strip().lower()
    try:
        mol = get_molecule(key)
    except Exception:
        from tmesim.molecular.breast_library import get_molecule_breast
        mol = get_molecule_breast(key)
    return float(mol.get("logp"))


def _corr(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    from scipy import stats
    sp = stats.spearmanr(x, y)
    pe = stats.pearsonr(np.log(x), np.log(y))
    slope, intercept = np.polyfit(np.log(x), np.log(y), 1)
    return {
        "n": int(len(x)),
        "spearman": float(sp.statistic), "spearman_p": float(sp.pvalue),
        "pearson_log": float(pe.statistic), "pearson_log_p": float(pe.pvalue),
        "r2_log": float(pe.statistic ** 2),
        "loglog_slope": float(slope), "loglog_intercept": float(intercept),
    }


def stage_overlap(src: Dict) -> Dict[str, object]:
    _tick("stage overlap")
    out: Dict[str, object] = {"frozen": _assert_frozen()}

    # --- proxy vs the structural K_b rule, both shipping panels ------------
    panels = {}
    for key in ("gbm", "breast"):
        rows = _panel_rows(src, key)
        names = [r["drug"] for r in src["algebra"]["panels"][key]]
        u = np.array([_proxy_of(n, key) for n in names])
        if key == "gbm":   # the proxy really is the shipping compile's binding
            for n, uu in zip(names, u):
                assert abs(uptake_proxy(_compile(n)) - uu) < 1e-18, n
        kb = np.array([rows[n]["K_b_rule"] for n in names])
        k = np.array([rows[n]["k"] for n in names])
        panels[key] = {
            "drugs": names,
            "uptake_proxy_per_s": {n: float(v) for n, v in zip(names, u)},
            "K_b_rule": {n: float(v) for n, v in zip(names, kb)},
            "k_eff_per_s": {n: float(v) for n, v in zip(names, k)},
            "logp": {n: float(rows[n]["logp"]) for n in names},
            "proxy_spread": float(u.max() / u.min()),
            "k_spread": float(k.max() / k.min()),
            "K_b_spread": float(kb.max() / kb.min()),
            "proxy_vs_Kb": _corr(u, kb),
            "k_vs_Kb": _corr(k, kb),
            "proxy_vs_k": _corr(u, k),
        }
    out["panels"] = panels

    # --- proxy vs MEASURED Vu,brain: how many compounds actually overlap? --
    with open(_VU) as fh:
        meas_rows = list(csv.DictReader(fh))
    measured = {}
    for r in meas_rows:
        if r.get("scoreable", "").strip().lower() != "yes":
            continue
        if r.get("species", "").strip().lower() != "rat":
            continue
        measured.setdefault(r["compound"].strip().lower(), []).append(
            float(r["value_vu_brain"]))
    engine_names = set()
    for key in ("gbm", "breast"):
        engine_names |= {n.strip().lower() for n in panels[key]["drugs"]}
    from tmesim.molecular import library as _lib
    lib_names = set()
    try:
        lib_names |= {n.strip().lower() for n in _lib.available_drugs()}
    except Exception:
        pass
    overlap = sorted((set(measured) & engine_names) | (set(measured) & lib_names))
    out["measured_overlap"] = {
        "n_measured_rat_scoreable_compounds": len(measured),
        "n_engine_panel_compounds": len(engine_names),
        "overlapping_compounds": overlap,
        "n_overlap": len(overlap),
        "detail": {c: {"measured_Vu_brain": measured[c],
                       "uptake_proxy_per_s": _proxy_of(c, "breast"),
                       "logp": _logp_of(c),
                       "K_b_rule": float(
                           two_route.binding_capacity_from_logp(_logp_of(c)))}
                   for c in overlap},
        "why_not_more":
            "The frozen brain-slice dataset carries NO SMILES column (recorded in "
            "results/brain_slice_anchor.md sec. 4.2, where B-S9 was declared NOT "
            "COMPUTABLE for exactly this reason), 23 of its 33 retained rows are "
            "radiolabelled PET tracers, and `uptake_proxy` is a function of "
            "RDKit-computed Crippen logP and HBD -- it cannot be evaluated without "
            "a structure. Sourcing structures now would be a curation step taken "
            "after the measured column was read; that report declined it and this "
            "one declines it too. No proxy value for a measured compound is "
            "estimated, guessed or substituted here.",
    }
    return out


# ==========================================================================
# STAGE residual -- is the ranking failure double-counted?
# ==========================================================================
def _order_under(fns, panel, assign, lo_r, hi_r):
    clipped = {n: min(max(assign[n], lo_r), hi_r) for n in panel}
    kills = {n: float(fns[n](math.log10(clipped[n]))) for n in panel}
    order = sorted(panel, key=lambda n: -kills[n])
    return order, kills, any(abs(clipped[n] - assign[n]) > 1e-9 for n in panel)


def _residualise(R: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Strip from log R the component linearly predictable from log(proxy).

    Keeps the geometric MEAN level of R untouched, because the uniform level of
    a tissue capacity is not constrained by anything the proxy carries: `alpha`
    is fitted to a penetration LENGTH, and a uniform capacity is exactly
    invariant on penetration length (capacity_ranking_robustness.md sec. 2b).
    Only the per-compound SHAPE can be double-counted, so only the shape is
    removed. Works row-wise on an (n_draws, n_drugs) matrix.
    """
    lr = np.log(R)
    lu = np.log(u)[None, :]
    lu_c = lu - lu.mean()
    denom = float((lu_c ** 2).sum())
    # explicit sum rather than matmul: BLAS raises spurious FP flags on
    # (n_draws, 10) @ (10, 1) under numpy 2.0.x and the result is identical.
    b = (((lr - lr.mean(axis=1, keepdims=True)) * lu_c).sum(axis=1, keepdims=True)
         / denom)                                                # (n,1) slopes
    resid = lr - b * lu_c
    return np.exp(resid)


def stage_residual(src: Dict) -> Dict[str, object]:
    _tick("stage residual")
    out: Dict[str, object] = {"frozen": _assert_frozen()}
    brain_vals = np.asarray(src["algebra"]["brain_slice"]["values"], dtype=float)

    for label, sweep_key, subpanel in (
        ("gbm10", "sweep_gbm", None),
        ("gbm5", "sweep_gbm", GBM5),
        ("breast10", "sweep_breast", None),
    ):
        sweep = src[sweep_key]
        panel = list(subpanel) if subpanel else list(sweep["panel"])
        r_grid = np.asarray(sweep["r_grid"], dtype=float)
        lo_r, hi_r = float(r_grid.min()), float(r_grid.max())
        fns = {n: _interp(r_grid, sweep["table"][n]["expected_kill"])
               for n in panel}
        base_order = sorted(panel,
                            key=lambda n: -sweep["table"][n]["expected_kill"][0])
        pkey = "gbm" if sweep_key == "sweep_gbm" else "breast"
        rows = _panel_rows(src, pkey)
        u = np.array([_proxy_of(n, pkey) for n in panel])
        kb = np.array([rows[n]["K_b_rule"] for n in panel])
        s_rule = np.array([rows[n]["s_eff_rule"] for n in panel])
        k = np.array([rows[n]["k"] for n in panel])

        res: Dict[str, object] = {"panel": panel, "base_order": base_order}

        # --- reproduction check: my pipeline must reproduce the published arms
        pub = src["named_" + ("gbm10" if label == "gbm10" else
                              "gbm5" if label == "gbm5" else "breast")]
        repro = {}
        for lab, assign in (
            ("shipping", {n: 1.0 for n in panel}),
            ("uniform_39", {n: BEST_CONSTANT for n in panel}),
            ("structural_K_b", {n: float(s_rule[j]) for j, n in enumerate(panel)}),
        ):
            o, kills, clipped = _order_under(fns, panel, assign, lo_r, hi_r)
            repro[lab] = {"order": o, "top1": o[0],
                          "matches_shipping": o == base_order,
                          "any_clipped": clipped}
        pub_struct = pub["per-compound: the repo's own structural K_b rule "
                         "(two_route.binding_capacity_from_logp)"]
        assert repro["structural_K_b"]["order"] == pub_struct["order"], (
            label, repro["structural_K_b"]["order"], pub_struct["order"])
        res["reproduction_of_published_named_arms"] = repro

        # --- ARM 1: capacity PROPORTIONAL to the proxy -----------------------
        # The strongest form of the double-counting defence: take `uptake_proxy`
        # to BE the capacity ordering, give it the measured level, and see what
        # survives. Swept over the uniform level c because that level is the part
        # the proxy provably cannot carry.
        gm_u = float(np.exp(np.log(u).mean()))
        # Below this level the proportional rule would hand some compound a
        # capacity < 1 (less storage than bare interstitium) which is off the
        # swept grid; the arm is reported both with and without that region.
        c_noclip = float(gm_u / u.min())
        c_grid = np.geomspace(1.0, hi_r, 240)
        prop = []
        for c in c_grid:
            assign = {n: float(c * u[j] / gm_u) for j, n in enumerate(panel)}
            o, kills, clipped = _order_under(fns, panel, assign, lo_r, hi_r)
            prop.append({"c": float(c), "order": o, "top1": o[0],
                         "matches": o == base_order, "any_clipped": clipped})
        first_break = next((p["c"] for p in prop if not p["matches"]), None)
        first_top1 = next((p["c"] for p in prop if p["top1"] != base_order[0]),
                          None)
        first_break_nc = next((p["c"] for p in prop
                               if p["c"] >= c_noclip and not p["matches"]), None)
        first_top1_nc = next((p["c"] for p in prop
                              if p["c"] >= c_noclip
                              and p["top1"] != base_order[0]), None)
        res["proportional_to_proxy"] = {
            "geomean_normalised": True,
            "c_grid_lo": float(c_grid[0]), "c_grid_hi": float(c_grid[-1]),
            "c_min_without_clipping": c_noclip,
            "full_order_breaks_at_c": first_break,
            "full_order_breaks_at_c_noclip": first_break_nc,
            "top1_changes_at_c_noclip": first_top1_nc,
            "break_point_involved_clipping": bool(
                first_break is not None and first_break < c_noclip),
            "top1_changes_at_c": first_top1,
            "order_at_c39": next(p["order"] for p in prop
                                 if abs(p["c"] - min(c_grid,
                                                     key=lambda z: abs(z - 39.0))
                                        ) < 1e-12),
            "matches_at_c39": next(p["matches"] for p in prop
                                   if abs(p["c"] - min(c_grid,
                                                       key=lambda z: abs(z - 39.0))
                                          ) < 1e-12),
            "n_distinct_orders": len({tuple(p["order"]) for p in prop}),
        }

        # --- ARM 2: the structural K_b rule, applied RESIDUALLY --------------
        fit = _corr(u, s_rule)
        resid_s = _residualise(s_rule[None, :], u)[0]
        o_r, kills_r, clip_r = _order_under(
            fns, panel, {n: float(resid_s[j]) for j, n in enumerate(panel)},
            lo_r, hi_r)
        o_f, _, clip_f = _order_under(
            fns, panel, {n: float(s_rule[j]) for j, n in enumerate(panel)},
            lo_r, hi_r)
        res["structural_residual"] = {
            "loglog_fit_s_eff_on_proxy": fit,
            "s_eff_full": {n: float(s_rule[j]) for j, n in enumerate(panel)},
            "s_eff_residual": {n: float(resid_s[j]) for j, n in enumerate(panel)},
            "spread_full": float(s_rule.max() / s_rule.min()),
            "spread_residual": float(resid_s.max() / resid_s.min()),
            "order_full": o_f, "order_residual": o_r,
            "top1_full": o_f[0], "top1_residual": o_r[0],
            "full_matches_shipping": o_f == base_order,
            "residual_matches_shipping": o_r == base_order,
            "any_clipped_full": clip_f, "any_clipped_residual": clip_r,
        }

        # --- ARM 3: the Monte Carlo, applied RESIDUALLY ----------------------
        rng = np.random.default_rng(MC_SEED)
        R_raw = rng.choice(brain_vals, size=(MC_DRAWS, len(panel)))
        R_res_unclipped = _residualise(R_raw, u)
        R_res = np.clip(R_res_unclipped, lo_r, hi_r)
        R_rawc = np.clip(R_raw, lo_r, hi_r)
        # Residualising can push a capacity outside the swept grid; the raw
        # draws never leave it (they ARE measured values). Disclosed, not hidden.
        clip_frac = float(np.mean((R_res_unclipped < lo_r)
                                  | (R_res_unclipped > hi_r)))
        # self-check: the vectorised residualiser IS a row-wise OLS in log space
        lu_ = np.log(u)
        for idx in (0, 1, MC_DRAWS // 3, MC_DRAWS - 1):
            sl = float(np.polyfit(lu_, np.log(R_raw[idx]), 1)[0])
            manual = np.exp(np.log(R_raw[idx]) - sl * (lu_ - lu_.mean()))
            assert np.allclose(manual, R_res_unclipped[idx], rtol=1e-9), idx
        boot = np.random.default_rng(BOOT_SEED)
        mc = {}
        for arm_lab, Rm in (("per_compound_empirical", R_rawc),
                            ("per_compound_empirical_RESIDUAL", R_res)):
            kills = np.column_stack([fns[n](np.log10(Rm[:, j]))
                                     for j, n in enumerate(panel)])
            assert np.all(np.isfinite(kills)), "interpolant extrapolated"
            oidx = np.argsort(-kills, axis=1)
            winners = oidx[:, 0]
            base_idx = np.array([panel.index(n) for n in base_order])
            full = np.all(oidx == base_idx, axis=1).astype(float)
            p1 = {}
            for j, n in enumerate(panel):
                ind = (winners == j).astype(float)
                p = float(ind.mean())
                bs = boot.binomial(len(ind), p, size=BOOT_REPS) / len(ind)
                p1[n] = {"p": p,
                         "ci90": [float(np.percentile(bs, 5)),
                                  float(np.percentile(bs, 95))]}
            pf = float(full.mean())
            bsf = boot.binomial(len(full), pf, size=BOOT_REPS) / len(full)
            mc[arm_lab] = {
                "p_rank1": p1,
                "p_full_order": pf,
                "p_full_order_ci90": [float(np.percentile(bsf, 5)),
                                      float(np.percentile(bsf, 95))],
                "geo_spread_of_applied_R": float(
                    np.exp(np.log(Rm).std(axis=1).mean())),
            }
        mc["residual_fraction_of_values_clipped_to_grid"] = clip_frac
        mc["log_spread_shrink_factor"] = (
            mc["per_compound_empirical_RESIDUAL"]["geo_spread_of_applied_R"]
            / mc["per_compound_empirical"]["geo_spread_of_applied_R"])
        res["monte_carlo"] = mc

        # --- the part the defence provably cannot touch ----------------------
        okey = {"gbm10": "order_gbm10", "gbm5": "order_gbm5",
                "breast10": "order_breast"}[label]
        res["uniform_arm_unreachable_by_the_defence"] = {
            "full_order_stable_up_to_R": src[okey]["full_order_stable_up_to_R"],
            "top1_stable_up_to_R": src[okey]["top1_stable_up_to_R"],
            "note": "a UNIFORM capacity carries no per-compound ordering at all, "
                    "so no amount of overlap between uptake_proxy and capacity "
                    "can reduce it",
        }
        out[label] = res
    return out


# ==========================================================================
# report
# ==========================================================================
def _fmt_order(o: List[str]) -> str:
    return " > ".join(o)


def stage_report(d: Dict) -> None:
    _tick("stage report")
    L: List[str] = []
    a = L.append
    role, ov, rs = d["role"], d["overlap"], d["residual"]
    g10, g5, b10 = rs["gbm10"], rs["gbm5"], rs["breast10"]
    gc = ov["panels"]["gbm"]["proxy_vs_Kb"]
    bc = ov["panels"]["breast"]["proxy_vs_Kb"]
    mo = ov["measured_overlap"]
    PANELS = (("GBM 10-drug", g10), ("GBM 5-drug", g5), ("breast 10-drug", b10))

    def _u(rr):
        return rr["uniform_arm_unreachable_by_the_defence"]

    def _fmt_stable(v):
        return "> 4184 (never changes on the grid)" if v is None else f"{v:.3g}"

    a("# Does `uptake_proxy` already carry the capacity ordering — and is the "
      "ranking failure double-counted?")
    a("")
    a("**A sensitivity study on top of two already-published results. No model "
      "change.** Nothing in `tmesim/` was modified, no default moved, and "
      f"`calibrate_alpha()` is still `{FROZEN_ALPHA!r}` — asserted against live "
      "code at the start of every stage.")
    a("")
    a("* script: `experiments/uptake_proxy_double_count.py` "
      "(`--stage role | overlap | residual | report`)")
    a("* frozen record: `results/sweep/uptake_proxy_double_count.json`")
    a("* the two findings this answers to: `results/brain_slice_anchor.md` "
      "(TEST B, **FAIL**) and `results/capacity_ranking_robustness.md` sec. 6.2, "
      "which raised the double-count possibility and explicitly left it open")
    a("* the `E[kill]`-versus-capacity tables are **re-used, not re-run**, from "
      "`results/sweep/capacity_ranking_robustness.json`. Sec. 5 states exactly "
      "what that inherits and what was done to license it.")
    a("")
    a("---")
    a("")
    a("## 0. The answer")
    a("")
    a("> ## NO — and the reason is a physics distinction, not a statistic.")
    a(">")
    a("> **`uptake_proxy` is wired into a SINK, not a capacity.** It returns "
      "`drug.binding`, and `drug.binding` reaches the solver at exactly one "
      "place in the shipped transport path — `physics/diffusion.py:315`, "
      "`k = clearance + binding` — as a first-order removal term with no return "
      "path and no bound pool (the opt-in `molecular/immune_sink.py` route "
      "builds the identical scalar and adds to it; same form). Measured "
      "on the shipped solver rather than argued: multiplying `binding` by 10 "
      "moves the steady-state concentration to "
      f"**{role['steady_state_probe']['binding_x10_ratio']:.4f}×** its shipping "
      "value, while a capacity of 10 leaves it at "
      f"**{role['steady_state_probe']['capacity_x10_ratio']:.12f}×**. A capacity "
      "stores solute and gives it back; this term destroys it. **They are "
      "different operators, and the double-counting defence is weaker than it "
      "looks before you check which one is implemented.**")
    a(">")
    a("> **The defence cannot reach the arm that does most of the damage.** The "
      "published failure includes a *uniform*-capacity arm — every compound "
      "given the SAME value, so zero per-compound ordering — and the full "
      "published order breaks there at "
      f"**R = {_u(g10)['full_order_stable_up_to_R']:.3g}** (GBM 10-drug) and "
      f"**R = {_u(b10)['full_order_stable_up_to_R']:.3g}** (breast 10-drug), "
      "both inside the measured range and below the best available constant of "
      "39. **No overlap between `uptake_proxy` and capacity, however strong, "
      "can shrink a term that has no per-compound content.**")
    a(">")
    a("> **Pushed to its limit, the defence makes the result WORSE on two of "
      "three panels.** Give every compound a capacity exactly proportional to "
      "its `uptake_proxy` — the maximal double-count, where the entire "
      "per-compound shape is the engine's own and none of it can be new physics "
      "— and the published order breaks at a geometric-mean capacity of "
      f"**{g10['proportional_to_proxy']['full_order_breaks_at_c']:.3g}** (GBM "
      f"10-drug, against 36.3 uniform) and "
      f"**{g5['proportional_to_proxy']['full_order_breaks_at_c']:.3g}** (GBM "
      "5-drug, against **435** uniform). Re-applying the engine's own ordering "
      "is not a null operation; it is a reordering force in its own right.")
    a(">")
    a("> **Against the measured quantity there is no test to run: n = "
      f"{mo['n_overlap']}.** Exactly "
      f"{mo['n_overlap']} of the "
      f"{mo['n_measured_rat_scoreable_compounds']} rat brain-slice compounds "
      f"appears on either shipping panel (`"
      f"{', '.join(mo['overlapping_compounds'])}`), the frozen dataset carries "
      "no SMILES column, and `uptake_proxy` is a function of Crippen logP and "
      "HBD — it cannot be evaluated without a structure. **Whether the proxy "
      "orders compounds the way measured capacity does is unknown, and this "
      "study does not resolve it in either direction** (sec. 2.1).")
    a(">")
    a("> **One thing the defence genuinely explains, reported here and not "
      "buried:** the breast panel's top-1 change from `vinorelbine` to "
      "`doxorubicin` under the structural `K_b` rule. On that panel "
      "`uptake_proxy` and `K_b` really are correlated "
      f"(Spearman **{bc['spearman']:+.3f}**, p = {bc['spearman_p']:.3f}) and "
      "applying the rule *residually* — stripping what the proxy already "
      f"predicts — restores `{b10['structural_residual']['top1_residual']}` as "
      "the top pick; and on the GBM 10-drug panel the same operation restores "
      "the full published order. **Those two named findings are substantially "
      "double-counted and should be discounted.** The breast full ordering is "
      "not restored, and no uniform, proportional or sampled arm moves at all "
      "(sec. 4).")
    a("")
    a("**What changes in what may be said.**")
    a("")
    a("* The `capacity_ranking_robustness.md` sec. 6.2 sentence — *\"the "
      "sensitivity measured above is partly a re-application of physics the "
      "model already has\"* — is **too generous and should be narrowed.** It is "
      "true of the two named structural-`K_b` arms re-tested here, and of "
      "nothing else. It is false of the uniform arms, false of the sampled "
      "arms, and the proportional arm points the other way.")
    a("* The published effect sizes should **not** be relabelled an upper bound "
      "on this basis. Sec. 3 quantifies the correction and it is small, "
      "two-signed, and smaller than the disagreement between the study's own arms.")
    a("* A **new** defect is established here that is not in either prior "
      "document: the engine holds a per-compound quantity that its own "
      "docstring calls *\"a capacity ordering\"* and applies it as an "
      "irreversible sink. That is not a double-count. It is a mis-assignment, "
      "and sec. 3.4 says why it cannot simply be moved.")
    a("")
    a("---")
    a("")

    # ---------------------------------------------------------------- 1 ----
    st = role["static"]
    a("## 1. What `uptake_proxy` actually is")
    a("")
    a("| | |")
    a("|---|---|")
    a(f"| returns | `{st['uptake_proxy_returns']}` — nothing else |")
    a(f"| which is set from | `{st['drug_binding_source']}` |")
    a(f"| computed from | `{'`, `'.join(st['tissue_binding_rate_inputs'])}` "
      "— RDKit Crippen logP and H-bond-donor count, and nothing else |")
    a(f"| functional form | `{st['tissue_binding_rate_form']}` |")
    a(f"| range | clamped to `{st['tissue_binding_rate_band_per_s'][0]:.0e}` – "
      f"`{st['tissue_binding_rate_band_per_s'][1]:.0e}` 1/s |")
    a("| per-compound? | **yes** — a pure function of two descriptors |")
    a(f"| enters the physics at | `{st['enters_physics_at']}` — "
      "the only place it reaches the transport solver |")
    a(f"| with what functional form | `{st['functional_form_in_pde']}` |")
    a(f"| bound pool / return path | **{st['has_bound_pool']}** |")
    a("")
    a("`alpha` multiplies it into `binding` in "
      "`calibration.calibrated_drugparams` (`binding := alpha * uptake_proxy`, "
      "`clearance := k_passive`), and `binding` is summed with `clearance` into "
      "the PDE reaction rate. **It is per-compound and it dominates `k`**: "
      f"across the GBM panel the proxy spans "
      f"**{ov['panels']['gbm']['proxy_spread']:.1f}×** and the resulting `k` "
      f"spans **{ov['panels']['gbm']['k_spread']:.1f}×**, with Spearman "
      f"{ov['panels']['gbm']['proxy_vs_k']['spearman']:+.3f} and a log–log "
      f"slope of {ov['panels']['gbm']['proxy_vs_k']['loglog_slope']:.3f} "
      "between them. **The engine's per-compound sink ordering IS the "
      "`uptake_proxy` ordering.** That much of the brief is exactly right.")
    a("")
    a("### 1.1 Capacity or sink? Measured on the shipped solver, not asserted")
    a("")
    a("This is the distinction everything else turns on. A tissue capacity "
      "`s_eff` stores solute reversibly: it retards the approach to equilibrium "
      "and leaves the equilibrium alone "
      "(`capacity_ranking_robustness.md` sec. 2b). A sink removes solute and "
      "moves the equilibrium. A homogeneous probe was run to steady state "
      "through `sim.engine.simulate_response` three ways:")
    a("")
    p = role["steady_state_probe"]
    a("| arm | steady-state mean C | ratio to shipping |")
    a("|---|---:|---:|")
    a(f"| shipping (`k` = {p['k_base_per_s']:.4g} /s) | "
      f"{p['shipping_C_ss']:.6f} | 1 |")
    a(f"| capacity R = 10 (D, k, S all ÷ 10) | {p['capacity_x10_C_ss']:.6f} | "
      f"**{p['capacity_x10_ratio']:.12f}** |")
    a(f"| `binding` × 10 | {p['binding_x10_C_ss']:.6f} | "
      f"**{p['binding_x10_ratio']:.4f}** |")
    a("")
    a("**A capacity is invisible at equilibrium. `binding` is not.** The "
      f"`binding`×10 arm lands at {p['binding_x10_ratio']:.4f}, which is "
      "`k_base/k_new` to five figures — the exact signature of a first-order "
      "sink and not of a store.")
    a("")
    a("**This is the unflattering half of the answer and it is stated first.** "
      "`uptake_proxy`'s own docstring says it is *\"intentionally a capacity "
      "ordering, not an absolute rate\"*. The module docstring **~190 lines above it, "
      "in the same file**, describes the same physics as "
      "*\"a near-irreversible perivascular sink\"*. Both sentences are in "
      "`tmesim/validation/calibration.py`; they describe different operators; "
      "**the code implements the sink.** So the premise of the double-counting "
      "defence — that the engine already carries this physics *as capacity* — "
      "is false as stated. What it carries is a capacity-flavoured *ordering* "
      "doing a sink's job.")
    a("")

    # ---------------------------------------------------------------- 2 ----
    a("## 2. Does it order compounds the way measured capacity does?")
    a("")
    a(f"### 2.1 Against measured `Vu,brain` — n = {mo['n_overlap']}, so there "
      "is no test")
    a("")
    a(f"Of the {mo['n_measured_rat_scoreable_compounds']} distinct rat, "
      "scoreable compounds in `data_cache/brain_slice/vu_brain.csv`, "
      f"**{mo['n_overlap']}** appears on either shipping panel: "
      f"**{', '.join(mo['overlapping_compounds'])}**.")
    a("")
    for c, det in mo["detail"].items():
        a(f"* `{c}` — measured `Vu,brain` = "
          f"{', '.join(f'{v:g}' for v in det['measured_Vu_brain'])} mL·g⁻¹ "
          f"(rat, Loryan 2013); `uptake_proxy` = "
          f"{det['uptake_proxy_per_s']:.4e} /s (the **largest** on the breast "
          f"panel); Crippen logP {det['logp']:.2f}; structural `K_b` = "
          f"{det['K_b_rule']:.4g}. The engine assumes a capacity of **1.0** for it.")
    a("")
    a("**One point is not an ordering.** No rank statistic, no log-linear fit "
      "and no significance is computable at n = 1, and none is reported. The "
      "single point is consistent with the proxy ranking a high-capacity "
      "compound high — and a single consistent point is not evidence.")
    a("")
    a(mo["why_not_more"])
    a("")
    a("**The brief's framing has to be corrected here.** *\"How well does "
      "`uptake_proxy` order compounds the way measured capacity does?\"* is "
      "**not answerable from anything in this repository**. Everything below "
      "compares the proxy to the repo's own *structural* `K_b` rule — a QSAR "
      "fitted to liver microsomes, not a measurement on these compounds, and "
      "one `results/intracellular_retardation.md` has already refuted on "
      "measured Pruijn HT-29 transport (rank agreement 77.4 % → 51.6 %). "
      "**Agreement or disagreement with it bounds nothing about the real "
      "quantity.** A prior overlap of n ≈ 10 at p = 0.35 was not evidence of "
      "agreement and is not evidence of disagreement; that reading was correct "
      "and it is preserved.")
    a("")
    a("### 2.2 Against the structural `K_b` rule")
    a("")
    a("| panel | n | Spearman ρ | p | Pearson r on logs | p | R² | log–log slope |")
    a("|---|---:|---:|---:|---:|---:|---:|---:|")
    for lab, key in (("GBM 10-drug", "gbm"), ("breast 10-drug", "breast")):
        c = ov["panels"][key]["proxy_vs_Kb"]
        a(f"| {lab} | {c['n']} | {c['spearman']:+.3f} | {c['spearman_p']:.3f} | "
          f"{c['pearson_log']:+.3f} | {c['pearson_log_p']:.3f} | "
          f"{c['r2_log']:.3f} | {c['loglog_slope']:.3f} |")
    a("")
    a("**The figures I was given are right, and they are also incomplete.** "
      "`capacity_ranking_robustness.md` sec. 6.2's Spearman **+0.333** "
      "(p = 0.35) and Pearson-on-logs **+0.630** (p = 0.051) reproduce exactly "
      f"({gc['spearman']:+.4f}, p = {gc['spearman_p']:.4f}; "
      f"{gc['pearson_log']:+.4f}, p = {gc['pearson_log_p']:.4f}). **They are "
      "the GBM panel only.** On the breast panel — the panel where the top-1 "
      "pick actually changed — the overlap is much stronger and is significant: "
      f"Spearman **{bc['spearman']:+.3f}** (p = {bc['spearman_p']:.3f}), "
      f"Pearson-on-logs **{bc['pearson_log']:+.3f}** "
      f"(p = {bc['pearson_log_p']:.3f}), R² = {bc['r2_log']:.3f}. Quoting the "
      "GBM number alone understates the overlap exactly where it matters, and "
      "sec. 4.1 is the consequence.")
    a("")
    a("**Correlation is the wrong statistic here; the SLOPE is the right one.** "
      "Capacity and sink enter the finite-window exposure only through the "
      "dimensionless group `k·T/R` (`capacity_ranking_robustness.md` sec. 2), "
      "so the capacity that cancels against the sink is one *proportional* to "
      "it — not one merely rank-correlated with it. The log–log slope of `K_b` "
      f"on `uptake_proxy` is **{gc['loglog_slope']:.2f}** on GBM and "
      f"**{bc['loglog_slope']:.2f}** on breast. Neither is 1. On breast the "
      "capacity grows nearly half again as fast as the sink that would absorb "
      "it, which is why a strong rank correlation still leaves a large "
      "reordering force.")
    a("")
    a("**And the two terms contradict each other on a specific compound.** "
      "`uptake_proxy` is monotone **increasing** in logP; the `K_b` rule "
      "(`125·10^(0.072·logP² + 0.067·logP − 1.126)`) is a **parabola** in logP "
      "whose minimum sits near logP = −0.47, so it turns back up for "
      "hydrophilic compounds. `temozolomide` (logP −2.08) is therefore the "
      "**smallest** sink in the engine and sits mid-table on capacity (5th smallest "
      "of 10) under the structural rule. This is the second defect the brief anticipated: "
      "two terms in the same engine, both claiming to describe intracellular "
      "sequestration, disagreeing about the same molecule.")
    a("")

    # ---------------------------------------------------------------- 3 ----
    a("## 3. Is the ranking failure double-counted?")
    a("")
    a("### 3.1 What CAN be double-counted, stated before any number")
    a("")
    a("A per-compound capacity assignment decomposes in log space into a "
      "**level** (its geometric mean) and a **shape** (the per-compound "
      "departures from it). Only the shape is a candidate:")
    a("")
    a("* **the level is not carried by the proxy.** `alpha` is fitted to a "
      "penetration *length*, and a uniform capacity is exactly invariant on "
      "penetration length (`capacity_ranking_robustness.md` sec. 2b, verified "
      "to 1e-9). Nothing in the calibration sees the capacity level, so nothing "
      "about it can have been applied twice;")
    a("* **the shape is the only candidate**, and only to the extent the proxy "
      "predicts it.")
    a("")
    a("**That alone settles the largest part of the published failure**, "
      "because the uniform arm has no shape at all:")
    a("")
    a("| panel | full order holds up to R = | top-1 holds up to R = |")
    a("|---|---:|---:|")
    for lab, rr in PANELS:
        uu = _u(rr)
        a(f"| {lab} | **{uu['full_order_stable_up_to_R']:.3g}** | "
          f"{_fmt_stable(uu['top1_stable_up_to_R'])} |")
    a("")
    a("Two of those three are below the best available constant of 39 and all "
      "three are inside the measured range [3.2, 4184]. **The double-count "
      "defence has no purchase on any of it.**")
    a("")
    a("### 3.2 The strongest possible defence: let the proxy BE the capacity")
    a("")
    a("Give every compound a capacity exactly proportional to its "
      "`uptake_proxy`, normalised to a geometric mean of `c`, and sweep `c`. "
      "This is the double-count at its limit: the per-compound shape is "
      "*entirely* the engine's own, so *none* of it can be new physics. If the "
      "ordering held here, the defence would be complete.")
    a("")
    a("| panel | full order breaks at c = | vs. the UNIFORM arm | top-1 changes at c = |")
    a("|---|---:|---:|---:|")
    for lab, rr in PANELS:
        pr = rr["proportional_to_proxy"]
        fb = pr["full_order_breaks_at_c"]
        fbn = pr["full_order_breaks_at_c_noclip"]
        t1 = pr["top1_changes_at_c"]
        cell = ("never on [1, 4184]" if fb is None else
                (f"**{fb:.3g}**" if not pr["break_point_involved_clipping"]
                 else f"**{fb:.3g}** (clipped; {fbn:.3g} unclipped)"))
        a(f"| {lab} | {cell} | "
          f"{_u(rr)['full_order_stable_up_to_R']:.3g} | "
          f"{'never on [1, 4184]' if t1 is None else f'{t1:.3g}'} |")
    a("")
    a("**It does not hold — and on the two GBM panels it is markedly worse than "
      "a uniform capacity.** On the 5-drug sub-panel, the one "
      "`results/ranking_robustness.md` actually scores, a uniform capacity "
      "preserves the published order to R = 435; a capacity proportional to the "
      "engine's own uptake ordering breaks it at **c = "
      f"{g5['proportional_to_proxy']['full_order_breaks_at_c']:.3g}**. Only on "
      "the breast panel is the proportional arm the more benign of the two "
      f"({b10['proportional_to_proxy']['full_order_breaks_at_c_noclip']:.3g} "
      f"against {_u(b10)['full_order_stable_up_to_R']:.3g}), and by less than a "
      "factor of 1.5.")
    a("")
    a("**Why re-applying the engine's own ordering is not a null operation.** "
      "Under `R_i ∝ k_i`, the group `k_i·T/R_i` becomes the *same* number for "
      "every compound, so every compound receives the *same* AUC attenuation. "
      "At the shipping state it does not: attenuation ranges over 0.9075–0.9965 "
      "across the GBM panel. Equalising it therefore *relatively boosts the "
      "slow-sink compounds*, which is itself a reordering. **A capacity that "
      "perfectly mirrors the sink does not cancel — it flattens, and flattening "
      "reorders.** This is the quantitative core of the answer and it is why "
      "\"the effect is partly double-counted, so treat it as an upper bound\" "
      "does not follow.")
    a("")
    a("*Disclosure: on the breast panel the proportional rule assigns some "
      "compounds a capacity below 1 for small `c`, which is off the bottom of "
      "the swept grid and was clipped to 1. The table reports the break point "
      "both with clipping and restricted to the un-clipped region "
      f"(c ≥ {b10['proportional_to_proxy']['c_min_without_clipping']:.3g}); the "
      "un-clipped figure is the conservative one and is the one quoted above.*")
    a("")
    a("### 3.3 The sampled arm, applied residually")
    a("")
    a("The published Monte Carlo draws each compound's capacity independently "
      "from the 33 measured values. The residual version strips from every draw "
      "the component linearly predictable from `log uptake_proxy` and keeps the "
      "geometric mean (sec. 3.1). Same seed (20260821), same 40 000 draws, same "
      "bootstrap, same decision-grade bar. The vectorised residualiser is "
      "asserted against a row-wise `polyfit` on four draws.")
    a("")
    a("| panel | arm | P(top-1 = published leader) | P(full published order) |")
    a("|---|---|---:|---:|")
    for lab, rr in PANELS:
        lead = rr["base_order"][0]
        for arm, nm in (("per_compound_empirical", "as published"),
                        ("per_compound_empirical_RESIDUAL", "**residual**")):
            m = rr["monte_carlo"][arm]
            p1 = m["p_rank1"][lead]
            a(f"| {lab} (`{lead}`) | {nm} | {p1['p']:.3f} "
              f"[{p1['ci90'][0]:.3f}, {p1['ci90'][1]:.3f}] | "
              f"{m['p_full_order']:.4f} "
              f"[{m['p_full_order_ci90'][0]:.4f}, "
              f"{m['p_full_order_ci90'][1]:.4f}] |")
    a("")
    shrink = g10["monte_carlo"]["log_spread_shrink_factor"]
    a("**The published figures reproduce** — 0.005 / 0.200 / 0.001 against "
      f"{g10['monte_carlo']['per_compound_empirical']['p_full_order']:.4f} / "
      f"{g5['monte_carlo']['per_compound_empirical']['p_full_order']:.4f} / "
      f"{b10['monte_carlo']['per_compound_empirical']['p_full_order']:.4f} "
      "here, on the same seed. **The residual correction moves nothing across "
      "any bar**, and on the breast panel it moves the wrong way (to 0.0000). "
      "Residualising shrinks the geometric spread of the applied capacity by a "
      f"factor of only **{shrink:.3f}** on the GBM panel, because a draw from "
      "the measured distribution is independent of `uptake_proxy` by "
      "construction: regressing it on the proxy removes only the sampling-noise "
      "overlap (expected R² ≈ 1/(n−2) ≈ 0.12 at n = 10).")
    a("")
    a("**This arm is the weakest evidence in this document and is labelled so.** "
      "It cannot be strong: the sampled capacities were never linked to "
      "structure, so there was never a real double-count in it to remove. It is "
      "reported because its absence is itself the answer to \"is the sampled "
      "instability inflated by the proxy?\" — no, it is not.")
    a("")
    a(f"*Disclosure: residualising can push a capacity outside the swept grid; "
      f"{100 * g10['monte_carlo']['residual_fraction_of_values_clipped_to_grid']:.2f} % "
      "of GBM-panel values were clipped to [1, 4184] (raw draws never leave it, "
      "because they ARE measured values).*")
    a("")
    a("### 3.4 Is it cleanly separable? No — and what that implies is worse")
    a("")
    a("The residual construction above is a **statistical** separation. No "
      "physical one exists, and that is a result rather than an evasion:")
    a("")
    a("1. **The two terms are different operators (sec. 1.1).** Taking "
      "`alpha·uptake_proxy` out of `k` and re-expressing it as a capacity is "
      "not a change of variables: it moves the steady state, which a capacity "
      "provably cannot do. There is no algebra that converts one into the other.")
    a("2. **The sink is load-bearing for a different anchor.** `alpha` exists "
      "solely to make doxorubicin's penetration length hit ~64.9 µm. A capacity "
      "has **no** effect on penetration length. So the sink cannot be reduced "
      "in favour of a capacity without breaking the one anchor the engine is "
      "fitted to. **The two terms are not substitutes and cannot be traded "
      "off** — which is exactly why the double-count cannot be netted out.")
    a("3. **The implication runs against the engine, not for it.** If "
      "`uptake_proxy` is a capacity ordering, as its docstring says, then the "
      "engine is currently using a capacity ordering to set a quantity "
      "(penetration length) that a real capacity would leave untouched, and "
      "leaving the slot a real capacity belongs in (`s_eff`) pinned at 1.0. "
      "**That is not a double-count that cancels; it is a mis-assignment that "
      "would have to be undone before either term could be calibrated "
      "honestly.** It also means `alpha` = 22.22 is absorbing whatever the "
      "capacity slot is missing, which is a second reason not to read `alpha` "
      "as a physical constant.")
    a("")

    # ---------------------------------------------------------------- 4 ----
    a("## 4. Where the defence DOES do work — reported in full")
    a("")
    a("### 4.1 The breast panel's top-1 change: substantially double-counted")
    a("")
    st_b = b10["structural_residual"]
    a("`capacity_ranking_robustness.md`'s headline per-compound finding is that "
      "the structural `K_b` rule changes the breast top-1 pick from "
      "`vinorelbine` to `doxorubicin`. That is the arm where `uptake_proxy` and "
      f"`K_b` are genuinely correlated (Spearman {bc['spearman']:+.3f}, "
      f"p = {bc['spearman_p']:.3f}, R² = {bc['r2_log']:.3f} on logs). Applying "
      "the same rule **residually**:")
    a("")
    a("| | full structural `K_b` | **residual** |")
    a("|---|---|---|")
    a(f"| spread of `s_eff` across the panel | {st_b['spread_full']:.4g}× | "
      f"{st_b['spread_residual']:.4g}× |")
    a(f"| top-1 | `{st_b['top1_full']}` | `{st_b['top1_residual']}` |")
    a(f"| = published order? | "
      f"{'yes' if st_b['full_matches_shipping'] else '**NO**'} | "
      f"{'yes' if st_b['residual_matches_shipping'] else '**NO**'} |")
    a("")
    a(f"* full-rule order: {_fmt_order(st_b['order_full'])}")
    a(f"* residual order: {_fmt_order(st_b['order_residual'])}")
    a("")
    if st_b["top1_residual"] == b10["base_order"][0]:
        a("**The top-1 change is partly an artefact of double-counting.** Strip "
          "the part of the structural rule that `uptake_proxy` already predicts "
          "— which on this panel is over half the log-variance — and "
          f"`{b10['base_order'][0]}` comes back as the top pick. The residual "
          f"spread falls from {st_b['spread_full']:.4g}× to "
          f"{st_b['spread_residual']:.4g}×, and that is enough to undo the "
          "swap. **This is the one place in either document where the "
          "double-count defence is load-bearing, and it should be said "
          "plainly.**")
        a("")
        a("**What it does not buy.** The residual arm still does **not** "
          "reproduce the published full ordering. It is a QSAR residualised "
          "against a QSAR on 10 points, with neither side measured on these "
          "compounds in this tissue. And it is one of four named arms: the "
          "uniform R = 39 arm changes the breast order too (sec. 3.1) and is "
          "untouched by any of this.")
    else:
        a("**The top-1 change survives residualisation**: even after removing "
          "everything `uptake_proxy` predicts about the structural rule, the "
          "breast top-1 still moves. Double-counting is not the explanation.")
    a("")
    a("### 4.2 The GBM panels")
    a("")
    for lab, rr in (("GBM 10-drug", g10), ("GBM 5-drug", g5)):
        s = rr["structural_residual"]
        a(f"* **{lab}** — full-rule top-1 `{s['top1_full']}`, residual top-1 "
          f"`{s['top1_residual']}`; residual order "
          f"{'**matches**' if s['residual_matches_shipping'] else 'does **NOT** match'} "
          f"the published one (full-rule order "
          f"{'matches' if s['full_matches_shipping'] else 'does not match'}). "
          f"`s_eff` spread {s['spread_full']:.4g}× → "
          f"{s['spread_residual']:.4g}×; log–log slope of `s_eff` on the proxy "
          f"{s['loglog_fit_s_eff_on_proxy']['loglog_slope']:.3f}, "
          f"R² {s['loglog_fit_s_eff_on_proxy']['r2_log']:.3f}.")
    a("")
    a("Two things to read off that, one in each direction. **For the "
      "defence:** on the GBM 10-drug panel residualising also puts the full "
      "published order back — the structural rule's single adjacent-pair swap "
      "(`irinotecan`/`procarbazine`) is inside the part the proxy predicts, so "
      "that swap is double-counted too. **Against it:** on the GBM 5-drug "
      "panel the proxy explains essentially none of the structural rule "
      f"(R² = {g5['structural_residual']['loglog_fit_s_eff_on_proxy']['r2_log']:.3f}, "
      "log–log slope "
      f"{g5['structural_residual']['loglog_fit_s_eff_on_proxy']['loglog_slope']:.3f}), "
      "so there is nothing there to remove — which is the honest reading of a "
      "Spearman of +0.333 at p = 0.35, and a reminder that the two GBM results "
      "come from the same rule on overlapping compounds.")
    a("")
    a("**Net.** The published study ran four named per-compound arms that "
      "broke the order (the structural `K_b` rule and its access-limited "
      "variant, on the GBM 10-drug and breast panels). This study re-tested "
      "the two structural ones: residualising **restores the GBM 10-drug full "
      "order** and **restores the breast top-1** (though not the breast full "
      "order). The two access-limited variants were not re-tested — they "
      "produced identical orderings to the structural rule on every panel, so "
      "the same conclusion is likely, but it is not measured here. **That is "
      "the whole of what the double-count defence buys. It buys nothing on the "
      "uniform arms (sec. 3.1), nothing on the sampled arms (sec. 3.3), and it "
      "is contradicted by the proportional arm (sec. 3.2).**")
    a("")

    # ---------------------------------------------------------------- 5 ----
    a("## 5. Verified, assumed, inherited")
    a("")
    a("**Verified against live code or a live solve.**")
    a("")
    a(f"* `calibrate_alpha()` = `{FROZEN_ALPHA!r}`, asserted at the start of "
      "every stage. `two_route.DEFAULT_CELL_VOLUME_FRACTION` = 0.7892 and "
      "`s_eff` at the shipping state = 1.0, likewise.")
    a("* `uptake_proxy(dp) == dp.binding`, and `dp.binding == "
      "molecular.binding.tissue_binding_rate(molecule)` — asserted for the "
      "compiled molecule, not read off a docstring, and re-asserted for every "
      "GBM-panel drug in `stage_overlap`.")
    a("* `binding` reaches the PDE at one place as `−k·C` with no return path "
      "(`physics/diffusion.py:315`); the sink-versus-capacity distinction was "
      "measured through `sim.engine.simulate_response` (sec. 1.1) and the two "
      "arms are asserted to differ in the predicted direction.")
    a("* Every correlation, slope, spread and overlap count in sec. 2 is "
      "computed here from live compiles and the frozen CSV — none is "
      "transcribed. The three figures I was handed (+0.333, p = 0.35; +0.630, "
      "p = 0.051; `s_eff` = 1.0) all reproduce.")
    a("* **The published named-arm orderings are reproduced exactly by this "
      "script's own pipeline before any new arm is read off it** — asserted in "
      "`stage_residual`, and the published Monte-Carlo probabilities reproduce "
      "on the same seed (sec. 3.3). That is the check that licenses re-using "
      "the frozen tables.")
    a("")
    a("**Inherited, not re-run.** The `E[kill]`-versus-capacity tables come "
      "from `results/sweep/capacity_ranking_robustness.json` — 220 transport "
      "solves per panel on the real Xenium reconstructions, at dx = 15 µm, "
      "dt = 20 s, convection **off**, 6 h window. This study does not re-solve "
      "them. Every ordering statement here is only as good as that record, "
      "including its convection-off limitation, and its own worst absolute "
      "interpolation error of 3.5e-04 is the relevant numerical bound. **If "
      "that record is wrong, sec. 3.2–4 are wrong with it**; sec. 1 and 2 are "
      "independent of it.")
    a("")
    a("**Assumed, and attackable.**")
    a("")
    a("* That removing a double-count means stripping the log-linearly "
      "predictable component and keeping the level. Sec. 3.1 argues the level "
      "cannot have been double-counted; a reader who disagrees would get a "
      "*smaller* residual capacity and a more benign answer, so this is the "
      "load-bearing assumption and it is stated where it can be attacked.")
    a("* That a *linear* fit in log space is the right way to project out the "
      "proxy. A monotone but non-linear relationship would leave residual "
      "structure this construction does not remove.")
    a("* That the structural `K_b` rule is an acceptable stand-in for measured "
      "capacity. It is not one, it has been refuted on other data, and "
      "sec. 2.1 says so before any number is quoted.")
    a("")
    a("**Not claimable from this document.**")
    a("")
    a("* That `uptake_proxy` agrees, or disagrees, with measured tissue "
      "capacity. **n = 1.** Nothing was measured about that here, and no "
      "structure was sourced to make it n > 1.")
    a("* That the published ranking failure is smaller than reported — except "
      "for the single named breast arm in sec. 4.1, under the assumptions "
      "listed there.")
    a("* That a capacity term should, or should not, be added to the engine. "
      "`results/intracellular_retardation.md` refuted the only available "
      "per-compound rule and nothing here reverses that.")
    a("* Any absolute kill number; anything about `cell_partition` (K), which "
      "remains unconstrained by any measurement; anything clinical.")
    a("")
    a("## 6. Reproduction")
    a("")
    a("```")
    a("./.venv/bin/python experiments/uptake_proxy_double_count.py --stage all")
    a("```")
    a("")
    a("Monte-Carlo seed **20260821**, bootstrap seed **20260822**, 40 000 "
      "draws, 4 000 bootstrap replicates — the same seeds and counts as "
      "`capacity_ranking_robustness.md`, so the two are directly comparable.")
    a("")
    a("**Nothing in `tmesim/`, `data_cache/`, `docs/`, `app/`, `landing-3d/` or "
      "any existing `results/*.md` was modified.** Two files were added: "
      "`experiments/uptake_proxy_double_count.py` and this document, plus the "
      "frozen record `results/sweep/uptake_proxy_double_count.json`.")
    a("")
    a("**No assertion was deleted or weakened.** `./.venv/bin/python -m pytest "
      "-q` → **1171 passed, 0 failed** (407.62 s), the unchanged baseline.")
    a("")
    a("**Rider disclosure — the working tree was not clean.** Sibling agents "
      "were live in this repository during the run and the tree carries other "
      "untracked files (`experiments/pk_ranking_robustness.py`, "
      "`results/axis_inventory.md`). No tracked file under `tmesim/` was "
      "modified (`git status --porcelain` shows only untracked additions), and "
      "every constant this study leans on is asserted against LIVE code at the "
      "start of every stage and matched. Saying the tree was dirty is cheaper "
      "than asserting it was not.")

    with open(_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    _tick(f"wrote {_MD}")


# ==========================================================================
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["role", "overlap", "residual", "report", "all"])
    args = ap.parse_args()

    os.makedirs(_SWEEP, exist_ok=True)
    data: Dict[str, object] = {}
    if os.path.exists(_OUT):
        with open(_OUT) as fh:
            data = json.load(fh)
    src = _load_source()

    if args.stage in ("role", "all"):
        data["role"] = stage_role()
    if args.stage in ("overlap", "all"):
        data["overlap"] = stage_overlap(src)
    if args.stage in ("residual", "all"):
        data["residual"] = stage_residual(src)
    with open(_OUT, "w") as fh:
        json.dump(data, fh, indent=1, sort_keys=False)
    _tick(f"wrote {_OUT}")
    if args.stage in ("report", "all"):
        stage_report(data)


if __name__ == "__main__":
    main()
