"""HOW OFTEN does the headline delivery/dose verdict flip?

Driver for ``results/verdict_inversion_census.md``.

THE QUESTION
============
``results/permeability_split.md`` §4.2 established that the product's headline
per-cell verdict — "was the drug never delivered, or delivered and too weak?" —
**inverts** for daunorubicin when the transvascular exchange term (which ships
OFF) is switched ON, at the shipped parameters with nothing scaled. It also
found doxorubicin sitting at 0.497 : 0.503, one cell from the same flip.

That was 4 compounds on 1 section, and three of the four are not in the shipped
library at all. The panel now tells a user the verdict CAN reverse. It cannot
tell them how often. **A caveat without a frequency is honest and not
actionable.** This file supplies the frequency.

WHAT IT MEASURES
================
For **every compound the product will actually run** — ``sorted(set(
molecular.library.GBM_DRUGS) | {"palbociclib"})``, which is exactly what
``app/server/api.py::_known_drugs`` offers before any user-compiled compound is
added — on the flagship section, exchange OFF versus exchange ON, at the
SHIPPED parameters with **nothing scaled**:

* the delivery : dose split and the dominant limiter under each arm;
* whether the dominant limiter INVERTS;
* the per-cell label churn (how many of the 40,887 cells change label, and in
  which direction);
* how far each non-inverting compound sits from 0.5, in SHARE and in CELLS —
  because a compound reporting 0.497 : 0.503 states a dominant limiter with
  full confidence while one cell from the opposite answer, which is arguably
  worse than an inversion because it does not announce itself.

STAGES
======
``reproduce`` — gate. Re-runs doxorubicin and daunorubicin through the SAME
compile path ``results/permeability_split.md`` §4.2 used
(``isomer_penetration_map.compile_all``, shared IC50) and requires the
delivery : dose splits to come back BIT-IDENTICAL to the published ones. If the
harness cannot reproduce the one result this file is extending, the census does
not run.

``census`` — the 13 shipped compounds, both arms.

RULES OBSERVED
==============
* **Nothing shipped is changed and no parameter is tuned.**
  ``PDEConfig.transvascular_exchange``, ``transvascular.SPLIT_DEFAULT`` and
  ``drugs.PS_SCALE_S_INV`` are asserted at their shipped values (``False``,
  ``False``, ``1.0``) at the TOP and the BOTTOM of every stage, and
  ``calibrate_alpha()`` is asserted live at ``22.223477711172563`` in the same
  two places. Nothing is scaled: the compiled ``DrugParams`` go to the solver
  exactly as the product compiles them.
* Configuration travels with every number (see ``_tissue_cfg``).
* Every solve goes through ``tmesim.sim.engine.simulate_response``, the ONE
  production place that turns ``PDEConfig.transvascular_exchange`` into an
  ``exchange_k`` argument. This file calls ``simulate_drug_transport`` nowhere
  and is therefore not a second wiring point.

BEVACIZUMAB IS THE DEGENERATE CASE AND IS SCORED AS ONE.
``bbb_penetrance == 0.0`` makes the delivery source identically zero, so no
drug enters under EITHER law, every surviving tumour cell is UNREACHED, and the
DELIVERY verdict is an **input, not a finding**. It is reported in full and
excluded from the inversion-rate denominator, with the exclusion stated. Any
compound whose compiled ``bbb_penetrance`` is exactly zero is detected the same
way, by value rather than by name.

THE FOUR CHECKS, AND THE DEFECT EACH ONE EXISTS TO CATCH
========================================================
Run with ``--check`` (they also run inline at the end of ``census``).

A ``check_split_is_not_inert`` — the split adds a back-flux ``-PS*C`` to a
  source that was ``+PS*Cp``, so with ``C >= 0`` the ON arm's exposure must be
  STRICTLY below the OFF arm's wherever any drug arrived. An exchange flag
  dropped at the engine boundary (the real defect
  ``tmesim/validation/propagation.py`` records) makes the two arms bit-identical
  and would report a **0 % inversion rate with every stage green**. This is the
  vacuous-green trap for this file and it is checked first.
B ``check_reproduces_published`` — the §4.2 gate above.
C ``check_churn_is_measured`` — the per-cell churn must be counted BETWEEN the
  arms and must reconcile with the transition matrix. Catches a churn diff
  taken against the wrong arm (which reports 0 changed cells, also green).
D ``check_degenerate_is_degenerate`` — a ``bbb_penetrance == 0`` compound must
  come back bit-identical across the two arms with delivery share exactly 1.0.
  Catches a degenerate compound quietly inflating the stable-verdict count.

Each of the four was planted as a defect and watched to bite; the verbatim
failures are in ``results/verdict_inversion_census.md`` §6.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from collections import Counter
from typing import Dict, List, Optional, Sequence

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from tmesim.core import PDEConfig, PKConfig, SimConfig                 # noqa: E402
from tmesim.molecular.binding_site_barrier import (                    # noqa: E402
    apply_barrier_to_drugparams,
)
from tmesim.molecular.breast_library import get_molecule_breast        # noqa: E402
from tmesim.molecular.library import GBM_DRUGS, get_molecule           # noqa: E402
from tmesim.sim.engine import simulate_response                        # noqa: E402
from tmesim.sim.mechanism import LABELS                                # noqa: E402
from tmesim.treatment.drugs import PS_SCALE_S_INV                      # noqa: E402
from tmesim.treatment.transvascular import SPLIT_DEFAULT               # noqa: E402
from tmesim.validation import vessel_profile as VP                     # noqa: E402
from tmesim.validation.calibration import (                            # noqa: E402
    calibrate_alpha,
    calibrated_drugparams,
)

ALPHA_LIVE = 22.223477711172563

#: Drugs curated in ``molecular.breast_library`` that are nonetheless valid
#: GBM-context names. Mirrors ``app/server/api.py::_CROSS_LIBRARY_GBM_DRUGS``,
#: which is the single source of truth for what the product will run.
CROSS_LIBRARY = ("palbociclib",)

#: Flagship section arm, taken UNCHANGED from
#: ``results/permeability_split.md`` §4.2 / ``results/kyle2004_penetration.md``
#: §7 so the numbers below are directly comparable with the published ones.
TISSUE_DX_UM = 12.0
TISSUE_REPAIR = 1
TISSUE_DT_S = 10.0
TISSUE_TOTAL_S = 3600.0 * 6
TISSUE_CONVECTION = True
TISSUE_SAVE_EVERY = 500
TISSUE_SEED = 0

#: ``results/permeability_split.md`` §4.2, verbatim. The gate, not a target.
PUBLISHED_SPLITS: Dict[str, Dict[str, Dict[str, float]]] = {
    "doxorubicin": {
        "OFF": {"delivery": 0.38444610712556665, "dose": 0.6155538928744333},
        "ON": {"delivery": 0.49685314685314685, "dose": 0.5031468531468531},
    },
    "daunorubicin": {
        "OFF": {"delivery": 0.42881624044292116, "dose": 0.5711837595570788},
        "ON": {"delivery": 0.5411667175705747, "dose": 0.4588332824294254},
    },
}

_DATA_PATH = os.path.join(_ROOT, "results", "verdict_inversion_census_data.json")


# --------------------------------------------------------------------------- #
# Frozen-shipped-state guard
# --------------------------------------------------------------------------- #
def _assert_frozen(where: str) -> float:
    """Re-assert that nothing shipped moved. Called at the TOP and the BOTTOM of
    every stage, as ``results/permeability_split.md`` does."""
    alpha = calibrate_alpha()
    if alpha != ALPHA_LIVE:
        raise AssertionError(
            f"[{where}] calibrate_alpha() is {alpha!r}, not the asserted live "
            f"value {ALPHA_LIVE!r}. Every number in this file would be on a "
            "different calibration and not comparable with results/.")
    if PDEConfig().transvascular_exchange is not False:
        raise AssertionError(
            f"[{where}] PDEConfig.transvascular_exchange is not False; the "
            "shipped default has been changed and the OFF arm is no longer the "
            "shipped arm.")
    if SPLIT_DEFAULT is not False:
        raise AssertionError(f"[{where}] transvascular.SPLIT_DEFAULT is not False")
    if PS_SCALE_S_INV != 1.0:
        raise AssertionError(
            f"[{where}] drugs.PS_SCALE_S_INV is {PS_SCALE_S_INV!r}, not the "
            "shipped 1.0; every absolute exposure below would be on a "
            "different scale and every threshold call with it.")
    return alpha


def shipped_library() -> List[str]:
    """Exactly what ``app/server/api.py::_known_drugs`` offers before any
    user-compiled compound is added."""
    return sorted(set(GBM_DRUGS) | set(CROSS_LIBRARY))


def _molecule_any(name: str):
    try:
        return get_molecule(name)
    except KeyError:
        return get_molecule_breast(name)


def compile_shipped(name: str, alpha: float):
    """Compile ONE compound the way the product compiles it.

    ``calibrated_drugparams`` then ``apply_barrier_to_drugparams`` is the path
    ``experiments/real_drug_panel.py::compile_real_drug`` and the Studio's own
    ``drug_rows`` both use. Nothing is scaled and nothing is overridden.
    """
    dp = calibrated_drugparams(name, alpha)
    return apply_barrier_to_drugparams(dp, _molecule_any(name))


def _tissue_cfg() -> str:
    return (f"10x Xenium human GBM flagship section, dx = {TISSUE_DX_UM:g} um, "
            f"repair = {TISSUE_REPAIR}, dt = {TISSUE_DT_S:g} s, "
            f"{TISSUE_TOTAL_S / 3600:.0f} h window, convection "
            f"{'ON' if TISSUE_CONVECTION else 'OFF'}, "
            f"PKConfig(plasma_peak={PKConfig().plasma_peak:g} uM, "
            f"half_life={PKConfig().half_life:g} s, "
            f"n_doses={PKConfig().n_doses})")


_MICRO_CACHE: Dict = {}


def _micro():
    if "m" not in _MICRO_CACHE:
        t0 = time.time()
        micro, _, _ = VP.build_arm_microenvironment(dx_um=TISSUE_DX_UM,
                                                    repair=TISSUE_REPAIR)
        _MICRO_CACHE["m"] = micro
        _MICRO_CACHE["build_s"] = time.time() - t0
    return _MICRO_CACHE["m"]


#: Indirection so the inert-flag defect can be planted at exactly the place the
#: real one lives (``tmesim/validation/propagation.py``): the flag's journey
#: from PDEConfig to the solver.
def _exchange_arm(exchange: bool) -> bool:
    return bool(exchange)


def _is_degenerate(dp) -> bool:
    """A compound whose delivery source is identically zero, detected BY VALUE
    rather than by name. Its own function so the 'degenerate compound quietly
    counted as a stable verdict' defect has somewhere to be planted."""
    return float(dp.bbb_penetrance) == 0.0


def _run(dp, exchange: bool) -> Dict:
    """One solve, one delivery law, through the single production wiring point."""
    pde = PDEConfig(transvascular_exchange=_exchange_arm(exchange),
                    dt=TISSUE_DT_S, total_time=TISSUE_TOTAL_S,
                    convection=TISSUE_CONVECTION, save_every=TISSUE_SAVE_EVERY)
    t0 = time.time()
    res = simulate_response(_micro(), SimConfig(drug=dp, pde=pde, seed=TISSUE_SEED),
                            attribute=True)
    counts = {k: int(v) for k, v in res.mechanism["counts"].items()}
    bd = {k: float(v) for k, v in res.mechanism["survival_breakdown"].items()}
    n_tumor = int(res.mechanism["n_tumor"])
    n_surv = int(res.mechanism["n_survived"])
    return {
        "kill_fraction": float(res.kill_fraction),
        "n_tumor": n_tumor,
        "n_survived": n_surv,
        "counts": counts,
        "survival_breakdown": bd,
        "dominant_limiter": max(bd, key=bd.get),
        "dominant_share": float(max(bd.values())),
        "unreached_fraction_of_tumor": counts["SURVIVED_UNREACHED"] / n_tumor
                                       if n_tumor else None,
        "mean_auc": float(np.asarray(res.fields["auc"]).mean()),
        "max_auc": float(np.asarray(res.fields["auc"]).max()),
        "solve_s": time.time() - t0,
        "labels": np.asarray(res.mechanism["labels"]),
    }


def _margin(arm: Dict) -> Dict:
    """How far this arm's verdict sits from an even split, in SHARE and in CELLS.

    TWO different questions are scored, because the panel reports two different
    things and they are not the same claim:

    * ``delivery_share_of_survivors`` — ``SURVIVED_UNREACHED / n_survived``,
      exactly what ``mechanism.attribute``'s ``survival_breakdown`` returns and
      exactly the number ``results/permeability_split.md`` §4.2 quotes. Its
      competitors include RESISTANCE and CHANCE.
    * ``delivery_share_of_exposure_limited`` — ``n_delivery / (n_delivery +
      n_dose)``. **This is the product's headline discrimination**: "the drug
      never arrived" versus "it arrived and there was not enough of it", asked
      only of the cells for which that question is even posed. A compound can
      have a rock-solid four-way dominant limiter and a coin-flip answer to
      this one.

    ``cells_to_flip`` is the smallest number of surviving cells that would have
    to move from the leading exposure limiter to the other for the reported
    delivery-vs-dose answer to reverse. Moving one cell moves the difference by
    two, so it is ``floor(|n_delivery - n_dose| / 2) + 1``.
    """
    c = arm["counts"]
    n_del = int(c["SURVIVED_UNREACHED"])
    n_dose = int(c["SURVIVED_UNDERDOSED"])
    n_exp = n_del + n_dose
    n_surv = int(arm["n_survived"])
    bd = arm["survival_breakdown"]
    share_surv = float(bd["delivery"]) if n_surv else float("nan")
    share_exp = (n_del / n_exp) if n_exp else float("nan")
    return {
        "n_delivery": n_del,
        "n_dose": n_dose,
        "n_exposure_limited": n_exp,
        "n_survived": n_surv,
        "delivery_share_of_survivors": share_surv,
        "delivery_share_of_exposure_limited": share_exp,
        "distance_from_even_share": abs(share_exp - 0.5) if n_exp
                                    else float("nan"),
        "delivery_dose_verdict": ("delivery" if n_del > n_dose else
                                  "dose" if n_dose > n_del else "tied")
                                 if n_exp else "undefined",
        "cells_to_flip": (abs(n_del - n_dose) // 2) + 1 if n_exp else None,
        "cells_to_flip_as_fraction_of_exposure_limited":
            (((abs(n_del - n_dose) // 2) + 1) / n_exp) if n_exp else None,
    }


def _pair(dp, name: str, verbose: bool = True) -> Dict:
    """Both arms for one compound, plus the per-cell churn between them."""
    arms = {}
    for exchange in (False, True):
        arm = _run(dp, exchange)
        arms["ON" if exchange else "OFF"] = arm
        if verbose:
            bd = arm["survival_breakdown"]
            print(f"  {name:<14} exchange {'ON ' if exchange else 'OFF'}  "
                  f"kill = {arm['kill_fraction']:.6f}  "
                  f"of survivors delivery {bd['delivery']:.4f} / dose "
                  f"{bd['dose']:.4f} / resistance {bd['resistance']:.4f} / "
                  f"chance {bd['chance']:.4f}  "
                  f"dominant = {arm['dominant_limiter']:<10} "
                  f"({arm['solve_s']:.0f}s)", flush=True)

    lab_off, lab_on = arms["OFF"]["labels"], arms["ON"]["labels"]
    changed = _churn_mask(lab_off, lab_on)
    trans = Counter(zip(lab_off[changed].tolist(), lab_on[changed].tolist()))
    rec = {
        "OFF": {k: v for k, v in arms["OFF"].items() if k != "labels"},
        "ON": {k: v for k, v in arms["ON"].items() if k != "labels"},
        "per_cell": {
            "n_cells": int(changed.size),
            "n_changed": int(changed.sum()),
            "fraction_changed": float(changed.mean()),
            "transitions": {f"{LABELS[u]}->{LABELS[v]}": int(c)
                            for (u, v), c in trans.items()},
        },
        "margin_OFF": _margin(arms["OFF"]),
        "margin_ON": _margin(arms["ON"]),
        "bbb_penetrance": float(dp.bbb_penetrance),
        "ic50_uM": float(dp.ic50),
        "diffusivity_um2_s": float(dp.diffusivity),
        "k_loss_s_inv": float(dp.binding) + float(dp.clearance),
        "degenerate_zero_source": _is_degenerate(dp),
        "arms_bit_identical": bool(np.array_equal(lab_off, lab_on)),
        "mean_auc_ratio_on_over_off":
            arms["ON"]["mean_auc"] / arms["OFF"]["mean_auc"]
            if arms["OFF"]["mean_auc"] > 0 else None,
    }
    rec["dominant_limiter_inverts"] = (
        arms["OFF"]["dominant_limiter"] != arms["ON"]["dominant_limiter"])
    v_off = rec["margin_OFF"]["delivery_dose_verdict"]
    v_on = rec["margin_ON"]["delivery_dose_verdict"]
    rec["delivery_dose_verdict_OFF"] = v_off
    rec["delivery_dose_verdict_ON"] = v_on
    rec["delivery_dose_verdict_inverts"] = (
        v_off in ("delivery", "dose") and v_on in ("delivery", "dose")
        and v_off != v_on)
    rec["delivery_dose_verdict_defined"] = (
        v_off in ("delivery", "dose") and v_on in ("delivery", "dose"))
    if verbose:
        top = ", ".join(f"{k}: {v:,}" for k, v in
                        sorted(rec["per_cell"]["transitions"].items(),
                               key=lambda kv: -kv[1])[:4])
        print(f"    -> {rec['per_cell']['n_changed']:,} of "
              f"{rec['per_cell']['n_cells']:,} cells change label "
              f"({100 * rec['per_cell']['fraction_changed']:.3f} %)"
              + (f": {top}" if top else ""))
        print(f"    -> dominant limiter (4-way) "
              f"{arms['OFF']['dominant_limiter']} -> "
              f"{arms['ON']['dominant_limiter']}"
              + ("  *** INVERTS ***" if rec["dominant_limiter_inverts"] else ""))
        print(f"    -> delivery-vs-dose verdict {v_off} -> {v_on}"
              + ("  *** INVERTS ***" if rec["delivery_dose_verdict_inverts"]
                 else "") + "  " + _closest_str(rec), flush=True)
    return rec


def _closest_str(rec: Dict) -> str:
    ds = [rec[k]["distance_from_even_share"] for k in ("margin_OFF", "margin_ON")]
    ds = [d for d in ds if isinstance(d, float) and math.isfinite(d)]
    cs = [rec[k]["cells_to_flip"] for k in ("margin_OFF", "margin_ON")]
    cs = [c for c in cs if c is not None]
    if not ds:
        return "[no exposure-limited survivors in one arm: verdict undefined]"
    return (f"[closest approach to even: {min(ds):.4f} from 0.5, "
            f"{min(cs):,} cells from the flip]")


def _churn_mask(lab_off: np.ndarray, lab_on: np.ndarray) -> np.ndarray:
    """Per-cell label churn BETWEEN the two arms. Its own function so the
    'diffed against the wrong arm' defect has somewhere to be planted."""
    return lab_off != lab_on


# --------------------------------------------------------------------------- #
# STAGE: reproduce — the gate
# --------------------------------------------------------------------------- #
def stage_reproduce(verbose: bool = True, tol: float = 1e-12) -> Dict:
    """Re-run §4.2's doxorubicin and daunorubicin through §4.2's own compile
    path and require the published delivery : dose splits back."""
    alpha = _assert_frozen("reproduce/top")
    import isomer_penetration_map as IPM  # noqa: E402  (heavy; only this stage)

    comp = IPM.compile_all(alpha, two_route=True)
    out: Dict = {"alpha": alpha, "config": _tissue_cfg(), "tol": tol,
                 "compounds": {}}
    if verbose:
        print("=" * 82)
        print("STAGE reproduce — does this harness reproduce "
              "results/permeability_split.md §4.2?")
        print("=" * 82)
        print(f"\n{_tissue_cfg()}")
        print("compile path: isomer_penetration_map.compile_all(two_route=True), "
              "shared IC50 — §4.2's own\n")

    for name, published in PUBLISHED_SPLITS.items():
        rec = _pair(comp[name]["drug_params"], name, verbose=verbose)
        deltas = {}
        for arm in ("OFF", "ON"):
            for key in ("delivery", "dose"):
                got = rec[arm]["survival_breakdown"][key]
                want = published[arm][key]
                deltas[f"{arm}.{key}"] = {"got": got, "published": want,
                                          "abs_delta": abs(got - want)}
        rec["reproduction"] = deltas
        out["compounds"][name] = rec
        worst = max(d["abs_delta"] for d in deltas.values())
        if verbose:
            print(f"    -> worst |delta| against the published §4.2 split: "
                  f"{worst:.3e}\n")
        if worst > tol:
            raise AssertionError(
                f"REPRODUCTION FAILED for {name}: worst |delta| against "
                f"results/permeability_split.md §4.2 is {worst:.6e} > {tol:.0e}. "
                f"{json.dumps({k: v['abs_delta'] for k, v in deltas.items()})}. "
                "This harness is not measuring the same thing the published "
                "inversion was measured on; the census is not run.")
    out["worst_abs_delta"] = max(
        d["abs_delta"] for rec in out["compounds"].values()
        for d in rec["reproduction"].values())
    _assert_frozen("reproduce/bottom")
    return out


# --------------------------------------------------------------------------- #
# STAGE: census — the whole shipped library
# --------------------------------------------------------------------------- #
def stage_census(verbose: bool = True,
                 names: Optional[Sequence[str]] = None) -> Dict:
    alpha = _assert_frozen("census/top")
    names = list(names) if names else shipped_library()
    out: Dict = {"alpha": alpha, "config": _tissue_cfg(),
                 "library_source": ("sorted(set(molecular.library.GBM_DRUGS) | "
                                    "{'palbociclib'}) — app/server/api.py::"
                                    "_known_drugs before user compounds"),
                 "n_compounds": len(names), "compounds": {}}
    if verbose:
        print("\n" + "=" * 82)
        print(f"STAGE census — {len(names)} shipped compounds, exchange OFF vs ON, "
              "NOTHING scaled")
        print("=" * 82)
        print(f"\n{_tissue_cfg()}")
        print("compile path: calibrated_drugparams(name, alpha) -> "
              "apply_barrier_to_drugparams  (the product's own)\n")

    for name in names:
        dp = compile_shipped(name, alpha)
        if verbose:
            print(f"{name}  (ic50 = {dp.ic50:g} uM, D = {dp.diffusivity:.4g} "
                  f"um^2/s, bbb_penetrance = {dp.bbb_penetrance:.6f})", flush=True)
        out["compounds"][name] = _pair(dp, name, verbose=verbose)

    out.update(_summarise(out["compounds"]))
    _assert_frozen("census/bottom")
    return out


#: A verdict this close to even is reported here as UNDECIDED rather than as a
#: finding. It is the doxorubicin row of ``results/permeability_split.md`` §4.2
#: — 0.497 : 0.503 — rounded OUT to a round number, i.e. the loosest band that
#: still contains the published near-miss. It is a reporting threshold for this
#: file only: no simulation parameter depends on it, no compound was inspected
#: before it was written down, and both the raw share and the raw cell margin
#: are tabulated beside every verdict so a reader can pick their own.
NEAR_EVEN_BAND = 0.05


def _summarise(compounds: Dict[str, Dict]) -> Dict:
    """The headline. Degenerate zero-source compounds are named and EXCLUDED
    from the denominator rather than allowed to inflate the stable count.

    TWO rates are reported because the panel makes two different claims:
    ``delivery_dose_inversion_rate`` is the product's headline discrimination
    (delivery vs dose among the exposure-limited survivors), and
    ``dominant_limiter_inversion_rate`` is the four-way limiter the report
    prints. A compound can be stable on one and not the other.
    """
    degenerate = sorted(n for n, r in compounds.items()
                        if r["degenerate_zero_source"])
    scored = sorted(n for n in compounds if n not in degenerate)
    dom_inv = sorted(n for n in scored
                     if compounds[n]["dominant_limiter_inverts"])
    dd_scored = sorted(n for n in scored
                       if compounds[n]["delivery_dose_verdict_defined"])
    dd_undefined = sorted(n for n in scored if n not in dd_scored)
    dd_inv = sorted(n for n in dd_scored
                    if compounds[n]["delivery_dose_verdict_inverts"])

    def closest(n: str) -> float:
        ds = [compounds[n][k]["distance_from_even_share"]
              for k in ("margin_OFF", "margin_ON")]
        ds = [d for d in ds if isinstance(d, float) and math.isfinite(d)]
        return min(ds) if ds else float("nan")

    def closest_cells(n: str):
        cs = [compounds[n][k]["cells_to_flip"]
              for k in ("margin_OFF", "margin_ON")]
        cs = [c for c in cs if c is not None]
        return min(cs) if cs else None

    near = {n: {"closest_distance_from_even": closest(n),
                "closest_cells_to_flip": closest_cells(n),
                "inverts": compounds[n]["delivery_dose_verdict_inverts"]}
            for n in dd_scored}
    near_miss = sorted(n for n in dd_scored
                       if not compounds[n]["delivery_dose_verdict_inverts"]
                       and closest(n) <= NEAR_EVEN_BAND)
    not_robust = sorted(set(dd_inv) | set(near_miss))
    return {
        "degenerate_zero_source": degenerate,
        "scored_compounds": scored,
        "n_scored": len(scored),
        # --- four-way dominant limiter -------------------------------------
        "dominant_limiter_inverting_compounds": dom_inv,
        "n_dominant_limiter_inverting": len(dom_inv),
        "dominant_limiter_inversion_rate":
            len(dom_inv) / len(scored) if scored else None,
        # --- the product's headline discrimination -------------------------
        "delivery_dose_scored_compounds": dd_scored,
        "delivery_dose_undefined_compounds": dd_undefined,
        "n_delivery_dose_scored": len(dd_scored),
        "delivery_dose_inverting_compounds": dd_inv,
        "n_delivery_dose_inverting": len(dd_inv),
        "delivery_dose_inversion_rate":
            len(dd_inv) / len(dd_scored) if dd_scored else None,
        # --- inverting OR one coin-flip away -------------------------------
        "near_even_band": NEAR_EVEN_BAND,
        "near_miss_compounds": near_miss,
        "not_robust_compounds": not_robust,
        "n_not_robust": len(not_robust),
        "not_robust_rate": len(not_robust) / len(dd_scored) if dd_scored else None,
        "near_even": near,
        "panel_churn": {
            n: compounds[n]["per_cell"]["fraction_changed"] for n in compounds},
        "max_churn_fraction": max(
            (compounds[n]["per_cell"]["fraction_changed"] for n in compounds),
            default=0.0),
        "total_cells_changed": sum(
            compounds[n]["per_cell"]["n_changed"] for n in compounds),
    }


# --------------------------------------------------------------------------- #
# THE CHECKS
# --------------------------------------------------------------------------- #
class CheckFailure(AssertionError):
    pass


def check_split_is_not_inert(compounds: Dict[str, Dict]) -> str:
    """A. The split adds a back-flux ``-PS*C`` to a source that was ``+PS*Cp``.
    With ``C >= 0`` the ON arm's exposure is STRICTLY below the OFF arm's
    wherever any drug arrived at all, and equal ONLY where the source is
    identically zero. An exchange flag dropped between PDEConfig and the solver
    makes the two arms bit-identical and this file would then report a 0 %
    inversion rate with every other stage green. That is the vacuous-green trap
    here, so it is checked as a PROPERTY of the law rather than as an
    arm-versus-arm difference that an inert flag satisfies trivially."""
    bad = []
    for n, r in compounds.items():
        off, on = r["OFF"]["mean_auc"], r["ON"]["mean_auc"]
        # Degeneracy is re-derived from the compiled value here, NOT read off
        # the record's own flag, so this check stays independent of the flag
        # check D owns. Two checks that share a predicate fail together and
        # neither of them is then covering its own subject.
        if float(r["bbb_penetrance"]) == 0.0:
            if not (off == 0.0 and on == 0.0):
                bad.append(f"{n}: zero source but mean AUC {off!r} / {on!r}")
        elif not (on < off):
            bad.append(f"{n}: split did not reduce exposure — "
                       f"mean AUC OFF {off!r}, ON {on!r} "
                       f"(ratio {on / off if off else float('nan')!r})")
    if bad:
        raise CheckFailure(
            "the exchange flag is INERT or the split does not saturate: "
            + "; ".join(bad))
    n_live = sum(1 for r in compounds.values()
                 if float(r["bbb_penetrance"]) != 0.0)
    return (f"A ok — the split strictly reduced exposure on all {n_live} "
            f"live-source compounds; an inert flag would have tied them.")


def check_reproduces_published(repro: Dict, tol: float = 1e-12) -> str:
    """B. The §4.2 gate, re-stated as a check so ``--check`` can rerun it."""
    worst = repro.get("worst_abs_delta")
    if worst is None or worst > tol:
        raise CheckFailure(
            f"harness does not reproduce results/permeability_split.md §4.2: "
            f"worst |delta| {worst!r} > {tol:.0e}")
    return (f"B ok — reproduces §4.2's doxorubicin and daunorubicin "
            f"delivery:dose splits to {worst:.3e}.")


def check_churn_is_measured(compounds: Dict[str, Dict]) -> str:
    """C. The churn must be counted BETWEEN the arms and must reconcile with
    the transition matrix. A diff taken against the wrong arm returns 0 changed
    cells for everything, which reads as 'the split changes nothing' and is
    green."""
    bad = []
    total = 0
    for n, r in compounds.items():
        pc = r["per_cell"]
        tsum = sum(pc["transitions"].values())
        if tsum != pc["n_changed"]:
            bad.append(f"{n}: transitions sum to {tsum} but n_changed is "
                       f"{pc['n_changed']}")
        if any(t.split("->")[0] == t.split("->")[1] for t in pc["transitions"]):
            bad.append(f"{n}: a transition has the same source and target label")
        if r["arms_bit_identical"] != (pc["n_changed"] == 0):
            bad.append(f"{n}: arms_bit_identical={r['arms_bit_identical']} "
                       f"disagrees with n_changed={pc['n_changed']}")
        if float(r["bbb_penetrance"]) == 0.0 and pc["n_changed"] != 0:
            bad.append(f"{n}: zero source but {pc['n_changed']} cells changed")
        total += pc["n_changed"]
    if total <= 0:
        bad.append("no cell anywhere on the panel changed label between the "
                   "arms — the churn is being diffed against the wrong arm, or "
                   "the flag is inert")
    if bad:
        raise CheckFailure("per-cell churn is not being measured: "
                           + "; ".join(bad))
    return (f"C ok — {total:,} label changes counted between the arms across "
            f"the panel, each reconciling with its transition matrix.")


def check_degenerate_is_degenerate(compounds: Dict[str, Dict]) -> str:
    """D. A ``bbb_penetrance == 0`` compound has an identically zero source, so
    it must come back bit-identical across the arms with a delivery share of
    exactly 1.0 — its DELIVERY verdict is an input. If it is not detected as
    degenerate it silently inflates the count of compounds with a stable
    verdict, which is the number this file exists to report."""
    deg = [n for n, r in compounds.items() if r["degenerate_zero_source"]]
    bad = []
    for n, r in compounds.items():
        zero_src = float(r["bbb_penetrance"]) == 0.0
        if zero_src != r["degenerate_zero_source"]:
            bad.append(f"{n}: bbb_penetrance {r['bbb_penetrance']!r} but "
                       f"degenerate_zero_source={r['degenerate_zero_source']}")
        if not zero_src:
            continue
        if not r["arms_bit_identical"]:
            bad.append(f"{n}: zero source but the arms differ")
        for arm in ("OFF", "ON"):
            share = r[arm]["survival_breakdown"]["delivery"]
            if share != 1.0:
                bad.append(f"{n}/{arm}: zero source but delivery share is "
                           f"{share!r}, not exactly 1.0")
    if bad:
        raise CheckFailure("a zero-source compound is not being scored as the "
                           "degenerate case it is: " + "; ".join(bad))
    return (f"D ok — {len(deg)} zero-source compound(s) {deg} detected by value, "
            "bit-identical across the arms, delivery share exactly 1.0, and "
            "excluded from the inversion denominator.")


def run_checks(census: Dict, repro: Optional[Dict] = None,
               verbose: bool = True) -> List[str]:
    comps = census["compounds"]
    lines = [check_split_is_not_inert(comps)]
    if repro is not None:
        lines.append(check_reproduces_published(repro))
    lines.append(check_churn_is_measured(comps))
    lines.append(check_degenerate_is_degenerate(comps))
    if verbose:
        print("\n" + "=" * 82)
        print("CHECKS")
        print("=" * 82)
        for ln in lines:
            print("  " + ln)
    return lines


# --------------------------------------------------------------------------- #
def _report(census: Dict) -> None:
    c = census["compounds"]
    print("\n" + "=" * 82)
    print("HEADLINE")
    print("=" * 82)
    print(f"\n{census['config']}\n")
    print("  delivery share = SURVIVED_UNREACHED / (SURVIVED_UNREACHED + "
          "SURVIVED_UNDERDOSED),\n  i.e. the product's own delivery-vs-dose "
          "question asked of the cells it is posed about.\n")
    print(f"  {'compound':<14}{'bbb':>9}{'OFF share':>11}{'ON share':>11}"
          f"{'delivery-vs-dose':>26}{'4-way limiter':>28}{'cells to flip':>15}")
    for n in sorted(c):
        r = c[n]
        so = r["margin_OFF"]["delivery_share_of_exposure_limited"]
        sn = r["margin_ON"]["delivery_share_of_exposure_limited"]
        if r["degenerate_zero_source"]:
            verdict = "DEGENERATE (input)"
        elif not r["delivery_dose_verdict_defined"]:
            verdict = "UNDEFINED (no exp-limited)"
        elif r["delivery_dose_verdict_inverts"]:
            verdict = (f"INVERTS {r['delivery_dose_verdict_OFF']}->"
                       f"{r['delivery_dose_verdict_ON']}")
        else:
            verdict = f"stable ({r['delivery_dose_verdict_OFF']})"
        dom = (f"INVERTS {r['OFF']['dominant_limiter']}->"
               f"{r['ON']['dominant_limiter']}" if r["dominant_limiter_inverts"]
               else f"stable ({r['OFF']['dominant_limiter']})")
        cs = [r[k]["cells_to_flip"] for k in ("margin_OFF", "margin_ON")]
        cs = [x for x in cs if x is not None]
        ctf = f"{min(cs):,}" if cs else "-"
        print(f"  {n:<14}{r['bbb_penetrance']:>9.5f}"
              f"{so:>11.4f}{sn:>11.4f}{verdict:>26}{dom:>28}{ctf:>15}")

    print(f"\n  DELIVERY-vs-DOSE INVERSION RATE: "
          f"{census['n_delivery_dose_inverting']} of "
          f"{census['n_delivery_dose_scored']} scored compounds = "
          f"{100 * (census['delivery_dose_inversion_rate'] or 0.0):.1f} %  "
          f"({census['delivery_dose_inverting_compounds']})")
    print(f"  4-WAY DOMINANT-LIMITER INVERSION RATE: "
          f"{census['n_dominant_limiter_inverting']} of {census['n_scored']} = "
          f"{100 * (census['dominant_limiter_inversion_rate'] or 0.0):.1f} %  "
          f"({census['dominant_limiter_inverting_compounds']})")
    print(f"  INVERTS OR WITHIN {census['near_even_band']:g} OF EVEN: "
          f"{census['n_not_robust']} of {census['n_delivery_dose_scored']} = "
          f"{100 * (census['not_robust_rate'] or 0.0):.1f} %  "
          f"({census['not_robust_compounds']})")
    print(f"  excluded as degenerate zero-source (verdict is an INPUT): "
          f"{census['degenerate_zero_source']}")
    print(f"  excluded as undefined (no exposure-limited survivors in one arm): "
          f"{census['delivery_dose_undefined_compounds']}")
    near = sorted(census["near_even"].items(),
                  key=lambda kv: kv[1]["closest_distance_from_even"])
    print("\n  closest approach to an even delivery:dose split, per compound "
          "(either arm):")
    for n, m in near:
        print(f"    {n:<14}{m['closest_distance_from_even']:>10.4f} from 0.5, "
              f"{m['closest_cells_to_flip']:>8,} cells from the flip"
              + ("   [INVERTS]" if m["inverts"] else ""))


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return str(o)
    return o


def main(argv: Optional[List[str]] = None) -> Dict:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", default="all",
                    choices=("all", "reproduce", "census"))
    ap.add_argument("--drug", action="append", default=None,
                    help="census: restrict to these compounds")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)

    t0 = time.time()
    res: Dict = {
        "alpha_asserted": ALPHA_LIVE,
        "shipped_defaults": {
            "PDEConfig.transvascular_exchange":
                PDEConfig().transvascular_exchange,
            "transvascular.SPLIT_DEFAULT": SPLIT_DEFAULT,
            "drugs.PS_SCALE_S_INV": PS_SCALE_S_INV,
        },
        "config": _tissue_cfg(),
    }
    repro = None
    if args.stage in ("all", "reproduce"):
        repro = stage_reproduce()
        res["reproduce"] = repro
    if args.stage in ("all", "census"):
        res["census"] = stage_census(names=args.drug)
        res["checks"] = run_checks(res["census"], repro)
        _report(res["census"])
    res["microenvironment_build_s"] = _MICRO_CACHE.get("build_s")
    res["runtime_s"] = time.time() - t0
    _assert_frozen("main/bottom")

    if not args.no_write and args.stage == "all" and not args.drug:
        with open(_DATA_PATH, "w") as fh:
            json.dump(_jsonable(res), fh, indent=1)
        print(f"\nwrote {os.path.relpath(_DATA_PATH, _ROOT)} "
              f"({res['runtime_s']:.0f}s wall clock)")
    return res


if __name__ == "__main__":
    main()
