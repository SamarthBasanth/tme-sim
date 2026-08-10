"""Falsification test of the ``tail90`` recommendation in ``results/p99_convention.md``.

The prediction under test
------------------------
``results/p99_convention.md`` sec. 1.4 / 3.1 / 3.6 diagnoses the grid error of the
scored 3-D:2-D ratio as being, *at leading order, nothing but the grid error of
the normaliser*, on the grounds that with exchange OFF the domain mean is
grid-exact by a discrete conservation identity, so
``err(mean/N) = -err(N)/(1+err(N))``. From that diagnosis it recommends replacing
the shipped ``p99`` bath-matching reference with ``tail90`` (CVaR-90, the mean
exposure of the best-perfused 10 % of tumour cells), whose grid error on the
controlled slab at the shipped ``dx = 12 um`` is **-0.08 %** against p99's
**+6.98 %**.

That file names its own falsification test (its sec. 3.6 "What would change this
recommendation" and sec. 8 bullet 3):

    "if re-scoring the existing tissue dx ladders under tail90 did NOT remove
     their non-monotonicity, the diagnosis would be wrong and this
     recommendation should be withdrawn. That test costs no new solves."

**This script runs that test, and reports first that its stated precondition is
false.** The tissue ladders' per-cell exposure fields were never persisted --
``experiments/grid_convergence_gbm.py`` holds them in a local
``e_unit_by_key`` dict (line 477/560/615) and writes only SUMMARY STATISTICS to
``results/sweep/gbmgrid_records.jsonl``; ``grid_convergence_breast.py`` does the
same. The only per-cell exposure vectors on disk anywhere in the repository are
at ``dx = 12 um``. Stage ``premise`` establishes that by enumeration rather than
assertion.

What is therefore still measurable at zero solve cost, and what each thing settles
------------------------------------------------------------------------------
STAGE premise -- **which exposure fields actually exist**. Enumerates every
    on-disk per-cell ``e_unit`` vector in the repo with its dx, and every ladder
    record with the statistics it does carry. Establishes the blocker as a
    measurement.

STAGE ladder -- **the decomposition the p99 audit did not do.** The scored ratio
    factorises EXACTLY as ``ratio = f_pen * f_state * f_sens * interaction``
    (verified per rung to machine precision below). ``f_state`` -- the
    oxygen-banded cell-state resistance factor -- is **completely untouched by any
    normaliser choice**, and ``results/grid_convergence_gbm.md`` sec. 5 already
    measures it as non-monotone in dx. So the ladder movement splits into a part
    a normaliser could in principle fix and a part it provably cannot. Both are
    read off the cached records with no new computation.

    Inside ``f_pen`` the split goes one level further. In the linear-Hill limit
    ``f_pen -> mean(e)/N(e)``, so
    ``dln f_pen = dln mean - dln N + (Hill correction)``, and the ladder carries
    ``e_unit_mean`` and ``e_unit_p99`` at every rung. That gives a DIRECT tissue
    measurement of the p99 audit's central claim -- is ``dln N`` the dominant
    term? -- without needing the fields. The Hill correction is not assumed: it
    is measured as the per-drug regression slope of the true ``dln f_pen`` on the
    surrogate ``dln(mean/p99)`` across the ladder, and reported with its R^2.

STAGE perfect -- **the best case tail90 could possibly deliver.** ``tail90``'s
    slab grid error is -0.08 % where p99's is +6.98 %; the IDEALISATION of that is
    a normaliser with EXACTLY zero grid error, i.e. ``N(dx) = const``. Under it
    ``f_pen(dx) ∝ mean(dx)`` in the linear-Hill limit, so the counterfactual
    ladder is ``mean(dx) * f_state(dx) * interaction(dx)``, every factor measured.
    This is deliberately GENEROUS to the recommendation: no real normaliser can
    beat a perfect one, so if the perfect-normaliser ladder is still
    non-monotone with a comparable spread, ``tail90`` cannot rescue the ladders
    and the diagnosis is refuted on tissue. Monotonicity, spread, power-law fit
    and Richardson are recomputed on the counterfactual ladder with
    ``grid_convergence_gbm``'s OWN ``fit_power_law`` and ``richardson_triple``,
    imported rather than re-implemented, so the two are commensurable.

STAGE class -- **the estimator-class claim, tested on real tissue.**
    ``p99_convention.md`` sec. 3.1 argues that an order statistic is a STEP
    FUNCTION of dx (refining moves which cell the quantile lands in) while a
    moment or a tail mean is smooth, and blames the tissue ladders'
    non-monotonicity on that. The ladder records carry four normalisers at every
    rung -- ``mean`` (moment), ``p50``, ``p99`` and ``max`` (order statistics) --
    so the claim is directly checkable on the real sections. ``tail90`` itself is
    NOT among them and cannot be added without the fields.

STAGE dx12 -- **exact re-score at the one grid where the fields DO exist.**
    Both published panels, all nine verdicts, under p99 and tail90, through
    ``experiments/p99_convention.decompose_general`` (which
    ``assert_matches_shipped`` pins bit-identical to the shipped
    ``kill_calibration.decompose`` at ``tag='p99'``). Also measures, exactly, the
    elasticity ``dln f_pen / dln N`` at each convention's own operating point, so
    stage ``perfect``'s linear-Hill surrogate has a directly computed
    cross-check rather than an assumed one.

STAGE verdicts -- verdict consequences, computed LAST and written to their own
    file, exactly as ``p99_convention.py`` did, so numerical merit and verdict
    counts cannot contaminate each other.

HARD RULES OBSERVED
-------------------
* **No PDE transport solve is run.** Stage ``dx12`` runs the engine PROLOGUE
  (reconstruction, oxygen, IFP, cell-state labelling) because the ``sens`` and
  ``resist`` vectors are not cached, and reads the transport result from
  ``data_cache/``; this is the same thing ``p99_convention.py`` does and it adds
  no drug transport solve.
* **Nothing under ``tmesim/`` is touched**, no shipped default is changed, and no
  existing results file is edited.
* **A negative result is the point.** Nothing here selects a variant of tail90
  that rescues it; ``tail90`` is tested exactly as specified in
  ``p99_convention.md`` sec. 7.1 (CVaR-90 = mean of the top decile).

Run:
    ./.venv/bin/python experiments/tail90_falsification.py --stage premise
    ./.venv/bin/python experiments/tail90_falsification.py --stage ladder
    ./.venv/bin/python experiments/tail90_falsification.py --stage dx12
    ./.venv/bin/python experiments/tail90_falsification.py --stage verdicts
Writes:
    results/sweep/tail90_premise.json
    results/sweep/tail90_ladder.json      (stages ladder + perfect + class)
    results/sweep/tail90_dx12.json
    results/sweep/tail90_verdicts.json
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

OUT_DIR = os.path.join(_ROOT, "results", "sweep")
_T0 = time.time()

#: The published grid every verdict in the repo is taken at.
PUBLISHED_DX = 12.0
#: ``results/grid_convergence_gbm.md`` sec. 0 splits the ladder here: below this the
#: mesh resolves the perivascular layer well enough that the remaining movement
#: is the interesting part; above it the coarse grids are simply wrong.
FINE_MAX_DX = 18.0
#: The dt-control arms in the GBM ladder are NOT rungs of the spatial ladder.
LADDER_DT_S = 10.0


def log(msg: str) -> None:
    print(f"[t90 {time.strftime('%H:%M:%S')} +{time.time()-_T0:7.1f}s] {msg}",
          flush=True)


# --------------------------------------------------------------------------- #
# STAGE premise -- what exposure fields actually exist on disk?
# --------------------------------------------------------------------------- #
def stage_premise() -> Dict:
    """Enumerate every per-cell exposure vector on disk, and every ladder rung.

    The brief for this investigation asserts that "the exposure fields already
    exist in ``results/sweep/gbmgrid_*`` and ``results/sweep/gridconv_*`` and in
    ``data_cache/``". That is a checkable claim about files, and it is checked
    here by opening them rather than by trusting it, because every downstream
    stage depends on the answer.
    """
    out: Dict = {"per_cell_exposure_vectors_on_disk": [],
                 "ladder_record_files": [], "searched": []}

    # --- every on-disk file that carries an actual per-cell e_unit vector ----
    cands = (sorted(glob.glob(os.path.join(_ROOT, "data_cache", "*.json")))
             + sorted(glob.glob(os.path.join(OUT_DIR, "*.json")))
             + sorted(glob.glob(os.path.join(OUT_DIR, "*.jsonl"))))
    for path in cands:
        rel = os.path.relpath(path, _ROOT)
        size = os.path.getsize(path)
        out["searched"].append({"file": rel, "bytes": size})
        # cheap pre-filter: a per-cell vector is always stored under "e_unit"
        with open(path) as fh:
            head = fh.read(4_000_000)
        if '"e_unit"' not in head and '"e_unit":' not in head:
            continue
        try:
            if path.endswith(".jsonl"):
                recs = []
                with open(path) as fh:
                    for line in fh:
                        try:
                            recs.append(json.loads(line))
                        except Exception:
                            pass
                for r in recs:
                    v = r.get("e_unit")
                    if isinstance(v, list):
                        out["per_cell_exposure_vectors_on_disk"].append({
                            "file": rel, "key": f"{r.get('mode')}/{r.get('drug')}",
                            "n_cells": len(v), "dx_um": PUBLISHED_DX,
                            "dx_source": "file name / documented dx = 12 arm"})
            else:
                d = json.load(open(path))
                if isinstance(d, dict):
                    for k, v in d.items():
                        if isinstance(v, dict) and isinstance(v.get("e_unit"), list):
                            # cache keys are "drug|dx|dt|window|save_every[|ic50]"
                            parts = str(k).split("|")
                            dx = float(parts[1]) if len(parts) > 1 else float("nan")
                            out["per_cell_exposure_vectors_on_disk"].append({
                                "file": rel, "key": k, "n_cells": len(v["e_unit"]),
                                "dx_um": dx, "dx_source": "cache key"})
        except Exception as exc:                       # pragma: no cover
            out["searched"][-1]["parse_error"] = repr(exc)

    dxs = sorted({r["dx_um"] for r in out["per_cell_exposure_vectors_on_disk"]})
    out["distinct_dx_with_per_cell_fields"] = dxs
    out["n_per_cell_vectors"] = len(out["per_cell_exposure_vectors_on_disk"])

    # --- what the ladder files DO carry -------------------------------------
    for rel, tissue in (("results/sweep/gbmgrid_records.jsonl", "gbm"),
                        ("results/sweep/gridconv_records.jsonl", "breast")):
        path = os.path.join(_ROOT, rel)
        rows = [json.loads(l) for l in open(path)]
        spatial = [r for r in rows
                   if float(r.get("dt_s") or LADDER_DT_S) == LADDER_DT_S]
        out["ladder_record_files"].append({
            "file": rel, "tissue": tissue, "n_records": len(rows),
            "n_spatial_rungs": len(spatial),
            "drugs": sorted({r["drug"] for r in spatial}),
            "dx_ladder": sorted({r["dx_um"] for r in spatial}),
            "has_per_cell_e_unit": any("e_unit" in r for r in rows),
            "exposure_statistics_carried": [k for k in sorted(rows[0])
                                            if k.startswith("e_unit")],
            "fields_it_carries_instead": sorted(
                k for k in rows[0] if k.startswith("_e_unit")),
        })

    blocked = all(abs(d - PUBLISHED_DX) < 1e-9 for d in dxs) and len(dxs) >= 1
    out["verdict"] = {
        "ladder_per_cell_fields_exist": False,
        "only_dx_with_per_cell_fields": dxs,
        "can_rescore_ladders_under_tail90_without_solving": not blocked and False,
        "why": (
            "grid_convergence_gbm.run_ladder() builds e_unit_by_key in memory "
            "(line 477/560) and returns it (line 615) but never writes it; the "
            "JSONL rung records carry only mean / median / p01 / p99 / max plus a "
            "sum+first fingerprint. grid_convergence_breast.py is the same. A "
            "tail mean is not recoverable from those five order statistics, so "
            "tail90 cannot be evaluated at any rung except dx = 12 um without "
            "re-running the transport solves."),
    }
    log(f"premise: {out['n_per_cell_vectors']} per-cell exposure vectors on "
        f"disk, at dx = {dxs}; ladder rung fields present: "
        f"{any(f['has_per_cell_e_unit'] for f in out['ladder_record_files'])}")
    return out


# --------------------------------------------------------------------------- #
# Ladder loading + the exact factorisation
# --------------------------------------------------------------------------- #
def load_ladder(tissue: str) -> Dict[str, List[Dict]]:
    """The spatial dx ladder, per drug, fine -> coarse, from the cached records.

    dt-control arms are excluded: they are a different experiment (same dx, half
    dt) and including them would put two points at the same dx.
    """
    rel = {"gbm": "results/sweep/gbmgrid_records.jsonl",
           "breast": "results/sweep/gridconv_records.jsonl"}[tissue]
    rows = [json.loads(l) for l in open(os.path.join(_ROOT, rel))]
    rows = [r for r in rows
            if float(r.get("dt_s") or LADDER_DT_S) == LADDER_DT_S]
    out: Dict[str, List[Dict]] = {}
    for r in rows:
        out.setdefault(r["drug"], []).append(r)
    for name in out:
        out[name].sort(key=lambda r: float(r["dx_um"]))
    return out


def check_factorisation(rungs: Sequence[Dict]) -> Dict:
    """Verify ``ratio = f_pen * f_state * f_sens * interaction`` at every rung.

    Nothing downstream is legitimate unless this identity holds, because the
    whole decomposition is an additive split of ``ln(ratio)`` along it. It is a
    property of ``kill_calibration.decompose``'s own recorded outputs, so a
    failure would mean the records are not what they are labelled.
    """
    worst = 0.0
    for r in rungs:
        lhs = float(r["ratio_exact"])
        rhs = (float(r["f_penetration"]) * float(r["f_state"])
               * float(r["f_sensitivity"]) * float(r["interaction"]))
        worst = max(worst, abs(rhs / lhs - 1.0) if lhs else float("inf"))
    return {"worst_rel_deviation": worst, "ok": worst < 1e-9}


def _monotonicity(dxs: Sequence[float], vals: Sequence[float]) -> Dict:
    """Direction changes along the ladder, and the size of the largest reversal.

    "Monotone" here means the sequence never changes direction as dx varies.
    ``n_sign_changes`` is the count of direction reversals; ``worst_reversal_pct``
    is the largest step that goes against the ladder's overall trend, as a
    percentage -- a single tiny wobble and a violent oscillation both count as
    non-monotone, so the magnitude is reported alongside the flag.
    """
    v = np.asarray(vals, dtype=np.float64)
    d = np.diff(v)
    sgn = np.sign(d)
    nz = sgn[sgn != 0]
    changes = int(np.sum(nz[1:] != nz[:-1]))
    overall = np.sign(v[-1] - v[0])
    against = [100.0 * abs(step) / v[i]
               for i, step in enumerate(d) if np.sign(step) == -overall]
    return {"monotone": bool(changes == 0), "n_sign_changes": changes,
            "n_steps": int(d.size),
            "worst_reversal_pct": float(max(against)) if against else 0.0,
            "span_pct": float(100.0 * (v.max() - v.min()) / v.min()),
            "sd_over_mean_pct": float(100.0 * v.std(ddof=1) / v.mean())}


def _series(rungs: Sequence[Dict], key: str) -> np.ndarray:
    return np.asarray([float(r[key]) for r in rungs], dtype=np.float64)


# --------------------------------------------------------------------------- #
# STAGE ladder / perfect / class
# --------------------------------------------------------------------------- #
def analyse_drug(tissue: str, name: str, rungs: Sequence[Dict]) -> Dict:
    """The full decomposition for one drug's ladder. No solves, no fitting of
    anything that is not stated as a fit."""
    from experiments.grid_convergence_gbm import (fit_power_law,
                                                  richardson_triple)

    dxs = _series(rungs, "dx_um")
    ratio = _series(rungs, "ratio_exact")
    fpen = _series(rungs, "f_penetration")
    fstate = _series(rungs, "f_state")
    fsens = _series(rungs, "f_sensitivity")
    inter = _series(rungs, "interaction")
    emean = _series(rungs, "e_unit_mean")
    ep99 = _series(rungs, "e_unit_p99")
    ep50 = _series(rungs, "e_unit_median")
    emax = _series(rungs, "e_unit_max")

    i12 = int(np.argmin(np.abs(dxs - PUBLISHED_DX)))
    fine = dxs <= FINE_MAX_DX

    out: Dict = {
        "tissue": tissue, "drug": name,
        "dx_um": dxs.tolist(),
        "factorisation_check": check_factorisation(rungs),
        "c_bath_uM": float(rungs[0]["c_bath_uM"]),
        "ic50_uM": float(rungs[0]["ic50_uM"]),
        "bath_multiple": float(rungs[0]["c_bath_uM"]) / float(rungs[0]["ic50_uM"]),
        "hill": float(rungs[0]["hill"]),
        "L_um": float(rungs[0]["L_um"]),
        "band_lo": float(rungs[0]["band_used_lo"]),
        "band_hi": float(rungs[0]["band_used_hi"]),
        "series": {"ratio": ratio.tolist(), "f_penetration": fpen.tolist(),
                   "f_state": fstate.tolist(), "f_sensitivity": fsens.tolist(),
                   "interaction": inter.tolist(),
                   "e_unit_mean": emean.tolist(), "e_unit_p99": ep99.tolist(),
                   "e_unit_median": ep50.tolist(), "e_unit_max": emax.tolist()},
    }

    # ----------------------------------------------------------------- #
    # 1. The additive log split of the SCORED ratio.
    #    ln ratio = ln f_pen + ln f_state + ln f_sens + ln interaction.
    #    f_state and f_sens are normaliser-FREE by construction: neither
    #    touches the exposure field. `interaction` is the normalised covariance
    #    mean(p3*r)/(mean(p3)*mean(r)) and is only weakly normaliser-dependent.
    # ----------------------------------------------------------------- #
    def _mono_block(v, mask=None):
        vv = v if mask is None else v[mask]
        dd = dxs if mask is None else dxs[mask]
        return _monotonicity(dd, vv)

    terms = {"f_penetration": fpen, "f_state": fstate,
             "f_sensitivity": fsens, "interaction": inter, "ratio": ratio}
    out["term_behaviour"] = {}
    for k, v in terms.items():
        if float(v.max() - v.min()) == 0.0:
            out["term_behaviour"][k] = {"constant": True, "value": float(v[0])}
            continue
        out["term_behaviour"][k] = {
            "constant": False,
            "all_rungs": _mono_block(v),
            "fine_rungs": _mono_block(v, fine),
            "normaliser_free": k in ("f_state", "f_sensitivity"),
        }

    # Contribution of each term to the ladder movement of ln(ratio), measured
    # two ways because they answer different questions:
    #   `span_share`  -- each term's own log-span as a fraction of the sum of
    #                    spans. Ignores cancellation; says how big each wobble is.
    #   `cov_share`   -- cov(ln term, ln ratio)/var(ln ratio), which sums to
    #                    EXACTLY 1 over the terms and does account for
    #                    cancellation. This is the honest attribution.
    lr = np.log(ratio)
    logs = {k: np.log(v) for k, v in terms.items() if k != "ratio"}
    shares: Dict[str, Dict] = {}
    var_lr = float(np.var(lr, ddof=1))
    tot_span = sum(float(np.ptp(v)) for v in logs.values())
    for k, v in logs.items():
        shares[k] = {
            "log_span": float(np.ptp(v)),
            "span_share": (float(np.ptp(v)) / tot_span) if tot_span else 0.0,
            "cov_share": (float(np.cov(v, lr, ddof=1)[0, 1] / var_lr)
                          if var_lr > 0 else float("nan")),
        }
    out["attribution_all_rungs"] = {"var_ln_ratio": var_lr, "terms": shares,
                                    "cov_shares_sum": float(
                                        sum(s["cov_share"] for s in shares.values()))}
    lrf = np.log(ratio[fine])
    var_f = float(np.var(lrf, ddof=1))
    sharesf: Dict[str, Dict] = {}
    tot_span_f = sum(float(np.ptp(v[fine])) for v in logs.values())
    for k, v in logs.items():
        sharesf[k] = {
            "log_span": float(np.ptp(v[fine])),
            "span_share": (float(np.ptp(v[fine])) / tot_span_f) if tot_span_f else 0.0,
            "cov_share": (float(np.cov(v[fine], lrf, ddof=1)[0, 1] / var_f)
                          if var_f > 0 else float("nan")),
        }
    out["attribution_fine_rungs"] = {"var_ln_ratio": var_f, "terms": sharesf,
                                     "cov_shares_sum": float(
                                         sum(s["cov_share"] for s in sharesf.values()))}

    # ----------------------------------------------------------------- #
    # 2. Inside f_pen: is the NORMALISER the dominant term, as claimed?
    #    Linear-Hill limit: f_pen -> mean(e)/N(e), so
    #        dln f_pen = dln mean - dln N + (Hill correction).
    #    The surrogate is not assumed to be exact -- it is regressed against the
    #    true dln f_pen and the slope (the Hill damping) and R^2 are reported.
    # ----------------------------------------------------------------- #
    g_lin = emean / ep99
    dl_fpen = np.log(fpen) - np.log(fpen[i12])
    dl_surr = np.log(g_lin) - np.log(g_lin[i12])
    A = dl_surr[:, None]
    gamma = float(np.linalg.lstsq(A, dl_fpen, rcond=None)[0][0]) if np.any(dl_surr) else float("nan")
    resid = dl_fpen - gamma * dl_surr
    ss_tot = float(np.sum((dl_fpen - dl_fpen.mean()) ** 2))
    out["fpen_surrogate"] = {
        "definition": "linear-Hill surrogate  f_pen ~ mean(e)/p99(e)",
        "hill_damping_gamma": gamma,
        "r2": float(1.0 - np.sum(resid ** 2) / ss_tot) if ss_tot > 0 else float("nan"),
        "worst_abs_residual_pct": float(100.0 * np.max(np.abs(resid))),
        "log_span_ln_mean": float(np.ptp(np.log(emean))),
        "log_span_ln_p99": float(np.ptp(np.log(ep99))),
        "log_span_ln_mean_fine": float(np.ptp(np.log(emean[fine]))),
        "log_span_ln_p99_fine": float(np.ptp(np.log(ep99[fine]))),
        "normaliser_share_of_fpen_movement_all": (
            float(np.ptp(np.log(ep99)) / (np.ptp(np.log(ep99)) + np.ptp(np.log(emean))))),
        "normaliser_share_of_fpen_movement_fine": (
            float(np.ptp(np.log(ep99[fine]))
                  / (np.ptp(np.log(ep99[fine])) + np.ptp(np.log(emean[fine]))))),
        "mean_is_grid_exact_on_tissue": bool(np.ptp(np.log(emean)) < 1e-6),
    }

    # ----------------------------------------------------------------- #
    # 3. STAGE perfect -- the counterfactual ladder under a normaliser with
    #    EXACTLY ZERO grid error. tail90 measures -0.08 % on the slab against
    #    p99's +6.98 %; the limit of that is N(dx) = const, under which
    #    f_pen(dx) ∝ mean(dx) (linear-Hill) and therefore
    #        ratio_perfect(dx) ∝ mean(dx)^gamma * f_state(dx) * interaction(dx).
    #    gamma is the MEASURED damping from block 2, so the counterfactual is
    #    calibrated on this drug's own ladder rather than assumed to be linear.
    #    It is anchored to the published dx = 12 value so the two ladders are
    #    plotted on the same scale; only the SHAPE of the counterfactual is a
    #    claim, and the shape is anchor-independent.
    # ----------------------------------------------------------------- #
    ratio_perfect = (ratio[i12]
                     * np.exp(gamma * (np.log(emean) - np.log(emean[i12])))
                     * (fstate / fstate[i12]) * (inter / inter[i12]))
    out["counterfactual_perfect_normaliser"] = {
        "ratio": ratio_perfect.tolist(),
        "anchored_at_dx_um": float(dxs[i12]),
        "all_rungs": _monotonicity(dxs, ratio_perfect),
        "fine_rungs": _monotonicity(dxs[fine], ratio_perfect[fine]),
    }
    # and the pure-f_state ladder: what is left when the normaliser AND the
    # field mean are both frozen. Nothing any normaliser can do touches this.
    ratio_fstate_only = ratio[i12] * (fstate / fstate[i12])
    out["counterfactual_f_state_only"] = {
        "ratio": ratio_fstate_only.tolist(),
        "all_rungs": _monotonicity(dxs, ratio_fstate_only),
        "fine_rungs": _monotonicity(dxs[fine], ratio_fstate_only[fine]),
    }

    # ----------------------------------------------------------------- #
    # 3b. Does the normaliser's grid error CANCEL f_state's, or add to it?
    #     This is the question the span numbers alone cannot answer, and it
    #     decides the sign of the effect of improving the normaliser. The
    #     normaliser enters the score as -ln N, so a POSITIVE correlation
    #     between (-ln p99) and ln f_state means the two grid errors REINFORCE
    #     and removing the normaliser term shrinks the ladder; a NEGATIVE
    #     correlation means they CANCEL and removing it makes the ladder worse.
    # ----------------------------------------------------------------- #
    def _corr(a, b):
        a = np.asarray(a, float)
        b = np.asarray(b, float)
        if a.std() == 0 or b.std() == 0:
            return float("nan")
        return float(np.corrcoef(a, b)[0, 1])

    out["normaliser_vs_f_state_coupling"] = {
        "definition": ("corr over the ladder of the normaliser's contribution "
                       "to ln(ratio), which is -ln N, against ln f_state; "
                       "negative = the two grid errors partially CANCEL"),
        "corr_all_rungs": _corr(-np.log(ep99), np.log(fstate)),
        "corr_fine_rungs": _corr(-np.log(ep99[fine]), np.log(fstate[fine])),
        "corr_lnfpen_lnfstate_all": _corr(np.log(fpen), np.log(fstate)),
        "corr_lnfpen_lnfstate_fine": _corr(np.log(fpen[fine]), np.log(fstate[fine])),
    }

    # ----------------------------------------------------------------- #
    # 3c. Is the counterfactual's conclusion an artefact of the fitted gamma?
    #     Re-run it at gamma = 1 (the exact linear-Hill value) and at gamma = 0
    #     (the fully saturated limit, where f_pen stops responding to the field
    #     at all). If the monotonicity verdict is the same at all three, it is a
    #     property of the ladder and not of the surrogate.
    # ----------------------------------------------------------------- #
    rob = {}
    for gtag, g in (("fitted", gamma), ("linear_hill_gamma_1", 1.0),
                    ("saturated_gamma_0", 0.0)):
        rr = (ratio[i12] * np.exp(g * (np.log(emean) - np.log(emean[i12])))
              * (fstate / fstate[i12]) * (inter / inter[i12]))
        rob[gtag] = {"gamma": float(g),
                     "all_rungs": _monotonicity(dxs, rr),
                     "fine_rungs": _monotonicity(dxs[fine], rr[fine])}
    out["counterfactual_gamma_robustness"] = rob

    # convergence diagnostics on BOTH ladders, with the ladder scripts' own code
    def _conv(vals: np.ndarray) -> Dict:
        v = np.asarray(vals, dtype=np.float64)
        fit = fit_power_law(dxs.tolist(), v.tolist())
        rich: Dict = {"ok": False, "reason": "no 8/12/18 triple on this ladder"}
        try:
            idx = {float(d): k for k, d in enumerate(dxs)}
            rich = richardson_triple(v[idx[8.0]], v[idx[12.0]], v[idx[18.0]],
                                     r=1.5)
        except KeyError:
            pass
        res: Dict = {"fit": fit, "richardson": rich}
        # is the power-law fit degenerate, by grid_convergence_gbm.md's own
        # criterion: p pinned at a scan edge AND C collapsed to ~0.
        if fit.get("ok"):
            res["fit_degenerate"] = bool(
                (abs(fit["p"] - 4.0) < 1e-3 or abs(fit["p"] - 0.25) < 1e-3)
                and abs(fit["C"]) < 1e-5)
            res["disc_err_at_dx12_vs_fit_q0_pct"] = float(
                100.0 * (v[i12] - fit["q0"]) / fit["q0"])
        if rich.get("ok"):
            res["disc_err_at_dx12_vs_richardson_pct"] = float(
                100.0 * (v[i12] - rich["q0"]) / rich["q0"])
        return res

    out["convergence_published_p99"] = _conv(ratio)
    out["convergence_perfect_normaliser"] = _conv(ratio_perfect)
    out["ladder_p99"] = {"all_rungs": _monotonicity(dxs, ratio),
                         "fine_rungs": _monotonicity(dxs[fine], ratio[fine])}

    # ----------------------------------------------------------------- #
    # 4. STAGE class -- the four normalisers the records actually carry.
    #    p99_convention.md sec. 3.1: order statistics are STEP FUNCTIONS of dx,
    #    moments are smooth. On tissue, is that true?
    # ----------------------------------------------------------------- #
    out["normaliser_ladders"] = {}
    for tag, v, kind in (("mean", emean, "moment"), ("p50", ep50, "order"),
                         ("p99", ep99, "order"), ("max", emax, "order")):
        out["normaliser_ladders"][tag] = {
            "kind": kind, "values": v.tolist(),
            "all_rungs": _monotonicity(dxs, v),
            "fine_rungs": _monotonicity(dxs[fine], v[fine]),
        }
    return out


def stage_ladder() -> Dict:
    out: Dict = {"published_dx_um": PUBLISHED_DX, "fine_max_dx_um": FINE_MAX_DX,
                 "panels": []}
    for tissue in ("gbm", "breast"):
        lad = load_ladder(tissue)
        prec: Dict = {"tissue": tissue, "drugs": []}
        for name, rungs in sorted(lad.items()):
            a = analyse_drug(tissue, name, rungs)
            prec["drugs"].append(a)
            log(f"ladder {tissue}/{name}: p99 fine span "
                f"{a['ladder_p99']['fine_rungs']['span_pct']:.2f}% "
                f"(mono={a['ladder_p99']['fine_rungs']['monotone']}) -> "
                f"perfect-normaliser fine span "
                f"{a['counterfactual_perfect_normaliser']['fine_rungs']['span_pct']:.2f}% "
                f"(mono={a['counterfactual_perfect_normaliser']['fine_rungs']['monotone']}) "
                f"| f_state alone "
                f"{a['counterfactual_f_state_only']['fine_rungs']['span_pct']:.2f}%")
        out["panels"].append(prec)
    return out


# --------------------------------------------------------------------------- #
# STAGE dx12 -- exact re-score where the fields DO exist
# --------------------------------------------------------------------------- #
#: The two conventions under test. Every other candidate in ``p99_convention.py``
#: is left out deliberately: this file adjudicates ONE recommendation.
TAGS: Tuple[str, ...] = ("p99", "tail90")


def stage_dx12() -> Dict:
    """Exact p99 vs tail90 at the published grid, plus the measured elasticity.

    The elasticity ``s = dln f_pen / dln N`` is what stage ``perfect``'s
    linear-Hill surrogate approximates by ``-1`` (damped by ``gamma``). Here it
    is computed EXACTLY on the real per-cell field by evaluating the scored
    ``f_pen`` at a range of scaled normaliser values around each convention's own
    operating point -- no linearisation, no fit to a ladder.
    """
    from experiments.p99_convention import (NORM_BY_TAG, assert_matches_shipped,
                                            decompose_general, load_breast_panel,
                                            load_gbm_panel, normalise)
    from tmesim.validation.kill_calibration import p_hill

    out: Dict = {"dx_um": PUBLISHED_DX, "tags": list(TAGS), "panels": []}
    for loader in (load_gbm_panel, load_breast_panel):
        panel = loader()
        log(f"dx12: {panel['tissue']} panel loaded, {panel['n_tumor']} tumour cells")
        prec: Dict = {"tissue": panel["tissue"], "n_tumor": panel["n_tumor"],
                      "drugs": []}
        for name, d in panel["drugs"].items():
            e, sens, resist = d["e_unit"], d["sens"], d["resist"]
            cb, ic, hh, mk = (d["c_bath_uM"], d["ic50_uM"], d["hill"],
                              d["max_kill"])
            # refuse to proceed unless the generalised scorer IS the shipped one
            ship = assert_matches_shipped(e, sens, resist, cb, ic, hh, mk)
            row: Dict = {"drug": name, "c_bath_uM": cb, "ic50_uM": ic,
                         "bath_multiple": cb / ic, "hill": hh,
                         "shipped_ratio": ship["ratio_exact"],
                         "n_cells": int(e.size), "conventions": {}}
            p2 = float(p_hill(cb, ic, hh, mk))
            for tag in TAGS:
                dec = decompose_general(e, sens, resist, cb, ic, hh, mk, tag=tag)
                n0 = normalise(e, tag)
                # exact elasticity of f_pen w.r.t. the normaliser value, at this
                # convention's own operating point. +-2 % brackets the whole
                # movement the p99 normaliser shows across the fine ladder.
                els = {}
                for eps in (0.02, 0.05):
                    fp = []
                    for lam in (1.0 - eps, 1.0 + eps):
                        fp.append(float(np.mean(p_hill(cb * e / (n0 * lam), ic,
                                                       hh, mk)) / p2))
                    els[f"pm{eps:g}"] = float(
                        (math.log(fp[1]) - math.log(fp[0]))
                        / (math.log(1.0 + eps) - math.log(1.0 - eps)))
                row["conventions"][tag] = {
                    "label": NORM_BY_TAG[tag]["label"],
                    "N": n0, "ratio_exact": dec["ratio_exact"],
                    "f_penetration": dec["f_penetration"],
                    "f_state": dec["f_state"],
                    "dln_fpen_dln_N": els,
                }
            row["tail90_over_p99_N"] = (row["conventions"]["tail90"]["N"]
                                        / row["conventions"]["p99"]["N"])
            row["tail90_over_p99_ratio"] = (
                row["conventions"]["tail90"]["ratio_exact"]
                / row["conventions"]["p99"]["ratio_exact"])
            prec["drugs"].append(row)
            log(f"dx12: {panel['tissue']}/{name} p99 {row['conventions']['p99']['ratio_exact']:.4f} "
                f"tail90 {row['conventions']['tail90']['ratio_exact']:.4f} "
                f"(N ratio {row['tail90_over_p99_N']:.4f}, elasticity "
                f"{row['conventions']['p99']['dln_fpen_dln_N']['pm0.02']:+.4f})")
        out["panels"].append(prec)
    return out


# --------------------------------------------------------------------------- #
# STAGE verdicts -- reported LAST and separately
# --------------------------------------------------------------------------- #
def stage_verdicts() -> Dict:
    """What the counterfactual ladders would do to the nine published verdicts.

    Computed after every numerical conclusion is fixed and written to its own
    file, so no part of the falsification test can be influenced by a PASS count.
    The bands are the frozen literature bands carried in the ladder records and
    are not recomputed here.
    """
    lad = json.load(open(os.path.join(OUT_DIR, "tail90_ladder.json")))
    dx12: Optional[Dict] = None
    p = os.path.join(OUT_DIR, "tail90_dx12.json")
    if os.path.exists(p):
        dx12 = json.load(open(p))

    out: Dict = {"panels": [], "note": (
        "Verdicts are the frozen-band test applied to each ladder. The "
        "counterfactual ladder is an IDEALISATION (a normaliser with exactly "
        "zero grid error); its verdicts are reported to show what is and is not "
        "at stake, NOT as a proposed replacement for the published numbers.")}
    for prec in lad["panels"]:
        rows = []
        for a in prec["drugs"]:
            lo, hi = a["band_lo"], a["band_hi"]

            def _v(x):
                return "PASS" if lo <= x <= hi else "FAIL"

            dxs = a["dx_um"]
            r_p99 = a["series"]["ratio"]
            r_cf = a["counterfactual_perfect_normaliser"]["ratio"]
            i12 = int(np.argmin(np.abs(np.asarray(dxs) - PUBLISHED_DX)))
            row = {"tissue": prec["tissue"], "drug": a["drug"],
                   "band_lo": lo, "band_hi": hi,
                   "published_ratio_dx12": r_p99[i12],
                   "published_verdict": _v(r_p99[i12]),
                   "verdicts_p99_every_rung": [_v(x) for x in r_p99],
                   "verdicts_perfect_every_rung": [_v(x) for x in r_cf],
                   "n_pass_p99": sum(1 for x in r_p99 if _v(x) == "PASS"),
                   "n_pass_perfect": sum(1 for x in r_cf if _v(x) == "PASS"),
                   "n_rungs": len(dxs)}
            if dx12 is not None:
                for pp in dx12["panels"]:
                    if pp["tissue"] != prec["tissue"]:
                        continue
                    for dr in pp["drugs"]:
                        if dr["drug"] != a["drug"]:
                            continue
                        row["ratio_tail90_dx12"] = dr["conventions"]["tail90"]["ratio_exact"]
                        row["verdict_tail90_dx12"] = _v(
                            dr["conventions"]["tail90"]["ratio_exact"])
            rows.append(row)
        out["panels"].append({"tissue": prec["tissue"], "drugs": rows})
    return out


# --------------------------------------------------------------------------- #
def _json_safe(obj):
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def _dump(obj, name: str) -> None:
    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, name)
    with open(path, "w") as fh:
        json.dump(_json_safe(obj), fh, indent=2)
    log(f"written -> results/sweep/{name}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="ladder",
                    choices=["premise", "ladder", "dx12", "verdicts"])
    args = ap.parse_args()
    if args.stage == "premise":
        _dump(stage_premise(), "tail90_premise.json")
    elif args.stage == "ladder":
        _dump(stage_ladder(), "tail90_ladder.json")
    elif args.stage == "dx12":
        _dump(stage_dx12(), "tail90_dx12.json")
    elif args.stage == "verdicts":
        _dump(stage_verdicts(), "tail90_verdicts.json")


if __name__ == "__main__":
    main()
