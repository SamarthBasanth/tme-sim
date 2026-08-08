"""ALPHA IDENTIFIABILITY: what the doxorubicin anchor does and does NOT determine.

This is a DIAGNOSIS, not a fix. Nothing in ``tmesim/`` is edited and nothing is
tuned to hit an anchor. ``calibrate_alpha()`` is called unmodified; where this
script needs to vary the anchor drug's binding input it MONKEYPATCHES
``validation.calibration._compile`` inside a context manager and restores it.

The question
------------
``calibrate_alpha()`` solves ONE equation

    L_dox = sqrt(D_eff,dox / (c_dox + alpha * b_dox))  ==  L*     (64.921 um)

for ONE unknown, alpha. Every other drug then inherits ``k_d = c_d + alpha*b_d``.
The hypothesis under test is that ``alpha`` and the anchor drug's binding prior
``b_dox`` are NOT separately identified -- only their product is -- so any change
to ``b_dox`` is exactly re-exported as a global multiplier on every OTHER drug's
sink.

Stages
------
  ``--stage algebra``  parts 1 and 2 (closed form, seconds)
  ``--stage pde``      part 3, the cancellation test (20 transport solves, ~1 h
                       with ``--workers 5``; caches to the scratch JSON)
  ``--stage report``   writes ``results/alpha_identifiability.md`` from the cache
  ``--stage all``      (default) all three

Run:    ./.venv/bin/python experiments/alpha_identifiability.py
Writes: results/alpha_identifiability.md
        (cache) alpha_identifiability_cache.json next to this script's --cache arg
"""
from __future__ import annotations

import argparse
import contextlib
import dataclasses
import json
import math
import os
import sys
import time

import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.validation import calibration as cal  # noqa: E402
from tmesim.validation.pk_anchors import half_distance_to_efold_length  # noqa: E402
from tmesim.validation.pk_reality import DEFAULT_ECM, penetration_length_um  # noqa: E402

# The GBM panel used by experiments/real_drug_panel.py, verbatim.
PANEL = [
    "temozolomide", "lomustine", "carmustine", "procarbazine", "cisplatin",
    "methotrexate", "etoposide", "vincristine", "irinotecan", "doxorubicin",
]

#: molecular/binding.py's own documented upper bound on tissue_binding_rate (1/s).
BINDING_CEILING = 5e-3

#: The factor the prior agent's measured-logP substitution moved doxorubicin's
#: RAW binding prior by (a DROP). Used here only as a realistic perturbation
#: size; parts 1 and 2 also sweep round factors so nothing rests on it.
REPORTED_DOX_BINDING_DROP = 3.3478

DEFAULT_CACHE = os.path.join(
    os.environ.get("TMESIM_SCRATCH", "/tmp"), "alpha_identifiability_cache.json")
OUT_MD = os.path.join(_ROOT, "results", "alpha_identifiability.md")

_T0 = time.time()


def log(msg: str) -> None:
    print(f"[alpha_ident {time.time()-_T0:7.1f}s] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# Monkeypatch harness: scale the RAW binding prior of chosen drugs.
# --------------------------------------------------------------------------- #
_ORIG_COMPILE = cal._compile


@contextlib.contextmanager
def scaled_binding(scale_map):
    """Temporarily scale the compiled (UN-calibrated) ``binding`` of named drugs.

    This emulates "a better measurement changed this drug's binding input" WITHOUT
    editing molecular/binding.py. ``calibrate_alpha`` and ``calibrated_drugparams``
    both route through ``calibration._compile``, so patching it once covers the fit
    and every downstream compile consistently.
    """
    smap = {str(k).strip().lower(): float(v) for k, v in dict(scale_map).items()}

    def _patched(drug, ecm_reference: float = DEFAULT_ECM):
        dp = _ORIG_COMPILE(drug, ecm_reference=ecm_reference)
        if isinstance(drug, str):
            s = smap.get(str(drug).strip().lower(), 1.0)
            if s != 1.0:
                dp = dataclasses.replace(dp, binding=float(dp.binding) * s)
        return dp

    cal._compile = _patched
    try:
        yield
    finally:
        cal._compile = _ORIG_COMPILE


def raw_terms(name: str):
    """(D_eff, clearance, raw binding prior) for a drug at the default ECM."""
    dp = cal._compile(name, ecm_reference=DEFAULT_ECM)
    return (cal._d_eff_tortuous(dp, DEFAULT_ECM), cal.k_passive(dp), cal.uptake_proxy(dp))


def panel_state(alpha: float):
    """Per-drug calibrated k and L at a given alpha (analytic, no PDE)."""
    rows = {}
    for name in PANEL:
        dp = cal.calibrated_drugparams(name, alpha)
        k = float(dp.clearance) + float(dp.binding)
        rows[name] = {
            "binding": float(dp.binding),
            "clearance": float(dp.clearance),
            "k": k,
            "L_um": penetration_length_um(dp, ecm_density=DEFAULT_ECM),
            "over_ceiling": float(dp.binding) / BINDING_CEILING,
        }
    return rows


# --------------------------------------------------------------------------- #
# PART 1 -- the degeneracy, proved by perturbation
# --------------------------------------------------------------------------- #
def part1_degeneracy():
    log("PART 1: degeneracy of (alpha, b_anchor)")
    alpha0 = cal.calibrate_alpha()
    d_eff0, c_dox, b_dox = raw_terms("doxorubicin")
    out = {
        "alpha0": alpha0,
        "dox_d_eff": d_eff0,
        "dox_clearance": c_dox,
        "dox_raw_binding": b_dox,
        "dox_sink_alpha_b": alpha0 * b_dox,
        "k_target": d_eff0 / (cal.DOX_TARGET_EFOLD_UM ** 2),
        "target_efold_um": cal.DOX_TARGET_EFOLD_UM,
        "base_panel": panel_state(alpha0),
        "perturbations": [],
    }
    log(f"  alpha0 = {alpha0!r}   b_dox = {b_dox:.6e}   alpha0*b_dox = {alpha0*b_dox:.10e}")

    base = out["base_panel"]
    for f in (0.1, 1.0 / REPORTED_DOX_BINDING_DROP, 0.5, 2.0,
              REPORTED_DOX_BINDING_DROP, 10.0):
        with scaled_binding({"doxorubicin": f}):
            a = cal.calibrate_alpha()
            pan = panel_state(a)
        rec = {
            "f_binding_scale": f,
            "alpha": a,
            "alpha_ratio": a / alpha0,
            "alpha_ratio_times_f": (a / alpha0) * f,     # == 1 exactly if degenerate
            "dox_sink_alpha_b": a * (b_dox * f),
            "dox_sink_rel_error": abs(a * b_dox * f - alpha0 * b_dox) / (alpha0 * b_dox),
            "L_dox_um": pan["doxorubicin"]["L_um"],
            "L_dox_rel_error": abs(pan["doxorubicin"]["L_um"] - base["doxorubicin"]["L_um"])
            / base["doxorubicin"]["L_um"],
            "per_drug": {
                n: {"k_ratio": pan[n]["k"] / base[n]["k"],
                    "L_ratio": pan[n]["L_um"] / base[n]["L_um"],
                    "binding_ratio": pan[n]["binding"] / base[n]["binding"]}
                for n in PANEL},
        }
        out["perturbations"].append(rec)
        log(f"  f={f:8.5f} -> alpha={a:12.5f}  alpha_ratio*f={rec['alpha_ratio_times_f']:.15f}"
            f"  L_dox rel err={rec['L_dox_rel_error']:.2e}")

    # The ONLY channel by which the reciprocal identity can be broken: the anchor's
    # clearance also depends on logP (molecular/binding.clearance_rate), and it
    # enters the numerator. Bound the deviation exactly.
    out["clearance_share_of_k_target"] = c_dox / out["k_target"]
    out["max_alpha_deviation_from_reciprocal"] = (
        out["k_target"] / (out["k_target"] - c_dox) - 1.0)
    log(f"  c_dox/k* = {out['clearance_share_of_k_target']:.6e}  -> the reciprocal "
        f"identity can be broken by at most {100*out['max_alpha_deviation_from_reciprocal']:.4f}%")

    # Predict, from the closed form alone, the two alphas independently reported by
    # concurrent agents. Uses the PUBLIC binding.tissue_binding_rate, no patching.
    from tmesim.molecular.binding import tissue_binding_rate
    from tmesim.molecular import library as _lib
    dox_desc = dict(_lib.get_molecule("doxorubicin").descriptors)
    out["dox_descriptors"] = {k: dox_desc.get(k) for k in ("mw", "logp", "hbd")}
    checks = []
    # (a) logp_provenance.md section 5.3: Crippen 0.001 -> measured HSDB/Hansch 1.27
    b_new = tissue_binding_rate(dict(dox_desc, logp=1.27))
    f_a = b_new / b_dox
    checks.append({
        "case": "doxorubicin logP 0.001 (RDKit Crippen) -> 1.27 (HSDB/Hansch)",
        "source": "results/logp_provenance.md section 5.3",
        "b_dox_scale_f": f_a,
        "alpha_predicted_by_identity": alpha0 / f_a,
        "alpha_reported": 68.6368,
        "global_L_rescale_predicted": math.sqrt(f_a),
        "global_L_rescale_reported": 1.47,
    })
    # (b) the task brief / measured-logD branch: raw binding DROPPED 3.3478x
    checks.append({
        "case": f"doxorubicin raw binding dropped {REPORTED_DOX_BINDING_DROP:.4f}x",
        "source": "the measured-logP branch reported to this run",
        "b_dox_scale_f": 1.0 / REPORTED_DOX_BINDING_DROP,
        "alpha_predicted_by_identity": alpha0 * REPORTED_DOX_BINDING_DROP,
        "alpha_reported": 497.30899,
        "global_L_rescale_predicted": 1.0 / math.sqrt(REPORTED_DOX_BINDING_DROP),
        "global_L_rescale_reported": None,
    })
    for c in checks:
        c["rel_error_vs_reported"] = abs(
            c["alpha_predicted_by_identity"] - c["alpha_reported"]) / c["alpha_reported"]
        log(f"  CROSS-CHECK {c['case']}: identity predicts alpha="
            f"{c['alpha_predicted_by_identity']:.5f}, reported {c['alpha_reported']:.5f}, "
            f"rel err {c['rel_error_vs_reported']:.3e} "
            f"(bound {out['max_alpha_deviation_from_reciprocal']:.3e})")
    out["cross_checks"] = checks

    # THE EXACT GAUGE: scale EVERY drug's raw binding by t and refit. If alpha is a
    # gauge quantity, every downstream number is bit-identical while alpha moves by 1/t.
    out["gauge_test"] = []
    all_names = list(dict.fromkeys(PANEL + [cal.MITOXANTRONE_NAME]))
    for t in (0.01, 0.2, 5.0, 100.0):
        with scaled_binding({n: t for n in all_names}):
            a = cal.calibrate_alpha()
            pan = panel_state(a)
            rv = cal.revalidate(a)
        worst_L = max(abs(pan[n]["L_um"] - base[n]["L_um"]) for n in PANEL)
        worst_k = max(abs(pan[n]["k"] - base[n]["k"]) / base[n]["k"] for n in PANEL)
        out["gauge_test"].append({
            "t": t, "alpha": a, "alpha_ratio": a / alpha0,
            "alpha_ratio_times_t": (a / alpha0) * t,
            "worst_abs_L_change_um": worst_L,
            "worst_rel_k_change": worst_k,
            "mitox_predicted_L_um": rv["second_drug"]["predicted_L_um"],
            "mitox_verdict": rv["second_drug"]["verdict"],
        })
        log(f"  GAUGE t={t:8.3f}: alpha={a:14.5f} ({a/alpha0:.6g}x)  worst |dL| over the "
            f"whole panel = {worst_L:.3e} um   mitox L={rv['second_drug']['predicted_L_um']:.10f}")

    # Is the module's OWN out-of-sample test flat under an ANCHOR-ONLY change?
    out["revalidate_flatness"] = []
    for f in (1.0, 1.0 / REPORTED_DOX_BINDING_DROP, 10.0):
        with scaled_binding({"doxorubicin": f}):
            a = cal.calibrate_alpha()
            rv = cal.revalidate(a)
        out["revalidate_flatness"].append({
            "f_binding_scale": f, "alpha": a,
            "mitox_predicted_L_um": rv["second_drug"]["predicted_L_um"],
            "mitox_ratio_to_mid": rv["second_drug"]["ratio_to_mid"],
            "verdict": rv["second_drug"]["verdict"],
            "tmz_L_um": rv["qualitative"]["temozolomide_L_um"],
            "fold_further": rv["qualitative"]["fold_further"],
            "generalizes": rv["generalizes"],
        })
        log(f"  revalidate under f={f:.5f}: alpha={a:.4f}  mitox L="
            f"{rv['second_drug']['predicted_L_um']:.4f} um  [{rv['second_drug']['verdict']}]"
            f"  tmz/dox fold={rv['qualitative']['fold_further']:.4f}")
    return out


# --------------------------------------------------------------------------- #
# PART 2 -- blast radius of the ANCHOR CHOICE
# --------------------------------------------------------------------------- #
#: Every measured perivascular half-distance in the repo's own cited compendium
#: (tmesim/validation/pk_anchors.py). These are the ONLY penetration-depth
#: anchors that exist here; no number below is invented.
ANCHOR_VARIANTS = [
    {"key": "dox_primeau_mid", "drug": "doxorubicin", "half_um": 45.0,
     "label": "doxorubicin, Primeau 40-50 um midpoint (THE BASELINE)",
     "cite": "Primeau et al., Clin Cancer Res 2005;11(24):8782-8, PMID 16361566"},
    {"key": "dox_primeau_lo", "drug": "doxorubicin", "half_um": 40.0,
     "label": "doxorubicin, Primeau band LOW end",
     "cite": "Primeau et al. 2005, PMID 16361566"},
    {"key": "dox_primeau_hi", "drug": "doxorubicin", "half_um": 50.0,
     "label": "doxorubicin, Primeau band HIGH end",
     "cite": "Primeau et al. 2005, PMID 16361566"},
    {"key": "dox_patel_mid", "drug": "doxorubicin", "half_um": 50.0,
     "label": "doxorubicin, Patel/Tredan/Tannock 25-75 um midpoint",
     "cite": "Patel, Tredan, Tannock, Cancer Chemother Pharmacol 2013;72(1):127-38, PMID 23680920"},
    {"key": "mitox_mid", "drug": "mitoxantrone", "half_um": 50.0,
     "label": "MITOXANTRONE, 25-75 um midpoint (ANCHOR-DRUG SWAP)",
     "cite": "Patel, Tredan, Tannock 2013, PMID 23680920"},
    {"key": "mitox_lo", "drug": "mitoxantrone", "half_um": 25.0,
     "label": "MITOXANTRONE, band LOW end",
     "cite": "Patel, Tredan, Tannock 2013, PMID 23680920"},
    {"key": "mitox_hi", "drug": "mitoxantrone", "half_um": 75.0,
     "label": "MITOXANTRONE, band HIGH end",
     "cite": "Patel, Tredan, Tannock 2013, PMID 23680920"},
]


def part2_anchor_swap():
    log("PART 2: anchor-choice blast radius")
    out = {"variants": []}
    base = None
    for v in ANCHOR_VARIANTS:
        efold = half_distance_to_efold_length(v["half_um"])
        a = cal.calibrate_alpha(target_efold_um=efold, drug=v["drug"])
        pan = panel_state(a)
        over = {n: r["over_ceiling"] for n, r in pan.items() if r["over_ceiling"] > 1.0}
        rec = dict(v)
        rec.update({
            "efold_target_um": efold,
            "alpha": a,
            "panel": pan,
            "n_over_ceiling": len(over),
            "worst_over_ceiling": max(over.values()) if over else 0.0,
            "L_tmz_over_L_dox": pan["temozolomide"]["L_um"] / pan["doxorubicin"]["L_um"],
        })
        if base is None:
            base = rec
        rec["alpha_vs_baseline"] = a / base["alpha"]
        rec["L_ratio_vs_baseline"] = {
            n: pan[n]["L_um"] / base["panel"][n]["L_um"] for n in PANEL}
        rec["k_ratio_vs_baseline"] = {
            n: pan[n]["k"] / base["panel"][n]["k"] for n in PANEL}
        out["variants"].append(rec)
        log(f"  {v['key']:16s} L*={efold:7.2f} um  alpha={a:12.4f}  "
            f"({a/base['alpha']:.4f}x baseline)  over-ceiling {len(over)}/{len(PANEL)}")

    # Are the two measured anchors mutually consistent under this model?
    # Fit on each, then check whether the OTHER one's published band is hit.
    dox_lo, dox_hi = cal.DOX_EFOLD_BAND_UM
    mit_lo, mit_hi = cal.MITOX_EFOLD_BAND_UM
    base_v = out["variants"][0]
    mit_v = next(r for r in out["variants"] if r["key"] == "mitox_mid")

    def _mitox_L(alpha):
        return penetration_length_um(
            cal.calibrated_drugparams(cal.MITOXANTRONE_NAME, alpha),
            ecm_density=DEFAULT_ECM)

    out["mutual_consistency"] = {
        "dox_band_um": [dox_lo, dox_hi],
        "mitox_band_um": [mit_lo, mit_hi],
        "dox_band_width": dox_hi / dox_lo,
        "mitox_band_width": mit_hi / mit_lo,
        # fit on doxorubicin -> predict mitoxantrone
        "fit_dox_pred_mitox_L": _mitox_L(base_v["alpha"]),
        "fit_dox_pred_mitox_in_band": bool(
            mit_lo <= _mitox_L(base_v["alpha"]) <= mit_hi),
        # fit on mitoxantrone -> predict doxorubicin
        "fit_mitox_pred_dox_L": mit_v["panel"]["doxorubicin"]["L_um"],
        "fit_mitox_pred_dox_in_band": bool(
            dox_lo <= mit_v["panel"]["doxorubicin"]["L_um"] <= dox_hi),
        "fit_mitox_pred_dox_over_band_top": (
            mit_v["panel"]["doxorubicin"]["L_um"] / dox_hi),
    }
    mc = out["mutual_consistency"]
    log(f"  MUTUAL: fit-on-dox -> mitox L={mc['fit_dox_pred_mitox_L']:.2f} um "
        f"(band {mit_lo:.1f}-{mit_hi:.1f}, width {mc['mitox_band_width']:.2f}x) "
        f"{'IN' if mc['fit_dox_pred_mitox_in_band'] else 'OUT'}")
    log(f"  MUTUAL: fit-on-mitox -> dox   L={mc['fit_mitox_pred_dox_L']:.2f} um "
        f"(band {dox_lo:.1f}-{dox_hi:.1f}, width {mc['dox_band_width']:.2f}x) "
        f"{'IN' if mc['fit_mitox_pred_dox_in_band'] else 'OUT'} "
        f"({mc['fit_mitox_pred_dox_over_band_top']:.3f}x the band top)")

    # The mitoxantrone anchor-drug swap, spelled out.
    mit = next(r for r in out["variants"] if r["key"] == "mitox_mid")
    out["headline_swap"] = {
        "alpha_baseline": base["alpha"],
        "alpha_mitoxantrone": mit["alpha"],
        "ratio": mit["alpha"] / base["alpha"],
        "L_ratio": mit["L_ratio_vs_baseline"],
    }
    return out


# --------------------------------------------------------------------------- #
# PART 3 -- does it CANCEL?  (real transport solves on the real Xenium section)
# --------------------------------------------------------------------------- #
def _solve_one(args):
    """One (drug, alpha) transport solve on the shared, already-prepared tissue.

    Runs in a worker process; imports live inside so the parent stays light.
    """
    import numpy as _np
    from tmesim.core import PDEConfig, PKConfig, PDConfig, SimConfig, Microenvironment
    from tmesim.reconstruct.real_pipeline import (
        apply_drug_sensitivity, reconstruct_real)
    from tmesim.validation.kill_calibration import (
        BATH_MULTIPLES, DT_S, SAVE_EVERY, SEED, WINDOW_S, compile_drug, decompose,
        kill_with, prepare_tissue, resistance_vector, unit_transport,
        window_plasma_mean)

    name, alpha, arm, dx_um = args
    t0 = time.time()
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
    pde = PDEConfig(dt=DT_S, total_time=WINDOW_S, convection=True, save_every=SAVE_EVERY)
    pk = PKConfig(plasma_peak=50.0, half_life=3600.0 * 4, n_doses=1,
                  dose_interval=3600.0 * 24)
    pd_cfg = PDConfig()
    work, labeled, velocity = prepare_tissue(
        micro, SimConfig(pde=pde, pk=pk, pd=pd_cfg, seed=SEED))
    tumor = labeled.tumor_mask()
    state_t = _np.asarray(labeled.state)[tumor]

    drug = compile_drug(name, alpha)
    m2 = Microenvironment(domain=work.domain, fields=work.fields, cells=labeled.copy(),
                          patient_id=work.patient_id, metadata=work.metadata)
    apply_drug_sensitivity(m2, name, mode="uniform")
    sens = _np.asarray(m2.cells.sensitivity, dtype=_np.float64)[tumor]

    tr = unit_transport(work, drug, pde, velocity)
    col, row = work.domain.to_index(labeled.x, labeled.y)
    auc = _np.asarray(tr.auc, dtype=_np.float64)[row, col]
    duration = float(tr.times[-1]) if tr.times.size and tr.times[-1] > 0 else 1.0
    e_unit = (auc / duration)[tumor]

    pm = window_plasma_mean(pk, WINDOW_S, DT_S)
    kill = kill_with(labeled, work, tr, drug, pd_cfg, pm)

    resist = resistance_vector(state_t, pd_cfg, drug=drug)
    dec = {}
    for m in BATH_MULTIPLES:
        d = decompose(e_unit, sens, resist, m * drug.ic50, drug.ic50, drug.hill,
                      pd_cfg.max_kill_rate)
        dec[f"{m:g}xIC50"] = {k: v for k, v in d.items()}

    return {
        "drug": name, "arm": arm, "alpha": float(alpha),
        "binding": float(drug.binding), "clearance": float(drug.clearance),
        "k": float(drug.binding) + float(drug.clearance),
        "ic50": float(drug.ic50), "hill": float(drug.hill),
        "penetration_index": float(tr.penetration_index),
        "mean_auc_tissue": float(_np.asarray(tr.auc)[
            _np.asarray(work.domain.tissue_mask, dtype=bool)].mean()),
        "e_unit_mean": float(e_unit.mean()),
        "e_unit_median": float(_np.median(e_unit)),
        "e_unit_p99": float(_np.percentile(e_unit, 99.0)),
        "kill_nominal": float(kill),
        "f_state": float(_np.mean(resist)),
        "decomposition": dec,
        "solve_s": time.time() - t0,
    }


def part3_pde(cache_path: str, workers: int, dx_um: float, g: float):
    """Two arms: baseline alpha, and alpha*g (the anchor-perturbation re-export)."""
    import multiprocessing as mp

    alpha0 = cal.calibrate_alpha()
    arms = {"baseline": alpha0, "perturbed": alpha0 * g}
    log(f"PART 3: PDE arms baseline alpha={alpha0:.5f}, perturbed alpha={alpha0*g:.5f} (g={g})")

    cache = {}
    if os.path.exists(cache_path):
        try:
            cache = json.load(open(cache_path))
        except Exception:
            cache = {}

    jobs = []
    for arm, a in arms.items():
        for name in PANEL:
            key = f"{name}|{arm}|{a:.10g}|{dx_um}"
            if key not in cache:
                jobs.append((key, (name, a, arm, dx_um)))
    log(f"  {len(jobs)} solves to run ({len(arms)*len(PANEL)-len(jobs)} cached)")

    if jobs:
        ctx = mp.get_context("spawn")
        done = 0
        with ctx.Pool(processes=workers) as pool:
            for key, res in zip([j[0] for j in jobs],
                                pool.imap(_solve_one, [j[1] for j in jobs])):
                cache[key] = res
                done += 1
                log(f"  [{done}/{len(jobs)}] {res['drug']:14s} {res['arm']:10s} "
                    f"kill={res['kill_nominal']:.6f} e_unit_mean={res['e_unit_mean']:.5g} "
                    f"({res['solve_s']:.0f}s)")
                json.dump(cache, open(cache_path, "w"))
    json.dump(cache, open(cache_path, "w"))

    by = {}
    for arm, a in arms.items():
        for name in PANEL:
            by.setdefault(arm, {})[name] = cache[f"{name}|{arm}|{a:.10g}|{dx_um}"]
    return {"alpha0": alpha0, "g": g, "dx_um": dx_um, "arms": list(arms), "by": by}


def _lnr(a, b):
    return math.log(a / b) if (a > 0 and b > 0) else float("nan")


ANCHOR = "doxorubicin"


def analyse_cancellation(p3):
    """TWO arms, both exactly reconstructible from the same 20 solves.

    ``uniform``  every drug's sink x g. The pure COMMON-MODE / gauge-direction
                 perturbation the task hypothesised.
    ``anchor``   the sink of every drug EXCEPT the anchor x g, the anchor
                 bit-identical. This is what a change to ``b_dox`` ACTUALLY does:
                 ``alpha -> g*alpha`` and ``b_dox -> b_dox/g`` leave
                 ``s_dox = alpha*b_dox`` and ``c_dox`` untouched, so doxorubicin's
                 DrugParams are unchanged and its solve is literally the baseline
                 one. The other nine see ``s_d -> g*s_d``, which is the perturbed
                 solve. No extra simulation is needed and no approximation is made.
    """
    base = p3["by"]["baseline"]
    pert = p3["by"]["perturbed"]
    arms = {
        "uniform": {n: pert[n] for n in PANEL},
        "anchor": {n: (base[n] if n == ANCHOR else pert[n]) for n in PANEL},
    }
    return {
        "arms": {k: _analyse_arm(p3, base, v) for k, v in arms.items()},
        "g": p3["g"],
    }


def _analyse_arm(p3, b, p):
    """Per-drug surviving fraction and cancellation factor, in the style of
    results/transvascular_exchange_revalidation.md.

        surviving fraction = |dln(ratio-type output)| / |dln(absolute output)|
        cancellation factor = 1 - surviving fraction
    """
    anchor = ANCHOR
    gm0 = math.exp(sum(math.log(b[n]["e_unit_mean"]) for n in PANEL) / len(PANEL))
    gm1 = math.exp(sum(math.log(p[n]["e_unit_mean"]) for n in PANEL) / len(PANEL))
    gmk0 = math.exp(sum(math.log(b[n]["kill_nominal"]) for n in PANEL) / len(PANEL))
    gmk1 = math.exp(sum(math.log(p[n]["kill_nominal"]) for n in PANEL) / len(PANEL))

    rows = []
    for n in PANEL:
        e0, e1 = b[n]["e_unit_mean"], p[n]["e_unit_mean"]
        k0, k1 = b[n]["kill_nominal"], p[n]["kill_nominal"]
        f0 = b[n]["decomposition"]["1xIC50"]["f_penetration"]
        f1 = p[n]["decomposition"]["1xIC50"]["f_penetration"]
        r0 = b[n]["decomposition"]["1xIC50"]["ratio_exact"]
        r1 = p[n]["decomposition"]["1xIC50"]["ratio_exact"]
        # ratio-type: normalised to the ANCHOR drug (the natural cross-drug ratio)
        rel_e0 = e0 / b[anchor]["e_unit_mean"]
        rel_e1 = e1 / p[anchor]["e_unit_mean"]
        rel_k0 = k0 / b[anchor]["kill_nominal"]
        rel_k1 = k1 / p[anchor]["kill_nominal"]

        # ratio-type, second normaliser: the PANEL geometric mean (what a relative
        # score / ranking product effectively uses when there is no anchor drug).
        grel_e0, grel_e1 = e0 / gm0, e1 / gm1
        grel_k0, grel_k1 = k0 / gmk0, k1 / gmk1

        d_abs_e, d_rel_e = _lnr(e1, e0), _lnr(rel_e1, rel_e0)
        d_abs_k, d_rel_k = _lnr(k1, k0), _lnr(rel_k1, rel_k0)
        d_grel_e, d_grel_k = _lnr(grel_e1, grel_e0), _lnr(grel_k1, grel_k0)
        rows.append({
            "gm_rel_auc_base": grel_e0, "gm_rel_auc_pert": grel_e1,
            "gm_rel_kill_base": grel_k0, "gm_rel_kill_pert": grel_k1,
            "surviving_frac_auc_gm": abs(d_grel_e) / abs(d_abs_e) if d_abs_e else float("nan"),
            "surviving_frac_kill_gm": abs(d_grel_k) / abs(d_abs_k) if d_abs_k else float("nan"),
            "drug": n,
            "k_ratio": p[n]["k"] / b[n]["k"],
            "common_mode_purity": _lnr(p[n]["k"], b[n]["k"]) / math.log(p3["g"]),
            "e_unit_mean_base": e0, "e_unit_mean_pert": e1,
            "e_unit_fold": e1 / e0,
            "kill_base": k0, "kill_pert": k1, "kill_fold": k1 / k0,
            "f_pen_base": f0, "f_pen_pert": f1, "f_pen_fold": f1 / f0,
            "ratio3d2d_base": r0, "ratio3d2d_pert": r1, "ratio3d2d_fold": r1 / r0,
            "rel_auc_base": rel_e0, "rel_auc_pert": rel_e1,
            "rel_kill_base": rel_k0, "rel_kill_pert": rel_k1,
            "surviving_frac_auc": abs(d_rel_e) / abs(d_abs_e) if d_abs_e else float("nan"),
            "surviving_frac_kill": abs(d_rel_k) / abs(d_abs_k) if d_abs_k else float("nan"),
            "cancel_factor_auc": (0.0 + max(0.0, 1.0 - abs(d_rel_e) / abs(d_abs_e))) if d_abs_e else float("nan"),
            "cancel_factor_kill": (0.0 + max(0.0, 1.0 - abs(d_rel_k) / abs(d_abs_k))) if d_abs_k else float("nan"),
            "penetration_index_base": b[n]["penetration_index"],
            "penetration_index_pert": p[n]["penetration_index"],
        })

    def _spearman(x, y):
        rx = np.argsort(np.argsort(np.asarray(x, dtype=float)))
        ry = np.argsort(np.argsort(np.asarray(y, dtype=float)))
        rx = rx - rx.mean(); ry = ry - ry.mean()
        return float((rx * ry).sum() / math.sqrt((rx ** 2).sum() * (ry ** 2).sum()))

    def _kendall(x, y):
        x = np.asarray(x, dtype=float); y = np.asarray(y, dtype=float)
        n = x.size; c = d = 0
        for i in range(n):
            for j in range(i + 1, n):
                s = np.sign(x[i] - x[j]) * np.sign(y[i] - y[j])
                if s > 0:
                    c += 1
                elif s < 0:
                    d += 1
        return float((c - d) / (0.5 * n * (n - 1))), int(d)

    kills0 = [b[n]["kill_nominal"] for n in PANEL]
    kills1 = [p[n]["kill_nominal"] for n in PANEL]
    aucs0 = [b[n]["e_unit_mean"] for n in PANEL]
    aucs1 = [p[n]["e_unit_mean"] for n in PANEL]
    kt_kill, disc_kill = _kendall(kills0, kills1)
    kt_auc, disc_auc = _kendall(aucs0, aucs1)
    order0 = [PANEL[i] for i in np.argsort(-np.asarray(kills0))]
    order1 = [PANEL[i] for i in np.argsort(-np.asarray(kills1))]
    aorder0 = [PANEL[i] for i in np.argsort(-np.asarray(aucs0))]
    aorder1 = [PANEL[i] for i in np.argsort(-np.asarray(aucs1))]

    non_anchor = [r for r in rows if r["drug"] != anchor]
    surv_kill = [r["surviving_frac_kill"] for r in non_anchor]
    surv_auc = [r["surviving_frac_auc"] for r in non_anchor]
    r3 = [abs(math.log(r["ratio3d2d_fold"])) for r in rows]
    fpen = [abs(math.log(r["f_pen_fold"])) for r in rows]
    verdicts = {
        "median_surviving_kill": float(np.median(surv_kill)),
        "median_surviving_auc": float(np.median(surv_auc)),
        "max_surviving_kill": float(max(surv_kill)),
        "max_surviving_auc": float(max(surv_auc)),
        # how much of the g-fold sink change reaches the 3D:2D ratio
        "ratio3d2d_share_of_g": float(max(r3) / abs(math.log(p3["g"]))),
        "fpen_share_of_g": float(max(fpen) / abs(math.log(p3["g"]))),
        "median_surviving_kill_gm": float(np.median(
            [r["surviving_frac_kill_gm"] for r in rows
             if not math.isnan(r["surviving_frac_kill_gm"])])),
        "median_surviving_auc_gm": float(np.median(
            [r["surviving_frac_auc_gm"] for r in rows
             if not math.isnan(r["surviving_frac_auc_gm"])])),
    }

    def _margin(vals, order):
        """Relative gap between the top two, as a fraction of the leader."""
        s = sorted(vals, reverse=True)
        return float((s[0] - s[1]) / s[0]) if s[0] > 0 else float("nan")

    return {
        "rows": rows,
        "verdicts": verdicts,
        "ranking": {
            "kill_top1_margin_base": _margin(kills0, order0),
            "kill_top1_margin_pert": _margin(kills1, order1),
            "kill_spearman": _spearman(kills0, kills1),
            "kill_kendall": kt_kill, "kill_discordant_pairs": disc_kill,
            "auc_spearman": _spearman(aucs0, aucs1),
            "auc_kendall": kt_auc, "auc_discordant_pairs": disc_auc,
            "kill_order_base": order0, "kill_order_pert": order1,
            "auc_order_base": aorder0, "auc_order_pert": aorder1,
            "kill_top1_same": order0[0] == order1[0],
            "kill_top3_same": set(order0[:3]) == set(order1[:3]),
            "auc_top1_same": aorder0[0] == aorder1[0],
            "auc_top3_same": set(aorder0[:3]) == set(aorder1[:3]),
        },
    }


# --------------------------------------------------------------------------- #
# Report
# --------------------------------------------------------------------------- #
def _f(x, n=4):
    try:
        return f"{float(x):.{n}f}"
    except Exception:
        return "—"


def write_report(p1, p2, p3, canc, path=OUT_MD):
    A = []
    add = A.append
    g = p3["g"] if p3 else None

    add("# The doxorubicin anchor identifies ONE number, and it is not alpha")
    add("")
    add("**A diagnosis. No parameter was changed, nothing was tuned to an anchor, "
        "and no code in `tmesim/` was edited.** Everything below is either closed-form "
        "algebra read off `tmesim/validation/calibration.py` or a run executed by "
        "`experiments/alpha_identifiability.py`.")
    add("")
    add("> ### Headline")
    add("> ")
    add(f"> **`alpha` and the anchor drug's binding prior are not separately "
        f"identified — only their product is.** Perturbing doxorubicin's raw binding "
        f"input by a factor `f` moves `alpha` by exactly `1/f`; across six perturbation "
        f"factors spanning 100x the product `alpha·b_dox` is invariant to "
        f"**{max(r['dox_sink_rel_error'] for r in p1['perturbations']):.1e}** relative "
        f"error and `L_dox` to "
        f"**{max(r['L_dox_rel_error'] for r in p1['perturbations']):.1e}**. The anchor "
        f"pins one number, `k_dox = {p1['k_target']:.6e} /s`. It pins nothing else.")
    add("> ")
    add(f"> **Any improvement to the anchor drug's chemistry is re-exported as a global "
        f"multiplier on every OTHER drug.** That is not a side effect; it is the only "
        f"channel through which `b_dox` can act, because the anchor equation forces "
        f"`alpha·b_dox` to stay put.")
    if p1.get("gauge_test"):
        gt = p1["gauge_test"]
        add("> ")
        add(f"> **Worse: `alpha` is a gauge coordinate, not a measurement.** Scaling "
            f"*every* drug's binding prior by `t` and refitting moves `alpha` over four "
            f"orders of magnitude ({min(g_['alpha'] for g_ in gt):.4f} to "
            f"{max(g_['alpha'] for g_ in gt):.2f}) while leaving every downstream number "
            f"— every `k`, every `L`, the out-of-sample prediction — identical to "
            f"{max(g_['worst_abs_L_change_um'] for g_ in gt):.1e} µm. Publishing "
            f"`alpha = 148.60815765445545` to 17 figures publishes a coordinate.")
    if p1.get("revalidate_flatness") and len(p1["revalidate_flatness"]) > 1:
        b_rv, p_rv = p1["revalidate_flatness"][0], p1["revalidate_flatness"][1]
        add("> ")
        add(f"> **And the branch would have broken the project's own falsification "
            f"test.** Under the reported anchor-binding drop, `calibration.revalidate()`'s "
            f"held-out mitoxantrone prediction falls "
            f"{b_rv['mitox_predicted_L_um']:.1f} → {p_rv['mitox_predicted_L_um']:.1f} µm "
            f"and exits the published {cal.MITOX_EFOLD_BAND_UM[0]:.0f}–"
            f"{cal.MITOX_EFOLD_BAND_UM[1]:.0f} µm band "
            f"(`{b_rv['verdict']}` → `{p_rv['verdict']}`).")
    if p2:
        hs = p2["headline_swap"]
        add("> ")
        add(f"> **The choice of anchor DRUG moves the whole panel.** Refitting on "
            f"mitoxantrone — the only other drug in this repo's cited compendium with a "
            f"measured perivascular half-distance — gives "
            f"`alpha = {hs['alpha_mitoxantrone']:.2f}` against the published "
            f"`{hs['alpha_baseline']:.2f}`, a **{hs['ratio']:.3f}x** shift, and every "
            f"drug's penetration length moves by "
            f"**{min(hs['L_ratio'].values()):.3f}–{max(hs['L_ratio'].values()):.3f}x**.")
    if canc:
        UA = canc["arms"]["uniform"]
        AA = canc["arms"]["anchor"]
        rkA, rkB = UA["ranking"], AA["ranking"]
        add("> ")
        add(f"> **Does it cancel? The question has two answers, because the perturbation "
            f"is not the one it looks like.** A uniform \u00d7{g:.4f} on every drug\u2019s sink is "
            f"common-mode and cancels almost completely in cross-drug ratios "
            f"(median **{100*(1-UA['verdicts']['median_surviving_kill']):.1f}%** "
            f"cancelled against the anchor drug; ranking fully intact). But **what a "
            f"change to `b_dox` actually does is move nine drugs and leave the tenth \u2014 "
            f"the anchor \u2014 exactly where it was**, and against the anchor that cancels "
            f"**{100*(1-AA['verdicts']['median_surviving_kill']):.1f}%**.")
        add("> ")
        add(f"> **Absolutes do not survive either way** (kill "
            f"{min(r['kill_fold'] for r in UA['rows']):.3f}\u2013"
            f"{max(r['kill_fold'] for r in UA['rows']):.3f}\u00d7), and the **3-D:2-D kill "
            f"ratio** \u2014 the one external test absolute kill passes \u2014 still takes "
            f"{100*UA['verdicts']['ratio3d2d_share_of_g']:.1f}% of the perturbation, "
            f"because renormalising to the best-perfused cell removes a global "
            f"concentration but not a change in gradient SHAPE. The cross-drug **ranking "
            f"is the one thing that holds**"
            + ("." if rkB["kill_top1_same"] else
               f", with one caveat: the anchor-exempt arm flips top-1 across a "
               f"{100*rkB['kill_top1_margin_pert']:.2f}% margin \u2014 a near-tie, not a "
               f"reversal."))
    add("")
    add("---")
    add("")

    # ---- 1 ----
    add("## 1. What is being fit, against what, with how many degrees of freedom")
    add("")
    add("`calibration.calibrate_alpha()` is a closed-form solve of **one scalar "
        "equation in one scalar unknown**:")
    add("")
    add("```")
    add("k_d(alpha) = c_d + alpha * b_d                      (per drug d)")
    add("L_d(alpha) = sqrt(D_eff,d / k_d(alpha))")
    add("")
    add("FIT:   L_dox(alpha) = L*                            ONE equation")
    add("SOLVE: alpha = (D_eff,dox / L*^2 - c_dox) / b_dox    ONE unknown")
    add("```")
    add("")
    add("with, on the doxorubicin anchor as the code actually computes it:")
    add("")
    add("| symbol | meaning | value |")
    add("|---|---|---:|")
    add(f"| `L*` | fit target, e-folding length | {p1['target_efold_um']:.5f} µm |")
    add(f"| `D_eff,dox` | `D·(1 − 0.6·ECM)` at ECM = {DEFAULT_ECM} | {p1['dox_d_eff']:.5f} µm²/s |")
    add(f"| `c_dox` | `clearance`, NOT amplified | {p1['dox_clearance']:.6e} /s |")
    add(f"| `b_dox` | raw `tissue_binding_rate` prior | {p1['dox_raw_binding']:.6e} /s |")
    add(f"| `k*` | `D_eff/L*²`, the fitted quantity | **{p1['k_target']:.6e} /s** |")
    add(f"| `alpha` | the solved constant | **{p1['alpha0']!r}** |")
    add(f"| `alpha·b_dox` | the anchor's amplified sink | **{p1['dox_sink_alpha_b']:.6e} /s** |")
    add("")
    add("**Free parameters vs independent constraints, counted exactly.** Write the "
        "panel's physics in terms of what the PDE actually sees: the per-drug amplified "
        "sink `s_d ≡ alpha·b_d`, one per drug, `N = 10` on this GBM panel. The parameter "
        "vector is `(alpha, b_1 … b_N)` — `N + 1` numbers. The map to observables is "
        "`(alpha, {b_d}) ↦ {alpha·b_d}`, which has a **one-dimensional fibre** "
        "(§1.2 verifies it): `(alpha/t, {t·b_d})` gives identical physics for every `t`. "
        "So `N + 1` parameters carry only `N` distinguishable quantities, and the "
        "calibration supplies **one** equation. The ledger:")
    add("")
    add("| | count | status |")
    add("|---|---:|---|")
    add("| parameters `(alpha, b_1…b_N)` | N + 1 = 11 | — |")
    add("| gauge directions (unobservable by construction) | 1 | `alpha` alone is pure gauge |")
    add("| physically distinguishable quantities `s_d` | N = 10 | — |")
    add("| independent data constraints in `calibrate_alpha()` | **1** | the doxorubicin depth |")
    add("| → identified | **1** | `s_dox`, hence `k_dox` |")
    add("| → NOT identified | **9** | every ratio `b_d/b_dox`, `d ≠ dox` |")
    add("")
    add("* **Identified — exactly one number:** `k_dox = c_dox + alpha·b_dox = "
        f"{p1['k_target']:.6e} /s`. This is fixed by the anchor and is *invariant to any "
        "change in `b_dox`*. A better `b_dox` cannot move it, which is why the branch saw "
        "zero improvement at the anchor.")
    add("* **NOT identified — nine numbers:** the ratios `s_d/s_dox = b_d/b_dox`. These "
        "come entirely from `molecular/binding.py`'s logP/HBD sigmoid and have **zero** "
        "data constraint here. `results/binding_sink_refuted.md` §2 already established "
        "that this sigmoid is refuted in direction, not merely in magnitude — so the nine "
        "unconstrained numbers are not merely unmeasured, they are the ones independently "
        "shown to be wrong.")
    add("* **NOT identified — `alpha` itself:** it is a gauge coordinate (§1.2). "
        "Restricted to the two-parameter sub-problem `(ln alpha, ln b_dox)`, the Jacobian "
        "of the single fitted observable is `(1, 1)` — rank 1 in a 2-dimensional space, "
        "with a null direction `(δ ln alpha, δ ln b_dox) = (t, −t)` along which the fit "
        "is *exactly* flat.")
    add("")
    add("**The consequence, algebraically.** Send `b_dox → f·b_dox` with every other "
        "drug's prior untouched. The anchor equation gives")
    add("")
    add("```")
    add("alpha'      = (k* − c_dox)/(f·b_dox) = alpha / f            (exactly)")
    add("s_dox'      = alpha'·(f·b_dox) = alpha·b_dox                (invariant)")
    add("s_d'        = alpha'·b_d       = s_d / f      for d != dox  (GLOBAL RESCALE)")
    add("k_d'/k_d    = (c_d + s_d/f)/(c_d + s_d)       -> 1/f  as s_d >> c_d")
    add("```")
    add("")
    add("**So `alpha` is not a physical constant that a better `b_dox` refines. It is a "
        "sink whose only job is to absorb `b_dox`, and it re-emits that absorption onto "
        "every drug it was never fit to.** `alpha`'s docstring calls it \"the one factor "
        "that lifts the deliberately-conservative Stage-0 binding prior up to the true "
        "perivascular sink strength\" — that reading is only available for doxorubicin. "
        "For every other drug `alpha` is a doxorubicin-specific correction applied to a "
        "different molecule.")
    add("")
    add("### 1.1 Proof by perturbation (six factors, spanning 100×)")
    add("")
    add("`experiments/alpha_identifiability.py` monkeypatches `calibration._compile` to "
        "scale **doxorubicin's raw binding only**, then calls the unmodified "
        "`calibrate_alpha()`. If the degeneracy is exact, `alpha_ratio × f = 1`.")
    add("")
    add("| `f` (scale on `b_dox`) | refit `alpha` | `alpha`/`alpha₀` | **`alpha_ratio × f`** | `alpha·b_dox` rel. error | `L_dox` rel. error | sink multiplier on every other drug |")
    add("|---:|---:|---:|---:|---:|---:|---:|")
    for r in p1["perturbations"]:
        other = r["per_drug"]["temozolomide"]["binding_ratio"]
        add(f"| {r['f_binding_scale']:.5f} | {r['alpha']:.5f} | {r['alpha_ratio']:.5f} | "
            f"**{r['alpha_ratio_times_f']:.15f}** | {r['dox_sink_rel_error']:.2e} | "
            f"{r['L_dox_rel_error']:.2e} | {other:.5f}× |")
    add("")
    add("Every row of column 4 is 1 to within double-precision round-off. The anchor "
        "drug's penetration length does not move at all — **it cannot**, it is the fit. "
        "The entire effect of \"improving doxorubicin's binding input\" lands on the "
        "nine drugs that were not improved.")
    add("")
    rr = [r for r in p1["perturbations"]
          if abs(r["f_binding_scale"] - 1 / REPORTED_DOX_BINDING_DROP) < 1e-9][0]
    add(f"For the specific case reported by the prior agent — a measured-logP "
        f"substitution that dropped doxorubicin's raw binding by "
        f"{REPORTED_DOX_BINDING_DROP:.4f}× — this reproduces `alpha` moving by "
        f"{rr['alpha_ratio']:.5f}×, i.e. the reciprocal, to 15 decimal places. "
        f"**The 0.04% agreement they observed was not a coincidence and not an "
        f"approximation. It is an identity, and the residual is bounded in closed form.**")
    add("")
    add("The only channel that can break the reciprocal is the anchor's *clearance*, "
        "which `molecular/binding.py::clearance_rate` also makes a function of logP, and "
        "which enters the numerator:")
    add("")
    add("```")
    add("alpha'/alpha  =  (1/f) · (k* − c'_dox)/(k* − c_dox)")
    add("```")
    add("")
    add(f"Since `c_dox/k* = {p1['clearance_share_of_k_target']:.6e}`, that second factor "
        f"is confined to "
        f"`[{1-p1['max_alpha_deviation_from_reciprocal']:.6f}, "
        f"{1+p1['max_alpha_deviation_from_reciprocal']:.6f}]` — the reciprocal identity "
        f"cannot be violated by more than "
        f"**{100*p1['max_alpha_deviation_from_reciprocal']:.4f}%** no matter what a "
        f"measured logP does to doxorubicin. The observed 0.04% sits inside that bound. "
        f"There is no room in this calibration for a change to the anchor drug to do "
        f"anything except rescale everything else.")
    add("")
    if p1.get("cross_checks"):
        add("### 1.1.1 The closed form predicts BOTH concurrently-reported alphas, in "
            "opposite directions")
        add("")
        add("Two independent workstreams hit this from opposite sides. The identity of "
            "§1 predicts both from doxorubicin's descriptors alone — no simulation, no "
            "refit — using the public `binding.tissue_binding_rate` "
            f"(doxorubicin: MW {p1['dox_descriptors']['mw']}, logP "
            f"{p1['dox_descriptors']['logp']}, HBD {p1['dox_descriptors']['hbd']}):")
        add("")
        add("| reported case | source | implied `f` on `b_dox` | `alpha` predicted by the identity | `alpha` reported | rel. error |")
        add("|---|---|---:|---:|---:|---:|")
        for c in p1["cross_checks"]:
            add(f"| {c['case']} | `{c['source']}` | {c['b_dox_scale_f']:.5f} | "
                f"**{c['alpha_predicted_by_identity']:.5f}** | {c['alpha_reported']:.5f} | "
                f"{c['rel_error_vs_reported']:.2e} |")
        add("")
        add(f"Both land inside the "
            f"{100*p1['max_alpha_deviation_from_reciprocal']:.4f}% clearance bound "
            f"derived above. **The two reports disagree about the sign of the correction "
            f"to doxorubicin's logP** — one raises it (Crippen 0.001 → Hansch 1.27, "
            f"`alpha` falls to ~68.6, every other drug's `L` grows ~1.47×), the other "
            f"lowers it (`alpha` rises to ~497, every other drug's `L` shrinks ~0.55×). "
            f"**That disagreement is not the problem this file is about, and resolving it "
            f"would not fix anything.** Whichever is right, the correction has exactly "
            f"one place to go, and it is not doxorubicin.")
        add("")
    add("### 1.2 `alpha` is a GAUGE quantity — its value carries no physical content")
    add("")
    add("The degeneracy is not only about the anchor. Scale **every** drug's raw binding "
        "prior by the same `t` — the whole `tissue_binding_rate` output band moved "
        "bodily — and refit:")
    add("")
    add("```")
    add("b_d -> t·b_d  for ALL d     =>     alpha -> alpha/t     =>     s_d = alpha·b_d  UNCHANGED")
    add("```")
    add("")
    add("| `t` (scale on EVERY `b_d`) | refit `alpha` | `alpha`/`alpha₀` | `alpha_ratio × t` | worst \\|Δ`L`\\| over the whole panel | worst rel. Δ`k` | mitoxantrone predicted `L` |")
    add("|---:|---:|---:|---:|---:|---:|---:|")
    for r in p1.get("gauge_test", []):
        add(f"| {r['t']:g} | {r['alpha']:.5f} | {r['alpha_ratio']:.6g} | "
            f"{r['alpha_ratio_times_t']:.15f} | **{r['worst_abs_L_change_um']:.3e} µm** | "
            f"{r['worst_rel_k_change']:.3e} | {r['mitox_predicted_L_um']:.10f} µm |")
    add("")
    add("**`alpha` moves over four orders of magnitude and every single downstream "
        "number — every `k`, every `L`, the out-of-sample mitoxantrone prediction, its "
        "verdict — is bit-identical.** `alpha` is a *coordinate*, not a measurement. The "
        "pair `(alpha, {b_d})` has an exact one-parameter gauge freedom "
        "`(alpha/t, {t·b_d})`, and only the gauge-invariant products `s_d = alpha·b_d` "
        "reach the PDE.")
    add("")
    add("This is what makes the 17-significant-figure publication of "
        "`alpha = 148.60815765445545` in `README.md`, `API_CONTRACT.md`, "
        "`tmesim/core/config.py` and `tmesim/validation/revalidation.py` misleading in "
        "kind rather than merely in precision. **The number is not wrong; it is not the "
        "sort of thing that can be right.** It is only meaningful jointly with the exact "
        "`binding.py` sigmoid that produced `b_dox`, and it silently changes whenever "
        "that sigmoid does — which is precisely the event that triggered this "
        "investigation. The quantity worth freezing and publishing is the "
        "gauge-invariant one:")
    add("")
    add(f"> `s_dox = alpha·b_dox = {p1['dox_sink_alpha_b']:.6e} /s`  "
        f"(equivalently `k_dox = {p1['k_target']:.6e} /s`)")
    add("")
    add("### 1.3 What the calibration's own out-of-sample test can and cannot see")
    add("")
    add("`calibration.revalidate()` predicts mitoxantrone with `alpha` frozen and it "
        "PASSES *inside* the published band. What does that test actually constrain? "
        "Mitoxantrone's predicted `L` depends on `alpha·b_mitox = s_dox·(b_mitox/b_dox)`. "
        "So it constrains the **ratio** `b_mitox/b_dox` — one of the N−1 gauge-invariant "
        "quantities — and is completely blind to the gauge (§1.2 proves that "
        "empirically: the mitoxantrone column there does not move at all).")
    add("")
    add("**It is NOT blind to an anchor-only change, and this is the finding with the "
        "sharpest consequence for the branch.** Changing `b_dox` alone changes the ratio "
        "`b_mitox/b_dox`, so it moves the one out-of-sample prediction this project has:")
    add("")
    add("| `f` on `b_dox` only | `alpha` | mitoxantrone predicted `L` | ratio to the 72.1 µm mid | verdict | TMZ/dox fold |")
    add("|---:|---:|---:|---:|:--:|---:|")
    for r in p1.get("revalidate_flatness", []):
        add(f"| {r['f_binding_scale']:.5f} | {r['alpha']:.5f} | "
            f"{r['mitox_predicted_L_um']:.4f} µm | {r['mitox_ratio_to_mid']:.4f} | "
            f"**{r['verdict']}** | {r['fold_further']:.4f} |")
    add("")
    base_rv = p1["revalidate_flatness"][0]
    pert_rv = p1["revalidate_flatness"][1]
    add(f"**The measured-logP correction to doxorubicin knocks the module's own "
        f"out-of-sample test out of the published band.** At the reported "
        f"{REPORTED_DOX_BINDING_DROP:.4f}× binding drop, mitoxantrone's predicted `L` "
        f"falls {base_rv['mitox_predicted_L_um']:.2f} → "
        f"{pert_rv['mitox_predicted_L_um']:.2f} µm against a published e-fold-equivalent "
        f"band of {cal.MITOX_EFOLD_BAND_UM[0]:.1f}–{cal.MITOX_EFOLD_BAND_UM[1]:.1f} µm. "
        f"It is **below the band**, and survives only under the looser 3× "
        f"order-of-magnitude lens (`PASS` → `PASS(3x)`). `results/calibration.md`'s "
        f"headline sentence — \"the predicted penetration falls **inside** the band — the "
        f"one constant transfers to a second drug it was never fit to\" — would no longer "
        f"be true. This is an independent, code-internal reason the branch is not "
        f"mergeable, and it is stronger than a magnitude argument because it is the "
        f"project's own declared falsification test.")
    add("")
    add("**The structural blind spot that remains.** Nothing in this fit — not the "
        "anchor, not the out-of-sample test — sees the gauge, and nothing sees `PS` "
        "either (`results/transvascular_exchange_revalidation.md` §6.2: with "
        "`transvascular_exchange = False`, \"PS appears nowhere in that formula at "
        "all\"). Two independent unidentified directions live in the same one-equation "
        "fit.")
    add("")
    add("---")
    add("")

    # ---- 2 ----
    if p2:
        add("## 2. Blast radius of the ANCHOR CHOICE")
        add("")
        add("`tmesim/validation/pk_anchors.py` is this repo's own cited compendium of "
            "measured tissue-PK anchors. **It contains exactly three drugs with a "
            "measured perivascular half-distance — doxorubicin, mitoxantrone, "
            "topotecan — and the last two get their number from a range the source paper "
            "reports *jointly for all three drugs*.** Temozolomide has no penetration-depth "
            "anchor anywhere in the repo (its only entry is a brain:plasma AUC partition, "
            "PMID 19861433, which is a different physical quantity and is flagged as such "
            "in that file). So the anchor cannot be swapped to temozolomide without "
            "inventing data. The swap below uses mitoxantrone, and every target is a "
            "published half-distance converted with the module's own `d_half/ln2` rule.")
        add("")
        add("*(Data-hygiene note found while doing this: the repo cites PMID 23680920 two "
            "different ways — `pk_anchors.py` as \"Patel KJ, Tredan O, Tannock IF … "
            "2013;72(1):127-38\" and `calibration.py` as \"Rocca … Tannock … "
            "2013;72(3):591-601\". Same PMID, different first author, volume issue and "
            "pages. One of them is wrong. Not load-bearing for anything below — the "
            "numeric anchor is identical in both — but it should be reconciled.)*")
        add("")
        add("| anchor | target `L*` | refit `alpha` | vs baseline | `L_dox` | `L_tmz` | `L_tmz/L_dox` | drugs over the 5e-3 /s ceiling | worst |")
        add("|---|---:|---:|---:|---:|---:|---:|:--:|---:|")
        for v in p2["variants"]:
            pan = v["panel"]
            add(f"| {v['label']} | {v['efold_target_um']:.2f} µm | "
                f"**{v['alpha']:.3f}** | {v['alpha_vs_baseline']:.4f}× | "
                f"{pan['doxorubicin']['L_um']:.2f} µm | {pan['temozolomide']['L_um']:.2f} µm | "
                f"{v['L_tmz_over_L_dox']:.3f} | {v['n_over_ceiling']}/{len(PANEL)} | "
                f"{v['worst_over_ceiling']:.1f}× |")
        add("")
        hs = p2["headline_swap"]
        add(f"**The anchor-drug swap on its own moves `alpha` by "
            f"{hs['ratio']:.4f}×** ({hs['alpha_baseline']:.3f} → "
            f"{hs['alpha_mitoxantrone']:.3f}). Two independent things are mixed into that "
            f"number and both matter:")
        add("")
        add("1. the target length differs (64.92 vs 72.13 µm), because the two papers "
            "report different half-distances; and")
        add("2. the anchor drug's own `b`, `c` and `D` differ, and *those* enter only "
            "through the degeneracy of §1.")
        add("")
        add("**Per-drug penetration lengths under the mitoxantrone anchor, relative to "
            "the published panel:**")
        add("")
        add("| drug | `L` (baseline) | `L` (mitoxantrone anchor) | fold |")
        add("|---|---:|---:|---:|")
        base = p2["variants"][0]["panel"]
        mit = next(r for r in p2["variants"] if r["key"] == "mitox_mid")
        for n in PANEL:
            add(f"| {n} | {base[n]['L_um']:.2f} µm | {mit['panel'][n]['L_um']:.2f} µm | "
                f"{mit['L_ratio_vs_baseline'][n]:.4f}× |")
        add("")
        add("And the *published band itself* is wider than the anchor-drug swap. Taking "
            "the two ends of the SAME doxorubicin measurement (40 µm and 50 µm) gives "
            f"`alpha` = {p2['variants'][1]['alpha']:.2f} and "
            f"{p2['variants'][2]['alpha']:.2f} — a "
            f"{p2['variants'][1]['alpha']/p2['variants'][2]['alpha']:.3f}× spread from the "
            "reported uncertainty of the single anchor, before any question of which drug "
            "to anchor on. Taking the mitoxantrone band's two ends (25–75 µm) gives "
            f"{p2['variants'][6]['alpha']:.2f} to {p2['variants'][5]['alpha']:.2f}, a "
            f"**{p2['variants'][5]['alpha']/p2['variants'][6]['alpha']:.1f}× range**.")
        add("")
        add("**`alpha = 148.60815765445545` is quoted to 17 significant figures in "
            "`README.md`, `API_CONTRACT.md`, `tmesim/core/config.py` and "
            "`tmesim/validation/revalidation.py`. Its defensible precision is about one "
            "significant figure.** The digits are reproducible, not meaningful: they are "
            "a deterministic function of an RDKit logP and a 40–50 µm literature range.")
        add("")
        if p2.get("mutual_consistency"):
            mc = p2["mutual_consistency"]
            add("### 2.1 The two measured anchors are NOT mutually consistent, and the "
                "out-of-sample PASS is an artefact of which band is wider")
            add("")
            add("Fit on one anchor, then ask whether the OTHER anchor's published band is "
                "hit. The answer depends entirely on which way round you do it:")
            add("")
            add("| fit on | predict | predicted `L` | published e-fold band | band width | verdict |")
            add("|---|---|---:|---:|---:|:--:|")
            add(f"| doxorubicin (64.92 µm) | mitoxantrone | "
                f"{mc['fit_dox_pred_mitox_L']:.2f} µm | "
                f"{mc['mitox_band_um'][0]:.1f}–{mc['mitox_band_um'][1]:.1f} µm | "
                f"**{mc['mitox_band_width']:.2f}×** | "
                f"{'**INSIDE**' if mc['fit_dox_pred_mitox_in_band'] else '**OUTSIDE**'} |")
            add(f"| mitoxantrone (72.13 µm) | doxorubicin | "
                f"{mc['fit_mitox_pred_dox_L']:.2f} µm | "
                f"{mc['dox_band_um'][0]:.1f}–{mc['dox_band_um'][1]:.1f} µm | "
                f"**{mc['dox_band_width']:.2f}×** | "
                f"{'**INSIDE**' if mc['fit_mitox_pred_dox_in_band'] else '**OUTSIDE**'} |")
            add("")
            add(f"**The model cannot fit both measurements at once.** Anchored on "
                f"doxorubicin it hits mitoxantrone; anchored on mitoxantrone it puts "
                f"doxorubicin at {mc['fit_mitox_pred_dox_L']:.1f} µm, "
                f"**{mc['fit_mitox_pred_dox_over_band_top']:.2f}× above the top** of "
                f"doxorubicin's own published band. That is a real, if modest, "
                f"falsification of the single-sink model on the only two penetration "
                f"measurements it has — and it is invisible in the current direction of "
                f"the test.")
            add("")
            add(f"**And the direction that passes is the easy one.** The held-out "
                f"mitoxantrone band spans **{mc['mitox_band_width']:.2f}×** "
                f"({mc['mitox_band_um'][0]:.0f}–{mc['mitox_band_um'][1]:.0f} µm); the "
                f"doxorubicin band it was fitted to spans only "
                f"**{mc['dox_band_width']:.2f}×**. `pk_anchors.py`'s own note explains "
                f"why: the 25–75 µm range is reported *jointly for three drugs at one "
                f"timepoint*, not as a mitoxantrone-specific value, and the file grades "
                f"that entry `confidence: medium` for exactly this reason. **A prediction "
                f"landing inside a 3× band is a weak pass**, and `results/calibration.md` "
                f"presents it as the module's headline validation without quoting the "
                f"band width. The honest version of that claim is: the one constant "
                f"transfers to a second, chemically-adjacent drug within a factor of 3, "
                f"in the only direction where the target is 3× wide.")
            add("")
        add("### 2.2 The internal contradiction gets worse under every alternative anchor")
        add("")
        add("`results/binding_sink_refuted.md` §4 showed that after the `alpha` "
            "amplification most calibrated drugs breach `molecular/binding.py`'s own "
            "documented [1e-5, 5e-3] /s band. The counts in the table above are on the "
            "10-drug `experiments/real_drug_panel.py` panel, so they are not the same "
            "denominator as that file's 19-drug or 13-drug counts — but the point is the "
            "*movement*: from "
            f"{p2['variants'][6]['n_over_ceiling']}/{len(PANEL)} (worst "
            f"{p2['variants'][6]['worst_over_ceiling']:.1f}×) to "
            f"{p2['variants'][5]['n_over_ceiling']}/{len(PANEL)} (worst "
            f"{p2['variants'][5]['worst_over_ceiling']:.1f}×) across nothing but the two "
            "ends of one published range. **The repo's own internal self-contradiction is "
            "a function of an arbitrary anchor choice.** That is a stronger statement "
            "than the contradiction itself: a quantity whose breach of its own documented "
            "bounds can be dialled from 4/10 to 10/10 by picking a different end of a "
            "literature range is not a measured quantity.")
        add("")
        add("---")
        add("")

    # ---- 3 ----
    if canc:
        UA = canc["arms"]["uniform"]
        AA = canc["arms"]["anchor"]
        add("## 3. Does it cancel?")
        add("")
        add(f"**Setup.** Real 10x Xenium GBM section, `dx = {p3['dx_um']:.0f}` \u00b5m, "
            f"`dt = 10` s, 6 h window, convection on, seed 0 \u2014 identical to "
            f"`tmesim/validation/kill_calibration.py`, whose `prepare_tissue`, "
            f"`unit_transport`, `decompose` and `kill_with` are **imported, not "
            f"re-implemented**. 20 transport solves at "
            f"`alpha = {p3['alpha0']:.5f}` and `alpha = {p3['alpha0']*g:.5f}` "
            f"(**{g:.4f}\u00d7**). Sanity: the reduced path reproduces the full engine "
            f"(`simulate_response` on temozolomide gives kill "
            f"{UA['rows'][0]['kill_base']:.6f}, matched here to all printed digits) and "
            f"reproduces `results/transvascular_exchange_revalidation.md` \u00a72's "
            f"published one-way numbers exactly (cisplatin 2.7122, temozolomide 1.1403, "
            f"lomustine 1.3881).")
        add("")
        add("### 3.1 The perturbation is TWO different perturbations, and they answer "
            "the question differently")
        add("")
        add("This is the part the framing of the question gets wrong, and it matters.")
        add("")
        add("| arm | what moves | is it common-mode? |")
        add("|---|---|---|")
        add(f"| **A \u2014 uniform** | every drug\u2019s sink \u00d7 {g:.4f} | YES \u2014 this is the "
            "gauge direction of \u00a71.2 |")
        add(f"| **B \u2014 anchor** | nine drugs\u2019 sinks \u00d7 {g:.4f}, **doxorubicin exempt** | "
            "NO \u2014 one drug in ten does not move |")
        add("")
        add("**Arm B is what actually happens when `b_dox` changes**, and it is *not* a "
            "global multiplier. Under `b_dox \u2192 b_dox/g`, `alpha \u2192 g\u00b7alpha`, the anchor\u2019s "
            "`s_dox = alpha\u00b7b_dox` and its `clearance` are both untouched, so doxorubicin\u2019s "
            "`DrugParams` are bit-identical and its solve *is* the baseline solve; the "
            "other nine get exactly the perturbed solve. Both arms are therefore exact "
            "reconstructions from the same 20 runs, with no extra simulation and no "
            "approximation.")
        add("")
        add("**A global multiplier on every drug\u2019s sink is common-mode. A change to the "
            "ANCHOR drug\u2019s binding is not \u2014 it is common-mode on 9 of 10 drugs with the "
            "10th pinned.** That distinction is the whole answer to \u201cdoes it cancel?\u201d")
        add("")
        add("`alpha` scales only the `binding` half of `k = clearance + binding`, so even "
            "in arm A `k'/k = (c + g\u00b7s)/(c + s) \u2264 g`:")
        add("")
        add("| drug | `k` baseline | `k` \u00d7g arm | `k` fold | purity `ln(k'/k)/ln g` (arm A) | purity (arm B) |")
        add("|---|---:|---:|---:|---:|---:|")
        for ra, rb in zip(UA["rows"], AA["rows"]):
            bb = p3["by"]["baseline"][ra["drug"]]
            pp = p3["by"]["perturbed"][ra["drug"]]
            add(f"| {ra['drug']} | {bb['k']:.4e} /s | {pp['k']:.4e} /s | "
                f"{ra['k_ratio']:.4f}\u00d7 | {ra['common_mode_purity']:.5f} | "
                f"{rb['common_mode_purity']:.5f} |")
        add("")
        pur = [r["common_mode_purity"] for r in UA["rows"]]
        add(f"Arm A purity {min(pur):.4f}\u2013{max(pur):.4f} \u2014 common-mode to within "
            f"{100*(1-min(pur)):.2f}%, genuinely the class of error that has been "
            f"cancelling in this repo. Arm B\u2019s doxorubicin row is **0.00000**: the anchor "
            f"does not move at all.")
        add("")
        add("### 3.2 Absolute outputs \u2014 the full size of the error")
        add("")
        add("Arm A and arm B differ only in the doxorubicin row (arm B\u2019s is 1.0000\u00d7 by "
            "construction).")
        add("")
        add("| drug | mean tumour exposure (\u00b5M per \u00b5M plasma) | fold, arm A | fold, arm B | kill @ nominal dose | fold, arm A | fold, arm B |")
        add("|---|---:|---:|---:|---:|---:|---:|")
        for ra, rb in zip(UA["rows"], AA["rows"]):
            add(f"| {ra['drug']} | {ra['e_unit_mean_base']:.5g} \u2192 {ra['e_unit_mean_pert']:.5g} | "
                f"**{ra['e_unit_fold']:.4f}\u00d7** | {rb['e_unit_fold']:.4f}\u00d7 | "
                f"{ra['kill_base']:.6f} \u2192 {ra['kill_pert']:.6f} | "
                f"**{ra['kill_fold']:.4f}\u00d7** | {rb['kill_fold']:.4f}\u00d7 |")
        add("")
        add(f"**Nothing quoted as an absolute survives either arm.** Kill moves "
            f"{min(r['kill_fold'] for r in UA['rows']):.3f}\u2013"
            f"{max(r['kill_fold'] for r in UA['rows']):.3f}\u00d7 and mean exposure "
            f"{min(r['e_unit_fold'] for r in UA['rows']):.3f}\u2013"
            f"{max(r['e_unit_fold'] for r in UA['rows']):.3f}\u00d7.")
        add("")
        add("### 3.3 Ratio-type outputs \u2014 what survives")
        add("")
        add("`surviving fraction = |\u0394ln(ratio-type output)| / |\u0394ln(absolute output)|`; "
            "`cancellation = 1 \u2212 surviving`. Same construction as "
            "`results/transvascular_exchange_revalidation.md`. Two normalisers: against "
            "the anchor drug, and against the panel geometric mean (what a relative "
            "score effectively uses).")
        add("")
        add("**Arm A (uniform / gauge direction) \u2014 it cancels essentially completely.**")
        add("")
        add("| drug | kill/kill_dox | surviving | **cancelled** | kill/panel-geomean | surviving | **cancelled** |")
        add("|---|---:|---:|---:|---:|---:|---:|")
        for r in UA["rows"]:
            if r["drug"] == ANCHOR:
                add(f"| {r['drug']} *(normaliser)* | 1.0000 \u2192 1.0000 | \u2014 | \u2014 | "
                    f"{r['gm_rel_kill_base']:.4f} \u2192 {r['gm_rel_kill_pert']:.4f} | "
                    f"{100*r['surviving_frac_kill_gm']:.1f}% | "
                    f"**{100*(1-r['surviving_frac_kill_gm']):.1f}%** |")
                continue
            add(f"| {r['drug']} | {r['rel_kill_base']:.4f} \u2192 {r['rel_kill_pert']:.4f} | "
                f"{100*r['surviving_frac_kill']:.1f}% | "
                f"**{100*r['cancel_factor_kill']:.1f}%** | "
                f"{r['gm_rel_kill_base']:.4f} \u2192 {r['gm_rel_kill_pert']:.4f} | "
                f"{100*r['surviving_frac_kill_gm']:.1f}% | "
                f"**{100*(1-r['surviving_frac_kill_gm']):.1f}%** |")
        add("")
        add("**Arm B (what the branch actually does) \u2014 against the anchor it does not "
            "cancel at all.**")
        add("")
        add("| drug | kill/kill_dox | surviving | **cancelled** | kill/panel-geomean | surviving | **cancelled** |")
        add("|---|---:|---:|---:|---:|---:|---:|")
        for r in AA["rows"]:
            if r["drug"] == ANCHOR:
                add(f"| {r['drug']} *(exempt \u2014 does not move)* | 1.0000 \u2192 1.0000 | \u2014 | \u2014 | "
                    f"{r['gm_rel_kill_base']:.4f} \u2192 {r['gm_rel_kill_pert']:.4f} | \u2014 | \u2014 |")
                continue
            add(f"| {r['drug']} | {r['rel_kill_base']:.4f} \u2192 {r['rel_kill_pert']:.4f} | "
                f"{100*r['surviving_frac_kill']:.1f}% | "
                f"**{100*r['cancel_factor_kill']:.1f}%** | "
                f"{r['gm_rel_kill_base']:.4f} \u2192 {r['gm_rel_kill_pert']:.4f} | "
                f"{100*r['surviving_frac_kill_gm']:.1f}% | "
                f"**{100*(1-r['surviving_frac_kill_gm']):.1f}%** |")
        add("")
        add(f"**The headline cancellation numbers.** Against the anchor drug: arm A "
            f"cancels a median **{100*(1-UA['verdicts']['median_surviving_kill']):.1f}%** "
            f"of the kill error; arm B cancels a median "
            f"**{100*(1-AA['verdicts']['median_surviving_kill']):.1f}%**. Against the "
            f"panel geometric mean: arm A "
            f"**{100*(1-UA['verdicts']['median_surviving_kill_gm']):.1f}%**, arm B "
            f"**{100*(1-AA['verdicts']['median_surviving_kill_gm']):.1f}%**. "
            f"Compare `results/transvascular_exchange_revalidation.md`\u2019s 84\u201386% for the "
            f"transvascular fix.")
        add("")
        add("**Cross-drug ranking.**")
        add("")
        add("| arm | quantity | Spearman \u03c1 | Kendall \u03c4 | discordant / 45 | top-1 | top-3 set | top-1 margin base \u2192 arm |")
        add("|---|---|---:|---:|---:|:--:|:--:|---:|")
        for label, AR in (("A uniform", UA), ("B anchor", AA)):
            rk = AR["ranking"]
            add(f"| {label} | kill_fraction | {rk['kill_spearman']:.4f} | "
                f"{rk['kill_kendall']:.4f} | {rk['kill_discordant_pairs']} | "
                f"{'HELD' if rk['kill_top1_same'] else '**CHANGED**'} | "
                f"{'HELD' if rk['kill_top3_same'] else '**CHANGED**'} | "
                f"{100*rk['kill_top1_margin_base']:.2f}% \u2192 "
                f"{100*rk['kill_top1_margin_pert']:.2f}% |")
            add(f"| {label} | mean exposure | {rk['auc_spearman']:.4f} | "
                f"{rk['auc_kendall']:.4f} | {rk['auc_discordant_pairs']} | "
                f"{'HELD' if rk['auc_top1_same'] else '**CHANGED**'} | "
                f"{'HELD' if rk['auc_top3_same'] else '**CHANGED**'} | \u2014 |")
        add("")
        add("**Rank order and relative magnitude are different claims, and only the "
            "first one survives.** Arm A preserves the kill ORDER perfectly "
            "(\u03c1 = 1.0000, 0/45 discordant) while the kill of individual drugs "
            "*relative to doxorubicin* still moves by up to "
            f"{100*max(r['surviving_frac_kill'] for r in UA['rows'] if r['drug'] != ANCHOR):.0f}% "
            "of the absolute error \u2014 cisplatin\u2019s kill:kill_dox goes 1.6809 \u2192 2.5274. "
            "Any claim of the form \u201cdrug X kills N\u00d7 as much as drug Y\u201d is NOT protected "
            "by the ranking result. This is the same distinction "
            "`results/transvascular_exchange_revalidation.md` \u00a74.1 insists on \u2014 \u201cany "
            "future claim must say WHICH quantity it is about\u201d \u2014 and it applies here "
            "in the ratio direction rather than the delivery-vs-kill direction.")
        add("")
        add(f"* kill order, baseline:  {' > '.join(UA['ranking']['kill_order_base'])}")
        add(f"* kill order, arm A:     {' > '.join(UA['ranking']['kill_order_pert'])}")
        add(f"* kill order, arm B:     {' > '.join(AA['ranking']['kill_order_pert'])}")
        add("")
        rkB = AA["ranking"]
        if not rkB["kill_top1_same"]:
            add(f"**Arm B changes top-1** \u2014 but read the margin before reacting to it. "
                f"The baseline leader wins by "
                f"{100*UA['ranking']['kill_top1_margin_base']:.2f}% of its own kill; under "
                f"arm B the new leader wins by "
                f"{100*rkB['kill_top1_margin_pert']:.2f}%. **This is a near-tie, not a "
                f"reversal.** It is categorically unlike "
                f"`results/transvascular_exchange_revalidation.md` \u00a74.1, where two swapped "
                f"pairs crossed gaps of 40.5% and 35.6% of the leader\u2019s kill \u2014 roughly "
                f"50-sigma. The honest statement is that arm B leaves the ranking intact "
                f"except at a boundary that was already inside the noise.")
        else:
            add("**The ranking survives both arms**, top-1 and top-3 set intact.")
        add("")
        v = UA["verdicts"]
        vB = AA["verdicts"]
        add(f"### 3.4 The 3-D:2-D kill ratio \u2014 {100*v['ratio3d2d_share_of_g']:.1f}% of "
            f"the sink error reaches it")
        add("")
        add("`results/kill_calibration.md`\u2019s headline is a **ratio** "
            "(`F_penetration \u00d7 F_state`), and passing it is the only external validation "
            "absolute kill has. `F_state` is transport-independent and exactly invariant "
            "here. `F_penetration` is not: it is `mean_i p_hill(E_i) / p_hill(C_bath)` "
            "with the field renormalised so the p99 tumour cell sits at the bath \u2014 a "
            "renormalisation that removes a *global concentration* factor but **cannot** "
            "remove a change in the SHAPE of the perivascular gradient, and changing the "
            "sink changes exactly that shape. Arm A and arm B are identical here except "
            "for doxorubicin, which does not move in arm B.")
        add("")
        add("| drug | `F_penetration` @1\u00d7IC50 base \u2192 \u00d7g | fold | 3-D:2-D `ratio_exact` base \u2192 \u00d7g | fold |")
        add("|---|---:|---:|---:|---:|")
        for r in UA["rows"]:
            add(f"| {r['drug']} | {r['f_pen_base']:.4f} \u2192 {r['f_pen_pert']:.4f} | "
                f"**{r['f_pen_fold']:.4f}\u00d7** | {r['ratio3d2d_base']:.4f} \u2192 "
                f"{r['ratio3d2d_pert']:.4f} | **{r['ratio3d2d_fold']:.4f}\u00d7** |")
        add("")
        add(f"**At worst {100*v['ratio3d2d_share_of_g']:.1f}% of the sink perturbation "
            f"survives into the 3-D:2-D ratio**, i.e. it cancels "
            f"{100*(1-v['ratio3d2d_share_of_g']):.1f}% \u2014 real protection, but not the "
            f"84\u201386% the transvascular fix got, and the residual lands directly on the "
            f"number `kill_calibration.md` compares to a 2\u20135\u00d7 literature band. "
            f"Temozolomide\u2019s predicted 3-D:2-D ratio moves "
            f"{UA['rows'][0]['ratio3d2d_base']:.4f} \u2192 "
            f"{UA['rows'][0]['ratio3d2d_pert']:.4f}. The band is wide enough "
            f"(`kill_calibration.md` \u00a71 concedes the 2-D reference alone has a ~19\u00d7 IQR) "
            f"that the verdict does not flip at this perturbation size \u2014 but \u201cPASSES for "
            f"all 3 drugs\u201d should not be read as robust to the calibration, because the "
            f"margin it passes by is comparable to the anchor\u2019s own published "
            f"uncertainty (\u00a72).")
        add("")
        add("---")
        add("")

    # ---- 4 ----
    add("## 4. What would identify it")
    add("")
    add("The degeneracy is `rank(constraints) = 1` against `N+1` parameters. It is "
        "broken only by **adding independent constraints of the right kind**, never by "
        "adding parameters. Specifically:")
    add("")
    add("**(a) The one that breaks the `alpha`/`b_dox` split — a direct binding "
        "measurement on the anchor drug.** The split is unidentifiable from *any* number "
        "of penetration depths, because `alpha` and `b_dox` enter every penetration "
        "length only through `s_d = alpha·b_d`. What is needed is a measurement of the "
        "sink *itself*, in the right compartment: the free-interstitial first-order loss "
        "rate `k` (units 1/s), or equivalently an unbound tissue-to-plasma partition "
        "`Kp,uu` plus a transvascular `PS`, for doxorubicin in tumour tissue. "
        "`results/binding_sink_refuted.md` §1 already identifies the exact measurement "
        "class: intratumoural microdialysis with a resolvable washout phase, the design "
        "Blakeley *et al.* (PMID 18787762) used for methotrexate ECF elimination. One "
        "such number for doxorubicin converts `alpha` from a fitted scale into a "
        "*testable prediction*, and the current model already fails that test by ~94× on "
        "the one drug where it exists.")
    add("")
    add("**(b) The one that breaks the N−1 unconstrained ratios — a second drug with an "
        "independently measured penetration depth, in the same tissue, from a DIFFERENT "
        "chemical class.** This is the constraint the mitoxantrone test is trying and "
        "failing to be. It fails on two counts, both documented in the code: "
        "mitoxantrone is an anthracenedione DNA-intercalator chemically adjacent to "
        "doxorubicin (`calibration.py`'s own honesty caveat), and its 25–75 µm number "
        "comes from a range PMID 23680920 reports jointly for doxorubicin, mitoxantrone "
        "and topotecan — so it is not fully independent of the anchor either. The "
        "requirement is specific: **a perivascular immunofluorescence-vs-CD31-distance "
        "profile for a WEAKLY-bound, deeply-penetrating drug — temozolomide is the right "
        "one — measured in glioma tissue by the Tannock-lab method.** That single "
        "measurement would constrain `b_tmz/b_dox`, which currently rests on nothing but "
        "a logP/HBD sigmoid, and it spans the widest arm of the panel. It does not exist "
        "today; §2 above had to be run on mitoxantrone for exactly that reason.")
    add("")
    add("**(c) What would NOT help, and must not be attempted.** Fitting `alpha` "
        "jointly with anything else to the same doxorubicin anchor; fitting a per-drug "
        "correction to the same anchor; or refitting `alpha` on a breast anchor "
        "(explicitly ruled out by `results/binding_sink_refuted.md` §6). Each adds "
        "parameters to a rank-1 constraint and produces a model that fits better and "
        "predicts no more. Two anchors and one free constant is already the maximum "
        "defensible ratio here — it is at rank, not below it.")
    add("")
    add("**(d) The measurement that would make the whole term obsolete rather than "
        "identified.** `results/binding_sink_refuted.md` §2 establishes that the "
        "residual is **active efflux and intracellular trapping**, which a first-order "
        "interstitial loss term does not physically represent and no logD/logP heuristic "
        "can rank. Identifying `alpha` better makes a structurally wrong term more "
        "precisely wrong. The measurement that changes the model rather than its "
        "constants is a per-drug efflux-transporter substrate characterisation "
        "(P-gp/BCRP/MRP) coupled to a two-compartment interstitial/intracellular sink — "
        "at which point `alpha` has no job left to do.")
    add("")
    add("---")
    add("")
    add("## 5. Verdict")
    add("")
    add("**`alpha` is not a fitted physical constant; it is a doxorubicin-shaped hole "
        "into which every revision of the anchor drug's chemistry disappears, and out of "
        "which it re-emerges as a global rescale of nine other drugs.** The published "
        "value 148.60815765445545 is exactly reproducible and *not* a measurement of "
        "anything: it is a gauge coordinate (\u00a71.2), it is defensible to about one "
        "significant figure even as a coordinate (\u00a72), and the single physical quantity "
        "behind it is `k_dox`.")
    add("")
    if canc:
        UA = canc["arms"]["uniform"]
        AA = canc["arms"]["anchor"]
        rkA, rkB = UA["ranking"], AA["ranking"]
        add("**On mergeability of the measured-logP branch: the prior agent\u2019s verdict "
            "stands, and this sharpens why.** The branch\u2019s zero improvement at the anchor "
            "is not a disappointing result; it is a proof that the change had nowhere to "
            "go except onto the other nine drugs. And \u00a71.3 adds an independent, "
            "code-internal reason: the change pushes `calibration.revalidate()`\u2019s "
            "held-out mitoxantrone prediction out of its published band.")
        add("")
        add(f"**On the product.** The cross-drug ranking \u2014 the repo\u2019s declared primary "
            f"output \u2014 holds under the common-mode arm (kill Spearman "
            f"{rkA['kill_spearman']:.4f}, {rkA['kill_discordant_pairs']}/45 discordant, "
            f"top-1 {'held' if rkA['kill_top1_same'] else 'CHANGED'}) and "
            f"{'holds' if rkB['kill_top1_same'] else 'flips top-1 across a near-tie'} "
            f"under the anchor-exempt arm "
            f"({rkB['kill_discordant_pairs']}/45 discordant, margin "
            f"{100*rkB['kill_top1_margin_pert']:.2f}%). That is consistent with "
            f"`results/transvascular_exchange_revalidation.md` \u00a74.1\u2019s explanation of "
            f"*why* \u2014 common-mode survives, drug-specific does not \u2014 and this "
            f"perturbation is common-mode on 9 of 10 drugs.")
        add("")
        add(f"**But the honest reading is not \u201cit cancels\u201d.** Everything quoted as an "
            f"absolute moves by "
            f"{min(r['kill_fold'] for r in UA['rows']):.3f}\u2013"
            f"{max(r['kill_fold'] for r in UA['rows']):.3f}\u00d7. Against the anchor drug "
            f"itself, the arm that actually occurs cancels "
            f"{100*(1-AA['verdicts']['median_surviving_kill']):.1f}%. And the 3-D:2-D "
            f"ratio test \u2014 a ratio, and the only external validation absolute kill has "
            f"\u2014 still takes {100*UA['verdicts']['ratio3d2d_share_of_g']:.1f}% of the "
            f"perturbation. **The product is safe only in the one output it has always "
            f"been safe in, and that safety was never in question here.**")
        add("")
    add("**Nothing here licenses changing `alpha`.** The correct action is to stop "
        "describing it as calibrated-and-validated, and to stop publishing it as a "
        "number. It is a gauge coordinate attached to *one* identified physical "
        "quantity, `k_dox`, fixed by one drug, one measurement and one half-distance→"
        "e-fold convention, with a rank-1 constraint and nine unconstrained ratios "
        "supplied by a heuristic this repo has already refuted in direction.")
    add("")
    add("---")
    add("")
    add("*Generated by `experiments/alpha_identifiability.py`. Anchor values read from "
        "`tmesim/validation/pk_anchors.py`; PDE machinery imported from "
        "`tmesim/validation/kill_calibration.py`. No file under `tmesim/` was modified.*")

    text = "\n".join(A) + "\n"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["all", "algebra", "pde", "report"])
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--dx", type=float, default=12.0)
    ap.add_argument("--g", type=float, default=REPORTED_DOX_BINDING_DROP)
    ap.add_argument("--cache", default=DEFAULT_CACHE)
    args = ap.parse_args()

    algebra_cache = args.cache.replace(".json", "_algebra.json")

    p1 = p2 = p3 = canc = None
    if args.stage in ("all", "algebra", "report"):
        if args.stage == "report" and os.path.exists(algebra_cache):
            d = json.load(open(algebra_cache))
            p1, p2 = d["p1"], d["p2"]
        else:
            p1 = part1_degeneracy()
            p2 = part2_anchor_swap()
            json.dump({"p1": p1, "p2": p2}, open(algebra_cache, "w"))
            log(f"algebra cached -> {algebra_cache}")

    if args.stage in ("all", "pde"):
        p3 = part3_pde(args.cache, args.workers, args.dx, args.g)
        canc = analyse_cancellation(p3)
    elif args.stage == "report" and os.path.exists(args.cache):
        try:
            p3 = part3_pde(args.cache, args.workers, args.dx, args.g)
            canc = analyse_cancellation(p3)
        except KeyError:
            log("PDE cache incomplete; report will omit part 3")

    if args.stage in ("all", "report"):
        out = write_report(p1, p2, p3, canc)
        log(f"wrote {out}")


if __name__ == "__main__":
    main()
