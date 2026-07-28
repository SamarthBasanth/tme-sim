"""Does adding active efflux improve the MEASURED scorecards? Scored, not assumed.

``./.venv/bin/python -m tmesim.validation.efflux_impact``

What this module is for
-----------------------
:mod:`tmesim.molecular.efflux` implements P-gp / BCRP efflux as a
literature-sourced, zero-free-parameter penetrance factor. Implementing a
mechanism is not evidence for it. This module re-scores the two measured human
datasets the repo already has, WITH the term active, and reports before/after:

1. **The six-drug human glioma partition scorecard** (``results/glioma_regional_pk.md``):
   predicted vs measured tissue:plasma ratio for eribulin, cefazolin,
   methotrexate, temozolomide, levetiracetam and adavosertib. Current shipped
   model: geometric-mean absolute fold error 1.99x, 4/6 within 3x, 6/6 within 10x.
2. **The GIST imatinib tumour:normal partition target** (``results/gist_imatinib.md``):
   n = 6 matched within-patient pairs, measured 0.023-0.716, current model 7.92x
   geometric-mean fold error and structurally floored.

It also runs one **external** check that neither scorecard can provide: the
measured brain-penetration ordering of the three PARP inhibitors, which is a
published contrast between drugs the efflux table covers and which was not used
to build anything here.

The result is a NEGATIVE, and it is reported as one
----------------------------------------------------
Efflux does not improve either scorecard. The term stays **default-OFF** in
``treatment.drugs.delivery_source``. The reasons are specific and are the useful
output of this run:

* On the glioma scorecard the term is dominated by **adavosertib**, whose
  in vitro efflux ratio (14.1, PMID 28928160) is contradicted by a direct human
  measurement of the same drug in the same tissue — human GBM tumour Kp,uu
  median 3.2, because OATP1A2-mediated UPTAKE dominates in the acidic tumour
  microenvironment (Li 2017, same paper). Adding efflux drives an already
  under-predicted drug further down.
* On the GIST target the term **cannot** help, for a structural reason exactly
  parallel to the one the Moncayo explant run established for scalar
  ion-trapping: a scalar per-drug factor cancels out of a within-drug
  tissue RATIO. This is demonstrated numerically below, not asserted.
* Both scorecards, plus the PARP ordering check, point at the same thing the
  GIST run named: the missing mechanism is **uptake** transport, not efflux.

Architecture
------------
Additive only. ``core``, ``physics`` and ``sim.engine`` are untouched. The one
change outside this module and ``molecular/efflux.py`` is a default-off ``efflux``
keyword on ``treatment.drugs.delivery_source``, whose bit-exact backward
compatibility is verified numerically by :func:`backward_compatibility_check`
below and reported in the output.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence

import numpy as np

from ..molecular.efflux import (
    EFFLUX_FORMS,
    EFFLUX_TABLE,
    FORM_KPUU,
    FORM_MONOLAYER,
    PREFER_IN_VITRO,
    PREFER_IN_VIVO,
    efflux_detail,
    efflux_factor,
    substrate_summary,
)
from ..reconstruct.barrier import L_PARACELLULAR_LEAK, effective_penetrance

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(REPO_ROOT, "results")
REPORT_MD = os.path.join(RESULTS_DIR, "efflux_impact.md")

BAND_TIGHT = 3.0
BAND_OOM = 10.0


# ===========================================================================
# 0. Backward compatibility of the one touched function
# ===========================================================================
def backward_compatibility_check() -> Dict:
    """Verify ``delivery_source(efflux=False)`` is BIT-IDENTICAL to the prior code.

    The brief permits a minimal addition to the delivery path only if the default
    reproduces the previous behaviour exactly. "Exactly" here means bit-for-bit
    on the float64 array, not "close" — so the comparison is ``==`` on every
    element, over a grid of penetrances and plasma concentrations, and over both
    the no-barrier and barrier branches.

    The reference is the expression as it stood before the ``efflux`` keyword
    existed, written out here rather than imported, so a future edit to the
    shipped function cannot silently redefine what "unchanged" means.
    """
    from ..core import (CellPopulation, CellState, CellType, Domain, DrugParams,
                        Microenvironment, TissueFields)
    from ..treatment.drugs import delivery_source

    rng = np.random.default_rng(20260727)
    domain = Domain(nx=29, ny=23, dx=11.0)
    fields = TissueFields(ecm_density=domain.ones() * 0.45,
                          vascular_density=rng.random(domain.shape))
    cells = CellPopulation(x=np.array([10.0]), y=np.array([10.0]),
                           cell_type=np.array([CellType.TUMOR]),
                           state=np.array([CellState.VIABLE]),
                           sensitivity=np.array([0.5]))
    micro = Microenvironment(domain=domain, fields=fields, cells=cells)
    vascular = np.asarray(fields.vascular_density, dtype=np.float64)
    disruption_field = rng.random(domain.shape)

    checks: List[Dict] = []
    # Include a drug name that IS in the efflux table with a strong factor, so a
    # leaked application would be impossible to miss.
    for name in ("temozolomide", "eribulin", "imatinib", "not_in_the_table"):
        for pen in (0.0, 0.035, 0.5, 0.9399, 1.0):
            drug = DrugParams(name=name, bbb_penetrance=pen)
            for plasma in (0.0, 1.0, 7.3):
                ref = vascular * pen * float(plasma)
                got = delivery_source(micro, drug, plasma)
                checks.append({"name": name, "pen": pen, "plasma": plasma,
                               "branch": "no_barrier",
                               "identical": bool(np.array_equal(got, ref))})
                got0 = delivery_source(micro, drug, plasma, efflux=False)
                checks.append({"name": name, "pen": pen, "plasma": plasma,
                               "branch": "no_barrier_explicit_false",
                               "identical": bool(np.array_equal(got0, ref))})
                ref_b = vascular * np.clip(pen + L_PARACELLULAR_LEAK * disruption_field,
                                           0.0, 1.0) * float(plasma)
                got_b = delivery_source(micro, drug, plasma,
                                        bbb_disruption=disruption_field)
                checks.append({"name": name, "pen": pen, "plasma": plasma,
                               "branch": "barrier",
                               "identical": bool(np.array_equal(got_b, ref_b))})

    # And confirm the switch DOES bite when asked, so the check is not vacuous.
    drug = DrugParams(name="eribulin", bbb_penetrance=0.5)
    on = delivery_source(micro, drug, 1.0, efflux=True)
    off = delivery_source(micro, drug, 1.0)
    bites = bool(np.all(on <= off) and on.sum() < off.sum() * 0.1)

    n = len(checks)
    n_ok = sum(c["identical"] for c in checks)
    return {"n_comparisons": n, "n_bit_identical": n_ok,
            "all_bit_identical": n_ok == n,
            "switch_has_an_effect_when_on": bites,
            "failures": [c for c in checks if not c["identical"]]}


# ===========================================================================
# 1. The substrate table, as data for the report
# ===========================================================================
def substrate_rows(form: str = FORM_KPUU,
                   prefer: str = PREFER_IN_VITRO) -> List[Dict]:
    """One row per drug in :data:`EFFLUX_TABLE`, with the resolved factor."""
    rows = []
    for name, rec in EFFLUX_TABLE.items():
        d = efflux_detail(name, form=form, prefer=prefer)
        rows.append({
            "drug": name,
            "pgp": rec.pgp, "bcrp": rec.bcrp,
            "er_in_vitro": rec.er_in_vitro,
            "er_pgp": rec.er_pgp, "er_bcrp": rec.er_bcrp,
            "fold_ko": rec.fold_ko,
            "er_used": d["efflux_ratio"], "factor": d["factor"],
            "source": d["source"], "reason": d["reason"],
            "pgp_citation": rec.pgp_citation, "bcrp_citation": rec.bcrp_citation,
            "er_citation": rec.er_citation, "er_assay": rec.er_assay,
            "fold_ko_citation": rec.fold_ko_citation,
            "fold_ko_assay": rec.fold_ko_assay,
            "notes": rec.notes,
            "consistency": rec.consistency(form),
        })
    return rows


# ===========================================================================
# 2. Scorecard 1 -- the six-drug measured human glioma partition
# ===========================================================================
#: Where the efflux factor multiplies.
#:
#: ``"transcellular"`` (the physical default) multiplies only the drug's own
#: penetrance, leaving the regional paracellular leak alone: a broken tight
#: junction is a shunt AROUND the pump, which is why efflux inhibition helps
#: least where the barrier is already open.
#:
#: ``"total"`` multiplies the whole bracket including the leak. It is physically
#: wrong and is run only to show that the verdict does not depend on the choice.
PLACEMENTS = ("transcellular", "total")


def _aggregate(rows: Sequence[Dict], drop: Optional[str] = None) -> Dict:
    """Geometric-mean signed and unsigned fold error, band counts, worst case.

    Mirrors ``glioma_regional_pk.run_scorecard``'s aggregation exactly so the
    "before" column here reproduces the published numbers rather than
    re-deriving them differently.
    """
    rs = [r for r in rows if r["drug"] != drop and math.isfinite(r["log10_fold"])]
    if not rs:
        return {}
    logs = [r["log10_fold"] for r in rs]
    alogs = [abs(x) for x in logs]
    worst = max((abs(r["fold_error"]) if r["fold_error"] >= 1.0
                 else 1.0 / r["fold_error"]) for r in rs)
    return {
        "n": len(rs),
        "gm_signed": 10.0 ** (sum(logs) / len(logs)),
        "gm_abs": 10.0 ** (sum(alogs) / len(alogs)),
        "worst_abs_fold": worst,
        "n_within_3x": sum(1 for r in rs if 1 / BAND_TIGHT <= r["fold_error"] <= BAND_TIGHT),
        "n_within_10x": sum(1 for r in rs if 1 / BAND_OOM <= r["fold_error"] <= BAND_OOM),
    }


def glioma_scorecard(form: str = FORM_KPUU,
                     prefer: str = PREFER_IN_VITRO,
                     placement: str = "transcellular") -> Dict:
    """Re-score the six measured human glioma drugs with and without efflux.

    The BEFORE column is not pasted: it is
    ``glioma_regional_pk.run_scorecard(aggregation="geometric", apply_leak=True)``
    called live, i.e. the model exactly as it currently ships (geometric-mean
    CNS-MPO veto plus the fitted additive paracellular leak, each drug evaluated
    at its own measured regional sample composition).

    The AFTER column applies the efflux factor to the structure-only penetrance
    and then re-adds the leak::

        P_eff = clip( bbb_penetrance * f_efflux + L * D_sample , 0, 1 )

    Nothing else changes: same measured values, same regional weights, same
    aggregation.
    """
    if placement not in PLACEMENTS:
        raise ValueError(f"placement must be one of {PLACEMENTS}")
    from .glioma_regional_pk import run_scorecard

    before = run_scorecard(aggregation="geometric", apply_leak=True)
    rows: List[Dict] = []
    for r in before["rows"]:
        drug = r["drug"]
        struct = float(r["structure_only_partition"])
        disr = float(r["sample_mean_bbb_disruption"])
        meas = float(r["measured_mean"])
        det = efflux_detail(drug, form=form, prefer=prefer)
        f = float(det["factor"])
        if placement == "transcellular":
            pred = float(np.clip(struct * f + L_PARACELLULAR_LEAK * disr, 0.0, 1.0))
        else:
            pred = float(effective_penetrance(struct, disr)) * f
        fold_before = float(r["fold_error"])
        fold_after = pred / meas if meas > 0 else float("inf")
        rows.append({
            "drug": drug, "measured": meas, "method": r["method"],
            "structure_only": struct, "sample_disruption": disr,
            "pred_before": float(r["predicted_partition"]), "fold_before": fold_before,
            "log10_fold_before": math.log10(fold_before) if fold_before > 0 else float("nan"),
            "efflux_factor": f, "efflux_reason": det["reason"],
            "efflux_ratio": det["efflux_ratio"], "efflux_source": det["source"],
            "pred_after": pred, "fold_error": fold_after,
            "log10_fold": math.log10(fold_after) if fold_after > 0 else float("nan"),
        })

    before_rows = [{"drug": r["drug"], "fold_error": r["fold_before"],
                    "log10_fold": r["log10_fold_before"]} for r in rows]
    agg_before = _aggregate(before_rows)
    agg_after = _aggregate(rows)
    # adavosertib is the drug whose measured value (2.46, tissue homogenate) lies
    # ABOVE the model's structural ceiling of 1.0, i.e. it is a known separate
    # defect. Reported ex-adavosertib too -- NOT as the headline.
    return {
        "form": form, "prefer": prefer, "placement": placement,
        "rows": rows,
        "before": agg_before, "after": agg_after,
        "before_ex_adavosertib": _aggregate(before_rows, drop="adavosertib"),
        "after_ex_adavosertib": _aggregate(rows, drop="adavosertib"),
        "n_drugs_moved": sum(1 for r in rows if r["efflux_factor"] < 1.0),
        "improved": agg_after["gm_abs"] < agg_before["gm_abs"],
    }


def glioma_paired_test(sc: Dict) -> Dict:
    """Is the before/after difference distinguishable from nothing at n = 6?

    A paired comparison of ``|log10 fold error|`` per drug, by an exact sign test
    over the drugs that actually move. With only a handful of movers this is
    reported to show that it is UNPOWERED, not to claim a result either way.
    """
    from scipy import stats
    b = np.array([abs(r["log10_fold_before"]) for r in sc["rows"]])
    a = np.array([abs(r["log10_fold"]) for r in sc["rows"]])
    d = a - b
    movers = int(np.count_nonzero(np.abs(d) > 1e-12))
    n_worse = int(np.sum(d > 1e-12))
    p_sign = (float(stats.binomtest(n_worse, movers, 0.5).pvalue)
              if movers else float("nan"))
    import warnings
    try:
        with warnings.catch_warnings():
            # scipy warns that n = 6 with 3 ties is too small for its normal
            # approximation. That warning IS the finding; it is carried in
            # `note` below rather than emitted as console noise.
            warnings.simplefilter("ignore")
            p_w = float(stats.wilcoxon(b, a, zero_method="wilcox").pvalue)
    except ValueError:
        p_w = float("nan")
    return {"n_drugs": len(b), "n_movers": movers, "n_worse": n_worse,
            "n_better": movers - n_worse,
            "mean_delta_log10": float(d.mean()), "p_sign_test": p_sign,
            "p_wilcoxon": p_w,
            "note": ("with %d of %d drugs moving at all, neither test can reach "
                     "significance at any effect size; both are reported to make "
                     "that explicit" % (movers, len(b)))}


# ===========================================================================
# 3. Scorecard 2 -- the GIST tumour:normal target
# ===========================================================================
def gist_impact() -> Dict:
    """Can efflux move the measured GIST tumour:normal imatinib ratio? No, and why.

    Four placements are evaluated against the same n = 6 matched pairs. Only the
    first is what a per-drug efflux property actually IS; the rest exist to
    quantify what a *tissue-differential* transporter term would have to look
    like, and to show that the required asymmetry runs opposite to the known
    biology.

    ==========================  ==================================================
    placement                   what it asserts
    ==========================  ==================================================
    ``both`` (the honest one)   imatinib is effluxed in tumour and normal alike
    ``normal_only``             efflux only in normal tissue -- which IS the real
                                hepatic biology (canalicular ABCB1/ABCG2 export
                                imatinib into bile) and moves the prediction the
                                WRONG WAY
    ``tumour_only``             efflux only in tumour -- would improve the number,
                                requires tumour P-gp >> normal liver P-gp, which
                                is backwards: liver is among the highest
                                ABCB1/ABCG2-expressing tissues in the body
    ``required``                solves for the tumour/normal factor ratio that
                                would reproduce the measurement
    ==========================  ==================================================
    """
    from .gist_imatinib import (HALF_SPACING_UM, compile_imatinib,
                                delivery_efficiency, measured_targets)

    targets = measured_targets()
    compiled = compile_imatinib()
    f = efflux_factor("imatinib")
    det = efflux_detail("imatinib")

    e_t = delivery_efficiency(compiled["L_um"], HALF_SPACING_UM["tumour"])
    e_n = {t: delivery_efficiency(compiled["L_um"], HALF_SPACING_UM[t])
           for t in ("liver", "intestine")}
    key = {"Liver": "liver", "Intestine": "intestine"}

    def _score(f_tumour: float, f_normal: float) -> Dict:
        folds, preds = [], []
        for _, r in targets.iterrows():
            p = (f_tumour * e_t) / (f_normal * e_n[key[r.organ]])
            preds.append(p)
            folds.append(max(p / r.ratio_mean, r.ratio_mean / p))
        folds = np.asarray(folds)
        return {"geomean_fold_error": float(10 ** np.mean(np.log10(folds))),
                "max_fold_error": float(folds.max()),
                "n_within_3x": int((folds <= 3.0).sum()),
                "pred_liver": float((f_tumour * e_t) / (f_normal * e_n["liver"])),
                "pred_intestine": float((f_tumour * e_t) / (f_normal * e_n["intestine"])),
                "per_patient": [{"patient": r.patient, "organ": r.organ,
                                 "measured": float(r.ratio_mean),
                                 "pred": float(p), "fold": float(fo)}
                                for (_, r), p, fo in zip(targets.iterrows(), preds, folds)]}

    baseline = _score(1.0, 1.0)                     # == M1, no efflux anywhere
    both = _score(f, f)                             # scalar per-drug property
    normal_only = _score(1.0, f)                    # the real hepatic biology
    tumour_only = _score(f, 1.0)                    # would need tumour >> liver P-gp

    # Exact cancellation, to machine precision, over every pair.
    cancels = all(abs(a["pred"] - b["pred"]) <= 1e-15 * max(1.0, abs(a["pred"]))
                  for a, b in zip(baseline["per_patient"], both["per_patient"]))

    liver = targets[targets.organ == "Liver"]
    meas_liver_gm = float(10 ** np.mean(np.log10(liver.ratio_mean.values)))
    required_ratio = meas_liver_gm / baseline["pred_liver"]
    # what combined ER asymmetry that implies, under FORM_KPUU (f = 1/ER)
    required_er_asymmetry = 1.0 / required_ratio

    return {
        "imatinib_efflux_factor": f,
        "imatinib_efflux_ratio": det["efflux_ratio"],
        "imatinib_efflux_source": det["source"],
        "M1_baseline": baseline,
        "scalar_both_sides": both,
        "scalar_cancels_exactly": bool(cancels),
        "normal_only": normal_only,
        "tumour_only": tumour_only,
        "measured_liver_geomean": meas_liver_gm,
        "required_tumour_over_normal_factor": required_ratio,
        "required_er_asymmetry": required_er_asymmetry,
        "n": int(len(targets)),
    }


# ===========================================================================
# 4. External check -- the PARP-inhibitor brain-penetration ordering
# ===========================================================================
#: Measured brain-to-plasma AUC ratios for the PARP inhibitors, from a source
#: used NOWHERE else in this repo and not used to build the efflux table.
#: Sun et al. 2018 Oncotarget, PMID 30647846: niraparib 0.29-0.33 vs olaparib
#: ~0.03 across two xenograft models. Rucaparib is not in that comparison, so it
#: carries no measured value here.
PARP_MEASURED_BRAIN_PLASMA: Dict[str, Optional[float]] = {
    "niraparib": 0.31,     # midpoint of the reported 0.29-0.33
    "olaparib": 0.03,
    "rucaparib": None,
}
PARP_MEASURED_CITATION = ("Sun et al. 2018 Oncotarget, PMID 30647846 "
                          "(brain:plasma AUC ratio, two xenograft models)")


def parp_ordering_check(form: str = FORM_KPUU) -> Dict:
    """Does the efflux term reproduce a measured between-drug brain-penetration gap?

    This is the only test here that is fully EXTERNAL: neither scorecard contains
    a PARP inhibitor, and the measured contrast (PMID 30647846) played no part in
    building :data:`EFFLUX_TABLE`. Niraparib's brain:plasma AUC ratio is ~10x
    olaparib's; the question is whether the efflux factors alone reproduce that.
    """
    facs = {d: efflux_factor(d, form=form) for d in PARP_MEASURED_BRAIN_PLASMA}
    n, o = facs["niraparib"], facs["olaparib"]
    meas_n = PARP_MEASURED_BRAIN_PLASMA["niraparib"]
    meas_o = PARP_MEASURED_BRAIN_PLASMA["olaparib"]
    return {
        "factors": facs,
        "measured": dict(PARP_MEASURED_BRAIN_PLASMA),
        "citation": PARP_MEASURED_CITATION,
        "measured_niraparib_over_olaparib": meas_n / meas_o,
        "predicted_niraparib_over_olaparib": n / o,
        "ordering_correct": bool(n > o),
        "separation_shortfall_fold": (meas_n / meas_o) / (n / o),
    }


# ===========================================================================
# 5. Sensitivity -- does the verdict depend on any choice we made?
# ===========================================================================
def sensitivity() -> List[Dict]:
    """Re-run the glioma scorecard across every discretionary choice in the term.

    Functional form (``kpuu`` / ``monolayer``), which published observable wins
    when a drug has both (``in_vitro`` / ``in_vivo``), and where the factor
    multiplies (``transcellular`` / ``total``). If the verdict were sensitive to
    any of these, the verdict would not be a verdict.
    """
    out = []
    for form in EFFLUX_FORMS:
        for prefer in (PREFER_IN_VITRO, PREFER_IN_VIVO):
            for placement in PLACEMENTS:
                sc = glioma_scorecard(form=form, prefer=prefer, placement=placement)
                out.append({"form": form, "prefer": prefer, "placement": placement,
                            "gm_abs_before": sc["before"]["gm_abs"],
                            "gm_abs_after": sc["after"]["gm_abs"],
                            "within_3x_before": sc["before"]["n_within_3x"],
                            "within_3x_after": sc["after"]["n_within_3x"],
                            "within_10x_before": sc["before"]["n_within_10x"],
                            "within_10x_after": sc["after"]["n_within_10x"],
                            "worst_before": sc["before"]["worst_abs_fold"],
                            "worst_after": sc["after"]["worst_abs_fold"],
                            "improved": sc["improved"]})
    return out


def per_drug_leave_one_out(sc: Dict) -> List[Dict]:
    """Which drug is responsible? Aggregate with each single drug removed."""
    rows = sc["rows"]
    before_rows = [{"drug": r["drug"], "fold_error": r["fold_before"],
                    "log10_fold": r["log10_fold_before"]} for r in rows]
    out = []
    for r in rows:
        b = _aggregate(before_rows, drop=r["drug"])
        a = _aggregate(rows, drop=r["drug"])
        out.append({"dropped": r["drug"], "gm_abs_before": b["gm_abs"],
                    "gm_abs_after": a["gm_abs"],
                    "improved_without_it": a["gm_abs"] < b["gm_abs"]})
    return out


# ===========================================================================
# 6. Report
# ===========================================================================
def _fmt(v, spec="%.3f", dash="--"):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return dash
    return spec % v


def write_report(compat: Dict, subs: List[Dict], sc: Dict, paired: Dict,
                 loo: List[Dict], gist: Dict, parp: Dict, sens: List[Dict],
                 path: str = REPORT_MD) -> str:
    A: List[str] = []
    a = A.append
    s = substrate_summary()
    b, af = sc["before"], sc["after"]
    bx, ax = sc["before_ex_adavosertib"], sc["after_ex_adavosertib"]

    a("# Active efflux (P-gp / BCRP) — implemented, cited, and it does NOT help")
    a("")
    a("**Verdict up front: implemented as a zero-free-parameter, fully cited mechanism, "
      "scored against both measured human datasets, and left DEFAULT-OFF because it "
      "makes the one scorecard it can act on worse and cannot act on the other at all.**")
    a("")
    a("| scorecard | before | after efflux | verdict |")
    a("|---|---|---|---|")
    a(f"| human glioma partition, n = 6 drugs (GM abs fold error) | "
      f"**{b['gm_abs']:.3f}x** | **{af['gm_abs']:.3f}x** | "
      f"{'improved' if sc['improved'] else 'WORSE'} |")
    a(f"| .. within 3x | {b['n_within_3x']}/6 | {af['n_within_3x']}/6 | "
      f"{'unchanged' if b['n_within_3x'] == af['n_within_3x'] else 'changed'} |")
    a(f"| .. within 10x | {b['n_within_10x']}/6 | {af['n_within_10x']}/6 | "
      f"{'unchanged' if b['n_within_10x'] == af['n_within_10x'] else 'WORSE'} |")
    a(f"| .. worst case | {b['worst_abs_fold']:.2f}x | {af['worst_abs_fold']:.2f}x | "
      f"{'WORSE' if af['worst_abs_fold'] > b['worst_abs_fold'] else 'better'} |")
    a(f"| GIST tumour:normal, n = 6 pairs (geomean fold error) | "
      f"**{gist['M1_baseline']['geomean_fold_error']:.2f}x** | "
      f"**{gist['scalar_both_sides']['geomean_fold_error']:.2f}x** | "
      f"IDENTICAL — the term cancels exactly |")
    a("")
    a("Three findings, in order of how much they are worth:")
    a("")
    a("1. **A scalar per-drug efflux factor is structurally incapable of moving a "
      "tumour:normal ratio.** It cancels: `R = (f·E_t)/(f·E_n) = E_t/E_n`. Verified to "
      "machine precision on all six GIST pairs (§5). This is the SAME class of "
      "structural limit the Moncayo explant run established for scalar ion-trapping, "
      "arrived at independently, and it means the GIST diagnosis cannot be answered by "
      "any per-drug property — only by a per-TISSUE transporter field.")
    a("2. **On the glioma scorecard the term is contradicted by a direct human "
      "measurement of the same drug in the same tissue.** Adavosertib carries the "
      "best-evidenced efflux ratio in the table (14.1, two assay types, two groups), "
      "and the same paper reports human GBM tumour Kp,uu of median 3.2 — adavosertib "
      "ACCUMULATES in human glioma because OATP1A2 UPTAKE dominates in the acidic "
      "tumour microenvironment. Applying the efflux term drives it from 0.24x to "
      f"{[r for r in sc['rows'] if r['drug'] == 'adavosertib'][0]['fold_error']:.3f}x.")
    a("3. **Literature coverage is the binding constraint, not the physics.** Of "
      f"{s['n_drugs']} drugs, {s['substrate_with_number']} carry a usable published "
      f"number and {s['substrate_without_number']} are CONFIRMED substrates with no "
      "usable ratio anywhere — including doxorubicin, the textbook P-gp substrate and "
      "the molecule this repo's perivascular constant was fitted to. A mechanism that "
      "silently does nothing for half the library is not ready to be on by default.")
    a("")

    # ---------------- 1. backward compatibility ----------------
    a("## 1. The one touched function is bit-identical when off")
    a("")
    a("`treatment.drugs.delivery_source` gained a single keyword, `efflux=False`. "
      "`core`, `physics` and `sim.engine` are untouched.")
    a("")
    a(f"- **{compat['n_bit_identical']}/{compat['n_comparisons']} comparisons "
      f"bit-for-bit identical** (`np.array_equal`, not `allclose`) across 4 drug "
      "names x 5 penetrances x 3 plasma concentrations x 3 branches "
      "(implicit default, explicit `efflux=False`, and the barrier branch).")
    a(f"- All identical: **{compat['all_bit_identical']}**.")
    a(f"- And the check is not vacuous: with `efflux=True` on a strong substrate the "
      f"source does change ({compat['switch_has_an_effect_when_on']}).")
    a("")

    # ---------------- 2. the table ----------------
    a("## 2. The literature-sourced substrate table")
    a("")
    a(f"{s['n_drugs']} drugs: all ten in `molecular.library.GBM_DRUGS`, the four "
      "glioma-scorecard drugs outside it, imatinib (GIST), and the three PARP "
      "inhibitors (Moncayo). Every non-`unknown` cell carries a PMID. `unknown` means "
      "searched and not found, and is treated downstream as a non-substrate — but is "
      "reported separately, because *measured not to be a substrate* and *nobody has "
      "looked* are different states of knowledge.")
    a("")
    a(f"- P-gp (ABCB1): **{s['pgp_yes']} yes / {s['pgp_no']} no / "
      f"{s['pgp_unknown']} unknown**")
    a(f"- BCRP (ABCG2): **{s['bcrp_yes']} yes / {s['bcrp_no']} no / "
      f"{s['bcrp_unknown']} unknown**")
    a(f"- Confirmed substrate WITH a usable published number: "
      f"**{s['substrate_with_number']}**")
    a(f"- Confirmed substrate WITHOUT one: **{s['substrate_without_number']}** "
      "(factor forced to exactly 1.0 — we do not invent ratios)")
    a(f"- Measured non-substrate: **{s['non_substrate']}**")
    a("")
    a("| drug | P-gp | BCRP | in vitro ER | in vivo KO fold | ER used | **factor** | source | status |")
    a("|---|---|---|---|---|---|---|---|---|")
    for r in subs:
        a(f"| {r['drug']} | {r['pgp']} | {r['bcrp']} | "
          f"{_fmt(r['er_in_vitro'], '%.1f')} | {_fmt(r['fold_ko'], '%.1f')} | "
          f"{_fmt(r['er_used'], '%.2f')} | **{r['factor']:.4f}** | "
          f"{r['source'] or '--'} | {r['reason']} |")
    a("")
    a("### Citations, per assignment")
    a("")
    for r in subs:
        a(f"**{r['drug']}** — P-gp `{r['pgp']}`"
          + (f": {r['pgp_citation']}" if r['pgp_citation'] else " (no source found)"))
        a(f"  BCRP `{r['bcrp']}`"
          + (f": {r['bcrp_citation']}" if r['bcrp_citation'] else " (no source found)"))
        if r["er_in_vitro"] is not None:
            a(f"  In vitro ER: P-gp {_fmt(r['er_pgp'], '%.1f')}, "
              f"BCRP {_fmt(r['er_bcrp'], '%.1f')}, combined "
              f"{r['er_in_vitro']:.2f}. {r['er_assay']}. {r['er_citation']}")
        if r["fold_ko"] is not None:
            a(f"  In vivo knockout fold: {r['fold_ko']:.1f}. {r['fold_ko_assay']}. "
              f"{r['fold_ko_citation']}")
        a(f"  *Notes:* {r['notes']}")
        a("")

    # ---------------- 3. the form ----------------
    a("## 3. The functional form, and why it is not a knob")
    a("")
    a("Both candidate forms parameterise the same quantity, the efflux burden "
      "`b = CL_efflux / P_passive`, and both give `ER = 1 + b`:")
    a("")
    a("- **`kpuu` (default): `f = 1/ER`.** Single-barrier unbound-partition relation — "
      "efflux clearance in parallel with passive permeability, which is the derivation "
      "behind the published definition `ER = (P_pass + Vmax/Km)/P_pass` itself "
      "(PMID 23421687). The `+1` is already inside ER, so dividing by ER does not "
      "double-count.")
    a("- **`monolayer`: `f = 2/(1+ER)`.** Resolves the cell interior (apical / cell / "
      "basolateral, pump at the apical face); exactly 2x more permissive at high ER.")
    a("")
    a("Both return exactly 1.0 at ER = 1, which is why the plain `1/(1+ER)` starting "
      "form was rejected: it would silently halve every molecule with no efflux at all.")
    a("")
    a("`kpuu` is the default on three grounds, all pre-committed: it is the form with "
      "published empirical support against measured brain `Kp,uu` (r^2 = 0.813 for a "
      "high-expressing MDCK-MDR1 line, PMID 24009309; transwell-to-knockout agreement "
      "over 32 CNS compounds, PMID 17962372; `Kp,uu` inversely tracking ER above "
      "ER ~ 3, PMID 35156331); it does not double-count the `+1`; and it wins the one "
      "internal consistency test available. **An in vivo knockout fold change is the "
      "same model read differently** — deleting the pump sets `b = 0`, so "
      "`f = 1/fold_KO` in EITHER form — which is why using whichever observable a drug "
      "happens to have is not form-shopping.")
    a("")
    a("### The consistency test (n = 2 — a tie-breaker, not a calibration)")
    a("")
    a("Two drugs carry both observables. Under `kpuu` an in vitro ER predicts "
      "`fold_KO = ER`; under `monolayer`, `(1+ER)/2`.")
    a("")
    a("| drug | in vitro ER | measured KO fold | `kpuu` predicts | `monolayer` predicts |")
    a("|---|---|---|---|---|")
    for r in subs:
        c_k = r["consistency"]
        if not c_k:
            continue
        rec = EFFLUX_TABLE[r["drug"]]
        c_m = rec.consistency(FORM_MONOLAYER)
        a(f"| {r['drug']} | {rec.er_in_vitro:.1f} | {rec.fold_ko:.1f} | "
          f"{c_k['fold_ko_predicted_from_er']:.1f} "
          f"({c_k['in_vivo_over_in_vitro']:.2f}x off) | "
          f"{c_m['fold_ko_predicted_from_er']:.1f} "
          f"({c_m['in_vivo_over_in_vitro']:.2f}x off) |")
    a("")
    a("Geometric-mean absolute discrepancy: **1.57x for `kpuu`, 2.15x for `monolayer`**. "
      "n = 2. Reported because it points the same way as the published correlations, "
      "not because two drugs settle anything.")
    a("")
    a("**A scaling factor was available and was refused.** The best-supported practical "
      "formulation in the literature is `Kp,uu = 1/(1 + SF·(ER-1))` with an assay- and "
      "species-specific `SF` calibrated to neuro-PK data (PMID 30683809; human BBB P-gp "
      "is 2.3x LOWER than mouse mdr1a and human BCRP 1.9x higher, PMID 21291474). There "
      "is exactly one dataset here we could calibrate `SF` on — the six-drug scorecard "
      "this term is about to be tested against — so fitting it would have made the test "
      "circular. The shipped term has **zero free parameters** and is correspondingly "
      "allowed to fail. It did.")
    a("")

    # ---------------- 4. glioma scorecard ----------------
    a("## 4. Scorecard 1 — the six-drug measured human glioma partition")
    a("")
    a("BEFORE is not pasted: it is `glioma_regional_pk.run_scorecard(aggregation="
      "\"geometric\", apply_leak=True)` called live, i.e. the model exactly as it ships. "
      "AFTER changes one thing:")
    a("")
    a("```")
    a("  before:  P_eff = clip( bbb_penetrance              + L * D_sample , 0, 1 )")
    a("  after :  P_eff = clip( bbb_penetrance * f_efflux   + L * D_sample , 0, 1 )")
    a("```")
    a("")
    a("Efflux multiplies the TRANSCELLULAR term only — a disrupted tight junction is a "
      "paracellular shunt around the pump. (The physically-wrong alternative, "
      "multiplying the whole bracket, is run in §7 and changes nothing.)")
    a("")
    a("| drug | measured | modality | struct. penetrance | disrupted-sample frac | "
      "**f_efflux** | why | pred before | fold before | pred after | **fold after** |")
    a("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sc["rows"]:
        a(f"| {r['drug']} | {r['measured']:.4f} | {r['method']} | "
          f"{r['structure_only']:.4f} | {r['sample_disruption']:.2f} | "
          f"**{r['efflux_factor']:.4f}** | {r['efflux_reason']} | "
          f"{r['pred_before']:.4f} | {r['fold_before']:.3f}x | "
          f"{r['pred_after']:.4f} | **{r['fold_error']:.3f}x** |")
    a("")
    a("| aggregate | before | after | direction |")
    a("|---|---|---|---|")
    a(f"| GM \\|fold error\\| | **{b['gm_abs']:.3f}x** | **{af['gm_abs']:.3f}x** | "
      f"{'better' if af['gm_abs'] < b['gm_abs'] else 'WORSE'} |")
    a(f"| GM signed bias | {b['gm_signed']:.3f}x | {af['gm_signed']:.3f}x | |")
    a(f"| within 3x | {b['n_within_3x']}/6 | {af['n_within_3x']}/6 | |")
    a(f"| within 10x | {b['n_within_10x']}/6 | {af['n_within_10x']}/6 | "
      f"{'WORSE' if af['n_within_10x'] < b['n_within_10x'] else ''} |")
    a(f"| worst \\|fold error\\| | {b['worst_abs_fold']:.2f}x | "
      f"{af['worst_abs_fold']:.2f}x | "
      f"{'WORSE' if af['worst_abs_fold'] > b['worst_abs_fold'] else 'better'} |")
    a(f"| GM \\|fold error\\|, ex-adavosertib | {bx['gm_abs']:.3f}x | "
      f"{ax['gm_abs']:.3f}x | "
      f"{'better' if ax['gm_abs'] < bx['gm_abs'] else 'WORSE'} |")
    a("")
    a(f"**{sc['n_drugs_moved']} of 6 drugs move at all.** The other "
      f"{6 - sc['n_drugs_moved']} get a factor of exactly 1.0 — not because they are "
      "non-substrates, but because the literature has no usable number for them "
      "(methotrexate and cefazolin are confirmed BCRP substrates with no applicable "
      "ratio; levetiracetam's P-gp status is genuinely disputed).")
    a("")
    a("### Which drug is doing this")
    a("")
    a("| drug dropped | GM \\|fold\\| before | GM \\|fold\\| after | efflux helps without it? |")
    a("|---|---|---|---|")
    for r in loo:
        a(f"| {r['dropped']} | {r['gm_abs_before']:.3f}x | {r['gm_abs_after']:.3f}x | "
          f"{'YES' if r['improved_without_it'] else 'no'} |")
    a("")
    a("The term only helps once adavosertib is removed — and even then it buys "
      f"{bx['gm_abs']:.3f}x -> {ax['gm_abs']:.3f}x, a "
      f"{100 * (1 - ax['gm_abs'] / bx['gm_abs']):.0f}% reduction in error on five drugs. "
      "**Dropping the drug that a mechanism gets wrong, and then reporting that the "
      "mechanism works, is exactly the move this repo exists not to make.** The "
      "ex-adavosertib figure is reported because adavosertib has an independently "
      "documented defect (its measured 2.46 is a tissue-homogenate ratio above the "
      "model's structural ceiling of 1.0), not because it licenses the headline.")
    a("")
    a("### Per drug, honestly")
    a("")
    for r in sc["rows"]:
        if r["efflux_factor"] >= 1.0:
            continue
        direction = ("improves" if abs(math.log10(r["fold_error"]))
                     < abs(r["log10_fold_before"]) else "WORSENS")
        a(f"- **{r['drug']}** ({direction}): factor {r['efflux_factor']:.4f} from "
          f"ER {r['efflux_ratio']:.1f} ({r['efflux_source']}). "
          f"{r['fold_before']:.2f}x -> {r['fold_error']:.3f}x.")
    a("")
    a("  - *eribulin* is the drug this exercise was for, and efflux does move it the "
      "right way — but past the target. It was over-predicted 6.94x; with a 30x mouse "
      "knockout factor it is now under-predicted, so the error shrinks without the "
      "prediction becoming good, and it still fails the 3x band.")
    a("  - *temozolomide* is a genuine surprise from the literature pass: TMZ IS a weak "
      "efflux substrate (a good Bcrp1 substrate, PMID 29852323), not a non-substrate as "
      "commonly assumed. But the measured combined-knockout effect is only 1.5x, and "
      "TMZ was already slightly under-predicted, so even that tiny factor hurts.")
    a("  - *adavosertib* is the whole story. See §4.1.")
    a("")
    a("### 4.1 Why adavosertib breaks it — and why that is a finding, not bad luck")
    a("")
    a("Adavosertib has the best-evidenced efflux data in the table: in vitro ER 10.6 "
      "(MDCKII-ABCB1) and 4.5 (MDCKII-ABCG2) with proper inhibitor controls and 85-115% "
      "mass recovery (PMID 28928160), independently confirmed with a 25x knockout brain "
      "fold change by a second group (PMID 29147815). If the efflux term is going to "
      "work anywhere, it should work here.")
    a("")
    a("It fails, and **the same paper says why**: Li et al. measured human GBM tumour "
      "`Kp,uu` at a **median of 3.2 (range 1.3-24.4)**. Adavosertib does not merely "
      "reach human glioma despite dual efflux — it *accumulates* there, because "
      "OATP1A2-mediated **uptake** dominates in the acidic tumour microenvironment. The "
      "authors state the mouse and human results disagree. Two further details point "
      "the same way: the ERs are pH-dependent and roughly halve as basolateral pH falls "
      "7.4 -> 6.5 (i.e. tumour pH suppresses the efflux the model is applying), and the "
      "measured value the scorecard grades against (2.46) is itself above 1.0.")
    a("")
    a("So the model over-corrects a drug whose real behaviour in the target tissue is "
      "the opposite of the mechanism being added. That is not a tuning failure; it is "
      "the efflux term being **the wrong mechanism for this drug in this tissue**.")
    a("")
    a("### 4.2 Is the change even distinguishable at n = 6?")
    a("")
    a(f"- drugs that move at all: **{paired['n_movers']} of {paired['n_drugs']}**; "
      f"of those, {paired['n_worse']} worsen and {paired['n_better']} improve")
    a(f"- mean change in `|log10 fold error|`: **{paired['mean_delta_log10']:+.4f}** "
      "(positive = worse)")
    a(f"- exact sign test over the movers: p = {_fmt(paired['p_sign_test'], '%.3f')}; "
      f"Wilcoxon signed-rank: p = {_fmt(paired['p_wilcoxon'], '%.3f')}")
    a("")
    a(f"**{paired['note']}.** The honest reading is that six drugs cannot establish "
      "that efflux helps OR that it hurts. What they can do is show that it does not "
      "deliver the large, obvious improvement that would justify enabling an untested "
      "term — and that the one drug with the best efflux data moves the furthest in the "
      "wrong direction for a reason the literature already documented.")
    a("")

    # ---------------- 5. GIST ----------------
    a("## 5. Scorecard 2 — GIST tumour:normal, where the term cannot even apply")
    a("")
    a("`results/gist_imatinib.md` diagnosed the missing mechanism as being on the "
      "**normal-liver** side: active hepatic uptake (OCT1/OATP) plus lysosomal trapping "
      "concentrating imatinib in normal liver. The question here is whether efflux, the "
      "tractable half of that diagnosis, moves the prediction toward the measured "
      "0.023-0.142.")
    a("")
    a(f"Imatinib's factor is **{gist['imatinib_efflux_factor']:.4f}** "
      f"(ER {gist['imatinib_efflux_ratio']:.0f}, {gist['imatinib_efflux_source']}, "
      "PMID 12604685) — one of the strongest in the table. It changes nothing.")
    a("")
    a("### 5.1 The structural result")
    a("")
    a("The prediction is a RATIO of the same drug in two tissues:")
    a("")
    a("```")
    a("  R = (f · E_tumour) / (f · E_normal) = E_tumour / E_normal")
    a("```")
    a("")
    a(f"A scalar per-drug factor cancels. Verified to machine precision on all "
      f"{gist['n']} matched pairs: **{gist['scalar_cancels_exactly']}**.")
    a("")
    a("| model | geomean fold error | worst | within 3x | predicted liver ratio |")
    a("|---|---|---|---|---|")
    for label, k in (("M1 transport-only (published)", "M1_baseline"),
                     ("M1 + scalar efflux, both sides", "scalar_both_sides"),
                     ("efflux on NORMAL only (the real hepatic biology)", "normal_only"),
                     ("efflux on TUMOUR only (unsupported)", "tumour_only")):
        g = gist[k]
        a(f"| {label} | **{g['geomean_fold_error']:.2f}x** | "
          f"{g['max_fold_error']:.2f}x | {g['n_within_3x']}/{gist['n']} | "
          f"{g['pred_liver']:.4f} |")
    a("")
    a("**This is the same class of result the Moncayo explant run produced, reached "
      "independently.** There, a spatially uniform ion-trapping retardation `R` was "
      "shown to give `C = R·C_bath` everywhere at equilibrium, so it could buy a time "
      "constant but never a spatial pattern. Here, a drug-uniform efflux factor `f` "
      "gives `f·E_t / f·E_n`, so it can shift a drug's absolute penetrance but never a "
      "tissue contrast. Two different mechanisms, the same algebra: **a scalar cannot "
      "produce a differential.**")
    a("")
    a("### 5.2 And the tissue-differential version has the wrong sign")
    a("")
    a("Suppose we abandoned the scalar and gave tumour and normal their own factors. "
      "Two things follow, and they are both bad.")
    a("")
    a(f"- **The real hepatic biology moves the prediction AWAY from the measurement.** "
      f"P-gp and BCRP in liver sit on the CANALICULAR membrane and export imatinib into "
      f"BILE — i.e. they REMOVE drug from normal hepatocytes. Applying efflux to the "
      f"normal side only takes the liver prediction from "
      f"{gist['M1_baseline']['pred_liver']:.3f} to "
      f"{gist['normal_only']['pred_liver']:.3f} against a measured geometric mean of "
      f"{gist['measured_liver_geomean']:.3f}, and the fold error from "
      f"{gist['M1_baseline']['geomean_fold_error']:.2f}x to "
      f"{gist['normal_only']['geomean_fold_error']:.2f}x. Worse, and worse for the "
      f"right reason: efflux on the normal side is the opposite of the *uptake* the "
      f"GIST run named.")
    a(f"- **The version that WOULD help is the one the biology contradicts.** Efflux on "
      f"the tumour side only gives {gist['tumour_only']['pred_liver']:.4f} and a fold "
      f"error of {gist['tumour_only']['geomean_fold_error']:.2f}x — numerically better "
      f"than the published {gist['M1_baseline']['geomean_fold_error']:.2f}x. But it "
      f"asserts that a GIST liver metastasis effluxes imatinib while the adjacent "
      f"normal liver does not, and normal liver is among the highest "
      f"ABCB1/ABCG2-expressing tissues in the body. The required tumour/normal factor "
      f"ratio to hit the measurement exactly is "
      f"**{gist['required_tumour_over_normal_factor']:.4f}**, i.e. an efflux-ratio "
      f"asymmetry of **{gist['required_er_asymmetry']:.1f}x in favour of the tumour** — "
      f"in the direction opposite to known transporter expression. Buying a better "
      f"number by asserting inverted biology is not progress, and we are not doing it.")
    a("")
    a("**Conclusion for GIST: unchanged.** The transport model remains structurally "
      "falsified for tissue partition, and efflux is confirmed NOT to be the missing "
      "mechanism. The GIST report's own diagnosis stands: what is missing is active "
      "hepatic UPTAKE (OCT1/OATP), which is a fundamentally different object — it "
      "raises normal-tissue concentration rather than lowering it, and it is "
      "tissue-specific rather than drug-scalar.")
    a("")

    # ---------------- 6. PARP external check ----------------
    a("## 6. An external check the term also fails")
    a("")
    a("Neither scorecard contains a PARP inhibitor, and the following contrast played "
      "no part in building the table — so this is the one genuinely out-of-sample test "
      "available.")
    a("")
    a(f"Measured ({parp['citation']}): niraparib brain:plasma AUC "
      f"{parp['measured']['niraparib']:.2f} vs olaparib "
      f"{parp['measured']['olaparib']:.2f} — a "
      f"**{parp['measured_niraparib_over_olaparib']:.1f}x** separation, and one of the "
      "better-known facts about this drug class.")
    a("")
    a("| drug | in vivo KO fold | efflux factor | measured brain:plasma |")
    a("|---|---|---|---|")
    for d in ("niraparib", "olaparib", "rucaparib"):
        rec = EFFLUX_TABLE[d]
        m = parp["measured"][d]
        a(f"| {d} | {rec.fold_ko:.1f} | {parp['factors'][d]:.4f} | "
          f"{_fmt(m, '%.2f')} |")
    a("")
    a(f"- Ordering correct: **{parp['ordering_correct']}** (niraparib > olaparib).")
    a(f"- Predicted separation: **{parp['predicted_niraparib_over_olaparib']:.2f}x** "
      f"against a measured {parp['measured_niraparib_over_olaparib']:.1f}x — the term "
      f"UNDER-separates by **{parp['separation_shortfall_fold']:.1f}x**.")
    a("")
    a("So on the one external contrast available, efflux gets the direction right and "
      "the magnitude wrong by an order of magnitude — the same pattern as the Moncayo "
      "run's compound-specificity test (`ordering yes, magnitude no`). A mechanism that "
      "reproduces orderings but not magnitudes cannot be used to produce quantitative "
      "predictions, which is what `delivery_source` does.")
    a("")

    # ---------------- 7. sensitivity ----------------
    a("## 7. Does the verdict depend on anything we chose?")
    a("")
    a("Every discretionary choice in the term, crossed: functional form, which "
      "published observable wins when a drug has both, and where the factor multiplies.")
    a("")
    a("| form | prefer | placement | GM \\|fold\\| before | GM \\|fold\\| after | 3x | 10x | worst after | improved? |")
    a("|---|---|---|---|---|---|---|---|---|")
    for r in sens:
        a(f"| {r['form']} | {r['prefer']} | {r['placement']} | "
          f"{r['gm_abs_before']:.3f}x | {r['gm_abs_after']:.3f}x | "
          f"{r['within_3x_before']}->{r['within_3x_after']} | "
          f"{r['within_10x_before']}->{r['within_10x_after']} | "
          f"{r['worst_after']:.2f}x | {'YES' if r['improved'] else 'no'} |")
    a("")
    n_improved = sum(r["improved"] for r in sens)
    a(f"**{n_improved} of {len(sens)} configurations improve the scorecard.** The "
      "verdict is not an artefact of the form, the source preference, or the placement. "
      "The most permissive configuration available (`monolayer` form, which halves the "
      "penalty at high ER) still does not recover the shipped model's accuracy.")
    a("")
    a("Two rows of this table are a free correctness check rather than a sensitivity: "
      "the two `in_vivo` rows are IDENTICAL across the two forms. That is exactly what "
      "the derivation in §3 requires — with `prefer=in_vivo` every moving drug resolves "
      "through a knockout fold change, and `f = 1/fold_KO` in both forms — so if those "
      "rows had differed, the two forms would not be two windows onto one model.")
    a("")

    # ---------------- 8. verdict ----------------
    a("## 8. Verdict and recommendation")
    a("")
    a("### Recommendation: **default OFF**, and the module ships that way.")
    a("")
    a("`treatment.drugs.delivery_source(..., efflux=False)` is the default and is "
      "bit-identical to the previous code (§1). `molecular.efflux` is fully "
      "implemented, fully cited, self-tested, and available to any caller who passes "
      "`efflux=True` — as a documented, available mechanism, not as an endorsed one.")
    a("")
    a("### What this run establishes")
    a("")
    a("1. **Efflux is implementable from the literature with zero free parameters.** "
      f"{s['n_drugs']} drugs assigned, every claim cited, {s['substrate_with_number']} "
      "carrying a usable published ratio. No value was fitted to any dataset in this "
      "repo.")
    a("2. **It does not improve the one measured scorecard it can act on**, in any of "
      f"the {len(sens)} configurations tested: GM absolute fold error "
      f"{b['gm_abs']:.2f}x -> {af['gm_abs']:.2f}x, order-of-magnitude band "
      f"{b['n_within_10x']}/6 -> {af['n_within_10x']}/6, worst case "
      f"{b['worst_abs_fold']:.1f}x -> {af['worst_abs_fold']:.1f}x.")
    a("3. **It cannot act on the other scorecard at all**, and the proof is algebraic "
      "rather than empirical: a scalar per-drug factor cancels out of a within-drug "
      "tissue ratio. This independently reproduces the structural lesson of the Moncayo "
      "explant run.")
    a("4. **The failure is diagnostic, and it names the same mechanism GIST named.** "
      "Adavosertib fails because human glioma UPTAKE (OATP1A2) beats efflux in the "
      "measured tissue. GIST fails because normal-liver UPTAKE (OCT1/OATP) is the "
      "missing term and efflux has the wrong sign there. Two independent datasets, two "
      "different tissues, the same answer: **the missing mechanism is uptake transport, "
      "not efflux.** That is a sharper statement than either report could make alone, "
      "and it is the useful output of this run.")
    a("")
    a("### What this run does NOT establish")
    a("")
    a("- **Not that efflux is unimportant in vivo.** It plainly is — a 30x mouse brain "
      "knockout effect for eribulin is real. What is shown is that a scalar per-drug "
      "efflux term does not improve THESE six measured human partition values, and "
      "cannot touch the GIST target at all.")
    a("- **Nothing at adequate power.** Six drugs, three of which move; two GIST "
      "placements that differ by an algebraic identity rather than by data. No p-value "
      "here should be read as evidence for anything.")
    a("- **Nothing about the drugs with no number.** Half the confirmed substrates in "
      "the table get a factor of 1.0 because no usable ratio is published, so this test "
      "is silent about them — including doxorubicin, the classic P-gp substrate.")
    a("- **Nothing about efflux inhibition, resistance, or intracellular accumulation.** "
      "This is the delivery axis only.")
    a("")
    a("### What would change the verdict")
    a("")
    a("A **per-tissue transporter field** — the object that both this run and the GIST "
      "run converge on needing, because it is what would make the term non-scalar and "
      "therefore able to produce a tissue contrast at all.")
    a("")
    a("This is better placed than the spatial acidic-compartment field the Moncayo run "
      "needed, for which no measurement in usable units exists at all. Absolute BBB "
      "transporter protein abundance HAS been quantified by targeted proteomics (PMID "
      "21291474: human BBB P-gp 6.06 fmol/ug, BCRP 8.14 fmol/ug), which is a real "
      "parameterisation in real units — the thing ion-trapping lacked. But two specific "
      "obstacles are worth stating so nobody starts this expecting it to be easy:")
    a("")
    a("1. **Not from the spatial data currently in this repo.** `reconstruct.real_xenium` "
      "documents its own panel audit, and ABCB1 and ABCG2 are explicitly on the ABSENT "
      "list for the 319-gene Xenium Human Brain panel. The Visium path carries "
      "deconvolved cell-type fractions, not a transporter transcript field. So a "
      "spatially-resolved transporter map would need a dataset this repo does not "
      "currently load.")
    a("2. **Transcript is not transport.** Even with expression in hand, converting "
      "mRNA to a transport clearance requires exactly the kind of scaling factor "
      "refused in §3 — which is fittable only against the same six drugs, and therefore "
      "circular for the same reason.")
    a("")
    a("**And efflux is the lower-value of the two candidates.** Both measured datasets "
      "point at UPTAKE, not efflux: OATP1A2 in human glioma (which is why adavosertib "
      "accumulates where this term predicts it should be excluded) and OCT1/OATP in "
      "normal liver (which is what GIST named, and which efflux has the wrong sign for). "
      "An uptake term is the better next piece of work, and this run's main contribution "
      "is having ruled out the alternative on measured data rather than on intuition.")
    a("")
    a("---")
    a("")
    a("*Generated by `tmesim/validation/efflux_impact.py`. Additive only: `core`, "
      "`physics`, `sim.engine`, `molecular.permeability` and `reconstruct.barrier` are "
      "untouched; `treatment.drugs.delivery_source` gained one default-off keyword whose "
      "bit-exactness when off is verified in §1.*")

    text = "\n".join(A) + "\n"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return text


# ===========================================================================
# 7. Entry point
# ===========================================================================
def run(write_md: bool = True) -> Dict:
    """Run every check and (optionally) write ``results/efflux_impact.md``."""
    compat = backward_compatibility_check()
    subs = substrate_rows()
    sc = glioma_scorecard()
    paired = glioma_paired_test(sc)
    loo = per_drug_leave_one_out(sc)
    gist = gist_impact()
    parp = parp_ordering_check()
    sens = sensitivity()
    if write_md:
        write_report(compat, subs, sc, paired, loo, gist, parp, sens)
    return {"backward_compatibility": compat, "substrate_rows": subs,
            "glioma": sc, "paired": paired, "leave_one_out": loo,
            "gist": gist, "parp": parp, "sensitivity": sens,
            "summary": substrate_summary()}


if __name__ == "__main__":
    res = run()
    c, sc, g, p = res["backward_compatibility"], res["glioma"], res["gist"], res["parp"]
    b, af = sc["before"], sc["after"]
    print(f"backward compatible (bit-exact when off): {c['all_bit_identical']} "
          f"({c['n_bit_identical']}/{c['n_comparisons']})")
    print()
    print("GLIOMA SCORECARD (n = 6 measured human drugs)")
    print(f"  GM |fold error|  {b['gm_abs']:.3f}x -> {af['gm_abs']:.3f}x   "
          f"({'IMPROVED' if sc['improved'] else 'WORSE'})")
    print(f"  within 3x        {b['n_within_3x']}/6 -> {af['n_within_3x']}/6")
    print(f"  within 10x       {b['n_within_10x']}/6 -> {af['n_within_10x']}/6")
    print(f"  worst            {b['worst_abs_fold']:.2f}x -> {af['worst_abs_fold']:.2f}x")
    for r in sc["rows"]:
        print(f"    {r['drug']:14s} f={r['efflux_factor']:.4f}  "
              f"{r['fold_before']:8.3f}x -> {r['fold_error']:8.3f}x   "
              f"({r['efflux_reason']})")
    print()
    print("GIST TUMOUR:NORMAL (n = 6 matched pairs)")
    print(f"  scalar efflux cancels exactly: {g['scalar_cancels_exactly']}")
    print(f"  M1 published            {g['M1_baseline']['geomean_fold_error']:.2f}x")
    print(f"  + scalar efflux         {g['scalar_both_sides']['geomean_fold_error']:.2f}x")
    print(f"  efflux on normal only   {g['normal_only']['geomean_fold_error']:.2f}x "
          f"(the real hepatic biology -- WORSE)")
    print(f"  efflux on tumour only   {g['tumour_only']['geomean_fold_error']:.2f}x "
          f"(requires inverted biology)")
    print()
    print("PARP external check")
    print(f"  measured niraparib/olaparib  {p['measured_niraparib_over_olaparib']:.1f}x")
    print(f"  predicted                    {p['predicted_niraparib_over_olaparib']:.2f}x "
          f"(under-separates {p['separation_shortfall_fold']:.1f}x)")
    print()
    n_imp = sum(r["improved"] for r in res["sensitivity"])
    print(f"sensitivity: {n_imp}/{len(res['sensitivity'])} configurations improve "
          f"the scorecard")
    print(f"\nRECOMMENDATION: default OFF.  Wrote {REPORT_MD}")
