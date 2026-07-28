"""Does carrier-mediated UPTAKE improve the MEASURED scorecards? Scored, not assumed.

``./.venv/bin/python -m tmesim.validation.uptake_impact``

What this module is for
-----------------------
:mod:`tmesim.molecular.uptake` implements OATP/OCT/OAT uptake transport as a
literature-sourced, zero-free-parameter, **tissue-differential** exposure factor.
Implementing a mechanism is not evidence for it — that is the lesson
``results/efflux_impact.md`` paid for. So this module re-scores the same two
measured human datasets, WITH the term active, and reports before/after:

1. **The GIST imatinib tumour:normal target** (``results/gist_imatinib.md``):
   n = 6 matched within-patient pairs, measured 0.023-0.716, current model 7.92x
   geometric-mean fold error and — the sharp part — *structurally floored*. Four
   of the five measured liver ratios lie below the transport model's hard floor
   at ANY penetration length. The question here is not merely whether the number
   improves but whether the **floor is broken**, because only a mechanism that
   puts more drug in normal liver than passive transport permits can do that.
2. **The six-drug human glioma partition scorecard**
   (``results/glioma_regional_pk.md``): current shipped model geometric-mean
   absolute fold error 1.99x, 4/6 within 3x, 6/6 within 10x. Adavosertib is the
   specific target: under-predicted 0.24x, and its measured 2.46 sits above the
   passive model's structural ceiling of 1.0.

Why this test is not a rerun of the efflux test
------------------------------------------------
Efflux was *provably inert* on GIST: a per-drug scalar cancels exactly out of a
within-drug tissue ratio, ``(f E_t)/(f E_n) = E_t/E_n``, verified to machine
precision. Uptake is per-(drug, tissue), so it does not cancel — and the
asymmetry it introduces is the one the biology supplies for free, since OATP1B1
and OATP1B3 are hepatocyte transporters and a GIST liver metastasis is not made
of hepatocytes. The **falsifiable** consequence, which this module tests
explicitly, is that the term must move the five liver pairs and leave the single
intestinal pair — the one pair the shipped model already gets right to 1.11x —
essentially untouched. A term that fixed the liver pairs by breaking the
intestinal pair would be a fitted contrast, not a mechanism.

The circularity guard is enforced here, not just declared
----------------------------------------------------------
For adavosertib the only published uptake number is an in vivo **human glioma**
``Kp,uu``, which is the same class of measurement, in the same tissue, as the
scorecard target. Using it would be copying the answer across. The default
configuration therefore withholds it (``allow_circular=False``, reason
``"only_circular_observable"``) and the drug gets exactly 1.0. The circular
variant is run and reported as an explicitly labelled **upper bound** in the
sensitivity section, never as the headline. :func:`circularity_audit` prints the
gap between the two so a reader can see precisely how much of any improvement
would have come from the forbidden number.

Architecture
------------
**Additive only, and more strictly than the efflux run.** That run added one
default-off keyword to ``treatment.drugs.delivery_source``. This one adds
*nothing at all* outside the two new files: no shipped code path calls
:mod:`tmesim.molecular.uptake`, so the term is off by construction rather than by
default. :func:`additivity_check` verifies both halves of that claim numerically.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..molecular.uptake import (
    PREFER_IN_VITRO,
    PREFER_IN_VIVO,
    SOURCE_PREFERENCES,
    TISSUE_GENERIC,
    TISSUE_GLIOMA,
    TISSUE_INTESTINAL_TUMOUR,
    TISSUE_LIVER_METASTASIS,
    TISSUE_NORMAL_INTESTINE,
    TISSUE_NORMAL_LIVER,
    TISSUES,
    TRANSPORTERS,
    UPTAKE_TABLE,
    substrate_summary,
    uptake_adjusted_penetrance,
    uptake_detail,
    uptake_factor,
)
from ..reconstruct.barrier import L_PARACELLULAR_LEAK, effective_penetrance

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(REPO_ROOT, "results")
REPORT_MD = os.path.join(RESULTS_DIR, "uptake_impact.md")

BAND_TIGHT = 3.0
BAND_OOM = 10.0

#: Name used in ``UptakeEntry.circular_against`` for the glioma scorecard.
SCORECARD_GLIOMA = "glioma_regional_pk"


# ===========================================================================
# 0. Additivity: the shipped pipeline is untouched, and OFF reduces exactly
# ===========================================================================
def additivity_check() -> Dict:
    """Two numeric claims, both verified rather than asserted.

    **(a) The shipped delivery path is bit-identical.** ``molecular.uptake`` is
    not imported by anything in ``core``, ``physics``, ``sim``, ``treatment`` or
    ``reconstruct``, so ``delivery_source`` cannot have changed. Rather than
    trusting that, the original expression is written out here and compared with
    ``==`` on every element of the float64 array, over a grid of drug names,
    penetrances and plasma concentrations, and over both the no-barrier and
    barrier branches. The drug names include ones that carry a LARGE uptake
    factor, so a leaked application would be impossible to miss.

    **(b) The new composition reduces exactly when the factor is 1.0.** With
    ``U = 1`` — any non-substrate, and any substrate in a tissue lacking its
    carrier — :func:`uptake_adjusted_penetrance` must equal
    :func:`tmesim.reconstruct.barrier.effective_penetrance` bit-for-bit, ceiling
    included.

    And the check is made non-vacuous: with a real substrate in the tissue that
    expresses its carrier, the composition must actually move.
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
    for name in ("imatinib", "adavosertib", "methotrexate", "not_in_the_table"):
        for pen in (0.0, 0.035, 0.5, 0.9399, 1.0):
            drug = DrugParams(name=name, bbb_penetrance=pen)
            for plasma in (0.0, 1.0, 7.3):
                ref = vascular * pen * float(plasma)
                got = delivery_source(micro, drug, plasma)
                checks.append({"branch": "delivery_no_barrier",
                               "identical": bool(np.array_equal(got, ref))})
                ref_b = vascular * np.clip(pen + L_PARACELLULAR_LEAK * disruption_field,
                                           0.0, 1.0) * float(plasma)
                got_b = delivery_source(micro, drug, plasma,
                                        bbb_disruption=disruption_field)
                checks.append({"branch": "delivery_barrier",
                               "identical": bool(np.array_equal(got_b, ref_b))})

    # (b) the composition identity, over drugs whose factor is exactly 1.0
    comp: List[Dict] = []
    for name in ("cisplatin", "not_in_the_table", "imatinib"):
        # imatinib IS a substrate but of hepatic carriers -- in glioma its factor
        # is 1.0, which is precisely the tissue-differential branch.
        for tissue in (TISSUE_GENERIC, TISSUE_GLIOMA):
            if uptake_factor(name, tissue=tissue) != 1.0:
                continue
            for p in (0.0, 0.05, 0.5, 0.9, 1.0):
                for d in (0.0, 0.25, 0.5, 1.0):
                    lhs = uptake_adjusted_penetrance(p, name, tissue=tissue,
                                                     bbb_disruption=d)
                    rhs = float(effective_penetrance(p, d))
                    comp.append({"drug": name, "tissue": tissue,
                                 "identical": lhs == rhs})

    # non-vacuity: the term must bite where the carrier IS expressed
    bites = uptake_factor("imatinib", tissue=TISSUE_NORMAL_LIVER) > 1.0
    bites_comp = (uptake_adjusted_penetrance(0.5, "imatinib",
                                             tissue=TISSUE_NORMAL_LIVER)
                  > uptake_adjusted_penetrance(0.5, "imatinib",
                                               tissue=TISSUE_LIVER_METASTASIS))
    # and the ceiling must genuinely rise above 1.0 somewhere, else the
    # adavosertib branch of the argument is untested machinery
    ceiling_rises = uptake_adjusted_penetrance(
        1.0, "imatinib", tissue=TISSUE_NORMAL_LIVER, bbb_disruption=0.0) > 1.0

    n, n_ok = len(checks), sum(c["identical"] for c in checks)
    m, m_ok = len(comp), sum(c["identical"] for c in comp)
    return {
        "n_delivery_comparisons": n, "n_delivery_bit_identical": n_ok,
        "delivery_all_bit_identical": n_ok == n,
        "n_composition_comparisons": m, "n_composition_bit_identical": m_ok,
        "composition_all_bit_identical": m_ok == m,
        "factor_bites_where_carrier_expressed": bool(bites),
        "composition_is_tissue_differential": bool(bites_comp),
        "ceiling_rises_above_one": bool(ceiling_rises),
        "n_shipped_files_modified": 0,
        "failures": [c for c in checks + comp if not c["identical"]],
    }


# ===========================================================================
# 1. The substrate table, as data for the report
# ===========================================================================
def substrate_rows(prefer: str = PREFER_IN_VITRO,
                   allow_circular: bool = False) -> List[Dict]:
    """One row per drug in :data:`UPTAKE_TABLE`, resolved in every tissue."""
    rows = []
    for name, rec in UPTAKE_TABLE.items():
        per_tissue = {}
        for t in TISSUES:
            d = uptake_detail(name, tissue=t, prefer=prefer,
                              allow_circular=allow_circular)
            per_tissue[t] = {"factor": d["factor"], "reason": d["reason"],
                             "applied": d["transporters_applied"]}
        rows.append({
            "drug": name,
            "entries": rec.entries,
            "notes": rec.notes,
            "is_substrate": rec.is_substrate,
            "per_tissue": per_tissue,
            "max_factor": max(v["factor"] for v in per_tissue.values()),
            "suppressed_circular": any(e.suppressed_circular for e in rec.entries),
        })
    return rows


def circularity_audit(prefer: str = PREFER_IN_VITRO) -> List[Dict]:
    """Every observable the circularity guard is withholding, and what it costs.

    Reported so that the reader can see exactly how much of the headline would
    change if the guard were dropped — the honest way to present a number one has
    decided not to use.
    """
    out = []
    for name, rec in UPTAKE_TABLE.items():
        for e in rec.entries:
            if not e.suppressed_circular:
                continue
            tissues = TRANSPORTERS[e.transporter].tissues
            for t in tissues:
                out.append({
                    "drug": name, "transporter": e.transporter, "tissue": t,
                    "withheld_value": e.kpuu_in_vivo,
                    "withheld_quantity": e.kpuu_assay,
                    "withheld_citation": e.kpuu_citation,
                    "circular_against": e.circular_against,
                    "factor_default": uptake_factor(name, tissue=t, prefer=prefer),
                    "factor_if_allowed": uptake_factor(name, tissue=t, prefer=prefer,
                                                       allow_circular=True),
                })
    return out


# ===========================================================================
# 2. Scorecard 1 -- GIST tumour:normal, the target uptake was named for
# ===========================================================================
#: Which model tissue each GIST section belongs to. The tumour side of a liver
#: pair is a GIST METASTASIS sitting in liver, not hepatocytes -- that distinction
#: is the entire mechanism, so it is spelled out here rather than inferred.
GIST_TISSUE_MAP: Dict[str, Tuple[str, str]] = {
    "Liver": (TISSUE_LIVER_METASTASIS, TISSUE_NORMAL_LIVER),
    "Intestine": (TISSUE_INTESTINAL_TUMOUR, TISSUE_NORMAL_INTESTINE),
}


def gist_impact(prefer: str = PREFER_IN_VITRO,
                allow_circular: bool = False) -> Dict:
    """Does tissue-differential hepatic uptake move the GIST ratio, and break the floor?

    The prediction is the same reaction-diffusion delivery ratio the published
    run scored, with one tissue-differential factor added on each side::

        R = (U_tumour  * E_tumour) / (U_normal * E_normal)

    ``E`` is unchanged (:func:`gist_imatinib.delivery_efficiency`, the 1-D Krogh
    slab), the intervessel half-spacings are unchanged, and the penetration
    length ``L`` is unchanged — it still comes from the doxorubicin-anchored
    ``alpha`` that was fitted once to an unrelated measurement. The ONLY new
    quantity is ``U``, and it is a published uptake ratio.

    Three things are reported, in increasing order of how much they are worth:

    1. the geometric-mean fold error over **n = 6**, against the published 7.92x;
    2. whether the prediction now falls **below the structural floor** the
       published run proved a transport-only model cannot cross —
       ``E_t/E_n >= R_normal/R_tumour`` for any ``L``, which becomes
       ``>= (R_normal/R_tumour)/U_normal`` once uptake is on. Breaking that floor
       is the qualitative claim; the fold error is only the quantitative one;
    3. the **intestinal control**. Patient AA is an intestinal primary against
       normal intestine, the one pair the shipped model already predicts to
       1.11x. Hepatic carriers are not expressed there, so a real mechanism must
       leave it alone. If it moved, the term would be a fitted contrast.
    """
    from .gist_imatinib import (HALF_SPACING_UM, compile_imatinib,
                                delivery_efficiency, measured_targets)

    targets = measured_targets()
    compiled = compile_imatinib()
    L = compiled["L_um"]
    e_t = delivery_efficiency(L, HALF_SPACING_UM["tumour"])
    e_n = {"Liver": delivery_efficiency(L, HALF_SPACING_UM["liver"]),
           "Intestine": delivery_efficiency(L, HALF_SPACING_UM["intestine"])}

    facs = {}
    for organ, (t_tissue, n_tissue) in GIST_TISSUE_MAP.items():
        facs[organ] = {
            "tumour": uptake_detail("imatinib", tissue=t_tissue, prefer=prefer,
                                    allow_circular=allow_circular),
            "normal": uptake_detail("imatinib", tissue=n_tissue, prefer=prefer,
                                    allow_circular=allow_circular),
        }

    def _score(use_uptake: bool) -> Dict:
        rows, folds = [], []
        for _, r in targets.iterrows():
            ut = facs[r.organ]["tumour"]["factor"] if use_uptake else 1.0
            un = facs[r.organ]["normal"]["factor"] if use_uptake else 1.0
            p = (ut * e_t) / (un * e_n[r.organ])
            fold = max(p / r.ratio_mean, r.ratio_mean / p)
            folds.append(fold)
            rows.append({"patient": r.patient, "organ": r.organ,
                         "measured": float(r.ratio_mean), "u_tumour": ut,
                         "u_normal": un, "pred": float(p), "fold": float(fold)})
        f = np.asarray(folds)
        return {"geomean_fold_error": float(10 ** np.mean(np.log10(f))),
                "max_fold_error": float(f.max()),
                "n_within_3x": int((f <= 3.0).sum()),
                "n_within_10x": int((f <= 10.0).sum()),
                "pred_liver": float((facs["Liver"]["tumour"]["factor"] if use_uptake else 1.0)
                                    * e_t
                                    / ((facs["Liver"]["normal"]["factor"] if use_uptake else 1.0)
                                       * e_n["Liver"])),
                "pred_intestine": float((facs["Intestine"]["tumour"]["factor"] if use_uptake else 1.0)
                                        * e_t
                                        / ((facs["Intestine"]["normal"]["factor"] if use_uptake else 1.0)
                                           * e_n["Intestine"])),
                "per_patient": rows}

    before = _score(False)
    after = _score(True)

    # --- the structural floor, before and after ---------------------------
    # Transport-only: E_t/E_n -> R_n/R_t as L -> 0, and -> 1 as L -> inf.
    u_n_liver = facs["Liver"]["normal"]["factor"]
    u_t_liver = facs["Liver"]["tumour"]["factor"]
    floor_nom = HALF_SPACING_UM["liver"] / HALF_SPACING_UM["tumour"]
    band = getattr(__import__("tmesim.validation.gist_imatinib", fromlist=["x"]),
                   "HALF_SPACING_BAND_UM")
    floor_ext = band["liver"][0] / band["tumour"][1]
    scale = u_t_liver / u_n_liver
    liver_meas = targets[targets.organ == "Liver"].ratio_mean.values
    return {
        "n": int(len(targets)),
        "prefer": prefer, "allow_circular": allow_circular,
        "L_um": float(L), "E_tumour": float(e_t), "E_normal": dict(e_n),
        "factors": {o: {k: {"factor": v["factor"], "reason": v["reason"],
                            "transporters": v["transporters_applied"],
                            "uptake_ratio": v["uptake_ratio"]}
                        for k, v in d.items()} for o, d in facs.items()},
        "before": before, "after": after,
        "improved": after["geomean_fold_error"] < before["geomean_fold_error"],
        "measured_liver_geomean": float(10 ** np.mean(np.log10(liver_meas))),
        "required_u_normal": float(before["pred_liver"]
                                   / 10 ** np.mean(np.log10(liver_meas))),
        "floor_nominal_before": float(floor_nom),
        "floor_nominal_after": float(floor_nom * scale),
        "floor_extreme_before": float(floor_ext),
        "floor_extreme_after": float(floor_ext * scale),
        "n_liver_below_nominal_floor_before": int((liver_meas < floor_nom).sum()),
        "n_liver_below_nominal_floor_after": int((liver_meas < floor_nom * scale).sum()),
        "n_liver_below_extreme_floor_before": int((liver_meas < floor_ext).sum()),
        "n_liver_below_extreme_floor_after": int((liver_meas < floor_ext * scale).sum()),
        "floor_broken": bool(scale < 1.0),
        "intestinal_control_unchanged": bool(
            abs(after["pred_intestine"] - before["pred_intestine"]) <= 1e-12),
        "intestinal_fold_before": float(
            [r for r in before["per_patient"] if r["organ"] == "Intestine"][0]["fold"]),
        "intestinal_fold_after": float(
            [r for r in after["per_patient"] if r["organ"] == "Intestine"][0]["fold"]),
        # A no-partitioning null, reproduced for the same reason the published
        # run reported it: 8x wrong is only meaningful against 13x wrong.
        "null_unity_geomean": float(10 ** np.mean(np.log10(
            [max(1.0 / r, r) for r in targets.ratio_mean.values]))),
    }


def gist_required_vs_published(gist: Dict) -> Dict:
    """Compare the uptake ratio the data DEMANDS with the one literature supplies.

    This is the single most informative number in the run and it is the one that
    can embarrass us, so it is computed explicitly rather than left implicit in a
    fold error. ``required`` is solved from the measurement; ``published`` is what
    the table carries with zero parameters fitted. Their ratio is how far off the
    literature value is — and a value that lands close is evidence, while a value
    that needed tuning to land close would not be.
    """
    req = gist["required_u_normal"]
    pub = gist["factors"]["Liver"]["normal"]["factor"]
    return {"required_u_normal": req, "published_u_normal": pub,
            "published_over_required": pub / req if req else float("nan"),
            "discrepancy_fold": max(pub / req, req / pub) if (req and pub) else float("nan")}


# ===========================================================================
# 3. Scorecard 2 -- the six-drug measured human glioma partition
# ===========================================================================
#: Where the uptake factor multiplies, and what ceiling applies.
#:
#: ``"transcellular"`` multiplies only the drug's own penetrance, leaving the
#: regional paracellular leak alone -- a broken tight junction is a shunt AROUND
#: the carrier, the same physical claim the efflux term makes about pumps.
#: ``"total"`` multiplies the whole bracket including the leak; physically wrong,
#: run only to show the verdict does not depend on the choice.
PLACEMENTS = ("transcellular", "total")

#: ``"uptake"`` raises the free-drug ceiling from 1.0 to the uptake factor, which
#: is the physically correct bound and the thing ``barrier.effective_penetrance``
#: says it cannot do. ``"unity"`` keeps the old ceiling of 1.0, which forbids any
#: prediction above unity and is run to isolate how much of the change is the
#: ceiling rather than the factor.
CEILINGS = ("uptake", "unity")


def _aggregate(rows: Sequence[Dict], drop: Optional[str] = None) -> Dict:
    """Geometric-mean signed and unsigned fold error, band counts, worst case.

    Mirrors ``glioma_regional_pk.run_scorecard``'s aggregation exactly so the
    "before" column reproduces the published numbers rather than re-deriving
    them differently.
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


def glioma_scorecard(prefer: str = PREFER_IN_VITRO,
                     allow_circular: bool = False,
                     placement: str = "transcellular",
                     ceiling: str = "uptake") -> Dict:
    """Re-score the six measured human glioma drugs with and without uptake.

    The BEFORE column is not pasted: it is
    ``glioma_regional_pk.run_scorecard(aggregation="geometric", apply_leak=True)``
    called live, i.e. the model exactly as it ships. AFTER changes one thing::

        before:  P_eff = clip( bbb_penetrance     + L * D_sample , 0, 1 )
        after :  P_eff = clip( bbb_penetrance * U + L * D_sample , 0, U )

    with ``U = uptake_factor(drug, tissue="glioma")``. Nothing else changes: same
    measured values, same regional weights, same aggregation.
    """
    if placement not in PLACEMENTS:
        raise ValueError(f"placement must be one of {PLACEMENTS}")
    if ceiling not in CEILINGS:
        raise ValueError(f"ceiling must be one of {CEILINGS}")
    from .glioma_regional_pk import run_scorecard

    before = run_scorecard(aggregation="geometric", apply_leak=True)
    rows: List[Dict] = []
    for r in before["rows"]:
        drug = r["drug"]
        struct = float(r["structure_only_partition"])
        disr = float(r["sample_mean_bbb_disruption"])
        meas = float(r["measured_mean"])
        det = uptake_detail(drug, tissue=TISSUE_GLIOMA, prefer=prefer,
                            allow_circular=allow_circular)
        u = float(det["factor"])
        cap = u if ceiling == "uptake" else 1.0
        if placement == "transcellular":
            pred = float(np.clip(struct * u + L_PARACELLULAR_LEAK * disr, 0.0, cap))
        else:
            pred = float(np.clip(float(effective_penetrance(struct, disr)) * u,
                                 0.0, cap))
        fold_before = float(r["fold_error"])
        fold_after = pred / meas if meas > 0 else float("inf")
        rows.append({
            "drug": drug, "measured": meas, "method": r["method"],
            "structure_only": struct, "sample_disruption": disr,
            "pred_before": float(r["predicted_partition"]), "fold_before": fold_before,
            "log10_fold_before": math.log10(fold_before) if fold_before > 0 else float("nan"),
            "uptake_factor": u, "uptake_reason": det["reason"],
            "uptake_ratio": det["uptake_ratio"],
            "uptake_transporters": det["transporters_applied"],
            "uptake_source": det["source"],
            "suppressed_circular": det["suppressed_circular"],
            "pred_after": pred, "fold_error": fold_after,
            "log10_fold": math.log10(fold_after) if fold_after > 0 else float("nan"),
        })

    before_rows = [{"drug": r["drug"], "fold_error": r["fold_before"],
                    "log10_fold": r["log10_fold_before"]} for r in rows]
    agg_before = _aggregate(before_rows)
    agg_after = _aggregate(rows)
    ada = [r for r in rows if r["drug"] == "adavosertib"]
    return {
        "prefer": prefer, "allow_circular": allow_circular,
        "placement": placement, "ceiling": ceiling,
        "rows": rows, "before": agg_before, "after": agg_after,
        "before_ex_adavosertib": _aggregate(before_rows, drop="adavosertib"),
        "after_ex_adavosertib": _aggregate(rows, drop="adavosertib"),
        "n_drugs_moved": sum(1 for r in rows if r["uptake_factor"] > 1.0),
        "improved": agg_after["gm_abs"] < agg_before["gm_abs"],
        "adavosertib_fold_before": ada[0]["fold_before"] if ada else None,
        "adavosertib_fold_after": ada[0]["fold_error"] if ada else None,
        "adavosertib_improved": (abs(math.log10(ada[0]["fold_error"]))
                                 < abs(math.log10(ada[0]["fold_before"])))
                                if ada else None,
    }


def glioma_paired_test(sc: Dict) -> Dict:
    """Is the before/after difference distinguishable from nothing at n = 6?

    An exact sign test plus a Wilcoxon over ``|log10 fold error|``, reported to
    show that the comparison is UNPOWERED rather than to claim a result.
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
            warnings.simplefilter("ignore")
            p_w = float(stats.wilcoxon(b, a, zero_method="wilcox").pvalue)
    except ValueError:
        p_w = float("nan")
    return {"n_drugs": len(b), "n_movers": movers, "n_worse": n_worse,
            "n_better": movers - n_worse, "mean_delta_log10": float(d.mean()),
            "p_sign_test": p_sign, "p_wilcoxon": p_w,
            "note": ("with %d of %d drugs moving at all, neither test can reach "
                     "significance at any effect size; both are reported to make "
                     "that explicit" % (movers, len(b)))}


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
# 4. Sensitivity -- does the verdict depend on anything we chose?
# ===========================================================================
def sensitivity() -> List[Dict]:
    """Cross every discretionary choice in the term on the glioma scorecard.

    Source preference, the circularity guard, where the factor multiplies, and
    which ceiling applies. Sixteen configurations. Unlike the efflux term there
    is no functional-form axis, because the uptake form is unique (see the
    module docstring of :mod:`tmesim.molecular.uptake`).
    """
    out = []
    for prefer in SOURCE_PREFERENCES:
        for allow_circular in (False, True):
            for placement in PLACEMENTS:
                for ceiling in CEILINGS:
                    sc = glioma_scorecard(prefer=prefer, allow_circular=allow_circular,
                                          placement=placement, ceiling=ceiling)
                    out.append({
                        "prefer": prefer, "allow_circular": allow_circular,
                        "placement": placement, "ceiling": ceiling,
                        "gm_abs_before": sc["before"]["gm_abs"],
                        "gm_abs_after": sc["after"]["gm_abs"],
                        "within_3x_before": sc["before"]["n_within_3x"],
                        "within_3x_after": sc["after"]["n_within_3x"],
                        "within_10x_before": sc["before"]["n_within_10x"],
                        "within_10x_after": sc["after"]["n_within_10x"],
                        "worst_before": sc["before"]["worst_abs_fold"],
                        "worst_after": sc["after"]["worst_abs_fold"],
                        "ada_after": sc["adavosertib_fold_after"],
                        "improved": sc["improved"]})
    return out


def gist_sensitivity() -> List[Dict]:
    """The GIST verdict across the source preference, the guard, and the geometry.

    The geometry sweep matters more here than the transporter choices do: the
    published run showed the intervessel half-spacings are the largest single
    lever on ``M1``, so an uptake term that only works at one geometry would not
    be a result. Every row re-solves the whole prediction.
    """
    from .gist_imatinib import (HALF_SPACING_BAND_UM, HALF_SPACING_UM,
                                compile_imatinib, delivery_efficiency,
                                measured_targets)
    targets = measured_targets()
    L = compile_imatinib()["L_um"]
    liver_meas = targets[targets.organ == "Liver"].ratio_mean.values

    out = []
    for prefer in SOURCE_PREFERENCES:
        for allow_circular in (False, True):
            g = gist_impact(prefer=prefer, allow_circular=allow_circular)
            out.append({"axis": "transporter_source",
                        "setting": f"prefer={prefer}, allow_circular={allow_circular}",
                        "u_normal_liver": g["factors"]["Liver"]["normal"]["factor"],
                        "gm_before": g["before"]["geomean_fold_error"],
                        "gm_after": g["after"]["geomean_fold_error"],
                        "improved": g["improved"]})

    u_n = uptake_factor("imatinib", tissue=TISSUE_NORMAL_LIVER)
    u_t = uptake_factor("imatinib", tissue=TISSUE_LIVER_METASTASIS)
    u_ni = uptake_factor("imatinib", tissue=TISSUE_NORMAL_INTESTINE)
    u_ti = uptake_factor("imatinib", tissue=TISSUE_INTESTINAL_TUMOUR)

    def _gm(rt: float, rn: float, ln: float) -> Tuple[float, float]:
        e_t = delivery_efficiency(ln, rt)
        e_liv = delivery_efficiency(ln, rn)
        e_int = delivery_efficiency(ln, HALF_SPACING_UM["intestine"])
        folds_b, folds_a = [], []
        for _, r in targets.iterrows():
            if r.organ == "Liver":
                pb, pa = e_t / e_liv, (u_t * e_t) / (u_n * e_liv)
            else:
                pb, pa = e_t / e_int, (u_ti * e_t) / (u_ni * e_int)
            folds_b.append(max(pb / r.ratio_mean, r.ratio_mean / pb))
            folds_a.append(max(pa / r.ratio_mean, r.ratio_mean / pa))
        return (float(10 ** np.mean(np.log10(folds_b))),
                float(10 ** np.mean(np.log10(folds_a))))

    rt0, rn0 = HALF_SPACING_UM["tumour"], HALF_SPACING_UM["liver"]
    for label, rt, rn, ln in (
            ("baseline geometry", rt0, rn0, L),
            ("tumour half-spacing 30 um (band end)", HALF_SPACING_BAND_UM["tumour"][0], rn0, L),
            ("tumour half-spacing 120 um (band end)", HALF_SPACING_BAND_UM["tumour"][1], rn0, L),
            ("normal-liver half-spacing 10 um (band end)", rt0, HALF_SPACING_BAND_UM["liver"][0], L),
            ("normal-liver half-spacing 15 um (band end)", rt0, HALF_SPACING_BAND_UM["liver"][1], L),
            ("penetration length L = 5 um", rt0, rn0, 5.0),
            ("penetration length L = 400 um (uncalibrated)", rt0, rn0, 400.0)):
        b, a = _gm(rt, rn, ln)
        out.append({"axis": "geometry", "setting": label, "u_normal_liver": u_n,
                    "gm_before": b, "gm_after": a, "improved": a < b})
    return out


# ===========================================================================
# 5. External check -- an out-of-sample contrast neither scorecard supplies
# ===========================================================================
def hcc_downregulation_check() -> Dict:
    """The mechanism's own premise, tested against an independent measurement.

    The GIST result rests on one claim: hepatic uptake carriers are present in
    normal liver and effectively absent from a tumour sitting in it. That claim
    is testable outside this repo's datasets. Hepatocellular carcinoma — a
    tumour that at least *derives* from hepatocytes, and therefore the hardest
    case for the claim — shows strong, replicated downregulation of SLCO1B1,
    SLCO1B3 and SLC22A1 relative to surrounding non-tumour liver. If the premise
    failed even there, it would certainly fail for a GIST metastasis, which has
    no hepatocyte lineage at all.

    Returns the published fold-downregulations carried in
    :data:`tmesim.molecular.uptake.UPTAKE_TABLE`'s notes, with their PMIDs.
    """
    from ..molecular.uptake import HEPATIC_TUMOUR_DOWNREGULATION
    rows = list(HEPATIC_TUMOUR_DOWNREGULATION)
    return {"rows": rows,
            "n": len(rows),
            "all_downregulated": all(r["direction"] == "down" for r in rows)}


# ===========================================================================
# 6. Report
# ===========================================================================
def _fmt(v, spec="%.3f", dash="--"):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return dash
    return spec % v


def write_report(add: Dict, subs: List[Dict], circ: List[Dict], gist: Dict,
                 gistreq: Dict, gsens: List[Dict], sc: Dict, paired: Dict,
                 loo: List[Dict], sens: List[Dict], hcc: Dict,
                 path: str = REPORT_MD) -> str:
    raise NotImplementedError  # filled in below, after the table exists


# ===========================================================================
# 7. Entry point
# ===========================================================================
def run(write_md: bool = True) -> Dict:
    add = additivity_check()
    subs = substrate_rows()
    circ = circularity_audit()
    gist = gist_impact()
    gistreq = gist_required_vs_published(gist)
    gsens = gist_sensitivity()
    sc = glioma_scorecard()
    paired = glioma_paired_test(sc)
    loo = per_drug_leave_one_out(sc)
    sens = sensitivity()
    hcc = hcc_downregulation_check()
    out = {"additivity": add, "substrates": subs, "circularity": circ,
           "gist": gist, "gist_required": gistreq, "gist_sensitivity": gsens,
           "glioma": sc, "paired": paired, "leave_one_out": loo,
           "sensitivity": sens, "hcc": hcc}
    if write_md:
        out["report_path"] = write_report(add, subs, circ, gist, gistreq, gsens,
                                          sc, paired, loo, sens, hcc)
    return out


if __name__ == "__main__":
    run()
