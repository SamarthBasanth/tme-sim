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

    # Non-vacuity. These MUST use a drug/tissue pair that actually carries a
    # number, otherwise the whole check passes trivially. Doxorubicin in glioma
    # (OATP1A2, 2.0-fold, PMID 24554572) is the strongest such pair in the table.
    bites = uptake_factor("doxorubicin", tissue=TISSUE_GLIOMA) > 1.0
    bites_comp = (uptake_adjusted_penetrance(0.5, "doxorubicin",
                                             tissue=TISSUE_GLIOMA)
                  > uptake_adjusted_penetrance(0.5, "doxorubicin",
                                               tissue=TISSUE_NORMAL_LIVER))
    # and the ceiling must genuinely rise above 1.0 somewhere, else the
    # ceiling-raising branch of the argument is untested machinery
    ceiling_rises = uptake_adjusted_penetrance(
        1.0, "doxorubicin", tissue=TISSUE_GLIOMA, bbb_disruption=0.0) > 1.0

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
    from .gist_imatinib import HALF_SPACING_BAND_UM
    floor_nom = HALF_SPACING_UM["liver"] / HALF_SPACING_UM["tumour"]
    floor_ext = HALF_SPACING_BAND_UM["liver"][0] / HALF_SPACING_BAND_UM["tumour"][1]
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


#: The largest hepatic uptake ratio for imatinib that ANY published source can
#: support, used as a falsifying UPPER BOUND and never as a value.
#:
#: Zimmerman et al. 2013 (PMID 23340295) publish imatinib's OATP1B1/OATP1B3
#: uptake only as bars in Figure 1A, with no number anywhere in the text, tables
#: or caption. What the figure DOES license is a bound: sorafenib was the highest
#: TKI in that panel at ~1.5-fold, so imatinib is necessarily below 1.5-fold. The
#: only numeric imatinib uptake ratio that exists anywhere — OATP1A2, 1.5-fold at
#: pH 7.4 (PMID 21508937) — happens to land at the same value, for a transporter
#: that is not hepatic.
#:
#: Using a bound to FALSIFY is legitimate; using it as a value would not be. It
#: is therefore never fed to :func:`uptake_factor` and appears only in
#: :func:`gist_upper_bound`.
IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND = 1.5
IMATINIB_UPPER_BOUND_CITATION = (
    "PMID 23340295 (Zimmerman 2013 Clin Cancer Res), Figure 1A: imatinib below the "
    "panel maximum of ~1.5-fold (sorafenib); concordant with the only numeric "
    "imatinib uptake ratio published for any carrier, 1.5-fold on OATP1A2 at "
    "pH 7.4, PMID 21508937 (Eechoute 2011 Clin Pharmacol Ther)")


def gist_upper_bound(gist: Dict) -> Dict:
    """Is the uptake hypothesis merely unimplemented, or quantitatively excluded?

    This is the analogue, for the uptake mechanism, of the structural-floor
    argument that falsified the transport model in ``results/gist_imatinib.md``
    — and it is the reason this run produces a *result* rather than only a
    coverage complaint.

    The table gives imatinib a hepatic uptake factor of exactly 1.0 because no
    number is published. That alone would leave open the possibility that the
    mechanism is right and merely unquantified. So we ask the sharper question:
    **what happens at the largest value any published source can support?**
    :data:`IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND` is that ceiling, and if even it
    falls far short of the value the measurement demands, then the hepatic-uptake
    hypothesis is not just unimplemented — it is *insufficient by the literature's
    own numbers*.
    """
    from .gist_imatinib import (HALF_SPACING_UM, compile_imatinib,
                                delivery_efficiency, measured_targets)
    targets = measured_targets()
    L = compile_imatinib()["L_um"]
    e_t = delivery_efficiency(L, HALF_SPACING_UM["tumour"])
    e_n = {"Liver": delivery_efficiency(L, HALF_SPACING_UM["liver"]),
           "Intestine": delivery_efficiency(L, HALF_SPACING_UM["intestine"])}
    required = gist["required_u_normal"]

    rows = []
    for u_n in (1.0, 1.2, IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND, 2.0, 3.0, 5.0,
                8.0, required, 15.0, 20.0):
        folds = []
        for _, r in targets.iterrows():
            un = u_n if r.organ == "Liver" else 1.0
            p = e_t / (un * e_n[r.organ])
            folds.append(max(p / r.ratio_mean, r.ratio_mean / p))
        f = np.asarray(folds)
        rows.append({"u_normal": float(u_n),
                     "pred_liver": float(e_t / (u_n * e_n["Liver"])),
                     "geomean_fold_error": float(10 ** np.mean(np.log10(f))),
                     "n_within_3x": int((f <= 3.0).sum()),
                     "is_upper_bound": abs(u_n - IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND) < 1e-9,
                     "is_required": abs(u_n - required) < 1e-9})
    at_bound = [r for r in rows if r["is_upper_bound"]][0]
    return {
        "upper_bound": IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND,
        "citation": IMATINIB_UPPER_BOUND_CITATION,
        "required_u_normal": required,
        "shortfall_fold": required / IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND,
        "gm_at_upper_bound": at_bound["geomean_fold_error"],
        "gm_before": gist["before"]["geomean_fold_error"],
        "pred_liver_at_upper_bound": at_bound["pred_liver"],
        "measured_liver_geomean": gist["measured_liver_geomean"],
        "floor_at_upper_bound": gist["floor_nominal_before"]
                                / IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND,
        "n_liver_below_floor_at_upper_bound": int(
            (np.asarray([0.023127, 0.141834, 0.043893, 0.043283, 0.047287])
             < gist["floor_nominal_before"] / IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND).sum()),
        "excluded_by_its_own_literature": bool(
            required / IMATINIB_HEPATIC_UPTAKE_UPPER_BOUND > 3.0),
        "sweep": rows,
    }


# ---------------------------------------------------------------------------
# The redirect: what the chain excluded too early
# ---------------------------------------------------------------------------
#: Imatinib's basic pKa, as already used by ``validation/gist_imatinib.py``'s
#: Rodgers-Rowland composition model. Not a new parameter.
IMATINIB_PKA_BASE = 8.07
PH_LYSOSOME = 4.7
PH_CYTOSOL = 7.2


def lysosomal_sizing() -> Dict:
    """Size the mechanism the chain ruled out, and show the ruling-out was invalid.

    ``results/moncayo_explant.md`` §4.3 proved that a spatially UNIFORM
    ion-trapping retardation ``R`` cannot generate spatial dispersion, because at
    equilibrium ``C = R * C_bath`` everywhere however large ``R`` becomes. That
    proof is correct and this module does not dispute it. **But it is a statement
    about dispersion WITHIN one tissue, and the GIST target is a ratio BETWEEN two
    tissues.** Two tissues with different acidic-compartment content have
    different ``R``, and

    ::

        C_tumour / C_normal = (R_t * E_t) / (R_n * E_n)

    does not cancel — exactly as a tissue-differential uptake term does not
    cancel, and for exactly the same algebraic reason. The chain that motivated
    this run therefore eliminated ion-trapping using an argument that does not
    apply to the target it was eliminated for.

    That matters because the surviving evidence points there. The group that
    dismantled the imatinib/OCT1 story concluded in the same programme that
    **lysosomal sequestration, not SLC22A1, determines intracellular imatinib**
    (Burger et al. 2015, PMID 26108972) — tested in HEK293, K562, SD-1 and, to
    the point, **GIST-T1**. The prazosin-sensitive uptake signal that the whole
    clinical "OCT-1 activity" literature rests on is, in that work, lysosomal.

    This function does not implement the mechanism. It SIZES it, using the
    Henderson-Hasselbalch trapping ratio already implemented in
    ``validation/moncayo_explant.py`` and imatinib's own pKa already used in
    ``validation/gist_imatinib.py``, and asks one question: is the lysosomal
    volume fraction that would reproduce the measurement physically plausible, or
    absurd? An absurd answer would rule the mechanism out cheaply. A plausible
    one identifies the next piece of work.
    """
    ratio = ((1.0 + 10.0 ** (IMATINIB_PKA_BASE - PH_LYSOSOME))
             / (1.0 + 10.0 ** (IMATINIB_PKA_BASE - PH_CYTOSOL)))
    gist = gist_impact()
    required_R = gist["required_u_normal"]
    # R = 1 + phi*(ratio - 1)  =>  phi = (R - 1)/(ratio - 1)
    phi_required = (required_R - 1.0) / (ratio - 1.0)
    # for reference, the value the Moncayo module uses for generic tissue
    phi_reference = 0.014
    R_at_reference = 1.0 + phi_reference * (ratio - 1.0)
    return {
        "pka_base": IMATINIB_PKA_BASE,
        "lysosome_cytosol_trapping_ratio": float(ratio),
        "required_R_normal_over_tumour": float(required_R),
        "phi_lys_required_in_hepatocytes": float(phi_required),
        "phi_lys_reference_generic_tissue": phi_reference,
        "R_at_reference_phi": float(R_at_reference),
        "moncayo_swept_max_phi": 0.05,
        "R_at_moncayo_swept_max": float(1.0 + 0.05 * (ratio - 1.0)),
        "plausible": bool(phi_required <= 0.05),
        "citation_mechanism": ("PMID 26108972 (Burger/Sparreboom 2015 Mol "
                               "Pharmacol): lysosomal sequestration, not SLC22A1, "
                               "determines intracellular imatinib; tested in "
                               "HEK293, K562, SD-1 and GIST-T1"),
        "citation_pka": ("basic pKa 8.07, the same literature value already used "
                         "by validation/gist_imatinib.py's composition model"),
        "note": ("SIZING ONLY. This is not a model, nothing here is fitted, and no "
                 "prediction is made. It asks whether the required acidic-"
                 "compartment asymmetry is physically reachable, because a "
                 "mechanism that needed an impossible value could be dismissed "
                 "without implementing it."),
    }


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
                 gistreq: Dict, gub: Dict, lys: Dict, gsens: List[Dict],
                 sc: Dict, paired: Dict, loo: List[Dict], sens: List[Dict],
                 hcc: Dict, path: str = REPORT_MD) -> str:
    A: List[str] = []
    a = A.append
    s = substrate_summary()
    b, af = sc["before"], sc["after"]
    gb, ga = gist["before"], gist["after"]

    def _dir(before: float, after: float) -> str:
        if abs(after - before) < 1e-9:
            return "IDENTICAL"
        return "BETTER" if after < before else "WORSE"

    # ---------------------------------------------------------------- header
    a("# Carrier-mediated uptake (OATP / OCT / OAT) — implemented, cited, and it "
      "does NOT help either")
    a("")
    a("**Verdict up front: this is the mechanism both prior runs named, it was "
      "implemented with zero free parameters and one unambiguous functional form, "
      "and it fails on both measured scorecards. It is left DEFAULT-OFF. But this "
      "is not the same negative twice — unlike the efflux run, this one produces a "
      "positive DISPROOF of the specific mechanism `results/gist_imatinib.md` "
      "named, a quantitative bound showing the whole hepatic-uptake class is "
      "insufficient by the literature's own numbers, and a redirect to a candidate "
      "the chain eliminated using an argument that does not apply to it.**")
    a("")
    a("| scorecard | before | after uptake | verdict |")
    a("|---|---|---|---|")
    a(f"| GIST tumour:normal, n = 6 pairs (geomean fold error) | "
      f"**{gb['geomean_fold_error']:.2f}x** | **{ga['geomean_fold_error']:.2f}x** | "
      f"{_dir(gb['geomean_fold_error'], ga['geomean_fold_error'])} — the factor "
      f"resolves to 1.0 |")
    a(f"| .. at the literature's UPPER BOUND on hepatic uptake ({gub['upper_bound']:.1f}x) "
      f"| {gb['geomean_fold_error']:.2f}x | {gub['gm_at_upper_bound']:.2f}x | "
      f"still {gub['shortfall_fold']:.1f}x short of what the data demands |")
    a(f"| human glioma partition, n = 6 drugs (GM abs fold error) | "
      f"**{b['gm_abs']:.3f}x** | **{af['gm_abs']:.3f}x** | "
      f"{_dir(b['gm_abs'], af['gm_abs'])} |")
    a(f"| .. within 3x | {b['n_within_3x']}/6 | {af['n_within_3x']}/6 | |")
    a(f"| .. within 10x | {b['n_within_10x']}/6 | {af['n_within_10x']}/6 | |")
    a(f"| .. adavosertib specifically | {sc['adavosertib_fold_before']:.3f}x | "
      f"{sc['adavosertib_fold_after']:.3f}x | unchanged — the number does not exist |")
    a("")
    a("Five findings, in order of how much they are worth:")
    a("")
    a("1. **The mechanism GIST named is not merely unquantified — it is falsified.** "
      "`results/gist_imatinib.md` named *\"active hepatic uptake (OCT1/OATP)\"*. "
      "Verified rather than assumed, the OCT1 half is contradicted by five "
      "OCT1-expressing cell lines and oocytes with a working positive control "
      "(MPP+ up 8-fold in the same oocytes), by shRNA knockdown, and — decisively "
      "— by **Oct1(-/-) mice whose HEPATIC imatinib is indistinguishable from "
      "wild-type** (PMID 24352644). Imatinib is a potent OCT1 *inhibitor* (IC50 "
      "0.095 uM), which is exactly the property that makes the "
      "inhibitor-difference assays underlying the clinical \"OCT-1 activity\" "
      "literature read as though it were a substrate.")
    a(f"2. **And the surviving hepatic-uptake candidates are quantitatively "
      f"excluded, by their own literature.** The measurement demands "
      f"`U_normal` = **{gistreq['required_u_normal']:.1f}x**. The largest value any "
      f"published source supports is **{gub['upper_bound']:.1f}x** — a "
      f"**{gub['shortfall_fold']:.1f}x shortfall**. Even granting that ceiling, the "
      f"geomean fold error only reaches {gub['gm_at_upper_bound']:.2f}x and "
      f"{gub['n_liver_below_floor_at_upper_bound']}/5 liver ratios still lie below "
      f"the reachable floor. This is the same *kind* of argument that falsified "
      "the transport model, now applied to its proposed replacement.")
    a("3. **On glioma the term is blocked on the one drug that mattered, twice "
      "over.** Adavosertib is a properly controlled OATP1A2 substrate (PMID "
      "28928160, HEK293 vs vector P = 0.001, abolished by 100 uM "
      "sulfobromophthalein P = 0.005, with OATP2B1 and OAT3 negative in the same "
      "experiment) — but the uptake fold appears **only as bars in Figure 4A**, "
      "with no value in the text, tables or abstract. And the one number that does "
      "exist, Li 2017's human GBM `Kp,uu` of 3.2, is circular against the "
      "scorecard target *and* is itself back-fitted: that paper's PBPK uptake "
      "clearance of 44 L/h was \"assigned based on sensitivity analysis\", i.e. "
      "tuned so the model reproduces 3.2. Passive plus efflux alone give 0.78-0.92.")
    a("4. **A widely repeated mechanistic story does not survive reading the "
      "source.** \"The acidic tumour microenvironment activates OATP1A2 uptake\" is "
      "asserted in `results/efflux_impact.md` on Li 2017's authority. In fact **no "
      "pH experiment was ever run on OATP1A2** in that paper — all of its "
      "pH-dependence work is on MDCKII efflux, and its own mechanism is *reduced "
      "efflux* at low pH (ER 9.2 -> 3.7 for ABCB1), not activated uptake. The PBPK "
      "model explicitly holds uptake constant. This run corrects the record.")
    a(f"5. **The redirect is the most useful output, and it is a correction to the "
      f"chain's own reasoning.** `results/moncayo_explant.md` eliminated "
      f"ion-trapping by proving a uniform retardation `R` gives `C = R*C_bath` "
      f"everywhere, so it cannot produce *spatial dispersion within a section*. "
      f"That proof is correct — **and it does not apply to a ratio BETWEEN two "
      f"tissues**, which is what GIST measures. Meanwhile the group that "
      f"dismantled the OCT1 story concluded that **lysosomal sequestration** "
      f"determines intracellular imatinib, tested in GIST-T1 among other lines "
      f"(PMID 26108972). Sizing it with imatinib's own pKa: the required "
      f"hepatocyte lysosomal volume fraction is "
      f"**{lys['phi_lys_required_in_hepatocytes']*100:.1f}%**, inside the range "
      f"`moncayo_explant.py` already sweeps. §10.")
    a("")

    # ---------------------------------------------------------- 1 additivity
    a("## 1. Nothing outside the two new files was touched")
    a("")
    a("The efflux run added one default-off keyword to "
      "`treatment.drugs.delivery_source`. **This run adds nothing at all** outside "
      "`tmesim/molecular/uptake.py` and `tmesim/validation/uptake_impact.py`: no "
      "shipped code path imports `molecular.uptake`, so the term is off by "
      "construction rather than by default. Both halves are verified numerically.")
    a("")
    a(f"- **{add['n_delivery_bit_identical']}/{add['n_delivery_comparisons']} "
      f"delivery-source comparisons bit-for-bit identical** (`np.array_equal`, not "
      f"`allclose`) across 4 drug names x 5 penetrances x 3 plasma concentrations "
      f"x 2 branches. All identical: **{add['delivery_all_bit_identical']}**.")
    a(f"- **{add['n_composition_bit_identical']}/"
      f"{add['n_composition_comparisons']} composition comparisons bit-for-bit "
      f"identical** to `reconstruct.barrier.effective_penetrance` wherever the "
      f"factor is exactly 1.0 — ceiling included. All identical: "
      f"**{add['composition_all_bit_identical']}**.")
    a(f"- Shipped files modified: **{add['n_shipped_files_modified']}**.")
    a("")
    a("And the checks are not vacuous. Each uses a drug/tissue pair that actually "
      "carries a number — doxorubicin in glioma, OATP1A2, 2.0-fold, PMID 24554572 "
      "— because a non-vacuity check run on a pair whose factor is 1.0 would pass "
      "trivially:")
    a("")
    a(f"- the factor bites where the carrier IS expressed: "
      f"**{add['factor_bites_where_carrier_expressed']}**")
    a(f"- the composition is genuinely tissue-differential (glioma > normal liver "
      f"for the same drug): **{add['composition_is_tissue_differential']}**")
    a(f"- the free-drug ceiling really does rise above 1.0: "
      f"**{add['ceiling_rises_above_one']}**")
    a("")
    a("That last one is a real capability even though nothing on the scorecards "
      "gets to use it. `reconstruct.barrier.effective_penetrance` clips penetrance "
      "at 1.0 and its docstring has always said why: *\"a partition above 1 "
      "requires tissue binding or active uptake, neither of which this term "
      "models.\"* This module models active uptake, so it raises the ceiling to the "
      "uptake factor — not a free choice but the steady-state bound "
      "`Kp,uu -> (P + CL_up)/P = U`. Adavosertib's measured 2.46 sits above the old "
      "ceiling and is unreachable without it. The machinery is in place; the "
      "literature will not supply the number to drive it.")
    a("")

    # ------------------------------------------------------------- 2 the form
    a("## 2. The functional form — and why there is only one of it")
    a("")
    a("Passive permeability `P` in both directions plus an active influx clearance "
      "`CL_up` in parallel gives a steady-state unbound partition")
    a("")
    a("```")
    a("  f_uptake = Kp,uu = (P + CL_up)/P = 1 + b_up = UR")
    a("```")
    a("")
    a("The efflux module had to choose between `1/ER` and `2/(1+ER)` — forms that "
      "differ by 2x at high ER — and justify the choice on n = 2. **That ambiguity "
      "does not arise here.** All three ways an uptake ratio is measured collapse "
      "to the same number:")
    a("")
    a("| observable | value |")
    a("|---|---|")
    a("| cell accumulation fold, transfected / vector control | `1 + CL_up/P` |")
    a("| directional `Papp(A->B) / Papp(B->A)` | `1 + CL_up/P` |")
    a("| `Papp(A->B)` vs the passive-only monolayer | `1 + CL_up/P` |")
    a("")
    a("Derived from the three-compartment mass balance and **verified numerically** "
      "in `molecular/uptake.py`'s self-test, not asserted. An in vivo tissue "
      "`Kp,uu` is the same quantity read in vivo, so it needs no conversion either. "
      "One form, no form-shopping — one fewer discretionary choice than the efflux "
      "term had, which is why the sensitivity table in §9 has no form axis.")
    a("")
    a("Carriers acting in parallel on one membrane have additive clearances, so "
      "`UR_combined = sum_i UR_i - (n - 1)`, the exact mirror of the efflux "
      "module's parallel-pump rule.")
    a("")
    a("**Zero free parameters.** Every ratio in the table is a published "
      "measurement with a PMID. Nothing is fitted to either scorecard. The one "
      "scaling factor that WAS available — converting transporter protein "
      "abundance in fmol/ug into a transport clearance — is refused for exactly "
      "the reason the efflux module refused its `SF`: the only data we could "
      "calibrate it on are the two scorecards under test.")
    a("")

    # ------------------------------------------------- 3 tissue expression map
    a("## 3. The tissue expression map — what makes this term non-scalar")
    a("")
    a("A transporter contributes in a tissue if and only if that tissue is listed "
      "here. This is the property efflux lacked: a per-drug scalar cancels out of "
      "a within-drug tissue ratio, `(f E_t)/(f E_n) = E_t/E_n`, whereas a "
      "per-(drug, tissue) factor scales it by `U_t/U_n`. Both identities are "
      "verified in the module's self-test.")
    a("")
    a("| transporter | gene | expressed in | source |")
    a("|---|---|---|---|")
    for name, t in TRANSPORTERS.items():
        a(f"| {name} | {t.gene} | {', '.join(t.tissues) if t.tissues else '**(none)**'} "
          f"| {t.tissue_citation or '*no verifiable localisation citation — see notes*'} |")
    a("")
    a("Four honest caveats, all of which cut against the result this map was built "
      "to produce:")
    a("")
    a("1. **Binary is a modelling choice.** The alternative — scaling transport by "
      "measured protein abundance — needs a transport-per-fmol conversion fittable "
      "only against the datasets under test. The binary map is a literature fact; "
      "a scaled map would be a fitted parameter.")
    a("2. **`OATP1B3` is NOT tumour-absent in general.** It is liver-restricted in "
      "*normal* tissue (PMID 10779507) but abundantly expressed in gastric, colon "
      "and pancreatic cancer (PMID 11375950). It is coded as absent from the GIST "
      "liver metastasis because GIST is a mesenchymal KIT-driven tumour rather "
      "than an adenocarcinoma — an inference from lineage, since **no measurement "
      "of any uptake carrier in GIST tissue exists** (PubMed `SLCO1A2 AND GIST` "
      "returns zero; the one GIST transporter review, PMID 21163869, covers only "
      "ABCB1/ABCC1).")
    a("3. **`CTR1` and `SLC22A16` are encoded with NO tissues**, and therefore "
      "contribute nothing. CTR1 is universally described as the cisplatin uptake "
      "route but no verifiable primary human tissue-distribution citation could be "
      "retrieved; SLC22A16 carries the cleanest doxorubicin uptake Km anywhere "
      "(5.2 uM, PMID 15963465) but was cloned as testis-specific (PMID 12089149). "
      "Both are recorded and neutralised rather than silently omitted.")
    a("4. **`OATP1A2` is deliberately NOT given the intestine**, and this "
      "exclusion is not neutral. OATP1A2 is widely asserted to sit in the "
      "enterocyte apical membrane and to mediate duodenal absorption of imatinib "
      "(PMID 21508937), but no verified human intestinal localisation citation "
      "could be retrieved, and the module forbids an uncited tissue claim. "
      "Including it would move the GIST intestinal control pair — the one pair the "
      "shipped model already gets right. Sources also disagree on brain: PMID "
      "16357150 places OATP1A2 in glioma BBB/blood-tumour-barrier endothelium and "
      "explicitly **not** in glioma cells (which is exactly where a penetrance "
      "multiplier belongs), while PMID 25132355 localises it to neurons and puts "
      "OATP2B1 in the capillaries instead.")
    a("")

    # ------------------------------------------------------------- 4 the table
    a("## 4. The literature-sourced substrate table")
    a("")
    a(f"{s['n_drugs']} drugs and {len(TRANSPORTERS)} carriers: all ten in "
      "`molecular.library.GBM_DRUGS`, the four glioma-scorecard drugs outside it, "
      "imatinib (GIST), and the three PARP inhibitors (Moncayo) — the same 18 the "
      "efflux table covers, so the two mechanisms are scored over an identical "
      "library. Every non-`unknown` cell carries a PMID.")
    a("")
    a(f"- confirmed substrate of at least one carrier: **{s['substrate_any']}**")
    a(f"- .. WITH a usable published number: **{s['substrate_with_number']}**")
    a(f"- .. WITHOUT one (factor forced to exactly 1.0): "
      f"**{s['substrate_without_number']}**")
    a(f"- .. whose only number is withheld as circular: "
      f"**{s['circular_suppressed']}**")
    a(f"- measured non-substrate: **{s['non_substrate']}**")
    a(f"- unknown status throughout: **{s['unknown_only']}**")
    a("")
    a("**This is the same binding constraint the efflux run hit, and it is worse "
      "here.** There, 7 of 18 drugs carried a usable published ratio. Here it is "
      f"{s['substrate_with_number']} of {s['n_drugs']} — and the field's reason is "
      "structural rather than accidental: uptake-transporter papers report "
      "**Km/Vmax and inhibitor IC50s**, which characterise the carrier, while a "
      "penetrance multiplier needs a **ratio against passive permeability in the "
      "same system**. Methotrexate is the extreme case: a confirmed substrate of "
      "six carriers spanning liver, kidney, BBB and every tissue, with abundant "
      "kinetics, and not one usable fold.")
    a("")
    hdr_t = [t for t in TISSUES if t != TISSUE_GENERIC]
    a("| drug | " + " | ".join(hdr_t) + " | reason (glioma) |")
    a("|---|" + "---|" * (len(hdr_t) + 1))
    for r in subs:
        cells = []
        for t in hdr_t:
            v = r["per_tissue"][t]
            cells.append(f"**{v['factor']:.2f}**" if v["factor"] > 1.0 else "1.00")
        a(f"| {r['drug']} | " + " | ".join(cells)
          + f" | {r['per_tissue'][TISSUE_GLIOMA]['reason']} |")
    a("")
    a("### Citations, per assignment")
    a("")
    for r in subs:
        a(f"**{r['drug']}**")
        if not r["entries"]:
            a("  - no carrier assignment: see notes.")
        for e in r["entries"]:
            line = f"  - `{e.transporter}` **{e.status}**"
            if e.citation:
                line += f": {e.citation}"
            a(line)
            if e.uptake_ratio is not None:
                a(f"    - in vitro uptake ratio **{e.uptake_ratio:.2f}** — "
                  f"{e.ratio_assay} — {e.ratio_citation}")
            if e.kpuu_in_vivo is not None:
                flag = (f" — **WITHHELD as circular against "
                        f"{', '.join(e.circular_against)}**" if e.circular_against else "")
                a(f"    - in vivo Kp,uu **{e.kpuu_in_vivo:.2f}** ({e.kpuu_tissue}) — "
                  f"{e.kpuu_assay} — {e.kpuu_citation}{flag}")
            if e.notes:
                a(f"    - *{e.notes}*")
        a(f"  *Notes:* {r['notes']}")
        a("")
    a("### Corrections this pass makes to earlier reports in this repo")
    a("")
    a("Three, all recorded because a repo that only ever adds findings is not "
      "checking itself:")
    a("")
    a("1. **`results/efflux_impact.md` states adavosertib \"ACCUMULATES in human "
      "glioma because OATP1A2-mediated uptake dominates in the acidic tumour "
      "microenvironment\".** The accumulation is real and well measured. The "
      "*explanation* is not: no pH experiment was ever run on OATP1A2 in the "
      "source paper, its pH work is entirely on MDCKII efflux, its own mechanism "
      "is reduced efflux at low pH, and its PBPK model holds uptake clearance "
      "constant at a value fitted to reproduce the answer.")
    a("2. **`results/efflux_impact.md` records cefazolin as \"a confirmed OAT1 "
      "substrate (PMID 15618660, Ki 72 uM vs PAH)\".** Two errors compounded: that "
      "PMID is **rat** Oat1, and human OAT1 does not transport cefazolin (uptake "
      "within two-fold of control, PMID 16098483); and a Ki is *inhibition* of PAH "
      "transport, not substrate status. Cefazolin's real carrier is OAT3 (PMID "
      "14984259).")
    a("3. **`results/gist_imatinib.md` names \"active hepatic uptake (OCT1/OATP)\" "
      "as the missing physics.** The OCT1 half is falsified (§5.1). The OATP half "
      "survives as substrate status but is quantitatively insufficient (§5.2).")
    a("")

    # ------------------------------------------------------- 5 GIST scorecard
    a("## 5. Scorecard 1 — GIST tumour:normal, the target uptake was named for")
    a("")
    a("The prediction is the published one with one factor added per side:")
    a("")
    a("```")
    a("  before:  R = E_tumour / E_normal")
    a("  after :  R = (U_tumour * E_tumour) / (U_normal * E_normal)")
    a("```")
    a("")
    a(f"`E` is unchanged (1-D Krogh slab, L = {gist['L_um']:.1f} um from the "
      "doxorubicin-anchored `alpha` fitted once to an unrelated measurement), the "
      "intervessel half-spacings are unchanged. The ONLY new quantity is `U`.")
    a("")
    a("| organ | side | tissue | U | carriers applied | reason |")
    a("|---|---|---|---|---|---|")
    for organ, d in gist["factors"].items():
        for side in ("tumour", "normal"):
            v = d[side]
            tis = GIST_TISSUE_MAP[organ][0 if side == "tumour" else 1]
            a(f"| {organ} | {side} | {tis} | **{v['factor']:.3f}** | "
              f"{', '.join(v['transporters']) if v['transporters'] else '--'} | "
              f"{v['reason']} |")
    a("")
    a(f"**`U_normal` for liver resolves to 1.000**, so the prediction is unchanged "
      f"at {ga['geomean_fold_error']:.2f}x. The reason is not that imatinib fails "
      "to be a hepatic uptake substrate — OATP1B3 is a confirmed one — but that "
      "the only place its uptake fold was ever published is Figure 1A of PMID "
      "23340295, as bars with asterisks and no number in the text, tables or "
      "caption. This is the doxorubicin situation from `results/efflux_impact.md` "
      "repeated exactly, on the anchor drug of the other scorecard.")
    a("")
    a("| patient | organ | measured | pred before | fold before | pred after | "
      "fold after |")
    a("|---|---|---|---|---|---|---|")
    for rb, ra in zip(gb["per_patient"], ga["per_patient"]):
        a(f"| {rb['patient']} | {rb['organ']} | **{rb['measured']:.4f}** | "
          f"{rb['pred']:.4f} | {rb['fold']:.2f}x | {ra['pred']:.4f} | "
          f"{ra['fold']:.2f}x |")
    a("")

    a("### 5.1 The named mechanism is falsified, not merely unquantified")
    a("")
    a("`results/gist_imatinib.md` concluded: *\"The dominant term is not the tumour "
      "being hard to reach — it is normal liver actively concentrating a weak base "
      "via hepatic uptake transporters and lysosomal trapping.\"* The uptake half "
      "of that sentence was checked rather than adopted, and it does not survive.")
    a("")
    a("| design | result | source |")
    a("|---|---|---|")
    a("| 5 OCT1-expressing lines (HEK, MDCK, V79, CHO) + oocytes | **no increase "
      "over vector control**; MPP+ positive control up 8-fold in the SAME oocytes "
      "| PMID 24352644 |")
    a("| decynium-22 and prazosin on OCT1 cells | **failed to inhibit** imatinib "
      "uptake | PMID 24352644 |")
    a("| shRNA OCT1 knockdown | **no effect** | PMID 24352644 |")
    a("| **Oct1(-/-) mice, plasma and HEPATIC imatinib** | **no change vs "
      "wild-type** | PMID 24352644 |")
    a("| oocytes, 0.2 uM | **not a substrate, P = 0.21**; HEK293 only 1.20-fold | "
      "PMID 18483382 |")
    a("| reduced-function SLC22A1 variants, oral clearance, n = 74 | **no effect, "
      "P = 0.99** | PMID 18483382 |")
    a("| oocytes + HEK | **not a substrate of OCT1, OCT2 or OCT3** | PMID 30734690 |")
    a("")
    a("The large clinical \"OCT-1 activity\" literature (PMID 16597591, 17761829, "
      "17568400; PMID 18559609 gives a prazosin-inhibitable IUR of 6.38 "
      "ng/200,000 cells) is real and genuinely predicts response — but every one "
      "of those studies measures an **inhibitor-difference in non-transfected "
      "cells**, never transport by OCT1. PMID 26108972 shows in HEK293, K562, SD-1 "
      "and **GIST-T1** that the prazosin-sensitive component is **lysosomal "
      "sequestration**, not SLC22A1. And imatinib inhibits OCT1 with IC50 0.095 uM "
      "— the precise property that makes an inhibitor-difference assay read like "
      "substrate evidence. **Inhibitor is not substrate**, and that single "
      "distinction is the difference between this row being \"no\" and being "
      "\"yes\".")
    a("")

    a("### 5.2 And the rest of the class is excluded by its own numbers")
    a("")
    a("A factor of 1.0 for want of a published number leaves open the possibility "
      "that the mechanism is right and merely unmeasured. So the sharper question "
      "is asked: **what happens at the largest value any published source can "
      "support?**")
    a("")
    a(f"- Upper bound used: **{gub['upper_bound']:.1f}x**. Source: "
      f"{gub['citation']}.")
    a(f"- Value the measurement DEMANDS: **{gub['required_u_normal']:.2f}x**.")
    a(f"- **Shortfall: {gub['shortfall_fold']:.1f}x.**")
    a("")
    a("| U_normal | predicted liver ratio | geomean fold error (n = 6) | within 3x |")
    a("|---|---|---|---|")
    for r in gub["sweep"]:
        tag = ""
        if r["is_upper_bound"]:
            tag = " **<- literature ceiling**"
        if r["is_required"]:
            tag = " **<- what the data demands**"
        a(f"| {r['u_normal']:.2f}{tag} | {r['pred_liver']:.4f} | "
          f"{r['geomean_fold_error']:.2f}x | {r['n_within_3x']}/6 |")
    a("")
    a(f"Measured liver geomean: **{gub['measured_liver_geomean']:.4f}**. At the "
      f"literature ceiling the prediction is "
      f"{gub['pred_liver_at_upper_bound']:.4f} — still "
      f"{gub['pred_liver_at_upper_bound']/gub['measured_liver_geomean']:.1f}x too "
      f"high — and the structural floor only falls from "
      f"{gist['floor_nominal_before']:.3f} to {gub['floor_at_upper_bound']:.3f}, "
      f"leaving {gub['n_liver_below_floor_at_upper_bound']}/5 measured liver ratios "
      "still below it.")
    a("")
    a(f"**Excluded by its own literature: {gub['excluded_by_its_own_literature']}.** "
      "This is the same *class* of argument that made the original GIST run "
      "valuable — there, no penetration length could reach the measurement; here, "
      "no published hepatic uptake ratio can. A mechanism can be ruled out by a "
      "bound without ever being implemented at a value, and doing so is cheaper "
      "and more honest than fitting one.")
    a("")

    a("### 5.3 Does the GIST verdict survive the choices we made?")
    a("")
    a("The published run showed the intervessel half-spacings are the largest "
      "single lever on the transport prediction. Every row re-solves the whole "
      "prediction.")
    a("")
    a("| axis | setting | U_normal | GM before | GM after | improved? |")
    a("|---|---|---|---|---|---|")
    for r in gsens:
        a(f"| {r['axis']} | {r['setting']} | {r['u_normal_liver']:.2f} | "
          f"{r['gm_before']:.2f}x | {r['gm_after']:.2f}x | "
          f"{'yes' if r['improved'] else 'no'} |")
    n_imp = sum(1 for r in gsens if r["improved"])
    a("")
    a(f"**{n_imp} of {len(gsens)} configurations improve the GIST prediction.**")
    a("")
    a("Read that table for what it is: **every row is identical before and after, "
      "because `U_normal` resolves to 1.000 in all of them.** It shows the verdict "
      "is not an artefact of the geometry, the source preference or the "
      "circularity guard, and it shows the geometry lever is real (the published "
      "7.92x moves between 3.67x and 12.87x across the band) — but it cannot say "
      "anything about the mechanism, because the mechanism never engages. The "
      "sweep that carries information is the upper-bound table in §5.2, where "
      "`U_normal` is forced through the whole range and cannot reach the "
      f"measurement below {gistreq['required_u_normal']:.1f}x.")
    a("")

    a("### 5.4 What this does NOT establish, at n = 6")
    a("")
    a("The published run's caveat applies unchanged: **n = 6 matched pairs, five "
      "of them the same tissue pairing**, i.e. roughly ONE independent "
      "discriminating comparison. That was enough to falsify a transport model and "
      "it is enough to bound a replacement, but:")
    a("")
    a("- The upper-bound argument inherits every assumption of the delivery model "
      "it sits on top of — the Krogh geometry, the literature half-spacings, and "
      "the doxorubicin-anchored `alpha`. It shows hepatic uptake cannot rescue "
      "*this* model, not that hepatic uptake is unimportant.")
    a("- Ion suppression differs between liver parenchyma and GIST stroma and does "
      "not cancel in the ratio. It is an uncontrolled confounder inflating the "
      "measured contrast, and therefore inflating the required `U_normal`.")
    a("- `U_tumour = 1.0` is an inference from lineage and from HCC "
      "downregulation (§8), not a measurement on these sections.")
    a("")

    # ----------------------------------------------------- 6 glioma scorecard
    a("## 6. Scorecard 2 — the six-drug measured human glioma partition")
    a("")
    a("BEFORE is not pasted: it is "
      "`glioma_regional_pk.run_scorecard(aggregation=\"geometric\", "
      "apply_leak=True)` called live. AFTER changes one thing:")
    a("")
    a("```")
    a("  before:  P_eff = clip( bbb_penetrance     + L * D_sample , 0, 1 )")
    a("  after :  P_eff = clip( bbb_penetrance * U + L * D_sample , 0, U )")
    a("```")
    a("")
    a("| drug | measured | struct. penetrance | disrupted frac | **U** | why | "
      "pred before | fold before | pred after | fold after |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for r in sc["rows"]:
        a(f"| {r['drug']} | {r['measured']:.4f} | {r['structure_only']:.4f} | "
          f"{r['sample_disruption']:.2f} | **{r['uptake_factor']:.4f}** | "
          f"{r['uptake_reason']} | {r['pred_before']:.4f} | "
          f"{r['fold_before']:.3f}x | {r['pred_after']:.4f} | "
          f"{r['fold_error']:.3f}x |")
    a("")
    a("| aggregate | before | after | direction |")
    a("|---|---|---|---|")
    a(f"| GM \\|fold error\\| | **{b['gm_abs']:.3f}x** | **{af['gm_abs']:.3f}x** | "
      f"{_dir(b['gm_abs'], af['gm_abs'])} |")
    a(f"| GM signed bias | {b['gm_signed']:.3f}x | {af['gm_signed']:.3f}x | |")
    a(f"| within 3x | {b['n_within_3x']}/6 | {af['n_within_3x']}/6 | |")
    a(f"| within 10x | {b['n_within_10x']}/6 | {af['n_within_10x']}/6 | |")
    a(f"| worst \\|fold error\\| | {b['worst_abs_fold']:.2f}x | "
      f"{af['worst_abs_fold']:.2f}x | |")
    a("")
    a(f"**{sc['n_drugs_moved']} of 6 drugs moves at all.** The `why` column gives "
      "which of the eight distinct reasons produced each 1.0, and the distinctions "
      "matter: `non_substrate` (measured not to be one) is a finding, "
      "`substrate_but_no_published_ratio` is a coverage gap, "
      "`no_transporter_in_this_tissue` is the mechanism working correctly, and "
      "`only_circular_observable` is a number we refused.")
    a("")

    a("### 6.1 Adavosertib — the pre-registered target, blocked twice")
    a("")
    a("The efflux run identified adavosertib as the case where efflux was most "
      "clearly the wrong mechanism and named uptake as the right one. This run was "
      "written to test that. It cannot.")
    a("")
    a("**The substrate status verifies.** Li et al. did not merely hypothesise it: "
      "uptake of AZD1775 in HEK293 overexpressing OATP1A2 was significantly higher "
      "than vector control (P = 0.001) and was reduced to control level by 100 uM "
      "sulfobromophthalein (P = 0.005), with an estrone-3-sulfate positive control "
      "and with OATP2B1 and OAT3 negative in the same experiment. That is a "
      "properly controlled positive, and PMID 16357150 puts OATP1A2 in the glioma "
      "blood-tumour-barrier endothelium — exactly where a penetrance multiplier "
      "acts.")
    a("")
    a("**Block 1: the number does not exist.** The uptake fold appears only as "
      "bars in Figure 4A. No value in the text, tables or abstract; no Km, no "
      "Vmax; PubMed `AZD1775 AND OATP1A2` returns exactly one record and it "
      "publishes only P-values. Under this repo's standing rule — the one that "
      "forced doxorubicin to 1.0 in `results/efflux_impact.md` — a figure is not a "
      "number.")
    a("")
    a("**Block 2: the only number that does exist is circular, and back-fitted.** "
      "See §7.")
    a("")
    a(f"So adavosertib's factor is 1.00 and its fold error is unchanged at "
      f"{sc['adavosertib_fold_after']:.3f}x. **The mechanism is right, the "
      "measurement is right, and the model can use neither.**")
    a("")

    a("### 6.2 Temozolomide is the only mover, and it moves the wrong way")
    a("")
    tmz = [r for r in sc["rows"] if r["drug"] == "temozolomide"][0]
    a(f"TMZ carries the table's other usable ratio — GLUT3, {tmz['uptake_factor']:.2f}x "
      "(PMID 40813374) — and it makes the prediction worse, "
      f"{tmz['fold_before']:.3f}x -> {tmz['fold_error']:.3f}x, taking a drug that "
      "was slightly under-predicted to slightly over-predicted. This is the "
      "temozolomide-and-efflux pattern in reverse: there, a genuine but tiny "
      "1.5-fold efflux term hurt because TMZ was already under-predicted; here a "
      "1.70-fold uptake term overshoots.")
    a("")
    a("It is also the weakest row in the table and is reported as such: a single "
      "unreplicated 2025 paper, whose companion affinity data are internally "
      "inconsistent by ~1000x (GLUT3 glucose Kd reported at 1-10 uM against an "
      "established Km of ~1.4 mM), and whose knockdown arm implies a larger ratio "
      "(2.19) than the overexpression arm used here (1.70). Leave-one-out below "
      "shows the scorecard verdict does not depend on it either way.")
    a("")

    a("### 6.3 Which drug is doing what")
    a("")
    a("| drug dropped | GM \\|fold\\| before | GM \\|fold\\| after | "
      "uptake helps without it? |")
    a("|---|---|---|---|")
    for r in loo:
        a(f"| {r['dropped']} | {r['gm_abs_before']:.3f}x | {r['gm_abs_after']:.3f}x "
          f"| {'YES' if r['improved_without_it'] else 'no'} |")
    a("")

    a("### 6.4 Is the change even distinguishable at n = 6?")
    a("")
    a(f"- drugs that move at all: **{paired['n_movers']} of {paired['n_drugs']}**; "
      f"of those, {paired['n_worse']} worsen and {paired['n_better']} improve")
    a(f"- mean change in `|log10 fold error|`: **{paired['mean_delta_log10']:+.4f}** "
      f"(positive = worse)")
    a(f"- exact sign test: p = {_fmt(paired['p_sign_test'])}; Wilcoxon signed-rank: "
      f"p = {_fmt(paired['p_wilcoxon'])}")
    a("")
    a(f"**{paired['note']}.** No p-value here should be read as evidence for "
      "anything. What the scorecard establishes is not that uptake hurts; it is "
      "that the term cannot be *exercised* on five of the six drugs, and why each "
      "one is blocked.")
    a("")

    # ------------------------------------------------------- 7 the guard cost
    a("## 7. The circularity guard, and exactly what it costs")
    a("")
    a("A number one has decided not to use should be printed, not hidden.")
    a("")
    if not circ:
        a("*(nothing withheld)*")
    else:
        a("| drug | carrier | tissue | withheld value | quantity | circular "
          "against | factor default | factor if allowed |")
        a("|---|---|---|---|---|---|---|---|")
        for r in circ:
            a(f"| {r['drug']} | {r['transporter']} | {r['tissue']} | "
              f"**{r['withheld_value']:.2f}** | {r['withheld_quantity'][:90]}... | "
              f"{', '.join(r['circular_against'])} | {r['factor_default']:.2f} | "
              f"{r['factor_if_allowed']:.2f} |")
    a("")
    a("**Why this is circular and not merely \"in vivo\".** The efflux module used "
      "in vivo mouse knockout fold changes without complaint, and rightly: a "
      "knockout contrast in a mouse is a mechanistic input, measured against a "
      "genetic control, in a different species from the target. Li 2017's `Kp,uu` "
      "is none of those things. It is a human glioblastoma tumour-to-plasma drug "
      "concentration ratio — the same quantity, in the same tissue, in the same "
      "disease as the scorecard target.")
    a("")
    a("The two studies ARE independent (Li 2017: Karmanos/Barrow, n = 20, unbound "
      "ratio, median 3.2 range 1.3-24.4; Lee 2025 PMID 39820282: Adult Brain Tumor "
      "Consortium, total homogenate ratio 4.18 +/- 3.36 enhancing and 0.74 +/- "
      "0.63 non-enhancing), so this is not double-counting one dataset. Enabling "
      "it would be a real cross-cohort replication. **It would not be a prediction "
      "from mechanism**, which is the only thing this repo scores.")
    a("")
    a("And there is a second, independent reason to refuse it. Li 2017's PBPK "
      "uptake clearance of **44 L/h was \"assigned based on sensitivity "
      "analysis\"** — i.e. tuned so the model reproduces 3.2. Passive plus efflux "
      "alone give Kp,uu 0.78-0.92 (max 1.03 at sensitivity extremes). So the value "
      "is not an independent measurement of uptake at all; it is the residual "
      "required to close the mass balance, and importing it would import a fitted "
      "parameter under the name of a literature value.")
    a("")
    sens_circ = [r for r in sens
                 if r["allow_circular"] and r["placement"] == "transcellular"
                 and r["ceiling"] == "uptake" and r["prefer"] == PREFER_IN_VITRO]
    if sens_circ:
        r = sens_circ[0]
        a(f"For completeness, the number it would have bought: GM abs fold error "
          f"{r['gm_abs_before']:.3f}x -> **{r['gm_abs_after']:.3f}x**, adavosertib "
          f"to {r['ada_after']:.3f}x. **That is an upper bound obtained by using "
          "the answer, and it is not claimed as a result.**")
        a("")
    a("Note what this implies about the headline. If the guard were dropped, the "
      "glioma scorecard would improve — and the improvement would be worthless. "
      "That is the whole reason the guard is a code path with a reason code rather "
      "than a paragraph of good intentions.")
    a("")

    # ------------------------------------------------------- 8 external check
    a("## 8. The mechanism's premise, tested out of sample")
    a("")
    a("The GIST argument rests on a claim the GIST dataset cannot test, because no "
      "transporter stain accompanies the deposit: **hepatic uptake carriers are "
      "present in normal liver and effectively absent from a tumour sitting in "
      "it.** The hardest possible case for that claim is hepatocellular carcinoma "
      "— a tumour that actually derives from hepatocytes. If the premise failed "
      "even there, it would certainly fail for a GIST metastasis.")
    a("")
    a("| gene | carrier | tumour vs adjacent non-tumour liver | source |")
    a("|---|---|---|---|")
    for r in hcc["rows"]:
        a(f"| {r['gene']} | {r['transporter']} | {r['finding']} | {r['citation']} |")
    a("")
    a(f"Uniformly downregulated across all rows: **{hcc['all_downregulated']}** — "
      "and the dissent is in the table rather than dropped: OATP1B1 was **not** "
      "significantly decreased in the same paired samples that found OATP1B3 down "
      "in 60% (PMID 14739090). The premise is directionally supported and not "
      "uniform. No number from these papers enters the model, and HCC is not GIST.")
    a("")

    # ---------------------------------------------------------- 9 sensitivity
    a("## 9. Does the glioma verdict depend on anything we chose?")
    a("")
    a("Every discretionary choice, crossed: source preference, the circularity "
      "guard, where the factor multiplies, and which ceiling applies. Sixteen "
      "configurations. There is no functional-form axis because the form is unique "
      "(§2) — the efflux run needed eight configurations partly to cover a form "
      "ambiguity that does not exist here.")
    a("")
    a("| prefer | circular | placement | ceiling | GM before | GM after | 3x | 10x "
      "| worst after | ada after | improved? |")
    a("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in sens:
        a(f"| {r['prefer']} | {r['allow_circular']} | {r['placement']} | "
          f"{r['ceiling']} | {r['gm_abs_before']:.3f}x | {r['gm_abs_after']:.3f}x | "
          f"{r['within_3x_before']}->{r['within_3x_after']} | "
          f"{r['within_10x_before']}->{r['within_10x_after']} | "
          f"{r['worst_after']:.2f}x | {_fmt(r['ada_after'])}x | "
          f"{'yes' if r['improved'] else 'no'} |")
    honest = [r for r in sens if not r["allow_circular"]]
    n_honest = sum(1 for r in honest if r["improved"])
    n_all = sum(1 for r in sens if r["improved"])
    a("")
    a(f"**{n_honest} of {len(honest)} non-circular configurations improve the "
      f"glioma scorecard** ({n_all} of {len(sens)} including the circular ones). "
      "The verdict is not an artefact of the placement, the ceiling or the source "
      "preference.")
    a("")

    # ------------------------------------------------------------ 10 redirect
    a("## 10. The redirect — and a correction to this chain's own reasoning")
    a("")
    a("Two mechanisms have now been implemented, cited and scored against these "
      "datasets, and neither worked. That would be a dead end if the elimination "
      "argument that got us here were sound. It is not, in one specific place, and "
      "that is the most useful thing this run found.")
    a("")
    a("**What Moncayo proved.** `results/moncayo_explant.md` §4.3: with `phi_lys` "
      "uniform in space, the retardation `R` is uniform, so at equilibrium "
      "`C_total = R * C_bath` **everywhere** and the dispersion goes to zero "
      "however large `R` becomes. Correct, and this module does not dispute it.")
    a("")
    a("**What it does not prove.** That is a statement about dispersion *within one "
      "tissue*. The GIST target is a ratio *between two tissues*, and two tissues "
      "with different acidic-compartment content have different `R`:")
    a("")
    a("```")
    a("  C_tumour / C_normal = (R_t * E_t) / (R_n * E_n)      does NOT cancel")
    a("```")
    a("")
    a("This is the same algebra that lets a tissue-differential uptake term move "
      "the ratio where a scalar efflux term could not. **The chain eliminated "
      "ion-trapping for GIST using a proof about a different quantity.** Both the "
      "efflux report (§5.1) and this task's framing inherited that error.")
    a("")
    a("**And the evidence points there.** The same group whose experiments "
      "falsified imatinib/OCT1 concluded that **lysosomal sequestration, not "
      "SLC22A1, determines intracellular imatinib** — in HEK293, K562, SD-1 and, "
      f"to the point, **GIST-T1** ({lys['citation_mechanism'].split(': ')[0]}). The "
      "prazosin-sensitive signal underpinning the entire clinical \"OCT-1 "
      "activity\" literature is, in that work, lysosomal. So the mechanism that "
      "explains why the uptake story looked true is the mechanism the chain "
      "discarded.")
    a("")
    a("### Sizing it — is the required asymmetry reachable, or absurd?")
    a("")
    a("Not a model, not a prediction, and nothing fitted. This uses the "
      "Henderson-Hasselbalch trapping ratio already implemented in "
      "`validation/moncayo_explant.py` and imatinib's own basic pKa already used "
      "by `validation/gist_imatinib.py`. An absurd answer would rule the mechanism "
      "out cheaply; a plausible one identifies the next piece of work.")
    a("")
    a("| quantity | value |")
    a("|---|---|")
    a(f"| imatinib basic pKa (already in the repo) | {lys['pka_base']} |")
    a(f"| lysosome:cytosol trapping ratio at pH 4.7 / 7.2 | "
      f"**{lys['lysosome_cytosol_trapping_ratio']:.1f}** |")
    a(f"| `R_normal/R_tumour` required by the measurement | "
      f"**{lys['required_R_normal_over_tumour']:.2f}** |")
    a(f"| implied hepatocyte lysosomal volume fraction (tumour ~0) | "
      f"**{lys['phi_lys_required_in_hepatocytes']*100:.2f}%** |")
    a(f"| `phi_lys` the Moncayo module uses for generic tissue | "
      f"{lys['phi_lys_reference_generic_tissue']*100:.1f}% |")
    a(f"| largest `phi_lys` the Moncayo module swept | "
      f"{lys['moncayo_swept_max_phi']*100:.0f}% |")
    a(f"| **reachable inside that swept range?** | "
      f"**{lys['plausible']}** |")
    a("")
    a(f"The required {lys['phi_lys_required_in_hepatocytes']*100:.1f}% sits inside "
      "the range `moncayo_explant.py` already sweeps and is high but not absurd "
      "for a hepatocyte, which is among the most lysosome-rich cells in the body. "
      "Note carefully what this is and is not: **it is a plausibility check on a "
      "single number, not evidence.** It says the candidate cannot be dismissed on "
      "magnitude, which is exactly the test the hepatic-uptake candidate just "
      "failed in §5.2. That asymmetry — one candidate excluded by a bound, the "
      "other not — is the entire content of this section.")
    a("")
    a("**Two things would be needed to turn it into a test**, and neither is "
      "available today: a measured acidic-compartment or lysosomal-content "
      "contrast between hepatocytes and GIST tissue (a literature pass for this "
      "project established that no measured mg/g or mol% acidic-phospholipid "
      "content exists for tumour tissue, and no published PBPK model carries "
      "tumour-specific lysosomal parameters), and an independent tumour-in-liver "
      "dataset so the test is not run on the same six pairs that generated the "
      "hypothesis.")
    a("")

    # ------------------------------------------------------------- 11 verdict
    a("## 11. Verdict and recommendation")
    a("")
    a("### Recommendation: **default OFF, and more strongly than efflux.**")
    a("")
    a("Efflux was left off after being tested and found wrong. Uptake is left off "
      "after being tested and found **inapplicable**: on 16 of 18 drugs the "
      "literature supplies no usable ratio, on the GIST anchor the named carrier "
      "is falsified and the class is bounded out, and on the glioma anchor the one "
      "usable number would be circular. Nothing in the shipped pipeline calls "
      "`molecular.uptake` and nothing was changed to accommodate it; it is "
      "available to any caller who imports it, as a documented mechanism rather "
      "than an endorsed default.")
    a("")
    a("**Do not read this as \"uptake transport is unimportant.\"** It plainly is "
      "important — a 4.1-fold liver-to-plasma change in Oatp1a/1b-/- mice for "
      "doxorubicin is real. What is shown is narrower and firmer: a "
      "literature-sourced, zero-free-parameter uptake term cannot be built for "
      "these drugs today, and where it can be built it does not improve these two "
      "measured human datasets.")
    a("")
    a("### What this run establishes")
    a("")
    a("1. **Uptake is implementable with zero free parameters and, unlike efflux, "
      "one unambiguous functional form** — all three published observables collapse "
      "to `UR = 1 + CL_up/P`, verified from the mass balance.")
    a("2. **A tissue-differential factor is structurally capable of moving a "
      "tumour:normal ratio**, where the efflux scalar provably was not. Both "
      "identities are verified to machine precision. The architecture is correct; "
      "the inputs are missing.")
    a("3. **The mechanism `results/gist_imatinib.md` named is falsified** by five "
      "cell lines, oocytes with a positive control, shRNA knockdown, and Oct1(-/-) "
      "mice with unchanged hepatic imatinib.")
    a(f"4. **The remaining hepatic-uptake class is excluded by a bound**: "
      f"{gub['required_u_normal']:.1f}x demanded, {gub['upper_bound']:.1f}x "
      f"available, {gub['shortfall_fold']:.1f}x short. Ruling a mechanism out by a "
      "literature ceiling, without ever fitting it, is the cheapest honest form of "
      "falsification available and it is the main methodological contribution here.")
    a("5. **Literature coverage is the binding constraint for the second run "
      "running, and now the reason is structural.** Uptake papers report Km/Vmax "
      "and inhibitor IC50s; a penetrance multiplier needs a ratio against passive "
      "permeability in the same system. Those are different measurements and the "
      "field mostly publishes the first.")
    a("6. **Three factual corrections to earlier reports in this repo** (§4), "
      "including one — the acidic-TME/OATP1A2 story — that this repo had been "
      "repeating on a misreading of its own cited source.")
    a("")
    a("### What this run does NOT establish")
    a("")
    a("- **Nothing at adequate power.** Six matched pairs with one independent "
      "tissue contrast; six drugs of which one moves.")
    a("- **Not that adavosertib's uptake mechanism is wrong.** It is documented, "
      "controlled, in the right direction, and unquantified. If the value behind "
      "Figure 4A of PMID 28928160 were published, the glioma half of this test "
      "could actually be run — one number, from one experiment already performed.")
    a("- **Nothing about the tumour side of the GIST contrast.** `U_tumour = 1.0` "
      "is an inference from lineage and HCC data, not a measurement.")
    a("- **Nothing about lysosomal sequestration.** §10 sizes it; it does not test "
      "it. A plausible magnitude is not evidence.")
    a("")
    a("### Where the chain stands now")
    a("")
    a("| candidate | status |")
    a("|---|---|")
    a("| passive reaction-diffusion transport | **falsified** for tissue partition "
      "— below the structural floor at any L (`gist_imatinib.md`) |")
    a("| scalar lysosomal ion-trapping, WITHIN a section | **falsified** — uniform "
      "R gives no dispersion at equilibrium (`moncayo_explant.md`) |")
    a("| active efflux (P-gp / BCRP) | **tested, does not help**; provably inert on "
      "a tissue ratio; wrong sign in liver (`efflux_impact.md`) |")
    a("| carrier-mediated uptake (OATP / OCT / OAT) | **tested, does not help**; "
      "named carrier falsified; class excluded by a literature bound (this run) |")
    a("| **tissue-differential acidic-compartment / lysosomal content** | "
      "**OPEN, and wrongly eliminated earlier.** Direct evidence in GIST-T1 (PMID "
      "26108972); required magnitude plausible (§10); no measurement exists |")
    a("| spatial transporter or acidic-compartment FIELD (within-tumour structure) "
      "| **OPEN**, unchanged from `moncayo_explant.md`; no measurement in usable "
      "units exists |")
    a("")
    a("The last two are different objects and should not be conflated again. A "
      "per-tissue *constant* is enough to move a tumour:normal ratio and cannot "
      "produce structure inside a tumour; a spatial *field* is needed for the "
      "Moncayo target and is not needed for the GIST one. This run's contribution "
      "to that distinction is that the GIST-shaped version of ion-trapping was "
      "never actually ruled out.")
    a("")
    a("---")
    a("")
    a("*Generated by `tmesim/validation/uptake_impact.py`. Additive only: no file "
      "under `core/`, `physics/`, `molecular/` (other than the new "
      "`molecular/uptake.py`), `reconstruct/`, `sim/` or `treatment/` was "
      "modified, and the bit-exactness of the untouched delivery path is verified "
      "in §1.*")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\n".join(A) + "\n")
    return path


# ===========================================================================
# 7. Entry point
# ===========================================================================
def run(write_md: bool = True) -> Dict:
    add = additivity_check()
    subs = substrate_rows()
    circ = circularity_audit()
    gist = gist_impact()
    gistreq = gist_required_vs_published(gist)
    gub = gist_upper_bound(gist)
    lys = lysosomal_sizing()
    gsens = gist_sensitivity()
    sc = glioma_scorecard()
    paired = glioma_paired_test(sc)
    loo = per_drug_leave_one_out(sc)
    sens = sensitivity()
    hcc = hcc_downregulation_check()
    out = {"additivity": add, "substrates": subs, "circularity": circ,
           "gist": gist, "gist_required": gistreq, "gist_upper_bound": gub,
           "lysosomal_sizing": lys, "gist_sensitivity": gsens,
           "glioma": sc, "paired": paired, "leave_one_out": loo,
           "sensitivity": sens, "hcc": hcc}
    if write_md:
        out["report_path"] = write_report(add, subs, circ, gist, gistreq, gub,
                                          lys, gsens, sc, paired, loo, sens, hcc)
    return out


if __name__ == "__main__":
    run()
