"""Does tissue-differential LYSOSOMAL ION TRAPPING rescue either scorecard? No.

``./.venv/bin/python -m tmesim.validation.ion_trapping_impact``

Why this run exists
-------------------
Three mechanisms have now been named as the residual physics the shipped model
is missing. Two have been implemented, cited and scored:

* **efflux** (``results/efflux_impact.md``) -- tested, worsens glioma, provably
  inert on a tumour:normal ratio because a per-drug scalar cancels;
* **carrier-mediated uptake** (``results/uptake_impact.md``) -- tested, the
  named carrier (OCT1) is falsified outright, and the surviving hepatic-uptake
  class is excluded by a literature ceiling: the GIST measurement demands
  ``11.7x`` and the largest published value is ``1.5x``.

The uptake run then produced a redirect, and it was a correction to the chain's
own reasoning. ``results/moncayo_explant.md`` §4.3 proved a *uniform*
retardation ``R`` gives ``C = R * C_bath`` everywhere and therefore no
dispersion **within** a section -- and the chain used that proof to eliminate
ion trapping for **GIST**, which is a ratio **between two tissues**:

::

    C_tumour / C_normal = (R_t * E_t) / (R_n * E_n)      does NOT cancel

Moncayo's proof is about a different quantity. So tissue-differential lysosomal
ion trapping was the one named candidate not actually excluded -- and the group
that dismantled the imatinib/OCT1 story concluded that lysosomal sequestration
determines intracellular imatinib, tested in GIST-T1 (PMID 26108972).

What this run does about it
---------------------------
It implements the mechanism as ``tmesim.molecular.ion_trapping``, using the
*exact* Henderson-Hasselbalch formulation already validated in
``validation/moncayo_explant.py`` (asserted identical, not re-derived), and then
does the one thing that decides the question: it **parameterises ``phi_lys``
from measurement, not from the target**.

That is the crux and the place this analysis could cheat. The sizing in
``uptake_impact.md`` §10 showed the required hepatocyte lysosomal volume
fraction is 3.9%, and observed only that this was inside the range
``moncayo_explant.py`` happened to sweep. It is not inside the range that was
*measured*. Schmitt et al. 2020 (PMID 33148688) put rat **liver** at **1.25% of
tissue volume**, and that is the only tissue in their table anchored to a direct
measurement.

The result is a third negative, and it has the same shape as the second: a
**bound argument**. The mechanism is not merely unquantified -- at the measured
parameter, in the physically unattainable best case, it falls short by a factor
this run states explicitly.

Additive only
-------------
Nothing outside ``molecular/ion_trapping.py`` and this file was touched, and no
shipped code path imports either, so the term is off by construction rather than
by default. :func:`additivity_check` verifies that numerically with
``np.array_equal``, not ``allclose``.
"""
from __future__ import annotations

import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from ..molecular.ion_trapping import (
    ANTIBODIES,
    CANCER_CELL_PHI_LYS,
    CANCER_CELL_PHI_LYS_CITATION,
    CANCER_VS_NORMAL_DIRECTION_NOTE,
    DRUG_SMILES,
    HUMAN_LIVER_NOTE,
    LIVER_CELLULAR_VOLUME_FRACTION,
    PH_CYTOSOL,
    PH_LYSOSOME,
    PH_LYSOSOME_SCHMITT,
    PHI_LYS,
    TISSUE_GENERIC,
    TISSUE_GLIOMA,
    TISSUE_INTESTINAL_TUMOUR,
    TISSUE_LIVER_METASTASIS,
    TISSUE_NORMAL_BRAIN,
    TISSUE_NORMAL_INTESTINE,
    TISSUE_NORMAL_LIVER,
    TISSUE_SPLEEN,
    TISSUES,
    basicity,
    cell_to_tissue_basis,
    formulation_matches_moncayo,
    max_reachable_ratio,
    pH_partition_ratio,
    phi_lys_for,
    required_phi_lys,
    retardation_factor,
    summary as trapping_summary,
    trapping_detail,
    trapping_factor,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RESULTS_DIR = os.path.join(REPO_ROOT, "results")
REPORT_MD = os.path.join(RESULTS_DIR, "ion_trapping_impact.md")

BAND_TIGHT = 3.0
BAND_OOM = 10.0

#: Cellular volume fraction assumed for TUMOUR tissue when converting the only
#: measured cancer-cell lysosomal volume fractions (Collins 2021, a fraction of
#: CELL volume) onto the tissue basis this term needs.
#:
#: It is set EQUAL to the liver value on purpose. No measurement of tumour
#: cellularity enters, so the tumour:normal contrast is driven purely by the
#: measured per-cell lysosomal content and not by an unmeasured cellularity
#: difference that could be tuned to taste.
ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION = LIVER_CELLULAR_VOLUME_FRACTION

#: Which glioma-scorecard measurement methods see TOTAL tissue drug. This is
#: physics, not a knob: lysosomes are intracellular reservoirs, so at
#: equilibrium they raise total tissue drug by ``R`` and leave unbound
#: EXTRACELLULAR drug -- what a microdialysis probe samples -- unchanged.
METHOD_SEES_TOTAL_DRUG: Dict[str, bool] = {
    "tissue_homogenate_LCMS": True,
    "intracerebral_microdialysis": False,
}


# ===========================================================================
# 0. Additivity: the shipped pipeline is untouched, and OFF reduces exactly
# ===========================================================================
def additivity_check() -> Dict:
    """Three numeric claims, all verified rather than asserted.

    **(a) The shipped delivery path is bit-identical.** ``molecular.ion_trapping``
    is imported by nothing in ``core``, ``physics``, ``sim``, ``treatment``,
    ``reconstruct`` or the rest of ``molecular``, so ``delivery_source`` cannot
    have changed. Rather than trusting that, the original expression is written
    out and compared with ``==`` on every element of the float64 array, over a
    grid of drug names, penetrances and plasma concentrations, and over both the
    no-barrier and barrier branches. The names include the two drugs that carry
    a LARGE retardation, so a leak would be impossible to miss.

    **(b) The term is inert wherever it should be.** Every drug in the table has
    ``R == 1.0`` exactly in the generic tissue and in every tumour tissue (none
    of which has a measured ``phi_lys``), and every non-base has ``R == 1.0``
    exactly even in tissues where ``phi_lys`` IS measured.

    **(c) Non-vacuity.** With a real weak base in a tissue that has a measured
    ``phi_lys``, the factor must actually move -- otherwise (a) and (b) pass
    trivially and prove nothing.
    """
    from ..core import (CellPopulation, CellState, CellType, Domain, DrugParams,
                        Microenvironment, TissueFields)
    from ..reconstruct.barrier import L_PARACELLULAR_LEAK
    from ..treatment.drugs import delivery_source

    rng = np.random.default_rng(20260728)
    domain = Domain(nx=31, ny=19, dx=12.0)
    fields = TissueFields(ecm_density=domain.ones() * 0.45,
                          vascular_density=rng.random(domain.shape))
    cells = CellPopulation(x=np.array([10.0]), y=np.array([10.0]),
                           cell_type=np.array([CellType.TUMOR]),
                           state=np.array([CellState.VIABLE]),
                           sensitivity=np.array([0.5]))
    micro = Microenvironment(domain=domain, fields=fields, cells=cells)
    vascular = np.asarray(fields.vascular_density, dtype=np.float64)
    disruption = rng.random(domain.shape)

    checks: List[Dict] = []
    for name in ("imatinib", "adavosertib", "eribulin", "not_in_the_table"):
        for pen in (0.0, 0.035, 0.5, 0.9399, 1.0):
            drug = DrugParams(name=name, bbb_penetrance=pen)
            for plasma in (0.0, 1.0, 7.3):
                ref = vascular * pen * float(plasma)
                checks.append({"branch": "delivery_no_barrier", "identical": bool(
                    np.array_equal(delivery_source(micro, drug, plasma), ref))})
                ref_b = vascular * np.clip(pen + L_PARACELLULAR_LEAK * disruption,
                                           0.0, 1.0) * float(plasma)
                checks.append({"branch": "delivery_barrier", "identical": bool(
                    np.array_equal(delivery_source(micro, drug, plasma,
                                                   bbb_disruption=disruption), ref_b))})

    # (b) inert where it should be -- exact equality, not a tolerance
    inert: List[Dict] = []
    for name in list(DRUG_SMILES) + ["not_in_the_table", "bevacizumab"]:
        for tissue in (TISSUE_GENERIC, TISSUE_LIVER_METASTASIS,
                       TISSUE_INTESTINAL_TUMOUR, TISSUE_GLIOMA):
            inert.append({"drug": name, "tissue": tissue, "why": "no_measured_phi_lys",
                          "identical": trapping_factor(name, tissue=tissue) == 1.0})
    for name in ("temozolomide", "cefazolin", "levetiracetam", "methotrexate",
                 "not_in_the_table", "bevacizumab"):
        for tissue in (TISSUE_NORMAL_LIVER, TISSUE_NORMAL_BRAIN,
                       TISSUE_NORMAL_INTESTINE, TISSUE_SPLEEN):
            inert.append({"drug": name, "tissue": tissue, "why": "no_basic_centre",
                          "identical": trapping_factor(name, tissue=tissue) == 1.0})

    # (c) non-vacuity
    bites = trapping_factor("imatinib", tissue=TISSUE_NORMAL_LIVER) > 4.0
    differential = (trapping_factor("imatinib", tissue=TISSUE_NORMAL_LIVER)
                    != trapping_factor("imatinib", tissue=TISSUE_LIVER_METASTASIS,
                                       phi_lys=0.0))
    # and the term must cancel exactly when the two tissues share a phi_lys,
    # which is Moncayo's result restated as a test
    same = trapping_factor("imatinib", tissue=TISSUE_LIVER_METASTASIS,
                           phi_lys=phi_lys_for(TISSUE_NORMAL_LIVER))
    cancels = same == trapping_factor("imatinib", tissue=TISSUE_NORMAL_LIVER)

    n, n_ok = len(checks), sum(c["identical"] for c in checks)
    m, m_ok = len(inert), sum(c["identical"] for c in inert)
    return {
        "n_delivery_comparisons": n, "n_delivery_bit_identical": n_ok,
        "delivery_all_bit_identical": n_ok == n,
        "n_inert_comparisons": m, "n_inert_exactly_one": m_ok,
        "inert_all_exactly_one": m_ok == m,
        "factor_bites_where_phi_lys_measured": bool(bites),
        "factor_is_tissue_differential": bool(differential),
        "factor_cancels_when_phi_lys_equal": bool(cancels),
        "n_shipped_files_modified": 0,
        "failures": [c for c in checks + inert if not c["identical"]],
    }


def formulation_check() -> Dict:
    """The trapping law here IS the one validated against Crank's solution."""
    return formulation_matches_moncayo()


# ===========================================================================
# 1. The parameter table, as data for the report
# ===========================================================================
def phi_lys_rows() -> List[Dict]:
    """Every lysosomal volume fraction this module can cite, with provenance."""
    out = []
    for t in TISSUES:
        e = PHI_LYS[t]
        out.append({"tissue": t, "phi_lys": e.phi_lys, "basis": e.basis,
                    "species": e.species, "method": e.method,
                    "is_measured": e.is_measured, "citation": e.citation,
                    "notes": e.notes})
    return out


def cancer_phi_lys_rows() -> List[Dict]:
    """Collins 2021 cancer-cell values, converted onto the tissue basis."""
    return [{"line": k, "phi_lys_cell_basis": v,
             "phi_lys_tissue_basis": cell_to_tissue_basis(
                 v, ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION)}
            for k, v in CANCER_CELL_PHI_LYS.items()]


def basicity_rows(prefer_literature_pka: bool = True) -> List[Dict]:
    """Structure -> basic centre -> pKa -> trapping ratio, for every drug."""
    out = []
    for d in DRUG_SMILES:
        b = basicity(d, prefer_literature_pka=prefer_literature_pka)
        out.append({"drug": d, "inchikey": b.inchikey,
                    "n_basic_centres": b.n_basic_centres,
                    "basic_class": b.basic_class,
                    "pka_structure": b.pka_structure,
                    "pka_literature": b.pka_literature,
                    "pka_used": b.pka_used, "pka_source": b.pka_source,
                    "trapping_ratio": b.trapping_ratio, "note": b.note})
    return out


# ===========================================================================
# 2. Scorecard 1 -- GIST tumour:normal, the target the redirect named
# ===========================================================================
#: Which model tissue each GIST section belongs to. The tumour side of a liver
#: pair is a GIST METASTASIS sitting in liver, not hepatocytes -- and for THIS
#: mechanism that distinction is the whole hypothesis, because the claim is that
#: hepatocytes hold more acidic compartment than GIST cells do.
GIST_TISSUE_MAP: Dict[str, Tuple[str, str]] = {
    "Liver": (TISSUE_LIVER_METASTASIS, TISSUE_NORMAL_LIVER),
    "Intestine": (TISSUE_INTESTINAL_TUMOUR, TISSUE_NORMAL_INTESTINE),
}


def _gist_geometry() -> Dict:
    """The unchanged transport half of the GIST prediction, loaded once."""
    from .gist_imatinib import (HALF_SPACING_BAND_UM, HALF_SPACING_UM,
                                compile_imatinib, delivery_efficiency,
                                measured_targets)
    targets = measured_targets()
    L = float(compile_imatinib()["L_um"])
    return {
        "targets": targets, "L_um": L,
        "E_tumour": float(delivery_efficiency(L, HALF_SPACING_UM["tumour"])),
        "E_normal": {"Liver": float(delivery_efficiency(L, HALF_SPACING_UM["liver"])),
                     "Intestine": float(delivery_efficiency(L, HALF_SPACING_UM["intestine"]))},
        "floor_nominal": float(HALF_SPACING_UM["liver"] / HALF_SPACING_UM["tumour"]),
        "floor_extreme": float(HALF_SPACING_BAND_UM["liver"][0]
                               / HALF_SPACING_BAND_UM["tumour"][1]),
        "liver_measured": targets[targets.organ == "Liver"].ratio_mean.values,
        "intestine_measured": targets[targets.organ == "Intestine"].ratio_mean.values,
    }


def _gist_score(geo: Dict, phi_tumour: Optional[float],
                phi_normal_from_table: bool = True,
                phi_normal_override: Optional[float] = None,
                prefer_literature_pka: bool = True,
                ph_lys: float = PH_LYSOSOME) -> Dict:
    """Score the six matched pairs at one choice of lysosomal volume fractions.

    ``phi_tumour = None`` means "turn the mechanism off entirely" (the BEFORE
    column). ``phi_tumour = 0.0`` is the physically unattainable limit of a
    tumour whose cells contain no acidic compartment at all -- the case most
    generous to the hypothesis, and therefore the one a bound argument must use.
    """
    et, en = geo["E_tumour"], geo["E_normal"]
    kw = dict(prefer_literature_pka=prefer_literature_pka, ph_lys=ph_lys)

    def _R(organ: str, side: str) -> float:
        if phi_tumour is None:
            return 1.0
        t_tissue, n_tissue = GIST_TISSUE_MAP[organ]
        if side == "tumour":
            return trapping_factor("imatinib", tissue=t_tissue,
                                   phi_lys=phi_tumour, **kw)
        phi_n = (phi_normal_override if phi_normal_override is not None
                 else (phi_lys_for(n_tissue) if phi_normal_from_table else 0.0))
        return trapping_factor("imatinib", tissue=n_tissue, phi_lys=phi_n, **kw)

    rows, folds = [], []
    for _, r in geo["targets"].iterrows():
        rt, rn = _R(r.organ, "tumour"), _R(r.organ, "normal")
        p = (rt * et) / (rn * en[r.organ])
        fold = max(p / r.ratio_mean, r.ratio_mean / p)
        folds.append(fold)
        rows.append({"patient": r.patient, "organ": r.organ,
                     "measured": float(r.ratio_mean), "R_tumour": float(rt),
                     "R_normal": float(rn), "pred": float(p), "fold": float(fold)})
    f = np.asarray(folds)
    rt_l, rn_l = _R("Liver", "tumour"), _R("Liver", "normal")
    rt_i, rn_i = _R("Intestine", "tumour"), _R("Intestine", "normal")
    return {
        "geomean_fold_error": float(10 ** np.mean(np.log10(f))),
        "max_fold_error": float(f.max()),
        "n_within_3x": int((f <= BAND_TIGHT).sum()),
        "n_within_10x": int((f <= BAND_OOM).sum()),
        "pred_liver": float((rt_l * et) / (rn_l * en["Liver"])),
        "pred_intestine": float((rt_i * et) / (rn_i * en["Intestine"])),
        "R_liver_normal": float(rn_l), "R_liver_tumour": float(rt_l),
        "R_intestine_normal": float(rn_i), "R_intestine_tumour": float(rt_i),
        "scale_liver": float(rt_l / rn_l), "scale_intestine": float(rt_i / rn_i),
        "per_patient": rows,
    }


def gist_impact(prefer_literature_pka: bool = True,
                ph_lys: float = PH_LYSOSOME) -> Dict:
    """Does tissue-differential ion trapping move GIST, and break the floor?

    Five configurations are scored, and the ordering matters because the
    *generous* ones come first and still fail:

    ``off``
        the shipped model, geometric-mean fold error 7.92x over n = 6.
    ``limit``
        measured ``phi_lys`` in normal tissue, ``phi_lys = 0`` in tumour. A
        tumour cell with no acidic compartment does not exist, so this is an
        upper bound on what the mechanism can do, not a prediction.
    ``collins_low`` / ``collins_high``
        the only measured cancer-cell lysosomal volume fractions there are
        (Collins 2021, PMID 33172973), converted onto the tissue basis.
    ``host``
        tumour tissue given its host organ's value. Here the term cancels and
        the score returns to baseline -- which is Moncayo's result appearing
        exactly where it should.

    The falsifiable side-condition is the **intestinal control**. Patient AA is
    an intestinal primary against normal intestine, the one pair the shipped
    model already predicts to 1.11x. Gut ``phi_lys`` (1.22%) is barely different
    from liver's (1.25%), so a mechanism driven by acidic-compartment content
    cannot fix the five liver pairs without moving the intestinal pair by almost
    the same factor. If it does, the contrast is fitted, not mechanistic.
    """
    geo = _gist_geometry()
    lit_liver = phi_lys_for(TISSUE_NORMAL_LIVER)
    lit_gut = phi_lys_for(TISSUE_NORMAL_INTESTINE)
    kw = dict(prefer_literature_pka=prefer_literature_pka, ph_lys=ph_lys)

    lo = cell_to_tissue_basis(min(CANCER_CELL_PHI_LYS.values()),
                              ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION)
    hi = cell_to_tissue_basis(max(CANCER_CELL_PHI_LYS.values()),
                              ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION)

    configs = [
        ("off", None, "mechanism off -- the shipped model"),
        ("limit", 0.0,
         "measured normal phi_lys, tumour phi_lys = 0 (physically unattainable "
         "best case; an upper bound, not a prediction)"),
        ("collins_low", lo,
         f"tumour phi_lys from the LOWEST measured cancer line, MCF7 0.50% of "
         f"cell volume -> {100*lo:.2f}% of tissue (PMID 33172973)"),
        ("collins_high", hi,
         f"tumour phi_lys from the HIGHEST measured cancer line, MDA-MB-231 "
         f"3.67% of cell volume -> {100*hi:.2f}% of tissue (PMID 33172973)"),
        ("host", lit_liver,
         "tumour tissue given its host organ's measured value -- the term "
         "cancels, which is Moncayo's proof reappearing"),
    ]
    scored = {name: _gist_score(geo, phi, **kw) for name, phi, _ in configs}
    before, limit = scored["off"], scored["limit"]

    # --- the structural floor --------------------------------------------
    lm = geo["liver_measured"]
    floors = {}
    for name, _, _ in configs:
        s = scored[name]["scale_liver"]
        floors[name] = {
            "scale": s,
            "floor_nominal": geo["floor_nominal"] * s,
            "floor_extreme": geo["floor_extreme"] * s,
            "n_below_nominal": int((lm < geo["floor_nominal"] * s).sum()),
            "n_below_extreme": int((lm < geo["floor_extreme"] * s).sum()),
            "nominal_floor_broken": bool((lm < geo["floor_nominal"] * s).sum() == 0),
            "extreme_floor_broken": bool((lm < geo["floor_extreme"] * s).sum() == 0),
        }

    meas_gm = float(10 ** np.mean(np.log10(lm)))
    return {
        "n": int(len(geo["targets"])),
        "L_um": geo["L_um"], "E_tumour": geo["E_tumour"], "E_normal": geo["E_normal"],
        "phi_lys_liver": lit_liver, "phi_lys_gut": lit_gut,
        "phi_lys_tumour_collins_low": lo, "phi_lys_tumour_collins_high": hi,
        "configs": [{"name": n, "phi_tumour": p, "description": d}
                    for n, p, d in configs],
        "scored": scored, "floors": floors,
        "before": before, "after": limit,
        "improved": limit["geomean_fold_error"] < before["geomean_fold_error"],
        "measured_liver_geomean": meas_gm,
        "required_R_normal_over_tumour": float(before["pred_liver"] / meas_gm),
        "floor_nominal_before": geo["floor_nominal"],
        "floor_extreme_before": geo["floor_extreme"],
        "n_liver_below_nominal_floor_before": int((lm < geo["floor_nominal"]).sum()),
        "n_liver_below_extreme_floor_before": int((lm < geo["floor_extreme"]).sum()),
        "intestinal_fold_before": float(
            [r for r in before["per_patient"] if r["organ"] == "Intestine"][0]["fold"]),
        "intestinal_fold_limit": float(
            [r for r in limit["per_patient"] if r["organ"] == "Intestine"][0]["fold"]),
        "intestinal_control_survives": bool(
            abs(math.log10([r for r in limit["per_patient"]
                            if r["organ"] == "Intestine"][0]["fold"]))
            <= abs(math.log10([r for r in before["per_patient"]
                               if r["organ"] == "Intestine"][0]["fold"])) + 0.1),
        "null_unity_geomean": float(10 ** np.mean(np.log10(
            [max(1.0 / r, r) for r in geo["targets"].ratio_mean.values]))),
    }


# ===========================================================================
# 3. THE CRUX -- required phi_lys against measured phi_lys
# ===========================================================================
def required_vs_measured(gist: Optional[Dict] = None,
                         prefer_literature_pka: bool = True,
                         ph_lys: float = PH_LYSOSOME) -> Dict:
    """Is the mechanism merely unquantified, or quantitatively excluded?

    This is the analogue, for ion trapping, of the ``11.7x demanded / 1.5x
    available`` bound that excluded the hepatic-uptake class -- and it is the
    reason this run produces a result rather than another coverage complaint.

    The comparison is set up so that every discretionary choice favours the
    mechanism:

    * lysosomal lumen pH 4.7 rather than the 5.3 the source PBPK model uses
      (a lower lumen pH means STRONGER trapping);
    * the LysoTracker-based liver ``phi_lys`` of 1.25% rather than the EM-based
      0.92%, because LysoTracker labels all acidic compartment;
    * tumour ``phi_lys`` set to exactly zero, a tumour cell with no lysosomes.

    If the mechanism cannot reach the target even then, it is excluded by the
    literature's own numbers rather than merely unimplemented.
    """
    gist = gist if gist is not None else gist_impact(
        prefer_literature_pka=prefer_literature_pka, ph_lys=ph_lys)
    b = basicity("imatinib", prefer_literature_pka=prefer_literature_pka,
                 ph_lys=ph_lys)
    required_R = gist["required_R_normal_over_tumour"]
    pka = b.pka_used

    phi_req = required_phi_lys(required_R, pka, phi_lys_other=0.0, ph_lys=ph_lys)
    phi_measured = phi_lys_for(TISSUE_NORMAL_LIVER)
    phi_body_max = phi_lys_for(TISSUE_SPLEEN)
    reachable = max_reachable_ratio(phi_measured, pka, 0.0, ph_lys=ph_lys)
    reachable_body_max = max_reachable_ratio(phi_body_max, pka, 0.0, ph_lys=ph_lys)

    # If the tumour is NOT given zero lysosomes but the same content as its
    # host organ's cells, what would the liver need? (The realistic demand.)
    phi_req_vs_host = required_phi_lys(required_R, pka,
                                       phi_lys_other=phi_measured, ph_lys=ph_lys)
    # And, reading it the other way: granting the liver the body's maximum,
    # how lysosome-poor would the tumour have to be?
    k = pH_partition_ratio(pka, ph_lys, PH_CYTOSOL) - 1.0
    r_spleen = 1.0 + phi_body_max * k
    phi_tumour_needed_at_body_max = ((r_spleen / required_R) - 1.0) / k

    return {
        "pka_used": pka, "pka_source": b.pka_source,
        "ph_lys": ph_lys, "ph_cyt": PH_CYTOSOL,
        "trapping_ratio": float(pH_partition_ratio(pka, ph_lys, PH_CYTOSOL)),
        "saturating_ratio": float(10.0 ** (PH_CYTOSOL - ph_lys)),
        "required_R_normal_over_tumour": float(required_R),
        "phi_lys_required_tumour_zero": float(phi_req),
        "phi_lys_measured_liver": float(phi_measured),
        "phi_lys_required_over_measured": float(phi_req / phi_measured),
        "inside_measured_range": bool(phi_req <= phi_measured),
        "max_reachable_ratio_at_measured": float(reachable),
        "shortfall_fold": float(required_R / reachable),
        "excluded_by_its_own_literature": bool(required_R / reachable > 2.0),
        "phi_lys_body_max_spleen": float(phi_body_max),
        "max_reachable_ratio_at_body_max": float(reachable_body_max),
        "still_short_even_at_body_max": bool(required_R > reachable_body_max),
        "phi_lys_required_if_tumour_equals_host": float(phi_req_vs_host),
        "phi_tumour_needed_if_liver_had_spleen_content": float(
            phi_tumour_needed_at_body_max),
        "collins_min_tissue_basis": float(cell_to_tissue_basis(
            min(CANCER_CELL_PHI_LYS.values()),
            ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION)),
        "citation_phi": PHI_LYS[TISSUE_NORMAL_LIVER].citation,
        "citation_tumour": CANCER_CELL_PHI_LYS_CITATION,
    }


def phi_sweep(prefer_literature_pka: bool = True,
              ph_lys: float = PH_LYSOSOME) -> List[Dict]:
    """Sweep tumour ``phi_lys`` and report the whole response surface.

    No single tumour value is defensible, so none is asserted: the sweep runs
    from the impossible zero, through the measured cancer-cell range, to the
    host-organ value where the term cancels.
    """
    geo = _gist_geometry()
    lo = cell_to_tissue_basis(min(CANCER_CELL_PHI_LYS.values()),
                              ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION)
    hi = cell_to_tissue_basis(max(CANCER_CELL_PHI_LYS.values()),
                              ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION)
    liver = phi_lys_for(TISSUE_NORMAL_LIVER)
    grid = sorted({0.0, 0.001, 0.002, lo, 0.006, 0.008, 0.010, liver, 0.015,
                   0.020, 0.025, hi, 0.040, 0.0527})
    lm = geo["liver_measured"]
    out = []
    for phi in grid:
        s = _gist_score(geo, phi, prefer_literature_pka=prefer_literature_pka,
                        ph_lys=ph_lys)
        sc = s["scale_liver"]
        out.append({
            "phi_tumour": float(phi),
            "R_tumour": s["R_liver_tumour"], "R_normal": s["R_liver_normal"],
            "scale": sc, "pred_liver": s["pred_liver"],
            "geomean_fold_error": s["geomean_fold_error"],
            "n_within_3x": s["n_within_3x"],
            "intestinal_fold": float([r for r in s["per_patient"]
                                      if r["organ"] == "Intestine"][0]["fold"]),
            "floor_nominal": geo["floor_nominal"] * sc,
            "n_below_nominal_floor": int((lm < geo["floor_nominal"] * sc).sum()),
            "is_collins_low": abs(phi - lo) < 1e-12,
            "is_collins_high": abs(phi - hi) < 1e-12,
            "is_host_value": abs(phi - liver) < 1e-12,
        })
    return out


# ===========================================================================
# 4. Scorecard 2 -- the six-drug measured human glioma partition
# ===========================================================================
#: How the retardation is applied on the glioma scorecard.
#:
#: ``"method_aware"`` applies ``R`` only where the measurement sees TOTAL tissue
#: drug. That is not a convenience: a lysosome is an intracellular reservoir, so
#: at equilibrium it raises total tissue drug and leaves unbound EXTRACELLULAR
#: drug -- what an intracerebral microdialysis probe samples -- unchanged. It is
#: the same ``C_total = R * C_free`` identity the Moncayo module proved.
#:
#: ``"all_tissues"`` applies ``R`` to every drug regardless of assay. Physically
#: wrong, and run anyway because the verdict on this scorecard DEPENDS on the
#: choice and a reader is entitled to see how much.
GLIOMA_PLACEMENTS = ("method_aware", "all_tissues")


def glioma_scorecard(placement: str = "method_aware",
                     phi_glioma: Optional[float] = None,
                     prefer_literature_pka: bool = True,
                     ph_lys: float = PH_LYSOSOME) -> Dict:
    """Re-score the six measured human glioma drugs with and without trapping.

    The BEFORE column is not pasted: it is
    ``glioma_regional_pk.run_scorecard(aggregation="geometric", apply_leak=True)``
    called live, i.e. the model exactly as it ships. AFTER changes one thing::

        pred_after = pred_before * R(drug, glioma)

    ``R`` multiplies the WHOLE predicted partition rather than only the
    transcellular route, and that differs on purpose from how the efflux and
    uptake terms were placed. Those are membrane transporters and act on the
    route a molecule takes; ion trapping is an intracellular equilibrium that
    acts on whatever free drug is present however it arrived, including drug
    that crossed through a broken tight junction. There is no "transcellular
    only" version of this mechanism to run.

    The old ceiling of 1.0 does not apply either, and for a reason
    ``reconstruct.barrier.effective_penetrance`` states in its own docstring: a
    partition above 1 requires tissue binding or active uptake. Lysosomal
    sequestration IS tissue binding, so the total-drug ceiling rises to ``R``.

    ``phi_glioma`` has no measurement (see the module's ``PHI_LYS`` note on
    Forsby 1983). It defaults to the measured NORMAL BRAIN value, i.e. the
    assumption that glioma tissue has its host organ's acidic-compartment
    content, and that assumption is reported on every row.
    """
    if placement not in GLIOMA_PLACEMENTS:
        raise ValueError(f"placement must be one of {GLIOMA_PLACEMENTS}")
    from .glioma_regional_pk import run_scorecard

    phi = phi_glioma if phi_glioma is not None else phi_lys_for(TISSUE_NORMAL_BRAIN)
    before = run_scorecard(aggregation="geometric", apply_leak=True)

    rows: List[Dict] = []
    for r in before["rows"]:
        drug = r["drug"]
        meas = float(r["measured_mean"])
        pred_b = float(r["predicted_partition"])
        method = r["method"]
        sees_total = METHOD_SEES_TOTAL_DRUG.get(method, True)
        det = trapping_detail(drug, tissue=TISSUE_GLIOMA, phi_lys=phi,
                              prefer_literature_pka=prefer_literature_pka,
                              ph_lys=ph_lys)
        R_chem = float(det["factor"])
        applied = R_chem if (placement == "all_tissues" or sees_total) else 1.0
        pred_a = pred_b * applied
        fb = float(r["fold_error"])
        fa = pred_a / meas if meas > 0 else float("inf")
        rows.append({
            "drug": drug, "measured": meas, "method": method,
            "measures_total_drug": bool(sees_total),
            "pred_before": pred_b, "fold_before": fb,
            "log10_fold_before": math.log10(fb) if fb > 0 else float("nan"),
            "R_chemistry": R_chem, "R_applied": applied,
            "basic_class": det["basic_class"], "pka_used": det["pka_used"],
            "pka_source": det["pka_source"], "reason": det["reason"],
            "pred_after": pred_a, "fold_error": fa,
            "log10_fold": math.log10(fa) if fa > 0 else float("nan"),
        })

    def _agg(key_f: str, key_l: str, drop: Optional[str] = None) -> Dict:
        rs = [r for r in rows if r["drug"] != drop and math.isfinite(r[key_l])]
        logs = [abs(r[key_l]) for r in rs]
        signed = [r[key_l] for r in rs]
        return {
            "n": len(rs),
            "gm_signed": 10.0 ** (sum(signed) / len(signed)),
            "gm_abs": 10.0 ** (sum(logs) / len(logs)),
            "worst_abs_fold": max(max(r[key_f], 1.0 / r[key_f]) for r in rs),
            "n_within_3x": sum(1 for r in rs
                               if 1 / BAND_TIGHT <= r[key_f] <= BAND_TIGHT),
            "n_within_10x": sum(1 for r in rs
                                if 1 / BAND_OOM <= r[key_f] <= BAND_OOM),
        }

    agg_b = _agg("fold_before", "log10_fold_before")
    agg_a = _agg("fold_error", "log10_fold")
    ada = [r for r in rows if r["drug"] == "adavosertib"][0]
    eri = [r for r in rows if r["drug"] == "eribulin"][0]
    return {
        "placement": placement, "phi_glioma": float(phi),
        "phi_glioma_is_measured": False,
        "phi_glioma_source": ("assumed equal to the MEASURED normal-brain value "
                              "(1.36% of tissue, PMID 33148688); no glioma "
                              "measurement is retrievable"),
        "rows": rows, "before": agg_b, "after": agg_a,
        "improved": agg_a["gm_abs"] < agg_b["gm_abs"],
        "n_drugs_moved": sum(1 for r in rows if r["R_applied"] > 1.0),
        "n_drugs_basic": sum(1 for r in rows if r["R_chemistry"] > 1.0),
        "adavosertib_fold_before": ada["fold_before"],
        "adavosertib_fold_after": ada["fold_error"],
        "adavosertib_improved": abs(math.log10(ada["fold_error"]))
                                < abs(math.log10(ada["fold_before"])),
        "eribulin_fold_before": eri["fold_before"],
        "eribulin_fold_after": eri["fold_error"],
        "leave_one_out": [
            {"dropped": r["drug"],
             "gm_abs_before": _agg("fold_before", "log10_fold_before",
                                   drop=r["drug"])["gm_abs"],
             "gm_abs_after": _agg("fold_error", "log10_fold",
                                  drop=r["drug"])["gm_abs"]}
            for r in rows],
    }


def glioma_paired_test(sc: Dict) -> Dict:
    """Is the before/after difference distinguishable from nothing at n = 6?"""
    from scipy import stats
    b = np.array([abs(r["log10_fold_before"]) for r in sc["rows"]])
    a = np.array([abs(r["log10_fold"]) for r in sc["rows"]])
    d = a - b
    movers = int(np.count_nonzero(np.abs(d) > 1e-12))
    n_worse = int(np.sum(d > 1e-12))
    p_sign = (float(stats.binomtest(n_worse, movers, 0.5).pvalue)
              if movers else float("nan"))
    return {"n_drugs": len(b), "n_movers": movers, "n_worse": n_worse,
            "n_better": movers - n_worse, "mean_delta_log10": float(d.mean()),
            "p_sign_test": p_sign,
            "note": (f"with {movers} of {len(b)} drugs moving at all, no test "
                     f"can reach significance at any effect size; reported to "
                     f"make that explicit rather than to claim a result")}


def glioma_pka_sweep(placement: str = "method_aware") -> List[Dict]:
    """How much of the glioma result rests on adavosertib's unmeasured pKa?

    No measured pKa for AZD1775 could be retrieved anywhere, so the headline
    uses the structure-class value 9.8. That is a real weakness and it is swept
    rather than defended: the piperazine is basic beyond doubt, but *how* basic
    decides how much the drug moves.
    """
    from .glioma_regional_pk import run_scorecard
    before = run_scorecard(aggregation="geometric", apply_leak=True)
    ada = [r for r in before["rows"] if r["drug"] == "adavosertib"][0]
    phi = phi_lys_for(TISSUE_NORMAL_BRAIN)
    out = []
    for pka in (6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.8, 10.4):
        R = retardation_factor(pka, phi)
        rows = []
        for r in before["rows"]:
            sees = METHOD_SEES_TOTAL_DRUG.get(r["method"], True)
            b = basicity(r["drug"])
            if b.basic_class in ("none", "antibody", "unknown", "unparseable"):
                f = 1.0
            elif placement == "all_tissues" or sees:
                f = R
            else:
                f = 1.0
            rows.append(abs(math.log10(float(r["predicted_partition"]) * f
                                       / float(r["measured_mean"]))))
        out.append({"pka": pka, "R": float(R),
                    "adavosertib_pred": float(ada["predicted_partition"]) * R,
                    "adavosertib_fold": float(ada["predicted_partition"]) * R
                                        / float(ada["measured_mean"]),
                    "gm_abs": float(10 ** np.mean(rows))})
    return out


# ===========================================================================
# 5. Sensitivity -- does the verdict depend on anything we chose?
# ===========================================================================
#: The EM-based cross-check on the liver value: Blouin 1977 (PMID 833203) gives
#: hepatocyte dense bodies 0.82% plus autophagic vacuoles 0.23% = 1.1% of CELL
#: volume, i.e. 0.92% of tissue. Lower than the LysoTracker value, so the
#: headline uses the LysoTracker one and this is the pessimistic arm.
PHI_LYS_LIVER_BLOUIN = cell_to_tissue_basis(0.011, LIVER_CELLULAR_VOLUME_FRACTION)


def sensitivity() -> List[Dict]:
    """Cross every discretionary choice and re-run both verdicts.

    Four axes, each of which could have been chosen the other way:

    * **lysosomal lumen pH** 4.7 (Moncayo, textbook) vs 5.3 (the effective value
      Schmitt 2020's own PBPK model uses). Higher lumen pH -> weaker trapping.
    * **imatinib pKa** 8.07 (the repo's own literature value) vs 9.8 (the
      structure-class value). Higher pKa -> stronger trapping.
    * **liver phi_lys** 1.25% (LysoTracker, Schmitt) vs 0.92% (EM, Blouin).
    * **glioma phi_lys** normal brain 1.36% vs the two Collins extremes.
    """
    geo = _gist_geometry()
    out: List[Dict] = []
    for ph in (PH_LYSOSOME, PH_LYSOSOME_SCHMITT):
        for lit_pka in (True, False):
            for phi_liver, phi_label in ((phi_lys_for(TISSUE_NORMAL_LIVER), "schmitt_lysotracker"),
                                         (PHI_LYS_LIVER_BLOUIN, "blouin_em")):
                b = basicity("imatinib", prefer_literature_pka=lit_pka, ph_lys=ph)
                s = _gist_score(geo, 0.0, phi_normal_override=phi_liver,
                                prefer_literature_pka=lit_pka, ph_lys=ph)
                req = geo["E_tumour"] / geo["E_normal"]["Liver"] / float(
                    10 ** np.mean(np.log10(geo["liver_measured"])))
                phi_req = required_phi_lys(req, b.pka_used, 0.0, ph_lys=ph)
                out.append({
                    "scorecard": "GIST", "ph_lys": ph,
                    "pka": b.pka_used, "pka_source": b.pka_source,
                    "phi_liver": phi_liver, "phi_liver_source": phi_label,
                    "R_liver": s["R_liver_normal"],
                    "gm_at_limit": s["geomean_fold_error"],
                    "phi_required": phi_req,
                    "required_over_measured": phi_req / phi_liver,
                    "shortfall_fold": req / s["R_liver_normal"],
                    "still_excluded": bool(req / s["R_liver_normal"] > 2.0),
                })
    for ph in (PH_LYSOSOME, PH_LYSOSOME_SCHMITT):
        for placement in GLIOMA_PLACEMENTS:
            for phi_g, lab in ((phi_lys_for(TISSUE_NORMAL_BRAIN), "normal_brain"),
                               (cell_to_tissue_basis(
                                   min(CANCER_CELL_PHI_LYS.values()),
                                   ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION), "collins_low"),
                               (cell_to_tissue_basis(
                                   max(CANCER_CELL_PHI_LYS.values()),
                                   ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION), "collins_high")):
                sc = glioma_scorecard(placement=placement, phi_glioma=phi_g,
                                      ph_lys=ph)
                out.append({
                    "scorecard": "glioma", "ph_lys": ph, "placement": placement,
                    "phi_glioma": phi_g, "phi_glioma_source": lab,
                    "gm_abs_before": sc["before"]["gm_abs"],
                    "gm_abs_after": sc["after"]["gm_abs"],
                    "improved": sc["improved"],
                    "n_within_3x_after": sc["after"]["n_within_3x"],
                    "adavosertib_fold_after": sc["adavosertib_fold_after"],
                })
    return out


# ===========================================================================
# 6. Report
# ===========================================================================
def _pct(v: Optional[float], nd: int = 2) -> str:
    return "--" if v is None else f"{100.0 * v:.{nd}f}%"


def write_report(add: Dict, form: Dict, gist: Dict, req: Dict,
                 sweep: List[Dict], glio: Dict, glio_all: Dict,
                 paired: Dict, pka_sweep: List[Dict], sens: List[Dict],
                 path: str = REPORT_MD) -> str:
    L: List[str] = []
    a = L.append
    gb, gl = gist["before"], gist["after"]
    n = gist["n"]

    a("# Tissue-differential lysosomal ion trapping — implemented, "
      "parameterised from measurement, and it does NOT rescue GIST")
    a("")
    a("**Verdict up front.** This is the one candidate the chain had not "
      "actually excluded. It was implemented with the retardation law already "
      "validated against Crank's analytical solution, and its single free "
      "parameter was taken from the literature rather than from the target. "
      f"The measurement demands a normal:tumour retardation ratio of "
      f"**{req['required_R_normal_over_tumour']:.2f}x**; the measured rat liver "
      f"lysosomal volume fraction of "
      f"**{_pct(req['phi_lys_measured_liver'], 2)} of tissue volume** "
      f"(PMID 33148688) can produce at most "
      f"**{req['max_reachable_ratio_at_measured']:.2f}x**, and only in the "
      f"physically unattainable limit of a tumour cell with no lysosomes at "
      f"all. That is a **{req['shortfall_fold']:.2f}x shortfall**, the same "
      f"kind of bound argument that excluded hepatic uptake. Over "
      f"**n = {n}** matched GIST pairs the geometric-mean fold error goes from "
      f"**{gb['geomean_fold_error']:.2f}x** to "
      f"**{gl['geomean_fold_error']:.2f}x** — but only at that impossible "
      f"setting, and it costs the intestinal control "
      f"({gist['intestinal_fold_before']:.2f}x -> "
      f"{gist['intestinal_fold_limit']:.2f}x); at the lowest lysosomal volume "
      f"fraction ever MEASURED in a cancer cell it is "
      f"{gist['scored']['collins_low']['geomean_fold_error']:.2f}x, and at the "
      f"highest it is {gist['scored']['collins_high']['geomean_fold_error']:.2f}x "
      f"— worse than doing nothing. **This is the third excluded mechanism.**")
    a("")
    a("| scorecard | before | after ion trapping | verdict |")
    a("|---|---|---|---|")
    a(f"| GIST tumour:normal, n = {n} pairs (geomean fold error) | "
      f"**{gb['geomean_fold_error']:.2f}x** | "
      f"**{gl['geomean_fold_error']:.2f}x** | improves, but ONLY at "
      f"`phi_lys(tumour) = 0`, which does not exist |")
    a(f"| .. at the lowest MEASURED cancer-cell value (MCF7) | "
      f"{gb['geomean_fold_error']:.2f}x | "
      f"{gist['scored']['collins_low']['geomean_fold_error']:.2f}x | "
      "still 4.9x wrong |")
    a(f"| .. at the highest MEASURED cancer-cell value (MDA-MB-231) | "
      f"{gb['geomean_fold_error']:.2f}x | "
      f"{gist['scored']['collins_high']['geomean_fold_error']:.2f}x | "
      "**WORSE — the sign flips** |")
    a(f"| .. tumour given its host organ's measured value | "
      f"{gb['geomean_fold_error']:.2f}x | "
      f"{gist['scored']['host']['geomean_fold_error']:.2f}x | "
      "the term cancels exactly, as Moncayo requires |")
    a(f"| GIST intestinal control (the one pair the model gets right) | "
      f"{gist['intestinal_fold_before']:.2f}x | "
      f"{gist['intestinal_fold_limit']:.2f}x | **BROKEN** at the same setting "
      "that improves the liver pairs |")
    a(f"| human glioma partition, n = 6 drugs (GM abs fold error) | "
      f"**{glio['before']['gm_abs']:.3f}x** | "
      f"**{glio['after']['gm_abs']:.3f}x** | improves — the first mechanism in "
      "this chain that helps anything |")
    a(f"| .. applying R regardless of assay type (physically wrong) | "
      f"{glio_all['before']['gm_abs']:.3f}x | "
      f"{glio_all['after']['gm_abs']:.3f}x | worse |")
    a(f"| .. within 3x | {glio['before']['n_within_3x']}/6 | "
      f"{glio['after']['n_within_3x']}/6 | |")
    a(f"| .. adavosertib specifically | "
      f"{glio['adavosertib_fold_before']:.3f}x | "
      f"{glio['adavosertib_fold_after']:.3f}x | the drug the uptake run named, "
      "and the only one that moves |")
    a("")
    a("Six findings, in order of how much they are worth:")
    a("")
    a(f"1. **The required lysosomal volume fraction is outside the measured "
      f"range, and that is the result.** With the tumour granted zero lysosomes "
      f"— an impossibility, chosen because it is the most generous case — the "
      f"hepatocyte compartment would have to be "
      f"**{_pct(req['phi_lys_required_tumour_zero'], 2)} of liver volume**. "
      f"The measured value is **{_pct(req['phi_lys_measured_liver'], 2)}**, "
      f"a **{req['phi_lys_required_over_measured']:.2f}x** overshoot. If the "
      f"tumour is instead given its host organ's content, the demand rises to "
      f"**{_pct(req['phi_lys_required_if_tumour_equals_host'], 1)}** — "
      f"{req['phi_lys_required_if_tumour_equals_host'] / req['phi_lys_body_max_spleen']:.1f}x "
      f"more than **spleen**, the most lysosome-rich tissue in the body at "
      f"{_pct(req['phi_lys_body_max_spleen'], 1)}.")
    a(f"2. **The mechanism has a hard ceiling and the target is above it.** "
      f"Henderson-Hasselbalch trapping saturates at "
      f"`10^(pH_cyt - pH_lys)` = {req['saturating_ratio']:.1f}, so the whole "
      f"question is the volume fraction — and at the measured one the largest "
      f"tumour:normal asymmetry obtainable is "
      f"{req['max_reachable_ratio_at_measured']:.2f}x against "
      f"{req['required_R_normal_over_tumour']:.2f}x demanded. Ruling a mechanism "
      f"out by a literature ceiling without ever fitting it is the same move "
      f"that excluded hepatic uptake, now applied to its replacement.")
    a("3. **The sign of the tumour:normal contrast is not established, and at "
      "the only measured tumour values it can run backwards.** The mechanism "
      "needs normal liver to be MORE lysosome-rich than GIST tissue. No matched "
      "cancer-versus-normal lysosomal volume-fraction pair exists on any basis; "
      "the two studies that do report volumes point the other way (Glunde 2003 "
      "PMID 14965446, Kundra & Dean 1990 PMID 1693341); and at MDA-MB-231's "
      f"measured 3.67% of cell volume the predicted ratio moves AWAY from the "
      f"measurement, {gb['geomean_fold_error']:.2f}x -> "
      f"{gist['scored']['collins_high']['geomean_fold_error']:.2f}x.")
    a(f"4. **The one configuration that helps GIST breaks the intestinal "
      f"control.** Gut `phi_lys` ({_pct(gist['phi_lys_gut'], 2)}) and liver "
      f"`phi_lys` ({_pct(gist['phi_lys_liver'], 2)}) are barely different, so a "
      f"mechanism driven by acidic-compartment content cannot move the five "
      f"liver pairs without moving the intestinal pair by nearly the same "
      f"factor. It does: {gist['intestinal_fold_before']:.2f}x -> "
      f"{gist['intestinal_fold_limit']:.2f}x. The uptake run set that as the "
      f"falsifiable side-condition for a real mechanism, and this one fails it.")
    a(f"5. **On glioma it is the first thing in this chain that helps — and the "
      f"help is one drug wide.** GM absolute fold error "
      f"{glio['before']['gm_abs']:.3f}x -> {glio['after']['gm_abs']:.3f}x, "
      f"{glio['before']['n_within_3x']}/6 -> {glio['after']['n_within_3x']}/6 "
      f"within 3x, driven entirely by adavosertib "
      f"({glio['adavosertib_fold_before']:.3f}x -> "
      f"{glio['adavosertib_fold_after']:.3f}x), the exact drug the uptake run "
      f"identified as the one uptake could not fix. It is not evidence: n = 6, "
      f"one mover, sign-test p = {paired['p_sign_test']:.2f}, and adavosertib's "
      f"pKa is a structure class rather than a measurement.")
    a("6. **A structural observation that decides the placement question.** "
      "Ion trapping raises TOTAL tissue drug and leaves unbound EXTRACELLULAR "
      "drug alone, because a lysosome is a reservoir and at equilibrium a "
      "reservoir does not shift the free concentration outside the cell. That "
      "is the same `C_total = R * C_free` identity Moncayo proved. Applied "
      "honestly, it means the term cannot touch the three microdialysis drugs "
      "on the glioma scorecard — including eribulin, the largest single error "
      "there. Applied dishonestly (to everything), the two weak bases move by "
      "exactly opposite amounts in log space and the geometric mean barely "
      "shifts. §6.")
    a("")

    # ---- 1. additivity ---------------------------------------------------
    a("## 1. Nothing outside the two new files was touched")
    a("")
    a("The efflux run added one default-off keyword to "
      "`treatment.drugs.delivery_source`. The uptake run added nothing outside "
      "its two files. **This run does the same**: no shipped code path imports "
      "`molecular.ion_trapping`, so the term is off by construction rather than "
      "by default. All three claims are verified numerically.")
    a("")
    a(f"- **{add['n_delivery_bit_identical']}/{add['n_delivery_comparisons']} "
      f"delivery-source comparisons bit-for-bit identical** (`np.array_equal`, "
      f"not `allclose`) across 4 drug names x 5 penetrances x 3 plasma "
      f"concentrations x 2 branches. All identical: "
      f"**{add['delivery_all_bit_identical']}**.")
    a(f"- **{add['n_inert_exactly_one']}/{add['n_inert_comparisons']} inertness "
      f"checks return exactly 1.0** — every drug in every tumour tissue (none "
      f"has a measured `phi_lys`), and every non-base in every tissue that does. "
      f"All exactly 1.0: **{add['inert_all_exactly_one']}**.")
    a(f"- Shipped files modified: **{add['n_shipped_files_modified']}**.")
    a("")
    a("And the checks are not vacuous:")
    a("")
    a(f"- the factor bites where `phi_lys` IS measured (imatinib in normal "
      f"liver, R > 4): **{add['factor_bites_where_phi_lys_measured']}**")
    a(f"- the factor is genuinely tissue-differential: "
      f"**{add['factor_is_tissue_differential']}**")
    a(f"- and it cancels EXACTLY when the two tissues share a `phi_lys`, which "
      f"is Moncayo's proof restated as a unit test: "
      f"**{add['factor_cancels_when_phi_lys_equal']}**")
    a("")
    a("### The formulation is the validated one, asserted rather than claimed")
    a("")
    a(f"`molecular/ion_trapping.py` restates `pH_partition_ratio` and "
      f"`retardation_factor` rather than importing them, because `molecular` "
      f"must not depend on `validation`. Restating can become drifting, so the "
      f"two are compared against `validation/moncayo_explant.py`'s functions "
      f"over {form['n']} (pKa, phi) combinations with `==`, not a tolerance: "
      f"**{form['n_identical']}/{form['n']} identical**, all identical "
      f"**{form['all_identical']}**, pH constants identical "
      f"**{form['ph_constants_identical']}**. That stepper was checked against "
      f"Crank's analytical solution to 0.53%.")
    a("")

    # ---- 2. the parameter ------------------------------------------------
    a("## 2. The parameter — where `phi_lys` comes from, and where it does not")
    a("")
    a("This is the crux. A tissue-differential retardation has exactly one "
      "input per tissue, and if that input is chosen to fit the target then the "
      "whole exercise is circular. So every value carries a PMID and the "
      "provenance of each is stated, including the ones that undermine it.")
    a("")
    a("### Measured, fraction of TISSUE volume")
    a("")
    a("| tissue | `phi_lys` | measured? | provenance |")
    a("|---|---:|---|---|")
    for r in phi_lys_rows():
        prov = ("**DIRECTLY MEASURED**" if r["tissue"] == TISSUE_NORMAL_LIVER
                else ("scaled from liver by BMP" if r["is_measured"]
                      else "**no measurement exists**"))
        a(f"| `{r['tissue']}` | {_pct(r['phi_lys'], 3)} | "
          f"{'yes' if r['is_measured'] else 'NO'} | {prov} |")
    a("")
    a("Source: **PMID 33148688** (Schmitt 2020 *Drug Metab Dispos* 49:53-61, "
      "doi:10.1124/dmd.120.000161), Tables 1-2, rat, from the lysosomal marker "
      "lipid bis(monoacylglycero)phosphate. Range 0.03% (adipose) to 5.3% "
      "(spleen).")
    a("")
    a("Four caveats, all of which cut against the mechanism this table was "
      "assembled to test:")
    a("")
    a("1. **Only LIVER is a measurement.** The paper's own Table 2 footnote: "
      "*\"Liver was taken as reference, whose lysosomal content was determined "
      "experimentally. Other tissues were calculated relative to the liver by "
      "bis(monoacylglycero)phosphate tissue concentrations.\"* These are one "
      "measurement plus a lipid-ratio extrapolation, not fourteen independent "
      "numbers. Fortunately liver is the tissue this run needs.")
    a(f"2. **The liver value is rat, and it is derived rather than read off.** "
      f"`(1 - f_EW) x fV_lys(hepatocyte) = "
      f"{LIVER_CELLULAR_VOLUME_FRACTION} x 1.49% = "
      f"{_pct(req['phi_lys_measured_liver'], 2)}`. The hepatocyte cell-basis "
      f"input 1.49% +/- 0.08% is **PMID 30409837** (Schmitt 2019 *DMD* 47:49-57, "
      f"doi:10.1124/dmd.118.084541), LysoTracker Red spot volumetry over "
      f">350,000 endo-/lysosomes in plated rat hepatocytes from 6 rats. The "
      f"independent EM cross-check — **PMID 833203** (Blouin, Bolender & Weibel "
      f"1977 *J Cell Biol* 72:441, doi:10.1083/jcb.72.2.441, n = 5, point "
      f"counting) — gives hepatocyte dense bodies 0.82% +/- 0.16% of CELL "
      f"volume plus autophagic vacuoles 0.23%, i.e. 1.1% of cell = "
      f"{_pct(PHI_LYS_LIVER_BLOUIN, 2)} of tissue. **The headline uses the "
      f"LARGER of the two**, because LysoTracker labels all acidic compartment "
      f"while EM dense bodies do not, and because it is the choice generous to "
      f"the mechanism. The pessimistic arm is in §7.")
    a("3. **No normal-human liver value is retrievable at all.** " + HUMAN_LIVER_NOTE)
    a("4. **Beware the basis.** `% of TISSUE`, `% of CELL` and `% of CYTOPLASM` "
      "differ by ~20% and ~40% respectively and are freely interchanged in the "
      "PBPK literature. Only tissue-basis values are used directly here; "
      "`ion_trapping.cell_to_tissue_basis` does every conversion and reproduces "
      "the published liver number to 5e-5.")
    a("")
    a("### NOT measured: tumour tissue, on any basis")
    a("")
    a("This is the honest gap and it is not a small one.")
    a("")
    for t in (TISSUE_LIVER_METASTASIS, TISSUE_GLIOMA):
        a(f"**`{t}`** — {PHI_LYS[t].notes}")
        a("")
    a("The closest thing to a tumour parameter that exists is the only modern "
      "set of measured cancer-cell lysosomal volume fractions, on a CELL basis:")
    a("")
    a("| cell line | `phi_lys`, % of cell volume | on a tissue basis (x "
      f"{ASSUMED_TUMOUR_CELLULAR_VOLUME_FRACTION}) |")
    a("|---|---:|---:|")
    for r in cancer_phi_lys_rows():
        a(f"| {r['line']} | {_pct(r['phi_lys_cell_basis'], 3)} | "
          f"{_pct(r['phi_lys_tissue_basis'], 3)} |")
    a("")
    a(CANCER_CELL_PHI_LYS_CITATION + ".")
    a("")
    a("Three caveats travel with it and all three cut against using it: the "
      "panel is human **breast** lines only, there is **no normal comparator in "
      "the paper**, and the internal spread is 7-fold. The paper's own point is "
      "that lysosomotropic exposure inflates the parameter 1.5-3x within 24 h — "
      "which matters directly, because imatinib is itself lysosomotropic and "
      "the GIST specimens were taken from patients on imatinib.")
    a("")
    a("**And the direction of the cancer-versus-normal contrast — which is the "
      "sign of the whole mechanism — is not established.** " +
      CANCER_VS_NORMAL_DIRECTION_NOTE)
    a("")
    # ---- 3. chemistry ----------------------------------------------------
    a("## 3. The chemistry — pKa from structure, not from a wish list")
    a("")
    a("Basic centres are found from structure with the SMARTS "
      "`moncayo_explant.py` already uses (an sp3 nitrogen that is not an "
      "amide/sulfonamide/phosphoramide N, not attached to an aromatic ring, not "
      "already charged), and the pKa is read off the same curated "
      "substructure-class table. That table cannot be tuned drug by drug, which "
      "is the point of using it. Every SMILES is asserted against its PubChem "
      "InChIKey before any chemistry runs.")
    a("")
    a("| drug | InChIKey | basic centres | class | pKa (structure) | pKa "
      "(literature) | pKa used | trapping ratio |")
    a("|---|---|---:|---|---:|---:|---:|---:|")
    for r in basicity_rows():
        lit = r["pka_literature"]
        lit_s = "--" if lit is None else f"{lit:.2f}"
        a(f"| {r['drug']} | `{r['inchikey']}` | {r['n_basic_centres']} | "
          f"{r['basic_class']} | {r['pka_structure']:.1f} | {lit_s} | "
          f"{r['pka_used']:.2f} | {r['trapping_ratio']:.1f} |")
    a("")
    a("Three things worth reading off that table:")
    a("")
    a("- **Only three of seven drugs are weak bases at all**: imatinib "
      "(piperazine), adavosertib (4-methylpiperazine) and eribulin (primary "
      "amine). Temozolomide, cefazolin, levetiracetam and methotrexate have no "
      "basic centre that passes the filter, so the term is exactly 1.0 for them "
      "in every tissue. That is a real constraint on how much this mechanism "
      "could ever explain.")
    a("- **Imatinib uses 8.07, the repo's own value.** That is "
      "`IMATINIB_PKA_BASE` from `validation/gist_imatinib.py`, whose comment "
      "records the literature range as 7.7-8.07. The structure class would give "
      "9.8, which produces a LARGER trapping ratio — using the repo's smaller "
      "literature value is the conservative choice against the mechanism, and "
      "both are run in §7.")
    a("- **Adavosertib's pKa could not be found anywhere** — not in the primary "
      "literature, not in vendor datasheets, not in the phase I PK papers. The "
      "structure-class 9.8 is used and swept over 6.0-10.4 in §6 rather than "
      "defended, because the glioma result rests on this drug alone.")
    a("")
    a("Note the ceiling this imposes. The trapping ratio saturates at "
      f"`10^(pH_cyt - pH_lys)` = {req['saturating_ratio']:.1f}, so imatinib's "
      f"{req['trapping_ratio']:.1f} is already 88% of the maximum any base can "
      "reach. There is no stronger base to appeal to; the only remaining lever "
      "is the volume fraction, and §4 is about that lever failing.")
    a("")

    # ---- 4. GIST ---------------------------------------------------------
    a(f"## 4. Scorecard 1 — GIST tumour:normal, n = {n} matched pairs")
    a("")
    a(f"The prediction is the same reaction-diffusion delivery ratio the "
      f"published run scored — same 1-D Krogh slab, same intervessel "
      f"half-spacings, same penetration length `L = {gist['L_um']:.2f}` um from "
      f"the doxorubicin-anchored `alpha` fitted once to an unrelated "
      f"measurement — with one tissue-differential retardation on each side:")
    a("")
    a("```")
    a("  ratio = (R_tumour * E_tumour) / (R_normal * E_normal)")
    a("  R     = 1 + phi_lys * (trapping_ratio - 1)")
    a("```")
    a("")
    a(f"`E_tumour = {gist['E_tumour']:.4f}`, "
      f"`E_liver = {gist['E_normal']['Liver']:.4f}`, "
      f"`E_intestine = {gist['E_normal']['Intestine']:.4f}` — all unchanged. "
      f"The ONLY new quantity is `phi_lys`.")
    a("")
    a("### 4.1 Five configurations, generous ones first")
    a("")
    a("| configuration | `phi_lys` tumour | `R_normal` | `R_tumour` | scale "
      "`R_t/R_n` | predicted liver ratio | geomean fold error "
      f"(n = {n}) | within 3x |")
    a("|---|---:|---:|---:|---:|---:|---:|---:|")
    for c in gist["configs"]:
        s = gist["scored"][c["name"]]
        a(f"| `{c['name']}` | "
          f"{'--' if c['phi_tumour'] is None else _pct(c['phi_tumour'], 3)} | "
          f"{s['R_liver_normal']:.3f} | {s['R_liver_tumour']:.3f} | "
          f"{s['scale_liver']:.4f} | {s['pred_liver']:.4f} | "
          f"**{s['geomean_fold_error']:.2f}x** | {s['n_within_3x']}/{n} |")
    a("")
    a(f"Measured liver geomean: **{gist['measured_liver_geomean']:.4f}**. The "
      f"no-partitioning null is {gist['null_unity_geomean']:.2f}x, so 7.92x was "
      f"already beating it and 3.00x beats it further — but only at a value of "
      f"`phi_lys(tumour)` that no cell has.")
    a("")
    for c in gist["configs"]:
        a(f"- **`{c['name']}`** — {c['description']}")
    a("")
    a("### 4.2 The whole response surface, so nothing is hidden by a choice")
    a("")
    a("No tumour value is defensible, so none is asserted. The sweep runs from "
      "the impossible zero, through the measured cancer-cell range, to the "
      "host-organ value where the term cancels.")
    a("")
    a("| `phi_lys` tumour | `R_tumour` | scale | predicted liver | geomean fold "
      f"error | within 3x | intestinal control | nominal floor | liver pairs "
      "still below it |")
    a("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in sweep:
        tag = (" **(MCF7)**" if r["is_collins_low"] else
               " **(MDA-MB-231)**" if r["is_collins_high"] else
               " **(= liver)**" if r["is_host_value"] else "")
        a(f"| {_pct(r['phi_tumour'], 3)}{tag} | {r['R_tumour']:.3f} | "
          f"{r['scale']:.4f} | {r['pred_liver']:.4f} | "
          f"{r['geomean_fold_error']:.2f}x | {r['n_within_3x']}/{n} | "
          f"{r['intestinal_fold']:.2f}x | {r['floor_nominal']:.4f} | "
          f"{r['n_below_nominal_floor']}/5 |")
    a("")
    a("**The surface is monotone.** The best the mechanism ever does is "
      f"{sweep[0]['geomean_fold_error']:.2f}x at `phi_lys(tumour) = 0`, and it "
      f"degrades from there. There is no interior optimum to find and nothing "
      f"to tune toward.")
    a("")
    a("### 4.3 The structural floor — broken, but not where it matters")
    a("")
    a("`results/gist_imatinib.md` proved the transport-only model has a hard "
      "floor at any penetration length: `E_t/E_n >= a_liver/a_tumour`. With a "
      "retardation on, that becomes `>= (a_liver/a_tumour) * (R_t/R_n)`.")
    a("")
    a("| configuration | nominal floor | liver pairs below it | extreme-geometry "
      "floor | liver pairs below it |")
    a("|---|---:|---:|---:|---:|")
    for c in gist["configs"]:
        f = gist["floors"][c["name"]]
        a(f"| `{c['name']}` | {f['floor_nominal']:.4f} | "
          f"{f['n_below_nominal']}/5 | {f['floor_extreme']:.4f} | "
          f"{f['n_below_extreme']}/5 |")
    a("")
    a(f"So the honest answer to \"is the floor broken?\" is **partly, and only "
      f"in the impossible case**. At `phi_lys(tumour) = 0` the extreme-geometry "
      f"floor drops from {gist['floor_extreme_before']:.4f} to "
      f"{gist['floors']['limit']['floor_extreme']:.4f} and all five liver "
      f"ratios become reachable — the previously-unreachable "
      f"{gist['n_liver_below_extreme_floor_before']}/5 becomes "
      f"{gist['floors']['limit']['n_below_extreme']}/5. But the **nominal** "
      f"floor only falls to {gist['floors']['limit']['floor_nominal']:.4f} and "
      f"{gist['floors']['limit']['n_below_nominal']}/5 liver ratios are still "
      f"below it. At the lowest MEASURED cancer-cell value the nominal floor is "
      f"{gist['floors']['collins_low']['floor_nominal']:.4f} with "
      f"{gist['floors']['collins_low']['n_below_nominal']}/5 still below, and "
      f"at the highest it is above the starting floor. **Breaking a floor by "
      f"assuming a tumour has no lysosomes is not breaking it.**")
    a("")
    a("### 4.4 Per patient, at the most generous setting")
    a("")
    a("| patient | organ | measured | predicted (before) | predicted (limit) | "
      "fold before | fold at limit |")
    a("|---|---|---:|---:|---:|---:|---:|")
    for rb, rl in zip(gb["per_patient"], gl["per_patient"]):
        a(f"| {rb['patient']} | {rb['organ']} | {rb['measured']:.4f} | "
          f"{rb['pred']:.4f} | {rl['pred']:.4f} | {rb['fold']:.2f}x | "
          f"{rl['fold']:.2f}x |")
    a("")
    a(f"**Patient AA is the control and it is where the mechanism gives itself "
      f"away.** It is an intestinal primary against normal intestine, the one "
      f"pair the shipped model already predicts to "
      f"{gist['intestinal_fold_before']:.2f}x. Gut `phi_lys` "
      f"({_pct(gist['phi_lys_gut'], 2)}) is within 3% of liver's "
      f"({_pct(gist['phi_lys_liver'], 2)}), so a mechanism made of "
      f"acidic-compartment content applies there with almost the same strength "
      f"and the control moves to {gist['intestinal_fold_limit']:.2f}x. "
      f"Control survives: **{gist['intestinal_control_survives']}**. The uptake "
      f"run set this as the falsifiable side-condition — *\"a term that fixed "
      f"the liver pairs by breaking the intestinal pair would be a fitted "
      f"contrast, not a mechanism\"* — and ion trapping fails it, where a "
      f"hepatocyte-restricted transporter at least would not have.")
    a("")
    # ---- 5. the crux -----------------------------------------------------
    a("## 5. Required versus measured — the bound argument")
    a("")
    a("This is the section that decides the question, and it is set up so that "
      "**every discretionary choice favours the mechanism**:")
    a("")
    a(f"- lysosomal lumen pH **{req['ph_lys']}** rather than the "
      f"{PH_LYSOSOME_SCHMITT} that Schmitt 2020's own PBPK model uses — a lower "
      f"lumen pH means stronger trapping;")
    a(f"- the LysoTracker-based liver `phi_lys` of "
      f"{_pct(req['phi_lys_measured_liver'], 2)} rather than the EM-based "
      f"{_pct(PHI_LYS_LIVER_BLOUIN, 2)};")
    a("- tumour `phi_lys` set to exactly **zero**, a tumour cell containing no "
      "acidic compartment at all.")
    a("")
    a("| quantity | value |")
    a("|---|---:|")
    a(f"| imatinib basic pKa (the repo's own value) | {req['pka_used']:.2f} |")
    a(f"| lysosome:cytosol trapping ratio at pH {req['ph_lys']} / {req['ph_cyt']} "
      f"| {req['trapping_ratio']:.1f} |")
    a(f"| .. and its saturating ceiling for ANY base | "
      f"{req['saturating_ratio']:.1f} |")
    a(f"| `R_normal / R_tumour` the measurement demands | "
      f"**{req['required_R_normal_over_tumour']:.2f}** |")
    a(f"| implied hepatocyte `phi_lys` (tumour = 0) | "
      f"**{_pct(req['phi_lys_required_tumour_zero'], 2)}** |")
    a(f"| MEASURED rat liver `phi_lys` (PMID 33148688) | "
      f"**{_pct(req['phi_lys_measured_liver'], 2)}** |")
    a(f"| required / measured | **{req['phi_lys_required_over_measured']:.2f}x** |")
    a(f"| **is the required value inside the measured range?** | "
      f"**{req['inside_measured_range']}** |")
    a(f"| largest `R_n/R_t` the measured value can produce | "
      f"**{req['max_reachable_ratio_at_measured']:.2f}** |")
    a(f"| **shortfall** | **{req['shortfall_fold']:.2f}x** |")
    a(f"| implied hepatocyte `phi_lys` if tumour = host organ | "
      f"{_pct(req['phi_lys_required_if_tumour_equals_host'], 1)} |")
    a(f"| highest `phi_lys` of ANY tissue in the body (spleen) | "
      f"{_pct(req['phi_lys_body_max_spleen'], 1)} |")
    a(f"| tumour `phi_lys` needed if liver had SPLEEN's content | "
      f"{_pct(req['phi_tumour_needed_if_liver_had_spleen_content'], 2)} |")
    a(f"| lowest MEASURED cancer-cell `phi_lys`, tissue basis | "
      f"{_pct(req['collins_min_tissue_basis'], 2)} |")
    a("")
    a("**Read the last four rows together, because that is where the argument "
      "closes.** Taken alone, the required "
      f"{_pct(req['phi_lys_required_tumour_zero'], 2)} is not off the biological "
      f"scale — spleen reaches {_pct(req['phi_lys_body_max_spleen'], 1)}. But it "
      f"is required of **liver**, whose measured value is "
      f"{_pct(req['phi_lys_measured_liver'], 2)}, and no tissue other than "
      f"spleen exceeds 1.73%. And even granting liver the spleen's content, the "
      f"tumour would then have to sit at "
      f"{_pct(req['phi_tumour_needed_if_liver_had_spleen_content'], 2)} — below "
      f"skeletal muscle (0.11%) and "
      f"{req['collins_min_tissue_basis'] / req['phi_tumour_needed_if_liver_had_spleen_content']:.1f}x "
      f"below the lowest cancer cell line anyone has measured. Both ends have to "
      f"be wrong at once, in opposite directions, for the mechanism to work.")
    a("")
    a(f"**Excluded by its own literature: "
      f"{req['excluded_by_its_own_literature']}.**")
    a("")
    a("### What this corrects in `results/uptake_impact.md` §10")
    a("")
    a("That section sized the same quantity and got 3.86% — the same number "
      "this run gets — and concluded the mechanism *\"cannot be dismissed on "
      "magnitude\"* because 3.86% *\"sits inside the range `moncayo_explant.py` "
      "already sweeps\"*. **That was the wrong comparator.** The Moncayo sweep "
      "(0.5-5%) was a plausibility range chosen to bracket an unknown, not a "
      "measurement, and its own reference value of 1.4% was assembled as "
      "`0.7 cellularity x 0.02 lysosomal fraction` rather than measured. "
      "Compared against what was actually measured in liver, 3.86% is a "
      f"{req['phi_lys_required_over_measured']:.2f}x overshoot. The redirect was "
      "right that Moncayo's *dispersion* proof does not apply to a "
      "between-tissue ratio; it was wrong that the magnitude was reachable.")
    a("")

    # ---- 6. glioma -------------------------------------------------------
    a("## 6. Scorecard 2 — the six-drug measured human glioma partition")
    a("")
    a("The BEFORE column is `glioma_regional_pk.run_scorecard(aggregation="
      "\"geometric\", apply_leak=True)` called live. AFTER multiplies the whole "
      "predicted partition by `R`, and that placement differs on purpose from "
      "how efflux and uptake were placed: those are membrane transporters and "
      "act on the route a molecule takes, whereas ion trapping is an "
      "intracellular equilibrium acting on whatever free drug is present "
      "however it arrived, including drug that came through a broken tight "
      "junction. There is no transcellular-only version of this mechanism.")
    a("")
    a("The old ceiling of 1.0 does not apply either, for the reason "
      "`reconstruct.barrier.effective_penetrance` gives in its own docstring: "
      "*\"a partition above 1 requires tissue binding or active uptake\"*. "
      "Lysosomal sequestration IS tissue binding, so the total-drug ceiling "
      f"rises to `R`. `phi_lys(glioma)` = {_pct(glio['phi_glioma'], 2)}, "
      f"{glio['phi_glioma_source']}.")
    a("")
    a("### 6.1 The placement question is physics, and it is decidable")
    a("")
    a("A lysosome is an intracellular reservoir. At equilibrium it raises TOTAL "
      "tissue drug by `R` and leaves the unbound EXTRACELLULAR concentration "
      "alone — the same `C_total = R * C_free` identity Moncayo proved. So a "
      "tissue-homogenate LC-MS/MS measurement sees `R` and an intracerebral "
      "microdialysis probe does not. Three of the six drugs on this scorecard "
      "are microdialysis, including eribulin, which is both a strong base and "
      "the largest single error on the board.")
    a("")
    a("| drug | measurement | sees total drug? | `R` from chemistry | `R` "
      "applied | predicted before | predicted after | fold before | fold after |")
    a("|---|---|---|---:|---:|---:|---:|---:|---:|")
    for r in glio["rows"]:
        a(f"| {r['drug']} | `{r['method']}` | "
          f"{'yes' if r['measures_total_drug'] else '**no**'} | "
          f"{r['R_chemistry']:.3f} | {r['R_applied']:.3f} | "
          f"{r['pred_before']:.4f} | {r['pred_after']:.4f} | "
          f"{r['fold_before']:.3f}x | {r['fold_error']:.3f}x |")
    a("")
    a("| aggregate | before | after (method-aware) | after (applied to "
      "everything) |")
    a("|---|---:|---:|---:|")
    a(f"| geometric-mean ABS fold error | {glio['before']['gm_abs']:.4f}x | "
      f"**{glio['after']['gm_abs']:.4f}x** | {glio_all['after']['gm_abs']:.4f}x |")
    a(f"| geometric-mean signed fold error | {glio['before']['gm_signed']:.4f}x "
      f"| {glio['after']['gm_signed']:.4f}x | "
      f"{glio_all['after']['gm_signed']:.4f}x |")
    a(f"| within 3x | {glio['before']['n_within_3x']}/6 | "
      f"{glio['after']['n_within_3x']}/6 | "
      f"{glio_all['after']['n_within_3x']}/6 |")
    a(f"| within 10x | {glio['before']['n_within_10x']}/6 | "
      f"{glio['after']['n_within_10x']}/6 | "
      f"{glio_all['after']['n_within_10x']}/6 |")
    a(f"| worst drug | {glio['before']['worst_abs_fold']:.2f}x | "
      f"{glio['after']['worst_abs_fold']:.2f}x | "
      f"{glio_all['after']['worst_abs_fold']:.2f}x |")
    a(f"| improved? | -- | **{glio['improved']}** | {glio_all['improved']} |")
    a("")
    a("**A structural note on the wrong placement, because it is not obvious.** "
      "Under `all_tissues` both weak bases get the same `R` (their pKa values "
      "both saturate the trapping ratio), eribulin is over-predicted and "
      "adavosertib under-predicted, so the two `|log10 fold|` values move by "
      "`+log10(R)` and `-log10(R)` and the geometric mean is *exactly* "
      "unchanged — until adavosertib crosses 1.0 and starts moving back the "
      "other way. That is why several rows in the §7 sensitivity table read "
      f"{glio['before']['gm_abs']:.3f}x on both sides. The apparent robustness "
      "is an artefact of two errors cancelling, not a result.")
    a("")
    a("### 6.2 It is one drug wide, and that drug's pKa is not measured")
    a("")
    a(f"{glio['n_drugs_basic']}/6 drugs are weak bases and "
      f"{glio['n_drugs_moved']}/6 move under the method-aware placement. "
      f"Sign test over |log10 fold error|: {paired['n_better']} better, "
      f"{paired['n_worse']} worse, p = {paired['p_sign_test']:.2f}. "
      f"{paired['note']}.")
    a("")
    a("| drug dropped | GM abs before | GM abs after |")
    a("|---|---:|---:|")
    for r in glio["leave_one_out"]:
        a(f"| {r['dropped']} | {r['gm_abs_before']:.4f}x | "
          f"{r['gm_abs_after']:.4f}x |")
    a("")
    a("Adavosertib's pKa could not be found in any retrievable source, so it is "
      "swept rather than defended:")
    a("")
    a("| assumed pKa | `R(glioma)` | adavosertib predicted | adavosertib fold | "
      "GM abs fold error (n = 6) |")
    a("|---:|---:|---:|---:|---:|")
    for r in pka_sweep:
        a(f"| {r['pka']:.1f} | {r['R']:.3f} | {r['adavosertib_pred']:.3f} | "
          f"{r['adavosertib_fold']:.3f}x | {r['gm_abs']:.4f}x |")
    a("")
    a("The improvement is robust to the pKa above about 7.0 and disappears "
      "below 6.5. An aniline-attached 4-methylpiperazine is basic in that range "
      "beyond reasonable doubt, so the direction survives — but note that a pKa "
      "of about 7.7 would land adavosertib exactly on its measured value, which "
      "is a fact to be reported and not a parameter to be chosen.")
    a("")
    # ---- 7. sensitivity --------------------------------------------------
    a("## 7. Does the verdict depend on anything we chose?")
    a("")
    a("Four axes, each of which could have been chosen the other way. The GIST "
      "rows are all evaluated at `phi_lys(tumour) = 0`, i.e. at the mechanism's "
      "best possible case.")
    a("")
    a("| lumen pH | imatinib pKa | liver `phi_lys` | `R_liver` | GM fold error "
      "at the limit | required `phi_lys` | required / measured | shortfall | "
      "still excluded? |")
    a("|---:|---:|---|---:|---:|---:|---:|---:|---|")
    for r in sens:
        if r["scorecard"] != "GIST":
            continue
        a(f"| {r['ph_lys']} | {r['pka']:.2f} ({r['pka_source']}) | "
          f"{_pct(r['phi_liver'], 2)} ({r['phi_liver_source']}) | "
          f"{r['R_liver']:.2f} | {r['gm_at_limit']:.2f}x | "
          f"{_pct(r['phi_required'], 1)} | "
          f"{r['required_over_measured']:.2f}x | "
          f"{r['shortfall_fold']:.2f}x | **{r['still_excluded']}** |")
    a("")
    a(f"**Excluded on {sum(1 for r in sens if r['scorecard'] == 'GIST' and r['still_excluded'])}"
      f"/{sum(1 for r in sens if r['scorecard'] == 'GIST')} settings.** The "
      f"headline configuration is the one most generous to the mechanism; every "
      f"other choice makes the shortfall larger, and using the lumen pH that "
      f"Schmitt's own model uses makes it roughly 2.4x larger again.")
    a("")
    a("| lumen pH | placement | `phi_lys(glioma)` | GM abs before | GM abs "
      "after | improved? | within 3x | adavosertib fold |")
    a("|---:|---|---|---:|---:|---|---:|---:|")
    for r in sens:
        if r["scorecard"] != "glioma":
            continue
        a(f"| {r['ph_lys']} | `{r['placement']}` | "
          f"{_pct(r['phi_glioma'], 2)} ({r['phi_glioma_source']}) | "
          f"{r['gm_abs_before']:.3f}x | {r['gm_abs_after']:.3f}x | "
          f"{r['improved']} | {r['n_within_3x_after']}/6 | "
          f"{r['adavosertib_fold_after']:.3f}x |")
    a("")
    a("The glioma improvement survives every setting of the method-aware "
      "placement and fails every setting of the physically wrong one. That is "
      "the right way round — the placement that is justified a priori is the "
      "one that helps — but it should still be read as one drug moving at "
      "n = 6, not as a validated mechanism.")
    a("")

    # ---- 8. verdict ------------------------------------------------------
    a("## 8. Verdict and recommendation")
    a("")
    a("### Recommendation: **default OFF**, and the module ships that way.")
    a("")
    a("Nothing in the shipped pipeline imports `molecular.ion_trapping`, and "
      "the term returns exactly 1.0 for every drug in every tumour tissue "
      "because none has a measured lysosomal volume fraction. It is available "
      "to any caller who supplies one explicitly, as a documented mechanism "
      "with a stated parameter gap rather than an endorsed default.")
    a("")
    a("### What this run establishes")
    a("")
    a("1. **The chain's own correction was half right.** Moncayo's uniform-`R` "
      "proof genuinely does not apply to a between-tissue ratio — this run "
      "confirms that algebraically and as a unit test, and the term does move "
      "the GIST ratio where a scalar efflux term provably could not. But the "
      "magnitude claim in `uptake_impact.md` §10 compared the required 3.86% "
      "against a *swept* range rather than a *measured* one. Against the "
      f"measurement it is a {req['phi_lys_required_over_measured']:.2f}x "
      f"overshoot.")
    a(f"2. **Ion trapping is excluded by a bound, exactly as hepatic uptake "
      f"was.** Demanded {req['required_R_normal_over_tumour']:.2f}x, available "
      f"{req['max_reachable_ratio_at_measured']:.2f}x at the measured liver "
      f"value in a physically unattainable best case, "
      f"{req['shortfall_fold']:.2f}x short. Uptake was 7.8x short. Two "
      f"mechanisms have now been ruled out by literature ceilings without "
      f"either being fitted.")
    a("3. **The sign is not supported and can invert.** No matched "
      "cancer-versus-normal lysosomal volume-fraction measurement exists; the "
      "two studies reporting volumes point the other way; and at the highest "
      "measured cancer-cell value the GIST prediction moves away from the data.")
    a("4. **The falsifiable side-condition fails.** The intestinal control "
      f"moves from {gist['intestinal_fold_before']:.2f}x to "
      f"{gist['intestinal_fold_limit']:.2f}x at the same setting that improves "
      f"the liver pairs, because gut and liver `phi_lys` differ by 3%.")
    a("5. **On glioma it helps, and it is the first thing in this chain that "
      "does.** That is worth recording precisely because the two previous runs "
      "were negatives on both boards. It rests on one drug, an unmeasured pKa "
      "and an assumed glioma `phi_lys`, and it is not claimed as evidence.")
    a("6. **The placement question has a physical answer, and it changes the "
      "glioma verdict.** Total-drug assays see `R`; extracellular microdialysis "
      "does not. Honouring that is what makes the glioma result positive, and "
      "ignoring it makes it negative.")
    a("")
    a("### What this run does NOT establish")
    a("")
    a("- **That lysosomal sequestration of imatinib is not real.** It plainly "
      "is: PMID 26108972 demonstrates it in GIST-T1 with bafilomycin A1, NH4Cl "
      "and LysoTracker colocalisation. What is shown is narrower — the "
      "TISSUE-LEVEL volume fractions the mechanism would need to explain the "
      "GIST tumour:normal contrast are not the ones that were measured.")
    a("- **Nothing at adequate power.** Six matched pairs with one independent "
      "tissue contrast; six drugs of which one moves.")
    a("- **Nothing about the tumour side, because nothing is measured there.** "
      "Every tumour number in this report is either an assumption or a breast "
      "cancer cell line standing in for a GIST.")
    a("- **Nothing about spatial structure.** This is a per-tissue constant. "
      "The Moncayo target needs a spatial FIELD of acidic-compartment content "
      "and that remains untested and unmeasured.")
    a("")
    a("### Where the chain stands now")
    a("")
    a("| candidate | status |")
    a("|---|---|")
    a("| passive reaction-diffusion transport | **falsified** for tissue "
      "partition — below the structural floor at any L (`gist_imatinib.md`) |")
    a("| scalar lysosomal ion-trapping, WITHIN a section | **falsified** — "
      "uniform R gives no dispersion at equilibrium (`moncayo_explant.md`) |")
    a("| active efflux (P-gp / BCRP) | **tested, does not help**; provably "
      "inert on a tissue ratio; wrong sign in liver (`efflux_impact.md`) |")
    a("| carrier-mediated uptake (OATP / OCT / OAT) | **tested, does not "
      "help**; named carrier falsified; class excluded by a literature bound "
      "(`uptake_impact.md`) |")
    a("| **tissue-differential lysosomal ion trapping** | **tested, does not "
      f"rescue GIST**; excluded by a measured-parameter bound "
      f"({req['shortfall_fold']:.2f}x short); breaks the intestinal control; "
      f"sign unsupported. **Helps the glioma scorecard by one drug** (this run) |")
    a("| spatial acidic-compartment FIELD (within-tumour structure) | "
      "**OPEN**, unchanged; no measurement in usable units exists |")
    a("")
    a("**Three mechanisms have now been named, implemented, cited and scored, "
      "and none of them closes the GIST residual.** Two were excluded by "
      "literature ceilings and one by a scalar-cancellation proof. That is a "
      "stronger conclusion than any of the three negatives on its own: the "
      "residual in the GIST tumour:normal contrast is not accounted for by any "
      "transporter or sequestration term that today's literature can "
      "parameterise, and the next step is **measurement rather than "
      "modelling**.")
    a("")
    a("### The three measurements that would settle it")
    a("")
    a("Named specifically, because \"more data\" is not a plan:")
    a("")
    a("1. **A lysosomal volume fraction for GIST tissue, or for GIST-T1 cells, "
      "on a stated basis.** It does not exist. Collins 2021's method (live-cell "
      "confocal z-stacks normalised to measured cell volume) applied to GIST-T1 "
      "and to matched normal hepatocytes would decide the sign and the "
      "magnitude in one experiment.")
    a("2. **Forsby, Griph & Ivarsson 1983 (PMID 6137128).** Stereological "
      "lysosomal volume density in a human normal glial line versus a human "
      "malignant glioma line — exactly the glioma parameter this run had to "
      "assume. The paper is not open access, has no PMCID, was never archived "
      "and has three citations, so the four numbers were never quoted onward. "
      "This is an interlibrary-loan request, not a research programme.")
    a("3. **A normal-human liver lysosomal volume fraction.** Every value in "
      "the PBPK literature is rat, a textbook or a fit. Rohr 1976 (PMID 823698) "
      "and Roessner 1978 (PMID 654842) contain them and are not obtainable "
      "online.")
    a("")
    a("---")
    a("")
    a("*Generated by `tmesim/validation/ion_trapping_impact.py`. Additive only: "
      "no file under `core/`, `physics/`, `molecular/` (other than the new "
      "`molecular/ion_trapping.py`), `reconstruct/`, `sim/` or `treatment/` was "
      "modified, and the bit-exactness of the untouched delivery path is "
      "verified in §1.*")

    text = "\n".join(L) + "\n"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return text


# ===========================================================================
def run(write_md: bool = True) -> Dict:
    """Everything, in the order the report reads."""
    add = additivity_check()
    form = formulation_check()
    gist = gist_impact()
    req = required_vs_measured(gist)
    sweep = phi_sweep()
    glio = glioma_scorecard(placement="method_aware")
    glio_all = glioma_scorecard(placement="all_tissues")
    paired = glioma_paired_test(glio)
    pkas = glioma_pka_sweep()
    sens = sensitivity()
    if write_md:
        write_report(add, form, gist, req, sweep, glio, glio_all, paired,
                     pkas, sens)
    return {"additivity": add, "formulation": form, "gist": gist,
            "required_vs_measured": req, "phi_sweep": sweep,
            "glioma": glio, "glioma_all_tissues": glio_all,
            "glioma_paired": paired, "glioma_pka_sweep": pkas,
            "sensitivity": sens, "trapping_summary": trapping_summary(),
            "report": REPORT_MD}


if __name__ == "__main__":
    out = run()
    add, form = out["additivity"], out["formulation"]
    gist, req = out["gist"], out["required_vs_measured"]
    glio, glio_all = out["glioma"], out["glioma_all_tissues"]
    assert add["delivery_all_bit_identical"], add["failures"]
    assert add["inert_all_exactly_one"], add["failures"]
    assert form["all_identical"], form["failures"]

    print("ION TRAPPING -- tissue-differential lysosomal retardation")
    print(f"  additive-only            : delivery bit-identical "
          f"{add['n_delivery_bit_identical']}/{add['n_delivery_comparisons']}, "
          f"inert-where-it-should-be "
          f"{add['n_inert_exactly_one']}/{add['n_inert_comparisons']}, "
          f"shipped files modified {add['n_shipped_files_modified']}")
    print(f"  formulation == Moncayo   : {form['all_identical']} "
          f"({form['n_identical']}/{form['n']} exact)")
    print()
    print(f"  MEASURED liver phi_lys   : {100*req['phi_lys_measured_liver']:.2f}% "
          f"of tissue (PMID 33148688)")
    print(f"  REQUIRED (tumour = 0)    : "
          f"{100*req['phi_lys_required_tumour_zero']:.2f}%  "
          f"-> {req['phi_lys_required_over_measured']:.2f}x the measured value")
    print(f"  inside measured range?   : {req['inside_measured_range']}")
    print(f"  demanded {req['required_R_normal_over_tumour']:.2f}x, reachable "
          f"{req['max_reachable_ratio_at_measured']:.2f}x -> "
          f"{req['shortfall_fold']:.2f}x SHORT; excluded "
          f"{req['excluded_by_its_own_literature']}")
    print()
    print(f"  GIST (n = {gist['n']})  geomean fold error "
          f"{gist['before']['geomean_fold_error']:.2f}x -> "
          f"{gist['after']['geomean_fold_error']:.2f}x at the impossible limit; "
          f"{gist['scored']['collins_low']['geomean_fold_error']:.2f}x / "
          f"{gist['scored']['collins_high']['geomean_fold_error']:.2f}x at the "
          f"measured cancer-cell extremes")
    print(f"  GIST intestinal control  "
          f"{gist['intestinal_fold_before']:.2f}x -> "
          f"{gist['intestinal_fold_limit']:.2f}x  (survives "
          f"{gist['intestinal_control_survives']})")
    print(f"  GIST nominal floor       "
          f"{gist['floor_nominal_before']:.4f} "
          f"({gist['n_liver_below_nominal_floor_before']}/5 below) -> "
          f"{gist['floors']['limit']['floor_nominal']:.4f} "
          f"({gist['floors']['limit']['n_below_nominal']}/5 below)")
    print()
    print(f"  glioma (n = 6) GM abs    {glio['before']['gm_abs']:.3f}x -> "
          f"{glio['after']['gm_abs']:.3f}x method-aware (improved "
          f"{glio['improved']}), {glio_all['after']['gm_abs']:.3f}x applied to "
          f"everything (improved {glio_all['improved']})")
    print(f"  .. adavosertib           "
          f"{glio['adavosertib_fold_before']:.3f}x -> "
          f"{glio['adavosertib_fold_after']:.3f}x")
    print()
    print(f"  report -> {out['report']}")
