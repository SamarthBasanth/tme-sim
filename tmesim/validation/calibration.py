"""CALIBRATION: turn the uncalibrated Stage-0 tissue-loss prior into a MECHANISTIC,
one-constant model of perivascular drug consumption -- then test whether that ONE
constant generalizes to a SECOND drug WITHOUT refitting.

Why this module exists
----------------------
The PK reality check (:mod:`tmesim.validation.pk_reality`) exposed a hard failure:
tmesim over-predicts doxorubicin tissue penetration by ~15x. The analytic
reaction-diffusion length ``L = sqrt(D_eff / k)`` comes out ~700 um versus the
published ~40-50 um (Primeau et al. 2005, PMID 16361566). The root cause is the
loss rate ``k = clearance + binding``: the structure-only priors give
``k ~ 3e-4 /s`` (a ~1-hour timescale), whereas reproducing a ~45 um boundary layer
needs ``k = D_eff / L^2 ~ 7e-2 /s`` (a ~15-second timescale) -- roughly **300x
faster**. The missing physics is not slow metabolic clearance; it is **rapid
cellular uptake plus high-capacity intracellular / DNA binding**, which acts as a
near-irreversible perivascular sink and is exactly what carves the steep gradient.

The mechanistic model (ONE free constant)
-----------------------------------------
We keep the frozen engine's first-order sink but split it into a slow passive part
and a fast, molecularly-scaled uptake part::

    k_eff = k_passive + alpha * uptake_proxy(drug)

    k_passive      = drug.clearance      (slow metabolic washout, ~1e-4 /s; the
                     genuinely-slow term, left UN-amplified)
    uptake_proxy   = drug.binding        (the structure-based tissue/intracellular
                     BINDING-CAPACITY rate from binding.tissue_binding_rate: a
                     monotone map of lipophilicity + H-bond-donor "stickiness",
                     i.e. how strongly the drug is sequestered into cell lipids,
                     protein and DNA -- the real driver of perivascular capture)
    alpha          = the SINGLE free scaling constant (dimensionless).

Physical reading of the proxy and alpha
    ``uptake_proxy`` is the module's tissue-binding rate, which we take as a
    structure-defensible *ordering* of high-capacity intracellular sequestration:
    a drug that sticks to tissue more (higher logP / more H-bond donors, e.g. an
    intercalator) is consumed into the perivascular cell layer faster. It is a
    ~30x-wide, well-behaved quantity (dox high, temozolomide low), unlike passive
    membrane permeability which spans ~5 orders of magnitude in the logP model and
    mis-ranks doxorubicin (a low-passive-permeability P-gp substrate whose steep
    gradient is binding-, not permeability-, driven). ``alpha`` is the one factor
    that lifts the deliberately-conservative Stage-0 binding prior up to the true
    perivascular sink strength -- i.e. it encodes the intracellular binding
    CAPACITY (bound:free ratio) the structure-only heuristic cannot know. We fit
    it on doxorubicin ONCE and never touch it again.

Because ``L = sqrt(D_eff / k_eff)`` and ``k_eff`` is *linear* in alpha, the
calibration is a single closed-form solve (no iteration, no optimizer): given the
doxorubicin anchor we invert to ``k_eff* = D_eff / L*^2`` and read off the one
alpha. That is the ONLY fit in this module.

UNIT CONVENTION (corrected): the published anchor is a HALF-distance (40-50 um to
50% of perivascular intensity), while ``L`` is an E-FOLDING length. They differ by
``ln2``: ``L = d_half / ln2``, so ``L* = 45/ln2 ~= 64.9 um``, NOT 45 um. Fitting
to the raw half-distance inflates alpha by ``1/ln^2(2) = 2.081x`` (it gave 309.6
instead of the correct 148.6). ``calibrate_alpha`` now takes an e-folding target
and every anchor band is converted before comparison. See the constant block below.

The honest test (out-of-sample)
-------------------------------
With alpha FROZEN from doxorubicin we PREDICT the penetration of a SECOND drug and
compare to an INDEPENDENT literature measurement we never fit to:

    Mitoxantrone: drug intensity falls to half within ~25-75 um of tumor vessels
    (same immunofluorescence-vs-distance methodology as the dox anchor).
        Rocca, ... Tannock et al. (Tannock lab), "Distribution of the anticancer
        drugs doxorubicin, mitoxantrone and topotecan in tumors and normal
        tissues." Cancer Chemother Pharmacol 2013; 72(3):591-601. PMID 23680920.

Mitoxantrone is a genuinely separate measurement, though (honestly) it is an
anthracenedione DNA-intercalator chemically related to doxorubicin, so it tests
MAGNITUDE TRANSFER within the strongly-bound class. We therefore ALSO run the
opposite-regime qualitative check: temozolomide (small, weakly bound) must be
predicted to penetrate MUCH further than doxorubicin (large L).

Nothing here modifies core, the engine, or the molecular modules; it only drives
them and returns a normal ``core.DrugParams`` the frozen engine consumes.

Units: length um, time s, diffusivity um^2/s, rates 1/s.
"""
from __future__ import annotations

import dataclasses
import math
from typing import Dict, Optional, Union

import numpy as np

from ..core import DrugParams
from ..core.config import TMZ_HYPOXIA_RESISTANCE
from ..molecular import descriptors as _descriptors
from ..molecular import library
from ..molecular.compile_drug import molecule_to_drugparams
from ..physics.diffusion import ALPHA_TORTUOSITY
from .pk_anchors import half_distance_to_efold_length
from .pk_reality import (
    DEFAULT_ECM,
    DOX_CITATION,
    DOX_PENETRATION_MID,
    DOX_PENETRATION_UM,
    ORDER_OF_MAGNITUDE_FACTOR,
    penetration_length_um,
)

# ---------------------------------------------------------------------------
# Anchors.
# ---------------------------------------------------------------------------
# THE fit anchor (doxorubicin, Primeau et al. 2005; re-used from pk_reality).
#
# UNIT CONVENTION (this is the one that used to be got wrong -- read it).
# The published anchor is a HALF-DISTANCE: the distance from a vessel at which
# doxorubicin immunofluorescence falls to 50% of its perivascular value, 40-50 um
# (midpoint 45 um). The simulator's ``penetration_length_um`` is an E-FOLDING
# length L, defined by C(r) ~ exp(-r/L). For an exponential profile the two are
# NOT the same number:
#
#     d_half = L * ln2        =>      L = d_half / ln2 = 1.4427 * d_half
#
# so the 40-50 um half-distance anchor is an L-equivalent of ~58-72 um, midpoint
# ~64.9 um (:func:`pk_anchors.half_distance_to_efold_length`). Calibrating against
# the raw 45 um treats a half-distance as if it were an e-folding length and
# inflates alpha by exactly 1/ln^2(2) = 2.081x (it produced alpha = 309.6 instead
# of 148.6). This module now standardises on the E-FOLD convention throughout:
# alpha is fit to DOX_TARGET_EFOLD_UM and every out-of-sample band is converted
# with the same rule before comparison.
DOX_TARGET_UM = DOX_PENETRATION_MID  # 45 um -- the raw published HALF-distance
#: The fit target actually used: the doxorubicin half-distance midpoint converted
#: to an e-folding length (45 / ln2 ~= 64.9 um), directly comparable to the sim's
#: ``penetration_length_um``.
DOX_TARGET_EFOLD_UM = half_distance_to_efold_length(DOX_PENETRATION_MID)
#: The 40-50 um published half-distance band as an L-equivalent (~58-72 um).
DOX_EFOLD_BAND_UM = (
    half_distance_to_efold_length(DOX_PENETRATION_UM[0]),
    half_distance_to_efold_length(DOX_PENETRATION_UM[1]),
)

# THE out-of-sample anchor (a SECOND drug we never fit to).
# Mitoxantrone is not in the GBM library, so we compile it from a verified SMILES
# (1,4-dihydroxy-5,8-bis{[2-(2-hydroxyethylamino)ethyl]amino}anthracene-9,10-dione,
# MW ~444), the same molecular pipeline the library drugs use.
MITOXANTRONE_NAME = "mitoxantrone"
MITOXANTRONE_SMILES = "OCCNCCNc1ccc(NCCNCCO)c2c1C(=O)c1c(O)ccc(O)c1C2=O"
MITOX_PENETRATION_UM = (25.0, 75.0)        # raw published HALF-distances
MITOX_PENETRATION_MID = 50.0
#: Same half -> e-fold conversion applied to the out-of-sample band (~36-108 um)
#: so the prediction is compared like-for-like with the sim's L.
MITOX_EFOLD_BAND_UM = (
    half_distance_to_efold_length(MITOX_PENETRATION_UM[0]),
    half_distance_to_efold_length(MITOX_PENETRATION_UM[1]),
)
MITOX_EFOLD_MID_UM = half_distance_to_efold_length(MITOX_PENETRATION_MID)
MITOX_CITATION = (
    "Tannock lab, 'Distribution of the anticancer drugs doxorubicin, mitoxantrone "
    "and topotecan in tumors and normal tissues,' Cancer Chemother Pharmacol 2013; "
    "72(3):591-601, PMID 23680920 (mitoxantrone intensity falls to half within "
    "~25-75 um of tumor blood vessels; e-fold equivalent ~36-108 um)"
)


# ---------------------------------------------------------------------------
# Compilation helper (library drug OR mitoxantrone-from-SMILES).
# ---------------------------------------------------------------------------
def _compile(drug: Union[str, DrugParams], ecm_reference: float = DEFAULT_ECM) -> DrugParams:
    """Resolve ``drug`` to a FRESHLY-compiled (un-calibrated) ``DrugParams``.

    A ``DrugParams`` passes through unchanged; a library name is compiled through
    the real molecular pipeline; the special name ``"mitoxantrone"`` is compiled
    from :data:`MITOXANTRONE_SMILES` (it is not in the GBM library).

    Applies the same TMZ-specific ``hypoxia_resistance`` literature override as
    ``core.config.DRUG_LIBRARY["temozolomide-like"]`` (see that constant's and
    ``PDConfig``'s docstrings for the operating-point derivation) -- this is the
    real, GDSC-calibrated molecular-compile pathway every real-data report
    (``kill_calibration.md``, ``real_drug_panel.md``, ``dose_schedule.md``,
    ``cell_state_impact.md``) actually uses for "temozolomide", so the override
    has to be set HERE too, not just on the synthetic library entry, or every
    real-drug-panel run would silently keep using the uncited global default.
    """
    if isinstance(drug, DrugParams):
        return drug
    name = str(drug).strip().lower()
    if name == MITOXANTRONE_NAME:
        mol = _descriptors.from_smiles(MITOXANTRONE_NAME, MITOXANTRONE_SMILES)
    elif name == "palbociclib":
        # Structure curated in molecular.breast_library, not here -- the two
        # drug libraries are kept disjoint by test_disjoint_from_gbm_library
        # (a name must never have two silently-diverging definitions). This is
        # deliberate GBM-context reuse of that SAME real molecule (identical
        # SMILES; see core.config.DRUG_LIBRARY's "palbociclib-like" comment),
        # the same cross-library pattern doxorubicin/cisplatin already use in
        # the opposite direction. potency.resolve_ic50 below is unaffected by
        # which library supplied the structure -- it looks the name up in the
        # glioma GDSC cache directly, so this correctly gets THIS tissue's own
        # measured IC50 (52.19 uM), not the breast one (59.98 uM).
        from ..molecular.breast_library import get_molecule_breast

        mol = get_molecule_breast(name)
    else:
        mol = library.get_molecule(name)
    dp = molecule_to_drugparams(mol, ecm_reference=ecm_reference)
    if name == "temozolomide":
        dp = dataclasses.replace(dp, hypoxia_resistance=TMZ_HYPOXIA_RESISTANCE)
    return dp


# ---------------------------------------------------------------------------
# Mechanistic loss-rate model: k_eff = k_passive + alpha * uptake_proxy.
# ---------------------------------------------------------------------------
def uptake_proxy(drug: DrugParams) -> float:
    """Molecular proxy for RAPID CELLULAR UPTAKE + HIGH-CAPACITY intracellular /
    DNA binding, the driver of steep perivascular gradients (1/s).

    We use the compiled drug's ``binding`` field, i.e. the structure-based tissue
    binding rate from :func:`molecular.binding.tissue_binding_rate` -- a monotone
    map of lipophilicity (logP) and H-bond-donor count representing how strongly a
    drug is sequestered into cell membranes, protein and DNA. It is intentionally
    a *capacity ordering*, not an absolute rate: doxorubicin (intercalator) ranks
    high, temozolomide (small, hydrophilic) ranks low, spanning a well-behaved
    ~30x range. The absolute scale it is missing is supplied by the single fitted
    ``alpha``. Pass a FRESHLY-compiled (un-calibrated) DrugParams.
    """
    return float(drug.binding)


def k_passive(drug: DrugParams) -> float:
    """Slow passive loss (metabolic washout/clearance), 1/s -- the genuinely-slow
    term (~1e-4 /s) that is NOT amplified by alpha. Pass an un-calibrated dp."""
    return float(drug.clearance)


def k_eff(drug: DrugParams, alpha: float) -> float:
    """Effective first-order tissue loss rate ``k_passive + alpha*uptake_proxy``
    (1/s). Pass a FRESHLY-compiled (un-calibrated) DrugParams."""
    return k_passive(drug) + float(alpha) * uptake_proxy(drug)


def _d_eff_tortuous(drug: DrugParams, ecm_density: float) -> float:
    """Effective diffusivity the analytic length uses: ``diffusivity *
    (1 - ALPHA_TORTUOSITY*ecm)`` (um^2/s), matching
    :func:`pk_reality.penetration_length_um` exactly so calibration is consistent
    with the metric it targets."""
    d = float(drug.diffusivity) * (
        1.0 - ALPHA_TORTUOSITY * float(np.clip(ecm_density, 0.0, 1.0))
    )
    return max(d, 1e-6)


# ---------------------------------------------------------------------------
# THE single fit: solve for alpha on doxorubicin.
# ---------------------------------------------------------------------------
def calibrate_alpha(
    target_efold_um: float = DOX_TARGET_EFOLD_UM,
    drug: str = "doxorubicin",
    ecm_density: float = DEFAULT_ECM,
) -> float:
    """Solve for the ONE scaling constant ``alpha`` so the compiled drug's analytic
    penetration length ``L = sqrt(D_eff / k_eff)`` equals ``target_efold_um``.

    ``target_efold_um`` is an **E-FOLDING length**, i.e. the r in ``C(r) ~ exp(-r/L)``
    -- the same quantity :func:`pk_reality.penetration_length_um` returns. It is NOT
    a published half-distance. The default is the doxorubicin 40-50 um half-distance
    anchor (Primeau et al. 2005, PMID 16361566) CONVERTED with
    :func:`pk_anchors.half_distance_to_efold_length`::

        L = d_half / ln2 = 45 / 0.6931 ~= 64.9 um     (band ~58-72 um)

    Passing the raw 45 um half-distance instead treats a half-distance as an
    e-folding length and inflates alpha by exactly ``1/ln^2(2) = 2.081x``; that is
    the historical inconsistency (alpha 309.6 vs the correct 148.6) this signature
    change removes. This is the ONLY fit in the module.

    Closed form (``k_eff`` is linear in alpha, ``L`` monotone in ``k_eff``)::

        k_eff*  = D_eff / target_efold_um^2
        alpha   = (k_eff* - k_passive) / uptake_proxy

    Returns the (dimensionless) alpha; raises if the proxy is non-positive or the
    required k_eff is below the passive floor (target too large to reach).
    """
    dp = _compile(drug, ecm_reference=ecm_density)
    d_eff = _d_eff_tortuous(dp, ecm_density)
    if target_efold_um <= 0:
        raise ValueError(f"target_efold_um must be > 0, got {target_efold_um}")
    k_target = d_eff / (float(target_efold_um) ** 2)
    proxy = uptake_proxy(dp)
    if proxy <= 0.0:
        raise ValueError(f"{drug}: uptake_proxy is non-positive ({proxy}); cannot calibrate")
    alpha = (k_target - k_passive(dp)) / proxy
    if alpha <= 0.0:
        raise ValueError(
            f"{drug}: required k_eff={k_target:.3e} is at/below the passive floor "
            f"k_passive={k_passive(dp):.3e}; target_efold_um={target_efold_um} unreachable"
        )
    return float(alpha)


# ---------------------------------------------------------------------------
# Build a calibrated, engine-ready DrugParams.
# ---------------------------------------------------------------------------
def calibrated_drugparams(drug_name: str, alpha: Optional[float] = None) -> DrugParams:
    """Compile ``drug_name`` and override its loss terms so that
    ``clearance + binding == k_eff = k_passive + alpha*uptake_proxy``.

    The passive metabolic term stays in ``clearance``; the amplified uptake sink
    goes into ``binding`` (``= alpha * uptake_proxy``). The result is a normal,
    frozen-engine ``core.DrugParams`` -- diffusivity, MW, ic50, hill, bbb are
    untouched. If ``alpha`` is None it is calibrated from doxorubicin first.
    """
    if alpha is None:
        alpha = calibrate_alpha()
    dp = _compile(drug_name)
    new_binding = float(alpha) * uptake_proxy(dp)
    return dataclasses.replace(dp, binding=new_binding, clearance=k_passive(dp))


# ---------------------------------------------------------------------------
# Out-of-sample revalidation (alpha FIXED from doxorubicin).
# ---------------------------------------------------------------------------
def revalidate(
    alpha: float,
    second_drug: str = MITOXANTRONE_NAME,
    ecm_density: float = DEFAULT_ECM,
) -> Dict:
    """With ``alpha`` FIXED from doxorubicin, PREDICT a second drug's penetration
    and compare to its independent literature anchor. Nothing is refit here.

    E-FOLD CONVENTION (see :func:`calibrate_alpha`): every published anchor here is
    a HALF-distance and is converted to an e-folding length (``d_half / ln2``)
    before it is compared to the sim's ``L``. The doxorubicin fit target is
    ~64.9 um (not the raw 45 um) and the mitoxantrone band is ~36-108 um (not the
    raw 25-75 um).

    Returns a dict with: the doxorubicin post-calibration check (should hit
    ~64.9 um by construction), the out-of-sample second-drug prediction vs. the
    mitoxantrone e-fold-equivalent band, and the opposite-regime qualitative check
    (temozolomide, small + weakly bound, must penetrate much further than
    doxorubicin). Reports the ratio and an honest pass/fail for each.
    """
    # --- doxorubicin sanity (by construction ~64.9 um e-folding length) ---
    dox_cal = calibrated_drugparams("doxorubicin", alpha)
    dox_L = penetration_length_um(dox_cal, ecm_density=ecm_density)

    # --- OUT-OF-SAMPLE second drug (mitoxantrone), compared on the SAME
    #     (e-fold-converted) footing the fit used. ---
    sec_cal = calibrated_drugparams(second_drug, alpha)
    sec_L = penetration_length_um(sec_cal, ecm_density=ecm_density)
    lo, hi = MITOX_EFOLD_BAND_UM
    ratio_mid = sec_L / MITOX_EFOLD_MID_UM
    in_band = bool(lo <= sec_L <= hi)
    # Two verdict lenses: strict (inside the published range) and the same
    # order-of-magnitude (3x) band pk_reality uses.
    within_3x = bool(
        (1.0 / ORDER_OF_MAGNITUDE_FACTOR) <= ratio_mid <= ORDER_OF_MAGNITUDE_FACTOR
    )
    second = {
        "drug": second_drug,
        "predicted_L_um": sec_L,
        # ``anchor_um`` is the e-fold-CONVERTED band (the like-for-like comparison);
        # the raw published half-distances are kept alongside for traceability.
        "anchor_um": MITOX_EFOLD_BAND_UM,
        "anchor_mid_um": MITOX_EFOLD_MID_UM,
        "half_distance_band_um": MITOX_PENETRATION_UM,
        "citation": MITOX_CITATION,
        "ratio_to_mid": ratio_mid,
        "inside_published_band": in_band,
        "within_3x": within_3x,
        "verdict": "PASS" if in_band else ("PASS(3x)" if within_3x else "FAIL"),
    }

    # --- opposite-regime qualitative check (temozolomide penetrates FURTHER) ---
    tmz_cal = calibrated_drugparams("temozolomide", alpha)
    tmz_L = penetration_length_um(tmz_cal, ecm_density=ecm_density)
    qualitative = {
        "claim": "small, weakly-bound temozolomide penetrates MUCH further than "
                 "doxorubicin (larger L)",
        "temozolomide_L_um": tmz_L,
        "doxorubicin_L_um": dox_L,
        "fold_further": (tmz_L / dox_L) if dox_L > 0 else math.inf,
        "pass": bool(tmz_L > dox_L),
    }

    generalizes = second["verdict"] in ("PASS", "PASS(3x)") and qualitative["pass"]
    return {
        "alpha": float(alpha),
        "ecm_density": ecm_density,
        "doxorubicin_post_cal_L_um": dox_L,
        "second_drug": second,
        "qualitative": qualitative,
        "generalizes": generalizes,
        "verdict": _verdict_text(alpha, dox_L, second, qualitative, generalizes),
    }


def _verdict_text(alpha, dox_L, second, qualitative, generalizes) -> str:
    parts = [
        f"Single fitted constant alpha={alpha:.1f} (the ~150x amplification of the "
        f"structure-only binding prior up to the true perivascular sink), fit on the "
        f"E-FOLD-converted doxorubicin anchor L=d_half/ln2={DOX_TARGET_EFOLD_UM:.1f} um.",
        f"In-sample: doxorubicin post-calibration L={dox_L:.1f} um vs the "
        f"{DOX_EFOLD_BAND_UM[0]:.0f}-{DOX_EFOLD_BAND_UM[1]:.0f} um e-fold-equivalent of "
        f"the published 40-50 um half-distance (hits the target by construction -- this "
        f"is the fit, not a test).",
        f"OUT-OF-SAMPLE: {second['drug']} PREDICTED L={second['predicted_L_um']:.1f} um "
        f"vs the e-fold-equivalent {second['anchor_um'][0]:.0f}-{second['anchor_um'][1]:.0f} um "
        f"(published half-distance {MITOX_PENETRATION_UM[0]:.0f}-{MITOX_PENETRATION_UM[1]:.0f} um; "
        f"{second['ratio_to_mid']:.2f}x the ~{MITOX_EFOLD_MID_UM:.0f} um midpoint) "
        f"=> {second['verdict']}"
        + (" (inside the published band)." if second["inside_published_band"]
           else " (outside band but within 3x)." if second["within_3x"]
           else " (FAILS the 3x band)."),
        f"Opposite regime: temozolomide L={qualitative['temozolomide_L_um']:.0f} um is "
        f"{qualitative['fold_further']:.1f}x further than doxorubicin "
        f"({qualitative['doxorubicin_L_um']:.0f} um) -- "
        + ("correct direction." if qualitative["pass"] else "WRONG direction."),
    ]
    if generalizes:
        parts.append(
            "VERDICT: the one mechanistic constant GENERALIZES -- an independent, "
            "same-method measurement of a second drug falls where the model predicts, "
            "and the weakly-bound regime is captured qualitatively. HONESTY CAVEAT: "
            "mitoxantrone is an anthracenedione intercalator chemically related to "
            "doxorubicin, so this is a strong magnitude-transfer test within the "
            "strongly-bound class rather than a maximally-independent one; the "
            "temozolomide contrast is what probes the opposite regime."
        )
    else:
        parts.append(
            "VERDICT: the constant does NOT cleanly generalize on this test -- see the "
            "out-of-sample ratio above; report it honestly rather than refit."
        )
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Results write-up.
# ---------------------------------------------------------------------------
def write_report(result: Dict, path: str = "results/calibration.md") -> str:
    """Write results/calibration.md from a :func:`revalidate` result dict."""
    import os

    sd = result["second_drug"]
    q = result["qualitative"]
    lo, hi = sd["anchor_um"]
    lines = [
        "# Calibration — one mechanistic constant for perivascular drug consumption",
        "",
        "**Module:** `tmesim/validation/calibration.py` "
        "(run: `python -m tmesim.validation.calibration`)",
        "",
        "**Scope:** fixes the calibration gap the PK reality check exposed "
        "(tmesim over-predicted doxorubicin tissue penetration ~15x: ~700 µm vs the "
        "published ~40–50 µm). We make the tissue loss rate *mechanistic* with exactly "
        "**one** free constant, fit it **once** on doxorubicin, then test whether it "
        "predicts a **second** drug with **no** refitting. Nothing in core, the engine, "
        "or the molecular modules is modified — this module only drives them and returns "
        "a normal `core.DrugParams`.",
        "",
        "---",
        "",
        "## The mechanism and the proxy",
        "",
        "The frozen transport PDE loses drug at a first-order rate "
        "`k = clearance + binding`, giving penetration length `L = sqrt(D_eff / k)`. "
        "The structure-only priors give `k ≈ 3e-4 /s` (a ~1-hour timescale) → `L ≈ 700 µm`. "
        f"Reproducing the {DOX_TARGET_EFOLD_UM:.0f} µm boundary layer (the e-fold "
        "equivalent of the published 40–50 µm half-distance — see the convention "
        "section below) needs `k = D_eff / L² ≈ 3e-2 /s` (a ~30-second timescale) — "
        "about **150× faster**. The missing physics is not "
        "slow metabolic clearance; it is **rapid cellular uptake + high-capacity "
        "intracellular / DNA binding**, a near-irreversible perivascular sink that is "
        "what actually carves the steep gradient.",
        "",
        "We keep the engine's first-order sink but split it:",
        "",
        "```",
        "k_eff = k_passive + alpha * uptake_proxy(drug)",
        "",
        "k_passive    = drug.clearance   (slow metabolic washout, ~1e-4 /s; NOT amplified)",
        "uptake_proxy = drug.binding     (structure-based tissue/DNA BINDING-CAPACITY rate,",
        "                                 binding.tissue_binding_rate: monotone in logP + HBD)",
        "alpha        = the ONE free scaling constant (dimensionless)",
        "```",
        "",
        "**Why this proxy.** `uptake_proxy` is the molecular tissue-binding rate — a "
        "structure-defensible *ordering* of how strongly a drug is sequestered into cell "
        "lipids, protein and DNA (intercalators rank high, small hydrophilic drugs low). "
        "It is a well-behaved ~30× range (doxorubicin high, temozolomide low), unlike "
        "passive membrane permeability, which spans ~5 orders of magnitude in the logP "
        "model and mis-ranks doxorubicin (a low-passive-permeability P-gp substrate whose "
        "gradient is binding-, not permeability-, driven). `alpha` is the single factor "
        "that lifts the deliberately-conservative Stage-0 binding prior up to the true "
        "perivascular sink strength — it encodes the intracellular binding **capacity** "
        "(bound:free ratio) the structure-only heuristic cannot know.",
        "",
        "Because `k_eff` is linear in `alpha`, calibration is one closed-form solve: "
        "`k_eff* = D_eff / L*²`, `alpha = (k_eff* − k_passive) / uptake_proxy`. No "
        "optimizer, no iteration.",
        "",
        "---",
        "",
        "## Unit convention: half-distance vs e-folding length (CORRECTED)",
        "",
        "The published anchors are **half-distances** `d_half` — the distance from a "
        "vessel at which drug intensity falls to 50 % of its perivascular value. The "
        "simulator's `penetration_length_um` is an **e-folding length** `L`, defined by "
        "`C(r) ~ exp(−r/L)`. For an exponential profile these are different numbers:",
        "",
        "```",
        "d_half = L · ln2        =>        L = d_half / ln2 = 1.4427 · d_half",
        "```",
        "",
        f"So doxorubicin's published **40–50 µm half-distance** is an L-equivalent of "
        f"**{DOX_EFOLD_BAND_UM[0]:.0f}–{DOX_EFOLD_BAND_UM[1]:.0f} µm**, midpoint "
        f"**{DOX_TARGET_EFOLD_UM:.1f} µm** "
        f"(`pk_anchors.half_distance_to_efold_length`). `calibrate_alpha` now takes its "
        "target as an **e-folding length** and defaults to that converted midpoint.",
        "",
        "> **Correction notice.** An earlier revision of this file reported "
        "**alpha = 309.6**, obtained by fitting to the raw 45 µm half-distance *as if it "
        "were* an e-folding length. That is the wrong convention, and it inflates alpha "
        "by exactly `1/ln²(2) = 2.081×`. The corrected, e-fold-converted fit gives "
        f"**alpha = {result['alpha']:.1f}** — the value `results/pk_multidrug.md` already "
        "used, so the two files are now consistent. Any figure or table still quoting "
        "α = 309.6 predates this correction.",
        "",
        "---",
        "",
        "## The single calibrated constant",
        "",
        f"Fit on doxorubicin against the **{DOX_TARGET_EFOLD_UM:.1f} µm e-fold-equivalent** "
        f"of the 40–50 µm half-distance anchor "
        f"(Primeau et al., Clin Cancer Res 2005; 11(24):8782–8788, PMID 16361566):",
        "",
        f"> **alpha = {result['alpha']:.1f}**   (dimensionless; a ~150× amplification of "
        f"the structure-only binding prior)",
        "",
        "This is the **only** fit. Everything below uses this same alpha.",
        "",
        "### In-sample check (doxorubicin) — this is the fit, not a test",
        "",
        f"Doxorubicin post-calibration `L = {result['doxorubicin_post_cal_L_um']:.1f} µm` "
        f"vs the {DOX_EFOLD_BAND_UM[0]:.0f}–{DOX_EFOLD_BAND_UM[1]:.0f} µm e-fold-equivalent "
        f"band → hits ~{DOX_TARGET_EFOLD_UM:.0f} µm **by construction**.",
        "",
        "---",
        "",
        "## The real test: out-of-sample second drug (NOT fit)",
        "",
        f"With alpha **fixed at {result['alpha']:.1f}**, predict **{sd['drug']}** and "
        f"compare to an independent, same-method literature measurement — converted to "
        f"an e-folding length with the same `d_half/ln2` rule used for the fit:",
        "",
        f"> {sd['citation']}",
        "",
        "| drug (out-of-sample) | predicted L | published half-distance | e-fold-equivalent band | ratio to mid | verdict |",
        "|---|---:|---:|---:|---:|:--:|",
        f"| {sd['drug']} | **{sd['predicted_L_um']:.1f} µm** | "
        f"{sd['half_distance_band_um'][0]:.0f}–{sd['half_distance_band_um'][1]:.0f} µm | "
        f"{lo:.0f}–{hi:.0f} µm | "
        f"{sd['ratio_to_mid']:.2f}× | **{sd['verdict']}** |",
        "",
        (f"The predicted penetration falls **inside** the {lo:.0f}–{hi:.0f} µm "
         "e-fold-equivalent band — the one constant transfers to a second drug it was "
         "never fit to."
         if sd["inside_published_band"] else
         "The prediction is **outside** the converted band; reported honestly, not refit."),
        "",
        "### Opposite-regime qualitative check (temozolomide)",
        "",
        f"Temozolomide (small, weakly bound) is predicted at `L = "
        f"{q['temozolomide_L_um']:.0f} µm`, **{q['fold_further']:.1f}× further** than "
        f"doxorubicin (`{q['doxorubicin_L_um']:.0f} µm`) — "
        + ("correct direction (the model captures the weakly-bound regime it was not fit on)."
           if q["pass"] else "WRONG direction."),
        "",
        "---",
        "",
        "## Honest verdict",
        "",
        result["verdict"],
        "",
        "### Caveats (the honesty line)",
        "",
        "- **One anchor sets the scale.** alpha is fit to a single doxorubicin number; "
        "the absolute penetration of *every* drug rides on that one measurement.",
        "- **The anchor is a half-distance, converted.** `L = d_half/ln2` assumes the "
        "perivascular profile is a single exponential. Real immunofluorescence profiles "
        "are only approximately exponential, so the conversion carries its own (modest) "
        "error — but applying it is strictly better than the previous convention, which "
        "silently equated the two lengths and mis-scaled alpha by 2.08×.",
        "- **Mitoxantrone is chemically related to doxorubicin** (both anthracenedione "
        "DNA-intercalators). The measurement is genuinely independent, but this is a "
        "strong *magnitude-transfer* test within the strongly-bound class, not a "
        "maximally-independent one. The temozolomide contrast is what probes the "
        "opposite (weakly-bound, deep-penetration) regime — and only qualitatively, "
        "since temozolomide has no clean perivascular-distance dataset.",
        "- **The proxy is an ordering, not a first-principles rate.** "
        "`tissue_binding_rate` is a low-confidence logP/HBD heuristic; alpha corrects "
        "its scale but inherits its ranking errors. A drug the heuristic mis-ranks "
        "(e.g. active transport, efflux, ion-trapping) would be mis-predicted.",
        "- **Not fit to the second drug.** No parameter here was tuned to mitoxantrone "
        "or temozolomide; both are pure predictions.",
        "",
    ]
    text = "\n".join(lines)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


# ---------------------------------------------------------------------------
def _print(result: Dict) -> None:
    sd = result["second_drug"]
    q = result["qualitative"]
    print("=" * 88)
    print("tmesim CALIBRATION -- one mechanistic constant for perivascular drug consumption")
    print("=" * 88)
    print(f"Mechanism: k_eff = k_passive(clearance) + alpha * uptake_proxy(tissue binding)")
    print(f"Convention: anchors are HALF-distances; fit target is the e-folding length "
          f"L = d_half/ln2")
    print(f"THE single fit (doxorubicin 40-50 um half-dist -> L = "
          f"{DOX_TARGET_EFOLD_UM:.1f} um):  alpha = {result['alpha']:.2f}")
    print("-" * 88)
    print(f"In-sample  doxorubicin  L = {result['doxorubicin_post_cal_L_um']:6.1f} um "
          f"(e-fold band {DOX_EFOLD_BAND_UM[0]:.0f}-{DOX_EFOLD_BAND_UM[1]:.0f} um; "
          f"by construction)")
    print("-" * 88)
    print("OUT-OF-SAMPLE (alpha fixed, NOT refit):")
    lo, hi = sd["anchor_um"]
    print(f"  {sd['drug']:13s} predicted L = {sd['predicted_L_um']:6.1f} um   "
          f"e-fold band {lo:.0f}-{hi:.0f} um   ratio={sd['ratio_to_mid']:.2f}x   "
          f"[{sd['verdict']}]")
    print(f"    anchor: {sd['citation']}")
    print(f"  temozolomide  L = {q['temozolomide_L_um']:.0f} um is "
          f"{q['fold_further']:.1f}x further than doxorubicin "
          f"({q['doxorubicin_L_um']:.0f} um)  [{'PASS' if q['pass'] else 'FAIL'}]")
    print("=" * 88)
    print("VERDICT")
    print("=" * 88)
    print(result["verdict"])
    print("=" * 88)


if __name__ == "__main__":
    import time as _time

    t0 = _time.time()
    alpha = calibrate_alpha(target_efold_um=DOX_TARGET_EFOLD_UM, drug="doxorubicin")
    result = revalidate(alpha)
    _print(result)
    out = write_report(result)
    print(f"\n[wrote {out}]  [runtime {_time.time() - t0:.1f} s]")
