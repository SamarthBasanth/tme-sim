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
doxorubicin anchor ``L* = 45 um`` we invert to ``k_eff* = D_eff / L*^2`` and read
off the one alpha. That is the ONLY fit in this module.

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
from ..molecular import descriptors as _descriptors
from ..molecular import library
from ..molecular.compile_drug import molecule_to_drugparams
from ..physics.diffusion import ALPHA_TORTUOSITY
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
DOX_TARGET_UM = DOX_PENETRATION_MID  # 45 um

# THE out-of-sample anchor (a SECOND drug we never fit to).
# Mitoxantrone is not in the GBM library, so we compile it from a verified SMILES
# (1,4-dihydroxy-5,8-bis{[2-(2-hydroxyethylamino)ethyl]amino}anthracene-9,10-dione,
# MW ~444), the same molecular pipeline the library drugs use.
MITOXANTRONE_NAME = "mitoxantrone"
MITOXANTRONE_SMILES = "OCCNCCNc1ccc(NCCNCCO)c2c1C(=O)c1c(O)ccc(O)c1C2=O"
MITOX_PENETRATION_UM = (25.0, 75.0)
MITOX_PENETRATION_MID = 50.0
MITOX_CITATION = (
    "Tannock lab, 'Distribution of the anticancer drugs doxorubicin, mitoxantrone "
    "and topotecan in tumors and normal tissues,' Cancer Chemother Pharmacol 2013; "
    "72(3):591-601, PMID 23680920 (mitoxantrone intensity falls to half within "
    "~25-75 um of tumor blood vessels)"
)


# ---------------------------------------------------------------------------
# Compilation helper (library drug OR mitoxantrone-from-SMILES).
# ---------------------------------------------------------------------------
def _compile(drug: Union[str, DrugParams], ecm_reference: float = DEFAULT_ECM) -> DrugParams:
    """Resolve ``drug`` to a FRESHLY-compiled (un-calibrated) ``DrugParams``.

    A ``DrugParams`` passes through unchanged; a library name is compiled through
    the real molecular pipeline; the special name ``"mitoxantrone"`` is compiled
    from :data:`MITOXANTRONE_SMILES` (it is not in the GBM library).
    """
    if isinstance(drug, DrugParams):
        return drug
    name = str(drug).strip().lower()
    if name == MITOXANTRONE_NAME:
        mol = _descriptors.from_smiles(MITOXANTRONE_NAME, MITOXANTRONE_SMILES)
    else:
        mol = library.get_molecule(name)
    return molecule_to_drugparams(mol, ecm_reference=ecm_reference)


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
    target_um: float = DOX_TARGET_UM,
    drug: str = "doxorubicin",
    ecm_density: float = DEFAULT_ECM,
) -> float:
    """Solve for the ONE scaling constant ``alpha`` so the compiled drug's analytic
    penetration length ``L = sqrt(D_eff / k_eff)`` equals ``target_um`` (default
    45 um, the doxorubicin anchor). This is the ONLY fit in the module.

    Closed form (``k_eff`` is linear in alpha, ``L`` monotone in ``k_eff``)::

        k_eff*  = D_eff / target_um^2
        alpha   = (k_eff* - k_passive) / uptake_proxy

    Returns the (dimensionless) alpha; raises if the proxy is non-positive or the
    required k_eff is below the passive floor (target too large to reach).
    """
    dp = _compile(drug, ecm_reference=ecm_density)
    d_eff = _d_eff_tortuous(dp, ecm_density)
    if target_um <= 0:
        raise ValueError(f"target_um must be > 0, got {target_um}")
    k_target = d_eff / (float(target_um) ** 2)
    proxy = uptake_proxy(dp)
    if proxy <= 0.0:
        raise ValueError(f"{drug}: uptake_proxy is non-positive ({proxy}); cannot calibrate")
    alpha = (k_target - k_passive(dp)) / proxy
    if alpha <= 0.0:
        raise ValueError(
            f"{drug}: required k_eff={k_target:.3e} is at/below the passive floor "
            f"k_passive={k_passive(dp):.3e}; target_um={target_um} unreachable"
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

    Returns a dict with: the doxorubicin post-calibration check (should hit ~45 um
    by construction), the out-of-sample second-drug prediction vs. the mitoxantrone
    ~25-75 um band, and the opposite-regime qualitative check (temozolomide, small
    + weakly bound, must penetrate much further than doxorubicin). Reports the
    ratio and an honest pass/fail for each.
    """
    # --- doxorubicin sanity (by construction ~45 um) ---
    dox_cal = calibrated_drugparams("doxorubicin", alpha)
    dox_L = penetration_length_um(dox_cal, ecm_density=ecm_density)

    # --- OUT-OF-SAMPLE second drug (mitoxantrone) ---
    sec_cal = calibrated_drugparams(second_drug, alpha)
    sec_L = penetration_length_um(sec_cal, ecm_density=ecm_density)
    lo, hi = MITOX_PENETRATION_UM
    ratio_mid = sec_L / MITOX_PENETRATION_MID
    in_band = bool(lo <= sec_L <= hi)
    # Two verdict lenses: strict (inside the published range) and the same
    # order-of-magnitude (3x) band pk_reality uses.
    within_3x = bool(
        (1.0 / ORDER_OF_MAGNITUDE_FACTOR) <= ratio_mid <= ORDER_OF_MAGNITUDE_FACTOR
    )
    second = {
        "drug": second_drug,
        "predicted_L_um": sec_L,
        "anchor_um": MITOX_PENETRATION_UM,
        "anchor_mid_um": MITOX_PENETRATION_MID,
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
        f"Single fitted constant alpha={alpha:.1f} (the ~300x amplification of the "
        f"structure-only binding prior up to the true perivascular sink).",
        f"In-sample: doxorubicin post-calibration L={dox_L:.1f} um vs 40-50 um anchor "
        f"(hits ~45 um by construction -- this is the fit, not a test).",
        f"OUT-OF-SAMPLE: {second['drug']} PREDICTED L={second['predicted_L_um']:.1f} um "
        f"vs published {second['anchor_um'][0]:.0f}-{second['anchor_um'][1]:.0f} um "
        f"({second['ratio_to_mid']:.2f}x the ~{MITOX_PENETRATION_MID:.0f} um midpoint) "
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
        "Reproducing a ~45 µm boundary layer needs `k = D_eff / L² ≈ 7e-2 /s` "
        "(a ~15-second timescale) — about **300× faster**. The missing physics is not "
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
        "## The single calibrated constant",
        "",
        f"Fit on doxorubicin against the **45 µm** anchor "
        f"(Primeau et al., Clin Cancer Res 2005; 11(24):8782–8788, PMID 16361566):",
        "",
        f"> **alpha = {result['alpha']:.1f}**   (dimensionless; the ~300× amplification "
        f"the diagnosis predicted)",
        "",
        "This is the **only** fit. Everything below uses this same alpha.",
        "",
        "### In-sample check (doxorubicin) — this is the fit, not a test",
        "",
        f"Doxorubicin post-calibration `L = {result['doxorubicin_post_cal_L_um']:.1f} µm` "
        f"vs the 40–50 µm anchor → hits ~45 µm **by construction**.",
        "",
        "---",
        "",
        "## The real test: out-of-sample second drug (NOT fit)",
        "",
        f"With alpha **fixed at {result['alpha']:.1f}**, predict **{sd['drug']}** and "
        f"compare to an independent, same-method literature measurement:",
        "",
        f"> {sd['citation']}",
        "",
        "| drug (out-of-sample) | predicted L | published band | ratio to mid | verdict |",
        "|---|---:|---:|---:|:--:|",
        f"| {sd['drug']} | **{sd['predicted_L_um']:.1f} µm** | {lo:.0f}–{hi:.0f} µm | "
        f"{sd['ratio_to_mid']:.2f}× | **{sd['verdict']}** |",
        "",
        ("The predicted penetration falls **inside** the published 25–75 µm band — the "
         "one constant transfers to a second drug it was never fit to."
         if sd["inside_published_band"] else
         "The prediction is **outside** the published band; reported honestly, not refit."),
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
    print(f"THE single fit (doxorubicin -> 45 um anchor):  alpha = {result['alpha']:.2f}")
    print("-" * 88)
    print(f"In-sample  doxorubicin  L = {result['doxorubicin_post_cal_L_um']:6.1f} um "
          f"(anchor 40-50 um; by construction)")
    print("-" * 88)
    print("OUT-OF-SAMPLE (alpha fixed, NOT refit):")
    lo, hi = sd["anchor_um"]
    print(f"  {sd['drug']:13s} predicted L = {sd['predicted_L_um']:6.1f} um   "
          f"published {lo:.0f}-{hi:.0f} um   ratio={sd['ratio_to_mid']:.2f}x   "
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
    alpha = calibrate_alpha(target_um=DOX_TARGET_UM, drug="doxorubicin")
    result = revalidate(alpha)
    _print(result)
    out = write_report(result)
    print(f"\n[wrote {out}]  [runtime {_time.time() - t0:.1f} s]")
