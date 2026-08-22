"""The transvascular delivery law, and the disclosure that has to travel with it.

Every delivery number this app renders — mean tumour AUC, therapeutic coverage,
penetration length, kill fraction — and therefore the SHAPE of every exposure
map built out of them, is produced under one of two source terms in the
transport solve. (This sentence used to say "every cross-drug RANKING built out
of them". Compound ordering was withdrawn on 2026-08-21 after it was measured —
results/pk_ranking_robustness.md, results/capacity_ranking_robustness.md — and
the law matters MORE under the replacement claim, not less: a plasma or a
capacity factor cancels out of a within-map ratio, and the source term does
not.) The two terms are:

* ``S = PS*Cp``          the historical ONE-WAY leak. No path from interstitium
                         back to blood, so the steady state is unbounded above.
* ``S = PS*(Cp - C)``    gradient-driven EXCHANGE. Same ``PS``, applied to both
                         halves of the gradient, so ``C_ss <= Cp`` by
                         construction.

Which one is in force is :attr:`tmesim.core.PDEConfig.transvascular_exchange`,
and the choice is NOT cosmetic: it re-weights delivery by permeability, which
spans 47x across the GBM drug list, and it changes the top-1 GBM pick
(cisplatin -> vincristine, ``results/transvascular_exchange_revalidation.md``
§4.1). That top-1 flip is now reported as EVIDENCE OF MAGNITUDE rather than as a
consequence for a product output, because compound ordering was withdrawn on
2026-08-21 (``results/pk_ranking_robustness.md``,
``results/capacity_ranking_robustness.md``). The consequence that remains is
sharper, not softer: the source term sets the SHAPE of the exposure field, and a
shape is precisely what a within-map ratio does NOT cancel. A user who reads a
map off this Studio without being told which law produced it has been handed a
gradient with a hidden branch in it.

**Everything here is READ FROM ``PDEConfig`` AT CALL TIME.** Nothing about the
law — not the formula, not the word "default", not the direction of the caveat —
is hard-coded to today's shipped value. If the default flips, the disclosure
rewrites itself rather than going quietly stale, which is the only version of
this that stays true without someone remembering to come back.

The external anchor and the magnitude of the law's effect are quoted verbatim
from the repo's own measurements:

* ``results/breast_kp_validation.md`` — against human breast microdialysis the
  one-way leak puts 5-fluorouracil at 7.80x plasma where 0.61x was measured, and
  methotrexate at 0.057x where 0.60x was measured. OPPOSITE directions, so no
  single constant absorbs both.
* ``results/transvascular_exchange_revalidation.md`` §4.1 — correcting the law
  changes top-1 on GBM (cisplatin -> vincristine) and swaps 4 of 45 pairs.
"""
from __future__ import annotations

from typing import Dict, Optional

__all__ = ["LEAK_FORM", "EXCHANGE_FORM", "REFERENCES", "delivery_law"]

LEAK_FORM = "S = PS*Cp"
EXCHANGE_FORM = "S = PS*(Cp - C)"

REFERENCES = (
    "results/breast_kp_validation.md, "
    "results/transvascular_exchange_revalidation.md"
)

# The measured miss, in both directions. Kept as one string because the two
# halves are only damning together: a single-direction error could be absorbed
# by rescaling PS, and this one cannot.
_ANCHOR_MISS = (
    "against human breast microdialysis the one-way leak puts 5-fluorouracil at "
    "7.80x plasma where 0.61x was measured, and methotrexate at 0.057x where "
    "0.60x was measured — OPPOSITE directions, so no single constant absorbs "
    "both"
)

# How much the law moves the answer, measured on GBM. This used to be presented
# as "the ranking consequence"; the ranking is withdrawn, and the measurement is
# kept because it is the only external calibration of the law's MAGNITUDE that
# this repo has. Read it as: the source term is large enough to reorder five
# drugs, therefore it is large enough to retilt one map.
_LAW_MAGNITUDE = (
    "the corrected law re-weights delivery by permeability (which spans 47x "
    "across this drug list) by enough to have CHANGED THE TOP-1 GBM DRUG, "
    "cisplatin -> vincristine, swapping 4 of 45 pairs, back when this repo "
    "still published an ordering"
)


def _shipped_default() -> bool:
    """The value ``PDEConfig`` ships with, read now, not remembered.

    Imported inside the function so that importing this module never pulls in
    numpy/scipy through ``tmesim.core`` — the API imports it at module scope.
    """
    from tmesim.core import PDEConfig

    return bool(getattr(PDEConfig(), "transvascular_exchange", False))


def delivery_law(in_effect: Optional[bool] = None) -> Dict:
    """Describe the delivery law that produced a set of numbers.

    :param in_effect: the value actually used by the run being described. When
        ``None`` the shipped default is reported, which is what a surface that
        is not scoped to one run (the report limits block, the cross-specimen
        comparison) can honestly say.
    :returns: ``{transvascular_exchange, shipped_default, is_shipped_default,
        law, alternative_law, title, headline, caveat, reference}``.

    ``caveat`` is written to match the shape ``app/web/js/ui/inspect.js``
    ``splitCaveat()`` parses: a ``<subject>: <body>`` prefix under 88 characters
    and a single trailing parenthesis carrying the repo citations. It contains
    "MUST NOT", so it is counted and styled as a HARD STOP alongside the other
    hard stops in the caveats panel — which is the correct severity. A map read
    as law-independent is a wrong map, not a slightly qualified one: the two
    laws put the drug in different places, not merely at different levels.
    """
    shipped = _shipped_default()
    on = shipped if in_effect is None else bool(in_effect)
    law = EXCHANGE_FORM if on else LEAK_FORM
    other = LEAK_FORM if on else EXCHANGE_FORM
    is_default = on == shipped

    # How this run's value relates to what the package ships. Stated explicitly
    # in every branch, so "the shipped default" is never an assumption the
    # reader has to make.
    if is_default:
        standing = (
            f"PDEConfig.transvascular_exchange={on}, which is the SHIPPED "
            "DEFAULT — nothing on this screen opted into it"
        )
    else:
        standing = (
            f"PDEConfig.transvascular_exchange={on}, which OVERRIDES the "
            f"shipped default of {shipped}"
        )

    if on:
        title = f"DELIVERY LAW {law} (exchange)"
        body = (
            f"every delivery number and the SHAPE of every exposure map here "
            f"was produced by the gradient-driven EXCHANGE source {law} "
            f"({standing}), which is bounded by the plasma level by "
            f"construction. The alternative one-way leak {other} is the form "
            f"every stored reference result, calibration and regression in this "
            f"repo was produced under, and the two laws DISAGREE about where "
            f"the drug goes, not merely how much of it there is: "
            f"{_LAW_MAGNITUDE}. A plasma or a tissue-capacity factor cancels "
            f"out of a within-map ratio; this does not. Numbers and maps here "
            f"MUST NOT be compared against, or read as continuous with, any "
            f"result produced under the leak"
        )
    else:
        title = f"DELIVERY LAW {law} (one-way leak)"
        body = (
            f"every delivery number and the SHAPE of every exposure map here "
            f"was produced by the ONE-WAY vascular leak {law} ({standing}). It "
            f"has no path from interstitium back to blood, so its steady state "
            f"is unbounded above, and it MISSES ITS EXTERNAL ANCHOR IN BOTH "
            f"DIRECTIONS: {_ANCHOR_MISS}. Correcting the form to {other} is not "
            f"a rescale — {_LAW_MAGNITUDE} — so it moves the gradient, which is "
            f"the one thing a within-map ratio does not cancel. The map shown "
            f"here MUST NOT be read as law-independent"
        )

    # Short enough for a table footnote or a one-line panel note, and it still
    # names the law, its status, and the consequence.
    if on:
        headline = (
            f"Delivery law: {law} (exchange, "
            + ("shipped default" if is_default else f"overriding the shipped {shipped}")
            + f"). The one-way leak {other} that every stored result in this "
            f"repo used gives a DIFFERENT top-1 GBM drug."
        )
    else:
        headline = (
            f"Delivery law: {law} (one-way leak, "
            + ("shipped default" if is_default else f"overriding the shipped {shipped}")
            + f"), which misses human microdialysis in both directions. "
            f"Correcting it to {other} changes the top-1 GBM drug "
            f"(cisplatin -> vincristine)."
        )

    return {
        "transvascular_exchange": on,
        "shipped_default": shipped,
        "is_shipped_default": is_default,
        "law": law,
        "alternative_law": other,
        "title": title,
        "headline": headline,
        "caveat": f"{title}: {body} ({REFERENCES}).",
        "reference": REFERENCES,
    }
