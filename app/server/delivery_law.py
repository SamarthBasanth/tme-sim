"""The transvascular delivery law, and the disclosure that has to travel with it.

Every delivery number this app renders — mean tumour AUC, therapeutic coverage,
penetration length, kill fraction — and therefore every cross-drug RANKING built
out of them, is produced under one of two source terms in the transport solve:

* ``S = PS*Cp``          the historical ONE-WAY leak. No path from interstitium
                         back to blood, so the steady state is unbounded above.
* ``S = PS*(Cp - C)``    gradient-driven EXCHANGE. Same ``PS``, applied to both
                         halves of the gradient, so ``C_ss <= Cp`` by
                         construction.

Which one is in force is :attr:`tmesim.core.PDEConfig.transvascular_exchange`,
and the choice is NOT cosmetic: it re-weights the panel by permeability, which
spans 47x across the GBM drug list, and it changes the top-1 GBM pick
(cisplatin -> vincristine, ``results/transvascular_exchange_revalidation.md``
§4.1). A user who reads a ranking off this Studio without being told which law
produced it has been handed a decision-shaped number with a hidden branch in it.

**Everything here is READ FROM ``PDEConfig`` AT CALL TIME.** Nothing about the
law — not the formula, not the word "default", not the direction of the caveat —
is hard-coded to today's shipped value. If the default flips, the disclosure
rewrites itself rather than going quietly stale, which is the only version of
this that stays true without someone remembering to come back.

The external anchor and the ranking consequence are quoted verbatim from the
repo's own measurements:

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

# The ranking consequence, measured on GBM.
_RANK_FLIP = (
    "the corrected law re-weights the panel by permeability (which spans 47x "
    "across this drug list) and CHANGES THE TOP-1 GBM DRUG, cisplatin -> "
    "vincristine, swapping 4 of 45 pairs"
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
        ranking) can honestly say.
    :returns: ``{transvascular_exchange, shipped_default, is_shipped_default,
        law, alternative_law, title, headline, caveat, reference}``.

    ``caveat`` is written to match the shape ``app/web/js/ui/inspect.js``
    ``splitCaveat()`` parses: a ``<subject>: <body>`` prefix under 88 characters
    and a single trailing parenthesis carrying the repo citations. It contains
    "MUST NOT", so it is counted and styled as a HARD STOP alongside the other
    hard stops in the caveats panel — which is the correct severity. A ranking
    read as law-independent is a wrong ranking, not a slightly qualified one.
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
            f"every delivery number and every ranking here was produced by the "
            f"gradient-driven EXCHANGE source {law} ({standing}), which is "
            f"bounded by the plasma level by construction. The alternative "
            f"one-way leak {other} is the form every stored reference result, "
            f"calibration and regression in this repo was produced under, and "
            f"the two laws DISAGREE on the ordering: {_RANK_FLIP}. Numbers here "
            f"MUST NOT be compared against, or read as continuous with, any "
            f"result produced under the leak"
        )
    else:
        title = f"DELIVERY LAW {law} (one-way leak)"
        body = (
            f"every delivery number and every ranking here was produced by the "
            f"ONE-WAY vascular leak {law} ({standing}). It has no path from "
            f"interstitium back to blood, so its steady state is unbounded "
            f"above, and it MISSES ITS EXTERNAL ANCHOR IN BOTH DIRECTIONS: "
            f"{_ANCHOR_MISS}. Correcting the form to {other} is not a rescale — "
            f"{_RANK_FLIP}. The ordering shown here MUST NOT be read as "
            f"law-independent"
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
