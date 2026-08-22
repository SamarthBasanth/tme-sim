"""The SINGLE-COMPOUND MAP DOCUMENT — the deliverable the pipeline produces.

A **report** is ONE compound simulated on ONE real tissue: where the drug
reaches, how much of the tumour it never reaches, and the provenance and limits
of every number that says so. It is the artefact a customer keeps and shows a
colleague.

THE RANKED PANEL WAS WITHDRAWN ON 2026-08-21
--------------------------------------------
Until that date this module produced "one ranked, mechanistically attributed
drug panel", and the ORDERING was named as its defensible output. That ordering
was measured this week and it does not hold:

* ``results/pk_ranking_robustness.md`` sec. 0 — per-compound plasma exposure
  across the compounds the GBM panel ranks spans **113.1x**, while the published
  GBM order first breaks at a per-compound spread of **1.53x**. A single
  compound's own licensed dose range (irinotecan, 1,660 vs 3,392 ng/mL Cmax on
  one FDA label table, a 2.04x span) contains a swap that needs only 1.961x. The
  ordering was conditional on a quantity that varies by two orders of magnitude
  between the very compounds being ordered.
* ``results/capacity_ranking_robustness.md`` sec. 0 — under this repo's OWN
  structural capacity rule the breast top-1 CHANGES (``vinorelbine`` ->
  ``doxorubicin``), and no sampled per-compound capacity arm reaches
  decision-grade on any panel. The measured spread of the capacity term the
  engine assumes to be 1.0 for every compound is 1307.5x (n = 33 rat brain
  slices, 3.2-4184 mL/g).

So a report no longer ranks anything. It takes ONE drug. What survives, and what
this document now says, is stated in the ``map_is_the_output`` limit and repeated
in every disclosure on the page:

    Within a single compound on a single tissue, the RELATIVE exposure between
    regions is the defensible output. The absolute level is not, and comparison
    BETWEEN compounds is not.

How it is built
---------------
Nothing here simulates anything. A report is assembled out of RUNS:

* every arm is a normal ``RunSpec`` submitted through the ordinary
  :class:`app.server.jobs.JobManager` queue, so it shares the single worker, the
  reconstruction cache, cancellation, on-disk persistence and ``/api/runs``;
* the one drug contributes TWO arms on the same reconstruction — the real run,
  and the ``resistance_off=true`` counterfactual whose kill is the denominator of
  ``rho``. That counterfactual is PER-COMPOUND: it compares this compound with
  itself on this tissue, so nothing about it depended on the withdrawn ranking
  and it is what MEASURES the resistance term rather than asserting it;
* when both arms are ``done`` the report is assembled by
  :func:`attribute_compound` below.

Because both arms share ``(section, dx, vascular_mode, oxygen_mode,
state_mode)``, the expensive reconstruction is built once and the second arm is a
cache hit. That is the whole reason the report goes through the queue instead of
around it.

Why the attribution is NOT ``tmesim.report.attribution``
--------------------------------------------------------
That module computes every term as a deviation from the PANEL MEAN
(``delivery_i = log10 C_i - mean_j log10 C_j`` and so on). **At N = 1 every one of
those terms is identically zero**, so ported unchanged it would print three
zeros and call them an explanation. It is a ranking artefact, not a per-compound
reading, and it could not be carried over. It is untouched and still used by
``experiments/patient_report.py``; :func:`attribute_compound` here is a new,
absolute, single-compound decomposition. See its docstring for why the
replacement is strictly better than the thing it replaces.

What it refuses to do
---------------------
* It will not run in ``sensitivity_mode='determinant'``: on a section that
  measures some compounds' determinants and not others, an estimate of intrinsic
  resistance for THIS compound would be a measurement-coverage artifact rather
  than a measurement (``results/real_drug_panel.md``).
* It will not present an absolute kill without the interval and the state-route
  span that ``results/kill_calibration.md`` and ``results/cell_state_impact.md``
  measured for it.
* It will not invent ``C_plasma``. If the run's ``plasma_mean_uM`` is missing or
  non-finite the decomposition REFUSES with a stated reason and prints no
  numbers, rather than substituting one.
* It carries :data:`NOT_FOR_CLINICAL_USE` at the top AND the bottom of every
  rendering, and the limits block is the first section of the document, never a
  footnote.
"""
from __future__ import annotations

import html as _html
import json
import math
import os
import threading
import time
from typing import Dict, List, Optional, Sequence

from .cache import (
    PRIVATE_FILE_MODE,
    describe_section,
    lookup_section,
    own_only,
    private_makedirs,
    redact_exception,
    repo_root,
    section_genes,
    state_route_support,
)
from .delivery_law import delivery_law
from .jobs import (PDE_DT_S, JobManager, claim_set_stamp, claim_set_status,
                   jsonable)

__all__ = [
    "MAP_CLAIM",
    "NOT_FOR_CLINICAL_USE",
    "REPORT_SENSITIVITY_MODE",
    "TERMS",
    "TERM_LABEL",
    "Report",
    "ReportManager",
    "attribute_compound",
    "build_limits",
    "render_html",
    "reports_dir",
]

#: The claim this document makes, stated once and quoted verbatim by the
#: ``map_is_the_output`` limit, by the hero, and by every rewritten disclosure
#: that used to talk about a ranking. One sentence, three clauses, and the two
#: negative clauses are not optional trimming.
MAP_CLAIM = (
    "Within a single compound on a single tissue, the RELATIVE exposure between "
    "regions is the defensible output. The absolute level is not, and comparison "
    "BETWEEN compounds is not."
)

#: Forced, not defaulted, and it stays forced now that there is no ranking to be
#: fair about. On a section carrying NO measured resistance determinant for this
#: compound, per-cell sensitivity 1.0 is the honest reading: it asserts no
#: intrinsic resistance, which makes every kill on the page an explicit
#: no-intrinsic-resistance UPPER BOUND — a bound a reader can reason about. The
#: alternative, 'determinant', would scale this compound by whichever determinant
#: the panel happens to carry and leave it at 1.0 where the panel does not, so
#: the number would move with MEASUREMENT COVERAGE rather than with
#: pharmacology. An upper bound that says so is defensible; an estimate that
#: silently encodes what the panel could see is not (results/real_drug_panel.md).
REPORT_SENSITIVITY_MODE = "uniform"

NOT_FOR_CLINICAL_USE = (
    "NOT FOR CLINICAL USE. This is a research prototype. Every number in this "
    "report is a PREDICTION of a physical model, not a measurement of a "
    "patient. Nothing here has been validated against a treatment outcome. Do "
    "not use it to select, rank, exclude, or dose any therapy for any patient."
)

#: Terminal report states.
TERMINAL = ("done", "failed", "cancelled")


def reports_dir() -> str:
    """Absolute path to ``app/reports/`` (created on demand)."""
    return private_makedirs(os.path.join(repo_root(), "app", "reports"))


# --------------------------------------------------------------------------- #
# The limits block — mandatory, prominent, and derived from this repo's own
# measurements rather than asserted.
# --------------------------------------------------------------------------- #
def build_limits(section_id: str, spec: Dict, drug: Optional[str]) -> List[Dict]:
    """The limits that travel with EVERY report, in the order they are shown.

    ``drug`` is the ONE compound this document is about. It used to be a panel;
    the ranked panel was withdrawn on 2026-08-21 (see the module docstring), so
    nothing here compares one compound against another any more.

    Each entry is ``{id, severity, title, body, source, applies}``. ``applies``
    is a property of this section and this spec, never a global switch: an item
    that does not apply is still shown, saying why, because a caveat that
    disappears cannot be audited.
    """
    support = state_route_support(section_id)
    genes_known = bool(support.get("genes_known"))
    hypoxia_on_panel = bool(support.get("hypoxia_present"))
    section = lookup_section(section_id) or {}
    dx = spec.get("dx_um")
    limits: List[Dict] = []

    # 1 -------------------------------------------------------------------- #
    limits.append({
        "id": "not_for_clinical_use",
        "severity": "critical",
        "applies": True,
        "title": "NOT FOR CLINICAL USE",
        "source": None,
        "body": NOT_FOR_CLINICAL_USE,
    })

    # 2 -------------------------------------------------------------------- #
    # THE claim of this document, and the reason the ranked panel it replaced was
    # withdrawn on 2026-08-21. Every other disclosure on this page that used to
    # say "the ranking" now says this instead; none of them was shortened.
    limits.append({
        "id": "map_is_the_output",
        "severity": "critical",
        "applies": True,
        "title": ("The relative exposure MAP is the defensible output. The "
                  "absolute level is not, and a comparison BETWEEN compounds is "
                  "not."),
        "source": ("results/cell_state_impact.md §3.2, §3.3, §5; "
                   "results/pk_ranking_robustness.md sec. 0, sec. 1.2; "
                   "results/capacity_ranking_robustness.md sec. 0, sec. 2b; "
                   "results/descriptor_falsification.md"),
        "body": (
            f"{MAP_CLAIM} "
            "WHY THE ABSOLUTE LEVEL IS NOT. With tissue, drug, transport and RNG "
            "held identical, kill moves by up to 0.3252 absolute (1.77x "
            "relative) purely with the choice of cell-state route "
            "('oxygen', 'vessel_distance', 'markers'). Read where the drug goes. "
            "Do not read how much of the tumour it kills. "
            "WHY A CROSS-COMPOUND COMPARISON IS NOT, AND WHY THIS PAGE NO LONGER "
            "CARRIES ONE. This report used to rank a panel and to name that "
            "ordering as its defensible output. THAT CLAIM IS WITHDRAWN as of "
            "2026-08-21. Per-compound plasma exposure across the compounds the "
            "GBM panel ranked spans 113.1x, while the published GBM order first "
            "breaks at a per-compound spread of only 1.53x — and one FDA label "
            "table for a single compound (irinotecan: Cmax 1,660 ng/mL at "
            "125 mg/m2, 3,392 ng/mL at 340 mg/m2, a 2.04x span) already contains "
            "more spread than the 1.961x a swap in that order needs. "
            "Independently, under this repo's own structural tissue-capacity "
            "rule the breast top-1 CHANGES (vinorelbine -> doxorubicin), and no "
            "sampled per-compound capacity arm reached decision-grade on any "
            "panel; the engine assumes that capacity term is exactly 1.0 for "
            "every compound, and its measured spread is 1307.5x (n = 33, "
            "3.2-4184 mL/g). "
            "WHAT THE RELATIVE READING IS PROTECTED AGAINST, AT THE STRENGTH "
            "EACH WAS ACTUALLY MEASURED — THE TWO ARE NOT THE SAME STRENGTH AND "
            "THE DIFFERENCE MATTERS. (1) PLASMA / DOSE SCALE: EXACT for the "
            "finite-window AUC field. The transport PDE is linear in C and "
            "starts from C = 0, so the AUC field is exactly proportional to the "
            "plasma level; measured at the shipping configuration (real Xenium "
            "GBM section, dx = 12 um, dt = 10 s, convection ON, 6 h window) over "
            "plasma factors 0.005 to 26 — a 5,200x span — the worst relative "
            "deviation of the field from an exact rescale anywhere was 1.41e-14, "
            "i.e. machine precision. A per-compound plasma factor therefore "
            "moves every point of this map together and cancels exactly from a "
            "ratio between two points. (2) TISSUE CAPACITY: EXACT FOR A "
            "CONVERGED GRADIENT SHAPE ONLY, verified to 1e-9 — and NOT "
            "unconditionally. In the source's own words, 'steady-state and "
            "equilibrium-referenced outputs — penetration length, "
            "interstitial:plasma ratio, converged gradient shape — are EXACTLY "
            "invariant to tissue capacity ... The finite-window kill ranking is "
            "not.' A capacity term R enters as dC/dt = (1/R)[...], a pure "
            "rescale of TIME, so over a FINITE exposure window "
            "AUC_R(x,T) = R * AUC_1(x,T/R): a region of this map that has NOT "
            "converged within the window carries a residual, and the invariance "
            "holds there only to the extent the field has converged. This report "
            "runs a finite window. That is not asserted to be converged anywhere. "
            "WHAT NEITHER INVARIANCE PROTECTS AGAINST. Error in the transport "
            "LAW itself. A wrong effective diffusivity changes the SHAPE of the "
            "gradient, not merely its level, so it is not a common factor and "
            "nothing above divides it out — and that law's descriptor inputs "
            "have never been scored "
            "(results/descriptor_falsification.md). "
            "THIS IS A LIST OF THE AXES THAT WERE SWEPT, NOT A COVERAGE CLAIM. "
            "Exactly two perturbations were measured: plasma/dose scale, and "
            "tissue capacity. An axis that does not appear in that list of two "
            "is UNTESTED, not covered; silence about it here is silence, not a "
            "pass."
        ),
    })

    # 2b ------------------------------------------------------------------- #
    # Limit 2 names two invariances the map DOES have (plasma scale exactly;
    # tissue capacity for a converged shape). This one sits immediately beside it
    # because those invariances are common-mode factors, and the delivery LAW is
    # not: it sets the source term at the vessel wall, so changing it changes the
    # SHAPE of the gradient this map is entirely made of, not merely its level.
    # The two sentences have to sit next to each other or the first reads as a
    # broader guarantee than it is. The law is read from PDEConfig at call time
    # (app/server/delivery_law.py), so this limit cannot describe a law the
    # engine is no longer using.
    law = delivery_law()
    limits.append({
        "id": "transvascular_delivery_law",
        "severity": "critical",
        "applies": True,
        "title": (law["title"] + " — the MAP above is LAW-DEPENDENT, and the law "
                  "is the one thing limit 2's invariances do NOT divide out"),
        "source": law["reference"],
        "body": (
            law["caveat"]
            + " HOW THIS BEARS ON THE MAP RATHER THAN ON A RANKING. This "
            "document no longer ranks compounds, so the law's power to reorder a "
            "panel is no longer the relevant harm; the relevant harm is narrower "
            "and closer to home. The source term is what injects drug at the "
            "vessel wall, so the choice of law sets the boundary condition of "
            "the very gradient this page reads regions off. It is therefore NOT "
            "a factor that multiplies every point of the map together, and "
            "neither of the invariances in limit 2 removes it: it can move the "
            "relative exposure between two regions of THIS tissue for THIS "
            "compound, which is precisely the quantity this document offers as "
            "defensible. Read every falloff, every reach ratio and every region "
            "comparison below as conditional on the law named in this limit's "
            "title."
        ),
    })

    # Limit 2c. The two limits above condition the map on things that are
    # uncertain. This one names what the transport model is measured to get
    # WRONG, and it has to sit with them rather than in the numerics block at
    # the bottom: D_eff is the single input that sets how far this compound
    # penetrates from a vessel, which is the shape of the map itself.
    # Scored blind, unfitted, on 71 compounds.
    limits.append({
        "id": "diffusivity_two_route_partial",
        "severity": "critical",
        "applies": True,
        "title": ("Transport reproduces about a QUARTER of the measured spread "
                  "between compounds — similar-mass analogues are separated "
                  "weakly, not reliably"),
        "source": ("results/two_route_promotion.md; "
                   "results/diffusivity_falsification.md; "
                   "results/two_route_transport.md (Pruijn 2008, "
                   "doi:10.1071/CH08240)"),
        "body": (
            "WHAT CHANGED, BECAUSE THIS LIMIT USED TO SAY SOMETHING STRONGER. "
            "Until results/two_route_promotion.md the D_eff handed to the "
            "transport solve was derived from molecular weight and the tissue "
            "ECM fraction and NOTHING ELSE, and that model was FALSIFIED: "
            "scored blind against 71 compounds measured by one method in one "
            "cell line it reached R^2 = 0.022 against a 29.7x measured spread, "
            "and at IDENTICAL molecular weight it was flat BY CONSTRUCTION -- "
            "exactly 1.00x -- where the measurement spans 3.7-3.8x. The engine "
            "now applies a two-route cell term that brings logP and "
            "hydrogen-bond donor count into D_eff through a membrane "
            "permeability. On the same 71 compounds, blind and unfitted, it "
            "reaches R^2 = 0.718 with Spearman +0.821, 99 % of compounds within "
            "3x of measurement, and 24 of 31 exact-mass pairs ordered correctly "
            "against the old model's 0 of 31. "
            "WHAT IS STILL NOT TRUE, AND IT IS WHY THIS LIMIT IS STILL "
            "CRITICAL. The model reproduces 8.3x of the measured 29.7x spread "
            "-- roughly a QUARTER -- so it UNDER-SEPARATES a congeneric series "
            "systematically rather than separating it correctly. The entire "
            "test is ONE cell line. And absolute transport remains wrong: the "
            "same drug's compiled D_eff misses its measured multicellular-layer "
            "band by 1.5x (49.3 um^2/s against a measured 74-130). "
            "CONSEQUENCE FOR THIS PAGE, RESTATED FOR A SINGLE-COMPOUND MAP. "
            "This limit used to end by telling a reader to treat the ORDER of a "
            "similar-mass series as a hypothesis rather than a measurement. "
            "There is no order on this page any more, and the consequence is "
            "narrower but it is not smaller: D_eff is the one compound-specific "
            "input that sets how far THIS compound travels from a vessel before "
            "it is cleared, so it fixes the length scale of the gradient every "
            "region comparison below is read off. The two-route model that "
            "produces it is scored, blind and unfitted, at R^2 = 0.718 across "
            "compounds, but the SAME compiled D_eff misses its own measured "
            "multicellular-layer band by 1.5x (49.3 um^2/s against a measured "
            "74-130), and the entire scoring set is ONE cell line. Read this "
            "compound's absolute penetration length L and absolute D_eff as "
            "carrying at least that 1.5x, and read the SHAPE of the falloff as "
            "conditional on a length scale that has never been measured on this "
            "tissue, this compound, or any 3-D specimen. This is exactly the "
            "error limit 2 names as the one its plasma and capacity invariances "
            "do NOT divide out: it changes the shape of the gradient, not just "
            "its level."
        ),
    })

    # 3 -------------------------------------------------------------------- #
    if genes_known and not hypoxia_on_panel:
        applies, body = True, (
            f"This section's panel carries NONE of the hypoxia-response "
            f"transcripts {support['hypoxia_missing']}, so hypoxia is "
            "UNMEASURED here — assigned to no cell, which is NOT the same "
            "statement as 'computed and found to be zero'. On the 319-gene "
            "flagship panel the proliferation markers MKI67 and TOP2A fire in "
            "only 1.42% of tumour cells against a published Ki-67 index of "
            "19-40%, i.e. 13x below the lower bound, so the per-cell "
            "proliferation call is refused there too and a neighbourhood "
            "enrichment fallback fires whose LEVEL moves 8.79x with the window "
            "radius. NO state route clears the bar on such a panel. An absolute "
            "kill fraction MUST NOT be quoted from this section under any "
            "state_mode."
        )
    elif genes_known:
        applies, body = False, (
            f"This section's panel does carry hypoxia-response transcripts "
            f"{support['hypoxia_present']}, so the disqualification that "
            "applies to the 319-gene flagship panel (no hypoxia gene at all, "
            "and MKI67/TOP2A firing in 1.42% of tumour cells against a "
            "published 19-40% Ki-67 index) does not apply verbatim here. The "
            "hypoxia band remains an uncalibrated ENRICHMENT label, so the "
            "absolute level is still not earned."
        )
    else:
        applies, body = True, (
            "This section's panel gene list could not be enumerated from the "
            "cached bytes, so whether it can support a cell-state call is "
            "UNKNOWN. Treat the absolute kill as unsupported until it is."
        )
    limits.append({
        "id": "cell_state_support",
        "severity": "critical" if applies else "info",
        "applies": applies,
        "title": "Cell state — can this panel support an absolute kill at all?",
        "source": "results/cell_state_impact.md §1.1, §3.1, §5",
        "body": body,
    })

    # 4 -------------------------------------------------------------------- #
    limits.append({
        "id": "kill_is_an_interval",
        "severity": "critical",
        "applies": True,
        "title": "Where absolute kill IS defensible, it is an INTERVAL, and a wide one",
        "source": "results/kill_calibration.md §5.1",
        "body": (
            "Propagating both open terms — the measured inter-culture potency "
            "spread and the measured cell-state route span — gives temozolomide "
            "a 3-D kill of 0.0048-0.1442 around a point estimate of 0.0326: a "
            "30x span. Cisplatin lands at 0.2667-0.8852 and doxorubicin at "
            "0.1666-0.8065. A number with a 30x error bar is not a number "
            "anyone doses from. Every single-valued kill printed in this report "
            "is the centre of a bar of that order."
        ),
    })

    # 5 -------------------------------------------------------------------- #
    limits.append({
        "id": "uniform_upper_bound",
        "severity": "critical",
        "applies": True,
        "title": f"sensitivity_mode='{REPORT_SENSITIVITY_MODE}' — every kill here is an UPPER BOUND",
        "source": "reconstruct.real_pipeline.apply_drug_sensitivity; results/real_drug_panel.md",
        "body": (
            f"{drug or 'This compound'} carries per-cell sensitivity 1.0 "
            "everywhere, so NO intrinsic resistance is asserted for it anywhere "
            "on this section. That is deliberate and it is FORCED, not "
            "defaulted. THE REASON IT IS FORCED HAS CHANGED, AND THE NEW REASON "
            "IS THE ONE THAT OUTLIVES THE RANKING. It used to be justified as "
            "'the only mode in which a cross-drug comparison is fair, because an "
            "identical per-cell factor cancels out of it' — an argument about a "
            "ranking, and this document no longer ranks. The reason it survives "
            "is simpler and it applies to a single compound on a single tissue: "
            "on a section that carries NO measured resistance determinant for "
            "this compound, a no-resistance UPPER BOUND is the honest reading, "
            "and an estimate would not be. A bound that declares itself a bound "
            "can be reasoned about — the true kill is at or below it, and the "
            "reader knows in which direction they are wrong. The alternative "
            "'determinant' mode produces no such guarantee: it scales a compound "
            "by whichever determinant this panel happens to carry (doxorubicin "
            "by TOP2A, cisplatin / lomustine / methotrexate by an MKI67 proxy, "
            "to ~0.585) and leaves a compound whose determinant is ABSENT from "
            "the panel at 1.0 — temozolomide, whose determinant MGMT is on "
            "neither Xenium panel, keeps 1.0 and is thereby flattered ~1.7x for "
            "a MEASUREMENT-COVERAGE reason and not a pharmacological one. A "
            "number that moves with what the panel could see, while looking like "
            "a number about the drug, is worse than a bound. The price of the "
            "bound is stated plainly: every ABSOLUTE kill below is a "
            "no-intrinsic-resistance UPPER BOUND, not a prediction of kill."
        ),
    })

    # 6 -------------------------------------------------------------------- #
    oxy_mode = spec.get("oxygen_mode") or "pde"
    # This limit previously read "The oxygen PDE is FALSIFIED as value-add",
    # asserting the solve was beaten by a distance transform. That was an
    # INVERTED INEQUALITY -- every predictor in the source table is sign-aligned
    # to oxygen, so the stronger one is the one with the larger |rho|, and the
    # solve has the larger |rho|. See results/oxygen_demotion.md, which is now
    # the record of the correction. The limit is NOT deleted: the honest reading
    # is still a warning, because the winning association is weak.
    limits.append({
        "id": "oxygen_pde_falsified",
        "severity": "warn",
        "applies": True,
        "title": "The oxygen field is WEAK — the best of several weak predictors",
        "source": "results/oxygen_demotion.md; results/state_validation.md §1.2",
        "body": (
            "Against measured hypoxia markers in 816,769 real cells the "
            "reaction-diffusion oxygen solve scores rho = -0.1523, "
            "out-predicting a one-line distance transform built from its own "
            "source cells (-0.1204) and its own undiffused source term "
            "(-0.1111). It is the default ON that evidence. But read the "
            "magnitude, not just the ordering: |rho| ~ 0.15 is a WEAK "
            "association, measured within one section rather than out of "
            "sample, and carried by VEGFA alone (HIF1A contributes +0.0140). "
            "The defensible statement is 'the best of several weak predictors', "
            "never 'the oxygen solve is validated'. Separately and still true: "
            "OxygenConfig.consumption has never been fitted to a measured "
            "oxygen tension, so the solve's effective diffusion length is "
            "several times the published ~150 um oxygen-diffusion distance — a "
            "rank correlation cannot detect that, being invariant to monotone "
            "rescaling, so absolute oxygen levels and band edges remain "
            f"uncalibrated. This report ran oxygen_mode='{oxy_mode}'. The "
            "oxygen route's downstream effect on predicted kill is small and "
            "therapeutic coverage is oxygen-invariant by construction — so the "
            "EXPOSURE MAP that is this document's output does not rest on the "
            "oxygen field at all, and no absolute hypoxia claim should rest on "
            "it either. Two things on this page DO read the oxygen field and "
            "must be read with this limit beside them: the RESISTANCE decade, "
            "whose rho is measured by re-simulating this compound with the "
            "hypoxic and quiescent kill penalties switched off, and the 'viable' "
            "tumour-cell band the reach curve is taken over, which is the "
            "engine's own oxygen banding. Both inherit the |rho| ~ 0.15 "
            "association above, and neither is a measurement of this specimen's "
            "oxygen tension. AN EARLIER VERSION OF THIS LIMIT "
            "DECLARED THE SOLVE FALSIFIED AS VALUE-ADD; that claim was "
            "backwards by its own quoted numbers and has been withdrawn."
        ),
    })

    # 7 -------------------------------------------------------------------- #
    limits.append({
        "id": "no_outcome_label",
        "severity": "critical",
        "applies": True,
        "title": "Nothing here is validated against a patient outcome",
        "source": "app/server/cache.py section registry (has_outcome_label)",
        "body": (
            f"The specimen {section_id!r} is an untreated resection "
            f"({section.get('accession') or 'open dataset'}). There is no "
            "paired, measured drug response for any drug on it, and there is NO "
            "outcome-labelled section anywhere in this repository. THE MAP BELOW "
            "HAS THEREFORE NEVER BEEN SCORED: not its shape, not its falloff "
            "with distance from vessel, not the fraction of tumour it says the "
            "drug never reached. Its accuracy is not 'good', not 'bad', it is "
            "UNKNOWN — and note that withdrawing the ranking did not weaken this "
            "limit by one word, because the thing that was missing was never a "
            "second compound, it was a MEASUREMENT to compare against. Potency "
            "additionally comes from GDSC "
            "cell-line viability screens or curated literature — dissociated "
            "cells in monolayer, carrying none of this tumour's 3-D "
            "architecture, stromal/immune context or heterogeneity."
        ),
    })

    # 8 -------------------------------------------------------------------- #
    # This limit used to assert, at the renderer's LOWEST severity, that
    # "kill_fraction is grid-converged for dx <= 20 um" and that the penetration
    # index is over-estimated on coarse grids. Both halves are retracted by the
    # sources below, so the limit shipped a reassurance the documentation had
    # already withdrawn. It is NOT deleted -- the dx a report ran at is exactly
    # the kind of thing a limits block exists to state -- but it now reports
    # what was actually measured, and it can no longer be filed as "info".
    if isinstance(dx, (int, float)):
        dx_f = float(dx)
        dx_body = (
            f"This report ran at dx = {dx_f:g} um. THERE IS NO dx AT WHICH THIS "
            "PIPELINE HAS BEEN SHOWN TO BE GRID-CONVERGED, and an earlier "
            "version of this limit said otherwise. It asserted 'kill_fraction "
            "is grid-converged for dx <= 20 um' and that the penetration index "
            "is 'over-estimated on coarse grids (0.285 at dx = 60 um falling to "
            "0.187 at dx = 10 um)'. BOTH HALVES ARE WITHDRAWN. (a) The kill "
            "plateau was measured on ONE Xenium GBM section for ONE metric and "
            "does not transfer: on breast, docetaxel scores 0.1484 at the "
            "published dx = 12 um against a Richardson dx -> 0 estimate of "
            "0.1537, and its verdict FLIPS FAIL -> PASS between dx = 12 and "
            "dx = 16 um; on the GBM panel the scored ratio is NON-MONOTONE in "
            "dx for all three drugs tested, which is not an asymptotic regime, "
            "so Richardson and power-law fits are being applied outside their "
            "assumptions (observed orders p = 0.96 and p = -3.33, GCI 5.44% and "
            "6.42% on the finest pairs, both power-law fits flagged "
            "DEGENERATE). (b) The penetration ladder RISES with refinement — "
            "0.1110 at dx = 60 um to 0.3062 at dx = 10 um — so a coarse mesh "
            "UNDER-states the index, the opposite direction to what was "
            "printed here, and the quoted 0.285 / 0.187 pair appears nowhere in "
            "that column. What survives, and is the reason this limit is a "
            "warning rather than a footnote: dx is the LARGEST single driver of "
            "the penetration index measured on that section (|elasticity| "
            "0.747, 3.1x the next knob), and it reaches the answer through the "
            "RECONSTRUCTION — the oxygen solve, the cell-state banding and the "
            "IFP/velocity solve all live on the grid — so a converged transport "
            "step does not imply a converged prediction. Read every number on "
            "this page as carrying an unquantified mesh term. This is not a "
            "claim that the solver is broken: mass is conserved to ~4e-14 on "
            "the real irregular mask and the prediction is dt-converged to "
            "<1e-2 for dt <~ 20 s, and neither of those is affected by any of "
            "the above."
        )
        if dx_f > 18.0:
            dx_body += (
                f" This run is COARSER STILL: at dx = {dx_f:g} um it sits on "
                "the part of both ladders where the coarse rungs (24 and 30 um) "
                "understate the scored ratio by the largest margin measured — "
                "up to 11.7% in a single refinement step."
            )
        limits.append({
            "id": "resolution",
            "severity": "warn",
            "applies": True,
            "title": "Resolution — the grid-convergence claim is WITHDRAWN",
            "source": ("results/numerical_rigor.md CORRECTION + §2; "
                       "results/grid_convergence_breast.md; "
                       "results/grid_convergence_gbm.md"),
            "body": dx_body,
        })

    # 9 -------------------------------------------------------------------- #
    limits.append({
        "id": "one_anchor",
        "severity": "warn",
        "applies": True,
        "title": "One measurement sets the absolute delivery scale",
        "source": "validation.calibration.calibrate_alpha (PMID 16361566)",
        "body": (
            "The perivascular-sink constant alpha is fit to a SINGLE "
            "doxorubicin penetration measurement (a 40-50 um half-distance, "
            "converted to an e-folding length L = d_half/ln2 ~ 65 um) and "
            "revalidated out-of-sample on exactly one other drug, "
            "mitoxantrone. This compound's absolute penetration rides on that "
            "one anchor, and the underlying binding proxy is a structure-based "
            "ORDERING of compounds, not a first-principles rate — so a compound "
            "that proxy places wrongly (active transport, efflux, ion trapping) "
            "gets the wrong alpha and therefore the wrong perivascular sink. "
            "THIS IS A MAP-SHAPE RISK, NOT ONLY A LEVEL RISK, WHICH IS WHY IT "
            "SURVIVES THE WITHDRAWAL OF THE RANKING UNCHANGED IN SEVERITY. "
            "alpha sets how steeply exposure falls away from a vessel, so it "
            "does not multiply the map by a constant — it re-shapes it, and "
            "limit 2's plasma and capacity invariances do not divide it out. "
            f"Read the falloff and the reach curve for "
            f"{drug or 'this compound'} as resting on a single doxorubicin "
            "measurement plus one out-of-sample compound."
        ),
    })

    # 10 ------------------------------------------------------------------- #
    limits.append({
        "id": "pk_simplification",
        "severity": "info",
        "applies": True,
        "title": "Systemic PK is the documented Stage-0 simplification",
        "source": "tmesim/sim/engine.py",
        "body": (
            "A single time-averaged sustained plasma level drives the transport "
            f"solve at dt = {PDE_DT_S:g} s. It preserves total delivered drug, "
            "and — this is the part that matters for a map rather than for the "
            "ranking this sentence used to be about — the plasma level enters "
            "the solve as ONE MULTIPLICATIVE FACTOR on the source term, so it "
            "scales every point of the exposure field together and cancels "
            "exactly from a ratio between two points (limit 2, invariance 1: "
            "worst relative deviation 1.41e-14 over plasma factors 0.005 to 26). "
            "What it does NOT preserve is anything absolute or anything "
            "time-dependent: it is not a time-resolved infusion, the absolute "
            "exposure level and therefore the IC50 coverage move with it "
            "one-for-one, and no schedule optimisation, no C_max reading and no "
            "dosing inference should be read off it."
        ),
    })
    return limits


# --------------------------------------------------------------------------- #
# Per-drug provenance
# --------------------------------------------------------------------------- #
def _mgmt_status(drug: str, spec: Dict, prov: Dict, section_id: str) -> Dict:
    """MGMT modelling status for one drug on this section.

    Says ``unmodelled - MGMT not on panel`` where that is the truth, rather than
    silently omitting a resistance determinant that is known to matter.
    """
    from tmesim.molecular.mgmt import MGMT_CALIBRATED_DRUGS, is_mgmt_calibrated

    if not is_mgmt_calibrated(drug):
        return {
            "state": "not applicable",
            "detail": (
                f"MGMT is not a calibrated resistance determinant for {drug}; "
                f"only {list(MGMT_CALIBRATED_DRUGS)} is "
                "(tmesim.molecular.mgmt.MGMT_CALIBRATED_DRUGS)."
            ),
            "value": None,
        }
    if prov.get("mgmt_applied"):
        return {
            "state": "modelled",
            "detail": (
                f"MGMT sensitivity scalar {prov.get('mgmt_scalar')} "
                f"(95% CI {prov.get('mgmt_scalar_ci95')}, calibration n = "
                f"{prov.get('mgmt_calibration_n')}, model "
                f"{prov.get('mgmt_model')!r}) applied from the supplied value "
                f"{prov.get('mgmt_input')!r}. {prov.get('mgmt_note') or ''}"
            ).strip(),
            "value": prov.get("mgmt_input"),
        }

    genes = section_genes(section_id)
    on_panel = "MGMT" in set(genes)
    if genes and not on_panel:
        return {
            "state": "unmodelled - MGMT not on panel",
            "detail": (
                f"MGMT is the resistance determinant for {drug} and it is NOT "
                f"on this section's {len(genes)}-gene panel, so no per-patient "
                "MGMT term could be measured or applied. MGMT is additionally a "
                "promoter-METHYLATION state that transcript counts cannot "
                "report even where the gene is on the panel, and it has no "
                "panel surrogate (0 of 396 candidates pass FDR). This drug's "
                "kill is therefore an upper bound with a KNOWN, NAMED missing "
                "resistance term."
            ),
            "value": None,
        }
    if not genes:
        return {
            "state": "unmodelled - panel gene list unknown",
            "detail": (
                "This section's gene list could not be enumerated from the "
                f"cached bytes, so whether MGMT — the determinant for {drug} — "
                "is measurable here is unknown. No MGMT term was applied."
            ),
            "value": None,
        }
    return {
        "state": "unmodelled - no MGMT value supplied",
        "detail": (
            f"MGMT is on this panel but no MGMT value was supplied in the "
            f"RunSpec, so no per-patient term was applied for {drug}. Supply "
            "spec_overrides.mgmt (log2 CPM, or 'methylated'/'unmethylated') to "
            "model it."
        ),
        "value": None,
    }


def _ic50_provenance(prov: Dict) -> Dict:
    conf = prov.get("ic50_confidence")
    n = prov.get("ic50_n_cell_lines")
    if conf == "measured":
        label = f"GDSC measured (n={n} glioma cell lines)"
        kind = "measured"
    elif conf == "literature":
        label = "curated literature"
        kind = "assumed"
    else:
        label = "flagged default (LOW confidence)"
        kind = "assumed"
    return {
        "ic50_uM": prov.get("ic50_uM"),
        "confidence": conf,
        "source": label,
        "n_cell_lines": n,
        "note": prov.get("ic50_note"),
        "provenance_kind": kind,
    }


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #
def _num(v, default=float("nan")) -> float:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    return f if math.isfinite(f) else default


def _compound_row(drug: str, real: Dict, cf: Dict, mean_ecm: float,
                  section_id: str) -> Dict:
    """The ONE compound's row, built only out of numbers the engine produced."""
    from tmesim.physics.diffusion import ALPHA_TORTUOSITY

    prov = real.get("provenance") or {}
    metrics = real.get("metrics") or {}
    extra = real.get("extra_metrics") or {}
    dp = real.get("drug_params") or {}
    cf_metrics = (cf or {}).get("metrics") or {}

    diffusivity = _num(dp.get("diffusivity_um2_s"))
    d_eff = (diffusivity * (1.0 - ALPHA_TORTUOSITY * min(max(mean_ecm, 0.0), 1.0))
             if math.isfinite(diffusivity) else float("nan"))

    return {
        "name": drug,
        # --- the inputs :func:`attribute_compound` reads --------------------- #
        "tumor_auc": _num(metrics.get("mean_tumor_auc"), 0.0),
        "ic50": _num(prov.get("ic50_uM"), float("nan")),
        "kill": _num(metrics.get("kill_fraction"), 0.0),
        "kill_no_resistance": _num(cf_metrics.get("kill_fraction"), 0.0),
        # C_plasma. Read from the run's own extra_metrics, NEVER recomputed and
        # never defaulted: if it is absent the decomposition refuses (see
        # attribute_compound), because inventing the denominator of the DELIVERY
        # term would invent the term.
        "c_plasma_uM": _num(extra.get("plasma_mean_uM")),
        "plasma_peak_uM": _num(extra.get("plasma_peak_uM")),
        # --- reported alongside --------------------------------------------- #
        "L_um": _num(extra.get("penetration_length_analytic_um")),
        "mw": _num(dp.get("molecular_weight")),
        "bbb": _num(dp.get("bbb_penetrance")),
        "diffusivity": diffusivity,
        "d_eff": d_eff,
        "hill": _num(dp.get("hill")),
        "binding": _num(dp.get("binding_1_per_s")),
        "clearance": _num(dp.get("clearance_1_per_s")),
        "pen": _num(extra.get("penetration_index")),
        "mean_auc": _num(extra.get("mean_tissue_auc")),
        "median_tumor_auc": _num(extra.get("median_tumor_auc")),
        # The IC50 coverage LADDER, not one number. 1x IC50 acquired an authority
        # it had not earned by being quoted alone; all three rungs were already
        # computed by every run and only one of them was ever rendered.
        "coverage_0_5x": _num(extra.get("therapeutic_coverage_0_5x")),
        "coverage_1x": _num(metrics.get("therapeutic_coverage_1x")),
        "coverage_2x": _num(extra.get("therapeutic_coverage_2x")),
        "E_p10_uM": _num(extra.get("E_p10_uM")),
        "E_p50_uM": _num(extra.get("E_p50_uM")),
        "E_p90_uM": _num(extra.get("E_p90_uM")),
        "viable_tumor_fraction": _num(extra.get("viable_tumor_fraction")),
        "hypoxic_tumor_fraction": _num(extra.get("hypoxic_tumor_fraction")),
        "penetration_length_fit_um": _num(metrics.get("penetration_length_um")),
        # --- provenance ------------------------------------------------------ #
        "provenance": {
            "run_id": real.get("run_id"),
            "counterfactual_run_id": (cf or {}).get("run_id"),
            "alpha": prov.get("alpha"),
            "alpha_source": prov.get("alpha_source"),
            "ic50": _ic50_provenance(prov),
            "sensitivity_mode": prov.get("sensitivity_mode"),
            "sensitivity_upper_bound": bool(prov.get("sensitivity_upper_bound")),
            "sensitivity_note": prov.get("sensitivity_note"),
            "sensitivity_determinant_gene": prov.get("sensitivity_determinant_gene"),
            "sensitivity_mean": prov.get("sensitivity_mean"),
            "mgmt": _mgmt_status(drug, real.get("spec") or {}, prov, section_id),
            "seed": prov.get("seed"),
            "pde_dt_s": prov.get("pde_dt_s"),
            "state_mode": prov.get("state_mode"),
            "caveats": prov.get("caveats") or [],
        },
    }


# --------------------------------------------------------------------------- #
# The single-compound decomposition
# --------------------------------------------------------------------------- #
#: The three decades, keyed and labelled in the SAME vocabulary the per-cell
#: mechanism attribution uses (``app/web/js/mechanism.js``: delivery / dose /
#: resistance / chance, from ``tmesim.sim.mechanism.LABEL_LIMITER``). The
#: document carries two attributions — this whole-tumour one and the per-cell
#: one — and they must not use different words for the same idea. ``chance`` has
#: no whole-tumour analogue and is deliberately absent rather than faked.
TERMS: Sequence[str] = ("delivery_decades", "dose_decades", "resistance_decades")

TERM_LABEL: Dict[str, str] = {
    "delivery_decades": "DELIVERY",
    "dose_decades": "DOSE",
    "resistance_decades": "RESISTANCE",
}

#: What each decade means, in one sentence, shipped in the payload beside it.
TERM_MEANING: Dict[str, str] = {
    "delivery_decades": (
        "decades LOST getting from blood into this tumour's interstitium: "
        "log10(C_tumor / C_plasma). Negative is the normal case and its "
        "magnitude is the delivery deficit."
    ),
    "dose_decades": (
        "decades of HEADROOM the plasma level has over this compound's own "
        "potency: log10(C_plasma / IC50). This is the term the dose uncertainty "
        "lands in, and it is quarantined here on purpose."
    ),
    "resistance_decades": (
        "decades of kill LOST to this tumour's oxygen-driven resistance state: "
        "log10(rho), where rho = kill / kill(resistance OFF) is MEASURED by the "
        "counterfactual arm, not asserted."
    ),
}

#: Written into the document beside the numbers. Every line of it is checkable
#: against the run it describes.
ATTRIBUTION_METHOD = (
    "C_tumor  = (mean AUC at this tumour's TUMOR-cell positions) / T   [uM]\n"
    "C_plasma = this run's own plasma mean (extra_metrics.plasma_mean_uM)  [uM]\n"
    "drive    = C_tumor / IC50                                         [-]\n"
    "rho      = kill / kill(resistance OFF)                            [-]\n"
    "\n"
    "DELIVERY  = log10(C_tumor  / C_plasma)   decades lost blood -> interstitium\n"
    "DOSE      = log10(C_plasma / IC50)       decades of plasma headroom over IC50\n"
    "RESISTANCE= log10(rho)                   decades of kill lost to O2 resistance\n"
    "\n"
    "  rho is MEASURED by a second simulation of THIS compound on THIS tissue\n"
    "  with hypoxia_resistance=1 and quiescent_resistance=1 -- the\n"
    "  resistance_off counterfactual arm -- never assumed.\n"
    "\n"
    "identity:  DELIVERY + DOSE  ==  log10(drive)\n"
    "\n"
    "Every term is an ABSOLUTE decade count referenced to a quantity from THIS\n"
    "SAME RUN. Nothing is centred on a panel mean."
)

#: Why this decomposition replaced the panel-centred one, in the payload.
ATTRIBUTION_RATIONALE = (
    "This replaced tmesim.report.attribution's PANEL-CENTRED decomposition on "
    "2026-08-21, and it is strictly better for three reasons that are worth "
    "stating rather than assuming. (1) IT EXISTS AT N = 1. The old terms were "
    "deviations from a panel mean (delivery_i = log10 C_i - mean_j log10 C_j and "
    "so on), so on ONE compound every one of them is identically zero: ported "
    "unchanged it would print three zeros and call them an explanation. (2) IT "
    "DOES NOT DEPEND ON THE COMPANY THE COMPOUND KEPT. Every term here is an "
    "absolute, physically interpretable decade count referenced to a quantity "
    "from the SAME run — the run's own plasma mean, the compound's own IC50, the "
    "compound's own resistance-off counterfactual. Panel-centring made a "
    "compound's explanation move when a DIFFERENT compound was added to or "
    "removed from the selection, which was never a defensible property of an "
    "explanation. (3) IT QUARANTINES THE DOSE UNCERTAINTY IN ONE NAMED TERM "
    "INSTEAD OF SMEARING IT ACROSS ALL THREE. DELIVERY is INVARIANT to the "
    "plasma level: C_tumor is exactly linear in the plasma source (worst "
    "relative deviation 1.41e-14 over plasma factors 0.005 to 26, a 5,200x span, "
    "measured at the shipping configuration — results/pk_ranking_robustness.md "
    "sec. 1.2), so rescaling the dose multiplies C_tumor and C_plasma by the "
    "same factor and it cancels from log10(C_tumor / C_plasma) exactly. "
    "RESISTANCE is a ratio of two kills of the same compound at the same dose, "
    "so it does not move either. The entire dose uncertainty therefore lands in "
    "DOSE and nowhere else — read that term as the uncertain one, and read "
    "DELIVERY as the term that survives not knowing the dose. Saying that out "
    "loud is the honest reading and it is a genuine improvement over a "
    "decomposition in which the same uncertainty was spread invisibly over "
    "three terms."
)


def attribute_compound(row: Dict, window_s: float) -> Dict:
    """Decompose ONE compound's standing on ONE tissue into three decades.

    Writes ``c_tumor_uM``, ``drive``, ``rho``, ``resistance_cost`` and the three
    :data:`TERMS` onto ``row`` in place, and returns the method block that ships
    with them.

    The exact identity, asserted in ``tests/test_app_reports.py``::

        DELIVERY + DOSE == log10(drive)

    REFUSAL, NOT DEFAULTING. ``C_plasma`` is read from the run's own
    ``extra_metrics.plasma_mean_uM`` (``app/server/jobs.py``) and is never
    recomputed. If it is missing or non-finite this returns a ``null_reason`` and
    writes NO decades, because the DELIVERY term is a ratio to that number and
    substituting one would invent the answer. The same rule applies to the
    exposure window, the IC50, ``C_tumor`` and the counterfactual kill; each
    refuses on its own terms and says which one failed.

    See :data:`ATTRIBUTION_RATIONALE` for why this is not
    ``tmesim.report.attribution`` (which is untouched, still panel-centred, and
    still what ``experiments/patient_report.py`` calls).
    """
    block = {
        "terms": list(TERMS),
        "term_labels": dict(TERM_LABEL),
        "term_meanings": dict(TERM_MEANING),
        "method": ATTRIBUTION_METHOD,
        "rationale": ATTRIBUTION_RATIONALE,
        "implementation": "app.server.report.attribute_compound",
        "panel_centred_implementation_not_used": (
            "tmesim.report.attribution.attribute — panel-centred, identically "
            "zero at N = 1, unchanged and still used by "
            "experiments/patient_report.py"
        ),
        "null_reason": None,
        "resistance_null_reason": None,
    }
    for key in ("c_tumor_uM", "drive", "rho", "resistance_cost",
                "delivery_decades", "dose_decades", "resistance_decades",
                "log10_drive"):
        row[key] = None

    T = _num(window_s)
    c_plasma = _num(row.get("c_plasma_uM"))
    ic50 = _num(row.get("ic50"))
    auc = _num(row.get("tumor_auc"))

    if not (math.isfinite(T) and T > 0):
        block["null_reason"] = (
            f"the exposure window is {window_s!r}: a time-average of the AUC "
            "field over a non-positive window is undefined, so C_tumor cannot "
            "be formed and no decade is reported."
        )
        return block
    if not (math.isfinite(c_plasma) and c_plasma > 0):
        block["null_reason"] = (
            f"this run reported plasma_mean_uM = {row.get('c_plasma_uM')!r}. "
            "C_plasma is the reference of the DELIVERY term and the numerator "
            "of the DOSE term, and it is read from the run's own extra_metrics "
            "(app/server/jobs.py) rather than recomputed. It is NOT substituted: "
            "a plausible plasma level chosen here would not be a measurement of "
            "this run, it would be an invention, and both decades would inherit "
            "it silently. No decade is reported."
        )
        return block
    if not (math.isfinite(auc) and auc > 0):
        block["null_reason"] = (
            f"mean tumour-cell AUC is {row.get('tumor_auc')!r}, so C_tumor is "
            "not a positive concentration and log10(C_tumor / C_plasma) is "
            "undefined. This is the 'no TUMOR cells in this reconstruction' "
            "case, and it refuses rather than reporting a floor."
        )
        return block
    if not (math.isfinite(ic50) and ic50 > 0):
        block["null_reason"] = (
            f"this compound's IC50 resolved to {row.get('ic50')!r}, so the DOSE "
            "decade log10(C_plasma / IC50) and the drive are undefined. No "
            "potency was substituted."
        )
        return block

    c_tumor = auc / T
    drive = c_tumor / ic50
    delivery = math.log10(c_tumor / c_plasma)
    dose = math.log10(c_plasma / ic50)

    row["c_tumor_uM"] = c_tumor
    row["drive"] = drive
    row["log10_drive"] = math.log10(drive)
    row["delivery_decades"] = delivery
    row["dose_decades"] = dose

    kill = _num(row.get("kill"))
    kill_off = _num(row.get("kill_no_resistance"))
    if not (math.isfinite(kill_off) and kill_off > 0):
        block["resistance_null_reason"] = (
            f"the resistance-off counterfactual arm returned kill = "
            f"{row.get('kill_no_resistance')!r}. rho = kill / kill(resistance "
            "OFF) has no denominator, so the RESISTANCE decade is not reported. "
            "It is MEASURED or it is absent; it is never assumed."
        )
        return block
    if not math.isfinite(kill) or kill <= 0:
        block["resistance_null_reason"] = (
            f"the real arm returned kill = {row.get('kill')!r}. log10(rho) is "
            "undefined at zero kill, so the RESISTANCE decade is not reported "
            "rather than clamped to a floor. DELIVERY and DOSE are unaffected: "
            "neither reads a kill."
        )
        return block

    rho = kill / kill_off
    row["rho"] = rho
    row["resistance_cost"] = 1.0 - rho
    row["resistance_decades"] = math.log10(rho)
    return block


def _mechanism_reading(real: Dict) -> Dict:
    """The PER-CELL answer: how much of this tumour the drug never reached.

    Read straight off ``tmesim.sim.mechanism``'s attribution as the run already
    computed it (``app/server/jobs.py``), which reconciles bit-exactly with the
    engine's own kill fraction or the run fails. Nothing is recomputed here.
    """
    mech = real.get("mechanism")
    if not isinstance(mech, dict):
        return {
            "null_reason": (
                "this run carries no per-cell mechanism attribution: "
                + str(real.get("mechanism_reason")
                      or "no reason was recorded by the run")
                + ". The headline of this document is a per-cell fraction, so it "
                "is not reported at all rather than approximated from a "
                "whole-tumour average."
            ),
            "unreached_fraction_of_tumor": None,
            "unreached_fraction_of_survivors": None,
            "survival_breakdown": None,
            "counts": None,
            "n_tumor": None,
            "n_survived": None,
            "killed_fraction_of_tumor": None,
            "limiter_vocabulary": dict(TERM_LABEL),
        }
    counts = mech.get("counts") or {}
    n_tumor = _num(mech.get("n_tumor"), 0.0)
    unreached = _num(counts.get("SURVIVED_UNREACHED"), float("nan"))
    frac_tumor = (unreached / n_tumor
                  if n_tumor > 0 and math.isfinite(unreached) else None)
    breakdown = mech.get("survival_breakdown") or {}
    return {
        "null_reason": None,
        "unreached_fraction_of_tumor": frac_tumor,
        "unreached_fraction_of_survivors": breakdown.get("delivery"),
        "survival_breakdown": breakdown,
        "survival_breakdown_note": mech.get("survival_breakdown_note"),
        "counts": counts,
        "n_tumor": mech.get("n_tumor"),
        "n_survived": mech.get("n_survived"),
        "killed_fraction_of_tumor": mech.get("killed_fraction_of_tumor"),
        "joint_limited": mech.get("joint_limited"),
        "reconciliation": mech.get("reconciliation"),
        "limiter_vocabulary": dict(TERM_LABEL),
        "definition": (
            "SURVIVED_UNREACHED is a TUMOR cell that survived because the drug "
            "never arrived in quantity at its position — the delivery limiter of "
            "tmesim.sim.mechanism, the same vocabulary the whole-tumour DELIVERY "
            "decade above uses. The attribution reconciles bit-exactly with the "
            "engine's own kill fraction or the run fails."
        ),
    }


def assemble(section_id: str, spec_base: Dict, drug: str,
             arms: Dict[str, Dict]) -> Dict:
    """Build the report body from the two completed arms of ONE compound.

    ``arms`` maps ``{'real': RunResult, 'counterfactual': RunResult}``. There is
    no ranking, no sort and no rank: the ranked panel was withdrawn on
    2026-08-21 (module docstring), and a rank over one row would be decoration.
    """
    real = arms["real"]
    cf = arms.get("counterfactual") or {}
    report = real.get("report") or {}
    mean_ecm = _num(report.get("mean_ecm_density"), 0.0)
    window_s = float(spec_base.get("total_time_s") or 0.0)

    row = _compound_row(drug, real, cf, mean_ecm, section_id)
    attribution = attribute_compound(row, window_s=window_s)

    # The scale-free under-reached curve, CONSUMED from the run, never
    # re-implemented here: app/server/reach.py owns it and app/server/jobs.py
    # attaches it to every drug run. A build whose runs predate it degrades with
    # a stated reason rather than with a silently missing section.
    reach = real.get("reach")
    if isinstance(reach, dict):
        reach_reason = reach.get("null_reason")
    else:
        reach = None
        reach_reason = (
            "this run carries no 'reach' block. The scale-free under-reached "
            "curve is produced by app.server.reach.compute_reach and attached by "
            "app/server/jobs.py; a run made before that landed has none, and no "
            "substitute is computed here — a second implementation of the same "
            "curve is exactly how two numbers that must agree stop agreeing."
        )

    # The union of both arms' own caveats, de-duplicated, order preserved.
    seen, caveats = set(), []
    for arm in ("real", "counterfactual"):
        run = arms.get(arm) or {}
        for c in ((run.get("provenance") or {}).get("caveats") or []):
            if c not in seen:
                seen.add(c)
                caveats.append(c)

    return {
        "compound": jsonable(row),
        "attribution": jsonable(attribution),
        "mechanism": jsonable(_mechanism_reading(real)),
        "reach": jsonable(reach),
        "reach_reason": reach_reason,
        "microenvironment": jsonable(report),
        "caveats": caveats,
        # The claim set the caveats above were written against, stamped beside
        # them for the same reason a run stamps its own: this document is
        # rendered to report.html and KEPT, and read months later. Without the
        # stamp a document assembled today is indistinguishable, to a reader,
        # from one assembled before the claim it repeats was withdrawn.
        "claim_set": claim_set_stamp(),
        "window_s": window_s,
        "map_claim": MAP_CLAIM,
    }


# --------------------------------------------------------------------------- #
# Report objects
# --------------------------------------------------------------------------- #
class Report:
    """One single-compound map document: its spec, its two arms, its status."""

    def __init__(self, report_id: str, section: str, drug: Optional[str],
                 spec_base: Dict, spec_overrides: Dict,
                 label: Optional[str] = None):
        self.report_id = report_id
        self.section = section
        self.drug = drug
        self.spec_base = spec_base
        self.spec_overrides = dict(spec_overrides or {})
        self.label = label
        self.created_at = time.time()
        self.finished_at: Optional[float] = None
        self.status = "queued"
        self.error: Optional[str] = None
        self.warnings: List[str] = []
        #: ``{'real': run_id, 'counterfactual': run_id}`` — flat, because there
        #: is one compound. It was ``drug -> {arm -> run_id}`` while this was a
        #: panel; :meth:`ReportManager._restore` flattens the old shape.
        self.arms: Dict[str, str] = {}
        self.body: Optional[Dict] = None

    # -- views ------------------------------------------------------------- #
    @property
    def run_ids(self) -> List[str]:
        return [rid for rid in (self.arms.get("real"),
                                self.arms.get("counterfactual")) if rid]

    def summary_view(self, progress: Optional[Dict] = None) -> Dict:
        # The list-row headline. It is deliberately NOT a kill and NOT a rank:
        # ``top_ranked`` was removed with the ranking on 2026-08-21. What a row
        # carries instead is the per-cell reading this document leads with — the
        # fraction of tumour cells the compound never reached — which is
        # per-compound, per-cell, and on the real tissue.
        headline = None
        if self.body:
            mech = self.body.get("mechanism") or {}
            headline = {
                "drug": self.drug,
                "unreached_fraction_of_tumor":
                    mech.get("unreached_fraction_of_tumor"),
                "unreached_null_reason": mech.get("null_reason"),
            }
        return {
            "report_id": self.report_id,
            "label": self.label,
            "status": self.status,
            "section": self.section,
            "drug": self.drug,
            "n_runs": len(self.run_ids),
            "created_at": self.created_at,
            "finished_at": self.finished_at,
            "runtime_s": (round(self.finished_at - self.created_at, 3)
                          if self.finished_at else None),
            "progress": progress,
            "headline": headline,
            # Every number in ``headline`` is read off an exposure field whose
            # SHAPE the delivery law sets, and a list row is read without opening
            # the document — so the law travels with the SUMMARY and not only
            # with the full view's limits block. It is re-parented here rather
            # than deleted with the ranking it used to qualify. Read from
            # PDEConfig at call time (app/server/delivery_law.py).
            "delivery_law": delivery_law(),
            # Likewise the claim itself: a summary consumer meets it too.
            "map_claim": MAP_CLAIM,
            "error": self.error,
            "not_for_clinical_use": True,
        }

    def full_view(self, progress: Optional[Dict] = None) -> Dict:
        out = self.summary_view(progress)
        out.update({
            "spec_base": self.spec_base,
            "spec_overrides": self.spec_overrides,
            "arms": self.arms,
            "run_ids": self.run_ids,
            "warnings": self.warnings,
            "sensitivity_mode": REPORT_SENSITIVITY_MODE,
            "not_for_clinical_use_text": NOT_FOR_CLINICAL_USE,
            "limits": build_limits(self.section, self.spec_base, self.drug),
            "section_detail": _safe_section(self.section),
            "generator": "app/server/report.py (POST /api/reports)",
            "attribution_implementation": "app.server.report.attribute_compound",
        })
        if self.body:
            out.update(self.body)
        else:
            out.setdefault("compound", None)
            out.setdefault("attribution", None)
            out.setdefault("mechanism", None)
            out.setdefault("reach", None)
            out.setdefault("reach_reason", None)
            out.setdefault("microenvironment", None)
            out.setdefault("caveats", [])
        return out


def _safe_section(section_id: str) -> Dict:
    try:
        return describe_section(section_id)
    except Exception as exc:  # noqa: BLE001 - a metadata read must not 500
        return {"id": section_id,
                "reason": redact_exception(
                    exc, f"section {section_id!r} could not be described")}


class ReportManager:
    """Submits report arms to the run queue and assembles the result.

    Deliberately has NO worker thread of its own: the arms are ordinary jobs on
    the one existing worker, and assembly is pure arithmetic over their finished
    results. Status is therefore derived from the jobs on read
    (:meth:`refresh`), which makes ``GET /api/reports/{id}`` the poll target with
    no second scheduler to keep consistent.
    """

    def __init__(self, jobs: JobManager, root: Optional[str] = None):
        self.jobs = jobs
        # Sibling of the queue's own directory, so a JobManager pointed at a
        # temporary root (tests) gets a temporary report store too, and the
        # production default lands on ``app/reports/`` next to ``app/runs/``.
        self.root = root or os.path.join(
            os.path.dirname(os.path.abspath(jobs.root)), "reports")
        private_makedirs(self.root)
        self._reports: "Dict[str, Report]" = {}
        self._order: List[str] = []
        self._counter = 0
        self._lock = threading.RLock()
        self._restore()

    # -- persistence -------------------------------------------------------- #
    def _dir(self, report_id: str) -> str:
        return os.path.join(self.root, report_id)

    def _restore(self) -> None:
        try:
            names = sorted(os.listdir(self.root))
        except OSError:
            return
        for name in names:
            path = os.path.join(self.root, name, "report.json")
            if not (name.startswith("rep_") and os.path.exists(path)):
                continue
            try:
                with open(path) as fh:
                    doc = json.load(fh)
            except Exception:  # noqa: BLE001 - skip an unreadable report
                continue
            # A document persisted before 2026-08-21 has ``drugs[]`` and a
            # nested ``arms`` map, because it was a panel. It is restored, not
            # skipped — a stored artefact does not stop existing when the
            # product's claim changes — but it is restored as what it is: the
            # first compound of that panel, with a warning that says so, so a
            # reader is never shown a withdrawn ranking as a current document.
            drug = doc.get("drug")
            legacy = None
            if not drug:
                panel = doc.get("drugs") or []
                drug = panel[0] if panel else None
                if panel:
                    legacy = (
                        f"RESTORED FROM A WITHDRAWN PANEL DOCUMENT. This report "
                        f"was created as a {len(panel)}-compound RANKING "
                        f"({list(panel)}) before the cross-compound ranking was "
                        "withdrawn on 2026-08-21 "
                        "(results/pk_ranking_robustness.md, "
                        "results/capacity_ranking_robustness.md). Its ranking is "
                        "NOT re-rendered and must not be read. Only its first "
                        f"compound ({drug!r}) is carried forward as this "
                        "document's subject; re-run it to obtain a "
                        "single-compound map."
                    )
            rep = Report(doc.get("report_id", name), doc.get("section", ""),
                         drug, doc.get("spec_base") or {},
                         doc.get("spec_overrides") or {}, doc.get("label"))
            rep.status = doc.get("status", "done")
            rep.created_at = doc.get("created_at") or rep.created_at
            rep.finished_at = doc.get("finished_at")
            rep.error = doc.get("error")
            rep.warnings = list(doc.get("warnings") or [])
            if legacy:
                rep.warnings.append(legacy)
            # The ``legacy`` warning above catches a PANEL document, because
            # that is the shape the withdrawal changed. It does not catch a
            # SINGLE-COMPOUND document from the same era: that one has no
            # ``drugs[]``, so it took no warning at all and re-rendered its
            # stored ``caveats`` verbatim — including, on four of the six
            # documents in ``app/reports/`` when this was written, a sentence
            # naming the cross-compound ranking as the defensible output. That
            # text goes into ``report.html``, which is the file a customer
            # keeps and forwards, so it is the copy of this defect that has
            # already left the building.
            #
            # This WIDENS the condition rather than replacing it: a panel
            # document still takes ``legacy`` (it is a statement about the
            # document's SHAPE that this stamp cannot make), and now every
            # document also gets the era/compliance verdict that
            # ``app/server/jobs.py`` computes for runs — same function, same
            # vocabulary, so the two fronts cannot drift into different words.
            # Nothing is rewritten: the stored ``caveats`` are restored whole
            # and the warning sits beside them.
            cs = claim_set_status(doc.get("claim_set"), doc.get("caveats"))
            if cs["warning"]:
                rep.warnings.append(cs["warning"])
            arms = doc.get("arms") or {}
            if drug and isinstance(arms.get(drug), dict):
                arms = arms[drug]          # flatten the old drug -> arm nesting
            rep.arms = {k: v for k, v in arms.items() if isinstance(v, str)}
            rep.body = {k: doc[k] for k in
                        ("compound", "attribution", "mechanism", "reach",
                         "reach_reason", "microenvironment", "caveats",
                         "claim_set", "window_s", "map_claim") if k in doc}
            with self._lock:
                self._reports[rep.report_id] = rep
                self._order.append(rep.report_id)
                try:
                    self._counter = max(self._counter,
                                        int(rep.report_id.split("_")[-1]))
                except ValueError:
                    pass

    def _persist(self, rep: Report) -> None:
        d = self._dir(rep.report_id)
        private_makedirs(d)
        doc = rep.full_view(self.progress(rep))
        tmp = os.path.join(d, "report.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(jsonable(doc), fh, indent=1)
        own_only(tmp, PRIVATE_FILE_MODE)
        os.replace(tmp, os.path.join(d, "report.json"))
        if rep.status == "done":
            tmp = os.path.join(d, "report.html.tmp")
            with open(tmp, "w") as fh:
                fh.write(render_html(doc))
            own_only(tmp, PRIVATE_FILE_MODE)
            os.replace(tmp, os.path.join(d, "report.html"))

    def html_path(self, report_id: str) -> str:
        return os.path.join(self._dir(report_id), "report.html")

    # -- submission --------------------------------------------------------- #
    def submit(self, section: str, drug: str, spec_overrides: Dict,
               base_spec: Dict, label: Optional[str] = None,
               warnings: Optional[List[str]] = None) -> Report:
        """Enqueue the compound's TWO arms and return the (queued) report."""
        with self._lock:
            self._counter += 1
            report_id = f"rep_{self._counter:05d}"
            rep = Report(report_id, section, drug, base_spec, spec_overrides,
                         label)
            rep.warnings = list(warnings or [])

        # BOTH arms are enqueued BEFORE the report becomes visible. A report that
        # is half-registered would look, to :meth:`refresh`, like a report whose
        # missing arm has nothing pending — i.e. finished — and would be
        # assembled without its counterfactual. Without the counterfactual there
        # is no rho, and the RESISTANCE decade would be missing rather than
        # measured.
        real = dict(base_spec, drug=drug, resistance_off=False)
        cf = dict(base_spec, drug=drug, resistance_off=True)
        rep.arms = {
            "real": self.jobs.submit(real).run_id,
            "counterfactual": self.jobs.submit(cf).run_id,
        }
        rep.status = "running"
        with self._lock:
            self._reports[report_id] = rep
            self._order.append(report_id)
        self._persist(rep)
        return rep

    # -- access ------------------------------------------------------------- #
    def get(self, report_id: str) -> Optional[Report]:
        with self._lock:
            return self._reports.get(report_id)

    def list(self, limit: Optional[int] = None) -> List[Report]:
        with self._lock:
            ids = list(reversed(self._order))
            reps = [self._reports[i] for i in ids if i in self._reports]
        return reps[:limit] if limit else reps

    def cancel(self, report_id: str) -> Optional[Report]:
        rep = self.get(report_id)
        if rep is None or rep.status in TERMINAL:
            return rep
        for rid in rep.run_ids:
            self.jobs.cancel(rid)
        rep.status = "cancelled"
        rep.error = "cancelled by request; every queued arm was cancelled"
        rep.finished_at = time.time()
        self._persist(rep)
        return rep

    # -- progress + assembly ------------------------------------------------- #
    def progress(self, rep: Report) -> Dict:
        states: Dict[str, int] = {}
        stage = None
        for rid in rep.run_ids:
            job = self.jobs.get(rid)
            st = job.status if job else "missing"
            states[st] = states.get(st, 0) + 1
            if job and job.status == "running":
                stage = f"{job.stage} ({rid})"
        total = max(len(rep.run_ids), 1)
        done = states.get("done", 0)
        return {
            "runs_total": len(rep.run_ids),
            "runs_done": done,
            "by_status": states,
            "fraction": round(done / total, 4),
            "stage": stage,
        }

    def refresh(self, rep: Report) -> Report:
        """Recompute the report's status from its arms; assemble when ready.

        Serialised per report: two concurrent polls must not both assemble and
        both persist. Cheap — assembly is pure arithmetic over results the
        worker has already produced.
        """
        if rep.status in TERMINAL:
            return rep
        with self._lock:
            return self._refresh_locked(rep)

    def _refresh_locked(self, rep: Report) -> Report:
        if rep.status in TERMINAL:      # another poll got there first
            return rep
        # Defensive: never assemble before BOTH arms are registered.
        drug = rep.drug
        if not drug or set(rep.arms or {}) != {"real", "counterfactual"}:
            return rep

        arms: Dict[str, Dict] = {}
        pending = 0
        for arm in ("real", "counterfactual"):
            rid = rep.arms[arm]
            job = self.jobs.get(rid)
            if job is None:
                rep.status = "failed"
                rep.error = (
                    f"arm {arm!r} for {drug!r} (run {rid}) is missing from "
                    "the job queue; the report cannot be assembled"
                )
                rep.finished_at = time.time()
                self._persist(rep)
                return rep
            if job.status == "failed":
                rep.status = "failed"
                rep.error = (
                    f"arm {arm!r} for {drug!r} (run {rid}) FAILED: "
                    f"{job.error}"
                )
                rep.finished_at = time.time()
                self._persist(rep)
                return rep
            if job.status == "cancelled":
                rep.status = "cancelled"
                rep.error = (
                    f"arm {arm!r} for {drug!r} (run {rid}) was cancelled. "
                    "Both arms are required and no report was assembled: "
                    "without the resistance-off counterfactual there is no "
                    "denominator for rho, so the RESISTANCE decade would be "
                    "asserted rather than measured, and this document does not "
                    "assert it."
                )
                rep.finished_at = time.time()
                self._persist(rep)
                return rep
            if job.status != "done":
                pending += 1
                continue
            arms[arm] = job.result_view()
        if pending:
            rep.status = "running"
            return rep

        try:
            rep.body = assemble(rep.section, rep.spec_base, drug, arms)
            rep.status = "done"
        except Exception as exc:  # noqa: BLE001 - surfaced, never swallowed
            rep.status = "failed"
            rep.error = redact_exception(
                exc, "this report's arms finished but could not be assembled")
        rep.finished_at = time.time()
        self._persist(rep)
        return rep

    def view(self, rep: Report, full: bool = True) -> Dict:
        self.refresh(rep)
        prog = self.progress(rep)
        return rep.full_view(prog) if full else rep.summary_view(prog)


# --------------------------------------------------------------------------- #
# HTML rendering — self-contained, product-styled, zero external requests
# --------------------------------------------------------------------------- #
def _tokens_css() -> str:
    """The Studio's OWN design tokens, inlined so the document is self-contained.

    Read from ``app/web/css/tokens.css`` rather than duplicated, so a token
    change in the product is a token change in every report rendered after it.
    """
    path = os.path.join(repo_root(), "app", "web", "css", "tokens.css")
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        # A report must render even from a checkout with no front end. These are
        # the handful of tokens the document actually consumes.
        return (
            ":root{--font-ui:'Helvetica Neue',Arimo,Arial,sans-serif;"
            "--font-mono:'SF Mono',Menlo,Consolas,monospace;--radius:2px;"
            "--bg-0:#0d1013;--bg-1:#171b1f;--bg-2:#1e2328;--bg-3:#272d33;"
            "--line-0:#232a30;--line-1:#2f363d;--line-2:#454e57;"
            "--fg-0:#e9edf1;--fg-1:#b3bac1;--fg-2:#828a92;--fg-3:#5c646c;"
            "--accent:#38bcd6;--accent-bg:rgba(56,188,214,.13);--accent-hi:#6fd6ea;"
            "--ok:#3fbf85;--warn:#efab30;--warn-bg:rgba(239,171,48,.14);"
            "--crit:#e8544c;--crit-bg:rgba(232,84,76,.15);--info:#7ba7d8;"
            "--info-bg:rgba(123,167,216,.13);--fs-micro:9.5px;--fs-xs:10.5px;"
            "--fs-sm:11.5px;--fs-md:12.5px;--fs-lg:15px;--fs-xl:19px;}"
        )


#: The light palette the printed document uses, as a last resort. The real
#: source is the product's own ``tokens.css`` — see :func:`_print_tokens`.
_FALLBACK_LIGHT = """
  color-scheme: light;
  --bg-0: hsl(0 0% 100%); --bg-1: hsl(210 24% 97%); --bg-2: hsl(0 0% 100%);
  --bg-3: hsl(210 22% 94%); --bg-4: hsl(210 20% 88%);
  --line-0: hsl(210 18% 91%); --line-1: hsl(210 16% 83%);
  --line-2: hsl(210 14% 68%);
  --fg-0: hsl(213 28% 13%); --fg-1: hsl(213 13% 33%);
  --fg-2: hsl(213 10% 47%); --fg-3: hsl(213 9% 63%);
  --accent: hsl(197 88% 30%); --accent-bg: hsl(197 78% 40% / 0.12);
  --ok: hsl(150 66% 26%); --warn: hsl(30 92% 32%);
  --warn-bg: hsl(34 92% 46% / 0.16); --crit: hsl(356 70% 43%);
  --crit-bg: hsl(356 70% 48% / 0.12); --info: hsl(214 52% 38%);
  --info-bg: hsl(214 52% 44% / 0.12);
"""


def _print_tokens() -> str:
    """The product's LIGHT token values, lifted for the printed page.

    A researcher prints this. The screen document is the Studio's instrument
    surface and is dark on purpose — but a dark ground on paper is either a
    solid black page of toner or, with backgrounds off (the browser default),
    near-white text on white and therefore nothing at all. So print re-declares
    the SAME token names with light values and every component follows, because
    every component is styled through the tokens rather than around them.

    The values are not retyped here: they are cut out of the ``tokens.css`` that
    :func:`_tokens_css` has already inlined, which is the same reason that
    function exists. A token the product re-picks is re-picked on paper too, and
    the two cannot drift.
    """
    css = _tokens_css()
    key = ':root[data-theme="light"]'
    i = css.find(key)
    if i < 0:
        return _FALLBACK_LIGHT
    open_at = css.find("{", i)
    if open_at < 0:
        return _FALLBACK_LIGHT
    depth, j = 0, open_at
    while j < len(css):
        if css[j] == "{":
            depth += 1
        elif css[j] == "}":
            depth -= 1
            if depth == 0:
                return css[open_at + 1:j]
        j += 1
    return _FALLBACK_LIGHT


_DOC_CSS = """
/* ---- report document ---------------------------------------------------
   Same rules as the Studio: square corners (<=2px), mono tabular numerics,
   one accent, semantic colour reserved for state. Nothing in this stylesheet
   references anything off this page - no stylesheet import, no resource
   reference, no webfont, no image, no script. Every graphic is a styled
   element, and the whole document is one file.

   The document's one structural device is a MARGIN RAIL. Section numbers,
   limit severities and field labels sit in it; prose sits beside it at a
   readable measure; tables and code break out across both columns. That single
   grid does four jobs: severity becomes scannable vertically without reading,
   prose stops running 180 characters wide, the section numbering the hero
   cross-references becomes visible, and the printed page gets a real gutter. */
:root {
  --rail: 132px;          /* the margin rail */
  --gut: 20px;            /* rail to content */
  --measure: 76ch;        /* running prose never exceeds this */
  --fs-2xl: 25px;         /* section 2's pick name */
  --fs-3xl: 33px;         /* document title */
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0; background: var(--bg-0); color: var(--fg-1);
  font-family: var(--font-ui); font-size: var(--fs-md); line-height: 1.55;
}
.doc { max-width: 1240px; margin: 0 auto; padding: 0 24px 64px; }
.mono, .num, code, pre, td.num, th.num {
  font-family: var(--font-mono); font-variant-numeric: tabular-nums;
}
h1, h2, h3 { color: var(--fg-0); font-weight: 600; margin: 0;
  text-wrap: balance; }
h1 { font-size: var(--fs-3xl); line-height: 1.15; letter-spacing: -0.012em; }
h3 { font-size: var(--fs-md); margin: 18px 0 6px; }
p { margin: 0 0 9px; max-width: var(--measure); }
a { color: var(--accent); }
a:focus-visible, [tabindex]:focus-visible {
  outline: var(--focus-w, 2px) solid var(--focus, currentColor);
  outline-offset: var(--focus-off, 1px);
}
b, strong { color: var(--fg-0); font-weight: 600; }
.lead { font-size: var(--fs-lg); line-height: 1.5; color: var(--fg-1);
  max-width: var(--measure); }

/* ---- the rail ---------------------------------------------------------- */
.sec { display: grid; grid-template-columns: var(--rail) minmax(0, 1fr);
  column-gap: var(--gut); margin: 40px 0 0; }
.sec > * { grid-column: 2; }
.sec > .full, .sec > h2, .sec > .ledger, .sec > .census { grid-column: 1 / -1; }
h2 {
  display: grid; grid-template-columns: var(--rail) minmax(0, 1fr);
  column-gap: var(--gut); align-items: baseline;
  border-bottom: 1px solid var(--line-2); padding-bottom: 6px;
  margin: 0 0 16px;
}
h2 .sn { font-family: var(--font-mono); font-size: var(--fs-xs);
  color: var(--fg-3); letter-spacing: 0.06em; }
h2 .st { font-size: var(--fs-lg); color: var(--fg-0); font-weight: 600;
  letter-spacing: -0.004em; }
/* An eyebrow in the rail, aligned with the section number above it. */
.railed { display: grid; grid-template-columns: var(--rail) minmax(0, 1fr);
  column-gap: var(--gut); }
.eyebrow { font-family: var(--font-mono); font-size: var(--fs-micro);
  text-transform: uppercase; letter-spacing: 0.09em; color: var(--fg-3); }

/* ---- the banner. Top and bottom, unmissable, never a footnote. --------- */
.nfcu {
  border: 2px solid var(--crit); background: var(--crit-bg);
  padding: 13px 16px; margin: 20px 0; border-radius: var(--radius);
  display: grid; grid-template-columns: var(--rail) minmax(0, 1fr);
  column-gap: var(--gut); align-items: start;
}
.nfcu .t {
  font-family: var(--font-mono); font-size: var(--fs-md); font-weight: 700;
  line-height: 1.25; letter-spacing: 0.08em; color: var(--crit); display: block;
}
/* Wider than the running measure on purpose: this is a stamp, not prose, and
   it should read as two or three emphatic lines rather than a narrow column. */
.nfcu .b { color: var(--fg-0); max-width: 96ch; }

/* ---- masthead ---------------------------------------------------------- */
header.hd { padding: 34px 0 0; }
.hd .kicker { font-family: var(--font-mono); font-size: var(--fs-xs);
  text-transform: uppercase; letter-spacing: 0.16em; color: var(--accent);
  margin-bottom: 10px; }
.hd .sub { color: var(--fg-2); font-size: var(--fs-md); margin-top: 8px;
  max-width: var(--measure); }
/* Three named groups, not one flat dump: what the SPECIMEN is, how the RUN was
   configured, and what identifies the DOCUMENT. A reader auditing a number
   needs to know which of those three a field belongs to. */
.metagrid { display: grid; gap: 0; margin: 24px 0 0;
  border-top: 1px solid var(--line-2); }
.metagrid .grp { display: grid; grid-template-columns: var(--rail) minmax(0, 1fr);
  column-gap: var(--gut); padding: 11px 0; border-bottom: 1px solid var(--line-0); }
.metagrid .gl { font-family: var(--font-mono); font-size: var(--fs-micro);
  text-transform: uppercase; letter-spacing: 0.09em; color: var(--fg-3);
  padding-top: 2px; }
.metagrid .gv { display: flex; flex-wrap: wrap; gap: 4px 28px; }
.meta .k { color: var(--fg-3); font-size: var(--fs-micro); text-transform: uppercase;
  letter-spacing: 0.06em; margin-right: 6px; }
.meta .v { color: var(--fg-0); font-family: var(--font-mono);
  font-variant-numeric: tabular-nums; font-size: var(--fs-sm); }

/* ---- limits ------------------------------------------------------------
   Severity is stated as a WORD in the rail as well as a colour, so it survives
   greyscale printing and colour-vision deficiency. The wall of near-identical
   red the old rendering produced is broken by that word plus the census above:
   a reader can count the criticals before reading any of them. */
/* The ledger's own columns are sized so its titles land on exactly the same
   axis as the section titles and the limit titles below it. One rail, one
   alignment, all the way down the page. */
.ledger { width: 100%; border-collapse: collapse; margin: 0 0 22px;
  table-layout: fixed; }
.ledger td { padding: 5px 10px 5px 0; border-bottom: 1px solid var(--line-0);
  vertical-align: baseline; }
.ledger .lsev { font-family: var(--font-mono); font-size: var(--fs-micro);
  letter-spacing: 0.09em; white-space: nowrap;
  width: calc(var(--rail) + var(--gut) - 30px); }
.ledger .lnum { font-family: var(--font-mono); font-size: var(--fs-micro);
  color: var(--fg-3); width: 30px; text-align: right; }
.ledger a { color: var(--fg-1); text-decoration: none; }
.ledger a:hover { color: var(--fg-0); text-decoration: underline; }
.sev-critical { color: var(--crit); }
.sev-warn { color: var(--warn); }
.sev-info { color: var(--info); }
.census { font-family: var(--font-mono); font-size: var(--fs-sm);
  color: var(--fg-2); margin: 0 0 12px; }
.census b { font-variant-numeric: tabular-nums; }

/* Spans the section grid so the severity word in its rail sits on the SAME
   axis as the section number above it, and stops at the measure so the tinted
   panel is the width of its own text rather than of the viewport. */
.limit {
  --lpad: 21px;                 /* 3px severity rule + 18px padding */
  grid-column: 1 / -1;
  display: grid;
  grid-template-columns: calc(var(--rail) - var(--lpad)) minmax(0, 1fr);
  column-gap: var(--gut);
  border-left: 3px solid var(--line-2); background: var(--bg-1);
  padding: 13px 18px 15px; margin: 0 0 4px;
  max-width: calc(var(--rail) + var(--gut) + var(--measure) + 18px);
}
.limit.critical { border-left-color: var(--crit); background: var(--crit-bg); }
.limit.warn { border-left-color: var(--warn); background: var(--warn-bg); }
.limit.info { border-left-color: var(--info); background: var(--info-bg); }
.limit .lm { font-family: var(--font-mono); font-size: var(--fs-micro);
  letter-spacing: 0.09em; line-height: 1.6; }
.limit .lm .lid { display: block; color: var(--fg-3); letter-spacing: 0.02em;
  margin-top: 3px; font-size: var(--fs-sm); }
.limit .lt { color: var(--fg-0); font-weight: 600; display: block;
  font-size: var(--fs-lg); line-height: 1.35; max-width: var(--measure);
  text-wrap: balance; }
.limit .lb { margin-top: 6px; max-width: var(--measure); }
.limit .ls {
  display: block; margin-top: 8px; font-family: var(--font-mono);
  font-size: var(--fs-micro); color: var(--fg-3); max-width: var(--measure);
  word-break: break-word;
}
.limit.na { opacity: 0.72; }
.limit.na .lt::after {
  content: " — DOES NOT APPLY TO THIS SECTION"; color: var(--fg-3);
  font-weight: 400; letter-spacing: 0.04em; font-size: var(--fs-micro);
}

/* ---- chips ------------------------------------------------------------- */
.chip {
  display: inline-block; font-family: var(--font-mono);
  font-size: var(--fs-micro); text-transform: uppercase; letter-spacing: 0.06em;
  padding: 1px 5px; border: 1px solid var(--line-2); border-radius: var(--radius);
  color: var(--fg-2); white-space: nowrap; line-height: 1.5;
}
.chip.measured { color: var(--ok); border-color: var(--ok); }
.chip.assumed { color: var(--warn); border-color: var(--warn); }
.chip.computed { color: var(--info); border-color: var(--info); }
.chip.crit { color: var(--crit); border-color: var(--crit); }
.chip.accent { color: var(--accent); border-color: var(--accent); }
/* The key for those chips, stated once where they first appear. A colour a
   reader has to infer is a colour that is not carrying information. */
.key { display: flex; flex-wrap: wrap; gap: 6px 22px; margin: 0 0 12px;
  font-size: var(--fs-sm); color: var(--fg-2); align-items: baseline; }
.key > span { display: flex; gap: 7px; align-items: baseline; }

/* ---- tables ------------------------------------------------------------ */
.tblwrap { overflow-x: auto; max-width: 100%; margin: 10px 0 14px; }
table.tbl { border-collapse: collapse; width: 100%; font-size: var(--fs-sm); }
table.tbl th {
  text-align: left; font-weight: 600; color: var(--fg-2);
  font-size: var(--fs-micro); text-transform: uppercase; letter-spacing: 0.06em;
  border-bottom: 1px solid var(--line-2); padding: 6px 9px; white-space: nowrap;
  vertical-align: bottom;
}
/* The column-group row. It is the answer to "which of these did the model
   COMPUTE and which was it GIVEN?", encoded in the table's own structure
   instead of left to the reader to work out drug by drug. */
table.tbl tr.grp th { border-bottom: 0; padding-bottom: 2px; color: var(--fg-3);
  letter-spacing: 0.09em; font-weight: 400; }
table.tbl tr.grp th.g { border-bottom: 1px solid var(--line-1); color: var(--fg-2); }
table.tbl td {
  padding: 6px 9px; border-bottom: 1px solid var(--line-0);
  vertical-align: baseline;
}
table.tbl td.num, table.tbl th.num { text-align: right; }
table.tbl tbody tr:nth-child(even) td { background: var(--bg-1); }
/* `tr.top` survives the ranking: it no longer marks rank 1 (there are no ranks)
   but the row of the reach ladder the headline is READ OFF, so a reader can see
   that the headline is one point on a curve rather than a constant. The two
   rules that styled a rank NUMBER (`td.rank`, `.rank`) are genuinely orphaned
   and are removed; nothing else lost styling with them. */
table.tbl tbody tr.top td { background: var(--accent-bg); }
table.tbl td .u { color: var(--fg-3); font-size: var(--fs-micro); margin-left: 3px; }
.drug { color: var(--fg-0); font-weight: 600; white-space: nowrap; }

/* ---- bars (inline, no images) ------------------------------------------
   `display: inline-block` is load-bearing, not cosmetic. As a plain inline box
   this element's height, min-width and background were all ignored, its
   absolutely-positioned children inherited a zero-size containing block, and
   every bar in the attribution table rendered 0x0 -- the document's only
   quantitative graphic, invisible in every report ever exported. */
.bar { position: relative; display: inline-block; vertical-align: middle;
  height: 9px; width: 128px; background: var(--bg-3); }
.bar > i { position: absolute; top: 0; bottom: 0; display: block; }
.bar > i.pos { background: var(--accent); }
/* fg-2, not fg-3: a negative attribution term is DATA, not a disabled state,
   and fg-3 is the placeholder neutral — it all but vanished on paper. */
.bar > i.neg { background: var(--fg-2); }
.bar > u { position: absolute; top: -2px; bottom: -2px; left: 50%; width: 1px;
  background: var(--line-2); }
/* The unsigned magnitude bar: zero-anchored, one direction. It carried the
   ranking's kill-as-a-fraction-of-top bar; it now carries the reach ladder's
   F(x) and the exposure falloff, both of which are fractions in [0, 1] and are
   the reason it is zero-anchored rather than centred. */
.mag { position: relative; display: inline-block; vertical-align: middle;
  height: 9px; width: 108px; background: var(--bg-3); }
.mag > i { position: absolute; top: 0; bottom: 0; left: 0; display: block;
  background: var(--accent-dim, var(--fg-3)); }
tr.top .mag > i { background: var(--accent); }

/* ---- misc -------------------------------------------------------------- */
pre.code {
  background: var(--bg-1); border: 1px solid var(--line-1); padding: 11px 13px;
  overflow-x: auto; font-size: var(--fs-xs); line-height: 1.5;
  color: var(--fg-1); margin: 10px 0 14px;
}
.kv { display: grid; grid-template-columns: 148px minmax(0, 1fr); gap: 0 16px;
  font-size: var(--fs-xs); margin: 8px 0 0; }
.kv dt { color: var(--fg-3); font-size: var(--fs-micro); text-transform: uppercase;
  letter-spacing: 0.06em; padding: 4px 0; border-bottom: 1px solid var(--line-0); }
.kv dd { margin: 0; color: var(--fg-0); font-family: var(--font-mono);
  line-height: 1.5; padding: 4px 0; border-bottom: 1px solid var(--line-0);
  word-break: break-word; }
.card { border: 1px solid var(--line-1); background: var(--bg-1);
  padding: 14px 16px; margin: 0 0 10px; }
.card .ch { display: flex; flex-wrap: wrap; align-items: baseline; gap: 8px;
  border-bottom: 1px solid var(--line-1); padding-bottom: 8px; margin-bottom: 10px; }
.card .ch .n { color: var(--fg-0); font-weight: 600; font-size: var(--fs-lg); }
ul.pts { margin: 0 0 9px; padding-left: 18px; max-width: var(--measure); }
ul.pts li { margin: 0 0 6px; }
ul.pts.wide, ul.pts.wide li { max-width: none; }
.foot { margin-top: 30px; padding-top: 12px; border-top: 1px solid var(--line-1);
  color: var(--fg-3); font-size: var(--fs-micro); font-family: var(--font-mono);
  line-height: 1.7; }
.hl { color: var(--fg-0); }

/* ---- the hero ----------------------------------------------------------
   Two halves, and the second is not a footnote to the first: the pick, then
   the condition the pick is only true under. They are given the same visual
   weight class deliberately -- a reader who takes the name and stops has taken
   half a sentence. */
.hero { border: 1px solid var(--line-2); background: var(--bg-1);
  padding: 0; margin: 0 0 18px; }
.hero .pick { display: grid; grid-template-columns: var(--rail) minmax(0, 1fr);
  column-gap: var(--gut); padding: 18px 18px 16px;
  border-bottom: 1px solid var(--line-1); }
.hero .pn { font-size: var(--fs-2xl); line-height: 1.15; color: var(--fg-0);
  font-weight: 600; letter-spacing: -0.012em; }
.hero .pn .under { display: block; font-size: var(--fs-sm); font-weight: 400;
  color: var(--fg-2); font-family: var(--font-mono); letter-spacing: 0;
  margin-top: 6px; }
.hero .pchips { display: flex; flex-wrap: wrap; gap: 6px; margin-top: 10px; }
.hero .pfacts { font-family: var(--font-mono); font-size: var(--fs-sm);
  color: var(--fg-1); margin-top: 10px; }
.hero .cond { display: grid; grid-template-columns: var(--rail) minmax(0, 1fr);
  column-gap: var(--gut); padding: 16px 18px; }
.hero .cond + .cond { border-top: 1px solid var(--line-0); }
.hero .cond p:last-child { margin-bottom: 0; }
.hero .cl { font-family: var(--font-mono); font-size: var(--fs-micro);
  text-transform: uppercase; letter-spacing: 0.09em; color: var(--crit);
  line-height: 1.6; }
.hero .cl.q { color: var(--fg-3); }

@media (prefers-reduced-motion: reduce) {
  * { animation: none !important; transition: none !important; }
}
"""


#: Everything in the printed page that is NOT a colour. The colours arrive by
#: re-declaring the tokens, so there is no second palette to keep in step.
_PRINT_RULES = """
  @page { margin: 14mm 12mm; }
  body { font-size: 9.5pt; line-height: 1.45; }
  .doc { max-width: none; padding: 0; }
  :root { --rail: 108px; --gut: 14px; --measure: 92ch; }
  h1 { font-size: 21pt; }
  h2 .st { font-size: 12pt; }
  .hd .kicker { color: var(--fg-2); }
  .sec { margin-top: 22px; }
  h2 { break-after: avoid; page-break-after: avoid; }
  .limit, .card, .hero, .nfcu, tr, .metagrid .grp {
    break-inside: avoid; page-break-inside: avoid;
  }
  thead { display: table-header-group; }
  .tblwrap { overflow-x: visible; }
  table.tbl { font-size: 8pt; }
  table.tbl tbody tr:nth-child(even) td { background: transparent; }
  table.tbl tbody tr.top td { background: var(--accent-bg); }
  /* Backgrounds are the browser's to honour or drop. Severity must survive
     either way, so it is carried by a rule and a word as well as by a wash. */
  .limit { border: 1px solid var(--line-1); border-left-width: 4px; }
  .nfcu { border-width: 2px; }
  .foot { break-inside: avoid; page-break-inside: avoid; }
"""


def _print_css() -> str:
    """The paper stylesheet.

    Re-declares the SAME token names with the product's own light values, so
    every rule in :data:`_DOC_CSS` follows onto paper without one
    print-specific colour being invented. Then it stops the page breaking
    through a limit, a drug card or a table header — a caveat split across a
    page fold is a caveat half its readers will not finish.
    """
    return ("@media print {\n"
            '  :root, :root[data-theme="dark"] {'
            f"{_print_tokens()}\n"
            "    --bg-0: hsl(0 0% 100%);\n"
            "  }\n"
            f"{_PRINT_RULES}"
            "}\n")


def _esc(v) -> str:
    return _html.escape("" if v is None else str(v), quote=True)


def _fmt(v, digits: int = 4) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "—"
    if not math.isfinite(f):
        return "—"
    if f != 0 and (abs(f) < 1e-3 or abs(f) >= 1e5):
        return f"{f:.3g}"
    return f"{f:.{digits}f}"


def _sig(v, digits: int = 3) -> str:
    try:
        f = float(v)
    except (TypeError, ValueError):
        return "—"
    return f"{f:.{digits}g}" if math.isfinite(f) else "—"


def _window_label(window_s: float) -> str:
    """Human exposure window. Never prints "0 h" for a 600 s window."""
    try:
        s = float(window_s)
    except (TypeError, ValueError):
        return "—"
    if not math.isfinite(s) or s <= 0:
        return "—"
    if s >= 3600.0:
        h = s / 3600.0
        return f"{h:.0f} h" if abs(h - round(h)) < 1e-9 else f"{h:.2f} h"
    return f"{s:.0f} s"


def _bar(value: float, span: float = 2.0) -> str:
    """A signed, centred bar. Inline CSS only — no <img>, no external request."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(v):
        return ""
    frac = max(-1.0, min(1.0, v / span)) * 50.0
    if frac >= 0:
        style = f"left:50%;width:{frac:.2f}%"
        cls = "pos"
    else:
        style = f"right:50%;width:{-frac:.2f}%"
        cls = "neg"
    return f'<span class="bar"><u></u><i class="{cls}" style="{style}"></i></span>'


# ``_md_bold`` used to live here. It existed for exactly one caller: the "Rank
# gap" bullet, which rendered ``tmesim.report.attribution.gap_reason``'s ``**x**``
# markers. That bullet was a statement about an adjacent RANK and went with the
# ranking on 2026-08-21, and no other string in this document contains ``**``.
# It is removed as genuinely orphaned code, not as a disclosure: nothing it
# rendered was a caveat, and no text lost its emphasis.


#: The severity words as the document prints them. Same three values
#: ``build_limits`` assigns — spelled out, because a colour alone does not
#: survive a greyscale printer or a red-green reader.
_SEV_WORD = {"critical": "CRITICAL", "warn": "WARNING", "info": "NOTE"}


def _limit_anchor(lim: Dict) -> str:
    """A stable in-page anchor for one limit. Never an external reference."""
    return "lim-" + _esc(lim.get("id") or "unknown")


def _limits_census(limits: Sequence[Dict]) -> str:
    """How conditioned this document is, as a count, before the reading starts.

    Twelve dense paragraphs all washed the same red is not a severity scale, it
    is a wall — and a wall is skimmed. The count is stated first so the weight
    of the block registers even on a reader who then skims it.
    """
    n = {"critical": 0, "warn": 0, "info": 0}
    for lim in limits:
        sev = lim.get("severity") or "info"
        n[sev] = n.get(sev, 0) + 1
    parts = [f'<b class="sev-{k}">{n[k]} {_SEV_WORD[k]}</b>'
             for k in ("critical", "warn", "info") if n.get(k)]
    return (f'<p class="census">{len(limits)} limits travel with this report — '
            + " &middot; ".join(parts) + ". None of them is optional.</p>")


def _limits_ledger(limits: Sequence[Dict]) -> str:
    """An index of the limits, above their full text.

    Two jobs. It lets a reader see the shape of the objection list before
    reading 4,000 words of it, and it gives every limit a name the rest of the
    document can point AT — the hero's condition cites the delivery-law limit,
    and until now that cross-reference resolved to nothing on the page.
    """
    rows = []
    for i, lim in enumerate(limits, start=1):
        sev = lim.get("severity") or "info"
        na = "" if lim.get("applies", True) else \
            ' <span class="eyebrow">does not apply</span>'
        rows.append(
            f'<tr><td class="lnum">{i}</td>'
            f'<td class="lsev sev-{_esc(sev)}">{_SEV_WORD.get(sev, sev.upper())}</td>'
            f'<td><a href="#{_limit_anchor(lim)}">{_esc(lim.get("title"))}</a>'
            f'{na}</td></tr>'
        )
    return (f'<table class="ledger"><tbody>{"".join(rows)}</tbody></table>')


def _limits_html(limits: Sequence[Dict]) -> str:
    out = []
    for i, lim in enumerate(limits, start=1):
        sev = lim.get("severity") or "info"
        cls = sev
        if not lim.get("applies", True):
            cls += " na"
        # The identifier belongs on the full-measure source line, not in the
        # rail: the longest of them is 29 characters and no rail narrow enough
        # to be a rail can set it without breaking it mid-word.
        src_bits = [f'limit id <b class="hl">{_esc(lim.get("id"))}</b>']
        if lim.get("source"):
            src_bits.append(_esc(lim["source"]))
        src = f'<span class="ls">{" &middot; ".join(src_bits)}</span>'
        out.append(
            f'<div class="limit {cls}" id="{_limit_anchor(lim)}">'
            f'<div class="lm sev-{_esc(sev)}">{_SEV_WORD.get(sev, sev.upper())}'
            f'<span class="lid">{i}</span></div>'
            f'<div><span class="lt">{_esc(lim.get("title"))}</span>'
            f'<div class="lb">{_esc(lim.get("body"))}</div>{src}</div></div>'
        )
    return "\n".join(out)


def _banner(where: str) -> str:
    return (
        f'<div class="nfcu" role="alert" data-where="{_esc(where)}">'
        f'<span class="t">NOT FOR CLINICAL USE</span>'
        f'<div class="b">{_esc(NOT_FOR_CLINICAL_USE)}</div></div>'
    )


def _exposure_table(row: Dict) -> str:
    """The exposure map's summary for the ONE compound.

    This replaced the ranked-panel table on 2026-08-21. It is a KEY-VALUE
    reading of one compound rather than a comparison of several, and every line
    of it is a number the run already produced.

    The ``pen. index`` column the ranked table used to carry is still GONE, and
    the reason is a methodology disclosure rather than a ranking artefact, so it
    is restated in the section prose below rather than lost with its host table:
    ``penetration_index = mean(AUC)/p99(AUC)`` is a scale-invariant SHAPE
    statistic — unchanged when the whole exposure field is multiplied by any
    constant — so on real tissue it reports the vascular geometry of the section
    and not the compound (results/delivery_metrics.md §4: it spans 0.3044 to
    0.5107 across a panel whose delivered exposure spans 158x, and temozolomide
    held the HIGHEST index of the GBM panel while killing least). It is retained
    below as a labelled tissue-geometry statistic, not as a delivery number.
    """
    # The column-group row answers, structurally, the question a reader of a
    # simulation output has to be able to answer before trusting any cell: which
    # of these did the model PRODUCE, and which was it HANDED? Left to prose it
    # is answered nowhere; left to the reader it is answered wrongly.
    prov = row.get("provenance") or {}
    ic50 = prov.get("ic50") or {}
    kind = ic50.get("provenance_kind") or "assumed"
    rows_spec = [
        ("computed", "mean tumour-cell exposure",
         f"{_sig(row.get('tumor_auc'))} <span class='u'>uM.s</span>",
         "mean AUC at this tumour's own TUMOR-cell positions"),
        ("computed", "time-averaged C_tumor",
         f"{_sig(row.get('c_tumor_uM'))} <span class='u'>uM</span>",
         "that exposure divided by the window — the concentration the tumour "
         "cells actually saw"),
        ("computed", "median tumour-cell exposure",
         f"{_sig(row.get('median_tumor_auc'))} <span class='u'>uM.s</span>",
         "half the tumour cells sat below this"),
        ("computed", "p10 / p50 / p90 exposure",
         f"{_sig(row.get('E_p10_uM'))} / {_sig(row.get('E_p50_uM'))} / "
         f"{_sig(row.get('E_p90_uM'))} <span class='u'>uM</span>",
         "the SPREAD across the tumour — this is the relative reading, and it "
         "is the defensible one"),
        ("computed", "fitted penetration length",
         f"{_fmt(row.get('penetration_length_fit_um'), 0)} <span class='u'>um</span>",
         "exp+floor fit of AUC against distance from the nearest vessel"),
        ("computed", "tissue-geometry index (retired as a delivery metric)",
         _fmt(row.get("pen"), 4),
         "penetration_index = mean(AUC)/p99(AUC): scale-INVARIANT, so it "
         "reports this section's vascular geometry and not this compound "
         "(results/delivery_metrics.md §4). Retained labelled, never read as "
         "delivery"),
        ("computed", "viable / hypoxic tumour fraction",
         f"{_fmt(row.get('viable_tumor_fraction'), 4)} / "
         f"{_fmt(row.get('hypoxic_tumor_fraction'), 4)}",
         "the engine's own oxygen banding; the reach curve is taken over the "
         "viable set"),
        ("assumed", "plasma mean C_plasma",
         f"{_sig(row.get('c_plasma_uM'))} <span class='u'>uM</span>",
         "the single sustained level the transport solve was driven at — one "
         "multiplicative factor on the source term"),
        ("assumed", "IC50",
         f"{_sig(row.get('ic50'))} <span class='u'>uM</span>",
         _esc(ic50.get("source") or "") + " — " + _esc(ic50.get("note") or "")),
        ("assumed", "D_eff / analytic L",
         f"{_fmt(row.get('d_eff'), 0)} <span class='u'>um2/s</span> / "
         f"{_fmt(row.get('L_um'), 0)} <span class='u'>um</span>",
         "the length scale of the falloff, compiled from descriptors that have "
         "never been scored on 3-D tissue"),
    ]
    head = ("<tr><th>provenance</th><th>quantity</th><th class='num'>value</th>"
            "<th>what it is</th></tr>")
    body = "".join(
        f"<tr><td><span class='chip {k}'>{k}</span></td>"
        f"<td class='drug'>{_esc(label)}</td>"
        f"<td class='num'>{value}</td><td>{what}</td></tr>"
        for k, label, value, what in rows_spec)
    ic_row = (
        f"<tr><td><span class='chip {_esc(kind)}'>{_esc(kind)}</span></td>"
        f"<td class='drug'>IC50 provenance</td>"
        f"<td class='num'>n = {_esc(ic50.get('n_cell_lines'))}</td>"
        f"<td>{_esc(ic50.get('confidence'))} — {_esc(ic50.get('source'))}</td></tr>"
    )
    return (f"<div class='tblwrap'><table class='tbl'><thead>{head}</thead>"
            f"<tbody>{body}{ic_row}</tbody></table></div>")


def _coverage_table(row: Dict) -> str:
    """The IC50 coverage LADDER, never one number.

    All three rungs were computed by every run and only ``1x`` was ever
    rendered. ``1x IC50`` acquired an authority it had not earned by being quoted
    alone; the ladder is printed so the threshold is visibly a choice — the same
    argument ``app/server/reach.py`` makes about its own ``x = 0.1``.
    """
    head = ("<tr><th>threshold</th><th class='num'>fraction of tumour cells "
            "covered</th><th>reading</th></tr>")
    rungs = [
        ("0.5x IC50", row.get("coverage_0_5x"),
         "a half-potency threshold — the most generous of the three"),
        ("1x IC50", row.get("coverage_1x"),
         "the rung this document used to quote alone; it is one point on a "
         "ladder, not a standard"),
        ("2x IC50", row.get("coverage_2x"),
         "a twice-potency threshold — the most demanding of the three"),
    ]
    body = "".join(
        f"<tr><td class='drug'>{_esc(label)}</td>"
        f"<td class='num'><b>{_fmt(v, 4)}</b></td><td>{_esc(note)}</td></tr>"
        for label, v, note in rungs)
    return (f"<div class='tblwrap'><table class='tbl'><thead>{head}</thead>"
            f"<tbody>{body}</tbody></table></div>")


def _reach_block(reach: Optional[Dict], reason: Optional[str]) -> str:
    """The scale-free under-reached reading, CONSUMED from ``app/server/reach.py``.

    Not re-implemented here. If the block is absent or refused, the refusal is
    printed with its reason: a missing section that says nothing is
    indistinguishable from a section whose answer was zero.
    """
    if not isinstance(reach, dict) or reach.get("null_reason") or reason:
        why = (reason or (reach or {}).get("null_reason")
               or "no reason was recorded")
        return (
            '<div class="limit warn"><div class="lm sev-warn">WARNING'
            '<span class="lid">reach</span></div><div>'
            '<span class="lt">The scale-free under-reached curve is NOT '
            'available for this run</span>'
            f'<div class="lb">{_esc(why)}</div>'
            '<span class="ls">app/server/reach.py &middot; no substitute '
            'number was computed, and no absolute exposure figure elsewhere on '
            'this page should be read in its place</span></div></div>'
        )
    curve = reach.get("curve") or {}
    xs = curve.get("x") or []
    fs = curve.get("under_reached_fraction") or []
    head = ("<tr><th class='num'>reach ratio x</th>"
            "<th class='num'>F(x) — fraction of viable tumour cells below it</th>"
            "<th>&nbsp;</th></tr>")
    hx = reach.get("headline_x")
    body = "".join(
        "<tr class='{cls}'><td class='num'>{x}</td><td class='num'><b>{f}</b></td>"
        "<td><span class='mag'><i style='width:{w:.2f}%'></i></span></td></tr>".format(
            cls="top" if (isinstance(hx, (int, float))
                          and abs(_num(x, -1) - float(hx)) < 1e-12) else "",
            x=_sig(x), f=_fmt(f, 4), w=max(0.0, min(1.0, _num(f, 0.0))) * 100.0)
        for x, f in zip(xs, fs))
    falloff = reach.get("falloff") or {}
    d = falloff.get("distance_um") or []
    n = falloff.get("normalised_mean_auc") or []
    c = falloff.get("count") or []
    fhead = ("<tr><th class='num'>distance from nearest vessel</th>"
             "<th class='num'>mean AUC / E_ref</th>"
             "<th class='num'>grid points in bin</th><th>&nbsp;</th></tr>")
    fbody = "".join(
        "<tr><td class='num'>{d}<span class='u'>um</span></td>"
        "<td class='num'><b>{v}</b></td><td class='num'>{c}</td>"
        "<td><span class='mag'><i style='width:{w:.2f}%'></i></span></td></tr>".format(
            d=_fmt(dd, 0), v=_fmt(vv, 4), c=_esc(cc),
            w=max(0.0, min(1.0, _num(vv, 0.0))) * 100.0)
        for dd, vv, cc in zip(d, n, c))
    return (
        f"<div class='tblwrap'><table class='tbl'><thead>{head}</thead>"
        f"<tbody>{body}</tbody></table></div>"
        f"<p><b>Exposure falloff with distance from the nearest vessel</b>, in the "
        f"SAME bins that produced the reference E_ref = "
        f"{_sig(reach.get('E_ref_uM_s'))} "
        # READ the unit, never hard-code it: E_ref is an AUC (uM*s) over a
        # finite window and a CONCENTRATION (uM) on the steady map, and
        # app/server/reach.py says so at its definition. A document that
        # printed the wrong one would be wrong about a physical quantity while
        # looking perfectly well-formed. A record with no E_ref_units predates
        # steady mode entirely, so uM*s is the correct reading for it -- that
        # is an inference about a era that could not have produced the other
        # value, not a default.
        f"{reach.get('E_ref_units') or 'uM*s'}, each divided by it. This is the "
        f"relative reading, and it is the one this document offers as "
        f"defensible.</p>"
        f"<div class='tblwrap'><table class='tbl'><thead>{fhead}</thead>"
        f"<tbody>{fbody}</tbody></table></div>"
        f"<p><b>The threshold is a convention.</b> {_esc(reach.get('convention'))}</p>"
        f"<p><b>What the ratio is and is not protected against.</b> "
        f"{_esc(reach.get('scale_invariance'))}</p>"
        f"<p><b>Which cells this is over.</b> {_esc(reach.get('viable_definition'))} "
        f"E_ref itself is {_esc(reach.get('E_ref_source'))}.</p>"
    )


def _provenance_key() -> str:
    """What the chips mean, stated where they first appear.

    The document has always coloured its provenance chips green / amber / blue
    and never once said what the three colours were. A colour a reader has to
    infer is a colour that is carrying no information.
    """
    return (
        '<div class="key">'
        '<span><span class="chip measured">measured</span>'
        'a number somebody observed in an experiment</span>'
        '<span><span class="chip computed">computed</span>'
        'a number this simulation produced</span>'
        '<span><span class="chip assumed">assumed</span>'
        'a number chosen, curated or defaulted — not measured on this tissue</span>'
        "</div>"
    )


def _decades(v) -> str:
    """A signed decade count, or an em dash where the term was REFUSED.

    Never ``+0.00`` for an absent term: a refusal that renders as a zero is a
    refusal a reader will read as a measurement.
    """
    f = _num(v)
    return f"{f:+.2f}" if math.isfinite(f) else "—"


def _attribution_table(row: Dict, attribution: Dict) -> str:
    """The three absolute decades for the ONE compound.

    Each row states what the term IS, so the number and its meaning arrive
    together. Nothing here is centred on a panel mean — see
    :func:`attribute_compound`.
    """
    head = ("<tr><th>term</th><th class='num'>decades</th>"
            "<th>magnitude</th><th>what it measures</th></tr>")
    meanings = attribution.get("term_meanings") or TERM_MEANING
    labels = attribution.get("term_labels") or TERM_LABEL
    body = []
    for key in (attribution.get("terms") or TERMS):
        v = row.get(key)
        body.append(
            "<tr><td class='drug'>{label}</td><td class='num'><b>{d}</b></td>"
            "<td>{bar}</td><td>{what}</td></tr>".format(
                label=_esc(labels.get(key, key)), d=_decades(v),
                bar=_bar(v) if math.isfinite(_num(v)) else "",
                what=_esc(meanings.get(key, ""))))
    ident = (
        "<tr><td class='drug'>DELIVERY + DOSE</td>"
        f"<td class='num'><b>{_decades((_num(row.get('delivery_decades'), 0.0)) + (_num(row.get('dose_decades'), 0.0)) if (math.isfinite(_num(row.get('delivery_decades'))) and math.isfinite(_num(row.get('dose_decades')))) else float('nan'))}</b></td>"
        f"<td>{_bar(row.get('log10_drive'))}</td>"
        f"<td>= log10(drive), exactly. drive = C_tumor / IC50 = "
        f"{_sig(row.get('drive'))}, i.e. log10 = {_decades(row.get('log10_drive'))}. "
        "This identity is asserted in tests/test_app_reports.py.</td></tr>"
    )
    rho_row = (
        "<tr><td class='drug'>rho (measured)</td>"
        f"<td class='num'><b>{_fmt(row.get('rho'), 4)}</b></td><td></td>"
        f"<td>kill {_fmt(row.get('kill'))} / kill(resistance OFF) "
        f"{_fmt(row.get('kill_no_resistance'))} — this tumour's oxygen-driven "
        f"resistance state costs this compound "
        f"{_fmt((_num(row.get('resistance_cost'), 0.0)) * 100, 1)}% of the kill "
        "it would otherwise achieve</td></tr>"
    )
    return (f"<div class='tblwrap'><table class='tbl'><thead>{head}</thead>"
            f"<tbody>{''.join(body)}{ident}{rho_row}</tbody></table></div>")


def _compound_card(r: Dict, window_s: float) -> str:
    """ONE card, for the one compound.

    ``#{rank}``, the "percentile of this panel" figures on the delivery and
    potency bullets, and the "Rank gap" bullet are GONE: each of them was a
    statement about other compounds, and there are none. Everything that was
    genuinely per-compound is kept — L_um, D_eff, the IC50 provenance, the rho
    counterfactual, MGMT status and the transport parameters — because none of it
    ever depended on the panel.
    """
    prov = r.get("provenance") or {}
    ic50 = prov.get("ic50") or {}
    mgmt = prov.get("mgmt") or {}
    mgmt_cls = "crit" if str(mgmt.get("state", "")).startswith("unmodelled") \
        else "measured" if mgmt.get("state") == "modelled" else "assumed"
    ub = prov.get("sensitivity_upper_bound")
    window = _window_label(window_s)
    L = _num(r.get("L_um"))
    l_ratio = (f"{L / 45.0:.1f}x the calibrated ~45 um perivascular boundary "
               "layer" if math.isfinite(L) else
               "no analytic penetration length was produced for this compound")
    return f"""
<div class="card">
  <div class="ch">
    <span class="n">{_esc(r.get('name'))}</span>
    <span class="chip {'crit' if ub else 'computed'}">
      {'kill = UPPER BOUND' if ub else 'kill = computed'}</span>
    <span class="chip crit">RELATIVE EXPOSURE ONLY</span>
    <span class="chip {mgmt_cls}">MGMT: {_esc(mgmt.get('state'))}</span>
  </div>
  <ul class="pts">
    <li><b>DELIVERY {_decades(r.get('delivery_decades'))} decades</b> —
      log10(C_tumor / C_plasma), the decades LOST getting from blood into this
      tumour's interstitium. Mean exposure at this tumour's own TUMOR-cell
      positions is <span class="mono">{_sig(r.get('tumor_auc'))} uM.s</span> over
      the {window} window, i.e. a time-averaged
      <span class="mono">C_tumor = {_sig(r.get('c_tumor_uM'))} uM</span> against
      a plasma mean of
      <span class="mono">C_plasma = {_sig(r.get('c_plasma_uM'))} uM</span>.
      This term does NOT move when the dose moves: both numbers are exactly
      linear in the plasma level, so the factor cancels from the ratio.
      D_eff = {_fmt(r.get('d_eff'), 0)} um2/s and BBB penetrance
      {_fmt(r.get('bbb'), 3)} give an analytic penetration length
      L = {_fmt(r.get('L_um'), 0)} um ({_esc(l_ratio)}).</li>
    <li><b>DOSE {_decades(r.get('dose_decades'))} decades</b> —
      log10(C_plasma / IC50), the decades of headroom the plasma level has over
      this compound's own potency. IC50 = {_sig(r.get('ic50'))} uM —
      <span class="chip {_esc(ic50.get('provenance_kind') or 'assumed')}">{_esc(ic50.get('provenance_kind'))}</span>
      {_esc(ic50.get('source'))}. {_esc(ic50.get('note'))}
      <b>This is where the dose uncertainty lands, and the decomposition puts it
      here on purpose.</b> The plasma level is a single sustained
      Stage-0 simplification and the IC50 is a monolayer cell-line number; both
      uncertainties are quarantined in this one term instead of being smeared
      across all three.</li>
    <li><b>RESISTANCE {_decades(r.get('resistance_decades'))} decades</b> —
      log10(rho). COUNTERFACTUAL, not asserted: re-simulating this same compound
      on this same tissue with the hypoxic AND quiescent kill penalties switched
      off (run <span class="mono">{_esc(prov.get('counterfactual_run_id'))}</span>)
      moves the predicted kill from {_fmt(r.get('kill'))} to
      {_fmt(r.get('kill_no_resistance'))}, so rho = {_fmt(r.get('rho'), 4)} and
      this tumour's oxygen-driven resistance state costs this compound
      {_fmt((_num(r.get('resistance_cost'), 0.0)) * 100, 1)}% of the kill it
      would otherwise achieve. That is a comparison of this compound with
      ITSELF, which is why it survived the withdrawal of the panel unchanged.</li>
    <li><b>The identity.</b> DELIVERY + DOSE = log10(drive) exactly, where
      drive = C_tumor / IC50 = {_sig(r.get('drive'))}. Both terms are absolute
      decade counts referenced to quantities from THIS run; neither is centred on
      a panel mean, and neither moves if a different compound is or is not
      simulated.</li>
  </ul>
  <dl class="kv">
    <dt>run</dt><dd>{_esc(prov.get('run_id'))} (real) ·
      {_esc(prov.get('counterfactual_run_id'))} (resistance-off)</dd>
    <dt>alpha</dt><dd>{_fmt(prov.get('alpha'), 1)} — {_esc(prov.get('alpha_source'))}</dd>
    <dt>IC50 provenance</dt><dd>{_esc(ic50.get('source'))} · confidence
      {_esc(ic50.get('confidence'))} · n = {_esc(ic50.get('n_cell_lines'))}</dd>
    <dt>sensitivity</dt><dd>mode={_esc(prov.get('sensitivity_mode'))} ·
      mean={_fmt(prov.get('sensitivity_mean'), 4)} ·
      upper_bound={_esc(bool(ub))} ·
      determinant_gene={_esc(prov.get('sensitivity_determinant_gene') or '— none on this panel')}</dd>
    <dt>MGMT</dt><dd>{_esc(mgmt.get('state'))} — {_esc(mgmt.get('detail'))}</dd>
    <dt>transport</dt><dd>MW {_fmt(r.get('mw'), 0)} Da · D {_fmt(r.get('diffusivity'), 0)}
      um2/s · binding {_sig(r.get('binding'))} 1/s ·
      clearance {_sig(r.get('clearance'))} 1/s · hill {_fmt(r.get('hill'), 2)}</dd>
    <dt>plasma</dt><dd>mean {_sig(r.get('c_plasma_uM'))} uM · peak
      {_sig(r.get('plasma_peak_uM'))} uM — read from this run's own
      extra_metrics, never recomputed here</dd>
  </dl>
</div>"""


def render_html(doc: Dict) -> str:
    """Render a report document as ONE self-contained HTML page.

    Self-contained is a hard requirement, not a nicety: the file must render
    identically off a USB stick, from an email attachment, or inside an air-
    gapped review. There is therefore NO ``<img>``, no ``@import``, no
    ``url(...)``, no webfont, no external script and no fetch of any kind — the
    Studio's own ``tokens.css`` is inlined and every graphic is a styled element.
    """
    row = doc.get("compound") or None
    limits = doc.get("limits") or []
    attribution = doc.get("attribution") or {}
    mech = doc.get("mechanism") or {}
    reach = doc.get("reach")
    reach_reason = doc.get("reach_reason")
    micro = doc.get("microenvironment") or {}
    spec = doc.get("spec_base") or {}
    section = doc.get("section_detail") or {}
    window_s = float(doc.get("window_s") or 0.0)
    drug_name = (row or {}).get("name") or doc.get("drug")

    # Three named groups, not sixteen undifferentiated fields. Auditing a number
    # on this page means knowing whether the field beside it describes the
    # TISSUE (fixed, measured, not ours), the RUN (chosen, and every choice is a
    # limit below), or the DOCUMENT (identity and provenance).
    meta_groups = [
        ("specimen", [
            ("section", doc.get("section")),
            ("panel", section.get("panel")),
            ("assay", section.get("assay")),
            ("accession", section.get("accession")),
            ("cells", micro.get("n_cells")),
            ("tumour cells", micro.get("n_tumor")),
            ("outcome label",
             "ABSENT" if not section.get("has_outcome_label") else "present"),
        ]),
        ("run settings", [
            ("grid", " x ".join(str(g) for g in (micro.get("grid") or []))),
            ("dx", f"{_fmt(spec.get('dx_um'), 1)} um"),
            ("window", _window_label(window_s)),
            ("sensitivity mode", doc.get("sensitivity_mode")),
            ("vascular mode", spec.get("vascular_mode")),
            ("oxygen mode", spec.get("oxygen_mode")),
            ("state mode", spec.get("state_mode")),
        ]),
        ("document", [
            ("report", doc.get("report_id")),
            ("status", doc.get("status")),
            ("generated", time.strftime("%Y-%m-%d %H:%M:%S",
                                        time.localtime(doc.get("finished_at")
                                                       or time.time()))),
        ]),
    ]
    meta = "".join(
        f'<div class="grp meta"><div class="gl">{_esc(label)}</div>'
        '<div class="gv">' + "".join(
            f'<span><span class="k">{_esc(k)}</span>'
            f'<span class="v">{_esc(v if v not in (None, "") else "—")}</span>'
            "</span>" for k, v in items) + "</div></div>"
        for label, items in meta_groups)

    # The hero used to name the top of a RANKING. That ranking was withdrawn on
    # 2026-08-21, so the hero now leads with the PER-CELL mechanism reading — the
    # fraction of tumour cells the compound never reached. That number is
    # per-compound, per-cell, computed on the real tissue, and it reconciles
    # bit-exactly with the engine's own kill fraction; it is the answer this
    # product has, not a supporting figure for one.
    #
    # The delivery-law disclosure is NOT deleted with the ranking it used to
    # qualify. It is re-parented: the law sets the source term at the vessel
    # wall, so it conditions the SHAPE of the map, which is a stronger reason to
    # keep it in the hero than the one it had before. The law is the one stored
    # with this document (read from PDEConfig when the runs were made), falling
    # back to the shipped default for a document that predates it. See
    # app/server/delivery_law.py.
    law = doc.get("delivery_law")
    if not isinstance(law, dict):
        law = delivery_law()
    # The law limit is on this same page; the hero's condition links to the
    # limit's own anchor — an in-page fragment, not an external reference.
    law_ref = "#lim-transvascular_delivery_law"
    map_ref = "#lim-map_is_the_output"
    headline = ""
    if row:
        unreached = mech.get("unreached_fraction_of_tumor")
        if unreached is None:
            lead = (
                '<div class="pn">How much of this tumour '
                f'{_esc(drug_name)} never reached: <b>NOT REPORTED</b>'
                '<span class="under">'
                + _esc(mech.get("null_reason") or
                       "the per-cell mechanism attribution is unavailable for "
                       "this run, and no whole-tumour average is substituted "
                       "for it")
                + "</span></div>")
        else:
            lead = (
                '<div class="pn">'
                f'<b>{_fmt(_num(unreached) * 100, 1)}% of this tumour\'s cells '
                f'were never reached</b> by {_esc(drug_name)}'
                '<span class="under">SURVIVED_UNREACHED — per-cell, on this '
                'measured tissue: cells that survived because the drug never '
                'arrived in quantity at their position, not because they '
                'resisted it. This is the answer. '
                f'{_fmt(_num(mech.get("unreached_fraction_of_survivors"), 0.0) * 100, 1)}% '
                'of all SURVIVING tumour cells are in this bucket.</span></div>')
        headline = f"""
<div class="hero">
  <div class="pick">
    <div class="eyebrow">the finding</div>
    <div>
      {lead}
      <div class="pchips">
        <span class="chip crit">RELATIVE EXPOSURE ONLY</span>
        <span class="chip crit">UPPER BOUND</span>
        <span class="chip crit">LAW-DEPENDENT</span>
        <span class="chip crit">NOT A RECOMMENDATION</span>
      </div>
      <div class="pfacts">DELIVERY {_decades(row.get('delivery_decades'))} dec
        &middot; DOSE {_decades(row.get('dose_decades'))} dec
        &middot; RESISTANCE {_decades(row.get('resistance_decades'))} dec
        &middot; L = {_fmt(row.get('L_um'), 0)} um
        &middot; IC50 {_sig(row.get('ic50'))} uM
        ({_esc(((row.get('provenance') or {}).get('ic50') or {}).get('source'))})
        &middot; tumour-cell exposure {_sig(row.get('tumor_auc'))} uM.s</div>
    </div>
  </div>
  <div class="cond">
    <div class="cl">what this document claims</div>
    <div>
      <p><b>{_esc(MAP_CLAIM)}</b> Read WHERE the drug went, relative to
         elsewhere in this same tissue. Do not read how much of the tumour it
         killed, and do not compare any number on this page with a number from a
         different compound (<a href="{map_ref}">see the map limit above</a>).</p>
      <p><b>This is a model output, not a recommendation.</b> It has never been
         checked against a measured response for this tissue or any other.
         Treat it as a hypothesis to test in the lab, nothing more.</p>
    </div>
  </div>
  <div class="cond">
    <div class="cl">what was measured, and how far</div>
    <div>
      <p><b>Two perturbations were swept, at two DIFFERENT strengths, and that
         list is not a coverage claim.</b> The plasma / dose scale is EXACT for
         this finite-window exposure field — worst relative deviation 1.41e-14
         over plasma factors 0.005 to 26, a 5,200x span, measured at the shipping
         configuration (results/pk_ranking_robustness.md sec. 1.2) — so a
         per-compound plasma factor moves every point of this map together and
         cancels from a ratio between points. Tissue capacity is EXACT FOR A
         CONVERGED GRADIENT SHAPE ONLY, verified to 1e-9
         (results/capacity_ranking_robustness.md sec. 2b); a capacity term acts
         as a pure rescale of time, so over a FINITE window a region that has not
         converged carries a residual. <b>Those are the only two axes that were
         swept. An axis absent from that list is UNTESTED, not covered.</b></p>
      <p><b>Neither protects against error in the transport law itself.</b>
         Every exposure number here was produced under
         <span class="mono">{_esc(law.get('law'))}</span>; the alternative source
         term <span class="mono">{_esc(law.get('alternative_law'))}</span>
         re-weights delivery by permeability, which spans <b>47x</b> across this
         engine's GBM compound list — enough that on the panel this document used
         to rank it moved the top-1 pick outright, <b>cisplatin &rarr;
         vincristine</b>, swapping 4 of 45 pairs. That swap is quoted here as
         EVIDENCE OF THE LAW'S LEVERAGE and not as a comparison: it is how large
         a re-shaping of delivery the choice of source term is worth. A wrong
         effective diffusivity does the same thing to the SHAPE of the gradient
         rather than to its level, which is exactly what the two invariances
         above cannot divide out (<a href="{law_ref}">see the delivery-law limit
         above</a>; results/transvascular_exchange_revalidation.md §4.1). That
         law's descriptor inputs have never been scored
         (results/descriptor_falsification.md).</p>
    </div>
  </div>
  <div class="cond">
    <div class="cl">the potency readout that is still missing</div>
    <div>
      <p><b>Until an ex-vivo potency readout exists on this specimen's own
         cells, everything on this page that reads an IC50 is not
         decision-grade.</b> That is the DOSE decade and the whole IC50 coverage
         ladder — not the exposure map, which reads no potency at all. This
         sentence used to appear here as a caveat on a RANKING, in
         results/ranking_robustness.md §8's words: <b>"the reportable output of
         this pipeline for the top of the panel is a SET, not a winner"</b>,
         with P(cisplatin = rank 1) at 0.931 under that report's declared priors
         falling to <b>0.587</b> in its stress scenario — and the measurement it
         named as the fix was a patient-specific potency readout. The ranking is
         gone; the missing measurement is not, and it lands on exactly the terms
         named above. A per-compound IC50 is still a monolayer cell-line number
         carrying none of this tumour's 3-D architecture.</p>
    </div>
  </div>
</div>"""

    caveat_items = "".join(f"<li>{_esc(c)}</li>" for c in (doc.get("caveats") or []))
    warn_items = "".join(f"<li>{_esc(c)}</li>" for c in (doc.get("warnings") or []))

    # Section 05 renders micro["patient_id"] under the label "specimen id". The
    # divergence is deliberate: the field name is load-bearing across the engine
    # and the API payloads and renaming it is a data-model change, but a field
    # labelled "patient id" on a rendered page invites a user to type one, and
    # this tool is scoped to de-identified tissue (docs/DEPLOYMENT.md §0). If the
    # field is ever renamed, the label here is already correct.
    body = f"""
<div class="doc">
<header class="hd">
  <div class="kicker">tmesim · single-compound exposure map</div>
  <h1>Where {_esc(drug_name or 'this compound')} reaches in
    {_esc(section.get('label') or doc.get('section'))} — and how much of the
    tumour it never reaches</h1>
  <div class="sub">ONE compound, ONE measured specimen: the exposure map, the
    fraction of tumour cells the drug never reached, and the provenance and
    limits of every number that says so. Every number below is a prediction; the
    sections are ordered so that the limits on it are read first. This document
    ranked a drug panel until 2026-08-21; that ranking is WITHDRAWN and the
    reason is limit 2.</div>
  <div class="metagrid">{meta}</div>
</header>

{_banner('top')}

<section class="sec">
<h2><span class="sn">01</span><span class="st">Limits — read these before any number on this page</span></h2>
{_limits_census(limits)}
{_limits_ledger(limits)}
{_limits_html(limits)}
</section>

<section class="sec">
<h2><span class="sn">02</span><span class="st">The single-compound reading — where the drug went</span></h2>
<div class="full">{headline}</div>
{_provenance_key()}
<div class="full">{_exposure_table(row or {})}</div>
<p>{_esc(drug_name or 'This compound')} was simulated on ONE reconstruction of
this specimen, over a {_window_label(window_s)} exposure window, at seed
{_esc(((row or {}).get('provenance') or {}).get('seed'))}, through the
standard compilation pipeline (<span class="mono">calibrate_alpha &rarr;
calibrated_drugparams &rarr; apply_barrier_to_drugparams</span>). The
counterfactual arm reused the SAME reconstruction, so the two arms differ in
exactly one setting and rho is a ratio rather than a comparison of two tissues.
Every level in the table above is the centre of a wide interval — see limit
"Where absolute kill IS defensible, it is an INTERVAL".</p>

<h3>How much of the tumour is UNDER-REACHED — the scale-free reading</h3>
<p>This is the reading that survives not knowing the dose. It is a RATIO between
two points of one map — each viable tumour cell's exposure over the perivascular
reference exposure <span class="mono">E_ref</span> — so a factor that multiplies
every point of the map divides out of it, at the two strengths the hero states
and the map limit cites. It is computed by
<span class="mono">app.server.reach.compute_reach</span> and consumed here; this
document does not re-implement it, because a second implementation of a number
that must agree is how two numbers stop agreeing.</p>
<div class="full">{_reach_block(reach, reach_reason)}</div>

<h3>IC50 coverage across the dose ladder</h3>
<p><b>Not one number.</b> <span class="mono">therapeutic_coverage</span> is the
fraction of this section's tumour cells whose time-averaged exposure exceeds a
multiple of this compound's own IC50 over the {_window_label(window_s)} window.
All three rungs were already computed by every run and only
<span class="mono">1x</span> was ever printed, which gave that rung an authority
it had not earned. <b>Coverage is an ABSOLUTE quantity and is therefore not
protected by anything in limit 2:</b> it reads the plasma level and the IC50,
which are exactly the two inputs the robustness studies measured as uncertain
(per-compound plasma exposure spans 113.1x across the compounds the GBM panel
ranked, results/pk_ranking_robustness.md sec. 0). Rescale the dose and every
number in this table moves; the reach curve above does not. Read the ladder, and
prefer the curve.</p>
<div class="full">{_coverage_table(row or {})}</div>
<p><b>On the retired penetration index.</b> The exposure table above still
carries <span class="mono">penetration_index = mean(AUC)/p99(AUC)</span>, but
labelled as a TISSUE-GEOMETRY statistic and never as a delivery number. It is a
scale-invariant SHAPE statistic: multiply the whole exposure field by any
constant and it does not move, so on real tissue it reports this section's
vascular geometry rather than this compound. It was retired as a delivery metric
for that reason (results/delivery_metrics.md §4) and the withdrawn GBM panel
demonstrated why — temozolomide held the HIGHEST penetration index of that panel
while killing least. Neither it nor coverage is exempt from anything else on this
page: both are produced under the delivery law named in the hero, and both carry
the same unquantified mesh term as everything else.</p>
</section>

<section class="sec">
<h2><span class="sn">03</span><span class="st">Why this compound did what it did on this tissue</span></h2>
<p>Three <b>absolute log10 terms</b> ("decades"), each referenced to a quantity
from <b>this same run</b> — never to a panel mean, and never to whichever other
compounds happened to be selected. <span class="mono">DELIVERY + DOSE</span> is
exactly <span class="mono">log10(drive)</span>. Nothing is asserted: the
RESISTANCE term comes from a real counterfactual simulation of this compound
against itself, which is why it was simulated twice.</p>
<div class="full"><pre class="code">{_esc(attribution.get('method'))}</pre>
{_attribution_table(row or {}, attribution)}</div>
{('<div class="limit critical"><div class="lm sev-critical">CRITICAL'
  '<span class="lid">attribution</span></div><div>'
  '<span class="lt">The decomposition REFUSED — no decade is reported</span>'
  f'<div class="lb">{_esc(attribution.get("null_reason"))}</div>'
  '<span class="ls">app.server.report.attribute_compound &middot; nothing was '
  'substituted for the missing quantity</span></div></div>')
 if attribution.get('null_reason') else ''}
{('<div class="limit warn"><div class="lm sev-warn">WARNING'
  '<span class="lid">rho</span></div><div>'
  '<span class="lt">The RESISTANCE decade was not reported</span>'
  f'<div class="lb">{_esc(attribution.get("resistance_null_reason"))}</div>'
  '<span class="ls">app.server.report.attribute_compound &middot; DELIVERY and '
  'DOSE are unaffected: neither reads a kill</span></div></div>')
 if attribution.get('resistance_null_reason') else ''}
<p><b>Why this decomposition, and not the one this report used to print.</b>
{_esc(attribution.get('rationale'))}</p>
</section>

<section class="sec">
<h2><span class="sn">04</span><span class="st">This compound's reading and provenance</span></h2>
<div class="full">{_compound_card(row, window_s) if row else
  '<p>No compound block was assembled for this document, so there is nothing to '
  'report here. The limits above still apply in full.</p>'}</div>
</section>

<section class="sec">
<h2><span class="sn">05</span><span class="st">Specimen</span></h2>
<dl class="kv">
  <dt>specimen id</dt><dd>{_esc(micro.get('patient_id'))}</dd>
  <dt>source</dt><dd>{_esc(micro.get('source'))}</dd>
  <dt>cells</dt><dd>{_esc(micro.get('n_cells'))} measured
    ({_esc(micro.get('n_tumor'))} TUMOR)</dd>
  <dt>extent</dt><dd>{_esc(' x '.join(_fmt(v, 0) for v in (micro.get('extent_um') or [])))} um</dd>
  <dt>mean vascular density</dt><dd>{_fmt(micro.get('mean_vascular_density'), 4)}</dd>
  <dt>mean ECM density</dt><dd>{_fmt(micro.get('mean_ecm_density'), 4)}</dd>
  <dt>mean oxygen</dt><dd>{_fmt(micro.get('mean_oxygen'), 4)}</dd>
  <dt>tissue mask</dt><dd>{_esc(micro.get('tissue_mask_mode'))} ·
    coverage {_fmt((micro.get('tissue_coverage_fraction') or 0) * 100, 1)}%</dd>
  <dt>cell state route</dt><dd>{_esc((micro.get('cell_state') or {}).get('mode'))} ·
    measured={_esc((micro.get('cell_state') or {}).get('is_measured'))}</dd>
  <dt>outcome label</dt><dd>ABSENT — untreated resection, nothing to score against</dd>
</dl>
</section>

<section class="sec">
<h2><span class="sn">06</span><span class="st">Caveats carried by the individual runs</span></h2>
<ul class="pts">{caveat_items or '<li>none reported</li>'}</ul>
{('<h3>Spec warnings at submission</h3><ul class="pts">' + warn_items + '</ul>') if warn_items else ''}
</section>

<section class="sec">
<h2><span class="sn">07</span><span class="st">What would have to be true for this to be trustworthy</span></h2>
<ul class="pts">
  <li><b>A SPATIALLY RESOLVED drug measurement on tissue like this — this is the
    experiment that would falsify the map claim, and it is not the one this
    section used to name.</b> It used to ask for a tumour with a known outcome so
    that "the predicted ranking is scored against the observed one". There is no
    predicted ranking any more, and an outcome label would not test this document
    anyway: the claim is about relative exposure BETWEEN REGIONS of one specimen,
    so the falsifying measurement is a measurement of drug concentration at known
    distances from known vessels in the same tissue — mass-spectrometry imaging,
    a labelled-compound penetration profile, or an equivalent. Score the measured
    falloff against the predicted falloff. If the predicted normalised exposure
    at 100 um from a vessel is materially wrong relative to the perivascular
    value, THE MAP CLAIM IS FALSE, and the two invariances in limit 2 do not
    save it — they protect the ratio against a dose factor and a capacity term,
    not against a wrong gradient.</li>
  <li><b>A scored transport law.</b> The delivery law and the effective
    diffusivity set the SHAPE of the gradient, which is the entire claim. The
    law's descriptor inputs have never been scored
    (results/descriptor_falsification.md), and the compiled D_eff misses its own
    measured multicellular-layer band by 1.5x. Scoring those against 3-D
    penetration data would test the map claim more directly than any outcome
    label could.</li>
  <li><b>A potency measurement in context.</b> Cell-line IC50 must be replaced by
    (or at minimum reconciled with) a 3-D / organoid / patient-derived
    measurement from tissue resembling this one. This falsifies the DOSE decade
    and the coverage ladder — the absolute parts — not the map.</li>
  <li><b>More than one penetration anchor.</b> alpha rides on a single
    doxorubicin measurement; independent per-compound perivascular-distance data
    would test whether the predicted falloff transfers outside the
    strongly-bound class. alpha re-shapes the gradient rather than scaling it, so
    this is a test OF the map claim and not merely of its level.</li>
  <li><b>A measured per-cell sensitivity</b>, replacing today's explicit
    assumption (uniform 1.0), which is what makes every kill on this page an
    upper bound rather than a prediction.</li>
  <li><b>A panel that can measure cell state.</b> Until hypoxia and
    proliferation are measurable on the section itself, the absolute kill has no
    defensible level under any state route — and the RESISTANCE decade and the
    "viable" band the reach curve is taken over both inherit that.</li>
  <li><b>A converged window, or a measurement of how far from converged this one
    is.</b> The tissue-capacity invariance in limit 2 is exact for a CONVERGED
    gradient shape and carries a residual over a finite window. Nothing on this
    page asserts that this window converged anywhere in this tissue, and
    measuring that residual would turn a conditional invariance into an
    unconditional one or refute it.</li>
</ul>
</section>

{_banner('bottom')}

<div class="foot">
  Generated by {_esc(doc.get('generator'))} ·
  attribution: {_esc(doc.get('attribution_implementation'))} ·
  runs: {_esc(', '.join(doc.get('run_ids') or []))} ·
  deterministic (seed
  {_esc(((row or {}).get('provenance') or {}).get('seed'))}).
  Self-contained document: no external stylesheet, font, image, script or
  network request.
</div>
</div>"""

    return (
        "<!doctype html>\n"
        '<html lang="en" data-theme="dark"><head><meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
        '<meta name="color-scheme" content="dark light">\n'
        '<meta name="robots" content="noindex">\n'
        # Matches the <h1>. It read "tmesim therapy report" until 2026-08-15 and
        # "simulated compound ranking" until 2026-08-21; naming a deliverable a
        # "therapy report" is an implication this product does not make, and
        # naming it a ranking is a claim this product has withdrawn. NOT FOR
        # CLINICAL USE stays IN the title on purpose: this document travels off a
        # USB stick, and the browser tab is the one part of it that survives a
        # screenshot of the window frame.
        f"<title>Codon exposure map {_esc(doc.get('report_id'))} — "
        f"{_esc(drug_name or 'no compound')} on {_esc(doc.get('section'))} — "
        f"NOT FOR CLINICAL USE</title>\n"
        f"<style>\n{_tokens_css()}\n{_DOC_CSS}\n{_print_css()}\n</style>\n"
        "</head><body>\n" + body + "\n</body></html>\n"
    )
