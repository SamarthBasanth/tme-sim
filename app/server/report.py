"""Therapy reports — the deliverable the pipeline exists to produce.

A **report** is one ranked, mechanistically attributed drug panel simulated on
ONE section, with the provenance of every input and the limits of the output
attached to it. It is the in-product form of ``experiments/patient_report.py``.

How it is built
---------------
Nothing here simulates anything. A report is assembled out of RUNS:

* every arm is a normal ``RunSpec`` submitted through the ordinary
  :class:`app.server.jobs.JobManager` queue, so it shares the single worker, the
  reconstruction cache, cancellation, on-disk persistence and ``/api/runs``;
* each drug contributes TWO arms on the same reconstruction — the real run, and
  the ``resistance_off=true`` counterfactual whose kill is the denominator of
  ``rho`` in the attribution;
* when every arm is ``done`` the report is assembled by
  :func:`tmesim.report.attribution.attribute` — the SAME function
  ``experiments/patient_report.py`` calls. The science is imported, not forked.

Because all 2N arms share ``(section, dx, vascular_mode, oxygen_mode,
state_mode)``, the expensive reconstruction is built once and every subsequent
arm is a cache hit. That is the whole reason the report goes through the queue
instead of around it.

What it refuses to do
---------------------
* It will not rank in ``sensitivity_mode='determinant'``. On a panel that
  measures some drugs' determinants and not others that mode is not a ranking,
  it is a measurement-coverage artifact (``results/real_drug_panel.md``).
* It will not present an absolute kill without the interval and the state-route
  span that ``results/kill_calibration.md`` and ``results/cell_state_impact.md``
  measured for it.
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
from .jobs import PDE_DT_S, JobManager, jsonable

__all__ = [
    "DEFAULT_PANEL",
    "NOT_FOR_CLINICAL_USE",
    "REPORT_SENSITIVITY_MODE",
    "Report",
    "ReportManager",
    "build_limits",
    "render_html",
    "reports_dir",
]

#: The five-drug GBM panel ``experiments/patient_report.py`` ranks, in its order.
DEFAULT_PANEL: Sequence[str] = (
    "temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin",
)

#: The ONLY per-cell sensitivity mode in which a cross-drug ranking is fair: it
#: is identical for every drug, so it cancels out of the comparison. Forced, not
#: defaulted — a report is a ranking.
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
def build_limits(section_id: str, spec: Dict, drugs: Sequence[str]) -> List[Dict]:
    """The limits that travel with EVERY report, in the order they are shown.

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
    limits.append({
        "id": "ranking_is_the_output",
        "severity": "critical",
        "applies": True,
        "title": "The RANKING is the defensible output. The absolute kill is not.",
        "source": "results/cell_state_impact.md §3.2, §3.3, §5",
        "body": (
            "The cross-drug ordering below is INVARIANT across all three "
            "cell-state routes ('oxygen', 'vessel_distance', 'markers') and "
            "across the whole measured potency band. The absolute kill fraction "
            "is not: with tissue, drug, transport and RNG held identical, kill "
            "moves by up to 0.3252 absolute (1.77x relative) purely with the "
            "choice of state route. Read the order. Do not read the level."
        ),
    })

    # 2b ------------------------------------------------------------------- #
    # Limit 2 says the ordering is invariant. It is — to the CELL-STATE route,
    # which is a common-mode perturbation. It is NOT invariant to the
    # transvascular delivery law, which is drug-specific and changes top-1 on
    # GBM. Those two sentences have to sit next to each other or the first one
    # reads as a broader guarantee than it is. The law is read from PDEConfig at
    # call time (app/server/delivery_law.py), so this limit cannot describe a
    # law the engine is no longer using.
    law = delivery_law()
    limits.append({
        "id": "transvascular_delivery_law",
        "severity": "critical",
        "applies": True,
        "title": law["title"] + " — the ranking above is LAW-DEPENDENT",
        "source": law["reference"],
        "body": law["caveat"],
    })

    # Limit 2c. The two limits above condition the ranking on things that are
    # uncertain. This one names the case where the ranking is measured to be
    # WRONG, and it has to sit with them rather than in the numerics block at
    # the bottom: a reader who is ranking a congeneric series is the reader most
    # likely to stop after limit 2, and that reader is the one this concerns.
    # The engine's only molecular input to diffusivity is molecular weight, so
    # at matched MW its transport prediction is flat BY CONSTRUCTION and the
    # ordering it produces over such a series is potency alone wearing a
    # transport label. Scored blind, unfitted, on 71 compounds.
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
            "CONSEQUENCE FOR THIS PAGE: if the drugs ranked below are a "
            "similar-mass series, read their ORDER as a hypothesis worth "
            "testing rather than as a measurement, and read no absolute "
            "delivery number from it at all. Across drugs that differ "
            "substantially in mass the delivery term discriminates strongly, "
            "and that remains the best-supported case."
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
            "Every drug carries per-cell sensitivity 1.0, so NO intrinsic "
            "resistance is asserted for any of them. That is deliberate and it "
            "is the only setting in which a cross-drug ranking is fair, because "
            "an identical factor cancels out of the comparison. The price is "
            "that every ABSOLUTE kill below is a no-intrinsic-resistance UPPER "
            "BOUND. The alternative 'determinant' mode is NOT used here: it "
            "scales doxorubicin (TOP2A) and cisplatin / lomustine / "
            "methotrexate (MKI67 proxy) to ~0.585 while temozolomide keeps 1.0 "
            "— because TMZ's determinant MGMT is absent from the panel — which "
            "flatters TMZ ~1.7x for a measurement-coverage reason, not a "
            "pharmacological one."
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
            "ranking does not rest on the oxygen field, and no absolute hypoxia "
            "claim should rest on it either. AN EARLIER VERSION OF THIS LIMIT "
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
            "outcome-labelled section anywhere in this repository. The ranking "
            "below has therefore never been scored: its accuracy is not 'good', "
            "not 'bad', it is UNKNOWN. Potency additionally comes from GDSC "
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
            "mitoxantrone. Every drug's absolute penetration rides on that one "
            "anchor, and the underlying binding proxy is a structure-based "
            "ORDERING, not a first-principles rate. A drug it mis-ranks (active "
            "transport, efflux, ion trapping) would be mis-predicted."
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
            f"solve at dt = {PDE_DT_S:g} s; it preserves total delivered drug "
            "and the relative ranking, but it is not a time-resolved infusion "
            "and no schedule optimisation should be read off it."
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


def _row_from_arms(drug: str, real: Dict, cf: Dict, mean_ecm: float,
                   section_id: str) -> Dict:
    """One panel row, built ONLY out of numbers the engine already produced."""
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
        # --- the four inputs tmesim.report.attribution needs ---------------- #
        "tumor_auc": _num(metrics.get("mean_tumor_auc"), 0.0),
        "ic50": _num(prov.get("ic50_uM"), float("nan")),
        "kill": _num(metrics.get("kill_fraction"), 0.0),
        "kill_no_resistance": _num(cf_metrics.get("kill_fraction"), 0.0),
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
        "coverage_1x": _num(metrics.get("therapeutic_coverage_1x")),
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


def assemble(section_id: str, spec_base: Dict, drugs: Sequence[str],
             arms: Dict[str, Dict[str, Dict]]) -> Dict:
    """Build the report body from completed run results.

    ``arms`` maps ``drug -> {'real': RunResult, 'counterfactual': RunResult}``.
    """
    from tmesim.report.attribution import TERM_LABEL, TERMS, attribute, gap_reason

    first = arms[drugs[0]]["real"]
    report = first.get("report") or {}
    mean_ecm = _num(report.get("mean_ecm_density"), 0.0)
    window_s = float(spec_base.get("total_time_s") or 0.0)

    rows = [
        _row_from_arms(d, arms[d]["real"], arms[d].get("counterfactual"),
                       mean_ecm, section_id)
        for d in drugs
    ]
    check = attribute(rows, window_s=window_s)
    rows.sort(key=lambda r: r["kill"], reverse=True)

    for i, r in enumerate(rows, start=1):
        r["rank"] = i
        if i < len(rows):
            r["gap_reason"] = gap_reason(r, rows[i])
        elif len(rows) > 1:
            r["gap_reason"] = gap_reason(rows[i - 2], r, from_below=True)
        else:
            r["gap_reason"] = None

    sp = check["spearman_total_vs_kill"]
    check["faithful"] = (
        "faithful" if (sp is not None and math.isfinite(sp) and sp >= 0.9)
        else "partially faithful" if (sp is not None and math.isfinite(sp) and sp >= 0.5)
        else "NOT faithful"
    )

    # The union of every arm's own caveats, de-duplicated, order preserved.
    seen, caveats = set(), []
    for d in drugs:
        for arm in ("real", "counterfactual"):
            run = arms[d].get(arm) or {}
            for c in ((run.get("provenance") or {}).get("caveats") or []):
                if c not in seen:
                    seen.add(c)
                    caveats.append(c)

    return {
        "ranking": jsonable(rows),
        "attribution": {
            "terms": list(TERMS),
            "term_labels": dict(TERM_LABEL),
            "check": jsonable(check),
            "method": (
                "C_tumor_i = (mean AUC at this tumour's TUMOR-cell positions) / T   [uM]\n"
                "drive_i   = C_tumor_i / IC50_i                                     [-]\n"
                "\n"
                "delivery_i = log10 C_tumor_i - mean_j log10 C_tumor_j\n"
                "potency_i  = -( log10 IC50_i - mean_j log10 IC50_j )\n"
                "o2resist_i = log10 rho_i     - mean_j log10 rho_j\n"
                "\n"
                "  rho_i = kill_i / kill_i(resistance OFF)   <- MEASURED by a second\n"
                "                                              simulation of the same drug\n"
                "                                              on the same tissue with\n"
                "                                              hypoxia_resistance=1 and\n"
                "                                              quiescent_resistance=1\n"
                "\n"
                "identity:  delivery_i + potency_i  ==  log10 drive_i - mean_j log10 drive_j"
            ),
            "implementation": "tmesim.report.attribution.attribute",
        },
        "microenvironment": jsonable(report),
        "caveats": caveats,
        "window_s": window_s,
    }


# --------------------------------------------------------------------------- #
# Report objects
# --------------------------------------------------------------------------- #
class Report:
    """One therapy report: its spec, its arms, its status, its document."""

    def __init__(self, report_id: str, section: str, drugs: Sequence[str],
                 spec_base: Dict, spec_overrides: Dict,
                 label: Optional[str] = None):
        self.report_id = report_id
        self.section = section
        self.drugs = list(drugs)
        self.spec_base = spec_base
        self.spec_overrides = dict(spec_overrides or {})
        self.label = label
        self.created_at = time.time()
        self.finished_at: Optional[float] = None
        self.status = "queued"
        self.error: Optional[str] = None
        self.warnings: List[str] = []
        #: drug -> {'real': run_id, 'counterfactual': run_id}
        self.arms: Dict[str, Dict[str, str]] = {}
        self.body: Optional[Dict] = None

    # -- views ------------------------------------------------------------- #
    @property
    def run_ids(self) -> List[str]:
        out: List[str] = []
        for d in self.drugs:
            for arm in ("real", "counterfactual"):
                rid = self.arms.get(d, {}).get(arm)
                if rid:
                    out.append(rid)
        return out

    def summary_view(self, progress: Optional[Dict] = None) -> Dict:
        best = None
        if self.body and self.body.get("ranking"):
            top = self.body["ranking"][0]
            best = {"drug": top["name"], "kill": top["kill"],
                    "dominant": top.get("dominant")}
        return {
            "report_id": self.report_id,
            "label": self.label,
            "status": self.status,
            "section": self.section,
            "drugs": self.drugs,
            "n_runs": len(self.run_ids),
            "created_at": self.created_at,
            "finished_at": self.finished_at,
            "runtime_s": (round(self.finished_at - self.created_at, 3)
                          if self.finished_at else None),
            "progress": progress,
            "top_ranked": best,
            # ``top_ranked`` is a RANKING, read straight off the list without
            # opening the document, so the law that produced it travels with the
            # SUMMARY and not only with the full view's limits block. Read from
            # PDEConfig at call time (app/server/delivery_law.py).
            "delivery_law": delivery_law(),
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
            "limits": build_limits(self.section, self.spec_base, self.drugs),
            "section_detail": _safe_section(self.section),
            "generator": "app/server/report.py (POST /api/reports)",
            "attribution_implementation": "tmesim.report.attribution.attribute",
        })
        if self.body:
            out.update(self.body)
        else:
            out.setdefault("ranking", [])
            out.setdefault("attribution", None)
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
            rep = Report(doc.get("report_id", name), doc.get("section", ""),
                         doc.get("drugs") or [], doc.get("spec_base") or {},
                         doc.get("spec_overrides") or {}, doc.get("label"))
            rep.status = doc.get("status", "done")
            rep.created_at = doc.get("created_at") or rep.created_at
            rep.finished_at = doc.get("finished_at")
            rep.error = doc.get("error")
            rep.warnings = doc.get("warnings") or []
            rep.arms = doc.get("arms") or {}
            rep.body = {k: doc[k] for k in
                        ("ranking", "attribution", "microenvironment",
                         "caveats", "window_s") if k in doc}
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
    def submit(self, section: str, drugs: Sequence[str], spec_overrides: Dict,
               base_spec: Dict, label: Optional[str] = None,
               warnings: Optional[List[str]] = None) -> Report:
        """Enqueue 2 x len(drugs) arms and return the (queued) report."""
        with self._lock:
            self._counter += 1
            report_id = f"rep_{self._counter:05d}"
            rep = Report(report_id, section, drugs, base_spec, spec_overrides,
                         label)
            rep.warnings = list(warnings or [])

        # Every arm is enqueued BEFORE the report becomes visible. A report that
        # is half-registered would look, to :meth:`refresh`, like a report whose
        # missing drugs have no pending arms — i.e. finished — and would be
        # assembled from an incomplete panel. A partial panel is not a ranking.
        for drug in rep.drugs:
            real = dict(base_spec, drug=drug, resistance_off=False)
            cf = dict(base_spec, drug=drug, resistance_off=True)
            rep.arms[drug] = {
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
        # Defensive: never assemble a panel whose arms are not all registered.
        expected = {(d, a) for d in rep.drugs
                    for a in ("real", "counterfactual")}
        have = {(d, a) for d, m in (rep.arms or {}).items() for a in m}
        if have != expected:
            return rep

        arms: Dict[str, Dict[str, Dict]] = {}
        pending = 0
        for drug in rep.drugs:
            arms[drug] = {}
            for arm, rid in (rep.arms.get(drug) or {}).items():
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
                        f"arm {arm!r} for {drug!r} (run {rid}) was cancelled; "
                        "a partial panel is not a ranking, so no report was "
                        "assembled"
                    )
                    rep.finished_at = time.time()
                    self._persist(rep)
                    return rep
                if job.status != "done":
                    pending += 1
                    continue
                arms[drug][arm] = job.result_view()
        if pending:
            rep.status = "running"
            return rep

        try:
            rep.body = assemble(rep.section, rep.spec_base, rep.drugs, arms)
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
table.tbl tbody tr.top td { background: var(--accent-bg); }
table.tbl tbody tr.top td.rank { color: var(--accent); }
table.tbl td .u { color: var(--fg-3); font-size: var(--fs-micro); margin-left: 3px; }
.rank { font-family: var(--font-mono); color: var(--fg-2); }
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
/* The ranking's own magnitude bar: unsigned, zero-anchored, one direction. */
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


def _md_bold(text: str) -> str:
    """Escape, then promote the attribution helpers' ``**x**`` to <b>x</b>."""
    out = _esc(text)
    parts = out.split("**")
    if len(parts) % 2 == 0:  # unbalanced — leave it escaped and literal
        return out
    return "".join(p if i % 2 == 0 else f"<b>{p}</b>"
                   for i, p in enumerate(parts))


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


def _ranking_table(rows: Sequence[Dict]) -> str:
    """The ranked panel.

    The ``pen. index`` column that used to sit here is GONE. ``penetration_index
    = mean(AUC)/p99(AUC)`` is a scale-invariant SHAPE statistic — it is
    unchanged when the whole exposure field is multiplied by any constant — so
    on real tissue it reports the vascular geometry of the section, which every
    drug in the panel shares (results/delivery_metrics.md §4: it spans 0.3044 to
    0.5107 across a panel whose delivered exposure spans 158x). Its failure mode
    used to be on display in this very table: temozolomide carries the HIGHEST
    penetration index of the GBM panel (0.5107) while ranking LAST. Its
    replacement, ``therapeutic_coverage`` at 1x IC50, was already computed and
    stored for every row and simply never rendered; it scales with delivered
    dose, it knows the drug's IC50, and on the same panel its separation
    (across-drug spread / grid drift) is 67.7 against 15.2 for the index.
    """
    # The column-group row answers, structurally, the question a reader of a
    # simulation output has to be able to answer before trusting any cell: which
    # of these did the model PRODUCE, and which was it HANDED? Left to prose it
    # is answered nowhere; left to the reader it is answered wrongly.
    group = (
        "<tr class='grp'><th colspan='2'></th>"
        "<th class='g' colspan='4'>computed by this run</th>"
        "<th class='g' colspan='5'>inputs it was given</th>"
        "<th class='g'>attribution</th></tr>"
    )
    head = (
        "<tr><th class='num'>#</th><th>drug</th><th class='num'>predicted kill</th>"
        "<th>ranked by kill</th>"
        "<th class='num'>tumour-cell AUC</th>"
        "<th class='num'>coverage @1x IC50</th>"
        "<th class='num'>IC50</th>"
        "<th>IC50 source</th><th class='num'>BBB</th><th class='num'>D_eff</th>"
        "<th class='num'>L</th>"
        "<th>dominant driver</th></tr>"
    )
    kills = [k for k in (_num(r.get("kill"), 0.0) for r in rows) if k > 0]
    top_kill = max(kills) if kills else 0.0
    body = []
    for r in rows:
        prov = r.get("provenance") or {}
        ic50 = prov.get("ic50") or {}
        kind = ic50.get("provenance_kind") or "assumed"
        k = _num(r.get("kill"), 0.0)
        frac = (max(0.0, k) / top_kill * 100.0) if top_kill > 0 else 0.0
        body.append(
            "<tr class='{cls}'>"
            "<td class='num rank'>{rank}</td>"
            "<td class='drug'>{name}</td>"
            "<td class='num'><b>{kill}</b></td>"
            "<td><span class='mag'><i style='width:{frac:.2f}%'></i></span></td>"
            "<td class='num'>{auc}<span class='u'>uM.s</span></td>"
            "<td class='num'>{cov}</td>"
            "<td class='num'>{ic}<span class='u'>uM</span></td>"
            "<td><span class='chip {kind}'>{kind}</span> {src}</td>"
            "<td class='num'>{bbb}</td>"
            "<td class='num'>{deff}<span class='u'>um2/s</span></td>"
            "<td class='num'>{L}<span class='u'>um</span></td>"
            "<td>{dom}</td></tr>".format(
                cls="top" if r.get("rank") == 1 else "",
                rank=r.get("rank"), name=_esc(r.get("name")),
                kill=_fmt(r.get("kill")), frac=frac,
                auc=_sig(r.get("tumor_auc")),
                ic=_sig(r.get("ic50")), kind=_esc(kind),
                src=_esc(ic50.get("source")), bbb=_fmt(r.get("bbb"), 3),
                deff=_fmt(r.get("d_eff"), 0), L=_fmt(r.get("L_um"), 0),
                cov=_fmt(r.get("coverage_1x"), 4), dom=_esc(r.get("dominant")),
            )
        )
    return (f"<div class='tblwrap'><table class='tbl'><thead>{group}{head}</thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>")


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


def _attribution_table(rows: Sequence[Dict]) -> str:
    head = (
        "<tr><th class='num'>#</th><th>drug</th><th class='num'>delivery</th>"
        "<th class='num'>potency</th><th class='num'>O2-resist.</th>"
        "<th class='num'>net</th><th>net (decades vs panel mean)</th>"
        "<th class='num'>rho</th><th class='num'>O2 cost</th></tr>"
    )
    body = []
    for r in rows:
        body.append(
            "<tr><td class='num rank'>{rank}</td><td class='drug'>{name}</td>"
            "<td class='num'>{d}</td><td class='num'>{p}</td>"
            "<td class='num'>{h}</td><td class='num'><b>{t}</b></td>"
            "<td>{bar}</td><td class='num'>{rho}</td>"
            "<td class='num'>{cost}<span class='u'>%</span></td></tr>".format(
                rank=r.get("rank"), name=_esc(r.get("name")),
                d=f"{r.get('delivery_decades', 0):+.2f}",
                p=f"{r.get('potency_decades', 0):+.2f}",
                h=f"{r.get('hypoxia_decades', 0):+.2f}",
                t=f"{r.get('total_decades', 0):+.2f}",
                bar=_bar(r.get("total_decades")),
                rho=_fmt(r.get("rho"), 4),
                cost=_fmt((r.get("hypoxia_cost") or 0) * 100, 1),
            )
        )
    return (f"<div class='tblwrap'><table class='tbl'><thead>{head}</thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>")


def _drug_cards(rows: Sequence[Dict], window_s: float) -> str:
    out = []
    window = _window_label(window_s)
    for r in rows:
        prov = r.get("provenance") or {}
        ic50 = prov.get("ic50") or {}
        mgmt = prov.get("mgmt") or {}
        mgmt_cls = "crit" if str(mgmt.get("state", "")).startswith("unmodelled") \
            else "measured" if mgmt.get("state") == "modelled" else "assumed"
        ub = prov.get("sensitivity_upper_bound")
        out.append(f"""
<div class="card">
  <div class="ch">
    <span class="rank">#{r.get('rank')}</span>
    <span class="n">{_esc(r.get('name'))}</span>
    <span class="chip accent">{_esc(r.get('dominant'))}</span>
    <span class="chip {'crit' if ub else 'computed'}">
      {'kill = UPPER BOUND' if ub else 'kill = computed'}</span>
    <span class="chip {mgmt_cls}">MGMT: {_esc(mgmt.get('state'))}</span>
  </div>
  <ul class="pts">
    <li><b>Delivery {r.get('delivery_decades', 0):+.2f} decades.</b>
      Mean exposure at this tumour's own TUMOR-cell positions is
      <span class="mono">{_sig(r.get('tumor_auc'))} uM.s</span> over the
      {window} window (time-averaged C = {_sig(r.get('c_tumor_uM'))} uM),
      the {_fmt(r.get('exposure_pct'), 0)}th percentile of this panel.
      D_eff = {_fmt(r.get('d_eff'), 0)} um2/s and BBB penetrance
      {_fmt(r.get('bbb'), 3)} give an analytic penetration length
      L = {_fmt(r.get('L_um'), 0)} um
      ({_fmt((r.get('L_um') or 0) / 45.0, 1)}x the calibrated ~45 um
      perivascular boundary layer).</li>
    <li><b>Potency {r.get('potency_decades', 0):+.2f} decades.</b>
      IC50 = {_sig(r.get('ic50'))} uM —
      <span class="chip {_esc(ic50.get('provenance_kind') or 'assumed')}">{_esc(ic50.get('provenance_kind'))}</span>
      {_esc(ic50.get('source'))} — the {_fmt(r.get('potency_pct'), 0)}th
      percentile of the panel. {_esc(ic50.get('note'))}</li>
    <li><b>O2-resistance {r.get('hypoxia_decades', 0):+.2f} decades.</b>
      COUNTERFACTUAL, not asserted: re-simulating this same drug on this same
      tissue with the hypoxic AND quiescent kill penalties switched off
      (run <span class="mono">{_esc(prov.get('counterfactual_run_id'))}</span>)
      moves the predicted kill from {_fmt(r.get('kill'))} to
      {_fmt(r.get('kill_no_resistance'))}, so
      rho = {_fmt(r.get('rho'), 4)} and this tumour's oxygen-driven resistance
      state costs this drug {_fmt((r.get('hypoxia_cost') or 0) * 100, 1)}% of
      the kill it would otherwise achieve.</li>
    <li><b>Rank gap.</b> {_md_bold(r.get('gap_reason') or 'single-drug panel: no adjacent rank to compare against.')}</li>
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
  </dl>
</div>""")
    return "\n".join(out)


def render_html(doc: Dict) -> str:
    """Render a report document as ONE self-contained HTML page.

    Self-contained is a hard requirement, not a nicety: the file must render
    identically off a USB stick, from an email attachment, or inside an air-
    gapped review. There is therefore NO ``<img>``, no ``@import``, no
    ``url(...)``, no webfont, no external script and no fetch of any kind — the
    Studio's own ``tokens.css`` is inlined and every graphic is a styled element.
    """
    rows = doc.get("ranking") or []
    limits = doc.get("limits") or []
    attribution = doc.get("attribution") or {}
    check = attribution.get("check") or {}
    micro = doc.get("microenvironment") or {}
    spec = doc.get("spec_base") or {}
    section = doc.get("section_detail") or {}
    window_s = float(doc.get("window_s") or 0.0)
    best = rows[0] if rows else None

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

    # The hero names the top of a RANKING, which is the defensible output — but
    # it may not name it as settled while limit 2b on this same page says the
    # ordering is law-dependent. The law is the one stored with this document
    # (read from PDEConfig when the runs were made), falling back to the shipped
    # default for a document that predates it. See app/server/delivery_law.py.
    law = doc.get("delivery_law")
    if not isinstance(law, dict):
        law = delivery_law()
    # The law limit is on this same page; the hero's condition used to cite it
    # as "limit 2b", which named nothing a reader could find. It now links to
    # the limit's own anchor — an in-page fragment, not an external reference.
    law_ref = "#lim-transvascular_delivery_law"
    headline = ""
    if best:
        headline = f"""
<div class="hero">
  <div class="pick">
    <div class="eyebrow">the finding</div>
    <div>
      <div class="pn">Top of this ranking, under {_esc(law.get('law'))}:
        <b>{_esc(best.get('name'))}</b>
        <span class="under">predicted kill {_fmt(best.get('kill'))} — read the
          ORDER, not this level</span></div>
      <div class="pchips">
        <span class="chip crit">UPPER BOUND</span>
        <span class="chip crit">LAW-DEPENDENT</span>
        <span class="chip crit">NOT A RECOMMENDATION</span>
        <span class="chip accent">{_esc(best.get('dominant'))}</span>
      </div>
      <div class="pfacts">L = {_fmt(best.get('L_um'), 0)} um &middot;
        IC50 {_sig(best.get('ic50'))} uM
        ({_esc(((best.get('provenance') or {}).get('ic50') or {}).get('source'))})
        &middot; tumour-cell exposure {_sig(best.get('tumor_auc'))} uM.s</div>
    </div>
  </div>
  <div class="cond">
    <div class="cl">what it is not</div>
    <div>
      <p><b>This is a model output, not a recommendation.</b> It has never been
         checked against a measured response for this tissue or any other. The
         ORDER is the claim; the LEVEL is not. Treat it as a hypothesis to test
         in the lab, nothing more.</p>
    </div>
  </div>
  <div class="cond">
    <div class="cl">the condition on the order</div>
    <div>
      <p><b>The order carries a condition — read this with the name above.</b>
         Every delivery number here was produced under
         <span class="mono">{_esc(law.get('law'))}</span>; under the alternative
         source term <span class="mono">{_esc(law.get('alternative_law'))}</span>
         the panel is re-weighted by permeability (which spans 47x across this
         drug list) and <b>the top-1 GBM pick changes, cisplatin &rarr;
         vincristine</b>, swapping 4 of 45 pairs (<a href="{law_ref}">see the
         delivery-law limit above</a>;
         results/transvascular_exchange_revalidation.md §4.1).</p>
      <p>Independently of the law, results/ranking_robustness.md §8 concludes
         that until an ex-vivo potency readout exists on this patient's own
         cells, <b>"the reportable output of this pipeline for the top of the
         panel is a SET, not a winner"</b> — P(cisplatin = rank 1) is 0.931
         under that report's declared priors and falls to 0.587 in its stress
         scenario. So: this is the top of the ranking one law and one set of
         priors produced, not this specimen's drug.</p>
    </div>
  </div>
</div>"""

    caveat_items = "".join(f"<li>{_esc(c)}</li>" for c in (doc.get("caveats") or []))
    warn_items = "".join(f"<li>{_esc(c)}</li>" for c in (doc.get("warnings") or []))

    sp = check.get("spearman_total_vs_kill")

    # Section 05 renders micro["patient_id"] under the label "specimen id". The
    # divergence is deliberate: the field name is load-bearing across the engine
    # and the API payloads and renaming it is a data-model change, but a field
    # labelled "patient id" on a rendered page invites a user to type one, and
    # this tool is scoped to de-identified tissue (docs/DEPLOYMENT.md §0). If the
    # field is ever renamed, the label here is already correct.
    body = f"""
<div class="doc">
<header class="hd">
  <div class="kicker">tmesim · simulated compound ranking</div>
  <h1>Simulated compound ranking — {_esc(section.get('label') or doc.get('section'))}</h1>
  <div class="sub">A ranked, mechanistically attributed drug panel simulated on
    one measured specimen. Every number below is a prediction; the sections are
    ordered so that the limits on it are read first.</div>
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
<h2><span class="sn">02</span><span class="st">Ranked panel</span></h2>
<div class="full">{headline}</div>
{_provenance_key()}
<div class="full">{_ranking_table(rows)}</div>
<p>All {len(rows)} drugs were simulated on the <b>same</b> reconstruction, with
the same {_window_label(window_s)} exposure window, the same seed
({_esc(((rows[0].get('provenance') if rows else {}) or {}).get('seed'))}) and the
same compilation pipeline (<span class="mono">calibrate_alpha &rarr;
calibrated_drugparams &rarr; apply_barrier_to_drugparams</span>), so every
difference is attributable to that drug's own physics and potency on identical
tissue. Ranked by predicted kill fraction; the bar is that kill as a fraction of
the top row's, and like every level on this page it is the centre of a wide
interval.</p>
<p><b>On the delivery column.</b>
<span class="mono">coverage @1x IC50</span> is
<span class="mono">therapeutic_coverage</span> — the fraction of this section's
tumour cells whose time-averaged exposure exceeds that drug's own IC50 over the
{_window_label(window_s)} window. It replaces the
<span class="mono">pen. index</span> column earlier versions of this report
printed. <span class="mono">penetration_index = mean(AUC)/p99(AUC)</span> is a
scale-invariant SHAPE statistic: multiply the whole exposure field by any
constant and it does not move, so on real tissue it reports the section's
vascular geometry, which every drug here shares. It was retired for that reason
(results/delivery_metrics.md §4) and the panel demonstrates why — temozolomide
holds the HIGHEST penetration index of the GBM panel while ranking LAST.
Coverage is not exempt from anything else on this page: it is a delivery number
like the others, so it is produced under the delivery law named above and it
carries the same mesh term as everything else.</p>
</section>

<section class="sec">
<h2><span class="sn">03</span><span class="st">Why each drug ranks where it does</span></h2>
<p>Each drug's standing is decomposed into three <b>exactly additive log10
terms</b> ("decades"), each centred on the panel mean so <span class="mono">+</span>
reads as an advantage over the panel and <span class="mono">−</span> as a penalty.
Nothing is asserted: the O2-resistance term comes from a real counterfactual
simulation, which is why every drug here was simulated twice.</p>
<div class="full"><pre class="code">{_esc(attribution.get('method'))}</pre>
{_attribution_table(rows)}</div>
<p><b>Is the attribution itself faithful?</b> Spearman rank correlation across
the panel between the net attribution and the simulated kill:
<span class="mono hl">rho = {_fmt(sp, 3)}</span> &rarr;
{_esc(check.get('faithful'))}. Component-wise, rho(delivery, kill) =
{_fmt(check.get('spearman_delivery_vs_kill'), 3)} and rho(potency, kill) =
{_fmt(check.get('spearman_potency_vs_kill'), 3)}. With only {len(rows)} drugs a
single rho is a weak statistic — it is reported to be checkable, not because it
is powerful.</p>
</section>

<section class="sec">
<h2><span class="sn">04</span><span class="st">Per-drug reading and provenance</span></h2>
<div class="full">{_drug_cards(rows, window_s)}</div>
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
  <li><b>A tissue with a paired, measured drug response.</b> Everything here is
    unfalsifiable until a tumour with a known outcome runs through the identical
    pipeline and the predicted ranking is scored against the observed one.</li>
  <li><b>A potency measurement in context.</b> Cell-line IC50 must be replaced by
    (or at minimum reconciled with) a 3-D / organoid / patient-derived
    measurement from tissue resembling this one.</li>
  <li><b>More than one penetration anchor.</b> alpha rides on a single
    doxorubicin measurement; independent per-drug perivascular-distance data
    would test whether the ordering transfers outside the strongly-bound
    class.</li>
  <li><b>A measured per-cell sensitivity</b>, replacing today's explicit
    assumption (uniform 1.0).</li>
  <li><b>A panel that can measure cell state.</b> Until hypoxia and
    proliferation are measurable on the section itself, the absolute kill has no
    defensible level under any state route.</li>
</ul>
</section>

{_banner('bottom')}

<div class="foot">
  Generated by {_esc(doc.get('generator'))} ·
  attribution: {_esc(doc.get('attribution_implementation'))} ·
  runs: {_esc(', '.join(doc.get('run_ids') or []))} ·
  deterministic (seed
  {_esc(((rows[0].get('provenance') if rows else {}) or {}).get('seed'))}).
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
        # Matches the <h1> ("Simulated compound ranking"). It read "tmesim
        # therapy report" until 2026-08-15 — two renames stale, and naming a
        # deliverable a "therapy report" is the implication this product does
        # not make. NOT FOR CLINICAL USE stays IN the title on purpose: this
        # document travels off a USB stick, and the browser tab is the one part
        # of it that survives a screenshot of the window frame.
        f"<title>Codon simulated compound ranking {_esc(doc.get('report_id'))} — "
        f"{_esc(doc.get('section'))} — NOT FOR CLINICAL USE</title>\n"
        f"<style>\n{_tokens_css()}\n{_DOC_CSS}\n{_print_css()}\n</style>\n"
        "</head><body>\n" + body + "\n</body></html>\n"
    )
