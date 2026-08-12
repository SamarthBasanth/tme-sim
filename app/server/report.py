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
    describe_section,
    lookup_section,
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
    path = os.path.join(repo_root(), "app", "reports")
    os.makedirs(path, exist_ok=True)
    return path


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
    if isinstance(dx, (int, float)):
        converged = float(dx) <= 15.0
        dx_body = (
            f"This report ran at dx = {float(dx):g} um. kill_fraction is "
            "grid-converged for dx <= 20 um; penetration_index is NOT converged "
            "above ~15 um and is over-estimated on coarse grids because a mesh "
            "coarser than the ~45 um perivascular boundary layer cannot resolve "
            "it (0.285 at dx = 60 um falling to 0.187 at dx = 10 um)."
        ) + ("" if converged else
             " THIS RUN IS OUTSIDE THE CONVERGED PENETRATION REGIME; read every "
             "penetration number as an upper bound." if float(dx) <= 20.0 else
             " THIS RUN IS OUTSIDE THE CONVERGED REGIME FOR BOTH kill_fraction "
             "AND penetration_index.")
        limits.append({
            "id": "resolution",
            "severity": "info" if converged else "warn",
            "applies": True,
            "title": "Resolution sensitivity",
            "source": "results/numerical_rigor.md §2",
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
                "reason": f"{type(exc).__name__}: {exc}"}


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
        os.makedirs(self.root, exist_ok=True)
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
        os.makedirs(d, exist_ok=True)
        doc = rep.full_view(self.progress(rep))
        tmp = os.path.join(d, "report.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(jsonable(doc), fh, indent=1)
        os.replace(tmp, os.path.join(d, "report.json"))
        if rep.status == "done":
            tmp = os.path.join(d, "report.html.tmp")
            with open(tmp, "w") as fh:
                fh.write(render_html(doc))
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
            rep.error = f"assembly failed: {type(exc).__name__}: {exc}"
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


_DOC_CSS = """
/* ---- report document ---------------------------------------------------
   Same rules as the Studio: square corners (<=2px), mono tabular numerics,
   one accent, semantic colour reserved for state. Nothing in this stylesheet
   references anything off this page - no stylesheet import, no resource
   reference, no webfont, no image, no script. Every graphic is a styled
   element, and the whole document is one file. */
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0; background: var(--bg-0); color: var(--fg-1);
  font-family: var(--font-ui); font-size: var(--fs-sm); line-height: 1.45;
}
.doc { max-width: 1180px; margin: 0 auto; padding: 0 18px 56px; }
.mono, .num, code, pre, td.num, th.num {
  font-family: var(--font-mono); font-variant-numeric: tabular-nums;
}
h1, h2, h3 { color: var(--fg-0); font-weight: 600; margin: 0; }
h1 { font-size: var(--fs-xl); letter-spacing: 0.01em; }
h2 {
  font-size: var(--fs-md); text-transform: uppercase; letter-spacing: 0.06em;
  color: var(--fg-2); border-bottom: 1px solid var(--line-2);
  padding-bottom: 5px; margin: 30px 0 12px;
}
h3 { font-size: var(--fs-md); margin: 18px 0 6px; }
p { margin: 6px 0; }
a { color: var(--accent); }
b, strong { color: var(--fg-0); font-weight: 600; }

/* ---- the banner. Top and bottom, unmissable, never a footnote. --------- */
.nfcu {
  border: 2px solid var(--crit); background: var(--crit-bg);
  padding: 11px 14px; margin: 16px 0; border-radius: var(--radius);
}
.nfcu .t {
  font-family: var(--font-mono); font-size: var(--fs-md); font-weight: 700;
  letter-spacing: 0.16em; color: var(--crit); display: block; margin-bottom: 5px;
}
.nfcu .b { color: var(--fg-0); }

header.hd { padding: 20px 0 0; border-bottom: 1px solid var(--line-1); }
.hd .sub { color: var(--fg-2); font-size: var(--fs-xs); margin-top: 4px; }
.meta {
  display: flex; flex-wrap: wrap; gap: 3px 20px; padding: 9px 0 14px;
  font-size: var(--fs-xs);
}
.meta .k { color: var(--fg-3); text-transform: uppercase; letter-spacing: 0.06em;
  margin-right: 5px; }
.meta .v { color: var(--fg-0); font-family: var(--font-mono); }

/* ---- limits ------------------------------------------------------------ */
.limit {
  border-left: 3px solid var(--line-2); background: var(--bg-1);
  padding: 8px 12px; margin: 0 0 3px;
}
.limit.critical { border-left-color: var(--crit); background: var(--crit-bg); }
.limit.warn { border-left-color: var(--warn); background: var(--warn-bg); }
.limit.info { border-left-color: var(--info); background: var(--info-bg); }
.limit .lt { color: var(--fg-0); font-weight: 600; display: block; }
.limit .lb { margin-top: 3px; }
.limit .ls {
  display: block; margin-top: 4px; font-family: var(--font-mono);
  font-size: var(--fs-micro); color: var(--fg-3);
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
  color: var(--fg-2); white-space: nowrap;
}
.chip.measured { color: var(--ok); border-color: var(--ok); }
.chip.assumed { color: var(--warn); border-color: var(--warn); }
.chip.computed { color: var(--info); border-color: var(--info); }
.chip.crit { color: var(--crit); border-color: var(--crit); }
.chip.accent { color: var(--accent); border-color: var(--accent); }

/* ---- tables ------------------------------------------------------------ */
.tblwrap { overflow-x: auto; max-width: 100%; margin: 8px 0; }
table.tbl { border-collapse: collapse; width: 100%; font-size: var(--fs-sm); }
table.tbl th {
  text-align: left; font-weight: 600; color: var(--fg-2);
  font-size: var(--fs-micro); text-transform: uppercase; letter-spacing: 0.06em;
  border-bottom: 1px solid var(--line-2); padding: 5px 8px; white-space: nowrap;
}
table.tbl td {
  padding: 4px 8px; border-bottom: 1px solid var(--line-0);
  vertical-align: top;
}
table.tbl td.num, table.tbl th.num { text-align: right; }
table.tbl tbody tr:nth-child(even) td { background: var(--bg-1); }
table.tbl tbody tr.top td { background: var(--accent-bg); }
table.tbl td .u { color: var(--fg-3); font-size: var(--fs-micro); margin-left: 3px; }
.rank { font-family: var(--font-mono); color: var(--fg-2); }
.drug { color: var(--fg-0); font-weight: 600; }

/* ---- attribution bars (inline, no images) ------------------------------ */
.bar { position: relative; height: 9px; background: var(--bg-3); min-width: 96px; }
.bar > i { position: absolute; top: 0; bottom: 0; display: block; }
.bar > i.pos { background: var(--accent); }
.bar > i.neg { background: var(--fg-3); }
.bar > u { position: absolute; top: -2px; bottom: -2px; left: 50%; width: 1px;
  background: var(--line-2); }

/* ---- misc -------------------------------------------------------------- */
pre.code {
  background: var(--bg-1); border: 1px solid var(--line-1); padding: 9px 11px;
  overflow-x: auto; font-size: var(--fs-xs); color: var(--fg-1); margin: 8px 0;
}
.kv { display: grid; grid-template-columns: auto 1fr; gap: 0 12px;
  font-size: var(--fs-xs); margin: 5px 0; }
.kv dt { color: var(--fg-3); font-size: var(--fs-micro); text-transform: uppercase;
  letter-spacing: 0.06em; padding: 2px 0; border-bottom: 1px solid var(--line-0); }
.kv dd { margin: 0; color: var(--fg-0); font-family: var(--font-mono);
  padding: 2px 0; border-bottom: 1px solid var(--line-0); word-break: break-word; }
.card { border: 1px solid var(--line-1); background: var(--bg-1); padding: 10px 12px;
  margin: 8px 0; }
.card .ch { display: flex; flex-wrap: wrap; align-items: baseline; gap: 8px;
  border-bottom: 1px solid var(--line-1); padding-bottom: 5px; margin-bottom: 6px; }
.card .ch .n { color: var(--fg-0); font-weight: 600; font-size: var(--fs-md); }
ul.pts { margin: 5px 0; padding-left: 17px; }
ul.pts li { margin: 3px 0; }
.foot { margin-top: 22px; padding-top: 10px; border-top: 1px solid var(--line-1);
  color: var(--fg-3); font-size: var(--fs-micro); font-family: var(--font-mono); }
.hl { color: var(--fg-0); }
"""


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


def _limits_html(limits: Sequence[Dict]) -> str:
    out = []
    for lim in limits:
        cls = lim.get("severity") or "info"
        if not lim.get("applies", True):
            cls += " na"
        src = (f'<span class="ls">{_esc(lim["source"])}</span>'
               if lim.get("source") else "")
        out.append(
            f'<div class="limit {cls}">'
            f'<span class="lt">{_esc(lim.get("title"))}</span>'
            f'<div class="lb">{_esc(lim.get("body"))}</div>{src}</div>'
        )
    return "\n".join(out)


def _banner(where: str) -> str:
    return (
        f'<div class="nfcu" role="alert" data-where="{_esc(where)}">'
        f'<span class="t">NOT FOR CLINICAL USE</span>'
        f'<div class="b">{_esc(NOT_FOR_CLINICAL_USE)}</div></div>'
    )


def _ranking_table(rows: Sequence[Dict]) -> str:
    head = (
        "<tr><th class='num'>#</th><th>drug</th><th class='num'>predicted kill</th>"
        "<th class='num'>tumour-cell AUC</th><th class='num'>IC50</th>"
        "<th>IC50 source</th><th class='num'>BBB</th><th class='num'>D_eff</th>"
        "<th class='num'>L</th><th class='num'>pen. index</th>"
        "<th>dominant driver</th></tr>"
    )
    body = []
    for r in rows:
        prov = r.get("provenance") or {}
        ic50 = prov.get("ic50") or {}
        kind = ic50.get("provenance_kind") or "assumed"
        body.append(
            "<tr class='{cls}'>"
            "<td class='num rank'>{rank}</td>"
            "<td class='drug'>{name}</td>"
            "<td class='num'><b>{kill}</b></td>"
            "<td class='num'>{auc}<span class='u'>uM.s</span></td>"
            "<td class='num'>{ic}<span class='u'>uM</span></td>"
            "<td><span class='chip {kind}'>{kind}</span> {src}</td>"
            "<td class='num'>{bbb}</td>"
            "<td class='num'>{deff}<span class='u'>um2/s</span></td>"
            "<td class='num'>{L}<span class='u'>um</span></td>"
            "<td class='num'>{pen}</td>"
            "<td>{dom}</td></tr>".format(
                cls="top" if r.get("rank") == 1 else "",
                rank=r.get("rank"), name=_esc(r.get("name")),
                kill=_fmt(r.get("kill")), auc=_sig(r.get("tumor_auc")),
                ic=_sig(r.get("ic50")), kind=_esc(kind),
                src=_esc(ic50.get("source")), bbb=_fmt(r.get("bbb"), 3),
                deff=_fmt(r.get("d_eff"), 0), L=_fmt(r.get("L_um"), 0),
                pen=_fmt(r.get("pen"), 3), dom=_esc(r.get("dominant")),
            )
        )
    return (f"<div class='tblwrap'><table class='tbl'><thead>{head}</thead>"
            f"<tbody>{''.join(body)}</tbody></table></div>")


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

    meta_items = [
        ("report", doc.get("report_id")),
        ("section", doc.get("section")),
        ("panel", section.get("panel")),
        ("assay", section.get("assay")),
        ("accession", section.get("accession")),
        ("cells", micro.get("n_cells")),
        ("tumour cells", micro.get("n_tumor")),
        ("grid", " x ".join(str(g) for g in (micro.get("grid") or []))),
        ("dx", f"{_fmt(spec.get('dx_um'), 1)} um"),
        ("window", _window_label(window_s)),
        ("sensitivity mode", doc.get("sensitivity_mode")),
        ("vascular mode", spec.get("vascular_mode")),
        ("oxygen mode", spec.get("oxygen_mode")),
        ("state mode", spec.get("state_mode")),
        ("outcome label", "ABSENT" if not section.get("has_outcome_label") else "present"),
        ("generated", time.strftime("%Y-%m-%d %H:%M:%S",
                                    time.localtime(doc.get("finished_at")
                                                   or time.time()))),
    ]
    meta = "".join(
        f'<span><span class="k">{_esc(k)}</span>'
        f'<span class="v">{_esc(v if v not in (None, "") else "—")}</span></span>'
        for k, v in meta_items)

    headline = ""
    if best:
        headline = f"""
<div class="card">
  <div class="ch">
    <span class="n">Highest predicted kill on this specimen: {_esc(best.get('name'))}</span>
    <span class="chip accent">{_fmt(best.get('kill'))}</span>
    <span class="chip crit">UPPER BOUND</span>
  </div>
  <p>Dominant driver <b>{_esc(best.get('dominant'))}</b> ·
     L = {_fmt(best.get('L_um'), 0)} um ·
     IC50 {_sig(best.get('ic50'))} uM
     ({_esc(((best.get('provenance') or {}).get('ic50') or {}).get('source'))}) ·
     tumour-cell exposure {_sig(best.get('tumor_auc'))} uM.s.</p>
  <p><b>This is a model output, not a recommendation.</b> It has never been
     checked against a measured response for this tissue or any other. The
     ORDER is the claim; the LEVEL is not. Treat it as a hypothesis to test in
     the lab, nothing more.</p>
</div>"""

    caveat_items = "".join(f"<li>{_esc(c)}</li>" for c in (doc.get("caveats") or []))
    warn_items = "".join(f"<li>{_esc(c)}</li>" for c in (doc.get("warnings") or []))

    sp = check.get("spearman_total_vs_kill")
    body = f"""
<div class="doc">
<header class="hd">
  <h1>Virtual therapy panel — {_esc(section.get('label') or doc.get('section'))}</h1>
  <div class="sub">Ranked, mechanistically attributed drug panel simulated on one
    measured specimen. Report {_esc(doc.get('report_id'))}, status
    {_esc(doc.get('status'))}.</div>
  <div class="meta">{meta}</div>
</header>

{_banner('top')}

<h2>1 · Limits — read these before any number on this page</h2>
{_limits_html(limits)}

<h2>2 · Ranked panel</h2>
{headline}
{_ranking_table(rows)}
<p>All {len(rows)} drugs were simulated on the <b>same</b> reconstruction, with
the same {_window_label(window_s)} exposure window, the same seed
({_esc(((rows[0].get('provenance') if rows else {}) or {}).get('seed'))}) and the
same compilation pipeline (<span class="mono">calibrate_alpha &rarr;
calibrated_drugparams &rarr; apply_barrier_to_drugparams</span>), so every
difference is attributable to that drug's own physics and potency on identical
tissue. Ranked by predicted kill fraction.</p>

<h2>3 · Why each drug ranks where it does</h2>
<p>Each drug's standing is decomposed into three <b>exactly additive log10
terms</b> ("decades"), each centred on the panel mean so <span class="mono">+</span>
reads as an advantage over the panel and <span class="mono">−</span> as a penalty.
Nothing is asserted: the O2-resistance term comes from a real counterfactual
simulation, which is why every drug here was simulated twice.</p>
<pre class="code">{_esc(attribution.get('method'))}</pre>
{_attribution_table(rows)}
<p><b>Is the attribution itself faithful?</b> Spearman rank correlation across
the panel between the net attribution and the simulated kill:
<span class="mono hl">rho = {_fmt(sp, 3)}</span> &rarr;
{_esc(check.get('faithful'))}. Component-wise, rho(delivery, kill) =
{_fmt(check.get('spearman_delivery_vs_kill'), 3)} and rho(potency, kill) =
{_fmt(check.get('spearman_potency_vs_kill'), 3)}. With only {len(rows)} drugs a
single rho is a weak statistic — it is reported to be checkable, not because it
is powerful.</p>

<h2>4 · Per-drug reading and provenance</h2>
{_drug_cards(rows, window_s)}

<h2>5 · Specimen</h2>
<dl class="kv">
  <dt>patient / sample id</dt><dd>{_esc(micro.get('patient_id'))}</dd>
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

<h2>6 · Caveats carried by the individual runs</h2>
<ul class="pts">{caveat_items or '<li>none reported</li>'}</ul>
{('<h3>Spec warnings at submission</h3><ul class="pts">' + warn_items + '</ul>') if warn_items else ''}

<h2>7 · What would have to be true for this to be trustworthy</h2>
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
        f"<title>tmesim therapy report {_esc(doc.get('report_id'))} — "
        f"{_esc(doc.get('section'))} — NOT FOR CLINICAL USE</title>\n"
        f"<style>\n{_tokens_css()}\n{_DOC_CSS}\n</style>\n"
        "</head><body>\n" + body + "\n</body></html>\n"
    )
