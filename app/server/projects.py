"""Projects, specimens, and the cross-specimen comparison.

The Studio was single-specimen and stateless between sessions: runs lived flat
under ``app/runs/<run_id>/`` with no notion of which patient they belonged to.
This module adds the missing layer, and only that layer:

* A **Project** is a study. It has a name, notes, and one or more **Specimens**.
* A **Specimen** binds a ``sections`` registry id to a *capability report* —
  captured at the moment it is added, from the section's own panel — plus the
  list of ``run_id``s that belong to it. Runs are **associated, never copied**:
  ``app/runs/`` remains the single store, and a specimen holds ids into it.
* ``compare_project`` compares the matched runs of a project **one compound at a
  time**: the same compound on specimen A against the same compound on specimen
  B, plus the tissue contrast between the two specimens. It emits no ordering.

Everything is persisted as one JSON document per project under
``app/projects/<project_id>/project.json`` (gitignored, like ``app/runs/``), so
a study survives a restart.

What this module compares, and what it withdrew
-----------------------------------------------

**Rank agreement was the primary result of this module and has been WITHDRAWN.**
Until 2026-08 ``compare_project`` returned a ``ranking`` block — a per-compound
ordering per specimen, a Kendall tau-b between specimens, and a
``RANKING AGREES / DISAGREES`` verdict — on the strength of
``results/heldout_patient.md`` §3.3 (the kill ordering replicated on two real GBM
sections) and ``results/cell_state_impact.md`` §3.3 (the ordering was invariant
across all three cell-state routes). Two later measurements took that away:

* ``results/pk_ranking_robustness.md`` — per-compound plasma exposure across the
  compounds the GBM panel ranks spans **113.1x**, and the published GBM ordering
  first breaks at a spread of only **1.53x**. A swap in the published order sits
  inside a single compound's own licensed dose range (irinotecan's FDA label
  publishes a 2.04x Cmax span; ``irinotecan > temozolomide`` needs 1.961x to
  reverse). The ordering the module presented was an ordering at roughly 9x
  supratherapeutic exposure.
* ``results/capacity_ranking_robustness.md`` — the engine assumes a tissue
  capacity of exactly 1.0 for every compound, against a measured rat brain-slice
  spread of 3.2 - 4184 mL/g. Under this repo's own structural per-compound
  capacity rule the **breast top-1 pick changes** (``vinorelbine`` ->
  ``doxorubicin``), and the full uniform-capacity order breaks at R = 36.3 on
  GBM and R = 5.28 on breast.

Ordering compounds by predicted effect is therefore gone from this module.
**Comparison of runs stays**, because it is a different question: the same
compound on two specimens is a legitimate comparison and nothing above touches
it — a per-compound plasma factor and a per-compound capacity term are the SAME
factor on both sides of a same-compound comparison. What survives is stated in
full by :data:`RELATIVE_EXPOSURE_CLAIM`, :data:`RELATIVE_EXPOSURE_BASIS` and
:data:`RELATIVE_EXPOSURE_LIMITS` below, and every surface that used to say "the
ranking is the defensible output" now says that instead.

The rest of the shape is unchanged and still earned:

* ``results/heldout_patient.md`` §6.5 — *"the response layer travels; the
  reconstruction layer does not"*. **Absolute kills do not travel.** Its §3.5
  names two mechanisms, neither of them a drug property: the absolute oxygen
  level is uncalibrated on both sections (``OxygenConfig.consumption`` was never
  fitted to a measured tension), and the ECM field is a max-normalised KDE over
  *whichever matrix markers the panel happens to carry* — 11 genes on one
  section, 3 on the other — so the two have no common absolute scale even on
  identical tissue. Its §5 table labels the ECM row "NOT comparable as biology".
  Both external reconstruction validations **failed** on the second patient: the
  ENDOTHELIAL fraction (9.74% -> 4.18% against a volume-EM 8.0% anchor) and the
  vascular area fraction (2.50% -> 5.79% against a measured 2.1 +/- 0.3%).
* ``results/cell_state_impact.md`` supplies the capability rule. Its §3.2
  measures absolute kill moving by up to **0.3252 absolute (1.77x relative)**
  purely with the choice of cell-state route, with tissue, compound, transport
  and RNG held identical. Its §5 conclusion is that a panel carrying no
  hypoxia-response transcript cannot support an absolute kill under *any* route
  — the default route's HYPOXIC band is empty for a bookkeeping reason, the
  marker route declares hypoxia UNMEASURED, and the vessel-distance route is a
  geometric proxy.

So this module refuses to present a cross-specimen absolute-kill delta as a
finding unless *both* specimens' panels can support a measured cell state *and*
both runs actually used a measured route. A measured kill is never silently
differenced against an assumed one. The primary result of a comparison is now
the specimen capability contrast and the composition contrast — what each panel
can support, and how the two tissues differ — not an ordering.
"""
from __future__ import annotations

import copy
import json
import math
import os
import re
import threading
import time
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from .cache import (
    PRIVATE_FILE_MODE,
    all_sections,
    describe_section,
    own_only,
    private_makedirs,
    redact_exception,
    repo_root,
    runs_dir,
    state_route_support,
)
from .delivery_law import delivery_law

__all__ = [
    "PROJECT_ID_RE",
    "RUN_ID_RE",
    "MAX_NAME_CHARS",
    "MAX_NOTES_CHARS",
    "ProjectError",
    "ProjectStore",
    "PROJECT_STORE",
    "projects_dir",
    "capability_report",
    "run_composition",
    "THUMBNAIL_NEVER_FABRICATES",
    "spec_signature",
    "compare_project",
    "RELATIVE_EXPOSURE_CLAIM",
    "RELATIVE_EXPOSURE_BASIS",
    "RELATIVE_EXPOSURE_LIMITS",
    "RELATIVE_EXPOSURE_FULL",
]

PROJECT_ID_RE = re.compile(r"^p_\d{5}$")
SPECIMEN_ID_RE = re.compile(r"^s_\d{5}$")

#: A run id, for the SAME reason ``PROJECT_ID_RE`` exists: the thumbnail reads
#: ``app/runs/<run_id>/result.json`` straight off disk, so a ``run_id`` is a path
#: segment there and is validated before it is ever joined to a root. The ids in
#: a project document arrive through ``attach_run``, which today refuses an id
#: the JobManager cannot resolve — but a document is a file, an operator can edit
#: it, and the traversal fixed in ``ProjectStore._dir`` was exactly a guard that
#: existed and was not wired to the boundary that needed it. ``\d{5,}`` and not
#: ``\d{5}``: ``JobManager.submit`` formats ``r_{counter:05d}``, which is five
#: digits only until the 100,000th run.
RUN_ID_RE = re.compile(r"^r_\d{5,}$")

# --------------------------------------------------------------------------- #
# FREE-TEXT LIMITS, and why a store that only checked for EMPTY was a DoS.
#
# ``create``/``update``/``duplicate``/``add_specimen``/``update_specimen``
# validated their text for emptiness and nothing else, so a name was refused at
# zero characters and accepted at twenty million. That is not a cosmetic gap,
# because of where the text GOES: ``_summary`` puts ``notes`` in every row of
# ``store.list()``, and ``store.list()`` is what ``GET /api/workspace`` returns
# -- the SINGLE boot call for the home page.
#
# Measured against a running server on 2026-09-20:
#
#     POST /api/projects  {"name": "pwn", "notes": "N" * 20_971_520}  -> 201
#     GET  /api/workspace                             163 B -> 20,972,002 B
#
# One accepted request, and the front door of the product serves a 21 MB boot
# payload to every visitor, from disk, across restarts, until someone finds the
# directory and deletes it by hand.
#
# THE BODY CAP DID NOT STOP IT, and that is the other half worth writing down:
# ``api.py::_cap_request_body`` refuses on ``Content-Length``, and a request
# sent ``Transfer-Encoding: chunked`` has none. That is a defect in a different
# file. These caps are deliberately NOT a restatement of it -- a request-size
# limit answers "how much may you send me", and this answers "how much am I
# willing to persist and then hand to every future reader". The second question
# is this module's to answer, and its answer must hold even when the first one
# is wrong.
#
# The numbers: a study name is a label, not a document -- 200 characters is
# longer than any real one and short enough that a thousand projects cost
# 200 KB. Notes are prose and people paste protocol text into them, so 20,000
# characters (~4 pages) is generous; 1,000 projects at the cap is 20 MB, which
# is the point at which the operator has a store problem rather than an
# attacker. Counted in CHARACTERS, not bytes: the limit has to be explicable to
# the person typing into the field, and ``len()`` on a ``str`` is the number a
# UI can enforce too.
# --------------------------------------------------------------------------- #
#: Longest project name / specimen label, in characters.
MAX_NAME_CHARS = 200

#: Longest project notes / specimen notes, in characters.
MAX_NOTES_CHARS = 20_000

#: The figures this module quotes, kept as named constants so the message and
#: the evidence cannot drift apart. Both are from results/cell_state_impact.md
#: §3.2 (kill across the three state routes, everything else held identical).
STATE_ROUTE_KILL_SWING_ABS = 0.3252
STATE_ROUTE_KILL_SWING_REL = 1.77

#: The claim that REPLACED "the cross-drug ranking is the defensible output"
#: when rank agreement was withdrawn (see the module docstring). One constant,
#: so the sentence cannot drift between the capability report, the comparison
#: headline and the browser.
RELATIVE_EXPOSURE_CLAIM = (
    "Within a single compound on a single tissue, the RELATIVE exposure between "
    "regions is the defensible output. The absolute level is not, and comparison "
    "BETWEEN compounds is not."
)

#: Why the first half of that sentence stands — stated at the two DIFFERENT
#: strengths the two measurements actually support. They are not interchangeable
#: and must not be written as if they were: one is exact for the finite-window
#: field, the other only for a converged one.
RELATIVE_EXPOSURE_BASIS = (
    "It rests on two measured invariances, and they are NOT equally strong. "
    "(1) PLASMA/DOSE SCALE — EXACT for the finite-window field. A per-compound "
    "plasma factor is an exact linear rescale of the whole AUC field: worst "
    "relative deviation 1.41e-14 over plasma factors 0.005 to 26, a 5,200x span, "
    "measured at the shipping configuration on real Xenium GBM tissue "
    "(results/pk_ranking_robustness.md sec. 1.2). "
    "(2) TISSUE CAPACITY — EXACT only for a CONVERGED field. The omitted "
    "capacity term divides the transport coefficient and the sink together, so "
    "the steady-state field, the penetration length and the interstitial:plasma "
    "ratio are exactly invariant to it, verified to 1e-9 "
    "(results/capacity_ranking_robustness.md sec. 2b). But it enters as a pure "
    "time rescale, so over a FINITE exposure window a between-region ratio is "
    "capacity-invariant only where the field has converged — the same document "
    "measures the finite-window kill ranking moving under it. So a within-map "
    "ratio is structurally blind to a per-compound dose factor always, and to a "
    "per-compound capacity term where the map has converged, and only that far."
)

#: What the claim does NOT cover. Not optional: an invariance to two SCALE
#: factors says nothing about the law that sets the SHAPE, and a ratio is a
#: shape quantity.
RELATIVE_EXPOSURE_LIMITS = (
    "NEITHER INVARIANCE PROTECTS AGAINST ERROR IN THE TRANSPORT LAW ITSELF. A "
    "wrong effective diffusivity changes the SHAPE of the gradient, not just its "
    "level, so a within-map ratio is fully exposed to it — and the descriptor "
    "path feeding that law is UNSCORED. The published R^2 = 0.718 was obtained "
    "by feeding the law the source paper's own measured logP and hydrogen-bond "
    "donor columns; the shipping engine feeds it RDKit's Crippen logP and "
    "Lipinski donor count instead, which miss by RMS 1.62 log units on the 5 of "
    "71 compounds the engine can even be run on, and propagating that error "
    "collapses the score to about 0.27 (results/descriptor_falsification.md "
    "sec. 0). Finally, the perturbations that have been swept against an "
    "ordering or a verdict are A LIST, and a list is all it is: per-compound "
    "IC50, transvascular delivery law, delivery source geometry, tissue "
    "capacity, anchor compound, section x cell-state route, calibration alpha, "
    "ECM tortuosity, KDE bandwidth, endothelial typing threshold, oxygen level, "
    "cell-state route and exposure unit (results/axis_inventory.md sec. 2). "
    "Every axis not named there is UNSWEPT, and silence about an unswept axis is "
    "not evidence about it."
)

#: The whole disclosure, in the order it has to be read.
RELATIVE_EXPOSURE_FULL = " ".join((
    RELATIVE_EXPOSURE_CLAIM,
    RELATIVE_EXPOSURE_BASIS,
    RELATIVE_EXPOSURE_LIMITS,
))


class ProjectError(Exception):
    """A bad request against the project store (mapped to 400/404/409)."""

    def __init__(self, message: str, status: int = 400, detail: Optional[Dict] = None):
        super().__init__(message)
        self.message = message
        self.status = int(status)
        self.detail = detail or {}


# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #
def projects_dir() -> str:
    """Absolute path to the project store (created on demand).

    A HOSTED INSTANCE KEEPS NOTHING WITHOUT THIS.

    Measured 2026-09-20 in a real container: this path and
    ``projects.projects_dir()`` were both hardcoded off ``repo_root()``, which
    puts them at ``/app/app/runs`` and ``/app/app/projects`` -- siblings under
    ``/app/app``, whose parent also holds ``server/`` and ``web/``. Render
    allows ONE persistent disk per service, and mounting it at the only
    directory that contains both would hide the application. So every run and
    every project a hosted user created was destroyed on the next deploy.

    That is survivable for runs, which are re-derivable from a RunSpec. It is
    not survivable for PROJECTS, which are the thing a user is asked to
    organise their work into -- a workspace that forgets is worse than no
    workspace, because it invites the investment and then loses it.

    ``TMESIM_STATE_DIR`` is one directory an operator can point a single disk
    at; ``runs/`` and ``projects/`` become children of it, so one mount covers
    both and neither shadows the code. Unset, the behaviour is byte-identical
    to before -- the local single-user workflow this was built for is
    untouched, and no existing installation moves.

    See ``app/server/cache.py::state_dir``, which this shares.
    """
    from .cache import state_dir

    root = state_dir()
    if root:
        return private_makedirs(os.path.join(root, "projects"))
    return private_makedirs(os.path.join(repo_root(), "app", "projects"))


# --------------------------------------------------------------------------- #
# Capability
# --------------------------------------------------------------------------- #
def capability_report(section_id: str) -> Dict:
    """What this specimen's PANEL can and cannot support, from its gene list.

    Captured when a specimen is added and refreshed on read, because it is cheap
    (``state_route_support`` reads the memoised ``features.tsv.gz`` symbols, never
    the count matrix) and because a stale capability claim is worse than none.

    The load-bearing field is ``absolute_kill_supported``. It is ``False`` unless
    the panel carries at least one hypoxia-response transcript, which is the rule
    ``results/cell_state_impact.md`` §5 arrives at: with no hypoxia gene, the
    'oxygen' default's band is empty for a bookkeeping reason, 'markers' declares
    hypoxia UNMEASURED, and 'vessel_distance' is a geometric proxy — so NO route
    clears the bar and the 0.3252-absolute swing (§3.2) is unresolved.

    What survives when the level does not, verbatim from
    :data:`RELATIVE_EXPOSURE_FULL` (this text replaced "the cross-drug RANKING is
    invariant across all three routes and is what survives", withdrawn 2026-08 —
    see the module docstring): *Within a single compound on a single tissue, the
    RELATIVE exposure between regions is the defensible output. The absolute
    level is not, and comparison BETWEEN compounds is not.* It rests on two
    measured invariances, and they are NOT equally strong. (1) PLASMA/DOSE SCALE
    — EXACT for the finite-window field: a per-compound plasma factor is an exact
    linear rescale of the whole AUC field, worst relative deviation 1.41e-14 over
    plasma factors 0.005 to 26, a 5,200x span, measured at the shipping
    configuration on real Xenium GBM tissue
    (``results/pk_ranking_robustness.md`` sec. 1.2). (2) TISSUE CAPACITY — EXACT
    only for a CONVERGED field: the omitted capacity term divides the transport
    coefficient and the sink together, so the steady-state field, the penetration
    length and the interstitial:plasma ratio are exactly invariant to it, verified
    to 1e-9 (``results/capacity_ranking_robustness.md`` sec. 2b) — but it enters
    as a pure time rescale, so over a FINITE exposure window a between-region
    ratio is capacity-invariant only where the field has converged, and that same
    document measures the finite-window kill ranking moving under it. NEITHER
    INVARIANCE PROTECTS AGAINST ERROR IN THE TRANSPORT LAW ITSELF: a wrong
    effective diffusivity changes the SHAPE of the gradient, not just its level,
    so a within-map ratio is fully exposed to it — and the descriptor path
    feeding that law is UNSCORED, missing the measured columns by RMS 1.62 log
    units on the 5 of 71 compounds the engine can even be run on and collapsing
    the published R^2 = 0.718 to about 0.27
    (``results/descriptor_falsification.md`` sec. 0). The perturbations swept
    against an ordering or a verdict are A LIST — per-compound IC50, transvascular
    delivery law, delivery source geometry, tissue capacity, anchor compound,
    section x cell-state route, calibration alpha, ECM tortuosity, KDE bandwidth,
    endothelial typing threshold, oxygen level, cell-state route and exposure
    unit (``results/axis_inventory.md`` sec. 2) — and every axis not named there
    is UNSWEPT, with silence about an unswept axis being no evidence about it.

    Note what this does NOT claim: a gene being ON the panel says nothing about
    its detection RATE, which is a property of the count matrix and is only
    measured at reconstruction time. That distinction is exactly what disqualifies
    the 319-gene flagship panel's proliferation call (MKI67 + TOP2A are present
    and fire in 1.42% of tumour cells, 13x below the published 19-40% Ki-67 range,
    §1.1). ``markers_route.per_cell_call_knowable`` says so explicitly.
    """
    if section_id not in all_sections():
        raise ProjectError(
            f"unknown section {section_id!r}; known ids: {sorted(all_sections())}",
            status=400,
        )
    section = describe_section(section_id)
    support = state_route_support(section_id)
    source = section.get("source")
    genes_known = bool(support.get("genes_known"))
    hypoxia = list(support.get("hypoxia_present") or [])
    prolif = list(support.get("proliferation_present") or [])
    cycle = list(support.get("cycle_specific_present") or [])

    # --- per-route support, mirroring app.server.api._validate_state_mode ---- #
    if source == "visium":
        markers = {
            "supported": False,
            "reason": (
                "state_mode='markers' needs per-cell expression and is Xenium "
                "only; a Visium section is spot-level and carries no per-cell "
                "transcripts, so reconstruct_real would raise."
            ),
        }
    elif not genes_known:
        markers = {
            "supported": None,
            "reason": (
                "this section's panel gene list could not be read from "
                "data_cache/, so no verdict on the marker route is evidence."
            ),
        }
    elif not prolif and not hypoxia:
        markers = {
            "supported": False,
            "reason": (
                "the panel carries none of the proliferation transcripts "
                f"{support['proliferation_missing']} and none of the "
                f"hypoxia-response transcripts {support['hypoxia_missing']}. "
                "NOTHING would be measured, so the route has no content here."
            ),
        }
    else:
        markers = {
            "supported": True,
            "reason": (
                f"proliferation from {prolif or 'nothing on this panel'}; "
                f"hypoxia from {hypoxia or 'nothing on this panel'}."
            ),
        }
    markers["per_cell_call_knowable"] = False
    markers["per_cell_call_reason"] = (
        "whether the marker route makes a PER-CELL call or falls back to a "
        "neighbourhood enrichment rule depends on the union DETECTION RATE "
        "clearing a 5% floor. That is a property of the count matrix, is not "
        "knowable from the panel, and is measured at reconstruction time and "
        "returned in provenance.cell_state.detail.coverage. Below the floor the "
        "fallback's ORDERING is measured but its LEVEL moves 8.79x with the "
        "window radius (results/cell_state_impact.md §1.1, §2.4)."
    )

    routes = {
        "oxygen": {
            "supported": True,
            "is_measured": False,
            "reason": (
                "the shipping default bands the SIMULATED oxygen field. It always "
                "runs, but nothing in it is measured: OxygenConfig.consumption has "
                "never been fitted to a measured oxygen tension, and the VIABLE "
                "band moves from 0.00% to 10.04% purely with the tissue-mask "
                "convention and the mesh (results/cell_state_impact.md §2.2)."
            ),
        },
        "vessel_distance": {
            "supported": source != "visium",
            "is_measured": False,
            "reason": (
                "bands exp(-d/L) in distance to the nearest MEASURED endothelial "
                "cell. The vessel positions are measured; the state is not."
                if source != "visium" else
                "a Visium section places cells by deconvolution and carries no "
                "measured endothelial positions to take a distance from."
            ),
        },
        "markers": dict(markers, is_measured=bool(markers["supported"])),
    }

    if not genes_known:
        absolute_ok = False
        absolute_reason = (
            f"section {section_id!r}: its panel gene list is not readable from "
            "this checkout, so whether ANY state route can measure hypoxia here "
            "is unknown. An absolute kill is not supported by an unknown panel."
        )
    elif not hypoxia:
        absolute_ok = False
        absolute_reason = (
            f"section {section_id!r} carries none of the hypoxia-response "
            f"transcripts {support['hypoxia_missing']}, so hypoxia is UNMEASURED "
            "here and NO state route clears the bar "
            "(results/cell_state_impact.md §5): the 'oxygen' default's HYPOXIC "
            "band is empty for a bookkeeping reason, 'markers' declares hypoxia "
            "UNMEASURED, and 'vessel_distance' is a geometric proxy. Kill moves "
            f"by up to {STATE_ROUTE_KILL_SWING_ABS} absolute "
            f"({STATE_ROUTE_KILL_SWING_REL}x relative) purely with the choice of "
            "route, with tissue, compound, transport and RNG held identical. "
            + RELATIVE_EXPOSURE_FULL
        )
    else:
        absolute_ok = True
        absolute_reason = (
            f"section {section_id!r} carries the hypoxia-response transcript(s) "
            f"{hypoxia}, so state_mode='markers' can put a MEASURED label on "
            "cells.state here. The HYPOXIC band is still an ENRICHMENT label "
            "whose SIZE is uncalibrated (VEGFA/ADM/SLC2A1 all have real "
            "non-hypoxic expression), so its magnitude must not be quoted as a "
            "hypoxic fraction — but the panel does clear the bar of "
            "results/cell_state_impact.md §5."
        )

    return {
        "section": section_id,
        "label": section.get("label"),
        "available": bool(section.get("available")),
        "unavailable_reason": section.get("reason"),
        "source": source,
        "panel": section.get("panel"),
        "assay": section.get("assay"),
        "accession": section.get("accession"),
        "n_cells": section.get("n_cells"),
        "n_genes": section.get("n_genes"),
        "n_spots": section.get("n_spots"),
        "extent_um": section.get("extent_um"),
        "has_outcome_label": bool(section.get("has_outcome_label")),
        # --- panel content ------------------------------------------------- #
        "genes_known": genes_known,
        "proliferation_present": prolif,
        "proliferation_missing": list(support.get("proliferation_missing") or []),
        "cycle_specific_present": cycle,
        "hypoxia_present": hypoxia,
        "hypoxia_missing": list(support.get("hypoxia_missing") or []),
        # --- verdicts ------------------------------------------------------ #
        "state_routes": routes,
        "supports_cell_state": bool(genes_known and hypoxia and markers["supported"]),
        "absolute_kill_supported": bool(absolute_ok),
        "absolute_kill_reason": absolute_reason,
        "captured_at": time.time(),
    }


# --------------------------------------------------------------------------- #
# Thumbnail — a project's card art, and the one thing it must never be
# --------------------------------------------------------------------------- #
#: What ``{"kind": "none"}`` is for. A thumbnail is the most glanceable thing on
#: the home page and therefore the easiest place in the product to publish a
#: number nobody computed: an empty stacked bar reads as "this tissue is all
#: one type", and zeros read as "measured, and it was zero". Neither is true
#: when no run has finished. So the ONLY two outputs are a composition lifted
#: verbatim out of a finished run's own reconstruction report, or a refusal with
#: a reason — never a placeholder shaped like data.
#:
#: The same rule ``_composition_rows`` already applies to the comparison table:
#: results/cell_state_impact.md §3.1 — "UNMEASURED and assigned to no cell" is
#: NOT the statement "computed and found to be zero" and must not be rendered as
#: if it were.
THUMBNAIL_NEVER_FABRICATES = (
    "A thumbnail is either a cell-type composition read verbatim from a finished "
    "run's reconstruction report, or kind='none' with the reason no composition "
    "exists. It is never an empty chart and never zeros: an empty bar reads as "
    "'all one cell type' and a zero reads as 'measured, and it was zero', and "
    "neither has been computed when no run has finished."
)


def run_composition(run: Dict) -> Optional[Dict]:
    """The thumbnail block of ONE run view, or ``None`` if it cannot supply one.

    ``run`` is a ``RunResult`` view — the dict ``Job.result_view()`` returns and
    the dict ``app/runs/<run_id>/result.json`` stores; they are the same shape,
    which is what lets the store read a thumbnail off disk without importing the
    JobManager.

    Composition lives in ``run["report"]``, written by the reconstruction and
    persisted by ``JobManager._persist``: ``celltype_fractions`` (the marker-gated
    per-type fractions), ``n_cells`` (the segmented cell count) and the section on
    ``run["spec"]``. Nothing here is recomputed and nothing is defaulted — a
    missing field returns ``None`` so the caller emits ``kind='none'`` with a
    reason, per :data:`THUMBNAIL_NEVER_FABRICATES`.

    ``status != "done"`` returns ``None`` unconditionally. A failed or cancelled
    run can still carry a ``report`` (the reconstruction succeeds before the solve
    fails), and showing that tissue on a card would present a run that did not
    finish as one that did.
    """
    if not isinstance(run, dict) or run.get("status") != "done":
        return None
    report = run.get("report")
    if not isinstance(report, dict):
        return None
    fractions = report.get("celltype_fractions")
    if not isinstance(fractions, dict) or not fractions:
        return None
    clean: Dict[str, float] = {}
    for label, value in fractions.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        v = float(value)
        if not math.isfinite(v):
            continue
        clean[str(label)] = v
    if not clean:
        return None
    n_cells = report.get("n_cells")
    if isinstance(n_cells, bool) or not isinstance(n_cells, (int, float)):
        return None
    section = (run.get("spec") or {}).get("section")
    if not isinstance(section, str) or not section:
        section = report.get("source")
    if not isinstance(section, str) or not section:
        return None
    return {
        "kind": "composition",
        "section": section,
        "n_cells": int(n_cells),
        "type_fractions": clean,
        "run_id": run.get("run_id"),
        "measured": (
            "Cell positions are segmented from the spatial assay; the cell-TYPE "
            "call is a marker-gated classification of measured transcripts. This "
            "is a chart of composition, not a picture of the tissue."
        ),
    }


def _no_thumbnail(reason: str) -> Dict:
    return {"kind": "none", "reason": reason}


# --------------------------------------------------------------------------- #
# The store
# --------------------------------------------------------------------------- #
class ProjectStore:
    """Projects on disk under ``app/projects/``, one JSON document each.

    Deliberately not a database: a study is a handful of ids and a name, the
    write rate is a click per minute, and a plain JSON file per project is
    inspectable, diffable and trivially recoverable. Writes go through a temp
    file + ``os.replace`` so a crashed write can never leave a half-parsed
    project behind, and every mutation is serialised by one lock.
    """

    def __init__(self, root: Optional[str] = None,
                 runs_root: Optional[str] = None):
        self.root = root or projects_dir()
        private_makedirs(self.root)
        self._lock = threading.RLock()
        # ``app/runs/`` is READ, never written, and only to build a thumbnail.
        # Left None it is resolved lazily by ``runs_dir()`` -- resolving it here
        # would create app/runs/ as a side effect of importing this module,
        # which is the JobManager's job and not the project store's.
        self.runs_root = runs_root
        # run_id -> ((st_mtime_ns, st_size), (run_stub | None, composition | None))
        # See ``_run_doc_composition`` for what this does and does not go stale
        # on. Values are small by construction: a composition and a two-field
        # stub, never a run body.
        self._composition_cache: Dict[
            str, Tuple[Tuple[int, int], Tuple[Optional[Dict], Optional[Dict]]]] = {}

    # -- paths -------------------------------------------------------------- #
    def _dir(self, project_id: str) -> str:
        # SECURITY: THE ID IS A PATH SEGMENT, SO IT IS VALIDATED HERE AND NOT
        # AT THE CALL SITES.
        #
        # Starlette percent-decodes before route matching, so `%2e%2e` arrives
        # as the literal `..` inside a single segment and the router's `[^/]+`
        # never sees a slash. Joined straight onto the store root it resolved
        # one level up. MEASURED 2026-09-03 with a planted canary:
        #   GET    /api/projects/%2e%2e  -> 200, served the file above the root
        #   DELETE /api/projects/%2e%2e  -> 200, DELETED it
        # A NUL byte in the same position reached open() and raised
        # ValueError("embedded null byte") -- a 500 out of a lookup.
        #
        # Bounded to one level and to the fixed filename, and on a real
        # deployment the reachable object does not exist -- so today's blast
        # radius is nil. It was still a read-and-delete primitive one root
        # reconfiguration away from mattering, and DELETE does not care that
        # the file it removed was uninteresting.
        #
        # `PROJECT_ID_RE` already existed and was applied only when LISTING the store
        # directory. The guard was written and then not wired to the boundary
        # that needed it.
        if not PROJECT_ID_RE.match(project_id or ""):
            raise ProjectError(f"no such project {project_id!r}",
                               status=404)
        return os.path.join(self.root, project_id)

    def _path(self, project_id: str) -> str:
        return os.path.join(self._dir(project_id), "project.json")

    def _run_path(self, run_id: str) -> Optional[str]:
        """Absolute path of a run's ``result.json``, or ``None`` for a bad id.

        SECURITY: ``run_id`` IS A PATH SEGMENT AND IS VALIDATED HERE, for the
        same reason and in the same shape as :meth:`_dir`. ``%2e%2e`` reached
        ``_dir`` as a literal ``..`` inside one segment and was a read-AND-DELETE
        primitive (fixed f9d443d); this method only ever reads, but an id from a
        project document that an operator hand-edited would otherwise choose
        which file the home page parses. A bad id is not an exception here: it
        is a run this store cannot find, which the caller reports as a true
        ``kind='none'`` reason.
        """
        if not RUN_ID_RE.match(run_id or ""):
            return None
        return os.path.join(self.runs_root or runs_dir(), run_id, "result.json")

    # -- io ----------------------------------------------------------------- #
    def _read(self, project_id: str) -> Optional[Dict]:
        try:
            with open(self._path(project_id)) as fh:
                doc = json.load(fh)
        except FileNotFoundError:
            return None
        except ProjectError:
            # Already a well-formed refusal with its own status -- the id
            # validation in _dir raises one. Re-raised UNTOUCHED so a 404 stays
            # a 404: the broad catch below would otherwise relabel it a 500 and
            # log "is on disk but could not be read" about a path that was
            # never opened and a project that does not exist. Closing the
            # traversal without this line turned the hole into a crash.
            raise
        except Exception as exc:  # noqa: BLE001 - a corrupt file is reported
            raise ProjectError(
                redact_exception(
                    exc,
                    f"project {project_id!r} is on disk but could not be read"),
                status=500,
            ) from exc
        if not isinstance(doc, dict):
            raise ProjectError(
                f"project {project_id!r} is not a JSON object", status=500)
        doc.setdefault("id", project_id)
        doc.setdefault("specimens", [])
        # MIGRATION, and the only shape it may take. Every project written
        # before 2026-09-20 has no ``last_opened_at`` key, and the field means
        # "a human looked at it" -- so the ONLY honest value for a document that
        # never recorded one is null. Not ``time.time()`` (that would claim an
        # open that did not happen, and would put every old project at the top
        # of a recency sort forever), and not ``created_at`` (creating a project
        # is not opening it, and the home page renders null as "not opened yet").
        #
        # setdefault and NOT a write: a read must not rewrite a document. Backing
        # this default into every project.json on read would rewrite the whole
        # store on the first home-page load, move every mtime, and -- since
        # _write stamps updated_at -- reorder a workspace sorted by it.
        doc.setdefault("last_opened_at", None)
        return doc

    def _write(self, doc: Dict, touch: bool = True) -> Dict:
        """Persist ``doc``. ``touch=False`` leaves ``updated_at`` exactly as read.

        Exactly one caller passes ``touch=False``: :meth:`mark_opened`. If
        stamping ``last_opened_at`` also moved ``updated_at``, the two fields
        would rise together on every open and the whole distinction the contract
        draws in §3 would be gone — a workspace sorted by ``updated_at`` is
        supposed to be able to say "a background job touched this", and it
        cannot say that if looking at a project counts as touching it.
        """
        if touch:
            doc["updated_at"] = time.time()
        d = self._dir(doc["id"])
        private_makedirs(d)
        tmp = os.path.join(d, "project.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(doc, fh, indent=1, sort_keys=False)
        own_only(tmp, PRIVATE_FILE_MODE)
        os.replace(tmp, os.path.join(d, "project.json"))
        return doc

    def _ids(self) -> List[str]:
        try:
            names = os.listdir(self.root)
        except OSError:
            return []
        return sorted(n for n in names if PROJECT_ID_RE.match(n)
                      and os.path.exists(self._path(n)))

    # -- projects ----------------------------------------------------------- #
    def list(self, thumbnails: bool = True,
             resolve_run: Optional[Callable[[str], Optional[Dict]]] = None,
             ) -> List[Dict]:
        """Every project, newest first, as a summary row.

        ``thumbnails=True`` (the default, because ``GET /api/projects`` is the
        home page's project feed and the contract's ``ProjectSummary`` carries
        the field) adds ``thumbnail`` to each row. See :meth:`thumbnail` for what
        that costs and what it caches; the short version is that it is one
        ``stat`` per attached run plus one JSON parse per run document that has
        changed since it was last read, and that it stops at the first finished
        run of a project rather than scanning them all.

        ``thumbnails=False`` is for callers that want the rows and not the art —
        ``len(store.list())`` in the health route and in ``account.py`` is the
        live example, and it should not pay for 96 stat calls to learn a count.
        """
        out = []
        with self._lock:
            for pid in self._ids():
                try:
                    doc = self._read(pid)
                except ProjectError as exc:
                    out.append({"id": pid, "name": None, "unreadable": exc.message})
                    continue
                if doc is None:
                    continue
                row = self._summary(doc)
                if thumbnails:
                    row["thumbnail"] = self._thumbnail(doc, resolve_run)
                out.append(row)
        out.sort(key=lambda p: p.get("created_at") or 0, reverse=True)
        return out

    @staticmethod
    def _summary(doc: Dict) -> Dict:
        specimens = doc.get("specimens") or []
        return {
            "id": doc["id"],
            "name": doc.get("name"),
            "notes": doc.get("notes"),
            "created_at": doc.get("created_at"),
            "updated_at": doc.get("updated_at"),
            # NEW (contract §3). A SEPARATE field from updated_at: null means
            # "not opened yet" -- never an old date, never now.
            "last_opened_at": doc.get("last_opened_at"),
            "n_specimens": len(specimens),
            "n_runs": sum(len(s.get("run_ids") or []) for s in specimens),
            "sections": [s.get("section") for s in specimens],
        }

    def get(self, project_id: str, refresh_capability: bool = True) -> Dict:
        with self._lock:
            doc = self._read(project_id)
        if doc is None:
            raise ProjectError(f"no such project {project_id!r}", status=404)
        if refresh_capability:
            for spec in doc.get("specimens") or []:
                sid = spec.get("section")
                if sid in all_sections():
                    try:
                        spec["capability"] = capability_report(sid)
                    except ProjectError:
                        pass
        return doc

    @staticmethod
    def _text(value, field: str, limit: int, *, required: bool = False,
              strip: bool = True, required_msg: Optional[str] = None) -> str:
        """One gate for every free-text field this store will persist.

        A ``staticmethod`` on the store rather than a module function because
        the rule it enforces is the store's: these limits exist to bound what
        ``list()`` -- and therefore ``GET /api/workspace`` -- can be made to
        return. See the MAX_* constants for the measurement.

        The refusal names the field, the limit and what was sent, because a
        400 that says only "too long" cannot be acted on by whoever is typing.

        EVERY message keeps the ``"<field>: <what is wrong>"`` shape, which is
        load-bearing: ``app/web/js/home/createProject.js::attachFieldError``
        parses that prefix to decide which input to mark red, and a message
        that loses it is silently demoted to a form-level error.
        """
        text = str(value or "")
        if strip:
            text = text.strip()
        if required and not text:
            raise ProjectError(
                required_msg or f"{field}: required, a non-empty {field}")
        if len(text) > limit:
            raise ProjectError(
                f"{field}: {len(text)} characters, above the {limit}-character "
                f"limit. This field is stored and returned in every workspace "
                f"listing, so its size is bounded on purpose.")
        return text

    def create(self, name: str, notes: str = "") -> Dict:
        name = self._text(name, "name", MAX_NAME_CHARS, required=True,
                          required_msg="name: required, a non-empty project name")
        notes = self._text(notes, "notes", MAX_NOTES_CHARS, strip=False)
        with self._lock:
            now = time.time()
            doc = {
                "id": self._next_id(),
                "name": name,
                "notes": str(notes or ""),
                "created_at": now,
                "updated_at": now,
                # Creating a project is NOT opening it.
                "last_opened_at": None,
                "specimens": [],
            }
            return self._write(doc)

    def _next_id(self) -> str:
        """Next free ``p_NNNNN``. Call under ``self._lock``."""
        n = 0
        for pid in self._ids():
            try:
                n = max(n, int(pid.split("_")[-1]))
            except ValueError:
                pass
        return f"p_{n + 1:05d}"

    def update(self, project_id: str, name: Optional[str] = None,
               notes: Optional[str] = None) -> Dict:
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            if name is not None:
                doc["name"] = self._text(
                    name, "name", MAX_NAME_CHARS, required=True,
                    required_msg="name: must not be empty")
            if notes is not None:
                doc["notes"] = self._text(notes, "notes", MAX_NOTES_CHARS,
                                          strip=False)
            return self._write(doc)

    def delete(self, project_id: str) -> Dict:
        """Delete the project document. Runs are NOT deleted — they are
        associated, not owned, and ``app/runs/`` is their only store."""
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            run_ids = [r for s in doc.get("specimens") or []
                       for r in (s.get("run_ids") or [])]
            path = self._path(project_id)
            try:
                os.remove(path)
                os.rmdir(self._dir(project_id))
            except OSError:
                pass
        return {
            "id": project_id, "deleted": True,
            "runs_retained": sorted(set(run_ids)),
            "reason": (
                "specimens hold run IDS into app/runs/, which remains the single "
                "store; deleting a project detaches runs, it never deletes them."
            ),
        }

    # -- recency ------------------------------------------------------------ #
    def mark_opened(self, project_id: str) -> Dict:
        """Stamp ``last_opened_at``. A HUMAN looked at this project.

        Idempotent in the only sense that matters for a recency stamp: calling it
        twice is not an error and leaves the document in the same shape, with the
        later timestamp.

        **It does not move ``updated_at``**, and that is the whole point of the
        field (contract §3). ``updated_at`` answers "did the project DOCUMENT
        change" — a specimen added, a run attached by a background job, the notes
        edited. ``last_opened_at`` answers "did someone look". A workspace sorted
        by ``updated_at`` puts a project on top because a job attached a run to
        it; sorted by ``last_opened_at`` it is on top because you were working on
        it. Bump both on an open and the first sentence stops being true.
        """
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            doc["last_opened_at"] = time.time()
            self._write(doc, touch=False)
        return {
            "ok": True,
            "id": project_id,
            "last_opened_at": doc["last_opened_at"],
            # Returned so a caller can SEE that it did not move.
            "updated_at": doc.get("updated_at"),
            "reason": (
                "last_opened_at means a human looked at this project; updated_at "
                "means the project document changed. Opening is not a change, so "
                "updated_at is deliberately untouched here."
            ),
        }

    # -- duplication -------------------------------------------------------- #
    def duplicate(self, project_id: str, name: Optional[str] = None) -> Dict:
        """Copy a project's specimens and their run REFERENCES into a new project.

        **Runs are neither re-executed nor copied.** ``app/runs/`` is the single
        store for run results and binary payloads (module docstring, contract
        §2), a specimen holds ``run_id``s into it, and a duplicate holds the same
        ids. Two specimens in two projects legitimately point at the same run —
        that is what makes this operation cheap, and it is also why deleting
        either project still deletes no run.

        What is copied verbatim, and why nothing is restamped: each specimen
        keeps its ``id``, ``label``, ``notes``, ``created_at``, ``run_ids`` and
        its captured ``capability`` block. ``created_at`` on a specimen records
        when that section was bound and its panel inspected; moving it to now
        would claim an inspection that did not happen. The new PROJECT gets a new
        ``created_at``, a null ``last_opened_at`` (nobody has opened it) and
        ``duplicated_from``, so the copy's provenance is on the document rather
        than only in someone's memory.

        The copy is a deep copy: mutating the duplicate's specimens can never
        reach the original's in-memory document, and neither project.json is
        shared.
        """
        with self._lock:
            src = self.get(project_id, refresh_capability=False)
            if name is None:
                # TRUNCATE rather than refuse: the caller supplied no name, so
                # a 400 here would be a refusal they cannot act on -- and
                # duplicating a project whose name is already AT the cap would
                # otherwise be impossible, because " (copy)" pushes it over.
                suffix = " (copy)"
                stem = str(src.get("name") or project_id)
                room = MAX_NAME_CHARS - len(suffix)
                new_name = f"{stem[:room]}{suffix}"
            else:
                new_name = self._text(
                    name, "name", MAX_NAME_CHARS, required=True,
                    required_msg="name: must not be empty")
            now = time.time()
            doc = {
                "id": self._next_id(),
                "name": new_name,
                "notes": str(src.get("notes") or ""),
                "created_at": now,
                "updated_at": now,
                "last_opened_at": None,
                "specimens": copy.deepcopy(src.get("specimens") or []),
                "duplicated_from": project_id,
                "duplicated_at": now,
                "duplicated_reason": (
                    "specimens and their run IDS were copied. No run was "
                    "re-executed and no run was copied: app/runs/ is the single "
                    "store and both projects now reference the same run ids."
                ),
            }
            return self._write(doc)

    # -- thumbnail ---------------------------------------------------------- #
    def thumbnail(self, project_id: str,
                  resolve_run: Optional[Callable[[str], Optional[Dict]]] = None,
                  ) -> Dict:
        """The card art for one project (contract §4), or a reasoned refusal.

        ``{"kind": "composition", "section", "n_cells", "type_fractions"}`` from
        a run that actually finished, or ``{"kind": "none", "reason"}``. Never
        anything else — see :data:`THUMBNAIL_NEVER_FABRICATES`.

        ``resolve_run`` maps a ``run_id`` to its ``RunResult`` view, exactly as
        :func:`compare_project` takes one. Passed in (the API has a live
        JobManager) the manager is the truth and nothing is cached here. Omitted,
        the store reads ``app/runs/<run_id>/result.json`` itself and memoises per
        run document — see :meth:`_run_doc_composition`.
        """
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            return self._thumbnail(doc, resolve_run)

    def _thumbnail(self, doc: Dict,
                   resolve_run: Optional[Callable[[str], Optional[Dict]]],
                   ) -> Dict:
        """``_thumbnail`` on an already-read document. Call under ``self._lock``.

        COST, because this runs on every home-page load and a home page that
        takes a second to paint is a worse product than one without card art:
        the scan STOPS at the first run that can supply a composition, so the
        common project pays one stat and (once per changed run document) one JSON
        parse. It is O(runs of this project) only in the case where NO run can
        supply one, which is exactly the case that returns ``kind='none'`` and so
        cannot be cached against the project document — a run finishing does not
        touch project.json, so a thumbnail cached on the project's mtime would
        stay ``none`` after the run that should have filled it finished.

        The FIRST such run, not the newest: a thumbnail is an identity, and one
        that changes every time a background job finishes is a card that flickers
        for no information. Picking the newest would also force a parse of every
        run document to find out which one that is.
        """
        specimens = doc.get("specimens") or []
        if not specimens:
            return _no_thumbnail(
                "this project has no specimens yet, so no tissue has been "
                "reconstructed. Add a specimen and run it.")
        run_ids = [rid for s in specimens for rid in (s.get("run_ids") or [])]
        if not run_ids:
            return _no_thumbnail(
                f"no run is attached to any of this project's {len(specimens)} "
                "specimen(s), so nothing has been reconstructed yet.")

        missing = 0
        statuses: Dict[str, int] = {}
        done_without_composition: List[str] = []
        for rid in run_ids:
            if resolve_run is not None:
                run = resolve_run(rid)
                comp = run_composition(run) if run is not None else None
            else:
                run, comp = self._run_doc_composition(rid)
            if comp is not None:
                return comp
            if run is None:
                missing += 1
                continue
            status = str(run.get("status") or "unknown")
            statuses[status] = statuses.get(status, 0) + 1
            if status == "done":
                done_without_composition.append(rid)

        n = len(run_ids)
        parts = [f"none of the {n} run(s) attached to this project can supply a "
                 "measured composition"]
        if statuses:
            parts.append("attached run status: " + ", ".join(
                f"{k} x{v}" for k, v in sorted(statuses.items())))
        if done_without_composition:
            parts.append(
                f"run(s) {sorted(done_without_composition)} finished but their "
                "report carries no celltype_fractions/n_cells, so there is no "
                "composition to draw")
        if missing:
            parts.append(
                f"{missing} attached run id(s) are not in app/runs/ — a run was "
                "removed or its id is malformed; the project still references it")
        parts.append("No chart is drawn from this: " + THUMBNAIL_NEVER_FABRICATES)
        return _no_thumbnail(". ".join(parts) + ".")

    def _run_doc_composition(self, run_id: str,
                             ) -> Tuple[Optional[Dict], Optional[Dict]]:
        """``(run_view, composition)`` read from ``app/runs/<id>/result.json``.

        Memoised on ``(st_mtime_ns, st_size)`` of that file, which is what makes
        the home page affordable: a run document is written ONCE, by
        ``JobManager._persist`` via ``os.replace``, so a cache hit is a single
        ``stat`` and a miss is the only time 20 kB of JSON is parsed.

        WHAT GOES STALE. Nothing that this store can observe: the key moves the
        moment the file does, a deleted run drops out of the cache through the
        ``OSError`` path, and a run id that was never on disk is not cached at
        all. The one theoretical stale case is a rewrite that lands on the SAME
        nanosecond mtime AND the same byte length, which ``os.replace`` does not
        produce. What is NOT cached here is the choice of run — see
        :meth:`_thumbnail`: the first run that can supply a composition wins and
        keeps winning, by design, so a later run finishing does not repaint the
        card. That is a stable identity, not a stale read.

        A STUB of the run — ``{"run_id", "status"}`` and nothing else — is
        returned alongside, so a refusal can state the true reason (``queued``,
        ``error``, finished-but-no-report) rather than a bare "no thumbnail".
        It is a stub deliberately: the cache must be able to answer the
        no-composition case WITHOUT reparsing, because that is the case that
        scans every run of the project, and holding whole 20 kB run bodies in a
        cache that exists to make a page load would be the wrong trade.
        """
        path = self._run_path(run_id)
        if path is None:
            self._composition_cache.pop(run_id, None)
            return None, None
        try:
            st = os.stat(path)
        except (OSError, ValueError):
            # ValueError: a NUL in the id reaches the syscall as "embedded null
            # byte" and would be a 500 out of a home-page load. RUN_ID_RE
            # already refuses it; this is the second layer, because the
            # traversal fixed in f9d443d was a guard that existed and was not
            # wired to the boundary that needed it.
            self._composition_cache.pop(run_id, None)
            return None, None
        key = (st.st_mtime_ns, st.st_size)
        cached = self._composition_cache.get(run_id)
        if cached is not None and cached[0] == key:
            stub, comp = cached[1]
            return stub, comp
        run = self._read_run_doc(path)
        if run is None:
            stub, comp = None, None
        else:
            stub = {"run_id": run.get("run_id") or run_id,
                    "status": run.get("status")}
            comp = run_composition(run)
        self._composition_cache[run_id] = (key, (stub, comp))
        return stub, comp

    @staticmethod
    def _read_run_doc(path: str) -> Optional[Dict]:
        """Parse one ``result.json``. An unreadable run is a run we cannot draw,
        not a 500: the home page must still load."""
        try:
            with open(path) as fh:
                doc = json.load(fh)
        except (OSError, ValueError):
            return None
        return doc if isinstance(doc, dict) else None

    # -- specimens ---------------------------------------------------------- #
    @staticmethod
    def _find_specimen(doc: Dict, specimen_id: str) -> Dict:
        for spec in doc.get("specimens") or []:
            if spec.get("id") == specimen_id:
                return spec
        raise ProjectError(
            f"no such specimen {specimen_id!r} in project {doc['id']!r}; "
            f"known: {[s.get('id') for s in doc.get('specimens') or []]}",
            status=404,
        )

    def add_specimen(self, project_id: str, section: str,
                     label: Optional[str] = None, notes: str = "") -> Dict:
        if section not in all_sections():
            raise ProjectError(
                f"section: unknown id {section!r}. Known ids: {sorted(all_sections())}")
        # BEFORE capability_report, which reads the section off disk: refusing
        # oversize text must not cost a disk read first.
        if label is not None:
            label = self._text(label, "label", MAX_NAME_CHARS)
        notes = self._text(notes, "notes", MAX_NOTES_CHARS, strip=False)
        cap = capability_report(section)
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            specimens = doc.setdefault("specimens", [])
            if any(s.get("section") == section for s in specimens):
                raise ProjectError(
                    f"section {section!r} is already a specimen of project "
                    f"{project_id!r}. One specimen per section keeps a matched "
                    "comparison unambiguous.",
                    status=409,
                )
            n = 0
            for s in specimens:
                try:
                    n = max(n, int(str(s.get("id", "s_00000")).split("_")[-1]))
                except ValueError:
                    pass
            specimen = {
                "id": f"s_{n + 1:05d}",
                "label": (label or cap.get("label") or section).strip(),
                "section": section,
                "notes": str(notes or ""),
                "created_at": time.time(),
                "run_ids": [],
                "capability": cap,
            }
            specimens.append(specimen)
            self._write(doc)
        return specimen

    def update_specimen(self, project_id: str, specimen_id: str,
                        label: Optional[str] = None,
                        notes: Optional[str] = None) -> Dict:
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            spec = self._find_specimen(doc, specimen_id)
            if label is not None:
                spec["label"] = self._text(
                    label, "label", MAX_NAME_CHARS, required=True,
                    required_msg="label: must not be empty")
            if notes is not None:
                spec["notes"] = self._text(notes, "notes", MAX_NOTES_CHARS,
                                           strip=False)
            self._write(doc)
        return spec

    def remove_specimen(self, project_id: str, specimen_id: str) -> Dict:
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            spec = self._find_specimen(doc, specimen_id)
            doc["specimens"] = [s for s in doc["specimens"]
                                if s.get("id") != specimen_id]
            self._write(doc)
        return {
            "project_id": project_id, "specimen_id": specimen_id, "removed": True,
            "runs_retained": list(spec.get("run_ids") or []),
            "reason": "runs are associated, not owned; app/runs/ is untouched.",
        }

    # -- run association ---------------------------------------------------- #
    def attach_run(self, project_id: str, specimen_id: str, run_id: str,
                   run: Optional[Dict]) -> Dict:
        """Associate an existing run with a specimen.

        ``run`` is the ``RunResult`` view from the JobManager, or ``None`` if no
        such run exists. The section on the run's own spec MUST equal the
        specimen's section: a specimen is a patient, and a run from another
        patient filed under it would silently corrupt every comparison built on
        top of this store. That is refused, not warned about.
        """
        if run is None:
            raise ProjectError(
                f"no such run {run_id!r}; nothing was attached", status=404)
        run_section = ((run.get("spec") or {}).get("section"))
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            spec = self._find_specimen(doc, specimen_id)
            if run_section != spec.get("section"):
                raise ProjectError(
                    f"run {run_id!r} was computed on section {run_section!r}, but "
                    f"specimen {specimen_id!r} is bound to section "
                    f"{spec.get('section')!r}. A specimen is a patient; filing "
                    "another patient's run under it would corrupt every "
                    "cross-specimen comparison built on this project.",
                    status=409,
                )
            run_ids = spec.setdefault("run_ids", [])
            already = run_id in run_ids
            if not already:
                run_ids.append(run_id)
                self._write(doc)
        return {
            "project_id": project_id, "specimen_id": specimen_id,
            "run_id": run_id, "attached": not already,
            "n_runs": len(spec.get("run_ids") or []),
            "reason": ("run was already attached to this specimen" if already
                       else None),
        }

    def detach_run(self, project_id: str, specimen_id: str,
                   run_id: str) -> Dict:
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            spec = self._find_specimen(doc, specimen_id)
            run_ids = spec.get("run_ids") or []
            if run_id not in run_ids:
                raise ProjectError(
                    f"run {run_id!r} is not attached to specimen "
                    f"{specimen_id!r}", status=404)
            spec["run_ids"] = [r for r in run_ids if r != run_id]
            self._write(doc)
        return {
            "project_id": project_id, "specimen_id": specimen_id,
            "run_id": run_id, "detached": True,
            "reason": ("the run itself is untouched in app/runs/; only the "
                       "association was removed"),
        }


#: Process-wide store used by the API.
PROJECT_STORE = ProjectStore()


# --------------------------------------------------------------------------- #
# Matching
# --------------------------------------------------------------------------- #
#: Every RunSpec field EXCEPT ``section`` and ``drug``. Two runs match when this
#: signature is equal: same mesh, same modes, same PK, same window — so the only
#: thing that differs between the two sides of a comparison is the SPECIMEN.
#: ``drug`` sits outside the signature and is carried alongside it as the row
#: key (see ``spec_signature`` and the ``{signature: {drug: {specimen: run}}}``
#: grouping in ``compare_project``), which is what makes a comparison
#: same-compound by construction.
_SPEC_MATCH_FIELDS: Sequence[Tuple[str, object]] = (
    ("dx_um", None),
    ("sensitivity_mode", "uniform"),
    ("vascular_mode", "discrete"),
    ("oxygen_mode", "pde"),
    ("state_mode", "oxygen"),
    ("mgmt", None),
    ("total_time_s", None),
    # The COUNTERFACTUAL arm (PDConfig resistance multipliers forced to 1). A
    # true arm is not a prediction for the tissue, so it must never be matched
    # against a normal arm of the same drug.
    ("resistance_off", False),
)
_PK_MATCH_FIELDS = ("plasma_peak", "half_life_s", "n_doses", "dose_interval_s")

#: The only RunSpec keys a signature deliberately drops. Everything else — named
#: above or not — is folded in, so a field added to the RunSpec later can never
#: silently let two different experiments be matched. ``resistance_off`` is the
#: live example: it is the COUNTERFACTUAL arm (resistance multipliers set to 1),
#: and matching it against a normal arm would difference a prediction against a
#: what-if.
#:
#: ``drug`` is in here for a reason that CHANGED when cross-compound ordering was
#: withdrawn, and the hashing did not change with it — see ``spec_signature``.
_SPEC_IGNORED = frozenset({"section", "drug", "pk"})


def _round(v):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return v
    return round(float(v), 6)


def spec_signature(spec: Dict) -> str:
    """Canonical string for everything in a RunSpec except section and drug.

    **What it is for now.** It identifies the EXPERIMENT — the mesh, the modes,
    the PK, the exposure window, the arm — and nothing about the tissue or the
    compound. That is what a comparison needs to be legitimate: anything else
    differing (a different mesh, a different state route, a different exposure
    window, a counterfactual arm) makes two runs incomparable, and grouping by
    this signature is what stops them being compared anyway.

    **Why ``drug`` is still excluded, restated.** It used to be excluded so runs
    differing only in patient+drug could be matched into one cross-compound
    ranking. That ranking is gone (see the module docstring), so this is no
    longer a matching requirement: ``compare_project`` groups
    ``{signature: {drug: {specimen: run}}}`` and every comparison it emits is
    within one ``drug`` key. Excluding ``drug`` from the signature is what lets
    the compound be that row key rather than being smeared into the match key —
    one signature names one experiment, and the whole compound set run under it
    hangs off it. **The hashing is unchanged and must stay unchanged**: this
    string is the stable identity of a matched set, it is returned to clients in
    ``matched_sets[].signature``, and moving ``drug`` into it would silently
    resignature every set in every existing study.

    Two runs with the same signature therefore differ only in *which specimen*
    and *which compound*, and the compound is pinned equal before anything is
    differenced.

    Fails CLOSED on an unrecognised field: any RunSpec key that is not
    ``section``, ``drug`` or ``pk`` enters the signature even if this module has
    never heard of it, so a field added to the RunSpec later splits a matched set
    — which surfaces as "not evaluable" and prompts a rerun — rather than
    silently pairing two runs that meant different things. An unknown key whose
    value is ``None`` or ``False`` is treated as absent, because "flag off" and
    "flag not sent" are the same experiment and must not split a set.
    """
    spec = spec or {}
    known = {k for k, _ in _SPEC_MATCH_FIELDS}
    parts = [f"{k}={_round(spec.get(k, default))}"
             for k, default in _SPEC_MATCH_FIELDS]
    for k in sorted(set(spec) - known - _SPEC_IGNORED):
        v = spec.get(k)
        if v is None or v is False:
            continue
        parts.append(f"{k}={_round(v)}")
    pk = spec.get("pk") or {}
    parts += [f"pk.{k}={_round(pk.get(k))}" for k in _PK_MATCH_FIELDS]
    for k in sorted(set(pk) - set(_PK_MATCH_FIELDS)):
        v = pk.get(k)
        if v is None or v is False:
            continue
        parts.append(f"pk.{k}={_round(v)}")
    return "; ".join(parts)


def _drug_of(run: Dict) -> Optional[str]:
    d = (run.get("spec") or {}).get("drug")
    return str(d).strip().lower() if isinstance(d, str) and d.strip() else None


# --------------------------------------------------------------------------- #
# Comparability of an ABSOLUTE number across specimens
# --------------------------------------------------------------------------- #
def _run_state_measured(run: Dict) -> bool:
    prov = run.get("provenance") or {}
    if "state_is_measured" in prov:
        return bool(prov.get("state_is_measured"))
    cs = (run.get("report") or {}).get("cell_state") or {}
    return bool(cs.get("is_measured"))


def _kill_comparability(entries: List[Dict]) -> Dict:
    """Is an absolute difference between these specimens readable at all?

    ``entries`` are ``{specimen, capability, run}`` for ONE compound — the caller
    builds them from a single ``drug`` key, so every entry here is the same
    compound on a different specimen and this function never sees two compounds
    at once. That is the only comparison the module makes: cross-compound
    ordering was withdrawn (module docstring), cross-specimen comparison of one
    compound was not.

    The verdict is ``False`` unless EVERY specimen's panel can support a measured
    cell state AND every run actually used a measured route — the rule from
    ``results/cell_state_impact.md`` §5, and the reason
    ``results/heldout_patient.md`` §3.5 gives for refusing to compare absolute
    kills between its two sections.
    """
    drugs = {_drug_of(e["run"]) for e in entries}
    if len(drugs) > 1:  # pragma: no cover - the caller groups by drug
        raise ProjectError(
            f"_kill_comparability was handed more than one compound {sorted(drugs)}. "
            "Only the SAME compound may be compared across specimens; ordering "
            "two different compounds by predicted effect is withdrawn "
            "(results/pk_ranking_robustness.md, "
            "results/capacity_ranking_robustness.md).",
            status=500,
        )
    blockers: List[str] = []
    caveats: List[str] = []

    unsupported = [e for e in entries
                   if not (e["capability"] or {}).get("absolute_kill_supported")]
    for e in unsupported:
        blockers.append(
            f"{e['specimen']['label']} ({e['specimen']['section']}): "
            + str((e["capability"] or {}).get("absolute_kill_reason"))
        )

    measured = {e["specimen"]["id"]: _run_state_measured(e["run"])
                for e in entries}
    if len(set(measured.values())) > 1:
        yes = [e["specimen"]["label"] for e in entries
               if measured[e["specimen"]["id"]]]
        no = [e["specimen"]["label"] for e in entries
              if not measured[e["specimen"]["id"]]]
        blockers.append(
            f"cells.state is MEASURED on {yes} and ASSUMED on {no}. Differencing "
            "a measured kill against an assumed one is not a patient difference; "
            "the state label alone moves kill by up to "
            f"{STATE_ROUTE_KILL_SWING_ABS} absolute "
            f"({STATE_ROUTE_KILL_SWING_REL}x relative) with everything else held "
            "identical (results/cell_state_impact.md §3.2)."
        )
    elif not any(measured.values()):
        blockers.append(
            "no specimen used a MEASURED cell-state route: every run here is on "
            "state_mode='oxygen' or 'vessel_distance', neither of which measures "
            "anything. Select state_mode='markers' on both specimens — where the "
            "panel supports it — before reading an absolute difference."
        )

    # Caveats that apply even when the panels DO clear the bar.
    panels = {e["specimen"]["section"]: (e["capability"] or {}).get("panel")
              for e in entries}
    if len(set(panels.values())) > 1:
        caveats.append(
            f"the specimens are on DIFFERENT panels ({panels}). Every absolute "
            "field difference is jointly a difference in biology and a difference "
            "in what is measurable: the ECM field is a max-normalised KDE over "
            "whichever matrix transcripts the panel carries, so an 11-marker and "
            "a 3-marker field have no common absolute scale even on identical "
            "tissue (results/heldout_patient.md §3.5)."
        )
    modes = {(e["run"].get("spec") or {}).get("sensitivity_mode", "uniform")
             for e in entries}
    if modes == {"uniform"}:
        caveats.append(
            "sensitivity_mode='uniform' on every side: each kill is a "
            "no-intrinsic-resistance UPPER BOUND, so the delta is a delta "
            "between two upper bounds, not a measured treatment effect."
        )
    caveats.append(
        "the absolute oxygen level is uncalibrated on every section "
        "(OxygenConfig.consumption has never been fitted to a measured oxygen "
        "tension), and it sets the band occupancies every absolute kill inherits "
        "(results/heldout_patient.md §3.5)."
    )
    caveats.append(
        "neither section carries a drug-response outcome label, so no kill here "
        "has ever been scored against a measured response. NOT FOR CLINICAL USE."
    )
    caveats.append(RELATIVE_EXPOSURE_FULL)

    return {
        "comparable": not blockers,
        "compound": sorted(drugs)[0] if drugs else None,
        "blockers": blockers,
        "caveats": caveats,
        "state_measured": measured,
    }


# --------------------------------------------------------------------------- #
# Composition
# --------------------------------------------------------------------------- #
def _tissue_area_mm2(report: Dict) -> Optional[float]:
    ext = report.get("extent_um")
    cov = report.get("tissue_coverage_fraction")
    if not ext or len(ext) < 2 or cov is None:
        return None
    try:
        return float(ext[0]) * float(ext[1]) * float(cov) / 1e6
    except (TypeError, ValueError):
        return None


def _composition_rows(entries: List[Dict]) -> List[Dict]:
    """Cross-specimen composition, each row with provenance and a verdict.

    This is where a user sees that patient B is more vascularised than patient A
    — and immediately sees which half of that sentence the data supports.
    ``results/heldout_patient.md`` §1.1 is the cautionary case: the endothelial
    FRACTION fell (9.74% -> 4.18%) while the endothelial AREAL DENSITY rose
    (143 -> 382 /mm2), because the second section is 6.2x more cellular. Both
    rows are shown, adjacent, for exactly that reason.
    """
    rows: List[Dict] = []

    def add(key: str, label: str, unit: str, getter: Callable[[Dict], object],
            kind: str, why: str, comparable: bool, reason: str,
            fmt: str = "num") -> None:
        values, provenance = {}, {}
        for e in entries:
            sid = e["specimen"]["id"]
            try:
                values[sid] = getter(e)
            except Exception:  # noqa: BLE001 - a missing field is a null, not a 500
                values[sid] = None
            provenance[sid] = {"kind": kind, "why": why}
        rows.append({
            "key": key, "label": label, "unit": unit, "fmt": fmt,
            "values": values, "provenance": provenance,
            "comparable": bool(comparable), "reason": reason,
        })

    frac = lambda name: (lambda e: (e["report"].get("celltype_fractions") or {}).get(name))  # noqa: E731
    _MEASURED_POSITIONS = (
        "Cell positions come from the segmented spatial assay; the cell-TYPE "
        "call is a marker-gated classification of measured transcripts, so it is "
        "computed from measured inputs, not read off an instrument."
    )
    _PANEL_WARN = (
        "the two panels use DIFFERENT marker sets for this gate (e.g. ENDOTHELIAL "
        "is PECAM1 on the 319-gene brain panel and the mean z of VWF/FLT1/PLVAP on "
        "the 480-gene IO panel), and a 3-marker gate over sparse probes behaves "
        "differently from a 1-marker gate. results/heldout_patient.md §1.2 "
        "concludes the endothelial cell-type FRACTION is not a portable quantity "
        "across panels."
    )
    panels = {(e["capability"] or {}).get("panel") for e in entries}
    same_panel = len(panels) == 1

    add("n_cells", "Cells measured", "n", lambda e: e["report"].get("n_cells"),
        "measured", "Segmented cell count from the assay.", True,
        "A cell count is a direct measurement on both sides.", fmt="int")
    add("tissue_area_mm2", "Tissue area", "mm2",
        lambda e: _tissue_area_mm2(e["report"]),
        "computed",
        "extent_um[0] x extent_um[1] x tissue_coverage_fraction, from the "
        "reconstruction's own tissue mask.",
        True, "Derived from the measured extent and the tissue mask.")
    add("cell_density_mm2", "Cell density", "/mm2",
        lambda e: (e["report"]["n_cells"] / _tissue_area_mm2(e["report"])
                   if _tissue_area_mm2(e["report"]) else None),
        "computed", "Measured cell count over the masked tissue area.", True,
        "Both terms are measured or masked-measured; density is comparable. It "
        "is also the denominator that makes every cell-type FRACTION below "
        "section-specific (results/heldout_patient.md §1.1).")

    for name, lab in (("TUMOR", "TUMOR"), ("STROMAL", "STROMAL"),
                      ("IMMUNE", "IMMUNE"), ("ENDOTHELIAL", "ENDOTHELIAL")):
        add(f"frac_{name.lower()}", f"{lab} fraction", "frac", frac(name),
            "computed", _MEASURED_POSITIONS,
            same_panel,
            ("Both specimens are on the same panel, so the same marker gate "
             "produced both numbers." if same_panel else _PANEL_WARN),
            fmt="pct")

    add("endothelial_per_mm2", "ENDOTHELIAL areal density", "/mm2",
        lambda e: (((e["report"].get("celltype_counts") or {}).get("ENDOTHELIAL"))
                   / _tissue_area_mm2(e["report"])
                   if _tissue_area_mm2(e["report"]) else None),
        "computed",
        "Measured endothelial cell count over the masked tissue area. This — not "
        "the fraction — is what the oxygen and drug SOURCE term is built from.",
        same_panel,
        ("Same panel on both sides." if same_panel else
         "The gate differs across panels (see the fraction rows), but note that "
         # 0.67x / 4.19x, NOT the 0.43x / 2.67x this used to cite. Those are
         # the superseded mean_z reading (results/endothelial_density.md:215);
         # heldout_patient.md §1.1's own table, under the matched_null_z gate
         # that ships, reads 10.50% -> 7.07% (0.67x) and 154 -> 645 cells/mm2
         # (4.19x). Same qualitative finding, wrong numbers against the
         # section being cited -- and the correct pair is STRONGER, since the
         # divergence is wider than the retired one claimed.
         "in results/heldout_patient.md §1.1 the fraction FELL 0.67x while this "
         "density ROSE 4.19x on the same pair of sections. They point in "
         "opposite directions; quoting only one of them is how that report's "
         "conclusion gets inverted."))
    add("mean_vascular_density", "Mean vascular density", "frac",
        lambda e: e["report"].get("mean_vascular_density"),
        "computed",
        "Field built by painting discrete segmented vessel objects (or a KDE "
        "under vascular_mode='kde') from the MEASURED endothelial positions.",
        same_panel,
        ("Same panel and the same vessel source on both sides." if same_panel else
         "The vessel painter uses a fixed 7.5 um radius, which at dx = 12 um is "
         "already near the single-cell equivalent radius, so on a denser section "
         "the painted area inflates. results/heldout_patient.md §2.1 measured "
         "exactly this: the vascular area fraction went 2.50% -> 5.79% against a "
         "measured 2.1 +/- 0.3%, and calls it the clearest single failure of "
         "generalisation in the report."))
    add("mean_ecm_density", "Mean ECM density", "frac",
        lambda e: e["report"].get("mean_ecm_density"),
        "computed",
        "Max-normalised KDE over whichever matrix transcripts the panel carries.",
        False,
        "NOT COMPARABLE ACROSS PANELS AS BIOLOGY. The field is max-normalised "
        "over the matrix markers that happen to be on each panel — 11 genes on "
        "the 319-gene brain panel against 3 on the IO panel, with both fibrillar "
        "collagen markers absent from the latter — so the two have no common "
        "absolute scale even on identical tissue (results/heldout_patient.md "
        "§3.5). The ECM markers actually used are listed per specimen.")
    add("ecm_markers_used", "ECM markers used", "genes",
        lambda e: e["report"].get("ecm_markers_used"),
        "measured", "The panel transcripts the ECM field was actually built from.",
        True, "Listed so the row above can be read.", fmt="list")
    add("mean_oxygen", "Mean oxygen (tissue)", "frac",
        lambda e: e["report"].get("mean_oxygen"),
        "assumed",
        "Steady-state reaction-diffusion solve whose consumption constant has "
        "never been fitted to a measured oxygen tension.",
        False,
        "NOT COMPARABLE AS BIOLOGY. The absolute level is uncalibrated on every "
        "section, and it moves with cell density because a denser tissue is a "
        "bigger sink — results/heldout_patient.md §3.5 calls the 0.2955 vs "
        "0.1232 gap a statement about OxygenConfig.consumption and cell density, "
        "not about the two tumours.")
    add("tissue_coverage_fraction", "Tissue coverage of bbox", "frac",
        lambda e: e["report"].get("tissue_coverage_fraction"),
        "computed", "Fraction of the bounding box inside the tissue mask.",
        True, "A geometric property of the mask on both sides.", fmt="pct")

    # --- measured cell state, where the panel supports it ------------------- #
    for band in ("VIABLE", "QUIESCENT", "HYPOXIC"):
        values, provenance = {}, {}
        comparable = True
        reasons = []
        for e in entries:
            sid = e["specimen"]["id"]
            cs = (e["report"].get("cell_state") or {})
            fr = cs.get("fractions_tumor") or {}
            values[sid] = fr.get(band)
            is_measured = bool(cs.get("is_measured"))
            unmeasured = list(cs.get("unmeasured_states") or [])
            if band in unmeasured:
                values[sid] = None
                provenance[sid] = {
                    "kind": "unknown",
                    "why": (
                        f"{band} is UNMEASURED on this panel and is assigned to "
                        "NO cell. That is NOT the same statement as 'computed "
                        "and found to be zero' and must not be reported as if it "
                        "were (results/cell_state_impact.md §3.1)."
                    ),
                }
                comparable = False
                reasons.append(
                    f"{e['specimen']['label']}: {band} is UNMEASURED on this "
                    "panel — the band is empty because nothing measured it.")
            elif is_measured:
                provenance[sid] = {
                    "kind": "measured",
                    "why": (
                        f"state_mode={cs.get('mode')!r}: {cs.get('source')}. The "
                        "HYPOXIC band remains an ENRICHMENT label whose SIZE is "
                        "uncalibrated."
                    ),
                }
            else:
                provenance[sid] = {
                    "kind": "assumed",
                    "why": (
                        f"state_mode={cs.get('mode')!r} is not a measurement. "
                        "Across the three routes the tumour VIABLE fraction on "
                        "the flagship section runs 0.00% / 10.04% / 29.80% / "
                        "91.08% depending on route, mask and mesh "
                        "(results/cell_state_impact.md §2.1, §3.1)."
                    ),
                }
                comparable = False
                reasons.append(
                    f"{e['specimen']['label']}: cells.state came from "
                    f"state_mode={cs.get('mode')!r}, which measures nothing.")
        rows.append({
            "key": f"state_{band.lower()}",
            "label": f"Tumour {band} fraction",
            "unit": "frac", "fmt": "pct",
            "values": values, "provenance": provenance,
            "comparable": comparable,
            "reason": ("; ".join(reasons) if reasons else
                       "Both specimens put a MEASURED label on cells.state, so "
                       "the band occupancies are comparable — subject to the "
                       "HYPOXIC band's uncalibrated size."),
        })

    # --- deltas against the first specimen ---------------------------------- #
    ref = entries[0]["specimen"]["id"] if entries else None
    for row in rows:
        row["reference"] = ref
        deltas, ratios = {}, {}
        base = row["values"].get(ref)
        for sid, v in row["values"].items():
            if sid == ref or not isinstance(v, (int, float)) or isinstance(v, bool):
                continue
            if not isinstance(base, (int, float)) or isinstance(base, bool):
                continue
            deltas[sid] = float(v) - float(base)
            ratios[sid] = (float(v) / float(base)) if base else None
        row["delta"] = deltas
        row["ratio"] = ratios
    return rows


# --------------------------------------------------------------------------- #
# The comparison
# --------------------------------------------------------------------------- #
#: The metrics a cross-specimen row may be read on. Every row carries all four;
#: ``metric`` picks which one is the FEATURED column. It used to pick which
#: metric a cross-compound ranking was taken over — there is no ranking now, so
#: it selects the metric the SAME-COMPOUND cross-specimen rows report, and
#: nothing else. It cannot change which rows exist or whether they are readable.
_HEADLINE_METRICS = (
    "kill_fraction", "mean_tumor_auc", "therapeutic_coverage_1x",
    "penetration_length_um",
)


def _headline_block(spec_views: List[Dict], kill_rows: List[Dict],
                    composition_rows: List[Dict], any_incomparable: bool,
                    metric: str) -> Dict:
    """What this comparison actually established — never an ordering.

    Returns ``{tag, level, text}``. The old headline was ``_rank_block``'s
    ``RANKING AGREES / DISAGREES`` verdict; that block is gone (module
    docstring), so this states the two things a cross-specimen comparison still
    supports: which compounds were matched SAME-COMPOUND across which specimens,
    and how much of the tissue contrast is flagged readable.
    """
    n_matched = len(kill_rows)
    compared = sorted({sid for row in kill_rows for sid in (row.get("runs") or {})})
    n_comp = len(composition_rows)
    n_readable = sum(1 for r in composition_rows if r.get("comparable"))
    with_runs = [s for s in spec_views if s.get("n_runs_done")]

    if not kill_rows:
        return {
            "tag": "NOT EVALUABLE",
            "level": "info",
            "text": (
                "SAME-COMPOUND COMPARISON NOT EVALUABLE: no compound is matched "
                "across two specimens at an identical RunSpec. This project has "
                f"{len(spec_views)} specimen(s), {len(with_runs)} of them with a "
                "completed run. Two runs are matched only when everything but "
                "the section and the compound is identical — same mesh, same "
                "modes, same PK, same exposure window, same arm. Run the SAME "
                "compound at the SAME spec on both specimens. "
                f"The tissue contrast below is still readable: {n_readable} of "
                f"{n_comp} composition row(s) are flagged comparable. "
                + RELATIVE_EXPOSURE_CLAIM
            ),
        }

    plural = "" if n_matched == 1 else "s"
    text = (
        f"{n_matched} compound{plural} matched across {len(compared)} specimens "
        f"at one identical RunSpec, and every comparison below is "
        f"SAME-COMPOUND: one compound's value on one specimen against THAT SAME "
        f"COMPOUND's value on the other, read on {metric}. No ordering of "
        "compounds is computed anywhere in this document. The primary result of "
        "a cross-specimen comparison is the specimen capability contrast — what "
        "each panel can and cannot support — together with the composition "
        f"contrast, of which {n_readable} of {n_comp} row(s) are flagged "
        "comparable. " + RELATIVE_EXPOSURE_CLAIM
    )
    if any_incomparable:
        return {
            "tag": "SAME COMPOUND · NOT COMPARABLE",
            "level": "crit",
            "text": text + (
                "  ABSOLUTE KILLS ARE FLAGGED NOT-COMPARABLE for at least one "
                "compound; see the blockers on each row."
            ),
        }
    return {"tag": "SAME COMPOUND · MATCHED", "level": "ok", "text": text}


def compare_project(project: Dict,
                    resolve_run: Callable[[str], Optional[Dict]],
                    metric: str = "kill_fraction") -> Dict:
    """Cross-specimen comparison of the matched runs in one project.

    ``resolve_run`` maps a ``run_id`` to its ``RunResult`` view (or ``None``).
    Runs are looked up, never duplicated.

    **Same compound, two specimens.** Every comparison this returns holds the
    compound fixed and varies the specimen. Ordering two different compounds by
    predicted effect was withdrawn (module docstring: ``pk_ranking_robustness``
    and ``capacity_ranking_robustness``), and there is no ``ranking`` key.

    The returned document is ordered the way it must be read:

    1. ``specimens`` — every specimen with its capability verdict and the routes
       its panel CANNOT support, so no number below is read without knowing what
       produced it.
    2. ``kill`` — one row per compound matched on two or more specimens: that
       compound's own absolute values side by side with their delta and ratio,
       each row carrying a ``comparable`` flag that is ``False`` unless both
       specimens' panels support a measured cell state and both runs used one.
    3. ``composition`` — cell-type fractions, vascular statistics and measured
       state fractions, each with provenance and its own comparability verdict.

    ``metric`` selects which of ``_HEADLINE_METRICS`` the rows in (2) are
    featured on. It does not decide which rows exist.
    """
    if metric not in _HEADLINE_METRICS:
        raise ProjectError(
            f"metric: {metric!r} is not one of {list(_HEADLINE_METRICS)}")

    specimens = project.get("specimens") or []
    spec_views: List[Dict] = []
    resolved: Dict[str, List[Dict]] = {}
    problems: List[str] = []

    for spec in specimens:
        cap = spec.get("capability") or {}
        runs: List[Dict] = []
        missing: List[str] = []
        not_done: List[Dict] = []
        for rid in spec.get("run_ids") or []:
            run = resolve_run(rid)
            if run is None:
                missing.append(rid)
                continue
            if run.get("status") != "done":
                not_done.append({"run_id": rid, "status": run.get("status"),
                                 "error": run.get("error")})
                continue
            runs.append(run)
        resolved[spec["id"]] = runs

        unsupported = [name for name, blk in (cap.get("state_routes") or {}).items()
                       if blk.get("supported") is False]
        labels: List[str] = []
        if not cap.get("absolute_kill_supported"):
            labels.append("ABSOLUTE KILL NOT SUPPORTED BY THIS PANEL")
        if unsupported:
            labels.append("UNSUPPORTED STATE ROUTE(S): " + ", ".join(sorted(unsupported)))
        if not cap.get("available", True):
            labels.append("SECTION NOT LOADABLE IN THIS CHECKOUT")
        if missing:
            problems.append(
                f"specimen {spec['id']} references run(s) {missing} that no "
                "longer exist in app/runs/")
        spec_views.append({
            "id": spec["id"],
            "label": spec.get("label"),
            "section": spec.get("section"),
            "notes": spec.get("notes"),
            "capability": cap,
            "supports_cell_state": bool(cap.get("supports_cell_state")),
            "absolute_kill_supported": bool(cap.get("absolute_kill_supported")),
            "absolute_kill_reason": cap.get("absolute_kill_reason"),
            "unsupported_state_routes": sorted(unsupported),
            "labels": labels,
            "n_runs_attached": len(spec.get("run_ids") or []),
            "n_runs_done": len(runs),
            "runs_not_done": not_done,
            "runs_missing": missing,
        })

    order = [s["id"] for s in spec_views]
    by_id = {s["id"]: s for s in spec_views}

    # --- group by (spec signature, drug) ------------------------------------ #
    # {signature: {drug: {specimen_id: run}}}
    groups: Dict[str, Dict[str, Dict[str, Dict]]] = {}
    for sid, runs in resolved.items():
        for run in runs:
            drug = _drug_of(run)
            if drug is None:
                continue  # a reconstruction-only run has no compound to compare
            sig = spec_signature(run.get("spec") or {})
            slot = groups.setdefault(sig, {}).setdefault(drug, {})
            prev = slot.get(sid)
            # Newest wins; the earlier one is reported, not silently dropped.
            if prev is None or (run.get("created_at") or 0) > (prev.get("created_at") or 0):
                slot[sid] = run

    def _coverage(sig: str) -> Tuple[int, int]:
        drugs = groups[sig]
        sids = {s for d in drugs.values() for s in d}
        matched = sum(1 for d in drugs.values() if len(d) >= 2)
        return (len(sids), matched)

    signatures = sorted(groups, key=lambda s: (_coverage(s), s), reverse=True)
    primary_sig = signatures[0] if signatures else None

    matched_sets = []
    for sig in signatures:
        drugs = groups[sig]
        matched_sets.append({
            "signature": sig,
            "is_primary": sig == primary_sig,
            "specimens": sorted({s for d in drugs.values() for s in d},
                                key=lambda s: order.index(s) if s in order else 99),
            "drugs": sorted(drugs),
            "n_matched_drugs": sum(1 for d in drugs.values() if len(d) >= 2),
        })

    # --- SAME COMPOUND, TWO SPECIMENS ---------------------------------------- #
    # One row per compound that BOTH (or all) specimens ran at the primary
    # signature. The compound is fixed inside a row and the specimen varies; two
    # different compounds are never placed on the same axis. Rows are still
    # flagged, because an absolute level does not travel between sections even
    # when the compound does (results/heldout_patient.md §3.5).
    kill_rows = []
    any_incomparable = False
    if primary_sig:
        for drug in sorted(groups[primary_sig]):
            per_spec = groups[primary_sig][drug]
            entries = [
                {"specimen": by_id[sid],
                 "capability": by_id[sid]["capability"],
                 "run": run}
                for sid, run in sorted(
                    per_spec.items(),
                    key=lambda kv: order.index(kv[0]) if kv[0] in order else 99)
                if sid in by_id
            ]
            if len(entries) < 2:
                continue
            verdict = _kill_comparability(entries)
            any_incomparable = any_incomparable or not verdict["comparable"]
            ref = entries[0]["specimen"]["id"]
            per_metric = {}
            for m in _HEADLINE_METRICS:
                vals = {}
                for e in entries:
                    v = (e["run"].get("metrics") or {}).get(m)
                    vals[e["specimen"]["id"]] = (
                        float(v) if isinstance(v, (int, float))
                        and not isinstance(v, bool) else None)
                base = vals.get(ref)
                per_metric[m] = {
                    "values": vals,
                    "delta": {k: (v - base) for k, v in vals.items()
                              if k != ref and v is not None and base is not None},
                    "ratio": {k: (v / base) for k, v in vals.items()
                              if k != ref and v is not None and base},
                }
            kill_rows.append({
                "drug": drug,
                "compound": drug,     # the axis held FIXED in this row
                "varies": "specimen",  # ...and the axis that varies
                "reference": ref,
                "runs": {e["specimen"]["id"]: e["run"]["run_id"] for e in entries},
                "metrics": per_metric,
                "comparable": verdict["comparable"],
                "blockers": verdict["blockers"],
                "caveats": verdict["caveats"],
                "state_measured": verdict["state_measured"],
            })

    # --- composition --------------------------------------------------------- #
    comp_entries = []
    for sv in spec_views:
        runs = resolved.get(sv["id"]) or []
        chosen = None
        if primary_sig:
            for run in runs:
                if spec_signature(run.get("spec") or {}) == primary_sig:
                    chosen = run
                    break
        if chosen is None and runs:
            chosen = max(runs, key=lambda r: r.get("created_at") or 0)
        if chosen is None or not (chosen.get("report") or {}):
            continue
        comp_entries.append({
            "specimen": sv, "capability": sv["capability"],
            "report": chosen.get("report") or {}, "run": chosen,
        })
    composition = {
        "specimens": [e["specimen"]["id"] for e in comp_entries],
        "from_run": {e["specimen"]["id"]: e["run"]["run_id"]
                     for e in comp_entries},
        "rows": _composition_rows(comp_entries) if comp_entries else [],
        "note": (
            "Composition is read off the reconstruction report of ONE run per "
            "specimen (named in from_run); every run of a specimen at the same "
            "RunSpec shares one cached reconstruction, so the choice does not "
            "change these numbers. results/heldout_patient.md §6.5: 'the "
            "response layer travels; the reconstruction layer does not' — the "
            "two quantities this project cites as external validations, the "
            "ENDOTHELIAL fraction and the vascular area fraction, are exactly "
            "the ones that failed to replicate on a second patient."
        ),
    }

    head = _headline_block(spec_views, kill_rows, composition["rows"],
                           any_incomparable, metric)

    # The transvascular delivery law that produced every delivery number above.
    # It used to hang off the ranking block; the ranking is gone and the
    # disclosure is not, so it moves to the top level of the document. Read from
    # PDEConfig at call time, never cached (app/server/delivery_law.py).
    law = delivery_law()

    return {
        "project_id": project.get("id"),
        "name": project.get("name"),
        "metric": metric,
        "metric_means": (
            "which metric the SAME-COMPOUND cross-specimen rows are featured on. "
            "It selects a column, not a conclusion: every row carries all of "
            f"{list(_HEADLINE_METRICS)} regardless, and no metric produces an "
            "ordering of compounds."
        ),
        "generated_at": time.time(),
        "specimens": spec_views,
        "matched_sets": matched_sets,
        "primary_signature": primary_sig,
        "delivery_law": law,
        "kill": kill_rows,
        "composition": composition,
        "headline": head["text"],
        "headline_tag": head["tag"],
        "headline_level": head["level"],
        "compares": (
            "the SAME compound across specimens, and the tissue between them. "
            "Ordering two different compounds by predicted effect was withdrawn "
            "(results/pk_ranking_robustness.md: per-compound plasma exposure "
            "spans 113.1x across the compounds the GBM panel ranks and the "
            "published ordering breaks at 1.53x; "
            "results/capacity_ranking_robustness.md: the breast top-1 pick "
            "changes under this repo's own structural capacity rule)."
        ),
        "problems": problems,
        "reading_order": [
            "specimens — what each panel can and cannot support",
            "kill — the SAME compound on two specimens: its own values side by "
            "side, only readable where comparable is true",
            "composition — tissue differences, each with its own verdict",
        ],
    }
