"""Projects, specimens, and the cross-specimen comparison.

The Studio was single-specimen and stateless between sessions: runs lived flat
under ``app/runs/<run_id>/`` with no notion of which patient they belonged to.
This module adds the missing layer, and only that layer:

* A **Project** is a study. It has a name, notes, and one or more **Specimens**.
* A **Specimen** binds a ``sections`` registry id to a *capability report* —
  captured at the moment it is added, from the section's own panel — plus the
  list of ``run_id``s that belong to it. Runs are **associated, never copied**:
  ``app/runs/`` remains the single store, and a specimen holds ids into it.
* ``compare_project`` answers the question the whole product is for: *do two
  patients agree?*

Everything is persisted as one JSON document per project under
``app/projects/<project_id>/project.json`` (gitignored, like ``app/runs/``), so
a study survives a restart.

Why the comparison is shaped the way it is
------------------------------------------

``results/heldout_patient.md`` ran the identical frozen pipeline on two real GBM
sections and found a hard asymmetry, restated in its §6.5 as *"the response layer
travels; the reconstruction layer does not"*:

* The **cross-drug kill RANKING replicates exactly** (cisplatin > doxorubicin >
  lomustine > methotrexate > temozolomide on both sections, §3.3).
* **Absolute kills do not travel.** §3.5 names two mechanisms, neither of them a
  drug property: the absolute oxygen level is uncalibrated on both sections
  (``OxygenConfig.consumption`` was never fitted to a measured tension), and the
  ECM field is a max-normalised KDE over *whichever matrix markers the panel
  happens to carry* — 11 genes on one section, 3 on the other — so the two have
  no common absolute scale even on identical tissue. Its §5 table therefore
  labels the ECM row "NOT comparable as biology".
* Both external reconstruction validations **failed** on the second patient: the
  ENDOTHELIAL fraction (9.74% -> 4.18% against a volume-EM 8.0% anchor) and the
  vascular area fraction (2.50% -> 5.79% against a measured 2.1 +/- 0.3%).

``results/cell_state_impact.md`` supplies the capability rule. Its §3.2 measures
absolute kill moving by up to **0.3252 absolute (1.77x relative, i.e. 77%)**
purely with the choice of cell-state route, with tissue, drug, transport and RNG
held identical, while §3.3 finds the ranking **INVARIANT** across all three
routes. Its §5 conclusion is that a panel carrying no hypoxia-response
transcript cannot support an absolute kill under *any* route — the default
route's HYPOXIC band is empty for a bookkeeping reason, the marker route
declares hypoxia UNMEASURED, and the vessel-distance route is a geometric proxy.

So this module makes **rank agreement the primary result**, and refuses to
present a cross-specimen absolute-kill delta as a finding unless *both*
specimens' panels can support a measured cell state *and* both runs actually
used a measured route. A measured kill is never silently differenced against an
assumed one.
"""
from __future__ import annotations

import json
import os
import re
import threading
import time
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from .cache import (
    all_sections,
    describe_section,
    repo_root,
    state_route_support,
)
from .delivery_law import delivery_law

__all__ = [
    "PROJECT_ID_RE",
    "ProjectError",
    "ProjectStore",
    "PROJECT_STORE",
    "projects_dir",
    "capability_report",
    "spec_signature",
    "compare_project",
]

PROJECT_ID_RE = re.compile(r"^p_\d{5}$")
SPECIMEN_ID_RE = re.compile(r"^s_\d{5}$")

#: The figures this module quotes, kept as named constants so the message and
#: the evidence cannot drift apart. Both are from results/cell_state_impact.md
#: §3.2 (kill across the three state routes, everything else held identical).
STATE_ROUTE_KILL_SWING_ABS = 0.3252
STATE_ROUTE_KILL_SWING_REL = 1.77


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
    """Absolute path to ``app/projects/`` (created on demand)."""
    path = os.path.join(repo_root(), "app", "projects")
    os.makedirs(path, exist_ok=True)
    return path


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
    clears the bar and the 0.3252-absolute swing (§3.2) is unresolved. The
    cross-drug RANKING is invariant across all three routes and is what survives.

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
            "route, with tissue, drug, transport and RNG held identical. The "
            "CROSS-DRUG RANKING is invariant across all three routes and is the "
            "defensible output; the level is not."
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

    def __init__(self, root: Optional[str] = None):
        self.root = root or projects_dir()
        os.makedirs(self.root, exist_ok=True)
        self._lock = threading.RLock()

    # -- paths -------------------------------------------------------------- #
    def _dir(self, project_id: str) -> str:
        return os.path.join(self.root, project_id)

    def _path(self, project_id: str) -> str:
        return os.path.join(self._dir(project_id), "project.json")

    # -- io ----------------------------------------------------------------- #
    def _read(self, project_id: str) -> Optional[Dict]:
        try:
            with open(self._path(project_id)) as fh:
                doc = json.load(fh)
        except FileNotFoundError:
            return None
        except Exception as exc:  # noqa: BLE001 - a corrupt file is reported
            raise ProjectError(
                f"project {project_id!r} is on disk but unreadable: "
                f"{type(exc).__name__}: {exc}",
                status=500,
            ) from exc
        if not isinstance(doc, dict):
            raise ProjectError(
                f"project {project_id!r} is not a JSON object", status=500)
        doc.setdefault("id", project_id)
        doc.setdefault("specimens", [])
        return doc

    def _write(self, doc: Dict) -> Dict:
        doc["updated_at"] = time.time()
        d = self._dir(doc["id"])
        os.makedirs(d, exist_ok=True)
        tmp = os.path.join(d, "project.json.tmp")
        with open(tmp, "w") as fh:
            json.dump(doc, fh, indent=1, sort_keys=False)
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
    def list(self) -> List[Dict]:
        """Every project, newest first, as a summary row."""
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
                out.append(self._summary(doc))
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

    def create(self, name: str, notes: str = "") -> Dict:
        name = (name or "").strip()
        if not name:
            raise ProjectError("name: required, a non-empty project name")
        with self._lock:
            existing = self._ids()
            n = 0
            for pid in existing:
                try:
                    n = max(n, int(pid.split("_")[-1]))
                except ValueError:
                    pass
            project_id = f"p_{n + 1:05d}"
            doc = {
                "id": project_id,
                "name": name,
                "notes": str(notes or ""),
                "created_at": time.time(),
                "updated_at": time.time(),
                "specimens": [],
            }
            return self._write(doc)

    def update(self, project_id: str, name: Optional[str] = None,
               notes: Optional[str] = None) -> Dict:
        with self._lock:
            doc = self.get(project_id, refresh_capability=False)
            if name is not None:
                if not str(name).strip():
                    raise ProjectError("name: must not be empty")
                doc["name"] = str(name).strip()
            if notes is not None:
                doc["notes"] = str(notes)
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
                if not str(label).strip():
                    raise ProjectError("label: must not be empty")
                spec["label"] = str(label).strip()
            if notes is not None:
                spec["notes"] = str(notes)
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
#: thing that differs between the two sides of a comparison is the patient.
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
_SPEC_IGNORED = frozenset({"section", "drug", "pk"})


def _round(v):
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return v
    return round(float(v), 6)


def spec_signature(spec: Dict) -> str:
    """Canonical string for everything in a RunSpec except section and drug.

    Two runs with the same signature differ only in *which patient* and *which
    drug* — which is exactly the pair of axes a cross-specimen therapy panel
    varies. Anything else differing (a different mesh, a different state route,
    a different exposure window, a counterfactual arm) makes the two runs
    incomparable, and grouping by this signature is what stops them being
    compared anyway.

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
# Ranking
# --------------------------------------------------------------------------- #
def _order(values: Dict[str, float]) -> List[str]:
    """Drugs ordered by value, descending; ties broken by name so the string
    form is stable. Ties are reported separately — they are not hidden here."""
    return sorted(values, key=lambda d: (-values[d], d))


def _kendall(a: Dict[str, float], b: Dict[str, float],
             drugs: Sequence[str]) -> Dict:
    """Kendall tau-b over the drugs common to both, plus the discordant pairs.

    Written out rather than imported so the discordant pairs — the *useful*
    output when a ranking disagrees, because they name which two drugs swapped —
    fall out of the same pass.
    """
    conc = disc = tie_a = tie_b = 0
    discordant: List[List[str]] = []
    for i in range(len(drugs)):
        for j in range(i + 1, len(drugs)):
            x, y = drugs[i], drugs[j]
            da = a[x] - a[y]
            db = b[x] - b[y]
            if da == 0 and db == 0:
                tie_a += 1
                tie_b += 1
                continue
            if da == 0:
                tie_a += 1
                continue
            if db == 0:
                tie_b += 1
                continue
            if (da > 0) == (db > 0):
                conc += 1
            else:
                disc += 1
                discordant.append([x, y])
    n0 = conc + disc
    denom = ((n0 + tie_a) * (n0 + tie_b)) ** 0.5
    tau = ((conc - disc) / denom) if denom > 0 else None
    return {
        "kendall_tau_b": tau,
        "concordant_pairs": conc,
        "discordant_pairs": disc,
        "ties": {"a": tie_a, "b": tie_b},
        "swapped": discordant,
    }


def _rank_block(specimen_values: Dict[str, Dict[str, float]],
                specimen_order: Sequence[str],
                metric: str) -> Dict:
    """The PRIMARY result: does the drug ranking agree across specimens?

    Restricted to the drugs every compared specimen has a run for — ranking a
    5-drug panel against a 3-drug one is not a disagreement, it is a different
    question, and mixing them would manufacture one.
    """
    sids = [s for s in specimen_order if specimen_values.get(s)]
    common = None
    for s in sids:
        keys = set(specimen_values[s])
        common = keys if common is None else (common & keys)
    common_drugs = sorted(common or set())

    by_specimen = {
        s: {
            "order": _order({d: specimen_values[s][d] for d in common_drugs}),
            "values": {d: specimen_values[s][d] for d in common_drugs},
        }
        for s in sids
    }

    # Which transvascular delivery law produced every value in this block. Read
    # from PDEConfig at call time, never cached, because the ordering below is
    # law-dependent: it is invariant to the COMMON-MODE perturbations the
    # verdict cites (cell-state route, section) and NOT to this one, which is
    # drug-specific and changes top-1 on GBM. See app/server/delivery_law.py.
    law = delivery_law()

    if len(sids) < 2 or len(common_drugs) < 2:
        return {
            "metric": metric,
            "specimens": sids,
            "drugs": common_drugs,
            "by_specimen": by_specimen,
            "delivery_law": law,
            "agrees": None,
            "pairwise": [],
            "verdict": (
                "RANKING NOT EVALUABLE: a rank comparison needs at least 2 "
                f"specimens with at least 2 drugs in common, and this project has "
                f"{len(sids)} specimen(s) with {len(common_drugs)} common drug(s) "
                "at a matched RunSpec. Run the same drugs at the same spec on "
                "both specimens."
            ),
            "reason": "insufficient matched runs",
        }

    pairwise = []
    all_agree = True
    for i in range(len(sids)):
        for j in range(i + 1, len(sids)):
            a, b = sids[i], sids[j]
            k = _kendall(by_specimen[a]["values"], by_specimen[b]["values"],
                         common_drugs)
            agrees = by_specimen[a]["order"] == by_specimen[b]["order"]
            all_agree = all_agree and agrees
            pairwise.append({
                "a": a, "b": b, "agrees": agrees,
                "order_a": by_specimen[a]["order"],
                "order_b": by_specimen[b]["order"],
                **k,
            })

    verdict = (
        "RANKING AGREES across every specimen pair. This is the primary result "
        "of a cross-patient comparison and the one the evidence supports: "
        "results/heldout_patient.md §3.3 found the kill ordering IDENTICAL on "
        "two real GBM sections, and results/cell_state_impact.md §3.3 found it "
        "INVARIANT across all three cell-state routes. The absolute levels did "
        "NOT travel in either study — read the kill table below with its "
        "comparability flags, not as a difference between two patients."
        if all_agree else
        "RANKING DISAGREES on at least one specimen pair. The swapped drug pairs "
        "are listed per pair. A disagreement here is a real finding — it is the "
        "one quantity results/heldout_patient.md §3.3 and "
        "results/cell_state_impact.md §3.3 both found to be invariant, so a "
        "swap means either a genuine patient difference in delivery x potency or "
        "a pair of drugs whose predicted kills are close enough that the order "
        "is not resolved. Check the margin between the swapped drugs before "
        "reading it as biology."
    )
    return {
        "metric": metric,
        "specimens": sids,
        "drugs": common_drugs,
        "by_specimen": by_specimen,
        "delivery_law": law,
        "agrees": all_agree,
        "pairwise": pairwise,
        "verdict": verdict,
        "reason": None,
    }


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
    """Is an absolute kill difference between these specimens readable at all?

    ``entries`` are ``{specimen, capability, run}`` for one matched drug. The
    verdict is ``False`` unless EVERY specimen's panel can support a measured
    cell state AND every run actually used a measured route — the rule from
    ``results/cell_state_impact.md`` §5, and the reason
    ``results/heldout_patient.md`` §3.5 gives for refusing to compare absolute
    kills between its two sections.
    """
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

    return {
        "comparable": not blockers,
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
         "in results/heldout_patient.md §1.1 the fraction FELL 0.43x while this "
         "density ROSE 2.67x on the same pair of sections. They point in "
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
_HEADLINE_METRICS = (
    "kill_fraction", "mean_tumor_auc", "therapeutic_coverage_1x",
    "penetration_length_um",
)


def compare_project(project: Dict,
                    resolve_run: Callable[[str], Optional[Dict]],
                    metric: str = "kill_fraction") -> Dict:
    """Cross-specimen comparison of the matched runs in one project.

    ``resolve_run`` maps a ``run_id`` to its ``RunResult`` view (or ``None``).
    Runs are looked up, never duplicated.

    The returned document is ordered the way it must be read:

    1. ``specimens`` — every specimen with its capability verdict and the routes
       its panel CANNOT support, so no number below is read without knowing what
       produced it.
    2. ``ranking`` — **the primary result.** Does the drug ordering agree?
    3. ``kill`` — the matched per-drug absolute values and their deltas, each
       carrying a ``comparable`` flag that is ``False`` unless both specimens'
       panels support a measured cell state and both runs used one.
    4. ``composition`` — cell-type fractions, vascular statistics and measured
       state fractions, each with provenance and its own comparability verdict.
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
                continue  # a reconstruction-only run has no drug to rank
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

    # --- ranking (PRIMARY) --------------------------------------------------- #
    values: Dict[str, Dict[str, float]] = {}
    if primary_sig:
        for drug, per_spec in groups[primary_sig].items():
            for sid, run in per_spec.items():
                v = (run.get("metrics") or {}).get(metric)
                if isinstance(v, (int, float)) and not isinstance(v, bool):
                    values.setdefault(sid, {})[drug] = float(v)
    ranking = _rank_block(values, order, metric)
    ranking["signature"] = primary_sig

    # --- absolute per-drug table (SECONDARY, flagged) ------------------------ #
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

    headline = ranking["verdict"]
    if any_incomparable:
        headline += (
            "  ABSOLUTE KILLS ARE FLAGGED NOT-COMPARABLE for at least one drug; "
            "see the blockers on each row."
        )

    return {
        "project_id": project.get("id"),
        "name": project.get("name"),
        "metric": metric,
        "generated_at": time.time(),
        "specimens": spec_views,
        "matched_sets": matched_sets,
        "primary_signature": primary_sig,
        "ranking": ranking,
        "kill": kill_rows,
        "composition": composition,
        "headline": headline,
        "problems": problems,
        "reading_order": [
            "specimens — what each panel can and cannot support",
            "ranking — THE PRIMARY RESULT: does the drug ordering agree?",
            "kill — absolute values, only readable where comparable is true",
            "composition — tissue differences, each with its own verdict",
        ],
    }
