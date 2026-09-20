"""FastAPI application implementing API_CONTRACT.md v1.

Route map (everything in the contract, nothing invented on top of it):

===========  ==================================  ==================================
Method       Path                                Returns
===========  ==================================  ==================================
GET          /api/sections                       Section[]
GET          /api/drugs                          drug rows with real IC50 + source
GET          /api/defaults                       a valid RunSpec (production values)
POST         /api/runs                           {run_id} (enqueues, returns at once)
GET          /api/runs                           RunResult[] (newest first)
GET          /api/runs/{id}                      RunResult
GET          /api/runs/{id}/status               {status, progress, stage,
                                                  elapsed_s, log_tail}
DELETE       /api/runs/{id}                      cancel a queued/running job
GET          /api/runs/{id}/geometry             TMSG binary cell payload
GET          /api/runs/{id}/field/{name}         TMSF binary field payload
GET          /api/runs/{id}/mechanism            TMSM per-cell mechanism labels
                                                 + the aggregate breakdown
GET          /api/runs/{id}/timelapse            TMSL drug-concentration frames
POST         /api/validate                       {ok, errors[], warnings[]}
POST         /api/sections/scan                  ingest ANY directory (or, with
                                                 TMESIM_STUDIO_DATA_ROOTS set,
                                                 any directory under an
                                                 allowlisted root): detect() +
                                                 capability report; optionally
                                                 register it as a section
DELETE       /api/sections/{id}                  remove an INGESTED section
GET/POST     /api/projects                       Project[] / create a Project
GET/PATCH/   /api/projects/{id}                  one Project (with its Specimens)
DELETE
POST         /api/projects/{id}/specimens        bind a section as a Specimen
GET/PATCH/   /api/projects/{id}/specimens/{sid}  one Specimen + its capability
DELETE
POST         /api/projects/{id}/specimens/       associate an existing run
             {sid}/runs
DELETE       /api/projects/{id}/specimens/       detach it (the run is kept)
             {sid}/runs/{run_id}
GET          /api/projects/{id}/compare          cross-specimen comparison
POST         /api/reports                        {report_id} (2 runs per drug, enqueued)
GET          /api/reports                        ReportSummary[] (newest first)
GET          /api/reports/{id}                   the full Report document (poll target)
GET          /api/reports/{id}/html              the rendered, self-contained document
DELETE       /api/reports/{id}                   cancel a queued/running report
===========  ==================================  ==================================

Two additive diagnostics that the contract does not forbid and the UI does not
have to use: ``GET /api/health`` (cache + queue counters) and
``GET /api/runs/{id}/traceback`` (the failure stack for a ``failed`` run).

One additive resource, ``/api/compounds`` (see ``app/server/compounds.py``):

===========  ==================================  ==================================
POST         /api/compounds/compile              dry-run compile of a user SMILES:
                                                 the DrugParams a run would
                                                 simulate + a provenance row per
                                                 field. 422 with a specific
                                                 refusal reason, never a stack
                                                 trace. Stores nothing.
GET/POST     /api/compounds                      stored user compounds / compile
                                                 and store one
GET/DELETE   /api/compounds/{id}                 recompile one from its stored
                                                 inputs / forget it
===========  ==================================  ==================================

A stored compound's name becomes a legal ``RunSpec.drug`` (see ``_known_drugs``),
which is the whole point: a researcher maps THEIR compound on THEIR section and
reads where it reaches. (This line said "ranks THEIR series" until 2026-08-21;
compound ordering was measured and withdrawn - see
``results/pk_ranking_robustness.md`` and ``results/capacity_ranking_robustness.md``.)

The static front end in ``app/web/`` is mounted last, at ``/``, so it can never
shadow an ``/api`` route.
"""
from __future__ import annotations

import math
import os
import secrets
import threading
import time
from contextlib import asynccontextmanager, contextmanager
from typing import Dict, List, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from starlette.middleware.trustedhost import TrustedHostMiddleware
from starlette.routing import get_route_path

from .cache import (
    OXYGEN_MODES,
    RECON_CACHE,
    SECTIONS,
    PathNotAllowed,
    all_sections,
    confine_scan_path,
    lookup_section,
    SENSITIVITY_MODES,
    STATE_MODES,
    VASCULAR_MODES,
    list_sections,
    make_key,
    redact_exception,
    repo_root,
    section_genes,
    state_route_support,
)
from .jobs import PDE_DT_S, JobManager, jsonable
from .payloads import FIELD_NAMES, MAX_TIMELAPSE_FRAMES
from .projects import (
    PROJECT_STORE,
    ProjectError,
    ProjectStore,
    capability_report,
    compare_project,
)
from .report import (
    MAP_CLAIM,
    REPORT_SENSITIVITY_MODE,
    ReportManager,
    render_html,
)

__all__ = ["create_app", "app", "validate_spec", "drug_rows", "default_spec",
           "validate_report_request"]

#: Grid resolutions this repo has ever run a convergence LADDER at. These are
#: NOT "converged below this" thresholds and must not be described as such: the
#: claim that kill_fraction is grid-converged for dx <= 20 um is withdrawn (see
#: the correction block in results/numerical_rigor.md, and the two ladders that
#: contradict it — results/grid_convergence_breast.md, .../grid_convergence_gbm.md).
#: They survive as the coarseness at which the measured drift becomes large
#: enough to warn about, which is all they were ever entitled to be.
DX_KILL_LADDER_COARSE_UM = 20.0
DX_PENETRATION_LADDER_COARSE_UM = 15.0

#: Kept as aliases so nothing importing the old names breaks; the names lie
#: about what was measured, so nothing in this file uses them.
DX_KILL_CONVERGED_UM = DX_KILL_LADDER_COARSE_UM
DX_PENETRATION_CONVERGED_UM = DX_PENETRATION_LADDER_COARSE_UM

#: Above this the grid is big enough that a run is measured in minutes; the UI
#: gets a warning rather than a refusal, because it is a cost, not an error.
_LARGE_GRID_CELLS = 1_500_000

#: SECURITY / capacity: hard REFUSAL ceiling on the grid, ~27x the warning
#: threshold. Grid size scales as 1/dx^2, and dx_um used to be bounded only by
#: "must be > 0" -- everything above that was a warning, and warnings do not
#: block (``validate_spec`` returns ``ok = not errors``). On a 5000 x 5000 um
#: section, dx_um = 0.1 asks for 2.5e9 voxels: one float64 field is 20 GB, and
#: several fields are allocated per solve. One POST could OOM the box.
#: 40 M voxels is ~320 MB per float64 field: expensive but survivable, and far
#: beyond any dx a real section needs.
#:
#: OVERRIDABLE, because 40 M is the "no request should ever ask for this"
#: ceiling and NOT the "this machine can afford it" one. Those are different
#: numbers and only the operator knows the second.
#:
#: Measured on a 512 MB instance (Render's free tier), flagship Xenium GBM,
#: real container, no bind mount: dx = 30 um is 46,070 voxels and peaks at
#: 220 MiB -- comfortable. dx = 12 um is 286,200 voxels, ~0.7 % of the ceiling
#: below, and the container is OOM-KILLED. Not the run failing: the PROCESS
#: dies, status=exited, and /healthz stops answering, so one ordinary request
#: takes the service down for everyone. A ceiling that only catches 20 GB
#: allocations does not catch that, because 286,200 voxels is a perfectly
#: reasonable ask on a machine with room for it.
#:
#: So a constrained host sets TMESIM_MAX_GRID_CELLS to what it can actually
#: afford and gets a 400 with a reason instead of a dead container. 100_000
#: permits dx >= 20.25 um on this section and refuses dx = 12 -- and
#: refuses dx = 20 TOO, which this line used to say it allowed. Measured
#: 2026-08-30 with this file's own extent formula on xenium_319:
#: dx = 20 is 102,419 voxels, 2.4 % OVER the ceiling. Break-even is
#: dx = 20.2405 um. An operator sizing a host off '~20' would have set a
#: ceiling that refuses the resolution they picked it for.
#: ``TMESIM_MAX_GRID_CELLS=auto`` adopts the largest grid MEASURED to survive
#: on a host of this size -- 46,070 voxels, the dx = 30 um point above, which
#: peaked at 220 MiB on the 512 MB instance. It applies a MEASUREMENT; it does
#: not fit one. The only two points available are "46,070 survived" and
#: "286,200 died", so the only honest automatic ceiling is the one that
#: survived, and it is deliberately conservative: it refuses grids between
#: those points that might well have been fine. An operator who wants that
#: range back sets a number instead.
#:
#: OPT-IN, and that is the whole design. The ceiling could have moved on its
#: own the moment a small cgroup limit was detected, and it does not: a server
#: that silently refuses a spec the operator expects to run, because it read a
#: number out of /sys, is a debugging session nobody asked for. Detection
#: drives the WARNING unconditionally (see account.oom_warning); it drives the
#: CEILING only when asked to, by name.
_DEFAULT_MAX_GRID_CELLS = 40_000_000


def _resolve_max_grid_cells():
    """(ceiling, source). See the note above for why ``auto`` is opt-in."""
    raw = (os.environ.get("TMESIM_MAX_GRID_CELLS") or "").strip()
    if not raw:
        return _DEFAULT_MAX_GRID_CELLS, "default"
    if raw.lower() == "auto":
        from .account import (MEASURED_INSTANCE_MB, MEASURED_OK_VOXELS,
                              memory_limit)
        mem = memory_limit()
        if mem["bytes"] and mem["bytes"] <= MEASURED_INSTANCE_MB * 1024 * 1024:
            return MEASURED_OK_VOXELS, f"auto ({mem['source']})"
        # `auto` on a host that is NOT small, or whose size is undeclared,
        # falls back to the default rather than guessing a middle value --
        # and reports that it did, so "auto" never hides which branch ran.
        return _DEFAULT_MAX_GRID_CELLS, "auto (no small-host limit detected)"
    try:
        return int(raw), "TMESIM_MAX_GRID_CELLS"
    except ValueError:
        # A typo must not silently restore the 40 M default on a small box.
        raise SystemExit(
            f"TMESIM_MAX_GRID_CELLS={raw!r} is neither an integer nor 'auto'. "
            "Refusing to start rather than fall back to the 40,000,000 "
            "default, which on a small host permits the spec measured to "
            "OOM-kill it.")


_MAX_GRID_CELLS, _MAX_GRID_CELLS_SOURCE = _resolve_max_grid_cells()

#: App-wide request-body ceiling, in bytes.
#:
#: `app/server/protocol_api.py` has had a 4 MB cap with a Content-Length
#: pre-check since a ~500 MB POST cost ~2.8 GB resident and a handful of them
#: OOM'd the process. Its docstring makes the argument in one line: "a limit
#: that fires after the allocation is not a limit."
#:
#: It applied to ONE router. Every other POST route -- /api/runs,
#: /api/validate, /api/reports, /api/projects, /api/compounds,
#: /api/sections/scan -- did a bare `await request.json()` and buffered the
#: lot. MEASURED 2026-09-03: a 64 MB body was buffered and parsed by all six,
#: in 0.07-0.45 s each, while the protocol router refused the same body in
#: 0.009 s.
#:
#: Same number as the protocol router's, deliberately: the largest legitimate
#: body here is a RunSpec or a compound definition, both kilobytes.
MAX_REQUEST_BODY_BYTES = 4 * 1024 * 1024

#: Ceiling on pk.n_doses. An ERROR, not a warning, for the reason the grid
#: ceiling above is one: ``ok = not errors``, so a warning does not block the
#: allocation it warns about.
#:
#: ``treatment/drugs.py::plasma_profile`` builds ``np.arange(n_doses)`` plus
#: broadcast arrays of shape ``(n_doses, len(times))``. MEASURED on this
#: machine with only THREE time points: n_doses = 2,000,000 peaks at 160 MB in
#: 56 ms, and it is linear -- so 1e9 doses is ~80 GB for three points, and far
#: more over a real time grid. ``POST /api/runs`` accepted n_doses = 10**15
#: with a 202.
#:
#: This is the same class as the OOM the grid ceiling was added for: one
#: ordinary POST that takes the process down for everyone. The grid was capped
#: and the dose count was not. 10,000 is far above any real schedule -- a full
#: Stupp course is 30-42 administrations.
MAX_N_DOSES = 10_000

#: Ceiling on total_time_s, in the same spirit. 100 years.
#: POST /api/runs accepted total_time_s = 1e308 with a 202 before this. The window already
#: warns above 72 h; this only refuses the values that are a memory request
#: rather than a protocol.
MAX_TOTAL_TIME_S = 100.0 * 365.25 * 24 * 3600

#: Hard floor on dx_um for the case where the section extent is unknown and
#: the voxel ceiling above cannot be evaluated. 0.5 um is well below a cell
#: diameter, so no legitimate tissue-scale run is refused by it.
MIN_DX_UM = 0.5

#: Wall-clock at import, used for the uptime figure on /api/account.
#:
#: At IMPORT, not at create_app(): the tests build several apps in one process
#: and a per-app boot time would report the age of the most recent test
#: fixture. What an operator wants from "uptime" is how long THIS PROCESS has
#: been answering, and that is what this measures. It is time.time(), so it is
#: a wall clock and it moves if the host clock is adjusted -- fine for an
#: order-of-magnitude "has this restarted?", not a duration to compute against.
_PROCESS_STARTED_AT = time.time()


# --------------------------------------------------------------------------- #
# Drugs
# --------------------------------------------------------------------------- #
_drug_rows_cache: Optional[List[Dict]] = None

#: SINGLE-FLIGHT. Compiling the library is one job, and it must be done by one
#: thread at a time no matter how many callers arrive.
#:
#: MEASURED ON THE LIVE DEPLOYMENT 2026-09-20. `GET /api/drugs` did not answer
#: inside 120 s, repeatedly, while `/healthz` answered the same container in
#: 0.2 s and the same image answered in 0.012 s on a laptop. The memo made
#: this look impossible, but the memo is only consulted at ENTRY: before the
#: first result exists, every concurrent caller ran the whole compile
#: independently. The startup warm thread is one of those callers.
#:
#: On a free tier with a fraction of a CPU, that is a pile-up that feeds
#: itself: a client times out, the server keeps computing, the client retries,
#: and the retry starts another full compile while the abandoned ones still
#: hold the CPU. The Studio's boot does
#: `Promise.all([getSections, getDrugs, getDefaults])`, so the whole
#: workstation waits on the slowest of them.
#:
#: With the lock, N callers cost one compile: the first computes, the rest
#: block and then read the memo. This changes no result -- `drug_rows()`
#: returns the same rows either way -- it changes only how many times the
#: machine is asked to produce them.
_drug_rows_lock = threading.Lock()


def drug_rows(force: bool = False) -> List[Dict]:
    """``GET /api/drugs`` payload: the real library, compiled once and memoised.

    ``ic50_uM`` and ``ic50_source`` come from
    ``molecular.potency.resolve_ic50`` — measured GDSC glioma median where the
    drug is screened, curated literature where it is not, and a flagged default
    otherwise. ``bbb_penetrance`` comes from the fully compiled
    ``calibrated_drugparams -> apply_barrier_to_drugparams`` chain, i.e. the
    exact ``DrugParams`` a run would simulate. A drug that fails to compile is
    returned with nulls and a ``reason``, never dropped and never faked.

    Rows for compounds the operator compiled from their own SMILES
    (``app/server/compounds.py``) are appended on EVERY call, never memoised:
    the library is fixed for the life of the process, a user's compound list is
    not. Those rows carry ``is_user_compound: true`` and their own
    ``predicted_fields`` count, because the picker must be able to say that a
    self-compiled compound is mostly predictions before it is chosen, not after.
    """
    global _drug_rows_cache
    if _drug_rows_cache is not None and not force:
        return _drug_rows_cache + _user_compound_rows()

    # Serialise the compile. The re-check inside the lock is what makes this
    # single-flight rather than merely serialised: everyone who queued while
    # the first caller worked finds the memo populated and returns it instead
    # of repeating work that is already done.
    #
    # `_user_compound_rows()` stays OUTSIDE the memo, exactly as before: the
    # library is fixed for the life of the process, an operator's own compiled
    # compounds are not.
    with _drug_rows_lock:
        if _drug_rows_cache is not None and not force:
            return _drug_rows_cache + _user_compound_rows()
        return _compile_drug_rows()


def _compile_drug_rows() -> List[Dict]:
    """The compile itself. Split out so the memo/lock dance above reads as one
    thing and the body below is unchanged from when it was inline."""
    global _drug_rows_cache

    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.breast_library import get_molecule_breast
    from tmesim.molecular.library import GBM_DRUGS, get_molecule
    from tmesim.molecular.potency import resolve_ic50
    from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

    def _get_molecule_any(name: str):
        try:
            return get_molecule(name)
        except KeyError:
            if name in _CROSS_LIBRARY_GBM_DRUGS:
                return get_molecule_breast(name)
            raise

    try:
        alpha = float(calibrate_alpha())
        alpha_error = None
    except Exception as exc:  # noqa: BLE001
        alpha, alpha_error = None, redact_exception(
            exc, "calibrate_alpha() failed", brief=True)

    rows: List[Dict] = []
    for name in sorted(set(GBM_DRUGS) | set(_CROSS_LIBRARY_GBM_DRUGS)):
        row: Dict = {
            "name": name,
            "mw": None,
            "logp": None,
            "ic50_uM": None,
            "ic50_source": None,
            "bbb_penetrance": None,
            # --- additive provenance --------------------------------------- #
            "ic50_confidence": None,
            "ic50_note": None,
            "ic50_n_cell_lines": None,
            "is_antibody": None,
            "smiles_source": None,
            "diffusivity_um2_s": None,
            "reason": None,
        }
        try:
            mol = _get_molecule_any(name)
            row["mw"] = _num(mol.get("mw"))
            row["logp"] = _num(mol.get("logp"))
            row["is_antibody"] = bool(mol.get("is_antibody"))
            row["smiles_source"] = mol.source
        except Exception as exc:  # noqa: BLE001
            row["reason"] = redact_exception(
                exc, "the molecular descriptors for this drug could not be built")

        try:
            pot = resolve_ic50(name)
            row["ic50_uM"] = _num(pot["ic50_uM"])
            row["ic50_confidence"] = pot["confidence"]
            row["ic50_note"] = pot["note"]
            row["ic50_n_cell_lines"] = pot.get("n_cell_lines")
            if pot["confidence"] == "measured":
                row["ic50_source"] = f"GDSC measured (n={pot.get('n_cell_lines')})"
            elif pot["confidence"] == "literature":
                row["ic50_source"] = "curated literature"
            else:
                row["ic50_source"] = "flagged default (LOW confidence)"
        except Exception as exc:  # noqa: BLE001
            row["reason"] = redact_exception(
                exc, "the IC50 for this drug could not be resolved")

        if alpha is None:
            row["reason"] = f"bbb_penetrance unavailable: {alpha_error}"
        else:
            try:
                dp = apply_barrier_to_drugparams(
                    calibrated_drugparams(name, alpha), _get_molecule_any(name))
                row["bbb_penetrance"] = _num(dp.bbb_penetrance)
                row["diffusivity_um2_s"] = _num(dp.diffusivity)
                # The compiled ic50 is the resolved one; report it from the
                # object the engine will actually consume.
                row["ic50_uM"] = _num(dp.ic50)
            except Exception as exc:  # noqa: BLE001
                row["reason"] = redact_exception(
                    exc, "this drug could not be compiled")
        rows.append(row)

    _drug_rows_cache = rows
    return rows + _user_compound_rows()


def _user_compound_rows() -> List[Dict]:
    """``GET /api/drugs`` rows for the operator's own compiled compounds.

    Same schema as a library row so the picker needs no second code path, plus
    the fields that make the difference visible at the point of choice:
    ``is_user_compound``, ``smiles``, ``predicted_fields`` / ``total_fields``,
    and ``ic50_supplied``. A compound that no longer compiles (a physics change
    moved a range bound under it) is returned with nulls and a ``reason`` —
    never dropped, so the operator can see why it went away.
    """
    from .compounds import COMPOUND_STORE

    rows: List[Dict] = []
    try:
        stored = COMPOUND_STORE.list()
    except Exception as exc:  # noqa: BLE001 - a broken store must not 500 /api/drugs
        return [{
            "name": "(user compounds unavailable)", "mw": None, "logp": None,
            "ic50_uM": None, "ic50_source": None, "bbb_penetrance": None,
            "is_user_compound": True, "unavailable": True,
            "reason": redact_exception(
                exc, "the user compound store could not be read"),
        }]
    for rec in stored:
        row: Dict = {
            "name": rec.get("name"),
            "mw": None, "logp": None, "ic50_uM": None, "ic50_source": None,
            "bbb_penetrance": None, "ic50_confidence": None, "ic50_note": None,
            "ic50_n_cell_lines": None, "is_antibody": False,
            "smiles_source": "user_smiles", "diffusivity_um2_s": None,
            "reason": None,
            "is_user_compound": True,
            "id": rec.get("id"),
            "smiles": rec.get("smiles"),
            "note": rec.get("note"),
            "created_at": rec.get("created_at"),
            "ic50_supplied": rec.get("ic50_uM") is not None,
        }
        try:
            doc = COMPOUND_STORE.compile_stored(rec.get("id") or rec.get("name"))
        except Exception as exc:  # noqa: BLE001
            row["reason"] = redact_exception(
                exc, "this compound no longer compiles")
            rows.append(row)
            continue
        dp = doc["drug_params"]
        row.update({
            "mw": _num(dp.get("molecular_weight")),
            "logp": _num(doc.get("descriptors", {}).get("logp")),
            "ic50_uM": _num(dp.get("ic50_uM")),
            "bbb_penetrance": _num(dp.get("bbb_penetrance")),
            "diffusivity_um2_s": _num(dp.get("diffusivity_um2_s")),
            "ic50_confidence": "user_supplied" if row["ic50_supplied"] else "default_LOW",
            "ic50_source": ("supplied by you (not measured here)"
                            if row["ic50_supplied"]
                            else "flagged default (LOW confidence) — you supplied none"),
            "ic50_note": next((f["method"] for f in doc.get("fields", [])
                               if f["key"] == "ic50"), None),
            "fields": doc.get("fields"),
            "caveats": doc.get("caveats"),
            "warnings": doc.get("warnings"),
            "can_conclude": doc.get("can_conclude"),
            "cannot_conclude": doc.get("cannot_conclude"),
            "predicted_fields": sum(1 for f in doc.get("fields", [])
                                    if f.get("is_prediction")),
            "total_fields": len(doc.get("fields", [])),
        })
        rows.append(row)
    return rows


async def _json_body(request: Request) -> Dict:
    """The request body as a JSON object, or a 400 saying exactly what was wrong."""
    try:
        body = await request.json()
    except Exception as exc:  # noqa: BLE001
        raise HTTPException(
            status_code=400,
            detail=f"request body is not valid JSON: {exc}") from exc
    if not isinstance(body, dict):
        raise HTTPException(
            status_code=400,
            detail=f"request body must be a JSON object, got "
                   f"{type(body).__name__}")
    return body


def _num(v):
    """float(v) unless it is missing or non-finite, in which case ``None``."""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError, OverflowError):
        # OverflowError was missing and is the one that mattered: a JSON
        # integer literal above sys.float_info.max is a perfectly good int, so
        # TypeError and ValueError never fire and float() raises instead --
        # straight out through validate_spec as a 500. Found while fixing the
        # same bug in _as_number; this call site is reached one line later, so
        # patching only that one moved the crash rather than removing it.
        return None
    return f if math.isfinite(f) else None


#: Drugs curated in molecular.breast_library, not molecular.library.GBM_DRUGS
#: (the two are kept disjoint by test_disjoint_from_gbm_library), that are
#: nonetheless valid GBM-context drugs -- calibration._compile special-cases
#: each of these the same way it special-cases mitoxantrone. Single source of
#: truth for every GBM-context drug lookup in this module (_known_drugs,
#: drug_rows) so a name can never appear in one and be rejected by the other.
_CROSS_LIBRARY_GBM_DRUGS = ("palbociclib",)


def _known_drugs() -> List[str]:
    """Every name a RunSpec's ``drug`` may take: the curated GBM library, the
    cross-library reuses, and any compound the operator compiled from their own
    SMILES (``app/server/compounds.py``). User compounds are included here rather
    than special-cased at submit time so ``POST /api/validate`` and the run
    pipeline agree about what exists."""
    from tmesim.molecular.library import GBM_DRUGS

    from .compounds import user_compound_names

    return sorted(set(GBM_DRUGS) | set(_CROSS_LIBRARY_GBM_DRUGS)
                  | set(user_compound_names()))


# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
#: The Studio's default exposure window, in seconds. Named rather than inlined
#: because it is ONE OF THREE different defaults shipping in this repo — the CLI
#: runs a fast 1,800 s window (``tmesim.cli._FAST_TOTAL_TIME``) and the library
#: default is 43,200 s (``core.config.PDEConfig.total_time``). Exposure, coverage
#: and kill all scale with the window, so numbers from two surfaces are not
#: comparable unless the window matches. They are not silently harmonised: each
#: surface states the window it used instead (the report's meta block; the
#: ``compare-drugs`` CLI command did the same from ``results/export_honesty.md``
#: until it was removed on 2026-08-21 with the rest of the ranking). This is 6 h,
#: the
#: window ``experiments/patient_report.py`` runs.
STUDIO_DEFAULT_TOTAL_TIME_S = 3600.0 * 6


def _default_dx_um(section: str, preferred: float = 12.0) -> float:
    """``preferred``, coarsened until THIS server will accept it.

    THE DEFECT THIS CLOSES, measured on the live deployment 2026-09-20.
    ``GET /api/defaults`` handed a first-time user ``dx_um = 12.0``. On a
    512 MB host, where ``TMESIM_MAX_GRID_CELLS=auto`` resolves to 100,000,
    ``POST /api/runs`` answered that very spec with 400:

        dx_um = 12 um on this section gives a grid of ~284,498 voxels, above
        the 100,000 ceiling this server is configured with ...

    So the hosted product's FIRST CLICK FAILED. The defaults route and the
    validator disagreed about what this server can run, and nothing coupled
    them -- ``default_spec`` was a literal and the ceiling is resolved from the
    environment at import.

    This is the residue of a good fix. The same POST used to OOM-kill the
    container in about a second (findings table, row 8); the ceiling turned a
    crash into a refusal, which is strictly better and is why the ceiling
    exists. But converting a crash into a refusal left the defaults still
    pointing at the thing that gets refused, and a refusal on the primary
    action is a product defect even when the refusal is correct.

    UNCONSTRAINED SERVERS DO NOT MOVE. With the shipped ceiling (40,000,000)
    every built-in section clears ``preferred`` by orders of magnitude, so this
    returns 12.0 unchanged and no existing installation sees a different
    default. It coarsens ONLY where the alternative is handing the user a spec
    the same server will reject.

    It is NOT a claim that the coarser grid is converged. 12 um is not
    converged either -- see :func:`default_spec` -- and every report carries
    that in its own limits block. This changes which un-converged resolution a
    constrained host offers first, nothing more.

    AND ON THE SMALLEST TIER IT CANNOT EVEN REACH THE LADDER. Measured for
    ``xenium_319`` (8087.8 x 5065.4 um) at a 100,000 ceiling: the smallest dx
    that fits is 20.24 um, so this returns 21. The FINEST rung this repo has
    run a convergence ladder on is 20 um, which needs ~102,433 voxels and does
    not fit. So on a 512 MB host there is NO dx that both satisfies the memory
    ceiling and sits on a rung anyone has studied -- the two constraints have
    no overlap on this section, by about 2.4%.

    Nothing here papers over that. ``validate_spec`` already answers the
    returned spec with "dx_um = 21 um is coarser than every rung this repo has
    run a convergence ladder on", the run proceeds carrying that warning, and
    the honest reading of a hosted demo is that it demonstrates the pipeline
    rather than producing a number anyone should quote. A ceiling of ~110,000
    would put 20 um in reach; whether the tier survives that is an unmeasured
    question and must not be assumed from this comment.
    """
    row = next((r for r in list_sections() if r["id"] == section), None)
    extent = (row or {}).get("extent_um")
    if not extent:
        return preferred
    try:
        ex, ey = float(extent[0]), float(extent[1])
    except (TypeError, ValueError, IndexError):
        return preferred
    if ex <= 0 or ey <= 0:
        return preferred

    # The validator's own arithmetic: n = (ex/dx) * (ey/dx). Inverted, the
    # smallest dx that fits is sqrt(ex*ey / ceiling). Derived from the same
    # expression rather than guessed, so the two cannot disagree about the
    # boundary -- which is the whole failure being fixed.
    if (ex / preferred) * (ey / preferred) <= _MAX_GRID_CELLS:
        return preferred

    import math
    dx = math.sqrt((ex * ey) / float(_MAX_GRID_CELLS))
    # Up to a whole micrometre: a defaults panel reading "20.24 um" implies a
    # precision this number does not have, and rounding DOWN would re-cross
    # the ceiling it was computed from.
    dx = float(math.ceil(dx))
    # Then confirm against the real check, stepping coarser if the rounding or
    # any future rule still refuses. Bounded, and it falls back to `preferred`
    # rather than looping or inventing an absurd grid: a defaults route that
    # cannot find a good answer must not become a defaults route that hangs.
    for _ in range(64):
        if (ex / dx) * (ey / dx) <= _MAX_GRID_CELLS:
            return dx
        dx += 1.0
    return preferred


def default_spec() -> Dict:
    """A valid ``RunSpec`` carrying the current production defaults.

    Sourced from the production entry points, not retyped: ``dx = 12 um`` and a
    6 h window from ``experiments/patient_report.py``, the ``'uniform'``
    sensitivity mode it ranks in, ``reconstruct.real_fields``'s default
    ``'discrete'`` vascular mode, ``real_pipeline.OXYGEN_MODE_DEFAULT``, and the
    ``core.PKConfig`` field defaults.

    ``dx = 12 um`` is what the production entry point runs. It is NOT a
    converged resolution — an earlier version of this docstring called it "the
    converged resolution per results/numerical_rigor.md section 2", and that
    claim is withdrawn by the correction block at the top of that file and by
    ``results/grid_convergence_breast.md`` / ``results/grid_convergence_gbm.md``.
    Every report says so in its own limits block.
    """
    from tmesim.core import PKConfig
    from tmesim.reconstruct.real_pipeline import (
        OXYGEN_MODE_DEFAULT,
        STATE_MODE_DEFAULT,
    )

    pk = PKConfig()
    section = "xenium_319"
    if not any(s["id"] == section and s["available"] for s in list_sections()):
        first = next((s["id"] for s in list_sections() if s["available"]), None)
        if first:
            section = first
    return {
        "section": section,
        "dx_um": _default_dx_um(section),
        "drug": "temozolomide",
        "sensitivity_mode": "uniform",
        "vascular_mode": "discrete",
        "oxygen_mode": OXYGEN_MODE_DEFAULT,
        "state_mode": STATE_MODE_DEFAULT,
        "mgmt": None,
        "pk": {
            "plasma_peak": float(pk.plasma_peak),
            "half_life_s": float(pk.half_life),
            "n_doses": int(pk.n_doses),
            "dose_interval_s": float(pk.dose_interval),
        },
        "total_time_s": STUDIO_DEFAULT_TOTAL_TIME_S,
    }


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #
_SPEC_KEYS = {
    "section", "dx_um", "drug", "sensitivity_mode", "vascular_mode",
    "oxygen_mode", "state_mode", "mgmt", "pk", "total_time_s",
    "resistance_off",
    # --- mechanism attribution + time-lapse (all optional, all defaulted) --- #
    "attribute_mechanism", "timelapse_frames", "time_resolved",
}
_PK_KEYS = {"plasma_peak", "half_life_s", "n_doses", "dose_interval_s"}


def _is_whole(value) -> bool:
    """Is this number an exact integer, without ever calling ``float()`` on a
    value that cannot survive it?

    ``float(x) != int(x)`` is the natural spelling and it raises OverflowError
    on a JSON integer literal above ``sys.float_info.max`` -- a 500 out of a
    validator. An ``int`` is whole by definition and needs no conversion at
    all; only a genuine float has a fractional part to check.
    """
    if isinstance(value, int):
        return True
    return float(value).is_integer()


def _as_number(value, name: str, errors: List[str]) -> Optional[float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append(f"{name}: expected a number, got {type(value).__name__} "
                      f"({value!r})")
        return None
    try:
        v = float(value)
    except OverflowError:
        # A JSON integer literal above sys.float_info.max. Python parses it
        # exactly, as an int, so the isinstance check above passes and float()
        # is where it fails. MEASURED: 309 nines is a 500 out of
        # POST /api/validate; 308 is a clean 200. Reproduced on dx_um,
        # total_time_s, mgmt and three pk fields, and through /api/runs and
        # /api/reports' spec_overrides.
        #
        # This is the failure the docstring below already argues against -- a
        # 500 out of the endpoint whose whole job is turning bad input into a
        # readable message.
        #
        # bit_length(), NOT the value: str() of an int over ~4,300 digits
        # raises ValueError on Python 3.11+, so repr-ing it into the message
        # would re-create the 500 inside the handler for it.
        errors.append(f"{name}: integer is too large to be a floating-point "
                      f"number ({value.bit_length()} bits)")
        return None
    if not math.isfinite(v):
        errors.append(f"{name}: must be finite, got {value!r}")
        return None
    return v


def validate_spec(spec: Dict) -> Dict:
    """``POST /api/validate``: ``{ok, errors[], warnings[]}``.

    Every message names the offending field and says what to do about it,
    because the UI's code editor shows these verbatim. Errors mean the run
    would fail or is meaningless; warnings mean it would run but the result
    needs reading with care.
    """
    errors: List[str] = []
    warnings: List[str] = []

    if not isinstance(spec, dict):
        return {"ok": False,
                "errors": [f"RunSpec must be a JSON object, got "
                           f"{type(spec).__name__}"],
                "warnings": []}

    unknown = sorted(set(spec) - _SPEC_KEYS)
    if unknown:
        warnings.append(
            f"unrecognised field(s) {unknown}: ignored by the backend. Valid "
            f"RunSpec fields are {sorted(_SPEC_KEYS)}."
        )

    # --- section ---------------------------------------------------------- #
    section_id = spec.get("section")
    section_row = None
    if not isinstance(section_id, str) or not section_id:
        errors.append("section: required, must be a section id string. "
                      f"Known ids: {sorted(all_sections())}")
    elif section_id not in all_sections():
        errors.append(f"section: unknown id {section_id!r}. Known ids: "
                      f"{sorted(all_sections())}")
    else:
        section_row = next(
            (s for s in list_sections() if s["id"] == section_id), None)
        if section_row is not None and not section_row["available"]:
            errors.append(
                f"section {section_id!r} is not loadable: "
                f"{section_row.get('reason')}"
            )

    # Only ever look a section up by a STRING id. The code editor can POST
    # "section": ["x"] or {"a": 1}, and ``lookup_section`` does a dict ``get``,
    # which raises ``TypeError: unhashable type`` — a 500 out of a validation
    # endpoint whose whole job is to turn bad input into a readable message. The
    # error naming the bad type is already recorded above.
    if not isinstance(section_id, str) or not section_id:
        section_id = None
    source = (lookup_section(section_id) or {}).get("source") if section_id else None

    # --- modes ------------------------------------------------------------- #
    sens_mode = spec.get("sensitivity_mode", "uniform")
    if sens_mode not in SENSITIVITY_MODES:
        errors.append(f"sensitivity_mode: {sens_mode!r} is not one of "
                      f"{list(SENSITIVITY_MODES)}")
    elif sens_mode == "determinant" and source == "visium":
        warnings.append(
            "sensitivity_mode='determinant' has no effect on a Visium section: "
            "spot-level data carries no per-cell expression, so sensitivity "
            "stays uniform 1.0 and the kill remains a no-resistance upper bound."
        )
    elif sens_mode == "determinant":
        warnings.append(
            "sensitivity_mode='determinant' scales a compound down when its "
            "determinant is on the panel and leaves it at 1.0 when the "
            "determinant is absent, so the kill it produces is conditional on "
            "measurement COVERAGE rather than on biology. Use 'uniform' unless "
            "you specifically want that conditioning, and read the result as an "
            "upper bound on this one map either way. (This warning used to end "
            "'use uniform for any cross-drug ranking'; compound ordering was "
            "withdrawn on 2026-08-21 - results/pk_ranking_robustness.md, "
            "results/capacity_ranking_robustness.md.)"
        )

    vasc_mode = spec.get("vascular_mode", "discrete")
    if vasc_mode not in VASCULAR_MODES:
        errors.append(f"vascular_mode: {vasc_mode!r} is not one of "
                      f"{list(VASCULAR_MODES)}")
    elif vasc_mode == "kde" and source == "visium":
        warnings.append(
            "vascular_mode is ignored on a Visium section: its vascular field "
            "is built from the measured ENDOTHELIAL spot fraction and the "
            "reconstructor exposes no mode switch."
        )
    elif vasc_mode == "kde":
        warnings.append(
            "vascular_mode='kde' is the historical smooth 40 um endothelial "
            "KDE, nonzero over ~99.7% of the tissue. The production default is "
            "'discrete' (results/vessel_promotion.md); the KDE field is also "
            "re-derived AFTER the production reconstruction, which costs one "
            "extra oxygen solve."
        )

    oxy_mode = spec.get("oxygen_mode", "pde")
    if oxy_mode not in OXYGEN_MODES:
        errors.append(f"oxygen_mode: {oxy_mode!r} is not one of "
                      f"{list(OXYGEN_MODES)}")
    elif oxy_mode == "markers":
        if source == "visium":
            errors.append(
                "oxygen_mode='markers' needs per-cell expression and is Xenium "
                "only; a Visium section carries no per-cell transcripts."
            )
        elif section_id:
            from tmesim.physics.oxygen_proxy import HYPOXIA_MARKER_GENES

            genes = set(section_genes(section_id))
            if genes:
                present = [g for g in HYPOXIA_MARKER_GENES if g in genes]
                if not present:
                    errors.append(
                        f"oxygen_mode='markers' is unavailable on section "
                        f"{section_id!r}: its panel carries none of the "
                        f"hypoxia-response transcripts "
                        f"{list(HYPOXIA_MARKER_GENES)}. On this section hypoxia "
                        "is UNMEASURED — use 'pde' or 'vessel_distance'."
                    )
                else:
                    warnings.append(
                        f"oxygen_mode='markers' will use {present} — a hypoxia-"
                        "RESPONSE readout, rank-normalised, NOT an oxygen "
                        "tension. Also note sim.engine re-solves the oxygen PDE "
                        "internally, so this field does not drive the kill."
                    )
    elif oxy_mode == "vessel_distance":
        warnings.append(
            "oxygen_mode='vessel_distance' is a geometric SURROGATE, not a "
            "solved oxygen tension, and sim.engine re-solves the PDE "
            "internally, so it does not change this run's kill fraction."
        )

    _validate_state_mode(spec, section_id, source, errors, warnings)

    # --- dx ---------------------------------------------------------------- #
    dx = None
    if "dx_um" not in spec:
        errors.append("dx_um: required (grid spacing in micrometres)")
    else:
        dx = _as_number(spec["dx_um"], "dx_um", errors)
        if dx is not None:
            if dx <= 0:
                errors.append(f"dx_um: must be > 0, got {dx}")
                dx = None
            elif dx < MIN_DX_UM:
                # A REFUSAL, not a warning: memory scales as 1/dx^2 and a
                # warning does not block the run (see _MAX_GRID_CELLS).
                errors.append(
                    f"dx_um: must be >= {MIN_DX_UM:g} um, got {dx:g}. Grid "
                    "memory scales as 1/dx^2, and no tissue-scale question "
                    "needs a mesh finer than half a micrometre.")
                dx = None
            elif dx > DX_KILL_LADDER_COARSE_UM:
                warnings.append(
                    f"dx_um = {dx:g} um is coarser than every rung this repo "
                    f"has run a convergence ladder on "
                    f"({DX_KILL_LADDER_COARSE_UM:g} um and finer). NO dx is "
                    "grid-converged here — that claim is withdrawn (see the "
                    "correction in results/numerical_rigor.md, and "
                    "results/grid_convergence_breast.md / "
                    "results/grid_convergence_gbm.md, which contradict it) — "
                    "but the coarse rungs are measurably the worst: on both "
                    "published ladders the 24 and 30 um rungs understate the "
                    "scored ratio by up to 11.7% in a single refinement step."
                )
            elif dx > DX_PENETRATION_LADDER_COARSE_UM:
                warnings.append(
                    f"dx_um = {dx:g} um: the exposure metrics still move here. "
                    "The measured penetration ladder RISES with refinement "
                    "(0.1110 at dx = 60 um to 0.3062 at dx = 10 um), so a "
                    "coarse mesh UNDER-states penetration. An earlier version "
                    "of this warning said 'kill_fraction is converged' and that "
                    "penetration is '15-25% high'; both are withdrawn."
                )
            # `dx` is cleared above when it is non-positive, so both of the
            # remaining checks have to re-test it: a dx of -1 is an ERROR, and
            # crashing on it here would turn a 400 into a 500.
            # THE REFUSAL THE GRID CEILING STRUCTURALLY CANNOT MAKE.
            # TMESIM_MAX_GRID_CELLS counts VOXELS; a Visium reconstruction
            # spends its memory on the spot matrix, the scRNA reference and a
            # signature densification, all of which are sized by the PANEL and
            # not by dx. Measured 2026-09-20: ~3.9-4.3 GB at dx=30 on a host
            # declaring 512 MB, and coarsening dx moves none of it. So a
            # visitor choosing a Visium section walked straight past the
            # ceiling and OOM-killed the container.
            #
            # Done HERE, at submit time, rather than left to raise inside the
            # job: this file's own rule is that a deliberate refusal is a 400
            # the user reads, not a started-then-failed run they have to open
            # a log to understand.
            #
            # It refuses only against a DECLARED limit (cgroup v2 -> v1 ->
            # TMESIM_HOST_MEMORY_MB). An ordinary workstation declares none
            # and is untouched.
            if source == "visium" and section_id:
                try:
                    from tmesim.reconstruct.deconvolve import (
                        check_visium_memory_budget,
                    )
                    check_visium_memory_budget(section_id)
                except MemoryError as exc:
                    errors.append(str(exc))
                except Exception:  # noqa: BLE001
                    # A budget check that cannot run must never be the reason
                    # a valid spec is refused. The job-time refusal still
                    # stands behind it.
                    pass

            if dx is not None and source == "visium" and dx < 50.0:
                warnings.append(
                    f"dx_um = {dx:g} um over-resolves a 100 um Visium spot "
                    "pitch; reconstruct_real falls back to dx = 100 um for "
                    "Visium below 50 um, so this value will NOT be used."
                )
            if (dx is not None and section_row
                    and section_row.get("extent_um") and dx > 0):
                ex, ey = section_row["extent_um"]
                n_cells_grid = (ex / dx) * (ey / dx)
                if n_cells_grid > _MAX_GRID_CELLS:
                    # An ERROR, not a warning. `ok = not errors`, and
                    # POST /api/runs only refuses when ok is False, so a
                    # warning here would have let the allocation proceed.
                    # Size the units to the grid. The original message always
                    # said GB, which reads as "0.0 GB" for the case an operator
                    # actually hits -- a lowered TMESIM_MAX_GRID_CELLS on a
                    # small instance, where the grid is megabytes and the
                    # refusal is about the WHOLE SOLVE not fitting, not about
                    # one array. A refusal that quotes an absurd number invites
                    # the reader to dismiss it.
                    field_b = n_cells_grid * 8
                    field = (f"{field_b / 1e9:,.1f} GB" if field_b >= 1e9
                             else f"{field_b / 1e6:,.1f} MB")
                    tuned = _MAX_GRID_CELLS != _DEFAULT_MAX_GRID_CELLS
                    errors.append(
                        f"dx_um = {dx:g} um on this section gives a grid of "
                        f"~{n_cells_grid:,.0f} voxels, above the "
                        f"{_MAX_GRID_CELLS:,} ceiling"
                        + (" this server is configured with "
                           "(TMESIM_MAX_GRID_CELLS)" if tuned else "")
                        + f" — one float64 field alone would be {field}, and a "
                        "solve allocates several plus the cell population. "
                        "Use a coarser dx_um."
                    )
                elif n_cells_grid > _LARGE_GRID_CELLS:
                    warnings.append(
                        f"dx_um = {dx:g} um on this section gives a grid of "
                        f"~{n_cells_grid:,.0f} voxels; the reconstruction and "
                        "each transport solve will take minutes."
                    )

    # --- total_time --------------------------------------------------------- #
    total_time = None
    if "total_time_s" not in spec:
        errors.append("total_time_s: required (exposure window in seconds)")
    else:
        total_time = _as_number(spec["total_time_s"], "total_time_s", errors)
        if total_time is not None:
            if total_time <= 0:
                errors.append(f"total_time_s: must be > 0, got {total_time}")
            elif total_time < PDE_DT_S:
                errors.append(
                    f"total_time_s = {total_time:g} s is shorter than the "
                    f"production solver step dt = {PDE_DT_S:g} s, so the "
                    "transport solve would take zero steps and every exposure "
                    f"metric would be 0. Use total_time_s >= {PDE_DT_S:g}."
                )
            elif total_time > MAX_TOTAL_TIME_S:
                errors.append(
                    f"total_time_s = {total_time:g} s is above the "
                    f"{MAX_TOTAL_TIME_S:g} s ceiling ({MAX_TOTAL_TIME_S / (365.25 * 86400):.0f} "
                    f"years). At dt = {PDE_DT_S:g} s that is "
                    f"{total_time / PDE_DT_S:.3g} timesteps; the window is a "
                    "compute request at this size, not an exposure window. "
                    "Use a window you would actually observe.")
            elif total_time > 3600.0 * 72:
                warnings.append(
                    f"total_time_s = {total_time:g} s is "
                    f"{total_time / 3600.0:.0f} h; at dt = {PDE_DT_S:g} s that "
                    f"is {total_time / PDE_DT_S:,.0f} implicit solves."
                )

    # --- advection CFL ------------------------------------------------------ #
    _validate_cfl(section_id, dx, vasc_mode, oxy_mode,
                  spec.get("state_mode", "oxygen"), errors, warnings)

    # --- drug + mgmt --------------------------------------------------------- #
    drug = spec.get("drug", None)
    if drug is not None:
        if not isinstance(drug, str):
            errors.append(f"drug: expected a drug name or null, got "
                          f"{type(drug).__name__}")
        else:
            known = _known_drugs()
            if drug.strip().lower() not in known:
                errors.append(
                    f"drug: unknown drug {drug!r}. Available: {known}")
            else:
                from .compounds import user_compound_names

                if drug.strip().lower() in set(user_compound_names()):
                    warnings.append(
                        f"drug {drug!r} is a SELF-COMPILED compound: every one "
                        "of its transport parameters except molecular weight "
                        "and TPSA is a prediction from the SMILES you supplied, "
                        "and its logP is an RDKit Crippen prediction that cannot "
                        "see stereochemistry. Its kill_fraction, mean_tumor_auc "
                        "and penetration_length_um are NOT defensible as "
                        "absolute numbers. DO NOT run it alongside other "
                        "compounds and read the ordering: the cross-compound "
                        "RANKING was measured and WITHDRAWN on 2026-08-21 "
                        "(results/pk_ranking_robustness.md, "
                        "results/capacity_ranking_robustness.md), and this "
                        "warning told you to do exactly that until "
                        "2026-08-22. What IS readable is this compound "
                        "against ITSELF on this one section: the RELATIVE "
                        "exposure between regions of the same map, because a "
                        "per-compound plasma factor is an exact linear "
                        "rescale of the whole field and cancels from a ratio "
                        "between points. On how far the structure-to-"
                        "transport step can be trusted at all: against 71 "
                        "measured compounds the shipping two-route model "
                        "reaches R^2 = 0.718 blind "
                        "(results/two_route_promotion.md), which is a large "
                        "improvement on the molecular-weight-only model it "
                        "replaced (R^2 = 0.022, and flat by construction at "
                        "matched mass) and is still only about a QUARTER of the "
                        "measured 29.7x spread, from ONE cell line. A "
                        "similar-mass series is separated weakly rather than "
                        "not at all."
                    )
    else:
        # PERFECTLY INVERTED until 2026-08-22: this sat INSIDE the
        # `if drug is not None:` branch above, so it fired on every run that
        # HAD a drug -- announcing "drug is null" as the second thing a
        # first-time user read, contradicting the RunSpec editor two panels
        # left -- and was silent in the one case it exists for. A warning that
        # is wrong whenever it appears and absent whenever it applies is worse
        # than no warning: it teaches the reader to discount the whole list.
        warnings.append(
            "drug is null: this is a reconstruction-only run. kill_fraction, "
            "mean_tumor_auc, therapeutic_coverage_1x and penetration_length_um "
            "will all be null with a reason, and no drug_auc field is produced."
        )

    # --- counterfactual arm --------------------------------------------------- #
    ro = spec.get("resistance_off", False)
    if not isinstance(ro, bool):
        errors.append(
            f"resistance_off: expected true/false, got {type(ro).__name__} "
            f"({ro!r}). It selects the counterfactual arm of the rank "
            "attribution — PDConfig(hypoxia_resistance=1, "
            "quiescent_resistance=1) — not a dose."
        )
    elif ro:
        if drug is None:
            warnings.append(
                "resistance_off=true with drug=null has no effect: there is no "
                "kill to switch resistance off from."
            )
        else:
            warnings.append(
                "resistance_off=true is a COUNTERFACTUAL: this run switches "
                "this tumour's hypoxic AND quiescent kill penalties off "
                "(PDConfig(hypoxia_resistance=1, quiescent_resistance=1)). Its "
                "kill_fraction is NOT a prediction for this tissue — it is the "
                "denominator of rho = kill / kill(resistance off) in the rank "
                "attribution. Do not quote it as a therapy result."
            )

    mgmt = spec.get("mgmt")
    if mgmt is not None:
        if isinstance(mgmt, str):
            if mgmt.strip().lower() not in ("methylated", "unmethylated"):
                errors.append(
                    f"mgmt: a string mgmt must be 'methylated' or "
                    f"'unmethylated', got {mgmt!r}. A number is read as "
                    "log2(CPM+1) expression."
                )
        else:
            v = _as_number(mgmt, "mgmt", errors)
            if v is not None and v < 0:
                errors.append(
                    f"mgmt: log2(CPM+1) expression cannot be negative, got {v}")
        if isinstance(drug, str):
            from tmesim.molecular.mgmt import MGMT_CALIBRATED_DRUGS, is_mgmt_calibrated

            if not is_mgmt_calibrated(drug):
                warnings.append(
                    f"mgmt was supplied but {drug!r} is not MGMT-calibrated "
                    f"(only {list(MGMT_CALIBRATED_DRUGS)} are). It will be "
                    "recorded as a no-op with a reason; sensitivity is unchanged."
                )

    # --- PK ------------------------------------------------------------------ #
    pk = spec.get("pk")
    if pk is None:
        errors.append("pk: required, an object with "
                      f"{sorted(_PK_KEYS)}")
    elif not isinstance(pk, dict):
        errors.append(f"pk: expected an object, got {type(pk).__name__}")
    else:
        unknown_pk = sorted(set(pk) - _PK_KEYS)
        if unknown_pk:
            warnings.append(f"pk: unrecognised key(s) {unknown_pk}: ignored. "
                            f"Valid keys are {sorted(_PK_KEYS)}.")
        peak = _as_number(pk.get("plasma_peak", 0.0), "pk.plasma_peak", errors) \
            if "plasma_peak" in pk else None
        if "plasma_peak" not in pk:
            warnings.append("pk.plasma_peak not set; the core.PKConfig default "
                            "will be used.")
        elif peak is not None:
            if peak <= 0:
                errors.append(
                    f"pk.plasma_peak: must be > 0 uM, got {peak}. A zero or "
                    "negative peak delivers no drug and every exposure metric "
                    "would be 0.")
            elif peak > 1e4:
                warnings.append(
                    f"pk.plasma_peak = {peak:g} uM is far above any clinically "
                    "attainable plasma level for these agents (TMZ peaks near "
                    "50 uM); the result is a hypothetical, not a dose.")

        if "half_life_s" in pk:
            hl = _as_number(pk["half_life_s"], "pk.half_life_s", errors)
            if hl is not None:
                if hl <= 0:
                    errors.append(
                        f"pk.half_life_s: must be > 0 s, got {hl}; "
                        "plasma_profile divides by it (lambda = ln2/half_life).")
                elif hl > 3600.0 * 24 * 30:
                    warnings.append(
                        f"pk.half_life_s = {hl:g} s ({hl / 86400.0:.1f} days) "
                        "makes elimination negligible over the window.")

        if "n_doses" in pk:
            nd = pk["n_doses"]
            if isinstance(nd, bool) or not isinstance(nd, (int, float)):
                errors.append(f"pk.n_doses: expected an integer, got "
                              f"{type(nd).__name__}")
            elif not _is_whole(nd):
                errors.append(f"pk.n_doses: must be a whole number, got {nd}")
            elif int(nd) < 1:
                errors.append(f"pk.n_doses: must be >= 1, got {int(nd)}")
            elif int(nd) > MAX_N_DOSES:
                errors.append(
                    f"pk.n_doses = {int(nd)}: above the {MAX_N_DOSES:,} "
                    "ceiling. plasma_profile allocates an "
                    "(n_doses x n_times) float64 array, so this is a memory "
                    "request rather than a schedule -- a full Stupp course is "
                    "30-42 administrations.")
            elif int(nd) > 100:
                warnings.append(
                    f"pk.n_doses = {int(nd)}: doses superpose linearly, so a "
                    "large count raises the effective plasma level well above "
                    "plasma_peak.")

        n_doses = _num(pk.get("n_doses", 1))
        if "dose_interval_s" in pk:
            di = _as_number(pk["dose_interval_s"], "pk.dose_interval_s", errors)
            if di is not None and di <= 0 and (n_doses or 1) > 1:
                errors.append(
                    f"pk.dose_interval_s: must be > 0 s when n_doses > 1, got "
                    f"{di}; all doses would land at t = 0.")

        if (peak is not None and total_time is not None
                and n_doses is not None and n_doses > 1):
            di = _num(pk.get("dose_interval_s"))
            if di and di >= total_time:
                warnings.append(
                    f"pk.dose_interval_s = {di:g} s >= total_time_s = "
                    f"{total_time:g} s, so only the first dose falls inside the "
                    "simulated window.")

    _validate_readouts(spec, section_row, dx, total_time, errors, warnings)

    return {"ok": not errors, "errors": errors, "warnings": warnings}


def _validate_readouts(spec: Dict, section_row, dx, total_time,
                       errors: List[str], warnings: List[str]) -> None:
    """Validate the mechanism-attribution and time-lapse knobs.

    All three fields are optional and defaulted, so absence is never an error.
    What IS worth saying out loud: a time-lapse on a static plasma source cannot
    show the drug clearing, and the payload is ``n_frames`` full ``f32`` grids,
    which on a fine mesh is tens of megabytes.
    """
    from .jobs import DEFAULT_TIMELAPSE_FRAMES, TIMELAPSE_BUDGET_BYTES

    if "attribute_mechanism" in spec and not isinstance(
            spec["attribute_mechanism"], bool):
        errors.append(
            f"attribute_mechanism: expected true/false, got "
            f"{type(spec['attribute_mechanism']).__name__}")
    if spec.get("attribute_mechanism") is False and spec.get("drug"):
        warnings.append(
            "attribute_mechanism=false: this run will report kill_fraction "
            "without saying WHY any cell lived or died. The attribution is a "
            "read-out of arrays the kill model already computes and leaves the "
            "run bit-identical, so there is rarely a reason to switch it off.")

    time_resolved = spec.get("time_resolved")
    if time_resolved is not None and not isinstance(time_resolved, bool):
        errors.append(f"time_resolved: expected true/false, got "
                      f"{type(time_resolved).__name__}")
        time_resolved = None
    if time_resolved:
        warnings.append(
            "time_resolved=true drives the transport source with the real "
            "plasma(t) curve instead of its window average. Both paths deliver "
            "the same TOTAL drug, so a difference against a stored result is a "
            "SCHEDULE effect — but every calibration and regression in this "
            "repo was produced on the static path, so this run's absolute "
            "numbers are not directly comparable with them "
            "(results/dose_schedule.md).")

    frames = spec.get("timelapse_frames")
    if frames is None:
        frames = DEFAULT_TIMELAPSE_FRAMES
    elif isinstance(frames, bool) or not isinstance(frames, (int, float)):
        errors.append(f"timelapse_frames: expected an integer, got "
                      f"{type(spec['timelapse_frames']).__name__}")
        return
    elif not _is_whole(frames):
        errors.append(f"timelapse_frames: must be a whole number, got {frames}")
        return
    else:
        frames = int(frames)
        if frames < 0:
            errors.append(f"timelapse_frames: must be >= 0 (0 disables the "
                          f"time-lapse), got {frames}")
            return
        if frames == 1:
            errors.append(
                "timelapse_frames = 1: a time-lapse needs at least a first and "
                "a last frame. Use 0 to disable it, or >= 2.")
            return
        if frames > MAX_TIMELAPSE_FRAMES:
            warnings.append(
                f"timelapse_frames = {frames} exceeds the "
                f"MAX_TIMELAPSE_FRAMES = {MAX_TIMELAPSE_FRAMES} ceiling and "
                f"will be capped; the run reports the count it actually used.")
            frames = MAX_TIMELAPSE_FRAMES

    if frames >= 2 and spec.get("drug") is None:
        warnings.append(
            "timelapse_frames is set but drug is null: a reconstruction-only "
            "run solves no transport, so there is no drug field to film and no "
            "time-lapse payload will be produced.")
        return

    if frames >= 2 and not time_resolved:
        warnings.append(
            "a time-lapse on the default static plasma source shows the drug "
            "ENTER and SPREAD toward a quasi-steady state, but never CLEAR: the "
            "delivery rate is held at its window average for the whole window "
            "by construction. Set time_resolved=true to watch concentration "
            "fall.")

    if (frames >= 2 and section_row and section_row.get("extent_um")
            and dx and dx > 0):
        ex, ey = section_row["extent_um"]
        rows, cols = round(ey / dx), round(ex / dx)
        per_frame = rows * cols * 4
        total = per_frame * frames
        if total > TIMELAPSE_BUDGET_BYTES:
            warnings.append(
                f"timelapse_frames = {frames} at dx = {dx:g} um is "
                f"~{rows}x{cols} = {per_frame / 1e6:.2f} MB/frame = "
                f"{total / 1e6:.0f} MB, over the "
                f"{TIMELAPSE_BUDGET_BYTES / 1e6:.0f} MB budget. The frame count "
                f"will be REDUCED to what fits and the run will say by how "
                f"much; nothing is truncated silently.")
        elif total > 8e6:
            warnings.append(
                f"time-lapse payload will be ~{total / 1e6:.0f} MB "
                f"({frames} x {rows}x{cols} f32 frames).")

    if frames >= 2 and total_time:
        warnings.append(
            f"time-lapse interval will be ~{total_time / (frames - 1):.0f} s "
            f"({frames} frames over {total_time:g} s). Capture points are "
            f"rounded onto the solver's own {PDE_DT_S:g} s step grid, so the "
            f"realised interval can differ by up to one step.")


def _validate_state_mode(spec: Dict, section_id, source,
                         errors: List[str], warnings: List[str]) -> None:
    """Validate ``state_mode`` — and refuse to let a panel produce a state it
    cannot measure.

    The rule comes from ``results/cell_state_impact.md``, and the point it turns
    on is that **the cell-state label is the largest single term in every absolute
    kill number this project reports**: §3.2 measures kill moving by up to 0.3252
    absolute (1.77x relative) purely with the choice of route, everything else
    held identical.

    * ``'markers'`` on a **Visium** section is an ERROR — spot-level data carries
      no per-cell transcripts and ``_install_states`` raises.
    * ``'markers'`` on a panel carrying **neither** proliferation nor
      hypoxia-response transcripts is an ERROR: nothing at all would be measured,
      so the route has no content on that section.
    * ``'markers'`` on a panel missing ONE axis is a WARNING naming the axis, the
      genes and the consequence, because that state is then assigned to NO cell
      and "UNMEASURED" must never be read as "computed and found to be zero".
      This is the 319-gene flagship case for hypoxia.
    * The **detection RATE** cannot be known without reading the count matrix, so
      it is never asserted here. What is stated instead is the floor the engine
      will apply and what it does when the panel misses it — the flagship panel's
      measured 1.42% is 13x below the published 19-40% Ki-67 range (§1.1), so the
      per-cell call is refused there and a neighbourhood ENRICHMENT rule is used
      whose LEVEL moves 8.79x with the window radius.
    """
    state_mode = spec.get("state_mode", "oxygen")
    if state_mode not in STATE_MODES:
        errors.append(f"state_mode: {state_mode!r} is not one of "
                      f"{list(STATE_MODES)}")
        return

    if state_mode != "oxygen":
        warnings.append(
            f"state_mode={state_mode!r} changes cells.state on the "
            "RECONSTRUCTION only. sim.engine.simulate_response re-bands every "
            "living cell from its own oxygen solve before applying the kill, so "
            "this setting moves the reported band occupancies and their "
            "provenance — NOT this run's kill_fraction. It is also part of the "
            "reconstruction cache key, so selecting it forces a rebuild."
        )

    support = state_route_support(section_id) if section_id else {}

    if state_mode == "markers":
        if source == "visium":
            errors.append(
                "state_mode='markers' needs per-cell expression and is Xenium "
                "only; a Visium section is spot-level and carries no per-cell "
                "transcripts, so reconstruct_real would raise."
            )
        elif support.get("genes_known"):
            prolif = support["proliferation_present"]
            hypox = support["hypoxia_present"]
            if not prolif and not hypox:
                errors.append(
                    f"state_mode='markers' is unavailable on section "
                    f"{section_id!r}: its panel carries none of the proliferation "
                    f"transcripts {support['proliferation_missing']} and none of "
                    f"the hypoxia-response transcripts "
                    f"{support['hypoxia_missing']}. NOTHING would be measured on "
                    "this section, so the route has no content here — use "
                    "'oxygen' or 'vessel_distance'."
                )
            else:
                if not hypox:
                    warnings.append(
                        f"state_mode='markers' on section {section_id!r}: the "
                        f"panel carries NO hypoxia-response transcript (none of "
                        f"{support['hypoxia_missing']}), so HYPOXIC is UNMEASURED "
                        "and is assigned to NO cell. That is NOT the same "
                        "statement as 'computed and found to be zero' and must "
                        "not be reported as if it were. Do not quote a hypoxic "
                        "fraction, or an absolute kill, from this section."
                    )
                if not support["cycle_specific_present"]:
                    warnings.append(
                        "state_mode='markers': this panel carries no MKI67 or "
                        "TOP2A, so the proliferation call would rest on "
                        "low-specificity markers only. PCNA is deliberately "
                        "excluded from the default call (long protein half-life, "
                        "transcript expressed in non-cycling and DNA-repairing "
                        "cells)."
                    )
                warnings.append(
                    "state_mode='markers' will make a PER-CELL call only if the "
                    "union detection rate of "
                    f"{support['cycle_specific_present'] or prolif} reaches the "
                    "5% floor; that rate is a property of the count matrix and "
                    "is NOT knowable before the run, so it is measured at "
                    "reconstruction time and returned in "
                    "provenance.cell_state.detail.coverage. Below the floor the "
                    "module falls back to a NEIGHBOURHOOD enrichment rule whose "
                    "ORDERING is measured but whose LEVEL is not — it moves "
                    "8.79x across the window radius. On the 319-gene flagship "
                    "panel the measured union is 1.42%, 13x below the published "
                    "19-40% Ki-67 range, so that fallback is what fires there "
                    "(results/cell_state_impact.md §1.1, §2.4)."
                )
                if hypox:
                    warnings.append(
                        f"state_mode='markers' hypoxia will use {hypox} — "
                        "hypoxia-RESPONSE transcripts, which have real "
                        "non-hypoxic expression. The HYPOXIC band is an "
                        "ENRICHMENT label whose SIZE is UNCALIBRATED; do not "
                        "quote its magnitude as a hypoxic fraction."
                    )
    elif state_mode == "vessel_distance":
        warnings.append(
            "state_mode='vessel_distance' bands a geometric SURROGATE "
            "(exp(-d/L) in distance to the nearest endothelial cell, L from "
            "Thomlinson & Gray, PMID 13304213) through the unmodified "
            "cells.agents.update_cell_states. Whatever the positions are, the "
            "STATE is not measured. On the three Xenium sections scored in "
            "results/cell_state_impact.md §2.1 it called 91-100% of tumour cells "
            "VIABLE, so its dx-stability was there partly the stability of a "
            "near-constant — that figure is a property of those sections, not a "
            "prediction for this one."
        )
        if source == "visium":
            warnings.append(
                "state_mode='vessel_distance' on a Visium section measures "
                "distance to SYNTHESISED endothelial pseudo-cells: "
                "reconstruct.real_visium places round(fraction * 12) cells per "
                "spot at the spot centroid with a seeded jitter. The distance is "
                "therefore to a deconvolved 100 um spot grid, NOT to a measured "
                "vessel, and the resulting band occupancies inherit the spot "
                "pitch. Xenium sections carry measured centroids; this one does "
                "not."
            )
        elif support.get("genes_known") and not support.get("endothelial_present"):
                errors.append(
                    "state_mode='vessel_distance' is unavailable on section "
                    f"{section_id!r}: its panel carries no endothelial marker, so "
                    "no cell can be typed ENDOTHELIAL and the distance transform "
                    "has no point to measure from. Every cell would be labelled "
                    "HYPOXIC by construction — an artefact of an empty input, not "
                    "a state. reconstruct_real raises with this reason."
                )

    # The conclusion of results/cell_state_impact.md §5, applied to whichever
    # route was chosen: a panel with no hypoxia-response gene cannot support an
    # absolute kill under ANY route.
    if (support.get("genes_known") and not support.get("hypoxia_present")
            and spec.get("drug") is not None):
        warnings.append(
            f"ABSOLUTE KILL SHOULD NOT BE QUOTED FROM SECTION {section_id!r} "
            "UNDER ANY state_mode. Its panel measures no hypoxia, so NO state "
            "route clears the bar here: the 'oxygen' default's VIABLE band moves "
            "from 0.00% to 10.04% purely with the tissue-mask convention and the "
            "mesh, 'markers' falls back to a neighbourhood rule whose level moves "
            "8.79x with the window radius, and 'vessel_distance' is a geometric "
            "proxy. Kill moves up to 0.3252 absolute (1.77x) across the three "
            "routes with everything else identical. WHAT SURVIVES IS THE MAP: "
            "within this one compound on this one section, the RELATIVE "
            "exposure between regions, because a per-compound plasma factor is "
            "an exact linear rescale of the whole field (1.41e-14 over a 5,200x "
            "span, results/pk_ranking_robustness.md §1.2) and so cancels from a "
            "ratio between points. The LEVEL does not survive, and neither does "
            "a comparison against another compound - that claim was withdrawn "
            "on 2026-08-21 (results/pk_ranking_robustness.md, "
            "results/capacity_ranking_robustness.md). "
            "(results/cell_state_impact.md §3.2, §3.3, §5)."
        )


def _validate_cfl(section_id, dx, vasc_mode, oxy_mode, state_mode,
                  errors: List[str], warnings: List[str]) -> None:
    """Check ``dt * max|v| / dx <= 1`` for the fixed production ``dt``.

    ``physics.diffusion.simulate_drug_transport`` integrates advection with an
    explicit upwind step and raises when the CFL is violated — tens of seconds
    into a run, after the whole reconstruction has been paid for. The peak Darcy
    speed is a property of a reconstruction, so it can only come from one that
    has actually been built. The cache records every measured speed on disk, so
    the evidence survives a restart.

    Evidence is used in three tiers, and the message always says which tier it
    came from:

    1. **This exact key** has been reconstructed -> the speed is this spec's own
       measurement. A violation is an ERROR.
    2. **The same section and the same ``vascular_mode``**, at some other
       ``dx``, has been reconstructed. ``solve_ifp`` reads the vascular field
       and the domain and never reads oxygen or ``cells.state``, so those keys
       share this spec's vascular source and differ only in mesh (and in the
       oxygen / state routes, neither of which the pressure solve touches). A
       violation is still an ERROR, labelled with the key it was measured on.
    3. **Only a different ``vascular_mode``**, or nothing at all -> a WARNING
       saying the CFL is unverified and why. No speed is ever assumed.
    """
    if dx is None or not section_id or lookup_section(section_id) is None:
        return

    key = make_key(section_id, dx, vasc_mode, oxy_mode, state_mode)
    speed = RECON_CACHE.known_speed(key)
    evidence = "measured on this exact reconstruction"
    if speed is None:
        same_vasc = RECON_CACHE.speeds_for_section(section_id, vasc_mode)
        if same_vasc:
            src = max(same_vasc, key=lambda k: same_vasc[k])
            speed = same_vasc[src]
            evidence = (
                f"measured on key [{src}] — the same section and the same "
                f"vascular_mode, at a different mesh; physics.pressure."
                f"solve_ifp reads the vascular field, not oxygen"
            )

    if speed is not None:
        cfl = speed * PDE_DT_S / dx
        if cfl > 1.0:
            errors.append(
                f"dx_um = {dx:g} um violates the advection CFL at the "
                f"production dt = {PDE_DT_S:g} s: dt*max|v|/dx = {cfl:.3f} > 1, "
                f"with peak Darcy speed {speed:.4g} um/s ({evidence}). "
                f"physics.diffusion.simulate_drug_transport will refuse the "
                f"solve. Use dx_um >= {speed * PDE_DT_S:.4g} um."
            )
        elif cfl > 0.8:
            warnings.append(
                f"advection CFL = {cfl:.3f} at dx = {dx:g} um and dt = "
                f"{PDE_DT_S:g} s — inside the limit but close to it "
                f"(peak Darcy speed {speed:.4g} um/s, {evidence})."
            )
        return

    other_modes = RECON_CACHE.speeds_for_section(section_id)
    detail = (
        f" The only speeds recorded for this section come from a different "
        f"vascular_mode ({sorted(other_modes)}), which drives a different "
        f"pressure field, so they are not used here."
        if other_modes else
        f" No reconstruction of section {section_id!r} has ever been built in "
        f"this checkout."
    )
    warnings.append(
        f"advection CFL is UNVERIFIED for this spec: the limit is "
        f"dt*max|v|/dx <= 1 at the production dt = {PDE_DT_S:g} s, and the peak "
        f"interstitial Darcy speed max|v| is a property of the reconstruction."
        f"{detail} It will be measured on the first run and checked exactly "
        f"from then on."
    )


# --------------------------------------------------------------------------- #
# Report requests
# --------------------------------------------------------------------------- #
_REPORT_KEYS = {"section", "drug", "spec_overrides", "label"}

#: RunSpec fields a report OWNS and will not let a caller override, with the
#: reason each one is fixed. ``drug`` is the report's single subject and is
#: selected by the top-level ``drug`` field; ``resistance_off`` is swept across
#: the two arms; ``sensitivity_mode`` is forced because an upper bound that says
#: so is defensible and an estimate is not.
_REPORT_FIXED = {
    "drug": ("a report is about ONE compound and the top-level 'drug' field "
             "selects it; spec_overrides.drug would let the two disagree"),
    "resistance_off": ("both arms are generated for that one compound — the "
                       "real run and the resistance-off counterfactual that "
                       "MEASURES rho for the RESISTANCE decade"),
    "sensitivity_mode": (
        f"forced to {REPORT_SENSITIVITY_MODE!r}: on a section carrying no "
        "measured resistance determinant for this compound, a no-resistance "
        "UPPER BOUND is the honest reading and an estimate is not. "
        "'determinant' scales a compound by whichever determinant the panel "
        "happens to carry and leaves it at 1.0 where the panel does not, so the "
        "number moves with measurement coverage rather than with pharmacology "
        "(results/real_drug_panel.md)."
    ),
}


def validate_report_request(body: Dict) -> Dict:
    """``POST /api/reports`` body -> ``{ok, errors[], warnings[], section,
    drug, base_spec}``.

    The body takes a single ``drug``. It took a ``drugs[]`` panel until
    2026-08-21, when the cross-compound ranking was withdrawn
    (``results/pk_ranking_robustness.md``,
    ``results/capacity_ranking_robustness.md``): one compound is now the ONLY
    legal shape, and a ``drugs[]`` array is rejected with that reason rather than
    silently truncated to its first element.

    ``base_spec`` is a full, already-``validate_spec``-clean RunSpec with
    ``drug=null``; the manager clones it for each of the two arms, which
    therefore share one reconstruction cache key by construction.
    """
    errors: List[str] = []
    warnings: List[str] = []
    if not isinstance(body, dict):
        return {"ok": False, "errors": [
            f"report request must be a JSON object, got {type(body).__name__}"],
            "warnings": [], "section": None, "drug": None, "base_spec": None}

    unknown = sorted(set(body) - _REPORT_KEYS)
    if unknown:
        warnings.append(
            f"unrecognised field(s) {unknown}: ignored. Valid fields are "
            f"{sorted(_REPORT_KEYS)}.")

    section = body.get("section")
    if not isinstance(section, str) or not section:
        errors.append("section: required, a section id string. Known ids: "
                      f"{sorted(all_sections())}")
    elif section not in all_sections():
        errors.append(f"section: unknown id {section!r}. Known ids: "
                      f"{sorted(all_sections())}")

    # ``drugs[]`` is REJECTED rather than truncated. A caller who sends a panel
    # is asking for a ranking, and this endpoint no longer produces one; quietly
    # keeping the first element would return a document that answers a different
    # question from the one that was asked, under the name of the one that was.
    drug = None
    raw = body.get("drug")
    if "drugs" in body:
        errors.append(
            "drugs[]: a report is about ONE compound and no longer takes a "
            "panel. The cross-compound ranking was withdrawn on 2026-08-21 — "
            "per-compound plasma exposure across the compounds the GBM panel "
            "ranked spans 113.1x while the order breaks at 1.53x "
            "(results/pk_ranking_robustness.md sec. 0), and the breast top-1 "
            "changes under this repo's own structural capacity rule "
            "(results/capacity_ranking_robustness.md sec. 0). Send "
            "{'drug': '<name>'}. Your list is not truncated to its first "
            "element, because that would answer a different question under the "
            "name of the one you asked.")
    if raw is None:
        if "drugs" not in body:
            errors.append(
                "drug: required, a single drug name string. Available: "
                f"{_known_drugs()}")
    elif not isinstance(raw, str) or not raw.strip():
        errors.append(f"drug: expected a single drug name string, got "
                      f"{type(raw).__name__} ({raw!r})")
    else:
        known = _known_drugs()
        name = raw.strip().lower()
        if name not in known:
            errors.append(f"drug: unknown drug {raw!r}. Available: {known}")
        else:
            drug = name

    overrides = body.get("spec_overrides") or {}
    if not isinstance(overrides, dict):
        errors.append(f"spec_overrides: expected an object, got "
                      f"{type(overrides).__name__}")
        overrides = {}

    base = default_spec()
    if isinstance(section, str):
        base["section"] = section
    for k, v in overrides.items():
        if k in _REPORT_FIXED:
            warnings.append(f"spec_overrides.{k} was IGNORED: {_REPORT_FIXED[k]}")
            continue
        base[k] = v
    base["sensitivity_mode"] = REPORT_SENSITIVITY_MODE
    base["drug"] = None
    base["resistance_off"] = False

    # Validate the spec the arms will actually run, with the real drug in it, so
    # a report can never enqueue two runs that /api/validate would have refused.
    if drug and not errors:
        probe = dict(base, drug=drug)
        verdict = validate_spec(probe)
        errors.extend(verdict["errors"])
        for w in verdict["warnings"]:
            if w not in warnings:
                warnings.append(w)

    return {"ok": not errors, "errors": errors, "warnings": warnings,
            "section": section, "drug": drug, "base_spec": base,
            "map_claim": MAP_CLAIM}


# --------------------------------------------------------------------------- #
# Section caveats
#
# A section record says what a specimen IS. This says what it CANNOT support,
# in the one place every picker in the UI already reads. Two sources feed it:
#
#   1. The PLATFORM, which is knowable for every section, built-in or ingested.
#      Spot transcriptomics is the load-bearing case: reconstruct.real_visium
#      does not measure where the cells inside a 55 um spot sit, it SYNTHESISES
#      12 pseudo-cells per spot at the centroid with uniform jitter. Any
#      per-cell spatial statistic on such a section measures that jitter rather
#      than the tissue, so results/immune_sink_impact.md §5.6 declares it
#      INVALID and the generator refuses to compute it. `spatial_spread_valid`
#      carries that refusal to the UI so the Studio cannot offer the number.
#   2. The stored CAPABILITY report, which only an ingested section carries
#      (a user's panel is the unknown one). Its NO / CAVEAT verdicts are
#      surfaced verbatim - they are already written for a reader.
#
# Additive only: no existing field is changed or removed.
# --------------------------------------------------------------------------- #
_SPOT_SOURCES = ("visium",)

#: Stated once here so the API, the picker and the report cannot drift apart.
SPOT_SPATIAL_REFUSAL = (
    "Per-cell spatial statistics are INVALID on this section. A Visium spot "
    "pools several cells over ~55 um and this pipeline does not measure where "
    "they sit - it synthesises 12 pseudo-cells per spot at the spot centroid "
    "with uniform jitter. The per-cell attenuation spread and the "
    "immune-neighbourhood structure ratio would measure that jitter, not the "
    "tissue (results/immune_sink_impact.md §5.6). THE DELIVERY MAP IS INVALID "
    "HERE TOO - BUT NOT BY THE ROUTE THIS SENTENCE USED TO NAME. The falloff "
    "curve bins by distance to the nearest endothelial cell "
    "(validation.vessel_profile.direct_vessel_distance_um), the perivascular "
    "reference is the first bin of that curve, and the delivery spread is a "
    "percentile ratio over the exposure field - and all three move under a "
    "pure change of the placement seed. MEASURED 2026-09-20 on all four "
    "built-in Visium sections: a reseed alone moves penetration_length_um by "
    "up to 19.4% (visium_gbm3: 301.2 -> 331.8 -> 242.7 um). A quantity that "
    "moves that far when nothing about the tissue changed is a statement "
    "about the synthesis, not about this tissue. Read NOTHING here as a "
    "statement about where its vessels are. "
    "THE MECHANISM, corrected 2026-09-20. This paragraph said the exposure "
    "field's SOURCE is a density over synthesised ENDOTHELIAL CELL POSITIONS. "
    "On a spot section that is FALSE: reconstruct.real_fields builds the "
    "vascular field as a KDE over the MEASURED endothelial SPOT fraction, "
    "vascular_density_from_cells is called ZERO times on these sections "
    "(against once on the xenium_319 control), and vascular_density, "
    "ecm_density and tissue_mask are all BITWISE IDENTICAL across placement "
    "seeds. The real route is CONVECTION: physics.pressure rasterises TUMOUR "
    "cell positions into the interstitial pressure field, which advects the "
    "drug. With convection ON - the shipped default - 100% of tissue pixels "
    "move under a reseed (p99 1.76%, max 2.82%); with convection OFF the "
    "field is bitwise identical. The refusal stands; only its reason was "
    "wrong, and a refusal defended by a false reason is one relaxation away "
    "from being dropped by someone who checks it. "
    "(An earlier version, written 2026-08-21, claimed the exposure map and "
    "its falloff DID survive this refusal because they are solved on the "
    "grid. That was wrong too: it confused where the PDE is discretised with "
    "where its source term comes from. The same sentence previously named the "
    "cross-drug RANKING as what survived, which was withdrawn the same day - "
    "results/pk_ranking_robustness.md, "
    "results/capacity_ranking_robustness.md.)"
)


def _with_caveats(section: Dict) -> Dict:
    """Add ``spatial_spread_valid`` + ``caveats`` to a section record.

    Never raises: a section that cannot be described is still listed, just
    without its caveats, because dropping it from the picker would hide a
    specimen rather than qualify it.
    """
    out = dict(section)
    caveats: List[Dict] = []
    try:
        platform = str(out.get("platform") or out.get("source") or "").lower()
        is_spot = platform in _SPOT_SOURCES
        # A section whose platform we cannot name gets None, not True: "we did
        # not establish this" is a different claim from "this is valid".
        out["spatial_spread_valid"] = (None if not platform else not is_spot)
        if is_spot:
            caveats.append({
                "capability": "spatial_spread",
                "status": "NO",
                "reason": SPOT_SPATIAL_REFUSAL,
            })
        for row in ((out.get("capability") or {}).get("headline") or []):
            if row.get("status") not in ("NO", "CAVEAT"):
                continue
            if row.get("capability") == "spatial_spread" and is_spot:
                continue  # already stated above, verbatim from one source
            caveats.append({
                "capability": row.get("capability"),
                "status": row.get("status"),
                "reason": row.get("reason"),
            })
        for issue in (out.get("detect_issues") or []):
            caveats.append({"capability": "detection", "status": "CAVEAT",
                            "reason": issue})
    except Exception:  # noqa: BLE001 - a caveat must never break the listing
        out.setdefault("spatial_spread_valid", None)
    out["caveats"] = caveats
    out["n_caveats_blocking"] = sum(1 for c in caveats if c["status"] == "NO")
    return out


# --------------------------------------------------------------------------- #
# App
# --------------------------------------------------------------------------- #
def create_app(manager: Optional[JobManager] = None,
               projects: Optional[ProjectStore] = None,
               reports: Optional[ReportManager] = None) -> FastAPI:
    """Build the FastAPI application (a factory, so tests get a fresh one)."""
    @asynccontextmanager
    async def lifespan(_app: FastAPI):
        """Compile the drug library off the request path.

        ``drug_rows()`` reads the GDSC cache and compiles ten molecules, which
        takes ~15 s the first time. Doing it in a background thread means the
        UI's first ``GET /api/drugs`` waits on the memo rather than on a cold
        read, and nothing else is delayed by it.
        """
        threading.Thread(target=_warm_drugs, name="tmesim-warm-drugs",
                         daemon=True).start()
        yield

    app = FastAPI(
        title="Vivum Studio API",
        version="1",
        description=__doc__,
        lifespan=lifespan,
    )
    # SECURITY: no cross-origin grant by default. The front end is served
    # same-origin from this very app (see _mount_web), so it needs none.
    # `allow_origins=["*"]` used to be here, and it meant any web page the
    # operator happened to visit could read AND delete their whole compound
    # library through 127.0.0.1 -- binding to loopback is no defence, because
    # a browser is a local process. Set TMESIM_STUDIO_ALLOWED_ORIGINS to a
    # comma-separated list of explicit origins (e.g. "http://localhost:5173")
    # only for the split dev-port workflow. "*" is rejected on purpose.
    _allowed_origins = [
        o.strip() for o in
        os.environ.get("TMESIM_STUDIO_ALLOWED_ORIGINS", "").split(",")
        if o.strip() and o.strip() != "*"
    ]
    if _allowed_origins:
        app.add_middleware(
            CORSMiddleware,
            allow_origins=_allowed_origins,
            allow_methods=["GET", "POST", "DELETE", "OPTIONS"],
            allow_headers=["content-type", "x-api-key"],
        )

    # SECURITY (defence in depth, NOT the fix for the Host-header auth bypass
    # -- that is fixed at the gate itself, below). Validating the Host header
    # blocks DNS rebinding: a hostile name that resolves to 127.0.0.1 makes a
    # hostile page same-origin with this server, which no CORS policy can
    # stop. Off unless TMESIM_STUDIO_ALLOWED_HOSTS is set; __main__.py sets a
    # loopback allowlist for the default local bind. Note starlette's
    # implementation splits the Host header on ":" and takes [0], so a
    # bracketed IPv6 literal cannot be allowlisted, and a value like
    # "127.0.0.1:8000/../" passes it -- which is exactly why it is not
    # load-bearing here.
    _allowed_hosts = [
        h.strip() for h in
        os.environ.get("TMESIM_STUDIO_ALLOWED_HOSTS", "").split(",")
        if h.strip()
    ]
    if _allowed_hosts and "*" not in _allowed_hosts:
        app.add_middleware(TrustedHostMiddleware, allowed_hosts=_allowed_hosts)

    class _CapRequestBodyASGI:
        """Enforce the body cap on the BYTES THAT ARRIVE, not on what the
        client says will arrive.

        WHY THIS EXISTS ALONGSIDE the Content-Length pre-check below. That
        check is honest about being a pre-check -- its own docstring says the
        header "is client-supplied and may be absent or a lie". A security
        pass on 2026-09-20 showed the lie is cheap to tell: the same
        multi-megabyte body that gets a 413 with a Content-Length sails
        through under ``Transfer-Encoding: chunked``, because no header
        declares a size. Only ``protocol_api`` used ``_capped_body``; the
        other ``await request.json()`` sites were uncapped in practice.

        So the cap is moved to where it cannot be lied to. This reads the body
        itself, stopping at ``MAX_REQUEST_BODY_BYTES + 1`` -- bounded by
        construction, so the defence cannot become the exhaustion it prevents
        -- and either refuses with 413 or REPLAYS the buffered bytes to the
        application so every downstream ``request.json()`` behaves exactly as
        before.

        Raw ASGI rather than ``@app.middleware("http")`` deliberately:
        BaseHTTPMiddleware cannot read the stream and still hand an intact
        stream onward, which is precisely the operation needed here.

        GET/HEAD/OPTIONS pass through untouched, for the reason the pre-check
        gives: refusing a body on a liveness probe is a new way to fail a
        health check, and this file has done that twice already.
        """

        def __init__(self, app, max_bytes: int):
            self.app = app
            self.max_bytes = max_bytes

        async def __call__(self, scope, receive, send):
            if (scope.get("type") != "http"
                    or scope.get("method") in ("GET", "HEAD", "OPTIONS")):
                return await self.app(scope, receive, send)

            chunks = []
            total = 0
            more = True
            while more:
                message = await receive()
                if message["type"] == "http.disconnect":
                    # Hand it on unchanged; the app decides what a disconnect
                    # means. Swallowing it here would hang the request.
                    chunks.append(message)
                    break
                body = message.get("body", b"")
                total += len(body)
                if total > self.max_bytes:
                    return await self._refuse(send, total)
                chunks.append(message)
                more = bool(message.get("more_body", False))

            sent = False

            async def replay():
                nonlocal sent
                if chunks:
                    return chunks.pop(0)
                # The stream is exhausted. Anything further is a disconnect --
                # returning a body message with more_body would deadlock a
                # handler waiting for the end.
                sent = True
                return {"type": "http.disconnect"}

            return await self.app(scope, replay, send)

        async def _refuse(self, send, seen: int):
            payload = json.dumps({"detail": {
                "message": (f"request body exceeds the "
                            f"{self.max_bytes} B limit"),
                "reason": ("a RunSpec, a report request and a compound "
                           "definition are all kilobytes. This was measured "
                           "on the bytes received, so declaring no length "
                           "does not raise the limit."),
            }}).encode("utf-8")
            # The size SEEN is deliberately not echoed: the read stops at the
            # cap, so it is not the body's real size, and printing a number
            # that is merely "where we stopped" invites the reader to treat it
            # as a measurement.
            await send({"type": "http.response.start", "status": 413,
                        "headers": [(b"content-type", b"application/json"),
                                    (b"content-length",
                                     str(len(payload)).encode())]})
            await send({"type": "http.response.body", "body": payload})

    app.add_middleware(_CapRequestBodyASGI, max_bytes=MAX_REQUEST_BODY_BYTES)

    @app.middleware("http")
    async def _cap_request_body(request: Request, call_next):
        """Refuse an oversized body BEFORE a handler buffers it.

        Content-Length first, because it lets us refuse without reading a byte.
        It is client-supplied and may be absent or a lie, so this is a cheap
        pre-check and not the whole defence -- the protocol router additionally
        counts its stream as it arrives, and a route needing that guarantee
        should use `_capped_body`. What this stops is the ordinary case: a
        large declared body reaching `await request.json()`.

        GET/HEAD/OPTIONS are skipped rather than trusted. They can technically
        carry a body, and refusing one on a liveness probe would be a new way
        to fail a health check -- something this file has already done twice.
        """
        if request.method not in ("GET", "HEAD", "OPTIONS"):
            declared = request.headers.get("content-length")
            if (declared is not None and declared.isdigit()
                    and int(declared) > MAX_REQUEST_BODY_BYTES):
                return JSONResponse(
                    {"detail": {
                        "message": (f"request body is {int(declared)} B, above "
                                    f"the {MAX_REQUEST_BODY_BYTES} B limit"),
                        "reason": ("a RunSpec, a report request and a compound "
                                   "definition are all kilobytes"),
                    }},
                    status_code=413)
        return await call_next(request)

    # Optional shared-secret gate on every /api/* route, OFF by default (empty
    # env var) so the local single-user workflow this app was built for is
    # completely unaffected. __main__.py sets TMESIM_STUDIO_API_KEY when a
    # non-loopback --host is requested without an explicit --api-key, so a
    # remote deployment can never come up silently unauthenticated. Static
    # files (the frontend itself) are deliberately NOT gated, so the page can
    # load and prompt for the key rather than failing to render at all.
    api_key = os.environ.get("TMESIM_STUDIO_API_KEY", "").strip()
    if api_key:
        @app.middleware("http")
        async def _require_api_key(request: Request, call_next):
            # SECURITY: decide on the SAME path string the router dispatches
            # on -- starlette.routing.Route.matches() calls get_route_path()
            # too, so the gate and the dispatch physically cannot disagree.
            # `request.url.path` must never be used here: starlette builds
            # request.url by concatenating the client-supplied Host header
            # with the path and re-parsing, so `Host: x/../` (or
            # `Host: 127.0.0.1:8000/../`, which passes host validation)
            # yields url.path='/..//api/compounds' -- the prefix test fails,
            # the gate is skipped, and the router still runs the endpoint.
            # That is GHSA-86qp-5c8j-p5mr and it was live here.
            route_path = get_route_path(request.scope)
            # Fail closed: anything that is not a plain read of a static
            # asset needs the key, so a future non-/api route that mutates
            # state is authenticated by default rather than exposed by
            # default. OPTIONS stays exempt -- a CORS preflight cannot carry
            # the header.
            # /healthz is the ONE exempt path, and it is exempt because a
            # platform health checker cannot carry a key. Render probes the
            # configured path every ~10 s with no headers; pointed at
            # /api/health it got 401 forever, concluded the service never came
            # up, and killed a container that was serving correctly on the port
            # it had already detected. The log read "401 Unauthorized" a hundred
            # times and then "Timed Out" -- a healthy process, failed by its own
            # health check.
            #
            # It is safe to exempt because of what it returns: a literal
            # {"ok": true} and nothing else. NOT /api/health, which reports the
            # runs directory, run counts and cache statistics -- that stays
            # behind the key. Liveness is not a data route, and this one is
            # incapable of becoming one.
            # GET/HEAD only, so the fail-closed rule below is untouched: a POST
            # to /healthz still needs the key like any other state-changing
            # request.
            # /api/health is exempt TOO, and only because reality insisted.
            # /healthz exists and is the right path to probe -- but a platform
            # health check is configured in a DASHBOARD, and the one already
            # pointed at /api/health kept 401ing and killing a container that
            # was serving correctly. Telling the operator to go and change a
            # setting is not a fix when the code can simply stop being wrong.
            #
            # Safe because of what the endpoint DOES when unauthenticated: see
            # the handler, which returns a bare {"ok": true} to an anonymous
            # caller and the runs directory / counts / cache statistics only to
            # one holding the key. The route is exempt from the GATE; the DATA
            # is not.
            is_liveness = (route_path in ("/healthz", "/api/health")
                           and request.method in ("GET", "HEAD"))
            needs_key = not is_liveness and (
                route_path.startswith("/api/")
                or request.method not in ("GET", "HEAD", "OPTIONS")
            )
            if needs_key and request.method != "OPTIONS":
                # constant-time: a `!=` on a secret leaks its prefix by timing
                if not secrets.compare_digest(
                        request.headers.get("x-api-key", ""), api_key):
                    return JSONResponse(
                        {"detail": "missing or incorrect X-API-Key header"},
                        status_code=401,
                    )
            return await call_next(request)
    jobs = manager if manager is not None else JobManager()
    store = projects if projects is not None else PROJECT_STORE
    reports = reports if reports is not None else ReportManager(jobs)
    app.state.jobs = jobs
    app.state.cache = RECON_CACHE
    app.state.projects = store
    app.state.reports = reports

    # ----------------------------------------------------------------- #
    # HEAD AS WELL AS GET, and this is not tidiness.
    #
    # FastAPI registers only the method you name, so HEAD /healthz returned 404
    # — while the API-key middleware had already been written to exempt
    # `request.method in ("GET", "HEAD")` on this path, i.e. the gate expected
    # HEAD probes that the router then refused. A platform health check
    # configured for HEAD would have marked this service down while it was
    # serving perfectly.
    #
    # That is the failure this route exists because of: a health check pointed
    # at a route it could not satisfy killed a healthy container here twice.
    # Found by an input sweep, not by a probe, because nothing in this repo
    # probes with HEAD.
    @app.api_route("/healthz", methods=["GET", "HEAD"])
    def healthz() -> Dict:
        """Liveness, for a platform health checker. NOT ``/api/health``.

        Deliberately outside ``/api/`` and deliberately empty of information:
        a platform probe cannot send an ``X-API-Key`` header, so this route is
        exempt from the gate, and the only safe thing for an ungated route to
        say is that the process is answering.

        ``/api/health`` stays behind the key and keeps reporting the runs
        directory, run counts and cache statistics. Point a health check at
        THIS one; pointing it at that one returns 401 to a checker that cannot
        authenticate, and the platform then kills a container that was serving
        perfectly well.
        """
        return {"ok": True}

    @app.get("/api/sections")
    def get_sections() -> List[Dict]:
        return jsonable([_with_caveats(s) for s in list_sections()])

    @app.post("/api/sections/scan")
    async def post_scan(request: Request) -> Dict:
        """Inspect ANY directory: what is it, and what can it support?

        Body ``{path, measure_detection_rates=false, register=false,
        section_id=null, label=null}``.

        The default is a READ-ONLY scan that touches only filenames and
        ``features.tsv.gz`` — milliseconds even on the 851k-cell section — and
        reports every RATE-dependent verdict as UNMEASURED rather than guessing
        it. ``measure_detection_rates`` performs a full ingest so the report can
        state the measured proliferation and hypoxia detection rates, which is
        what turns "MKI67 is on the panel" into "MKI67 fires in 1.34% of cells,
        below the 5% floor". ``register`` adds the sample to the section
        registry (implies a full ingest, since a section is only registered once
        it has loaded and passed validation).

        Never 500s on a bad directory: an unrecognised or broken sample is an
        ANSWER, returned in ``detection.issues`` / ``error``.

        SECURITY: "ANY directory" is the literal truth and the API key does not
        change it — the key authenticates this route, it does not confine it,
        so a key handed to a field tester is read access to this host as this
        process's uid. Setting ``TMESIM_STUDIO_DATA_ROOTS`` confines ``path``
        to an allowlist of roots, checked against the RESOLVED path (``~``,
        ``..`` and symlinks all followed first) rather than against the text.
        Unset, scanning stays unconfined — the bring-your-own-tissue workflow
        points at wherever the operator's export lives — and the first scan in
        the process says so on stderr. See ``docs/DEPLOYMENT.md`` §4.

        A refusal is a 403 that names the caller's own path, the variable and
        how to add a root, and never the roots themselves. It does leak one
        bit per query about whether a given directory is inside the allowlist;
        that is the price of a diagnostic that tells an operator what to do,
        and this route is behind the API key whenever one is set.
        """
        try:
            body = await request.json()
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(
                status_code=400,
                detail=f"request body is not valid JSON: {exc}") from exc
        # isinstance(str), not truthiness. The test was `not body.get("path")`
        # and the next line was `str(body["path"])`, so `{"path": {"a": 1}}`
        # became the literal path `{'a': 1}`, was resolved against the server's
        # working directory, and came back in the 200 body -- turning a type
        # error into a disclosure of the server's CWD. No new read capability,
        # but no reason to answer a dict with a filesystem path either.
        if (not isinstance(body, dict)
                or not isinstance(body.get("path"), str)
                or not body["path"].strip()):
            raise HTTPException(
                status_code=400,
                detail="body must be an object with a 'path' string naming the "
                       "sample directory to scan")

        from . import ingest_registry

        path = body["path"]
        try:
            confined = confine_scan_path(path)
        except PathNotAllowed as exc:
            raise HTTPException(status_code=403, detail=exc.as_detail()) from None
        if confined is not None:
            # Confinement is in force: hand the ingester the path that was
            # actually checked, so the object scanned and the object allowed
            # cannot be two different directories.
            path = confined
        deep = bool(body.get("measure_detection_rates", False))
        if bool(body.get("register", False)):
            out = ingest_registry.register(
                path,
                section_id=body.get("section_id"),
                label=body.get("label"),
            )
            if out["ok"]:
                out["detection"] = {
                    k: v for k, v in out["section"].items()
                    if k in ("platform", "layout", "detect_confidence",
                             "detect_issues")
                }
            return jsonable(out)
        return jsonable(ingest_registry.scan(path, measure_detection_rates=deep))

    @app.delete("/api/sections/{section_id}")
    def delete_section(section_id: str) -> Dict:
        """Remove an INGESTED section. Built-in sections are never removable."""
        from . import ingest_registry

        out = ingest_registry.unregister(section_id)
        if not out["ok"]:
            raise HTTPException(status_code=400, detail=out["error"])
        return jsonable(out)

    @app.get("/ingest", response_class=HTMLResponse, include_in_schema=False)
    def ingest_page() -> str:
        """Server-rendered form over ``POST /api/sections/scan``.

        **Why this exists.** ``results/first_user_readiness.md`` found that
        bring-your-own-tissue works end to end through the API and the CLI but
        appears NOWHERE in ``app/web/`` — the Studio's section picker is a
        read-only list from ``GET /api/sections``. So a researcher who will not
        hand-POST JSON cannot load their own sample, which is the entire value
        proposition. This page closes that with no front-end change: it is one
        self-contained document served from the backend, and it calls the same
        endpoint ``curl`` would.

        It adds NO capability and NO new trust boundary — every request it makes
        is one ``POST /api/sections/scan`` that any client could already make,
        under the same ``X-API-Key`` gate. It renders ``capability_text``
        verbatim (the exact report ``vivum ingest`` prints) rather than
        re-deriving a summary, so this surface can never disagree with the CLI
        about what a sample can support.

        This is a stopgap, not the destination: the Studio picker itself should
        learn to ingest. Registered sections appear in ``GET /api/sections``
        immediately, so a Studio reload picks them up.
        """
        return _INGEST_PAGE

    @app.get("/api/drugs")
    def get_drugs() -> List[Dict]:
        return jsonable(drug_rows())

    # --- bring your own compound ------------------------------------------ #
    # A researcher's own series is the point: they want to rank THEIR molecules.
    # POST .../compile is a dry run (nothing is stored) so the form can show the
    # per-field provenance BEFORE the operator commits to the name.

    @app.post("/api/compounds/compile")
    async def post_compound_compile(request: Request) -> Dict:
        from .compounds import CompoundRefused, compile_preview

        body = await _json_body(request)
        try:
            return jsonable(compile_preview(
                body.get("name"), body.get("smiles"), body.get("ic50_uM")))
        except CompoundRefused as exc:
            raise HTTPException(status_code=422, detail=exc.as_dict())

    @app.get("/api/compounds")
    def get_compounds() -> List[Dict]:
        from .compounds import COMPOUND_STORE

        return jsonable(COMPOUND_STORE.list())

    @app.post("/api/compounds", status_code=201)
    async def post_compound(request: Request) -> Dict:
        from .compounds import COMPOUND_STORE, CompoundError, CompoundRefused

        body = await _json_body(request)
        try:
            return jsonable(COMPOUND_STORE.create(
                body.get("name"), body.get("smiles"), body.get("ic50_uM"),
                note=body.get("note")))
        except CompoundRefused as exc:
            raise HTTPException(status_code=422, detail=exc.as_dict())
        except CompoundError as exc:
            raise HTTPException(status_code=exc.status,
                                detail={"code": exc.code,
                                        "message": exc.message,
                                        "field": "name"})

    @app.get("/api/compounds/{compound_id}")
    def get_compound(compound_id: str) -> Dict:
        from .compounds import COMPOUND_STORE, CompoundError, CompoundRefused

        try:
            return jsonable(COMPOUND_STORE.compile_stored(compound_id))
        except CompoundError as exc:
            raise HTTPException(status_code=exc.status,
                                detail={"code": exc.code,
                                        "message": exc.message})
        except CompoundRefused as exc:
            raise HTTPException(status_code=422, detail=exc.as_dict())

    @app.delete("/api/compounds/{compound_id}")
    def delete_compound(compound_id: str) -> Dict:
        from .compounds import COMPOUND_STORE, CompoundError

        # CompoundError is CAUGHT here, unlike every sibling route above which
        # already did. The store now validates the id at its own boundary (see
        # compounds.py::_dir and the traversal it closed), so a malformed id
        # raises here instead of being joined onto a path -- and an uncaught
        # store error out of a DELETE is a 500 where the store meant a 404.
        try:
            deleted = COMPOUND_STORE.delete(compound_id)
        except CompoundError as exc:
            raise HTTPException(status_code=exc.status,
                                detail={"code": exc.code,
                                        "message": exc.message}) from exc
        if not deleted:
            raise HTTPException(
                status_code=404,
                detail={"code": "compound_not_found",
                        "message": f"no user compound {compound_id!r}"})
        return {"deleted": compound_id}

    @app.get("/api/defaults")
    def get_defaults() -> Dict:
        return jsonable(default_spec())

    @app.post("/api/validate")
    async def post_validate(request: Request) -> Dict:
        try:
            spec = await request.json()
        except Exception as exc:  # noqa: BLE001 - the editor sends broken JSON
            return {"ok": False,
                    "errors": [f"request body is not valid JSON: {exc}"],
                    "warnings": []}
        return jsonable(validate_spec(spec))

    @app.post("/api/runs", status_code=202)
    async def post_run(request: Request) -> Dict:
        try:
            spec = await request.json()
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(
                status_code=400,
                detail=f"request body is not valid JSON: {exc}") from exc
        verdict = validate_spec(spec)
        if not verdict["ok"]:
            raise HTTPException(status_code=400, detail={
                "message": "RunSpec rejected by /api/validate",
                **verdict,
            })
        job = jobs.submit(spec)
        return {"run_id": job.run_id, "status": job.status,
                "warnings": verdict["warnings"]}

    @app.get("/api/runs")
    def get_runs(limit: int = Query(200, ge=1, le=5000)) -> List[Dict]:
        return jsonable([j.result_view() for j in jobs.list(limit=limit)])

    @app.get("/api/runs/{run_id}")
    def get_run(run_id: str) -> Dict:
        job = _require(jobs, run_id)
        return jsonable(job.result_view())

    @app.get("/api/runs/{run_id}/status")
    def get_run_status(run_id: str,
                       log_lines: int = Query(40, ge=1, le=500)) -> Dict:
        job = _require(jobs, run_id)
        return jsonable(job.status_view(log_lines))

    @app.get("/api/runs/{run_id}/traceback")
    def get_run_traceback(run_id: str) -> Dict:
        """The failure stack for a ``failed`` run — opt-in, and path-scrubbed.

        SECURITY: ``traceback.format_exc()`` always carries absolute host
        paths, so serving it disclosed the repo root, the operator's OS
        username and the internal module layout to any client. The typed
        diagnosis (``job.error``) is what the product promises and is always
        returned; the stack itself needs
        ``TMESIM_STUDIO_EXPOSE_TRACEBACKS=1``, and even then the repo root is
        replaced with ``<repo>``. The unscrubbed stack is always on the
        server's stderr for whoever is at the console.
        """
        job = _require(jobs, run_id)
        out = {"run_id": run_id, "status": job.status, "error": job.error,
               "traceback": None}
        if os.environ.get("TMESIM_STUDIO_EXPOSE_TRACEBACKS", "") != "1":
            out["traceback_withheld"] = (
                "stack traces are not served by default because they contain "
                "absolute host paths. Set TMESIM_STUDIO_EXPOSE_TRACEBACKS=1 "
                "on the server to enable, or read the server's stderr. The "
                "typed error is in `error`.")
            return out
        tb = job.traceback
        if tb:
            tb = tb.replace(repo_root(), "<repo>")
        out["traceback"] = tb
        return out

    @app.delete("/api/runs/{run_id}")
    def delete_run(run_id: str) -> Dict:
        job = _require(jobs, run_id)
        if job.status in ("done", "failed", "cancelled"):
            return {
                "run_id": run_id, "status": job.status, "cancelled": False,
                "reason": f"run is already {job.status}; nothing to cancel",
            }
        jobs.cancel(run_id)
        return {
            "run_id": run_id, "status": job.status, "cancelled": True,
            "reason": (
                "cancel flag set. A queued run is cancelled immediately; a "
                "running run stops at the next stage boundary, so one already "
                "inside simulate_response finishes that solve first."
            ),
        }

    @app.get("/api/runs/{run_id}/geometry")
    def get_geometry(run_id: str):
        return _binary(jobs, run_id, "geometry", "cell geometry")

    @app.get("/api/runs/{run_id}/field/{name}")
    def get_field(run_id: str, name: str):
        if name not in FIELD_NAMES:
            raise HTTPException(
                status_code=404,
                detail=f"unknown field {name!r}; expected one of "
                       f"{list(FIELD_NAMES)}")
        return _binary(jobs, run_id, f"field_{name}", f"field {name!r}")

    @app.get("/api/runs/{run_id}/mechanism")
    def get_mechanism(run_id: str):
        """Per-cell mechanism attribution (TMSM).

        One ``uint8`` label per cell in the geometry payload's order, plus the
        aggregate breakdown — counts, the delivery/dose/resistance/chance split
        of the survivors, the thresholds that produced them and the bit-exact
        reconciliation against ``kill_fraction`` — carried in the same response
        as a short JSON block. See ``tmesim.sim.mechanism``.
        """
        return _binary(jobs, run_id, "mechanism", "mechanism attribution")

    @app.get("/api/runs/{run_id}/timelapse")
    def get_timelapse(run_id: str):
        """Drug-concentration time-lapse (TMSL).

        ``n_frames`` snapshots of the CONCENTRATION field (uM) through the
        exposure window, each a complete TMSF block prefixed by its frame index
        and time stamp. Frame count and interval are in the container header and
        repeated in ``RunResult.timelapse``.
        """
        return _binary(jobs, run_id, "timelapse", "drug time-lapse")

    # ------------------------------------------------------------------ #
    # Projects / specimens — the multi-patient layer.
    #
    # A project owns nothing but ids: specimens hold ``run_id``s into
    # ``app/runs/``, which stays the single store for run results and binary
    # payloads. Attaching, detaching and deleting therefore never touch a run.
    # ------------------------------------------------------------------ #
    def _run_view(run_id: str) -> Optional[Dict]:
        job = jobs.get(run_id)
        return jsonable(job.result_view()) if job is not None else None

    async def _body(request: Request) -> Dict:
        try:
            doc = await request.json()
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(
                status_code=400,
                detail=f"request body is not valid JSON: {exc}") from exc
        if not isinstance(doc, dict):
            raise HTTPException(
                status_code=400,
                detail=f"request body must be a JSON object, got "
                       f"{type(doc).__name__}")
        return doc

    @app.get("/api/projects")
    def get_projects() -> List[Dict]:
        return jsonable(store.list())

    @app.post("/api/projects", status_code=201)
    async def post_project(request: Request) -> Dict:
        doc = await _body(request)
        with _project_errors():
            return jsonable(store.create(doc.get("name"), doc.get("notes", "")))

    @app.get("/api/projects/{project_id}")
    def get_project(project_id: str) -> Dict:
        with _project_errors():
            return jsonable(store.get(project_id))

    @app.patch("/api/projects/{project_id}")
    async def patch_project(project_id: str, request: Request) -> Dict:
        doc = await _body(request)
        with _project_errors():
            return jsonable(store.update(project_id, doc.get("name"),
                                         doc.get("notes")))

    @app.delete("/api/projects/{project_id}")
    def delete_project(project_id: str) -> Dict:
        with _project_errors():
            return jsonable(store.delete(project_id))

    @app.post("/api/projects/{project_id}/specimens", status_code=201)
    async def post_specimen(project_id: str, request: Request) -> Dict:
        doc = await _body(request)
        with _project_errors():
            return jsonable(store.add_specimen(
                project_id, doc.get("section"), doc.get("label"),
                doc.get("notes", "")))

    @app.get("/api/projects/{project_id}/specimens/{specimen_id}")
    def get_specimen(project_id: str, specimen_id: str) -> Dict:
        with _project_errors():
            project = store.get(project_id)
            spec = ProjectStore._find_specimen(project, specimen_id)
        out = dict(spec)
        out["runs"] = [r for r in
                       (_run_view(rid) for rid in spec.get("run_ids") or [])
                       if r is not None]
        out["runs_missing"] = [rid for rid in spec.get("run_ids") or []
                               if _run_view(rid) is None]
        return jsonable(out)

    @app.patch("/api/projects/{project_id}/specimens/{specimen_id}")
    async def patch_specimen(project_id: str, specimen_id: str,
                             request: Request) -> Dict:
        doc = await _body(request)
        with _project_errors():
            return jsonable(store.update_specimen(
                project_id, specimen_id, doc.get("label"), doc.get("notes")))

    @app.delete("/api/projects/{project_id}/specimens/{specimen_id}")
    def delete_specimen(project_id: str, specimen_id: str) -> Dict:
        with _project_errors():
            return jsonable(store.remove_specimen(project_id, specimen_id))

    @app.post("/api/projects/{project_id}/specimens/{specimen_id}/runs",
              status_code=201)
    async def post_specimen_run(project_id: str, specimen_id: str,
                                request: Request) -> Dict:
        doc = await _body(request)
        run_id = doc.get("run_id")
        if not isinstance(run_id, str) or not run_id:
            raise HTTPException(
                status_code=400,
                detail="run_id: required, the id of an existing run")
        with _project_errors():
            return jsonable(store.attach_run(project_id, specimen_id, run_id,
                                             _run_view(run_id)))

    @app.delete(
        "/api/projects/{project_id}/specimens/{specimen_id}/runs/{run_id}")
    def delete_specimen_run(project_id: str, specimen_id: str,
                            run_id: str) -> Dict:
        with _project_errors():
            return jsonable(store.detach_run(project_id, specimen_id, run_id))

    @app.get("/api/projects/{project_id}/compare")
    def get_project_compare(
        project_id: str,
        metric: str = Query("kill_fraction"),
    ) -> Dict:
        """Compare the matched runs of a project — the SAME compound across
        specimens, plus the tissue contrast between them.

        ``metric`` selects which of ``projects._HEADLINE_METRICS`` the
        same-compound cross-specimen rows are FEATURED on. It used to choose
        which metric a cross-compound ranking was taken over; that ranking was
        withdrawn (``results/pk_ranking_robustness.md``,
        ``results/capacity_ranking_robustness.md``) and the response carries no
        ``ranking`` key. Every row in ``kill`` carries all four metrics whatever
        is passed here, so this parameter changes which column is emphasised and
        can never change which rows exist, whether a row is comparable, or what
        the headline concludes. An unknown value is a 400 rather than a silent
        fallback: a caller asking for a metric this build does not have is
        asking about numbers that are not on screen.
        """
        with _project_errors():
            project = store.get(project_id)
            return jsonable(compare_project(project, _run_view, metric=metric))

    @app.get("/api/sections/{section_id}/capability")
    def get_section_capability(section_id: str) -> Dict:
        """What a section's PANEL can support, without creating a project.

        The same block a specimen captures, so the browser can show what a
        section will and will not be able to say BEFORE it is added to a study.
        """
        with _project_errors():
            return jsonable(capability_report(section_id))

    # ------------------------------------------------------------------ #
    # Therapy reports
    # ------------------------------------------------------------------ #
    @app.post("/api/reports", status_code=202)
    async def post_report(request: Request) -> Dict:
        """``{section, drug, spec_overrides}`` -> ``{report_id}``.

        Enqueues exactly TWO ordinary runs of that one compound — the real arm
        and the ``resistance_off`` counterfactual that measures rho — on the SAME
        reconstruction key, so the expensive reconstruction is built once and the
        second arm is a cache hit. Returns immediately; poll
        ``GET /api/reports/{id}``.

        The body took a ``drugs[]`` panel until 2026-08-21. The cross-compound
        ranking was withdrawn (``results/pk_ranking_robustness.md``,
        ``results/capacity_ranking_robustness.md``) and one compound is now the
        only legal shape; a panel is a 400 with that reason, never a truncation.
        """
        body = await _body(request)
        verdict = validate_report_request(body)
        if not verdict["ok"]:
            raise HTTPException(status_code=400, detail={
                "message": "report request rejected",
                **{k: v for k, v in verdict.items() if k != "base_spec"},
            })
        rep = reports.submit(
            section=verdict["section"], drug=verdict["drug"],
            spec_overrides=body.get("spec_overrides") or {},
            base_spec=verdict["base_spec"], label=body.get("label"),
            warnings=verdict["warnings"],
        )
        return {"report_id": rep.report_id, "status": rep.status,
                "run_ids": rep.run_ids, "warnings": verdict["warnings"],
                # The claim the document will make, on the acknowledgement, so a
                # caller that never opens the HTML still meets it.
                "map_claim": MAP_CLAIM}

    @app.get("/api/reports")
    def get_reports(limit: int = Query(200, ge=1, le=5000)) -> List[Dict]:
        return jsonable([reports.view(r, full=False)
                         for r in reports.list(limit=limit)])

    @app.get("/api/reports/{report_id}")
    def get_report(report_id: str) -> Dict:
        return jsonable(reports.view(_require_report(reports, report_id)))

    @app.get("/api/reports/{report_id}/html", response_class=HTMLResponse)
    def get_report_html(report_id: str) -> HTMLResponse:
        rep = _require_report(reports, report_id)
        doc = reports.view(rep)
        if rep.status != "done":
            raise HTTPException(status_code=409, detail={
                "message": "the rendered document is not available yet",
                "report_id": report_id,
                "status": rep.status,
                "progress": doc.get("progress"),
                "reason": (
                    f"a report renders when every arm has completed; this one "
                    f"is {rep.status}."
                    + (f" Error: {rep.error}" if rep.error else "")
                ),
            })
        return HTMLResponse(
            content=render_html(doc),
            headers={"X-Report-Id": report_id,
                     "X-Not-For-Clinical-Use": "true"},
        )

    @app.delete("/api/reports/{report_id}")
    def delete_report(report_id: str) -> Dict:
        rep = _require_report(reports, report_id)
        if rep.status in ("done", "failed", "cancelled"):
            return {"report_id": report_id, "status": rep.status,
                    "cancelled": False,
                    "reason": f"report is already {rep.status}"}
        reports.cancel(report_id)
        return {"report_id": report_id, "status": rep.status, "cancelled": True,
                "reason": ("every queued arm was cancelled; an arm already "
                           "inside simulate_response finishes that solve first")}

    @app.api_route("/api/health", methods=["GET", "HEAD"])
    def health(request: Request) -> Dict:
        """Liveness to anyone; diagnostics only to a caller holding the key.

        THE GATE EXEMPTS THIS ROUTE, SO THE HANDLER HAS TO DO THE WITHHOLDING.
        A platform health checker cannot send ``X-API-Key``, and pointing one
        at a 401 means the platform kills a container that is serving fine --
        measured, twice, on this project. So an anonymous GET gets exactly
        ``{"ok": true}``: enough to prove the process answers, and nothing a
        stranger can learn from.

        The rest -- the runs directory, run counts, the projects root, cache
        statistics, the solver timestep -- is real operational detail about a
        machine, and it still requires the key. ``/healthz`` remains the route
        to prefer for a probe; this one is exempt because a dashboard somewhere
        is already pointed at it and being right in the repo beats being right
        in an argument.
        """
        api_key = os.environ.get("TMESIM_STUDIO_API_KEY", "").strip()
        if api_key and not secrets.compare_digest(
                request.headers.get("x-api-key", ""), api_key):
            return {"ok": True}

        by_status: Dict[str, int] = {}
        for j in jobs.list():
            by_status[j.status] = by_status.get(j.status, 0) + 1
        return jsonable({
            "ok": True,
            "runs_dir": jobs.root,
            "runs": by_status,
            "projects_dir": store.root,
            # thumbnails=False -- a COUNT must not pay for a picture.
            # list() now computes a composition thumbnail per project by
            # reading run documents off disk. That is 0.078 ms on this store
            # and 176 ms cold at 1,250 run documents that cannot supply one.
            # /api/health is polled by a platform health checker; a liveness
            # route whose cost grows with the run store is a liveness route
            # that eventually reports death for a healthy process. This file
            # has already killed a container twice by getting health wrong.
            "projects": len(store.list(thumbnails=False)),
            "reconstruction_cache": RECON_CACHE.stats(),
            "pde_dt_s": PDE_DT_S,
        })

    # /api/protocol/* — the TPL editor's routes. Bound to the SAME job manager,
    # so a protocol sweep shares this process's single worker and the one
    # reconstruction cache with the /api/runs path.
    from .protocol_api import build_protocol_router

    app.include_router(build_protocol_router(jobs))

    # /api/uncertainty/* and /api/runs/{id}/uncertainty -- the interval this
    # project's own findings have been asking for. results/kill_calibration.md
    # says it plainly: "Absolute kill is defensible as an INTERVAL, not as a
    # number, and the interval is wide." That sentence lived in a prose caveat
    # inside a results document, while every metric the engine emitted was a
    # single value.
    #
    # It re-reads a finished run's stored fields and re-propagates two MEASURED
    # spans through that run's own Hill. Nothing is re-simulated. The two kinds
    # of uncertainty are kept structurally separate and a merged or unlabelled
    # interval is refused at the boundary, because conflating them would make
    # the feature worse than no feature: one of them shrinks with a better
    # model and the other only shrinks by measuring this patient.
    from .uncertainty_api import build_uncertainty_router

    app.include_router(build_uncertainty_router(jobs))

    # /api/cross-section — ONE compound, MANY sections. The only cross-run
    # comparison this product draws.
    #
    # Ordering two different COMPOUNDS was withdrawn 2026-08-21 and stays
    # withdrawn. The same compound on two specimens is a different question and
    # a legitimate one, and the interesting output is the DISAGREEMENT: a
    # compound whose reach verdict flips between two patients' tissue is the
    # most informative thing this engine currently says.
    #
    # The module audits its own finished PAYLOAD for a second compound at any
    # depth, rather than trusting its input filter -- the filter is what a
    # future edit would relax, and one grid on one colour scale is a ranking
    # however the compounds got there.
    from .cross_section import build_cross_section_router

    app.include_router(build_cross_section_router(jobs=jobs))

    # /api/limits and /api/account — the deployment's own surface. See
    # app/server/account.py for why these are a separate module, and why the
    # limits are PASSED IN rather than imported there.
    #
    # `_limits` is built here, from the module-level constants the refusal
    # path itself reads, so the number the spec editor shows before a run and
    # the number validate_spec refuses with are the same number. The failure
    # this prevents is the one this project keeps having: a second copy of a
    # value, correct on the day it was written.
    from .account import build_account_router, oom_exposure

    _limits = {
        "max_grid_cells": _MAX_GRID_CELLS,
        "max_grid_cells_default": _DEFAULT_MAX_GRID_CELLS,
        # The RESOLVED source, not one inferred from the value. The inferred
        # version said "default" for an explicit TMESIM_MAX_GRID_CELLS=40000000
        # and could not describe `auto` at all.
        "max_grid_cells_source": _MAX_GRID_CELLS_SOURCE,
        # The point above which a spec is allowed but slow, distinct from the
        # ceiling above which it is refused. A client that shows only the
        # ceiling lets a user walk into a ten-minute reconstruction with no
        # warning.
        "large_grid_cells": _LARGE_GRID_CELLS,
        "min_dx_um": MIN_DX_UM,
        # Published for the same reason the grid ceiling is: a client that
        # cannot see a limit finds it by tripping it.
        "max_n_doses": MAX_N_DOSES,
        "max_total_time_s": MAX_TOTAL_TIME_S,
        "dx_kill_ladder_coarse_um": DX_KILL_LADDER_COARSE_UM,
        "dx_penetration_ladder_coarse_um": DX_PENETRATION_LADDER_COARSE_UM,
        "pde_dt_s": PDE_DT_S,
        # Stated so a client can reimplement it without reading this file, and
        # so the reimplementation has something to be checked AGAINST. The
        # extent comes from the section row's `extent_um`; a section without
        # one cannot be checked client-side at all, and the client must say so
        # rather than assume the spec is fine.
        "grid_cells_formula": (
            "(extent_um[0] / dx_um) * (extent_um[1] / dx_um)"),
        "field_bytes_per_grid_cell": 8,
        "grid_cells_requires": "section extent_um",
        # Whether this server's ceiling still permits the spec that was
        # MEASURED to OOM-kill a container of this host's size. A comparison
        # between a configured number and a recorded measurement -- not a
        # prediction, and not a derived ceiling. See account.oom_exposure.
        "oom_exposure": oom_exposure(
            _MAX_GRID_CELLS, _MAX_GRID_CELLS_SOURCE == "default"),
    }
    app.include_router(build_account_router(
        jobs=jobs, store=store, reports=reports, limits=_limits,
        cache=RECON_CACHE, started_at=_PROCESS_STARTED_AT))

    # The workspace router (docs/WORKSPACE_CONTRACT.md). Bound to the SAME job
    # manager and the SAME project store as the routes above, so the workspace
    # is a view of this process's state and not a second reading of the disk.
    # It mounts THREE routes, and two of them are project routes:
    #
    #   GET  /api/workspace                     the home page's boot call
    #   POST /api/projects/{id}/open            stamps last_opened_at
    #   POST /api/projects/{id}/duplicate       copies specimens + run ids
    #
    # The last two are not beside the other eleven project routes above
    # because they exist only because the home page exists -- nothing in the
    # Studio opens or duplicates a project. This comment is here so the next
    # reader grepping this file for a project route finds them anyway.
    #
    # Gated by the api_key middleware for exactly one reason: they are under
    # /api/. Nothing about this router opts in, and nothing could forget to.
    from .workspace import build_workspace_router

    app.include_router(build_workspace_router(jobs=jobs, store=store))

    # The UI, mounted LAST so /api always wins.
    _mount_web(app)
    return app


def _warm_drugs() -> None:
    try:
        drug_rows()
    except Exception:  # noqa: BLE001 - warming is best-effort
        pass


@contextmanager
def _project_errors():
    """Translate a :class:`ProjectError` into its HTTP status, verbatim.

    The store raises with the status it means (400 bad request, 404 unknown id,
    409 a conflict such as attaching another patient's run to a specimen) and a
    message written for the operator. Re-wording it here would only make the two
    drift apart, so the message travels unchanged.
    """
    try:
        yield
    except ProjectError as exc:
        raise HTTPException(status_code=exc.status, detail={
            "message": exc.message, **exc.detail}) from exc


def _require(jobs: JobManager, run_id: str):
    job = jobs.get(run_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"no such run {run_id!r}")
    return job


def _require_report(reports: ReportManager, report_id: str):
    rep = reports.get(report_id)
    if rep is None:
        raise HTTPException(status_code=404,
                            detail=f"no such report {report_id!r}")
    return rep


def _binary(jobs: JobManager, run_id: str, payload: str, label: str):
    """Serve a persisted binary payload, or explain precisely why it is absent."""
    job = _require(jobs, run_id)
    path = jobs.payload_path(run_id, payload)
    if os.path.exists(path):
        return FileResponse(
            path,
            media_type="application/octet-stream",
            headers={
                # A run's payloads are written once and never change.
                "Cache-Control": "public, max-age=31536000, immutable",
                "X-Run-Id": run_id,
            },
        )
    if job.status != "done":
        raise HTTPException(status_code=409, detail={
            "message": f"{label} is not available yet",
            "run_id": run_id,
            "status": job.status,
            "stage": job.stage,
            "reason": (
                f"payloads are written when a run completes; this run is "
                f"{job.status}." + (f" Error: {job.error}" if job.error else "")
            ),
        })
    body = job.body or {}
    # A payload deleted to reclaim disk is NOT a payload that was never
    # produced, and the branches below would say the latter -- naming a cause
    # ("timelapse_frames was 0, or this run predates the feature") that is
    # simply untrue for a run that captured one and had it pruned. That is the
    # exact failure this repository keeps finding: a confident explanation
    # pointed at the wrong object. So a recorded prune answers first, and says
    # what was removed, when, and how to get it back.
    pruned = (body.get("pruned_payloads") or {}).get(payload)
    if isinstance(pruned, dict):
        raise HTTPException(status_code=410, detail={
            "message": f"{label} was pruned from this run",
            "run_id": run_id,
            "reason": (
                f"this run DID produce {label}; the file was deleted on "
                f"{pruned.get('at', 'an unrecorded date')} to reclaim disk "
                f"({pruned.get('bytes', 0):,} bytes). "
                + str(pruned.get("reason") or "")
            ).strip(),
            "recoverable": "re-run the RunSpec in this record",
        })
    if payload == "field_drug_auc" and not job.spec.get("drug"):
        reason = ("the drug_auc field only exists for a run with a drug; this "
                  "RunSpec had drug=null (reconstruction-only).")
    elif payload == "mechanism":
        # The run itself recorded why it has no attribution — serve that, never
        # a generic apology.
        reason = body.get("mechanism_reason") or (
            "no mechanism attribution was produced by this run "
            "(attribute_mechanism was off, or this run predates the feature).")
    elif payload == "timelapse":
        tl = body.get("timelapse") or {}
        reason = tl.get("reason") or (
            "no drug time-lapse was captured for this run "
            "(timelapse_frames was 0, or this run predates the feature).")
    else:
        reason = (f"{label} was not produced by this run (the underlying field "
                  f"was None in the reconstruction).")
    raise HTTPException(status_code=404, detail={
        "message": f"{label} unavailable", "run_id": run_id, "reason": reason,
    })


#: The bring-your-own-tissue page served at ``GET /ingest``. Self-contained on
#: purpose: it is served by the backend and must work whether or not
#: ``app/web/`` has been built, so it carries its own CSS and script and pulls
#: in nothing external. The only endpoint it touches is
#: ``POST /api/sections/scan``.
_INGEST_PAGE = """<!doctype html><meta charset="utf-8">
<title>Vivum Studio &mdash; bring your own tissue</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
 body{background:#101214;color:#d8dadd;font:13px/1.6 ui-monospace,Menlo,monospace;
      margin:0;padding:32px}
 h1{font-size:15px;font-weight:600;letter-spacing:.04em;text-transform:uppercase;
    margin:0 0 6px;color:#f0f2f4}
 .sub{color:#8b9199;margin:0 0 20px;max-width:76ch}
 .b{border:1px solid #2a2e33;padding:16px;max-width:76ch;margin:0 0 16px}
 label{display:block;margin:10px 0 3px;color:#9aa1a9}
 input[type=text]{width:100%;box-sizing:border-box;background:#0a0c0e;
   border:1px solid #2a2e33;color:#d8dadd;padding:7px 8px;font:inherit}
 .row{display:flex;gap:10px;align-items:center;margin-top:14px;flex-wrap:wrap}
 button{background:#1b2733;border:1px solid #33506e;color:#cfe2f6;padding:7px 14px;
   font:inherit;cursor:pointer}
 button:hover:not([disabled]){background:#22344a}
 button[disabled]{opacity:.5;cursor:default}
 .chk{color:#9aa1a9;font-size:12px}
 pre{background:#0a0c0e;border:1px solid #2a2e33;padding:12px;overflow-x:auto;
   white-space:pre;margin:0}
 .err{border-color:#6e3333;color:#f0b8b8}
 .ok{color:#9ad7a0}
 code{color:#7fb6e8} a{color:#7fb6e8}
 .note{color:#8b9199;font-size:12px;max-width:76ch}
</style>
<h1>Bring your own tissue</h1>
<p class="sub">Point this at a directory holding a 10x Xenium or Visium export
(or a bare MEX triple plus a coordinates table) and it will tell you what the
sample <em>is</em> and what it <strong>cannot</strong> support. Registering adds
it to the section list, so it becomes selectable in Studio after a reload.
Same endpoint as <code>POST /api/sections/scan</code> and the same report
<code>vivum ingest</code> prints &mdash; this page just saves you the
<code>curl</code>.</p>

<div class="b">
 <label for="path">Sample directory (absolute path on the machine running this
 server)</label>
 <input type="text" id="path" spellcheck="false"
        placeholder="/data/xenium_run_2026_08/output-XETG00001">

 <label for="sid">Section id (optional &mdash; defaults to the directory name)</label>
 <input type="text" id="sid" spellcheck="false">

 <label for="label">Label shown in the picker (optional)</label>
 <input type="text" id="label" spellcheck="false">

 <label for="key">X-API-Key (only if this server was started with one)</label>
 <input type="text" id="key" spellcheck="false" autocomplete="off">

 <div class="row">
  <button id="scan">Scan &mdash; read only</button>
  <button id="reg">Scan + register</button>
  <label class="chk"><input type="checkbox" id="deep"> measure detection rates
   (full ingest &mdash; slower, but turns &ldquo;the gene is on the panel&rdquo;
   into a measured firing rate)</label>
 </div>
</div>

<div id="out"></div>
<p class="note">Registering never shadows a built-in section, and a directory
that is not ingestable is answered with a reason rather than a broken section.
Remove an ingested section with
<code>DELETE /api/sections/&lt;id&gt;</code>.</p>

<script>
(function () {
  var out = document.getElementById("out");
  var els = ["path", "sid", "label", "key", "deep", "scan", "reg"]
    .reduce(function (a, k) { a[k] = document.getElementById(k); return a; }, {});

  function esc(s) {
    return String(s).replace(/[&<>]/g, function (c) {
      return { "&": "&amp;", "<": "&lt;", ">": "&gt;" }[c];
    });
  }
  function box(cls, html) { return '<div class="b ' + cls + '">' + html + "</div>"; }
  function show(html) { out.innerHTML = html; }

  function render(res, registered) {
    var html = "";
    if (res.error) {
      html += box("err", "<strong>Not ingested.</strong> " + esc(res.error));
    } else if (registered && res.ok) {
      html += box("", '<span class="ok">Registered as <code>' +
        esc(res.section_id) + "</code>.</span> Reload Studio and it will be in " +
        "the section picker." +
        (res.section && res.section.runnable === false
          ? " <strong>Inspectable only:</strong> the reconstruction pipeline " +
            "has no route for this platform."
          : ""));
    }
    var det = res.detection || (res.section ? {
      platform: res.section.platform, layout: res.section.layout,
      confidence: res.section.detect_confidence,
      issues: res.section.detect_issues
    } : null);
    if (det) {
      var bits = [];
      if (det.platform) bits.push("platform " + esc(det.platform));
      if (det.layout) bits.push("layout " + esc(det.layout));
      if (det.confidence !== undefined && det.confidence !== null) {
        bits.push("confidence " + esc(det.confidence));
      }
      if (bits.length) html += box("", bits.join(" &middot; "));
      var iss = det.issues || det.detect_issues;
      if (iss && iss.length) {
        html += box("err", "<strong>Issues</strong><pre>" +
          esc(iss.join("\\n")) + "</pre>");
      }
    }
    if (res.capability_text) {
      html += "<pre>" + esc(res.capability_text) + "</pre>";
    }
    if (!html) {
      html = box("", "The server returned no capability report:<pre>" +
        esc(JSON.stringify(res, null, 1)) + "</pre>");
    }
    show(html);
  }

  function go(registered) {
    var path = els.path.value.trim();
    if (!path) { show(box("err", "Enter a sample directory path.")); return; }
    var body = {
      path: path,
      measure_detection_rates: !!els.deep.checked,
      register: !!registered
    };
    if (els.sid.value.trim()) body.section_id = els.sid.value.trim();
    if (els.label.value.trim()) body.label = els.label.value.trim();
    var headers = { "Content-Type": "application/json" };
    if (els.key.value.trim()) headers["X-API-Key"] = els.key.value.trim();

    els.scan.disabled = els.reg.disabled = true;
    show(box("", (registered ? "Ingesting" : "Scanning") + " &hellip; " +
      (body.measure_detection_rates || registered
        ? "a full ingest of a large section can take a minute."
        : "this reads filenames and the panel only.")));

    fetch("/api/sections/scan", {
      method: "POST", headers: headers, body: JSON.stringify(body)
    }).then(function (r) {
      return r.json().then(function (j) { return { status: r.status, body: j }; });
    }).then(function (r) {
      if (r.status === 401) {
        show(box("err", "401 &mdash; this server is running with an API key. " +
          "Paste it into the X-API-Key field above."));
        return;
      }
      if (r.status === 403 && r.body && r.body.detail && r.body.detail.message) {
        // Path confinement (TMESIM_STUDIO_DATA_ROOTS). The server already
        // wrote a sentence telling the operator what to do; show THAT, not a
        // JSON blob of it.
        show(box("err", "403 &mdash; " + esc(r.body.detail.message)));
        return;
      }
      if (r.status >= 400) {
        show(box("err", esc(r.status + " " +
          (r.body && r.body.detail ? JSON.stringify(r.body.detail) : ""))));
        return;
      }
      render(r.body, registered);
    }).catch(function (e) {
      show(box("err", "Request failed: " + esc(e)));
    }).then(function () {
      els.scan.disabled = els.reg.disabled = false;
    });
  }

  els.scan.addEventListener("click", function () { go(false); });
  els.reg.addEventListener("click", function () { go(true); });
})();
</script>
"""


_PLACEHOLDER = """<!doctype html><meta charset="utf-8">
<title>Vivum Studio</title>
<style>
 body{background:#101214;color:#d8dadd;font:13px/1.5 ui-monospace,Menlo,monospace;
      margin:0;padding:32px}
 h1{font-size:15px;font-weight:600;letter-spacing:.04em;text-transform:uppercase;
    margin:0 0 16px;color:#f0f2f4}
 code{color:#7fb6e8} a{color:#7fb6e8} li{margin:2px 0}
 .b{border:1px solid #2a2e33;padding:14px 16px;max-width:760px}
</style>
<h1>Vivum Studio &mdash; backend up, front end not built</h1>
<div class="b">
<p>The API is serving. <code>app/web/index.html</code> does not exist yet, so this
placeholder is shown instead.</p>
<p>Live endpoints:</p>
<ul>
 <li><a href="/api/sections">/api/sections</a></li>
 <li><a href="/api/drugs">/api/drugs</a></li>
 <li><a href="/api/defaults">/api/defaults</a></li>
 <li><a href="/api/runs">/api/runs</a></li>
 <li><a href="/api/health">/api/health</a></li>
 <li><a href="/ingest">/ingest</a> &mdash; load YOUR OWN 10x export</li>
 <li><a href="/docs">/docs</a> &mdash; interactive OpenAPI</li>
</ul>
</div>
"""


def _mount_web(app: FastAPI) -> None:
    """Serve ``app/web/`` — the workspace home at ``/``, the Studio at
    ``/studio``, and everything else as a static asset.

    THE ROUTE TABLE IS ``docs/WORKSPACE_CONTRACT.md`` §1.

    Until 2026-09-20 this mounted ``StaticFiles(html=True)`` at ``/`` and that
    was the whole routing story: ``/`` served ``index.html``, the Studio, which
    is the product's HARDEST screen. A first-time visitor landed on a
    three-rail workstation that assumes you already know which section you
    want, which compound, and what a RunSpec is.

    Now ``/`` is the workspace and the Studio is one click in. Three
    properties, each load-bearing:

    * **``/index.html`` keeps working, byte for byte.** The static mount still
      serves it, so every existing bookmark and every link anyone has shared
      survives. That is the contract's compatibility bar and it is checked by
      ``tests/test_home_page.py``.
    * **``/studio`` is an explicit route**, not a redirect, so the Studio has a
      URL that reads like what it is and can carry ``?project=``.
    * **Both are ungated**, like every static asset. A page must be able to
      load in order to prompt for a key; a key gate in front of the login
      screen is a loop.

    The explicit routes are declared BEFORE the catch-all mount because
    starlette matches in registration order — a ``StaticFiles`` mount at ``/``
    swallows every path under it, so anything registered after it is
    unreachable. That ordering is the entire mechanism here, and reversing it
    fails silently: the routes still exist, they are simply never matched.
    """
    web_dir = os.path.join(repo_root(), "app", "web")
    index = os.path.join(web_dir, "index.html")
    home = os.path.join(web_dir, "home.html")
    if not (os.path.isdir(web_dir) and os.path.exists(index)):
        # No front end built. Keep the original behaviour exactly: a single
        # placeholder page at `/` explaining what is missing, rather than a
        # bare 404 that looks like the server is broken.
        @app.get("/", response_class=HTMLResponse, include_in_schema=False)
        def _placeholder() -> str:
            return _PLACEHOLDER

        return

    from fastapi.staticfiles import StaticFiles

    def _page(path: str):
        return FileResponse(path, media_type="text/html",
                            headers={"Cache-Control": "no-cache"})

    # The Studio, at a name of its own. Registered whether or not home.html
    # exists, so this route does not silently depend on the other page.
    @app.get("/studio", response_class=HTMLResponse, include_in_schema=False)
    def studio_page():
        return _page(index)

    if os.path.exists(home):
        # The workspace. `/` only becomes the home page when the home page is
        # actually on disk -- a half-deployed image serves the Studio at `/`
        # exactly as it always did, rather than 500ing on the front door.
        @app.get("/", response_class=HTMLResponse, include_in_schema=False)
        def home_page():
            return _page(home)

    app.mount("/", StaticFiles(directory=web_dir, html=True), name="web")


_app: Optional[FastAPI] = None


def __getattr__(name: str):
    """Lazily build the module-level ``app`` for ``uvicorn app.server.api:app``.

    PEP 562 module ``__getattr__`` so that merely importing this module does not
    spin up a job worker or scan ``app/runs/`` — tests build their own
    :class:`~app.server.jobs.JobManager` through :func:`create_app`.
    """
    if name == "app":
        global _app
        if _app is None:
            _app = create_app()
        return _app
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
