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
POST         /api/sections/scan                  ingest ANY directory: detect() +
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

The static front end in ``app/web/`` is mounted last, at ``/``, so it can never
shadow an ``/api`` route.
"""
from __future__ import annotations

import math
import os
import threading
from contextlib import asynccontextmanager, contextmanager
from typing import Dict, List, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse

from .cache import (
    OXYGEN_MODES,
    RECON_CACHE,
    SECTIONS,
    all_sections,
    lookup_section,
    SENSITIVITY_MODES,
    STATE_MODES,
    VASCULAR_MODES,
    list_sections,
    make_key,
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
    DEFAULT_PANEL,
    REPORT_SENSITIVITY_MODE,
    ReportManager,
    render_html,
)

__all__ = ["create_app", "app", "validate_spec", "drug_rows", "default_spec",
           "validate_report_request"]

#: Grid resolutions the repo has convergence evidence for
#: (results/numerical_rigor.md section 2).
DX_KILL_CONVERGED_UM = 20.0
DX_PENETRATION_CONVERGED_UM = 15.0

#: Above this the grid is big enough that a run is measured in minutes; the UI
#: gets a warning rather than a refusal, because it is a cost, not an error.
_LARGE_GRID_CELLS = 1_500_000


# --------------------------------------------------------------------------- #
# Drugs
# --------------------------------------------------------------------------- #
_drug_rows_cache: Optional[List[Dict]] = None


def drug_rows(force: bool = False) -> List[Dict]:
    """``GET /api/drugs`` payload: the real library, compiled once and memoised.

    ``ic50_uM`` and ``ic50_source`` come from
    ``molecular.potency.resolve_ic50`` — measured GDSC glioma median where the
    drug is screened, curated literature where it is not, and a flagged default
    otherwise. ``bbb_penetrance`` comes from the fully compiled
    ``calibrated_drugparams -> apply_barrier_to_drugparams`` chain, i.e. the
    exact ``DrugParams`` a run would simulate. A drug that fails to compile is
    returned with nulls and a ``reason``, never dropped and never faked.
    """
    global _drug_rows_cache
    if _drug_rows_cache is not None and not force:
        return _drug_rows_cache

    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.library import GBM_DRUGS, get_molecule
    from tmesim.molecular.potency import resolve_ic50
    from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

    try:
        alpha = float(calibrate_alpha())
        alpha_error = None
    except Exception as exc:  # noqa: BLE001
        alpha, alpha_error = None, f"{type(exc).__name__}: {exc}"

    rows: List[Dict] = []
    for name in sorted(GBM_DRUGS):
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
            mol = get_molecule(name)
            row["mw"] = _num(mol.get("mw"))
            row["logp"] = _num(mol.get("logp"))
            row["is_antibody"] = bool(mol.get("is_antibody"))
            row["smiles_source"] = mol.source
        except Exception as exc:  # noqa: BLE001
            row["reason"] = f"descriptor build failed: {type(exc).__name__}: {exc}"

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
            row["reason"] = f"IC50 unresolved: {type(exc).__name__}: {exc}"

        if alpha is None:
            row["reason"] = (
                f"bbb_penetrance unavailable: calibrate_alpha() failed "
                f"({alpha_error})"
            )
        else:
            try:
                dp = apply_barrier_to_drugparams(
                    calibrated_drugparams(name, alpha), get_molecule(name))
                row["bbb_penetrance"] = _num(dp.bbb_penetrance)
                row["diffusivity_um2_s"] = _num(dp.diffusivity)
                # The compiled ic50 is the resolved one; report it from the
                # object the engine will actually consume.
                row["ic50_uM"] = _num(dp.ic50)
            except Exception as exc:  # noqa: BLE001
                row["reason"] = (
                    f"compilation failed: {type(exc).__name__}: {exc}")
        rows.append(row)

    _drug_rows_cache = rows
    return rows


def _num(v):
    """float(v) unless it is missing or non-finite, in which case ``None``."""
    if v is None:
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def _known_drugs() -> List[str]:
    from tmesim.molecular.library import GBM_DRUGS

    return sorted(GBM_DRUGS)


# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
def default_spec() -> Dict:
    """A valid ``RunSpec`` carrying the current production defaults.

    Sourced from the production entry points, not retyped:
    ``dx = 12 um`` and a 6 h window from ``experiments/patient_report.py``
    (the converged resolution per ``results/numerical_rigor.md`` section 2), the
    ``'uniform'`` sensitivity mode it ranks in, ``reconstruct.real_fields``'s
    default ``'discrete'`` vascular mode, ``real_pipeline.OXYGEN_MODE_DEFAULT``,
    and the ``core.PKConfig`` field defaults.
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
        "dx_um": 12.0,
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
        "total_time_s": 3600.0 * 6,
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


def _as_number(value, name: str, errors: List[str]) -> Optional[float]:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append(f"{name}: expected a number, got {type(value).__name__} "
                      f"({value!r})")
        return None
    v = float(value)
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
            "sensitivity_mode='determinant' is NOT comparable across drugs: a "
            "drug whose determinant is on the panel is scaled down while a drug "
            "whose determinant is absent keeps 1.0. Use 'uniform' for any "
            "cross-drug ranking."
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
            elif dx > DX_KILL_CONVERGED_UM:
                warnings.append(
                    f"dx_um = {dx:g} um is coarser than the {DX_KILL_CONVERGED_UM:g} "
                    "um at which kill_fraction is grid-converged "
                    "(results/numerical_rigor.md section 2). Both kill_fraction "
                    "and penetration_index should be read as upper bounds."
                )
            elif dx > DX_PENETRATION_CONVERGED_UM:
                warnings.append(
                    f"dx_um = {dx:g} um: kill_fraction is converged, but "
                    "penetration_index is not below ~15 um and is typically "
                    "15-25% high."
                )
            # `dx` is cleared above when it is non-positive, so both of the
            # remaining checks have to re-test it: a dx of -1 is an ERROR, and
            # crashing on it here would turn a 400 into a 500.
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
                if n_cells_grid > _LARGE_GRID_CELLS:
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
            elif float(nd) != int(nd):
                errors.append(f"pk.n_doses: must be a whole number, got {nd}")
            elif int(nd) < 1:
                errors.append(f"pk.n_doses: must be >= 1, got {int(nd)}")
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
    elif float(frames) != int(frames):
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
            "routes with everything else identical. The CROSS-DRUG RANKING is "
            "INVARIANT across all three and is the defensible output; the level "
            "is not (results/cell_state_impact.md §3.2, §3.3, §5)."
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
_REPORT_KEYS = {"section", "drugs", "spec_overrides", "label"}

#: RunSpec fields a report OWNS and will not let a caller override, with the
#: reason each one is fixed. ``drug`` and ``resistance_off`` are swept per arm;
#: ``sensitivity_mode`` is what makes the ranking a ranking.
_REPORT_FIXED = {
    "drug": "a report sweeps the whole panel; drugs[] selects it",
    "resistance_off": ("both arms are generated per drug — the real run and the "
                       "resistance-off counterfactual the attribution needs"),
    "sensitivity_mode": (
        f"forced to {REPORT_SENSITIVITY_MODE!r}: it is the only mode in which a "
        "cross-drug ranking is fair, because an identical per-cell factor "
        "cancels out of the comparison. 'determinant' scales some drugs by "
        "their measured determinant and leaves others at 1.0, which is a "
        "measurement-coverage artifact, not a ranking "
        "(results/real_drug_panel.md)."
    ),
}


def validate_report_request(body: Dict) -> Dict:
    """``POST /api/reports`` body -> ``{ok, errors[], warnings[], section,
    drugs[], base_spec}``.

    ``base_spec`` is a full, already-``validate_spec``-clean RunSpec with
    ``drug=null``; the manager clones it per arm. Every arm therefore shares one
    reconstruction cache key by construction.
    """
    errors: List[str] = []
    warnings: List[str] = []
    if not isinstance(body, dict):
        return {"ok": False, "errors": [
            f"report request must be a JSON object, got {type(body).__name__}"],
            "warnings": [], "section": None, "drugs": [], "base_spec": None}

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

    raw_drugs = body.get("drugs")
    if raw_drugs is None:
        drugs = list(DEFAULT_PANEL)
        warnings.append(
            f"drugs[] not supplied: defaulting to the {len(DEFAULT_PANEL)}-drug "
            f"panel {list(DEFAULT_PANEL)} that experiments/patient_report.py "
            "ranks.")
    elif not isinstance(raw_drugs, list) or not raw_drugs:
        errors.append("drugs: expected a non-empty array of drug names")
        drugs = []
    else:
        known = _known_drugs()
        drugs, seen = [], set()
        for d in raw_drugs:
            if not isinstance(d, str):
                errors.append(f"drugs: expected drug names, got "
                              f"{type(d).__name__} ({d!r})")
                continue
            name = d.strip().lower()
            if name not in known:
                errors.append(f"drugs: unknown drug {d!r}. Available: {known}")
                continue
            if name in seen:
                warnings.append(f"drugs: {name!r} listed more than once; the "
                                "duplicate was dropped.")
                continue
            seen.add(name)
            drugs.append(name)
        if len(drugs) == 1:
            warnings.append(
                "a one-drug panel has no ranking and no panel mean, so every "
                "attribution term is 0 by construction. The attribution is a "
                "CROSS-drug decomposition; use two or more drugs.")

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

    # Validate the spec the arms will actually run, with a real drug in it, so
    # a report can never enqueue 2N runs that /api/validate would have refused.
    if drugs and not errors:
        probe = dict(base, drug=drugs[0])
        verdict = validate_spec(probe)
        errors.extend(verdict["errors"])
        for w in verdict["warnings"]:
            if w not in warnings:
                warnings.append(w)

    return {"ok": not errors, "errors": errors, "warnings": warnings,
            "section": section, "drugs": drugs, "base_spec": base}


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
        title="tmesim Studio API",
        version="1",
        description=__doc__,
        lifespan=lifespan,
    )
    # Bound to 127.0.0.1 by __main__; permissive CORS so the front end can be
    # served from a different dev port without a proxy.
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_methods=["*"],
        allow_headers=["*"],
    )

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
            if request.url.path.startswith("/api/") and request.method != "OPTIONS":
                if request.headers.get("x-api-key", "") != api_key:
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
    @app.get("/api/sections")
    def get_sections() -> List[Dict]:
        return jsonable(list_sections())

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
        """
        try:
            body = await request.json()
        except Exception as exc:  # noqa: BLE001
            raise HTTPException(
                status_code=400,
                detail=f"request body is not valid JSON: {exc}") from exc
        if not isinstance(body, dict) or not body.get("path"):
            raise HTTPException(
                status_code=400,
                detail="body must be an object with a 'path' string naming the "
                       "sample directory to scan")

        from . import ingest_registry

        path = str(body["path"])
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

    @app.get("/api/drugs")
    def get_drugs() -> List[Dict]:
        return jsonable(drug_rows())

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
        job = _require(jobs, run_id)
        return {"run_id": run_id, "status": job.status, "error": job.error,
                "traceback": job.traceback}

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
        """``{section, drugs[], spec_overrides}`` -> ``{report_id}``.

        Enqueues 2 x len(drugs) ordinary runs (the real arm and the
        resistance-off counterfactual for each drug) on the SAME reconstruction
        key, so the expensive reconstruction is built once and every later arm
        is a cache hit. Returns immediately; poll ``GET /api/reports/{id}``.
        """
        body = await _body(request)
        verdict = validate_report_request(body)
        if not verdict["ok"]:
            raise HTTPException(status_code=400, detail={
                "message": "report request rejected",
                **{k: v for k, v in verdict.items() if k != "base_spec"},
            })
        rep = reports.submit(
            section=verdict["section"], drugs=verdict["drugs"],
            spec_overrides=body.get("spec_overrides") or {},
            base_spec=verdict["base_spec"], label=body.get("label"),
            warnings=verdict["warnings"],
        )
        return {"report_id": rep.report_id, "status": rep.status,
                "run_ids": rep.run_ids, "warnings": verdict["warnings"]}

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

    @app.get("/api/health")
    def health() -> Dict:
        by_status: Dict[str, int] = {}
        for j in jobs.list():
            by_status[j.status] = by_status.get(j.status, 0) + 1
        return jsonable({
            "ok": True,
            "runs_dir": jobs.root,
            "runs": by_status,
            "projects_dir": store.root,
            "projects": len(store.list()),
            "reconstruction_cache": RECON_CACHE.stats(),
            "pde_dt_s": PDE_DT_S,
        })

    # /api/protocol/* — the TPL editor's routes. Bound to the SAME job manager,
    # so a protocol sweep shares this process's single worker and the one
    # reconstruction cache with the /api/runs path.
    from .protocol_api import build_protocol_router

    app.include_router(build_protocol_router(jobs))

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


_PLACEHOLDER = """<!doctype html><meta charset="utf-8">
<title>tmesim Studio</title>
<style>
 body{background:#101214;color:#d8dadd;font:13px/1.5 ui-monospace,Menlo,monospace;
      margin:0;padding:32px}
 h1{font-size:15px;font-weight:600;letter-spacing:.04em;text-transform:uppercase;
    margin:0 0 16px;color:#f0f2f4}
 code{color:#7fb6e8} a{color:#7fb6e8} li{margin:2px 0}
 .b{border:1px solid #2a2e33;padding:14px 16px;max-width:760px}
</style>
<h1>tmesim Studio &mdash; backend up, front end not built</h1>
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
 <li><a href="/docs">/docs</a> &mdash; interactive OpenAPI</li>
</ul>
</div>
"""


def _mount_web(app: FastAPI) -> None:
    """Serve ``app/web/`` at ``/`` if the front end exists."""
    web_dir = os.path.join(repo_root(), "app", "web")
    index = os.path.join(web_dir, "index.html")
    if os.path.isdir(web_dir) and os.path.exists(index):
        from fastapi.staticfiles import StaticFiles

        app.mount("/", StaticFiles(directory=web_dir, html=True), name="web")
        return

    @app.get("/", response_class=HTMLResponse, include_in_schema=False)
    def _placeholder() -> str:
        return _PLACEHOLDER


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
