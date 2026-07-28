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
POST         /api/validate                       {ok, errors[], warnings[]}
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
from contextlib import asynccontextmanager
from typing import Dict, List, Optional

from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse

from .cache import (
    OXYGEN_MODES,
    RECON_CACHE,
    SECTIONS,
    SENSITIVITY_MODES,
    VASCULAR_MODES,
    list_sections,
    make_key,
    repo_root,
    section_genes,
)
from .jobs import PDE_DT_S, JobManager, jsonable
from .payloads import FIELD_NAMES

__all__ = ["create_app", "app", "validate_spec", "drug_rows", "default_spec"]

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
    from tmesim.reconstruct.real_pipeline import OXYGEN_MODE_DEFAULT

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
    "oxygen_mode", "mgmt", "pk", "total_time_s",
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
                      f"Known ids: {sorted(SECTIONS)}")
    elif section_id not in SECTIONS:
        errors.append(f"section: unknown id {section_id!r}. Known ids: "
                      f"{sorted(SECTIONS)}")
    else:
        section_row = next(
            (s for s in list_sections() if s["id"] == section_id), None)
        if section_row is not None and not section_row["available"]:
            errors.append(
                f"section {section_id!r} is not loadable: "
                f"{section_row.get('reason')}"
            )

    source = SECTIONS.get(section_id, {}).get("source") if section_id else None

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
            if source == "visium" and dx < 50.0:
                warnings.append(
                    f"dx_um = {dx:g} um over-resolves a 100 um Visium spot "
                    "pitch; reconstruct_real falls back to dx = 100 um for "
                    "Visium below 50 um, so this value will NOT be used."
                )
            if section_row and section_row.get("extent_um") and dx > 0:
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
    _validate_cfl(section_id, dx, vasc_mode, oxy_mode, errors, warnings)

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

    return {"ok": not errors, "errors": errors, "warnings": warnings}


def _validate_cfl(section_id, dx, vasc_mode, oxy_mode,
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
       and the domain and never reads oxygen, so those keys share this spec's
       vascular source and differ only in mesh. A violation is still an ERROR,
       labelled with the key it was measured on.
    3. **Only a different ``vascular_mode``**, or nothing at all -> a WARNING
       saying the CFL is unverified and why. No speed is ever assumed.
    """
    if dx is None or not section_id or section_id not in SECTIONS:
        return

    key = make_key(section_id, dx, vasc_mode, oxy_mode)
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
# App
# --------------------------------------------------------------------------- #
def create_app(manager: Optional[JobManager] = None) -> FastAPI:
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
    jobs = manager if manager is not None else JobManager()
    app.state.jobs = jobs
    app.state.cache = RECON_CACHE

    # ----------------------------------------------------------------- #
    @app.get("/api/sections")
    def get_sections() -> List[Dict]:
        return jsonable(list_sections())

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

    @app.get("/api/health")
    def health() -> Dict:
        by_status: Dict[str, int] = {}
        for j in jobs.list():
            by_status[j.status] = by_status.get(j.status, 0) + 1
        return jsonable({
            "ok": True,
            "runs_dir": jobs.root,
            "runs": by_status,
            "reconstruction_cache": RECON_CACHE.stats(),
            "pde_dt_s": PDE_DT_S,
        })

    # The UI, mounted LAST so /api always wins.
    _mount_web(app)
    return app


def _warm_drugs() -> None:
    try:
        drug_rows()
    except Exception:  # noqa: BLE001 - warming is best-effort
        pass


def _require(jobs: JobManager, run_id: str):
    job = jobs.get(run_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"no such run {run_id!r}")
    return job


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
    reason = (
        "the drug_auc field only exists for a run with a drug; this RunSpec "
        "had drug=null (reconstruction-only)."
        if payload == "field_drug_auc" and not job.spec.get("drug")
        else f"{label} was not produced by this run (the underlying field was "
             f"None in the reconstruction)."
    )
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
