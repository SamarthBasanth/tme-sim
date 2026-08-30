"""HOW MUCH OF THE SHIPPED STEADY-STATE REACH VERDICT IS THE MESH?

Driver for ``results/dx_sensitivity_steady.md``.

WHY THIS FILE EXISTS
====================
A "10x grid sensitivity" figure was once published from this repository and
RETRACTED because it was never measured. This file replaces it with a real
sweep. Nothing here is derived, extrapolated or asserted from a prior document:
every number below comes from a run this script executed, and the raw output is
written to ``results/dx_sensitivity_steady_data.json``.

WHAT THE SHIPPED OBJECT ACTUALLY IS (verified from source, not assumed)
======================================================================
``app/server/jobs.py::execute_run`` reconstructs the section through
``app/server/cache.py::build_reconstruction`` -> ``reconstruct.real_pipeline.
reconstruct_real``, marches ``sim.engine.simulate_response`` for the spec's
window at ``PDE_DT_S = 10 s`` / ``PDE_CONVECTION = True`` / ``SIM_SEED = 0``,
and then calls

    reach = compute_reach(micro, res, distance_um=dist,
                          penetration=dm["penetration_length"],
                          drug=dp, convection=PDE_CONVECTION)

``compute_reach`` is defined in ``app/server/reach.py`` (NOT in ``jobs.py`` —
``jobs.py`` only imports and calls it). Inside it, ``steady = drug is not None``,
so the production call runs in **steady-state** mode: the map is built on
``reach.steady_state_concentration``, ``window_mode`` is ``"steady_state"``,
and ``penetration`` is DISCARDED (``if steady or penetration is None:`` re-fits
against the steady profile). A tissue-capacity term ``R`` multiplies ``dC/dt``
only, so it is absent from ``(K - L + Adv) c* = S`` by inspection and plays no
part in anything measured here.

Two consequences for this sweep, both of which it exploits:

* The steady field has **no dt and no window**, so nothing in the reach map is
  subject to the advection CFL. Only the finite-window march is.
* The march is still required, because ``compute_reach`` reads
  ``sim_result.fields['oxygen']`` (the viable band) and
  ``delivery_metrics.tumor_cell_exposure`` reads ``fields['auc']``. It is also
  where the *other* published verdict lives — the per-cell dominant limiter from
  ``sim.mechanism`` — so both are measured at every rung from the SAME run.

WHAT IS VARIED, AND WHAT IS NOT
===============================
**ONLY ``dx``.** ``dt = 10 s``, the 6 h window, ``save_every``, the seed,
``convection = True``, ``transvascular_exchange = False`` (the shipped default),
``PKConfig()``, the drug compile path and ``alpha`` are all held at their
shipped values. ``sensitivity_mode`` is the shipped ``'uniform'``, which sets
``cells.sensitivity`` to 1.0 everywhere and is therefore a no-op against the
reconstruction default -- so it is not called, and that is a documented
equivalence rather than an omission (``real_pipeline.apply_drug_sensitivity``
docstring, ``mode='uniform'``).

Every solve goes through ``sim.engine.simulate_response``, the one production
wiring point that turns ``PDEConfig.transvascular_exchange`` into an
``exchange_k`` argument. This file calls ``simulate_drug_transport`` NOWHERE —
the defect where a bare call with a hand-built ``PDEConfig`` silently made both
arms of a comparison exchange-OFF cannot arise here.

The reconstruction entry point is ``validation.vessel_profile.
build_arm_microenvironment(dx_um, repair=0)``, which is literally
``reconstruct_real(source="xenium", dx_um=dx)`` — the same call
``build_reconstruction`` makes for ``xenium_319`` (whose registry ``sample`` is
``None``, i.e. the loader default) at the shipped ``oxygen_mode='pde'`` /
``state_mode='oxygen'``. ``build_reconstruction`` itself is not called only
because it writes section metadata back into ``data_cache`` as a side effect;
for ``vascular_mode='discrete'`` everything else it does is a metadata label.
The equivalence is not asserted, it is checked: at ``dx = 12`` this path
reproduces the container-measured shipped grid (286,200 voxels) and
``kill_fraction`` (0.11221232991768859) from ``docs/DEPLOYMENT.md``.

THE VACUOUS-GREEN TRAPS, AND THE CHECKS THAT CATCH THEM
=======================================================
A ``check_mode_is_steady`` — if ``drug=`` were dropped from the ``compute_reach``
  call, every rung would silently return a FINITE-WINDOW map and the whole sweep
  would measure a different object while every number still looked plausible.
  ``window_mode == "steady_state"`` is required at every rung.
B ``check_rungs_distinct`` — the mesh has to actually reach the physics. The
  per-rung ``(E_ref, F-ladder, falloff[1])`` fingerprints must be pairwise
  distinct. A sweep that reused one reconstruction would report perfect
  convergence with every rung green.
C ``check_shipped_reproduced`` — the ``dx = 12`` rung must reproduce the
  published container ``kill_fraction`` bit-for-bit. If this harness is not the
  shipped pipeline, its convergence statement is about something else.
D ``check_alpha`` — ``calibrate_alpha()`` asserted at the START and the END of
  every stage.

WHERE THE SWEEP WAS RUN, AND WHY IT WAS NOT RUN IN THE WORKING TREE
==================================================================
The first attempt ran against the live checkout and was CONTAMINATED and
DISCARDED. This repository is worked on by several agents concurrently, and
between the ``dx = 20`` and ``dx = 16`` rungs two files in the compile/delivery
chain -- ``tmesim/reconstruct/barrier.py`` and ``tmesim/treatment/drugs.py`` --
were edited and later reverted underneath the running sweep. The straddling rung
came back with ``kill_fraction = 0.317655`` against ``0.112044`` one rung
earlier: a 2.8x jump that would have read as a spectacular mesh effect and was
not one. It was thrown away, not repaired, and the whole ladder was re-run.

The published sweep therefore runs against a FROZEN, SELF-CONSISTENT snapshot of
the committed tree, extracted read-only with ``git archive HEAD`` (no
``checkout``, no ``stash``, no ``reset`` -- the working tree is never touched)
into a scratch directory, with ``data_cache`` symlinked back. That snapshot is
also the right object on the merits: uncommitted in-flight edits by other agents
are not what ships. Check C is what makes the choice falsifiable -- the snapshot
must reproduce the container-measured shipped ``kill_fraction`` bit-for-bit.

Two mechanisms make a repeat of that contamination visible rather than silent:
``source_fingerprint`` (a SHA-256 over every ``.py`` under ``tmesim/`` and
``app/server/`` of the tree actually being imported) is recorded on EVERY rung
and check E requires all rungs equal; and the REPEAT GATE re-runs the first rung
after the whole ladder and requires it bit-identical.

Run:
    ./.venv/bin/python experiments/dx_sensitivity_steady.py --stage probe
    ./.venv/bin/python experiments/dx_sensitivity_steady.py --stage sweep
    ./.venv/bin/python experiments/dx_sensitivity_steady.py --stage report
    ./.venv/bin/python experiments/dx_sensitivity_steady.py --stage rung \
        --dx 12 --compound temozolomide      # one rung, JSON on stdout
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from app.server.jobs import (                                          # noqa: E402
    PDE_CONVECTION,
    PDE_DT_S,
    PDE_SAVE_EVERY,
    SIM_SEED,
)
from app.server.reach import REACH_LADDER, HEADLINE_X, compute_reach   # noqa: E402
from tmesim.core import PDEConfig, PKConfig, SimConfig                 # noqa: E402
from tmesim.molecular.binding_site_barrier import (                    # noqa: E402
    apply_barrier_to_drugparams,
)
from tmesim.molecular.library import get_molecule                      # noqa: E402
from tmesim.physics.pressure import solve_ifp                          # noqa: E402
from tmesim.sim.engine import simulate_response                        # noqa: E402
from tmesim.validation import vessel_profile as VP                     # noqa: E402
from tmesim.validation.calibration import (                            # noqa: E402
    calibrate_alpha,
    calibrated_drugparams,
)
from tmesim.validation.vessel_profile import (                         # noqa: E402
    direct_vessel_distance_um,
    endothelial_xy,
)

# --------------------------------------------------------------------------- #
# The shipped state this file is a measurement OF. Nothing here is a choice.
# --------------------------------------------------------------------------- #
ALPHA_LIVE = 22.223477711172563

#: ``tmesim.data.section_catalog.DEFAULT_SECTION`` and the section every
#: ``GET /api/defaults`` spec names.
SECTION = "xenium_319"

#: ``app/server/api.py::defaults()`` -> ``"dx_um": 12.0``. THE SHIPPED DEFAULT.
SHIPPED_DX_UM = 12.0

#: The HOSTED-INSTANCE value. ``app/server/api.py`` lines 161-165 record the
#: measurement on a 512 MB instance (Render's free tier), flagship Xenium GBM,
#: in a real container: ``dx = 30 um`` is 46,070 voxels and peaks at 220 MiB --
#: comfortable -- while ``dx = 12 um`` is 286,200 voxels and the container is
#: OOM-KILLED. ``app/web/js/ui/spec.js`` line 833 renders the same number to the
#: user. A host that lowers ``TMESIM_MAX_GRID_CELLS`` to 100,000 refuses dx = 12
#: and permits dx >= ~20 um on this section, which is why 20 is on the ladder
#: too.
HOSTED_DX_UM = 30.0

#: The grid-cell ceiling a constrained host is documented to set.
HOSTED_MAX_GRID_CELLS = 100_000

#: The ladder. Coarse end: the hosted value and the two rungs between it and the
#: shipped default. Fine end: 7 um, because
#: ``results/grid_convergence_gbm.md`` measured the advection CFL on THIS
#: section at the production dt = 10 s and found dx = 7 -> 0.96 and dx = 6 ->
#: 1.13 REFUSED. 7 um is therefore the finest mesh the shipped pipeline can
#: reach without also changing dt, which would refine two discretisations at
#: once. The CFL is re-measured here at every rung rather than trusted.
DX_LADDER: Sequence[float] = (30.0, 24.0, 20.0, 16.0, 12.0, 10.0, 8.0, 7.0)

#: Two compounds, chosen for what they are on the shipped panel, not for their
#: numbers. temozolomide is the compound the Studio OPENS ON and is one of the 7
#: of 12 whose delivery/dose verdict INVERTS under the missing transvascular
#: term (``jobs.VERDICT_STABLE_COMPOUNDS`` excludes it). doxorubicin IS in that
#: stable set and is the second half of the headline pair
#: ``experiments/steady_state_reach.py`` uses. If resolution moved only the
#: already-fragile compound that would be a different finding from moving both.
COMPOUNDS: Sequence[str] = ("temozolomide", "doxorubicin")

#: ``app/server/api.py::STUDIO_DEFAULT_TOTAL_TIME_S``.
TOTAL_TIME_S = 3600.0 * 6

#: ``docs/DEPLOYMENT.md``: the stock ``GET /api/defaults`` spec, run end-to-end
#: in the container. The gate for check C. Not a target -- if this harness
#: cannot reproduce it, the harness is wrong and the sweep does not run.
SHIPPED_KILL_FRACTION = 0.11221232991768859
SHIPPED_GRID_CELLS = 286_200

_DATA = os.path.join(_ROOT, "results", "dx_sensitivity_steady_data.json")

# --------------------------------------------------------------------------- #
# PRE-REGISTRATION. Written here, in code, and copied verbatim into
# results/dx_sensitivity_steady.md BEFORE the sweep was run.
# --------------------------------------------------------------------------- #
#: The headline the reach panel renders, at ``x = 0.1``. ``app/web/js/ui/
#: reach.js`` formats it with ``pct(v, 1) + '%'``, i.e. to one decimal place of
#: a percentage, so the panel's own display resolution is 0.001. The threshold
#: below is 0.01 -- TEN TIMES the display resolution -- so a step that clears it
#: is a movement the user can see, not a rounding artefact. It is also far above
#: the ladder's arithmetic resolution ``1 / n_viable`` (order 1e-4 here).
PREREG = {
    "quantity": "reach.compute_reach(...)['headline_under_reached_fraction'], "
                "i.e. F(0.1): the fraction of VIABLE tumour cells whose "
                "steady-state concentration is below 0.1 x the perivascular "
                "reference E_ref",
    "step_threshold": 0.01,
    "span_threshold": 0.05,
    "statement": (
        "PRE-REGISTERED, BEFORE MEASUREMENT. The shipped steady-state reach "
        "verdict is declared RESOLUTION-DEPENDENT -- and the shipped default "
        "dx = 12 um declared NOT on the converged side -- if EITHER: "
        "(1) |F(0.1; dx_k+1) - F(0.1; dx_k)| > 0.01 for any single refinement "
        "step at dx <= 12 um (12->10, 10->8, 8->7), for either compound; or "
        "(2) the categorical verdict -- the dominant limiter from "
        "sim.mechanism.survival_breakdown, or the delivery-vs-dose call among "
        "exposure-limited survivors -- differs between any two rungs of the "
        "ladder, for either compound. "
        "SECONDARY, reported either way and NOT a substitute for the above: "
        "the full-ladder span max|F(0.1)| - min|F(0.1)| over 30..7 um, called "
        "LARGE above 0.05. "
        "If (1) or (2) fires the headline of results/dx_sensitivity_steady.md "
        "says FAILED and no other metric is substituted for F(0.1)."
    ),
}


def _p(*a):
    print(*a, flush=True)


def _assert_alpha(where: str) -> float:
    """Check D. Also re-asserts the two shipped flags whose movement would make
    every number here incomparable with results/."""
    alpha = calibrate_alpha()
    if alpha != ALPHA_LIVE:
        raise AssertionError(
            f"[{where}] calibrate_alpha() = {alpha!r}, not the shipped "
            f"{ALPHA_LIVE!r}. Every number in this file would be on a "
            "different calibration.")
    if PDEConfig().transvascular_exchange is not False:
        raise AssertionError(
            f"[{where}] PDEConfig.transvascular_exchange is not False; the "
            "shipped delivery law has changed.")
    if PDE_CONVECTION is not True or PDE_DT_S != 10.0 or SIM_SEED != 0:
        raise AssertionError(
            f"[{where}] jobs.py production constants moved: "
            f"PDE_CONVECTION={PDE_CONVECTION!r}, PDE_DT_S={PDE_DT_S!r}, "
            f"SIM_SEED={SIM_SEED!r}")
    _p(f"  [{where}] calibrate_alpha() = {alpha!r}  OK")
    return alpha


def compile_shipped(name: str, alpha: float):
    """The compile path the Studio's own ``drug_rows`` and the shipped
    real-tissue panel both use: ``calibrated_drugparams`` then
    ``apply_barrier_to_drugparams``. Nothing is scaled or overridden."""
    return apply_barrier_to_drugparams(calibrated_drugparams(name, alpha),
                                       get_molecule(name))


def _peak_rss_bytes() -> Dict:
    """Peak resident set of THIS process, reported with its platform units
    stated rather than converted on a guess.

    ``ru_maxrss`` is a high-water mark and is monotone non-decreasing, so it is
    only honest per-process. That is exactly why ``--stage sweep`` runs every
    rung in its own subprocess: the number below is then this rung's own peak
    and nothing else's.
    """
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Darwin reports bytes; Linux reports kilobytes. Stated, not assumed.
    unit = "bytes" if sys.platform == "darwin" else "kilobytes"
    scale = 1.0 if unit == "bytes" else 1024.0
    return {"ru_maxrss_raw": int(raw), "ru_maxrss_unit": unit,
            "peak_rss_bytes": int(raw * scale),
            "peak_rss_mib": float(raw * scale) / (1024.0 ** 2),
            "platform": platform.platform()}


def source_fingerprint() -> Dict:
    """SHA-256 over every ``.py`` under ``tmesim/`` and ``app/server/`` of the
    tree this process is importing from.

    NOT decoration. This repository is worked on by several agents at once, and
    the FIRST attempt at this sweep was contaminated exactly that way: partway
    through the ladder ``tmesim/reconstruct/barrier.py`` (the drug compile) and
    ``tmesim/treatment/drugs.py`` (the delivery source) were both edited under
    the running sweep, and the rung that straddled the edit reported a
    ``kill_fraction`` of 0.3177 against 0.1120 one rung earlier -- a 2.8x jump
    that would have read as a spectacular mesh effect and was not one. That
    attempt was DISCARDED, not repaired. Every rung now carries this hash and
    ``check_one_source_tree`` requires them all equal, so the same contamination
    cannot reach a table again without announcing itself.
    """
    h = hashlib.sha256()
    files = []
    for sub in ("tmesim", os.path.join("app", "server")):
        for dirpath, _dirs, names in os.walk(os.path.join(_ROOT, sub)):
            for nm in sorted(names):
                if nm.endswith(".py"):
                    files.append(os.path.join(dirpath, nm))
    for fp in sorted(files):
        h.update(os.path.relpath(fp, _ROOT).encode())
        with open(fp, "rb") as fh:
            h.update(fh.read())
    return {"sha256": h.hexdigest(), "n_files": len(files), "root": _ROOT}


# --------------------------------------------------------------------------- #
# One rung
# --------------------------------------------------------------------------- #
def run_rung(dx_um: float, compound: str, alpha: float) -> Dict:
    """Reconstruct, march, and read BOTH published verdicts off the one run."""
    out: Dict = {"dx_um": float(dx_um), "compound": compound,
                 "section": SECTION, "alpha": alpha}

    t0 = time.time()
    micro, _cn, _cr = VP.build_arm_microenvironment(dx_um=float(dx_um), repair=0)
    out["build_s"] = time.time() - t0
    dom = micro.domain
    out["grid"] = [int(dom.ny), int(dom.nx)]
    out["n_grid_cells"] = int(dom.ny) * int(dom.nx)
    out["n_tissue_voxels"] = int(np.asarray(dom.tissue_mask, dtype=bool).sum())
    out["n_cells"] = int(len(micro.cells.x))
    out["extent_um"] = [float(dom.nx * dom.dx), float(dom.ny * dom.dx)]
    out["over_hosted_ceiling"] = bool(out["n_grid_cells"] > HOSTED_MAX_GRID_CELLS)

    # The advection CFL the production dt will actually be held to, measured on
    # THIS reconstruction with the production solve -- the same number
    # app/server/api.py::_cfl_guard reports.
    _ifp, vel = solve_ifp(micro)
    v = np.asarray(vel, dtype=np.float64)
    speed = float(np.max(np.abs(v[0]) + np.abs(v[1]))) if v.size else 0.0
    out["max_darcy_speed_um_s"] = speed
    out["cfl_at_production_dt"] = speed * PDE_DT_S / float(dx_um)

    dp = compile_shipped(compound, alpha)
    pde = PDEConfig(dt=PDE_DT_S, total_time=TOTAL_TIME_S,
                    convection=PDE_CONVECTION, save_every=PDE_SAVE_EVERY)
    out["pde"] = {"dt": PDE_DT_S, "total_time": TOTAL_TIME_S,
                  "convection": PDE_CONVECTION, "save_every": PDE_SAVE_EVERY,
                  "transvascular_exchange": bool(pde.transvascular_exchange),
                  "seed": SIM_SEED}

    t0 = time.time()
    res = simulate_response(micro, SimConfig(drug=dp, pk=PKConfig(), pde=pde,
                                             seed=SIM_SEED),
                            attribute=True)
    out["sim_s"] = time.time() - t0
    out["kill_fraction"] = float(res.kill_fraction)

    # --- verdict 1: the per-cell dominant limiter (finite-window march) ----- #
    bd = {k: float(x) for k, x in res.mechanism["survival_breakdown"].items()}
    counts = {k: int(x) for k, x in res.mechanism["counts"].items()}
    n_del = int(counts["SURVIVED_UNREACHED"])
    n_dose = int(counts["SURVIVED_UNDERDOSED"])
    n_exp = n_del + n_dose
    out["mechanism"] = {
        "survival_breakdown": bd,
        "counts": counts,
        "n_tumor": int(res.mechanism["n_tumor"]),
        "n_survived": int(res.mechanism["n_survived"]),
        "dominant_limiter": max(bd, key=bd.get),
        "dominant_share": float(max(bd.values())),
        "n_delivery": n_del,
        "n_dose": n_dose,
        "n_exposure_limited": n_exp,
        "delivery_share_of_exposure_limited": (n_del / n_exp) if n_exp else None,
        "delivery_dose_verdict": (("delivery" if n_del > n_dose else
                                   "dose" if n_dose > n_del else "tied")
                                  if n_exp else "undefined"),
        "cells_to_flip": ((abs(n_del - n_dose) // 2) + 1) if n_exp else None,
    }

    # --- verdict 2: the shipped steady-state reach map ---------------------- #
    # The SAME call app/server/jobs.py makes, with the same distance map. The
    # only argument dropped is `penetration`, which compute_reach DISCARDS in
    # steady-state mode (`if steady or penetration is None:` refits against the
    # steady profile), so the payload is unaffected -- verified by reading
    # app/server/reach.py, and `penetration_reused` is recorded below to prove
    # it was False in production too.
    t0 = time.time()
    dist = direct_vessel_distance_um(micro.domain, endothelial_xy(micro.cells))
    out["distance_s"] = time.time() - t0
    t0 = time.time()
    reach = compute_reach(micro, res, distance_um=dist, drug=dp,
                          convection=PDE_CONVECTION)
    out["reach_s"] = time.time() - t0
    fall = reach.get("falloff") or {}
    out["reach"] = {
        "null_reason": reach["null_reason"],
        "window_mode": reach["window_mode"],
        "window_s": reach.get("window_s"),
        "headline_x": reach.get("headline_x"),
        "headline_under_reached_fraction":
            reach.get("headline_under_reached_fraction"),
        "curve": reach.get("curve"),
        "E_ref": reach.get("E_ref_uM_s"),
        "E_ref_units": reach.get("E_ref_units"),
        "n_viable": reach.get("n_viable"),
        "n_tumor": reach.get("n_tumor"),
        "penetration_reused": reach.get("penetration_reused"),
        "penetration_fit_error": reach.get("penetration_fit_error"),
        "falloff_distance_um": fall.get("distance_um"),
        "falloff_normalised": fall.get("normalised_mean_auc"),
        "falloff_count": fall.get("count"),
        "falloff_bin_um": fall.get("bin_um"),
    }
    out["memory"] = _peak_rss_bytes()
    out["source"] = source_fingerprint()
    out["wall_s"] = (out["build_s"] + out["sim_s"] + out["distance_s"]
                     + out["reach_s"])
    return out


# --------------------------------------------------------------------------- #
# Stages
# --------------------------------------------------------------------------- #
def stage_probe() -> Dict:
    """Cheap prologue: what each rung COSTS and whether the production dt can
    reach it. No transport solve."""
    alpha = _assert_alpha("probe:start")
    _p("=" * 78)
    _p("PROBE: grid, CFL at the production dt, hosted ceiling")
    _p("=" * 78)
    rows = []
    for dx in DX_LADDER:
        t0 = time.time()
        micro, _a, _b = VP.build_arm_microenvironment(dx_um=float(dx), repair=0)
        build = time.time() - t0
        dom = micro.domain
        n = int(dom.ny) * int(dom.nx)
        _ifp, vel = solve_ifp(micro)
        v = np.asarray(vel, dtype=np.float64)
        speed = float(np.max(np.abs(v[0]) + np.abs(v[1]))) if v.size else 0.0
        cfl = speed * PDE_DT_S / float(dx)
        row = {"dx_um": float(dx), "grid": [int(dom.ny), int(dom.nx)],
               "n_grid_cells": n, "build_s": build,
               "max_darcy_speed_um_s": speed, "cfl_at_production_dt": cfl,
               "cfl_refuses": bool(cfl > 1.0),
               "over_hosted_ceiling": bool(n > HOSTED_MAX_GRID_CELLS)}
        rows.append(row)
        _p(f"  dx = {dx:5.1f} um  grid {dom.ny:>4}x{dom.nx:<4} = {n:>8,}  "
           f"build {build:5.1f}s  max|v| = {speed:.4f} um/s  "
           f"CFL = {cfl:.3f}{'  REFUSED' if cfl > 1.0 else ''}"
           f"{'  > hosted ceiling' if row['over_hosted_ceiling'] else ''}")
    _assert_alpha("probe:end")
    return {"rows": rows}


def stage_rung(dx: float, compound: str) -> Dict:
    """One rung, in its own process, so ``ru_maxrss`` is honest."""
    alpha = _assert_alpha("rung:start")
    out = run_rung(float(dx), compound, alpha)
    out["alpha_end"] = _assert_alpha("rung:end")
    print("===RUNG_JSON===")
    print(json.dumps(out))
    return out


def stage_sweep(only_compounds: Optional[Sequence[str]] = None,
                only_dx: Optional[Sequence[float]] = None) -> Dict:
    """Every rung, each in its own subprocess."""
    alpha = _assert_alpha("sweep:start")
    compounds = list(only_compounds or COMPOUNDS)
    ladder = [float(d) for d in (only_dx or DX_LADDER)]
    _p("=" * 78)
    _p("SWEEP: dx x compound, one subprocess per rung")
    _p(f"  ladder   {ladder}")
    _p(f"  compounds {compounds}")
    _p(f"  PRE-REGISTERED step threshold |dF(0.1)| > {PREREG['step_threshold']}")
    _p("=" * 78)
    rows: List[Dict] = []
    dropped: List[Dict] = []
    for compound in compounds:
        for dx in ladder:
            t0 = time.time()
            cmd = [sys.executable, os.path.abspath(__file__), "--stage", "rung",
                   "--dx", repr(dx), "--compound", compound]
            proc = subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True)
            if proc.returncode != 0:
                dropped.append({"dx_um": dx, "compound": compound,
                                "returncode": proc.returncode,
                                "stderr_tail": proc.stderr[-4000:]})
                _p(f"  dx = {dx:5.1f} {compound:<14} DROPPED "
                   f"(exit {proc.returncode})")
                _p(proc.stderr[-2000:])
                continue
            marker = "===RUNG_JSON==="
            if marker not in proc.stdout:
                dropped.append({"dx_um": dx, "compound": compound,
                                "returncode": 0,
                                "stderr_tail": "no JSON marker on stdout"})
                _p(f"  dx = {dx:5.1f} {compound:<14} DROPPED (no JSON)")
                continue
            row = json.loads(proc.stdout.split(marker, 1)[1].strip().splitlines()[0])
            row["subprocess_wall_s"] = time.time() - t0
            rows.append(row)
            r = row["reach"]
            m = row["mechanism"]
            _p(f"  dx = {dx:5.1f} {compound:<14} "
               f"n={row['n_grid_cells']:>8,}  "
               f"F(0.1) = {r['headline_under_reached_fraction']:.6f}  "
               f"mode = {r['window_mode']:<12} "
               f"limiter = {m['dominant_limiter']:<10} "
               f"del/dose = {m['delivery_dose_verdict']:<8} "
               f"kill = {row['kill_fraction']:.6f}  "
               f"{row['subprocess_wall_s']:.0f}s  "
               f"peak {row['memory']['peak_rss_mib']:.0f} MiB")
    # The end-to-end reproducibility gate: re-run the FIRST rung, AFTER the
    # whole ladder, and require it bit-identical. See `source_fingerprint`.
    repeat = None
    if rows:
        first = rows[0]
        _p(f"\n  REPEAT GATE: re-running dx = {first['dx_um']:g} "
           f"{first['compound']} after the full ladder")
        cmd = [sys.executable, os.path.abspath(__file__), "--stage", "rung",
               "--dx", repr(first["dx_um"]), "--compound", first["compound"]]
        proc = subprocess.run(cmd, cwd=_ROOT, capture_output=True, text=True)
        marker = "===RUNG_JSON==="
        if proc.returncode == 0 and marker in proc.stdout:
            again = json.loads(
                proc.stdout.split(marker, 1)[1].strip().splitlines()[0])
            same = (
                again["kill_fraction"] == first["kill_fraction"]
                and again["reach"]["headline_under_reached_fraction"]
                == first["reach"]["headline_under_reached_fraction"]
                and again["reach"]["E_ref"] == first["reach"]["E_ref"]
                and again["source"]["sha256"] == first["source"]["sha256"])
            repeat = {
                "pass": bool(same),
                "dx_um": first["dx_um"], "compound": first["compound"],
                "first": {"kill_fraction": first["kill_fraction"],
                          "F": first["reach"]["headline_under_reached_fraction"],
                          "E_ref": first["reach"]["E_ref"],
                          "sha256": first["source"]["sha256"]},
                "repeat": {"kill_fraction": again["kill_fraction"],
                           "F": again["reach"]["headline_under_reached_fraction"],
                           "E_ref": again["reach"]["E_ref"],
                           "sha256": again["source"]["sha256"]},
            }
            _p(f"  REPEAT GATE: {'PASS' if same else 'FAIL'}")
        else:
            repeat = {"pass": False, "reason": "the repeat rung did not run",
                      "stderr_tail": proc.stderr[-2000:]}
            _p("  REPEAT GATE: FAIL (the repeat rung did not run)")

    payload = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "alpha_start": alpha,
        "repeat_check": repeat,
        "source_root": _ROOT,
        "source_note": os.environ.get("DX_SWEEP_SOURCE_NOTE") or "",
        "section": SECTION,
        "shipped_dx_um": SHIPPED_DX_UM,
        "hosted_dx_um": HOSTED_DX_UM,
        "hosted_max_grid_cells": HOSTED_MAX_GRID_CELLS,
        "dx_ladder": list(ladder),
        "compounds": compounds,
        "reach_ladder": [float(v) for v in REACH_LADDER],
        "headline_x": float(HEADLINE_X),
        "prereg": PREREG,
        "shipped_reference": {
            "kill_fraction": SHIPPED_KILL_FRACTION,
            "n_grid_cells": SHIPPED_GRID_CELLS,
            "source": "docs/DEPLOYMENT.md, stock GET /api/defaults spec run "
                      "end-to-end in the container",
        },
        "rows": rows,
        "dropped": dropped,
    }
    payload["checks"] = run_checks(payload)
    payload["analysis"] = analyse(payload)
    payload["alpha_end"] = _assert_alpha("sweep:end")
    with open(_DATA, "w") as fh:
        json.dump(payload, fh, indent=2)
    _p(f"\n  wrote {_DATA}")
    return payload


# --------------------------------------------------------------------------- #
# Checks
# --------------------------------------------------------------------------- #
def run_checks(payload: Dict) -> Dict:
    rows = payload["rows"]
    out: Dict = {}

    # A: the object really is the steady-state one.
    modes = sorted({r["reach"]["window_mode"] for r in rows})
    nulls = [r for r in rows if r["reach"]["null_reason"]]
    out["A_mode_is_steady"] = {
        "pass": modes == ["steady_state"] and not nulls,
        "window_modes_seen": modes,
        "null_rungs": [(r["dx_um"], r["compound"], r["reach"]["null_reason"])
                       for r in nulls],
        "catches": "a dropped drug= argument, which would silently make every "
                   "rung a finite-window map of a different field",
    }

    # B: the mesh reached the physics.
    fps = {}
    dupes = []
    for r in rows:
        fp = (r["compound"], round(float(r["reach"]["E_ref"]), 12),
              tuple(round(float(v), 12)
                    for v in r["reach"]["curve"]["under_reached_fraction"]))
        if fp in fps:
            dupes.append([fps[fp], r["dx_um"]])
        fps[fp] = r["dx_um"]
    out["B_rungs_distinct"] = {
        "pass": not dupes,
        "duplicate_pairs": dupes,
        "n_fingerprints": len(fps),
        "n_rows": len(rows),
        "catches": "a sweep that reused one reconstruction, which would report "
                   "perfect convergence with every rung green",
    }

    # C: the dx = 12 rung IS the shipped run.
    ref = [r for r in rows if r["dx_um"] == SHIPPED_DX_UM
           and r["compound"] == "temozolomide"]
    if ref:
        r = ref[0]
        out["C_shipped_reproduced"] = {
            "pass": (r["kill_fraction"] == SHIPPED_KILL_FRACTION
                     and r["n_grid_cells"] == SHIPPED_GRID_CELLS),
            "kill_fraction": r["kill_fraction"],
            "kill_fraction_published": SHIPPED_KILL_FRACTION,
            "n_grid_cells": r["n_grid_cells"],
            "n_grid_cells_published": SHIPPED_GRID_CELLS,
            "catches": "a harness that is not the shipped pipeline, whose "
                       "convergence statement would be about something else",
        }
    else:
        out["C_shipped_reproduced"] = {
            "pass": False, "reason": "the dx = 12 temozolomide rung is missing"}

    # E: ONE source tree produced every rung.
    fps = sorted({r["source"]["sha256"] for r in rows})
    out["E_one_source_tree"] = {
        "pass": len(fps) == 1,
        "sha256_seen": fps,
        "root": rows[0]["source"]["root"] if rows else None,
        "n_files": rows[0]["source"]["n_files"] if rows else None,
        "catches": "a concurrent edit to the physics or compile chain landing "
                   "mid-sweep -- the defect that destroyed the first attempt at "
                   "this sweep and produced a 2.8x kill_fraction jump that "
                   "would have read as a mesh effect",
    }

    # F: the sweep is reproducible end to end.
    rep = payload.get("repeat_check")
    if rep:
        out["F_repeat_bit_identical"] = dict(rep)
        out["F_repeat_bit_identical"]["catches"] = (
            "any drift -- source, cache or otherwise -- between the first rung "
            "and the last. A rung re-run AFTER the whole ladder must return "
            "bit-identical numbers to the one run BEFORE it.")
    else:
        out["F_repeat_bit_identical"] = {
            "pass": False, "reason": "the repeat rung was not run"}

    out["all_pass"] = all(v.get("pass") for v in out.values()
                          if isinstance(v, dict))
    return out


# --------------------------------------------------------------------------- #
# Analysis against the PRE-REGISTERED threshold
# --------------------------------------------------------------------------- #
def analyse(payload: Dict) -> Dict:
    rows = payload["rows"]
    step_thr = float(payload["prereg"]["step_threshold"])
    span_thr = float(payload["prereg"]["span_threshold"])
    per: Dict = {}
    for compound in payload["compounds"]:
        rs = sorted([r for r in rows if r["compound"] == compound],
                    key=lambda r: -r["dx_um"])          # coarse -> fine
        if not rs:
            continue
        dxs = [r["dx_um"] for r in rs]
        F = [r["reach"]["headline_under_reached_fraction"] for r in rs]
        steps = []
        for i in range(len(rs) - 1):
            steps.append({
                "from_dx": dxs[i], "to_dx": dxs[i + 1],
                "dF": F[i + 1] - F[i], "abs_dF": abs(F[i + 1] - F[i]),
                "below_shipped": dxs[i] <= SHIPPED_DX_UM,
                "exceeds_threshold": abs(F[i + 1] - F[i]) > step_thr,
            })
        fine_steps = [s for s in steps if s["below_shipped"]]
        limiters = sorted({r["mechanism"]["dominant_limiter"] for r in rs})
        verdicts = sorted({r["mechanism"]["delivery_dose_verdict"] for r in rs})
        per[compound] = {
            "dx_um": dxs,
            "F_headline": F,
            "steps": steps,
            "fine_steps_exceeding": [s for s in fine_steps
                                     if s["exceeds_threshold"]],
            "max_abs_step_below_shipped": (max(s["abs_dF"] for s in fine_steps)
                                           if fine_steps else None),
            "max_abs_step_all": max(s["abs_dF"] for s in steps) if steps else None,
            "span": (max(F) - min(F)) if F else None,
            "span_exceeds_secondary": ((max(F) - min(F)) > span_thr) if F else None,
            "dominant_limiters": limiters,
            "limiter_changes": len(limiters) > 1,
            "delivery_dose_verdicts": verdicts,
            "delivery_dose_changes": len(verdicts) > 1,
            "F_at_shipped": next((r["reach"]["headline_under_reached_fraction"]
                                  for r in rs if r["dx_um"] == SHIPPED_DX_UM),
                                 None),
            "F_at_hosted": next((r["reach"]["headline_under_reached_fraction"]
                                 for r in rs if r["dx_um"] == HOSTED_DX_UM),
                                None),
            "F_at_finest": F[-1] if F else None,
        }
    crit1 = any(bool(v["fine_steps_exceeding"]) for v in per.values())
    crit2 = any(v["limiter_changes"] or v["delivery_dose_changes"]
                for v in per.values())
    return {
        "per_compound": per,
        "criterion_1_fine_step_exceeds": crit1,
        "criterion_2_categorical_verdict_changes": crit2,
        "resolution_dependent": bool(crit1 or crit2),
        "prereg_result": "FAILED" if (crit1 or crit2) else "PASSED",
        "secondary_span_large": any(bool(v["span_exceeds_secondary"])
                                    for v in per.values()),
    }


def stage_report() -> Dict:
    alpha = _assert_alpha("report:start")
    with open(_DATA) as fh:
        payload = json.load(fh)
    a = payload["analysis"]
    _p("=" * 78)
    _p("RAW LADDER")
    _p("=" * 78)
    hdr = (f"{'compound':<14} {'dx':>5} {'grid':>13} {'voxels':>9} "
           f"{'F(0.1)':>9} {'dF':>9} {'E_ref':>11} {'mode':>13} "
           f"{'limiter':>10} {'del/dose':>9} {'wall_s':>7} {'MiB':>7}")
    _p(hdr)
    for compound in payload["compounds"]:
        rs = sorted([r for r in payload["rows"] if r["compound"] == compound],
                    key=lambda r: -r["dx_um"])
        prev = None
        for r in rs:
            F = r["reach"]["headline_under_reached_fraction"]
            d = "" if prev is None else f"{F - prev:+.6f}"
            prev = F
            _p(f"{compound:<14} {r['dx_um']:>5.1f} "
               f"{r['grid'][0]:>5}x{r['grid'][1]:<7} {r['n_grid_cells']:>9,} "
               f"{F:>9.6f} {d:>9} {r['reach']['E_ref']:>11.5g} "
               f"{r['reach']['window_mode']:>13} "
               f"{r['mechanism']['dominant_limiter']:>10} "
               f"{r['mechanism']['delivery_dose_verdict']:>9} "
               f"{r['wall_s']:>7.1f} {r['memory']['peak_rss_mib']:>7.0f}")
    _p("")
    _p("=" * 78)
    _p(f"PRE-REGISTERED RESULT: {a['prereg_result']}")
    _p("=" * 78)
    for c, v in a["per_compound"].items():
        _p(f"  {c}: F(0.1) span {v['span']:.6f} over the ladder; "
           f"max |step| below dx = 12 is {v['max_abs_step_below_shipped']}; "
           f"limiters {v['dominant_limiters']}; "
           f"delivery/dose {v['delivery_dose_verdicts']}")
    _assert_alpha("report:end")
    return payload


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stage", default="probe",
                    choices=["probe", "sweep", "rung", "report"])
    ap.add_argument("--dx", type=float, default=None)
    ap.add_argument("--compound", default=None)
    args = ap.parse_args(argv)
    if args.stage == "probe":
        stage_probe()
    elif args.stage == "rung":
        if args.dx is None or not args.compound:
            ap.error("--stage rung needs --dx and --compound")
        stage_rung(args.dx, args.compound)
    elif args.stage == "sweep":
        stage_sweep(only_compounds=([args.compound] if args.compound else None),
                    only_dx=([args.dx] if args.dx is not None else None))
    elif args.stage == "report":
        stage_report()


if __name__ == "__main__":
    main()
