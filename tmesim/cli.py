"""Command-line interface for **tmesim** / **Vivum** — the glioblastoma +
breast-cancer tumor-microenvironment drug-response simulator (REAL data only).

Installed as the ``vivum`` console script (``pip install -e .``), or run as a
module directly — both are equivalent::

    vivum reconstruct --source xenium --dx 12
    python -m tmesim.cli reconstruct --source xenium --dx 12

    vivum compare-drugs --drugs temozolomide,cisplatin,doxorubicin
    vivum studio

See ``docs/QUICKSTART.md`` for a full copy-pasteable walkthrough (install,
first run, and the underlying Python API for notebook use).

Subcommands
-----------
``reconstruct``
    Reconstruct a whole-niche REAL microenvironment from open spatial data
    (10x Xenium single-cell or 10x Visium spot), print the reconstruction report,
    and save a microenvironment figure. No synthetic data anywhere.
``simulate``
    Reconstruct a real microenvironment, compile a REAL calibrated drug
    (calibrated_drugparams -> binding-site barrier), run the mechanistic
    transport+kill simulation, print kill fraction / penetration, save a figure.
``compare-drugs``
    Reconstruct a real microenvironment ONCE and simulate several drugs on the
    SAME reconstruction, printing (and optionally saving) a comparison table.
    Works for any drug in the molecular library — GBM
    (:data:`tmesim.molecular.library.GBM_DRUGS`) by default, or breast
    (:data:`tmesim.molecular.breast_library.BREAST_DRUGS`) via ``--library breast``
    — so a newly added drug is picked up automatically, with no CLI change.
``molecules``
    List the library GBM drugs with structure-derived descriptors (no simulation).
``molecular-run``
    Compile a library molecule to DrugParams (Stokes-Einstein -> Brownian
    hindrance -> BBB/binding) and simulate it on a REAL reconstructed patient.
``acquire``
    Best-effort live download of open glioblastoma datasets into a cache dir.
``studio``
    Launch the tmesim Studio web backend (thin delegate to
    ``python -m app.server``; needs the optional ``pip install -e '.[app]'``
    extra).

Design notes
------------
* Every heavy module (``reconstruct.real_pipeline``, ``sim.engine``,
  ``validation.calibration``, ``molecular.*``, ``viz.plots``, ``data.acquire``) is
  imported **lazily inside each command** so ``--help`` stays fast and never
  crashes if an optional dependency is missing.
* ``--seed`` is global and threaded into the simulation config for determinism.
* The commands use a *reduced* ``PDEConfig.total_time`` so the CLI finishes
  quickly; the full-fidelity window lives in the library defaults.
* The repository root is put on ``sys.path`` unconditionally at startup (see
  ``_ensure_repo_root_on_path``) so ``vivum studio`` can find the ``app``
  package regardless of the caller's current working directory. Commands that
  read or write real data (``data_cache/``, ``results/``) still resolve those
  paths relative to the current working directory, so run ``vivum`` from the
  repository root (exactly as ``python -m tmesim.cli`` already required).
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from typing import List, Optional

# Fast default simulated-exposure window (seconds) for the CLI. The library
# default (see core.config.PDEConfig) is 12 hours, which is far too slow for an
# interactive run; a short window still exercises the whole transport+kill
# pipeline and produces meaningful figures.
_FAST_TOTAL_TIME = 1800.0     # 30 simulated minutes
_FAST_SAVE_EVERY = 200        # snapshot cadence (PDE steps)
_DEFAULT_DRUG = "temozolomide"

# The open breast Xenium section (Janesick et al. 2023, PMID 38114474) used by
# ``experiments/breast_drug_panel.py``. Only used as the ``--sample`` default
# when ``--library breast`` is selected and the caller did not pass ``--sample``
# explicitly — the Xenium loader's own default sample is the GBM section.
_BREAST_SAMPLE_PATH = os.path.join(
    "data_cache", "xenium_breast_janesick", "rep1"
)


def _ensure_repo_root_on_path() -> str:
    """Put the repository root on ``sys.path`` so ``import app`` (Studio) works
    regardless of the caller's current working directory.

    Mirrors ``app/server/__main__.py``'s own ``_ensure_repo_on_path`` so the two
    entry points behave identically. Data/results paths used elsewhere in this
    CLI are still resolved relative to the CWD (as they always have been) —
    this only affects whether the optional ``app`` package can be imported.
    """
    root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    if root not in sys.path:
        sys.path.insert(0, root)
    return root


# --------------------------------------------------------------------------- #
# Config helpers
# --------------------------------------------------------------------------- #
def _build_sim_config(seed: int, total_time: Optional[float] = None):
    """Construct a :class:`~tmesim.core.SimConfig` with a CLI-friendly (reduced)
    PDE exposure window. The drug is set separately by the caller."""
    from tmesim.core import SimConfig

    cfg = SimConfig(seed=int(seed))
    cfg.pde.total_time = float(total_time) if total_time is not None else _FAST_TOTAL_TIME
    cfg.pde.save_every = _FAST_SAVE_EVERY
    return cfg


def _print_report(report: dict) -> None:
    """Pretty-print a reconstruction_report dict."""
    print("\n=== REAL whole-niche microenvironment ===")
    for k, v in report.items():
        if isinstance(v, dict):
            print(f"  {k}:")
            for kk, vv in v.items():
                vv_s = f"{vv:.4f}" if isinstance(vv, float) else str(vv)
                print(f"      {kk:14s} {vv_s}")
        elif isinstance(v, float):
            print(f"  {k:22s} {v:.4f}")
        else:
            print(f"  {k:22s} {v}")


# --------------------------------------------------------------------------- #
# reconstruct
# --------------------------------------------------------------------------- #
def cmd_reconstruct(args: argparse.Namespace) -> int:
    """Reconstruct a REAL microenvironment, print the report, save a figure."""
    from tmesim.reconstruct.real_pipeline import (
        reconstruct_real,
        reconstruction_report,
    )
    from tmesim.viz.plots import plot_microenvironment

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    print(f"[reconstruct] assembling REAL microenvironment "
          f"(source={args.source}, sample={args.sample}, dx={args.dx} um) ...")
    micro = reconstruct_real(source=args.source, sample=args.sample, dx_um=args.dx)
    report = reconstruction_report(micro)
    _print_report(report)

    fig = os.path.join(out_dir, f"reconstruct_{args.source}.png")
    plot_microenvironment(micro, fig)
    print(f"\n[reconstruct] figure written: {os.path.abspath(fig)}")
    return 0


# --------------------------------------------------------------------------- #
# simulate
# --------------------------------------------------------------------------- #
def _compile_real_drug(drug_name: str):
    """Compile a REAL, engine-ready drug: calibrated transport/loss params with
    the binding-site-barrier retention folded in. Returns ``(molecule, DrugParams)``.
    """
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.library import get_molecule
    from tmesim.validation.calibration import calibrated_drugparams

    dp = calibrated_drugparams(drug_name)
    molecule = get_molecule(drug_name)
    dp = apply_barrier_to_drugparams(dp, molecule)
    return molecule, dp


def cmd_simulate(args: argparse.Namespace) -> int:
    """Reconstruct a real microenvironment, compile a real drug, simulate, plot."""
    from tmesim.reconstruct.real_pipeline import (
        reconstruct_real,
        reconstruction_report,
    )
    from tmesim.sim.engine import simulate_response
    from tmesim.viz.plots import plot_microenvironment, plot_simulation

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)
    drug_name = args.drug or _DEFAULT_DRUG

    print(f"[simulate] reconstructing REAL microenvironment "
          f"(source={args.source}, sample={args.sample}, dx={args.dx} um) ...")
    micro = reconstruct_real(source=args.source, sample=args.sample, dx_um=args.dx)
    report = reconstruction_report(micro)
    _print_report(report)

    print(f"\n[simulate] compiling REAL calibrated drug {drug_name!r} "
          "(calibrated_drugparams -> binding-site barrier) ...")
    molecule, dp = _compile_real_drug(drug_name)
    print(f"  drug             : {dp.name}")
    print(f"  effective_D      : {dp.diffusivity:12.3f} um^2/s")
    print(f"  ic50 / hill      : {dp.ic50:.3f} uM / {dp.hill:.3f}")
    print(f"  binding+clearance: {dp.binding:.3e} + {dp.clearance:.3e} 1/s")
    print(f"  bbb_penetrance   : {dp.bbb_penetrance:.3f}")

    cfg = _build_sim_config(args.seed, total_time=args.total_time)
    cfg.drug = dp
    if args.plasma is not None:
        cfg.pk.plasma_peak = float(args.plasma)

    print(f"\n[simulate] running transport+kill simulation "
          f"(total_time={cfg.pde.total_time:.0f}s) ...")
    result = simulate_response(micro, cfg)
    pen = result.extras.get("penetration_index", float("nan"))

    micro_png = os.path.join(out_dir, f"simulate_{args.source}_microenvironment.png")
    sim_png = os.path.join(out_dir, f"simulate_{args.source}_{drug_name}.png")
    plot_microenvironment(micro, micro_png)
    plot_simulation(result, sim_png)

    print("\n=== simulate results ===")
    print(f"  patient_id          : {micro.patient_id}")
    print(f"  drug (molecule)     : {molecule.name}")
    print(f"  predicted_response  : {result.predicted_response:.3f}  (= kill fraction)")
    print(f"  kill_fraction       : {result.kill_fraction:.3f}")
    print(f"  viable_tumor_frac   : {result.viable_tumor_fraction:.3f}")
    print(f"  penetration_index   : {pen:.3f}")
    print("\n[simulate] figures written:")
    print(f"  {os.path.abspath(micro_png)}")
    print(f"  {os.path.abspath(sim_png)}")
    return 0


# --------------------------------------------------------------------------- #
# compare-drugs
# --------------------------------------------------------------------------- #
def _compile_drug_for_comparison(name: str, alpha: float, library: str):
    """Compile ONE drug through the calibrated real-data pipeline, for either
    library.

    ``library='gbm'`` — ``tmesim.molecular.library.get_molecule`` +
    ``validation.calibration.calibrated_drugparams(name, alpha)``: the exact
    calls ``experiments/real_drug_panel.py`` and this module's own
    ``_compile_real_drug`` make.

    ``library='breast'`` — doxorubicin is REUSED VERBATIM from the GBM branch
    (same molecule, same calibrated physics; it is not breast-specific).
    Every other breast drug is compiled from
    ``tmesim.molecular.breast_library.get_molecule_breast`` with a REAL
    GDSC2-breast IC50 (``tmesim.molecular.potency_breast.resolve_ic50``), then
    calibrated through the SAME ``calibrated_drugparams`` call (which accepts
    an already-built ``DrugParams`` and applies the identical alpha-scaled
    perivascular-sink calibration a library-name lookup would). This is the
    identical composition ``experiments/breast_drug_panel.py::compile_breast_drug``
    uses; no calibration code is duplicated here.

    Returns ``(molecule, DrugParams)``. Raises on an unknown / uncompilable
    drug name — the caller decides whether that is fatal for the whole run.
    """
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.library import get_molecule as get_molecule_gbm
    from tmesim.validation.calibration import calibrated_drugparams

    key = name.strip().lower()
    if library == "gbm" or key == "doxorubicin":
        dp = calibrated_drugparams(key, alpha)
        mol = get_molecule_gbm(key)
    else:
        from tmesim.molecular.breast_library import get_molecule_breast
        from tmesim.molecular.compile_drug import molecule_to_drugparams
        from tmesim.molecular.potency_breast import resolve_ic50 as resolve_ic50_breast
        from tmesim.validation.pk_reality import DEFAULT_ECM

        mol = get_molecule_breast(key)
        ic50_info = resolve_ic50_breast(key)
        dp_uncal = molecule_to_drugparams(
            mol, ecm_reference=DEFAULT_ECM, ic50=ic50_info["ic50_uM"]
        )
        dp = calibrated_drugparams(dp_uncal, alpha)
    dp = apply_barrier_to_drugparams(dp, mol)
    return mol, dp


def cmd_compare_drugs(args: argparse.Namespace) -> int:
    """Reconstruct a real microenvironment ONCE and simulate several drugs on
    that SAME reconstruction, so every kill/exposure difference is attributable
    to drug physics rather than to re-reconstructing the tissue. Prints a
    comparison table, sorted by predicted kill fraction, and optionally saves
    it as CSV.

    This is a thin orchestration wrapper: no science is reimplemented here —
    every number comes from the same ``reconstruct_real`` / ``calibrate_alpha``
    / ``calibrated_drugparams`` / ``simulate_response`` calls the underlying
    experiment scripts (``experiments/real_drug_panel.py``,
    ``experiments/breast_drug_panel.py``) use.
    """
    from tmesim.core import SimConfig
    from tmesim.reconstruct.real_pipeline import (
        reconstruct_real,
        reconstruction_report,
    )
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.calibration import calibrate_alpha

    library = args.library
    sample = args.sample
    if library == "breast" and sample is None:
        sample = _BREAST_SAMPLE_PATH
        if not os.path.isdir(sample):
            print(f"[compare-drugs] --library breast needs the Janesick breast "
                  f"Xenium cache at {sample!r}, which is not present. Pass "
                  f"--sample <path> to point at your own cache, or see "
                  f"docs/QUICKSTART.md for how this one was downloaded.")
            return 1

    print(f"[compare-drugs] reconstructing REAL microenvironment "
          f"(source={args.source}, sample={sample}, dx={args.dx} um) — once, "
          f"reused for every drug ...")
    micro = reconstruct_real(source=args.source, sample=sample, dx_um=args.dx)
    report = reconstruction_report(micro)
    _print_report(report)

    if args.drugs:
        drug_names = [d.strip().lower() for d in args.drugs.split(",") if d.strip()]
    else:
        if library == "breast":
            from tmesim.molecular.breast_library import BREAST_DRUGS as _LIB
        else:
            from tmesim.molecular.library import GBM_DRUGS as _LIB
        drug_names = sorted(_LIB)
    if not drug_names:
        print("[compare-drugs] no drugs to compare (empty --drugs list).")
        return 1

    alpha = calibrate_alpha()
    print(f"\n[compare-drugs] alpha = {alpha:.4f} "
          "(perivascular-sink constant, fit once on doxorubicin, shared by "
          "every drug below) ...")

    rows: List[dict] = []
    for name in drug_names:
        print(f"  compiling + simulating {name:16s} ...", end="", flush=True)
        try:
            mol, dp = _compile_drug_for_comparison(name, alpha, library)
            cfg = _build_sim_config(args.seed, total_time=args.total_time)
            cfg.drug = dp
            if args.plasma is not None:
                cfg.pk.plasma_peak = float(args.plasma)
            result = simulate_response(micro, cfg)
        except Exception as exc:
            # A drug library is actively growing (see docs/QUICKSTART.md); one
            # bad or not-yet-wired entry should not sink the whole comparison.
            print(f" SKIPPED ({type(exc).__name__}: {exc})")
            continue
        rows.append(dict(
            drug=name,
            mw=float(mol.get("mw")) if mol.get("mw") is not None else float("nan"),
            ic50_uM=float(dp.ic50),
            bbb_penetrance=float(dp.bbb_penetrance),
            diffusivity_um2_s=float(dp.diffusivity),
            kill_fraction=float(result.kill_fraction),
            penetration_index=float(result.extras.get("penetration_index", float("nan"))),
            mean_auc_uM_s=float(result.extras.get("mean_auc", float("nan"))),
        ))
        print(f" kill={rows[-1]['kill_fraction']:.4f}  "
              f"ic50={rows[-1]['ic50_uM']:.3g} uM  "
              f"mean_auc={rows[-1]['mean_auc_uM_s']:.3g}")

    if not rows:
        print("\n[compare-drugs] every requested drug failed to compile/simulate; "
              "nothing to report.")
        return 1

    rows.sort(key=lambda r: r["kill_fraction"], reverse=True)
    cols = ["drug", "mw", "ic50_uM", "bbb_penetrance", "diffusivity_um2_s",
            "kill_fraction", "penetration_index", "mean_auc_uM_s"]
    print("\n=== drug comparison (sorted by predicted kill_fraction, best first) ===")
    print("PREDICTIONS on a real reconstructed section, not validated clinical "
          "outcomes — see docs/QUICKSTART.md.")
    hdr = f"{'drug':16s} {'MW':>8s} {'IC50(uM)':>9s} {'bbb':>6s} " \
          f"{'D_eff':>8s} {'pen_idx':>8s} {'mean_auc':>10s} {'kill':>7s}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['drug']:16s} {r['mw']:8.0f} {r['ic50_uM']:9.3g} "
              f"{r['bbb_penetrance']:6.3f} {r['diffusivity_um2_s']:8.1f} "
              f"{r['penetration_index']:8.3f} {r['mean_auc_uM_s']:10.3g} "
              f"{r['kill_fraction']:7.4f}")

    if args.csv:
        os.makedirs(os.path.dirname(os.path.abspath(args.csv)) or ".", exist_ok=True)
        with open(args.csv, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=cols)
            writer.writeheader()
            writer.writerows(rows)
        print(f"\n[compare-drugs] comparison table written: {os.path.abspath(args.csv)}")

    return 0


# --------------------------------------------------------------------------- #
# studio
# --------------------------------------------------------------------------- #
def cmd_studio(args: argparse.Namespace) -> int:
    """Launch the tmesim Studio web backend.

    Thin delegate to ``app.server.__main__`` (the existing, already-working
    Studio entry point, normally run as ``python -m app.server``) — nothing
    about the server is reimplemented here. Needs the optional ``app`` extra
    (``pip install -e '.[app]'``); a clear message is printed and this command
    exits 1 rather than crashing with an ImportError traceback if that extra
    is not installed.
    """
    _ensure_repo_root_on_path()
    try:
        from app.server.__main__ import main as _studio_main
    except ImportError as exc:
        print("[studio] the Studio web backend needs the optional 'app' extra "
              f"(fastapi / uvicorn / httpx), which is not installed ({exc}).")
        print("[studio] install it with:  pip install -e '.[app]'")
        print("[studio] then re-run:      vivum studio")
        return 1

    argv: List[str] = ["--host", args.host, "--port", str(args.port),
                       "--log-level", args.log_level]
    if args.reload:
        argv.append("--reload")
    if args.cache_size is not None:
        argv += ["--cache-size", str(args.cache_size)]
    return int(_studio_main(argv))


# --------------------------------------------------------------------------- #
# acquire
# --------------------------------------------------------------------------- #
def cmd_acquire(args: argparse.Namespace) -> int:
    """Lazily import and run the best-effort real-data acquisition."""
    try:
        from tmesim.data.acquire import acquire_all
    except Exception as exc:
        print(f"[acquire] data-acquisition module unavailable: {exc}")
        print("[acquire] nothing downloaded. (This step is optional.)")
        return 0  # optional module missing is not a hard failure

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)
    print(f"[acquire] attempting live download into {out_dir!r} "
          f"(timeout={args.timeout}s) ...")
    try:
        status = acquire_all(cache_dir=out_dir, timeout=args.timeout)
    except Exception as exc:
        print(f"[acquire] acquisition failed cleanly: {exc}")
        return 0

    print("\n=== acquire status ===")
    if isinstance(status, dict):
        for source, info in status.items():
            if isinstance(info, dict):
                st = info.get("status", "?")
                note = info.get("note", "")
                path = info.get("path", "")
                line = f"  {source:20s} : {st}"
                if path:
                    line += f"  -> {path}"
                if note:
                    line += f"  ({note})"
                print(line)
            else:
                print(f"  {source:20s} : {info}")
    else:
        print(f"  {status}")
    return 0


# --------------------------------------------------------------------------- #
# molecular layer: molecules / molecular-run
# --------------------------------------------------------------------------- #
def _import_molecular_library():
    """Lazily import the molecular library, degrading gracefully.

    Returns ``(GBM_DRUGS, get_molecule, all_molecules)`` or ``None`` (with a
    printed notice) if the optional ``tmesim.molecular`` layer is unavailable.
    """
    try:
        from tmesim.molecular.library import GBM_DRUGS, all_molecules, get_molecule
    except Exception as exc:  # module not built / RDKit / import error
        print(f"[molecular] molecular layer unavailable: {exc}")
        print("[molecular] (this is optional; the other commands work without it.)")
        return None
    return GBM_DRUGS, get_molecule, all_molecules


def cmd_molecules(args: argparse.Namespace) -> int:
    """List the library drugs (GBM_DRUGS) with structure-derived key descriptors.

    Prints molecular weight (Da), logP, and a quick physicochemical
    ``bbb_penetrance`` estimate for each molecule. No heavy simulation is run —
    all values are cheap cheminformatics readouts from structure.
    """
    mod = _import_molecular_library()
    if mod is None:
        return 1
    GBM_DRUGS, get_molecule, _all = mod

    # bbb_penetrance is a cheap physicochemical rule (no sim); import lazily and
    # degrade to NaN per-drug if it (or its module) is missing.
    try:
        from tmesim.molecular.permeability import bbb_penetrance as _bbb
    except Exception as exc:
        print(f"[molecules] bbb estimator unavailable ({exc}); showing MW/logP only.")
        _bbb = None

    print(f"[molecules] {len(GBM_DRUGS)} library GBM drugs "
          "(structure-derived descriptors; no simulation):\n")
    header = (f"{'drug':16s} {'MW (Da)':>12s} {'logP':>8s} "
              f"{'bbb_penetrance':>16s}  {'source':<14s}")
    print(header)
    print("-" * len(header))
    for name in sorted(GBM_DRUGS):
        try:
            mol = get_molecule(name)
        except Exception as exc:
            print(f"{name:16s} {'ERR':>12s}  ({type(exc).__name__}: {exc})")
            continue
        mw = mol.get("mw")
        logp = mol.get("logp")
        bbb = float("nan")
        if _bbb is not None:
            try:
                bbb = float(_bbb(mol.descriptors))
            except Exception:
                bbb = float("nan")
        mw_s = f"{mw:12.2f}" if isinstance(mw, (int, float)) else f"{str(mw):>12s}"
        logp_s = f"{logp:8.2f}" if isinstance(logp, (int, float)) else f"{str(logp):>8s}"
        print(f"{name:16s} {mw_s} {logp_s} {bbb:16.3f}  {mol.source:<14s}")
    return 0


def _print_compiled_drugparams(dp) -> None:
    """Print the compiled DrugParams so the molecule->params bridge is visible."""
    print("\n=== compiled DrugParams (molecule -> engine params bridge) ===")
    print(f"  name             : {dp.name}")
    print(f"  effective_D      : {dp.diffusivity:12.3f} um^2/s  (hindered tissue D)")
    print(f"  molecular_weight : {dp.molecular_weight:12.2f} Da")
    print(f"  ic50             : {dp.ic50:12.3f} uM")
    print(f"  hill             : {dp.hill:12.3f}")
    print(f"  binding          : {dp.binding:12.3e} 1/s")
    print(f"  clearance        : {dp.clearance:12.3e} 1/s")
    print(f"  bbb_penetrance   : {dp.bbb_penetrance:12.3f}")


def cmd_molecular_run(args: argparse.Namespace) -> int:
    """Compile a library molecule to DrugParams, then simulate it on a REAL patient.

    The named GBM library molecule is compiled to a ``core.DrugParams`` via the
    molecular layer's ``molecule_to_drugparams`` bridge (Stokes-Einstein free D ->
    Brownian-dynamics hindrance -> permeability/BBB -> binding/clearance). The
    compiled params are printed, then the drug is simulated on a REAL reconstructed
    microenvironment (:func:`reconstruct_real`), writing figures.
    """
    # Optional molecular layer — lazy import so the rest of the CLI still works.
    try:
        from tmesim.molecular.compile_drug import molecule_to_drugparams
        from tmesim.molecular.library import get_molecule
    except Exception as exc:
        print(f"[molecular-run] molecular layer unavailable: {exc}")
        print("[molecular-run] (this is optional; the other commands work without it.)")
        return 1

    from tmesim.reconstruct.real_pipeline import (
        reconstruct_real,
        reconstruction_report,
    )
    from tmesim.sim.engine import simulate_response
    from tmesim.viz.plots import plot_microenvironment, plot_simulation

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    print(f"[molecular-run] compiling library molecule {args.drug!r} "
          "-> DrugParams (Stokes-Einstein -> Brownian hindrance -> BBB/binding) ...")
    try:
        molecule = get_molecule(args.drug)
    except KeyError as exc:
        raise SystemExit(f"[molecular-run] {exc}")
    dp = molecule_to_drugparams(molecule)
    _print_compiled_drugparams(dp)

    print(f"\n[molecular-run] reconstructing REAL microenvironment "
          f"(source={args.source}, sample={args.sample}, dx={args.dx} um) ...")
    micro = reconstruct_real(source=args.source, sample=args.sample, dx_um=args.dx)
    report = reconstruction_report(micro)
    _print_report(report)

    # Build a CLI-fast config and inject the molecularly-compiled drug directly
    # (this IS the molecule->params bridge output).
    cfg = _build_sim_config(args.seed, total_time=args.total_time)
    cfg.drug = dp
    if args.plasma is not None:
        cfg.pk.plasma_peak = float(args.plasma)

    print(f"\n[molecular-run] simulating (drug={cfg.drug.name}, "
          f"total_time={cfg.pde.total_time:.0f}s, "
          f"plasma_peak={cfg.pk.plasma_peak:.1f} uM) ...")
    result = simulate_response(micro, cfg)

    micro_png = os.path.join(out_dir, "molecular_microenvironment.png")
    sim_png = os.path.join(out_dir, f"molecular_{args.drug}.png")
    plot_microenvironment(micro, micro_png)
    plot_simulation(result, sim_png)

    pen = result.extras.get("penetration_index", float("nan"))
    print("\n=== molecular-run results ===")
    print(f"  drug (molecule)     : {molecule.name}")
    print(f"  patient_id          : {micro.patient_id}")
    print(f"  predicted_response  : {result.predicted_response:.3f}  (= kill fraction)")
    print(f"  kill_fraction       : {result.kill_fraction:.3f}")
    print(f"  penetration_index   : {pen:.3f}")
    print("\n[molecular-run] figures written:")
    print(f"  {os.path.abspath(micro_png)}")
    print(f"  {os.path.abspath(sim_png)}")
    return 0


# --------------------------------------------------------------------------- #
# argparse plumbing
# --------------------------------------------------------------------------- #
def _add_real_source_args(p: argparse.ArgumentParser) -> None:
    """Common real-reconstruction args: --source / --sample / --dx."""
    p.add_argument("--source", choices=["xenium", "visium"], default="xenium",
                   help="real spatial data source")
    p.add_argument("--sample", default=None,
                   help="sample id / cache path (default: loader default)")
    p.add_argument("--dx", type=float, default=30.0,
                   help="continuum grid spacing in micrometers")


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="vivum",
        description="Vivum (tmesim) — glioblastoma + breast-cancer tumor-"
                    "microenvironment drug-response simulator (REAL data "
                    "only). Also runnable as `python -m tmesim.cli`.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--seed", type=int, default=0,
                        help="global reproducibility seed")

    sub = parser.add_subparsers(
        dest="command",
        metavar="{reconstruct,simulate,compare-drugs,molecules,molecular-run,"
                 "acquire,studio}")
    sub.required = True

    # reconstruct
    p_rec = sub.add_parser(
        "reconstruct",
        help="reconstruct a REAL microenvironment and print the report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    _add_real_source_args(p_rec)
    p_rec.add_argument("--out", default="results", help="output directory")
    p_rec.set_defaults(func=cmd_reconstruct)

    # simulate
    p_sim = sub.add_parser(
        "simulate",
        help="reconstruct + compile a REAL drug + simulate transport/kill",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    _add_real_source_args(p_sim)
    p_sim.add_argument("--drug", default=None,
                       help=f"library GBM drug name (default: {_DEFAULT_DRUG})")
    p_sim.add_argument("--out", default="results", help="output directory")
    p_sim.add_argument("--plasma", dest="plasma", type=float, default=None,
                       help="override PK plasma_peak (uM) for the run")
    p_sim.add_argument("--total-time", dest="total_time", type=float, default=None,
                       help="PDE exposure window in seconds (default: fast CLI window)")
    p_sim.set_defaults(func=cmd_simulate)

    # compare-drugs
    p_cmp = sub.add_parser(
        "compare-drugs",
        help="reconstruct ONE real section and compare several drugs on it",
        description="Reconstruct a real microenvironment once and simulate "
                    "several drugs on that SAME reconstruction, printing (and "
                    "optionally saving) a comparison table sorted by predicted "
                    "kill fraction. Works for any drug in the molecular "
                    "library, so a newly added drug needs no CLI change.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    _add_real_source_args(p_cmp)
    p_cmp.add_argument("--library", choices=["gbm", "breast"], default="gbm",
                       help="which molecular drug library to draw the default "
                            "drug list from, and how to compile a requested "
                            "drug that is not in GBM_DRUGS")
    p_cmp.add_argument("--drugs", default=None,
                       help="comma-separated drug names to compare (case-"
                            "insensitive). Default: every drug in the selected "
                            "--library (GBM_DRUGS or BREAST_DRUGS)")
    p_cmp.add_argument("--plasma", dest="plasma", type=float, default=None,
                       help="override PK plasma_peak (uM), shared by every drug")
    p_cmp.add_argument("--total-time", dest="total_time", type=float, default=None,
                       help="PDE exposure window in seconds (default: fast CLI window)")
    p_cmp.add_argument("--csv", default=None,
                       help="path to write the comparison table as CSV "
                            "(default: print only, do not save)")
    p_cmp.set_defaults(func=cmd_compare_drugs)

    # studio
    p_studio = sub.add_parser(
        "studio",
        help="launch the tmesim Studio web backend (needs the 'app' extra)",
        description="Launch the tmesim Studio web backend — a thin delegate "
                    "to 'python -m app.server'. Needs 'pip install -e .[app]'.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_studio.add_argument("--host", default="127.0.0.1",
                          help="bind address (local diagnostic workstation "
                               "app, not a network service — pass explicitly "
                               "to expose it beyond localhost)")
    p_studio.add_argument("--port", type=int, default=8000, help="bind port")
    p_studio.add_argument("--reload", action="store_true",
                          help="uvicorn autoreload (development)")
    p_studio.add_argument("--log-level", default="info",
                          choices=("critical", "error", "warning", "info", "debug"))
    p_studio.add_argument("--cache-size", type=int, default=None,
                          help="max cached reconstructions (default 4)")
    p_studio.set_defaults(func=cmd_studio)

    # acquire
    p_acq = sub.add_parser(
        "acquire", help="best-effort live download of open GBM datasets",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_acq.add_argument("--out", default="data_cache", help="cache directory")
    p_acq.add_argument("--timeout", type=int, default=60,
                       help="per-source download timeout (seconds)")
    p_acq.set_defaults(func=cmd_acquire)

    # molecules
    p_mols = sub.add_parser(
        "molecules",
        help="list library GBM drugs with structure-derived descriptors "
             "(MW, logP, bbb_penetrance); no simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_mols.set_defaults(func=cmd_molecules)

    # molecular-run
    p_mrun = sub.add_parser(
        "molecular-run",
        help="compile a library molecule to DrugParams and simulate on a REAL patient",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_mrun.add_argument("--drug", required=True,
                        help="library molecule name from GBM_DRUGS "
                             "(e.g. temozolomide, doxorubicin)")
    _add_real_source_args(p_mrun)
    p_mrun.add_argument("--out", default="results", help="output directory")
    p_mrun.add_argument("--plasma", dest="plasma", type=float, default=None,
                        help="override PK plasma_peak (uM) for the run")
    p_mrun.add_argument("--total-time", dest="total_time", type=float, default=None,
                        help="PDE exposure window in seconds (default: fast CLI window)")
    p_mrun.set_defaults(func=cmd_molecular_run)

    return parser


def main(argv: Optional[List[str]] = None) -> int:
    """CLI entry point. Returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except SystemExit:
        raise
    except KeyboardInterrupt:
        print("\n[tmesim] interrupted.", file=sys.stderr)
        return 130
    except Exception as exc:  # pragma: no cover - top-level guard
        print(f"[tmesim] command '{args.command}' failed: "
              f"{type(exc).__name__}: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
