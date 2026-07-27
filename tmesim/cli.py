"""Command-line interface for **tmesim** — the glioblastoma tumor-microenvironment
drug-response simulator (REAL data only).

Run as a module::

    python -m tmesim.cli reconstruct --source xenium --dx 30
    python -m tmesim.cli simulate --source xenium --drug temozolomide
    python -m tmesim.cli molecules
    python -m tmesim.cli molecular-run --drug doxorubicin --source xenium
    python -m tmesim.cli acquire --out data_cache

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
``molecules``
    List the library GBM drugs with structure-derived descriptors (no simulation).
``molecular-run``
    Compile a library molecule to DrugParams (Stokes-Einstein -> Brownian
    hindrance -> BBB/binding) and simulate it on a REAL reconstructed patient.
``acquire``
    Best-effort live download of open glioblastoma datasets into a cache dir.

Design notes
------------
* Every heavy module (``reconstruct.real_pipeline``, ``sim.engine``,
  ``validation.calibration``, ``molecular.*``, ``viz.plots``, ``data.acquire``) is
  imported **lazily inside each command** so ``--help`` stays fast and never
  crashes if an optional dependency is missing.
* ``--seed`` is global and threaded into the simulation config for determinism.
* The commands use a *reduced* ``PDEConfig.total_time`` so the CLI finishes
  quickly; the full-fidelity window lives in the library defaults.
"""
from __future__ import annotations

import argparse
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
        prog="tmesim",
        description="tmesim — glioblastoma tumor-microenvironment drug-response "
                    "simulator (REAL data only).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--seed", type=int, default=0,
                        help="global reproducibility seed")

    sub = parser.add_subparsers(
        dest="command",
        metavar="{reconstruct,simulate,molecules,molecular-run,acquire}")
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
