"""Command-line interface for **tmesim** — the glioblastoma tumor-microenvironment
drug-response simulator.

Run as a module::

    python -m tmesim.cli demo
    python -m tmesim.cli run --n 40 --data synthetic --out results
    python -m tmesim.cli acquire --out data_cache

Subcommands
-----------
``demo``
    Generate ONE synthetic GBM patient, run the full mechanistic simulation, write
    the microenvironment + simulation figures, and print the headline predictions.
``run``
    Build a cohort (synthetic / real / auto), run the honest Stage-0 kill-test
    experiment (sim vs. cheap baselines), and write a markdown report.
``acquire``
    Best-effort live download of open glioblastoma datasets into a cache dir.

Design notes
------------
* The heavy pipeline modules (``sim.engine``, ``validation.experiment``,
  ``viz.report``, ``viz.plots``) and the real-data modules (``data.acquire``,
  ``data.loaders``) are imported **lazily inside each command**. This keeps
  ``--help`` responsive and lets the CLI degrade gracefully with a clear message
  if an optional module (e.g. real-data acquisition) is not present yet.
* ``--seed`` is global and threaded through cohort generation and the simulation
  config so every run is deterministic.
* The ``run`` / ``demo`` paths use a *reduced* ``PDEConfig.total_time`` so the CLI
  finishes quickly; the full-fidelity window lives in the library defaults.
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional

# Fast default simulated-exposure window (seconds) for the CLI. The library
# default (see core.config.PDEConfig) is 12 hours, which is far too slow for an
# interactive demo; a short window still exercises the whole transport+kill
# pipeline and produces meaningful figures.
_FAST_TOTAL_TIME = 1800.0     # 30 simulated minutes
_FAST_SAVE_EVERY = 200        # snapshot cadence (PDE steps)


# --------------------------------------------------------------------------- #
# Config helpers
# --------------------------------------------------------------------------- #
def _build_sim_config(seed: int, total_time: Optional[float] = None,
                      drug_name: Optional[str] = None):
    """Construct a :class:`~tmesim.core.SimConfig` with a CLI-friendly (reduced)
    PDE exposure window.

    Parameters
    ----------
    seed : reproducibility seed threaded into the config.
    total_time : override for ``PDEConfig.total_time`` (seconds); defaults to the
        fast CLI window so runs complete quickly.
    drug_name : optional key into ``DRUG_LIBRARY`` selecting the treating agent.

    Returns
    -------
    SimConfig
    """
    from tmesim.core import DRUG_LIBRARY, SimConfig

    cfg = SimConfig(seed=int(seed))
    cfg.pde.total_time = float(total_time) if total_time is not None else _FAST_TOTAL_TIME
    cfg.pde.save_every = _FAST_SAVE_EVERY
    if drug_name:
        if drug_name not in DRUG_LIBRARY:
            raise SystemExit(
                f"unknown drug {drug_name!r}; choose from {sorted(DRUG_LIBRARY)}")
        cfg.drug = DRUG_LIBRARY[drug_name]
    return cfg


# --------------------------------------------------------------------------- #
# demo
# --------------------------------------------------------------------------- #
def cmd_demo(args: argparse.Namespace) -> int:
    """Generate one synthetic GBM patient, simulate drug response, and plot it.

    Writes ``<out>/demo_microenvironment.png`` and ``<out>/demo_simulation.png``
    and prints ``predicted_response``, ``kill_fraction`` and ``penetration_index``.
    """
    from tmesim.data.synthetic import generate_synthetic_patient
    from tmesim.sim.engine import simulate_response
    from tmesim.viz.plots import plot_microenvironment, plot_simulation

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    print(f"[demo] generating one synthetic GBM patient (seed={args.seed}) ...")
    record = generate_synthetic_patient(patient_id="DEMO-000", seed=args.seed)
    micro = record.micro

    cfg = _build_sim_config(args.seed, total_time=args.total_time,
                            drug_name=args.drug)
    print(f"[demo] simulating drug response (drug={cfg.drug.name}, "
          f"total_time={cfg.pde.total_time:.0f}s) ...")
    result = simulate_response(micro, cfg)

    micro_png = os.path.join(out_dir, "demo_microenvironment.png")
    sim_png = os.path.join(out_dir, "demo_simulation.png")
    plot_microenvironment(micro, micro_png)
    plot_simulation(result, sim_png)

    # penetration_index lives on the TransportResult that drove the kill.
    pen = getattr(getattr(result, "transport", None), "penetration_index", float("nan"))

    print("\n=== demo results ===")
    print(f"  patient_id          : {record.patient_id}")
    print(f"  observed_response   : {record.observed_response:.3f}  (ground truth)")
    print(f"  predicted_response  : {result.predicted_response:.3f}")
    print(f"  kill_fraction       : {result.kill_fraction:.3f}")
    print(f"  penetration_index   : {pen:.3f}")
    print("\n[demo] figures written:")
    print(f"  {os.path.abspath(micro_png)}")
    print(f"  {os.path.abspath(sim_png)}")
    return 0


# --------------------------------------------------------------------------- #
# run
# --------------------------------------------------------------------------- #
def _load_real_cohort(cache_dir: str) -> List:
    """Lazily import and call the real-data loader, degrading gracefully.

    Returns an empty list (never raises) if ``tmesim.data.loaders`` is missing or
    fails, so callers can fall back to synthetic data.
    """
    try:
        from tmesim.data.loaders import load_real_cohort
    except Exception as exc:  # module not built yet / import error
        print(f"[run] real-data loader unavailable ({exc}); no real cohort.")
        return []
    try:
        return load_real_cohort(cache_dir=cache_dir) or []
    except Exception as exc:
        print(f"[run] real-data load failed ({exc}); no real cohort.")
        return []


def _build_cohort(data: str, n: int, seed: int, cache_dir: str) -> List:
    """Build the patient cohort per the ``--data`` mode.

    ``synthetic`` : always synthetic.
    ``real``      : real only (may be empty if no data cached).
    ``auto``      : try real; fall back to synthetic with a printed notice.
    """
    from tmesim.data.synthetic import generate_synthetic_cohort

    if data == "synthetic":
        print(f"[run] building synthetic cohort (n={n}, seed={seed}) ...")
        return generate_synthetic_cohort(n=n, seed=seed)

    if data == "real":
        print(f"[run] loading real cohort from {cache_dir!r} ...")
        records = _load_real_cohort(cache_dir)
        if not records:
            print("[run] no real records available; cohort is empty.")
        return records

    # auto
    print(f"[run] auto mode: attempting real cohort from {cache_dir!r} ...")
    records = _load_real_cohort(cache_dir)
    if records:
        print(f"[run] using {len(records)} real records.")
        return records
    print("[run] NOTICE: real data unavailable/empty — falling back to "
          f"synthetic (n={n}, seed={seed}).")
    return generate_synthetic_cohort(n=n, seed=seed)


def cmd_run(args: argparse.Namespace) -> int:
    """Build a cohort, run the Stage-0 experiment, and write a report."""
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.experiment import run_stage0_experiment
    from tmesim.viz.report import build_report

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    records = _build_cohort(args.data, args.n, args.seed, args.cache_dir)
    if not records:
        print("[run] ERROR: empty cohort — nothing to evaluate. "
              "Try --data synthetic or run `acquire` first.")
        return 1

    cfg = _build_sim_config(args.seed, total_time=args.total_time,
                            drug_name=args.drug)

    print(f"[run] running Stage-0 experiment on {len(records)} patients "
          f"(total_time={cfg.pde.total_time:.0f}s) ...")
    exp = run_stage0_experiment(records, cfg, out_dir=out_dir)

    # Build the report using the first record's simulation + microenvironment as
    # the worked example figure set.
    example = records[0]
    print(f"[run] simulating example patient {example.patient_id} for report ...")
    sim_example = simulate_response(example.micro, cfg)
    report_path = build_report(exp, sim_example, example.micro, out_dir=out_dir)

    print("\n=== run results ===")
    print(f"  cohort size    : {len(records)}")
    print(f"  primary metric : {getattr(exp, 'primary_metric', 'roc_auc')}")
    print(f"  beat_baseline  : {getattr(exp, 'beat_baseline', False)}")
    print(f"  VERDICT        : {getattr(exp, 'verdict', '(none)')}")
    print(f"\n[run] report written: {os.path.abspath(report_path)}")
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
        print("[acquire] nothing downloaded. (This is optional; synthetic data "
              "works without it.)")
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
# molecular layer (optional): molecules / molecular-run
# --------------------------------------------------------------------------- #
def _import_molecular_library():
    """Lazily import the molecular library, degrading gracefully.

    Returns ``(GBM_DRUGS, get_molecule, all_molecules)`` or ``None`` (with a
    printed notice) if the optional ``tmesim.molecular`` layer is unavailable, so
    the rest of the CLI keeps working.
    """
    try:
        from tmesim.molecular.library import GBM_DRUGS, all_molecules, get_molecule
    except Exception as exc:  # module not built / RDKit / import error
        print(f"[molecular] molecular layer unavailable: {exc}")
        print("[molecular] (this is optional; demo/run/acquire work without it.)")
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


def _compile_drug(get_molecule, molecule_to_drugparams, name: str):
    """Resolve a library molecule by ``name`` and compile it to ``DrugParams``.

    Returns ``(molecule, drugparams)``. Raises ``SystemExit`` with a clear
    message on an unknown drug name.
    """
    try:
        molecule = get_molecule(name)
    except KeyError as exc:
        raise SystemExit(f"[molecular-run] {exc}")
    dp = molecule_to_drugparams(molecule)
    return molecule, dp


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
    """Compile a library molecule to DrugParams, then simulate with it.

    The named GBM library molecule is compiled to a ``core.DrugParams`` via the
    molecular layer's ``molecule_to_drugparams`` bridge (Stokes-Einstein free D ->
    Brownian-dynamics hindrance -> permeability/BBB -> binding/clearance). The
    compiled params are printed, then either a single-patient ``demo``-style
    simulation (``--n 1``) or the Stage-0 cohort experiment (``--n > 1``) is run
    with that molecularly-compiled drug, writing figures / a report.
    """
    # Optional molecular layer — lazy import so the rest of the CLI still works.
    try:
        from tmesim.molecular.compile_drug import molecule_to_drugparams
        from tmesim.molecular.library import get_molecule
    except Exception as exc:
        print(f"[molecular-run] molecular layer unavailable: {exc}")
        print("[molecular-run] (this is optional; demo/run/acquire work without it.)")
        return 1

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)

    print(f"[molecular-run] compiling library molecule {args.drug!r} "
          "-> DrugParams (Stokes-Einstein -> Brownian hindrance -> BBB/binding) ...")
    molecule, dp = _compile_drug(get_molecule, molecule_to_drugparams, args.drug)
    _print_compiled_drugparams(dp)

    # Build a CLI-fast config and inject the molecularly-compiled drug directly
    # (bypassing DRUG_LIBRARY — this IS the molecule->params bridge output).
    cfg = _build_sim_config(args.seed, total_time=args.total_time)
    cfg.drug = dp
    if args.plasma is not None:
        cfg.pk.plasma_peak = float(args.plasma)

    if args.n <= 1:
        # ---- single-patient, demo-style simulation ----
        from tmesim.data.synthetic import generate_synthetic_patient
        from tmesim.sim.engine import simulate_response
        from tmesim.viz.plots import plot_microenvironment, plot_simulation

        print(f"\n[molecular-run] single-patient simulation "
              f"(seed={args.seed}, drug={cfg.drug.name}, "
              f"total_time={cfg.pde.total_time:.0f}s, "
              f"plasma_peak={cfg.pk.plasma_peak:.1f} uM) ...")
        record = generate_synthetic_patient(patient_id="MOL-000", seed=args.seed)
        micro = record.micro
        result = simulate_response(micro, cfg)

        micro_png = os.path.join(out_dir, "molecular_microenvironment.png")
        sim_png = os.path.join(out_dir, "molecular_simulation.png")
        plot_microenvironment(micro, micro_png)
        plot_simulation(result, sim_png)

        pen = getattr(getattr(result, "transport", None),
                      "penetration_index", float("nan"))
        print("\n=== molecular-run results (single patient) ===")
        print(f"  drug (molecule)     : {molecule.name}")
        print(f"  patient_id          : {record.patient_id}")
        print(f"  observed_response   : {record.observed_response:.3f}  (ground truth)")
        print(f"  predicted_response  : {result.predicted_response:.3f}")
        print(f"  kill_fraction       : {result.kill_fraction:.3f}")
        print(f"  penetration_index   : {pen:.3f}")
        print("\n[molecular-run] figures written:")
        print(f"  {os.path.abspath(micro_png)}")
        print(f"  {os.path.abspath(sim_png)}")
        return 0

    # ---- cohort Stage-0 experiment ----
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.experiment import run_stage0_experiment
    from tmesim.viz.report import build_report

    records = _build_cohort(args.data, args.n, args.seed, "data_cache")
    if not records:
        print("[molecular-run] ERROR: empty cohort — nothing to evaluate. "
              "Try --data synthetic.")
        return 1

    print(f"\n[molecular-run] Stage-0 experiment on {len(records)} patients "
          f"(drug={cfg.drug.name}, total_time={cfg.pde.total_time:.0f}s) ...")
    exp = run_stage0_experiment(records, cfg, out_dir=out_dir)

    example = records[0]
    print(f"[molecular-run] simulating example patient {example.patient_id} "
          "for report ...")
    sim_example = simulate_response(example.micro, cfg)
    report_path = build_report(exp, sim_example, example.micro, out_dir=out_dir)

    print("\n=== molecular-run results (cohort) ===")
    print(f"  drug (molecule): {molecule.name}")
    print(f"  cohort size    : {len(records)}")
    print(f"  primary metric : {getattr(exp, 'primary_metric', 'roc_auc')}")
    print(f"  beat_baseline  : {getattr(exp, 'beat_baseline', False)}")
    print(f"  VERDICT        : {getattr(exp, 'verdict', '(none)')}")
    print(f"\n[molecular-run] report written: {os.path.abspath(report_path)}")
    return 0


# --------------------------------------------------------------------------- #
# argparse plumbing
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="tmesim",
        description="tmesim — glioblastoma tumor-microenvironment drug-response "
                    "simulator and honest Stage-0 kill-test.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--seed", type=int, default=0,
                        help="global reproducibility seed")

    sub = parser.add_subparsers(
        dest="command",
        metavar="{demo,run,acquire,molecules,molecular-run}")
    sub.required = True

    # demo
    p_demo = sub.add_parser(
        "demo", help="simulate ONE synthetic patient and write figures",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_demo.add_argument("--out", default="results", help="output directory")
    p_demo.add_argument("--drug", default=None,
                        help="drug name from DRUG_LIBRARY (default: config default)")
    p_demo.add_argument("--total-time", dest="total_time", type=float, default=None,
                        help="PDE exposure window in seconds (default: fast CLI window)")
    p_demo.set_defaults(func=cmd_demo)

    # run
    p_run = sub.add_parser(
        "run", help="run the Stage-0 cohort experiment and write a report",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_run.add_argument("--n", type=int, default=40, help="cohort size (synthetic)")
    p_run.add_argument("--data", choices=["synthetic", "real", "auto"],
                       default="synthetic", help="data source")
    p_run.add_argument("--out", default="results", help="output directory")
    p_run.add_argument("--cache-dir", dest="cache_dir", default="data_cache",
                       help="where real data is cached (for real/auto)")
    p_run.add_argument("--drug", default=None,
                       help="drug name from DRUG_LIBRARY (default: config default)")
    p_run.add_argument("--total-time", dest="total_time", type=float, default=None,
                       help="PDE exposure window in seconds (default: fast CLI window)")
    p_run.set_defaults(func=cmd_run)

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
        help="compile a library molecule to DrugParams and simulate with it",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_mrun.add_argument("--drug", required=True,
                        help="library molecule name from GBM_DRUGS "
                             "(e.g. temozolomide, doxorubicin)")
    p_mrun.add_argument("--n", type=int, default=1,
                        help="patients: 1 => demo-style single patient; "
                             ">1 => Stage-0 cohort experiment")
    p_mrun.add_argument("--data", choices=["synthetic", "real", "auto"],
                        default="synthetic", help="data source (cohort mode)")
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
