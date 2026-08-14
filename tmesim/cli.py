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
``validate``
    Run the EXTERNAL-ANCHOR suite and print a scorecard: every place this model
    is scored against an independent MEASURED number from the literature, with
    a PMID and a PASS / FAIL / NOT TESTABLE verdict per row. ``--quick`` runs
    only the analytic anchors (seconds, no ``data_cache/``); the full run needs
    real tissue and takes minutes. **Exits non-zero if any anchor FAILS**, so it
    is usable as a CI gate. See :mod:`tmesim.validation.scorecard`.
``ingest``
    Point Vivum at YOUR OWN raw 10x Xenium / Visium export directory. Detects
    the platform, validates the bytes, and reports what the section IS (cells,
    extent, panel) and what it CANNOT support — a missing endothelial marker
    means no oxygen or delivery result at all, and on spot data the per-cell
    spatial statistics are refused outright. ``--register`` makes it selectable
    in the Studio. See :mod:`tmesim.data.ingest` and
    ``results/bring_your_own_tissue.md``.
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
import json
import os
import sys
import textwrap
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


def _wrap(text: str, indent: str = "  ") -> str:
    """Wrap a long disclosure for a terminal without truncating a word of it."""
    return textwrap.fill(str(text), width=86, subsequent_indent=indent)


def _delivery_law_note(transvascular_exchange: bool) -> Optional[dict]:
    """The transvascular delivery-law disclosure for a run, or ``None``.

    ``compare-drugs`` prints a RANKING, and the ranking is not law-independent:
    the two source terms disagree on top-1 for GBM. The Studio has carried this
    disclosure on every one of its surfaces for a while; this command — the one
    the README calls the core use case — did not.

    There is exactly ONE implementation of the disclosure
    (``app/server/delivery_law.py``) and it is a RUNTIME read of ``PDEConfig``,
    so it cannot go stale. It is imported rather than restated here precisely so
    that no second copy of the claim can drift from the first. That module
    imports nothing but ``typing`` at module scope, so this does not pull in the
    optional ``[app]`` extra; if the ``app`` package is genuinely absent, the
    caller says so loudly instead of quietly printing nothing.
    """
    _ensure_repo_root_on_path()
    try:
        from app.server.delivery_law import delivery_law
    except Exception:  # noqa: BLE001 — absence is reported, never swallowed
        return None
    return delivery_law(bool(transvascular_exchange))


def _coverage_at_1x_ic50(micro, result, dp) -> float:
    """``therapeutic_coverage`` at 1x IC50 — the headline delivery metric.

    Same code path ``validation.delivery_metrics.compute_all`` uses (it calls
    exactly these three functions); the expensive distance transform behind
    ``penetration_length_fit`` is not needed for coverage and is not run.

    This is what replaced ``penetration_index`` in this command's output.
    ``penetration_index = mean(AUC)/p99(AUC)`` is a scale-invariant SHAPE
    statistic — multiply the exposure field by any constant and it is unchanged
    — so on real tissue it reports the section's vascular geometry, shared by
    every drug (``results/delivery_metrics.md`` §4). It was retired there, and
    this table used to demonstrate exactly why: temozolomide carried the highest
    penetration index of the panel while ranking last.
    """
    from tmesim.validation.delivery_metrics import (
        exposure_window_s,
        therapeutic_coverage,
        tumor_cell_exposure,
    )

    exp = tumor_cell_exposure(micro, result)
    cov = therapeutic_coverage(exp["auc_uM_s"], float(dp.ic50),
                               exposure_window_s(result))
    return float(cov["coverage_1x"])


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
    from tmesim.core import PDEConfig, SimConfig
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

    # --- the two conditions every number below is contingent on ------------- #
    # Printed BEFORE the reconstruction, not after the table, because they are
    # conditions on the whole run rather than footnotes to it.
    #
    # (1) The exposure window. Three surfaces ship three different defaults —
    #     this CLI's fast window, the Studio's, and the library's PDEConfig —
    #     and AUC, coverage and kill all scale with it, so two numbers taken
    #     from two surfaces are not comparable. They are deliberately NOT
    #     silently harmonised here: raising this command's window to the library
    #     default would multiply every quickstart run's cost ~24x, and lowering
    #     the others would change what the Studio has already published. Instead
    #     each surface states the window it used, which is the honest version of
    #     the same fix (``cmd_simulate`` and ``cmd_molecular_run`` already did;
    #     this command did not).
    # (2) The transvascular delivery law, read from PDEConfig at call time.
    probe = _build_sim_config(args.seed, total_time=args.total_time)
    window_s = float(probe.pde.total_time)
    exchange = bool(getattr(probe.pde, "transvascular_exchange", False))
    library_window_s = float(PDEConfig().total_time)
    law = _delivery_law_note(exchange)

    print(f"\n[compare-drugs] EXPOSURE WINDOW = {window_s:.0f} s "
          f"({window_s / 3600.0:.2f} h) — "
          + ("this CLI's fast default (--total-time to change it)"
             if args.total_time is None else "set by --total-time")
          + f". The library default (core.config.PDEConfig.total_time) is "
            f"{library_window_s:.0f} s and the Studio submits a different one "
            f"again; exposure, coverage and kill ALL scale with the window, so "
            f"numbers taken from two surfaces are not comparable unless the "
            f"window matches.")

    if law is not None:
        print("\n[compare-drugs] " + _wrap(law["headline"]))
        print("[compare-drugs] " + _wrap(law["caveat"]))
    else:
        print(f"\n[compare-drugs] WARNING: the delivery-law disclosure could "
              f"not be loaded (app/server/delivery_law.py is not importable "
              f"from here). This run uses "
              f"PDEConfig.transvascular_exchange={exchange!r}. The two source "
              f"terms disagree on the top-1 GBM drug, so the ranking below is "
              f"LAW-DEPENDENT and you are reading it without its disclosure.")

    print(f"\n[compare-drugs] reconstructing REAL microenvironment "
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
            # therapeutic_coverage @1x IC50 REPLACES penetration_index here, in
            # the printed table and in the CSV. See _coverage_at_1x_ic50 for
            # why the index is not carried under any name: it is retired, and
            # this is the export most likely to end up in a slide.
            therapeutic_coverage_1x=_coverage_at_1x_ic50(micro, result, dp),
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
            "kill_fraction", "therapeutic_coverage_1x", "mean_auc_uM_s"]
    print("\n=== drug comparison (sorted by predicted kill_fraction, best first) ===")
    print("PREDICTIONS on a real reconstructed section, not validated clinical "
          "outcomes — see docs/QUICKSTART.md.")
    hdr = f"{'drug':16s} {'MW':>8s} {'IC50(uM)':>9s} {'bbb':>6s} " \
          f"{'D_eff':>8s} {'cov@IC50':>9s} {'mean_auc':>10s} {'kill':>7s}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['drug']:16s} {r['mw']:8.0f} {r['ic50_uM']:9.3g} "
              f"{r['bbb_penetrance']:6.3f} {r['diffusivity_um2_s']:8.1f} "
              f"{r['therapeutic_coverage_1x']:9.4f} {r['mean_auc_uM_s']:10.3g} "
              f"{r['kill_fraction']:7.4f}")
    print(f"\ncov@IC50 = therapeutic_coverage at 1x IC50: the fraction of this "
          f"section's tumour cells whose time-averaged exposure exceeds that "
          f"drug's own IC50 over the {window_s:.0f} s window. It replaces the "
          f"pen_idx column this table used to print — penetration_index is a "
          f"scale-invariant shape statistic that reports the section's vascular "
          f"geometry rather than the drug's delivery, and it is retired "
          f"(results/delivery_metrics.md §4).")
    if law is not None:
        print("\n" + _wrap("RANKING CONDITION — " + law["headline"], indent=""))

    if args.csv:
        prov = _comparison_provenance(
            args=args, sample=sample, library=library, alpha=alpha,
            window_s=window_s, report=report, law=law, exchange=exchange,
            n_drugs=len(rows))
        os.makedirs(os.path.dirname(os.path.abspath(args.csv)) or ".", exist_ok=True)
        with open(args.csv, "w", newline="") as fh:
            # A comment block FIRST, so the provenance cannot be separated from
            # the numbers by copying the file. `pandas.read_csv(..., comment='#')`
            # and every `#`-aware reader skip it; the sidecar below exists for
            # readers that do not.
            for line in _provenance_comment_lines(prov):
                fh.write(f"# {line}\n" if line else "#\n")
            writer = csv.DictWriter(fh, fieldnames=cols)
            writer.writeheader()
            writer.writerows(rows)
        side = os.path.abspath(args.csv) + ".provenance.json"
        with open(side, "w") as fh:
            json.dump(prov, fh, indent=1, sort_keys=True)
        print(f"\n[compare-drugs] comparison table written: {os.path.abspath(args.csv)}")
        print(f"[compare-drugs] it opens with a `#` provenance header — section, "
              f"dx, window, alpha, delivery law, sensitivity mode and NOT FOR "
              f"CLINICAL USE — so read it with "
              f"pandas.read_csv(path, comment='#'). The same block, "
              f"machine-readable, is in: {side}")

    return 0


def _comparison_provenance(*, args, sample, library, alpha, window_s, report,
                           law, exchange, n_drugs) -> dict:
    """Everything a reader needs to know what a comparison CSV actually is.

    The exported CSV used to carry eight columns and nothing else: no section,
    no ``dx``, no exposure window, no ``alpha``, no delivery law, no sensitivity
    mode and no "not for clinical use". It is the only machine-readable export
    in the product and therefore the artefact most likely to become a slide, so
    everything the run knows about its own conditions travels inside it.

    Every value here is READ from the run — the reconstruction report, the
    config actually used, the calibrated alpha — never restated from memory.
    """
    sens = (report or {}).get("sensitivity") or {}
    return {
        "not_for_clinical_use": (
            "NOT FOR CLINICAL USE. Every number in this file is a PREDICTION of "
            "a physical model, not a measurement of a patient. Nothing here has "
            "been validated against a treatment outcome. Do not use it to "
            "select, rank, exclude or dose any therapy for any patient."
        ),
        "generator": "tmesim.cli cmd_compare_drugs (vivum compare-drugs)",
        "section_source": (report or {}).get("source", args.source),
        "sample_path": sample,
        "patient_id": (report or {}).get("patient_id"),
        "library": library,
        "n_drugs": int(n_drugs),
        "dx_um": (report or {}).get("dx_um", args.dx),
        "grid": list((report or {}).get("grid") or []),
        "n_cells": (report or {}).get("n_cells"),
        "n_tumor": (report or {}).get("n_tumor"),
        "exposure_window_s": float(window_s),
        "plasma_peak_uM_override": (float(args.plasma)
                                    if args.plasma is not None else None),
        "seed": int(args.seed),
        "alpha": float(alpha),
        "alpha_source": ("validation.calibration.calibrate_alpha — ONE "
                         "doxorubicin penetration measurement (PMID 16361566), "
                         "shared by every drug in this file"),
        "oxygen_mode": (report or {}).get("oxygen_mode"),
        "state_mode": (report or {}).get("state_mode"),
        "sensitivity_mode": sens.get("mode") or "uniform 1.0 (reconstruction default)",
        "sensitivity_is_upper_bound": bool(sens.get("upper_bound", True)),
        "transvascular_exchange": bool(exchange),
        "delivery_law": (law or {}).get("law"),
        "delivery_law_caveat": (
            (law or {}).get("caveat")
            or (f"UNAVAILABLE — app/server/delivery_law.py was not importable "
                f"when this file was written. The run used "
                f"PDEConfig.transvascular_exchange={bool(exchange)}, and the two "
                f"source terms disagree on the top-1 GBM drug, so this ranking "
                f"is LAW-DEPENDENT.")
        ),
        "retired_metric_note": (
            "penetration_index is NOT in this file. It is mean(AUC)/p99(AUC), a "
            "scale-invariant shape statistic that reports the section's "
            "vascular geometry rather than a drug's delivery, and it is retired "
            "(results/delivery_metrics.md section 4). The column "
            "therapeutic_coverage_1x replaces it."
        ),
        "what_is_defensible": (
            "The ORDER, under the delivery law named above. It is strongest "
            "between compounds that differ substantially in mass and WEAKER "
            "between similar-mass analogues — see diffusivity_model below, "
            "which now separates such compounds rather than being flat, but "
            "reproduces only about a quarter of the measured spread. The "
            "absolute kill_fraction is not defensible: it is a "
            "no-intrinsic-resistance UPPER BOUND and, where it has been "
            "interval-propagated at all, the interval spans ~30x "
            "(results/kill_calibration.md section 5.1)."
        ),
        # The single most load-bearing caveat on this export, because the
        # comparison a real compound series asks for -- rank these analogues --
        # is exactly the one the transport model has been measured to fail.
        "diffusivity_model": (
            "The diffusivity_um2_s column comes from a two-route tissue "
            "transport model — paracellular between cells plus transcellular "
            "through membranes — so it consumes molecular weight, logP AND "
            "hydrogen-bond donor count. Measured blind against 71 compounds in "
            "multicellular layers (Pruijn 2008, doi 10.1071/CH08240) it reaches "
            "R^2 = 0.718 with 99 % of compounds within 3x of measurement and 24 "
            "of 31 exact-mass pairs ordered correctly. The MOLECULAR-WEIGHT-ONLY "
            "model this replaced reached R^2 = 0.022, was flat by construction "
            "at matched mass, and ordered 0 of 31 "
            "(results/diffusivity_falsification.md). WHAT IS STILL NOT TRUE: it "
            "reproduces 8.3x of the measured 29.7x spread — about a QUARTER — "
            "so it under-separates a congeneric series systematically; the whole "
            "test is ONE cell line; and the ABSOLUTE value in this column is not "
            "defensible, the same drug's compiled D_eff missing its measured "
            "multicellular-layer band by 1.5x. Full result: "
            "results/two_route_promotion.md."
        ),
    }


def _provenance_comment_lines(prov: dict) -> List[str]:
    """The provenance block as `#`-comment lines, wrapped for a terminal."""
    out: List[str] = ["NOT FOR CLINICAL USE — see not_for_clinical_use below.",
                      ""]
    for key, value in prov.items():
        if value is None or value == []:
            value = "—"
        text = f"{key}: {value}"
        out.extend(textwrap.wrap(text, width=94, subsequent_indent="    ")
                   or [text])
        out.append("")
    out.extend(textwrap.wrap(
        "--- data below: one row per drug, header row next. Every line above "
        "starts with '#'; read with pandas.read_csv(path, comment='#') or skip "
        "'#' lines. The same block is in <this file>.provenance.json.",
        width=94, subsequent_indent="    "))
    return out


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
# validate — the external-anchor scorecard
# --------------------------------------------------------------------------- #
#: Column widths for the scorecard table. Kept here (not in the scorecard
#: module) because they are presentation, not science.
_V_STATUS_W = 12
_V_SCOPE_W = 8
_V_ANCHOR_W = 34


def _v_rule(char: str = "-", width: int = 100) -> str:
    return char * width


def _print_scorecard(result: dict) -> None:
    """Render :func:`tmesim.validation.scorecard.run_suite`'s result.

    Deliberately prints failures THREE times — in the table, again in a
    dedicated FAILURES block with the full quantity/PMID/conditionals, and once
    more in the summary line. A scorecard that made the model look good by
    omission would be worse than no scorecard at all.
    """
    from tmesim.validation import scorecard as SC

    rows = result["rows"]
    mode = result["mode"]

    print()
    print(_v_rule("="))
    print("VIVUM EXTERNAL-ANCHOR SCORECARD"
          f"    mode={mode}    "
          f"{result['seconds']:.1f}s")
    print(_v_rule("="))
    print("Every row compares a MODEL-PREDICTED number to an INDEPENDENT")
    print("MEASURED number from the cited publication. Rows are never")
    print("suppressed and failures are never softened.")

    # Physics provenance, up front. The breast anchor moves ~10x between the
    # two delivery laws, so a scorecard that does not say which one it ran is
    # ambiguous by an order of magnitude.
    phys = result.get("physics", {})
    if phys:
        print()
        print(f"PHYSICS   {phys['law']}")
        if phys["is_shipping_default"]:
            print("          This IS the law the engine ships with "
                  "(PDEConfig.transvascular_exchange).")
        else:
            print("          !! NOT THE SHIPPING DEFAULT. The engine ships "
                  f"transvascular_exchange="
                  f"{phys['transvascular_exchange_repo_default']}; this run "
                  "overrode it.")
            print("          !! Numbers below are NOT the numbers the shipping "
                  "engine produces.")
        print(f"          {phys['note']}")
    if mode == "quick":
        print()
        print("!! QUICK MODE. Only the analytic anchors ran. Anchors needing real")
        print("!! tissue are DEFERRED below -- DEFERRED IS NOT A PASS, and some of")
        print("!! them FAIL in the full run. A green --quick is not a green suite.")
    print()

    # One block per row rather than a fixed-width table: every field here
    # (the compared quantity, the band, the reason a row has no verdict) is a
    # sentence, and truncating any of them to fit a column is how a scorecard
    # starts lying by abbreviation.
    pad = " " * (_V_STATUS_W + _V_SCOPE_W + 2)
    last_group = None
    for r in rows:
        if r["group"] != last_group:
            last_group = r["group"]
            grp = next((g for g in result["groups"] if g["key"] == last_group), {})
            print(_v_rule("-"))
            print(f"[{last_group}]  {grp.get('title', '')}")
            print(f"    source: {grp.get('source', '')}")
            print(f"    report: {grp.get('report', '')}")
        pm = ", ".join(r["pmids"]) or "n/a"
        print(f"\n{r['status']:{_V_STATUS_W}s} {r['scope']:{_V_SCOPE_W}s} "
              f"{r['label']}   [PMID {pm}]")
        print(f"{pad}compares  : {r['quantity']}")
        print(f"{pad}predicted : {r['predicted']}")
        print(f"{pad}measured  : {r['measured']}")
        if r.get("physics"):
            print(f"{pad}physics   : {r['physics']}")
        if r["bar"]:
            print(f"{pad}bar       : {r['bar']}")
        if r["reason"]:
            print(f"{pad}NO VERDICT: {r['reason']}")
        if r["note"]:
            print(f"{pad}note      : {r['note']}")
        if r["conditional"]:
            print(f"{pad}conditional on: {', '.join(r['conditional'])}")

    # --- failures, repeated in full ---------------------------------------- #
    fails = [r for r in rows if r["status"] == SC.FAIL]
    print()
    print(_v_rule("="))
    if fails:
        print(f"FAILURES ({len(fails)}) -- the model misses these MEASURED values")
        print(_v_rule("="))
        for r in fails:
            print(f"\n  FAIL  {r['group']} / {r['label']}   [PMID "
                  f"{', '.join(r['pmids']) or 'n/a'}]")
            print(f"        compares  : {r['quantity']}")
            print(f"        predicted : {r['predicted']}")
            print(f"        measured  : {r['measured']}")
            if r.get("physics"):
                print(f"        physics   : {r['physics']}")
            print(f"        bar       : {r['bar']}")
            if r["note"]:
                print(f"        note      : {r['note']}")
    else:
        print("FAILURES: none among the rows that carried a verdict.")
        print(_v_rule("="))

    # --- everything without a verdict -------------------------------------- #
    nt = [r for r in rows if r["status"] in (SC.NOT_TESTABLE, SC.DEFERRED)]
    if nt:
        print()
        print(_v_rule("="))
        print(f"NO VERDICT ({len(nt)}) -- NOT a pass. Each row states why.")
        print(_v_rule("="))
        for r in nt:
            print(f"  {r['status']:12s} {r['group']}/{r['label']}: {r['reason']}")

    # --- conditionals ------------------------------------------------------ #
    if result["conditionals"]:
        print()
        print(_v_rule("="))
        print("CONDITIONALS cited by the rows above")
        print(_v_rule("="))
        for key, text in sorted(result["conditionals"].items()):
            print(f"  [{key}] {text}")

    # --- summary ----------------------------------------------------------- #
    c = result["counts"]
    print()
    print(_v_rule("="))
    print(f"SUMMARY   PASS {c[SC.PASS]}   FAIL {c[SC.FAIL]}   "
          f"NOT TESTABLE {c[SC.NOT_TESTABLE]}   DEFERRED {c[SC.DEFERRED]}   "
          f"(of {result['n_rows']} rows)")
    if result["n_scored"]:
        print(f"          {result['n_scored_not_pde']} of the "
              f"{result['n_scored']} scored rows never invoke the transport PDE. "
              f"They cannot be cited as evidence about the solver.")
    if c[SC.FAIL]:
        print(f"VERDICT   {c[SC.FAIL]} anchor(s) FAIL against measured data. "
              f"Exiting 1.")
    else:
        print("VERDICT   no scored anchor fails. Read the NO VERDICT block "
              "before calling this validated.")
    print(_v_rule("="))


def _print_anchor_list() -> None:
    """``--list``: the registry and the deliberate omissions, no computation."""
    from tmesim.validation import scorecard as SC

    print("\nRegistered external anchors (vivum validate):\n")
    print(f"{'KEY':28s} {'TIER':6s} {'SCOPE':9s} {'TISSUE':8s} EST")
    print(_v_rule("-", 92))
    for g in SC.build_registry():
        print(f"{g.key:28s} {g.tier:6s} {g.scope:9s} {g.tissue:8s} {g.est}")
        print(f"{'':28s} {g.title}")
        print(f"{'':28s} source: {g.source}")
        print(f"{'':28s} report: {g.report}")
        if g.needs:
            print(f"{'':28s} needs : {', '.join(g.needs)}")
    print("\nInventoried and DELIBERATELY NOT registered, with the reason:\n")
    for what, why in SC.OUT_OF_SCOPE:
        print(f"  - {what}")
        print(f"      {why}")
    print()


def cmd_validate(args: argparse.Namespace) -> int:
    """Run the external-anchor suite and print the scorecard.

    Thin orchestration only: every comparison is computed by the validation
    module that already owns it (see ``tmesim/validation/scorecard.py``'s
    registry), and every verdict is that module's own declared bar. Nothing
    under ``results/`` is written.

    Returns 1 if any anchor FAILS, so this is usable in CI.
    """
    import json

    from tmesim.validation import scorecard as SC

    if args.list_anchors:
        _print_anchor_list()
        return 0

    only = [k.strip() for k in args.only.split(",")] if args.only else None

    def _progress(g) -> None:
        if args.json:
            return
        if args.quick and g.tier != SC.TIER_QUICK:
            print(f"[validate] {g.key:28s} DEFERRED (--quick)", flush=True)
        else:
            print(f"[validate] {g.key:28s} running ({g.tier}, {g.est or 'n/a'}) ...",
                  flush=True)

    try:
        result = SC.run_suite(quick=args.quick, only=only,
                              exchange=args.exchange, on_progress=_progress)
    except KeyError as exc:
        print(f"[validate] {exc}")
        return 2

    if args.json:
        print(json.dumps(result, indent=2, sort_keys=False, default=str))
    else:
        _print_scorecard(result)
    return int(result["exit_code"])


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
# ingest — bring your own tissue
# --------------------------------------------------------------------------- #
def _ingest_print_detection(det: dict) -> None:
    """Print what the directory was recognised AS, and why."""
    print("\n=== detection ===")
    print(f"  platform      : {det.get('platform', '?')}")
    if det.get("layout"):
        print(f"  layout        : {det['layout']}")
    conf = det.get("confidence")
    if conf is not None:
        print(f"  confidence    : {float(conf):.2f}")
    print(f"  ingestable    : {'yes' if det.get('ingestable') else 'NO'}")
    for issue in det.get("issues") or []:
        print(f"  ! {issue}")


def _ingest_print_sample(sample: dict) -> None:
    """Print what the ingested section IS: size, extent, and what was checked."""
    print("\n=== sample ===")
    print(f"  cells         : {sample.get('n_cells'):,}")
    print(f"  panel genes   : {sample.get('n_genes'):,}")
    ext = sample.get("extent_um") or [float("nan")] * 2
    print(f"  extent        : {ext[0]:,.1f} x {ext[1]:,.1f} um")
    print(f"  coordinates   : {sample.get('coordinate_units', '?')}")
    val = sample.get("validation") or {}
    meas = val.get("measured") or {}
    if "cell_density_per_mm2" in meas:
        print(f"  cell density  : {meas['cell_density_per_mm2']:,.0f} cells/mm^2")
    if "median_counts_per_cell" in meas:
        print(f"  median counts : {meas['median_counts_per_cell']:,.1f} per cell")
    for w in val.get("warnings") or []:
        print(f"  ! WARNING: {w}")


def cmd_ingest(args: argparse.Namespace) -> int:
    """Scan (and optionally register) a raw 10x export directory.

    The default is READ-ONLY: it detects the platform, validates the bytes, and
    prints the capability report — what this specimen can and CANNOT support —
    without touching the section registry. ``--register`` adds it to the table
    the Studio's section picker reads.

    Returns 1 when the directory is not ingestable or registration is refused,
    so this is usable as a gate in a data-prep script. The refusal is always a
    named file, column or bound, never a stack trace.
    """
    import json

    from tmesim.data import ingest as ing

    # --- registry-only modes (no path needed) ------------------------------ #
    if args.list_sections or args.remove:
        reg = _import_ingest_registry()
        if reg is None:
            return 2
        if args.remove:
            out = reg.unregister(args.remove)
            if not out["ok"]:
                print(f"[ingest] REFUSED: {out['error']}")
                return 1
            print(f"[ingest] removed ingested section {args.remove!r}")
            return 0
        sections = reg.ingested_sections()
        if args.json:
            print(json.dumps(sections, indent=2, default=str))
            return 0
        if not sections:
            print("[ingest] no ingested sections. "
                  "Add one with: vivum ingest <dir> --register")
            return 0
        print(f"=== ingested sections ({len(sections)}) ===")
        for sid, rec in sections.items():
            runnable = "runnable" if rec.get("runnable", True) else "INSPECT-ONLY"
            print(f"  {sid:24s} {rec.get('platform', '?'):8s} "
                  f"{rec.get('n_cells_ingested', 0):>9,} cells  {runnable}")
            print(f"  {'':24s} {rec.get('cache_dir', '')}")
        return 0

    if not args.path:
        print("[ingest] a sample directory is required "
              "(or pass --list / --remove ID)")
        return 2

    path = os.path.abspath(os.path.expanduser(args.path))
    if not os.path.isdir(path):
        print(f"[ingest] REFUSED: {path!r} is not a directory. Point this at "
              "the sample folder itself — the one holding cells.csv.gz + "
              "cell_feature_matrix/ (Xenium), or spatial/ + a MEX matrix "
              "directory (Visium).")
        return 1

    # --json must be PURE JSON on stdout: a progress line would make the
    # output unparseable for the data-prep script this flag exists for.
    if not args.json:
        # --register always loads the matrix, whatever --deep says.
        mode = "full load" if (args.deep or args.register) else "panel-only, fast"
        print(f"[ingest] scanning {path!r} ({mode}) ...")

    # --- register: the full load, refused as a whole if anything fails ----- #
    if args.register:
        reg = _import_ingest_registry()
        if reg is None:
            return 2
        out = reg.register(path, section_id=args.id, label=args.label)
        if args.json:
            print(json.dumps(out, indent=2, default=str))
            return 0 if out["ok"] else 1
        if not out["ok"]:
            print(f"[ingest] REFUSED — nothing was registered:\n  {out['error']}")
            return 1
        sec = out["section"]
        _ingest_print_detection({
            "platform": sec.get("platform"), "layout": sec.get("layout"),
            "confidence": sec.get("detect_confidence"), "ingestable": True,
            "issues": sec.get("detect_issues"),
        })
        print("\n=== sample ===")
        print(f"  cells         : {sec.get('n_cells_ingested', 0):,}")
        print(f"  panel genes   : {sec.get('n_genes_ingested', 0):,}")
        ext = sec.get("extent_um_ingested") or [float("nan")] * 2
        print(f"  extent        : {ext[0]:,.1f} x {ext[1]:,.1f} um")
        print()
        print(out["capability_text"])
        print(f"\n[ingest] registered as section id {out['section_id']!r}")
        if not sec.get("runnable", True):
            print("[ingest] NOTE: this section is INSPECT-ONLY — "
                  "reconstruct_real has no route for its platform.")
        else:
            print(f"[ingest] run it:  vivum compare-drugs --source "
                  f"{sec.get('platform')} --sample {sec.get('sample')}")
        print("[ingest] the Studio picks up sections registered by another "
              "process on RESTART, not live.")
        return 0

    # --- scan: read-only ---------------------------------------------------- #
    report = ing.scan_report(path, measure_detection_rates=args.deep)
    if args.json:
        print(json.dumps(report, indent=2, default=str))
        return 0 if report.get("error") is None else 1

    _ingest_print_detection(report.get("detection") or {})
    if report.get("error"):
        print(f"\n[ingest] REFUSED: {report['error']}")
        return 1
    if report.get("sample"):
        _ingest_print_sample(report["sample"])
    print()
    print(report.get("capability_text") or "(no capability report)")
    if not args.deep:
        print("\n[ingest] panel PRESENCE only — no detection rate was measured. "
              "Re-run with --deep to load the matrix and upgrade every "
              "rate-dependent verdict.")
    print("\n[ingest] nothing was registered (this was a read-only scan). "
          "Add --register to make it selectable in the Studio.")
    return 0


def _import_ingest_registry():
    """Lazily import the section registry, degrading with a printed reason."""
    try:
        from app.server import ingest_registry
    except Exception as exc:  # noqa: BLE001 - optional 'app' extra
        print(f"[ingest] the section registry is unavailable: {exc}")
        print("[ingest] registering a section needs the Studio extra: "
              "pip install -e '.[app]'. Scanning (without --register) "
              "works without it.")
        return None
    return ingest_registry


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
        metavar="{reconstruct,simulate,compare-drugs,validate,ingest,"
                 "molecules,molecular-run,acquire,studio}")
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

    # validate
    p_val = sub.add_parser(
        "validate",
        help="score the model against every INDEPENDENT MEASURED literature "
             "anchor and print a PASS/FAIL/NOT-TESTABLE scorecard",
        description="Run the external-anchor suite: every place this model is "
                    "scored against an independent MEASURED number, each row "
                    "citing its PMID and stating exactly what quantity is being "
                    "compared. Failures are reported prominently and the "
                    "command EXITS NON-ZERO if any anchor fails, so it is "
                    "usable as a CI gate. NOT TESTABLE is a distinct status "
                    "from PASS and is never pooled with it.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_val.add_argument("--quick", action="store_true",
                       help="run only the analytic / literature-table anchors "
                            "(seconds, no data_cache/ tissue needed). Anchors "
                            "that need real tissue are reported DEFERRED, which "
                            "is NOT a pass")
    p_val.add_argument("--json", action="store_true",
                       help="emit the scorecard as JSON on stdout instead of "
                            "the human-readable table")
    p_val.add_argument("--only", default=None,
                       help="comma-separated anchor keys to run "
                            "(see --list for the registry)")
    p_val.add_argument("--list", dest="list_anchors", action="store_true",
                       help="list the registered anchors and the validations "
                            "deliberately left out (with reasons), then exit; "
                            "runs nothing")
    p_val.add_argument("--exchange", action="store_true",
                       help="run the breast transvascular anchor with the "
                            "opt-in gradient-driven delivery law "
                            "(PDEConfig.transvascular_exchange) instead of the "
                            "shipping one-way leak. Every number produced this "
                            "way is conditional on the uncalibrated "
                            "PS_SCALE_S_INV -- see the CONDITIONALS block")
    p_val.set_defaults(func=cmd_validate)

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

    # ingest
    p_ing = sub.add_parser(
        "ingest",
        help="scan (and optionally register) YOUR OWN raw 10x Xenium / Visium "
             "export directory",
        description="Point this at a raw 10x export directory and it says what "
                    "the sample IS — platform, cell count, physical extent, "
                    "panel — and, more importantly, what it CANNOT support: "
                    "which cell types have no marker on the panel, whether "
                    "hypoxia is measured at all, whether the vascular field "
                    "has any endothelial marker behind it, and whether "
                    "per-cell spatial statistics are valid on that platform. "
                    "The default is a READ-ONLY scan. --register adds the "
                    "section to the table the Studio's picker reads; "
                    "registration is REFUSED as a whole, never degraded, when "
                    "the sample does not validate.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_ing.add_argument("path", nargs="?", default=None,
                       help="the sample directory (Xenium: cells.csv.gz + "
                            "cell_feature_matrix/; Visium: spatial/ + a MEX "
                            "matrix directory; or a bare MEX triple beside a "
                            "coordinates table)")
    p_ing.add_argument("--deep", action="store_true",
                       help="load the count matrix and measure real DETECTION "
                            "RATES. Without it only panel PRESENCE is read "
                            "(fast), and every rate-dependent verdict is "
                            "reported UNMEASURED rather than guessed")
    p_ing.add_argument("--register", action="store_true",
                       help="add the section to the registry so the Studio and "
                            "the other subcommands can select it (always a "
                            "full load, whatever --deep says)")
    p_ing.add_argument("--id", default=None,
                       help="section id to register under "
                            "(default: slug of the directory name). A BUILT-IN "
                            "section id is always refused")
    p_ing.add_argument("--label", default=None,
                       help="human-readable label shown in the Studio picker")
    p_ing.add_argument("--list", dest="list_sections", action="store_true",
                       help="list the sections already ingested, then exit")
    p_ing.add_argument("--remove", default=None, metavar="ID",
                       help="unregister an ingested section by id "
                            "(built-in sections are never removable)")
    p_ing.add_argument("--json", action="store_true",
                       help="emit the machine-readable report on stdout "
                            "instead of the human-readable one")
    p_ing.set_defaults(func=cmd_ingest)

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
