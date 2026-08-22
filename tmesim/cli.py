"""Command-line interface for **tmesim** / **Codon** — the glioblastoma +
breast-cancer tumor-microenvironment drug-response simulator (REAL data only).

Installed as the ``codon`` console script (``pip install -e .``), or run as a
module directly — both are equivalent::

    codon reconstruct --source xenium --dx 12
    python -m tmesim.cli reconstruct --source xenium --dx 12

    codon simulate --source xenium --dx 12 --drug temozolomide
    codon studio

The pre-rename name ``vivum`` is still registered as an alias console script and
runs exactly the same entry point, so older copy-pasted commands keep working.

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
    ONE tissue, ONE compound — this is the product. Within a single map the
    per-compound scale factors are common to every voxel and cancel exactly, so
    the RELATIVE exposure between regions of that section survives; ordering one
    compound against another does not, and there is no longer a command that
    does it (see the ``compare-drugs`` removal note in ``docs/QUICKSTART.md``,
    ``results/pk_ranking_robustness.md`` and
    ``results/capacity_ranking_robustness.md``).
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
    Point Codon at YOUR OWN raw 10x Xenium / Visium export directory. Detects
    the platform, validates the bytes, and reports what the section IS (cells,
    extent, panel) and what it CANNOT support — a missing endothelial marker
    means no oxygen or delivery result at all, and on spot data the per-cell
    spatial statistics are refused outright. ``--register`` makes it selectable
    in the Studio. See :mod:`tmesim.data.ingest` and
    ``results/bring_your_own_tissue.md``.
``acquire``
    Best-effort live download of open glioblastoma datasets into a cache dir:
    four tabular sources (potency, expression, niche context) plus — this is the
    part that decides whether anything is runnable — one or more *spatial tissue
    sections*. ``--sections`` selects them, ``--list`` prints the plan without
    downloading. See :mod:`tmesim.data.section_catalog`.
``sections``
    What tissue sections exist, which are runnable right now, and for each one
    that is not, **the exact command that fixes it**. This is the command to run
    when anything says a section is unavailable.
``studio``
    Launch the Codon Studio web backend (thin delegate to
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
  ``_ensure_repo_root_on_path``) so ``codon studio`` can find the ``app``
  package regardless of the caller's current working directory. Commands that
  read or write real data (``data_cache/``, ``results/``) still resolve those
  paths relative to the current working directory, so run ``codon`` from the
  repository root (exactly as ``python -m tmesim.cli`` already required).
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

# The open breast Xenium section (Janesick et al. 2023, PMID 38114474) used by
# ``experiments/breast_drug_panel.py``. There is no longer a ``--library``
# switch that resolves it implicitly (it belonged to the removed
# ``compare-drugs`` command), so it is named here as the path to pass to
# ``--sample`` when you want the breast section instead of the GBM default —
# the Xenium loader's own default sample is the GBM section.
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
# studio
# --------------------------------------------------------------------------- #
def cmd_studio(args: argparse.Namespace) -> int:
    """Launch the Codon Studio web backend.

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
        print("[studio] then re-run:      codon studio")
        return 1

    argv: List[str] = ["--host", args.host, "--port", str(args.port),
                       "--log-level", args.log_level]
    if args.reload:
        argv.append("--reload")
    if args.cache_size is not None:
        argv += ["--cache-size", str(args.cache_size)]
    # THREE STATES, and the difference is load-bearing -- `is not None`, never
    # a truthiness test. app.server.__main__ treats them differently:
    #   absent            -> respect TMESIM_STUDIO_API_KEY, and GENERATE one
    #                        off loopback rather than serve unauthenticated
    #   --api-key SECRET  -> use it, overriding the environment
    #   --api-key ''      -> explicitly DISABLE the gate
    # `if args.api_key:` would silently collapse the third into the first and
    # turn "serve with no authentication, I mean it" into "generate a key",
    # which is the opposite of what the operator asked for.
    if args.api_key is not None:
        argv += ["--api-key", args.api_key]
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
    print("CODON EXTERNAL-ANCHOR SCORECARD"
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

    # A first-time user reads a non-zero exit as "my install is broken". It is
    # not: these rows fail on a clean checkout, on every machine, by design.
    # Say so BEFORE the wall of failing rows rather than after it.
    n_fail = result["counts"].get(SC.FAIL, 0)
    if n_fail:
        print()
        print(_v_rule("!"))
        print(f"!! {n_fail} anchor(s) below FAIL, and this command will EXIT 1. "
              "THAT IS THE")
        print("!! DESIGNED OUTCOME OF A CORRECT INSTALL. These are PUBLISHED "
              "LIMITATIONS of")
        print("!! the model -- each one is written up in results/ and named in "
              "README.md's")
        print("!! 'What has failed' section. Nothing about your machine, your "
              "data or your")
        print("!! install produced them, and a build that reported them as "
              "passing would be")
        print("!! the broken one. Do not treat exit 1 here as a setup problem.")
        print(_v_rule("!"))
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
        print("Published limitations, not install problems. Each is written up "
              "in results/.")
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
        print("          Exit 1 is the CORRECT result on a clean install: these "
              "anchors are")
        print("          published limitations of the model, not symptoms of a "
              "bad setup.")
        print("          `codon validate --list` is the registry, with what was "
              "deliberately")
        print("          left out and why. To gate CI on a SUBSET, use "
              "`--only <key,key>`.")
    else:
        print("VERDICT   no scored anchor fails. Read the NO VERDICT block "
              "before calling this validated.")
    print(_v_rule("="))


def _print_anchor_list() -> None:
    """``--list``: the registry and the deliberate omissions, no computation."""
    from tmesim.validation import scorecard as SC

    print("\nRegistered external anchors (codon validate):\n")
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
#: The four non-spatial sources :func:`tmesim.data.acquire.acquire_all` always
#: attempts, and — the part that used to be invisible — what they are FOR.
_TABULAR_SOURCES = (
    ("tcga_gbm_clinical", "TCGA-GBM clinical records (GDC REST API)"),
    ("tcga_gbm_expression", "open STAR gene-count files for those same patients"),
    ("ivy_gap", "Ivy GAP anatomic-structure RNA-seq atlas"),
    ("gdsc_glioma", "GDSC2 glioma cell-line dose-response (drug potency)"),
)


def _acquire_print_plan(section_ids, out_dir: str) -> None:
    """State, BEFORE downloading, exactly what this run will and will not do.

    The failure this replaces: ``acquire`` fetched four tabular sources, said
    ``4/4 ok``, and left every tissue section unavailable, with nothing in the
    output connecting those two facts. A user ran the command the README
    recommends, waited, and was returned to the start.
    """
    from tmesim.data import section_catalog as SC

    print(f"[acquire] plan for {out_dir!r}:")
    print("[acquire]   4 tabular sources — potency, expression, niche context:")
    for name, what in _TABULAR_SOURCES:
        print(f"[acquire]       {name:22s} {what}")
    print("[acquire]     NONE of those four makes a tissue section runnable. "
          "They are")
    print("[acquire]     inputs to a drug, not the tissue a drug is simulated ON.")

    if section_ids:
        total = SC.total_auto_bytes(section_ids)
        print(f"[acquire]   {len(section_ids)} spatial tissue section(s) — "
              f"this is what makes")
        print("[acquire]   `codon reconstruct` / `simulate` / the Studio work:")
        for sid in section_ids:
            rec = SC.SECTION_RECIPES.get(sid, {})
            print(f"[acquire]       {sid:22s} {rec.get('label', '?')}  "
                  f"({SC.human_bytes(rec.get('bytes'))})")
        print(f"[acquire]     total to download: ~{SC.human_bytes(total)} "
              "(already-present files are skipped)")
    else:
        print("[acquire]   NO spatial section (--sections none). This run will "
              "not make")
        print("[acquire]   any section runnable. Drop --sections to fetch "
              f"{SC.DEFAULT_SECTION}.")

    not_fetched = [s for s in SC.SECTION_RECIPES if s not in section_ids]
    if not_fetched:
        print("[acquire]   NOT fetched by this run: " + ", ".join(not_fetched))
        print("[acquire]     `codon sections` prints the exact command for each.")
    print()


def _acquire_print_next_steps(out_dir: str) -> None:
    """Where the user actually is now, and the one command to run next."""
    from tmesim.data import section_catalog as SC

    try:
        st = SC.all_section_status()
    except Exception as exc:  # noqa: BLE001 - a summary must never be the crash
        print(f"[acquire] could not check section availability: {exc}")
        return
    avail = [s for s in st if s["available"]]
    print()
    print(f"=== built-in tissue sections: {len(avail)} of {len(st)} runnable ===")
    for s in st:
        mark = "OK " if s["available"] else "-- "
        print(f"  {mark}{s['id']:22s} {s['label']}")
        if not s["available"] and s["reason"]:
            print(f"       {s['reason']}")
    if avail:
        pick = (SC.DEFAULT_SECTION
                if any(s["id"] == SC.DEFAULT_SECTION for s in avail)
                else avail[0]["id"])
        print()
        print("[acquire] you can run a simulation now:")
        print(f"[acquire]   {SC.run_command(pick)}")
    else:
        print()
        print("[acquire] NO section is runnable yet — nothing can be simulated.")
    print("[acquire] what to do about the unavailable ones: codon sections")
    if os.path.abspath(out_dir) != os.path.abspath("data_cache"):
        print(f"[acquire] NOTE: you downloaded into {out_dir!r}, but the "
              "loaders read ./data_cache/. The")
        print("[acquire] availability list above is about ./data_cache/, "
              "not about your --out directory.")


def cmd_acquire(args: argparse.Namespace) -> int:
    """Lazily import and run the best-effort real-data acquisition."""
    try:
        from tmesim.data.acquire import acquire_all
        from tmesim.data import section_catalog as SC
    except Exception as exc:
        print(f"[acquire] data-acquisition module unavailable: {exc}")
        print("[acquire] nothing downloaded. (This step is optional.)")
        return 0  # optional module missing is not a hard failure

    section_ids, unknown = SC.resolve_section_arg(args.sections)
    if unknown:
        print(f"[acquire] unknown section id(s): {', '.join(unknown)}")
        print(f"[acquire] known: {', '.join(SC.SECTION_RECIPES)} "
              "(or 'all' / 'none')")
        return 2

    if args.list_plan:
        _acquire_print_plan(section_ids, args.out)
        print("[acquire] --list: nothing was downloaded.")
        return 0

    out_dir = args.out
    os.makedirs(out_dir, exist_ok=True)
    _acquire_print_plan(section_ids, out_dir)
    print(f"[acquire] attempting live download into {out_dir!r} "
          f"(timeout={args.timeout}s) ...")
    try:
        status = acquire_all(cache_dir=out_dir, timeout=args.timeout,
                             sections=section_ids)
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
                line = f"  {source:22s} : {st}"
                if path:
                    line += f"  -> {path}"
                if note:
                    line += f"  ({note})"
                print(line)
            else:
                print(f"  {source:22s} : {info}")
    else:
        print(f"  {status}")

    _acquire_print_next_steps(out_dir)
    return 0


# --------------------------------------------------------------------------- #
# sections — what is runnable, and what to do about what is not
# --------------------------------------------------------------------------- #
def cmd_sections(args: argparse.Namespace) -> int:
    """List every built-in section with its availability AND its remedy.

    The point of this command is the remedy. ``GET /api/sections`` and the
    reconstruct error path both say *what is wrong* ("cache directory
    data_cache/xenium_gbm does not exist"); a name of a missing directory is not
    an instruction. Every unavailable row here carries the literal command that
    fixes it — one ``codon acquire --sections`` line where the data is a plain
    HTTP GET, and the actual shell block where it is not.

    Exits 0 when at least one section is runnable, 1 when none is, so a setup
    script can gate on it.
    """
    from tmesim.data import section_catalog as SC

    wanted = args.section_ids or list(SC.SECTION_RECIPES)
    unknown = [s for s in wanted if s not in SC.SECTION_RECIPES]
    if unknown:
        print(f"[sections] unknown section id(s): {', '.join(unknown)}")
        print(f"[sections] known: {', '.join(SC.SECTION_RECIPES)}")
        return 2

    rows = [SC.section_status(sid) for sid in wanted]
    n_avail = sum(1 for r in rows if r["available"])

    if args.json:
        import json

        detail = []
        for r in rows:
            rec = SC.SECTION_RECIPES[r["id"]]
            detail.append({**r,
                           "platform": rec["platform"],
                           "bytes": rec.get("bytes"),
                           "verified": rec.get("verified")})
        print(json.dumps({"available": n_avail, "total": len(rows),
                          "sections": detail}, indent=2))
        return 0 if n_avail else 1

    print()
    print(f"=== built-in tissue sections: {n_avail} of {len(rows)} runnable ===")
    print("A section is a directory of measured tissue under data_cache/, "
          "registered by")
    print("PATH — nothing is uploaded anywhere, and nothing here is "
          "identifiable patient data.")
    for r in rows:
        rec = SC.SECTION_RECIPES[r["id"]]
        print()
        print(f"  [{'RUNNABLE   ' if r['available'] else 'UNAVAILABLE'}] "
              f"{r['id']}  —  {r['label']}")
        print(f"      platform  : {rec['platform']}")
        print(f"      directory : {rec['cache_dir']}")
        print(f"      size      : {SC.human_bytes(rec.get('bytes'))}"
              f"   ({'auto-fetchable' if rec['mode'] == SC.AUTO else 'manual'})")
        if rec.get("note"):
            print(f"      note      : {rec['note']}")
        if args.verbose and rec.get("verified"):
            print(f"      verified  : {rec['verified']}")
        if r["available"]:
            continue
        if r["reason"]:
            print(f"      why not   : {r['reason']}")
        print("      TO FIX IT:")
        for ln in r["remedy"]:
            print(f"        {ln}")

    print()
    if n_avail:
        pick = (SC.DEFAULT_SECTION
                if any(r["id"] == SC.DEFAULT_SECTION and r["available"]
                       for r in rows)
                else next(r["id"] for r in rows if r["available"]))
        print("[sections] at least one section is runnable. Try:")
        print(f"[sections]   {SC.run_command(pick)}")
    else:
        print("[sections] nothing is runnable yet. The shortest path to a first "
              "result:")
        print(f"[sections]   codon acquire --sections {SC.DEFAULT_SECTION}"
              f"   ({SC.human_bytes(SC.total_auto_bytes([SC.DEFAULT_SECTION]))})")
        print(f"[sections]   {SC.run_command(SC.DEFAULT_SECTION)}")
    print("[sections] sections YOU registered are listed separately: "
          "codon ingest --list")
    print("[sections] your own tissue is a first-class option too: "
          "codon ingest <dir> --register")
    return 0 if n_avail else 1


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
                  "Add one with: codon ingest <dir> --register")
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
            print(f"[ingest] run it:  codon simulate --source "
                  f"{sec.get('platform')} --sample {sec.get('sample')} "
                  f"--drug {_DEFAULT_DRUG}")
            print("[ingest] that maps ONE compound on THIS section: where it "
                  "reaches, and how much of the tumour it never reaches. "
                  "Change --drug to map a different compound on the same "
                  "tissue; the two maps are not a ranking of the two "
                  "compounds.")
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
        prog="codon",
        description="Codon (tmesim) — glioblastoma + breast-cancer tumor-"
                    "microenvironment drug-response simulator (REAL data "
                    "only). Also runnable as `python -m tmesim.cli`.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--seed", type=int, default=0,
                        help="global reproducibility seed")

    sub = parser.add_subparsers(
        dest="command",
        metavar="{reconstruct,simulate,validate,ingest,"
                 "molecules,molecular-run,acquire,sections,studio}")
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

    # NOTE (2026-08-21): `compare-drugs` was REMOVED, not renamed. It
    # reconstructed one section and printed several compounds ranked against
    # each other. That cross-compound ordering was measured and does not hold:
    # results/pk_ranking_robustness.md (per-compound plasma exposure spans
    # 113.1x across the retained compounds and the published GBM order first
    # breaks at a spread of 1.53x) and results/capacity_ranking_robustness.md
    # (the breast top-1 pick changes from vinorelbine to doxorubicin under this
    # repo's own structural binding-capacity rule). Within ONE map the
    # per-compound scale factors are common to every voxel and cancel exactly,
    # so `simulate` — one tissue, one compound — is what remains.

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
        help="launch the Codon Studio web backend (needs the 'app' extra)",
        description="Launch the Codon Studio web backend — a thin delegate "
                    "to 'python -m app.server'. Needs 'pip install -e .[app]'.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_studio.add_argument("--host", default="127.0.0.1",
                          help="bind address (local research workstation "
                               "app, not a network service — pass explicitly "
                               "to expose it beyond localhost)")
    p_studio.add_argument("--port", type=int, default=8000, help="bind port")
    p_studio.add_argument("--reload", action="store_true",
                          help="uvicorn autoreload (development)")
    p_studio.add_argument("--log-level", default="info",
                          choices=("critical", "error", "warning", "info", "debug"))
    p_studio.add_argument("--cache-size", type=int, default=None,
                          help="max cached reconstructions (default 4)")
    # Added 2026-08-22. docs/DEPLOYMENT.md §1b's headline command has always
    # been `codon studio --host 0.0.0.0 --api-key <secret>`, which exited 2 --
    # the flag existed on `python -m app.server` and was never forwarded here.
    # default=None is deliberate: it keeps "--api-key ''" distinguishable from
    # not passing the flag at all. See cmd_studio.
    p_studio.add_argument("--api-key", default=None,
                          help="require this key in the X-API-Key header on "
                               "/api/* routes. Omitted and --host non-loopback: "
                               "one is generated and printed once. Pass '' to "
                               "explicitly serve with NO authentication.")
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
        "acquire",
        help="best-effort live download of open GBM datasets, including the "
             "spatial tissue section a first simulation needs",
        description="Downloads four TABULAR sources (TCGA clinical, TCGA "
                    "expression, Ivy GAP, GDSC2 potency) and, by default, the "
                    "SPATIAL SECTION that makes reconstruct/simulate work. The "
                    "tabular four on their own leave every section unavailable "
                    "— run with --list to see the plan before downloading.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_acq.add_argument("--out", default="data_cache", help="cache directory")
    p_acq.add_argument("--timeout", type=int, default=60,
                       help="per-source download timeout (seconds)")
    p_acq.add_argument("--sections", default=None, metavar="IDS",
                       help="comma-separated spatial sections to fetch, or "
                            "'all' (every auto-fetchable one) or 'none'. "
                            "Default: the smallest runnable section "
                            "(xenium_319). `codon sections` lists them")
    p_acq.add_argument("--list", dest="list_plan", action="store_true",
                       help="print what this run WOULD download, then exit "
                            "without downloading anything")
    p_acq.set_defaults(func=cmd_acquire)

    # sections
    p_sec = sub.add_parser(
        "sections",
        help="which tissue sections are runnable, and the exact command to fix "
             "each one that is not",
        description="Every built-in section, whether it is loadable right now, "
                    "and for each unavailable one the literal command that "
                    "makes it available. Exits 1 if no section is runnable.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p_sec.add_argument("section_ids", nargs="*", metavar="ID",
                       help="restrict to these section ids (default: all)")
    p_sec.add_argument("--json", action="store_true",
                       help="machine-readable output")
    p_sec.add_argument("--verbose", action="store_true",
                       help="also print how each recipe's URLs were verified")
    p_sec.set_defaults(func=cmd_sections)

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


def _is_missing_data_refusal(exc: BaseException) -> bool:
    """Is this the loaders' deliberate "I have no data for that" refusal?

    ``IngestError`` is imported lazily: this predicate runs inside the
    top-level handler, where an import failure would replace a useful error
    with a confusing one.
    """
    if isinstance(exc, FileNotFoundError):
        return True
    try:
        from tmesim.data.ingest import IngestError

        return isinstance(exc, IngestError)
    except Exception:  # noqa: BLE001 - fall back to the generic path
        return False


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
        # A missing-data refusal is not a crash. IngestError is raised on
        # purpose, with a sentence written for the person reading it; wrapping
        # that sentence in twenty lines of our own stack frames buries the one
        # part of the output that tells them what to do. Everything else still
        # gets a full traceback, because everything else is a real bug.
        if _is_missing_data_refusal(exc):
            print(f"[{args.command}] REFUSED: {exc}", file=sys.stderr)
            print(f"[{args.command}] this is a missing-data refusal, not a "
                  "crash — nothing was computed and no synthetic data was "
                  "substituted.", file=sys.stderr)
            print(f"[{args.command}] `codon sections` lists every section and "
                  "the exact command that makes each one runnable.",
                  file=sys.stderr)
            return 1
        print(f"[tmesim] command '{args.command}' failed: "
              f"{type(exc).__name__}: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
