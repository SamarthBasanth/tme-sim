"""FLAGSHIP INTEGRATION: run the validated physics engine on a REAL microenvironment.

This assembles the whole-niche REAL glioblastoma microenvironment from the open
10x Xenium single-cell section (40,887 measured cells; threshold-gated marker
typing; REAL vascular field from endothelial cells; REAL ECM field from the
panel's collagen/matrix transcripts; real-perfusion oxygen), compiles a REAL drug
(temozolomide, calibrated perivascular-consumption physics + binding-site barrier),
and runs the transport -> kill engine on that real tissue.

HONESTY (read this): the predicted kill on this section is a PREDICTION with NO
outcome label to check it against. This open Xenium section is an untreated
resection — there is no paired measured drug response. What is real here is the
GEOMETRY (real cells, real vasculature, real matrix) and the DRUG PHYSICS
(structure-compiled transport + a doxorubicin-calibrated perivascular sink). The
kill number is NOT a validated response; it is what the validated physics
predicts on real tissue. Do not read it as an outcome.

RESOLUTION (read this too): the grid-refinement audit in
``results/numerical_rigor.md`` shows that ``kill_fraction`` is grid-converged by
dx <= 20 um, but ``penetration_index`` is NOT converged at dx = 30 um: the
calibrated ~45 um perivascular boundary layer is unresolved on coarse grids, so
penetration is over-estimated by ~15-25%. This experiment therefore runs at the
CONVERGED resolution dx = 12 um (override with ``--dx``).

Run:  ./.venv/bin/python experiments/real_simulation.py
Writes: results/real_simulation.png
"""
from __future__ import annotations

import argparse
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

from tmesim.cells.agents import update_cell_states
from tmesim.cells.pkpd import apply_drug_kill
from tmesim.core import (
    CellState,
    CellType,
    Microenvironment,
    PDEConfig,
    SimConfig,
    TissueFields,
)
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.reconstruct.real_pipeline import (
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "real_simulation.png")

# CONVERGED grid resolution. penetration_index only settles once dx resolves the
# ~45 um perivascular boundary layer (see results/numerical_rigor.md section 2);
# dx = 12 um is inside the converged regime, dx = 30 um is NOT.
DX_UM = 12.0

# CellType -> colour (matches real_xenium.py).
_CT_COLORS = {
    CellType.TUMOR: "#d62728",
    CellType.STROMAL: "#2ca02c",
    CellType.IMMUNE: "#1f77b4",
    CellType.ENDOTHELIAL: "#ff7f0e",
    CellType.NECROTIC: "#7f7f7f",
}


def compile_real_temozolomide():
    """Compile the REAL temozolomide drug: calibrated perivascular-consumption
    physics (alpha fit ONCE on doxorubicin) then the binding-site barrier.

    The barrier is a no-op for a small molecule (temozolomide has no saturable
    high-affinity target), which is the correct, honest behaviour: TMZ is not
    binding-site-barrier-limited. We still route it through the same pipeline an
    antibody would use.
    """
    alpha = calibrate_alpha()  # single fit on doxorubicin (see calibration.py)
    drug = calibrated_drugparams("temozolomide", alpha)
    tmz_mol = get_molecule("temozolomide")
    drug = apply_barrier_to_drugparams(drug, tmz_mol)
    return drug, alpha


def reproduce_kill(micro, result, cfg, drug):
    """Reproduce the engine's EXACT per-cell kill to build a faithful kill map.

    ``simulate_response`` runs on a private copy and does not hand back the final
    cell states, so we replay the deterministic kill on a copy of the same cells
    using the SAME oxygen field (from the result) and the SAME seeded RNG the
    engine used (``default_rng(cfg.seed)``). This yields the identical per-cell
    kill probability and dead set the reported ``kill_fraction`` came from.
    """
    oxygen = result.fields["oxygen"]
    work = Microenvironment(
        domain=micro.domain,
        fields=TissueFields(
            ecm_density=micro.fields.ecm_density,
            vascular_density=micro.fields.vascular_density,
            oxygen=oxygen,
        ),
        cells=micro.cells.copy(),
        patient_id=micro.patient_id,
        metadata=micro.metadata,
    )
    update_cell_states(work, oxygen, cfg.oxygen)
    p_kill = apply_drug_kill(
        work, result.transport, cfg.pd, drug, rng=np.random.default_rng(cfg.seed)
    )
    dead = np.asarray(work.cells.state) == int(CellState.DEAD)
    return p_kill, dead


def main(dx_um: float = DX_UM):
    t0 = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    # --- 1. REAL whole-niche microenvironment (Xenium) ----------------------
    print(f"Reconstructing the REAL Xenium GBM microenvironment (dx={dx_um} um) ...")
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))

    # Per-cell sensitivity for THIS drug. Temozolomide's determinant is MGMT, which
    # is NOT on the 319-gene panel, so this leaves sensitivity at uniform 1.0 and
    # says so in the report rather than inventing heterogeneity.
    # mode='determinant' is stated explicitly: this is a SINGLE-drug run, so the
    # cross-drug comparability problem the 'uniform' default guards against does
    # not arise here. (For temozolomide the two modes coincide anyway.)
    sens_info = apply_drug_sensitivity(micro, "temozolomide", mode="determinant")
    report = reconstruction_report(micro)

    # --- 2. REAL drug: temozolomide (calibrated + barrier) ------------------
    drug, alpha = compile_real_temozolomide()
    print(
        f"\nREAL drug compiled: temozolomide  (calibration alpha={alpha:.1f})\n"
        f"  diffusivity   = {drug.diffusivity:.1f} um^2/s\n"
        f"  ic50          = {drug.ic50:.3g} uM\n"
        f"  binding       = {drug.binding:.3e} /s   clearance = {drug.clearance:.3e} /s\n"
        f"  bbb_penetrance= {drug.bbb_penetrance:.3f}"
    )

    # --- 3. Run the validated physics engine on the real tissue -------------
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)
    cfg = SimConfig(drug=drug, pde=pde, seed=0)
    print("\nRunning simulate_response on the real tissue ...")
    result = simulate_response(micro, cfg)
    p_kill, dead = reproduce_kill(micro, result, cfg, drug)

    # --- 4. Report (stdout) -------------------------------------------------
    print("\n" + "=" * 78)
    print("RECONSTRUCTION REPORT (real Xenium GBM microenvironment)")
    print("=" * 78)
    print(f"  patient_id : {report['patient_id']}   source: {report['source']}")
    print(f"  n_cells    : {report['n_cells']:,}  (n_tumor {report['n_tumor']:,})")
    print(f"  grid       : {report['grid']}  dx={report['dx_um']} um  "
          f"extent={report['extent_um']} um")
    print("  cell-type fractions (CORRECTED, threshold-gated typing):")
    for name, frac in report["celltype_fractions"].items():
        print(f"      {name:12s} {frac * 100:5.2f}%  ({report['celltype_counts'][name]:,})")
    print(f"  mean vascular_density : {report['mean_vascular_density']:.4f}")
    print(f"  mean ecm_density      : {report['mean_ecm_density']:.4f}  "
          f"(REAL matrix transcripts: {', '.join(report['ecm_markers_used'])})")
    print(f"  mean oxygen           : {report['mean_oxygen']:.4f}")
    sb = report["sensitivity"]
    print(f"  per-cell sensitivity  : {sb['model']}")
    print(f"      kind      : {sb['kind']}")
    print(f"      rationale : {sb['rationale']}")
    print(f"      tumor mean={sb['mean_tumor']:.4f} sd={sb['std_tumor']:.4f} "
          f"range=[{sb['min']:.3f}, {sb['max']:.3f}]  uniform={sb['is_uniform']}")
    print("-" * 78)
    print("PREDICTED DRUG RESPONSE (temozolomide on the real tissue)")
    print(f"  kill_fraction     : {result.kill_fraction:.4f}")
    print(f"  penetration_index : {result.extras['penetration_index']:.4f}")
    print(f"  mean_auc          : {result.extras['mean_auc']:.4g} uM*s")
    print(f"  hypoxic_tumor_frac: {result.extras['hypoxic_tumor_fraction']:.4f}")
    print("=" * 78)
    print(
        "HONESTY: this kill_fraction is a PREDICTION with NO outcome label to check\n"
        "against. This open Xenium section is an untreated resection (no paired\n"
        "measured drug response). REAL geometry + REAL drug physics -> a prediction,\n"
        "NOT a validated response. Do not read the number as an outcome.\n"
        "It is also CONDITIONAL on the per-cell sensitivity assumption printed above:\n"
        "temozolomide's determinant is MGMT, which this 319-gene panel does not carry,\n"
        "so sensitivity is uniform 1.0 (no intrinsic resistance modelled) and the kill\n"
        "fraction is the no-intrinsic-resistance UPPER BOUND."
    )
    print("=" * 78)

    # --- 5. FLAGSHIP multi-panel figure -------------------------------------
    W, H = micro.domain.extent_um
    ext = [0, W, H, 0]  # image convention: y increases downward
    cells = micro.cells
    ct = np.asarray(cells.cell_type)

    fig, axes = plt.subplots(2, 3, figsize=(19, 11), constrained_layout=True)

    # (a) real cells by CellType
    ax = axes[0, 0]
    for c in CellType:
        m = ct == int(c)
        if m.any():
            ax.scatter(cells.x[m], cells.y[m], s=2, c=_CT_COLORS[c], linewidths=0,
                       label=f"{c.name} ({int(m.sum()):,})")
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_title("REAL cells by CellType (threshold-gated typing)")
    ax.legend(markerscale=4, fontsize=7, loc="upper right")

    # (b) real vascular density
    ax = axes[0, 1]
    im = ax.imshow(micro.fields.vascular_density, cmap="inferno", extent=ext,
                   aspect="equal", vmin=0, vmax=1)
    ax.set_title("REAL vascular density (endothelial KDE)")
    fig.colorbar(im, ax=ax, shrink=0.8)

    # (c) real oxygen
    ax = axes[0, 2]
    im = ax.imshow(micro.fields.oxygen, cmap="RdYlBu", extent=ext, aspect="equal",
                   vmin=0, vmax=1)
    ax.set_title("REAL oxygen (physics solve on real perfusion)")
    fig.colorbar(im, ax=ax, shrink=0.8)

    # (d) drug penetration (AUC)
    ax = axes[1, 0]
    auc = result.fields["auc"]
    im = ax.imshow(auc, cmap="viridis", extent=ext, aspect="equal")
    ax.set_title(f"Temozolomide penetration (AUC, uM*s)\n"
                 f"penetration_index = {result.extras['penetration_index']:.3f}")
    fig.colorbar(im, ax=ax, shrink=0.8)

    # (e) kill map: tumor cells coloured by predicted per-cell kill probability
    ax = axes[1, 1]
    tum = ct == int(CellType.TUMOR)
    # non-tumor cells as faint grey context
    ax.scatter(cells.x[~tum], cells.y[~tum], s=1, c="#dddddd", linewidths=0)
    sc = ax.scatter(cells.x[tum], cells.y[tum], s=3, c=p_kill[tum], cmap="magma",
                    vmin=0, vmax=1, linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_title(f"Predicted KILL map (tumor cells)\n"
                 f"kill_fraction = {result.kill_fraction:.3f}  "
                 f"({int(dead[tum].sum()):,}/{int(tum.sum()):,} tumor killed)")
    fig.colorbar(sc, ax=ax, shrink=0.8, label="per-cell kill probability")

    # (f) text panel: report + honesty note
    ax = axes[1, 2]
    ax.axis("off")
    frac = report["celltype_fractions"]
    lines = [
        "REAL microenvironment + REAL drug physics",
        "",
        f"source: 10x Xenium GBM  |  {report['n_cells']:,} real cells",
        f"grid {report['grid']}  dx={report['dx_um']:.0f} um (grid-CONVERGED)",
        "",
        "Cell-type fractions (threshold-gated):",
        f"  TUMOR {frac['TUMOR']*100:.1f}%   ENDOTHELIAL {frac['ENDOTHELIAL']*100:.1f}%",
        f"  IMMUNE {frac['IMMUNE']*100:.1f}%   STROMAL {frac['STROMAL']*100:.1f}%",
        f"  (endothelial corrected: 20.6% -> {frac['ENDOTHELIAL']*100:.1f}%,",
        "   consistent with Visium ~9.8%)",
        "",
        f"mean vascular {report['mean_vascular_density']:.3f}  "
        f"ecm {report['mean_ecm_density']:.3f}  o2 {report['mean_oxygen']:.3f}",
        "ECM from REAL collagen/matrix transcripts.",
        "",
        f"Drug: temozolomide (alpha={alpha:.1f}, +barrier)",
        f"  kill_fraction     = {result.kill_fraction:.3f}",
        f"  penetration_index = {result.extras['penetration_index']:.3f}",
        "",
        "alpha fit on the E-FOLD-converted doxorubicin",
        "anchor (40-50um half-distance -> L=d_half/ln2",
        "~65um). Earlier runs used 309.6, which treated",
        "the half-distance as an e-folding length.",
        "",
        "Per-cell sensitivity = UNIFORM 1.0: TMZ's",
        "determinant is MGMT, absent from this panel,",
        "so no intrinsic resistance is modelled and the",
        "kill above is a no-resistance upper bound.",
        "",
        f"Run at CONVERGED dx={report['dx_um']:.0f} um: penetration",
        "needs dx<=12um to resolve the ~45um",
        "perivascular layer (numerical_rigor.md).",
        "Earlier dx=30um runs over-stated it.",
        "",
        "PREDICTION ONLY -- no outcome label.",
        "Untreated resection: no paired measured",
        "response. Real geometry + real drug",
        "physics, NOT a validated response.",
    ]
    ax.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", fontsize=10.5,
            family="monospace", transform=ax.transAxes)

    fig.suptitle(
        "FLAGSHIP: validated drug-transport physics on a REAL glioblastoma "
        "microenvironment (temozolomide)  —  PREDICTION, not a validated outcome",
        fontsize=14,
    )
    fig.savefig(_OUT_PNG, dpi=130)
    print(f"\nflagship figure written: {_OUT_PNG}")
    print(f"[runtime {time.time() - t0:.1f} s]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dx", type=float, default=DX_UM,
                    help="grid spacing in um (default %(default)s, grid-converged)")
    main(dx_um=ap.parse_args().dx)
