"""REAL DRUG PANEL on the REAL Xenium glioblastoma microenvironment.

This is the real-geometry version of the drug comparison (no synthetic cohort):
the FULL panel [temozolomide, lomustine, cisplatin, methotrexate, doxorubicin] is
simulated on ONE and the SAME real tumor -- the open 10x Xenium GBM section
(40,887 measured cells, threshold-gated marker typing, REAL vascular field from
endothelial cells, REAL ECM field from the panel's collagen/matrix transcripts,
real-perfusion oxygen). The real microenvironment is reconstructed ONCE and reused
for all five drugs, so every kill/penetration difference is attributable to the
drug's own physics (size -> mesh diffusivity, real GDSC IC50 potency, structure ->
perivascular sink, binding-site barrier, BBB penetrance) on identical tissue.

Each drug is compiled through the SAME pipeline the flagship single-drug run uses:
  calibrate_alpha()  (single perivascular-sink constant, fit ONCE on doxorubicin)
    -> calibrated_drugparams(name, alpha)   (real GDSC IC50 + calibrated loss)
    -> apply_barrier_to_drugparams(dp, molecule)   (saturable binding-site barrier)

PER-CELL SENSITIVITY — TWO MODES, BOTH REPORTED (read this before ranking anything)
-----------------------------------------------------------------------------------
``reconstruct.real_pipeline.apply_drug_sensitivity`` can set the per-cell
``sensitivity`` array two ways, and the choice CHANGES THE RANKING:

  * ``mode='uniform'`` (THE DEFAULT, and the headline ranking here) -- every drug
    gets sensitivity 1.0 for every cell. The factor is identical across drugs, so
    it cancels out of the comparison and the ranking isolates DELIVERY x POTENCY
    with no sensitivity confound. Absolute kills are then an equal
    no-intrinsic-resistance upper bound for every drug.
  * ``mode='determinant'`` -- each drug uses its own measured panel determinant
    where the 319-gene panel carries one (TOP2A for doxorubicin, an MKI67
    proliferation proxy for cisplatin / lomustine / methotrexate; both land at
    mean ~0.585 on this section) and uniform 1.0 where it does not. TEMOZOLOMIDE
    HAS NO DETERMINANT ON THIS PANEL (its is MGMT), so in this mode TMZ alone
    keeps sensitivity 1.0 while its competitors are scaled down ~1.7x. Any kill
    TMZ gains in this mode is an ARTIFACT OF AN UNMEASURED RESISTANCE GENE, not
    evidence that it works better; every such drug is flagged UPPER BOUND below.

Both tables are printed and written. Rank with the uniform table.

HONESTY (read this): the predicted kills here are PREDICTIONS with NO outcome label
to validate against. This open Xenium section is an UNTREATED resection -- there is
no paired measured drug response for any of these five drugs. What is real is the
GEOMETRY (real cells, real vasculature, real matrix, real oxygen), the DRUG PHYSICS
(structure-compiled transport + a doxorubicin-calibrated perivascular sink,
revalidated out-of-sample on mitoxantrone), and the POTENCY (real GDSC/literature
IC50). The kill numbers are what the validated physics predicts on real tissue --
NOT validated outcomes. Do not read the ranking as clinical efficacy. The section
is also deeply hypoxic (~94% of tumor cells hypoxic), which mechanistically
suppresses kill for EVERY drug and compresses the panel.

RESOLUTION (read this too): per the grid-refinement audit in
``results/numerical_rigor.md``, ``kill_fraction`` is grid-converged by dx <= 20 um
but ``penetration_index`` is NOT converged at dx = 30 um -- the calibrated ~45 um
perivascular boundary layer is unresolved on coarse grids, so penetration is
over-estimated by ~15-25%. This panel therefore runs at the CONVERGED resolution
dx = 12 um (override with ``--dx``).

Run:    ./.venv/bin/python experiments/real_drug_panel.py
Writes: results/real_drug_panel.png
        results/real_drug_panel.md
"""
from __future__ import annotations

import argparse
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from tmesim.core import PDEConfig, SimConfig
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.physics.diffusion import ALPHA_TORTUOSITY
from tmesim.reconstruct.real_pipeline import (
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

PANEL = ["temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin"]

# Per-cell sensitivity modes. 'uniform' is the HEADLINE (fair cross-drug) mode;
# 'determinant' is reported alongside it so the confound is visible, never hidden.
RANKING_MODE = "uniform"
SENSITIVITY_MODES = ("uniform", "determinant")

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "real_drug_panel.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "real_drug_panel.md")

# CONVERGED grid resolution. penetration_index only settles once dx resolves the
# ~45 um perivascular boundary layer (see results/numerical_rigor.md section 2);
# dx = 12 um is inside the converged regime, dx = 30 um is NOT.
DX_UM = 12.0


def compile_real_drug(name: str, alpha: float):
    """Compile ONE real drug through the full pipeline: calibrated perivascular
    consumption (alpha fit once on doxorubicin) + saturable binding-site barrier.
    The barrier is a no-op for small molecules with no saturable high-affinity
    target -- the correct, honest behaviour -- but every drug is routed through the
    identical pipeline an antibody would use."""
    dp = calibrated_drugparams(name, alpha)
    mol = get_molecule(name)
    dp = apply_barrier_to_drugparams(dp, mol)
    return dp, mol


def effective_diffusivity(diffusivity: float, ecm_density: float) -> float:
    """Effective tissue diffusivity the PDE actually experiences after ECM
    tortuosity: D_eff = D * (1 - ALPHA_TORTUOSITY * ecm) (um^2/s)."""
    return float(diffusivity) * (1.0 - ALPHA_TORTUOSITY * float(np.clip(ecm_density, 0.0, 1.0)))


def main(dx_um: float = DX_UM):
    t0 = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    # --- 1. REAL whole-niche microenvironment (Xenium), reconstructed ONCE -----
    print(f"Reconstructing the REAL Xenium GBM microenvironment (once, dx={dx_um} um) ...")
    micro = reconstruct_real(source="xenium", dx_um=float(dx_um))
    report = reconstruction_report(micro)
    mean_ecm = report["mean_ecm_density"]
    print(
        f"  {report['n_cells']:,} real cells  (tumor {report['n_tumor']:,})  "
        f"grid {report['grid']}  dx={report['dx_um']:.0f} um\n"
        f"  mean vascular {report['mean_vascular_density']:.3f}  "
        f"ecm {mean_ecm:.3f}  o2 {report['mean_oxygen']:.3f}"
    )

    # --- 2. Single calibration constant (fit ONCE on doxorubicin) -------------
    alpha = calibrate_alpha()
    print(f"\nperivascular-sink calibration alpha = {alpha:.1f} (fit once on doxorubicin)")

    # --- 3. Simulate every drug on the SAME real tissue, under BOTH modes ------
    # Sensitivity is the ONLY thing that differs between the two passes: same
    # reconstruction, same drugs, same PDE, same seed.
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)
    by_mode = {}
    for mode in SENSITIVITY_MODES:
        print(f"\n--- per-cell sensitivity mode = '{mode}' "
              f"{'(HEADLINE, fair cross-drug)' if mode == RANKING_MODE else '(comparison only)'}")
        rows = []
        for name in PANEL:
            dp, mol = compile_real_drug(name, alpha)
            info = apply_drug_sensitivity(micro, name, mode=mode)
            cfg = SimConfig(drug=dp, pde=pde, seed=0)
            print(f"  simulating {name:14s} ...", end="", flush=True)
            res = simulate_response(micro, cfg)
            d_eff = effective_diffusivity(dp.diffusivity, mean_ecm)
            rows.append(
                dict(
                    name=name,
                    mw=float(mol.get("mw")),
                    ic50=float(dp.ic50),
                    bbb=float(dp.bbb_penetrance),
                    diffusivity=float(dp.diffusivity),
                    d_eff=d_eff,
                    binding=float(dp.binding),
                    clearance=float(dp.clearance),
                    kill=float(res.kill_fraction),
                    pen=float(res.extras["penetration_index"]),
                    mean_auc=float(res.extras["mean_auc"]),
                    hypoxic=float(res.extras["hypoxic_tumor_fraction"]),
                    auc_field=np.asarray(res.fields["auc"]),
                    sens_mode=mode,
                    sens_gene=info.get("determinant_gene"),
                    sens_mean=float(info.get("mean", 1.0)),
                    upper_bound=bool(info.get("upper_bound", True)),
                )
            )
            flag = "  [UPPER BOUND: no resistance modelled]" if (
                mode == "determinant" and rows[-1]["upper_bound"]) else ""
            print(f" kill={rows[-1]['kill']:.4f}  pen={rows[-1]['pen']:.3f}  "
                  f"mean_auc={rows[-1]['mean_auc']:.3g}  "
                  f"sens={rows[-1]['sens_mean']:.3f}"
                  f"({rows[-1]['sens_gene'] or 'none'}){flag}")
        by_mode[mode] = rows

    # Restore the headline (fair) sensitivity on the shared reconstruction.
    apply_drug_sensitivity(micro, PANEL[0], mode=RANKING_MODE)

    # --- 4. Rankings (HEADLINE = the fair uniform mode) -----------------------
    rows = by_mode[RANKING_MODE]
    best_kill = max(rows, key=lambda r: r["kill"])
    best_pen = max(rows, key=lambda r: r["pen"])
    hypoxic_frac = rows[0]["hypoxic"]  # identical tissue for all drugs

    # --- 5. Figure: per-drug AUC heatmaps (shared scale) + bar chart ----------
    W, H = micro.domain.extent_um
    ext = [0, W, H, 0]
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    # Shared colour scale from the 99th percentile of tissue AUC across ALL drugs.
    all_tissue_auc = np.concatenate([r["auc_field"][tissue] for r in rows])
    vmax = float(np.percentile(all_tissue_auc, 99.0)) or 1.0

    fig = plt.figure(figsize=(19, 9))
    gs = GridSpec(2, len(PANEL), figure=fig, height_ratios=[1.15, 1.0],
                  hspace=0.28, wspace=0.12)

    im = None
    for j, r in enumerate(rows):
        ax = fig.add_subplot(gs[0, j])
        im = ax.imshow(r["auc_field"], cmap="viridis", extent=ext, aspect="equal",
                       vmin=0, vmax=vmax)
        ax.set_title(f"{r['name']}\npen={r['pen']:.3f}  kill={r['kill']:.3f}",
                     fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
    cbar = fig.colorbar(im, ax=fig.axes[:len(PANEL)], shrink=0.85, pad=0.01,
                        location="right")
    cbar.set_label("drug exposure AUC (uM*s), shared scale")

    # Bottom: grouped bar chart of kill_fraction & penetration_index per drug.
    ax_bar = fig.add_subplot(gs[1, :])
    names = [r["name"] for r in rows]
    kills = [r["kill"] for r in rows]
    pens = [r["pen"] for r in rows]
    x = np.arange(len(names))
    w = 0.38
    b1 = ax_bar.bar(x - w / 2, kills, w, label="predicted kill_fraction", color="#d62728")
    b2 = ax_bar.bar(x + w / 2, pens, w, label="penetration_index", color="#1f77b4")
    for b in list(b1) + list(b2):
        ax_bar.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.005,
                    f"{b.get_height():.3f}", ha="center", va="bottom", fontsize=8)
    ax_bar.set_xticks(x); ax_bar.set_xticklabels(names)
    ax_bar.set_ylabel("fraction / index")
    ax_bar.set_ylim(0, max(max(kills), max(pens)) * 1.25 + 1e-6)
    ax_bar.legend(loc="upper right", fontsize=9)
    pen_lo, pen_hi = min(pens), max(pens)
    ax_bar.set_title(
        f"Predicted kill & penetration on the SAME real tumor  "
        f"(best kill: {best_kill['name']}, kill={best_kill['kill']:.3f};  "
        f"tumor hypoxic fraction {hypoxic_frac*100:.0f}%)\n"
        f"penetration_index is FLAT across the panel at converged dx "
        f"({pen_lo:.3f}–{pen_hi:.3f}) — it is a shape statistic set by vascular "
        f"geometry, NOT a drug ranking; compare drugs by mean AUC / kill instead",
        fontsize=9.5,
    )

    fig.suptitle(
        "REAL drug panel on the REAL Xenium glioblastoma microenvironment  "
        "(one tumor, five drugs)  —  PREDICTIONS, not validated outcomes\n"
        f"grid-CONVERGED dx = {report['dx_um']:.0f} µm · "
        f"α = {alpha:.1f} (e-fold convention) · "
        f"per-cell sensitivity mode = '{RANKING_MODE}' (identical for every drug "
        "→ fair cross-drug ranking) · "
        "GEOMETRIC (dx-independent) tissue mask",
        fontsize=13,
    )
    fig.savefig(_OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"\nfigure written: {_OUT_PNG}")

    # --- 6. Markdown report ---------------------------------------------------
    md = build_markdown(report, alpha, rows, best_kill, best_pen, hypoxic_frac,
                        mean_ecm, by_mode)
    with open(_OUT_MD, "w") as fh:
        fh.write(md)
    print(f"report written: {_OUT_MD}")

    # --- 7. stdout summary ----------------------------------------------------
    print("\n" + "=" * 92)
    print("REAL DRUG PANEL on the REAL Xenium GBM microenvironment (predictions only)")
    print(f"HEADLINE sensitivity mode = '{RANKING_MODE}' (identical across drugs "
          "-> fair cross-drug ranking)")
    print("=" * 92)
    hdr = (f"{'drug':14s} {'MW':>7} {'IC50(uM)':>9} {'bbb':>6} {'D_eff':>8} "
           f"{'pen_idx':>8} {'mean_auc':>10} {'kill':>7}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['name']:14s} {r['mw']:7.0f} {r['ic50']:9.3g} {r['bbb']:6.3f} "
              f"{r['d_eff']:8.1f} {r['pen']:8.3f} {r['mean_auc']:10.3g} {r['kill']:7.4f}")
    print("=" * 92)
    print("RANKING UNDER BOTH SENSITIVITY MODES (kill_fraction, best first)")
    for mode in SENSITIVITY_MODES:
        ordered = sorted(by_mode[mode], key=lambda r: r["kill"], reverse=True)
        tag = "FAIR / headline" if mode == RANKING_MODE else "CONFOUNDED — do not rank"
        print(f"  mode='{mode}' [{tag}]")
        for i, r in enumerate(ordered, 1):
            ub = " <- UPPER BOUND (determinant unmeasured)" if (
                mode == "determinant" and r["upper_bound"]) else ""
            print(f"    {i}. {r['name']:14s} kill={r['kill']:.4f}  "
                  f"sens={r['sens_mean']:.3f} ({r['sens_gene'] or 'none'}){ub}")
    print("=" * 92)
    print(f"best predicted KILL       : {best_kill['name']} "
          f"(kill={best_kill['kill']:.4f}, pen={best_kill['pen']:.3f}, "
          f"IC50={best_kill['ic50']:.3g} uM, bbb={best_kill['bbb']:.3f})")
    print(f"best predicted PENETRATION: {best_pen['name']} "
          f"(pen={best_pen['pen']:.3f}, D_eff={best_pen['d_eff']:.0f} um^2/s)")
    print(f"tumor hypoxic fraction    : {hypoxic_frac*100:.1f}%")
    print(f"grid resolution           : dx={report['dx_um']:.0f} um {report['grid']}  "
          f"(CONVERGED; penetration_index needs dx<=12um -- numerical_rigor.md)")
    print(f"tissue mask               : {report['tissue_mask_mode']} "
          f"(r={report['tissue_disc_radius_um']} um), coverage "
          f"{report['tissue_coverage_fraction']*100:.1f}% -- dx-INDEPENDENT")
    print("=" * 92)
    print("HONESTY: untreated resection -> NO outcome label. Real geometry + real drug\n"
          "physics + real potency -> PREDICTIONS, not validated responses. Do not read the\n"
          "ranking as clinical efficacy.")
    print("=" * 92)
    print(f"[runtime {time.time() - t0:.1f} s]")


def _ranking_line(rows) -> str:
    ordered = sorted(rows, key=lambda r: r["kill"], reverse=True)
    return " > ".join(f"**{r['name']}** ({r['kill']:.4f})" for r in ordered)


def build_markdown(report, alpha, rows, best_kill, best_pen, hypoxic_frac,
                   mean_ecm, by_mode) -> str:
    frac = report["celltype_fractions"]
    lines = [
        "# Real drug panel on the real Xenium glioblastoma microenvironment",
        "",
        "**Experiment:** `experiments/real_drug_panel.py` "
        "(run: `./.venv/bin/python experiments/real_drug_panel.py`)",
        "",
        "**Figure:** `results/real_drug_panel.png`",
        "",
        f"> **Provenance of this run.** Perivascular-sink calibration "
        f"**α = {alpha:.1f}** (the e-folding-length convention, "
        f"`L = d_half/ln 2`; the retired α = 309.6 came from fitting the raw "
        f"half-distance as if it were an e-fold length and is superseded). "
        f"Per-cell sensitivity mode = **`{RANKING_MODE}`** for the headline "
        f"ranking — identical for every drug, so the comparison is fair; the "
        f"drug-specific `determinant` mode is reported separately in §*Sensitivity "
        f"mode* below. Tissue mask = **geometric** (union of "
        f"{report['tissue_disc_radius_um']:.0f} µm discs around the measured "
        f"cell centroids), which is **independent of the grid spacing** — this "
        f"replaces the old \"grid cells containing ≥ 1 cell\" mask and changes "
        f"the oxygen field and the hypoxic fraction substantially (see §*Tissue "
        f"mask*).",
        "",
        "The full panel "
        "**[temozolomide, lomustine, cisplatin, methotrexate, doxorubicin]** is "
        "simulated on **one and the same real tumor** — the open 10x Xenium GBM "
        f"section ({report['n_cells']:,} measured cells; tumor {report['n_tumor']:,}). "
        "The real microenvironment is reconstructed **once** and reused for every "
        "drug, so all differences come from each drug's own physics on identical "
        "tissue (real vasculature, real ECM from collagen/matrix transcripts, "
        "real-perfusion oxygen).",
        "",
        "Each drug is compiled through the same pipeline: "
        "`calibrate_alpha()` (single perivascular-sink constant, fit **once** on "
        f"doxorubicin → α={alpha:.1f}) → `calibrated_drugparams(name, α)` (real GDSC "
        "IC50 + calibrated loss) → `apply_barrier_to_drugparams` (saturable "
        "binding-site barrier).",
        "",
        "## Real microenvironment (reconstructed once)",
        "",
        f"- Source: 10x Xenium GBM — {report['n_cells']:,} real cells, "
        f"grid {report['grid']}, dx={report['dx_um']:.0f} µm, "
        f"extent {report['extent_um']} µm.",
        f"- Cell-type fractions: TUMOR {frac['TUMOR']*100:.1f}%, "
        f"ENDOTHELIAL {frac['ENDOTHELIAL']*100:.1f}%, "
        f"IMMUNE {frac['IMMUNE']*100:.1f}%, STROMAL {frac['STROMAL']*100:.1f}%.",
        f"- Mean vascular density {report['mean_vascular_density']:.3f}, "
        f"mean ECM {mean_ecm:.3f}, mean oxygen {report['mean_oxygen']:.3f}.",
        f"- Tissue mask: **geometric**, radius "
        f"{report['tissue_disc_radius_um']:.0f} µm, covering "
        f"{report['tissue_coverage_fraction']*100:.1f}% of the bounding box "
        "(dx-independent — see §*Tissue mask*).",
        f"- **Tumor hypoxic fraction ≈ {hypoxic_frac*100:.1f}%** at this "
        "resolution under the grid-independent tissue mask. See §*Tissue mask* — "
        "this number moved a great deal when the mask was fixed and the honest "
        "reading has changed.",
        "",
        "## Tissue mask — the geometry fix, and what it changed",
        "",
        "`Domain.tissue_mask` used to be *\"grid cells that contain at least one "
        "segmented cell\"*. That is a property of the **grid**, not of the "
        "tissue: as `dx` shrinks each cell can still only mark one tile, so the "
        "covered fraction collapsed — **50.6% / 32.7% / 13.9%** of the bounding "
        "box at dx = 30 / 20 / 12 µm. `real_fields._kde_on_grid` multiplies the "
        "vascular field by that mask, so the **total O₂ supply shrank as the mesh "
        "was refined**, which is a discretisation artifact, not biology.",
        "",
        "The mask is now **geometric**: a grid cell is tissue when its centre lies "
        f"within **{report['tissue_disc_radius_um']:.0f} µm** of at least one "
        "measured cell centroid (the union of fixed physical discs). The radius "
        "is a declared constant justified against this section's own geometry — "
        "the nearest-neighbour spacing between the 40,887 measured centroids is "
        "median **16.7 µm** (p10 10.4, p90 24.2, p99 34.2 µm), so 30 µm is ~1.8× "
        "the median and just above the p90 spacing: the union of discs bridges "
        ">90% of real inter-cell gaps (giving a *connected* domain, which the "
        "diffusion operator needs) while extending only ~1 cell diameter beyond "
        "the outermost measured cell. It is not fitted to anything.",
        "",
        "Measured effect (same section, same physics, only the mask definition "
        "changed):",
        "",
        "| dx (µm) | coverage — legacy | coverage — geometric | mean O₂ @cells — "
        "legacy | mean O₂ @cells — geometric | hypoxic tumor frac — legacy | "
        "hypoxic tumor frac — geometric |",
        "|---:|---:|---:|---:|---:|---:|---:|",
        "| 30 | 0.506 | 0.672 | 0.522 | 0.570 | 0.000 | 0.000 |",
        "| 20 | 0.327 | 0.674 | 0.426 | 0.565 | 0.000 | 0.000 |",
        "| 12 | 0.139 | 0.677 | 0.251 | 0.558 | 0.134 | 0.000 |",
        "| **spread (max/min)** | **3.64×** | **1.007×** | **2.08×** | "
        "**1.02×** | — | — |",
        "",
        "So the geometric mask does what it was built to do: tissue coverage and "
        "mean oxygen are now invariant to the grid (1.007× and 1.02× across a "
        "2.5× range of dx), where before they varied by 3.6× and 2.1×.",
        "",
        "**Honest consequence — the hypoxia story changes again, and not in the "
        "model's favour.** With the O₂ supply no longer artificially starved by a "
        "shrinking mask, mean oxygen at the measured cell positions rises to "
        "≈0.56–0.57 and **no tumor cell falls below the 0.2 hypoxia threshold at "
        "any dx: the predicted hypoxic fraction is 0.0%.** For a glioblastoma "
        "section — a tumor type defined partly by pseudopalisading necrosis — "
        "that is not credible as biology. It is best read as evidence about the "
        "*oxygen model*, not about this tumor: `OxygenConfig.consumption` has "
        "never been calibrated against a measured oxygen tension, and the "
        "previously-reported hypoxia (94% at dx = 30, 13% at dx = 12) was being "
        "manufactured by the grid-dependent supply term rather than by the "
        "tissue. **Both the old numbers and this 0% are un-calibrated model "
        "output; the correct statement is that this simulator cannot currently "
        "predict this section's hypoxic fraction, and no conclusion in this "
        "report should rest on it.** What the fix does establish is that the "
        "answer no longer depends on the mesh.",
        "",
        "## Sensitivity mode — why the headline ranking uses `uniform`",
        "",
        "Per-cell drug sensitivity can be set two ways "
        "(`reconstruct.real_pipeline.apply_drug_sensitivity`), and the choice "
        "changes the ranking:",
        "",
        "- **`mode='uniform'` (default, and the headline table above).** Every "
        "drug gets sensitivity 1.0 at every cell. Because the factor is "
        "*identical across drugs* it cancels out of the comparison, so the "
        "ranking isolates **delivery × potency** with no sensitivity confound. "
        "Absolute kills are an equal no-intrinsic-resistance **upper bound for "
        "every drug**.",
        "- **`mode='determinant'`.** Each drug uses its own measured determinant "
        "where the 319-gene Xenium panel carries one — **TOP2A** for doxorubicin "
        "(its literal molecular target) and an **MKI67** proliferation proxy for "
        "cisplatin / lomustine / methotrexate — and uniform 1.0 where it does "
        "not. On this section both determinants land at mean sensitivity ≈0.585. "
        "**Temozolomide has no determinant on this panel** (its is MGMT, which is "
        "absent and is a methylation state that transcript counts cannot report), "
        "so TMZ *alone* keeps 1.0 while every competitor is scaled down ~1.7×.",
        "",
        "**That asymmetry is an artifact of measurement coverage, not of drug "
        "biology.** Any kill temozolomide gains in `determinant` mode is bought "
        "purely by its resistance gene being unmeasured. It is reported here, "
        "flagged **UPPER BOUND**, and it is NOT the ranking:",
        "",
        "| drug | sensitivity (uniform) | sensitivity (determinant) | determinant "
        "gene | kill (uniform) | kill (determinant) | comparable across drugs? |",
        "|---|---:|---:|:--|---:|---:|:--|",
    ]
    _u = {r["name"]: r for r in by_mode["uniform"]}
    _d = {r["name"]: r for r in by_mode["determinant"]}
    for name in PANEL:
        u, d = _u[name], _d[name]
        flag = ("❌ **UPPER BOUND** — no resistance modelled"
                if d["upper_bound"] else "scaled by its own determinant")
        lines.append(
            f"| {name} | {u['sens_mean']:.3f} | {d['sens_mean']:.3f} | "
            f"{d['sens_gene'] or '— (none on panel)'} | {u['kill']:.4f} | "
            f"{d['kill']:.4f} | {flag} |")
    lines += [
        "",
        f"**Ranking under `mode='uniform'` (FAIR — this is the ranking):** "
        f"{_ranking_line(by_mode['uniform'])}",
        "",
        f"**Ranking under `mode='determinant'` (CONFOUNDED — reported for "
        f"comparison only, do not cite as a ranking):** "
        f"{_ranking_line(by_mode['determinant'])}",
        "",
        "Read the difference plainly: in `determinant` mode every drug except "
        "temozolomide is multiplied by ≈0.585, so their kills drop by roughly "
        "that factor while TMZ's does not move at all. Any place TMZ climbs is "
        "**entirely** attributable to that, and would reverse the moment MGMT "
        "were measurable. `determinant` mode is legitimate for a *within-drug* "
        "spatial question (\"where on this section does doxorubicin work?\"); it "
        "is not legitimate for ranking drugs against each other on a panel that "
        "measures some of their determinants and not others.",
        "",
        "## Predicted panel on this real tumor",
        "",
        f"*(sensitivity mode = `{RANKING_MODE}` — identical for every drug.)*",
        "",
        "| drug | MW (Da) | real IC50 (µM) | bbb_penetrance | D_eff (µm²/s) | "
        "penetration_index | mean AUC (µM·s) | predicted kill |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in rows:
        lines.append(
            f"| {r['name']} | {r['mw']:.0f} | {r['ic50']:.3g} | {r['bbb']:.3f} | "
            f"{r['d_eff']:.1f} | {r['pen']:.3f} | {r['mean_auc']:.3g} | "
            f"{r['kill']:.4f} |"
        )
    lines += [
        "",
        "IC50 values are the real GDSC/literature potencies that flow in through the "
        "molecular layer; MW is from the verified molecular structure; `bbb_penetrance` "
        "and `D_eff` are structure-derived (`D_eff = D · (1 − "
        f"{ALPHA_TORTUOSITY}·ECM)` at the tumor mean ECM = {mean_ecm:.3f}).",
        "",
        "## Which drug wins on THIS real tumor, and why",
        "",
        f"- **Penetration does NOT separate these drugs at converged resolution.** "
        f"`penetration_index` spans only {min(r['pen'] for r in rows):.3f}–"
        f"{max(r['pen'] for r in rows):.3f} across the whole panel — a spread of "
        f"{max(r['pen'] for r in rows) - min(r['pen'] for r in rows):.3f}, which is "
        "not a meaningful ranking. This is a **correction**: an earlier, "
        "unconverged dx = 30 µm run of the same panel appeared to span 0.266–0.396 "
        "and was reported as 'temozolomide penetrates best, doxorubicin barely "
        "reaches the tumor bulk'. **That ranking was a grid artifact and is "
        "retracted.** The "
        "reason is definitional: `penetration_index = mean(AUC) / p99(AUC)` over "
        "tissue is a *scale-invariant shape* statistic — it measures how uniform the "
        "exposure field is, not how much drug arrives. Once the perivascular boundary "
        "layer is actually resolved, that shape is set by the **vascular geometry**, "
        "which is identical for all five drugs, and these drugs' D_eff values differ "
        "by only ~1.5× (175–263 µm²/s) — far too little to reshape the field.",
        f"- **Delivery differences are real, but they live in the AUC, not the index.** "
        f"Mean exposure spans {min(r['mean_auc'] for r in rows):.3g}–"
        f"{max(r['mean_auc'] for r in rows):.3g} µM·s "
        f"(~{max(r['mean_auc'] for r in rows) / max(min(r['mean_auc'] for r in rows), 1e-30):.0f}×), "
        "driven by BBB penetrance and the calibrated perivascular sink. Doxorubicin "
        "really is delivery-starved on this section — its mean AUC is ~120× below "
        "temozolomide's — but that shows up as a *lower absolute exposure*, not as a "
        "flatter spatial profile. **Quote mean AUC, not penetration_index, when "
        "comparing drugs.**",
        f"- **Best predicted kill (mode = `{RANKING_MODE}`):** "
        f"**{best_kill['name']}** "
        f"(kill_fraction {best_kill['kill']:.4f}, penetration_index {best_kill['pen']:.3f}, "
        f"IC50 {best_kill['ic50']:.3g} µM, mean AUC {best_kill['mean_auc']:.3g} µM·s, "
        f"bbb {best_kill['bbb']:.3f}). Kill is **delivery × potency**, discounted by "
        "whatever oxygen-driven resistance the tissue carries: a drug wins only when "
        "it both reaches the cells AND is potent enough at the exposure it achieves. "
        f"Full ranking: {_ranking_line(rows)}. Because every drug carries the same "
        "sensitivity in this mode, none of this ordering is a sensitivity artifact — "
        "it is delivery and potency only.",
        "",
        "## Honesty and caveats",
        "",
        "- **No outcome label.** This open Xenium section is an **untreated "
        "resection**. There is no paired measured drug response for any of these five "
        "drugs. The kill/penetration numbers are **predictions of the validated "
        "physics on real tissue — not validated outcomes.** Do not read the ranking as "
        "clinical efficacy.",
        "- **What is genuinely real here:** the geometry (real cells, real "
        "vasculature, real matrix, real oxygen), the transport physics (structure → "
        "mesh diffusivity + a doxorubicin-calibrated perivascular sink, revalidated "
        "out-of-sample on mitoxantrone), and the potency (real GDSC/literature IC50).",
        f"- **Hypoxia is now predicted at {hypoxic_frac*100:.1f}%, and that is a "
        "statement about the model, not about this tumor.** Every earlier hypoxia "
        "number from this pipeline (≈94% at dx = 30 µm, ≈13% at dx = 12 µm) was "
        "produced under the grid-dependent occupancy tissue mask, which starved "
        "the O₂ supply more and more as the mesh was refined; those numbers are "
        "**retracted**. Under the fixed geometric mask the hypoxic fraction is "
        "0.0% at every dx — which is *also* not credible for a GBM section. "
        "`OxygenConfig.consumption` has never been calibrated against a measured "
        "oxygen tension, so the absolute level is un-anchored in both directions. "
        "**No claim in this report depends on the hypoxic fraction, and none "
        "should.**",
        "- **Absolute kill fractions are an upper bound.** In the headline "
        "`uniform` sensitivity mode no intrinsic per-cell resistance is modelled "
        "for any drug. That is what makes the *ranking* fair; it also means every "
        "*absolute* kill here is the best case. See §*Sensitivity mode*.",
        "- **One anchor sets the scale.** The single α is fit to one doxorubicin "
        "penetration number; the absolute penetration of every drug rides on that "
        "measurement, and the binding proxy is an ordering, not a first-principles "
        "rate.",
        "",
        "## Resolution & convergence",
        "",
        f"This panel was run at **dx = {report['dx_um']:.0f} µm** on grid "
        f"{report['grid']}. That choice is not cosmetic — see "
        "[`results/numerical_rigor.md`](numerical_rigor.md) §2 (grid refinement on "
        "this same real section):",
        "",
        "- **`kill_fraction` is grid-converged; `penetration_index` requires "
        "dx ≲ 12 µm.** The calibrated perivascular boundary layer is ~45 µm "
        "thick, so any coarser grid cannot resolve the steep near-vessel gradient "
        "and reports an artificially flat — i.e. too high — index. That is why "
        "this panel runs at dx = 12 µm. dx = 12 µm is also the finest grid that "
        "satisfies the upwind-advection CFL condition at the published dt = 10 s.",
        "- **An earlier dx = 30 µm run of this panel OVER-ESTIMATED "
        "`penetration_index`** by 1.4–2.0×, drug-dependently, which manufactured "
        "a spurious penetration ranking (see above). Those figures are superseded; "
        "do not cite the dx = 30 µm penetration numbers.",
        "- **⚠️ The per-drug convergence ladder previously printed here is "
        "RETRACTED, not merely updated.** It was measured under the legacy "
        "grid-dependent occupancy tissue mask, so its dx sweep confounded genuine "
        "mesh refinement with a shrinking simulated tissue domain — exactly the "
        "defect §*Tissue mask* fixes. Re-deriving a 5-drug × 6-grid ladder under "
        "the geometric mask is a separate compute job that was not run for this "
        "report. **For the current, geometric-mask grid study see "
        "[`results/numerical_rigor.md`](numerical_rigor.md) §2**, which is "
        "regenerated alongside this file on the same geometry (doxorubicin "
        "reference drug, advection off).",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Real drug panel on the real Xenium GBM section")
    ap.add_argument("--dx", type=float, default=DX_UM,
                    help="grid spacing in um (default %(default)s, grid-converged)")
    main(dx_um=ap.parse_args().dx)
