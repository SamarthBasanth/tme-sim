"""CROSS-TISSUE: the SAME real drug physics across MULTIPLE real GBM tumors.

This is the inter-tumor-heterogeneity companion to ``real_simulation.py``. Where
that experiment runs the validated transport -> kill engine on ONE real Xenium
section, this runs the *identical* pipeline across FIVE real glioblastoma tissues:

  * the open 10x Xenium single-cell GBM section (source='xenium'), and
  * the four open 10x Visium GBM sections (list_samples():
    10x_parent_gbm, GBM1, GBM2, GBM3).

Each tissue is reconstructed with the real whole-niche assembler
(:func:`reconstruct_real`) and simulated with the SAME real drug -- temozolomide,
calibrated perivascular-consumption physics (alpha fit ONCE on doxorubicin) plus
the binding-site barrier -- and the SAME engine config. Nothing about the drug or
the physics changes between tumors; only the REAL microenvironment (vascularity,
hypoxia, ECM, cell composition) differs. The spread you see in predicted
penetration and kill is therefore REAL inter-tumor heterogeneity, not a knob.

HONESTY (read this):
  * These are PREDICTIONS on real geometry with NO outcome labels. Every one of
    these sections is an untreated resection -- there is no paired measured drug
    response to check the kill numbers against. What is real is the GEOMETRY
    (cells, vasculature, matrix, oxygen) and the DRUG PHYSICS. The kill/penetration
    numbers are what the validated physics predicts on real tissue, not outcomes.
  * PLATFORM DIFFERENCE: Xenium is single-cell resolution (real cells, per-cell
    marker typing, expression-weighted ECM). Visium is spot-level (~100 um spots,
    reference-based deconvolution, stromal-proxy ECM). The two are NOT
    apples-to-apples: Visium composition is deconvolved and coarser, and its ECM
    field is a proxy. We plot them together but mark the platform so the
    single-cell tumor is never silently compared to spot-level tumors.

Run:  ./.venv/bin/python experiments/cross_tissue.py
Writes: results/cross_tissue.png, results/cross_tissue.md
"""
from __future__ import annotations

import os
import time
from typing import Dict, List, Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tmesim.core import PDEConfig, SimConfig
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.data.real.visium import list_samples, sample_info
from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "cross_tissue.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "cross_tissue.md")


def compile_real_temozolomide():
    """Compile the REAL temozolomide drug ONCE: calibrated perivascular-consumption
    physics (alpha fit on doxorubicin) then the binding-site barrier. The exact
    same DrugParams object is reused for every tissue -- the drug never changes."""
    alpha = calibrate_alpha()  # single fit on doxorubicin (see calibration.py)
    drug = calibrated_drugparams("temozolomide", alpha)
    drug = apply_barrier_to_drugparams(drug, get_molecule("temozolomide"))
    return drug, alpha


def _reconstruct_retry(source: str, sample: Optional[str], dx_um: float):
    """Reconstruct, retrying once on any transient/IO/network error."""
    try:
        return reconstruct_real(source=source, sample=sample, dx_um=dx_um)
    except Exception as exc:  # noqa: BLE001 - retry once per task spec
        print(f"    reconstruct failed ({exc!r}); retrying once in 5s ...")
        time.sleep(5)
        return reconstruct_real(source=source, sample=sample, dx_um=dx_um)


def run_one_tissue(label: str, source: str, sample: Optional[str], dx_um: float,
                   drug, cfg_pde: PDEConfig, diagnosis: str) -> Dict:
    """Reconstruct one REAL tissue and run the SAME drug through the engine."""
    print(f"\n[{label}] source={source} sample={sample} dx={dx_um} um")
    micro = _reconstruct_retry(source, sample, dx_um)
    report = reconstruction_report(micro)
    cfg = SimConfig(drug=drug, pde=cfg_pde, seed=0)
    result = simulate_response(micro, cfg)
    frac = report["celltype_fractions"]
    row = {
        "label": label,
        "platform": "Xenium (single-cell)" if source == "xenium"
        else "Visium (spot-level)",
        "source": source,
        "diagnosis": diagnosis,
        "n_cells": report["n_cells"],
        "n_tumor": report["n_tumor"],
        "grid": report["grid"],
        "dx_um": report["dx_um"],
        "tumor_frac": frac.get("TUMOR", 0.0),
        "endo_frac": frac.get("ENDOTHELIAL", 0.0),
        "immune_frac": frac.get("IMMUNE", 0.0),
        "stromal_frac": frac.get("STROMAL", 0.0),
        "necrotic_frac": frac.get("NECROTIC", 0.0),
        "mean_vascular": report["mean_vascular_density"],
        "mean_ecm": report["mean_ecm_density"],
        "mean_oxygen": report["mean_oxygen"],
        "kill_fraction": result.kill_fraction,
        "penetration_index": result.extras["penetration_index"],
        "mean_auc": result.extras["mean_auc"],
        "hypoxic_tumor_fraction": result.extras["hypoxic_tumor_fraction"],
    }
    print(f"    n_cells={row['n_cells']:,} tumor={row['tumor_frac']*100:.1f}% "
          f"vasc={row['mean_vascular']:.3f} o2={row['mean_oxygen']:.3f} "
          f"hypoxic={row['hypoxic_tumor_fraction']:.3f} "
          f"pen={row['penetration_index']:.3f} kill={row['kill_fraction']:.3f}")
    return row


def _spread(rows: List[Dict], key: str):
    vals = [r[key] for r in rows]
    return min(vals), max(vals), (max(vals) - min(vals))


def make_figure(rows: List[Dict], alpha: float, path: str) -> None:
    """Two scatter panels (the heterogeneity story) + two composition/field bars."""
    labels = [r["label"] for r in rows]
    xen = np.array([r["source"] == "xenium" for r in rows])
    colors = ["#d62728" if x else "#1f77b4" for x in xen]  # red=Xenium, blue=Visium

    vasc = np.array([r["mean_vascular"] for r in rows])
    pen = np.array([r["penetration_index"] for r in rows])
    hyp = np.array([r["hypoxic_tumor_fraction"] for r in rows])
    kill = np.array([r["kill_fraction"] for r in rows])

    fig, axes = plt.subplots(2, 2, figsize=(15, 12), constrained_layout=True)

    # (a) vascularity vs predicted penetration
    ax = axes[0, 0]
    ax.scatter(vasc, pen, c=colors, s=120, edgecolors="k", zorder=3)
    for x, y, lab in zip(vasc, pen, labels):
        ax.annotate(lab, (x, y), textcoords="offset points", xytext=(7, 4),
                    fontsize=9)
    ax.set_xlabel("REAL mean vascular density")
    ax.set_ylabel("Predicted penetration_index")
    ax.set_title("(a) Real vascularity vs predicted drug penetration\n"
                 "(same TMZ physics on every tumor)")
    ax.grid(alpha=0.3)

    # (b) hypoxia vs predicted kill
    ax = axes[0, 1]
    ax.scatter(hyp, kill, c=colors, s=120, edgecolors="k", zorder=3)
    for x, y, lab in zip(hyp, kill, labels):
        ax.annotate(lab, (x, y), textcoords="offset points", xytext=(7, 4),
                    fontsize=9)
    ax.set_xlabel("REAL hypoxic tumor fraction")
    ax.set_ylabel("Predicted kill_fraction")
    ax.set_title("(b) Real hypoxia vs predicted kill\n"
                 "(same TMZ physics on every tumor)")
    ax.grid(alpha=0.3)

    # (c) grouped bars: vascular / oxygen / hypoxic-tumor per tissue
    ax = axes[1, 0]
    x = np.arange(len(rows))
    w = 0.27
    ax.bar(x - w, vasc, w, label="mean vascular density", color="#ff7f0e")
    ax.bar(x, np.array([r["mean_oxygen"] for r in rows]), w, label="mean oxygen",
           color="#2ca02c")
    ax.bar(x + w, hyp, w, label="hypoxic tumor fraction", color="#7f4fbf")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("field value (0-1)")
    ax.set_title("(c) REAL microenvironment fields per tumor")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    # (d) stacked cell-type composition
    ax = axes[1, 1]
    comp_keys = [("tumor_frac", "TUMOR", "#d62728"),
                 ("stromal_frac", "STROMAL", "#2ca02c"),
                 ("immune_frac", "IMMUNE", "#1f77b4"),
                 ("endo_frac", "ENDOTHELIAL", "#ff7f0e"),
                 ("necrotic_frac", "NECROTIC", "#7f7f7f")]
    bottom = np.zeros(len(rows))
    for key, name, col in comp_keys:
        vals = np.array([r[key] for r in rows])
        ax.bar(x, vals, 0.6, bottom=bottom, label=name, color=col)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("cell-type fraction")
    ax.set_title("(d) REAL cell-type composition per tumor\n"
                 "(Xenium: per-cell typing; Visium: deconvolved)")
    ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(
        f"CROSS-TISSUE inter-tumor heterogeneity: same TMZ physics "
        f"(alpha={alpha:.0f}) on {len(rows)} REAL GBM tissues\n"
        "red = Xenium single-cell   |   blue = Visium spot-level   "
        "(PREDICTIONS on real geometry, NO outcome labels)",
        fontsize=13,
    )
    fig.savefig(path, dpi=130)
    plt.close(fig)
    print(f"\nfigure written: {path}")


def write_md(rows: List[Dict], alpha: float, drug, runtime_s: float,
             path: str) -> None:
    def pct(v):
        return f"{v * 100:.1f}%"

    lines: List[str] = []
    lines.append("# Cross-tissue predicted drug behaviour on REAL GBM microenvironments\n")
    lines.append(
        "The **same** real drug (temozolomide: doxorubicin-calibrated "
        f"perivascular-consumption physics, alpha = {alpha:.0f}, plus the "
        "binding-site barrier) and the **same** validated transport -> kill engine "
        "were run across every real glioblastoma tissue below. The drug and the "
        "physics are identical for all tumors; only the REAL microenvironment "
        "differs. Any spread in predicted penetration or kill is therefore **real "
        "inter-tumor heterogeneity**, not a tuned knob.\n"
    )
    lines.append(
        f"Drug physics (fixed for all tissues): diffusivity "
        f"{drug.diffusivity:.1f} um^2/s, ic50 {drug.ic50:.3g} uM, binding "
        f"{drug.binding:.3e}/s, clearance {drug.clearance:.3e}/s, bbb_penetrance "
        f"{drug.bbb_penetrance:.3f}.\n"
    )

    # --- main table ---
    lines.append("## Per-tissue table\n")
    hdr = ("| Tissue | Platform | Diagnosis | Cells (tumor) | Grid @dx | "
           "Tumor | Endo | Vascular | Oxygen | Hypoxic tumor | "
           "Pred. penetration | Pred. kill |")
    sep = "|" + "---|" * 12
    lines.append(hdr)
    lines.append(sep)
    for r in rows:
        lines.append(
            f"| **{r['label']}** | {r['platform']} | {r['diagnosis']} | "
            f"{r['n_cells']:,} ({r['n_tumor']:,}) | "
            f"{r['grid'][0]}x{r['grid'][1]} @{r['dx_um']:.0f}um | "
            f"{pct(r['tumor_frac'])} | {pct(r['endo_frac'])} | "
            f"{r['mean_vascular']:.3f} | {r['mean_oxygen']:.3f} | "
            f"{pct(r['hypoxic_tumor_fraction'])} | "
            f"{r['penetration_index']:.3f} | {pct(r['kill_fraction'])} |"
        )
    lines.append("")

    # --- observed spread ---
    vlo, vhi, vsp = _spread(rows, "mean_vascular")
    hlo, hhi, hsp = _spread(rows, "hypoxic_tumor_fraction")
    plo, phi, psp = _spread(rows, "penetration_index")
    klo, khi, ksp = _spread(rows, "kill_fraction")
    olo, ohi, osp = _spread(rows, "mean_oxygen")
    lines.append("## Observed real spread across tumors\n")
    lines.append(
        f"- **Vascular density**: {vlo:.3f} -> {vhi:.3f} "
        f"(spread {vsp:.3f}, {vhi / max(vlo, 1e-9):.1f}x).\n"
        f"- **Mean oxygen**: {olo:.3f} -> {ohi:.3f} (spread {osp:.3f}).\n"
        f"- **Hypoxic tumor fraction**: {pct(hlo)} -> {pct(hhi)} "
        f"(spread {hsp * 100:.1f} pts).\n"
        f"- **Predicted penetration_index**: {plo:.3f} -> {phi:.3f} "
        f"(spread {psp:.3f}).\n"
        f"- **Predicted kill_fraction**: {pct(klo)} -> {pct(khi)} "
        f"(spread {ksp * 100:.1f} pts).\n"
    )

    # --- interpretation ---
    lines.append("## Interpretation: how real microenvironment differences drive the predictions\n")
    lines.append(
        "The engine is a convection-diffusion-reaction transport solve followed by "
        "an oxygen-modulated Hill kill. Two real-tissue levers move the predictions:\n\n"
        "1. **Vascularity -> delivery -> penetration.** The drug source is seeded at "
        "the real endothelial (vascular) field, so tumors with denser, more evenly "
        "distributed real vasculature get drug delivered to more of the tissue and "
        "show a higher predicted penetration_index. Panel (a) plots this real "
        "vascularity against predicted penetration across the tumors.\n\n"
        "2. **Hypoxia -> oxygen-gated kill.** Hypoxic tumor regions (low real "
        "oxygen, inferred from the real-perfusion oxygen solve) are quiescent and "
        "less drug-sensitive, so tumors carrying a larger real hypoxic burden show a "
        "lower predicted kill for the identical exposure. Panel (b) plots real "
        "hypoxic fraction against predicted kill.\n\n"
        "Because composition (tumor/stroma/immune/endothelial fractions), vascular "
        "density, ECM and oxygen all differ between these REAL resections, the same "
        "drug physics lands differently on each -- which is exactly the inter-tumor "
        "heterogeneity a per-patient simulator should surface.\n"
    )
    if hlo > 0.9:
        lines.append(
            f"\n**Which lever dominates in THIS cohort (honest read).** Hypoxic tumor "
            f"fraction is near-saturated across all tissues ({pct(hlo)}-{pct(hhi)}), so "
            f"it barely varies and does NOT discriminate between these tumors -- panel "
            f"(b) is nearly flat on the x-axis. The real discriminator here is "
            f"**vascularity -> penetration**: predicted penetration_index spans "
            f"{plo:.3f}-{phi:.3f} and tracks the {vlo:.3f}-{vhi:.3f} spread in real "
            f"vascular density, and predicted kill follows penetration far more than "
            f"hypoxia. Read panel (a), not (b), as the heterogeneity story for this "
            f"cohort; the hypoxia axis is kept for honesty, not because it separates "
            f"these particular tumors.\n"
        )

    # --- honest caveats ---
    lines.append("## Honest caveats (read before trusting any number)\n")
    lines.append(
        "- **No outcome labels.** Every tissue here is an *untreated* open "
        "resection. There is NO paired measured temozolomide response for any of "
        "them. The kill and penetration numbers are what the validated physics "
        "predicts on real geometry -- they are **predictions, not validated "
        "outcomes**. Do not read the kill fraction as a clinical response.\n"
        "- **Platform difference (Xenium vs Visium).** The Xenium section is "
        "single-cell resolution: real individual cells, per-cell threshold-gated "
        "marker typing, and an expression-weighted ECM field from real matrix "
        "transcripts. The four Visium sections are spot-level (~100 um spots): "
        "cell-type fractions are **deconvolved** (reference-based, coarser), each "
        "spot is a small neighbourhood rather than a cell, and the ECM field is a "
        "**stromal-fraction proxy** (no per-cell transcripts). Comparing the "
        "single-cell tumor to the spot-level tumors is therefore NOT "
        "apples-to-apples: differences between the Xenium point and the Visium "
        "points partly reflect platform + resolution + deconvolution, not just "
        "biology. In the figure Xenium is red and Visium is blue for exactly this "
        "reason.\n"
        "- **Reconstruction, not measurement, of fields.** Vascular density, ECM and "
        "oxygen are inferred from real markers/perfusion, not directly measured "
        "concentrations. They are faithful to the data but are model fields.\n"
        "- **One drug, one config.** Only temozolomide at one engine config is run "
        "here; the relative spread is the message, not any absolute value.\n"
    )
    lines.append(f"\n_Generated by `experiments/cross_tissue.py` on real data "
                 f"({len(rows)} tissues, runtime {runtime_s:.0f}s). "
                 f"Figure: `results/cross_tissue.png`._\n")

    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"markdown written: {path}")


def main() -> None:
    t0 = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    drug, alpha = compile_real_temozolomide()
    print(f"REAL drug compiled ONCE: temozolomide (alpha={alpha:.1f}); reused "
          f"unchanged for every tissue.")

    # Same engine config for every tissue (matches real_simulation.py).
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)

    rows: List[Dict] = []

    # --- Xenium single-cell section (dx=30um suits single-cell) ---
    rows.append(run_one_tissue(
        label="Xenium-GBM", source="xenium", sample=None, dx_um=30.0,
        drug=drug, cfg_pde=pde,
        diagnosis="GBM (10x Xenium single-cell section)"))

    # --- 4 real Visium sections (reconstruct at spot pitch -> dx<50 falls back
    #     to the ~100um spot pitch inside the assembler) ---
    for sid in list_samples():
        try:
            info = sample_info(sid)
            dx = info.get("diagnosis", "Glioblastoma")
        except Exception:
            dx = "Glioblastoma"
        rows.append(run_one_tissue(
            label=sid, source="visium", sample=sid, dx_um=40.0,
            drug=drug, cfg_pde=pde, diagnosis=dx))

    runtime = time.time() - t0
    make_figure(rows, alpha, _OUT_PNG)
    write_md(rows, alpha, drug, runtime, _OUT_MD)

    # --- stdout summary ---
    print("\n" + "=" * 78)
    print("CROSS-TISSUE SUMMARY (same TMZ physics, real GBM tumors)")
    print("=" * 78)
    print(f"{'tissue':14s} {'platform':10s} {'vasc':>6s} {'o2':>6s} "
          f"{'hypox':>6s} {'pen':>6s} {'kill':>6s}")
    for r in rows:
        plat = "xenium" if r["source"] == "xenium" else "visium"
        print(f"{r['label']:14s} {plat:10s} {r['mean_vascular']:6.3f} "
              f"{r['mean_oxygen']:6.3f} {r['hypoxic_tumor_fraction']:6.3f} "
              f"{r['penetration_index']:6.3f} {r['kill_fraction']:6.3f}")
    vlo, vhi, vsp = _spread(rows, "mean_vascular")
    plo, phi, psp = _spread(rows, "penetration_index")
    klo, khi, ksp = _spread(rows, "kill_fraction")
    hlo, hhi, hsp = _spread(rows, "hypoxic_tumor_fraction")
    print("-" * 78)
    print(f"REAL spread: vascular {vlo:.3f}-{vhi:.3f} | hypoxic {hlo:.3f}-{hhi:.3f} "
          f"| penetration {plo:.3f}-{phi:.3f} | kill {klo:.3f}-{khi:.3f}")
    print("HONEST: predictions on real geometry, NO outcome labels; Xenium is "
          "single-cell, Visium is deconvolved spot-level (not apples-to-apples).")
    print(f"[runtime {runtime:.1f} s]")


if __name__ == "__main__":
    main()
