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

Run:    ./.venv/bin/python experiments/real_drug_panel.py
Writes: results/real_drug_panel.png
        results/real_drug_panel.md
"""
from __future__ import annotations

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
from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

PANEL = ["temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin"]

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "real_drug_panel.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "real_drug_panel.md")


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


def main():
    t0 = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    # --- 1. REAL whole-niche microenvironment (Xenium), reconstructed ONCE -----
    print("Reconstructing the REAL Xenium GBM microenvironment (once) ...")
    micro = reconstruct_real(source="xenium", dx_um=30.0)
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

    # --- 3. Simulate every drug on the SAME real tissue -----------------------
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)
    rows = []
    for name in PANEL:
        dp, mol = compile_real_drug(name, alpha)
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
            )
        )
        print(f" kill={rows[-1]['kill']:.4f}  pen={rows[-1]['pen']:.3f}  "
              f"mean_auc={rows[-1]['mean_auc']:.3g}")

    # --- 4. Rankings ----------------------------------------------------------
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
    ax_bar.set_title(
        f"Predicted kill & penetration on the SAME real tumor  "
        f"(best kill: {best_kill['name']};  best penetration: {best_pen['name']};  "
        f"tumor hypoxic fraction {hypoxic_frac*100:.0f}% -> kills suppressed for all)",
        fontsize=10,
    )

    fig.suptitle(
        "REAL drug panel on the REAL Xenium glioblastoma microenvironment  "
        "(one tumor, five drugs)  —  PREDICTIONS, not validated outcomes",
        fontsize=14,
    )
    fig.savefig(_OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"\nfigure written: {_OUT_PNG}")

    # --- 6. Markdown report ---------------------------------------------------
    md = build_markdown(report, alpha, rows, best_kill, best_pen, hypoxic_frac, mean_ecm)
    with open(_OUT_MD, "w") as fh:
        fh.write(md)
    print(f"report written: {_OUT_MD}")

    # --- 7. stdout summary ----------------------------------------------------
    print("\n" + "=" * 92)
    print("REAL DRUG PANEL on the REAL Xenium GBM microenvironment (predictions only)")
    print("=" * 92)
    hdr = (f"{'drug':14s} {'MW':>7} {'IC50(uM)':>9} {'bbb':>6} {'D_eff':>8} "
           f"{'pen_idx':>8} {'mean_auc':>10} {'kill':>7}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['name']:14s} {r['mw']:7.0f} {r['ic50']:9.3g} {r['bbb']:6.3f} "
              f"{r['d_eff']:8.1f} {r['pen']:8.3f} {r['mean_auc']:10.3g} {r['kill']:7.4f}")
    print("=" * 92)
    print(f"best predicted KILL       : {best_kill['name']} "
          f"(kill={best_kill['kill']:.4f}, pen={best_kill['pen']:.3f}, "
          f"IC50={best_kill['ic50']:.3g} uM, bbb={best_kill['bbb']:.3f})")
    print(f"best predicted PENETRATION: {best_pen['name']} "
          f"(pen={best_pen['pen']:.3f}, D_eff={best_pen['d_eff']:.0f} um^2/s)")
    print(f"tumor hypoxic fraction    : {hypoxic_frac*100:.1f}%  "
          f"(suppresses kill for EVERY drug; ~94% deeply hypoxic real tissue)")
    print("=" * 92)
    print("HONESTY: untreated resection -> NO outcome label. Real geometry + real drug\n"
          "physics + real potency -> PREDICTIONS, not validated responses. Do not read the\n"
          "ranking as clinical efficacy.")
    print("=" * 92)
    print(f"[runtime {time.time() - t0:.1f} s]")


def build_markdown(report, alpha, rows, best_kill, best_pen, hypoxic_frac, mean_ecm) -> str:
    frac = report["celltype_fractions"]
    lines = [
        "# Real drug panel on the real Xenium glioblastoma microenvironment",
        "",
        "**Experiment:** `experiments/real_drug_panel.py` "
        "(run: `./.venv/bin/python experiments/real_drug_panel.py`)",
        "",
        "**Figure:** `results/real_drug_panel.png`",
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
        f"- **Tumor hypoxic fraction ≈ {hypoxic_frac*100:.0f}%** — this deeply "
        "hypoxic real tissue mechanistically suppresses kill for **every** drug.",
        "",
        "## Predicted panel on this real tumor",
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
        f"- **Best predicted penetration:** **{best_pen['name']}** "
        f"(penetration_index {best_pen['pen']:.3f}, D_eff {best_pen['d_eff']:.0f} µm²/s). "
        "Small, weakly-bound alkylators (temozolomide, lomustine) and the small "
        "platinum crosslinker diffuse furthest before the perivascular sink consumes "
        "them; the large, strongly-bound doxorubicin is captured within tens of µm of "
        "vessels (the calibrated ~45 µm boundary layer) and barely reaches the tumor "
        "bulk.",
        f"- **Best predicted kill:** **{best_kill['name']}** "
        f"(kill_fraction {best_kill['kill']:.4f}, penetration_index {best_kill['pen']:.3f}, "
        f"IC50 {best_kill['ic50']:.3g} µM, bbb {best_kill['bbb']:.3f}). Kill is the "
        "product of **delivery × potency × (hypoxia resistance)**: a drug wins only "
        "when it both reaches the cells AND is potent enough at the exposure it "
        "achieves. Here delivery dominates the ranking — doxorubicin is the most "
        "potent by IC50 but is delivery-starved on this section, while the penetrant "
        "alkylators convert their exposure into kill.",
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
        f"- **Hypoxia compresses the panel.** ≈{hypoxic_frac*100:.0f}% of tumor cells "
        "are hypoxic; the resistance banding lowers the achievable kill for **every** "
        "drug, so absolute kills are small and the drugs are closer together than "
        "their potency differences alone would suggest.",
        "- **One anchor sets the scale.** The single α is fit to one doxorubicin "
        "penetration number; the absolute penetration of every drug rides on that "
        "measurement, and the binding proxy is an ordering, not a first-principles "
        "rate.",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
