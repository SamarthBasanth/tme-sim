"""REAL DRUG PANEL on the REAL Xenium breast-cancer microenvironment.

Breast-cancer analogue of ``experiments/real_drug_panel.py`` -- SAME engine,
SAME reconstruction entry point (``tmesim.reconstruct.real_pipeline.
reconstruct_real``), SAME calibration machinery (``tmesim.validation.
calibration.calibrate_alpha`` / ``calibrated_drugparams``), SAME binding-site
barrier (``tmesim.molecular.binding_site_barrier.apply_barrier_to_drugparams``).
Nothing in ``tmesim.core``, ``tmesim.physics``, ``tmesim.sim``, or the
GBM-specific calibration constants is touched -- this script only DRIVES the
existing, frozen engine on a different real tissue.

Data: Janesick et al. 2023 (Nat Commun 14, 8353; PMID 38114474) human
breast-cancer FFPE Xenium Replicate 1 (167,780 real cells, 313-gene targeted
panel), already cached at ``data_cache/xenium_breast_janesick/rep1`` by an
earlier task (``results/registration_3d_breast_poc.md``) -- read-only here, no
re-download.

Panel: [doxorubicin, paclitaxel, cyclophosphamide, tamoxifen, docetaxel,
epirubicin, 5-fluorouracil, vinorelbine, cisplatin, palbociclib] -- expanded
from the original 4-drug panel (see ``results/breast_cancer_flagship.md`` for
that panel's own writeup) to a 10-drug panel credible for pharma/biotech drug
comparison. docetaxel / epirubicin / 5-fluorouracil / vinorelbine are new
molecules in ``tmesim.molecular.breast_library`` with real, measured GDSC2
breast (BRCA) IC50s; cisplatin reuses the GBM library's molecule verbatim but
with a real, measured breast-specific IC50 (a real-data substitute for
carboplatin, which is absent from GDSC2 entirely -- see
``tmesim/molecular/breast_library.py`` for the full rejection rationale, also
covering capecitabine); palbociclib (CDK4/6 inhibitor) is a genuine
mechanism-category contrast case -- cytostatic, not cytotoxic -- included
deliberately to expose that boundary to a pharma audience rather than hidden.
Trastuzumab was RE-EVALUATED and remains excluded (no cytotoxic IC50 exists
for a biologic in GDSC or the literature this pipeline's Hill-kill model can
use).

  * doxorubicin -- REUSED VERBATIM from the GBM library
    (``tmesim.molecular.library.get_molecule`` +
    ``tmesim.validation.calibration.calibrated_drugparams("doxorubicin", alpha)``,
    the EXACT function call the GBM panel makes). Same molecule, same IC50,
    same calibrated physics -- this is the cross-cancer sanity-check drug: any
    difference in its behaviour here versus on GBM tissue comes ONLY from the
    breast tissue's own geometry/ECM/vasculature, not from a re-fit parameter.
  * paclitaxel / cyclophosphamide / tamoxifen -- compiled via
    ``tmesim.molecular.breast_library.get_molecule_breast`` (real, RDKit-
    verified SMILES) with a REAL measured IC50 from the live GDSC2 breast
    (BRCA) cell-line dose-response screen
    (``tmesim.molecular.potency_breast``, 50 cell lines per drug), run through
    the IDENTICAL calibration call the GBM panel uses
    (``calibrated_drugparams``, alpha fit ONCE on doxorubicin, never refit
    here) and the identical binding-site barrier.

Per-cell sensitivity: 'uniform' ONLY (every drug at 1.0, identical across
drugs -- the fair, delivery x potency comparison). The 'determinant' mode
(TOP2A / MKI67 mechanistic proxies,
``tmesim.reconstruct.real_xenium.SENSITIVITY_DETERMINANTS``) is GBM-panel
machinery this script does not extend to the breast drugs -- extending it
would mean asserting a mechanistic determinant for paclitaxel /
cyclophosphamide / tamoxifen that has not been checked here, so it is left
out rather than silently applied. Doxorubicin's own TOP2A determinant IS
still resolvable if a caller asks for 'determinant' mode (real_xenium.
SENSITIVITY_DETERMINANTS is untouched), but this script's headline table uses
'uniform' throughout, exactly matching the GBM panel's own headline choice.

Resolution: dx = 12 um, the SAME grid-converged resolution
``experiments/real_drug_panel.py`` uses (``results/numerical_rigor.md``): the
~45 um calibrated perivascular boundary layer needs dx <= 20 um to resolve, and
12 um is also the finest grid the upwind-advection CFL condition allows at
dt = 10 s. Running the breast panel at the SAME dx as the GBM panel (rather
than a coarser, cheaper one) is what makes the cross-cancer comparison in
Sec. "penetration vs GBM" a fair one.

Run:    ./.venv/bin/python experiments/breast_drug_panel.py
Writes: results/breast_drug_panel.png
        results/breast_drug_panel_data.json  (raw numbers, consumed by
                                               results/breast_cancer_flagship.md)
"""
from __future__ import annotations

import argparse
import json
import os
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from tmesim.core import PDEConfig, SimConfig
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.breast_library import get_molecule_breast
from tmesim.molecular.compile_drug import molecule_to_drugparams
from tmesim.molecular.library import get_molecule as get_molecule_gbm
from tmesim.molecular.potency_breast import resolve_ic50 as resolve_ic50_breast
from tmesim.physics.diffusion import ALPHA_TORTUOSITY
from tmesim.reconstruct.marker_panels import select_marker_map
from tmesim.reconstruct.real_pipeline import (
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import (
    DOX_TARGET_EFOLD_UM,
    calibrate_alpha,
    calibrated_drugparams,
)
from tmesim.validation.pk_reality import DEFAULT_ECM

PANEL = [
    "doxorubicin", "paclitaxel", "cyclophosphamide", "tamoxifen",
    "docetaxel", "epirubicin", "5-fluorouracil", "vinorelbine",
    "cisplatin", "palbociclib",
]

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
_RESULTS_DIR = os.path.join(_REPO, "results")
_OUT_PNG = os.path.join(_RESULTS_DIR, "breast_drug_panel.png")
_OUT_JSON = os.path.join(_RESULTS_DIR, "breast_drug_panel_data.json")

BREAST_SAMPLE_PATH = os.path.join(_REPO, "data_cache", "xenium_breast_janesick", "rep1")

# SAME grid-converged resolution the GBM panel uses (results/numerical_rigor.md).
DX_UM = 12.0


def compile_breast_drug(name: str, alpha: float):
    """Compile ONE drug for the breast panel through the SAME pipeline the GBM
    panel uses (calibrate perivascular sink -> apply binding-site barrier).

    doxorubicin: reuses ``calibrated_drugparams("doxorubicin", alpha)`` and
    ``library.get_molecule("doxorubicin")`` VERBATIM -- the identical function
    calls ``experiments/real_drug_panel.py`` makes.

    cisplatin: molecule REUSED VERBATIM from ``tmesim.molecular.library``
    (the GBM panel's own compiled cisplatin -- same SMILES, same structure),
    but compiled with a REAL, MEASURED GDSC2-BREAST (BRCA) IC50
    (``potency_breast.resolve_ic50``), not the GBM literature IC50 -- unlike
    doxorubicin, cisplatin IS screened in GDSC2's breast panel (n=50 cell
    lines), so the measured-tissue-specific number is used in preference to
    reusing GBM's own literature-tier value. This is the real-data
    substitute for carboplatin (evaluated and rejected -- absent from GDSC2
    entirely, across every tissue; see ``tmesim/molecular/breast_library.py``
    for the full rationale).

    paclitaxel / cyclophosphamide / tamoxifen / docetaxel / epirubicin /
    5-fluorouracil / vinorelbine / palbociclib: a fresh Molecule from
    ``breast_library`` compiled with a REAL GDSC2-breast IC50
    (``potency_breast.resolve_ic50``), then passed as an already-built
    ``DrugParams`` into ``calibrated_drugparams`` -- which recognises an
    already-compiled ``DrugParams`` and applies the SAME alpha-scaled
    perivascular-sink calibration to it that a library-name lookup would (see
    ``tmesim.validation.calibration._compile``'s passthrough branch). No
    calibration code is duplicated or reimplemented.

    palbociclib is a CDK4/6 inhibitor, predominantly CYTOSTATIC rather than
    cytotoxic -- included as a genuine category-mismatch contrast case (see
    ``tmesim/molecular/breast_library.py``'s module docstring); it goes
    through the IDENTICAL compilation and simulation pipeline as every other
    drug here (nothing is special-cased in the physics/kill model), but its
    kill_fraction output should be read with that mechanism caveat attached.
    """
    key = name.strip().lower()
    if key == "doxorubicin":
        dp = calibrated_drugparams("doxorubicin", alpha)
        mol = get_molecule_gbm("doxorubicin")
        ic50_info = {
            "ic50_uM": float(dp.ic50), "source": "literature",
            "confidence": "literature",
            "note": (
                "REUSED VERBATIM from the GBM panel "
                "(tmesim.molecular.potency.LITERATURE_IC50['doxorubicin']); "
                "doxorubicin is not screened in GDSC2 at all (neither the "
                "glioma nor the breast BRCA subset carries it)."
            ),
        }
    elif key == "cisplatin":
        mol = get_molecule_gbm("cisplatin")
        ic50_info = resolve_ic50_breast(key)
        dp_uncal = molecule_to_drugparams(
            mol, ecm_reference=DEFAULT_ECM, ic50=ic50_info["ic50_uM"]
        )
        dp = calibrated_drugparams(dp_uncal, alpha)
    else:
        mol = get_molecule_breast(key)
        ic50_info = resolve_ic50_breast(key)
        dp_uncal = molecule_to_drugparams(
            mol, ecm_reference=DEFAULT_ECM, ic50=ic50_info["ic50_uM"]
        )
        dp = calibrated_drugparams(dp_uncal, alpha)
    dp = apply_barrier_to_drugparams(dp, mol)
    return dp, mol, ic50_info


def effective_diffusivity(diffusivity: float, ecm_density: float) -> float:
    return float(diffusivity) * (1.0 - ALPHA_TORTUOSITY * float(np.clip(ecm_density, 0.0, 1.0)))


def main(dx_um: float = DX_UM):
    t0 = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    # --- 1. REAL breast microenvironment (Xenium Rep1), reconstructed ONCE ---
    print(f"Reconstructing the REAL Xenium breast-cancer microenvironment "
          f"(once, dx={dx_um} um) ...")
    micro = reconstruct_real(source="xenium", sample=BREAST_SAMPLE_PATH, dx_um=float(dx_um))
    report = reconstruction_report(micro)
    mean_ecm = report["mean_ecm_density"]
    _, panel_name, panel_cov = select_marker_map(
        # Re-derive which marker map was selected (reconstruction doesn't
        # surface it directly in the report); genes come from the same load
        # reconstruct_real already did internally, so re-load once here for
        # the provenance line only (cheap: ~2s, not in the simulated hot path).
        __import__("tmesim.data.real.singlecell_spatial", fromlist=["load"])
        .load(path=BREAST_SAMPLE_PATH)["genes"]
    )
    print(
        f"  {report['n_cells']:,} real cells  (tumor {report['n_tumor']:,})  "
        f"grid {report['grid']}  dx={report['dx_um']:.0f} um\n"
        f"  marker panel selected: {panel_name}  coverage={panel_cov}\n"
        f"  celltype fractions: {report['celltype_fractions']}\n"
        f"  mean vascular {report['mean_vascular_density']:.4f}  "
        f"ecm {mean_ecm:.4f}  o2 {report['mean_oxygen']:.4f}\n"
        f"  ecm markers used: {report['ecm_markers_used']}"
    )

    # --- 2. Single calibration constant, fit ONCE on doxorubicin (unchanged) ---
    alpha = calibrate_alpha()
    print(f"\nperivascular-sink calibration alpha = {alpha:.2f} "
          f"(fit once on doxorubicin, target L={DOX_TARGET_EFOLD_UM:.1f} um; "
          "REUSED verbatim, not refit on breast tissue)")

    # --- 3. Simulate every drug on the SAME real breast tissue -----------------
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)
    rows = []
    for name in PANEL:
        dp, mol, ic50_info = compile_breast_drug(name, alpha)
        info = apply_drug_sensitivity(micro, name, mode="uniform")
        cfg = SimConfig(drug=dp, pde=pde, seed=0)
        print(f"  simulating {name:16s} ...", end="", flush=True)
        t1 = time.time()
        res = simulate_response(micro, cfg, attribute=True)
        dt_run = time.time() - t1
        d_eff = effective_diffusivity(dp.diffusivity, mean_ecm)
        mech = res.mechanism  # dict: labels/counts/survival_breakdown/reconciliation
        rows.append(
            dict(
                name=name,
                mw=float(mol.get("mw")),
                ic50=float(dp.ic50),
                ic50_source=ic50_info["source"],
                ic50_confidence=ic50_info["confidence"],
                ic50_n_cell_lines=ic50_info.get("n_cell_lines"),
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
                mechanism_counts=(mech.get("counts") if mech is not None else None),
                mechanism_survival_breakdown=(
                    mech.get("survival_breakdown") if mech is not None else None
                ),
                mechanism_reconciled=(
                    bool(mech["reconciliation"]["bit_exact"]) if mech is not None else None
                ),
                runtime_s=dt_run,
            )
        )
        print(f" kill={rows[-1]['kill']:.4f}  pen={rows[-1]['pen']:.4f}  "
              f"mean_auc={rows[-1]['mean_auc']:.4g}  ic50={rows[-1]['ic50']:.4g} uM "
              f"({ic50_info['source']})  mech_ok={rows[-1]['mechanism_reconciled']}  "
              f"[{dt_run:.1f}s]")

    # --- 4. Figure: per-drug AUC heatmaps (shared scale) + bar chart ----------
    W, H = micro.domain.extent_um
    ext = [0, W, H, 0]
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    all_tissue_auc = np.concatenate([r["auc_field"][tissue] for r in rows])
    vmax = float(np.percentile(all_tissue_auc, 99.0)) or 1.0

    # Width scales with panel size so a 10-drug panel doesn't cram 10 AUC
    # heatmaps into a figure sized for the original 4-drug panel.
    fig = plt.figure(figsize=(max(17, 2.0 * len(PANEL)), 9))
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
        "Predicted kill & penetration on the SAME real breast tumor (Rep1)  "
        "-- PREDICTIONS, not validated outcomes (untreated resection)",
        fontsize=9.5,
    )
    fig.suptitle(
        "REAL drug panel on the REAL Xenium breast-cancer microenvironment "
        "(Janesick et al. 2023 Rep1)\n"
        f"grid dx = {report['dx_um']:.0f} µm (matches the GBM panel's converged "
        f"resolution) · α = {alpha:.1f} (SAME constant, fit on GBM doxorubicin, "
        "never refit) · sensitivity mode = 'uniform'",
        fontsize=12,
    )
    fig.savefig(_OUT_PNG, dpi=130, bbox_inches="tight")
    print(f"\nfigure written: {_OUT_PNG}")

    # --- 5. Raw data dump (consumed by results/breast_cancer_flagship.md) -----
    dump_rows = []
    for r in rows:
        d = {k: v for k, v in r.items() if k != "auc_field"}
        dump_rows.append(d)
    dump = {
        "reconstruction_report": {
            k: v for k, v in report.items()
            if k not in ("sensitivity", "oxygen", "cell_state")
        },
        "marker_panel_selected": panel_name,
        "marker_panel_coverage": panel_cov,
        "alpha": alpha,
        "dx_um": dx_um,
        "rows": dump_rows,
        "runtime_s": time.time() - t0,
    }
    with open(_OUT_JSON, "w") as fh:
        json.dump(dump, fh, indent=2, default=float)
    print(f"data written: {_OUT_JSON}")

    print("\n" + "=" * 92)
    print("REAL DRUG PANEL on the REAL Xenium breast-cancer microenvironment")
    print("=" * 92)
    hdr = (f"{'drug':16s} {'MW':>8} {'IC50(uM)':>10} {'src':>11} {'bbb':>6} "
           f"{'D_eff':>8} {'pen_idx':>8} {'mean_auc':>10} {'kill':>7}")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['name']:16s} {r['mw']:8.1f} {r['ic50']:10.4g} "
              f"{r['ic50_source']:>11} {r['bbb']:6.3f} {r['d_eff']:8.1f} "
              f"{r['pen']:8.4f} {r['mean_auc']:10.4g} {r['kill']:7.4f}")
    print("=" * 92)
    print(f"[runtime {time.time() - t0:.1f} s]")
    return dump


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Real drug panel on the real Xenium breast section")
    ap.add_argument("--dx", type=float, default=DX_UM,
                     help="grid spacing in um (default %(default)s, matches GBM panel)")
    main(dx_um=ap.parse_args().dx)
