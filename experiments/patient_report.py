"""THE USE CASE — a per-tumor VIRTUAL DRUG PANEL report.

Given ONE real patient's tissue, produce the decision-support artifact the whole
simulator exists to generate: a measured description of that patient's tumor
microenvironment, a ranked five-drug virtual panel simulated on it, a COMPUTED
attribution of why each drug lands where it does, and a decision figure — plus a
loud, unmissable statement of what this is and is not.

What is REAL here
-----------------
* the **geometry** — real measured cells at real positions, a real vascular
  field (KDE of measured ENDOTHELIAL cells), a real ECM field (expression-
  weighted KDE of the panel's collagen / matrix transcripts), and a
  physics-solved oxygen field driven by that real vasculature
  (``reconstruct.real_pipeline.reconstruct_real``);
* the **niche regions** — each cell's microenvironment zone assigned by
  correlating its local measured expression against the Ivy GAP anatomic-region
  signatures (``reconstruct.region_label``);
* the **drug physics** — structure-compiled transport (Stokes–Einstein radius →
  ECM-mesh diffusivity) plus a perivascular sink whose single constant α was fit
  ONCE on doxorubicin and revalidated out-of-sample on mitoxantrone
  (``validation.calibration``), then routed through the saturable binding-site
  barrier (``molecular.binding_site_barrier``);
* the **potency** — real measured GDSC glioma IC50 where it exists, curated
  literature IC50 where it does not, each labelled with its source
  (``molecular.potency``).

What is NOT real
----------------
The **outcomes**. This is an untreated open resection: there is no paired
measured drug response for any of these five drugs, on this tissue or any other
tissue in this repository. Every kill / exposure / ranking number below is a
PREDICTION of the physics, never a validated result. See the LIMITATIONS block
at the top and bottom of the generated report.

The attribution (computed, not asserted)
----------------------------------------
Predicted kill is driven by exposure relative to potency, then discounted by the
tumor's own oxygen-driven resistance (hypoxic AND quiescent cells). We decompose
each drug's standing into three **exactly additive log10 terms** ("decades"),
each centred on the panel mean so they read as advantage (+) / penalty (−)
versus the panel:

    C_tumor_i      = mean AUC sampled at this tumor's TUMOR-cell positions / T
    drive_i        = C_tumor_i / IC50_i                    (dimensionless)

    delivery_i     = log10 C_tumor_i - mean_j log10 C_tumor_j
    potency_i      = -( log10 IC50_i - mean_j log10 IC50_j )
    o2resist_i     = log10 rho_i     - mean_j log10 rho_j

      where rho_i = kill_i / kill_i(resistance OFF) is measured by a SECOND,
      otherwise-identical simulation of the same drug on the same tissue with
      ``PDConfig(hypoxia_resistance=1, quiescent_resistance=1)`` — i.e. the
      resistance term is a counterfactual simulation, not a hand-wave. Because
      BOTH penalties are switched off together, the term is named for the whole
      oxygen-driven resistance state rather than for hypoxia alone.

    delivery_i + potency_i  ==  log10 drive_i - mean_j log10 drive_j   (exact)

The "dominant driver" for a drug is simply the term of largest magnitude, and
the reason for each ADJACENT RANK GAP is the term that differs most between the
two drugs. We also report the Spearman rank correlation between the summed
attribution and the actual simulated kill ranking, so the attribution is itself
falsifiable rather than decorative.

Per-cell sensitivity mode (states which ranking this is)
--------------------------------------------------------
This report ranks drugs, so it runs the ONLY mode in which a cross-drug ranking
is fair: ``apply_drug_sensitivity(..., mode='uniform')`` — every drug gets
per-cell ``sensitivity`` 1.0, so the factor is identical across drugs, cancels
out of the comparison, and the ranking isolates DELIVERY x POTENCY (discounted
by the tissue's own oxygen-driven resistance, which is measured by the
counterfactual run). Absolute kills are therefore an equal
no-intrinsic-resistance UPPER BOUND for every drug.

The alternative ``mode='determinant'`` (each drug scaled by its own measured
panel determinant where one exists) is NOT used for ranking here, because the
319-gene panel carries determinants for some of these drugs and not others:
TOP2A for doxorubicin and an MKI67 proliferation proxy for cisplatin /
lomustine / methotrexate (both ~0.585 mean on this section), but NOTHING for
temozolomide, whose determinant MGMT is absent. In that mode TMZ alone keeps
1.0 and is flattered ~1.7x for a purely technical reason. Both rankings are
reported side by side in ``results/real_drug_panel.md``.

Tissue mask
-----------
``Domain.tissue_mask`` is now GEOMETRIC — the union of 30 µm discs around the
measured cell centroids — which makes the simulated tissue domain independent of
``dx``. The legacy "grid cells containing >= 1 segmented cell" mask shrank with
the mesh (50.6% / 32.7% / 13.9% coverage at dx = 30/20/12 µm) and, because the
vascular KDE is masked to tissue, made the total O2 supply a function of the
grid. Every hypoxia number this pipeline reported before that fix is retracted.

Resolution
----------
Per ``results/numerical_rigor.md``, ``kill_fraction`` is grid-converged for
dx <= 20 µm but ``penetration_index`` is NOT converged at dx = 30 µm (the
calibrated ~45 µm perivascular boundary layer is unresolved on coarse grids, so
penetration is over-estimated by ~15–25%). This report therefore runs at
dx = 12 µm, inside the converged regime.

Run:    ./.venv/bin/python experiments/patient_report.py
        ./.venv/bin/python experiments/patient_report.py --source visium
Writes: results/patient_report_<tissue>.md
        results/patient_report_<tissue>.png

ADD-ONLY: nothing in core / physics / molecular / reconstruct is modified.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:  # allow `python experiments/patient_report.py`
    sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..")))

from tmesim.core import CellType, OxygenConfig, PDConfig, PDEConfig, SimConfig
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.molecular.potency import resolve_ic50
from tmesim.physics.diffusion import ALPHA_TORTUOSITY
from tmesim.reconstruct.real_pipeline import (
    apply_drug_sensitivity,
    reconstruct_real,
    reconstruction_report,
)
from tmesim.reconstruct.region_label import (
    REGION_ORDER,
    label_regions,
    region_field,
    region_signatures,
)
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.validation.pk_reality import penetration_length_um
from tmesim.viz.patient_report import plot_patient_report

PANEL: List[str] = [
    "temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin",
]

_RESULTS_DIR = os.path.abspath(os.path.join(_HERE, "..", "results"))

# CONVERGED grid resolution (see results/numerical_rigor.md section 2).
DX_UM = 12.0
# Exposure window used for every drug (identical tissue, identical schedule).
WINDOW_S = 3600.0 * 6
# Neighbourhood pitch for aggregating single-cell Xenium counts before region
# scoring. Per-cell Xenium counts are far too sparse to correlate against a
# bulk (Ivy GAP) signature; ~100 um matches the Visium spot pitch the region
# labeller was designed around.
NEIGHBOURHOOD_UM = 100.0
_MIN_CELLS_PER_NEIGHBOURHOOD = 5

_TISSUE_LABEL = {
    "xenium": "10x Xenium single-cell glioblastoma section",
    "visium": "10x Visium spot-level glioblastoma section",
}

_LIMITATIONS_TITLE = "LIMITATIONS — READ BEFORE ANY OTHER NUMBER ON THIS PAGE"

# The ONLY per-cell-sensitivity mode in which a cross-drug ranking is fair: it is
# identical for every drug, so it cancels out of the comparison. See the module
# docstring and `reconstruct.real_pipeline.SENSITIVITY_MODES`.
SENSITIVITY_MODE = "uniform"


# --------------------------------------------------------------------------- #
# LIMITATIONS (rendered verbatim at the TOP and the BOTTOM of the report)
# --------------------------------------------------------------------------- #
def limitations_block(tissue_label: str, dx_um: float) -> str:
    # Do not claim convergence we did not run at: numerical_rigor.md puts the
    # penetration_index asymptote below ~15 um.
    if dx_um <= 15.0:
        dx_clause = (
            f"This report runs at **dx = {dx_um:.0f} µm**, inside the converged "
            "regime; any penetration number quoted at dx ≳ 20 µm should be read "
            "as an upper bound.")
    else:
        dx_clause = (
            f"⚠️ **This report was run at dx = {dx_um:.0f} µm, which is NOT in "
            "the converged regime** — its `penetration_index` is an upper bound "
            "and is likely ~15–25% high. Re-run at dx ≤ 12 µm for a converged "
            "penetration number.")
    return "\n".join([
        f"> ## ⚠️ {_LIMITATIONS_TITLE}",
        ">",
        "> **This is a research prototype. Every number in this report is a "
        "PREDICTION of a physical model. NOT FOR CLINICAL USE. Do not use it to "
        "select, rank, exclude, or dose any therapy for any patient.**",
        ">",
        "> 1. **The predictions are NOT clinically validated.** No prediction in "
        "this report has been compared against a measured treatment response, "
        "because none exists to compare against.",
        f"> 2. **No outcome-labelled data exists for this tissue.** The "
        f"{tissue_label} is an **untreated resection**. There is no paired drug "
        "response for any of the five drugs, so the ranking below has never been "
        "scored and its accuracy is **unknown** — not 'good', not 'bad', unknown.",
        "> 3. **Potency is measured in dissociated cells.** The IC50 values come "
        "from GDSC cell-line viability screens (dissociated glioma lines in "
        "monolayer) or from curated literature ranges. They carry none of this "
        "tumor's 3-D architecture, stromal/immune context, or intra-tumoral "
        "heterogeneity, and cell-line IC50 is a weak surrogate for in-patient "
        "sensitivity.",
        f"> 4. **`penetration_index` is resolution-sensitive.** Per "
        "`results/numerical_rigor.md` §2 it is *over-estimated on coarse grids* "
        "(≈0.285 at dx = 60 µm falling to ≈0.187 at dx = 10 µm) because a grid "
        "coarser than the ~45 µm perivascular boundary layer cannot resolve it. "
        f"{dx_clause} `kill_fraction` *is* grid-converged for dx ≤ 20 µm.",
        "> 5. **One measurement sets the absolute scale.** The perivascular sink "
        "constant α is fit to a *single* doxorubicin penetration measurement "
        "(a 40–50 µm half-distance, converted to the e-folding length "
        "L = d_half/ln2 ≈ 65 µm); every drug's absolute penetration rides on that "
        "one anchor, "
        "and the underlying binding proxy is a structure-based *ordering*, not a "
        "first-principles rate. A drug it mis-ranks (active transport, efflux, "
        "ion trapping) would be mis-predicted.",
        "> 6. **Per-cell drug sensitivity is an assumption, not a measurement — "
        f"and this ranking was produced in `mode='{SENSITIVITY_MODE}'`.** Every "
        "drug carries uniform `sensitivity` = 1.0, so no intrinsic resistance is "
        "asserted for any of them. That is deliberate: it is the only setting in "
        "which a *cross-drug ranking* is fair, because an identical factor "
        "cancels out of the comparison. The price is that every **absolute** kill "
        "below is a no-resistance **upper bound**. The alternative "
        "`mode='determinant'` scales each drug by its own measured panel "
        "determinant where one exists (TOP2A for doxorubicin, an MKI67 "
        "proliferation proxy for cisplatin / lomustine / methotrexate, both "
        "≈0.585 here) — but temozolomide's determinant **MGMT is absent from this "
        "panel**, so in that mode TMZ alone stays at 1.0 and gains ~1.7× on its "
        "competitors for a purely technical reason. That is an artifact, not "
        "efficacy; both rankings are printed side by side in "
        "`results/real_drug_panel.md`. No absolute kill fraction here is a "
        "measured per-cell drug response.",
        "> 8. **The tissue mask was grid-dependent and has been fixed; every "
        "earlier hypoxia number from this pipeline is retracted.** "
        "`Domain.tissue_mask` used to be \"grid cells containing ≥ 1 segmented "
        "cell\", whose coverage collapsed 50.6% → 32.7% → 13.9% at "
        "dx = 30 → 20 → 12 µm; because the vascular KDE is masked to tissue, the "
        "total O₂ *supply* shrank with the mesh. It is now the union of 30 µm "
        "discs around the measured centroids (coverage 67.2/67.4/67.7% over the "
        "same grids), so mean O₂ is grid-invariant to ~2%. The consequence is "
        "reported honestly in §1.3: the predicted hypoxic fraction is now 0%, "
        "which for a glioblastoma is *also* not credible — the un-calibrated "
        "`OxygenConfig.consumption` means the absolute oxygen level is unanchored "
        "in both directions. Nothing in this report should rest on it.",
        "> 7. **Systemic PK is simplified.** A single time-averaged sustained "
        "plasma level drives the transport solve (documented Stage-0 "
        "simplification in `tmesim/sim/engine.py`); it preserves total delivered "
        "drug and relative ranking but is not a time-resolved infusion.",
    ])


# --------------------------------------------------------------------------- #
# 1. microenvironment profile
# --------------------------------------------------------------------------- #
def _percentiles(field: np.ndarray, mask: np.ndarray) -> Dict[str, float]:
    v = np.asarray(field, dtype=np.float64)[np.asarray(mask, dtype=bool)]
    if v.size == 0:
        return {k: float("nan") for k in ("p10", "p50", "p90", "mean", "max")}
    return {
        "p10": float(np.percentile(v, 10)),
        "p50": float(np.percentile(v, 50)),
        "p90": float(np.percentile(v, 90)),
        "mean": float(v.mean()),
        "max": float(v.max()),
    }


def _neighbourhood_expression(xy: np.ndarray, expression, pitch_um: float):
    """Aggregate single-cell counts into square neighbourhoods of ``pitch_um``.

    Returns ``(agg_expr (n_bins, n_genes), centroids (n_bins, 2), n_cells_per
    bin, cell -> bin index)``. Aggregation is a plain sparse row-sum, so a bin's
    profile is the raw pooled count vector the region labeller expects.
    """
    from scipy import sparse

    xy = np.asarray(xy, dtype=np.float64)
    ix = np.floor((xy[:, 0] - xy[:, 0].min()) / pitch_um).astype(np.int64)
    iy = np.floor((xy[:, 1] - xy[:, 1].min()) / pitch_um).astype(np.int64)
    key = iy * (ix.max() + 1) + ix
    _uniq, inv = np.unique(key, return_inverse=True)
    n_bins = int(inv.max()) + 1

    picker = sparse.csr_matrix(
        (np.ones(inv.size), (inv, np.arange(inv.size))),
        shape=(n_bins, inv.size),
    )
    expr = expression if sparse.issparse(expression) else sparse.csr_matrix(expression)
    agg = (picker @ expr).tocsr()

    counts = np.bincount(inv, minlength=n_bins).astype(np.float64)
    safe = np.maximum(counts, 1.0)
    cx = np.bincount(inv, weights=xy[:, 0], minlength=n_bins) / safe
    cy = np.bincount(inv, weights=xy[:, 1], minlength=n_bins) / safe
    return agg, np.column_stack([cx, cy]), counts, inv


def _location_expression(source: str, sample: Optional[str]):
    """Measured expression at spatial LOCATIONS for the given real source.

    Xenium -> ~100 um pooled neighbourhoods of the measured single cells (per-cell
    counts are far too sparse to correlate against a bulk signature).
    Visium  -> the spots themselves, which are already the right granularity.

    Returns ``(expr, coords_um, genes, weights, granularity_note)`` or ``None``
    if the source's raw counts are unavailable.
    """
    src = source.strip().lower()
    try:
        if src == "xenium":
            from tmesim.data.real.singlecell_spatial import load as load_sc

            kwargs = {"dense": False}
            if sample is not None:
                kwargs["path"] = sample
            try:
                data = load_sc(**kwargs)
            except Exception:  # retry once on a transient IO/network hiccup
                time.sleep(3)
                data = load_sc(**kwargs)
            agg, cent, counts, _inv = _neighbourhood_expression(
                np.asarray(data["cell_xy_um"], dtype=np.float64),
                data["expression"], NEIGHBOURHOOD_UM,
            )
            keep = counts >= _MIN_CELLS_PER_NEIGHBOURHOOD
            note = (
                f"pooled {NEIGHBOURHOOD_UM:.0f} µm × {NEIGHBOURHOOD_UM:.0f} µm "
                f"neighbourhoods of ≥ {_MIN_CELLS_PER_NEIGHBOURHOOD} measured "
                "cells"
            )
            return (agg[keep], cent[keep], np.asarray(data["genes"], dtype=str),
                    counts[keep], note)

        if src == "visium":
            from tmesim.data.real.visium import load_sample

            sid = sample or "10x_parent_gbm"
            try:
                data = load_sample(sid)
            except Exception:
                time.sleep(3)
                data = load_sample(sid)
            expr = data["expression"]
            n = expr.shape[0]
            note = "measured in-tissue Visium spots (55 µm, 100 µm pitch)"
            return (expr, np.asarray(data["coords_um"], dtype=np.float64),
                    np.asarray(data["genes"], dtype=str), np.ones(n), note)
    except Exception as exc:  # noqa: BLE001 — the report degrades, never crashes
        print(f"  [niche regions unavailable: {type(exc).__name__}: {exc}]")
    return None


def niche_regions(micro, source: str, sample: Optional[str]) -> Dict:
    """Assign every measured location — and hence every reconstructed cell — to
    an Ivy GAP anatomic niche (margin / core / perivascular / necrotic).

    Locations are matched into the reconstruction's coordinate frame by nearest
    neighbour after removing each cloud's own origin offset, so this does not
    assume the loader and the reconstructor share a row order or an origin.

    Returns ``{}`` (and the report degrades gracefully) if the raw counts or the
    Ivy GAP reference are unavailable.
    """
    loaded = _location_expression(source, sample)
    if loaded is None:
        return {}
    expr, coords, genes, _w, granularity = loaded

    try:
        sigs = region_signatures()
        labels, _scores, _names, shared = label_regions(
            expr, coords, genes, signatures=sigs, return_scores=True)
    except Exception as exc:  # noqa: BLE001
        print(f"  [niche regions unavailable: {type(exc).__name__}: {exc}]")
        return {}

    # Put the locations in the domain frame (both clouds are origin-shifted).
    loc = np.asarray(coords, dtype=np.float64) - np.asarray(coords, dtype=np.float64).min(axis=0)
    domain = micro.domain
    W, H = domain.extent_um
    loc[:, 0] = np.clip(loc[:, 0], 0.0, W - 1e-6)
    loc[:, 1] = np.clip(loc[:, 1], 0.0, H - 1e-6)

    from scipy.spatial import cKDTree

    cell_xy = np.column_stack([np.asarray(micro.cells.x, dtype=np.float64),
                               np.asarray(micro.cells.y, dtype=np.float64)])
    _d, nn = cKDTree(loc).query(cell_xy, k=1)
    cell_region = np.asarray(labels, dtype=object).astype(str)[nn]

    tumor = np.asarray(micro.cells.cell_type, dtype=int) == int(CellType.TUMOR)
    present = [r for r in REGION_ORDER if (labels == r).any()]
    n_all = int(cell_region.size)
    n_tum = int(tumor.sum())
    fields = region_field(domain, loc, labels, bandwidth_um=NEIGHBOURHOOD_UM)

    return {
        "granularity": granularity,
        "n_locations": int(len(labels)),
        "n_shared_genes": int(len(shared)),
        "regions_present": present,
        "location_counts": {r: int((labels == r).sum()) for r in present},
        "cell_fractions": {
            r: (float((cell_region == r).sum()) / n_all if n_all else 0.0)
            for r in present},
        "tumor_fractions": {
            r: (float((tumor & (cell_region == r)).sum()) / n_tum if n_tum else 0.0)
            for r in present},
        "region_fields": fields,
        "cell_region": cell_region,
    }


def resistance_state_profile(micro) -> Dict[str, float]:
    """Fraction of tumor cells in each oxygen-banded resistance state.

    Reproduces exactly what the engine does at its oxygen step — solve the
    steady-state O2 field, then band every living cell VIABLE / QUIESCENT /
    HYPOXIC — on a private copy, using the public physics and cell modules. This
    matters because the resistance counterfactual in §3 switches off *both* the
    hypoxic and the quiescent penalty, so the report must state both fractions
    rather than attributing the whole discount to hypoxia.
    """
    from tmesim.cells.agents import state_fractions, update_cell_states
    from tmesim.core import Microenvironment, TissueFields
    from tmesim.physics.oxygen import solve_oxygen

    # Private cells AND a private field container: never mutate the caller's.
    work = Microenvironment(
        domain=micro.domain,
        fields=TissueFields(ecm_density=micro.fields.ecm_density,
                            vascular_density=micro.fields.vascular_density),
        cells=micro.cells.copy(), patient_id=micro.patient_id,
        metadata=micro.metadata)
    update_cell_states(work, solve_oxygen(work, OxygenConfig()), OxygenConfig())
    return state_fractions(work)


def microenvironment_profile(micro, report: Dict, regions: Dict) -> Dict:
    """Measured description of THIS patient's tumor (no drug involved yet)."""
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    f = micro.fields
    return {
        "states": resistance_state_profile(micro),
        "vascular": _percentiles(f.vascular_density, tissue),
        "ecm": _percentiles(f.ecm_density, tissue),
        "oxygen": _percentiles(f.oxygen, tissue) if f.oxygen is not None else None,
        "n_tissue_voxels": int(tissue.sum()),
        "tissue_fraction": float(tissue.mean()),
        "regions": regions,
    }


# --------------------------------------------------------------------------- #
# 2. simulate the panel on THIS tumor
# --------------------------------------------------------------------------- #
def compile_real_drug(name: str, alpha: float):
    """Compile ONE real drug through the full, unmodified pipeline: real GDSC /
    literature potency + calibrated perivascular sink + saturable binding-site
    barrier. Identical to the pipeline every other tmesim experiment uses."""
    dp = calibrated_drugparams(name, alpha)
    mol = get_molecule(name)
    return apply_barrier_to_drugparams(dp, mol), mol


def effective_diffusivity(diffusivity: float, ecm_density: float) -> float:
    """Tissue diffusivity the PDE actually experiences after ECM tortuosity:
    ``D_eff = D · (1 − ALPHA_TORTUOSITY · ecm)`` (µm²/s)."""
    return float(diffusivity) * (
        1.0 - ALPHA_TORTUOSITY * float(np.clip(ecm_density, 0.0, 1.0)))


def _sample_field_at_cells(field: np.ndarray, micro, mask: np.ndarray) -> float:
    """Mean of a grid field sampled at the positions of the selected cells."""
    if not mask.any():
        return float("nan")
    col, row = micro.domain.to_index(np.asarray(micro.cells.x)[mask],
                                    np.asarray(micro.cells.y)[mask])
    return float(np.asarray(field, dtype=np.float64)[row, col].mean())


def simulate_panel(micro, alpha: float, mean_ecm: float,
                   drugs: Sequence[str] = PANEL, dx_um: float = DX_UM) -> List[Dict]:
    """Simulate every drug on the SAME reconstructed tumor, twice each:

    1. the real run (this tumor's hypoxia + quiescence resistance active), and
    2. a COUNTERFACTUAL run with ``PDConfig(hypoxia_resistance=1,
       quiescent_resistance=1)`` — identical geometry, identical transport,
       identical seed — which isolates how much kill this tumor's hypoxic burden
       is costing each drug.
    """
    pde = PDEConfig(dt=10.0, total_time=WINDOW_S, convection=True, save_every=500)
    no_resist = PDConfig(hypoxia_resistance=1.0, quiescent_resistance=1.0)
    tumor_mask = np.asarray(micro.cells.cell_type, dtype=int) == int(CellType.TUMOR)

    rows: List[Dict] = []
    for name in drugs:
        dp, mol = compile_real_drug(name, alpha)
        pot = resolve_ic50(name)
        # EXPLICIT sensitivity mode. 'uniform' is identical for every drug, so it
        # cancels out of the ranking (see the module docstring); it is set per
        # drug anyway so the choice is visible at the call site, never implicit.
        sens_info = apply_drug_sensitivity(micro, name, mode=SENSITIVITY_MODE)
        print(f"  simulating {name:14s} ...", end="", flush=True)

        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
        cf = simulate_response(
            micro, SimConfig(drug=dp, pde=pde, seed=0, pd=no_resist))

        auc = np.asarray(res.fields["auc"], dtype=np.float64)
        tumor_auc = _sample_field_at_cells(auc, micro, tumor_mask)
        rows.append(dict(
            name=name,
            mw=float(mol.get("mw")),
            ic50=float(dp.ic50),
            ic50_source=str(pot["source"]),
            ic50_confidence=str(pot["confidence"]),
            ic50_note=str(pot["note"]),
            ic50_n_lines=pot.get("n_cell_lines"),
            bbb=float(dp.bbb_penetrance),
            diffusivity=float(dp.diffusivity),
            d_eff=effective_diffusivity(dp.diffusivity, mean_ecm),
            binding=float(dp.binding),
            clearance=float(dp.clearance),
            hill=float(dp.hill),
            # Patient-specific analytic penetration length at THIS tumor's ECM.
            L_um=float(penetration_length_um(dp, ecm_density=mean_ecm)),
            kill=float(res.kill_fraction),
            kill_no_resistance=float(cf.kill_fraction),
            pen=float(res.extras["penetration_index"]),
            mean_auc=float(res.extras["mean_auc"]),
            tumor_auc=float(tumor_auc),
            hypoxic=float(res.extras["hypoxic_tumor_fraction"]),
            auc_field=auc,
            oxygen_field=np.asarray(res.fields["oxygen"], dtype=np.float64),
            sens_mode=str(sens_info["mode"]),
            sens_mean=float(sens_info["mean"]),
            sens_gene=sens_info["determinant_gene"],
            sens_upper_bound=bool(sens_info["upper_bound"]),
        ))
        r = rows[-1]
        print(f" kill={r['kill']:.4f} (no-resist {r['kill_no_resistance']:.4f})"
              f"  L={r['L_um']:7.1f} µm  pen={r['pen']:.3f}"
              f"  tumor AUC={r['tumor_auc']:.3g}")
    return rows


# --------------------------------------------------------------------------- #
# 3. attribution — computed, additive, falsifiable
# --------------------------------------------------------------------------- #
_TERMS = ("delivery_decades", "potency_decades", "hypoxia_decades")
_TERM_LABEL = {
    "delivery_decades": "DELIVERY",
    "potency_decades": "POTENCY",
    # The counterfactual switches off BOTH the hypoxic and the quiescent kill
    # penalty, so the honest name for this term is the whole oxygen-driven
    # resistance state, not hypoxia alone.
    "hypoxia_decades": "O2-RESISTANCE",
}


def _pct_rank(values: Sequence[float], i: int, higher_is_better: bool = True) -> float:
    """Percentile of item ``i`` within the panel (0 = worst, 100 = best)."""
    v = np.asarray(values, dtype=np.float64)
    n = v.size
    if n < 2:
        return 100.0
    order = np.argsort(v if higher_is_better else -v)
    rank = int(np.where(order == i)[0][0])
    return 100.0 * rank / (n - 1)


def _spearman(a: Sequence[float], b: Sequence[float]) -> float:
    from scipy.stats import rankdata

    ra, rb = rankdata(a), rankdata(b)
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def attribute(rows: List[Dict], window_s: float = WINDOW_S) -> Dict:
    """Attach the exact additive log10 attribution to every row, in place.

    ``delivery + potency`` reproduces the drug's log10 exposure-to-potency drive
    relative to the panel mean *exactly*; ``hypoxia`` is measured by the
    counterfactual resistance-off simulation. Returns the method's own
    faithfulness check.
    """
    eps = 1e-30
    c_tumor = np.array([max(r["tumor_auc"], eps) / window_s for r in rows])
    ic50 = np.array([max(r["ic50"], eps) for r in rows])
    rho = np.array([
        (r["kill"] / r["kill_no_resistance"]) if r["kill_no_resistance"] > 0 else 1.0
        for r in rows])
    rho = np.clip(rho, eps, None)

    log_c, log_ic50, log_rho = np.log10(c_tumor), np.log10(ic50), np.log10(rho)
    delivery = log_c - log_c.mean()
    potency = -(log_ic50 - log_ic50.mean())
    hypoxia = log_rho - log_rho.mean()

    kills = [r["kill"] for r in rows]
    for i, r in enumerate(rows):
        r["c_tumor_uM"] = float(c_tumor[i])
        r["drive"] = float(c_tumor[i] / ic50[i])
        r["rho"] = float(rho[i])
        r["hypoxia_cost"] = float(1.0 - rho[i])
        r["delivery_decades"] = float(delivery[i])
        r["potency_decades"] = float(potency[i])
        r["hypoxia_decades"] = float(hypoxia[i])
        r["total_decades"] = float(delivery[i] + potency[i] + hypoxia[i])
        r["exposure_pct"] = _pct_rank(c_tumor, i, higher_is_better=True)
        r["potency_pct"] = _pct_rank(ic50, i, higher_is_better=False)
        r["hypoxia_pct"] = _pct_rank(rho, i, higher_is_better=True)
        r["penetration_pct"] = _pct_rank([x["L_um"] for x in rows], i, True)
        dom = max(_TERMS, key=lambda k: abs(r[k]))
        # Name the term AND its direction — "DELIVERY" alone is ambiguous
        # between "wins on delivery" and "is crippled by delivery".
        r["dominant"] = (f"{_TERM_LABEL[dom]} "
                         f"({'advantage' if r[dom] >= 0 else 'penalty'})")
        r["dominant_term"] = _TERM_LABEL[dom]
        r["dominant_key"] = dom

    totals = [r["total_decades"] for r in rows]
    return {
        "spearman_total_vs_kill": _spearman(totals, kills),
        "spearman_delivery_vs_kill": _spearman(
            [r["delivery_decades"] for r in rows], kills),
        "spearman_potency_vs_kill": _spearman(
            [r["potency_decades"] for r in rows], kills),
        "panel_mean_hypoxia_cost": float(np.mean([r["hypoxia_cost"] for r in rows])),
        "window_s": float(window_s),
    }


def _gap_reason(better: Dict, worse: Dict, *, from_below: bool = False) -> str:
    """Which attribution term actually explains an adjacent rank gap.

    Deltas are always computed as ``better - worse`` (the higher-ranked drug
    minus the lower-ranked one); ``from_below`` only changes how the sentence is
    phrased, so the last-placed drug is described relative to the drug ABOVE it
    rather than nonsensically against itself.
    """
    deltas = {k: better[k] - worse[k] for k in _TERMS}
    key = max(deltas, key=lambda k: abs(deltas[k]))
    fold = 10.0 ** abs(deltas[key])
    if from_below:
        # Restate every delta from the LOWER-ranked drug's point of view, so the
        # signs in one sentence all share a single reference drug.
        others = ", ".join(f"{_TERM_LABEL[k].lower()} {-deltas[k]:+.2f}"
                           for k in _TERMS if k != key)
        return (
            f"Versus **{better['name']}** immediately above it, the term that "
            f"differs most is **{_TERM_LABEL[key]}** ({-deltas[key]:+.2f} "
            f"decades for this drug ≈ {fold:.1f}× "
            f"{'against it' if deltas[key] > 0 else 'in its favour'}); its "
            f"other terms stand at {others} relative to that drug."
        )
    others = ", ".join(
        f"{_TERM_LABEL[k].lower()} {deltas[k]:+.2f}" for k in _TERMS if k != key)
    return (
        f"Versus **{worse['name']}** immediately below it, the term that differs "
        f"most is **{_TERM_LABEL[key]}** ({deltas[key]:+.2f} decades ≈ "
        f"{fold:.1f}× {'in its favour' if deltas[key] > 0 else 'against it'}); "
        f"the other terms differ by {others}."
    )


# --------------------------------------------------------------------------- #
# 4. the markdown report
# --------------------------------------------------------------------------- #
def _ic50_provenance(r: Dict) -> str:
    if r["ic50_confidence"] == "measured":
        n = r.get("ic50_n_lines")
        return f"GDSC **measured** (n={n} glioma lines)"
    if r["ic50_confidence"] == "literature":
        return "curated **literature**"
    return "flagged **default (LOW)**"


def _fmt_pct(v: float) -> str:
    return f"{v:.0f}ᵗʰ" if not math.isnan(v) else "—"


def _drug_narrative(rank: int, r: Dict, rows: List[Dict], hypoxic_frac: float,
                    mean_ecm: float) -> List[str]:
    n = len(rows)
    lines = [
        f"### #{rank} · {r['name']} — dominant driver: **{r['dominant']}**",
        "",
        f"`predicted kill {r['kill']:.4f}` · "
        f"`tumor-cell exposure {r['tumor_auc']:.3g} µM·s "
        f"(C̄ = {r['c_tumor_uM']:.3g} µM)` · "
        f"`L = {r['L_um']:.0f} µm` · "
        f"`IC50 = {r['ic50']:.3g} µM` · "
        f"`BBB {r['bbb']:.3f}` · "
        f"`D_eff {r['d_eff']:.0f} µm²/s`",
        "",
        f"- **Delivery {r['delivery_decades']:+.2f} decades.** The mean drug "
        f"exposure actually sampled at this tumor's own tumor-cell positions is "
        f"{r['tumor_auc']:.3g} µM·s over the {WINDOW_S / 3600:.0f} h window "
        f"(time-averaged C̄ = {r['c_tumor_uM']:.3g} µM), the "
        f"{_fmt_pct(r['exposure_pct'])} percentile of the {n}-drug panel. That "
        f"follows from D_eff = {r['d_eff']:.0f} µm²/s at this tumor's mean ECM "
        f"({mean_ecm:.3f}) and BBB penetrance {r['bbb']:.3f}, giving a "
        f"penetration length L = {r['L_um']:.0f} µm "
        f"({r['L_um'] / 45.0:.1f}× the calibrated ~45 µm perivascular boundary "
        f"layer, {_fmt_pct(r['penetration_pct'])} percentile).",
        f"- **Potency {r['potency_decades']:+.2f} decades.** IC50 = "
        f"{r['ic50']:.3g} µM — {_ic50_provenance(r)} — the "
        f"{_fmt_pct(r['potency_pct'])} percentile of the panel.",
        f"- **O₂-resistance {r['hypoxia_decades']:+.2f} decades.** "
        f"Counterfactual check: re-simulating this same drug on this same tissue "
        f"with the hypoxic *and* quiescent kill penalties switched OFF raises "
        f"the predicted kill from {r['kill']:.4f} to "
        f"{r['kill_no_resistance']:.4f} — so this tumor's oxygen-driven "
        f"resistance state (§1.3) costs this drug "
        f"**{r['hypoxia_cost'] * 100:.1f}%** of the kill it would otherwise "
        f"achieve ({_fmt_pct(r['hypoxia_pct'])} percentile — higher is less "
        f"penalised). Note this term is nearly identical across the panel: it "
        f"depresses every drug rather than re-ordering them.",
    ]

    ordered = sorted(_TERMS, key=lambda k: abs(r[k]), reverse=True)
    a, b = ordered[0], ordered[1]
    lines.append(
        f"- **Reading.** Its position is set mainly by **{_TERM_LABEL[a]}** "
        f"({r[a]:+.2f} decades), with {_TERM_LABEL[b].lower()} "
        f"({r[b]:+.2f}) second; the three terms net to "
        f"{r['total_decades']:+.2f} decades versus the panel average."
    )
    if rank < n:
        lines.append(f"- **Why it ranks above the next drug.** "
                     f"{_gap_reason(r, rows[rank])}")
    else:
        lines.append(f"- **Why it ranks last.** "
                     f"{_gap_reason(rows[rank - 2], r, from_below=True)}")
    lines.append("")
    return lines


def build_markdown(*, tissue_key: str, tissue_label: str, report: Dict,
                   profile: Dict, rows: List[Dict], check: Dict, alpha: float,
                   hypoxic_frac: float, mean_ecm: float, dx_um: float,
                   png_name: str, runtime_s: float) -> str:
    frac = report["celltype_fractions"]
    regions = profile.get("regions") or {}
    states = profile.get("states") or {}
    thr = OxygenConfig().hypoxia_threshold
    best = rows[0]
    lim = limitations_block(tissue_label, dx_um)

    L: List[str] = [
        f"# Virtual drug panel — one real {tissue_label}",
        "",
        f"**Patient / sample id:** `{report['patient_id']}` · "
        f"**source:** `{report['source']}` · "
        f"**grid:** {report['grid'][0]}×{report['grid'][1]} at "
        f"dx = {report['dx_um']:.0f} µm · "
        f"**extent:** {report['extent_um'][0]:.0f} × "
        f"{report['extent_um'][1]:.0f} µm",
        "",
        f"**Generated by:** `experiments/patient_report.py` "
        f"(`./.venv/bin/python experiments/patient_report.py "
        f"--source {tissue_key}`) · **figure:** `results/{png_name}`",
        "",
        f"**Run provenance:** perivascular-sink calibration **α = {alpha:.1f}** "
        "(e-folding-length convention `L = d_half/ln 2`; the retired α = 309.6 "
        "fitted the raw half-distance as if it were an e-fold length and is "
        f"superseded) · per-cell sensitivity **`mode='{SENSITIVITY_MODE}'`** "
        "(identical for every drug → the cross-drug ranking is fair; absolute "
        "kills are a no-resistance upper bound) · tissue mask **geometric** "
        f"(union of {report.get('tissue_disc_radius_um') or 30:.0f} µm discs "
        "around the measured centroids, "
        f"{(report.get('tissue_coverage_fraction') or 0.0) * 100:.1f}% coverage) "
        "— **dx-independent**, replacing the old grid-dependent occupancy mask.",
        "",
        lim,
        "",
        "---",
        "",
        "## 0. What this artifact is",
        "",
        "One real patient's tumor tissue in, one **ranked, mechanistically "
        "attributed drug panel** out. The microenvironment is reconstructed "
        "**once** from measured spatial transcriptomics and reused for every "
        "drug, so every difference in the table below is attributable to that "
        "drug's own physics and potency acting on *identical* tissue — not to a "
        "different patient, a different geometry, or a different random draw.",
        "",
        "| what | is it real? | provenance |",
        "|---|:--:|---|",
        f"| cell positions & phenotypes | **real** | {report['n_cells']:,} "
        f"measured cells, threshold-gated marker typing |",
        "| vascular field | **real** | KDE of the measured ENDOTHELIAL cells |",
        f"| ECM field | **real** | expression-weighted KDE of measured matrix "
        f"transcripts ({', '.join(report['ecm_markers_used']) or 'stromal proxy'}) |",
        "| oxygen field | **real-driven** | physics solver driven by the real "
        "vascular field |",
        "| niche regions | **real** | Ivy GAP anatomic-region signature "
        "correlation |",
        "| drug transport | **real physics** | structure-compiled diffusivity + "
        f"α-calibrated perivascular sink (α = {alpha:.1f}, fit once on "
        "doxorubicin, revalidated out-of-sample on mitoxantrone) |",
        "| potency (IC50) | **real, but ex-vivo** | GDSC glioma cell-line "
        "screen, else curated literature — labelled per drug below |",
        "| **treatment outcome** | ❌ **absent** | this tissue is an untreated "
        "resection; there is **no** measured response to score against |",
        "",
        "---",
        "",
        "## 1. Tumor microenvironment profile — the measured description of THIS tumor",
        "",
        "### 1.1 Cell-type composition",
        "",
        "| cell type | cells | fraction |",
        "|---|---:|---:|",
    ]
    counts = report["celltype_counts"]
    for ct in CellType:
        if counts.get(ct.name, 0) == 0:
            continue
        L.append(f"| {ct.name.title()} | {counts[ct.name]:,} | "
                 f"{frac[ct.name] * 100:.2f}% |")
    L += [
        f"| **total** | **{report['n_cells']:,}** | 100% |",
        "",
        f"Tumor cells: **{report['n_tumor']:,}** "
        f"({frac['TUMOR'] * 100:.1f}% of all measured cells).",
        "",
        "### 1.2 Physical fields (over the tissue mask)",
        "",
        "| field | mean | p10 | median | p90 | max |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for key, label in (("vascular", "vascular density"),
                       ("ecm", "ECM volume fraction"),
                       ("oxygen", "oxygen (relative)")):
        p = profile.get(key)
        if not p:
            continue
        L.append(f"| {label} | **{p['mean']:.3f}** | {p['p10']:.3f} | "
                 f"{p['p50']:.3f} | {p['p90']:.3f} | {p['max']:.3f} |")
    L += [
        "",
        f"Tissue mask covers {profile['n_tissue_voxels']:,} grid voxels "
        f"({profile['tissue_fraction'] * 100:.1f}% of the bounding box).",
        "",
        "### 1.3 Hypoxia and the resistance state of the tumor cells",
        "",
        f"> **{hypoxic_frac * 100:.1f}% of this tumor's cells are hypoxic** "
        f"(O₂ below the {thr:g} threshold in the real-perfusion oxygen solve).",
        "",
        "The oxygen solve bands every living tumor cell, and *both* the hypoxic "
        "and the marginal-oxygen quiescent band carry a kill penalty — so the "
        "full resistance picture, not just hypoxia, is what suppresses the panel:",
        "",
        "| tumor-cell state | O₂ band | fraction | kill penalty |",
        "|---|---|---:|---|",
        f"| VIABLE | O₂ ≥ {2 * thr:g} | {states.get('viable', 0.0) * 100:.1f}% "
        "| none |",
        f"| QUIESCENT | {thr:g} ≤ O₂ < {2 * thr:g} | "
        f"{states.get('quiescent', 0.0) * 100:.1f}% | ×"
        f"{PDConfig().quiescent_resistance:g} |",
        f"| HYPOXIC | O₂ < {thr:g} | {states.get('hypoxic', 0.0) * 100:.1f}% "
        f"| ×{PDConfig().hypoxia_resistance:g} |",
        "",
        f"So although only {hypoxic_frac * 100:.1f}% of tumor cells are "
        f"outright hypoxic, "
        f"{(states.get('hypoxic', 0.0) + states.get('quiescent', 0.0)) * 100:.1f}%"
        " carry *some* oxygen-driven resistance. §3 measures what that costs "
        "each drug by re-simulating with both penalties switched off.",
        "",
        "> **Honest caveat on this number — read it before using it.** The "
        "hypoxic fraction is a threshold-crossing count sitting close to the O₂ "
        "cutoff, and this simulator does not currently predict it reliably.",
        ">",
        "> *What changed.* `Domain.tissue_mask` used to be \"grid cells "
        "containing ≥ 1 segmented cell\", whose coverage collapsed 50.6% → 32.7% "
        "→ 13.9% at dx = 30 → 20 → 12 µm. Because `real_fields._kde_on_grid` "
        "masks the vascular field to tissue, the total O₂ **supply** shrank as "
        "the mesh was refined, and the hypoxic fraction moved with it "
        "(≈94% → 36% → 16% → 13% → 11% → 7% over dx = 30 → 8 µm). **That entire "
        "ladder is a discretisation artifact and is retracted.** The mask is now "
        "*geometric* — the union of 30 µm discs around the measured centroids — "
        "covering 67.2/67.4/67.7% at dx = 30/20/12 µm, so mean O₂ at the measured "
        "cell positions is now invariant (0.570 / 0.565 / 0.558, a 1.02× spread) "
        "where it previously ran 0.522 / 0.426 / 0.251 (2.08×).",
        ">",
        "> *What that leaves.* With the supply no longer artificially starved, "
        "**no tumor cell falls below the 0.2 threshold at any dx — the predicted "
        "hypoxic fraction is 0.0%.** For a glioblastoma, a tumor type partly "
        "*defined* by pseudopalisading necrosis, that is not credible either. "
        "`OxygenConfig.consumption` has never been calibrated against a measured "
        "oxygen tension, so the absolute level is unanchored in both directions: "
        "the old numbers were too high for a numerical reason, and this one is "
        "plausibly too low for a calibration reason. **The honest position is "
        "that the hypoxic fraction of this section is currently unknown to this "
        "model.** What the fix does establish is that the answer no longer "
        "depends on the mesh. `kill_fraction` is grid-converged; this is not "
        "calibrated.",
        "",
        "### 1.4 Niche-region breakdown",
        "",
    ]
    if regions:
        L += [
            f"Each of {regions['n_locations']:,} measured locations "
            f"({regions['granularity']}) is scored against the four Ivy GAP "
            f"anatomic-region signatures over "
            f"{regions['n_shared_genes']} shared genes and assigned to its "
            "best-correlating niche; every reconstructed cell inherits its "
            "nearest location's niche.",
            "",
            "| niche region | locations | % of all cells | % of TUMOR cells |",
            "|---|---:|---:|---:|",
        ]
        for reg in regions["regions_present"]:
            L.append(
                f"| {reg} | {regions['location_counts'][reg]:,} | "
                f"{regions['cell_fractions'][reg] * 100:.1f}% | "
                f"{regions['tumor_fractions'][reg] * 100:.1f}% |")
        L += [
            "",
            "*Read this honestly:* the niche assignment is a correlation to a "
            "bulk reference atlas (Ivy GAP) over the small set of genes shared "
            "with this assay's panel — it is a regional *prior*, not a "
            "histological read of this section.",
            "",
        ]
    else:
        L += ["Niche-region labelling was unavailable for this sample "
              "(raw counts or the Ivy GAP reference could not be loaded).", ""]

    # ---- section 2: the panel ------------------------------------------- #
    L += [
        "---",
        "",
        "## 2. Ranked virtual drug panel",
        "",
        f"> ### 🏁 Highest predicted kill on this tumor: **{best['name']}** "
        f"({best['kill']:.4f})",
        f"> Dominant driver **{best['dominant']}** · L = {best['L_um']:.0f} µm · "
        f"IC50 {best['ic50']:.3g} µM — {_ic50_provenance(best)} · "
        f"tumor-cell exposure {best['tumor_auc']:.3g} µM·s.",
        "> ",
        "> **This is a model output, not a recommendation.** It has never been "
        "checked against a measured response for this tissue or any other. "
        "Treat it as a hypothesis to test in the lab, nothing more.",
        "",
        "All five drugs are simulated on the **same** reconstruction with the "
        f"same {WINDOW_S / 3600:.0f} h exposure window, the same seed, and the "
        "same compilation pipeline "
        "(`calibrate_alpha` → `calibrated_drugparams` → "
        "`apply_barrier_to_drugparams`). Ranked by predicted kill fraction:",
        "",
        f"> **This ranking was produced in per-cell sensitivity "
        f"`mode='{SENSITIVITY_MODE}'`** — every drug carries sensitivity 1.0, so "
        "the factor is identical across drugs and cancels out of the comparison. "
        "That is what makes the ordering a fair delivery × potency contest; it "
        "also makes each **absolute** kill a no-resistance upper bound. The "
        "drug-specific `determinant` mode would scale doxorubicin (TOP2A) and "
        "cisplatin / lomustine / methotrexate (MKI67 proxy) down to ≈0.585 while "
        "leaving temozolomide at 1.0 — because TMZ's determinant MGMT is not on "
        "this panel — which flatters TMZ ~1.7× for a measurement-coverage reason "
        "rather than a pharmacological one. Both rankings are tabulated in "
        "[`results/real_drug_panel.md`](real_drug_panel.md).",
        "",
        "| # | drug | MW (Da) | real IC50 (µM) | IC50 source | BBB penetrance | "
        "D_eff (µm²/s) | penetration L (µm) | tumor-cell exposure AUC (µM·s) | "
        "tissue exposure AUC (µM·s) | penetration_index | **predicted kill** |",
        "|---:|---|---:|---:|:--|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for i, r in enumerate(rows, start=1):
        src = ("GDSC (measured)" if r["ic50_confidence"] == "measured"
               else "literature" if r["ic50_confidence"] == "literature"
               else "default (LOW)")
        L.append(
            f"| {i} | **{r['name']}** | {r['mw']:.0f} | {r['ic50']:.3g} | {src} | "
            f"{r['bbb']:.3f} | {r['d_eff']:.0f} | {r['L_um']:.0f} | "
            f"{r['tumor_auc']:.3g} | {r['mean_auc']:.3g} | {r['pen']:.3f} | "
            f"**{r['kill']:.4f}** |")

    L += [
        "",
        "**Column definitions.**",
        "",
        "- **real IC50 / IC50 source** — measured median across GDSC glioma cell "
        "lines where the drug is screened, otherwise a curated literature "
        "anchor. Per-drug provenance:",
        "",
    ]
    for r in rows:
        if r["ic50_confidence"] == "measured":
            tag = f"**GDSC measured**, n={r['ic50_n_lines']} glioma cell lines"
        elif r["ic50_confidence"] == "literature":
            tag = "**curated literature** (not measured in GDSC's glioma screen)"
        else:
            tag = "**flagged default — LOW confidence**"
        L.append(f"  - `{r['name']}` — {tag}. {r['ic50_note']}")
    L += [
        "",
        "- **BBB penetrance** — structure-derived fraction of the vascular "
        "concentration that crosses into tissue.",
        f"- **D_eff** — `D · (1 − {ALPHA_TORTUOSITY}·ECM)` evaluated at *this "
        f"tumor's* mean ECM = {mean_ecm:.3f}, i.e. the diffusivity the PDE "
        "actually experiences here.",
        "- **penetration L** — analytic reaction–diffusion length "
        "`L = √(D_eff / k_eff)` (µm) at this tumor's ECM; the grid-independent "
        "delivery metric. Compare against the calibrated ~45 µm perivascular "
        "boundary layer.",
        "- **tumor-cell exposure AUC** — the exposure field sampled at this "
        "tumor's actual tumor-cell coordinates (what the cells that matter "
        "really see). **tissue exposure AUC** is the plain mean over the tissue "
        "mask.",
        "- **penetration_index** — the solver's spatial-uniformity index. "
        f"Resolution-sensitive; quoted here at the converged dx = {dx_um:.0f} µm "
        "(see limitations).",
        "- **predicted kill** — fraction of tumor cells killed. **A model "
        "output with no validation label.**",
        "",
        "---",
        "",
    ]

    # ---- section 3: attribution ----------------------------------------- #
    sp = check["spearman_total_vs_kill"]
    L += [
        "## 3. Why each drug ranks where it does — computed attribution",
        "",
        "### 3.1 The method (an exact additive decomposition, not a story)",
        "",
        "Predicted kill is driven by how much drug the tumor cells see relative "
        "to how much they need, then discounted by this tumor's own resistance "
        "state. Each drug's standing is decomposed into three **exactly "
        "additive log10 terms** ('decades'), each centred on the panel mean so "
        "`+` reads as advantage over the panel and `−` as penalty:",
        "",
        "```",
        "C_tumor_i = (mean AUC at this tumor's TUMOR-cell positions) / T_window   [µM]",
        "drive_i   = C_tumor_i / IC50_i                                           [-]",
        "",
        "delivery_i = log10 C_tumor_i - mean_j log10 C_tumor_j",
        "potency_i  = -( log10 IC50_i - mean_j log10 IC50_j )",
        "o2resist_i = log10 rho_i     - mean_j log10 rho_j",
        "",
        "  rho_i = kill_i / kill_i(resistance OFF)   <- measured by a SECOND simulation",
        "                                              of the same drug on the same",
        "                                              tissue with hypoxia_resistance=1",
        "                                              and quiescent_resistance=1",
        "",
        "identity:  delivery_i + potency_i  ==  log10 drive_i - mean_j log10 drive_j",
        "```",
        "",
        "The **dominant driver** for a drug is the term of largest magnitude "
        "(labelled *advantage* or *penalty* by its sign). The reason for each "
        "**adjacent rank gap** is the term that differs most between the two "
        "drugs. Nothing here is asserted — the O₂-resistance term in particular "
        "comes from a real counterfactual simulation, which is why this "
        "experiment runs each drug twice. `rho_i` bundles the hypoxic *and* the "
        "quiescent kill penalty (both are switched off together), so it is "
        "named for the oxygen-driven resistance state as a whole, not hypoxia "
        "alone.",
        "",
        "### 3.2 Attribution table (decades vs the panel mean)",
        "",
        "| # | drug | delivery | potency | O₂-resist. | net | dominant "
        "driver | exposure pct | potency pct | O₂-resistance cost |",
        "|---:|---|---:|---:|---:|---:|:--|---:|---:|---:|",
    ]
    for i, r in enumerate(rows, start=1):
        L.append(
            f"| {i} | **{r['name']}** | {r['delivery_decades']:+.2f} | "
            f"{r['potency_decades']:+.2f} | {r['hypoxia_decades']:+.2f} | "
            f"{r['total_decades']:+.2f} | {r['dominant']} | "
            f"{r['exposure_pct']:.0f} | {r['potency_pct']:.0f} | "
            f"{r['hypoxia_cost'] * 100:.1f}% |")

    faith = ("faithful — the decomposition reproduces the simulator's own "
             "ranking" if sp >= 0.9 else
             "partially faithful — the decomposition captures the trend but not "
             "the exact ordering" if sp >= 0.5 else
             "NOT faithful — report this as a failure of the attribution, not of "
             "the simulation")
    L += [
        "",
        "### 3.3 Is the attribution itself faithful? (falsifiability check)",
        "",
        "An attribution that does not reproduce the thing it claims to explain "
        "is decoration. Spearman rank correlation, across the panel, between the "
        "**net attribution** and the **simulated kill**:",
        "",
        f"> **ρ(net attribution, predicted kill) = {sp:+.3f}** → {faith}.",
        "",
        f"Component-wise: ρ(delivery, kill) = "
        f"{check['spearman_delivery_vs_kill']:+.3f}, ρ(potency, kill) = "
        f"{check['spearman_potency_vs_kill']:+.3f}. With only "
        f"{len(rows)} drugs a single ρ is a weak statistic — it is reported to "
        "be checkable, not because it is powerful.",
        "",
        "### 3.4 Per-drug reading",
        "",
    ]
    for i, r in enumerate(rows, start=1):
        L += _drug_narrative(i, r, rows, hypoxic_frac, mean_ecm)

    # ---- section 4: figure + what would falsify -------------------------- #
    L += [
        "---",
        "",
        "## 4. Decision panel figure",
        "",
        f"![virtual drug panel]({png_name})",
        "",
        f"`results/{png_name}` — row 1: the measured tumor (cell types, "
        "vascular density, ECM density, oxygen with the hypoxia threshold "
        "contoured, niche regions). Row 2: the predicted exposure (AUC) map for "
        "each drug on that same tissue, in ranked order on one shared colour "
        "scale, with the top-ranked drug ringed. Row 3: the decision plane "
        "(penetration length vs IC50, marker area ∝ predicted kill) and the "
        "ranked kill bars. Row 4: the attribution decomposition above.",
        "",
        "---",
        "",
        "## 5. What would have to be true for this to be trustworthy",
        "",
        "State this plainly, because the report is otherwise easy to over-read:",
        "",
        "1. **A tissue with a paired, measured drug response.** Everything here "
        "is unfalsifiable until a tumor with a known outcome is run through the "
        "identical pipeline and the predicted ranking is scored against the "
        "observed one.",
        "2. **A potency measurement in context.** Cell-line IC50 must be "
        "replaced by (or at minimum reconciled with) a 3-D / organoid / "
        "patient-derived measurement from tissue resembling this one.",
        "3. **More than one penetration anchor.** α currently rides on a single "
        "doxorubicin measurement; independent per-drug perivascular-distance "
        "data would test whether the ordering transfers outside the "
        "strongly-bound class.",
        "4. **A measured per-cell sensitivity**, replacing today's explicit "
        "assumption (uniform 1.0, or a panel-determinant proxy where one exists).",
        "",
        "Until then the correct reading of this artifact is: *this is what "
        "validated transport physics and measured ex-vivo potency predict on "
        "this specific measured geometry* — a hypothesis generator, not a "
        "decision.",
        "",
        "---",
        "",
        lim,
        "",
        f"*Generated {time.strftime('%Y-%m-%d %H:%M:%S')} by "
        f"`experiments/patient_report.py` in {runtime_s:.0f} s. "
        f"Deterministic (seed 0) — re-running reproduces this file.*",
        "",
    ]
    return "\n".join(L)


# --------------------------------------------------------------------------- #
# 5. driver
# --------------------------------------------------------------------------- #
def main(source: str = "xenium", sample: Optional[str] = None,
         dx_um: float = DX_UM, drugs: Sequence[str] = PANEL) -> Dict:
    t0 = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)
    key = source.strip().lower()
    tissue_label = _TISSUE_LABEL.get(key, key)

    print(f"[1/5] reconstructing the REAL {key} microenvironment "
          f"(once, dx={dx_um} µm) ...")
    try:
        micro = reconstruct_real(source=key, sample=sample, dx_um=float(dx_um))
    except (OSError, ConnectionError, TimeoutError) as exc:  # retry once
        print(f"      transient failure ({exc}); retrying once ...")
        time.sleep(5)
        micro = reconstruct_real(source=key, sample=sample, dx_um=float(dx_um))
    report = reconstruction_report(micro)
    mean_ecm = float(report["mean_ecm_density"])
    print(f"      {report['n_cells']:,} real cells (tumor {report['n_tumor']:,}) "
          f"grid {report['grid']} dx={report['dx_um']:.0f} µm | "
          f"vasc {report['mean_vascular_density']:.3f} ecm {mean_ecm:.3f} "
          f"o2 {report['mean_oxygen']:.3f}")

    print("[2/5] labelling niche regions against the Ivy GAP atlas ...")
    regions = niche_regions(micro, key, sample)
    if regions:
        print("      " + "  ".join(
            f"{r}={regions['tumor_fractions'][r] * 100:.1f}%"
            for r in regions["regions_present"]) + "  (% of tumor cells)")
    profile = microenvironment_profile(micro, report, regions)

    alpha = calibrate_alpha()
    print(f"[3/5] perivascular-sink calibration α = {alpha:.1f} "
          "(fit ONCE on doxorubicin, not refit here)")

    print(f"[4/5] simulating {len(drugs)} drugs × 2 runs "
          "(real + resistance-off counterfactual) on the SAME tumor ...")
    rows = simulate_panel(micro, alpha, mean_ecm, drugs=drugs, dx_um=dx_um)
    hypoxic_frac = float(rows[0]["hypoxic"])  # identical tissue for every drug

    check = attribute(rows, window_s=WINDOW_S)
    rows.sort(key=lambda r: r["kill"], reverse=True)

    print("[5/5] rendering the decision panel + report ...")
    png_name = f"patient_report_{key}.png"
    md_name = f"patient_report_{key}.md"
    out_png = os.path.join(_RESULTS_DIR, png_name)
    out_md = os.path.join(_RESULTS_DIR, md_name)

    panel = {
        "tissue_name": key,
        "tissue_label": tissue_label,
        "report": report,
        "extent_um": report["extent_um"],
        "dx_um": float(report["dx_um"]),
        "window_h": WINDOW_S / 3600.0,
        "tissue_mask": np.asarray(micro.domain.tissue_mask, dtype=bool),
        "cells": {"x": np.asarray(micro.cells.x),
                  "y": np.asarray(micro.cells.y),
                  "cell_type": np.asarray(micro.cells.cell_type)},
        "fields": {"vascular": micro.fields.vascular_density,
                   "ecm": micro.fields.ecm_density,
                   "oxygen": rows[0]["oxygen_field"]},
        "hypoxic_fraction": hypoxic_frac,
        "hypoxia_threshold": OxygenConfig().hypoxia_threshold,
        "region_fields": regions.get("region_fields", {}),
        "region_tumor_fractions": regions.get("tumor_fractions", {}),
        "rows": rows,
        "recommended": rows[0]["name"],
    }
    plot_patient_report(panel, out_png)
    print(f"      figure written: {out_png}")

    runtime = time.time() - t0
    md = build_markdown(
        tissue_key=key, tissue_label=tissue_label, report=report,
        profile=profile, rows=rows, check=check, alpha=alpha,
        hypoxic_frac=hypoxic_frac, mean_ecm=mean_ecm, dx_um=float(report["dx_um"]),
        png_name=png_name, runtime_s=runtime)
    with open(out_md, "w") as fh:
        fh.write(md)
    print(f"      report written: {out_md}")

    _print_summary(key, report, rows, check, hypoxic_frac,
                   profile.get('states', {}), runtime)
    return {"rows": rows, "check": check, "png": out_png, "md": out_md}


def _print_summary(key, report, rows, check, hypoxic_frac, states, runtime) -> None:
    W = 128
    print("\n" + "=" * W)
    print(f"VIRTUAL DRUG PANEL — real {key} tumor "
          f"({report['n_cells']:,} cells, dx={report['dx_um']:.0f} µm) "
          "— PREDICTIONS ONLY, NOT CLINICALLY VALIDATED")
    print(f"sensitivity mode = '{SENSITIVITY_MODE}' (identical for every drug "
          "-> fair cross-drug ranking; absolute kills are upper bounds)  |  "
          f"tissue mask = {report.get('tissue_mask_mode')} "
          f"({(report.get('tissue_coverage_fraction') or 0.0)*100:.1f}% coverage, "
          "dx-independent)")
    print("=" * W)
    hdr = (f"{'#':>2} {'drug':14s} {'IC50(uM)':>9} {'src':>10} {'bbb':>6} "
           f"{'D_eff':>7} {'L(um)':>8} {'tumorAUC':>10} {'pen':>6} {'kill':>8} "
           f"  {'dominant driver':29}")
    print(hdr)
    print("-" * len(hdr))
    for i, r in enumerate(rows, start=1):
        src = "GDSC" if r["ic50_confidence"] == "measured" else r["ic50_source"]
        print(f"{i:>2} {r['name']:14s} {r['ic50']:9.3g} {src:>10} {r['bbb']:6.3f} "
              f"{r['d_eff']:7.0f} {r['L_um']:8.1f} {r['tumor_auc']:10.3g} "
              f"{r['pen']:6.3f} {r['kill']:8.4f}   {r['dominant']:29}")
    print("-" * len(hdr))
    print("attribution (decades vs panel mean):")
    for i, r in enumerate(rows, start=1):
        print(f"  #{i} {r['name']:14s} delivery {r['delivery_decades']:+6.2f}  "
              f"potency {r['potency_decades']:+6.2f}  "
              f"o2-resist {r['hypoxia_decades']:+6.2f}  "
              f"net {r['total_decades']:+6.2f}  "
              f"(O2 resistance costs {r['hypoxia_cost'] * 100:5.1f}% of its kill)")
    print(f"attribution faithfulness: rho(net, kill) = "
          f"{check['spearman_total_vs_kill']:+.3f}")
    st = states or {}
    print(f"tumor resistance state  : hypoxic {hypoxic_frac * 100:.1f}%, "
          f"quiescent {st.get('quiescent', 0.0) * 100:.1f}%, "
          f"viable {st.get('viable', 0.0) * 100:.1f}% "
          "(suppresses kill for EVERY drug; hypoxic fraction NOT grid-converged)")
    print("=" * W)
    print("LIMITATIONS: research prototype. Predictions are NOT clinically "
          "validated; no outcome-labelled\ndata exists for this tissue; potency "
          "is measured in dissociated cells; penetration_index is\n"
          "resolution-sensitive (see results/numerical_rigor.md). NOT FOR "
          "CLINICAL USE.")
    print("=" * W)
    print(f"[runtime {runtime:.1f} s]")


def _cli() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--source", default="xenium", choices=("xenium", "visium"),
                   help="real spatial-transcriptomics source (default: xenium)")
    p.add_argument("--sample", default=None,
                   help="sample id / cache path (default: the loader default)")
    p.add_argument("--dx", type=float, default=DX_UM,
                   help=f"grid spacing in um (default {DX_UM:g}, the converged "
                        "resolution for penetration_index)")
    a = p.parse_args()
    main(source=a.source, sample=a.sample, dx_um=a.dx)


if __name__ == "__main__":
    _cli()
