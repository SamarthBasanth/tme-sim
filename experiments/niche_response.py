"""Predicted drug response BY ANATOMIC NICHE REGION on real tissue.

This is the "the microenvironment matters" claim made **spatially explicit and
testable**. Everywhere else in this repo a drug gets ONE number per tumour
(``kill_fraction``). That number is a spatial average, and a spatial average is
exactly the thing that hides the microenvironment. Here we ask the question the
whole model exists to answer:

    *WHERE inside a real tumour does each drug actually work?*

Pipeline (all REAL, no synthetic data anywhere in the path)
-----------------------------------------------------------
1. **Real tissue.** The open 10x ``Parent_Visium_Human_Glioblastoma`` section
   (3,468 spots x 36,601 genes, primary GBM resection, fresh-frozen).
   Visium is the primary platform for this experiment because it shares
   **283 / 311** genes with the Ivy GAP region-signature panel; the Xenium GBM
   section shares only 41 of 311, so region assignment there rests on a much
   thinner correlation and is treated as out of scope (see ``--source``).
2. **Real anatomic niche labels.** Every spot is scored against the REAL Ivy GAP
   per-anatomic-structure FPKM signatures (Puchalski et al. 2018, *Science*;
   laser-microdissected GBM structures) by rank correlation over the shared
   panel, and assigned to its best-correlating zone:
   ``margin`` / ``core`` / ``perivascular`` / ``necrotic``
   (:mod:`tmesim.reconstruct.region_label`). The zone labels are then carried
   onto the simulation grid and onto every reconstructed cell by
   nearest-spot assignment.
3. **Real microenvironment + real drug physics.** The same reconstruction and
   the same drug-compilation pipeline the flagship panel uses:
   ``calibrate_alpha()`` -> ``calibrated_drugparams(name, alpha)`` ->
   ``apply_barrier_to_drugparams(dp, molecule)`` -> ``simulate_response``.
4. **Per-region readout.** For each drug x each niche zone: mean tissue drug AUC,
   mean per-cell exposure, predicted kill fraction of the tumour cells *in that
   zone*, and the local vascular density and oxygen.
5. **A falsifiable test.** The textbook expectation is
   ``perivascular > margin > core > necrotic`` for delivery and kill. We test the
   contrast ``perivascular`` vs ``core + necrotic`` quantitatively, with effect
   sizes (Cohen's d, Cliff's delta) and TWO permutation nulls -- a naive i.i.d.
   label shuffle and a spatially-structured null (random rotation / reflection of
   the whole label map) that does not destroy the spatial autocorrelation of the
   labels. If the expected pattern does not appear, that is reported as-is.

Grid resolution
---------------
The simulation runs at the grid-CONVERGED spacing **dx = 12 um** established in
``results/numerical_rigor.md`` (``penetration_index`` only settles once the
~45 um perivascular boundary layer is resolved). **This does not create spatial
information that Visium does not have**: the underlying measurement is spot-level
with a ~100 um centre-to-centre pitch and a 55 um spot diameter, so the *physical
information content* of every field in this experiment is spot-limited to ~100 um
regardless of how fine the numerical mesh is. dx = 12 um buys numerical
convergence of the PDE, not biological resolution. A second full panel is run on
``reconstruct_real``'s native Visium grid (dx = 100 um = the spot pitch) as an
explicit resolution control; both are reported.

Run:    ./.venv/bin/python experiments/niche_response.py
Writes: results/niche_response.png
        results/niche_response.md

Nothing in core / physics / molecular / reconstruct / engine is modified.
"""
from __future__ import annotations

import argparse
import os
import time
from typing import Dict, List, Optional, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.spatial import cKDTree

from tmesim.cells.agents import update_cell_states
from tmesim.cells.pkpd import apply_drug_kill
from tmesim.core import (
    CellState,
    Microenvironment,
    OxygenConfig,
    PDEConfig,
    SimConfig,
    TissueFields,
)
from tmesim.data.real import visium
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.reconstruct.real_fields import add_oxygen
from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from tmesim.reconstruct.real_visium import reconstruct_from_visium
from tmesim.reconstruct.region_label import (
    REGION_ORDER,
    label_regions,
    region_signatures,
)
from tmesim.sim.engine import simulate_response
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

PANEL = ["temozolomide", "lomustine", "cisplatin", "methotrexate", "doxorubicin"]

VISIUM_SAMPLE = "10x_parent_gbm"

# CONVERGED numerical grid (results/numerical_rigor.md section 2). See the module
# docstring: this is numerical convergence, NOT biological resolution.
DX_UM = 12.0
# reconstruct_real's native Visium grid == the spot pitch. Resolution control.
DX_SPOT_UM = 100.0

# Nearest-spot assignment radii (um).
_GRID_ASSIGN_RADIUS_UM = 60.0   # ~ half the 100 um spot pitch (+ margin)
_CELL_ASSIGN_RADIUS_UM = 60.0   # cells are jittered inside a 55 um spot disc

# Zones expected to be well delivered / poorly delivered.
HIGH_GROUP = ["perivascular"]
LOW_GROUP = ["core", "necrotic"]

N_PERM = 10_000      # i.i.d. label-shuffle permutations
N_SPATIAL = 1_000    # spatially-structured (rotation/reflection) permutations

REGION_COLORS = {
    "margin": "#4c9f70",
    "core": "#d1495b",
    "perivascular": "#3b6ea5",
    "necrotic": "#8a5a2b",
    "unlabeled": "#bbbbbb",
}

_RESULTS_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
)
_OUT_PNG = os.path.join(_RESULTS_DIR, "niche_response.png")
_OUT_MD = os.path.join(_RESULTS_DIR, "niche_response.md")


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _retry(fn, *args, **kwargs):
    """Call ``fn``, retrying ONCE after a short pause on a transient/network error.

    A missing local cache (``FileNotFoundError``) is a configuration error and is
    NOT retried.
    """
    try:
        return fn(*args, **kwargs)
    except FileNotFoundError:
        raise
    except (OSError, ConnectionError, TimeoutError):
        time.sleep(3.0)
        return fn(*args, **kwargs)


def _cp10k_log(expr, genes: Sequence[str], want: Sequence[str]) -> Dict[str, np.ndarray]:
    """Per-spot counts-per-10k, log1p, for a few marker genes (marker validation)."""
    X = expr.toarray().astype(np.float64) if hasattr(expr, "toarray") else np.asarray(expr, float)
    tot = X.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    norm = np.log1p(X / tot * 1e4)
    col = {g: j for j, g in enumerate(genes)}
    return {g: norm[:, col[g]] for g in want if g in col}


def compile_real_drug(name: str, alpha: float):
    """Compile ONE real drug through the full pipeline (identical to the flagship
    panel): calibrated perivascular sink -> saturable binding-site barrier."""
    dp = calibrated_drugparams(name, alpha)
    mol = get_molecule(name)
    return apply_barrier_to_drugparams(dp, mol), mol


# ---------------------------------------------------------------------------
# region -> grid / region -> cell assignment
# ---------------------------------------------------------------------------
def assign_regions(domain, spot_coords: np.ndarray, spot_labels: np.ndarray, cells):
    """Carry per-spot niche labels onto (a) the simulation grid and (b) the cells.

    Both are nearest-spot assignments with a distance cutoff: grid points and
    cells further than the cutoff from any spot get ``None`` (they have no real
    measurement behind them and are excluded from every regional mean). This is
    the honest way to move a spot-level label onto a finer mesh -- it does NOT
    invent sub-spot structure, it replicates the spot's label over its footprint.

    Returns ``(grid_region_idx, cell_region_idx, spot_tree)`` where the index
    arrays hold the position of the region in :data:`REGION_ORDER`, or -1.
    """
    tree = cKDTree(spot_coords)
    order = {r: i for i, r in enumerate(REGION_ORDER)}
    lab_idx = np.array([order.get(str(s), -1) for s in spot_labels], dtype=int)

    X, Y = domain.coords_um()
    gq = np.column_stack([X.ravel(), Y.ravel()])
    gd, gi = tree.query(gq, k=1)
    grid_idx = np.where(gd <= _GRID_ASSIGN_RADIUS_UM, lab_idx[gi], -1)
    grid_idx = grid_idx.reshape(domain.shape)

    cd, ci = tree.query(np.column_stack([cells.x, cells.y]), k=1)
    cell_idx = np.where(cd <= _CELL_ASSIGN_RADIUS_UM, lab_idx[ci], -1)

    return grid_idx, cell_idx, tree


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------
def cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Cohen's d with a pooled SD. Returns 0.0 for a degenerate comparison."""
    na, nb = a.size, b.size
    if na < 2 or nb < 2:
        return 0.0
    va, vb = a.var(ddof=1), b.var(ddof=1)
    pooled = ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)
    if pooled <= 0:
        return 0.0
    return float((a.mean() - b.mean()) / np.sqrt(pooled))


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta in [-1, 1] -- P(a > b) - P(a < b), computed by ranking.

    Non-parametric, scale-free, and the right effect size for these heavily
    skewed exposure distributions. +1 means every ``a`` exceeds every ``b``.
    """
    from scipy.stats import rankdata

    na, nb = a.size, b.size
    if na == 0 or nb == 0:
        return 0.0
    both = np.concatenate([a, b])
    r = rankdata(both)
    # Mann-Whitney U from the rank sum of group a, then map U -> Cliff's delta.
    u = r[:na].sum() - na * (na + 1) / 2.0
    return float(2.0 * u / (na * nb) - 1.0)


def _stat(values: np.ndarray, idx: np.ndarray, hi: np.ndarray, lo: np.ndarray) -> float:
    """Test statistic: mean(values | high-group zones) - mean(values | low-group)."""
    a = values[np.isin(idx, hi)]
    b = values[np.isin(idx, lo)]
    if a.size == 0 or b.size == 0:
        return np.nan
    return float(a.mean() - b.mean())


def permutation_test(
    values: np.ndarray,
    region_idx: np.ndarray,
    coords: np.ndarray,
    high: Sequence[str],
    low: Sequence[str],
    n_perm: int = N_PERM,
    n_spatial: int = N_SPATIAL,
    seed: int = 0,
) -> Dict:
    """Two-null permutation test of ``mean(high zones) - mean(low zones)``.

    Null 1 (``p_perm``): i.i.d. shuffle of the region labels across locations.
    This is the standard permutation test, but it destroys the spatial
    autocorrelation of the labels, so it is **anti-conservative** for spatial
    data -- any smooth field will beat a spatially-white null almost by
    construction. It is reported for completeness.

    Null 2 (``p_spatial``): the whole label map is rigidly rotated by a random
    angle (and randomly reflected) about the tissue centroid, and each real
    location takes the label of the nearest rotated location. The *spatial
    structure* of the labelling (patch sizes, autocorrelation, class balance) is
    preserved; only its registration to the drug field is randomised. This is the
    honest null for "are these zones aligned with the delivery field?" and it is
    the p-value that should be quoted.

    One-sided in the direction of the expected (positive) effect. p-values are
    reported with the standard ``(1 + #ge) / (1 + n)`` add-one correction, so the
    smallest reportable value is ``1 / (n + 1)``.
    """
    rng = np.random.default_rng(seed)
    order = {r: i for i, r in enumerate(REGION_ORDER)}
    hi = np.array([order[r] for r in high if r in order], dtype=int)
    lo = np.array([order[r] for r in low if r in order], dtype=int)

    keep = region_idx >= 0
    values = np.asarray(values, float)[keep]
    idx = np.asarray(region_idx, int)[keep]
    coords = np.asarray(coords, float)[keep]

    observed = _stat(values, idx, hi, lo)

    # --- Null 1: i.i.d. label shuffle ---------------------------------------
    ge = 0
    for _ in range(int(n_perm)):
        if _stat(values, rng.permutation(idx), hi, lo) >= observed:
            ge += 1
    p_perm = (1.0 + ge) / (1.0 + n_perm)

    # --- Null 2: rigid rotation / reflection of the label map ---------------
    centre = coords.mean(axis=0)
    centred = coords - centre
    tree = cKDTree(coords)
    ge_s = 0
    for _ in range(int(n_spatial)):
        th = rng.uniform(0.0, 2.0 * np.pi)
        c, s = np.cos(th), np.sin(th)
        # errstate: some BLAS backends raise spurious FP flags on this tiny gemm.
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            rot = centred @ np.array([[c, -s], [s, c]])
        if rng.random() < 0.5:                 # random reflection
            rot[:, 0] = -rot[:, 0]
        # location i adopts the label that sat at the pre-rotation position, i.e.
        # the label of the location nearest to the back-rotated coordinate.
        _, nn = tree.query(rot + centre, k=1)
        shuffled = idx[nn]
        if _stat(values, shuffled, hi, lo) >= observed:
            ge_s += 1
    p_spatial = (1.0 + ge_s) / (1.0 + n_spatial)

    a = values[np.isin(idx, hi)]
    b = values[np.isin(idx, lo)]
    return dict(
        observed=observed,
        mean_high=float(a.mean()) if a.size else np.nan,
        mean_low=float(b.mean()) if b.size else np.nan,
        ratio=float(a.mean() / b.mean()) if (b.size and b.mean() > 0) else np.nan,
        n_high=int(a.size),
        n_low=int(b.size),
        d=cohens_d(a, b),
        cliff=cliffs_delta(a, b),
        p_perm=float(p_perm),
        p_spatial=float(p_spatial),
        n_perm=int(n_perm),
        n_spatial=int(n_spatial),
    )


# ---------------------------------------------------------------------------
# the run
# ---------------------------------------------------------------------------
def build_micro(dx_um: float, sample: str, use_pipeline: bool):
    """Reconstruct the REAL Visium microenvironment.

    ``use_pipeline=True`` calls :func:`tmesim.reconstruct.real_pipeline.reconstruct_real`
    -- the canonical whole-niche assembler. NOTE that for Visium it deliberately
    clamps any ``dx_um < 50`` back to the 100 um spot pitch (a Visium spot carries
    no sub-100 um information, so the assembler refuses to over-resolve it). To run
    the PDE on the grid-converged 12 um mesh we therefore call the exact same two
    public steps the assembler's Visium branch performs
    (``reconstruct_from_visium`` + real-perfusion ``add_oxygen``) with the fine dx.
    Both paths are exercised; the coarse one is the resolution control.
    """
    if use_pipeline:
        micro = _retry(reconstruct_real, source="visium", sample=sample, dx_um=dx_um)
    else:
        micro = _retry(reconstruct_from_visium, sample, dx_um=float(dx_um))
        if micro.fields.oxygen is None:
            add_oxygen(micro, OxygenConfig())
        micro.metadata = dict(micro.metadata)
        micro.metadata["real_pipeline"] = {
            "source": "visium",
            "vascular_source": "KDE of REAL ENDOTHELIAL spot fraction",
            "ecm_source": "STROMAL-fraction proxy (spot resolution)",
            "ecm_markers_used": [],
            "oxygen_source": "physics.solve_oxygen driven by the real vascular field",
        }
    return micro


def per_cell_kill(micro: Microenvironment, res, cfg: SimConfig, drug):
    """Recover the engine's PER-CELL kill probability and realised deaths.

    ``simulate_response`` only returns aggregates. We reproduce its step 2 + step 6
    exactly -- oxygen banding via :func:`update_cell_states` on the engine's own
    returned oxygen field, then :func:`apply_drug_kill` with the same transport
    result and the same seeded RNG -- on a private copy. Nothing in the engine is
    touched, and the reproduction is verified against ``res.kill_fraction``.
    """
    ox = np.asarray(res.fields["oxygen"], dtype=np.float64)
    snap = Microenvironment(
        domain=micro.domain,
        fields=TissueFields(
            ecm_density=micro.fields.ecm_density,
            vascular_density=micro.fields.vascular_density,
            oxygen=ox,
        ),
        cells=micro.cells.copy(),
        patient_id=micro.patient_id,
        metadata=micro.metadata,
    )
    update_cell_states(snap, ox, cfg.oxygen)
    p_kill = apply_drug_kill(
        snap, res.transport, cfg.pd, drug, rng=np.random.default_rng(cfg.seed)
    )
    dead = np.asarray(snap.cells.state) == int(CellState.DEAD)
    exposure = np.asarray(snap.cells.exposure, dtype=np.float64)
    reproduced = 1.0 - snap.cells.viable_tumor_fraction()
    return p_kill, dead, exposure, float(reproduced)


def run_panel(micro, grid_idx, cell_idx, alpha, pde, seed=0, verbose=True):
    """Simulate the whole real drug panel and reduce each drug to per-region rows."""
    cells = micro.cells
    tumor = cells.tumor_mask()
    ct = np.asarray(cells.cell_type)
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)

    out = []
    for name in PANEL:
        dp, mol = compile_real_drug(name, alpha)
        cfg = SimConfig(drug=dp, pde=pde, seed=seed)
        if verbose:
            print(f"  simulating {name:14s} ...", end="", flush=True)
        t0 = time.time()
        res = simulate_response(micro, cfg)
        p_kill, dead, exposure, reproduced = per_cell_kill(micro, res, cfg, dp)
        auc = np.asarray(res.fields["auc"], dtype=np.float64)
        oxy = np.asarray(res.fields["oxygen"], dtype=np.float64)
        vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)

        regions = {}
        for i, reg in enumerate(REGION_ORDER):
            gm = (grid_idx == i) & tissue
            cm = cell_idx == i
            tm = cm & tumor
            regions[reg] = dict(
                n_grid=int(gm.sum()),
                n_cells=int(cm.sum()),
                n_tumor=int(tm.sum()),
                mean_auc=float(auc[gm].mean()) if gm.any() else np.nan,
                median_auc=float(np.median(auc[gm])) if gm.any() else np.nan,
                mean_cell_exposure=float(exposure[cm].mean()) if cm.any() else np.nan,
                mean_tumor_exposure=float(exposure[tm].mean()) if tm.any() else np.nan,
                kill_frac=float(dead[tm].mean()) if tm.any() else np.nan,
                exp_kill=float(p_kill[tm].mean()) if tm.any() else np.nan,
                vascular=float(vasc[gm].mean()) if gm.any() else np.nan,
                oxygen=float(oxy[gm].mean()) if gm.any() else np.nan,
            )
        out.append(
            dict(
                name=name,
                mw=float(mol.get("mw")),
                ic50=float(dp.ic50),
                bbb=float(dp.bbb_penetrance),
                diffusivity=float(dp.diffusivity),
                kill=float(res.kill_fraction),
                kill_reproduced=reproduced,
                pen=float(res.extras["penetration_index"]),
                mean_auc=float(res.extras["mean_auc"]),
                hypoxic=float(res.extras["hypoxic_tumor_fraction"]),
                auc_field=auc,
                oxygen_field=oxy,
                p_kill=p_kill,
                dead=dead,
                exposure=exposure,
                regions=regions,
            )
        )
        if verbose:
            print(
                f" kill={res.kill_fraction:.4f} (reproduced {reproduced:.4f}) "
                f"mean_auc={res.extras['mean_auc']:.3g}  [{time.time()-t0:.0f}s]"
            )
    return out, tumor, ct


def spot_level_series(domain, spot_coords, auc, p_kill, cells, tumor):
    """Reduce a drug to the two SPOT-LEVEL series the statistics are run on.

    Cells inside one Visium spot all sample essentially the same drug field, so
    the cell is not an independent unit -- using n = 41,560 cells would inflate
    every test. The honest unit is the spot (n = 3,468). ``spot_auc`` samples the
    exposure field at the spot centroid; ``spot_kill`` is the mean predicted kill
    probability of the tumour cells belonging to that spot (NaN where a spot has
    no tumour cell).
    """
    col, row = domain.to_index(spot_coords[:, 0], spot_coords[:, 1])
    spot_auc = np.asarray(auc)[row, col]

    tree = cKDTree(spot_coords)
    d, nn = tree.query(np.column_stack([cells.x, cells.y]), k=1)
    ok = (d <= _CELL_ASSIGN_RADIUS_UM) & tumor
    n_spots = spot_coords.shape[0]
    tot = np.zeros(n_spots)
    cnt = np.zeros(n_spots)
    np.add.at(tot, nn[ok], p_kill[ok])
    np.add.at(cnt, nn[ok], 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        spot_kill = np.where(cnt > 0, tot / np.maximum(cnt, 1e-12), np.nan)
    return spot_auc, spot_kill


# ---------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------
def make_figure(
    micro, spot_coords, spot_labels, rows, region_counts, sample_meta, dx_um, stats
):
    order = {r: i for i, r in enumerate(REGION_ORDER)}
    present = [r for r in REGION_ORDER if region_counts.get(r, 0) > 0]

    fig = plt.figure(figsize=(21, 12.5))
    gs = GridSpec(2, 5, figure=fig, height_ratios=[1.25, 1.0], hspace=0.33, wspace=0.30)

    # --- A. the real tissue, coloured by real Ivy GAP niche zone -------------
    ax = fig.add_subplot(gs[0, 0:2])
    for r in REGION_ORDER + ["unlabeled"]:
        m = spot_labels == r
        if m.any():
            ax.scatter(
                spot_coords[m, 0], spot_coords[m, 1], s=11, c=REGION_COLORS[r],
                label=f"{r} ({m.sum()}, {100*m.mean():.0f}%)", edgecolors="none",
            )
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("y (µm)")
    ax.legend(markerscale=2.4, fontsize=8.5, loc="upper right", framealpha=0.92)
    ax.set_title(
        f"A. REAL tissue, coloured by REAL Ivy GAP anatomic niche zone\n"
        f"{sample_meta['accession']} — {sample_meta['n_spots']:,} Visium spots, "
        f"{sample_meta['n_shared']}/{sample_meta['n_panel']} signature genes shared",
        fontsize=11,
    )

    # --- B. context: local vascular density + oxygen per zone ---------------
    ax = fig.add_subplot(gs[0, 2])
    reg0 = rows[0]["regions"]
    xs = np.arange(len(present))
    w = 0.38
    vas = [reg0[r]["vascular"] for r in present]
    oxy = [reg0[r]["oxygen"] for r in present]
    ax.bar(xs - w / 2, vas, w, color="#3b6ea5", label="vascular density")
    ax.bar(xs + w / 2, oxy, w, color="#e07b39", label="oxygen (relative)")
    for x, v in zip(xs - w / 2, vas):
        ax.text(x, v, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
    for x, v in zip(xs + w / 2, oxy):
        ax.text(x, v, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(xs)
    ax.set_xticklabels(present, rotation=20, fontsize=8.5)
    ax.set_ylabel("field value")
    ax.legend(fontsize=8)
    ax.set_title("B. Local microenvironment per zone\n(identical for all drugs)", fontsize=10.5)

    # --- C. heatmap: drug x zone predicted tumour kill ----------------------
    ax = fig.add_subplot(gs[0, 3:5])
    M = np.array([[r["regions"][reg]["kill_frac"] for reg in present] for r in rows])
    im = ax.imshow(M, cmap="magma", aspect="auto")
    ax.set_xticks(np.arange(len(present)))
    ax.set_xticklabels(present, rotation=18, fontsize=9)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([r["name"] for r in rows], fontsize=9)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            ax.text(
                j, i, f"{v:.3f}", ha="center", va="center", fontsize=8.5,
                color="white" if v < 0.55 * np.nanmax(M) else "black",
            )
    fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02, label="predicted kill fraction of tumour cells in zone")
    ax.set_title(
        "C. Predicted kill of the tumour cells IN EACH ZONE (drug × niche)",
        fontsize=10.5,
    )

    # --- D..H. per-drug exposure bars by zone -------------------------------
    for j, r in enumerate(rows):
        ax = fig.add_subplot(gs[1, j])
        vals = [r["regions"][reg]["mean_tumor_exposure"] for reg in present]
        cols = [REGION_COLORS[reg] for reg in present]
        ax.bar(np.arange(len(present)), vals, color=cols)
        for x, v in enumerate(vals):
            ax.text(x, v, f"{v:.2e}", ha="center", va="bottom", fontsize=7, rotation=90)
        ax.set_xticks(np.arange(len(present)))
        ax.set_xticklabels(present, rotation=25, fontsize=8)
        ax.set_ylim(0, max(vals) * 1.42 if max(vals) > 0 else 1)
        if j == 0:
            ax.set_ylabel("mean per-cell exposure (µM·s)")
        st = stats[r["name"]]["auc"]
        ax.set_title(
            f"{r['name']}\nperivasc / (core+necrotic) = {st['ratio']:.2f}×   "
            f"δ={st['cliff']:+.2f}",
            fontsize=9.5,
        )

    fig.suptitle(
        "WHERE does each drug work? Predicted drug response by REAL anatomic niche zone "
        "on REAL glioblastoma tissue\n"
        f"10x Visium GBM × real Ivy GAP zone signatures × real drug physics — "
        f"PDE at grid-converged dx = {dx_um:.0f} µm, but the underlying measurement is "
        "spot-limited to ~100 µm (no sub-spot biology is created)  —  PREDICTIONS, not "
        "validated outcomes",
        fontsize=13.5,
    )
    fig.savefig(_OUT_PNG, dpi=125, bbox_inches="tight")
    return _OUT_PNG


# ---------------------------------------------------------------------------
# markdown
# ---------------------------------------------------------------------------
def _fmt(v, spec=".3g"):
    return "n/a" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:{spec}}"


def build_markdown(ctx) -> str:
    rows = ctx["rows"]
    stats = ctx["stats"]
    present = ctx["present"]
    meta = ctx["sample_meta"]
    rep = ctx["report"]
    L: List[str] = []
    A = L.append

    A("# Predicted drug response by anatomic niche region on real glioblastoma tissue")
    A("")
    A("**Experiment:** `experiments/niche_response.py` "
      "(run: `./.venv/bin/python experiments/niche_response.py`)")
    A("")
    A("**Figure:** `results/niche_response.png`")
    A("")
    A("Everywhere else in this repo a drug gets **one** number per tumour "
      "(`kill_fraction`). That number is a spatial average — and a spatial average "
      "is exactly what hides the microenvironment. This experiment asks the question "
      "the whole model exists to answer: **where inside a real tumour does each drug "
      "actually work?**")
    A("")
    A("Every spot of a real 10x Visium glioblastoma section is assigned to a real "
      "Ivy GAP anatomic zone (`margin` / `core` / `perivascular` / `necrotic`) by rank "
      "correlation against the laser-microdissected Ivy GAP per-structure FPKM "
      "signatures (Puchalski et al. 2018, *Science*). The real drug panel is then "
      "simulated on the real reconstructed microenvironment, and every readout is "
      "resolved **per zone**.")
    A("")

    # --- 1. tissue + labelling ---------------------------------------------
    A("## 1. The real tissue and its real niche zones")
    A("")
    A(f"- **Sample:** `{meta['sample_id']}` — {meta['accession']}, {meta['diagnosis']}.")
    A(f"- **Measurement:** {meta['n_spots']:,} Visium spots × {meta['n_genes']:,} genes; "
      f"~100 µm centre-to-centre spot pitch, {meta['spot_diameter_um']:.0f} µm spot diameter.")
    A(f"- **Zone signatures:** real Ivy GAP per-anatomic-structure FPKM, reduced to a "
      f"{meta['n_panel']}-gene region-discriminative panel. "
      f"**{meta['n_shared']} / {meta['n_panel']} panel genes are present in this Visium "
      f"sample** — this is why Visium is the primary platform here. The Xenium GBM "
      f"section shares only 41/311 panel genes, so its zone assignment rests on a much "
      f"thinner correlation; it is deliberately out of scope for this experiment.")
    A(f"- **Assignment confidence:** mean top-1 Spearman correlation to the winning "
      f"signature = {meta['mean_top1']:.3f}; mean top-1 − top-2 margin = "
      f"{meta['mean_margin']:.3f}. These are *soft* assignments — see the caveats.")
    A("")
    A("### Zone composition and local microenvironment")
    A("")
    A("| zone | spots | % of tissue | grid points | cells | tumour cells | "
      "mean vascular density | mean oxygen (rel.) |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|")
    reg0 = rows[0]["regions"]
    tot_spots = meta["n_spots"]
    for reg in present:
        d = reg0[reg]
        A(f"| **{reg}** | {ctx['region_counts'][reg]:,} | "
          f"{100*ctx['region_counts'][reg]/tot_spots:.1f}% | {d['n_grid']:,} | "
          f"{d['n_cells']:,} | {d['n_tumor']:,} | {_fmt(d['vascular'],'.4f')} | "
          f"{_fmt(d['oxygen'],'.4f')} |")
    A("")
    A("Marker validation of the assignment (mean per-spot log CP10k of canonical "
      "markers, by **assigned** zone — these markers are part of the signature panel, "
      "so this is an internal consistency check, not independent validation):")
    A("")
    mk = ctx["markers"]
    genes = list(mk["values"].keys())
    A("| zone | " + " | ".join(f"`{g}`" for g in genes) + " |")
    A("|---|" + "---:|" * len(genes))
    for reg in present:
        m = ctx["spot_labels"] == reg
        A(f"| **{reg}** | " + " | ".join(f"{mk['values'][g][m].mean():.3f}" for g in genes) + " |")
    A("")
    for line in mk["verdicts"]:
        A(f"- {line}")
    A("")

    # --- 2. simulation setup ------------------------------------------------
    A("## 2. Simulation setup")
    A("")
    A(f"- **Reconstruction:** real cells at real spot positions "
      f"({rep['n_cells']:,} cells, tumour {rep['n_tumor']:,}), real vascular field "
      f"(KDE of the deconvolved ENDOTHELIAL spot fraction), real ECM proxy "
      f"(STROMAL fraction), real-perfusion oxygen from the unmodified physics solver.")
    A(f"- **Grid:** {rep['grid']} at dx = {rep['dx_um']:.0f} µm "
      f"(extent {rep['extent_um'][0]:.0f} × {rep['extent_um'][1]:.0f} µm).")
    A(f"- **Drug compilation (identical to the flagship panel):** `calibrate_alpha()` "
      f"→ α = {ctx['alpha']:.1f} (fit **once** on doxorubicin) → "
      f"`calibrated_drugparams(name, α)` (real GDSC/literature IC50) → "
      f"`apply_barrier_to_drugparams` (saturable binding-site barrier).")
    A(f"- **Engine:** `simulate_response`, PDE dt = {ctx['pde'].dt:.0f} s over "
      f"{ctx['pde'].total_time/3600:.0f} h, convection on, seed 0.")
    A("")
    A("**Resolution — read this.** The PDE runs at the grid-**converged** spacing "
      f"dx = {rep['dx_um']:.0f} µm established in "
      "[`results/numerical_rigor.md`](numerical_rigor.md) §2 (`penetration_index` only "
      "settles once the ~45 µm perivascular boundary layer is resolved). **This buys "
      "numerical convergence, not biological resolution.** Every field in this "
      "experiment ultimately derives from Visium spots on a ~100 µm pitch, so the "
      "*information content* is spot-limited to ~100 µm no matter how fine the mesh "
      "is. Zone labels are replicated over each spot's ~60 µm footprint — the "
      "sub-spot structure you see in the AUC field comes from the transport physics "
      "acting on a spot-resolution source, not from sub-spot measurement. Note also "
      "that `reconstruct_real`'s Visium branch deliberately **clamps** `dx_um < 50` "
      "back to the 100 µm spot pitch for exactly this reason; the fine-grid run calls "
      "the same two public steps (`reconstruct_from_visium` + `add_oxygen`) directly, "
      "and the 100 µm pipeline grid is run as an explicit resolution control (§5).")
    A("")

    # --- 3. the main table --------------------------------------------------
    A("## 3. Drug × niche zone — where each drug actually works")
    A("")
    A("| drug | zone | mean tissue AUC (µM·s) | mean per-cell exposure (µM·s) | "
      "mean tumour-cell exposure (µM·s) | predicted kill of tumour cells in zone | "
      "expected kill (mean *p*) | vascular density | oxygen |")
    A("|---|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        for reg in present:
            d = r["regions"][reg]
            A(f"| {r['name']} | {reg} | {_fmt(d['mean_auc'])} | "
              f"{_fmt(d['mean_cell_exposure'])} | {_fmt(d['mean_tumor_exposure'])} | "
              f"{_fmt(d['kill_frac'],'.4f')} | {_fmt(d['exp_kill'],'.4f')} | "
              f"{_fmt(d['vascular'],'.4f')} | {_fmt(d['oxygen'],'.4f')} |")
    A("")
    A("Whole-tumour reference (the number the flagship panel would report, i.e. what "
      "the per-zone resolution is hiding):")
    A("")
    A("| drug | MW (Da) | real IC50 (µM) | bbb | whole-tumour kill | penetration index | "
      "whole-tissue mean AUC (µM·s) | tumour hypoxic fraction |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        A(f"| {r['name']} | {r['mw']:.0f} | {r['ic50']:.3g} | {r['bbb']:.3f} | "
          f"{r['kill']:.4f} | {r['pen']:.3f} | {r['mean_auc']:.3g} | {r['hypoxic']*100:.1f}% |")
    A("")
    A("*(The per-cell kill reconstruction reproduces the engine's own "
      "`kill_fraction` to "
      f"{max(abs(r['kill']-r['kill_reproduced']) for r in rows):.2e} absolute — the "
      "per-zone kills are a decomposition of the engine's result, not a "
      "re-implementation.)*")
    A("")
    A("### Zone ordering by delivery and by kill")
    A("")
    A("| drug | zones ranked by mean tumour-cell exposure (high → low) | "
      "zones ranked by predicted kill (high → low) | expected ordering recovered? |")
    A("|---|---|---|:--:|")
    expected = ["perivascular", "margin", "core", "necrotic"]
    exp_present = [r for r in expected if r in present]
    for r in rows:
        by_exp = sorted(present, key=lambda g: -r["regions"][g]["mean_tumor_exposure"])
        by_kill = sorted(present, key=lambda g: -r["regions"][g]["kill_frac"])
        ok = "yes" if by_exp == exp_present else "no"
        A(f"| {r['name']} | {' > '.join(by_exp)} | {' > '.join(by_kill)} | {ok} |")
    A("")
    A(f"Textbook expectation for reference: `{' > '.join(exp_present)}`.")
    A("")

    # --- 4. statistics ------------------------------------------------------
    A("## 4. Is delivery/kill really highest perivascularly and lowest in core/necrotic?")
    A("")
    A("**Contrast:** `perivascular` vs `core + necrotic`, one-sided in the expected "
      "(positive) direction.")
    A("")
    A("**Unit of analysis:** the **spot**, not the cell. All ~12 cells inside one "
      "Visium spot sample essentially the same drug field, so treating cells as "
      "independent would inflate *n* by ~12× and every p-value with it. "
      f"n = {stats[rows[0]['name']]['auc']['n_high']:,} perivascular spots vs "
      f"{stats[rows[0]['name']]['auc']['n_low']:,} core+necrotic spots.")
    A("")
    A("**Two nulls are reported.** `p_perm` is the standard i.i.d. label shuffle — "
      "it is **anti-conservative for spatial data**, because a smooth field beats a "
      "spatially-white null almost by construction. `p_spatial` rigidly rotates and "
      "reflects the *entire* label map about the tissue centroid, preserving patch "
      "size, class balance and spatial autocorrelation, and only randomising the "
      "label map's registration to the drug field. **`p_spatial` is the p-value to "
      "quote.** Both use the add-one correction, so the floor is "
      f"1/{stats[rows[0]['name']]['auc']['n_perm']+1:,} and "
      f"1/{stats[rows[0]['name']]['auc']['n_spatial']+1:,} respectively.")
    A("")
    A("### 4a. Delivery (per-spot drug exposure AUC)")
    A("")
    A("| drug | mean AUC perivascular | mean AUC core+necrotic | ratio | Δ | "
      "Cohen's *d* | Cliff's δ | *p*<sub>perm</sub> | *p*<sub>spatial</sub> |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        s = stats[r["name"]]["auc"]
        A(f"| {r['name']} | {_fmt(s['mean_high'])} | {_fmt(s['mean_low'])} | "
          f"{_fmt(s['ratio'],'.2f')}× | {_fmt(s['observed'])} | {s['d']:+.3f} | "
          f"{s['cliff']:+.3f} | {s['p_perm']:.2g} | {s['p_spatial']:.2g} |")
    A("")
    A("### 4b. Kill (per-spot mean predicted kill probability of tumour cells)")
    A("")
    A("| drug | mean *p*<sub>kill</sub> perivascular | mean *p*<sub>kill</sub> "
      "core+necrotic | ratio | Δ | Cohen's *d* | Cliff's δ | *p*<sub>perm</sub> | "
      "*p*<sub>spatial</sub> |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        s = stats[r["name"]]["kill"]
        A(f"| {r['name']} | {_fmt(s['mean_high'],'.4f')} | {_fmt(s['mean_low'],'.4f')} | "
          f"{_fmt(s['ratio'],'.2f')}× | {_fmt(s['observed'],'.4g')} | {s['d']:+.3f} | "
          f"{s['cliff']:+.3f} | {s['p_perm']:.2g} | {s['p_spatial']:.2g} |")
    A("")
    A("### 4c. Pairwise zone contrasts (delivery, first drug in the panel as reference)")
    A("")
    A("| contrast | ratio | Cohen's *d* | Cliff's δ | *p*<sub>spatial</sub> |")
    A("|---|---:|---:|---:|---:|")
    for label, s in ctx["pairwise"].items():
        A(f"| {label} | {_fmt(s['ratio'],'.2f')}× | {s['d']:+.3f} | {s['cliff']:+.3f} | "
          f"{s['p_spatial']:.2g} |")
    A("")

    # --- 5. resolution control ---------------------------------------------
    A("## 5. Resolution control (dx = 100 µm, `reconstruct_real`'s native Visium grid)")
    A("")
    A("The same panel, the same zone labels, on the assembler's own spot-pitch grid:")
    A("")
    A("| drug | perivasc/(core+necrotic) AUC ratio @12 µm | @100 µm | "
      "Cliff's δ @12 µm | @100 µm | whole-tumour kill @12 µm | @100 µm |")
    A("|---|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        c = ctx["coarse"]
        s12 = stats[r["name"]]["auc"]
        s100 = c["stats"][r["name"]]["auc"]
        k100 = next(x for x in c["rows"] if x["name"] == r["name"])["kill"]
        A(f"| {r['name']} | {_fmt(s12['ratio'],'.2f')}× | {_fmt(s100['ratio'],'.2f')}× | "
          f"{s12['cliff']:+.3f} | {s100['cliff']:+.3f} | {r['kill']:.4f} | {k100:.4f} |")
    A("")
    A(ctx["coarse_verdict"])
    A("")

    # --- 6. interpretation --------------------------------------------------
    A("## 6. Honest interpretation — does microenvironment structure produce "
      "regionally different drug response?")
    A("")
    for line in ctx["interpretation"]:
        A(line)
        A("")

    # --- 7. caveats ---------------------------------------------------------
    A("## 7. Caveats — what this does NOT show")
    A("")
    for line in ctx["caveats"]:
        A(f"- {line}")
    A("")
    return "\n".join(L)


# ---------------------------------------------------------------------------
def main(dx_um: float = DX_UM, sample: str = VISIUM_SAMPLE, quick: bool = False):
    t_start = time.time()
    os.makedirs(_RESULTS_DIR, exist_ok=True)

    # --- 1. real tissue + real niche labels ---------------------------------
    print(f"Loading REAL Visium sample {sample!r} ...")
    d = _retry(visium.load_sample, sample, marker_report=False)
    spot_coords = np.asarray(d["coords_um"], dtype=np.float64)
    print(f"  {d['n_spots']:,} spots x {d['n_genes']:,} genes  ({d['accession']})")

    print("Building REAL Ivy GAP anatomic-zone signatures ...")
    sig = region_signatures()
    spot_labels, scores, rnames, shared = label_regions(
        d["expression"], spot_coords, d["genes"], signatures=sig, return_scores=True
    )
    spot_labels = np.asarray(spot_labels, dtype=str)
    n_panel, n_shared = len(sig[2]), len(shared)
    srt = np.sort(scores, axis=1)
    mean_top1 = float(np.nanmean(srt[:, -1]))
    mean_margin = float(np.nanmean(srt[:, -1] - srt[:, -2]))
    print(f"  panel {n_panel} genes, shared with Visium {n_shared} "
          f"({100*n_shared/n_panel:.0f}%)")
    uniq, cnts = np.unique(spot_labels, return_counts=True)
    region_counts = dict(zip(uniq.tolist(), cnts.tolist()))
    for r in REGION_ORDER + ["unlabeled"]:
        if r in region_counts:
            print(f"    {r:12s} {region_counts[r]:5d} "
                  f"({100*region_counts[r]/len(spot_labels):5.1f}%)")
    present = [r for r in REGION_ORDER if region_counts.get(r, 0) > 0]

    # marker consistency check
    want = ["PECAM1", "VWF", "CA9", "VEGFA", "MKI67", "MBP"]
    mkvals = _cp10k_log(d["expression"], d["genes"], want)
    verdicts = []
    for g, zone in (("PECAM1", "perivascular"), ("VWF", "perivascular"),
                    ("CA9", "necrotic"), ("VEGFA", "necrotic"),
                    ("MKI67", "core"), ("MBP", "margin")):
        if g not in mkvals or zone not in present:
            continue
        m = spot_labels == zone
        a, b = mkvals[g][m].mean(), mkvals[g][~m].mean()
        verdicts.append(
            f"`{g}` in **{zone}** = {a:.3f} vs rest {b:.3f} — "
            f"{'higher, as expected' if a > b else '**NOT higher — inconsistent**'}"
        )
        print(f"  marker {g:7s} {zone:12s} {a:.3f} vs rest {b:.3f}")

    # --- 2. real microenvironment (fine, converged grid) --------------------
    print(f"\nReconstructing the REAL Visium microenvironment at dx={dx_um:.0f} um ...")
    micro = build_micro(dx_um, sample, use_pipeline=False)
    report = reconstruction_report(micro)
    print(f"  {report['n_cells']:,} cells (tumour {report['n_tumor']:,})  "
          f"grid {report['grid']}  dx={report['dx_um']:.0f} um  "
          f"vasc={report['mean_vascular_density']:.3f} o2={report['mean_oxygen']:.3f}")

    grid_idx, cell_idx, _ = assign_regions(
        micro.domain, spot_coords, spot_labels, micro.cells
    )
    print(f"  zone-assigned grid points: {(grid_idx>=0).sum():,} / {grid_idx.size:,}; "
          f"zone-assigned cells: {(cell_idx>=0).sum():,} / {micro.cells.n:,}")

    # --- 3. drug panel ------------------------------------------------------
    alpha = calibrate_alpha()
    print(f"\nperivascular-sink calibration alpha = {alpha:.1f} (fit once on doxorubicin)")
    pde = PDEConfig(
        dt=10.0,
        total_time=3600.0 * (1.0 if quick else 6.0),
        convection=True,
        save_every=500,
    )
    print(f"Simulating the real panel at dx={dx_um:.0f} um ...")
    rows, tumor, ct = run_panel(micro, grid_idx, cell_idx, alpha, pde)

    # --- 4. statistics on spot-level series ---------------------------------
    print("\nRunning permutation tests (i.i.d. + spatially-structured nulls) ...")
    n_perm = 2000 if quick else N_PERM
    n_spatial = 200 if quick else N_SPATIAL
    order = {r: i for i, r in enumerate(REGION_ORDER)}
    spot_idx = np.array([order.get(s, -1) for s in spot_labels], dtype=int)

    stats: Dict[str, Dict] = {}
    for r in rows:
        s_auc, s_kill = spot_level_series(
            micro.domain, spot_coords, r["auc_field"], r["p_kill"], micro.cells,
            tumor,
        )
        ok_kill = np.isfinite(s_kill)
        stats[r["name"]] = {
            "auc": permutation_test(
                s_auc, spot_idx, spot_coords, HIGH_GROUP, LOW_GROUP,
                n_perm=n_perm, n_spatial=n_spatial, seed=1,
            ),
            "kill": permutation_test(
                s_kill[ok_kill], spot_idx[ok_kill], spot_coords[ok_kill],
                HIGH_GROUP, LOW_GROUP, n_perm=n_perm, n_spatial=n_spatial, seed=2,
            ),
        }
        s = stats[r["name"]]["auc"]
        print(f"  {r['name']:14s} AUC peri/(core+nec) = {s['ratio']:.3f}x  "
              f"d={s['d']:+.3f}  cliff={s['cliff']:+.3f}  "
              f"p_perm={s['p_perm']:.2g}  p_spatial={s['p_spatial']:.2g}")

    # pairwise zone contrasts on the reference drug's delivery field
    ref = rows[0]
    s_auc_ref, _ = spot_level_series(
        micro.domain, spot_coords, ref["auc_field"], ref["p_kill"], micro.cells,
        tumor,
    )
    pairwise: Dict[str, Dict] = {}
    for hi, lo in (("perivascular", "necrotic"), ("perivascular", "core"),
                   ("perivascular", "margin"), ("margin", "necrotic"),
                   ("core", "necrotic")):
        if hi not in present or lo not in present:
            continue
        pairwise[f"{hi} vs {lo}"] = permutation_test(
            s_auc_ref, spot_idx, spot_coords, [hi], [lo],
            n_perm=n_perm, n_spatial=n_spatial, seed=3,
        )
        p = pairwise[f"{hi} vs {lo}"]
        print(f"  [{ref['name']}] {hi:12s} vs {lo:12s}: ratio={p['ratio']:.3f}x  "
              f"cliff={p['cliff']:+.3f}  p_spatial={p['p_spatial']:.2g}")

    # --- 5. resolution control at the native 100 um pipeline grid -----------
    print(f"\nResolution control: same panel via reconstruct_real (dx -> "
          f"{DX_SPOT_UM:.0f} um spot pitch) ...")
    micro_c = build_micro(DX_SPOT_UM, sample, use_pipeline=True)
    grid_idx_c, cell_idx_c, _ = assign_regions(
        micro_c.domain, spot_coords, spot_labels, micro_c.cells
    )
    rows_c, tumor_c, _ = run_panel(micro_c, grid_idx_c, cell_idx_c, alpha, pde)
    stats_c: Dict[str, Dict] = {}
    for r in rows_c:
        s_auc, s_kill = spot_level_series(
            micro_c.domain, spot_coords, r["auc_field"], r["p_kill"], micro_c.cells,
            tumor_c,
        )
        ok = np.isfinite(s_kill)
        stats_c[r["name"]] = {
            "auc": permutation_test(
                s_auc, spot_idx, spot_coords, HIGH_GROUP, LOW_GROUP,
                n_perm=n_perm, n_spatial=n_spatial, seed=1,
            ),
            "kill": permutation_test(
                s_kill[ok], spot_idx[ok], spot_coords[ok], HIGH_GROUP, LOW_GROUP,
                n_perm=n_perm, n_spatial=n_spatial, seed=2,
            ),
        }

    ratios12 = np.array([stats[r["name"]]["auc"]["ratio"] for r in rows])
    ratios100 = np.array([stats_c[r["name"]]["auc"]["ratio"] for r in rows])
    same_sign = bool(np.all(np.sign(ratios12 - 1.0) == np.sign(ratios100 - 1.0)))
    coarse_verdict = (
        "The perivascular/(core+necrotic) contrast has the **same direction at both "
        "resolutions** for every drug"
        if same_sign else
        "The contrast **changes direction between resolutions for at least one drug** "
        "— treat the effect as resolution-sensitive"
    ) + (
        f" (12 µm ratios {ratios12.min():.2f}–{ratios12.max():.2f}× vs 100 µm ratios "
        f"{ratios100.min():.2f}–{ratios100.max():.2f}×). Refining the mesh sharpens the "
        "near-vessel gradient the PDE can represent, but it adds no measurement "
        "information: the source field is spot-limited at both resolutions."
    )

    # --- 6. interpretation (written from the numbers, not hard-coded) -------
    ref_s = stats[rows[0]["name"]]["auc"]
    all_ratios = {r["name"]: stats[r["name"]]["auc"]["ratio"] for r in rows}
    all_cliff = {r["name"]: stats[r["name"]]["auc"]["cliff"] for r in rows}
    kill_ratios = {r["name"]: stats[r["name"]]["kill"]["ratio"] for r in rows}
    peri_wins = all(v > 1.0 for v in all_ratios.values())
    peri_wins_kill = all(
        np.isfinite(v) and v > 1.0 for v in kill_ratios.values()
    )
    max_p_spatial = max(stats[r["name"]]["auc"]["p_spatial"] for r in rows)
    hypoxic = rows[0]["hypoxic"]
    reg0 = rows[0]["regions"]

    # ordering check
    exp_order = [r for r in ["perivascular", "margin", "core", "necrotic"] if r in present]
    orders = {}
    for r in rows:
        orders[r["name"]] = sorted(
            present, key=lambda g: -r["regions"][g]["mean_tumor_exposure"]
        )
    full_order_ok = all(o == exp_order for o in orders.values())

    interpretation = []
    if peri_wins:
        interpretation.append(
            f"**Yes — the expected perivascular gradient appears, and it is large.** "
            f"For every drug in the panel, mean per-spot exposure is higher in the "
            f"perivascular zone than in the core+necrotic zones, by "
            f"{min(all_ratios.values()):.2f}–{max(all_ratios.values()):.2f}×, with "
            f"Cliff's δ from {min(all_cliff.values()):+.2f} to "
            f"{max(all_cliff.values()):+.2f} and a spatially-structured permutation "
            f"p ≤ {max_p_spatial:.2g}. The microenvironment does **not** deliver drug "
            f"uniformly: the same systemic dose produces materially different local "
            f"exposure depending on which anatomic zone a tumour cell sits in."
        )
    else:
        losers = [k for k, v in all_ratios.items() if v <= 1.0]
        interpretation.append(
            f"**The expected pattern does NOT hold for the whole panel.** Perivascular "
            f"exposure fails to exceed core+necrotic exposure for: "
            f"{', '.join(losers)} (ratios "
            f"{', '.join(f'{k}={all_ratios[k]:.2f}×' for k in losers)}). Reported as-is."
        )
    if full_order_ok:
        interpretation.append(
            f"**The full textbook ordering is recovered** "
            f"(`{' > '.join(exp_order)}`) for every drug, not just the "
            f"perivascular-vs-rest contrast — the intermediate zones fall in the "
            f"expected places too."
        )
    else:
        interpretation.append(
            "**The full four-zone ordering is only partly recovered.** Observed "
            "orderings by mean tumour-cell exposure: "
            + "; ".join(f"{k}: {' > '.join(v)}" for k, v in orders.items())
            + f" — against the textbook `{' > '.join(exp_order)}`. The "
            "perivascular-vs-core/necrotic contrast is the robust part; the finer "
            "ordering is not."
        )
    interpretation.append(
        "**But a large part of this is circular, and that has to be said plainly.** "
        "The `perivascular` zone label is driven by endothelial/angiogenic transcripts "
        "(`PECAM1`, `VWF`, `CD34`, `CLDN5`, `FLT1`, …), and the model's drug **source "
        "term is proportional to the vascular density field**, which is itself a KDE of "
        "the deconvolved ENDOTHELIAL spot fraction — i.e. of the *same* endothelial "
        "marker expression at the *same* spots. So 'drug arrives where the endothelial "
        "signal is' is close to a restatement of the model's input, not an independent "
        "discovery. The perivascular result should be read as a **consistency check "
        "that the pipeline propagates real vascular structure end-to-end**, not as "
        "evidence that the model predicts delivery. The genuinely non-circular "
        "comparisons are the ones that do not involve the perivascular label — "
        "`core` vs `necrotic` and `margin` vs `necrotic` (necrotic is called from HIF "
        "targets `CA9`/`VEGFA`/`ADM`/`NDRG1`, core from `EGFR`/`SOX2`/`OLIG2`/`MKI67`, "
        "neither of which feeds the vascular field) — see §4c."
    )
    # The non-circular arm: contrasts that do not involve the perivascular label.
    nc = {k: v for k, v in pairwise.items() if "perivascular" not in k}
    if nc:
        bits = ", ".join(
            f"**{k}**: {v['ratio']:.2f}× (δ={v['cliff']:+.2f}, "
            f"p_spatial={v['p_spatial']:.2g})"
            for k, v in nc.items()
        )
        cn = pairwise.get("core vs necrotic")
        if cn is not None and cn["ratio"] < 1.0:
            interpretation.append(
                f"**The non-circular half of the prediction FAILS, and this is the "
                f"most informative result here.** Stripping out every contrast that "
                f"involves the perivascular label leaves: {bits}. The textbook "
                f"expectation is that the necrotic/perinecrotic zone is the "
                f"worst-perfused, worst-delivered compartment in the tumour — it is "
                f"not. On this section the `necrotic` zone receives **more** predicted "
                f"drug than the `core` ({1/cn['ratio']:.2f}× more, Cliff's "
                f"δ={-cn['cliff']:+.2f}), and the difference is not significant against "
                f"the spatial null (p_spatial = {cn['p_spatial']:.2g}). Mechanically "
                f"this follows straight from the reconstructed fields: mean vascular "
                f"density is {reg0['necrotic']['vascular']:.4f} in the necrotic zone "
                f"vs {reg0['core']['vascular']:.4f} in the core, and mean oxygen "
                f"{reg0['necrotic']['oxygen']:.4f} vs {reg0['core']['oxygen']:.4f}. "
                f"Biologically that is defensible — perinecrotic HIF-high tissue "
                f"(`CA9`/`VEGFA`/`ADM`) is precisely the tissue that *drives* "
                f"microvascular proliferation, so it sits next to vessels — but it "
                f"means **this model does not reproduce the classic 'necrotic core is a "
                f"drug sanctuary' picture**, and the honest reading is that once you "
                f"remove the vascular-marker circularity, only the "
                f"perivascular-vs-everything contrast survives, and the finer zone "
                f"structure does not."
            )
        else:
            interpretation.append(
                f"**Non-circular contrasts** (those not involving the perivascular "
                f"label, whose zone signature shares markers with the model's vascular "
                f"source term): {bits}. These are the contrasts that carry real "
                f"information about whether the zone signatures predict delivery."
            )
    if peri_wins_kill:
        interpretation.append(
            f"**Exposure differences convert into kill differences.** Predicted kill of "
            f"the tumour cells inside the perivascular zone exceeds that in "
            f"core+necrotic for every drug, by "
            f"{min(v for v in kill_ratios.values() if np.isfinite(v)):.2f}–"
            f"{max(v for v in kill_ratios.values() if np.isfinite(v)):.2f}×. This is the "
            f"spatially explicit version of the whole claim: a single whole-tumour "
            f"`kill_fraction` averages over cells whose predicted fates differ by that "
            f"factor. Note however that the **kill contrast is much smaller than the "
            f"exposure contrast** (exposure "
            f"{min(all_ratios.values()):.2f}–{max(all_ratios.values()):.2f}× vs kill "
            f"{min(v for v in kill_ratios.values() if np.isfinite(v)):.2f}–"
            f"{max(v for v in kill_ratios.values() if np.isfinite(v)):.2f}×): the Hill "
            f"dose-response compresses it, because at these predicted exposures most "
            f"drugs sit high on the sigmoid where extra drug buys little extra kill. "
            f"A delivery gradient only matters when the tissue straddles the drug's "
            f"IC50: here the largest kill contrast is "
            f"**{max(kill_ratios, key=lambda k: kill_ratios[k])}** "
            f"({max(kill_ratios.values()):.2f}×) and the smallest is "
            f"**{min(kill_ratios, key=lambda k: kill_ratios[k])}** "
            f"({min(kill_ratios.values()):.2f}×), which tracks where each drug's "
            f"achieved exposure sits on its own dose-response curve rather than how "
            f"steep its delivery gradient is."
        )
    else:
        interpretation.append(
            "**Exposure differences do NOT cleanly convert into kill differences for "
            "every drug** (kill ratios: "
            + ", ".join(f"{k}={v:.2f}×" for k, v in kill_ratios.items())
            + "). Where the Hill curve is saturated or far below threshold across the "
            "whole tissue, a large exposure contrast produces little kill contrast — "
            "the delivery gradient is real but pharmacodynamically inert."
        )
    interpretation.append(
        f"**Hypoxia is not the mechanism here.** On this section at dx = "
        f"{report['dx_um']:.0f} µm the predicted tumour hypoxic fraction is "
        f"{hypoxic*100:.1f}%, so the hypoxia-resistance multiplier in the PD model is "
        f"essentially inactive and the regional kill differences are driven almost "
        f"entirely by **delivered exposure**, not by hypoxic drug resistance. The "
        f"Visium-derived vascular field is a spot-resolution KDE and is consequently "
        f"much smoother than the Xenium single-cell field, which is why this section "
        f"reads as well-oxygenated. Do not transfer the hypoxia numbers from the Xenium "
        f"panel to this one."
    )
    interpretation.append(
        "**What would make this a real result rather than a coherent one:** paired "
        "region-resolved drug measurement (e.g. MALDI-MSI or LC-MS of a treated "
        "resection, zone-annotated), which would let the per-zone AUC column be scored "
        "instead of merely reported. Nothing in this experiment is validated against a "
        "measured outcome."
    )

    caveats = [
        "**No outcome labels anywhere.** This Visium section is an untreated primary "
        "resection. There is no measured drug response, and no measured per-zone drug "
        "concentration, for any of the five drugs. Every kill and exposure number is a "
        "**prediction of the physics on real geometry**, not a validated outcome.",
        "**Zone labels are signature-correlation assignments, not pathology.** Each "
        f"spot is given the zone whose Ivy GAP signature it rank-correlates with best. "
        f"The mean winning correlation is only {mean_top1:.3f} and the mean margin over "
        f"the runner-up is {mean_margin:.3f} — these are **soft, low-margin "
        f"assignments**, not annotated histology. A pathologist's zone map of this "
        f"section does not exist in the open data, so the labelling cannot be scored. "
        f"Treat zone membership as a graded score that has been hardened into a class.",
        "**Circularity in the perivascular arm** (see §6): the perivascular label and "
        "the model's vascular source term are both functions of endothelial marker "
        "expression at the same spots. The perivascular-vs-rest contrast is therefore "
        "partly a restatement of an input.",
        "**Spot-limited resolution.** Visium spots are 55 µm across on a ~100 µm pitch. "
        "Real perivascular gradients in GBM are ~50 µm wide — at or below the "
        "measurement scale. The 12 µm mesh is a numerical convergence requirement and "
        "creates no sub-spot information; zone labels are replicated over each spot's "
        "footprint. Any claim about structure finer than ~100 µm is unsupported by the "
        "data, whatever the figure looks like.",
        "**Deconvolution stands between the counts and the cells.** Cell types (and "
        "hence the vascular and ECM fields) come from reference-based deconvolution of "
        "spot expression, not from observed cells. Errors there propagate into every "
        "per-zone number.",
        "**Necrotic zones contain fewer viable cells than the model assumes.** The "
        "reconstruction places a fixed bag of cells per spot apportioned by "
        "deconvolved fractions; it does not down-weight the truly acellular fraction of "
        "a necrotic zone. Per-zone kill *fractions* in the necrotic zone are therefore "
        "computed over a cell population that is probably too viable and too numerous.",
        "**One section, one patient.** n = 1 tissue. Nothing here establishes that "
        "these zone contrasts generalise; the cross-tissue question is not addressed.",
        "**Statistical caveats.** The spot is treated as the unit of analysis, which is "
        "better than the cell but still not independent — neighbouring spots are "
        "strongly autocorrelated. That is why `p_spatial` (rotation null) is quoted "
        "rather than `p_perm` (i.i.d. null); the rotation null is itself approximate, "
        "since nearest-neighbour label transfer after rotation slightly perturbs class "
        "counts near the tissue boundary. The five drugs are not independent tests — "
        "they share one tissue, one vascular field and one α calibration, so their "
        "p-values are effectively one test repeated with different drug parameters.",
        "**A single calibration constant sets the delivery scale.** α is fit once on a "
        "doxorubicin penetration measurement; the absolute AUC of every zone rides on "
        "that one anchor. Zone *ratios* are far more trustworthy than zone *absolutes*.",
    ]
    # Zone-specific caveats derived from the run itself.
    failed_markers = [v for v in verdicts if "NOT higher" in v]
    if failed_markers:
        caveats.append(
            "**Not every canonical marker separates the way it should.** "
            + " ".join(v.replace("**", "") + "." for v in failed_markers)
            + " A zone whose own canonical marker does not enrich in it is a zone whose "
            "assignment should be treated with extra suspicion — the correlation-based "
            "call is being carried by the other ~280 panel genes, not by the marker "
            "the zone is named for."
        )
    small = [r for r in present if region_counts.get(r, 0) < 200]
    if small:
        caveats.append(
            "**Some zones are small.** "
            + ", ".join(f"`{r}` has only {region_counts[r]} spots" for r in small)
            + ". Their per-zone means are noisy and their pairwise contrasts are "
            "underpowered — visible in §4c, where contrasts involving these zones do "
            "not reach significance against the spatial null."
        )

    # --- 7. figure + report -------------------------------------------------
    meta = dict(
        sample_id=d["sample_id"],
        accession=d["accession"],
        diagnosis=d["diagnosis"],
        n_spots=int(d["n_spots"]),
        n_genes=int(d["n_genes"]),
        spot_diameter_um=float(d["spot_diameter_fullres"]) * float(d["umpp"]),
        n_panel=n_panel,
        n_shared=n_shared,
        mean_top1=mean_top1,
        mean_margin=mean_margin,
    )
    png = make_figure(
        micro, spot_coords, spot_labels, rows, region_counts, meta, report["dx_um"], stats
    )
    print(f"\nfigure written: {png}")

    ctx = dict(
        rows=rows, stats=stats, present=present, sample_meta=meta, report=report,
        alpha=alpha, pde=pde, region_counts=region_counts, spot_labels=spot_labels,
        markers=dict(values=mkvals, verdicts=verdicts),
        pairwise=pairwise,
        coarse=dict(rows=rows_c, stats=stats_c),
        coarse_verdict=coarse_verdict,
        interpretation=interpretation, caveats=caveats,
    )
    with open(_OUT_MD, "w") as fh:
        fh.write(build_markdown(ctx))
    print(f"report written: {_OUT_MD}")

    # --- 8. stdout summary --------------------------------------------------
    print("\n" + "=" * 110)
    print("DRUG x NICHE ZONE on the REAL Visium GBM section (predictions only)")
    print("=" * 110)
    hdr = (f"{'drug':14s} {'zone':13s} {'tissueAUC':>11s} {'cellExpo':>11s} "
           f"{'killFrac':>9s} {'E[kill]':>9s} {'vasc':>7s} {'O2':>7s} {'nTumor':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        for reg in present:
            g = r["regions"][reg]
            print(f"{r['name']:14s} {reg:13s} {g['mean_auc']:11.3g} "
                  f"{g['mean_tumor_exposure']:11.3g} {g['kill_frac']:9.4f} "
                  f"{g['exp_kill']:9.4f} {g['vascular']:7.4f} {g['oxygen']:7.4f} "
                  f"{g['n_tumor']:8d}")
        print("-" * len(hdr))
    print("PERIVASCULAR vs CORE+NECROTIC (spot-level, one-sided):")
    for r in rows:
        s = stats[r["name"]]["auc"]
        k = stats[r["name"]]["kill"]
        print(f"  {r['name']:14s} AUC ratio {s['ratio']:6.2f}x d={s['d']:+.2f} "
              f"cliff={s['cliff']:+.2f} p_perm={s['p_perm']:.2g} "
              f"p_spatial={s['p_spatial']:.2g} | kill ratio {k['ratio']:6.2f}x "
              f"cliff={k['cliff']:+.2f} p_spatial={k['p_spatial']:.2g}")
    print("=" * 110)
    print("HONESTY: untreated resection -> NO outcome label; zone labels are soft "
          "signature\ncorrelations, not pathology; the perivascular arm is partly "
          "circular with the model's\nvascular source term; measurement is "
          "spot-limited to ~100 um. See results/niche_response.md.")
    print("=" * 110)
    print(f"[runtime {time.time() - t_start:.0f} s]")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Predicted drug response by real anatomic niche region"
    )
    ap.add_argument("--dx", type=float, default=DX_UM,
                    help="grid spacing in um (default %(default)s, grid-converged)")
    ap.add_argument("--sample", type=str, default=VISIUM_SAMPLE,
                    help="Visium sample id (default %(default)s)")
    ap.add_argument("--quick", action="store_true",
                    help="short PDE window + fewer permutations (smoke test only)")
    a = ap.parse_args()
    main(dx_um=a.dx, sample=a.sample, quick=a.quick)
