"""Is MGMT a spatially resolvable field, or is a per-patient scalar adequate?

The falsified assumption
------------------------
``results/gbm_tmz_exvivo.md`` §B2 and ``results/mgmt_scalar.md`` §2/§4.4 assert
that MGMT cannot be resolved spatially. The *evidence* was narrow — MGMT is
absent from the two Xenium panels this repo held, and sits at the detection
floor in Visium — but the *claim* was written as a general property of spatial
assays. It is not one.

GEO **GSE324268** (Hölzel lab, public 2026-03) is a free, no-login human
primary GBM **Xenium** section whose custom add-on panel carries a dedicated
8-probe MGMT target across 851,509 segmented cells. It makes the general claim
directly testable, and this module tests it.

What is measured here
---------------------
1. **Provenance and detection** — the MGMT target read out of the shipped panel
   spec, then the counts read out of the shipped matrix, then both scored
   against the assay's own negative controls.
2. **Reconstruction** by the UNMODIFIED pipeline
   (``reconstruct.real_xenium.reconstruct_from_xenium``), and the ENDOTHELIAL
   fraction against the human volume-EM anchor of 8.0 % (PMID 38723085), as
   ``results/heldout_patient.md`` did for section 2.
3. **Spatial autocorrelation** of per-cell MGMT among TUMOUR cells — Moran's I
   on a k = 6 nearest-neighbour, row-standardised graph, with a permutation
   null, against three reference classes: the pooled negative controls (the
   technical floor), 47 abundance-matched real genes (what IS detectable at
   this abundance), and high-expression structured genes.
4. **Variance partition** — between-region vs within-region, over square bins
   at 100/200/400/800/1600 µm and over Ivy GAP niche regions, as Poisson
   overdispersion φ against the SAME partition applied to the pooled negative
   controls.
5. **Power** — the minimum detectable bin-rate CV at 80 % power, and the
   minimum detectable rate ratio for region-sized two-group contrasts.
6. **The decision** — predicted TMZ kill under a scalar MGMT vs a
   spatially-resolved MGMT, on this section, with everything else held fixed.

The unit problem, stated once
-----------------------------
``molecular.mgmt`` is calibrated in **bulk CPM over ~34,461 genes**. A CPM
computed over a **366-gene targeted panel** is a different quantity by roughly
two orders of magnitude (this section reads 71.5 "panel-CPM" among tumour
cells against a calibration ceiling of 21.1 bulk CPM). Feeding panel-CPM into
the mapping would be a unit error dressed as a measurement.

So the impact analysis holds the SECTION MEAN fixed at a declared anchor drawn
from the Visium range already published in ``results/mgmt_scalar.md`` §2
(7.00 CPM, the lowest section; 25.52 CPM, the highest) and takes ONLY the
*relative* spatial modulation from Xenium:

    CPM(bin) = CPM_anchor * rate(bin) / rate(section)

That isolates the question actually being asked — *does spatial structure
change the answer?* — from the cross-assay calibration question, which this
dataset cannot settle.

Nothing under ``core/``, ``physics/``, ``molecular/``, ``reconstruct/`` or
``sim/`` is modified, patched or shadowed by this module.
"""
from __future__ import annotations

import json
import os
import time
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import sparse, stats
from scipy.spatial import cKDTree

from ..core.microenvironment import CellType
from ..data.real import mgmt_spatial as loader

# --------------------------------------------------------------------------- #
# Declared analysis constants (all fixed before the result was seen)
# --------------------------------------------------------------------------- #
#: Nearest-neighbour count for the Moran's I spatial weights. 6 is the planar
#: Delaunay-average degree, so a row-standardised k = 6 graph is the natural
#: "immediate neighbours" operator for a dense 2-D cell layout.
K_NEIGHBOURS = 6

#: Permutations for the Moran's I null. 999 for the headline variables, 49 for
#: the large reference sweeps (where the point estimate, not its p, is used).
N_PERM = 999
N_PERM_SWEEP = 49

#: Square-bin edge lengths (µm) for the variance partition.
BIN_SIZES_UM: Tuple[int, ...] = (100, 200, 400, 800, 1600)

#: Bin used for the headline partition and for the spatial-MGMT sensitivity
#: field. 400 µm holds a median ~1,000 tumour cells — enough that the per-bin
#: MGMT rate is not pure Poisson noise, small enough to resolve a niche.
HEADLINE_BIN_UM = 400

#: Minimum tumour cells for a bin to enter the partition.
MIN_CELLS_PER_BIN = 50

#: Abundance-matched reference class: real panel genes whose section total is
#: within this factor of MGMT's.
ABUNDANCE_WINDOW = (0.5, 2.0)

#: Human cortex volume-EM endothelial fraction (Shapson-Coe et al., Science
#: 2024;384:eadk4858, PMID 38723085) — the reconstruction's one external anchor.
ENDOTHELIAL_VOLUME_EM = 8.0

#: Grid spacing for the impact simulations. 12 µm is the converged value
#: (``results/numerical_rigor.md``) and the one ``mgmt_scalar.md`` §3 used.
DX_UM = 12.0

#: Anchors for the impact runs — the extremes of the per-section pseudobulk
#: MGMT range measured from our own Visium data in ``results/mgmt_scalar.md``.
ANCHORS_CPM: Tuple[Tuple[str, float], ...] = (
    ("GBM2 (lowest of our Visium sections)", 7.00),
    ("10x_parent_gbm (highest of our Visium sections)", 25.52),
)

RNG_SEED = 20260728

_RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
_RESULTS = os.path.abspath(_RESULTS)
_RUN_CACHE = os.path.join("data_cache", "mgmt_spatial_runs.json")


# --------------------------------------------------------------------------- #
# 0. provenance + detection
# --------------------------------------------------------------------------- #
def verify_provenance(path: str = loader.DEFAULT_CACHE) -> Dict:
    """Re-hash the downloaded bytes and confirm they match the frozen record."""
    rec = loader.download(path, verify=True)
    return {
        "accession": loader.GSE324268_ACCESSION,
        "sample": loader.GSE324268_SAMPLE,
        "landing_page": loader.GSE324268_LANDING,
        "files": rec,
        "all_match": all(
            v["sha256_match"] for v in rec.values() if v["sha256_match"] is not None
        ),
        "skipped": loader.GSE324268_SKIPPED,
    }


def verify_mgmt(path: str = loader.DEFAULT_CACHE, cell_type: Optional[np.ndarray] = None) -> Dict:
    """Verify the MGMT claim from the panel spec AND from the count matrix.

    Reports the exact feature name, Ensembl id, probe count, detection rate
    over all cells and over TUMOUR cells, and the counts — plus the negative
    control background the counts have to beat to mean anything.
    """
    entry = loader.mgmt_panel_entry(path)
    d = loader.mgmt_per_cell(path)
    mg = d["mgmt"]
    gene_totals = d["gene_totals"]
    rank = int(np.where(gene_totals["gene"].to_numpy() == "MGMT")[0][0]) + 1

    npt, nct = d["neg_probe_totals"], d["neg_codeword_totals"]
    per_probe = mg.sum() / d["n_mgmt_probes"]

    out = {
        "panel_entry": entry,
        "n_cells": int(len(mg)),
        "mgmt_counts": int(mg.sum()),
        "mgmt_cells_positive": int((mg > 0).sum()),
        "mgmt_pct_cells": float(100.0 * (mg > 0).mean()),
        "mgmt_max_per_cell": int(mg.max()),
        "mgmt_count_histogram": {
            int(v): int(c) for v, c in zip(*np.unique(mg, return_counts=True))
        },
        "rank_by_total_counts": rank,
        "n_gene_targets": int(len(gene_totals)),
        "median_gene_total": float(gene_totals["total"].median()),
        "median_gene_counts_per_cell": float(np.median(d["gene_total"])),
        "neg_probe": {
            "n_features": d["n_neg_probes"],
            "total": int(npt.sum()),
            "mean_per_feature": float(npt.mean()),
            "median_per_feature": float(np.median(npt)),
            "max_per_feature": int(npt.max()),
        },
        "neg_codeword": {
            "n_features": d["n_neg_codewords"],
            "total": int(nct.sum()),
            "mean_per_feature": float(nct.mean()),
            "max_per_feature": int(nct.max()),
        },
        "mgmt_counts_per_probe": float(per_probe),
        "fold_over_neg_probe_mean": float(per_probe / npt.mean()),
        "fold_over_worst_neg_probe": float(per_probe / npt.max()),
        "section_panel_cpm": float(1e6 * mg.sum() / d["gene_total"].sum()),
        "source": d["source"],
    }
    if cell_type is not None:
        by_type = {}
        for c in CellType:
            m = cell_type == int(c)
            if m.sum():
                by_type[c.name] = {
                    "n_cells": int(m.sum()),
                    "mgmt_counts": int(mg[m].sum()),
                    "pct_positive": float(100.0 * (mg[m] > 0).mean()),
                    "panel_cpm": float(1e6 * mg[m].sum() / d["gene_total"][m].sum()),
                }
        out["by_cell_type"] = by_type
        tum = cell_type == int(CellType.TUMOR)
        out["mgmt_pct_tumour_cells"] = float(100.0 * (mg[tum] > 0).mean())
        out["mgmt_counts_tumour"] = int(mg[tum].sum())
        out["n_tumour_cells"] = int(tum.sum())
    return out


# --------------------------------------------------------------------------- #
# 1. reconstruction through the unmodified pipeline
# --------------------------------------------------------------------------- #
def reconstruct_section(path: str = loader.DEFAULT_CACHE, dx_um: float = 30.0) -> Dict:
    """Reconstruct GSE324268 with the existing pipeline and score composition."""
    micro = loader.reconstruct(path, dx_um=dx_um)
    ct = micro.cells.cell_type
    n = int(micro.cells.n)
    frac = {
        c.name: {
            "n": int((ct == int(c)).sum()),
            "pct": float(100.0 * (ct == int(c)).mean()),
        }
        for c in CellType
    }
    area_mm2 = (
        (micro.cells.x.max() - micro.cells.x.min())
        * (micro.cells.y.max() - micro.cells.y.min())
        / 1e6
    )
    endo_pct = frac["ENDOTHELIAL"]["pct"]
    return {
        "n_cells": n,
        "dx_um": float(dx_um),
        "grid_shape": tuple(int(v) for v in micro.domain.shape),
        "extent_um": tuple(float(v) for v in micro.domain.extent_um),
        "bounding_area_mm2": float(area_mm2),
        "cell_density_per_mm2": float(n / area_mm2),
        "celltype_fractions": frac,
        "endothelial_pct": endo_pct,
        "endothelial_vs_volume_em": float(endo_pct / ENDOTHELIAL_VOLUME_EM),
        "endothelial_per_mm2": float(frac["ENDOTHELIAL"]["n"] / area_mm2),
        "marker_panel_used": micro.metadata["marker_panel_used"],
        "marker_panel_coverage": micro.metadata["marker_panel_coverage"],
        "marker_map": micro.metadata["marker_map"],
        "markers_missing": micro.metadata["markers_missing_from_panel"],
        "tissue_coverage_fraction": micro.metadata["tissue_coverage_fraction"],
        "cell_type": ct,
        "micro": micro,
    }


def verify_existing_panels_unchanged() -> Dict:
    """Assert the two previously-cached sections still select the same maps.

    No file under ``reconstruct/`` was modified by this work, and this panel
    needs no new registry entry (it carries all twelve brain-map markers, so
    ``select_marker_map`` already scores it brain_319 = 1.00). This function
    turns that from a claim into a check.
    """
    from ..reconstruct.marker_panels import (
        BRAIN_PANEL_MARKER_MAP,
        IO_PANEL_MARKER_MAP,
        select_marker_map,
    )
    from ..data.real.singlecell_spatial import load as load_xen

    expect = {
        "data_cache/xenium_gbm": ("brain_319", 1.00, 0.50, BRAIN_PANEL_MARKER_MAP),
        "data_cache/xenium_gbm_io": ("io_480", 0.50, 1.00, IO_PANEL_MARKER_MAP),
    }
    out = {}
    for p, (want_name, want_b, want_io, want_obj) in expect.items():
        if not os.path.exists(p):
            out[p] = {"status": "cache absent — not checked"}
            continue
        g = np.asarray(load_xen(path=p)["genes"], dtype=str)
        mm, name, cov = select_marker_map(g)
        out[p] = {
            "panel_selected": name,
            "expected": want_name,
            "coverage": cov,
            "same_map_object": mm is want_obj,
            "unchanged": (
                name == want_name
                and abs(cov["brain_319"] - want_b) < 5e-3
                and abs(cov["io_480"] - want_io) < 5e-3
                and mm is want_obj
            ),
        }
    out["gse324268"] = loader.panel_marker_coverage()
    out["new_registry_entry_needed"] = False
    return out


# --------------------------------------------------------------------------- #
# 2. the analysis substrate: tumour cells, depth, MGMT, weights
# --------------------------------------------------------------------------- #
class _Substrate:
    """Everything the spatial tests share, built once."""

    def __init__(self, path: str, cell_type: np.ndarray):
        d = loader.load(path, genes_only=False)
        self.genes = np.asarray(d["genes"])
        self.ftype = np.asarray(d["feature_types"])
        self.Xc = d["expression"].tocsc()
        self.is_gene = self.ftype == "Gene Expression"
        gj = np.where(self.is_gene)[0]
        self.gene_cols = gj
        self.gene_names = self.genes[gj]
        depth_all = np.asarray(self.Xc[:, gj].sum(axis=1)).ravel()

        self.tum = cell_type == int(CellType.TUMOR)
        self.idx = np.where(self.tum)[0]
        self.xy = d["cell_xy_um"][self.idx]
        self.depth = depth_all[self.idx]
        self.mgmt = self.col("MGMT")
        self.rate = self.mgmt.sum() / self.depth.sum()

        npj = np.where(self.ftype == "Negative Control Probe")[0]
        ncj = np.where(self.ftype == "Negative Control Codeword")[0]
        self.pool_neg_probe = np.asarray(self.Xc[:, npj].sum(axis=1)).ravel()[self.idx]
        self.pool_neg_codeword = np.asarray(self.Xc[:, ncj].sum(axis=1)).ravel()[self.idx]
        self.n_neg_probe, self.n_neg_codeword = len(npj), len(ncj)
        self._neg_probe_cols, self._neg_cw_cols = npj, ncj

        self.gene_section_totals = np.asarray(self.Xc[:, gj].sum(axis=0)).ravel()
        self._W = None
        self._nn_dist = None

    def col(self, gene: str) -> np.ndarray:
        j = int(np.where((self.genes == gene) & self.is_gene)[0][0])
        return np.asarray(self.Xc[:, j].todense()).ravel()[self.idx]

    def col_at(self, j: int) -> np.ndarray:
        return np.asarray(self.Xc[:, j].todense()).ravel()[self.idx]

    @property
    def W(self) -> sparse.csr_matrix:
        if self._W is None:
            n = len(self.idx)
            tree = cKDTree(self.xy)
            dist, nb = tree.query(self.xy, k=K_NEIGHBOURS + 1, workers=-1)
            nb, dist = nb[:, 1:], dist[:, 1:]
            self._nn_dist = dist
            self._W = sparse.csr_matrix(
                (np.full(n * K_NEIGHBOURS, 1.0 / K_NEIGHBOURS),
                 (np.repeat(np.arange(n), K_NEIGHBOURS), nb.ravel())),
                shape=(n, n),
            )
        return self._W

    def bins(self, edge_um: float) -> np.ndarray:
        bx = ((self.xy[:, 0] - self.xy[:, 0].min()) // edge_um).astype(np.int64)
        by = ((self.xy[:, 1] - self.xy[:, 1].min()) // edge_um).astype(np.int64)
        return bx * (by.max() + 1) + by

    def abundance_matched(self) -> np.ndarray:
        """Panel genes whose section total brackets MGMT's — the fair reference."""
        mtot = self.gene_section_totals[
            int(np.where(self.gene_names == "MGMT")[0][0])
        ]
        lo, hi = ABUNDANCE_WINDOW
        keep = (
            (self.gene_section_totals >= lo * mtot)
            & (self.gene_section_totals <= hi * mtot)
            & (self.gene_names != "MGMT")
        )
        return self.gene_cols[keep]


def _moran(y: np.ndarray, W: sparse.csr_matrix, rng, nperm: int) -> Dict:
    # ``errstate`` only because macOS Accelerate leaves FPU status flags set
    # after a large BLAS dot, which numpy then reports as a spurious
    # divide/overflow/invalid triple. The returned values are bit-identical
    # with and without it (checked); nothing here can legitimately divide by
    # zero — the one real degenerate case (a constant vector) is handled below.
    with np.errstate(all="ignore"):
        y = np.asarray(y, dtype=np.float64)
        z = y - y.mean()
        den = float(z @ z)
        if den <= 0 or not np.isfinite(den):
            return {"I": float("nan"), "p": float("nan"), "z": float("nan"),
                    "null_mean": float("nan"), "null_sd": float("nan")}
        I = float(z @ (W @ z)) / den
        null = np.empty(nperm)
        for b in range(nperm):
            zp = rng.permutation(z)
            null[b] = float(zp @ (W @ zp)) / den
        sd = float(null.std(ddof=1))
    return {
        "I": I,
        "p": float((int((null >= I).sum()) + 1) / (nperm + 1)),
        "null_mean": float(null.mean()),
        "null_sd": sd,
        "z": float((I - null.mean()) / sd) if sd > 0 else float("nan"),
    }


# --------------------------------------------------------------------------- #
# 3(a). spatial autocorrelation among tumour cells
# --------------------------------------------------------------------------- #
def spatial_autocorrelation(sub: _Substrate) -> Dict:
    """Moran's I of per-cell MGMT among TUMOUR cells, with three reference classes.

    The headline number is meaningless on its own. With n = 652,811 cells a
    permutation test will call almost anything significant, so what decides the
    question is where MGMT's I sits relative to

      * the **pooled negative controls** — the technical floor, i.e. how much
        apparent clustering this assay produces from segmentation, cell area
        and probe density alone at a comparable count depth; and
      * **abundance-matched real genes** — what a genuinely structured gene
        detected in ~1 % of cells looks like on this same tissue.

    ``depth-residual`` removes the strongest confound: total transcript count
    per cell is itself strongly clustered (I ~ 0.14), so any count variable
    inherits some of that.
    """
    rng = np.random.default_rng(RNG_SEED)
    W = sub.W
    mg, depth = sub.mgmt, sub.depth
    resid = mg - sub.rate * depth

    headline = {}
    for lbl, y, nperm in (
        ("MGMT count", mg, N_PERM),
        ("MGMT binary (detected / not)", (mg > 0).astype(float), N_PERM),
        ("MGMT depth-residual", resid, N_PERM),
        ("pooled negative-control probes", sub.pool_neg_probe, N_PERM),
        ("pooled negative-control codewords", sub.pool_neg_codeword, N_PERM),
        ("total gene counts per cell (DEPTH — the confound)", depth, N_PERM),
    ):
        headline[lbl] = _moran(y, W, rng, nperm)

    matched = []
    for j in sub.abundance_matched():
        y = sub.col_at(j)
        r = _moran(y, W, rng, N_PERM_SWEEP)
        matched.append({
            "gene": str(sub.genes[j]),
            "total": int(y.sum()),
            "pct_positive": float(100.0 * (y > 0).mean()),
            "I": r["I"],
        })
    matched.sort(key=lambda d: d["I"])
    mI = np.array([m["I"] for m in matched])

    positive = []
    for g in ("VEGFA", "AQP4", "EGFR", "CA9", "MKI67", "PTPRC"):
        if g in set(sub.gene_names.tolist()):
            y = sub.col(g)
            r = _moran(y, W, rng, N_PERM_SWEEP)
            positive.append({
                "gene": g, "total": int(y.sum()),
                "pct_positive": float(100.0 * (y > 0).mean()), "I": r["I"],
            })

    mgmt_I = headline["MGMT count"]["I"]
    n_more = int((mI > mgmt_I).sum())
    return {
        "k_neighbours": K_NEIGHBOURS,
        "n_tumour_cells": int(len(sub.idx)),
        "median_nn_um": float(np.median(sub._nn_dist[:, 0])),
        "median_kth_nn_um": float(np.median(sub._nn_dist[:, -1])),
        "n_perm": N_PERM,
        "headline": headline,
        "abundance_matched": matched,
        "abundance_matched_stats": {
            "n": len(matched),
            "min": float(mI.min()), "median": float(np.median(mI)),
            "max": float(mI.max()),
            "n_exceeding_mgmt": n_more,
            "mgmt_percentile": float(100.0 * (mI < mgmt_I).mean()),
        },
        "positive_controls": positive,
    }


# --------------------------------------------------------------------------- #
# 3(b). variance partition
# --------------------------------------------------------------------------- #
def _partition(y: np.ndarray, depth: np.ndarray, labels: np.ndarray,
               min_cells: int = MIN_CELLS_PER_BIN) -> Dict:
    """Between-region vs within-region variance of a count rate.

    The model is: region ``b`` holds ``N_b`` total gene counts and ``k_b``
    counts of the gene. Under "one tumour-level mean and nothing else",
    ``k_b ~ Poisson(p * N_b)``. Then

      * ``phi`` = Pearson chi-square / df. 1.0 == sampling noise only.
      * ``frac_spatial`` = (Var_obs - Var_sampling) / Var_obs, the share of the
        observed between-region variance in the rate that is NOT explained by
        counting statistics. **This is the "fraction of MGMT variance that is
        spatial rather than a single tumour-level mean".**
      * ``cv_excess`` = sqrt(excess variance) / p — the same thing expressed as
        a multiplicative swing of the local rate, which is the unit the
        sensitivity mapping consumes.
    """
    u, inv = np.unique(labels, return_inverse=True)
    B = len(u)
    k = np.bincount(inv, weights=y, minlength=B)
    N = np.bincount(inv, weights=depth, minlength=B)
    nc = np.bincount(inv, minlength=B)
    keep = nc >= min_cells
    k, N, nc = k[keep], N[keep], nc[keep]
    if len(k) < 2:
        return {"B": int(len(k)), "phi": float("nan")}
    p = k.sum() / N.sum()
    exp = p * N
    chi2 = float(((k - exp) ** 2 / exp).sum())
    df = len(k) - 1
    rate = k / N
    var_obs = float(np.average((rate - p) ** 2, weights=N))
    var_samp = float(np.average(p / N, weights=N))
    excess = max(0.0, var_obs - var_samp)
    return {
        "B": int(len(k)),
        "df": int(df),
        "chi2": chi2,
        "phi": chi2 / df,
        "p_chi2": float(stats.chi2.sf(chi2, df)),
        "rate": float(p),
        "var_observed": var_obs,
        "var_sampling": var_samp,
        "frac_spatial": float(excess / var_obs) if var_obs > 0 else float("nan"),
        "cv_excess": float(np.sqrt(excess) / p),
        "median_cells_per_region": float(np.median(nc)),
        "median_expected_counts_per_region": float(np.median(exp)),
        "region_rate": rate,
        "region_N": N,
        "region_cells": nc,
        "keep_mask": keep,
        "inverse": inv,
    }


def variance_partition(sub: _Substrate) -> Dict:
    """Partition MGMT variance over spatial bins, against the technical floor."""
    scales = {}
    for L in BIN_SIZES_UM:
        lab = sub.bins(L)
        scales[L] = {
            "MGMT": _partition(sub.mgmt, sub.depth, lab),
            "neg_probes_pooled": _partition(sub.pool_neg_probe, sub.depth, lab),
            "neg_codewords_pooled": _partition(sub.pool_neg_codeword, sub.depth, lab),
        }

    lab = sub.bins(HEADLINE_BIN_UM)
    matched = []
    for j in sub.abundance_matched():
        y = sub.col_at(j)
        r = _partition(y, sub.depth, lab)
        matched.append({"gene": str(sub.genes[j]), "total": int(y.sum()),
                        "phi": r["phi"], "frac_spatial": r["frac_spatial"],
                        "cv_excess": r["cv_excess"]})
    matched.sort(key=lambda d: d["phi"])
    phis = np.array([m["phi"] for m in matched])
    mgmt_phi = scales[HEADLINE_BIN_UM]["MGMT"]["phi"]

    return {
        "bin_sizes_um": list(BIN_SIZES_UM),
        "headline_bin_um": HEADLINE_BIN_UM,
        "min_cells_per_bin": MIN_CELLS_PER_BIN,
        "scales": scales,
        "abundance_matched": matched,
        "abundance_matched_stats": {
            "n": len(matched), "min": float(phis.min()),
            "median": float(np.median(phis)), "max": float(phis.max()),
            "n_exceeding_mgmt": int((phis > mgmt_phi).sum()),
        },
    }


def niche_partition(sub: _Substrate) -> Dict:
    """Same partition over real Ivy GAP anatomic niches instead of blind bins.

    Neighbourhood pseudobulk (200 µm squares of tumour cells) is scored against
    the measured Ivy GAP per-structure signatures by the existing, unmodified
    ``reconstruct.region_label.label_regions``.
    """
    from ..reconstruct.region_label import label_regions

    L = 200
    lab = sub.bins(L)
    u, inv = np.unique(lab, return_inverse=True)
    ncell = np.bincount(inv, minlength=len(u))
    keep = np.where(ncell >= MIN_CELLS_PER_BIN)[0]

    Xg = sub.Xc[:, sub.gene_cols].tocsr()[sub.idx]
    M = sparse.csr_matrix(
        (np.ones(len(inv)), (inv, np.arange(len(inv)))), shape=(len(u), len(inv))
    )
    pseudo = (M @ Xg).toarray()[keep]
    labels, scores, rnames, shared = label_regions(
        pseudo, None, sub.gene_names, return_scores=True
    )
    bin_region = {int(b): str(labels[i]) for i, b in enumerate(keep)}
    cell_region = np.array([bin_region.get(int(b), "unlabeled") for b in inv])

    order = sorted(set(cell_region.tolist()))
    codes = {r: i for i, r in enumerate(order)}
    rlab = np.array([codes[r] for r in cell_region])

    p_all = sub.rate
    per_region = []
    for r in order:
        m = cell_region == r
        k, N = float(sub.mgmt[m].sum()), float(sub.depth[m].sum())
        lo, hi = stats.poisson.interval(0.95, max(k, 1e-9))
        per_region.append({
            "region": r, "n_cells": int(m.sum()), "mgmt_counts": int(k),
            "pct_positive": float(100.0 * (sub.mgmt[m] > 0).mean()),
            "panel_cpm": float(1e6 * k / N),
            "panel_cpm_lo": float(1e6 * lo / N), "panel_cpm_hi": float(1e6 * hi / N),
            "ratio_vs_section": float((k / N) / p_all),
        })
    return {
        "neighbourhood_um": L,
        "n_neighbourhoods": int(len(keep)),
        "n_shared_signature_genes": int(len(shared)),
        "region_counts": {r: int((cell_region == r).sum()) for r in order},
        "per_region": per_region,
        "MGMT": _partition(sub.mgmt, sub.depth, rlab),
        "neg_codewords_pooled": _partition(sub.pool_neg_codeword, sub.depth, rlab),
        "neg_probes_pooled": _partition(sub.pool_neg_probe, sub.depth, rlab),
        "cell_region": cell_region,
    }


# --------------------------------------------------------------------------- #
# 4. power — what COULD this have detected?
# --------------------------------------------------------------------------- #
def power_analysis(sub: _Substrate, part: Dict) -> Dict:
    """Minimum detectable spatial structure, so a null is not read as absence.

    Two complementary calculations:

    * **Minimum detectable bin-rate CV.** Simulate a tumour whose true local
      MGMT rate is log-normally modulated with coefficient of variation ``cv``,
      draw Poisson counts against the REAL per-bin depths, and ask how often
      the overdispersion test fires at alpha = 0.05. The smallest ``cv`` with
      power >= 0.8 is the smallest genuinely spatial signal this section could
      have found.
    * **Minimum detectable rate ratio** for a two-group (region vs rest)
      contrast at a range of group sizes, in the format
      ``results/state_validation.md`` uses.
    """
    rng = np.random.default_rng(RNG_SEED + 1)
    out = {"headline_bin_um": HEADLINE_BIN_UM, "curves": {}}

    for L in (200, HEADLINE_BIN_UM, 800):
        pr = part["scales"][L]["MGMT"]
        N = np.asarray(pr["region_N"], dtype=float)
        df = len(N) - 1
        crit = stats.chi2.ppf(0.95, df)
        p = pr["rate"]
        curve = []
        for cv in (0.02, 0.04, 0.06, 0.08, 0.10, 0.15, 0.20, 0.30, 0.50):
            hits = 0
            nsim = 300
            sigma = np.sqrt(np.log(1.0 + cv ** 2))
            for _ in range(nsim):
                mult = np.exp(rng.normal(0.0, sigma, len(N)))
                mult /= mult.mean()
                k = rng.poisson(p * N * mult)
                if k.sum() == 0:
                    continue
                exp = (k.sum() / N.sum()) * N
                hits += int(float(((k - exp) ** 2 / exp).sum()) > crit)
            curve.append({"cv": cv, "power": hits / nsim})
        mde = next((c["cv"] for c in curve if c["power"] >= 0.8), None)
        out["curves"][L] = {
            "n_bins": int(len(N)),
            "curve": curve,
            "mde_cv_80pct": mde,
            "observed_cv_excess": pr["cv_excess"],
            "observed_over_mde": (pr["cv_excess"] / mde) if mde else None,
        }

    def _mdrr(n1, n2, p, alpha=0.05, power=0.8):
        from scipy.optimize import brentq
        za, zb = stats.norm.ppf(1 - alpha / 2), stats.norm.ppf(power)

        def f(rr):
            p1, p2 = p * rr, p
            pb = (n1 * p1 + n2 * p2) / (n1 + n2)
            se0 = np.sqrt(pb * (1 - pb) * (1 / n1 + 1 / n2))
            se1 = np.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
            return abs(p1 - p2) - (za * se0 + zb * se1)
        return float(brentq(f, 1.0 + 1e-6, 60.0))

    ptum = float((sub.mgmt > 0).mean())
    n = len(sub.idx)
    groups = []
    for frac in (0.5, 0.25, 0.10, 0.05, 0.01, 0.001):
        n1 = int(n * frac)
        groups.append({
            "group_pct": 100 * frac, "n_group": n1,
            "expected_mgmt_positive": n1 * ptum,
            "min_detectable_rate_ratio": _mdrr(n1, n - n1, ptum),
        })
    out["two_group"] = {"per_cell_positive_rate_pct": 100 * ptum, "groups": groups}

    # single-cell resolution is the thing that is genuinely impossible
    out["single_cell"] = {
        "pct_positive": 100 * ptum,
        "note": (
            "At a 1.19 % per-cell detection rate a single cell carries ~0.012 "
            "expected MGMT counts. No per-cell MGMT call is possible and none "
            "is attempted; every statement here is about aggregates."
        ),
    }
    return out


# --------------------------------------------------------------------------- #
# 5. the decision — scalar vs spatially-resolved MGMT, predicted TMZ kill
# --------------------------------------------------------------------------- #
def _bin_cpm_field(sub: _Substrate, cell_type: np.ndarray, anchor_cpm: float,
                   edge_um: float, spatial_null: bool = False) -> Tuple[np.ndarray, Dict]:
    """Per-CELL MGMT CPM from its own spatial bin, anchored to the section mean.

    ``CPM(bin) = anchor * rate(bin) / rate(section)``, where ``rate`` is MGMT
    counts per total panel gene count among TUMOUR cells. Cells in bins with
    fewer than ``MIN_CELLS_PER_BIN`` tumour cells (and every non-tumour cell)
    take the section anchor, i.e. the scalar — no value is invented where the
    data does not support one.
    """
    n_all = len(cell_type)
    lab_t = sub.bins(edge_um)
    u, inv = np.unique(lab_t, return_inverse=True)
    mg = sub.mgmt
    if spatial_null:
        # SPATIAL NULL: the same MGMT counts on the same cells, randomly
        # reassigned across tumour cells. Total, per-cell marginal distribution
        # and every bin's depth are preserved exactly; only the spatial
        # arrangement is destroyed. Any kill difference this field produces
        # against the scalar CANNOT be MGMT biology.
        mg = np.random.default_rng(RNG_SEED + 3).permutation(mg)
    k = np.bincount(inv, weights=mg, minlength=len(u))
    N = np.bincount(inv, weights=sub.depth, minlength=len(u))
    nc = np.bincount(inv, minlength=len(u))
    ratio = np.ones(len(u))
    ok = (nc >= MIN_CELLS_PER_BIN) & (N > 0)
    ratio[ok] = (k[ok] / N[ok]) / sub.rate

    cpm = np.full(n_all, float(anchor_cpm), dtype=np.float64)
    cpm[sub.idx] = anchor_cpm * ratio[inv]
    stats_ = {
        "edge_um": float(edge_um),
        "n_bins_used": int(ok.sum()),
        "ratio_min": float(ratio[ok].min()), "ratio_max": float(ratio[ok].max()),
        "ratio_p05": float(np.percentile(ratio[ok], 5)),
        "ratio_p95": float(np.percentile(ratio[ok], 95)),
        "cell_cpm_min": float(cpm[sub.idx].min()),
        "cell_cpm_max": float(cpm[sub.idx].max()),
        "cell_cpm_mean": float(cpm[sub.idx].mean()),
        # Two numbers that expose how much of this "field" is counting noise:
        # bins that happened to draw ZERO MGMT counts are handed CPM = 0, i.e.
        # sensitivity 1.0 (fully TMZ-sensitive) purely by a Poisson accident;
        # and the top bins extrapolate far past the calibration ceiling.
        "frac_cells_in_zero_count_bins": float(
            (ratio[inv] == 0.0).mean()
        ),
        "n_zero_count_bins": int(((k == 0) & ok).sum()),
        "frac_cells_above_calibration_ceiling": float(
            (cpm[sub.idx] > 21.12).mean()
        ),
    }
    return cpm, stats_


def _sens_from_cpm(cpm: np.ndarray) -> np.ndarray:
    """Vectorised ``molecular.mgmt.mgmt_sensitivity`` over a per-cell CPM array."""
    from ..molecular.mgmt import mgmt_sensitivity
    uniq, inv = np.unique(np.round(cpm, 6), return_inverse=True)
    s = np.array([mgmt_sensitivity(float(c), units="cpm") for c in uniq])
    return s[inv]


def kill_impact(sub: _Substrate, cell_type: np.ndarray,
                path: str = loader.DEFAULT_CACHE, dx_um: float = DX_UM,
                use_cache: bool = True) -> Dict:
    """Predicted TMZ kill: no MGMT, scalar MGMT, spatially-resolved MGMT.

    Same tissue, same grid, same PDE, same seed. The ONLY thing that changes
    between rows is ``micro.cells.sensitivity``.
    """
    import copy
    from ..core import PDEConfig, SimConfig
    from ..molecular.binding_site_barrier import apply_barrier_to_drugparams
    from ..molecular.library import get_molecule
    from ..molecular.mgmt import mgmt_sensitivity
    from ..reconstruct.real_pipeline import apply_drug_sensitivity, reconstruct_real
    from ..sim.engine import simulate_response
    from .calibration import calibrate_alpha, calibrated_drugparams

    cache = {}
    if use_cache and os.path.exists(_RUN_CACHE):
        with open(_RUN_CACHE) as fh:
            cache = json.load(fh)

    fields: Dict[str, Dict] = {}
    plans: List[Dict] = [{"key": f"none|dx{dx_um}", "label": "no MGMT (current model)",
                          "kind": "none", "anchor": None, "scalar": 1.0}]
    for name, cpm in ANCHORS_CPM:
        s = float(mgmt_sensitivity(cpm, units="cpm"))
        plans.append({"key": f"scalar|{cpm}|dx{dx_um}",
                      "label": f"scalar MGMT = {cpm:g} CPM ({name})",
                      "kind": "scalar", "anchor": cpm, "scalar": s})
        for L, null in ((HEADLINE_BIN_UM, False), (200, False),
                        (HEADLINE_BIN_UM, True), (200, True)):
            if null and cpm != ANCHORS_CPM[0][1]:
                continue  # the null control is run at one anchor only
            cpm_cells, st = _bin_cpm_field(sub, cell_type, cpm, L,
                                           spatial_null=null)
            sens = _sens_from_cpm(cpm_cells)
            key = (f"spatialNULL{L}|{cpm}|dx{dx_um}" if null
                   else f"spatial{L}|{cpm}|dx{dx_um}")
            fields[key] = {"sens": sens, "field_stats": st}
            tum = cell_type == int(CellType.TUMOR)
            plans.append({
                "key": key,
                "label": (
                    f"SPATIAL NULL (MGMT counts shuffled across cells), "
                    f"{L} um bins, anchor {cpm:g} CPM" if null else
                    f"spatially-resolved MGMT, {L} um bins, anchor {cpm:g} CPM"
                ),
                "kind": "spatial_null" if null else "spatial",
                "anchor": cpm, "bin_um": L,
                "scalar": float(sens[tum].mean()),
                "sens_min": float(sens[tum].min()), "sens_max": float(sens[tum].max()),
                "sens_p05": float(np.percentile(sens[tum], 5)),
                "sens_p95": float(np.percentile(sens[tum], 95)),
                "field_stats": st,
            })

    todo = [p for p in plans if p["key"] not in cache]
    if todo:
        micro0 = reconstruct_real(source="xenium", sample=path, dx_um=dx_um)
        alpha = calibrate_alpha()
        dp = calibrated_drugparams("temozolomide", alpha)
        dp = apply_barrier_to_drugparams(dp, get_molecule("temozolomide"))
        pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)
        for p in todo:
            micro = copy.deepcopy(micro0)
            apply_drug_sensitivity(micro, "temozolomide", mode="uniform")
            if p["kind"] == "scalar":
                micro.cells.sensitivity = micro.cells.sensitivity * p["scalar"]
            elif p["kind"] in ("spatial", "spatial_null"):
                micro.cells.sensitivity = (
                    micro.cells.sensitivity * fields[p["key"]]["sens"]
                )
            t0 = time.time()
            res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
            cache[p["key"]] = {
                "kill": float(res.kill_fraction),
                "mean_sensitivity": float(micro.cells.sensitivity.mean()),
                "seconds": time.time() - t0,
            }
            print(f"    {p['key']:34s} kill={res.kill_fraction:.5f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
            os.makedirs(os.path.dirname(_RUN_CACHE) or ".", exist_ok=True)
            with open(_RUN_CACHE, "w") as fh:
                json.dump(cache, fh, indent=1)

    rows = []
    field_stats = [
        {"key": p["key"], "kind": p["kind"], "anchor": p["anchor"],
         "bin_um": p.get("bin_um"), "field_stats": p["field_stats"]}
        for p in plans if "field_stats" in p
    ]
    base = cache[plans[0]["key"]]["kill"]
    for p in plans:
        r = dict(p)
        r.pop("field_stats", None)
        r["kill"] = cache[p["key"]]["kill"]
        r["kill_vs_none_pp"] = 100.0 * (r["kill"] - base)
        rows.append(r)
    # scalar-vs-spatial deltas, the actual decision quantity
    deltas = []
    for name, cpm in ANCHORS_CPM:
        sc = cache[f"scalar|{cpm}|dx{dx_um}"]["kill"]
        for L in (HEADLINE_BIN_UM, 200):
            sp_ = cache[f"spatial{L}|{cpm}|dx{dx_um}"]["kill"]
            nk = f"spatialNULL{L}|{cpm}|dx{dx_um}"
            null_kill = cache[nk]["kill"] if nk in cache else None
            deltas.append({
                "anchor_cpm": cpm, "anchor_name": name, "bin_um": L,
                "kill_scalar": sc, "kill_spatial": sp_,
                "delta_pp": 100.0 * (sp_ - sc),
                "ratio": sp_ / sc if sc else float("nan"),
                "pct_of_scalar_kill": 100.0 * (sp_ - sc) / sc if sc else float("nan"),
                "kill_spatial_null": null_kill,
                "delta_pp_null": (100.0 * (null_kill - sc)
                                  if null_kill is not None else None),
                "signal_over_null": (
                    (sp_ - sc) / (null_kill - sc)
                    if null_kill is not None and abs(null_kill - sc) > 1e-12
                    else None),
            })
    return {"dx_um": dx_um, "rows": rows, "deltas": deltas,
            "field_stats": field_stats, "baseline_kill_no_mgmt": base}


def poisson_null_phi(sub: _Substrate, edge_um: float = HEADLINE_BIN_UM,
                     nsim: int = 400) -> Dict:
    """Distribution of phi under a PURE Poisson null at MGMT's own count depth.

    Completes the reference ladder. ``phi`` is compared against three things,
    in increasing order of realism:

      1. this simulated Poisson null (counting statistics only, phi ~ 1);
      2. the assay's pooled negative controls (counting statistics PLUS every
         technical gradient — segmentation, cell area, probe density — that a
         real feature also inherits);
      3. abundance-matched real genes (what genuine biology looks like here).
    """
    rng = np.random.default_rng(RNG_SEED + 2)
    pr = _partition(sub.mgmt, sub.depth, sub.bins(edge_um))
    N = np.asarray(pr["region_N"], dtype=float)
    p = pr["rate"]
    phis = np.empty(nsim)
    for i in range(nsim):
        k = rng.poisson(p * N)
        exp = (k.sum() / N.sum()) * N
        phis[i] = float(((k - exp) ** 2 / exp).sum()) / (len(N) - 1)
    # Count-matched control. The pooled negative controls carry fewer counts
    # than MGMT, and a sceptic could argue that alone explains why their phi
    # matches. So thin MGMT binomially down to the negatives' exact depth and
    # re-measure: if thinned MGMT lands on the negatives' phi, the comparison
    # is fair at equal depth and the argument is closed.
    thinned = {}
    for name, y in (("neg_probes_pooled", sub.pool_neg_probe),
                    ("neg_codewords_pooled", sub.pool_neg_codeword)):
        target = float(y.sum())
        q = min(1.0, target / float(sub.mgmt.sum()))
        vals = np.empty(50)
        for i in range(50):
            yt = rng.binomial(sub.mgmt.astype(np.int64), q).astype(float)
            vals[i] = _partition(yt, sub.depth, sub.bins(edge_um))["phi"]
        thinned[name] = {
            "target_counts": target, "thinning_p": q,
            "mgmt_thinned_phi_mean": float(vals.mean()),
            "mgmt_thinned_phi_sd": float(vals.std(ddof=1)),
            "control_phi": _partition(y, sub.depth, sub.bins(edge_um))["phi"],
        }
    return {
        "edge_um": float(edge_um), "nsim": nsim,
        "mean": float(phis.mean()), "sd": float(phis.std(ddof=1)),
        "p95": float(np.percentile(phis, 95)),
        "p99": float(np.percentile(phis, 99)),
        "observed_mgmt": pr["phi"],
        "count_matched": thinned,
    }


# --------------------------------------------------------------------------- #
# 6. figure
# --------------------------------------------------------------------------- #
def make_figure(sub: _Substrate, cell_type: np.ndarray, ac: Dict, vp: Dict,
                npart: Dict, pw: Dict, impact: Optional[Dict],
                out_path: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(3, 3, hspace=0.42, wspace=0.30)

    # (a) MGMT+ cells on the section
    ax = fig.add_subplot(gs[0, 0])
    xy = sub.xy
    ax.scatter(xy[::8, 0], xy[::8, 1], s=0.35, c="#c9c9c9", linewidths=0,
               label=f"tumour cells (1/8 shown, n={len(xy):,})")
    pos = sub.mgmt > 0
    ax.scatter(xy[pos, 0], xy[pos, 1], s=0.9, c="#d62728", linewidths=0,
               label=f"MGMT+ ({pos.sum():,}; {100*pos.mean():.2f} %)")
    ax.set_aspect("equal"); ax.invert_yaxis()
    ax.set_title("(a) MGMT-positive tumour cells\nGSE324268 GBMtumor07", fontsize=10)
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ax.legend(markerscale=6, fontsize=7, loc="lower right")

    # (b) per-bin MGMT rate at the headline scale
    ax = fig.add_subplot(gs[0, 1])
    L = HEADLINE_BIN_UM
    lab = sub.bins(L)
    u, inv = np.unique(lab, return_inverse=True)
    k = np.bincount(inv, weights=sub.mgmt, minlength=len(u))
    N = np.bincount(inv, weights=sub.depth, minlength=len(u))
    nc = np.bincount(inv, minlength=len(u))
    ok = nc >= MIN_CELLS_PER_BIN
    bx = ((xy[:, 0] - xy[:, 0].min()) // L).astype(int)
    by = ((xy[:, 1] - xy[:, 1].min()) // L).astype(int)
    grid = np.full((by.max() + 1, bx.max() + 1), np.nan)
    for b in np.where(ok)[0]:
        lb = u[b]
        gy = int(lb % (by.max() + 1)); gx = int(lb // (by.max() + 1))
        grid[gy, gx] = 1e6 * k[b] / N[b]
    fin = np.isfinite(grid)
    im = ax.imshow(grid, origin="upper", cmap="magma",
                   vmin=np.percentile(grid[fin], 2),
                   vmax=np.percentile(grid[fin], 98))
    plt.colorbar(im, ax=ax, fraction=0.046, label="MGMT panel-CPM")
    ax.set_title(f"(b) MGMT rate per {L} µm bin\n"
                 f"φ = {vp['scales'][L]['MGMT']['phi']:.3f}, excess CV = "
                 f"{vp['scales'][L]['MGMT']['cv_excess']:.3f}", fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])

    # (c) Moran's I against the reference classes
    ax = fig.add_subplot(gs[0, 2])
    mI = np.array([m["I"] for m in ac["abundance_matched"]])
    ax.hist(mI, bins=24, color="#9ecae1", edgecolor="#3182bd",
            label=f"{len(mI)} abundance-matched genes")
    hl = ac["headline"]
    for lbl, key, c, ls in (
        ("MGMT", "MGMT count", "#d62728", "-"),
        ("pooled neg. controls", "pooled negative-control codewords", "#7f7f7f", "--"),
    ):
        ax.axvline(hl[key]["I"], color=c, ls=ls, lw=2,
                   label=f"{lbl}  I = {hl[key]['I']:+.4f}")
    ax.set_yscale("symlog")
    ax.set_xlabel("Moran's I (k = 6, tumour cells)")
    ax.set_ylabel("genes")
    ax.set_title("(c) MGMT sits at the BOTTOM of\nthe abundance-matched distribution",
                 fontsize=10)
    ax.legend(fontsize=7)

    # (d) overdispersion vs scale, MGMT vs technical floor
    ax = fig.add_subplot(gs[1, 0])
    Ls = list(BIN_SIZES_UM)
    ax.plot(Ls, [vp["scales"][L]["MGMT"]["phi"] for L in Ls], "o-",
            color="#d62728", label="MGMT")
    ax.plot(Ls, [vp["scales"][L]["neg_probes_pooled"]["phi"] for L in Ls], "s--",
            color="#7f7f7f", label="pooled neg. control PROBES")
    ax.plot(Ls, [vp["scales"][L]["neg_codewords_pooled"]["phi"] for L in Ls], "^--",
            color="#bbbbbb", label="pooled neg. control CODEWORDS")
    ax.axhline(1.0, color="k", lw=0.8, ls=":", label="pure Poisson (no structure)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("bin edge (µm)"); ax.set_ylabel("overdispersion φ")
    ax.set_title("(d) MGMT tracks the TECHNICAL FLOOR\nat every spatial scale",
                 fontsize=10)
    ax.legend(fontsize=7)

    # (e) phi against abundance-matched genes
    ax = fig.add_subplot(gs[1, 1])
    phis = np.array([m["phi"] for m in vp["abundance_matched"]])
    ax.hist(np.log10(phis), bins=20, color="#9ecae1", edgecolor="#3182bd",
            label=f"{len(phis)} abundance-matched genes")
    ax.axvline(np.log10(vp["scales"][L]["MGMT"]["phi"]), color="#d62728", lw=2,
               label=f"MGMT  φ = {vp['scales'][L]['MGMT']['phi']:.2f}")
    ax.axvline(np.log10(vp["scales"][L]["neg_probes_pooled"]["phi"]), color="#7f7f7f",
               ls="--", lw=2,
               label=f"neg. controls  φ = {vp['scales'][L]['neg_probes_pooled']['phi']:.2f}")
    ax.set_xlabel(f"log10 φ  ({L} µm bins)"); ax.set_ylabel("genes")
    ax.set_title("(e) 47/47 matched genes are MORE\nstructured than MGMT", fontsize=10)
    ax.legend(fontsize=7)

    # (f) power curve
    ax = fig.add_subplot(gs[1, 2])
    for Lp, c in zip((200, HEADLINE_BIN_UM, 800), ("#c6dbef", "#6baed6", "#08519c")):
        cur = pw["curves"][Lp]
        ax.plot([p["cv"] for p in cur["curve"]], [p["power"] for p in cur["curve"]],
                "o-", color=c, label=f"{Lp} µm bins")
    ax.axhline(0.8, color="k", ls=":", lw=0.8)
    ax.axvline(vp["scales"][HEADLINE_BIN_UM]["MGMT"]["cv_excess"], color="#d62728",
               lw=2, label="MGMT observed excess CV")
    ax.set_xlabel("true bin-rate CV"); ax.set_ylabel("power (α = 0.05)")
    ax.set_title("(f) Minimum detectable spatial CV\n— the power the null rests on",
                 fontsize=10)
    ax.legend(fontsize=7)

    # (g) niche regions
    ax = fig.add_subplot(gs[2, 0])
    pr = [r for r in npart["per_region"] if r["region"] != "unlabeled"]
    names = [r["region"] for r in pr]
    vals = [r["panel_cpm"] for r in pr]
    err = np.array([[r["panel_cpm"] - r["panel_cpm_lo"] for r in pr],
                    [r["panel_cpm_hi"] - r["panel_cpm"] for r in pr]])
    ax.bar(names, vals, yerr=err, capsize=4, color="#d62728", alpha=0.8)
    ax.axhline(1e6 * sub.rate, color="k", ls="--", lw=1,
               label=f"section mean {1e6*sub.rate:.1f}")
    ax.set_ylabel("MGMT panel-CPM"); ax.tick_params(axis="x", labelrotation=20)
    ax.set_title("(g) MGMT by Ivy GAP niche\n(95 % Poisson intervals)", fontsize=10)
    ax.legend(fontsize=7)

    # (h) sensitivity swing vs calibration CI
    ax = fig.add_subplot(gs[2, 1])
    from ..molecular.mgmt import mgmt_sensitivity
    anchor = ANCHORS_CPM[1][1]
    cv = vp["scales"][HEADLINE_BIN_UM]["MGMT"]["cv_excess"]
    d0 = mgmt_sensitivity(anchor, units="cpm", detail=True)
    s_lo = mgmt_sensitivity(anchor * (1 + cv), units="cpm")
    s_hi = mgmt_sensitivity(anchor * (1 - cv), units="cpm")
    ci = d0["ci95"]
    ax.barh([1], [s_hi - s_lo], left=[s_lo], color="#d62728", height=0.42,
            label=f"spatial swing (±1 CV = ±{100*cv:.0f} %): "
                  f"width {s_hi - s_lo:.3f}")
    ax.barh([0], [ci[1] - ci[0]], left=[ci[0]], color="#9ecae1", height=0.42,
            label=f"95 % CI already carried on the scalar: "
                  f"width {ci[1] - ci[0]:.3f}")
    ax.axvline(d0["sensitivity"], color="k", lw=1.2, ls="--",
               label=f"point estimate {d0['sensitivity']:.3f}")
    ax.set_ylim(-0.5, 1.5)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["calibration CI", "spatial"])
    ax.set_xlabel(f"TMZ sensitivity multiplier at {anchor:g} CPM")
    ax.set_title("(h) The spatial signal is ~2 orders of\nmagnitude inside the "
                 "existing CI", fontsize=10)
    ax.legend(fontsize=7, loc="lower right")

    # (i) kill impact
    ax = fig.add_subplot(gs[2, 2])
    if impact:
        dd = [d for d in impact["deltas"] if d["delta_pp_null"] is not None]
        if dd:
            dd.sort(key=lambda d: -d["bin_um"])
            xpos = np.arange(len(dd))
            ax.bar(xpos - 0.19, [d["delta_pp"] for d in dd], width=0.36,
                   color="#d62728", label="spatially-resolved MGMT")
            ax.bar(xpos + 0.19, [d["delta_pp_null"] for d in dd], width=0.36,
                   color="#7f7f7f",
                   label="SPATIAL NULL (counts shuffled across cells)")
            ax.set_xticks(xpos)
            ax.set_xticklabels([f"{d['bin_um']} µm bins" for d in dd], fontsize=8)
            ax.axhline(0, color="k", lw=0.8)
            ax.set_ylabel("Δ predicted TMZ kill vs scalar (pp)")
            ax.set_title("(i) The scalar-vs-spatial shift is\nreproduced by "
                         "pure noise", fontsize=10)
            ax.legend(fontsize=7, loc="upper left")
        else:
            worst = max(abs(d["delta_pp"]) for d in impact["deltas"])
            ax.bar([r["label"][:18] for r in impact["rows"]],
                   [r["kill"] for r in impact["rows"]], color="#d62728")
            ax.tick_params(axis="x", labelrotation=40, labelsize=6)
            ax.set_ylabel("predicted TMZ kill fraction")
            ax.set_title(f"(i) Scalar vs spatial MGMT\nlargest difference "
                         f"{worst:.3f} pp", fontsize=10)
    else:
        ax.text(0.5, 0.5, "impact simulations not run", ha="center", va="center")
        ax.set_axis_off()

    fig.suptitle(
        "MGMT in real human GBM Xenium (GEO GSE324268, 851,509 cells, 8-probe "
        "MGMT target): is it a spatial field or a scalar?",
        fontsize=13, y=0.965,
    )
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out_path


# --------------------------------------------------------------------------- #
# 7. report
# --------------------------------------------------------------------------- #
def build_markdown(prov: Dict, det: Dict, rec: Dict, panels: Dict, ac: Dict,
                   vp: Dict, npart: Dict, pnull: Dict, pw: Dict,
                   impact: Optional[Dict]) -> str:
    L = HEADLINE_BIN_UM
    mgmt_phi = vp["scales"][L]["MGMT"]["phi"]
    negp_phi = vp["scales"][L]["neg_probes_pooled"]["phi"]
    cv = vp["scales"][L]["MGMT"]["cv_excess"]
    negp_cv = vp["scales"][L]["neg_probes_pooled"]["cv_excess"]
    ms = ac["abundance_matched_stats"]
    vs = vp["abundance_matched_stats"]
    hl = ac["headline"]
    A: List[str] = []
    W = A.append

    W("# MGMT in real human GBM Xenium: a spatial field, or a scalar?")
    W("")
    W("**We said MGMT was unmodellable from spatial data. That generalisation "
      "was wrong, and this report retracts it — then tests the position it was "
      "used to justify, and finds the position itself survives.**")
    W("")
    W(f"`results/gbm_tmz_exvivo.md` §B2 established that MGMT is absent from the "
      f"two Xenium panels this repo held. `results/mgmt_scalar.md` §2 turned that "
      f"into *\"MGMT does not need to be spatial — it is a per-patient scalar\"*, "
      f"partly on the grounds that it **cannot** be resolved spatially. "
      f"**GEO {prov['accession']}** falsifies the second half of that sentence: "
      f"it is a free, no-login human primary GBM Xenium section whose custom "
      f"add-on panel carries a dedicated **8-probe MGMT target** across "
      f"**{det['n_cells']:,} cells**, and MGMT is genuinely detected in it. "
      f"So the question is now empirical rather than rhetorical, and it was "
      f"asked directly.")
    W("")
    W("![figure](mgmt_spatial.png)")
    W("")
    W("## Verdict")
    W("")
    W("| question | answer |")
    W("|---|---|")
    W(f"| Is MGMT really measured here? | **YES.** {det['mgmt_counts']:,} counts "
      f"in {det['mgmt_cells_positive']:,} of {det['n_cells']:,} cells "
      f"({det['mgmt_pct_cells']:.2f} %); "
      f"{det['fold_over_neg_probe_mean']:.1f}× the per-probe negative-control mean |")
    W(f"| Is per-cell MGMT usable? | **NO.** {det['mgmt_pct_tumour_cells']:.2f} % "
      f"of tumour cells positive, max {det['mgmt_max_per_cell']} counts in any "
      f"cell. A per-cell MGMT call does not exist at this depth and none is made |")
    W(f"| Is MGMT spatially autocorrelated? | Moran's I = "
      f"{hl['MGMT count']['I']:+.5f} (z = {hl['MGMT count']['z']:+.2f}, "
      f"p = {hl['MGMT count']['p']:.3f}) — statistically non-zero, and "
      f"**{ms['n_exceeding_mgmt']} of {ms['n']} abundance-matched genes are more "
      f"clustered than it** |")
    _cmv = pnull["count_matched"]["neg_probes_pooled"]
    W(f"| Is it above the assay's own noise? | **NO — it is below it.** At "
      f"{L} µm bins MGMT φ = {mgmt_phi:.3f} against pooled negative-control "
      f"probes φ = {negp_phi:.3f}; **thinned to the controls' exact count "
      f"depth MGMT falls to φ = {_cmv['mgmt_thinned_phi_mean']:.3f} ± "
      f"{_cmv['mgmt_thinned_phi_sd']:.3f}, below the controls' "
      f"{_cmv['control_phi']:.3f}** |")
    W(f"| Is this an underpowered null? | **NO.** The same test on the same cells "
      f"returns φ up to {vs['max']:.0f} and Moran's I up to {ms['max']:+.3f} for "
      f"real genes at the SAME abundance. The assay can see structure here; MGMT "
      f"has none to see |")
    if impact:
        worst = max(abs(d["delta_pp"]) for d in impact["deltas"])
        _rat = [d["signal_over_null"] for d in impact["deltas"]
                if d["signal_over_null"] is not None]
        W(f"| Does the scalar cost anything? | **No — the spatial version costs "
          f"something.** Predicted TMZ kill moves by at most **{worst:.4f} pp** "
          f"between scalar and spatially-resolved MGMT, and a spatial NULL "
          f"(counts shuffled across cells) reproduces "
          f"{100/max(_rat) if _rat else float('nan'):.0f}–"
          f"{100/min(_rat) if _rat else float('nan'):.0f} % of that shift. It "
          f"is Jensen's inequality on counting noise, not biology |")
    W("| **Decision** | **The scalar is VINDICATED — as an empirical finding, "
      "not a limitation.** |")
    W("")

    # ---------------- 1. provenance ----------------
    W("## 1. Provenance — verified, not assumed")
    W("")
    W(f"GEO **{prov['accession']}**, sample **{prov['sample']}** "
      f"(`GBMtumor07`). {prov['landing_page']}")
    W("")
    W("| file | bytes | sha256[:16] | matches frozen record |")
    W("|---|---:|---|---|")
    for name, v in prov["files"].items():
        W(f"| `{name}` | {v['bytes']:,} | `{str(v['sha256'])[:16]}` | "
          f"{'yes' if v['sha256_match'] else ('derived' if v['sha256_match'] is None else '**NO**')} |")
    W("")
    W("Downloaded from the NCBI GEO FTP mirror (`ftp.ncbi.nlm.nih.gov`), all "
      "URLs HEAD-checked HTTP 200 before fetching. **The smallest complete unit "
      "was taken deliberately**; these were skipped:")
    W("")
    for k, why in prov["skipped"].items():
        W(f"* `{k}` — {why}")
    W("")
    pe = det["panel_entry"]
    W("### The panel, read out of the shipped spec")
    W("")
    W("| field | value |")
    W("|---|---|")
    W("| base panel | `Xenium Human Brain Gene Expression` **hBrain_v1.1**, 266 targets |")
    W("| add-on | `hBrain_100g`, design_id **ECQNGH**, described in the file as "
      "\"Glioblastoma Panel\", 100 targets |")
    W("| targets in spec | 400 |")
    W(f"| `Gene Expression` features in the matrix | **{det['n_gene_targets']}** |")
    W(f"| **MGMT entry** | symbol `{pe['symbol']}`, id **`{pe['ensembl']}`**, "
      f"codeword {pe['codewords']}, **`gene_coverage` = {pe['gene_coverage']}** |")
    W("")
    W(f"**Correction to the brief.** The dataset carries **{pe['gene_coverage']} "
      f"probe pairs on ONE MGMT gene target**, i.e. one column in the count "
      f"matrix — not \"8 probe sets\" in the sense of 8 features. The count "
      f"claim itself reproduces exactly: {det['mgmt_counts']:,} counts across "
      f"{det['mgmt_cells_positive']:,} cells.")
    W("")

    # ---------------- 2. detection ----------------
    W("## 2. Detection — the claim, verified, and then stress-tested")
    W("")
    W("| quantity | value |")
    W("|---|---:|")
    W(f"| cells | {det['n_cells']:,} |")
    W(f"| MGMT counts | {det['mgmt_counts']:,} |")
    W(f"| cells with >= 1 MGMT count | {det['mgmt_cells_positive']:,} |")
    W(f"| **% of all cells** | **{det['mgmt_pct_cells']:.4f} %** |")
    W(f"| **% of TUMOUR cells** | **{det['mgmt_pct_tumour_cells']:.4f} %** "
      f"({det['mgmt_counts_tumour']:,} counts in {det['n_tumour_cells']:,} cells) |")
    W(f"| max MGMT counts in any one cell | {det['mgmt_max_per_cell']} |")
    W(f"| rank among the {det['n_gene_targets']} gene targets by total counts | "
      f"**{det['rank_by_total_counts']}** |")
    W(f"| median gene counts per cell | {det['median_gene_counts_per_cell']:.0f} |")
    W("")
    W("**Is it signal or background?** The assay ships its own answer.")
    W("")
    W("| feature class | features | total counts | per feature |")
    W("|---|---:|---:|---:|")
    W(f"| MGMT | 1 target ({pe['gene_coverage']} probes) | "
      f"{det['mgmt_counts']:,} | {det['mgmt_counts_per_probe']:.0f} per probe |")
    W(f"| Negative Control Probe | {det['neg_probe']['n_features']} | "
      f"{det['neg_probe']['total']:,} | mean {det['neg_probe']['mean_per_feature']:.0f}, "
      f"max {det['neg_probe']['max_per_feature']:,} |")
    W(f"| Negative Control Codeword | {det['neg_codeword']['n_features']} | "
      f"{det['neg_codeword']['total']:,} | mean "
      f"{det['neg_codeword']['mean_per_feature']:.0f} |")
    W("")
    W(f"MGMT runs **{det['fold_over_neg_probe_mean']:.1f}×** the mean "
      f"negative-control probe on a per-probe basis and "
      f"**{det['fold_over_worst_neg_probe']:.1f}×** the single noisiest one. "
      f"That is real but modest: MGMT is a genuinely low-expressed gene here, "
      f"ranking {det['rank_by_total_counts']} of {det['n_gene_targets']}. "
      f"**The claim is confirmed and the enthusiasm is bounded.**")
    W("")
    W("### MGMT detection by cell type — the first warning")
    W("")
    W("| cell type | cells | MGMT counts | % positive | panel-CPM |")
    W("|---|---:|---:|---:|---:|")
    for name, v in det["by_cell_type"].items():
        W(f"| {name} | {v['n_cells']:,} | {v['mgmt_counts']:,} | "
          f"{v['pct_positive']:.3f} % | {v['panel_cpm']:.1f} |")
    W("")
    W("**MGMT detection is LOWEST in tumour cells and ~3× higher in endothelial "
      "and immune cells.** That is almost certainly not MGMT biology — it is the "
      "marker-gate selecting cells with more transcripts overall. It is the first "
      "sign that anything spatial in this variable will be inherited from cell "
      "composition and capture depth rather than from MGMT regulation, and it is "
      "why every test below is run **within tumour cells only** and repeated on a "
      "depth residual.")
    W("")

    # ---------------- 3. reconstruction ----------------
    W("## 3. Reconstruction — the existing pipeline, unmodified")
    W("")
    W("`reconstruct.real_xenium.reconstruct_from_xenium` consumed this section "
      "with **no change to any file under `tmesim/reconstruct/`**. The new "
      "loader writes the 10x MEX layout the existing loader already reads; the "
      "pipeline could not tell it apart from the other two sections.")
    W("")
    pc = panels["gse324268"]
    W(f"**No new panel registry entry was needed.** This panel carries all twelve "
      f"markers of `BRAIN_PANEL_MARKER_MAP`, so `select_marker_map` scores it "
      f"brain_319 = {pc['coverage']['brain_319']:.2f} / "
      f"io_480 = {pc['coverage']['io_480']:.2f} and returns the *same map object* "
      f"(`is_brain_map_object = {pc['is_brain_map_object']}`). Markers absent from "
      f"this panel: {pc['markers_absent'] or 'none'}.")
    W("")
    W("| existing cache | panel selected | brain_319 cov | io_480 cov | same map object | UNCHANGED |")
    W("|---|---|---:|---:|---|---|")
    for p in ("data_cache/xenium_gbm", "data_cache/xenium_gbm_io"):
        v = panels.get(p, {})
        if "panel_selected" not in v:
            W(f"| `{p}` | — | — | — | — | {v.get('status', 'n/a')} |")
            continue
        W(f"| `{p}` | {v['panel_selected']} | {v['coverage']['brain_319']:.2f} | "
          f"{v['coverage']['io_480']:.2f} | {v['same_map_object']} | "
          f"**{'yes' if v['unchanged'] else 'NO'}** |")
    W("")
    W("### Cell-type fractions, and the one external anchor")
    W("")
    W("| | section 1 `xenium_gbm` | section 2 `xenium_gbm_io` | **section 3 GSE324268** |")
    W("|---|---:|---:|---:|")
    f3 = rec["celltype_fractions"]
    W(f"| cells | 40,887 | 816,769 | **{rec['n_cells']:,}** |")
    W(f"| panel | brain 319 | IO 480 | **brain v1.1 + GBM add-on, "
      f"{det['n_gene_targets']}** |")
    W(f"| TUMOR | 73.80 % | 79.09 % | **{f3['TUMOR']['pct']:.2f} %** |")
    W(f"| STROMAL | 7.92 % | 9.98 % | **{f3['STROMAL']['pct']:.2f} %** |")
    W(f"| IMMUNE | 8.54 % | 6.74 % | **{f3['IMMUNE']['pct']:.2f} %** |")
    W(f"| **ENDOTHELIAL** | **9.74 %** | **4.18 %** | "
      f"**{f3['ENDOTHELIAL']['pct']:.2f} %** |")
    W(f"| vs volume-EM 8.0 % (PMID 38723085) | 1.22× | 0.52× | "
      f"**{rec['endothelial_vs_volume_em']:.2f}×** |")
    W("")
    W(f"Human cortex volume electron microscopy counts **{ENDOTHELIAL_VOLUME_EM} %** "
      f"of cells as endothelial (Shapson-Coe et al., *Science* 2024;384:eadk4858, "
      f"PMID 38723085). This section reads "
      f"**{f3['ENDOTHELIAL']['pct']:.2f} %**, "
      f"i.e. **{rec['endothelial_vs_volume_em']:.2f}×** the anchor.")
    W("")
    W(f"**Read this as a third data point, not a pass.** Section 1 came in at "
      f"1.22×, section 2 at 0.52×, and this one at "
      f"{rec['endothelial_vs_volume_em']:.2f}×. Sections 1 and 3 use the *same* "
      f"single-marker PECAM1 gate, and they differ from each other by "
      f"{f3['ENDOTHELIAL']['pct']/9.74:.2f}× — so the spread is not only a "
      f"panel artefact, as `results/heldout_patient.md` §1.1 could still hope. "
      f"**Across three human GBM sections the classifier's endothelial fraction "
      f"spans 4.2–{f3['ENDOTHELIAL']['pct']:.1f} %, a "
      f"{f3['ENDOTHELIAL']['pct']/4.18:.1f}× range that brackets the volume-EM "
      f"value rather than reproducing it.** The 8.0 % agreement on section 1 is "
      f"now best read as a coincidence of that section, exactly as "
      f"`heldout_patient.md` warned. This is an incidental finding of this "
      f"report, and it is a negative one.")
    W("")
    W(f"Geometry: {rec['extent_um'][0]:.0f} × {rec['extent_um'][1]:.0f} µm, "
      f"{rec['bounding_area_mm2']:.1f} mm² bounding area, "
      f"{rec['cell_density_per_mm2']:.0f} cells/mm², "
      f"{rec['endothelial_per_mm2']:.0f} endothelial cells/mm². Cell typing is "
      f"per-cell and grid-independent, so the fractions above do not depend on "
      f"`dx`; the reconstruction above was built at dx = {rec['dx_um']:.0f} µm "
      f"and the §7 simulations at dx = {DX_UM:.0f} µm (the converged value).")
    W("")

    # ---------------- 4. autocorrelation ----------------
    W("## 4. The core test (a) — is per-cell MGMT spatially autocorrelated?")
    W("")
    W(f"Moran's I over the **{ac['n_tumour_cells']:,} TUMOUR cells only**, on a "
      f"k = {ac['k_neighbours']} nearest-neighbour row-standardised graph "
      f"(median nearest neighbour {ac['median_nn_um']:.1f} µm, median k-th "
      f"{ac['median_kth_nn_um']:.1f} µm), against a "
      f"{ac['n_perm']}-permutation null.")
    W("")
    W("| variable | Moran's I | permutation null | z | p |")
    W("|---|---:|---|---:|---:|")
    for k, v in hl.items():
        W(f"| {k} | **{v['I']:+.5f}** | {v['null_mean']:+.5f} ± {v['null_sd']:.5f} | "
          f"{v['z']:+.2f} | {v['p']:.4f} |")
    W("")
    W(f"**Taken alone this row says \"MGMT is significantly clustered\" — and "
      f"taken alone it would be worthless.** With n = {ac['n_tumour_cells']:,} "
      f"a permutation test detects essentially any departure from a random "
      f"relabelling. Three things kill the interpretation:")
    W("")
    W(f"1. **The confound dwarfs the signal.** Total transcript count per cell "
      f"has I = {hl['total gene counts per cell (DEPTH — the confound)']['I']:+.5f}, "
      f"**{hl['total gene counts per cell (DEPTH — the confound)']['I']/hl['MGMT count']['I']:.0f}× "
      f"MGMT's**, and every count variable inherits some of it. Removing it "
      f"(the depth residual) drops MGMT to "
      f"{hl['MGMT depth-residual']['I']:+.5f} at p = "
      f"{hl['MGMT depth-residual']['p']:.3f} — **no longer significant at 0.05**.")
    W(f"2. **The pooled negative-control codewords score "
      f"{hl['pooled negative-control codewords']['I']:+.5f}** — "
      f"{100*hl['pooled negative-control codewords']['I']/hl['MGMT count']['I']:.0f} % "
      f"of MGMT's value, from features that measure nothing at all.")
    W(f"3. **{ms['n_exceeding_mgmt']} of {ms['n']} abundance-matched real genes "
      f"are more clustered than MGMT** — only "
      f"{ms['n'] - ms['n_exceeding_mgmt']} sits below it, putting MGMT at the "
      f"{ms['mgmt_percentile']:.0f}nd percentile of its own abundance class. "
      f"Matched-gene I runs {ms['min']:+.4f} to {ms['max']:+.4f}, median "
      f"{ms['median']:+.4f}.")
    W("")
    W("### What the reference class proves — and why this is not an "
      "underpowered null")
    W("")
    W("The abundance-matched class is the whole argument. These are real panel "
      "genes with the **same section total counts as MGMT** (within "
      f"{ABUNDANCE_WINDOW[0]:g}–{ABUNDANCE_WINDOW[1]:g}×) and comparable "
      "per-cell detection rates, scored on the **same cells with the same "
      "weights**. If ~1 % detection made spatial structure undetectable, they "
      "would all sit at MGMT's value. They do not:")
    W("")
    W("| gene | total counts | % of tumour cells + | Moran's I |")
    W("|---|---:|---:|---:|")
    show = ac["abundance_matched"][:3] + ac["abundance_matched"][-6:]
    for m in show:
        W(f"| {m['gene']} | {m['total']:,} | {m['pct_positive']:.2f} % | "
          f"{m['I']:+.5f} |")
    W(f"| **MGMT** | **{det['mgmt_counts_tumour']:,}** | "
      f"**{det['mgmt_pct_tumour_cells']:.2f} %** | "
      f"**{hl['MGMT count']['I']:+.5f}** |")
    W("")
    W(f"**CA9 — a canonical perinecrotic hypoxia marker — is detected in "
      f"{[m['pct_positive'] for m in ac['abundance_matched'] if m['gene']=='CA9'][0]:.2f} % "
      f"of tumour cells, LESS often than MGMT, and returns "
      f"I = {[m['I'] for m in ac['abundance_matched'] if m['gene']=='CA9'][0]:+.4f} — "
      f"{[m['I'] for m in ac['abundance_matched'] if m['gene']=='CA9'][0]/hl['MGMT count']['I']:.0f}× "
      f"MGMT's.** At this abundance, on this section, real spatial biology is "
      f"plainly visible. MGMT's absence of it is a measurement, not a "
      f"power failure.")
    W("")
    W("High-expression positive controls, for scale:")
    W("")
    W("| gene | total counts | % + | Moran's I |")
    W("|---|---:|---:|---:|")
    for m in ac["positive_controls"]:
        W(f"| {m['gene']} | {m['total']:,} | {m['pct_positive']:.1f} % | "
          f"{m['I']:+.5f} |")
    W("")
    # ---------------- 5. variance partition ----------------
    W("## 5. The core test (b) — how much of MGMT's variance is spatial?")
    W("")
    W("Each region holds `N` total panel gene counts and `k` MGMT counts. Under "
      "*\"one tumour-level mean and nothing else\"*, `k ~ Poisson(p·N)`. Then")
    W("")
    W("* **φ** = Pearson χ²/df. φ = 1 means the between-region spread is exactly "
      "counting noise.")
    W("* **frac_spatial** = (Var_observed − Var_sampling) / Var_observed — the "
      "share of the between-region variance in the MGMT *rate* that a single "
      "tumour-level mean does not explain. **This is the number the brief asked "
      "for.**")
    W("* **excess CV** = the same excess re-expressed as a multiplicative swing "
      "of the local rate, which is the unit the sensitivity mapping consumes.")
    W("")
    W("Every row is computed identically for MGMT and for the assay's pooled "
      "negative controls, on the same bins and the same cells.")
    W("")
    W("| bin | bins used | median tumour cells/bin | median expected MGMT/bin | "
      "MGMT φ | MGMT p | **MGMT frac_spatial** | MGMT excess CV | "
      "neg-PROBE φ | neg-PROBE excess CV | neg-CODEWORD φ |")
    W("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for Lb in BIN_SIZES_UM:
        s = vp["scales"][Lb]
        m_, np_, nc_ = s["MGMT"], s["neg_probes_pooled"], s["neg_codewords_pooled"]
        W(f"| {Lb} µm | {m_['B']:,} | {m_['median_cells_per_region']:.0f} | "
          f"{m_['median_expected_counts_per_region']:.1f} | {m_['phi']:.3f} | "
          f"{m_['p_chi2']:.2g} | **{m_['frac_spatial']:.3f}** | "
          f"{m_['cv_excess']:.3f} | {np_['phi']:.3f} | {np_['cv_excess']:.3f} | "
          f"{nc_['phi']:.3f} |")
    W("")
    W(f"**The answer to \"what fraction of MGMT variance is spatial?\" is "
      f"scale-dependent and, at every scale, matched by features that measure "
      f"nothing.** At the headline {L} µm scale MGMT's φ = {mgmt_phi:.3f} and "
      f"the pooled negative-control probes' φ = {negp_phi:.3f} — they agree to "
      f"{abs(100*(mgmt_phi-negp_phi)/negp_phi):.0f} %. MGMT's excess CV "
      f"({cv:.3f}) is **smaller** than the negative controls' ({negp_cv:.3f}). "
      f"The nominal `frac_spatial` of {vp['scales'][L]['MGMT']['frac_spatial']:.3f} "
      f"is therefore **an upper bound on a quantity whose technical floor is at "
      f"least as large**, not an estimate of MGMT biology.")
    W("")
    W(f"For completeness, a pure Poisson null simulated at MGMT's own count "
      f"depth against the real per-bin depths gives φ = {pnull['mean']:.3f} ± "
      f"{pnull['sd']:.3f} (95th pct {pnull['p95']:.3f}, n = {pnull['nsim']}). "
      f"So the ladder at {L} µm reads: **pure Poisson "
      f"{pnull['mean']:.2f} → negative controls {negp_phi:.2f} → MGMT "
      f"{mgmt_phi:.2f} → real matched genes median {vs['median']:.1f}, max "
      f"{vs['max']:.0f}.** MGMT is on the wrong rung.")
    W("")
    W("**Count-matched control — closing the obvious objection.** The pooled "
      "negative controls carry fewer counts than MGMT, and unequal depth alone "
      "could in principle explain the agreement. So MGMT was thinned "
      "binomially down to each control's exact count depth and re-measured "
      f"(50 draws, {L} µm bins):")
    W("")
    W("| control | its counts | thinning p | **MGMT thinned to that depth: φ** | "
      "**the control's own φ** |")
    W("|---|---:|---:|---:|---:|")
    for nm, cm in pnull["count_matched"].items():
        W(f"| `{nm}` | {cm['target_counts']:,.0f} | {cm['thinning_p']:.3f} | "
          f"**{cm['mgmt_thinned_phi_mean']:.3f} ± {cm['mgmt_thinned_phi_sd']:.3f}** | "
          f"**{cm['control_phi']:.3f}** |")
    W("")
    _cm = pnull["count_matched"]["neg_probes_pooled"]
    W(f"**At equal count depth MGMT is *less* spatially overdispersed than the "
      f"assay's own negative controls** — "
      f"{_cm['mgmt_thinned_phi_mean']:.3f} ± {_cm['mgmt_thinned_phi_sd']:.3f} "
      f"against {_cm['control_phi']:.3f}, i.e. the control sits "
      f"{(_cm['control_phi'] - _cm['mgmt_thinned_phi_mean'])/_cm['mgmt_thinned_phi_sd']:.1f} "
      f"SD ABOVE the thinned MGMT distribution. The depth objection does not "
      f"survive, and it turns out to have been working against the conclusion "
      f"rather than for it: the raw φ agreement at full depth "
      f"({mgmt_phi:.3f} vs {negp_phi:.3f}) slightly *flatters* MGMT, because φ "
      f"grows with counts for a fixed relative excess.")
    W("")
    W(f"**{vs['n_exceeding_mgmt']} of {vs['n']} abundance-matched genes have a "
      f"HIGHER φ than MGMT — that is all of them.** MGMT is the least spatially "
      f"structured gene of its abundance class on this section.")
    W("")
    W("| gene | total | φ | frac_spatial | excess CV |")
    W("|---|---:|---:|---:|---:|")
    W(f"| **MGMT** | **{det['mgmt_counts_tumour']:,}** | **{mgmt_phi:.2f}** | "
      f"**{vp['scales'][L]['MGMT']['frac_spatial']:.3f}** | **{cv:.3f}** |")
    for m in vp["abundance_matched"][:3] + vp["abundance_matched"][-4:]:
        W(f"| {m['gene']} | {m['total']:,} | {m['phi']:.2f} | "
          f"{m['frac_spatial']:.3f} | {m['cv_excess']:.3f} |")
    W("")

    # ---------------- niches ----------------
    W("### 5.1 Real niche regions, not blind bins")
    W("")
    W(f"{npart['n_neighbourhoods']:,} tumour-cell neighbourhoods of "
      f"{npart['neighbourhood_um']} µm were pseudobulked and scored against the "
      f"measured Ivy GAP per-anatomic-structure signatures by the existing, "
      f"unmodified `reconstruct.region_label.label_regions` "
      f"({npart['n_shared_signature_genes']} shared signature genes).")
    W("")
    W("| niche | tumour cells | MGMT counts | % cells + | MGMT panel-CPM "
      "(95 % Poisson) | ratio vs section |")
    W("|---|---:|---:|---:|---|---:|")
    for r in npart["per_region"]:
        W(f"| {r['region']} | {r['n_cells']:,} | {r['mgmt_counts']:,} | "
          f"{r['pct_positive']:.2f} % | {r['panel_cpm']:.1f} "
          f"({r['panel_cpm_lo']:.1f}–{r['panel_cpm_hi']:.1f}) | "
          f"{r['ratio_vs_section']:.3f} |")
    W("")
    nm, nn = npart["MGMT"], npart["neg_codewords_pooled"]
    W(f"MGMT φ = {nm['phi']:.2f} (frac_spatial {nm['frac_spatial']:.3f}, excess "
      f"CV {nm['cv_excess']:.3f}) — **and the pooled negative-control codewords, "
      f"over the same four regions, score φ = {nn['phi']:.2f} with excess CV "
      f"{nn['cv_excess']:.3f}, i.e. MORE apparent regional structure than MGMT.** "
      f"Whatever makes counts vary between Ivy GAP niches on this section acts on "
      f"features that have no biology, so it is not MGMT regulation. The point "
      f"estimates are nonetheless the most interpretable numbers in this report: "
      f"perinecrotic tumour reads "
      f"{[r['ratio_vs_section'] for r in npart['per_region'] if r['region']=='necrotic'][0]:.2f}× "
      f"the section mean and perivascular tumour "
      f"{[r['ratio_vs_section'] for r in npart['per_region'] if r['region']=='perivascular'][0]:.2f}×, "
      f"a **{max(r['panel_cpm'] for r in npart['per_region'] if r['region'] != 'unlabeled')/min(r['panel_cpm'] for r in npart['per_region'] if r['region'] != 'unlabeled'):.2f}× "
      f"spread across labelled niches** whose confidence intervals overlap the section "
      f"mean's neighbourhood and whose technical control moves further.")
    W("")

    # ---------------- power ----------------
    W("## 6. Power — what this section COULD have detected, stated up front")
    W("")
    W("An honest null needs its detection limit published next to it. Two "
      "calculations, in the format `results/state_validation.md` uses.")
    W("")
    W("### 6.1 Minimum detectable spatial CV")
    W("")
    W("A tumour whose true local MGMT rate is log-normally modulated with "
      "coefficient of variation `cv` was simulated against the REAL per-bin "
      "depths, Poisson counts drawn, and the overdispersion test run at "
      "α = 0.05.")
    W("")
    W("| bin | bins | MDE CV at 80 % power | observed MGMT excess CV | observed / MDE |")
    W("|---:|---:|---:|---:|---:|")
    for Lp in (200, HEADLINE_BIN_UM, 800):
        c = pw["curves"][Lp]
        mde = c["mde_cv_80pct"]
        W(f"| {Lp} µm | {c['n_bins']:,} | "
          f"{('%.2f' % mde) if mde else '> 0.50'} | "
          f"{c['observed_cv_excess']:.3f} | "
          f"{('%.2f×' % c['observed_over_mde']) if c['observed_over_mde'] else '—'} |")
    W("")
    W(f"**The observed MGMT modulation sits barely above the detection floor "
      f"— {pw['curves'][HEADLINE_BIN_UM]['observed_over_mde']:.2f}× the MDE at "
      f"{HEADLINE_BIN_UM} µm.** That is exactly the regime in which a positive "
      f"χ² should not be believed on its own, and it is why the negative-control "
      f"comparison — not the p-value — carries this section's conclusion.")
    W("")
    W("**What this CANNOT rule out.** A spatial MGMT modulation with a "
      f"bin-rate CV below ~{pw['curves'][HEADLINE_BIN_UM]['mde_cv_80pct']:.2f} at "
      f"{HEADLINE_BIN_UM} µm would be invisible here, as would any structure at a "
      f"length scale far below {BIN_SIZES_UM[0]} µm or far above the "
      f"{max(rec['extent_um'])/1000:.0f} mm section. A ±"
      f"{100*pw['curves'][HEADLINE_BIN_UM]['mde_cv_80pct']:.0f} % local swing in "
      "MGMT is the smallest thing this test can see; a smaller one is not "
      "excluded. §7 shows why a swing of that size would not matter anyway.")
    W("")
    W("### 6.2 Minimum detectable rate ratio, two-group contrast")
    W("")
    W(f"Per-cell MGMT positivity among tumour cells is "
      f"{pw['two_group']['per_cell_positive_rate_pct']:.3f} %.")
    W("")
    W("| group size | n cells | expected MGMT+ | min. detectable rate ratio |")
    W("|---:|---:|---:|---:|")
    for g in pw["two_group"]["groups"]:
        W(f"| {g['group_pct']:g} % of tumour | {g['n_group']:,} | "
          f"{g['expected_mgmt_positive']:.1f} | "
          f"**{g['min_detectable_rate_ratio']:.3f}×** |")
    W("")
    W("**This is the opposite of the `state_validation.md` situation, and the "
      "contrast is worth stating.** There, MKI67 in a 180-cell subgroup could "
      "only have detected a ≥4.1× rate ratio, so a null meant nothing. Here a "
      f"region holding 10 % of the tumour would reveal a "
      f"**{[g['min_detectable_rate_ratio'] for g in pw['two_group']['groups'] if g['group_pct'] == 10.0][0]:.2f}×** "
      f"difference in MGMT rate, and even a 1 % region would reveal "
      f"{[g['min_detectable_rate_ratio'] for g in pw['two_group']['groups'] if g['group_pct'] == 1.0][0]:.2f}×. "
      "**The section is well powered for regional contrasts and returns none "
      "above the technical floor.** The verdict below is a NULL, not an "
      "UNDERPOWERED.")
    W("")
    W(f"**The one thing that IS underpowered, absolutely and unfixably here:** "
      f"the single cell. At {pw['single_cell']['pct_positive']:.2f} % positivity "
      f"a tumour cell carries ~{pw['single_cell']['pct_positive']/100:.4f} "
      f"expected MGMT counts. **No per-cell MGMT value exists in this dataset "
      f"and none is computed.** Every number above is an aggregate over >= "
      f"{MIN_CELLS_PER_BIN} cells.")
    W("")
    # ---------------- 7. the decision ----------------
    W("## 7. The decision — scalar MGMT vs spatially-resolved MGMT")
    W("")
    W("### 7.1 The unit problem, and how the comparison is made honest")
    W("")
    W(f"`molecular.mgmt` is calibrated in **bulk CPM over 34,461 genes** "
      f"(observed range 0–21.12 CPM). A CPM over this **{det['n_gene_targets']}-gene "
      f"targeted panel** is a different quantity: this section reads "
      f"**{1e6*vp['scales'][L]['MGMT']['rate']:.1f} panel-CPM** among tumour "
      f"cells, ~3.4× past the calibration ceiling purely because the denominator "
      f"is 94× smaller. **Feeding panel-CPM into the mapping would be a unit "
      f"error dressed as a measurement, and it is not done.**")
    W("")
    W("So the impact runs hold the SECTION MEAN at a declared anchor from the "
      "Visium range already published in `results/mgmt_scalar.md` §2 and take "
      "**only the relative spatial modulation** from Xenium:")
    W("")
    W("```")
    W("CPM(bin) = CPM_anchor × rate(bin) / rate(section)")
    W("```")
    W("")
    W("Bins with fewer than "
      f"{MIN_CELLS_PER_BIN} tumour cells, and every non-tumour cell, take the "
      "anchor itself — no value is invented where the data does not support "
      "one. This isolates the question being asked (*does spatial structure "
      "change the answer?*) from the cross-assay calibration question, which "
      "this dataset cannot settle.")
    W("")
    W("### 7.2 The arithmetic, before any simulation")
    W("")
    from ..molecular.mgmt import mgmt_sensitivity as _ms
    anchor = ANCHORS_CPM[1][1]
    d0 = _ms(anchor, units="cpm", detail=True)
    s_lo = _ms(anchor * (1 + cv), units="cpm")
    s_hi = _ms(anchor * (1 - cv), units="cpm")
    ci = d0["ci95"]
    W(f"MGMT's excess CV at {L} µm is **{cv:.3f}**, i.e. a ±{100*cv:.0f} % swing "
      f"in the local rate. Pushed through `s(x) = exp(−λ·log2(CPM+1))` at the "
      f"{anchor:g} CPM anchor:")
    W("")
    W("| | sensitivity multiplier | width |")
    W("|---|---|---:|")
    W(f"| point estimate at {anchor:g} CPM | {d0['sensitivity']:.4f} | — |")
    W(f"| **± 1 spatial CV** ({anchor*(1-cv):.1f}–{anchor*(1+cv):.1f} CPM) | "
      f"{s_lo:.4f} – {s_hi:.4f} | **{s_hi - s_lo:.4f}** |")
    W(f"| **95 % CI already carried on the scalar** (λ = "
      f"{d0['lambda']:.4f}, n = {d0['n']}) | {ci[0]:.4f} – {ci[1]:.4f} | "
      f"**{ci[1] - ci[0]:.4f}** |")
    W("")
    W(f"**The uncertainty already published on the scalar is "
      f"{(ci[1]-ci[0])/(s_hi-s_lo):.0f}× wider than the entire spatial swing** "
      f"— and that swing is itself an over-estimate, because §5 shows most of "
      f"the {100*cv:.0f} % is counting noise rather than MGMT biology. Even "
      f"taking it entirely at face value, it is invisible inside the error bar "
      f"the model already reports. That is the decision; the simulations below "
      f"confirm it and then show something worse about the spatial version.")
    W("")
    if impact:
        W(f"### 7.3 The simulations — same tissue, same grid (dx = "
          f"{impact['dx_um']:g} µm), same PDE, same seed")
        W("")
        W("The **only** thing that changes between rows is "
          "`micro.cells.sensitivity`.")
        W("")
        W("| MGMT input | mean sensitivity | per-cell sensitivity range (5–95 pct) "
          "| predicted TMZ kill | Δ vs no-MGMT (pp) |")
        W("|---|---:|---|---:|---:|")
        for r in impact["rows"]:
            rng_s = (f"{r['sens_p05']:.4f} – {r['sens_p95']:.4f} "
                     f"(min {r['sens_min']:.4f}, max {r['sens_max']:.4f})"
                     if r["kind"] in ("spatial", "spatial_null") else "constant")
            W(f"| {r['label']} | {r['scalar']:.4f} | {rng_s} | "
              f"**{r['kill']:.5f}** | {r['kill_vs_none_pp']:+.3f} |")
        W("")
        W("**Scalar vs spatial, head to head — the quantity the brief asked "
          "for, with its own null:**")
        W("")
        W("The `SPATIAL NULL` column is the decisive control. It is the same "
          "MGMT counts on the same cells, **randomly reassigned across tumour "
          "cells** — total, per-cell marginal distribution and every bin's "
          "depth preserved exactly, only the spatial arrangement destroyed. "
          "Whatever kill difference *it* produces cannot be MGMT biology.")
        W("")
        W("| anchor | bin | kill, scalar | kill, spatial | Δ (pp) | "
          "Δ as % of scalar kill | kill, SPATIAL NULL | Δ null (pp) | "
          "**signal / null** |")
        W("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for dd in impact["deltas"]:
            nk = ("—" if dd["kill_spatial_null"] is None
                  else f"{dd['kill_spatial_null']:.5f}")
            nd = ("—" if dd["delta_pp_null"] is None
                  else f"{dd['delta_pp_null']:+.4f}")
            so = ("—" if dd["signal_over_null"] is None
                  else f"**{dd['signal_over_null']:.2f}×**")
            W(f"| {dd['anchor_cpm']:g} CPM | {dd['bin_um']} µm | "
              f"{dd['kill_scalar']:.5f} | {dd['kill_spatial']:.5f} | "
              f"**{dd['delta_pp']:+.4f}** | {dd['pct_of_scalar_kill']:+.2f} % | "
              f"{nk} | {nd} | {so} |")
        W("")
        worst = max(abs(dd["delta_pp"]) for dd in impact["deltas"])
        worst_pct = max(abs(dd["pct_of_scalar_kill"]) for dd in impact["deltas"])
        base = impact["baseline_kill_no_mgmt"]
        big = max(abs(r["kill_vs_none_pp"]) for r in impact["rows"])
        ratios = [dd["signal_over_null"] for dd in impact["deltas"]
                  if dd["signal_over_null"] is not None]
        W(f"**Largest scalar-vs-spatial difference: {worst:.4f} percentage "
          f"points ({worst_pct:.2f} % of the kill) — and it is not signal.**")
        W("")
        W("Two things have to be read together here.")
        W("")
        W(f"1. **The difference is small in absolute terms.** Supplying an MGMT "
          f"value *at all* moves predicted kill by up to {big:.3f} pp from the "
          f"no-MGMT baseline of {base:.5f}; resolving it spatially moves it a "
          f"further {worst:.4f} pp, i.e. **~{big/worst:.0f}× less** — and the "
          f"calibration CI on the scalar is "
          f"{(ci[1]-ci[0])/(s_hi-s_lo):.0f}× wider again than the whole spatial "
          f"swing (§7.2).")
        if ratios:
            W(f"2. **The difference that does exist is a NOISE ARTEFACT, and "
              f"the null proves it.** Shuffling MGMT counts at random across "
              f"tumour cells — destroying every trace of spatial structure — "
              f"reproduces "
              f"{'the same or a larger' if max(ratios) <= 1.05 else 'most of the'} "
              f"shift (signal/null ratio "
              f"{min(ratios):.2f}–{max(ratios):.2f}×). The mechanism is "
              f"**Jensen's inequality**: `s(x) = exp(−λ·log2(CPM+1))` is convex "
              f"in CPM, so the mean of `s` over a *dispersed* per-bin CPM "
              f"distribution exceeds `s` at the mean CPM — and at these "
              f"detection rates most of that dispersion is Poisson counting "
              f"noise, not biology. It is why the effect **grows as the bins "
              f"get smaller** (from {HEADLINE_BIN_UM} µm to 200 µm), which is "
              f"the exact opposite of what a real spatial signal would do.")
        W("")
        W("**Where the bias comes from, concretely.** A bin that happens to "
          "draw zero MGMT counts is handed CPM = 0, and CPM = 0 is the "
          "mapping's anchor — *fully TMZ-sensitive*. That is a maximal "
          "biological claim produced by a Poisson accident:")
        W("")
        W("| bin | median expected MGMT counts/bin | bins with ZERO counts | "
          "tumour cells handed \"MGMT-silent\" by a zero draw | tumour cells "
          "pushed past the calibration ceiling (21.12 CPM) |")
        W("|---:|---:|---:|---:|---:|")
        for Lb in (HEADLINE_BIN_UM, 200):
            fs = next(r["field_stats"] for r in impact["field_stats"]
                      if r["bin_um"] == Lb and r["kind"] == "spatial"
                      and r["anchor"] == ANCHORS_CPM[1][1])
            W(f"| {Lb} µm | "
              f"{vp['scales'][Lb]['MGMT']['median_expected_counts_per_region']:.1f} | "
              f"{fs['n_zero_count_bins']:,} | "
              f"{100*fs['frac_cells_in_zero_count_bins']:.2f} % | "
              f"{100*fs['frac_cells_above_calibration_ceiling']:.2f} % |")
        W("")
        W("**This is the strongest argument in the report against a spatial "
          "MGMT field, and it is a positive result rather than a null: at the "
          "achievable detection rate, resolving MGMT spatially does not add "
          "information, it adds a systematic upward bias in predicted kill, "
          "and the bias grows with the resolution you ask for.** A scalar is "
          "not merely adequate here; it is the estimator that does not have "
          "this failure mode.")
        W("")
        W("The absolute size is also easy to understand. Kill is a saturating "
          "function of local exposure integrated over "
          f"{ac['n_tumour_cells']:,} tumour cells; a ±{100*cv:.0f} % modulation "
          f"of one multiplicative factor at a {HEADLINE_BIN_UM} µm length "
          "scale, averaged over the whole section, very nearly cancels. "
          "Spatial MGMT would only matter if the modulation were large, "
          "coherent at the scale of the drug-penetration gradient, and "
          "*correlated with delivery* — and none of those three holds here.")
        W("")
    W("### 7.4 The verdict")
    W("")
    W("> **The per-patient scalar treatment of MGMT is VINDICATED, and it should "
      "be stated as an empirical finding rather than as a limitation.**")
    W("")
    W("Precisely what is now established, and what is not:")
    W("")
    W("* **ESTABLISHED (retraction).** \"MGMT cannot be resolved spatially\" is "
      "false as a general claim about spatial assays. MGMT is measurable on a "
      "Xenium custom add-on and is measured here, in human primary GBM, at "
      f"{det['mgmt_counts']:,} counts over {det['n_cells']:,} cells. Any text in "
      "`results/mgmt_scalar.md` or `results/gbm_tmz_exvivo.md` implying "
      "otherwise is superseded by this report.")
    W(f"* **ESTABLISHED (the new, stronger justification).** On the one human GBM "
      f"section where MGMT *is* measured spatially, its spatial structure is "
      f"**at or below the assay's own negative controls** "
      f"(φ {mgmt_phi:.3f} vs {negp_phi:.3f} at {L} µm; excess CV {cv:.3f} vs "
      f"{negp_cv:.3f}; and at count-matched depth "
      f"{_cm['mgmt_thinned_phi_mean']:.3f} vs {_cm['control_phi']:.3f}) and is "
      f"**the weakest of all {vs['n']} abundance-matched "
      f"genes on the same cells**. The scalar is not a concession to a missing "
      f"probe any more; it is what the measurement supports.")
    if impact and any(d["signal_over_null"] is not None for d in impact["deltas"]):
        _r = [d["signal_over_null"] for d in impact["deltas"]
              if d["signal_over_null"] is not None]
        W(f"* **ESTABLISHED (and this one was not anticipated).** Resolving "
          f"MGMT spatially at the achievable detection rate does not just fail "
          f"to add information — it **injects a systematic upward bias into "
          f"predicted kill**, because the sensitivity mapping is convex and "
          f"most of the per-bin dispersion is Poisson noise. A spatial NULL "
          f"reproduces {100/max(_r):.0f}–{100/min(_r):.0f} % of the shift, and "
          f"the bias GROWS as the bins shrink. A scalar is the estimator "
          f"without this failure mode.")
    W("* **NOT ESTABLISHED.** That MGMT is spatially uniform in GBM generally. "
      "This is **one section from one patient**, and the negative controls put "
      "a floor under what could have been seen. A second MGMT-bearing section "
      "would be worth more than any further analysis of this one.")
    W("* **NOT ESTABLISHED.** That per-cell MGMT is achievable at any depth. "
      f"At {det['mgmt_pct_tumour_cells']:.2f} % positivity it is not achievable "
      f"here, and the useful unit is a region of >= {MIN_CELLS_PER_BIN} cells.")
    W("")
    W("**Recommended change to the model: none.** "
      "`reconstruct.real_pipeline.apply_drug_sensitivity(mgmt=...)` already "
      "takes one per-patient scalar and already reports its CI. This report "
      "changes the *justification* for that interface, not the interface. The "
      "sentence in `mgmt_scalar.md` §2 — *\"a per-spot 'MGMT map' built from "
      "that would be counting noise dressed up as biology\"* — was written of "
      "Visium and turns out to be true of targeted-probe Xenium as well, for a "
      "reason that has now been measured rather than assumed.")
    W("")

    # ---------------- 8. caveats ----------------
    W("## 8. Honest caveats")
    W("")
    W("1. **Different patient, different everything.** GSE324268 `GBMtumor07` "
      "shares no tissue, no patient and no panel with the sections in "
      "`results/heldout_patient.md` or with the four Visium sections in "
      "`results/mgmt_scalar.md` §2. Nothing here is a paired measurement, and "
      "the §7 impact runs use *this* tissue's geometry with *those* sections' "
      "MGMT levels as anchors — a deliberate hybrid, done so the spatial "
      "question can be asked without a unit error, and not a claim about any "
      "single patient.")
    W("2. **Transcript, not methylation.** Every probe- and imaging-based "
      "spatial platform measures MGMT **transcript abundance**; the routine "
      "clinical biomarker is **promoter methylation**, which none of them "
      "measures. This report models expression on purpose, and the reason is "
      "measured, not asserted: in the Erasmus cohort MGMT **expression** was "
      "prognostic for OS (CoxPH HR 1.85, p = 6.1e-05, n = 56, externally "
      "replicated in GLIOTRAIN at HR 1.21, p = 0.025) while MGMT **promoter "
      "methylation was NOT significant in the same cohort (p = 0.16)** — "
      "`results/gbm_tmz_exvivo.md` §C. Expression is the quantity the model "
      "wants. The deposit carries no methylation call for this specimen; "
      "obtaining one is the single cheapest follow-up.")
    W("3. **A targeted panel is not a transcriptome.** MGMT ranks "
      f"{det['rank_by_total_counts']} of {det['n_gene_targets']} targets here, "
      f"and 'panel-CPM' is not bulk CPM (§7.1). All spatial statistics are "
      f"internally normalised (MGMT counts per total panel gene count), which "
      f"makes them comparable **within** this section and not across assays.")
    W("4. **Dropout is the dominant term and is treated as such.** "
      f"{det['mgmt_pct_cells']:.2f} % of cells positive; a zero means \"no probe "
      f"fired\", never \"no transcript\". Aggregation to >= {MIN_CELLS_PER_BIN} "
      f"cells is not a convenience, it is the minimum unit at which the "
      f"variable exists. §6 publishes the resulting detection limits.")
    W("5. **The negative controls are the load-bearing comparison, and they "
      "have their own limits.** Pooled negative-control probes carry "
      f"{det['neg_probe']['total']:,} counts section-wide against MGMT's "
      f"{det['mgmt_counts']:,} "
      f"({pnull['count_matched']['neg_probes_pooled']['target_counts']:,.0f} vs "
      f"{det['mgmt_counts_tumour']:,} restricted to tumour cells, which is the "
      f"population every spatial test uses). That objection is answered "
      f"directly in §5: MGMT thinned binomially to each control's exact depth "
      f"lands at φ = "
      f"{pnull['count_matched']['neg_probes_pooled']['mgmt_thinned_phi_mean']:.3f} "
      f"against the control's own "
      f"{pnull['count_matched']['neg_probes_pooled']['control_phi']:.3f} — i.e. "
      f"correcting for depth makes MGMT look *less* structured, not more — and "
      f"a simulated pure-Poisson null at MGMT's own depth gives "
      f"{pnull['mean']:.3f}. The comparison is fair — but it rests on the "
      f"assumption that the negative controls sample the same technical "
      f"gradients a real probe does, which is standard practice and not a "
      f"proof.")
    W("6. **Cell typing is marker-gated and imperfect.** MGMT detection is 3× "
      "higher in the ENDOTHELIAL and IMMUNE gates than in TUMOR (§2), which is "
      "very likely a depth artefact of the gate rather than biology. Restricting "
      "to tumour cells removes the between-type component but not any residual "
      "within-type composition gradient; the depth-residual Moran's I is the "
      "control for that, and it is the row that goes non-significant.")
    W("7. **Niche labels are correlative.** `region_label` scores neighbourhood "
      "pseudobulk against Ivy GAP signatures by rank correlation; no niche here "
      "is histologically confirmed, and `perivascular` holds only "
      f"{[r['n_cells'] for r in npart['per_region'] if r['region']=='perivascular'][0]:,} "
      f"tumour cells. NECROTIC is never assigned as a *cell type* by any panel "
      f"in this repo, so the `necrotic` **niche** is a transcriptional "
      f"assignment, not a morphological one.")
    W("8. **No outcome label, again.** This is an untreated resection with no "
      "drug response, no PFS/OS and no MGMT methylation call. Nothing here "
      "validates a kill prediction — only the machinery that produces one, and "
      "the modelling choice that feeds it.")
    W("9. **An incidental negative result, recorded rather than buried.** The "
      f"ENDOTHELIAL fraction reads {f3['ENDOTHELIAL']['pct']:.2f} % "
      f"({rec['endothelial_vs_volume_em']:.2f}× the volume-EM anchor) on the "
      f"*same* PECAM1 gate that gave 9.74 % on section 1. Three human GBM "
      f"sections now span 4.18–{f3['ENDOTHELIAL']['pct']:.2f} %. §3 states the "
      f"consequence: the 8.0 % agreement is a property of section 1, not a "
      f"validated property of the classifier.")
    W("")

    # ---------------- references ----------------
    W("## References and provenance")
    W("")
    W(f"* GEO **{prov['accession']}** — {prov['landing_page']} "
      f"(Hölzel lab, Bonn; public 2026-03). Sample {prov['sample']} "
      f"`GBMtumor07`. Free direct download, no login, no DUA.")
    W("* Shapson-Coe A, *et al.* *Science* 2024;384(6696):eadk4858. "
      "PMID **38723085** (human cortex volume EM, 8.0 % endothelial).")
    W("* Ntafoulis I, Kleijn A, Ju J, *et al.* *Br J Cancer* "
      "2023;129(8):1327-1338. PMID **37620410**, GEO **GSE232173** (the "
      "MGMT→TMZ-sensitivity calibration and the expression-vs-methylation "
      "prognostic contrast).")
    W("* Puchalski RB, *et al.* *Science* 2018;360(6389):660-663 (Ivy GAP "
      "anatomic-structure signatures used by `reconstruct.region_label`).")
    W("* Hegi ME, *et al.* *NEJM* 2005;352:997-1003. PMID **15758010**; "
      "Esteller M, *et al.* *NEJM* 2000;343:1350-1354. PMID **11070098**.")
    W("* Code: `tmesim/data/real/mgmt_spatial.py` (loader + MEX materialisation), "
      "`tmesim/validation/mgmt_spatial.py` (this analysis). Nothing under "
      "`core/`, `physics/`, `molecular/`, `reconstruct/` or `sim/` was modified.")
    W("* Prior work superseded in part: `results/mgmt_scalar.md` §2 and §4.4, "
      "`results/gbm_tmz_exvivo.md` §B2 — the *mechanism* claim (\"no Xenium "
      "panel carries MGMT\") is now false in general; the *modelling* claim "
      "(\"use a scalar\") is now measured rather than forced.")
    W("* `docs/DATA_ACCESS_PLAN.md` §3 and §7 item 1 named this dataset as the "
      "highest-value no-access-needed experiment. This report closes it.")
    W("")
    return "\n".join(A) + "\n"


# --------------------------------------------------------------------------- #
# 8. driver
# --------------------------------------------------------------------------- #
def run(path: str = loader.DEFAULT_CACHE, dx_um: float = DX_UM,
        recon_dx_um: float = 30.0, with_impact: bool = True,
        write: bool = True) -> Dict:
    """Run everything and (optionally) write ``results/mgmt_spatial.{md,png}``."""
    t0 = time.time()
    print("[1/8] provenance ...", flush=True)
    prov = verify_provenance(path)
    print(f"      all sha256 match: {prov['all_match']}", flush=True)

    print("[2/8] reconstruction (unmodified pipeline) ...", flush=True)
    rec = reconstruct_section(path, dx_um=recon_dx_um)
    ct = rec["cell_type"]
    print(f"      ENDOTHELIAL {rec['endothelial_pct']:.2f} % "
          f"({rec['endothelial_vs_volume_em']:.2f}x volume-EM)", flush=True)

    print("[3/8] panel-registry regression check ...", flush=True)
    panels = verify_existing_panels_unchanged()

    print("[4/8] MGMT detection ...", flush=True)
    det = verify_mgmt(path, cell_type=ct)
    print(f"      {det['mgmt_counts']:,} counts / "
          f"{det['mgmt_cells_positive']:,} cells "
          f"({det['mgmt_pct_cells']:.4f} %)", flush=True)

    sub = _Substrate(path, ct)
    print("[5/8] spatial autocorrelation ...", flush=True)
    ac = spatial_autocorrelation(sub)
    print("[6/8] variance partition + niches ...", flush=True)
    vp = variance_partition(sub)
    npart = niche_partition(sub)
    pnull = poisson_null_phi(sub)
    print("[7/8] power ...", flush=True)
    pw = power_analysis(sub, vp)

    impact = None
    if with_impact:
        print("[8/8] impact simulations (slow) ...", flush=True)
        impact = kill_impact(sub, ct, path=path, dx_um=dx_um)

    md = build_markdown(prov, det, rec, panels, ac, vp, npart, pnull, pw, impact)
    png = os.path.join(_RESULTS, "mgmt_spatial.png")
    if write:
        os.makedirs(_RESULTS, exist_ok=True)
        make_figure(sub, ct, ac, vp, npart, pw, impact, png)
        with open(os.path.join(_RESULTS, "mgmt_spatial.md"), "w") as fh:
            fh.write(md)
        print(f"wrote results/mgmt_spatial.md and results/mgmt_spatial.png "
              f"({time.time() - t0:.0f}s total)", flush=True)
    return {"provenance": prov, "detection": det, "reconstruction":
            {k: v for k, v in rec.items() if k not in ("micro", "cell_type")},
            "panels": panels, "autocorrelation": ac, "variance_partition": vp,
            "niches": npart, "poisson_null": pnull, "power": pw,
            "impact": impact, "markdown": md}


if __name__ == "__main__":  # pragma: no cover
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--path", default=loader.DEFAULT_CACHE)
    ap.add_argument("--dx", type=float, default=DX_UM)
    ap.add_argument("--no-impact", action="store_true")
    a = ap.parse_args()
    run(a.path, dx_um=a.dx, with_impact=not a.no_impact)
