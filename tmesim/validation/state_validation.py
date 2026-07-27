"""TWO COMBINED INTERNAL VALIDATIONS of the oxygen -> cell-state -> resistance chain,
each tested against MEASURED biology in real tissue.

Both components under test are already known to be fragile, and both are
load-bearing for every absolute kill number this project reports.

TEST 1 — is the SIMULATED OXYGEN FIELD real?
--------------------------------------------
``physics.solve_oxygen`` produces the oxygen field that
``cells.agents.update_cell_states`` bands into HYPOXIC / QUIESCENT / VIABLE, and
those states multiply drug kill through ``PDConfig.hypoxia_resistance`` (0.35)
and ``PDConfig.quiescent_resistance`` (0.50). The field has never been compared
against a measured hypoxia readout, because the 319-gene Xenium Human Brain panel
carries no hypoxia-response gene worth the name.

The second cached section does. ``data_cache/xenium_gbm_io`` is the 10x Xenium
Human Brain Cancer FFPE section on the Immuno-Oncology panel — 816,769 cells,
480 gene targets — and it measures **VEGFA** (63.4% of cells) and **HIF1A**
(54.9%). VEGFA and HIF1A are expressed by the perinecrotic pseudopalisading
glioma cells that define GBM hypoxia (Fischer et al., Brain Pathol 2005, PMID
16389942), so a simulated oxygen field that is describing real tissue oxygen must
be NEGATIVELY correlated with them.

THE CONTROL THAT DECIDES THE TEST. A correlation with the measured hypoxia
markers is worthless on its own, because the oxygen PDE is *driven* by the
measured endothelial cells: its source term is a KDE of them. Anything that
falls off with distance from a vessel will correlate with a hypoxia marker. So
this module always reports the simulated field ALONGSIDE three baselines that
use the same inputs and no physics at all:

  * ``dist_endothelial`` — Euclidean distance to the nearest model-ENDOTHELIAL
    cell (one cKDTree query, no PDE);
  * ``vascular_kde`` — the PDE's OWN SOURCE TERM sampled at each cell, i.e. the
    field before it is diffused. If the solve degrades this, the solve is
    subtracting information;
  * ``local_cell_density`` — cells within 50 um (pure crowding, no vessels).

**If a baseline matches or beats the solved field, the PDE adds nothing** and
this module says so in those words.

TEST 2 — is the CELL-STATE LABELLING real?
------------------------------------------
On the 319-gene section the engine labels the great majority of tumour cells
QUIESCENT purely from oxygen banding, and ``results/dose_schedule.md`` shows
~89% of the temozolomide kill ceiling comes from that quiescent band alone.
QUIESCENT is an empirical claim — "alive but not proliferating" — and ``MKI67``
and ``TOP2A`` are measured per cell on that panel, so it is falsifiable.

THE POWER PROBLEM, HANDLED UP FRONT. MKI67 is detected in only ~0.35% and TOP2A
in ~1.08% of cells on this targeted panel, against published GBM Ki-67 labelling
indices of ~19-40% (:data:`LITERATURE_KI67`). The measurement is therefore
dropout-dominated by one to two orders of magnitude, and a null result is very
easily a power failure rather than a refutation. Every state contrast here
carries an explicit minimum-detectable-effect calculation, and the verdict is
allowed to be **UNDERPOWERED** — which is a real answer, not a hedge.

Caveats that apply to everything below
--------------------------------------
1. **Different patients.** The IO section and the 319-gene section are different
   tumours. Test 1 and Test 2 are two experiments, not two views of one.
2. **VEGFA/HIF1A are hypoxia-RESPONSE genes, not an oxygen tension.** They are
   the transcriptional output of a signalling pathway with its own kinetics, and
   both are induced by non-hypoxic stimuli. HIF1A is the weaker of the two by
   construction: HIF-1alpha is regulated mainly by oxygen-dependent prolyl
   hydroxylation and VHL-mediated proteolysis of the PROTEIN (Ivan et al.,
   Science 2001, PMID 11292862; Jaakkola et al., Science 2001, PMID 11292861),
   so its mRNA need not track oxygen at all.
3. **Targeted-panel dropout.** A zero is "no probe fired", not "no transcript".
4. **Typing circularity.** ``FLT1`` and ``PLVAP`` are IO-panel ENDOTHELIAL
   markers and therefore help *build* the oxygen source; correlations with them
   are circular and are reported as such, never as corroboration.

ADDS ONLY. This module drives the frozen stack (``reconstruct.real_pipeline`` ->
``physics.solve_oxygen`` -> ``cells.agents.update_cell_states``) and reuses the
statistics kernels already written in ``validation.cellstate_validation``. It
modifies nothing in ``core`` / ``physics`` / ``molecular`` / ``sim.engine`` /
``reconstruct``.

Outputs ``results/state_validation.md`` and ``results/state_validation.png``.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import rankdata

from ..cells.agents import update_cell_states
from ..core.config import OxygenConfig, PDConfig
from ..core.microenvironment import CellState, CellType
from ..data.real.singlecell_spatial import load as load_singlecell_spatial
from ..reconstruct.marker_panels import (
    BRAIN_PANEL_MARKER_MAP,
    HYPOXIA_RESPONSE_GENES,
    select_marker_map,
)
from ..reconstruct.real_pipeline import reconstruct_real
# Statistics kernels are REUSED, not re-implemented: identical estimators to
# ``results/cellstate_validation.png`` so the two reports are directly comparable.
from .cellstate_validation import (
    _column,
    _zscore,
    cliffs_delta,
    cohens_d,
    depth_deciles,
    normalized_log1p,
    partial_spearman,
    permutation_test_means,
    spearman_permutation,
    stratified_permutation_test,
)

# ---------------------------------------------------------------------------
# Declared constants. Every modelling choice is named here, not buried.
# ---------------------------------------------------------------------------

#: The two cached real sections. DIFFERENT PATIENTS, different panels.
IO_SAMPLE = "data_cache/xenium_gbm_io"      # 816,769 cells, 480-gene IO panel
BRAIN_SAMPLE = "data_cache/xenium_gbm"      # 40,887 cells, 319-gene brain panel

#: Grid spacings swept in both tests. 12 um is the grid-converged production
#: spacing (``results/numerical_rigor.md``); 30 and 20 expose grid sensitivity.
DX_CONVERGED = 12.0
DX_SWEEP: Tuple[float, ...] = (30.0, 20.0, DX_CONVERGED)

#: Tissue-mask conventions. The mask gates the vascular field (the oxygen
#: SOURCE), so it sets the absolute oxygen level and therefore every state
#: fraction — yet it is a bookkeeping convention, not a measurement.
#:   ``geometric`` — union of fixed-radius discs around measured centroids; the
#:                   current shipping default, grid-independent.
#:   ``occupancy`` — legacy "grid cells containing >= 1 cell"; grid-dependent,
#:                   and the convention under which ``results/dose_schedule.md``
#:                   quoted its 86.0% QUIESCENT / 13.4% HYPOXIC headline.
TISSUE_MASK_MODES: Tuple[str, ...] = ("geometric", "occupancy")
DEFAULT_MASK_MODE = "geometric"

#: The mask+dx configuration ``results/dose_schedule.md`` actually ran under.
#: Test 2's headline contrast is computed there, because that is the
#: configuration whose kill ceilings are at stake.
DOSE_SCHEDULE_CONFIG = ("occupancy", DX_CONVERGED)

#: Permutation draws. 2,000 on the IO section (816k cells — every effect of
#: interest is either p < 5e-4 or indistinguishable from zero, so finer
#: resolution buys nothing) and 10,000 on the 319-gene section, matching
#: ``cellstate_validation``.
N_PERM_IO = 2_000
N_PERM_BRAIN = 10_000

#: DOCUMENTED SUBSAMPLE. Permutation kernels are O(n x draws); at 816,769 cells
#: that is 1.6e9 element-shuffles per test. Permutation p-values on the IO
#: section are therefore computed on a fixed random subsample of this many cells
#: (seed :data:`PERM_SEED`). Every EFFECT SIZE and every Spearman rho is computed
#: on the FULL population — only the null distribution is subsampled, which can
#: only make a p-value more conservative.
IO_PERM_SUBSAMPLE = 50_000

#: Deterministic seed. This module contains no other randomness.
PERM_SEED = 0

#: Measured hypoxia-RESPONSE genes on the IO panel. Not an oxygen tension —
#: see caveat 2 in the module docstring.
HYPOXIA_GENES: Tuple[str, ...] = ("VEGFA", "HIF1A")

#: Measured proliferation markers on the 319-gene panel. TOP2A is NOT on the IO
#: panel and PCNA/CDK1 are not either, which is why Test 2 cannot be run there.
PROLIFERATION_MARKERS: Tuple[str, ...] = ("MKI67", "TOP2A")

#: Radius (um) for the ``local_cell_density`` trivial baseline. ~3x the measured
#: median nearest-neighbour spacing, so the count is a stable local crowding
#: estimate rather than a nearest-neighbour coin flip.
DENSITY_RADIUS_UM = 50.0

#: Percentile positivity rule for Test 2: "a cell is proliferating if it is in
#: the top ``100 - PERCENTILE_RULE_Q`` % of depth-normalised marker expression
#: among tumour cells". 80 is chosen to BRACKET the published Ki-67 index
#: (:data:`LITERATURE_KI67`) — if the panel measured proliferation the way IHC
#: does, the top 20% would be roughly the Ki-67-positive population.
#:
#: HONEST NOTE: with >99% of cells at exactly zero this threshold falls inside a
#: tie mass at 0 and the rule is DEGENERATE. :func:`positivity_report` detects
#: that and says so instead of silently returning a meaningless number.
PERCENTILE_RULE_Q = 80.0

#: Minimum cells per tile for the coarse-tile spatial (grid-level) correlation.
#: A tile mean over 1-2 cells is mostly dropout noise.
TILE_MIN_CELLS = 5

#: Fixed physical tile size (um) for the grid-level correlation, held constant
#: across the dx sweep so the three numbers are comparable to each other rather
#: than to the mesh.
TILE_UM = 60.0

#: Two-sided alpha and target power for the minimum-detectable-effect statements.
MDE_ALPHA = 0.05
MDE_POWER = 0.80

# ---------------------------------------------------------------------------
# Literature anchors. Every entry is PubMed-verified; nothing is estimated.
# ---------------------------------------------------------------------------

#: Published Ki-67 labelling indices in glioblastoma. The task's "~15-30%"
#: sits at the low end of this range; the published central values run ~19-40%.
#: Either way the panel's ~0.44% MKI67 detection rate is 40-90x below them.
LITERATURE_KI67: Tuple[Dict[str, object], ...] = (
    {
        "label": "GBM, n=35 (prognostic cutoff)",
        "value_pct": 19.0,
        "citation": "Fudaba et al., World Neurosurg 2021;149:e469-e480",
        "pmid": "33567368",
        "doi": "10.1016/j.wneu.2021.02.002",
        "note": "Ki-67 LI cutoff 19% was an independent predictor of overall survival",
    },
    {
        "label": "GBM, n=44 (cohort split point)",
        "value_pct": 27.0,
        "citation": "Alkhaibary et al., Hematol Oncol Stem Cell Ther 2018;12(2):82-88",
        "pmid": "30552865",
        "doi": "10.1016/j.hemonc.2018.11.001",
        "note": "cohort dichotomised at a Ki-67 LI of 27%; LI alone did not predict OS",
    },
    {
        "label": "epithelioid GBM, n=16 (median, range 10-50%)",
        "value_pct": 30.0,
        "citation": "Zheng et al., Zhonghua Bing Li Xue Za Zhi 2021;50(3):229-235",
        "pmid": "33677887",
        "doi": "10.3760/cma.j.cn112151-20200617-00479",
        "note": "median Ki-67 index 30% (range 10-50%)",
    },
    {
        "label": "GBM, n=70 (median)",
        "value_pct": 40.0,
        "citation": "Caramanti et al., Cureus 2022;14(4):e24246",
        "pmid": "35602791",
        "doi": "10.7759/cureus.24246",
        "note": "median Ki-67 index 40% across 70 GBM resections",
    },
)

#: The directional prediction the state model must reproduce, measured in tissue.
LITERATURE_HYPOXIA_PROLIFERATION: Dict[str, object] = {
    "claim": (
        "hypoxic (pseudopalisading, HIF-1alpha-positive) glioblastoma cells are "
        "5-50% LESS proliferative and 6-20x more apoptotic than adjacent astrocytoma"
    ),
    "citation": "Brat et al., Cancer Res 2004;64(3):920-927",
    "pmid": "14871821",
    "doi": "10.1158/0008-5472.CAN-03-2073",
}

#: Why VEGFA is the right measured readout of hypoxia in THIS tumour type.
LITERATURE_VEGF_HYPOXIA: Dict[str, object] = {
    "claim": (
        "in GBM, vascular involution produces necrosis and hypoxia, which induces "
        "HIF-1alpha and VEGF expression specifically in the perinecrotic "
        "pseudopalisading glioma cells"
    ),
    "citation": "Fischer et al., Brain Pathol 2005;15(4):297-310",
    "pmid": "16389942",
    "doi": "10.1111/j.1750-3639.2005.tb00115.x",
}

#: Why HIF1A *mRNA* is a weak hypoxia readout — the mechanism is proteolytic.
LITERATURE_HIF1A_POSTTRANSLATIONAL: Tuple[Dict[str, object], ...] = (
    {
        "claim": (
            "HIF-1alpha is targeted for VHL-mediated destruction by O2-dependent "
            "prolyl hydroxylation — oxygen sensing acts on the PROTEIN, not the mRNA"
        ),
        "citation": "Ivan et al., Science 2001;292(5516):464-468",
        "pmid": "11292862",
        "doi": "10.1126/science.1059817",
    },
    {
        "claim": (
            "HIF-alpha is targeted to the VHL ubiquitylation complex by O2-regulated "
            "prolyl hydroxylation; HIF-PH requires dioxygen as cosubstrate"
        ),
        "citation": "Jaakkola et al., Science 2001;292(5516):468-472",
        "pmid": "11292861",
        "doi": "10.1126/science.1059796",
    },
)

#: Ki-67 interobserver variability — the reference standard is itself imprecise.
LITERATURE_KI67_VARIABILITY: Dict[str, object] = {
    "claim": (
        "Ki-67 LI scoring in gliomas agrees poorly between observers: kappa 0.32 "
        "for whole-slide and as low as 0.04 for hot-spot scoring (n=235 gliomas)"
    ),
    "citation": "Nielsen et al., Diagn Pathol 2018;13:38",
    "pmid": "29885671",
    "doi": "10.1186/s13000-018-0711-2",
}

# ---------------------------------------------------------------------------
# Palette. Validated with the dataviz skill's checker on the light surface
# (#fcfcfb): categorical slots {red 6, yellow 3, blue 1} pass the lightness
# band, the chroma floor and CVD separation (worst adjacent dE 37.7, deutan);
# yellow is 2.11:1 against the surface, sub-3:1, so the RELIEF RULE applies and
# every mark drawn in yellow also carries a direct value label. The
# model-vs-baseline pair {orange 8, blue 1} passes all four checks (dE 96.7).
# ---------------------------------------------------------------------------
STATE_COLORS: Dict[str, str] = {
    "HYPOXIC": "#e34948",    # slot 6, red
    "QUIESCENT": "#eda100",  # slot 3, yellow  (sub-3:1 -> always direct-labelled)
    "VIABLE": "#2a78d6",     # slot 1, blue
}
MODEL_COLOR = "#eb6834"      # slot 8, orange — the simulated (PDE) quantity
BASELINE_COLOR = "#2a78d6"   # slot 1, blue   — the trivial baselines
MEASURED_COLOR = "#1baf7a"   # slot 2, aqua   — measured quantities
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8984"
GRID_COLOR = "#e4e3df"
SURFACE = "#fcfcfb"

#: Sequential blue ramp (magnitude, light->dark) for the simulated O2 map, and
#: the aqua hue for the second simultaneous sequential context (measured score).
SEQ_BLUE = ("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#104281")
SEQ_AQUA = ("#d6f2e6", "#a5e0c8", "#6fcaa5", "#31b183", "#1baf7a", "#0d6b48")

#: Ordered low-oxygen -> high-oxygen. Used for every state axis in the report.
STATE_ORDER: Tuple[str, ...] = ("HYPOXIC", "QUIESCENT", "VIABLE")

_STATE_ENUM = {
    "HYPOXIC": CellState.HYPOXIC,
    "QUIESCENT": CellState.QUIESCENT,
    "VIABLE": CellState.VIABLE,
}


def _results_dir() -> str:
    d = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results")
    )
    os.makedirs(d, exist_ok=True)
    return d


def load_section(path: str) -> Dict:
    """Load a cached real Xenium section, retrying once on a transient failure."""
    try:
        return load_singlecell_spatial(path=path, dense=False)
    except Exception:
        return load_singlecell_spatial(path=path, dense=False)


# ---------------------------------------------------------------------------
# Driving the FROZEN pipeline: reconstruct -> physics.solve_oxygen -> banding.
# ---------------------------------------------------------------------------
def band_states(
    sample: str, dx_um: float, tissue_mask_mode: str = DEFAULT_MASK_MODE
) -> Dict:
    """Reconstruct one section at one ``dx`` and apply the shipping state banding.

    Calls the unmodified :func:`reconstruct.real_pipeline.reconstruct_real` (which
    itself runs ``physics.solve_oxygen`` on the real vascular field) and then the
    unmodified :func:`cells.agents.update_cell_states`. NOTHING here re-implements
    the banding — the object under test is the production code path.

    Returns the per-cell oxygen, states, types, positions, the raw fields and the
    grid, plus the tumour-cell state fractions. Rows are in the loader's
    deterministic order, i.e. index-aligned to the expression matrix (asserted by
    the callers, never assumed).
    """
    cfg = OxygenConfig()
    micro = reconstruct_real(
        source="xenium", sample=sample, dx_um=float(dx_um),
        tissue_mask_mode=str(tissue_mask_mode),
    )
    oxygen = np.asarray(micro.fields.oxygen, dtype=np.float64)
    update_cell_states(micro, oxygen, cfg)

    cells = micro.cells
    col, row = micro.domain.to_index(cells.x, cells.y)
    local_o2 = oxygen[row, col]
    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)[row, col]

    cell_type = np.asarray(cells.cell_type)
    tumor = cell_type == int(CellType.TUMOR)
    n_tumor = int(tumor.sum())
    fractions = {
        name: (float((cells.state[tumor] == int(_STATE_ENUM[name])).sum()) / n_tumor
               if n_tumor else 0.0)
        for name in STATE_ORDER
    }
    fractions_all = {
        name: float((cells.state == int(_STATE_ENUM[name])).mean())
        for name in STATE_ORDER
    }

    return {
        "sample": sample,
        "dx_um": float(dx_um),
        "tissue_mask_mode": str(tissue_mask_mode),
        "grid": tuple(int(v) for v in micro.domain.shape),
        "state": cells.state.copy(),
        "cell_type": cell_type.copy(),
        "local_o2": local_o2,
        "vascular_at_cell": vasc,
        "x": cells.x.copy(),
        "y": cells.y.copy(),
        "row": row,
        "col": col,
        "tumor_mask": tumor,
        "n_cells": int(cells.n),
        "n_tumor": n_tumor,
        "fractions_tumor": fractions,
        "fractions_all": fractions_all,
        "hypoxia_threshold": float(cfg.hypoxia_threshold),
        "mean_oxygen_tissue": float(oxygen[micro.domain.tissue_mask].mean()),
        "mean_oxygen_at_cells": float(local_o2.mean()),
        "oxygen_field": oxygen,
        "tissue_mask": np.asarray(micro.domain.tissue_mask, dtype=bool),
    }


# ---------------------------------------------------------------------------
# THE TRIVIAL BASELINES — the control that decides Test 1.
# ---------------------------------------------------------------------------
def trivial_baselines(
    x: np.ndarray, y: np.ndarray, cell_type: np.ndarray,
    vascular_at_cell: np.ndarray, radius_um: float = DENSITY_RADIUS_UM,
) -> Dict[str, Dict[str, object]]:
    """Physics-free predictors of hypoxia built from the SAME inputs as the PDE.

    The oxygen PDE consumes exactly two spatial inputs: the ENDOTHELIAL cell
    positions (source) and the TUMOR cell positions (sink). If a one-line function
    of those same inputs predicts the measured hypoxia markers as well as the
    solved field does, then the reaction-diffusion solve is not contributing.

    Each baseline is returned with the SIGN CONVENTION ALIGNED TO OXYGEN — i.e.
    larger means "more oxygenated" — so its Spearman rho against a hypoxia score
    is directly comparable to the simulated field's, no mental sign-flipping.

    ``dist_endothelial``
        Negated Euclidean distance (um) to the nearest model-ENDOTHELIAL cell, via
        one ``cKDTree`` query. This is the crudest possible perfusion proxy.
    ``vascular_kde``
        The PDE's own source term sampled at each cell — the Gaussian KDE of the
        endothelial cells BEFORE the diffusion solve. Comparing against this
        isolates the contribution of the solve itself: if the solved field is
        worse, the solve destroyed information that its input already had.
    ``local_cell_density``
        Negated count of cells within ``radius_um``. Pure crowding, no vessels; a
        sanity floor for "any spatial variable at all correlates with anything".
    """
    xy = np.column_stack([np.asarray(x, float), np.asarray(y, float)])
    endo = np.asarray(cell_type) == int(CellType.ENDOTHELIAL)
    n_endo = int(endo.sum())

    if n_endo:
        dist_endo, _ = cKDTree(xy[endo]).query(xy, k=1)
    else:
        dist_endo = np.full(xy.shape[0], np.nan)

    tree_all = cKDTree(xy)
    n_within = np.asarray(
        tree_all.query_ball_point(xy, r=float(radius_um), return_length=True),
        dtype=np.float64,
    )

    return {
        "dist_endothelial": {
            "values": -dist_endo,
            "label": f"-(distance to nearest ENDOTHELIAL cell), n_endo={n_endo:,}",
            "what": "one cKDTree query on the PDE's own source cells; no physics",
            "median_um": float(np.median(dist_endo)) if n_endo else float("nan"),
        },
        "vascular_kde": {
            "values": np.asarray(vascular_at_cell, float),
            "label": "vascular-density KDE at the cell (the PDE's SOURCE TERM)",
            "what": "the solver's own input field, undiffused",
            "median_um": float("nan"),
        },
        "local_cell_density": {
            "values": -n_within,
            "label": f"-(cells within {radius_um:.0f} um)",
            "what": "pure local crowding; ignores vessels entirely",
            "median_um": float("nan"),
        },
    }


# ---------------------------------------------------------------------------
# Small statistics helpers not already in cellstate_validation.
# ---------------------------------------------------------------------------
def spearman_full(x: np.ndarray, y: np.ndarray) -> float:
    """Tie-corrected Spearman rho on the FULL population (no permutation).

    ``rankdata`` averages ties, which is essential here: the marker vectors are
    heavily zero-inflated, so a naive rank would fabricate an ordering inside the
    tie mass.
    """
    rx = rankdata(np.asarray(x, float))
    ry = rankdata(np.asarray(y, float))
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denom = np.sqrt(float((rx ** 2).sum()) * float((ry ** 2).sum()))
    return float((rx * ry).sum() / denom) if denom > 0 else 0.0


def subsample_index(n: int, k: int, seed: int = PERM_SEED) -> np.ndarray:
    """Fixed, seeded random subsample of ``k`` row indices out of ``n`` (or all)."""
    if k >= n:
        return np.arange(n)
    return np.sort(np.random.default_rng(seed).choice(n, size=int(k), replace=False))


def min_detectable_rate_ratio(
    p0: float, n_group: int, n_rest: int,
    alpha: float = MDE_ALPHA, power: float = MDE_POWER,
) -> Dict[str, float]:
    """Smallest marker-positive RATE RATIO the group contrast could have detected.

    A null result only means something if the test COULD have rejected. Using the
    normal approximation to a two-proportion test, the minimum detectable absolute
    difference at two-sided ``alpha`` and ``power`` is

        delta = (z_{alpha/2} + z_{power}) * sqrt(p0 (1 - p0) (1/n1 + 1/n0))

    and the reported ``rate_ratio`` is ``(p0 + delta) / p0`` — "the group would have
    had to be at least this many times more marker-positive for us to have seen
    it". A rate ratio far above any biologically plausible enrichment (Ki-67 in
    cycling vs quiescent cells differs by well under an order of magnitude) means
    the test is UNDERPOWERED and its null carries no information.

    Also returns ``expected_positives``, the count the group would contain under
    the null — below ~5 no rate estimate from it is meaningful at all.
    """
    z_a, z_b = 1.959963985, 0.8416212336  # alpha=0.05 two-sided, power=0.80
    if not (0.0 < p0 < 1.0) or n_group <= 0 or n_rest <= 0:
        return {"rate_ratio": float("inf"), "abs_difference": float("nan"),
                "expected_positives": 0.0, "n_group": int(n_group)}
    se = np.sqrt(p0 * (1.0 - p0) * (1.0 / n_group + 1.0 / n_rest))
    delta = (z_a + z_b) * se
    return {
        "rate_ratio": float((p0 + delta) / p0),
        "abs_difference": float(delta),
        "expected_positives": float(p0 * n_group),
        "n_group": int(n_group),
        "alpha": float(alpha),
        "power": float(power),
    }


def positivity_report(
    counts: np.ndarray, totals: np.ndarray, subset: np.ndarray,
    q: float = PERCENTILE_RULE_Q,
) -> Dict[str, object]:
    """Two DOCUMENTED positivity rules for one marker over ``subset`` cells.

    ``nonzero``
        ``raw count >= 1``. The permissive rule: "this probe fired in this cell".
        It is the closest thing a targeted in-situ panel has to a detection call.
    ``percentile``
        depth-normalised log1p expression strictly above the ``q``-th percentile of
        the subset — i.e. "the top ``100-q``%", chosen to bracket the published
        Ki-67 index so the two are on the same footing.

    The percentile rule is reported as **degenerate** when the ``q``-th percentile
    is 0, because then the threshold sits inside the tie mass at zero and the rule
    silently collapses to ``nonzero`` (or to nothing). That degeneracy IS the
    result — it is the dropout problem stated quantitatively — so it is surfaced
    rather than papered over.
    """
    counts = np.asarray(counts, float)
    subset = np.asarray(subset, bool)
    c = counts[subset]
    n = int(c.size)
    norm = normalized_log1p(counts, totals)[subset]

    thresh = float(np.percentile(norm, q)) if n else 0.0
    degenerate = bool(thresh <= 0.0)
    pct_positive = norm > thresh if not degenerate else (c >= 1)

    return {
        "n_cells": n,
        "zero_fraction": float((c == 0).mean()) if n else 1.0,
        "mean_count": float(c.mean()) if n else 0.0,
        "max_count": float(c.max()) if n else 0.0,
        "nonzero": {
            "rule": "raw count >= 1",
            "n_positive": int((c >= 1).sum()),
            "fraction": float((c >= 1).mean()) if n else 0.0,
            "mask": c >= 1,
        },
        "percentile": {
            "rule": (f"depth-normalised log1p > p{q:.0f} of tumour cells"
                     + (" [DEGENERATE: p%.0f == 0, fell back to count >= 1]" % q
                        if degenerate else "")),
            "threshold": thresh,
            "degenerate": degenerate,
            "n_positive": int(pct_positive.sum()),
            "fraction": float(pct_positive.mean()) if n else 0.0,
            "mask": pct_positive,
        },
    }


def rate_by_state(positive: np.ndarray, states: np.ndarray) -> Dict[str, Dict[str, float]]:
    """Marker-positive rate per model state with a Wilson 95% interval.

    Wilson rather than the normal approximation because several of these rates are
    a handful of cells out of a few hundred, where the normal interval is nonsense.
    """
    out: Dict[str, Dict[str, float]] = {}
    for name in STATE_ORDER:
        m = states == int(_STATE_ENUM[name])
        n = int(m.sum())
        if n == 0:
            out[name] = {"n": 0, "k": 0, "rate": float("nan"),
                         "lo": float("nan"), "hi": float("nan")}
            continue
        k = int(np.asarray(positive, float)[m].sum())
        p = k / n
        z = 1.959963985
        denom = 1.0 + z * z / n
        centre = (p + z * z / (2 * n)) / denom
        half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
        out[name] = {"n": n, "k": k, "rate": float(p),
                     "lo": float(max(0.0, centre - half)),
                     "hi": float(min(1.0, centre + half))}
    return out


# ---------------------------------------------------------------------------
# TEST 1 — simulated oxygen vs MEASURED hypoxia markers on the IO section.
# ---------------------------------------------------------------------------
def hypoxia_score(
    expr, genes: np.ndarray, totals: np.ndarray, gene_set: Sequence[str] = HYPOXIA_GENES,
) -> Dict[str, object]:
    """Measured per-cell hypoxia-RESPONSE score: mean z of depth-normalised markers.

    Each gene is counts-per-median normalised then log1p'd (:func:`normalized_log1p`,
    the repo's standard depth correction) and z-scored across all cells, so the two
    genes contribute on a comparable scale. The composite is their mean.

    The per-gene vectors are returned alongside the composite because averaging is
    only defensible if the genes actually co-vary — and here they barely do (the
    measured VEGFA/HIF1A Spearman rho is reported as ``gene_gene_rho``). A near-zero
    inter-gene correlation means the composite is not measuring one latent
    "hypoxia" variable, and the report says so.
    """
    used: List[str] = []
    missing: List[str] = []
    per_gene: Dict[str, np.ndarray] = {}
    raw_counts: Dict[str, np.ndarray] = {}
    for g in gene_set:
        raw = _column(expr, genes, g)
        if raw is None:
            missing.append(g)
            continue
        used.append(g)
        raw_counts[g] = raw
        per_gene[g] = _zscore(normalized_log1p(raw, totals))

    if not per_gene:
        n = int(expr.shape[0])
        return {"score": np.zeros(n), "per_gene": {}, "raw_counts": {},
                "used": used, "missing": missing, "gene_gene_rho": float("nan"),
                "detection": {}}

    score = np.mean(np.vstack([per_gene[g] for g in used]), axis=0)
    gg = (spearman_full(per_gene[used[0]], per_gene[used[1]])
          if len(used) == 2 else float("nan"))
    return {
        "score": score,
        "per_gene": per_gene,
        "raw_counts": raw_counts,
        "used": used,
        "missing": missing,
        "gene_gene_rho": gg,
        "detection": {g: float((raw_counts[g] > 0).mean()) for g in used},
    }


def all_gene_ranking(
    expr, genes: np.ndarray, totals: np.ndarray, oxygen_at_cell: np.ndarray,
) -> List[Tuple[str, float, float]]:
    """Spearman rho of EVERY panel gene against the simulated oxygen, sorted.

    This is the check that turns "VEGFA correlates" into a falsifiable statement.
    A single gene correlating at rho = -0.15 is unremarkable; VEGFA being the most
    anti-correlated of all 480 measured targets is not. It also exposes the
    circularity honestly — the genes at the POSITIVE end are the vascular and
    myeloid programmes, several of which are the ENDOTHELIAL markers that build
    the oxygen source in the first place.

    Returns ``[(gene, rho, detection_rate)]`` sorted most-negative first. Genes are
    depth-normalised before ranking so the ordering is not a depth ordering.
    """
    genes = np.asarray(genes, dtype=str)
    totals = np.asarray(totals, float)
    scale = float(np.median(totals[totals > 0])) if np.any(totals > 0) else 1.0
    safe = np.where(totals > 0, totals, 1.0)

    ro = rankdata(np.asarray(oxygen_at_cell, float))
    ro = ro - ro.mean()
    no = float(np.sqrt((ro ** 2).sum()))

    csc = expr.tocsc() if sparse.issparse(expr) else expr
    out: List[Tuple[str, float, float]] = []
    for j in range(len(genes)):
        if sparse.issparse(csc):
            cnt = np.asarray(csc[:, [j]].todense()).ravel().astype(np.float64)
        else:
            cnt = np.asarray(csc)[:, j].astype(np.float64)
        v = np.log1p(cnt / safe * scale)
        rv = rankdata(v)
        rv = rv - rv.mean()
        nv = float(np.sqrt((rv ** 2).sum()))
        rho = float((ro * rv).sum() / (no * nv)) if nv > 0 and no > 0 else 0.0
        out.append((str(genes[j]), rho, float((cnt > 0).mean())))
    out.sort(key=lambda r: r[1])
    return out


def tile_correlation(
    row: np.ndarray, col: np.ndarray, dx_um: float, grid: Tuple[int, int],
    oxygen_field: np.ndarray, score: np.ndarray,
    tile_um: float = TILE_UM, min_cells: int = TILE_MIN_CELLS,
) -> Dict[str, float]:
    """FIELD-vs-FIELD correlation on fixed physical tiles (not on the mesh).

    The per-cell test asks "does a cell in low simulated O2 express more VEGFA?".
    This asks the coarser and arguably fairer question: "does a REGION the model
    calls hypoxic express more VEGFA on average?" Averaging over a tile suppresses
    per-cell dropout, which is the dominant noise source, so this is the version of
    the test most favourable to the model.

    Tiles are a FIXED ``tile_um`` square regardless of the solver's ``dx``, so the
    three dx values produce comparable numbers instead of numbers that drift with
    the mesh. Tiles holding fewer than ``min_cells`` measured cells are dropped —
    a tile mean over one or two cells is mostly dropout.
    """
    ny, nx = grid
    oxygen_field = np.asarray(oxygen_field, float)
    score = np.asarray(score, float)

    # Cell centre in um, then binned into fixed physical tiles.
    cx = (np.asarray(col, float) + 0.5) * float(dx_um)
    cy = (np.asarray(row, float) + 0.5) * float(dx_um)
    tx = (cx / float(tile_um)).astype(np.int64)
    ty = (cy / float(tile_um)).astype(np.int64)
    ntx, nty = int(tx.max()) + 1, int(ty.max()) + 1
    flat = ty * ntx + tx

    n_tiles = nty * ntx
    cnt = np.bincount(flat, minlength=n_tiles).astype(np.float64)
    ssum = np.bincount(flat, weights=score, minlength=n_tiles)
    osum = np.bincount(flat, weights=oxygen_field[row, col], minlength=n_tiles)

    keep = cnt >= int(min_cells)
    if int(keep.sum()) < 10:
        return {"rho": float("nan"), "pearson": float("nan"), "n_tiles": int(keep.sum()),
                "tile_um": float(tile_um), "min_cells": int(min_cells)}
    mean_score = ssum[keep] / cnt[keep]
    mean_o2 = osum[keep] / cnt[keep]
    return {
        "rho": spearman_full(mean_o2, mean_score),
        "pearson": float(np.corrcoef(mean_o2, mean_score)[0, 1]),
        "n_tiles": int(keep.sum()),
        "tile_um": float(tile_um),
        "min_cells": int(min_cells),
        "_mean_o2": mean_o2,
        "_mean_score": mean_score,
    }


def run_test1(
    sample: str = IO_SAMPLE, n_perm: int = N_PERM_IO,
    dx_sweep: Sequence[float] = DX_SWEEP,
) -> Dict:
    """Test 1 in full: (a) correlation, (b) state enrichment, (c) spatial, (d) dx."""
    print("\n" + "=" * 78)
    print("TEST 1 — SIMULATED OXYGEN vs MEASURED HYPOXIA MARKERS (IO section)")
    print("=" * 78)

    data = load_section(sample)
    expr = data["expression"]
    genes = np.asarray(data["genes"], dtype=str)
    totals = np.asarray(expr.sum(axis=1)).ravel().astype(np.float64)
    n_cells = int(expr.shape[0])
    _, panel_name, panel_cov = select_marker_map(genes)
    print(f"  {n_cells:,} cells x {len(genes)} genes — panel selected: {panel_name} "
          f"(coverage {panel_cov})")

    hs = hypoxia_score(expr, genes, totals)
    score = hs["score"]
    print(f"  measured hypoxia genes: {hs['used']} "
          + "  ".join(f"{g} det={100*hs['detection'][g]:.2f}%" for g in hs["used"]))
    print(f"  measured VEGFA-vs-HIF1A Spearman rho = {hs['gene_gene_rho']:+.4f} "
          "(if ~0 they are NOT one latent hypoxia variable)")

    log_depth = np.log(np.maximum(totals, 1.0))

    # ---- (d) the dx sweep drives everything else --------------------------
    sweep: Dict[float, Dict] = {}
    for dx in dx_sweep:
        print(f"\n  --- dx = {dx:.0f} um ---")
        st = band_states(sample, dx, DEFAULT_MASK_MODE)
        if st["n_cells"] != n_cells:
            raise RuntimeError(
                f"row misalignment at dx={dx}: {st['n_cells']} reconstructed cells "
                f"vs {n_cells} expression rows — no per-cell comparison is valid"
            )
        o2 = st["local_o2"]
        base = trivial_baselines(st["x"], st["y"], st["cell_type"],
                                 st["vascular_at_cell"])

        preds: Dict[str, Dict[str, object]] = {
            "simulated_o2": {
                "values": o2,
                "label": "SIMULATED oxygen (physics.solve_oxygen)",
                "what": "the reaction-diffusion solve under test",
                "median_um": float("nan"),
            }
        }
        preds.update(base)

        per_pred: Dict[str, Dict] = {}
        for key, spec in preds.items():
            v = np.asarray(spec["values"], float)
            entry = {
                "label": spec["label"],
                "what": spec["what"],
                "rho_score": spearman_full(v, score),
                "partial_rho_depth_removed": partial_spearman(v, score, log_depth),
            }
            for g in hs["used"]:
                entry[f"rho_{g}"] = spearman_full(v, hs["per_gene"][g])
            per_pred[key] = entry
            print(f"    {key:<20} rho(score)={entry['rho_score']:+.4f}  "
                  + "  ".join(f"rho({g})={entry['rho_'+g]:+.4f}" for g in hs["used"]))

        # Nested comparison: does the PDE add anything the raw geometry lacks?
        de = np.asarray(base["dist_endothelial"]["values"], float)
        nested = {
            "score_vs_simo2_given_dist": partial_spearman(score, o2, de),
            "score_vs_dist_given_simo2": partial_spearman(score, de, o2),
            "simo2_vs_dist": spearman_full(o2, de),
        }
        print(f"    NESTED: rho(score, simO2 | dist_endo) = "
              f"{nested['score_vs_simo2_given_dist']:+.4f}   vs   "
              f"rho(score, dist_endo | simO2) = {nested['score_vs_dist_given_simo2']:+.4f}")

        # (c) field-vs-field on fixed physical tiles
        tiles = tile_correlation(st["row"], st["col"], dx, st["grid"],
                                 st["oxygen_field"], score)
        print(f"    TILE-LEVEL ({tiles['tile_um']:.0f} um tiles, >={tiles['min_cells']} "
              f"cells, n={tiles['n_tiles']:,}): rho={tiles['rho']:+.4f} "
              f"pearson={tiles['pearson']:+.4f}")

        sweep[dx] = {
            "dx_um": float(dx),
            "grid": st["grid"],
            "mean_oxygen_tissue": st["mean_oxygen_tissue"],
            "mean_oxygen_at_cells": st["mean_oxygen_at_cells"],
            "fractions_all": st["fractions_all"],
            "fractions_tumor": st["fractions_tumor"],
            "n_tumor": st["n_tumor"],
            "predictors": per_pred,
            "nested": nested,
            "tiles": tiles,
            "_state": st,
        }

    ref = sweep[DX_CONVERGED]
    st = ref["_state"]
    o2 = st["local_o2"]
    states = st["state"]
    tumor = st["tumor_mask"]

    # ---- (a) permutation p-values on a DOCUMENTED subsample ----------------
    idx = subsample_index(n_cells, IO_PERM_SUBSAMPLE)
    print(f"\n  (a) permutation p-values on a seeded subsample of {idx.size:,} cells "
          f"({n_perm:,} draws); effect sizes are on all {n_cells:,}")
    perm: Dict[str, Dict] = {}
    for key in ("simulated_o2", "dist_endothelial", "vascular_kde"):
        v = (o2 if key == "simulated_o2"
             else np.asarray(trivial_baselines(
                 st["x"], st["y"], st["cell_type"], st["vascular_at_cell"]
             )[key]["values"], float))
        perm[key] = spearman_permutation(v[idx], score[idx], n_perm, PERM_SEED)
        print(f"    {key:<20} rho(sub)={perm[key]['rho']:+.4f}  p={perm[key]['p_value']:.5f}")

    # ---- (b) model-HYPOXIC vs model-VIABLE enrichment ----------------------
    print("\n  (b) model-HYPOXIC vs model-VIABLE, measured marker enrichment")
    hyp = states == int(CellState.HYPOXIC)
    via = states == int(CellState.VIABLE)
    qui = states == int(CellState.QUIESCENT)
    strata = depth_deciles(totals)

    enrichment: Dict[str, Dict] = {}
    for label, vals in [("hypoxia_score", score)] + [
            (g, hs["per_gene"][g]) for g in hs["used"]]:
        block: Dict[str, Dict] = {}
        for cname, ma, mb in (("hypoxic_vs_viable", hyp, via),
                              ("hypoxic_vs_quiescent", hyp, qui),
                              ("hypoxic_vs_viable_TUMOUR", hyp & tumor, via & tumor)):
            a, b = vals[ma], vals[mb]
            pair = ma | mb
            # Permutation on a subsample of the LARGE group, keeping the small
            # group whole (documented; can only make the p-value conservative).
            sub = subsample_index(int(pair.sum()), IO_PERM_SUBSAMPLE)
            pv = permutation_test_means(vals[pair][sub], ma[pair][sub], n_perm, PERM_SEED)
            block[cname] = {
                "n_a": int(ma.sum()), "n_b": int(mb.sum()),
                "mean_a": float(a.mean()) if a.size else float("nan"),
                "mean_b": float(b.mean()) if b.size else float("nan"),
                "cohens_d": cohens_d(a, b),
                "cliffs_delta": cliffs_delta(a, b),
                "p_value": pv["p_value"],
                "n_permuted": int(sub.size),
            }
            print(f"    {label:<14} {cname:<26} d={block[cname]['cohens_d']:+.4f} "
                  f"cliff={block[cname]['cliffs_delta']:+.4f} p={pv['p_value']:.5f} "
                  f"(n={block[cname]['n_a']:,} vs {block[cname]['n_b']:,})")
        # Depth-stratified control: same contrast, labels permuted within decile.
        pair = hyp | via
        sub = subsample_index(int(pair.sum()), IO_PERM_SUBSAMPLE)
        block["hypoxic_vs_viable_depth_stratified"] = stratified_permutation_test(
            vals[pair][sub], hyp[pair][sub], strata[pair][sub], n_perm, PERM_SEED)
        print(f"    {label:<14} depth-stratified HYPOXIC vs VIABLE: "
              f"diff={block['hypoxic_vs_viable_depth_stratified']['observed']:+.5f} "
              f"p={block['hypoxic_vs_viable_depth_stratified']['p_value']:.5f}")
        enrichment[label] = block

    # ---- the all-480-gene ranking (confirm / refute the interrupted run) ---
    print("\n  all-480-gene Spearman ranking against the simulated oxygen ...")
    ranking = all_gene_ranking(expr, genes, totals, o2)
    rank_of = {g: i + 1 for i, (g, _, _) in enumerate(ranking)}
    rho_of = {g: r for g, r, _ in ranking}
    for g in hs["used"]:
        print(f"    {g}: rank {rank_of[g]}/{len(ranking)} most-negative, "
              f"rho={rho_of[g]:+.5f}")

    # Depth is a technical variable; how strongly does the field track it?
    o2_vs_depth = spearman_full(o2, totals)
    print(f"    (control) rho(simulated O2, total transcript count) = {o2_vs_depth:+.4f}")

    return {
        "sample": sample,
        "dataset": data["source"]["dataset"],
        "panel": data["source"]["panel"],
        "panel_name": panel_name,
        "panel_coverage": panel_cov,
        "n_cells": n_cells,
        "n_genes": int(len(genes)),
        "hypoxia_genes_used": hs["used"],
        "hypoxia_genes_missing": hs["missing"],
        "detection": hs["detection"],
        "gene_gene_rho": hs["gene_gene_rho"],
        "hypoxia_threshold": st["hypoxia_threshold"],
        "sweep": {dx: {k: v for k, v in sweep[dx].items() if k != "_state"}
                  for dx in dx_sweep},
        "perm": perm,
        "perm_subsample_n": int(idx.size),
        "n_perm": int(n_perm),
        "enrichment": enrichment,
        "ranking": ranking,
        "rank_of": {g: rank_of[g] for g in hs["used"]},
        "o2_vs_depth_rho": o2_vs_depth,
        "_ref": {
            "o2": o2, "score": score, "states": states, "tumor": tumor,
            "x": st["x"], "y": st["y"], "row": st["row"], "col": st["col"],
            "grid": st["grid"], "dx_um": DX_CONVERGED,
            "oxygen_field": st["oxygen_field"], "tissue_mask": st["tissue_mask"],
            "per_gene": hs["per_gene"], "totals": totals,
            "dist_endo": np.asarray(trivial_baselines(
                st["x"], st["y"], st["cell_type"], st["vascular_at_cell"]
            )["dist_endothelial"]["values"], float),
            "tiles": ref["tiles"],
        },
    }


# ---------------------------------------------------------------------------
# PANEL REGRESSION GUARD — the IO map must not have moved the brain-panel result.
# ---------------------------------------------------------------------------
#: The 319-gene brain-panel reconstruction as it stood BEFORE any panel-aware
#: marker map existed. :func:`check_brain_panel_unchanged` re-derives these and
#: fails loudly on any drift, so adding the IO panel cannot silently perturb the
#: section every other result in ``results/`` is computed on.
BRAIN_PANEL_REFERENCE_PCT: Dict[str, float] = {
    "TUMOR": 73.80, "STROMAL": 7.92, "IMMUNE": 8.54, "ENDOTHELIAL": 9.74,
}
BRAIN_PANEL_TOLERANCE_PCT = 0.02


def check_brain_panel_unchanged(sample: str = BRAIN_SAMPLE) -> Dict:
    """Assert the 319-gene reconstruction still gives the historical fractions.

    The IO marker map is only admissible if it is inert on the other section. This
    re-runs the real typing path on the brain panel and compares against
    :data:`BRAIN_PANEL_REFERENCE_PCT` to 0.02 percentage points, and additionally
    checks that :func:`marker_panels.select_marker_map` returns the *same object*
    as :data:`BRAIN_PANEL_MARKER_MAP` for that panel (identity, not just equality).
    """
    data = load_section(sample)
    genes = np.asarray(data["genes"], dtype=str)
    mm, panel_name, cov = select_marker_map(genes)
    st = band_states(sample, DX_CONVERGED, DEFAULT_MASK_MODE)
    ct = st["cell_type"]
    n = int(ct.size)
    pct = {c.name: 100.0 * float((ct == int(c)).sum()) / n for c in CellType}
    drift = {k: abs(pct[k] - v) for k, v in BRAIN_PANEL_REFERENCE_PCT.items()}
    ok = (panel_name == "brain_319" and mm is BRAIN_PANEL_MARKER_MAP
          and max(drift.values()) <= BRAIN_PANEL_TOLERANCE_PCT)
    print(f"\n  PANEL REGRESSION GUARD (319-gene brain section): "
          f"{'PASS' if ok else 'FAIL'}")
    print("    " + "  ".join(f"{k}={pct[k]:.2f}%" for k in BRAIN_PANEL_REFERENCE_PCT))
    print(f"    map identity is BRAIN_PANEL_MARKER_MAP: {mm is BRAIN_PANEL_MARKER_MAP}"
          f"   selected panel: {panel_name}  coverage {cov}")
    if not ok:
        raise RuntimeError(
            "the 319-gene brain-panel reconstruction CHANGED: "
            f"{pct} vs reference {BRAIN_PANEL_REFERENCE_PCT}"
        )
    return {"panel_name": panel_name, "coverage": cov, "percent": pct,
            "reference_percent": dict(BRAIN_PANEL_REFERENCE_PCT),
            "max_drift_pct": float(max(drift.values())),
            "is_same_map_object": bool(mm is BRAIN_PANEL_MARKER_MAP), "pass": ok}


# ---------------------------------------------------------------------------
# TEST 2 — cell-state model vs MEASURED proliferation on the 319-gene section.
# ---------------------------------------------------------------------------
def run_test2(sample: str = BRAIN_SAMPLE, n_perm: int = N_PERM_BRAIN) -> Dict:
    """Test 2 in full: measured proliferating fraction, state contrast, power, dx."""
    print("\n" + "=" * 78)
    print("TEST 2 — CELL-STATE MODEL vs MEASURED PROLIFERATION (319-gene section)")
    print("=" * 78)

    data = load_section(sample)
    expr = data["expression"]
    genes = np.asarray(data["genes"], dtype=str)
    totals = np.asarray(expr.sum(axis=1)).ravel().astype(np.float64)
    n_cells = int(expr.shape[0])
    print(f"  {n_cells:,} cells x {len(genes)} genes")

    # ---- the full mask x dx sweep of the model's state vector --------------
    sweep: Dict[str, Dict[float, Dict]] = {m: {} for m in TISSUE_MASK_MODES}
    for mode in TISSUE_MASK_MODES:
        for dx in DX_SWEEP:
            st = band_states(sample, dx, mode)
            if st["n_cells"] != n_cells:
                raise RuntimeError("row misalignment on the 319-gene section")
            sweep[mode][dx] = st
            f = st["fractions_tumor"]
            print(f"    mask={mode:<10} dx={dx:>2.0f}  mean O2={st['mean_oxygen_tissue']:.4f}  "
                  + "  ".join(f"{k}={100*f[k]:6.2f}%" for k in STATE_ORDER))

    mode0, dx0 = DOSE_SCHEDULE_CONFIG
    ref = sweep[mode0][dx0]
    tumor = ref["tumor_mask"]
    n_tumor = int(tumor.sum())
    states_t = ref["state"][tumor]
    o2_t = ref["local_o2"][tumor]
    totals_t = totals[tumor]
    log_depth = np.log(np.maximum(totals_t, 1.0))
    strata = depth_deciles(totals_t)
    print(f"\n  headline configuration = mask '{mode0}', dx={dx0:.0f} um "
          f"(the one results/dose_schedule.md ran under); {n_tumor:,} tumour cells")

    # ---- measured proliferating fraction, two documented rules -------------
    print("\n  MEASURED proliferating fraction of tumour cells:")
    markers: Dict[str, Dict] = {}
    for g in PROLIFERATION_MARKERS:
        raw = _column(expr, genes, g)
        if raw is None:
            markers[g] = {"on_panel": False}
            print(f"    {g}: NOT ON PANEL")
            continue
        rep = positivity_report(raw, totals, tumor)
        markers[g] = {"on_panel": True, "report": rep, "counts_tumor": raw[tumor]}
        print(f"    {g}: nonzero {100*rep['nonzero']['fraction']:.3f}% "
              f"({rep['nonzero']['n_positive']}/{rep['n_cells']})   "
              f"percentile-p{PERCENTILE_RULE_Q:.0f} "
              f"{100*rep['percentile']['fraction']:.3f}% "
              f"[{'DEGENERATE' if rep['percentile']['degenerate'] else 'ok'}]   "
              f"zeros {100*rep['zero_fraction']:.2f}%")

    lit_lo = min(float(e["value_pct"]) for e in LITERATURE_KI67)
    lit_hi = max(float(e["value_pct"]) for e in LITERATURE_KI67)
    mki = markers["MKI67"]["report"]["nonzero"]["fraction"] if markers["MKI67"]["on_panel"] else float("nan")
    print(f"    published GBM Ki-67 index {lit_lo:.0f}-{lit_hi:.0f}%  ->  the panel "
          f"detects MKI67 in {100*mki:.3f}% of tumour cells, "
          f"{lit_lo/(100*mki):.0f}-{lit_hi/(100*mki):.0f}x lower")

    # ---- model VIABLE fraction vs the measured fraction --------------------
    viable_range = [100 * sweep[m][dx]["fractions_tumor"]["VIABLE"]
                    for m in TISSUE_MASK_MODES for dx in DX_SWEEP]
    print(f"    model VIABLE fraction across the mask x dx grid: "
          f"{min(viable_range):.2f}% .. {max(viable_range):.2f}% "
          f"({max(viable_range)/max(min(viable_range), 1e-9):.0f}x swing)")

    # ---- state contrast + power -------------------------------------------
    print("\n  ARE model-VIABLE cells more proliferative than QUIESCENT/HYPOXIC?")
    viable = states_t == int(CellState.VIABLE)
    hypoxic = states_t == int(CellState.HYPOXIC)
    quiescent = states_t == int(CellState.QUIESCENT)

    contrasts: Dict[str, Dict] = {}
    for g in PROLIFERATION_MARKERS:
        if not markers[g]["on_panel"]:
            continue
        counts = markers[g]["counts_tumor"]
        pos = (counts >= 1).astype(np.float64)
        overall = float(pos.mean())
        norm = normalized_log1p(_column(expr, genes, g), totals)[tumor]
        rates = rate_by_state(pos, states_t)

        hq = hypoxic | quiescent
        block = {
            "rates": rates,
            "overall_rate": overall,
            "viable_vs_rest": {
                **permutation_test_means(pos, viable, n_perm, PERM_SEED),
                "cohens_d": cohens_d(norm[viable], norm[~viable]),
                "cliffs_delta": cliffs_delta(norm[viable], norm[~viable]),
                "power": min_detectable_rate_ratio(
                    overall, int(viable.sum()), int((~viable).sum())),
            },
            "hypoxic_vs_quiescent": {
                **permutation_test_means(pos[hq], hypoxic[hq], n_perm, PERM_SEED),
                "cohens_d": cohens_d(norm[hypoxic], norm[quiescent]),
                "cliffs_delta": cliffs_delta(norm[hypoxic], norm[quiescent]),
                "power": min_detectable_rate_ratio(
                    overall, int(hypoxic.sum()), int(quiescent.sum())),
            },
            "viable_vs_rest_depth_stratified": stratified_permutation_test(
                norm, viable, strata, n_perm, PERM_SEED),
            "hypoxic_vs_quiescent_depth_stratified": stratified_permutation_test(
                norm[hq], hypoxic[hq], strata[hq], n_perm, PERM_SEED),
            "o2_continuous": {
                **spearman_permutation(o2_t, norm, n_perm, PERM_SEED),
                "partial_rho_depth_removed": partial_spearman(o2_t, norm, log_depth),
            },
        }
        contrasts[g] = block
        r = rates
        print(f"    {g}+ by state: " + "  ".join(
            f"{k}={100*r[k]['rate']:.3f}% ({r[k]['k']}/{r[k]['n']:,})"
            for k in STATE_ORDER if r[k]["n"]))
        vv = block["viable_vs_rest"]
        print(f"      VIABLE vs rest        d={vv['cohens_d']:+.4f} p={vv['p_value']:.4f}"
              f" | POWER: {vv['power']['expected_positives']:.1f} expected positives, "
              f"needs a {vv['power']['rate_ratio']:.1f}x rate ratio to detect")
        hv = block["hypoxic_vs_quiescent"]
        print(f"      HYPOXIC vs QUIESCENT  d={hv['cohens_d']:+.4f} p={hv['p_value']:.4f}"
              f" | POWER: needs {hv['power']['rate_ratio']:.2f}x")
        print(f"      depth-stratified: VIABLE p="
              f"{block['viable_vs_rest_depth_stratified']['p_value']:.4f}  "
              f"HYP/QUI p={block['hypoxic_vs_quiescent_depth_stratified']['p_value']:.4f}")
        print(f"      continuous rho(simO2, {g}) = {block['o2_continuous']['rho']:+.4f} "
              f"p={block['o2_continuous']['p_value']:.4f}  "
              f"partial(depth removed)={block['o2_continuous']['partial_rho_depth_removed']:+.4f}")

    # ---- dx dependence of the same contrast --------------------------------
    print("\n  dx / mask dependence of the state contrast:")
    dx_dep: Dict[str, Dict[float, Dict]] = {m: {} for m in TISSUE_MASK_MODES}
    for mode in TISSUE_MASK_MODES:
        for dx in DX_SWEEP:
            s = sweep[mode][dx]
            tm = s["tumor_mask"]
            stt = s["state"][tm]
            row = {"fractions": s["fractions_tumor"],
                   "mean_oxygen_tissue": s["mean_oxygen_tissue"],
                   "n_tumor": int(tm.sum()), "grid": s["grid"]}
            for g in PROLIFERATION_MARKERS:
                if not markers[g]["on_panel"]:
                    continue
                pos = (markers[g]["counts_tumor"] >= 1).astype(np.float64)
                rr = rate_by_state(pos, stt)
                v = stt == int(CellState.VIABLE)
                row[g] = {
                    "rates": rr,
                    "viable_rate": rr["VIABLE"]["rate"],
                    "nonviable_rate": (float(pos[~v].mean()) if (~v).any() else float("nan")),
                    "power": min_detectable_rate_ratio(
                        float(pos.mean()), int(v.sum()), int((~v).sum())),
                }
            dx_dep[mode][dx] = row
            g0 = "MKI67"
            if g0 in row:
                print(f"    mask={mode:<10} dx={dx:>2.0f}  VIABLE={100*row['fractions']['VIABLE']:6.2f}%  "
                      f"{g0}+ in VIABLE={100*row[g0]['viable_rate']:.3f}%  "
                      f"elsewhere={100*row[g0]['nonviable_rate']:.3f}%  "
                      f"(needs {row[g0]['power']['rate_ratio']:.1f}x to detect)")

    return {
        "sample": sample,
        "dataset": data["source"]["dataset"],
        "panel": data["source"]["panel"],
        "n_cells": n_cells,
        "n_tumor": n_tumor,
        "headline_mask": mode0,
        "headline_dx": float(dx0),
        "hypoxia_threshold": ref["hypoxia_threshold"],
        "pd_config": {"hypoxia_resistance": PDConfig().hypoxia_resistance,
                      "quiescent_resistance": PDConfig().quiescent_resistance},
        "markers": {g: (markers[g]["report"] if markers[g]["on_panel"]
                        else {"on_panel": False}) for g in PROLIFERATION_MARKERS},
        "contrasts": contrasts,
        "dx_dependence": dx_dep,
        "viable_fraction_range_pct": [float(min(viable_range)), float(max(viable_range))],
        "literature_ki67_range_pct": [lit_lo, lit_hi],
        "n_perm": int(n_perm),
        "_ref": {
            "states_t": states_t, "o2_t": o2_t, "totals_t": totals_t,
            "x_t": ref["x"][tumor], "y_t": ref["y"][tumor],
            "counts": {g: markers[g]["counts_tumor"] for g in PROLIFERATION_MARKERS
                       if markers[g]["on_panel"]},
        },
    }


# ---------------------------------------------------------------------------
# VERDICTS. Mechanical rules, applied to the numbers — not written by hand.
# ---------------------------------------------------------------------------
#: Verdict vocabulary, fixed in advance so the conclusion cannot drift to fit the
#: result: SUPPORTED / PARTIALLY SUPPORTED / CONTRADICTED / UNDERPOWERED.
VERDICTS = ("SUPPORTED", "PARTIALLY SUPPORTED", "CONTRADICTED", "UNDERPOWERED")

#: A predictor must beat the best trivial baseline's |rho| by at least this
#: relative margin to count as "the physics is adding something". 5% is
#: deliberately generous to the model.
BASELINE_MARGIN = 0.05

#: Below this |rho| the correlation is treated as absent regardless of p — at
#: n = 816,769 a rho of 0.01 is "significant" and meaningless.
RHO_FLOOR = 0.02

#: A group contrast whose minimum detectable rate ratio exceeds this is declared
#: UNDERPOWERED: no plausible proliferation biology produces a 3x difference in
#: Ki-67 positivity between perfusion bands (Brat et al. measured 5-50% LESS,
#: i.e. a rate ratio between ~0.5 and ~0.95).
MAX_PLAUSIBLE_RATE_RATIO = 3.0


def verdict_test1(t1: Dict) -> Dict[str, object]:
    """Mechanical verdict for Test 1 from the measured numbers."""
    ref = t1["sweep"][DX_CONVERGED]["predictors"]
    sim = ref["simulated_o2"]["rho_score"]
    p = t1["perm"]["simulated_o2"]["p_value"]
    baselines = {k: v["rho_score"] for k, v in ref.items() if k != "simulated_o2"}
    best_key = max(baselines, key=lambda k: abs(baselines[k]))
    best = baselines[best_key]
    beats = abs(sim) > abs(best) * (1.0 + BASELINE_MARGIN)

    reasons: List[str] = []
    if abs(sim) < RHO_FLOOR:
        verdict = "CONTRADICTED"
        reasons.append(
            f"the simulated oxygen field has essentially no relationship to the "
            f"measured hypoxia response (|rho| = {abs(sim):.4f} < {RHO_FLOOR})")
    elif sim > 0:
        verdict = "CONTRADICTED"
        reasons.append(
            f"the correlation has the WRONG SIGN (rho = {sim:+.4f}); the model's "
            "low-oxygen regions express LESS hypoxia response, not more")
    elif not beats:
        verdict = "PARTIALLY SUPPORTED"
        reasons.append(
            f"the sign and significance are right (rho = {sim:+.4f}, p = {p:.5f}) "
            f"but a trivial physics-free baseline ({best_key}, rho = {best:+.4f}) "
            f"matches or BEATS the solved field, so the reaction-diffusion solve "
            "is not what is carrying the signal")
    else:
        verdict = "SUPPORTED"
        reasons.append(
            f"rho = {sim:+.4f} (p = {p:.5f}) and the solved field beats every "
            f"trivial baseline (best: {best_key} at {best:+.4f})")

    # Per-gene split — the interrupted run's claim, adjudicated.
    per_gene = {g: ref["simulated_o2"].get(f"rho_{g}") for g in t1["hypoxia_genes_used"]}
    null_genes = [g for g, r in per_gene.items()
                  if r is not None and abs(r) < RHO_FLOOR]
    if null_genes:
        reasons.append(
            "the composite is carried by ONE gene: "
            + ", ".join(f"{g} rho={per_gene[g]:+.4f}" for g in per_gene)
            + f" — {', '.join(null_genes)} is null")
    return {
        "verdict": verdict,
        "reasons": reasons,
        "rho_simulated": sim,
        "p_value": p,
        "best_baseline": best_key,
        "best_baseline_rho": best,
        "beats_baselines": bool(beats),
        "per_gene_rho": per_gene,
    }


def verdict_test2(t2: Dict) -> Dict[str, object]:
    """Mechanical verdict for Test 2, with UNDERPOWERED as a first-class outcome."""
    reasons: List[str] = []
    underpowered: List[str] = []
    wrong_direction: List[str] = []
    supported: List[str] = []

    for g, block in t2["contrasts"].items():
        vv = block["viable_vs_rest"]
        rr = vv["power"]["rate_ratio"]
        exp = vv["power"]["expected_positives"]
        if rr > MAX_PLAUSIBLE_RATE_RATIO or exp < 5.0:
            underpowered.append(
                f"{g} VIABLE-vs-rest could only have detected a {rr:.1f}x rate "
                f"ratio ({exp:.1f} positives expected in VIABLE under the null)")
        hq = block["hypoxic_vs_quiescent"]
        rates = block["rates"]
        if rates["HYPOXIC"]["n"] and rates["QUIESCENT"]["n"]:
            if rates["HYPOXIC"]["rate"] > rates["QUIESCENT"]["rate"]:
                wrong_direction.append(
                    f"{g}+ is HIGHER in model-HYPOXIC ({100*rates['HYPOXIC']['rate']:.3f}%) "
                    f"than in model-QUIESCENT ({100*rates['QUIESCENT']['rate']:.3f}%), "
                    "the opposite of the measured biology")
            elif hq["p_value"] < 0.05:
                supported.append(f"{g} HYPOXIC < QUIESCENT at p={hq['p_value']:.4f}")

    if underpowered and not supported:
        verdict = "UNDERPOWERED"
        reasons.extend(underpowered)
        reasons.append(
            "the decisive VIABLE-vs-rest comparison cannot distinguish 'the state "
            "labelling is right' from 'the state labelling is noise'")
    elif wrong_direction and not supported:
        verdict = "CONTRADICTED"
        reasons.extend(wrong_direction)
    elif supported and not wrong_direction:
        verdict = "SUPPORTED"
        reasons.extend(supported)
    else:
        verdict = "PARTIALLY SUPPORTED"
        reasons.extend(supported + wrong_direction + underpowered)

    lo, hi = t2["viable_fraction_range_pct"]
    reasons.append(
        f"independently of power, the model's VIABLE fraction is not a measurement: "
        f"it runs {lo:.2f}% to {hi:.2f}% across two tissue-mask conventions and "
        "three grid spacings, all of which are bookkeeping choices")
    return {
        "verdict": verdict,
        "reasons": reasons,
        "underpowered_contrasts": underpowered,
        "wrong_direction": wrong_direction,
        "supported_contrasts": supported,
    }


# ---------------------------------------------------------------------------
# The experiment.
# ---------------------------------------------------------------------------
def run(
    io_sample: str = IO_SAMPLE, brain_sample: str = BRAIN_SAMPLE,
    n_perm_io: int = N_PERM_IO, n_perm_brain: int = N_PERM_BRAIN,
) -> Dict:
    """Run both validations and return every number the report quotes."""
    print("=== COMBINED STATE VALIDATION: simulated O2 and the cell-state model ===")
    guard = check_brain_panel_unchanged(brain_sample)
    t1 = run_test1(io_sample, n_perm_io)
    t2 = run_test2(brain_sample, n_perm_brain)
    v1 = verdict_test1(t1)
    v2 = verdict_test2(t2)
    print(f"\nVERDICT Test 1 (simulated oxygen field): {v1['verdict']}")
    for r in v1["reasons"]:
        print(f"  - {r}")
    print(f"\nVERDICT Test 2 (cell-state labelling): {v2['verdict']}")
    for r in v2["reasons"]:
        print(f"  - {r}")
    return {"guard": guard, "test1": t1, "test2": t2,
            "verdict1": v1, "verdict2": v2}
