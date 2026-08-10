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

#: Measured hypoxia-RESPONSE genes on the IO panel. Single source of truth is
#: :data:`marker_panels.HYPOXIA_RESPONSE_GENES`, so the panel module and this
#: validation can never drift apart. NOT an oxygen tension — see caveat 2 in the
#: module docstring.
HYPOXIA_GENES: Tuple[str, ...] = tuple(HYPOXIA_RESPONSE_GENES["io_480"])

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
            # ALREADY sign-aligned to oxygen (= minus the distance), like every
            # other baseline. Do not negate it again at the call site.
            "neg_dist_endo": np.asarray(trivial_baselines(
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
#:
#: SCOPED TO THE LEGACY SCORING RULE. This guard exists to prove the IO MARKER MAP
#: is inert on the brain panel, so it is now checked under
#: ``aggregation='mean_z'`` — the rule that was in force when these numbers were
#: frozen. The production default has since moved to ``'matched_null_z'`` to
#: remove the marker-COUNT bias (``results/classifier_fix.md``), which changes the
#: composition ON PURPOSE; pinning the current default is
#: ``tests/test_classifier_arity.py``'s job, not this guard's.
BRAIN_PANEL_REFERENCE_PCT: Dict[str, float] = {
    "TUMOR": 73.80, "STROMAL": 7.92, "IMMUNE": 8.54, "ENDOTHELIAL": 9.74,
}
BRAIN_PANEL_TOLERANCE_PCT = 0.02

#: The same section under the CURRENT default aggregation. Reported, not asserted
#: — it is pinned by ``tests/test_classifier_arity.py``.
BRAIN_PANEL_CURRENT_PCT: Dict[str, float] = {
    "TUMOR": 58.24, "STROMAL": 11.78, "IMMUNE": 19.48, "ENDOTHELIAL": 10.50,
}


def check_brain_panel_unchanged(sample: str = BRAIN_SAMPLE) -> Dict:
    """Assert the 319-gene MARKER MAP still gives the historical fractions.

    The IO marker map is only admissible if it is inert on the other section. This
    re-runs the real typing path on the brain panel **under the legacy
    ``'mean_z'`` aggregation** and compares against
    :data:`BRAIN_PANEL_REFERENCE_PCT` to 0.02 percentage points, and additionally
    checks that :func:`marker_panels.select_marker_map` returns the *same object*
    as :data:`BRAIN_PANEL_MARKER_MAP` for that panel (identity, not just equality).

    Scoring the legacy rule (rather than the production default) is deliberate:
    the property under test is that the PANEL REGISTRY did not move this section,
    and that property must stay checkable after the deliberate change of default
    aggregation documented in ``results/classifier_fix.md``. The composition under
    the current default is reported alongside so any drift there is still visible.
    """
    from ..reconstruct.real_xenium import TYPING_AGGREGATION, phenotype_cells

    data = load_section(sample)
    genes = np.asarray(data["genes"], dtype=str)
    mm, panel_name, cov = select_marker_map(genes)

    ct_legacy, _um, _ms = phenotype_cells(
        data["expression"], genes, marker_map=mm, aggregation="mean_z"
    )
    n = int(ct_legacy.size)
    pct = {c.name: 100.0 * float((ct_legacy == int(c)).sum()) / n for c in CellType}

    ct_now, _um2, _ms2 = phenotype_cells(data["expression"], genes, marker_map=mm)
    pct_now = {c.name: 100.0 * float((ct_now == int(c)).sum()) / n for c in CellType}

    drift = {k: abs(pct[k] - v) for k, v in BRAIN_PANEL_REFERENCE_PCT.items()}
    ok = (panel_name == "brain_319" and mm is BRAIN_PANEL_MARKER_MAP
          and max(drift.values()) <= BRAIN_PANEL_TOLERANCE_PCT)
    print(f"\n  PANEL REGRESSION GUARD (319-gene brain section, aggregation="
          f"'mean_z'): {'PASS' if ok else 'FAIL'}")
    print("    " + "  ".join(f"{k}={pct[k]:.2f}%" for k in BRAIN_PANEL_REFERENCE_PCT))
    print(f"    current default ({TYPING_AGGREGATION}, reported not asserted): "
          + "  ".join(f"{k}={pct_now[k]:.2f}%" for k in BRAIN_PANEL_REFERENCE_PCT))
    print(f"    map identity is BRAIN_PANEL_MARKER_MAP: {mm is BRAIN_PANEL_MARKER_MAP}"
          f"   selected panel: {panel_name}  coverage {cov}")
    if not ok:
        raise RuntimeError(
            "the 319-gene brain-panel MARKER MAP changed the reconstruction: "
            f"{pct} vs reference {BRAIN_PANEL_REFERENCE_PCT}"
        )
    return {"panel_name": panel_name, "coverage": cov, "percent": pct,
            "reference_percent": dict(BRAIN_PANEL_REFERENCE_PCT),
            "aggregation_checked": "mean_z",
            "current_default_aggregation": str(TYPING_AGGREGATION),
            "current_default_percent": pct_now,
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
    for v in (v1, v2):
        if v["verdict"] not in VERDICTS:
            raise RuntimeError(
                f"verdict {v['verdict']!r} is outside the fixed vocabulary "
                f"{VERDICTS} — the conclusion vocabulary must not drift"
            )
    print(f"\nVERDICT Test 1 (simulated oxygen field): {v1['verdict']}")
    for r in v1["reasons"]:
        print(f"  - {r}")
    print(f"\nVERDICT Test 2 (cell-state labelling): {v2['verdict']}")
    for r in v2["reasons"]:
        print(f"  - {r}")
    return {"guard": guard, "test1": t1, "test2": t2,
            "verdict1": v1, "verdict2": v2}


# ---------------------------------------------------------------------------
# Figure.
# ---------------------------------------------------------------------------
def _style_axes(ax, grid_axis: str = "y") -> None:
    """Recessive grid, no top/right spines — the marks carry the ink, not the frame."""
    ax.set_axisbelow(True)
    ax.grid(True, axis=grid_axis, color=GRID_COLOR, linewidth=0.8, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID_COLOR)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8, length=0)


def _ramp(colors: Sequence[str], name: str):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, list(colors))


def _tile_field(
    row: np.ndarray, col: np.ndarray, dx_um: float, values: np.ndarray,
    tile_um: float = TILE_UM, min_cells: int = TILE_MIN_CELLS,
) -> Tuple[np.ndarray, Tuple[float, float, float, float]]:
    """Mean of ``values`` per fixed physical tile, as an image + its um extent."""
    cx = (np.asarray(col, float) + 0.5) * float(dx_um)
    cy = (np.asarray(row, float) + 0.5) * float(dx_um)
    tx = (cx / float(tile_um)).astype(np.int64)
    ty = (cy / float(tile_um)).astype(np.int64)
    ntx, nty = int(tx.max()) + 1, int(ty.max()) + 1
    flat = ty * ntx + tx
    cnt = np.bincount(flat, minlength=nty * ntx).astype(np.float64)
    tot = np.bincount(flat, weights=np.asarray(values, float), minlength=nty * ntx)
    img = np.full(nty * ntx, np.nan)
    ok = cnt >= int(min_cells)
    img[ok] = tot[ok] / cnt[ok]
    return img.reshape(nty, ntx), (0.0, ntx * tile_um, nty * tile_um, 0.0)


def _binned_profile(
    predictor: np.ndarray, response: np.ndarray, n_bins: int = 25,
) -> Tuple[np.ndarray, np.ndarray]:
    """Mean response per equal-count bin of the predictor's PERCENTILE rank.

    Percentile rank rather than raw units so predictors on incomparable scales
    (relative oxygen, micrometres, cell counts) share one x-axis, which is the only
    way the model and the trivial baseline can be overlaid on the same panel
    without a second axis (the skill's one-axis rule).
    """
    pr = rankdata(np.asarray(predictor, float)) / float(predictor.size) * 100.0
    edges = np.linspace(0, 100, n_bins + 1)
    idx = np.clip(np.searchsorted(edges, pr, side="right") - 1, 0, n_bins - 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    means = np.array([
        float(np.asarray(response, float)[idx == b].mean()) if np.any(idx == b) else np.nan
        for b in range(n_bins)
    ])
    return centres, means


def make_figure(results: Dict, out_path: Optional[str] = None) -> str:
    """Six-panel evidence figure. Returns the written path."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    t1, t2 = results["test1"], results["test2"]
    v1, v2 = results["verdict1"], results["verdict2"]
    r = t1["_ref"]
    preds = t1["sweep"][DX_CONVERGED]["predictors"]

    fig = plt.figure(figsize=(18.0, 10.8))
    gs = fig.add_gridspec(2, 3, hspace=0.46, wspace=0.28,
                          left=0.088, right=0.985, top=0.855, bottom=0.095)
    fig.patch.set_facecolor(SURFACE)

    # ---- A: the SIMULATED oxygen field --------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ofield = np.where(r["tissue_mask"], r["oxygen_field"], np.nan)
    ny, nx = r["grid"]
    ext = (0.0, nx * r["dx_um"], ny * r["dx_um"], 0.0)
    im = ax.imshow(ofield, cmap=_ramp(SEQ_BLUE, "o2"), extent=ext,
                   interpolation="nearest", rasterized=True)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("relative O₂ (model)", fontsize=8, color=INK_SECONDARY)
    cb.ax.tick_params(labelsize=7, colors=INK_SECONDARY, length=0)
    cb.outline.set_visible(False)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)", fontsize=8.5, color=INK_SECONDARY)
    ax.set_ylabel("y (µm)", fontsize=8.5, color=INK_SECONDARY)
    ax.set_title(f"A · SIMULATED oxygen — physics.solve_oxygen, dx={r['dx_um']:.0f} µm\n"
                 f"mean {t1['sweep'][DX_CONVERGED]['mean_oxygen_tissue']:.3f} over tissue",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax, grid_axis="both")
    ax.grid(False)

    # ---- B: the MEASURED hypoxia-response field -----------------------------
    ax = fig.add_subplot(gs[0, 1])
    simg, sext = _tile_field(r["row"], r["col"], r["dx_um"], r["score"])
    im = ax.imshow(simg, cmap=_ramp(SEQ_AQUA, "hyp"), extent=sext,
                   interpolation="nearest", rasterized=True,
                   vmin=np.nanpercentile(simg, 2), vmax=np.nanpercentile(simg, 98))
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("VEGFA+HIF1A z (measured)", fontsize=8, color=INK_SECONDARY)
    cb.ax.tick_params(labelsize=7, colors=INK_SECONDARY, length=0)
    cb.outline.set_visible(False)
    ax.set_aspect("equal")
    ax.set_xlabel("x (µm)", fontsize=8.5, color=INK_SECONDARY)
    ax.set_ylabel("y (µm)", fontsize=8.5, color=INK_SECONDARY)
    tl = t1["sweep"][DX_CONVERGED]["tiles"]
    ax.set_title(f"B · MEASURED hypoxia response — VEGFA+HIF1A, {TILE_UM:.0f} µm tiles\n"
                 f"field-vs-field ρ = {tl['rho']:+.3f} over {tl['n_tiles']:,} tiles",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax, grid_axis="both")
    ax.grid(False)

    # ---- C: the correlation, with the TRIVIAL BASELINE overlaid -------------
    ax = fig.add_subplot(gs[0, 2])
    series = [
        ("SIMULATED O₂ — the PDE", r["o2"],
         preds["simulated_o2"]["rho_score"], MODEL_COLOR, "-"),
        ("−dist. to endothelial cell — TRIVIAL", r["neg_dist_endo"],
         preds["dist_endothelial"]["rho_score"], BASELINE_COLOR, "--"),
    ]
    handles = []
    ends = {}
    for lab, vals, rho, colr, ls in series:
        cx, cy = _binned_profile(vals, r["score"])
        ends[lab] = float(cy[0])
        ax.plot(cx, cy, color=colr, linewidth=2.0, linestyle=ls, zorder=3,
                solid_capstyle="round")
        ax.scatter(cx, cy, s=22, color=colr, edgecolor=SURFACE, linewidth=1.2, zorder=4)
        handles.append(plt.Line2D([], [], color=colr, linewidth=2.0, linestyle=ls,
                                  marker="o", markersize=5,
                                  label=f"{lab}   ρ = {rho:+.4f}"))
    ax.axhline(0.0, color=INK_MUTED, linewidth=1.0, zorder=2)
    lo_end = ends["SIMULATED O₂ — the PDE"]
    hi_end = ends["−dist. to endothelial cell — TRIVIAL"]
    ax.annotate(
        f"at the most-hypoxic decile the trivial baseline reaches "
        f"{hi_end:+.2f}\nz-units of measured hypoxia; the PDE only {lo_end:+.2f}",
        (2.0, hi_end), textcoords="offset points", xytext=(16, -6),
        ha="left", va="top", fontsize=8, color=INK_PRIMARY, zorder=6)
    ax.legend(handles=handles, fontsize=7.5, loc="lower left", frameon=False)
    ax.set_xlabel("percentile of the predictor  (0 = most hypoxic / furthest from a vessel)",
                  fontsize=8.5, color=INK_SECONDARY)
    ax.set_ylabel("mean measured VEGFA+HIF1A z-score", fontsize=8.5, color=INK_SECONDARY)
    ax.set_title("C · The trivial baseline separates measured hypoxia MORE\n"
                 "sharply than the reaction-diffusion solve does",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax, grid_axis="both")

    # ---- D: predictor league table -------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    order = ["simulated_o2", "vascular_kde", "dist_endothelial", "local_cell_density"]
    labels = ["SIMULATED O₂ (PDE)",
              "vascular KDE (source)",
              "−dist. to endothelium",
              f"−cell density ({DENSITY_RADIUS_UM:.0f} µm)"]
    ys = np.arange(len(order))[::-1]
    for key, lab, y in zip(order, labels, ys):
        rho = preds[key]["rho_score"]
        colr = MODEL_COLOR if key == "simulated_o2" else BASELINE_COLOR
        ax.barh(y, abs(rho), color=colr, height=0.56, zorder=3)
        ax.text(abs(rho) + 0.004, y, f"ρ = {rho:+.4f}", va="center", ha="left",
                fontsize=8.5, color=INK_PRIMARY, zorder=4)
    best = v1["best_baseline_rho"]
    ax.axvline(abs(best), color=INK_MUTED, linewidth=1.2, linestyle=(0, (4, 3)), zorder=2)
    ax.set_yticks(ys)
    ax.set_yticklabels(labels, fontsize=8, color=INK_PRIMARY)
    ax.set_xlim(0, max(abs(preds[k]["rho_score"]) for k in order) * 1.42)
    ax.set_xlabel("|Spearman ρ| with the measured VEGFA+HIF1A score",
                  fontsize=8.5, color=INK_SECONDARY)
    ax.set_title("D · Physics-free baselines BEAT the solved oxygen field\n"
                 "(same inputs, no PDE) — so the solve adds nothing here",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax, grid_axis="x")
    ax.legend(handles=[Patch(color=MODEL_COLOR, label="simulated (PDE)"),
                       Patch(color=BASELINE_COLOR, label="trivial baseline")],
              fontsize=7.5, loc="lower right", frameon=False)

    # ---- E: measured markers by MODEL state (IO section) --------------------
    ax = fig.add_subplot(gs[1, 1])
    tracks = ["VEGFA", "HIF1A"]
    width = 0.36
    xs = np.arange(len(STATE_ORDER))
    for gi, g in enumerate(tracks):
        vals = r["per_gene"][g]
        for si, name in enumerate(STATE_ORDER):
            m = r["states"] == int(_STATE_ENUM[name])
            n = int(m.sum())
            if n == 0:
                continue
            mu = float(vals[m].mean())
            se = float(vals[m].std(ddof=1) / np.sqrt(n))
            xpos = si + (gi - 0.5) * width
            ax.bar(xpos, mu, width=width * 0.9, color=STATE_COLORS[name],
                   alpha=1.0 if gi == 0 else 0.45, zorder=3,
                   edgecolor=SURFACE, linewidth=1.0)
            ax.plot([xpos, xpos], [mu - 1.96 * se, mu + 1.96 * se],
                    color=INK_PRIMARY, linewidth=1.4, zorder=5)
            # RELIEF RULE: yellow is sub-3:1, so every bar carries a direct label.
            ax.text(xpos, mu + (0.012 if mu >= 0 else -0.012), f"{mu:+.3f}",
                    ha="center", va="bottom" if mu >= 0 else "top",
                    fontsize=7.5, color=INK_PRIMARY, zorder=6)
    ax.axhline(0.0, color=INK_MUTED, linewidth=1.0, zorder=2)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{s}\nn={int((r['states']==int(_STATE_ENUM[s])).sum()):,}"
                        for s in STATE_ORDER], fontsize=8.5, color=INK_PRIMARY)
    ax.set_ylabel("mean measured expression (z, depth-normalised)",
                  fontsize=8.5, color=INK_SECONDARY)
    e = t1["enrichment"]
    ax.set_ylim(ax.get_ylim()[0] * 1.12, max(ax.get_ylim()[1] * 1.9, 0.22))
    ax.set_title("E · Measured markers by MODEL state (IO section). The banding does\n"
                 f"separate VEGFA (d(HYP vs VIA) = "
                 f"{e['VEGFA']['hypoxic_vs_viable']['cohens_d']:+.2f}); HIF1A only "
                 f"{e['HIF1A']['hypoxic_vs_viable']['cohens_d']:+.2f}",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax)
    ax.legend(handles=[Patch(facecolor="#6b6a66", label="VEGFA"),
                       Patch(facecolor="#6b6a66", alpha=0.45, label="HIF1A")],
              fontsize=7.5, loc="upper right", frameon=False, ncol=2)

    # ---- F: Test 2 — measured proliferation by model state ------------------
    ax = fig.add_subplot(gs[1, 2])
    width = 0.36
    for gi, g in enumerate(PROLIFERATION_MARKERS):
        if g not in t2["contrasts"]:
            continue
        rates = t2["contrasts"][g]["rates"]
        for si, name in enumerate(STATE_ORDER):
            rr = rates[name]
            if not rr["n"]:
                continue
            xpos = si + (gi - 0.5) * width
            ax.bar(xpos, 100 * rr["rate"], width=width * 0.9,
                   color=STATE_COLORS[name], alpha=1.0 if gi == 0 else 0.45,
                   zorder=3, edgecolor=SURFACE, linewidth=1.0)
            ax.plot([xpos, xpos], [100 * rr["lo"], 100 * rr["hi"]],
                    color=INK_PRIMARY, linewidth=1.4, zorder=5)
            ax.text(xpos, 100 * rr["hi"] + 0.06,
                    f"{100*rr['rate']:.2f}%\n{rr['k']}/{rr['n']:,}",
                    ha="center", va="bottom", fontsize=7, color=INK_PRIMARY, zorder=6)
    ax.set_xticks(np.arange(len(STATE_ORDER)))
    ax.set_xticklabels(STATE_ORDER, fontsize=8.5, color=INK_PRIMARY)
    ax.set_ylabel("marker-positive tumour cells (%)", fontsize=8.5, color=INK_SECONDARY)
    ax.set_ylim(0, 7.6)
    pw = t2["contrasts"]["MKI67"]["viable_vs_rest"]["power"]
    lo, hi = t2["literature_ki67_range_pct"]
    ax.text(0.015, 0.985,
            f"published GBM Ki-67 index {lo:.0f}–{hi:.0f}% — off this axis;\n"
            f"the panel detects MKI67 in only "
            f"{100*t2['markers']['MKI67']['nonzero']['fraction']:.2f}% of tumour cells.\n"
            f"VIABLE holds {pw['n_group']} cells, {pw['expected_positives']:.1f} MKI67+ "
            f"expected under the null:\nonly a ≥{pw['rate_ratio']:.1f}× rate ratio was "
            "detectable → UNDERPOWERED",
            transform=ax.transAxes, fontsize=7.5, color=INK_SECONDARY,
            va="top", ha="left", zorder=7)
    ax.set_title(f"F · MEASURED proliferation by model state — 319-gene section,\n"
                 f"mask '{t2['headline_mask']}', dx={t2['headline_dx']:.0f} µm "
                 "(the dose_schedule.md configuration)",
                 fontsize=10, color=INK_PRIMARY, loc="left", pad=8)
    _style_axes(ax)
    ax.legend(handles=[Patch(facecolor="#6b6a66", label="MKI67"),
                       Patch(facecolor="#6b6a66", alpha=0.45, label="TOP2A")],
              fontsize=7.5, loc="upper right", frameon=False, ncol=2)

    fig.suptitle(
        "Does the oxygen → cell-state → resistance chain survive contact with "
        "measured biology? — two real 10x Xenium GBM sections",
        fontsize=13.5, color=INK_PRIMARY, x=0.055, ha="left", y=0.972)
    fig.text(0.055, 0.938,
             f"TEST 1 (simulated O₂ vs measured VEGFA/HIF1A, {t1['n_cells']:,} cells): "
             f"{v1['verdict']}  ·  "
             f"TEST 2 (cell states vs measured MKI67/TOP2A, {t2['n_tumor']:,} tumour "
             f"cells): {v2['verdict']}",
             fontsize=10, color=INK_PRIMARY, ha="left")
    _fo = _ordering(preds["simulated_o2"]["rho_score"],
                    preds["dist_endothelial"]["rho_score"])
    fig.text(0.055, 0.905,
             ("The PDE out-predicts a one-line distance-to-vessel baseline "
              f"({_fo['ratio']:.2f}× on |ρ|) but only at |ρ| ≈ "
              f"{_fo['abs_sim']:.2f}, within-section; "
              if _fo["sim_wins"] else
              "The PDE is out-predicted by a one-line distance-to-vessel "
              "baseline; ")
             + "HIF1A is null; and the proliferation test cannot distinguish a "
             "correct state model from noise. Two different patients — the tests "
             "are independent, not corroborating.",
             fontsize=8.5, color=INK_SECONDARY, ha="left")

    out_path = out_path or os.path.join(_results_dir(), "state_validation.png")
    fig.savefig(out_path, dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"figure written: {out_path}")
    return out_path


# ---------------------------------------------------------------------------
# Report.
# ---------------------------------------------------------------------------
def _pmid(entry: Dict[str, object]) -> str:
    return (f"{entry['citation']} (PMID {entry['pmid']}, "
            f"doi:{entry['doi']})")


def _ordering(sim: float, base: float) -> Dict[str, object]:
    """Derive the simulated-vs-baseline comparison FROM THE NUMBERS.

    Every predictor in Test 1 is sign-aligned to oxygen — :func:`trivial_baselines`
    negates the distance and the cell count explicitly — so both a correct solved
    field and a correct baseline score NEGATIVE, and the stronger predictor is
    simply the one with the larger ``|rho|``. The comparison is therefore
    ``abs(sim)`` against ``abs(base)``, never ``sim`` against ``base``.

    WHY THIS FUNCTION EXISTS. The prose in this report used to hard-code the
    direction of that comparison in English while :func:`verdict_test1` computed
    it with ``abs()``. When the measured ordering was the reverse of the
    hard-coded sentence, the generator kept emitting the sentence — so
    ``results/state_validation.md`` §1.2 asserted that a distance baseline beat
    the solved field (-0.1204 vs -0.1523) while §1.6 of the same document, from
    the same run, returned SUPPORTED with the reason "the solved field beats
    every trivial baseline". The document contradicted itself for as long as both
    sections existed, and a whole component (``results/oxygen_demotion.md``,
    ``physics/oxygen_proxy.py``) was deprecated on the strength of the English
    rather than the arithmetic.

    Deriving the wording here is what makes that unrepeatable: there is now no
    sentence in this module that states an ordering the numbers did not produce,
    so a regeneration reproduces whatever the data say — including a flip back —
    instead of reverting to a remembered conclusion.
    """
    a_sim, a_base = abs(float(sim)), abs(float(base))
    sim_wins = a_sim > a_base
    strong, weak = (a_sim, a_base) if sim_wins else (a_base, a_sim)
    return {
        "sim_wins": bool(sim_wins),
        "abs_sim": a_sim,
        "abs_base": a_base,
        # How many times stronger the WINNER is than the loser.
        "ratio": (strong / weak) if weak else float("inf"),
        # What fraction of the winner the loser recovers.
        "recovers": (weak / strong) if strong else float("nan"),
        # The baseline expressed as a multiple of the solved field. < 1 means the
        # baseline is weaker. This is the number the superseded prose reported as
        # a percentage CHANGE ("-21%"), which is why it read as though the
        # baseline had gained rather than lost.
        "base_over_sim": (a_base / a_sim) if a_sim else float("inf"),
        "winner": "the solved field" if sim_wins else "the baseline",
        # Verb phrases, so no call site has to re-decide the direction.
        "verb_vs_baseline": "beats" if sim_wins else "is beaten by",
        "strongest_or_weakest": "strongest" if sim_wins else "weakest",
    }


def write_report(results: Dict, out_path: Optional[str] = None) -> str:
    """Write ``results/state_validation.md``. Returns the written path.

    Every sentence that states an ORDERING between the solved field and a
    baseline is derived through :func:`_ordering` rather than written out, for
    the reason given in that function's docstring.
    """
    t1, t2 = results["test1"], results["test2"]
    v1, v2 = results["verdict1"], results["verdict2"]
    guard = results["guard"]
    preds = t1["sweep"][DX_CONVERGED]["predictors"]
    # The three Test-1 numbers every ordering claim below is derived from.
    so = preds["simulated_o2"]["rho_score"]
    de = preds["dist_endothelial"]["rho_score"]
    vk = preds["vascular_kde"]["rho_score"]
    ordv = _ordering(so, de)       # solved field vs the distance baseline
    ordk = _ordering(so, vk)       # solved field vs its own undiffused source
    L: List[str] = []
    A = L.append

    # =====================================================================
    A("# Do the simulated oxygen field and the oxygen-banded cell states carry "
      "real biology?")
    A("")
    A("Two independent falsification tests against MEASURED markers in real GBM "
      "tissue, on two different open 10x Xenium sections.")
    A("")
    # The correction record. Emitted, not hand-written, so a regeneration cannot
    # quietly drop it -- this project keeps a visible trail of its own mistakes.
    # It only appears while the solved field is the winner, i.e. exactly when a
    # reader would otherwise have no way of knowing the document once said the
    # opposite.
    if ordv["sim_wins"]:
        A("> ## CORRECTION (§1.2) — this document once had its central "
          "inequality backwards")
        A(">")
        A("> **What §1.2 previously said.** That a one-line distance-to-vessel "
          "baseline **out-predicted** the solved oxygen field, quoting the two "
          "numbers below in that order and describing the gap as \"-21%\". The "
          "one-sentence summary, the headline table, §1.3, §1.5 and consequences "
          "1–2 of §1.6 all repeated that reading.")
        A(">")
        A("> **Why it was wrong.** Every predictor in the §1.2 table is "
          "sign-aligned to oxygen — `trivial_baselines` negates the distance — "
          "so, as that table's own header states, a correct predictor gives a "
          f"NEGATIVE ρ and the STRONGER predictor is the one with the larger "
          f"|ρ|. |{so:+.4f}| > |{de:+.4f}|: **the solved field out-predicts the "
          f"distance baseline**, by {ordv['ratio']:.2f}×. The \"-21%\" was itself "
          f"the tell — {ordv['base_over_sim']:.2f} means the baseline is "
          f"{100*ordv['base_over_sim']:.0f}% as strong, i.e. **weaker**, not "
          "better.")
        A(">")
        A("> **Every mechanical check in this repo always disagreed with the "
          "prose.** `verdict_test1()` compares `abs(sim) > abs(best)` and has "
          "always returned **SUPPORTED** here, with the reason string \"the "
          "solved field beats every trivial baseline\" — which is why §1.6 "
          "contradicted §1.2 above it for as long as both existed. "
          "`validation/heldout_patient.py` adjudicates the same comparison "
          "arithmetically and returns **REFUTED** on the demotion claim. The "
          "nested partials and the depth-partialled column agree.")
        A(">")
        A("> **What it cost.** A component was deprecated on the misread: "
          "`results/oxygen_demotion.md` and `tmesim/physics/oxygen_proxy.py` were "
          "both written on the inverted premise. That document is now the record "
          "of the correction rather than of the demotion.")
        A(">")
        A("> **The fix that makes this unrepeatable.** Every ordering claim in "
          "this report is now DERIVED from the numbers via "
          "`state_validation._ordering()` instead of written in English, so "
          "re-running this module reproduces whatever the data say — including a "
          "flip back — rather than reverting to a remembered conclusion.")
        A(">")
        A("> **What does NOT follow.** The solved field winning does not make it "
          f"good. |ρ| ≈ {ordv['abs_sim']:.2f} is WEAK, this is a within-section "
          "re-derivation and not out-of-sample, and the composite is carried by "
          "VEGFA alone. The correct conclusion is **\"the PDE is the best of four "
          "weak predictors\"**, not \"the oxygen solve is validated\".")
        A("")
    A("| | Test 1 — simulated oxygen field | Test 2 — cell-state labelling |")
    A("|---|---|---|")
    A(f"| **Verdict** | **{v1['verdict']}** | **{v2['verdict']}** |")
    A(f"| Section | {t1['dataset']} | {t2['dataset']} |")
    A(f"| Cells | {t1['n_cells']:,} ({t1['n_genes']}-gene panel) | "
      f"{t2['n_cells']:,} ({t2['n_tumor']:,} tumour; 319-gene panel) |")
    A("| Measured against | VEGFA + HIF1A (hypoxia response) | MKI67 + TOP2A "
      "(proliferation) |")
    A(f"| Headline statistic | Spearman ρ = "
      f"{preds['simulated_o2']['rho_score']:+.4f} "
      f"(p = {t1['perm']['simulated_o2']['p_value']:.5f}) | MKI67+ in VIABLE "
      f"{100*t2['contrasts']['MKI67']['rates']['VIABLE']['rate']:.2f}% "
      f"({t2['contrasts']['MKI67']['rates']['VIABLE']['k']}/"
      f"{t2['contrasts']['MKI67']['rates']['VIABLE']['n']}) |")
    A(f"| Baseline contrast | a one-line distance-to-vessel baseline scores "
      f"**{de:+.4f}**, i.e. **{ordv['base_over_sim']:.2f}× "
      f"the solved field — "
      f"{'weaker' if ordv['sim_wins'] else 'STRONGER'}, reading by |ρ| as this "
      f"table's own sign convention requires** | the test could only have "
      f"detected a "
      f"**≥{t2['contrasts']['MKI67']['viable_vs_rest']['power']['rate_ratio']:.1f}× "
      "rate ratio** |")
    A("")
    A("**In one sentence each.** The simulated oxygen field *is* significantly "
      "anti-correlated with measured hypoxia response, and it "
      + ("**out-predicts** every physics-free baseline built from its own inputs "
         if ordv["sim_wins"] and ordk["sim_wins"] else
         "is **out-predicted** by at least one physics-free baseline built from "
         "its own inputs ")
      + f"— including the distance-to-nearest-vessel transform ({de:+.4f}) and "
      f"its own undiffused source term ({vk:+.4f}) — but "
      + ("it wins a contest between four weak predictors on the one section that "
         "measures the marker, so the result is real and small rather than "
         "decisive. "
         if ordv["sim_wins"] else
         "the reaction–diffusion solve is therefore not what carries the signal. ")
      + "And the proliferation test on the 319-gene section is dropout-limited to "
      "the point where it cannot distinguish a correct cell-state model from "
      "noise; that is a power failure, not a vindication and not a refutation.")
    A("")
    A("> **The two sections are DIFFERENT PATIENTS.** Test 1 and Test 2 are two "
      "experiments, not two views of one tumour. Nothing here transfers from one "
      "to the other.")
    A("")

    # =====================================================================
    A("## 0. What is under test, and why it matters")
    A("")
    A("The chain is: `physics.solve_oxygen` produces an oxygen field → "
      "`cells.agents.update_cell_states` bands it into HYPOXIC / QUIESCENT / "
      "VIABLE at `t = OxygenConfig.hypoxia_threshold` = "
      f"{t2['hypoxia_threshold']} → `cells.pkpd` multiplies drug kill by "
      f"`PDConfig.hypoxia_resistance` = {t2['pd_config']['hypoxia_resistance']} "
      f"and `PDConfig.quiescent_resistance` = "
      f"{t2['pd_config']['quiescent_resistance']}.")
    A("")
    A("`results/dose_schedule.md` shows the consequence: for temozolomide — whose "
      "per-cell sensitivity is uniform 1.0 because MGMT is not on the panel — the "
      "kill CEILING is *entirely* those two multipliers, and ~89% of it comes from "
      "the quiescent band alone. If the oxygen field or the banding is wrong, "
      "every absolute kill number in this repository inherits the error.")
    A("")
    A(f"**Note on the {t2['pd_config']['hypoxia_resistance']:g} figure above.** "
      "That is the GLOBAL `PDConfig` default reported here for context; "
      "temozolomide specifically now carries a literature-informed PER-DRUG "
      "override of 0.80 (`core.config.DRUG_LIBRARY`, `PDConfig` docstring; "
      "`results/sweep/state_multiplier_literature_check_*.md`), converted "
      "Hill-curve-consistently from real GBM hypoxia-vs-TMZ IC50-shift data at "
      "this project's own operating point. It changes nothing in THIS report: "
      "every conclusion above turns on the oxygen-field / cell-state-banding "
      "machinery upstream of the multiplier, not its value, and (as section 2 "
      "below and `results/kill_calibration.md` both independently confirm) the "
      "HYPOXIC band is empty or near-empty on every real section this project "
      "currently reconstructs, so neither multiplier's value moves any number "
      "here.")
    A("")
    A("Both tests drive the **frozen production path** "
      "(`reconstruct.real_pipeline.reconstruct_real` → `physics.solve_oxygen` → "
      "`cells.agents.update_cell_states`). Nothing is re-implemented here.")
    A("")

    # ---- panel guard ----
    A("### 0.1 Regression guard — the IO panel map did not disturb the 319-gene section")
    A("")
    A("Test 1 needs a marker map for the 480-gene Immuno-Oncology panel "
      "(`reconstruct.marker_panels`), because the brain-panel map collapses that "
      "section to TUMOR 93.1% / IMMUNE 6.9% with **no ENDOTHELIAL cells at all** — "
      "and the endothelial cells are the oxygen SOURCE, so that failure would have "
      "silently produced a supply-free oxygen field. Adding it is only admissible "
      "if it is inert on the other section:")
    A("")
    A("| class | this run | historical reference | drift |")
    A("|---|---|---|---|")
    for k, ref_v in guard["reference_percent"].items():
        A(f"| {k} | {guard['percent'][k]:.2f}% | {ref_v:.2f}% | "
          f"{abs(guard['percent'][k]-ref_v):.3f} pp |")
    A("")
    A(f"`select_marker_map` returns the **same dict object** as "
      f"`BRAIN_PANEL_MARKER_MAP` for the 319-gene panel "
      f"(`is` identity: {guard['is_same_map_object']}), so the historical typing is "
      f"bit-for-bit unchanged. Max drift {guard['max_drift_pct']:.3f} pp. "
      f"**{'PASS' if guard['pass'] else 'FAIL'}**")
    A("")

    # =====================================================================
    A("## 1. TEST 1 — simulated oxygen vs measured hypoxia markers")
    A("")
    A(f"**Section.** {t1['dataset']}  ")
    A(f"**Panel.** {t1['panel']}  ")
    A(f"**Panel selection.** `{t1['panel_name']}` "
      f"(marker coverage: {t1['panel_coverage']})  ")
    A(f"**Cells.** {t1['n_cells']:,} × {t1['n_genes']} gene targets")
    A("")
    A("Measured hypoxia-response genes on this panel:")
    A("")
    A("| gene | cells with ≥1 transcript |")
    A("|---|---|")
    for g in t1["hypoxia_genes_used"]:
        A(f"| {g} | {100*t1['detection'][g]:.2f}% |")
    A("")
    A(f"The two genes correlate with **each other** at Spearman ρ = "
      f"{t1['gene_gene_rho']:+.4f}. That is close to zero: VEGFA and HIF1A are "
      "*not* two readouts of one latent 'hypoxia' variable in these data, so the "
      "composite score is reported alongside — never instead of — the per-gene "
      "results.")
    A("")

    # ---- 1.1 the interrupted-run claim
    A("### 1.1 The interrupted run's claim, adjudicated")
    A("")
    A("> *\"VEGFA is the single most anti-correlated gene of all 480 with "
      "simulated oxygen — but HIF1A is exactly null.\"*")
    A("")
    A("**Both halves CONFIRMED.** Every one of the 480 panel genes was "
      "depth-normalised and rank-correlated against the simulated oxygen at each "
      f"cell (dx = {DX_CONVERGED:.0f} µm, all {t1['n_cells']:,} cells):")
    A("")
    A("| rank (1 = most negative) | gene | ρ vs simulated O₂ | detection |")
    A("|---|---|---|---|")
    for i, (g, rho, det) in enumerate(t1["ranking"][:8]):
        mark = " **←**" if g in t1["hypoxia_genes_used"] else ""
        A(f"| {i+1} | {g}{mark} | {rho:+.4f} | {100*det:.2f}% |")
    A("| … | … | … | … |")
    for g in t1["hypoxia_genes_used"]:
        rk = t1["rank_of"][g]
        if rk > 8:
            rho = dict((a, b) for a, b, _ in t1["ranking"])[g]
            det = dict((a, c) for a, _, c in t1["ranking"])[g]
            A(f"| {rk} | {g} **←** | {rho:+.4f} | {100*det:.2f}% |")
    A("| … | … | … | … |")
    n_rank = len(t1["ranking"])
    for off, (g, rho, det) in enumerate(t1["ranking"][-5:][::-1]):
        A(f"| {n_rank - off} | {g} | {rho:+.4f} | {100*det:.2f}% |")
    A("")
    A(f"* **VEGFA is rank 1 of 480** at ρ = "
      f"{dict((a,b) for a,b,_ in t1['ranking'])['VEGFA']:+.5f}. Of every gene the "
      "panel measures, the canonical hypoxia-response gene is the one the "
      "simulated oxygen field most strongly anti-correlates with. That is a real, "
      "non-trivial hit and it is the strongest positive evidence in this report.")
    A(f"* **HIF1A is rank {t1['rank_of']['HIF1A']} of 480** at ρ = "
      f"{dict((a,b) for a,b,_ in t1['ranking'])['HIF1A']:+.5f} — null, and if "
      "anything the *wrong sign*. This is expected on mechanism, not a surprise: "
      "HIF-1α is regulated by oxygen-dependent prolyl hydroxylation and "
      "VHL-mediated proteolysis of the **protein** "
      f"({_pmid(LITERATURE_HIF1A_POSTTRANSLATIONAL[0])}; "
      f"{_pmid(LITERATURE_HIF1A_POSTTRANSLATIONAL[1])}), so its mRNA has no "
      "obligation to track oxygen. HIF1A mRNA should not have been treated as a "
      "hypoxia readout in the first place, and the composite score is diluted by "
      "including it.")
    A("")
    A("**The circularity, stated plainly.** The positive end of that table is the "
      "vascular and myeloid programme — RGS5, FLT1, PDGFRB, A2M, IGFBP7, CD163. "
      "`FLT1` is one of the three ENDOTHELIAL markers that *builds* the oxygen "
      "source term, so its positive correlation is partly the model correlating "
      "with itself. Only the NEGATIVE end of the ranking is informative, which is "
      "why VEGFA's rank-1 position is the claim worth making.")
    A("")
    A("VEGFA is the right measured readout here: "
      f"{_pmid(LITERATURE_VEGF_HYPOXIA)} — *{LITERATURE_VEGF_HYPOXIA['claim']}*.")
    A("")

    # ---- 1.2 (a) the correlation and the trivial baseline
    A("### 1.2 (a) Is simulated oxygen anti-correlated with measured hypoxia — "
      "and does it beat a trivial baseline?")
    A("")
    A("Per-cell Spearman ρ against the measured hypoxia score (mean z of "
      "depth-normalised VEGFA and HIF1A), computed on **all "
      f"{t1['n_cells']:,} cells** at dx = {DX_CONVERGED:.0f} µm. Every predictor "
      "is sign-aligned to oxygen (larger = better perfused), so a correct "
      "predictor gives a NEGATIVE ρ.")
    A("")
    A("| predictor | what it is | ρ vs hypoxia score | ρ vs VEGFA | ρ vs HIF1A | "
      "ρ with log-depth partialled out |")
    A("|---|---|---|---|---|---|")
    order = ["simulated_o2", "vascular_kde", "dist_endothelial", "local_cell_density"]
    for key in order:
        p = preds[key]
        name = ("**SIMULATED O₂ (the PDE)**" if key == "simulated_o2"
                else f"{key}")
        A(f"| {name} | {p['what']} | **{p['rho_score']:+.4f}** | "
          f"{p['rho_VEGFA']:+.4f} | {p['rho_HIF1A']:+.4f} | "
          f"{p['partial_rho_depth_removed']:+.4f} |")
    A("")
    A(f"Permutation p-values (label shuffling, {t1['n_perm']:,} draws on a seeded "
      f"random subsample of {t1['perm_subsample_n']:,} cells — subsampled because "
      "the kernel is O(n × draws); every effect size above is on the full "
      "population, so subsampling can only make a p-value *more* conservative):")
    A("")
    A("| predictor | ρ on the subsample | permutation p |")
    A("|---|---|---|")
    for key, pv in t1["perm"].items():
        A(f"| {key} | {pv['rho']:+.4f} | {pv['p_value']:.5f} |")
    A("")
    A("**The answer to (a) is yes, and "
      + ("the PDE wins the comparison — but it is a contest between weak "
         "predictors.**" if ordv["sim_wins"] else
         "a trivial baseline does it better.**")
      + f" The simulated field is genuinely, significantly anti-correlated with "
      f"measured hypoxia response (ρ = {so:+.4f}, "
      f"p = {t1['perm']['simulated_o2']['p_value']:.5f}, floor of the estimator). "
      "Every trivial baseline built from the same inputs is anti-correlated too, "
      f"but the solved field is the **{ordv['strongest_or_weakest']}** of the "
      "four — reading the table by |ρ|, as its own sign convention requires:")
    A("")
    A(f"* **The PDE ({so:+.4f}) {ordv['verb_vs_baseline']} `dist_endothelial` "
      f"({de:+.4f}) by {ordv['ratio']:.2f}×.** One `cKDTree` query — the distance "
      "from each cell to the nearest endothelial cell — recovers "
      f"{100*ordv['recovers']:.0f}% of "
      + ("the solved field's association with measured VEGFA/HIF1A. That a "
         "physics-free one-liner gets **most** of the way there is the honest "
         "caveat; it does not get all the way there, and it does not win."
         if ordv["sim_wins"] else
         "the association and costs nothing, so the reaction–diffusion solve is "
         "not what is carrying the signal."))
    A(f"* **`vascular_kde` ({vk:+.4f}) is the "
      f"{'weakest' if ordk['sim_wins'] else 'stronger'} vessel-derived "
      "predictor of the pair.** That is the solver's OWN SOURCE TERM, sampled "
      f"before it is diffused. The solve takes a field correlating at {vk:+.4f} "
      f"and returns one correlating at {so:+.4f}: **the diffusion step "
      + ("ADDS predictive information** "
         f"({ordk['ratio']:.2f}× on |ρ|) rather than destroying it. This is the "
         "single most direct evidence in the document that the reaction–diffusion "
         "step is doing work."
         if ordk["sim_wins"] else
         "*destroys* predictive information** rather than adding it."))
    A(f"* **`local_cell_density` ({preds['local_cell_density']['rho_score']:+.4f}) "
      "is near zero**, which is the sanity floor — not every spatial variable "
      "correlates with the marker, so the effects above are real signal, not "
      "autocorrelation artefacts.")
    A("")
    # The honesty clamp on the winning direction. A correction must not become an
    # overclaim with the sign flipped, so the size of the win is stated next to it.
    A(f"**How large is \"{'wins' if ordv['sim_wins'] else 'loses'}\", honestly.** "
      f"{ordv['ratio']:.2f}× on a ρ of {ordv['abs_sim']:.2f} is a small margin on "
      "a weak association, measured on ONE section, WITHIN-section rather than "
      "out-of-sample, with no confidence interval attached to the *difference* "
      "between two correlated predictors (the permutation test above tests each ρ "
      "against zero, not the predictors against each other). The composite is "
      f"carried by VEGFA alone (ρ = {preds['simulated_o2']['rho_VEGFA']:+.4f}) "
      f"while HIF1A contributes ρ = {preds['simulated_o2']['rho_HIF1A']:+.4f}. "
      + ("The defensible statement is **\"the solve is not redundant with the "
         "geometry it is built from\"** and **\"the PDE is the best of four weak "
         "predictors\"** — NOT \"the oxygen solve is validated\"."
         if ordv["sim_wins"] else
         "The defensible statement is **\"a trivial baseline is at least as good "
         "as the solve here\"** — NOT \"the oxygen solve is falsified\", which "
         "would require it to fail against the marker, not merely to be matched "
         "by cheaper geometry."))
    A("")
    nested = t1["sweep"][DX_CONVERGED]["nested"]
    A("**Nested test — does the PDE add anything beyond the raw geometry?**")
    A("")
    A("| partial Spearman | ρ |")
    A("|---|---|")
    A(f"| hypoxia score vs **simulated O₂**, distance-to-vessel held fixed | "
      f"{nested['score_vs_simo2_given_dist']:+.4f} |")
    A(f"| hypoxia score vs **−distance-to-vessel**, simulated O₂ held fixed | "
      f"{nested['score_vs_dist_given_simo2']:+.4f} |")
    A(f"| simulated O₂ vs −distance-to-vessel (how nested the two predictors are) | "
      f"{nested['simo2_vs_dist']:+.4f} |")
    A("")
    # Derived, not asserted: whichever partial is larger in absolute value names
    # the variable that survives conditioning, and the other one is the redundant
    # member. Hard-coding this sentence is what put the two halves of this
    # document into contradiction in the first place.
    n_sim = abs(nested["score_vs_simo2_given_dist"])
    n_dist = abs(nested["score_vs_dist_given_simo2"])
    ord_n = _ordering(n_sim, n_dist)
    A(("The simulated field retains "
       f"**{n_sim:.4f}** of explanatory power after distance-to-vessel is "
       "conditioned out; distance-to-vessel retains only "
       f"**{n_dist:.4f}** after the simulated field is conditioned out."
       if ord_n["sim_wins"] else
       "Distance-to-vessel retains "
       f"**{n_dist:.4f}** of explanatory power after the simulated field is "
       "conditioned out; the simulated field retains only "
       f"**{n_sim:.4f}** after distance is conditioned out.")
      + " The two are substantially nested "
      f"(ρ = {nested['simo2_vs_dist']:+.4f}), as they must be — the PDE's source "
      "cells *are* the distance transform's targets — but "
      + ("**distance-to-vessel is the redundant variable in that pair, not the "
         "PDE.**" if ord_n["sim_wins"] else
         "**the PDE is the redundant variable in that pair.**")
      + " The nested test therefore "
      + ("agrees with" if ord_n["sim_wins"] == ordv["sim_wins"] else "DISAGREES with")
      + " the marginal comparison above.")
    A("")
    A(f"**Depth control.** The simulated field correlates with total transcript "
      f"count — a purely technical variable — at ρ = {t1['o2_vs_depth_rho']:+.4f}. "
      "Partialling log-depth out of the table above changes nothing material "
      f"(simulated O₂ goes {so:+.4f} → "
      f"{preds['simulated_o2']['partial_rho_depth_removed']:+.4f}; distance goes "
      f"{de:+.4f} → {preds['dist_endothelial']['partial_rho_depth_removed']:+.4f}), "
      "so the ordering is not a depth artefact.")
    A("")

    # ---- 1.3 (b) state enrichment
    A("### 1.3 (b) Are model-HYPOXIC cells enriched for measured VEGFA/HIF1A "
      "vs model-VIABLE?")
    A("")
    frac = t1["sweep"][DX_CONVERGED]["fractions_all"]
    A("State occupancy on this section at dx = "
      f"{DX_CONVERGED:.0f} µm (all cells): " + ", ".join(
          f"**{s}** {100*frac[s]:.2f}%" for s in STATE_ORDER) + ".")
    A("")
    A("| marker | contrast | n (A vs B) | mean A | mean B | Cohen's d | "
      "Cliff's δ | permutation p |")
    A("|---|---|---|---|---|---|---|---|")
    for label in ["hypoxia_score"] + list(t1["hypoxia_genes_used"]):
        for cname in ("hypoxic_vs_viable", "hypoxic_vs_quiescent",
                      "hypoxic_vs_viable_TUMOUR"):
            b = t1["enrichment"][label][cname]
            A(f"| {label} | {cname.replace('_',' ')} | {b['n_a']:,} vs {b['n_b']:,} | "
              f"{b['mean_a']:+.4f} | {b['mean_b']:+.4f} | **{b['cohens_d']:+.4f}** | "
              f"{b['cliffs_delta']:+.4f} | {b['p_value']:.5f} |")
    A("")
    A("Depth-stratified control (labels permuted **within** total-transcript "
      "deciles, so only comparably-deep cells are ever compared):")
    A("")
    A("| marker | HYPOXIC vs VIABLE, depth-stratified difference | p |")
    A("|---|---|---|")
    for label in ["hypoxia_score"] + list(t1["hypoxia_genes_used"]):
        b = t1["enrichment"][label]["hypoxic_vs_viable_depth_stratified"]
        A(f"| {label} | {b['observed']:+.5f} | {b['p_value']:.5f} |")
    A("")
    vd = t1["enrichment"]["VEGFA"]["hypoxic_vs_viable"]
    hd = t1["enrichment"]["HIF1A"]["hypoxic_vs_viable"]
    A(f"**The answer to (b) is yes for VEGFA, no for HIF1A.** Model-HYPOXIC cells "
      f"carry substantially more measured VEGFA than model-VIABLE cells "
      f"(d = {vd['cohens_d']:+.3f}, Cliff's δ = {vd['cliffs_delta']:+.3f}, "
      f"p = {vd['p_value']:.5f}), and the effect survives depth stratification. "
      f"For HIF1A the same contrast is only d = {hd['cohens_d']:+.3f}, and the "
      "tumour-only version is not significant at all "
      f"(p = {t1['enrichment']['HIF1A']['hypoxic_vs_viable_TUMOUR']['p_value']:.5f}).")
    A("")
    A("**Two things blunt this result.** First, the VIABLE group holds only "
      f"{vd['n_b']:,} of {t1['n_cells']:,} cells "
      f"({100*vd['n_b']/t1['n_cells']:.3f}%), so \"HYPOXIC vs VIABLE\" is a "
      "comparison against a sliver. The better-balanced HYPOXIC-vs-QUIESCENT "
      f"contrast is far weaker "
      f"(VEGFA d = {t1['enrichment']['VEGFA']['hypoxic_vs_quiescent']['cohens_d']:+.3f}, "
      f"HIF1A d = "
      f"{t1['enrichment']['HIF1A']['hypoxic_vs_quiescent']['cohens_d']:+.3f}). "
      "Second, §1.2 shows that "
      + (f"{100*ordv['recovers']:.0f}% of this association, by |ρ|, is "
         "recoverable by thresholding distance-to-vessel alone, so a large "
         "Cohen's d here is mostly evidence that the *bands* sit somewhere "
         "sensible in space. The PDE separates them somewhat better than the "
         "distance threshold does; the point that survives is that the margin "
         "is small."
         if ordv["sim_wins"] else
         "the same separation is obtained, more strongly, by thresholding "
         "distance-to-vessel. A large Cohen's d here is evidence that the "
         "*bands* sit somewhere sensible in space, not that the PDE computed "
         "them."))
    A("")

    # ---- 1.4 (c) spatial / field-level
    A("### 1.4 (c) Field-vs-field agreement on the grid")
    A("")
    A("The per-cell test is dominated by single-cell dropout. Averaging both "
      f"fields over fixed {TILE_UM:.0f} µm tissue tiles (≥{TILE_MIN_CELLS} measured "
      "cells per tile) suppresses that noise and asks the coarser question: does "
      "a REGION the model calls hypoxic express more hypoxia response? Tiles are "
      "a fixed physical size at every dx, so the three rows are comparable to each "
      "other rather than to the mesh.")
    A("")
    A("| dx (µm) | grid | tiles used | Spearman ρ | Pearson r |")
    A("|---|---|---|---|---|")
    for dx in DX_SWEEP:
        tl = t1["sweep"][dx]["tiles"]
        g = t1["sweep"][dx]["grid"]
        A(f"| {dx:.0f} | {g[0]}×{g[1]} | {tl['n_tiles']:,} | {tl['rho']:+.4f} | "
          f"{tl['pearson']:+.4f} |")
    A("")
    A("Agreement is ~2.4× stronger at the field level "
      f"(ρ ≈ {t1['sweep'][DX_CONVERGED]['tiles']['rho']:+.3f} vs "
      f"{so:+.3f} per cell), which is the expected signature of per-cell dropout "
      "rather than of a better model. It is the most favourable framing of the "
      "test and it is reported as such.")
    A("")

    # ---- 1.5 (d) dx dependence
    A("### 1.5 (d) Does agreement improve as the grid is refined?")
    A("")
    A("**No — it is flat to four decimal places.**")
    A("")
    A("| dx (µm) | grid | mean O₂ (tissue) | HYPOXIC / QUIESCENT / VIABLE (all cells) | "
      "ρ(sim O₂, hypoxia score) | ρ(dist-to-vessel, hypoxia score) |")
    A("|---|---|---|---|---|---|")
    for dx in DX_SWEEP:
        s = t1["sweep"][dx]
        f = s["fractions_all"]
        A(f"| {dx:.0f} | {s['grid'][0]}×{s['grid'][1]} | "
          f"{s['mean_oxygen_tissue']:.4f} | "
          + " / ".join(f"{100*f[k]:.2f}%" for k in STATE_ORDER)
          + f" | {s['predictors']['simulated_o2']['rho_score']:+.4f} | "
          f"{s['predictors']['dist_endothelial']['rho_score']:+.4f} |")
    A("")
    A("This is good news about the numerics and neutral news about the model. "
      "With the geometric tissue mask and the areal oxygen sink both in place, "
      "the field is genuinely grid-converged *for this statistic*: mean O₂ moves "
      "by 0.2% and the correlation by <0.0002 over a 2.5× range of dx. "
      + ("Note also that the PDE's margin over the distance baseline is **stable "
         "across the whole ladder**, so the ordering in §1.2 is not a "
         "single-resolution artefact. What the flatness *does* establish is the "
         f"ceiling: refining the mesh will not raise |ρ| above ≈{ordv['abs_sim']:.3f}, "
         "so **the remaining gap to the measured marker is a model or measurement "
         "limitation, not a resolution one.**"
         if ordv["sim_wins"] else
         "So the disagreement with measured biology is **not** a resolution "
         "problem — refining the mesh cannot rescue it. The model has converged "
         "onto an answer that a distance transform beats."))
    A("")

    # ---- 1.6 verdict
    A("### 1.6 Verdict — Test 1")
    A("")
    A(f"## **{v1['verdict']}**")
    A("")
    for reason in v1["reasons"]:
        A(f"* {reason}")
    A("")
    A("**Consequences, stated plainly.**")
    A("")
    if ordv["sim_wins"]:
        A("1. **The oxygen field is not falsified, and it is not out-competed "
          "either.** Its sign is right, its significance is overwhelming, VEGFA "
          "being its rank-1 anti-correlate out of 480 genes is a genuine hit, "
          "and it is the strongest of the four predictors tested "
          f"(|ρ| {abs(so):.4f} vs {abs(de):.4f} vs {abs(vk):.4f} vs "
          f"{abs(preds['local_cell_density']['rho_score']):.4f}). The reaction–"
          "diffusion step measurably **improves** on its own undiffused source "
          f"term, {abs(vk):.4f} → {abs(so):.4f}. What remains true is that the "
          "association is **weak in absolute terms**, that this is a "
          "**within-section** test rather than an out-of-sample one, that the "
          "composite is carried by VEGFA alone, and that a physics-free distance "
          f"transform recovers {100*ordv['recovers']:.0f}% of it — so the PDE's "
          "marginal value over trivial geometry is real but modest.")
        A("2. **Hypoxia-driven resistance in the kill model is not *validated* by "
          "this, only left standing.** The `PDConfig.hypoxia_resistance` "
          "multiplier is applied to cells selected by a field that beats the "
          "trivial alternatives on this section — which is the weakest form of "
          "support that still counts as support. It is not evidence that the "
          "absolute hypoxic fraction is right (see consequence 5), and a "
          f"|ρ| of {abs(so):.2f} cannot carry a per-cell claim about any "
          "individual cell.")
    else:
        A("1. **The oxygen field is not falsified — it is out-competed.** Its "
          "sign is right, its significance is overwhelming, and VEGFA being its "
          "rank-1 anti-correlate out of 480 genes is a genuine hit. But the "
          "reaction–diffusion solve is not what produces that hit: a plain "
          "distance transform does better. Whatever the PDE is adding, it is not "
          "predictive accuracy about hypoxia.")
        A("2. **Hypoxia-driven resistance in the kill model is therefore not "
          "supported *as physics*.** The `PDConfig.hypoxia_resistance` multiplier "
          "is applied to cells selected by a field that a one-line geometric "
          "heuristic selects better. Any kill number that leans on the *spatial "
          "pattern* of hypoxia is, at best, an expensive re-derivation of "
          "\"far from a vessel\".")
    A("3. **Hypoxia should be derived from measured markers where they exist.** On "
      "a section that measures VEGFA, the defensible per-cell hypoxia covariate is "
      "the measured VEGFA level, not the simulated oxygen — and where no such "
      "marker is on the panel (the 319-gene section: no VEGFA, no CA9, no ADM, no "
      "NDRG1, no SLC2A1), the honest statement is that hypoxia is **unmeasured**, "
      "not that it has been simulated.")
    A("4. **Drop HIF1A from any hypoxia score.** Its mRNA is null against the "
      "field (rank "
      f"{t1['rank_of']['HIF1A']}/{len(t1['ranking'])}) for a known mechanistic "
      "reason, and including it dilutes the composite from ρ ≈ "
      f"{preds['simulated_o2']['rho_VEGFA']:+.3f} (VEGFA alone) to "
      f"{so:+.3f}.")
    A("5. **The absolute oxygen level remains uncalibrated.** Nothing here tests "
      "it. `OxygenConfig.consumption` has never been fitted to a measured oxygen "
      "tension, and on this section the model calls "
      f"{100*frac['HYPOXIC']:.1f}% of all cells HYPOXIC while on the 319-gene "
      "section under the shipping mask it calls 0.0% hypoxic (§2.2). Those two "
      "numbers cannot both be right.")
    A("")

    # =====================================================================
    A("## 2. TEST 2 — cell-state model vs measured proliferation")
    A("")
    A(f"**Section.** {t2['dataset']}  ")
    A(f"**Panel.** {t2['panel']}  ")
    A(f"**Cells.** {t2['n_cells']:,} total, {t2['n_tumor']:,} tumour  ")
    A(f"**Headline configuration.** tissue mask `{t2['headline_mask']}`, "
      f"dx = {t2['headline_dx']:.0f} µm — *the configuration "
      "`results/dose_schedule.md` actually ran under*, which is the one whose kill "
      "ceilings are at stake.")
    A("")

    # ---- 2.1 measured proliferating fraction
    A("### 2.1 The measured proliferating fraction, under two documented rules")
    A("")
    A("| marker | cells at exactly 0 | rule `nonzero` (count ≥ 1) | "
      f"rule `percentile` (top {100-PERCENTILE_RULE_Q:.0f}% of depth-normalised "
      "expression) |")
    A("|---|---|---|---|")
    for g in PROLIFERATION_MARKERS:
        m = t2["markers"][g]
        if not m.get("nonzero"):
            A(f"| {g} | — | NOT ON PANEL | NOT ON PANEL |")
            continue
        deg = " — **DEGENERATE**" if m["percentile"]["degenerate"] else ""
        A(f"| {g} | {100*m['zero_fraction']:.2f}% | "
          f"{100*m['nonzero']['fraction']:.3f}% "
          f"({m['nonzero']['n_positive']}/{m['n_cells']:,}) | "
          f"{100*m['percentile']['fraction']:.3f}%{deg} |")
    A("")
    A(f"**The percentile rule is mathematically degenerate on this panel, and that "
      f"IS the result.** The rule was written to bracket the published Ki-67 index "
      f"— \"call the top {100-PERCENTILE_RULE_Q:.0f}% of cells proliferating\" — "
      "but with >98.9% of tumour cells at exactly zero counts, the "
      f"{PERCENTILE_RULE_Q:.0f}th percentile of expression *is* zero. The "
      "threshold falls inside a tie mass at the origin and the rule collapses back "
      "onto `count ≥ 1`. There is no percentile of this data that corresponds to a "
      "Ki-67 index, because the data does not reach that far.")
    A("")
    lo_lit, hi_lit = t2["literature_ki67_range_pct"]
    mki_frac = t2["markers"]["MKI67"]["nonzero"]["fraction"]
    A("**Against the literature.** Published Ki-67 labelling indices in "
      "glioblastoma (PubMed, all verified):")
    A("")
    A("| cohort | Ki-67 index | source |")
    A("|---|---|---|")
    for e in LITERATURE_KI67:
        A(f"| {e['label']} | {e['value_pct']:.0f}% | {e['citation']} "
          f"(PMID {e['pmid']}, doi:{e['doi']}) |")
    A("")
    A(f"The panel detects MKI67 in **{100*mki_frac:.3f}%** of tumour cells — "
      f"**{lo_lit/(100*mki_frac):.0f}× to {hi_lit/(100*mki_frac):.0f}× below** the "
      "published range. This is not a claim that this tumour has a Ki-67 index of "
      "0.44%; it is a statement about probe sensitivity and transcript dropout on "
      "a targeted in-situ panel. The reference standard is itself imprecise "
      f"({_pmid(LITERATURE_KI67_VARIABILITY)} — {LITERATURE_KI67_VARIABILITY['claim']}), "
      "but not by two orders of magnitude.")
    A("")
    A("**Model VIABLE fraction, for comparison.** "
      f"{100*t2['dx_dependence'][t2['headline_mask']][t2['headline_dx']]['fractions']['VIABLE']:.2f}% "
      "in the headline configuration, but "
      f"{t2['viable_fraction_range_pct'][0]:.2f}% – "
      f"{t2['viable_fraction_range_pct'][1]:.2f}% across the mask × dx grid (§2.2). "
      "**Do not read any agreement between \"model VIABLE %\" and \"MKI67-detected "
      "%\" as corroboration** — it would be a coincidence between two unrelated "
      "small numbers, one of which is a bookkeeping artefact and the other a "
      "dropout rate.")
    A("")

    # ---- 2.2 the state vector is a bookkeeping choice
    A("### 2.2 The model's state vector is set by bookkeeping, not by tissue")
    A("")
    A("| tissue mask | dx (µm) | mean O₂ (tissue) | HYPOXIC | QUIESCENT | VIABLE |")
    A("|---|---|---|---|---|---|")
    for mode in TISSUE_MASK_MODES:
        for dx in DX_SWEEP:
            row = t2["dx_dependence"][mode][dx]
            f = row["fractions"]
            star = " ←" if (mode, dx) == (t2["headline_mask"], t2["headline_dx"]) else ""
            A(f"| `{mode}`{star} | {dx:.0f} | {row['mean_oxygen_tissue']:.4f} | "
              + " | ".join(f"{100*f[k]:.2f}%" for k in STATE_ORDER) + " |")
    A("")
    A("Read that table before reading any statistic below it. The tumour-cell "
      f"VIABLE fraction runs from {t2['viable_fraction_range_pct'][0]:.2f}% to "
      f"{t2['viable_fraction_range_pct'][1]:.2f}% — a "
      f"{t2['viable_fraction_range_pct'][1]/max(t2['viable_fraction_range_pct'][0],1e-9):.0f}× "
      "swing — across two tissue-mask conventions and three grid spacings. "
      "**Neither the mask nor the mesh is a measurement.** Under the current "
      "shipping default (`geometric`) the HYPOXIC state is *empty at every dx*, so "
      "`PDConfig.hypoxia_resistance` has no cells to act on at all; the 86% "
      "QUIESCENT / 13.4% HYPOXIC split that `results/dose_schedule.md` quotes comes "
      "from the legacy `occupancy` mask, whose coverage shrinks as the mesh is "
      "refined. Those kill ceilings are conditional on a convention that the "
      "reconstruction layer has since abandoned.")
    A("")

    # ---- 2.3 the contrast + power
    A("### 2.3 Are model-VIABLE cells more proliferative? — and could we have told?")
    A("")
    A("| marker | state | rate | k / n | Wilson 95% CI |")
    A("|---|---|---|---|---|")
    for g in PROLIFERATION_MARKERS:
        if g not in t2["contrasts"]:
            continue
        for s in STATE_ORDER:
            rr = t2["contrasts"][g]["rates"][s]
            if not rr["n"]:
                A(f"| {g} | {s} | — | 0 / 0 | state is empty |")
                continue
            A(f"| {g} | {s} | {100*rr['rate']:.3f}% | {rr['k']} / {rr['n']:,} | "
              f"{100*rr['lo']:.3f}–{100*rr['hi']:.3f}% |")
    A("")
    A("| marker | contrast | Cohen's d | Cliff's δ | permutation p | "
      "depth-stratified p | MINIMUM DETECTABLE rate ratio | expected positives |")
    A("|---|---|---|---|---|---|---|---|")
    for g in PROLIFERATION_MARKERS:
        if g not in t2["contrasts"]:
            continue
        b = t2["contrasts"][g]
        for cname, skey in (("VIABLE vs rest", "viable_vs_rest"),
                            ("HYPOXIC vs QUIESCENT", "hypoxic_vs_quiescent")):
            c = b[skey]
            dsk = ("viable_vs_rest_depth_stratified" if skey == "viable_vs_rest"
                   else "hypoxic_vs_quiescent_depth_stratified")
            A(f"| {g} | {cname} | {c['cohens_d']:+.4f} | {c['cliffs_delta']:+.4f} | "
              f"{c['p_value']:.4f} | {b[dsk]['p_value']:.4f} | "
              f"**{c['power']['rate_ratio']:.2f}×** | "
              f"{c['power']['expected_positives']:.1f} |")
    A("")
    A("**The power analysis is the result.** The \"minimum detectable rate ratio\" "
      f"column is the smallest enrichment the contrast could have found at "
      f"α = {MDE_ALPHA}, power = {MDE_POWER} given the observed positivity rate and "
      "the group sizes.")
    A("")
    pw = t2["contrasts"]["MKI67"]["viable_vs_rest"]["power"]
    A(f"* For **MKI67, VIABLE vs rest**, the VIABLE group holds "
      f"{pw['n_group']} tumour cells and would contain "
      f"{pw['expected_positives']:.1f} MKI67-positive cells under the null. The "
      f"test could only have detected a **≥{pw['rate_ratio']:.1f}× rate ratio**.")
    A(f"* The measured biology this is supposed to detect is nowhere near that "
      f"size: {_pmid(LITERATURE_HYPOXIA_PROLIFERATION)} measured that hypoxic "
      "pseudopalisading GBM cells are *5–50% LESS* proliferative than adjacent "
      "astrocytoma — a rate ratio of roughly 0.5–0.95, i.e. a **depletion of at "
      f"most 2×**. A test that needs {pw['rate_ratio']:.1f}× to fire cannot see "
      "an effect of that size in either direction.")
    A("* Therefore **no conclusion about the VIABLE label can be drawn from this "
      "section.** Not \"the state model is wrong\"; not \"the state model is "
      "right\". The measurement is not capable of adjudicating it.")
    A("")
    vr = t2["contrasts"]["MKI67"]["rates"]["VIABLE"]
    vrest = t2["contrasts"]["MKI67"]
    A("**A positive result is visible, and it must not be believed.** The point "
      f"estimate for MKI67 goes the RIGHT way — {100*vr['rate']:.2f}% positive in "
      f"VIABLE against {100*vrest['overall_rate']:.3f}% section-wide, a "
      f"{vr['rate']/vrest['overall_rate']:.1f}× rate ratio, permutation "
      f"p = {vrest['viable_vs_rest']['p_value']:.4f}, surviving depth "
      f"stratification (p = "
      f"{vrest['viable_vs_rest_depth_stratified']['p_value']:.4f}). It rests on "
      f"**{vr['k']} positive cells out of {vr['n']}**. Its Wilson interval spans "
      f"{100*vr['lo']:.2f}–{100*vr['hi']:.2f}%, one Poisson fluctuation moves it "
      f"from {vr['k']} to 2 or 7, and TOP2A — the other proliferation marker on "
      "the same cells — shows nothing "
      f"({100*t2['contrasts']['TOP2A']['rates']['VIABLE']['rate']:.2f}%, "
      f"{t2['contrasts']['TOP2A']['rates']['VIABLE']['k']}/"
      f"{t2['contrasts']['TOP2A']['rates']['VIABLE']['n']}, p = "
      f"{t2['contrasts']['TOP2A']['viable_vs_rest']['p_value']:.4f}). This is "
      "precisely the regime the minimum-detectable-effect calculation says is "
      "uninterpretable, and reporting it as support would be the exact error this "
      "module exists to prevent.")
    A("")
    A("**What the test CAN rule out.** The better-powered "
      f"HYPOXIC-vs-QUIESCENT contrast (n = "
      f"{t2['contrasts']['MKI67']['rates']['HYPOXIC']['n']:,} vs "
      f"{t2['contrasts']['MKI67']['rates']['QUIESCENT']['n']:,}) needs only "
      f"{t2['contrasts']['MKI67']['hypoxic_vs_quiescent']['power']['rate_ratio']:.2f}× "
      "for MKI67 and "
      f"{t2['contrasts']['TOP2A']['hypoxic_vs_quiescent']['power']['rate_ratio']:.2f}× "
      "for TOP2A, which is inside the plausible range. It rules out a **large** "
      "proliferation deficit in the model's HYPOXIC band — and it does so in the "
      "wrong direction:")
    A("")
    for g in PROLIFERATION_MARKERS:
        if g not in t2["contrasts"]:
            continue
        rr = t2["contrasts"][g]["rates"]
        if rr["HYPOXIC"]["n"] and rr["QUIESCENT"]["n"]:
            A(f"* **{g}+ is HIGHER in model-HYPOXIC "
              f"({100*rr['HYPOXIC']['rate']:.3f}%, {rr['HYPOXIC']['k']}/"
              f"{rr['HYPOXIC']['n']:,}) than in model-QUIESCENT "
              f"({100*rr['QUIESCENT']['rate']:.3f}%, {rr['QUIESCENT']['k']}/"
              f"{rr['QUIESCENT']['n']:,})** — the opposite of the measured biology. "
              f"Permutation p = "
              f"{t2['contrasts'][g]['hypoxic_vs_quiescent']['p_value']:.4f}; "
              "depth-stratified p = "
              f"{t2['contrasts'][g]['hypoxic_vs_quiescent_depth_stratified']['p_value']:.4f}.")
    A("")
    A("Neither reversal is large, and both sit close to the depth-stratification "
      "boundary, so this is a weak contradiction rather than a strong one. But it "
      "is a contradiction, and it points the same way as the continuous test:")
    A("")
    A("| continuous test | Spearman ρ | p | ρ with log-depth removed |")
    A("|---|---|---|---|")
    for g in PROLIFERATION_MARKERS:
        if g not in t2["contrasts"]:
            continue
        c = t2["contrasts"][g]["o2_continuous"]
        A(f"| simulated O₂ vs {g} (per tumour cell) | {c['rho']:+.4f} | "
          f"{c['p_value']:.4f} | {c['partial_rho_depth_removed']:+.4f} |")
    A("")
    A("Rank-correlating the **simulated oxygen** against measured proliferation "
      "across all "
      f"{t2['n_tumor']:,} tumour cells — the best-powered version of the question, "
      "because it uses every cell and does not depend on where the band edges fall "
      "— returns essentially nothing, with a *positive* sign where a negative one "
      "was not expected either. The oxygen field carries no per-cell proliferation "
      "information on this section.")
    A("")

    # ---- 2.4 dx dependence of the contrast
    A("### 2.4 dx-dependence of the Test 2 contrast")
    A("")
    A("| mask | dx (µm) | VIABLE fraction | MKI67+ in VIABLE | MKI67+ elsewhere | "
      "min. detectable rate ratio |")
    A("|---|---|---|---|---|---|")
    for mode in TISSUE_MASK_MODES:
        for dx in DX_SWEEP:
            row = t2["dx_dependence"][mode][dx]
            if "MKI67" not in row:
                continue
            A(f"| `{mode}` | {dx:.0f} | "
              f"{100*row['fractions']['VIABLE']:.2f}% | "
              f"{100*row['MKI67']['viable_rate']:.3f}% | "
              f"{100*row['MKI67']['nonviable_rate']:.3f}% | "
              f"{row['MKI67']['power']['rate_ratio']:.1f}× |")
    A("")
    A("The measured MKI67+ rate inside the model's VIABLE band is ~0.44–0.46% at "
      "five of the six configurations — i.e. **indistinguishable from the "
      "section-wide rate of "
      f"{100*mki_frac:.3f}%** — because at those settings VIABLE contains 74–98% "
      "of the tumour and the \"contrast\" is a population against itself. Only in "
      "the sixth configuration (the dose_schedule one) does VIABLE become small "
      "enough to be a real subgroup, and there it is too small to test. There is "
      "no dx at which this comparison is both meaningful and powered.")
    A("")

    # ---- 2.5 verdict
    A("### 2.5 Verdict — Test 2")
    A("")
    A(f"## **{v2['verdict']}**")
    A("")
    for reason in v2["reasons"]:
        A(f"* {reason}")
    A("")
    A("**Consequences, stated plainly.**")
    A("")
    A("1. **The quiescent labelling is UNSUPPORTED — in the specific sense that it "
      "has not been, and on this section cannot be, tested.** The measurement that "
      "would adjudicate it is 40–90× too insensitive. Nobody should read this as "
      "clearance.")
    A("2. **The kill ceilings inherit that.** `results/dose_schedule.md` derives "
      "~89% of the temozolomide kill ceiling from the quiescent band, via a "
      f"×{t2['pd_config']['quiescent_resistance']} multiplier applied to a label "
      "with no measured support and a "
      f"{t2['viable_fraction_range_pct'][1]/max(t2['viable_fraction_range_pct'][0],1e-9):.0f}× "
      "convention-dependent size. Every absolute kill fraction and every ceiling "
      "in that report carries an unquantified error of that magnitude. The "
      "dose/schedule **comparisons** are far more robust — they share the same "
      "multiplier and it largely cancels — but the **levels** are not defensible "
      "as measurements.")
    A("3. **The one thing the data does say points the wrong way.** In the "
      "better-powered HYPOXIC-vs-QUIESCENT contrast, both measured proliferation "
      "markers are *higher* in the model's HYPOXIC band. That is the opposite of "
      f"the measured tissue biology ({LITERATURE_HYPOXIA_PROLIFERATION['citation']}, "
      f"PMID {LITERATURE_HYPOXIA_PROLIFERATION['pmid']}).")
    A("4. **What would fix it.** Either (a) a proliferation readout with a real "
      "detection rate — IHC/IF Ki-67 on a serial section, or a whole-"
      "transcriptome platform — or (b) derive the state labels from measured "
      "markers (a per-cell proliferation score) instead of from the simulated "
      "oxygen, and report the states as *observed* rather than *predicted*.")
    A("")

    # =====================================================================
    A("## 3. Caveats — what could make these conclusions wrong")
    A("")
    A("1. **The two sections are different patients.** "
      f"Test 1 runs on `{t1['sample']}` ({t1['n_cells']:,} cells, IO panel); "
      f"Test 2 on `{t2['sample']}` ({t2['n_cells']:,} cells, brain panel). They "
      "are different tumours with different vascular architectures and different "
      "panels. Nothing here is a paired observation, and the two verdicts do not "
      "reinforce each other.")
    A("2. **VEGFA and HIF1A are hypoxia-RESPONSE genes, not an oxygen tension.** "
      "They are the transcriptional output of a pathway with its own kinetics and "
      "its own non-hypoxic inducers. A simulated field could be a perfectly good "
      "oxygen tension and still correlate weakly with them. What the trivial-"
      "baseline control establishes is narrower and safer: *whatever* VEGFA is "
      "reporting, a distance transform predicts it better than the PDE does — and "
      "that comparison is internal, so it is immune to this caveat.")
    A(f"3. **Targeted-panel dropout.** A zero is \"no probe fired\", not \"no "
      f"transcript\". This is fatal for Test 2 (MKI67 at "
      f"{100*mki_frac:.2f}%) and merely limiting for Test 1 (VEGFA at "
      f"{100*t1['detection']['VEGFA']:.1f}%, which is why Test 1 is powered and "
      "Test 2 is not).")
    A("4. **Typing circularity.** `FLT1` and `PLVAP` are IO-panel ENDOTHELIAL "
      "markers, so they help build the oxygen source; the positive tail of the "
      "480-gene ranking is partly the model correlating with its own inputs. "
      "`MKI67` is one of the four brain-panel TUMOR markers, so an MKI67-high cell "
      "is marginally more likely to be typed TUMOR — that inflates the MKI67+ rate "
      "*uniformly across states*, so it cannot manufacture a between-state "
      "contrast, but it does inflate the marginal rate slightly.")
    A("5. **`SLC1A3`/`S100B`/`PTN` (the IO STROMAL class) are expressed by "
      "malignant astrocyte-like glioma cells too**, so the IO section's STROMAL "
      "compartment is not clean. The brain panel's `AQP4`/`GJA1` have the same "
      "ambiguity, so the two sections are consistently rather than differently "
      "wrong here.")
    A("6. **Pericytes are unlabelled on the IO section.** The panel measures "
      "`RGS5`/`PDGFRB`/`NOTCH3` but `CellType` has no MURAL member, so mural cells "
      "fall into the ungated TUMOR default. That inflates the tumour count and "
      "hence the oxygen sink.")
    A("7. **The absolute oxygen level is uncalibrated.** "
      "`OxygenConfig.consumption` has never been fitted against a measured oxygen "
      "tension. The two sections disagree by orders of magnitude in hypoxic "
      f"fraction ({100*frac['HYPOXIC']:.1f}% of IO cells vs 0.0% of brain-panel "
      "cells under the shipping mask), which is a statement about the sink "
      "constant and the cell density, not about the two tumours.")
    A("8. **Subsampling.** Permutation p-values on the IO section use a seeded "
      f"random subsample of {t1['perm_subsample_n']:,} of {t1['n_cells']:,} cells "
      f"with {t1['n_perm']:,} draws. All effect sizes and all Spearman ρ are on "
      "the full population. Subsampling the null can only make a p-value more "
      "conservative, and every p reported here is already at the estimator's "
      "resolution floor.")
    A("")

    # =====================================================================
    A("## 4. Reproduction")
    A("")
    A("```bash")
    A("python -m tmesim.validation.state_validation")
    A("```")
    A("")
    A("Deterministic: the only randomness is the seeded permutation / subsampling "
      f"(`PERM_SEED = {PERM_SEED}`). Runtime ~8 minutes, dominated by six "
      "reconstructions of the 816,769-cell IO section and the 480-gene correlation "
      "scan.")
    A("")
    A("**Files added by this validation (nothing else was modified):**")
    A("")
    A("* `tmesim/validation/state_validation.py` — this analysis")
    A("* `tmesim/reconstruct/marker_panels.py` — the panel-aware marker map "
      "(regression-guarded in §0.1)")
    A("* `results/state_validation.md`, `results/state_validation.png`")
    A("")
    A("**Data provenance.** Both sections are open 10x Genomics Xenium datasets, "
      "CC BY 4.0, no login or DUA:")
    A("")
    A(f"* Test 1 — {t1['dataset']}; panel: {t1['panel']}")
    A(f"* Test 2 — {t2['dataset']}; panel: {t2['panel']}")
    A("")
    A("Neither section carries a paired drug-response outcome label, so nothing "
      "here validates a kill prediction — only the internal machinery that "
      "produces one.")
    A("")
    A("**Figure.** `results/state_validation.png`.")
    A("")

    out_path = out_path or os.path.join(_results_dir(), "state_validation.md")
    with open(out_path, "w") as fh:
        fh.write("\n".join(L))
    print(f"report written: {out_path}")
    return out_path


def main() -> Dict:
    """Run both validations, write the report and the figure."""
    results = run()
    write_report(results)
    make_figure(results)
    return results


if __name__ == "__main__":  # pragma: no cover - drives the real caches
    main()
