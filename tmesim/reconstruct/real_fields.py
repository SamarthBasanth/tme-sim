"""REAL tissue physical fields (vasculature, ECM, oxygen) from REAL cells.

This module replaces the *synthetic* field generators (``reconstruct.vasculature``
grows invented random-walk vessel trees; the oxygen self-test there fabricates a
radial tumor). Here every field is derived from the measured positions and types
of real segmented cells in a reconstructed :class:`~tmesim.core.Microenvironment`.

Design contract
---------------
These builders operate on the CORE interface only — a
:class:`~tmesim.core.Domain` plus a :class:`~tmesim.core.CellPopulation` (parallel
``x``/``y``/``cell_type`` arrays). They do NOT read per-cell marker expression:
cell *typing* (which cell is ENDOTHELIAL / STROMAL / TUMOR) is the job of the
reconstruction step (``reconstruct/real_xenium.py``) that produces the population.
That keeps this file independent of how the cells were classified — it consumes a
finished population regardless of who built it.

What is real vs. proxy (honest accounting)
-------------------------------------------
* ``vascular_density_from_cells`` — REAL. Vessel locations are the measured
  positions of endothelial (PECAM1+ / FLT1+) cells. Since
  ``results/vessel_promotion.md`` the DEFAULT construction is ``mode='discrete'``:
  those cells are clustered into discrete vessel objects and painted as a source
  that is exactly zero away from a vessel. The historical 40 um-bandwidth
  Gaussian KDE is retained as ``mode='kde'``. Vessels are real cell locations in
  either mode; no vessel geometry is invented.
* ``ecm_density_from_cells`` — PROXY. ECM is not itself a cell, so we cannot bin
  it directly. We use the density of STROMAL / reactive-mesenchymal cells (the
  cells that DEPOSIT matrix: CAF-like fibroblasts, reactive astrocytes) as a
  spatial proxy for ECM deposition. This is a proxy, clearly labeled as such.

  A HIGHER-FIDELITY ECM field is possible on THIS panel: the Xenium Human Brain
  panel used here DOES carry real matrix genes — collagens ``COL1A2``,
  ``COL12A1``, ``COL25A1``, plus ``POSTN`` (periostin), ``VCAN`` (versican),
  ``DCN`` (decorin), ``FBLN1`` (fibulin-1), ``LAMA2`` (laminin) and the
  matricellular ``THBS1``/``SPON1``/``TGFBI``. So true ECM could be a
  KDE of ECM-transcript abundance rather than of stromal cell counts. We do NOT
  do that here because the core ``CellPopulation`` interface carries no
  expression vector — only positions and types. When marker expression is
  available upstream (as it is in ``real_xenium.py``), prefer an expression-
  weighted ECM field; ``ecm_density_from_positions`` accepts optional per-cell
  weights precisely so a caller can pass ECM-marker abundance through.

The oxygen field is then produced by the EXISTING, unmodified physics solver
(``physics.oxygen.solve_oxygen``): steady-state reaction-diffusion with the real
vascular field as the O2 source and the real tumor-cell density as the O2 sink.
So oxygen is real-perfusion-driven, not painted on.
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy.ndimage import gaussian_filter

from ..core import (
    CellPopulation,
    CellType,
    Domain,
    Microenvironment,
    OxygenConfig,
    TissueFields,
)
from ..physics.oxygen import solve_oxygen

# Default KDE bandwidth (Gaussian sigma) in micrometers. Endothelium marks vessel
# walls a few cell-diameters apart; ~40 um smooths discrete cell hits into a
# continuous perfusion / deposition field without erasing real heterogeneity.
_DEFAULT_BANDWIDTH_UM = 40.0


def _kde_on_grid(
    domain: Domain,
    x: np.ndarray,
    y: np.ndarray,
    bandwidth_um: float = _DEFAULT_BANDWIDTH_UM,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Gaussian KDE of point positions onto the domain grid, normalized to [0,1].

    Cells are binned to their grid index (``Domain.to_index``, accumulated with
    ``np.add.at``), the count field is Gaussian-smoothed with sigma
    ``bandwidth_um / dx`` cells, masked to tissue, then divided by its own max so
    the densest region reads 1.0. Fully deterministic (no RNG).

    Parameters
    ----------
    x, y : (n,) float — cell positions in micrometers.
    bandwidth_um : Gaussian bandwidth in micrometers.
    weights : optional (n,) float — per-cell weight (e.g. marker abundance). If
        None every selected cell contributes 1.0.
    """
    density = domain.zeros()
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size == 0:
        return density

    col, row = domain.to_index(x, y)
    if weights is None:
        np.add.at(density, (row, col), 1.0)
    else:
        w = np.asarray(weights, dtype=np.float64)
        np.add.at(density, (row, col), w)

    sigma_cells = max(bandwidth_um / domain.dx, 1e-6)
    density = gaussian_filter(density, sigma=sigma_cells, mode="nearest")

    # Respect tissue extent (no signal outside simulated tissue).
    density *= domain.tissue_mask.astype(np.float64)

    dmax = float(density.max())
    if dmax > 0.0:
        density /= dmax
    return np.clip(density, 0.0, 1.0)


#: The vascular-source construction to use.
#:
#: ``'discrete'`` is the DEFAULT as of ``results/vessel_promotion.md``: the
#: segmented-vessel field of :mod:`reconstruct.vessel_network`, exactly zero away
#: from a vessel. ``'kde'`` is the historical 40 um-bandwidth endothelial KDE and
#: is retained, unchanged, as an explicit option.
#:
#: The promotion settled the two blockers ``vessel_network_impact.md`` §8.2 left
#: open, and both failed to survive measurement:
#:
#: * **alpha does not need re-fitting.** ``validation.calibration.calibrate_alpha``
#:   is a closed-form solve on the compiled drug plus the published anchor; no
#:   vascular field enters it. It returns 148.608... under either mode, and
#:   ``validation.pk_multidrug.run()`` is byte-identical. Forcing an *empirical*
#:   section-scale re-fit instead BREAKS the out-of-sample check under BOTH
#:   fields, so 148.6 stands.
#: * **grid convergence holds.** Coverage @1x IC50, mean tumour AUC and kill are
#:   converged for the discrete source by dx = 12 um (dx 8 -> 6 changes 0.0010
#:   abs / 0.18% / 0.0002 abs), and the source's own area fraction is ~4x MORE
#:   dx-stable than the max-normalised KDE's.
#:
#: TWO REAL COSTS, both quantified in ``results/vessel_promotion.md``:
#:
#: 1. **Distance-resolved claims need dx <= 8 um.** The fitted penetration length
#:    of the shortest-L drug runs 47.2 / 39.6 / 51.5 / 51.9 um at dx =
#:    15 / 12 / 8 / 6, so the production dx = 12 um value sits ~24% below the
#:    dx <= 8 asymptote. Bulk metrics are fine at dx = 12; quote ``L`` at
#:    dx <= 8. (Under the KDE, ``L`` was already converged at dx = 12.)
#: 2. **Tighter advection CFL** -- but the numbers below are MEASURED ON A
#:    SUPERSEDED RECONSTRUCTION and overstate the constraint on the shipped one.
#:    A concentrated source does steepen the IFP gradient, so the max
#:    interstitial speed is ~1.6-2x the KDE's; that part holds.
#:
#:    The figures originally recorded here -- dt_max 20.5 / 15.4 / 9.5 / 6.9 s at
#:    dx = 15 / 12 / 8 / 6, and "dt = 10 s VIOLATES the CFL at dx <= 8" -- were
#:    taken before ``4ad221f`` changed ``TYPING_AGGREGATION`` from ``mean_z`` to
#:    ``matched_null_z``, which moved this section from 3,984 endothelial /
#:    30,174 tumour cells to 4,295 / 23,812. On the SHIPPED field, re-measured:
#:
#:        dx (um)        15      12       8       7       6
#:        max speed    0.7077  0.6336  0.7513  0.6684  0.6773
#:        dt_max (s)   21.20   18.94   10.65   10.47    8.86
#:        CFL @ dt=10   0.472   0.528   0.939   0.955   1.129 (solver refuses)
#:
#:    So dt = 10 s is admissible down to dx = 7, with little margin at dx <= 8,
#:    and dx = 6 is refused. Prefer dt <= 5 s below dx = 8 for headroom rather
#:    than for stability. See ``results/cfl_resolution.md``; the solver gates on
#:    this itself and raises rather than degrading silently
#:    (``physics/diffusion.py``), so a rung that produced a number was stable.
#:
#: NORMALISATION WARNING. The two fields integrate to different totals over
#: tissue (mean 0.0250 vs 0.0840, i.e. 3.36x), because a max-normalised KDE
#: spread over ~100% of the tissue integrates to more than a unit source on 2.5%
#: of it. Neither amplitude is a validated absolute, but the referents differ:
#: KDE 1.0 = "the densest endothelial-KDE spot in THIS section" (section-
#: dependent); discrete 1.0 = "a perfused vessel wall leaking at the reference
#: rate" (section-independent). Flipping the default therefore changes the UNIT
#: of every absolute exposure by ~3.4x. Any absolute number in ``results/``
#: generated before this change is stale; see the list in
#: ``results/vessel_promotion.md`` §7.
VASCULAR_MODES = ("kde", "discrete")
_DEFAULT_VASCULAR_MODE = "discrete"


def vascular_density_from_cells(
    domain: Domain,
    cells: CellPopulation,
    bandwidth_um: float = _DEFAULT_BANDWIDTH_UM,
    mode: str = _DEFAULT_VASCULAR_MODE,
    **discrete_kwargs,
) -> np.ndarray:
    """Perfused-vessel-fraction field (ny, nx) in [0,1] from REAL endothelium.

    ``mode='discrete'`` (**DEFAULT** since ``results/vessel_promotion.md``) — the
    measured endothelial cells clustered into DISCRETE putative vessel objects
    and painted as a source that is nonzero only within a vessel radius and
    EXACTLY zero elsewhere
    (:func:`reconstruct.vessel_network.segment_vessels` +
    :func:`~reconstruct.vessel_network.discrete_vascular_field`). Extra keyword
    arguments (``link_distance_um``, ``radius_um``, ``amplitude``,
    ``match_total_to``) are forwarded. ``bandwidth_um`` is ignored in this mode.
    See :data:`VASCULAR_MODES` for the evidence, the two costs (dx <= 8 um for
    distance-resolved claims; dt <= 5 s for the advection CFL) and the 3.36x
    normalisation change.

    ``mode='kde'`` (the historical default, byte-for-byte unchanged) — KDE of the
    measured positions of
    ``CellType.ENDOTHELIAL`` cells (the PECAM1+/FLT1+ vessel-lining cells
    identified upstream). The value at each grid cell is the relative local
    density of endothelium — a real perfused-vessel fraction, high exactly where
    vessel cells actually cluster. No vessel geometry is invented. This field is
    SMOOTH: at the default 40 um bandwidth it is nonzero over ~99.7% of the
    tissue, which is the structural confound documented in
    ``results/SESSION_AUDIT.md`` §3, and the reason the default moved away from
    it. Passing ``mode='kde'`` reproduces the historical field exactly.

    Returns an all-zero field if the population contains no endothelial cells.
    """
    m = str(mode).strip().lower()
    if m == "kde":
        if discrete_kwargs:
            raise TypeError(
                f"unexpected keyword(s) {sorted(discrete_kwargs)} for mode='kde'")
        endo = cells.cell_type == int(CellType.ENDOTHELIAL)
        return _kde_on_grid(domain, cells.x[endo], cells.y[endo], bandwidth_um)
    if m == "discrete":
        # Imported lazily so the KDE path takes no new dependency and there is no
        # import cycle (vessel_network imports only from ..core).
        from .vessel_network import discrete_vascular_field, segment_vessels

        seg_kw = {}
        if "link_distance_um" in discrete_kwargs:
            seg_kw["link_distance_um"] = float(discrete_kwargs.pop("link_distance_um"))
        vessels = segment_vessels(cells, domain, **seg_kw)
        return discrete_vascular_field(domain, vessels, **discrete_kwargs)
    raise ValueError(f"unknown vascular mode {mode!r}; expected one of {VASCULAR_MODES}")


def ecm_density_from_cells(
    domain: Domain,
    cells: CellPopulation,
    bandwidth_um: float = _DEFAULT_BANDWIDTH_UM,
) -> np.ndarray:
    """ECM-density PROXY field (ny, nx) in [0,1] from REAL stromal cells.

    KDE of the measured positions of ``CellType.STROMAL`` cells (reactive /
    mesenchymal, CAF-like — the matrix-DEPOSITING population) used as a spatial
    proxy for extracellular-matrix density. Higher stromal density => more matrix
    deposition => higher effective ECM (raising tortuosity / lowering drug
    diffusivity downstream).

    PROXY caveat: ECM is acellular, so it cannot be counted directly from a cell
    population. On the Xenium Human Brain panel used in this project, REAL ECM
    transcripts ARE available (collagens COL1A2/COL12A1/COL25A1, POSTN, VCAN, DCN,
    FBLN1, LAMA2). A truer ECM field would KDE those transcript abundances; the
    core ``CellPopulation`` interface carries no expression, so we approximate
    with stromal-cell density here. Use :func:`ecm_density_from_positions` with
    ECM-marker weights when expression is available upstream.

    Returns an all-zero field if the population contains no stromal cells.
    """
    stro = cells.cell_type == int(CellType.STROMAL)
    return _kde_on_grid(domain, cells.x[stro], cells.y[stro], bandwidth_um)


def ecm_density_from_positions(
    domain: Domain,
    x: np.ndarray,
    y: np.ndarray,
    weights: Optional[np.ndarray] = None,
    bandwidth_um: float = _DEFAULT_BANDWIDTH_UM,
) -> np.ndarray:
    """Escape hatch for a higher-fidelity, expression-weighted ECM field.

    Same KDE as :func:`ecm_density_from_cells` but on arbitrary positions with
    optional per-point ``weights``. Pass the positions of ECM-marker-positive
    cells weighted by their summed collagen/POSTN/VCAN abundance to build ECM
    from REAL matrix transcripts rather than from a stromal-cell proxy.
    """
    return _kde_on_grid(domain, x, y, bandwidth_um, weights=weights)


def build_real_fields(
    domain: Domain,
    cells: CellPopulation,
    bandwidth_um: float = _DEFAULT_BANDWIDTH_UM,
    vascular_mode: str = _DEFAULT_VASCULAR_MODE,
) -> TissueFields:
    """Assemble the real ECM + vascular fields; leave oxygen/ifp unset.

    Oxygen is intentionally left ``None`` so the caller can run the physics solver
    on the assembled vascular field via :func:`add_oxygen` (keeping physics in the
    physics package). IFP is left for ``physics.pressure`` similarly.

    ``vascular_mode`` is forwarded to :func:`vascular_density_from_cells` and
    defaults to ``'discrete'`` (see :data:`VASCULAR_MODES`); pass ``'kde'`` for
    the historical smooth field.
    """
    ecm = ecm_density_from_cells(domain, cells, bandwidth_um)
    vascular = vascular_density_from_cells(domain, cells, bandwidth_um,
                                           mode=vascular_mode)
    return TissueFields(ecm_density=ecm, vascular_density=vascular,
                        oxygen=None, ifp=None)


def add_oxygen(micro: Microenvironment, cfg: OxygenConfig = OxygenConfig()) -> None:
    """Solve the real-perfusion-driven oxygen field and store it in place.

    Calls the EXISTING, unmodified ``physics.oxygen.solve_oxygen`` on ``micro``
    (which reads ``micro.fields.vascular_density`` as the O2 source and the real
    tumor-cell density as the O2 sink) and writes the result into
    ``micro.fields.oxygen``. Returns nothing.
    """
    oxygen = solve_oxygen(micro, cfg)
    micro.fields.oxygen = oxygen


# ---------------------------------------------------------------------------
# Self-test on REAL Xenium GBM data.
# ---------------------------------------------------------------------------
def _classify_cells_by_markers(expr, genes, cell_xy):
    """COARSE marker rule to get ENDOTHELIAL / STROMAL / TUMOR cells for the test.

    This is NOT the real reconstruction (that is deconvolution in real_xenium.py);
    it is only a simple, transparent rule so the field builders can be exercised
    on real coordinates. Each cell is counts-normalized into a lineage score;
    a cell is called ENDOTHELIAL or STROMAL only when that lineage's score clears
    a per-lineage percentile GATE and is the strongest of the specialised
    lineages, otherwise it defaults to TUMOR (the malignant background). Gating
    avoids the over-assignment a bare 3-way arg-max produces on rare markers.
    """
    genes = list(genes)
    E = expr.tocsc()
    n = E.shape[0]

    def score(marker_list):
        idx = [genes.index(g) for g in marker_list if g in genes]
        if not idx:
            return np.zeros(n)
        return np.asarray(E[:, idx].sum(axis=1)).ravel()

    total = np.asarray(E.sum(axis=1)).ravel()
    total = np.where(total <= 0, 1.0, total)

    endo = score(["PECAM1", "FLT1"]) / total
    stromal = score(["PDGFRA", "DCN", "VCAN", "COL1A2", "POSTN", "FBLN1"]) / total

    # Gate each specialised lineage at its 85th percentile among cells expressing
    # it at all (so gates track real marker abundance, not the zero-inflation).
    def gate(scr, q=0.85):
        pos = scr[scr > 0]
        return np.quantile(pos, q) if pos.size else np.inf

    endo_gate, stromal_gate = gate(endo), gate(stromal)
    is_endo = endo >= endo_gate
    is_stromal = (stromal >= stromal_gate) & ~is_endo  # endothelium wins ties

    ctype = np.full(n, int(CellType.TUMOR), dtype=np.int8)
    ctype[is_stromal] = int(CellType.STROMAL)
    ctype[is_endo] = int(CellType.ENDOTHELIAL)
    return ctype


if __name__ == "__main__":
    import os
    import time

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from tmesim.data.real.singlecell_spatial import load

    # --- load REAL Xenium GBM cells (retry once on transient/network error) ----
    try:
        data = load(dense=False)
    except Exception as exc:  # pragma: no cover - network/IO robustness
        print(f"load failed ({exc!r}); retrying once ...")
        time.sleep(3)
        data = load(dense=False)

    xy = data["cell_xy_um"]
    genes = data["genes"]
    expr = data["expression"]
    print(f"REAL cells loaded: {xy.shape[0]:,}  extent "
          f"x[{xy[:,0].min():.0f},{xy[:,0].max():.0f}] "
          f"y[{xy[:,1].min():.0f},{xy[:,1].max():.0f}] um")

    # --- report real ECM markers actually present on the panel -----------------
    ecm_panel = [g for g in ("COL1A2", "COL12A1", "COL25A1", "POSTN", "VCAN",
                             "DCN", "FBLN1", "LAMA2", "THBS1", "SPON1", "TGFBI")
                 if g in set(genes)]
    print(f"REAL ECM markers on panel ({len(ecm_panel)}): {ecm_panel}")

    # --- classify cells with a simple marker rule (test scaffolding only) ------
    ctype = _classify_cells_by_markers(expr, genes, xy)
    for t in (CellType.TUMOR, CellType.STROMAL, CellType.ENDOTHELIAL):
        print(f"  {t.name:12s}: {(ctype == int(t)).sum():6d} cells")

    # --- build a Domain covering the real extent -------------------------------
    dx = 20.0
    nx = int(np.ceil(xy[:, 0].max() / dx)) + 1
    ny = int(np.ceil(xy[:, 1].max() / dx)) + 1
    domain = Domain(nx=nx, ny=ny, dx=dx)

    rng = np.random.default_rng(0)
    cells = CellPopulation(
        x=xy[:, 0], y=xy[:, 1], cell_type=ctype,
        state=np.zeros(xy.shape[0], dtype=np.int8),
        sensitivity=rng.uniform(0, 1, xy.shape[0]),
    )

    # --- build real fields + real-perfusion oxygen -----------------------------
    fields = build_real_fields(domain, cells)
    micro = Microenvironment(domain=domain, fields=fields, cells=cells,
                             patient_id="xenium-gbm-real")
    assert fields.oxygen is None and fields.ifp is None, "build should leave O2/ifp unset"
    add_oxygen(micro, OxygenConfig())
    assert micro.fields.oxygen is not None, "add_oxygen must fill oxygen"

    vasc = micro.fields.vascular_density
    ecm = micro.fields.ecm_density
    oxy = micro.fields.oxygen
    for name, f in (("vascular", vasc), ("ecm", ecm), ("oxygen", oxy)):
        assert f.shape == domain.shape, f"{name} shape mismatch"
        assert f.min() >= 0.0 and f.max() <= 1.0, f"{name} out of [0,1]"

    # --- ASSERTION 1: vascular density HIGH where endothelial cells cluster ----
    endo = cells.cell_type == int(CellType.ENDOTHELIAL)
    ecol, erow = domain.to_index(cells.x[endo], cells.y[endo])
    vasc_at_endo = float(vasc[erow, ecol].mean())
    # compare against a matched sample of non-endothelial cell locations
    other = ~endo
    ocol, orow = domain.to_index(cells.x[other], cells.y[other])
    vasc_at_other = float(vasc[orow, ocol].mean())
    print(f"\nvascular @ endothelial cells = {vasc_at_endo:.4f}")
    print(f"vascular @ other cells       = {vasc_at_other:.4f}")
    assert vasc_at_endo > vasc_at_other, (
        "vascular density is not higher where endothelial cells cluster"
    )

    # --- ASSERTION 2: oxygen LOWER in dense-tumor / low-vessel regions ---------
    # tumor-cell density on the grid
    tum = cells.cell_type == int(CellType.TUMOR)
    tden = domain.zeros()
    tcol, trow = domain.to_index(cells.x[tum], cells.y[tum])
    np.add.at(tden, (trow, tcol), 1.0)
    tden = gaussian_filter(tden, sigma=_DEFAULT_BANDWIDTH_UM / dx, mode="nearest")

    tissue = tden > 0  # only compare where there is actual tumor tissue
    hi_tum_lo_vasc = tissue & (tden > np.quantile(tden[tissue], 0.7)) & \
        (vasc < np.quantile(vasc[tissue], 0.3))
    hi_vasc = tissue & (vasc > np.quantile(vasc[tissue], 0.7))
    o2_dense_tumor = float(oxy[hi_tum_lo_vasc].mean())
    o2_well_perfused = float(oxy[hi_vasc].mean())
    print(f"oxygen in dense-tumor/low-vessel = {o2_dense_tumor:.4f}")
    print(f"oxygen in well-perfused regions  = {o2_well_perfused:.4f}")
    assert o2_dense_tumor < o2_well_perfused, (
        "oxygen is not lower in dense-tumor/low-vessel regions"
    )

    # --- save figure: vascular + oxygen panels ---------------------------------
    os.makedirs("results", exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.2), constrained_layout=True)
    ext = [0, domain.extent_um[0], domain.extent_um[1], 0]
    im0 = ax[0].imshow(vasc, cmap="inferno", extent=ext, aspect="equal")
    ax[0].scatter(cells.x[endo], cells.y[endo], s=1, c="cyan", alpha=0.3,
                  linewidths=0, label="endothelial")
    ax[0].set_title("REAL vascular density (endothelial KDE)")
    ax[0].legend(loc="upper right", markerscale=6)
    fig.colorbar(im0, ax=ax[0], shrink=0.8)

    im1 = ax[1].imshow(ecm, cmap="viridis", extent=ext, aspect="equal")
    ax[1].set_title("ECM density PROXY (stromal KDE)")
    fig.colorbar(im1, ax=ax[1], shrink=0.8)

    im2 = ax[2].imshow(oxy, cmap="RdYlBu", extent=ext, aspect="equal", vmin=0, vmax=1)
    ax[2].set_title("REAL-perfusion oxygen (physics solve)")
    fig.colorbar(im2, ax=ax[2], shrink=0.8)
    for a in ax:
        a.set_xlabel("x (um)")
        a.set_ylabel("y (um)")
    fig.suptitle("Real Xenium GBM: physical fields from real cells", fontsize=14)
    out = os.path.join("results", "real_fields.png")
    fig.savefig(out, dpi=110)
    print(f"\nsaved {out}")
    print("\nSELF-TEST PASSED (real Xenium GBM data)")
