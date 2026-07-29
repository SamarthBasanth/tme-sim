"""Steady-state oxygen reaction-diffusion field and hypoxia labeling.

Oxygen matters for two reasons in the Stage 0 kill-test:

1. Tumor cores are poorly vascularized (see ``reconstruct.vasculature``), so O2
   diffuses in from vessels but is consumed by dense tumor tissue. The balance
   produces the classic hypoxic core.
2. Hypoxic cells are drug/therapy resistant, which the pharmacodynamic module
   (``cells.pkpd``) uses to modulate kill. This module supplies the O2 field and
   the per-cell HYPOXIC label that resistance keys off of.

We solve the *steady-state* reaction-diffusion balance

    D * lap(O) + supply * vasc * (1 - O) - consumption * tumor_density * O = 0

on the (ny, nx) grid with no-flux (Neumann) boundaries. Rearranged into a linear
system ``A @ O = b`` it is symmetric, sparse and diagonally dominant, so a single
sparse direct solve gives the exact steady state (no iteration/CFL concerns).

    diffusion:  D * lap(O)
    source:     supply * vasc          (turned off as O -> 1 via the (1-O) term)
    sink:       consumption * rho_tum  (proportional to LOCAL tumor-cell density)

Units follow the project convention: length um, diffusivity um^2/s. The oxygen
field is dimensionless relative saturation in [0, 1] (1 = fully oxygenated,
0 = anoxia).

MMHG REFERENCE POINT (read this before comparing O to any measured pO2)
-------------------------------------------------------------------------
``O`` had no defined physical unit until now. ``O = 1`` is DEFINED to represent
:data:`REFERENCE_PO2_MMHG` (35 mmHg, normal -- non-tumour-bearing -- human brain
TISSUE pO2, Erecinska & Silver 2001, PMID 11718758), not arterial/plasma pO2
(~100 mmHg): the PDE's ``vasc`` field is a vascular-density source term feeding
a tissue compartment, not a blood compartment, so the tissue anchor is the
correct physical match, and it is the anchor at which the EXISTING
``cells.agents.update_cell_states`` bands (``t=0.2`` HYPOXIC, ``[0.2,0.4)``
QUIESCENT, ``>=0.4`` VIABLE) convert to a biologically plausible range
(~7 / [7,14) / >=14 mmHg) rather than the arterial anchor's absurd one (~20 /
[20,40) / >=40 mmHg, which would call normal brain tissue "hypoxic"). See
``results/sweep/oxygen_mmhg_reference_proposal_2026-07-29T05-54-00Z.md`` for
the full derivation and both anchors' threshold tables. This conversion is a
DECLARATION about what the existing [0,1] number has always implicitly meant --
it changes no solved value of ``O``, no threshold, and no downstream label.
:func:`to_mmhg` is the one place to apply it.

GRID INDEPENDENCE OF THE SINK (fixed; read this)
------------------------------------------------
``rho_tum`` must be a *density* (cells per unit AREA), not a raw per-grid-cell
COUNT. An earlier revision binned tumor cells and fed the raw count straight into
the sink; a count scales as ``dx^2``, so halving the grid spacing quartered the
consumption term and the hypoxic fraction moved by almost an order of magnitude
(measured on the real Xenium reconstruction: hypoxic tumor fraction 0.937 / 0.364
/ 0.134 at dx = 30 / 20 / 12 um). :func:`_tumor_density` now divides the binned
count by the grid-cell area and re-expresses it against a fixed reference tile
(:data:`TUMOR_REFERENCE_AREA_UM2`), so the sink is an areal density and no longer
carries a ``dx^2`` factor.

THE SECOND HALF -- the tissue mask was also grid-dependent (now FIXED elsewhere)
---------------------------------------------------------------------------------
Fixing the sink was necessary but not sufficient for a grid-invariant oxygen field
on the real reconstruction. ``reconstruct.real_xenium`` used to define
``Domain.tissue_mask`` as "grid cells that contain at least one segmented cell",
which is itself a function of ``dx``: that mask covers 50.6% / 32.7% / 13.9% of the
rectangle at dx = 30 / 20 / 12 um. ``real_fields._kde_on_grid`` multiplies the
vascular field by that mask, so the total O2 SUPPLY shrank as the grid was
refined. Measured mean O2 at the (fixed, dx-independent) real cell positions:

    dx (um)                    30      20      12     max/min
    legacy count sink        0.154   0.214   0.251     1.63x
    areal sink, legacy mask  0.522   0.426   0.251     2.08x
    areal sink + GEOMETRIC mask (both fixes, current default)
                             0.570   0.565   0.558     1.02x

i.e. the two grid errors used to partially CANCEL, which is why the legacy spread
looked smaller. ``reconstruct.real_xenium.tissue_mask_from_cells`` now defines the
mask geometrically by default (a grid cell is tissue if its centre lies within a
fixed 30 um of a measured cell, so the mask is a union of fixed physical discs;
coverage 67.2 / 67.4 / 67.7% over the same grids), and with both fixes in place the
oxygen field on the REAL section is grid-invariant to ~2%. Nothing in THIS module
changed for that fix -- the mask lives in the reconstruction layer.

HONEST CAVEAT ON THE ABSOLUTE LEVEL (now the dominant uncertainty)
-------------------------------------------------------------------
``OxygenConfig.consumption`` has never been calibrated against a measured oxygen
tension. With the supply no longer artificially starved by a shrinking mask, mean
O2 on the real GBM section is ~0.56 and the predicted HYPOXIC TUMOR FRACTION IS 0.0
at every dx -- which, for a tumor type partly defined by pseudopalisading necrosis,
is no more credible than the 0.94 the coarse grid used to report. Both the old
ladder (0.94 -> 0.07 over dx = 30 -> 8 um) and the new 0.0 are un-calibrated model
output. What these fixes establish is that the answer no longer depends on the
numerical grid -- NOT that its level is right.

A LITERATURE PASS SHARPENED THIS UNCERTAINTY RATHER THAN RESOLVING IT
------------------------------------------------------------------------
Two follow-up literature checks each produced a real, unit-matched estimate for
``consumption``, and they disagree by roughly an order of magnitude IN OPPOSITE
DIRECTIONS relative to the current default (0.02 /s):

- Watabe et al. 2017 (rat C6 glioma, 15O-PET tumoral metabolic rate of oxygen,
  PMID 28210996, https://doi.org/10.1186/s13550-017-0263-6): TMRO2 = 1.9-2.8 mL
  O2/100 mL tissue/min, a real glioma-TISSUE-level (not cell-line, not a
  different tumor type) measurement in directly comparable units. Converting via
  the :data:`REFERENCE_PO2_MMHG` anchor above implies ``consumption`` should be
  roughly **12-25x HIGHER** than 0.02.
- A different reaction-diffusion-model estimate (Diepart et al., via Milotti
  2017) using the same linear-in-O functional form this module uses implies
  ``consumption`` should be roughly **4-100x LOWER** than 0.02.

Neither estimate is obviously wrong -- both are real, sourced, unit-matched
numbers -- and they bracket the current default from opposite sides by up to
two orders of magnitude combined. This is a SHARPER, more specific form of
"uncalibrated" than before (a quantified disagreement range with real sources
on both ends), not a resolution of it, and no value has been picked. See
``results/sweep/oxygen_ocr_literature_second_pass_2026-07-29T05-54-00Z.md``.

THE SINK'S FUNCTIONAL FORM IS A SEPARATE, ALSO-UNRESOLVED AXIS
-------------------------------------------------------------------
This module's sink is first-order in O (``consumption * rho_tum * O``); most of
the cited oxygen-consumption literature (e.g. Grimes, Fletcher & Partridge 2014,
PMID 26064525) instead uses a Michaelis-Menten form, ``a_max * O/(O + Km)``,
which does not vanish as O -> 0. A standalone prototype (unmodified
``solve_oxygen``, Picard-iterated MM sink, real flagship reconstruction) found
this is NOT a minor correction: sweeping the literature-plausible half-saturation
range swings the predicted hypoxic tumor FRACTION from the current 0.0 to
0.69-1.0, the opposite extreme. Neither 0.0 nor 69-100% is obviously correct
relative to the real human-GBM pO2 literature. This axis and the
``consumption``-magnitude uncertainty above are COUPLED, not independently
fixable: recalibrating ``consumption`` alone, without addressing the functional-
form mismatch, is unlikely to reproduce a realistic pO2 profile. See
``results/sweep/oxygen_functional_form_prototype_2026-07-29T05-54-00Z.md``.
Neither the functional form nor ``consumption`` has been changed in this module;
picking either without the other would be false precision, not a fix.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from ..core import CellState, Microenvironment, OxygenConfig

#: Reference tissue area (um^2) against which the tumor-cell areal density is
#: expressed, making ``rho_tum`` the dimensionless number of tumor cells per
#: reference tile instead of a grid-dependent per-cell count.
#:
#: Value = 144 um^2 = (12 um)^2. Two reasons:
#:   1. 12 um is the grid-converged production spacing the flagship runs at
#:      (``results/numerical_rigor.md``), so re-expressing the sink this way is a
#:      pure re-parameterisation there -- it changes nothing at dx = 12 um and
#:      removes the spurious dx^2 scaling at every other resolution. That keeps
#:      this a defect fix rather than a silent re-tuning of the model.
#:   2. It is close to the cross-sectional footprint of a single glioma cell in a
#:      2-D section (~100 um^2 for an ~11 um diameter cell), so ``rho_tum ~ 1``
#:      reads as "roughly one tumor cell per cell-sized tile" -- a physically
#:      interpretable occupancy.
#:
#: HONEST CAVEAT: ``OxygenConfig.consumption`` was never calibrated against a
#: measured oxygen tension, so the product ``consumption * TUMOR_REFERENCE_AREA_UM2``
#: is a single un-calibrated lumped constant. The ABSOLUTE hypoxic fraction is very
#: sensitive to it (a +/-40% change swings the real-section hypoxic tumor fraction
#: between ~0.0 and ~0.6). What this fix establishes is that the answer no longer
#: depends on the numerical grid -- not that the absolute level is calibrated.
TUMOR_REFERENCE_AREA_UM2 = 144.0

#: mmHg partial pressure that ``O = 1`` is DEFINED to represent. Anchored to
#: NORMAL (non-tumour-bearing) human brain TISSUE pO2 -- not arterial/plasma
#: pO2 -- because ``solve_oxygen``'s source term represents flux from vessels
#: INTO tissue, so ``O`` is a tissue-compartment quantity throughout; see the
#: module docstring's "MMHG REFERENCE POINT" section for why the arterial
#: anchor (~100 mmHg) was rejected.
#:
#: Value = 35.0 mmHg. Source: Erecinska M, Silver IA. "Tissue oxygen tension
#: and brain sensitivity to hypoxia." Respir Physiol 2001;128(3):263-276.
#: PMID 11718758, DOI 10.1016/s0034-5687(01)00306-1. Corroborated by (a) the
#: clinical PbtO2 (Licox-probe) monitoring literature's normal range of
#: 20-40 mmHg and (b) this project's own Beppu et al. 2002 GBM readings
#: (peritumoral ~59.8 mmHg, arterial ~109.2 mmHg, both above this
#: normal-tissue anchor, the expected ordering) -- see
#: ``results/oxygen_calibration_landscape.md`` §1 and
#: ``results/sweep/oxygen_mmhg_reference_proposal_2026-07-29T05-54-00Z.md``.
#:
#: HONEST CAVEAT: like :data:`TUMOR_REFERENCE_AREA_UM2`, this is a SINGLE POINT
#: anchor for a quantity that varies severalfold across individuals, brain
#: regions and measurement technique (Nozue et al. 1997, PMID 9290690, found
#: 1.5-5.6 mmHg tumour-pO2 spread across labs on an identical protocol). It
#: fixes what O=1 MEANS; it does NOT calibrate ``OxygenConfig.consumption``,
#: which independently sets where the SOLVED O actually lands on this new
#: scale -- see that field's own uncalibrated-level caveat above.
REFERENCE_PO2_MMHG = 35.0


def to_mmhg(oxygen: np.ndarray) -> np.ndarray:
    """Convert a dimensionless oxygen field/array to mmHg via :data:`REFERENCE_PO2_MMHG`.

    Pure linear rescale (``oxygen * REFERENCE_PO2_MMHG``). Touches nothing in
    ``solve_oxygen`` or ``label_hypoxia`` -- no threshold, no solved value, and
    no test-observable behavior changes. Exists so a caller who wants to
    compare the model's O to a measured pO2 (e.g. the four glioma studies in
    ``results/oxygen_calibration_landscape.md``) has one canonical place to do
    it instead of each re-deriving the anchor.
    """
    return np.asarray(oxygen, dtype=np.float64) * REFERENCE_PO2_MMHG


def _tumor_density(micro: Microenvironment) -> np.ndarray:
    """Local tumor-cell AREAL density on the grid (ny, nx), grid-spacing independent.

    Uses ``Domain.to_index`` to map each tumor cell's (x, y) in micrometers to its
    (col, row) grid index, accumulates counts with ``np.add.at`` (vectorized), then
    converts the count to a density::

        rho_tum = (count per grid cell) / dx^2 * TUMOR_REFERENCE_AREA_UM2
                = tumor cells per TUMOR_REFERENCE_AREA_UM2 of tissue

    The division by ``dx^2`` is the fix: a raw count scales as ``dx^2`` and made the
    oxygen sink -- hence the hypoxic fraction and every downstream hypoxia-gated
    resistance term -- depend on the numerical grid. An areal density does not.
    Multiplying by the fixed reference area keeps ``rho_tum`` a dimensionless
    occupancy of order 1 so ``OxygenConfig.consumption`` retains its meaning and
    magnitude (see :data:`TUMOR_REFERENCE_AREA_UM2`).

    Note this is a *point* density estimate (no smoothing): the oxygen diffusion
    length ``sqrt(D / (consumption * rho_tum))`` is O(100 um) >> dx, so the solver
    itself averages the field over many grid cells.
    """
    domain = micro.domain
    cells = micro.cells
    density = domain.zeros()
    tum = cells.tumor_mask()
    if not np.any(tum):
        return density
    col, row = domain.to_index(cells.x[tum], cells.y[tum])
    np.add.at(density, (row, col), 1.0)
    # count -> cells per um^2 -> cells per reference tile (grid-independent).
    density *= TUMOR_REFERENCE_AREA_UM2 / (float(domain.dx) ** 2)
    return density


def _build_laplacian(ny: int, nx: int, dx: float) -> sp.csr_matrix:
    """5-point finite-difference Laplacian on (ny, nx) with no-flux boundaries.

    No-flux (Neumann) boundaries are imposed by ghost cells equal to the interior
    cell, i.e. a missing neighbor simply drops out of the stencil. The returned
    operator ``L`` approximates ``lap(O)`` such that ``(L @ O.ravel())`` has the
    same sign convention as ``(neighbors - k*center)/dx^2``.
    """
    n = ny * nx
    idx = np.arange(n).reshape(ny, nx)
    inv = 1.0 / (dx * dx)

    rows = []
    cols = []
    vals = []
    diag = np.zeros((ny, nx), dtype=np.float64)

    # Four cardinal neighbors; each existing neighbor contributes +inv off-diagonal
    # and -inv to the center. Missing (boundary) neighbors are dropped => no flux.
    for src, dst, dslice in (
        (idx[1:, :], idx[:-1, :], (slice(1, None), slice(None))),   # north
        (idx[:-1, :], idx[1:, :], (slice(0, -1), slice(None))),     # south
        (idx[:, 1:], idx[:, :-1], (slice(None), slice(1, None))),   # west
        (idx[:, :-1], idx[:, 1:], (slice(None), slice(0, -1))),     # east
    ):
        rows.append(src.ravel())
        cols.append(dst.ravel())
        vals.append(np.full(src.size, inv))
        diag[dslice] -= inv

    rows.append(idx.ravel())
    cols.append(idx.ravel())
    vals.append(diag.ravel())

    L = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
        shape=(n, n),
    ).tocsr()
    return L


def solve_oxygen(micro: Microenvironment, cfg: OxygenConfig) -> np.ndarray:
    """Solve the steady-state reaction-diffusion oxygen field.

    Solves ``D*lap(O) + supply*vasc*(1-O) - consumption*rho_tum*O = 0`` on the
    domain grid with no-flux boundaries, where ``vasc`` is
    ``micro.fields.vascular_density`` (the O2 source) and ``rho_tum`` is the local
    tumor-cell AREAL density from :func:`_tumor_density` (the O2 sink) -- cells per
    :data:`TUMOR_REFERENCE_AREA_UM2`, which is independent of the grid spacing.
    The linear system

        (-D*L + diag(supply*vasc + consumption*rho_tum)) @ O = supply*vasc

    is symmetric positive-definite and solved directly with a sparse solver.

    Parameters
    ----------
    micro : Microenvironment
        Provides the domain, ``fields.vascular_density`` and the cell population.
    cfg : OxygenConfig
        ``diffusivity`` (um^2/s), ``supply`` and ``consumption`` rate constants.

    Returns
    -------
    np.ndarray
        Oxygen saturation field, shape (ny, nx), clipped to [0, 1]. The caller is
        responsible for storing it into ``micro.fields.oxygen``.
    """
    domain = micro.domain
    ny, nx = domain.shape
    n = ny * nx

    vasc = np.asarray(micro.fields.vascular_density, dtype=np.float64)
    if vasc.shape != (ny, nx):
        raise ValueError(
            f"vascular_density shape {vasc.shape} != domain shape {(ny, nx)}"
        )
    vasc = np.clip(vasc, 0.0, 1.0)
    rho_tum = _tumor_density(micro)

    D = float(cfg.diffusivity)
    supply = float(cfg.supply)
    consumption = float(cfg.consumption)

    L = _build_laplacian(ny, nx, domain.dx)

    source = (supply * vasc).ravel()                          # RHS b
    reaction = (supply * vasc + consumption * rho_tum).ravel()  # O-proportional loss

    # A = -D * L + diag(reaction). (-L is positive-definite; reaction >= 0.)
    A = (-D * L) + sp.diags(reaction, format="csr")

    # If there is no source and no reaction anywhere the system is singular; guard.
    if not np.any(source) and not np.any(reaction):
        return domain.zeros()

    oxygen = spsolve(A.tocsc(), source)
    oxygen = np.asarray(oxygen, dtype=np.float64).reshape(ny, nx)
    np.clip(oxygen, 0.0, 1.0, out=oxygen)
    return oxygen


def label_hypoxia(
    micro: Microenvironment, oxygen: np.ndarray, cfg: OxygenConfig
) -> None:
    """Mark non-dead cells sitting in low-oxygen grid cells as HYPOXIC (in place).

    Each cell's local oxygen is sampled at its (col, row) grid index. Cells that
    are not DEAD and sit where ``oxygen < cfg.hypoxia_threshold`` have their state
    set to ``CellState.HYPOXIC``. Modifies ``micro.cells.state`` in place; returns
    nothing.

    Parameters
    ----------
    micro : Microenvironment
        Population whose ``cells.state`` is updated.
    oxygen : np.ndarray
        Oxygen field (ny, nx), e.g. from :func:`solve_oxygen`.
    cfg : OxygenConfig
        Uses ``hypoxia_threshold``.
    """
    domain = micro.domain
    cells = micro.cells
    oxygen = np.asarray(oxygen, dtype=np.float64)
    if oxygen.shape != domain.shape:
        raise ValueError(
            f"oxygen shape {oxygen.shape} != domain shape {domain.shape}"
        )

    col, row = domain.to_index(cells.x, cells.y)
    local_o2 = oxygen[row, col]

    alive = cells.state != int(CellState.DEAD)
    hypoxic = alive & (local_o2 < float(cfg.hypoxia_threshold))
    cells.state[hypoxic] = np.int8(CellState.HYPOXIC)


if __name__ == "__main__":
    from ..core import (
        CellPopulation,
        CellType,
        Domain,
        TissueFields,
    )

    # ------------------------------------------------------------------
    # Build a synthetic Domain: vasculature is sparse in a central tumor
    # region (well-perfused rim, avascular core) with tumor cells packed
    # into that central region. Expect a hypoxic core.
    # ------------------------------------------------------------------
    rng = np.random.default_rng(0)
    ny = nx = 80
    dx = 15.0
    domain = Domain(nx=nx, ny=ny, dx=dx)

    X, Y = domain.coords_um()
    cx, cy = X.mean(), Y.mean()
    r = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    core_radius = 0.30 * min(domain.extent_um)

    # Background vasculature ~0.6, suppressed toward the core (avascular center).
    vasc = 0.6 - 0.59 * np.exp(-(r ** 2) / (2 * (0.5 * core_radius) ** 2))
    vasc = np.clip(vasc, 0.005, 1.0)

    ecm = domain.ones() * 0.4
    fields = TissueFields(ecm_density=ecm, vascular_density=vasc)

    # Tumor cells packed densely into the central core (strong O2 sink); a few
    # scattered stromal cells elsewhere. The count is set so the core reaches an
    # areal density of ~8 tumor cells per TUMOR_REFERENCE_AREA_UM2, which is what
    # it takes to drive a genuinely hypoxic core at the (un-calibrated) default
    # OxygenConfig.consumption -- see the note on TUMOR_REFERENCE_AREA_UM2.
    n_tumor = 24000
    ang = rng.uniform(0, 2 * np.pi, n_tumor)
    rad = core_radius * np.sqrt(rng.uniform(0, 1, n_tumor))
    tx = np.clip(cx + rad * np.cos(ang), 0, domain.extent_um[0] - 1e-3)
    ty = np.clip(cy + rad * np.sin(ang), 0, domain.extent_um[1] - 1e-3)

    n_stroma = 300
    sx = rng.uniform(0, domain.extent_um[0], n_stroma)
    sy = rng.uniform(0, domain.extent_um[1], n_stroma)

    x = np.concatenate([tx, sx])
    y = np.concatenate([ty, sy])
    ctype = np.concatenate([
        np.full(n_tumor, int(CellType.TUMOR)),
        np.full(n_stroma, int(CellType.STROMAL)),
    ])
    state = np.zeros(x.size, dtype=np.int8)
    sensitivity = rng.uniform(0, 1, x.size)
    cells = CellPopulation(x=x, y=y, cell_type=ctype, state=state,
                           sensitivity=sensitivity)

    micro = Microenvironment(domain=domain, fields=fields, cells=cells,
                             patient_id="oxygen-selftest")

    cfg = OxygenConfig()
    oxygen = solve_oxygen(micro, cfg)

    assert oxygen.shape == domain.shape
    assert oxygen.min() >= 0.0 and oxygen.max() <= 1.0, "oxygen out of [0,1]"

    # Determinism: solving again gives an identical field.
    oxygen2 = solve_oxygen(micro, cfg)
    assert np.array_equal(oxygen, oxygen2), "solve_oxygen is not deterministic"

    # Core vs periphery: mean O2 in the central core disc vs an outer annulus.
    core_mask = r < 0.5 * core_radius
    periph_mask = r > 1.3 * core_radius
    core_o2 = float(oxygen[core_mask].mean())
    periph_o2 = float(oxygen[periph_mask].mean())

    print(f"grid = {domain.shape}, dx = {dx} um")
    print(f"oxygen range = [{oxygen.min():.4f}, {oxygen.max():.4f}]")
    print(f"tumor-core mean O2  = {core_o2:.4f}")
    print(f"periphery mean O2   = {periph_o2:.4f}")
    assert core_o2 < periph_o2, (
        f"tumor core ({core_o2:.4f}) is not more hypoxic than periphery "
        f"({periph_o2:.4f})"
    )

    # Hypoxia labeling marks some (core) cells.
    label_hypoxia(micro, oxygen, cfg)
    n_hypoxic = int((micro.cells.state == int(CellState.HYPOXIC)).sum())
    print(f"hypoxia_threshold   = {cfg.hypoxia_threshold}")
    print(f"n hypoxic cells     = {n_hypoxic} / {micro.cells.n}")
    assert n_hypoxic > 0, "label_hypoxia marked no cells"

    # ------------------------------------------------------------------
    # GRID INDEPENDENCE of the sink: the SAME physical tissue solved on
    # grids spanning a 4x range of dx must give the same oxygen field.
    # (With the legacy raw-COUNT sink the hypoxic fraction ran 1.00 -> 0.00
    # over this same sweep.)
    # ------------------------------------------------------------------
    print("\ngrid-independence sweep (same physical tissue, same cells):")
    means = []
    for dx_test in (30.0, 20.0, 15.0, 7.5):
        n_test = int(round(domain.extent_um[0] / dx_test))
        dom_t = Domain(nx=n_test, ny=n_test, dx=dx_test)
        Xt, Yt = dom_t.coords_um()
        rt = np.sqrt((Xt - cx) ** 2 + (Yt - cy) ** 2)
        vt = np.clip(0.6 - 0.59 * np.exp(-(rt ** 2) / (2 * (0.5 * core_radius) ** 2)),
                     0.005, 1.0)
        micro_t = Microenvironment(
            domain=dom_t,
            fields=TissueFields(ecm_density=dom_t.ones() * 0.4, vascular_density=vt),
            cells=cells, patient_id="oxygen-dx-sweep",
        )
        o_t = solve_oxygen(micro_t, cfg)
        ct, rt_i = dom_t.to_index(cells.x, cells.y)
        local = o_t[rt_i, ct]
        means.append(float(local.mean()))
        print(f"  dx = {dx_test:5.1f} um  grid {dom_t.shape}  "
              f"mean O2 @ cells = {local.mean():.4f}  "
              f"hypoxic frac = {float((local < cfg.hypoxia_threshold).mean()):.4f}")
    spread = max(means) - min(means)
    print(f"  mean-O2 spread across a 4x dx range = {spread:.4f}")
    assert spread < 0.01, f"oxygen field is still grid-dependent (spread {spread:.4f})"

    print("\nself-test PASSED")
