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
on both ends). See
``results/sweep/oxygen_ocr_literature_second_pass_2026-07-29T05-54-00Z.md``.

**Resolution adopted (see ``results/oxygen_sensitivity_band.md``):** rather than
pick a single point in this range -- which would be false precision, not a fix
-- ``consumption`` remains at its current default (0.02, UNCHANGED) and the full
literature-defensible band this section derives (roughly 0.0002-0.5, combining
both estimates' brackets) has been swept end-to-end on the real flagship
reconstruction, for both sink forms, and reported as a sensitivity table. That
turns "we don't know" into "here is the quantified shape of the disagreement" --
a real resolution of the uncertainty's SHAPE, even though no single point in the
band is established as correct without real human GBM pO2 or OCR data this
project does not have.

THE SINK'S FUNCTIONAL FORM IS A SEPARATE, RELATED AXIS -- NOW AN AVAILABLE OPTION
-----------------------------------------------------------------------------------
This module's DEFAULT sink is still first-order in O
(``consumption * rho_tum * O``, ``OxygenConfig.sink_form = "linear"``); most of
the cited oxygen-consumption literature (e.g. Grimes, Fletcher & Partridge 2014,
PMID 26064525) instead uses a Michaelis-Menten form,
``consumption * rho_tum * O / (O + mm_half_saturation)``, which does not vanish
as O -> 0. A standalone prototype (Picard-iterated MM sink, real flagship
reconstruction) found this is NOT a minor correction: sweeping the literature-
plausible half-saturation range swings the predicted hypoxic tumor FRACTION from
0.0 to 0.69-1.0, the opposite extreme, holding ``consumption`` fixed. See
``results/sweep/oxygen_functional_form_prototype_2026-07-29T05-54-00Z.md`` for
that prototype and ``results/oxygen_sensitivity_band.md`` for the production
version's full grid over both this axis and ``consumption`` together.

**What changed in this module:** the Michaelis-Menten form is now an available,
documented, OPT-IN option -- ``OxygenConfig(sink_form="michaelis_menten")`` --
solved with a proper Newton-Raphson iteration (see :func:`_solve_michaelis_menten`),
not the prototype's throwaway Picard script. ``sink_form`` defaults to
``"linear"``, so every existing call site, test and published ``results/`` number
is UNCHANGED (:func:`solve_oxygen`'s linear branch is bit-identical to the
pre-existing code).

**What did NOT change:** which functional form, and which ``consumption`` value,
is CORRECT for this tissue. Neither 0.0 nor 69-100% hypoxic fraction is obviously
right relative to the real human-GBM pO2 literature (``results/
oxygen_calibration_landscape.md`` Sec.1), and the two axes are coupled (Sec.5 of
the prototype report). Resolving that needs real human GBM pO2 or tissue-level
OCR data this project does not have -- seeing the disagreement quantified across
a defensible grid (``results/oxygen_sensitivity_band.md``) is the honest
resolution available without inventing a number; picking a single point on that
grid and presenting it as calibrated would be false precision.
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

#: Valid values for ``OxygenConfig.sink_form``. "linear" is the ORIGINAL,
#: default sink (``consumption * rho_tum * O``); "michaelis_menten" is the
#: opt-in alternative (``consumption * rho_tum * O / (O + mm_half_saturation)``,
#: Grimes, Fletcher & Partridge 2014, PMID 26064525). See the module docstring
#: section "THE SINK'S FUNCTIONAL FORM..." and :func:`solve_oxygen`.
SINK_FORMS = ("linear", "michaelis_menten")


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


def _solve_michaelis_menten(
    L: sp.csr_matrix,
    D: float,
    source: np.ndarray,
    reaction_bulk: np.ndarray,
    K: float,
    max_iter: int,
    tol: float,
) -> tuple[np.ndarray, int, float]:
    """Newton-Raphson solve of the Michaelis-Menten oxygen sink variant.

    Solves, on the flattened (ravel'd) grid::

        D*lap(O) + supply*vasc*(1-O) - consumption*rho_tum*O/(O+K) = 0

    where ``source = supply*vasc`` and ``reaction_bulk = consumption*rho_tum``
    (both precomputed by the caller). Written as a residual in the same sign
    convention :func:`solve_oxygen`'s linear branch uses for its ``A @ O = b``::

        R(O) = (-D*L + diag(source)) @ O + reaction_bulk*O/(O+K) - source

    The first term is a FIXED linear operator (``M``, below) -- it does not
    depend on ``O`` -- so only the sink's derivative changes between iterations.
    Newton's method linearizes ``R`` about the current iterate::

        J(O) = M + diag(reaction_bulk * K / (O+K)^2)

    since ``d/dO [O/(O+K)] = K/(O+K)^2 > 0``, so ``J`` is ``M`` (already
    diagonally dominant / SPD-like, same operator the linear sink factors)
    plus a nonnegative diagonal -- a well-conditioned linear solve at every
    step, not new numerical machinery. Each step solves ``J @ delta = -R(O)``
    with the SAME sparse direct solver (`spsolve`) the linear sink uses, then
    updates ``O <- clip(O + step*delta, 0, 1)`` with a backtracking line
    search (halve the damping ``step`` until the residual norm decreases) so
    the iteration cannot diverge even from a poor start; clipping to [0, 1]
    keeps every intermediate iterate in the domain the O/(O+K) linearization
    is valid on. Seeded from ``O = 1`` (fully oxygenated; matches the source-
    only path the linear solve's RHS represents before any sink acts).

    Converges when ``max|delta| < tol``. This is standard Newton-Raphson for a
    smooth, monotone saturating nonlinearity -- quadratic local convergence, so
    it is expected (and, on the real flagship reconstruction, measured) to take
    a handful of iterations rather than the 15-40 Picard iterations the
    prototype needed for the same problem (``results/sweep/
    oxygen_functional_form_prototype_2026-07-29T05-54-00Z.md``). Raises
    ``RuntimeError`` rather than silently returning an unconverged field if
    ``max_iter`` is exhausted -- an unconverged nonlinear oxygen field is worse
    than an explicit failure.

    Returns
    -------
    (oxygen, n_iter, resid) : the flattened oxygen field, iterations used, and
    the final max|delta| (for diagnostics / the sensitivity-band report).
    """
    if K <= 0.0:
        raise ValueError(f"mm_half_saturation must be > 0, got {K}")
    if max_iter < 1:
        raise ValueError(f"mm_max_iter must be >= 1, got {max_iter}")

    # Fixed linear part of the operator: -D*L + diag(supply*vasc). Assembled
    # once and reused every Newton step -- only the sink's diagonal changes.
    M = (-D * L) + sp.diags(source, format="csr")

    def residual(o: np.ndarray) -> np.ndarray:
        return M @ o + reaction_bulk * o / (o + K) - source

    oxygen = np.ones_like(source)  # warm start: fully oxygenated
    resid = residual(oxygen)
    resid_norm = float(np.max(np.abs(resid))) if resid.size else 0.0

    n_iter = 0
    delta_norm = np.inf
    for n_iter in range(1, int(max_iter) + 1):
        deriv = reaction_bulk * K / (oxygen + K) ** 2
        J = M + sp.diags(deriv, format="csr")
        delta = spsolve(J.tocsc(), -resid)
        delta_norm = float(np.max(np.abs(delta))) if delta.size else 0.0

        # Backtracking line search: halve the step until the residual norm
        # actually decreases (or the step is negligibly small).
        step = 1.0
        while True:
            trial = np.clip(oxygen + step * delta, 0.0, 1.0)
            trial_resid = residual(trial)
            trial_norm = float(np.max(np.abs(trial_resid))) if trial_resid.size else 0.0
            if trial_norm <= resid_norm or step < 1e-4:
                break
            step *= 0.5

        oxygen, resid, resid_norm = trial, trial_resid, trial_norm

        if delta_norm * step < tol:
            break
    else:
        raise RuntimeError(
            f"Michaelis-Menten oxygen sink did not converge in {max_iter} Newton "
            f"iterations (K={K}, final max|delta|={delta_norm:.3e}, "
            f"max|residual|={resid_norm:.3e}); increase cfg.mm_max_iter or "
            f"loosen cfg.mm_tol."
        )

    return oxygen, n_iter, delta_norm


def solve_oxygen(micro: Microenvironment, cfg: OxygenConfig) -> np.ndarray:
    """Solve the steady-state reaction-diffusion oxygen field.

    Solves ``D*lap(O) + supply*vasc*(1-O) - sink(O) = 0`` on the domain grid
    with no-flux boundaries, where ``vasc`` is ``micro.fields.vascular_density``
    (the O2 source) and the sink is one of two functional forms selected by
    ``cfg.sink_form``:

    - ``"linear"`` (DEFAULT, unchanged behavior): ``sink(O) = consumption*rho_tum*O``,
      proportional to O everywhere. The resulting linear system

          (-D*L + diag(supply*vasc + consumption*rho_tum)) @ O = supply*vasc

      is symmetric positive-definite and solved directly with a single sparse
      solve -- no iteration, no CFL concern. This branch is BIT-IDENTICAL to the
      pre-existing code (same operations, same order) for any ``cfg`` that
      leaves ``sink_form`` at its default.
    - ``"michaelis_menten"`` (opt-in): ``sink(O) = consumption*rho_tum*O/(O+K)``
      where ``K = cfg.mm_half_saturation`` (Grimes, Fletcher & Partridge 2014,
      PMID 26064525). Nonlinear in O, so it is solved by Newton-Raphson
      iteration (:func:`_solve_michaelis_menten`) reusing the same sparse
      direct solver at every step. See the module docstring's "THE SINK'S
      FUNCTIONAL FORM..." section and ``results/oxygen_sensitivity_band.md``
      for why this is not a minor correction on the real flagship section.

    ``rho_tum`` is the local tumor-cell AREAL density from :func:`_tumor_density`
    -- cells per :data:`TUMOR_REFERENCE_AREA_UM2`, which is independent of the
    grid spacing -- for both sink forms.

    Parameters
    ----------
    micro : Microenvironment
        Provides the domain, ``fields.vascular_density`` and the cell population.
    cfg : OxygenConfig
        ``diffusivity`` (um^2/s), ``supply`` and ``consumption`` rate constants,
        plus ``sink_form`` / ``mm_half_saturation`` / ``mm_max_iter`` / ``mm_tol``
        (only consulted when ``sink_form == "michaelis_menten"``).

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
    sink_form = str(getattr(cfg, "sink_form", "linear")).strip().lower()
    if sink_form not in SINK_FORMS:
        raise ValueError(f"cfg.sink_form={sink_form!r} must be one of {SINK_FORMS}")

    L = _build_laplacian(ny, nx, domain.dx)

    source = (supply * vasc).ravel()                    # RHS b == supply*vasc
    reaction_bulk = (consumption * rho_tum).ravel()      # consumption*rho_tum

    # If there is no source and no tumor sink anywhere the system is singular
    # for EITHER sink form (both reduce to sink(O) == 0 identically); guard
    # once, identically, for both. (Equivalent, for "linear", to the original
    # guard on `not any(source + reaction_bulk)`, since both terms are >= 0.)
    if not np.any(source) and not np.any(reaction_bulk):
        return domain.zeros()

    if sink_form == "linear":
        reaction = source + reaction_bulk  # O-proportional loss
        # A = -D * L + diag(reaction). (-L is positive-definite; reaction >= 0.)
        A = (-D * L) + sp.diags(reaction, format="csr")
        oxygen = spsolve(A.tocsc(), source)
    else:
        K = float(cfg.mm_half_saturation)
        max_iter = int(cfg.mm_max_iter)
        tol = float(cfg.mm_tol)
        oxygen, _n_iter, _resid = _solve_michaelis_menten(
            L, D, source, reaction_bulk, K, max_iter, tol
        )

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
