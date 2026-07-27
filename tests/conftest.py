"""Shared fixtures for the tmesim test suite — REAL-ONLY (no synthetic data).

The suite is deliberately fast. Two kinds of fixtures live here:

* **Real pipeline fixtures** — ``real_micro`` (a whole-niche microenvironment
  reconstructed from an open 10x Xenium GBM section via
  :func:`tmesim.reconstruct.real_pipeline.reconstruct_real`) and ``real_drug``
  (a calibrated + binding-site-barrier ``core.DrugParams`` for temozolomide).
  ``real_micro`` is session-scoped and built at a coarse grid so it assembles in
  a fraction of a second; ``simulate_response`` never mutates it (it works on a
  private copy), so it is safe to share across tests.

* **Physics-numerics builders** — ``flat_micro_factory`` / ``radial_micro_factory``
  and the ``inert_drug`` tracer. These construct ``Microenvironment`` objects
  directly from the frozen ``core`` containers (no synthetic patient generator)
  purely to exercise the transport / oxygen / pressure solvers on controlled
  geometry. They contain no invented "ground truth" or outcome labels.
"""
from __future__ import annotations

import numpy as np
import pytest

from tmesim.core import (
    CellPopulation,
    CellState,
    CellType,
    Domain,
    DrugParams,
    Microenvironment,
    TissueFields,
)
from tmesim.molecular import library
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.reconstruct.real_pipeline import reconstruct_real
from tmesim.validation.calibration import calibrated_drugparams

# Coarse grid keeps the real reconstruction (and any simulation on it) fast.
REAL_SOURCE = "xenium"
REAL_DX_UM = 80.0


# --------------------------------------------------------------------------- #
# Real pipeline fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="session")
def real_micro() -> Microenvironment:
    """A REAL whole-niche microenvironment from an open 10x Xenium GBM section.

    Built once per session at a coarse grid (dx=80 um) so it assembles in a
    fraction of a second. Real cells at real positions with REAL vascular / ECM /
    oxygen fields — no synthetic data anywhere in the path. Retries once on a
    transient load error.
    """
    try:
        return reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM)
    except Exception:
        import time

        time.sleep(3)
        return reconstruct_real(source=REAL_SOURCE, dx_um=REAL_DX_UM)


@pytest.fixture(scope="session")
def real_drug() -> DrugParams:
    """A REAL, engine-ready drug: temozolomide, calibrated (perivascular-sink
    alpha fit on the doxorubicin anchor) with the binding-site barrier applied.
    """
    dp = calibrated_drugparams("temozolomide")
    tmz = library.get_molecule("temozolomide")
    return apply_barrier_to_drugparams(dp, tmz)


# --------------------------------------------------------------------------- #
# Physics-numerics builders (controlled geometry, no synthetic ground truth)
# --------------------------------------------------------------------------- #
def make_flat_micro(ny: int = 40, nx: int = 40, dx: float = 15.0,
                    ecm_value: float = 0.3, tissue_mask=None) -> Microenvironment:
    """Minimal homogeneous microenvironment for transport numerics tests."""
    domain = Domain(nx=nx, ny=ny, dx=dx, tissue_mask=tissue_mask)
    fields = TissueFields(
        ecm_density=np.full((ny, nx), float(ecm_value)),
        vascular_density=np.zeros((ny, nx)),
    )
    cells = CellPopulation(
        x=np.array([nx * dx / 2]), y=np.array([ny * dx / 2]),
        cell_type=np.array([int(CellType.TUMOR)]),
        state=np.array([int(CellState.VIABLE)]),
        sensitivity=np.array([0.5]),
    )
    return Microenvironment(domain=domain, fields=fields, cells=cells)


def make_radial_tumor_micro(ny: int = 60, nx: int = 60, dx: float = 15.0,
                            n_tumor: int = 4000, seed: int = 0):
    """A central dense tumor with sparse-core vasculature + dense-core ECM.

    Mirrors the generative structure the physics self-tests rely on: sparse
    vessels in the core (drives hypoxia), leaky/elevated coupling toward the
    center (drives elevated core IFP). Returns (micro, r_field, core_radius).
    """
    rng = np.random.default_rng(seed)
    domain = Domain(nx=nx, ny=ny, dx=dx)
    X, Y = domain.coords_um()
    cx, cy = X.mean(), Y.mean()
    r = np.hypot(X - cx, Y - cy)
    core_radius = 0.30 * min(domain.extent_um)

    # Vasculature: perfused rim, suppressed avascular core.
    vasc = 0.6 - 0.59 * np.exp(-(r ** 2) / (2 * (0.5 * core_radius) ** 2))
    vasc = np.clip(vasc, 0.005, 1.0)
    # ECM: denser in the core.
    ecm = 0.2 + 0.5 * np.exp(-(r ** 2) / (2 * (0.7 * core_radius) ** 2))
    fields = TissueFields(ecm_density=ecm, vascular_density=vasc)

    ang = rng.uniform(0, 2 * np.pi, n_tumor)
    rad = core_radius * np.sqrt(rng.uniform(0, 1, n_tumor))
    tx = np.clip(cx + rad * np.cos(ang), 0, domain.extent_um[0] - 1e-3)
    ty = np.clip(cy + rad * np.sin(ang), 0, domain.extent_um[1] - 1e-3)
    cells = CellPopulation(
        x=tx, y=ty,
        cell_type=np.full(n_tumor, int(CellType.TUMOR), dtype=np.int8),
        state=np.full(n_tumor, int(CellState.VIABLE), dtype=np.int8),
        sensitivity=rng.uniform(0.2, 0.9, n_tumor),
    )
    micro = Microenvironment(domain=domain, fields=fields, cells=cells,
                             patient_id="radial-test")
    return micro, r, core_radius


@pytest.fixture
def flat_micro_factory():
    return make_flat_micro


@pytest.fixture
def radial_micro_factory():
    return make_radial_tumor_micro


@pytest.fixture
def inert_drug():
    """A conservative tracer: no clearance, no binding (mass must be conserved)."""
    return DrugParams(name="inert", diffusivity=300.0, clearance=0.0, binding=0.0)
