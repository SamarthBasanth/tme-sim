"""Shared fixtures and small-microenvironment builders for the tmesim test suite.

The whole suite is deliberately kept fast (small grids, short PDE windows). These
helpers build *real* pipeline objects (via ``data.synthetic`` and the ``core``
containers) so the tests exercise genuine signatures, not mocks.
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
    PDEConfig,
    TissueFields,
)
from tmesim.data.synthetic import generate_synthetic_patient

# Small, fast defaults used across the suite.
SMALL_GRID = (48, 48)
SMALL_DX = 15.0


def make_small_patient(seed: int = 1, grid=SMALL_GRID, dx: float = SMALL_DX,
                       subtype=None):
    """A small but fully-populated synthetic ``PatientRecord`` (fast to build)."""
    return generate_synthetic_patient(
        patient_id=f"SYN-T{seed:03d}", seed=seed, grid=grid, dx=dx, subtype=subtype)


def make_small_cohort(n: int = 6, seed: int = 0, grid=(40, 40), dx: float = SMALL_DX):
    """A small deterministic cohort spanning subtypes, on a small grid for speed.

    We build it directly (rather than via ``generate_synthetic_cohort``, which is
    hard-wired to a 120x120 grid) so the end-to-end pipeline tests stay well under
    the time budget.
    """
    subtypes = ("proneural", "classical", "mesenchymal")
    master = np.random.default_rng(seed)
    records = []
    for i in range(n):
        pseed = int(master.integers(0, 2**31 - 1))
        rec = generate_synthetic_patient(
            patient_id=f"SYN-{i:03d}", seed=pseed, grid=grid, dx=dx,
            subtype=subtypes[i % len(subtypes)])
        rec.micro.metadata["drug"] = "temozolomide-like"
        records.append(rec)
    return records


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
def small_patient():
    return make_small_patient(seed=1)


@pytest.fixture
def small_cohort():
    return make_small_cohort(n=6, seed=0)


# --------------------------------------------------------------------------- #
# Factory fixtures (import-mode independent access to the builders above)
# --------------------------------------------------------------------------- #
@pytest.fixture
def flat_micro_factory():
    return make_flat_micro


@pytest.fixture
def radial_micro_factory():
    return make_radial_tumor_micro


@pytest.fixture
def cohort_factory():
    return make_small_cohort


@pytest.fixture
def patient_factory():
    return make_small_patient


@pytest.fixture
def inert_drug():
    """A conservative tracer: no clearance, no binding (mass must be conserved)."""
    return DrugParams(name="inert", diffusivity=300.0, clearance=0.0, binding=0.0)


@pytest.fixture
def short_pde():
    return PDEConfig(dt=10.0, total_time=600.0, convection=False, save_every=10)
