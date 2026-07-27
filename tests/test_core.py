"""Tests for the frozen ``tmesim.core`` data model.

Covers Domain geometry/indexing, CellPopulation readouts (tumor_mask /
viable_tumor_fraction), Microenvironment.summary, and import stability of the
public core surface.
"""
from __future__ import annotations

import numpy as np
import pytest

import tmesim.core as core
from tmesim.core import (
    CellPopulation,
    CellState,
    CellType,
    Domain,
    Microenvironment,
    TissueFields,
)


# --------------------------------------------------------------------------- #
# Import stability
# --------------------------------------------------------------------------- #
def test_core_public_api_stable():
    """The documented core symbols must remain importable from tmesim.core."""
    expected = {
        "Domain",
        "CellPopulation", "CellState", "CellType", "Microenvironment",
        "TissueFields",
        "DRUG_LIBRARY", "DrugParams", "OxygenConfig", "PDConfig", "PDEConfig",
        "PKConfig", "SimConfig",
        "ExperimentResult", "PatientRecord", "SimulationResult", "TransportResult",
    }
    missing = {name for name in expected if not hasattr(core, name)}
    assert not missing, f"core is missing public symbols: {missing}"
    # __all__ should advertise the same surface.
    assert expected.issubset(set(core.__all__))


def test_celltype_cellstate_enum_values():
    """Enum ordinals are load-bearing (arrays store raw ints); pin them."""
    assert int(CellType.TUMOR) == 0
    assert int(CellState.VIABLE) == 0
    assert int(CellState.DEAD) == 3
    # Distinct values.
    assert len({int(c) for c in CellType}) == len(list(CellType))


# --------------------------------------------------------------------------- #
# Domain
# --------------------------------------------------------------------------- #
def test_domain_shape_and_extent():
    dom = Domain(nx=7, ny=5, dx=10.0)
    assert dom.shape == (5, 7)              # (ny, nx)
    assert dom.zeros().shape == (5, 7)
    assert dom.ones().shape == (5, 7)
    assert np.all(dom.zeros() == 0.0)
    assert np.all(dom.ones() == 1.0)
    assert dom.extent_um == (7 * 10.0, 5 * 10.0)


def test_domain_default_tissue_mask_all_true():
    dom = Domain(nx=4, ny=3)
    assert dom.tissue_mask.shape == (3, 4)
    assert dom.tissue_mask.dtype == bool
    assert dom.tissue_mask.all()


def test_domain_tissue_mask_shape_validation():
    with pytest.raises(ValueError):
        Domain(nx=4, ny=3, tissue_mask=np.ones((3, 3), dtype=bool))


def test_domain_coords_um_cell_centers():
    dom = Domain(nx=3, ny=2, dx=10.0)
    X, Y = dom.coords_um()
    assert X.shape == (2, 3) and Y.shape == (2, 3)
    # Cell centers are offset by half a cell.
    np.testing.assert_allclose(X[0], [5.0, 15.0, 25.0])
    np.testing.assert_allclose(Y[:, 0], [5.0, 15.0])


def test_domain_to_index_roundtrip_and_clamp():
    dom = Domain(nx=5, ny=4, dx=10.0)
    # A point inside cell (col=2, row=1): x in [20,30), y in [10,20).
    col, row = dom.to_index(np.array([25.0]), np.array([15.0]))
    assert (int(col[0]), int(row[0])) == (2, 1)
    # Out-of-range coordinates clamp into the grid.
    col, row = dom.to_index(np.array([-100.0, 1e6]), np.array([-100.0, 1e6]))
    assert list(col) == [0, dom.nx - 1]
    assert list(row) == [0, dom.ny - 1]


def test_domain_to_index_matches_coords():
    """coords_um cell centers must map back to their own indices."""
    dom = Domain(nx=6, ny=4, dx=15.0)
    X, Y = dom.coords_um()
    col, row = dom.to_index(X.ravel(), Y.ravel())
    exp_col = np.tile(np.arange(6), 4)
    exp_row = np.repeat(np.arange(4), 6)
    np.testing.assert_array_equal(col, exp_col)
    np.testing.assert_array_equal(row, exp_row)


# --------------------------------------------------------------------------- #
# CellPopulation
# --------------------------------------------------------------------------- #
def _make_pop():
    # 3 tumor, 2 stromal.
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    ctype = np.array([CellType.TUMOR, CellType.TUMOR, CellType.TUMOR,
                      CellType.STROMAL, CellType.STROMAL])
    state = np.full(5, int(CellState.VIABLE))
    sens = np.array([0.1, 0.5, 0.9, 0.2, 0.8])
    return CellPopulation(x=x, y=y, cell_type=ctype, state=state, sensitivity=sens)


def test_cellpopulation_basic_readouts():
    pop = _make_pop()
    assert pop.n == 5
    np.testing.assert_array_equal(pop.tumor_mask(), [True, True, True, False, False])
    assert pop.alive_mask().all()
    # exposure auto-initializes to zeros of length n.
    assert pop.exposure.shape == (5,)
    assert np.all(pop.exposure == 0.0)


def test_cellpopulation_length_validation():
    with pytest.raises(ValueError):
        CellPopulation(
            x=np.array([1.0, 2.0]), y=np.array([1.0]),
            cell_type=np.array([0, 0]), state=np.array([0, 0]),
            sensitivity=np.array([0.5, 0.5]),
        )


def test_viable_tumor_fraction_all_alive():
    pop = _make_pop()
    assert pop.viable_tumor_fraction() == 1.0


def test_viable_tumor_fraction_partial_kill():
    pop = _make_pop()
    # Kill one of the three tumor cells; killing a stromal cell must not count.
    pop.state[0] = int(CellState.DEAD)   # tumor
    pop.state[3] = int(CellState.DEAD)   # stromal (irrelevant to tumor fraction)
    assert pop.viable_tumor_fraction() == pytest.approx(2.0 / 3.0)


def test_viable_tumor_fraction_no_tumor_is_zero():
    pop = CellPopulation(
        x=np.array([1.0]), y=np.array([1.0]),
        cell_type=np.array([int(CellType.STROMAL)]),
        state=np.array([int(CellState.VIABLE)]),
        sensitivity=np.array([0.5]),
    )
    assert pop.viable_tumor_fraction() == 0.0


def test_cellpopulation_copy_is_independent():
    pop = _make_pop()
    dup = pop.copy()
    dup.state[0] = int(CellState.DEAD)
    # Mutating the copy must not touch the original.
    assert pop.state[0] == int(CellState.VIABLE)
    assert dup.state[0] == int(CellState.DEAD)


# --------------------------------------------------------------------------- #
# Microenvironment.summary
# --------------------------------------------------------------------------- #
def test_microenvironment_summary_fields():
    dom = Domain(nx=8, ny=6, dx=12.0)
    pop = _make_pop()
    fields = TissueFields(ecm_density=dom.ones() * 0.4,
                          vascular_density=dom.ones() * 0.25)
    micro = Microenvironment(domain=dom, fields=fields, cells=pop,
                             patient_id="P42", metadata={"subtype": "classical"})
    s = micro.summary()
    assert s["patient_id"] == "P42"
    assert s["grid"] == (6, 8)
    assert s["dx_um"] == 12.0
    assert s["n_cells"] == 5
    assert s["n_tumor"] == 3
    assert s["mean_ecm"] == pytest.approx(0.4)
    assert s["mean_vascular"] == pytest.approx(0.25)
    assert s["metadata_keys"] == ["subtype"]


def test_summary_on_real_microenvironment(real_micro):
    """summary() on a genuine REAL reconstructed microenvironment is well-formed."""
    s = real_micro.summary()
    assert s["n_cells"] > 0
    assert s["n_tumor"] > 0
    assert s["n_tumor"] <= s["n_cells"]
    assert 0.0 <= s["mean_ecm"] <= 1.0
    assert 0.0 <= s["mean_vascular"] <= 1.0
    assert s["grid"] == real_micro.domain.shape
