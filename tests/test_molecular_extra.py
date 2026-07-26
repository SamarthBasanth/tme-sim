"""Regression tests for the potency + ECM-mesh layers and their wiring into the
compile bridge (added after the initial molecular test suite)."""
import numpy as np
import pytest

from tmesim.molecular import ecm_mesh, potency, library
from tmesim.molecular.compile_drug import molecule_to_drugparams
from tmesim.core import DrugParams


# ----------------------------- potency (real IC50) -----------------------------

def test_resolve_ic50_shape_and_default():
    r = potency.resolve_ic50("temozolomide")
    assert set(["ic50_uM", "source", "confidence"]).issubset(r)
    assert r["ic50_uM"] > 0 and np.isfinite(r["ic50_uM"])
    # an unknown drug must degrade to a flagged default, never crash
    ru = potency.resolve_ic50("not_a_real_drug_xyz")
    assert ru["ic50_uM"] > 0
    assert "default" in ru["confidence"].lower() or "literature" in ru["confidence"].lower()


def test_temozolomide_ic50_is_measured_not_placeholder():
    # GDSC glioma median for TMZ is large (weak in vitro); must NOT be the old 25 uM stub
    r = potency.resolve_ic50("temozolomide")
    if r["source"] == "GDSC":
        assert r["ic50_uM"] > 100  # genuinely the measured weak-methylator value
        assert r.get("n_cell_lines", 0) >= 10


# ----------------------------- ECM mesh size gate ------------------------------

def test_ogston_partition_bounds_and_monotonic():
    small = ecm_mesh.ogston_renkin_hindrance(0.4)
    big = ecm_mesh.ogston_renkin_hindrance(5.5)
    for d in (small, big):
        assert 0.0 <= d["partition"] <= 1.0
    # larger solute -> lower partition and lower effective hindrance
    assert big["partition"] < small["partition"]
    assert big["effective_hindrance"] < small["effective_hindrance"]


def test_mesh_size_gate_is_strong():
    small = ecm_mesh.ogston_renkin_hindrance(0.4)["effective_hindrance"]
    antibody = ecm_mesh.ogston_renkin_hindrance(5.5)["effective_hindrance"]
    assert small > 5 * antibody  # strong size dependence (antibody excluded)


# ------------------------- integrated compile bridge ---------------------------

def test_compiled_drug_uses_measured_ic50():
    dp = molecule_to_drugparams(library.get_molecule("temozolomide"))
    assert isinstance(dp, DrugParams)
    # real potency flows through the bridge (measured/literature), not the 25 uM stub
    assert dp.ic50 > 100


def test_mesh_gate_lowers_antibody_diffusivity_far_below_small_molecule():
    tmz = molecule_to_drugparams(library.get_molecule("temozolomide"))
    bev = molecule_to_drugparams(library.get_molecule("bevacizumab"))
    # after the nm-scale mesh gate, the antibody's effective diffusivity collapses
    assert bev.diffusivity < tmz.diffusivity / 10
