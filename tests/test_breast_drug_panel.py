"""Tests for the breast-cancer drug panel: the new marker map, the new drug
library / potency modules, and the compiled-drug pipeline
(``experiments/breast_drug_panel.py``).

REAL-ONLY, like the rest of the suite: everything here either needs no data
(the molecule library, MW self-checks) or skips cleanly when a real cache is
absent (the breast Xenium section, the breast GDSC cache) rather than
fabricating a substitute.
"""
from __future__ import annotations

import importlib.util
import math
import os

import numpy as np
import pytest

from tmesim.data.real.singlecell_spatial import load as load_xenium
from tmesim.molecular.breast_library import (
    BREAST_DRUGS,
    all_molecules_breast,
    get_molecule_breast,
)
from tmesim.molecular.potency_breast import resolve_ic50 as resolve_ic50_breast
from tmesim.reconstruct.marker_panels import (
    BRAIN_PANEL_MARKER_MAP,
    BREAST_PANEL_MARKER_MAP,
    IO_PANEL_MARKER_MAP,
    PANEL_REGISTRY,
    select_marker_map,
)

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))

BREAST_REP1 = os.path.join(_REPO, "data_cache", "xenium_breast_janesick", "rep1")
BRAIN_319 = os.path.join(_REPO, "data_cache", "xenium_gbm")
IO_480 = os.path.join(_REPO, "data_cache", "xenium_gbm_io")
BREAST_GDSC_CACHE = os.path.join(_REPO, "data_cache", "gdsc_breast_dose_response.csv")


def _skip_if_missing(path: str) -> None:
    if not os.path.exists(path):
        pytest.skip(f"real data not present at {path}")


def _load_breast_panel_module():
    path = os.path.join(_REPO, "experiments", "breast_drug_panel.py")
    spec = importlib.util.spec_from_file_location("_breast_drug_panel", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- #
# marker_panels.py: the additive breast registry entry
# --------------------------------------------------------------------------- #
class TestBreastMarkerPanel:
    def test_registered_in_panel_registry(self):
        names = [n for n, _ in PANEL_REGISTRY]
        assert "breast_313" in names
        assert names.index("brain_319") < names.index("breast_313"), (
            "brain_319 must stay first so it keeps winning ties against an "
            "unrecognised panel (the pre-existing behaviour)"
        )

    def test_breast_map_disjoint_content_from_brain_and_io(self):
        # Not a strict requirement, but the deliberate design: breast markers
        # are tissue-appropriate genes, not a copy of the brain/IO maps.
        assert BREAST_PANEL_MARKER_MAP != BRAIN_PANEL_MARKER_MAP
        assert BREAST_PANEL_MARKER_MAP != IO_PANEL_MARKER_MAP

    def test_gbm_panels_still_select_their_own_map_unaffected(self):
        """The load-bearing regression: adding a THIRD registry entry must not
        change which map wins on the two EXISTING real GBM panels."""
        _skip_if_missing(BRAIN_319)
        d = load_xenium(path=BRAIN_319)
        mm, name, cov = select_marker_map(d["genes"])
        assert name == "brain_319"
        assert mm is BRAIN_PANEL_MARKER_MAP
        assert cov["brain_319"] == 1.0

    def test_io_panel_still_selects_its_own_map_unaffected(self):
        _skip_if_missing(IO_480)
        d = load_xenium(path=IO_480)
        mm, name, cov = select_marker_map(d["genes"])
        assert name == "io_480"
        assert cov["io_480"] == 1.0

    def test_breast_panel_selects_the_breast_map_with_full_coverage(self):
        """The actual capability check: on the real 313-gene breast panel the
        new map should win, unambiguously (coverage 1.0), rather than falling
        back to the brain map by chance overlap (the pre-fix behaviour: brain
        0.583 / io 0.571 / no breast entry)."""
        _skip_if_missing(BREAST_REP1)
        d = load_xenium(path=BREAST_REP1)
        mm, name, cov = select_marker_map(d["genes"])
        assert name == "breast_313"
        assert cov["breast_313"] == 1.0
        assert cov["breast_313"] > cov["brain_319"]
        assert cov["breast_313"] > cov["io_480"]
        # Every gene in the map really is on this panel (not just claimed).
        genes = set(np.asarray(d["genes"], dtype=str).tolist())
        for markers in mm.values():
            for g in markers:
                assert g in genes, f"{g} claimed present but missing from panel"

    def test_breast_reconstruction_gives_nonzero_stromal_fraction(self):
        """The measured effect of the fix: under the brain-map fallback,
        STROMAL collapses to 0.0% on this panel (MOG/AQP4/GJA1 absent) and
        every non-immune, non-endothelial support cell silently defaults to
        TUMOR. With the breast map selected, STROMAL should be a real,
        nonzero fraction (fibroblast/myoepithelial markers ARE on this
        panel)."""
        _skip_if_missing(BREAST_REP1)
        from tmesim.reconstruct.real_xenium import phenotype_cells
        from tmesim.core import CellType

        d = load_xenium(path=BREAST_REP1)
        mm, name, _ = select_marker_map(d["genes"])
        assert name == "breast_313"
        cell_type, _, _ = phenotype_cells(d["expression"], d["genes"], marker_map=mm)
        n = len(cell_type)
        stromal_frac = float((np.asarray(cell_type) == int(CellType.STROMAL)).sum()) / n
        assert stromal_frac > 0.05, (
            f"STROMAL fraction {stromal_frac:.4f} should be a real double-digit-"
            "percent-scale population on breast tissue, not near-zero"
        )


# --------------------------------------------------------------------------- #
# breast_library.py: real, RDKit-verified molecules
# --------------------------------------------------------------------------- #
class TestBreastLibrary:
    _LIT_MW = {
        "paclitaxel": (850.0, 858.0),
        "cyclophosphamide": (259.0, 263.0),
        "tamoxifen": (369.0, 374.0),
    }

    @pytest.mark.parametrize("name", sorted(BREAST_DRUGS))
    def test_molecule_parses_and_matches_literature_mw(self, name):
        mol = get_molecule_breast(name)
        assert mol.source == "rdkit", f"{name} should parse via RDKit"
        lo, hi = self._LIT_MW[name]
        mw = mol.get("mw")
        assert lo <= mw <= hi, f"{name} MW {mw} outside literature band [{lo}, {hi}]"
        assert mol.get("is_antibody") is False

    def test_unknown_drug_raises_keyerror(self):
        with pytest.raises(KeyError):
            get_molecule_breast("trastuzumab")

    def test_disjoint_from_gbm_library(self):
        from tmesim.molecular.library import GBM_DRUGS

        assert set(BREAST_DRUGS).isdisjoint(set(GBM_DRUGS)), (
            "breast and GBM libraries should not silently overlap; doxorubicin "
            "is intentionally reused from the GBM library directly by name, "
            "not duplicated here"
        )

    def test_all_molecules_breast_builds_every_entry(self):
        mols = all_molecules_breast()
        assert {m.name for m in mols} == set(BREAST_DRUGS)


# --------------------------------------------------------------------------- #
# potency_breast.py: real GDSC2 breast (BRCA) IC50 resolution
# --------------------------------------------------------------------------- #
class TestPotencyBreast:
    @pytest.mark.parametrize("name", ["paclitaxel", "cyclophosphamide", "tamoxifen"])
    def test_measured_ic50_from_real_gdsc_breast_cache(self, name):
        _skip_if_missing(BREAST_GDSC_CACHE)
        res = resolve_ic50_breast(name)
        assert res["source"] == "GDSC"
        assert res["confidence"] == "measured"
        assert res["n_cell_lines"] >= 10
        assert res["ic50_uM"] > 0.0 and math.isfinite(res["ic50_uM"])

    def test_unknown_drug_falls_back_to_flagged_default(self):
        res = resolve_ic50_breast("totally_unknown_drug_xyz")
        assert res["confidence"] == "default_LOW"
        assert res["ic50_uM"] == 25.0

    def test_missing_cache_dir_degrades_to_literature_or_default(self):
        # A cache_dir with nothing in it must never raise.
        res = resolve_ic50_breast("paclitaxel", cache_dir="/nonexistent_dir_xyz")
        assert res["source"] in ("literature", "default")
        assert res["ic50_uM"] > 0.0


# --------------------------------------------------------------------------- #
# experiments/breast_drug_panel.py: the compiled-drug pipeline end to end
# --------------------------------------------------------------------------- #
class TestCompileBreastDrug:
    @pytest.fixture(scope="class")
    def mod(self):
        return _load_breast_panel_module()

    @pytest.fixture(scope="class")
    def alpha(self, mod):
        from tmesim.validation.calibration import calibrate_alpha

        return calibrate_alpha()

    def test_doxorubicin_reused_verbatim_from_gbm(self, mod, alpha):
        """doxorubicin must go through calibrated_drugparams("doxorubicin", ...)
        exactly as the GBM panel does -- same IC50, same everything except
        whatever the caller does with it downstream."""
        from tmesim.validation.calibration import calibrated_drugparams

        dp, mol, ic50_info = mod.compile_breast_drug("doxorubicin", alpha)
        dp_gbm = calibrated_drugparams("doxorubicin", alpha)
        assert dp.ic50 == dp_gbm.ic50
        assert dp.molecular_weight == dp_gbm.molecular_weight
        assert ic50_info["source"] == "literature"

    @pytest.mark.parametrize("name", ["paclitaxel", "cyclophosphamide", "tamoxifen"])
    def test_new_breast_drug_compiles_to_finite_sane_drugparams(self, mod, alpha, name):
        _skip_if_missing(BREAST_GDSC_CACHE)
        dp, mol, ic50_info = mod.compile_breast_drug(name, alpha)
        for fld in ("diffusivity", "molecular_weight", "ic50", "hill",
                    "binding", "clearance", "bbb_penetrance"):
            v = getattr(dp, fld)
            assert isinstance(v, (int, float)) and math.isfinite(v), f"{name}.{fld}={v}"
        assert dp.ic50 > 0.0
        assert 0.0 <= dp.bbb_penetrance <= 1.0
        assert ic50_info["confidence"] == "measured"
        assert ic50_info["source"] == "GDSC"

    def test_compiled_drugs_differ_in_ic50_ordering_as_expected(self, mod, alpha):
        """Sanity check on the real, not invented, IC50 ordering: paclitaxel
        (a highly potent taxane) should be far more potent (lower IC50) than
        cyclophosphamide (a weakly-cytotoxic-in-vitro pro-drug)."""
        _skip_if_missing(BREAST_GDSC_CACHE)
        dp_ptx, _, _ = mod.compile_breast_drug("paclitaxel", alpha)
        dp_cyc, _, _ = mod.compile_breast_drug("cyclophosphamide", alpha)
        assert dp_ptx.ic50 < dp_cyc.ic50


# --------------------------------------------------------------------------- #
# Full pipeline smoke test (skips cleanly without the real breast cache)
# --------------------------------------------------------------------------- #
class TestBreastPanelEndToEnd:
    def test_reconstruct_and_simulate_one_drug(self):
        """The step-3 de-risking check the task asked for: verify the EXISTING
        engine runs on breast tissue with zero core changes, and that
        kill_fraction / mechanism attribution reconcile exactly as they do for
        GBM. Runs at a coarse dx to stay fast in the test suite (the published
        panel numbers in results/breast_cancer_flagship.md use the converged
        dx=12um; this test only checks the pipeline WORKS, not the converged
        numbers)."""
        _skip_if_missing(BREAST_REP1)
        _skip_if_missing(BREAST_GDSC_CACHE)
        from tmesim.core import PDEConfig, SimConfig
        from tmesim.reconstruct.real_pipeline import (
            apply_drug_sensitivity,
            reconstruct_real,
        )
        from tmesim.sim.engine import simulate_response
        from tmesim.validation.calibration import calibrate_alpha

        mod = _load_breast_panel_module()
        micro = reconstruct_real(source="xenium", sample=BREAST_REP1, dx_um=40.0)
        assert micro.cells.n > 0
        assert micro.fields.vascular_density is not None
        assert micro.fields.ecm_density is not None

        alpha = calibrate_alpha()
        dp, mol, _ = mod.compile_breast_drug("paclitaxel", alpha)
        apply_drug_sensitivity(micro, "paclitaxel", mode="uniform")
        pde = PDEConfig(dt=10.0, total_time=3600.0, convection=True, save_every=500)
        cfg = SimConfig(drug=dp, pde=pde, seed=0)
        res = simulate_response(micro, cfg, attribute=True)

        assert 0.0 <= res.kill_fraction <= 1.0
        assert math.isfinite(res.kill_fraction)
        assert res.mechanism is not None
        assert res.mechanism["reconciliation"]["bit_exact"], (
            "mechanism attribution must reconcile bit-exact with kill_fraction "
            "on breast tissue exactly as it does for GBM -- same engine, zero "
            "core changes"
        )
