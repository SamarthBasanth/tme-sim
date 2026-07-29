"""GBM-Space adapter: ``.h5ad`` -> Space Ranger dir -> the existing Visium
pipeline, UNCHANGED. REAL-ONLY, like the rest of the suite (see
``tests/test_ingest.py``'s module docstring for what that means here).

Two tiers of real data:

* :data:`_RAW_H5AD` — the ONE raw GBM-Space section (AT10-BRA-5-FO-3_1,
  ~154 MB) kept under ``data_cache/gbmspace_visium/_raw_h5ad/`` purely so
  :func:`inspect_section`'s h5py-decoding logic (categorical columns, the
  on-disk CSR, the obsm/spatial axis convention) is tested against genuine
  bytes rather than a hand-built fixture. Tests that need it are skipped if
  it is not present in a given checkout.
* :data:`_MATERIALIZED` — the 9 sections ingested and registered across two
  passes (3 in the first pass, 6 more in a follow-on pass that finished off
  every ``.h5ad`` member of the "selected" tarball), under
  ``data_cache/gbmspace_visium/<SAMPLE>/``. These are what
  ``results/gbmspace_adapter.md`` reports numbers from; the tests here assert
  those numbers are reproducible, not just that the code runs.

See ``results/gbmspace_ingestion.md`` (the format-mismatch discovery) and
``results/gbmspace_adapter.md`` (what was built on top of it, §1-8 the first
pass and §9 the follow-on) for the full account this adapter is the code for.
"""
from __future__ import annotations

import os

import numpy as np
import pytest
from scipy import sparse

from tmesim.data import ingest as ing
from tmesim.data.adapters import gbmspace
from tmesim.data.real import visium

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))

_RAW_H5AD = os.path.join(
    _REPO, "data_cache", "gbmspace_visium", "_raw_h5ad", "AT10-BRA-5-FO-3_1.h5ad"
)

#: (sample_name, materialized dir) for the 3 sections ingested in the first
#: pass, plus the 6 more added in the follow-on pass that finished off the
#: "selected" tarball (9 of 9 .h5ad members - see results/gbmspace_adapter.md
#: §9). Three of these (the AT15-BRA-4-FO-* ones) turned out to carry a
#: smaller ~18k-gene probe-based panel (Visium V5/CytAssist chemistry) rather
#: than the ~36.6k-gene whole-transcriptome panel the first 3 sections use -
#: the marker assertions below only check for specific genes, not panel size,
#: so this list is not split by chemistry.
_MATERIALIZED = [
    (sid, os.path.join(_REPO, "data_cache", "gbmspace_visium", sid))
    for sid in (
        "AT10-BRA-5-FO-3_1", "AT10-BRA-5-FO-4_2", "AT10-BRA-5-FO-2_1",
        "AT10-BRA-5-FO-1_1", "AT10-BRA-5-FO-1_2", "AT14-BRA-4-FO-2_1",
        "AT15-BRA-4-FO-C2-S27", "AT15-BRA-4-FO-C3-S34", "AT15-BRA-4-FO-B2-S25",
        # results/gbmspace_adapter.md §10: four more patients added to the
        # atlas (AT5, AT7, AT9, AT12), one representative section each.
        "AT12-BRA-4-FO-4_1", "AT9-BRA-4-FO-1_1", "AT7-BRA-3-FO-3_1",
        "AT5-BRA-5-FO-1_2",
    )
]

#: The whole point: markers the incumbent Xenium panels (319/480-gene) carry
#: none or almost none of (results/cell_state_impact.md), confirmed present on
#: GBM-Space's whole-transcriptome panel.
_HYPOXIA_MGMT_MARKERS = ["MGMT", "VEGFA", "CA9", "HIF1A", "ABCB1", "ABCG2"]


def _skip_missing(path: str) -> None:
    if not os.path.exists(path):
        pytest.skip(f"{path} not present in this checkout")


# --------------------------------------------------------------------------- #
# inspect_section(): the h5py-only .h5ad decoder
# --------------------------------------------------------------------------- #
def test_inspect_section_schema_and_markers():
    _skip_missing(_RAW_H5AD)
    section = gbmspace.inspect_section(_RAW_H5AD)

    assert section["sample_name"] == "AT10-BRA-5-FO-3_1"
    n_spots = len(section["barcodes"])
    assert n_spots > 1000
    assert len(set(section["barcodes"])) == n_spots  # unique after suffix strip

    # Whole-transcriptome panel, "Gene Expression" rows only (cell2location /
    # NMF columns dropped per the GBM-Space README's own instruction).
    assert section["n_genes_gene_expression"] == len(section["genes"])
    assert 36000 < section["n_genes_gene_expression"] < 36700
    assert section["n_genes_gene_expression"] < section["n_genes_total"]

    gene_set = set(section["genes"])
    for marker in _HYPOXIA_MGMT_MARKERS:
        assert marker in gene_set, f"{marker} missing from GBM-Space panel"

    expr = section["expression"]
    assert sparse.issparse(expr)
    assert expr.shape == (n_spots, len(section["genes"]))
    # Raw integer UMI counts, not the continuous cell2location/NMF columns.
    assert np.allclose(expr.data, np.round(expr.data))
    assert expr.data.min() >= 1.0  # CSR stores only explicit (nonzero) entries

    sf = section["scalefactors"]
    assert sf["spot_diameter_fullres"] > 0
    assert section["software_version"] == "spaceranger-2.0.0"

    # obsm/spatial axis convention: confirmed (not assumed) by correlation
    # against the independently-recorded array_row/array_col grid indices.
    corr_col = np.corrcoef(section["coords_px"][:, 0], section["array_col"])[0, 1]
    corr_row = np.corrcoef(section["coords_px"][:, 1], section["array_row"])[0, 1]
    assert abs(corr_col) > 0.999
    assert abs(corr_row) > 0.999


def test_inspect_section_rejects_a_non_gbmspace_h5ad(tmp_path):
    """A file with no var/feature_types group is refused, not silently misread."""
    import h5py

    bad = tmp_path / "not_gbmspace.h5ad"
    with h5py.File(bad, "w") as f:
        f.create_group("obs")
        f.create_group("var")
        f["var/_index"] = np.array([b"GENE1", b"GENE2"], dtype="S10")
    with pytest.raises(Exception):
        gbmspace.inspect_section(str(bad))


# --------------------------------------------------------------------------- #
# materialize_spaceranger_dir(): the adapter's actual output
# --------------------------------------------------------------------------- #
def test_materialize_then_detect_reaches_full_confidence(tmp_path):
    _skip_missing(_RAW_H5AD)
    section = gbmspace.inspect_section(_RAW_H5AD)
    dest = str(tmp_path / "materialized")
    gbmspace.materialize_spaceranger_dir(section, dest)

    assert os.path.isfile(os.path.join(dest, "spatial", "tissue_positions.csv"))
    assert os.path.isfile(os.path.join(dest, "spatial", "scalefactors_json.json"))
    matrix_dir = os.path.join(dest, "filtered_feature_bc_matrix")
    for name in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
        assert os.path.isfile(os.path.join(matrix_dir, name))

    # This is the actual claim: ingest.py's Visium detector needs ZERO new
    # code to recognise a materialized GBM-Space section.
    det = ing.detect(dest)
    assert det["platform"] == "visium"
    assert det["layout"] == "spaceranger_outs"
    assert det["confidence"] == 1.0, det["issues"]
    assert det["ingestable"]


def test_materialize_is_idempotent(tmp_path):
    _skip_missing(_RAW_H5AD)
    section = gbmspace.inspect_section(_RAW_H5AD)
    dest = str(tmp_path / "materialized")
    gbmspace.materialize_spaceranger_dir(section, dest)
    marker = os.path.join(dest, ".gbmspace_materialized")
    assert os.path.isfile(marker)
    mtime_before = os.path.getmtime(marker)

    gbmspace.materialize_spaceranger_dir(section, dest)  # no overwrite=True
    assert os.path.getmtime(marker) == mtime_before


def test_register_end_to_end_matches_inspect_section(tmp_path):
    """register() parity: the ingested SampleRecord's arrays trace back to
    exactly what inspect_section() read from the .h5ad, through the stock
    Visium loader - no second implementation to drift."""
    _skip_missing(_RAW_H5AD)
    section = gbmspace.inspect_section(_RAW_H5AD)

    out = gbmspace.register(
        _RAW_H5AD, sample_id="test_gbmspace_section",
        cache_root=str(tmp_path / "cache_root"),
    )
    rec = out["record"]
    assert rec.platform == "visium"
    assert rec.n_cells == len(section["barcodes"])
    assert rec.n_genes == len(section["genes"])
    assert set(rec.genes.tolist()) == set(section["genes"])
    assert rec.coordinate_units == "micrometer"
    assert rec.validation["ok"], rec.validation

    # Source dict is honest about the missing Ensembl ids (see module
    # docstring): symbol is written into both MEX id and symbol columns.
    assert "Ensembl" not in rec.source.get("panel", "") + str(rec.source)
    # sample_id resolves through the registry reconstruct_from_visium uses.
    assert "test_gbmspace_section" in visium.list_samples()
    spatial_dir, matrix_dir = visium._locate("test_gbmspace_section")
    assert os.path.abspath(spatial_dir) == os.path.join(out["dest_dir"], "spatial")
    visium.unregister_sample("test_gbmspace_section")


# --------------------------------------------------------------------------- #
# Capability report: hypoxia + MGMT measurable, per-cell routes still off
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("sample_name,path", _MATERIALIZED)
def test_materialized_sections_ingest_and_validate(sample_name, path):
    _skip_missing(path)
    det = ing.detect(path)
    assert det["platform"] == "visium"
    assert det["confidence"] == 1.0, det["issues"]

    rec = ing.ingest(path, sample_id=sample_name, detection=det)
    assert rec.validation["ok"], rec.validation
    # Threshold set against the actual spot-count distribution across all 20
    # sections on disk (780-12,120; results/gbmspace_adapter.md §10): AT5 and
    # AT9's smallest sections have ~800-900 spots, still a real tissue section,
    # not near-empty. 500 catches a genuinely broken/near-empty ingest without
    # excluding real small sections.
    assert rec.n_cells > 500
    # Two chemistries among the materialized sections (see the _MATERIALIZED
    # docstring above): ~36.6k-gene whole-transcriptome Visium V1, or the
    # ~18k-gene probe-based Visium V5/CytAssist panel. Either is fine; a
    # near-empty panel would not be.
    assert rec.n_genes > 15000

    cap = rec.capability
    # The point of this whole adapter: hypoxia AND MGMT are measurable on a
    # real GBM section for the first time in this repo, simultaneously.
    assert cap["hypoxia"]["present"], "expected VEGFA/CA9/... on the panel"
    assert cap["supports"]["hypoxia_measured"]["available"] is True
    assert cap["mgmt"]["present"] == ["MGMT"]
    assert cap["supports"]["mgmt_from_panel"]["available"] is True

    # Per-cell routes are still disabled BY CONSTRUCTION (spot-level data),
    # reusing the exact logic the 4 built-in Visium sections already hit -
    # a richer panel does not unlock what only per-cell resolution can.
    for name in ("state_mode:markers", "oxygen_mode:markers",
                 "sensitivity_mode:determinant"):
        assert cap["supports"][name]["available"] is False
        assert cap["supports"][name]["evidence"] == "platform"


def test_capability_confirms_full_marker_panel_on_first_section():
    sample_name, path = _MATERIALIZED[0]
    _skip_missing(path)
    genes = ing.panel_genes(path)
    assert len(genes) > 36000
    gene_set = set(genes)
    for marker in _HYPOXIA_MGMT_MARKERS + ["PECAM1", "MKI67", "TOP2A", "COL1A2"]:
        assert marker in gene_set


# --------------------------------------------------------------------------- #
# Reconstruction + a cheap simulate_response, end to end
# --------------------------------------------------------------------------- #
def test_reconstruct_from_visium_and_cheap_kill_run():
    sample_name, path = _MATERIALIZED[0]
    _skip_missing(path)

    from tmesim.core import (
        Microenvironment,
        OxygenConfig,
        PDEConfig,
        SimConfig,
        TissueFields,
    )
    from tmesim.cells.agents import state_fractions, update_cell_states
    from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
    from tmesim.molecular.library import get_molecule
    from tmesim.physics.oxygen import solve_oxygen
    from tmesim.reconstruct.real_pipeline import apply_drug_sensitivity
    from tmesim.reconstruct.real_visium import reconstruct_from_visium
    from tmesim.sim.engine import simulate_response
    from tmesim.validation.calibration import calibrated_drugparams

    reg_id = "test_reconstruct_gbmspace"
    visium.register_sample(reg_id, path, source="test", accession=sample_name,
                           diagnosis="test", open_access=True)
    try:
        micro = reconstruct_from_visium(reg_id)
        assert micro.cells.x.size > 1000

        # Same geometric-oxygen state-banding computation
        # results/real_drug_panel.md uses (private copy, no mutation).
        work = Microenvironment(
            domain=micro.domain,
            fields=TissueFields(ecm_density=micro.fields.ecm_density,
                                vascular_density=micro.fields.vascular_density),
            cells=micro.cells.copy(), patient_id=micro.patient_id,
            metadata=micro.metadata)
        update_cell_states(work, solve_oxygen(work, OxygenConfig()), OxygenConfig())
        states = state_fractions(work)
        assert set(states) == {"viable", "hypoxic", "quiescent", "dead"}
        total = sum(states.values())
        assert abs(total - 1.0) < 1e-9

        # A cheap drug run proves the whole transport -> kill engine accepts
        # this source, not just that the loader runs.
        dp = calibrated_drugparams("temozolomide", 148.6)
        mol = get_molecule("temozolomide")
        dp = apply_barrier_to_drugparams(dp, mol)
        apply_drug_sensitivity(micro, "temozolomide", mode="uniform")
        pde = PDEConfig(dt=10.0, total_time=600.0, convection=True, save_every=500)
        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
        assert 0.0 <= res.predicted_response <= 1.0
        assert 0.0 <= res.viable_tumor_fraction <= 1.0
    finally:
        visium.unregister_sample(reg_id)


# --------------------------------------------------------------------------- #
# Server wiring: the Studio section catalogue
# --------------------------------------------------------------------------- #
@pytest.fixture
def client(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient

    from app.server import cache as cache_mod
    from app.server import ingest_registry
    from app.server.api import create_app
    from app.server.jobs import JobManager

    runs = tmp_path / "runs"
    runs.mkdir()
    monkeypatch.setattr(cache_mod, "runs_dir", lambda: str(runs))
    monkeypatch.setattr(ingest_registry, "_store", None)
    monkeypatch.setattr(cache_mod, "_meta_cache", None)
    app = create_app(JobManager(root=str(runs / "jobs")))
    with TestClient(app) as c:
        yield c
    monkeypatch.setattr(ingest_registry, "_store", None)


def test_register_materialized_gbmspace_section_via_studio_api(client):
    """A materialized GBM-Space section registers exactly like any other
    Visium directory - ``ingest_registry`` does not know GBM-Space exists;
    it only sees a Space Ranger-shaped directory."""
    sample_name, path = _MATERIALIZED[0]
    _skip_missing(path)

    resp = client.post("/api/sections/scan", json={
        "path": path, "register": True,
        "section_id": "gbmspace_studio_test", "label": "GBM-Space (test)",
    })
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["ok"], body.get("error")
    assert body["section"]["platform"] == "visium"
    assert body["section"]["runnable"]
    assert body["capability"]["supports"]["hypoxia_measured"]["available"] is True
    assert body["capability"]["supports"]["mgmt_from_panel"]["available"] is True

    rows = client.get("/api/sections").json()
    ids = [r["id"] for r in rows]
    assert ids[:7] == ["xenium_319", "xenium_io", "gse324268",
                       "visium_10x_parent_gbm", "visium_gbm1", "visium_gbm2",
                       "visium_gbm3"]
    assert "gbmspace_studio_test" in ids

    # Reachable by reconstruct_from_visium, exactly like GBM1/2/3.
    assert "gbmspace_studio_test" in visium.list_samples()

    client.delete("/api/sections/gbmspace_studio_test")
    assert "gbmspace_studio_test" not in visium.list_samples()
