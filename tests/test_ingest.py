"""Generic sample ingestion: detection, loader parity, capability, API wiring.

REAL-ONLY, like the rest of the suite. The two "generic layout" fixtures
(:func:`stock_xenium_bundle`, :func:`bare_mex_bundle`) contain no invented
numbers: they SYMLINK the real cached bytes into the directory shape a stock 10x
bundle uses, so what is being tested is the resolver and the generic reader
against real data in a layout this repo has never had on disk. That is the only
way to test "a new patient sample loads without editing code" without a new
patient sample.

The load-bearing test is :func:`test_ingest_parity_*`: the generic path must
reproduce ``singlecell_spatial.load`` and ``visium.load_sample`` bit-for-bit on
all seven cached specimens. If it cannot reproduce the specific loaders it is not
ready, and the built-in sections would silently move.
"""
from __future__ import annotations

import gzip
import json
import os

import numpy as np
import pytest
from scipy import sparse

from tmesim.data import ingest as ing
from tmesim.data.real import visium
from tmesim.data.real.singlecell_spatial import load as load_xenium

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(_HERE, ".."))


def _cache(*parts: str) -> str:
    return os.path.join(REPO, "data_cache", *parts)


#: (section id, cache dir) for the three cached Xenium sections.
XENIUM_SAMPLES = [
    ("xenium_319", _cache("xenium_gbm")),
    ("xenium_io", _cache("xenium_gbm_io")),
    ("gse324268", _cache("gse324268_xenium")),
]

#: (section id, registered visium sample id, cache dir) for the four Visium ones.
VISIUM_SAMPLES = [
    ("visium_10x_parent_gbm", "10x_parent_gbm", _cache("visium_gbm")),
    ("visium_gbm1", "GBM1", _cache("visium_gbm", "GBM1")),
    ("visium_gbm2", "GBM2", _cache("visium_gbm", "GBM2")),
    ("visium_gbm3", "GBM3", _cache("visium_gbm", "GBM3")),
]

ALL_PATHS = [p for _, p in XENIUM_SAMPLES] + [p for _, _, p in VISIUM_SAMPLES]


def _skip_missing(path: str) -> None:
    if not os.path.isdir(path):
        pytest.skip(f"cache {path} not present in this checkout")


def _same_matrix(a, b) -> bool:
    if sparse.issparse(a) != sparse.issparse(b):
        return False
    if not sparse.issparse(a):
        return np.array_equal(np.asarray(a), np.asarray(b))
    return a.shape == b.shape and (a != b).nnz == 0


# --------------------------------------------------------------------------- #
# detect()
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("sid,path", XENIUM_SAMPLES)
def test_detect_xenium_caches(sid, path):
    """All three Xenium caches are recognised at full confidence."""
    _skip_missing(path)
    det = ing.detect(path)
    assert det["platform"] == "xenium"
    assert det["layout"] in ("xenium_bundle_extracted", "xenium_bundle_tar")
    assert det["confidence"] == 1.0, det["issues"]
    assert det["ingestable"]
    assert det["files"]["cells"] and det["files"]["panel"]


@pytest.mark.parametrize("sid,vsample,path", VISIUM_SAMPLES)
def test_detect_visium_caches(sid, vsample, path):
    """All four Visium samples are recognised, with spatial/ and matrix/ paired."""
    _skip_missing(path)
    det = ing.detect(path)
    assert det["platform"] == "visium"
    assert det["layout"] == "spaceranger_outs"
    assert det["confidence"] == 1.0, det["issues"]
    assert det["ingestable"]
    # The two directories must be SIBLINGS, never pulled from different samples.
    assert (os.path.dirname(det["files"]["spatial_dir"])
            == os.path.dirname(det["files"]["matrix_dir"]))


def test_detect_resolves_repo_canonical_names_first():
    """gse324268 holds BOTH the canonical copy and the original GSM* spelling.

    The canonical copy is what every published number on that section was
    measured from, so the resolver must pick it, not the GSM one.
    """
    path = _cache("gse324268_xenium")
    _skip_missing(path)
    found = ing.resolve_xenium_files(path)
    assert os.path.basename(found["cells"]) == "Xenium_GBM_cells.csv.gz"
    assert os.path.basename(found["panel"]) == "Xenium_GBM_gene_panel.json"
    # Both spellings really are present, so this is a genuine choice.
    names = set(os.listdir(path))
    assert any(n.startswith("GSM") and n.endswith("_cells.csv.gz") for n in names)


def test_detect_container_is_not_a_sample():
    """A directory of samples is reported as a container, never guessed at.

    ``data_cache/`` has a ``spatial/`` under visium_gbm and a MEX triple under
    gse324268_xenium; pairing those two unrelated specimens into one "sample" is
    exactly the failure this rule exists to prevent.
    """
    path = _cache()
    _skip_missing(path)
    det = ing.detect(path)
    assert det["platform"] == "multi_sample"
    assert det["layout"] == "container"
    assert not det["ingestable"]
    names = {c["name"] for c in det["contained_samples"]}
    assert {"xenium_gbm", "visium_gbm"} <= names


def test_detect_visium_root_is_both_sample_and_container():
    """``data_cache/visium_gbm`` IS the 10x parent sample and holds GBM1/2/3."""
    path = _cache("visium_gbm")
    _skip_missing(path)
    det = ing.detect(path)
    assert det["platform"] == "visium" and det["ingestable"]
    assert {c["name"] for c in det["contained_samples"]} == {"GBM1", "GBM2", "GBM3"}


def test_detect_unknown_directory_explains_itself(tmp_path):
    det = ing.detect(str(tmp_path))
    assert det["platform"] == "unknown"
    assert det["confidence"] == 0.0
    assert not det["ingestable"]
    assert det["issues"] and "Xenium" in det["issues"][0]


def test_detect_nonexistent_path():
    det = ing.detect("/nonexistent/definitely/not/here")
    assert det["platform"] == "unknown"
    assert "does not exist" in det["issues"][0]


# --------------------------------------------------------------------------- #
# Parity with the specific loaders  -- THE test
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("sid,path", XENIUM_SAMPLES)
def test_ingest_parity_xenium(sid, path):
    """``ingest(dir)`` == ``singlecell_spatial.load(dir)``, array for array."""
    _skip_missing(path)
    ref = load_xenium(path=path, dense=False)
    rec = ing.ingest(path, sample_id=sid)

    assert rec.platform == "xenium"
    assert np.array_equal(ref["cell_xy_um"], rec.cell_xy_um)
    assert np.array_equal(np.asarray(ref["genes"], dtype=str), rec.genes)
    assert np.array_equal(np.asarray(ref["cell_ids"], dtype=str), rec.cell_ids)
    assert _same_matrix(ref["expression"], rec.expression)
    assert rec.coordinate_units == ref["coordinate_units"] == "micrometer"
    # Provenance follows the bytes, not the caller.
    for key in ("dataset", "accession", "sample", "assay", "panel",
                "n_panel_genes"):
        assert rec.source[key] == ref["source"][key], key


@pytest.mark.parametrize("sid,vsample,path", VISIUM_SAMPLES)
def test_ingest_parity_visium(sid, vsample, path):
    """``ingest(dir)`` == ``visium.load_sample(id)``, array for array."""
    _skip_missing(path)
    ref = visium.load_sample(vsample)
    rec = ing.ingest(path, sample_id=sid)

    assert rec.platform == "visium"
    assert np.array_equal(ref["coords_um"], rec.cell_xy_um)
    assert np.array_equal(np.asarray(ref["genes"], dtype=str), rec.genes)
    assert np.array_equal(np.asarray(ref["barcodes"], dtype=str), rec.cell_ids)
    assert _same_matrix(ref["expression"], rec.expression)
    # The extras the Visium loader carries survive under their original names.
    assert ref["umpp"] == rec.native["umpp"]
    assert ref["n_spots"] == rec.native["n_spots"]
    assert np.array_equal(ref["coords_px"], rec.native["coords_px"])
    assert np.array_equal(ref["total_counts"], rec.native["total_counts"])
    assert ref["markers_present"] == rec.native["markers_present"]
    assert os.path.abspath(ref["spatial_dir"]) == os.path.abspath(
        rec.native["spatial_dir"])


def test_visium_load_sample_still_delegates_to_load_from_dirs():
    """The registry loader IS the generic reader plus provenance, not a copy."""
    path = _cache("visium_gbm", "GBM3")
    _skip_missing(path)
    spatial_dir, matrix_dir = visium._locate("GBM3")
    a = visium.load_from_dirs(spatial_dir, matrix_dir, sample_id="GBM3")
    b = visium.load_sample("GBM3")
    assert np.array_equal(a["coords_um"], b["coords_um"])
    assert _same_matrix(a["expression"], b["expression"])
    # load_sample adds exactly the three registry-only fields.
    assert set(b) - set(a) == {"source", "accession", "diagnosis"}


# --------------------------------------------------------------------------- #
# Generic layouts: real bytes, in a shape this repo has never had on disk
# --------------------------------------------------------------------------- #
@pytest.fixture
def stock_xenium_bundle(tmp_path):
    """The real 319-panel section, renamed to STOCK 10x basenames.

    ``cells.csv.gz`` / ``cell_feature_matrix/`` / ``gene_panel.json`` is what a
    user downloading from 10x actually has. Before this module the loader
    required ``Xenium_GBM_*`` and this directory was unloadable.
    """
    src = _cache("xenium_gbm")
    _skip_missing(src)
    dest = tmp_path / "Patient_042_Xenium"
    dest.mkdir()
    os.symlink(os.path.join(src, "Xenium_GBM_cells.csv.gz"),
               dest / "cells.csv.gz")
    os.symlink(os.path.join(src, "Xenium_GBM_gene_panel.json"),
               dest / "gene_panel.json")
    os.symlink(os.path.join(src, "cell_feature_matrix"),
               dest / "cell_feature_matrix")
    return str(dest)


def test_stock_named_xenium_bundle_ingests(stock_xenium_bundle):
    det = ing.detect(stock_xenium_bundle)
    assert det["platform"] == "xenium"
    assert det["confidence"] == 1.0, det["issues"]

    rec = ing.ingest(stock_xenium_bundle)
    ref = load_xenium(path=_cache("xenium_gbm"), dense=False)
    assert np.array_equal(ref["cell_xy_um"], rec.cell_xy_um)
    assert _same_matrix(ref["expression"], rec.expression)
    # Provenance is read from the bundle's own gene_panel.json, so a renamed
    # copy of a known dataset is still identified as that dataset.
    assert rec.source["panel"] == ref["source"]["panel"]


@pytest.fixture
def bare_mex_bundle(tmp_path):
    """A real Visium GBM3 matrix laid out as a BARE MEX triple + coordinates.

    No ``spatial/``, no Space Ranger structure — the escape-hatch layout, built
    from real bytes so the generic reader is tested against real counts.
    """
    src = _cache("visium_gbm", "GBM3")
    _skip_missing(src)
    ref = visium.load_sample("GBM3")
    dest = tmp_path / "unknown_platform_sample"
    dest.mkdir()
    for name in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
        os.symlink(os.path.join(src, "filtered_feature_bc_matrix", name),
                   dest / name)
    # Coordinates for EVERY barcode in the matrix (an incomplete join is
    # refused, and that refusal is tested separately).
    with gzip.open(os.path.join(src, "filtered_feature_bc_matrix",
                                "barcodes.tsv.gz"), "rt") as fh:
        all_bc = [ln.strip() for ln in fh if ln.strip()]
    pos = {b: xy for b, xy in zip(ref["barcodes"], ref["coords_um"])}
    with open(dest / "coordinates.csv", "w") as fh:
        fh.write("cell_id,x_um,y_um\n")
        for b in all_bc:
            x, y = pos.get(b, (0.0, 0.0))
            fh.write(f"{b},{x},{y}\n")
    return str(dest), ref


def test_bare_mex_plus_coordinates_ingests(bare_mex_bundle):
    path, ref = bare_mex_bundle
    det = ing.detect(path)
    assert det["platform"] == "mex"
    assert det["ingestable"]
    assert det["files"]["coordinate_columns"][:2] == ["x_um", "y_um"]

    rec = ing.ingest(path)
    assert rec.platform == "mex"
    assert rec.n_genes == len(ref["genes"])
    # Coordinates are joined BY BARCODE, so the shared barcodes must land on
    # the same positions the Space Ranger reader put them at.
    index = {b: i for i, b in enumerate(rec.cell_ids)}
    rows = np.array([index[b] for b in ref["barcodes"]], dtype=int)
    assert np.allclose(rec.cell_xy_um[rows], ref["coords_um"])
    # The unit assumption is RECORDED, not silently made.
    assert "ASSUMED" in rec.coordinate_units
    assert rec.source["assumptions"]


def test_mex_refuses_incomplete_coordinate_join(bare_mex_bundle, tmp_path):
    """A partial coordinate table must be refused, not silently truncated."""
    path, ref = bare_mex_bundle
    dest = tmp_path / "partial"
    dest.mkdir()
    for name in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
        os.symlink(os.path.join(path, name), dest / name)
    with open(dest / "coordinates.csv", "w") as fh:
        fh.write("cell_id,x_um,y_um\n")
        for b, (x, y) in list(zip(ref["barcodes"], ref["coords_um"]))[:10]:
            fh.write(f"{b},{x},{y}\n")
    with pytest.raises(ing.IngestError, match="absent from"):
        ing.ingest(str(dest))


# --------------------------------------------------------------------------- #
# Validation
# --------------------------------------------------------------------------- #
def test_validation_passes_on_every_real_sample():
    for path in ALL_PATHS:
        if not os.path.isdir(path):
            continue
        det = ing.detect(path)
        assert det["ingestable"], path


def test_validation_refuses_pixel_coordinates(bare_mex_bundle, tmp_path):
    """Coordinates left in PIXELS give an absurd extent and must be refused."""
    path, ref = bare_mex_bundle
    dest = tmp_path / "pixels"
    dest.mkdir()
    for name in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
        os.symlink(os.path.join(path, name), dest / name)
    import csv as _csv

    with open(os.path.join(path, "coordinates.csv")) as fh:
        rows = list(_csv.DictReader(fh))
    with open(dest / "coordinates.csv", "w") as fh:
        fh.write("cell_id,x_um,y_um\n")
        for r in rows:  # ~50x, the um-per-pixel scale of a Visium fullres image
            fh.write(f"{r['cell_id']},{float(r['x_um']) * 50},"
                     f"{float(r['y_um']) * 50}\n")
    with pytest.raises(ing.IngestError, match="PIXELS or nanometres"):
        ing.ingest(str(dest))


def test_validation_refuses_a_degenerate_cell_count(bare_mex_bundle, tmp_path):
    path, ref = bare_mex_bundle
    rec = ing.ingest(path)
    tiny = ing.SampleRecord(
        sample_id="tiny", platform="mex", path=path,
        cell_xy_um=rec.cell_xy_um[:5], expression=rec.expression[:5],
        genes=rec.genes, cell_ids=rec.cell_ids[:5],
        coordinate_units="micrometer", source={},
    )
    verdict = ing.validate_record(tiny)
    assert not verdict["ok"]
    assert any("floor" in e for e in verdict["errors"])


def test_validation_catches_matrix_misalignment(bare_mex_bundle):
    path, _ref = bare_mex_bundle
    rec = ing.ingest(path)
    broken = ing.SampleRecord(
        sample_id="broken", platform="mex", path=path,
        cell_xy_um=rec.cell_xy_um, expression=rec.expression[:-1],
        genes=rec.genes, cell_ids=rec.cell_ids,
        coordinate_units="micrometer", source={},
    )
    verdict = ing.validate_record(broken)
    assert not verdict["ok"]
    assert any("not aligned to the coordinates" in e for e in verdict["errors"])


def test_validation_measures_what_it_checked():
    path = _cache("xenium_gbm")
    _skip_missing(path)
    rec = ing.ingest(path, sample_id="xenium_319")
    m = rec.validation["measured"]
    assert rec.validation["ok"]
    assert m["n_cells"] == 40887
    assert m["n_genes"] == 319
    assert 8000 < m["extent_um"][0] < 8200
    assert m["total_counts"] > 0


# --------------------------------------------------------------------------- #
# Capability reporting
# --------------------------------------------------------------------------- #
def test_capability_319_panel_reproduces_the_published_verdicts():
    """The 319-gene flagship panel fails BOTH state axes, and says so.

    This is ``results/cell_state_impact.md`` §1.1 and §5, stated BEFORE a run
    rather than discovered after quoting a number.
    """
    path = _cache("xenium_gbm")
    _skip_missing(path)
    rec = ing.ingest(path, sample_id="xenium_319")
    cap = rec.capability

    assert cap["panel"]["selected_map"] == "brain_319"
    assert cap["panel"]["selected_coverage"] == 1.0

    # Hypoxia: not a single gene of the set is on the panel.
    assert cap["hypoxia"]["present"] == []
    assert cap["supports"]["hypoxia_measured"]["available"] is False
    assert "UNMEASURED" in cap["supports"]["hypoxia_measured"]["reason"]

    # Proliferation: the markers ARE present, and fire far below the floor.
    assert set(cap["proliferation"]["cycle_specific_present"]) == {"MKI67", "TOP2A"}
    rate = cap["proliferation"]["detection_rate"]
    assert 0.005 < rate < 0.05, rate
    markers_route = cap["supports"]["state_mode:markers"]
    assert markers_route["available"] is False
    assert "below the 5.0% per-cell floor" in markers_route["reason"]
    assert f"{rate:.2%}" in markers_route["reason"]

    # Therefore no absolute kill, under any route.
    assert cap["supports"]["absolute_kill"]["available"] is False
    assert cap["supports"]["oxygen_mode:markers"]["available"] is False
    # MGMT is on neither Xenium panel.
    assert cap["supports"]["mgmt_from_panel"]["available"] is False
    # What DOES work: vasculature, typing, a real ECM field.
    assert cap["supports"]["vascular_field"]["available"] is True
    assert cap["supports"]["cell_typing"]["available"] is True
    assert cap["supports"]["ecm_field"]["available"] is True


def test_capability_io_panel_measures_hypoxia_but_loses_the_brain_markers():
    path = _cache("xenium_gbm_io")
    _skip_missing(path)
    genes = ing.panel_genes(path)
    assert len(genes) == 480
    cap = ing.capability_report(genes, platform="xenium", sample_id="xenium_io")

    assert cap["panel"]["selected_map"] == "io_480"
    assert cap["hypoxia"]["present"] == ["VEGFA"]
    assert cap["supports"]["hypoxia_measured"]["available"] is True
    # VEGFA is not CA9, and the report says so rather than implying specificity.
    assert "CA9" in cap["supports"]["hypoxia_measured"]["reason"]
    assert cap["supports"]["oxygen_mode:markers"]["available"] is True
    # The brain panel's endothelial and neural markers are gone; the IO map's
    # substitutes carry the classes instead.
    assert "PECAM1" not in cap["endothelial"]["present"]
    assert cap["cell_types"]["ENDOTHELIAL"]["present"] == ["VWF", "FLT1", "PLVAP"]
    assert not cap["cell_types"]["ENDOTHELIAL"]["missing"]


def test_capability_gse324268_is_the_only_panel_carrying_mgmt():
    path = _cache("gse324268_xenium")
    _skip_missing(path)
    cap = ing.capability_report(ing.panel_genes(path), platform="xenium",
                                sample_id="gse324268")
    assert cap["mgmt"]["present"] == ["MGMT"]
    assert cap["supports"]["mgmt_from_panel"]["available"] is True
    # And it is the only one of the three: the other two carry no MGMT.
    for other in (_cache("xenium_gbm"), _cache("xenium_gbm_io")):
        if os.path.isdir(other):
            assert "MGMT" not in ing.panel_genes(other)


def test_capability_visium_disables_per_cell_routes_by_construction():
    path = _cache("visium_gbm", "GBM1")
    _skip_missing(path)
    cap = ing.capability_report(ing.panel_genes(path), platform="visium",
                                sample_id="visium_gbm1")
    for name in ("state_mode:markers", "oxygen_mode:markers",
                 "sensitivity_mode:determinant"):
        assert cap["supports"][name]["available"] is False
        assert cap["supports"][name]["evidence"] == "platform"
    # The reason is the platform, NOT the panel: a whole transcriptome has the
    # genes; a spot is just not a cell.
    assert "SLC2A1" in cap["hypoxia"]["present"]


def test_capability_without_a_matrix_says_unmeasured_rather_than_guessing():
    path = _cache("xenium_gbm")
    _skip_missing(path)
    cap = ing.capability_report(ing.panel_genes(path), expression=None,
                                platform="xenium")
    assert cap["detection_rates_measured"] is False
    assert cap["proliferation"]["detection_rate"] is None
    route = cap["supports"]["state_mode:markers"]
    assert route["available"] is None           # unknown, not False
    assert "UNVERIFIED" in route["reason"]
    # Panel-level verdicts are still definite, because presence IS knowable.
    assert cap["supports"]["hypoxia_measured"]["available"] is False


def test_capability_flags_a_panel_no_map_fits():
    """A panel matching nothing must SAY nothing fits, not silently use brain."""
    cap = ing.capability_report(["ACTB", "GAPDH", "TUBB", "RPL13A"],
                                platform="xenium")
    assert cap["panel"]["fits"] is False
    assert "NO REGISTERED MARKER MAP FITS" in cap["panel"]["verdict"]
    # Every class is empty, including the one the oxygen source is built from.
    assert cap["supports"]["vascular_field"]["available"] is False
    assert "oxygen SOURCE" in cap["supports"]["cell_typing"]["reason"]


def test_capability_report_text_renders_the_no_lines():
    path = _cache("xenium_gbm")
    _skip_missing(path)
    rec = ing.ingest(path, sample_id="xenium_319")
    text = ing.format_capability_report(rec.capability)
    assert "CAPABILITY REPORT" in text
    assert "[ NO ] hypoxia_measured" in text
    assert "[ NO ] state_mode:markers" in text
    assert "[ OK ] vascular_field" in text


def test_scan_report_is_cheap_and_never_raises(tmp_path):
    out = ing.scan_report(str(tmp_path))
    assert out["capability"] is None
    assert out["error"]
    good = _cache("xenium_gbm")
    if os.path.isdir(good):
        out = ing.scan_report(good)
        assert out["capability"] is not None
        assert out["detection"]["platform"] == "xenium"
        assert out["capability"]["detection_rates_measured"] is False


# --------------------------------------------------------------------------- #
# Server wiring
# --------------------------------------------------------------------------- #
@pytest.fixture
def client(tmp_path, monkeypatch):
    from fastapi.testclient import TestClient

    from app.server import cache as cache_mod
    from app.server import ingest_registry
    from app.server.api import create_app
    from app.server.jobs import JobManager

    # Isolate the ingested-section store AND the metadata memo from the repo.
    runs = tmp_path / "runs"
    runs.mkdir()
    monkeypatch.setattr(cache_mod, "runs_dir", lambda: str(runs))
    monkeypatch.setattr(ingest_registry, "_store", None)
    monkeypatch.setattr(cache_mod, "_meta_cache", None)
    app = create_app(JobManager(root=str(runs / "jobs")))
    with TestClient(app) as c:
        yield c
    monkeypatch.setattr(ingest_registry, "_store", None)


def test_builtin_sections_are_unchanged_by_the_ingest_wiring(client):
    """The seven built-in sections, same ids, same order, still available."""
    rows = client.get("/api/sections").json()
    assert [r["id"] for r in rows] == [
        "xenium_319", "xenium_io", "gse324268", "visium_10x_parent_gbm",
        "visium_gbm1", "visium_gbm2", "visium_gbm3",
    ]
    for r in rows:
        assert r["available"], (r["id"], r.get("reason"))
        assert "ingested" not in r
    by_id = {r["id"]: r for r in rows}
    assert by_id["xenium_319"]["n_cells"] == 40887
    assert by_id["xenium_319"]["n_genes"] == 319


def test_scan_endpoint_reports_detection_and_capability(client):
    path = _cache("xenium_gbm")
    _skip_missing(path)
    resp = client.post("/api/sections/scan", json={"path": path})
    assert resp.status_code == 200
    body = resp.json()
    assert body["detection"]["platform"] == "xenium"
    assert body["detection"]["confidence"] == 1.0
    assert body["capability"]["supports"]["hypoxia_measured"]["available"] is False
    assert "[ NO ] hypoxia_measured" in body["capability_text"]


def test_scan_endpoint_on_a_bad_path_answers_rather_than_500s(client, tmp_path):
    resp = client.post("/api/sections/scan", json={"path": str(tmp_path)})
    assert resp.status_code == 200
    body = resp.json()
    assert body["detection"]["platform"] == "unknown"
    assert body["capability"] is None
    assert body["error"]


def test_scan_endpoint_requires_a_path(client):
    assert client.post("/api/sections/scan", json={}).status_code == 400


def test_register_ingested_section_and_run_the_registry(client, stock_xenium_bundle):
    """A directory that is in no table anywhere becomes a first-class section."""
    resp = client.post("/api/sections/scan", json={
        "path": stock_xenium_bundle, "register": True,
        "section_id": "patient_042", "label": "Patient 042 (ingested)",
    })
    assert resp.status_code == 200, resp.text
    body = resp.json()
    assert body["ok"], body.get("error")
    assert body["section_id"] == "patient_042"

    rows = client.get("/api/sections").json()
    ids = [r["id"] for r in rows]
    assert ids[:7] == ["xenium_319", "xenium_io", "gse324268",
                       "visium_10x_parent_gbm", "visium_gbm1", "visium_gbm2",
                       "visium_gbm3"]
    assert "patient_042" in ids

    row = next(r for r in rows if r["id"] == "patient_042")
    assert row["available"] and row["ingested"] and row["runnable"]
    assert row["n_cells"] == 40887 and row["n_genes"] == 319
    assert row["capability"]["supports"]["hypoxia_measured"]["available"] is False

    # It validates like any other section, with the SAME panel-derived refusals.
    verdict = client.post("/api/validate", json={
        "section": "patient_042", "dx_um": 30.0, "drug": "temozolomide",
        "oxygen_mode": "markers", "state_mode": "oxygen",
        "pk": {"plasma_peak": 1.0, "half_life_s": 14400.0, "n_doses": 1,
               "dose_interval_s": 86400.0},
        "total_time_s": 600.0,
    }).json()
    assert not verdict["ok"]
    assert any("hypoxia" in e for e in verdict["errors"])

    # And it can be removed.
    assert client.delete("/api/sections/patient_042").status_code == 200
    assert "patient_042" not in [r["id"] for r in client.get("/api/sections").json()]


def test_register_refuses_to_shadow_a_builtin_section(client, stock_xenium_bundle):
    resp = client.post("/api/sections/scan", json={
        "path": stock_xenium_bundle, "register": True,
        "section_id": "xenium_319",
    })
    body = resp.json()
    assert not body["ok"]
    assert "BUILT-IN" in body["error"]
    assert client.delete("/api/sections/xenium_319").status_code == 400


def test_register_refuses_an_uningestable_directory(client, tmp_path):
    body = client.post("/api/sections/scan",
                       json={"path": str(tmp_path), "register": True}).json()
    assert not body["ok"]
    assert body["error"]
    assert body["section"] is None


def test_register_ingested_visium_section(client):
    """An ingested Visium directory is reachable by reconstruct_from_visium."""
    path = _cache("visium_gbm", "GBM3")
    _skip_missing(path)
    body = client.post("/api/sections/scan", json={
        "path": path, "register": True, "section_id": "ingested_gbm3",
    }).json()
    assert body["ok"], body.get("error")
    assert body["section"]["source"] == "visium"
    # reconstruct_from_visium takes a sample id, so the id must now resolve.
    assert "ingested_gbm3" in visium.list_samples()
    spatial_dir, _m = visium._locate("ingested_gbm3")
    assert os.path.abspath(spatial_dir) == os.path.join(path, "spatial")
    client.delete("/api/sections/ingested_gbm3")
    assert "ingested_gbm3" not in visium.list_samples()


def test_ingested_section_survives_a_registry_reload(client, stock_xenium_bundle,
                                                     monkeypatch):
    from app.server import ingest_registry

    body = client.post("/api/sections/scan", json={
        "path": stock_xenium_bundle, "register": True, "section_id": "p43",
    }).json()
    assert body["ok"]
    monkeypatch.setattr(ingest_registry, "_store", None)  # forget in-memory
    assert "p43" in ingest_registry.ingested_sections()
    assert "p43" in [r["id"] for r in client.get("/api/sections").json()]
    client.delete("/api/sections/p43")
