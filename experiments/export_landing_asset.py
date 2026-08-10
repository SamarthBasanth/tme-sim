"""Export the REAL reconstructed GBM tissue section for the landing page's hero
canvas: the cells as a glTF-Binary point cloud
(``landing-3d/public/models/tumor.glb``) and the SOLVED interstitial fluid
pressure field over that same section as iso-contours
(``landing-3d/public/models/ifp_contours.json``).

Both outputs come from ONE reconstruction and share ONE coordinate transform,
which is the whole reason the field registers with the cells rather than merely
sitting near them.

WHY THIS EXISTS
---------------
The landing page shipped a procedurally-generated tumour (a distorted sphere,
``landing-3d/components/canvas/tumorGeometry.ts``) while its own copy claimed
the reconstruction uses no synthetic geometry. That claim was true of the
PRODUCT and false of the PICTURE. This script closes that gap: it writes the
actual reconstruction the rest of the repo runs on, so the thing a visitor
looks at is the thing the product measures.

WHAT IS REAL HERE
-----------------
Everything in the geometry. Every point is one segmented cell from a real
human glioblastoma section, at its MEASURED micrometre centroid:

    Dataset   10x Genomics Xenium In Situ, FFPE Human Brain Glioblastoma with
              Add-on (Xenium Analyzer v1.3.0)
    Sample id Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon
    Licence   CC BY 4.0, open direct download
    Contents  40,887 segmented cells, 319 gene targets (hBrain_254g + 65-gene
              add-on), 8.1 x 5.1 mm section
    Loader    tmesim.data.real.singlecell_spatial.load()
    Cache     data_cache/xenium_gbm  (see docs/DATA_INVENTORY.md section A1)

Positions come straight from ``cells.x`` / ``cells.y`` of the assembled
:class:`~tmesim.core.Microenvironment` produced by
:func:`tmesim.reconstruct.real_pipeline.reconstruct_real` — the same call every
headline result in ``results/`` is built on. No point is moved, jittered,
smoothed, or invented.

Per-vertex colour encodes CELL TYPE, which is a DERIVED measured quantity: the
threshold-gated marker typing the reconstruction performs on the 319-gene
panel (``matched_null_z`` aggregation, z-threshold 1.0 against a detection-rate
matched null; see ``micro.metadata['typing_note']``). The marker sets are

    TUMOR        EGFR, SOX2, OLIG2, MKI67
    IMMUNE       PTPRC, CD68, AIF1, CD163
    ENDOTHELIAL  PECAM1
    STROMAL      MOG, AQP4, GJA1

so a colour on the page means "this cell's transcript counts for that marker
set scored above the null", not "an artist decided this looked tumour-like".
NECROTIC is a fifth type in :class:`~tmesim.core.CellType` and is legitimately
absent from this section (0 cells), so it does not appear in the asset.

WHAT IS *NOT* REAL / HONEST CAVEATS
-----------------------------------
1. THE SECTION IS PLANAR. Xenium measures a single 2D tissue section; there is
   no measured z per cell. Every point is written at y = 0 EXACTLY. The
   vertical axis of the asset carries no data — it is the sectioning axis, and
   it is empty because nothing was measured along it. Nothing was extruded or
   given fake depth. Any 3D impression on the page comes from camera motion
   around a flat sheet of real cells, not from invented volume.
2. COLOUR HUES are the project's own canonical cell-type hues
   (``tmesim.viz.explainer.CELLTYPE_COLORS``, the same hues as every figure in
   ``results/``), but two display-only channels are adjusted for the dark hero:
   lightness is normalised to a single value so the four classes read at equal
   weight, and saturation is capped. HUE — i.e. which class is which — is
   untouched, so the asset and the figures agree.
3. THIS IS THE UNTREATED RECONSTRUCTION. ``cells.state`` is uniformly VIABLE
   (``state_mode='oxygen'`` leaves state alone), so no colour here encodes a
   drug response, a kill, or an outcome. This section has no outcome label.

COORDINATE TRANSFORM (documented so the page can state it)
----------------------------------------------------------
Micrometre section coordinates -> three.js scene units:

    s      = FIT_EXTENT / max(width_um, height_um)      # uniform, isotropic
    X      = (x_um - x_mid) * s
    Y      = 0
    Z      = (y_um - y_mid) * s

with FIT_EXTENT = 2.6 scene units, chosen to match the footprint of the
procedural placeholder it replaces (whose outermost shell radius was 1.34).
Uniform scale + centring only: no rotation, no shear, no non-uniform stretch,
so every distance in the asset is a true measured distance times one constant.
The exact micrometres-per-unit is written into the GLB's ``asset.extras`` so a
viewer can put a scale bar on it.

THE SOLVED IFP FIELD (second output)
------------------------------------
The landing page's section 02 is the beat that claims the simulator actually
solves something ("Then run physics on it."), and it used to be anchored by a
decorative contour ring whose own annotation read "Illustrative, not solved".
The honest fix is not to delete that disclaimer but to make the picture real,
which is what this half of the script does.

``tmesim.physics.pressure.solve_ifp`` is a genuine steady-state solve of

    -div( K grad p ) + (a + d) * p = a * p_v

on the SAME reconstructed microenvironment the point cloud is built from —
leaky-vessel Starling sources driven by the real measured ENDOTHELIAL cell
positions, lymphatic drainage suppressed inside the tumour, Dirichlet p = 0 at
the domain border and outside the tissue mask. Nothing about the field is
authored; it is a sparse direct solve on the 675 x 424 grid.

WHAT IS DISPLAYED, AND THE ONE DISPLAY CHOICE IN IT
---------------------------------------------------
The solved field's structure is much FINER than a pixel at the zoom the landing
page draws the section at. With the discrete vessel source, the screening
length sqrt(K/(a+d)) is of order 60 um, so p is a superposition of thousands of
small per-vessel haloes: thresholding the raw field yields ~1,400-3,000
disconnected components, each a fraction of a screen pixel across when 8.1 mm
of tissue is fitted into ~2.6 scene units. Drawn literally it is invisible
speckle, not a field.

So what is contoured is the field's LOCAL MEAN over a
:data:`IFP_COARSE_GRAIN_UM` micrometre Gaussian window, computed only over
tissue nodes (the mask is normalised out, so the average is not diluted by the
p = 0 exterior). This is a resolution choice about the picture, not a change to
the physics: no term is added, removed or retuned, and the averaging length is
written into the output and stated on the page. Coarse-graining is the ONLY
transformation applied between the solver and the contour vertices.

What survives that averaging is real regional structure — a single dominant
high-pressure lobe in the lower-left of the section, which sits on top of the
section's measured endothelial-density maximum (independently derived as
``ENTRY_POINTS[0]`` in the landing page's own shader constants). High IFP where
the vessels actually are is the model's central mechanism, and here it falls
out of the measured geometry rather than being placed by hand.

Levels are fixed fractions of the coarse-grained field's in-tissue maximum, so
each contour's colour is ``viridis(level / max)`` under the product's own ramp
table (``app/web/js/gl/ramps.js``, ``STOPS.viridis``) — the same ramp the Studio
paints scalar overlays with. The absolute mmHg value of every level is carried
in the JSON so the page can legend it in units.

OUTPUT
------
``landing-3d/public/models/tumor.glb`` — a valid binary glTF 2.0 containing one
root node with one child node per cell type, each a POINTS-mode primitive with
POSITION (float32 VEC3) and COLOR_0 (normalised uint8 VEC4). Written by hand
with numpy + struct + json; no glTF dependency is added to the project.
Provenance is embedded in ``asset.extras`` so the file is self-describing.

``landing-3d/public/models/ifp_contours.json`` — the solved IFP iso-contours,
outermost (lowest) level first, as polylines already in the page's scene
coordinates, plus the provenance needed to state what they are.

USAGE
-----
    .venv/bin/python experiments/export_landing_asset.py
    .venv/bin/python experiments/export_landing_asset.py --max-cells 20000
    .venv/bin/python experiments/export_landing_asset.py --no-glb   # IFP only
"""
import warnings; warnings.filterwarnings("ignore")

import argparse
import colorsys
import json
import os
import struct
import sys

import numpy as np
from scipy import ndimage
from skimage import measure

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:  # allow `python experiments/export_landing_asset.py`
    sys.path.insert(0, _ROOT)

from tmesim.core import CellState, CellType                      # noqa: E402
from tmesim.physics.pressure import solve_ifp                    # noqa: E402
from tmesim.reconstruct.real_pipeline import (                   # noqa: E402
    reconstruct_real,
    reconstruction_report,
)
from tmesim.viz.explainer import CELLTYPE_COLORS                 # noqa: E402

# --- reconstruction settings -------------------------------------------------
SAMPLE = "data_cache/xenium_gbm"
DX_UM = 12.0            # continuum grid spacing; affects the fields, not the cells
OXYGEN_MODE = "pde"     # historical default, bit-identical to every shipped result
STATE_MODE = "oxygen"   # leaves cells.state untouched -> all VIABLE

# --- asset settings ----------------------------------------------------------
OUT_REL = os.path.join("landing-3d", "public", "models", "tumor.glb")
FIT_EXTENT = 2.6        # scene units across the section's longest axis
DISPLAY_LIGHTNESS = 0.60  # HSL lightness all four cell-type hues are set to
DISPLAY_MAX_SATURATION = 0.60  # ceiling, so a fully-saturated print hue is not neon

CITATION = (
    "10x Genomics Xenium In Situ, FFPE Human Brain Glioblastoma with Add-on "
    "(Xenium Analyzer v1.3.0), sample Xenium_V1_FFPE_Human_Brain_Glioblastoma_"
    "With_Addon, CC BY 4.0."
)

# --- solved-IFP contour settings ---------------------------------------------
IFP_OUT_REL = os.path.join("landing-3d", "public", "models", "ifp_contours.json")

#: Gaussian averaging length applied to the SOLVED field before contouring, in
#: micrometres. The only transformation between solver and contour vertices;
#: see the module docstring for why it is required (the raw field's ~60 um
#: per-vessel structure is far below one screen pixel at the page's zoom).
IFP_COARSE_GRAIN_UM = 250.0

#: Iso-levels as fractions of the coarse-grained field's in-tissue maximum,
#: outermost (lowest) first. Fractions rather than round mmHg values so that
#: each level's colour is exactly ``viridis(fraction)`` — the ramp position and
#: the level are then the same number, and the legend cannot drift from the
#: picture. Spread up the ramp rather than bunched at the top so the drawn set
#: spans a visible range of viridis instead of four shades of yellow.
IFP_LEVEL_FRACTIONS = (0.30, 0.50, 0.70, 0.90)

#: Drop contours shorter than this many marching-squares vertices. At dx = 12 um
#: a 40-vertex loop is roughly a 150 um feature — below the averaging length, so
#: anything smaller is a ripple in the window, not a resolved structure.
IFP_MIN_CONTOUR_VERTICES = 40

#: Douglas-Peucker tolerance for polyline simplification, in SCENE units.
#: 0.0035 units is ~11 um, i.e. under one grid cell, so simplification cannot
#: move a vertex by as much as the grid the field was solved on.
IFP_SIMPLIFY_TOLERANCE = 0.0035

#: matplotlib viridis control points, copied from the PRODUCT's own ramp table
#: (``app/web/js/gl/ramps.js``, ``STOPS.viridis``) so a level's colour here is
#: the colour the Studio would paint that same normalised value. Kept as the
#: full 17 stops; the landing page's GLSL copy subsamples 9 of them, which
#: agrees with this to within a couple of 1/255 steps.
VIRIDIS_STOPS_255 = (
    (68, 1, 84), (72, 24, 106), (71, 45, 123), (66, 64, 134), (59, 82, 139),
    (51, 99, 141), (44, 114, 142), (38, 130, 142), (33, 145, 140),
    (31, 160, 136), (40, 174, 128), (63, 188, 115), (94, 201, 98),
    (132, 212, 75), (173, 220, 48), (216, 226, 25), (253, 231, 37),
)

# glTF / GL constants
GL_FLOAT = 5126
GL_UNSIGNED_BYTE = 5121
GL_ARRAY_BUFFER = 34962
GL_POINTS = 0


# ---------------------------------------------------------------- palette ----
def _hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _srgb_to_linear(c):
    return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4


def display_palette():
    """Canonical cell-type hues, lightness-normalised for the dark hero.

    Returns ``{cell_type_int: (hex_srgb, (r, g, b) linear-light floats)}``.
    HUE is taken verbatim from :data:`tmesim.viz.explainer.CELLTYPE_COLORS`, so
    the landing page and every figure in ``results/`` agree on which colour
    means which cell type. Only the two display-only channels move: HSL
    lightness is set to :data:`DISPLAY_LIGHTNESS` (equal weight — no class
    shouts louder than another against ``--bg: #07080a``) and saturation is
    capped at :data:`DISPLAY_MAX_SATURATION` (STROMAL's print hue #008300 is
    fully saturated and would otherwise lift to neon green).

    The linear-light triple is what goes into the GLB: glTF ``COLOR_0`` is
    defined in linear space, and three.js (r152+) treats it as already being in
    the renderer's linear working space.
    """
    out = {}
    for ct, hex_src in CELLTYPE_COLORS.items():
        r, g, b = _hex_to_rgb(hex_src)
        h, _l, s = colorsys.rgb_to_hls(r, g, b)
        r2, g2, b2 = colorsys.hls_to_rgb(h, DISPLAY_LIGHTNESS,
                                         min(s, DISPLAY_MAX_SATURATION))
        hex_disp = "#%02x%02x%02x" % tuple(round(c * 255) for c in (r2, g2, b2))
        out[int(ct)] = (hex_disp, tuple(_srgb_to_linear(c) for c in (r2, g2, b2)))
    return out


# -------------------------------------------------------------- glb writer ---
def _pad4(buf, fill=b"\x00"):
    rem = (-len(buf)) % 4
    return buf + fill * rem


def write_glb(path, groups, extras):
    """Write a POINTS-mode binary glTF 2.0.

    ``groups`` is a list of ``(name, positions (n,3) float32, colors (n,4)
    uint8, node_extras)``. Each becomes its own node/mesh so the page can build
    a legend or toggle classes without re-parsing the buffer. One shared
    material; three.js's GLTFLoader converts it to a ``PointsMaterial`` and
    enables ``vertexColors`` because ``COLOR_0`` is present.

    Hand-rolled rather than pulling in trimesh/pygltflib: a point cloud needs
    two accessors per group and no scene graph to speak of, and this is a
    one-off asset build, not a runtime dependency.
    """
    bin_parts, offset = [], 0
    buffer_views, accessors, meshes, nodes = [], [], [], []

    for name, pos, col in ((g[0], g[1], g[2]) for g in groups):
        pos = np.ascontiguousarray(pos, dtype="<f4")
        col = np.ascontiguousarray(col, dtype=np.uint8)
        n = int(pos.shape[0])

        # POSITION — 12 B/element, inherently 4-byte aligned.
        pos_bytes = pos.tobytes()
        buffer_views.append({"buffer": 0, "byteOffset": offset,
                             "byteLength": len(pos_bytes), "target": GL_ARRAY_BUFFER})
        bin_parts.append(pos_bytes)
        offset += len(pos_bytes)
        acc_pos = len(accessors)
        accessors.append({
            "bufferView": len(buffer_views) - 1, "componentType": GL_FLOAT,
            "count": n, "type": "VEC3",
            "min": [float(v) for v in pos.min(axis=0)],
            "max": [float(v) for v in pos.max(axis=0)],
        })

        # COLOR_0 — VEC4/UNSIGNED_BYTE is 4 B/element. VEC3 would be 3 B and
        # glTF requires every vertex-attribute element to be 4-byte aligned,
        # so the alpha channel is not optional here.
        col_bytes = col.tobytes()
        buffer_views.append({"buffer": 0, "byteOffset": offset,
                             "byteLength": len(col_bytes), "target": GL_ARRAY_BUFFER})
        bin_parts.append(col_bytes)
        offset += len(col_bytes)
        acc_col = len(accessors)
        accessors.append({
            "bufferView": len(buffer_views) - 1, "componentType": GL_UNSIGNED_BYTE,
            "count": n, "type": "VEC4", "normalized": True,
        })

        meshes.append({"name": name, "primitives": [{
            "attributes": {"POSITION": acc_pos, "COLOR_0": acc_col},
            "mode": GL_POINTS, "material": 0,
        }]})

    # node 0 is the root; children follow. Deliberately avoids the substrings
    # "core"/"edema"/"margin"/"tissue", which GLTFTumor.tsx matches on to drive
    # the placeholder's split animation — these are cell types, not regions.
    nodes.append({"name": "xenium_gbm_section",
                  "children": list(range(1, len(groups) + 1))})
    for i, g in enumerate(groups):
        nodes.append({"name": g[0], "mesh": i, "extras": g[3]})

    bin_blob = _pad4(b"".join(bin_parts))
    gltf = {
        "asset": {"version": "2.0", "generator": "tme-sim experiments/export_landing_asset.py",
                  "extras": extras},
        "scene": 0,
        "scenes": [{"nodes": [0]}],
        "nodes": nodes,
        "meshes": meshes,
        "materials": [{
            "name": "cell_points",
            "pbrMetallicRoughness": {"baseColorFactor": [1.0, 1.0, 1.0, 1.0],
                                     "metallicFactor": 0.0, "roughnessFactor": 1.0},
            "doubleSided": True,
        }],
        "accessors": accessors,
        "bufferViews": buffer_views,
        "buffers": [{"byteLength": len(bin_blob)}],
    }

    json_blob = _pad4(json.dumps(gltf, separators=(",", ":")).encode("utf-8"), b" ")
    total = 12 + 8 + len(json_blob) + 8 + len(bin_blob)
    with open(path, "wb") as fh:
        fh.write(struct.pack("<III", 0x46546C67, 2, total))       # 'glTF', v2
        fh.write(struct.pack("<II", len(json_blob), 0x4E4F534A))  # 'JSON'
        fh.write(json_blob)
        fh.write(struct.pack("<II", len(bin_blob), 0x004E4942))   # 'BIN\0'
        fh.write(bin_blob)
    return total


# ------------------------------------------------------------ solved IFP ----
def viridis_hex(t):
    """``viridis(t)`` as an sRGB hex string, from the product's own ramp table.

    Piecewise-linear between :data:`VIRIDIS_STOPS_255`, which is the same
    interpolation ``app/web/js/gl/ramps.js`` does when it bakes its 256-entry
    LUT, so a level drawn on the landing page and the same normalised value
    painted in the Studio come out the same colour.
    """
    t = min(max(float(t), 0.0), 1.0)
    n = len(VIRIDIS_STOPS_255) - 1
    u = t * n
    i0 = min(n, int(np.floor(u)))
    i1 = min(n, i0 + 1)
    f = u - i0
    c = [VIRIDIS_STOPS_255[i0][k] * (1 - f) + VIRIDIS_STOPS_255[i1][k] * f
         for k in range(3)]
    return "#%02x%02x%02x" % tuple(int(round(v)) for v in c)


def coarse_grain(field, mask, sigma_px):
    """Local mean of ``field`` over tissue nodes only, Gaussian window.

    A plain ``gaussian_filter`` would average the p = 0 exterior in and pull
    every near-boundary value down, inventing a pressure falloff at the section
    edge that the solve does not have. Dividing the smoothed field by the
    smoothed mask (Nadaraya-Watson / "normalised convolution") makes this a
    weighted mean over tissue nodes alone, so the value at a boundary node is
    the mean of the tissue actually around it.
    """
    num = ndimage.gaussian_filter(np.where(mask, field, 0.0), sigma_px, mode="nearest")
    den = ndimage.gaussian_filter(mask.astype(np.float64), sigma_px, mode="nearest")
    return np.where(den > 1e-3, num / np.maximum(den, 1e-9), 0.0)


def export_ifp_contours(path, micro, x_mid, y_mid, scale, sample, dx_um):
    """Solve interstitial fluid pressure and write its iso-contours as JSON.

    ``x_mid`` / ``y_mid`` / ``scale`` are the SAME transform constants the point
    cloud was written with — passed in rather than recomputed so the two
    outputs cannot drift apart. Contour vertices come out in the page's
    standing-up frame: the exporter lays the section in the glTF XZ plane with
    y == 0 and the page turns it -90 degrees about X, which sends the asset's z
    to world y, so a scene-space ``(x, y)`` here is ``(X, Z)`` there.
    """
    p, _velocity = solve_ifp(micro)
    mask = np.asarray(micro.domain.tissue_mask, dtype=bool)
    dx = float(micro.domain.dx)

    raw_max = float(p[mask].max())
    sigma_px = IFP_COARSE_GRAIN_UM / dx
    ps = coarse_grain(p, mask, sigma_px)
    disp_max = float(ps[mask].max())
    if not np.isfinite(disp_max) or disp_max <= 0.0:
        raise RuntimeError("solved IFP field is empty in tissue; refusing to export")

    levels = []
    for frac in IFP_LEVEL_FRACTIONS:
        value = frac * disp_max
        paths = []
        for contour in measure.find_contours(ps, value):
            if len(contour) < IFP_MIN_CONTOUR_VERTICES:
                continue
            # find_contours returns (row, col) in array-index space; array
            # element (r, c) is the grid cell CENTRE, at ((c+0.5)*dx, (r+0.5)*dx)
            # micrometres, which is Domain.coords_um's convention.
            x_um = (contour[:, 1] + 0.5) * dx
            y_um = (contour[:, 0] + 0.5) * dx
            xy = np.stack([(x_um - x_mid) * scale, (y_um - y_mid) * scale], axis=1)
            closed = bool(np.allclose(xy[0], xy[-1]))
            xy = measure.approximate_polygon(xy, IFP_SIMPLIFY_TOLERANCE)
            if len(xy) < 3:
                continue
            if closed:  # approximate_polygon can drop the repeated last vertex
                xy = np.vstack([xy, xy[0]])
            paths.append([[round(float(a), 5), round(float(b), 5)] for a, b in xy])
        if not paths:
            continue
        levels.append({
            "fraction_of_max": round(float(frac), 4),
            "value_mmhg": round(value, 5),
            # colour == viridis(fraction_of_max) by construction, so the ramp
            # position and the iso-level are the same number.
            "color_srgb_hex": viridis_hex(frac),
            "n_paths": len(paths),
            "n_vertices": int(sum(len(q) for q in paths)),
            "paths": paths,
        })
    if not levels:
        raise RuntimeError("no IFP contour survived; refusing to export")

    doc = {
        "what_this_is": (
            "Iso-contours of a SOLVED interstitial fluid pressure field over the "
            "same real Xenium GBM section as tumor.glb. Not illustrative: the "
            "field is the output of tmesim.physics.pressure.solve_ifp on the "
            "reconstruction, and the only transformation applied before "
            "contouring is a stated local mean (see display_coarse_grain_um)."
        ),
        "quantity": "interstitial fluid pressure (IFP)",
        "units": "mmHg",
        "solver": "tmesim.physics.pressure.solve_ifp",
        "model": ("steady-state screened Poisson  -div(K grad p) + (a+d) p = a p_v  "
                  "with leaky-vessel Starling sources from the measured "
                  "ENDOTHELIAL cell positions, lymphatic drainage suppressed "
                  "inside the tumour, Dirichlet p=0 at the border and outside "
                  "the tissue mask; 5-point finite volume, harmonic-mean face "
                  "conductances, sparse direct solve"),
        "is_solved": True,
        "is_calibrated_to_patient": False,
        "parameter_note": ("Transport parameters are order-of-magnitude "
                           "biological values, not patient-calibrated (see "
                           "tmesim/physics/pressure.py). The GEOMETRY driving "
                           "the solve is measured; the coefficients are not "
                           "fitted to this section."),
        "dataset": "10x Genomics Xenium Human Brain — Glioblastoma (v1.3.0)",
        "citation": CITATION,
        "licence": "CC BY 4.0",
        "sample_id": "Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon",
        "pipeline": ("tmesim.reconstruct.real_pipeline.reconstruct_real"
                     f"(source='xenium', sample='{sample}', dx_um={dx_um}, "
                     f"oxygen_mode='{OXYGEN_MODE}', state_mode='{STATE_MODE}') "
                     "-> tmesim.physics.pressure.solve_ifp"),
        "patient_id": micro.patient_id,
        "grid": [int(micro.domain.ny), int(micro.domain.nx)],
        "dx_um": dx,
        "raw_field_max_mmhg": round(raw_max, 5),
        "raw_field_mean_in_tissue_mmhg": round(float(p[mask].mean()), 6),
        "display_coarse_grain_um": IFP_COARSE_GRAIN_UM,
        "display_coarse_grain_note": (
            "Contours are of the solved field's local mean over a "
            f"{IFP_COARSE_GRAIN_UM:g} um Gaussian window taken across tissue "
            "nodes only. The raw field's structure is per-vessel (screening "
            "length ~60 um) and lies below one screen pixel at the zoom the "
            "page draws the whole 8.1 mm section at, so it is averaged for "
            "display. No solver term was changed. Because this is a MEAN, the "
            "level values below are much lower than raw_field_max_mmhg."
        ),
        "displayed_field_max_mmhg": round(disp_max, 5),
        "ramp": "viridis",
        "ramp_source": "app/web/js/gl/ramps.js STOPS.viridis (the product's own table)",
        "ramp_mapping": "color = viridis(value_mmhg / displayed_field_max_mmhg)",
        "um_per_scene_unit": round(1.0 / scale, 2),
        "transform": ("Identical to tumor.glb's, then expressed in the page's "
                      "standing-up frame: scene x=(x_um-x_mid)*s, "
                      f"scene y=(y_um-y_mid)*s with s={scale:.8g} units/um. "
                      "Uniform, isotropic, centred — no rotation or stretch."),
        "level_order": "outermost (lowest pressure) first",
        "simplify_tolerance_scene_units": IFP_SIMPLIFY_TOLERANCE,
        "levels": levels,
    }
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    return doc


# --------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=os.path.join(_ROOT, OUT_REL))
    ap.add_argument("--ifp-out", default=os.path.join(_ROOT, IFP_OUT_REL))
    ap.add_argument("--no-glb", action="store_true",
                    help="skip the point cloud and write only the solved-IFP "
                         "contours (both are built from the same reconstruction "
                         "and the same transform either way)")
    ap.add_argument("--no-ifp", action="store_true",
                    help="skip the solved-IFP contours and write only the point "
                         "cloud")
    ap.add_argument("--sample", default=SAMPLE)
    ap.add_argument("--dx-um", type=float, default=DX_UM)
    ap.add_argument("--max-cells", type=int, default=0,
                    help="0 = keep every measured cell (default; the full section "
                         "already fits the web budget). If set, cells are chosen by "
                         "UNIFORM RANDOM subsample without replacement (seed 0), "
                         "which thins every region by the same factor and so "
                         "preserves relative spatial structure. Density-aware or "
                         "region-of-interest sampling is deliberately not offered: "
                         "it would misrepresent the section.")
    args = ap.parse_args()

    os.chdir(_ROOT)  # the loaders resolve data_cache/ relative to the repo root
    print("reconstructing REAL Xenium GBM microenvironment ...", flush=True)
    micro = reconstruct_real(source="xenium", sample=args.sample, dx_um=args.dx_um,
                             oxygen_mode=OXYGEN_MODE, state_mode=STATE_MODE)
    rep = reconstruction_report(micro)
    cells = micro.cells
    n_total = cells.n
    print(f"  {n_total} measured cells, patient_id={micro.patient_id}")

    # -- subsample (default: none) --------------------------------------------
    if args.max_cells and n_total > args.max_cells:
        rng = np.random.default_rng(0)
        idx = np.sort(rng.choice(n_total, args.max_cells, replace=False))
        subsample = f"uniform random without replacement, seed 0, {len(idx)}/{n_total}"
    else:
        idx = np.arange(n_total)
        subsample = f"none — all {n_total} measured cells retained"
    x_um, y_um = cells.x[idx], cells.y[idx]
    ctype = cells.cell_type[idx].astype(int)

    # -- coordinate transform (see module docstring) --------------------------
    x_lo, x_hi = float(x_um.min()), float(x_um.max())
    y_lo, y_hi = float(y_um.min()), float(y_um.max())
    w_um, h_um = x_hi - x_lo, y_hi - y_lo
    scale = FIT_EXTENT / max(w_um, h_um)
    x_mid, y_mid = 0.5 * (x_lo + x_hi), 0.5 * (y_lo + y_hi)
    X = (x_um - x_mid) * scale
    Z = (y_um - y_mid) * scale
    Y = np.zeros_like(X)                       # planar section: no measured depth
    P = np.stack([X, Y, Z], axis=1).astype(np.float32)
    um_per_unit = 1.0 / scale

    # -- colour ---------------------------------------------------------------
    pal = display_palette()
    legend = {}
    groups = []
    for ct in sorted(set(ctype.tolist())):
        m = ctype == ct
        n = int(m.sum())
        name = CellType(ct).name
        hex_disp, lin = pal[ct]
        rgba = np.empty((n, 4), dtype=np.uint8)
        rgba[:, :3] = np.round(np.clip(np.array(lin), 0, 1) * 255).astype(np.uint8)
        rgba[:, 3] = 255
        legend[name] = {"hex": hex_disp, "n_cells": n,
                        "fraction": round(n / len(idx), 4)}
        groups.append((name, P[m], rgba, {
            "cell_type": name, "cell_type_id": int(ct), "n_cells": n,
            "display_hex_srgb": hex_disp,
            "canonical_hex_srgb": CELLTYPE_COLORS[int(ct)],
        }))
        print(f"  {name:<12s} {n:>7d}  {hex_disp}")

    extras = {
        "what_this_is": (
            "Every point is one segmented cell from a real human glioblastoma "
            "section at its measured micrometre centroid. No synthetic geometry."
        ),
        "dataset": "10x Genomics Xenium Human Brain — Glioblastoma (v1.3.0)",
        "citation": CITATION,
        "licence": "CC BY 4.0",
        "sample_id": "Xenium_V1_FFPE_Human_Brain_Glioblastoma_With_Addon",
        "pipeline": ("tmesim.reconstruct.real_pipeline.reconstruct_real"
                     f"(source='xenium', sample='{args.sample}', dx_um={args.dx_um}, "
                     f"oxygen_mode='{OXYGEN_MODE}', state_mode='{STATE_MODE}')"),
        "patient_id": micro.patient_id,
        "n_cells_measured": int(n_total),
        "n_points": int(len(idx)),
        "subsample": subsample,
        "section_extent_um": [round(w_um, 1), round(h_um, 1)],
        "um_per_scene_unit": round(um_per_unit, 2),
        "transform": ("X=(x_um-x_mid)*s, Y=0, Z=(y_um-y_mid)*s with "
                      f"s={scale:.8g} units/um (uniform, isotropic, centred). "
                      "Y is exactly 0 for every point: Xenium measures a single "
                      "2D section, so the vertical axis carries no data and no "
                      "depth was invented."),
        "colour_means": ("cell type, from threshold-gated marker typing on the "
                         "319-gene Xenium panel (matched_null_z aggregation, "
                         "z-threshold 1.0 vs a detection-rate-matched null)"),
        "marker_map": micro.metadata.get("marker_map"),
        "typing_note": micro.metadata.get("typing_note"),
        "legend": legend,
        "cell_state": (f"all cells are {CellState.VIABLE.name} — this is the "
                       "untreated reconstruction; no colour encodes a drug "
                       "response, a kill, or an outcome"),
        "has_outcome_label": False,
        "colour_note": (
            "Hue is tmesim.viz.explainer.CELLTYPE_COLORS verbatim (same as "
            "every figure in results/); HSL lightness is set to "
            f"{DISPLAY_LIGHTNESS} and saturation capped at "
            f"{DISPLAY_MAX_SATURATION}, for equal weight on the dark "
            "background. COLOR_0 is stored linear-light per the glTF spec."
        ),
        "reconstruction_report": {k: (round(v, 6) if isinstance(v, float) else v)
                                  for k, v in rep.items()
                                  if isinstance(v, (int, float, str, bool))},
    }

    if not args.no_glb:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        total = write_glb(args.out, groups, extras)
        print(f"\nwrote {args.out}  ({total / 1e6:.2f} MB)")
        print(f"  points      : {len(idx)} ({subsample})")
        print(f"  section     : {w_um:.0f} x {h_um:.0f} um -> "
              f"{w_um * scale:.3f} x {h_um * scale:.3f} scene units, y=0")
        print(f"  scale       : {um_per_unit:.1f} um per scene unit")
        print(f"  colour      : cell type ({', '.join(legend)})")
    else:
        print("\nskipping tumor.glb (--no-glb)")

    # -- solved interstitial fluid pressure, same transform -------------------
    if not args.no_ifp:
        print("\nsolving interstitial fluid pressure ...", flush=True)
        doc = export_ifp_contours(args.ifp_out, micro, x_mid, y_mid, scale,
                                  args.sample, args.dx_um)
        size = os.path.getsize(args.ifp_out)
        print(f"wrote {args.ifp_out}  ({size / 1e3:.1f} kB)")
        print(f"  raw solved p: 0 .. {doc['raw_field_max_mmhg']:.3f} mmHg "
              f"(mean in tissue {doc['raw_field_mean_in_tissue_mmhg']:.4f})")
        print(f"  displayed   : {IFP_COARSE_GRAIN_UM:g} um local mean, "
              f"max {doc['displayed_field_max_mmhg']:.4f} mmHg")
        for lv in doc["levels"]:
            print(f"  level {lv['fraction_of_max']:.2f} = "
                  f"{lv['value_mmhg']:.4f} mmHg  {lv['color_srgb_hex']}  "
                  f"{lv['n_paths']} path(s), {lv['n_vertices']} vertices")


if __name__ == "__main__":
    main()
