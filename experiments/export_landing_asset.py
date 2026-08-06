"""Export the REAL reconstructed GBM tissue section to a glTF-Binary point cloud
for the landing page's hero canvas (``landing-3d/public/models/tumor.glb``).

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

OUTPUT
------
``landing-3d/public/models/tumor.glb`` — a valid binary glTF 2.0 containing one
root node with one child node per cell type, each a POINTS-mode primitive with
POSITION (float32 VEC3) and COLOR_0 (normalised uint8 VEC4). Written by hand
with numpy + struct + json; no glTF dependency is added to the project.
Provenance is embedded in ``asset.extras`` so the file is self-describing.

USAGE
-----
    .venv/bin/python experiments/export_landing_asset.py
    .venv/bin/python experiments/export_landing_asset.py --max-cells 20000
"""
import warnings; warnings.filterwarnings("ignore")

import argparse
import colorsys
import json
import os
import struct
import sys

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:  # allow `python experiments/export_landing_asset.py`
    sys.path.insert(0, _ROOT)

from tmesim.core import CellState, CellType                      # noqa: E402
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


# --------------------------------------------------------------------- main --
def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=os.path.join(_ROOT, OUT_REL))
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

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    total = write_glb(args.out, groups, extras)
    print(f"\nwrote {args.out}  ({total / 1e6:.2f} MB)")
    print(f"  points      : {len(idx)} ({subsample})")
    print(f"  section     : {w_um:.0f} x {h_um:.0f} um -> "
          f"{w_um * scale:.3f} x {h_um * scale:.3f} scene units, y=0")
    print(f"  scale       : {um_per_unit:.1f} um per scene unit")
    print(f"  colour      : cell type ({', '.join(legend)})")


if __name__ == "__main__":
    main()
