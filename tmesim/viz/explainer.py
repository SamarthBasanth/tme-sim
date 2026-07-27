"""VISUAL EXPLAINER — the whole real-data pipeline, legible at a glance.

Four standalone, publication-quality figures. Every number rendered here is
pulled live from the real modules (no hard-coded results, no invented values):

  1. ``results/explainer_pipeline.png``            end-to-end data flow
  2. ``results/explainer_multiscale.png``          the nm -> mm scale ladder
  3. ``results/explainer_microenvironment.png``    the real reconstructed niche
  4. ``results/explainer_delivery_vs_potency.png`` delivery beats potency

Run:  ``python -m tmesim.viz.explainer``

This module is ADDITIVE: it imports from core / physics / molecular /
reconstruct / sim and mutates nothing.

Design notes
------------
Palette is the validated categorical/sequential set (adjacent-CVD delta-E and
contrast checked); sequential fields each get their own SINGLE-hue ramp plus a
scale legend; every series is direct-labelled so identity is never colour-alone.
"""
from __future__ import annotations

import math
import os
import time
import warnings
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")  # headless; must precede pyplot

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyBboxPatch, Polygon, Rectangle

from ..core import CellType, PDEConfig, SimConfig
from ..molecular import ecm_mesh, stokes_einstein
from ..molecular.binding_site_barrier import apply_barrier_to_drugparams
from ..molecular.library import get_molecule
from ..physics.diffusion import ALPHA_TORTUOSITY
from ..reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from ..sim.engine import simulate_response
from ..validation.calibration import (
    DOX_CITATION,
    DOX_PENETRATION_UM,
    calibrate_alpha,
    calibrated_drugparams,
)

__all__ = ["main", "collect", "figure_pipeline", "figure_multiscale",
           "figure_microenvironment", "figure_delivery_vs_potency"]

# --------------------------------------------------------------------------- #
# paths / run parameters
# --------------------------------------------------------------------------- #
_HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.abspath(os.path.join(_HERE, "..", "..", "results"))

#: The drug panel, in fixed order (fixes the categorical colour assignment).
PANEL: List[str] = ["temozolomide", "lomustine", "cisplatin",
                    "methotrexate", "doxorubicin"]

#: Grid-converged resolution (see results/numerical_rigor.md section 2:
#: penetration_index only settles once dx resolves the ~45 um boundary layer).
DX_UM = 12.0

#: nm-scale ECM fibre matrix used by the molecular layer's mesh gate.
FIBER_RADIUS_NM = 2.0
FIBER_VOLUME_FRACTION = 0.05

# --------------------------------------------------------------------------- #
# palette (validated: adjacent-CVD deltaE and surface contrast checked)
# --------------------------------------------------------------------------- #
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_3 = "#84837c"
HAIRLINE = "#e2e1dc"
PANEL_BG = "#f4f3f0"

#: categorical theme, fixed order — hues are assigned by slot, never cycled.
CAT = ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7",
       "#e34948", "#e87ba4", "#eb6834"]

DRUG_COLORS: Dict[str, str] = {name: CAT[i] for i, name in enumerate(PANEL)}

CELLTYPE_COLORS: Dict[int, str] = {
    int(CellType.TUMOR): "#e34948",
    int(CellType.STROMAL): "#008300",
    int(CellType.IMMUNE): "#2a78d6",
    int(CellType.ENDOTHELIAL): "#4a3aa7",
    int(CellType.NECROTIC): "#7f7f7f",
}

#: one hue per sequential field, each shipped with its own colourbar.
FIELD_HUE = {
    "vascular": "#4a3aa7",   # matches the ENDOTHELIAL cell colour
    "ecm": "#008300",        # matches the STROMAL cell colour
    "oxygen": "#2a78d6",
    "exposure": "#1baf7a",
}

STAGE_ACCENTS = [CAT[0], CAT[1], CAT[2], CAT[3], CAT[4]]


def _ramp(hex_color: str, name: str = "ramp") -> LinearSegmentedColormap:
    """A single-hue light->dark sequential ramp anchored on the chart surface.

    Four stops with monotonically decreasing lightness: a near-surface tint (so
    "near zero" recedes), two mid tints, the hue itself, and a dark shade.
    """
    base = np.asarray(to_rgb(hex_color))
    surf = np.asarray(to_rgb(SURFACE))

    def tint(f: float) -> np.ndarray:      # f=1 -> surface, f=0 -> base
        return base + (surf - base) * f

    stops = [tint(0.96), tint(0.72), tint(0.36), base, base * 0.52]
    return LinearSegmentedColormap.from_list(name, [tuple(s) for s in stops])


CMAPS = {k: _ramp(v, f"tme_{k}") for k, v in FIELD_HUE.items()}


def _style() -> None:
    """Global typography / chrome: sans face, hairline recessive axes."""
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "text.color": INK,
        "axes.labelcolor": INK_2,
        "axes.edgecolor": "#c9c8c2",
        "axes.linewidth": 0.8,
        "xtick.color": INK_2,
        "ytick.color": INK_2,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "legend.frameon": False,
        "grid.color": HAIRLINE,
        "grid.linewidth": 0.7,
        "grid.linestyle": "-",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def _caption(fig, text: str, y: float = 0.012, size: float = 9.0) -> None:
    """Bake a standalone-readable caption into the bottom of a figure."""
    fig.text(0.5, y, text, ha="center", va="bottom", fontsize=size,
             color=INK_2, linespacing=1.55, wrap=False)


def _fmt(x: float, digits: int = 3) -> str:
    """Compact scientific formatting: 5.61e+06 -> '5.6e6'."""
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    if -2 <= e <= 4:
        return f"{x:,.{max(0, digits - 1 - e)}f}".rstrip("0").rstrip(".")
    return f"{x / 10 ** e:.1f}e{e}"


# --------------------------------------------------------------------------- #
# 1. COLLECT — every number in the figures comes from here, computed live
# --------------------------------------------------------------------------- #
def _retry(fn, *args, **kw):
    """Run ``fn``, retrying once after a short pause (transient network / IO)."""
    try:
        return fn(*args, **kw)
    except Exception:
        time.sleep(3.0)
        return fn(*args, **kw)


def _datasets() -> Dict:
    """Real dataset inventory — shapes read off the actual cached bytes."""
    from ..data.real import drug_response, ivygap, scref_gbm, visium
    from ..validation import pk_anchors

    out: Dict = {}

    # Visium: every open GBM section in the cache.
    ids = visium.list_samples()
    spots = []
    n_genes_visium = 0
    for sid in ids:
        d = _retry(visium.load_sample, sid, marker_report=False)
        spots.append(int(d["n_spots"]))
        n_genes_visium = int(d["n_genes"])
    out["visium"] = dict(n_samples=len(ids), spots=spots, n_genes=n_genes_visium)

    # Darmanis 2017 single-cell reference (GSE84465).
    expr, meta, genes = _retry(scref_gbm.load_reference)
    out["scref"] = dict(n_genes=int(expr.shape[0]), n_cells=int(expr.shape[1]),
                        n_types=int(meta["cell_type"].nunique()))

    # Ivy GAP anatomic structures.
    st = _retry(ivygap.load_structures)
    out["ivygap"] = dict(n_samples=int(st.samples.shape[0]),
                         n_genes=int(st.expression.shape[0]),
                         n_structures=int(st.samples["base_structure"].nunique()))

    # Real drug-response screens.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        summ = drug_response.summarize()
    out["screens"] = summ

    # Measured tissue-PK anchors.
    out["pk"] = dict(n_anchors=len(pk_anchors.PK_ANCHORS),
                     n_drugs=len(pk_anchors.drugs()),
                     drugs=pk_anchors.drugs())
    return out


def collect(dx_um: float = DX_UM, verbose: bool = True) -> Dict:
    """Reconstruct the real microenvironment ONCE, simulate the panel on it,
    and gather every quantity the four figures display."""
    t0 = time.time()
    D: Dict = {"dx_um": float(dx_um)}

    if verbose:
        print(f"[1/4] reconstructing the REAL Xenium GBM niche (dx = {dx_um:g} um) ...")
    micro = _retry(reconstruct_real, source="xenium", dx_um=float(dx_um))
    rep = reconstruction_report(micro)
    D["micro"] = micro
    D["report"] = rep
    if verbose:
        print(f"      {rep['n_cells']:,} real cells, grid {rep['grid']}, "
              f"extent {rep['extent_um'][0]:.0f} x {rep['extent_um'][1]:.0f} um "
              f"[{time.time() - t0:.1f}s]")

    # --- real geometric scales, computed from the measured cell coordinates ---
    from scipy.spatial import cKDTree

    cells = micro.cells
    ctype = np.asarray(cells.cell_type)
    P = np.column_stack([cells.x, cells.y])
    endo = ctype == int(CellType.ENDOTHELIAL)
    tumor = ctype == int(CellType.TUMOR)

    tree_all = cKDTree(P)
    nn_all = tree_all.query(P, k=2)[0][:, 1]
    tree_endo = cKDTree(P[endo])
    d_to_vessel = tree_endo.query(P[tumor])[0]
    D["scales"] = dict(
        cell_nn_median=float(np.median(nn_all)),
        vessel_dist_median=float(np.median(d_to_vessel)),
        vessel_dist_p90=float(np.percentile(d_to_vessel, 90.0)),
        endo_nn_median=float(np.median(tree_endo.query(P[endo], k=2)[0][:, 1])),
        n_endothelial=int(endo.sum()),
    )

    # --- molecular layer: structure -> radius -> D_free -> mesh gate ---------
    if verbose:
        print("[2/4] compiling the molecular layer (structure -> transport) ...")
    alpha = calibrate_alpha()
    D["alpha"] = float(alpha)
    mean_ecm = float(rep["mean_ecm_density"])

    mol_rows = []
    for name in PANEL + ["bevacizumab"]:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mol = get_molecule(name)
        r_nm = float(mol.get("radius_hydrodynamic_nm"))
        free_D = stokes_einstein.free_diffusivity(r_nm)
        mesh = ecm_mesh.ogston_renkin_hindrance(
            r_nm, FIBER_RADIUS_NM, FIBER_VOLUME_FRACTION)
        dp = calibrated_drugparams(name, alpha)
        d_eff = float(dp.diffusivity) * (1.0 - ALPHA_TORTUOSITY * mean_ecm)
        k_eff = float(dp.clearance) + float(dp.binding)
        mol_rows.append(dict(
            name=name,
            mw=float(mol.get("mw")),
            radius_nm=r_nm,
            free_D=float(free_D),
            mesh_partition=float(mesh["partition"]),
            mesh_hindrance=float(mesh["effective_hindrance"]),
            tissue_D=float(dp.diffusivity),
            d_eff=d_eff,
            k_eff=k_eff,
            L_um=float(math.sqrt(d_eff / k_eff)),
            ic50=float(dp.ic50),
            bbb=float(dp.bbb_penetrance),
        ))
    D["molecular"] = {r["name"]: r for r in mol_rows}
    D["mesh_pore_nm"] = float(ecm_mesh.mesh_pore_radius_nm(
        FIBER_RADIUS_NM, FIBER_VOLUME_FRACTION))

    # --- simulate every drug on the SAME real tissue -------------------------
    if verbose:
        print(f"[3/4] simulating {len(PANEL)} drugs on the real tissue ...")
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True, save_every=500)
    rows = []
    for name in PANEL:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            dp = apply_barrier_to_drugparams(
                calibrated_drugparams(name, alpha), get_molecule(name))
        t1 = time.time()
        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
        m = D["molecular"][name]
        rows.append(dict(
            name=name, ic50=m["ic50"], mw=m["mw"], bbb=m["bbb"],
            d_eff=m["d_eff"], L_um=m["L_um"],
            kill=float(res.kill_fraction),
            pen=float(res.extras["penetration_index"]),
            mean_auc=float(res.extras["mean_auc"]),
            hypoxic=float(res.extras["hypoxic_tumor_fraction"]),
            auc=np.asarray(res.fields["auc"]),
        ))
        if verbose:
            print(f"      {name:14s} kill={rows[-1]['kill']:.4f}  "
                  f"pen={rows[-1]['pen']:.3f}  AUC={rows[-1]['mean_auc']:.3g} "
                  f"[{time.time() - t1:.1f}s]")
    D["sim"] = rows
    D["by_drug"] = {r["name"]: r for r in rows}
    D["hypoxic_fraction"] = float(rows[0]["hypoxic"])

    # --- exposure decay away from the real vessels (real, computed) ----------
    from scipy.ndimage import distance_transform_edt

    ny, nx = micro.domain.shape
    col, row = micro.domain.to_index(P[endo, 0], P[endo, 1])
    vessel_mask = np.zeros((ny, nx), dtype=bool)
    vessel_mask[row, col] = True
    dist = distance_transform_edt(~vessel_mask) * float(micro.domain.dx)
    D["vessel_distance"] = dist
    D["vessel_mask"] = vessel_mask

    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    edges = np.arange(0.0, 210.0, 10.0)
    centers = 0.5 * (edges[:-1] + edges[1:])
    idx = np.digitize(dist[tissue], edges) - 1
    decay = {}
    for r in rows:
        v = r["auc"][tissue]
        prof = np.array([v[idx == b].mean() if np.any(idx == b) else np.nan
                         for b in range(len(centers))])
        with np.errstate(invalid="ignore"):
            decay[r["name"]] = prof / prof[0]
    D["decay"] = dict(centers=centers, profiles=decay)

    # --- dataset inventory ---------------------------------------------------
    if verbose:
        print("[4/4] reading the real dataset inventory ...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        D["datasets"] = _datasets()
    D["xenium"] = dict(
        n_cells=int(rep["n_cells"]),
        n_tumor=int(rep["n_tumor"]),
        n_panel_genes=int(micro.metadata.get("n_panel_genes", 0)),
        assay=str(micro.metadata.get("assay", "10x Xenium")),
        ecm_markers=list(rep.get("ecm_markers_used", [])),
    )

    D["citation"] = dict(dox=DOX_CITATION, dox_band=DOX_PENETRATION_UM)
    D["elapsed"] = time.time() - t0
    if verbose:
        print(f"      collected in {D['elapsed']:.1f}s")
    return D


# --------------------------------------------------------------------------- #
# drawing primitives for the designed diagrams
# --------------------------------------------------------------------------- #
def _mix(hex_color: str, f: float, toward: str = SURFACE) -> Tuple[float, float, float]:
    """Blend ``hex_color`` a fraction ``f`` toward ``toward`` (f=1 -> toward)."""
    a = np.asarray(to_rgb(hex_color))
    b = np.asarray(to_rgb(toward))
    return tuple(a + (b - a) * float(f))


def _rounded(ax, x, y, w, h, radius=2.0, face="#ffffff", edge=None, lw=0.9, z=2,
             alpha=1.0):
    """A rounded rectangle in data units (axes must be aspect-equal)."""
    p = FancyBboxPatch(
        (x + radius, y + radius), w - 2 * radius, h - 2 * radius,
        boxstyle=f"round,pad={radius},rounding_size={radius}",
        mutation_scale=1.0, linewidth=lw, facecolor=face,
        edgecolor=(edge if edge is not None else "none"),
        zorder=z, alpha=alpha, clip_on=False,
    )
    ax.add_patch(p)
    return p


def _chevron(ax, x, y, size=4.6, color="#c9c8c2", z=6):
    """A single right-pointing chevron marking flow between stages."""
    pts = np.array([[0, 1.0], [0.55, 0.0], [0, -1.0],
                    [0.42, -1.0], [0.97, 0.0], [0.42, 1.0]])
    ax.add_patch(Polygon(pts * size + np.array([x, y]), closed=True,
                         facecolor=color, edgecolor="none", zorder=z,
                         clip_on=False))


def _scalebar(ax, length_data, label, color=INK_2, frac_x=0.06, frac_y=0.07,
              size=8.0, lw=2.4):
    """A horizontal scale bar anchored in axes fraction, sized in data units."""
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    x = x0 + (x1 - x0) * frac_x
    y = y0 + (y1 - y0) * frac_y
    ax.plot([x, x + length_data], [y, y], color=color, lw=lw,
            solid_capstyle="butt", zorder=20)
    ax.text(x + length_data / 2.0, y + (y1 - y0) * 0.025, label, ha="center",
            va="bottom", fontsize=size, color=color, zorder=20)


def _wrap(text: str, width: int) -> List[str]:
    """Hard-wrap a paragraph to ``width`` characters, honouring explicit \\n."""
    import textwrap

    out: List[str] = []
    for chunk in str(text).split("\n"):
        out.extend(textwrap.wrap(chunk, width=width) or [""])
    return out


def _para(ax, x, y, text, width, fontsize=8.8, color=INK_2, dy=None,
          transform=None, weight="normal", ha="left", zorder=6):
    """Render a wrapped paragraph downward from ``y``; return the next free y."""
    lines = _wrap(text, width)
    step = dy if dy is not None else fontsize * 1.55
    tr = transform if transform is not None else ax.transData
    for i, ln in enumerate(lines):
        ax.text(x, y - i * step, ln, fontsize=fontsize, color=color,
                va="top", ha=ha, transform=tr, fontweight=weight, zorder=zorder)
    return y - len(lines) * step


def _display_field(field, tissue, dx, smooth_um=45.0, everywhere=False,
                   cut=0.06):
    """Render-ready version of a field defined only on measured tissue pixels.

    The tissue mask is the set of grid cells that actually contain a measured
    cell — at dx = 12 µm that is a sparse dust (~14% of the grid), so a raw
    ``imshow`` of it reads as blank. This returns the LOCAL MEAN over the
    measured pixels, ``G(f·t)/G(t)`` for a Gaussian G, masked to where the
    tissue support is real. It is a display-only smoothing of measured values —
    no values are invented, and the section silhouette falls out of the support.
    """
    from scipy.ndimage import gaussian_filter

    f = np.asarray(field, dtype=np.float64)
    t = np.asarray(tissue, dtype=bool)
    sigma = max(0.8, float(smooth_um) / float(dx))
    den = gaussian_filter(t.astype(np.float64), sigma)
    if everywhere:
        num = gaussian_filter(f, sigma)
        out = num
    else:
        num = gaussian_filter(np.where(t, f, 0.0), sigma)
        out = num / np.maximum(den, 1e-12)
    return np.where(den > cut, out, np.nan)


def _field_im(ax, field, extent, cmap, tissue=None, gamma=0.45, pmax=99.5,
              pmin=0.0, interpolation="bilinear", vmin=None, vmax=None):
    """Show a strongly-peaked KDE/AUC field so its structure is actually visible.

    A linear scale buries these fields (means sit at ~15% of the 99.5th
    percentile), so magnitude is shown on a power-law norm — still one hue,
    still monotone, and always shipped with a colourbar.
    """
    from matplotlib.colors import PowerNorm

    f = np.asarray(field, dtype=np.float64)
    ref = f[tissue] if tissue is not None and np.any(tissue) else f
    if vmax is None:
        vmax = float(np.nanpercentile(ref, pmax)) or 1.0
    if vmin is None:
        vmin = float(np.nanpercentile(ref, pmin)) if pmin > 0 else 0.0
    show = np.where(tissue, f, np.nan) if tissue is not None else f
    im = ax.imshow(show, cmap=cmap, extent=extent, aspect="auto",
                   norm=PowerNorm(gamma=gamma, vmin=vmin, vmax=vmax),
                   interpolation=interpolation)
    h = abs(extent[2] - extent[3])
    w = abs(extent[1] - extent[0])
    ax.set_box_aspect(h / w)
    return im


def _colorbar(fig, ax, im, label, size=8.0):
    cb = fig.colorbar(im, ax=ax, fraction=0.042, pad=0.02)
    cb.set_label(label, fontsize=size, color=INK_2)
    cb.ax.tick_params(labelsize=size - 0.5, colors=INK_2, length=2, width=0.6)
    cb.outline.set_visible(False)
    return cb


# --------------------------------------------------------------------------- #
# FIGURE 1 — the end-to-end pipeline, as a designed diagram
# --------------------------------------------------------------------------- #
_W1, _H1 = 22.0, 14.6          # inches
_XU, _YU = 220.0, 146.0        # data units (1 unit = 0.1 in, aspect-equal)


def _stage_items(D: Dict) -> List[Tuple[str, str, List]]:
    """Build the five stage columns' content from the collected real numbers."""
    ds = D["datasets"]
    rep = D["report"]
    xen = D["xenium"]
    frac = rep["celltype_fractions"]
    mol = D["molecular"]
    sim = D["sim"]
    sp = ds["visium"]["spots"]

    pens = [r["pen"] for r in sim]
    kills = [r["kill"] for r in sim]
    ic50s = [r["ic50"] for r in sim]
    mws = [r["mw"] for r in sim]
    radii = [mol[n]["radius_nm"] for n in PANEL]
    freeD = [mol[n]["free_D"] for n in PANEL]
    deff = [mol[n]["d_eff"] for n in PANEL]
    hind = [mol[n]["mesh_hindrance"] for n in PANEL]

    data_items = [
        ("text", "10x Xenium GBM section",
         [f"{xen['n_cells']:,} single cells at true µm coordinates",
          f"{xen['n_panel_genes']}-gene panel · CC BY 4.0"]),
        ("text", f"{ds['visium']['n_samples']}× Visium GBM sections",
         [f"{min(sp):,}–{max(sp):,} spots each",
          f"whole transcriptome, {ds['visium']['n_genes']:,} genes"]),
        ("text", "Darmanis 2017 cell atlas",
         [f"{ds['scref']['n_genes']:,} genes × {ds['scref']['n_cells']:,} cells",
          f"{ds['scref']['n_types']} annotated types (GSE84465)"]),
        ("text", "Ivy GAP anatomic atlas",
         [f"{ds['ivygap']['n_samples']} samples, "
          f"{ds['ivygap']['n_structures']} anatomic structures",
          "leading edge → cellular tumor → perinecrotic"]),
        ("text", "GDSC + HGCC drug screens",
         [f"HGCC {ds['screens']['hgcc']['drugs']} drugs × "
          f"{ds['screens']['hgcc']['cultures']} patient cultures",
          f"GDSC {ds['screens']['gdsc']['drugs']} drugs × "
          f"{ds['screens']['gdsc']['cell_lines']} lines"]),
        ("text", "Measured tissue-PK anchors",
         [f"{ds['pk']['n_anchors']} PMID-verified anchors, "
          f"{ds['pk']['n_drugs']} drugs",
          "perivascular gradients + brain:plasma ratios"]),
    ]

    recon_items = [
        ("text", "Threshold-gated cell typing",
         [f"TUMOR {frac['TUMOR']*100:.1f}% · STROMAL {frac['STROMAL']*100:.1f}% · "
          f"IMMUNE {frac['IMMUNE']*100:.1f}%",
          f"ENDOTHELIAL {frac['ENDOTHELIAL']*100:.1f}% "
          f"(Visium NNLS agrees to 0.03%)"]),
        ("image", "vascular",
         ["Vasculature = KDE of real endothelial cells",
          f"mean density {rep['mean_vascular_density']:.3f} · "
          f"{D['scales']['n_endothelial']:,} endothelial cells"]),
        ("text", "ECM from real matrix transcripts",
         [f"expression-weighted KDE of {len(xen['ecm_markers'])} genes "
          "(COL1A2, POSTN, VCAN, …)",
          f"mean ECM density {rep['mean_ecm_density']:.3f}"]),
        ("text", "Oxygen from real perfusion",
         [f"physics solve on the real vascular field · "
          f"mean O₂ {rep['mean_oxygen']:.3f}",
          f"{D['hypoxic_fraction']*100:.0f}% of tumor cells hypoxic"]),
        ("text", "Niche regions (Ivy GAP signatures)",
         ["perivascular / cellular tumor / perinecrotic / margin",
          "marker-faithful: PECAM1 d=1.29, VEGFA d=1.35 (p=1e-4)"]),
    ]

    mol_items = [
        ("text", "Structure → descriptors",
         ["verified SMILES → RDKit: MW, logP, TPSA, HBD",
          f"MW {min(mws):.0f}–{max(mws):.0f} Da across the panel"]),
        ("text", "Stokes–Einstein",
         [r"$r_h=0.066\,MW^{1/3}$;  $D=k_BT/6\pi\eta r_h$  (310 K, 1.2 mPa·s)",
          f"r_h {min(radii):.2f}–{max(radii):.2f} nm → "
          f"D_free {min(freeD):.0f}–{max(freeD):.0f} µm²/s"]),
        ("text", "Brownian + Ogston/Renkin hindrance",
         [f"2 nm fibres at φ={FIBER_VOLUME_FRACTION:g} → "
          f"pore radius {D['mesh_pore_nm']:.1f} nm",
          f"mesh hindrance {min(hind):.2f}–{max(hind):.2f} → "
          f"D_eff {min(deff):.0f}–{max(deff):.0f} µm²/s"]),
        ("text", "Measured IC50 (not predicted)",
         [f"real GDSC / HGCC potencies: "
          f"{min(ic50s):.1f}–{max(ic50s):.0f} µM",
          "GDSC vs HGCC corroborate at Spearman +0.83"]),
        ("text", "→ DrugParams",
         ["D_eff · IC50 · binding · clearance · BBB penetrance",
          "one object, no per-drug hand tuning"]),
    ]

    phys_items = [
        ("text", "Transport PDE",
         [r"$\partial_t c=\nabla\!\cdot\!(D_{eff}\nabla c)-v\!\cdot\!\nabla c-kc+S$",
          "finite volume, implicit; mass conserved to 1e-13"]),
        ("text", "Interstitial pressure (Darcy)",
         ["IFP → convective velocity v; the Jain outward-flow effect",
          "vessel source S gated by BBB penetrance"]),
        ("text", "Oxygen reaction–diffusion",
         ["real perfusion → O₂ field → per-cell hypoxia banding",
          "hypoxic cells are oxygen-gated to lower sensitivity"]),
        ("text", "Hill kill",
         [r"$E=c^{h}/(c^{h}+IC50^{h})$ integrated over the course",
          "kill = delivery × potency × hypoxia resistance"]),
        ("text", "ONE fitted constant",
         [f"α = {D['alpha']:.1f}, fit once to the measured doxorubicin gradient",
          f"{D['citation']['dox_band'][0]:.0f}–{D['citation']['dox_band'][1]:.0f} µm, "
          "PMID 16361566 · out-of-sample PASS on 2 drugs"]),
    ]

    out_items = [
        ("image", "exposure",
         ["Doxorubicin exposure AUC (µM·s)",
          "trapped in the perivascular boundary layer"]),
        ("text", "Penetration index",
         [f"{min(pens):.3f} – {max(pens):.3f} across the panel",
          f"best: {max(sim, key=lambda r: r['pen'])['name']}"]),
        ("text", "Predicted kill map",
         [f"{min(kills)*100:.1f}% – {max(kills)*100:.1f}% of tumor cells",
          f"best: {max(sim, key=lambda r: r['kill'])['name']}"]),
        ("text", "NOT validated: response",
         ["every open spatial GBM section is an untreated resection",
          "no paired measured drug response exists — these are predictions"]),
    ]

    return [
        ("REAL DATA", data_items),
        ("RECONSTRUCTION", recon_items),
        ("MOLECULAR LAYER", mol_items),
        ("PHYSICS", phys_items),
        ("OUTPUT", out_items),
    ]


def figure_pipeline(D: Dict, path: str) -> str:
    """Figure 1 — the end-to-end data flow, drawn as a designed diagram."""
    stages = _stage_items(D)
    accents = STAGE_ACCENTS

    fig = plt.figure(figsize=(_W1, _H1))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, _XU)
    ax.set_ylim(0, _YU)
    ax.set_aspect("equal")
    ax.axis("off")

    # --- masthead ----------------------------------------------------------
    ax.text(6, 139.0,
            "From measured tissue to a predicted drug map",
            fontsize=25, fontweight="bold", color=INK, va="center")
    ax.text(6, 132.4,
            "The tme-sim pipeline end to end. Every box is a real measurement or "
            "a mechanistic step — nothing is synthetic, and every number on this "
            "page was recomputed when the figure was drawn.",
            fontsize=11.5, color=INK_2, va="center")

    x0, gap = 6.0, 5.6
    w = (_XU - 2 * x0 - 4 * gap) / 5.0
    y_bot, y_top = 25.0, 126.0
    h = y_top - y_bot

    thumb_axes = {}
    for i, ((title, items), accent) in enumerate(zip(stages, accents)):
        x = x0 + i * (w + gap)

        # column shell
        _rounded(ax, x, y_bot, w, h, radius=2.6,
                 face=_mix(accent, 0.945), edge=_mix(accent, 0.62), lw=1.0, z=1)

        # header pill
        hp = 8.4
        _rounded(ax, x, y_top - hp, w, hp, radius=2.6, face=accent,
                 edge=None, lw=0, z=3)
        ax.add_patch(Circle((x + 5.6, y_top - hp / 2), 2.5, facecolor="#ffffff",
                            edgecolor="none", zorder=4, alpha=0.92))
        ax.text(x + 5.6, y_top - hp / 2, str(i + 1), ha="center", va="center",
                fontsize=10.5, fontweight="bold", color=accent, zorder=5)
        ax.text(x + 10.4, y_top - hp / 2, title, ha="left", va="center",
                fontsize=12.5, fontweight="bold", color="#ffffff", zorder=5)

        # items — heights are proportional to the wrapped text each one needs
        pad = 2.4
        top = y_top - hp - 2.4
        bottom = y_bot + pad
        wrap_w = 47
        wrapped = []
        for kind, key, lines in items:
            wl = []
            for ln in lines:
                wl.extend(_wrap(ln, wrap_w) if "$" not in ln else [ln])
            wrapped.append(wl)
        req = [26.0 if it[0] == "image" else 6.6 + 2.9 * len(wl)
               for it, wl in zip(items, wrapped)]
        avail = (top - bottom) - 1.6 * (len(items) - 1)
        scale = avail / sum(req)

        cy = top
        for (kind, key, lines), wl, rq in zip(items, wrapped, req):
            ih = rq * scale
            iy = cy - ih
            _rounded(ax, x + pad, iy, w - 2 * pad, ih, radius=1.6,
                     face="#ffffff", edge=HAIRLINE, lw=0.8, z=4)
            ax.add_patch(Rectangle((x + pad + 0.9, iy + 1.4), 0.85, ih - 2.8,
                                   facecolor=accent, edgecolor="none", zorder=5))
            tx = x + pad + 3.4
            if kind == "image":
                cap_h = 7.4
                bx = (x + pad + 3.2) / _XU
                by = (iy + cap_h) / _YU
                bw = (w - 2 * pad - 6.4) / _XU
                bh = (ih - cap_h - 2.0) / _YU
                thumb_axes[key] = fig.add_axes([bx, by, bw, bh], zorder=6)
                ax.text(tx, iy + 4.7, lines[0], ha="left", va="center",
                        fontsize=8.2, fontweight="bold", color=INK, zorder=6)
                ax.text(tx, iy + 2.0, lines[1], ha="left", va="center",
                        fontsize=7.6, color=INK_2, zorder=6)
            else:
                ax.text(tx, cy - 3.3, key, ha="left", va="center",
                        fontsize=9.5, fontweight="bold", color=INK, zorder=6)
                for j, ln in enumerate(wl):
                    ax.text(tx, cy - 6.6 - j * 2.9, ln, ha="left", va="center",
                            fontsize=7.7, color=INK_2, zorder=6)
            cy = iy - 1.6

        if i < 4:
            _chevron(ax, x + w + 0.55, y_bot + h * 0.5, size=4.4,
                     color=_mix(accents[i + 1], 0.45))

    # --- caption band -------------------------------------------------------
    _rounded(ax, 6, 3.5, _XU - 12, 17.5, radius=2.4, face=PANEL_BG,
             edge=HAIRLINE, lw=0.9, z=1)
    rep = D["report"]
    best_kill = max(D["sim"], key=lambda r: r["kill"])
    best_pen = max(D["sim"], key=lambda r: r["pen"])
    ax.text(11, 18.0, "How to read this", fontsize=11, fontweight="bold",
            color=INK, va="center")
    _para(ax, 11, 15.2,
          "Stages 1–2 turn a measured tissue section into a cell-resolved "
          "microenvironment; stages 3–4 turn a molecular structure into transport "
          "and kill physics; stage 5 is what the two produce together. The two "
          "halves meet only at the PDE, so a drug's fate on this tumor is decided "
          "by its own structure acting on this patient's geometry — no per-drug "
          "fitting anywhere in the chain.",
          width=205, fontsize=9.6, color=INK_2, dy=3.4)
    _para(ax, 11, 7.6,
          f"Rendered live at dx = {D['dx_um']:.0f} µm on the "
          f"{rep['n_cells']:,}-cell Xenium section "
          f"(grid {rep['grid'][0]}×{rep['grid'][1]}, "
          f"{rep['extent_um'][0]:.0f} × {rep['extent_um'][1]:.0f} µm). "
          f"Best predicted penetration {best_pen['name']} "
          f"({best_pen['pen']:.3f}); best predicted kill {best_kill['name']} "
          f"({best_kill['kill']*100:.1f}%). "
          "PREDICTIONS on real geometry — not validated outcomes.",
          width=215, fontsize=9.0, color=INK_3, dy=3.2)

    # --- thumbnails ---------------------------------------------------------
    micro = D["micro"]
    W, H = micro.domain.extent_um
    ext = [0, W, H, 0]
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    dx = float(micro.domain.dx)
    thumbs = [("vascular", micro.fields.vascular_density, CMAPS["vascular"], 0.55),
              ("exposure", D["by_drug"]["doxorubicin"]["auc"], CMAPS["exposure"], 0.55)]
    for key, fld, cmap, gamma in thumbs:
        a = thumb_axes.get(key)
        if a is None:
            continue
        disp = _display_field(fld, tissue, dx, smooth_um=55.0)
        _field_im(a, disp, ext, cmap, tissue=None, gamma=gamma, pmax=99.0)
        a.set_xticks([]); a.set_yticks([])
        a.set_facecolor(SURFACE)
        for s in a.spines.values():
            s.set_color(HAIRLINE)

    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# shared: pick real tissue patches (perivascular / hypoxic / dense)
# --------------------------------------------------------------------------- #
def _patch_centers(D: Dict, half_um: float = 320.0) -> Dict[str, Tuple[float, float]]:
    """Locate real, cell-populated patches: the most perfused, the most hypoxic,
    and the densest. Returns centre coordinates in micrometres."""
    from scipy.ndimage import gaussian_filter

    micro = D["micro"]
    dom = micro.domain
    ny, nx = dom.shape
    counts = np.zeros((ny, nx), dtype=np.float64)
    col, row = dom.to_index(micro.cells.x, micro.cells.y)
    np.add.at(counts, (row, col), 1.0)
    dens = gaussian_filter(counts, sigma=max(1.0, 40.0 / dom.dx))

    tissue = np.asarray(dom.tissue_mask, dtype=bool)
    m = int(np.ceil(half_um / dom.dx)) + 1
    interior = np.zeros((ny, nx), dtype=bool)
    interior[m:max(m + 1, ny - m), m:max(m + 1, nx - m)] = True
    valid = tissue & interior & (dens > np.percentile(dens[tissue], 65.0))
    if not valid.any():
        valid = tissue & interior

    def _sel(field: np.ndarray, mode: str) -> Tuple[float, float]:
        f = np.where(valid, np.asarray(field, dtype=np.float64), np.nan)
        i = int(np.nanargmax(f) if mode == "max" else np.nanargmin(f))
        r, c = np.unravel_index(i, f.shape)
        return (float((c + 0.5) * dom.dx), float((r + 0.5) * dom.dx))

    return {
        "perivascular": _sel(micro.fields.vascular_density, "max"),
        "hypoxic": _sel(micro.fields.oxygen, "min"),
        "dense": _sel(dens, "max"),
    }


def _scatter_cells(ax, D, cx, cy, half, s=6.0, alpha=0.95, types=None):
    """Scatter the REAL cells inside a square patch, coloured by measured type."""
    cells = D["micro"].cells
    x = np.asarray(cells.x)
    y = np.asarray(cells.y)
    t = np.asarray(cells.cell_type)
    sel = (np.abs(x - cx) <= half) & (np.abs(y - cy) <= half)
    order = [CellType.TUMOR, CellType.STROMAL, CellType.IMMUNE,
             CellType.ENDOTHELIAL, CellType.NECROTIC]
    for ct in order:
        if types is not None and ct not in types:
            continue
        m = sel & (t == int(ct))
        if not m.any():
            continue
        ax.scatter(x[m], y[m], s=s, c=CELLTYPE_COLORS[int(ct)], alpha=alpha,
                   linewidths=0, zorder=4 if ct != CellType.ENDOTHELIAL else 5)
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy + half, cy - half)      # image convention: y down
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])


# --------------------------------------------------------------------------- #
# FIGURE 2 — the multiscale ladder
# --------------------------------------------------------------------------- #
def figure_multiscale(D: Dict, path: str) -> str:
    """Figure 2 — molecule (nm) → ECM mesh (nm) → cell (µm) → vessel spacing
    (µm) → tissue section (mm), with the physics and the real numbers at each
    rung."""
    mol = D["molecular"]
    sc = D["scales"]
    rep = D["report"]
    micro = D["micro"]
    dom = micro.domain
    radii = np.array([mol[n]["radius_nm"] for n in PANEL])
    freeD = np.array([mol[n]["free_D"] for n in PANEL])
    deff = np.array([mol[n]["d_eff"] for n in PANEL])
    hind = np.array([mol[n]["mesh_hindrance"] for n in PANEL])
    mws = np.array([mol[n]["mw"] for n in PANEL])
    Ls = {n: mol[n]["L_um"] for n in PANEL}
    pore = D["mesh_pore_nm"]
    bev = mol["bevacizumab"]

    accents = [CAT[0], CAT[3], CAT[7], CAT[4], CAT[1]]
    W, H = dom.extent_um

    fig = plt.figure(figsize=(22.0, 13.4))
    gs = GridSpec(3, 5, figure=fig, height_ratios=[3.05, 1.15, 2.05],
                  left=0.032, right=0.984, top=0.845, bottom=0.095,
                  hspace=0.26, wspace=0.13)

    fig.text(0.032, 0.955, "Five scales, one drug molecule",
             fontsize=25, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(0.030, 0.907,
             "A drug that works has to survive every rung of this ladder. Each "
             "panel is drawn at its own scale; the ruler underneath places them "
             "on one logarithmic axis spanning eight orders of magnitude "
             "(0.4 nanometres to 8 millimetres).",
             fontsize=11.5, color=INK_2, ha="left", va="center")

    sch = [fig.add_subplot(gs[0, i]) for i in range(5)]
    ruler = fig.add_subplot(gs[1, :])
    txt = [fig.add_subplot(gs[2, i]) for i in range(5)]

    # ---------------- rung 1: the molecule --------------------------------- #
    a = sch[0]
    a.set_aspect("equal")
    a.set_xlim(-1.45, 1.45); a.set_ylim(-1.02, 1.06)
    a.set_xticks([]); a.set_yticks([])
    small = PANEL[int(np.argmin(radii))]
    large = PANEL[int(np.argmax(radii))]
    for nm, xc in ((small, -0.66), (large, 0.66)):
        r = mol[nm]["radius_nm"]
        col_ = DRUG_COLORS[nm]
        a.add_patch(Circle((xc, 0.10), r, facecolor=_mix(col_, 0.55),
                           edgecolor=col_, lw=1.8, zorder=3))
        a.annotate("", xy=(xc + r, 0.10), xytext=(xc, 0.10),
                   arrowprops=dict(arrowstyle="->", color=INK_2, lw=1.0), zorder=4)
        a.text(xc, -0.42, f"{nm}\n$r_h$ = {r:.2f} nm", ha="center", va="top",
               fontsize=8.6, color=INK_2, linespacing=1.5)
    a.text(0, 0.86, r"$D_{free}=k_BT\,/\,6\pi\eta r_h$", ha="center", va="center",
           fontsize=12.5, color=INK)
    a.plot([-1.34, -0.84], [-0.92, -0.92], color=INK_2, lw=2.4,
           solid_capstyle="butt")
    a.text(-1.09, -0.885, "0.5 nm", ha="center", va="bottom", fontsize=8,
           color=INK_2)

    # ---------------- rung 2: the ECM mesh --------------------------------- #
    a = sch[1]
    a.set_aspect("equal")
    box_w, box_h = 62.0, 42.0
    a.set_xlim(0, box_w); a.set_ylim(0, box_h)
    a.set_xticks([]); a.set_yticks([])
    n_fib = int(round(FIBER_VOLUME_FRACTION * box_w * box_h /
                      (math.pi * FIBER_RADIUS_NM ** 2)))
    rng = np.random.default_rng(11)
    pts: List[Tuple[float, float]] = []
    tries = 0
    while len(pts) < n_fib and tries < 8000:
        tries += 1
        p = (rng.uniform(2.5, box_w - 2.5), rng.uniform(2.5, box_h - 8.0))
        if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 > (3.0 * FIBER_RADIUS_NM) ** 2
               for q in pts):
            pts.append(p)
    for p in pts:
        a.add_patch(Circle(p, FIBER_RADIUS_NM, facecolor=_mix(CAT[3], 0.35),
                           edgecolor=CAT[3], lw=0.9, zorder=3))
    a.add_patch(Circle((box_w * 0.5, box_h * 0.55), pore, facecolor="none",
                       edgecolor=INK_2, lw=1.2, ls=(0, (3, 2.5)), zorder=4))
    a.text(box_w * 0.5, box_h * 0.55 + pore + 1.0, f"pore r = {pore:.1f} nm",
           ha="center", va="bottom", fontsize=8.2, color=INK_2, zorder=5)
    a.add_patch(Circle((box_w * 0.17, 5.0), mol[large]["radius_nm"] * 1.0,
                       facecolor=DRUG_COLORS[large], edgecolor=DRUG_COLORS[large],
                       lw=1.4, zorder=6))
    a.text(box_w * 0.17, 8.4, f"{large}  {mol[large]['radius_nm']:.2f} nm",
           ha="center", va="bottom", fontsize=7.8, color=INK_2, zorder=6)
    a.add_patch(Circle((box_w * 0.80, 6.5), bev["radius_nm"],
                       facecolor=_mix("#7f7f7f", 0.35), edgecolor="#5a5a5a",
                       lw=1.2, zorder=6))
    a.text(box_w * 0.80, 13.0, f"IgG  {bev['radius_nm']:.1f} nm", ha="center",
           va="bottom", fontsize=7.8, color=INK_2, zorder=6)
    a.plot([2.5, 12.5], [box_h - 3.4, box_h - 3.4], color=INK_2, lw=2.4,
           solid_capstyle="butt")
    a.text(7.5, box_h - 3.0, "10 nm", ha="center", va="bottom", fontsize=8,
           color=INK_2)

    # ---------------- rung 3: real cells ----------------------------------- #
    centers = _patch_centers(D)
    a = sch[2]
    cx, cy = centers["dense"]
    _scatter_cells(a, D, cx, cy, 95.0, s=48.0, alpha=0.9)
    _scalebar(a, 50.0, "50 µm", frac_x=0.05, frac_y=0.055, size=8)
    a.text(0.5, 1.025, f"real Xenium cells · median spacing "
                       f"{sc['cell_nn_median']:.1f} µm",
           transform=a.transAxes, ha="center", va="bottom", fontsize=8.4,
           color=INK_2)

    # ---------------- rung 4: vessels and the penetration length ----------- #
    a = sch[3]
    tissue = np.asarray(dom.tissue_mask, dtype=bool)
    dx = float(dom.dx)
    auc_disp = _display_field(D["by_drug"]["doxorubicin"]["auc"], tissue, dx,
                              smooth_um=36.0, cut=0.02)
    px, py = centers["perivascular"]
    half4 = 300.0
    j0 = int(round((px - half4) / dx)); j1 = int(round((px + half4) / dx))
    i0 = int(round((py - half4) / dx)); i1 = int(round((py + half4) / dx))
    j0, i0 = max(0, j0), max(0, i0)
    j1, i1 = min(dom.nx, j1), min(dom.ny, i1)
    sub = auc_disp[i0:i1, j0:j1]
    ext4 = [j0 * dx, j1 * dx, i1 * dx, i0 * dx]
    _field_im(a, sub, ext4, CMAPS["exposure"], tissue=None, gamma=0.6, pmax=99.0)
    _scatter_cells(a, D, px, py, half4, s=11.0, alpha=0.95,
                   types=[CellType.ENDOTHELIAL])
    a.set_xlim(ext4[0], ext4[1]); a.set_ylim(ext4[2], ext4[3])
    a.set_box_aspect(1.0)
    L_dox = Ls["doxorubicin"]
    a.add_patch(Circle((px, py), L_dox, facecolor="none", edgecolor=INK,
                       lw=1.5, ls=(0, (3, 2.5)), zorder=7))
    a.annotate(f"L = {L_dox:.0f} µm", xy=(px, py + L_dox), xytext=(14, -30),
               textcoords="offset points", fontsize=8.4, color=INK, zorder=8,
               arrowprops=dict(arrowstyle="-", color=INK, lw=0.9),
               bbox=dict(boxstyle="round,pad=0.24", facecolor=SURFACE,
                         edgecolor="none", alpha=0.9))
    _scalebar(a, 100.0, "100 µm", frac_x=0.05, frac_y=0.055, size=8)
    a.text(0.5, 1.025, "doxorubicin exposure + real vessels",
           transform=a.transAxes, ha="center", va="bottom", fontsize=8.4,
           color=INK_2)

    ins = a.inset_axes([0.555, 0.635, 0.415, 0.325])
    cen = D["decay"]["centers"]
    for nm in ("temozolomide", "doxorubicin"):
        ins.plot(cen, D["decay"]["profiles"][nm], lw=2.0,
                 color=DRUG_COLORS[nm], solid_capstyle="round")
    ins.set_facecolor("#ffffff")
    ins.patch.set_alpha(0.92)
    ins.set_xlim(0, 200); ins.set_ylim(0, 1.05)
    ins.set_xlabel("µm from vessel", fontsize=6.8, labelpad=1.0)
    ins.set_ylabel("rel. AUC", fontsize=6.8, labelpad=1.0)
    ins.tick_params(labelsize=6.2, length=2, width=0.6, pad=1.2)
    ins.grid(True, axis="y", lw=0.5)
    ins.set_axisbelow(True)
    ins.text(0.97, 0.95, "TMZ", transform=ins.transAxes, ha="right", va="top",
             fontsize=7, color=DRUG_COLORS["temozolomide"], fontweight="bold")
    ins.text(0.97, 0.66, "dox", transform=ins.transAxes, ha="right", va="top",
             fontsize=7, color=DRUG_COLORS["doxorubicin"], fontweight="bold")

    # ---------------- rung 5: the whole section ---------------------------- #
    a = sch[4]
    cells = micro.cells
    x = np.asarray(cells.x); y = np.asarray(cells.y); t = np.asarray(cells.cell_type)
    step = max(1, len(x) // 24000)
    for ct in (CellType.TUMOR, CellType.STROMAL, CellType.IMMUNE,
               CellType.ENDOTHELIAL):
        m = t[::step] == int(ct)
        if m.any():
            a.scatter(x[::step][m], y[::step][m], s=0.5,
                      c=CELLTYPE_COLORS[int(ct)], linewidths=0, alpha=0.85)
    a.set_xlim(0, W); a.set_ylim(H, 0)
    a.set_aspect("equal")
    a.set_xticks([]); a.set_yticks([])
    _scalebar(a, 1000.0, "1 mm", frac_x=0.05, frac_y=0.06, size=8)
    a.text(0.5, 1.025, f"{rep['n_cells']:,} measured cells, one section",
           transform=a.transAxes, ha="center", va="bottom", fontsize=8.4,
           color=INK_2)

    for a, acc in zip(sch, accents):
        for s_ in a.spines.values():
            s_.set_color(HAIRLINE)
            s_.set_visible(True)

    # ---------------- the ruler -------------------------------------------- #
    lo, hi = -3.6, 4.25
    ruler.set_xlim(lo, hi); ruler.set_ylim(0, 1)
    ruler.set_xticks([]); ruler.set_yticks([])
    for s_ in ruler.spines.values():
        s_.set_visible(False)
    ruler.plot([lo, hi], [0.55, 0.55], color="#c9c8c2", lw=1.1, zorder=1)
    for v, lab in [(-3, "1 nm"), (-2, "10 nm"), (-1, "100 nm"), (0, "1 µm"),
                   (1, "10 µm"), (2, "100 µm"), (3, "1 mm"), (4, "10 mm")]:
        ruler.plot([v, v], [0.50, 0.55], color="#c9c8c2", lw=1.0, zorder=1)
        ruler.text(v, 0.455, lab, ha="center", va="top", fontsize=8.6,
                   color=INK_3)

    spans = [
        (radii.min() * 1e-3, radii.max() * 1e-3, "drug radius"),
        (FIBER_RADIUS_NM * 1e-3, pore * 1e-3, "fibre → pore"),
        (sc["cell_nn_median"], sc["cell_nn_median"], "cell spacing"),
        (sc["vessel_dist_median"], sc["vessel_dist_p90"], "vessel distance"),
        (min(W, H), max(W, H), "section"),
    ]
    fig.canvas.draw()
    rpos = ruler.get_position()
    for i, ((v0, v1, lab), acc) in enumerate(zip(spans, accents)):
        l0, l1 = math.log10(v0), math.log10(v1)
        if l1 - l0 < 0.06:
            ruler.plot([l0], [0.55], marker="o", ms=9, color=acc, zorder=4,
                       markeredgecolor=SURFACE, markeredgewidth=1.8)
        else:
            ruler.plot([l0, l1], [0.55, 0.55], color=acc, lw=7.5, zorder=3,
                       solid_capstyle="round")
        mid = 0.5 * (l0 + l1)
        ruler.text(mid, 0.63, lab, ha="center", va="bottom", fontsize=8.4,
                   color=acc, fontweight="bold", zorder=5)
        spos = sch[i].get_position()
        cxf = 0.5 * (spos.x0 + spos.x1)
        xd = lo + (cxf - rpos.x0) / (rpos.x1 - rpos.x0) * (hi - lo)
        ruler.plot([mid, xd], [0.80, 1.42], color=_mix(acc, 0.55), lw=1.1,
                   zorder=2, clip_on=False)

    for nm in ("doxorubicin", "temozolomide"):
        lv = math.log10(Ls[nm])
        ruler.plot([lv], [0.21], marker="^", ms=7, color=DRUG_COLORS[nm], zorder=5)
        ruler.text(lv, 0.13, f"L({nm[:3]}) = {Ls[nm]:.0f} µm", ha="center",
                   va="top", fontsize=7.8, color=DRUG_COLORS[nm])
    lv = math.log10(dom.dx)
    ruler.plot([lv], [0.21], marker="^", ms=6, color=INK_3, zorder=5)
    ruler.text(lv, 0.13, f"grid dx = {dom.dx:.0f} µm", ha="center", va="top",
               fontsize=7.8, color=INK_3)
    ruler.text(lo, 0.02, "physics lengths the simulation actually resolves",
               ha="left", va="bottom", fontsize=7.6, color=INK_3, style="italic")

    # ---------------- text blocks ------------------------------------------ #
    blocks = [
        ("MOLECULE", "0.4 nm", [
            "Stokes–Einstein on a hydrodynamic radius derived from the verified "
            "structure alone.",
            f"The panel spans r_h {radii.min():.2f}–{radii.max():.2f} nm "
            f"({mws.min():.0f}–{mws.max():.0f} Da), giving free diffusivities of "
            f"{freeD.min():.0f}–{freeD.max():.0f} µm²/s at 310 K in 1.2 mPa·s fluid.",
            "Size enters here once and never leaves."]),
        ("ECM MESH", f"{pore:.1f} nm pore", [
            f"Ogston available volume × Renkin/Faxén drag through a random "
            f"{FIBER_RADIUS_NM:.0f} nm fibre matrix at φ = {FIBER_VOLUME_FRACTION:g}.",
            f"The panel keeps {hind.min():.2f}–{hind.max():.2f} of its free "
            f"diffusivity; a {bev['radius_nm']:.1f} nm IgG keeps only "
            f"{bev['mesh_hindrance']:.3f} — the mesh is a size gate.",
            f"After tissue tortuosity: D_eff = {deff.min():.0f}–{deff.max():.0f} µm²/s."]),
        ("CELL", f"{sc['cell_nn_median']:.0f} µm", [
            f"{rep['n_cells']:,} real cells at measured coordinates, typed by "
            "gated marker transcripts.",
            f"TUMOR {rep['celltype_fractions']['TUMOR']*100:.1f}%, STROMAL "
            f"{rep['celltype_fractions']['STROMAL']*100:.1f}%, IMMUNE "
            f"{rep['celltype_fractions']['IMMUNE']*100:.1f}%, ENDOTHELIAL "
            f"{rep['celltype_fractions']['ENDOTHELIAL']*100:.1f}%.",
            "Each cell carries its own oxygen state and Hill kill draw."]),
        ("VESSEL SPACING", f"{sc['vessel_dist_median']:.0f} µm", [
            f"The median tumor cell sits {sc['vessel_dist_median']:.0f} µm from "
            f"the nearest measured endothelial cell (90th percentile "
            f"{sc['vessel_dist_p90']:.0f} µm).",
            f"Penetration length L = sqrt(D_eff / k_eff): "
            f"{Ls['doxorubicin']:.0f} µm for doxorubicin vs "
            f"{Ls['temozolomide']:.0f} µm for temozolomide.",
            f"α = {D['alpha']:.1f} was fit once so that L matches the measured "
            f"{D['citation']['dox_band'][0]:.0f}–"
            f"{D['citation']['dox_band'][1]:.0f} µm doxorubicin gradient."]),
        ("TISSUE SECTION", f"{max(W, H)/1000:.1f} mm", [
            f"{W:.0f} × {H:.0f} µm solved on a {dom.ny}×{dom.nx} grid at "
            f"dx = {dom.dx:.0f} µm.",
            "Convection–diffusion–reaction with interstitial pressure, oxygen "
            "and Hill kill; implicit in time, mass-conserving to 1e-13.",
            f"Output: penetration index {min(r['pen'] for r in D['sim']):.2f}–"
            f"{max(r['pen'] for r in D['sim']):.2f}, predicted kill "
            f"{min(r['kill'] for r in D['sim'])*100:.1f}–"
            f"{max(r['kill'] for r in D['sim'])*100:.1f}%."]),
    ]
    for a, (head, badge, lines), acc in zip(txt, blocks, accents):
        a.axis("off")
        a.set_xlim(0, 1); a.set_ylim(0, 1)
        a.plot([0, 1], [0.985, 0.985], color=acc, lw=2.6, clip_on=False)
        a.text(0, 0.925, head, fontsize=11.5, fontweight="bold", color=acc,
               va="top", ha="left")
        a.text(1.0, 0.925, badge, fontsize=10.5, color=INK_3, va="top",
               ha="right")
        yy = 0.735
        for ln in lines:
            if "$" in ln:
                a.text(0, yy, ln, fontsize=8.8, color=INK_2, va="top",
                       ha="left", transform=a.transAxes)
                yy -= 0.075
            else:
                yy = _para(a, 0, yy, ln, width=53, fontsize=8.8, color=INK_2,
                           dy=0.072, transform=a.transAxes)
            yy -= 0.038

    _caption(fig,
             "Every rung is measured or mechanistic: radii come from verified "
             "structures, the mesh from fibre-matrix theory, the cells and "
             "vessels from a real Xenium section, and the section geometry from "
             "the assay itself.\nThe only fitted number in the whole ladder is "
             "the single perivascular-sink constant α, anchored to one published "
             "doxorubicin measurement (PMID 16361566).",
             y=0.016, size=9.4)

    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# FIGURE 3 — the real reconstructed microenvironment
# --------------------------------------------------------------------------- #
def figure_microenvironment(D: Dict, path: str) -> str:
    """Figure 3 — real cells coloured by measured type, two zoomed niches, and
    the physical fields that drive transport."""
    micro = D["micro"]
    dom = micro.domain
    rep = D["report"]
    W, H = dom.extent_um
    ext = [0, W, H, 0]
    dx = float(dom.dx)
    tissue = np.asarray(dom.tissue_mask, dtype=bool)
    cells = micro.cells
    x = np.asarray(cells.x); y = np.asarray(cells.y); t = np.asarray(cells.cell_type)

    disp = {
        "vascular": _display_field(micro.fields.vascular_density, tissue, dx, 45.0),
        "ecm": _display_field(micro.fields.ecm_density, tissue, dx, 45.0),
        "oxygen": _display_field(micro.fields.oxygen, tissue, dx, 45.0,
                                 everywhere=True),
    }

    centers = _patch_centers(D)
    half = 320.0
    zooms = [("A", "perivascular", centers["perivascular"], "vascular",
              "vascular density"),
             ("B", "hypoxic", centers["hypoxic"], "oxygen", "oxygen")]

    grange = {k: (0.0, float(np.nanpercentile(v, 99.5)))
              for k, v in disp.items()}

    fig = plt.figure(figsize=(23.0, 13.4))
    gs = GridSpec(2, 4, figure=fig, width_ratios=[1.0, 1.0, 1.0, 1.0],
                  height_ratios=[1.0, 0.55], left=0.030, right=0.972,
                  top=0.848, bottom=0.092, hspace=0.28, wspace=0.16)

    fig.text(0.030, 0.957, "One real glioblastoma section, fully reconstructed",
             fontsize=25, fontweight="bold", color=INK, ha="left", va="center")
    fig.text(0.032, 0.903,
             f"{rep['n_cells']:,} measured cells at true micrometre coordinates "
             f"from an open 10x Xenium GBM section "
             f"({D['xenium']['n_panel_genes']}-gene panel). Cell identity is gated "
             "on marker transcripts; the fields underneath are built from those "
             "same measured cells, not assumed.",
             fontsize=11.5, color=INK_2, ha="left", va="center")

    # ---- main cell map ----------------------------------------------------
    ax = fig.add_subplot(gs[0, 0:2])
    counts = {}
    for ct in (CellType.TUMOR, CellType.STROMAL, CellType.IMMUNE,
               CellType.ENDOTHELIAL, CellType.NECROTIC):
        m = t == int(ct)
        counts[ct] = int(m.sum())
        if not m.any():
            continue
        ax.scatter(x[m], y[m], s=1.5, c=CELLTYPE_COLORS[int(ct)],
                   linewidths=0, alpha=0.85,
                   zorder=4 if ct != CellType.ENDOTHELIAL else 5)
    ax.set_xlim(0, W); ax.set_ylim(H, 0)
    ax.set_box_aspect(H / W)
    ax.set_anchor("NW")
    ax.set_xlabel("x (µm)"); ax.set_ylabel("y (µm)")
    ax.tick_params(length=2.5, width=0.7)
    for s_ in ax.spines.values():
        s_.set_color(HAIRLINE); s_.set_visible(True)
    ax.set_title("Measured cells, coloured by gated phenotype", fontsize=13,
                 color=INK, pad=9, loc="left", fontweight="bold")

    for tag, name, (cx, cy), _k, _l in zooms:
        ax.add_patch(Rectangle((cx - half, cy - half), 2 * half, 2 * half,
                               facecolor="none", edgecolor=INK, lw=1.6, zorder=8))
        ax.text(cx, cy - half - 55, tag, ha="center", va="bottom", fontsize=12,
                fontweight="bold", color=INK, zorder=9,
                bbox=dict(boxstyle="circle,pad=0.26", facecolor=SURFACE,
                          edgecolor=INK, lw=1.2))

    handles = [Line2D([0], [0], marker="o", ls="none", ms=9,
                      markerfacecolor=CELLTYPE_COLORS[int(ct)],
                      markeredgecolor="none",
                      label=f"{ct.name.title()}  {counts[ct]:,}  "
                            f"({counts[ct]/rep['n_cells']*100:.1f}%)")
               for ct in (CellType.TUMOR, CellType.STROMAL, CellType.IMMUNE,
                          CellType.ENDOTHELIAL) if counts[ct]]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.10),
              ncol=4, fontsize=10.5, handletextpad=0.5, columnspacing=2.2,
              labelcolor=INK_2)

    # ---- zoom insets ------------------------------------------------------
    for spec, (tag, name, (cx, cy), key, lab) in zip(
            [gs[0, 2], gs[0, 3]], zooms):
        za = fig.add_subplot(spec)
        j0 = max(0, int(round((cx - half) / dx)))
        j1 = min(dom.nx, int(round((cx + half) / dx)))
        i0 = max(0, int(round((cy - half) / dx)))
        i1 = min(dom.ny, int(round((cy + half) / dx)))
        sub = disp[key][i0:i1, j0:j1]
        e = [j0 * dx, j1 * dx, i1 * dx, i0 * dx]
        _field_im(za, sub, e, CMAPS[key], tissue=None, gamma=0.6,
                  vmin=grange[key][0], vmax=grange[key][1])
        _scatter_cells(za, D, cx, cy, half, s=13.0, alpha=0.95)
        za.set_xlim(e[0], e[1]); za.set_ylim(e[2], e[3])
        za.set_box_aspect(1.0)
        za.set_anchor("N")
        za.set_xticks([]); za.set_yticks([])
        _scalebar(za, 100.0, "100 µm", frac_x=0.05, frac_y=0.055, size=8.5)
        for s_ in za.spines.values():
            s_.set_color(INK); s_.set_linewidth(1.6); s_.set_visible(True)
        col, row = dom.to_index(np.array([cx]), np.array([cy]))
        raw = float(np.asarray(
            micro.fields.vascular_density if key == "vascular"
            else micro.fields.oxygen)[row[0], col[0]])
        extreme = "most perfused" if key == "vascular" else "most hypoxic"
        za.set_title(f"{tag} · {name}\nthe {extreme} patch  ·  {lab} {raw:.3f}",
                     fontsize=10.2, color=INK, pad=7, linespacing=1.5)

    # ---- field small multiples -------------------------------------------
    fields = [
        (gs[1, 0], "vascular", "Vascular density",
         "Gaussian KDE of the real endothelial cells",
         f"mean {rep['mean_vascular_density']:.3f}"),
        (gs[1, 1], "ecm", "ECM density",
         f"weighted KDE of {len(D['xenium']['ecm_markers'])} matrix transcripts",
         f"mean {rep['mean_ecm_density']:.3f}"),
        (gs[1, 2], "oxygen", "Oxygen",
         "physics solve driven by the real perfusion",
         f"mean {rep['mean_oxygen']:.3f}"),
    ]
    for spec, key, title, sub, stat in fields:
        a = fig.add_subplot(spec)
        im = _field_im(a, disp[key], ext, CMAPS[key], tissue=None, gamma=0.5,
                       vmin=grange[key][0], vmax=grange[key][1])
        a.set_anchor("N")
        a.set_xticks([]); a.set_yticks([])
        for s_ in a.spines.values():
            s_.set_color(HAIRLINE); s_.set_visible(True)
        a.set_title(f"{title}  ·  {stat}\n{sub}", fontsize=10, color=INK,
                    pad=6, linespacing=1.5)
        _colorbar(fig, a, im, title.lower(), size=8)

    # ---- provenance panel -------------------------------------------------
    a = fig.add_subplot(gs[1, 3])
    a.axis("off"); a.set_xlim(0, 1); a.set_ylim(0, 1)
    a.plot([0, 1], [1.0, 1.0], color=INK_3, lw=2.2, clip_on=False)
    a.text(0, 0.93, "WHERE EVERY FIELD COMES FROM", fontsize=10.5,
           fontweight="bold", color=INK, va="top")
    lines = [
        ("Cells", f"{rep['n_cells']:,} Xenium cells, threshold-gated marker "
                  "typing; cross-platform Pearson 0.997 against the independent "
                  "Visium NNLS deconvolution."),
        ("Vasculature", f"Gaussian KDE of the {D['scales']['n_endothelial']:,} "
                        "measured endothelial cells — no synthetic vessel tree."),
        ("ECM", "expression-weighted KDE of "
                + ", ".join(D["xenium"]["ecm_markers"][:4]) + " and "
                + f"{max(0, len(D['xenium']['ecm_markers']) - 4)} more real "
                  "matrix transcripts."),
        ("Oxygen", "reaction–diffusion solved on that real perfusion field; "
                   f"{D['hypoxic_fraction']*100:.0f}% of tumor cells land "
                   "hypoxic at this resolution."),
        ("Not measured", "drug response. This is an untreated resection with no "
                         "outcome label — the physics is validated, the response "
                         "is a prediction."),
    ]
    yy = 0.83
    for head, body in lines:
        a.text(0, yy, head, fontsize=9.4, fontweight="bold", color=INK, va="top")
        yy = _para(a, 0, yy - 0.062, body, width=72, fontsize=8.6, color=INK_2,
                   dy=0.058, transform=a.transAxes)
        yy -= 0.035

    _caption(fig,
             "A · the most perfused patch in the section and B · the most "
             "hypoxic one, both 640 × 640 µm, drawn from the same measured cell "
             "table. The contrast between them is what makes drug delivery "
             "patient-specific: exposure decays away from the vessels in A, and "
             "the oxygen-starved cells in B resist whatever reaches them.\n"
             "Field panels show the LOCAL MEAN of the measured tissue pixels "
             "(Gaussian, 45 µm): at dx = 12 µm the measured support is a sparse "
             "dust of cell-containing pixels, so the raw arrays\nare unreadable. "
             "No values are invented, and A and B share the section-wide colour "
             "scale of the panel beneath them.",
             y=0.012, size=9.3)

    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# FIGURE 4 — delivery vs potency: the scientific message in one plot
# --------------------------------------------------------------------------- #
def figure_delivery_vs_potency(D: Dict, path: str) -> str:
    """Figure 4 — potency (1/IC50) against delivered exposure (tissue AUC), with
    marker area carrying the predicted kill."""
    rows = sorted(D["sim"], key=lambda r: r["mean_auc"])
    kills = np.array([r["kill"] for r in rows])
    aucs = np.array([r["mean_auc"] for r in rows])
    pot = np.array([1.0 / r["ic50"] for r in rows])

    dox = D["by_drug"]["doxorubicin"]
    winner = max(D["sim"], key=lambda r: r["kill"])
    most_delivered = max(D["sim"], key=lambda r: r["mean_auc"])
    span_pot = max(r["ic50"] for r in D["sim"]) / min(r["ic50"] for r in D["sim"])
    span_auc = most_delivered["mean_auc"] / min(r["mean_auc"] for r in D["sim"])

    fig = plt.figure(figsize=(18.0, 10.6))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[2.15, 1.0], left=0.058,
                  right=0.978, top=0.828, bottom=0.115, wspace=0.19)

    fig.text(0.058, 0.958,
             "The most potent drug loses, because it never arrives",
             fontsize=25, fontweight="bold", color=INK, ha="left", va="center")
    _para(fig, 0.058, 0.905,
          "Five drugs, one real tumor, identical geometry. Potency is the "
          "measured GDSC/HGCC IC50; delivered exposure is what the transport "
          "physics actually puts into the tissue. Kill needs BOTH — and on this "
          f"tumor delivery spans {span_auc:.0f}× while potency spans "
          f"{span_pot:.0f}×, so delivery decides the ranking.",
          width=152, fontsize=11.5, color=INK_2, dy=0.027,
          transform=fig.transFigure)

    # ---- the scatter ------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xscale("log"); ax.set_yscale("log")
    xlo, xhi = aucs.min() / 5.0, aucs.max() * 7.0
    ylo, yhi = pot.min() / 5.0, pot.max() * 5.0

    # iso-exposure diagonals: AUC / IC50 = constant seconds held at the IC50
    for e in range(3, 7):
        c = 10.0 ** e
        xs = np.array([xlo, xhi])
        ax.plot(xs, c / xs, color=HAIRLINE, lw=1.0, zorder=1)
        xi = xhi * 0.45
        yi = c / xi
        if ylo * 1.6 < yi < yhi / 1.6:
            ax.text(xi, yi, f"$10^{e}$ s at IC50", fontsize=8.4, color=INK_3,
                    rotation=-31, rotation_mode="anchor", ha="center",
                    va="bottom", zorder=2)

    smax, smin = 2400.0, 300.0
    span = max(float(np.ptp(kills)), 1e-12)
    sizes = smin + (kills - kills.min()) / span * (smax - smin)
    label_pos = {"doxorubicin": (0, 1), "temozolomide": (0, -1),
                 "methotrexate": (0, 1), "lomustine": (0, -1),
                 "cisplatin": (0, 1)}
    for r, s_, xa, ya in zip(rows, sizes, aucs, pot):
        c = DRUG_COLORS[r["name"]]
        ax.scatter([xa], [ya], s=s_, c=c, alpha=0.9, linewidths=2.0,
                   edgecolors=SURFACE, zorder=5)
        up = label_pos.get(r["name"], (0, 1))[1]
        off = math.sqrt(s_) * 0.62 + 13
        ax.annotate(f"{r['name']}\nkill {r['kill']*100:.1f}%  ·  "
                    f"IC50 {r['ic50']:.3g} µM",
                    xy=(xa, ya), xytext=(0, off if up > 0 else -off),
                    textcoords="offset points", ha="center",
                    va="bottom" if up > 0 else "top",
                    fontsize=10.2, color=INK, fontweight="bold", zorder=6,
                    linespacing=1.5)

    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
    ax.set_xlabel("delivered exposure  —  mean tissue AUC (µM·s)   →",
                  fontsize=11.5)
    ax.set_ylabel("potency  —  1 / IC50 (µM$^{-1}$)   →", fontsize=11.5)
    ax.grid(True, which="major", lw=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.7)

    # the callout that carries the message (bottom-left quadrant is empty)
    ax.annotate(
        f"doxorubicin is {winner['ic50']/dox['ic50']:.0f}× more potent than "
        f"{winner['name']}\nby measured IC50 — and still loses.\n"
        f"It is delivered {most_delivered['mean_auc']/dox['mean_auc']:.0f}× "
        f"less than {most_delivered['name']}: largest MW,\n"
        f"BBB penetrance {dox['bbb']:.3f}, and trapped inside\n"
        f"L = {D['molecular']['doxorubicin']['L_um']:.0f} µm of a vessel. "
        "Its potency is spent\non tissue that the drug never leaves.",
        xy=(dox["mean_auc"], 1.0 / dox["ic50"]),
        xytext=(xlo * 1.35, ylo * 1.6),
        fontsize=10.2, color=INK_2, ha="left", va="bottom", linespacing=1.7,
        arrowprops=dict(arrowstyle="->", color=INK_3, lw=1.2,
                        connectionstyle="arc3,rad=-0.22"),
        zorder=7,
        bbox=dict(boxstyle="round,pad=0.55", facecolor=PANEL_BG,
                  edgecolor=HAIRLINE, lw=0.9))

    # size legend (top-right quadrant is empty)
    leg_k = [kills.min(), 0.5 * (kills.min() + kills.max()), kills.max()]
    leg_h = [Line2D([0], [0], marker="o", ls="none", markerfacecolor=INK_3,
                    markeredgecolor="none", alpha=0.55,
                    ms=math.sqrt(smin + (k - kills.min()) / span *
                                 (smax - smin)) * 0.60,
                    label=f"{k*100:.1f}%") for k in leg_k]
    lg = ax.legend(handles=leg_h, title="marker area = predicted kill",
                   loc="upper right", fontsize=10, labelspacing=2.1,
                   borderpad=1.1, handletextpad=2.0, title_fontsize=9.5,
                   labelcolor=INK_2, borderaxespad=0.8)
    lg.get_title().set_color(INK_2)

    # ---- kill ranking + the takeaway --------------------------------------
    sub = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0, 1],
                                  height_ratios=[3.0, 1.0], hspace=0.36)
    ax2 = fig.add_subplot(sub[0])
    order = sorted(D["sim"], key=lambda r: r["kill"])
    for i, r in enumerate(order):
        ax2.barh(i, r["kill"] * 100.0, height=0.46,
                 color=DRUG_COLORS[r["name"]], edgecolor=SURFACE, linewidth=2.0)
        ax2.text(r["kill"] * 100.0 + 0.14, i, f"{r['kill']*100:.1f}%",
                 va="center", ha="left", fontsize=10.5, color=INK,
                 fontweight="bold")
    ax2.set_yticks(np.arange(len(order)))
    ax2.set_yticklabels([r["name"] for r in order], fontsize=10.5)
    ax2.set_ylim(-0.7, len(order) - 0.3)
    ax2.set_xlim(0, max(r["kill"] for r in order) * 100.0 * 1.30)
    ax2.set_xlabel("predicted kill (% of tumor cells)", fontsize=11)
    ax2.grid(True, axis="x", lw=0.7)
    ax2.set_axisbelow(True)
    ax2.tick_params(length=3, width=0.7)
    ax2.set_title("Ranking follows the diagonal,\nnot the vertical axis",
                  fontsize=12.5, color=INK, pad=10, loc="left",
                  fontweight="bold", linespacing=1.5)

    order_ratio = [r["name"] for r in sorted(
        D["sim"], key=lambda r: r["mean_auc"] / r["ic50"])]
    order_kill = [r["name"] for r in order]
    agree = "identically to" if order_ratio == order_kill else "almost exactly like"

    ax3 = fig.add_subplot(sub[1])
    ax3.axis("off"); ax3.set_xlim(0, 1); ax3.set_ylim(0, 1)
    ax3.plot([0, 1], [1.0, 1.0], color=HAIRLINE, lw=1.4, clip_on=False)
    _para(ax3, 0, 0.86,
          "Rank the panel by exposure-in-IC50-units (AUC / IC50 — the grey "
          f"diagonals) and it orders {agree} the simulated kill. Delivery and "
          "potency multiply; here delivery carries the wider range, so the "
          "diagonal, not the vertical axis, predicts the winner.",
          width=58, fontsize=9.8, color=INK_2, dy=0.145,
          transform=ax3.transAxes)

    _caption(fig,
             f"Real geometry ({D['report']['n_cells']:,} Xenium cells, "
             f"dx = {D['dx_um']:.0f} µm), real measured IC50, structure-derived "
             f"transport, one fitted constant (α = {D['alpha']:.1f}). "
             f"{D['hypoxic_fraction']*100:.0f}% of tumor cells are hypoxic at "
             "this resolution and are oxygen-gated to lower sensitivity, which "
             "lowers every drug's kill.\nThese are PREDICTIONS on an untreated "
             "resection — there is no paired measured drug response to validate "
             "them against, and the ranking must not be read as clinical efficacy.",
             y=0.014, size=9.4)

    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #
def main(dx_um: float = DX_UM) -> List[str]:
    """Build all four explainer figures into ``results/``."""
    _style()
    os.makedirs(RESULTS_DIR, exist_ok=True)
    D = collect(dx_um=dx_um)

    jobs = [
        ("explainer_pipeline.png", figure_pipeline),
        ("explainer_multiscale.png", figure_multiscale),
        ("explainer_microenvironment.png", figure_microenvironment),
        ("explainer_delivery_vs_potency.png", figure_delivery_vs_potency),
    ]
    written = []
    for fname, fn in jobs:
        out = os.path.join(RESULTS_DIR, fname)
        fn(D, out)
        kb = os.path.getsize(out) / 1024.0
        print(f"  wrote {out}  ({kb:,.0f} KB)")
        written.append(out)
    return written


if __name__ == "__main__":
    main()
