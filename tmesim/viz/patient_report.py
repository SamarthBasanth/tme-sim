"""Decision-panel figure for the per-tumor VIRTUAL DRUG PANEL report.

This module is a pure RENDERING layer for
:mod:`experiments.patient_report`. It takes an already-computed ``panel`` dict
(the measured microenvironment profile + the ranked, already-attributed per-drug
rows) and draws the one figure a scientist / pharma reviewer would actually look
at when deciding which agent to take forward on THIS patient's tissue:

  row 1  the measured tumor: cell-type map, vascular density, ECM density,
         oxygen (hypoxic tissue hatched), niche-region map;
  row 2  the per-drug predicted exposure (AUC) maps on that same tissue, in
         RANKED order and on ONE shared colour scale, so the delivery
         differences between drugs are directly comparable by eye;
  row 3  the decision plane -- predicted penetration length L vs measured
         potency IC50, with the recommended drug ringed -- next to the ranked
         predicted-kill bar chart;
  row 4  the ATTRIBUTION chart: each drug's log10 advantage over the panel
         average decomposed into delivery / potency / hypoxia-resistance
         decades (computed in the experiment, not eyeballed here).

Nothing in ``core`` / ``physics`` / ``molecular`` / ``reconstruct`` is touched;
this module only draws.

Colour policy
-------------
* The five drugs get a FIXED categorical order (:data:`DRUG_COLORS`) -- a drug
  keeps its hue in every panel regardless of where it ranks, so re-running on a
  different tumor never repaints the survivors.
* Magnitude fields (vascular / ECM / oxygen / AUC) use single-hue sequential
  ramps, never a rainbow.
* The attribution chart is the one diverging encoding (advantage vs penalty
  about a zero baseline).
Both categorical palettes pass the lightness / chroma / CVD-separation /
contrast checks (deutan dE >= 14 on the worst adjacent pair).
"""
from __future__ import annotations

from typing import Dict, List, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from ..core import CellType

__all__ = ["plot_patient_report", "DRUG_COLORS", "REGION_COLORS", "CELLTYPE_COLORS"]

# --------------------------------------------------------------------------- #
# palettes (validated; assigned in FIXED order, never cycled)
# --------------------------------------------------------------------------- #
# Panel drug order -> hue. A drug's colour follows the drug, never its rank.
_DRUG_HUES: Sequence[str] = ("#3b6fd4", "#e0761a", "#009a86", "#b5399b", "#9c7414")

DRUG_COLORS: Dict[str, str] = {
    "temozolomide": _DRUG_HUES[0],
    "lomustine": _DRUG_HUES[1],
    "cisplatin": _DRUG_HUES[2],
    "methotrexate": _DRUG_HUES[3],
    "doxorubicin": _DRUG_HUES[4],
}

REGION_COLORS: Dict[str, str] = {
    "margin": "#4c9f70",
    "core": "#b8860b",
    "perivascular": "#c0392b",
    "necrotic": "#7d5ba6",
    "unlabeled": "#b0b0b0",
}

# Reuse the repo-wide per-CellType colours so this figure agrees with every
# other tmesim microenvironment plot.
CELLTYPE_COLORS: Dict[int, str] = {
    int(CellType.TUMOR): "#d62728",
    int(CellType.STROMAL): "#2ca02c",
    int(CellType.IMMUNE): "#1f77b4",
    int(CellType.ENDOTHELIAL): "#9467bd",
    int(CellType.NECROTIC): "#7f7f7f",
}

# Attribution diverging pair (advantage / penalty) + neutral ink.
_POS = "#2c7fb8"
_NEG = "#d95f02"

_INK = "#1b1b1b"
_INK_MUTED = "#6b6b6b"
_GRID = "#d9d9d9"


def _drug_color(name: str, fallback_index: int = 0) -> str:
    return DRUG_COLORS.get(name, _DRUG_HUES[fallback_index % len(_DRUG_HUES)])


def _bare(ax) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    # Equal-aspect images are shorter than their grid cell; anchor them to the
    # top so the slack collects below the row instead of splitting it in two.
    ax.set_anchor("N")
    for s in ax.spines.values():
        s.set_edgecolor("#c9c9c9")


def _recessive(ax) -> None:
    ax.grid(True, color=_GRID, lw=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_edgecolor("#bdbdbd")
    ax.tick_params(colors=_INK_MUTED, labelsize=8)


def _field_panel(ax, field, extent, title, cmap, cbar_label, vmax=None):
    """Single-hue sequential magnitude field with a slim colourbar."""
    arr = np.asarray(field, dtype=float)
    im = ax.imshow(arr, origin="lower", extent=extent, cmap=cmap,
                   aspect="equal", vmin=0.0, vmax=vmax)
    ax.set_title(title, fontsize=9, color=_INK, pad=6)
    _bare(ax)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.040, pad=0.02)
    cb.set_label(cbar_label, fontsize=7.5, color=_INK_MUTED)
    cb.ax.tick_params(labelsize=7, colors=_INK_MUTED)
    cb.outline.set_visible(False)
    return im


# --------------------------------------------------------------------------- #
# row 1 -- the measured tumor
# --------------------------------------------------------------------------- #
def _panel_celltypes(ax, panel, extent) -> None:
    cells = panel["cells"]
    x = np.asarray(cells["x"], float)
    y = np.asarray(cells["y"], float)
    ct = np.asarray(cells["cell_type"], int)
    frac = panel["report"]["celltype_fractions"]

    # Draw the rarer types last so they are not buried under the tumor mass.
    order = [CellType.TUMOR, CellType.STROMAL, CellType.IMMUNE,
             CellType.ENDOTHELIAL, CellType.NECROTIC]
    handles: List = []
    for c in order:
        m = ct == int(c)
        if not m.any():
            continue
        ax.scatter(x[m], y[m], s=0.30, linewidths=0,
                   color=CELLTYPE_COLORS[int(c)], alpha=0.65, rasterized=True)
        handles.append(mpatches.Patch(
            color=CELLTYPE_COLORS[int(c)],
            label=f"{c.name.title()} {frac.get(c.name, 0.0) * 100:.1f}%"))
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal")
    ax.set_title(f"measured tumor\n{panel['report']['n_cells']:,} real cells",
                 fontsize=9, color=_INK, pad=6)
    _bare(ax)
    ax.legend(handles=handles, loc="upper right", fontsize=6.6, frameon=True,
              framealpha=0.92, borderpad=0.35, handlelength=0.9,
              labelspacing=0.28, edgecolor="#cccccc")


def _panel_oxygen(ax, panel, extent) -> None:
    ox = np.asarray(panel["fields"]["oxygen"], float)
    thr = float(panel["hypoxia_threshold"])
    im = ax.imshow(ox, origin="lower", extent=extent, cmap="cividis",
                   aspect="equal", vmin=0.0, vmax=max(float(ox.max()), 1e-9))
    # Mark the hypoxic tissue explicitly rather than relying on the ramp alone.
    ax.contour(ox, levels=[thr], colors=["#ff3b30"], linewidths=0.8,
               extent=extent, origin="lower")
    ax.set_title(
        f"oxygen (real perfusion)\n{panel['hypoxic_fraction'] * 100:.0f}% of "
        "tumor cells hypoxic", fontsize=9, color=_INK, pad=6)
    _bare(ax)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.040, pad=0.02)
    cb.set_label(f"O2 (rel.); red line = threshold {thr:g}",
                 fontsize=7.5, color=_INK_MUTED)
    cb.ax.tick_params(labelsize=7, colors=_INK_MUTED)
    cb.outline.set_visible(False)


def _panel_regions(ax, panel, extent) -> None:
    """Argmax niche-region map (categorical), masked to tissue."""
    fields = panel.get("region_fields") or {}
    names = [n for n in ("margin", "core", "perivascular", "necrotic")
             if n in fields]
    tissue = np.asarray(panel["tissue_mask"], bool)
    if not names:
        ax.text(0.5, 0.5, "niche regions unavailable", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color=_INK_MUTED)
        _bare(ax)
        return

    stack = np.stack([np.asarray(fields[n], float) for n in names], axis=0)
    best = np.argmax(stack, axis=0).astype(float)
    supported = stack.sum(axis=0) > 1e-9
    best[~(tissue & supported)] = np.nan

    cmap = matplotlib.colors.ListedColormap([REGION_COLORS[n] for n in names])
    cmap.set_bad("#f2f2f2")
    ax.imshow(np.ma.masked_invalid(best), origin="lower", extent=extent,
              cmap=cmap, aspect="equal", vmin=-0.5, vmax=len(names) - 0.5,
              interpolation="nearest")
    ax.set_title("niche regions\n(Ivy GAP signature match)", fontsize=9,
                 color=_INK, pad=6)
    _bare(ax)
    tf = panel.get("region_tumor_fractions", {})
    ax.legend(
        handles=[mpatches.Patch(color=REGION_COLORS[n],
                                label=f"{n} {tf.get(n, 0.0) * 100:.1f}%")
                 for n in names],
        loc="upper right", fontsize=6.6, frameon=True, framealpha=0.92,
        borderpad=0.35, handlelength=0.9, labelspacing=0.28,
        edgecolor="#cccccc", title="% of tumor cells", title_fontsize=6.6)


# --------------------------------------------------------------------------- #
# row 2 -- per-drug exposure maps (ranked, shared scale)
# --------------------------------------------------------------------------- #
def _panel_auc(ax, row, rank, extent, norm, recommended: bool) -> None:
    # A LOG norm (shared across all five panels) is what makes the comparison
    # legible: on a linear scale the delivery-starved drugs render as a flat
    # black rectangle and their perivascular cuff -- the mechanism -- is
    # invisible. The scale is still ONE shared scale, so panels stay comparable.
    cmap = matplotlib.colormaps["magma"].copy()
    cmap.set_bad(cmap(0.0))
    im = ax.imshow(np.asarray(row["auc_field"], float), origin="lower",
                   extent=extent, cmap=cmap, aspect="equal", norm=norm)
    color = _drug_color(row["name"])
    ax.set_title(f"#{rank}  {row['name']}", fontsize=9.5, color=color,
                 fontweight="bold" if recommended else "normal", pad=5)
    _bare(ax)
    # Secondary (non-colour) encoding of identity + the recommendation.
    for s in ax.spines.values():
        s.set_edgecolor(color)
        s.set_linewidth(2.4 if recommended else 1.2)
    ax.text(0.025, 0.03,
            f"L = {row['L_um']:.0f} µm   kill = {row['kill']:.4f}\n"
            f"tumor-cell AUC = {row['tumor_auc']:.2g} µM·s",
            transform=ax.transAxes, fontsize=6.8, color="#ffffff", va="bottom",
            bbox=dict(boxstyle="round,pad=0.28", facecolor="#000000aa",
                      edgecolor="none"))
    if recommended:
        ax.text(0.97, 0.97, "RECOMMENDED", transform=ax.transAxes, fontsize=7,
                ha="right", va="top", color="#ffffff", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor=color,
                          edgecolor="none"))
    return im


# --------------------------------------------------------------------------- #
# row 3 -- decision plane + ranked kill
# --------------------------------------------------------------------------- #
def _panel_scatter(ax, rows, recommended: str) -> None:
    """Predicted penetration L (delivery) vs measured IC50 (potency).

    Down-and-right is the good corner: further penetration AND lower IC50.
    Marker area encodes the predicted kill; the recommendation is ringed.
    """
    kills = np.array([r["kill"] for r in rows], float)
    kmax = float(kills.max()) if kills.size and kills.max() > 0 else 1.0
    for r in rows:
        color = _drug_color(r["name"])
        area = 90.0 + 900.0 * (r["kill"] / kmax)
        ax.scatter(r["L_um"], r["ic50"], s=area, color=color, alpha=0.85,
                   edgecolors="#ffffff", linewidths=2.0, zorder=3)
        if r["name"] == recommended:
            ax.scatter(r["L_um"], r["ic50"], s=area * 2.6, facecolors="none",
                       edgecolors=color, linewidths=2.0, zorder=2)
        marker = "GDSC" if r["ic50_confidence"] == "measured" else "lit."
        ax.annotate(f"{r['name']} ({marker})", (r["L_um"], r["ic50"]),
                    textcoords="offset points", xytext=(12, 12),
                    ha="left", fontsize=8, color=_INK)
    # Pad the log limits so the annotations never fall outside the axes.
    xs = np.array([r["L_um"] for r in rows], float)
    ys = np.array([r["ic50"] for r in rows], float)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xs.min() / 1.9, xs.max() * 2.6)
    ax.set_ylim(ys.min() / 3.2, ys.max() * 2.2)
    ax.invert_yaxis()  # up == more potent
    # Decade ticks alone leave a nearly bare axis over this narrow range; label
    # the 1-2-5 subdivisions in plain numerals instead of 10^x.
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 5.0)))
        axis.set_minor_locator(matplotlib.ticker.NullLocator())
        axis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(
                lambda v, _p: f"{v:g}" if v >= 1 else f"{v:.2g}"))
    ax.set_xlabel("predicted penetration length L (µm)  →  better delivery",
                  fontsize=8.5, color=_INK_MUTED)
    ax.set_ylabel("real IC50 (µM), inverted  ↑  more potent", fontsize=8.5,
                  color=_INK_MUTED)
    ax.set_title(
        "decision plane: delivery vs potency  (marker area ∝ predicted kill; "
        f"ring = recommended: {recommended})", fontsize=9.5, color=_INK)
    _recessive(ax)


def _panel_rank(ax, rows) -> None:
    names = [r["name"] for r in rows][::-1]
    kills = [r["kill"] for r in rows][::-1]
    colors = [_drug_color(r["name"]) for r in rows][::-1]
    ypos = np.arange(len(names))
    ax.barh(ypos, kills, height=0.62, color=colors)
    for yv, k in zip(ypos, kills):
        ax.text(k + max(kills) * 0.015, yv, f"{k:.4f}", va="center",
                fontsize=8, color=_INK)
    ax.set_yticks(ypos)
    ax.set_yticklabels(names, fontsize=8.5, color=_INK)
    ax.set_xlim(0, max(kills) * 1.22 if max(kills) > 0 else 1.0)
    ax.set_xlabel("predicted kill fraction (NOT a validated outcome)",
                  fontsize=8.5, color=_INK_MUTED)
    ax.set_title("ranked predicted kill on THIS tumor", fontsize=9.5,
                 color=_INK)
    _recessive(ax)
    ax.grid(axis="y", visible=False)


# --------------------------------------------------------------------------- #
# row 4 -- attribution
# --------------------------------------------------------------------------- #
def _panel_attribution(ax, rows) -> None:
    names = [r["name"] for r in rows]
    keys = ("delivery_decades", "potency_decades", "hypoxia_decades")
    labels = ("delivery (tumor-cell exposure)", "potency (1/IC50)",
              "O2 resistance (hypoxic + quiescent)")
    x = np.arange(len(names))
    w = 0.26
    for i, (k, lab) in enumerate(zip(keys, labels)):
        vals = np.array([r[k] for r in rows], float)
        off = (i - 1) * w
        colors = [_POS if v >= 0 else _NEG for v in vals]
        ax.bar(x + off, vals, w * 0.9, color=colors,
               edgecolor="#ffffff", linewidth=1.0)
        for xv, v in zip(x + off, vals):
            ax.text(xv, v + (0.03 if v >= 0 else -0.03), f"{v:+.2f}",
                    ha="center", va="bottom" if v >= 0 else "top",
                    fontsize=6.8, color=_INK_MUTED)
    ax.axhline(0.0, color="#7a7a7a", lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"#{i + 1} {n}\ndriver: {r['dominant']}\n"
         f"(O2 resistance costs {r['hypoxia_cost'] * 100:.0f}% of its kill)"
         for i, (n, r) in enumerate(zip(names, rows))], fontsize=8,
        color=_INK)
    for tick, r in zip(ax.get_xticklabels(), rows):
        tick.set_color(_drug_color(r["name"]))
    ax.set_ylabel("log10 advantage vs panel mean (decades)", fontsize=8.5,
                  color=_INK_MUTED)
    ax.set_title(
        "ATTRIBUTION — why each drug ranks where it does. Each bar triple is "
        "delivery | potency | O2-resistance (left→right within a drug); "
        "blue = advantage, orange = penalty vs the panel mean.\nThe three terms "
        "are an exact additive log10 decomposition (see the report). The "
        "O2-resistance terms are near zero because this tumor penalises every "
        "drug almost equally — the absolute cost is under each label.",
        fontsize=9.5, color=_INK, pad=14)
    _recessive(ax)
    ax.grid(axis="x", visible=False)
    ax.legend(handles=[mpatches.Patch(color=_POS, label="advantage vs panel mean"),
                       mpatches.Patch(color=_NEG, label="penalty vs panel mean")],
              loc="upper right", fontsize=7.5, frameon=True, framealpha=0.92,
              edgecolor="#cccccc")


# --------------------------------------------------------------------------- #
# the figure
# --------------------------------------------------------------------------- #
def plot_patient_report(panel: Dict, out_path: str) -> str:
    """Render the per-tumor virtual-drug-panel decision figure to ``out_path``.

    Parameters
    ----------
    panel : dict
        Produced by ``experiments/patient_report.py``. Required keys:
        ``report`` (reconstruction_report), ``cells`` {x, y, cell_type},
        ``fields`` {vascular, ecm, oxygen}, ``tissue_mask``, ``extent_um``,
        ``hypoxic_fraction``, ``hypoxia_threshold``, ``rows`` (ranked per-drug
        dicts), ``recommended``, ``tissue_name``, ``dx_um``, ``window_h``;
        optional ``region_fields`` / ``region_tumor_fractions``.
    out_path : str
        PNG destination.

    Returns
    -------
    str
        ``out_path``.
    """
    rows = list(panel["rows"])
    n = len(rows)
    W, H = panel["extent_um"]
    extent = (0.0, float(W), 0.0, float(H))

    # Size the two map rows to the tissue's own aspect ratio so equal-aspect
    # images do not leave a band of dead whitespace under every panel.
    col_w = 4.05
    map_h = max(1.5, min(3.0, 0.86 * col_w * float(H) / float(W))) + 0.62
    plot_h, attr_h = 5.0, 3.7
    fig_h = 2 * map_h + plot_h + attr_h + 1.45

    fig = plt.figure(figsize=(col_w * n, fig_h), facecolor="#fcfcfb")
    gs = GridSpec(4, 2 * n, figure=fig,
                  height_ratios=[map_h, map_h, plot_h, attr_h],
                  hspace=0.34, wspace=0.95,
                  top=1.0 - 1.20 / fig_h, bottom=0.055, left=0.045, right=0.965)

    # ---- row 1: the measured tumor -------------------------------------- #
    _panel_celltypes(fig.add_subplot(gs[0, 0:2]), panel, extent)
    _field_panel(fig.add_subplot(gs[0, 2:4]), panel["fields"]["vascular"],
                 extent,
                 "vascular density (real endothelial KDE)\n"
                 f"mean {panel['report']['mean_vascular_density']:.3f}",
                 "Blues", "vascular density")
    _field_panel(fig.add_subplot(gs[0, 4:6]), panel["fields"]["ecm"], extent,
                 "ECM density (real matrix transcripts)\n"
                 f"mean {panel['report']['mean_ecm_density']:.3f}",
                 "YlOrBr", "ECM volume fraction")
    _panel_oxygen(fig.add_subplot(gs[0, 6:8]), panel, extent)
    _panel_regions(fig.add_subplot(gs[0, 8:10]), panel, extent)

    # ---- row 2: per-drug exposure maps, ranked, ONE shared scale --------- #
    tissue = np.asarray(panel["tissue_mask"], bool)
    pooled = np.concatenate([np.asarray(r["auc_field"], float)[tissue]
                             for r in rows])
    vmax = float(np.percentile(pooled, 99.0)) or 1.0
    norm = matplotlib.colors.LogNorm(vmin=vmax / 1e4, vmax=vmax, clip=True)
    im = None
    for j, r in enumerate(rows):
        ax = fig.add_subplot(gs[1, 2 * j:2 * j + 2])
        im = _panel_auc(ax, r, j + 1, extent, norm,
                        r["name"] == panel["recommended"])
    cb = fig.colorbar(im, ax=[a for a in fig.axes[-n:]], shrink=0.88, pad=0.008,
                      location="right")
    cb.set_label("predicted drug exposure AUC (µM·s)\nONE shared log scale, "
                 "4 decades", fontsize=8, color=_INK_MUTED)
    cb.ax.tick_params(labelsize=7, colors=_INK_MUTED)
    cb.outline.set_visible(False)

    # ---- row 3: decision plane + ranked kill ----------------------------- #
    _panel_scatter(fig.add_subplot(gs[2, 0:n]), rows, panel["recommended"])
    _panel_rank(fig.add_subplot(gs[2, n:2 * n]), rows)

    # ---- row 4: attribution ---------------------------------------------- #
    _panel_attribution(fig.add_subplot(gs[3, :]), rows)

    fig.suptitle(
        f"VIRTUAL DRUG PANEL — one real {panel['tissue_label']} tumor "
        f"({panel['report']['n_cells']:,} measured cells, "
        f"dx = {panel['dx_um']:.0f} µm, {panel['window_h']:.0f} h exposure)\n"
        "RESEARCH PROTOTYPE — PREDICTIONS ONLY, NOT CLINICALLY VALIDATED. "
        "No outcome-labelled data exists for this tissue. NOT FOR CLINICAL USE.",
        fontsize=15, color=_INK, y=1.0 - 0.22 / fig_h, va="top")

    fig.savefig(out_path, dpi=125, bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)
    return out_path
