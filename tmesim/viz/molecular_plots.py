"""Figures for the MOLECULAR layer demo.

Two public entry points that visualize how *different real drugs* — with every
transport parameter derived from molecular structure — penetrate and act on the
SAME patient tumor microenvironment differently:

* :func:`plot_drug_comparison` — per-drug drug-exposure (AUC) heatmaps on a shared
  colour scale (small lipophilics reaching the core vs large/polar drugs and the
  antibody failing to), plus a penetration/kill bar chart and a structure-derived
  parameter table.
* :func:`plot_hindrance_curves` — effective diffusivity vs ECM obstacle fraction
  for each drug (the Brownian-dynamics hindrance curve), exposing size-dependent
  hindrance.

Rendered headless (Agg) and written to disk; nothing is shown interactively.
Every function degrades gracefully if an optional field is missing.
"""
from __future__ import annotations

from typing import List, Sequence

import matplotlib

matplotlib.use("Agg")  # headless: must be set before pyplot is imported

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

__all__ = [
    "plot_drug_comparison",
    "plot_hindrance_curves",
]

# A consistent, colour-blind-friendly per-drug palette (assigned by order).
_DRUG_COLORS = [
    "#1f77b4", "#d62728", "#2ca02c", "#9467bd",
    "#ff7f0e", "#17becf", "#8c564b", "#e377c2",
    "#7f7f7f", "#bcbd22",
]

# Colormap for the AUC (drug-exposure) heatmaps.
_AUC_CMAP = "inferno"


def _color(i: int) -> str:
    return _DRUG_COLORS[i % len(_DRUG_COLORS)]


def _shared_vmax(fields: Sequence[np.ndarray]) -> float:
    """A robust shared upper limit for the AUC colour scale.

    Uses the 99.5th percentile across every drug's exposure field (ignoring a
    lone hot pixel) so the shared scale is dominated by the best-penetrating
    drug's tissue level — which is exactly what makes weakly-penetrating drugs
    and the antibody read as (correctly) dark / unreached.
    """
    vals = []
    for f in fields:
        a = np.asarray(f, dtype=float)
        if a.size:
            vals.append(a[np.isfinite(a)])
    if not vals:
        return 1.0
    allv = np.concatenate(vals)
    vmax = float(np.percentile(allv, 99.5)) if allv.size else 1.0
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(allv.max()) if allv.size and allv.max() > 0 else 1.0
    return vmax


def plot_drug_comparison(results: List[dict], out_path: str) -> None:
    """Multi-panel comparison of how different drugs penetrate the SAME tumor.

    Parameters
    ----------
    results : list of dict
        One entry per drug, each with keys:
          * ``name`` (str)
          * ``profile_summary`` (dict from ``MolecularProfile.summary()``: needs
            ``free_D_um2s``, ``eff_D_um2s``, ``bbb_penetrance``)
          * ``penetration_index`` (float in [0,1])
          * ``kill_fraction`` (float in [0,1])
          * ``auc_field`` ((ny, nx) drug-exposure AUC array, uM*s)
    out_path : str
        Destination PNG path.

    Panels
    ------
    * Top row: each drug's AUC heatmap on ONE shared colour scale (single
      colorbar), so the spatial reach is directly comparable.
    * Bottom-left: grouped bar chart of penetration_index and kill_fraction.
    * Bottom-right: a table of structure-derived free_D / effective_D / BBB.
    """
    results = [r for r in results if r is not None]
    if not results:
        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "no drug results to plot", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        return

    n = len(results)
    names = [r.get("name", f"drug{i}") for i, r in enumerate(results)]
    fields = [np.asarray(r.get("auc_field"), dtype=float) for r in results]
    vmax = _shared_vmax(fields)

    # Layout: top row = n heatmaps (+ a thin colorbar column); bottom row split
    # into a bar-chart panel and a table panel.
    fig = plt.figure(figsize=(max(3.2 * n, 12.0), 9.2))
    gs = GridSpec(
        2, n + 1, figure=fig,
        height_ratios=[1.05, 0.9],
        width_ratios=[1.0] * n + [0.06],
        hspace=0.42, wspace=0.28,
    )

    # --- top row: AUC heatmaps on a shared scale --------------------------
    im = None
    for i, (r, field) in enumerate(zip(results, fields)):
        ax = fig.add_subplot(gs[0, i])
        if field is None or field.size == 0 or not np.isfinite(field).any():
            ax.text(0.5, 0.5, "n/a", ha="center", va="center",
                    transform=ax.transAxes, color="0.5")
            ax.set_xticks([]); ax.set_yticks([])
        else:
            im = ax.imshow(field, origin="lower", cmap=_AUC_CMAP,
                           vmin=0.0, vmax=vmax, aspect="equal")
            ax.set_xticks([]); ax.set_yticks([])
        pen = r.get("penetration_index", float("nan"))
        ax.set_title(f"{names[i]}\npen. index = {pen:.2f}", fontsize=10)

    # single shared colorbar in the thin last column
    if im is not None:
        cax = fig.add_subplot(gs[0, n])
        cbar = fig.colorbar(im, cax=cax)
        cbar.set_label("drug exposure AUC (uM*s)  — shared scale", fontsize=9)

    # --- bottom-left: penetration_index + kill_fraction bar chart ----------
    ax_bar = fig.add_subplot(gs[1, : max(1, (n + 1) // 2)])
    pen_vals = [float(r.get("penetration_index", np.nan)) for r in results]
    kill_vals = [float(r.get("kill_fraction", np.nan)) for r in results]
    x = np.arange(n)
    w = 0.38
    b1 = ax_bar.bar(x - w / 2, pen_vals, w, label="penetration index",
                    color="#4c78a8", edgecolor="0.3")
    b2 = ax_bar.bar(x + w / 2, kill_vals, w, label="kill fraction",
                    color="#e45756", edgecolor="0.3")
    for bars in (b1, b2):
        for rect in bars:
            h = rect.get_height()
            if np.isfinite(h):
                ax_bar.text(rect.get_x() + rect.get_width() / 2, h + 0.01,
                            f"{h:.2f}", ha="center", va="bottom", fontsize=7)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
    ax_bar.set_ylabel("fraction / index [0,1]")
    ax_bar.set_ylim(0, 1.08)
    ax_bar.set_title("Penetration index vs kill fraction", fontsize=11)
    ax_bar.grid(True, axis="y", alpha=0.3)
    ax_bar.legend(fontsize=8, loc="upper right")

    # --- bottom-right: structure-derived parameter table -------------------
    ax_tab = fig.add_subplot(gs[1, max(1, (n + 1) // 2):])
    ax_tab.axis("off")
    col_labels = ["drug", "free_D\n(um^2/s)", "effective_D\n(um^2/s)", "BBB\npenetrance"]
    rows = []
    for r in results:
        ps = r.get("profile_summary", {}) or {}
        free_d = ps.get("free_D_um2s", float("nan"))
        eff_d = ps.get("eff_D_um2s", float("nan"))
        bbb = ps.get("bbb_penetrance", float("nan"))
        rows.append([
            r.get("name", "?"),
            f"{free_d:.1f}" if _finite(free_d) else "n/a",
            f"{eff_d:.1f}" if _finite(eff_d) else "n/a",
            f"{bbb:.3f}" if _finite(bbb) else "n/a",
        ])
    table = ax_tab.table(cellText=rows, colLabels=col_labels,
                         loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.55)
    # colour the header row and tag each drug row with its palette colour
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_facecolor("#31456b")
            cell.set_text_props(color="white", fontweight="bold")
        elif col == 0:
            cell.set_text_props(color=_color(row - 1), fontweight="bold")
    ax_tab.set_title("Structure-derived transport parameters", fontsize=11, pad=12)

    fig.suptitle(
        "Different real drugs, same GBM microenvironment — "
        "structure-derived transport & kill",
        fontsize=15, fontweight="bold",
    )
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_hindrance_curves(profiles: list, out_path: str) -> None:
    """Effective diffusivity vs ECM obstacle fraction, per drug.

    Draws each drug's Brownian-dynamics ``hindrance_curve`` (``{ecm_fraction:
    D_eff}``) as a line, exposing size-dependent hindrance: small molecules keep
    a high effective D as the ECM crowds; large drugs / the antibody collapse.

    Parameters
    ----------
    profiles : list of MolecularProfile
        Each must expose ``.molecule.name`` and ``.hindrance_curve`` (dict) and
        ``.free_diffusivity``.
    out_path : str
        Destination PNG path.
    """
    profiles = [p for p in profiles if p is not None]
    fig, ax = plt.subplots(figsize=(9.0, 6.2))

    if not profiles:
        ax.text(0.5, 0.5, "no profiles to plot", ha="center", va="center",
                transform=ax.transAxes)
        ax.axis("off")
        fig.savefig(out_path, dpi=130, bbox_inches="tight")
        plt.close(fig)
        return

    # order legend by molecule size (MW) so it reads big->small sensibly
    def _mw(p):
        try:
            return float(p.molecule.get("mw", 0.0) or 0.0)
        except Exception:
            return 0.0

    for i, p in enumerate(sorted(profiles, key=_mw)):
        curve = getattr(p, "hindrance_curve", {}) or {}
        if not curve:
            continue
        fracs = sorted(curve.keys())
        d_eff = [curve[f] for f in fracs]
        name = getattr(getattr(p, "molecule", None), "name", f"drug{i}")
        mw = _mw(p)
        ax.plot(fracs, d_eff, "-o", color=_color(i), lw=2.0, markersize=6,
                label=f"{name} (MW {mw:.0f})")

    ax.set_xlabel("ECM / obstacle volume fraction  φ", fontsize=11)
    ax.set_ylabel("effective diffusivity  D_eff (um^2/s)", fontsize=11)
    ax.set_title(
        "Size-dependent ECM hindrance (Brownian dynamics)\n"
        "larger drugs lose more effective diffusivity as ECM crowds",
        fontsize=12, fontweight="bold",
    )
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    ax.legend(fontsize=9, loc="upper right", title="drug (structure-derived)")

    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _finite(x) -> bool:
    try:
        return isinstance(x, (int, float)) and np.isfinite(x)
    except Exception:
        return False
