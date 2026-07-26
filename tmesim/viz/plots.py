"""Figure generation for tmesim.

Three public entry points, one per stage of the pipeline:

* :func:`plot_microenvironment` — the reconstructed TME (cells + tissue fields).
* :func:`plot_simulation` — the outcome of one drug-response simulation.
* :func:`plot_validation` — the honest Stage-0 verdict (sim vs. baselines).

All figures are rendered with the non-interactive ``Agg`` backend and written to
disk; nothing is ever shown interactively. Every function tolerates missing
optional fields (oxygen, ifp, per-patient scores, ...) and simply omits the
panels it cannot draw.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless: must be set before pyplot is imported

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from tmesim.core import (
    CellType,
    ExperimentResult,
    Microenvironment,
    SimulationResult,
)

__all__ = [
    "plot_microenvironment",
    "plot_simulation",
    "plot_validation",
]

# Consistent per-CellType colors used across every figure.
_CELLTYPE_COLORS = {
    CellType.TUMOR: "#d62728",        # red
    CellType.STROMAL: "#2ca02c",      # green
    CellType.IMMUNE: "#1f77b4",       # blue
    CellType.ENDOTHELIAL: "#9467bd",  # purple
    CellType.NECROTIC: "#7f7f7f",     # grey
}


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _imshow_field(ax: plt.Axes, field: np.ndarray, extent: tuple,
                  title: str, cmap: str, label: str) -> None:
    """Draw a (ny, nx) field as a heatmap with an attached colorbar.

    ``extent`` is (0, width_um, 0, height_um). ``origin='lower'`` keeps the
    physical y-axis pointing up while matching the numpy (row=y) convention.
    """
    im = ax.imshow(np.asarray(field), origin="lower", extent=extent,
                   cmap=cmap, aspect="equal")
    ax.set_title(title)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(label)


def _blank_panel(ax: plt.Axes, message: str) -> None:
    """Turn an axis into a centered 'not available' placeholder."""
    ax.text(0.5, 0.5, message, ha="center", va="center",
            transform=ax.transAxes, fontsize=10, color="0.4", wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])


def _save(fig: Figure, out_path: str) -> None:
    """Tight-layout and write a figure to ``out_path`` as PNG, then close it."""
    fig.tight_layout()
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# 1. microenvironment
# --------------------------------------------------------------------------- #
def plot_microenvironment(micro: Microenvironment, out_path: str) -> None:
    """Render a multi-panel view of a reconstructed microenvironment.

    Panels
    ------
    * Cell scatter in micrometers, colored by :class:`CellType` (legend).
    * Vascular-density heatmap (drug/oxygen source).
    * ECM-density heatmap (tortuosity / diffusion barrier).
    * Oxygen and IFP heatmaps, drawn only when those fields are populated.

    Parameters
    ----------
    micro : Microenvironment
        The reconstructed TME to visualize.
    out_path : str
        Destination PNG path.
    """
    domain = micro.domain
    fields = micro.fields
    cells = micro.cells
    width_um, height_um = domain.extent_um
    extent = (0.0, width_um, 0.0, height_um)

    # Decide which optional fields exist so we can size the grid.
    optional = []
    if getattr(fields, "oxygen", None) is not None:
        optional.append(("oxygen", fields.oxygen, "Oxygen (relative)", "viridis"))
    if getattr(fields, "ifp", None) is not None:
        optional.append(("ifp", fields.ifp, "Interstitial fluid pressure (mmHg)", "magma"))

    n_panels = 3 + len(optional)  # cells + vascular + ecm + optionals
    ncols = 3
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.6 * nrows))
    axes = np.atleast_1d(axes).ravel()

    # -- panel 0: cells colored by type -----------------------------------
    ax = axes[0]
    ctype = np.asarray(cells.cell_type)
    for ct in CellType:
        sel = ctype == int(ct)
        if not np.any(sel):
            continue
        ax.scatter(cells.x[sel], cells.y[sel], s=6, alpha=0.7,
                   color=_CELLTYPE_COLORS.get(ct, "#000000"), label=ct.name.title(),
                   edgecolors="none")
    ax.set_title(f"Cells by type (n={cells.n})")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_xlim(0, width_um)
    ax.set_ylim(0, height_um)
    ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=7, framealpha=0.85, markerscale=1.5)

    # -- panel 1: vascular density ----------------------------------------
    _imshow_field(axes[1], fields.vascular_density, extent,
                  "Vascular density", "Reds", "fraction perfused")

    # -- panel 2: ECM density ---------------------------------------------
    _imshow_field(axes[2], fields.ecm_density, extent,
                  "ECM density", "YlGn", "relative density")

    # -- optional panels ---------------------------------------------------
    for i, (_key, field, title, cmap) in enumerate(optional, start=3):
        _imshow_field(axes[i], field, extent, title, cmap, title.split("(")[0].strip())

    # blank out any unused axes
    for j in range(n_panels, len(axes)):
        axes[j].axis("off")

    fig.suptitle(f"Microenvironment — {micro.patient_id}", fontsize=14,
                 fontweight="bold")
    _save(fig, out_path)


# --------------------------------------------------------------------------- #
# 2. simulation
# --------------------------------------------------------------------------- #
def plot_simulation(result: SimulationResult, out_path: str) -> None:
    """Render the outcome of one drug-response simulation.

    Panels
    ------
    * Any spatial fields present in ``result.fields`` (oxygen, ifp, auc /
      drug penetration, vascular, ...) as heatmaps.
    * The viable-tumor-fraction-vs-time curve, annotated with the predicted
      response and kill fraction.

    Missing fields are skipped gracefully. Time is shown in hours.

    Parameters
    ----------
    result : SimulationResult
        The simulation output to visualize.
    out_path : str
        Destination PNG path.
    """
    # Preferred display order + nice labels/colormaps for known field names.
    field_specs = [
        ("oxygen", "Oxygen (relative)", "viridis", "O2"),
        ("ifp", "Interstitial fluid pressure (mmHg)", "magma", "mmHg"),
        ("auc", "Drug exposure AUC (uM*s)", "inferno", "AUC"),
        ("drug_penetration", "Drug penetration", "inferno", "penetration"),
        ("vascular", "Vascular density", "Reds", "fraction"),
    ]
    fields = dict(result.fields or {})
    ordered = [(k, t, c, lbl) for (k, t, c, lbl) in field_specs if k in fields]
    # include any extra fields not in our known list
    known = {k for k, *_ in field_specs}
    for k, v in fields.items():
        if k not in known and v is not None:
            ordered.append((k, k.replace("_", " ").title(), "viridis", k))

    n_panels = len(ordered) + 1  # + viable curve
    ncols = min(3, max(2, n_panels))
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.6 * nrows))
    axes = np.atleast_1d(axes).ravel()

    idx = 0
    for key, title, cmap, label in ordered:
        field = fields.get(key)
        if field is None:
            _blank_panel(axes[idx], f"{title}\n(not available)")
            axes[idx].set_title(title)
        else:
            # extent unknown here (no domain on result); use array index extent.
            arr = np.asarray(field)
            im = axes[idx].imshow(arr, origin="lower", cmap=cmap, aspect="equal")
            axes[idx].set_title(title)
            axes[idx].set_xlabel("x (grid)")
            axes[idx].set_ylabel("y (grid)")
            cbar = fig.colorbar(im, ax=axes[idx], fraction=0.046, pad=0.04)
            cbar.set_label(label)
        idx += 1

    # -- viable tumor fraction vs time ------------------------------------
    ax = axes[idx]
    hist = np.asarray(result.viable_history) if result.viable_history is not None else np.array([])
    times = np.asarray(result.times) if result.times is not None else np.array([])
    if hist.size and times.size == hist.size:
        t_hours = times / 3600.0
        ax.plot(t_hours, hist, "-o", color="#d62728", markersize=3, lw=1.8)
        ax.set_xlabel("time (hours)")
    elif hist.size:
        ax.plot(np.arange(hist.size), hist, "-o", color="#d62728", markersize=3, lw=1.8)
        ax.set_xlabel("step")
    else:
        _blank_panel(ax, "no viable-fraction history")
    ax.set_ylabel("viable tumor fraction")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Viable tumor fraction over treatment")
    ax.grid(True, alpha=0.3)

    # annotation box with headline numbers
    txt = (f"predicted response: {result.predicted_response:.3f}\n"
           f"kill fraction: {result.kill_fraction:.3f}\n"
           f"final viable: {result.viable_tumor_fraction:.3f}")
    ax.text(0.97, 0.95, txt, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, bbox=dict(boxstyle="round", facecolor="#fff4e6",
                                  edgecolor="0.6", alpha=0.95))
    idx += 1

    for j in range(idx, len(axes)):
        axes[j].axis("off")

    fig.suptitle("Simulation result", fontsize=14, fontweight="bold")
    _save(fig, out_path)


# --------------------------------------------------------------------------- #
# 3. validation
# --------------------------------------------------------------------------- #
def _extract_true_score(per_patient: list) -> tuple:
    """Best-effort pull of (y_true, sim_score) arrays from per_patient dicts.

    Looks for common key spellings; returns (None, None) if unavailable.
    """
    if not per_patient:
        return None, None
    true_keys = ("observed_response", "observed", "y_true", "label", "truth")
    score_keys = ("predicted_response", "sim", "sim_score", "predicted", "y_score", "score")
    y_true, y_score = [], []
    for row in per_patient:
        if not isinstance(row, dict):
            return None, None
        tk = next((k for k in true_keys if k in row), None)
        sk = next((k for k in score_keys if k in row), None)
        if tk is None or sk is None:
            return None, None
        y_true.append(row[tk])
        y_score.append(row[sk])
    return np.asarray(y_true, dtype=float), np.asarray(y_score, dtype=float)


def _roc_points(y_true: np.ndarray, y_score: np.ndarray) -> tuple:
    """Compute ROC (fpr, tpr) for a binary label with a simple threshold sweep."""
    y_true = np.asarray(y_true)
    # binarize continuous labels at their median so a curve is still drawable
    if not np.all(np.isin(np.unique(y_true), (0, 1))):
        y_true = (y_true >= np.median(y_true)).astype(int)
    pos = y_true == 1
    n_pos, n_neg = int(pos.sum()), int((~pos).sum())
    if n_pos == 0 or n_neg == 0:
        return None, None
    order = np.argsort(-y_score)
    ys = y_true[order]
    tpr = np.concatenate([[0.0], np.cumsum(ys == 1) / n_pos])
    fpr = np.concatenate([[0.0], np.cumsum(ys == 0) / n_neg])
    return fpr, tpr


def _calibration_points(y_true: np.ndarray, y_score: np.ndarray, bins: int = 8) -> tuple:
    """Return (mean_pred, mean_true) per score bin for a reliability diagram."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    edges = np.linspace(0.0, 1.0, bins + 1)
    idx = np.clip(np.digitize(y_score, edges[1:-1]), 0, bins - 1)
    mp, mt = [], []
    for b in range(bins):
        sel = idx == b
        if np.any(sel):
            mp.append(float(y_score[sel].mean()))
            mt.append(float(y_true[sel].mean()))
    return np.asarray(mp), np.asarray(mt)


def plot_validation(exp: ExperimentResult, out_path: str) -> None:
    """Render the Stage-0 verdict: mechanistic sim vs. the cheap baselines.

    Panels
    ------
    * Bar chart of the primary metric for the sim model and every baseline,
      with the sim bar highlighted and a value on each bar.
    * ROC curve for the sim (and baselines when per-patient scores are present),
      or a placeholder when raw scores are unavailable.
    * Calibration / reliability diagram for the sim predictions.

    The suptitle states plainly whether the sim beat the baselines.

    Parameters
    ----------
    exp : ExperimentResult
        The experiment output to visualize.
    out_path : str
        Destination PNG path.
    """
    primary = exp.primary_metric or "roc_auc"
    model_scores = dict(exp.model_scores or {})
    baseline_scores = dict(exp.baseline_scores or {})

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
    ax_bar, ax_roc, ax_cal = axes

    # -- panel A: primary-metric bar chart --------------------------------
    names = ["sim"] + list(baseline_scores.keys())
    values = [model_scores.get(primary, np.nan)]
    for bname in baseline_scores:
        values.append(baseline_scores[bname].get(primary, np.nan))
    colors = ["#d62728"] + ["#7f9cb5"] * len(baseline_scores)
    bars = ax_bar.bar(range(len(names)), values, color=colors, edgecolor="0.3")
    ax_bar.set_xticks(range(len(names)))
    ax_bar.set_xticklabels(names, rotation=25, ha="right", fontsize=8)
    ax_bar.set_ylabel(primary)
    ax_bar.set_title(f"Primary metric: {primary}")
    ax_bar.grid(True, axis="y", alpha=0.3)
    for rect, val in zip(bars, values):
        if np.isfinite(val):
            ax_bar.text(rect.get_x() + rect.get_width() / 2, rect.get_height(),
                        f"{val:.3f}", ha="center", va="bottom", fontsize=8)
    # reference line at best baseline
    base_vals = [v for v in values[1:] if np.isfinite(v)]
    if base_vals:
        ax_bar.axhline(max(base_vals), color="0.4", ls="--", lw=1,
                       label="best baseline")
        ax_bar.legend(fontsize=8, loc="lower right")

    # -- panel B: ROC ------------------------------------------------------
    y_true, y_score = _extract_true_score(exp.per_patient)
    drew_roc = False
    if y_true is not None and y_score is not None and y_true.size >= 2:
        fpr, tpr = _roc_points(y_true, y_score)
        if fpr is not None:
            _trapz = getattr(np, "trapezoid", np.trapz)
            auc = float(model_scores.get("roc_auc", _trapz(tpr, fpr)))
            ax_roc.plot(fpr, tpr, "-", color="#d62728", lw=2,
                        label=f"sim (AUC={auc:.3f})")
            drew_roc = True
    if drew_roc:
        ax_roc.plot([0, 1], [0, 1], "--", color="0.6", lw=1, label="chance")
        ax_roc.set_xlabel("false positive rate")
        ax_roc.set_ylabel("true positive rate")
        ax_roc.set_xlim(0, 1)
        ax_roc.set_ylim(0, 1)
        ax_roc.set_title("ROC — sim")
        ax_roc.legend(fontsize=8, loc="lower right")
    else:
        _blank_panel(ax_roc, "per-patient scores unavailable\n(ROC from metric only)")
        ax_roc.set_title("ROC")

    # -- panel C: calibration ---------------------------------------------
    drew_cal = False
    if y_true is not None and y_score is not None and y_true.size >= 2:
        # only meaningful if scores look like probabilities in [0,1]
        if np.nanmin(y_score) >= -0.01 and np.nanmax(y_score) <= 1.01:
            yt = y_true.copy()
            if not np.all(np.isin(np.unique(yt), (0, 1))):
                yt = (yt >= np.median(yt)).astype(float)
            mp, mt = _calibration_points(yt, y_score)
            if mp.size:
                ax_cal.plot(mp, mt, "-o", color="#1f77b4", lw=1.8, markersize=5,
                            label="sim")
                drew_cal = True
    if drew_cal:
        ax_cal.plot([0, 1], [0, 1], "--", color="0.6", lw=1, label="perfect")
        ax_cal.set_xlabel("mean predicted response")
        ax_cal.set_ylabel("observed responder fraction")
        ax_cal.set_xlim(0, 1)
        ax_cal.set_ylim(0, 1)
        ax_cal.set_title("Calibration")
        ax_cal.legend(fontsize=8, loc="upper left")
    else:
        _blank_panel(ax_cal, "calibration unavailable\n(need per-patient probabilities)")
        ax_cal.set_title("Calibration")

    # -- headline verdict --------------------------------------------------
    if exp.beat_baseline:
        banner = "SIM BEAT THE BASELINES"
        color = "#1a7f37"
    else:
        banner = "SIM DID NOT BEAT THE BASELINES"
        color = "#b42318"
    verdict = exp.verdict or ""
    fig.suptitle(f"Stage-0 validation — {banner}\n{verdict}",
                 fontsize=13, fontweight="bold", color=color)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# self-test
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    import os

    from tmesim.core import (
        CellPopulation,
        CellState,
        Domain,
        TissueFields,
    )

    rng = np.random.default_rng(0)
    out_dir = "/tmp"

    # --- fabricate a small Microenvironment ------------------------------
    ny, nx, dx = 40, 40, 15.0
    domain = Domain(nx=nx, ny=ny, dx=dx)
    Xg, Yg = domain.coords_um()
    cx, cy = nx * dx / 2, ny * dx / 2
    r = np.hypot(Xg - cx, Yg - cy)
    # vascular: sparser in the core
    vascular = np.clip(0.15 + 0.5 * (r / r.max()) + 0.1 * rng.standard_normal((ny, nx)), 0, 1)
    ecm = np.clip(0.3 + 0.4 * np.exp(-((r / (0.4 * r.max())) ** 2)) + 0.05 * rng.standard_normal((ny, nx)), 0, 1)
    oxygen = np.clip(vascular * 0.9 + 0.1, 0, 1)
    ifp = np.clip(12.0 * np.exp(-((r / (0.5 * r.max())) ** 2)), 0, None)
    fields = TissueFields(ecm_density=ecm, vascular_density=vascular,
                          oxygen=oxygen, ifp=ifp)

    n_cells = 400
    ang = rng.uniform(0, 2 * np.pi, n_cells)
    rad = rng.uniform(0, min(cx, cy), n_cells)
    xs = np.clip(cx + rad * np.cos(ang), 0, nx * dx)
    ys = np.clip(cy + rad * np.sin(ang), 0, ny * dx)
    ctype = rng.choice([0, 1, 2, 3, 4], size=n_cells, p=[0.55, 0.2, 0.1, 0.1, 0.05])
    state = np.full(n_cells, int(CellState.VIABLE), dtype=int)
    sens = rng.uniform(0, 1, n_cells)
    cells = CellPopulation(x=xs, y=ys, cell_type=ctype, state=state, sensitivity=sens)
    micro = Microenvironment(domain=domain, fields=fields, cells=cells,
                             patient_id="selftest-001",
                             metadata={"subtype": "mesenchymal"})

    # --- fabricate a SimulationResult ------------------------------------
    times = np.linspace(0, 72 * 3600, 12)
    viable_hist = np.linspace(1.0, 0.35, 12) + 0.02 * rng.standard_normal(12)
    sim = SimulationResult(
        predicted_response=0.62,
        viable_tumor_fraction=float(viable_hist[-1]),
        kill_fraction=float(1 - viable_hist[-1]),
        viable_history=viable_hist,
        times=times,
        transport=None,
        fields={"oxygen": oxygen, "ifp": ifp,
                "auc": np.abs(rng.standard_normal((ny, nx))) * vascular,
                "vascular": vascular},
        extras={},
    )

    # --- fabricate an ExperimentResult -----------------------------------
    n_pat = 30
    obs = rng.integers(0, 2, n_pat).astype(float)
    pred = np.clip(0.5 + 0.25 * (obs - 0.5) * 2 + 0.2 * rng.standard_normal(n_pat), 0, 1)
    per_patient = [{"patient_id": f"p{i}", "observed_response": float(obs[i]),
                    "predicted_response": float(pred[i])} for i in range(n_pat)]
    exp = ExperimentResult(
        model_scores={"roc_auc": 0.78, "brier": 0.17},
        baseline_scores={
            "mean_prior": {"roc_auc": 0.50, "brier": 0.25},
            "intrinsic_sensitivity": {"roc_auc": 0.66, "brier": 0.21},
            "expression_logreg": {"roc_auc": 0.71, "brier": 0.19},
        },
        per_patient=per_patient,
        primary_metric="roc_auc",
        verdict="sim roc_auc 0.78 > best baseline 0.71 (+0.07)",
        beat_baseline=True,
    )

    paths = {
        "micro": os.path.join(out_dir, "selftest_microenvironment.png"),
        "sim": os.path.join(out_dir, "selftest_simulation.png"),
        "val": os.path.join(out_dir, "selftest_validation.png"),
    }
    plot_microenvironment(micro, paths["micro"])
    plot_simulation(sim, paths["sim"])
    plot_validation(exp, paths["val"])

    for name, p in paths.items():
        assert os.path.exists(p), f"{name} figure not written: {p}"
        size = os.path.getsize(p)
        assert size > 0, f"{name} figure is empty: {p}"
        print(f"OK  {name:6s} -> {p}  ({size} bytes)")

    print("self-test passed")
