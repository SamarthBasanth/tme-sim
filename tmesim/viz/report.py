"""Stage-0 report generator for tmesim.

This module is pure *formatting and orchestration*: it runs the three figure
functions from :mod:`tmesim.viz.plots` and stitches their output together with
the numbers already computed by the validation harness into a single, honest
Markdown report. It introduces **no new science** — every metric, verdict, and
prediction it prints comes verbatim from the :class:`ExperimentResult` and
:class:`SimulationResult` it is handed.

Public entry point:

* :func:`build_report` — write ``<out_dir>/report.md`` plus the three PNG
  figures it embeds, and return the path to the Markdown file.
"""
from __future__ import annotations

import os
from typing import List, Optional

import numpy as np

from tmesim.core import ExperimentResult, Microenvironment, SimulationResult
from tmesim.viz.plots import (
    plot_microenvironment,
    plot_simulation,
    plot_validation,
)

__all__ = ["build_report"]

# Metrics shown, in order, in the comparison table. Each is (key, direction)
# where direction is "up" if higher is better and "down" if lower is better.
_TABLE_METRICS = [
    ("roc_auc", "up"),
    ("brier", "down"),
    ("spearman", "up"),
]


# --------------------------------------------------------------------------- #
# small formatting helpers
# --------------------------------------------------------------------------- #
def _fmt(value: Optional[float]) -> str:
    """Format a metric value for a table cell, or an em dash when absent."""
    if value is None:
        return "—"
    try:
        v = float(value)
    except (TypeError, ValueError):
        return "—"
    if not np.isfinite(v):
        return "—"
    return f"{v:.3f}"


def _best_baseline_name(baseline_scores: dict, metric: str, direction: str) -> Optional[str]:
    """Return the baseline name with the best value on ``metric`` (or None)."""
    best_name: Optional[str] = None
    best_val: Optional[float] = None
    for name, scores in baseline_scores.items():
        val = scores.get(metric) if isinstance(scores, dict) else None
        if val is None:
            continue
        try:
            val = float(val)
        except (TypeError, ValueError):
            continue
        if not np.isfinite(val):
            continue
        if best_val is None:
            best_name, best_val = name, val
            continue
        if direction == "up" and val > best_val:
            best_name, best_val = name, val
        elif direction == "down" and val < best_val:
            best_name, best_val = name, val
    return best_name


def _metrics_table(exp: ExperimentResult) -> str:
    """Build the sim-vs-baselines Markdown table.

    One row per model (sim first, then every baseline), one column per metric in
    :data:`_TABLE_METRICS`. The best *baseline* value in each metric column is
    bolded so the reader can see at a glance what the sim had to beat. The
    primary metric column header is annotated as such.
    """
    model_scores = dict(exp.model_scores or {})
    baseline_scores = dict(exp.baseline_scores or {})
    primary = exp.primary_metric or "roc_auc"

    # header
    headers = ["model"]
    for key, direction in _TABLE_METRICS:
        arrow = "higher better" if direction == "up" else "lower better"
        tag = ", primary" if key == primary else ""
        headers.append(f"{key} ({arrow}{tag})")
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")

    # precompute which baseline is best per metric (for bolding)
    best_per_metric = {
        key: _best_baseline_name(baseline_scores, key, direction)
        for key, direction in _TABLE_METRICS
    }

    # sim row
    sim_cells = ["**sim (mechanistic)**"]
    for key, _direction in _TABLE_METRICS:
        sim_cells.append(_fmt(model_scores.get(key)))
    lines.append("| " + " | ".join(sim_cells) + " |")

    # baseline rows
    for name, scores in baseline_scores.items():
        cells = [f"`{name}`"]
        for key, _direction in _TABLE_METRICS:
            val = scores.get(key) if isinstance(scores, dict) else None
            text = _fmt(val)
            if best_per_metric.get(key) == name and text != "—":
                text = f"**{text}**"
            cells.append(text)
        lines.append("| " + " | ".join(cells) + " |")

    footnote = (
        "\n*Best baseline per metric is shown in **bold**. The Stage-0 bar is that "
        "the mechanistic sim must beat the best baseline on the primary metric "
        f"(`{primary}`).*"
    )
    return "\n".join(lines) + "\n" + footnote


def _verdict_section(exp: ExperimentResult) -> str:
    """Build the headline verdict block with an explicit PASS/FAIL banner."""
    if exp.beat_baseline:
        banner = "> ## ✅ PASS — the mechanistic sim beat the baselines"
    else:
        banner = "> ## ❌ FAIL — the mechanistic sim did NOT beat the baselines"
    verdict = (exp.verdict or "").strip() or "(no verdict string provided)"
    # Present the verdict string verbatim, inside the same blockquote banner.
    return f"{banner}\n>\n> **Verdict:** {verdict}\n"


# --------------------------------------------------------------------------- #
# public entry point
# --------------------------------------------------------------------------- #
def build_report(
    exp: ExperimentResult,
    sim_example: SimulationResult,
    micro_example: Microenvironment,
    out_dir: str = "results",
) -> str:
    """Generate the Stage-0 Markdown report and its three figures.

    Renders (into ``out_dir``) ``microenvironment.png``, ``simulation_example.png``,
    and ``validation.png`` via the :mod:`tmesim.viz.plots` functions, then writes
    ``report.md`` embedding them (by relative path) around the honest verdict and
    metrics already contained in ``exp``.

    The report presents, in order:

    1. the headline verdict (verbatim ``exp.verdict``) with a PASS/FAIL banner
       driven by ``exp.beat_baseline``;
    2. a metrics table (sim vs. every baseline on roc_auc / brier / spearman)
       with the best baseline highlighted;
    3. an example reconstructed microenvironment figure + caption;
    4. an example simulation figure (drug penetration / oxygen / IFP / viable
       curve) + caption;
    5. the validation figure (ROC / calibration / metric bars) + caption;
    6. an honest interpretation section, including the Stage-0 scope caveat.

    Parameters
    ----------
    exp : ExperimentResult
        The scored Stage-0 experiment (source of the verdict, metrics, banner).
    sim_example : SimulationResult
        One representative simulation to illustrate the mechanistic output.
    micro_example : Microenvironment
        One representative reconstructed microenvironment to illustrate the input.
    out_dir : str, optional
        Directory to create (if missing) and write the report + figures into.

    Returns
    -------
    str
        Absolute path to the written ``report.md``.
    """
    os.makedirs(out_dir, exist_ok=True)

    # --- 1. render the three figures -------------------------------------
    micro_png = "microenvironment.png"
    sim_png = "simulation_example.png"
    val_png = "validation.png"
    plot_microenvironment(micro_example, os.path.join(out_dir, micro_png))
    plot_simulation(sim_example, os.path.join(out_dir, sim_png))
    plot_validation(exp, os.path.join(out_dir, val_png))

    # --- 2. assemble the markdown ----------------------------------------
    primary = exp.primary_metric or "roc_auc"
    n_patients = len(exp.per_patient or [])

    parts: List[str] = []

    parts.append("# tmesim — Stage 0 kill-test report\n")
    parts.append(
        "*Does a mechanistic tumor-microenvironment + drug-transport simulation "
        "predict observed glioblastoma drug response better than a cheap, "
        "non-spatial baseline?* This report answers that question for the cohort "
        f"that was run ({n_patients} patient-sample(s)), and shows the "
        "reconstruction, one simulation, and the head-to-head scoring behind the "
        "verdict.\n"
    )

    # (1) headline verdict
    parts.append("## 1. Headline verdict\n")
    parts.append(_verdict_section(exp))

    # (2) metrics table
    parts.append("\n## 2. Metrics — sim vs. baselines\n")
    parts.append(
        f"Primary metric: **`{primary}`**. Scores below are cohort-level "
        "(cross-validated where a baseline needs training). A blank cell "
        "means that metric was not computed for that model.\n"
    )
    parts.append(_metrics_table(exp))

    # (3) microenvironment figure
    parts.append("\n## 3. Example reconstructed microenvironment\n")
    parts.append(f"![Reconstructed microenvironment]({micro_png})\n")
    parts.append(
        f"**Figure 1.** Reconstructed microenvironment for `{micro_example.patient_id}` "
        f"({micro_example.cells.n} cells on a {micro_example.domain.shape[0]}×"
        f"{micro_example.domain.shape[1]} grid, dx = {micro_example.domain.dx:g} µm). "
        "Panels: cells colored by phenotype (tumor / stromal / immune / endothelial "
        "/ necrotic); vascular density (the drug and oxygen source, sparser in the "
        "tumor core); ECM density (raises tortuosity and lowers effective drug "
        "diffusivity); and, when present, the derived oxygen and interstitial-fluid-"
        "pressure fields. This is the spatial input the mechanistic model uses and "
        "the non-spatial baselines throw away.\n"
    )

    # (4) simulation figure
    parts.append("\n## 4. Example simulation\n")
    parts.append(f"![Example simulation]({sim_png})\n")
    parts.append(
        "**Figure 2.** One drug-response simulation. Spatial panels show oxygen, "
        "interstitial fluid pressure, and the drug-exposure AUC (drug penetration) "
        "over the treatment course, alongside vascular density. The final panel is "
        "the viable-tumor-fraction-vs-time curve. Headline numbers for this run: "
        f"predicted response = {sim_example.predicted_response:.3f}, kill fraction = "
        f"{sim_example.kill_fraction:.3f}, final viable tumor fraction = "
        f"{sim_example.viable_tumor_fraction:.3f}. Poor penetration into the "
        "hypoxic, high-IFP core is exactly the effect a non-spatial model cannot "
        "represent.\n"
    )

    # (5) validation figure
    parts.append("\n## 5. Validation\n")
    parts.append(f"![Stage-0 validation]({val_png})\n")
    parts.append(
        "**Figure 3.** Stage-0 head-to-head. Left: the primary metric for the sim "
        "and every baseline, with the best-baseline reference line — the bar the "
        "sim must clear. Middle: the ROC curve for the sim predictions (drawn when "
        "per-patient scores are available). Right: the calibration / reliability "
        "diagram, i.e. whether predicted responses match observed responder "
        "fractions. The figure's title states plainly whether the sim beat the "
        "baselines.\n"
    )

    # (6) honest interpretation
    parts.append("\n## 6. Interpretation (read this before believing anything)\n")
    if exp.beat_baseline:
        parts.append(
            f"On this cohort the mechanistic sim **beat** the best baseline on the "
            f"primary metric (`{primary}`). What this *does* suggest: the spatial, "
            "transport-aware structure carried signal about drug response that the "
            "cheap non-spatial baselines did not capture, at least on these data. "
            "What it **does not** establish: clinical utility, generalization to "
            "unseen or real patients, or that the mechanism is correct rather than "
            "merely correlated — a favorable synthetic result can reflect the "
            "generator's own assumptions.\n"
        )
    else:
        parts.append(
            f"On this cohort the mechanistic sim **did not beat** the best baseline "
            f"on the primary metric (`{primary}`). Per the Stage-0 contract this is "
            "the honest, non-flattering outcome: the added mechanistic machinery did "
            "not earn its complexity here. That does **not** prove the mechanism is "
            "wrong — it may reflect the data, the cohort size, or unmodeled biology "
            "— but it does mean the idea is **not ready** on this evidence and "
            "should not be oversold.\n"
        )
    parts.append(
        "\n### Stage-0 scope caveat\n"
        "This is a **Stage-0 kill-test on synthetic / best-effort data**, not a "
        "clinical validation. The microenvironments and observed responses are "
        "synthetic (or best-effort open-data proxies); labels are a *noisy* "
        "function of the modeled biology, chosen so the test is fair rather than "
        "trivially solvable. Nothing here has been checked against real patient "
        "outcomes, no regulatory or clinical claim is implied, and the numbers "
        "above should be read as evidence about **whether the machinery is worth "
        "building further**, not as a measure of real-world predictive accuracy. A "
        "genuine validation would require controlled-access spatial and drug-screen "
        "data and prospective outcomes.\n"
    )

    report_md = "\n".join(parts)

    # --- 3. write and return ---------------------------------------------
    report_path = os.path.join(out_dir, "report.md")
    with open(report_path, "w", encoding="utf-8") as fh:
        fh.write(report_md)

    return os.path.abspath(report_path)


# --------------------------------------------------------------------------- #
# self-test
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    from tmesim.core import (
        CellPopulation,
        CellState,
        Domain,
        TissueFields,
    )

    rng = np.random.default_rng(0)
    test_out = "/tmp/tmesim_report_test"

    # --- fabricate a small Microenvironment ------------------------------
    ny, nx, dx = 40, 40, 15.0
    domain = Domain(nx=nx, ny=ny, dx=dx)
    Xg, Yg = domain.coords_um()
    cx, cy = nx * dx / 2, ny * dx / 2
    r = np.hypot(Xg - cx, Yg - cy)
    vascular = np.clip(0.15 + 0.5 * (r / r.max()) + 0.1 * rng.standard_normal((ny, nx)), 0, 1)
    ecm = np.clip(0.3 + 0.4 * np.exp(-((r / (0.4 * r.max())) ** 2)), 0, 1)
    oxygen = np.clip(vascular * 0.9 + 0.1, 0, 1)
    ifp = np.clip(12.0 * np.exp(-((r / (0.5 * r.max())) ** 2)), 0, None)
    fields = TissueFields(ecm_density=ecm, vascular_density=vascular,
                          oxygen=oxygen, ifp=ifp)

    n_cells = 300
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
    viable_hist = np.linspace(1.0, 0.35, 12)
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
    verdict_text = "sim roc_auc 0.780 > best baseline 0.710 (+0.070) — PASS"
    exp = ExperimentResult(
        model_scores={"roc_auc": 0.78, "brier": 0.17, "spearman": 0.55},
        baseline_scores={
            "mean_prior": {"roc_auc": 0.50, "brier": 0.25, "spearman": 0.00},
            "intrinsic_sensitivity": {"roc_auc": 0.66, "brier": 0.21, "spearman": 0.41},
            "expression_logreg": {"roc_auc": 0.71, "brier": 0.19, "spearman": 0.48},
        },
        per_patient=per_patient,
        primary_metric="roc_auc",
        verdict=verdict_text,
        beat_baseline=True,
    )

    # --- build the report ------------------------------------------------
    path = build_report(exp, sim, micro, out_dir=test_out)

    # --- assertions ------------------------------------------------------
    assert os.path.exists(path), f"report.md not written: {path}"
    for png in ("microenvironment.png", "simulation_example.png", "validation.png"):
        p = os.path.join(test_out, png)
        assert os.path.exists(p), f"figure not written: {p}"
        assert os.path.getsize(p) > 0, f"figure is empty: {p}"
        print(f"OK  figure -> {p}  ({os.path.getsize(p)} bytes)")

    with open(path, encoding="utf-8") as fh:
        content = fh.read()
    assert verdict_text in content, "verdict text not found verbatim in report.md"
    assert "PASS" in content, "PASS banner missing"
    assert "Stage-0 scope caveat" in content, "scope caveat missing"
    assert "microenvironment.png" in content, "micro figure not embedded"
    assert "simulation_example.png" in content, "sim figure not embedded"
    assert "validation.png" in content, "validation figure not embedded"
    print(f"OK  report -> {path}  ({len(content)} chars)")
    print("self-test passed")
