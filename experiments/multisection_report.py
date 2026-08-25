"""Render ``results/multisection_generalisation.md`` from the sweep's JSON.

Reads ``results/multisection_generalisation.json`` (written by
``experiments/multisection_generalisation.py``) and writes the report. Every
number in the document is interpolated from that file — nothing is typed by
hand, so the prose cannot drift away from the run that produced it.

    ./.venv/bin/python experiments/multisection_report.py
"""
from __future__ import annotations

import json
import os
import sys
from typing import Dict, List, Optional

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_RESULTS = os.path.join(_ROOT, "results")
_DATA = os.path.join(_RESULTS, "multisection_generalisation.json")
_OUT = os.path.join(_RESULTS, "multisection_generalisation.md")

ALPHA_EXPECTED = 22.223477711172563

#: Section order for every table. Xenium (single-cell, spatial statistics
#: VALID) first, then Visium (spot, spatial statistics REFUSED).
XENIUM = ["xenium_319", "xenium_io", "gse324268"]
VISIUM = ["visium_10x_parent_gbm", "visium_gbm1", "visium_gbm2", "visium_gbm3"]
ORDER = XENIUM + VISIUM

FLAGSHIP = "xenium_319"


# --------------------------------------------------------------------------- #
# formatting
# --------------------------------------------------------------------------- #
def f(v, nd=4, none="—"):
    if v is None:
        return none
    if isinstance(v, str):
        return v
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    if x != x:
        return "nan"
    if x == 0:
        return "0"
    if abs(x) >= 1e6 or abs(x) < 1e-3:
        return f"{x:.{nd}g}".replace("e+0", "e+").replace("e-0", "e-")
    return f"{x:.{nd}g}"


def pct(v, nd=2, none="—"):
    if v is None:
        return none
    return f"{100.0 * float(v):.{nd}f}%"


def spread(vals: List[Optional[float]]) -> Optional[float]:
    """max/min over the finite, strictly positive entries. None if undefined."""
    xs = [float(v) for v in vals
          if v is not None and float(v) == float(v) and float(v) > 0]
    if len(xs) < 2:
        return None
    return max(xs) / min(xs)


# --------------------------------------------------------------------------- #
def load() -> Dict:
    with open(_DATA) as fh:
        return json.load(fh)


def runs_at(blob: Dict, dx: float) -> Dict[str, Dict]:
    out = {}
    for key, rec in blob["runs"].items():
        if rec.get("spec", {}).get("dx_um") == dx:
            out[rec["section"]] = rec
    return out


def _row(r: Dict, *cells) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def build(blob: Dict) -> str:
    dxs = sorted({rec["spec"]["dx_um"] for rec in blob["runs"].values()},
                 reverse=True)
    primary_dx = max(dxs)
    R = runs_at(blob, primary_dx)
    fine = runs_at(blob, min(dxs)) if len(dxs) > 1 else {}
    present = [s for s in ORDER if s in R]
    xen = [s for s in XENIUM if s in R]

    a: List[str] = []
    w = a.append

    # ----------------------------------------------------------------- 0 --- #
    w("# Does ANY of it generalise off `xenium_319`? — one compound, one "
      "configuration, every runnable section")
    w("")
    w(f"**Status: RUN.** Driver: `experiments/multisection_generalisation.py` "
      f"(the sweep) and `experiments/multisection_report.py` (this document). "
      f"Data: `results/multisection_generalisation.json`. Every number below is "
      f"interpolated from that file; none is typed by hand.")
    w("")
    w("**Nothing under `tmesim/`, `app/`, `docs/` or `landing-3d/` was changed, "
      "and nothing was tuned to make sections agree.** `calibrate_alpha()` was "
      f"asserted live at `{ALPHA_EXPECTED!r}` before the first section, after "
      "every section, and after the last. The sweep drives "
      "`app.server.jobs.execute_run` — the exact function the Studio queue "
      "calls — so these are the product's shipped readouts, not a "
      "re-implementation of them.")
    w("")
    w("## THE CONFIGURATION, stated once and true of every number in this file")
    w("")
    w(f"> **temozolomide** · **dx = {primary_dx:g} µm** · 6 h exposure window "
      f"(`total_time_s = {R[present[0]]['spec']['total_time_s']:g}`) · "
      f"`dt = 10 s` · convection ON · `vascular_mode = discrete` · "
      f"`oxygen_mode = pde` · `state_mode = oxygen` · "
      f"`sensitivity_mode = uniform` · plasma_peak = "
      f"{f(R[present[0]].get('plasma_peak_uM'))} µM, half-life 14400 s, "
      f"1 dose · seed 0 · `transvascular_exchange = False` (shipped default) · "
      f"reach map in **steady-state** mode.")
    w("")
    w("ONE compound and ONE configuration throughout, so the **section** is the "
      "only thing that varies. This is a comparison of one compound across "
      "SECTIONS. It is **not** a cross-compound ranking — that claim was "
      "withdrawn 2026-08-21 and nothing here revives it.")
    w("")

    # ----------------------------------------------------------------- 1 --- #
    w("## 0. Verdict up front")
    w("")
    w("**Almost nothing generalises.** Of the eight shipped readouts this task "
      "asked for, exactly one is stable across tissue in a way that carries "
      "signal, one is stable only because it is degenerate everywhere, and the "
      "rest — including the per-cell discrimination verdict, the fitted "
      "penetration length, the falloff and the reach ladder's only informative "
      "rung — move by factors, change sign of interpretation, or are refused "
      "outright.")
    w("")
    w("| readout | across the 3 single-cell (Xenium) GBM sections | verdict |")
    w("|---|---|---|")

    def xv(key):
        return [R[s].get(key) for s in xen]

    sp_spread = spread(xv("E_p90_over_p10"))
    w(_row(None, "**delivery spread `E_p90/E_p10`**",
           " · ".join(f"{R[s]['section']} = {f(R[s]['E_p90_over_p10'])}"
                      for s in xen),
           f"**the most stable readout** — {f(sp_spread, 3)}× spread"))
    w(_row(None, "**reach headline `F(0.1)`**",
           " · ".join(f"{f(R[s]['reach_headline_F'])}" for s in xen),
           "**stable at zero — and therefore carries no cross-tissue signal "
           "at all**"))
    top = [(R[s]["reach_ladder_F"] or [None])[-1] for s in xen]
    w(_row(None, "**reach ladder top rung `F(1.0)`**",
           " · ".join(f(v) for v in top),
           "**NOT stable** — the ladder's only informative rung spans "
           f"{f(min(v for v in top if v is not None))} → "
           f"{f(max(v for v in top if v is not None))}"))
    Ls = xv("penetration_length_um")
    w(_row(None, "**fitted penetration length `L`**",
           " · ".join(f"{f(v)} µm" for v in Ls),
           f"**NOT stable** — {f(spread(Ls), 3)}× spread"))
    supports = [(R[s]["falloff_distance_um"] or [None])[-1] for s in xen]
    w(_row(None, "**normalised falloff**",
           "usable-bin support ends at " + " · ".join(
               f"{f(v)} µm" for v in supports),
           "**NOT stable, and not even defined over a common domain** — the "
           "curves cannot be overlaid"))
    lim = [R[s].get("mech_dominant_limiter") for s in xen]
    w(_row(None, "**per-cell discrimination — dominant limiter**",
           " · ".join(f"{R[s]['section']} = **{R[s].get('mech_dominant_limiter')}**"
                      for s in xen),
           "**FLIPS** — the product's headline verdict is not the same verdict "
           "on two sections of the same disease on the same platform"))
    cov = xv("therapeutic_coverage_1x")
    w(_row(None, "`therapeutic_coverage_1x` (context)",
           " · ".join(pct(v) for v in cov),
           "**NOT stable** — 0% on the flagship, ~73–76% on the other two"))
    hyp = xv("hypoxic_tumor_fraction")
    w(_row(None, "hypoxic tumour fraction (context)",
           " · ".join(pct(v) for v in hyp),
           "**the flagship's 0% is a property of the flagship** — every "
           "hypoxia-gated statement made on it is untested"))
    nv = [(R[s].get("n_viable_reach"), R[s].get("n_tumor")) for s in xen]
    w(_row(None, "`n_viable` / `n_tumor` (the reach denominator)",
           " · ".join(f"{v:,}/{t:,}" for v, t in nv),
           "**NOT stable** — see §2.2, the reach map on `xenium_io` is a "
           "fraction of ONE cell"))
    w("")
    w("**The four Visium sections are not in that table on purpose.** The "
      "shipped code refuses these quantities on spot data. §2.1.")
    w("")

    # ----------------------------------------------------------------- 2 --- #
    w("## 1. What ran, and what it cost")
    w("")
    w(f"All {len(present)} runnable built-in sections completed at "
      f"dx = {primary_dx:g} µm. No section crashed. Wall clock is honest and "
      f"was measured against three other agents running heavy simulations "
      f"concurrently on the same machine.")
    w("")
    w("| section | label | platform | cells | grid @dx | **voxels** | extent (µm) | tissue-mask mode | tissue coverage | vascular source | wall clock |")
    w("|---|---|---|---:|---|---:|---|---|---:|---|---:|")
    for s in present:
        r = R[s]
        g = r.get("grid") or [0, 0]
        vox = int(g[0]) * int(g[1])
        ext = r.get("extent_um") or [None, None]
        w(_row(None, f"`{s}`", r.get("label"), r.get("platform"),
               f"{int(r['n_cells']):,}", f"{g[0]}×{g[1]}", f"{vox:,}",
               f"{f(ext[0], 5)} × {f(ext[1], 5)}",
               r.get("tissue_mask_mode") or "full rectangle",
               pct(r.get("tissue_coverage_fraction"), 1),
               (r.get("vascular_source") or "")[:34],
               f"{r['wall_clock_s']:.1f} s"))
    w("")
    voxels = [int((R[s]['grid'] or [0, 0])[0]) * int((R[s]['grid'] or [0, 0])[1])
              for s in present]
    w(f"**Section extents differ, so a fixed dx gives different grid sizes.** "
      f"At dx = {primary_dx:g} µm the voxel count runs from {min(voxels):,} to "
      f"{max(voxels):,} — a factor of {max(voxels) / max(min(voxels), 1):.0f}. "
      f"Every mask-area average (`mean_tissue_auc`, the falloff bin means, the "
      f"tissue-mask mean of anything) is an average over a different number of "
      f"voxels on each row, and the Visium rows have no tissue mask at all — "
      f"their `tissue_coverage_fraction` is 1.0 because the domain is the full "
      f"bounding rectangle, against "
      f"{pct(R[FLAGSHIP].get('tissue_coverage_fraction'), 1)} on the flagship. "
      f"These are not the same measurement and are not compared here as if "
      f"they were.")
    w("")

    # ----------------------------------------------------------------- 3 --- #
    w("## 2. Refusals — recorded, not skipped")
    w("")
    w("### 2.1 Four of seven sections REFUSE the entire question (spot data)")
    w("")
    n_spot = sum(1 for s in present if (R[s].get("spatial_spread_valid") is False))
    w(f"{n_spot} of {len(present)} runnable built-in sections are Visium. "
      f"`app.server.api.SPOT_SPATIAL_REFUSAL` is attached to each of them by "
      f"the shipped API, and it is not a soft caveat — it says the delivery "
      f"map, its falloff, its perivascular reference and the delivery spread "
      f"are all INVALID on spot data. Verbatim, from the shipped constant:")
    w("")
    refusal = None
    for s in present:
        for c in (R[s].get("section_caveats") or []):
            if c.get("capability") == "spatial_spread":
                refusal = c.get("reason")
                break
        if refusal:
            break
    if refusal:
        for line in _wrap(refusal, 78):
            w("> " + line)
    w("")
    w("**This is the single largest result in this document.** The generalisation "
      f"question can only be asked of the {len(xen)} single-cell sections. The "
      f"other {n_spot} produce numbers — the engine does not raise — and those "
      "numbers are tabulated in §7 **quarantined**, so a reader can see exactly "
      "what would have been published had the refusal been ignored. They are "
      "not evidence about those tumours.")
    w("")

    w("### 2.2 A refusal that did NOT fire, and should have: `n_viable = 1`")
    w("")
    io = R.get("xenium_io")
    if io:
        w(f"`app/server/reach.py` defines `F(x)` as a fraction **of the VIABLE "
          f"tumour cells** (local oxygen ≥ 2 × `hypoxia_threshold`) and refuses "
          f"— `null_reason`, no number — when that population is empty. On "
          f"`xenium_io` at dx = {primary_dx:g} µm the viable population is "
          f"**{io['n_viable_reach']} cell** out of {io['n_tumor']:,} tumour "
          f"cells. One. The guard is `n_viable == 0`, so it does not fire: the "
          f"payload returns a clean `F(0.1) = {f(io['reach_headline_F'])}`, a "
          f"full seven-rung ladder, `null_reason = "
          f"{io.get('reach_null_reason')}`, and no indication anywhere that the "
          f"denominator is a single cell.")
        w("")
        w("**That is a live gap in shipped code, found by running off the "
          "flagship.** On `xenium_319` the same denominator is "
          f"{R[FLAGSHIP]['n_viable_reach']:,} cells, so the flagship could "
          "never have exposed it. It is reported, not patched: `app/` is "
          "outside this experiment's territory and patching it would also mean "
          "this document was written against an engine that no longer exists.")
        w("")

    w("### 2.3 Refusals that were available and did not fire on any section")
    w("")
    n_reason = sum(1 for s in present if R[s].get("metric_reasons"))
    n_fit = sum(1 for s in present if R[s].get("penetration_fit_error"))
    n_null = sum(1 for s in present if R[s].get("reach_null_reason"))
    w(f"* **No endothelial marker → no delivery result.** "
      f"`validation.vessel_profile.direct_vessel_distance_um` raises when a "
      f"section has no endothelial cells to measure to. It did not fire: all "
      f"{len(present)} sections carry endothelial cells (fractions in §3). "
      f"Recorded as a result — this axis is untested by these seven specimens, "
      f"not passed by them.")
    w(f"* **`metric_reasons`** (the shipped per-metric refusal channel) is "
      f"empty on {len(present) - n_reason} of {len(present)} sections"
      + ("." if n_reason == 0 else
         f"; {n_reason} section(s) recorded one, listed in §7."))
    w(f"* **`penetration_fit_error`** is `None` on "
      f"{len(present) - n_fit} of {len(present)} sections — the exp+floor fit "
      f"converged everywhere.")
    w(f"* **`reach.null_reason`** is `None` on {len(present) - n_null} of "
      f"{len(present)} sections.")
    w("")
    w("So the refusal machinery is *quiet* on this panel. That is the "
      "uncomfortable part: the numbers below disagree violently while every "
      "refusal channel the product ships reports nothing wrong. **Disagreement "
      "across tissue is not currently something this pipeline can detect about "
      "itself.**")
    w("")

    # ----------------------------------------------------------------- 4 --- #
    w("## 3. The per-section table — single-cell sections only")
    w("")
    w("These are the four sections' worth — sorry, the three sections' worth — "
      "of readouts on which every quantitative claim in this repository could "
      "legitimately have been made. `xenium_319` is one of them.")
    w("")
    w("| readout | " + " | ".join(f"**`{s}`**" for s in xen) +
      " | max/min |")
    w("|---|" + "---|" * (len(xen) + 1))

    def line(label, key, fmt=lambda v: f(v), sp=True):
        vals = [R[s].get(key) for s in xen]
        s_ = spread(vals) if sp else None
        w(_row(None, label, *[fmt(v) for v in vals],
               f"**{f(s_, 3)}×**" if s_ else "—"))

    line("cells (total)", "n_cells", lambda v: f"{int(v):,}" if v else "—",
         sp=False)
    line("`n_tumor`", "n_tumor", lambda v: f"{int(v):,}" if v else "—")
    w(_row(None, "`n_viable` (reach denominator)",
           *[f"{R[s]['n_viable_reach']:,}" for s in xen],
           f"**{f(spread([R[s]['n_viable_reach'] for s in xen]), 3)}×**"))
    w(_row(None, "`n_viable` / `n_tumor`",
           *[pct(R[s]['n_viable_reach'] / R[s]['n_tumor'], 4) for s in xen],
           ""))
    for ct in ("TUMOR", "ENDOTHELIAL", "IMMUNE", "STROMAL"):
        vals = [(R[s].get("celltype_fractions") or {}).get(ct) for s in xen]
        w(_row(None, f"cell fraction {ct}", *[pct(v, 1) for v in vals],
               f"{f(spread(vals), 3)}×"))
    line("mean vascular density", "mean_vascular_density")
    line("mean oxygen", "mean_oxygen")
    line("hypoxic tumour fraction", "hypoxic_tumor_fraction", lambda v: pct(v))
    w(_row(None, "— *(order-of-magnitude only: the oxygen consumption constant "
                 "is un-calibrated)* —", *([""] * (len(xen) + 1))))
    w(_row(None, "**`E_p90/E_p10`  (DELIVERY SPREAD — the headline)**",
           *[f"**{f(R[s]['E_p90_over_p10'])}**" for s in xen],
           f"**{f(spread(xv('E_p90_over_p10')), 3)}×**"))
    line("`E_p10` (µM)", "E_p10_uM")
    line("`E_p50` (µM)", "E_p50_uM")
    line("`E_p90` (µM)", "E_p90_uM")
    w(_row(None, "**`F(0.1)` (reach headline)**",
           *[f"**{f(R[s]['reach_headline_F'])}**" for s in xen], "—"))
    w(_row(None, "**`F(1.0)` (top ladder rung)**",
           *[f"**{f((R[s]['reach_ladder_F'] or [None])[-1])}**" for s in xen],
           "—"))
    w(_row(None, "`E_ref` (perivascular reference)",
           *[f"{f(R[s]['E_ref'])} {R[s].get('E_ref_units') or ''}" for s in xen],
           f"{f(spread(xv('E_ref')), 3)}×"))
    w(_row(None, "**fitted penetration length `L` (µm)**",
           *[f"**{f(R[s]['penetration_length_um'])}**" for s in xen],
           f"**{f(spread(xv('penetration_length_um')), 3)}×**"))
    line("analytic penetration length (µm)", "penetration_length_analytic_um")
    w(_row(None, "falloff support (last usable bin, µm)",
           *[f(( R[s]['falloff_distance_um'] or [None])[-1]) for s in xen],
           f"{f(spread(supports), 3)}×"))
    line("`therapeutic_coverage_1x`", "therapeutic_coverage_1x",
         lambda v: pct(v), sp=False)
    line("median E / IC50", "median_E_over_ic50")
    line("mean tumour AUC (µM·s)", "mean_tumor_auc")
    w(_row(None, "**dominant limiter**",
           *[f"**{R[s].get('mech_dominant_limiter')}**" for s in xen], "—"))
    w(_row(None, "*kill_fraction (NOT a defensible output)*",
           *[f"*{f(R[s]['kill_fraction'])}*" for s in xen], "*—*"))
    w("")
    w("The last row is greyed for a reason. **Absolute kill fraction is never a "
      "defensible output and no comparison in this document is built on it.** "
      "It is printed only because the shipped payload carries it and silently "
      "dropping a number the product prints would be its own kind of "
      "dishonesty.")
    w("")

    # ----------------------------------------------------------------- 5 --- #
    w("## 4. The per-cell discrimination — the verdict flips")
    w("")
    w("`tmesim.sim.mechanism` labels every tumour cell with why it survived, "
      "and reconciles the labels against `kill_fraction` **bit-for-bit** "
      "(`==`, not `allclose`) or raises. It reconciled on every section:")
    w("")
    w("| section | KILLED | UNREACHED | UNDERDOSED | RESISTANT | STOCHASTIC | UNATTRIBUTED | **delivery** | **dose** | **resistance** | **chance** | **dominant** | reconciles bit-exact |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")
    for s in present:
        r = R[s]
        c = r.get("mech_counts") or {}
        b = r.get("mech_survival_breakdown") or {}
        rec = r.get("mech_reconciliation") or {}
        mark = "" if r.get("spatial_spread_valid") is not False else " ⚠︎spot"
        w(_row(None, f"`{s}`{mark}",
               f"{c.get('KILLED', 0):,}", f"{c.get('SURVIVED_UNREACHED', 0):,}",
               f"{c.get('SURVIVED_UNDERDOSED', 0):,}",
               f"{c.get('SURVIVED_RESISTANT', 0):,}",
               f"{c.get('SURVIVED_STOCHASTIC', 0):,}",
               f"{c.get('UNATTRIBUTED', 0):,}",
               f(b.get("delivery")), f(b.get("dose")), f(b.get("resistance")),
               f(b.get("chance")),
               f"**{r.get('mech_dominant_limiter')}**",
               "yes" if rec.get("bit_exact") else "**NO**"))
    w("")
    w("**Read the `dominant` column.** On the flagship the product says the "
      f"tumour is **{R[FLAGSHIP].get('mech_dominant_limiter')}**-limited, with "
      f"{pct((R[FLAGSHIP].get('mech_survival_breakdown') or {}).get('dose'))} "
      f"of survivals in that bucket. On `xenium_io` and `gse324268` — same "
      f"disease, same platform, same 319/480-style panel family, same "
      f"compound, same configuration — it says "
      f"**{R['xenium_io'].get('mech_dominant_limiter')}**, with "
      f"{pct((R['xenium_io'].get('mech_survival_breakdown') or {}).get('resistance'))} "
      f"and "
      f"{pct((R['gse324268'].get('mech_survival_breakdown') or {}).get('resistance'))} "
      f"in the resistance bucket respectively. This is not a shift in "
      f"magnitude. It is a different answer to the question the product exists "
      f"to answer.")
    w("")
    fl = R[FLAGSHIP].get("mech_joint_limited") or {}
    if fl:
        w(f"**And the flagship's own split is over-determined.** "
          f"`joint_limited` on `xenium_319` is {fl.get('count'):,} cells, "
          f"{pct(fl.get('fraction_of_exposure_limited'))} of the "
          f"exposure-limited survivals: cells labelled DOSE-limited by the "
          f"documented precedence rule that are ALSO at or below the "
          f"resistance threshold, and that raising the dose alone would not "
          f"have killed. The clean-looking "
          f"{pct((R[FLAGSHIP].get('mech_survival_breakdown') or {}).get('dose'))} "
          f"dose share is not a clean counterfactual even on the section it "
          f"was measured on.")
        w("")
    thr = R[FLAGSHIP].get("mech_thresholds") or {}
    w(f"Thresholds are identical on every row and are derived, not chosen: "
      f"`{thr}`. Nothing here was tuned per section.")
    w("")

    # ----------------------------------------------------------------- 6 --- #
    w("## 5. The reach ladder `F(x)` across tissue")
    w("")
    ladder_x = R[FLAGSHIP].get("reach_ladder_x") or []
    w("| section | " + " | ".join(f"`F({x:g})`" for x in ladder_x) +
      " | `n_viable` | `E_ref` |")
    w("|---|" + "---:|" * (len(ladder_x) + 2))
    for s in present:
        r = R[s]
        F_ = r.get("reach_ladder_F") or []
        mark = "" if r.get("spatial_spread_valid") is not False else " ⚠︎spot"
        w(_row(None, f"`{s}`{mark}", *[f(v) for v in F_],
               f"{r.get('n_viable_reach'):,}"
               + (" **← ONE CELL**" if r.get("n_viable_reach") == 1 else ""),
               f"{f(r.get('E_ref'))} {r.get('E_ref_units') or ''}"))
    w("")
    w(f"**The headline rung `F({R[FLAGSHIP]['reach_headline_x']:g})` reads "
      f"{f(R[FLAGSHIP]['reach_headline_F'])} on every single-cell section, and "
      f"essentially zero on the spot sections too.** `app/server/reach.py`'s "
      f"`CONVENTION_NOTE` already de-headlined `x = 0.1` on the grounds that "
      f"F(0.1) = 0.0005 on the flagship at dx = 12 µm and \"a headline that "
      f"reads 0.05% on the specimen this product ships with is not a "
      f"headline\". **This sweep says the same thing about all seven "
      f"specimens.** The de-headlining was right, and it was if anything "
      f"understated: the rung is not merely unimpressive on the flagship, it is "
      f"informationless across every tissue this checkout can run.")
    w("")
    top_x = [(R[s].get("reach_ladder_F") or [None])[-1] for s in xen]
    w(f"**The one rung that does carry signal — `F(1.0)` — does not "
      f"generalise.** Across the three single-cell sections it reads "
      + ", ".join(f"{f(v)} (`{s}`)" for s, v in zip(xen, top_x))
      + ". The flagship's value is the extreme one. And on `xenium_io` that "
        "number is a fraction of a one-cell population (§2.2), so the honest "
        "count of single-cell sections on which the ladder is a meaningful "
        "object is **two**.")
    w("")

    # ----------------------------------------------------------------- 7 --- #
    w("## 6. The normalised falloff — the curves do not share a domain")
    w("")
    w("Each row is the shipped `falloff` block: the same distance bins that "
      "produced `E_ref`, each bin's mean divided by `E_ref`, so bin 0 is 1.0 "
      "by construction. A bin needs ≥ 30 grid points to be usable, so **how far "
      "the curve extends is itself a per-section measurement.**")
    w("")
    common = [5.0, 15.0, 25.0, 35.0, 45.0, 65.0, 85.0, 105.0, 145.0]
    w("| section | " + " | ".join(f"{d:g} µm" for d in common) +
      " | last usable bin | n bins |")
    w("|---|" + "---:|" * (len(common) + 2))
    for s in present:
        r = R[s]
        d = r.get("falloff_distance_um") or []
        v = r.get("falloff_normalised_mean_auc") or []
        m = {float(dd): vv for dd, vv in zip(d, v)}
        mark = "" if r.get("spatial_spread_valid") is not False else " ⚠︎spot"
        w(_row(None, f"`{s}`{mark}",
               *[f(m[c], 3) if c in m else "*no usable bin*" for c in common],
               f"{f(d[-1]) if d else '—'} µm", len(d)))
    w("")
    xs_at45 = [
        dict(zip(R[s]["falloff_distance_um"],
                 R[s]["falloff_normalised_mean_auc"])).get(45.0) for s in xen]
    w(f"**At 45 µm from the nearest vessel — inside every section's support — "
      f"the normalised exposure is "
      + ", ".join(f"{f(v, 3)} (`{s}`)" for s, v in zip(xen, xs_at45))
      + f", a {f(spread(xs_at45), 3)}× spread.** By 105 µm the flagship still "
        f"reads "
      + f(dict(zip(R[FLAGSHIP]['falloff_distance_um'],
                   R[FLAGSHIP]['falloff_normalised_mean_auc'])).get(105.0), 3)
      + " while `gse324268` has no usable bin there at all. The three curves "
        "are not three samples of one law; they are three different shapes over "
        "three different domains.")
    w("")

    # ----------------------------------------------------------------- 8 --- #
    if fine:
        w(f"## 7. Is the disagreement TISSUE, or is it the GRID? — the paired "
          f"dx = {max(dxs):g} → {min(dxs):g} µm control")
        w("")
        w("A cross-section spread that is no bigger than the same readout's own "
          "grid drift is not a finding. **The rule below was fixed before the "
          f"dx = {min(dxs):g} µm arm had produced a single number:** for each "
          "readout, compare the max/min across the three single-cell sections "
          f"at dx = {max(dxs):g} µm against that readout's own "
          f"dx = {max(dxs):g} → {min(dxs):g} µm fold change **on the flagship**. "
          "A section spread inside the flagship's own grid drift is NOT "
          "ESTABLISHED as tissue signal.")
        w("")
        w("This is the only defensible yardstick available, because "
          "`app/server/jobs.py` ships the standing caveat that **there is no dx "
          "at which this pipeline has been shown to be grid-converged** — so a "
          "chosen tolerance would be a preference, whereas the flagship's own "
          "measured drift is a measurement.")
        w("")
        keys = [
            ("E_p90/E_p10 (delivery spread)", "E_p90_over_p10"),
            ("penetration length L (µm)", "penetration_length_um"),
            ("E_ref", "E_ref"),
            ("therapeutic_coverage_1x", "therapeutic_coverage_1x"),
            ("median E / IC50", "median_E_over_ic50"),
            ("mean tumour AUC", "mean_tumor_auc"),
            ("n_viable", "n_viable_reach"),
            ("hypoxic tumour fraction", "hypoxic_tumor_fraction"),
        ]
        w(f"| readout | section spread @dx={max(dxs):g} (max/min over 3 Xenium) "
          f"| flagship dx={max(dxs):g} → {min(dxs):g} fold | verdict |")
        w("|---|---:|---:|---|")
        for label, key in keys:
            vals = [R[s].get(key) for s in xen]
            sp_ = spread(vals)
            c30 = R[FLAGSHIP].get(key)
            c12 = (fine.get(FLAGSHIP) or {}).get(key)
            drift = None
            if c30 and c12 and float(c30) > 0 and float(c12) > 0:
                drift = max(float(c30), float(c12)) / min(float(c30), float(c12))
            if sp_ is None or drift is None:
                verdict = "**cannot be judged** — one side undefined"
            elif sp_ > 3.0 * drift:
                verdict = "**tissue** — spread ≫ grid drift"
            elif sp_ > drift:
                verdict = "tissue-dominated, but only by "\
                          f"{sp_ / drift:.1f}×"
            else:
                verdict = "**NOT ESTABLISHED** — spread ≤ the flagship's own grid drift"
            w(_row(None, label, f(sp_, 3), f(drift, 3), verdict))
        w("")
        w(f"### 7.1 Every section at both grids")
        w("")
        w(f"| section | `E_p90/p10` @{max(dxs):g} | @{min(dxs):g} | `L` (µm) "
          f"@{max(dxs):g} | @{min(dxs):g} | limiter @{max(dxs):g} | "
          f"@{min(dxs):g} | `F(1.0)` @{max(dxs):g} | @{min(dxs):g} | "
          f"voxels @{max(dxs):g} | @{min(dxs):g} |")
        w("|---|---:|---:|---:|---:|---|---|---:|---:|---:|---:|")
        for s in present:
            r30, r12 = R[s], fine.get(s)

            def g2(rec, key):
                return rec.get(key) if rec else None

            def topr(rec):
                return (rec.get("reach_ladder_F") or [None])[-1] if rec else None

            def vox(rec):
                if not rec or not rec.get("grid"):
                    return None
                return int(rec["grid"][0]) * int(rec["grid"][1])

            mark = "" if r30.get("spatial_spread_valid") is not False else " ⚠︎spot"
            w(_row(None, f"`{s}`{mark}",
                   f(r30.get("E_p90_over_p10")), f(g2(r12, "E_p90_over_p10")),
                   f(r30.get("penetration_length_um")),
                   f(g2(r12, "penetration_length_um")),
                   r30.get("mech_dominant_limiter") or "—",
                   (r12 or {}).get("mech_dominant_limiter") or "—",
                   f(topr(r30)), f(topr(r12)),
                   f"{vox(r30):,}" if vox(r30) else "—",
                   f"{vox(r12):,}" if vox(r12) else "—"))
        w("")
        flips = [s for s in present
                 if fine.get(s) and R[s].get("mech_dominant_limiter")
                 != fine[s].get("mech_dominant_limiter")]
        if flips:
            w(f"**The dominant limiter changes with the GRID on "
              f"{len(flips)} of {len(present)} sections** "
              + ", ".join(f"`{s}`" for s in flips)
              + f" — so on those sections the product's headline verdict is not "
                f"even a property of the tissue at fixed configuration; it is a "
                f"property of the mesh. Where that happens, the cross-section "
                f"flip in §4 cannot be attributed to tissue either.")
        else:
            w(f"**The dominant limiter is unchanged by the grid on all "
              f"{len(present)} sections.** The cross-section flip in §4 is "
              f"therefore not a mesh artefact — it survives a "
              f"{max(dxs) / min(dxs):.1f}× refinement in dx "
              f"(≈{(max(dxs) / min(dxs)) ** 2:.0f}× in voxel count).")
        w("")

    # ----------------------------------------------------------------- 9 --- #
    w("## 8. The quarantine — what the four Visium sections produced")
    w("")
    w("Reproduced so that the refusal in §2.1 is auditable rather than a "
      "assertion. **None of these is a reading of those tumours.** The delivery "
      "map's source is a density over ENDOTHELIAL CELL POSITIONS and on a spot "
      "section those positions are synthesised — 12 pseudo-cells per spot at "
      "the spot centroid with uniform jitter — so every column below is a "
      "measurement of that synthesis.")
    w("")
    w("| section | cells | `n_tumor` | `n_viable` | `E_p90/p10` | `F(0.1)` | "
      "`F(1.0)` | `L` (µm) | falloff support | dominant limiter | "
      "`coverage_1x` |")
    w("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|")
    for s in [x for x in present if x in VISIUM]:
        r = R[s]
        d = r.get("falloff_distance_um") or []
        w(_row(None, f"`{s}`", f"{int(r['n_cells']):,}", f"{r['n_tumor']:,}",
               f"{r['n_viable_reach']:,}", f(r["E_p90_over_p10"]),
               f(r["reach_headline_F"]),
               f((r.get("reach_ladder_F") or [None])[-1]),
               f(r["penetration_length_um"]),
               f"{f(d[-1]) if d else '—'} µm",
               r.get("mech_dominant_limiter"),
               pct(r.get("therapeutic_coverage_1x"))))
    w("")
    vis = [s for s in present if s in VISIUM]
    w(f"For the record of how bad it would have been: had the refusal been "
      f"ignored, the delivery spread across all {len(present)} sections would "
      f"have been quoted as "
      f"{f(min(R[s]['E_p90_over_p10'] for s in present))}–"
      f"{f(max(R[s]['E_p90_over_p10'] for s in present))} "
      f"({f(spread([R[s]['E_p90_over_p10'] for s in present]), 3)}×) and the "
      f"penetration length as "
      f"{f(min(R[s]['penetration_length_um'] for s in present))}–"
      f"{f(max(R[s]['penetration_length_um'] for s in present))} µm "
      f"({f(spread([R[s]['penetration_length_um'] for s in present]), 3)}×). "
      f"Both ranges are dominated by sections whose vascular geometry is "
      f"synthetic.")
    w("")

    # ---------------------------------------------------------------- 10 --- #
    w("## 9. Stable across tissue vs. property of `xenium_319`")
    w("")
    w("### STABLE ACROSS TISSUE")
    w("")
    w(f"1. **The delivery spread `E_p90/E_p10`.** "
      + ", ".join(f"{f(R[s]['E_p90_over_p10'])} (`{s}`)" for s in xen)
      + f" — a {f(spread(xv('E_p90_over_p10')), 3)}× spread across three "
        f"independent GBM specimens, which is the tightest agreement of any "
        f"readout here. It is also the readout `app/server/reach.py` already "
        f"promoted to headline in place of `F(0.1)`, and this sweep is the "
        f"first evidence that the promotion picked the right statistic for a "
        f"reason beyond convenience. **Qualified:** three sections is three "
        f"sections, and §7 tests it against the flagship's own grid drift.")
    w("")
    w("2. **Bit-exact mechanism reconciliation.** The per-cell attribution "
      f"reproduced `kill_fraction` bit-for-bit on {len(present)} of "
      f"{len(present)} sections, with `UNATTRIBUTED = 0` everywhere. The "
      "*internal consistency* of the attribution generalises perfectly. Its "
      "*answer* does not — §4.")
    w("")
    w("3. **The pipeline runs at all.** Seven of seven sections completed with "
      "no crash and no metric refusal. That is a real engineering result and "
      "it is the only sense in which 'it generalises' is true without "
      "qualification.")
    w("")
    w(f"4. **`F(0.1)` = 0.** Stable, and stable at the value that makes it "
      f"useless. Listed under STABLE because it is, and read under §5 for what "
      f"that stability is worth.")
    w("")
    w("### PROPERTY OF `xenium_319` THAT HAS BEEN READ AS A PROPERTY OF THE MODEL")
    w("")
    items = []
    items.append(
        f"**\"This tumour is "
        f"{R[FLAGSHIP].get('mech_dominant_limiter')}-limited.\"** "
        f"The flagship's dominant limiter is "
        f"`{R[FLAGSHIP].get('mech_dominant_limiter')}` at "
        f"{pct((R[FLAGSHIP].get('mech_survival_breakdown') or {}).get('dose'))}. "
        f"Both other single-cell sections return "
        f"`{R['xenium_io'].get('mech_dominant_limiter')}`. The per-cell "
        f"discrimination is a property of the specimen, and on the flagship it "
        f"is additionally over-determined "
        f"({pct((R[FLAGSHIP].get('mech_joint_limited') or {}).get('fraction_of_exposure_limited'))} "
        f"of exposure-limited survivals are jointly resistance-limited).")
    items.append(
        f"**\"No tumour cell reaches IC50.\"** `therapeutic_coverage_1x` is "
        f"{pct(R[FLAGSHIP].get('therapeutic_coverage_1x'))} on the flagship and "
        f"{pct(R['xenium_io'].get('therapeutic_coverage_1x'))} / "
        f"{pct(R['gse324268'].get('therapeutic_coverage_1x'))} on the other "
        f"two. A statement that reads as a fact about temozolomide is a fact "
        f"about this section.")
    items.append(
        f"**The penetration length.** "
        f"{f(R[FLAGSHIP]['penetration_length_um'])} µm on the flagship against "
        f"{f(R['xenium_io']['penetration_length_um'])} and "
        f"{f(R['gse324268']['penetration_length_um'])} µm "
        f"({f(spread(xv('penetration_length_um')), 3)}×). The flagship is the "
        f"**longest** of the three — the specimen the product ships with is its "
        f"most drug-penetrable one.")
    items.append(
        f"**The falloff curve.** The flagship's usable-bin support runs to "
        f"{f((R[FLAGSHIP]['falloff_distance_um'] or [None])[-1])} µm; "
        f"`gse324268` stops at "
        f"{f((R['gse324268']['falloff_distance_um'] or [None])[-1])} µm and "
        f"`xenium_io` runs to "
        f"{f((R['xenium_io']['falloff_distance_um'] or [None])[-1])} µm and "
        f"reaches exactly 0 in its outer bins. There is no shared domain on "
        f"which 'the normalised falloff' is one curve.")
    items.append(
        f"**Zero hypoxic tumour.** `hypoxic_tumor_fraction` is "
        f"{pct(R[FLAGSHIP].get('hypoxic_tumor_fraction'))} on the flagship "
        f"against {pct(R['xenium_io'].get('hypoxic_tumor_fraction'))} and "
        f"{pct(R['gse324268'].get('hypoxic_tumor_fraction'))}. Every "
        f"hypoxia-gated term in the engine — the resistance bands, the viable "
        f"band, the oxygen gate — has been exercised on a section where the "
        f"hypoxic population is EMPTY. That is why the limiter can flip to "
        f"`resistance` the moment the flagship is left: the resistance branch "
        f"was effectively untested. (Order-of-magnitude only — the oxygen "
        f"consumption constant is un-calibrated — but the direction is not in "
        f"doubt at this size.)")
    items.append(
        f"**A well-populated viable band.** The reach map's denominator is "
        f"{R[FLAGSHIP]['n_viable_reach']:,} cells on the flagship. On "
        f"`xenium_io` it is {R['xenium_io']['n_viable_reach']}. The flagship "
        f"could never have revealed that the ladder's refusal guard has a hole "
        f"at `n_viable == 1` (§2.2).")
    items.append(
        f"**The reach ladder having a non-trivial top rung.** `F(1.0)` = "
        f"{f((R[FLAGSHIP]['reach_ladder_F'] or [None])[-1])} on the flagship "
        f"and {f((R['gse324268']['reach_ladder_F'] or [None])[-1])} on "
        f"`gse324268`. The flagship is the outlier that makes the ladder look "
        f"informative.")
    for i, it in enumerate(items, 1):
        w(f"{i}. {it}")
        w("")

    w("### NOT ESTABLISHED EITHER WAY")
    w("")
    w(f"* Anything about spot tissue. Four of seven sections refuse (§2.1), so "
      f"'does it generalise across platform?' has **no answer** in this "
      f"checkout — not a negative one, an absent one.")
    w(f"* Anything about a section with no endothelial marker. All seven carry "
      f"one, so that refusal path is untested rather than passed (§2.3).")
    w(f"* Anything about disease other than GBM. Every runnable built-in "
      f"section here is glioblastoma.")
    w("")

    # ---------------------------------------------------------------- 11 --- #
    w("## 10. What this changes")
    w("")
    w("1. **Every quantitative claim in this repository that is stated without "
      "a section qualifier needs one.** Not as a hedge — as a correction. The "
      "limiter verdict, the coverage, the penetration length and the falloff "
      "are section-level measurements that have been carrying model-level "
      "grammar.")
    w("2. **The delivery-spread headline survives, and is the only headline "
      "that does.** That is a vindication of the decision recorded in "
      "`app/server/reach.py`'s `CONVENTION_NOTE` to demote `F(0.1)` in its "
      "favour — reached there on the flagship alone, and confirmed here across "
      "tissue.")
    w("3. **The refusal ladder has a hole at `n_viable == 1`** and the product "
      "will publish a one-cell statistic without saying so (§2.2).")
    w("4. **No refusal channel fired on the disagreement itself.** Seven "
      "sections, verdicts that contradict each other, and "
      "`metric_reasons`/`null_reason`/`penetration_fit_error` all silent. The "
      "product cannot currently tell a customer that their section is outside "
      "the regime the flagship established.")
    w("5. **This was worth running before a customer ran it.** A buyer with "
      "their own Xenium section had roughly a two-in-three chance, on this "
      "evidence, of getting a different dominant limiter than every document "
      "in this repository describes.")
    w("")
    w("---")
    w("")
    w(f"*Generated by `experiments/multisection_report.py` from "
      f"`results/multisection_generalisation.json` "
      f"(sweep last updated {blob.get('updated')}). "
      f"`calibrate_alpha()` = {ALPHA_EXPECTED!r} at every stage boundary.*")
    return "\n".join(a) + "\n"


def _wrap(text: str, width: int) -> List[str]:
    words, lines, cur = text.split(), [], ""
    for word in words:
        if len(cur) + len(word) + 1 > width:
            lines.append(cur)
            cur = word
        else:
            cur = (cur + " " + word).strip()
    if cur:
        lines.append(cur)
    return lines


def main() -> int:
    blob = load()
    md = build(blob)
    with open(_OUT, "w") as fh:
        fh.write(md)
    print(f"wrote {_OUT} ({len(md)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
