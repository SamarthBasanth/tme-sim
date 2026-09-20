"""Does the oxygen sink's FUNCTIONAL FORM explain the missing hypoxia?

Controlled A/B of ``OxygenConfig.sink_form`` on the flagship GBM section, scored
against MEASURED intraoperative human-glioma pO2, with the decision rule frozen
in the PRE_REGISTRATION block below BEFORE any solve is run.

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
``results/oxygen_calibration_landscape.md`` Sec.4 argued ON PAPER (2026-07-28) that
the engine's first-order-everywhere sink ``-c*rho*O`` "will systematically
under-predict hypoxia", because real tissue consumption saturates (Michaelis-Menten)
rather than throttling off linearly as oxygen runs out. The engine reports
``HYPOXIC = 0.0000`` on the flagship section, which is the predicted signature.

**THE HYPOXIC-FRACTION HALF OF THIS A/B IS ALREADY SPENT.** See the
PRIOR_INFORMATION block: ``results/oxygen_sensitivity_band.md`` (2026-07-29) ran
both sink forms on this exact section at this exact ``consumption`` and published
the hypoxic fractions; ``results/oxygen_joint_grid_test.md`` (2026-08-15) published
the linear arm's median pO2 in mmHg. Those numbers were READ BEFORE THIS SCRIPT WAS
WRITTEN and are therefore NOT pre-registerable here. They are re-derived below as a
REPRODUCTION CHECK ONLY.

What is genuinely unmeasured -- and is what this script pre-registers -- is the
DISTRIBUTION: the Michaelis-Menten arm's median, quantile ladder and minimum in
mmHg, the sub-2.5 mmHg tail fraction under either form, and the scoring of both
against the measured human-glioma anchors.

Run: ``./.venv/bin/python experiments/oxygen_sink_form.py``
Writes: ``results/oxygen_sink_form_data.json``. Touches no module under ``tmesim/``.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from tmesim.core import CellType, OxygenConfig
from tmesim.physics.oxygen import REFERENCE_PO2_MMHG, solve_oxygen, to_mmhg
from tmesim.reconstruct.real_pipeline import reconstruct_real
from tmesim.validation.calibration import calibrate_alpha

REPO = Path(__file__).resolve().parent.parent
OUT_JSON = REPO / "results" / "oxygen_sink_form_data.json"

ALPHA_PINNED = 22.223477711172563

# ---------------------------------------------------------------------------
# PRIOR INFORMATION -- what was ALREADY PUBLISHED and was read BEFORE writing
# this script. Declared so that nothing below can be passed off as a blind
# prediction. Every one of these is a number I knew in advance.
# ---------------------------------------------------------------------------
PRIOR_INFORMATION = {
    "results/oxygen_sensitivity_band.md": {
        "date": "2026-07-29",
        "what_it_already_ran": (
            "BOTH sink forms, on reconstruct_real('xenium', dx_um=12.0), across an "
            "11-point consumption grid, at K=1/35. The row at the shipping default "
            "consumption=0.02 is the exact A/B this script re-runs."
        ),
        "known_in_advance": {
            "linear_hypoxic_fraction_at_c0.02": 0.0000,
            "linear_mean_O2_tumor_at_c0.02": 0.3486,
            "michaelis_menten_hypoxic_fraction_at_c0.02": 1.0000,
            "michaelis_menten_mean_O2_tumor_at_c0.02": 0.0347,
        },
    },
    "results/oxygen_joint_grid_test.md": {
        "date": "2026-08-15",
        "what_it_already_ran": (
            "The SHIPPED (linear) configuration's median pO2 in mmHg on this "
            "section at dx=12, plus p10/p90 for two hypothetical non-shipping "
            "configurations."
        ),
        "known_in_advance": {
            "shipped_linear_mean_O2_tumor": 0.348628,
            "shipped_linear_median_pO2_mmHg": 12.107,
            "shipped_linear_hypoxic_fraction": 0.0000,
        },
    },
    "results/sweep/oxsupply_sections.json": {
        "date": "2026-08-10",
        "what_it_already_ran": (
            "The LINEAR arm's ENTIRE per-tumour-cell pO2 percentile ladder in mmHg on "
            "this exact section at dx=12, plus the sub-2.5 and sub-10 mmHg fractions. "
            "Stored in raw JSON under sections/gbm/shipped_0.05 and NEVER surfaced in "
            "any results/*.md. Found by a prior-art sweep in this pass and read before "
            "measuring."
        ),
        "known_in_advance": {
            "tumour_o2_mmhg_percentiles": {
                "1": 10.21880049736665, "5": 10.59759610495022,
                "10": 10.859829506732337, "25": 11.523922427523368,
                "50": 12.106600909723596, "75": 12.630937942385584,
                "90": 14.004533189996282, "99": 14.707176344941557,
            },
            "tumour_frac_below_2p5mmhg": 0.0,
            "tumour_frac_below_10mmhg": 0.0006719301192675962,
            "mean_o2_tumour_mmhg": 12.201971805280552,
            "n_tumour_cells": 23812,
        },
    },
    "consequence": (
        "The hypoxic fraction under BOTH forms was known before measuring. Because "
        "the linear arm's hypoxic fraction (O < 0.2, i.e. < 7 mmHg) was already "
        "published as exactly 0.0, its sub-2.5 and sub-5 mmHg fractions are "
        "DERIVABLE as 0.0 and are likewise not pre-registerable. The linear arm's "
        "lower tail is entirely spent. What remains genuinely open is the "
        "Michaelis-Menten arm's DISTRIBUTION SHAPE (median, quantiles, minimum, "
        "and the split of its 100% sub-7 mmHg mass across the 2.5 and 5 mmHg "
        "thresholds), and the scoring of both arms against measured pO2."
    ),
    "the_linear_arm_is_fully_spent": (
        "With oxsupply_sections.json in hand, EVERY quantity this script computes for "
        "the LINEAR arm was already on disk before it ran. The linear arm is therefore "
        "reported here as a REPRODUCTION, not as a measurement. The genuinely new "
        "measurement is the MICHAELIS-MENTEN arm's mmHg distribution on a real "
        "reconstructed section, which the prior-art sweep confirmed does not exist "
        "anywhere in the repository -- the only both-forms mmHg execution "
        "(stale_blockers_2026-09-20.md Sec.2.3) was on a synthetic 64x64 lattice and "
        "says so."
    ),
}

# ---------------------------------------------------------------------------
# KNOWN CONFOUNDS -- declared BEFORE measuring, because each one bounds what any
# result below can mean. None of these is discovered after the fact.
# ---------------------------------------------------------------------------
CONFOUNDS = {
    "C1_the_source_constant_is_off_by_~4000x": (
        "results/oxygen_supply_calibration.md derives OxygenConfig.supply ~ 194 /s "
        "against the shipped 0.05 -- a factor of ~3,880 -- and states that 'Every "
        "published F_state, absolute kill and kill-derived ratio on the GBM section "
        "should be read as conditional on an oxygen source constant that is ~4,000x "
        "from its measured value.' BOTH arms of this A/B sit on that same source. "
        "The A/B is therefore a controlled comparison of two SINKS on a SOURCE that "
        "is independently known to be wrong, so no arm's ABSOLUTE agreement or "
        "disagreement with a measured pO2 can be attributed to the sink alone. What "
        "the A/B does isolate is the DIFFERENCE between the two forms, because the "
        "source is byte-identical across them."
    ),
    "C2_the_mmHg_anchor_is_itself_contested": (
        "REFERENCE_PO2_MMHG = 35.0 is a normal-brain-TISSUE anchor. "
        "results/oxygen_supply_calibration.md Sec.2 argues from Secomb's reference "
        "cortical values (P_cap = 58.75 mmHg, P_tiss = 47.65 mmHg) that 'a "
        "compartment-correct anchor for O = 1 would be ~55-60 mmHg, not 35, and "
        "every mmHg figure in this report would scale up by ~1.7x'. The anchor is "
        "FROZEN at the shipped 35.0 for the primary verdict. A declared-in-advance "
        "robustness pass re-scores both arms at 57.5 mmHg so the verdict's dependence "
        "on the anchor is MEASURED rather than assumed."
    ),
    "C3_the_median_statistic_is_mesh_sensitive": (
        "results/oxygen_joint_grid_test.md Sec.1.4 shows the median-pO2 statistic "
        "crossing the Lally value as dx changes (2.19 -> 9.00 mmHg over dx = 30 -> 8 "
        "um for the jointly calibrated configuration). dx is therefore FROZEN at 12.0 "
        "um -- the grid-converged production spacing every published results/ number "
        "is computed at -- BEFORE looking at any output, and is not varied."
    ),
    "C4_cell_weighting_vs_needle_weighting": (
        "An Eppendorf needle samples readings along a 22 mm track through space "
        "(Rampling 1994: 'Six different tracks, each with a path length of 22 mm, "
        "were recorded per patient representing 192 readings'). The engine's per-cell "
        "distribution weights by CELL DENSITY instead. These are different sampling "
        "measures. Both are reported -- per-cell AND area-weighted over the tissue "
        "grid -- rather than choosing the more favourable one."
    ),
}

# Declared in advance: the alternative mmHg anchor the robustness pass uses.
ROBUSTNESS_PREF_MMHG = 57.5  # midpoint of the contested 55-60 mmHg capillary anchor

# ---------------------------------------------------------------------------
# MEASURED ANCHORS -- every one fetched by me in this pass from the primary
# PubMed record and quoted verbatim. Nothing here is inherited from another
# results/ document without re-fetching.
# ---------------------------------------------------------------------------
ANCHORS = {
    "A1_lally2006_cohort_median": {
        "pmid": "17207315",
        "url": "https://pubmed.ncbi.nlm.nih.gov/17207315/",
        "citation": (
            "Lally BE, Rockwell S, Fischer DB, Collingridge DR, Piepmeier JM, "
            "Knisely JP. Cancer J 2006;12(6):461-6."
        ),
        "quote": (
            "Tissue oxygen tension was measured intraoperatively using an Eppendorf "
            "pO2 Histograph. ... for analysis we included 13 patients with high "
            "grade gliomas (HGG) and 10 with low grade gliomas (LGG). The median "
            "tumor oxygen pressure for the entire cohort was 5.1 mmHg."
        ),
        "value_mmHg": 5.1,
        "IMPORTANT_CORRECTION": (
            "The abstract dichotomises 'the median pO2' per patient and stratifies "
            "'patients ... with a median tumor oxygen tensions either above or "
            "below' the cut-points. 5.1 mmHg is therefore a COHORT-LEVEL median OF "
            "PER-PATIENT MEDIANS, not the median of a pooled per-reading histogram. "
            "That is a different object from the engine's per-cell distribution "
            "over one section, and this script does not pretend otherwise."
        ),
        "tissue": "human patient glioma, intraoperative",
    },
    "A2_collingridge1999_tail_exists": {
        "pmid": "10665789",
        "url": "https://pubmed.ncbi.nlm.nih.gov/10665789/",
        "citation": (
            "Collingridge DR, Piepmeier JM, Rockwell S, Knisely JP. "
            "Radiother Oncol 1999;53(2):127-31."
        ),
        "quote": (
            "This study quantifies the spatial distribution of pO2 in glioma and in "
            "the surrounding brain tissue. Both glioma and peritumoural brain "
            "contain regions at oxygen tensions less than 2.5 mmHg."
        ),
        "constraint": "P(pO2 < 2.5 mmHg) > 0 in human glioma tissue.",
        "tissue": "human patient glioma, intraoperative",
        "note": (
            "The published abstract gives NO median, NO n and NO percentage. It "
            "supports a TAIL-EXISTENCE constraint and nothing stronger."
        ),
    },
    "A3_nozue1997_interlab_floor": {
        "pmid": "9290690",
        "url": "https://pubmed.ncbi.nlm.nih.gov/9290690/",
        "citation": (
            "Nozue M, Lee I, Yuan F, ... Jain RK, Suit HD. J Surg Oncol 1997;66(1):30-8."
        ),
        "quote": (
            "the median values of tumor pO2 varied from 1.5 mmHg to 5.6 mmHg among "
            "the laboratories, and pO2 of normal subcutis also varied from 28 mmHg "
            "to 38 mmHg. ... Caution needs to be exercised in using the absolute, "
            "median, or distribution of pO2 values measured by the Eppendorf "
            "'Histograph' to compare the data between laboratories"
        ),
        "band_mmHg": [1.5, 5.6],
        "tissue": "FSaII murine fibrosarcoma -- NOT glioma, NOT human",
        "role": (
            "Used ONLY as the measurement-noise floor: the span a single median may "
            "occupy for an identical protocol on identical tissue across six labs. "
            "It brackets Lally's 5.1 mmHg, which is why it is the acceptance band."
        ),
    },
    "A4_xenograft_tail_fractions_SUPPORTING_ONLY": {
        "entries": [
            {
                "pmid": "15093922",
                "url": "https://pubmed.ncbi.nlm.nih.gov/15093922/",
                "citation": "Pinel S, et al. Int J Radiat Oncol Biol Phys 2004;59(1):250-9.",
                "quote": (
                    "Percentage of pO(2) values <2.5 mm Hg was reduced in "
                    "rHuEPO-treated tumors as compared with control groups (37.1 "
                    "+/- 19.1% vs. 58.5 +/- 27.0%; p = 0.009 for GBM Nan1; 81.6 "
                    "+/- 13.4% vs. 91.5 +/- 8.3%; p = 0.035 for U87)."
                ),
                "control_frac_below_2p5": [0.585, 0.915],
            },
            {
                "pmid": "9783885",
                "url": "https://pubmed.ncbi.nlm.nih.gov/9783885/",
                "citation": "Stuben G, et al. Radiother Oncol 1998;48(2):143-8.",
                "quote": (
                    "All three xenografts were found to be poorly oxygenated (about "
                    "80% of all pO2 values were < or =2.5 mmHg)."
                ),
                "control_frac_below_2p5": [0.80],
            },
        ],
        "tissue": "human glioma cell lines XENOGRAFTED in nude mice",
        "role": (
            "SUPPORTING ONLY, never decisive. The host vasculature is murine and "
            "the tumour is a cell line, not a patient's tumour microenvironment. "
            "These are the ONLY published per-reading fraction-below-threshold "
            "numbers for glioma I could find; no such number exists in the open "
            "abstracts for human PATIENT glioma."
        ),
    },
}

# ---------------------------------------------------------------------------
# PRE-REGISTRATION -- frozen before any solve. A failed threshold is reported
# FAILED; no threshold below may be moved after seeing a number.
# ---------------------------------------------------------------------------
PRE_REGISTRATION = {
    "section": "reconstruct_real('xenium', dx_um=12.0) -- the flagship GBM Xenium "
               "section at the grid-converged production spacing",
    "runspec": {
        "arm_A": "OxygenConfig()  [sink_form='linear', the shipping default]",
        "arm_B": "OxygenConfig(sink_form='michaelis_menten')  [K=1/35, the shipped default K]",
        "held_identical": (
            "diffusivity, supply, consumption, hypoxia_threshold, mm_half_saturation, "
            "mm_max_iter, mm_tol, the reconstruction object, the domain, the vascular "
            "field, the cell positions and the tumour mask. The two OxygenConfig "
            "dataclasses must differ in EXACTLY ONE field: sink_form. This is "
            "asserted at runtime, not claimed."
        ),
        "no_parameter_is_tuned": (
            "consumption stays at its shipped 0.02 and K at its shipped 1/35 in both "
            "arms. Nothing is fitted to any anchor. If an arm misses the band it "
            "misses it."
        ),
    },
    "quantities": [
        "per-tumour-cell pO2 in mmHg (O sampled at each tumour cell's grid index, "
        "exactly as physics.oxygen.label_hypoxia samples it, then to_mmhg)",
        "the full quantile ladder: min, p1, p5, p10, p25, p50, p75, p90, p95, p99, max",
        "P(pO2 < 2.5), P(pO2 < 5), P(pO2 < 7 == the engine's own hypoxia band), P(pO2 < 10)",
        "the same ladder over the tissue-masked GRID (area-weighted), because an "
        "Eppendorf needle track samples space, not cells -- reported alongside the "
        "per-cell ladder rather than instead of it",
    ],
    "thresholds": {
        "T1_median_band": {
            "rule": "An arm PASSES iff its per-tumour-cell median pO2 lies within "
                    "[1.5, 5.6] mmHg.",
            "provenance": "Nozue 1997's measured inter-laboratory span of median "
                          "tumour pO2 under an identical protocol (A3). It contains "
                          "Lally 2006's measured human-glioma cohort median of 5.1 "
                          "mmHg (A1). The band is the measurement's own noise floor, "
                          "so an arm outside it is outside by more than laboratory "
                          "variability can absorb.",
            "direction_is_recorded": True,
        },
        "T2_tail_exists": {
            "rule": "An arm PASSES iff P(pO2 < 2.5 mmHg) > 0 among tumour cells.",
            "provenance": "Collingridge 1999 (A2): human glioma contains regions "
                          "below 2.5 mmHg. Deliberately the WEAKEST possible form of "
                          "that constraint -- a single cell suffices -- so that "
                          "failing it is unambiguous and cannot be argued down.",
        },
        "T3_form_explains_the_gap": {
            "rule": "The functional form EXPLAINS the missing hypoxia iff arm A "
                    "(linear) fails T1 or T2, AND arm B (michaelis_menten) passes "
                    "BOTH T1 and T2.",
        },
        "T4_null_no_effect": {
            "rule": "The two forms are INDISTINGUISHABLE iff |median_A - median_B| < "
                    "0.5 mmHg AND |P_A(<2.5) - P_B(<2.5)| < 0.05.",
            "provenance": "0.5 mmHg is ~12% of the width of the T1 acceptance band "
                          "(4.1 mmHg); a median difference smaller than that is "
                          "smaller than the measurement can resolve. Declared so "
                          "that a NULL is as interpretable as a hit.",
        },
        "T5_refuted_by_overshoot": {
            "rule": "The functional form DOES NOT explain the gap, by OVERSHOOT, iff "
                    "arm B passes T2 but fails T1 LOW (median < 1.5 mmHg) while arm A "
                    "fails T1 HIGH (median > 5.6 mmHg) -- i.e. the two forms BRACKET "
                    "the measurement instead of one of them matching it.",
            "why_declared_in_advance": (
                "This is the outcome most likely to be rationalised after the fact as "
                "'the MM form is directionally right'. Naming it in advance as a "
                "REFUTATION, with its own label, removes that escape. If T5 fires, it "
                "leads the report."
            ),
        },
    },
    "what_would_be_reported_as_a_failure_of_this_experiment_itself": (
        "If the Newton solve fails to converge, or if the degeneracy control "
        "(CONTROL_2 below) does not reproduce the linear arm, the A/B is void and "
        "is reported void rather than interpreted."
    ),
    "controls": {
        "CONTROL_1_determinism": "Arm A solved twice must be bit-identical.",
        "CONTROL_2_degeneracy": (
            "michaelis_menten with K >> 1 and consumption scaled by the same K "
            "degenerates analytically to the linear sink, since c*K*O/(O+K) -> c*O as "
            "K -> inf. Running it must reproduce arm A to tight tolerance. This "
            "proves the difference between the arms is the FUNCTIONAL FORM and not a "
            "solver artefact."
        ),
        "CONTROL_3_input_identity": (
            "SHA-256 of every array fed to the solver (vascular density, cell x/y/"
            "type, the derived tumour areal density, the domain parameters) taken "
            "before arm A, between the arms and after arm B. All three must match."
        ),
    },
}

QUANTILES = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]
THRESHOLDS_MMHG = [2.5, 5.0, 7.0, 10.0]
NOZUE_BAND = (1.5, 5.6)


def _h(*arrays) -> str:
    m = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        m.update(str(a.dtype).encode())
        m.update(str(a.shape).encode())
        m.update(a.tobytes())
    return m.hexdigest()[:16]


def input_fingerprint(micro) -> str:
    """SHA-256 over everything solve_oxygen reads. CONTROL_3."""
    from tmesim.physics.oxygen import _tumor_density
    d = micro.domain
    return _h(
        np.array([d.nx, d.ny], dtype=np.int64),
        np.array([float(d.dx)], dtype=np.float64),
        np.asarray(micro.fields.vascular_density, dtype=np.float64),
        np.asarray(micro.cells.x, dtype=np.float64),
        np.asarray(micro.cells.y, dtype=np.float64),
        np.asarray(micro.cells.cell_type, dtype=np.int64),
        _tumor_density(micro),
    )


def describe(po2: np.ndarray) -> dict:
    """Full distributional description of a pO2 sample, in mmHg."""
    q = np.percentile(po2, QUANTILES)
    out = {
        "n": int(po2.size),
        "mean_mmHg": float(po2.mean()),
        "std_mmHg": float(po2.std()),
        "quantiles_mmHg": {f"p{p}": float(v) for p, v in zip(QUANTILES, q)},
        "min_mmHg": float(po2.min()),
        "max_mmHg": float(po2.max()),
    }
    out["median_mmHg"] = out["quantiles_mmHg"]["p50"]
    out["frac_below"] = {
        f"{t}": float((po2 < t).mean()) for t in THRESHOLDS_MMHG
    }
    # Empirical CDF on a fixed mmHg grid, so the SHAPE can be compared by eye
    # and by a reader who has a measured histogram we do not.
    grid = [0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 7.0, 10.0, 14.0, 20.0, 28.0, 35.0]
    out["ecdf"] = {f"{g}": float((po2 < g).mean()) for g in grid}
    return out


def score(arm_name: str, d: dict) -> dict:
    med = d["median_mmHg"]
    tail = d["frac_below"]["2.5"]
    lo, hi = NOZUE_BAND
    t1 = lo <= med <= hi
    direction = "in_band" if t1 else ("LOW" if med < lo else "HIGH")
    t2 = tail > 0.0
    return {
        "arm": arm_name,
        "median_mmHg": med,
        "T1_median_in_1.5_5.6": bool(t1),
        "T1_direction": direction,
        "T1_verdict": "PASS" if t1 else f"FAIL ({direction})",
        "frac_below_2.5_mmHg": tail,
        "T2_tail_exists": bool(t2),
        "T2_verdict": "PASS" if t2 else "FAIL",
    }


def main() -> None:
    print("=" * 78)
    a0 = calibrate_alpha()
    print(f"START calibrate_alpha() = {a0!r}")
    assert a0 == ALPHA_PINNED, f"CALIBRATION DRIFT at START: {a0!r} != {ALPHA_PINNED!r}"
    print("ASSERT OK")
    print("=" * 78)

    print("\nPRE-REGISTRATION (frozen in this file before any solve):")
    print(f"  T1 median band      : [{NOZUE_BAND[0]}, {NOZUE_BAND[1]}] mmHg (Nozue 1997)")
    print(f"  T2 tail existence   : P(pO2 < 2.5 mmHg) > 0 (Collingridge 1999)")
    print(f"  T3 form explains gap: A fails, B passes both")
    print(f"  T4 null             : |dmedian| < 0.5 mmHg and |dtail| < 0.05")
    print(f"  T5 refuted-by-overshoot: B fails LOW while A fails HIGH")
    print(f"  REFERENCE_PO2_MMHG  = {REFERENCE_PO2_MMHG}  (O=1 anchor)")

    print("\nReconstructing flagship section (xenium, dx_um=12.0)...")
    t0 = time.time()
    micro = reconstruct_real("xenium", dx_um=12.0)
    print(f"  done in {time.time()-t0:.1f}s  grid={micro.domain.shape}  "
          f"dx={micro.domain.dx}  n_cells={micro.cells.n}")

    tum = micro.cells.cell_type == int(CellType.TUMOR)
    n_tumor = int(tum.sum())
    col, row = micro.domain.to_index(micro.cells.x[tum], micro.cells.y[tum])
    tmask = np.asarray(micro.domain.tissue_mask, dtype=bool)
    print(f"  n_tumor_cells={n_tumor}  tissue grid points={int(tmask.sum())}")

    fp_before = input_fingerprint(micro)
    print(f"\nCONTROL_3 input fingerprint BEFORE arm A: {fp_before}")

    cfg_A = OxygenConfig()
    cfg_B = OxygenConfig(sink_form="michaelis_menten")

    dA, dB = dataclasses.asdict(cfg_A), dataclasses.asdict(cfg_B)
    diff = {k: (dA[k], dB[k]) for k in dA if dA[k] != dB[k]}
    print(f"\nCONFIG DIFF (must be exactly one field): {diff}")
    assert set(diff) == {"sink_form"}, f"arms differ in more than sink_form: {diff}"
    print("  ASSERT OK -- the two arms differ in exactly one field.")
    print(f"  shared: consumption={cfg_A.consumption}  supply={cfg_A.supply}  "
          f"D={cfg_A.diffusivity}  hypoxia_threshold={cfg_A.hypoxia_threshold}  "
          f"K={cfg_A.mm_half_saturation!r}")

    # ---------------- ARM A ----------------
    print("\n--- ARM A: sink_form='linear' (shipping default) ---")
    t0 = time.time()
    oxy_A = solve_oxygen(micro, cfg_A)
    tA = time.time() - t0
    print(f"  solved in {tA:.2f}s")

    # CONTROL_1
    oxy_A2 = solve_oxygen(micro, cfg_A)
    det = bool(np.array_equal(oxy_A, oxy_A2))
    print(f"  CONTROL_1 determinism (bit-identical re-solve): {det}")
    assert det, "arm A is not deterministic"

    fp_mid = input_fingerprint(micro)
    print(f"  CONTROL_3 fingerprint BETWEEN arms: {fp_mid}")

    # ---------------- ARM B ----------------
    print("\n--- ARM B: sink_form='michaelis_menten' (K=1/35) ---")
    t0 = time.time()
    oxy_B = solve_oxygen(micro, cfg_B)
    tB = time.time() - t0
    print(f"  solved in {tB:.2f}s")

    fp_after = input_fingerprint(micro)
    print(f"  CONTROL_3 fingerprint AFTER arm B : {fp_after}")
    identical = (fp_before == fp_mid == fp_after)
    print(f"  CONTROL_3 all three identical: {identical}")
    assert identical, "solver inputs changed between arms -- A/B is not controlled"

    # ---------------- CONTROL_2: degeneracy ----------------
    print("\n--- CONTROL_2: MM degenerates to linear as K -> inf ---")
    K_big = 1.0e4
    cfg_C = OxygenConfig(sink_form="michaelis_menten",
                         consumption=cfg_A.consumption * K_big,
                         mm_half_saturation=K_big)
    oxy_C = solve_oxygen(micro, cfg_C)
    deg_err = float(np.max(np.abs(oxy_C - oxy_A)))
    print(f"  K={K_big:g}, consumption={cfg_A.consumption*K_big:g}  "
          f"max|O_MM - O_linear| = {deg_err:.3e}")
    degen_ok = deg_err < 1e-3
    print(f"  CONTROL_2 {'PASS' if degen_ok else 'FAIL'} "
          f"(the arms differ by the FORM, not by the solver)")

    # ---------------- distributions ----------------
    po2_A_cells = to_mmhg(oxy_A[row, col])
    po2_B_cells = to_mmhg(oxy_B[row, col])
    po2_A_grid = to_mmhg(oxy_A[tmask])
    po2_B_grid = to_mmhg(oxy_B[tmask])

    res = {
        "linear": {"cells": describe(po2_A_cells), "grid": describe(po2_A_grid),
                   "solve_time_s": tA},
        "michaelis_menten": {"cells": describe(po2_B_cells), "grid": describe(po2_B_grid),
                             "solve_time_s": tB},
    }

    def show(tag: str, d: dict) -> None:
        q = d["quantiles_mmHg"]
        print(f"\n  [{tag}]  n={d['n']}")
        print(f"    min {d['min_mmHg']:.4f} | p1 {q['p1']:.4f} | p5 {q['p5']:.4f} | "
              f"p10 {q['p10']:.4f} | p25 {q['p25']:.4f}")
        print(f"    MEDIAN {q['p50']:.4f} mmHg | p75 {q['p75']:.4f} | "
              f"p90 {q['p90']:.4f} | p95 {q['p95']:.4f} | p99 {q['p99']:.4f} | "
              f"max {d['max_mmHg']:.4f}")
        print(f"    mean {d['mean_mmHg']:.4f} +/- {d['std_mmHg']:.4f}")
        fb = d["frac_below"]
        print(f"    P(<2.5) {fb['2.5']:.4f} | P(<5) {fb['5.0']:.4f} | "
              f"P(<7) {fb['7.0']:.4f} | P(<10) {fb['10.0']:.4f}")

    print("\n" + "=" * 78)
    print("MEASURED DISTRIBUTIONS (pO2 in mmHg via REFERENCE_PO2_MMHG = 35)")
    print("=" * 78)
    show("LINEAR / tumour cells", res["linear"]["cells"])
    show("LINEAR / tissue grid", res["linear"]["grid"])
    show("MICHAELIS-MENTEN / tumour cells", res["michaelis_menten"]["cells"])
    show("MICHAELIS-MENTEN / tissue grid", res["michaelis_menten"]["grid"])

    # ---------------- reproduction check against the spent prior ----------------
    prior = PRIOR_INFORMATION["results/oxygen_sensitivity_band.md"]["known_in_advance"]
    repro = {
        "linear_hypoxic_fraction": {
            "published": prior["linear_hypoxic_fraction_at_c0.02"],
            "measured_now": res["linear"]["cells"]["frac_below"]["7.0"],
        },
        "michaelis_menten_hypoxic_fraction": {
            "published": prior["michaelis_menten_hypoxic_fraction_at_c0.02"],
            "measured_now": res["michaelis_menten"]["cells"]["frac_below"]["7.0"],
        },
        "linear_mean_O2_tumor": {
            "published": prior["linear_mean_O2_tumor_at_c0.02"],
            "measured_now": float(po2_A_cells.mean() / REFERENCE_PO2_MMHG),
        },
        "michaelis_menten_mean_O2_tumor": {
            "published": prior["michaelis_menten_mean_O2_tumor_at_c0.02"],
            "measured_now": float(po2_B_cells.mean() / REFERENCE_PO2_MMHG),
        },
        "linear_median_pO2_mmHg_joint_grid_test": {
            "published": PRIOR_INFORMATION["results/oxygen_joint_grid_test.md"
                                           ]["known_in_advance"]["shipped_linear_median_pO2_mmHg"],
            "measured_now": res["linear"]["cells"]["median_mmHg"],
        },
    }
    print("\n" + "=" * 78)
    print("REPRODUCTION CHECK vs the ALREADY-PUBLISHED numbers (NOT a new result)")
    print("=" * 78)
    for k, v in repro.items():
        print(f"  {k:48s} published {v['published']!r:12} now {v['measured_now']:.6f}")

    # ---------------- scoring ----------------
    sA = score("linear", res["linear"]["cells"])
    sB = score("michaelis_menten", res["michaelis_menten"]["cells"])

    dmed = abs(sA["median_mmHg"] - sB["median_mmHg"])
    dtail = abs(sA["frac_below_2.5_mmHg"] - sB["frac_below_2.5_mmHg"])
    t4_null = (dmed < 0.5) and (dtail < 0.05)
    t3 = ((not (sA["T1_median_in_1.5_5.6"] and sA["T2_tail_exists"]))
          and sB["T1_median_in_1.5_5.6"] and sB["T2_tail_exists"])
    t5 = (sB["T2_tail_exists"] and sB["T1_direction"] == "LOW"
          and sA["T1_direction"] == "HIGH")

    if t3:
        verdict = "T3 CONFIRMED -- the sink's functional form EXPLAINS the missing hypoxia"
    elif t5:
        verdict = ("T5 REFUTED-BY-OVERSHOOT -- the functional form does NOT explain the "
                   "gap; the two forms BRACKET the measurement")
    elif t4_null:
        verdict = "T4 NULL -- the two sink forms are indistinguishable on this section"
    else:
        verdict = "NEITHER T3 NOR T5 NOR T4 -- see the per-threshold table"

    print("\n" + "=" * 78)
    print("PRE-REGISTERED SCORING")
    print("=" * 78)
    for s in (sA, sB):
        print(f"  {s['arm']:18s} median {s['median_mmHg']:9.4f} mmHg  "
              f"T1 {s['T1_verdict']:14s}  P(<2.5)={s['frac_below_2.5_mmHg']:.4f}  "
              f"T2 {s['T2_verdict']}")
    print(f"\n  |dmedian| = {dmed:.4f} mmHg   |dtail| = {dtail:.4f}")
    print(f"  T3 (form explains gap)   : {t3}")
    print(f"  T4 (null, no effect)     : {t4_null}")
    print(f"  T5 (refuted by overshoot): {t5}")
    print(f"\n  VERDICT: {verdict}")

    # ---------------- ROBUSTNESS: does the verdict survive the contested anchor? ----
    scale = ROBUSTNESS_PREF_MMHG / REFERENCE_PO2_MMHG
    rob = {}
    print("\n" + "=" * 78)
    print(f"ROBUSTNESS (declared in advance) -- re-anchored at p_ref = "
          f"{ROBUSTNESS_PREF_MMHG} mmHg (x{scale:.3f})")
    print("=" * 78)
    for name, d in (("linear", res["linear"]["cells"]),
                    ("michaelis_menten", res["michaelis_menten"]["cells"])):
        med2 = d["median_mmHg"] * scale
        lo, hi = NOZUE_BAND
        ok = lo <= med2 <= hi
        direction = "in_band" if ok else ("LOW" if med2 < lo else "HIGH")
        # What p_ref would this arm NEED to land its median on Lally's 5.1 mmHg?
        med_dimensionless = d["median_mmHg"] / REFERENCE_PO2_MMHG
        implied = (5.1 / med_dimensionless) if med_dimensionless > 0 else float("inf")
        rob[name] = {
            "median_mmHg_at_57.5_anchor": med2,
            "T1_verdict_at_57.5_anchor": "PASS" if ok else f"FAIL ({direction})",
            "p_ref_mmHg_required_to_hit_Lally_5.1": implied,
        }
        print(f"  {name:18s} median {med2:9.4f} mmHg -> "
              f"{'PASS' if ok else 'FAIL (' + direction + ')'}")
        print(f"  {'':18s} p_ref it would NEED to hit Lally's 5.1 mmHg: "
              f"{implied:.1f} mmHg")
    print("\n  (Arterial pO2 is ~100 mmHg. A required p_ref far above that is not a "
          "\n   rescuable anchor choice -- it is a statement that the field is wrong.)")
    verdict_robust = (rob["linear"]["T1_verdict_at_57.5_anchor"].startswith("FAIL")
                      and rob["michaelis_menten"]["T1_verdict_at_57.5_anchor"].startswith("FAIL"))
    print(f"\n  Both arms still FAIL T1 under the contested anchor: {verdict_robust}")

    # ---------------- POST-HOC (NOT pre-registered) -------------------------
    # Labelled post-hoc because it was devised AFTER seeing the two distributions.
    # It is reported as a diagnostic, never as a registered test.
    print("\n" + "=" * 78)
    print("POST-HOC DIAGNOSTIC -- NOT PRE-REGISTERED, devised after seeing the output")
    print("=" * 78)
    cA, cB = res["linear"]["cells"], res["michaelis_menten"]["cells"]
    # (a) Support disjointness.
    gap_lo, gap_hi = cB["max_mmHg"], cA["min_mmHg"]
    disjoint = gap_lo < gap_hi
    lally = ANCHORS["A1_lally2006_cohort_median"]["value_mmHg"]
    in_gap = gap_lo < lally < gap_hi
    print(f"  (a) SUPPORT DISJOINTNESS")
    print(f"      michaelis_menten MAX = {gap_lo:.4f} mmHg")
    print(f"      linear           MIN = {gap_hi:.4f} mmHg")
    print(f"      supports disjoint: {disjoint}   gap = "
          f"[{gap_lo:.4f}, {gap_hi:.4f}] mmHg")
    print(f"      Lally's measured cohort median ({lally} mmHg) lies IN that gap: "
          f"{in_gap}")
    if in_gap:
        print(f"      => at the shipped consumption, NEITHER sink form places a "
              f"single\n         tumour cell at the measured median.")
    # (b) Anchor-INVARIANT shape test. to_mmhg is a pure rescale, so min/median
    #     is invariant under ANY choice of REFERENCE_PO2_MMHG.
    print(f"\n  (b) ANCHOR-INVARIANT DYNAMIC RANGE  R = min / median")
    print(f"      Lally (median 5.1) + Collingridge (regions < 2.5) together require")
    print(f"      the tissue to reach 2.5/5.1 = {2.5/lally:.4f} of its own median.")
    shape = {}
    for name, d in (("linear", cA), ("michaelis_menten", cB)):
        R = d["min_mmHg"] / d["median_mmHg"]
        ok = R <= (2.5 / lally)
        shape[name] = {"R_min_over_median": R, "can_reach_tail_at_any_anchor": bool(ok)}
        print(f"      {name:18s} R = {R:.4f}  -> can reach a sub-2.5 mmHg region at "
              f"SOME anchor: {ok}")
    print("      This is the one statement in this file that no choice of")
    print("      REFERENCE_PO2_MMHG can change.")

    posthoc = {
        "IS_POST_HOC_NOT_PREREGISTERED": True,
        "support_disjointness": {
            "michaelis_menten_max_mmHg": gap_lo,
            "linear_min_mmHg": gap_hi,
            "supports_disjoint": bool(disjoint),
            "gap_mmHg": [gap_lo, gap_hi],
            "lally_median_mmHg": lally,
            "lally_median_lies_in_gap": bool(in_gap),
        },
        "anchor_invariant_shape": {
            "required_R": 2.5 / lally,
            "basis": "Lally 2006 median 5.1 mmHg with Collingridge 1999 regions "
                     "below 2.5 mmHg; R = min/median is invariant under any pure "
                     "rescale of REFERENCE_PO2_MMHG.",
            "per_arm": shape,
        },
    }

    a1 = calibrate_alpha()
    print("\n" + "=" * 78)
    print(f"END calibrate_alpha() = {a1!r}")
    assert a1 == ALPHA_PINNED, f"CALIBRATION DRIFT at END: {a1!r} != {ALPHA_PINNED!r}"
    print("ASSERT OK")
    print("=" * 78)

    payload = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "script": "experiments/oxygen_sink_form.py",
        "git_head": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO,
                                   capture_output=True, text=True).stdout.strip(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "calibrate_alpha_start": a0,
        "calibrate_alpha_end": a1,
        "REFERENCE_PO2_MMHG": REFERENCE_PO2_MMHG,
        "prior_information": PRIOR_INFORMATION,
        "confounds": CONFOUNDS,
        "anchors": ANCHORS,
        "pre_registration": PRE_REGISTRATION,
        "section": {
            "call": "reconstruct_real('xenium', dx_um=12.0)",
            "grid": list(micro.domain.shape),
            "dx_um": float(micro.domain.dx),
            "n_cells": int(micro.cells.n),
            "n_tumor_cells": n_tumor,
            "n_tissue_grid_points": int(tmask.sum()),
        },
        "config_diff": {k: {"linear": v[0], "michaelis_menten": v[1]}
                        for k, v in diff.items()},
        "config_shared": {k: dA[k] for k in dA if k not in diff},
        "controls": {
            "CONTROL_1_determinism": det,
            "CONTROL_2_degeneracy_max_abs_err": deg_err,
            "CONTROL_2_pass": degen_ok,
            "CONTROL_3_fingerprints": {"before": fp_before, "between": fp_mid,
                                       "after": fp_after, "identical": identical},
        },
        "distributions": res,
        "reproduction_check": repro,
        "scoring": {
            "linear": sA, "michaelis_menten": sB,
            "delta_median_mmHg": dmed, "delta_tail_frac": dtail,
            "T3_form_explains_gap": bool(t3),
            "T4_null": bool(t4_null),
            "T5_refuted_by_overshoot": bool(t5),
            "verdict": verdict,
        },
        "robustness_alternative_anchor": {
            "p_ref_mmHg": ROBUSTNESS_PREF_MMHG,
            "scale_vs_shipped": scale,
            "per_arm": rob,
            "both_arms_still_fail_T1": bool(verdict_robust),
        },
        "post_hoc_diagnostics": posthoc,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {OUT_JSON.relative_to(REPO)}")


if __name__ == "__main__":
    main()
