"""DOES THE MAP REPRODUCE A MEASURED PENETRATION DIFFERENCE?
Doxorubicin / epirubicin / daunorubicin / mitoxantrone against Kyle 2004 Table 1.

THE QUESTION
------------
Codon's surviving claim is a **map**: given one real tumour tissue section and one
compound, where does the drug reach and how much of the tumour does it never
reach. Within a single map the uncertain per-compound absolutes cancel, so "this
region got 10x less than that region" needs no absolute calibration.

**Nobody has ever checked that the map is right.** This script is that check.

THE MEASUREMENT
---------------
Kyle AH, Huxham LA, Chiam ASJ, Sim DH, Minchinton AI. "Direct assessment of drug
penetration into tissue using a novel application of three-dimensional cell
culture." Cancer Res 2004;64(17):6304-9. PMID 15342419,
doi:10.1158/0008-5472.CAN-04-1099.

**Table 1**, verbatim header: *"Drug exposure to the near side relative to the far
side of MCCs as estimated from the effect-based assay or drug derived
fluorescence"*:

    drug   effect-based assay   drug fluorescence
    DOX          ~12               11 +- 5 (SD)
    EPI          ~10               12 +- 3 (SD)
    DAU          ~30               50 +- 15 (SD)
    MIT          >30               (not measurable)

HCT-116 multilayered cell culture, ~150 um thick (140 +- 15 um SD when frozen
immediately after exposure), 1 h exposure to ONE side with the other side sealed.

WHY THIS SOURCE AND NOT THE OTHERS — the two gates
--------------------------------------------------
An earlier candidate (Trédan 2009, mitoxantrone vs AQ4N) was **rejected**, and the
two reasons are now standing criteria applied to every candidate:

**GATE 1 — DOSE PARITY. Kyle 2004 PASSES, verbatim:** *"The ratio of fluorescence
on each side is calculated from data shown in Fig. 6 using the MCCs frozen
immediately after exposure to **10 μmol/L** drug over 1 hour."* All three
fluorescence-measured drugs at the SAME 10 umol/L, and the effect-based column is
likewise referenced to *"the effect seen on the far side at 10 μmol/L drug"*.
Trédan compared 1 vs 10 mg/mL in vitro and 40 vs 200 mg/kg in vivo, doses chosen
*"to facilitate detection"* — a 5-10x asymmetry that no transport comparison
survives.

**GATE 2 — IS THE MEASURED QUANTITY ACTUALLY TRANSPORT? Kyle 2004 PASSES.** The
fluorescence column is a spatially resolved profile of the **parent molecule**:
*"Only DOX, EPI, and DAU, which all possess the same fluorophore, could be
examined. Emission intensities of the three drugs in aqueous solution were within
10% of each other."* One species, one fluorophore, cross-calibrated to 10%.
Trédan's AQ4N signal was AQ4N **plus** its reduced metabolite AQ4, which
*"cannot be separated"*, with hypoxic bioreductive trapping inflating apparent
depth — transport, activation and binding inseparable in one number.

**What Kyle 2004's number does include, and this is disclosed everywhere it is
quoted:** cellular sequestration as well as diffusion. The authors are explicit —
*"DOX appears to accumulate the least over the 1-hour exposure with EPI and DAU
reaching ~1.5 and ~7 times higher levels, respectively, on the near sides."*
That is legitimate here because the engine models exactly this composite: a
transport coefficient ``D_eff`` against a loss/uptake sink ``k``, giving
``L = sqrt(D_eff/k)``. It would NOT be legitimate against a pure-diffusion model.

WHY THIS PAIR IS A RISKY TEST AND NOT A RIGGED ONE
--------------------------------------------------
Daunorubicin is doxorubicin minus the 14-hydroxyl:

    doxorubicin    PubChem CID 31703   C27H29NO11  MW 543.5  HBD 6  XLogP 1.3
    daunorubicin   PubChem CID 30323   C27H29NO10  MW 527.5  HBD 5  XLogP 1.8

Daunorubicin is **more lipophilic and has one fewer hydrogen-bond donor**. The
shipping two-route law (``results/two_route_promotion.md``) opens the
transcellular route as a function of exactly those two descriptors, with an HBD
decrement of 1.34 log units per donor. So the transport law, on its own, predicts
daunorubicin's membrane is the **more open** of the two and its ``D_eff`` the
**larger** — i.e. that daunorubicin should penetrate **further**.

**The measurement says the opposite, by 4.5-fold.** Whether the engine's uptake
sink ``k`` overturns its own diffusivity term hard enough to recover the measured
direction is genuinely unknown before the run, and it is the content of this test.
The sink is a structure-only prior that ``results/alpha_identifiability.md``
proves is **unidentified** for every drug except doxorubicin.

And the pair is nearly matched in mass (543.5 vs 527.5, 3% apart), so the
superseded molecular-weight-only law of
``results/diffusivity_falsification.md`` is almost flat between them: it predicts
a 1.01x separation against a measured 4.5x.

A FREE POSITIVE CONTROL
-----------------------
Epirubicin is doxorubicin's 4'-epimer — same formula, same mass, same graph, only
stereochemistry differs. The engine's descriptors cannot see stereochemistry, so
it **must** predict them identical. Kyle Table 1 measures them as identical too
(11 +- 5 vs 12 +- 3; ~12 vs ~10). A model that separated them would be predicting
a difference the measurement excludes. This costs one extra solve and it checks
the readout in the direction where the right answer is "no difference".

IN-SAMPLE DISCLOSURE — READ BEFORE QUOTING ANY NUMBER
-----------------------------------------------------
**Doxorubicin is the drug ``calibrate_alpha()`` is fit on.** Its absolute
penetration length is pinned to 64.92 um by construction (Primeau 2005,
PMID 16361566) and cannot fail. ``scorecard.py::_run_perivascular_penetration``
already declines to score it for this reason and this file inherits that
judgement: **doxorubicin's absolute agreement with Kyle 2004 is NOT evidence.**

What IS out-of-sample is the **ratio**. Nothing was ever fit to daunorubicin,
epirubicin or mitoxantrone; alpha is frozen; so ``ratio(DAU)/ratio(DOX)`` is
entirely a prediction about daunorubicin, in the same structure the repository
already accepts for mitoxantrone and topotecan in ``pk_multidrug.run_penetration``.

Daunorubicin is new to this repository. None of the four compounds is among the
71 Pruijn 2008 HT-29 compounds that ``results/transport_anchor_prereg.md`` §0.1
declares **spent**; stage ``compounds`` checks that against the table rather than
trusting prose.

WHAT IS CHANGED
---------------
Nothing. No default, no constant, no shipped parameter. ``calibrate_alpha()`` is
asserted live at ``22.223477711172563`` in every stage. Nothing under ``tmesim/``
is written to; this script drives it and reuses
``tmesim.validation.vessel_profile``'s profile/fit machinery unmodified so the
readout cannot drift from the one the published perivascular-gradient report uses.

THE BAR IS FIXED BEFORE THE MEASUREMENT
---------------------------------------
``--stage prereg`` writes ``results/sweep/kyle2004_penetration_prereg.json`` the
first time and **asserts equality** on every later run. Every other stage refuses
to start until that file exists.

Run:
    ./.venv/bin/python experiments/isomer_penetration_map.py --stage prereg
    ./.venv/bin/python experiments/isomer_penetration_map.py --stage compounds
    ./.venv/bin/python experiments/isomer_penetration_map.py --stage mcc
    ./.venv/bin/python experiments/isomer_penetration_map.py --stage all
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from tmesim.core import (
    CellPopulation,
    Domain,
    Microenvironment,
    PDEConfig,
    SimConfig,
    TissueFields,
)
from tmesim.molecular import descriptors as _descriptors
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.molecular.compile_drug import build_profile, molecule_to_drugparams
from tmesim.physics.diffusion import ALPHA_TORTUOSITY
from tmesim.sim.engine import simulate_response
from tmesim.validation import vessel_profile as VP
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams

# --------------------------------------------------------------------------- #
# 0. Constants
# --------------------------------------------------------------------------- #

#: The shipped calibration constant. Asserted live in every stage.
SHIPPED_ALPHA = 22.223477711172563

DOX, EPI, DAU, MIT = "doxorubicin", "epirubicin", "daunorubicin", "mitoxantrone"
ORDER: Tuple[str, ...] = (DOX, EPI, DAU, MIT)

#: Structures. Three are REUSED verbatim from the repository so this file cannot
#: introduce a second, silently-diverging definition of a molecule the codebase
#: already carries. Only daunorubicin is new, and its provenance is recorded.
SMILES: Dict[str, str] = {
    DOX: ("COc1cccc2c1C(=O)c1c(O)c3c(c(O)c1C2=O)CC(O)(C(=O)CO)CC3OC1CC(N)"
          "C(O)C(C)O1"),                       # tmesim/molecular/library.py
    EPI: ("CC1C(C(CC(O1)OC2CC(CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)"
          "(C(=O)CO)O)N)O"),                   # tmesim/molecular/breast_library.py
    DAU: ("CC1C(C(CC(O1)OC2CC(CC3=C2C(=C4C(=C3O)C(=O)C5=C(C4=O)C(=CC=C5)OC)O)"
          "(C(=O)C)O)N)O"),                    # PubChem CID 30323 (see PROVENANCE)
    MIT: "OCCNCCNc1ccc(NCCNCCO)c2c1C(=O)c1c(O)ccc(O)c1C2=O",  # calibration.py
}

PROVENANCE: Dict[str, str] = {
    DOX: "REUSED verbatim from tmesim/molecular/library.py::GBM_DRUGS.",
    EPI: "REUSED verbatim from tmesim/molecular/breast_library.py::BREAST_DRUGS.",
    DAU: ("NEW to this repository. PubChem PUG-REST property endpoint, name "
          "lookup 'daunorubicin' -> CID 30323, ConnectivitySMILES. Returned "
          "MolecularFormula C27H29NO10, MolecularWeight 527.5, XLogP 1.8, "
          "HBondDonorCount 5, HBondAcceptorCount 11, TPSA 186."),
    MIT: ("REUSED verbatim from tmesim/validation/calibration.py::"
          "MITOXANTRONE_SMILES."),
}

#: Kyle 2004 Table 1, transcribed. `fluor` is the primary column (a direct,
#: equimolar, single-fluorophore measurement of the parent drug); `effect` is the
#: independent effect-based estimate. `None` where the paper reports none.
KYLE_TABLE1: Dict[str, Dict] = {
    DOX: {"effect": 12.0, "effect_kind": "~", "fluor": 11.0, "fluor_sd": 5.0},
    EPI: {"effect": 10.0, "effect_kind": "~", "fluor": 12.0, "fluor_sd": 3.0},
    DAU: {"effect": 30.0, "effect_kind": "~", "fluor": 50.0, "fluor_sd": 15.0},
    MIT: {"effect": 30.0, "effect_kind": ">", "fluor": None, "fluor_sd": None},
}

MEASUREMENT = {
    "citation": ("Kyle AH, Huxham LA, Chiam ASJ, Sim DH, Minchinton AI. Direct "
                 "assessment of drug penetration into tissue using a novel "
                 "application of three-dimensional cell culture. Cancer Res "
                 "2004;64(17):6304-9."),
    "pmid": "15342419",
    "doi": "10.1158/0008-5472.CAN-04-1099",
    "table": "Table 1",
    "quantity": ("ratio of near-side to far-side drug exposure across a "
                 "multilayered cell culture"),
    "system": ("HCT-116 human colorectal carcinoma multilayered cell culture "
               "(MCC), one-sided 1 h exposure with the far side sealed"),
    "thickness_um": 140.0,
    "thickness_sd_um": 15.0,
    "thickness_quote": ("MCCs frozen immediately after exposure measured "
                        "140 +- 15 (SD) um in thickness"),
    "edge_window_um": 30.0,
    "edge_window_quote": ("the first 30 um of tissue on the near side and far "
                          "side of MCCs"),
    "exposure_h": 1.0,
    "concentration_uM": 10.0,
    "gate1_dose_parity_quote": (
        "The ratio of fluorescence on each side is calculated from data shown "
        "in Fig. 6 using the MCCs frozen immediately after exposure to "
        "10 umol/L drug over 1 hour."),
    "gate2_species_quote": (
        "Only DOX, EPI, and DAU, which all possess the same fluorophore, could "
        "be examined. Emission intensities of the three drugs in aqueous "
        "solution were within 10% of each other."),
    "sequestration_quote": (
        "DOX appears to accumulate the least over the 1-hour exposure with EPI "
        "and DAU reaching ~1.5 and ~7 times higher levels, respectively, on "
        "the near sides."),
    "near_side_accumulation_vs_dox": {DOX: 1.0, EPI: 1.5, DAU: 7.0, MIT: None},
    "far_side_indistinguishable_quote": (
        "Differences in the effects of the drugs on the far sides of the "
        "cultures were found to be statistically indistinguishable, P > 0.05, "
        "except at the 10 umol/L x hour exposure, P < 0.005 (ANOVA), where the "
        "effect of EPI was significantly greater than DOX."),
    "near_side_dominance_quote": (
        "in both the effect-based assay and direct fluorescence assay, the "
        "drugs performed similarly on the far sides of the MCCs regardless of "
        "their differing lipophilicity, and it was only on the near sides that "
        "their effects could be differentiated."),
    "n_underlying": "n = 4-6 (Fig. 3, SE bars); n for the fluorescence column "
                    "is not stated in the paper",
    "table1": KYLE_TABLE1,
}

#: Sources considered and REJECTED, with the gate each failed. Recorded so a
#: later reader can see the choice was made on stated criteria and not on which
#: dataset flattered the model.
REJECTED_SOURCES = {
    "Tredan 2009 Cancer Res 69:940 (PMID 19176397), mitoxantrone vs AQ4N": {
        "gate1_dose_parity": ("FAIL. 1 vs 10 mg/mL in the multicellular "
                              "culture and 40 vs 200 mg/kg in vivo, chosen "
                              "verbatim 'to facilitate detection' — a 5-10x "
                              "asymmetry."),
        "gate2_quantity": ("FAIL. The signal is AQ4N plus its reduced "
                           "metabolite AQ4, which 'cannot be separated', and "
                           "hypoxic bioreductive trapping plus stable AQ4-DNA "
                           "intercalation inflate the apparent depth. "
                           "Transport, activation and binding are one number."),
        "also": "No penetration number of any kind is tabulated in the paper.",
    },
    "Grantab 2006 Cancer Res 66:1033 (PMID 16424039), Table 2": {
        "gate1_dose_parity": ("FAIL, and stated outright: '25 umol/L "
                              "paclitaxel, 5 umol/L doxorubicin, 10 umol/L "
                              "methotrexate, and 77 umol/L 5-FU ... permit "
                              "sensitive detection of the drug in compartment "
                              "2'. A 15.4-fold spread across the compared "
                              "compounds."),
        "gate2_quantity": ("PASS in form (radiolabel in a receiving "
                           "compartment, method-consistent across drugs) but "
                           "the authors cap it: radiolabel 'might be "
                           "associated with a metabolite; thus, this provides "
                           "an upper limit'."),
        "note": ("Its dose-matched axis is CELL LINE, not drug. The "
                 "epithelioid-vs-round contrast at fixed drug and dose is "
                 "clean; the drug-vs-drug contrast is not dose-controlled. It "
                 "is the right dataset for a different experiment."),
    },
    "Tunggal 1999 Clin Cancer Res 5:1583 (PMID 10389947)": {
        "gate1_dose_parity": ("FAIL. ~9 uM doxorubicin, 10 uM methotrexate, "
                              "50 uM mitoxantrone, ~78 uM 5-FU — an 8.6-fold "
                              "spread with no rationale given. The "
                              "methotrexate concentration is also "
                              "self-contradictory between Methods (10 uM) and "
                              "the Fig. 4 legend (50 uM)."),
        "gate2_quantity": ("FAIL for the headline pair. Doxorubicin was "
                           "counted as radiolabel (parent plus any labelled "
                           "metabolite) while mitoxantrone was quantified by "
                           "HPLC at 600 nm (parent only) — two different "
                           "measurands presented as one comparison."),
        "also": ("No numeric SD or SEM is printed anywhere in the paper, so it "
                 "cannot supply a measurement-level noise floor."),
    },
}

#: The IC50 handed to every compound. Potency plays NO part in transport (the
#: kill model is strictly downstream of the PDE) and every compound gets the same
#: number, so it cannot differentiate them. It exists because the compile
#: signature requires one.
SHARED_IC50_UM = 25.0

#: The repository's own loosest declared agreement bar, reused rather than
#: re-invented: tmesim/validation/pk_reality.py::ORDER_OF_MAGNITUDE_FACTOR.
ORDER_OF_MAGNITUDE_FACTOR = 3.0

#: A compound's real-section half-distance must beat the SOURCE NULL by at least
#: this factor for the real-section arm to count as measuring transport rather
#: than vessel/KDE geometry. results/vessel_profile.md found the null alone buys
#: d_half = 10.0 um with no drug physics at all.
SOURCE_NULL_MARGIN = 2.0

#: MCC-analogue geometry. dx is FINER than the tissue grid because this is a
#: 1-D uniform slab, not the reconstructed section: 140 um at dx = 12 um is under
#: 12 cells and the 30 um edge windows would be 2.5 cells each. dx is swept.
MCC_DX_UM = 2.0
MCC_DX_SWEEP: Sequence[float] = (1.0, 2.0, 4.0)
MCC_EXPOSURE_S_SWEEP: Sequence[float] = (3600.0, 3600.0 * 6)

#: Explicit-diffusion stability margin. dt = SAFETY * dx^2 / D_eff, with
#: SAFETY well inside the 2-D limit of 0.25. This is a property of THIS script's
#: slab geometry, not an engine parameter, and it is swept through MCC_DX_SWEEP.
MCC_CFL_SAFETY = 0.15

_RESULTS = os.path.join(_ROOT, "results")
_SWEEP = os.path.join(_RESULTS, "sweep")
_PREREG_PATH = os.path.join(_SWEEP, "kyle2004_penetration_prereg.json")
_DATA_PATH = os.path.join(_RESULTS, "kyle2004_penetration_data.json")
_FIG_PATH = os.path.join(_RESULTS, "kyle2004_penetration_map.png")

#: Nuisance axes for the real-section noise floor. Each varies between two runs
#: for reasons that have NOTHING to do with the compound.
NUISANCE_BINS_UM: Sequence[float] = (5.0, 10.0, 20.0)
NUISANCE_ARMS: Sequence[str] = ("dx12_repaired", "dx12_native", "dx30_native")
NUISANCE_FIELDS: Sequence[str] = ("auc", "conc")
NUISANCE_METRICS: Sequence[str] = ("primary", "legacy")


# --------------------------------------------------------------------------- #
# 1. Pre-registration
# --------------------------------------------------------------------------- #
def _measured_ratio_bounds() -> Dict:
    """The measured fold-separation and its own uncertainty, from Table 1 SDs.

    This is the experiment's noise floor, and it is what decides whether a
    difference the model reproduces is bigger than the measurement's own scatter.
    """
    d, u = KYLE_TABLE1[DAU], KYLE_TABLE1[DOX]
    central = d["fluor"] / u["fluor"]
    lo = (d["fluor"] - d["fluor_sd"]) / (u["fluor"] + u["fluor_sd"])
    hi = (d["fluor"] + d["fluor_sd"]) / (u["fluor"] - u["fluor_sd"])
    return {
        "central_fluorescence": central,
        "sd_interval": [lo, hi],
        "central_effect_based": d["effect"] / u["effect"],
        "derivation": ("central = 50/11; interval = (50-15)/(11+5) to "
                       "(50+15)/(11-5), i.e. Table 1 fluorescence +- 1 SD "
                       "combined in the worst-case direction. The effect-based "
                       "column gives an independent 30/12."),
    }


def build_prereg() -> Dict:
    return {
        "_what": ("PRE-REGISTRATION. The bar for 'does the Codon map reproduce "
                  "a MEASURED tissue-penetration difference between two "
                  "compounds'. Frozen before any concentration field was "
                  "solved. Scoring stages assert live code against this file."),
        "_rule": ("Nothing is fitted to the measurement. No default, constant "
                  "or shipped parameter is changed; calibrate_alpha() is "
                  "asserted live at the shipped value. The source was chosen "
                  "for dose parity and species purity (the two gates below) "
                  "BEFORE any model output was inspected."),
        "frozen_utc": "pre-solve",
        "measurement": MEASUREMENT,
        "sources_rejected_and_why": REJECTED_SOURCES,
        "gates_applied_to_choose_the_source": {
            "gate1_dose_parity": ("PASS. All compared drugs at 10 umol/L. An "
                                  "earlier candidate (Trédan 2009 Cancer Res "
                                  "69:940, mitoxantrone vs AQ4N, PMID 19176397) "
                                  "was REJECTED on this gate: 1 vs 10 mg/mL in "
                                  "vitro and 40 vs 200 mg/kg in vivo, chosen "
                                  "'to facilitate detection'."),
            "gate2_quantity_is_transport": (
                "PASS. Parent-drug fluorescence, one shared fluorophore, "
                "emissions cross-calibrated within 10%. The same earlier "
                "candidate was REJECTED on this gate too: its signal was AQ4N "
                "plus its reduced metabolite AQ4, which 'cannot be separated', "
                "with hypoxic bioreductive trapping inflating apparent depth. "
                "Kyle's number DOES include cellular sequestration as well as "
                "diffusion; that is admissible here only because the engine "
                "models the same composite, L = sqrt(D_eff/k), and it is "
                "disclosed wherever the result is quoted."),
        },
        "compounds": {n: {"smiles": SMILES[n], "provenance": PROVENANCE[n],
                          "measured": KYLE_TABLE1[n]} for n in ORDER},
        "in_sample_disclosure": {
            DOX: ("IN-SAMPLE. calibrate_alpha() is fit on doxorubicin, pinning "
                  "its penetration length to 64.92 um by construction. Its "
                  "absolute agreement with Kyle 2004 is NOT evidence and is "
                  "reported NOT TESTABLE, following "
                  "scorecard.py::_run_perivascular_penetration."),
            EPI: ("Out-of-sample as a fit, but the engine cannot distinguish it "
                  "from doxorubicin (4'-epimer; descriptors are "
                  "stereochemistry-blind), so its numbers are expected to be "
                  "doxorubicin's. That is the point of criterion C0."),
            DAU: "OUT-OF-SAMPLE. New to this repository. Nothing fit to it.",
            MIT: ("OUT-OF-SAMPLE for this measurement. It carries a separate "
                  "perivascular anchor in pk_anchors against Patel 2013's band, "
                  "which is reported JOINTLY for three drugs and therefore "
                  "contains no compound-to-compound contrast; nothing was fit "
                  "to it."),
        },
        "engine_state": {
            "calibrate_alpha": SHIPPED_ALPHA,
            "alpha_fit_on": "doxorubicin (Primeau 2005, PMID 16361566)",
            "two_route": "ON (shipping default)",
            "efflux": "OFF (shipping default)",
            "intracellular_retardation": "OFF (shipping default)",
            "section": ("10x Xenium human GBM, hBrain 319-gene panel "
                        "(data_cache/xenium_gbm), the repository's default"),
            "section_dx_um": VP.DX_UM,
            "mcc_dx_um": MCC_DX_UM,
        },
        "arms": {
            "A_mcc_analogue": (
                "The apples-to-apples arm. A uniform slab of the measured "
                "thickness (140 um) with a one-sided planar source and a sealed "
                "far side — the same geometry Kyle's apparatus imposes — solved "
                "by the SAME frozen transport solver. Readout: the ratio of "
                "mean exposure over the near 30 um to mean exposure over the "
                "far 30 um, which is Table 1's quantity with no conversion and "
                "no assumption of exponentiality."),
            "B_real_section_map": (
                "The product's actual output. Both compounds on ONE "
                "reconstructed real tissue section, reduced to exposure vs "
                "distance-from-nearest-measured-endothelial-cell, with the "
                "model-free empirical half-distance and the C0*exp(-d/L)+floor "
                "fit."),
            "C_noise_floor": (
                "Every choice that differs between two runs for reasons "
                "unrelated to the compound: grid, tissue-mask repair, distance "
                "bin width, exposure vs final concentration, distance metric."),
        },
        "criteria": {
            "C0_epimer_control": {
                "statement": ("Doxorubicin and epirubicin are 4'-epimers. The "
                              "engine's descriptors are stereochemistry-blind, "
                              "so it MUST predict them identical — and Kyle "
                              "Table 1 measures them as identical "
                              "(11+-5 vs 12+-3; ~12 vs ~10). A model that "
                              "separated them would predict a difference the "
                              "measurement excludes."),
                "bar": ("relative difference in the engine's near:far ratio "
                        "< 1e-6, AND the measured intervals overlap"),
                "note": ("This is a control on the READOUT, not evidence of "
                         "predictive power: agreeing by being blind is not the "
                         "same as agreeing by being right, and the report must "
                         "say so."),
            },
            "C1_direction": {
                "statement": ("PRIMARY. The engine must place DAUNORUBICIN as "
                              "the WORSE penetrator of the pair, which is what "
                              "the measurement says by 4.5x on the primary "
                              "column and 2.5x on the independent effect-based "
                              "column."),
                "bar": ("near:far ratio(daunorubicin) > near:far "
                        "ratio(doxorubicin) in arm A, AND empirical half-"
                        "distance(daunorubicin) < half-distance(doxorubicin) "
                        "in arm B. Both required."),
                "why_this_can_fail": (
                    "The shipping two-route law consumes logP and HBD. "
                    "Daunorubicin is MORE lipophilic (XLogP 1.8 vs 1.3) and has "
                    "FEWER donors (5 vs 6), so the diffusivity term alone "
                    "predicts the OPPOSITE of the measurement. Recovering the "
                    "measured direction requires the uptake sink to overturn "
                    "the engine's own transport term — and that sink is a "
                    "structure-only prior that alpha_identifiability.md proves "
                    "is unidentified for every drug but doxorubicin."),
            },
            "C2_noise_floor": {
                "statement": ("The separation must be larger than the span "
                              "produced by choices that differ between the two "
                              "runs for reasons unrelated to the compound."),
                "nuisance_axes": {
                    "arm_A": {"dx_um": list(MCC_DX_SWEEP),
                              "exposure_s": list(MCC_EXPOSURE_S_SWEEP)},
                    "arm_B": {"distance_bin_width_um": list(NUISANCE_BINS_UM),
                              "arm": list(NUISANCE_ARMS),
                              "field": list(NUISANCE_FIELDS),
                              "distance_metric": list(NUISANCE_METRICS)},
                },
                "bar": ("the [min, max] interval of the readout over the full "
                        "nuisance grid must be DISJOINT between doxorubicin "
                        "and daunorubicin, in both arms"),
            },
            "C3_magnitude": {
                "statement": ("The modelled fold-separation must agree with the "
                              "measured fold-separation. Two bars are "
                              "registered: a STRICT one set by the "
                              "measurement's own SDs, and a LOOSE one set by "
                              "the repository's declared order-of-magnitude "
                              "factor."),
                "measured": _measured_ratio_bounds(),
                "bar_strict": ("R_sim inside the measurement's own +-1 SD "
                               "interval"),
                "bar_loose": (f"R_sim within {ORDER_OF_MAGNITUDE_FACTOR}x of "
                              f"the central measured fold"),
                "note": ("The loose bar is deliberately reported as weak: a 3x "
                         "band around a 4.5x effect admits R_sim below 1.0, "
                         "i.e. the wrong direction. C1 therefore carries the "
                         "verdict and C3 qualifies it."),
            },
            "C4_absolute_ratios": {
                "statement": ("Secondary and heavily qualified. Each compound's "
                              "own modelled near:far ratio against its Table 1 "
                              "value."),
                "bar": "inside the measured +- 1 SD interval",
                "not_evidence_for": [DOX, EPI],
                "reason": ("doxorubicin is the alpha fit and epirubicin is "
                           "numerically doxorubicin; only daunorubicin's row "
                           "is a prediction"),
            },
            "C6_far_side_disclosure": {
                "statement": ("The paper's own ANOVA finds the far-side "
                              "(post-transport) signal statistically "
                              "INDISTINGUISHABLE across the four drugs, and "
                              "attributes the whole near:far spread to "
                              "near-side accumulation. The engine's far-side "
                              "exposures are therefore REPORTED against that "
                              "finding."),
                "bar": ("NOT SCOREABLE — the paper prints no numeric far-side "
                        "values, only the ANOVA verdict. The engine's far-side "
                        "spread across the four compounds is reported and the "
                        "qualitative agreement stated. No number will be read "
                        "off a figure to manufacture a bar."),
                "why_it_still_matters": (
                    "If the engine reproduces the near:far ratio by predicting "
                    "a large far-side DIFFERENCE, it is getting the right "
                    "answer for a mechanism the measurement excludes."),
            },
            "C7_near_side_accumulation": {
                "statement": ("The paper measures near-side accumulation "
                              "relative to doxorubicin: EPI ~1.5x, DAU ~7x. "
                              "The engine's near-side field is reported "
                              "against those."),
                "bar": ("REGISTERED AS EXPECTED TO FAIL, before the run. The "
                        "engine solves for FREE INTERSTITIAL concentration and "
                        "represents cellular uptake as a LOSS term "
                        "(k = clearance + alpha*uptake_proxy), not as a bound "
                        "reservoir — intracellular retardation was built and "
                        "REFUTED in results/intracellular_retardation.md and "
                        "ships OFF. A stronger sink therefore LOWERS the "
                        "engine's near-side concentration where the "
                        "measurement RAISES total tissue drug. The engine has "
                        "no mechanism to accumulate, so this row is expected "
                        "to come out with the wrong sign."),
                "why_the_ratio_test_survives_it": (
                    "A near:far RATIO is large for a high-uptake compound "
                    "under BOTH descriptions — a store inflates the numerator, "
                    "a sink steepens the gradient. So C1 and C3 remain "
                    "legitimate even though C7 is expected to fail. This is "
                    "registered in advance precisely so it cannot be "
                    "discovered afterwards and presented as a nuance."),
            },
            "C5_real_section_informativeness": {
                "statement": ("results/vessel_profile.md found the drug field "
                              "on the reconstructed section largely inherits "
                              "the shape of the vascular SOURCE field. If that "
                              "confound is total here, arm B measures vessel "
                              "geometry rather than transport and must be "
                              "declared uninformative."),
                "bar": (f"each compound's real-section half-distance must "
                        f"exceed the SOURCE NULL half-distance by at least "
                        f"{SOURCE_NULL_MARGIN}x"),
                "if_it_fails": ("arm B is reported UNINFORMATIVE, arm A carries "
                                "C1, and the limitation is stated as a "
                                "limitation of the map's spatial readout."),
            },
        },
        "registered_predictions_before_solving": {
            "superseded_law": ("The molecular-weight-only law of "
                               "diffusivity_falsification.md predicts a 1.01x "
                               "separation between doxorubicin and "
                               "daunorubicin (543.5 vs 527.5 Da) against a "
                               "measured 4.5x. It is registered as failing this "
                               "test before it is run."),
            "two_route_diffusivity_alone": ("predicts daunorubicin FASTER, "
                                            "i.e. the WRONG direction"),
            "unknown": ("whether the uptake sink overturns it, and whether any "
                        "separation clears the noise floor"),
        },
        "what_this_cannot_establish": [
            "Anything about clinical response. Out of scope.",
            ("Absolute penetration depths. alpha is a gauge fit on doxorubicin "
             "and two_route_promotion.md §8.3 records two external anchors now "
             "disagreeing about where the engine sits."),
            ("That the map is right for compounds unlike these. n = 4, one "
             "structural class, one section, one tumour type."),
            ("That Kyle's ratio is pure diffusion. It is diffusion against "
             "cellular sequestration, and the authors attribute the ranking "
             "mainly to the latter."),
            ("Anything from the doxorubicin or epirubicin rows taken alone."),
        ],
    }


def stage_prereg(verbose: bool = True) -> Dict:
    _assert_alpha()
    prereg = build_prereg()
    os.makedirs(_SWEEP, exist_ok=True)
    if os.path.exists(_PREREG_PATH):
        with open(_PREREG_PATH) as fh:
            frozen = json.load(fh)
        for key in ("measurement", "compounds", "criteria", "engine_state",
                    "arms", "gates_applied_to_choose_the_source",
                    "in_sample_disclosure",
                    "registered_predictions_before_solving"):
            if json.dumps(frozen.get(key), sort_keys=True) != \
                    json.dumps(prereg.get(key), sort_keys=True):
                raise AssertionError(
                    f"PRE-REGISTRATION DRIFT in '{key}'. The frozen file "
                    f"{_PREREG_PATH} disagrees with live code. This is exactly "
                    f"what the registration exists to catch. Resolve it "
                    f"deliberately; do not delete the file.")
        if verbose:
            print(f"PREREG {_PREREG_PATH}\n  exists and MATCHES live code.  OK")
    else:
        with open(_PREREG_PATH, "w") as fh:
            json.dump(prereg, fh, indent=1)
            fh.write("\n")
        if verbose:
            print(f"PREREG written to {_PREREG_PATH}")
    if verbose:
        mb = _measured_ratio_bounds()
        print(f"  calibrate_alpha() = {calibrate_alpha()!r}  "
              f"(asserted == {SHIPPED_ALPHA!r})")
        print(f"  measurement: PMID {MEASUREMENT['pmid']} Table 1, "
              f"all drugs at {MEASUREMENT['concentration_uM']} uM")
        print(f"  measured fold DAU/DOX = "
              f"{mb['central_fluorescence']:.2f}x  "
              f"(+-1SD interval {mb['sd_interval'][0]:.2f}-"
              f"{mb['sd_interval'][1]:.2f}x; effect-based "
              f"{mb['central_effect_based']:.2f}x)")
        print(f"  REGISTERED RISK: the two-route diffusivity term alone "
              f"predicts the WRONG direction.")
    return prereg


def _require_prereg() -> Dict:
    if not os.path.exists(_PREREG_PATH):
        raise SystemExit(
            "REFUSING TO RUN: the pre-registration does not exist. Run "
            "`--stage prereg` first. The bar is fixed before the measurement.")
    return stage_prereg(verbose=False)


def _assert_alpha() -> float:
    a = calibrate_alpha()
    if a != SHIPPED_ALPHA:
        raise AssertionError(
            f"calibrate_alpha() = {a!r}, expected the shipped {SHIPPED_ALPHA!r}. "
            f"A shipped constant has moved; this experiment is invalid until "
            f"that is explained.")
    return a


# --------------------------------------------------------------------------- #
# 2. Compiling — ONE path, so only the molecule differs
# --------------------------------------------------------------------------- #
def compile_all(alpha: float, two_route: bool = True) -> Dict[str, Dict]:
    out: Dict[str, Dict] = {}
    for name in ORDER:
        mol = _descriptors.from_smiles(name, SMILES[name])
        prof = build_profile(mol, ic50=SHARED_IC50_UM, two_route=two_route)
        dp = molecule_to_drugparams(prof, ic50=SHARED_IC50_UM)
        dp = calibrated_drugparams(dp, alpha)
        dp = apply_barrier_to_drugparams(dp, mol)
        tr = prof.two_route or {}
        out[name] = {
            "name": name, "source": mol.source,
            "mw": float(mol.get("mw")), "logp": float(mol.get("logp")),
            "hbd": float(mol.get("hbd")), "hba": float(mol.get("hba")),
            "tpsa": float(mol.get("tpsa")),
            "formal_charge": float(mol.get("formal_charge")),
            "radius_nm": float(mol.get("radius_hydrodynamic_nm")),
            "interstitial_D": float(prof.interstitial_diffusivity),
            "D_eff": float(dp.diffusivity),
            "two_route_hindrance": float(tr.get("hindrance", float("nan"))),
            "membrane_P_um_s": float(
                tr.get("membrane_permeability_um_s", float("nan"))),
            "binding": float(dp.binding), "clearance": float(dp.clearance),
            "k_loss": float(dp.binding + dp.clearance),
            "bbb": float(dp.bbb_penetrance),
            "drug_params": dp, "molecule": mol,
        }
    return out


def stage_compounds(verbose: bool = True) -> Dict:
    _require_prereg()
    alpha = _assert_alpha()
    ship = compile_all(alpha, two_route=True)
    old = compile_all(alpha, two_route=False)

    def L_an(v, ecm=0.4):
        return math.sqrt(v["D_eff"] * max(1.0 - ALPHA_TORTUOSITY * ecm, 1e-6)
                         / v["k_loss"])

    res = {
        "alpha": alpha,
        "shipping": {k: {kk: vv for kk, vv in v.items()
                         if kk not in ("drug_params", "molecule")}
                     for k, v in ship.items()},
        "superseded_mw_only_D": {k: v["interstitial_D"] for k, v in old.items()},
        "L_analytic_um": {k: L_an(v) for k, v in ship.items()},
    }
    res["superseded_fold_DAU_over_DOX"] = (
        res["superseded_mw_only_D"][DAU] / res["superseded_mw_only_D"][DOX])
    res["epimer_D_relative_difference"] = abs(
        ship[EPI]["D_eff"] - ship[DOX]["D_eff"]) / max(ship[DOX]["D_eff"], 1e-30)
    res["epimer_k_relative_difference"] = abs(
        ship[EPI]["k_loss"] - ship[DOX]["k_loss"]) / max(ship[DOX]["k_loss"],
                                                         1e-30)

    if verbose:
        print("\n" + "=" * 82)
        print("STAGE: compounds")
        print("=" * 82)
        print(f"\n{'':<15}{'MW':>9}{'logP':>8}{'HBD':>5}{'TPSA':>7}"
              f"{'D_inter':>10}{'hind':>8}{'P_mem':>11}{'D_eff':>10}"
              f"{'k_loss':>11}{'L':>8}")
        for n in ORDER:
            v = ship[n]
            print(f"{n:<15}{v['mw']:>9.2f}{v['logp']:>8.3f}{v['hbd']:>5.0f}"
                  f"{v['tpsa']:>7.1f}{v['interstitial_D']:>10.3f}"
                  f"{v['two_route_hindrance']:>8.4f}"
                  f"{v['membrane_P_um_s']:>11.3g}{v['D_eff']:>10.3f}"
                  f"{v['k_loss']:>11.3e}{res['L_analytic_um'][n]:>8.2f}")
        print(f"\nSUPERSEDED MW-only law, DAU/DOX fold: "
              f"{res['superseded_fold_DAU_over_DOX']:.4f}x "
              f"(against a measured 4.55x)")
        print(f"EPIMER DEGENERACY  dox vs epi: D rel-diff "
              f"{res['epimer_D_relative_difference']:.3e}, k rel-diff "
              f"{res['epimer_k_relative_difference']:.3e}")
        print(f"ANALYTIC L fold DAU/DOX: "
              f"{res['L_analytic_um'][DAU] / res['L_analytic_um'][DOX]:.4f}x "
              f"(>1 means the engine says daunorubicin goes FURTHER, which is "
              f"the WRONG direction)")
    return res


# --------------------------------------------------------------------------- #
# 3. Arm A — the MCC analogue
# --------------------------------------------------------------------------- #
def _slab(ecm: float, thickness_um: float, dx_um: float, ny: int = 5):
    """A uniform slab: one-sided planar source at column 0, sealed far side.

    This is Kyle's apparatus as the engine's frozen solver sees it. ``ny`` is
    small and the field is uniform along y, so the solve is effectively 1-D; the
    engine's zero-flux boundaries are exactly the sealed far side.
    """
    nx = int(round(float(thickness_um) / float(dx_um))) + 1
    dom = Domain(nx=nx, ny=ny, dx=float(dx_um))
    vasc = np.zeros((ny, nx), dtype=np.float64)
    vasc[:, 0] = 1.0
    fields = TissueFields(ecm_density=np.full((ny, nx), float(ecm)),
                          vascular_density=vasc)
    cells = CellPopulation(x=np.array([0.0]), y=np.array([0.0]),
                           cell_type=np.array([0]), state=np.array([0]),
                           sensitivity=np.array([1.0]))
    return Microenvironment(domain=dom, fields=fields, cells=cells,
                            patient_id="mcc_slab"), nx


def mcc_ratio(dp, ecm: float, dx_um: float = MCC_DX_UM,
              thickness_um: Optional[float] = None,
              exposure_s: float = 3600.0,
              window_um: Optional[float] = None) -> Dict:
    """The engine's near:far exposure ratio across a slab — Table 1's quantity."""
    thickness_um = (MEASUREMENT["thickness_um"] if thickness_um is None
                    else thickness_um)
    window_um = (MEASUREMENT["edge_window_um"] if window_um is None
                 else window_um)
    micro, nx = _slab(ecm, thickness_um, dx_um)
    d_eff = float(dp.diffusivity) * max(1.0 - ALPHA_TORTUOSITY * ecm, 1e-6)
    dt = min(2.0, MCC_CFL_SAFETY * dx_um * dx_um / max(d_eff, 1e-9))
    pde = PDEConfig(convection=False, dt=dt,
                    total_time=float(exposure_s), save_every=100000)
    res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
    auc = np.asarray(res.fields["auc"], dtype=np.float64).mean(axis=0)
    x = np.arange(nx) * float(dx_um)
    near = auc[x <= window_um]
    far = auc[x >= (thickness_um - window_um)]
    ratio = float(near.mean() / far.mean()) if far.mean() > 0 else float("inf")
    return {"ratio": ratio, "x_um": x.tolist(), "auc": auc.tolist(),
            "near_over_far": ratio,
            "near_mean": float(near.mean()), "far_mean": float(far.mean()),
            "dx_um": float(dx_um), "exposure_s": float(exposure_s),
            "thickness_um": float(thickness_um)}


def run_mcc(ecm: float, verbose: bool = True) -> Dict:
    alpha = _assert_alpha()
    compiled = compile_all(alpha, two_route=True)
    out: Dict = {"ecm_density": ecm, "drugs": {}, "sweep": []}
    for n in ORDER:
        dp = compiled[n]["drug_params"]
        out["drugs"][n] = mcc_ratio(dp, ecm)
        for dxs in MCC_DX_SWEEP:
            for te in MCC_EXPOSURE_S_SWEEP:
                r = mcc_ratio(dp, ecm, dx_um=dxs, exposure_s=te)
                out["sweep"].append({"drug": n, "dx_um": dxs, "exposure_s": te,
                                     "ratio": r["ratio"]})
        if verbose:
            meas = KYLE_TABLE1[n]
            ms = (f"{meas['fluor']:.0f} +- {meas['fluor_sd']:.0f}"
                  if meas["fluor"] else f"{meas['effect_kind']}"
                                        f"{meas['effect']:.0f} (effect only)")
            print(f"  {n:<15} engine near:far = {out['drugs'][n]['ratio']:8.2f}"
                  f"     measured {ms}")
    return out


# --------------------------------------------------------------------------- #
# 4. Arm B — the map on real tissue
# --------------------------------------------------------------------------- #
def _profile_all(field: np.ndarray, dist: np.ndarray, tissue: np.ndarray,
                 bin_um: float) -> Dict:
    prof = VP.distance_profile(field, dist, tissue, bin_um=bin_um)
    d, y = np.asarray(prof["distance_um"]), np.asarray(prof["mean"])
    cnt = np.asarray(prof["count"])
    keep = cnt >= VP.MIN_BIN_COUNT
    d, y, w = d[keep], y[keep], cnt[keep].astype(np.float64)
    out = {"distance_um": d.tolist(), "mean": y.tolist(), "n_bins": int(d.size),
           "max_um": float(d.max()) if d.size else float("nan"),
           "empirical_d_half_um": VP.empirical_half_distance_um(d, y)}
    fit = VP.fit_exponential_plus_floor(d, y, w) if d.size >= 4 else {}
    out["fit_L_um"] = float(fit.get("L_um", float("nan")))
    out["fit_r2_linear"] = float(fit.get("r2_linear_space", float("nan")))
    out["fit_floor_fraction"] = float(fit.get("floor_fraction", float("nan")))
    out["fit_error"] = fit.get("error")
    return out


def run_map(arm_key: str = "dx12_repaired", bin_um: float = VP.BIN_UM,
            verbose: bool = True) -> Dict:
    alpha = _assert_alpha()
    arm = next(a for a in VP.ARMS if a["key"] == arm_key)
    t0 = time.time()
    micro, _cn, _cr = VP.build_arm_microenvironment(dx_um=arm["dx_um"],
                                                    repair=arm["repair"])
    tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
    endo = VP.endothelial_xy(micro.cells)
    dist_p = VP.direct_vessel_distance_um(micro.domain, endo)
    legacy = VP.vessel_mask_from_density(micro.fields.vascular_density, tissue)
    dist_l = VP.distance_from_vessel_um(legacy["mask"], float(arm["dx_um"]))
    ecm = float(np.asarray(micro.fields.ecm_density)[tissue].mean())
    if verbose:
        print(f"  tissue rebuilt in {time.time() - t0:.1f}s: "
              f"{micro.domain.ny}x{micro.domain.nx} @ dx={arm['dx_um']} um, "
              f"{int(tissue.sum())} tissue px, {endo.shape[0]} endothelial "
              f"cells, mean ECM {ecm:.4f}")

    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True,
                    save_every=500)
    compiled = compile_all(alpha, two_route=True)
    out: Dict = {"arm": arm_key, "dx_um": float(arm["dx_um"]),
                 "bin_um": float(bin_um), "ecm_density": ecm, "drugs": {},
                 "n_endothelial": int(endo.shape[0]),
                 "n_tissue_px": int(tissue.sum())}
    for n in ORDER:
        t1 = time.time()
        dp = compiled[n]["drug_params"]
        res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
        auc = np.asarray(res.fields["auc"], dtype=np.float64)
        conc = np.asarray(res.transport.final_concentration, dtype=np.float64)
        rec = {"auc_primary": _profile_all(auc, dist_p, tissue, bin_um),
               "conc_primary": _profile_all(conc, dist_p, tissue, bin_um),
               "auc_legacy": _profile_all(auc, dist_l, tissue, bin_um),
               "mean_auc": float(res.extras["mean_auc"]),
               "solve_s": time.time() - t1, "_auc_field": auc}
        out["drugs"][n] = rec
        if verbose:
            print(f"  {n:<15} solved {rec['solve_s']:5.1f}s   d_half = "
                  f"{_fmt(rec['auc_primary']['empirical_d_half_um'])}   "
                  f"L_fit = {rec['auc_primary']['fit_L_um']:.1f} um")
    out["source_null"] = _profile_all(
        np.asarray(micro.fields.vascular_density, dtype=np.float64),
        dist_p, tissue, bin_um)
    out["_micro"], out["_tissue"] = micro, tissue
    if verbose:
        print(f"  SOURCE NULL     d_half = "
              f"{_fmt(out['source_null']['empirical_d_half_um'])}  "
              f"(no drug physics at all)")
    return out


def run_noise(verbose: bool = True) -> Dict:
    alpha = _assert_alpha()
    compiled = compile_all(alpha, two_route=True)
    pde = PDEConfig(dt=10.0, total_time=3600.0 * 6, convection=True,
                    save_every=500)
    rows: List[Dict] = []
    for arm_key in NUISANCE_ARMS:
        arm = next(a for a in VP.ARMS if a["key"] == arm_key)
        micro, _cn, _cr = VP.build_arm_microenvironment(dx_um=arm["dx_um"],
                                                        repair=arm["repair"])
        tissue = np.asarray(micro.domain.tissue_mask, dtype=bool)
        legacy = VP.vessel_mask_from_density(micro.fields.vascular_density,
                                             tissue)
        dists = {"primary": VP.direct_vessel_distance_um(
                     micro.domain, VP.endothelial_xy(micro.cells)),
                 "legacy": VP.distance_from_vessel_um(legacy["mask"],
                                                      float(arm["dx_um"]))}
        for n in (DOX, DAU):
            dp = compiled[n]["drug_params"]
            res = simulate_response(micro, SimConfig(drug=dp, pde=pde, seed=0))
            flds = {"auc": np.asarray(res.fields["auc"], dtype=np.float64),
                    "conc": np.asarray(res.transport.final_concentration,
                                       dtype=np.float64)}
            for fk, fv in flds.items():
                for mk, dv in dists.items():
                    for b in NUISANCE_BINS_UM:
                        p = _profile_all(fv, dv, tissue, b)
                        rows.append({"drug": n, "arm": arm_key, "field": fk,
                                     "metric": mk, "bin_um": b,
                                     "d_half_um": p["empirical_d_half_um"],
                                     "fit_L_um": p["fit_L_um"]})
        if verbose:
            print(f"  {arm_key:<15} done")
    return {"rows": rows}


def _span(rows: List[Dict], drug: str, key: str = "d_half_um"):
    vals = [r[key] for r in rows if r["drug"] == drug]
    fin = [float(v) for v in vals if v is not None and np.isfinite(v)]
    never = sum(1 for v in vals if v is None or not np.isfinite(v))
    if not fin:
        return None, None, 0, never
    return min(fin), max(fin), len(fin), never


# --------------------------------------------------------------------------- #
# 5. Scoring
# --------------------------------------------------------------------------- #
def score(comp: Dict, mcc: Dict, mapping: Dict, noise: Optional[Dict]) -> Dict:
    v: Dict[str, Dict] = {}
    mb = _measured_ratio_bounds()
    r = {n: mcc["drugs"][n]["ratio"] for n in ORDER}

    # C0 --------------------------------------------------------------------
    rel = abs(r[EPI] - r[DOX]) / max(abs(r[DOX]), 1e-30)
    v["C0_epimer_control"] = {
        "pass": bool(rel < 1e-6),
        "engine_ratio": {DOX: r[DOX], EPI: r[EPI]},
        "relative_difference": rel,
        "measured": {DOX: "11 +- 5", EPI: "12 +- 3"},
        "detail": ("the engine must NOT separate a measured-identical epimer "
                   "pair (bar: relative difference < 1e-6)"),
    }

    # C1 --------------------------------------------------------------------
    dir_A = r[DAU] > r[DOX]
    h = {n: mapping["drugs"][n]["auc_primary"]["empirical_d_half_um"]
         for n in ORDER}

    def _shallower(x, y):
        xi = x is None or not np.isfinite(x)
        yi = y is None or not np.isfinite(y)
        if xi and yi:
            return None
        if yi:
            return True          # y never halves => x is shallower
        if xi:
            return False
        return float(x) < float(y)

    dir_B = _shallower(h[DAU], h[DOX])
    v["C1_direction"] = {
        "pass": bool(dir_A and dir_B is True),
        "arm_A_mcc": {"ratio": {DOX: r[DOX], DAU: r[DAU]},
                      "daunorubicin_worse": bool(dir_A)},
        "arm_B_map": {"d_half_um": h, "daunorubicin_shallower": dir_B},
        "measured_direction": ("daunorubicin penetrates WORSE: 50 +- 15 vs "
                               "11 +- 5 (fluorescence), ~30 vs ~12 (effect)"),
        "detail": "both arms must place daunorubicin as the worse penetrator",
    }

    # C2 --------------------------------------------------------------------
    A = {n: [s["ratio"] for s in mcc["sweep"] if s["drug"] == n]
         for n in (DOX, DAU)}
    a_disjoint = bool(min(A[DAU]) > max(A[DOX]) or min(A[DOX]) > max(A[DAU]))
    if noise is None:
        b_disjoint = None
        spans = None
    else:
        sd, sa = _span(noise["rows"], DOX), _span(noise["rows"], DAU)
        spans = {DOX: sd, DAU: sa}
        b_disjoint = (bool(sa[0] > sd[1] or sd[0] > sa[1])
                      if (sd[0] is not None and sa[0] is not None) else None)
    v["C2_noise_floor"] = {
        "pass": (bool(a_disjoint and b_disjoint) if b_disjoint is not None
                 else None),
        "arm_A_span": {n: [min(A[n]), max(A[n])] for n in (DOX, DAU)},
        "arm_A_disjoint": a_disjoint,
        "arm_B_span": ({n: {"min": s[0], "max": s[1], "n_finite": s[2],
                            "n_never_halves": s[3]}
                        for n, s in spans.items()} if spans else None),
        "arm_B_disjoint": b_disjoint,
        "detail": "the two compounds' readout intervals must be disjoint",
    }

    # C3 --------------------------------------------------------------------
    R_sim = r[DAU] / r[DOX] if r[DOX] > 0 else float("nan")
    lo, hi = mb["sd_interval"]
    strict = bool(np.isfinite(R_sim) and lo <= R_sim <= hi)
    loose = bool(np.isfinite(R_sim) and R_sim > 0 and
                 abs(math.log(R_sim / mb["central_fluorescence"]))
                 <= math.log(ORDER_OF_MAGNITUDE_FACTOR))
    v["C3_magnitude"] = {
        "pass": strict, "pass_loose": loose,
        "R_sim": R_sim, "R_measured": mb,
        "detail": (f"strict bar: R_sim in [{lo:.2f}, {hi:.2f}] "
                   f"(measurement's own +-1 SD); loose bar: within "
                   f"{ORDER_OF_MAGNITUDE_FACTOR}x of "
                   f"{mb['central_fluorescence']:.2f}"),
    }

    # C4 --------------------------------------------------------------------
    rows = {}
    for n in ORDER:
        m = KYLE_TABLE1[n]
        if m["fluor"] is None:
            rows[n] = {"engine": r[n], "measured": f">{m['effect']:.0f} "
                                                   f"(effect-based only)",
                       "in_band": None}
            continue
        blo, bhi = m["fluor"] - m["fluor_sd"], m["fluor"] + m["fluor_sd"]
        rows[n] = {"engine": r[n],
                   "measured": f"{m['fluor']:.0f} +- {m['fluor_sd']:.0f}",
                   "band": [blo, bhi], "in_band": bool(blo <= r[n] <= bhi)}
    v["C4_absolute_ratios"] = {
        "pass": rows[DAU]["in_band"], "rows": rows,
        "detail": ("only the daunorubicin row is evidence; doxorubicin is the "
                   "alpha fit and epirubicin is numerically doxorubicin"),
    }

    # C6 / C7 ---------------------------------------------------------------
    near = {n: mcc["drugs"][n]["near_mean"] for n in ORDER}
    far = {n: mcc["drugs"][n]["far_mean"] for n in ORDER}
    far_rel = {n: (far[n] / far[DOX] if far[DOX] > 0 else float("nan"))
               for n in ORDER}
    near_rel = {n: (near[n] / near[DOX] if near[DOX] > 0 else float("nan"))
                for n in ORDER}
    fvals = [x for x in far_rel.values() if np.isfinite(x) and x > 0]
    v["C6_far_side_disclosure"] = {
        "pass": None,
        "engine_far_side_relative_to_dox": far_rel,
        "engine_far_side_spread_fold": (max(fvals) / min(fvals) if fvals
                                        else float("nan")),
        "measured": ("statistically indistinguishable across drugs, P > 0.05 "
                     "(ANOVA), except EPI > DOX at 10 uM x h"),
        "detail": ("NOT SCOREABLE — the paper prints no numeric far-side "
                   "values. Reported so a right answer cannot be reached by a "
                   "mechanism the measurement excludes."),
    }
    meas_near = MEASUREMENT["near_side_accumulation_vs_dox"]
    signs = {n: (None if meas_near[n] is None or not np.isfinite(near_rel[n])
                 else bool((near_rel[n] - 1.0) * (meas_near[n] - 1.0) >= 0))
             for n in ORDER}
    v["C7_near_side_accumulation"] = {
        "pass": None,
        "engine_near_side_relative_to_dox": near_rel,
        "measured_relative_to_dox": meas_near,
        "sign_agrees": signs,
        "detail": ("REGISTERED AS EXPECTED TO FAIL before the run: the engine "
                   "solves for FREE interstitial drug and models uptake as a "
                   "LOSS, not a reservoir, so it cannot accumulate near the "
                   "source. Reported, not scored."),
    }

    # C5 --------------------------------------------------------------------
    sn = mapping["source_null"]["empirical_d_half_um"]
    ratios = {}
    for n in ORDER:
        hh = h[n]
        if hh is None or not np.isfinite(hh):
            ratios[n] = float("inf")
        elif sn and np.isfinite(sn) and sn > 0:
            ratios[n] = float(hh) / float(sn)
        else:
            ratios[n] = float("nan")
    v["C5_real_section_informativeness"] = {
        "pass": bool(all((not np.isfinite(x)) or x >= SOURCE_NULL_MARGIN
                         for x in ratios.values())),
        "source_null_d_half_um": sn, "ratio_vs_null": ratios,
        "detail": (f"each compound's half-distance must exceed the source null "
                   f"({_fmt(sn)}) by >= {SOURCE_NULL_MARGIN}x"),
    }

    scored = [x["pass"] for x in v.values() if x["pass"] is not None]
    v["_overall"] = {"all_scored_pass": bool(all(scored)),
                     "n_scored": len(scored),
                     "n_not_scoreable": sum(1 for x in v.values()
                                            if x.get("pass") is None)}
    return v


def _fmt(x, nd: int = 1) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "never halves"
    return f"{float(x):.{nd}f} um"


# --------------------------------------------------------------------------- #
# 6. Figure
# --------------------------------------------------------------------------- #
def make_figure(comp: Dict, mcc: Dict, mapping: Dict,
                noise: Optional[Dict]) -> str:
    micro, tissue = mapping["_micro"], mapping["_tissue"]
    dx = mapping["dx_um"]
    ext = [0, micro.domain.nx * dx, 0, micro.domain.ny * dx]
    colour = {DOX: "#B2182B", EPI: "#D6604D", DAU: "#2166AC", MIT: "#4D4D4D"}

    fig = plt.figure(figsize=(18.0, 9.8))
    gs = GridSpec(2, 4, figure=fig, height_ratios=[1.0, 0.9], hspace=0.34,
                  wspace=0.42)

    norm = {}
    for n in ORDER:
        f = np.asarray(mapping["drugs"][n]["_auc_field"], dtype=np.float64)
        p99 = np.percentile(f[tissue], 99) if tissue.any() else 1.0
        norm[n] = np.clip(np.where(tissue, f / max(p99, 1e-30), np.nan), 0, 1)

    for i, n in enumerate((DOX, DAU)):
        ax = fig.add_subplot(gs[0, i])
        im = ax.imshow(norm[n], origin="lower", cmap="magma", vmin=0, vmax=1,
                       extent=ext)
        ax.set_title(f"{n}\nAUC / p99(AUC) — same tissue, same solver",
                     fontsize=10)
        ax.set_xlabel("µm")
        if i == 0:
            ax.set_ylabel("µm")
        plt.colorbar(im, ax=ax, fraction=0.040, pad=0.03,
                     label="fraction of p99")

    ax = fig.add_subplot(gs[0, 2])
    diff = norm[DAU] - norm[DOX]
    vm = float(np.nanmax(np.abs(diff))) if np.isfinite(diff).any() else 1.0
    vm = max(vm, 1e-6)
    im = ax.imshow(diff, origin="lower", cmap="RdBu_r", vmin=-vm, vmax=vm,
                   extent=ext)
    ax.set_title("daunorubicin − doxorubicin\n(each on its own p99 scale)",
                 fontsize=10)
    ax.set_xlabel("µm")
    plt.colorbar(im, ax=ax, fraction=0.040, pad=0.03,
                     label="Δ fraction of p99")

    # measured vs modelled near:far
    ax = fig.add_subplot(gs[0, 3])
    xs = np.arange(len(ORDER))
    for i, n in enumerate(ORDER):
        m = KYLE_TABLE1[n]
        if m["fluor"] is not None:
            ax.errorbar(i - 0.13, m["fluor"], yerr=m["fluor_sd"], fmt="s",
                        color="0.25", ms=7, capsize=4,
                        label="measured (Kyle 2004 Table 1)" if i == 0 else None)
        else:
            ax.plot(i - 0.13, 30, "s", color="0.25", ms=7, mfc="none")
            ax.annotate(">30", (i - 0.13, 33), color="0.25", fontsize=8,
                        ha="center", va="bottom")
        ax.plot(i + 0.13, mcc["drugs"][n]["ratio"], "o", color=colour[n], ms=9,
                label="Codon (MCC-analogue arm)" if i == 0 else None)
    ax.set_xticks(xs)
    ax.set_xticklabels([n[:4] for n in ORDER], fontsize=8)
    ax.set_yscale("log")
    ax.set_ylabel("near : far exposure ratio")
    ax.set_title("The test — measured vs modelled\n(140 µm slab, 30 µm edge "
                 "windows)", fontsize=9.5)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.25, which="both")

    ax = fig.add_subplot(gs[1, 0])
    for n in ORDER:
        d = mcc["drugs"][n]
        y = np.asarray(d["auc"], dtype=float)
        y = y / y[0] if y.size and y[0] > 0 else y
        ax.semilogy(d["x_um"], np.clip(y, 1e-8, None), color=colour[n], lw=1.8,
                    label=n)
    ax.set_xlabel("depth into the slab (µm)")
    ax.set_ylabel("exposure / near-surface exposure")
    ax.set_title("Arm A — the MCC analogue", fontsize=10)
    ax.legend(fontsize=7.5)

    ax = fig.add_subplot(gs[1, 1])
    for n in ORDER:
        p = mapping["drugs"][n]["auc_primary"]
        y = np.asarray(p["mean"], dtype=float)
        y = y / y[0] if y.size and y[0] > 0 else y
        ax.plot(p["distance_um"], y, "o-", color=colour[n], ms=3.5, lw=1.6,
                label=n)
    sn = mapping["source_null"]
    ys = np.asarray(sn["mean"], dtype=float)
    if ys.size and ys[0] > 0:
        ax.plot(sn["distance_um"], ys / ys[0], "k--", lw=1.2,
                label="SOURCE NULL (no drug physics)")
    ax.axhline(0.5, color="0.5", lw=0.8, ls=":")
    ax.set_xlabel("distance from nearest endothelial cell (µm)")
    ax.set_ylabel("AUC / perivascular AUC")
    ax.set_title("Arm B — the map's own readout", fontsize=10)
    ax.legend(fontsize=7)
    ax.set_ylim(0, 1.05)

    ax = fig.add_subplot(gs[1, 2])
    for j, n in enumerate((DOX, DAU)):
        vals = [s["ratio"] for s in mcc["sweep"] if s["drug"] == n]
        rng = np.random.default_rng(0).uniform(-0.1, 0.1, len(vals))
        ax.scatter(np.full(len(vals), j) + rng, vals, s=22, color=colour[n],
                   alpha=0.8)
        ax.hlines([min(vals), max(vals)], j - 0.22, j + 0.22, color=colour[n],
                  lw=2)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([DOX, DAU], fontsize=8)
    ax.set_ylabel("near : far ratio")
    ax.set_title("Noise floor, arm A\n(grid, exposure window)", fontsize=9.5)

    ax = fig.add_subplot(gs[1, 3])
    if noise is not None:
        for j, n in enumerate((DOX, DAU)):
            vals = [r["d_half_um"] for r in noise["rows"] if r["drug"] == n]
            fin = [float(x) for x in vals if x is not None and np.isfinite(x)]
            nev = len(vals) - len(fin)
            if fin:
                rng = np.random.default_rng(1).uniform(-0.1, 0.1, len(fin))
                ax.scatter(np.full(len(fin), j) + rng, fin, s=20,
                           color=colour[n], alpha=0.8)
                ax.hlines([min(fin), max(fin)], j - 0.22, j + 0.22,
                          color=colour[n], lw=2)
            ax.annotate(f"{nev}/{len(vals)} never halve", (j, 0.02),
                        xycoords=("data", "axes fraction"), ha="center",
                        fontsize=7, color=colour[n])
        ax.set_xticks([0, 1])
        ax.set_xticklabels([DOX, DAU], fontsize=8)
        ax.set_ylabel("empirical d½ (µm)")
        ax.set_title("Noise floor, arm B\n(grid, mask, bin, field, metric)",
                     fontsize=9.5)
    else:
        ax.axis("off")

    fig.suptitle(
        "Does the map reproduce a MEASURED penetration difference?   "
        "doxorubicin / epirubicin / daunorubicin / mitoxantrone vs Kyle 2004 "
        "Table 1 (PMID 15342419), all drugs at 10 µM",
        fontsize=11.5, y=0.985)
    fig.savefig(_FIG_PATH, dpi=150)
    plt.close(fig)
    return _FIG_PATH


# --------------------------------------------------------------------------- #
# 7. Driver
# --------------------------------------------------------------------------- #
def _jsonable(o):
    if isinstance(o, dict):
        return {k: _jsonable(x) for k, x in o.items() if not k.startswith("_")}
    if isinstance(o, (list, tuple)):
        return [_jsonable(x) for x in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def main() -> Dict:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=["prereg", "compounds", "mcc", "all"])
    ap.add_argument("--arm", default="dx12_repaired")
    ap.add_argument("--skip-noise", action="store_true")
    args = ap.parse_args()

    t0 = time.time()
    if args.stage == "prereg":
        stage_prereg()
        return {}
    comp = stage_compounds()
    if args.stage == "compounds":
        return {}

    print("\n" + "=" * 82)
    print("ARM B: the map — every compound on ONE real tissue section")
    print("=" * 82)
    mapping = run_map(arm_key=args.arm)

    print("\n" + "=" * 82)
    print("ARM A: the MCC analogue — Kyle 2004's geometry, Kyle 2004's quantity")
    print("=" * 82)
    mcc = run_mcc(mapping["ecm_density"])
    if args.stage == "mcc":
        return {}

    noise = None
    if not args.skip_noise:
        print("\n" + "=" * 82)
        print("ARM C: the noise floor")
        print("=" * 82)
        noise = run_noise()

    verdicts = score(comp, mcc, mapping, noise)
    print("\n" + "=" * 82)
    print("VERDICTS, against the frozen pre-registration")
    print("=" * 82)
    for k, val in verdicts.items():
        if k.startswith("_"):
            continue
        p = val.get("pass")
        tag = "PASS" if p is True else ("FAIL" if p is False
                                        else "NOT SCOREABLE")
        print(f"  {k:<36} {tag}")
        print(f"      {val.get('detail', '')}")
    print(f"\n  overall: {verdicts['_overall']}")

    out = {"compounds": comp, "mcc": mcc, "map": mapping, "noise": noise,
           "verdicts": verdicts, "measurement": MEASUREMENT,
           "alpha": calibrate_alpha(), "elapsed_s": time.time() - t0}
    fig = make_figure(comp, mcc, mapping, noise)
    with open(_DATA_PATH, "w") as fh:
        json.dump(_jsonable(out), fh, indent=1)
    print(f"\nwrote {_DATA_PATH}\nwrote {fig}\ntotal {time.time() - t0:.1f}s")
    return out


if __name__ == "__main__":
    main()
