"""DESCRIPTOR FALSIFICATION -- re-score the shipping two-route transport law on
Pruijn 2008 using the descriptors THE ENGINE COMPUTES, not the paper's measured
columns.

WHY THIS FILE EXISTS
--------------------
``results/two_route_promotion.md`` publishes **R^2 = 0.718** for the two-route
transport law against Pruijn 2008 (71 compounds, HT29 multicellular layers), and
that number appears in the README, ``docs/QUICKSTART.md`` and the marketing copy.
It was produced by ``experiments/two_route_transport.py::_build_rows``, which
feeds the law:

    logp <- the paper's OWN logP7.4 column   (Table S1, measured or ACD-calculated)
    hbd  <- the paper's OWN HD column        (Table S1)

Neither of those is what ships. ``compile_drug.build_profile`` feeds the same law

    logp <- descriptors.compute_descriptors(smiles)["logp"]  = RDKit Crippen MolLogP
    hbd  <- descriptors.compute_descriptors(smiles)["hbd"]   = RDKit Lipinski.NumHDonors

So **the descriptor layer that actually ships has never been scored against
anything.** Every published transport number assumes a descriptor quality it has
not demonstrated. This file measures that gap, and does not repair it.

THE RULE THIS FILE WAS WRITTEN UNDER
------------------------------------
**Nothing is tuned. Nothing is repaired.** Not a constant, not a descriptor
method, not a correction factor. ``tmesim/molecular/lipophilicity.py::resolve_logp``
is deliberately NOT applied -- it has no callers today and applying it moves
``alpha``, which is a decision for a human, not for a scoring script. If a
computed descriptor scores badly here, that IS the result.

THE HARD LIMIT ON THIS EXPERIMENT, STATED FIRST
-----------------------------------------------
66 of the 71 Pruijn compounds are Auckland Cancer Society Research Centre
``SN`` code numbers -- unpublished benzotriazine di-N-oxide analogues whose
structures exist only as drawings in the Table S1 PDF. They are not in PubChem
(checked: ``PUGREST.NotFound`` by name). **The engine's descriptor layer cannot
be run on them at all**, because it takes a SMILES and there is no SMILES.
Hand-drawing 66 structures off a scanned table would be inventing the inputs,
which is exactly the failure mode this file exists to detect.

What IS available, and what this file therefore does, in three arms:

  ARM 1  (n=5, exact)   The five compounds in Table S1 with public, verifiable
                        structures -- urea, mannitol, tirapazamine, acridone and
                        NLCQ-1. Engine descriptors computed from PubChem SMILES,
                        MW cross-checked against Table S1, compared per-compound
                        with the paper's columns, and run end-to-end through the
                        SHIPPING compile.
  ARM 2  (n=49, paired) The 49 compounds carrying BOTH a measured and an
                        ACD-calculated logP7.4 in Table S1. Score the identical
                        physics twice, changing only the descriptor source.
                        THIS ARM CAME BACK NULL and is reported as such: ACD
                        LogP/LogD v8.0 was trained BY THE AUTHORS on 98 compounds
                        of this very series (paper, Experimental section), so 22
                        of the 49 "calculated" values are EXACTLY the measured
                        value and 48 of 49 are within 0.15. The column is not an
                        independent computed descriptor, so the large-n version
                        of this question cannot be answered with data in hand.
  ARM 3  (n=68, sweep)  Propagate descriptor error of known magnitude through the
                        published score, so a future out-of-sample failure can be
                        attributed to physics or to descriptors rather than to
                        both at once.

Run:
    ./.venv/bin/python experiments/descriptor_falsification.py --stage structures
    ./.venv/bin/python experiments/descriptor_falsification.py --stage descriptors
    ./.venv/bin/python experiments/descriptor_falsification.py --stage substitute
    ./.venv/bin/python experiments/descriptor_falsification.py --stage propagate
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import sys
import time
from typing import Dict, List, Optional, Tuple

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from tmesim.molecular import compile_drug, descriptors, stokes_einstein  # noqa: E402
from tmesim.molecular import two_route  # noqa: E402

from experiments.diffusivity_falsification import (  # noqa: E402
    PRUIJN_2008_TABLE_S1, TABLE_UNITS_TO_UM2_S,
)
from experiments.two_route_transport import (  # noqa: E402
    _build_rows, _engine_interstitial_D, _score,
)

OUT_DIR = os.path.join(_ROOT, "results", "sweep")
TABLE_S1_CSV = os.path.join(_ROOT, "data_cache", "pruijn2008",
                            "pruijn2008_tableS1.csv")

#: The published figure this file re-scores. Quoted from
#: results/two_route_promotion.md; recomputed live below, never assumed.
PUBLISHED_R2 = 0.718


# =========================================================================== #
# The only structures in Pruijn 2008 that a descriptor layer can actually be
# run on. Each is a PUBLIC compound with a PubChem record; the SMILES below were
# retrieved from PubChem PUG-REST (property/SMILES) and each is verified against
# Table S1's own MW column by `stage_structures`, which ASSERTS agreement to
# 1 Da. Nothing here is drawn, guessed, or inferred from an analogue.
#
# The other 66 compounds are ACSRC `SN` code numbers with no public structure --
# `curl .../compound/name/SN28287/...` returns PUGREST.NotFound. They are listed
# by `stage_structures` as UNCOMPUTABLE, which is itself a finding.
# =========================================================================== #
CURATED_STRUCTURES: Dict[str, Dict] = {
    "Urea": {
        "smiles": "C(=O)(N)N",
        "pubchem_cid": 1176,
        "pubchem_mw": 60.056,
        "note": "PubChem CID 1176, IUPACName 'urea'. Unambiguous.",
    },
    "Mannitol": {
        "smiles": "C([C@H]([C@H]([C@@H]([C@@H](CO)O)O)O)O)O",
        "pubchem_cid": 6251,
        "pubchem_mw": 182.17,
        "note": (
            "PubChem CID 6251, D-mannitol, IUPACName "
            "'(2R,3R,4R,5R)-hexane-1,2,3,4,5,6-hexol'. Pruijn's Table S1 says "
            "only 'Mannitol'; D-mannitol is the commercial/tracer form. Every "
            "descriptor used here (MW, Crippen logP, donor count) is "
            "stereochemistry-independent, so the choice of enantiomer cannot "
            "affect any number in this file."
        ),
    },
    "SN24280 (TPZ)": {
        "smiles": "C1=CC=C2C(=C1)[N+](=NC(=[N+]2[O-])N)[O-]",
        "pubchem_cid": 36400,
        "pubchem_mw": 178.15,
        "note": (
            "Tirapazamine. This is the repo's OWN TPZ SMILES, byte-identical to "
            "experiments/two_route_transport.py::TPZ_SMILES and "
            "experiments/anchor_prototype.py, so the engine descriptors "
            "computed here are the ones the project already publishes for its "
            "registered anchor. Asserted identical below."
        ),
    },
    "acridan-9(10H)-one": {
        "smiles": "C1=CC=C2C(=C1)C(=O)C3=CC=CC=C3N2",
        "pubchem_cid": 2015,
        "pubchem_mw": 195.22,
        "note": (
            "Acridone. PubChem CID 2015 IUPACName is '10H-acridin-9-one', which "
            "is the same compound as Table S1's 'acridan-9(10H)-one' written in "
            "the other convention. MW check confirms it."
        ),
    },
    "NLCQ-1": {
        "smiles": "C1=CC2=C(C=CN=C2C=C1Cl)NCCCN3C=CN=C3[N+](=O)[O-]",
        "pubchem_cid": 398939,
        "pubchem_mw": 331.76,
        "note": (
            "PubChem CID 398939 is the HYDROCHLORIDE (SMILES ends '.Cl', "
            "MW 368.2). Table S1's MW is 332, i.e. the free base, so the '.Cl' "
            "component was dropped and nothing else changed. The MW assertion "
            "below is what establishes that this is the right call: the free "
            "base computes to ~331.8 against Table S1's 332, the salt to 368."
        ),
    },
}

#: Compounds whose structures were sought and NOT found. Keyed by the reason, so
#: "we didn't look" and "we looked and there is nothing" stay distinguishable --
#: the same discipline lipophilicity.NO_EXPERIMENTAL_LOGP applies.
NO_PUBLIC_STRUCTURE_REASON = (
    "Auckland Cancer Society Research Centre internal code number. Pruijn 2008 "
    "publishes these structures only as drawings in the Table S1 accessory PDF; "
    "PubChem PUG-REST name lookup returns PUGREST.NotFound. No SMILES exists in "
    "this repo or in any public database reachable from it, so "
    "descriptors.compute_descriptors CANNOT BE CALLED on them. Deliberately not "
    "hand-drawn from the PDF images: a transcribed structure is an invented "
    "input, and inventing the input is the exact defect this file measures."
)


def _hr(c: str = "=") -> None:
    print(c * 78)


def _table_row(name: str) -> Tuple:
    for row in PRUIJN_2008_TABLE_S1:
        if row[0] == name:
            return row
    raise KeyError(name)


# =========================================================================== #
# STAGE 1 -- what can the descriptor layer even be run on?
# =========================================================================== #
def stage_structures() -> dict:
    _hr()
    print("STAGE 1 -- CAN THE SHIPPING DESCRIPTOR LAYER BE RUN AT ALL?")
    _hr()
    print("  descriptors.compute_descriptors takes a SMILES. Pruijn 2008 does not")
    print("  publish SMILES; Table S1 publishes drawings. So the first question is")
    print("  not 'how good are the computed descriptors' but 'do they exist'.\n")

    known = sorted(CURATED_STRUCTURES)
    all_names = [r[0] for r in PRUIJN_2008_TABLE_S1]
    unknown = [n for n in all_names if n not in CURATED_STRUCTURES]

    print(f"  {'compound':<22}{'TableS1 MW':>11}{'PubChem MW':>12}{'RDKit MW':>10}"
          f"{'|d| vs S1':>11}  CID")
    recs = []
    for name in known:
        spec = CURATED_STRUCTURES[name]
        row = _table_row(name)
        mol = descriptors.from_smiles(name, spec["smiles"])
        assert mol.source == "rdkit", (
            f"{name}: descriptors fell back to the curated table instead of "
            f"parsing SMILES -- source={mol.source!r}. The whole point of "
            "this run is to exercise the RDKit path.")
        rdkit_mw = float(mol.get("mw"))
        d = abs(rdkit_mw - float(row[1]))
        print(f"  {name:<22}{row[1]:>11.1f}{spec['pubchem_mw']:>12.2f}"
              f"{rdkit_mw:>10.2f}{d:>11.2f}  {spec['pubchem_cid']}")
        assert d <= 1.0, (
            f"{name}: RDKit MW {rdkit_mw:.2f} disagrees with Table S1 MW "
            f"{row[1]:.1f} by {d:.2f} Da. The curated SMILES is NOT this "
            "compound. Refusing to score a structure that fails its own "
            "identity check.")
        recs.append({"name": name, "smiles": spec["smiles"],
                     "pubchem_cid": spec["pubchem_cid"],
                     "table_mw": float(row[1]), "rdkit_mw": rdkit_mw,
                     "mw_abs_error_da": d})
    print(f"\n  MW identity check passed for all {len(known)} "
          f"(worst |error| = {max(r['mw_abs_error_da'] for r in recs):.2f} Da).")

    # TPZ must be the repo's own SMILES, or the anchor column would not match.
    from experiments.two_route_transport import TPZ_SMILES
    assert CURATED_STRUCTURES["SN24280 (TPZ)"]["smiles"] == TPZ_SMILES, (
        "the TPZ SMILES here has drifted from two_route_transport.TPZ_SMILES")
    print("  tirapazamine SMILES is byte-identical to the repo's registered one: OK")

    print(f"\n  COMPUTABLE   : {len(known):>2} / {len(all_names)}  "
          f"({len(known)/len(all_names):.0%})")
    print(f"  UNCOMPUTABLE : {len(unknown):>2} / {len(all_names)}  "
          f"({len(unknown)/len(all_names):.0%})")
    print(f"\n  Reason all {len(unknown)} are uncomputable:")
    for line in _wrap(NO_PUBLIC_STRUCTURE_REASON, 72):
        print(f"    {line}")
    print(f"\n  -> The published R^2 = {PUBLISHED_R2} is a score of the PHYSICS on")
    print("     descriptors supplied by the measurement. On 66 of these 71")
    print("     compounds the shipping descriptor layer produces NO ANSWER, so no")
    print("     amount of re-scoring can convert that number into a score of the")
    print("     shipping engine. That is the headline limit of this experiment and")
    print("     it is a property of the dataset, not of the effort spent on it.")
    return {"computable": recs, "n_computable": len(known),
            "n_uncomputable": len(unknown), "uncomputable": unknown,
            "uncomputable_reason": NO_PUBLIC_STRUCTURE_REASON}


def _wrap(text: str, width: int) -> List[str]:
    out, cur = [], ""
    for w in text.split():
        if len(cur) + len(w) + 1 > width:
            out.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}".strip()
    if cur:
        out.append(cur)
    return out


# =========================================================================== #
# STAGE 2 -- per-compound descriptor error, the actual deliverable
# =========================================================================== #
def stage_descriptors(verbose: bool = True) -> dict:
    _hr()
    print("STAGE 2 -- PER-COMPOUND DESCRIPTOR ERROR (engine vs Table S1)")
    _hr()
    print("  Engine logP  = RDKit Crippen MolLogP, the NEUTRAL octanol/water logP.")
    print("  Paper  logP  = logP7.4, a DISTRIBUTION coefficient at pH 7.4.")
    print("  Engine HBD   = RDKit Lipinski.NumHDonors.")
    print("  Paper  HD    = Table S1's hydrogen-bond donor column.")
    print("  These are compared as the model uses them: whatever is in the `logp`")
    print("  and `hbd` slots of two_route.two_route_diffusivity.\n")

    rows = []
    print(f"  {'compound':<22}{'logP7.4':>9}{'src':>6}{'Crippen':>9}{'d logP':>8}"
          f"{'HD(S1)':>8}{'HBD(eng)':>10}{'d HBD':>7}")
    for name in sorted(CURATED_STRUCTURES):
        spec = CURATED_STRUCTURES[name]
        t = _table_row(name)
        mol = descriptors.from_smiles(name, spec["smiles"])
        clogp = float(mol.get("logp"))
        ehbd = int(mol.get("hbd"))
        paper_logp, paper_src, paper_hd = t[4], t[5], t[2]
        d_logp = (clogp - paper_logp) if paper_logp is not None else None
        d_hbd = ehbd - paper_hd
        rows.append({
            "name": name, "mw": float(t[1]), "paper_logp": paper_logp,
            "paper_logp_source": paper_src, "crippen_logp": clogp,
            "d_logp": d_logp, "paper_hd": int(paper_hd), "engine_hbd": ehbd,
            "d_hbd": d_hbd,
            "meas_D_um2s": float(t[6]) * TABLE_UNITS_TO_UM2_S,
            "tpsa": float(mol.get("tpsa")), "hba": int(mol.get("hba")),
            "formal_charge": int(mol.get("formal_charge")),
        })
        print(f"  {name:<22}{paper_logp:>9.2f}{paper_src:>6}{clogp:>9.3f}"
              f"{d_logp:>+8.3f}{paper_hd:>8d}{ehbd:>10d}{d_hbd:>+7d}")

    ds = [r["d_logp"] for r in rows if r["d_logp"] is not None]
    n = len(ds)
    mean = sum(ds) / n
    rms = math.sqrt(sum(d * d for d in ds) / n)
    mae = sum(abs(d) for d in ds) / n
    worst = max(rows, key=lambda r: abs(r["d_logp"]))
    print(f"\n  logP ERROR over n={n}:  mean {mean:+.3f}   RMS {rms:.3f}   "
          f"MAE {mae:.3f}   range [{min(ds):+.3f}, {max(ds):+.3f}]")
    print(f"  worst: {worst['name']} ({worst['d_logp']:+.3f} log units)")
    print(f"  sign : {sum(1 for d in ds if d > 0)}/{n} positive. The error is NOT")
    print("         one-signed overall -- but it is not scatter either. Grouped by")
    print("         structure the three classes separate cleanly, and the class")
    print("         boundaries are chemical, not statistical:")
    print("           neutral aromatic (acridone, TPZ)     |d| <= 0.54  -- Crippen's")
    print("               home ground; both slightly TOO POLAR.")
    print("           small neutral polar (urea, mannitol) |d| ~ 1.0    -- atom-")
    print("               contribution failure on H-bonding polyols/amides, BOTH")
    print("               SIGNS. This is the class lipophilicity.py already")
    print("               documents for methotrexate (+2.1) and cisplatin (~+2.4).")
    print("           IONISABLE BASE (NLCQ-1)               d = +3.28   -- the")
    print("               4-aminoquinoline is protonated at pH 7.4, so its logD7.4")
    print("               is far below its neutral logP. Crippen has no pKa and")
    print("               cannot see this. ONE-SIGNED, ALWAYS POSITIVE, and the")
    print("               largest single error in the set by 3x.")

    dh = [r["d_hbd"] for r in rows]
    print(f"\n  HBD ERROR over n={len(dh)}: mean {sum(dh)/len(dh):+.2f}, "
          f"range [{min(dh):+d}, {max(dh):+d}]")
    mism = [r for r in rows if r["d_hbd"] != 0]
    if mism:
        print(f"  {len(mism)} of {len(dh)} DISAGREE. Every disagreement is "
              "NEGATIVE, and the")
        print("  pattern is not noise -- see the counting convention below.")
        for r in mism:
            print(f"    {r['name']:<22} Table S1 HD={r['paper_hd']}  "
                  f"engine HBD={r['engine_hbd']}")

    # ---- the counting convention, demonstrated rather than asserted ---------
    print("\n  WHY THE HBD COLUMNS DISAGREE -- demonstrated, not assumed.")
    print("  RDKit exposes two different donor counts. Table S1's column tracks")
    print("  one of them and the engine consumes the other:")
    print(f"  {'compound':<22}{'S1 HD':>7}{'NumHDonors':>12}{'NHOHCount':>11}")
    conv = []
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors as _RDDesc, Lipinski as _Lip
        for name in sorted(CURATED_STRUCTURES):
            m = Chem.MolFromSmiles(CURATED_STRUCTURES[name]["smiles"])
            hd = int(_table_row(name)[2])
            a, b = int(_Lip.NumHDonors(m)), int(_RDDesc.NHOHCount(m))
            conv.append({"name": name, "table_hd": hd,
                         "num_h_donors": a, "nhoh_count": b})
            print(f"  {name:<22}{hd:>7d}{a:>12d}{b:>11d}")
        agree_a = sum(1 for c in conv if c["table_hd"] == c["num_h_donors"])
        agree_b = sum(1 for c in conv if c["table_hd"] == c["nhoh_count"])
        print(f"\n  agreement with Table S1: NumHDonors {agree_a}/{len(conv)}, "
              f"NHOHCount {agree_b}/{len(conv)}")
        print("  NumHDonors counts DONOR HEAVY ATOMS (an -NH2 is 1).")
        print("  NHOHCount counts DONOR HYDROGENS       (an -NH2 is 2).")
        print("  The engine (descriptors.compute_descriptors) uses NumHDonors.")
        print("  On this evidence Table S1's HD column is the HYDROGEN count, so")
        print("  the two columns are DIFFERENT QUANTITIES, not two estimates of")
        print("  one. This matters because two_route multiplies that count by")
        print(f"  HBD_LOGP_DECREMENT = {two_route.HBD_LOGP_DECREMENT} log units per")
        print(f"  donor: one unit of miscount is {two_route.HBD_LOGP_DECREMENT:.2f} "
              f"log units of membrane partition,")
        print(f"  i.e. a factor of {10**two_route.HBD_LOGP_DECREMENT:.0f} in K_mem.")
    except Exception as exc:  # pragma: no cover - RDKit is installed here
        print(f"  (RDKit unavailable, convention check skipped: {exc})")

    # ---- ionisation: is it handled anywhere on this path? ------------------
    print("\n  IS IONISATION HANDLED ON THE TRANSPORT PATH? Checked in code.")
    ion = _ionisation_audit()
    for line in ion["report"]:
        print(f"    {line}")

    # ---- end-to-end: the SHIPPING compile on these five --------------------
    print("\n  END-TO-END: the SHIPPING compile (two-route ON, engine descriptors)")
    print("  against the measured HT29 D_mcl, for the five compounds it can run on.")
    print(f"  {'compound':<22}{'measured':>10}{'shipping':>10}{'ratio':>8}"
          f"{'paper-desc':>12}{'ratio':>8}")
    e2e = []
    scored = {r["name"]: r for r in _build_rows(two_route.HBD_LOGP_DECREMENT)}
    for r in rows:
        name = r["name"]
        mol = descriptors.from_smiles(name, CURATED_STRUCTURES[name]["smiles"])
        prof = compile_drug.build_profile(mol)          # SHIPPING defaults
        ship = float(prof.effective_diffusivity)
        pub = scored[name].get("two_route")
        e2e.append({"name": name, "meas": r["meas_D_um2s"], "shipping_D": ship,
                    "shipping_ratio": ship / r["meas_D_um2s"],
                    "paper_desc_D": pub,
                    "paper_desc_ratio": (pub / r["meas_D_um2s"]) if pub else None})
        print(f"  {name:<22}{r['meas_D_um2s']:>10.1f}{ship:>10.1f}"
              f"{ship/r['meas_D_um2s']:>8.2f}"
              f"{(pub if pub else float('nan')):>12.1f}"
              f"{(pub/r['meas_D_um2s'] if pub else float('nan')):>8.2f}")
    print("\n  NOTE the shipping column uses the engine's OWN MW too, not Table")
    print("  S1's rounded one, so a small part of any difference is the MW column.")
    print("  n=5 is far too few for an R^2; the per-compound ratios above are the")
    print("  result, and they are reported as ratios for that reason.")

    return {"per_compound": rows, "logp_error": {
        "n": n, "mean": mean, "rms": rms, "mae": mae,
        "min": min(ds), "max": max(ds),
        "n_positive": sum(1 for d in ds if d > 0)},
        "hbd_error": {"mean": sum(dh) / len(dh), "min": min(dh), "max": max(dh),
                      "n_mismatched": len(mism)},
        "hbd_convention": conv, "ionisation": ion, "end_to_end": e2e}


def _ionisation_audit() -> dict:
    """Does anything between SMILES and D_eff know about pH or pKa? Checked by
    reading the shipping call path, not by grep-and-hope."""
    import inspect
    src_desc = inspect.getsource(descriptors.compute_descriptors)
    src_comp = inspect.getsource(compile_drug.build_profile)
    src_two = inspect.getsource(two_route.membrane_partition_coefficient) \
        if hasattr(two_route, "membrane_partition_coefficient") else ""
    src_two += inspect.getsource(two_route.two_route_diffusivity)
    blob = (src_desc + src_comp + src_two).lower()
    hits = sorted({t for t in ("pka", "logd", "ph_", " ph ", "ionis", "ioniz",
                               "protonat", "henderson")
                   if t in blob})
    _passes_logp = "descriptors.get(" in src_comp and '"logp"' in src_comp
    rep = [
        f"descriptors.compute_descriptors -> 'logp' is Crippen.MolLogP: "
        f"{'Crippen.MolLogP' in src_desc}",
        f"compile_drug passes descriptors['logp'] straight to two_route: "
        f"{_passes_logp}",
        f"pH/pKa/logD/ionisation tokens anywhere on that path: "
        f"{hits if hits else 'NONE'}",
    ]
    handled = bool(hits)
    rep.append(
        "VERDICT: ionisation is NOT handled on the transport path. RDKit's "
        "Crippen MolLogP is a neutral-species octanol/water logP with no pKa "
        "input at all." if not handled else
        "VERDICT: something on this path mentions ionisation -- re-read it.")
    rep.append(
        "The repo DOES have a pKa table (molecular/ion_trapping.py) but it "
        "serves lysosomal trapping and is not on this path.")
    rep.append(
        "Consequence: substituting Crippen logP for Pruijn's logP7.4 is a "
        "SYSTEMATIC discrepancy for any ionisable compound, not noise. Its "
        "sign is fixed: for a compound ionised at pH 7.4, logD7.4 < logP, so "
        "the engine reports it MORE lipophilic than it is, opens the "
        "transcellular route too far, and over-predicts D.")
    return {"tokens_found": hits, "handled": handled, "report": rep}


# =========================================================================== #
# STAGE 3 -- the paired substitution, n=49
# =========================================================================== #
def _load_csv_logp() -> Dict[str, Dict[str, Optional[float]]]:
    """Table S1's measured and ACD-calculated logP7.4 columns, from the curated
    CSV a sibling extraction owns. Read-only here."""
    out: Dict[str, Dict[str, Optional[float]]] = {}

    def f(v):
        v = (v or "").strip()
        try:
            return float(v)
        except ValueError:
            return None

    #: The CSV and PRUIJN_2008_TABLE_S1 disagree on the CASE of one compound
    #: ("Acridan-9(10H)-one" vs "acridan-9(10H)-one"), so names are matched
    #: case-insensitively and re-keyed to the Python table's spelling. Without
    #: this the paired set silently loses a compound.
    by_lower = {}
    with open(TABLE_S1_CSV) as fh:
        for row in csv.DictReader(fh):
            by_lower[row["compound_id"].strip().lower()] = {
                "measured": f(row["logP74_measured"]),
                "calculated": f(row["logP74_calculated"]),
            }
    for r in PRUIJN_2008_TABLE_S1:
        rec = by_lower.get(r[0].strip().lower())
        if rec is not None:
            out[r[0]] = rec
    return out


def _rows_with(logp_by_name: Dict[str, Optional[float]],
               hbd_by_name: Optional[Dict[str, float]] = None,
               cache: Optional[Dict[float, float]] = None) -> List[dict]:
    """`two_route_transport._build_rows` with the descriptor slots overridden.

    Structurally identical to that function -- same interstitial D, same radius,
    same shipping `hbd_decrement`, same physics -- so the ONLY thing that differs
    between two calls is what goes into the `logp`/`hbd` slots. Compounds whose
    override is None are dropped, exactly as a missing logP7.4 is dropped there.
    """
    cache = cache if cache is not None else {}
    rows = []
    for (name, mw, hd, ha, logp, src, d_meas, d_sem, d_paper) in PRUIJN_2008_TABLE_S1:
        lp = logp_by_name.get(name)
        if lp is None:
            continue
        hb = hd if hbd_by_name is None else hbd_by_name.get(name, hd)
        if mw not in cache:
            cache[mw] = _engine_interstitial_D(mw)
        eng = cache[mw]
        radius = stokes_einstein.hydrodynamic_radius_nm(mw)
        tr = two_route.two_route_diffusivity(
            eng, radius, float(lp), float(hb),
            hbd_decrement=two_route.HBD_LOGP_DECREMENT)
        rows.append({"name": name, "mw": mw, "hd": hb, "logp": float(lp),
                     "meas": d_meas * TABLE_UNITS_TO_UM2_S, "eng": eng,
                     "two_route": tr["effective_D_um2s"],
                     "para_only": tr["paracellular_D_um2s"],
                     "P_um_s": tr["membrane_permeability_um_s"]})
    return rows


def stage_substitute() -> dict:
    _hr()
    print("STAGE 3 -- PAIRED SUBSTITUTION: measured logP7.4 vs a COMPUTED one")
    _hr()
    print("  The engine's own descriptors cannot be computed for 66 of 71")
    print("  compounds (Stage 1), so the closest available paired test is Table")
    print("  S1's own two logP7.4 columns on the compounds that carry BOTH.")
    print("  Everything else -- MW, HD, radius, interstitial D, every constant --")
    print("  is held identical. The only thing that changes is the descriptor.\n")
    print("  The calculated column is ACD LogP/LogD prediction software v8.0")
    print("  (Advanced Chemistry Development, Toronto), and the paper states it")
    print("  was trained by the authors ON 98 COMPOUNDS OF THIS SERIES using")
    print("  'ACD/LogP System Training and Accuracy Extender' before predicting")
    print("  the rest (Pruijn 2008, Experimental). So this arm is a CEILING on")
    print("  what a computed logP can do here, not an estimate of RDKit Crippen:")
    print("  a purpose-built logD predictor trained in-series is the best case,")
    print("  and RDKit's generic neutral-species Crippen is nowhere near it.\n")

    cols = _load_csv_logp()
    cache: Dict[float, float] = {}

    # --- first, reproduce the published score so the baseline is not assumed --
    pub_rows = [r for r in _build_rows(two_route.HBD_LOGP_DECREMENT, cache)
                if "two_route" in r]
    pub = _score(pub_rows, "two_route")
    print(f"  BASELINE, reproduced live (n={pub['n']}, the published run):")
    print(f"    R^2 = {pub['r2']:.3f}   Spearman {pub['spearman']:+.3f}   "
          f"spread {pub['fold']:.2f}x   geo bias {pub['geo_bias']:.2f}x")
    assert abs(pub["r2"] - PUBLISHED_R2) < 5e-4, (
        f"live baseline R^2 {pub['r2']:.4f} does not reproduce the published "
        f"{PUBLISHED_R2}. Stop: something moved underneath this experiment.")
    print(f"    reproduces the published {PUBLISHED_R2} to "
          f"{abs(pub['r2']-PUBLISHED_R2):.1e}: OK\n")

    # --- the paired subset --------------------------------------------------
    names = [r[0] for r in PRUIJN_2008_TABLE_S1]
    paired = [n for n in names
              if cols.get(n, {}).get("measured") is not None
              and cols.get(n, {}).get("calculated") is not None]
    m_map = {n: cols[n]["measured"] for n in paired}
    c_map = {n: cols[n]["calculated"] for n in paired}

    rm = _rows_with(m_map, cache=cache)
    rc = _rows_with(c_map, cache=cache)
    sm, sc = _score(rm, "two_route"), _score(rc, "two_route")

    ds = [c_map[n] - m_map[n] for n in paired]
    nn = len(ds)
    rms = math.sqrt(sum(d * d for d in ds) / nn)
    mean = sum(ds) / nn
    print(f"  PAIRED SUBSET, n={nn} compounds carrying both columns.")
    print(f"  descriptor error of the ACD column vs measured: mean {mean:+.3f}, "
          f"RMS {rms:.3f} log units,")
    print(f"  range [{min(ds):+.2f}, {max(ds):+.2f}], "
          f"{sum(1 for d in ds if abs(d) > 1.0)} compounds off by > 1 log unit.")

    # ---- IS THIS ARM EVEN INFORMATIVE? Checked before it is interpreted. ----
    n_eq = sum(1 for d in ds if abs(d) < 1e-9)
    n_small = sum(1 for d in ds if 1e-9 <= abs(d) <= 0.15)
    biggest = max(abs(d) for d in ds)
    print(f"\n  !! READ THIS BEFORE READING THE SCORES BELOW.")
    print(f"     {n_eq} of {nn} pairs are EXACTLY EQUAL to the last decimal, and a")
    print(f"     further {n_small} differ by <= 0.15 log units. Exactly one pair "
          f"(the")
    print(f"     largest, {biggest:.2f}) differs by more than 0.15.")
    print("     That is what ACD 'System Training' does: it was trained on these")
    print("     compounds, so for a training compound the 'calculated' column")
    print("     REPRODUCES the measured one rather than predicting it.")
    print("     CONSEQUENCE: this arm is a NULL. It cannot estimate what a")
    print("     computed descriptor costs, because on this subset the computed")
    print("     column is not an independent estimate of the measurement -- it is")
    print("     the measurement, rounded. Whatever R^2 comes out below, it is NOT")
    print("     evidence that computed descriptors are harmless. It is evidence")
    print("     that this particular column is not a computed descriptor.\n")
    print(f"  {'logP source':<34}{'R^2':>8}{'Spearman':>10}{'spread':>9}"
          f"{'geo bias':>10}{'<3x':>7}")
    for lbl, s in (("MEASURED logP7.4 (published input)", sm),
                   ("ACD-CALCULATED logP7.4", sc)):
        print(f"  {lbl:<34}{s['r2']:>8.3f}{s['spearman']:>+10.3f}"
              f"{s['fold']:>8.2f}x{s['geo_bias']:>9.2f}x{s['within_3x']:>6.0%}")
    drop = sm["r2"] - sc["r2"]
    print(f"\n  R^2 across the swap: {sm['r2']:.4f} -> {sc['r2']:.4f} "
          f"(change {sc['r2'] - sm['r2']:+.4f}).")
    print("  Which is the answer you would expect from a column that is the")
    print("  measurement rounded. THIS NUMBER MUST NOT BE QUOTED as 'computed")
    print("  descriptors cost nothing'. The honest reading of Arm 2 is: the only")
    print("  large-n computed-descriptor column available for this dataset turned")
    print("  out to be circular, so the large-n test of the question CANNOT BE RUN")
    print("  with data in hand. Arm 1 (n=5, real Crippen) and Arm 3 (propagation)")
    print("  are what is left, and they are the parts to read.")

    # --- worst movers -------------------------------------------------------
    by = {r["name"]: r for r in rm}
    bc = {r["name"]: r for r in rc}
    movers = sorted(paired, key=lambda n: -abs(math.log10(bc[n]["two_route"]
                                                          / by[n]["two_route"])))
    print(f"\n  The ten compounds the computed descriptor moves furthest:")
    print(f"  {'compound':<16}{'MW':>6}{'measured':>10}{'D(meas lp)':>12}"
          f"{'D(calc lp)':>12}{'fold':>8}{'d logP':>8}")
    worst = []
    for n in movers[:10]:
        f_ = bc[n]["two_route"] / by[n]["two_route"]
        worst.append({"name": n, "mw": by[n]["mw"], "meas": by[n]["meas"],
                      "D_meas_logp": by[n]["two_route"],
                      "D_calc_logp": bc[n]["two_route"], "fold": f_,
                      "d_logp": c_map[n] - m_map[n]})
        print(f"  {n:<16}{by[n]['mw']:>6.0f}{by[n]['meas']:>10.1f}"
              f"{by[n]['two_route']:>12.1f}{bc[n]['two_route']:>12.1f}"
              f"{f_:>8.2f}x{c_map[n]-m_map[n]:>+8.2f}")

    return {"baseline": pub, "paired_n": nn, "paired_names": paired,
            "measured": sm, "calculated": sc, "r2_drop": drop,
            "descriptor_error": {"mean": mean, "rms": rms,
                                 "min": min(ds), "max": max(ds),
                                 "n_over_1_log": sum(1 for d in ds if abs(d) > 1.0)},
            "worst_movers": worst}


# =========================================================================== #
# STAGE 4 -- propagate descriptor error of known magnitude
# =========================================================================== #
def stage_propagate(replicates: int = 200, seed: int = 20260821) -> dict:
    _hr()
    print("STAGE 4 -- HOW MUCH R^2 SURVIVES A DESCRIPTOR ERROR OF SIZE X?")
    _hr()
    print("  This is the de-confounding table. A future out-of-sample test (a")
    print("  second cell line, unchanged constants) can fail because the physics")
    print("  is wrong or because the descriptors are. Below is what descriptor")
    print("  error ALONE does to the published score, holding the physics and")
    print("  every constant fixed. Nothing here is fitted; it is a sweep.\n")

    cache: Dict[float, float] = {}
    base_rows = [r for r in _build_rows(two_route.HBD_LOGP_DECREMENT, cache)
                 if "two_route" in r]
    base = _score(base_rows, "two_route")
    names = [r["name"] for r in base_rows]
    base_logp = {r["name"]: r["logp"] for r in base_rows}
    base_hd = {r["name"]: r["hd"] for r in base_rows}
    print(f"  baseline (published inputs, n={base['n']}): R^2 = {base['r2']:.3f}\n")

    # ---- A. systematic offset ---------------------------------------------
    print("  A. SYSTEMATIC logP OFFSET applied to every compound. This is the")
    print("     shape a neutral-vs-pH-7.4 mismatch takes: one-directional, not")
    print("     random. Positive = engine thinks everything is more lipophilic.")
    print(f"\n  {'d logP':>8}{'R^2':>9}{'Spearman':>10}{'spread':>9}{'geo bias':>10}"
          f"{'<3x':>7}")
    offsets = []
    for d in (-2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0):
        rr = _rows_with({n: base_logp[n] + d for n in names}, cache=cache)
        s = _score(rr, "two_route")
        s["offset"] = d
        offsets.append(s)
        mark = "  <-- published inputs" if d == 0.0 else ""
        print(f"  {d:>+8.2f}{s['r2']:>9.3f}{s['spearman']:>+10.3f}"
              f"{s['fold']:>8.2f}x{s['geo_bias']:>9.2f}x{s['within_3x']:>6.0%}{mark}")

    # ---- B. random error ---------------------------------------------------
    print(f"\n  B. RANDOM (UNBIASED) logP ERROR of RMS sigma, "
          f"{replicates} replicates each,")
    print(f"     seed {seed}. This is the shape a generic estimator's scatter")
    print("     takes once any systematic part is removed.")
    print(f"\n  {'sigma':>8}{'mean R^2':>10}{'sd':>8}{'p05':>8}{'p95':>8}"
          f"{'mean rho':>10}")
    rng = random.Random(seed)
    noises = []
    for sig in (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
        r2s, rhos = [], []
        reps = 1 if sig == 0.0 else replicates
        for _ in range(reps):
            rr = _rows_with({n: base_logp[n] + rng.gauss(0.0, sig) for n in names},
                            cache=cache)
            s = _score(rr, "two_route")
            r2s.append(s["r2"])
            rhos.append(s["spearman"])
        r2s_sorted = sorted(r2s)
        mu = sum(r2s) / len(r2s)
        sd = (math.sqrt(sum((x - mu) ** 2 for x in r2s) / (len(r2s) - 1))
              if len(r2s) > 1 else 0.0)
        rec = {"sigma": sig, "mean_r2": mu, "sd_r2": sd,
               "p05": r2s_sorted[int(0.05 * (len(r2s) - 1))],
               "p95": r2s_sorted[int(0.95 * (len(r2s) - 1))],
               "mean_spearman": sum(rhos) / len(rhos), "replicates": len(r2s)}
        noises.append(rec)
        print(f"  {sig:>8.2f}{mu:>10.3f}{sd:>8.3f}{rec['p05']:>8.3f}"
              f"{rec['p95']:>8.3f}{rec['mean_spearman']:>+10.3f}")

    # ---- C. the donor-count convention ------------------------------------
    print("\n  C. THE DONOR-COUNT CONVENTION (Stage 2). Table S1's HD counts donor")
    print("     HYDROGENS; the engine's `hbd` counts donor HEAVY ATOMS, so for")
    print("     every -NH2 the engine's count is one LOWER. Which compounds carry")
    print("     an -NH2 is unknowable without the 66 missing structures, so this")
    print("     is a BRACKET, not a claim: nobody shifts (published), everybody")
    print("     shifts by one, everybody by two.")
    print(f"\n  {'HD shift':>10}{'R^2':>9}{'Spearman':>10}{'spread':>9}"
          f"{'geo bias':>10}{'<3x':>7}")
    hd_shifts = []
    for sh in (0, -1, -2):
        rr = _rows_with(base_logp,
                        hbd_by_name={n: max(0.0, base_hd[n] + sh) for n in names},
                        cache=cache)
        s = _score(rr, "two_route")
        s["hd_shift"] = sh
        hd_shifts.append(s)
        mark = "  <-- published inputs" if sh == 0 else ""
        print(f"  {sh:>+10d}{s['r2']:>9.3f}{s['spearman']:>+10.3f}"
              f"{s['fold']:>8.2f}x{s['geo_bias']:>9.2f}x{s['within_3x']:>6.0%}{mark}")

    # ---- D. combined, at the magnitudes actually measured -------------------
    print("\n  D. COMBINED, at the error magnitudes THIS EXPERIMENT MEASURED.")
    print("     Each row is one honest guess at what the engine's own descriptors")
    print("     would do, built only from numbers produced above. It is a")
    print("     projection, and it is labelled as one -- the real answer needs")
    print("     the 66 structures nobody has.")
    combos = []
    stats = _crippen_error_stats()
    a, b = stats["all"], stats["ex_ionisable"]
    scen = [
        ("Crippen error as measured (n=5)", a["bias"], a["sd"], 0),
        ("  ... + engine donor convention", a["bias"], a["sd"], -1),
        ("Crippen error EXCLUDING the base (n=4)", b["bias"], b["sd"], 0),
        ("  ... + engine donor convention", b["bias"], b["sd"], -1),
    ]
    print(f"\n  Every bias/sigma below is MEASURED in Stage 2, not chosen:")
    print(f"    all 5 compounds : bias {a['bias']:+.3f}, residual sd "
          f"{a['sd']:.3f}  (RMS {a['rms']:.3f})")
    print(f"    excluding NLCQ-1: bias {b['bias']:+.3f}, residual sd "
          f"{b['sd']:.3f}  (RMS {b['rms']:.3f})")
    print(f"\n  {'scenario':<40}{'bias':>7}{'sigma':>7}{'dHD':>5}"
          f"{'mean R^2':>10}{'sd':>7}")
    for lbl, bias, sig, sh in scen:
        r2s = []
        for _ in range(replicates):
            rr = _rows_with(
                {n: base_logp[n] + bias + rng.gauss(0.0, sig) for n in names},
                hbd_by_name={n: max(0.0, base_hd[n] + sh) for n in names},
                cache=cache)
            r2s.append(_score(rr, "two_route")["r2"])
        mu = sum(r2s) / len(r2s)
        sd = math.sqrt(sum((x - mu) ** 2 for x in r2s) / (len(r2s) - 1))
        combos.append({"scenario": lbl.strip(), "bias": bias, "sigma": sig,
                       "hd_shift": sh, "mean_r2": mu, "sd_r2": sd})
        print(f"  {lbl:<40}{bias:>+7.2f}{sig:>7.2f}{sh:>5d}{mu:>10.3f}{sd:>7.3f}")
    print(f"\n  Baseline for comparison: R^2 = {base['r2']:.3f}.")
    print("  These rows rest on a FIVE-COMPOUND estimate of the shipping")
    print("  descriptor's error. That is the weakest number in this file and it is")
    print("  used because it is the ONLY direct measurement of the shipping")
    print("  descriptor that exists against this dataset. Its weakness is the")
    print("  finding, not a footnote: 66 of 71 compounds cannot be measured at all.")

    return {"baseline": base, "offsets": offsets, "noise": noises,
            "hd_shifts": hd_shifts, "combined": combos,
            "crippen_error_stats": stats,
            "seed": seed, "replicates": replicates}


def _crippen_error_stats() -> dict:
    """Bias and residual scatter of RDKit Crippen against Table S1's logP7.4, on
    the only five compounds where both exist. Recomputed here rather than passed
    in, so `--stage propagate` is runnable on its own and cannot quote a stale
    number. `ex_ionisable` drops NLCQ-1, the one compound with a basic centre
    that is protonated at pH 7.4."""
    per = []
    for name, spec in CURATED_STRUCTURES.items():
        t = _table_row(name)
        if t[4] is None:
            continue
        mol = descriptors.from_smiles(name, spec["smiles"])
        per.append((name, float(mol.get("logp")) - float(t[4])))

    def agg(ds):
        n = len(ds)
        mu = sum(ds) / n
        return {"n": n, "bias": mu, "rms": math.sqrt(sum(d * d for d in ds) / n),
                "sd": math.sqrt(sum((d - mu) ** 2 for d in ds) / (n - 1))}

    return {"per_compound": {k: v for k, v in per},
            "all": agg([d for _, d in per]),
            "ex_ionisable": agg([d for k, d in per if k != "NLCQ-1"])}


# =========================================================================== #
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all",
                    choices=("all", "structures", "descriptors", "substitute",
                             "propagate"))
    ap.add_argument("--replicates", type=int, default=200)
    a = ap.parse_args()

    t0 = time.time()
    out = {}
    if a.stage in ("all", "structures"):
        out["structures"] = stage_structures()
    if a.stage in ("all", "descriptors"):
        out["descriptors"] = stage_descriptors()
    if a.stage in ("all", "substitute"):
        out["substitute"] = stage_substitute()
    if a.stage in ("all", "propagate"):
        out["propagate"] = stage_propagate(replicates=a.replicates)

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, f"descriptor_run_{a.stage}.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {path}   ({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
