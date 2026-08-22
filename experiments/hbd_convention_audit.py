"""HBD CONVENTION AUDIT -- donor HEAVY ATOMS vs donor HYDROGENS.

MEASUREMENT ONLY. This script changes NOTHING in ``tmesim/``. It monkey-patches
``descriptors.compute_descriptors`` IN-PROCESS for the "after" arm and restores
it, so the repo default (``Lipinski.NumHDonors``) is never modified on disk.

The question
------------
``two_route.membrane_permeability_um_s`` applies
``HBD_LOGP_DECREMENT = 1.34`` log units per hydrogen-bond donor. Two different
integer counts are in play:

  * ``Lipinski.NumHDonors``  -- donor HEAVY ATOMS  (-NH2 counts 1). THE ENGINE.
  * ``Descriptors.NHOHCount`` -- donor HYDROGENS   (-NH2 counts 2). PRUIJN's HD.

The law is SCORED against Pruijn 2008 Table S1's ``HD`` column
(``experiments/diffusivity_falsification.py::PRUIJN_2008_TABLE_S1``, field 3,
which is fed straight into ``two_route`` by
``experiments/two_route_transport.py::_build_rows``) and is FED
``Lipinski.NumHDonors`` at runtime (``tmesim/molecular/descriptors.py:87``).

Stages
------
``--stage convention``  Pruijn's own stated definition + per-compound check on
                        every Table S1 row whose structure is resolvable.
``--stage library``     ``hbd`` under both conventions for every built-in drug.
``--stage alpha``       ``calibrate_alpha()`` before and after.
``--stage pruijn``      Pruijn R^2 baseline and the HD-shift bracket.
``--stage validate``    ``codon validate`` verdicts before and after.

Run:
    ./.venv/bin/python experiments/hbd_convention_audit.py --stage convention
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Tuple

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from rdkit import Chem                                        # noqa: E402
from rdkit.Chem import Descriptors, Lipinski                  # noqa: E402

from tmesim.molecular import descriptors as D                 # noqa: E402


# --------------------------------------------------------------------------- #
# The only structures in Pruijn 2008 Table S1 that resolve to a SMILES.
#
# Table S1 carries NO SMILES and NO CAS numbers: 67 of the 71 rows are Auckland
# Cancer Society Research Centre `SN` codes whose structures appear only as
# embedded 2-D drawings (data_cache/pruijn2008/PROVENANCE.md, "Structural
# identifiers -- a real limitation"). The four non-SN rows are named compounds;
# TPZ is named in its own row label. Nothing here was guessed: each SMILES is a
# textbook structure for a named compound.
#
# name-in-table -> (SMILES, why we trust the structure)
# --------------------------------------------------------------------------- #
RESOLVABLE: Dict[str, Tuple[str, str]] = {
    "Urea": ("NC(N)=O", "urea, CH4N2O"),
    "Mannitol": ("OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", "D-mannitol, C6H14O6"),
    "SN24280 (TPZ)": (
        "Nc1n[n+]([O-])c2ccccc2[n+]1[O-]",
        "tirapazamine = 3-amino-1,2,4-benzotriazine 1,4-dioxide, C7H6N4O2"),
    "Acridan-9(10H)-one": (
        "O=C1c2ccccc2Nc2ccccc21", "acridin-9(10H)-one, C13H9NO"),
    "NLCQ-1": (
        "Clc1ccc2c(NCCCn3ccnc3[N+](=O)[O-])ccnc2c1",
        "NLCQ-1 = 4-[3-(2-nitro-1-imidazolyl)propylamino]-7-chloroquinoline, "
        "C15H14ClN5O2"),
}

#: Pruijn's own definition, verbatim from the article's "Physicochemical
#: Parameters" section (Aust J Chem 2008;61(9):687-693, doi:10.1071/CH08240;
#: retrieved from the Internet Archive copy of the publisher's own open-access
#: PDF, the same route recorded in data_cache/pruijn2008/PROVENANCE.md).
PRUIJN_DEFINITION = (
    "MW, HD (count of -NH and -OH groups), HA (count of N and O atoms), and "
    "apparent (macroscopic) pKa values for the side chain were calculated "
    "using ACD/PhysChem package (version 8.0, Advanced Chemistry Development "
    "Inc., Toronto, Canada)."
)


# --------------------------------------------------------------------------- #
# The patch. Donor HYDROGENS instead of donor HEAVY ATOMS.
# --------------------------------------------------------------------------- #
_ORIG_COMPUTE = D.compute_descriptors


def _compute_descriptors_nhoh(smiles: str) -> Dict[str, float]:
    """``compute_descriptors`` with ``hbd = NHOHCount`` instead of NumHDonors."""
    out = _ORIG_COMPUTE(smiles)
    mol = Chem.MolFromSmiles(smiles)
    out["hbd"] = int(Descriptors.NHOHCount(mol))
    return out


#: Fallback-table rows whose ``hbd`` would also have to move, computed from the
#: library SMILES rather than typed. Filled by :func:`_fallback_overrides`.
def _fallback_overrides() -> Dict[str, int]:
    from tmesim.molecular.library import GBM_DRUGS
    from tmesim.molecular.breast_library import BREAST_DRUGS
    smi = dict(GBM_DRUGS)
    smi.update(BREAST_DRUGS)
    out = {}
    for name in D._FALLBACK_TABLE:
        s = smi.get(name)
        if not s:
            continue
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        out[name] = int(Descriptors.NHOHCount(m))
    return out


class nhoh_convention:
    """Context manager: engine-wide ``hbd`` = donor hydrogens, then restore."""

    def __enter__(self):
        self._saved_table = {k: dict(v) for k, v in D._FALLBACK_TABLE.items()}
        D.compute_descriptors = _compute_descriptors_nhoh
        for name, n in _fallback_overrides().items():
            D._FALLBACK_TABLE[name]["hbd"] = n
        return self

    def __exit__(self, *exc):
        D.compute_descriptors = _ORIG_COMPUTE
        for k, v in self._saved_table.items():
            D._FALLBACK_TABLE[k] = v
        return False


# --------------------------------------------------------------------------- #
# Stage 1: the convention
# --------------------------------------------------------------------------- #
def stage_convention() -> None:
    import csv

    print("=" * 78)
    print("STAGE 1 -- WHAT IS PRUIJN'S `HD`?")
    print("=" * 78)
    print("\nThe paper's own Methods, VERBATIM:\n")
    print("  \"" + PRUIJN_DEFINITION + "\"\n")
    print("  Aust J Chem 2008;61(9):687-693, section 'Physicochemical")
    print("  Parameters'. doi:10.1071/CH08240.")
    print("\nThe HA half of that sentence is DIRECTLY CHECKABLE ('count of N and")
    print("O atoms'), which is why it is quoted here: if HA reproduces on the")
    print("resolvable rows, the same sentence's HD clause is trustworthy.\n")

    path = os.path.join(_ROOT, "data_cache", "pruijn2008",
                        "pruijn2008_tableS1.csv")
    with open(path) as fh:
        table = list(csv.DictReader(fh))
    print(f"Table S1 rows: {len(table)}")
    resolvable = [r for r in table if r["compound_id"] in RESOLVABLE]
    print(f"Rows whose structure is resolvable to a SMILES: {len(resolvable)}")
    print(f"Rows that are ACSRC `SN` codes with drawings only: "
          f"{len(table) - len(resolvable)}")
    print("\n*** THE CONCLUSION BELOW RESTS ON "
          f"{len(resolvable)} OF {len(table)} ROWS. ***\n")

    hdr = (f"{'compound':22s} {'paper HD':>8s} {'NumHDonors':>11s} "
           f"{'NHOHCount':>10s} {'paper HA':>8s} {'N+O atoms':>10s}  verdict")
    print(hdr)
    print("-" * len(hdr))
    agree_atoms = agree_h = agree_ha = 0
    discriminating = []
    for r in resolvable:
        name = r["compound_id"]
        smi = RESOLVABLE[name][0]
        m = Chem.MolFromSmiles(smi)
        nd = int(Lipinski.NumHDonors(m))
        nh = int(Descriptors.NHOHCount(m))
        no = sum(1 for a in m.GetAtoms() if a.GetSymbol() in ("N", "O"))
        phd = int(r["h_bond_donors"])
        pha = int(r["h_bond_acceptors"])
        agree_atoms += (nd == phd)
        agree_h += (nh == phd)
        agree_ha += (no == pha)
        if nd != nh:
            discriminating.append(name)
            verdict = ("HYDROGENS" if phd == nh else
                       "HEAVY ATOMS" if phd == nd else "NEITHER")
        else:
            verdict = "(cannot discriminate)"
        print(f"{name:22s} {phd:8d} {nd:11d} {nh:10d} {pha:8d} {no:10d}  "
              f"{verdict}")
    print("-" * len(hdr))
    n = len(resolvable)
    print(f"paper HD == NumHDonors (heavy atoms): {agree_atoms}/{n}")
    print(f"paper HD == NHOHCount  (hydrogens)  : {agree_h}/{n}")
    print(f"paper HA == count of N and O atoms  : {agree_ha}/{n}   "
          "<- checks the quoted sentence itself")
    print(f"\nRows where the two conventions DIFFER at all (i.e. rows that carry")
    print(f"any information about which convention it is): "
          f"{len(discriminating)}  {discriminating}")

    print("\nWhat the paper's own text says about the HD RANGE (a second, "
          "independent check):")
    print("  Study A: 'with 0-3 HD and 5-11 HA'")
    print("  Study B: 'range of 60-347, HA of 2-7, and HD of 1-6'")
    hd_all = sorted({int(r["h_bond_donors"]) for r in table})
    print(f"  Table S1 HD values actually present: {hd_all}")


# --------------------------------------------------------------------------- #
# Stage 2: the library
# --------------------------------------------------------------------------- #
def stage_library() -> None:
    from tmesim.molecular.library import GBM_DRUGS
    from tmesim.molecular.breast_library import BREAST_DRUGS

    print("=" * 78)
    print("STAGE 2 -- BUILT-IN LIBRARY: WHICH COMPOUNDS MOVE?")
    print("=" * 78)
    smi: Dict[str, str] = {}
    tag: Dict[str, str] = {}
    for t, dd in (("GBM", GBM_DRUGS), ("breast", BREAST_DRUGS)):
        for k, v in dd.items():
            if k in smi:
                tag[k] = "both"
            else:
                smi[k], tag[k] = v, t
    print(f"\n{len(smi)} unique compounds across GBM_DRUGS + BREAST_DRUGS.\n")
    hdr = (f"{'compound':20s} {'panel':7s} {'NumHDonors':>10s} {'NHOHCount':>9s} "
           f"{'delta':>6s}   {'d(logK_mem)':>11s}  {'x K_mem':>9s}")
    print(hdr)
    print("-" * len(hdr))
    from tmesim.molecular import two_route
    dec = two_route.HBD_LOGP_DECREMENT
    moved = unmoved = unparseable = 0
    for name in sorted(smi):
        m = Chem.MolFromSmiles(smi[name])
        if m is None:
            print(f"{name:20s} {tag[name]:7s} {'n/a':>10s} {'n/a':>9s} "
                  f"{'--':>6s}   (no SMILES; fallback table)")
            unparseable += 1
            continue
        nd, nh = int(Lipinski.NumHDonors(m)), int(Descriptors.NHOHCount(m))
        d = nh - nd
        if d:
            moved += 1
        else:
            unmoved += 1
        print(f"{name:20s} {tag[name]:7s} {nd:10d} {nh:9d} {d:+6d}   "
              f"{-dec * d:11.2f}  {10 ** (-dec * d):9.4f}")
    print("-" * len(hdr))
    print(f"CHANGED: {moved}    UNCHANGED: {unmoved}    "
          f"NO SMILES: {unparseable}")
    print(f"\n(d(logK_mem) is the shift in log10 membrane/water partition "
          f"coefficient\n at HBD_LOGP_DECREMENT = {dec}, BEFORE the "
          f"logk_max = {two_route.LOGK_MEMBRANE_MAX} clamp.)")


# --------------------------------------------------------------------------- #
# Stage 3: alpha
# --------------------------------------------------------------------------- #
def stage_alpha() -> None:
    print("=" * 78)
    print("STAGE 3 -- calibrate_alpha() BEFORE AND AFTER")
    print("=" * 78)
    from tmesim.validation.calibration import calibrate_alpha
    import tmesim.validation.calibration as CAL

    def _fresh():
        # calibrate_alpha may memoise; clear anything that looks like a cache.
        for attr in dir(CAL):
            o = getattr(CAL, attr)
            if hasattr(o, "cache_clear"):
                o.cache_clear()

    _fresh()
    before = float(calibrate_alpha())
    print(f"\nbefore (Lipinski.NumHDonors, the shipping default): {before!r}")
    with nhoh_convention():
        _fresh()
        after = float(calibrate_alpha())
    _fresh()
    again = float(calibrate_alpha())
    print(f"after  (Descriptors.NHOHCount)                    : {after!r}")
    print(f"restored default re-run (must equal 'before')     : {again!r}")
    assert again == before, "monkeypatch did not restore cleanly"
    print(f"\ndelta = {after - before:+.6f}   ratio = {after / before:.6f}")
    print("\ncalibrate_alpha() is fit on DOXORUBICIN, whose donor count moves.")


# --------------------------------------------------------------------------- #
# Stage 4: Pruijn
# --------------------------------------------------------------------------- #
def stage_pruijn() -> None:
    print("=" * 78)
    print("STAGE 4 -- THE PRUIJN SCORE, AND THE HD-SHIFT BRACKET")
    print("=" * 78)
    import experiments.two_route_transport as TRT
    from experiments.diffusivity_falsification import PRUIJN_2008_TABLE_S1
    from tmesim.molecular import two_route

    print("\nNote the direction of the question. `_build_rows` feeds the PAPER'S")
    print("HD column into the law, so switching the ENGINE's convention does not")
    print("move this number by itself. What moves it is asking what the score")
    print("WOULD have been had the law been scored on the count the engine")
    print("actually computes -- i.e. HD shifted down wherever an -NH2/-NH3+ is")
    print("present. Structures are unavailable for 66/71 rows, so this is a")
    print("BRACKET, not a measurement.\n")

    cache: Dict[float, float] = {}
    base = [r for r in TRT._build_rows(two_route.HBD_LOGP_DECREMENT, cache)
            if "two_route" in r]
    r2_base = TRT._r2(base, "two_route") if hasattr(TRT, "_r2") else None
    print(f"baseline rows scored: {len(base)}")
    print(f"baseline R^2 (paper HD, as published): {_score(TRT, base)}")

    orig = list(PRUIJN_2008_TABLE_S1)
    try:
        for shift in (-1, -2):
            PRUIJN_2008_TABLE_S1[:] = [
                (t[0], t[1], max(0, t[2] + shift)) + tuple(t[3:]) for t in orig]
            rows = [r for r in TRT._build_rows(two_route.HBD_LOGP_DECREMENT,
                                               dict(cache))
                    if "two_route" in r]
            print(f"HD {shift:+d} on every row (floored at 0): "
                  f"R^2 = {_score(TRT, rows)}")
        # Structure-aware: shift ONLY the rows whose structure we can resolve
        # and whose two counts actually differ.
        known = {}
        for name, (s, _why) in RESOLVABLE.items():
            m = Chem.MolFromSmiles(s)
            known[name] = int(Lipinski.NumHDonors(m)) - int(
                Descriptors.NHOHCount(m))

        print("\nWhich rows CAN move at all, with no structures needed:")
        n_hd01 = sum(1 for t in orig if t[2] <= 1)
        print(f"  HD <= 1 rows: {n_hd01} of {len(orig)}. A row with HD = 1 has")
        print("  exactly ONE N-H/O-H hydrogen, so its donor group holds one")
        print("  hydrogen, so NumHDonors == NHOHCount == 1. HD = 0 likewise.")
        print(f"  These {n_hd01} rows are PROVABLY UNAFFECTED by the convention.")
        print(f"  Only the {len(orig) - n_hd01} rows with HD >= 2 can move, and")
        print("  only if a donor group carries more than one hydrogen.")

        def scenario(label, fn):
            PRUIJN_2008_TABLE_S1[:] = [
                (t[0], t[1], max(0, t[2] + (known[t[0]] if t[0] in known
                                            else fn(t[2])))) + tuple(t[3:])
                for t in orig]
            rr = [r for r in TRT._build_rows(two_route.HBD_LOGP_DECREMENT,
                                             dict(cache))
                  if "two_route" in r]
            print(f"{label}\n  {_score(TRT, rr)}")

        scenario("A) ONLY the 2 resolvable rows that differ (Urea -2, TPZ -1); "
                 "every\n   unresolved row assumed UNAFFECTED -- the "
                 "least-damage end:",
                 lambda hd: 0)
        scenario("B) every HD >= 2 row carries ONE -NH2 (TPZ's own 3-amino "
                 "group,\n   retained by its analogues) -> -1; resolvable rows "
                 "pinned to truth:",
                 lambda hd: -1 if hd >= 2 else 0)
        scenario("C) every HD >= 2 row is ALL -NH2 -> -(HD // 2); resolvable "
                 "rows\n   pinned to truth -- the most-damage end:",
                 lambda hd: -(hd // 2))
    finally:
        PRUIJN_2008_TABLE_S1[:] = orig
    rows = [r for r in TRT._build_rows(two_route.HBD_LOGP_DECREMENT, dict(cache))
            if "two_route" in r]
    print(f"\nrestored baseline re-run (must equal baseline): "
          f"{_score(TRT, rows)}")


def _score(TRT, rows) -> str:
    """The PUBLISHED score. Delegates to ``two_route_transport._score``, which
    is the function that produced 0.718 -- squared log10 Pearson, plus the
    geometric-bias and within-3x figures the promotion decision also used."""
    s = TRT._score(rows, "two_route")
    return (f"R^2 = {s['r2']:.4f}  (n={s['n']}, "
            f"geo bias {s['geo_bias']:.2f}x, "
            f"within-3x {s['within_3x']:.0%}, "
            f"spread {s['fold']:.0f}x)")


# --------------------------------------------------------------------------- #
# Stage 5: codon validate
# --------------------------------------------------------------------------- #
def stage_validate(quick: bool, only: Optional[str]) -> None:
    print("=" * 78)
    print("STAGE 5 -- EVERY VERDICT IN `codon validate`, BEFORE AND AFTER")
    print("=" * 78)
    from tmesim.validation import scorecard as SC

    keys = [k.strip() for k in only.split(",")] if only else None

    def run(label):
        print(f"\n[{label}] running suite (quick={quick}) ...", flush=True)
        return SC.run_suite(quick=quick, only=keys,
                            on_progress=lambda g: print(
                                f"  [{label}] {g.key}", flush=True))

    before = run("before")
    with nhoh_convention():
        after = run("after")

    def index(res):
        out = {}
        for i, r in enumerate(res["rows"]):
            out[(r["group"], r.get("label", ""), i)] = r
        return out
    b, a = index(before), index(after)
    print()
    hdr = (f"{'group':30s} {'row':22s} {'before':13s} {'after':13s} "
           f"{'pred before':>22s} {'pred after':>22s}  flipped?")
    print(hdr)
    print("-" * len(hdr))
    flips = []
    moved = []
    for k in b:
        rb, ra = b[k], a.get(k, {})
        vb, va = str(rb.get("status")), str(ra.get("status"))
        pb, pa = str(rb.get("predicted")), str(ra.get("predicted"))
        f = "*** FLIPPED ***" if vb != va else ("value moved" if pb != pa else "")
        if vb != va:
            flips.append((k[0], k[1], vb, va))
        elif pb != pa:
            moved.append((k[0], k[1], pb, pa))
        print(f"{k[0][:30]:30s} {str(k[1])[:22]:22s} {vb:13s} {va:13s} "
              f"{pb[:22]:>22s} {pa[:22]:>22s}  {f}")
    print("-" * len(hdr))
    print(f"rows whose PREDICTED VALUE moved but verdict held: {len(moved)}")
    for g, l, pb, pa in moved:
        print(f"  {g} / {l}: {pb}  ->  {pa}")
    print(f"before counts: {before['counts']}   exit={before['exit_code']}")
    print(f"after  counts: {after['counts']}   exit={after['exit_code']}")
    print(f"\nVERDICT FLIPS: {len(flips)}")
    for g, l, vb, va in flips:
        print(f"  {g} / {l}: {vb} -> {va}")

    out = os.path.join(_ROOT, "results", "hbd_convention_validate.json")
    with open(out, "w") as fh:
        json.dump({"before": before, "after": after,
                   "flips": flips, "quick": quick}, fh, indent=2, default=str)
    print(f"\nfull JSON: {out}")


# --------------------------------------------------------------------------- #
# Stage 6: is the score RECOMPUTABLE under the engine's own convention?
# --------------------------------------------------------------------------- #
def stage_recompute(trials: int) -> None:
    """Bound the engine-convention Pruijn R^2 WITHOUT transcribing a structure.

    The claim being tested is "0.718 is not recomputable, because 66 of 71
    structures are unavailable". That claim is too strong, and this stage says
    by how much. Two facts do the work:

      1. `HD <= 1` pins the heavy-atom count EXACTLY (a group with one hydrogen
         has one heavy atom), with no structure.
      2. `HD = h` bounds the heavy-atom count to `ceil(h/2) .. h` for the rest,
         because no organic donor group carries more than two donor hydrogens
         on one heavy atom (-NH2, -NH3+ being the outlier, handled separately).

    So every row's heavy-atom count lies in a SMALL known set. We enumerate the
    prediction at each candidate, then sample the joint space for min/max R^2.
    """
    import random

    import experiments.two_route_transport as TRT
    from experiments.diffusivity_falsification import PRUIJN_2008_TABLE_S1
    from tmesim.molecular import two_route

    print("=" * 78)
    print("STAGE 6 -- IS THE PRUIJN SCORE RECOMPUTABLE ON THE ENGINE'S "
          "CONVENTION?")
    print("=" * 78)

    orig = list(PRUIJN_2008_TABLE_S1)
    resolved = {}
    for name, (s, _why) in RESOLVABLE.items():
        m = Chem.MolFromSmiles(s)
        resolved[name] = int(Lipinski.NumHDonors(m))

    # Candidate heavy-atom counts per row.
    cands: List[Tuple[str, List[int]]] = []
    for t in orig:
        name, hd = t[0], t[2]
        if name in resolved:
            cands.append((name, [resolved[name]]))
        elif hd <= 1:
            cands.append((name, [hd]))
        else:
            lo = (hd + 1) // 2          # every donor group carries two H
            cands.append((name, list(range(lo, hd + 1))))
    n_pinned = sum(1 for _, c in cands if len(c) == 1)
    print(f"\nrows whose heavy-atom count is PINNED with no structure: "
          f"{n_pinned} of {len(orig)}")
    print(f"  ({sum(1 for t in orig if t[2] <= 1)} by the HD <= 1 argument, "
          f"{len(resolved)} by a resolvable structure)")
    free = [(n, c) for n, c in cands if len(c) > 1]
    print(f"rows left genuinely free: {len(free)}, each over a set of size "
          f"{sorted({len(c) for _, c in free})}")
    total = 1
    for _n, c in free:
        total *= len(c)
    print(f"joint space: {total:.3g} configurations")

    # Cache the prediction for every (row, candidate hbd) pair ONCE.
    print("\nprecomputing predictions for every (row, candidate) pair ...")
    pred: Dict[Tuple[str, int], float] = {}
    meas: Dict[str, float] = {}
    dcache: Dict[float, float] = {}
    allc = sorted({h for _n, c in cands for h in c})
    for h in allc:
        PRUIJN_2008_TABLE_S1[:] = [
            (t[0], t[1], h) + tuple(t[3:]) for t in orig]
        for r in TRT._build_rows(two_route.HBD_LOGP_DECREMENT, dcache):
            if "two_route" not in r:
                continue
            pred[(r["name"], h)] = r["two_route"]
            meas[r["name"]] = r["meas"]
    PRUIJN_2008_TABLE_S1[:] = orig
    names = [n for n, _c in cands if n in meas]
    print(f"  {len(pred)} predictions cached over {len(names)} scorable rows")

    import math

    def r2(assign: Dict[str, int]) -> float:
        lp = [math.log10(pred[(n, assign[n])]) for n in names]
        lm = [math.log10(meas[n]) for n in names]
        return TRT._pearson(lp, lm) ** 2

    base = {n: dict(cands)[n][0] for n in names}   # not meaningful alone
    lo_assign = {n: max(dict(cands)[n]) for n in names}   # fewest donors removed
    hi_assign = {n: min(dict(cands)[n]) for n in names}
    a_hi, a_lo = r2(lo_assign), r2(hi_assign)
    print(f"\nextreme configuration 'every free row keeps HD heavy atoms': "
          f"R^2 = {a_hi:.4f}")
    print(f"extreme configuration 'every free row is all -NH2'          : "
          f"R^2 = {a_lo:.4f}")

    rng = random.Random(20260821)
    lo = hi = None
    for _ in range(trials):
        asn = {n: rng.choice(dict(cands)[n]) for n in names}
        v = r2(asn)
        lo = v if lo is None else min(lo, v)
        hi = v if hi is None else max(hi, v)
    print(f"\n{trials} random configurations sampled from the joint space:")
    print(f"  min R^2 = {lo:.4f}   max R^2 = {hi:.4f}")
    print("  (R^2 is NOT monotone in the per-row shift, so the two 'extreme'")
    print("   configurations above are NOT the extremes of the interval.)")

    # Coordinate ascent/descent from random restarts. 26 binary choices, so
    # this finds the true extremes to all practical purposes.
    cmap = dict(cands)

    def climb(sign: int) -> float:
        best = None
        for _ in range(40):
            asn = {n: rng.choice(cmap[n]) for n in names}
            cur = r2(asn)
            improved = True
            while improved:
                improved = False
                for n in names:
                    if len(cmap[n]) == 1:
                        continue
                    keep = asn[n]
                    for alt in cmap[n]:
                        if alt == keep:
                            continue
                        asn[n] = alt
                        v = r2(asn)
                        if sign * v > sign * cur:
                            cur, keep, improved = v, alt, True
                    asn[n] = keep
            best = cur if best is None else (max(best, cur) if sign > 0
                                             else min(best, cur))
        return best

    hi_opt, lo_opt = climb(+1), climb(-1)
    print(f"\ncoordinate ascent/descent, 40 restarts each:")
    print(f"  TRUE-ish min R^2 = {lo_opt:.4f}   TRUE-ish max R^2 = {hi_opt:.4f}")
    lo, hi = min(lo, lo_opt), max(hi, hi_opt)
    print(f"\nPUBLISHED (paper HD, i.e. donor HYDROGENS): "
          f"{_score(TRT, [r for r in TRT._build_rows(two_route.HBD_LOGP_DECREMENT, dcache) if 'two_route' in r])}")
    print("\nCONCLUSION: the score is NOT 'unrecomputable'. It is recomputable")
    print(f"to an INTERVAL of width {abs(hi - lo):.3f} in R^2, without")
    print("transcribing a single structure from a PDF drawing.")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=("convention", "library", "alpha", "pruijn",
                             "validate", "recompute"))
    ap.add_argument("--trials", type=int, default=20000)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    {"convention": stage_convention,
     "library": stage_library,
     "alpha": stage_alpha,
     "pruijn": stage_pruijn,
     "validate": lambda: stage_validate(args.quick, args.only),
     "recompute": lambda: stage_recompute(args.trials)}[args.stage]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
