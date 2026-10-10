"""Check the written record against the data it claims to come from.

The rule this study is held to is "never invent a number". That is easy to
assert and easy to violate by transcription, so it is checked mechanically:

1.  No unfilled placeholder survives (``TODO``, ``TBD``, ``XXX``, ``<...>``,
    ``FIXME``, ``???``).
2.  Every decimal number of 4 or more significant digits in the report appears
    verbatim in one of the study's data files, or in a repository file the
    report cites. Four digits is the threshold because short numbers (``38``,
    ``1.27``, section counts, years, section numbers) are structural or rounded
    and are checked by eye; a long decimal is a measurement and must be
    traceable.

Run::

    ./.venv/bin/python experiments/surrogate_feasibility/check_report.py \
        results/surrogate_feasibility_2026-10-09.md

Exits non-zero on any violation and names it.
"""
from __future__ import annotations

import os
import re
import sys
from typing import List, Set

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, "..", ".."))

#: Where a number in the report is allowed to have come from.
SOURCES = [
    os.path.join(_HERE, "sweep_data.json"),
    os.path.join(_HERE, "fit_results.json"),
    os.path.join(_HERE, "analyse_results.json"),
    os.path.join(_HERE, "derived.json"),
    os.path.join(_HERE, "coarse_baseline.json"),
    os.path.join(_HERE, "ladder_data.json"),
    os.path.join(_HERE, "sweep.log"),
    os.path.join(_HERE, "ladder.log"),
    os.path.join(_HERE, "coarse.log"),
    os.path.join(_HERE, "fit.log"),
    # Repository records the report quotes rather than re-measures.
    os.path.join(_ROOT, "README.md"),
    os.path.join(_ROOT, "app", "server", "reach.py"),
    os.path.join(_ROOT, "results", "multisection_generalisation.json"),
    os.path.join(_ROOT, "results", "dx_sensitivity_steady_data.json"),
    os.path.join(_ROOT, "results", "delivery_spread_noise.md"),
    os.path.join(_ROOT, "results", "performance_profile_2026-09-20.md"),
    os.path.join(_ROOT, "results", "subgrid_closure.md"),
    os.path.join(_ROOT, "results", "CLAIM_INVENTORY.md"),
    os.path.join(_ROOT, "results", "verdict_inversion_census.md"),
    os.path.join(_ROOT, "tmesim", "ml", "panel_imputation.py"),
    os.path.join(_ROOT, "results", "hypoxia_imputation_model.md"),
]

PLACEHOLDERS = ("TODO", "TBD", "XXX", "FIXME", "???", "<FILL", "<NUM",
                "PLACEHOLDER")

#: A decimal with >=4 significant digits, e.g. 1.0004, 22.2234, 0.0037313.
_NUM = re.compile(r"\d+\.\d{3,}")

#: Any numeric literal in a source file, for the rounding check.
_SRC_NUM = re.compile(r"\d+\.\d+")


def _rounds_from_source(tok: str, src_nums) -> bool:
    """Is ``tok`` a correct rounding of some number in the sources?

    A report may quote a measured ``0.7132559626060712`` as ``0.7133``. That is
    traceable and must pass; a typo of ``0.7134`` must not. Truncation is also
    accepted (handled by the prefix test in the caller), so this only has to
    catch the round-half cases.
    """
    try:
        want = float(tok)
    except ValueError:
        return False
    places = len(tok.split(".", 1)[1])
    for s in src_nums:
        try:
            if round(float(s), places) == want:
                return True
        except (ValueError, OverflowError):
            continue
    return False


def load_sources() -> str:
    blobs: List[str] = []
    for p in SOURCES:
        if os.path.exists(p):
            with open(p, errors="replace") as fh:
                blobs.append(fh.read())
    return "\n".join(blobs)


def main(argv: List[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    report_path = argv[1]
    with open(report_path) as fh:
        report = fh.read()

    problems: List[str] = []

    for ph in PLACEHOLDERS:
        for i, line in enumerate(report.splitlines(), 1):
            if ph in line:
                problems.append(f"line {i}: placeholder {ph!r}: "
                                f"{line.strip()[:90]}")

    hay = load_sources()
    # Candidate source numbers, deduplicated, for the rounding check. Only
    # those sharing the token's integer part are worth testing, so the scan
    # stays cheap even on a few MB of JSON.
    src_nums = set(_SRC_NUM.findall(hay))
    by_int: dict = {}
    for s in src_nums:
        by_int.setdefault(s.split(".", 1)[0], []).append(s)

    seen: Set[str] = set()
    untraceable: List[str] = []
    for m in _NUM.finditer(report):
        tok = m.group(0)
        if tok in seen:
            continue
        seen.add(tok)
        if tok in hay:
            continue
        # A TRUNCATED quote of a longer measured value is legitimate: accept a
        # token that is a prefix of something in the sources.
        if re.search(re.escape(tok) + r"\d", hay):
            continue
        # So is a correctly ROUNDED quote.
        if _rounds_from_source(tok, by_int.get(tok.split(".", 1)[0], [])):
            continue
        untraceable.append(tok)
    if untraceable:
        problems.append(
            f"{len(untraceable)} long decimals do not appear in any data "
            f"source: {sorted(untraceable)[:20]}")

    print(f"[check] {report_path}")
    print(f"[check] distinct long decimals checked: {len(seen)}")
    print(f"[check] sources present: "
          f"{sum(1 for p in SOURCES if os.path.exists(p))}/{len(SOURCES)}")
    if problems:
        print("\nPROBLEMS:")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("[check] OK — no placeholders, every long decimal traces to data")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
