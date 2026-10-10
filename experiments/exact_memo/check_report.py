#!/usr/bin/env python3
"""Re-verify results/exact_memo_2026-10-10.md against this study's data files.

Same convention as ``experiments/surrogate_feasibility/check_report.py``: every
decimal of 4 or more significant digits in the record must appear VERBATIM in
one of this study's JSON files, in a source file in the repository, or in a
repository record the text cites. Placeholders fail.

A record whose numbers nobody re-reads is how
``results/codon_claim-surface-drift`` happens: the engine stays sound and the
document stops being true.
"""
from __future__ import annotations

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))

DATA = ("probe_scaling.json", "probe_classification.json",
        "axis_exponents.json", "speedup_dx30.json", "speedup_dx12.json",
        "verify_against_solve.json", "plants.json")

#: Files whose live text may also substantiate a figure.
SOURCES = ("app/server/memo.py", "app/server/memo_api.py", "app/server/jobs.py",
           "app/server/reach.py", "tests/test_memo.py",
           "results/surrogate_feasibility_2026-10-09.md",
           "results/performance_profile_2026-09-20.md")

PLACEHOLDERS = ("TODO", "TBD", "XXX", "FIXME", "<fill", "???")

#: A decimal with >= 4 significant digits, as written.
NUM = re.compile(r"\d[\d.eE+\-]*")

#: ISO dates are not measurements. Stripped before the scan so a dated
#: filename or a withdrawal date is not demanded of a data file.
DATE = re.compile(r"(?<!\d)(19|20)\d\d-\d\d-\d\d(?!\d)")


def significant(token: str) -> int:
    body = token.lower().split("e")[0].lstrip("-+")
    digits = body.replace(".", "").lstrip("0")
    return len(digits.rstrip("0")) if "." in body else len(digits)


def haystack() -> str:
    parts = []
    for name in DATA:
        path = os.path.join(HERE, name)
        if not os.path.isfile(path):
            print(f"MISSING DATA FILE {name}", file=sys.stderr)
            continue
        with open(path) as fh:
            blob = json.load(fh)
        parts.append(json.dumps(blob))
        # Also the repr of every float, so "1.676e-15" matches a stored
        # 1.6757...e-15 only through the rounded forms below, never by luck.
    for rel in SOURCES:
        path = os.path.join(ROOT, rel)
        if os.path.isfile(path):
            with open(path, encoding="utf-8") as fh:
                parts.append(fh.read())
    return "\n".join(parts)


def rounded_forms(hay: str) -> set:
    """Every rounding of every float in the data, so a record may quote
    ``1.676e-15`` for a stored ``1.6757...e-15`` and still be checkable."""
    out = set()
    for token in re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?", hay):
        try:
            value = float(token)
        except ValueError:
            continue
        for places in range(1, 10):
            out.add(f"{value:.{places}f}")
            out.add(f"{value:.{places}e}")
            out.add(f"{value:.{places}g}")
    return out


def main():
    record = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        ROOT, "results", "exact_memo_2026-10-10.md")
    with open(record, encoding="utf-8") as fh:
        text = fh.read()

    bad = [p for p in PLACEHOLDERS if p in text]
    if bad:
        print(f"PLACEHOLDERS SURVIVE: {bad}", file=sys.stderr)
        return 1

    text_scanned = DATE.sub(" ", text)
    hay = haystack()
    forms = rounded_forms(hay)
    unverified = []
    for token in NUM.findall(text_scanned):
        token = token.rstrip(".").rstrip("x").rstrip("+-")
        if not token or significant(token) < 4:
            continue
        if token in hay or token in forms:
            continue
        unverified.append(token)
    unverified = sorted(set(unverified))
    if unverified:
        print("NOT FOUND IN ANY DATA FILE OR CITED SOURCE:", file=sys.stderr)
        for token in unverified:
            print(f"  {token}", file=sys.stderr)
        return 1
    print(f"OK: every >=4-significant-digit decimal in "
          f"{os.path.basename(record)} is substantiated")
    return 0


if __name__ == "__main__":
    sys.exit(main())
