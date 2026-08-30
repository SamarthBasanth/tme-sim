#!/usr/bin/env python3
"""Reclaim disk from ``app/runs/`` WITHOUT making the server lie about it.

WHY THIS IS A SCRIPT AND NOT ``rm``
-----------------------------------
``app/runs/*/timelapse.bin`` is by far the largest thing this repository
generates -- measured 2026-08-30 at 1,164 MB across 73 runs, against 55 MB of
geometry and 2 MB of results. Deleting those files with ``rm`` reclaims the
disk and costs almost nothing in capability: a time-lapse is re-derivable by
re-running the RunSpec, which every run record still carries in full.

What ``rm`` DOES cost is the truth. ``app/server/api.py::_binary`` answers a
request for a missing payload with a REASON, and for a missing time-lapse that
reason reads "no drug time-lapse was captured for this run (timelapse_frames
was 0, or this run predates the feature)". For a run whose time-lapse we just
deleted, every clause of that sentence is false, and it is delivered with the
same confidence as a true one. This project's most persistent defect is
correct reasoning pointed at the wrong object; a bare ``rm`` here would
manufacture one on purpose.

So this script records what it removed, in the run's own ``result.json``, under
``pruned_payloads``. ``_binary`` reads that first and answers 410 Gone with the
real story: this run DID produce it, it was deleted on this date to reclaim
this many bytes, re-run the spec to get it back.

WHAT IS PROTECTED, AND WHY
--------------------------
1. Any run id named anywhere in ``results/``, ``archive/``, ``docs/`` or
   ``README.md``. Those are dated published records; a reader following a
   citation must find what the citation points at.
2. The ``--keep-recent`` most recent runs, so the Studio still demonstrates the
   time-lapse to anyone opening it, rather than presenting a feature that is
   uniformly 410.

Default is a DRY RUN. Nothing is deleted without ``--apply``.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import re
import sys

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUN_ID_RE = re.compile(rb"r_\d{5}")  # bytes: results/ holds PNGs and
#: pickles alongside the markdown, and a text-mode read would raise on them.

#: Where a citation could live. A run named in any of these is off limits.
CITATION_ROOTS = ("results", "archive", "docs", "README.md", "API_CONTRACT.md")


def cited_run_ids(repo: str = _REPO) -> set:
    """Every ``r_NNNNN`` mentioned in a published record.

    Deliberately a TEXT scan over whole trees rather than a curated list: a
    curated list is a second copy of the truth and would go stale the first
    time someone cited a run in a new file. A false positive here costs some
    disk; a false negative breaks a citation in a dated record.
    """
    found = set()
    for root in CITATION_ROOTS:
        path = os.path.join(repo, root)
        if os.path.isfile(path):
            files = [path]
        elif os.path.isdir(path):
            files = [os.path.join(d, f)
                     for d, _, fs in os.walk(path) for f in fs]
        else:
            continue
        for f in files:
            try:
                with open(f, "rb") as fh:
                    blob = fh.read()
            except OSError:
                continue
            found.update(m.decode() for m in RUN_ID_RE.findall(blob))
    return found


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--runs-dir", default=os.path.join(_REPO, "app", "runs"))
    ap.add_argument("--payload", action="append", default=None,
                    help="payload base name to prune (repeatable); "
                         "default: timelapse")
    ap.add_argument("--keep-recent", type=int, default=10,
                    help="never prune the N most recent runs (default 10)")
    ap.add_argument("--apply", action="store_true",
                    help="actually delete; without it this is a dry run")
    args = ap.parse_args(argv)
    payloads = args.payload or ["timelapse"]

    runs_dir = args.runs_dir
    if not os.path.isdir(runs_dir):
        print(f"no runs directory at {runs_dir}", file=sys.stderr)
        return 1

    names = sorted(n for n in os.listdir(runs_dir)
                   if n.startswith("r_")
                   and os.path.isdir(os.path.join(runs_dir, n)))
    cited = cited_run_ids()
    recent = set(names[-args.keep_recent:]) if args.keep_recent > 0 else set()
    protected = cited | recent

    stamp = _dt.datetime.now().strftime("%Y-%m-%d")
    total = 0
    pruned_n = 0
    skipped_cited = sorted(cited & set(names))

    for name in names:
        if name in protected:
            continue
        d = os.path.join(runs_dir, name)
        doc_path = os.path.join(d, "result.json")
        if not os.path.exists(doc_path):
            continue
        removed = {}
        for payload in payloads:
            f = os.path.join(d, f"{payload}.bin")
            if not os.path.exists(f):
                continue
            size = os.path.getsize(f)
            removed[payload] = {
                "at": stamp,
                "bytes": size,
                "reason": ("Pruned by scripts/prune_run_payloads.py to reclaim "
                           "disk. The RunSpec in this record is complete, so "
                           "re-running it reproduces this payload."),
            }
            total += size
            if args.apply:
                os.remove(f)
        if not removed:
            continue
        pruned_n += 1
        if args.apply:
            # Read-modify-write the run document, atomically, so a crash
            # mid-prune cannot leave a run whose file is gone and whose record
            # does not say so -- which is the exact state this script exists to
            # avoid producing.
            with open(doc_path) as fh:
                doc = json.load(fh)
            doc.setdefault("pruned_payloads", {}).update(removed)
            tmp = doc_path + ".tmp"
            with open(tmp, "w") as fh:
                json.dump(doc, fh, indent=1)
            os.replace(tmp, doc_path)

    mode = "PRUNED" if args.apply else "WOULD PRUNE (dry run; pass --apply)"
    print(f"{mode}: {pruned_n} runs, {total:,} bytes "
          f"({total / 1024 / 1024:.0f} MB) of {', '.join(payloads)}")
    print(f"protected: {len(cited & set(names))} cited by a published record "
          f"({', '.join(skipped_cited) or 'none'}), "
          f"{len(recent)} most recent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
