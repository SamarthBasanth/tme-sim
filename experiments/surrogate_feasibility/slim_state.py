"""Shrink the ingested-section registry copy to the fields a run needs.

``_state/runs/_ingested_sections.json`` is a copy of the Studio's registry so
the 31 ingested sections resolve from this worktree. The original is ~15 MB
because each record carries a full per-gene capability report, which only the
UI reads. ``cache.lookup_section`` and ``ingest_registry._wire_visium`` need
the identity and path fields only.

This drops the rest and verifies afterwards that all 38 sections still
enumerate. Run once; it is idempotent.
"""
from __future__ import annotations

import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_STORE = os.path.join(_HERE, "_state", "runs", "_ingested_sections.json")

#: Everything a reconstruction and a section record need.
KEEP = ("id", "label", "source", "sample", "cache_dir", "panel", "assay",
        "accession", "has_outcome_label", "platform", "diagnosis")


def main() -> int:
    if not os.path.exists(_STORE):
        print(f"no store at {_STORE}")
        return 1
    before = os.path.getsize(_STORE)
    with open(_STORE) as fh:
        raw = json.load(fh)
    slim = {k: {f: v[f] for f in KEEP if f in v}
            for k, v in raw.items() if isinstance(v, dict)}
    tmp = _STORE + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(slim, fh, indent=1)
    os.replace(tmp, _STORE)
    after = os.path.getsize(_STORE)
    print(f"{len(slim)} records, {before} -> {after} bytes")

    sys.path.insert(0, os.path.abspath(os.path.join(_HERE, "..", "..")))
    from app.server.cache import all_sections

    n = len(all_sections())
    print(f"sections enumerating: {n}")
    if n != 38:
        print("REFUSING to call this a success: expected 38 sections")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
