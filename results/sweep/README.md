
## xplatform_deviations_amd64_run31963874951.jsonl

The measurement behind the cross-platform tolerance in
`tests/test_time_resolved_pk.py`. 59 field-level comparisons of a linux/amd64
GitHub runner against the macOS/arm64 stored references, from CI run
31963874951 — the first run that reached every field instead of dying at the
first mismatch.

Kept here because the numbers quoted in that module's comment (worst
2.1488e-14, 89.3x headroom, 33 of 59 fields bit-identical, mean_hist rising
monotonically with checkpoint index) are derived from this file, and the CI
artifact it came from expires 2026-11-14. A number in a comment whose source
has been garbage-collected is a number nobody can check.
