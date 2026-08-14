# GBM-Space ingestion — Batch B

Follow-on to `results/gbmspace_adapter.md` (12 sections across 5 patients: AT10,
AT11, AT14, AT15, AT6). This batch ran the coordination protocol against
sibling batches A/C/D pulling from the same 97-section full atlas
(`spatial_data_visium.tar.gz`) concurrently, claimed and ingested 4 more
sections, and confirmed reconstruction on one of them.

## Coordination

- `data_cache/gbmspace_visium/` at start: 12 sections (5 patients), matching
  `gbmspace_adapter.md`.
- `data_cache/.gbmspace_claims/` created (didn't yet exist) and was empty at
  first check.
- Full-atlas archive listing via streamed `tar tzvf -` (no full download,
  ~3.9 GB compressed stream read once for the TOC only) enumerated all 97
  members. **Correction to my own first attempt**: an initial listing script
  used `awk '{print $3, ...}'` to extract file size, which is actually the
  `tar tzvf` **uid** column (always `20944`), not the size (`$5`). Caught this
  before selecting anything, discarded that listing, and re-ran with the
  correct field — flagging it here in case another batch made the same
  mistake, since a wrong "smallest first" ordering would defeat the
  coordination heuristic.
- Waited 14s (randomized), re-checked claims: 8 sections already claimed by
  sibling batches (`AT12-BRA-4-FO-4_1`, `AT13-BRA-3-FO-1_2`, `AT4-BRA5-FO-3_1`,
  `AT5-BRA-5-FO-1_2`, `AT5-BRA-5-FO-2_1`, `AT5-BRA-5-FO-3_1`,
  `AT7-BRA-3-FO-1_1`, `AT9-BRA-4-FO-1_1`).
- Computed available = 97 − 12 ingested − 8 claimed = 77. Applied "skip the
  first 4 smallest, take the next 4": smallest-available were
  `AT3-BRA5-FO-1_2`, `AT3-BRA5-FO-4_1`, `AT5-BRA-5-FO-1_1`,
  `AT7-BRA-3-FO-3_2` (skipped); claimed the next 4: `AT7-BRA-3-FO-3_1`,
  `AT11-BRA-6-FO-4_2`, `AT7-BRA-3-FO-1_2`, `AT15-BRA-4-FO-F1-S18`.
- **A genuine tooling bug caught before it caused damage**: the claim step
  first used `for name in $CANDIDATES` with `CANDIDATES` an unquoted
  space-separated string. This shell is **zsh**, and unlike bash, zsh does
  **not** word-split an unquoted parameter expansion by default — the loop ran
  once with `name` bound to the entire 4-name string, creating a single
  bogus claim file literally named
  `AT7-BRA-3-FO-3_1 AT11-BRA-6-FO-4_2 AT7-BRA-3-FO-1_2 AT15-BRA-4-FO-F1-S18.claimed`
  instead of 4 real per-section claims. Caught immediately by re-listing the
  claims directory, deleted the bogus file, and redid the claim with a zsh
  array (`candidates=(...); for name in "${candidates[@]}"`), which claimed
  all 4 correctly with no collisions. **Flagging this for any other batch
  using a similar shell loop pattern for claim files**: an unquoted
  `for x in $VAR` loop in zsh silently fails to iterate as a bash script
  author would expect, and the failure mode here (one wrong marker file, real
  sections never actually reserved) is exactly the kind of thing that could
  cause two batches to both download the same section without a visible
  error.

## Sections claimed and ingested

All 4 extracted from `anndata/<name>.h5ad` in a single streamed `tar xzf -`
pass over the compressed archive (only the 4 named members pulled, not the
full 3.9 GB), then run through the existing adapter (`gbmspace.register()` →
materialize → `visium.register_sample` → `ingest.detect`/`ingest.ingest` →
`ingest.validate_record`), then registered into the Studio catalogue via
`app.server.ingest_registry.register(path, section_id="gbmspace_<slug>")` to
match the established `gbmspace_<slug>` id convention (the first registration
attempt omitted `section_id` and got an un-prefixed id; caught by checking
`app/runs/_ingested_sections.json` against the other batches' entries,
unregistered, and re-registered correctly). Raw `.h5ad` deleted immediately
after each materialization.

| Section | Patient | Raw `.h5ad` | Spots (in-tissue) | Genes | Chemistry | Studio id | `validate_record` |
|---|---|---:|---:|---:|---|---|---|
| `AT7-BRA-3-FO-3_1` | AT7 (new) | 39.4 MB | 2,854 | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | `gbmspace_at7_bra_3_fo_3_1` | ok=True |
| `AT11-BRA-6-FO-4_2` | AT11 | 43.0 MB | 2,704 | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | `gbmspace_at11_bra_6_fo_4_2` | ok=True |
| `AT7-BRA-3-FO-1_2` | AT7 | 43.4 MB | 1,721 | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | `gbmspace_at7_bra_3_fo_1_2` | ok=True |
| `AT15-BRA-4-FO-F1-S18` | AT15 | 43.5 MB | 12,120 | **18,085** | **Visium V5 Slide / spaceranger-2.1.0** (probe-based) | `gbmspace_at15_bra_4_fo_f1_s18` | ok=True |

`AT15-BRA-4-FO-F1-S18` re-confirms the two-chemistry finding already
documented in `gbmspace_adapter.md` §9.3 (probe-based Visium V5, 18,085-gene
panel, distinct from the whole-transcriptome 36,601-gene chemistry the other
3 sections here use) — not a new format wrinkle, just another data point for
it. `capability_report` confirms both hypoxia markers
(`VEGFA`/`CA9`/`ADM`/`NDRG1`/`SLC2A1`, union detection 5.1%–58.7% across the
4) and MGMT (`MGMT`, union detection 0.37%–1.45%) are present and measurable
on all 4, same pattern as every GBM-Space section ingested so far.

Repo now carries **20 real GBM-Space sections** in `data_cache/gbmspace_visium/`
(the pre-existing 12 + this batch's 4 + 4 more already landed from sibling
batches by the time this batch finished: `AT12-BRA-4-FO-4_1`,
`AT5-BRA-5-FO-1_2`, `AT7-BRA-3-FO-1_1`, `AT9-BRA-4-FO-1_1`).

## Reconstruction check

`app.server.ingest_registry.load_store()` (re-wires all persisted
`gbmspace_*` ids into `tmesim.data.real.visium`'s in-process registry — needed
in a fresh process, same mechanism `gbmspace_adapter.md` §5 documents) then
`reconstruct_from_visium("gbmspace_at7_bra_3_fo_3_1")` — unmodified — ran
cleanly:

```
n_cells (placed) = 34,170   n_tumor = 16,547 (48.4%)
mean_ecm = 0.2014   mean_vascular = 0.0410
```

Same geometric-route state-banding computation `gbmspace_adapter.md` §6/§9.4
uses (steady-state O2 field from the reconstructed real ENDOTHELIAL-derived
`vascular_density`, then `update_cell_states`):

| VIABLE | QUIESCENT | HYPOXIC |
|---:|---:|---:|
| 78.55% | 21.45% | **0.0%** |

Then a cheap `simulate_response` (temozolomide, `sensitivity_mode='uniform'`,
10-minute window, `dt=10s`) completed end to end:
`predicted_response (kill_fraction) = 0.0601`, `viable_tumor_fraction = 0.9399`.

0.0% hypoxic is consistent with `gbmspace_adapter.md` §9.4's own finding that
the non-zero hypoxic band is **not universal** — this section's
`mean_vascular` (0.0410) sits in the middle of the range where §9.4 already
noted the trend is "mostly monotonic but not a clean single predictor" (e.g.
`AT6-BRA-3-FO-2_2` at `mean_vascular=0.024` also showed 0.0% against the
general trend). Not a new finding, not calibrated — same standing caveat
(`OxygenConfig.consumption` never calibrated against measured oxygen tension)
as every prior entry.

## Format issues

None in the adapter or data itself. Two issues found and fixed were entirely
on my side, not the codebase's:

1. My first listing script read the wrong `tar tzvf` column for file size
   (uid instead of size) — caught before any section was picked, re-derived
   correctly.
2. My first capability-report check called
   `ing.capability_report(rec)` (passing the whole `SampleRecord`, whose
   default `platform` param is `"xenium"`) instead of reading the
   already-computed `rec.capability` the way `gbmspace.register()`'s
   `__main__` block and `gbmspace_adapter.md` do. This produced a nonsensical
   report (`platform: xenium`, `panel: 910539 genes`) that looked like a real
   bug at first; reproduced it in isolation, traced it to the wrong call
   signature, and redid it correctly (`rec.capability` via
   `ing.format_capability_report(rec.capability)`), which gave the expected,
   sane Visium capability report for all 4 sections. Recorded here so it
   isn't mistaken for a genuine format wrinkle in any later read of this
   report.

Smoke import: `./.venv/bin/python -c "import tmesim.data.adapters.gbmspace"` →
OK. No file under `tmesim/core`, `tmesim/physics`, `tmesim/sim`,
`tmesim/molecular`, or the adapter itself was touched.

## Disk state

- Before this batch started: 21 GiB free.
- Fluctuated down to 17 GiB during peak concurrent sibling-batch downloads
  (observed mid-pass, before this batch's own download).
- Before this batch's own extraction: 19 GiB free.
- After ingesting all 4 sections and deleting all 4 raw `.h5ad` files: 19 GiB
  free (materialized MEX directories are ~5–6 MB each, 22 MB total for all 4 —
  far smaller than the ~170 MB combined raw `.h5ad` size, deleted immediately).
- Never approached the 10 GiB floor at any point; peak transient use was all
  4 raw files at once during the single streamed extraction (~170 MB total),
  well within budget.
- No leftover raw `.h5ad` or temp extraction directory
  (`data_cache/_batchB_tmp`) remains — removed after ingestion.

## Not done / out of scope for this batch

- Only this batch's 4 sections were pulled; dozens of the 97-section atlas
  remain untouched (left for other batches / future passes, per the
  coordination protocol's intent).
- No pytest suite run (smoke import only, per instructions).
- No git add/commit — left for the central collection pass.
- No calibration work — `OxygenConfig.consumption`'s un-anchored status is
  unchanged, same caveat as every prior GBM-Space entry.
