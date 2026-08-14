# GBM-Space ingestion — Batch A

Follow-on to `results/gbmspace_adapter.md` (12 sections across 5 patients: AT10,
AT11, AT14, AT15, AT6). This batch ran the coordination protocol against
sibling batches (B/C/D, running concurrently on the same 97-section full
atlas, `spatial_data_visium.tar.gz`), claimed 4 sections, ingested them, and
confirmed reconstruction on one.

## Coordination

- `data_cache/gbmspace_visium/` at start: 12 sections (5 patients), matching
  `gbmspace_adapter.md`. `data_cache/.gbmspace_claims/` did not exist yet;
  created it, empty at first check.
- **First listing attempt used `curl | tar tzvf -` and was unreliable**: the
  output only enumerated 65 of the archive's 97 members, and the first 36
  entries were missing their size/permission columns entirely (just a bare
  `20944 anndata/<name>.h5ad` — `20944` is the uid field, not a size), likely a
  `bsdtar` listing-format quirk under a non-tty redirect on this stream.
  Discarded that listing rather than risk picking "smallest" from bad data,
  and re-listed using Python's `tarfile` module streaming over the same URL
  (`tarfile.open(fileobj=urllib.request.urlopen(url), mode="r|gz")`), which
  parses tar headers directly rather than relying on a CLI's display
  formatting. That run correctly enumerated all 97 `.h5ad` members with
  reliable integer sizes.
- Re-checked `data_cache/gbmspace_visium/` and `.gbmspace_claims/` immediately
  before claiming (time had passed during the two listing attempts): sibling
  batches had already ingested 8 more sections and placed 8 claim markers by
  that point (`AT11-BRA-6-FO-4_2`, `AT12-BRA-4-FO-4_1`, `AT15-BRA-4-FO-F1-S18`,
  `AT5-BRA-5-FO-1_2`, `AT7-BRA-3-FO-1_1`, `AT7-BRA-3-FO-1_2`,
  `AT7-BRA-3-FO-3_1`, `AT9-BRA-4-FO-1_1`).
- Walked the size-sorted candidate list, skipping anything already present or
  claimed, and claimed via an atomic exclusive create (`set -C` / noclobber
  redirect — fails if another process already created the same marker,
  avoiding a benign-touch race). Claimed, in order: `AT3-BRA5-FO-1_2`
  (23.7 MB), `AT3-BRA5-FO-4_1` (28.2 MB), `AT5-BRA-5-FO-1_1` (29.6 MB),
  `AT7-BRA-3-FO-3_2` (37.2 MB) — the 4 smallest unclaimed, unheld members at
  that moment. No collisions encountered during claiming.
- (Confirmed after the fact via Batch B's own report: it independently
  reached the same 4 names as its "smallest available" candidates and
  correctly skipped all 4 because my claim markers were already there — the
  claim mechanism worked as intended.)

## Sections claimed and ingested

All 4 extracted from `anndata/<name>.h5ad` in a single streamed `tar xzf -`
pass over the compressed archive (only the 4 named members pulled in one pass,
not the full 3.9 GB written to disk), then run through the existing adapter
exactly as documented in `gbmspace_adapter.md` §3: `gbmspace.register()` →
`inspect_section` (h5py) → `materialize_spaceranger_dir` →
`visium.register_sample` → `ingest.detect`/`ingest.ingest` →
`ingest.validate_record`. Raw `.h5ad` deleted immediately after each
materialization. Then registered into the Studio catalogue via
`app.server.ingest_registry.register(path, section_id="gbmspace_<slug>")`,
matching the established `gbmspace_<slug>` id convention (verified against
`app/runs/_ingested_sections.json` after each call, with up to 3 retries per
section, to guard against the registry's read-modify-write persistence race
under concurrent sibling-batch writers — all 4 persisted on the first
attempt).

| Section | Patient | Raw `.h5ad` | Spots (in-tissue) | Genes | Chemistry | Studio id | `validate_record` |
|---|---|---:|---:|---:|---|---|---|
| `AT3-BRA5-FO-1_2` | AT3 (new) | 23.7 MB | 321 | 36,601 | Visium (whole-transcriptome) | `gbmspace_at3_bra5_fo_1_2` | ok=True |
| `AT3-BRA5-FO-4_1` | AT3 (new) | 28.2 MB | 2,040 | 36,601 | Visium (whole-transcriptome) | `gbmspace_at3_bra5_fo_4_1` | ok=True |
| `AT5-BRA-5-FO-1_1` | AT5 (new) | 29.6 MB | 813 | 36,601 | Visium (whole-transcriptome) | `gbmspace_at5_bra_5_fo_1_1` | ok=True |
| `AT7-BRA-3-FO-3_2` | AT7 (new) | 37.2 MB | 2,374 | 36,601 | Visium (whole-transcriptome) | `gbmspace_at7_bra_3_fo_3_2` | ok=True |

Two new patients this batch adds to the repo's GBM-Space collection: **AT3**
and **AT5** (previously only AT10, AT11, AT14, AT15, AT6 were represented).
`AT3-BRA5-FO-1_2` is notably small (321 in-tissue spots — the smallest
GBM-Space section ingested anywhere in this repo so far), consistent with it
being the smallest raw file in the whole 97-member atlas by a wide margin
(next smallest untouched-by-any-batch member was ~24 MB). `capability_report`
confirms hypoxia markers (`VEGFA`/`CA9`/`ADM`/`NDRG1`/`SLC2A1`, union
detection 32.9%–73.4% across the 4) and MGMT (`MGMT`, union detection
0.20%–1.72%) present and measurable on all 4 — same pattern every prior
GBM-Space section has shown. All 4 use the whole-transcriptome
(spaceranger-2.0.0, 36,601-gene) chemistry — no new probe-based (V5) section
in this batch.

Repo now carries several more GBM-Space sections than the 12 this batch
started from (this batch's 4, plus others landed concurrently by sibling
batches B/C/D by the time this batch finished).

## Reconstruction check

`app.server.ingest_registry.load_store()` (re-wires all persisted `gbmspace_*`
ids into `tmesim.data.real.visium`'s in-process registry, needed in a fresh
process — same mechanism `gbmspace_adapter.md` §5 documents) then
`reconstruct_from_visium("gbmspace_at3_bra5_fo_4_1")` — unmodified — ran
cleanly:

```
grid = (50, 56), dx_um = 100.0
n_cells (placed) = 24,437   n_tumor = 14,807 (60.6%)
mean_ecm = 0.0872   mean_vascular = 0.1358
```

Same geometric-route state-banding computation `gbmspace_adapter.md` §6/§9.4
uses (`experiments.patient_report.resistance_state_profile`: steady-state O2
field from the reconstructed real ENDOTHELIAL-derived `vascular_density`,
then `update_cell_states`):

| VIABLE | QUIESCENT | HYPOXIC |
|---:|---:|---:|
| 100.0% | 0.0% | **0.0%** |

`mean_vascular` here (0.136) is the highest seen for any GBM-Space section in
this repo so far (previous range 0.015–0.351, with `AT14-BRA-4-FO-2_1` at
0.351 the only higher value) — well above the ~0.11–0.13 range the two
built-in Visium sections show, and the 0% hypoxic result is exactly what
`gbmspace_adapter.md` §9.4's own trend (hypoxic fraction falls as
`mean_vascular` rises) would predict at this vascular density. Consistent
with, not a departure from, the standing finding — and the standing caveat
applies unchanged (`OxygenConfig.consumption` never calibrated against
measured oxygen tension; this is not a marker-based hypoxia measurement).

## Format issues

None. No new format wrinkle in the `.h5ad` schema itself — all 4 sections
decoded and materialized with zero adapter changes. The only issue this batch
hit was tooling-side (my own first listing attempt, see Coordination above),
not a data-format issue, and was caught before any section was picked.

Smoke import: `./.venv/bin/python -c "import tmesim.data.adapters.gbmspace"` →
OK. No file under `tmesim/core`, `tmesim/physics`, `tmesim/sim`,
`tmesim/molecular`, or the adapter itself was touched.

## Disk state

- Before this batch started: 21 GiB free (per initial `df -h`).
- Before this batch's own extraction (after sibling batches' concurrent
  activity): 19 GiB free.
- After ingesting all 4 sections, deleting all 4 raw `.h5ad`, and removing the
  scratch extraction directory: 19 GiB free (materialized MEX directories are
  2.3–5.4 MB each, ~15 MB total for all 4 — far smaller than the ~119 MB
  combined raw `.h5ad` size, deleted immediately after each materialize).
- Never approached the 10 GiB floor at any point; peak transient use was all 4
  raw files at once during the single streamed extraction (~119 MB total).
- No leftover raw `.h5ad` or scratch directory remains
  (`data_cache/_scratch_batchA` removed after ingestion).

## Not done / out of scope for this batch

- Only this batch's 4 sections were pulled; most of the 97-section atlas
  remains untouched (left for other batches / future passes, per the
  coordination protocol's intent).
- No pytest suite run (smoke import only, per instructions).
- No git add/commit — left for the central collection pass.
- No calibration work — `OxygenConfig.consumption`'s un-anchored status is
  unchanged, same caveat as every prior GBM-Space entry.
