# GBM-Space ingestion — Batch D

Follow-on to `results/gbmspace_adapter.md` (12 sections across 5 patients: AT10,
AT11, AT14, AT15, AT6). This batch ran the coordination protocol against
sibling batches A/B/C pulling from the same 97-section full atlas
(`spatial_data_visium.tar.gz`) concurrently, then claimed and ingested 4 more
sections — all from **patients not previously represented in this repo**.

## Coordination

- `data_cache/gbmspace_visium/` at start: 12 sections (5 patients), matching
  `gbmspace_adapter.md`.
- `data_cache/.gbmspace_claims/` at start: empty.
- Full-atlas archive listing via streamed `tar tzvf -` (no full download)
  turned up patient prefixes never before ingested: AT3, AT4, AT5, AT7, AT9,
  AT12, AT13 (in addition to the already-done AT10, AT11, AT14, AT15, AT6).
- Waited 35s, re-checked claims (still empty), then claimed 4 sections
  immediately via marker files before downloading: `AT12-BRA-4-FO-4_1.claimed`,
  `AT5-BRA-5-FO-1_2.claimed`, `AT7-BRA-3-FO-1_1.claimed`,
  `AT9-BRA-4-FO-1_1.claimed`. No collisions encountered.
- By the time this batch finished, sibling agents had added their own claims
  for `AT13-BRA-3-FO-1_2`, `AT4-BRA5-FO-3_1`, `AT5-BRA-5-FO-2_1`,
  `AT5-BRA-5-FO-3_1` — confirming the protocol prevented overlap (none of
  those match this batch's picks).

## Sections claimed and ingested

All 4 extracted from `anndata/<name>.h5ad` in a single streamed `tar xzf -`
pass over the compressed archive (not a full 3.9GB download — only the 4
named members were pulled), then run through the existing adapter
(`gbmspace.register()` → materialize → `visium.register_sample` →
`ingest.detect`/`ingest.ingest`), then registered into the Studio catalogue
via `app.server.ingest_registry.register()`, exactly matching the method in
`results/gbmspace_adapter.md` §3–5. Raw `.h5ad` deleted immediately after each
materialization.

| Section | Patient (new) | Raw `.h5ad` | Spots (in-tissue) | Genes | Chemistry | Studio id |
|---|---|---:|---:|---:|---|---|
| `AT12-BRA-4-FO-4_1` | AT12 | 24.4 MB | 3,401 | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | `gbmspace_at12_bra_4_fo_4_1` |
| `AT5-BRA-5-FO-1_2` | AT5 | 33.6 MB | 876 | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | `gbmspace_at5_bra_5_fo_1_2` |
| `AT7-BRA-3-FO-1_1` | AT7 | 47.6 MB | 1,712 | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | `gbmspace_at7_bra_3_fo_1_1` |
| `AT9-BRA-4-FO-1_1` | AT9 | 62.9 MB | 780 | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | `gbmspace_at9_bra_4_fo_1_1` |

All 4 registered `ok=True` with no errors (`ingest_registry.register()`
output), same as the 12 already ingested. `visium.register_sample`'s
in-process `sample_id` for `reconstruct_from_visium` is the raw section name
(e.g. `AT12-BRA-4-FO-4_1`) inside `gbmspace.register()` itself, but only the
Studio slug (`gbmspace_<slug>`) survives across a fresh process via
`ingest_registry.load_store()`'s `_wire_visium` re-wiring — same mechanism
`gbmspace_adapter.md` §5 already documents; not a new format wrinkle, just
noted here since it tripped the reconstruction script on the first attempt in
this pass (fixed by calling `load_store()` before
`reconstruct_from_visium(<studio id>)` rather than the raw section id).

Repo now carries **16 real GBM-Space sections across 9 patients** (AT5, AT6,
AT7, AT9, AT10, AT11, AT12, AT14, AT15).

## Reconstruction check

`reconstruct_from_visium("gbmspace_at12_bra_4_fo_4_1")` — unmodified — ran
cleanly:

```
grid = (45, 54), dx_um = 100.0
n_cells (placed) = 40,756   n_tumor = 18,224 (44.7%)
mean_ecm = 0.4016   mean_vascular = 0.1477
```

`resistance_state_profile()` (`experiments/patient_report.py`, same
geometric-route computation `gbmspace_adapter.md` §6/§9.4 uses — steady-state
O2 field from the reconstructed real ENDOTHELIAL-derived `vascular_density`,
then `update_cell_states` banding):

| VIABLE | QUIESCENT | HYPOXIC |
|---:|---:|---:|
| 100.0% | 0.0% | 0.0% |

Consistent with `gbmspace_adapter.md` §9.4's own finding that the non-zero
hypoxic band is **not universal** across GBM-Space sections — this one's
`mean_vascular` (0.1477) is the highest of any GBM-Space section reconstructed
so far in this repo (higher than the two built-in Visium sections, 0.110–0.132),
which tracks with the same "hypoxic fraction falls as mean_vascular rises"
pattern already flagged there, now extended to a 13th section (0.0% hypoxic at
the highest vascular density seen yet in this collection). Not a new finding,
not calibrated — same standing caveat as every prior entry.

## Format issues

None. Both chemistries previously seen (Visium V1 whole-transcriptome,
36,601 genes) — all 4 of this batch's sections match the more common
whole-transcriptome chemistry; no probe-based (V5) sections were among this
batch's picks. No adapter code changes needed.

## Disk state

- Before this batch's downloads: 18 GiB free (dropped from 21 GiB seen at the
  very start of this pass, due to concurrent sibling-batch downloads).
- After ingesting all 4 sections and deleting all 4 raw `.h5ad` files: 19 GiB
  free.
- Never approached the 10 GiB floor; peak transient use was one raw file at a
  time (largest 62.9 MB) during the single streamed extraction pass.
- No leftover raw `.h5ad` or temp extraction directory
  (`data_cache/_batchD_tmp`) remains — cleaned up after ingestion.

## Not done / out of scope for this batch

- Only 4 of the ~7 newly-visible patient prefixes were pulled (AT12, AT5, AT7,
  AT9 done; AT3, AT4, AT13 left for a future batch — the archive listing this
  pass captured was itself truncated to the first ~48 of 97 members by the
  streamed `tar tzvf | head -200` pipeline ending early under concurrent
  network load from sibling agents, so the remaining ~49 members' contents,
  and any patients beyond those named above, were not enumerated in this
  pass).
- No pytest suite run (smoke import only, per instructions):
  `./.venv/bin/python -c "import tmesim.data.adapters.gbmspace"` → OK.
- No git add/commit — left for the central collection pass.
