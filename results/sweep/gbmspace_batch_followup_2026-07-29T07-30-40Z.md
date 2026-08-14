# GBM-Space batch follow-up: 24 -> 30 sections, 10 -> 12 patients (AT4, AT13)

This is a short, single-agent follow-on ingestion pass on top of
`results/gbmspace_adapter.md` (§1-11, the full history and protocol). No
other batch was running concurrently (checked `data_cache/.gbmspace_claims/`
before claiming — its 12 existing markers matched exactly the 12
most-recently-ingested sections from §10/§11, no orphaned or foreign claims).
Same adapter (`tmesim/data/adapters/gbmspace.py`, untouched), same
registration path (`app.server.ingest_registry.register`), same
reconstruction call (`reconstruct_from_visium`). No file under
`tmesim/core`, `tmesim/physics`, `tmesim/sim`, or `tmesim/molecular` was
touched. No git command was run in this pass (per instruction — a human
reviews and commits separately).

## 1. Listing the full atlas: 97 sections confirmed, Python `tarfile` streaming used throughout

Per §11's explicit warning that `curl | tar tzvf -` piped to a non-tty
silently under-enumerates and mis-parses size columns, this pass used
Python's `tarfile` module streaming directly over the URL
(`urllib.request.urlopen` -> `tarfile.open(fileobj=..., mode="r|gz")`) for
both the listing and the extraction, never `curl | tar`.

```
TOTAL_MEMBERS=98   (97 .h5ad sections + 1 anndata/README.md)
```

This matches §11's corrected count exactly (97 sections). Cross-referencing
the full listing against the 24 already-ingested sections and their 10
patients (AT3, AT5, AT6, AT7, AT9, AT10, AT11, AT12, AT14, AT15) showed
**AT4 and AT13 are the only two patients in the entire 97-section atlas not
yet represented in this repo** — every other patient code appearing in the
tarball (AT3/5/6/7/9/10/11/12/14/15) already has at least one section here.
No `curl`-based listing was used at any point in this pass.

## 2. What was claimed and ingested: 6 sections, both new patients, 3 each

Per the task's instruction to maximize patient diversity, this pass claimed
3 sections from AT4 and 3 from AT13 — the smallest available `.h5ad` per
patient, to keep per-section download/materialize time and disk churn low —
rather than spreading the 6 across more patients (none of which would have
added a *new* patient, since AT4/AT13 were the only ones left).

All 6 claim markers were written to `data_cache/.gbmspace_claims/` *before*
any download began, and none conflicted with an existing marker (confirmed
by listing the directory immediately beforehand).

A single streamed pass over the compressed `spatial_data_visium.tar.gz`
(3.9 GB) extracted all 6 named members as they were encountered in the
archive's internal ordering — not 6 separate downloads, and the full
archive was never written to disk. Materialize -> register (Studio) ->
`ingest.validate_record` -> delete-raw-immediately ran for each section the
moment its bytes were extracted, before the stream continued to the next
target member.

| Section | Patient | Spots | Genes | Chemistry | Software version | Raw `.h5ad` |
|---|---|---:|---:|---|---|---:|
| AT13-BRA-3-FO-1_2 | AT13 (new) | 3,821 | 36,601 | Visium V1 Slide | spaceranger-2.0.0 | 153.0 MB |
| AT4-BRA5-FO-2_1   | AT4 (new)  | 1,440 | 36,601 | **Spatial 3' v1** | **spaceranger-1.3.0** | 69.2 MB |
| AT13-BRA-3-FO-4_2 | AT13 (new) | 2,604 | 36,601 | **Spatial 3' v1** | **spaceranger-1.3.0** | 51.9 MB |
| AT4-BRA5-FO-4_1   | AT4 (new)  | 2,944 | 36,601 | **Spatial 3' v1** | **spaceranger-1.3.0** | 116.5 MB |
| AT4-BRA5-FO-1_1   | AT4 (new)  | 2,099 | 36,601 | **Spatial 3' v1** | **spaceranger-1.3.0** | 60.9 MB |
| AT13-BRA-3-FO-1_1 | AT13 (new) | 3,039 | 36,601 | Visium V1 Slide | spaceranger-2.0.0 | 101.7 MB |

Every one of the 6 sections: `ing.validate_record(rec)["ok"] == True`, no
errors, and `ing.ingest(..., strict=True)` (the default) completed without
raising — i.e. each section also passed strict validation internally, not
just the separate explicit check logged above.

### A genuine, previously-undocumented finding: a *third* chemistry/pipeline label in this source

§9.3 found two chemistries in GBM-Space (`"Visium V1 Slide"`
spaceranger-2.0.0, 36,601-gene whole-transcriptome; `"Visium V5 Slide"`
spaceranger-2.1.0, 18,085-gene probe-based). **This pass found a third:**
4 of the 6 new sections (3 of 3 AT4 sections, 1 of 3 AT13 sections) report
`chemistry = "Spatial 3' v1"`, `software_version = "spaceranger-1.3.0"` —
an older Space Ranger version and a different chemistry label than either of
the two already documented. Confirmed directly from
`uns/spatial/*/metadata` via `inspect_section()`, not inferred from gene
count: **the gene panel is still the full 36,601-gene whole-transcriptome
set** for these sections, identical to the "Visium V1 Slide" sections — so
this is a metadata/processing-lineage difference (an earlier Space Ranger
run, on the older "Spatial 3' v1" chemistry name 10x used before rebranding
that assay "Visium"), not a different probe panel. `inspect_section()`
required no code change to handle this — it does not branch on chemistry
string, only on `feature_types == "Gene Expression"` — so this is reported
as a new fact about the source data, not a bug or a gap in the adapter.

## 3. A process bug this pass introduced, caught, and fixed: missing `gbmspace_` prefix

`ingest_registry.register(dest_dir)` was initially called without an
explicit `section_id`, unlike every prior pass. `ingest_registry._make_id`
falls back to slugifying the bare directory basename when `section_id` is
`None`, so the 6 new sections were first registered into the Studio
catalogue as `at13_bra_3_fo_1_1`, `at4_bra5_fo_2_1`, etc — **without** the
`gbmspace_` prefix every one of the other 24 sections in this repo's
catalogue uses. This was caught by cross-checking the resulting
`app/runs/_ingested_sections.json` against the existing naming convention
before writing this report, not left in the store.

Fixed: all 6 were unregistered (`ingest_registry.unregister`, which also
correctly removed their entries from `tmesim.data.real.visium`'s in-memory
registry) and re-registered with an explicit `section_id=f"gbmspace_{slug}"`,
matching the pattern every prior pass used. Confirmed in a **fresh Python
process** (not the process that did the fix) that `ingest_registry.load_store()`
re-wires all 6 corrected ids into `data.real.visium`'s registry and that
`reconstruct_from_visium("gbmspace_<slug>")` resolves each one — see §4.

**A second, unrelated instance of the same bug was found, not fixed.** Four
stray entries already present in `app/runs/_ingested_sections.json` before
this pass began — `at3_bra5_fo_1_2`, `at3_bra5_fo_4_1`, `at5_bra_5_fo_1_1`,
`at7_bra_3_fo_3_2` — are the same failure mode (missing `gbmspace_` prefix),
left over from an earlier session's independent-verification re-ingest step
(§10.1/§11 both mention re-ingesting one section per new patient "myself").
Confirmed each stray id points to the exact same `cache_dir` as its correctly
prefixed `gbmspace_*` counterpart (true duplicates, not distinct sections).
An attempt to clean these up in this pass was **blocked by the permission
classifier** as out-of-scope work on a prior, unrelated session's artifacts;
per that guidance this was not worked around. They are harmless (both ids
for the same directory resolve identically; nothing downstream is broken by
the duplication) but are flagged here honestly rather than silently left
undocumented, and are a ready, low-risk cleanup target for a future pass —
or for the four sections newly registered here (`gbmspace_at13_bra_3_fo_1_1`,
`gbmspace_at13_bra_3_fo_1_2`, `gbmspace_at13_bra_3_fo_4_2`,
`gbmspace_at4_bra5_fo_2_1`, `gbmspace_at4_bra5_fo_4_1`,
`gbmspace_at4_bra5_fo_1_1` — all 6 confirmed *clean*, no stray duplicate).

## 4. Reconstruction + smoke test: all 6 new sections, not just one

The task asked for at least one; all 6 were run since the marginal cost is
small once the pipeline is confirmed working (same judgment call §9.4 and
§11 made). `reconstruct_from_visium(...)` — unmodified — followed by the
same geometric-route state-banding computation
(`experiments.patient_report.resistance_state_profile`, steady-state O2 from
the reconstructed `vascular_density` field, then `update_cell_states`) used
in every prior pass:

| Section | Patient | Grid | n_cells (placed) | mean_vascular | mean_ecm | VIABLE | QUIESCENT | HYPOXIC |
|---|---|---|---:|---:|---:|---:|---:|---:|
| AT4-BRA5-FO-4_1   | AT4  | 57x54 | 35,225 | 0.0526 | 0.0675 | 68.7% | 29.7% | **1.5%** |
| AT13-BRA-3-FO-1_1 | AT13 | 54x54 | 36,421 | 0.0121 | 0.4717 | 100.0% | 0.0% | **0.0%** |
| AT13-BRA-3-FO-1_2 | AT13 | 57x47 | 45,730 | 0.0211 | 0.4367 | 17.4% | 82.6% | **0.0%** |
| AT13-BRA-3-FO-4_2 | AT13 | 57x54 | 31,166 | 0.0865 | 0.3156 | 100.0% | 0.0% | **0.0%** |
| AT4-BRA5-FO-1_1   | AT4  | 49x49 | 25,173 | 0.1516 | 0.1806 | 100.0% | 0.0% | **0.0%** |
| AT4-BRA5-FO-2_1   | AT4  | 44x47 | 17,245 | 0.1326 | 0.1257 | 100.0% | 0.0% | **0.0%** |

**Honest read, consistent with §9.4's standing finding, not a new
contradiction of it:** 5 of these 6 sections show the same 0.0% geometric
hypoxic fraction the built-in Visium/Xenium sections show; only
`AT4-BRA5-FO-4_1` shows a non-zero band (1.5%, at the low end of the
2.8%-80.6% range §9.4 already reported). AT4 and AT13 join AT14, AT11 and
AT6 as patients whose sections mostly report 0.0% under this uncalibrated
geometric route — §9.4's "it generalizes partially, not universally"
conclusion is reinforced, not overturned, by these 6 new data points: 8 of
the repo's now-18 whole-transcriptome-chemistry sections show a non-zero
hypoxic band, 10 do not. This is not a marker-based hypoxia measurement
(the per-cell marker route is off by construction for all spot-level
GBM-Space data, unchanged from every prior pass's caveat), and no
calibration work was done or implied by this table.

Capability report (`ing.format_capability_report`, unmodified) confirmed on
both representative sections checked directly
(`AT4-BRA5-FO-4_1`, `AT13-BRA-3-FO-1_1`): hypoxia panel
(VEGFA/CA9/ADM/NDRG1/SLC2A1, union detection 89-96%) and MGMT (detection
15-16%) both measurable, matching every other GBM-Space section in this
repo; per-cell routes (`state_mode:markers`, `oxygen_mode:markers`,
`sensitivity_mode:determinant`) all correctly disabled with
`evidence: platform`, the same by-construction behaviour every prior pass
observed.

## 5. Studio registration and tests

All 6 sections registered into the Studio catalogue under
`gbmspace_at13_bra_3_fo_1_1`, `gbmspace_at13_bra_3_fo_1_2`,
`gbmspace_at13_bra_3_fo_4_2`, `gbmspace_at4_bra5_fo_2_1`,
`gbmspace_at4_bra5_fo_4_1`, `gbmspace_at4_bra5_fo_1_1` (after the §3 fix).
Confirmed in a fresh process: `ingest_registry.load_store()` re-wires all 6
into `tmesim.data.real.visium`'s registry and `reconstruct_from_visium(...)`
resolves each one with no manual setup — same guarantee §5/§9.5/§10.1/§11
all confirmed for their own new sections.

`tests/test_gbmspace_adapter.py`'s `_MATERIALIZED` list was extended with
one representative section per new patient — `AT4-BRA5-FO-4_1` and
`AT13-BRA-3-FO-1_1` — following the exact precedent of §10.3/§11 (one
section per new patient, not all of them, to keep the parametrized suite's
runtime proportional to patient count rather than section count). No
existing assertion needed loosening (both new sections are
whole-transcriptome, 36,601 genes, 2,944 and 3,039 spots — comfortably
inside the existing `n_cells > 500` / `n_genes > 15000` thresholds set in
§10.3 against the full 780-12,120 spot-count range already on disk).

```
./.venv/bin/python -m pytest tests/test_gbmspace_adapter.py -q
-> 24 passed in 194.89s
```

+2 tests (parametrization grew from 22 to 24 entries), 0 regressions, 0
skips (all required caches present in this checkout). Per the task's
instruction, the full repo-wide `pytest tests/` suite was **not** run in
this pass — only `test_gbmspace_adapter.py`.

## 6. Disk discipline

| Checkpoint | Free space |
|---|---:|
| Before this pass (claims-dir check + `df`) | 19.45 GiB |
| After each of the 6 raw `.h5ad` deletions | 19.28-19.42 GiB (never dipped below 19.28 GiB) |
| After this pass (final `df`) | 19.51 GiB |

The 10 GiB floor was never approached — the closest point was 19.28 GiB,
9.28 GiB of headroom above the floor. Peak transient disk use was one raw
file at a time (largest: 153.0 MB for `AT13-BRA-3-FO-1_2`); each raw
`.h5ad` was deleted immediately after `gbmspace.register()` materialized it,
before the stream moved on to the next target member — same discipline as
every prior pass.

## 7. Repo total: 24 -> 30 sections, 10 -> 12 patients

Confirmed by listing `data_cache/gbmspace_visium/` directly (not by trusting
this pass's own running count):

```
AT10 x5   AT11 x3   AT12 x1   AT13 x3   AT14 x1   AT15 x4
AT3  x2   AT4  x3   AT5  x2   AT6  x1   AT7  x4   AT9  x1
```

**30 sections across 12 patients** — every patient code present anywhere in
the full 97-section atlas is now represented by at least one section in
this repo.

## 8. What was NOT done, honestly

- **67 of the atlas's 97 sections remain untouched** (97 - 30 = 67). All 12
  patients now have at least one representative section, so any further
  ingestion adds sections-per-patient depth, not new-patient diversity —
  a different, and lower-priority per the task's own framing, kind of
  progress than this pass made.
- **The stray duplicate Studio-registry entries flagged in §3
  (`at3_bra5_fo_1_2`, `at3_bra5_fo_4_1`, `at5_bra_5_fo_1_1`,
  `at7_bra_3_fo_3_2`) were not cleaned up** — an attempt was blocked by the
  permission classifier as out of scope for this task, and that boundary was
  respected rather than worked around. They are harmless duplicates, not a
  correctness problem, but are a concrete, low-risk item for whichever pass
  next touches the Studio registry.
- **No calibration work was done.** `OxygenConfig.consumption`'s
  un-anchored status (flagged in every prior pass, ultimately traced to
  `results/real_drug_panel.md`'s standing caveat) is unchanged by this pass
  and applies identically to the §4 table above.
- **The chemistry/version heterogeneity newly found in §2 (`"Spatial 3'
  v1"` / spaceranger-1.3.0) was documented, not investigated further** —
  e.g. whether it correlates with anything about the AT4/AT13 tumors
  themselves, or is purely an artifact of when each specimen happened to be
  processed, was not examined here.
- **H&E images, Ensembl gene ids**: unchanged from every prior pass — still
  not materialized, still a symbol-as-placeholder respectively.
