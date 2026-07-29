# GBM-Space adapter: built, ingested, reconstructed — the follow-on to the format-mismatch report

**Relationship to `results/gbmspace_ingestion.md`.** That report found GBM-Space's public
deposit is a per-section AnnData `.h5ad`, not the Space Ranger `outs/` layout
`tmesim/data/ingest.py`'s Visium detector expects, and stopped there deliberately — building
a reader was out of scope for that pass. This report is that next step: an adapter was
written, 3 sections were downloaded and ingested through it, registered into the Studio
section catalogue, and run end to end through reconstruction and a drug simulation. Nothing
in the earlier report is superseded; §1 below adds one correction of scope (which tarball
the earlier report's probe file actually came from) and one new fact (the selected tarball's
full manifest), neither of which changes its conclusions.

---

## 1. Two small corrections/additions to the earlier report

- The earlier report's 65 MB probe file (`anndata/AT11-BRA-6-FO-2_1.h5ad`) was extracted from
  the **full** `spatial_data_visium.tar.gz` (its path prefix is `anndata/`, not
  `anndata_selected/`) — consistent with that report's own account, just noted here because
  this pass's listing of the **selected** tarball (`spatial_data_visium_selected.tar.gz`,
  streamed with `tar tzvf -`, zero bytes retained) turned up a **different, smaller set of
  section names**:

  | Member | Size |
  |---|---|
  | `anndata_selected/AT10-BRA-5-FO-1_1.h5ad` | 252,670,404 |
  | `anndata_selected/AT15-BRA-4-FO-C2-S27.h5ad` | 855,749,642 |
  | `anndata_selected/AT14-BRA-4-FO-2_1.h5ad` | 222,841,017 |
  | `anndata_selected/AT15-BRA-4-FO-C3-S34.h5ad` | 811,663,185 |
  | `anndata_selected/AT15-BRA-4-FO-B2-S25.h5ad` | 648,490,559 |
  | `anndata_selected/AT10-BRA-5-FO-4_2.h5ad` | 165,163,809 |
  | `anndata_selected/AT10-BRA-5-FO-1_2.h5ad` | 289,873,485 |
  | `anndata_selected/AT10-BRA-5-FO-3_1.h5ad` | 161,885,173 |
  | `anndata_selected/AT10-BRA-5-FO-2_1.h5ad` | 184,232,653 |

  **10 sections total** in the "selected" (manuscript-figure) tarball, not previously
  enumerated. The 3 smallest (`AT10-BRA-5-FO-3_1`, `AT10-BRA-5-FO-4_2`, `AT10-BRA-5-FO-2_1`,
  ~155–185 MB each) were the ones downloaded and ingested in this pass — see §3.

- The schema description in the earlier report (§3) is confirmed correct on every point this
  pass re-derived independently, plus two things it did not check that turned out to matter
  for building an adapter (both below, §2).

---

## 2. What was confirmed by direct inspection, beyond what §3 of the earlier report checked

Both were verified against the actual downloaded bytes with `h5py`, not assumed from the
README quoted in the earlier report:

1. **`obsm/spatial` axis order and units.** It is `(x, y)` = `(pxl_col_in_fullres,
   pxl_row_in_fullres)` in FULLRES PIXEL units — the same convention and units Space Ranger's
   own `tissue_positions.csv` uses. This was confirmed, not assumed, by correlating
   `obsm/spatial[:, 0]` against the independently-recorded `obs/array_col` (`corr = -0.99999999`)
   and `obsm/spatial[:, 1]` against `obs/array_row` (`corr = 0.99999998`) for
   `AT10-BRA-5-FO-3_1` — see `tests/test_gbmspace_adapter.py::test_inspect_section_schema_and_markers`
   for the same check, run on every CI pass against the kept raw `.h5ad`.
2. **`X` is a plain `scipy`-compatible on-disk CSR.** `X/data` + `X/indices` + `X/indptr` plus
   an `attrs["shape"]` on the `X` group is exactly what
   `scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)` expects — **no `anndata`
   sparse-matrix decoder was needed.** `anndata` remains uninstalled in `.venv` (still only
   `h5py` 3.14.0); the trade-off flagged as open in the earlier report's §8 is resolved in
   favour of hand-rolled `h5py` decoding, and it was small: the categorical `obs`/`var`
   columns (`obs/sample`, `obs/sample_name`, `var/feature_types`) are a `categories` array +
   an integer `codes` array, about 8 lines to decode (`gbmspace._read_categorical`).
3. **The `"Gene Expression"` feature-type block is genuinely raw integer counts**, and the
   other 3 categories (`"Cell state abundances"`, `"Spatial niche abundances"`,
   `"Histopath annotation overlap"`) are genuinely continuous — checked directly (`min=1.0`,
   `max=873.0`, 0% non-integer values among the Gene Expression nonzeros for
   `AT10-BRA-5-FO-3_1`; the other blocks' sampled values were fractional, e.g. `12.65`,
   `4.17`, `14.23`). The README's instruction to filter to `"Gene Expression"` before treating
   the matrix as counts is not just cautionary — the other columns really are not counts.

---

## 3. The adapter: what was built, and why it looks the way it does

New files (additive only — nothing under `tmesim/core`, `tmesim/physics`, `tmesim/sim`,
`tmesim/molecular`, or `tmesim/reconstruct/registration_3d.py` was touched):

- `tmesim/data/adapters/__init__.py` — new package for format adapters.
- `tmesim/data/adapters/gbmspace.py` — the adapter itself.
- `tests/test_gbmspace_adapter.py` — 11 new tests (all passing; see §7).

### Design decision: materialize a Space Ranger directory, don't hand-build a `SampleRecord`

The task briefing asked that ingested sections "flow through `reconstruct_from_visium`
UNCHANGED." Tracing that call chain (`tmesim/reconstruct/real_visium.py` →
`tmesim/reconstruct/deconvolve.py` → `tmesim/data/real/visium.py`) shows
`reconstruct_from_visium(sample_id)` resolves a **string id** through
`data.real.visium`'s registry (`register_sample` / `_locate` / `load_from_dirs`), which reads
**literal files off disk** (`tissue_positions.csv`, `scalefactors_json.json`,
`matrix.mtx.gz`+`features.tsv.gz`+`barcodes.tsv.gz`). There is no in-memory entry point that
takes a `SampleRecord` or a dict directly — the whole chain, including
`tmesim/data/ingest.py`'s `detect`/`ingest` and `app/server/ingest_registry.py`'s
`register`, is coupled to "a directory with these exact filenames in it."

So rather than writing a second `SampleRecord`-construction path that duplicates
`ingest._ingest_visium`'s logic (and would still need a *new* reconstruction entry point to
consume it, since `reconstruct_from_visium` doesn't take a record), **the adapter reads the
`.h5ad` with `h5py` and writes out a genuine, ordinary Space Ranger `outs/`-shaped directory**
from those arrays (`gbmspace.materialize_spaceranger_dir`). That directory then flows through
**every** downstream module — `ingest.py`, `visium.py`, `deconvolve.py`, `real_visium.py`,
`ingest_registry.py` — with **zero code changes to any of them**. `ingest.detect()` on a
materialized directory returns confidence `1.0`, platform `"visium"`, layout
`"spaceranger_outs"` — genuinely indistinguishable, at that layer, from `data_cache/visium_gbm/GBM1`.
This is the strongest available sense of "flows through unchanged" given how those modules
are actually wired, and it means every future improvement to the Visium reconstruction path
(deconvolution reference, KDE bandwidth, cell-placement density) applies to GBM-Space
sections automatically, with no adapter-side maintenance.

**Trade-off, stated plainly**: one extra disk-materialization pass per section (a few seconds;
measured 31–36 s per section end to end including `ingest.ingest()`'s detection-rate
measurement over the full 36,601-gene matrix). The materialized MEX triple is *smaller* than
the source `.h5ad` — 32–37 MB per section vs. 155–184 MB raw — because the non-"Gene
Expression" `var` columns (cell2location, NMF) are dropped and gzip compresses the sparse
integer matrix well.

**What is NOT carried over**, and why that's fine here: the H&E histology image
(`uns/spatial/*/images/{hires,lowres}`) is not materialized — nothing in
`tmesim/data/real/visium.py`'s reader touches Space Ranger's image files, so there was
nothing downstream that would use it. GBM-Space's `var/_index` carries HGNC gene **symbols
only, no Ensembl ids** anywhere in the source file; the materialized `features.tsv.gz` writes
the symbol into both the id and symbol columns (real Space Ranger output never does this).
This is flagged in the adapter's module docstring and is a genuine, if minor, divergence from
a real Space Ranger bundle — nothing downstream reads `ensembl_ids` for anything load-bearing
(checked: it is carried through to `native["ensembl_ids"]` and not otherwise consumed).

---

## 4. What was ingested

3 sections (the 3 smallest in the selected tarball, per the task's "prove the pipeline before
scaling up" instruction), downloaded by extracting only those named members from the
streamed `.tar.gz` (no full-archive download):

| Section | Spots (in-tissue) | Genes (Gene Expression) | Materialized size |
|---|---:|---:|---:|
| `AT10-BRA-5-FO-3_1` | 3,272 | 36,601 | 33 MB |
| `AT10-BRA-5-FO-4_2` | 2,736 | 36,601 | 32 MB |
| `AT10-BRA-5-FO-2_1` | 3,144 | 36,601 | 37 MB |

All three: `chemistry = "Visium V1 Slide"`, `software_version = "spaceranger-2.0.0"`,
`obs/sample` provenance strings of the form
`spaceranger200_count_<N>_GBM_SPA<id>_GRCh38-2020-A` (Space Ranger 2.0.0, standard
GRCh38-2020-A reference — same processing lineage the earlier report found for its probe
file). Disk: peak usage ~511 MB raw `.h5ad` (all deleted after materialization except one
154 MB copy kept for adapter tests — see §7), 103 MB materialized MEX output for the 3
sections, well inside the ~17–18 GiB free at the start of this pass (16 GiB free at the end).

Every one of the 3 sections was ingested through **`tmesim.data.ingest.ingest()`, the exact
same function `app/server/ingest_registry.py` and the other 7 real samples use** — no
GBM-Space-specific validation code exists; `validate_record` passed on all 3 with `ok: True`.

### Marker-panel confirmation (the actual point of this whole exercise)

`ing.capability_report()` — unmodified — on all 3 sections confirms:

```
CELL TYPES
  TUMOR        EGFR, SOX2, OLIG2, MKI67  detected in 85–96% of cells
  IMMUNE       PTPRC, CD68, AIF1, CD163  detected in 85–91% of cells
  ENDOTHELIAL  PECAM1  detected in 21–25% of cells
  STROMAL      MOG, AQP4, GJA1  detected in 92–99% of cells

MARKER AXES
  hypoxia   VEGFA, CA9, ADM, NDRG1, SLC2A1   union detection 91–95%
  MGMT      MGMT   union detection 5.9–11.0%
  proliferation  MKI67, TOP2A, PCNA   union detection 21–47%

VERDICTS
  [ OK ] hypoxia_measured   — panel presence + detection rate 91-95%
  [ OK ] mgmt_from_panel    — panel presence + detection 5.9-11.0%
  [ NO ] state_mode:markers        — spot-level, by construction
  [ NO ] oxygen_mode:markers       — spot-level, by construction
  [ NO ] sensitivity_mode:determinant — spot-level, by construction
```

This is the first real GBM section in this repo where **hypoxia and MGMT are both
measurable simultaneously** — the earlier finding (§4 of `gbmspace_ingestion.md`) that "only
`gse324268_xenium` has MGMT, and only `xenium_gbm_io` has any hypoxia marker, and neither has
both" is now closed for this source. Per-cell routes remain unavailable **by construction**
(spot-level data), reusing the exact same disable logic the 4 built-in Visium sections already
hit — nothing new was written for that; the capability report's platform-level checks fire
identically. A richer panel does not unlock per-cell verdicts; it unlocks *spot-level*
hypoxia/MGMT signal that no other real section in this repo carries at all.

---

## 5. Studio registration

All 3 sections were registered into the section catalogue via
`app.server.ingest_registry.register()` — unmodified, same function `/api/sections/scan`
calls:

```
gbmspace_at10_bra_5_fo_3_1   platform=visium  runnable=True  n_cells=3272  n_genes=36601
gbmspace_at10_bra_5_fo_4_2   platform=visium  runnable=True  n_cells=2736  n_genes=36601
gbmspace_at10_bra_5_fo_2_1   platform=visium  runnable=True  n_cells=3144  n_genes=36601
```

Persisted to `app/runs/_ingested_sections.json` (previously `{}`, empty — nothing shadowed).
Confirmed in a **fresh Python process** (not just the process that ran `register()`) that
`ingest_registry.load_store()` re-wires each section into `data.real.visium`'s registry
(`_wire_visium`) and `reconstruct_from_visium("gbmspace_at10_bra_5_fo_3_1")` succeeds with no
manual setup — exactly the guarantee `_wire_visium`'s docstring promises for any ingested
Visium section. The Studio API path (`POST /api/sections/scan` → `register: true`) was
exercised end to end in `tests/test_gbmspace_adapter.py::test_register_materialized_gbmspace_section_via_studio_api`,
mirroring `tests/test_ingest.py::test_register_ingested_visium_section`.

---

## 6. Reconstruction + a cheap drug run, end to end

`reconstruct_from_visium("<registered sample>")` — unmodified — was run on
`AT10-BRA-5-FO-3_1`:

```
grid = (51, 54), dx_um = 100.0
n_cells (placed, 12/spot) = 39,248   n_tumor = 20,146
mean_ecm = 0.0927   mean_vascular = 0.0204
```

then a cheap `simulate_response` (temozolomide, `sensitivity_mode='uniform'`, 1-hour window,
`dt=10s`) completed in 0.1 s: `predicted_response (kill_fraction) = 0.031`,
`viable_tumor_fraction = 0.969`. Full pipeline — h5ad → materialize → detect → ingest →
register → deconvolve → reconstruct → transport → kill — runs end to end on real GBM-Space
bytes.

### Cell-state fractions: geometric route, not marker route — read this before the table

The per-cell `oxygen_mode='markers'` route is **unavailable** on this source (§4) — spot-level
data has no per-cell transcripts to read CA9/VEGFA/HIF1A from directly, regardless of panel
content. The VIABLE/QUIESCENT/HYPOXIC fractions below come from the **same geometric route**
`results/real_drug_panel.md` uses for its Xenium flagship number: solve the steady-state O₂
field from the reconstructed `vascular_density` field (itself built from the real Darmanis-deconvolved
ENDOTHELIAL fraction — real composition, not a marker readout of oxygen), then band every
living cell by local O₂ (`tmesim.cells.agents.update_cell_states`). Reported exactly as
`resistance_state_profile()` in `experiments/patient_report.py` computes it.

| Section | mean_vascular | VIABLE | QUIESCENT | HYPOXIC |
|---|---:|---:|---:|---:|
| `visium_10x_parent_gbm` (built-in) | 0.132 | 94.6% | 5.4% | **0.0%** |
| `GBM1` (built-in) | 0.110 | 71.1% | 28.9% | **0.0%** |
| `AT10-BRA-5-FO-3_1` (GBM-Space) | 0.020 | 11.1% | 64.7% | **24.2%** |
| `AT10-BRA-5-FO-4_2` (GBM-Space) | 0.015 | 4.5% | 34.5% | **61.1%** |
| `AT10-BRA-5-FO-2_1` (GBM-Space) | 0.045 | 36.1% | 58.8% | **5.1%** |

The flagship Xenium section's geometric hypoxic fraction is separately reported at **0.0%**
in `results/real_drug_panel.md`, which calls that "not credible" for a tumor type "defined
partly by pseudopalisading necrosis." Both built-in Visium sections checked here reproduce
that same 0.0% under the identical computation. **All 3 GBM-Space sections instead show
substantial, nonzero hypoxic bands (5–61%)** — qualitatively far more plausible for GBM
biology than 0.0%, and a genuinely different empirical result from every other real section
in this repo.

**Be precise about what this is and is not a claim of.** This is *not* "GBM-Space measures
hypoxia via its markers and the flagship doesn't" — the marker route is off for both, by
construction. What differs is the **reconstructed vascular geometry**: GBM-Space's real
per-spot expression, run through the same Darmanis-based deconvolution every Visium section
here uses, happens to yield a markedly lower mean ENDOTHELIAL fraction (`mean_vascular`
0.015–0.045 vs. 0.110–0.132 for the two built-in sections checked), which starves the same
uncalibrated oxygen-diffusion model of supply and produces the non-zero hypoxic band. Per
`results/real_drug_panel.md`'s own standing caveat, `OxygenConfig.consumption` "has never
been calibrated against a measured oxygen tension" — that caveat applies here unchanged. What
this result *does* establish, honestly: (a) the GBM-Space sections are not degenerate or
pathological inputs — the reconstruction and physics run cleanly and produce a materially
different, more biologically plausible output than the sections already in this repo; (b)
whether that difference is real tumor biology (these sections genuinely are more
hypovascular) or an artifact of Visium spot-pooling / this deconvolution reference / this
uncalibrated oxygen model cannot be determined from this pass alone — it is a question for
calibration work, not something this adapter resolves. This is reported as a qualitative,
directional finding, exactly as the task asked, not as a validated hypoxia measurement.

---

## 7. Tests

`tests/test_gbmspace_adapter.py`, 11 tests, all passing:

- `test_inspect_section_schema_and_markers` — decodes the raw `.h5ad` fixture (see below),
  checks shape, marker presence, integer-valued counts, and re-derives the `obsm/spatial`
  axis-order correlation check from §2 on every run.
- `test_inspect_section_rejects_a_non_gbmspace_h5ad` — a malformed file raises rather than
  silently misreading.
- `test_materialize_then_detect_reaches_full_confidence` — the actual claim of §3: a
  materialized directory hits `ing.detect()` confidence `1.0` with zero new detector code.
- `test_materialize_is_idempotent`.
- `test_register_end_to_end_matches_inspect_section` — parity between `inspect_section()`'s
  raw arrays and what the stock Visium loader reports after the materialize→detect→ingest
  round-trip, mirroring `test_ingest.py`'s parity tests for Xenium/Visium.
- `test_materialized_sections_ingest_and_validate` (parametrized over the 3 shipped sections)
  and `test_capability_confirms_full_marker_panel_on_first_section` — capability-report
  assertions matching `test_ingest.py::test_capability_visium_disables_per_cell_routes_by_construction`'s
  pattern.
- `test_reconstruct_from_visium_and_cheap_kill_run` — the §6 pipeline, asserted mechanically
  (state fractions sum to 1, kill fraction in `[0, 1]`) rather than pinned to exact numbers,
  so it doesn't break on an unrelated physics-tuning change.
- `test_register_materialized_gbmspace_section_via_studio_api` — the §5 Studio path, mirroring
  `test_ingest.py::test_register_ingested_visium_section`.

**One raw `.h5ad` was kept** (`data_cache/gbmspace_visium/_raw_h5ad/AT10-BRA-5-FO-3_1.h5ad`,
154 MB) specifically so `inspect_section()`'s `h5py`-decoding logic is tested against real
bytes rather than a hand-built fixture, consistent with `test_ingest.py`'s "REAL-ONLY, no
invented numbers" stance. All other tests run against the 3 already-materialized directories
under `data_cache/gbmspace_visium/<SAMPLE>/`, mirroring how `test_ingest.py` tests the 7
existing real samples out of `data_cache/`. Every test that depends on a specific cache
`skip`s (not fails) if that cache is absent from a given checkout, matching the existing
`_skip_missing` convention.

**Full suite**:

```
before: ./.venv/bin/python -m pytest tests/ -q  ->  518 passed, 1 warning in 111.17s
after:  ./.venv/bin/python -m pytest tests/ -q  ->  529 passed, 1 warning in 259.83s
```

+11 tests, 0 regressions, 0 skips observed in this checkout (all 3 caches present).

---

## 8. What was NOT done, honestly

- **Only 3 of the 10 sections in the "selected" tarball were ingested**, and neither of the
  other 87 sections in the full `spatial_data_visium.tar.gz` was touched. The task asked for
  "2-3 sections... to prove the pipeline before scaling up" — scaling up to the rest is a
  mechanical repeat of `gbmspace.register(h5ad_path)` per section, not new engineering; it
  was not done here because it was out of the stated scope for this pass, not because
  anything failed.
- **The EGA raw-CRAM deposit (§1(a) of the earlier report) was not pursued.** Still blocked by
  DAC approval, as before.
- **H&E images are not materialized or used anywhere.** Nothing downstream currently reads
  Space Ranger image files, so this was a deliberate scope cut, not an oversight — flagged
  again here in case a future histology-overlay feature needs it.
- **Ensembl gene ids are a placeholder (= symbol)**, because GBM-Space's source file carries
  none. Documented in the adapter's docstring and in every `source` dict it produces; nothing
  currently consumes `ensembl_ids` for anything load-bearing, but this would matter if that
  ever changed.
- **The hypoxic-fraction comparison in §6 is a geometric-route result, explicitly not a
  marker-based hypoxia measurement** — restated here because it is the single easiest part of
  this report to over-claim, and the honest reading is more specific than "GBM-Space measures
  hypoxia."
- **No calibration work was done.** `OxygenConfig.consumption`'s un-anchored status, noted in
  `results/real_drug_panel.md`, is unchanged by this pass and applies identically to every
  number in §6.

---

## 9. Follow-on pass: the rest of the "selected" tarball, plus a taste of the full atlas

This section reports a later pass whose job was exactly what §8 above flagged as not done:
scale past the 3 proof-of-concept sections now that disk space allowed it (24 GiB free at the
start of this pass, checked directly; 20 GiB free at the end). Same adapter
(`tmesim/data/adapters/gbmspace.py`), same registration path (`ingest_registry.register`), same
reconstruction call (`reconstruct_from_visium`) — nothing in §1–8 needed correcting, and no file
under `tmesim/core`, `tmesim/physics`, `tmesim/sim`, or `tmesim/molecular` was touched. The only
tracked-file change is `tests/test_gbmspace_adapter.py` (extended, not rewritten — see §9.5).

### 9.1 One correction to §1: the "selected" tarball has 9 sections, not 10

§1 above said "10 sections total" for `spatial_data_visium_selected.tar.gz` but its own table
lists 9 `.h5ad` members. Re-listing the archive directly in this pass (`curl | tar tzvf -`,
zero bytes retained) confirms **9** is correct — the same 9 rows §1's table already had; "10"
there was simply a miscount (plausibly counting `anndata_selected/README.md` as a section).
This changes nothing about what was already ingested; it just means this pass's 6 new sections
bring the "selected" tarball to **9 of 9 — full coverage**, not 9 of 10.

### 9.2 What was ingested: 6 more sections from "selected", 3 more from the full atlas

**"Selected" tarball, the remaining 6 members** (downloaded via one streamed `tar xzf -` pulling
all 6 named members out of a single pass over the compressed stream — not 6 separate downloads,
and not the full archive written to disk):

| Section | Raw `.h5ad` | Genes (Gene Expression) | Chemistry | Spots |
|---|---:|---:|---|---:|
| `AT10-BRA-5-FO-1_1` | 253 MB | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | 2,956 |
| `AT10-BRA-5-FO-1_2` | 290 MB | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | 3,999 |
| `AT14-BRA-4-FO-2_1` | 223 MB | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | 3,534 |
| `AT15-BRA-4-FO-C2-S27` | 856 MB | **18,085** | **Visium V5 Slide / spaceranger-2.1.0** | 9,697 |
| `AT15-BRA-4-FO-C3-S34` | 812 MB | **18,085** | **Visium V5 Slide / spaceranger-2.1.0** | 11,709 |
| `AT15-BRA-4-FO-B2-S25` | 648 MB | **18,085** | **Visium V5 Slide / spaceranger-2.1.0** | 10,356 |

**Full atlas tarball (`spatial_data_visium.tar.gz`, 97 sections, 3.9 GB), 3 sections** —
pursued per the task's step 3, within a stated, small budget: only the first 7 members of the
archive's internal ordering were in scope (confirmed by streaming `tar tzvf -` first), so
extracting the 3 wanted members required streaming through at most that much of the compressed
file, not the full 3.9 GB. Picked deliberately from **patients not otherwise represented**
(the "selected" tarball's 9 sections span only 3 patients — AT10, AT14, AT15):

| Section | Raw `.h5ad` | Genes | Chemistry | Spots | Patient |
|---|---:|---:|---|---:|---|
| `AT11-BRA-6-FO-1_2` | 114 MB | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | 3,609 | AT11 (new) |
| `AT11-BRA-6-FO-2_1` | 65 MB | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | 3,753 | AT11 (new) |
| `AT6-BRA-3-FO-2_2` | 87 MB | 36,601 | Visium V1 Slide / spaceranger-2.0.0 | 2,297 | AT6 (new) |

**`AT11-BRA-6-FO-2_1` is the exact same file `results/gbmspace_ingestion.md` probed months
earlier** (65,048,148 bytes, 3,753 spots — both numbers match that report's §2/§3 exactly). It
was extracted to a scratch directory then, inspected with `h5py`, and deleted without being
ingested, because at that point no adapter existed yet. This pass is the first time that
specific file has actually gone through the pipeline.

**Disk discipline**: each raw `.h5ad` was deleted immediately after `materialize_spaceranger_dir`
converted it (same pattern as the first 3 sections — the materialized MEX triple is far smaller
than the source file). Peak transient disk use was one raw file at a time (largest: 856 MB for
`AT15-BRA-4-FO-C2-S27`) plus whatever had already been materialized; free space never dropped
below **17.5 GiB** at any point this pass's ingestion scripts logged it. Final state: **9
sections** materialized under `data_cache/gbmspace_visium/` (978 MB total, including the one
154 MB raw `.h5ad` kept for tests since the first pass) **+ 3 more** (total 12 real GBM-Space
sections now in this repo), **20 GiB free** at the end. The full 97-section tarball itself was
never downloaded in bulk — only 3 named members were ever pulled from its compressed stream.

**Stopping point, stated honestly**: this pass stopped at 9 (all of "selected") + 3 (from the
full atlas) not because disk ran out — 20 GiB free remained comfortable headroom — but because
the task's own framing treated the full 97-section atlas as a "few sections, worth evaluating"
add-on once the selected tarball was exhausted, not a second full scale-up pass. 12 real
sections across 5 patients was judged enough to answer the generalization question (§9.4) without
further engineering; pulling more of the full atlas remains a mechanical repeat of the same
`gbmspace.register(h5ad_path)` call per member, exactly as §8 already said about the selected
tarball.

### 9.3 A finding that wasn't anticipated: two different Visium chemistries in this one source

The first pass's report (and its adapter docstring) implicitly assumed one GBM-Space chemistry —
Visium V1 Slide, spaceranger-2.0.0, 36,601-gene whole-transcriptome panel — because all 3
sections it ingested happened to share it. **This pass's 3 `AT15-BRA-4-FO-*` sections do not**:
`inspect_section()` reports `chemistry = "Visium V5 Slide"`, `software_version =
"spaceranger-2.1.0"`, and a 18,085-gene panel (confirmed by direct re-download-and-inspect of
one of the three, `AT15-BRA-4-FO-B2-S25`, specifically to check this rather than guess from the
gene count alone — `obs/sample` provenance string
`spaceranger210_count_47450_GBM_SPA13740723_GRCh38-2020-A`, same `GRCh38-2020-A` reference,
different Space Ranger version and slide chemistry).

18,085 genes is consistent with 10x Genomics' probe-based Visium (CytAssist / FFPE) Human
Transcriptome Probe Set, as distinct from the poly-A-capture whole-transcriptome chemistry the
other 9 sections in this repo's GBM-Space collection use. This adapter handles both without any
code change — `inspect_section()` already filters to `feature_types == "Gene Expression"`
generically, whatever that block's width turns out to be — but it is a real, previously-unknown
heterogeneity in the source data, not a bug. **MGMT and the full hypoxia marker panel
(`VEGFA`, `CA9`, `ADM`, `NDRG1`, `SLC2A1`) are present on both panels**, confirmed for every one
of the 6 new "selected"-tarball sections via the same `capability_report()` path as the first 3
— the smaller probe panel does not drop either axis. `tests/test_gbmspace_adapter.py`'s
`n_genes > 36000` assertion (which would have failed on the 3 `AT15` sections) was loosened to
`n_genes > 15000` with a comment explaining why, rather than silently narrowing the parametrized
section list to dodge the difference.

### 9.4 Reconstruction on all 9 new sections: does the earlier hypoxia finding generalize?

`reconstruct_from_visium` — unmodified — ran cleanly on all 9 new sections (the task asked for
"at least 2"; all 9 were run, since the marginal cost per section is small once the pipeline is
confirmed working). Same geometric-route state-banding computation as §6 above (steady-state O₂
from the reconstructed `vascular_density` field, then `update_cell_states`), reported the same
way:

| Section | Patient | Chemistry | mean_vascular | n_cells (placed) | tumor % | VIABLE | QUIESCENT | HYPOXIC |
|---|---|---|---:|---:|---:|---:|---:|---:|
| `AT10-BRA-5-FO-3_1` *(1st pass)* | AT10 | whole-transcriptome | 0.020 | 39,248 | 51.3% | 11.1% | 64.7% | **24.2%** |
| `AT10-BRA-5-FO-4_2` *(1st pass)* | AT10 | whole-transcriptome | 0.015 | — | — | 4.5% | 34.5% | **61.1%** |
| `AT10-BRA-5-FO-2_1` *(1st pass)* | AT10 | whole-transcriptome | 0.045 | — | — | 36.1% | 58.8% | **5.1%** |
| `AT10-BRA-5-FO-1_1` | AT10 | whole-transcriptome | 0.039 | 35,451 | 64.3% | 47.4% | 49.8% | **2.8%** |
| `AT10-BRA-5-FO-1_2` | AT10 | whole-transcriptome | 0.035 | 47,956 | 61.1% | 17.8% | 82.2% | **0.0%** |
| `AT14-BRA-4-FO-2_1` | AT14 | whole-transcriptome | 0.351 | 42,456 | 41.3% | 100.0% | 0.0% | **0.0%** |
| `AT11-BRA-6-FO-1_2` | AT11 | whole-transcriptome | 0.050 | 43,203 | 41.4% | 80.5% | 19.5% | **0.0%** |
| `AT11-BRA-6-FO-2_1` | AT11 | whole-transcriptome | 0.055 | 44,958 | 9.4% | 100.0% | 0.0% | **0.0%** |
| `AT6-BRA-3-FO-2_2` | AT6 | whole-transcriptome | 0.024 | 27,521 | 30.1% | 79.2% | 20.8% | **0.0%** |
| `AT15-BRA-4-FO-C2-S27` | AT15 | probe-based (V5) | 0.019 | 116,286 | 67.1% | 6.5% | 12.9% | **80.6%** |
| `AT15-BRA-4-FO-C3-S34` | AT15 | probe-based (V5) | 0.046 | 140,435 | 88.1% | 15.0% | 32.9% | **52.1%** |
| `AT15-BRA-4-FO-B2-S25` | AT15 | probe-based (V5) | 0.038 | 124,216 | 57.6% | 36.7% | 22.6% | **40.7%** |

**Honest answer to the question this section was asked: it generalizes partially, not
universally, and it correlates with something that wasn't controlled for.**

- **Non-zero geometric hypoxic fractions recur**: 7 of 12 real GBM-Space sections now show a
  non-zero hypoxic band (2.8%–80.6%), so the original finding — that this source produces
  materially different, non-degenerate hypoxic fractions unlike the 0.0% every built-in
  Visium/Xenium section in this repo reports — is not a fluke of the original 3 sections. The
  single highest value seen anywhere in this repo, **80.6%**, comes from this expanded set, not
  the original 3.
- **It is not universal**: 5 of 12 sections show the same 0.0% the built-in sections show —
  `AT10-BRA-5-FO-1_2` (same patient as three non-zero sections), all of `AT14`, both of `AT11`,
  and `AT6`. Three of the repo's now-5 GBM-Space patients (AT14, AT11, AT6) show **zero** hypoxic
  signal across all their sections (4 sections, 3 patients) — every non-zero-or-higher-than-5%
  section in this table comes from just 2 patients, AT10 and AT15.
- **A confound this pass cannot rule out**: all 3 of the highest hypoxic values (40.7–80.6%) come
  from the 3 sections that also share the different, probe-based Visium chemistry (§9.3) —
  a chemistry no other section in this comparison uses. Whether that reflects genuine tumor
  biology (AT15's tumor really is more hypovascular) or a processing artifact (probe-based
  capture producing a systematically different deconvolved ENDOTHELIAL fraction than poly-A
  capture, which is what `mean_vascular` — and therefore the O₂ field and the hypoxic band — is
  computed from) cannot be distinguished from this pass alone. This is flagged, not resolved;
  resolving it would need either same-chemistry sections from an AT15-like hypovascular patient,
  or a same-tissue Visium V1/V5 comparison, neither of which exists in this dataset.
  Restricting to the 9 whole-transcriptome sections alone (removing that confound) still shows a
  non-zero hypoxic band in 4 of 9 (2.8%–61.1%), so the qualitative finding survives even with the
  probe-based sections set aside entirely — just at lower magnitudes.
- **`mean_vascular` is not a clean single predictor even within one chemistry**: restricted to
  the 9 whole-transcriptome sections, hypoxic fraction falls roughly as `mean_vascular` rises
  (0.015→61.1%, 0.020→24.2%, 0.035→0.0%, 0.039→2.8%, 0.045→5.1%, 0.050→0.0%, 0.055→0.0%) — a
  real, mostly-monotonic trend — but `AT6-BRA-3-FO-2_2` breaks it (`mean_vascular` 0.024, among
  the lowest in the set, yet 0.0% hypoxic, not the ~50%+ the trend would predict). One exception
  in nine points is not disqualifying, but it means "just look at vascular density" is not a
  reliable summary of what this uncalibrated model does section-to-section — consistent with
  `results/real_drug_panel.md`'s standing caveat that `OxygenConfig.consumption` has never been
  calibrated against measured oxygen tension, restated here because this expanded set makes the
  model's sensitivity to something other than mean vascular density more visible, not less.

**What this is still not**: not a marker-based hypoxia measurement (§6's caveat applies
unchanged — GBM-Space's per-cell marker route is off by construction for all 12 sections, spot-
level data), not a calibrated result, and not evidence either way on the chemistry confound. It
is a wider, still-qualitative replication: real GBM-Space tissue, across 5 real patients now
instead of 1, keeps producing non-zero and sometimes very large geometric hypoxic fractions that
every other real section in this repo (Xenium flagship, both built-in Visium sections) reports as
exactly 0.0% — while also, in the majority of new sections, reproducing that same 0.0%. Both
halves of that sentence are the honest finding; reporting only the high numbers would overclaim.

### 9.5 Studio registration, tests, and the full suite

All 9 new sections registered via `app.server.ingest_registry.register()` — unmodified — under
ids `gbmspace_at10_bra_5_fo_1_1` … `gbmspace_at6_bra_3_fo_2_2` (same `gbmspace_<slug>` pattern as
the first 3). Confirmed in a fresh process, same check as §5: `ingest_registry.load_store()`
re-wires all **12** `gbmspace_*` ids into `tmesim.data.real.visium`'s registry, and
`reconstruct_from_visium(...)` resolves each one with no manual setup. The Studio catalogue now
carries 12 GBM-Space sections alongside the 7 built-in ones.

`tests/test_gbmspace_adapter.py` was extended, not rewritten: `_MATERIALIZED` now lists all 9
sections that survive a fresh checkout (the original 3 plus this pass's 6 "selected"-tarball
ones; the 3 full-atlas ones are not added to this list since, per the existing convention, only
sections actually shipped/reproducible from a clean `data_cache/` checkout belong there — same
reasoning as the single kept raw `.h5ad`). The one assertion that would have broken on the new
probe-based sections (`rec.n_genes > 36000`) was loosened to `> 15000` with a comment explaining
the two-chemistry finding, rather than quietly dropped or narrowed to dodge it.

```
before this pass: ./.venv/bin/python -m pytest tests/ -q  ->  577 passed, 1 warning in 278.55s
after this pass:  ./.venv/bin/python -m pytest tests/ -q  ->  583 passed, 1 warning in 275.60s
```

+6 tests (the `_MATERIALIZED` parametrization grew from 3 to 9 entries — 3 tests × 6 new
sections), 0 regressions, 0 skips (all caches present in this checkout). `git status` before and
after confirms the only tracked-file change is `tests/test_gbmspace_adapter.py`; no file under
`tmesim/core`, `tmesim/physics`, `tmesim/sim`, `tmesim/molecular`, `app/web/`, `marketing/`,
`landing-3d/`, `tmesim/ml/`, or `tmesim/data/adapters/codex.py` was touched.

### 9.6 What was NOT done, honestly (this pass)

- **Only 3 of the full atlas's ~97 sections were pulled**, deliberately, per the task's framing
  of that tarball as a "worth evaluating, within a clear budget" add-on rather than a second
  full scale-up. 94 sections remain untouched; pulling more is the same mechanical
  `gbmspace.register(h5ad_path)` call §8 already described, now proven to work across two
  different Visium chemistries as well.
- **The chemistry confound in §9.4 was flagged, not resolved.** Distinguishing "AT15 patient
  biology" from "probe-based capture artifact" needs data this pass does not have (either a
  same-chemistry hypovascular patient, or the same tissue run on both chemistries) and was not
  pursued further here.
- **No calibration work was done**, same caveat as §8, now applying to 12 sections instead of 3.
- **H&E images, Ensembl ids**: unchanged from §8 — still not materialized, still a
  symbol-as-placeholder respectively.
