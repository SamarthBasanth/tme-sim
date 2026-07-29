# GBM-Space ingestion attempt: found, inspected, blocked by format — not ingested

**Purpose.** `results/landscape_cancer_census.md` named GBM-Space (97 Visium sections, 12 GBM
tumours) as public spatial data this repo has not touched, and speculated that
`tmesim/data/ingest.py`'s existing Visium path — which reads a stock Space Ranger `outs/` layout
(`spatial/` + `tissue_positions[_list].csv` + `scalefactors_json.json`, plus a
`filtered_feature_bc_matrix/` MEX triple) — "should be usable with zero new ingestion code if
GBM-Space's public deposit uses the standard Space Ranger output shape." This document reports
what was actually found when that assumption was tested.

**Headline finding, stated first because it changes the shape of this whole report: it does
not.** GBM-Space's publicly downloadable spatial data is a per-section **AnnData `.h5ad` object**,
not a Space Ranger `outs/` directory. `tmesim/data/ingest.py`'s Visium detector looks for a
`spatial/` folder and a MEX triple on disk — an `.h5ad` is a single HDF5 file with none of that
directory structure, so `ingest.detect()` would return `platform="unknown"` on it, unconditionally,
regardless of section count. **Nothing was ingested.** No code was written, and no section was
registered. Per the task's explicit instruction — "do not build a whole new ingester in this pass
if the format turns out to need one; that would be new scope, report it instead and stop" — this
report stops at that line and documents precisely what was found instead, including a positive
result on the one question that didn't require ingestion to answer (marker-panel coverage).

---

## 1. What GBM-Space actually is, and where its data lives (two separate deposits)

GBM-Space turns out to be **two separate, differently-gated data deposits**, not one:

**(a) Raw sequencing data — EGA, controlled access.**
[EGA EGAD00001015527](https://ega-archive.org/datasets/EGAD00001015527) ("GBM-Space: Spatial
Transcriptomic Profiling of Glioblastoma (10x Genomics – Visium)") holds 123 files totalling
1.3 TB in **CRAM format** (aligned sequencing reads, not spot-level matrices). Access requires a
formal application and approval from Data Access Committee **EGAC00001000205**, under the
Wellcome Trust Sanger Institute Data Sharing Policy. This is not a small-friction "click to
download" gate — it is a DAC review process with an unknown turnaround time, and even after
approval it hands back raw reads that would need to be re-run through Space Ranger to produce the
directory shape `ingest.py` expects. This deposit was **not** pursued further in this pass; it
does not fit "download a small sample first."

**(b) Processed spatial data — gbmspace.org, genuinely open, no login.**
The project's own portal ([gbmspace.org](https://www.gbmspace.org/), part of the Wellcome LEAP
Delta Tissue programme; see also the
[Sanger Institute collaboration page](https://www.sanger.ac.uk/collaboration/gbm-space-spatial-genomic-atlas-of-glioblastoma/))
serves two direct, unauthenticated downloads, both confirmed live by direct fetch in this pass:

| File | URL | Size |
|---|---|---|
| Full spatial dataset | `https://gbmspace.cog.sanger.ac.uk/spatial_data_visium.tar.gz` | 3.9 GB (per portal; matches a `.tar.gz` whose members alone summed well past 1 GB in the first ~50 of its ~97 sections in this pass) |
| Manuscript-figures subset | `https://gbmspace.cog.sanger.ac.uk/spatial_data_visium_selected.tar.gz` | 770 MB, sections referenced directly in the paper |

Both were confirmed reachable and their **internal layout was inspected without downloading either
file in full** — see §2 for how, and why that matters for the "check disk space, do not commit to
a large download blind" instruction.

**Sources**: EGA dataset page (fetched directly); gbmspace.org homepage (fetched directly, all
hyperlinks enumerated); Sanger Institute collaboration page (surfaced via search);
`README.md` inside `spatial_data_visium_selected.tar.gz` (extracted and read directly, quoted in
§3); the two GBM-Space bioRxiv preprints already cited in `results/landscape_cancer_census.md`
([2025.05.13.653495](https://www.biorxiv.org/content/10.1101/2025.05.13.653495v1),
[2025.05.13.653733](https://www.biorxiv.org/content/10.1101/2025.05.13.653733v1) — bioRxiv blocked
direct fetch of these two in this pass with HTTP 403; their content was not independently
re-verified beyond what the census report and search-surfaced summaries already established).

---

## 2. How the deposit was inspected without a blind large download

Disk was tight (`df -h .` showed **19 GiB free** at the start of this session, falling to **18 GiB**
by the end — consistent with concurrent agents also writing to this machine, not with anything
this task wrote to the repo; see §6). Given that, the two `.tar.gz` files were **never downloaded
to completion or written to the repo**. Instead:

1. `curl -s -L <url> | tar tzvf - | head -N` — streams the gzip-compressed tar directly into
   `tar`'s list mode, discarding bytes after `head` has enough lines. This reveals every file name,
   size, and permission bit inside the archive without writing anything to disk, and without
   downloading the full 3.9 GB (the pipe closes once `head` is satisfied). This was run against
   **both** `spatial_data_visium.tar.gz` and `spatial_data_visium_selected.tar.gz`.
2. Once the layout inside was known (see §3), one small early member
   (`anndata/AT11-BRA-6-FO-2_1.h5ad`, 65,048,148 bytes per the tar listing) was extracted to a
   **scratch directory outside the repo** (`/private/tmp/.../scratchpad/gbmspace_probe/`, not
   `data_cache/`) using a byte-range-capped `curl --range 0-370000000` (370 MB of the compressed
   stream — enough to reach that member) piped through `tar xzOf`. This produced the file at
   exactly its listed size (65,048,148 bytes), confirming a clean extraction.
3. That one `.h5ad` was opened with `h5py` (already installed in `.venv`, version 3.14.0;
   `anndata` itself is **not** installed) and its internal HDF5 group structure was read directly —
   this is the basis for §3 and §4 below.
4. The scratch directory was deleted immediately after inspection
   (`rm -rf .../scratchpad/gbmspace_probe`). **Nothing was written to `data_cache/` or anywhere
   inside the `tme-sim` repo.** `data_cache/` size was checked before and after
   (6.6 GB, unchanged) to confirm this.

This is the "download a small sample first" instruction satisfied in spirit, but the actual small
sample (one `.h5ad`, ~65 MB, peak disk use) was used only to **determine ingestability**, not to
ingest — because the answer to "does the layout match" turned out to be no before any full-section
ingestion was attempted, so there was nothing to feed `tmesim/data/ingest.py`.

---

## 3. The layout, exactly, and exactly why it does not match `ingest.py`'s Visium detector

`tmesim/data/ingest.py`'s Visium path (`resolve_visium_dirs`, `_detect_visium` — read directly in
this pass) requires, as siblings under one parent: a `spatial/` directory containing
`tissue_positions[_list].csv` and `scalefactors_json.json`, plus a `filtered_feature_bc_matrix/` (or
`raw_...`) directory holding a `matrix.mtx[.gz]` + `features.tsv[.gz]` + `barcodes.tsv[.gz]` MEX
triple. That is the Space Ranger `outs/` shape the four existing real Visium samples in
`data_cache/visium_gbm/` all use.

Both GBM-Space tarballs instead contain a **flat directory of one `.h5ad` file per Visium section**
(`anndata/AT11-BRA-6-FO-2_1.h5ad`, `anndata/AT12-BRA-4-FO-2_1.h5ad`, … — confirmed by the streamed
tar listing) plus one `README.md` describing the schema. There is no `spatial/` directory, no
`tissue_positions*.csv`, no `scalefactors_json.json`, and no MEX triple anywhere in either archive.
The README (extracted and read directly — quoted below) confirms this is deliberate, not an
extraction artefact:

> *"This archive contains anndata representations of specific raw spatial transcriptomic data (10X
> Visium)… Otherwise anndata conforms to standard squidpy anndata format… obs: index contains the
> barcode ID and the sample ID… obsm: spatial — Spatial spot coordinates… uns: spatial — H&E image
> and associated metadata."*

Opening the extracted probe file with `h5py` confirmed the concrete internal shape:

```
X/            [csr_matrix]  data(4,181,283) indices(4,181,283) indptr(3,754)
obs/          [dataframe]   _index (barcode_sampleID strings), array_row, array_col,
                             in_tissue, sample (categorical), sample_name (categorical)
obsm/spatial  (3753, 2) float64      <- spot coordinates, but NOT the Space Ranger
                                          pixel/scalefactor convention ingest.py expects
var/          [dataframe]   _index (36,698 feature names), feature_types (categorical)
uns/spatial/  {'AT11-BRA-6-FO-2_1': ...}   <- H&E image + scalefactors, but keyed by
                                              HDF5 group, not a scalefactors_json.json file
```

`obs/sample` for this section reads `spaceranger200_count_GBM_SPA13219366_GRCh38-2020-A` — so the
underlying processing genuinely was Space Ranger 2.0.0 against the standard GRCh38-2020-A
reference. **The raw Space Ranger `outs/` this was built from is not what was redistributed**; what
gbmspace.org serves is the AnnData object built from it. The scale-factor and pixel-coordinate
information `ingest.py`'s Visium loader reads from `scalefactors_json.json` is present here too,
but inside `uns/spatial` as a nested HDF5 group with a different key structure, not as the
standalone JSON file the loader opens by name.

**This is why the census report's "should be usable with zero new code if the layout matches"
premise does not hold here**: the layout does not match, and the mismatch is not a filename
alias `resolve_visium_dirs` could absorb (the existing resolver already handles renamed
`spatial/`-adjacent siblings, tarball-not-yet-extracted cases, and nested `outs/` — see
`tmesim/data/ingest.py` docstring) — it is a **different file format entirely** (one HDF5 file
combining what Space Ranger splits across a directory tree, with `anndata`'s own internal
sparse-matrix and categorical-column encoding). Reading it correctly would require either the
`anndata` package (not installed) or hand-rolled HDF5 parsing of `anndata`'s on-disk schema
(CSR component groups, categorical `categories`/`codes` pairs, the `obsm`/`uns` HDF5-group
conventions) — a genuinely new reader, not a small adaptation of the MEX escape hatch already in
`ingest.py`. Per the task's constraint, that is new-ingester scope and was not built in this pass.

---

## 4. The one question answerable without ingesting: does this panel beat the incumbent flagship?

Even though nothing was ingested through the pipeline, the extracted probe `.h5ad` directly answers
the marker-coverage question the task asked about, using the same gene list the panel actually
ships (`var/_index`, read via `h5py`, no `ingest.py` code involved):

- **Panel breadth**: 36,601 `var` rows are tagged `feature_type == "Gene Expression"` — this is a
  **whole-transcriptome** Visium panel (36,601 matches the standard 10x GRCh38-2020-A reference
  gene count used by Space Ranger 2.0.0, consistent with the `obs/sample` provenance string above),
  not a targeted panel. The remaining 97 `var` rows are non-expression derived features the
  GBM-Space team added — 85 "Cell state abundances" (cell2location deconvolution results) and 12
  "Spatial niche abundances" (NMF-derived) — which the README explicitly warns must be filtered
  out before treating the matrix as raw counts.
- **Every marker this project's incumbent flagship (the 319-gene and 480-gene Xenium panels) is
  missing is present here**: `MGMT` ✅, `VEGFA` ✅, `CA9` ✅, `HIF1A` ✅ — all four are on this
  panel, confirmed by direct membership check against `var/_index`, not inferred. `PECAM1`
  (vasculature), `COL1A1`/`COL1A2`/`COL4A1` (ECM/collagen), and the existing typing markers
  (`MKI67`, `OLIG2`, `SOX2`, `EGFR`, `GFAP`, `PTPRC`, `CD68`, `AIF1`) are all present too — every
  marker checked was found. This matches `results/landscape_spatial_platforms.md`'s general
  finding that whole-transcriptome platforms capture hypoxia/efflux genes "for free," now
  confirmed gene-by-gene for GBM-Space specifically rather than assumed by category.
- **Caveat that matters**: this is spot-level Visium (3,753 spots for this one section, all
  `in_tissue`), exactly like the four Visium sections already in this repo. Per `ingest.py`'s own
  `capability_report` logic (confirmed by reading it directly), Visium **disables the per-cell
  routes by construction** (`state_mode='markers'`, `oxygen_mode='markers'`,
  `sensitivity_mode='determinant'` all report `available: False, evidence: 'platform'`) regardless
  of panel content — a spot pools many cells, so even a panel this complete would not unlock the
  per-cell verdicts the existing 319-gene Xenium panel fails. What it *would* unlock, if ingested,
  is a real MGMT and hypoxia signal at spot resolution — something none of the seven existing real
  sections currently have simultaneously (only `gse324268_xenium` has MGMT, and only
  `xenium_gbm_io` has any hypoxia marker, and neither has both).

---

## 5. What was NOT done, and why (explicitly, per the task's stop condition)

- **No sections were downloaded to `data_cache/`.** Only tar-listing peeks (zero bytes retained)
  and one 65 MB probe file to a scratch directory outside the repo (deleted after inspection).
- **No new ingestion or parsing code was written.** `tmesim/data/ingest.py` was read but not
  modified. Building an `.h5ad` reader — handling `anndata`'s CSR-in-HDF5-groups encoding,
  categorical `obs` columns, filtering `var/feature_types` down to `"Gene Expression"` before
  treating the matrix as counts, and pulling spot coordinates + scale factors out of `obsm/spatial`
  and `uns/spatial` instead of `scalefactors_json.json` — is a genuinely new reader path, not the
  "small, well-understood adaptation" the task authorized. This is exactly the "report it and stop"
  branch of the task's own instructions.
- **No section was registered** via `app/server/ingest_registry.py` — there is nothing ingested to
  register. (For the record: that registration path was read and confirmed to work exactly as
  described in the task — `POST /api/sections/scan {"path": ..., "register": true}` ingests,
  validates, and adds a section with the same shape as the seven built-in ones, refusing to shadow
  a built-in id; see `app/server/ingest_registry.py::register`. It is ready to use the moment a
  real GBM-Space section can be ingested.)
- **Disk space used**: peak ~435 MB in `/private/tmp/.../scratchpad/` (370 MB bounded tar prefix +
  65 MB extracted `.h5ad`), all deleted before this report was written. `data_cache/` is unchanged
  at 6.6 GB. Free space went from 19 GiB to 18 GiB over the session, consistent with other
  concurrent activity on this machine (per `git status`, unrelated files under `app/web/` and
  three `results/*.md` files this task was explicitly told not to touch were modified by a
  concurrent agent during this session), not with anything this task wrote.

---

## 6. Test suite

Run before any investigation, per the task's instruction:

```
./.venv/bin/python -m pytest tests/ -q
518 passed, 1 warning in 156.91s
```

No code in the repo was changed by this task (only this report file was added), so a second run
would be expected to reproduce the same result; it was not re-run given zero production code was
touched and the first run already establishes the pre-existing baseline this task did not disturb.

---

## 7. Honest bottom line

The census report's characterization of GBM-Space as "public data we have not touched" is correct
as far as it goes, but its implicit assumption — that the deposit would hand back Space Ranger
`outs/` directories the existing ingester could swallow whole — does not survive contact with the
actual files. What gbmspace.org actually serves is **openly downloadable (no DAC, no login) but in
a processed AnnData format that requires new reader code**, while the format that *would* match
`ingest.py` zero-code (raw Space Ranger output, or the CRAM this was built from) sits behind an EGA
Data Access Committee. Neither path was a "download 2-3 sections and ingest them" afternoon in the
end — one is blocked by format, the other by access control. That is a useful, if deflating,
finding in its own right: it means the actionable next step for GBM-Space specifically is not
"ingest more data" but "decide whether an `.h5ad` reader is worth building," which is a scoped,
well-understood follow-on now that the exact schema (§3) is documented, not a fresh investigation.

## 8. Next steps, ranked

1. **Write a dedicated `.h5ad` → `SampleRecord` adapter** (new file, e.g.
   `tmesim/data/real/gbmspace_visium.py`, additive — does not touch `ingest.py`'s existing
   detectors) that: opens the file with `h5py` directly (no new dependency — `h5py` is already
   installed; `anndata` is not and does not need to be added just to read this one schema), filters
   `var/feature_types == "Gene Expression"` to 36,601 columns, reads `obsm/spatial` for coordinates,
   and reads `uns/spatial/<sample>/scalefactors` for the pixel-to-micron conversion in place of
   `scalefactors_json.json`. This is the concrete, scoped version of "build a new ingester" that
   this task explicitly deferred rather than started.
2. **Apply for EGA DAC approval (EGAC00001000205)** in parallel — it is slow, but it is the only
   path to the raw Space Ranger `outs/` shape that needs zero new code, and unlocks re-processing
   flexibility (e.g. re-running Space Ranger with a custom reference) that the AnnData redistribution
   does not offer.
3. Once either path lands, this task's §4 marker-coverage finding (MGMT + full hypoxia panel,
   present and confirmed) means GBM-Space is worth prioritizing over other census candidates
   specifically **because** it is the first real GBM section in this repo that would carry both
   axes `results/cell_state_impact.md` found the flagship panel missing, simultaneously.
