# GBM 3D registration: the real Visium/CODEX pair, ingested and registered

**One-line summary.** Greenwald, Darnell, Hoefflin et al. (Tirosh/Suva labs), *Cell* 2024
(PMID 38653236) — the dataset `results/landscape_gbm_serial_sections.md` identified as the
best available near-adjacent, cross-modality, same-block GBM registration target — was
fetched for real (one tumour, `ZH1007nec`), ingested, and registered with
`tmesim/reconstruct/registration_3d.py`'s existing, unmodified multi-start ICP method. The
Visium layer ingests through `tmesim/data/ingest.py` at full confidence with **zero new
ingestion code** (a directory-restructuring step, not a code change — same pattern as the
breast Xenium tar-format quirk). The CODEX layer needed a genuinely new adapter
(`tmesim/data/adapters/codex.py`), and here the honest finding is sharper than the landscape
survey could tell from the outside: **the OME-TIFF's own embedded channel metadata carries
dye names, not marker names** — so per-cell CODEX **marker identity** is unresolved by this
pass, even though per-cell CODEX **position** is not. Registration on real Visium-spot vs
CODEX-nucleus point clouds genuinely improves tissue-frame overlap (Chamfer −40.2%,
Hausdorff-95 −36.8%, density correlation +10.8% relative), by the same honestly-computed
before/after metrics the breast pass used. The paper's own r=0.73 registration-quality number
is **not directly comparable** to either of these metrics — stated precisely below, not forced.

## What this is, precisely

- **Continues, does not replace,** `results/registration_3d_breast_poc.md` (breast Xenium
  Rep1/Rep2, same-modality, same-platform) — this is the flagship-cancer, cross-platform
  extension `results/landscape_gbm_serial_sections.md` recommended as the concrete next step.
- **New code is additive and narrowly scoped**: one new file
  (`tmesim/data/adapters/codex.py`), one new test file
  (`tests/test_codex_adapter.py`), one dependency added (`imagecodecs`, required to decode
  this real file's LZW compression — `pyproject.toml`), and one small, backward-compatible
  parameter added to `register_serial_sections` (`z_spacing_caveat`, see "Registration"
  below — the algorithm itself was not touched). `tmesim/core`, `tmesim/physics`,
  `tmesim/sim`, and `registration_3d.py`'s registration logic are unmodified.
- **Real, live-downloaded data only.** One tumour (`ZH1007`, necrotic region — the smallest
  matched Visium+CODEX pair by combined download size, ~943 MB total), chosen deliberately
  small per the task's disk-budget instruction, out of the paper's 13 IDH-wildtype GBM
  tumours. `df -h` showed 24 GiB free before any download; ~960 MB was used, confirmed after.

## The data

### Source

Zenodo [10.5281/zenodo.12624860](https://zenodo.org/records/12624860) (the current version
of the concept DOI 10.5281/zenodo.8105466 `results/landscape_gbm_serial_sections.md` cites —
resolved by following the redirect, not assumed), "Glioma spatialomics dataset", Greenwald et
al. 2024. 38 files total; full listing fetched via the Zenodo API before downloading anything,
confirming file sizes precisely matched the landscape survey's estimates (Visium per-sample
tarballs 10.2–49.9 MB; CODEX OME-TIFFs 718 MB–4.5 GB; two QuPath project `.zip` archives,
1.8 GB and 2.2 GB, not downloaded — see "What was deliberately not downloaded" below).

### Matching a Visium section to its near-adjacent CODEX section

`visium_dataset_description.csv` and `readme_CODEX.txt` (both <1 KB, downloaded first) confirm
Visium sample names and CODEX filenames share a `<patient>_<region>` naming convention
(`inf`=infiltrating, `T1`=contrast-enhancing, `nec`/`NEC`=necrotic). Three patients
(`ZH1007`, `ZH1019`, `ZH916`) have a clean same-patient, same-region name match across both
file sets; patient `ZH811`/`ZH881`/`ZH8811` does not (three different-looking IDs across the
two file sets for what may or may not be the same patient — **not resolved in this pass,
flagged rather than guessed**). Among the three clean matches, `ZH1007nec` (Visium 22.2 MB +
CODEX 920.9 MB = 943.1 MB) was the smallest combined download, so it was the one fetched.

### Visium layer — ingests through the existing pipeline, zero new ingestion code

The Zenodo Visium tarball is a **genuine Space Ranger 1.1.0 output** (confirmed from
`matrix.mtx`'s own header: `spaceranger-1.1.0`) but laid out **flat**, with every file
prefixed `zh1007nec_` instead of the standard names, no `spatial/` subdirectory, and
`matrix.mtx`/`features.tsv`/`barcodes.tsv` un-gzipped — none of which `resolve_visium_dirs()`
recognises by exact filename. This is the same class of quirk `results/gbmspace_adapter.md`
solved for GBM-Space's `.h5ad`, but much shallower here: the files are ALREADY valid Space
Ranger components, just misnamed/misplaced, so the fix is **symlinks into the expected
`outs/spatial/` and `outs/filtered_feature_bc_matrix/` shape** — a directory-restructuring
step done outside Python (the same "one manual step, not a code change" precedent
`results/registration_3d_breast_poc.md` set for the breast tarball's gzip quirk), not a new
adapter module.

```
detect('data_cache/greenwald_gbm_visium/GBM_ZH1007nec/outs')
  -> platform=visium, layout=spaceranger_outs, confidence=1.0, ingestable=True

ingest(...) -> n_cells (spots) = 955   n_genes = 33,538   extent_um = (2961.0, 3007.5)
  validation.ok = True (1 warning: 24 duplicate gene symbols, expected/handled)
```

This is a **whole-transcriptome** panel (33,538 genes) — far richer than either 10x Xenium
panel this project has used (319/480 genes). The capability report on it, run through the
exact same `ingest.capability_report()` every other sample in this repo goes through:

```
CELL TYPES     TUMOR (EGFR,SOX2,OLIG2,MKI67) 97.07% | IMMUNE 96.65% | ENDOTHELIAL (PECAM1) 36.54% | STROMAL 93.40%
MARKER AXES    proliferation 74.45%  hypoxia 94.03%  MGMT 54.45%  ECM 99.48%  endothelial 36.54%
VERDICTS       [OK] vascular_field  [OK] cell_typing (coverage 1.00)  [OK] hypoxia_measured
               [OK] mgmt_from_panel  [OK] ecm_field
               [NO] state_mode:markers / oxygen_mode:markers / sensitivity_mode:determinant
                    (SPOT-level data — by construction, same as every other Visium sample)
```

**MGMT, hypoxia, and vasculature are all measurable on this real GBM Visium section** —
matching what `results/gbmspace_adapter.md` found for GBM-Space's whole-transcriptome panel,
now confirmed on a second, independent, real GBM whole-transcriptome dataset. Per-cell routes
(`state_mode='markers'`, etc.) are unavailable for the structural reason every Visium sample in
this repo already has: a spot pools several cells, not because of anything specific to this
dataset.

### CODEX layer — ingests through a new, scoped adapter; marker identity is genuinely blocked

`tmesim/data/adapters/codex.py` is new. What it establishes, by direct inspection with
`tifffile` (not assumed from the README):

- `ZH1007_NEC_B.ome.tif`: pyramidal, LZW-compressed, `uint8`, **43 channels**, 7330×8205 px,
  `PhysicalSizeX`/`PhysicalSizeY` = **0.50670 µm/px** (read from the OME-XML), 3 pyramid
  levels. Decoding it required installing `imagecodecs` (tifffile's LZW backend) — added to
  `pyproject.toml`'s core dependencies, since without it `tifffile` raises a plain
  `ValueError` on this real file's actual compression, not at import time.
- **The channel names embedded in the OME-XML are dye names, not marker names.** Channel 0 is
  `"DAPI"`; the other 42 cycle through `"ATTO 550"`, `"FITC"`, `"Cy5"` — the three non-DAPI
  CODEX fluorescence channels, reused every cycle for a different antibody. This is standard
  for a raw-per-cycle CODEX/PhenoCycler export; the (cycle, channel) → marker mapping lives in
  a separate experiment manifest that **this file does not embed**, and that this pass did
  **not** find as a standalone small file in the Zenodo listing (only the two multi-GB QuPath
  project `.zip` archives are large enough to plausibly carry it — QuPath renames channels on
  import — and downloading either was out of this pass's disk budget). **This is the one place
  this task's own anticipated risk ("harder to parse than expected") turned out to be real**,
  and it is reported precisely rather than worked around by guessing a channel order.
- **What IS unambiguous, and is what this adapter uses**: channel 0 is always DAPI (the
  universal nuclear counterstain, imaged every cycle for drift correction) — no manifest
  needed. `extract_nucleus_positions()` Gaussian-smooths the DAPI channel and takes local
  intensity maxima (`skimage.feature.peak_local_max`) as approximate nucleus centres — the
  same class of "cheap, honest stand-in" `registration_3d.tissue_boundary_points` already uses
  for a tissue outline, one level down (per-nucleus, not per-silhouette). **This is explicitly
  not a validated segmentation** — no nucleus boundary, no cell-body assignment, no
  doublet handling; it will both miss real nuclei and produce some noise peaks. It is real
  enough, and distributed correctly enough across the tissue's actual footprint, to drive the
  same density/boundary-based registration method the breast pass used — which is exactly what
  this dataset was fetched to test.
- On the real file, at full resolution (level 0): **40,340 detected nucleus positions**,
  spanning real micron coordinates (79.6–4156.9 µm × 40.5–3573.2 µm), consistent with the
  image's own physical field of view (8205×0.5067 = 4157 µm, 7330×0.5067 = 3714 µm). Took
  34.4 s (dominated by decoding all 43 LZW-compressed channels, ~2.6 GB in memory).
- **One more genuine data-serving quirk, recorded exactly as found**: the OME-XML's own
  internal `TiffData`/`UUID` blocks self-reference the filename `"ZH1017_INF_B.ome.tif"`
  throughout — not `"ZH1007_NEC_B.ome.tif"`, the name Zenodo actually serves this file under
  (patient 1017 vs 1007, region INF vs NEC). Looks like a stale metadata artefact from an
  earlier pipeline/renaming step on the depositors' side. `read_ome_metadata()` detects and
  reports this mismatch rather than silently ignoring it; it does not affect channel 0 = DAPI
  or anything this adapter does with the file.
- **Declared-panel capability check** — the 39-antibody CODEX panel, read directly from the
  paper's own `CODEX_antibody_validations.pdf` (Zenodo, 1.4 MB, fetched and read directly, not
  inferred) and aliased to canonical HGNC gene symbols
  (`codex.CODEX_MARKER_ALIAS_TO_GENE`), run through `ingest.capability_report()` with
  `expression=None` (panel presence only — no per-cell value is claimed):

  ```
  vascular_field:    OK  — CD31 IS on the panel (aliased to PECAM1)
  hypoxia_measured:  OK  — CA9, NDRG1, GLUT1 (aliased to SLC2A1) are all on the panel
  mgmt_from_panel:   NO  — MGMT is genuinely not on this protein panel (honest negative)
  ```

  This directly answers the task's own question ("does this panel carry CD31/PECAM1, hypoxia
  markers?") — **yes, by panel design** — while being explicit that no per-cell measured
  detection RATE exists for any of them in this pass, because that would require the
  marker-identity mapping this pass could not recover.

### What was deliberately not downloaded, and why

Per the task's disk-budget instruction (download the smallest usable subset first): the raw
GEO `GSE237183_RAW.tar` (500.9 MB, all 19 samples bundled) was not fetched — the Zenodo
per-sample tarballs are smaller and already split by tumour. The two QuPath project archives
(`CODEX_qupath_project_v2.zip` 2.2 GB, `CODEX_organoid_qupath_project.zip` 1.8 GB) — which
very likely carry the channel-to-marker manifest this pass could not otherwise recover — were
not fetched; at 22 GiB free before this pass started, either one alone would have used roughly
10% of remaining disk for a single file whose actual content (a Java-serialized QuPath
project) `results/landscape_pathology_3d.md` already found is "likely to be prohibitively
difficult/impossible" to read outside QuPath itself. The other 42 non-DAPI, non-`ZH1007`
OME-TIFF images (718 MB–4.5 GB each) were not fetched — one tumour was the explicit scope for
this pass.

## Registration

`tmesim.reconstruct.registration_3d.register_serial_sections(xy_ref, xy_mov, z_spacing_um, ...)`
was called unmodified in its core algorithm — multi-start (identity, centroid-only, top-5
coarse rotation/scale grid search) each refined by trimmed ICP, final candidate chosen by
post-ICP boundary Chamfer distance, exactly as `results/registration_3d_breast_poc.md`
describes and validated. `xy_ref` = the 955 Visium spot centroids (already real micrometres,
via `scalefactors_json.json`'s `spot_diameter_fullres`); `xy_mov` = the 40,340 CODEX nucleus
positions (already real micrometres, via the OME-XML's `PhysicalSizeX`).

**One small, additive, backward-compatible change was made to `registration_3d.py`**: the
returned `SerialSectionRegistration.z_spacing_caveat` was hardcoded to the module's
breast-specific `Z_SPACING_CAVEAT` constant (about Janesick et al. 2023's Xenium Rep1/Rep2
pair specifically). Reusing the function on a genuinely different dataset (a different paper,
different platforms, different z-spacing provenance) without a way to override that string
would have made the returned object state something false about where its z-spacing
uncertainty comes from. Fix: `z_spacing_caveat` is now an optional parameter defaulting to the
original `Z_SPACING_CAVEAT` constant — every existing caller (including the breast PoC) is
byte-identical in behaviour (`tests/test_registration_3d.py`'s 9 tests still pass unchanged;
`tests/test_codex_adapter.py::test_z_spacing_caveat_defaults_to_module_constant_when_not_passed`
pins this explicitly). This is a documentation-string parameter, not a change to the
registration algorithm.

```python
caveat = (
    "z-spacing between the Visium and CODEX sections of ZH1007nec is UNVERIFIED: "
    "Greenwald et al. 2024 (PMID 38653236) states CODEX was performed on sections "
    "'near-adjacent' to the Visium sections of the same tumour block, but does not "
    "publish the exact micron offset. The value used here is an assumed placeholder, "
    "not a measured quantity."
)
reg = register_serial_sections(
    vis_rec.cell_xy_um, sec.nucleus_xy_um, z_spacing_um=10.0, z_spacing_caveat=caveat,
    ref_sample_id="ZH1007nec_visium", mov_sample_id="ZH1007nec_codex", seed=0,
)
```

### Result

```
transform:  rotation = 0.088°   scale = 0.7314   translation = (366.3, 446.7) um
z_spacing_um = 10.0  (ASSUMED — see caveat above; same class of placeholder the breast
                       pass used, now with dataset-specific prose, not reused verbatim)
n_points = 41,295  (955 Visium spots @ z=0  +  40,340 CODEX nuclei @ z=10.0)
```

| metric | before | after | change |
|---|---|---|---|
| density correlation (60 µm grid) | 0.2457 | 0.2721 | **+10.8%** relative |
| Chamfer distance, tissue outline (µm) | 325.6 | 194.7 | **−40.2%** (closer) |
| Hausdorff-95, tissue outline (µm) | 945.5 | 597.5 | **−36.8%** (closer) |

Winning seed was `seed_index=0` (identity rotation start) — final rotation only 0.088° from
identity, but a **substantial scale correction (0.731)**, i.e. the CODEX nucleus cloud was
shrunk ~27% to best overlap the Visium spot footprint. Full 7-candidate table available in
`reg.candidates`; the runner-up (seed_index=5, a coarse-search seed) converged to a
plausible-looking −115° rotation with WORSE Chamfer (219.4 µm vs the winner's 194.7 µm) — the
same "a wrong candidate can look locally competitive" pattern the breast pass's Chamfer-based
selection exists to catch, working correctly here too.

**On the scale correction, stated precisely rather than left unexplained.** The two platforms'
own pixel calibrations are close (Visium: 55 µm / 112.139 px = 0.4905 µm/px; CODEX: 0.5067
µm/px — a 3.2% difference, nowhere near 27%), and the Visium spot set is already filtered to
`in_tissue==1` only (955 of however many total spots the capture area covers — not background
inflating the apparent footprint). So the 27% scale gap is not a units/calibration bug this
pass introduced; the two independently-scanned near-adjacent sections' apparent tissue
footprints genuinely differ by roughly that much (bounding-box area ratio of the raw point
clouds, 0.618, is in the same ballpark as scale²=0.535). Two candidate explanations — real
tissue shrinkage/deformation between near-adjacent cuts during separate freezing/mounting
steps (a documented effect in serial histology; `results/landscape_pathology_3d.md`'s CODA
summary discusses comparable cross-section effects), or the CODEX scan capturing a smaller ROI
within the block than Visium's — are both plausible and **neither is resolved in this pass**.

**Visual check** (`results/gbm_3d_registration_overlay.png`, generated and inspected before
trusting the numbers above, per this project's standing verification discipline): before
registration the two point clouds already occupy a broadly overlapping region (not offset to
different quadrants — the two platforms' raw coordinate conventions are compatible enough to
start from); after registration the CODEX cloud is visibly rescaled into tighter alignment
with the Visium spot grid; the tissue-outline panel shows the Visium outline (spot-grid
boundary, smooth) and the CODEX outline (real tissue silhouette, irregular/branching) occupying
mostly the same region without a perfect match — expected, since a similarity transform cannot
fix genuine local tissue deformation, the same limitation the breast pass's "Honest next step"
already flagged.

### Comparison against the paper's own r=0.73 — not apples to apples, stated precisely

`results/landscape_gbm_serial_sections.md` flagged the paper's own reported
"r=0.73, p=0.011" CODEX/Visium registration-quality number as an external sanity check to
compare against. Read closely (direct fetch, this pass's earlier research pass, PMC11088502):
that r=0.73 is the **Pearson correlation, across near-adjacent sample pairs, of the relative
abundance of "structured" vs "disorganized" spots/cells** between the two modalities **after**
the authors' own STalign diffeomorphic registration — a downstream biological-organization
concordance statistic, not a spatial point-cloud alignment-quality metric. This pass's density
correlation (0.272) and Chamfer/Hausdorff distances are a **different kind of number**
(geometric tissue-outline overlap, not cross-modality tissue-organization agreement), computed
with a different, non-diffeomorphic method (rigid+scale ICP, not LDDMM), on cell/spot POSITION
alone (this pass could not recover per-cell CODEX marker identity, so no "structured vs
disorganized" call is even possible here — see "CODEX layer" above). **No direct numeric
comparison is offered, because none would be honest** — the two numbers answer different
questions. What can be said: this pass's registration produces a real, direction-correct
geometric improvement on the same real data the paper's own STalign step operated on, using a
cheaper and more limited method (rigid/similarity ICP vs the paper's diffeomorphic warp), which
is consistent with (not proof of) the paper's own finding that these two near-adjacent sections
of ZH1007 genuinely can be brought into better spatial agreement.

## What this DOES and does NOT establish

**Does establish:**
- The Greenwald/Tirosh 2024 GBM Visium layer ingests through `tmesim/data/ingest.py`
  completely unchanged, at full (1.0) detector confidence, once its flat Zenodo tarball layout
  is restructured (symlinks, no new Python) into a Space Ranger `outs/` shape — confirming
  `results/landscape_gbm_serial_sections.md`'s "plausible, not verified" prediction as TRUE.
- A whole-transcriptome (33,538-gene) real human GBM Visium section, with MGMT, hypoxia, and
  vasculature all independently measurable — a second, independent confirmation of
  `results/gbmspace_adapter.md`'s finding that a richer-than-Xenium panel unlocks these axes,
  now on data from a completely different lab/pipeline (Tirosh/Suva vs Sanger/GBM-Space).
- A new, scoped, real ingester for CODEX/PhenoCycler OME-TIFF output
  (`tmesim/data/adapters/codex.py`) that correctly recovers real physical pixel scale,
  real approximate nucleus positions, and the paper's own declared 39-marker antibody panel
  (confirming CD31/PECAM1 and 3 of 5 hypoxia markers are on it), with 11 passing tests.
  Registration-3d.py's existing, UNMODIFIED multi-start ICP method genuinely improves
  cross-platform (Visium-spot vs CODEX-nucleus) tissue-frame overlap on this real GBM data by
  the same class of honest, before/after metrics the breast pass used (Chamfer −40.2%,
  Hausdorff-95 −36.8%, density correlation +10.8%) — larger relative improvements than the
  breast same-modality pair saw, because this cross-platform pair had a real scale mismatch to
  correct that the breast Rep1/Rep2 pair did not.
- One small, justified, additive fix to `registration_3d.py` (`z_spacing_caveat` as a
  parameter) so the function's documentation string is honest when reused on a second,
  different dataset — the registration algorithm itself is byte-identical to the breast pass.

**Does NOT establish:**
- **Per-cell CODEX marker identity.** This is the pass's sharpest honest limitation: the
  OME-TIFF's 43 raw channels carry dye names, not antibody names, and the manifest that would
  resolve that mapping was not found as a small standalone file in this deposit and was not
  fetched from the multi-GB QuPath projects (disk-budget scoped out). No CD31+ vessel, no
  hypoxic-CA9+ region, no per-cell phenotype call can be made from this pass's CODEX
  extraction — only DECLARED-panel presence (yes, the panel covers these markers) is
  established, never a measured per-cell value.
- **A validated CODEX cell segmentation.** `extract_nucleus_positions()` is DAPI local-maxima
  detection, explicitly not a validated segmentation — no cell boundaries, no doublet
  handling, real over/under-counting expected. It is real enough for the geometric
  registration this dataset was fetched to test, and no further claim is made for it.
  Genuine consequence for the `capability_report`-style reporting on the CODEX layer: because
  no per-cell marker-identity-resolved expression matrix exists, EVERY detection-RATE-dependent
  verdict reads UNMEASURED, by construction, not by an oversight.
- **Cell-level or true anatomical correspondence** between the Visium and CODEX sections,
  exactly as `results/registration_3d_breast_poc.md` already establishes for its own pair —
  "near-adjacent" sections of the same block do not contain the same cells.
- **The exact z-spacing** between the Visium and CODEX sections. Genuinely unpublished (the
  paper states "near-adjacent," not a micron figure) — `z_spacing_um=10.0` here is the SAME
  class of labelled placeholder the breast pass used, with its own dataset-specific caveat
  text (not the breast pass's reused verbatim), per the `z_spacing_caveat` fix above.
- **A resolved explanation for the 27% scale mismatch** between the two sections' apparent
  tissue footprints — flagged, with two plausible candidate explanations, neither confirmed.
- **A full drug-transport run on this data.** Out of scope for this pass, as the task allowed;
  the Visium layer's capability report above (MGMT/hypoxia/vasculature all OK) is what a
  future `reconstruct_from_visium`-style run on this specific section could build on.
- **Multi-tumour registration.** One tumour (`ZH1007nec`) of the 13 available, per the
  explicit "prove the pipeline on 1-2 first" disk-budget instruction. The `ZH811`/`ZH881`
  naming ambiguity noted above would need resolving (via `visium_dataset_description.csv`
  cross-referenced against the CODEX filenames more carefully, or the QuPath project's own
  sample manifest) before a second tumour is attempted.

## Honest next step

In priority order:

1. **Resolve the CODEX channel-to-marker manifest.** Either download one QuPath project
   archive (2.2 GB — needs a disk-budget decision beyond this pass's scope) and extract its
   channel-naming metadata (QuPath projects store this in a readable `project.qpproj` JSON
   file alongside the unreadable per-image `.qpdata`, per `results/landscape_pathology_3d.md`
   — worth checking whether `project.qpproj` alone is extractable without the whole archive),
   or contact the paper's authors (`github.com/tiroshlab/Spatial_Glioma` — check whether the
   channel/cycle sheet is published there; not checked in this pass) for the raw manifest.
   Once available, `codex.py`'s `channel_intensities` array (already extracted, already kept)
   becomes a real per-cell marker expression table with almost no further work.
2. **A real cell segmentation**, not DAPI local-maxima — e.g. StarDist or Cellpose on the DAPI
   channel, both already common CODEX-analysis choices, would materially improve the nucleus
   position estimate this pass's registration already benefits from.
3. **Register more tumours** once the `ZH811`/`ZH881` naming ambiguity is resolved, to see
   whether the ~27% scale-mismatch finding is a one-off for this particular necrotic-region
   pair or a systematic Visium-vs-CODEX effect worth understanding before trusting future
   cross-platform registrations at face value.
4. **A non-rigid registration step** (thin-plate-spline or STalign-style diffeomorphic warp)
   after the current rigid/similarity stage — the same "Honest next step" #4 the breast pass
   already recommended, now doubly motivated: this cross-platform pair's tissue-outline panel
   (`results/gbm_3d_registration_overlay.png`, right panel) visibly shows local deformation a
   similarity transform cannot fix.

## Files

- `tmesim/data/adapters/codex.py` — the new CODEX OME-TIFF adapter (additive, new file).
- `tmesim/reconstruct/registration_3d.py` — one additive, backward-compatible parameter
  (`z_spacing_caveat`) added to `register_serial_sections`; registration algorithm unchanged.
- `pyproject.toml` — `imagecodecs` added to core dependencies (required to decode this real
  file's LZW compression).
- `tests/test_codex_adapter.py` — 11 tests: OME-XML metadata reading (including the real
  filename-mismatch quirk), nucleus-position extraction (determinism, physical bounds, a
  pinned real-data regression value), the declared-panel capability report, and one real
  end-to-end Visium-vs-CODEX registration test. All pass.
- `tests/test_registration_3d.py` — unchanged, all 9 tests still pass (confirms the
  `z_spacing_caveat` change is behaviourally invisible to existing callers).
- Full suite (`./.venv/bin/python -m pytest tests/ -q`): **552 passed, 0 failed** run before
  `tests/test_codex_adapter.py` existed (confirms the `registration_3d.py`/`pyproject.toml`
  changes alone broke nothing); **577 passed, 0 failed** run after it was added. The +25 delta
  is larger than this pass's own 11 new tests — this repo had other, unrelated concurrent work
  landing in the same working tree during this pass (`git status` showed in-progress edits to
  `app/web/js/*`, `tmesim/reconstruct/cell_state.py`, `tmesim/reconstruct/real_pipeline.py`,
  and a new `tmesim/ml/` module, none of which this pass touched or reviewed); both full-suite
  numbers are reported exactly as measured rather than adjusted to isolate this pass's own
  contribution, since a git-stash-based clean isolation was not attempted precisely so as not
  to disturb that concurrent work. This pass's own 11 tests were also confirmed green in
  isolation (`pytest tests/test_codex_adapter.py -q` → 11 passed).
- `data_cache/greenwald_gbm_visium/GBM_ZH1007nec/` — the real, live-downloaded Visium tarball
  contents plus the symlinked `outs/` restructuring (~22 MB).
- `data_cache/greenwald_gbm_codex/ZH1007_NEC_B.ome.tif` — the real, live-downloaded CODEX
  OME-TIFF (~921 MB).
- `results/gbm_3d_registration_overlay.png` — the verification plot referenced above.

---

## Second case: ZH1007's infiltrative (inf) region — does the first result generalize?

**One-line summary.** Same tumour (`ZH1007`) as above, its OTHER region: the infiltrative
margin (`inf`) instead of the necrotic core (`nec`). Chosen deliberately over an unrelated
tumour, per the task's own preference, because it is the one comparison in this deposit that
isolates *region* as the variable while holding *patient* constant. Both `tmesim/data/ingest.py`
and `tmesim/data/adapters/codex.py` handled this second real file completely unmodified — no
new code, no format surprises. **The direction of the finding repeats: registration measurably
improves cross-platform tissue-frame overlap on this second real pair too** (Chamfer −38.5%,
Hausdorff-95 −50.3%, density correlation +31.6%). The *magnitude* does not repeat exactly —
Hausdorff-95 and especially density correlation improve substantially MORE here than on `nec`,
and the absolute (post-registration) overlap quality is worse in absolute density-correlation
terms, consistent with the infiltrative margin genuinely being a more diffuse, less
sharply-bounded tissue structure than the necrotic core. Both are honest, real findings — stated
precisely below, not smoothed into one number.

### What is confirmed real, on disk, unchanged from the earlier checkpoint
- Chose `GBM_ZH1007inf.tar.gz` (Visium, 20.8 MB) + `ZH1007_INF.ome.tif` (CODEX, 1548.8 MB) —
  the SAME tumour (`ZH1007`) already registered in its necrotic (`nec`) region above, now its
  infiltrative (`inf`) margin, per the task's preferred "same tumour, other region" design
  (directly comparable pair, not a new unrelated tumour). Confirmed via
  `visium_dataset_description.csv` (fetched fresh from the Zenodo API, not assumed) that
  `GBM_ZH1007inf.tar.gz` is a real, distinct Visium sample entry (`zh1007inf_*` files),
  same naming convention as `zh1007nec_*`.
- Visium tarball downloaded, extracted, and restructured into the same symlinked
  `outs/{filtered_feature_bc_matrix,spatial}/` shape used for `ZH1007nec`
  (`data_cache/greenwald_gbm_visium/GBM_ZH1007inf/outs`). `ingest.detect()` returns
  `platform=visium, layout=spaceranger_outs, confidence=1.0` — same full-confidence result as
  the first pair, zero new code. `ingest.ingest()` gives **1,436 spots, 33,538 genes** (a
  richer spot count than `nec`'s 955). `capability_report` confirms the same whole-transcriptome
  verdicts as `nec`: vascular_field OK, cell_typing OK (coverage 1.00), hypoxia_measured OK,
  mgmt_from_panel OK, ecm_field OK — all four `[OK]` axes reproduce on this second section.
  One genuine format note: this file's `matrix.mtx` header reads
  `{"software_version": "v16.10-2-g378ed53-dirty", ...}`, NOT `"spaceranger-1.1.0"` like `nec`'s
  — a different processing-pipeline version string across the two regions of the same tumour,
  recorded here in case it matters later; it did not affect ingestion (still a valid
  MatrixMarket file, same 33,538×N shape convention, `resolve_visium_dirs`/`detect` do not key
  off this string).
- CODEX OME-TIFF fully downloaded (1,548,785,153 bytes, matches the Zenodo listing's 1548.8 MB
  to within rounding) and read with the EXISTING, unmodified `codex.read_ome_metadata()`:
  43 channels, 10584×9874 px (a different pixel footprint than `nec`'s 8205×7330 — a
  DIFFERENT physical ROI, expected for a different region of the same block), 0.50670 µm/px
  (IDENTICAL calibration to `nec`, both regions scanned on the same instrument/settings), 3
  pyramid levels — same format shape as `nec`, adapter needed zero changes. Channel 0 is DAPI,
  channels 1–42 cycle `ATTO 550`/`FITC`/`Cy5` — same dye-name-not-marker-name situation as
  `nec`, confirming this is a per-file property of the deposit's raw CODEX export, not a
  one-off. The self-referencing-filename OME-XML quirk reproduces here too, but points to a
  DIFFERENT stale name than expected: `ZH1017_NEC.ome.tif` (patient 1017, region NEC) embedded
  inside a file Zenodo actually serves as `ZH1007_INF.ome.tif` (patient 1007, region INF) — i.e.
  BOTH the patient-ID digit swap (1007↔1017) AND a region swap (INF↔NEC) are present in the
  stale metadata, consistent with the depositors reusing/cross-referencing an internal naming
  scheme sloppily across multiple files, not just the one `nec` pass 1 already flagged.
- Disk: 22 GiB free before this second pair was started, 18–21 GiB free during the ~1.57 GB
  combined download (other concurrent work in this shared repo is also using disk — confirmed
  via `git status` before starting, in-progress edits outside this pass's scope, same as pass
  1 noted). No cleanup needed yet.
- Baseline full suite (before any new code/config change for this pass — none has been made
  yet, since no adapter changes were needed): `577 passed, 0 failed`
  (`./.venv/bin/python -m pytest tests/ -q`, matches pass 1's post-pass count exactly, confirming
  no drift from other concurrent work between passes).

### Nucleus extraction and registration — completed

Full-resolution (`level=0`) DAPI local-maxima extraction on `ZH1007_INF.ome.tif`
(10584×9874 px, ~1.7× the pixel count of `nec`'s 8205×7330) took **29.4 s** (faster in
wall-clock terms than `nec`'s 34.4 s despite more pixels — machine load from other concurrent
work in this shared repo varies run to run; not a meaningful comparison) and found
**39,100 nucleus positions**, spanning real micron coordinates (52.2–4936.8 µm ×
46.1–5345.7 µm) — consistent with the image's own physical field of view
(9874×0.50670 = 5003 µm, 10584×0.50670 = 5363 µm). A comparable order of magnitude to `nec`'s
40,340, on a larger imaged footprint (bounding-box area 25.9 mm² vs `nec`'s ~14.4 mm²) —
i.e. a genuinely LOWER nucleus areal density on this section (39,100 / 25.9 mm² ≈
1,511 nuclei/mm² vs `nec`'s 40,340 / 14.4 mm² ≈ 2,802 nuclei/mm²), which is the expected
direction for an infiltrative margin (tumour cells more sparsely interspersed with normal
brain parenchyma) versus a necrotic core rim — recorded as a plausible, not proven,
biological read of a purely geometric measurement (no cell typing was done; see limitations).

Same real quirk reproduces, with a twist: this file's OME-XML `TiffData`/`UUID` blocks
self-reference `"ZH1017_NEC.ome.tif"` — NOT `"ZH1007_INF.ome.tif"`, the name Zenodo actually
serves it under. That is a DIFFERENT stale filename than `nec`'s own self-reference
(`"ZH1017_INF_B.ome.tif"`), and notably **the patient-digit swap (1007↔1017) is now paired
with a REGION swap too** (this INF file's embedded name says NEC; `nec`'s embedded name said
INF) — i.e. the two files' stale internal names look like they may be cross-referencing EACH
OTHER's region, not just carrying a generic typo. Recorded exactly as found; not chased
further (does not affect channel 0 = DAPI or anything this adapter computes).

`register_serial_sections()` — called with the exact same signature, defaults, and
`z_spacing_um=10.0` placeholder convention as the `nec` pass, only `z_spacing_caveat`'s text
changed to name this second, different pair explicitly (not reused verbatim from `nec`):

```
transform:  rotation = 4.279°   scale = 0.6156   translation = (484.0, 892.0) um
z_spacing_um = 10.0  (ASSUMED — same placeholder class as nec; not a measured quantity)
n_points = 40,536  (1,436 Visium spots @ z=0  +  39,100 CODEX nuclei @ z=10.0)
```

| metric | before | after | change |
|---|---|---|---|
| density correlation (60 µm grid) | 0.1370 | 0.1803 | **+31.6%** relative |
| Chamfer distance, tissue outline (µm) | 288.9 | 177.7 | **−38.5%** (closer) |
| Hausdorff-95, tissue outline (µm) | 1342.3 | 667.1 | **−50.3%** (closer) |

Winning candidate was again `seed_index=0` (identity rotation start), but unlike `nec`
(0.088° rotation, near-identity), this pair's winning transform has a real **4.28° rotation**
AND a **larger scale correction (0.616, i.e. ~38.4% linear shrink of the CODEX cloud)** than
`nec`'s 0.731 (~26.9% shrink). Six runner-up candidates (seed indices 2–6, all coarse-search
seeds converging to −105° to −120° rotations) produced WORSE Chamfer (219.3–231.8 µm vs the
winner's 177.7 µm) despite comparable or even higher raw density correlation (0.428–0.450 vs
the winner's 0.180) — the same "a wrong candidate can look locally competitive on density
correlation alone" pattern `nec`'s writeup already documented, Chamfer-based final selection
catching it correctly here too. Full table in the script output below.

**On the scale correction, checked the same way `nec`'s was, not asserted by analogy.** Visium
pixel calibration for this section: `55 µm / 112.3366 px = 0.4897 µm/px` (vs `nec`'s
`0.4905 µm/px` — 0.16% apart, trivial). CODEX pixel calibration: `0.50670 µm/px`, IDENTICAL to
`nec`'s to 5 decimal places (both regions scanned on the same instrument/settings, as
expected for one deposit). Visium-vs-CODEX calibration gap here: 3.5% (`nec`: 3.2%) — nowhere
near the 38% ICP-found scale correction, so, exactly as `nec` found, **this is not a units bug**.
Raw point-cloud bounding-box area ratio (Visium / CODEX) is 0.510 here (linear √ratio 0.714) —
in the same rough ballpark as the ICP scale (0.616) without matching it exactly, the same
"same ballpark, not exact" relationship `nec`'s writeup reported for its own 0.618 bbox ratio
vs 0.731 ICP scale. Two candidate explanations carried over unchanged from `nec` (real tissue
shrinkage/deformation between near-adjacent cuts vs a smaller CODEX-scanned ROI within the
block) — **still not resolved**, and this second data point does not resolve which one
dominates; if anything, a LARGER scale correction on `inf` than `nec` is consistent with either
explanation still being live (more deformation OR a differently-sized scan ROI on this second
section, no way to tell from geometry alone).

**Declared-panel capability check**: identical to `nec`'s, because it is the same antibody
panel/deposit, not a per-region property — `codex.declared_panel_capability_report()` on this
section confirms the same `vascular_field: OK` (CD31→PECAM1) and `hypoxia_measured: OK`
(CA9, NDRG1, GLUT1→SLC2A1) verdicts, `mgmt_from_panel: NO` (same honest negative). Not
re-derived from scratch as a new finding — recorded to confirm it reproduces, which it does
trivially (same panel PDF, same alias table, same file-independent computation).

**Visual check** (`results/gbm_3d_registration_zh1007inf_overlay.png`, generated and inspected
before trusting the numbers above, same discipline as `nec`): before registration, the raw
CODEX cloud is visibly LARGER than and offset to the right of the Visium spot footprint — a
much more visible pre-registration mismatch than `nec`'s already-broadly-overlapping start;
after registration, the CODEX cloud is rescaled and shifted into much closer alignment with
the Visium spot grid's shape; the tissue-outline panel shows the CODEX outline (irregular,
branching, genuinely infiltrative-looking) occupying mostly the right/interior portion of the
Visium footprint, a real but imperfect overlap — the same similarity-transform limitation
(cannot fix local deformation) `nec` already flagged, visibly present here too.

### Does the first result generalize? Direct comparison

| metric | ZH1007 **nec** (pass 1) | ZH1007 **inf** (pass 2) | same direction? |
|---|---|---|---|
| Chamfer, relative change | −40.2% | −38.5% | **yes, close** |
| Hausdorff-95, relative change | −36.8% | −50.3% | yes, but **substantially larger** here |
| Density correlation, relative change | +10.8% | +31.6% | yes, but **~3× larger** here |
| Density correlation, absolute (after) | 0.272 | 0.180 | **inf is LOWER in absolute terms** |
| ICP scale correction | 0.731 (27% shrink) | 0.616 (38% shrink) | same direction, **larger** on inf |
| Winning-candidate rotation | 0.088° (~identity) | 4.28° (small but real) | inf needed more rotation too |

**Honest read: the DIRECTION generalizes cleanly (all three metrics improve on both real
sections, by a comparable order of magnitude), but the MAGNITUDE does not repeat exactly, and
one number moves enough to be worth flagging rather than averaging away.** Chamfer is close
between the two (−40.2% vs −38.5%) — the strongest case for "this is a stable effect of the
method, not a lucky first case." Hausdorff-95 and density correlation both improve MORE on
`inf` than on `nec` in RELATIVE terms, but `inf`'s density correlation remains LOWER than
`nec`'s in ABSOLUTE terms even after registration (0.180 vs 0.272) — i.e. `inf` starts from a
much worse baseline (0.137 vs 0.246) and closes more of the gap, without ever reaching `nec`'s
absolute overlap quality. This is consistent with a real, non-arbitrary biological difference
between the two regions this task asked about: **the infiltrative margin's tissue structure is
genuinely less spatially concentrated/coherent than the necrotic core's**, in exactly the
direction a pathologist would expect (necrotic cores are compact, sharply demarcated
structures; infiltrative margins are, by definition, diffuse zones where tumour cells
interdigitate with normal brain) — the lower absolute density correlation and the visibly more
irregular CODEX outline in the overlay figure are two independent pieces of evidence pointing
the same way. **This pass does not claim to have proven that biological explanation** — only
two regions of one patient were measured, not enough to separate "biology" from "this
particular pair of scans" — but it is the most parsimonious reading of what was actually
measured, stated as a plausible interpretation, not a conclusion.

### What this second case does and does not add

**Does establish:**
- `tmesim/data/ingest.py` and `tmesim/data/adapters/codex.py` handle a SECOND real file pair
  from this same deposit completely unmodified — the `nec` pass's "zero new ingestion code"
  and "genuinely new CODEX format, but the adapter is now written" findings both hold on a
  second real file, not just the one they were built against. Confirms the adapter generalizes
  within-deposit, which was not guaranteed a priori.
- The `nec` pass's cross-platform registration improvement (Chamfer/Hausdorff/density
  correlation, all improving) is not a one-off lucky result on a single tumour/region — a
  second real section from the same patient shows the same directional effect.
- A genuine, geometrically-measured difference between necrotic-core and infiltrative-margin
  tissue structure within the same patient — lower absolute post-registration density
  correlation and a visibly more irregular tissue outline on the infiltrative section — that a
  future pass with real cell typing (not just DAPI positions) could investigate further.
- The channel-dye-not-marker-name problem, the OME-XML stale-filename quirk, and the ~part-in-30
  Visium-vs-CODEX pixel-calibration gap all reproduce on a second, independently-decoded file —
  each is a property of this DEPOSIT's export pipeline, not a one-off artefact of the first file
  fetched.

**Does NOT establish:**
- **Per-cell CODEX marker identity** — unresolved here for the exact same reason as `nec` (dye
  names, not antibody names, and the manifest was not fetched this pass either — same disk-budget
  reasoning as pass 1, now doubly true with two large OME-TIFFs already on disk).
- **A validated CODEX cell segmentation** — same DAPI-local-maxima caveat as `nec`, unchanged.
- **Multi-tumour generalization** — this is still ONE patient (`ZH1007`), two regions of it.
  The `ZH811`/`ZH881`/`ZH8811` naming ambiguity `nec`'s writeup flagged is still unresolved; a
  genuinely different tumour (not just a different region of the same one) has not yet been
  registered. This pass answers "does the SAME patient's other region look similar?" (yes,
  directionally) but not "does a DIFFERENT patient look similar?" (still open).
- **Why `inf`'s scale correction and absolute density correlation differ from `nec`'s** — two
  plausible explanations offered (real infiltrative-vs-necrotic tissue structure difference,
  or a scanning/ROI artefact independent of biology), neither ruled out.
- **z-spacing** — same unverified placeholder (`z_spacing_um=10.0`) as `nec`, for the same
  reason (paper states "near-adjacent," not a micron figure).

### Honest next step (updates `nec`'s list, does not replace it)

1. Register a genuinely different tumour (not another region of `ZH1007`) — `ZH916` or
   `ZH1019` both have clean same-patient Visium/CODEX filename matches per
   `visium_dataset_description.csv`/`readme_CODEX.txt` (unlike the `ZH811`/`ZH881`/`ZH8811`
   ambiguity) and would separate "does registration generalize across REGIONS of one patient"
   (this pass's question, now answered directionally-yes) from "does it generalize across
   PATIENTS" (still open).
2. Resolve the CODEX channel-to-marker manifest — unchanged priority from `nec`'s list; now
   twice as valuable, since it would unlock a REAL comparison of e.g. CD31+ vascular density or
   CA9+ hypoxic fraction between the necrotic and infiltrative regions of the same tumour,
   instead of only a geometric one.
3. A real cell segmentation (StarDist/Cellpose on DAPI) — unchanged priority; would also let
   the "lower nucleus areal density on `inf`" observation above be checked against an actual
   segmentation instead of a local-maxima proxy that could itself behave differently on more
   diffuse tissue.

### Files (this second case only; `nec`'s files list above is unchanged)

- `data_cache/greenwald_gbm_visium/GBM_ZH1007inf/` — the real, live-downloaded Visium tarball
  contents (20.8 MB) plus the same symlinked `outs/` restructuring pattern used for `nec`.
- `data_cache/greenwald_gbm_codex/ZH1007_INF.ome.tif` — the real, live-downloaded CODEX
  OME-TIFF (1,548,785,153 bytes ≈ 1548.8 MB, matches the Zenodo listing).
- `results/gbm_3d_registration_zh1007inf_overlay.png` — the verification plot for this pair.
- **No changes to `tmesim/data/adapters/codex.py`, `tmesim/reconstruct/registration_3d.py`, or
  any other source file** — both existing, unmodified modules from the `nec` pass handled this
  second real file with zero new code. No new test file was added for this pass (no new code
  exists to pin); the existing `tests/test_codex_adapter.py` and `tests/test_registration_3d.py`
  continue to cover the adapter/registration logic this pass exercised on new data.
- Test suite: the **577 passed, 0 failed** baseline recorded in the checkpoint above (run
  before this second pair's ingestion/registration script was written, confirming the starting
  state) still applies — no source file changed during this pass, so nothing invalidates it. A
  second full-suite run was deliberately NOT performed for this pass, per explicit
  coordinator instruction (a central pass verifies and commits across all concurrently-running
  agents' work rather than each pass re-running the full ~5-minute suite independently).
