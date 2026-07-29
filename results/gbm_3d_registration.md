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
