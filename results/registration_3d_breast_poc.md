# Serial-section 3D registration: proof-of-concept on real breast-cancer Xenium data

**One-line summary.** Two real, independently-reconstructed 5 µm Xenium sections of
the same human breast-cancer FFPE block (Janesick et al. 2023, Rep1/Rep2) were
registered into one shared (x, y, z) frame with a coarse-search-seeded, multi-start,
trimmed-ICP similarity-transform alignment. It works, ingests with zero new
ingestion code, and genuinely improves tissue-frame overlap by real (modest, honestly
reported) margins — but along the way the first version of the algorithm was caught
converging to a **worse**, plausible-looking alignment, and that failure (and the fix)
is reported here as carefully as the final numbers, per this task's own instructions.

## What this is, precisely

- **Product framing.** `tmesim.core.Domain` is (ny, nx) — no z-axis anywhere in the
  simulation engine. tme-sim's stated 3D strategy is to build the third dimension
  itself by registering multiple 2D sections cut at known z-spacing through the same
  tissue block (the PASTE2/STalign approach), not to wait for a vendor 3D format. This
  is the first concrete implementation of that strategy.
- **New module, fully additive.** `tmesim/reconstruct/registration_3d.py`. Nothing in
  `tmesim/core`, `tmesim/physics`, `tmesim/sim`, `tmesim/molecular`, or the existing
  2D `reconstruct` path (`assemble.py`, `spatial.py`, `vasculature.py`, ...) was
  touched. The GBM flagship pipeline and its cached results are untouched.
- **What "registration" means here — read this before trusting any number below.**
  Two adjacent 5 µm Xenium sections of the same block do **not** contain the same
  cells (a cell in Rep1 is physically absent from Rep2 — it was sectioned through or
  sits one slice away). There is no cell-to-cell ground truth to recover, and nothing
  in this module claims to have found one. "Registration" means finding the
  **rigid-plus-uniform-scale (similarity) transform of the (x, y) tissue coordinate
  frame** of Rep2 that best overlaps its **tissue geometry** (cell density / outline)
  with Rep1's — the same class of alignment PASTE2 calls its "spatial" step and
  STalign does with diffeomorphic image registration, done here with a simpler and
  more honestly-scoped tool.
- **What the z-axis is.** Rep1 is placed at z=0, Rep2 at z=z_spacing_um. This is a
  **2.5D reconstruction** (stacked, registered 2D sections with an assumed constant
  depth offset), not a volumetric one — there is no interpolation, no continuous
  volume, no third section to constrain a real z-gradient.

## The data

Both sections are real, live-downloaded 10x Xenium output from Janesick et al. 2023
(*Nat Commun* 14, 8353; PMID 38114474) — the human breast-cancer FFPE block that
STalign itself used to validate registration. Per the task's disk-budget constraint,
only the lightweight component files were downloaded (~100 MB total, not the
9-10 GB `outs.zip` bundles with morphology imagery):

```
https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_cells.csv.gz
https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_cell_feature_matrix.tar.gz
https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_gene_panel.json
(swap Rep1 -> Rep2 for the second section)
```

Cached at `data_cache/xenium_breast_janesick/{rep1,rep2}/`.

**One real wrinkle, worth recording plainly.** The `..._cell_feature_matrix.tar.gz`
files 10x serves for this specific sample are, despite the filename and despite
`content-type: application/x-tar` (checked directly with `curl -I`), **not actually
gzip-compressed** — they are plain tar archives. The existing loader
(`data/real/singlecell_spatial.py`) opens the tarball with `tarfile.open(path,
"r:gz")`, which fails on this real file (`BadGzipFile`). This is a genuine
data-serving quirk, not a bug introduced here. It was worked around **without writing
or modifying any ingestion code**: the tarball was extracted in place with the system
`tar` (which auto-detects the format), producing an already-extracted
`cell_feature_matrix/` directory — the exact layout `resolve_xenium_files()` already
prefers over the tarball path. So the "zero new ingestion code" claim in the task
brief holds, but not quite trivially — it took one manual extraction step outside
Python, not a code change, to get there.

### Ingestion (existing `tmesim.data.ingest` module, unmodified)

```
rep1: detect() -> platform=xenium, layout=xenium_bundle_extracted, confidence=1.0
      ingest() -> n_cells=167,780  n_genes=313  extent_um=(7521.0, 5474.9)  validation.ok=True
rep2: detect() -> platform=xenium, layout=xenium_bundle_extracted, confidence=1.0
      ingest() -> n_cells=118,752  n_genes=313  extent_um=(7519.8, 5471.0)  validation.ok=True
```

167,780 cells for Rep1 matches the paper's reported 167,885 (main-text figure is for
Rep1; the small difference is expected — cell counts/QC can differ slightly by
pipeline version). No warnings, no validation failures, no code changes to
`tmesim/data/ingest.py` or its loaders.

## Method

`tmesim.reconstruct.registration_3d.register_serial_sections(xy_ref, xy_mov,
z_spacing_um, ...)`:

1. **Quality metrics, measured BEFORE any transform** (see "Quality metrics" below).
2. **Multi-start seeding.** Two "do-nothing" seeds (identity; centroid-only
   translation) plus the top-5 candidates from a coarse rotation/scale grid search
   (72 angles × 3 scale factors, ranked by Pearson correlation of the two point
   clouds' density fields on a shared 100 µm grid).
3. **Trimmed ICP refinement of every seed** (Besl & McKay 1992): each iteration finds
   nearest-neighbour correspondences for a fixed random cell subsample under the
   current transform, keeps the closest 80% of pairs (trimming is the concession to
   partial section-to-section overlap — the same problem PASTE2's "partial" addresses
   with optimal transport instead), and refits a similarity transform in closed form
   via the Umeyama (1991) SVD solution.
4. **Final candidate selection by the post-ICP boundary Chamfer distance** — **not**
   by the coarse correlation score. This is the single most important design decision
   in this module, and it exists because of a failure caught during verification (see
   next section).
5. **Quality metrics, measured AFTER**, identical definitions, so before/after is a
   fair comparison.
6. Assemble `(x, y, z)`: Rep1's own (x, y) at z=0, Rep2's registered (x, y) at
   z=z_spacing_um.

## The failure this task's honesty requirement caught — and the fix

The very first version of this module picked **only the single best-scoring** coarse
rotation candidate (density correlation on the coarse grid) and fed it straight to
ICP, no second-guessing. On this real dataset that produced a **20.0° rotation** with
density correlation apparently *improving* from 0.242 to 0.425 — a plausible-looking
"success."

Plotting it (`results/registration_3d_breast_poc_overlay.png`, left two panels of the
first draft, reproduced below in the current form) showed the opposite: the "before"
panel, with **no transform applied at all**, already had Rep1's and Rep2's lobular
tumour structures roughly co-located, and the 20°-rotated "after" panel visibly
**broke** that alignment. The boundary metrics agreed with the eye and disagreed with
the correlation score:

| | density correlation | Chamfer (µm) | Hausdorff-95 (µm) |
|---|---|---|---|
| before (identity) | 0.2416 | 211.5 | 768.7 |
| 20°-rotation candidate (first version's answer) | 0.4247 | 254.1 | 1099.2 |
| centroid-translation-only | 0.2448 | 208.6 | 634.9 |

**Diagnosis.** Breast tissue in this section is dense with similarly-sized, roughly
circular glandular/lobular structures repeated across the field of view. A *coarse*
density grid (100 µm bins) can alias on that repetition: a wrong rotation can score a
higher coarse correlation than the true alignment simply because it phase-matches one
lobule to a different, similar-looking lobule elsewhere in the section. The coarse
correlation metric was the wrong thing to trust as a final answer — it is fine as a
cheap *seed generator* (which is still what it is used for) but not as a *selector*.

**Fix (now in the module, and load-bearing).** `register_serial_sections` no longer
trusts a single seed. It runs ICP from seven candidates (identity, centroid-only, and
the top-5 coarse-search rotations) and picks the winner by the **post-ICP boundary
Chamfer distance**, a metric computed on much finer (30 µm) tissue-outline geometry
and empirically far less alias-prone than the coarse correlation grid. On the real
data this correctly demotes the aliased 20° candidate and selects the near-identity
one (full candidate table in `SerialSectionRegistration.candidates`, reproduced
below).

This is exactly the check the task asked for — *"does the registered Rep2 tissue
outline actually overlap sensibly with Rep1's when plotted? ... do not force a
good-looking number"* — and it changed the shipped algorithm, not just the writeup.

## Final result on the real Rep1/Rep2 pair

```
transform:  rotation = 0.010°   scale = 0.99995   translation = (65.2, 622.2) um
z_spacing_um = 10.0  (ASSUMED — see "z-spacing" below)
n_points = 286,532  (167,780 Rep1 @ z=0  +  118,752 Rep2 @ z=10.0)
```

| metric | before | after | change |
|---|---|---|---|
| density correlation (60 µm grid) | 0.2416 | 0.2448 | **+1.35%** relative |
| Chamfer distance, tissue outline (µm) | 211.5 | 205.7 | **−2.73%** (closer) |
| Hausdorff-95, tissue outline (µm) | 768.7 | 634.5 | **−17.5%** (closer) |

Full multi-start candidate table (7 seeds tried; winner is seed_index=1):

| seed | seed rotation | final rotation | final scale | final Chamfer (µm) | final density corr |
|---|---|---|---|---|---|
| 0 identity | 0.0° | −0.004° | 1.0000 | 209.5 | 0.2416 |
| **1 centroid-only (WINNER)** | **0.0°** | **0.010°** | **0.99995** | **205.7** | **0.2448** |
| 2 coarse-search top-1 | 20.0° | 19.997° | 0.99993 | 254.1 | 0.4247 |
| 3 coarse-search #2 | 15.0° | 14.992° | 0.99992 | 226.0 | 0.4001 |
| 4 coarse-search #3 | 25.0° | 24.992° | 0.94988 | 284.9 | 0.4241 |
| 5 coarse-search #4 | 20.0° | 20.011° | 0.94977 | 255.5 | 0.3983 |
| 6 coarse-search #5 | 25.0° | 24.999° | 0.99987 | 285.4 | 0.4359 |

Visual check: `results/registration_3d_breast_poc_overlay.png` — left panel (raw,
unregistered coordinate frames), middle panel (after the winning near-identity
transform), right panel (tissue outlines, boundary points from
`tissue_boundary_points`). Rep1's cell coordinates were left untouched (they define
the shared frame); only Rep2 was transformed.

**Reading the result honestly.** The improvement is real but modest, and it is mostly
a translation, not a rotation or a scale correction: the winning transform is a 622 µm
shift in y with essentially zero rotation (0.01°) and essentially zero scale
correction (0.99995). The most plausible explanation is that Xenium's own per-run
coordinate export already places both replicates in a broadly consistent frame (same
slide/imaging-area convention), so most of the "registration problem" for *this
particular pair* was already solved before this module ran, and what remained to find
was a modest translational offset (plus a large chunk of genuinely un-overlapping
tissue at the top/bottom of each section's field of view — visible in the right panel
— which no rigid/similarity transform can fix, because it reflects the two sections
literally not covering the same area of tissue, not a misalignment). The boundary
metric used for candidate selection is also partly measuring **capture-window
overlap**, not purely biological tissue silhouette: both sections' cells fill almost
their entire rectangular scan field of view, so `tissue_boundary_points` (occupancy
threshold + morphological closing) traces something close to the scan rectangle for
much of the perimeter, not a fine anatomical contour. That is a real limitation of the
metric, stated plainly rather than left implicit.

## A second honest finding: the ICP basin of convergence is narrow on this data

Verifying the core ICP/Umeyama machinery on real cell geometry (not the ambiguous
cross-section case — a KNOWN synthetic rotation/scale/translation was applied to a
subsample of Rep1's actual centroids, so an exact correspondence genuinely exists)
turned up a second limitation worth recording. Seeded from the true answer, ICP
converges to it almost exactly (median residual ~1e-11 µm after 1-2 iterations,
`tests/test_registration_3d.py::test_icp_recovers_known_transform_given_a_reasonable_seed`).
But seeded only ~2-3° off the true rotation, on this dense, roughly homogeneous cell
point cloud (no distinguishing per-point features, fairly uniform nearest-neighbour
spacing), trimmed ICP can converge to a **self-consistent but wrong** local optimum
rather than the true correspondence — confirmed empirically: a true injected rotation
of 37° with a seed 3° off converged to −40° rather than the true −37°, and did not
improve with more iterations or a larger sample (it is a genuine local optimum, not
slow convergence). This is the classic ICP failure mode for feature-poor point
clouds, and it is why step 4 above (multi-start plus a robust final-candidate metric)
matters — the current coarse grid (5° steps in the default config) is not fine enough
to guarantee landing inside the true optimum's basin for an *arbitrary* large
rotation. It happened not to bite on the real Rep1/Rep2 pair only because that pair's
true offset turned out to be close to identity, comfortably inside the basin from the
identity/centroid seeds. **A future pair with a genuinely large rotation between
sections is not yet safely handled by this implementation** — see "Honest next step."

## What this DOES and does NOT establish

**Does establish:**
- The Janesick et al. 2023 breast-cancer Xenium serial-section pair ingests cleanly
  through tme-sim's existing generic ingester, with no new ingestion code (modulo the
  one-time manual tar-extraction workaround for a real upstream data-serving quirk,
  described above).
- A first, working, additive implementation of tme-sim's stated 3D strategy — 2D
  tissue-frame registration via coarse-search-seeded, multi-start, trimmed ICP — runs
  end to end on two real, independent cancer tissue sections and produces a genuine
  (x, y, z) point cloud.
- The registration is measurably, if modestly, better than doing nothing on real
  tissue geometry by two independent, honestly-defined metrics (density-field
  correlation, boundary Chamfer/Hausdorff distance), computed identically before and
  after.
- The verification step the task demanded — plot it, check an objective metric, don't
  trust the algorithm blind — caught a real failure mode (coarse-grid aliasing on
  repetitive tissue) before it shipped, and the fix is now load-bearing in the
  algorithm, not just noted in prose.

**Does NOT establish:**
- **Cell-level 3D correspondence.** No claim is or can be made that any specific cell
  in Rep1 corresponds to any specific cell in Rep2. The (x, y, z) point cloud is two
  independently-measured cell populations placed in a shared frame, not a tracked
  3D cell.
- **A precise z-spacing.** `z_spacing_um` is an explicit, labelled placeholder
  (`Z_SPACING_UM_ASSUMED = 10.0`, plausible range roughly 5-15 µm) because Janesick et
  al. 2023 documents 5 µm section thickness and that the two replicates are serial,
  but not the exact number of microns separating Rep1 from Rep2 specifically
  (`results/landscape_serial_sections.md`). Every call site must pass this
  explicitly — there is no silent default presented as measured.
- **A volumetric reconstruction.** This is 2.5D: two stacked, registered planes, not
  a continuous interpolated 3D density/expression field. There is no third section to
  even check whether a linear z-interpolation would be reasonable.
- **Robustness to a large true rotation between sections**, per the ICP basin-of-
  convergence finding above — untested and, based on the synthetic-recovery
  experiment, not currently reliable for this dataset's point-cloud density/structure.
- **A non-rigid or optimal-transport alignment.** No PASTE2-style partial optimal
  transport (spatial proximity + gene expression similarity jointly) and no
  diffeomorphic/non-rigid warp (STalign) was attempted. The trimming step in the ICP
  is a cheap partial concession to non-full overlap, not a substitute for either.

## Honest next step

In priority order, if this were continued:

1. **Fix the coarse-search granularity properly**, rather than relying on multi-start
   to paper over it: a finer angle grid (or a proper multi-resolution/pyramid search)
   so the true optimum's basin is reliably found even for a large rotation, and
   re-verify with the same synthetic-recovery test at several rotation magnitudes
   (5°, 15°, 30°, 60°) to characterize exactly where the current basin boundary is,
   rather than reporting the one 3°-off failure found here.
2. **Use gene expression, not just spatial density, in the correspondence step** —
   the actual core idea PASTE2 adds beyond pure geometric ICP: for each pair of
   nearby cells across sections, weight or gate correspondence by transcriptomic
   similarity (both sections share the identical 313-gene panel, so this is directly
   available with no new data). This would move the method from pure tissue-outline
   alignment toward PASTE2's actual (partial optimal transport) approach.
3. **Get a real z-spacing.** Either from a co-author/collaborator with access to the
   original Janesick et al. cutting log, or by switching to the DLPFC/spatialLIBD
   12-section dataset (`results/landscape_serial_sections.md`, Tier 2), whose z-
   spacing (10 µm / 300 µm) is exactly and publicly documented, to stress-test the
   pipeline on a real, precisely-known z-range even though that dataset is normal
   brain tissue, not cancer.
4. **A true non-rigid step.** A similarity transform cannot fix local tissue
   deformation (tearing, differential stretching between adjacent sections during
   mounting). A thin-plate-spline or diffeomorphic (STalign-style) refinement after
   the current rigid/affine stage would be the natural next layer, evaluated with the
   same before/after metrics already built here so any future non-rigid step is held
   to the same honesty standard as this one.

## Files

- `tmesim/reconstruct/registration_3d.py` — the module (additive, new file).
- `tests/test_registration_3d.py` — 9 tests: existing-ingester sanity on both real
  sections, ICP/Umeyama correctness given a reasonable seed, determinism given a fixed
  seed, seed-robustness, output shape/z-assignment contract, the before/after quality
  check on real data, and three unit-level checks on the metric primitives. All pass;
  full suite (`./.venv/bin/python -m pytest tests/ -q`) is 518 passed, 0 failed.
- `data_cache/xenium_breast_janesick/{rep1,rep2}/` — the real, live-downloaded data
  (~100 MB total).
- `results/registration_3d_breast_poc_overlay.png` — the verification plot referenced
  above.
