# Landscape survey: spatial omics / spatial proteomics platforms as ingest sources

**Purpose.** tme-sim is downstream of tissue characterization: it consumes whatever a spatial
platform already outputs (cell/spot positions + expression/protein counts + optionally
segmentation) and runs drug-transport physics on top. This document surveys what a customer
actually has on disk after a run, for every platform worth writing an ingester for, so the next
ingester can be scoped correctly the first time.

**Method.** Vendor documentation, primary papers, and two independent 2025 cross-platform
benchmarking studies were used as primary sources; every non-obvious claim below carries a
citation. Where a claim could not be pinned to a primary source, it is marked **UNVERIFIED**
rather than stated as fact. Two figures were verified directly against files already sitting in
this repo's `data_cache/` (marked "verified locally" below) rather than against a web source.

Today's date for currency purposes: 2026-07-28.

---

## 1. Answers to the four critical questions (read this first)

### Q1 — Does anyone ship a genuine 3D / serial-section product TODAY, off the shelf?

**No commercial spatial transcriptomics platform ships routine multi-section 3D reconstruction
as a standard analysis output.** What several platforms *do* ship is **optical sectioning within
a single physical section** (a z-stack a few microns deep, matching the section's own thickness),
which is not the same thing as a 3D tissue volume:

- **10x Xenium**: `morphology.ome.tif` is described by 10x as "a 3D Z-stack of the DAPI image,"
  captured at a **0.75 µm z-step** ([10x Genomics, Image Processing Algorithms](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/algorithms-overview/image-processing)),
  but the tissue itself is a single **5 µm (FFPE) or 10 µm (fresh-frozen)** section
  ([10x Genomics tissue-prep guidance](https://kb.10xgenomics.com/s/article/360036305271-What-fresh-frozen-tissue-section-thickness-should-I-place-on-Visium-Slides); [Xenium FF tissue prep guide](https://cdn.10xgenomics.com/image/upload/v1671645454/support-documents/CG000579_Demonstrated_Protocol_XeniumInSituProtocolsFF_TissuePreparationGuide.pdf)).
  That yields on the order of 7–14 optical planes spanning ~5–10 µm of depth — a thin
  optical stack through one section, not a reconstructed tissue volume. Per-transcript
  `z_location` is reported in `transcripts.csv.gz` ([Xenium file format doc](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html)), and the final `cells.csv.gz`/`cell_feature_matrix` is still a flat, 2D per-cell table (one x_centroid, y_centroid per cell — no z_centroid).
- **NanoString/Bruker CosMx SMI**: imaging uses an **~800 nm z-step** during acquisition (cited
  from CosMx-derived secondary literature — [NoButter, bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.11.25.625243v1) — the exact plane count per FOV is **UNVERIFIED**), and the per-transcript flat file carries a
  z-plane index and a `CellComp` subcellular-location tag
  ([AtoMx flat-file comparison](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/flat-file-exports/flat-files-compare.html)). Same story: sub-section optical resolution, not multi-section 3D. Standard CosMx sections are thin FFPE/FF sections (a few µm), consistent with the benchmark below.
- **Vizgen MERSCOPE/MERFISH**: explicitly documented — for a standard **10 µm** section, the
  instrument acquires **7 z-planes at 1.5 µm spacing**, from 1.5 µm to 10.5 µm depth
  ([Vizgen MERSCOPE Instrument User Guide, Rev H](https://genomics.uci.edu/wp-content/uploads/sites/30/91600001_MERSCOPE-Instrument-User-Guide_Rev-H.pdf)). `detected_transcripts.csv` carries a `global_z` plane index per transcript
  ([Vizgen post-processing docs](https://vizgen.github.io/vizgen-postprocessing/output_data_formats/detected_transcripts_format.html)). Again: one section's depth, not a tissue volume.
- **10x Visium / Visium HD, Stereo-seq, Curio Seeker, GeoMx**: no z dimension at all — these are
  flat capture surfaces (spot/bin arrays, bead pucks, ROI collection) applied to a single section.
- **Akoya PhenoCycler-Fusion**: the commercial whole-slide instrument is a **2D widefield**
  system (QPTIFF output, no native z-stack). Genuine volumetric CODEX/PhenoCycler-chemistry
  imaging exists but only as a **research protocol on confocal microscopes with thick (30–50 µm,
  i.e. 8–10× the standard 4–5 µm) cleared tissue sections** — published as *Highly multiplexed 3D
  profiling of cell states and immune niches in human tumors*, Nirmal/Yapp et al., *Nature
  Methods*, 2025 ([PMC12510885](https://pmc.ncbi.nlm.nih.gov/articles/PMC12510885/); [journal](https://www.nature.com/articles/s41592-025-02824-x)). This is not something you order from Akoya as a standard product run — it is a bespoke confocal workflow built on PhenoCycler chemistry.

**Bottom line for Q1: nobody sells "z-resolved tissue volume" as a menu item.** The z-stacks that
do exist (Xenium, CosMx, MERSCOPE) are an artifact of optical sectioning through one thin slide,
useful mainly for segmentation quality / transcript-to-nucleus assignment, not for reconstructing
tissue architecture in the z-direction. tme-sim's current 2D-only domain model is therefore not
missing a routinely-available data dimension — it is behind exactly one field (`z_location`) that
these three platforms already put in their transcript-level files but which the field is
essentially never used to build a true volumetric microenvironment from a single run.

### Q2 — Given that, what's the standard practice for getting real 3D?

**Serial sectioning + computational registration**, done by the research group, not the vendor.
Practice, per the literature surveyed:

- **Section spacing in practice**: reported ranges vary widely by tissue and platform. One
  volumetric-ST synthesis paper cites 16 µm sections at 32 µm center-to-center spacing as a
  common design ([arXiv 2511.14613](https://arxiv.org/pdf/2511.14613)); Open-ST (below) demonstrates alignment across sections spanning up to 50–100 µm; the BGI/STOmics **MOSTA** mouse-embryo atlas used **4–6 serial sagittal
  sections per stage** (17 at the latest stage, E16.5) as its standard depth of serial sampling ([STOmics MOSTA](https://db.cngb.org/stomics/mosta/)). There is no single industry-standard spacing — it is chosen per study, and specific counts/spacings beyond these examples are **UNVERIFIED** as a general rule.
- **Standard registration tools**, confirmed from primary sources:
  - **PASTE / PASTE2** — Fused Gromov-Wasserstein optimal transport alignment of consecutive
    slices (full or partial overlap); outputs a pairwise mapping matrix used to stack slices into
    a 3D reconstruction. PASTE2: Liu et al., *bioRxiv* 2023 / *RECOMB* 2023,
    [PMC9881963](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9881963/), [PubMed 37553263](https://pubmed.ncbi.nlm.nih.gov/37553263/), [GitHub](https://github.com/raphael-group/paste2).
  - **STalign** — converts ST measurements to raster images and registers them via
    diffeomorphic (LDDMM-style) metric mapping; explicitly built to align across sections,
    samples, *and technologies* (e.g. Visium ↔ Xenium), and to model local non-linear tissue
    distortion. (Cited via benchmarking review, [PMC11312151](https://pmc.ncbi.nlm.nih.gov/articles/PMC11312151/).)
  - **GPSA** — a two-layer Gaussian-process model that maps observed coordinates into a common
    coordinate system, then maps that common system to expression; probabilistic, handles
    non-rigid warps. (Same source.)
  - **SpaceFlow** — worth a correction here: SpaceFlow (Ren et al., *Nature Communications*
    2022) is a spatially-regularized deep-learning method for **domain identification / a
    pseudo-spatiotemporal map within one section**, not a multi-slice registration tool. It is
    sometimes cited alongside registration methods in reviews, but does not itself perform
    serial-section alignment; if it appeared on a list of registration tools, that is imprecise.
    Flagging this explicitly since the task named it as an example.
  - **Open-ST** (Schott et al., *Cell* 2024, [DOI](https://www.cell.com/cell/fulltext/S0092-8674(24)00636-6)) is the most complete public
    end-to-end example of true 3D spatial transcriptomics: serial sectioning + their own
    alignment tool (STIM) reconstructs a genuine tissue volume from consecutive Open-ST
    (sequencing-based, Stereo-seq-like) sections. This is a research protocol, not a commercial
    product.
- **Takeaway for tme-sim**: if we ever want a real z-axis, the path is "ingest N consecutive
  single-section outputs from any one of our existing 2D ingesters, then bolt on a registration
  step (PASTE2 or STalign are the best-evidenced choices)" — not "wait for a vendor to ship 3D."
  That is a new pipeline stage, not a new ingester.

### Q3 — Which platforms measure anything relevant to drug transport (vasculature, ECM, hypoxia, efflux transporters, proliferation)?

This splits sharply by **whole-transcriptome vs. targeted-panel** platforms:

- **Whole-transcriptome / very-high-plex platforms capture these genes for free**, no panel
  curation decision required: **10x Visium / Visium HD** (whole transcriptome by design), **Stereo-seq**
  (whole transcriptome), **CosMx Whole Transcriptome (WTx)** prototype (~19,000 genes, first
  shown at whole-transcriptome imaging in 2024 — [NanoString press release](https://nanostring.com/pressrelease/news-news-details-2024-nanostring-launches-the-cosmx-6000-plex-rna-assay-and-showcases-whole-transcriptome-imaging-capability/)), and **Xenium Prime 5K** (5,001 genes, broad
  "pan-tissue and pathways" curation — [10x Genomics](https://www.10xgenomics.com/products/xenium-5k-panel)) come closest to guaranteeing coverage of PECAM1/CD31, VEGFA, HIF1A/CA9, ABCB1/ABCG2, MKI67, and collagen genes without a custom add-on.
- **Targeted panels are hit-or-miss and must be checked gene-by-gene.** Verified **directly from
  the actual panel file already cached in this repo**
  (`data_cache/xenium_gbm/Xenium_GBM_gene_panel.json`, the 10x Xenium Human Brain 254-gene base
  panel + 65-gene GBM add-on, 319 genes total): **PECAM1 present, MKI67 present; VEGFA, CA9,
  HIF1A, ABCB1, ABCG2, COL1A1 all absent.** So even the panel tme-sim already ingests does **not**
  carry a hypoxia or efflux-transporter signal — this is a real, verified gap in the current
  pipeline's molecular inputs, not a hypothetical one.
- **Protein panels (Akoya PhenoCycler)**: the PhenoCode Discovery IO60 panel's "Tissue
  Architecture" category includes **CD31 (PECAM1 protein), Collagen IV, Podoplanin, CD34** —
  i.e. vasculature and basement-membrane/ECM markers are present by default in a stock protein
  panel ([Akoya PhenoCode IO60 brochure](https://www.akoyabio.com/wp-content/uploads/AkoyaBio-PhenoCode-IO60-Panel-Brochure.pdf)); it also carries Ki67 (proliferation). No hypoxia or efflux-transporter antibody was
  found in that panel — **UNVERIFIED whether any stock Akoya panel includes CA9 or ABCB1/ABCG2**,
  though both exist as off-the-shelf antibodies Akoya sells for custom panels.
- General pattern worth stating plainly: **efflux transporters (ABCB1/ABCG2) and hypoxia markers
  (CA9, HIF1A target genes) are not reliably in any platform's default targeted panel** — they
  have to be deliberately added as custom genes/antibodies, or you need a whole-transcriptome
  platform. This is a real constraint on which datasets can feed a drug-transport model without a
  custom panel order.

### Q4 — Is there a converging interchange format? Would one ingester buy multiple platforms?

**Yes, converging around `AnnData`/`.h5ad` at the table level and `SpatialData`/OME-Zarr at the
multi-modal level**, but it is not yet universal enough to replace platform-specific ingest.

- **SpatialData** (Marconato et al., *Nature Methods* 2025, [DOI 10.1038/s41592-024-02212-x](https://www.nature.com/articles/s41592-024-02212-x)) is
  explicitly designed as "an open and universal data framework for spatial omics." Its on-disk
  format is a **Zarr store that extends OME-NGFF**, combining OME-Zarr (images/labels),
  Parquet (points/shapes), and **AnnData-in-Zarr (tables)** into one container with five element
  types: images, labels, points, shapes, tables ([SpatialData docs](https://spatialdata.scverse.org/en/latest/design_doc.html)). The project states an explicit
  goal to converge further with OME-Zarr/NGFF over time.
  Readers already exist for Xenium, Visium/Visium HD, CosMx, MERSCOPE, and CODEX/PhenoCycler
  output into SpatialData — i.e. the *conversion* problem is largely solved by the ecosystem, but
  someone still has to write (or reuse) a platform-specific reader to get raw vendor files into
  that container in the first place.
  - This means: writing tme-sim's ingester to consume **AnnData (.h5ad) with `obsm['spatial']`
    coordinates** as its internal normalized form, rather than any one vendor's raw layout, would
    let us reuse the community's existing Xenium/Visium/CosMx/MERSCOPE/CODEX → SpatialData/AnnData
    readers (scanpy, squidpy, spatialdata-io) as a thin translation shim, instead of hand-rolling
    six parsers. That is the single highest-leverage format decision available.
- Bioconductor's parallel answer is **SpatialExperiment** (R), interoperable with the same
  Python-side AnnData/Zarr tables per the SpatialData design doc's stated goal of R↔Python
  interop (`anndataR`), but this only matters if tme-sim ever needs R interop — **UNVERIFIED**
  as a priority for this project.
- Practical caveat: none of Xenium/Visium/CosMx/MERSCOPE ship AnnData/SpatialData as their
  *native* customer-facing output — they all ship their own bespoke CSV/MEX/H5/QPTIFF layouts
  (documented per-platform below), and SpatialData is what you get *after* running a converter
  (`spatialdata-io`) on top of that native output. So a from-scratch ingester still has to know
  each vendor's real file layout; SpatialData mainly pays off if we standardize tme-sim's
  *internal* representation rather than trying to make vendors emit it directly.

---

## 2. Comparison table

| Platform (vendor) | Current instrument gen. | Molecule | Panel / plex | Whole-transcriptome? | Resolution | Single-cell or spot/bin? | Section thickness & z | Segmentation provided? | Native output container |
|---|---|---|---|---|---|---|---|---|---|
| **Xenium In Situ** (10x Genomics) | Xenium Analyzer, Xenium Prime chemistry | RNA (+optional protein add-on) | v1: 289/377/480-gene pre-designed or custom; **Prime 5K**: 5,001 genes ([10x](https://www.10xgenomics.com/products/xenium-5k-panel)) | No (targeted; 5K is broad but not WT) | Single-molecule, subcellular | Single-cell | 5 µm FFPE / 10 µm FF; optical z-stack 0.75 µm step within that section ([10x](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/algorithms-overview/image-processing)) — not multi-section 3D | Yes — nucleus-expansion or multimodal boundary-stain segmentation, shipped | CSV/Parquet + MEX/HDF5/Zarr + OME-TIFF |
| **Visium** (10x Genomics) | Visium v1/v2, CytAssist | RNA | Whole transcriptome | Yes | 55 µm spots, ~1–10 cells/spot | Spot (not single-cell) | No z; single section | No (spot pooling; `spaceranger segment` optional add-on) | MEX/HDF5 + spatial CSV/JSON/PNG |
| **Visium HD** (10x Genomics) | Visium HD w/ CytAssist | RNA | Whole transcriptome | Yes | 2 µm native bins, 8/16 µm analysis bins, native H&E-segmented cell output since SpaceRanger v4 ([Bin2cell](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11419951/); [10x outputs](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview)) | Bin-based, near-single-cell at 2 µm | No z; single section | Yes as of SpaceRanger v4 (`segmented_outputs/`) | Hierarchical `binned_outputs/{square_002um,008um,016um}` + `segmented_outputs/` + `barcode_mappings.parquet` |
| **CosMx SMI** (NanoString / now Bruker Spatial Biology) | CosMx SMI (5000-plate / 1000-plate); AtoMx SIP for reprocessing | RNA (+protein) | 1,000-plex UCC ([panel](https://nanostring.com/products/cosmx-spatial-molecular-imager/cosmx-rna-assays/cosmx-universal-cell-characterization-rna-panel/)), 6,000-plex, **WTx** ~19,000-gene prototype ([press release](https://nanostring.com/pressrelease/news-news-details-2024-nanostring-launches-the-cosmx-6000-plex-rna-assay-and-showcases-whole-transcriptome-imaging-capability/)) | 6K/WTx effectively whole-transcriptome-like coverage; 1K is targeted | Subcellular, single-molecule | Single-cell | Thin FFPE/FF section; ~800 nm z-step during acquisition ([NoButter, bioRxiv](https://www.biorxiv.org/content/10.1101/2024.11.25.625243v1)), z-plane index kept per transcript | Yes — multi-channel (nuclear/membrane/cytoplasm) stain-based segmentation shipped | AtoMx flat files (`exprMat`, `metadata`, `fov_positions` CSV, optional `polygons`/`transcripts`), or TileDB/Seurat exports |
| **GeoMx DSP** (NanoString / Bruker) — *not requested but relevant, added since shipping* | GeoMx DSP | RNA and/or protein | Up to whole transcriptome (WTA) or curated protein panels | Yes (WTA option) | **ROI-based, NOT single-cell** — user-drawn or AI-segmented regions | Region-of-interest | No z; single section | User-defined ROI, optionally with AOI cell-type masks; not per-cell segmentation | `.dcc` (per-ROI counts) + `.pkc` (probe config) ([GeoMx docs](https://brukerspatialbiology.com/products/geomx-digital-spatial-profiler/geomx-dsp-configuration-files/)) |
| **MERSCOPE / MERFISH** (Vizgen) | MERSCOPE, **MERSCOPE Ultra** (2024, 3.0 cm² imaging area) | RNA (+protein add-on) | 140 / 300 / 500 / up to 1,000-gene panels ([Vizgen](https://spatial.lifestrandsgx.com/vizgen-merscope-ultra/)) | No (targeted, up to 1,000 genes) | Subcellular, single-molecule | Single-cell | 10 µm section imaged as **7 z-planes at 1.5 µm spacing** ([Vizgen Instrument User Guide](https://genomics.uci.edu/wp-content/uploads/sites/30/91600001_MERSCOPE-Instrument-User-Guide_Rev-H.pdf)) — one section's optical stack, not 3D volume | Cellpose-based segmentation shipped by default; polygons provided | `detected_transcripts.csv`, `cell_by_gene.csv`, `cell_metadata.csv`, cell-boundary HDF5/Parquet |
| **PhenoCycler-Fusion (CODEX)** (Akoya Biosciences) | PhenoCycler-Fusion 2.0 | Protein (cyclic IF) | 10–60+ plex, curated PhenoCode Discovery panels ([Akoya](https://www.akoyabio.com/phenocycler/instrument/)) | N/A (protein) | Whole-slide, single-cell / subcellular | Single-cell | 2D widefield by default (QPTIFF); genuine z-stack 3D only via bespoke confocal + 30–50 µm thick-section research protocol, not the stock Fusion workflow ([Nirmal/Yapp, Nat Methods 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12510885/)) | No — segmentation done downstream (QuPath, HALO, Mesmer, etc.) | QPTIFF (multi-channel pyramidal TIFF) |
| **Curio Seeker (Slide-seq)** (Curio Bioscience) | Seeker kit, 8 tiles/run | RNA | Whole transcriptome | Yes | 10 µm beads, ~1–2 cells/bead | Bead (near-single-cell, not segmented) | No z; single section | No (bead pooling) | `.rds`, `.h5ad`, `.mtx` bead-by-gene + bead-location CSV ([Curio Seeker](https://curiobioscience.com/seeker/)) |
| **Curio Trekker (Slide-tags)** (Curio Bioscience, launched Oct 2024) | Trekker kit | RNA | Whole transcriptome (runs on standard scRNA-seq chemistry) | Yes | True single-cell (nuclei individually barcoded in situ, then sequenced) | Single-cell, **no image-based segmentation needed** | No z; single section | N/A — bypasses segmentation entirely by barcoding nuclei directly ([BusinessWire launch](https://www.businesswire.com/news/home/20241008842253/en/Curio-Bioscience-Launches-Trekker-Single-Cell-Spatial-Mapping-Kit)) | Standard single-cell matrix (h5ad/mtx) + per-cell spatial coordinates |
| **Stereo-seq** (BGI / STOmics) | DNBSEQ-T7RS / G400RS + Stereo-seq chip | RNA | Whole transcriptome | Yes | DNB spots ~500 nm pitch; commonly binned (Bin20 ≈ cell-sized, Bin50/100 also used) | Sub-single-cell native resolution, cell-bin (`.cellbin.gef`) available after segmentation | No z; single section (MOSTA atlas used 4–17 serial 2D sections per stage, not a fused 3D volume) ([MOSTA](https://db.cngb.org/stomics/mosta/)) | Optional — `.cellbin.gef` cell-bin output requires a segmentation step (image-based) | `.gem`/`.gem.gz` (text) and `.gef` (HDF5-based) — square-bin or cell-bin variants ([Stereopy docs](https://stereopy.readthedocs.io/en/v1.5.1/Tutorials/IO.html)) |
| **Element Biosciences (AVITI24 spatial / "Direct In Sample Sequencing")** | AVITI24, spatial capability announced for **2026** roadmap | RNA (+ expressed variants) | Not yet a fixed panel product — announced, not GA at time of writing | UNVERIFIED — described as high-plex, not confirmed WT | UNVERIFIED | UNVERIFIED | UNVERIFIED | UNVERIFIED | UNVERIFIED — **not yet shipping to customers as of this survey**; treat as pipeline/announced, not orderable ([Element 2026 roadmap](https://www.elementbiosciences.com/2026-innovation-roadmap)) |
| **Singular Genomics G4X** | G4X Spatial Sequencer | RNA + protein + fH&E, simultaneous | 500-plex RNA + 18-plex protein ([Singular Genomics](https://www.singulargenomics.com/g4x)) | No | Subcellular | Single-cell (on-instrument segmentation + clustering) | UNVERIFIED (no z/section-thickness spec found) | Yes — on-instrument cell segmentation and clustering shipped | "Standard spatial file outputs" per vendor; exact file/container spec **UNVERIFIED** from public sources — would need a vendor NDA/manual to ingest today |
| **Ultivue InSituPlex** (now part of Vizgen, merged Oct 2024) | InSituPlex + STARVUE | Protein | Multiplex IF, exact plex UNVERIFIED | N/A | Single-cell | Single-cell | UNVERIFIED | UNVERIFIED | Migrating under the Vizgen brand ([GenomeWeb](https://www.genengnews.com/topics/omics/vizgen-ultivue-merge-further-consolidating-spatial-biology/)); treat as a Vizgen-family protein product going forward, not a standalone vendor |
| **Lunaphore COMET** | COMET (4-slide, fully automated cyclic IF) | Protein | Up to 40-plex/cycle via sequential IF ([Lunaphore](https://lunaphore.com/products/comet/)) | N/A | Single-cell / subcellular, 20× | 2D, single optical section, DAPI/TRITC/Cy5 channels | 2D only (no z-stack found in public docs) | No — downstream tools (QuPath/HALO/HORIZON) | Background-subtracted OME-TIFF, channels 8-bit ([Lunaphore docs](https://forum.image.sc/t/cropping-wsi-fl-channels-with-python/102047)) |
| **Miltenyi MACSima (MICS)** | MACSima System | RNA + protein (sequential cycles, 3 antibodies/cycle) | Hundreds of targets across cycles | N/A / not WT for RNA mode | Subcellular, widefield | Single-cell after segmentation | UNVERIFIED (image "stacks" mentioned but not confirmed as true z or just multi-cycle) | No — MACS iQ View does segmentation as a downstream analysis step, not shipped raw | Proprietary image stacks + MACS iQ View exports; exact file/container spec **UNVERIFIED** from public sources |

---

## 3. Per-platform detail (file-by-file, for ingester scoping)

### 10x Xenium In Situ (incl. Xenium Prime 5K)

Two documentation sources cover this in detail and largely agree:
[Xenium File Format Documentation](https://cf.10xgenomics.com/supp/xenium/xenium_documentation.html) and
[Understanding Xenium Outputs](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/analysis/xoa-output-understanding-outputs).

Files a customer has after a run (Xenium Onboard Analysis, current generation):

- `experiment.xenium` — JSON manifest (run/sample metadata, links other files).
- `cells.csv.gz` (also `cells.parquet`) — one row per cell: `cell_id`, `x_centroid`,
  `y_centroid` (µm), `transcript_counts`, `control_probe_counts`, `control_codeword_counts`,
  `total_counts`, `cell_area`, `nucleus_area`. **No z_centroid** — this file is flat 2D.
- `transcripts.csv.gz` / `transcripts.parquet` — one row per called transcript: `transcript_id`,
  `cell_id`, `overlaps_nucleus`, `feature_name`, `x_location`, `y_location`, `z_location` (µm),
  `qv` (Phred quality score).
- `cell_boundaries.csv.gz`, `nucleus_boundaries.csv.gz` — polygon vertices (`cell_id`, `vertex_x`,
  `vertex_y`), clockwise, closed.
- `cell_feature_matrix/` — sparse count matrix, available as MEX (`matrix.mtx.gz` +
  `features.tsv.gz` + `barcodes.tsv.gz`), HDF5 (`.h5`), and Zarr.
- `morphology.ome.tif` — full 3D DAPI z-stack (0.75 µm step); `morphology_focus/` — 2D best-focus
  autofocus projection, pyramidal OME-TIFF, 16-bit, 0.2125 µm/px, JPEG-2000 compressed.
- `gene_panel.json` — panel metadata + full target list (this is exactly the file already cached
  at `data_cache/xenium_gbm/Xenium_GBM_gene_panel.json` in this repo).
- `metrics_summary.csv`, `analysis_summary.html`, `analysis/` (clustering/PCA/UMAP), `aux_outputs/`.

This is the file set tme-sim's `tmesim/reconstruct/real_xenium.py` and
`tmesim/data/spatial_ingest.py` already consume — confirms the existing ingester targets the
correct, current file layout.

Panel sizes: pre-designed v1 panels are ~289/377/480-gene ([10x panel page](https://www.10xgenomics.com/products/xenium-v1-panel)); **Xenium Prime 5K**
covers 5,001 genes, customizable +100 ([10x Prime 5K panels](https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/pre-designed-xenium-prime-5k)); a "dual chemistry" approach combining a v1
custom panel with Prime 5K on one section has also been demonstrated in the literature ([PMC12803277](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12803277/)).

Public datasets: `www.10xgenomics.com/datasets` hosts dozens of Xenium runs (breast, lung, skin,
heart, brain cancer, etc.), no login required ([10x datasets](https://www.10xgenomics.com/datasets)); this repo already pulls the GBM one.

### 10x Visium and Visium HD

[Space Ranger output overview](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/output-overview).

Visium v1/v2: `outs/spatial/` holds `tissue_positions.csv` (barcode→pixel/array coords),
`scalefactors_json.json` (px-to-µm/hires/lowres scale factors), `tissue_hires_image.png`,
`tissue_lowres_image.png`. Counts as `filtered_feature_bc_matrix/` (MEX) and `.h5`. 55 µm spot
diameter (verified against this repo's own `tmesim/data/spatial_ingest.py` docstring, which
already documents this correctly).

Visium HD adds a hierarchical structure: `binned_outputs/square_{002,008,016}um/` (each with its
own filtered/raw matrix + spatial folder; 8/16 µm bins also get `analysis/`), plus
`barcode_mappings.parquet` (tracks square↔nucleus↔cell↔bin relationships), `feature_slice.h5`
(fast 2 µm image retrieval), and — new as of SpaceRanger v4 (mid-2025) — `segmented_outputs/`
with native H&E-based cell segmentation, eliminating the need for third-party tools like Bin2cell
for basic cell-level output ([Bin2cell, PMC11419951](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11419951/); [EPFL Visium HD protocol](https://www.epfl.ch/research/facilities/gene-expression-core-facility/wp-content/uploads/2025/09/Guidelines-10XG-Visium-HD-Cytassist-by-GECF-vA.06.pdf)).

### NanoString/Bruker CosMx SMI

Primary paper: He et al., *Nature Biotechnology* 2022, DOI [10.1038/s41587-022-01483-z](https://doi.org/10.1038/s41587-022-01483-z).
File-format detail: [AtoMx flat-file comparison blog](https://nanostring-biostats.github.io/CosMx-Analysis-Scratch-Space/posts/flat-file-exports/flat-files-compare.html).

AtoMx SIP exports data as **flat files** (human-readable CSV), TileDB, or Seurat objects. Only
three flat files are strictly required for reanalysis:

- `exprMat_file.csv` — per-cell counts: `fov`, `cell_ID`, one column per gene target, plus
  `Negative(#)` (negative-probe controls) and `SystemControl(#)`.
- `metadata_file.csv` — per-cell morphology + position: `Area`, `Width`, `Height`,
  `AspectRatio`, `Perimeter`, `Circularity`, `Eccentricity`, `Solidity`, `CenterX/Y` (local +
  global, px), `Area.um²`, `NucArea`, plus a study-wide `cell_id` and FOV-scoped `cell_ID`.
  FOV-level QC stats (`median_RNA`, `median_negprobes`, etc.) also live here.
- `fov_positions_file.csv` — `FOV`, `x_global_px/mm`, `y_global_px/mm`.

Two more flat files are commonly exported but optional: a **transcript coordinates file** (`fov`,
`cell_ID`, `cell`, local/global x/y, **z-plane**, `target`, `CellComp` ∈ {Membrane, Nuclear,
Cytoplasm, None}) and a **polygons file** (per-cell boundary vertices, local + global px).

Panels: 1,000-plex Universal Cell Characterization ([panel page](https://nanostring.com/products/cosmx-spatial-molecular-imager/cosmx-rna-assays/cosmx-universal-cell-characterization-rna-panel/)), 6,000-plex RNA assay, and a
**whole-transcriptome (WTx) prototype at ~19,000 genes** shown in 2024 ([press release](https://nanostring.com/pressrelease/news-news-details-2024-nanostring-launches-the-cosmx-6000-plex-rna-assay-and-showcases-whole-transcriptome-imaging-capability/)).

Public dataset: **CosMx SMI NSCLC FFPE Dataset** — 8 samples / 5 patients, ~800,000 cells, ~260M
transcripts, 960-gene panel, freely downloadable ([Bruker Spatial Biology](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/nsclc-ffpe-dataset/)).

Note on corporate status: **Bruker acquired substantially all of NanoString's assets (including
CosMx, AtoMx, GeoMx, nCounter) for ~$392.6M in an asset deal that closed May 2024**; the business
now operates as **Bruker Spatial Biology** ([Bruker press release](https://brukerspatialbiology.com/pressrelease/news-news-details-2024-nanostring-technologies-to-be-acquired-by-bruker-corporation/)). Documentation/URLs exist under both
`nanostring.com` and `brukerspatialbiology.com` — both currently resolve, but expect the former
to eventually redirect/retire.

### Vizgen MERSCOPE / MERFISH (incl. MERSCOPE Ultra)

[Detected transcripts format doc](https://vizgen.github.io/vizgen-postprocessing/output_data_formats/detected_transcripts_format.html); [MERSCOPE Instrument User Guide](https://genomics.uci.edu/wp-content/uploads/sites/30/91600001_MERSCOPE-Instrument-User-Guide_Rev-H.pdf).

Core output files:

- `detected_transcripts.csv` — one row per called molecule: `barcode_id` (codebook row index),
  `global_x`/`global_y` (µm), **`global_z`** (index into the z-plane stack), `x`/`y` (px, FOV
  frame), `fov`, `gene`, `transcript_id`, optional `EntityID` (post-segmentation cell assignment,
  -1 if unassigned).
- `cell_by_gene.csv` — counts matrix, one row per segmented cell.
- `cell_metadata.csv` — per-cell spatial coordinates + segmented cell area.
- Cell-boundary geometry shipped as HDF5 or Parquet (per Seurat's `ReadVizgen` and squidpy's
  Vizgen tutorial conventions — exact schema **UNVERIFIED** beyond that container choice).

Imaging geometry (primary source): for a **10 µm** section, MERSCOPE acquires **7 z-planes at
1.5 µm spacing**, starting 1.5 µm above the coverslip ([Instrument User Guide, Rev H](https://genomics.uci.edu/wp-content/uploads/sites/30/91600001_MERSCOPE-Instrument-User-Guide_Rev-H.pdf), p. as cited
by search index — direct page number **UNVERIFIED**, but the guide is the primary vendor manual).

Panels: 140/300/500/up to 1,000-gene pre-designed and custom panels (18-bit/21-bit/27-bit
codebooks) ([Vizgen predesigned panels](https://vizgen.com/wp-content/uploads/2023/12/MERSCOPE-Predesigned-Panel-Flyer_Rev-B.pdf); [MERSCOPE Ultra](https://spatial.lifestrandsgx.com/vizgen-merscope-ultra/)). Not whole-transcriptome. **MERSCOPE Ultra** (first
shipments announced Sept 2024) triples imaging area to 3.0 cm²/slide and runs MERFISH 2.0
chemistry ([Vizgen shipment announcement](https://vizgen.com/vizgen-first-worldwide-shipments-merscope-ultra/)).

Public data: **Vizgen Data Release Program** — 26 samples, 10M+ cells, ~6B transcripts, free
download, largest public single-cell spatial genomics dataset per Vizgen's own framing
([data release program](https://vizgen.com/data-release-program/)).

Corporate note: **Vizgen and Ultivue merged (announced Oct 9, 2024)**; Ultivue's InSituPlex assays
and STARVUE platform are migrating under the Vizgen brand, name retained "for now" per Vizgen's
CEO ([GenomeWeb](https://www.genengnews.com/topics/omics/vizgen-ultivue-merge-further-consolidating-spatial-biology/)).

### Akoya PhenoCycler-Fusion (CODEX)

Output is **QPTIFF only** — a multi-channel pyramidal TIFF ([Akoya PhenoCycler instrument page](https://www.akoyabio.com/phenocycler/instrument/)).
Fusion is a whole-slide widefield 2D imager; there is no native per-cell table shipped by the
instrument — segmentation and quantification (QuPath, HALO, Mesmer/deepcell, or Akoya's own
tools) are a separate downstream step the customer runs. Protein plex 10–60+, curated
"PhenoCode Discovery" panels (e.g. IO60 includes CD31/PECAM1, Collagen IV, Podoplanin, CD34 under
"Tissue Architecture," plus Ki67 for proliferation — [panel brochure](https://www.akoyabio.com/wp-content/uploads/AkoyaBio-PhenoCode-IO60-Panel-Brochure.pdf)).

Public datasets: **HuBMAP Consortium data portal** hosts multiple CODEX/PhenoCycler FOV datasets,
combinable into a SpatialFeatureExperiment ([Voyager CODEX vignette](https://pachterlab.github.io/voyager/articles/codex_landing.html)); Akoya/Quanterix have also
released a 2.5-million-cell PhenoCycler-Fusion example dataset ([Quanterix blog](https://www.quanterix.com/blog/analyzing-2-5-million-cells-with-phenocycler-fusion-explore-the-dataset/)).

Genuine 3D CODEX-chemistry imaging exists only as a confocal, thick-section (30–50 µm) research
protocol, not a stock Fusion-instrument product — see Q1 above ([Nirmal/Yapp, Nat Methods 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12510885/)).

### Curio Seeker (Slide-seq) and Curio Trekker (Slide-tags)

**Seeker**: bead-based (10 µm beads, ~1–2 cells/bead), whole transcriptome. Pipeline outputs
bead-by-gene matrices in `.rds`, `.h5ad`, and `.mtx`, plus a spatial-location file and an HTML QC
report ([Curio Seeker product page](https://curiobioscience.com/seeker/)). Legacy Slide-seqV2 convention: raw DGE matrix
(`.txt.gz`) + bead-barcode-location `.csv` per puck ([Broad WARP Slide-seq pipeline](https://broadinstitute.github.io/warp/docs/Pipelines/SlideSeq_Pipeline/README)). No cell
segmentation — one bead is treated as the unit, not a segmented cell. No z.

**Trekker** (launched Oct 2024, based on Broad Institute "Slide-tags" method): fundamentally
different — barcodes tag **nuclei in situ**, which are then dissociated and run through standard
single-cell sequencing; spatial position is recovered computationally per nucleus. This sidesteps
image-based segmentation entirely and gives **true single-cell resolution** without a microscope
in the loop ([BusinessWire launch](https://www.businesswire.com/news/home/20241008842253/en/Curio-Bioscience-Launches-Trekker-Single-Cell-Spatial-Mapping-Kit)). Output is a standard single-cell matrix (h5ad/mtx-like) with a
spatial-coordinate column per cell — closer in shape to what tme-sim's Xenium ingester expects
than to a Visium-style spot file.

### Stereo-seq (BGI / STOmics / Complete Genomics)

DNA-nanoball (DNB) patterned array, ~500 nm nanoball pitch, giving effectively sub-cellular
native resolution before binning ([STOmics](https://en.stomics.tech/why-stereo-seq/index.html)). Output formats:

- `.gem`/`.gem.gz` — plain text, columns `geneID`, `x`, `y`, `MIDCount`, `ExonCount`.
- `.gef` — HDF5-based binary; two flavors: **square-bin** (`.raw.gef` / `.gef` / `.tissue.gef`,
  user-chosen bin size, e.g. Bin20 ≈ one cell, Bin50/Bin100 for lower-resolution analysis) and
  **cell-bin** (`.cellbin.gef`, requires an image-based cell-segmentation step to produce) ([Stereopy IO docs](https://stereopy.readthedocs.io/en/v1.5.1/Tutorials/IO.html)).

Whole transcriptome. No z-dimension per section. Public example: **MOSTA** (Mouse Organogenesis
Spatiotemporal Transcriptomic Atlas), 53 sagittal sections across E9.5–E16.5, freely downloadable
with bin50 and segmented single-cell data per section ([MOSTA download](https://db.cngb.org/stomics/mosta/download/)) — this is the best public
example of BGI's own serial-sectioning practice (4–6 sections/stage, 17 at E16.5), though sections
are analyzed and released as independent 2D slices, not fused into one 3D volume.

### Element Biosciences and Singular Genomics ("Element/Singular")

These are two distinct, unrelated companies, both sequencing-instrument makers entering spatial:

- **Element Biosciences** (AVITI24 sequencer): announced **"Direct In Sample Sequencing"** for
  tissue as part of a **2026 roadmap** — spatial expressed-variant detection + high-plex cell
  typing directly in FFPE/fresh-frozen tissue, up to 20 cm² imageable area per run, <1 day
  tissue-prep ([Element 2026 roadmap](https://www.elementbiosciences.com/2026-innovation-roadmap)). Protein detection is planned for 2027. **As of this survey this
  is an announced capability, not a shipping product with a documented customer file layout** —
  do not scope an ingester against it yet; revisit when GA.
- **Singular Genomics G4X**: this one **is shipping** — described as having processed ~20,000
  FFPE samples and profiled ~4B spatially resolved single cells across "dozens of sample types"
  and 25 technology-access partners ([Singular Genomics G4X](https://www.singulargenomics.com/g4x)). Simultaneous 500-plex RNA + 18-plex
  protein + formalin-H&E at subcellular resolution, 128 samples / 40 cm² per run (~10× competitors'
  throughput per the vendor). On-instrument base-calling, cell segmentation, and clustering are
  provided. **Exact output file names/formats are not published anywhere I could find** — vendor
  states "standard spatial file outputs" for downstream open-source/commercial tools, but the
  concrete file layout is **UNVERIFIED** without a vendor manual or a real sample dataset in hand.
  This is a real, shipping, high-throughput multiomic platform we currently know the least about
  mechanically.

### Ultivue InSituPlex

Antibody-oligo-barcode based multiplex IF with in-situ signal amplification ([overview](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9389457)). As of
Oct 2024, Ultivue is part of Vizgen (see MERSCOPE section above) — treat future Ultivue-derived
data as arriving under Vizgen's product umbrella; exact InSituPlex/STARVUE file format is
**UNVERIFIED** from public sources at the level of detail needed to write an ingester.

### Lunaphore COMET

Fully automated 4-slide sequential IF instrument, up to 40-plex per cycle, 20× objective, DAPI/
TRITC/Cy5 channels ([Lunaphore product page](https://lunaphore.com/products/comet/)). Output: background-subtracted **OME-TIFF**, each
channel downsized to 8-bit, viewable/exportable via the Lunaphore Viewer, compatible with QuPath/
HALO/HORIZON downstream ([image.sc forum thread on COMET OME-TIFF](https://forum.image.sc/t/cropping-wsi-fl-channels-with-python/102047)). 2D only — no evidence in public docs of a
z-stack acquisition mode. No cell segmentation shipped — this is an image-only product;
segmentation and quantification happen entirely downstream.

### Miltenyi MACSima (MICS)

Microfluidic cyclic staining: stain → image → erase, three fluorochrome-conjugated antibodies
per cycle, RNA cycles then protein cycles on the same section, widefield fluorescence, subcellular
resolution ([Miltenyi MICS technology page](https://www.miltenyibiotec.com/IE-en/products/macs-imaging-and-spatial-biology/macsima-platform/mics-technology.html); [Scientific Reports 2022](https://www.nature.com/articles/s41598-022-05841-4)). Analysis (segmentation, clustering) is
done in vendor software **MACS iQ View**, working on the raw "image stacks" the instrument
produces. **Exact raw output file format/container is UNVERIFIED from public sources** — public
documentation describes the workflow and the downstream software, not a customer-facing file
spec. This, like Singular Genomics G4X and Ultivue, would need direct vendor engagement or a
real sample in hand before an ingester could be written with confidence.

---

## 4. Interchange format landscape (detail beyond Q4 answer)

- **SpatialData** (`spatialdata` Python package, scverse ecosystem): Zarr store extending
  OME-NGFF; images/labels as OME-Zarr, points/shapes largely Parquet-backed, tables as
  AnnData-in-Zarr ([Marconato et al., Nat Methods 2025](https://www.nature.com/articles/s41592-024-02212-x); [design doc](https://spatialdata.scverse.org/en/latest/design_doc.html)). Has readers (`spatialdata-io`) for Xenium,
  Visium/Visium HD, CosMx, MERSCOPE, CODEX and others — this is the most-converged point in the
  ecosystem today, but it's a *post-hoc* conversion layer, not what any vendor ships natively.
- **OME-NGFF / OME-Zarr**: cloud-optimized bioimaging format that SpatialData explicitly builds
  on and aims to track ([Zarr paper, Histochem Cell Biol 2023](https://link.springer.com/article/10.1007/s00418-023-02209-1)). Covers the *imaging* side (morphology
  TIFFs, QPTIFF-equivalent) more than the *tabular* side.
  10x's own Xenium output already includes `analysis.zarr.zip` and `cells.zarr.zip`, and CosMx
  offers a TileDB export — so vendors are independently drifting toward Zarr-family containers
  even without adopting SpatialData wholesale.
- **AnnData / `.h5ad`**: the de facto Python single-cell/spatial table format (`obs` = per-cell
  metadata, `obsm['spatial']` = coordinates, `X`/`layers` = counts). Not a vendor-native output
  for any platform surveyed here, but it is what **every** platform's community-maintained reader
  converts to (Seurat's readers convert to Seurat objects; scanpy/squidpy readers convert to
  AnnData) — this is the most broadly-supported target if tme-sim wants one internal schema.
- **SpatialExperiment** (Bioconductor/R): the R-ecosystem analog; relevant only if tme-sim ever
  needs R interop (**UNVERIFIED** as a priority).

**Recommendation for tme-sim specifically**: don't try to get vendors to hand you SpatialData —
none of them do. Instead, standardize tme-sim's *internal* ingest contract on an AnnData-shaped
object (cell/spot table with `obsm['spatial']`, `X` counts, optional `obs['cell_type']`), and
write one thin per-vendor reader for each platform's *native* files into that shape. That gets
90% of the benefit of "one format buys many platforms" without depending on any vendor to change
what they ship.

---

## 5. Ranked recommendation: next 3 ingesters to write

Scored on **(customers who actually have this data) × (how much of tme-sim's model it feeds) ×
(implementation cost)** — lower cost and higher reach both count in favor.

### #1 — CosMx SMI (NanoString/Bruker) flat-file ingester

- **Customers**: CosMx is the #2 imaging-based platform by adoption after Xenium, with a large
  installed base and, critically, a **free public NSCLC FFPE dataset** (~800K cells, 960 genes)
  ready to test against today ([Bruker NSCLC dataset](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/nsclc-ffpe-dataset/)) — no waiting on a customer to hand us data.
- **Model fit**: single-cell resolution, shipped segmentation (nuclear/membrane/cytoplasm-based),
  and a per-transcript z-plane + subcellular-compartment tag — a strict superset of what our
  Xenium ingester already consumes. The file shapes (`exprMat`/`metadata`/`fov_positions` CSVs)
  are structurally close enough to Xenium's `cells.csv.gz`/`cell_feature_matrix` that most of
  `tmesim/reconstruct/real_xenium.py`'s downstream logic (phenotyping, field construction) should
  port with only the parsing layer rewritten.
  Caveat we should not gloss over: CosMx's Universal Cell Characterization panel composition
  relative to our drug-transport gene list (PECAM1/VEGFA/CA9/ABCB1/etc.) has **not** been checked
  gene-by-gene in this survey — verify before committing to it as a drug-transport data source,
  the same way we found gaps in the Xenium brain panel.
- **Cost**: low-medium. File formats are well-documented, a public dataset exists for immediate
  testing, and it reuses most of the existing Xenium reconstruction pipeline.

### #2 — MERSCOPE / MERFISH (Vizgen) ingester

- **Customers**: second-largest imaging-ST installed base; **26-sample, 10M+-cell public data
  release program** available with no login ([Vizgen data release program](https://vizgen.com/data-release-program/)) — again, immediately
  testable without waiting on a partner.
- **Model fit**: single-cell resolution with shipped Cellpose-based segmentation and boundary
  polygons — same shape of information as Xenium (cell table + boundaries + per-transcript
  coordinates with a z-plane index). MERSCOPE's PanCancer Pathways panel is explicitly cancer/
  TME-curated, which is closer to tme-sim's GBM/tumor use case than a generic panel — worth
  verifying its drug-transport gene coverage specifically (not done in this survey) since a
  cancer-pathway curation is more likely than a generic panel to include VEGFA/hypoxia genes.
- **Cost**: low-medium, similar shape to CosMx — CSV-based (`detected_transcripts.csv`,
  `cell_by_gene.csv`, `cell_metadata.csv`), well-documented, no proprietary binary container to
  reverse-engineer.

### #3 — Akoya PhenoCycler-Fusion (CODEX) protein ingester

- **Customers**: the dominant high-plex **spatial proteomics** platform, large installed base,
  public data via HuBMAP and Akoya/Quanterix example releases.
- **Model fit**: this is the highest-leverage pick for the physics model specifically, not just
  for "more data" — it is the one platform surveyed whose **default curated panel already
  includes CD31/PECAM1 (vasculature), Collagen IV/Podoplanin/CD34 (ECM/basement membrane), and
  Ki67 (proliferation)** without any custom add-on ([PhenoCode IO60 brochure](https://www.akoyabio.com/wp-content/uploads/AkoyaBio-PhenoCode-IO60-Panel-Brochure.pdf)) — i.e. it directly
  measures three of the five drug-transport-relevant axes tme-sim cares about (vasculature, ECM,
  proliferation) at the protein level, which is arguably more mechanistically direct than RNA
  proxies for these particular fields (a vessel is defined by CD31 protein, not PECAM1 mRNA).
  It does **not** natively provide hypoxia or efflux-transporter markers in the panels checked.
- **Cost**: higher than #1/#2. Native output is QPTIFF (a multi-channel image), not a per-cell
  table — this ingester has to do image-based cell segmentation (or consume someone else's, e.g.
  a HuBMAP-provided Mesmer/DeepCell segmentation) before it produces the cell/coordinate table
  our `reconstruct` layer expects. This is a materially different, harder ingest path than
  CSV-parsing Xenium/CosMx/MERSCOPE, which is why it's ranked #3 rather than #1 despite arguably
  the best drug-transport-relevant panel content of anything surveyed.

**Explicitly not recommended next, with reasons**: Visium HD (whole-transcriptome coverage is
attractive, but it's spot/bin-pooled, not single-cell, so it would need a genuinely new
reconstruction path rather than reusing the Xenium-style single-cell pipeline — bigger lift for
uncertain payoff given Visium classic is already ingested); Stereo-seq and Curio (compelling
whole-transcriptome coverage and, for Curio Trekker, a genuinely novel true-single-cell modality,
but currently no evidence of a customer/dataset pipeline pulling us toward either — revisit if
that changes); Singular Genomics G4X and Miltenyi MACSima (both promising on paper — G4X for
simultaneous RNA+protein+H&E throughput, MACSima for RNA+protein on one section — but neither has
a publicly documented file format, so step one for either would be reverse-engineering a format
from a vendor-provided sample, not a straightforward ingester write); Element Biosciences
(announced, not GA — nothing to ingest yet).

---

## 6. Sources not otherwise inlined above

- Ren et al., *Systematic benchmarking of high-throughput subcellular spatial transcriptomics
  platforms across human tumors*, Nature Communications 2025, DOI
  [10.1038/s41467-025-64292-3](https://doi.org/10.1038/s41467-025-64292-3), open access at
  [PMC12534522](https://pmc.ncbi.nlm.nih.gov/articles/PMC12534522/). Directly compared Stereo-seq
  v1.3, Visium HD FFPE, CosMx 6K, and Xenium 5K against CODEX protein and scRNA-seq on the same
  serial sections (colon adenocarcinoma, hepatocellular carcinoma, ovarian cancer). Confirmed:
  8 µm as the common cross-platform analysis bin size used for comparison; **5 µm section
  thickness for the FFPE platforms (Visium HD, CosMx, Xenium)**, 10 µm for fresh-frozen
  Stereo-seq and CODEX; no 3D/z-stack analysis performed by any platform in this benchmark.
  Xenium 5K showed the best segmentation accuracy and CODEX-protein concordance; CosMx 6K showed
  elevated negative-control signal.
- *Systematic benchmarking of imaging spatial transcriptomics platforms in FFPE tissues*
  (Vizgen/Broad-affiliated comparison of Xenium, MERSCOPE, CosMx on FFPE tumor TMAs), bioRxiv
  2023 / Nature Communications 2025, [PMC10723440](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10723440/); [Vizgen's own summary](https://vizgen.com/resources/comparison-of-imaging-based-single-cell-resolution-spatial-transcriptomics-profiling-platforms-using-formalin-fixed-paraffin-embedded-tumor-samples/).
  Confirms platform chemistry differences (Xenium: padlock probes + RCA; CosMx: probes + branched
  hybridization; MERSCOPE: direct hybridization + combinatorial barcoding) and that Xenium
  produced higher per-gene transcript counts on matched genes without loss of specificity in
  this comparison.

## 7. What this means concretely for tme-sim's current code

- The Xenium ingester (`tmesim/reconstruct/real_xenium.py`, `tmesim/data/spatial_ingest.py`) is
  targeting the correct, current file layout — verified against both 10x's live documentation and
  the actual cached panel JSON in `data_cache/xenium_gbm/`.
- The molecular gap flagged in Q3 is real and specific: the 319-gene GBM panel already in
  `data_cache/xenium_gbm/Xenium_GBM_gene_panel.json` has **no hypoxia gene (CA9, HIF1A) and no
  efflux transporter (ABCB1, ABCG2)** — any drug-transport claim keyed to hypoxia or efflux
  resistance from this specific dataset would be unsupported by the underlying panel, independent
  of anything about the physics model itself.
- The z-axis tme-sim doesn't currently model is not something any surveyed platform delivers as
  routine 3D tissue architecture — so "add a z dimension" is correctly scoped as a future
  registration-pipeline project (PASTE2/STalign on multiple 2D ingests), not as "wait for a
  richer single-run file format."
