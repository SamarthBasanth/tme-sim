# Closing the open threads in `data_scouting_2026-10-09.md` — 2026-10-09

**Status: VERIFICATION PASS. Nothing here is a result, a test, or a pre-registration.**
No code was changed. **No existing file in `results/` or `archive/` was edited, including
the record this file checks.** Where that record is wrong, the correction is recorded
here and only here. Nothing below is registered in `tmesim/validation/scorecard.py`.

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

*(Closing guard in §7.)*

**What this pass was for.** `results/data_scouting_2026-10-09.md` §7 separates what its
author verified against a primary record from what was reported to them, and says plainly
that nothing in the second category may be quoted onward. This pass closes those threads:
each one is now **verified**, **refuted**, or **recorded as a verified negative with the
exact blocker**. Four threads, in the order the brief gave them.

**Provenance rule, unchanged.** Every number below was retrieved live on 2026-10-09 from a
primary record. Nothing is inferred or reconstructed. Where a conclusion is mine rather
than a source's statement — the platform argument in §1.2, the footprint measurement in
§1.4 — it is labelled as such.

**Instruments that worked, recorded because two of them are not obvious.**

| what | how |
|---|---|
| **PMC full text including all tables** | `eutils efetch.fcgi?db=pmc&id=<numeric PMCID>&retmode=xml`. **This works where `europepmc.org/.../fullTextXML` returns HTTP 500**, which is what blocked §2.5 of the record being checked. |
| PMC supplementary binaries | `europepmc.org/webservices/rest/<PMCID>/supplementaryFiles` (OA articles only). `pmc.ncbi.nlm.nih.gov/articles/instance/<id>/bin/<file>` is behind a proof-of-work challenge and returns an HTML page, not the file. |
| bioRxiv full text | Cloudflare-rate-limited (error 1015); `curl --retry 4 --retry-delay 25 --retry-all-errors` with a browser UA got through. |
| GEO | `acc.cgi?targ=self|gsm|all&form=text&view=brief|full`, and `suppl/filelist.txt` (columns: Archive/File, Name, Time, **Size**, Type). |
| CELLxGENE | `api.cellxgene.cziscience.com/curation/v1/collections/<id>` — note it exposes per-dataset `cell_count` and asset `filesize` but **no per-cell-type counts**. |
| per-cell-type counts in a multi-GB remote H5AD | HTTP `Range` requests into an `h5py` file-like object: **2.36 MB transferred out of 8.13 GB** to read `obs` of core GBmap. |

---

## 1. Thread 1 — the matched-pair ambiguity. **The named source cannot settle it, and that is now a verified fact rather than an open question.**

The record's §2.3 says the load-bearing unverified item is whether a Visium↔Xenium matched
pair is *the same physical section* or *adjacent sections of the same block*, and names the
bioRxiv methods preprint as the thing that would settle it. Four findings.

### 1.1 The preprint cannot settle it, because the matched sections are not in it. **VERIFIED NEGATIVE.**

bioRxiv **2025.04.04.647096** (Khan, Wang, … Petti). Two versions, v1 2025-04-04 and
**v2 2025-04-22**; the API reports no published version. It also has **PMID 40235981 /
PMC11996482**, which the record did not have — but the PMC copy is marked
*"The publisher of this article does not allow downloading of the full text in XML form"*
and serves the abstract only. Full text was retrieved from bioRxiv directly.

- The word **"adjacent"** occurs 9 times in the full text, **every one of them inside the
  authors' `Cluster Adjacency` metric or a figure caption about it**. `serial`,
  `consecutive`, `same section`, `same block`: **0 occurrences each.** There is no
  sectioning-relationship statement in the preprint at all.
- More decisively, the preprint's Xenium cohort is a different set of patients. Verbatim:
  *"We generated Xenium data for two patient cohorts … the primary cohort contained five
  samples from four different patients; the secondary cohort contained four samples from
  one patient."* Nine Xenium samples, five patients.
- **Supplementary Table S1** (`media-15.docx` from the v2 supplementary material), titled
  *"Patients and samples for Visium and Xenium analysis"*, lists them explicitly:

| assay | patients in the preprint |
|---|---|
| Visium | TWBK-B172, B176, B178, B185, B186, **TWBK-B189 (aka ST073021)**, TWKI-WU1220, TWKI-WU1221, TWKI-WU1227 |
| Xenium | GBM024 (4 sections), GBM025, GBM026 (2), GBM030, GBM034 |

  **The two sets are disjoint. The published preprint contains no matched Visium↔Xenium
  pair of any kind.** The nine Xenium samples are GEO **GSE316917**.

- The six matched pairs come from **GSE335275**, whose GEO title ends `[xenium_cohort3]`
  and which was **submitted 2026-06-12, fourteen months after preprint v2**. The preprint
  knows only cohorts 1 and 2 (Fig. S7 is captioned *"Analysis of Cohort 1 and Cohort 2
  samples"*). **Cohort 3 postdates the public manuscript**, so no version of it describes
  how those six sections relate to the Visium sections.
- GEO is silent too: all six GSE335275 samples carry
  `!Sample_extract_protocol_ch1 = n/a`, and the series-level design pastes the same
  Primary/Secondary cohort text as GSE316917 without a cohort-3 preparation paragraph.
- The Zenodo deposit for this study (**10.5281/zenodo.20709378**, dated 2026-06-15) does
  carry cohort 3 — it ships `Combined_Xenium_both_Cohorts_Aggregated_GBM.rds` and
  `paired_xenium_object_list_with_broadcelltype_annotations.rds`, the authors' own use of
  the word *paired* — but its description contains no methods text. Its
  `HE_overview.pdf` (5.9 MB) is **Visium H&E only**, all 22 sections on their fiducial
  frames; there are no Xenium images in it.

**Conclusion: the source the record named cannot answer the question, and no other public
record does either.** The only remaining route is the corresponding author
(`apetti@mgh.harvard.edu`, in the GEO contact block).

### 1.2 "The same physical section" is ruled out on platform grounds. **My inference, from two verbatim protocols.**

- Xenium, verbatim (GSE335275 / GSE316917 design, and the preprint methods): *"Formalin-fixed,
  paraffin-embedded (FFPE) tumor samples were sectioned into **5-um sections** using a
  microtome **onto Xenium slides** according to Protocol CG000578."*
- Visium, verbatim (GSE318562 extract protocol): *"The next nine samples were preserved in
  FFPE, and cDNA libraries were generated according to the **Visium Spatial Gene Expression
  for FFPE User Guide Rev. E**."*

A section is mounted on one slide. A Visium FFPE capture slide and a Xenium slide are
different consumables, so **a Visium↔Xenium pair is necessarily two different physical
sections.** The record's binary ("same section *or* adjacent sections") therefore collapses
to a single open quantity: **how far apart the two sections are.** No paper states it.

### 1.3 The chemistry worry the record raised is retired. **VERIFIED NEGATIVE — all six pairs are FFPE↔FFPE.**

The record warns that *"the GSE318562 cohort is mixed fresh-frozen/FFPE while the Xenium
sections are all FFPE, so the chemistry of each specific pair needs checking too."*
Checked. GSE318562's extract protocol says verbatim: *"The first 13 samples were
fresh-frozen and OCT-embedded … The next nine samples were preserved in FFPE."* The nine
FFPE samples are **GSM9497624–GSM9497632**, carrying
`!Sample_description = Spatial transcriptomics of glioblastoma using 10x Visium (V1) FFPE`.
All six Visium partners of the matched pairs — ST073021-02 (GSM9497624), ST073021-04
(…625), WU1220-01 (…626), WU1220-03 (…628), WU1221-Core2 (…629), WU1221-Edge1 (…630) —
are inside that block. **Every matched pair is FFPE Visium against FFPE Xenium. The
mixed-chemistry confound does not touch any of the six.**

Also settled in passing: the **B189 = ST073021** identity now has three independent
confirmations rather than one — the GEO Visium characteristic `batch: B189` on
GSM9497624/5, the GEO Xenium characteristic `patient id: B189` on GSM9809293/4, and the
preprint's own Table S1 row **"TWBK-B189 (aka ST073021)"**. The Zenodo filenames
`ST073021.B189.2.outs.zip` / `ST073021.B189.4.outs.zip` are a fourth.

### 1.4 A measurement, since no document answers: the two assays do not cover the same tissue extent, and the worst-agreeing donor is the three-way-matched one.

Public records being silent, I measured the only thing that is cheap and decisive-adjacent:
the **tissue footprint** of each half of each pair. Visium in-tissue spot positions
(`tissue_positions_list.csv` + `scalefactors_json.json`, converted to µm by
`55 / spot_diameter_fullres`) against Xenium cell centroids (`cells.parquet`, already in
µm). Extents are PCA-aligned so they do not depend on slide rotation; Visium area is
in-tissue spots × the 100 µm hex-lattice cell area (8,660 µm²); Xenium area is
100 µm grid-cell occupancy of the centroids. **Total download for all six pairs: ~70 MB.**

| pair (Visium ↔ Xenium) | Visium spots | Visium extent µm | Visium area mm² | Xenium cells | Xenium extent µm | Xenium area mm² | **area ×** | long × | short × |
|---|---:|---|---:|---:|---|---:|---:|---:|---:|
| **B189** ST073021-02 ↔ ST073021_2 | 1,787 | 6507 × 4386 | 15.48 | 52,234 | 10706 × 8652 | 55.32 | **3.57** | 1.65 | 1.97 |
| **B189** ST073021-04 ↔ ST073021_4 | 2,850 | 5787 × 5529 | 24.68 | 185,662 | 11218 × 8813 | 50.18 | **2.03** | 1.94 | 1.59 |
| WU1220-01 ↔ WU1220_2 (see §1.5) | 2,915 | 7019 × 5004 | 25.24 | 56,669 | 8984 × 5901 | 34.01 | **1.35** | 1.28 | 1.18 |
| *(WU1220-02, the alternative partner)* | 3,194 | 6391 × 5451 | 27.66 | 56,669 | 8984 × 5901 | 34.01 | *1.23* | *1.41* | *1.08* |
| WU1220-03 ↔ WU1220_3 | 3,995 | 7077 × 6232 | 34.60 | 296,613 | 7552 × 7401 | 38.77 | **1.12** | 1.07 | 1.19 |
| WU1221-Core2 ↔ WU1221_Core2 | 1,919 | 4570 × 4130 | 16.62 | 103,688 | 6123 × 6359 | 22.27 | **1.34** | 1.34 | 1.54 |
| WU1221-Edge1 ↔ WU1221_Edge | 2,045 | 4922 × 3467 | 17.71 | 20,859 | 6222 × 4597 | 18.06 | **1.02** | 1.26 | 1.33 |

Three readings, in descending confidence:

1. **The Xenium half is always the larger field** — area ratio 1.02–3.57, never below 1.
   Partly expected: the Visium `in_tissue` mask is an algorithmic tissue call inside a
   6.5 × 6.5 mm capture window, while the Xenium ROI is operator-chosen on a much larger
   slide (*"regions of interest (ROIs) were initially selected using the on-instrument
   software"*, verbatim). **Consequence for the proposed test: the two halves of a pair do
   not cover the same tissue. A composition comparison pooled over each whole section is a
   sub-region against a superset region, and the difference will contain real spatial
   heterogeneity on top of unmixing error.** The comparison has to be restricted to the
   co-registered overlap, or the confound is built in. This is a second, independent
   confound from the adjacency question and nobody had named it.
2. **WU1220 and WU1221 look like the same piece of tissue; B189 does not.** Four of the
   five WU pairs sit at 1.02–1.35× in area and 1.07–1.54× in each linear dimension —
   what you would expect from two sections of one block with a conservative Visium mask.
   **Both B189 pairs are outliers at 2.03× and 3.57× area, 1.59–1.97× linear.** Either the
   B189 Visium sections were trimmed to a much smaller piece than the Xenium sections, or
   they are not the same piece of tissue. **B189 is precisely the donor the record
   identifies as the uniquely valuable three-way match**, so the one design element it
   ranks highest is the one whose two halves agree worst on extent.
3. **Segmented cell density varies ~8× across the six Xenium sections and ~4× within each
   donor** — 944 /mm² (ST073021_2) to 7,651 /mm² (WU1220_3); within-donor B189 944 vs
   3,700, WU1220 1,666 vs 7,651, WU1221 4,656 vs 1,155. Any Xenium cell-fraction is a
   ratio whose denominator moves by 4× between two sections of the same tumour. That is
   not fatal to a *fraction* comparison, but it should be looked at before the run, not
   after.

**What this measurement is not.** It is a footprint-size comparison, not a registration.
Equal area does not prove adjacency; unequal area does not prove non-adjacency. I did not
align the two coordinate systems or compare shapes. **Adjacency remains UNVERIFIED for all
six pairs.** What is now verified is that the extents differ, by how much, and which donor
is the outlier.

### 1.5 An error in GSE335275 that would silently mis-pair one of the six. **New, and load-bearing.**

**GSM9809295 contradicts itself.** Its GEO title is
`!Sample_title = HumanBrain-Tumor-WU1220-01`, but **all nine of its deposited
supplementary files are named `GSM9809295_WU1220_2_*`** (`…_cells.parquet.gz`,
`…_cell_feature_matrix.zarr.zip`, `…_morphology.ome.tif.gz`, and the rest), listed as such
in both `filelist.txt` and the sample's own `!Sample_supplementary_file` lines. The title
says section 01; the files say section 2.

Every other cohort-3 sample is internally consistent (ST073021-02 ↔ `ST073021_2`;
WU1220-3 ↔ `WU1220_3`; WU1221-Core2 ↔ `WU1221_Core2`; WU1221-Edge1 ↔ `WU1221_Edge`).

This decides which Visium sample GSM9809295 pairs with: **GSM9497626 `TWKI-WU1220-1`** (if
the title is right) or **GSM9497627 `TWKI-WU1220-2`** (if the filenames are right). The
record's §2.3 table asserts the former without noticing the conflict. The footprint
measurement does **not** resolve it — the Xenium field is 1.35× the WU1220-1 Visium area
and 1.23× the WU1220-2 Visium area, both plausible. **One of the six pairs is of unknown
identity until the authors say which label is correct.** Treat the matched set as **five
confirmed pairs plus one ambiguous**, not six.

### 1.6 What thread 1 costs the proposed test

The §2.3 design is still the most valuable item in the record being checked, and it is
still the only one that can change a verdict. But its cost sheet is now explicit:

- it **cannot** be written up as a same-section comparison — ruled out on platform grounds;
- it **cannot** be written up as an adjacent-section comparison either — no source says so,
  and the measurement is equivocal for one donor and adverse for the other;
- it has a **region-mismatch confound** (§1.4 reading 1) that is independent of adjacency
  and must be handled by registration-and-restriction, not by pooling;
- it is **five pairs plus one ambiguous**, not six (§1.5);
- the **three-way B189 design is the weakest-supported part of it**, not the strongest.

The honest framing for a pre-registration is: *two sections of the same FFPE block, spatial
separation unstated by any source, non-identical tissue extents measured at 1.0–3.6× in
area.* That is a weaker instrument than "six matched pairs" sounds, and it should be
declared that way **before** the Xenium half is downloaded, which remains possible because
the halves are separate accessions.

---

## 2. Thread 2 — do the microdialysis patients overlap the Xenium patients? **No public record says. The honest answer is: cohort-level pairing only, pending one email.**

### 2.1 What is now verified

| | |
|---|---|
| **Same trial** | Confirmed, both sides, verbatim. Habashy: *"four patients who underwent intraoperative LIPU/MB … in the NCT04528680 clinical trial."* GSE339006 treatment protocol: *"Patients underwent intraoperative low-intensity pulsed ultrasound with microbubbles (LIPU/MB) using the SonoCloud-9 device as part of clinical trial NCT04528680."* |
| **Same institution** | Habashy: *"Patients were enrolled from a single tertiary care center from 2023 to 2024"*, IRB of **Northwestern University**. GSE339006 `!Series_contact_institute = Northwestern University`. |
| **Same people** | **Karl J Habashy is himself a listed contributor on GSE339006**, alongside Adam M Sonabend and Roger Stupp. |
| **Xenium patient IDs** | All five, not the three the record names: **01-01-110, 01-01-121, 01-01-124, 01-01-133, 01-01-135**. The sibling Visium series writes them `NU01-01-124`, confirming `NU` = Northwestern and that these are trial subject IDs. |
| **Microdialysis patient IDs** | **There are none.** Habashy's full text identifies patients only as *"Patient 1"* to *"Patient 4"*, in Tables 1–3. No trial subject ID appears anywhere in the article. |

### 2.2 Why it cannot be closed from public records

1. **No companion paper exists.** The GSE339006 family is titled *"Ultrasound-Mediated BBB
   Opening Promotes Microglial Reprogramming and Enhanced Endothelial Interactions in the
   Human Peritumoral Brain"*. PubMed and EuropePMC (including the preprint index) return
   **0 hits** for that title and for `"microglial reprogramming" AND ultrasound`. The data
   are public; the manuscript is not. **The GEO record is the only authority, and it says
   nothing about microdialysis or about drug dosing beyond fluorescein.**
2. **The three sibling series use three incompatible ID schemes**, so the deposit cannot
   even be joined to itself:

| series | assay | n samples | patient IDs |
|---|---|---:|---|
| GSE339006 | Xenium Prime 5K | 10 | `01-01-110`, `-121`, `-124`, `-133`, `-135` |
| **GSE338999** | **Visium** | **4** | `NU01-01-124`, `NU01-01-126` |
| GSE339010 | scRNA-seq | 14 | `P004033`, `P104`, `P105`, `P107`, `P108`, `P109`, `P111` |
| GSE338996 | Multiome ATAC+RNA | 8 | `P8`, `P10` |

   The `P1xx` and `P8/P10` schemes do not map onto the `01-01-1xx` scheme by any rule
   visible in the deposit.
3. **The two sub-studies of NCT04528680 used different drugs.** Habashy's microdialysis
   arm: *"intravenous infusion of **carboplatin** (over 30 minutes, at an area under the
   curve [AUC] of 3.5) and **fluorescein** (500 mg intravenous bolus)"*, catheters
   implanted after resection, sampling 2–24 h. The biopsy arm described in the group's
   earlier paper (Arrieta et al., *JCI Insight* 2024, PMID 39724015, GEO GSE208074):
   *"we infused fluorescein 500 mg IV and initiated a 45-minute IV infusion of
   **nab-paclitaxel** chemotherapy (Abraxane)"*, with biopsies at 4–15 min and ~45 min, and
   **no mention of microdialysis anywhere in that paper**. GSE339006's biopsies were taken
   *"approximately 45 minutes after sonication"* — the biopsy-arm timing — and its protocol
   names only microbubbles and fluorescein. The two arms are compatible within one surgery
   but are not the same procedure, and nothing states that any one patient did both.
4. Full-text searches of EuropePMC for the literal subject IDs return only noise (the
   indexer splits on hyphens), so that is **not** a usable negative and is not claimed as one.
5. Habashy's Supplementary Table 1 (*"Individual patient characteristics"*) could not be
   retrieved: `NIHMS2090842-supplement-1.docx` sits behind the PMC proof-of-work challenge,
   and this is an NIHPA author manuscript so EuropePMC serves no supplementary package.
   **It is unlikely to help anyway** — the paper's own convention is Patient 1–4, and
   GSE339006 carries no age or sex to cross-check against (`genotype: N/A`, nothing else).

### 2.3 What each resolution would support, stated so the claim cannot drift

- **If they overlap (per-patient):** GSE339006 would be reconstructable tissue of the exact
  compartment in which that patient's own free-interstitial `k` was measured, under the same
  perturbation. That is the unprecedented thing the record hopes for.
- **If they do not (cohort-level, which is all that is currently supported):** the pairing is
  *same trial, same institution, same compartment, same perturbation, different patients*.
  The tissue then informs a **structural prior** — vessel/endothelial geometry,
  `bbb_penetrance`, the peritumoral compartment the Portnow 2009 anchor was measured in —
  and the Habashy `k` remains a **population-level** loss rate applied to a reconstruction
  of someone else's tissue. That is still useful and still better than anything held, but
  it is a cohort anchor and must be written as one.
- **The engine-facing claim does not change between the two cases.** Microdialysis is a
  time course at a point either way; the spatial prediction is unscored either way.

**Recommended wording until settled: "same trial and same compartment; patient overlap
unknown."** Do not write "the same patients". One email to Habashy or Sonabend settles it,
and Habashy is an author on both halves.

### 2.4 One thing the record missed, and it is good news

**GSE338999 is a Visium series on the same peritumoral brain**, listed in the record under
"not inspected". It carries **patient `NU01-01-124`**, who is also in the Xenium series
(`01-01-124`), with paired sonicated / non-sonicated sections on both platforms. **That is a
matched Visium↔Xenium pair in human peritumoral brain, same patient, same perturbation** —
a second, independent instance of the §2.3 design, in the compartment the TMZ partition
anchor came from. It is 4 samples and, like its siblings, ships processed matrices only
(*"Raw sequencing files are not included … due to patient privacy considerations"*, verbatim).
It deserves a size check and a look before anyone budgets for the §2.3 test.

---

## 3. Thread 3 — extractability of the top-ranked anchors. **Three of four are better than the record assumed; one is figure-only.**

Standard applied: the record's own, from its §4 rejection of Cribaro 2021 — the numbers must
be in a table or in text with a usable dispersion, not in a figure, and SEM is not a
population spread.

### 3.1 Habashy 2025 (PMID 40495421, PMC12310118) — **FULLY EXTRACTABLE, per patient, mean ± SD.** The record's §2.0 warning can be lifted.

Citation verified from the JATS front matter: *Clin Cancer Res* 2025 Aug 14; **31(16):
3562–3570**, DOI 10.1158/1078-0432.CCR-25-0080. Reference 18 verified as
**Portnow J, Badie B, Chen M, Liu A, Blanchard S, Synold TW**, *Clin Cancer Res* 2009;
15(22):7092–8, **PMID 19861433** — exactly the anchor the record said it was, and
*"Published data on the peritumoral brain concentration of TMZ for **7 patients** were used
for comparison"*.

Tables 1, 2 and 3 are **machine-readable in the PMC XML, with per-patient rows and
`mean ± sd`** (the paper writes "sd" explicitly in the Table 1 and Table 2 summary rows).

**Table 2 — carboplatin, the row that matters for `alpha_identifiability.md` §4(a):**

| | plasma t½ h | plasma AUC µM·h | **non-son brain t½ h** | non-son brain AUC | **son brain t½ h** | son brain AUC | brain:plasma AUC % non-son | son | son:non-son ratio |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Patient 1 | 2.2 | 175.3 | **4.0** | 24.7 | **2.9** | 83.6 | 14.1 | 47.7 | 3.4 |
| Patient 2 | 4.7 | 118.0 | **7.4** | 9.9 | **5.3** | 37.5 | 8.4 | 31.8 | 3.8 |
| Patient 3 | 2.9 | 90.1 | **5.8** | 25.1 | **4.9** | 70.4 | 27.9 | 78.1 | 2.8 |
| Patient 4 | 3.4 | 88.3 | **3.2** | 12.8 | **6.1** | 30.5 | 14.5 | 34.5 | 2.4 |
| mean ± sd | 3.3 ± 1.1 | 117.9 ± 40.6 | **5.1 ± 1.9** | 18.1 ± 7.9 | **4.8 ± 1.4** | 55.5 ± 25.6 | 16.2 ± 8.3 | 48.0 ± 21.2 | **3.1 ± 0.6** |

Table 3 is the same shape for fluorescein (non-son brain t½ 9.0 / 30.0 / 6.8 / 8.4 h,
mean 13.6 ± 11.0; son 7.7 ± 3.0; ratio 3.1 ± 1.7). Table 1 gives per-patient C<sub>max</sub>
in plasma / non-sonicated / sonicated brain for both compounds.

Four things that must travel with these numbers:

1. **The "3.1-fold" is the sonicated:non-sonicated *brain-to-plasma AUC* ratio**
   (3.1 ± 0.6 carboplatin, 3.1 ± 1.7 fluorescein, P = 0.03 for both), not an influx rate.
   The record's phrase "influx up 3.1-fold" is close but imprecise.
2. **Use t½, not the reported clearance.** The paper's *"carboplatin's clearance from
   non-sonicated brain averaged 10.4 L/hour"* is an apparent dose/AUC clearance in **L/h**;
   it is not a tissue first-order loss rate and it is not in any table. The usable quantity
   is `k = ln2 / t½`, per patient, from Table 2 / Table 3.
3. **The terminal brain slope is not input-limited, which is what makes it usable.** Brain
   t½ exceeds plasma t½ for both compounds (carboplatin 5.1 vs 3.3 h; fluorescein 13.6 vs
   9.6 h), so the brain decay is not simply tracking plasma washout. Had it been the other
   way round, the fitted `k` would have been the plasma rate wearing a tissue label.
4. **Method of the t½:** *"Half-lives (t1/2) were calculated from the elimination rate
   constant derived from the last six measured concentrations as described on Portnow et
   al."* — i.e. the same estimator as the anchor already held, which is convenient and also
   means the two are not independent in method.

Also verified: 6 patients selected, **2 excluded for catheter problems, 4 analysed (8
catheters)**; microdialysate samples every 2 h for 24 h with the first at 4 h (averaging
2–4 h, plotted at 3 h); *"In the brain, microdialysate collections sampled **unbound (free)**
drug levels"*; fluorescein plasma corrected for 85% protein binding **from DrugBank, not
measured** (*"we have not performed an equilibrium dialysis"*) — so the fluorescein
brain:plasma ratios carry a literature correction factor, while carboplatin does not.

### 3.2 Johnson/Sanai 2026 (PMID 41206763, PMC13070491) — **FULLY EXTRACTABLE per patient, both compounds — and the record points at the wrong table.**

CC-BY. Supplementary retrieved via EuropePMC (`noaf257_supplementary_data.zip`, 5.8 MB →
`RiboEve_Supplementary_Tables_FINAL_20251217.xlsx`, 9 sheets).

- **Correction to the record.** It says *"Per-patient values are reported by the strand to
  be in Supplementary Table S3."* **S3 is not per-patient** — it is median / maximum /
  minimum by dose group. **S4 is per-patient but is plasma only.**
- **The per-patient tumour data is in Supplementary Table 7**, *"Individual subject
  characteristics, data availability, pharmacokinetic, and pharmacodynamic results"*, **29
  subjects × 40 columns**, including `Ribociclib.Unbound.enhancing`,
  `Ribociclib.Unbound.non-enhancing`, `Ribociclib.Unbound.CSF`, `Ribociclib.Unbound.plasma`
  and the matching six everolimus columns. Units µmol/L (cross-checked: patient 100-0003's
  `Ribociclib.Unbound.non-enhancing` = 1.77 = the 1,770 nM maximum in S3's 400 mg column).
- **S3 verified against the abstract:** unbound ribociclib in non-enhancing tumour, median
  (max, min) = **170 (1770, 65) nM at 400 mg** and **693 (2345, 68) nM at 600 mg**. The
  900 mg column is flagged in the sheet itself as *"obtained from our previously published
  study (Tien A, Clin Cancer Res 2019) for comparison"* — it is **not** new data.

**The one thing that changes the proposed test.** Counting S7 directly:

| column | numeric | **BLQ** | NA |
|---|---:|---:|---:|
| `Ribociclib.Unbound.non-enhancing` | **23** | 0 | 6 |
| `Ribociclib.Unbound.enhancing` | **24** | 0 | 5 |
| `Everolimus.Unbound.non-enhancing` | **0** | **23** | 6 |
| `Everolimus.Unbound.enhancing` | **0** | **24** | 5 |
| `Everolimus.Unbound.plasma` | 9 | 15 | 5 |

**Everolimus unbound tumour concentration is below the limit of quantitation in every single
patient, in both regions.** So the record's *"two drugs, the same patients … separated by
≥3 orders of magnitude"* is a **censored, one-sided** comparison, not two measured numbers.
Stated exactly: 23 patients have a measured ribociclib unbound non-enhancing concentration
spanning **0.065–2.345 µmol/L (36-fold)**, and in the same 23 patients everolimus is
**< 0.1 nM**. The separation is therefore **≥ 650-fold at the lowest patient** (2.8 orders),
**≥ 1,700-fold at the 400 mg median**, **≥ 6,930-fold at the 600 mg median**. "≥3 orders"
is right at the median and wrong at the minimum.

That still makes a sharp falsification test — a solver must land inside a 36-fold measured
per-patient ribociclib spread *and* predict everolimus below 0.1 nM — but it is a
**one-sided bound on one of the two compounds** and has to be pre-registered as such.

**Scope warning on this file.** Supplementary Tables 1 and 7 also carry `PFS.months`,
`OS.months` and `Vital.status`. Those columns are a patient-outcome dataset and are **out of
scope by the standing prohibition** (`data_scouting_2026-10-09.md` §6). Only the PK columns
are usable. Anyone who ingests this workbook should drop the outcome columns at load time
rather than leave them sitting next to the PK columns.

### 3.3 Schiff 2026 (PMID 41099363, PMC12979031) — **PARTIALLY EXTRACTABLE: medians and ranges in text, per-patient figure-only. The record's unit warning is verified.**

NIHPA author manuscript, no OA supplementary. Four in-text tables — patient characteristics,
best response, and two adverse-event tables. **There is no PK table.** The PK is running text:

- *"Median unbound drug concentrations in non-enhancing and enhancing tumors were **0.188
  (range, 0.135–1.133)** and **0.209 (0.117–0.731) μmol/L (or μmol/kg)**."*
- *"the respective median Kp of **0.23 (range, 0.17–0.49)** and **0.39 (range, 0.19–0.74)**,
  and the respective median **Kp,uu of 0.92 (range, 0.38–1.52)** and **0.98 (range,
  0.42–1.77)**."*
- *"median unbound fraction in plasma [Fu, plasma], **0.031**"*; Fu,tumor **0.083** and
  **0.076**.
- **n = 7** patients for non-enhancing tumour, **n = 5** for enhancing. Per-patient values
  exist only as plotted symbols: *"The symbols (○▲×) represent data from ind[ividual
  patients]"* (Figure 2A–C). Digitisable, not tabulated.

**The record's flag is verified verbatim: the paper writes `μmol/L (or μmol/kg)` and never
commits to a volume-vs-mass basis.** Resolve it with the authors before use; it is the
convention-not-travelling failure the project already catalogues.

**A refinement the record missed.** It presents 0.92 / 0.98 as if they were point values
"at the very top of the 0.013–0.178 range" of the held partition anchors. They are
**medians of 7 and 5 patients spanning 0.38–1.52 and 0.42–1.77** — a ~4.7-fold spread that
still sits entirely above 0.178, but a free partition of "~1" with a lower bound of 0.38 is
a different object from a free partition of 0.92. Carry the range.

Also verified: 66 subjects enrolled; a single 60 mg pamiparib dose on the day of surgery
with blood and tumour sampled 2–8 h post-dose; and the record's claim that no TMZ tissue
concentration is reported — confirmed, TMZ appears only as the 20 mg/day co-dose.

### 3.4 Evans 2008 (PMID 19138031, PMC2629387) — **FULLY EXTRACTABLE. The record's UNVERIFIED flag was a fetch failure, not a property of the paper.**

The record could not reach the full text because EuropePMC returned HTTP 500. The
`efetch?db=pmc` route returns it, tables included. **Five tables, all machine-readable:**

| table | contents |
|---|---|
| **1** | *Approximate Conversion of EF5 Binding to Tissue pO₂* — the calibration ladder: physiological 1–3 %CRB = 100–20 mmHg; mild hypoxia 3–10 %CRB = 20–4 mmHg; moderate 10–30 %CRB = 4–0.75; severe 30–100 %CRB = 0.75–0; anoxia 100 %CRB = 0. **This is what converts the binding axis to pO₂, and it is tabulated.** |
| **2** | *Summary of Tumor Grade and Vessel Findings* — per patient A–F: histology, grade, slides, **vessels/mm²** (244, 70, 16, 127, 24, 11), **vessel area µm² ± SD** (416 ± 878 … 3385 ± 5929), and vessel counts by oxygenation class. **SD, not SEM.** |
| **3** | *Summary of Parameters for the Lines Describing Average EF5 Binding per Band* — **per patient, per oxygenation class, the intercept and slope of the binding-vs-distance line**, in %CRB and %CRB per µm × 1000, for bands 1–5 and 2–5. |
| **4** | *Quantitative Data Supporting Fig. 5* — **counts of vessels with positive vs negative slope**, per section, per oxygenation class (18 vs 13, 6 vs 6, 16 vs 10 in totals). |
| **5** | Calibration correction-factor ranges across the six patients (optics 2.8×, AUC 1.9×, CRB 7.2×). |

The spatial axis is defined verbatim: *"EF5 binding information (percentage cube reference
binding, CRB) was collected as a series of **20 10-μm-wide bands** between adjacent
contours."* The dose is verified as *"given to patients intravenously (**21 mg/kg**)
approximately **24 h** preceding surgery"*, and the six tumours as *"a hemangiopericytoma, a
WHO Grade 2 and WHO Grade 3 glial brain tumor"* plus three glioblastomas (D, E, F, F
recurrent) — exactly as the record describes.

**The one real limit, and it is new.** Table 3 tabulates only **bands 1–5 and 2–5**, i.e.
the first **0–50 µm** from the vessel. The full 0–200 µm profile exists only in the figures.
So the extractable numeric gradient covers **the first 50 µm**, which is still the
perivascular scale the engine cares about, but it is not the whole published profile.

**Conclusion for thread 3's ranking: Evans 2008 moves from "extractability UNVERIFIED" to
"fully tabulated, SD-bearing, free, and the only distance-resolved human in-tumour dataset
in the whole record."** Table 4's sign counts in particular turn *"gradients averaged to near
zero"* from a sentence into a countable statistic.

---

## 4. Thread 4 — the reported-to-me numbers. **Two refuted, the rest verified.**

### 4.1 §2.9, Siletti — verified except one label list. **REFUTED: `fibroblast` is not in the Vascular object.**

CELLxGENE collection `283d65eb-dd53-496d-adb7-7570c7caa443` resolves to *"Human Brain Cell
Atlas v1.0"*, DOI 10.1126/science.add7046, 138 datasets, every one `disease: normal`,
`Homo sapiens`, snRNA-seq (`suspension_type: nucleus`, 10x 3' v3, 4 donors).

| claim | verdict | measured |
|---|---|---|
| `Supercluster: Vascular` = 9,932 cells / 47 MB | **VERIFIED** | `cell_count` **9,932**; H5AD **49,520,426 B = 47.23 MiB** |
| labels incl. `fibroblast` | **REFUTED** | **three** labels only: `endothelial cell`, `pericyte`, `vascular associated smooth muscle cell`. Counts 5,165 / 3,693 / 1,074. |
| `Fibroblast` supercluster 9,156 | **VERIFIED** | 9,156 cells, a **separate** 49,676,931 B file |
| `All non-neuronal cells` 888,263 / 4.3 GB | **VERIFIED** | 888,263; **4,666,925,886 B = 4.35 GiB**; 13 labels, including all four |
| normal brain, not tumour | **VERIFIED** | `normal` (PATO:0000461) across all 138 datasets |

**Practical consequence:** vasculature *and* fibroblasts cannot be had from one 47 MiB file.
It is two ~47 MiB downloads, or the 4.35 GiB non-neuronal object. The record implies one.

### 4.2 §2.9, GBmap — the load-bearing caveat is **VERIFIED**; the vascular count is **REFUTED**.

Collection `999f2a15-3d7e-440b-96ae-2c806799c08c`. Paper identified as PMID **40312969** /
PMC **12526130**, *Neuro-Oncology* 27(9):2281–2295, 2025.

- **Core GBmap 338,564 cells / 110 donors / 16 studies / 7.6 GB — VERIFIED, all four**
  (H5AD 8,127,057,404 B = 7.57 GiB; `donor_id` array has 110 unique entries).
- **The non-independence caveat — VERIFIED from the paper's own Data Availability**, quoted
  verbatim: *"The published datasets that were included in the core GBmap can be accessed
  under GEO/EGA/SRA accession numbers: GSE103224, **GSE131928**, GSE138794 …, GSE148842,
  GSE139448, …, **GSE84465**, …, GSE163108, and GSE163120. Datasets included in the extended
  GBmap can be accessed under GEO/links: …, **GSE182109**, …"*
  **GSE84465 is Darmanis — the reference this project already uses — and it is inside core
  GBmap.** A "GBmap vs Darmanis" comparison is a superset against one of its own subsets.
  The record said this was structurally sound but unverified; it is now verified.
- **Vascular share — REFUTED.** The record says *"the vascular share is ~1% of cells, i.e.
  roughly 3,400 endothelial + mural combined."* Read directly out of `obs` of the core H5AD
  by HTTP range request (2.36 MB of 8.13 GB transferred):
  **endothelial cell 673 + mural cell 1,418 = 2,091**, which is the authors' own
  `annotation_level_2 = Vascular` total exactly, and **0.618%** of 338,564.
  **3,400 is 338,564 × the paper's rounded "1%"** — it is a round-trip through a rounded
  percentage, not a count. Use **2,091**, or **673 endothelial**.
- Two things that make GBmap worse as a vessel reference than the record already says:
  the entire core endothelial compartment is **673 cells across 110 donors** (~6 per donor),
  and the paper disclaims the proportions outright — *"the proportions of cell types may not
  accurately reflect the reality of the tumor."*
- Minor, but it will bite a pipeline: **the CELLxGENE collection's own DOI is the 2022
  bioRxiv preprint** `10.1101/2022.08.27.505439`, not the 2025 Neuro-Oncology DOI the record
  gives. Cite the 2025 paper explicitly.
- Verified in passing: endothelial subtypes are real and tiny — `Tip-like` **175**,
  `Endo capilar` **271**, `Endo arterial` **227** (the misspelling is in the data).

**This does not change the §2.9 recommendation; it sharpens it.** Wälchli's 49,999
GBM-specific endothelial cells for 118 MB remains ~74× the entire core-GBmap endothelial
compartment, and Siletti's 5,165 normal-brain endothelial cells for 47 MiB beat GBmap's 673
at 1/160th the download.

### 4.3 §2.9, Wälchli's CELLxGENE collection — **VERIFIED private and empty.**

Collection ID found from the paper (PMC11324530), not guessed:
`c95ca269-68a7-47c5-82db-da227f31b598` → `"visibility": "PRIVATE"`, `"datasets": []`,
`published_at: null`, created 2024-05-01. The record's "irrelevant, because GEO and UCSC both
serve it" stands. There is no CELLxGENE route to Wälchli.

### 4.4 §2.10 #11, the "25–40× contradiction" — **the numbers are right; the contradiction is not one.**

Both papers retrieved in full. Every number the record flagged as unverified checks out:

| number | verdict | note |
|---|---|---|
| Salavati K: PM(CRC/SB) **7.3 ± 0.9 × 10⁻¹⁵**, PM(OC) **6.1 ± 1.1**, pancreatic **4.1 ± 1.3** m²/(Pa·s) | VERIFIED | **dispersion is SE, labelled once and then dropped** |
| normal peritoneum **2.8 ± 0.7 × 10⁻¹³** | VERIFIED | but the paper says *"one order of magnitude higher"*, never 40×. The ~40× is a derived ratio. |
| *"The variation of K within a single tumor ranged **up to a factor of 4**"* | VERIFIED verbatim | the Discussion restates intra-*type* variation as 11–72%, a different grouping |
| inverse correlation with collagen fibre density | VERIFIED | R² = 0.79, p < 0.0001 pooled; regression in Excel, multiple samples per patient treated as independent |
| **no brain or glioma tissue** | **VERIFIED** | Table 1 is pancreatic, peritoneal/ovarian/small-bowel metastases and normal peritoneum. "brain"/"glioma" occur **only in the reference list**. |
| Elmghirbi TIFP **6.0 ± 3.7 mmHg, range 1.6–14.9** | VERIFIED | **SD** (*"All statistics are reported as mean ± standard deviation"*), and it is the **control group, n = 13**, not all 36 animals |
| IFP transition length **1.0–2.5 mm, mean 1.66 ± 0.48 mm** | VERIFIED verbatim | the paper's own term is *"the length of the pressure gradient"*, symbol `l`; **"boundary layer" is not the paper's word**; the WIN needle advanced *"in steps of 0.5 mm"*, so 1.66 ± 0.48 mm is quantised on a 0.5 mm grid — an 0.48 mm SD is one-step scatter; **`n` for this figure is not stated** |

**On the conversion.** Both papers do report a K, both define it as `v = −K∇p`, and both use
μ_IF = 3.5 × 10⁻³ Pa·s, so the comparison is convention-clean.
1 mm²/(mmHg·s) = 7.500616 × 10⁻⁹ m²/(Pa·s), so Elmghirbi's **(2.3 ± 3.1) × 10⁻⁵
mm²/(mmHg·s)** = **1.73 × 10⁻¹³ m²/(Pa·s)**, giving **23.6× to 42.1×** against Salavati's
4.1–7.3 × 10⁻¹⁵ band. **The arithmetic behind "25–40×" is correct.** But:

- **It compares a measurement to a calibrated model output.** Elmghirbi's MRI-K comes from
  Kozeny–Carman with the specific surface area `S` **fitted by regression against the
  WIN-TIFP of the same 13 control rats**, and the authors state the direction of the bias
  themselves: *"will underestimate S, and thus **bias high** the estimated hydraulic
  conductivity."*
- **The gap is inside each study's own spread.** Elmghirbi's control range is
  2.1 × 10⁻¹⁴ to 9.0 × 10⁻¹³ m²/(Pa·s), a **43× span across 13 rats with SD > mean**.
  Salavati's tumour samples span 9.5 × 10⁻¹⁶ to 1.6 × 10⁻¹⁴, **17×**, and in SE.
- **The denominator does the work.** Against Salavati's own highest tumour mean
  (ovarian metastasis from small bowel, 6.3 × 10⁻¹⁴) the gap is **2.7×**; against normal
  peritoneum (2.8 × 10⁻¹³) Elmghirbi's rat brain tumour is **lower**.
- **Species and organ are confounded and cannot be separated**, because Salavati measured no
  brain tissue at all.

**Recommendation: do not pre-register this as a contradiction.** Two distributions that each
span 20–40× internally, whose central values differ by 24–42×, across different species and
organs, one of them a high-biased model estimate, are not in demonstrated conflict. The
defensible finding is the negative: **no published human brain-tumour hydraulic conductivity
exists in this pair.** Elmghirbi cites Netti 2000 (*Cancer Res* 60:2497–503) for a U87 glioma
K but does not print its value — that is the lead if a glioma K is wanted, and it is
unverified here.

The **1.66 ± 0.48 mm gradient length** survives as the more useful of the two and keeps the
record's framing: a shape constraint independent of absolute magnitude. Carry the 0.5 mm
quantisation and the missing `n` with it.

### 4.5 Byte-unit convention, applying to the whole record being checked

Every size in `data_scouting_2026-10-09.md` written as "MB"/"GB" matches on a **binary**
reading. 47 MB = 47.23 MiB; 4.3 GB = 4.35 GiB; 7.6 GB = 7.57 GiB; 118 MB = 123,811,242 B =
118.07 MiB. The numbers are right; the unit labels are MiB/GiB.

---

## 5. Errors found in `data_scouting_2026-10-09.md`

Recorded here, not corrected there. Four substantive, three cosmetic.

**Substantive:**

1. **§2.2 — "GSE316917 … → 6 donors" is wrong; it is 5.** The nine sections carry
   `patient id:` GBM024 ×4, GBM025, GBM026 ×2, GBM030, GBM034. The preprint's own
   description (5 samples / 4 patients + 4 samples / 1 patient) and Table S1 agree: five
   patients. The record's own sample list implies 5; the count is a slip.
2. **§2.3 — the WU1220 pair may be mis-identified** (§1.5). GSM9809295's title and its nine
   filenames disagree about the section number. Five confirmed pairs plus one ambiguous.
3. **§2.9 — "roughly 3,400 endothelial + mural combined" for GBmap is 2,091** (§4.2). The
   3,400 is the paper's rounded 1% multiplied back out.
4. **§2.9 — the Siletti Vascular object does not carry a `fibroblast` label** (§4.1).
   Fibroblasts are a separate supercluster.
5. **§2.10 — "Per-patient values … in Supplementary Table S3"** for Johnson 2026 is wrong
   (§3.2). S3 is median/max/min by dose group; the per-patient tumour values are in S7.

**Cosmetic, but they would be noticed by anyone checking:**

6. **§2.1 — "BioProject PRJNA1405549"** is attached to the wrong series. GSE318562's is
   **PRJNA1372229**; **PRJNA1405549 belongs to GSE316917**; GSE335275's is PRJNA1477838.
7. **§2.1 — "No PubMed ID is linked in the GEO record"** is true of GEO, but the preprint
   does have **PMID 40235981 / PMC11996482**. The record's conclusion (treat as a preprint,
   not peer-reviewed) is unaffected and correct — the PMC record is abstract-only.
8. **§2.0 — "influx up 3.1-fold"** is imprecisely stated; the 3.1 is the
   sonicated:non-sonicated brain-to-plasma AUC ratio (§3.1).

**Checked and found correct**, since the record's method statements were themselves under
audit: `suppl/filelist.txt` does carry a Size column (I mis-parsed it before re-checking);
`GSE318562_RAW.tar` = 706,416,640 B, `GSE335275_RAW.tar` = 59,217,838,080 B,
`GSE316917_RAW.tar` = 96,007,546,880 B, `GSE339006_RAW.tar` = 42,638,120,960 B — all four as
published; the 22/6/9/13 sample counts; the B189↔ST073021 identity; Habashy's citation,
trial, design, compartment and reference 18 = Portnow 2009.

One labelling caution, not an error in the record: the Zenodo `HE_overview.pdf` shows **23**
Visium panels including a `B185-04` and a `B176-02-02` that **are not in GEO** (GEO has
B185-01/02/03 and B176-02-04), and a `WU1200-3` typo. Do not join Zenodo filenames to GEO
samples by label alone.

---

## 6. Does the ranking change?

**Yes, in three places.**

**Evans 2008 (#5) moves up, and should probably be acted on first.** Its only recorded
weakness was that extractability was unverified; it is now verified as **five tables, SD-bearing,
per-patient, with a tabulated %CRB→pO₂ conversion and a tabulated per-band slope/intercept
for 0–50 µm from the vessel**. It is free, it needs no download pipeline, it is human
in-tumour tissue, and §3.1 of the record being checked establishes that **nothing else like
it has been published since**. It is also the only item in either file that is adverse to a
claim the engine currently makes. Cheapest falsification available.

**The §2.3 matched-pair test (#3) drops — not out, but down.** It is still the only item
that can move a verdict, and it still has the rare pre-registrable property. But it is now
five pairs plus one ambiguous; it cannot be described as same-section or (on the evidence)
as adjacent-section; it carries a region-mismatch confound that needs registration before
any pooling; and the three-way B189 design the record calls uniquely valuable is the one
whose halves agree worst on extent (2.0× and 3.6× in area). **A pre-registration written
today would have to declare all four of those, and after declaring them the test is weaker
than "six matched pairs on three donors" implies.** It is still worth running. It is no
longer obviously worth running *first*.

**§2.0 splits into a strong half and a weak half.** The **Habashy `k`** is confirmed as
everything the record hoped: per-patient, tabulated, mean ± SD, free-interstitial, human
brain, with the terminal brain slope slower than plasma so it is a genuine tissue retention
rate. That half is ready to use today and is the single highest-value number in either file.
The **tissue half is demoted to a cohort-level pairing** until the authors confirm patient
overlap — and the record's framing of GSE339006 as *"reconstructable tissue of the exact
compartment the `k` was measured in"* is only supportable at cohort level right now. The
partial compensation is **GSE338999**, the uninspected Visium sibling, which gives a
within-patient Visium↔Xenium pair (`01-01-124`) in peritumoral brain and should be costed.

**Unchanged:** §2.1 GSE318562 is still the best single find and still changes no verdict.
§2.9's recommendation is unchanged and better supported — Wälchli first, Siletti second,
GBmap for malignant states only and never as an independent reference, now verified rather
than inferred. §2.4, §2.6, §2.7, §2.8 untouched by this pass.

**Dropped:** the §2.10 #11 "25–40× contradiction" should **not** be pre-registered (§4.4).
What survives from that pair is the 1.66 ± 0.48 mm shape constraint and a verified negative:
no published human brain-tumour hydraulic conductivity.

**Still open after this pass, named so the gap stays visible:**
adjacency of the six pairs (authors only); which WU1220 section GSM9809295 is (authors only);
whether Habashy's 4 patients are among GSE339006's 5 (authors only); Habashy Supplementary
Tables 1–2 and Schiff's supplementary (PMC proof-of-work); Evans' 50–200 µm band profile
(figure-only); GSE338999's size and format; and Netti 2000's U87 glioma K.

The record's own structural warning is unchanged and, if anything, reinforced: everything
cheap in both files improves the simulator's **inputs**, and the two items that produce a
**verification result** — Evans 2008 and the per-compound magnitude check — are the ones
that needed this pass to become usable.

---

## 7. Closing calibration guard

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

Unchanged from the opening guard. No code, data loader, scorecard entry or existing dated
record was touched by this pass. The only file this pass created or modified is this one.
Downloads were made to `/tmp` only; nothing was written to `data_cache/`.

**Nothing in this file is registered, pre-registered, or scored.**
