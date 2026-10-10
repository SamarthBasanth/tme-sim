# PRE-REGISTRATION — matched-pair, own-reference deconvolution test
**Written 2026-10-10, BEFORE any agreement number was computed.**
Frozen copy; `results/matched_pair_deconvolution_2026-10-10.md` quotes this file
verbatim and carries its SHA-256. Nothing below may be edited after the first
number in §4 of that file exists.

---

## 0. Calibration guard, start of pass

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

## 1. The question

`results/revalidation.md` §4 retracted a cross-platform composition
concordance: Visium-NNLS-deconvolution against Xenium single-cell fractions read
Pearson **0.997** on the 4-vector and **0.945** on the non-tumour subvector, and
after the classifier arity fix reads **0.980** and **0.016**. §4c concludes the
question is OPEN and names the suspected cause:

> *"Visium NNLS deconvolution against the Darmanis atlas has its own, different
> tumour-inflating tendency (82.66% TUMOR, 1.04% STROMAL) which the old Xenium
> rule happened to lean toward."*

That is a **reference-mismatch** hypothesis: the Darmanis reference (GSE84465,
4 donors, none of them the scored tumour) is wrong for the scored tissue. This
test attacks it directly. Donor **B189** has, in three separate GEO accessions,
its own scRNA-seq, its own Visium sections and its own Xenium sections. Build the
signature from B189's own cells, deconvolve B189's own spots, score against
B189's own single cells. If the §4c attribution is right, agreement improves. If
it does not improve, the attribution is wrong and the collapse is unmixing or
platform — a worse position than the project currently occupies.

## 2. The two confounds, registered before anything is pooled

### 2.1 Section separation — UNVERIFIED and unverifiable from public records

Taken as established by `results/data_scouting_threads_closed_2026-10-09.md` §1,
not re-derived here:

- A Visium↔Xenium pair is **necessarily two different physical sections** (§1.2:
  a 5 µm section is mounted on one slide and the two assays need different
  slides). "Same physical section" is ruled out on platform grounds.
- **How far apart those two sections are is stated by no source.** The bioRxiv
  methods preprint (2025.04.04.647096) covers a disjoint patient set; the matched
  deposit GSE335275 postdates preprint v2 by fourteen months and carries
  `extract_protocol = n/a` on all six samples; the Zenodo deposit has no methods
  text.
- All six pairs are **FFPE↔FFPE** (§1.3), so the mixed-chemistry worry does not
  apply.

**This confound is NOT controlled by anything below. It is registered as an
uncontrolled term and every result in this file inherits it.** Section
separation can only inflate disagreement, never deflate it, so it cannot
manufacture a pass — but it can manufacture a failure, and a failure in C1 must
be reported with that caveat attached.

### 2.2 Region mismatch — the Xenium half is always the larger field

Measured in the same verification pass (§1.4): the Xenium field is 1.02× to
3.57× the Visium field by area, every pair, never below 1. **The two halves of a
pair do not cover the same tissue.** Pooling composition over whole sections
compares a sub-region against a superset region, so real spatial heterogeneity
enters the unmixing error.

**Three options were considered. Two are rejected on recorded grounds.**

| option | verdict | why |
|---|---|---|
| **(a) restrict to the geometric overlap of the two footprints** | **REJECTED — unaffordable and, if afforded, unverifiable** | Co-registration needs a shared physical referent. Visium ships an H&E `tissue_hires_image.png`; the Xenium morphology image is `*_morphology.ome.tif.gz` at **3.98–11.76 GB per section** (GEO `filelist.txt`, measured 2026-10-10). The two B189 sections alone are **10.89 GB + 11.76 GB = 22.65 GB** against **2.3 GiB** of free disk on this host — short by ~10×. Separately: the sections are different physical cut planes (§2.1), so even a perfect H&E↔DAPI registration aligns two different pieces of tissue and its residual cannot be bounded without a third reference. |
| **(b) stratify by region** | **REJECTED — reintroduces free parameters** | Region labels in a common frame need the same registration as (a). An intrinsic conditioning variable (hypoxia score, distance-to-vessel) is itself estimated differently on the two platforms, so the stratification variable would carry the confound it is meant to remove. |
| **(c) report the result as a function of how much of the field is shared** | **CHOSEN** | Implementable from the files that are affordable, with no registration and no new tunable parameter. Specified in §2.3. |

### 2.3 The chosen control: a measured region-mismatch band, in the units of the outcome

The Xenium half resolves single cells over its whole field, so the confound can
be **measured rather than argued about**. For each pair:

1. PCA-align the Xenium cell centroids (same procedure as §1.4, so the areas
   below are an independent reproduction of that table — a check on the record
   this test is built on).
2. `A_vis` = (in-tissue Visium spots) × 8660.25 µm² (the 100 µm hex-lattice cell
   area). `A_xen` = (occupied 100 µm grid cells of the Xenium centroids) × 10000 µm².
3. Draw random axis-aligned sub-windows of the PCA-aligned Xenium field with the
   **aspect ratio of the PCA-aligned Visium extent** and area `A_vis`. Accept a
   draw when its *occupied* tissue area is within **±10%** of `A_vis` and it
   holds **≥ 200 cells**. Target **2000 accepted** draws, cap **200 000**
   attempts. Seed `numpy.random.default_rng(20261010)`.

Each accepted sub-window gives a 4-vector composition. That ensemble **is** the
region-mismatch confound, measured, in composition units, at this pair's area
ratio. From it:

| statistic | definition |
|---|---|
| `L1_full` | L1 distance, Visium-deconvolved composition vs **whole-field** Xenium composition. The naive pooled comparison the verification pass warned against. |
| `L1_best` | **minimum** L1 distance, Visium composition vs any accepted sub-window. The most charitable possible placement of the Visium footprint inside the Xenium field. |
| `L1_wf_50`, `L1_wf_95` | median and 95th percentile of L1 between a random accepted sub-window and the whole field. How much composition moves for region reasons alone, at this area ratio. |

`mean_abs_diff = L1 / 4`, in percentage points, so every number is directly
comparable to the 4.44% / 12.21% already published in `results/revalidation.md`.

**`L1_best` is selection-biased downward** (it is a minimum over 2000 draws).
That bias is why it is used **only** to declare FAILURE and **never** to declare
a pass. If even the best-placed window cannot be matched, the error is real.
This asymmetry is the whole point of the construction and is not negotiable
after the fact.

## 3. Criteria, thresholds and predicted directions

### 3.1 Thresholds, taken from already-published numbers

Both thresholds are arithmetic from `results/revalidation.md` §4a/§4b and contain
**no information from the data downloaded for this test**:

- Retracted-era agreement: mean abs. difference **4.44 pp** → **τ_retracted = L1 0.1778**.
  (Xenium 73.80/7.92/8.54/9.74 vs Visium 82.66/1.04/6.53/9.77 → |Δ| 8.86/6.88/2.01/0.03, mean 4.445.)
- Current published agreement: mean abs. difference **12.21 pp** → **τ_current = L1 0.4884**.
  (Xenium 58.24/11.78/19.48/10.50 vs the same Visium vector → |Δ| 24.42/10.74/12.95/0.73, mean 12.21.)

### 3.2 The criteria

| id | criterion | what each outcome means |
|---|---|---|
| **C1** | `L1_best` (own-reference) vs `τ_current = 0.4884` | **`L1_best > τ_current` ⇒ FAILURE, and the strongest available negative.** The matched-donor, own-reference, most-charitably-placed comparison would be worse than the project's current unmatched, foreign-reference comparison. Removing the reference-mismatch confound would have made agreement worse. |
| **C2** | `L1_best` (own-reference) vs `τ_retracted = 0.1778` | `L1_best ≤ τ_retracted` ⇒ the own-reference deconvolution reaches retracted-era agreement **without** the retracted claim's donor- and reference-mismatch. Not a pass by itself (C2 uses the downward-biased statistic); it is a *necessary* condition for any positive reading. |
| **C3** | `L1_full` vs `L1_wf_95` | `L1_full ≤ L1_wf_95` ⇒ the pooled comparison sits inside the measured region band: the deconvolution agrees with the Xenium truth about as well as a random same-size piece of the same tissue does. **PASS, but weak, and the band width must be quoted with it** — at an area ratio of 3.57 the band is wide by construction. `L1_full > L1_wf_95` while `L1_best ≤ τ_current` ⇒ **UNINFORMATIVE**: region mismatch and unmixing error cannot be separated at this area ratio, and no number may be claimed. |
| **C4** | C1–C3 recomputed on the renormalised **non-tumour 3-vector** (STROMAL, IMMUNE, ENDOTHELIAL) | This is the subvector whose Pearson collapsed 0.945 → 0.016, so it gets its own criterion instead of being averaged into the 4-vector. A 4-vector pass with a 3-vector failure is reported as a failure. |
| **C5** | Spearman ρ of the Visium 4-vector against the whole-field Xenium 4-vector, and the ρ distribution over sub-windows | **REPORTED, NOT A GATE.** With K=4 the achievable ρ values are coarse. Declared non-gating now so a favourable ρ cannot be promoted later. |
| **C6** | **The headline contrast.** `L1_best(own B189 reference)` vs `L1_best(Darmanis reference)`, same pairs, same windows, same everything else | **Predicted direction: own-reference is CLOSER, `L1_best(own) < L1_best(Darmanis)`.** If the prediction holds, the §4c reference-mismatch attribution is supported. **If it fails — if the patient's own cells do not improve agreement — the §4c attribution is wrong**, reference mismatch is eliminated as the explanation, the collapse is unmixing or platform, and the project's position is worse than it is today. |
| **C7** | **The unmixing-only arm.** Xenium cells aggregated into Visium-geometry pseudo-spots (hex lattice, 100 µm pitch, cells within 27.5 µm of a lattice point, ≥5 cells per pseudo-spot), deconvolved with the **same** B189 signature restricted to the 542 Xenium panel genes, scored against each pseudo-spot's known composition | Region mismatch is identically zero (truth and mixture are the same cells) and platform is held fixed. **Small C7 error with large C1 error ⇒ the failure is platform/region, not unmixing. Large C7 error ⇒ the unmixing itself is broken, a stronger negative than C1.** Declared asymmetry: C7's signature has only 542 genes against the Visium arm's whole-transcriptome intersection, so a *large* C7 error is partly attributable to the gene restriction and a *small* one is strong evidence. |

### 3.3 Declared direction for the reference-composition change

The brief records a prior prediction that swapping the vascular reference drives
the vascular fraction **down**, because Darmanis's single small `Vascular` column
absorbs pericyte signal with nowhere else to put it. This test swaps the entire
reference, so the direction is declared now:

> **Swapping Darmanis → B189's own cells will move the Visium-deconvolved
> ENDOTHELIAL fraction DOWN.** Darmanis's `Vascular` label is the original
> authors' mixed mural+endothelial compartment; the B189 reference's ENDOTHELIAL
> class is gated by `PECAM1` alone through `reconstruct.real_xenium.phenotype_cells`,
> which is narrower, and pericytes will fall to that function's ungated TUMOR
> default. Recorded before computing.

## 4. Method, fixed

**Nothing below has a tunable knob, and no quantity is chosen after seeing an
agreement number.** Every component is the repo's shipped default, called
unmodified; `tmesim/**` is not edited.

1. **Labels, both modalities, identical procedure.**
   `reconstruct.real_xenium.phenotype_cells(..., threshold=1.0,
   aggregation='matched_null_z', pool_size=25)` — the shipped defaults — with the
   marker map resolved by `marker_panels.select_marker_map` on the **542-gene
   GSE335275 panel**, applied to the Xenium cells and to the B189 scRNA cells
   **restricted to those same 542 genes**. Same classifier, same markers, same
   gene space on both sides; the only difference between the two sides is the
   platform. (`matched_null_z` is analytic and RNG-free, so this is
   deterministic.)
2. **Signature.** `reconstruct.deconvolve.build_signatures(shared_genes,
   reference=(expr_df, labels_df, gene_names))` — the repo's own function, with
   only the reference tuple swapped. `shared_genes` = Visium symbols present in
   the reference. Labels are emitted as the Darmanis vocabulary
   (`Neoplastic` / `Immune cell` / `Vascular` / `Astrocyte`) so the repo's
   `map_label` collapse is used verbatim.
3. **Deconvolution.** `reconstruct.deconvolve.deconvolve_spot` per in-tissue
   spot, unmodified NNLS. Section composition = the mean of the per-spot
   fraction vectors.
4. **Units.** Visium spot coordinates to µm by `55 / spot_diameter_fullres` from
   `scalefactors_json.json`; Xenium centroids are already µm.
5. **Reported either way.** Ungated-to-TUMOR fraction on both sides, cells and
   spots used, genes shared, and the reproduced §1.4 area table.

## 5. Pairs: used, and excluded with reasons

| id | Visium | Xenium | donor | area × (§1.4) | arm |
|---|---|---|---|---|---|
| **P1** | GSM9497624 `TWBK-ST073021-2` | GSM9809293 `ST073021_2` | B189 | 3.57 | **own-reference (headline)** + Darmanis |
| **P2** | GSM9497625 `TWBK-ST073021-4` | GSM9809294 `ST073021_4` | B189 | 2.03 | **own-reference (headline)** + Darmanis |
| P3 | GSM9497628 `TWKI-WU1220-3` | GSM9809296 `WU1220_3` | WU1220 | 1.12 | foreign-reference only |
| P4 | GSM9497629 `TWKI-WU1221-Core2` | GSM9809298 `WU1221_Core2` | WU1221 | 1.34 | foreign-reference only |
| P5 | GSM9497630 `TWKI-WU1221-Edge1` | GSM9809297 `WU1221_Edge` | WU1221 | 1.02 | foreign-reference only |

**Excluded:**

- **GSM9809295** — GEO title `HumanBrain-Tumor-WU1220-01`, all nine deposited
  files named `WU1220_2_*` (§1.5). Pair identity unknown; the footprint
  measurement does not resolve it (1.35× against WU1220-1, 1.23× against
  WU1220-2, both plausible). **Excluded entirely rather than assigned.**
- **WU1220 and WU1221 from the headline arm** — GSE318564's scRNA donors are
  B178, B183, B189, WU-1225, WU-1226, GBM034. Neither WU1220 nor WU1221 is
  among them, so the own-reference design is **structurally impossible** for
  them. They appear only in the foreign-reference arm, where they supply the
  area-ratio gradient 1.02 → 1.34.
- Every other GSE318562 Visium section and all GSE316917 Xenium samples — no
  matched partner.

**The structural finding this table forces, registered before any result:** the
two confounds are **anti-correlated across the available data**. The only donor
with its own scRNA (B189, the one the scouting record ranks highest) is also the
**worst** region-mismatch case, at 3.57× and 2.03×. The region-cleanest pairs
(1.02×, 1.12×) have **no** own-donor scRNA. **No pair in this dataset has both
confounds controlled at once**, and no analysis can create one.

## 6. What this test cannot establish, whatever it returns

- Nothing about **section separation** (§2.1): unverified, unverifiable from
  public records, uncontrolled here.
- Nothing about **drug transport, kill fraction or clinical response**. This is a
  composition-estimator check and nothing else.
- No **cross-compound ranking** (withdrawn 2026-08-21) and no absolute kill
  fraction.
- Nothing about **Xenium ground truth**: the Xenium side is the repo's
  marker-panel classifier, an estimator, not an annotation. The authors' own
  broad-cell-type labels exist (Zenodo `paired_xenium_object_list_with_broadcelltype_annotations.rds`,
  **4.49 GB**) but are unreadable here — the file exceeds the 2.3 GiB of free
  disk and no R interpreter is installed, and a Seurat S4 object is not readable
  by the available Python stack. **Recorded as unafforded, not as unavailable.**
- Nothing about any donor other than B189 in the headline arm, and nothing about
  any tumour type other than GBM.
- A pass does **not** reinstate the retracted 0.997. That number stays retired
  whatever this returns.
