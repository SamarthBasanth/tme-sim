# The matched-pair, own-reference deconvolution test — 2026-10-10

**Status: PRE-REGISTERED TEST, RUN AND SCORED. It FAILED, in the pre-registered
adverse direction, and the pre-registered prediction was wrong.**
No file in `results/` or `archive/` was edited; this file is the only one added.
Nothing in `tmesim/**`, `app/**`, `docs/**` or `landing-3d/**` was modified.
Nothing below is registered in `tmesim/validation/scorecard.py`.

**Calibration guard, start and end of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563      <- start
alpha OK 22.223477711172563      <- end
```

---

## 0. The one-paragraph answer

Donor **B189** has its own scRNA-seq, its own Visium sections and its own Xenium
sections in three separate GEO accessions, which makes it the only tissue in
either scouting record where the deconvolution reference mismatch can be
removed. It was removed. **Agreement got worse, not better** — on four of the
five pairs the patient's own cells are 2.5× to 4× *further* from the Xenium
composition than the foreign Darmanis atlas is. The pre-registered prediction
(§3.2 C6: own-reference closer) is **refuted**. The region-mismatch confound
that the 2026-10-09 verification pass identified as load-bearing was measured
rather than argued about, and it is **real but an order of magnitude too small
to matter**: at the worst pair it can move the composition by at most 2.6
percentage points, against observed disagreements of 10–41 points. Three
mechanisms were isolated, none of them the reference-mismatch hypothesis that
`results/revalidation.md` §4c put forward.

---

## 1. Pre-registration

The protocol was written in full, saved and hashed **before any agreement number
existed**. The frozen copy is
`experiments/matched_pair_deconv/PREREGISTRATION.md`:

```
sha256  a2433a78054d3ced578a00614f244cd2ef4be938b2237bd39b90a6b533748904
written 2026-10-10T09:42:48Z (experiments/matched_pair_deconv/PREREG_TIMESTAMP_UTC.txt)
```

The same hash is recorded inside the run's own output
(`stage_b.json: prereg_sha256`).

**What that hash is and is not worth.** It proves the file scored below is the
file that was hashed. It does **not** independently prove the ordering, because
the hash and the timestamp were both produced by the same agent that wrote the
file. The strong evidence would have been a commit made before the run; no
commit was made, because committing was not requested. A reader who wants the
ordering guaranteed rather than asserted should require that of the next
pre-registration.

**Thresholds were taken from numbers already published in this repo**, so they
could not have been contaminated by the data fetched for this test. Both are
arithmetic from `results/revalidation.md` §4a/§4b:

| | source vectors | mean abs. diff | L1 threshold |
|---|---|---:|---:|
| `τ_retracted` | Xenium 73.80/7.92/8.54/9.74 vs Visium 82.66/1.04/6.53/9.77 | 4.44 pp | **0.1778** |
| `τ_current` | Xenium 58.24/11.78/19.48/10.50 vs the same Visium vector | 12.21 pp | **0.4884** |

---

## 2. The region-mismatch confound: registered first, then measured

### 2.1 Two options were rejected before any data was scored

The 2026-10-09 verification pass established (§1.4) that the Xenium half is
always the larger field, 1.02× to 3.57× by area, so the two halves of a pair do
not cover the same tissue and a pooled whole-section comparison mixes spatial
heterogeneity into the unmixing error.

**(a) Restrict to the geometric overlap — REJECTED, unaffordable, and
unverifiable even if afforded.** Co-registering the two frames needs a shared
physical referent. The Visium half ships an H&E image; the Xenium morphology
image is `*_morphology.ome.tif.gz`, and the GEO `filelist.txt` gives its size
per section:

| section | morphology image |
|---|---:|
| ST073021_2 | **10.89 GB** |
| ST073021_4 | **11.76 GB** |
| WU1220_3 | 9.22 GB |
| WU1221_Core2 | 5.81 GB |
| WU1220_2 | 7.77 GB |
| WU1221_Edge | 3.98 GB |

The two B189 images alone are **22.65 GB** against **2.3 GiB** of free disk on
this host — short by roughly tenfold, before the Visium side or any working
copy. **This is reported as unafforded, not as unavailable.** Separately, and
independently of budget: the two sections of a pair are different physical cut
planes (established on platform grounds in the verification pass §1.2), so even
a perfect H&E↔DAPI registration would be aligning two different pieces of
tissue, and the residual could not be bounded without a third reference.

**(b) Stratify by region — REJECTED, reintroduces free parameters.** Region
labels in a common frame need the same registration as (a). An intrinsic
conditioning variable is estimated differently on the two platforms, so the
stratifier would carry the confound it is meant to remove.

### 2.2 What was done instead, and what it cost

**(c) Report the result as a function of how much of the field is shared.** The
Xenium half resolves single cells over its whole field, so the confound can be
*measured in the units of the outcome*. For each pair, 2,000 random sub-windows
of the PCA-aligned Xenium field were drawn with the aspect ratio of the
PCA-aligned Visium extent, each grown about its centre until the **occupied**
tissue inside it matched the Visium footprint area to within ±10%, requiring
≥200 cells. Seed `default_rng(20261010)`. From the ensemble:

- `L1_full` — Visium composition vs **whole-field** Xenium composition (the naive pooled comparison);
- `L1_best` — **minimum** over sub-windows (the most charitable possible placement of the Visium footprint);
- `L1_wf_50/95` — median and 95th percentile of a random sub-window against the whole field, i.e. **the confound itself**.

`L1_best` is selection-biased downward, so it was pre-registered for declaring
**failure only, never a pass**. That asymmetry is what makes the construction
honest: if even the best-placed window cannot be matched, the error is real.

**The cost.** Two things were given up. First, no claim can be made about where
in the Xenium field the Visium section actually came from — only about the range
of compositions available anywhere in it. Second, `L1_best` cannot support a
positive result, so a *pass* on this design would have been inconclusive. The
design can only return "fails", "uninformative", or "weak pass"; it was accepted
on those terms before it was run.

### 2.3 The measured confound — real, and far too small to matter

2,000 windows accepted out of 2,000 attempts on every pair.

| pair | area × | `L1_wf` p50 | `L1_wf` p95 | `L1_wf` max | p95 as mean abs. diff |
|---|---:|---:|---:|---:|---:|
| P1 ST073021_2 | 3.56 | 0.0712 | **0.1033** | 0.1101 | **2.58 pp** |
| P2 ST073021_4 | 2.04 | 0.0393 | **0.0592** | 0.0668 | 1.48 pp |
| P3 WU1220_3 | 1.12 | 0.0067 | **0.0154** | 0.0162 | 0.39 pp |
| P4 WU1221_Core2 | 1.34 | 0.0122 | **0.0320** | 0.0363 | 0.80 pp |
| P5 WU1221_Edge | 1.02 | 0.0089 | **0.0183** | 0.0185 | 0.46 pp |

The band widens with the area ratio exactly as it should — 2.58 pp at 3.56×,
0.39 pp at 1.12× — which is the "as a function of how much of the field is
shared" reading the design was built to produce.

**Robustness.** A window clipped by the field boundary is grown further and
therefore over-covers, which would narrow the band spuriously. Excluding every
clipped window does **not** widen it — P1 0.1033 → 0.0933, P3 0.0154 → 0.0124,
P5 0.0183 → 0.0123; only P4 moves up, 0.0320 → 0.0339. The band survives its own
most obvious failure mode.

**Caveat, in the conservative direction.** The band is measured on the repo's
Xenium classifier's labels and the windows are constrained to ≥90% occupied
tissue, so it is a **lower bound** on true region-mismatch variability. It would
have to be wrong by a factor of 9 to 101 to explain the disagreements in §4.

### 2.4 The footprint table of 2026-10-09 is independently reproduced

Computed here from the files, not copied:

| pair | spots | Visium mm² | Xenium cells | Xenium mm² | area × here | area × in §1.4 |
|---|---:|---:|---:|---:|---:|---:|
| P1 ST073021_2 | 1,787 | 15.48 | 52,234 | 55.14 | **3.56** | 3.57 |
| P2 ST073021_4 | 2,850 | 24.68 | 185,662 | 50.44 | **2.04** | 2.03 |
| P3 WU1220_3 | 3,995 | 34.60 | 296,613 | 38.76 | **1.12** | 1.12 |
| P4 WU1221_Core2 | 1,919 | 16.62 | 103,688 | 22.24 | **1.34** | 1.34 |
| P5 WU1221_Edge | 2,045 | 17.71 | 20,859 | 18.03 | **1.02** | 1.02 |

Spot counts and cell counts match exactly; areas differ only in the grid origin.
**`results/data_scouting_threads_closed_2026-10-09.md` §1.4 is verified.**

---

## 3. Pairs used and excluded

| id | Visium | Xenium | donor | area × | arm |
|---|---|---|---|---:|---|
| **P1** | GSM9497624 `TWBK-ST073021-2` | GSM9809293 `ST073021_2` | B189 | 3.56 | own-reference + Darmanis |
| **P2** | GSM9497625 `TWBK-ST073021-4` | GSM9809294 `ST073021_4` | B189 | 2.04 | own-reference + Darmanis |
| P3 | GSM9497628 `TWKI-WU1220-3` | GSM9809296 `WU1220_3` | WU1220 | 1.12 | Darmanis + B189-as-foreign |
| P4 | GSM9497629 `TWKI-WU1221-Core2` | GSM9809298 `WU1221_Core2` | WU1221 | 1.34 | Darmanis + B189-as-foreign |
| P5 | GSM9497630 `TWKI-WU1221-Edge1` | GSM9809297 `WU1221_Edge` | WU1221 | 1.02 | Darmanis + B189-as-foreign |

**Excluded:**

- **GSM9809295** — GEO title `HumanBrain-Tumor-WU1220-01`, all nine deposited
  files named `WU1220_2_*`. The pair identity is unknown and the footprint
  measurement does not resolve it. **Excluded entirely rather than assigned to
  the more convenient partner.**
- **WU1220 and WU1221 from the own-reference arm** — GSE318564's scRNA donors
  are B178, B183, B189, WU-1225, WU-1226, GBM034. Neither is among them, so the
  own-reference design is *structurally impossible* for those pairs.
- All other GSE318562 Visium sections and all GSE316917 Xenium samples — no
  matched partner.

### The structural fact this forces

**The two confounds are anti-correlated across the available data.** The only
donor with its own scRNA — B189, the one the scouting record ranks highest — is
also the worst region-mismatch case, at 3.56× and 2.04×. The region-cleanest
pairs, 1.02× and 1.12×, have no own-donor scRNA at all. **No pair in this dataset
has both confounds controlled at once, and no analysis can create one.** This was
registered before the run, not discovered in it.

---

## 4. Results against each pre-registered criterion

### 4.1 The compositions

| pair | side | TUMOR | STROMAL | IMMUNE | ENDOTHELIAL |
|---|---|---:|---:|---:|---:|
| P1 | **Xenium whole field (truth)** | 67.97 | 12.11 | 11.46 | 8.46 |
| P1 | Visium NNLS, **B189 own reference** | 10.68 | **85.03** | 0.40 | 3.90 |
| P1 | Visium NNLS, Darmanis reference | 2.40 | **96.18** | 1.28 | 0.14 |
| P2 | **Xenium whole field (truth)** | 60.82 | 15.00 | 18.47 | 5.71 |
| P2 | Visium NNLS, **B189 own reference** | 0.54 | **90.07** | 2.46 | 6.93 |
| P2 | Visium NNLS, Darmanis reference | 87.52 | 0.99 | 9.75 | 1.75 |
| P3 | **Xenium whole field (truth)** | 69.19 | 18.28 | 10.02 | 2.50 |
| P3 | Visium NNLS, B189 reference | 1.51 | **95.94** | 0.53 | 2.01 |
| P3 | Visium NNLS, Darmanis reference | 92.18 | 4.37 | 1.29 | 2.16 |
| P4 | **Xenium whole field (truth)** | 65.83 | 15.92 | 13.39 | 4.86 |
| P4 | Visium NNLS, B189 reference | 0.55 | **98.25** | 0.25 | 0.95 |
| P4 | Visium NNLS, Darmanis reference | **98.76** | 0.09 | 0.50 | 0.65 |
| P5 | **Xenium whole field (truth)** | 60.83 | 18.19 | 12.10 | 8.87 |
| P5 | Visium NNLS, B189 reference | 0.66 | **98.97** | 0.00 | 0.36 |
| P5 | Visium NNLS, Darmanis reference | 65.71 | 34.08 | 0.18 | 0.03 |

### 4.2 C1, C2, C3 — the 4-vector

| pair | area × | band p95 | reference | `L1_full` | `L1_best` | mean abs. diff at best | **C1** | **C2** | **C3** |
|---|---:|---:|---|---:|---:|---:|---|---|---|
| P1 | 3.56 | 0.1033 | **own B189** | 1.4583 | **1.3910** | **34.77 pp** | **FAIL** | no | outside band |
| P1 | 3.56 | 0.1033 | Darmanis | 1.6815 | 1.6141 | 40.35 pp | **FAIL** | no | outside band |
| P2 | 2.04 | 0.0592 | **own B189** | 1.5258 | **1.4645** | **36.61 pp** | **FAIL** | no | outside band |
| P2 | 2.04 | 0.0592 | Darmanis | 0.5340 | 0.4909 | 12.27 pp | **FAIL** | no | outside band |
| P3 | 1.12 | 0.0154 | B189 (foreign) | 1.5533 | 1.5460 | 38.65 pp | **FAIL** | no | outside band |
| P3 | 1.12 | 0.0154 | Darmanis | 0.4597 | 0.4435 | 11.09 pp | pass | no | outside band |
| P4 | 1.34 | 0.0320 | B189 (foreign) | 1.6466 | 1.6399 | 41.00 pp | **FAIL** | no | outside band |
| P4 | 1.34 | 0.0320 | Darmanis | 0.6587 | 0.6240 | 15.60 pp | **FAIL** | no | outside band |
| P5 | 1.02 | 0.0183 | B189 (foreign) | 1.6156 | 1.6085 | 40.21 pp | **FAIL** | no | outside band |
| P5 | 1.02 | 0.0183 | Darmanis | 0.4153 | 0.4030 | 10.08 pp | pass | no | outside band |

- **C1 (`L1_best` > `τ_current` = 0.4884 ⇒ failure): FAILS on both headline
  pairs with the own reference**, and on 8 of the 10 pair×reference runs. On the
  own-reference arm the most charitably placed comparison is **2.85× and 3.00×
  worse** than the project's current unmatched, foreign-reference figure.
- **C2 (`L1_best` ≤ `τ_retracted` = 0.1778): not met anywhere.** The closest any
  run comes is Darmanis on P5 at 0.4030, still 2.3× above it.
- **C3: `L1_full` is outside the measured region band on every single pair and
  every reference**, by factors of 9.0 (Darmanis, P2) to 100.7 (own reference, P3).
  There is **no "uninformative" verdict available**: region mismatch cannot
  account for any of this.

### 4.3 C4 — the non-tumour 3-vector, the subvector that collapsed

| pair | band p95 | reference | `L1nt_full` | `L1nt_best` |
|---|---:|---|---:|---:|
| P1 | 0.2288 | own B189 | 1.1475 | **0.9824** |
| P1 | 0.2288 | Darmanis | 1.2147 | 1.0495 |
| P2 | 0.1314 | own B189 | 1.0454 | **0.9230** |
| P2 | 0.1314 | Darmanis | 0.6192 | 0.4765 |
| P3 | 0.0234 | B189 | 0.7615 | 0.7350 |
| P3 | 0.0234 | Darmanis | 0.3899 | 0.3819 |
| P4 | 0.0222 | B189 | 1.0441 | 1.0280 |
| P4 | 0.0222 | Darmanis | 0.7882 | 0.7754 |
| P5 | 0.0235 | B189 | 1.0637 | 1.0431 |
| P5 | 0.0235 | Darmanis | 1.0587 | 1.0381 |

**C4 fails everywhere, by 3.6× to 46× over the band.** A maximum possible L1 on a
renormalised 3-vector is 2.0; the own-reference arm reaches 0.92–1.04, i.e. half
of the worst value attainable. **The 0.016 Pearson collapse recorded in
`results/revalidation.md` §4 is reproduced in a matched-donor, own-reference,
region-controlled design, and removing the reference mismatch did not touch it.**

### 4.4 C5 — Spearman ρ, reported and NOT a gate

P1 own +0.60 / Darmanis +0.80 · P2 own **−0.80** / Darmanis +0.80 · P3 own 0.00 /
Darmanis +0.80 · P4 own 0.00 / Darmanis +0.20 · P5 own +0.60 / Darmanis +1.00.

Declared non-gating in advance, and it is worth seeing why that mattered:
Darmanis on P5 scores ρ = **+1.00**, a perfect ordering match, on a pair whose
composition is wrong by 10.08 pp per class. With K = 4 the statistic is too
coarse to carry a verdict, and had it not been fenced off beforehand the
ρ = +1.00 would have been available to quote.

### 4.5 C6 — the headline contrast. **The pre-registered prediction is REFUTED.**

Predicted: `L1_best(own B189) < L1_best(Darmanis)` on the same pairs and windows.

| pair | `L1_best` own | `L1_best` Darmanis | own better? |
|---|---:|---:|---|
| P1 | 1.3910 | 1.6141 | yes |
| P2 | **1.4645** | **0.4909** | **NO — 3.0× worse** |
| P3 | **1.5460** | **0.4435** | **NO — 3.5× worse** |
| P4 | **1.6399** | **0.6240** | **NO — 2.6× worse** |
| P5 | **1.6085** | **0.4030** | **NO — 4.0× worse** |

**Refuted on 4 of 5 pairs, including one of the two headline B189 pairs.** The
single pair where own-reference wins (P1) is the one with the worst area ratio
and the lowest Visium depth, and both of its numbers — 1.39 and 1.61 — are
catastrophic in absolute terms, so P1 is a contest between two failures.

**What this does to the standing hypothesis.** `results/revalidation.md` §4c
attributed the composition disagreement to Visium NNLS against the Darmanis
atlas having "its own, different tumour-inflating tendency". That is a
reference-mismatch hypothesis, and it is the one this test was built to
examine. **Giving the deconvolution the patient's own cells makes the
disagreement substantially worse.** Reference mismatch is therefore not the
explanation, and the project's position is worse than it was yesterday: the open
question stays open, and one of the two candidate explanations has been
eliminated in the unfavourable direction.

### 4.6 C7 — the unmixing-only arm, where region mismatch is identically zero

Xenium cells aggregated into Visium-geometry pseudo-spots (hex lattice, 100 µm
pitch, cells within 27.5 µm of a lattice point, ≥5 cells), deconvolved with the
same B189 signature restricted to the 362 panel genes, scored against each
pseudo-spot's own known composition. Truth and mixture are the **same cells**, so
region mismatch, section separation and platform are all exactly zero.

| pair | pseudo-spots | mean per-spot L1 | pooled L1 | pooled mean abs. diff |
|---|---:|---:|---:|---:|
| P1 | 761 | 0.7924 | 0.4601 | **11.50 pp** |
| P2 | 4,559 | 0.9226 | 0.6798 | **16.99 pp** |

| pair | side | TUMOR | STROMAL | IMMUNE | ENDOTHELIAL |
|---|---|---:|---:|---:|---:|
| P1 | truth | 67.85 | 11.13 | 11.42 | 9.61 |
| P1 | predicted | 63.87 | **34.13** | **1.25** | **0.74** |
| P2 | truth | 60.38 | 15.67 | 19.06 | 4.89 |
| P2 | predicted | 50.03 | **49.66** | **0.29** | **0.02** |

**The unmixing fails on its own, with every other confound removed.** It drives
IMMUNE and ENDOTHELIAL to near-zero — 19.06% → 0.29% and 4.89% → 0.02% on P2 —
and dumps the mass into STROMAL. Per the pre-registered reading of C7, a large
C7 error is partly attributable to the 362-gene restriction; but a 19-point
immune fraction collapsing to 0.3 points is not a gene-count effect, and the
direction is identical to the full Visium arm, which had 17,941 genes.

### 4.7 Within-donor reproducibility — no Xenium truth needed for the Visium rows

P1 and P2 are two sections of the same tumour. Comparing each estimator with
itself:

| quantity | L1 | mean abs. diff |
|---|---:|---:|
| **Xenium truth**, P1 vs P2 | 0.1980 | 4.95 pp |
| Visium NNLS, **B189 own reference**, P1 vs P2 | 0.2028 | 5.07 pp |
| Visium NNLS, **Darmanis reference**, P1 vs P2 | **1.9040** | **47.60 pp** |

Three distinct readings, and they separate cleanly:

1. The **Xenium side is internally reproducible** across two sections of one
   tumour (4.95 pp). It behaves like a measurement.
2. The **own-reference Visium side is equally reproducible (5.07 pp) and
   systematically wrong** — STROMAL 85–99% against a truth of 12–15%. A stable
   estimator of the wrong quantity.
3. The **Darmanis Visium side is not even internally reproducible**: the same
   reference and the same pipeline call one B189 section 96.18% STROMAL and the
   other 87.52% TUMOR. **This needs no cross-platform comparison and no Xenium
   half to be a problem.**

---

## 5. Three mechanisms, isolated

### 5.1 The marker-panel classifier does not transfer to scRNA: the "ENDOTHELIAL" class is PECAM1-positive myeloid cells

The reference was labelled with the repo's shipped
`real_xenium.phenotype_cells` (threshold 1.0, `matched_null_z`, pool 25), the
same function, markers and gene space used on the Xenium side. Detection rate
(% of cells with count > 0) by assigned class, B189 scRNA:

| assigned class | PECAM1 | PTPRC | CD68 | AIF1 | CD163 | AQP4 | GJA1 | RGS5 | PDGFRB | n |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| TUMOR | 0.8% | 19.2% | 24.6% | 25.3% | 8.8% | 34.5% | 12.7% | 1.4% | 0.7% | 7,494 |
| STROMAL | 1.4% | 5.3% | 10.6% | 10.3% | 3.5% | 92.1% | 71.4% | 2.7% | 2.6% | 1,487 |
| IMMUNE | 15.9% | 88.9% | 90.6% | 92.4% | 82.4% | 4.4% | 1.0% | 0.6% | 0.1% | 2,296 |
| **ENDOTHELIAL** | **100.0%** | **93.5%** | **96.0%** | **97.3%** | **67.1%** | 4.5% | 1.4% | **0.6%** | **0.0%** | 1,718 |

The ENDOTHELIAL class is 100% PECAM1-positive **and** 93–97% positive for all
three myeloid markers. It is the PECAM1-positive subset of the immune
compartment. It is not pericytes either — RGS5 0.6%, PDGFRB 0.0%.

**Why the classifier behaves differently here.** On a 362-gene Xenium panel,
PECAM1 detection is sparse and discriminative. In deep whole-transcriptome
polyA scRNA, PECAM1 is detected in essentially every macrophage, so a
single-marker PECAM1 gate selects immune cells rather than endothelium. The
classifier was built and audited for Xenium (`results/classifier_fix.md`) and
**nothing in this repo claimed it transfers to scRNA** — this test is the first
thing to have tried it, and it does not.

### 5.2 That makes the own-donor signature ill-conditioned

Cosine similarity between signature columns, over the 17,941 genes shared with
the P1 Visium section:

| | own B189 | Darmanis |
|---|---:|---:|
| cos(IMMUNE, ENDOTHELIAL) | **0.9777** | 0.4802 |
| cos(TUMOR, STROMAL) | 0.6699 | 0.7551 |
| condition number | **10.20** | 4.46 |

Two columns at cosine 0.978 are very nearly the same vector, so NNLS cannot
apportion weight between them and both lose to whichever column is better
conditioned. That is exactly the observed failure: IMMUNE and ENDOTHELIAL both
collapse while STROMAL absorbs the mass.

**This also refutes the direction declared in advance.** §3.3 of the
pre-registration predicted that swapping Darmanis → B189's own cells would move
the Visium ENDOTHELIAL fraction **DOWN**, on the hypothesis that Darmanis's
single `Vascular` column absorbs pericyte signal. Measured: 0.14 → 3.90 (P1),
1.75 → 6.93 (P2), 2.16 → 2.01 (P3), 0.65 → 0.95 (P4), 0.03 → 0.36 (P5).
**It went UP in four of five pairs — the declared direction was wrong**, and
wrong for a reason opposite to the one hypothesised: the own-reference
endothelial column is not narrower than Darmanis's, it is an immune-contaminated
column collinear with IMMUNE.

### 5.3 The Visium NNLS returns an arg-max, not a mixture — and much more so on FFPE

Per-spot solutions from the repo's shipped pipeline, Darmanis reference
(uniform-vector entropy = 1.386):

| section | chemistry | counts/spot | spots with max fraction > 0.99 | mean entropy |
|---|---|---:|---:|---:|
| `10x_parent_gbm` (published) | fresh-frozen | — | 22.7% | 0.408 |
| `GBM1` | fresh-frozen | — | 1.2% | 0.566 |
| `GBM2` | fresh-frozen | — | 0.1% | 0.960 |
| `GBM3` | fresh-frozen | — | 0.1% | 0.911 |
| P1 ST073021-2 | **FFPE** | 1,767 | **73.9%** | **0.091** |
| P2 ST073021-4 | **FFPE** | 5,734 | 27.8% | 0.283 |
| P3 WU1220-3 | **FFPE** | 3,762 | **75.7%** | **0.125** |
| P4 WU1221-Core2 | **FFPE** | 20,904 | **86.0%** | **0.040** |
| P5 WU1221-Edge1 | **FFPE** | 12,892 | 9.6% | 0.567 |

Two things follow. First, **the degeneracy is not a depth effect** — P4 has the
deepest spots (20,904 counts) and is the most degenerate (86.0%). Second, **the
estimator is markedly more degenerate on FFPE probe-based Visium than on the
fresh-frozen sections it was developed against.** All six matched pairs are FFPE
on both halves; the repo's Visium deconvolution has, until now, only ever been
run on fresh-frozen sections. When 74–86% of spots are assigned a single class at
>0.99, the section-level mean is an arg-max histogram, not a composition, and
comparing it with a Xenium composition compares objects of different kinds.

The published 82.66% TUMOR Visium vector is itself the most degenerate of the
four fresh-frozen sections (22.7% of spots above 0.99). That is recorded here as
an observation about an existing published number, not a claim about it.

---

## 6. Verification of the machinery, so the failures cannot be blamed on this pass

| check | result |
|---|---|
| Repo's own `deconvolve_sample("10x_parent_gbm")` under this environment | **TUMOR 82.66 / STROMAL 1.04 / IMMUNE 6.53 / ENDOTHELIAL 9.77** — the published vector, to every digit |
| Sparse signature build vs the repo's own `build_signatures` | max abs. difference **5.551e-17** |
| Visium in-tissue spot counts vs the 2026-10-09 table | 1,787 / 2,850 / 3,995 / 1,919 / 2,045 — exact |
| Xenium cell counts vs the 2026-10-09 table | 52,234 / 185,662 / 296,613 / 103,688 / 20,859 — exact |
| All 21 downloaded files vs the GEO `filelist.txt` byte counts | **21 OK, 0 problems** |
| Classifier determinism | `matched_null_z` computes its null moments analytically; RNG-free |

Every scoring component is the repo's shipped default, imported and called
unmodified. **Nothing was tuned, and no estimator variant was explored to find
one that agreed better.** The failures above are the shipped pipeline's.

---

## 7. Deviations from the pre-registration, declared

1. **The panel is 362 genes, not 542.** The pre-registration says "the 542-gene
   GSE335275 panel". 542 is the *feature* count; `feature_types` resolves to 362
   `gene`, 111 unassigned codewords, 41 negative-control codewords, 20
   negative-control probes, 7 deprecated codewords and 1 `aggregate_gene`
   pseudo-feature. A factual correction, not a protocol change.
2. **The first stage-B run labelled Xenium cells using all 542 columns**,
   including the 180 control features, which broke the pre-registration's own
   requirement of the same gene space on both sides (the scRNA side had used
   exactly the 362 genes). It also departed from the repo's loader default
   (`singlecell_spatial.load(genes_only=True)` drops control and codeword
   features). It was re-run restricted to the 362 gene features. The rule for
   choosing between them was fixed by the code and the pre-registration, not by
   the outcome: **match the repo's loader convention.** Both runs are kept —
   `stage_b.json` (corrected, primary) and
   `stage_b_WITH_CONTROL_PROBES.json` — and the difference is immaterial: the
   Xenium truth compositions move by at most **0.25 pp**, no verdict changes,
   and the largest shift in any `L1_best` is 0.0034.
3. **The WU1221_Edge Xenium matrix was downloaded before the pre-registration
   was written**, to establish the archive format and the panel gene list. No
   composition, label or agreement number was computed from it before the
   pre-registration was saved, and it belongs to P5, a foreign-reference pair,
   not to either headline pair.
4. **Sub-window sizing.** The pre-registration says windows have "area `A_vis`"
   and are accepted on occupied area within ±10%. Implemented as: a rectangle of
   the Visium aspect ratio, grown about its centre until the occupied tissue
   inside it meets the ±10% criterion. This is the mechanical realisation of the
   acceptance rule, not a change to it.

---

## 8. What this does NOT establish

- **Nothing about section separation.** It remains UNVERIFIED and unverifiable
  from public records, and uncontrolled here. It can only inflate disagreement,
  never deflate it, so it cannot have manufactured the failures above — but it
  is a live term in every number, and no result here is a statement about how far
  apart the two sections are.
- **Nothing about the authors' own cell-type annotations.** The Xenium truth used
  here is the repo's marker-panel classifier, which is an estimator, not an
  annotation. The authors' own broad cell-type labels exist
  (Zenodo 10.5281/zenodo.20709378,
  `paired_xenium_object_list_with_broadcelltype_annotations.rds`, **4.49 GB**) and
  were **not read**: the file alone is nearly twice the free disk on this host,
  no R interpreter is installed, and a Seurat S4 object is not readable by the
  available Python stack. **Recorded as unafforded, not as unavailable** — it is
  the single highest-value next download, and it would convert the truth side
  from an estimator into an annotation.
- **Nothing about the correct composition of these tumours.** This test says the
  Visium deconvolution disagrees with the Xenium classifier far beyond what
  region mismatch can explain. It does not say which of the two is right, and
  §5.1 gives a concrete reason to doubt the classifier on at least one modality.
- **Nothing about drug transport, kill fraction, or clinical response.** This is
  a composition-estimator check and nothing else. No cross-compound ranking
  (withdrawn 2026-08-21), no absolute kill fraction, no clinical, diagnostic or
  AI/ML framing.
- **Nothing about any donor but B189 in the own-reference arm**, nothing about
  any tumour type but GBM, and nothing about fresh-frozen Visium — all five pairs
  are FFPE on both halves.
- **It does not reinstate the retracted Pearson 0.997.** That number stays
  retired. Nor does it produce a replacement number: the correct reading is not
  "the agreement is X", it is **"the Visium composition estimate is not fit to be
  compared, and the matched pairs made that visible rather than fixing it."**
- **A pass was not available from this design.** `L1_best` was pre-registered for
  declaring failure only, so the strongest positive outcome reachable was "weak
  pass, band width quoted". The failure margins — 4× to 110× over the band — are
  far outside any reading that would have depended on that asymmetry.

---

## 9. What would move this next, in order

1. **The authors' own Xenium annotations** (4.49 GB Zenodo RDS). Needs ~6 GB of
   free disk and an R interpreter, neither available on this host. It replaces
   the estimator on the truth side with an annotation and would settle whether
   §5.1's classifier problem is one-sided.
2. **A reference labelled by something other than the Xenium marker panel.** The
   own-donor reference failed because the labels were wrong, not because the
   donor was wrong. The design is still the right design; the labelling step is
   what broke, and §5.1 says precisely where.
3. **One email to `apetti@mgh.harvard.edu`.** Section separation is the one
   confound no analysis here can touch, and the corresponding author is the only
   remaining route to it.

---

## 10. Files

| path | what |
|---|---|
| `experiments/matched_pair_deconv/PREREGISTRATION.md` | the frozen protocol, sha256 `a2433a78…8904` |
| `experiments/matched_pair_deconv/PREREG_SHA256.txt`, `PREREG_TIMESTAMP_UTC.txt` | hash and write time |
| `experiments/matched_pair_deconv/run_test.py` | the two-stage run |
| `experiments/matched_pair_deconv/io_loaders.py` | GEO readers (zarr, MEX, parquet) |
| `experiments/matched_pair_deconv/control_checks.py`, `control_degeneracy.py` | §5.3 and §6 |
| `experiments/matched_pair_deconv/diagnostics.py`, `check_controls.py`, `check_feature_types.py` | §5.1, §2.3 robustness, §7.1 |
| `experiments/matched_pair_deconv/verify_bytes.sh`, `fetch.sh`, `extract_visium.sh` | download and byte verification |
| `experiments/matched_pair_deconv/summarise.py` | the tables in §4 |

Downloads live under `/tmp/vmp/` (not in the repo): 752 MB, 21 files, all byte-verified.
`data_cache/{scref_gbm,visium_gbm,xenium_gbm}` in this worktree are symlinks to
the main checkout's caches.
