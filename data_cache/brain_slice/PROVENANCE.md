# Brain-slice anchor dataset — provenance and access record

Assembled 2026-08-21. Governed by `INCLUSION_RULE.md` in this directory, which was written and
frozen **before** any value outside Loryan 2013 Table 6 was extracted.

Purpose: an out-of-sample constraint on the cell/matrix partition coefficient `K` and the
intracellular accumulation term in Codon's two-route transport law. Neither has ever been
independently constrained.

**All values are transcribed verbatim. Nothing has been converted, rescaled, harmonised,
log-transformed or unit-normalised.** `Vu,brain` is left in mL·g brain^-1 even though the model
may want something dimensionless. That conversion is a decision for the model owner, not for
this dataset.

Bibliographic data in this file was retrieved from **PubMed / Europe PMC**.

---

## 1. Files

| File | Contents | Rows |
|---|---|---|
| `vu_brain.csv` | `Vu,brain`, mL·g brain^-1 | 77 (52 scoreable) |
| `kp_uu_cell.csv` | `Kp,uu,cell`, dimensionless | 9 (4 scoreable, brain parenchyma) |
| `kp_uu_brain_DO_NOT_SCORE.csv` | `Kp,uu,brain` / `Kp,uu,WB` / `Kp,uu,SC` / `Kp,uu,CSF` | 20 |
| `INCLUSION_RULE.md` | Frozen collection rule | — |

`Vu,brain` and `Kp,uu,cell` are kept as **two separate datasets** and are two separate tests.
They are not merged, even where one paper reports both for the same compound.

---

## 2. THE HEADLINE NUMBERS

### `Vu,brain` — the test is viable

| Subset | n | min | max | spread |
|---|---|---|---|---|
| Rat, scoreable | 35 | 0.436 | 4184 | **9596x** |
| Rat, excluding cyclic peptides | 33 | 3.2 | 4184 | **1308x** |
| Mouse, scoreable | 17 | 1.2 | 1930 | 1608x |
| Loryan 2013 SD rat alone | 8 native | 4.49 | 2650 | 590x |

For comparison, the dataset the transport law was originally built on spans **29.7x**. Even the
most conservative subset here is **44x wider**.

### `Kp,uu,cell` — the test is NOT viable on open-access data alone

**Only 4 scoreable brain-parenchyma values were obtainable**, spanning 0.18 to 9.73 = **54x**.

Four points is not enough to discriminate a partition coefficient. The ~40-compound
`Kp,uu,cell` table that would make this test work exists, and is paywalled (section 4).

---

## 3. Sources included

All retrieved with **no institutional credentials**.

| Short name | Citation | PMID / PMCID | DOI | Access route |
|---|---|---|---|---|
| Loryan 2013 | Loryan I, Fridén M, Hammarlund-Udenaes M. *Fluids Barriers CNS* 10:6 | 23336814 / PMC3602653 | [10.1186/2045-8118-10-6](https://doi.org/10.1186/2045-8118-10-6) | Europe PMC `fullTextXML` |
| Fridén 2014 | Fridén M et al. *EJNMMI Res* | 26116114 / PMC4452637 | [10.1186/s13550-014-0050-6](https://doi.org/10.1186/s13550-014-0050-6) | Europe PMC + SpringerOpen ESM `MOESM1` .docx |
| Cyclic peptides 2023 | *Pharmaceutics* 15:1507 | 37242750 / PMC10222203 | [10.3390/pharmaceutics15051507](https://doi.org/10.3390/pharmaceutics15051507) | Europe PMC `fullTextXML` |
| CE-123 2023 | *Int J Mol Sci* 24:16956 | 38069277 / PMC10707468 | [10.3390/ijms242316956](https://doi.org/10.3390/ijms242316956) | Europe PMC `fullTextXML` |
| hMDR1 mouse 2015 | *PLoS ONE* | 25932627 / PMC4416786 | [10.1371/journal.pone.0118638](https://doi.org/10.1371/journal.pone.0118638) | Europe PMC `fullTextXML` |
| SLOH/SLM 2017 | *Sci Rep* | 29180684 / PMC5703966 | [10.1038/s41598-017-16635-4](https://doi.org/10.1038/s41598-017-16635-4) | Europe PMC `fullTextXML` |
| H+/OC antiporter 2024 | *Front Cell Neurosci* | 39534684 / PMC11554538 | [10.3389/fncel.2024.1493644](https://doi.org/10.3389/fncel.2024.1493644) | Europe PMC — `Kp,uu` only, DO-NOT-SCORE file |

---

## 4. NOT obtainable — and exactly why

These are the papers the original brief targeted. **Three of the four are hard-paywalled.**

| Paper | Status | Barrier |
|---|---|---|
| Loryan 2013 (PMC3602653) | **OBTAINED** | none |
| Fridén 2009, *DMD* 37:1226, PMID 19299522, [10.1124/dmd.108.026377](https://doi.org/10.1124/dmd.108.026377) | **NOT obtained** | Europe PMC: `isOpenAccess=N`, not in PMC. Listed "free after 12 months" at ASPET but `dmd.aspetjournals.org` returns **HTTP 403** to non-browser clients. DMD has since migrated to Elsevier; the ScienceDirect PII route also returns **403**. |
| Fridén 2009, *J Med Chem* 52:6233, PMID 19764786, [10.1021/jm901036q](https://doi.org/10.1021/jm901036q) | **NOT obtained** | ACS. `isOpenAccess=N`. SI path 302-redirects to the article landing page behind Cloudflare; no `suppl_file` anchor is rendered to an unauthenticated client. |
| Loryan 2014, *Pharm Res* 31:2203, PMID 24623476, [10.1007/s11095-014-1319-1](https://doi.org/10.1007/s11095-014-1319-1) | **NOT obtained** | Springer. `isOpenAccess=N`, `inPMC=N`, "Subscription required". |
| **Loryan 2015, *Mol Pharm* 12:520, PMID 25496026, [10.1021/mp5005965](https://doi.org/10.1021/mp5005965)** | **NOT obtained** | ACS, same Cloudflare/SI behaviour. **Not in the original brief — found during this task.** This is the single highest-value missing item: it reports `Kp,uu,brain`, `Kp,uu,cell` and `Vu,brain` together for the same 40 compounds, in Supporting Information Table S1. |

Also checked and rejected as a source: Uppsala's DiVA open repository (hosts the Loryan 2013 PDF
but holds no deposited copy of the 2014/2015 papers and no relevant thesis).

**What it would take:** an institutional subscription (ACS + Springer), per-article purchase, or
an email to the corresponding authors (Irena Loryan / Margareta Hammarlund-Udenaes, Uppsala
University; Markus Fridén). Author contact is cheap and these are 10+ year-old datasets.

### Range statements seen but deliberately NOT recorded as data

The ungated Loryan 2015 front matter states the 40-compound set spans `Kp,uu,cell` 0.15–24 and
`Vu,brain` 2.8–624 mL·g brain^-1. These are recorded here as **access metadata only**. They are
not in any CSV. Note the `Vu,brain` range there is *narrower* than what we already hold.

---

## 5. THE TRAP IN FRIDÉN 2014 — read before using that file

The EJNMMI paper's Table S1 lists 47 `Vu,brain` values. **They are not all measurements.**

Table S1 footnote **h** states verbatim:

> "Measured value in brain slices for tracers with compound material available. **Calculated from
> reported fu,brain for tracers with values of fu,brain tabulated.**"

So the rows carrying an `fu,brain` value have a `Vu,brain` that was **computed from a
pH-partition model**, not measured. Scoring Codon against those would be scoring Codon against
somebody else's model — the exact hazard flagged for ML compilations, sitting inside a primary
paper.

Applying footnote h mechanically:
* **23 rows measured** in rat brain slices → `scoreable = yes`
* **24 rows model-derived** → `scoreable = NO`, retained and flagged, never deleted

**Unresolved discrepancy, flagged not papered over:** the body text says `Vu,brain` was determined
for **31** tracers by slice and **16** calculated. Footnote h applied to Table S1 gives **23 / 24**.
These do not reconcile from the files available. I have taken the **conservative** reading — only
the 23 unambiguous rows are marked scoreable. If the body text is right, up to 8 further rows are
real measurements being needlessly withheld. Resolving this needs author contact.

---

## 6. Duplicate compounds = a free measurement of assay reproducibility

Per rule R4, duplicates are kept as separate rows and never averaged. Each overlap is an
empirical estimate of how well this assay reproduces itself, and therefore **a floor on the
accuracy any model can fairly be held to.**

### Genuine independent cross-paper replicates (rat, same method, different labs/papers)

| Compound | Loryan 2013 (SD rat) | Fridén 2014 (rat) | ratio |
|---|---|---|---|
| Paroxetine | 714 | 876 | **1.23x** |
| Diazepam | 17.8 | 20 | **1.12x** |

**These two are the most useful numbers in this file.** Two independent papers, same method,
agree to within ~12–23%. A model scored at a tolerance finer than roughly **1.25x** on rat
`Vu,brain` is being scored against noise.

### Cross-strain mouse replicates (Loryan 2013 NMRI vs hMDR1 paper C57BL/6)

| Compound | NMRI | C57BL/6 | ratio |
|---|---|---|---|
| Digoxin | 44.9 | 44.8 | 1.002 |
| Verapamil | 47.3 | 48.2 | 1.019 |
| Docetaxel | 796 | 789 | 1.009 |
| Oxycodone | 3.75 | 3.75 | 1.000 |

**Treat this agreement with suspicion, do not celebrate it.** Sub-2% agreement across different
mouse strains and different papers is implausibly good for a 5-hour tissue assay, and the
oxycodone values are *bit-identical*. The likeliest explanation is that the later paper reused
earlier values rather than measuring independently — though its Methods section explicitly
describes measuring `Vu,brain` in slices from four C57BL/6 WT mice. **Do not use this set as a
reproducibility estimate.** Use the rat cross-paper pair above.

### Cross-species, within Loryan 2013 (SD rat vs NMRI mouse)

Ratios range 1.02x (verapamil) to 1.37x (digoxin, thioridazine). Rows are kept separate and
**must not be pooled**; this is recorded only to show species is a small effect relative to the
1300x dynamic range.

### Not an independent replicate

Loryan 2013's SD-rat **oxycodone** value (4.20) is footnoted "Data from Fridén et al., 2011". It
is `value_origin = imported` and `scoreable = no`. It must never be used as a check against its
own source.

---

## 7. Other flags carried in the CSVs

* **Internal inconsistencies found in the sources themselves**, recorded rather than silently fixed:
  * Fridén 2014: `[11C]AZD2995` = **7** in main-text Table 1, **6.8** in Table S1.
  * Fridén 2014: `FE-PE2I` labelled `[18F]` in Table 1, `[11C]` in Table S1.
  * hMDR1 2015: oxycodone `Vu,brain` = **3.75** in Table 1, **4.0** in the body text.
* **Units differ cosmetically between papers** (`mL·g brain^-1`, `mL/g brain`, `ml/g_brain`).
  Recorded verbatim per rule R6. They appear to denote the same quantity — **confirm before
  pooling; do not assume.**
* **Mechanism tags, not exclusions** (rule R3): cyclic peptides, PET tracers, H+/OC antiporter
  substrates, P-gp/BCRP substrates are all **included** and tagged. Whether they are in scope for
  Codon is a pre-registration decision to be made against a written rule, *before* seeing scores.
* **SLOH/SLM**: the AD-transgenic arm is a disease model, not a healthy control. Flagged.
* **Cyclic peptide neuron/astrocyte `Kp,uu,cell`**: isolated cell culture, a different system from
  the brain slice. Marked not scoreable.
* **`Kp,uu,cell` SDs for S-CE-123 and R-modafinil** are error-propagated, not replicate SDs.

---

## 8. Why `Kp,uu,brain` is captured but must NEVER be scored

`kp_uu_brain_DO_NOT_SCORE.csv` exists so these values are not lost, and is named so nobody scores
them by accident.

**A brain slice has no intact blood-brain barrier.** `Vu,brain` and `Kp,uu,cell` are therefore
clean in-vitro measurements of tissue and cellular partitioning — which is what Codon's transport
law actually describes.

`Kp,uu,brain` is a different quantity entirely: an **in-vivo** unbound brain-to-plasma ratio,
confounded by P-gp and BCRP efflux at the BBB. **Codon does not model those transporters.**
Scoring against `Kp,uu,brain` would score Codon against machinery it does not contain, and would
produce a failure that says nothing about the partition coefficient.

The same applies to `Kp,uu,WB`, `Kp,uu,SC` and `Kp,uu,CSF` in that file.

---

## 9. Honest bottom line

* The **`Vu,brain` test is real and worth building.** 52 scoreable values, 35 in rat, spanning
  three orders of magnitude, from 6 independent open-access papers, with a measured
  reproducibility floor of ~1.25x.
* The **`Kp,uu,cell` test is not currently viable.** 4 open-access values is not a test. The
  intracellular accumulation term therefore **remains unconstrained**, and that is the honest
  status until Loryan 2015's SI is obtained.
* `Vu,brain` constrains tissue partitioning. It does **not** independently constrain the
  intracellular term. Do not let a good result on the first test be reported as validating both.
