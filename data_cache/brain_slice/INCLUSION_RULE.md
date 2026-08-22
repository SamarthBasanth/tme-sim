# Inclusion rule — rat/mouse brain-slice anchor dataset

**Written 2026-08-21, BEFORE extracting values from any source other than Loryan 2013 Table 6.**
Frozen. Any change to this file after extraction must be dated and justified in place, not silently edited.

Purpose: this dataset is intended as an out-of-sample constraint on the cell/matrix partition
coefficient `K` and the intracellular accumulation term in Codon's two-route transport law.
Because we are assembling it from N papers, **we are the curators**, and curation is the main
threat to its validity. These rules exist to make the assembly mechanical.

---

## R1 — Which sources qualify

A source qualifies if ALL of the following hold:

1. It is **open access** and retrievable without institutional credentials (PMC / Europe PMC
   full text, or publisher OA). Paywalled sources are recorded in `PROVENANCE.md` as
   *known-missing*, never partially transcribed from an abstract.
2. It reports **per-compound** `Vu,brain` and/or `Kp,uu,cell` in a table or figure-with-values.
   Range statements in an abstract ("Vu,brain ranged from 2.8 to 624") are **access metadata,
   not data**, and are never recorded as values.
3. The measurement is the **rat or mouse fresh brain-slice method** of the Uppsala
   (Hammarlund-Udenaes / Fridén / Loryan) lineage — 300 um slices, oxygenated aECF, ~5 h
   incubation, Vi correction. Brain **homogenate** (equilibrium dialysis) `fu,brain` values do
   NOT qualify and are not to be inverted into `Vu,brain`.

## R2 — Species and strain

Species and strain are **recorded per row and never pooled**. Established already:
Sprague-Dawley rat and NMRI mouse are separate columns in Loryan 2013 and are separate rows here.
Pooling decisions are deferred to scoring, not made at collection time.

## R3 — No per-compound selection

**If a paper qualifies under R1, every compound in its qualifying table is taken.** No compound
is dropped for being uninteresting, mechanistically awkward, out of Codon's modelling scope,
an outlier, or a suspected active-transport substrate.

Explicitly: cyclic peptides and H+/OC antiporter substrates are **included**. Their mechanism is
tagged in a `mechanism_flag` column. Whether they are scoreable is a pre-registration decision
made later against a written rule — never a decision made after seeing how they score.

An exclusion made now on stated physical grounds is legitimate. An exclusion made after scoring
is not. This rule is what keeps those two separable.

## R4 — Duplicates (same compound in more than one source)

Duplicates are **kept as separate rows**, never averaged, never de-duplicated, never
"best-value" selected.

Rationale: the same compound measured by the same method in two independent papers is a direct
empirical estimate of assay reproducibility. Every such overlap is reported in `PROVENANCE.md`
with its ratio (larger/smaller). That ratio sets a **floor on the accuracy any model can be
expected to achieve** — a model cannot be scored as wrong at a tolerance finer than the assay
disagrees with itself.

Known overlap at time of writing: **oxycodone** appears in Loryan 2013 Table 6 (SD rat, itself
flagged as imported from Friden et al. 2011) and in the 2025 Fluids Barriers CNS oxycodone paper.

## R5 — Native vs imported values

Each row records `value_origin`:
- `native` — measured in the experiments reported by that paper;
- `imported:<citation>` — reproduced by that paper from an earlier one.

Imported values are **not** independent replicates and must not be used as an overlap check
against their own source. Known case: Loryan 2013 Table 6 oxycodone SD-rat value is
`imported:Friden 2011`.

## R6 — Units verbatim, no conversion

Units are recorded **exactly as printed**. No conversion, no rescaling, no harmonisation, no
log transformation, no unit normalisation at collection time.

- `Vu,brain` is reported in the literature as **mL·g brain^-1** and stays that way, even though
  the model may want something dimensionless.
- `Kp,uu,cell` is **dimensionless** and stays that way.

If a source prints different units, the row is **kept and flagged**, not converted and not dropped.

## R7 — Vu,brain and Kp,uu,cell are separate datasets

They constrain different things and will be scored as two different tests. They are stored in
two separate files and are never merged into one table, even where the same paper reports both
for the same compound.

## R8 — Kp,uu,brain is captured but NEVER scored

`Kp,uu,brain` may be captured when it appears in the same table, and is stored in a file whose
name marks it as not-to-be-scored.

Reason, recorded here so nobody re-derives it wrongly later: **a brain slice has no intact
blood-brain barrier.** `Vu,brain` and `Kp,uu,cell` are therefore clean in-vitro measurements of
tissue and cellular partitioning. `Kp,uu,brain` is a different thing entirely — it is an in-vivo
brain-to-plasma ratio confounded by P-gp / BCRP efflux at the BBB, which Codon does not model.
Scoring against it would be scoring Codon against a transporter it does not contain.

## R9 — Poolability with Loryan 2013 must be checked, not assumed

Before any `Vu,brain` row from another paper is treated as poolable with the Loryan 2013 Table 6
values, method, species and units must be confirmed to match. Where a paper deviates, the row is
**kept and marked**, not dropped.
