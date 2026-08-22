# Erratum: the Abraham LSER coefficients are right, their citation is not

**Scope.** Four citation defects attached to the derivation of
`HBD_LOGP_DECREMENT = 1.34`. **No number changes. `3.621`, `0.37` and `1.34` all
stand**, and the arithmetic `0.034 − (−3.587) = 3.621` reproduces exactly against
a source I retrieved and read myself. What is wrong is *where the repo says the
numbers came from* — in a provenance chain whose entire purpose is to let a
reader check them.

**Nothing was edited to produce this file.** Three of the affected sites are
**frozen pre-registrations**, whose value depends on their never being revised
after the fact; the correct remedy for those is an erratum beside them, which is
what this is. The remaining sites are existing files under active edit by other
work, and this investigation was scoped to new files. **Every correction is given
below as exact replacement text so it can be applied in one pass**, by whoever
holds those files.

Found while verifying §2 of `results/hbd_convention_audit.md`.

## Status — 2026-08-21

**The four non-frozen sites have now been corrected**, with the replacement text
given below applied as written:

| file | site | state |
|---|---|---|
| `experiments/two_route_efflux.py` | the `# VERIFIED, verbatim, in the open-access coefficient tables of` comment block | **corrected** |
| `experiments/two_route_efflux.py` | the `THE SOURCE, AND EXACTLY HOW SOURCED IT IS` printed block | **corrected** |
| `results/two_route_efflux.md` | the `> Both read verbatim from…` block quote | **corrected**, and the withdrawn hexane/heptane cross-check is now stated inline rather than silently dropped |
| `tmesim/molecular/two_route.py` | the `#: Both coefficients read VERBATIM…` comment above `HBD_LOGP_DECREMENT` | **corrected** |

**The three frozen pre-registrations were not touched, and must not be.** This
file is their remedy. `HBD_LOGP_DECREMENT` itself, the value `1.34`, the
`Lipinski.NumHDonors` pairing and every assertion in `tests/test_two_route.py`
are unchanged — no test pins a citation string, so none moved.

A repo-wide search confirms the four wrong strings now survive **only** in the
three frozen registrations and in the two audit documents that quote them as
evidence.

---

## The four defects

### D1 — the ACS Omega volume, issue and pages are wrong

Cited everywhere as *"Abraham MH, Acree WE Jr. ACS Omega 2019;**4(7):12920-12926**,
doi:10.1021/acsomega.8b03242, PMC6648601"*.

Resolved against the Crossref API directly:

```
title  : Solvation Descriptors for Zwitterionic α-Aminoacids; Estimation of
         Water–Solvent Partition Coefficients, Solubilities, and Hydrogen-Bond
         Acidity and Hydrogen-Bond Basicity
journal: ACS Omega        published: 2019-02-28
volume : 4    issue: 2    pages: 2883-2892     authors: Abraham, Acree
```

Confirmed a second way against PMC6648601 itself, which returns the same record.

**The DOI and the PMC id are correct and point at a real Abraham & Acree paper.
Only the volume/issue/page string is wrong.** Correct: **ACS Omega
2019;4(2):2883–2892**.

### D2 — that paper's Table 1 does not contain the number attributed to it

Every site says both coefficients were *"read VERBATIM from the open-access
**Table 1**"* of that paper. Its Table 1 is captioned *"Coefficients in Equation 3
for Water–Solvent Partitions; SP = log P"* and **contains no hexadecane row** —
it covers methanol, ethanol, the propanols, butanols, hexanol, acetonitrile,
DMSO and others. So `a(water → hexadecane) = −3.587` **cannot have come from
there.** Verified directly against the article.

This is the serious one. D1 is a transcription slip; D2 means the stated source
does not contain the stated number.

### D3 — the cross-check reference is attributed to the wrong journal

Cited as *"cross-checked against Table 1 of **J Comput Aided Mol Des 2017**
(PMC5579141)"*. PMC5579141 is:

> Abraham MH, Acree WE Jr. "Descriptors for Pentane-2,4-dione and Its
> Derivatives." ***J Solution Chem*** 2017;**46(8):1625–1638**.
> doi:10.1007/s10953-017-0667-y

**But the cross-check itself is sound.** That paper's Table 1 — captioned
"Coefficients in Eqs. [1] and [2] for partition of solutes from water and the gas
phase to **wet** organic solvents at 298 K" — **does** carry both numbers:
octan-1-ol `a = 0.034` and hexadecane `a = −3.587`. Only the journal name
attached to the PMC id is wrong.

### D4 — a cross-check that does not reproduce

Three sites additionally cross-check against *"the open-access Acree/Abraham
compilation (hexadecane −3.587, cyclohexane −3.740, **hexane −3.578, heptane
−3.571**)"*. Read off that compilation's own table (below):

| solvent | quoted in the repo | actually printed | reproduces? |
|---|---:|---:|---|
| hexadecane | −3.587 | **−3.587** | ✔ |
| cyclohexane | −3.740 | **−3.740** | ✔ |
| hexane | −3.578 | **−3.599** | ✘ |
| heptane | −3.571 | **−3.317** | ✘ |

**Nothing numeric depends on hexane or heptane** — only hexadecane and wet
octanol enter the calculation — which is precisely why this would have survived
indefinitely. A cross-check that does not hold is worse than no cross-check,
because it reads as corroboration.

**Two readings, and both are defects.** Either (a) the two values were
mis-transcribed from the chapter identified below, or (b) *"the open-access
Acree/Abraham compilation"* means some other compilation in which hexane really
is −3.578. Reading (b) cannot be checked, because the phrase names no title, no
year and no identifier — **an unnamed source is not a source**. The chapter
identified below is the only Acree/Abraham compilation I could retrieve in which
hexadecane is −3.587 and cyclohexane is −3.740, i.e. in which the two values that
*do* reproduce, reproduce. Its table appears twice in the retrieved chapter text
with identical figures.

---

## The source that does support both numbers, and is open

**Acree WE Jr, Grubbs LM, Abraham MH. "Prediction of Partition Coefficients and
Permeability of Drug Molecules in Biological Systems with Abraham Model Solute
Descriptors Derived from Measured Solubilities and Water-to-Organic Solvent
Partition Coefficients." In: *Toxicity and Drug Testing*, IntechOpen, 2012.**
Open access.

Its wet-solvent water-to-organic-solvent log P coefficient table
(`log P = c + eE + sS + aA + bB + vV`), read off the retrieved chapter text:

```
Wet Solvent                    c        e        s        a        b       v
Octan-1-ol (wet)             0.088    0.562   -1.054    0.034   -3.460   3.814
Hexadecane                   0.087    0.667   -1.617   -3.587   -4.869   4.433
Cyclohexane                  0.159    0.784   -1.678   -3.740   -4.929   4.577
Hexane                       0.361    0.579   -1.723   -3.599   -4.764   4.344
Heptane                      0.325    0.670   -2.061   -3.317   -4.733   4.543
```

`0.034 − (−3.587) = 3.621` ✔, and `3.621 × 0.37 = 1.340` ✔.

**Two independent open sources now support both coefficients**: this chapter, and
PMC5579141 correctly named (D3). Either is a valid citation; quoting both is
better, since that is what the original text was trying to do.

---

## Exact replacement text, per site

### FROZEN — do not edit. Erratum only.

Every wrong string is quoted here **verbatim and in full**, with its exact JSON
key path, so that this erratum stands in for the edit that must not happen. Six
strings across three files, not three — the efflux registration carries the
defect in three separate list entries, and quoting only the first would leave
D3 and D4 unrecorded at the site where they occur.

| file | JSON key path | wrong string, verbatim | defects |
|---|---|---|---|
| `results/sweep/tworoute_promotion_prereg.json` | `parameters_changed.hbd_logp_decrement.source` | `Abraham MH, Acree WE Jr. ACS Omega 2019;4(7):12920-12926, doi:10.1021/acsomega.8b03242, PMC6648601 -- Table 1, read VERBATIM. Cross-checked against Table 1 of J Comput Aided Mol Des 2017 (PMC5579141) and the open-access Acree/Abraham compilation (hexadecane -3.587, cyclohexane -3.740, hexane -3.578, heptane -3.571). Measured A per donor group from the open-access tables of PMC5830486 and PMC3477586.` | D1, D2, D3, D4 |
| `results/sweep/tworoute_run_prereg.json` | `prereg.promotion_prereg.parameters_changed.hbd_logp_decrement.source` | byte-identical to the string above (verified by comparison, not by assumption) | D1, D2, D3, D4 |
| `results/sweep/tworoute_efflux_prereg.json` | `quantities_added.hbd_logp_decrement_sourced.coefficient_sources[0]` | `Abraham MH, Acree WE Jr. ACS Omega 2019;4(7):12920-12926, doi:10.1021/acsomega.8b03242, PMC6648601, Table 1 -- read VERBATIM, open access.` | D1, D2 |
| `results/sweep/tworoute_efflux_prereg.json` | `quantities_added.hbd_logp_decrement_sourced.coefficient_sources[1]` | `J Comput Aided Mol Des 2017, PMC5579141, Table 1 -- cross-check, read VERBATIM, open access.` | D3 |
| `results/sweep/tworoute_efflux_prereg.json` | `quantities_added.hbd_logp_decrement_sourced.coefficient_sources[2]` | `Acree/Abraham open-access compilation -- further cross-check: hexadecane -3.587, cyclohexane -3.740, hexane -3.578, heptane -3.571.` | D4 |

**These must not be rewritten.** A pre-registration's whole value is that it
cannot be revised after the fact; editing one to fix a citation destroys that
property for a gain this file already provides. The record is: the *coefficients
those files registered are correct and independently verifiable*; the *source
string attached to them is wrong in the four ways above*.

**Not a defect, and listed so it is not swept up by a search.**
`results/sweep/tworoute_efflux_prereg.json`, key
`explicitly_NOT_used_because_UNVERIFIED.clogpalk_group_increments.citation`,
reads `Kenny PW, Montanari CA, Prokopczyk IM. J Comput Aided Mol Des 2013;27:389-402.`
That one **is** *J Comput Aided Mol Des*, it is a different paper, and it is
correct. The same is true of `results/two_route_efflux.md` line 471. D3 concerns
only the 2017 / PMC5579141 attribution.

### NOT frozen — safe to correct, and here is the text

**Correct citation block** (the replacement, in prose form):

> Both coefficients read verbatim from the open-access coefficient table of
> Acree WE Jr, Grubbs LM, Abraham MH, "Prediction of Partition Coefficients and
> Permeability of Drug Molecules in Biological Systems with Abraham Model Solute
> Descriptors…", in *Toxicity and Drug Testing*, IntechOpen 2012 (wet-solvent
> water-to-organic log P table). Independently confirmed in Table 1 of Abraham MH,
> Acree WE Jr, "Descriptors for Pentane-2,4-dione and Its Derivatives",
> *J Solution Chem* 2017;46(8):1625–1638, doi:10.1007/s10953-017-0667-y,
> PMC5579141. Measured `A` per donor group from the open-access tables of
> PMC5830486 and PMC3477586.
>
> NOT the source: Abraham MH, Acree WE Jr, *ACS Omega* 2019;4(2):2883–2892,
> doi:10.1021/acsomega.8b03242, PMC6648601. Earlier revisions of this file cited
> it, with the wrong volume/issue/pages, as the source of both coefficients. Its
> Table 1 contains no hexadecane row. See `results/abraham_citation_erratum.md`.

Line numbers below are **as they stood before the correction** (2026-08-21) and
are recorded for provenance only; they have since shifted, and the block each
one names is the durable identifier. All four have been applied — see **Status**
above.

| file | lines (pre-correction) | what to replace |
|---|---|---|
| `experiments/two_route_efflux.py` | 85–88 | the `# VERIFIED, verbatim, in the open-access coefficient tables of / # Abraham MH, Acree WE Jr. ACS Omega 2019;4(7):12920-12926, / # doi:10.1021/acsomega.8b03242, PMC6648601, Table 1; / # cross-checked against the Table 1 of J Comput Aided Mol Des 2017, PMC5579141.` block |
| `experiments/two_route_efflux.py` | 396–399 | `Both coefficients read VERBATIM from the open-access Table 1 of Abraham & / Acree, ACS Omega 2019;4(7):12920-12926, doi:10.1021/acsomega.8b03242 / (PMC6648601), cross-checked against Table 1 of J Comput Aided Mol Des 2017 / (PMC5579141).` |
| `results/two_route_efflux.md` | 415–421 | the block quote beginning `> Both read verbatim from the open-access Table 1 of **Abraham MH, Acree WE Jr.**` — note this file was under concurrent edit by another agent at the time of writing, so re-read before patching |

### ROUTED AND APPLIED — was outside the investigation's territory

`tmesim/molecular/two_route.py`, **lines 206–210 as they stood before the
correction**. Text that was there, verbatim:

```
#: Both coefficients read VERBATIM from the open-access Table 1 of Abraham MH,
#: Acree WE Jr, ACS Omega 2019;4(7):12920-12926, doi:10.1021/acsomega.8b03242
#: (PMC6648601); cross-checked against Table 1 of J Comput Aided Mol Des 2017
#: (PMC5579141) and the open-access Acree/Abraham compilation (hexadecane
#: -3.587, cyclohexane -3.740, hexane -3.578, heptane -3.571).
```

Replacement:

```
#: Both coefficients read VERBATIM from the open-access wet-solvent
#: water-to-organic log P coefficient table of Acree WE Jr, Grubbs LM,
#: Abraham MH, in Toxicity and Drug Testing, IntechOpen 2012 (octan-1-ol wet
#: a = 0.034; hexadecane a = -3.587; cyclohexane a = -3.740). Independently
#: confirmed in Table 1 of Abraham MH, Acree WE Jr, "Descriptors for
#: Pentane-2,4-dione and Its Derivatives", J Solution Chem 2017;46(8):1625-1638,
#: doi:10.1007/s10953-017-0667-y, PMC5579141, which carries both values.
#: NOT the source, despite earlier revisions of this file: ACS Omega
#: 2019;4(2):2883-2892, doi:10.1021/acsomega.8b03242, PMC6648601 -- its Table 1
#: has no hexadecane row. See results/abraham_citation_erratum.md.
```

**No assertion is affected.** `tests/test_two_route.py` pins the *value* `1.34`
and its arithmetic (`3.621 × 0.37`), never the citation string, so this
correction moves no test.

---

## What this does and does not change

* **`3.621` stands.** Verified against two independent open sources.
* **`1.34` stands.** `3.621 × 0.37 = 1.340`, and 0.37 is ethanol's measured `A`.
* **`HBD_LOGP_DECREMENT` is unchanged and was not touched.**
* **The `NumHDonors` pairing is unchanged and is correct** — see
  `results/hbd_convention_audit.md` §2, which is what this file fell out of.
* What changes is that a reader can now actually find the numbers. Before this
  file, following the citation as written led to a paper whose Table 1 does not
  contain one of the two coefficients being attributed to it.
