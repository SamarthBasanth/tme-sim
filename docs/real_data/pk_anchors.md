# Real measured tissue-PK anchors (multi-drug compendium)

**Purpose.** The physics engine's penetration/tissue-concentration predictions
were previously checked against a *single* drug (doxorubicin ~40-50 µm,
`results/pk_reality.md`). This compendium collects **real, published,
quantitative** tissue-pharmacokinetics measurements for **multiple drugs** so the
transport engine can be **scored** against reality more broadly. (Scored, not
"validated": these anchors are a test the engine can fail, and some of them it
does — `results/CLAIM_INVENTORY.md` §1e.) Every number is
traced to a peer-reviewed paper and its PubMed ID was verified against PubMed.

**Attribution.** The numbers below were retrieved and verified using PubMed.
Machine-readable form: `tmesim/validation/pk_anchors.py` (`PK_ANCHORS`,
`anchors_for(drug)`, plus `penetration_anchors()` / `partition_anchors()`
accessors). That module is pure data — it does not import or run the engine.

---

## Two distinct quantities (do not conflate)

The engine produces two kinds of output that map to two *different* classes of
real measurement:

| Sim output | Real measurement class | Anchors here |
|---|---|---|
| `penetration_length_um` (L) and single-vessel `1/e` distance | **spatial perivascular gradient** — how far a drug reaches from a vessel | doxorubicin, mitoxantrone, topotecan (perivascular half-distance) |
| `DrugParams.bbb_penetrance` | **partition ratio** — tissue-or-CSF : plasma concentration | temozolomide, methotrexate, cisplatin (brain/CSF:plasma) |

A brain:plasma ratio is **not** a penetration distance. Only the perivascular
anchors are directly comparable to L / the `1/e` distance; the ratio anchors are
comparable to the structural `bbb_penetrance` prior.

### Half-distance vs 1/e distance (required conversion)

The perivascular anchors report the distance `d_half` at which drug intensity
falls to **half** its perivascular value. The sim reports the **e-folding**
length L (fall to `1/e` = 37%) and the `1/e` distance. For an exponential
profile `C(r) ~ exp(-r/L)`:

```
d_half   = L * ln2  ≈ 0.693 * L        =>   L = d_half / ln2 ≈ 1.443 * d_half
d_(1/e)  = L
```

So compare the sim's L / `1/e` numbers to the **L-equivalent** column below, not
to the raw half-distances. (`pk_reality.py` currently compares to the raw
40-50 µm, which is slightly conservative — it treats the half-distance as the
1/e distance, understating the L-equivalent anchor a bit; the correct
L-equivalent ~58-72 µm does not change the verdict it recorded at the time.)

**Superseded, 2026-09-03.** The parenthesis above carried a large predicted `L` for
doxorubicin and an order-of-magnitude FAIL. That is no longer what the engine produces:
`results/pk_multidrug.md` is the later dated record and reports doxorubicin inside its
cited band, with the two out-of-sample drugs 2/2 inside theirs. Read that file for the
current numbers; they are deliberately not copied here, because a second copy is how this
one went stale.

---

## The anchors

### Perivascular penetration gradients (comparable to sim L / `1/e`)

| Drug | Half-distance (µm) | L-equivalent (µm) | Model | Method | Citation | PMID | Conf. |
|---|---|---|---|---|---|---|---|
| doxorubicin | 40–50 | 58–72 | murine 16C, EMT6 + human PC-3 xenografts (mice) | quantitative immunofluorescence vs CD31+ vessels | Primeau et al., Clin Cancer Res 2005;11(24 Pt 1):8782–8 | 16361566 | high |
| doxorubicin | 25–75 (10 min) | 36–108 | human tumor xenografts (mice) | quantitative immunofluorescence | Patel, Trédan, Tannock, Cancer Chemother Pharmacol 2013;72(1):127–38 | 23680920 | high |
| mitoxantrone | 25–75 (10 min) | 36–108 | human tumor xenografts (mice) | quantitative immunofluorescence | Patel, Trédan, Tannock 2013 | 23680920 | medium |
| topotecan | 25–75 (10 min) | 36–108 | human tumor xenografts (mice) | quantitative immunofluorescence | Patel, Trédan, Tannock 2013 | 23680920 | medium |

DOIs: [10.1158/1078-0432.CCR-05-1664](https://doi.org/10.1158/1078-0432.CCR-05-1664),
[10.1007/s00280-013-2176-z](https://doi.org/10.1007/s00280-013-2176-z).

Notes:
- Primeau 2005 is the canonical anchor already used in `pk_reality.py`:
  doxorubicin fluorescence decreases exponentially with distance from tumor
  vessels, falling to half its perivascular concentration at ~40–50 µm; regions
  of hypoxia lie 90–140 µm from vessels, i.e. beyond the drug's reach.
- Patel 2013 independently reproduces the doxorubicin gradient and reports the
  ~25–75 µm half-distance as a **shared range for all three drugs at 10 min**
  (distribution improves by 3–24 h but a significant gradient persists). The
  mitoxantrone/topotecan values are therefore order-of-magnitude, not
  drug-specific point values — hence `medium` confidence.

### Brain/CSF : plasma partition ratios (comparable to `bbb_penetrance`)

| Drug | Ratio | Compartment | Model | Method | Citation | PMID | Conf. |
|---|---|---|---|---|---|---|---|
| temozolomide | 0.178 | brain interstitium : plasma (AUC) | human peritumoral brain, 150 mg/m² oral | intracerebral microdialysis + LC-MS/MS | Portnow et al., Clin Cancer Res 2009;15(22):7092–8 | 19861433 | high |
| methotrexate | 0.013 | CSF : plasma (steady state) | children with ALL, no CNS disease, high-dose IV | paired CSF/plasma antifolate assay | Morse et al., J Clin Oncol 1985;3(1):19–24 | 3855310 | high |
| cisplatin | 0.025–0.029 | CSF : plasma (peak total Pt) | solid metastatic brain tumor, IV 60/100 mg/m² | atomic absorption spectroscopy, Ommaya CSF | Nakagawa et al., Cancer Chemother Pharmacol 1996;37(4):317–26 | 8548876 | medium |

DOIs: [10.1158/1078-0432.CCR-09-1349](https://doi.org/10.1158/1078-0432.CCR-09-1349),
[10.1200/JCO.1985.3.1.19](https://doi.org/10.1200/JCO.1985.3.1.19),
[10.1007/s002800050391](https://doi.org/10.1007/s002800050391).

Notes:
- **Temozolomide** — mean brain-ECF:plasma AUC ratio **17.8%** (the widely-cited
  "~0.2" TMZ brain-penetration figure); mean peak brain ECF 0.6 ± 0.3 µg/mL,
  time-to-peak ~2 h. Directly measured in human peritumoral brain by
  microdialysis. High confidence.
- **Methotrexate** — mean steady-state CSF:plasma **0.013** (SD 0.01) in patients
  *without* meningeal disease; rose ~12× to 0.157 with overt CNS leukemia. Use
  the 0.013 value as the intact-barrier partition anchor. A classic example of a
  poorly-CNS-penetrant small molecule.
- **Cisplatin** — intact blood-CSF barrier (IV, solid metastases): peak
  CSF:plasma total-Pt ~**2.5%** (60 mg/m²) to ~**2.9%** (100 mg/m²). The same
  paper reports much higher ratios (up to ~40%) for disrupted-barrier /
  intra-arterial / post-operative-cavity conditions — those are **not** the
  systemic-partition anchor and are excluded. `medium` confidence (total Pt,
  Ommaya sampling, small n).

---

## Supporting references (context / directional, not point anchors)

- **Minchinton & Tannock, Nat Rev Cancer 2006;6(8):583–92** (PMID 16862189,
  [DOI](https://doi.org/10.1038/nrc1893)) — review establishing limited tissue
  penetration as a general chemoresistance mechanism.
- **Tannock, Lee, Tunggal, Cowan, Egorin, Clin Cancer Res 2002;8(3):878–84**
  (PMID 11895922) — multicellular-layer (~200 µm) method: cisplatin, etoposide,
  gemcitabine, paclitaxel and vinblastine all penetrate slowly. Directional
  support (slow penetration), not a clean distance.
- **Fujimori, Covell, Fletcher, Weinstein, J Nucl Med 1990;31(7):1191–8**
  (PMID 2362198) — the "binding-site barrier": antibodies penetrate poorly. The
  directional anchor for the sim's small-molecule-vs-antibody ordering check.

---

## Which anchors validate what

**Directly usable for the sim's penetration/tissue-concentration check (multi-drug):**

1. **doxorubicin** — half-distance 40–50 µm → L-equivalent **58–72 µm**
   (Primeau); independently 25–75 µm → 36–108 µm (Patel). The primary
   quantitative penetration-length target.
2. **mitoxantrone**, **topotecan** — L-equivalent **36–108 µm** (Patel). Extend
   the penetration-length check to two more drugs (order-of-magnitude
   confidence). Also enable a **same-decade** check: all three are large,
   strongly-bound anthracycline-class/topoisomerase agents and should each show
   short (tens-of-µm) penetration — none should come out at millimetres. (This was
   labelled a relative-ordering check until 2026-09-03. Its content was never an
   ordering — it is a magnitude band applied to each compound separately — and the
   label is corrected rather than the check.)

These four rows give a real, multi-drug penetration-length band (~36–108 µm,
tens of µm) that the engine's `penetration_length_um` / single-vessel `1/e`
distance must reproduce. This paragraph recorded a large over-prediction as the
then-current state; `results/pk_multidrug.md` supersedes it, and the compendium's
value now is breadth — the same check on more compounds — rather than sharpening a
failure that has since moved.

**Usable for the sim's `bbb_penetrance` / partition check (not penetration length):**

3. **temozolomide** 0.178, **methotrexate** 0.013, **cisplatin** ~0.025–0.029 —
   real brain/CSF:plasma ratios spanning ~1.3%→18%, i.e. a >10× spread. Compare
   each drug's measured ratio to its own compiled `DrugParams.bbb_penetrance`,
   one compound at a time, against the ~0.01–0.2 magnitude band. These do
   **not** validate penetration length.

   **This item used to print the three compounds as an explicit sequence and
   pre-register that sequence as the acceptance criterion.** Both halves are gone
   as of 2026-09-03. The sequence is the output withdrawn on 2026-08-21 and a
   pre-registration of it is a scheduled re-assertion; the prohibition covers a
   criterion as squarely as a result. The criterion had also **failed on its own
   terms** — `results/pk_multidrug.md:42` scores that concordance at Spearman
   ρ = 0.50 (Kendall τ = 0.33, 2 of 3 pairs concordant), a partial agreement this
   file never recorded. Per-compound magnitude against a published band is the
   check that survives, and it is the one stated above.

## Caveats (honest)

- **Half-distance ≠ 1/e distance.** Apply the `L = d_half / ln2` conversion
  (`pk_anchors.half_distance_to_efold_length`) before comparing to the sim's L;
  the raw µm numbers are half-distances.
- **Partition ratio ≠ penetration distance.** The temozolomide/methotrexate/
  cisplatin numbers are whole-compartment concentration ratios (in-vivo AUC or
  peak), reflecting BBB/BCB transport + clearance, not a spatial gradient. They
  belong to the `bbb_penetrance` check only.
- **Model heterogeneity.** Perivascular anchors are mouse xenografts (immuno-
  fluorescence); partition anchors mix human brain-tumor patients (TMZ, cisplatin)
  and pediatric leukemia CSF (MTX). Absolute values depend on drug, dose, tumor
  type, and barrier integrity — treat ranges as order-of-magnitude targets.
- **Confidence flags.** `high` = a drug-specific measured point/range from a
  primary study; `medium` = a value reported jointly for several drugs, or with
  small n / surrogate compartment (Ommaya CSF, total Pt). None are estimates —
  all are measured — but `medium` entries should not be over-interpreted as
  precise single-drug numbers.
