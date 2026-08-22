# Pruijn 2008 — multicellular-layer diffusion coefficients (`Dmcl`)

## Headline finding, read this first

**The paper does NOT publish per-compound `Dmcl` values for SiHa, HCT8-Ea or
HCT8-Ra.** For the three non-HT29 cell lines the only measured, per-compound,
per-line number in the entire article + Accessory Publication is **urea**
(`pruijn2008_per_cellline_urea_and_thickness.csv`). Everything else about those
lines is published only as **two fitted coefficients** (`a` and `c` of Eqn 1) in
Table 1.

The nine-compound Study B cassette *was* run in all four lines, but only the
**HT29 arm** is tabulated (the `Dmcl_measured_studyB_*` column of Table S1). The
SiHa / HCT8-Ea / HCT8-Ra arms appear nowhere as numbers — not in Table S1, not in
a figure of raw points. Fig. 4 is HT29-only (its caption says so). Fig. 5 is
*simulations* from Eqn 1, not data. Fig. 3a plots urea **flux traces** for the
four lines but is explicitly "only one of several replicates ... for clarity".

**Consequence for the out-of-sample plan:** a per-compound, cross-cell-line
prediction test cannot be built from this paper as published. What *can* be built
is (i) a 4-point urea test, and (ii) a test against the published `a`/`c`
coefficients. Both are described under "What test is actually possible" below.
The underlying per-compound data may exist with the authors — the Accessory
Publication note says the compound set "is available from the author".

---

## Full citation

> Pruijn FB, Patel K, Hay MP, Wilson WR, Hicks KO. "Prediction of Tumour Tissue
> Diffusion Coefficients of Hypoxia-Activated Prodrugs from Physicochemical
> Parameters." *Australian Journal of Chemistry* 2008;61(9):687–693.
> doi:10.1071/CH08240

Corresponding author email as printed on the Accessory Publication:
`f.pruijn@auckland.ac.nz`.

## Access route

The publisher (CSIRO Publishing) returns HTTP 403 to automated retrieval — a
Cloudflare bot challenge, not a paywall. Both PDFs were retrieved **HTTP 200**
from the Internet Archive Wayback Machine, which holds snapshots of the
publisher's own open-access PDF URLs:

| file | URL | result |
|---|---|---|
| Article (7 pp.) | `http://web.archive.org/web/20240530035922id_/https://www.publish.csiro.au/CH/pdf/CH08240` | 200, 442 258 B |
| Accessory Publication / Table S1 (13 pp.) | `http://web.archive.org/web/2015/https://www.publish.csiro.au/?act=view_file&file_id=CH08240_AC.pdf` | 200, 181 692 B |

These are the publisher's own files, archived. **This is a primary source, not a
reproduction in a later paper.**

Extraction: `pdfplumber` 0.11.8 word-coordinate column assignment (not text-flow
parsing, which mis-assigns the wrapped cells), cross-checked against an
independent `pdftotext -layout` rendering of the same PDFs.

PubMed/PMC were not the route: this record is not indexed in PubMed or PMC.

---

## File 1 — `pruijn2008_tableS1.csv` (71 compounds, HT29 only)

Source: **Accessory Publication, Table S1**, PDF pages 1–12 (page 13 is
footnotes). Caption, verbatim:

> "Table S1. Structures, physicochemical parameters, calculated Dmcl (from Eq. 1
> using the coefficients for HT29 in Study A as listed in Table 1) and measured
> Dmcl for HT29 MCLs. Errors are SEM and units of Dmcl are cm2s-1. The
> calculated Dmcl values were based on measured logP7.4 when available."

**Every row of this table is HT29.** Both the "Study A" and "Study B" measured
columns are HT29 measurements; Study A and Study B are two *experimental
campaigns*, not two cell lines.

### What the measured quantity physically is

`Dmcl` = the **effective diffusion coefficient of free (unbound) drug through a
multicellular layer**. MCLs are grown on a collagen-coated porous Teflon support
membrane and mounted in a diffusion chamber between two **stirred** compartments.
Drug enters the donor compartment; both compartments are sampled through the gas
outlet. Concentration–time data in *both* compartments are fitted to
reaction–diffusion models (Fick), and `Dmcl` is the fitted transport coefficient.
It is **not** a penetration depth, a % penetrated, or an AUC ratio.

The authors are explicit that **Eqn 1 (and hence `Dmcl`) applies to free
(unbound) drug only**; where high-affinity binding sites exist (e.g. DNA
intercalators) flux is slowed during an extended pre-steady-state phase that is
modelled separately, not folded into `Dmcl`.

### Units — NOT converted, read this before using the numbers

Table S1 tabulates the column as **`Dmcl x 10^6`**, in **cm² s⁻¹**. A tabulated
value *v* therefore means `Dmcl = v x 10^-6 cm² s⁻¹`.

The CSV stores *v* **exactly as printed**. No conversion has been applied.

If you want µm² s⁻¹ the conversion would be `x 100` (since 1 cm² s⁻¹ = 1e8 µm²
s⁻¹, so `v x 1e-6 cm²/s = v x 100 µm²/s`). **That conversion has deliberately not
been applied here.** The repo's existing `TABLE_UNITS_TO_UM2_S = 100.0` in
`experiments/diffusivity_falsification.py` encodes the same factor, and it is
correct.

> **Unit inconsistency inside the paper — flagged, not harmonised.** The
> Experimental section reports the per-cell-line urea values in
> **cm² s⁻¹ × 10⁷**, whereas Table S1 uses **× 10⁶**. They do agree for HT29:
> `4.50e-7` (Experimental) = `0.450e-6` (Table S1, urea Study A). File 3 stores
> the Experimental values in the paper's own `× 10⁷` units. **Do not join files 1
> and 3 without applying the factor of 10 yourself, consciously.**

### Experimental conditions (HT29, Study A)

* Growth medium: αMEM + 10% fetal bovine serum. Donor compartment 5.5 mL αMEM,
  **with** 10% FBS in Study A, **without** in Study B.
* Gas-phase O₂ raised to **95%** during flux experiments, to extinguish central
  hypoxia in the MCL and so minimise reductive metabolism of the prodrugs.
* Initial donor concentration: **(10–50) × 10⁻⁶ M** (Study A);
  **10 × 10⁻⁶ M** each (Study B cassette, with stated exceptions).
* `~1 × 10⁻⁶ M` ¹⁴C-urea co-added to **every** experiment as an internal standard
  for the thickness of each individual MCL.
* Timepoint: flux is a time course, not a single timepoint. Fig. 3a runs
  **0–6 h**.
* Support membrane porosity 11%; flux through bare support membranes was measured
  for each compound; the Teflon support measurement accounts for the unstirred
  water layers, **assumed** similar in the presence of MCLs.
* Assay: HPLC with photodiode-array detection (Study A TPZ analogues); ion-trap
  MS (Study B, seven compounds); scintillation counting (³H-mannitol, ¹⁴C-urea).
* **Temperature is not stated in the retrieved article text.** MCL growth and
  flux methodology is cited out to refs 35 and 38 (Hicks et al.). *No temperature
  has been assumed or recorded here.*
* MCL thickness: measured on H&E-stained frozen sections. Per-line estimates in
  file 3 are for the specific MCLs of Fig. 3a only.

### Columns

`compound_id`, `ref_footnote` (Table S1 "Ref" column → synthesis paper, footnote
list on PDF p.13), `mw_daltons`, `h_bond_donors`, `h_bond_acceptors`,
`logP74_calculated` (+`_sem`), `logP74_measured` (+`_sem`),
`Dmcl_calculated_eq1_x1e6_cm2_s`, `Dmcl_measured_studyA_x1e6_cm2_s` (+`_sem`),
`Dmcl_measured_studyB_x1e6_cm2_s`, `source_pdf_page`, `notes`.

* `logP74_measured` / `logP74_calculated` are **the paper's own columns**, kept
  separate. The repo's existing literal collapses them into one field
  (measured-if-available, else calculated); this CSV does not. **These are
  measured-by-the-paper values — the shipping model computes its own descriptors,
  so do not silently mix the two.**
* `Dmcl_calculated_eq1_*` is the **authors' own Eqn-1 prediction**, not a
  measurement. Never mix it with the measured columns.
* `Dmcl_measured_studyB_*` is populated for **9 of 71** rows.

### Coverage and what is uncertain

* **71 of 71 compounds obtained.** Nothing is missing.
* **9 of 9** Study B (HT29 cassette) values obtained.
* `logP74_calculated` is **empty and flagged UNCERTAIN for SN28860 and SN28864**.
  Those two cells render in the source PDF as a single stray glyph (`3` and `6`
  respectively) with no recoverable remainder. Both compounds *do* have a
  measured logP7.4 (2.15 ± 0.02 and 0.001 ± 0.012), so nothing load-bearing is
  lost. **No value was guessed.**
* `logP74_calculated` is empty for 4 rows printed `np` ("not published",
  footnote e): SN28870, SN28871, SN28874, SN28875. **3 compounds (SN28871,
  SN28874, SN28875) carry no logP7.4 at all**, measured or calculated.
* 19 of 71 rows have no *measured* logP7.4.
* 5 cells were **line-wrapped in the PDF** (minus sign or SEM on an adjacent
  line) and were reassembled from the layout, then confirmed against the
  `pdftotext -layout` rendering: SN28400 (−0.187 ± 0.027), SN28565
  (0.148 ± 0.004), SN28620 (−0.471 ± 0.011), SN28823 (−0.854 ± 0.015), SN29143
  (−0.071 ± 0.021), all in `logP74_measured`. Flagged in `notes`.
* 3 rows carry the source's own `***` marker, whose legend on PDF p.12 reads
  "*** in table but not in graph": SN28823, SN29650, SN29662. Flagged in `notes`.
* SN28246, SN28287 and NLCQ-1 have a measured logP7.4 printed **without a SEM**
  (footnotes f / g). Flagged in `notes`; the SEM field is left empty, not zero.
* Compound identity: TPZ is printed as `24280 (TPZ)` and is recorded
  `SN24280 (TPZ)`.

### Compound-count discrepancy in the source — unresolved, reported as-is

The paper's own counts do not fully reconcile and it does **not** say which
compounds were dropped:

| where | count |
|---|---|
| Table S1 (rows actually listed) | **71** |
| Body text, Study A | "67 TPZ analogues" (+ NLCQ-1 + 3 flux markers = 71) |
| Table 1 header | "Study A (**69** individual compounds)" |
| Fig. 4 caption | "65 TPZ analogues and four other compounds" = **69** |
| Rows with `***` (in table, not in graph) | **3** |
| Rows with no logP7.4 at all | **3** |

71 − 3 = 68, not 69. **This is left unresolved.** Do not assume which two
compounds the n = 69 regression excluded.

### Structural identifiers — a real limitation

**Table S1 contains no SMILES and no CAS numbers.** Structures are supplied only
as **embedded 2-D drawings**. The machine-readable identifiers available are:

1. the Auckland Cancer Society Research Centre **`SN` code** (e.g. `SN28682`);
2. the **`ref_footnote`** number, pointing to the synthesis/characterisation
   paper on PDF p.13 (12 references, captured in the CSV as `ref_footnote`;
   `np` = not published).

Only the four non-TPZ compounds (urea, mannitol, acridan-9(10H)-one, NLCQ-1) have
names that resolve directly to a structure. **The 67 SN-coded compounds cannot be
resolved to SMILES from this file alone** — that requires going to the cited
synthesis papers or OCR-ing the drawings. Neither was done, and no structure has
been guessed.

---

## File 2 — `pruijn2008_table1_eqn1_coefficients.csv`

Source: **Table 1, article p.690.** Caption: "Coefficients in Eqn 1, which
describes the dependence of diffusion coefficient in MCLs (Dmcl) on
physicochemical parameters. All coefficients were significant (P < 0.002)."

### Eqn 1, as printed

```
log(Dmcl) = a + b·log(MW)  +  ───────────────────────────────────────────────
                                        ⎛ logP7.4 − x + y·HD + z·HA ⎞
                                1 + exp ⎜ ───────────────────────── ⎟
                                        ⎝             w             ⎠
```
with numerator `c` on the second term. `a` and `c` are the coefficients that
determine the **low** and **high** logP7.4 asymptotes respectively.

Study A fit quality (HT29): **R² = 0.907**, **Q² = 0.879** (7-fold
leave-group-out cross-validation).

### The per-line parameters — this is the crux

The authors fitted the other cell lines by **holding every coefficient at the
HT29 Study A value except `a` and `c`**. Verbatim:

> "We, therefore, fitted Eqn 1 to each dataset, holding all coefficients at the
> values determined for HT29 in Study A, except for a and c which are the
> coefficients that determine the high and low log P7.4 asymptotes."

So exactly **two** numbers differ per cell line, and **both are fitted, not
measured.**

| coefficient | Study A HT29 | Study B HT29 | SiHa | HCT8-Ea | HCT8-Ra |
|---|---|---|---|---|---|
| `a` | −5.20 ± 0.33 | −5.35 ± 0.37 | −4.80 ± 0.37 | −5.74 ± 0.43 | −4.39 ± 0.43 |
| `c` | 1.149 ± 0.092 | 1.149 ± 0.092 | 1.149 ± 0.092 | 1.99 ± 0.24 | 0.611 ± 0.24 |
| `b`,`w`,`x`,`y`,`z` | fitted | – | – | – | – |

`–` is printed in the source and means "held at the Study A HT29 value".

> **Observation, flagged as inference and NOT as a claim by the authors:** the
> `c` printed for Study B HT29 and for SiHa is *identical to the Study A HT29
> value including its standard error* (1.149 ± 0.092). The arithmetic strongly
> suggests `c` was **not** re-fitted for those two, i.e. only `a` moved — despite
> the prose saying `a` and `c` were both free. The paper does not say this. The
> CSV marks those two cells `printed_in_table` rather than `fitted`; the other
> six are marked `fitted`. **Do not treat SiHa as having an independently fitted
> high-logP asymptote.**

### The "up-to-100-fold" figure — sourced precisely

The 100-fold is the authors' own, and it refers to **HCT8-Ea at the high-logP7.4
asymptote**, not to a spread in `a`:

> "The HCT8-Ra line gave much higher Dmcl values, especially at the low log P7.4
> asymptote, which presumably reflects the predominantly paracellular transport
> in these loosely packed MCLs. The epitheloid HCT8-Ea line gave a low log P7.4
> asymptote similar to HT29, but the largest increase at high log P7.4
> (100-fold)."

For calibration against that quote: the spread in `a` across the four Study B
lines is −5.74 → −4.39, i.e. **1.35 log₁₀ units ≈ 22×**, and the spread in `c` is
0.611 → 1.99. The "100-fold" is a statement about the fitted curve at high logP,
not about `a`. **Do not restate it as "the asymptote coefficients move 100-fold".**

### Why the authors say the lines differ

* Ranking of flux: **HCT8-Ra > SiHa > HT29 > HCT8-Ea**, "despite similar average
  MCL thicknesses".
* HCT8-Ea and HCT8-Ra are two clones of the **same parental HCT8 line**, with
  epithelial and rounded morphology respectively. HCT8 is hemizygous for
  α-E-catenin (`CTNNA1`) and segregates at high frequency to α-E-catenin-null
  round variants lacking adherens junctions. The authors note it is "surprising"
  that HCT8-Ra nonetheless expresses normal α-E-catenin yet still fails to form
  tight intercellular junctions.
* Across the four lines there is an **inverse correlation between urea flux and
  cellular packing density** (observed by microscopy, Fig. 2b), which the authors
  read as evidence that hydrophilic transport is **paracellular**, not
  transcellular.
* Proposed interpretation for HCT8-Ea: efficient adherens junctions occlude the
  paracellular route for hydrophilic compounds, while its mean packing density is
  relatively low compared with HT29 and SiHa.
* The authors' own caveat, verbatim: *"rigorous physiological interpretation of
  these apparent differences would require examination of larger numbers of
  compounds and quantitative assessment of cell packing, intercellular matrix
  composition, and intercellular junctions in these MCLs."*

**Note that packing density is reported only qualitatively (microscopy, Fig. 2b).
There is no tabulated per-line extracellular volume fraction or packing-density
number in this paper.**

---

## File 3 — `pruijn2008_per_cellline_urea_and_thickness.csv`

Source: **Experimental / "MCL Flux Studies", article p.692** (urea `Dmcl`), and
the **Fig. 3 caption, p.689** (thickness).

Urea `Dmcl`, verbatim: "The urea Dmcl values (mean ± s.e.m., [n], cm² s⁻¹ × 10⁷)
were 4.50 ± 0.23 [6] for HT29, 9.62 ± 0.85 [9] for SiHa, 1.31 ± 0.37 [15] for
HCT8-Ea, and 26.8 ± 0.72 [13] for HCT8-Ra."

**Units are × 10⁷ cm² s⁻¹ here — different from Table S1's × 10⁶.** Stored as
printed. Not converted. See the unit warning above.

Thickness caveat, verbatim from the Fig. 3 caption: "The thickness of these MCLs,
determined by fitting the independently measured Dmcl for urea (see Experimental)
were estimated as 163 µm (HT29), 245 µm (SiHa), 182 µm (HCT8-Ea), and 206 µm
(HCT8-Ra)." These are **estimates for the specific MCLs plotted in Fig. 3a**,
back-computed from the urea `Dmcl`. They are **not** independent per-line
thickness measurements and **must not** be used as such — in particular, using
them to predict urea `Dmcl` would be circular.

Measured spread in urea `Dmcl` across the four lines: **26.8 / 1.31 ≈ 20.5×**.

---

## What test is actually possible

Stated plainly so the scope is not overclaimed:

1. **A 4-point urea test.** Predict urea `Dmcl` in HT29, SiHa, HCT8-Ea, HCT8-Ra
   with unchanged constants. Real, measured, out-of-sample for 3 of 4 lines — but
   it is **one compound and four points**, and urea is a small hydrophilic
   paracellular marker, the least discriminating end of the range. Circularity
   warning: the Fig. 3a thicknesses are derived *from* these urea values.
2. **A coefficient-level test.** Fit the model's own predicted `Dmcl` vs logP7.4
   curve per line and compare the recovered `a`/`c` to Table 1. This tests the
   *shape*, not per-compound accuracy, and compares against **fitted** rather
   than measured quantities.
3. **The genuinely clean per-compound cross-line test is not available from the
   published record.** It needs the SiHa / HCT8-Ea / HCT8-Ra arms of the
   nine-compound cassette. Those were measured — Table 1's per-line fits could
   not exist otherwise — but they are not published. Obtaining them means
   **emailing the corresponding author** (`f.pruijn@auckland.ac.nz`, Auckland
   Cancer Society Research Centre); the Accessory Publication states the compound
   set "is available from the author". That is the single action that would
   unblock the real test.

---

## Verification performed

* **Independent re-extraction vs. the repo's existing literal.** The 71-row
  table already hardcoded at `experiments/diffusivity_falsification.py:92`
  (`PRUIJN_2008_TABLE_S1`) was cross-checked field-by-field against this
  extraction: MW, HD, HA, measured `Dmcl` (Study A), its SEM, the paper's
  calculated `Dmcl`, and logP7.4. **71 rows, 0 discrepancies.** The repo's
  literal is confirmed correct. It is, however, **missing the Study B column, the
  logP SEMs, the separate calculated/measured logP columns, and the `ref`
  column** — all of which are present here.
* **Study B membership cross-check.** The 9 Study B rows were identified purely
  from PDF word coordinates, then compared against the article's prose list
  ("six TPZ analogues (SN28205, SN28682, SN28870, SN29102, SN29347, and SN29366
  ...) and the lipophilic acridan-9(10H)-one ... two radiolabelled hydrophilic
  compounds (³H-mannitol, ¹⁴C-urea)"). **Exact match, 9/9** — an independent
  confirmation that the column was assigned to the right rows.
* **Cross-source unit check.** Table S1 urea Study A (`0.450 × 10⁻⁶ cm²/s`)
  equals the Experimental section's HT29 urea value (`4.50 × 10⁻⁷ cm²/s`),
  confirming both the extraction and the ×10⁶/×10⁷ relationship between the two
  parts of the paper.
* Every numeric cell was additionally spot-checked against an independent
  `pdftotext -layout` rendering.

## Standing rule

Nothing in the model has ever been fitted to the urea or coefficient data in
files 2 and 3, and nothing may be.

**The HT29 Table S1 data in file 1 is NOT an out-of-sample anchor.** The shipping
two-route transport law was selected against it (`results/two_route_promotion.md`,
`results/two_route_transport.md`). File 1 is retained for joining and provenance,
not as a validation target.
