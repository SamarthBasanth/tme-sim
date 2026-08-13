# Diffusivity falsification test — does the MW-only `D_eff` model survive contact with 71 measured tissue diffusion coefficients?

**Verdict up front: NO. The model is falsified, and not marginally.**

Measured `D_mcl` spans **29.7-fold** across 71 compounds measured by one method
in one cell line. The engine spans **2.4-fold** over the same compounds and
correlates with the measurement at **R² = 0.028**. Molecular weight — the
engine's only molecular input to diffusivity — explains **2.8 %** of the
variance in the measured values. Among four compounds of *identical* molecular
weight (263 Da) the measurement spans **3.7-fold** while the engine is flat by
construction.

Driver: `experiments/diffusivity_falsification.py`. **Nothing is registered in
`tmesim/validation/scorecard.py`, nothing in `tmesim/` was modified, and nothing
was fitted.** This is a one-shot falsification test.

---

## 1. What was being tested, and why it is decisive

`MolecularProfile.effective_diffusivity` — the `D_eff` that feeds the transport
PDE — is produced by this chain:

```
mw -> descriptors.hydrodynamic_radius_nm  = 0.066 * mw**(1/3)     [MW ONLY]
   -> stokes_einstein.free_diffusivity(radius, T)                 [radius ONLY]
   -> brownian.simulate_effective_diffusivity(free_D, ecm_frac)   [radius + ECM]
   -> ecm_mesh.ogston_renkin_hindrance(radius)                    [radius ONLY]
```

Every input is either molecular weight or the tissue ECM volume fraction.
`logP`, HBD, HBA, TPSA and formal charge are all computed by
`descriptors.compute_descriptors` and **none of them reaches
`effective_diffusivity`**.

This was not assumed — the script verifies it as a pre-flight, and the check is
in the run output:

| tirapazamine compiled via | `D_eff` |
|---|---:|
| `descriptors.from_smiles(SMILES)` → `build_profile` | 287.3216 µm²/s |
| `Molecule(descriptors={"mw": 178.151})` → `build_profile` | 287.3483 µm²/s |
| relative difference | **0.01 %** |

So at a fixed ECM fraction the engine **must** predict a near-constant `D_eff`
across compounds of similar MW. A measured dataset that spans a wide range at
matched MW therefore falsifies the model in one shot. That is the whole test.

---

## 2. The dataset, and exactly how it was obtained

> Pruijn FB, Patel K, Hay MP, Wilson WR, Hicks KO. "Prediction of Tumour Tissue
> Diffusion Coefficients of Hypoxia-Activated Prodrugs from Physicochemical
> Parameters." *Aust J Chem* 2008;61(9):687–693.
> [DOI](https://doi.org/10.1071/CH08240)

67 benzotriazine di-N-oxides (tirapazamine analogues) plus 4 miscellaneous
compounds (urea, mannitol, acridan-9(10*H*)-one, NLCQ-1) = **71**, all measured
in **HT29** multicellular layers, **one method, one cell line**. `D_mcl` is
obtained by fitting flux between two stirred compartments across the MCL to
Fick's second law — a tissue-scale effective diffusion coefficient in exactly
the engine's functional role. Commensurability is therefore not in question;
this is the same quantity in the same equation.

### 2.1 How the "publisher-blocked" paper was actually retrieved

`results/external_anchor_survey.md` §2.5 recorded this paper as **UNVERIFIED —
publisher-blocked (HTTP 403)** and listed it as the single highest-value
obtainable measurement. That assessment of its value was right. The 403
diagnosis was incomplete, and this is the part worth recording:

* CSIRO Publishing has **migrated** `publish.csiro.au` → `connectsci.au`.
* The 403 body is a **Cloudflare bot challenge** — literally
  *"Just a moment… Enable JavaScript and cookies to continue"* — **not a
  paywall**. It defeats any headless fetch regardless of headers.
* Unpaywall and OpenAlex both record this DOI as **bronze open access**, with
  `oa_status: "bronze"`, `is_oa: true`, and the publisher's own PDF as the OA
  location. The paper was never paywalled.

The content being open, the Internet Archive holds HTTP-200 snapshots of the
publisher's own OA PDF URLs. Both files were retrieved from there:

| file | source |
|---|---|
| article (7 pp.) | `web.archive.org/web/20240530035922id_/https://www.publish.csiro.au/CH/pdf/CH08240` |
| **Accessory Publication, Table S1** (13 pp.) | `web.archive.org/web/2015/https://www.publish.csiro.au/?act=view_file&file_id=CH08240_AC.pdf` |

Table S1 is the per-compound table. Its own header reads, verbatim:

> *"Table S1. Structures, physicochemical parameters, calculated Dmcl (from Eq. 1
> using the coefficients for HT29 in Study A as listed in Table 1) and measured
> Dmcl for HT29 MCLs. Errors are SEM and units of Dmcl are cm²s⁻¹."*

The paper states Table S1 was available from the journal *"until September
2012"*, so the Archive copy is now the accessible route to it.

### 2.2 Extraction and its validation

The Table S1 PDF interleaves chemical-structure drawings with the data columns
and wraps several cells across lines, so naive text-flow parsing mis-assigns
values. Extraction was done instead by **pdfplumber word-coordinate column
assignment**, with the column boundaries re-detected on every page (they drift
from x≈356 to x≈397 across the 13 pages), and cross-checked against an
independent `pdftotext -layout` rendering.

**All 71 rows were recovered.** Every `d_meas` in
`experiments/diffusivity_falsification.py` is read from Table S1. None is
derived, inferred, imputed or predicted.

There is also an **independent check that the extraction is correct**. Table S1
carries the authors' own Eqn-1 prediction alongside the measurement. Regressing
the extracted predicted column against the extracted measured column gives
**R² = 0.843** (Spearman +0.897), against the **R² = 0.907** the paper reports
for its own fit. Those agree to the extent expected — the paper's R² is the
regression's own goodness-of-fit, while this is a recomputation on the two
published columns. Had the columns been mis-paired by the parser, this number
would have collapsed. It did not.

---

## 3. Result 1 — the spread. This is the decisive statistic.

| | min | max | **fold range** |
|---|---:|---:|---:|
| **measured `D_mcl`** | 12.8 | 380.0 µm²/s | **29.7×** |
| **engine `D_eff`** | 183.2 | 430.8 µm²/s | **2.4×** |
| the paper's own Eqn 1 | 15.0 | 339.0 µm²/s | 22.6× |

**The engine reproduces 7.9 % of the measured spread.** The 2.4× it does produce
is entirely the MW term, and most of it comes from urea (60 Da) — across the
drug-like part of the set the engine is far flatter still.

The authors say the same thing in their own text, independently of this
analysis:

> *"MW changes in the range typical of small molecule drugs (200–500 daltons)
> have little effect."* — Pruijn 2008, Sensitivity Analysis

## 4. Result 2 — correlation with measurement

| | Pearson r (log₁₀) | R² | Spearman |
|---|---:|---:|---:|
| **engine vs measured** | **+0.168** | **0.028** | +0.099 |
| paper's Eqn 1 vs measured | +0.918 | 0.843 | +0.897 |

The engine has **no predictive relationship** with the measurement. This is not
an accuracy problem that a scale factor could fix; there is essentially nothing
to rescale.

## 5. Result 3 — the molecular-weight exponent

| | value |
|---|---:|
| d log₁₀(`D_mcl`) / d log₁₀(MW), measured, marginal | −0.383 |
| d log₁₀(`D_eff`) / d log₁₀(MW), engine | −0.390 |
| Pruijn 2008 Table 1 coefficient **b** (multivariate fit, n = 69) | **−0.62 ± 0.14** |
| variance in measured log₁₀(`D_mcl`) explained by MW alone | **0.028** |

Read this carefully, because it is the one place the engine looks good and it is
misleading. The engine's marginal MW slope (−0.390) is almost exactly the
measured marginal slope (−0.383). **That agreement is worthless**, because the
MW relationship explains 2.8 % of the variance. The engine has approximately the
right slope on approximately the wrong axis. The paper's own multivariate
coefficient is −0.62 ± 0.14 — steeper, and estimated *after* controlling for
lipophilicity and hydrogen bonding, which is why it differs from the marginal
slope. These two exponents are not the same quantity and are not pooled here.

## 6. Result 4 — the kill shot: spread at matched molecular weight

A MW-only model cannot produce *any* spread among compounds of equal MW.

| MW band | n | MW spread | measured range (µm²/s) | **measured fold** | **engine fold** |
|---|---:|---|---|---:|---:|
| 160–200 | 10 | 162–195 | 21.1 – 380.0 | **18.0×** | 1.07× |
| 200–225 | 7 | 207–221 | 35.4 – 209.0 | **5.9×** | 1.03× |
| 240–265 | 12 | 246–263 | 13.9 – 270.0 | **19.4×** | 1.02× |
| 285–310 | 12 | 288–305 | 36.3 – 171.0 | **4.7×** | 1.02× |
| 310–350 | 11 | 312–348 | 60.6 – 240.0 | **4.0×** | 1.05× |
| 450–500 | 8 | 455–498 | 12.8 – 147.0 | **11.5×** | 1.04× |

And at **exactly identical** molecular weight, where the engine is flat by
construction:

| MW | n | measured range (µm²/s) | measured fold | engine fold | compounds |
|---:|---:|---|---:|---:|---|
| 192 | 3 | 54.3 – 75.0 | 1.4× | 1.00× | SN27758, SN27759, SN27769 |
| 213 | 3 | 35.4 – 87.0 | 2.5× | 1.00× | SN27791, SN27800, SN27891 |
| **263** | **4** | **13.9 – 50.9** | **3.7×** | **1.00×** | SN28402, SN28427, SN28433, SN28682 |
| 289 | 3 | 59.0 – 83.5 | 1.4× | 1.00× | SN28620, SN29143, SN29147 |
| **305** | **4** | **45.2 – 171.0** | **3.8×** | **1.00×** | SN28643, SN29102, SN29267, SN29662 |

This is the cleanest possible form of the test: same method, same cell line,
same laboratory, same molecular weight, and the measurement still moves nearly
4-fold. No parameter choice inside the current model structure can reproduce
that, because MW is the only thing the model is allowed to look at.

## 7. Result 5 — what the spread actually tracks

| predictor of measured log₁₀(`D_mcl`) | n | Pearson r | R² |
|---|---:|---:|---:|
| **H-bond donors** | 71 | **−0.553** | **0.305** |
| **logP₇.₄** | 68 | **+0.538** | **0.290** |
| log₁₀(MW) — *the engine's only input* | 71 | −0.166 | 0.028 |

Both dominant predictors are **already computed** by
`descriptors.compute_descriptors` and are simply never passed to the diffusivity
calculation. The gap is structural, not informational.

The mechanism is stated by the authors: hydrophilic compounds are restricted to
the **paracellular** route (between cells), while lipophilic compounds also
access the **transcellular** route (through cell membranes). The engine's only
obstacle population is "ECM" with a Brownian tortuosity and an Ogston–Renkin
fibre-mesh gate — it has **no membrane, and therefore no transcellular route at
all**. Whether a molecule can cross a lipid bilayer is invisible to it.

This is the same structural gap `results/external_anchor_survey.md` §3.2 found
from a completely different direction — the 12-fold measured MCL hindrance
against the engine's 1.77-fold, diagnosed there as "cellular volume exclusion is
not in the hindrance law". Two independent measurements, one cause: **cells are
not in the transport model.**

## 8. Result 6 — absolute accuracy (secondary)

| | value |
|---|---:|
| engine/measured ratio, min → median → max | 0.76× → **3.68×** → 17.84× |
| within 3× of measurement | **27 / 71 = 38 %** |
| geometric-mean bias | **3.49× too high** |

The engine over-predicts tissue diffusivity by ~3.5× on average. This is
consistent in sign and magnitude with §3.1 of `external_anchor_survey.md`
(287.3 vs 74–130 µm²/s for tirapazamine, 2.2–3.9× over) — and note the present
dataset contains tirapazamine itself, measured at **40.4 µm²/s** in HT29,
against the engine's 287.3, a **7.1×** over-prediction in a third cell line.

The bias is the lesser finding. A systematic 3.5× could in principle be absorbed
by a calibration constant. **The 29.7× spread the model cannot represent at all
could not be.**

---

## 9. What this does and does not license saying

**Established, and defensible:**
* The engine's tissue diffusivity is a function of molecular weight and ECM
  fraction only. Verified in code and empirically, 0.01 % agreement.
* Against 71 same-method same-cell-line measurements it captures 7.9 % of the
  observed spread at R² = 0.028, and is flat where the measurement moves 3.7×
  at fixed MW.
* The dominant physicochemical determinants of tissue diffusivity — hydrogen
  bonding and lipophilicity — are computed by the codebase and discarded before
  the transport calculation.

**Not established, and must not be claimed:**
* That fixing this fixes the engine's predicted concentration *fields*. This
  test scores one PDE coefficient, not a solved field.
* That the Pruijn Eqn-1 QSAR should simply be dropped into `compile_drug.py`.
  It is fitted to HT29 MCLs — packed cells with little stroma — and the paper
  itself shows the asymptote coefficients change materially across SiHa,
  HCT8-Ea and HCT8-Ra (up to 100-fold at high logP₇.₄ for HCT8-Ea). Adopting it
  wholesale would import a cell-line-specific fit as if it were physics, and
  would make the engine's diffusivity in-sample against this dataset, destroying
  the one clean out-of-sample transport test the project now has.
* Anything about clinical response. Out of scope, as in the parent survey.

**The honest one-line statement.** *Vivum's tissue diffusivity model is a
molecular-weight scaling. Measured against 71 tissue diffusion coefficients from
one method and one cell line, it reproduces 8 % of the observed 30-fold spread
and is uncorrelated with the measurement (R² = 0.03), because the two properties
that actually govern extravascular transport — hydrogen-bond donor count and
lipophilicity — are computed by the code and then discarded. The model is
falsified as a predictor of compound-to-compound differences in tissue
diffusivity.*

---

## 10. Reproducing this

```
./.venv/bin/python experiments/diffusivity_falsification.py
```

Runs in ~45 s. Compiles all 71 molecular weights through the unmodified
`compile_drug.build_profile` path at the repo-default ECM reference of 0.4,
prints the pre-flight MW-only proof, the full per-compound table, and Results
1–6 above. The 71 measured values are embedded in the script with the
provenance of every one recorded in its module docstring.
