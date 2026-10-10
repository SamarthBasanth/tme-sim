# Scoring the ECM/diffusivity coupling against external measurement — 2026-10-10

**Status: PRE-REGISTERED SCORING PASS. No parameter was tuned. Nothing under
`tmesim/**` was edited.** The pre-registration was written and frozen before any
engine number was computed, and is preserved verbatim at
`experiments/ecm_coupling_score/PREREGISTRATION.md`.

**Calibration guard, start of pass:**

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

*(Closing guard in §9.)*

**Commission.** Score the ECM/diffusivity coupling
`D_eff = D·(1 − ALPHA_TORTUOSITY·ecm)` against the two external anchors that
`results/transport_falsifier_scouting_2026-10-09.md` ranked #1 and #4, as a
pre-registered test: Shafieyan 2014 (a FORM-level refutation candidate) and
Shi 2014 (a BAND-level refutation candidate).

---

## 0. HEADLINE — what this pass found

**Both candidate refutations fail, and they fail for the same underlying reason:
each compared a dimensionless engine quantity against a measured quantity with a
different denominator.** Neither is a refutation of the shipped law. Both are
corrections to `transport_falsifier_scouting_2026-10-09.md`.

1. **Candidate 1 (Shafieyan, FORM) is a NON-RESULT, by the pre-registered rule.**
   The implied coefficient scales as exactly `1/A` in the unit-mapping scale `A`,
   and `A` is a free property of the mapping, not of the engine. Three defensible
   mappings give three *qualitatively different* verdicts — inadmissible,
   admissible, and saturated. The `a > 1` inadmissibility is **an artefact of one
   mapping choice**, and that mapping is the one *least* faithful to the code.
2. **The one mapping-robust statement:** run the shipped `α = 0.6` forward, and
   the predicted-vs-measured endpoint ratio is off by **1.22× to 1.58×** in all
   nine solute × mapping cells — but the **SIGN flips with the mapping** (too
   flat under M1/M3, too steep under M2). The measurement constrains the coupling
   to within ~1.6×, and **cannot say in which direction**.
3. **Candidate 2 (Shi, BAND) inverts in sign once the denominator is fixed.** The
   scouting record compared Shi's `D_eff/D_free` against the *spatial factor
   alone*. Shi's denominator is free water; the engine's band denominator is the
   **already-hindered compiled** diffusivity, after a 12.58× pre-spatial chain.
   Against the **full chain** the engine's `D_eff/D_free` is **[0.0318, 0.0795]**,
   which is **below all six** measured ratios. The engine is **6.06× too SLOW**,
   not "~2.1× too fast".
4. **Therefore the scouting record's error partition is void.** It derived
   "~2.1× in `D_eff`, ~2.9× in `k`". With the sign corrected, a `D_eff` that is
   6.06× too *small* against a `D_eff/k` that is ~6× too *large*
   (`kyle2004_penetration.md`) requires `k` too small by **~36×**, not 2.9×.
   The partition claim should not be used.
5. **A new structural defect, found by a forward-substitution check on my own
   arithmetic:** `_effective_diffusivity` clips `ecm` to `[0,1]`. Any mapping that
   puts a measured point above `ecm = 1` gives the engine **zero dynamic range**
   over the measured interval — it predicts a diffusivity ratio of exactly
   **1.0000** against measured 1.41× and 1.57×. Under the engine's own declared
   matrix volume fraction (`φ_f = 0.05`) that is what happens for 2 of 3 solutes.

---

## 1. The three structural facts, re-verified by running the shipped code

Not read off a docstring — driven through `_effective_diffusivity` itself
(`experiments/ecm_coupling_score/score_ecm_coupling.py` §0).

| fact | `file:line` | verbatim |
|---|---|---|
| ECM field is MAX-NORMALISED per section | `tmesim/reconstruct/real_fields.py:112-114` | `dmax = float(density.max())` / `if dmax > 0.0:` / `density /= dmax` |
| the coupling | `tmesim/physics/diffusion.py:88` | `d = drug.diffusivity * (1.0 - ALPHA_TORTUOSITY * ecm)` |
| the constant | `tmesim/physics/diffusion.py:75` | `ALPHA_TORTUOSITY = 0.6` |
| the floor | `tmesim/physics/diffusion.py:89` | `floor = max(1e-6, 1e-3 * drug.diffusivity)` |
| the clip | `tmesim/physics/diffusion.py:87` | `ecm = np.clip(np.asarray(ecm_density, ...), 0.0, 1.0)` |

Measured output of the shipped function:

```
 ecm   D_eff/D
 0.00    1.0000
 0.50    0.7000
 1.00    0.4000
  -> reachable band D_eff/D = [0.4000, 1.0000]
  -> ecm = [-5, 0, 1, 5] clips to D_eff/D = [1.0, 1.0, 0.4, 0.4]
```

**All three facts CONFIRMED.** `D_eff/D` is confined to `[0.40, 1.00]` by
construction, independent of any tissue, and the band is *hard* — the internal
`np.clip` means even an out-of-range field cannot escape it.

**One addition the scouting record did not state:** the engine also hard-wires the
matrix itself. `ecm_mesh.ogston_renkin_hindrance` defaults to
`fiber_radius_nm = 2.0, fiber_volume_fraction = 0.05`, and `compile_drug.py:174`
calls it as `ogston_renkin_hindrance(radius_nm)` — radius only. So **φ_f = 0.05
for every compound and every tissue**, and that number matters in §3 below.

---

## 2. Anchor verification — what I confirmed and what I CORRECTED

Every anchor was re-fetched from a primary source by me, not taken from the
scouting record.

### 2.1 Shafieyan 2014 — CONFIRMED verbatim, with three corrections

> Shafieyan Y, Khosravi N, Moeini M, Quinn TM. "Diffusion of MRI and CT Contrast
> Agents in Articular Cartilage under Static Compression." *Biophys J*
> 107(2):485–492. **PMID 25028890, PMCID PMC4104049**, DOI 10.1016/j.bpj.2014.04.041.

Authors, journal, date, volume, issue and pages **re-confirmed independently by
PubMed `efetch` and `esummary`**: `sortfirstauthor = Shafieyan Y`, author list
`Shafieyan, Khosravi, Moeini, Quinn`, `source = Biophys J`, `pubdate = 2014 Jul
15`, `volume = 107`, `issue = 2`, `pages = 485-492`. Every field matches the
scouting record. Its warning that this paper circulates under a *wrong*
attribution ("Pouran et al., 107(9)") is therefore **upheld**.

**The key sentence, re-read by me from the PMC full text, character for character:**

> "Increased diffusivity significantly correlated with increasing fluid volume
> fraction for all three contrast agents. Diffusivity increased from 400 ± 123 to
> 778 ± 111, 110 ± 19 to 155 ± 14, and 124 ± 28 to 195 ± 26 μm²/s when fluid
> volume fraction increased from 0.66 to 0.97, 0.78 to 0.94, and 0.80 to 0.95 for
> sodium iodide, sodium diatrizoate, and Gd-DTPA, respectively (Fig. 6)."

**CONFIRMED exactly as the scouting record quoted it.** It is **body text, not a
figure read by eye** — no digitisation was performed by anyone.

**CORRECTION 1 — the paper reports NO tortuosity coefficient at all.** The string
`tortuos` occurs **zero times** in the PMC full text (counted
programmatically). The task framing given to me stated that "the implied
tortuosity coefficients are **reported as** 1.50 / 1.64 / 2.16". **That is wrong.**
Those three numbers are **derived by inverting the engine's own linear law**
against the paper's endpoints. They are a property of `tmesim`, not a measurement.
Any future text calling them *reported* or *measured* tortuosity coefficients
must be corrected. (Pre-registered as **C1c**; this is a FAIL for the framing.)

**CORRECTION 2 — the `p = 0.061` caveat is on the wrong axis.** The scouting
record attaches "sodium diatrizoate's trend was **not** significant (*p* = 0.061)"
to its fluid-volume-fraction table. The full text puts that *p* on the
**compression** axis (their Fig. 4B):

> "Although sodium diatrizoate approached a trend of decreasing diffusivity from
> 145 ± 19 to 116 ± 19 μm²/s when compression changed from 0% to 52% (Fig. 4 B),
> this was not statistically significant (p = 0.061)."

On the **fluid-fraction** axis — the axis this test actually uses — the paper says:

> "Increasing trends were found to be significant as confirmed by analysis of
> variance and p value was equal or<0.005 for all contrast agents."

So all three fluid-fraction slopes **are** significant, *including* sodium
diatrizoate. This correction is in the **anchor's favour**.

**CORRECTION 3 — a charge confound the scouting record did not record.** The
paper's own partition coefficients are **sodium iodide 0.740 ± 0.0415, Gd-DTPA
0.469 ± 0.0142, sodium diatrizoate 0.385 ± 0.229**, and it attributes them to
electrostatics, not geometry:

> "As these molecules carry electric charges (iodide: −1, diatrizoate: −1, and
> Gd-DTPA: −2) identical to the negatively charged proteoglycans of the matrix
> solid network, partition coefficients of less than unity were expected…"

The steepest implied coefficient (Gd-DTPA, 2.16) belongs to the only **−2**
solute. With *n = 3* that is **suggestive, not a test**, and I am not treating it
as one — but a purely geometric `(1 − α·ecm)` has **no charge term at all**, so
part of this measured slope is a physics the engine does not model. This weakens
the anchor as a test of a *tortuosity* law independently of any mapping.

**The fluid-fraction definition, needed for §3** (their Eq. 5, verbatim):

> "The total explant fluid volume fraction under compression was calculated as
> the ratio between the explant fluid volume under compression (V_f) and the sum
> of the fluid volume and the solid volume (V_s) under compression (total
> volume): ϕ = V_f / (V_f + V_s)."

Gravimetric, from wet/dry weights after lyophilisation, assuming ρ_f = 1 and
ρ_s = 1.4 g/mL. **This is a true volume fraction.**

### 2.2 Shi 2014 — CONFIRMED verbatim, with one counting correction

> Shi L, Zeng M, Sun Y, Fu BM. "Quantification of blood-brain barrier solute
> permeability and brain transport by multiphoton microscopy." *J Biomech Eng*
> 2014;**136(3):031005**. **PMID 24193698**. Record and full abstract retrieved
> by me via PubMed `efetch`.

> "The mean Deff values were 33.2, 4.4, 1.3, 0.89, 0.59, and 0.47 × 10-7 cm2/s,
> respectively, for sodium fluorescein, dextran-4k, -20k, -40k, -70k, and IgG.
> The corresponding mean ratio of Deff to the free diffusion coefficient Dfree,
> Deff/Dfree, were 0.46, 0.19, 0.12, 0.12, 0.11, and 0.11 for these solutes."

**All six ratios CONFIRMED verbatim.** Tissue: "rat brain tissue 100-250 μm below
the pia mater". Distance axis CONFIRMED: "P and Deff were determined from the
rate of tissue solute accumulation and the radial concentration gradient around
individual microvessels in the brain tissue."

**CORRECTION 4 — the count is 5, not 4.** The scouting record says "**four** of
its six measured `D_eff/D_free` ratios (0.19 → 0.11) lie outside the shipped
reachable band". There are **five** non-small solutes (dextran-4k, -20k, -40k,
-70k, IgG) and five values outside `[0.40, 1.00]`: 0.19, 0.12, 0.12, 0.11, 0.11.
Only sodium fluorescein's 0.46 is inside. Pre-registered criterion **C2a is
therefore recorded as FAIL** — the claim as published is off by one. The error is
in the *conservative* direction, so the band claim is slightly **stronger** than
published, but a wrong count is a wrong count.

**One half-quote to flag.** The scouting record quotes only "there is no
significant difference in Deff/Dfree between solutes with molecular weights from
20,000 to 160,000 Da". The full sentence begins: "**While there is a significant
difference in Deff/Dfree between small (e.g., sodium fluorescein) and larger
solutes,** there is no significant difference…". The size-dependence *is*
significant across the small/large boundary — which is the half that matters for
a law with no size term in it.

### 2.3 Counterweights — CONFIRMED verbatim (pre-registered C3)

Both re-fetched by me via PubMed `efetch`.

* **Nugent & Jain 1984** (PMID 6197161, *Cancer Res* 44(1):238–44), verbatim:
  "**Macromolecular transport in tumor tissue was hindered to a lesser extent
  than in normal tissue**, which is consistent with reports of reduced contents
  of glycosaminoglycans, and markedly large interstitial space in tumors."
* **Chary & Jain 1989** (PMID 2748592, *PNAS* 86(14):5385–9), verbatim: "The
  average fluid velocity was about 0.6 microns/s, and albumin diffusion
  coefficients were **5.8 +/- 1.3 x 10(-7) cm2/s and 6.3 +/- 1.9 x 10(-7) cm2/s
  in normal and neoplastic tissues**, respectively. The interstitial fluid flow,
  in general, was directed into postcapillary venules."

**Both report tumour hindered LESS than, or the same as, normal tissue — opposite
in sign to Shi/Fu.** Chary & Jain's two values overlap within one SD: albumin is
**not** differentially hindered in tumour at all. The engine raises tortuosity
where matrix markers are dense, which is the Shi/Fu direction.

**C3 PASS**, with the pre-committed statement made explicitly: **the literature
does not agree on the sign, and `ALPHA_TORTUOSITY` cannot be validated by picking
one study.** A pass quoting only Shi would manufacture a verdict. The scouting
record's own impeachment also stands: Kutuzov 2018 criticises Shi's
radial-symmetry assumption, so Candidates 2 and 9 there partly impeach each other
and must travel together.

---

## 3. The unit mapping — the crux, declared before the run

`ecm_density` and `φ_fluid` are **not the same quantity**:

| | Shafieyan `φ = V_f/(V_f+V_s)` | engine `ecm_density` |
|---|---|---|
| nature | true volume fraction | dimensionless **relative** index |
| scale | absolute, gravimetric | **max-normalised per section** (`density /= dmax`) |
| what `1.0` means | zero solid | **the densest pixel of whatever section this is** |
| units | dimensionless, physical | dimensionless, **unanchored** |

Writing `x = 1 − φ_fluid` (solid volume fraction) and `ecm = A·x + B`:
`B = 0` is **forced**, because the shipped law gives `D_eff = D` exactly at
`ecm = 0`. Only the scale `A` is free — and **`A` is a property of the mapping,
not of the engine.** Three mappings were declared in advance:

| mapping | `ecm = 1` denotes | `A` | what it assumes |
|---|---|---:|---|
| **M1** | pure solid, `φ_solid = 1` | 1 | that a section's densest pixel is water-free solid |
| **M2** | densest explant in that solute's own series | `1/x_max` | that the densest cartilage explant and the densest GBM pixel are the same object |
| **M3** | the engine's own `φ_f = 0.05` | 20 | that `φ_f = 0.05` is a real matrix volume fraction and `ecm = 1` saturates there |

**M2 is the mapping faithful to the code**, because `real_fields.py:114`
max-normalises. **M1 is the scouting record's mapping**, and it is the one
contradicted by the engine's own construction. `results/heldout_patient.md` §3.5
refutes M1 empirically: mean ECM **0.076 vs 0.457 on the same tumour type**,
driven by which matrix genes were on the panel — "panel content, not tissue
biology". A field whose mean moves 6× with the gene panel is not a solid volume
fraction.

**This was already repo doctrine, and the scouting record quoted it and then
violated it.** `results/reconstruction_validation_roadmap.md:98`: "**No published
measurement can validate an amplitude that is defined only relative to the
densest point in one section.**"

---

## 4. CRITERION 1 — result

### C1a (confirmatory): PASS

The scouting record's M1 arithmetic is correct.

| solute | `φ_fluid` | `x = φ_solid` | measured D ratio | implied `α` (M1) | published |
|---|---|---|---:|---:|---:|
| sodium iodide | 0.66→0.97 | 0.34→0.03 | 1.9450 | **1.4969** | 1.50 |
| sodium diatrizoate | 0.78→0.94 | 0.22→0.06 | 1.4091 | **1.6364** | 1.64 |
| Gd-DTPA | 0.80→0.95 | 0.20→0.05 | 1.5726 | **2.1646** | 2.16 |

All three match to ±0.02. Its "too flat by" column (1.58 / 1.27 / 1.43) also
reproduces exactly (§5). **This confirms arithmetic, not physics.**

### C1b (the real, pre-registered test): MAPPING-DEPENDENT → NON-RESULT

| solute | M1 (`A`=1) | M2 (`A`=1/x_max) | M3 (`A`=20) |
|---|---:|---:|---:|
| sodium iodide | 1.4969 | 0.5089 | 0.0748 |
| sodium diatrizoate | 1.6364 | 0.3600 | 0.0818 |
| Gd-DTPA | 2.1646 | 0.4329 | 0.1082 |
| **verdict** | **INADMISSIBLE (α>1)** | **ADMISSIBLE (α<1)** | **OUT OF DOMAIN** |

The implied `α` scales as **exactly `1/A`** (verified over `A` ∈ [1, 20]). The
three verdicts are qualitatively different, so by the **pre-declared decision
rule** this outcome is recorded as a **NON-RESULT for the form claim** — not as a
pass for either side.

**The flip point is the whole test.** `α = 1` at `x_ref = 1/α_M1`:

| solute | solid fraction that `ecm = 1` must denote for the slope to be inadmissible |
|---|---:|
| sodium iodide | **> 0.668** |
| sodium diatrizoate | **> 0.611** |
| Gd-DTPA | **> 0.462** |

So Candidate 1 refutes the form **if and only if** the engine's `ecm = 1` pixel
denotes a solid volume fraction above ~0.46–0.67. **It does not.** The engine's
own hard-wired matrix volume fraction is **0.05** — an order of magnitude below
the flip point. Under the code's own number the slope is nowhere near
inadmissible.

### The out-of-domain finding (not pre-registered; found by self-checking)

A forward-substitution check (`experiments/ecm_coupling_score/verify_inverse.py`)
was written to catch errors in my own inverse algebra by pushing each implied `α`
back through the **shipped** `_effective_diffusivity`. It caught one — in my own
M3 analysis, not in the engine:

```
sodium iodide      M3  a=0.0748  ecm=(0.6000,6.8000)  fwd ratio=1.0324 vs measured 1.9450  OUT OF DOMAIN
sodium diatrizoate M3  a=0.0818  ecm=(1.2000,4.4000)  fwd ratio=1.0000 vs measured 1.4091  OUT OF DOMAIN
Gd-DTPA            M3  a=0.1082  ecm=(1.0000,4.0000)  fwd ratio=1.0000 vs measured 1.5726  OUT OF DOMAIN
```

`ecm` is `np.clip`'d to `[0,1]` **inside** the shipped function. Under M3 both
endpoints clip to `ecm = 1` for two of three solutes, and the engine predicts a
diffusivity ratio of **exactly 1.0000** against measured 1.41× and 1.57× — **zero
dynamic range**. M3 is therefore reported as **SATURATED, not "admissible with a
small α"**; the implied-`α` solve is invalid out of domain. M1 and M2 are
in domain and reproduce the measured ratio exactly through the shipped code.

**This is a real structural limitation worth recording on its own:** if a tissue's
matrix content ever genuinely saturates the normalised field, the coupling loses
all spatial gradient, silently.

### What the shipped floor actually does (correcting the "goes NEGATIVE" claim)

The framing given to me, and the scouting record, say that `a > 1` "sends `D_eff`
**negative** above `ecm = 1/a`". **Driven through the shipped code, that is not
what happens.** The bare expression does go negative; the shipped function
**floors** at `max(1e-6, 1e-3·D)`:

```
alpha = 1.4969,  D = 100,  floor = 0.1
    ecm   bare D(1-a*ecm)   shipped D_eff   shipped/D
 0.0000           100.0000        100.0000     1.00000
 0.4000            40.1240         40.1240     0.40124
 0.6680             0.0000          0.1000     0.00100
 0.8000           -19.7520          0.1000     0.00100
 1.0000           -49.6900          0.1000     0.00100
  -> shipped min = 0.1 (FLOORED, never negative)
```

So the solver stays well-posed; it does not break. What `a > 1` actually produces
is a **1000× cliff followed by a dead flat region** — above `ecm = 1/a` the law
carries **no gradient at all**. That is still not a tortuosity law, but it is
inadmissible for a *different reason* than negativity, and the published reason is
wrong.

### The one mapping-robust statement — the headline

Drop the inverse solve entirely. Run the **shipped `α = 0.6`** forward and compare
predicted against measured. Nothing is fitted:

| solute | mapping | `ecm_lo` | `ecm_hi` | pred ratio | measured | error | sign |
|---|---|---:|---:|---:|---:|---:|---|
| sodium iodide | M1 | 0.0300 | 0.3400 | 1.2337 | 1.9450 | 1.577× | **FLAT** |
| sodium iodide | M2 | 0.0882 | 1.0000 | 2.3676 | 1.9450 | 1.217× | **STEEP** |
| sodium iodide | M3 | 0.6000 | 1.0000 | 1.6000 | 1.9450 | 1.216× | **FLAT** |
| sodium diatrizoate | M1 | 0.0600 | 0.2200 | 1.1106 | 1.4091 | 1.269× | **FLAT** |
| sodium diatrizoate | M2 | 0.2727 | 1.0000 | 2.0909 | 1.4091 | 1.484× | **STEEP** |
| sodium diatrizoate | M3 | 1.0000 | 1.0000 | 1.0000 | 1.4091 | 1.409× | **FLAT** |
| Gd-DTPA | M1 | 0.0500 | 0.2000 | 1.1023 | 1.5726 | 1.427× | **FLAT** |
| Gd-DTPA | M2 | 0.2500 | 1.0000 | 2.1250 | 1.5726 | 1.351× | **STEEP** |
| Gd-DTPA | M3 | 1.0000 | 1.0000 | 1.0000 | 1.5726 | 1.573× | **FLAT** |

* **Magnitude: 1.216×–1.577× in all nine cells.** Remarkably tight.
* **Sign: FLAT under M1 and M3, STEEP under M2.** Set entirely by the mapping.

**This is the defensible result of Candidate 1.** The measurement pins the
coupling to within about a **factor of 1.6** and is **silent on direction**. That
is a real, if modest, external constraint on `ALPHA_TORTUOSITY = 0.6` — and it is
the *opposite* of a form refutation.

---

## 5. CRITERION 2 — result: the sign inverts

### C2a (confirmatory): FAIL — see Correction 4. Five of six, not four.

### C2b (the real, pre-registered test): the denominator was wrong

Shi's `D_eff/D_free` has **free water** in the denominator. The engine's
`[0.40, 1.00]` band has `drug.diffusivity` — the **already-hindered compiled**
value — in the denominator. Comparing them directly is a units error. The
comparable engine quantity is the **full chain**, read from the engine's own pure
functions (read-only; `compile_drug.build_profile`, not a registered run):

| | doxorubicin | temozolomide |
|---|---:|---:|
| `free_D` (Stokes–Einstein, T = 310 K) | 351.31 µm²/s | 495.08 µm²/s |
| compiled `effective_D` | 27.93 µm²/s | 42.42 µm²/s |
| `compiled_D / free_D` | **0.07949** (12.58× hindrance) | **0.08569** (11.67×) |
| **FULL-CHAIN `D_eff/D_free` band** | **[0.03180, 0.07949]** | [0.03428, 0.08569] |
| at flagship mean `ecm = 0.076` | **0.07586** | 0.08179 |

(The scouting record's §0.2 table said 351.31 / 27.93 / 12.7× — confirmed; the
hindrance is 12.58×, so "12.7×" is a rounding of a slightly different quantity.)

**C2b-i: NONE of the six measured ratios falls inside the full-chain band.** All
six are **above** it. Even the most hindered measurement (0.11, IgG at ~160 kDa)
is **1.38× faster** than the engine's maximum (0.0795).

**C2b-ii: engine 6.06× too SLOW, not 2.1× too fast.** Shi's small-solute 0.46
against the engine's full chain at the flagship `ecm = 0.076` (0.07586) gives
**6.06×**, with the engine on the **hindered** side. The scouting record's
"~2.1× too fast" compared 0.46 against the *spatial factor alone* (0.954) and so
measured the wrong ratio. **Pre-registered sign rule fires: the published sign is
wrong.**

### Consequence — the error partition is void

`results/kyle2004_penetration.md` reports `λ` 2.45× too long, i.e. `D_eff/k` ~6×
too large. The scouting record split that as "~2.1× in `D_eff`, ~2.9× in `k`".
With the sign corrected, a `D_eff` that is **6.06× too small** against a
`D_eff/k` that is ~6× too **large** requires `k` too small by **~36×** — not
2.9×. **The published partition should not be used.** Note this is not a new
positive attribution either: it is a different *arithmetic consequence* of a
cross-species anchor, and §6 explains why I am not promoting it.

### C2c: the caveats that must travel with any use of this

Shi's tissue is **normal rat brain**, not human GBM. His small solute is **sodium
fluorescein (376 Da)**, not a repo compound. His denominator is free water at an
unstated temperature while the engine's `free_D` is Stokes–Einstein at 310 K.
Kutuzov 2018 impeaches his radial-symmetry assumption. This is a
cross-species, cross-tissue, cross-compound comparison and **cannot carry a
quantitative correction on its own.**

---

## 6. Honest verdict

**On the FORM: the linear form SURVIVES. Both refutation candidates fail.**

* Candidate 1 does not refute the form. Its `a > 1` depends entirely on a mapping
  (`ecm ≡ φ_solid`) that contradicts the engine's own max-normalisation and is
  refuted empirically by `heldout_patient.md` §3.5. Under the code-faithful
  mapping the implied coefficients are **0.36–0.51**, all admissible, all within
  1.7× of the shipped 0.6 — and under the engine's own `φ_f = 0.05` the field
  saturates instead. **Three mappings, three verdicts.**
* Candidate 2 does not refute the form either. Its four-ratios-outside-the-band
  argument (really five) compares against the wrong denominator; against the full
  chain the engine is slower than *every* measured value, which is an amplitude
  failure, not a reachability failure.

**So, as commissioned: can the comparison be made defensibly?**

**At the FORM level — NO.** The quantity the engine's `ecm = 1` denotes is
undefined, the verdict scales as `1/A` in an unidentifiable mapping constant, and
the three defensible mappings disagree qualitatively. **The comparison cannot be
made defensibly, and that is the finding.** This is not a new discovery so much as
a confirmation of `reconstruction_validation_roadmap.md:98`, which already
forbade exactly this comparison — and which the scouting record quoted, in the
same file, before performing it anyway. **That is the claim-surface-drift pattern
reproducing itself inside a single document.**

**At the AMPLITUDE level — PARTIALLY, and only as a bound.** Two mapping-robust
numbers survive:

1. The shipped `α = 0.6` mispredicts Shafieyan's endpoint ratio by **1.22×–1.58×**
   under every mapping, **direction unknown**. `ALPHA_TORTUOSITY = 0.6` is
   therefore externally constrained to roughly a factor of 1.6, which is **better
   than "constrained by nothing"** (`verification_routes_2026-09-20.md` R5) and
   **consistent with** `uncertainty.md`'s Sykova & Nicholson sanity check.
2. The engine's full-chain `D_eff/D_free` **[0.0318, 0.0795]** lies **below all
   six** of Shi's measured ratios. Direction here is unambiguous: the engine's
   total hindrance chain is **too aggressive**. But this scores the *whole chain*
   — dominated by the two-route cell factor 0.1512, not by the spatial coupling —
   so it is **not** an indictment of `ALPHA_TORTUOSITY` specifically.

**What I will not claim.** The literature sign is unresolved (§2.3: Chary & Jain
and Nugent & Jain both point the other way, and Chary & Jain find *no* tumour
difference at all). A charge confound sits inside Candidate 1 that the engine does
not model at all (§2.1, Correction 3). One anchor impeaches the other. **No new
value of `ALPHA_TORTUOSITY` was computed, suggested, or installed, and none
should be on this evidence.**

**Net effect of this pass: two published refutation claims withdrawn, four
anchor-number corrections, one new structural limitation recorded, and one modest
but genuine external bound (~1.6×, direction unknown) established on a constant
previously described as uncited.**

---

## 7. Corrections this pass makes to `transport_falsifier_scouting_2026-10-09.md`

That file is a dated published record and **was not edited.** Corrections are
recorded here only.

| # | its claim | correction |
|---|---|---|
| 1 | §2.1 implied coefficients are "tortuosity coefficients" | The paper never uses the word "tortuosity" (0 occurrences). They are derived by inverting `tmesim`'s own law. |
| 2 | §2.1 "sodium diatrizoate's trend was not significant (p = 0.061)", in the fluid-fraction table | That *p* is on the **compression** axis. On the fluid-fraction axis the paper reports *p* ≤ 0.005 for **all three**. |
| 3 | §2.1 "`D_eff` goes **negative** above `ecm = 1/a`" | The shipped function **floors** at `1e-3·D`; it never goes negative. The real effect is a 1000× cliff and loss of all gradient. |
| 4 | §2.1 "the measured slope … breaks the linear form" | Mapping-dependent. Admissible under the code-faithful max-normalised mapping (α = 0.36–0.51), saturated under the engine's own φ_f = 0.05. |
| 5 | §2.4 "**four** of its six ratios lie outside" | **Five** of six (0.19, 0.12, 0.12, 0.11, 0.11). Off by one. |
| 6 | §2.4 "the engine is **~2.1× too fast** in diffusivity" | **Sign inverted.** Against the full chain the engine is **6.06× too SLOW**. The 2.1× compared Shi's free-water-denominated ratio against the spatial factor alone. |
| 7 | §2.4 "partitions that 6× — ~2.1× in `D_eff`, ~2.9× in `k`" | **Void.** With the sign corrected, `k` would have to be ~36× too small. Do not use the partition. |
| 8 | §2.4 quotes only "no significant difference … 20,000 to 160,000 Da" | The sentence begins "**While there is a significant difference** … between small … and larger solutes". |

Credit where due: its provenance warning about the Pouran misattribution is
upheld, its M1 arithmetic is exactly right, its §0.1 analysis of the
max-normalisation problem is correct and is what defeats its own §2.1, and its
§2.11/§5.2 dissent (Chary & Jain, Nugent & Jain) is the reason this pass could
report both sides.

---

## 8. Artefacts

| path | what |
|---|---|
| `experiments/ecm_coupling_score/PREREGISTRATION.md` | frozen pre-registration, written before any engine number |
| `experiments/ecm_coupling_score/score_ecm_coupling.py` | the scoring run (§1, §4, §5) |
| `experiments/ecm_coupling_score/verify_inverse.py` | forward-substitution, floor and monotonicity checks (§4) |

Not registered in `tmesim/validation/scorecard.py`. No file under `tmesim/**`,
`app/**`, `docs/**`, `landing-3d/**` or `.github/**` was touched. No existing file
in `results/` or `archive/` was modified.

---

## 9. Calibration guard, end of pass

```
$ .venv/bin/python -c "from tmesim.validation.calibration import calibrate_alpha; a=calibrate_alpha(); assert a==22.223477711172563, a; print('alpha OK',a)"
alpha OK 22.223477711172563
```

`ALPHA_TORTUOSITY` on exit: **0.6**, unchanged, never written to.
