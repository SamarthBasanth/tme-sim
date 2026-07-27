# Calibration — one mechanistic constant for perivascular drug consumption

**Module:** `tmesim/validation/calibration.py` (run: `python -m tmesim.validation.calibration`)

**Scope:** fixes the calibration gap the PK reality check exposed (tmesim over-predicted doxorubicin tissue penetration ~15x: ~700 µm vs the published ~40–50 µm). We make the tissue loss rate *mechanistic* with exactly **one** free constant, fit it **once** on doxorubicin, then test whether it predicts a **second** drug with **no** refitting. Nothing in core, the engine, or the molecular modules is modified — this module only drives them and returns a normal `core.DrugParams`.

---

## The mechanism and the proxy

The frozen transport PDE loses drug at a first-order rate `k = clearance + binding`, giving penetration length `L = sqrt(D_eff / k)`. The structure-only priors give `k ≈ 3e-4 /s` (a ~1-hour timescale) → `L ≈ 700 µm`. Reproducing the 65 µm boundary layer (the e-fold equivalent of the published 40–50 µm half-distance — see the convention section below) needs `k = D_eff / L² ≈ 3e-2 /s` (a ~30-second timescale) — about **150× faster**. The missing physics is not slow metabolic clearance; it is **rapid cellular uptake + high-capacity intracellular / DNA binding**, a near-irreversible perivascular sink that is what actually carves the steep gradient.

We keep the engine's first-order sink but split it:

```
k_eff = k_passive + alpha * uptake_proxy(drug)

k_passive    = drug.clearance   (slow metabolic washout, ~1e-4 /s; NOT amplified)
uptake_proxy = drug.binding     (structure-based tissue/DNA BINDING-CAPACITY rate,
                                 binding.tissue_binding_rate: monotone in logP + HBD)
alpha        = the ONE free scaling constant (dimensionless)
```

**Why this proxy.** `uptake_proxy` is the molecular tissue-binding rate — a structure-defensible *ordering* of how strongly a drug is sequestered into cell lipids, protein and DNA (intercalators rank high, small hydrophilic drugs low). It is a well-behaved ~30× range (doxorubicin high, temozolomide low), unlike passive membrane permeability, which spans ~5 orders of magnitude in the logP model and mis-ranks doxorubicin (a low-passive-permeability P-gp substrate whose gradient is binding-, not permeability-, driven). `alpha` is the single factor that lifts the deliberately-conservative Stage-0 binding prior up to the true perivascular sink strength — it encodes the intracellular binding **capacity** (bound:free ratio) the structure-only heuristic cannot know.

Because `k_eff` is linear in `alpha`, calibration is one closed-form solve: `k_eff* = D_eff / L*²`, `alpha = (k_eff* − k_passive) / uptake_proxy`. No optimizer, no iteration.

---

## Unit convention: half-distance vs e-folding length (CORRECTED)

The published anchors are **half-distances** `d_half` — the distance from a vessel at which drug intensity falls to 50 % of its perivascular value. The simulator's `penetration_length_um` is an **e-folding length** `L`, defined by `C(r) ~ exp(−r/L)`. For an exponential profile these are different numbers:

```
d_half = L · ln2        =>        L = d_half / ln2 = 1.4427 · d_half
```

So doxorubicin's published **40–50 µm half-distance** is an L-equivalent of **58–72 µm**, midpoint **64.9 µm** (`pk_anchors.half_distance_to_efold_length`). `calibrate_alpha` now takes its target as an **e-folding length** and defaults to that converted midpoint.

> **Correction notice.** An earlier revision of this file reported **alpha = 309.6**, obtained by fitting to the raw 45 µm half-distance *as if it were* an e-folding length. That is the wrong convention, and it inflates alpha by exactly `1/ln²(2) = 2.081×`. The corrected, e-fold-converted fit gives **alpha = 148.6** — the value `results/pk_multidrug.md` already used, so the two files are now consistent. Any figure or table still quoting α = 309.6 predates this correction.

---

## The single calibrated constant

Fit on doxorubicin against the **64.9 µm e-fold-equivalent** of the 40–50 µm half-distance anchor (Primeau et al., Clin Cancer Res 2005; 11(24):8782–8788, PMID 16361566):

> **alpha = 148.6**   (dimensionless; a ~150× amplification of the structure-only binding prior)

This is the **only** fit. Everything below uses this same alpha.

### In-sample check (doxorubicin) — this is the fit, not a test

Doxorubicin post-calibration `L = 64.9 µm` vs the 58–72 µm e-fold-equivalent band → hits ~65 µm **by construction**.

---

## The real test: out-of-sample second drug (NOT fit)

With alpha **fixed at 148.6**, predict **mitoxantrone** and compare to an independent, same-method literature measurement — converted to an e-folding length with the same `d_half/ln2` rule used for the fit:

> Tannock lab, 'Distribution of the anticancer drugs doxorubicin, mitoxantrone and topotecan in tumors and normal tissues,' Cancer Chemother Pharmacol 2013; 72(3):591-601, PMID 23680920 (mitoxantrone intensity falls to half within ~25-75 um of tumor blood vessels; e-fold equivalent ~36-108 um)

| drug (out-of-sample) | predicted L | published half-distance | e-fold-equivalent band | ratio to mid | verdict |
|---|---:|---:|---:|---:|:--:|
| mitoxantrone | **48.2 µm** | 25–75 µm | 36–108 µm | 0.67× | **PASS** |

The predicted penetration falls **inside** the 36–108 µm e-fold-equivalent band — the one constant transfers to a second drug it was never fit to.

### Opposite-regime qualitative check (temozolomide)

Temozolomide (small, weakly bound) is predicted at `L = 265 µm`, **4.1× further** than doxorubicin (`65 µm`) — correct direction (the model captures the weakly-bound regime it was not fit on).

---

## Honest verdict

Single fitted constant alpha=148.6 (the ~150x amplification of the structure-only binding prior up to the true perivascular sink), fit on the E-FOLD-converted doxorubicin anchor L=d_half/ln2=64.9 um. In-sample: doxorubicin post-calibration L=64.9 um vs the 58-72 um e-fold-equivalent of the published 40-50 um half-distance (hits the target by construction -- this is the fit, not a test). OUT-OF-SAMPLE: mitoxantrone PREDICTED L=48.2 um vs the e-fold-equivalent 36-108 um (published half-distance 25-75 um; 0.67x the ~72 um midpoint) => PASS (inside the published band). Opposite regime: temozolomide L=265 um is 4.1x further than doxorubicin (65 um) -- correct direction. VERDICT: the one mechanistic constant GENERALIZES -- an independent, same-method measurement of a second drug falls where the model predicts, and the weakly-bound regime is captured qualitatively. HONESTY CAVEAT: mitoxantrone is an anthracenedione intercalator chemically related to doxorubicin, so this is a strong magnitude-transfer test within the strongly-bound class rather than a maximally-independent one; the temozolomide contrast is what probes the opposite regime.

### Caveats (the honesty line)

- **One anchor sets the scale.** alpha is fit to a single doxorubicin number; the absolute penetration of *every* drug rides on that one measurement.
- **The anchor is a half-distance, converted.** `L = d_half/ln2` assumes the perivascular profile is a single exponential. Real immunofluorescence profiles are only approximately exponential, so the conversion carries its own (modest) error — but applying it is strictly better than the previous convention, which silently equated the two lengths and mis-scaled alpha by 2.08×.
- **Mitoxantrone is chemically related to doxorubicin** (both anthracenedione DNA-intercalators). The measurement is genuinely independent, but this is a strong *magnitude-transfer* test within the strongly-bound class, not a maximally-independent one. The temozolomide contrast is what probes the opposite (weakly-bound, deep-penetration) regime — and only qualitatively, since temozolomide has no clean perivascular-distance dataset.
- **The proxy is an ordering, not a first-principles rate.** `tissue_binding_rate` is a low-confidence logP/HBD heuristic; alpha corrects its scale but inherits its ranking errors. A drug the heuristic mis-ranks (e.g. active transport, efflux, ion-trapping) would be mis-predicted.
- **Not fit to the second drug.** No parameter here was tuned to mitoxantrone or temozolomide; both are pure predictions.
