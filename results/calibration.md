# Calibration — one mechanistic constant for perivascular drug consumption

**Module:** `tmesim/validation/calibration.py` (run: `python -m tmesim.validation.calibration`)

**Scope:** fixes the calibration gap the PK reality check exposed (tmesim over-predicted doxorubicin tissue penetration ~15x: ~700 µm vs the published ~40–50 µm). We make the tissue loss rate *mechanistic* with exactly **one** free constant, fit it **once** on doxorubicin, then test whether it predicts a **second** drug with **no** refitting. Nothing in core, the engine, or the molecular modules is modified — this module only drives them and returns a normal `core.DrugParams`.

---

## The mechanism and the proxy

The frozen transport PDE loses drug at a first-order rate `k = clearance + binding`, giving penetration length `L = sqrt(D_eff / k)`. The structure-only priors give `k ≈ 3e-4 /s` (a ~1-hour timescale) → `L ≈ 700 µm`. Reproducing a ~45 µm boundary layer needs `k = D_eff / L² ≈ 7e-2 /s` (a ~15-second timescale) — about **300× faster**. The missing physics is not slow metabolic clearance; it is **rapid cellular uptake + high-capacity intracellular / DNA binding**, a near-irreversible perivascular sink that is what actually carves the steep gradient.

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

## The single calibrated constant

Fit on doxorubicin against the **45 µm** anchor (Primeau et al., Clin Cancer Res 2005; 11(24):8782–8788, PMID 16361566):

> **alpha = 309.6**   (dimensionless; the ~300× amplification the diagnosis predicted)

This is the **only** fit. Everything below uses this same alpha.

### In-sample check (doxorubicin) — this is the fit, not a test

Doxorubicin post-calibration `L = 45.0 µm` vs the 40–50 µm anchor → hits ~45 µm **by construction**.

---

## The real test: out-of-sample second drug (NOT fit)

With alpha **fixed at 309.6**, predict **mitoxantrone** and compare to an independent, same-method literature measurement:

> Tannock lab, 'Distribution of the anticancer drugs doxorubicin, mitoxantrone and topotecan in tumors and normal tissues,' Cancer Chemother Pharmacol 2013; 72(3):591-601, PMID 23680920 (mitoxantrone intensity falls to half within ~25-75 um of tumor blood vessels)

| drug (out-of-sample) | predicted L | published band | ratio to mid | verdict |
|---|---:|---:|---:|:--:|
| mitoxantrone | **33.4 µm** | 25–75 µm | 0.67× | **PASS** |

The predicted penetration falls **inside** the published 25–75 µm band — the one constant transfers to a second drug it was never fit to.

### Opposite-regime qualitative check (temozolomide)

Temozolomide (small, weakly bound) is predicted at `L = 184 µm`, **4.1× further** than doxorubicin (`45 µm`) — correct direction (the model captures the weakly-bound regime it was not fit on).

---

## Honest verdict

Single fitted constant alpha=309.6 (the ~300x amplification of the structure-only binding prior up to the true perivascular sink). In-sample: doxorubicin post-calibration L=45.0 um vs 40-50 um anchor (hits ~45 um by construction -- this is the fit, not a test). OUT-OF-SAMPLE: mitoxantrone PREDICTED L=33.4 um vs published 25-75 um (0.67x the ~50 um midpoint) => PASS (inside the published band). Opposite regime: temozolomide L=184 um is 4.1x further than doxorubicin (45 um) -- correct direction. VERDICT: the one mechanistic constant GENERALIZES -- an independent, same-method measurement of a second drug falls where the model predicts, and the weakly-bound regime is captured qualitatively. HONESTY CAVEAT: mitoxantrone is an anthracenedione intercalator chemically related to doxorubicin, so this is a strong magnitude-transfer test within the strongly-bound class rather than a maximally-independent one; the temozolomide contrast is what probes the opposite regime.

### Caveats (the honesty line)

- **One anchor sets the scale.** alpha is fit to a single doxorubicin number; the absolute penetration of *every* drug rides on that one measurement.
- **Mitoxantrone is chemically related to doxorubicin** (both anthracenedione DNA-intercalators). The measurement is genuinely independent, but this is a strong *magnitude-transfer* test within the strongly-bound class, not a maximally-independent one. The temozolomide contrast is what probes the opposite (weakly-bound, deep-penetration) regime — and only qualitatively, since temozolomide has no clean perivascular-distance dataset.
- **The proxy is an ordering, not a first-principles rate.** `tissue_binding_rate` is a low-confidence logP/HBD heuristic; alpha corrects its scale but inherits its ranking errors. A drug the heuristic mis-ranks (e.g. active transport, efflux, ion-trapping) would be mis-predicted.
- **Not fit to the second drug.** No parameter here was tuned to mitoxantrone or temozolomide; both are pure predictions.
