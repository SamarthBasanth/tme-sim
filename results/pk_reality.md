# PK Reality Check — does tmesim's predicted drug penetration match published tissue PK?

**Date:** 2026-07-26
**Module:** `tmesim/validation/pk_reality.py` (run: `python -m tmesim.validation.pk_reality`)
**Scope:** first real credibility test of the molecular → tissue pipeline. It does **not**
modify any frozen core/engine/molecular module; it only *drives* them (compiles real drugs
with `molecular.compile_drug.molecule_to_drugparams`, runs the frozen transport PDE) and
compares to the literature.

---

## The literature anchor (with citation)

**Doxorubicin penetrates only ~40–50 µm from tumor blood vessels.** A famous, repeatedly
reproduced result: doxorubicin autofluorescence falls to background within ~40–50 µm of the
nearest perfused vessel, producing steep perivascular gradients and leaving cells more than a
few cell-diameters from a vessel effectively untreated.

> Primeau AJ, Rendon A, Hedley D, Lilge L, Tannock IF. *The distribution of the anticancer
> drug doxorubicin in relation to blood vessels in solid tumors.* **Clin Cancer Res 2005;
> 11(24):8782–8788.** PMID **16361566** (confirmed via PubMed).

Supporting anchors (used for the **ordering** check, not exact magnitude):

- **Reaction-diffusion penetration length** `L = sqrt(D_eff / k)`, where `k = clearance + binding`
  is the first-order tissue loss rate. This is the steady-state e-folding length of the
  perivascular profile `C(r) ~ K₀(r/L)` for a vessel source (Thiele-modulus / boundary-layer
  theory; e.g. Thurber, Schmidt, Wittrup, *Adv Drug Deliv Rev* 2008;60:1421–1434).
- **Small molecules reach further than large / strongly-bound / rapidly-consumed molecules.**
- **Antibodies penetrate poorly** — the "binding-site barrier" (Fujimori et al., *J Nucl Med*
  1990;31:1191–1198).

---

## What the sim predicts

Drugs compiled through the real molecular pipeline (`molecule_to_drugparams`, ECM fraction 0.4,
so `D_eff = D·(1 − 0.6·ECM)`). Two independent penetration metrics:

- **`L_analytic`** — closed-form `sqrt(D_eff/k)`; grid-independent, the principled number.
- **`prof 1/e`** — single-vessel solve through the **frozen** transport engine; distance
  outward from the vessel edge at which steady-state `C` falls to `1/e` of its near-vessel value.
- **`efold fit`** — exponential length fit to the far-field tail of that same solve; recovers `L`.

| drug         |     MW (Da) | D_eff (µm²/s) |    k (1/s) | L_analytic (µm) | prof 1/e (µm) | efold fit (µm) |   bbb |
|--------------|------------:|--------------:|-----------:|----------------:|--------------:|---------------:|------:|
| doxorubicin  |         544 |         181.7 |  2.87e-04  |         **796** |       **210** |            866 | 0.051 |
| temozolomide |         194 |         256.4 |  7.79e-05  |            1815 |           710 |           2816 | 0.558 |
| bevacizumab  |     149 000 |          19.2 |  3.24e-05  |             770 |           190 |            831 | 0.000 |

The far-field fit for doxorubicin (866 µm) matches `L_analytic` (796 µm) — confirming `L` is the
true decay length, and the shorter `1/e` distance (210 µm) is the 2D geometric-spreading
(`K₀ ~ exp(−r/L)/√r`) near-field effect, not a different physics.

---

## Comparison to the anchor (doxorubicin, mid ≈ 45 µm)

| metric                    | predicted | anchor | ratio  | band (±3×) |
|---------------------------|----------:|-------:|-------:|:----------:|
| analytic penetration `L`  |    796 µm |  45 µm | **17.7×** | **FAIL** |
| single-vessel `1/e`       |    210 µm |  45 µm |  **4.7×** | **FAIL** |

Both metrics land **outside** the order-of-magnitude (3×) band. Even the most conservative
reading (the `1/e` distance) over-predicts penetration by ~5×; the principled `L` over-predicts
by ~18×.

### Ordering checks (directional, structure-driven) — both PASS

- ✅ **temozolomide `L` (1815 µm) > doxorubicin `L` (796 µm)** — the smaller, less-bound drug
  reaches further. Correct direction.
- ✅ **small-molecule D_eff (TMZ 256 µm²/s) ≫ antibody D_eff (bevacizumab 19 µm²/s), ~13×** —
  the antibody diffuses far slower. Correct direction.

---

## Honest verdict

**MAGNITUDE MISMATCH (fails), directional physics correct (passes).**

The sim **over-predicts doxorubicin tissue penetration by roughly an order of magnitude**
(≈18× on the principled reaction-diffusion length, ≈5× even on the more forgiving single-vessel
`1/e` distance). It does **not** reproduce the ~40–50 µm perivascular gradient that is the
headline experimental fact about doxorubicin.

**Root cause (diagnosed, not fudged):** the penetration length is `sqrt(D_eff/k)`, and the
compiled loss rate `k = clearance + binding ≈ 2.9e-4 /s` corresponds to a ~1-hour timescale.
That is far too slow. Doxorubicin's real perivascular gradient is set by **rapid, high-capacity
DNA intercalation / tissue binding and cellular uptake** — an effective consumption orders of
magnitude faster than the Stage-0 structure-only priors in `molecular/binding.py` produce
(those priors are explicitly documented there as *low-confidence, order-of-magnitude*). To match
~45 µm at `D_eff ≈ 182 µm²/s` would require `k ≈ D_eff / L² ≈ 182 / 45² ≈ 0.09 /s` — about **300×**
larger than the compiled value. The mismatch is a *binding/clearance-calibration* gap, not a
diffusion or numerics error (the far-field fit cleanly recovers the analytic `L`, and the
transport engine's mass-conservation self-tests pass).

**What works:** the *shape* of the physics is right. Larger/stickier drugs get lower `D_eff`
and shorter `L`; antibodies diffuse ~13× slower than small molecules. The relative ordering
tmesim relies on for its Stage-0 ranking is preserved. What is missing is any absolute
calibration of tissue loss rates against measured PK.

### Caveats

- **Uncalibrated parameters.** Binding & clearance are structure-only heuristics flagged
  low-confidence in `molecular/binding.py`; this check *is* the evidence that they need
  PK calibration. IC50/potency (not tested here) is even weaker.
- **Antibody penetration is under-modeled.** The sim gets the antibody right *via low diffusivity*
  but its `L` (770 µm) is ~ doxorubicin's, because the **binding-site barrier is not encoded**
  (the compiled antibody binding rate is tiny). Real antibodies penetrate *worse* than
  doxorubicin, not similarly. Ordering by `D_eff` passes; ordering by `L` would not capture this.
- **2D, single-vessel idealization.** The profile solve is a homogeneous 2D patch with one vessel,
  no interstitial convection, uniform ECM — an idealization of the real 3D, heterogeneous,
  pressure-driven microenvironment. Near a small source, 2D geometric spreading shortens the
  `1/e` distance relative to `L`.
- **Anchor is one drug/model.** 40–50 µm is doxorubicin in specific xenograft/tumor models;
  penetration varies with drug, dose, and tumor. The 3× band is deliberately generous.

**Bottom line:** tmesim's transport *machinery* is sound and its qualitative drug ordering is
correct, but its absolute perivascular penetration for doxorubicin is ~1 order of magnitude too
optimistic because the tissue loss (binding/clearance) rates are uncalibrated. This is a real,
useful, honestly-reported failure — the first quantitative target for PK calibration.
