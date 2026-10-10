# PRE-REGISTRATION — scoring the ECM/diffusivity coupling against external measurement

Written 2026-10-10, **BEFORE any engine number was computed.** Frozen; not edited after the run.

## Honesty declaration about blindness (read first)

This is **not a blind test.** `results/transport_falsifier_scouting_2026-10-09.md`
§2.1 already published the Mapping-1 implied coefficients (1.50 / 1.64 / 2.16)
and §2.4 already published the `[0.40, 1.00]` band claim. Their *direction* was
therefore known to me before I wrote this file. What is genuinely pre-registered
here, and had NOT been computed by anyone when this file was written, is:

* **C1b** — the mapping-sensitivity surface (does the Candidate-1 verdict survive
  a change of unit mapping, or invert?).
* **C2b** — the comparison of Shi's ratios against the engine's **full hindrance
  chain** rather than against the spatial factor alone.

Criteria C1a and C2a are confirmatory re-derivations of a published claim, and
are labelled as such. A confirmatory re-derivation is not a pass; it is a check
on someone else's arithmetic.

## Hard constraint, restated

No parameter is tuned. `ALPHA_TORTUOSITY` stays 0.6. Nothing under `tmesim/**` is
edited. If the law fails, the output is a recorded failure. If an implied
coefficient comes out at 1.50, that number is **reported, never installed**.

---

## The shipped law, as read from code (verified by me, 2026-10-10)

| fact | file:line | verbatim |
|---|---|---|
| ECM field is max-normalised per section | `tmesim/reconstruct/real_fields.py:112-114` | `dmax = float(density.max())` / `if dmax > 0.0:` / `density /= dmax` |
| the coupling | `tmesim/physics/diffusion.py:88` | `d = drug.diffusivity * (1.0 - ALPHA_TORTUOSITY * ecm)` |
| the constant | `tmesim/physics/diffusion.py:75` | `ALPHA_TORTUOSITY = 0.6` |
| the floor | `tmesim/physics/diffusion.py:89` | `floor = max(1e-6, 1e-3 * drug.diffusivity)` |
| matrix volume fraction, hard-wired | `tmesim/molecular/ecm_mesh.py:153-156` | `fiber_radius_nm: float = 2.0, fiber_volume_fraction: float = 0.05` |
| called with radius only | `tmesim/molecular/compile_drug.py:174` | `mesh = ecm_mesh.ogston_renkin_hindrance(radius_nm)` |

Consequence, re-derived: `ecm` is clipped to `[0,1]` and attains exactly `1.0` at
the densest pixel of *every* section, so **`D_eff/D` spans `[0.40, 1.00]` by
construction**, independent of tissue.

---

## CRITERION 1 — Candidate 1 (Shafieyan 2014), FORM level

**Anchor (verified verbatim from PMC4104049 full text before this file was written):**
> "Diffusivity increased from 400 ± 123 to 778 ± 111, 110 ± 19 to 155 ± 14, and
> 124 ± 28 to 195 ± 26 μm²/s when fluid volume fraction increased from 0.66 to
> 0.97, 0.78 to 0.94, and 0.80 to 0.95 for sodium iodide, sodium diatrizoate, and
> Gd-DTPA, respectively (Fig. 6)."

### The unit mapping is the crux, and it is declared BEFORE the run

Shafieyan's `φ = V_f/(V_f + V_s)` (their Eq. 5) is a **true, gravimetrically
measured volume fraction**. The engine's `ecm_density` is a **max-normalised,
dimensionless, section-relative KDE index over matrix transcript counts**. These
are not the same quantity and no measurement fixes the correspondence. I
therefore do not pick one mapping; I declare **three**, each defensible, and
pre-register the sensitivity of the verdict to the choice.

Let `x = 1 - φ_fluid` (solid volume fraction) and `ecm = A·x`. The intercept is
**not** free: the shipped law gives `D_eff = D` exactly at `ecm = 0`, so `B = 0`
is forced. Only the scale `A` is free, and `A` is a property of the *mapping*,
not of the engine.

| mapping | `ecm = 1` denotes | `A` |
|---|---|---|
| **M1** | pure solid, `φ_solid = 1` | `1` |
| **M2** | the densest explant in that solute's own series (mirrors `density /= dmax`) | `1/x_max` |
| **M3** | the engine's own declared matrix volume fraction, `φ_f = 0.05` | `1/0.05 = 20` |

**M2 is the mapping that is faithful to the code**, because `real_fields.py:114`
max-normalises. M1 is the scouting record's mapping. M3 takes the engine's own
hard-wired `φ_f = 0.05` at its word.

**What each mapping assumes.** M1 assumes the engine's index is the solid volume
fraction itself, i.e. that a section's densest pixel is solid free of water —
which `results/heldout_patient.md` §3.5 refutes empirically (mean ECM 0.076 vs
0.457 on the *same* tumour type, driven by panel gene content). M2 assumes the
densest cartilage explant and the densest GBM pixel are the same physical
object, which is unverifiable. M3 assumes the engine's `φ_f = 0.05` is a real
matrix volume fraction rather than a convenient default, and that `ecm = 1`
saturates there.

### Pre-declared computation

For each solute, solve for the coefficient `α` the shipped form would need to
reproduce the **measured ratio of the two reported endpoints**:

    D(x_lo)/D(x_hi) = (1 - α·A·x_lo) / (1 - α·A·x_hi)

The ratio is used, not the absolute values, because absolute `D` requires a
free-water reference I have not verified from primary sources.

### Pre-declared criteria and what each outcome means

* **C1a (confirmatory).** Under **M1**, are the implied coefficients 1.50 / 1.64
  / 2.16 (±0.02)? **PASS** = the scouting record's arithmetic is right.
  This checks arithmetic only, not the physics.
* **C1b (the real test, never computed before).** Compute the implied `α` under
  M1, M2 and M3. **Pre-declared decision rule:**
  * If **all three** mappings give `α > 1` → the measured slope is inadmissible
    in the shipped functional form under every defensible mapping. **The linear
    form is REFUTED.**
  * If **all three** give `α < 1` → the slope is admissible; the test degrades
    from a form refutation to an amplitude comparison, and I report the factor by
    which the shipped `0.6` is off.
  * If the mappings **disagree in verdict** → **the refutation is an artefact of
    the unit mapping**, the comparison cannot be made defensibly at the form
    level, and *that* is the recorded finding. I pre-commit to reporting this
    outcome as a NON-RESULT for the form claim, not as a pass for either side.
* **C1c.** Does the paper itself report a tortuosity coefficient? **Pre-declared:**
  if the word "tortuosity" does not appear in the full text, then 1.50/1.64/2.16
  are **derived by inverting the engine's own law**, not measured, and any text
  calling them "reported" tortuosity coefficients must be corrected.

**Predicted direction, stated before running:** I expect C1b to **disagree across
mappings** (M1 `α > 1`, M2 `α < 1`), because `α` scales linearly with `A` and the
three `A` values span a factor of 20. If so, Candidate 1 does not refute the form.

---

## CRITERION 2 — Candidate 2 (Shi 2014), BAND level

**Anchor (verified verbatim from the PubMed abstract before this file was written):**
> "The mean Deff values were 33.2, 4.4, 1.3, 0.89, 0.59, and 0.47 × 10-7 cm2/s,
> respectively, for sodium fluorescein, dextran-4k, -20k, -40k, -70k, and IgG.
> The corresponding mean ratio of Deff to the free diffusion coefficient Dfree,
> Deff/Dfree, were 0.46, 0.19, 0.12, 0.12, 0.11, and 0.11 for these solutes."

### The denominator is the crux here, and it is declared BEFORE the run

Shi's denominator is **free-water `D_free`**. The engine's `[0.40, 1.00]` band is
`D_eff / drug.diffusivity`, where `drug.diffusivity` is the **already-hindered
compiled** value — the output of Stokes–Einstein × Brownian micro-tortuosity ×
Ogston–Renkin gate × two-route cell term. The two ratios have **different
denominators**, so comparing Shi's `D_eff/D_free` against the spatial band alone
is a units error. The comparable engine quantity is the **full chain**:

    (D_eff/D_free)_engine = (compiled_D / free_D) × (1 - 0.6·ecm)

* **C2a (confirmatory).** Are four of Shi's six ratios outside `[0.40, 1.00]`?
  **PASS** = yes. Checks the band claim against the spatial factor only.
* **C2b (the real test, never computed before).** Compute the engine's full-chain
  `D_eff/D_free` interval from the engine's own pure functions, read-only.
  **Pre-declared criteria:**
  * **C2b-i:** does any Shi ratio fall inside the engine's full-chain interval?
  * **C2b-ii:** what is the signed factor between the engine's full-chain value
    at the flagship mean `ecm = 0.076` and Shi's small-solute 0.46?
  * **Pre-declared decision rule on SIGN:** the scouting record claims the engine
    is "~2.1× too **fast**". If the full chain puts the engine **below** 0.46 by
    more than 1.2×, the published sign is **wrong** and the scouting record's
    error-partition ("~2.1× in D_eff, ~2.9× in k") is void, because the repo's own
    `D_eff/k` ≈ 6× too large (`results/kyle2004_penetration.md`) would then
    require `k` to be too **small** by much more than 2.9×, not too large.
* **C2c.** Shi's tissue is **normal rat brain**, not tumour; his small solute is
  sodium fluorescein (376 Da), not a repo compound. Any pass is cross-species and
  cross-tissue and must say so.

**Predicted direction, stated before running:** I expect the full chain to put the
engine's `D_eff/D_free` at **~0.03–0.08** (the scouting record's §0.2 table
implies a 12.7× total hindrance for doxorubicin, so `1/12.7 ≈ 0.079` times the
spatial factor). If so, the engine is **too SLOW**, not too fast, and the
published sign inverts.

---

## CRITERION 3 — counterweights (mandatory, not optional)

The literature does not agree on the sign of the tumour-vs-normal hindrance
effect. Pre-committed: the record must state, with verified verbatim quotes, that
Chary & Jain 1989 and Nugent & Jain 1984 both report tumour hindered **LESS**
than normal tissue — opposite in sign to Shi/Fu — and that Kutuzov 2018
impeaches Shi's radial-symmetry assumption. **A verdict quoting only one side is
declared invalid in advance.**

* **C3 PASS** = both sides are in the record with verified quotes, and the record
  states explicitly that the constant cannot be validated by picking one study.

---

## What CANNOT come out of this, by pre-commitment

* No new value of `ALPHA_TORTUOSITY`. Not computed, not suggested, not installed.
* No edit to `tmesim/**`.
* No digitisation of a figure by eye. If a number is figure-only, it is reported
  as figure-only and not used.
