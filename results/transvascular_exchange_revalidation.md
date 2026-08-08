# Transvascular exchange: the full re-validation

**One physics change, eleven independent re-validations.** This is the
consolidated record of what a corrected vascular delivery law does to every
published result in this repo. It is the evidence base for whether branch
`transvascular-exchange` can be merged, and for what has to be rewritten if
its flag is ever defaulted on.

**Nothing here is speculative.** Every number below was measured by re-running
the module in question both ways on real tissue. Every re-validation
reproduced its module's published baseline first, and no comparison was
trusted until the two arms were proven to genuinely differ (see §6 — the
verification trap fired for real, repeatedly).

---

## 1. The change

`tmesim/treatment/drugs.py::delivery_source` is a ONE-WAY leak::

    S = PS * Cp                      (no path back from interstitium to blood)

so the steady state is `PS*Cp / (clearance + binding)` — unbounded above the
plasma level. The physically correct form is gradient-driven exchange::

    S = PS * (Cp - C)  =  PS*Cp  -  PS*C

which is the SAME source minus a first-order loss at the SAME rate, so it
needs no new parameter. `PDEConfig(transvascular_exchange=True)` adds that
second half via `exchange_coefficient` and `simulate_drug_transport`'s
`exchange_k`. Steady state becomes `PS*Cp/(PS + clearance + binding) <= Cp`,
bounded by construction. Default is off and bit-for-bit identical.

**Why it was made:** `results/breast_kp_validation.md` — the model missed a
measured human microdialysis anchor by 12.8x too high (5-FU) and 10.5x too low
(methotrexate), in OPPOSITE directions, so no calibration constant could
absorb it.

## 2. How wrong the one-way law was

Not a breast-specific defect. On the GBM flagship, **5 of 11 drugs were
predicted above plasma concentration** — which passive delivery cannot do:

| drug | one-way | exchange |
|---|---:|---:|
| cisplatin | **2.7122** | 0.6184 |
| carmustine | **2.5360** | 0.6041 |
| lomustine | **1.3881** | 0.4584 |
| procarbazine | **1.2272** | 0.4374 |
| temozolomide | **1.1403** | 0.4979 |

Panel spread 386x → 90x on GBM; 970x → 104x on breast. Max cellwise
AUC/plasma-AUC **72.6 → 1.02**.

## 3. What SURVIVES

| result | verdict |
|---|---|
| **Numerics** | Closed-form steady state reproduced to **6.3e-15** worst case across 9 configs, `dt*PS` to 973. Mass budget closes to round-off. No new timestep limit (PS enters the implicit operator). Advective CFL trip point bit-identical. Non-negativity **observed** — the zero-clamp fired 0 times in 720 steps. |
| **GBM calibration chain** | **alpha unchanged (148.60815765445545).** Doxorubicin anchor holds; held-out mitoxantrone PASSES for every representative-PS choice, including a degenerate one. Correct analytic treatment is PS = 0, because PS and the source are the same field and the analytic length is defined only where it vanishes. Two agents reached this independently. |
| **GBM 3-D:2-D kill calibration** | **3/3 PASS**, none near a boundary. |
| **Breast 3-D:2-D kill calibration** | **5/6 PASS + docetaxel FAIL, zero flips.** Docetaxel moves the wrong way to be rescued. |
| **Immune-sink spatial spread** | RAW ratio 2.0025x → 1.9794x (−1.15%), so the raw statistic is robust to the transvascular term. But the raw ratio is not the claim: `results/immune_sink_impact.md` §5.3 scores it against spatially-aware nulls that keep the macrophages' clumpiness — clumpiness-matched permutation 1.6581, 400 um block shuffle 1.8522 — leaving an excess of **1.08-1.21x**, above the entire range of both in 20 draws (p = 0.0476, the floor). The 20-permutation LABEL null this row originally cited is demoted as too weak: it under-clumps the sink field it is compared against (CV 0.83 vs the observed 1.34). Monotonicity across the quintiles is NOT evidence and is dropped — every control reproduces it, including uniform-in-tissue scatter. Those nulls were run with exchange OFF and have not been recomputed with it ON. The mesh swing is larger than the physics swing. |
| **Dose-schedule structure** | Every schedule conclusion survives with its number roughly halved. Bound-vs-measured agreement holds (0.3387% vs 0.3792% ON, same relation as OFF). |
| **Mechanism reconciliation** | Bit-exact against `kill_fraction` in all 10 arms, 0 unattributed cells. |

### 3.1 Two places it makes the model BETTER

* **GBM 3-D:2-D delivery/state split moves toward the literature.** The report
  itself flagged as an unresolved defect that the model over-attributes the
  3-D penalty to delivery. Under exchange: temozolomide 56% → 38%, cisplatin
  54% → 27%, doxorubicin 62% → 62%. Delivery stays dominant **only for
  doxorubicin** — precisely Durand 1986's published ordering (PMID 2425117,
  "penetration was a problem only for doxorubicin").
* **Partition scorecard improves on the modality-matched rows.** Scoring only
  the microdialysis (free-drug) rows, which are the ones commensurable with
  what the engine computes: GM |fold error| **5.235x → 3.969x**, 3x band 0/3 →
  1/3, 10x band 3/3 → 3/3. No regression. The apparent regression in the
  pooled score is entirely one tissue-homogenate row (adavosertib, measured
  2.46) — homogenate measures TOTAL drug and can exceed 1 by definition, while
  the engine computes FREE interstitial, which cannot. `glioma_regional_pk.md`
  §42 already states that rule. Scoring against it penalises the model for
  being right.

## 4. What BREAKS

### 4.1 "Ranking is invariant" is FALSE as written — the headline casualty

| GBM, kill_fraction | Spearman | Kendall |
|---|---:|---:|
| off vs on | +0.9273 | +0.8222 |

High correlation, and it does not save the claim. **Top-1 changes
(cisplatin → vincristine). Top-3 set changes.** Four of 45 pairs swap, and two
of them cross gaps of **40.5% and 35.6% of the leader's kill** against a
Bernoulli sd of ~0.003 — roughly 50-sigma reversals, not near-ties.

The Monte-Carlo ranking panel is worse: **P(cisplatin rank 1) 0.9277 →
0.3965**, P(doxorubicin rank 1) 0.0430 → 0.5586, modal order reverses. The
module's own pre-declared "decision-grade" bar (leader above 0.90, runner-up
below 0.10) is met OFF and destroyed ON.

On breast, top-1 survives (vinorelbine) but the top-3 set changes and 7 of 45
pairs are discordant.

**Why no constant can absorb it.** The attenuation is
`(clearance+binding)/(PS + clearance+binding)`, and PS is driven by
`bbb_penetrance`, which spans **47x** across the GBM panel. AUC ON/OFF runs
from 0.228 (cisplatin) to 0.975 (vincristine). This is a **differential
re-weighting by permeability**, not a global rescale. The ranking's historical
robustness — to alpha, to the 3.36x KDE→discrete source change, to the
cell-state route — does not generalise, because every one of those
perturbations was **common-mode** and this one is **drug-specific**.

**Delivery-rank invariance does not imply kill-rank invariance.** On breast,
AUC ranking is perfectly preserved (rho = 1.000) while kill ranking scrambles
(rho = 0.830), because the Hill map saturates: drugs on the plateau lose
almost no kill when AUC drops, drugs on the steep limb lose most of theirs.
Any future claim must say WHICH quantity it is about.

### 4.2 Cisplatin's mechanism verdict inverts — in four independent code paths

Published in two places as delivery-solved: `mechanism_attribution.md` §3
(RESISTANCE 96.42%) and `dose_schedule.md` (SENSITIVITY-limited, on 99.55%
IC50 coverage). Both rest on tissue sitting at 2.71x plasma. Bounded:

| cisplatin | one-way | exchange |
|---|---:|---:|
| survival: RESISTANCE | 96.42% | 3.59% |
| survival: DOSE | 0.44% | **95.71%** |
| IC50 coverage at 1x | 99.55% | **4.67%** |

**"Cisplatin is a biology problem" becomes "cisplatin is a dose problem."**
Reproduced independently by the mechanism, 3-D:2-D, ranking, and dose-schedule
re-runs. Both reports now carry conditional-retraction notices.

### 4.3 Invalidated rather than rescaled

* Every absolute AUC and kill (temozolomide −56%, cisplatin −77%,
  doxorubicin −3.5%).
* `dose_schedule.md` §3's whole delivery-geometry table. Exchange removes drug
  fastest exactly where the source is strongest, so it FLATTENS the field:
  temozolomide penetration index 0.5112 → 0.7140 (**+40%**). The qualitative
  claims survive; the numbers do not.
* "At high enough dose >95% of tumour cells clear IC50 for both drugs" fails
  for temozolomide: 30x coverage 97.80% → **3.43%**.
* The landing page's quoted `delivery 0.172, dose 0.825, resistance 0.003`
  (run `r_00067`) becomes `0.423 / 0.577 / 0.000`.

## 5. New limits the fix does not remove, and one it creates

* **It fixes the FORM of the delivery law, not the VALUE of the
  permeability.** PS is still `bbb_penetrance`, a blood-BRAIN-barrier term.
  Methotrexate does not improve on breast (0.10x → 0.09x of measured) because
  its defect is the wrong permeability, not the wrong law. Two independent
  defects; one fixed.
* **PS and `bbb_penetrance` are perfectly confounded.** Scaling
  `vascular_density` by f gives a BIT-IDENTICAL PS field to scaling
  `bbb_penetrance` by f. The model cannot distinguish a wrong vessel map from
  a wrong permeability. Separating them needs an independent absolute
  measurement this project does not have.
* **Sensitivity improves, unevenly, and "robust" can mean "degenerate".**
  Elasticity of AUC to PS is exactly **1.0000** under the one-way form for
  every drug, and drops to 0.085–0.966 under exchange — never rises. But
  temozolomide and cisplatin become robust only because they SATURATE
  (Kp → 0.91–0.96, i.e. ~= Cp): the model is now saying "equilibrates with
  blood", which is stable and carries little spatial information. Trust the
  robustness, not the informativeness. The fragile half (palbociclib,
  irinotecan, etoposide, methotrexate, vincristine, elasticity 0.74–0.99) gets
  almost no benefit.
* **~~NEW RESOLUTION REQUIREMENT~~ — MEASURED, AND THIS ENTRY WAS WRONG.**
  Exchange does introduce a perivascular length scale
  `sqrt(D_eff/(PS + clearance + binding))` that collapses to ~15-17 um for
  `bbb_penetrance` >~ 0.6, and this section originally predicted, from that,
  that dx = 12 um would be inadequate and that finer grids would push the
  ON-arm values **DOWN** (reasoning: an under-resolved boundary layer
  over-estimates tissue concentration, which would DEEPEN the §4.1
  reordering).

  **Both halves of that prediction are false.** Swept dx = 12 / 10 / 8 um on
  the real GBM section, exchange on, for exactly the flagged drugs, with
  doxorubicin as a low-permeability control:

  | drug | bbb_pen | ratio@12 | ratio@10 | ratio@8 | 12→8 |
  |---|---:|---:|---:|---:|---:|
  | lomustine | 0.973 | 0.4584 | 0.4623 | 0.4671 | **+1.9%** |
  | carmustine | 0.955 | 0.6041 | 0.6080 | 0.6123 | **+1.4%** |
  | procarbazine | 0.885 | 0.4374 | 0.4407 | 0.4453 | **+1.8%** |
  | ribociclib | 0.800 | 0.2609 | 0.2636 | 0.2677 | **+2.6%** |
  | cisplatin | 0.858 | 0.6184 | 0.6223 | 0.6266 | **+1.3%** |
  | doxorubicin (control) | 0.021 | 0.0165 | 0.0163 | 0.0165 | −0.2% |

  Values move **UP**, not down, and by only **1.3-2.6%** — so dx = 12 um is
  adequate after all and this entry overstated the risk. The doxorubicin
  control is flat (−0.2% ratio, −0.0% kill), confirming the sweep measures the
  boundary layer rather than a generic grid artefact.

  **The §4.1 reordering is unaffected**, though not for the reason originally
  given. Kills move +1.0% to +2.0%, and the decisive cisplatin/doxorubicin
  swap persists at dx = 8 (cisplatin 0.2335 vs doxorubicin 0.2547, versus
  0.2311 vs 0.2547 at dx = 12). The reordering is not a resolution artefact.

  Retained rather than deleted because the original reasoning was recorded as
  a falsifiable prediction and it failed; deleting it would erase the check.
* **Attribution degenerates to two categories.** When exposure drops below
  IC50 population-wide, RESISTANCE and CHANCE go to exactly zero cells because
  the exposure bands take precedence. Under bounded physics that is most of
  the panel.
* **A quantisation hazard.** With the flag on, temozolomide's schedule spread
  prints `0.00e+00` — cosmetically identical to a solver artefact
  `dose_schedule.md` was rewritten to retract. It is not that artefact
  (continuous spread 0.3792%, monotone), but kill is quantised at 1/23,812 and
  the level fell 2.6x. Keep the up-dose ladder or a spurious zero gets
  republished.

## 6. The verification trap, which fired for real

Several modules call `simulate_drug_transport` DIRECTLY rather than through
`simulate_response`, so `PDEConfig(transvascular_exchange=True)` is a **silent
no-op** for them. Confirmed by explicit control, the cleanest being
`dose_schedule`:

| run | auc_sum |
|---|---:|
| exchange OFF | 1.554447793e+09 |
| **flag set, `exchange_k` NOT passed** | **1.554447793e+09 — bit-identical** |
| honest ON (`exchange_k` passed) | 6.86608671e+08 |

Affected: `pk_reality`, `kill_calibration`, `immune_sink_impact`,
`vessel_profile`, `numerical_rigor`, `experiments/regimen`,
`experiments/dose_schedule`. Had this not been guarded against, five of the
eleven re-validations would have silently reported "nothing changed".

### 6.1 A THIRD category: modules that touch no PDE at all

Beyond `simulate_response` (flag works) and direct `simulate_drug_transport`
(flag is a no-op, pass `exchange_k` instead), there is a third class that
never invokes the transport solver in any form, so there is no injection
point and `exchange_k` is not an available remedy either.

Confirmed for **`gist_imatinib`** (and its transitive dependencies
`pk_reality.penetration_length_um` and `calibration.calibrate_alpha`, both
closed-form). Its prediction is an analytic 1-D Krogh slab,
`E = (L/R)*tanh(R/L)`. Verified by forcing `PDEConfig.__init__.__defaults__`
so every config anywhere in the call tree carries the flag, and instrumenting
call counters: `simulate_drug_transport` 0 calls, `delivery_source` 0,
`exchange_coefficient` 0, and the entire result payload **byte-identical**
between arms. A positive control under the same patch on a
`simulate_response` path moved mean AUC 14.1x, proving the forcing mechanism
was live rather than inert.

**These modules must be reported as UN-TESTABLE, not as unaffected.** "Robust
under exchange" and "never evaluates the delivery law" are very different
statements and must not be pooled in a coverage count. `glioma_regional_pk`
is also analytic (assumption A1 is a closed-form structure score), so it
belongs here too.

**Why GIST specifically cannot see the defect**, which generalises: its target
is a RATIO, `C_tumour / C_normal`. The plasma term cancels exactly, so any
breach of `C <= Cp` present in both numerator and denominator divides out.
That is precisely why 5 of 11 GBM drugs and 3 of 10 breast drugs were caught
and imatinib was not. Computed off-module, imatinib **would** show the same
pathology — 11.71x plasma uncalibrated, 1.96x calibrated, at the vascular
densities the repo's own fixtures use.

**And a reframe worth keeping.** Imatinib's `bbb_penetrance` is 0.5237,
ranking 8th of 13 — above temozolomide, cisplatin and doxorubicin, and 23x
methotrexate. As a blood-brain-barrier number that is wrong by one to two
orders of magnitude (imatinib's CSF:plasma is ~0.01-0.1; it is a P-gp
substrate and the score is structure-only). But for GIST, where imatinib
genuinely does distribute well into peripheral tissue (Vd ~295 L), a high
value is arguably closer to right. The defect is that the parameter is
**MISLABELLED** — a brain-specific quantity used as a general tissue
permeability — not simply mis-valued. That distinction matters for how it
should be fixed: the fix is a tissue-aware permeability, not a global rescale
of `bbb_penetrance`.

## 6.2 A THIRD DEFECT, found while trying to verify a lead that was itself wrong

A follow-up asked whether `bbb_penetrance` is miscalibrated in absolute
magnitude in brain — the engine's volume-averaged methotrexate PS is
6.039e-4 s^-1, against a directly measured intact-BBB value of 3.3e-5 s^-1
(Ohno K, **Fredericks WR**, Rapoport SI, *Surg Neurol* 1979;12(4):323-8, PMID
524251 — note the author list, an earlier sweep of ours misattributed this to
"Ohno, Pettigrew & Rapoport", which is a different paper, PMID 696840).

**That claim is REFUTED, and the refutation is more interesting than the
claim.** `bbb_penetrance` is NOT a permeability. It is a dimensionless
brain:plasma PARTITION coefficient, and it says so in three places
(`core/config.py:25`, `reconstruct/barrier.py:66` — "It is NOT a permeability
with physical units" — and `molecular/permeability.py`, whose docstring
calibrates it against measured tissue:plasma ratios). Judged in the units it
was built in, methotrexate at 0.0225 against a measured 0.178 is if anything
**low**, not 18x high.

**But the engine uses it as a first-order RATE.** `treatment/drugs.py:228`
declares `exchange_coefficient` returns "the permeability-surface-area field
PS(x, y) in 1/s"; `sim/engine.py:44` declares
`source = vascular_density * bbb_penetrance * plasma_mean (uM/s)`. The
dimensional audit:

    vascular_density [dimensionless] x bbb_penetrance [dimensionless]
        x plasma [uM]  =  uM,  NOT uM/s

The engine obtains `s^-1` by an **implicit, undocumented, uncalibrated
conversion factor of exactly 1 s^-1**. Nothing in the repo derives or
justifies that scale.

**The repo already contradicts itself about which quantity this is.**
`results/kill_calibration.md:119` insists the term is "a **rate**, not a
concentration... **must not be read as** a physical tissue:plasma partition
coefficient". `results/glioma_regional_pk.md:34` reads it as exactly that
partition. Both cannot be true.

**The diagnostic that settles it.** Read as RATES, the drugs are mutually
incoherent: methotrexate is 18.3x above measured, temozolomide 2.10x — an
**8.7x spread** between two drugs against the same normal rat brain. Read as
PARTITIONS with the engine's own barrier term on, methotrexate is 0.79x and
temozolomide 0.82x of measured — a **1.04x spread**. The drug-to-drug
incoherence exists ONLY under the rate reading. That is the signature of a
units error, not a calibration error.

**It does NOT cancel from rankings.** Measured temozolomide:methotrexate
permeability ratio is 48.5x; the engine's is 5.57x, so it **under-separates by
8.7x** (5.1x after correcting for protein binding). Note this runs OPPOSITE to
the breast finding, where the engine OVER-separates methotrexate:5-FU by 28.1x
against a measured 1.02x. The parameter is unreliable in **both** directions,
which no single rescale can fix. And because `Kp = PS/(PS + k_loss)` saturates
in PS, no scale factor can convert one reading into the other — the distortion
is necessarily drug-dependent.

**This defect is structurally invisible to the existing GBM validation, and
more completely than §5 claimed.** The penetration-length anchor is
`L = sqrt(D_eff / k)` with `k = clearance + binding`. With
`transvascular_exchange = False` (the default), **PS appears nowhere in that
formula at all.** The doxorubicin-fit → mitoxantrone-predicted chain could not
have caught this at any magnitude of error.

**Consequence for the exchange term added in this branch.** It inherits the
same implicit scale: `exchange_k` is that dimensionless product, added
directly to a reaction rate in `1/s`. The term is at least now *consistent* —
the same quantity drives influx and back-flux, so whatever the true scale
`k_0` is, it multiplies both halves — but the steady state
`k_0*PS'*Cp / (k_0*PS' + k)` depends on `k_0`, and `k_0 = 1 s^-1` is
unjustified. **A calibrated `k_0` would change every exchange-on number in
this document.** They should be read as conditional on that assumption.

Corroboration for the measured side is solid — three groups, three methods, a
decade apart: Ohno 1979 (PA 3.3e-5 s^-1 intact, 28e-5 opened), Shapiro 1988
(PMID 3335031, methotrexate K1 1.3-1.9 uL/g/min cortex, 1.3-8.2 in C6 glioma),
Groothuis 1983 (PMID 6850641). The "but glioma is leaky" defence fails on
direct measurement: Shapiro measured methotrexate itself in a rat glioma at
4.4-28x BELOW the engine.

**Recommended fix, in order:** (a) decide which quantity `bbb_penetrance` is
and make `config.py`, `barrier.py`, `drugs.py`, `engine.py`,
`kill_calibration.md` and `glioma_regional_pk.md` agree — they currently do
not; (b) if the rate reading is kept, introduce an explicit calibrated
`PS = k_0 * vascular_density * bbb_penetrance` with `k_0` fit to measured PS
values, not left at an implicit 1.

## 7. Merge position

**The physics is right and the implementation is correct.** It is verified to
machine precision, it removes a defect that put half the GBM panel at
physically unavailable concentrations, it leaves the calibration chain intact,
and it moves two independent benchmarks toward measured reality.

**It is not merge-ready as a DEFAULT**, for one reason: it changes the top-1
drug recommendation, which is the decision this product exists to support, and
the claim "ranking is invariant" has to be rewritten across the repo and the
landing page before that ships. Regenerating every affected results doc is a
prerequisite, not a follow-up.

**Recommended sequencing:** keep the flag opt-in and default-off (as it is);
rewrite the ranking claim now (done for the landing page in this branch);
re-run the high-permeability drugs at dx = 8-10 um to converge the §4.1
magnitudes; then regenerate the affected reports and flip the default in a
single commit that also updates the marketing copy.

## 8. Provenance

Eleven parallel re-validations, each reproducing its module's published
baseline before comparing, each proving its two arms differed. Modules
covered: calibration chain, glioma partition scorecard, GBM 3-D:2-D, breast
3-D:2-D, mechanism attribution, numerical rigour, ranking robustness,
dose schedule + time-resolved PK, PS sensitivity, immune-sink, and a
cross-panel level/ranking sweep. Sources: `experiments/breast_kp_anchor.py`,
`results/breast_kp_validation.md`, `results/breast_transport_anchors.md`, and
the conditional-retraction notices now carried in
`results/mechanism_attribution.md` and `results/dose_schedule.md`.
