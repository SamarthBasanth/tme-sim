# Does the compound ORDERING survive a tissue-capacity term?

**A sensitivity study on the shipping panels. No model change.** Nothing in
`tmesim/` was modified, no default moved, `calibrate_alpha()` is still
`22.223477711172563` (asserted against live code at every stage), and every
arm below is the unmodified engine driven with a per-run drug object. At
`R = 1` the arm IS the shipping arm, bit-for-bit.

* script: `experiments/capacity_ranking_robustness.py`
  (`--stage algebra | sweep | breast | analyse | gridcheck | report`)
* frozen execution record: `results/sweep/capacity_ranking_robustness.json`
* the finding this answers: `results/brain_slice_anchor.md` (TEST B, **FAIL**),
  registered in `results/transport_anchor_prereg.md` sec. 3.5/3.5a and in the
  `brain-slice-tissue-capacity` row of `tmesim/validation/scorecard.py`.
* ranking-stability conventions are `results/ranking_robustness.md`'s, reused
  unchanged; the one deviation is named in sec. 5.1.

---

## 0. The answer

> ## NO.
>
> **THE TOP-1 PICK CHANGES — say this first.** On the shipping **breast 10-drug**
> panel, giving each compound the capacity implied by this repo's OWN
> structural `K_b` rule (`two_route.binding_capacity_from_logp`) replaces the
> top pick `vinorelbine` with `doxorubicin`. It changes under the
> access-limited variant of the same rule as well. That is 2 of the
> named per-compound assignments tested, not a tail of a sampling
> distribution.
>
> **The full published ordering is not robust, and it fails in the weakest
> possible case.** Under a **uniform** capacity — every compound given the
> SAME value, which is exactly what the "a constant error cancels out of a
> ranking" defence covers — the published order breaks at:
>
> * **GBM 10-drug: R = 36.3** — below the best-available constant of 39.
> * **GBM 5-drug: R = 435** — above the best-available constant of 39.
> * **breast 10-drug: R = 5.28** — below the best-available constant of 39, and only 1.7x above the bottom of the measured range.
>
> **Per-compound capacity breaks the top-1 claim too.** Resampling each
> compound's capacity independently from the 33 measured values, NO top-1
> claim on ANY of the panels meets `results/ranking_robustness.md`'s
> decision-grade bar (leader's 90 % interval entirely above 0.90, runner-up's
> entirely below 0.10).
>
> **The one result that cuts the other way, reported at the top rather than
> buried:** on the **five-drug GBM sub-panel** — the one
> `results/ranking_robustness.md` actually scores — EVERY named assignment
> (uniform R = 39, the structural rule, the access-limited rule) reproduces
> the published order **exactly**. And under a uniform capacity the GBM
> 10-drug top-1 (`vincristine`) never changes anywhere on the measured
> range. Those are genuine, and sec. 5.3 explains the mechanism rather than
> waving at it.
>
> | | top-1, UNIFORM capacity | top-1, PER-COMPOUND (sampled) | top-1, PER-COMPOUND (structural rule) | full order, UNIFORM |
> |---|---|---|---|---|
> | **GBM 10-drug** | `vincristine` holds to R = **&gt; 4184** (no change) | P(`vincristine`) = **0.730** [0.726, 0.734] — **not decision-grade** | `vincristine` (unchanged) | holds to R = **36.3** |
> | **GBM 5-drug** | `cisplatin` holds to R = **435** | P(`cisplatin`) = **0.584** [0.580, 0.588] — **not decision-grade** | `cisplatin` (unchanged) | holds to R = **435** |
> | **breast 10-drug** | `vinorelbine` holds to R = **&gt; 4184** (no change) | P(`vinorelbine`) = **0.797** [0.794, 0.801] — **not decision-grade** | `doxorubicin` (**CHANGED**) | holds to R = **5.28** |
>
> *Measured rat brain-slice `Vu,brain`: 3.2 – 4184 mL/g, n = 33, best available constant 39. The engine assumes 1.0 for every compound.*

**What changes on the marketing page.**

* **"The compound ranking is robust" cannot be said.** The published order
  breaks under a UNIFORM capacity on all 3 panels tested, and on the breast panel it breaks at R = 5.28, only 1.7x above the bottom of the measured range.
* **"Our top pick is robust" cannot be said unqualified either.** It survives a
  uniform capacity well, but a named per-compound assignment changes it on the
  breast panel and no sampled per-compound arm reaches decision-grade.
* **What CAN be said, and it is worth saying:** *"Codon's steady-state and
  equilibrium-referenced outputs — penetration length, interstitial:plasma
  ratio, converged gradient shape — are EXACTLY invariant to tissue capacity
  (sec. 2b, verified to 1e-9). The finite-window kill ranking is not, and the
  range over which it is stable is stated in sec. 9."* That is a narrower
  claim than the one being made today, and unlike it, it is demonstrated.

**What is NOT claimed here.** This is not evidence that a capacity term
*belongs* in the engine. `results/intracellular_retardation.md` refuted exactly
that term on a different dataset, and sec. 6.1 says so before anything else.
The finding is narrower and harder to escape: **the ranking claim is
conditional on a term whose correct value is unknown and whose measured spread
is 1308x.** A defence that holds only if a 1308x-spread quantity happens to
be exactly 1.0 for every compound is not a defence — it is an assumption that
has now been priced.

---

## 1. Numbers I was given, and what they actually are

The brief supplied one figure that was already wrong, and asked me to check
rather than inherit. Every measured-set figure below is recomputed from
`data_cache/brain_slice/vu_brain.csv` by this script, not transcribed.

| quantity | supplied | **recomputed** | verdict |
|---|---|---|---|
| measured spread, scored set | 631x (withdrawn), then 1308x | **1307.5x** (n = 33, 3.2–4184 mL/g) | the corrected figure is right; **631x is not** |
| spread across the 10 non-PET named drugs | 710x | **710.5x** (3.73–2650) | confirmed |
| `GMFE_floor` | 4.77x | **4.774x** at constant **39** | confirmed |
| `GMFE(engine)` / ratio | 57.50x / 12.05 | **57.50x** / **12.045** | confirmed |
| `s_eff` at the shipping state | 1.0 for every compound | **1.000000** exactly | confirmed |

**Two supplied figures I could not reproduce, and the reason.** The brief
quotes an AUC attenuation factor spanning **0.71–0.97** today and **0.042–0.30**
at a uniform R = 39, and cites `PDEConfig.total_time = 43200 s`. Two things:

1. **No shipping panel uses 43200 s.** `PDEConfig.total_time` *defaults* to
   43200 s, but `experiments/real_drug_panel.py` and
   `experiments/breast_drug_panel.py` both construct
   `PDEConfig(dt=10.0, total_time=3600.0*6, ...)`, i.e. **21600 s**. Every arm
   here uses the panels' own window. A shorter window makes capacity matter
   MORE, not less, so this correction is against my own conclusion, not for it.
2. On that window the fraction of the steady-state AUC the shipping GBM panel
   actually accumulates is **0.9075–0.9965** (1.10x spread), and at a
   uniform R = 39 it is **0.1267–0.8654** (6.8x). The direction and the
   order of magnitude of the brief's reasoning survive; the digits do not.
   **The conclusion below rests on the recomputed numbers.**

---

## 2. The mechanism, and the proof that it is implemented exactly

With a capacity term `R = s_eff` the homogenised problem is

```
s_eff * dC/dt = div(D grad C) - div(v C) + S - k*C        k = clearance + binding
=>      dC/dt = div((D/R) grad C) - div((v/R) C) + S/R - (k/R)*C
```

so each run divides `diffusivity`, `binding`, `clearance` and `bbb_penetrance`
by R and nothing else. `bbb_penetrance` enters the engine in exactly one place
(`treatment.drugs.delivery_source` returns `vascular_density * bbb_penetrance *
plasma`), so scaling it is precisely `S -> S/R`. **`ic50` and `hill` are never
touched**, so the potency axis is identical in every arm and any reordering is
a delivery-side result.

**(a) The retarded solve reproduces the closed form.** A 9x9 homogeneous domain,
`R dC/dt = S - kC`, against `AUC = (S/k)[T - (R/k)(1 - exp(-kT/R))]`:

| R | numeric AUC | analytic AUC | rel. error |
|---:|---:|---:|---:|
| 1 | 2.20376e+08 | 2.20381e+08 | 2.4e-05 |
| 39 | 1.14692e+08 | 1.14695e+08 | 2.2e-05 |
| 777 | 9.07951e+06 | 9.07953e+06 | 1.9e-06 |
| 4184 | 1.72442e+06 | 1.72442e+06 | 3.6e-07 |

**(b) The STEADY STATE is capacity-invariant, measured not assumed.** Run far
past equilibrium, `C_ss = S/k` must not move with R:

| R | final mean C | analytic `S/k` |
|---:|---:|---:|
| 1 | 20.6912331591 | 20.6912331591 |
| 39 | 20.6912331591 | 20.6912331591 |
| 4184 | 20.6911060258 | 20.6912331591 |

So capacity **cannot** change a converged spatial pattern, the penetration
length `L = sqrt(D/k)`, or any steady-state number. It acts on the ranking
**only** through the finite exposure window — and the attenuation depends on
`k*T/R`, i.e. on **k as well as R**. That is the whole reason a uniform R can
reorder a panel: the panel's drugs differ in `k` by
**27x** (0.000501 – 0.0134 1/s).

**(c) The kill mapping is the frozen one, not a re-implementation.**
`cells.pkpd.apply_drug_kill`'s own mean per-tumour-cell kill probability for
`temozolomide` is **0.1112551495**; this script's `E[kill]` is
**0.1112551495**; absolute difference **0.0e+00**.

**(d) Convection is OFF in every arm, and that is a limitation, not a choice.**
The velocity field is built inside `sim.engine.simulate_response` and cannot be
divided by R from outside, so an exactly-retarded advection term is
unavailable. `results/ranking_robustness.md` runs its Monte Carlo
convection-off for the same class of reason; sec. 7 re-measures the
convection on/off difference here rather than inheriting it.

---

## 3. UNIFORM capacity — the weakest possible version

Every compound gets the SAME R. Real Xenium GBM section, 40,887 measured cells (23,812 tumour), grid (339, 541), dx = 15 um, dt = 20 s, window 21600 s, sensitivity mode `uniform`. One reconstruction, 220 transport solves.

### 3.1 `E[kill]` across the capacity grid (10-drug GBM panel)

| drug | R=1 | R=3.2 | R=8 | R=20 | R=39 | R=120 | R=410 | R=1420 | R=4184 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **vincristine** | 0.5292 | 0.5291 | 0.5288 | 0.5280 | 0.5266 | 0.5177 | 0.4692 | 0.3134 | 0.1534 |
| **cisplatin** | 0.5260 | 0.5244 | 0.5198 | 0.5047 | 0.4770 | 0.3661 | 0.1801 | 0.0598 | 0.0187 |
| **carmustine** | 0.5218 | 0.5202 | 0.5158 | 0.5017 | 0.4765 | 0.3779 | 0.1978 | 0.0667 | 0.0214 |
| **doxorubicin** | 0.4781 | 0.4766 | 0.4729 | 0.4619 | 0.4404 | 0.3469 | 0.1796 | 0.0623 | 0.0199 |
| **lomustine** | 0.4753 | 0.4728 | 0.4666 | 0.4476 | 0.4138 | 0.2972 | 0.1345 | 0.0411 | 0.0127 |
| **etoposide** | 0.3441 | 0.3397 | 0.3294 | 0.3005 | 0.2555 | 0.1442 | 0.0515 | 0.0139 | 0.0038 |
| **irinotecan** | 0.2933 | 0.2915 | 0.2875 | 0.2769 | 0.2588 | 0.1883 | 0.0847 | 0.0259 | 0.0079 |
| **procarbazine** | 0.1569 | 0.1525 | 0.1426 | 0.1178 | 0.0877 | 0.0381 | 0.0112 | 0.0027 | 0.0007 |
| **methotrexate** | 0.1295 | 0.1266 | 0.1200 | 0.1035 | 0.0810 | 0.0381 | 0.0117 | 0.0029 | 0.0007 |
| **temozolomide** | 0.1113 | 0.0883 | 0.0552 | 0.0261 | 0.0134 | 0.0039 | 0.0010 | 0.0002 | 0.0001 |

### 3.2 The ordering, drug 1 first

| R | ordering (best predicted kill first) | = shipping? |
|---:|---|:--:|
| 1 | vincristine > cisplatin > carmustine > doxorubicin > lomustine > etoposide > irinotecan > procarbazine > methotrexate > temozolomide | **yes** |
| 39 | vincristine > cisplatin > carmustine > doxorubicin > lomustine > irinotecan > etoposide > procarbazine > methotrexate > temozolomide | **NO** |
| 55 | vincristine > carmustine > cisplatin > doxorubicin > lomustine > irinotecan > etoposide > procarbazine > methotrexate > temozolomide | **NO** |
| 180 | vincristine > carmustine > cisplatin > doxorubicin > lomustine > irinotecan > etoposide > methotrexate > procarbazine > temozolomide | **NO** |
| 620 | vincristine > carmustine > doxorubicin > cisplatin > lomustine > irinotecan > etoposide > methotrexate > procarbazine > temozolomide | **NO** |

**5 distinct orderings** appear across the 22-point capacity grid.

### 3.3 The five-drug sub-panel `results/ranking_robustness.md` scores

| R | ordering | = published? |
|---:|---|:--:|
| 1 | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | **yes** |
| 620 | doxorubicin > cisplatin > lomustine > methotrexate > temozolomide | **NO** |

---

## 4. WHICH pairs swap, and at what R

Crossing R located by bisection on a monotone (PCHIP) interpolant of
`E[kill]` against log R, built on the 22 solved points; interpolation error
is measured in sec. 7.

### 4.1 10-drug GBM panel — every pair that swaps inside the measured range

| pair | shipping leader | **swaps at R =** | inside measured range [3.2, 4184]? |
|---|---|---:|:--:|
| etoposide vs irinotecan | etoposide | **36.25** | **yes** |
| carmustine vs cisplatin | cisplatin | **42.00** | **yes** |
| procarbazine vs methotrexate | procarbazine | **122.66** | **yes** |
| cisplatin vs doxorubicin | cisplatin | **434.22** | **yes** |

**4 of 45 pairs swap** somewhere on the grid; **4** of them swap
at a capacity inside the measured brain-slice range.

### 4.2 Five-drug sub-panel

| pair | shipping leader | **swaps at R =** | inside measured range? |
|---|---|---:|:--:|
| cisplatin vs doxorubicin | cisplatin | **434.22** | **yes** |

---

## 5. PER-COMPOUND capacity

### 5.1 Monte Carlo over the measured distribution

40,000 draws. Each compound is given an INDEPENDENT capacity and the
panel is re-ranked inside every draw — the same re-rank-inside-every-draw
design `results/ranking_robustness.md` sec. 1 uses, and the same
decision-grade bar, quoted verbatim from that file:

> *the recommended drug holds rank 1 with probability whose 90% credible
> interval lies entirely above **0.90**, and the runner-up's lies entirely
> below **0.10**.*

**Declared deviation.** That report uses a *cluster* bootstrap because its
design is nested (256 IC50 draws inside 40 tissue draws). Here the tissue is
fixed and every capacity draw is independent, so there is no cluster to
resample and the interval is a plain 4,000-replicate percentile bootstrap over
draws (computed through the exact `Binomial(n, p)/n` identity for the
bootstrap of an indicator mean, not by materialising the resample matrix).
The bar is unchanged. **A second, larger deviation is stated rather than
hidden: independent per-compound draws almost certainly OVERSTATE the
instability**, because real capacity is a smooth function of lipophilicity and
neighbouring compounds do not draw independently. Sec. 5.2 gives the
correlated, mechanistic counterpart, **and it is substantially more benign than
this arm** — which is stated here, before the numbers, because it cuts against
the conclusion this section reaches.

**P(rank 1), five-drug sub-panel.** Shipping order: cisplatin > doxorubicin > lomustine > methotrexate > temozolomide.

| drug | A. per-compound, resampled from the 33 measured values | B. per-compound, log-uniform on [3.2, 4184] | C. UNIFORM (one draw shared by all compounds) |
|---|:--:|:--:|:--:|
| cisplatin | 0.584 [0.580, 0.588] | 0.471 [0.467, 0.475] | 0.818 [0.815, 0.822] |
| doxorubicin | 0.247 [0.243, 0.250] | 0.298 [0.294, 0.301] | 0.182 [0.178, 0.185] |
| lomustine | 0.167 [0.164, 0.170] | 0.225 [0.221, 0.228] | 0.000 [0.000, 0.000] |
| methotrexate | 0.002 [0.002, 0.003] | 0.006 [0.006, 0.007] | 0.000 [0.000, 0.000] |
| temozolomide | 0.000 [0.000, 0.000] | 0.001 [0.000, 0.001] | 0.000 [0.000, 0.000] |
| **P(full shipping order holds)** | 0.200 [0.197, 0.204] | 0.135 [0.132, 0.138] | 0.818 [0.815, 0.822] |
| **verdict at the 0.90 / 0.10 bar** | **not decision-grade** | **not decision-grade** | **not decision-grade** |

**P(rank 1), ten-drug panel.** Shipping order: vincristine > cisplatin > carmustine > doxorubicin > lomustine > etoposide > irinotecan > procarbazine > methotrexate > temozolomide.

| drug | A. per-compound, resampled from the 33 measured values | B. per-compound, log-uniform on [3.2, 4184] | C. UNIFORM (one draw shared by all compounds) |
|---|:--:|:--:|:--:|
| vincristine | 0.730 [0.726, 0.734] | 0.614 [0.610, 0.618] | 1.000 [1.000, 1.000] |
| cisplatin | 0.132 [0.130, 0.135] | 0.162 [0.159, 0.165] | 0.000 [0.000, 0.000] |
| carmustine | 0.106 [0.104, 0.109] | 0.139 [0.136, 0.142] | 0.000 [0.000, 0.000] |
| doxorubicin | 0.019 [0.018, 0.020] | 0.048 [0.046, 0.050] | 0.000 [0.000, 0.000] |
| lomustine | 0.012 [0.012, 0.013] | 0.034 [0.032, 0.035] | 0.000 [0.000, 0.000] |
| etoposide | 0.000 [0.000, 0.000] | 0.002 [0.002, 0.003] | 0.000 [0.000, 0.000] |
| irinotecan | 0.000 [0.000, 0.000] | 0.001 [0.001, 0.001] | 0.000 [0.000, 0.000] |
| procarbazine | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| methotrexate | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| temozolomide | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] |
| **P(full shipping order holds)** | 0.005 [0.004, 0.005] | 0.001 [0.000, 0.001] | 0.483 [0.479, 0.487] |
| **verdict at the 0.90 / 0.10 bar** | **not decision-grade** | **not decision-grade** | **DECISION-GRADE** |

Read column **C** first. It is the *uniform* case — one capacity value shared
by every compound, drawn from the measured distribution. That is the exact
case the "constant error cancels" defence covers, and the ordering does not
survive it.

### 5.2 Named, mechanistic per-compound assignments

Rather than only sampling, four specific capacity assignments are evaluated by
name. The third and fourth use the repo's OWN structural rule
(`two_route.binding_capacity_from_logp`, the Kilford/Hallifax-Houston logP
QSAR) — i.e. the term `results/intracellular_retardation.md` built and
refuted. Read sec. 6.1 before reading these as a recommendation.

**Ten-drug GBM panel.**

| capacity assignment | resulting order | = shipping? | top-1 |
|---|---|:--:|---|
| shipping (s_eff = 1 for every compound) | vincristine > cisplatin > carmustine > doxorubicin > lomustine > etoposide > irinotecan > procarbazine > methotrexate > temozolomide | **yes** | vincristine |
| uniform R = 39 (the GMFE-minimising constant, B-S1) | vincristine > cisplatin > carmustine > doxorubicin > lomustine > irinotecan > etoposide > procarbazine > methotrexate > temozolomide | **NO** | vincristine |
| per-compound: the repo's own structural K_b rule (two_route.binding_capacity_from_logp) | vincristine > cisplatin > carmustine > doxorubicin > lomustine > etoposide > procarbazine > irinotecan > methotrexate > temozolomide | **NO** | vincristine |
| per-compound: the SAME rule, access-limited to the 6 h window (two_route.equilibration_time_s) | vincristine > cisplatin > carmustine > doxorubicin > lomustine > etoposide > procarbazine > irinotecan > methotrexate > temozolomide | **NO** | vincristine |

**Five-drug sub-panel.**

| capacity assignment | resulting order | = shipping? | top-1 |
|---|---|:--:|---|
| shipping (s_eff = 1 for every compound) | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | **yes** | cisplatin |
| uniform R = 39 (the GMFE-minimising constant, B-S1) | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | **yes** | cisplatin |
| per-compound: the repo's own structural K_b rule (two_route.binding_capacity_from_logp) | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | **yes** | cisplatin |
| per-compound: the SAME rule, access-limited to the 6 h window (two_route.equilibration_time_s) | cisplatin > doxorubicin > lomustine > methotrexate > temozolomide | **yes** | cisplatin |

The per-compound `s_eff` those two rules produce, beside the measured range:

| drug | logP | `K_b` (rule) | `s_eff` (rule) | membrane equilibration tau (s) | tau vs the 6 h window | `s_eff` ACCESS-LIMITED |
|---|---:|---:|---:|---:|:--:|---:|
| temozolomide | -2.08 | 13.9 | 11.96 | 687 | tau < T | 11.96 |
| lomustine | 2.25 | 30.7 | 25.2 | 0.0728 | tau < T | 25.2 |
| carmustine | 1.16 | 14 | 12.02 | 0.415 | tau < T | 12.02 |
| procarbazine | 1.05 | 13.2 | 11.41 | 244 | tau < T | 11.41 |
| cisplatin | 0.20 | 9.7 | 8.654 | 66.4 | tau < T | 8.654 |
| methotrexate | 0.27 | 9.86 | 8.784 | 6.87e+05 | **tau >> T** | 1.241 |
| etoposide | 1.34 | 15.5 | 13.22 | 202 | tau < T | 13.22 |
| vincristine | 3.52 | 125 | 99.83 | 11.4 | tau < T | 99.83 |
| irinotecan | 4.09 | 282 | 223.4 | 0.0128 | tau < T | 223.4 |
| doxorubicin | 0.00 | 9.35 | 8.382 | 2.81e+07 | **tau >> T** | 1.006 |

For reference the measured rat brain-slice range is **3.2 – 4184**,
the best constant is **39**, and the engine ships **1.000**.

#### 5.3 The result in 5.2 cuts AGAINST sec. 5.1, and it is reported first

On the GBM panels the mechanistic per-compound assignments are **far more
benign than the sampled ones**, and pretending otherwise would be exactly the
shaped conclusion this exercise exists to prevent:

* On the **five-drug GBM sub-panel** every named assignment — uniform R = 39,
  the structural rule, and the access-limited structural rule — reproduces
  the **published order exactly**. Not just the top-1: the whole order.
* On the **ten-drug GBM panel** the top-1 (`vincristine`) is unchanged
  under every named assignment, and the structural rule moves exactly one
  adjacent pair.

**Why**, and it is not luck: the structural rule makes capacity an increasing
function of logP, and on the GBM panel the drugs with the largest capacity
(`irinotecan` 223, `vincristine` 100) are also the drugs with the fastest
sink `k`, which is the term that makes a compound INSENSITIVE to capacity. The
two per-compound effects partly cancel. Independent draws from the measured
distribution destroy that correlation by construction, which is why sec. 5.1
looks so much worse.

**And the cancellation does NOT hold on the breast panel — which is why
sec. 5.2 must not be read as a general reassurance.** There the structural
rule gives `tamoxifen` `s_eff` = 7218 and `vinorelbine` 652 against
`doxorubicin`'s 8.4, the correlation with `k` runs the other way, and
**the top-1 pick moves from `vinorelbine` to
`doxorubicin`**. Same rule, same engine, opposite conclusion.

**Neither arm is a measurement, and they disagree — with each other AND across
tissues.** The sampled arm assumes nothing about structure and says no ordering
is decision-grade anywhere. The mechanistic arm assumes a QSAR this repo has
already REFUTED on other data (sec. 6.1) and says the GBM ordering barely
moves while the breast TOP-1 changes. **The honest reading is that the
ordering's status depends on an unmeasured modelling choice, and that the
answer is not even the same sign on the two shipping tissues.** That is a
different failure from "the ranking is wrong": it is "the ranking is not
determined by what has been measured".

*Disclosure: one or more named assignments produced an `s_eff` outside the
swept grid [1, 4184] and were CLIPPED to it before being read off the
interpolant — `tamoxifen`'s structural-rule `s_eff` of 7218 exceeds the largest measured `Vu,brain`. Clipping
moves that compound UP the ordering relative to an unclipped value, i.e.
**against** the reordering reported here.*

---

## 6. Three things that would make this result mean less — all of them true

### 6.1 The direction is unknown, and the only RANKING evidence points the other way

**This is the first thing to read, and it is not a hedge.**
`results/intracellular_retardation.md` built a per-compound capacity term
`K_b` from the same structural rule used in sec. 5.2 and scored it on measured
Pruijn HT-29 multicellular-layer transport. It got **worse on every axis**, and
the axis it destroyed is the one this document is about:

| | matched-MW within-tie rank agreement |
|---|:--:|
| retardation **OFF** (shipping) | **77.4 %** (24 concordant / 7 discordant) |
| structural `K_b` rule **ON** | **51.6 %** (16 concordant / 15 discordant) |

(R^2 also fell 0.718 -> 0.606.) **The term this FAIL points at and the term that
experiment refuted are the same term.** Two measured datasets disagree: rat
brain-slice `Vu,brain` says the capacity slot is 39x too small and cannot be
fixed by a constant; Pruijn HT-29 says filling that slot with the only
available per-compound rule takes rank agreement to a coin flip.

**Two things about that refutation, both of which cut against reading it as a
clean contradiction of what is measured here.**

1. **It scores a different object, and on THAT object a uniform capacity
   provably cannot move a ranking.** That report's own control arm D says so:
   *"A compound-independent `K_b` is a rigid shift in log space, so it cannot
   change R^2 or Spearman — only bias"*, and the table confirms it at three
   magnitudes (R^2 pinned at 0.718 for uniform `K_b` = 1, 10, 100). The scored
   quantity there is a per-compound *transport coefficient*, which a uniform R
   divides identically. The quantity scored HERE is `E[kill]` from a
   **finite-window AUC through a Hill curve**, where R enters as `k*T/R` — a
   nonlinear function of each compound's own `k`. **That is the entire reason
   a uniform capacity is rank-preserving there and rank-breaking here**, and it
   is not a disagreement between the two results.
2. **Its frozen-state guard would now trip.** `experiments/intracellular_
   retardation.py:494` asserts `calibrate_alpha() == 148.60815765445545`; live
   it returns `22.223477711172563`, so that stage would refuse to run today.
   Checked rather than assumed: `alpha` is a state guard there and does not
   enter the scored quantity (`two_route` diffusivity / `retardation_factor`),
   so **the 77.4 % -> 51.6 % figure itself is unaffected**. Flagged because a
   refutation recorded in a superseded engine state should be labelled as one.

**I am not resolving the disagreement by preferring the dataset I measured.**
The correct reading is the weaker and more uncomfortable one:

> The ranking's stability is **conditional on a term whose value is unknown**.
> Setting it to 1.0 is not a defensible default that happens to be wrong in
> the absolutes — it is a load-bearing choice that the ordering depends on.
> Both datasets agree on that much, and neither licenses shipping a fix.

### 6.2 Possible double-counting — and it is worse than the brief suggests

The brief raises GDSC potency: IC50 is measured on cells in medium, where the
same intracellular partitioning occurs, so some of `Vu,brain`'s capacity may
already sit inside the potency reference. That is real, and it cuts a specific
way: a GDSC IC50 is referenced to the **medium** concentration, and the engine
compares it to the **free interstitial** concentration, which is the
commensurable pair — so at EQUILIBRIUM the comparison is already consistent
and capacity should not be applied twice. The engine's window is not at
equilibrium, which is why the effect above is non-zero; but the size of the
effect is inflated by however much of the partitioning the assay already
absorbed, and that fraction is not known.

**There is a second, larger double-count inside the engine itself, and this one
is checkable.** `validation.calibration.uptake_proxy` — the quantity `alpha`
multiplies to produce every compiled drug's `binding` rate, which is the
dominant part of the sink `k` — is documented in the repo's own words as:

> *"Molecular proxy for RAPID CELLULAR UPTAKE + HIGH-CAPACITY intracellular /
> DNA binding ... It is intentionally a **capacity ordering**, not an absolute
> rate."*

So **the engine already carries a per-compound intracellular-capacity ordering**
— wired into an irreversible first-order SINK rather than into the reversible
capacity slot `s_eff`. That is this repo's recurring "wrong term doing the
right work" pattern, and it means the sensitivity measured above is partly a
re-application of physics the model already has.

Measured rather than asserted, across the 10-drug GBM panel: Spearman rho
between `uptake_proxy` and the structural `K_b` is **+0.333** (n = 10, p = 0.35); Pearson r on logs **+0.630** (p = 0.051). **The overlap is real but
weak and not significant at n = 10** — so this neither dismisses the
sensitivity nor lets it be taken at face value. It is a reason the effect
size above should be read as an upper bound.

### 6.3 Access, not just capacity

`two_route.retardation_factor`'s own docstring:

> *"A cell behind a shut membrane cannot store solute any faster than the
> membrane admits it, so for large `membrane_resistance_ratio` the full
> capacity is not accessed within a finite experiment and this factor
> OVER-states the retardation."*

`two_route.equilibration_time_s` returns that timescale, and it is decisive
here. Against the panels' 6 h window (21,600 s):

* **2 of 10** GBM-panel drugs have tau > T: `methotrexate` (tau = 6.87e+05 s), `doxorubicin` (tau = 2.81e+07 s).
* For those compounds the accessible capacity over the window collapses toward the shipping value of 1.0 — see the last column of sec. 5.2.

**Access-limiting is itself strongly per-compound**, so it does not act as a
softening of the correction: it nearly cancels it for the slowest permeators
while leaving it untouched for the fastest, WIDENING the per-compound spread
of `s_eff` on the GBM panel from **27x** to **222x**.

**And yet, measured rather than predicted: the access-limited arm produces
the IDENTICAL ordering to the un-limited structural rule on all 3 panels.** A
wider spread in `s_eff` did not translate into a different order, because
the compounds it moves (on the GBM panel, `methotrexate` and `doxorubicin`)
were not near a crossing. This is reported because the obvious expectation
— bigger spread,
more reordering — is simply not what happened, and the correct response to
that is to print it, not to explain it away.

---

## 7. Is any of this a numerical artifact?

**Interpolation.** R = 39 was held out of the
interpolant and predicted:

| drug | predicted `E[kill]` | computed `E[kill]` | abs. error |
|---|---:|---:|---:|
| temozolomide | 0.013473 | 0.013384 | 8.9e-05 |
| lomustine | 0.413501 | 0.413824 | 0.00032 |
| carmustine | 0.476240 | 0.476531 | 0.00029 |
| procarbazine | 0.087736 | 0.087697 | 3.9e-05 |
| cisplatin | 0.476652 | 0.477005 | 0.00035 |
| methotrexate | 0.081026 | 0.081044 | 1.8e-05 |
| etoposide | 0.255258 | 0.255463 | 0.00021 |
| vincristine | 0.526529 | 0.526557 | 2.8e-05 |
| irinotecan | 0.258583 | 0.258842 | 0.00026 |
| doxorubicin | 0.440054 | 0.440365 | 0.00031 |

Worst absolute interpolation error **0.00035**.

**Production grid.** The sweep runs at dx = 15 um / dt = 20 s
(the `results/ranking_robustness.md` Monte-Carlo configuration); the panel's
published configuration is dx = 12 um / dt = 10 s. Re-run there:

| R (dx = 12, dt = 10, convection off) | ordering | = the dx = 15 ordering? |
|---:|---|:--:|
| 1 | vincristine > cisplatin > carmustine > doxorubicin > lomustine > etoposide > irinotecan > procarbazine > methotrexate > temozolomide | **yes** |
| 13 | vincristine > cisplatin > carmustine > doxorubicin > lomustine > etoposide > irinotecan > procarbazine > methotrexate > temozolomide | **yes** |
| 39 | vincristine > cisplatin > carmustine > doxorubicin > lomustine > irinotecan > etoposide > procarbazine > methotrexate > temozolomide | **yes** |
| 120 | vincristine > carmustine > cisplatin > doxorubicin > lomustine > irinotecan > etoposide > procarbazine > methotrexate > temozolomide | **yes** |

**Convection.** At dx = 12 the on/off difference in `E[kill]` is at most
**0.00041** at R = 1 and **0.0045** at R = 39, and the ordering is
identical with it on and off at R = 1.

---

## 8. The breast panel — the same test, plus one directly measured compound

Real Xenium breast section (Janesick rep1), 167,780 cells,
grid (366, 503), same configuration.

| R | ordering | = shipping? |
|---:|---|:--:|
| 1 | vinorelbine > cisplatin > paclitaxel > doxorubicin > docetaxel > epirubicin > 5-fluorouracil > cyclophosphamide > palbociclib > tamoxifen | **yes** |
| 8 | vinorelbine > paclitaxel > doxorubicin > cisplatin > docetaxel > epirubicin > 5-fluorouracil > cyclophosphamide > palbociclib > tamoxifen | **NO** |
| 13 | vinorelbine > paclitaxel > doxorubicin > docetaxel > cisplatin > epirubicin > cyclophosphamide > 5-fluorouracil > palbociclib > tamoxifen | **NO** |
| 20 | vinorelbine > paclitaxel > docetaxel > doxorubicin > cisplatin > epirubicin > cyclophosphamide > 5-fluorouracil > palbociclib > tamoxifen | **NO** |
| 30 | vinorelbine > paclitaxel > docetaxel > doxorubicin > epirubicin > cisplatin > cyclophosphamide > palbociclib > 5-fluorouracil > tamoxifen | **NO** |
| 55 | vinorelbine > paclitaxel > docetaxel > doxorubicin > epirubicin > cisplatin > palbociclib > cyclophosphamide > tamoxifen > 5-fluorouracil | **NO** |
| 120 | vinorelbine > paclitaxel > docetaxel > doxorubicin > epirubicin > cisplatin > palbociclib > tamoxifen > cyclophosphamide > 5-fluorouracil | **NO** |
| 270 | vinorelbine > paclitaxel > docetaxel > doxorubicin > epirubicin > tamoxifen > cisplatin > palbociclib > cyclophosphamide > 5-fluorouracil | **NO** |
| 620 | vinorelbine > paclitaxel > docetaxel > doxorubicin > epirubicin > tamoxifen > palbociclib > cisplatin > cyclophosphamide > 5-fluorouracil | **NO** |

* full order stable up to **R = 5.28**; top-1 stable up to **R = &gt; 4184 (no change anywhere on the swept grid)**.
* **13 of 45** pairs swap on the grid.

**One compound on this panel has a directly measured value.** `Docetaxel`
appears in the scored brain-slice set at **Vu,brain = 777 mL/g** (Loryan
2013, PMID 23336814, rat, ion class Base). That is a rat-brain measurement
applied to a human breast tumour and it constrains nothing about this
tissue — but it is the only place on either shipping panel where the
capacity of a panel drug has been measured at all, and it is **777x** the
value the engine assumes for it.

At a uniform R = 777 the breast ordering is:

> vinorelbine > paclitaxel > docetaxel > doxorubicin > epirubicin > tamoxifen > palbociclib > cisplatin > cyclophosphamide > 5-fluorouracil

against the shipping order

> vinorelbine > cisplatin > paclitaxel > doxorubicin > docetaxel > epirubicin > 5-fluorouracil > cyclophosphamide > palbociclib > tamoxifen

**Named per-compound arms, breast panel:**

| capacity assignment | resulting order | = shipping? | top-1 |
|---|---|:--:|---|
| shipping (s_eff = 1 for every compound) | vinorelbine > cisplatin > paclitaxel > doxorubicin > docetaxel > epirubicin > 5-fluorouracil > cyclophosphamide > palbociclib > tamoxifen | **yes** | vinorelbine |
| uniform R = 39 (the GMFE-minimising constant, B-S1) | vinorelbine > paclitaxel > docetaxel > doxorubicin > epirubicin > cisplatin > cyclophosphamide > palbociclib > 5-fluorouracil > tamoxifen | **NO** | vinorelbine |
| per-compound: the repo's own structural K_b rule (two_route.binding_capacity_from_logp) | doxorubicin > cisplatin > docetaxel > paclitaxel > epirubicin > vinorelbine > 5-fluorouracil > cyclophosphamide > palbociclib > tamoxifen | **NO** | doxorubicin |
| per-compound: the SAME rule, access-limited to the 6 h window (two_route.equilibration_time_s) | doxorubicin > cisplatin > docetaxel > paclitaxel > epirubicin > vinorelbine > 5-fluorouracil > cyclophosphamide > palbociclib > tamoxifen | **NO** | doxorubicin |

Monte Carlo, per-compound empirical: P(top-1 = `vinorelbine`) = **0.797**
[0.794, 0.801]; P(full shipping order) = **0.001**.

---

## 9. Is there a regime where the ordering IS stable?

Yes, and it is worth stating precisely, because it is the only defensible
version of the claim that survives.

| panel | full order stable up to R = | top-1 stable up to R = |
|---|---:|---:|
| GBM, 10 drugs | **36.3** | **&gt; 4184 (no change anywhere on the swept grid)** |
| GBM, 5-drug sub-panel | **435** | **435** |
| breast, 10 drugs | **5.28** | **&gt; 4184 (no change anywhere on the swept grid)** |

*"&gt; 4184" means the ordering never changed anywhere on the swept capacity grid.*

> **The condition the claim could honestly carry:**
>
> *The predicted ordering is stable against a uniform tissue-capacity term up
> to R ~ 5.28 for the full order and R ~ 435 for the top-1 pick. The engine
> assumes R = 1. Measured rat brain-slice `Vu,brain` spans 3.2 – 4184 with a
> best constant of 39.*
>
> **The stable regime does not contain the measured range**, and for the full order it does not contain even the best available constant (39).

There is a second stable regime and it is worth naming because it is a
genuine defence of the physics: **capacity cannot touch a converged field at
all** (sec. 2b). Any Codon output that is a steady-state or
equilibrium-referenced quantity — penetration length, interstitial:plasma
ratio, the shape of a converged gradient — is exactly capacity-invariant and
this entire failure mode does not reach it. The failure is specific to
**finite-window AUC**, which is what the kill model and therefore the ranking
is built on.

---

## 10. What may and may not be claimed from this document

**Claimable.**

* The engine's predicted FULL compound ordering is NOT robust to a
  tissue-capacity term over the range that term has been measured to span:
  on the 10-drug GBM panel it breaks at a UNIFORM R = 36.3, below the
  best available constant of 39 and inside [3.2, 4184].
* It is therefore not robust even to a UNIFORM capacity, so the "a constant
  error cancels out of a ranking" defence does not apply to this axis. The
  reason is mechanical and is stated in sec. 2: the AUC attenuation depends
  on `k*T/R`, and the panel's `k` spans 27x.
* Under a PER-COMPOUND capacity drawn from the measured distribution, no
  ordering claim on any of the 3 panels tested clears this repo's
  own decision-grade bar (sec. 5.1).
* The 10-drug top-1 pick (`vincristine`) IS stable against a uniform
  capacity across the whole measured range, and the five-drug top-1
  (`cisplatin`) to R = 435. That is a real result and is not
  softened here.
* On the **breast 10-drug** panel the TOP-1 PICK CHANGES from `vinorelbine` to
  `doxorubicin` under this repo's own structural `K_b` rule and under its
  access-limited variant — a named assignment, not a sampling tail.
* Under this repo's own structural `K_b` rule — the only mechanistic
  per-compound capacity assignment available — the FIVE-DRUG published
  order is reproduced exactly and the ten-drug top-1 is unchanged
  (sec. 5.2/5.3). This is reported as claimable because it is true, not
  because it is convenient; sec. 6.1 explains why it cannot be leaned on.
* Steady-state and equilibrium-referenced outputs are exactly
  capacity-invariant (sec. 2b, measured to 1e-9).

**NOT claimable.**

* That a capacity term should be added to the engine. `results/intracellular_
  retardation.md` refuted the only available per-compound rule on measured
  Pruijn HT-29 rank agreement (77.4 % -> 51.6 %). Nothing here reverses that.
* That any ordering in sec. 3-5 is more correct than the shipping one. Every
  alternative ordering here is produced by a capacity value that is a
  hypothesis, not a measurement on these compounds in this tissue.
* Any absolute kill number in this document. All arms run sensitivity mode
  `uniform`, convection off, and one section with no outcome label.
* Anything about `cell_partition` (K), which `results/transport_anchor_prereg.md`
  sec. 3.2 / 3.8 leaves unconstrained by any measurement. This study does not
  constrain it either.
* That the effect size is right. Sec. 6.2 gives a concrete reason it is an
  upper bound, and sec. 5.3 gives a second: independent per-compound draws
  destroy a structure-driven correlation that the real quantity has.
* That the sampled arm and the mechanistic arm can be averaged, or that one
  supersedes the other. They disagree, neither is measured on these
  compounds in this tissue, and the disagreement IS the result.

---

## 11. Reproduction

```
./.venv/bin/python experiments/capacity_ranking_robustness.py --stage algebra
./.venv/bin/python experiments/capacity_ranking_robustness.py --stage sweep
./.venv/bin/python experiments/capacity_ranking_robustness.py --stage breast
./.venv/bin/python experiments/capacity_ranking_robustness.py --stage analyse
./.venv/bin/python experiments/capacity_ranking_robustness.py --stage gridcheck
./.venv/bin/python experiments/capacity_ranking_robustness.py --stage report
```

Monte-Carlo seed **20260821** (the brain-slice anchor's own), bootstrap seed
**20260822**, 40,000 draws, 4,000 bootstrap replicates.

**Nothing in `tmesim/`, `data_cache/`, `docs/`, `app/`, `landing-3d/`,
`results/transport_anchor_prereg.md`, `results/brain_slice_anchor.md`,
`results/CLAIM_INVENTORY.md`, `results/kill_calibration*.md` or `README.md`
was modified by this run.** Two files were added and nothing else:
`experiments/capacity_ranking_robustness.py` and this document, plus the
frozen record `results/sweep/capacity_ranking_robustness.json`. The `--stage
sweep` / `--stage breast` runs additionally write a per-drug resume
checkpoint beside that record while they are running and it is not an output.

**No assertion was deleted or weakened.** `./.venv/bin/python -m pytest -q`
-> **1171 passed, 0 failed** (379.20 s), the unchanged baseline. One test did
fail on an earlier draft of the script and the fix was to change the script,
not the test: `tests/test_propagation.py::test_no_undeclared_paths[
transvascular_exchange]` correctly caught the sec. 2a verification calling
`physics.diffusion.simulate_drug_transport` directly, which would have made
`PDEConfig.transvascular_exchange` a silent no-op on this file. The
verification was rerouted through `sim.engine.simulate_response`, the shipped
entry point, and the numbers are unchanged (identical relative errors).

**Rider disclosure — the working tree was not clean.** Sibling agents were live
in this repository during the run and the tree carries uncommitted edits,
including one to `tmesim/molecular/two_route.py`, which IS a module this study
reads. Inspected rather than waved past: that edit is a docstring-only citation
correction (the Abraham HBD-coefficient provenance note, see
`results/abraham_citation_erratum.md`); it changes no constant and no code
path. Every constant this study leans on
(`calibrate_alpha()`, `two_route.DEFAULT_CELL_VOLUME_FRACTION`, the four
`PDEConfig` defaults, `s_eff` at the shipping state) is asserted against LIVE
code at the start of every stage and matched; the sweep, the grid check and
the algebra stage each re-assert them. Saying the tree was dirty is cheaper
than asserting it was not.

