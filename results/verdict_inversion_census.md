# How often does the headline verdict reverse? — 5 of 12, and temozolomide is one

> **CORRECTION (2026-08-29, no number changed).** This file measures TWO
> reversals per compound and the summary lines below used to quote them
> interchangeably: the **delivery/dose verdict** (the panel's own headline
> sentence, **5 of 12**) and the **four-way dominant limiter** printed under
> it (**7 of 12**). §1's table always carried both and no measured value here
> has moved — but the conflation propagated: `app/server/jobs.py` shipped one
> boolean built from §4's four names and told `irinotecan`, `palbociclib` and
> `ribociclib` that their delivery/dose verdict reverses, while citing this
> file, which measures it holding. The corrected sentences are marked below,
> and the tiers are now asserted against `verdict_inversion_census_data.json`
> by `tests/test_verdict_stability.py`.

**Status: RUN AND MEASURED.** Driver: `experiments/verdict_inversion_census.py`.
Data: `results/verdict_inversion_census_data.json`. Runtime **25,502 s** (7.1 h)
under heavy contention.

**Nothing shipped was changed.** `PDEConfig.transvascular_exchange`,
`transvascular.SPLIT_DEFAULT` and `drugs.PS_SCALE_S_INV` are asserted at their
shipped values; `calibrate_alpha()` asserted live at `22.223477711172563`. No
parameter was fitted and no default was moved. The `ON` arm is the exchange term
switched on at the **shipped** `PS`, with nothing scaled.

---

## 0. Verdict up front

`results/permeability_split.md` established that the per-cell delivery/dose
verdict **inverts** for daunorubicin under a transvascular exchange term the
model is known to be missing and ships in the OFF position. That was 4 compounds
on 1 section, and the panel was updated to say the verdict "can reverse". It
could not say how often.

**The delivery/dose verdict reverses for 5 of the 12 scored compounds —
41.7 %.** The coarser four-way dominant limiter moves for more: **7 of 12 —
58.3 %.** They are different claims, the table in §1 scores them in separate
columns, and the panel's headline sentence is the first one.

| measure | value |
|---|---:|
| compounds in the library | 13 |
| scored (excluding one degenerate zero-source) | 12 |
| **delivery↔dose verdict inverts** (the panel's own sentence) | **5 / 12 = 41.7 %** |
| four-way dominant limiter inverts (the label under it) | 7 / 12 = 58.3 % |
| worst per-cell relabelling | **57.2 %** of tumour cells |
| total cell labels changed across the panel | 143,899 |

**And temozolomide is one of them.** The flagship compound, on the flagship
section, in the configuration the landing page and every results file are
written around, flips from **dose** to **delivery** with 57.2 % of its cells
changing label. The single most-quoted run in this repository has a headline
verdict that is not robust to a term the model does not have.

---

## 1. The census

Flagship 10x Xenium GBM section, dx = 12 µm, repair 1, dt = 10 s, 6 h window,
convection ON, `plasma_peak` 50 µM. Exchange OFF vs ON at the shipped `PS`.

| compound | verdict OFF | verdict ON | limiter inverts | cells relabelled | mean AUC ON/OFF |
|---|---|---|---:|---:|---:|
| `bevacizumab` *(degenerate)* | delivery | delivery | — | 0.0% | — |
| `carmustine` | delivery | dose | **YES** | 42.6% | 0.062 |
| `cisplatin` | delivery | dose | — | 21.4% | 0.054 |
| `doxorubicin` | delivery | delivery | — | 0.6% | 0.837 |
| `etoposide` | dose | dose | — | 11.2% | 0.693 |
| `irinotecan` | dose | dose | **YES** | 22.5% | 0.633 |
| `lomustine` | delivery | dose | **YES** | 48.2% | 0.103 |
| `methotrexate` | dose | dose | — | 2.3% | 0.806 |
| `palbociclib` | dose | dose | **YES** | 39.8% | 0.272 |
| `procarbazine` | dose | delivery | **YES** | 57.1% | 0.097 |
| `ribociclib` | dose | dose | **YES** | 49.1% | 0.200 |
| `temozolomide` | dose | delivery | **YES** | 57.2% | 0.140 |
| `vincristine` | delivery | delivery | — | 0.1% | 0.927 |

`bevacizumab` is excluded from every denominator: `bbb_penetrance = 0` makes its
source identically zero, both arms are bit-identical, and its DELIVERY share is
exactly 1.0. That is an **input**, not a finding, and counting it would have
flattered the inversion rate by diluting it.

## 2. What the exchange term does, and why the verdict moves

The `mean AUC ON/OFF` column is the mechanism. The shipped one-way law has no
saturation — entry elasticity is exactly 1 — so switching on a wall resistance
cuts tissue exposure, by between **7 %** (`vincristine`, 0.927) and **95 %**
(`cisplatin`, 0.054). A compound whose exposure falls by an order of magnitude
crosses the IC50 threshold in many more cells, and cells move from
SURVIVED_UNDERDOSED to SURVIVED_UNREACHED.

> **CORRECTION (2026-08-29).** This paragraph previously read "every compound
> whose verdict inverts loses most of its exposure (ratios 0.062–0.633), and
> every compound that is stable keeps most of it (0.693–0.927)". Those two bands
> partition the panel only if `cisplatin` is dropped from it: its dominant
> limiter does NOT invert and it retains **0.054** — the lowest ratio in the
> table. The measured structure is three bands, not two, and it is stated below.
> No ratio changed.

The correlation is the readable part, and it separates the panel into **three
bands that do not overlap**:

| band | mean AUC ON/OFF | compounds |
|---|---:|---|
| delivery/dose verdict **reverses** | 0.054–0.140 | `carmustine`, `cisplatin`, `lomustine`, `procarbazine`, `temozolomide` |
| verdict holds, **four-way limiter inverts** | 0.200–0.633 | `irinotecan`, `palbociclib`, `ribociclib` |
| **both hold** | 0.693–0.927 | `doxorubicin`, `etoposide`, `methotrexate`, `vincristine` |

`doxorubicin` at 0.837 relabels 0.6 % of cells; `temozolomide` at 0.140 relabels
57.2 %. The verdict is stable exactly where the missing term would not have
mattered much — and the middle band is the one a two-state summary has nowhere
to put, which is how three compounds ended up labelled as reversing.

## 3. Checks

* A ok — the split strictly reduced exposure on all 12 live-source compounds; an inert flag would have tied them.
* B ok — reproduces §4.2's doxorubicin and daunorubicin delivery:dose splits to 0.000e+00.
* C ok — 143,899 label changes counted between the arms across the panel, each reconciling with its transition matrix.
* D ok — 1 zero-source compound(s) ['bevacizumab'] detected by value, bit-identical across the arms, delivery share exactly 1.0, and excluded from the inversion denominator.

## 4. What this does and does not license

**It does not say the ON arm is right.** `results/permeability_split.md` §4 is
explicit that the split buys **falsifiability, not accuracy**: making it honest
needs a per-compound wall permeability and a tissue surface-to-volume ratio this
repository does not have, and the implied `P` spans 0.0206–0.102 µm/s entirely
on an unmeasured vessel radius. Neither arm is the truth.

**It says the verdict is not determined by the data.** Two defensible treatments
of the same missing physics give opposite answers for more than half the panel,
so for those compounds the delivery/dose call is a statement about a modelling
choice rather than about the tissue.

**The operational consequence, and the panel should say it:** a delivery/dose
verdict is only trustworthy where the compound retains most of its exposure
under the split. On this section that is `vincristine`, `doxorubicin`,
`methotrexate` and `etoposide`. It is not trustworthy for temozolomide, which is
the compound the product demonstrates itself with.

That is a recommendation about **exposure retention**, and it is not the same
question as whether a given compound's verdict was measured reversing — the
distinction §2's three bands make. `irinotecan`, `palbociclib` and `ribociclib`
kept their delivery/dose verdict and lost their four-way limiter; reading this
paragraph as the reversal column is precisely the error the `CORRECTION` at the
top of this file records.
