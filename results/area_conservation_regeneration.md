# Regenerating the published results under `geometry='area'` — and what the sub-grid closure does once the field stops being binary

> ### One-line answer
>
> **The nine published verdicts and the nine published ratios survive the
> promotion.** Re-solved live at `dx = 12 µm` — both arms, one process, one
> reconstruction, one `alpha` — the `binary` control reproduces all nine
> published ratios to the precision they are published at (worst deviation
> **4.92e-05** against values quoted to four decimals), **zero of nine verdicts
> flip**, and the largest ratio move is **−2.622 %** (GBM cisplatin). The
> project's only falsification (breast docetaxel **FAIL**, 0.148430 → 0.151639)
> survives and so does the flagship (GBM temozolomide **PASS**, 0.262649 →
> 0.258059). This independently reproduces `vessel_source_conservation.md` §7's
> "largest ratio move 2.62 %" from a re-solve rather than a replay.
>
> **Nothing in `results/p99_convention.md`'s convention ranking moves either.**
> Scored on the same solves, all **thirteen** homogeneous normalisers move by at
> most **6.096 %** (`max`, GBM doxorubicin — the single-cell order statistic that
> report already names as the fragile one). The shipped `p99` moves **2.622 %**;
> the `tail90` convention that report recommends over it moves **2.287 %**, so
> `tail90` stays the better-behaved of the two under a source change as well as
> under a grid change.
>
> **AND THERE IS SOMETHING THAT MOVES, WHICH THE ORIGINAL BLAST RADIUS COULD NOT
> SEE.** The nine verdicts all run under the **shipped** oxygen constants. Under
> the **jointly calibrated** constants of `results/oxygen_joint_grid_test.md` the
> area painting is **not** uniformly an improvement: it removes **99.0 %** of the
> `F_state` grid span on the shipped configuration (**0.04609 → 0.00044**) and
> **52.9 %** on joint-lo (0.15746 → 0.07418), but **increases it by 52.7 % on
> joint-hi** (0.05418 → **0.08271**). The published `dx = 12 µm` values in that
> report move by **+0.0863** (joint-lo) and **+0.0533** (joint-hi) against
> **−0.0080** on the shipped configuration — an order of magnitude larger.
> **`results/oxygen_joint_grid_test.md` §1 would be stale under the promotion**
> and is not regenerated here (§5). **This is a bookkeeping consequence, not a
> verdict move**, and §6 says plainly what it does and does not license.
>
> **The sub-grid closure is no longer a no-op, and that is a hard finding.**
> `results/subgrid_closure.md` §5.0's control — *"on the shipped field the
> closure is an exact no-op"* — is a property of the **binary** painting, not of
> the closure. Measured directly: the binary field has **exactly 0** partial
> cells at **every** rung of the ladder, so `beta ≡ 1` by construction. The
> `geometry='area'` field at `dx = 12 µm` has **16,228** partial cells —
> **99.92 %** of all source cells — with `beta` to **144.0** against a validated
> maximum of **5.677**. On that field the closure **destroys** the very
> grid-independence the source fix buys: `shipped` `F_state` span
> **0.00044 → 0.14320**, a **325-fold** degradation, against a baseline that is
> **not** degenerate (§4.2). It helps **one of three** configurations, exactly as
> it did on `'subgrid'`. **`results/subgrid_closure.md` is corrected (§7); it is
> hand-written, not generated, and the correction is a hand edit for that
> reason.**
>
> **VERDICT: the promotion is safe on the evidence, and default-off is now
> load-bearing for the closure.** See §6, including the one condition.

**Status: a measurement. NOTHING SHIPPED WAS CHANGED.** No file in `tmesim/` was
touched by this work. `OxygenConfig.supply = 0.05`, `consumption = 0.02`,
`hypoxia_threshold = 0.2`, the band edges, `VESSEL_RADIUS_UM = 7.5`,
`LINK_DISTANCE_UM = 20.0`, `real_fields._DEFAULT_VASCULAR_MODE = 'discrete'`,
`discrete_vascular_field(geometry='binary')` and the production `dx = 12 µm` are
all exactly as they were. `calibrate_alpha()` returns
**148.60815765445545** — checked before and after (§8). **The default was NOT
flipped; this file is evidence, not a decision.**

* Generator: `experiments/area_conservation_regeneration.py`
  (`--stage alpha | closure | ratios`)
* Also modified: `experiments/subgrid_closure.py` — `_tissue_arm` gained a
  `painting='area'` branch so the closure could be measured on the promotion
  candidate's field. `'binary'` and `'subgrid'` behaviour is unchanged and the
  existing `subgrid_transfer_partial.json` checkpoint stays valid.
* Raw: `results/sweep/areagen_alpha.json`, `areagen_closure.json`
  (+ `areagen_closure_partial.json`), `areagen_ratios_dx12.json`
  (+ `areagen_ratios_dx12.jsonl`)
* Section: `data_cache/xenium_gbm` (40,887 cells / 23,812 tumour) and the breast
  sample used by `experiments/kill_calibration_breast.py`
* Predecessor: `results/vessel_source_conservation.md` §10, which named items 1
  and 2 below as preconditions for a promotion

---

> ## STATUS — EVERY NUMBER IN THIS FILE IS REPORTED IN TWO ARMS
>
> **This file follows the two-arm convention `results/mechanism_attribution.md`
> and `results/vessel_profile.md` established for the two delivery laws.** The
> source painting is an **OPEN DECISION**, so no single number below is
> presented as "the" answer:
>
> | arm | painting | painted vessel area | status as of this writing |
> |---|---|---|---|
> | **`binary`** | centre-in-disc YES/NO, deposits a whole `dx²` cell | **grid-dependent**: 713,700 → 748,656 µm² over `dx = 30…5 µm`, a **4.7312 %** span | **the SHIPPED default** (`discrete_vascular_field(geometry='binary')`) |
> | **`area`** | the same discs, same radius, same endothelial cells, measured as covered-area FRACTION | **742,713.0 µm² at every rung**, span **6.3e-14 %** | **opt-in, default-off** (`geometry='area'`, `results/vessel_source_conservation.md`) |
>
> **The binary numbers are NOT overwritten anywhere in this file.** Every table
> carries both. Where only one arm is shown it is labelled.
>
> **What to quote.** If you quote a number from this file you must name the
> painting it came from. "The `F_state` grid span is 0.00044" is not a quotable
> statement about this engine; "under the opt-in area-conserving painting the
> `F_state` grid span is 0.00044, against 0.04609 under the shipped binary one"
> is.
>
> **`'area'` is not `'subgrid'`.** `'subgrid'` paints one footprint per vessel
> *object* at a data-derived radius and deposits ~347k µm² against ~744k;
> `results/delivery_source_geometry.md` rejected it for the delivery source and
> it flipped three of seven scored rows. `'area'` keeps the shipped object, the
> shipped radius and the shipped per-endothelial-cell weighting and changes only
> the rasterisation. **The argument against one is not an argument against the
> other, and — §4 — a result measured on one does not transfer to the other.**

---

## 0. What was re-solved and what was re-scored

The brief asks for cached fields to be re-scored rather than re-solved wherever
the quantity permits it. Here is the ledger, because two of the four quantities
do **not** permit it:

| quantity | re-solved or re-scored? | why |
|---|---|---|
| the nine published ratios (§1) | **RE-SOLVED**, both arms | the source field enters the transport PDE. There is no cached exposure field to re-score. |
| the thirteen normaliser conventions (§2) | **RE-SCORED**, off §1's solves | `decompose` reads only the per-tumour-cell exposure vector `e_unit`. All thirteen conventions come from **one** set of solves, at no extra solve cost, by capturing `e_unit` at the point the shipped `decompose` is called. |
| the `F_state` / mean-O2 ladder (§3) | **RE-SOLVED** | the vascular field enters the steady oxygen solve. |
| the sub-grid closure (§4) | **RE-SOLVED** | `beta` enters the *reaction diagonal* of the oxygen operator, not the post-solve scoring, so a closed arm cannot be recovered from a raw arm's field. |

**The `binary` control was re-solved rather than replayed, and this is a
deliberate departure from `vessel_source_conservation.md` §7.** That report
replayed the published `binary` panels from
`delivery_source_geometry`'s checkpoint, on the correct reasoning that a replay
cannot drift. The cost is that the reproduction becomes a tautology. Here both
arms are solved in the same process, off the same reconstruction, with the same
`alpha`, so §1's agreement with `PUBLISHED_RATIO` is **a live check that the
pipeline still produces the published numbers** — and it does, to the precision
they are published at. It also produces `e_unit` for the `binary` arm, without
which §2 could not exist.

`experiments/delivery_source_geometry.py` was **not edited** (it is owned
elsewhere). Its `_source_field` is patched **in this process only**, and the
patch delegates every arm it does not own back to the original — the same
technique `experiments/vessel_source_conservation.py` already uses. Likewise
`kill_calibration.decompose` is wrapped by a shim that **records its inputs and
returns the shipped function's return value unchanged**; it computes nothing.

---

## 1. The nine published ratios, regenerated

`dx = 12 µm`, `alpha = 148.60815765445545`, both arms re-solved.
Raw: `results/sweep/areagen_ratios_dx12.json`.

| tissue / drug | published | **`binary`** (control) | **`area`** | move | verdict `binary` | verdict `area` |
|---|---:|---:|---:|---:|---|---|
| gbm / temozolomide | 0.2626 | **0.262649** | **0.258059** | **−1.748 %** | PASS | **PASS** |
| gbm / cisplatin | 0.2791 | **0.279118** | **0.271801** | **−2.622 %** | PASS | **PASS** |
| gbm / doxorubicin | 0.2077 | **0.207676** | **0.207578** | −0.047 % | PASS | **PASS** |
| breast / doxorubicin | 0.2301 | **0.230133** | **0.230267** | +0.058 % | PASS | **PASS** |
| breast / paclitaxel | 0.1503 | **0.150287** | **0.153418** | +2.083 % | NOT TESTABLE | **NOT TESTABLE** |
| breast / docetaxel | 0.1484 | **0.148430** | **0.151639** | **+2.162 %** | **FAIL** | **FAIL** |
| breast / epirubicin | 0.2301 | **0.230133** | **0.230267** | +0.058 % | PASS | **PASS** |
| breast / 5-fluorouracil | 0.3379 | **0.337882** | **0.335825** | −0.609 % | NOT TESTABLE | **NOT TESTABLE** |
| breast / cisplatin | 0.3085 | **0.308463** | **0.305113** | −1.086 % | PASS | **PASS** |

**`scored_verdict_flips = 0`. `worst_abs_ratio_pct_move = 2.6216 %`.**

**The control reproduces.** `binary_control_reproduces_published = True`, worst
deviation **4.9205646857120033e-05** (GBM temozolomide, 0.262649 against a
published 0.2626). The published values are quoted to four decimal places, so a
worst deviation of 4.9e-05 **is** exact agreement at the published precision.
This is a live re-solve, not a replay.

**The two rows that matter most both hold.** `breast/docetaxel` is this
project's only falsification and it stays **FAIL** — indeed it moves *away* from
the band, 0.148430 → 0.151639. `gbm/temozolomide` is the flagship and stays
**PASS**.

### 1.1 Where the ratio move comes from

The move is not a single mechanism, and the decomposition says so:

| tissue / drug | `f_penetration` `binary` → `area` | `f_state` `binary` → `area` | mean `e_unit` move | implied IC50 fold `binary` → `area` |
|---|---|---|---:|---|
| gbm / temozolomide | 0.467835 → 0.468040 (+0.044 %) | 0.550206 → 0.542164 | −0.716 % | 3.2608 → 3.3104 |
| gbm / cisplatin | 0.497955 → 0.493279 (−0.939 %) | 0.550206 → 0.542164 | −0.704 % | 4.5531 → 4.6715 |
| gbm / doxorubicin | 0.368825 → 0.374914 (+1.651 %) | 0.550206 → 0.542164 | −0.620 % | 6.0259 → 6.0286 |
| breast / doxorubicin | 0.446970 → 0.447939 (+0.217 %) | 0.506308 → 0.505875 | +0.052 % | 5.4740 → 5.4710 |
| breast / paclitaxel | 0.288086 → 0.294786 (+2.326 %) | 0.506308 → 0.505875 | +0.405 % | 8.1000 → 7.9508 |
| breast / docetaxel | 0.284436 → 0.291283 (+2.407 %) | 0.506308 → 0.505875 | +0.422 % | 8.1912 → 8.0349 |
| breast / epirubicin | 0.446970 → 0.447939 (+0.217 %) | 0.506308 → 0.505875 | +0.052 % | 5.4740 → 5.4710 |
| breast / 5-fluorouracil | 0.668597 → 0.665338 (−0.487 %) | 0.506308 → 0.505875 | −0.067 % | 3.7721 → 3.7952 |
| breast / cisplatin | 0.606925 → 0.601109 (−0.958 %) | 0.506308 → 0.505875 | −0.042 % | 4.1295 → 4.1741 |

Three things are visible and worth naming:

1. **`f_state` moves by the same amount for every drug in a tissue**, because it
   is a property of the oxygen field, not of the drug: GBM **0.550206 →
   0.542164** (**−0.0080**, the value `vessel_source_conservation.md` §10
   quotes), breast **0.506308 → 0.505875** (**−0.00043**). **The GBM `f_state`
   change is 19x the breast one**, which is why the GBM rows move down and most
   breast rows move up.
2. **`f_penetration` moves in BOTH directions** (−0.958 % to +2.407 %) while the
   mean exposure moves almost uniformly **down on GBM** and **up on breast**.
   The two effects partly cancel, which is why no ratio moves more than 2.62 %
   even though the source total moved 0.794 %.
3. **The largest `f_penetration` moves are the taxanes** (+2.33 %, +2.41 %) —
   the two shortest-penetration drugs in the panel, i.e. the ones whose exposure
   field is most sensitive to *where* the source sits rather than how much of it
   there is. That is the expected signature of a rasterisation fix and it is
   reassuring rather than alarming.

---

## 2. The normaliser conventions of `results/p99_convention.md`

`results/p99_convention.md` establishes that the scored ratio is invariant to
multiplicative *scale* error under any degree-1-homogeneous normaliser, and
recommends `tail90` over the shipped `p99`. **Neither claim is about source
geometry**, so what follows is a **new** measurement rather than a regeneration
of that file's published numbers — see §5 for what of it is stale and why.

All thirteen homogeneous conventions, scored on §1's solves. Worst
`|ratio move|` across the nine rows:

| convention | kind | worst ratio move, `binary` → `area` |
|---|---|---:|
| `max` | order | **+6.096 %** (gbm/doxorubicin) |
| `p99.9` | order | −2.746 % |
| `p95` | order | +2.730 % |
| **`p99` (SHIPPED)** | order | **−2.622 %** |
| `tail95` | tail | +2.495 % |
| `tail99` | tail | +2.307 % |
| **`tail90` (recommended)** | tail | **+2.287 %** |
| `sauter` | moment | +2.060 % |
| `p90` | order | +1.988 % |
| `L4` | moment | +1.944 % |
| `p50` | order | −1.922 % |
| `mean` | moment | −1.765 % |
| `rms` | moment | −1.561 % |

**Three findings.**

1. **No convention is destabilised by the source change.** The whole family sits
   in a **±6.1 %** band and twelve of the thirteen sit inside **±2.8 %**. There
   is no convention under which the area painting produces a qualitatively
   different picture from the shipped one.
2. **`max` is the outlier, and it is the outlier `p99_convention.md` predicts.**
   That report's only written justification for `p99` is "so a single grid cell
   cannot set the reference". `max` — the convention that *does* let a single
   cell set the reference — is the single worst-behaved entry here (+6.096 %),
   at more than twice the shipped convention's move. **The one argument the repo
   actually had for `p99` is independently confirmed by a perturbation it was
   never tested against.**
3. **`tail90` remains the better convention under this perturbation too**
   (+2.287 % against `p99`'s −2.622 %), as do all three conditional tail means
   relative to their corresponding order statistics (`tail95` 2.495 % vs `p95`
   2.730 %; `tail99` 2.307 % vs `p99.9` 2.746 %). **The integral-over-a-shell
   estimators beat the single-order-statistic estimators on source-geometry
   sensitivity, which is the same direction `p99_convention.md` measured for
   grid sensitivity.** That is a second, independent reason for its
   recommendation and it did not exist before this run.

---

## 3. The `F_state` and mean-O2 ladders

Real GBM section, steady oxygen solve, `dx = 30/18/12/8` (the published ladder of
`results/grid_convergence_gbm.md` §5) plus a fine probe at 6 and 5 µm. Both
paintings, both arms, all re-solved. Raw: `results/sweep/areagen_closure.json`.

### 3.1 Under the SHIPPED oxygen constants — the fix reproduces, exactly

`F_state`, `supply = 0.05`, `consumption = 0.02`:

| arm | dx 30 | dx 18 | dx 12 | dx 8 | dx 6 | dx 5 | **span (30–8)** | span (all) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **`binary`** (shipped) | 0.5041 | 0.5401 | **0.5502** | 0.5488 | 0.5390 | 0.5435 | **0.046090** | 0.046090 |
| **`area`** | 0.5424 | 0.5422 | **0.5422** | 0.5420 | 0.5420 | 0.5419 | **0.000441** | 0.000483 |

Mean tumour O2, same configuration:

| arm | dx 30 | dx 18 | dx 12 | dx 8 | dx 6 | dx 5 | **span (30–8)** |
|---|---:|---:|---:|---:|---:|---:|---:|
| **`binary`** | 0.3364 | 0.3482 | **0.3486** | 0.3482 | 0.3476 | 0.3476 | **0.012269** |
| **`area`** | 0.3474 | 0.3473 | **0.3473** | 0.3473 | 0.3472 | 0.3472 | **0.000100** |

Painted vessel area over tissue, µm²:

| arm | dx 30 | dx 18 | dx 12 | dx 8 | dx 6 | dx 5 |
|---|---:|---:|---:|---:|---:|---:|
| **`binary`** | 713,700 | 746,496 | **748,656** | 746,816 | 743,904 | 744,000 |
| **`area`** | 742,713 | 742,713 | **742,713** | 742,713 | 742,713 | 742,713 |

**This is an independent reproduction of `results/vessel_source_conservation.md`,
from a separate driver.** The `binary` `F_state` span comes out at
**0.04609020661851171** against that report's **0.04609020661851171** — the same
floating-point number — and its mean-O2 span at **0.012269196110666014** against
**0.012269196110666014**. The `area` arm's extreme values agree too: this run's
minimum **0.5419116411893163** at `dx = 5` and maximum **0.5424575844112213** at
`dx = 30` are that report's `f_state_min` and `f_state_max` exactly. The span
numbers differ (**0.000483** here against **0.000546** there) for a stated
reason: that report's ladder has ten rungs (30/24/18/16/12/10/8/7/6/5) and this
one has six, and six rungs of a ten-rung ladder cannot span more than the ten.
**There is no disagreement between the two runs.**

**Band composition at the production grid**, shipped constants, `dx = 12 µm`:

| arm | HYPOXIC | QUIESCENT | VIABLE |
|---|---:|---:|---:|
| **`binary`** | 0.0000 | 0.8996 | 0.1004 |
| **`area`** | 0.0000 | 0.9157 | 0.0843 |

The whole of the **−0.0080** `F_state` change is **1.61 % of tumour cells moving
VIABLE → QUIESCENT**. **No cell enters HYPOXIC under either painting**, so the
promotion does not create a hypoxic population that was not there.

### 3.2 Under the JOINTLY CALIBRATED constants — the fix is NOT uniformly an improvement

This is the part the nine-verdict blast radius cannot see, because all nine
verdicts run under the shipped constants. `results/oxygen_joint_grid_test.md` §1
publishes two further configurations (`supply = 194`, `consumption = 10` and
`20`), and they behave differently:

| configuration | `binary` span (30–8) | `area` span (30–8) | change |
|---|---:|---:|---|
| shipped (`supply = 0.05`, `consumption = 0.02`) | **0.046090** | **0.000441** | **−99.0 %** |
| joint-lo (`194 / 10`) | **0.157458** | **0.074181** | **−52.9 %** |
| **joint-hi (`194 / 20`)** | **0.054179** | **0.082715** | **+52.7 % — WORSE** |

`F_state` at the published `dx = 12 µm`:

| configuration | `binary` | `area` | move |
|---|---:|---:|---:|
| shipped | 0.550206 | 0.542164 | **−0.0080** |
| joint-lo | 0.524853 | 0.611228 | **+0.0863** |
| joint-hi | 0.416744 | 0.469952 | **+0.0533** |

**Stated plainly: the area painting removes essentially all of the `F_state` grid
drift on the shipped configuration, removes half of it on joint-lo, and adds half
again to it on joint-hi.** And the *level* it reports at the production grid
moves by **+0.086 / +0.053** on the two jointly calibrated configurations against
**−0.008** on the shipped one — **an order of magnitude larger**.

Two qualifications, so this is neither over- nor under-read:

* **joint-hi's `binary` span was already the smallest of the three** (0.0542,
  against 0.1575 for joint-lo and 0.0461 for shipped), and the `area` arm's
  absolute span of **0.0827** is comparable to joint-lo's **0.0742**. This is a
  good baseline being made worse, not a good result being ruined.
* **Neither jointly calibrated configuration ships.** They exist because
  `oxygen_joint_grid_test.md` calibrated them against literature oxygen data and
  found the shipped `supply` to be 4,570x off. **No published verdict uses
  them.** What moves is that report's own table.

**The consequence is a staleness, not a verdict flip** — see §5 and §6.

---

## 4. THE SUB-GRID CLOSURE ON A CONTINUOUS FIELD

`results/subgrid_closure.md` §5.0 records, as a control: *"On the shipped binary
vascular field, `binary | closed` returns `F_state` = 0.5502 / 0.5249 / 0.4167 at
`dx = 12 µm` … identical to `binary | raw` to four decimals … `beta == 1`
everywhere because the field is exactly 0 or 1."*

**That control is reproduced here and then shown to be conditional on a default
that is under review.**

### 4.1 The no-op reproduces, and its mechanism is now measured rather than argued

`binary | closed` at `dx = 12 µm`, this run: **0.5502057785990258 /
0.5248530152864102 / 0.41674365865949936** for shipped / joint-lo / joint-hi —
**bit-identical to `binary | raw`**, not merely equal to four decimals.

The mechanism, counted directly on the field at every rung of the ladder:

| painting | cells with `0 < v < 1` at dx = 30 / 18 / 12 / 8 / 6 / 5 |
|---|---|
| **`binary`** | **0 / 0 / 0 / 0 / 0 / 0** |
| **`area`** | 6,520 / 10,468 / **16,228** / 25,013 / 32,380 / 37,772 |

**There is not one partial cell anywhere in the binary field at any resolution.**
`beta` is therefore identically 1 by construction, and the no-op is a property of
the *painting*, not of the closure, the solver or the tissue. **The moment the
default stops being binary, the control is void.**

### 4.2 On the promotion candidate's field the closure is destructive

`F_state`, `geometry='area'`, raw versus closed on the same field:

| configuration | arm | dx 30 | dx 18 | dx 12 | dx 8 | dx 6 | dx 5 | **span (30–8)** |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| shipped | **`area` \| raw** | 0.5424 | 0.5422 | 0.5422 | 0.5420 | 0.5420 | 0.5419 | **0.000441** |
| shipped | **`area` \| closed** | 0.3568 | 0.3574 | 0.4379 | 0.5000 | 0.5000 | 0.5000 | **0.143203** |
| joint-lo | **`area` \| raw** | 0.6666 | 0.6326 | 0.6112 | 0.5924 | 0.5830 | 0.5777 | **0.074181** |
| joint-lo | **`area` \| closed** | 0.4769 | 0.4080 | 0.3848 | 0.3798 | 0.3877 | 0.3965 | **0.097075** |
| joint-hi | **`area` \| raw** | 0.5361 | 0.4910 | 0.4700 | 0.4534 | 0.4437 | 0.4394 | **0.082715** |
| joint-hi | **`area` \| closed** | 0.4372 | 0.3904 | 0.3720 | 0.3644 | 0.3637 | 0.3661 | **0.072865** |

**The closure reduces the `F_state` grid span on ONE of the three
configurations and increases it on the other two** — the same one-of-three
result `subgrid_closure.md` §5.1 reported on the `'subgrid'` field. But two
things are **stronger** here:

* **The shipped arm's degradation is 325-fold and the baseline is NOT
  degenerate.** §5.1 had to caveat its shipped/subgrid raw arm as "degenerate,
  not converged — 100 % of tumour cells QUIESCENT, `F_state` pinned at exactly
  0.5, no dynamic range at all", so its zero span was "not evidence of
  anything". **That caveat does not apply here.** The `area | raw` shipped arm at
  `dx = 12 µm` is **HYPOXIC 0.0000 / QUIESCENT 0.9157 / VIABLE 0.0843** — it has
  a live VIABLE population and real dynamic range — and its span of **0.000441**
  is a genuine convergence result. The closure turns it into **0.143203**.
  **This is the first un-caveated demonstration that the closure actively
  destroys grid-independence on real tissue.**
* **It moves the answer by more than it could ever fix.** At `dx = 12 µm` the
  closure changes mean tumour O2 from 0.3473 to 0.2046 on shipped (**−41.1 %**),
  0.3266 → 0.1225 on joint-lo (**−62.5 %**), 0.2112 → 0.0832 on joint-hi
  (**−60.6 %**). On the controlled 1-D problem the same closure moves the mean by
  about 6 %.

`F_state` moves at `dx = 12 µm`, closed against raw: **−0.1043** (shipped),
**−0.2264** (joint-lo), **−0.0980** (joint-hi).

### 4.3 Why — the same mechanism, measured on the new field

`beta` degenerates to the purely geometric `1/v` in the thin-saturated-source
limit. On the `area` field the median relative deviation of `beta` from `1/v` is
**0.0** at `dx = 5, 6, 8, 12` and `18 µm` — it is *exactly* `1/v`, carrying no
information about `lambda_source` or `lambda_sink` at all. Confirming that: at
every one of those five rungs `beta` is **identical across the shipped, joint-lo
and joint-hi configurations to every digit recorded** (median 2.500000 / 2.571429
/ 2.909091 / 3.428571 / 3.724138 and max 25.0 / 36.0 / 64.0 / 144.0 / 324.0 in
all three), **whose rate constants differ by up to 4,000-fold.**

(At `dx = 30 µm` the deviation is nonzero — **0.484 / 0.554 / 0.735** for
joint-hi / joint-lo / shipped — so the thin-source limit breaks down at the
coarsest rung and `beta` does become configuration-dependent there. It is the
only rung where that is true.)

The `beta` distribution on the `area` field, against the controlled sweep's own
validated range. **The `controlled sweep` column is quoted from
`results/subgrid_closure.md` §5.2 and was not re-measured here**; every other
column is from this run.

| | controlled sweep | `area`, dx = 30 | dx = 18 | **dx = 12** | dx = 8 | dx = 6 | dx = 5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| partial cells | — | 6,520 | 10,468 | **16,228** | 25,013 | 32,380 | 37,772 |
| **% of source cells partial** | — | **100.0 %** | **100.0 %** | **99.92 %** | 93.6 % | 82.8 % | 74.5 % |
| median `v` | — | 0.1089 | 0.1698 | **0.2500** | 0.3281 | 0.3889 | 0.4000 |
| % of partial cells with `v < 0.1` | — | 47.0 % | 35.5 % | **26.9 %** | 22.1 % | 16.5 % | 13.9 % |
| median `beta` (shipped) | 1.038 | 1.041 | 3.724 | **3.429** | 2.909 | 2.571 | 2.500 |
| p90 `beta` | 1.536 | 28.13 | 36.00 | **28.80** | 21.33 | 18.00 | 12.50 |
| **max `beta`** | **5.677** | **900.0** | **324.0** | **144.0** | 64.0 | 36.0 | 25.0 |
| **% of partial cells beyond the validated max** | — | 32.7 % | 38.6 % | **35.8 %** | 31.4 % | 29.1 % | 25.5 % |
| source mass beyond the validated max | — | 19.7 % | 12.9 % | **8.7 %** | 5.5 % | 4.2 % | 3.1 % |

**At the production grid, 99.92 % of source cells are partial and 35.8 % of them
demand a `beta` larger than anything the closure was ever validated at, carrying
8.7 % of the total source mass.** These are the same magnitudes
`subgrid_closure.md` §5.2 measured on the `'subgrid'` field (median `v` 0.2292
against **0.2500** here; 38.7 % beyond the validated range against **35.8 %**;
max `beta` 144.0 against **144.0**). **The failure is the same failure, on a
different continuous field, for the same measured geometric reason: a 1-D slab
cell problem applied to 2-D sub-cell discs over-corrects, and it over-corrects
hardest where `v` is smallest.**

**And it gets worse as the mesh coarsens, not better** — max `beta` runs 25 → 36
→ 64 → 144 → 324 → **900** from `dx = 5` to `dx = 30 µm`, because `v` falls as
`dx` grows. That is the opposite of what a boundary-layer closure should do.

### 4.4 Restricting the closure to its validated range does not rescue it either

`subgrid_closure.md` §5.3 tested the diagnosis by capping `beta` at the
controlled sweep's own maximum of **5.677** and letting the closure decline to
act beyond it. Repeated here on the `area` field, published ladder:

| configuration | `area` \| raw | `area` \| closed | **`area` \| closed, `beta` capped at 5.677** |
|---|---:|---:|---:|
| shipped | **0.000441** | 0.143203 | **0.000000** † |
| joint-lo | **0.074181** | 0.097075 | 0.111897 |
| joint-hi | 0.082715 | **0.072865** | 0.081650 |

† **The capped shipped arm is degenerate and its zero span is not evidence of
anything.** Its band composition at `dx = 12 µm` is **HYPOXIC 0.0000 / QUIESCENT
1.0000 / VIABLE 0.0000** — 100 % of tumour cells in one band, `F_state` pinned at
exactly 0.5000, no dynamic range. This is precisely the degeneracy
`subgrid_closure.md` §5.1 flagged for the uncapped shipped/subgrid arm, and it
must be read the same way. **Capping removes the catastrophe by removing the
signal.**

On the two configurations where the score is non-degenerate the capped closure is
**worse than doing nothing at all** on joint-lo (0.1119 against a raw 0.0742) and
**worse than not capping** on joint-hi (0.0817 against 0.0729). **There is no
version of this closure — full, or restricted to its validated range — that
improves grid-independence on more than one of the three configurations, on
either of this project's two continuous fields.**

### 4.5 The corrected conclusion

`results/subgrid_closure.md`'s headline finding — *"it does not transfer to
tissue and it stays off"* — is **not** overturned. It is **strengthened and made
unconditional**:

| claim | as written (binary default) | corrected (continuous default) |
|---|---|---|
| the closure is an exact no-op on the shipped field | **TRUE**, and reproduced bit-identically here | **VOID under `geometry='area'`**: 99.92 % of source cells partial, `F_state` moves up to **−0.2264** at `dx = 12 µm` |
| it reduces the `F_state` grid span on 1 of 3 configurations | measured on `'subgrid'` | **reproduces on `'area'`** — 1 of 3, same sign pattern |
| the shipped arm's flat ladder is not evidence (degenerate) | true on `'subgrid'` | **no longer needed** — the `'area'` raw arm is non-degenerate and the closure destroys a real 0.000441 span into 0.143203 |
| 38.7 % of partial cells beyond the validated `beta` | on `'subgrid'` | **35.8 %** on `'area'`, max `beta` **144.0** either way |
| capping `beta` rescues the shipped arm | true on `'subgrid'` | **true but degenerate** on `'area'` — capped arm is 100 % QUIESCENT |
| **operational consequence** | harmless by default AND by no-op | **harmless by default ONLY** |

**The last row is the finding.** The closure used to be safe for two independent
reasons; the promotion removes one of them. **`tmesim/physics/subgrid.py` must
not be enabled alongside the source promotion, and "it does nothing anyway" is no
longer a valid reason to be relaxed about it.**

---

## 5. What is NOT regenerated here, and is therefore stale under a promotion

Stated as a banner rather than passed over silently, following the precedent of
the fallback-rule work.

| document | what would move | regenerated? |
|---|---|---|
| `results/kill_calibration.md` (the nine ratios) | ratios by up to **2.62 %**; **no verdict** | **YES — §1.** Both arms, re-solved. |
| `results/p99_convention.md` — the convention *ranking* | up to **6.10 %** (`max`); ranking unchanged; `tail90` still beats `p99` | **YES — §2**, but as a NEW measurement (source sensitivity), not a re-run of that file's numbers. |
| `results/p99_convention.md` — stages A–G (the slab grid-error budget, `err(mean/N) = −err(N)/(1+err(N))`, the +8.3 % `mean/p99` statistic, the bootstrap stability, the perivascular percentile) | **UNKNOWN**. Those are properties of the controlled slab and of the exposure *distribution*; the source painting changes the distribution. | **NO.** Stages A–G need the controlled-slab sweep and the `dtissue` sweep, which are a different and much larger set of solves. **STALE under a promotion.** |
| `results/grid_convergence_gbm.md` §5 (`F_state` ladder, shipped constants) | span **0.046090 → 0.000441**; `dx = 12` value **0.5502 → 0.5422** | **YES — §3.1**, for the four published rungs plus two fine ones. The rest of that file (dt/CFL, §7) is untouched by the source. |
| `results/oxygen_joint_grid_test.md` §1 (joint-lo / joint-hi ladders) | `dx = 12` values move **+0.0863 / +0.0533**; joint-hi's grid span gets **52.7 % WORSE** | **PARTIALLY — §3.2** gives the four published rungs plus two fine ones for `F_state` and mean O2. That file's other quantities (its band-change counts, its per-cell statistics) are **NOT** regenerated. **STALE under a promotion, and this is the largest staleness this work found.** |
| `results/subgrid_closure.md` §5.0 | the no-op control becomes false | **YES — §4**, and the file is corrected (§7). |
| `results/tail90_falsification.md` §6.1 (3.72 % fine-rung `F_state` span) | almost certainly shrinks, since it is the same drift | **NO.** Not run. **STALE under a promotion.** |
| everything downstream of `vascular_density` not in the nine verdicts — `cell_state_impact.md`, `delivery_metrics.md`, `mgmt_spatial.md`, `heldout_patient.md`, `ranking_robustness.md`, the breast flagship files | unmeasured | **NO.** Out of scope for this brief and not attempted. **STALE under a promotion.** |

**None of these is a reason the promotion is unsound. All of them are re-runs the
promotion would owe.** The distinction matters: a stale number is a bookkeeping
debt; a moved verdict would be a scientific one. **This work found seven of the
former and zero of the latter.**

---

## 6. Is the promotion safe? — a direct answer

**Yes, on the evidence gathered here — with one condition and one caveat.**

**The bar the brief set was: zero verdicts moved in the original blast radius, so
the question is whether anything else does.** Four things were checked that the
original blast radius could not see. Here is what each returned:

1. **The nine published ratios, re-solved rather than replayed (§1).** Zero
   verdict flips, worst move 2.622 %, the `binary` control reproducing the
   published values to their published precision from a live solve. **This
   strengthens `vessel_source_conservation.md` §7 rather than qualifying it.**
   **Nothing moved.**
2. **The thirteen normaliser conventions of `p99_convention.md` (§2).** The whole
   family inside ±6.1 %, twelve of thirteen inside ±2.8 %, the ranking unchanged,
   and `tail90` still better than `p99`. **Nothing moved.**
3. **The `F_state` ladder under the shipped constants (§3.1).** The fix
   reproduces to the floating-point digit, and the band change is 1.61 % of
   tumour cells moving VIABLE → QUIESCENT with **no** new hypoxic population.
   **Nothing moved that was not already reported.**
4. **The `F_state` ladder under the jointly calibrated constants (§3.2).**
   **SOMETHING MOVED.** joint-hi's grid span gets **52.7 % worse**, and the
   `dx = 12 µm` values move **+0.0863** and **+0.0533** — an order of magnitude
   more than the shipped configuration's **−0.0080**.

**Why finding 4 does not stop the promotion.** It is a real result and it is
reported as one, but it is a **staleness**, not a verdict move, for three
reasons that are measured rather than asserted: neither jointly calibrated
configuration ships; **no published verdict is scored under either of them** (all
nine run at `supply = 0.05`); and the configuration that *does* ship shows the
fix working almost perfectly (99.0 % of the drift removed, against a
non-degenerate baseline). **What it does establish is that the promotion's
benefit is configuration-specific rather than universal, and any claim that
`geometry='area'` "fixes the grid dependence" must be stated for the shipped
constants specifically.** `results/oxygen_joint_grid_test.md` §1 would need
re-running, and until it is it is stale by up to 0.086 in `F_state`.

**THE CONDITION — and it is not optional.** `tmesim/physics/subgrid.py` must
stay disabled, and the promotion must not be bundled with enabling it. Before
this run the closure was safe twice over: default-off, *and* a provable no-op on
the shipped field. **The promotion removes the second guarantee entirely.** On
the `geometry='area'` field the closure inflates the shipped `F_state` grid span
from **0.000441** to **0.143203** — it would destroy considerably more than the
promotion buys. `results/subgrid_closure.md` has been corrected so nobody reads
its no-op control as a general property (§7).

**THE CAVEAT.** §5 lists seven documents that would be stale under a promotion
and were not regenerated, including all of `p99_convention.md` stages A–G,
`tail90_falsification.md` §6.1, and everything downstream of `vascular_density`
outside the nine verdicts. **"Safe" here means "no verdict moves and no
convention is destabilised", not "the repository is consistent the morning
after".** The promotion is a recalibration and still needs to be run as one, as
`vessel_source_conservation.md` §10 says.

**On `OxygenConfig.supply`: it was not touched, and it must not be.** It remains
**0.05**. The −0.794 % source change is visible in every table above rather than
absorbed into a constant already documented as 4,570x wrong.

**On the default: it was NOT flipped.** `discrete_vascular_field`'s `geometry`
default is still `'binary'` and `_DEFAULT_VASCULAR_MODE` is still `'discrete'`,
verified programmatically in §8. **That decision is not this file's to make.**

---

## 7. Was `results/subgrid_closure.md` generated, and how was it corrected?

**It was checked first, because a hand edit to a generated file has been silently
deleted five times in this project.**

**It is NOT generated. It is hand-written, and the correction is therefore a hand
edit.** The evidence:

* `experiments/subgrid_closure.py` writes **only** JSON. Its single output
  function is `save(stage, payload)`, which writes
  `results/sweep/subgrid_<stage>.json` and nothing else. Its `main()` calls
  `save` for the five stages and returns; there is no markdown path in the file.
* A repository-wide search for Python that writes a results markdown file
  (`_OUT_MD` / `OUT_MD` / an `os.path.join(..., "*.md")` write target) returns a
  list of generators — `grid_convergence_gbm.py`, `grid_convergence_breast.py`,
  `cfl_resolution.py`, `alpha_identifiability.py`, `kill_calibration_breast.py`,
  `real_drug_panel.py`, `niche_response.py`, `cross_tissue.py`, `combination.py`,
  `regimen.py`, `flagship_markers.py`, `dtissue_sensitivity.py`,
  `vascular_validation_breast.py`, `dose_schedule.py`, and the
  `tmesim/validation/` writers — **and `subgrid_closure.md` is not among their
  targets.**
* `scripts/reproduce_real.sh` does not reference it. Nothing in `scripts/` or
  `pyproject.toml` does.
* Every reference to the string `subgrid_closure` in `experiments/` and
  `tmesim/` is a **docstring citation**, not a write.

This is the opposite of `results/fstate_mechanism.md`, which commit `65ac848`
("Emit the corrected F_state mechanism from the generator, not as a hand edit")
moved into its generator precisely because it *is* generated.

**What was changed in `results/subgrid_closure.md`** — four edits, no number
deleted, no published value overwritten:

1. **A banner in the header block**, stating that one claim in the file is
   conditional on a default under review, naming which, and pointing at §4 here.
2. **§5.0** — the no-op paragraph is left exactly as written and a block is added
   after it recording that the no-op is a property of the binary painting, that
   the binary field has **0** partial cells while `geometry='area'` has 99.92 %,
   and what the closure does on the promotion candidate's field.
3. **§6** — a fourth numbered finding, stating that the file's third finding is
   now unconditional (it reproduces on both continuous fields) while its first
   control is now conditional, and that default-off is the only remaining
   guarantee.
4. **§8 Reproduce** — the new command.

**The uncommitted `pytest` correction already present in that file's §7 was left
untouched.** It is someone else's edit and is not mine to revert.

---

## 8. Preconditions and hygiene

Every constant this work was forbidden to move, **checked programmatically after
the fact** (`results/sweep/areagen_alpha.json`, `--stage alpha`):

| constant | required | measured | ok |
|---|---|---|---|
| `calibrate_alpha()` | 148.60815765445545 | **148.60815765445545** | **YES** |
| `OxygenConfig.supply` | 0.05 | **0.05** | **YES** |
| `OxygenConfig.consumption` | 0.02 | **0.02** | **YES** |
| `OxygenConfig.hypoxia_threshold` | 0.2 (band edge) | **0.2** | **YES** |
| `real_fields._DEFAULT_VASCULAR_MODE` | `'discrete'` | **`'discrete'`** | **YES** |
| `discrete_vascular_field(geometry=…)` default | `'binary'` | **`'binary'`** | **YES** |
| `VESSEL_RADIUS_UM` | 7.5 | **7.5** | **YES** |
| `LINK_DISTANCE_UM` | 20.0 | **20.0** | **YES** |
| production `DX_UM` | 12.0 | **12.0** (`heldout_patient`, `mgmt_scalar`) | **YES** |

The band edges `(0.2, 0.4)` and the `PDConfig` multipliers `(0.35, 0.5, 1.0)`
were not modified; they are read from `OxygenConfig` / `PDConfig` by the shared
`update_cell_states` path in every arm above.

**`pytest`: 0 failures, before and after.**

| when | result |
|---|---|
| **before** this work started | **950 passed**, 8 warnings, 349.65 s — the stated baseline, matched exactly |
| **after** the edits to `experiments/subgrid_closure.py` and `results/subgrid_closure.md` | **961 passed**, 8 warnings, 578.48 s |

**The count rose by 11 and the reason is not this work.** Another agent was live
on `tmesim/molecular/two_route.py` and `tests/test_two_route.py` throughout;
`git diff --stat` shows that test file gained **131 lines** during the session,
and it now reports **64 passed** in isolation. **No test was added, removed,
skipped or modified by this work — `tests/` is untouched here** — and **zero
tests failed in either run**, so nothing is attributable to this work in either
direction. Both runs report the same 8 pre-existing, unrelated warnings: one
matplotlib `PyparsingDeprecationWarning` and seven `DeprecationWarning`s from
`experiments/regimen.py:423`.

`tests/test_propagation.py`, which scans `experiments/` for undeclared solver
callers and reimplemented solvers, passes with the new
`experiments/area_conservation_regeneration.py` present and with the
`painting='area'` branch added to `experiments/subgrid_closure.py`.

### 8.1 What this work does NOT establish

* **One section, one tumour type, for §3 and §4.** Every oxygen number is
  `data_cache/xenium_gbm`. §1 and §2 cover both tissues.
* **One grid for §1 and §2.** The ratios and conventions were regenerated at the
  production `dx = 12 µm` only. Their grid dependence under `geometry='area'` is
  not measured here.
* **§2 is a source-sensitivity measurement, not a re-run of `p99_convention.md`.**
  That file's stages A–G are untouched and stale (§5).
* **Nothing here measures whether `geometry='area'` is more physically correct.**
  It measures that it conserves painted area, that the published verdicts survive
  it, and what it does to the closure. `results/vessel_source_conservation.md`
  §0–§2 carries the physical argument; this file carries the blast radius.
* **The `area_matched` arm was not re-run.** `vessel_source_conservation.md` §4.1
  establishes it recovers only 24.10 % of the drift; nothing here revisits that.

---

## 9. Reproduce

```
./.venv/bin/python -u -m experiments.area_conservation_regeneration --stage alpha
./.venv/bin/python -u -m experiments.area_conservation_regeneration --stage closure
./.venv/bin/python -u -m experiments.area_conservation_regeneration --stage ratios
```

`alpha` runs in under a second and needs no tissue. `closure` checkpoints to
`results/sweep/areagen_closure_partial.json` after every rung and is restart-safe;
it re-solves **69** steady oxygen problems and takes about 8 minutes, the
`dx = 5 µm` rung dominating. `ratios` checkpoints per panel to
`results/sweep/areagen_ratios_dx12.jsonl`, re-solves 18 transport problems
(2 tissues x 2 arms x 3 or 6 drugs) and takes about 20 minutes; the normaliser
family in §2 is scored off those same solves and costs nothing extra.

Restarting `ratios` from a checkpoint replays the panels but **not** the captured
exposure fields, so §2 needs a clean run of any panel it covers. The log says so
when it happens.
