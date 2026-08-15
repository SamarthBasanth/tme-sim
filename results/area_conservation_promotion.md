# The area-conserving vessel painting is now the default — re-verified on the post-two-route engine first, and zero verdicts moved

> ### One-line answer
>
> **`discrete_vascular_field(geometry='area')` is the shipping default as of
> 2026-08-14.** The binary centre-in-disc painting deposited a whole `dx²` grid
> cell of vessel per centre-in-disc vote, making the total painted vessel area a
> function of the mesh — **4.7312 %** span, non-monotone, over `dx = 30…5 µm` on
> the real GBM section. The area painting paints the same discs at the same
> radius around the same endothelial cells and conserves the total **exactly**:
> **742,713.0 µm² at every rung, span 6.3e-14 %**, which is exact rather than
> convergent because the sub-pixel lattice is anchored in physical space.
>
> **EVERY PIECE OF THE PRE-FLIP EVIDENCE REPRODUCED ON THE CURRENT ENGINE, TO
> THE PUBLISHED DIGIT, BEFORE ANYTHING WAS CHANGED.** That was the gating
> condition and it is not a formality: `results/vessel_source_conservation.md`
> and `results/area_conservation_regeneration.md` were both measured at
> `alpha = 148.60815765445545`, and commit `3de97a1` moved `alpha` to
> **22.223477711172563**, changing every drug's compiled diffusivity and
> therefore every transport solve behind the nine verdicts. Re-measured here on
> the post-flip tree: painted-area span **4.731213 %** vs a published 4.7312 %;
> conserved value **742,713.0 µm²** at all ten rungs; `F_state` span
> **0.046090 → 0.000546**, **98.82 %** of the drift removed; `area_matched`
> **24.10 %**; `beta` max **144.0** on **16,228** partial cells (**99.92 %**).
> Not one of those numbers moved.
>
> **ZERO OF THE NINE PUBLISHED VERDICTS MOVED, re-solved live at the new
> alpha.** Worst ratio move **−2.603 %** (GBM cisplatin), against **−2.622 %**
> measured pre-flip — the prediction was zero and it held. The project's only
> falsification (breast docetaxel **FAIL**) survives and goes **deeper**: 3.94 %
> outside its band → **5.84 %** outside. The flagship (GBM temozolomide
> **PASS**) survives. All thirteen normaliser conventions of
> `results/p99_convention.md` move by at most **3.839 %** and the ranking is
> unchanged — `tail90` (2.054 %) still beats the shipped `p99` (2.603 %).
>
> **`alpha` DID NOT MOVE.** `calibrate_alpha()` returns
> **22.223477711172563** before and after, bit-exactly. It is a closed-form
> solve on the compiled anchor drug and the Primeau half-distance; **no vascular
> field enters it**. Checked programmatically (§7), not assumed.
>
> **WHAT IT COST, stated rather than buried.** The production source total moves
> **−0.794 %** (748,656 → 742,713 µm² at `dx = 12 µm`) and production `F_state`
> moves **0.5502 → 0.5422**. This is unavoidable, not a design choice: the binary
> total is a *different number at every rung*, so "conserve the shipped total" is
> not a well-posed constraint, and the `area_matched` arm that does conserve it
> recovers only **24.10 %** of the drift because **the total IS the defect**.
> `OxygenConfig.supply` was **not** re-tuned to absorb it.
>
> **THE CONDITION, NOW ENFORCED IN CODE.** On the binary field
> `tmesim/physics/subgrid.py` was safe *twice over* — default-off **and** a
> provable no-op, because there was not one partial cell in the field at any
> resolution. The promotion removes the second guarantee: 99.92 % of source cells
> are now partial and `beta` reaches **144.0** against a validated maximum of
> **5.677**. Enabling the closure on the field that now ships would inflate the
> `F_state` grid span **0.000441 → 0.143203, 325-fold** — three times worse than
> the original defect. Its default-off status is **load-bearing**, and §6 records
> where that is now written in the code rather than only in a results file.

**Status: a DEFAULT CHANGE, plus the re-verification that gated it.** Changed
under `tmesim/`: exactly two files — `reconstruct/vessel_network.py` (the
`geometry` default, plus docstrings) and `physics/subgrid.py` (**docstrings
only — its `enabled=False` default is untouched**). No band edge, no
`OxygenConfig.supply` / `consumption` / `hypoxia_threshold`, no
`VESSEL_RADIUS_UM`, no `LINK_DISTANCE_UM`, no `SOURCE_SUBPIXEL_PITCH_UM`, no
`DX_UM`, no `_DEFAULT_VASCULAR_MODE` and no alpha was changed.

* Generator: `experiments/area_conservation_promotion.py`
  (`--stage hygiene | area | fstate | beta | ratios`)
* Code: `tmesim/reconstruct/vessel_network.py`, `tmesim/physics/subgrid.py`
* Tests: `tests/test_vessel_source.py`, `tests/test_vessel_network.py`,
  `tests/test_time_resolved_pk.py`, `tests/test_app_server.py`,
  `tests/data/static_source_reference.json`
* Raw: `results/sweep/acprom_hygiene.json`, `acprom_area_ladder.json`,
  `acprom_fstate_ladder.json`, `acprom_beta.json`, `acprom_ratios_dx12.json`
  (+ `.jsonl`)
* Section: `data_cache/xenium_gbm`, 40,887 cells / 23,812 tumour / 4,295
  endothelial → 3,091 vessel objects. Breast rows from
  `data_cache/xenium_breast_janesick/rep1`.
* Predecessors: `results/vessel_source_conservation.md` (landed it opt-in),
  `results/area_conservation_regeneration.md` (measured the blast radius),
  `results/two_route_promotion.md` (the flip this one was deliberately held back
  from)

---

> ## STATUS — EVERY NUMBER IN THIS FILE IS REPORTED IN TWO ARMS
>
> **This file follows the two-arm convention `results/mechanism_attribution.md`
> and `results/vessel_profile.md` established.** The difference from its two
> predecessors is that the decision is now **made**, so the arms have swapped
> roles:
>
> | arm | painting | painted vessel area | status |
> |---|---|---|---|
> | **`binary`** | centre-in-disc YES/NO, deposits a whole `dx²` cell | **grid-dependent**: 713,700 → 748,800 µm² over `dx = 30…5 µm`, **4.7312 %** span | **the SUPERSEDED default** — still reachable by name, still bit-identical |
> | **`area`** | the same discs, radius and endothelial cells, measured as covered-area FRACTION | **742,713.0 µm² at every rung**, span **6.3e-14 %** | **THE SHIPPING DEFAULT** |
>
> **No binary number in this file is overwritten and none is deleted.** Every
> table carries both arms. History is preserved the same way in the code (the
> binary branch is byte-for-byte unchanged and reachable by keyword), in the
> tests (assertions inverted, never removed) and in
> `tests/data/static_source_reference.json` (superseded values kept under a dated
> key).
>
> **What to quote.** Name the painting. "The `F_state` grid span is 0.000546" is
> now a statement about the shipped engine — but say so, because nine of the ten
> published files in `results/` still quote the binary value.
>
> **`'area'` is not `'subgrid'`.** `'subgrid'` paints one footprint per vessel
> *object* at a data-derived radius and deposits ~347k µm² against ~744k;
> `results/delivery_source_geometry.md` rejected it after it flipped three of
> seven scored rows. It remains opt-in and was **not** promoted. The argument for
> one is not an argument for the other.

---

## 0. Why this was a separate change, and what gated it

`results/two_route_promotion.md` records the hold-back verbatim: *"vessel
geometry='area' (separately evidenced; two consequential flips in one change
makes attribution impossible)."* Two-route landed at `3de97a1` and settled. This
is the second flip, run on its own.

The cost of holding back is that **all the supporting evidence was measured on
an engine that no longer exists.** `alpha` went
`148.60815765445545 → 22.223477711172563`; temozolomide's compiled diffusivity
went `277.985 → 42.4245 µm²/s`. Every transport solve behind the nine verdicts
is a different solve. **Evidence measured on a superseded engine is not evidence
about this one**, so the first thing this work did was re-measure it, on the
unmodified tree, before touching the default.

### 0.1 The trap that re-verification walks into, and how it was avoided

`experiments/delivery_source_geometry._source_field` and
`experiments/vessel_source_conservation.source_field` both build their `binary`
control by calling `discrete_vascular_field(domain, vessels)` **with no geometry
keyword** — i.e. by the default. That is correct while the default is `'binary'`
and silently catastrophic the moment it is not: *the control arm would become the
treatment arm*, every measurement would read zero change, and the promotion would
appear to be free.

Every arm in `experiments/area_conservation_promotion.py` passes `geometry=`
explicitly, including the binary control, so the identical script measures the
identical thing on both sides of the flip. `stage_hygiene` asserts it
(`binary_arm_survives_flip: true`). **Three of the eleven test failures the flip
produced were exactly this bug sitting in the test suite** (§5).

---

## 1. The pre-flip evidence, re-verified on the current engine

Run on the unmodified tree at `alpha = 22.223477711172563`, before the default
moved. **The gate was: if this does not reproduce, stop.**

| claim | published (pre-two-route) | **re-measured (post-two-route)** | reproduced? |
|---|---:|---:|:--:|
| binary painted-area span, `dx = 30…5 µm` | 4.7312 % | **4.731213 %** | ✓ |
| area painted-area span | 6.3e-14 % | **6.27e-14 %** | ✓ |
| conserved painted area, every rung | 742,713.0 µm² | **742,713.0 µm²** | ✓ |
| production total change at `dx = 12 µm` | −0.794 % | **−0.7938 %** | ✓ |
| `binary` `F_state` span | 0.04609 | **0.046090** | ✓ |
| `area` `F_state` span | 0.000546 | **0.000546** | ✓ |
| drift removed | 98.82 % | **98.82 %** | ✓ |
| `area_matched` drift removed | 24.10 % | **24.10 %** | ✓ |
| mean-O₂ span, binary → area | 0.012269 → 0.000117 | **0.012269 → 0.000117** | ✓ |
| partial cells at `dx = 12 µm` | 16,228 (99.92 %) | **16,228 (99.92 %)** | ✓ |
| `beta` max on the area field | 144.0 | **144.0** | ✓ |
| partial cells beyond validated `beta` | 35.8 % | **35.77 %** | ✓ |

**Nothing moved.** This is expected and it is worth saying why rather than
treating the agreement as luck: the painted-area ladder is **pure geometry** (no
PDE at all), and the `F_state` ladder is the **steady oxygen solve**, which reads
`vascular_density` and `rho_tum` and never touches a drug. The two-route flip
changed compiled *diffusivities*, which enter the *transport* solve. The only
quantity that could move is the nine verdicts — and that is the one thing §2
re-solves rather than assumes.

### 1.1 The painted-area ladder

| `dx` (µm) | disc / cell | **binary** (µm²) | **area** (µm²) | binary support | area support |
|---:|---:|---:|---:|---:|---:|
| 30 | 0.196 | 713,700 | **742,713.0** | 0.0256 | 0.2105 |
| 24 | 0.307 | 748,800 | **742,713.0** | 0.0268 | 0.1629 |
| 18 | 0.545 | 746,496 | **742,713.0** | 0.0268 | 0.1216 |
| 16 | 0.690 | 745,728 | **742,713.0** | 0.0267 | 0.1090 |
| **12** | **1.227** | **748,656** | **742,713.0** | **0.0268** | **0.0838** |
| 10 | 1.767 | 739,200 | **742,713.0** | 0.0265 | 0.0726 |
| 8 | 2.761 | 746,816 | **742,713.0** | 0.0268 | 0.0613 |
| 7 | 3.606 | 741,517 | **742,713.0** | 0.0266 | 0.0558 |
| 6 | 4.909 | 743,904 | **742,713.0** | 0.0267 | 0.0505 |
| 5 | 7.069 | 744,000 | **742,713.0** | 0.0267 | 0.0454 |

> **binary 4.731213 % — area 6.27e-14 %.** Ladder mean of the binary arm is
> 741,881.7 µm², so the conserved value sits **+0.112 %** from the mean and
> **−0.794 %** from the production rung. The production rung is near the **top**
> of the quantisation lottery; that is the entire origin of the −0.794 %.

**The support is NOT conserved and cannot be** — that is the honest cost, and it
is a property of what a finite-volume representation of a sub-grid object
*means*, not a defect of the fix. A cell containing any part of a vessel is
nonzero, so at `dx = 30 µm` the source spreads over 21 % of tissue against 4.5 %
at `dx = 5 µm`. **The fix conserves the source's integral, not its placement**,
so distance-resolved perivascular claims still need `dx ≤ 8 µm` exactly as
`results/vessel_promotion.md` requires today.

### 1.2 The `F_state` ladder, three arms

| `dx` (µm) | **binary** | **area** | **area_matched** |
|---:|---:|---:|---:|
| 30 | 0.5041 | 0.5424 | 0.5124 |
| 24 | 0.5392 | 0.5425 | 0.5474 |
| 18 | 0.5401 | 0.5422 | 0.5455 |
| 16 | 0.5333 | 0.5421 | 0.5449 |
| **12** | **0.5502** | **0.5422** | 0.5471 |
| 10 | 0.5305 | 0.5420 | 0.5387 |
| 8 | 0.5488 | 0.5420 | 0.5454 |
| 7 | 0.5407 | 0.5419 | 0.5407 |
| 6 | 0.5390 | 0.5420 | 0.5428 |
| 5 | 0.5435 | 0.5419 | 0.5427 |
| **span** | **0.046090** | **0.000546** | **0.034982** |
| **drift removed** | — | **98.82 %** | **24.10 %** |

Band mix at the production `dx = 12 µm`: HYPOXIC 0.0000 / QUIESCENT 0.8996 /
VIABLE 0.1004 (binary) → HYPOXIC 0.0000 / QUIESCENT **0.9157** / VIABLE
**0.0843** (area). **1.61 % of tumour cells move VIABLE → QUIESCENT and no new
hypoxic population appears.**

**`area_matched` is why the −0.794 % cannot be finessed away.** It holds the
tissue-integrated total EXACTLY at the shipped per-rung value and lets only the
distribution move — the strictest possible reading of "conserve the shipped
total while removing the grid dependence" — and it recovers **24.10 %** against
the raw fix's 98.82 %. Rescaling to a total that itself spans 4.73 % re-imports
the whole defect through the gauge. **Conserving the total and removing the
drift are in direct conflict, because the total is the defect.**

### 1.3 The remaining 1.18 % is not reachable by refining anything

The drift is a **product**, not a sum:

```
F_state drift  =  (field movement)  ×  (band-edge gain)
     0.046090  =      0.012269      ×        3.76
```

The fix attacks the left factor and removes **99.05 %** of it (mean-O₂ span
0.012269 → 0.000117). It removes **0 %** of the right factor: the band-edge
amplifier is a hard 0.4 ceiling cutting through a tumour population that spans
only 0.13 in O₂, and no rasterisation can address it. **The residual 0.000546 is
that undiminished gain acting on the irreducible field movement the PDE
discretisation and the cell-to-grid sampler still contribute.** The response to
the amplifier is still the one `results/fstate_mechanism.md` §5 named — stop
using hard bands — and **nothing here implements that.**

---

## 2. THE BLAST RADIUS, re-solved on the post-flip engine

The nine published 3-D:2-D rows, both arms, **re-solved live** in one process off
one reconstruction with one alpha, through
`experiments/delivery_source_geometry.py`'s own `_run_panel` /
`score_ratio_row` path. That file was **not edited** (it is owned elsewhere); its
`_source_field` is patched in-process only, and the patch names every geometry
explicitly (§0.1).

> **ZERO of the nine rows change verdict. Zero of the seven scored rows flip.**
> Worst ratio move **−2.603 %** (GBM cisplatin), against **−2.622 %** measured
> pre-flip. **The prediction was zero and it held.**

| row | published (pre-2R) | shipping (post-2R) | **`binary`** | **`area`** | move | `binary` verdict | **`area` verdict** | moved? |
|---|---:|---:|---:|---:|---:|---|---|:--:|
| GBM temozolomide | 0.2626 | — | 0.256491 | **0.251786** | −1.834 % | PASS | **PASS** | no |
| GBM cisplatin | 0.2791 | — | 0.281212 | **0.273892** | **−2.603 %** | PASS | **PASS** | no |
| GBM doxorubicin | 0.2077 | — | 0.208952 | **0.208109** | −0.403 % | PASS | **PASS** | no |
| breast doxorubicin | 0.2301 | 0.2330 | 0.233006 | **0.233291** | +0.122 % | PASS | **PASS** | no |
| breast paclitaxel | 0.1503 | 0.1747 | 0.174742 | **0.177020** | +1.304 % | NOT TESTABLE | NOT TESTABLE | no |
| breast docetaxel | 0.1484 | 0.1520 | 0.152033 | **0.154806** | +1.824 % | **FAIL** | **FAIL** | no |
| breast epirubicin | 0.2301 | 0.2330 | 0.233006 | **0.233291** | +0.122 % | PASS | **PASS** | no |
| breast 5-fluorouracil | 0.3379 | 0.3355 | 0.335529 | **0.333269** | −0.674 % | NOT TESTABLE | NOT TESTABLE | no |
| breast cisplatin | 0.3085 | 0.3105 | 0.310521 | **0.306914** | −1.162 % | PASS | **PASS** | no |

### 2.1 The control, the false red, and how it was made to explain itself

`delivery_source_geometry.PUBLISHED_RATIO` holds **pre-two-route** values, so a
**correct** run on this engine must disagree with it: the harness reports
`binary_control_reproduces_published: false`, worst deviation 2.44e-02. **That is
commit `3de97a1`'s documentation debt, not this change's**, and reading it as a
control failure would be wrong.

**A red flag that needs a footnote to be dismissed is a trap**, so it was not
left as one. `experiments/area_conservation_promotion.py` now prints a control
summary that names both comparisons and says which one is the control:

```
=== CONTROL SUMMARY =======================================
vs PRE-two-route PUBLISHED_RATIO : False (worst dev 2.44e-02)
    -- EXPECTED False, see binary_control_vs_pre_two_route_note
vs CURRENT shipping (two_route_promotion.md 5.2 Arm B) : True
    (worst dev 4.21e-05 over 6 rows)
    breast|doxorubicin      shipping 0.2330 | binary 0.233006 | dev 5.53e-06
    breast|paclitaxel       shipping 0.1747 | binary 0.174742 | dev 4.21e-05
    breast|docetaxel        shipping 0.1520 | binary 0.152033 | dev 3.27e-05
    breast|epirubicin       shipping 0.2330 | binary 0.233006 | dev 5.53e-06
    breast|5-fluorouracil   shipping 0.3355 | binary 0.335529 | dev 2.93e-05
    breast|cisplatin        shipping 0.3105 | binary 0.310521 | dev 2.08e-05
===========================================================
```

**The live control is green and it is an assertion, not a print.** The six breast
rows as re-published at the current alpha in `results/two_route_promotion.md`
§5.2 Arm B are reproduced by this run's independent `binary` arm to **worst
deviation 4.21e-05**, inside the 5e-5 their four-decimal precision supports.
They are **re-solved here, not replayed from a cache**, so this is a live check
that the pipeline still produces the shipping numbers. If it ever stops, the
stage raises `AssertionError` and says explicitly that this is *not* the benign
pre-two-route staleness.

The three GBM rows are **deliberately excluded** from the control:
`two_route_promotion.md` scored only the breast panel, so no post-flip published
GBM value exists, and inventing one would manufacture a control rather than apply
one.

**Why `PUBLISHED_RATIO` was not simply updated.** It would be the cleaner fix if
it were one constant, but it is not. The same pre-flip values are duplicated in
`experiments/delivery_source_geometry.py`, `experiments/grid_convergence_gbm.py`
and `experiments/grid_convergence_breast.py` — three files owned elsewhere —
so updating the one this harness reads would leave the other two disagreeing,
which is worse than a documented staleness. And several reports compare against
that constant *precisely because it is the historical published record*; mutating
it would destroy the comparison they exist to make. **Updating it is the
two-route flip's debt to settle, across all three files at once, and it is
recorded in §8 as outstanding.**

### 2.2 The falsification goes deeper, again

`breast docetaxel` is the only FAIL in the published record, and
`results/delivery_source_geometry.md` rejected the sub-grid geometry partly
because it *removed* it. Against a band top edge of **0.1462692**:

| | ratio | distance outside the band |
|---|---:|---:|
| `binary` (superseded default) | 0.152033 | 3.94 % |
| **`area` (ships)** | **0.154806** | **5.84 %** |

**This is the fifth independent physics correction in this project to deepen its
only falsification rather than rescue it.** A change that flattered the panel
would be the thing to distrust.

### 2.3 The normaliser family — nothing destabilised

All thirteen homogeneous conventions of `results/p99_convention.md`, re-scored on
these same solves (`decompose` reads only `e_unit`, so they cost no extra solve):

| convention | worst ratio move | | convention | worst ratio move |
|---|---:|---|---|---:|
| `max` | **+3.839 %** | | `sauter` | +1.823 % |
| `p99.9` | −3.011 % | | `mean` | −1.761 % |
| **`p99` (shipped)** | **−2.603 %** | | `tail99` | +1.753 % |
| `p95` | +2.544 % | | `L4` | +1.623 % |
| `tail95` | +2.164 % | | `rms` | −1.531 % |
| **`tail90`** | **+2.054 %** | | `p90` | +1.956 % |
| `p50` | −1.826 % | | | |

The whole family sits inside **±3.9 %**, twelve of thirteen inside ±3.1 %, and
the ranking is unchanged. `max` is worst, which is the single-cell order
statistic that report already names as the fragile one. **`tail90` (2.054 %)
still moves less than the shipped `p99` (2.603 %)**, so that report's
recommendation survives a source change as well as a grid change.

---

## 3. What the promotion cost, at production

| quantity, `dx = 12 µm` | `binary` (superseded) | **`area` (ships)** | change |
|---|---:|---:|---:|
| painted vessel area | 748,656 µm² | **742,713.0 µm²** | **−0.794 %** |
| source Σ over tissue, GBM | 5199.0 | **5157.7** | −0.794 % |
| source Σ over tissue, breast | 15191.0 | **15164.5** | −0.174 % |
| source support, GBM | 0.0268 | **0.0838** | **×3.1** |
| source support, breast | 0.0567 | **0.1533** | ×2.7 |
| mean tumour O₂, GBM | 0.348628 | **0.347281** | −0.39 % |
| `F_state`, GBM | 0.5502 | **0.5422** | **−0.0080** |
| IFP max in tissue | 12.932829 mmHg | **13.464377 mmHg** | +4.110 % |
| IFP mean in tissue | 0.128387 mmHg | **0.133624 mmHg** | +4.079 % |
| peak Darcy speed | 0.479831 µm/s | **0.446766 µm/s** | **−6.891 %** |

The dose change is **section-dependent** — −0.794 % on GBM, −0.174 % on breast —
which is what it should be: each is the distance from that section's own `dx =
12 µm` quantisation draw to its own conserved value, and those draws are
independent.

**The IFP moving in two directions at once is expected, not a puzzle.** A
smoother source builds a slightly *higher, flatter* pressure plateau but with
*gentler gradients*, and `v = -K grad p` reads the gradient. That −6.891 % speed
drop is what re-pinned the CFL guard in §5.

### 3.1 `OxygenConfig.supply` was not touched, and must not be

`results/oxygen_supply_calibration.md` establishes that `supply` is **4,570×**
too small, and that it is a *compensating* error against a sink ~140× too small.
The source term the PDE sees is `supply × vasc`, and `vasc`'s tissue integral is
exactly what this fix conserves. So this promotion moves `supply × ∫vasc` by a
factor of **0.99206**, against a standing documented miscalibration of a factor
of **4,570** in the same term — on a log scale, ~1,060× smaller than the error
already sitting there.

That is **not** an argument that 0.794 % is negligible; it changes production
`F_state` by 0.0080 and is reported as a change in every table above. It is an
argument that **the `supply` question this promotion "re-opens" was never
closed.** What must not happen is `supply` being quietly re-tuned to absorb the
0.794 %, which would recreate exactly the compensating-error structure that
report documents. **`supply = 0.05`, unchanged, verified in §7.**

---

## 4. `alpha` did not move, and that is a checked fact

`calibrate_alpha()` returns **22.223477711172563** before and after the flip,
bit-exactly (`results/sweep/acprom_hygiene.json`, `alpha_ok: true`).

This is not luck and it is not a coincidence worth being nervous about:
`calibrate_alpha` is a **closed-form** solve, `k_eff* = D_eff / L*²`, on the
compiled anchor drug and the published Primeau half-distance. **No vascular field
enters it.** `results/vessel_promotion.md` §1 makes the same point about the
`kde → discrete` promotion, which also left alpha bit-identical.

`experiments/area_conservation_promotion.py --stage ratios` **hard-asserts** this
before it will solve anything, so a promotion that silently moved alpha would
abort rather than publish.

**There is therefore no new alpha to pin.** The two values that matter are both
already pinned in `tests/test_two_route.py`: `ALPHA_SUPERSEDED =
148.60815765445545` (pre-two-route) and `ALPHA_CURRENT = 22.223477711172563`
(current). This work adds neither.

---

## 5. Every test that fired, and what was done to it

Baseline before the flip: **967 passed, 0 failed.** The flip produced **11**
failures. Ten were mine and every one was a guard doing its job; the eleventh was
not.

### 5.1 The guard that existed to catch exactly this

`tests/test_vessel_source.py::test_binary_default_is_untouched` — its own
docstring said that `vascular_density_from_cells`, hence every reconstruction and
every published result, routes through this default, *"so a silent flip of this
default would silently move all of them."* **It fired.** The flip is not silent —
it is evidenced, measured and reported here — so the guard is **inverted, not
deleted**, following the precedent of
`test_compile_drug_does_not_import_two_route` in the two-route promotion. It is
now `test_area_painting_is_the_default_and_binary_is_reachable_by_name` and it
pins that the default **is** the area painting and **is not** the binary one, so
an unannounced flip **back** would fire just as loudly.

### 5.2 The full ledger

| test | why it fired | what was done |
|---|---|---|
| `test_vessel_source.py::test_binary_default_is_untouched` | asserted `default == binary` | **INVERTED** → `test_area_painting_is_the_default_and_binary_is_reachable_by_name`. Every binary assertion kept, reached by name; new assertions pin the new default. |
| `test_vessel_source.py::test_binary_painting_does_not_conserve_painted_area` | **§0.1's bug, in the suite**: no `geometry=`, so it silently became a test of the FIX and reported the defect as cured | names `geometry="binary"`. Assertion unchanged. |
| `test_vessel_source.py::test_match_total_to_restores_the_shipped_total_exactly` | same: matched the area field to **itself**, making the test vacuous | names `geometry="binary"`. Assertions unchanged. |
| `test_vessel_network.py::test_subgrid_field_is_opt_in_and_default_is_unchanged` | asserted `default == binary` | **INVERTED** → `..._and_default_is_area_not_subgrid`. Adds an assertion that `'subgrid'` is still not the default — the claim the test was really written to protect. |
| `test_vessel_network.py::test_discrete_field_is_exactly_zero_away_from_vessels` | asserted the default is 0/1-valued | names `geometry="binary"`, **and additionally asserts the no-floor property on the new default** — the property the test is named for is not binary-specific, so the test got wider, not narrower. |
| `test_vessel_network.py::test_subgrid_paints_vessels_smaller_than_a_grid_cell` | its "legacy binary field" comparison used the default | names `geometry="binary"`. Assertions unchanged. |
| `test_time_resolved_pk.py` ×3 (bitwise stored references) | the painted total moved | **RE-PINNED** — §5.3. |
| `test_app_server.py::test_cfl_is_flagged_once_a_speed_has_been_measured` | the peak Darcy speed fell 6.891 % | **RE-PINNED** — §5.4. |
| `test_cli.py::TestArgparsePlumbing::test_prog_name_is_vivum` | `assert 'codon' == 'vivum'` | **NOT MINE** — §5.5. |

**No assertion was deleted anywhere.**

### 5.3 The bitwise re-pin, and why the new values are RIGHT

`tests/data/static_source_reference.json` states its own rule: re-pinning must
cost real work and a written justification, superseded values must be kept, and
regeneration is deliberately not scripted. Followed exactly — the two moved cases
were rebuilt by hand, old values preserved under
`_superseded_2026-08-14_binary_vessel_painting`, and a full `_repin_history`
entry added, the way `results/tumor_sigma_units.md` established.

**Why these are a re-pin and not a regression is checkable arithmetic, not a
judgement.** The transport integrals are **linear in the source total**, so they
must move by exactly the source total's change. They do:

| grid | tissue-integrated source | `mean_auc` | `auc_sum` | `final_c_sum` |
|---|---:|---:|---:|---:|
| `dx = 12 µm` (production) | 5199.0000 → 5157.7292 (**−0.7938 %**) | **−0.794 %** | **−0.794 %** | **−0.794 %** |
| `dx = 80 µm` | 104.0000 → 113.9567 (**+9.5738 %**) | **+9.574 %** | **+9.574 %** | **+9.574 %** |

Three independent integrals reproducing the source change to three decimals, on
two grids, **with opposite signs**, is not something a regression does. The signs
differ because the binary total is a different number at every rung: at `dx =
80 µm` a 7.5 µm disc is 1/36 of a grid cell so the centre-vote misses most
vessels, while `dx = 12 µm` sits near the top of the lottery.

**The failure pattern was predicted and observed exactly.** The painting reaches
a result only through a reconstruction, so only real-tissue cases can move.
`real_xenium_dx80_tmz_3dose` and `production_xenium_dx12_tmz_3dose` moved;
`default_1dose`, `weak_1dose`, `weak_3dose`, `weak_10dose_noconv` (hand-built
synthetic) and `solver_flat` (raw solver, no reconstruction) are **byte-for-byte
unchanged and were not regenerated**. `plasma_mean` — pure PK, no spatial
dependence — is **bit-identical in both moved cases**, which is the positive
control.

**What moved non-proportionally, named so it is not read as drift.**
`kill_fraction` moves −3.361 % at `dx = 12` and +19.614 % at `dx = 80` because
kill is a saturating Hill function of exposure, not a linear one.
`penetration_index` rises +1.639 % because the source is spread over 3.1× the
support, so the exposure field is less peaked.

### 5.4 The CFL guard — READ THIS BEFORE CONCLUDING IT WAS BENT TO FIT

> **THIS TEST WAS ALREADY BROKEN BEFORE THIS PROMOTION EXISTED, AND IT WAS
> BROKEN IN A WAY THAT HAD NOTHING TO DO WITH PHYSICS.** It did not test the CFL
> rule. It tested the contents of the author's `app/runs/` directory. The
> promotion did not break it; the promotion is what finally made it fire.

`test_cfl_is_flagged_once_a_speed_has_been_measured` derived both of its probe
meshes from **this** cache key's own measured speed. But when the probe `dx` is
not itself a cached key — which it never is, both probes are derived values —
`validate_spec` falls through to **tier 2** and quotes the **maximum** speed
recorded for that section at **any** mesh. So the test's "safe" probe passed only
when

```
2 × (this key's speed)  >  (max speed over every cached mesh for the section)
```

`app/runs/_cfl_speeds.json` is **git-ignored, machine-local state that
accumulates across sessions**. On this machine it already held a `dx = 12 µm`
key at **0.6336 µm/s** against this test's `dx = 30 µm` key at **0.2762 µm/s** —
a margin of **1.09×**. **On a machine that had ever run this section at a finer
mesh, or in a different order, this test was red for reasons entirely unrelated
to any code change.** It was one small perturbation from failing for anybody.

The promotion supplied that perturbation: the area painting spreads the same
source total over 3.1× the support, which smooths the pressure field and drops
the peak Darcy speed **−6.891 %**, shrinking the left-hand side below the
threshold.

**Established, not assumed** — the same controlled comparison used elsewhere in
this report: reverting the default to `'binary'` passes **3/3**, restoring
`'area'` fails **3/3**, deterministically, with the persisted speeds file
byte-identical across all six runs (it is git-ignored and was backed up and
compared).

**The fix repoints the test at the rule rather than at the cache**: both probe
meshes are now derived from the speed tier 2 will *actually quote*, so the test
asserts the CFL contract on any machine in any order. **All four assertions are
kept** — the error case, its message content, the safe case, and the new
consistency check that this key's own speed is among those recorded. **Nothing
was loosened to accommodate the promotion**; if anything the test is now stricter,
because it can no longer pass by accident of cache contents.

### 5.5 One transient failure was not mine, and is resolved

Mid-session, `tests/test_cli.py::TestArgparsePlumbing::test_prog_name_is_vivum`
failed with `assert 'codon' == 'vivum'`. That was a **user-directed product
rename** landing concurrently from another agent (`tmesim/cli.py`,
`tests/test_cli.py`, `pyproject.toml`, `README.md`, `app/` and `landing-3d/` are
all part of it). **Nothing in this change touches any CLI surface.** It resolved
itself once that work completed and does not appear in the final run.

It is recorded here only because it briefly made the flip look like it had landed
red when it had not — the same confusion §5.6 covers at larger scale.

### 5.6 The 58-failure run was an artefact, and it was bisected rather than assumed

An intermediate full-suite run returned **58 failed / 915 passed** across
`test_protocol.py`, `test_ingest.py`, `test_byo_compound.py`, `test_projects.py`,
`test_mechanism.py` and `test_gbmspace_adapter.py` — none of them modules this
work touches, and the errors were `JSONDecodeError` rather than anything physical.

**Two independent causes, neither this promotion's, and both established rather
than asserted:**

1. **A leaked environment variable.**
   `tests/test_app_server.py::TestAutoGeneratedApiKey` was leaking
   `TMESIM_STUDIO_ALLOWED_HOSTS` into the process environment. Every app
   constructed by a *later* module then picked up `TrustedHostMiddleware` and
   rejected `TestClient`'s `Host: testserver` with `400 Invalid host header` —
   which is why the symptom was a JSON decode error rather than a `401`, and why
   it appeared only when modules shared one process. Bisected and fixed by the
   agent that owns that file.
2. **Siblings rewriting `app/server/` and the CLI while pytest was mid-run.**

The controlled comparison that separated this from the promotion is the same one
used in §5.4: those six modules gave **361 passed / 0 failed under the `'area'`
default** against **354 passed / 7 errors under the reverted `'binary'`
default** — i.e. *cleaner with this change in*.

**No failure outside the ten in §5.2 was ever attributable to this work, and none
was dismissed on anybody's say-so.** The final clean run (§7) is 1006/0.

---

## 6. The condition: `subgrid.py`'s default-off is now load-bearing

**Before this promotion the sub-grid closure was safe for two independent
reasons. It is now safe for one.**

The shipped field used to be exactly 0 or 1 everywhere — **not one partial cell
at any resolution** — so `beta ≡ 1` by construction and the closure was an exact
no-op *as well as* being off by default. Anyone who flipped `enabled=True` got a
bit-identical answer.

Measured on the field that now ships, at production `dx = 12 µm`:

| | `binary` | **`area` (ships)** |
|---|---:|---:|
| source cells | 5,199 | 16,241 |
| **partial cells (`0 < v < 1`)** | **0** | **16,228 (99.92 %)** |
| `beta` median / p90 / **max** | 1 / 1 / **1** | 3.43 / 28.8 / **144.0** |
| validated `beta` maximum | — | **5.677** |
| partial cells beyond validated max | — | **35.77 %** |
| source mass beyond validated max | — | 8.71 % |

`results/area_conservation_regeneration.md` §4.2 measured what enabling it would
do to **exactly this field**: the shipped `F_state` grid span goes **0.000441 →
0.143203, a 325-fold degradation**, and mean tumour O₂ at `dx = 12 µm` moves
**−41.1 %**. The baseline it destroys is **not** degenerate (HYPOXIC 0.0000 /
QUIESCENT 0.9157 / VIABLE 0.0843 — a live VIABLE population with real dynamic
range). **Enabling the closure would destroy roughly three times more
grid-independence than this promotion buys.**

**Where that is now written, in the code, where someone might flip it** — the
brief's requirement, and `tmesim/physics/subgrid.py`'s `enabled=False` default is
itself **untouched**:

1. **The module docstring banner (line 3)** now says the default is load-bearing
   and points at the precondition section, instead of only saying the module is
   opt-in.
2. **The "STRUCTURAL PRECONDITION" section** is retitled
   `DO NOT FLIP solve_oxygen_subgrid(enabled=...) TO True. THE DEFAULT IS
   LOAD-BEARING, AND IT HAS BEEN LOAD-BEARING ONLY SINCE 2026-08-14.` It records
   that the field used to be binary, that the no-op guarantee is **gone**, the
   measured 99.92 % / 144.0 / 5.677 / 35.77 % numbers, the 325-fold consequence,
   and that **anyone changing the default owes a re-run and a written
   justification** — closing with: *"It does nothing anyway" was true until
   2026-08-14 and is no longer a valid reason to be relaxed about it.*
3. **`solve_oxygen_subgrid`'s own docstring** — the function whose keyword would
   be flipped — carries the warning inline, so it is visible at the call site
   without reading the module header.
4. **`source_reweight`'s `vascular_fraction` parameter doc** now reads *"A binary
   field gives `beta == 1` everywhere; **the field that ships is no longer
   binary**."*
5. **`vessel_network.py`'s module docstring** names the consequence from the
   other side, so the coupling is discoverable from either file.

---

## 7. Preconditions and hygiene

Every constant this work was forbidden to move, **checked rather than claimed**
(`results/sweep/acprom_hygiene.json`, measured *after* the flip):

| | value | required |
|---|---|---|
| `calibrate_alpha()` | **22.223477711172563** | unchanged ✓ |
| `OxygenConfig.supply` | 0.05 | unchanged ✓ |
| `OxygenConfig.consumption` | 0.02 | unchanged ✓ |
| `OxygenConfig.hypoxia_threshold` | 0.2 | unchanged ✓ |
| QUIESCENT ceiling | 0.4 | unchanged ✓ |
| `VESSEL_RADIUS_UM` | 7.5 | unchanged ✓ |
| `LINK_DISTANCE_UM` | 20.0 | unchanged ✓ |
| `SOURCE_SUBPIXEL_PITCH_UM` | 1.0 | unchanged ✓ |
| `DX_UM` (heldout_patient, mgmt_scalar) | 12.0 / 12.0 | unchanged ✓ |
| `real_fields._DEFAULT_VASCULAR_MODE` | `discrete` | unchanged ✓ |
| `solve_oxygen_subgrid(enabled=)` default | **`False`** | **unchanged ✓** |
| `discrete_vascular_field(geometry=)` default | **`area`** | **CHANGED — this work** |
| shipped field == `geometry='area'` | **true** | ✓ |
| shipped field == `geometry='binary'` | **false** | ✓ |
| binary control arm still distinguishable | **true** | ✓ |

`tmesim/physics/oxygen.py` was read and **not modified**. Nothing under
`tmesim/molecular/`, `tmesim/validation/`, `app/`, `landing-3d/` or `docs/` was
touched by this work.

**Test suite.** `./.venv/bin/python -m pytest -q` → **1006 passed, 0 failed.**

Baseline before this work, same command, same machine: **967 passed, 0 failed**.
The suite grew to 1006 during this session because sibling agents added tests
(the Codon rename, the Studio auth fix, the protocol routes); **this work added
no test and removed none.** Ten existing tests fired on the flip and were
re-pinned or inverted per §5.2, and **no assertion was deleted anywhere.**

Intermediate runs of 956/11 and 915/58 appear in the session history; both are
explained and disposed of in §5.4–§5.6, and neither reflects the state of this
change.

---

## 8. What is now stale, stated as a banner rather than passed over

`results/area_conservation_regeneration.md` §5 listed seven documents that a
promotion would make stale. **The promotion has happened, so they are stale
now**, not hypothetically.

| document | what is stale | regenerated? |
|---|---|---|
| **`results/oxygen_joint_grid_test.md`** | §1's joint-lo / joint-hi ladders. `dx = 12` values move **+0.0863 / +0.0533**; joint-hi's grid span gets **52.7 % WORSE** | **NO — banner-flagged only** (§8.1). The largest staleness this work found. |
| `results/kill_calibration.md` (the nine ratios) | ratios by up to 2.603 %; **no verdict** | **§2 supersedes them in two-arm form.** Both arms carried; the file itself is not edited. |
| `results/p99_convention.md` — convention *ranking* | up to 3.839 % (`max`); ranking unchanged | **§2.3**, as a new source-sensitivity measurement. |
| `results/p99_convention.md` — stages A–G | UNKNOWN. Properties of the controlled slab and of the exposure *distribution*, which the painting changes | **NO.** Needs the slab and `dtissue` sweeps. **STALE.** |
| `results/grid_convergence_gbm.md` §5 | `F_state` span 0.046090 → 0.000546; `dx = 12` value 0.5502 → 0.5422 | **§1.2**, all ten rungs. |
| `results/subgrid_closure.md` §5.0 | its no-op control is now false | corrected already by `area_conservation_regeneration.md` §7; **§6 here makes it load-bearing in code.** |
| `results/tail90_falsification.md` §6.1 | the 3.72 % fine-rung `F_state` span — almost certainly shrinks, same drift | **NO.** Not run. **STALE.** |
| everything downstream of `vascular_density` outside the nine verdicts — `cell_state_impact.md`, `delivery_metrics.md`, `mgmt_spatial.md`, `heldout_patient.md`, `ranking_robustness.md`, the breast flagship files | unmeasured | **NO.** Out of scope. **STALE.** |
| **`PUBLISHED_RATIO`, duplicated in `delivery_source_geometry.py`, `grid_convergence_gbm.py`, `grid_convergence_breast.py`** | holds **pre-two-route** values, so it reads as a red control on any current-engine run | **NO — and it is `3de97a1`'s debt, not this one's** (§2.1). All three copies must move together or not at all. This harness names the reason in its own output instead. |

**None of these is a reason the promotion is unsound; all of them are re-runs the
promotion owes.** A stale number is a bookkeeping debt; a moved verdict would be
a scientific one. **This work found seven of the former and zero of the latter.**

### 8.1 The joint-grid staleness, and what it does and does not license

Under the **jointly calibrated** constants the area painting is **not uniformly
an improvement**: it removes 99.0 % of the `F_state` grid span on the shipped
configuration and **increases it by 52.7 % on joint-hi**.

This is **bookkeeping, not science**, for three measured reasons: *neither
jointly calibrated configuration ships*; **no published verdict is scored under
either of them** (all nine run at `supply = 0.05`); and the configuration that
*does* ship shows the fix working almost perfectly against a non-degenerate
baseline.

**What it does establish, and what this report will not overclaim:** the
promotion's benefit is **configuration-specific rather than universal.** Any
claim that `geometry='area'` "fixes the grid dependence" is a claim about **the
shipped constants specifically**. `results/oxygen_joint_grid_test.md` now carries
a banner saying exactly this, at the top, as the only edit made to that file.

---

## 9. Honest limits

1. **The promotion does not make the engine grid-converged.** It removes the
   *source's* contribution to the drift. The band-edge amplifier is untouched
   (§1.3) and is the dominant remaining term.
2. **The painted support is not conserved and cannot be** (§1.1). Any
   distance-resolved perivascular claim still needs `dx ≤ 8 µm`.
3. **The −0.794 % is a real recalibration of the production source**, not a bug
   fix with no consequences. It moves production `F_state` by 0.0080, moves three
   bitwise references, and is reported everywhere rather than absorbed.
4. **Seven documents are stale** (§8) and only two were regenerated.
5. **The absolute conserved value carries the estimator's own Riemann error.**
   The shipped `pitch = 1.0 µm` sits +0.0093 % from the `pitch = 0.25 µm` value —
   ~120× smaller than the movement removed, but not zero.
6. **Conservation is exact only where `dx` is an integer multiple of the pitch.**
   True at every production rung; off-lattice the guarantee degrades to a Riemann
   error (span 0.1295 % against binary's 8.391 % — still 65× better, but not
   exact).
7. **This changes where drug enters the tissue, not how much is right.**
   `OxygenConfig.supply` remains documented as 4,570× wrong (§3.1). This
   promotion makes the source *mesh-independent*, not *correct*.

---

## 10. Reproduce

```bash
./.venv/bin/python -u -m experiments.area_conservation_promotion --stage hygiene
./.venv/bin/python -u -m experiments.area_conservation_promotion --stage area
./.venv/bin/python -u -m experiments.area_conservation_promotion --stage fstate
./.venv/bin/python -u -m experiments.area_conservation_promotion --stage beta
./.venv/bin/python -u -m experiments.area_conservation_promotion --stage ratios
./.venv/bin/python -m pytest -q
```

`--stage ratios` takes ~15 min (four panels, nine drugs, both arms, re-solved)
and checkpoints per panel to `results/sweep/acprom_ratios_dx12.jsonl`, keyed by
alpha so a stale checkpoint cannot be silently reused.

**To measure the binary arm after the flip, pass `geometry='binary'`
explicitly** — `experiments/vessel_source_conservation.py` and
`experiments/delivery_source_geometry.py` both build their control from the
default and will now silently return the treatment arm instead (§0.1).
