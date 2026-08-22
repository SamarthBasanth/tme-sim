# The HBD convention: the engine is right and the validation is not

**Verdict, one sentence.** `descriptors.py:87` is **correct as it ships** and must
not be changed; the defect is in the *scoring*, which fed the two-route law
Pruijn's donor-**hydrogen** column into a coefficient defined per donor **group**,
so the published `R² = 0.718` describes neither the shipping engine nor the
correct physics.

Nothing in `tmesim/`, `app/`, `docs/`, `landing-3d/` or `data_cache/` was changed
by this investigation. `HBD_LOGP_DECREMENT` was not touched. No assertion was
moved, no pin re-pinned, no default flipped. All measurements were taken by
monkey-patching `descriptors.compute_descriptors` in-process and restoring it,
verified by re-running the "before" arm afterwards and asserting bitwise
equality. `data_cache/` was checked SHA-256-identical before and after. Suite
after this work: **1098 passed, 0 failed** — the stated 1076 baseline plus tests
added by concurrent work; none of them mine, and none broken.

Script: `experiments/hbd_convention_audit.py` (six stages, each reproducible
standalone). Row-level `codon validate` output:
`results/hbd_convention_validate.json`.

---

## 0. The four questions, answered up front

| question | answer |
|---|---|
| Is Pruijn's `HD` the count of donor hydrogens? | **Yes.** Settled from the paper's own Methods sentence, not from RDKit's naming. §1 |
| Which count is the dimensionally correct partner for `1.34`? | **Donor heavy atoms — `Lipinski.NumHDonors`, what the engine already uses.** The project's own frozen registration says so in as many words; the measured Abraham `A` values confirm it at 1.6× lower RMSE and near-zero bias. §2 |
| Bug, approximation, or coincidence of naming? | **A coincidence of naming that produced a real provenance defect in one published number.** Not a bug in the engine. §2.3 |
| Does anything flip? | **No.** 0 flips across the 40 `codon validate` rows that can be tested (11 predicted values move, none crosses a bar); 10 breast kill rows are stale-cached and remain unmeasured. §6 |

Two things that are *not* what was expected going in, and both matter:

* **The suspected direction was backwards.** The runtime path is the correct one.
  It is the validation that is mis-conditioned.
* **The affected drugs include doxorubicin, methotrexate, cisplatin and
  temozolomide** — i.e. the calibration drug and the flagship GBM drug. The
  expectation that "doxorubicin, methotrexate and similar are not affected" is
  wrong. §5.

And one thing found in passing that is arguably more consequential than the
question asked, because it is live on the current default: **`codon validate` is
publishing stale kill ratios**. The `kill_calibration` cache key carries no drug
parameter, so 13 of the 50 rows are served from `data_cache` regardless of what
changed about the drug. Re-run cold, the shipping default's own published
temozolomide ratio is off by **2.6 %** and doxorubicin's by **1.8 %**. No verdict
moves, but a scorecard whose purpose is checkable numbers is printing numbers
that do not recompute. §6.2, §6.3.

---

## 1. Pruijn's `HD` is the count of donor HYDROGENS. Settled from the paper.

The three-compound check that started this is correct, and it is now settled by
the paper's own words rather than by inference. From the "Physicochemical
Parameters" section of Pruijn FB, Patel K, Hay MP, Wilson WR, Hicks KO, *Aust J
Chem* 2008;61(9):687–693, doi:10.1071/CH08240 — retrieved by the route already
recorded in `data_cache/pruijn2008/PROVENANCE.md` (Internet Archive snapshot of
the publisher's own open-access PDF), **verbatim**:

> "MW, **HD (count of –NH and –OH groups)**, **HA (count of N and O atoms)**, and
> apparent (macroscopic) pKa values for the side chain were calculated using
> ACD/PhysChem package (version 8.0, Advanced Chemistry Development Inc.,
> Toronto, Canada)."

Two things follow.

**First, the sentence is self-verifying.** The `HA` clause — "count of N and O
atoms" — is a claim checkable directly against the table, and it holds on every
row whose structure is resolvable, **5 of 5**. A sentence whose checkable half
reproduces exactly is a sentence whose other half can be trusted.

**Second, read literally, "–NH" is one nitrogen-bound hydrogen**, so a primary
amine contributes 2. That is exactly what the table shows.

```
compound               paper HD  NumHDonors  NHOHCount paper HA  N+O atoms  verdict
SN24280 (TPZ)                 2           1          2        6          6  HYDROGENS
NLCQ-1                        1           1          1        7          7  (cannot discriminate)
Urea                          4           2          4        3          3  HYDROGENS
Mannitol                      6           6          6        6          6  (cannot discriminate)
Acridan-9(10H)-one            1           1          1        2          2  (cannot discriminate)

paper HD == NumHDonors (heavy atoms): 3/5
paper HD == NHOHCount  (hydrogens)  : 5/5
paper HA == count of N and O atoms  : 5/5   <- checks the quoted sentence itself
```

**The caveat, stated plainly: it is a 5-row check, and only 2 of those 5 rows
carry any information.** Table S1 has 71 rows; 66 are ACSRC `SN` codes whose
structures exist only as embedded 2-D drawings
(`data_cache/pruijn2008/PROVENANCE.md`, "Structural identifiers — a real
limitation"). Of the 5 resolvable rows, `NLCQ-1`, `Mannitol` and
`Acridan-9(10H)-one` have donor groups carrying exactly one hydrogen each, so
both conventions give the same integer and they discriminate nothing. Only
**Urea (4 vs 2)** and **TPZ (2 vs 1)** discriminate, and both say hydrogens.

I nonetheless call this **settled rather than provisional**, because the decisive
evidence is the Methods sentence and its self-check, not the row count. The rows
are corroboration.

Third corroboration, from the article body: the text gives the Study A descriptor
ranges as *"with 0–3 HD and 5–11 HA"* and Study B as *"range of 60–347, HA of
2–7, and HD of 1–6"*. Table S1's HD values are `{0, 1, 2, 3, 4, 6}` — the 4 and 6
belong to urea and mannitol, which are Study B, exactly as the text implies.

---

## 2. Which count is physically right — and it is not the one Pruijn used

Decided by what `HBD_LOGP_DECREMENT = 1.34` is *defined against*, not by which
one scores better.

### 2.1 The project already registered the answer, in advance, in writing

`results/sweep/tworoute_promotion_prereg.json` — the frozen registration for the
promotion of `1.34` — states the derivation **verbatim**:

> "Abraham's LSER is log P = c + e·E + s·S + a·A + b·B + v·V with A the solute's
> hydrogen-bond ACIDITY. … giving a Δ log P donor coefficient of 3.621 per unit
> A. **The model's input is an integer donor COUNT (Lipinski.NumHDonors), not A,
> so a typical A per donor GROUP is required**; the central aliphatic-OH value
> A = 0.37 gives 3.621 × 0.37 = 1.340."

`1.34` was registered as "**A per donor GROUP**", explicitly paired with
`Lipinski.NumHDonors`, before the model was ever scored. Pairing it with
`NHOHCount` would multiply a per-group quantity by a per-hydrogen count.
**Switching the default would introduce the dimensional error, not remove it.**

### 2.2 The measured `A` values confirm it, and I checked them myself

`A` (= Σα₂ᴴ) is not a count of anything. It is a continuous, whole-molecule
descriptor back-derived by regression from measured partition and solubility
data. No integer represents it faithfully; the question is only which integer is
the least-wrong partner for a constant pinned to an aliphatic hydroxyl.

**The calibration anchor cannot discriminate.** `A = 0.37` is ethanol, which has
exactly one donor heavy atom *and* exactly one donor hydrogen, so
`1 × 1.34 = 1.34` under either convention. Everything turns on how the two
conventions extrapolate to multi-hydrogen groups.

Descriptor values pulled from the **UFZ-LSER database v3.2.1** (Helmholtz Centre
for Environmental Research, open), each carrying its own primary-source tag;
`[A1994]` means the row is attributed to **Abraham MH, Chadha HS, Whiting GS,
Mitchell RC (1994) *J Pharm Sci* 83:1085–1100** — the canonical Δlog P paper,
which is itself closed-access, so its *descriptors* are recoverable even though
its regression text is not. Verified directly against the raw query dump:

| compound | donor heavy atoms | donor H | `A` | source |
|---|---:|---:|---:|---|
| ethanol (the anchor) | 1 | 1 | **0.37** | [A1994] |
| methanol | 1 | 1 | 0.43 | [A1994] |
| **methylamine** | 1 | **2** | **0.16** | [A1994] |
| **dimethylamine** | 1 | **1** | **0.08** | [A1994] |
| trimethylamine | 0 | 0 | 0.00 | [A1994] |
| **aniline** | 1 | **2** | **0.26** | [A1994] |
| N-methylaniline | 1 | 1 | 0.17 | [A1994] |
| N,N-dimethylaniline | 0 | 0 | 0.00 | [A1994] |
| acetamide | 1 | 2 | 0.54 | [A1994] |
| benzamide | 1 | 2 | 0.49 | [A1994] |
| N-methylacetamide | 1 | 1 | 0.40 | [A1994] |
| acetic acid | 1 | 1 | 0.61 | [A1994] |
| phenol | 1 | 1 | 0.60 | [A1994] |

**Aniline's 0.26 is the whole-molecule value and already covers both N–H
hydrogens.** A per-hydrogen reading would require 0.52, which is not the
published number. The same holds for benzamide (0.49, not 0.98). The `A` table
already quoted in `two_route.py`'s docstring is therefore a **per-group** table,
and consuming it one increment per group is what the code does.

**The honest complication, which does not change the answer.** The aliphatic
amine series *does* scale roughly with hydrogen count — 0.16 → 0.08 → 0.00 as
N–H hydrogens go 2 → 1 → 0. Taken alone that looks like support for the hydrogen
convention. What kills it is the **magnitude**: the per-hydrogen increment for an
amine is ≈ 0.08, nowhere near the anchor's 0.37. Charging an –NH₂ two full
hydroxyl-equivalents is 4.6× the measured value; charging it one is 2.3× — still
wrong, but half as wrong.

Quantitatively, error in the applied penalty against the true `3.621 × A`, over
**24 compounds** — the 13 above plus water, ammonia, ethylamine, diethylamine,
triethylamine, propylamine, dipropylamine, formamide, benzoic acid, urea and
N,N-dimethylacetamide (computed here, reproducible from the same dump):

| convention | RMSE | MAE | bias |
|---|---:|---:|---:|
| **heavy atoms — `Lipinski.NumHDonors` (shipping)** | **0.90** | **0.67** | **+0.02** |
| hydrogens — `Descriptors.NHOHCount` | 1.43 | 1.01 | **+0.75** |

The shipping convention is **1.6× better on RMSE and essentially unbiased**; the
hydrogen convention carries a **+0.75 log-unit systematic over-penalty**, i.e. it
would make every amine-bearing compound roughly 5× less membrane-permeable than
the measured acidity supports. (A sibling's independent pull over a slightly
different 21-compound set gave 0.72 vs 1.31 — same direction, same magnitude.)

Two honesty notes on that table. **Water is scored `NumHDonors = 0`** — RDKit's
SMARTS requires exactly one H on the oxygen, so `O` does not match — which
*penalises* the heavy-atom convention in the RMSE above; the 0.90 is therefore
conservative against the shipping choice, not flattering to it. And **urea's `A`
is taken as 0.50** (Gunatilleka & Poole 2000, the value carried in the dump);
Abraham (2010) reports 0.83 for the same compound, and the dump records both.
Neither substitution changes the ordering.

**Residual error is large either way, and that should be said.** Even with the
correct convention, one flat `1.34` carries RMSE ≈ 0.9 log units, because `A`
varies ~7× across donor *types*: 0.08 for an amine N–H, 0.37 for an aliphatic OH,
0.60–0.61 for a phenol or carboxylic acid. The flat constant **over-penalises
amines by ~0.8 log units and under-penalises acids, phenols and amides by
~0.5–0.9**. That is a known limitation of using any integer count, and it is not
an argument for the other integer.

### 2.3 The consequence: the published `R² = 0.718` has mixed units

`experiments/two_route_transport.py::_build_rows` reads field 3 of
`experiments/diffusivity_falsification.py::PRUIJN_2008_TABLE_S1` — the paper's
`HD` column, i.e. donor hydrogens — and hands it to
`two_route.membrane_permeability_um_s`, where it is multiplied by a per-group
`1.34`. The runtime path (`descriptors.compute_descriptors` →
`compile_drug.build_profile` → the same function) hands it
`Lipinski.NumHDonors`.

**The law was scored on one quantity and is fed another.** That much of the
original finding is confirmed exactly as stated. What is new is the *direction*:
the scoring is the side that is dimensionally wrong.

So: **not a bug, and not a documented approximation. A coincidence of naming**
— `NumHDonors` and `NHOHCount` agree for hydroxyls, which is most of medicinal
chemistry, and diverge for primary amines and amides — **which produced a real
provenance defect in exactly one place: the number `0.718`.**

---

## 3. Can `0.718` be recomputed on the engine's own convention? Partly — and the residue is too wide to be useful

`--stage recompute`. This is the sharpest form of the question, and the answer is
neither "yes" nor "unrecomputable".

**45 of 71 rows can be pinned exactly with no structure at all.** Two arguments,
neither of which transcribes anything from a drawing:

1. **`HD ≤ 1` pins it.** A row with `HD = 1` has exactly one N–H/O–H hydrogen; a
   donor group holding one hydrogen holds one donor heavy atom; so
   `NumHDonors == NHOHCount == 1`. `HD = 0` likewise. That is **42 rows** (34 at
   `HD = 1`, 8 at `HD = 0`).
2. **5 rows have a resolvable named structure** (§1) — of which `NLCQ-1` and
   `Acridan-9(10H)-one` are already covered by argument 1, so the union is
   42 + 5 − 2 = **45**.

That leaves **26 free rows**, all at `HD ∈ {2, 3}`, each with a heavy-atom count
in a 2-element set (`⌈HD/2⌉ … HD`) — 6.7 × 10⁷ admissible configurations.
Enumerating the prediction at every (row, candidate) pair and searching that
space:

| configuration | R² | geo bias | within 3× |
|---|---:|---:|---:|
| **as published** (paper HD = donor hydrogens) | **0.7183** | 1.03× | 99 % |
| every free row assumed unaffected (no `–NH₂`) | 0.7151 | 1.04× | 99 % |
| every free row carries one `–NH₂` | 0.6376 | 1.28× | 91 % |
| **searched minimum** over the 6.7 × 10⁷ space | **0.4978** | — | — |
| **searched maximum** | **0.8388** | — | — |

(20 000 random draws gave 0.5014–0.8295; 40-restart coordinate ascent/descent
tightened it to 0.4978–0.8388. **R² is not monotone in the per-row shift**, so
the two named scenarios are *plausible cases*, not bounds — a mistake worth
flagging, because it is the natural thing to assume.)

**Conclusion.** The score is not unrecomputable — 45 of 71 rows, nearly two
thirds, are pinned for free. But the residual 26 rows leave R² anywhere in **0.498–0.839**, an interval
0.34 wide that *contains* the published 0.718. So `0.718` is **neither refuted
nor confirmed** on the engine's own convention, and the interval is far too wide
to serve as a validation. The earlier crude "shift everything by −1 / −2" bracket
(0.5695 / 0.3622) is not a bound either — it silently shifts the 42 rows that
provably cannot shift.

**What would actually close it:** the 26 structures, from the `ref_footnote`
synthesis papers already captured in
`data_cache/pruijn2008/pruijn2008_tableS1.csv`, or from the corresponding author
(`f.pruijn@auckland.ac.nz`), whom `PROVENANCE.md` already names as the route to
the unpublished arms. **Not** by OCR-ing the Table S1 drawings — that is
inventing inputs.

One tempting shortcut, and why I am not taking it: the free rows are TPZ
analogues, and TPZ's only donor is its 3-amino group, so it is natural to guess
that every `HD = 2` row is one `–NH₂` and land on **0.638**. **That guess is not
supported.** A `HD = 2` row is equally consistent with a *secondary* 3-alkylamino
`–NH` plus a hydroxyl — and 3-alkylamino substitution is a standard variation on
this exact scaffold. Without the structures the two cases are indistinguishable,
and 0.638 is a scenario, not an estimate. It is not quotable.

**Where `0.718` is published and would need re-labelling** — reported, not
edited, since all of these are outside this investigation's territory and several
are pinned by tests: `README.md` (lines 31 and 372), `tmesim/cli.py:632`,
`tmesim/molecular/byo.py` (3 places), `app/server/api.py:729`,
`app/server/report.py:195`, `experiments/descriptor_falsification.py:104`
(`PUBLISHED_R2`), and asserted in `tests/test_byo_compound.py` and
`tests/test_time_resolved_pk.py`. `results/descriptor_falsification.md` §7
already recommended the re-label; this file is the quantitative case for it.

---

## 4. Blast radius, measured — `calibrate_alpha` moves 32 %

`--stage alpha`:

| | value |
|---|---|
| before — `Lipinski.NumHDonors`, the shipping default | `22.223477711172563` |
| after — `Descriptors.NHOHCount` | `15.100926279598324` |
| restored-default re-run (asserted equal to "before") | `22.223477711172563` |
| ratio | **0.6795** |

The pin `22.223477711172563` appears in `results/two_route_promotion.md` §0 and
throughout the anchor prereg. **It was not re-pinned and must not be** — the
number above is what a switch *would* produce, recorded so nobody has to
rediscover it.

The mechanism is not the membrane leg. Doxorubicin's `D_eff` is essentially
invariant (`27.9255826233 → 27.9255760326`): at logP ≈ 0 with 6 donors the
transcellular route is already closed, so a seventh changes nothing that survives
the effective-medium average. Alpha moves because `binding.tissue_binding_rate` —
a separate heuristic, `score = (logP−2.5)/2.5 + (hbd−2)/4` — moves 1.47×, and
alpha is re-solved to put the same measured doxorubicin perivascular anchor back
(`L = 64.9 µm`, unchanged in both arms of §6). Note that `alpha × binding` is
then almost exactly conserved for doxorubicin (`4.971e-3 → 4.973e-3`): **the
calibration absorbs the change at the anchor drug by construction**, and what
propagates to every other drug is only its *difference* from doxorubicin.

That is also why compounds whose own count does not move at all still shift in
§6. Mitoxantrone and topotecan carry only single-hydrogen donors, so
`NumHDonors == NHOHCount` for both and their descriptors are identical in the two
arms; their `L` moves purely because `alpha` fell, lowering
`k_eff = clearance + alpha·binding` and lengthening `L = sqrt(D_eff/k_eff)`.
**`alpha` is a global gauge, so a change confined to doxorubicin's descriptors
reaches every drug in the engine.**

### 4.1 `hbd` has three consumers and only one is the Abraham law

| consumer | what it does with `hbd` | provenance |
|---|---|---|
| `two_route.membrane_permeability_um_s` | `logK = logP − 1.34·hbd` | Abraham LSER, registered |
| `permeability.bbb_penetrance` (twin in `glioma_regional_pk`) | `d_hbd = 1/(1+exp((hbd−3)/0.8))` | declared heuristic |
| `binding.tissue_binding_rate` | `score = (logP−2.5)/2.5 + (hbd−2)/4` | declared heuristic, "CONFIDENCE: low" |

The two heuristics have sigmoids centred at `hbd = 3` and `hbd = 2` with widths
`0.8` and `4`. A ±1 shift is a large fraction of those widths, which is why they
move more than the physics leg does. **Neither heuristic has any Abraham
provenance, so §2's argument does not license either count for them** — they are
declared conventions and remain so under both.

---

## 5. Which library drugs are affected, named

`--stage library`. 20 unique compounds across `GBM_DRUGS` + `BREAST_DRUGS`; 19
have a parseable SMILES (bevacizumab is the antibody marker and has none).

**5 of 19 change. 14 do not.**

| compound | panel | `NumHDonors` | `NHOHCount` | Δ | Δ log K_mem |
|---|---|---:|---:|---:|---:|
| **cisplatin** | GBM | 2 | 4 | **+2** | −2.68 |
| **methotrexate** | GBM | 5 | 7 | **+2** | −2.68 |
| **doxorubicin** | GBM | 6 | 7 | **+1** | −1.34 |
| **epirubicin** | breast | 6 | 7 | **+1** | −1.34 |
| **temozolomide** | GBM | 1 | 2 | **+1** | −1.34 |
| 5-fluorouracil, carmustine, cyclophosphamide, docetaxel, etoposide, irinotecan, lomustine, paclitaxel, palbociclib, procarbazine, ribociclib, tamoxifen, vincristine, vinorelbine | both | — | — | **0** | 0 |

**The affected set is small but it is not the unimportant part of the set.**
Doxorubicin is the single compound `calibrate_alpha()` is fit on. Temozolomide is
the GBM drug the product is built around. Cisplatin moves by two.

Compiled consequences:

| drug | `hbd` | `D_eff` µm²/s | `bbb_penetrance` | `binding` 1/s |
|---|---|---|---|---|
| cisplatin | 2→4 | 38.67 → 35.71 (**0.923×**) | 0.858 → 0.560 (**0.653×**) | ×**2.00** |
| temozolomide | 1→2 | 42.42 → 42.04 (0.991×) | 0.125 → 0.116 (0.926×) | ×1.18 |
| doxorubicin | 6→7 | 27.9256 → 27.9256 (1.000×) | 0.0206 → 0.0199 (0.965×) | ×**1.47** |
| epirubicin | 6→7 | 27.9256 → 27.9256 (1.000×) | 0.0206 → 0.0199 (0.965×) | ×**1.47** |
| methotrexate | 5→7 | 30.0605 → 30.0602 (1.000×) | 0.0225 → 0.0203 (0.903×) | ×**2.16** |

The membrane leg itself moves 22–1000× (`P_mem` for methotrexate goes
`4.14e-05 → 8.64e-08`) and almost none of it reaches `D_eff`, because for these
compounds the transcellular route is already negligible. **The two-route law is
the least sensitive of the three consumers.** Cisplatin is the only compound
whose `D_eff` moves materially, by 7.7 %.

Beyond the library: any user compound submitted through `tmesim/molecular/byo.py`
carrying a primary amine or primary amide is affected on the same terms.

---

## 6. Every verdict in `codon validate` — **0 of 50 flip**

Full suite (not `--quick`), both arms in one process, `before` first and
unpatched. (The restore itself is proved in `--stage alpha` and `--stage pruijn`,
each of which re-runs its baseline after the patch and asserts bitwise equality;
this stage does not re-run the suite a third time.)

| | FAIL | NOT TESTABLE | PASS | exit |
|---|---:|---:|---:|---:|
| **before** | 17 | 19 | 14 | 1 |
| **after** | 17 | 19 | 14 | 1 |

**VERDICT FLIPS: 0.** Eleven rows' *predicted values* moved; not one crossed a
bar.

**How much of that 50 is genuinely tested: 40.** 18 rows carry no drug and cannot
respond; 13 were served from a stale cache, of which the 3 GBM ones were re-run
cold and still do not flip (§6.2) and the **10 breast ones remain untested**
(§6.1). So the honest claim is *0 flips across 40 rows, with 10 unmeasured* — not
0 across 50. Nothing in the 40 comes close to a bar.

| anchor / row | before | after |
|---|---|---|
| `gbm-penetration` / mitoxantrone | L = 48.3 µm (0.67× band mid) | L = 58.5 µm (0.81×) |
| `gbm-penetration` / topotecan | L = 148.3 µm (2.06×) | L = 178.8 µm (2.48×) |
| `glioma-partition` / eribulin | 0.0353 (6.94× measured) | 0.0329 (6.48×) |
| `glioma-partition` / methotrexate | 0.1415 (0.79×) | 0.1393 (0.78×) |
| `glioma-partition` / temozolomide | 0.1594 (0.82×) | 0.1501 (0.77×) |
| `glioma-partition` / levetiracetam | 0.9466 (1.36×) | 0.9117 (1.31×) |
| `gist-imatinib` / imatinib | GM \|fold error\| 11.38× | 11.83× |
| `mcl-transport-coefficients` / `D_eff` (TPZ) | 49.3 µm²/s (0.67× lowest measured line) | 43.7 µm²/s (0.59×) |
| `mcl-transport-coefficients` / hindrance (TPZ) | 10.33-fold (0.86× measured) | 11.66-fold (**0.97×**) |
| `mcl-transport-coefficients` / `k` (TPZ) | k_cal 6.918e-04 /s at α = 22.2235 | 6.033e-04 /s at α = 15.1009 |
| `breast-transvascular` / methotrexate | 0.0537 (0.09× measured) | 0.0215 (**0.04×**) |

Two deserve naming. The **TPZ hindrance factor improves** (0.86× → 0.97× of
measured) while the **TPZ `D_eff` degrades** (0.67× → 0.59×) — same compound,
same run, opposite directions on two readings. That is exactly why "it scored
better" is not an argument. And **`breast-transvascular` / methotrexate degrades
2.5×**, 0.09× → 0.04× of the measured human microdialysis value; already failing,
now failing worse.

### 6.1 A hole in that answer, found and then closed

**13 of the 50 rows were not initially a valid test.** `gbm-kill-3d2d` (3 rows)
and `breast-kill-3d2d` (10 rows) are served from
`data_cache/kill_calibration_runs.json` / `..._breast_runs.json`, whose cache key
is built at `tmesim/validation/kill_calibration.py:1001` as

```python
ckey = f"{name}|{dx_um}|{DT_S}|{WINDOW_S}|{SAVE_EVERY}"
```

— **drug name and grid only. No drug parameter enters the key.** The cached
object is `e_unit`, the per-cell unit-concentration exposure field, i.e. the
output of a PDE solve that used the drug's diffusivity, binding and clearance. A
descriptor change is therefore silently ignored on a warm cache. Both caches were
fully warm (verified: SHA-256 byte-identical before and after, so no key missed
and nothing recomputed), so both arms read the *same* pre-change fields.
Cisplatin's `D_eff` moves 7.7 % and its cached kill ratio is bit-identical at
0.2746 in both arms — the masking, visible.

**This is a latent cache-invalidation defect independent of the HBD question**,
reported and not fixed (`tmesim/**` is outside territory). Anyone who changes a
descriptor, a partition constant or `HBD_LOGP_DECREMENT` and reads a green
`codon validate` off a warm cache is reading a stale answer for those 13 rows.

A further **18 rows** (`vascular-morphometry`, `breast-vascular-morphometry`,
`endothelial-cellularity`, `heldout-patient-geometry`, `ifp-wick-in-needle`,
`ps-scale`) involve no drug and are unaffected by construction. That leaves 19
live drug-dependent rows, of which 11 moved and 0 flipped. The 3 GBM cached rows
are closed cold in §6.2; the 10 breast ones are not.

### 6.2 The 13 cached rows, re-run COLD — still 0 flips, but a stale published number falls out

`gbm-kill-3d2d` was re-run with `use_cache=False` in both arms, driving
`kill_calibration.run` → `analyse` → `scorecard._kill_rows` exactly as the
scorecard does, so these are the real verdict rows. `use_cache=False` also
**writes nothing** (`kill_calibration.py:1080` guards the write with
`if use_cache:`), confirmed: both cache files are SHA-256 unchanged.

| drug | before COLD (α = 22.223478) | after COLD (α = 15.100926) | verdict |
|---|---|---|---|
| cisplatin | ratio 0.2739 (IC50 fold 4.64×) | ratio 0.2603 (4.87×) | PASS → PASS |
| doxorubicin | ratio 0.2081 (6.01×) | ratio 0.2081 (6.01×) | PASS → PASS |
| temozolomide | ratio 0.2518 (3.38×) | ratio 0.2577 (3.31×) | PASS → PASS |

**0 flips, cold.** So the §6 headline holds on the full 50 rows, not just the 37
that were live. Doxorubicin is bit-identical, which is §4's alpha-absorption
argument showing up exactly where predicted: `alpha × binding` is conserved at
the anchor drug by construction.

**The blast radius that falls out of this, and it has nothing to do with HBD.**
Compare the cold `before` arm — the unmodified shipping default — against what
`codon validate` actually printed off the warm cache minutes earlier:

| drug | cache-served (what the scorecard publishes) | cold-recomputed (correct) | error |
|---|---:|---:|---:|
| doxorubicin | 0.2044 | **0.2081** | **−1.8 %** |
| temozolomide | 0.2583 | **0.2518** | **+2.6 %** |
| cisplatin | 0.2746 | **0.2739** | +0.26 % |

Same code, same alpha, same grid, `verify_engine=False` in both — the only
difference is cache-served versus recomputed. **`codon validate` has been
publishing stale 3-D:2-D kill ratios for these rows on the current default**, by
up to 2.6 %. The cached `e_unit` fields were solved under drug parameters that no
longer hold, and the key at `kill_calibration.py:1001` cannot notice. No verdict
is affected — the bands are far wider than 2.6 % — but these are published
numbers on a scorecard whose whole purpose is that its numbers are checkable.

**Not closed: the 10 `breast-kill-3d2d` rows.** They have the same defect through
`experiments/kill_calibration_breast.py:315` and were not re-run cold here.
Whether *they* are publishing stale values is unmeasured, and that is stated
rather than assumed.

### 6.3 The cache-key defect, characterised — and it is worse than under-keying

Reported, **not fixed**: `tmesim/**` is outside this investigation's territory.
The exact patch and its test are below so it can be applied in one pass.

The GBM key omits every drug parameter. The **breast twin is stranger**:
`experiments/kill_calibration_breast.py:315` keys on

```python
ckey = f"{name}|{dx_um}|{DT_S}|{WINDOW_S}|{SAVE_EVERY}|ic50={drug.ic50:.6g}"
```

`ic50` is the **one** drug parameter that does *not* determine the cached object.
`e_unit` is the unit-concentration exposure field from
`unit_transport(work, drug, pde, velocity)`, which reads `diffusivity`,
`binding`, `clearance` and `bbb_penetrance`; `ic50` enters only afterwards, in
the kill computation. So the breast key includes the irrelevant parameter and
omits all four determining ones. That is not a missing field — it is a key that
looks parameterised and is not, which is why it survived.

**Adding the four fields would be the wrong fix**, because it goes stale the next
time a parameter is added. The structural version cannot:

```python
import dataclasses, hashlib, json

def _ckey(name: str, drug: DrugParams, dx_um: float) -> str:
    """Cache key covering EVERY field of DrugParams, by construction.

    Deliberately over-keys: `ic50` does not affect `e_unit`, and including it
    costs an occasional needless resolve. That is the safe direction. A key
    derived from the parameter set cannot silently stop covering a parameter
    the way a hand-written f-string can.
    """
    blob = json.dumps(dataclasses.asdict(drug), sort_keys=True, default=repr)
    digest = hashlib.sha256(blob.encode()).hexdigest()[:16]
    return f"{name}|{dx_um}|{DT_S}|{WINDOW_S}|{SAVE_EVERY}|{digest}"
```

And the test that makes it structural rather than remembered — it fails if *any*
field stops being covered, including one added years from now:

```python
def test_kill_cache_key_covers_every_drug_parameter():
    """A cache key that under-keys is a check that has stopped checking.

    Same failure class this repo has hit before: the guard still runs, still
    passes, and no longer constrains anything. Perturb each DrugParams field in
    turn; every perturbation MUST move the key.
    """
    base = DRUG_LIBRARY["doxorubicin"]
    k0 = KC._ckey("doxorubicin", base, 12.0)
    unmoved = []
    for f in dataclasses.fields(DrugParams):
        v = getattr(base, f.name)
        alt = (v * 1.5 + 1.0) if isinstance(v, (int, float)) else "PERTURBED"
        if KC._ckey("doxorubicin", dataclasses.replace(
                base, **{f.name: alt}), 12.0) == k0:
            unmoved.append(f.name)
    assert not unmoved, f"cache key ignores DrugParams fields: {unmoved}"
```

Note this test would have **failed on the shipping key** for every field, and on
the breast key for every field except `ic50` — i.e. it detects the defect that
was actually present, which is the bar a regression test has to clear.

**One consequence to plan for, not to avoid.** A hashed key invalidates every
existing cache entry, so the first run after the fix is fully cold — minutes, not
seconds. Given §6.2 that is not a cost, it is the *point*: the existing entries
are the stale ones. Whoever lands this should expect the published kill ratios to
move by ~1–3 % and should record that movement rather than treat it as a
regression.

---

## 7. Citation defects in this repo, verified first-hand

Four of them, found while checking §2. **No number moves** — `3.621`, `0.37` and
`1.34` all stand and reproduce against sources I retrieved and read. What is
wrong is where the repo says they came from. In summary:

| # | defect |
|---|---|
| D1 | `doi:10.1021/acsomega.8b03242` is **ACS Omega 2019;4(2):2883–2892**, not the cited `4(7):12920–12926`. Confirmed against Crossref and against PMC6648601. |
| D2 | That paper's Table 1 **has no hexadecane row**, so `a = −3.587` cannot have come from the source it is attributed to. |
| D3 | The cross-check "`Table 1 of J Comput Aided Mol Des 2017 (PMC5579141)`" is really ***J Solution Chem* 2017;46(8):1625–1638**, doi:10.1007/s10953-017-0667-y — but its Table 1 **does** carry both `0.034` and `−3.587`, so the cross-check is sound and only the journal name is wrong. |
| D4 | The compilation cross-check quotes **hexane −3.578, heptane −3.571**; the table prints **−3.599 and −3.317**. Neither reproduces. Nothing numeric depends on them, which is exactly why it would have survived. |

Full write-up, the open-access source that *does* support both coefficients, and
exact replacement text for all seven affected sites:
**`results/abraham_citation_erratum.md`**. Three of the sites are frozen
pre-registrations and are handled by erratum rather than edit; one is
`tmesim/molecular/two_route.py:206–210` and needs routing.

---

## 8. What I would do, and why

**1. Do not change `descriptors.py`. Do not touch `HBD_LOGP_DECREMENT`.** The
engine is the correct side of this discrepancy. Switching to `NHOHCount` would
pair a per-group coefficient with a per-hydrogen count, add a +0.75 log-unit
systematic over-penalty on every amine, move `alpha` 32 %, and break a frozen
registration that names `Lipinski.NumHDonors` explicitly. `tests/test_two_route.py:558`
is correct and should stay.

**2. Re-label `0.718`, and say what the honest alternative is.** Everywhere it
appears it should read *"scored on the paper's own descriptors (`HD` = donor
hydrogens), which is not the convention the engine computes."* And the follow-on
sentence should be §3's, not a substitute number: on the engine's own convention
the score is pinned for 45 of 71 compounds and lies somewhere in **0.498–0.839**
for the rest — an interval that contains 0.718 and is too wide to validate
anything. **Do not publish a replacement point estimate.** §3 lists every call
site, including the two test files that pin the string.

**3. Fix the cache key, structurally, and add the test that keeps it fixed.**
`ckey` at `kill_calibration.py:1001` omits every drug parameter, and the breast
twin at `kill_calibration_breast.py:315` includes the one parameter that does not
matter and omits the four that do. 13 of 50 `codon validate` rows are therefore a
regression test for the cache, not for the model. **Do not just add the four
fields** — derive the key from `dataclasses.asdict(drug)` so it cannot go stale
when a parameter is added, and add the test in §6.3 that fails if any field stops
being covered. Both are written out ready to apply; they are in `tmesim/**` and
`experiments/kill_calibration_breast.py`, so they need routing.

**This is not hypothetical and it is not about HBD.** §6.2 measured it: on the
unmodified shipping default, the published temozolomide 3-D:2-D ratio is
**0.2583 cache-served against 0.2518 recomputed (2.6 % off)** and doxorubicin's
is **0.2044 against 0.2081 (1.8 % off)**. Those are claim-surface numbers.
Wherever a kill-ratio row is quoted it is currently quoting a cached value that
does not reproduce, and the fix is to invalidate the cache as part of the same
change. The 10 breast rows have the same defect and were **not** measured.

**4. Record the citation erratum.** Written:
`results/abraham_citation_erratum.md`, with exact replacement text for all seven
sites. A wrong source-to-number mapping in a provenance chain that exists
precisely to make numbers checkable is a defect even when the number is right —
and D2 is not a formatting slip: the cited table does not contain the cited
coefficient. **The three frozen pre-registrations must not be edited**; the
erratum beside them is the correct form, because a pre-registration that can be
revised after the fact is not one.

**5. What I would *not* do, and why.** If one wanted the law to be more right
rather than merely self-consistent, the move is **not** a different integer — it
is a per-group `A` lookup (aliphatic OH 0.37, amine N–H 0.08, aniline NH₂ 0.26,
amide NH₂ 0.49–0.54, phenol 0.60, acid 0.61), which cuts the ~0.9 log-unit RMSE
substantially and costs nothing at runtime. That is a larger change with its own
registration and its own blast radius, and **this finding does not justify it**.
It belongs in a separate pre-registered experiment, scored against something
other than Pruijn.

### The one-line summary for the "is it immaterial?" question

**No — "immaterial" is too strong, and "material" is too strong.** The difference
is real, it is dimensionally settled in the engine's favour, it moves `alpha` by
32 % and five library drugs including the calibration drug — and it flips **zero
of fifty** external-anchor verdicts. The thing genuinely damaged is one published
number, and the damage is that it cannot be cleanly recomputed rather than that
it is wrong.
