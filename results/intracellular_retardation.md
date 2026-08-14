# The `K > 1` intracellular retardation term — built, verified, and refuted

**What this file is.** `results/two_route_efflux.md` §5.3 named one thing as the
last item standing between a measured R² = 0.718–0.725 transport model and its
becoming the default: a `K > 1` intracellular retardation term, with a measured
same-cell-line target (Choi 2011). This file builds that term, verifies its
composition against a direct transient PDE solve, and scores it blind against
Pruijn. **It does not flip anything, and its conclusion is that the term should
not be flipped on either.**

Driver: `experiments/intracellular_retardation.py`. Code changed: three additive
functions and one default-`None` keyword in `tmesim/molecular/two_route.py`, plus
ten tests in `tests/test_two_route.py`. Nothing is wired into `compile_drug.py`,
`DrugParams.diffusivity` is bit-identical, and `calibrate_alpha()` still returns
**148.60815765445545** (asserted by `--stage pruijn`, which refuses to score
otherwise).

> **The rule this file was written under, unchanged from its parents.** Nothing
> is fitted to Pruijn. `results/sweep/retard_prereg.json` was written to disk
> **before the term was scored even once**, and it registers *in advance* both the
> rule and the five specific ways the author expected it to fail. All five
> happened. The parent `results/sweep/tworoute_prereg.json` is untouched,
> `HBD_LOGP_DECREMENT` is still 1.0, and the parent's checking stage still passes.

---

## 0. Bottom line

**The physics is right, the composition is verified, the term is real — and the
best-sourced version of it makes the model decisively worse. Retardation does not
explain paclitaxel, and the anomaly stands.**

Five findings, in order of weight:

1. **The commissioning brief's proposed form is wrong physics, and this is
   demonstrable rather than arguable.** `D_eff = D_free/(1+K)` applied to the
   *inclusion conductivity* fails a direct transient solve by **22–35 %**, growing
   with `K` (§3.3). Reversible binding does not change a steady-state conductance
   at all; it changes **capacity**, which appears only in a transient. The correct
   composition is `D_app = σ_eff / s_eff` with `s_eff` an exact volume average,
   and it tracks the direct solve to **0.01 %**.

2. **It follows that paclitaxel's transcellular gain *cannot* fall — under any
   value of `K`, for reasons that have nothing to do with whether retardation is
   real.** The gain is a ratio of two *conductances*. The brief nominated a fall
   in the gain as the test of this physics and a failure to fall as a refutation.
   That test is not diagnostic and **it was registered as non-diagnostic before
   scoring** (`retard_prereg.json`, registered prediction 4). The diagnostic
   quantity is the apparent D, and it is reported instead.

3. **`K` cannot be sourced per compound — not for one of the 71 anchor compounds,
   and not even for paclitaxel.** The taxane route the brief suggested is blocked
   by an **UNVERIFIED** input: no primary mammalian intracellular tubulin
   concentration could be retrieved. It brackets paclitaxel's `K_b` at **9–1852**,
   which cannot discriminate anything (§2.4).

4. **Sourced as a STRUCTURAL rule, it scores worse on every measure that matters,
   and catastrophically on two.** R² **0.718 → 0.606**, geometric bias **1.03× →
   0.098×**, within-3× coverage **99 % → 0 %**, and — the sharpest one — the
   matched-MW within-tie rank agreement collapses from **77.4 % to 51.6 %, i.e. to
   a coin flip** (§4). The disclosure sweep puts the best multiplier at **0.01–0.03
   of the sourced value**, i.e. the anchor data prefer *essentially no retardation*.

5. **Pruijn's own paper says why, and it was found by reading the paper rather
   than by rationalising the score.** Their Dmcl is explicitly a *free-drug*
   coefficient, and compounds that did show sequestration had it removed into
   separately fitted reaction terms rather than left in `D` (§4.4). The
   measurement and the model disagree because **the measurement is not the
   quantity a retardation term predicts.**

**And it makes the blast radius worse, not better — over-determinedly so.** The
panel differential goes **1.60× → 10.48×** (§5): retardation moves paclitaxel
from 1.6× *faster* than doxorubicin to **10× slower**. It assigns an 8.4×
retardation to doxorubicin, the α calibration anchor, whose membrane is shut and
whose cells would need **325 days** to load that bound pool (§5). It drops every
drug's boundary layer to **0.1–0.3 grid cells** at the published `dx = 12 µm`
(§6). **2 of 12 verdicts move — and the new one is cisplatin PASS → FAIL**, a
passing verdict destroyed. And in `recal`, the mode the codebase actually runs,
**α collapses from the parent's 22–36 to 0.03–2.5 and goes NEGATIVE for
paclitaxel (−0.0131)** — not a worse calibration but a meaningless one (§6).

**Recommendation: §7. Short version — the flip should proceed, WITHOUT this term,
and §5.3 item 1 of the parent should be closed as ANSWERED rather than left
pending.**

---

## 1. What was asked, and what each answer rests on

| question | answer | rests on |
|---|---|---|
| Does Choi 2011 verify? | **Yes, cell line and signature — but the brief's numbers are not time-matched** | §2.1 |
| Where does `K` come from? | **Not per compound. A structural rule, secondary-sourced** | §2.2–§2.5 |
| Is the composition verified? | **Yes, against a direct transient solve; and the brief's form is falsified by it** | §3 |
| What does it do to Pruijn? | **Worse on every axis; rank agreement collapses to chance** | §4 |
| What happens to paclitaxel? | **Its gain cannot move. Its apparent D falls 10× below the anchor** | §5 |
| Blast radius? | §6 |
| Flip? | **Yes — without this term** | §7 |

---

## 2. Sourcing `K`

### 2.1 Choi 2011 verified — and the motivating numbers partly corrected

The brief stated: *"Choi 2011 (PMID 21225235) measured paclitaxel at 38 %
penetration through an HT-29 multicellular layer at 1 µM against 100 % at 10 µM."*
The abstract was retrieved and read. Verbatim:

> **Choi MS, Kim SH, Kuh HJ.** *Oncol Rep* 2011;25(3):863–870, PMID **21225235**,
> [DOI](https://doi.org/10.3892/or.2011.1138):
> *"The penetration rate measured after 24 h by cytotoxicity of the conditioned
> media was 40 and 38% in DLD-1 (20 µM) and HT-29 MCLs (1 µM), respectively"*
> and
> *"In HT-29 MCL, full penetration of PTX was obtained at 10 µM after 48 h,
> whereas only 80% was obtained at 1 µM."*

**The cell line is right and the saturable signature is real.** But the 38 % and
the 100 % are **not the same experiment**: the 38 % is at **24 h by the
cytotoxicity method**, the 100 % at 10 µM is at **48 h by [¹⁴C]-paclitaxel**. The
**time- and method-matched pair is 80 % at 1 µM against 100 % at 10 µM.**

**This weakens the motivating evidence and it is stated here rather than buried.**
The saturable restriction in the anchor cell line, matched, is about **1.25-fold**,
not the ~2.6-fold the brief's pairing implies. The qualitative claim — saturable
penetration in a P-gp-negative line, therefore binding rather than a pump —
survives intact. The magnitude that made it look like a large effect does not.
(The unambiguously saturable case in that paper is **DLD-1**, which is not the
anchor line: 10 % at 24 h and 20 % at 96 h at 1 µM.)

`results/two_route_efflux.md` §2.5 and §6 quote the 38 %/100 % pairing. **Those
passages overstate the effect and should be read against this correction.** That
file is not owned here and was not edited.

The P-gp-negative premise re-verified independently, verbatim:

> **Wils P et al.** *Biochem Pharmacol* 1994;48(7):1528–1530, PMID **7945455**,
> [DOI](https://doi.org/10.1016/0006-2952(94)90580-0) — *"No functional expression
> of P-glycoprotein was found in the differentiated HT29-18-C1 cell line."*

### 2.2 Per-compound `K` for the anchor set: not available

The Pruijn anchor is 67 benzotriazine di-N-oxide analogues of tirapazamine plus
NLCQ-1, urea, mannitol and acridan-9(10H)-one. **No intracellular bound fraction,
cell-to-medium partition or unbound fraction could be retrieved for any of them.**
A per-compound `K_b` is unavailable and the house rule forbids inventing one. The
brief's item 2 provides for exactly this case, so a **structural** rule was
registered instead.

### 2.3 The structural rule, and its one-line derivation

Kilford 2008 predicts the unbound fraction in a **hepatocyte** incubation from
lipophilicity alone:

```
fu_hep = 1 / (1 + 125 * VR * X),   X = 10^(0.072*logP² + 0.067*logP − 1.126)
```

with `VR` the volume ratio of **cells to incubation**. The bound-to-free ratio in
the incubation is `125·VR·X`. All that bound material sits **inside the cells**,
which occupy a fraction `VR` of the incubation, so per unit **cell volume**

```
K_b = (125 · VR · X) / VR = 125 · X
```

**The volume ratio cancels exactly.** That is what makes the number usable here:
the result is a property of a millilitre of cells, not of Kilford's assay
geometry, so it transfers to a cell sitting in a tissue. It uses `logP`, which the
engine already computes. This is a genuine structural rule from an existing
descriptor, registered before scoring.

**Provenance, stated exactly, because it is weaker than the parents':**

| | status |
|---|---|
| **Kilford PJ, Gertz M, Houston JB, Galetin A.** *Drug Metab Dispos* 2008;36(7):1194–1197, PMID **18411401**, [DOI](https://doi.org/10.1124/dmd.108.020834) | **PAYWALLED, NOT READ.** HTTP 403 at ASPET, no PMC record, Unpaywall `is_oa: false`. Its **abstract is verified** and confirms the substitution architecture in the authors' own words — *"further extended to allow direct prediction of fu(hep) from drug lipophilicity data by substituting the fu(mic) term by the Hallifax and Houston predictive equation"* — but carries **no coefficients**. |
| **Hallifax D, Houston JB.** *Drug Metab Dispos* 2006;34(4):724–726, PMID **16552024**, [DOI](https://doi.org/10.1124/dmd.105.007658) | **PAYWALLED, NOT READ.** Letter; PubMed records *"[Abstract not available]"*. |
| **Eur J Pharm Sci** 2018, PMID 30098391, PMC6365006, [DOI](https://doi.org/10.1016/j.ejps.2018.08.008) | **READ VERBATIM, open access.** Gives both equations; defines `C` as *"the microsomal protein concentration in mg/mL"* and `VR` as *"the ratio of cell volume to incubation volume and is 0.005 at a cell concentration of 10⁶ cells/mL"*. |
| **Yim DS, Bae SH, Choi S.** *Transl Clin Pharmacol* 2021;29(2), PMID 34235120, PMC8255549, [DOI](https://doi.org/10.12793/tcp.2021.29.e12) | **READ VERBATIM, open access.** Same equation, same four constants, same `VR` definition and value. |

**So the four constants (125, 0.072, 0.067, −1.126) are SECONDARY-SOURCED from two
independent witnesses that agree exactly, not read from either primary.** That is
a weaker provenance than the Abraham coefficients the parent used, and it is
labelled as such rather than dressed up.

**Two things that are interpretation, not source, and are recorded as such:**

* The reading of `125` as *mg protein-equivalent per mL of cell volume* follows
  from the algebra of the substitution (`125·VR` occupies exactly the slot `C` in
  mg/mL occupies, and `VR` is dimensionless) and nothing retrieved contradicts it.
  But **no source states the authors intended it as a physical protein
  concentration**; Kilford's own abstract calls it *"A new nonlinear empirical
  equation"* fitted on 39 drugs. It is dimensionally interpretable that way, not
  measured that way.
* **A correction to the pre-registration, recorded rather than silently fixed.**
  `retard_prereg.json` as written said the rule uses *"logP for neutrals and bases
  and logD7.4 for acids"*. **That is wrong.** As published it is **logP for BASES,
  logD7.4 for ACIDS AND NEUTRALS** (Yim 2021, verbatim: *"In both equations, logD
  is used for acidic or neutral drugs and logP for basic drugs"*; PMC6365006
  attributes the convention to Austin 2002). The original text is left in place in
  the pre-registration and corrected in its `_corrections_after_writing` block, so
  the audit trail shows what was believed when. **No number moves**: the engine
  holds one lipophilicity per compound and could not have applied the distinction
  either way.

### 2.4 The rule's own domain says it does not cover the compounds that matter

| logP | `K_b` | bound % | `R = 1 + φ·K_b` |
|---:|---:|---:|---:|
| −2.62 | 19.48 | 95.1 % | 16.37 |
| −0.80 | 9.19 | 90.2 % | 8.25 |
| **−0.47** (vertex) | **9.02** | 90.0 % | 8.12 |
| 0.00 | 9.35 | 90.3 % | 8.38 |
| 2.00 | 24.71 | 96.1 % | 20.50 |
| 2.97 | 63.83 | 98.5 % | 51.37 |
| **3.74** (paclitaxel) | **169.28** | 99.4 % | **134.59** |

**Read the logP = 0 row.** The rule asserts that a compound with no lipophilicity
at all is **90 % bound** inside a cell. **Urea and mannitol are in the Pruijn set,
are classical non-binding flux markers, and Pruijn uses urea as its internal
standard.** The rule is also **non-monotone** — its vertex is at logP = −0.465, so
it predicts *more* binding for the most hydrophilic compounds in the set than for
neutral ones, an artefact of extrapolating a quadratic fitted over drug-like logP
down to −2.62. Both facts are registered in `retard_prereg.json` as expected
out-of-domain failures **before** the rule was scored.

Kilford's own abstract, verbatim, agrees: prediction accuracy is *"high for
relatively hydrophilic drugs (logP/D ≤ 2.5)"* and *"less accurate ... for
lipophilic drugs (logP > 3)"*. **Paclitaxel (3.74) and docetaxel (3.26) are both
outside the range the source vouches for.**

### 2.5 The per-compound taxane route: bracketed, not pinned

The brief named tubulin stoichiometry as a candidate. At tracer concentrations
`K_b = [T]/Kd`.

**Verified:** stoichiometry *"0.99 ± 0.04 ligand per alpha beta dimer"* (Díaz &
Andreu, *Biochemistry* 1993;32(11):2747–2755, PMID **8096151**,
[DOI](https://doi.org/10.1021/bi00062a003)). Polymer `Kd` *"8.7 x 10-7 M"*
(Parness & Horwitz, *J Cell Biol* 1981;91(2):479–487, PMID **6118377**, PMC2111958,
[DOI](https://doi.org/10.1083/jcb.91.2.479), **full text read**). A second primary
gives `Kb = 3.7 x 10^7 M^-1` ≈ 27 nM (Andreu & Barasoain, *Biochemistry*
2001;40(40):11975–11984, PMID **11580273**) — a **~30-fold disagreement, not
averaged**.

> *Housekeeping, because it would propagate:* the **PubMed abstract string for
> PMID 6118377 contains a typo, "8.7 x 19(-7) M"**. The PMC full text reads 10⁻⁷.

**UNVERIFIED, and it is the blocking input: the intracellular tubulin
concentration.** No primary mammalian µM or mg/mL value could be retrieved. The
only primary number is a **percentage** — *"approximately 2.5--3.3% of their total
protein"* (Hiller & Weber, *Cell* 1978;14(4):795–804, PMID **688394**). The
ubiquitous "~20 µM" is a secondary restatement of that paper and does not
reconstruct cleanly from its own percentage. A published honest bound is
*"anywhere between 8 µM and 50 µM"* (eLife 2022, PMID 35029146).

| Kd source | `K_b` bracket at [T] = 8–50 µM |
|---|---:|
| Parness & Horwitz 1981 (0.87 µM) | 9.2 – 57.5 |
| Andreu & Barasoain 2001 (0.027 µM) | 296 – 1852 |
| **structural rule, for comparison** | **169.3** |

The bracket spans **9 to 1852** and contains the rule's value. That is a
consistency check and nothing more — a range covering two orders of magnitude
discriminates nothing, and it exists only for taxanes, so it applies to **zero**
compounds in the anchor set. Also note: **a Kd for unpolymerised tubulin dimer
does not exist** — *"binding to the unassembled protein is not detected within
ligand solubility limits"* (Díaz, Menéndez & Andreu, PMID **8104479**), and
taxane-site ligands *"show a loss of affinity of at least four orders of
magnitude"* on unassembled tubulin (Prota et al., eLife 2023, PMID **36876916**).

**So the honest answer to "can `K` be sourced per compound?" is NO — not for the
anchor set, and not even for the one drug the whole question is about.**

---

## 3. The composition, verified against a direct transient solve

### 3.1 Why the brief's form is wrong

The brief specified `D_eff = D_free/(1+K)` applied to the inclusion conductivity.
**Reversible binding at local equilibrium does not change a steady-state
conductance at all.** At steady state the bound pool is time-invariant, carries no
flux, and drops out of the flux balance; only the free profile drives transport,
and it obeys the *unretarded* conductivity. Binding changes **capacity**, and
capacity appears only in a **transient**.

Writing the homogenised problem in matrix-equivalent potential `u`, two
coefficients appear and they homogenise differently:

```
conductance  σ_eff : matrix σ_m, inclusion K·D_cyto
                     -> Maxwell / Hasselman & Johnson, ALREADY here, UNAFFECTED
capacity     s_eff = (1−φ)·1 + φ·K·(1+K_b)     -> an EXACT volume average
apparent D         = σ_eff / s_eff
```

The minimal additive form adopted is the **ratio** of the capacity with binding to
the capacity without, which isolates the new physics and does not relitigate the
capacity convention the module already ships against:

```
R = s_eff(K_b) / s_eff(0) = 1 + φ·K·K_b / [(1−φ) + φ·K]
```

which at the pre-registered `K = 1` is exactly `R = 1 + φ·K_b`. `K_b = 0` gives
`R = 1` exactly, so the term is **bitwise off by default**.

### 3.2 The instrument, and the identity

`--stage algebra` builds a 1-D finite-volume laminate
`[matrix | membrane | cell | membrane] × N` carrying the **real 4 nm bilayer as its
own control volume**, and runs it as a transient from an empty start. It knows
nothing about Maxwell or about retardation factors. Daynes–Barrer: the cumulative
flux tends to `J_ss·(t − τ)`; **`J_ss` tests the conductance and `τ` tests the
capacity**, so the leg this work adds is tested separately from the leg already
verified.

| N periods | `K_b` | `J` meas/pred | **`τ` meas/pred** | `D_app` meas/pred |
|---:|---:|---:|---:|---:|
| 10 | 30 | 1.000000 | 1.004694 | 0.995328 |
| 20 | 30 | 1.000000 | 1.001173 | 0.998828 |
| 40 | 30 | 1.000000 | 1.000293 | 0.999707 |
| **80** | **30** | **1.000000** | **1.000073** | **0.999927** |

The steady flux matches to machine precision at every N. The **time lag converges
onto the volume-averaged capacity** as the laminate is refined — the residual is
the laminate's own finite-size correction and falls with N as it must.

`R = s_eff(K_b)/s_eff(0)` holds as an identity to **1.86e-16** over 72
`(φ, K, K_b)` points (asserted again in
`tests/test_two_route.py::TestIntracellularRetardation`).

### 3.3 The discriminating test — and the brief's form fails it

| `K_b` | `D_app` direct | **capacity form** | err | **brief's form** | err |
|---:|---:|---:|---:|---:|---:|
| 0 | 148.1641 | 148.1696 | 0.00 % | 148.1696 | 0.00 % |
| 3 | 44.0082 | 44.0110 | **0.01 %** | 53.8148 | **22.28 %** |
| 30 | 6.0065 | 6.0069 | **0.01 %** | 7.9948 | **33.10 %** |
| 300 | 0.6234 | 0.6234 | **0.01 %** | 0.8403 | **34.79 %** |

**The capacity form tracks the direct solve; the brief's form does not, and it
fails in a way that grows with `K_b`** — because it applies a capacity effect to a
conductance, and because at φ = 0.79 it entirely misses the retardation of
*matrix*-phase transport by cell storage, which is the dominant term.

### 3.4 Where local equilibrium fails, measured rather than asserted

The long-time asymptote is correct at **every finite `P`** — a slow membrane
delays equilibrium, it does not abolish it — so the composition is not wrong in
the shut-membrane limit. What is wrong is applying it to an experiment **shorter
than the cell's loading time** `τ_cell = (a/3)·K·(1+K_b)/P`:

| `P` (µm/s) | `τ_cell` (s) | note |
|---:|---:|---|
| 1e−4 | 8.1e5 | doxorubicin-like: **9 days** at `K_b`=30 |
| 1.0 | 81 | |
| 1e4 | 0.008 | |

**For the shut-membrane drugs the retardation factor is an UPPER BOUND, not an
estimate.** `two_route.equilibration_time_s` exposes this, and §5 shows it biting
on the α anchor.

---

## 4. Scored against Pruijn — blind, once

`--stage pruijn` asserts the live constants against `retard_prereg.json`, asserts
default-off is bitwise identical, and asserts `calibrate_alpha() ==
148.60815765445545` **before** it will score anything.

### 4.1 The score, at the sourced decrement 1.34

| arm | spread | **R²** | Spearman | geo bias | within 3× |
|---|---:|---:|---:|---:|---:|
| **A. retardation OFF (control)** | 8.3× | **0.718** | +0.821 | **1.03×** | **99 %** |
| C. structural rule × 0.01 | 8.3× | 0.729 | +0.825 | 0.936× | 97 % |
| C. structural rule × 0.03 | 8.2× | **0.731** | +0.822 | 0.792× | 96 % |
| C. structural rule × 0.1 | 9.0× | 0.702 | +0.808 | 0.521× | 82 % |
| C. structural rule × 0.3 | 10.1× | 0.648 | +0.783 | 0.265× | 34 % |
| **C. structural rule × 1 — THE SOURCED VALUE** | 11.2× | **0.606** | +0.732 | **0.098×** | **0 %** |
| C. structural rule × 3 | 11.6× | 0.590 | +0.720 | 0.035× | 0 % |
| D. CONTROL: uniform `K_b` = 1 | 8.3× | 0.718 | +0.821 | 0.577× | 96 % |
| D. CONTROL: uniform `K_b` = 10 | 8.3× | 0.718 | +0.821 | 0.116× | 1 % |
| D. CONTROL: uniform `K_b` = 100 | 8.3× | 0.718 | +0.821 | 0.013× | 0 % |

At the pre-registered decrement 1.0 the same pattern holds: **0.725 → 0.584**,
best multiplier 0.01.

**Read control arm D first.** A compound-independent `K_b` is a rigid shift in log
space, so it **cannot** change R² or Spearman — only bias. That is exactly what
the table shows, at three magnitudes. It was registered in advance as the expected
result, and it means **every** change in R² in arm C is attributable to the logP
*dependence* of the structural rule and to nothing else. If arm D had moved R²,
there would be a bug.

**At the sourced value every one of the 68 scored compounds is under-predicted**
(worst ratio 0.0149, best 0.293 — the whole distribution is below 1/3), which is
why coverage is 0 %.

**The disclosure that matters:** the multiplier that maximises R² is **0.01–0.03**,
i.e. the anchor data prefer a binding capacity **30–100× smaller than the sourced
one** — which is to say, essentially none. The adopted value is emphatically not
the sweep optimum, and the gap is not marginal.

### 4.2 The matched-MW test — the sharpest refutation

This is the test the two-route model earned its keep on: where a MW-only model is
flat by construction, does the model reproduce the measured within-tie ordering?

| | within-tie rank agreement |
|---|---:|
| **retardation OFF** | **77.4 %** (24 concordant / 7 discordant) |
| **structural rule** | **51.6 %** (16 concordant / 15 discordant) |

**Retardation destroys the matched-MW rank agreement, taking it to a coin flip.**
That is a cleaner refutation than R², because it is exactly the axis on which the
transcellular route was originally justified, and it is immune to any uniform
scaling.

| MW | measured spread | model OFF | model ON |
|---:|---:|---:|---:|
| 192 | 1.4× | 1.03× | 1.00× |
| 213 | 2.5× | 1.10× | 1.01× |
| 263 | 3.7× | 1.22× | 2.44× |
| 289 | 1.4× | 1.68× | 1.58× |
| 305 | 3.8× | 2.69× | 2.24× |

### 4.3 Every registered failure mode happened

| registered before scoring | outcome |
|---|---|
| 1. R² expected to **fall** | **Yes.** 0.718 → 0.606. |
| 2. Rule expected **out of domain on the flux markers** | **Yes.** Urea/mannitol carry `K_b` ≈ 9–19. |
| 3. Geometric bias expected to **collapse** | **Yes.** 1.03× → 0.098×; coverage 99 % → 0 %. |
| 4. Transcellular gain **cannot move** | **Yes.** Invariant at every `K_b` (§5). |
| 5. Local equilibrium expected to **fail for shut membranes** | **Yes.** τ_cell = 325 days for doxorubicin (§5). |

**None of these is a discovery of this run. All five were written down first.**

### 4.4 Why the measurement disagrees — from Pruijn's own paper, not from the residuals

The Pruijn PDF was retrieved (publisher's own OA file, via the Internet Archive
snapshot already recorded in `experiments/diffusivity_falsification.py`) and its
Methods and Discussion read. Two passages are load-bearing and neither was known
to the parent files.

**Methods, verbatim:**

> *"The diffusion coefficient of each compound was determined by fitting the
> concentration–time data for both the donor and receiver compartment to Fick's
> second law, **including reaction terms in the MCL when appropriate**, as
> previously described."* [ref 36 = Hicks KO et al., *Br J Cancer* 1997;76:894]

**Results, verbatim:**

> *"Under these conditions the flux curves were well modelled as Fickian
> diffusion, generally with little or no metabolic consumption. For subsets of
> compounds with DNA intercalating functionality (e.g., SN26955) ... **the mass
> balance showed significant loss during the flux experiments, which was
> accommodated by addition of reaction terms in the MCL.**"*

**Discussion, verbatim — the sentence that decides it:**

> *"The second caution is that **Eqn 1 applies to free (unbound) drug only.** If
> high affinity drug binding sites are present in high concentrations (e.g., DNA
> interbase-pair sites in the case of DNA intercalators), flux can be greatly
> slowed during the extended pre-steady-state phase."*

**Read together: Pruijn's `Dmcl` is intended and fitted as a FREE-drug diffusion
coefficient.** Sequestration, where it was large enough to show in the mass
balance, was **taken out into separately fitted reaction terms** rather than left
inside `D`. So a retardation term predicts a quantity **Pruijn deliberately did
not report**.

This is not a rescue of a bad score — it is an independent, documentary
explanation that **agrees with the score**: the data prefer a multiplier of
0.01–0.03, i.e. **no retardation**, which is precisely what a free-drug
coefficient should show. Two independent lines of evidence, one documentary and
one numerical, reach the same conclusion.

> **The honest counter-argument, stated because it is real.** Linear reversible
> binding at local equilibrium does *not* break Fickian behaviour — it merely
> rescales `D` — so "the flux curves were well modelled as Fickian" does **not**
> by itself exclude retardation being absorbed into the fitted `D`. The
> case rests on the authors' explicit framing plus their removal of sequestration
> into reaction terms, not on the Fickian statement. It is strong but it is not a
> proof, and `results/two_route_transport.md` §6's "one cell line, one method"
> caveat still binds.

---

## 5. The panel, and what actually happens to paclitaxel

At the sourced decrement 1.34, with the structural rule on:

| drug | logP | `K_b` | `R` | **gain** | scale (conductance) | **apparent scale** | vs dox | `τ_cell` (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| doxorubicin | 0.00 | 9.35 | 8.38 | 1.0000 | 0.1512 | 0.01803 | 1.0000 | **2.8e7** |
| **paclitaxel** | 3.74 | 168.34 | 133.85 | **1.5985** | 0.2416 | **0.00181** | **0.1001** | 205.7 |
| docetaxel | 3.26 | 90.06 | 72.07 | 1.0102 | 0.1527 | 0.00212 | 0.1175 | 7108 |
| epirubicin | 0.00 | 9.35 | 8.38 | 1.0000 | 0.1512 | 0.01803 | 1.0000 | 2.8e7 |
| 5-fluorouracil | −0.80 | 9.19 | 8.25 | 1.0083 | 0.1524 | 0.01847 | 1.0241 | 472.3 |
| cisplatin | 0.20 | 9.70 | 8.65 | 1.0832 | 0.1637 | 0.01892 | 1.0491 | 66.4 |
| **panel differential** | | | | **1.599×** | **1.599×** | **10.480×** | | |

**Three things, and each is bad for the term.**

1. **The `gain` column is bit-identical to the un-retarded model's, for every drug,
   at every `K_b`.** This is asserted as a test
   (`test_retardation_cannot_move_a_conductance`) and in the module self-test. The
   brief's nominated test — "paclitaxel's gain should fall; if it does not, that
   is a refutation" — **cannot be passed by any value of `K_b`**, and the reason is
   arithmetic, not physical. It was registered as non-diagnostic before scoring.
   **The refutation is real, but it is not this test that delivers it; §4.1 and
   §4.2 are.**

2. **The blast radius gets six times worse, not better.** The panel differential
   goes **1.599× → 10.480×**. Retardation does not tame paclitaxel — it **reverses
   its sign**, taking it from 1.60× *faster* than doxorubicin to **10× slower**.
   The parent's central argument for flipping was that sourcing the HBD decrement
   shrank the differential 4.94× → 1.60×. **This term throws that away.**

3. **It does unphysical things to the α calibration anchor.** Doxorubicin's
   membrane is shut, so its `τ_cell` is **2.8e7 s ≈ 325 days**. It cannot load an
   8.4× bound pool in any experiment, in any tissue, on any relevant timescale —
   yet the term applies the full retardation to it, and doxorubicin is the drug α
   is calibrated on. **This is the same failure mode that condemned the
   parameter-free HBD variant in `results/two_route_efflux.md` §4.1: wrong physics
   on the anchor drug.**

---

## 6. Blast radius

`--stage blast`, published `dx = 12 µm` breast grid, 458×628, 167,780 cells /
89,785 tumour, same `dt`, same 6 h window, same seed, same `PDConfig`, same p99
bath-matching convention, same frozen literature bands, both α arms. α is the
published **148.6082** in `frozen` and re-solved in `recal`.

**Scales, retardation ON at hbd_decrement 1.34:**

| drug | scale | fold | vs dox |
|---|---:|---:|---:|
| doxorubicin | 0.018034 | 55.5× | 1.0000× |
| paclitaxel | 0.001805 | 554.0× | 0.1001× |
| docetaxel | 0.002119 | 472.0× | 0.1175× |
| epirubicin | 0.018034 | 55.5× | 1.0000× |
| 5-fluorouracil | 0.018468 | 54.1× | 1.0241× |
| cisplatin | 0.018919 | 52.9× | 1.0491× |

**The magnitude alone should stop this.** The un-retarded correction was a 6.6×
reduction in `D` (scale 0.151); this is a **55× to 554× reduction**. Two
consequences are visible before any verdict is read:

* **Every drug's boundary layer collapses far below one grid cell** at
  `dx = 12 µm`. `results/two_route_efflux.md` §4.3 already flagged docetaxel at
  `L/dx = 0.94` as a resolution hazard at a 6.6× reduction. Scaling its published
  `frozen`-arm figures by `L ∝ √D`:

  | drug | `L/dx` at the flip config (scale) | **`L/dx` with retardation ON** |
  |---|---:|---:|
  | doxorubicin | 0.93 (0.1512) | **0.32** (0.018034) |
  | paclitaxel | 1.18 (0.2416) | **0.10** (0.001805) |
  | docetaxel | 0.94 (0.1527) | **0.11** (0.002119) |

  At 55–554× the `dx = 12 µm` grid does not resolve the gradient for **any**
  drug — not even the α anchor — so the `frozen` numbers below are quoted on a
  grid that cannot support them. This is a much stronger version of the parent's
  §5.3 item 4 caveat, and on its own it disqualifies shipping this term at this
  grid.
* **`recal` α must move a long way** to absorb a 55× change in `D`, which is what
  `recal` is for; `frozen` cannot, so `frozen` ratios fall through the floor.

**Solved.** Four rows are published **NOT TESTABLE** (paclitaxel and
5-fluorouracil, both arms) — single-study bands — leaving 8 scored verdicts.

| mode | drug | scale | ratio ship | ratio corr | verdict → verdict | α | moved? |
|---|---|---:|---:|---:|---|---:|:--:|
| frozen | doxorubicin | 0.018034 | 0.2301 | 0.0681 | PASS → PASS | 148.6082 | no |
| frozen | paclitaxel | 0.001805 | 0.1503 | 0.0211 | NOT TESTABLE | 148.6082 | no |
| **frozen** | **docetaxel** | 0.002119 | 0.1484 | **0.0210** | **FAIL → PASS** | 148.6082 | **YES** |
| frozen | epirubicin | 0.018034 | 0.2301 | 0.0681 | PASS → PASS | 148.6082 | no |
| frozen | 5-fluorouracil | 0.018468 | 0.3379 | 0.1620 | NOT TESTABLE | 148.6082 | no |
| **frozen** | **cisplatin** | 0.018919 | 0.3085 | **0.1312** | **PASS → FAIL** | 148.6082 | **YES** |
| recal | doxorubicin | 0.018034 | 0.2301 | 0.2395 | PASS → PASS | **2.4031** | no |
| recal | paclitaxel | 0.001805 | 0.1503 | 0.1520 | NOT TESTABLE | **−0.0131** | no |
| recal | docetaxel | 0.002119 | 0.1484 | 0.1595 | FAIL → FAIL | **0.0335** | no |
| recal | epirubicin | 0.018034 | 0.2301 | 0.2395 | PASS → PASS | **2.4031** | no |
| recal | 5-fluorouracil | 0.018468 | 0.3379 | 0.3091 | NOT TESTABLE | **2.4677** | no |
| recal | cisplatin | 0.018919 | 0.3085 | 0.2887 | PASS → PASS | **2.5349** | no |

**2 of 12 verdicts move, both in `frozen`. In `recal`, 0 of 12 move.**

**But the verdict count is the least of it, and reading only that column would
miss the finding.** Two things in this table are worse than a moved verdict.

**(i) `recal` α collapses, and for paclitaxel it goes NEGATIVE.** The parent's
`recal` α values at the flip config span **22.22 to 35.69** and it reports that
tightening as evidence the panel is approaching a single α. With retardation on:

| | parent, flip config | **this term** |
|---|---:|---:|
| doxorubicin / epirubicin | 22.2234 | **2.4031** |
| cisplatin | 24.0953 | 2.5349 |
| 5-fluorouracil | 22.4093 | 2.4677 |
| docetaxel | 22.4521 | **0.0335** |
| **paclitaxel** | 35.6934 | **−0.0131** |

**A negative α is not a worse fit, it is a meaningless one.** α is a
non-negative transport-calibration constant; the solver returning −0.0131 means
the retardation-scaled `D` has pushed the calibration outside the domain in which
it has a solution at all. Docetaxel's 0.0335 is three orders of magnitude below
the anchor's and is the same failure in a less obvious form. **The mode the
codebase actually runs does not survive this term** — it reports "0 verdicts
moved" only because it has already broken.

**(ii) The one `frozen` movement in the parent becomes two here, and the new one
goes the wrong way.** The parent's sole movement was docetaxel FAIL → PASS. That
persists, but **cisplatin also moves, PASS → FAIL** — a previously passing,
testable verdict destroyed. So this term does not merely fail to help; it costs a
verdict that the flip config keeps.

**Taken with §5 and the `L/dx` table above, the blast radius argument against the
term is over-determined:** it makes the panel differential 6.5× worse, it applies
an unphysical retardation to the α anchor, it drops every drug's boundary layer to
0.1–0.3 grid cells, it breaks a passing verdict, and it drives the calibration
constant negative.

---

## 7. Recommendation

**The decision is the reader's. This section produces the evidence and a
recommendation, and flips nothing.**

### 7.1 What I recommend

**Recommendation: FLIP THE TWO-ROUTE MODEL AS `results/two_route_efflux.md` §5.1
SPECIFIES — and DO NOT include this retardation term.**

| | value | why |
|---|---|---|
| `two_route` | **ON** | parent §5.2, unchanged by this work |
| `hbd_decrement` | **1.34**, re-registered | parent §3 |
| efflux | **OFF** | parent §2 |
| **intracellular retardation** | **OFF** | **this file** |
| α mode | **`recal`** | parent §5.2 |

**The parent made the flip conditional on this term (§5.3 item 1: *"Add the `K > 1`
retardation basis first, or flip with paclitaxel's row explicitly carrying that
caveat"*). That condition should now be discharged as ANSWERED, not as pending.**
The term was built, its composition was verified to 0.01 % against a direct
transient solve, and the best-sourced version of it makes the model worse on R²,
on bias, on coverage, and — decisively — on the matched-MW rank agreement that
justified the transcellular route in the first place. The blocker is not
outstanding; it has been investigated and closed negative.

**Paclitaxel's row should carry the caveat**, which is the parent's own stated
alternative. The caveat is now better specified than it was: paclitaxel's anomaly
is **not** explained by intracellular retardation, and §4.4 supplies a documented
reason why the anchor dataset cannot be used to argue otherwise.

### 7.2 What this work changes about the flip decision

| parent §5.3 reason not to flip | status after this work |
|---|---|
| 1. Model known-wrong for taxanes; add the `K > 1` term first | **ANSWERED NEGATIVE.** The term is built and verified; the sourced version scores materially worse, makes the blast radius 6.5× larger, destroys a passing verdict, and drives `recal` α negative. It should not be added. The motivating measurement is also weaker than quoted (§2.1). |
| 2. One cell line, still | **OPEN, unchanged**, and now the most valuable outstanding item by a distance. |
| 3. Descriptor layer wrong by ~2.4 log units for cisplatin | **OPEN, unchanged, owned elsewhere.** |
| 4. Taxane grid-resolution hazard | **OPEN, unchanged** for the flip config; §6 shows it would become disqualifying if this term were included. |

### 7.3 What I recommend *against*

* **Against shipping the retardation term at any multiplier**, including the
  0.01–0.03 that maximises R². That multiplier is **fitted to Pruijn** and adopting
  it would do exactly what every pre-registration in this line of work exists to
  prevent. Its improvement (0.718 → 0.731) is within noise and it is worth
  nothing set against the cost of a fitted constant in the one clean anchor.
* **Against using the Choi 2011 38 %/100 % pairing** in any further argument. Use
  the time-matched 80 %/100 % (§2.1).
* **Against treating "paclitaxel's gain did not fall" as the refutation.** It is
  arithmetically incapable of falling. The refutation is §4.1 and §4.2.
* **Against deleting the code.** The term is correct physics, verified, additive
  and default-off. It becomes relevant the moment this project scores against a
  measurement that reports an **apparent** rather than a free-drug diffusivity —
  and `results/two_route_efflux.md` §2.2 (iii) notes that Pruijn's authors invoke
  exactly such a transient/capacity effect when discussing efflux. The instrument
  should stay; only its application to this anchor is refused.

### 7.4 Follow-ups this produces, in value order

1. **Score against a second cell line.** Unchanged from both parents, and now
   unambiguously the top item.
2. **Find a transport dataset that reports an APPARENT diffusivity.** That is the
   measurement this term predicts, and it is the only way to test it. Pruijn
   cannot.
3. **Retrieve Hicks KO et al., *Br J Cancer* 1997;76:894** — Pruijn's "as
   previously described" MCL model. It would settle §4.4's remaining ambiguity
   definitively by showing whether the fitted model carries an equilibrium binding
   term or only irreversible reaction terms.
4. **Retrieve the two DMD primaries** (Kilford 2008, Hallifax & Houston 2006) to
   upgrade the four constants from secondary- to primary-sourced, if the rule is
   ever wanted for another purpose.
5. **A primary intracellular tubulin concentration.** It is the blocking input for
   any per-compound taxane binding calculation, and the widely cited "~20 µM" does
   not reconcile with its own source (§2.5).

---

## 8. What this licenses saying, and what it does not

**Established:**

* Reversible intracellular binding retards the **transient apparent** diffusivity
  and leaves the **steady-state conductance** exactly unchanged. The composition
  `D_app = σ_eff/s_eff` matches a direct transient PDE solve to 0.01 %; the form
  the commissioning brief proposed fails the same test by 22–35 %.
* `R = 1 + φ·K·K_b/[(1−φ)+φ·K]` is an exact capacity ratio, verified to 1.86e-16.
* `K` **cannot** be sourced per compound for the anchor set, nor for paclitaxel.
* Sourced as a structural rule, retardation makes the blind Pruijn score worse on
  R², Spearman, bias and coverage, and collapses matched-MW rank agreement to
  chance. The anchor data prefer a capacity 30–100× smaller than the sourced one.
* Pruijn's own Methods and Discussion say their `Dmcl` is a **free-drug**
  coefficient with sequestration removed into separate reaction terms.
* Choi 2011's HT-29 saturation, time- and method-matched, is ~1.25-fold, not the
  ~2.6-fold the brief's pairing implies.

**Not established, and must not be claimed:**

* **That intracellular retardation does not occur.** It certainly does. What is
  refuted is that it is measurable in, or applicable to, the Pruijn anchor, and
  that the Kilford-derived structural rule is the right magnitude for it.
* **That the structural rule is wrong in its own domain.** It is a published,
  twice-independently-quoted DMPK tool. What is shown is that it is **out of
  domain** here — on non-binding flux markers, and on the two taxanes its own
  source declines to vouch for.
* **That paclitaxel is explained.** It is not. `results/two_route_efflux.md` §2.5
  relocated the diagnosis from efflux to binding; this file shows binding does not
  carry it either, on the evidence available. **Paclitaxel remains the model's most
  questionable row, and its band is single-study and published NOT TESTABLE.**
* **That the four Kilford/Hallifax constants are primary-sourced.** They are
  secondary-sourced from two independent open-access witnesses. Neither primary
  was read.
* **That `125` is a measured intracellular protein concentration.** It is an
  empirical fitted constant that is dimensionally interpretable that way.

**Known missing physics, updated:**

1. **The `K > 1` retardation basis is no longer missing.** It is present, verified,
   default-off, and shown not to apply to the anchor. `results/two_route_efflux.md`
   §6 item 1 should be re-read against this file.
2. **A capacity-commensurable transport dataset is now the missing item**, in the
   sense that this project can model apparent diffusivity but has nothing to score
   it against.
3. **A BCRP/ABCG2 term** remains unmodelled (parent §6 item 2), unchanged.
4. Everything else in `results/two_route_transport.md` §6 and
   `results/cellular_exclusion.md` §7 carries over unchanged.

---

## 9. Reproducing this

```
./.venv/bin/python experiments/intracellular_retardation.py --stage sources   # <1 s
./.venv/bin/python experiments/intracellular_retardation.py --stage algebra   # ~25 s
./.venv/bin/python experiments/intracellular_retardation.py --stage pruijn    # ~30 s
./.venv/bin/python experiments/intracellular_retardation.py --stage panel     # ~5 s
./.venv/bin/python experiments/intracellular_retardation.py --stage blast     # long
./.venv/bin/python -m tmesim.molecular.two_route                              # self-test
./.venv/bin/python -m pytest tests/test_two_route.py -q
```
