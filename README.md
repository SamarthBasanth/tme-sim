# Vivum

**A mechanistic drug-transport simulator that runs on real, measured tumour tissue.**

Vivum reconstructs the microenvironment of a real tumour section from spatial
transcriptomics — cell types and states, vasculature, ECM, oxygen — then solves
convection–diffusion–reaction transport, interstitial fluid pressure and a
steady-state oxygen field on that geometry, and maps the concentration each cell
actually sees onto a Hill dose–response. Two tissues ship as flagships: a
glioblastoma section (**40,887 measured cells**) and a breast carcinoma section
(**167,780 measured cells**).

**Who it is for.** Pharma, biotech and academic *research* users asking a
transport question: on this tissue architecture, which of these compounds
actually reaches the cells, and is the limit delivery or potency? A compound that
never arrives and a compound that arrives too weak look identical in a survival
curve and need opposite fixes; Vivum separates them.

**The single output this project will defend is the *relative* comparison of
compounds on *one* tissue section.** Everything absolute has been tested against
external measurement and essentially all of it has failed. Both halves of that
sentence are documented below, with a pointer to the measurement behind each.

> ### Read this before you rank a compound series
>
> The shipped transport model derives a compound's tissue diffusivity from its
> **molecular weight, its logP and its hydrogen-bond donor count**, through a
> two-route law — paracellular between cells, transcellular through their
> membranes. Scored **blind and unfitted** against 71 compounds measured by one
> method in one cell line (Pruijn 2008, doi:10.1071/CH08240) it reaches
> **R² = 0.718, Spearman +0.821**, with **99 %** of compounds within 3× of
> measurement, a geometric bias of **1.03×**, and **24 of 31** exact-mass pairs
> ordered correctly.
>
> **That is a change of claim, and the previous one is worth stating because it
> was worse.** Until 2026-08-13 diffusivity was a function of **molecular weight
> alone** and that model was **FALSIFIED**: R² = 0.022 on those same compounds,
> and flat **by construction — exactly 1.00×** — at identical molecular weight
> where the measurement spans **3.7×** (four compounds at 263 Da) and **3.8×**
> (four at 305 Da), so it ordered **0 of 31** mass-matched pairs.
> → [`results/diffusivity_falsification.md`](results/diffusivity_falsification.md)
>
> **What is still not true, and it is narrower than "fixed".** The model
> reproduces **8.3×** of the measured **29.7×** spread — about **a quarter** — so
> it *under-separates* a congeneric series systematically rather than separating
> it correctly. It has been scored against exactly **one cell line**. And it
> bought nothing for absolute numbers: the same drug's compiled `D_eff` misses
> its measured multicellular-layer band by 1.5× (49.3 µm²/s against a measured
> 74–130), having previously missed it by 2.2× in the other direction.
>
> **So: read the order over a similar-mass series as a hypothesis worth testing,
> not as a measurement.** Ranking compounds that differ substantially in mass
> *and* in the properties below remains the best-supported case, and no absolute
> transport number is defensible in either.
> → [`results/two_route_promotion.md`](results/two_route_promotion.md)

> ### What Vivum does not do
>
> * It does **not** predict clinical response or patient outcome. No open dataset
>   anywhere pairs measured spatial drug distribution with measured treatment
>   response in the same tumour, and there is no glioma + temozolomide spatial
>   drug map at all — a 39-agent, 164-dataset survey established this.
>   → [`results/DATA_HUNT.md`](results/DATA_HUNT.md)
> * It is **not** a diagnostic, treatment-planning or surgical tool, and it is not
>   a "digital twin" or personalised-medicine product.
> * **Absolute kill fraction is not a defensible output**, and neither is any other
>   absolute quantity the engine produces. See *What has failed* below.
> * The spatial-attribution statistic is **refused in code** on spot data, where
>   immune positions are synthesised rather than measured.

---

## Run it

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e .
pip install rdkit openpyxl statsmodels xlrd pyxlsb rdata   # see docs/REPRODUCE.md §1
```

`pip install -e .` registers a `vivum` console script; `vivum <cmd>` and
`python -m tmesim.cli <cmd>` are the same tool. The second `pip` line is not in
`pyproject.toml` and is needed for the molecular descriptors and the
external-validation readers. Python ≥ 3.9; pin with
`pip install -c constraints.txt -e .` if you intend to reproduce a number rather
than merely run the code.

**If a file under `results/` tells you to run `codon ...`, run it as written.**
The product was named Vivum, renamed to Codon on 2026-08-14, and renamed back to
Vivum on 2026-08-22. The dated files written in between quote their reproduce
instructions as `codon ...`; they are a published record and are not edited, so
`codon` is still registered as an alias console script pointing at the same
entry point.

**A fresh clone has no data.** `data_cache/` is gitignored, so on a clone every
tissue section is unavailable and `reconstruct` refuses. One command fixes that:

```bash
vivum sections                  # what tissue exists, what is runnable, how to fix what is not
vivum acquire                   # 4 tabular sources + the xenium_319 tissue section
vivum simulate --source xenium --dx 30 --drug temozolomide   # first result
```

`vivum acquire` fetches two different kinds of thing and the difference matters:
four **tabular** sources (TCGA clinical, TCGA expression, Ivy GAP, GDSC2 potency)
that parameterise *drugs*, and one **spatial tissue section** — the thing a drug
is simulated *on*. Only the second makes `reconstruct`, `simulate` and the Studio
work. `vivum acquire --list` prints the plan without downloading; `--sections
all` adds the larger auto-fetchable sections; `--sections none` is the old
tabular-only behaviour. Sections that need extraction, a conversion or a
companion dataset are **not** auto-fetched, and `vivum sections` prints the exact
shell block for each one rather than pointing at a document.

```bash
vivum --help                    # every subcommand
vivum molecules                 # library drugs + structure-derived descriptors, no data
vivum reconstruct --source xenium --dx 12          # real GBM niche, 40,887 cells
vivum simulate --source xenium --dx 12 --drug temozolomide

# The core use case: ONE reconstruction, ONE compound, and the question
# "when this failed, did the drug not arrive, or arrive and stay too weak?"
vivum sections                  # which tissue is runnable, and how to fix the rest
vivum studio                    # the same thing with the per-cell map, in a browser

vivum validate --quick          # the scorecard, analytic tier only (seconds)
vivum validate                  # full scorecard, real tissue + PDE solves (~15-25 min)
vivum validate --list           # the anchor registry, plus what was deliberately left out
```

**`vivum validate` exits 1 on a correct install, and that is the point.** Several
anchors FAIL against measured data; they are the *published limitations* listed
under *What has failed* below, they fail identically on every machine, and a
build that reported them as passing would be the broken one. Read exit 1 here as
"the scorecard ran", not "the install is wrong". The command says so in its own
output.

> **`compare-drugs` WAS REMOVED ON 2026-08-21, WITH THE CLAIM IT EXISTED TO
> MAKE.** This block used to open on it — *"the core use case: one
> reconstruction, several compounds, ranked"* — and the cross-compound ordering
> was measured and does not hold (`results/pk_ranking_robustness.md`,
> `results/capacity_ranking_robustness.md`: per-compound plasma exposure spans
> 113.1× across the retained compounds, and there is no per-compound
> tissue-capacity term at all). The command is gone, so those two lines errored
> rather than misled — but the sentence above them still called ranking the core
> use case, which is the part that mattered. Both are replaced above.
>
> Dated files under `results/` still show `codon compare-drugs …`. They are a
> published record and are not edited; the `codon` console script is retained
> for exactly that reason.

**Pass `--dx 12` explicitly.** The CLI's fast defaults use a coarser grid, and
values quoted at dx = 30 are superseded — a coarse grid cannot resolve the
perivascular boundary layer. dx = 12 is the production resolution, and it is
**not** a converged one either; see *What is retracted*.

The two Xenium sections are fetched by `vivum acquire` (`--sections all`). The
Visium sections, GSE84465 and HGCC are **not**: they need extraction, a
conversion or a companion dataset, so no single command is honest about them.
`vivum sections` prints the literal shell block for each; the same commands are
in [`docs/REPRODUCE.md`](docs/REPRODUCE.md) §2b and the complete byte-level
record is [`docs/DATA_INVENTORY.md`](docs/DATA_INVENTORY.md).

**Studio** — the browser UI (reconstruct, compare, protocol editor, caveats
panel):

```bash
pip install -e '.[app]'
vivum studio                    # or: python -m app.server
```

**Your own tissue** — nothing here is limited to the sections this project
prepared. Point the ingester at a raw 10x export (Xenium or Visium, or a bare
MEX triple plus a coordinates table as the escape hatch for an unseen platform):

```bash
vivum ingest /path/to/your/export              # what is it, and what can it support?
vivum ingest /path/to/your/export --register   # + make it selectable in Studio
```

or open **`/ingest`** in a running Studio for the same thing in a browser. It
answers **before** anything runs, and the useful half of the answer is what your
sample *cannot* support — a panel with no endothelial marker means no oxygen or
delivery result at all, and it says so rather than producing one. A directory it
cannot use is refused with a reason instead of a broken section.

**De-identified tissue only.** Vivum is not built for identifiable patient data
and must not be used with it — no accounts, no audit log, no encryption at rest.
The sample directory is read where it lies and never uploaded anywhere, but its
absolute path, the section id and the label all become visible to anyone who can
reach the server, so keep identifiers out of those too.
[`docs/DEPLOYMENT.md` §0](docs/DEPLOYMENT.md) says where to look before you
ingest.

**`vivum validate`** is the honesty layer, and it is meant to be read as much as
run. It scores the model against independent published measurements, prints
PASS / FAIL / NOT TESTABLE / DEFERRED per row with the reasoning inline, refuses
to score a solver number whose delivery law is undeclared, and **exits 1 on any
FAIL**. It currently reports more failures than passes. That is the point of it.

**Landing page** — `cd landing-3d && npm install && npm run dev`. The guardrail
block at the top of [`landing-3d/lib/copy.ts`](landing-3d/lib/copy.ts) is the
best-calibrated claim text on any surface in this repo.

**Reproduce everything** — `./scripts/reproduce_real.sh` (add `--with-deep`,
`--list`, `--dry-run`, `--only <stage>`). Deterministic, seed 0. Setup, expected
outputs and measured runtimes: [`docs/REPRODUCE.md`](docs/REPRODUCE.md).

---

## The data

All open, all cited, all traceable. **No synthetic data.** `tmesim/data/synthetic.py`
and the fabricated-geometry loader are deleted from the package and the
reproducibility harness fails loudly if either reappears; the synthetic era is
archived under [`archive/synthetic_era/`](archive/synthetic_era/README.md).

| Layer | Dataset | Scale |
|---|---|---|
| GBM section (flagship) | 10x Xenium GBM FFPE, Human Brain 319-gene panel (CC BY 4.0) | **40,887 cells**, true µm coordinates |
| Breast section (flagship) | 10x Xenium breast, Janesick 2023 (PMID 38114474), Rep1 | **167,780 cells** (89,785 tumour) |
| GBM, held-out patient | 10x Xenium GBM FFPE, Immuno-Oncology 480-gene panel | **816,769 cells**, 6.2× the cell density |
| Spot spatial | 4 Visium GBM sections (10x public + GSE194329) | 1,775–4,764 spots each |
| Cell reference | Darmanis 2017, GEO GSE84465 | 23,460 genes × 3,589 annotated cells |
| Niche architecture | Ivy GAP (Puchalski 2018) | 270 RNA-seq samples, 7 anatomic structures |
| Drug potency | GDSC2 + HGCC (Johansson 2020, PMID 32668248) | 262 drugs × 85 patient-derived cultures |
| Physics anchors | PMID-verified tissue-PK measurements | penetration lengths, brain:plasma ratios |
| Validation corpus | `data_cache/hunt/` (~5.3 GB) | 164 datasets surveyed, 143 open — [`results/DATA_HUNT.md`](results/DATA_HUNT.md) |

Per-dataset provenance: [`docs/real_data/`](docs/real_data/).

---

## What is defensible

These survived attacks designed to break them. The wording below is the wording
that survived; it is deliberately narrow.

**1. Ratio-type outputs are substantially insulated from absolute parameter
error.** The delivery law was replaced with a corrected one (`S = PS·(Cp − C)`
instead of the one-way `S = PS·Cp`) and every breast panel drug re-run.
**Absolute delivery moved by up to 9.58×** (5-fluorouracil); **the scored 3-D:2-D
ratio moved by at most 1.45×**, and **zero verdicts flipped**. For the three
worst-affected drugs 84–86 % of the delivery error cancels in the ratio.
→ [`results/breast_kill_revalidation.md`](results/breast_kill_revalidation.md) §5

*The limit of this.* The residual is never zero: 13.8 % at best (5-FU), ~48 % for
the five low-permeability drugs. The cancellation is of *multiplicative scale*
error, which the `plasma = c_bath / p99(e_unit)` renormalisation removes partly by
construction. A *shape* perturbation is not covered — a sub-grid delivery source
moves the score materially ([`results/delivery_source_geometry.md`](results/delivery_source_geometry.md)).
The source file's own summary: *"'insulated' is not 'immune'."*

**2. The scored ratio is essentially insensitive to the timestep.** Halving `dt`
across six drug × grid combinations moves it by at most **0.0003 %**, with the
exposure field itself moving 4–5 orders of magnitude more — which is what the
renormalisation is supposed to do.
→ [`results/grid_convergence_gbm.md`](results/grid_convergence_gbm.md) §7.3.
That file states the limit itself: this licenses **nothing** about absolute AUC or
absolute kill, only about the ratio.

**3. Therapeutic coverage is *exactly* invariant to the oxygen field** — verified
with `==`, not `isclose`, across both vascular modes and three drugs, and
structurally guaranteed because oxygen is not an input to the transport solve at
all. → [`results/oxygen_level_shift.md`](results/oxygen_level_shift.md) §4

**4. Cross-drug ranking is invariant to common-mode perturbations** — all three
cell-state routes, the oxygen route, α, and the exposure unit (a 3.36× change in
the exposure unit moves decision probabilities by ≤ 0.004).
→ [`results/cell_state_impact.md`](results/cell_state_impact.md) §3.3,
[`results/oxygen_level_shift.md`](results/oxygen_level_shift.md) §4,
[`results/ranking_robustness.md`](results/ranking_robustness.md).
It is **not** invariant to drug-specific transport parameters — see *retracted*.

**5. The solver's internal numerics.** Mass is conserved on the real irregular
tissue mask to **~4e-14** (the round-off floor), concentrations never go negative,
the implicit scheme never destabilises across a 32× `dt` range, and the explicit
advection step computes its own CFL number and *raises* rather than degrading
silently. → [`results/numerical_rigor.md`](results/numerical_rigor.md) §1, §3 —
the sections that file's own correction block explicitly leaves standing.

**6. The tissue is real.** Every reconstruction above runs on measured cell
coordinates from a real section. Nothing upstream of the prediction is synthetic.

**7. Two independent real potency screens agree.** GDSC2 and HGCC, no model in the
loop, Spearman **0.827** (p = 2.35e-05) over the 18 shared drugs — which is why
the molecular layer uses measured IC50 rather than a structure→potency regression.
That same file records why the regression was rejected: a descriptor model reaches
Spearman 0.42 on real sensitivity, but the dominant feature is molecular weight and
the signal is a drug-**class** confound, not a within-class QSAR. Note this is a
check on the **data**, not on any Vivum model.
→ [`results/pharmacology_real.md`](results/pharmacology_real.md)

---

## What is contingent

Real measurements whose verdict depends on a configuration choice that is itself
unvalidated, or that rests on one section.

**The GBM 3-D:2-D kill panel passes 3/3 — under the shipped configuration only.**
Temozolomide 0.2518 (band 0.1590–0.4622), cisplatin 0.2739 (0.2086–1.5004),
doxorubicin 0.2081 (0.0575–1.0000), from a model never fitted to any 3-D-vs-2-D
measurement. Across five measured alternative configurations the panel counts are
**2/3, 2/3, 3/3, 3/3, 1/3**, and temozolomide — the matched anchor — is among the
failures in all three arms that change the count. The two 2/3 arms are independent
in provenance but land on the same lever (oxygen level) and give bit-identical
temozolomide ratios, so they are not two independent confirmations.
→ [`results/area_conservation_promotion.md`](results/area_conservation_promotion.md) §2
for the current ratios, and [`results/kill_calibration.md`](results/kill_calibration.md)
§3 / **§3.1** for the panel's construction and its alternative-configuration arms.
The three ratios above were **0.2626 / 0.2791 / 0.2077 until 2026-08-21**: that is
what `kill_calibration.md` §3 measured on 2026-08-12, before `d77d169`
(area-conserving vessel source) and `fd2a9e7` (two-route transport ON) both moved
`e_unit`. Neither commit regenerated the document, so this line quoted a superseded
table for a week. The bands and the 3/3 are unaffected; the citation now points at
the post-change measurement, which two independent computations agree on to four
decimals.
Separately, the 3/3 survives on the *width of the literature bands*, not on
numerics: `dx = 12` sits at the maximum of `F_state` over the whole grid ladder,
biasing all three ratios high by +1.2 / +2.4 / +2.7 %
([`results/grid_convergence_gbm.md`](results/grid_convergence_gbm.md)).

**The breast kill panel is 3 PASS / 1 FAIL / 6 NOT TESTABLE — settled, and
only four of ten drugs are testable at all.** This reconciles two earlier counts
that were **both wrong**. `kill_calibration_breast.md` §3 published **5 PASS /
1 FAIL** by falling back to a broad cross-tumour band for two drugs whose own
literature is a single row with no spread; the shipped `vivum validate` did not
fall back and printed **3 PASS / 3 FAIL**, scoring 5-fluorouracil and paclitaxel
against **zero-width** bands (`0.1636–0.1636`, `0.0412–0.0412`) that no
prediction can land inside — a bar nothing can fail versus a bar nothing can
pass, the mirror image of the defect for which the GBM IFP row was declared NOT
TESTABLE. A single-study band is now **NOT TESTABLE** in both paths, which can no
longer disagree. Live `vivum validate --only breast-kill-3d2d` prints
`PASS 3 FAIL 1 NOT TESTABLE 6 (of 10 rows)`. Docetaxel is the one FAIL.
→ [`results/zero_width_bands.md`](results/zero_width_bands.md)

**Docetaxel on breast is the one verdict that is a function of the grid.** Its
ratio 0.1484 sits **1.5 % outside** its band ceiling of 0.1463 — a margin smaller
than every available estimate of the grid movement of the statistic producing it,
and those estimates disagree in sign. The project's own position:
*"UNDETERMINED at dx = 12 µm — not FAIL, and not PASS"*
([`results/p99_convention.md`](results/p99_convention.md) §6.1). It survives the
corrected delivery law (1.5 % → 2.6 % outside), and it survives a cellular-exclusion
correction that raises hindrance from 1.77× to 6.61× in one arm while flipping to
PASS in the other; it does not survive a 2× grid coarsening, and a sub-grid
delivery source removes it entirely (0.1484 → 0.0752).

**Ranking under oxygen perturbation is *near*-invariant, not invariant.** Relative
kill changes span **71.4–77.0 %** against a 74 % common shift — a 5.6 pp spread
across drugs. `kill_fraction` is a covariance, not a scalar multiple, so the
"drug-independent per-cell multiplier" premise is false in general. Worse, the
oxygen → resistance map is genuinely **drug-dependent** in code:
`TMZ_HYPOXIA_RESISTANCE = 0.80` (`tmesim/core/config.py:264`) against a global
`hypoxia_resistance = 0.35` (`:238`). That drug-dependence is inert **today only
because the GBM section's hypoxic fraction is 0.0000 at every grid spacing
tested** — a non-empty hypoxic band would activate it and the near-invariance
above would have to be re-measured.
→ [`results/oxygen_level_shift.md`](results/oxygen_level_shift.md) §4,
[`results/grid_convergence_gbm.md`](results/grid_convergence_gbm.md) §2.1

**The top-1 pick is not decision-grade.** Correcting the delivery law changes the
top-1 GBM pick (cisplatin → vincristine), changes the top-3 set, and swaps 4 of 45
pairs — two of them across gaps of 40.5 % and 35.6 % of the leader's kill.
P(cisplatin rank 1) falls **0.9277 → 0.3965**.
→ [`results/transvascular_exchange_revalidation.md`](results/transvascular_exchange_revalidation.md) §4.1.
`results/ranking_robustness.md` §8: *"the reportable output of this pipeline for
the top of the panel is a SET, not a winner."*

**Every absolute this engine currently ships is computed under a delivery law the
codebase documents as measurably wrong.** `transvascular_exchange: bool = False`
(`tmesim/core/config.py:137`) selects the one-way leak `S = PS·Cp`, which has no
path from interstitium back to blood and whose steady state is unbounded above
plasma; five of eleven GBM drugs sit above plasma. The corrected law introduces no
new parameter and is bounded by plasma by construction. It is off by default for a
stated engineering reason — turning it on invalidates every published absolute and
is a re-validation exercise, not a flag flip — but that reasoning does not extend
to the claim surface, which is why it is stated here.

**The spatial-attribution statistic is a study, not a product output.** Tumour
cells in the macrophage-densest quintile retain **0.4291** of exposure against
**0.8592** in the sparsest — a raw 2.0025× spread that reproduces to the digit.
But against spatially-aware nulls that preserve the macrophages' own clustering
the excess is only **1.08–1.21×** (nulls return 1.6581 and 1.8522), **82 % of the
raw ratio is reproduced by distance-to-tissue-edge alone with zero immune input**
(1.6394), p = 0.0476 is the 20-draw floor and must not be read as a resolved
p-value, there is no confidence interval, and it is one section, one drug, one
sink strength. The module ships **default-off** and no shipped code path imports
it. → [`results/immune_sink_impact.md`](results/immune_sink_impact.md) §5

---

## What has failed

Every absolute quantity this project has scored against an independent published
measurement. These are listed so nobody has to rediscover them.

| Quantity | Predicted | Measured | Error | Source |
|---|---|---|---|---|
| Interstitial fluid pressure, breast | 0.370 mmHg | 15 ± 9 mmHg (PMID 1423283) | **41× low** | [`ifp_validation.md`](results/ifp_validation.md) §4.2 |
| IFP, GBM | 0.462 mmHg | 2.0 ± 2.5 mmHg | 4.3× low, but registered **NOT TESTABLE** — a model outputting zero would pass | [`ifp_validation.md`](results/ifp_validation.md) §4.1 |
| Oxygen `supply` constant | 0.05 | ≈194–251 s⁻¹ (PMID 31423930) | **~4,570× low** on GBM (2,160× on breast), partly cancelled by a consumption sink **112–166× low** | [`oxygen_supply_calibration.md`](results/oxygen_supply_calibration.md) §3.3, §5 |
| Effective diffusivity, multicellular layer | 287.3 µm²/s | 74–130 µm²/s | FAIL | [`anchor_registration.md`](results/anchor_registration.md) §1 |
| **Tissue diffusivity `D_eff` as a *law*, 71 compounds** | 8.3× spread, **R² = 0.718** as shipped (was 2.4×, R² = 0.022) | 29.7× spread (Pruijn 2008, doi:10.1071/CH08240) | **PARTIALLY SUPPORTED** — blind and unfitted, 99 % within 3×, and **24 of 31** exact-mass pairs ordered correctly against the superseded model's **0 of 31**. Still only ~**24 %** of the measured spread, and ONE cell line. This is the row that bears on ranking a compound series | [`two_route_promotion.md`](results/two_route_promotion.md) |
| Hindrance `D_free/D_eff` | 1.77× (superseded) → **10.33×** as shipped | 12× (PMID 9806526) | **PASS** — inside the band, 0.86× the measured value. The one external transport anchor the two-route promotion moved from FAIL to PASS | [`two_route_promotion.md`](results/two_route_promotion.md) |
| Effective diffusivity, tirapazamine | 287.3 (superseded) → **49.3 µm²/s** as shipped | 74–130 µm²/s (PMID 9806526, 17869669) | **FAIL, direction flipped** — was 2.21× above the band, now 1.50× below it | [`two_route_promotion.md`](results/two_route_promotion.md) |
| Breast Kp, 5-fluorouracil | 7.7956 | 0.61 ± 0.11 (PMID 9205062) | **12.8× high** | [`breast_kp_validation.md`](results/breast_kp_validation.md) |
| Breast Kp, methotrexate | 0.0571 | 0.60 ± 0.20 (PMID 9679959) | **10.5× low** — opposite direction | [`breast_kp_validation.md`](results/breast_kp_validation.md) |
| Vascular area fraction, held-out patient | 10.982 % | 4.70 ± 2.38 % (PMID 27583442) | **2.34×** — but see below | [`vascular_segmentation_audit.md`](results/vascular_segmentation_audit.md) §2.3 |
| Microvessel density, held-out patient | 233.27 /mm² | 65 (35–95) /mm² | **3.59×** (1.70× on section 1 — the merge factor is not a constant) | [`heldout_patient.md`](results/heldout_patient.md) §2 |
| Regional CE vs non-enhancing separation | exactly 1.00× | 1.45–4.26× (p < 0.0001) | FAIL — there is no barrier-integrity term in the shipped source | [`glioma_regional_pk.md`](results/glioma_regional_pk.md) |
| Passive transport as a sufficient explanation of tissue partition | 7.9× geometric-mean fold error | vs 13.0× for a no-partitioning null | **FALSIFIED, diagnostically** — beats the null by only ~1.6×, and 4 of 5 liver pairs lie below the model's structural floor | [`gist_imatinib.md`](results/gist_imatinib.md) |

**The vascular-area row is a trap and it cuts against us.** The commensurable
anchor (Marien 2016 areal fraction, 4.70 %) was **added** to the scorecard, not
substituted for the incommensurable one (Tamma 2022 CD34 stain-positivity 2.1 %,
where the lumen contributes no stained pixels but the model paints the whole vessel
profile). `scorecard.py` still reads the 2.1 % key, so **`vivum validate` still
prints 5.23×**. Do not quote 2.34× as the product's number until the scorecard
scores it. Either way it is a FAIL.

**Three scored rows compare a *solved field* to an *external measurement*, and all
three FAIL — in opposite directions**: breast IFP 41× low, 5-FU partition 12.8×
high, methotrexate partition 10.5× low. No single constant can absorb them.
Refinement makes IFP *worse*, not better: breast IFP falls monotonically
0.621 → 0.424 → 0.370 → 0.286 mmHg across dx = 24/18/12/8, a 2.17× spread
([`ifp_validation.md`](results/ifp_validation.md) §5).

**And the scorecard mostly does not test the solver at all.** Live
`vivum validate --quick` on the shipping physics returns
**PASS 7 · FAIL 6 · NOT TESTABLE 4 · DEFERRED 8** of 25 rows, and prints its own
disclosure: ***13 of the 13 scored rows never invoke the transport PDE.*** The full
run does reach the solver — 12 PDE-scoped scored rows — but every PASS among them
is a 3-D:2-D ratio conditional on the one fitted constant α. The full-run tally
after the zero-width-band fix is **PASS 13 · FAIL 18 · NOT TESTABLE 19** of 50
rows (it was FAIL 20 / NOT TESTABLE 17 before; two rows moved FAIL → NOT
TESTABLE and **none moved to PASS**).
→ [`results/zero_width_bands.md`](results/zero_width_bands.md) §5.
Roughly **18 distinct FAIL rows** exist across `results/`; they are enumerated in
[`results/CLAIM_INVENTORY.md`](results/CLAIM_INVENTORY.md) §1e. Docetaxel is one
falsification among many — this project's honesty record is better than "one
instance", and should be stated that way.

---

## What is retracted

Each of these was published by this project and is now measured false or
superseded. They are kept visible because the retractions are load-bearing.

| Retracted claim | What replaced it |
|---|---|
| *Ranking is unconditionally invariant.* | Invariant to **common-mode** perturbations only. A drug-specific transport parameter changes the top-1 GBM pick. → [`transvascular_exchange_revalidation.md`](results/transvascular_exchange_revalidation.md) §4.1 |
| *dx = 12 µm is grid-converged.* | **Neither panel is in an asymptotic regime.** Breast: docetaxel is −3.4 % from its Richardson limit, its boundary layer is 28.2 µm = 2.4 cells at dx = 12, ladder spread 49.7 % ([`grid_convergence_breast.md`](results/grid_convergence_breast.md)). GBM: the ladder is **non-monotone**, observed orders 0.96 and −3.33, GCI 5.44–6.42 %, power-law fits flagged DEGENERATE ([`grid_convergence_gbm.md`](results/grid_convergence_gbm.md) §0). |
| *dx = 12 is the finest grid satisfying the advection CFL at dt = 10 s.* | **Measured false.** dx = 10 runs at CFL 0.663; the GBM floor is **dx = 7** (CFL 0.955), with dx = 6 the first grid the solver refuses. Re-measured independently on the shipped configuration: dx = 8 at 0.939, dx = 7 at 0.955. → [`grid_convergence_gbm.md`](results/grid_convergence_gbm.md) §2.1, [`cfl_resolution.md`](results/cfl_resolution.md) §0. The present-tense comment at `tmesim/reconstruct/real_fields.py:146-149` describes a superseded reconstruction. |
| *The mesh is not the bottleneck once dx ≲ 12 µm* / **Numerics: VALIDATED**. | Withdrawn by [`numerical_rigor.md`](results/numerical_rigor.md)'s own correction block. The elasticities were mis-attributed (|E_pen| 0.747 belongs to `grid_dx_um`, not the endothelial threshold; |E_kill| 0.610 belongs to α) and the §2 grid trend had the sign backwards. Conservation, positivity and dt behaviour still stand; the headline does not. |
| *`penetration_index` is a delivery number.* | **Retired.** It is `mean(AUC)/p99(AUC)` — a **scale-invariant shape statistic**: multiply the whole exposure field by any constant and it does not move, so on real tissue it reports the *vascular geometry of the section*, which is shared by every drug. At dx = 12 µm it spans 0.304–0.511 across the panel while those drugs' mean delivered exposure differs by **158×**. Use `therapeutic_coverage` at 1× IC50 (spread 0.995, separation 67.7× its own grid drift, against 15.2× for `penetration_index`), quoted with its 0.5×/2× band. → [`delivery_metrics.md`](results/delivery_metrics.md) §4 |
| *The oxygen solve is out-predicted by a distance-to-vessel baseline, so demote it.* | **Reversed — the premise was an inverted inequality.** The solve out-predicts its baseline by 1.26× (ρ = −0.1523 vs −0.1204). But ρ ≈ 0.15 is weak, within-section, one patient, and the composite is carried by VEGFA alone — HIF1A's ρ is **+0.0140, the wrong sign**. → [`oxygen_demotion.md`](results/oxygen_demotion.md), [`heldout_patient_validation.md`](results/heldout_patient_validation.md) |
| *The vascular field is a smooth KDE, so every vessel-relative claim is confounded.* | **Superseded configuration.** `_DEFAULT_VASCULAR_MODE = "discrete"` (`tmesim/reconstruct/real_fields.py:179`) — segmented vessels, not a blur, since commit `ceb96c0`. The promotion also grew the oxygen route's downstream effect ~12×. Neither vascular mode is grid-converged at dx = 12. |
| *Cross-platform composition concordance, Pearson 0.997.* | After a classifier-arity fix the 4-vector is 0.980 and the **non-tumour subvector collapses 0.945 → 0.016** — indistinguishable from no agreement. → [`revalidation.md`](results/revalidation.md) |
| *Transport physics validated out-of-sample against measured tissue PK.* | Those are `analytic`-scope rows that **never invoke the solver**. `scorecard.py`'s own warning: *"'Robust under exchange' and 'never evaluates the delivery law' are very different statements."* They can falsify a coefficient law; they cannot confirm or refute the solver. |
| *A barrier-integrity fix is specified but not implemented.* | It **is** implemented — `L_PARACELLULAR_LEAK = 0.238` in `tmesim/reconstruct/barrier.py`, reachable via `treatment.drugs.delivery_source(..., bbb_disruption=...)`. It is **opt-in and OFF by default**, so every shipped result still predicts a 1.00× regional separation. Note it was *fitted* to the same contrasts it is scored against, so leak-on scores are not out-of-sample. |

One recovery, in the other direction: the held-out patient's endothelial fraction
now reads **7.07 % against a volume-EM measured 8.0 %** (PMID 38723085) after the
classifier fix and the discrete-vessel promotion, where it previously did not
replicate. The owning file qualifies even this: the agreement should be read as a
property of that section and that marker gate, not as a validated property of the
classifier ([`heldout_patient.md`](results/heldout_patient.md) §2).

---

## The honest summary

Vivum is a research instrument whose defensible output is the **relative
comparison of compounds on one tissue section**, with essentially **no external
validation of the transport physics itself**. The narrow claim has been attacked
from several directions and held — with one measured exception that a user must
know before running anything: **the diffusivity law reproduces about a quarter
of the measured compound-to-compound spread, so a ranking over similar-mass
analogues is weak evidence rather than a measurement** (see the box at the top of
this page). The broad claim — that the transport physics is validated, or that
any absolute number it produces is trustworthy — has been tested and is false.

The sharpest available statement of that gap is
[`results/external_anchor_survey.md`](results/external_anchor_survey.md) §5, and
it is worth reading in full. Its verdict: the solver's **input coefficients**
(`D_eff`, `k`, and the pressure field) *are* externally testable against public
absolute measurements — and on first contact the engine misses all of them. The
solver's **output fields** are a different matter: nothing public was found that
measures a spatially resolved free-interstitial concentration field in tumour
tissue that this engine could reconstruct. So *"our transport coefficients are
testable and currently fail in identifiable, structural ways"* is defensible;
*"our predicted concentration fields are validated"* is not, and no public data
currently exists that would make it so.

`vivum validate` reports more failures than passes and exits non-zero. That is
not a defect being disclosed reluctantly; it is the product's most reliable
component. The modules that could flatter a result and the modules that check it
are separate, and the checking modules are allowed to fail loudly. Every
retraction on this page was found that way, not by the test suite, which was green
throughout.

**What would most improve this, in order:** (0) score the transport law against
a **second cell line** — the two-route promotion (R² 0.022 → 0.718, blind) rests
on Pruijn's HT-29 measurements alone, and Pruijn's own paper reports its
asymptote coefficients moving up to 100-fold across SiHa, HCT8-Ea and HCT8-Ra, so
this is the cheapest unbought evidence in the project and the only item on this
list that changes what a *ranking* means; (1) a **second section of the same
tumour type** —
almost every defensible claim above is *n* = 1 in tissue, and this is the cheapest
large improvement available; (2) any external anchor that scores the solver on the
strength of a solve **and passes** — there are three such anchors and all three
fail, and [`results/external_anchor_survey.md`](results/external_anchor_survey.md)
§5 reports that nothing public measures a spatially resolved free-interstitial
concentration field in tumour tissue, so this may require generating data rather
than finding it; (3) flip `transvascular_exchange` and re-validate, or publish the
caveat on every absolute; (4) resolve `dx`; (5) re-derive the oxygen `supply`
constant for the discrete-vessel referent.

---

## Standing limitations

* **The response axis cannot currently be validated.** A data-access gap, not an
  engineering gap. → [`results/DATA_HUNT.md`](results/DATA_HUNT.md)
* **2D.** A whole-slide section is naturally 2D, but out-of-plane vascular
  connectivity is invisible: a vessel just above the plane makes a point look
  "far" when it is not. Distances are Euclidean, not geodesic.
* **One anchor sets the delivery scale.** α is fit once, in closed form, to a
  single doxorubicin perivascular measurement (Primeau 2005, PMID 16361566) whose
  own uncertainty is plausibly ±50–100 %. It is the **largest** local kill
  elasticity (|E_kill| = 0.610) and the second-largest penetration elasticity.
  Everything downstream inherits that one fit.
* **`D_tissue` and the perivascular sink strength are not separately identified**
  ([`dtissue_sensitivity.md`](results/dtissue_sensitivity.md) §8), so part of the
  observed robustness to a diffusivity change is a gauge artefact — α is re-fit to
  the same anchor and re-exports the change as a global multiplier.
* **Fields are inferred, not measured.** Vasculature, ECM and oxygen are
  reconstructed from real markers. `vascular_density` is max-normalised: it carries
  geometry, not an absolute perfused fraction, and is not comparable across tissues.
* **Temozolomide's resistance determinant is not spatially resolvable.** MGMT is
  absent from **both** flagship Xenium panels, so a determinant-mode comparison
  would quietly advantage TMZ; the reconstruction default is therefore uniform
  sensitivity 1.0 for every drug — the only fair cross-drug mode, and a documented
  no-resistance upper bound rather than a measurement. The stronger claim that
  MGMT is *unmodellable from spatial data* was **retracted**: a third public GBM
  Xenium section (GSE324268) does carry an MGMT probe and detects it. It is still
  not usable — 1.19 % of tumour cells positive, max 4 counts in any cell, and its
  spatial clustering falls **below** the assay's own negative-control probes when
  thinned to matched depth. The defensible position survived the test: MGMT is a
  per-patient **scalar**, and is available as one.
  → [`results/mgmt_spatial.md`](results/mgmt_spatial.md), [`results/mgmt_scalar.md`](results/mgmt_scalar.md)
* **The oxygen field's constants are uncalibrated, and the hypoxic band is
  currently empty.** GBM's hypoxic fraction is 0.0000 at every grid spacing
  tested under the shipped configuration; under a literature-derived `supply` it
  is 40–71 %. Treat any hypoxic fraction as a configuration output, not a
  measurement.
* **Reduced biology.** No adaptive resistance, no immune dynamics, and systemic PK
  is a simplified sustained plasma level. Active efflux (P-gp / BCRP) **is**
  implemented from literature efflux ratios and is opt-in — but default-off is the
  *validated* recommendation, because turning it on does not improve the measured
  anchor. ABCB1 and ABCG2 are on the explicit ABSENT list for both Xenium panels
  (`tmesim/reconstruct/real_xenium.py:693`), so nothing spatial constrains it here.
  → [`results/efflux_impact.md`](results/efflux_impact.md)
* **One section per tumour, no replicates.** Where a claim rests on a single
  section, the write-ups say so.
* **Research prototype.** No diagnostic, regulatory or treatment-selection use.

---

## Module map

```
data/real/  ─►  reconstruct/real_*  ─►  Microenvironment  ─►  physics/ + cells/  ─►  sim/  ─►  validation/  ─►  viz/
(real bytes)     (real niche)            (core object)         (transport + kill)     (engine)   (vs real data)   (figures)
                                              ▲
                                    molecular/ (structure → DrugParams)
```

| Package | Responsibility |
|---|---|
| `tmesim.core` | Frozen data model: `Domain`, `Microenvironment`, `CellPopulation`, `TissueFields`, `DrugParams`, the config dataclasses and the `*Result` containers. |
| `tmesim.data.real` | Loaders for the real caches: `singlecell_spatial` (both Xenium panels), `visium`, `scref_gbm` (GSE84465), `ivygap`, `drug_response` (HGCC + GDSC2). |
| `tmesim.molecular` | Structure → physics: verified SMILES, descriptors, Stokes–Einstein radius, Ogston–Renkin mesh hindrance, permeability/BBB score, measured GDSC potency, `compile_drug`. |
| `tmesim.reconstruct` | Real tissue → real niche. `real_xenium` (marker typing, geometric tissue mask), `marker_panels`, `real_visium` (NNLS deconvolution), `real_fields` (vascular/ECM/oxygen fields; **discrete** vessels by default), `vessel_network` (the segmentation), `barrier` (the opt-in paracellular leak), `real_pipeline`. |
| `tmesim.physics` | `diffusion` (convection–diffusion–reaction transport), `oxygen` (steady-state O₂ + hypoxia banding), `pressure` (IFP and the Jain-effect velocity). |
| `tmesim.cells` | `agents` (oxygen → cell state), `pkpd` (local AUC → Hill kill probability). |
| `tmesim.treatment` | `drugs` — lookup, plasma PK profile, spatial delivery source. |
| `tmesim.sim` | `engine.simulate_response` — the deterministic per-tissue integrator. |
| `tmesim.validation` | The honesty layer, and where this project spends most of its code. `scorecard` (what `vivum validate` runs), `calibration` (the single α), `numerical_rigor`, `kill_calibration`, `delivery_metrics`, `ranking_robustness`, plus the external-measurement modules: `pk_multidrug`, `glioma_regional_pk`, `gist_imatinib`, `gbm_tmz_exvivo`, `moncayo_explant`, `pharmacology_real`, `heldout_patient`. |
| `tmesim.cli` | `reconstruct`, `simulate`, `sections`, `validate`, `molecules`, `molecular-run`, `acquire`, `studio`. |
| `experiments/` | Runnable real-data studies: flagship simulation, drug panel, cross-tissue, niche response, combination, dose/schedule, per-patient report. |

The frozen public interface: [`ARCHITECTURE_CONTRACT.md`](ARCHITECTURE_CONTRACT.md);
where it and [`docs/architecture.md`](docs/architecture.md) disagree, the contract
wins. The molecular layer's own contract: [`MOLECULAR_CONTRACT.md`](MOLECULAR_CONTRACT.md).

---

## Where to look next

| Question | File |
|---|---|
| **What does this project claim, and what backs each claim?** | [`results/CLAIM_INVENTORY.md`](results/CLAIM_INVENTORY.md) — the single settled answer, and the source for most of this page |
| New here — install to first result | [`docs/QUICKSTART.md`](docs/QUICKSTART.md) |
| How do I reproduce all of this? | [`docs/REPRODUCE.md`](docs/REPRODUCE.md) |
| What is the model, mathematically? | [`docs/MATHEMATICAL_MODEL.md`](docs/MATHEMATICAL_MODEL.md) |
| Where did every byte of data come from? | [`docs/DATA_INVENTORY.md`](docs/DATA_INVENTORY.md) |
| What open data exists to validate this, and what does not? | [`results/DATA_HUNT.md`](results/DATA_HUNT.md) |
| How should I compare drugs? | [`results/delivery_metrics.md`](results/delivery_metrics.md) |
| Is the solver trustworthy? | [`results/numerical_rigor.md`](results/numerical_rigor.md) + its correction block |
| Is the grid fine enough? | [`results/grid_convergence_gbm.md`](results/grid_convergence_gbm.md), [`results/grid_convergence_breast.md`](results/grid_convergence_breast.md) |
| Scripting drug/dose/sweep experiments | [`PROTOCOL_LANGUAGE.md`](PROTOCOL_LANGUAGE.md) + Studio's protocol tab |
| Running Studio for a team | [`docs/DEPLOYMENT.md`](docs/DEPLOYMENT.md) |

[`results/SESSION_AUDIT.md`](results/SESSION_AUDIT.md) is a historical record of an
earlier audit pass. Several of its statements have since been superseded; where it
and [`results/CLAIM_INVENTORY.md`](results/CLAIM_INVENTORY.md) differ, the
inventory is later and wins — and where the inventory and a primary results file
differ, the primary file wins.
