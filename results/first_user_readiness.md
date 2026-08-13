# First-user readiness

**Question.** What actually stops a competent cancer-biology or DMPK researcher
from running Vivum on their own tissue and their own compounds this week and
getting a result they would put in a deck?

**Method.** The journey was walked, not modelled. Install read fresh from the
README, Studio started, a section that this project did not prepare ingested and
run, a compound in neither library compiled, a report generated and read as its
reader rather than its author, and the CLI ranking path run to completion. Every
item below is something that was hit. Where something works, it says so.

One caveat on the evidence: the Studio reachable on :8000 was a long-running
process started 2026-08-05, so the generated report's *numbers* come from an
older build than the CLI's. Every claim about **code** below was checked against
current source at `HEAD` and cites file and line; only the report's numeric
values carry that caveat, and it is flagged again where it matters.

**Scope note.** Two build agents are live on *bring-your-own-tissue* and
*bring-your-own-compound*. Those two questions are answered directly below
because the brief asked for a yes/no, but the ranked findings deliberately
concentrate on what is **outside** both lanes.

---

## The two direct answers

### Can a user bring their own tissue today? **Yes — but only through `curl`.**

This was the finding expected to be the headline, and it is not. The ingestion
path is real, generic and good.

A stock 10x Xenium Onboard Analysis bundle was assembled outside `data_cache/`
with canonical vendor filenames (`cells.csv.gz`, `cell_feature_matrix/`,
`gene_panel.json`) and offered to the ingester cold:

```
detect()  ->  platform "xenium", layout "xenium_bundle_extracted",
              confidence 1.0, issues []
POST /api/sections/scan {register:true}  ->  ok:true in 0.7 s
              118,752 cells, 313 genes, extent 7519.8 x 5471.0 um
POST /api/runs  ->  completed in 45.6 s
```

`tmesim/data/ingest.py` recognises `xenium`, `visium`, a bare `mex` triple plus a
coordinates table as the escape hatch for an unseen platform, and
`multi_sample` containers. It refuses with a reason rather than building a broken
`Microenvironment`, and it will not shadow a built-in section id. The capability
report it returns before anything runs is the best thing in the product: it told
me, unprompted, that this panel carries no hypoxia transcript, that
`oxygen_mode='markers'` was therefore unavailable, that MGMT was absent, and that
**absolute kill was not quotable from this sample under any state mode** — and
that same text then reappeared as a warning on the run I submitted.

The gap is exposure, not capability. `sections/scan` appears **nowhere** in
`app/web/` — the Studio's section picker is a read-only list from `GET
/api/sections`. So the answer to "can a researcher bring their own tissue" is
yes for one who will hand-POST JSON, and no for everyone else. *(Owned by the
BYO-tissue build.)*

### Can a user bring their own compound today? **Yes in Python, yes in the protocol editor, no in the main flow.**

Osimertinib — a real EGFR TKI in neither library — compiled from SMILES cleanly:

```
descriptors.from_smiles -> mw 486.58, logp 3.69, tpsa 100.44, r_h 0.5191 nm
build_profile           -> D_free 364.5, D_eff 193.5 um2/s, hindrance 0.531
molecule_to_drugparams  -> DrugParams(ic50=10.342, bbb_penetrance=0.503, ...)
```

It resolved a **measured** GDSC IC50 (median over 50 glioma cell lines) rather
than inventing one, and every derived field carried its own provenance string
with an explicit confidence — `free_diffusivity` "high (first-principles physics
given the radius)", `bbb_penetrance` "low (structure-only; known to under-predict
temozolomide)". That is exactly right.

But `GET /api/drugs` returns the 13-drug library and nothing else, and
`POST /api/runs` takes `drug` as a library name. The SMILES path reaches the
Studio only through the TPL protocol editor (`drug X from smiles "..."`). A user
in the ordinary Compare flow cannot add a compound. *(Owned by the BYO-compound
build.)*

---

## Ranked findings

### Blocks a first user

**1. Three different default exposure windows are in force in one product, and
the ranking command discloses none of them.**

Measured in current source:

| surface | default exposure window | where |
|---|---|---|
| `vivum compare-drugs` / `simulate` | **1,800 s** (30 min) | `tmesim/cli.py:88` `_FAST_TOTAL_TIME` |
| Studio + every report | **21,600 s** (6 h) | `app/server/api.py:399` `3600.0 * 6` |
| the library itself | **43,200 s** (12 h) | `PDEConfig().total_time` |

That is a **24x** spread end to end, and a **12x** difference between the
README's headline command and the Studio. The CLI's own docstring
(`tmesim/cli.py:67-68`) says it uses a reduced window and that "the full-fidelity
window lives in the library defaults" — but the Studio, the surface that produces
the exportable report, uses neither the CLI's nor the library's value.

Both surfaces were run for real on the same section (`xenium_319`), same drugs,
same `dx = 12`, and every drug disagreed:

| drug | `vivum compare-drugs` (1,800 s, window undisclosed) | Studio report `rep_00006` (21,600 s, window printed) |
|---|---|---|
| cisplatin | 0.4711 | 0.4362 |
| doxorubicin | 0.3199 | 0.2595 |
| temozolomide | 0.0240 | 0.0184 |

**Read that table as illustrative, not as a measurement of the window effect.**
The Studio available on :8000 was a process started 2026-08-05, and the engine
changed after it (`5d6e8e4` `_TUMOR_SIGMA` units fix, `d365e24` cellular volume
exclusion in the hindrance law, among others), so these deltas confound the
window with eight days of engine changes and cannot be attributed to either
alone. A clean measurement would need both surfaces on one commit.

The *finding* does not depend on that table: the three default windows are read
directly from current source, and they differ. And the confound itself is the
point of finding 9 — nothing on either surface would have let a user notice
that they were comparing two different builds.

`cmd_compare_drugs` never prints `total_time` — verified by inspection, the
substring occurs twice in the function and both are the `_build_sim_config` call.
`cmd_simulate` and `cmd_molecular_run` *do* print their window
(`tmesim/cli.py:222`, `:997`); the comparison command, the one whose entire
purpose is a ranking, does not. The only statement of the reduction lives in the
module docstring at `tmesim/cli.py:67-68`, which no user reads.

A first user who runs the README command and then opens the Studio gets two
different numbers for the same thing with no explanation on either surface. That
is the fastest available way to lose them.

**2. The CSV a user exports carries no provenance at all.**

This is the complete file `--csv` produced on the run above — nothing is elided:

```
drug,mw,ic50_uM,bbb_penetrance,diffusivity_um2_s,kill_fraction,penetration_index,mean_auc_uM_s
cisplatin,298.03,24.550285905059777,0.8576798963736021,236.18442715353106,0.47106500923903916,0.42381404756834784,211508.63959439844
doxorubicin,543.525,0.5,0.020567870235594195,184.74861430191334,0.3199227280362842,0.3030674593779195,1404.256710236273
temozolomide,194.154,606.0466722391783,0.12537158358322562,277.9851085523925,0.024021501763816566,0.4959460250277222,78452.08216346582
```

Eight columns, no sidecar, no comment header (`tmesim/cli.py`, the `if args.csv:`
block in `cmd_compare_drugs`). Absent: the section, `dx`, the exposure window,
`alpha`, the delivery law, the sensitivity mode, and every word of "not for
clinical use". Present, at full float64 precision: `kill_fraction` — the quantity
the README says is not defensible — and `penetration_index`, which is retired.

This is the only machine-readable export in the product, and it is the artefact
most likely to become a slide. Everything the audit built is stripped at exactly
the moment the numbers leave the tool.

### Makes them distrust it

**3. The report asserts a grid-convergence claim the README explicitly retracts,
and files it as reassurance.**

`app/server/report.py:303-324`, current code, ships this limit on every report:

> "This report ran at dx = 12 um. kill_fraction is **grid-converged for dx <= 20
> um**; penetration_index is NOT converged above ~15 um..."

with `"source": "results/numerical_rigor.md §2"` and, because `dx = 12` satisfies
`converged = float(dx) <= 15.0`, `"severity": "info"`.

The README's *What is retracted* table retracts both halves:

* "*dx = 12 µm is grid-converged.*" → "**Neither panel is in an asymptotic
  regime.**" GBM's ladder is non-monotone, observed orders 0.96 and −3.33, GCI
  5.44–6.42 %, power-law fits flagged DEGENERATE.
* "*The mesh is not the bottleneck once dx ≲ 12 µm* / **Numerics: VALIDATED**" →
  withdrawn, and "**the §2 grid trend had the sign backwards**" — §2 being the
  exact section this limit cites as its authority.

So the one document a researcher would export tells them, in the section headed
*"Limits — read these before any number on this page"*, at the lowest severity
the renderer has, that the production resolution is converged. It is the only
place in the product where a caveat surface is *less* honest than the README.

**4. The retired delivery metric has a column; its recommended replacement is
computed and never shown.**

`penetration_index` is under *What is retracted* in the README — "**Retired.** It
is `mean(AUC)/p99(AUC)` — a **scale-invariant shape statistic**... on real tissue
it reports the *vascular geometry of the section*, which is shared by every
drug", with the instruction "Use `therapeutic_coverage` at 1× IC50".

In current code:

* `app/server/report.py:499` puts `penetration_index` in the row dict and `:1217`
  renders it as the **"pen. index"** column of the ranked panel table.
* `app/server/report.py:502` computes `coverage_1x` from
  `therapeutic_coverage_1x` into the row dict, and it is **stored** in the report
  document (present in `app/reports/rep_00006/report.json`) — but it is not
  referenced by any render function in `report.py` and does not appear anywhere
  in the rendered HTML. It is carried and then dropped at the last step.
  *(To be precise: `therapeutic_coverage_1x` is well used elsewhere — the Studio
  project-comparison view labels it "Coverage @ 1x IC50"
  (`web/js/ui/projects.js:39`), the protocol layer exposes it, and the validation
  modules treat it as the headline delivery metric. The defect is specific to the
  report and the CLI, which are the two exportable surfaces.)*
* `tmesim/cli.py` shows `pen_idx` as a column of the comparison table and in the
  CSV; `therapeutic_coverage` appears in neither.

Both generated outputs demonstrate the retired metric's exact failure mode on
their own front page — temozolomide carries the **highest** penetration index of
the panel while ranking **last**:

| surface | TMZ pen. index | TMZ rank | TMZ kill |
|---|---|---|---|
| Studio report `rep_00006` | 0.511 (highest of 3) | 3 of 3 | 0.0184 |
| `vivum compare-drugs` | 0.496 (highest of 3) | 3 of 3 | 0.0240 |

A reader who does not know the metric was retired reads that column as "TMZ
penetrates best" — which is precisely the misreading
`results/delivery_metrics.md` §4 retired it to prevent.

**5. The report's headline names a single winner, and its own page 1 contradicts
it.**

`app/server/report.py:1363-1381` renders a hero card: *"Highest predicted kill on
this specimen: cisplatin"*, followed by *"The ORDER is the claim; the LEVEL is
not."*

On the same page, limit 2b (`report.py:143-159`, added 2026-08-12) correctly
states the ranking is **law-dependent**. And the README is blunter: "**The top-1
pick is not decision-grade.**" Correcting the delivery law changes the top-1 GBM
pick — cisplatin → vincristine — and `P(cisplatin rank 1)` falls 0.9277 → 0.3965;
`results/ranking_robustness.md` §8 says "the reportable output of this pipeline
for the top of the panel is a **SET**, not a winner."

The hero card names, as the product's single most prominent output, the exact
pick the project has measured to be unstable — and defends it with "the ORDER is
the claim", which is the one part of the ranking that is *least* defensible at
the top. The correct caveat exists 200 pixels below the incorrect one.

**6. The delivery-law disclosure reaches every Studio surface and misses the CLI
ranking.**

`app/server/delivery_law.py` is genuinely well built — it reads
`PDEConfig.transvascular_exchange` at call time, so the disclosure rewrites
itself rather than going stale. Verified present at: run provenance and caveats
(`jobs.py:563-567`, inserted at index 0), project rankings
(`projects.py:750,758,810`), report limits and document (`report.py:151,660`),
the project UI (`web/js/ui/projects.js:481`), and `vivum validate`
(`cli.py:480-492`).

It is absent from `vivum compare-drugs` — zero occurrences in the function — and
from `experiments/real_drug_panel.py` and `experiments/breast_drug_panel.py`
(zero mentions each). The README calls `compare-drugs` "**The core use case: one
reconstruction, several compounds, ranked.**" The one surface named as the core
use case is the one surface the disclosure did not reach.

Also worth a line: the Studio's always-on `STANDING_LIMITS`
(`web/js/ui/provenance.js:81`) has three entries and the delivery law is not one
of them. In practice a run's own caveats cover it, so this is redundancy rather
than a hole.

### Annoying

**7. The README's install block is wrong about its own `pyproject.toml`.**

README line 45 tells the user to run
`pip install rdkit openpyxl statsmodels xlrd pyxlsb rdata`, and lines 49-51 say
"The second `pip` line is **not in `pyproject.toml`**". Checking:

* `rdkit>=2025.9.2` is a **core dependency** (`pyproject.toml:62`) — already
  installed by `pip install -e .`.
* `statsmodels` and `openpyxl` are the **`[validation]` extra**
  (`pyproject.toml:105-113`).
* Only `xlrd`, `pyxlsb`, `rdata` are genuinely undeclared.

`docs/REPRODUCE.md` is correct and current — §1 gives
`pip install -c constraints.txt -e '.[validation]'` then
`pip install -c constraints.txt xlrd pyxlsb rdata`. The README never mentions the
`[validation]` extra at all, and `docs/QUICKSTART.md:44` repeats the same stale
line in shorter form. Nothing breaks — it over-installs — but the README states
something checkably false about a file two directories away, on a page whose
entire credibility rests on being checkable.

**8. `python -m app.server` reports success when it fails to bind.**

`app/server/__main__.py:97-99` calls `uvicorn.run(...)` and unconditionally
`return 0`. With port 8000 occupied the process printed the bind error and then
printed its normal banner — `serving http://127.0.0.1:8000` — and exited **0**.
The banner is printed before `uvicorn.run` (lines 77-88), so the last thing on
screen is a success message for a server that is not running.

**9. Nothing identifies which build of the code a Studio is running.**

The server on :8000 had been up since **Aug 5**; the delivery-law disclosure
landed **Aug 12**. Runs submitted to it silently lacked the disclosure. This cost
real time here — the missing caveat looked exactly like a product bug and was
only excluded by comparing the process start time to `git log`. `GET /api/health`
returns cache and queue counters and no version; the only `version` in `api.py`
(`:1139`) is the static API-contract version. A user running a team Studio
(`docs/DEPLOYMENT.md`) has no way to tell a current server from a stale one, and
the failure mode is silently weaker disclosure.

---

## What works, and should not be re-litigated

Stated because a short honest list is worth more than a long one, and because
these were tested rather than assumed:

* Install is clean. `pip install -e .` registers `vivum`; all six README-listed
  extras import; `vivum --help` and every subcommand resolve.
* `vivum validate --quick` reproduces the README's tally **exactly** — `PASS 7
  FAIL 6 NOT TESTABLE 4 DEFERRED 8 (of 25 rows)`, prints "13 of the 13 scored
  rows never invoke the transport PDE", and exits 1. The README is accurate here
  to the digit.
* The ingest capability report is the strongest asset in the product and it
  propagates: panel gaps measured at ingest reappeared verbatim as run warnings.
* The report is otherwise a serious document — a `NOT FOR CLINICAL USE` banner
  top and bottom, limits **before** any number, an exactly additive
  delivery/potency/O2-resistance attribution checked for faithfulness against the
  simulated kill, a real counterfactual simulation behind the resistance term,
  per-drug IC50 provenance with `n` cell lines, and a self-contained single-file
  HTML export with no external requests. Findings 3, 4 and 5 are defects *in a
  good document*, not reasons to rebuild it.

---

## The single highest-value thing to build next

**Make one run's disclosure travel with its numbers, and make the two ranking
surfaces agree.**

Concretely, and in this order:

1. Delete the retracted convergence limit (`report.py:303-324`) and render
   `coverage_1x`, which is already computed, in place of the retired `pen. index`
   column. This is subtraction, not construction, and it removes the only place
   where the product claims more than the README.
2. Give `compare-drugs` the same three lines the Studio already has: print the
   exposure window, print the delivery law via `delivery_law()`, and stop
   defaulting to a window no other surface uses.
3. Write a provenance header into the CSV — section, `dx`, window, `alpha`, law,
   sensitivity mode, and the not-for-clinical-use line.

The reason this outranks anything else outside the two build lanes: the project
has just spent two days establishing precisely what it may and may not claim, and
that work currently stops at the README. The two artefacts a first user actually
leaves with — the exported report and the exported CSV — are the two places the
audit has not yet reached. A researcher who compares the README to the report
will find the report is the more confident of the two, and that is the exact
inversion this project has spent its credibility avoiding.
