# Quickstart

A fast, copy-pasteable path to your first result with **Codon** (package name
`tmesim`) — the glioblastoma + breast-cancer tumor-microenvironment
drug-transport simulator. Every command and every code snippet on this page was
actually run against this repository before being written down; none of it is
aspirational.

This page assumes no prior familiarity with the codebase's internal module
layout. If you want the full scientific picture — what is real, what is
validated, and what is not — read [`README.md`](../README.md) and
[`results/SESSION_AUDIT.md`](../results/SESSION_AUDIT.md) next; this page is
only about getting the tool running.

---

## 1. Install

Python ≥ 3.9 (tested on CPython 3.9.6).

```bash
git clone <this-repo> tme-sim
cd tme-sim
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install --upgrade pip          # not optional — see below
pip install -e .
```

That installs the core simulator (numpy, scipy, matplotlib, pandas,
scikit-image, networkx, scikit-learn, requests, tifffile, imagecodecs) and
registers a `codon` command on your `PATH`.

**Why the `pip install --upgrade pip` line is there.** `pyproject.toml` has no
`setup.py`, so an editable install needs PEP 660 support, which arrived in pip
21.3. macOS 14's system Python 3.9.6 ships pip 21.2.4, and on that pip the line
above fails with:

```
ERROR: File "setup.py" or "setup.cfg" not found. Directory cannot be installed
in editable mode ... (A "pyproject.toml" file was found, but editable mode
currently requires a setuptools-based build.)
```

Measured on a clean clone, 2026-08-15. The upgrade takes about a second.

**If you want Studio** (the local web UI for exploring a reconstruction
interactively), install the additive `app` extra instead:

```bash
pip install -e '.[app]'
```

**Optional, for the full molecular layer** (real structure-derived molecular
descriptors instead of a fallback table) and some external-validation readers:

```bash
pip install rdkit openpyxl statsmodels
```

Nothing above is required for the two examples below — they work with just
`pip install -e .`. See [`docs/REPRODUCE.md`](REPRODUCE.md) §1 if you also
want the readers for `data_cache/hunt/`.

**Run everything from the repository root** (`cd tme-sim`, as above). The
tool reads real tissue data from `data_cache/` and writes figures/reports to
`results/`, both resolved relative to your current directory — exactly like
any other data-driven research CLI.

Confirm the install:

```bash
codon --help
```

```
usage: codon [-h] [--seed SEED]
             {reconstruct,simulate,compare-drugs,validate,ingest,molecules,molecular-run,acquire,sections,studio}
             ...
```

(`python -m tmesim.cli --help` is the exact same tool, if you prefer running
it as a module instead of the installed script.)

---

## 1b. Get the tissue (a clone has none)

`data_cache/` is gitignored. **A fresh clone contains no tissue at all**, so
every section reads `UNAVAILABLE` and `reconstruct` refuses with a named
missing directory. Confirm that for yourself, then fix it:

```bash
codon sections                 # 0 of 7 runnable on a clone, with the fix for each
codon acquire                  # fetches the tabular sources + the xenium_319 section
```

`acquire` downloads two different kinds of thing, and conflating them is what
used to make a first run dead-end:

| kind | what | makes a section runnable? |
|---|---|---|
| 4 **tabular** sources | TCGA-GBM clinical, TCGA-GBM expression, Ivy GAP, GDSC2 glioma potency | **no** — these parameterise *drugs* |
| 1 **spatial section** | `xenium_319`: 10x Xenium GBM, 40,887 measured cells, 16.6 MB, 3 files | **yes** — this is the tissue |

Useful variants:

```bash
codon acquire --list                 # print the plan, download nothing
codon acquire --sections all         # + the 480-panel IO section (287.6 MB)
codon acquire --sections none        # tabular sources only (the pre-2026-08 behaviour)
codon sections --verbose             # + how each recipe's URLs were verified
```

The Visium sections and `gse324268` are **not** auto-fetchable — they need
extraction, a format conversion or a companion dataset (the Darmanis
deconvolution reference), so no one-liner would be truthful about them.
`codon sections` prints the literal shell block for each, and
[`REPRODUCE.md`](REPRODUCE.md) §2b carries the same commands.

**The whole path, measured end to end on a clean clone** (macOS 14, CPython
3.9.6, `pip` cache warm, 2026-08-15):

| step | wall clock |
|---|---:|
| `git clone` (from a local path) | 2 s |
| `python3 -m venv .venv` | 1 s |
| `pip install --upgrade pip` | 1 s |
| `pip install -e .` | 13 s |
| `codon sections` → 0 of 7 runnable, exit 1 | 9 s |
| `codon acquire` → 4 of 5 sources ok, 1 of 7 sections runnable | 3 m 52 s |
| `codon simulate --source xenium --dx 30 --drug temozolomide` | 35 s |
| **clone → first completed simulation** | **4 m 53 s** |

`acquire` dominates that, and its time is dominated by the GDC endpoints: three
timed runs of the same command on the same day took **1 m 50 s, 3 m 40 s and
3 m 52 s**. Treat it as minutes, not as a number to plan around. The `pip` steps
were measured with a warm pip cache; a first-ever install downloads wheels and
takes longer.

In the run tabulated above the Ivy GAP host returned **HTTP 504** and that
source failed. `acquire` reported it as `failed` on its own line and carried on,
the section still landed, and the simulation still ran — a tabular source going
down does not cost you your first result. That is the intended behaviour, not a
lucky escape.

If you have your own de-identified export, you do not need any of this: register
it by path with `codon ingest /path/to/export --register` and select it the same
way. Nothing is uploaded anywhere — the directory is read where it lies.

---

## 2. First example: the CLI

### 2a. Reconstruct a real tumor section and see basic stats

`codon acquire` has now cached the open **10x Xenium glioblastoma** section
(40,887 real, measured cells — no synthetic data anywhere in this path).
Reconstruct it and print a report:

```bash
codon reconstruct --source xenium --dx 30
```

This prints cell-type fractions, mean vascular/ECM/oxygen density, the grid
size, and the per-cell drug-sensitivity assumption in force (documented
inline — reconstructing tissue is drug-agnostic, so this is always "uniform
1.0" at this stage), and saves a figure to `results/reconstruct_xenium.png`.
`--dx` is the grid spacing in µm; smaller is slower but higher-resolution —
30 µm finishes in a few seconds, the production value used throughout `results/`
is 12 µm. **Do not read that as "converged"** — the claim *"dx = 12 µm is
grid-converged"* is retracted (README "What is retracted";
`results/grid_convergence_gbm.md` §0, `results/grid_convergence_breast.md`).
Neither panel is in an asymptotic regime; 12 µm is a production choice.

### 2b. Compare several drugs on that same section

```bash
codon compare-drugs --dx 30 --drugs temozolomide,cisplatin,doxorubicin --csv results/my_comparison.csv
```

```
=== drug comparison (sorted by predicted kill_fraction, best first) ===
PREDICTIONS on a real reconstructed section, not validated clinical outcomes — see docs/QUICKSTART.md.
drug                   MW  IC50(uM)    bbb    D_eff  pen_idx   mean_auc    kill
-------------------------------------------------------------------------------
cisplatin             298      24.6  0.858    236.2    0.375   2.02e+05  0.4174
doxorubicin           544       0.5  0.021    184.7    0.225   1.34e+03  0.2573
temozolomide          194       606  0.125    278.0    0.472   7.48e+04  0.0210
```

The tissue section is reconstructed **once** and reused for every drug, so
every difference in the table is attributable to that drug's own parameters
(diffusivity, real GDSC/literature IC50, BBB penetrance, perivascular-sink
transport calibrated once on doxorubicin) on identical tissue — the same
methodology `experiments/real_drug_panel.py` uses, just wrapped in a CLI.

> **The `D_eff` column comes from a two-route transport model, and it recovers
> about a quarter of the measured compound-to-compound spread.** `logP` and
> hydrogen-bond donor count now reach transport through a membrane permeability.
> Scored blind and unfitted against 71 compounds measured in one cell line
> (Pruijn 2008): **R² = 0.718**, 99 % within 3× of measurement, **24 of 31**
> exact-mass pairs ordered correctly. The molecular-weight-only model this
> replaced was FALSIFIED on the same data — R² = 0.022, flat by construction at
> matched mass (exactly 1.00× against a measured 3.7×), **0 of 31** pairs
> ordered. **What is still not true:** it reproduces 8.3× of the measured 29.7×
> spread, so it *under-separates* similar-mass analogues rather than resolving
> them, and the whole test is one cell line. **Read a ranking over such a series
> as a hypothesis worth testing, not as a transport measurement.** Ranking
> compounds that differ substantially in mass remains the best-supported case.
> → [`../results/two_route_promotion.md`](../results/two_route_promotion.md)
> → [`../results/diffusivity_falsification.md`](../results/diffusivity_falsification.md),
> [`../results/two_route_transport.md`](../results/two_route_transport.md)
> (the two-route fix reaches R² = 0.725 blind, but is **default-OFF and wired
> into nothing**, so no number here benefits from it).

Drop `--drugs` and it compares **every** drug currently in the molecular
library instead of a fixed list:

```bash
codon compare-drugs --dx 30
```

(this runs one simulation per library drug on the same reconstructed section
— about 2 minutes for the 11 drugs in the library at the time of writing; add
`--drugs` to run only the ones you care about, as in 2b, if you want it
faster.)

`tmesim.molecular.library.GBM_DRUGS` is actively growing (new agents are
being added on an ongoing basis), which is exactly why `compare-drugs` reads
the list at run time rather than the CLI hardcoding today's names — a newly
added drug is comparable with no CLI change. For breast tissue, pass
`--library breast`, which draws its default list from
`tmesim.molecular.breast_library.BREAST_DRUGS` and reconstructs the open
Janesick et al. 2023 breast Xenium section instead:

```bash
codon compare-drugs --library breast --dx 30
```

A drug that fails to compile or simulate (for instance, one whose potency
data is still being wired up) is skipped with a printed reason rather than
aborting the whole comparison — useful precisely because these libraries
change over time.

Read the honesty caveats before quoting any number from this table: it is a
**prediction** on a real, **untreated** tissue section, not a validated clinical
outcome — there is no paired measured drug response to check it against. The
transport physics is **not** validated either: its input coefficients are
externally testable and the engine currently misses them (see the `D_eff` box
above, and `codon validate`). What is defensible here is the **order**, not the
level, and not over a similar-mass series. See
[`results/real_drug_panel.md`](../results/real_drug_panel.md) and the
README's "Five corrections" section for what is and is not defensible to
quote.

### 2c. Ask how the model currently scores against reality

```bash
codon validate --quick
```

This runs the **external-anchor suite**: every place in this repository where a
model-produced number is scored against an *independent measured* number from
the published literature. One row per comparison, each citing its PMID and
stating exactly what quantity is being compared against what.

`--quick` runs only the analytic anchors — closed-form penetration lengths,
descriptor-derived partition coefficients, literature-table comparisons — which
need no tissue section and finish in a few seconds. Drop it to run the full
suite, which reconstructs real sections from `data_cache/` and solves the
transport PDE (minutes, and it is the only mode that tests the solver at all).

**This command exits 1, and that is the correct result.** Anchors that FAIL are
published limitations of the model — they fail identically on a clean clone, on
every machine — and the scorecard prints a block saying exactly that before the
failing rows. An exit code of 1 here means the scorecard ran and reported
honestly; it is not a broken install, and a run in which those rows passed would
be the thing to distrust.

```
NOT TESTABLE analytic doxorubicin   [PMID 16361566]
                      compares  : analytic penetration length L = sqrt(D_eff/k_eff) vs measured
                                  perivascular half-distance
                      predicted : L = 64.9 um
                      measured  : 40-50 um half-distance (= 58-72 um e-fold, d_half/ln2)
                      bar       : n/a — this IS the fit
                      NO VERDICT: IN-SAMPLE. alpha is solved in closed form so that L equals this target
                                  exactly; the model reproduces it by construction and could not have
                                  failed. Reporting it as a PASS would be circular.

PASS         analytic topotecan   [PMID 23680920]
                      compares  : OUT-OF-SAMPLE analytic penetration length L = sqrt(D_eff/k_eff),
                                  alpha frozen, vs measured perivascular half-distance
                      predicted : L = 83.6 um  (1.16x band mid)
                      measured  : 25-75 um half-distance (= 36-108 um e-fold, d_half/ln2)
                      bar       : inside the published e-fold band, else within 3x of its midpoint
                      note      : camptothecin topoisomerase-I inhibitor (chemically DISTINCT from the
                                  anthracycline fit -- a more independent test)
                      conditional on: alpha-fit, analytic-only

FAIL         analytic adavosertib   [PMID 39820282]
                      compares  : predicted tissue:plasma partition vs MEASURED human glioma
                                  tissue:plasma AUC ratio
                      predicted : 0.5785  (0.24x measured)
                      measured  : 2.4600 (range 0.7400-4.1800)
                      bar       : fold error within 3x (module BAND_TIGHT); 10x reported alongside
```

At the time of writing `--quick` covers 8 of the 16 registered anchor groups and
ends on **PASS 7 / FAIL 6 / NOT TESTABLE 4 / DEFERRED 8** of 25 rows, exiting 1.
The full run adds the PDE and geometry anchors and ends deeper in the red — see
the failure list below. These counts move as anchors land; regenerate them from
live output rather than trusting this page.

**Four statuses, and the distinction between them is the point.**

| status | meaning |
|---|---|
| `PASS` | the model ran and landed inside the band the owning module declares (the bar is named on every row — different modules declare different bars) |
| `FAIL` | the model ran and missed. Reported three times over: in the row, again in a dedicated `FAILURES` block with full detail, and in the summary |
| `NOT TESTABLE` | **no verdict is possible**, for a stated reason — and this is never a soft pass |
| `DEFERRED` | skipped by `--quick`. Also not a pass: some anchors that are DEFERRED here **FAIL** in the full run, which is why `--quick` prints a banner saying so |

`NOT TESTABLE` is load-bearing and occurs for six genuinely different reasons,
each printed inline:

* **the number is the *fit*, not a test of it** — doxorubicin's penetration
  anchor is solved for in closed form, so the model reproduces it by
  construction and could not have failed;
* **no measured target exists** — nobody has published a free-drug macrophage
  uptake ratio, so the immune sink's magnitude has nothing to be scored
  against, and the same is true of every breast intercapillary-distance row
  (GBM has Yoshii & Sugiyama; breast provably does not);
* **the two sides are not the same quantity** — a measured *compartment*
  contrast is not the model's tile-to-tile spread;
* **the prediction is threshold-chosen** — the KDE vascular area fraction lands
  wherever the chosen threshold puts it, so agreeing with the measurement would
  be a choice rather than a result;
* **the model cannot express the measured quantity** — the `bbb_penetrance`
  partition-vs-rate units defect;
* **required data absent** — the anchor needs bytes under `data_cache/` that
  this checkout does not have.

Conflating any of those with `PASS` is exactly the error this project retracts
results for.

**Failures are the reason the command exists.** The full run currently ends on
**PASS 13 / FAIL 18 / NOT TESTABLE 19** across 50 rows — more failures than
passes, and all of it printed rather than buried:

```
====================================================================================================
FAILURES (18) -- the model misses these MEASURED values
====================================================================================================

  FAIL  breast-transvascular / 5-fluorouracil   [PMID 9205062]
        compares  : mean_tumor_auc / (plasma_mean * total_time) — free interstitial drug AUC at
                    tumour cells over plasma AUC — vs MEASURED in-vivo microdialysis
                    interstitial:plasma AUC in primary human breast tumour
        predicted : 7.7956  (12.78x measured)
        measured  : 0.61 +/- 0.11 (mean +/- SE, n=10)
        physics   : transvascular_exchange=False — one-way leak, S = PS*Cp
        bar       : inside measured mean +/- SE

  FAIL  breast-transvascular / methotrexate   [PMID 9679959]
        predicted : 0.0571  (0.10x measured)
        measured  : 0.60 +/- 0.20 (mean +/- SE, n=9)
        physics   : transvascular_exchange=False — one-way leak, S = PS*Cp
```

The 18 span every layer: the two breast microdialysis anchors, two of six
glioma partition drugs, the `bbb_penetrance` rank ordering, the GIST imatinib
falsification, the two MCL transport-coefficient rows, **one** of the breast
3-D:2-D kill ratios (docetaxel), all four GBM vascular-morphometry anchors,
breast microvessel density against Weidner, and section cellularity against
measured H&E nuclear counts.

The breast 3-D:2-D panel is **3 PASS / 1 FAIL / 6 NOT TESTABLE** — only four of
ten drugs are testable at all. Two earlier counts, "5 PASS / 1 FAIL" and "3 PASS
/ 3 FAIL", were **both wrong**: the second scored two drugs against **zero-width**
bands that no prediction can land inside. A single-study band is now NOT TESTABLE
in both the tool and the document generator.
→ [`../results/zero_width_bands.md`](../results/zero_width_bands.md)

`codon validate` **exits non-zero if any anchor fails**, so it drops straight
into CI:

```bash
codon validate --quick || echo "an external anchor regressed"
```

Note that `--quick` is a *weaker* gate than the full run, not a faster
equivalent: **13 of the 13 rows it scores never invoke the transport PDE**, so a
green `--quick` is not evidence about the solver.

#### Every row names the physics it ran under

The breast interstitial:plasma anchor moves by roughly **10x** depending on
which transvascular delivery law is active — 5-FU lands at 12.78x of measured
under the shipping one-way leak `S = PS*Cp`, and near 1.33x under the opt-in
gradient-driven `S = PS*(Cp - C)`. A bare number for that anchor is therefore
ambiguous by an order of magnitude, so the scorecard prints the law it ran
under at the top of the run *and* on every row:

```
PHYSICS   transvascular_exchange=False — one-way leak, S = PS*Cp (source never depends
          on local tissue concentration)
          This IS the law the engine ships with (PDEConfig.transvascular_exchange).

FAIL         PDE      5-fluorouracil   [PMID 9205062]
                      predicted : 7.7956  (12.78x measured)
                      measured  : 0.61 +/- 0.11 (mean +/- SE, n=10)
                      physics   : transvascular_exchange=False — one-way leak, S = PS*Cp
```

The rules behind that line:

* The shipping default is **read from `PDEConfig`**, never hardcoded, so a
  change to what the engine ships cannot silently desynchronise from what the
  scorecard claims the default is.
* Running `--exchange` prints a loud `!! NOT THE SHIPPING DEFAULT` banner and
  appends the same marker to every affected row, because those numbers are not
  the ones the shipping engine produces — and each is additionally conditional
  on the uncalibrated `PS_SCALE_S_INV`.
* The breast anchor reads the law back **off the config that was actually
  solved**, not off the flag it was asked for. A mismatch is reported as a
  `FAIL` provenance row rather than scored.
* A row whose `scope` is `PDE` and that carries a verdict but declares **no**
  law is downgraded to `NOT TESTABLE` with `PHYSICS PROVENANCE MISSING`. An
  unlabelled solver number is not a weaker result, it is an uninterpretable
  one, and it must not be allowed to pass quietly.
* Rows that never invoked the solver read `n/a — transport solver never invoked
  for this row`. That says the question does not arise, **not** that the two
  laws agreed.
* `--json` carries the same provenance twice: a top-level `physics` object for
  the run, and a `physics` string on every row — so two stored scorecards
  cannot be diffed without the delivery law being visible on both sides.

Three more things worth knowing before you quote a row:

* **Rows carry a `scope`.** `analytic` means the prediction is a closed-form
  expression over molecular descriptors and the transport solver was never
  invoked — such a row scores the molecular layer and *cannot* be cited as
  evidence about the PDE. `geometry` means the row scores the reconstructed
  *tissue* (vessel density, spacing, area fraction, cellularity) with no drug
  and no solver anywhere in it — a real claim, but about what the transport
  runs **on**, not about the transport. The summary line prints how many of the
  scored rows never touch the solver. Only the full (non-`--quick`) run
  produces `PDE` rows at all.
* **Rows list the conditionals they inherit**, and every cited conditional is
  printed in full at the bottom. The two with the longest reach are `alpha-fit`
  (everything downstream of the one fitted constant) and `ps-scale` (every
  transvascular-exchange number is conditional on
  `tmesim/treatment/drugs.py::PS_SCALE_S_INV = 1.0`, an implicit and underived
  dimensionless-to-rate conversion — see
  [`results/transvascular_exchange_revalidation.md`](../results/transvascular_exchange_revalidation.md)
  §6.2).
* **The omissions are auditable too.** `codon validate --list` prints the
  registry *and* the validations that were inventoried and deliberately left
  out, with the reason for each (ablations that would double-count the same six
  measurements; real-vs-real corroboration with no model in the loop; the
  `vessel_profile` anchor that `alpha` was itself fitted to).

Useful options:

```bash
codon validate --list                          # what is checked, and what is not
codon validate --quick --json                  # machine-readable, for CI
codon validate --only breast-transvascular     # one anchor
codon validate --only vascular-morphometry,breast-vascular-morphometry  # several
codon validate --only breast-transvascular --exchange   # with the opt-in delivery law
```

`--json` emits the whole scorecard — every row's predicted/measured values,
PMIDs, bar, status, reason and raw numbers — as a single parseable document
with no log lines around it. The registry itself lives in
`tmesim/validation/scorecard.py`; adding an anchor is one entry there and needs
no CLI change.

### 2d. Launch Studio (the interactive web UI)

Needs the `app` extra installed (`pip install -e '.[app]'` above):

```bash
codon studio
```

```
Codon Studio backend
  repo        /path/to/tme-sim
  run history /path/to/tme-sim/app/server/.runs
  recon cache maxsize=4
  serving     http://127.0.0.1:8000
  your tissue http://127.0.0.1:8000/ingest
```

Then open `http://127.0.0.1:8000` in a browser. `codon studio --port 8080` to
use a different port; `codon studio --help` for the full option list. This is
a thin delegate to the existing `python -m app.server` entry point — nothing
about the server itself changed.

**To load YOUR OWN tissue, open `/ingest`.** Studio's section picker currently
lists only the sections this project prepared; `/ingest` is a page served by the
backend that takes the path to your own 10x export, prints the same capability
report `codon ingest` does — including what your panel *cannot* support — and
registers the sample so it appears in the picker on reload. It is the same
`POST /api/sections/scan` endpoint, without the `curl` — and `codon ingest
<dir> --register` does the same thing from the command line.
[`results/rebrand_backend.md`](../results/rebrand_backend.md) §5 lists the
front-end work that would fold this into the picker itself.

**De-identified tissue only.** Codon is not built for identifiable patient data
and must not be used with it. The sample is read from a directory on the
server's own machine and never uploaded, but the absolute path, the section id
and the label are all served by `GET /api/sections` and shown in the picker, so
identifiers hidden in a directory name leak just as readily as ones in a column.
[`DEPLOYMENT.md` §0](DEPLOYMENT.md) lists the four places to check first.

### 2e. TPL — a small scripting language for defining and sweeping experiments

Studio's **"protocol"** tab (top nav) is a real code editor, not just a form:
syntax highlighting, line-numbered diagnostics, `Cmd+Enter` to plan+run,
`Cmd+S` to validate, `F8`/`Shift+F8` to jump between warnings. Its own spec
describes the design stance directly: **"MATLAB-adjacent in feel — assignment,
ranges, loops over parameters, an interactive run — but domain-specific... no
arbitrary control flow, no user functions, no I/O. Every construct maps onto a
real quantity in the model."** Full grammar and reference:
[`PROTOCOL_LANGUAGE.md`](../PROTOCOL_LANGUAGE.md).

A complete program — define a drug, define a dosing protocol, sweep a
parameter, run it, assert a hypothesis:

```
drug TMZ from library "temozolomide"
  ic50 from gdsc                  # measured; errors if TMZ is absent from the screen
end

protocol single_dose
  route oral
  dose 150 mg/m2
  every 24 h
  repeat 1
end

sweep dx from 12 to 60 um step 12  # 5 runs, one axis

run TMZ with single_dose on xenium_319
  sensitivity uniform              # NO-RESISTANCE UPPER BOUND, not an expected response
  report coverage, kill, mean_auc, penetration_length
end

expect coverage > 0.2              # recorded, never enforced — state the hypothesis first
```

Paste that into the editor (or start from one of the 6 built-in templates —
"single drug", "dose sweep", "Stupp-like multi-cycle", and others) and hit
`Cmd+Enter`. `plan` shows the expanded run count and estimated cost and waits
for explicit confirmation before anything queues — nothing runs by accident.
Every refusal is a parse error with a line and column, never a silent
fallback; the same caveats this project attaches to every `results/*.md`
report (e.g. "ABSOLUTE KILL MUST NOT BE QUOTED FROM this section under any
state route") appear as inline diagnostics at the exact line that would
trigger them, so they're visible at the point of use, not just in a
separately-read report.

The same three routes the editor uses are also a plain HTTP API, if you want
to script protocol runs without the browser:

```bash
curl -s -X POST http://127.0.0.1:8000/api/protocol/validate \
  -H "Content-Type: application/json" \
  -d '{"source": "drug TMZ from library \"temozolomide\"\n  ic50 from gdsc\nend\n\nprotocol single_dose\n  route oral\n  dose 150 mg/m2\n  every 24 h\n  repeat 1\nend\n\nrun TMZ with single_dose on xenium_319\n  dx 12 um\n  sensitivity uniform\n  report coverage, kill\nend\n"}'
```

`/api/protocol/plan` (expand the sweep, show run count + cost, launch
nothing) and `/api/protocol/run` (submit for real) follow the same shape —
see `app/server/protocol_api.py`'s module docstring for the full route table.

---

## 3. Second example: the Python API (for a notebook)

Every CLI command above is itself a thin wrapper around a small, stable set
of top-level functions. If you want to script against the simulator directly
— in a Jupyter notebook, in your own analysis pipeline — this is the same
real pipeline, called directly:

```python
from tmesim.reconstruct.real_pipeline import reconstruct_real, reconstruction_report
from tmesim.validation.calibration import calibrate_alpha, calibrated_drugparams
from tmesim.molecular.library import get_molecule
from tmesim.molecular.binding_site_barrier import apply_barrier_to_drugparams
from tmesim.core import SimConfig
from tmesim.sim.engine import simulate_response

# 1. Reconstruct a real tumor section (coarse grid here for a fast demo;
#    drop dx_um, or set it to 12.0, for the production value — production, NOT
#    converged; see the retraction note above and results/grid_convergence_gbm.md).
micro = reconstruct_real(source="xenium", dx_um=80.0)
report = reconstruction_report(micro)
print("cells:", report["n_cells"], "tumor:", report["n_tumor"], "grid:", report["grid"])
# cells: 40887 tumor: 23812 grid: (65, 103)

# 2. Compile a real, calibrated drug: structure -> Stokes-Einstein diffusivity
#    -> perivascular-sink calibration (alpha, fit once on doxorubicin) ->
#    binding-site barrier. Any name in tmesim.molecular.library.GBM_DRUGS works.
alpha = calibrate_alpha()
drug_params = calibrated_drugparams("temozolomide", alpha)
mol = get_molecule("temozolomide")
drug_params = apply_barrier_to_drugparams(drug_params, mol)
print("IC50 (uM):", drug_params.ic50, "BBB penetrance:", drug_params.bbb_penetrance)
# IC50 (uM): 606.0466722391783 BBB penetrance: 0.12537158358322562

# 3. Simulate transport + kill on the reconstructed section.
cfg = SimConfig(drug=drug_params, seed=0)
cfg.pde.total_time = 1800.0   # fast demo window (seconds); library default is 12 h
result = simulate_response(micro, cfg)
print("kill_fraction:", round(result.kill_fraction, 4))
# kill_fraction: 0.0198
```

`result` (a `tmesim.sim.engine.SimulationResult`) also carries
`result.viable_tumor_fraction`, `result.transport` (the raw exposure field),
`result.fields` (oxygen / IFP / AUC / vascular density, for your own
plotting), and `result.extras` (`penetration_index`, `mean_auc`,
`hypoxic_tumor_fraction`, and more) — the same object every experiment script
under `experiments/` and every `results/*.md` report is built from.

`reconstruct_real` and `reconstruction_report` are documented in
`tmesim/reconstruct/real_pipeline.py`; `calibrated_drugparams` and
`calibrate_alpha` in `tmesim/validation/calibration.py`; `simulate_response`
(including what each pipeline stage does and what `attribute=True` /
`capture_frames=` add) in `tmesim/sim/engine.py`. Each docstring is the
authoritative reference — this page only gets you to the point of being able
to read them productively.

For a breast-tissue reconstruction from the Python API, point `sample=` at
the cached Janesick et al. 2023 section and compile through
`tmesim.molecular.breast_library` for breast-specific drugs:

```python
micro = reconstruct_real(
    source="xenium", sample="data_cache/xenium_breast_janesick/rep1", dx_um=80.0,
)
```

---

## 4. Known friction points (read before filing an issue)

* **A clone has no data, and nothing downloads itself on first use.**
  `data_cache/` is gitignored, so a fresh clone has zero runnable sections.
  `codon acquire` (§1b) fetches the tabular sources *and* the `xenium_319`
  section, which is the one every example here uses; `--sections all` adds the
  other auto-fetchable section. The Visium sections, `gse324268`, the breast
  section and HGCC are **not** auto-fetchable and `codon sections` prints the
  exact commands for each. Until then, `reconstruct` / `simulate` /
  `compare-drugs` raise a clear file-not-found error naming the missing path
  rather than silently substituting synthetic data — there is no synthetic
  fallback anywhere in this codebase.
* **`codon validate` exits 1 on a healthy install.** Several external anchors
  FAIL against measured data on every machine, deliberately; they are this
  project's published limitations, not a symptom of a bad setup, and the
  command prints a block saying so above the failing rows. Do not "fix" it.
  `codon validate --list` is the registry of what is scored and what was
  deliberately left out.
* **`--sample` is a cache directory path, not a dataset name**, and its
  default depends on `--source`: for `xenium` it is the 319-gene GBM section;
  for `visium` it is `10x_parent_gbm`. Passing `--library breast` to
  `compare-drugs` without `--sample` now resolves the breast Xenium cache
  automatically; every other command needs `--sample <path>` spelled out for
  a non-default section (see `codon reconstruct --help`).
* **The molecular drug libraries are moving targets.** `GBM_DRUGS` and
  `BREAST_DRUGS` are both under active extension. `codon molecules` and
  `codon compare-drugs` (no `--drugs`) always read the current list rather
  than a value baked into this document — trust the command's own output over
  any specific drug name mentioned here.
* **Absolute kill fractions are predictions, not validated outcomes.** Every
  open spatial section this simulator ships with is an *untreated* resection
  — there is no paired measured drug response to check a predicted kill
  fraction against. Read `results/SESSION_AUDIT.md` before quoting one.
* **Run commands from the repository root.** Data (`data_cache/`) and output
  (`results/`) paths are resolved relative to your current working directory.

---

## Where to look next

| Question | File |
|---|---|
| How does the model score against measured reality, right now? | `codon validate` (§2c), then [`results/SESSION_AUDIT.md`](../results/SESSION_AUDIT.md) |
| What is the current state, honestly? | [`results/SESSION_AUDIT.md`](../results/SESSION_AUDIT.md) |
| How do I regenerate every published result from scratch? | [`docs/REPRODUCE.md`](REPRODUCE.md) |
| Where did every byte of data come from? | [`docs/DATA_INVENTORY.md`](DATA_INVENTORY.md) |
| What is the model, mathematically? | [`docs/MATHEMATICAL_MODEL.md`](MATHEMATICAL_MODEL.md) |
| How do the pieces fit together? | [`docs/architecture.md`](architecture.md) |
| How should I compare drugs correctly? | [`results/delivery_metrics.md`](../results/delivery_metrics.md) |
| My whole team needs to reach Studio, not just my machine | [`docs/DEPLOYMENT.md`](DEPLOYMENT.md) |
