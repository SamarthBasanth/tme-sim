# Quickstart

A fast, copy-pasteable path to your first result with **Vivum** (package name
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
pip install -e .
```

That installs the core simulator (numpy, scipy, matplotlib, pandas,
scikit-image, networkx, scikit-learn, requests, tifffile, imagecodecs) and
registers a `vivum` command on your `PATH`.

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
vivum --help
```

```
usage: vivum [-h] [--seed SEED]
             {reconstruct,simulate,compare-drugs,molecules,molecular-run,acquire,studio}
             ...
```

(`python -m tmesim.cli --help` is the exact same tool, if you prefer running
it as a module instead of the installed script.)

---

## 2. First example: the CLI

### 2a. Reconstruct a real tumor section and see basic stats

The repository ships an already-cached open **10x Xenium glioblastoma**
section (40,887 real, measured cells — no synthetic data anywhere in this
path). Reconstruct it and print a report:

```bash
vivum reconstruct --source xenium --dx 30
```

This prints cell-type fractions, mean vascular/ECM/oxygen density, the grid
size, and the per-cell drug-sensitivity assumption in force (documented
inline — reconstructing tissue is drug-agnostic, so this is always "uniform
1.0" at this stage), and saves a figure to `results/reconstruct_xenium.png`.
`--dx` is the grid spacing in µm; smaller is slower but higher-resolution —
30 µm finishes in a few seconds, the grid-converged production value used
throughout `results/` is 12 µm (see the README's "Five corrections" section
for why).

### 2b. Compare several drugs on that same section

```bash
vivum compare-drugs --dx 30 --drugs temozolomide,cisplatin,doxorubicin --csv results/my_comparison.csv
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
every difference in the table is attributable to that drug's own physics
(structure-derived diffusivity, real GDSC/literature IC50, BBB penetrance,
perivascular-sink transport calibrated once on doxorubicin) on identical
tissue — the same methodology `experiments/real_drug_panel.py` uses, just
wrapped in a CLI.

Drop `--drugs` and it compares **every** drug currently in the molecular
library instead of a fixed list:

```bash
vivum compare-drugs --dx 30
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
vivum compare-drugs --library breast --dx 30
```

A drug that fails to compile or simulate (for instance, one whose potency
data is still being wired up) is skipped with a printed reason rather than
aborting the whole comparison — useful precisely because these libraries
change over time.

Read the honesty caveats before quoting any number from this table: it is a
**prediction** of validated transport + potency physics on a real,
**untreated** tissue section, not a validated clinical outcome — there is no
paired measured drug response to check it against. See
[`results/real_drug_panel.md`](../results/real_drug_panel.md) and the
README's "Five corrections" section for what is and is not defensible to
quote.

### 2c. Launch Studio (the interactive web UI)

Needs the `app` extra installed (`pip install -e '.[app]'` above):

```bash
vivum studio
```

```
tmesim Studio backend
  repo        /path/to/tme-sim
  run history /path/to/tme-sim/app/server/.runs
  recon cache maxsize=4
  serving     http://127.0.0.1:8000
```

Then open `http://127.0.0.1:8000` in a browser. `vivum studio --port 8080` to
use a different port; `vivum studio --help` for the full option list. This is
a thin delegate to the existing `python -m app.server` entry point — nothing
about the server itself changed.

### 2d. TPL — a small scripting language for defining and sweeping experiments

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
#    drop dx_um, or set it to 12.0, for the grid-converged production value).
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

* **Data is pre-cached, not auto-downloaded on first use.** The Xenium /
  Visium sections these examples run on already live under `data_cache/` in
  this checkout. `vivum acquire` best-effort downloads the *smaller* open
  datasets; the large single-cell caches (Xenium GBM/breast, Visium, HGCC)
  are too big for that command and are fetched with the explicit `curl`
  commands in [`docs/REPRODUCE.md`](REPRODUCE.md) §2b. If you `git clone`
  this repository without its data cache, `reconstruct` / `simulate` /
  `compare-drugs` will raise a clear file-not-found error naming the missing
  path rather than silently substituting synthetic data — there is no
  synthetic fallback anywhere in this codebase.
* **`--sample` is a cache directory path, not a dataset name**, and its
  default depends on `--source`: for `xenium` it is the 319-gene GBM section;
  for `visium` it is `10x_parent_gbm`. Passing `--library breast` to
  `compare-drugs` without `--sample` now resolves the breast Xenium cache
  automatically; every other command needs `--sample <path>` spelled out for
  a non-default section (see `vivum reconstruct --help`).
* **The molecular drug libraries are moving targets.** `GBM_DRUGS` and
  `BREAST_DRUGS` are both under active extension. `vivum molecules` and
  `vivum compare-drugs` (no `--drugs`) always read the current list rather
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
| What is the current state, honestly? | [`results/SESSION_AUDIT.md`](../results/SESSION_AUDIT.md) |
| How do I regenerate every published result from scratch? | [`docs/REPRODUCE.md`](REPRODUCE.md) |
| Where did every byte of data come from? | [`docs/DATA_INVENTORY.md`](DATA_INVENTORY.md) |
| What is the model, mathematically? | [`docs/MATHEMATICAL_MODEL.md`](MATHEMATICAL_MODEL.md) |
| How do the pieces fit together? | [`docs/architecture.md`](architecture.md) |
| How should I compare drugs correctly? | [`results/delivery_metrics.md`](../results/delivery_metrics.md) |
