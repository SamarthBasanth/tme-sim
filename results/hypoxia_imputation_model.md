# ML-imputed hypoxia state — honest cross-section validation

*2026-07-28. This is the hypoxia-specific deliverable for the ML microenvironment-
reconstruction task. It was subsequently folded into, and cross-checked
against, a wider 19-target sweep across hypoxia / MGMT / proliferation / ECM /
drug efflux / vasculature — see `results/panel_gap_imputation_sweep.md` for the
full picture and the numbers this file's design choices are validated against.
This file stays focused on the one target and the one production route this
session actually wired: an opt-in ML-imputed HYPOXIC call for panels missing
the real hypoxia markers.*

## Verdict

**Real, honest, moderate signal — real enough to wire as an explicit opt-in
route, not real enough to trust at face value.** A `RandomForestRegressor`
trained on 3 real GBM-Space Visium sections predicts a continuous hypoxia
score from the 319-gene Xenium Human Brain panel's OTHER genes with:

| check | R² | Pearson r | what it measures |
|---|---|---|---|
| leave-one-section-out CV (GBM-Space, Visium, same assay) | **0.692** | **0.832** | can the model recover hypoxia score on a WHOLE tissue block it never trained on |
| cross-panel: independent Xenium single-cell section (`gse_366`, 4/5 hypoxia genes as ground truth) | **0.238** | **0.557** | can the SAME model, applied to a different assay/resolution/patient, still say something real |
| cross-panel: second independent Xenium single-cell section (`io_480`, VEGFA-only ground truth) | 0.051 | 0.270 | the honest floor — weak, but not the number that should gate the decision (only 60 of 319 feature genes overlap this panel, and ground truth is a single gene) |

**Read the degradation as the finding, not as noise to explain away.** Going
from R²=0.69 (same assay) to R²=0.24 (cross-assay, richer ground truth) to
R²=0.05 (cross-assay, thinner ground truth and thinner feature overlap) is
exactly the pattern expected when spot-pooled Visium co-expression is asked
to generalise to single-cell Xenium counts of genes measured by a completely
different probe chemistry. It is real signal, genuinely degraded by a real
domain shift, not a fluke in either direction.

## Exactly what was trained on, and what was held out

* **Target**: `cell_state.HYPOXIA_MARKERS` = VEGFA, CA9, ADM, NDRG1, SLC2A1
  (the SAME gene set `reconstruct.cell_state.py`'s own measured-marker
  HYPOXIC call uses — not a redefinition) pooled into one continuous score:
  the mean, per cell/spot, of each present gene's log1p count z-scored across
  that section's own population.
* **Training / primary validation — GBM-Space Visium**
  (`data_cache/gbmspace_visium/{AT10-BRA-5-FO-2_1,AT10-BRA-5-FO-3_1,AT10-BRA-5-FO-4_2}`),
  whole transcriptome (36,591 genes), 3 real, distinct GBM tissue blocks
  (De Jong et al., bioRxiv 2025). **3-fold leave-one-section-out**: each fold
  trains on 2 of the 3 sections and tests on the ENTIRE third, held out
  completely — never a split of cells/spots within one section. Features:
  the 319-gene Xenium Human Brain panel's own gene list, minus the 5 hypoxia
  genes (319 features; confirmed all present in GBM-Space's whole
  transcriptome).
* **Held out entirely from training — real single-cell Xenium sections.**
  `gse324268_xenium` (366-gene panel, carries VEGFA/CA9/ADM/SLC2A1 — 4 of the
  5 target genes, NDRG1 absent even here) and `xenium_gbm_io` (480-gene IO
  panel, carries VEGFA only). A model trained on ALL 3 GBM-Space sections
  pooled is scored against each section's own real measured hypoxia genes,
  using only the subset of the 319 feature genes that section also measures
  (269 of 319 for `gse_366`; 60 of 319 for `io_480`).
* **Model**: `sklearn.ensemble.RandomForestRegressor(n_estimators=200,
  max_depth=10, min_samples_leaf=4, random_state=0)` — declared, not tuned
  against these numbers. Every fold beat the best single-gene cheap baseline
  (fit on training sections only, scored unselected on the held-out one); the
  best single gene was `HILPDA` — a real, independently-known direct HIF-1
  target (`real_pipeline.HIF_TARGET_GENES`) — in 2 of the 3 folds, which is
  independent evidence the forest is picking up real hypoxia biology and not
  fitting noise.

## Why this specific target, this specific model

Hypoxia state gates the largest absolute-kill term this project reports
(`cell_state.py`'s own module docstring: ~89% of the temozolomide kill
ceiling on the 319-gene panel comes from the QUIESCENT/HYPOXIC banding
alone), and the 319-gene panel — the one this project's flagship
reconstruction actually runs on — measures NONE of the 5 canonical hypoxia
genes, so that state is currently `UNMEASURED` on it by construction
(`cell_state.py`'s `state_from_markers`). A RandomForestRegressor (over
gradient boosting: a coin's-width choice either way; over any deep-learning
approach: 3 independent tissue sections is not enough data to fit or honestly
validate a network, and PyTorch is not an installed dependency of this
project) was trained specifically because the honest validation protocol —
section-held-out, then cross-assay — needed a model that would not silently
overfit 3 sections' worth of data.

## How a caller opts in

```python
from tmesim.reconstruct.real_pipeline import reconstruct_real, install_ml_imputed_hypoxia
from tmesim.data.real.singlecell_spatial import load

micro = reconstruct_real("xenium", state_mode="markers")   # unchanged default route
data = load(path="data_cache/xenium_gbm", dense=False)

result = install_ml_imputed_hypoxia(micro, data["expression"], data["genes"])
```

* **Not a `state_mode=` value.** `cell_state.STATE_MODES` is unchanged
  (`("oxygen", "vessel_distance", "markers")`) — `reconstruct_real(...,
  state_mode="ml_imputed")` raises `ValueError` exactly as any unrecognised
  string would (`tests/test_panel_imputation.py::test_reconstruct_real_refuses_ml_imputed_state_mode`).
  `STATE_MODES` is mirrored and asserted equal across three files
  (`cell_state.py`, `app/server/cache.py`, `tmesim/protocol/ast.py`);
  extending that shared DSL/server surface for a route validated on 3 Visium
  sections was judged out of scope for this session — `install_ml_imputed_hypoxia`
  is a function a caller applies explicitly to an already-built
  `Microenvironment` instead.
* **Refuses on any panel this sweep did not validate.** Only the 319-gene
  brain panel (`brain_319`) is accepted; any other panel raises
  `ValueError` naming the reason
  (`tests/test_panel_imputation.py::test_ml_imputed_refuses_an_unvalidated_panel`).
* **Never silently overwrites a measured call.** If the panel already
  measures a hypoxia marker, `cells.state` is left exactly as the measured
  route produced it, and the function returns `applied=False` with an
  explanation — the one guarantee every other test in this file exists to
  keep true
  (`tests/test_panel_imputation.py::test_ml_imputed_never_overwrites_a_measured_hypoxia_call`).
* **Provenance tag `kind='ml_imputed'`** — distinct from `measured` /
  `computed` / `literature` / `assumed` (`app/web/js/ui/provenance.js`'s
  taxonomy, deliberately not reused) — carrying this file's own R²/r numbers
  in the returned `caveat` string, so nothing downstream can quote the
  resulting HYPOXIC fraction without also seeing exactly how much to trust
  it.
* **Threshold, not a calibrated cutoff.** A cell is called HYPOXIC iff its
  predicted score exceeds the SECTION'S OWN mean — the identical "enrichment
  above the section rate" convention `cell_state.neighbourhood_positive_fraction`
  already uses for the exact same reason: no absolute cutoff here is a
  measurement either.

## If this doesn't clear your bar

It didn't clear a strict one: R²=0.238 on the more relevant cross-panel check
means roughly a quarter of the variance in a real section's hypoxia score is
recovered from a completely different assay's co-expression structure — real,
worth having as a labelled, opt-in estimate instead of a bare "UNMEASURED",
but nowhere near a substitute for actually measuring VEGFA/CA9/ADM/NDRG1/
SLC2A1 on the tissue. What would move this number: more real, distinct GBM
tissue sections in training (3 is the minimum defensible number for a
leave-one-section-out claim, not a generous one); a matched single-cell
training resource (this entire result is trained on spot-pooled Visium and
tested on single-cell Xenium — a same-resolution training set would likely
close much of the 0.69 → 0.24 gap); and, per §4.2 of the wider sweep, treating
individual hypoxia genes separately from the composite, since CA9 — the most
biologically SPECIFIC of the five — is also the worst individually predicted
(R²=0.094 same-assay), a genuinely different question from "does the pooled
composite work."

## See also

`results/panel_gap_imputation_sweep.md` — the same protocol applied to 18
further targets (MGMT, TOP2A/MKI67, 11 ECM genes, ABCB1/ABCG2, PECAM1 and 8
further vascular genes, RGS5/PDGFRB/NOTCH3/ACTA2), including two clean
negative results (MGMT, POSTN) and the strongest result in the whole sweep
(`vascular_composite`, R²=0.652 cross-panel) — reported with the same
discipline and not wired past the generic `impute_gene_score()` API because
neither has a `CellState`/`CellType` slot to occupy in this project's schema.
