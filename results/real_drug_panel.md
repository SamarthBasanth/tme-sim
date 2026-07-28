# Real drug panel on the real Xenium glioblastoma microenvironment

> ## ⚠ STALE — the absolute numbers in this file no longer match the engine
>
> **Do not quote an absolute kill, AUC or coverage from this file.** The engine has moved
> since it was generated. The *ranking* is unaffected and still reproduces exactly.
>
> | | this file | current engine |
> |---|---|---|
> | tumour cells typed | 30,174 | **23,812** |
> | mean O₂ | 0.295 | **0.352** |
> | cisplatin `kill_fraction` | 0.3699 | **0.4361666386695784** (~18% higher) |
>
> **Cause:** the endothelial classifier fix (`matched_null_z` aggregation, which removed an
> arity bias in marker-gated phenotyping). Fewer cells are typed as tumour and the survivors
> sit in better-perfused tissue, so mean O₂ rises and kill rises with it.
>
> **This is not a report-layer bug.** The current value was reproduced *bit-identically* by
> two independent paths — the static-export builder and `app/server/report.py` — both landing
> on `0.4361666386695784`, and by re-running the CLI chain by hand on the current tree.
>
> This is exactly the failure mode this project's own limits describe: **the ranking is
> invariant, the level is not.** Regeneration is pending a tree in which
> `tmesim/reconstruct/real_pipeline.py` is not being actively edited; regenerating against a
> half-applied engine change would bake in a state that never existed.
>
> The same staleness applies to [`results/patient_report_xenium.md`](patient_report_xenium.md).

> **Regenerated under the discrete vessel source.** `vascular_mode = discrete`
> (`reconstruct.real_fields._DEFAULT_VASCULAR_MODE`, promoted in
> [`results/vessel_promotion.md`](vessel_promotion.md)) · **alpha = 148.6**, the
> e-folding-length convention `L = d_half / ln 2` · per-cell sensitivity mode =
> **`uniform` (headline ranking); `determinant` reported separately** · **dx = 12 µm**.
>
> Every absolute exposure / AUC / coverage / kill below is on the **discrete-source
> unit** and is NOT comparable with the KDE-era numbers this file carried before this
> regeneration (the tissue-integrated source total is 3.36x smaller; that is a change
> of *unit*, not of physics — the dose-matched, geometry-only effect on coverage is
> <= 0.0166). To recover the old numbers, re-run with `mode='kde'`.
>
> Headline shift vs the KDE-era file: mean tumour AUC 5.48e+06 -> 1.63e+06 µM·s (cisplatin) and kill 0.9125 -> 0.3699; the **kill ordering is unchanged** — cisplatin > doxorubicin > lomustine > methotrexate > temozolomide.

**Experiment:** `experiments/real_drug_panel.py` (run: `./.venv/bin/python experiments/real_drug_panel.py`)

**Figure:** `results/real_drug_panel.png`

> **Provenance of this run.** Perivascular-sink calibration **α = 148.6** (the e-folding-length convention, `L = d_half/ln 2`; the retired α = 309.6 came from fitting the raw half-distance as if it were an e-fold length and is superseded). Per-cell sensitivity mode = **`uniform`** for the headline ranking — identical for every drug, so the comparison is fair; the drug-specific `determinant` mode is reported separately in §*Sensitivity mode* below. Tissue mask = **geometric** (union of 30 µm discs around the measured cell centroids), which is **independent of the grid spacing** — this replaces the old "grid cells containing ≥ 1 cell" mask and changes the oxygen field and the hypoxic fraction substantially (see §*Tissue mask*).

The full panel **[temozolomide, lomustine, cisplatin, methotrexate, doxorubicin]** is simulated on **one and the same real tumor** — the open 10x Xenium GBM section (40,887 measured cells; tumor 30,174). The real microenvironment is reconstructed **once** and reused for every drug, so all differences come from each drug's own physics on identical tissue (real vasculature, real ECM from collagen/matrix transcripts, real-perfusion oxygen).

Each drug is compiled through the same pipeline: `calibrate_alpha()` (single perivascular-sink constant, fit **once** on doxorubicin → α=148.6) → `calibrated_drugparams(name, α)` (real GDSC IC50 + calibrated loss) → `apply_barrier_to_drugparams` (saturable binding-site barrier).

## Real microenvironment (reconstructed once)

- Source: 10x Xenium GBM — 40,887 real cells, grid (424, 675), dx=12 µm, extent (8100.0, 5088.0) µm.
- Cell-type fractions: TUMOR 73.8%, ENDOTHELIAL 9.7%, IMMUNE 8.5%, STROMAL 7.9%.
- Mean vascular density 0.025, mean ECM 0.076, mean oxygen 0.295.
- Tissue mask: **geometric**, radius 30 µm, covering 67.7% of the bounding box (dx-independent — see §*Tissue mask*).
- **Tumor hypoxic fraction ≈ 12.7%** at this resolution under the grid-independent tissue mask. See §*Tissue mask* — this number moved a great deal when the mask was fixed and the honest reading has changed.

## Tissue mask — the geometry fix, and what it changed

`Domain.tissue_mask` used to be *"grid cells that contain at least one segmented cell"*. That is a property of the **grid**, not of the tissue: as `dx` shrinks each cell can still only mark one tile, so the covered fraction collapsed — **50.6% / 32.7% / 13.9%** of the bounding box at dx = 30 / 20 / 12 µm. `real_fields._kde_on_grid` multiplies the vascular field by that mask, so the **total O₂ supply shrank as the mesh was refined**, which is a discretisation artifact, not biology.

The mask is now **geometric**: a grid cell is tissue when its centre lies within **30 µm** of at least one measured cell centroid (the union of fixed physical discs). The radius is a declared constant justified against this section's own geometry — the nearest-neighbour spacing between the 40,887 measured centroids is median **16.7 µm** (p10 10.4, p90 24.2, p99 34.2 µm), so 30 µm is ~1.8× the median and just above the p90 spacing: the union of discs bridges >90% of real inter-cell gaps (giving a *connected* domain, which the diffusion operator needs) while extending only ~1 cell diameter beyond the outermost measured cell. It is not fitted to anything.

Measured effect (same section, same physics, only the mask definition changed):

| dx (µm) | coverage — legacy | coverage — geometric | mean O₂ @cells — legacy | mean O₂ @cells — geometric | hypoxic tumor frac — legacy | hypoxic tumor frac — geometric |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 0.506 | 0.672 | 0.522 | 0.570 | 0.000 | 0.000 |
| 20 | 0.327 | 0.674 | 0.426 | 0.565 | 0.000 | 0.000 |
| 12 | 0.139 | 0.677 | 0.251 | 0.558 | 0.134 | 0.000 |
| **spread (max/min)** | **3.64×** | **1.007×** | **2.08×** | **1.02×** | — | — |

So the geometric mask does what it was built to do: tissue coverage and mean oxygen are now invariant to the grid (1.007× and 1.02× across a 2.5× range of dx), where before they varied by 3.6× and 2.1×.

**Honest consequence — the hypoxia story changes again, and not in the model's favour.** With the O₂ supply no longer artificially starved by a shrinking mask, mean oxygen at the measured cell positions rises to ≈0.56–0.57 and **no tumor cell falls below the 0.2 hypoxia threshold at any dx: the predicted hypoxic fraction is 0.0%.** For a glioblastoma section — a tumor type defined partly by pseudopalisading necrosis — that is not credible as biology. It is best read as evidence about the *oxygen model*, not about this tumor: `OxygenConfig.consumption` has never been calibrated against a measured oxygen tension, and the previously-reported hypoxia (94% at dx = 30, 13% at dx = 12) was being manufactured by the grid-dependent supply term rather than by the tissue. **Both the old numbers and this 0% are un-calibrated model output; the correct statement is that this simulator cannot currently predict this section's hypoxic fraction, and no conclusion in this report should rest on it.** What the fix does establish is that the answer no longer depends on the mesh.

## Sensitivity mode — why the headline ranking uses `uniform`

Per-cell drug sensitivity can be set two ways (`reconstruct.real_pipeline.apply_drug_sensitivity`), and the choice changes the ranking:

- **`mode='uniform'` (default, and the headline table above).** Every drug gets sensitivity 1.0 at every cell. Because the factor is *identical across drugs* it cancels out of the comparison, so the ranking isolates **delivery × potency** with no sensitivity confound. Absolute kills are an equal no-intrinsic-resistance **upper bound for every drug**.
- **`mode='determinant'`.** Each drug uses its own measured determinant where the 319-gene Xenium panel carries one — **TOP2A** for doxorubicin (its literal molecular target) and an **MKI67** proliferation proxy for cisplatin / lomustine / methotrexate — and uniform 1.0 where it does not. On this section both determinants land at mean sensitivity ≈0.585. **Temozolomide has no determinant on this panel** (its is MGMT, which is absent and is a methylation state that transcript counts cannot report), so TMZ *alone* keeps 1.0 while every competitor is scaled down ~1.7×.

**That asymmetry is an artifact of measurement coverage, not of drug biology.** Any kill temozolomide gains in `determinant` mode is bought purely by its resistance gene being unmeasured. It is reported here, flagged **UPPER BOUND**, and it is NOT the ranking:

| drug | sensitivity (uniform) | sensitivity (determinant) | determinant gene | kill (uniform) | kill (determinant) | comparable across drugs? |
|---|---:|---:|:--|---:|---:|:--|
| temozolomide | 1.000 | 1.000 | — (none on panel) | 0.0147 | 0.0147 | ❌ **UPPER BOUND** — no resistance modelled |
| lomustine | 1.000 | 0.590 | MKI67 | 0.1975 | 0.1190 | scaled by its own determinant |
| cisplatin | 1.000 | 0.590 | MKI67 | 0.3699 | 0.2197 | scaled by its own determinant |
| methotrexate | 1.000 | 0.590 | MKI67 | 0.0162 | 0.0096 | scaled by its own determinant |
| doxorubicin | 1.000 | 0.585 | TOP2A | 0.2165 | 0.1301 | scaled by its own determinant |

**Ranking under `mode='uniform'` (FAIR — this is the ranking):** **cisplatin** (0.3699) > **doxorubicin** (0.2165) > **lomustine** (0.1975) > **methotrexate** (0.0162) > **temozolomide** (0.0147)

**Ranking under `mode='determinant'` (CONFOUNDED — reported for comparison only, do not cite as a ranking):** **cisplatin** (0.2197) > **doxorubicin** (0.1301) > **lomustine** (0.1190) > **temozolomide** (0.0147) > **methotrexate** (0.0096)

Read the difference plainly: in `determinant` mode every drug except temozolomide is multiplied by ≈0.585, so their kills drop by roughly that factor while TMZ's does not move at all. Any place TMZ climbs is **entirely** attributable to that, and would reverse the moment MGMT were measurable. `determinant` mode is legitimate for a *within-drug* spatial question ("where on this section does doxorubicin work?"); it is not legitimate for ranking drugs against each other on a panel that measures some of their determinants and not others.

## Predicted panel on this real tumor

*(sensitivity mode = `uniform` — identical for every drug.)*

| drug | MW (Da) | real IC50 (µM) | bbb_penetrance | D_eff (µm²/s) | penetration_index | mean AUC (µM·s) | predicted kill |
|---|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194 | 606 | 0.125 | 265.3 | 0.470 | 6.86e+05 | 0.0147 |
| lomustine | 234 | 50 | 0.973 | 247.4 | 0.370 | 8.35e+05 | 0.1975 |
| cisplatin | 298 | 24.6 | 0.858 | 225.4 | 0.420 | 1.63e+06 | 0.3699 |
| methotrexate | 454 | 12 | 0.023 | 189.8 | 0.331 | 1.41e+04 | 0.0162 |
| doxorubicin | 544 | 0.5 | 0.021 | 176.3 | 0.310 | 1.03e+04 | 0.2165 |

IC50 values are the real GDSC/literature potencies that flow in through the molecular layer; MW is from the verified molecular structure; `bbb_penetrance` and `D_eff` are structure-derived (`D_eff = D · (1 − 0.6·ECM)` at the tumor mean ECM = 0.076).

## Which drug wins on THIS real tumor, and why

- **Penetration does NOT separate these drugs at converged resolution.** `penetration_index` spans only 0.310–0.470 across the whole panel — a spread of 0.161, which is not a meaningful ranking. This is a **correction**: an earlier, unconverged dx = 30 µm run of the same panel appeared to span 0.266–0.396 and was reported as 'temozolomide penetrates best, doxorubicin barely reaches the tumor bulk'. **That ranking was a grid artifact and is retracted.** The reason is definitional: `penetration_index = mean(AUC) / p99(AUC)` over tissue is a *scale-invariant shape* statistic — it measures how uniform the exposure field is, not how much drug arrives. Once the perivascular boundary layer is actually resolved, that shape is set by the **vascular geometry**, which is identical for all five drugs, and these drugs' D_eff values differ by only ~1.5× (175–263 µm²/s) — far too little to reshape the field.
- **Delivery differences are real, but they live in the AUC, not the index.** Mean exposure spans 1.03e+04–1.63e+06 µM·s (~157×), driven by BBB penetrance and the calibrated perivascular sink. Doxorubicin really is delivery-starved on this section — its mean AUC is ~120× below temozolomide's — but that shows up as a *lower absolute exposure*, not as a flatter spatial profile. **Quote mean AUC, not penetration_index, when comparing drugs.**
- **Best predicted kill (mode = `uniform`):** **cisplatin** (kill_fraction 0.3699, penetration_index 0.420, IC50 24.6 µM, mean AUC 1.63e+06 µM·s, bbb 0.858). Kill is **delivery × potency**, discounted by whatever oxygen-driven resistance the tissue carries: a drug wins only when it both reaches the cells AND is potent enough at the exposure it achieves. Full ranking: **cisplatin** (0.3699) > **doxorubicin** (0.2165) > **lomustine** (0.1975) > **methotrexate** (0.0162) > **temozolomide** (0.0147). Because every drug carries the same sensitivity in this mode, none of this ordering is a sensitivity artifact — it is delivery and potency only.

## Honesty and caveats

- **No outcome label.** This open Xenium section is an **untreated resection**. There is no paired measured drug response for any of these five drugs. The kill/penetration numbers are **predictions of the validated physics on real tissue — not validated outcomes.** Do not read the ranking as clinical efficacy.
- **What is genuinely real here:** the geometry (real cells, real vasculature, real matrix, real oxygen), the transport physics (structure → mesh diffusivity + a doxorubicin-calibrated perivascular sink, revalidated out-of-sample on mitoxantrone), and the potency (real GDSC/literature IC50).
- **Hypoxia is now predicted at 12.7%, and that is a statement about the model, not about this tumor.** Every earlier hypoxia number from this pipeline (≈94% at dx = 30 µm, ≈13% at dx = 12 µm) was produced under the grid-dependent occupancy tissue mask, which starved the O₂ supply more and more as the mesh was refined; those numbers are **retracted**. Under the fixed geometric mask the hypoxic fraction is 0.0% at every dx — which is *also* not credible for a GBM section. `OxygenConfig.consumption` has never been calibrated against a measured oxygen tension, so the absolute level is un-anchored in both directions. **No claim in this report depends on the hypoxic fraction, and none should.**
- **Absolute kill fractions are an upper bound.** In the headline `uniform` sensitivity mode no intrinsic per-cell resistance is modelled for any drug. That is what makes the *ranking* fair; it also means every *absolute* kill here is the best case. See §*Sensitivity mode*.
- **One anchor sets the scale.** The single α is fit to one doxorubicin penetration number; the absolute penetration of every drug rides on that measurement, and the binding proxy is an ordering, not a first-principles rate.

## Resolution & convergence

This panel was run at **dx = 12 µm** on grid (424, 675). That choice is not cosmetic — see [`results/numerical_rigor.md`](numerical_rigor.md) §2 (grid refinement on this same real section):

- **`kill_fraction` is grid-converged; `penetration_index` requires dx ≲ 12 µm.** The calibrated perivascular boundary layer is ~45 µm thick, so any coarser grid cannot resolve the steep near-vessel gradient and reports an artificially flat — i.e. too high — index. That is why this panel runs at dx = 12 µm. dx = 12 µm is also the finest grid that satisfies the upwind-advection CFL condition at the published dt = 10 s.
- **An earlier dx = 30 µm run of this panel OVER-ESTIMATED `penetration_index`** by 1.4–2.0×, drug-dependently, which manufactured a spurious penetration ranking (see above). Those figures are superseded; do not cite the dx = 30 µm penetration numbers.
- **⚠️ The per-drug convergence ladder previously printed here is RETRACTED, not merely updated.** It was measured under the legacy grid-dependent occupancy tissue mask, so its dx sweep confounded genuine mesh refinement with a shrinking simulated tissue domain — exactly the defect §*Tissue mask* fixes. Re-deriving a 5-drug × 6-grid ladder under the geometric mask is a separate compute job that was not run for this report. **For the current, geometric-mask grid study see [`results/numerical_rigor.md`](numerical_rigor.md) §2**, which is regenerated alongside this file on the same geometry (doxorubicin reference drug, advection off).
