# Does richer endothelial marker gating fix breast's near-CSR vasculature? (lever (a) from `results/vascular_validation_breast.md` §6)

Standalone prototype, not committed to `tmesim/`: calls
`tmesim.reconstruct.real_xenium.phenotype_cells()` (unmodified, its documented
`marker_map` override) and `tmesim.validation.vascular_validation`'s free
functions (`cluster_vessels`, `distance_to_nearest`, `tissue_footprint`,
`poisson_distance_stats`, `local_heterogeneity`) verbatim and unmodified, on
`data_cache/xenium_breast_janesick/rep1` at `dx_um=12` (the exact section and
grid `results/breast_cancer_flagship.md` and `results/vascular_validation_breast.md`
use). Script: `endothelial_marker_richness_prototype.py` (scratchpad, not under
`tmesim/`). Nothing under `tmesim/` was edited to produce this report,
including `tmesim/reconstruct/marker_panels.py`, which is out of scope per this
session's "prototype first, zero core-engine changes" discipline.

## 1. The panel gene list, read directly from the data, not from memory

`data_cache/xenium_breast_janesick/rep1/cell_feature_matrix/features.tsv.gz`
(cross-checked against the raw panel manifest,
`Xenium_FFPE_Human_Breast_Cancer_Rep1_gene_panel.json`, design `hBreast_320g`)
carries **313** `Gene Expression` targets. `marker_panels.BREAST_PANEL_MARKER_MAP`
types `ENDOTHELIAL` from 3 of them (`PECAM1`, `VWF`, `CLDN5`).

**Correction to prior documentation**: `marker_panels.py`'s own docstring and
both `results/vascular_validation_breast.md` and
`results/breast_cancer_flagship.md` say the panel carries "13
endothelial-associated genes" and name 10 more beyond the 3 used (`CLEC14A,
EGFL7, KDR, RAMP2, MMRN2, SOX17, SOX18, NOSTRIN, AQP1, ESM1, ANGPT2`). Counting
that named list directly: it has **11** entries, not 10, making **14** total,
not 13. Grepping the panel additionally turns up a **12th unused candidate**
those reports never named: **`CD93`** (pan-endothelial + monocyte/platelet
marker, `AA4`/`C1QR1`), present at 18.04% detection — higher than 9 of the 11
previously-named candidates. So the panel actually carries **15** plausible
endothelial-associated genes (3 used, 12 unused), not 13 (3 used, 10 unused).
This is a minor bookkeeping correction, not a substantive finding, but the task
asked for the list read from data rather than trusted from memory, and it
was off by one gene (plus one gene missing entirely) in every prior write-up
that quoted it.

Detection rate of every candidate (matched-null-z's own abundance statistic,
`real_xenium.detection_rates`, fraction of the 167,780 cells with ≥1 transcript):

| gene | detect % | used in production? |
|---|---:|---|
| PECAM1 | 25.49 | yes |
| VWF | 16.22 | yes |
| CLDN5 | 1.50 | yes |
| CD93 | 18.04 | no |
| AQP1 | 22.48 | no |
| KDR | 8.54 | no |
| RAMP2 | 8.67 | no |
| CLEC14A | 8.50 | no |
| NOSTRIN | 6.92 | no |
| MMRN2 | 6.95 | no |
| ANGPT2 | 5.18 | no |
| EGFL7 | 3.85 | no |
| SOX17 | 3.56 | no |
| ESM1 | 2.86 | no |
| SOX18 | 2.24 | no |

## 2. Which of the 12 unused genes are genuinely well-established, low-noise endothelial markers

Judged against standard vascular-biology / pan-endothelial-atlas literature
(Kalucka et al. 2020 Cell PMID 32059776 murine EC atlas panel; individual
discovery/characterisation papers), not thrown in uncritically:

**Included in the "selective" composite (8 additions, well-justified):**
- `CLEC14A` — endothelial-restricted C-type lectin, an established tumour
  endothelial marker (non-CD31 pan-EC surface marker, used in vascular biology
  FACS/IHC panels).
- `EGFL7` — endothelial-specific secreted matrix protein, vascular-development
  literature since its 2004 characterisation; near-exclusively EC-expressed.
- `KDR` (VEGFR2/FLK1) — the canonical VEGF receptor, gold-standard endothelial
  identity gene, including in tumour vasculature.
- `RAMP2` — adrenomedullin receptor accessory protein, strongly endothelium-
  enriched (Human Protein Atlas: "enhanced in endothelial cells").
- `MMRN2` (multimerin-2) — endothelial-specific ECM glycoprotein, discovered
  and used specifically as a tumour endothelial marker (Lorenzon et al. 2012).
- `NOSTRIN` — eNOS-trafficking protein, discovered as, and largely restricted
  to, vascular endothelium.
- `ESM1` (endocan) — endothelial-cell-specific secreted proteoglycan, an
  established tumour-angiogenesis marker.
- `ANGPT2` — stored with VWF in endothelial Weibel-Palade bodies, a
  well-established (though more activation/state-dependent than the others —
  see the sensitivity variant below) marker of tumour vasculature.

**Deliberately excluded — genuine noise risk, not just conservatism:**
- `AQP1` — broadly expressed outside endothelium (erythrocytes, renal
  proximal tubule, choroid plexus); its 22.48% detection rate (2nd-highest of
  all 12 candidates, well above every other exclusively-endothelial gene here)
  is itself circumstantial evidence of non-endothelial signal in this tissue,
  most plausibly red-cell contamination in vascular-adjacent regions.
- `CD93` — expressed on endothelium AND myeloid cells (monocytes, some
  granulocytes) AND platelets; adding it risks pulling `IMMUNE`-lineage cells
  (which this panel already types via `PTPRC/CD68/AIF1/CD163`) into
  `ENDOTHELIAL`, the opposite of a cleaner signal.
- `SOX17` — a real endothelial (particularly arterial) transcription factor,
  but also a primitive-endoderm/gut-epithelium lineage gene with documented
  aberrant re-expression in some carcinomas; breast ductal/glandular
  epithelium carries some risk of cross-reactivity a purely-vascular gene
  would not.
- `SOX18` — endothelial/lymphatic developmental TF, but also expressed in
  hair-follicle and skin appendage lineages; lower risk than `SOX17` in this
  tissue but still a TF with documented non-vascular expression, not a
  structural/surface identity gene.

A second, narrower **10-gene** variant drops `ANGPT2` from the selective set
(it is the most state/activation-dependent of the eight, rather than a stable
identity marker) as a sensitivity check on whether it is load-bearing. A
deliberately uncritical **15-gene "kitchen sink"** variant (selective 11 +
`SOX17`/`SOX18`/`AQP1`/`CD93`) is included specifically to test whether
throwing in everything the panel offers, noise and all, does something
different from the more careful subset — i.e., to make the "don't just throw
all 10 (12) in" caution an empirical comparison rather than an assertion.

## 3. Result: MVD, endothelial count, and heterogeneity all move; the CSR ratio does not improve, and on the module's own headline diagnostic it gets worse

All four variants use the identical `phenotype_cells(..., aggregation=
'matched_null_z', threshold=1.0)` production scoring rule and the identical
`TUMOR`/`IMMUNE`/`STROMAL` marker lists — only the `ENDOTHELIAL` gene list
changes, isolating lever (a) from every other typing choice.

| variant | ENDOTHELIAL genes (n) | n_endo | endo % of cells | endo density /mm² |
|---|---|---:|---:|---:|
| baseline (production) | 3 | 13,677 | 8.15% | 354.6 |
| selective | 11 | 16,693 | 9.95% | 432.8 |
| selective, no ANGPT2 | 10 | 16,298 | 9.71% | 422.6 |
| kitchen sink | 15 | 16,770 | 10.00% | 434.8 |

Richer typing is **not** a no-op: it adds ~22% more endothelial cells
(13,677 → ~16,700) and raises the raw endothelial density from 355 to
~423-435/mm², confirming the panel-richness lever does real, measurable work
on cell identity.

**The self-referential CSR test, at the module's own 20 µm-linking
convention** (median tumour→vessel distance vs. a Poisson null at the
reconstruction's own vessel density — the "single most important row" both the
GBM and breast validation reports call this):

| variant | MVD /mm² (20µm) | csr_median_ratio | csr_tail_ratio |
|---|---:|---:|---:|
| baseline (production) | 149.2 | 0.907× | **4.801×** |
| selective (11 genes) | 158.3 | 0.798× | **2.468×** |
| selective, no ANGPT2 (10 genes) | 159.6 | 0.795× | **2.215×** |
| kitchen sink (15 genes) | 151.3 | 0.814× | **2.383×** |

The tail ratio — the number that says "this reconstruction has a heavier far
tail than a random vessel field would, i.e. it carries real vascular
heterogeneity" — is **cut roughly in half by richer typing, in every variant,
including the careful 8/10-gene selective one**: 4.80× → 2.2-2.5×. That moves
**toward** 1.0 (toward indistinguishable-from-random), not away from it. The
median ratio also moves further below 1.0 (0.91× → ~0.80×), i.e. further from
random but in the "shorter than random" direction that, per
`vascular_validation_breast.md` §3.4, was already the harder-to-interpret of
the two signs (real vascular clustering should mean a *heavier tail*, not
necessarily a shorter median).

**This instability is worse, not better, across the full 15/20/25/30 µm
linking sweep** (the same robustness check `vascular_validation_breast.md`
§3.3 ran):

| link (µm) | baseline tail ratio | selective (11g) tail ratio | selective (10g, no ANGPT2) tail ratio | kitchen sink tail ratio |
|---:|---:|---:|---:|---:|
| 15 | 21.53× | 16.70× | 14.56× | 15.58× |
| 20 | 4.80× | 2.47× | 2.22× | 2.38× |
| 25 | 2.12× | **0.95×** | **0.83×** | **0.97×** |
| 30 | 1.20× | **0.50×** | **0.43×** | **0.55×** |

At the two wider linking distances (25/30 µm) every richer-typed variant's
tail ratio drops **below 1.0** — meaning the reconstruction now shows a
*thinner* far tail than a matched random field, the opposite of the
real-vasculature clustering signature this whole exercise is looking for. The
production 3-gene baseline never does this (it stays ≥1.2× everywhere in the
sweep). Richer typing makes the module's headline diagnostic *more* sensitive
to an arbitrary free parameter (linking distance) in a way that more often
points the wrong direction, not less.

**Decomposition A (raw endothelial-cell density, no clustering correction)**
tells a matching story: `csr_median_ratio` drops from 1.252× (baseline) to
1.11-1.16× (all three richer variants) — again moving toward 1.0, i.e. toward
CSR, not away from it.

**The one diagnostic that does move in the direction hoped for is the coarse
500 µm-tile heterogeneity ratio** (endothelial-cell density p90/p10 across
tiles holding ≥100 cells):

| variant | p90/p10 |
|---|---:|
| baseline (production) | 6.56× |
| selective (11 genes) | **9.09×** |
| selective, no ANGPT2 | 8.92× |
| kitchen sink | 8.58× |

This rises by ~40%, moving further from CSR-uniform and closer to (though
still well short of) Haisan et al.'s 4.6× within-tumour compartment spread
cited in `vascular_validation_breast.md` as the only breast-specific
heterogeneity anchor available (that anchor is *lower* than even the
baseline's 6.56×, so "higher" here is moving away from, not toward, the one
number that exists to check against — read this row cautiously, not as a
clean win).

## 4. Verdict

**Richer marker gating is a real, measurable lever — it changes the
endothelial cell count, the endothelial density, and the MVD by double-digit
percentages — but it does not fix, and on the module's own headline
vessel-clustering diagnostic (the CSR tail ratio at 20 µm linking, and
increasingly at 25-30 µm) it makes the reconstruction MORE indistinguishable
from — or even more uniform than — a random Poisson vessel field, not less.**
Only one secondary diagnostic (coarse tile-level heterogeneity) moves in the
hoped-for direction, and even that moves away from the one external anchor
available for it. This holds for the careful, literature-justified 8/10-gene
addition exactly as much as for the deliberately uncritical 15-gene kitchen
sink — the three richer variants are far more similar to each other (tail
ratio 2.2-2.5× at 20 µm) than any of them is to the 3-gene baseline (4.8×), so
gene *selectivity* was not the deciding factor here; gene *count/coverage*
was, and it moved the wrong way.

This is an honest null in the same family as the bandwidth sweep and the
adaptive-kNN field result already on record (`vascular_validation.md` §6a):
another concrete, previously-untried lever for the "vasculature looks too
close to random" finding was tried in full and did not resolve it — if
anything it modestly sharpened the negative result on the module's primary
diagnostic. Combined with those results, this strengthens rather than weakens
the case that the near-CSR finding reflects a real property of this section's
measured cell positions (or a limit of DBSCAN single-linkage clustering as the
cells-to-vessels object model — lever (b), being tested separately), not an
artefact of which 3 of the panel's endothelial genes happen to be used for
typing.

## 5. Caveats

- **One section, one patient**, exactly as `vascular_validation_breast.md`
  already flags — this does not test whether the direction of this result
  would replicate on Rep2 or a different breast section.
- **`matched_null_z` aggregation and `threshold=1.0` were held fixed** across
  every variant (the production defaults), isolating the marker-list choice
  from the scoring-rule choice; a different aggregation rule was not swept
  here.
- **`TUMOR`/`IMMUNE`/`STROMAL` cell counts move slightly as a side effect**
  (`n_tumor` falls from 89,785 to ~87,500-87,800 as more cells clear the
  richer `ENDOTHELIAL` gate instead of defaulting to `TUMOR`) — an honest,
  expected consequence of `phenotype_cells`' shared arg-max across types, not
  a separate finding, but it means `tumor_to_vessel` distances in this report
  are computed against a very slightly different tumour-cell set per variant,
  matching how `measure()` itself would behave under a swapped marker map.
- **No vascular-density field (KDE/discrete) was rebuilt or validated here** —
  this report tests only the CSR/vessel-clustering diagnostic on cell
  positions, exactly the piece `vascular_validation_breast.md` §6a identifies
  as untouchable by field-smoothing choices; it says nothing new about the
  field's own values.
- **The CD93 finding (a 12th unused candidate the prior reports missed) is a
  minor bookkeeping correction**, not independently validated as "should be
  excluded" beyond the standard literature reasoning given in §2 — it was
  never in the "selective" or "no ANGPT2" variants and only entered the
  kitchen-sink variant, where it did not change the overall conclusion.
- **Not run through the repo test suite or committed** — a standalone
  prototype script per this session's discipline; `tmesim/reconstruct/
  marker_panels.py` is unmodified and no change to it is recommended by this
  report (see §4).

---
*Reproduce:* `endothelial_marker_richness_prototype.py` (scratchpad, not
committed). Requires `data_cache/xenium_breast_janesick/rep1` (already
cached). No file under `tmesim/` is imported-and-monkeypatched or edited;
every number above comes from calling `phenotype_cells()` and
`tmesim.validation.vascular_validation`'s functions with a custom
`marker_map` argument, the extension point those functions already expose.
