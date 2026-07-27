# Rigor Audit — what the verification pass found, and what still stands

This session put the real-data build through an adversarial verification pass. It found
**nine distinct defects**, several of which were mutually cancelling, and produced the
project's **first validation against measured drug concentrations in human tissue**.

Every headline number reported before this audit was wrong in at least one way. The
engine and the reconstruction were sound; the *claims* were not.

---

## 1. The first real external validation — human glioma drug exposure

Scored against **125 human glioma tissue samples across 42 patients** plus 27
PMID-anchored regional PK measurements (`results/glioma_regional_pk.md`).

| Axis | Result |
|---|---|
| **Partition magnitude** | 5/6 drugs within 10x, 4/6 within 3x. Aggregate bias 2.02x — **but 1.01x (essentially unbiased) excluding one outlier** |
| **Cross-drug ordering** | Measured penetrant-vs-excluded separation 4.66x; **model predicts 3.65x from structure alone** |
| **Regional ordering** | **FAIL.** Measured CE-vs-non-enhancing 1.45x-4.26x (p<0.0001, within-patient permutation); model predicts exactly 1.00x |

The regional failure has a precise cause — there is no barrier-integrity term anywhere in
the model — and the validation **derived the fix from the data**: a multiplicative barrier
factor cannot work (3.2x spread across drugs), but an **additive leak is drug-independent
to within 1.3x** (L = 0.211 / 0.230 / 0.272 across three independent studies), giving
`source = vascular_density * (bbb_penetrance + L * bbb_disruption) * plasma`, L ~ 0.24.
That is the signature of size-nonselective paracellular opening. **Not yet implemented.**

It also traced a 66x outlier (eribulin) to a genuine code defect: `bbb_penetrance` takes
the **mean** of four desirability scores, so a molecule with two disqualifying properties
is rescued by two favourable ones. *A mean cannot express a veto.* **Not yet fixed.**

---

## 2. Defects found and fixed

| # | Defect | Status |
|---|---|---|
| 1 | Endothelial fraction 2x too high (single-marker argmax bias): 20.6% -> 9.7% | fixed |
| 2 | `penetration_index` not grid-converged **and definitionally unable to discriminate drugs** | retired, replaced |
| 3 | Hypoxic fraction a grid artifact (94% at dx=30 -> 13% at dx=12 -> 7% at dx=8, not converged) | retracted |
| 4 | alpha off by exactly 1/ln^2(2) — half-distance vs e-fold convention (309.6 -> **148.6**) | fixed + guarded by test |
| 5 | Per-cell sensitivity a placeholder suppressing all kills ~5x | fixed (drug-specific, real markers) |
| 6 | Grid-dependent oxygen sink (raw count, scales as dx^2) | fixed |
| 7 | Grid-dependent tissue mask — **was partially cancelling #6** | diagnosed, prototype validated, **not applied** |
| 8 | Unfair cross-drug comparison: TMZ advantaged because its resistance determinant is *unmeasured* | diagnosed, **not fixed** |
| 9 | `bbb_penetrance` mean cannot express a veto | diagnosed, **not fixed** |

**#7 is the most instructive.** Two grid-dependence bugs pulled in opposite directions, so
the output looked *more* stable than it was. Compensating errors that produce plausible
results are invisible to a passing test suite — 55 tests were green throughout. Only
writing the math down (`docs/MATHEMATICAL_MODEL.md`) and checking invariants found them.

---

## 3. The structural confound: our vasculature is a blur, not vessels

Three independent analyses converged on one root cause. The vascular source field is a
**smooth Gaussian KDE of endothelial cell positions** (bandwidth ~40 um), so:

- fitted penetration length is **26-28 um for ALL five drugs** against a **source-null of
  23.7 um** — the decay measures KDE bandwidth, not drug transport (`delivery_metrics.md`);
- the exposure profile has a **flat floor at ~12% of peak**, so no tissue is truly
  drug-starved and *delivery is never the binding constraint* (`dose_schedule.md`);
- "perivascular regions get more drug" is **circular** — the region label and the drug
  source share the endothelial signal (`niche_response.md`).

Real tumour vasculature is **discrete**: distinct vessels with genuinely avascular tissue
between them. Replacing the KDE with a segmented vessel network is the single
highest-value structural fix outstanding. Literature basis is established:
**~20 endothelial cells per mm of vessel, 1-2 EC per capillary cross-section**
(Shapson-Coe 2024 *Science* volume EM; Sargent 2023) — i.e. one EC ~ 50 um of vessel.

---

## 4. What survives the audit

- **Reconstruction fidelity**: cross-platform Pearson 0.997; endothelial agrees to 0.03%
  between two independent platforms. And an external check vindicates it — our **9.7%**
  endothelial matches human cortex volume-EM tissue counting (**8.0%**, Shapson-Coe 2024),
  *not* the 1-2% from dissociation-based scRNA-seq. Correcting toward scRNA-seq would have
  introduced a 5x error into the delivery source term.
- **Numerics**: mass conserved to ~1e-13 on the real irregular mask; dt-converged; kill
  fraction grid-converged beyond dx<=20 um. Production resolution is **dx = 12 um** (also
  the finest grid satisfying the advection CFL at dt = 10 s).
- **Physics vs real PK**: one constant fit to doxorubicin predicts mitoxantrone (48 um) and
  chemically-distinct topotecan (84 um) inside their real cited bands, out-of-sample.
- **Potency**: GDSC vs HGCC agree at Spearman +0.83 (p=2.4e-5), real-vs-real.
- **The delivery x potency result**: ranking by **AUC/IC50** reproduces simulated kill
  ordering exactly. Delivery spans 120x, potency 1212x; the product decides. Doxorubicin
  is 49x more potent than cisplatin, 120x less delivered, and loses.
- **Recommended headline metric**: `therapeutic_coverage @ 1x IC50` — the fraction of tumour
  cells whose exposure exceeds the drug's measured IC50. Dimensionally sound,
  grid-converged (+-0.03), discriminating (0.311-0.925), and it recovers the kill ordering
  without touching the Hill curve, sensitivities, hypoxia factors or the Bernoulli draw.
  It is computed over the *fixed set of real cells*, which is why it does not inherit the
  KDE geometry confound.

---

## 5. The data landscape, settled

A 39-agent sweep (164 datasets, 143 open) established (`results/DATA_HUNT.md`):

- **No open dataset anywhere pairs measured spatial drug distribution with measured
  treatment response in the same tumour.** Not one.
- **No glioma + temozolomide spatial drug map exists**, open or otherwise (METASPACE: zero).
- **MGMT is absent from every Xenium panel** — so temozolomide sensitivity is not merely
  unmodelled but *unmodellable* from this modality.
- The experiment that would validate this project end-to-end already exists — the
  Agar/White adavosertib GBM platform, drug concentration per pixel in phase-1 patient
  tissue with molecular response — and has **no accession**. It is a months-long
  collaboration request, not a download.

What *is* open and downloaded: GIST imatinib MSI (6 matched human tumour/normal pairs, with
a quantitative calibration series), Moncayo HGSOC explants (114 ROIs pairing measured drug
with 11,224 genes), the human glioma regional PK compendium, and GSE232173 (glioma + TMZ +
measured ex-vivo response).

---

## 6. Outstanding work (interrupted by the spend limit)

Highest value first:
1. Implement the BBB-disruption additive leak (L~0.24) and the veto-shaped `bbb_penetrance`.
2. Apply the geometric tissue mask (fix #7) and re-run the affected results.
3. Replace the KDE vascular field with a discrete segmented vessel network; re-test whether
   the flat floor disappears and drugs finally separate spatially.
4. Finish the three interrupted external validations: GIST imatinib, Moncayo explants,
   GSE232173 TMZ ex-vivo.
5. Finish the two interrupted internal validations: cell-state vs measured MKI67, and
   simulated oxygen vs measured VEGFA/HIF1A. **Partial result already observed: VEGFA is the
   single most anti-correlated gene of all 480 with simulated oxygen — encouraging — but
   HIF1A is exactly null.** Both need completing before either is claimed.
6. Cross-drug comparability modes (uniform vs determinant sensitivity).
7. Ranking robustness: cisplatin (0.241) and doxorubicin (0.231) differ by ~4%, well inside
   plausible parameter uncertainty. **The top-1 call is probably not decision-grade** and
   should be reported as a tie until shown otherwise.

---

## 7. The honest position

The physics engine is numerically sound and the reconstruction is externally validated. The
delivery/partition axis is now validated in human glioma with a quantified error budget. The
regional axis has a specified, data-derived fix that is not yet implemented. The spatial axis
is confounded by our own vascular smoothing. The response axis cannot be validated with open
data at all.

That is a real model with a real error budget — a fundamentally different object from a
simulator that produces confident numbers. The defects above were found because the model was
attacked, not because it failed; a version of this project that skipped the audit would still
be reporting nine wrong numbers with beautiful figures attached.
