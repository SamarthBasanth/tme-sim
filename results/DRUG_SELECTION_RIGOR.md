# Drug-Selection Rigor Study: Does Spatial-Transport Physics Help?

**Question.** Across 15 independently-authored ground-truth (GT) kill models — 12 biological
motifs plus 3 adversarial probes — does the mechanistic PDE sim select the right per-patient
drug better than (a) a potency-only shortcut and (b) a cross-validated feature-ML baseline?
Each GT was authored without reading `K_sim`, scores 30 patients × 5 drugs, and reports
top-1 accuracy and mean regret (higher acc / lower regret is better; random ≈ 0.20).

**Bottom line up front:** the sim reliably beats the naive potency shortcut (10/15), but it
almost never beats a cheap ML classifier on observable features (1/15, and that one win is
fragile and circular). The sim's value is confined to delivery-limited regimes, and even there
a logistic model on TME features usually recovers the same signal. Whether real GBM lives in
the regime where the sim adds value is **not answerable synthetically** and needs clinical data.

---

## 1. Headline counts

| Comparison | Sim wins | Mean acc gap | Median acc gap |
|---|---|---|---|
| Sim vs **potency** shortcut | **10 / 15** | **+0.207** | **+0.300** |
| Sim vs **feature-ML** baseline | **1 / 15** | **−0.429** | **−0.433** |

- **Beats potency (10/15):** multiplicative, liebig_min, threshold, hypoxia_resist,
  ecm_sieving, bbb_gated, subtype_pref, saturating, noisy_realistic, delivery×sens-gate.
- **Beats feature-ML (1/15):** delivery×sens-gate only — and its author flags that win as a
  small-sample seed artifact (see §2).
- The five motifs where the sim does **not** beat potency are the two potency-dominant nulls,
  anti_delivery, delivery_dominant (ties potency at acc 0), and subtype_blindspot (ties at 0.333).

The asymmetry is the headline: **transport physics clears the dumbest baseline but is dominated
by a cheap supervised model that sees the same observable TME features.**

---

## 2. The circularity test

We ask two things. (i) Does the sim's advantage track *how delivery-limited* the GT is
(`frac_best_not_most_potent`)? (ii) Does the sim only win when the GT happens to *agree with the
sim's own physics* (`sim_vs_gt_rank_corr`)? A win that only appears at high rank-corr is
**circular** — the GT was built from the same biophysical ingredients the sim integrates.

**Pearson correlations (n=15):**

| | vs `frac_best_not_most_potent` | vs `sim_vs_gt_rank_corr` |
|---|---|---|
| Sim − potency accuracy gap | **+0.894** | **+0.210** |
| Sim − ML accuracy gap | +0.685 | +0.493 |

**Reading it honestly:**

- **The sim's win over potency is genuine, not circular.** The advantage is strongly explained
  by delivery-limitedness (r = **0.894**) and only weakly by rank agreement with the sim's own
  physics (r = **0.21**). Concretely, the sim beats potency across the *full spread* of rank-corr,
  including plainly moderate values — threshold (corr 0.42), hypoxia_resist (0.42), liebig_min
  (0.48), ecm_sieving (0.56), noisy_realistic (0.50). It does not need the GT to echo its ranking;
  it needs the GT to be a regime where the most-potent drug is undeliverable. In the 11
  delivery-limited motifs (`frac ≥ 0.9`) the sim beats potency in **10/11**, mean acc gap +0.494.
  This is the one solid, defensible claim in the study.

- **The sim's single win over ML *is* the circular one.** The only ML-beating motif,
  delivery×sens-gate, has the **highest rank-corr in the whole set (0.882)** — its GT was
  admittedly built from the same delivery + Hill + apoptosis-dependence ingredients the sim
  integrates. Its own author ran a 20-seed robustness probe: sim beats both baselines in only
  9/20 seeds (mean acc 0.753 vs ML 0.747 — a statistical tie), and traced the seed-0 "win" to
  the sim being a *degenerate constant classifier* (K_sim.argmax = cisplatin for all 30 patients)
  that happens to match the cohort plurality. So the lone ML win is exactly the case the
  circularity test is designed to flag, and it does not survive reseeding.

**Verdict:** delivery-limitedness — not rank agreement — drives the sim's real advantage, so the
potency-shortcut win is not an artifact of building the GT in the sim's image. The ML win is the
opposite: high rank-corr, contrived, and fragile. The circularity test passes for the potency
comparison and fails (correctly caught) for the ML comparison.

---

## 3. Where the sim loses, and why — its blind spots

**A. Potency-dominant regimes (null cases): sim is the *worst* policy.**
`potency_dominant` and `potency_dominant2` both collapse to "always pick doxorubicin" (the
globally most-potent drug is best for all 30 patients; `frac_best_not_most_potent = 0`). The
potency shortcut is oracle-perfect (acc 1.0); the sim scores **acc 0.0**, because its
transport weighting systematically *down*-ranks doxorubicin (worst BBB, largest MW) — the one
drug that wins. This is honest and expected: when the intrinsic potency gap between the top two
drugs exceeds the delivery multiplier's dynamic range, delivery physics cannot reorder the top
and is strictly dominated by a one-line rule. **Blind spot: the sim over-penalizes potent-but-
poorly-delivered agents when delivery barely matters.**

**B. Pharmacodynamic / target-dependence signals the sim structurally can't see.**
`subtype_pref`, `subtype_blindspot`, hypoxia_resist, and the subtype term in `noisy_realistic`
put the deciding signal in an **observable biomarker** (GBM transcriptional subtype, hypoxia
gate, drug-specific apoptosis dependence) that a delivery/PK model has no representation for.
Feature-ML reads the one-hot subtype / hypoxic / sens features directly and wins big:
subtype_blindspot ML acc 1.0 vs sim 0.333; subtype_pref ML 0.90 vs sim 0.43; hypoxia_resist
ML 0.63 vs sim 0.33. **Blind spot: anything driven by pharmacodynamics or molecular target
dependence rather than pharmacokinetics/transport.**

**C. Delivery *is* the signal but it's legible from features, so physics adds nothing.**
The most sobering group: bbb_gated, ecm_sieving, saturating, delivery_dominant, threshold.
These are genuinely delivery-limited (so the sim beats potency) yet the delivery limit is a
near-deterministic function of one or two observable features (`vasc`, `ecm`). A logistic
classifier learns it as well or better: bbb_gated is an *exact tie* (0.867 = 0.867, ML collapses
the BBB gate to a function of `vasc`); ecm_sieving ML 0.73 vs sim 0.30; saturating ML 0.80 vs
sim 0.33. `delivery_dominant` is the sharpest embarrassment — a regime *labeled* "sim should
help," where the cached sim scores **acc 0.0** and **rank-corr −0.11** (anti-correlated) because
its residual potency weighting steers it to cisplatin while the delivery-dominant GT makes
lomustine best for 25/30. **Blind spot: even in delivery-limited worlds, if the limit is a simple
function of measurable TME features, cheap ML captures it and the PDE physics is redundant.**

**D. Adversarial delivery-paradox: physics is actively harmful.**
`anti_delivery` inverts the rule (least-deliverable drug wins via efflux saturation / size
trapping). The sim posts its worst result anywhere: **acc 0.0, regret 0.66 — worse than random
(0.45)**, rank-corr 0.024. It confidently picks exactly the wrong drugs. ML still wins (0.97)
because the paradox is a deterministic function of the observable `ecm` feature. **Blind spot:
the sim assumes more-delivered = better; when biology inverts that, delivery-awareness is a
liability.**

**Recurring degeneracy caveat (applies to nearly all motifs).** On this 30-patient, 5-drug
cohort the realized best-drug distribution is highly concentrated — most motifs have only 2 of 5
drugs ever winning (typically lomustine vs cisplatin), because ECM spans only ~0.57–0.63 and the
BBB panel is dominated by a few penetrant agents. That narrowness (i) makes ML's majority-class
baseline unusually strong, and (ii) means several "wins" and "ties" reflect base rates more than
per-patient discrimination. The cached sim is itself frequently degenerate (constant
K_sim.argmax = cisplatin), which caps its accuracy at the cisplatin base rate and carries zero
per-patient information in those motifs.

---

## 4. Honest bottom line

**When does spatial-transport modeling add value to drug selection?** On this evidence, in a
narrow band: the true best drug must be **delivery-limited** (the most-potent agent is not
deliverable — `frac_best_not_most_potent` high) **and** the delivery limit must **not be a simple
function of observable TME features** (`vasc`, `ecm`, BBB constants). Only then can transport
physics contribute signal a cheap classifier can't recover. In this 15-GT battery, essentially
zero motifs sit cleanly in that band: every delivery-limited GT here had its delivery limit
legible from one or two features, so ML matched or beat the sim in 14/15. The sim's one robust,
non-circular achievement is beating the *potency shortcut* (10/15, mean gap +0.21), i.e.
demonstrating it is "not fooled by raw potency in delivery-limited tumors." That is real but
modest — it clears a floor, it does not clear the cheap-ML bar.

**Is this answerable synthetically at all?** No — not decisively. The battery is self-defeating
by construction: to test whether physics helps, we hand-author a GT, but any GT whose delivery
term we can write down is a GT whose delivery term feature-ML can learn from the same inputs. The
only motif where the sim beat ML was the one built from the sim's own ingredients (rank-corr
0.882) — precisely the circular case. Synthetic GTs can establish *necessary* conditions
(delivery must dominate) and can falsify overclaims (physics is harmful in potency-dominant and
anti-delivery regimes), but they cannot establish that physics adds *unique* value, because a
synthetic ground truth exposes its own generative features to the ML baseline.

**What real data would settle it.** A cohort with (i) measured per-patient/per-region drug
*exposure* or a validated PK surrogate (contrast-enhancement / DCE-MRI permeability, ECM/stromal
density, hypoxia imaging), (ii) *observed* per-patient response to specific agents, and (iii)
enough patients that a feature-ML baseline has a fair shot. The decisive test is whether the
mechanistic sim beats a feature-ML model trained on the *same* clinically measurable inputs, on
held-out patients — because in the real world the true kill function is *not* a known function of
the recorded features, which is the exact condition this synthetic study cannot create. Until
then, the defensible claim is only: **the transport sim helps rank-order drugs in delivery-
limited tumors and correctly refuses to be seduced by potency — but whether real GBM drug
selection lives in the regime where that beats a cheap feature model is an open empirical
question that requires clinical PK/response data to answer.**

---

### Appendix — per-motif scoreboard

| Motif | acc_sim | acc_pot | acc_ml | frac_delivlim | rank_corr | beats_pot | beats_ml |
|---|---|---|---|---|---|---|---|
| multiplicative* | 1.00 | 0.00 | 1.00 | 1.00 | 0.40 | ✓ | ✗ (tie, ML degenerate) |
| liebig_min | 0.57 | 0.00 | 0.77 | 1.00 | 0.48 | ✓ | ✗ |
| threshold | 0.77 | 0.00 | 0.80 | 1.00 | 0.42 | ✓ | ✗ |
| hypoxia_resist | 0.33 | 0.03 | 0.63 | 0.97 | 0.42 | ✓ | ✗ |
| ecm_sieving | 0.30 | 0.00 | 0.73 | 1.00 | 0.56 | ✓ | ✗ |
| bbb_gated | 0.87 | 0.03 | 0.87 | 0.97 | 0.79 | ✓ | ✗ (exact tie) |
| subtype_pref | 0.43 | 0.00 | 0.90 | 1.00 | 0.38 | ✓ | ✗ |
| saturating | 0.33 | 0.00 | 0.80 | 1.00 | 0.77 | ✓ | ✗ |
| potency_dominant | 0.00 | 1.00 | 1.00 | 0.00 | 0.51 | ✗ | ✗ |
| delivery_dominant | 0.00 | 0.00 | 0.87 | 1.00 | −0.11 | ✗ (tie@0) | ✗ |
| noisy_realistic | 0.37 | 0.10 | 0.50 | 0.90 | 0.50 | ✓ | ✗ |
| anti_delivery (adv) | 0.00 | 0.33 | 0.97 | 0.67 | 0.02 | ✗ | ✗ |
| delivery×sens-gate (adv) | 0.73 | 0.10 | 0.63 | 0.90 | 0.88 | ✓ | ✓ (fragile/circular) |
| subtype_blindspot (adv) | 0.33 | 0.33 | 1.00 | 0.67 | 0.35 | ✗ (tie) | ✗ |
| potency_dominant2 (adv) | 0.00 | 1.00 | 1.00 | 0.00 | 0.63 | ✗ | ✗ |

\* multiplicative is degenerate (one drug best for all 30); both sim=1.0 and ML=1.0 are
majority-class artifacts, not evidence of personalization.
