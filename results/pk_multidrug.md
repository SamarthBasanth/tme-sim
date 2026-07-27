# Multi-drug physics reality check — real cited PK anchors (no synthetic data)

**Module:** `tmesim/validation/pk_multidrug.py` (run: `python -m tmesim.validation.pk_multidrug`)

**Scope.** Extends the single-drug doxorubicin check (`pk_reality.py`) and the one-constant calibration (`calibration.py`) to **six real drugs**, using only the PMID-verified measurements curated in `pk_anchors.py`. Two physically distinct tests are run and **never mixed**: a penetration-**length** test (compared to converted perivascular half-distances) and a BBB/partition-**ratio** ordering test (rank concordance only — a dimensionless `bbb_penetrance` score is not a partition ratio). Nothing in core, the engine, or the molecular layer is modified.

---

## Test 1 — penetration length generalizes out-of-sample

The single mechanistic constant `alpha` (`k_eff = k_passive + alpha·uptake_proxy`) is fit **once** on doxorubicin's **e-fold-equivalent** anchor — the Primeau half-distance 40–50 µm converted to an e-folding length `L = d_half/ln2 ≈ 58–72 µm` (Primeau et al., *Clin Cancer Res* 2005;11(24):8782–8, PMID 16361566, DOI 10.1158/1078-0432.CCR-05-1664). Result: **alpha = 148.6**. With alpha **frozen**, every other drug's `L = sqrt(D_eff / k_eff)` is a pure prediction. Sim `L` is compared to the **converted (L-equivalent)** band `d_half/ln2`, per the `pk_anchors` convention — not to the raw half-distance.

| drug | role | predicted L | real half-dist | L-equivalent band | ratio to mid | verdict | PMID |
|---|---|---:|---:|---:|---:|:--:|---|
| doxorubicin | in-sample (fit) | **64.9 µm** | 40–50 µm | 58–72 µm | 1.00× | **FIT** | 16361566 |
| mitoxantrone | out-of-sample | **48.2 µm** | 25–75 µm | 36–108 µm | 0.67× | **PASS** | 23680920 |
| topotecan | out-of-sample | **83.6 µm** | 25–75 µm | 36–108 µm | 1.16× | **PASS** | 23680920 |

Out-of-sample anchors (both from the same immunofluorescence-vs-distance study): Patel KJ, Tredan O, Tannock IF, *Cancer Chemother Pharmacol* 2013;72(1):127–38, PMID 23680920, DOI 10.1007/s00280-013-2176-z.

**Per-drug verdict.**
- **doxorubicin** — *in-sample fit* (defines alpha). Post-fit `L = 64.9 µm`, on the 58–72 µm target by construction. Not a test.
- **mitoxantrone** — PASS. Predicted `L = 48.2 µm` vs L-equivalent 36–108 µm (0.67× the midpoint). Class: anthracenedione DNA-intercalator (SAME strongly-bound class as doxorubicin -- a magnitude-transfer test, not a maximally independent one). Anchor confidence: medium — The ~25-75 um half-distance is reported as a shared range for doxorubicin/mitoxantrone/topotecan at 10 min, not as a mitoxantrone-specific single value; treat as order-of-magnitude. Perivascular tissue gradient (comparable to sim penetration length).
- **topotecan** — PASS. Predicted `L = 83.6 µm` vs L-equivalent 36–108 µm (1.16× the midpoint). Class: camptothecin topoisomerase-I inhibitor (chemically DISTINCT from the anthracycline fit -- a more independent test). Anchor confidence: medium — Shared ~25-75 um half-distance range for the three drugs at 10 min; not a topotecan-specific point value. Perivascular tissue gradient.

**Out-of-sample tally: 2/2 drugs land in (or within 3× of) their real cited bands.**

---

## Test 2 — bbb_penetrance ordering vs real partition ratios

The molecular layer's structure-only `bbb_penetrance` (CNS-MPO-style rules on MW/logP/TPSA/HBD) is a dimensionless 0–1 score, **not** a partition ratio, so only its **rank ordering** is compared to the real measured brain/CSF:plasma ratios. Real ordering (ascending penetration): **methotrexate < cisplatin < temozolomide**.

| drug | real ratio | quantity (measured) | bbb_penetrance | real rank | bbb rank | PMID |
|---|---:|---|---:|:--:|:--:|---|
| methotrexate | 0.013 | csf_to_plasma_ratio | 0.115 | 1 | 1 | 3855310 |
| cisplatin | 0.025–0.029 | csf_to_plasma_ratio | 0.893 | 2 | 3 | 8548876 |
| temozolomide | 0.178 | tissue_to_plasma_ratio | 0.558 | 3 | 2 | 19861433 |

Real partition anchors: temozolomide — Portnow et al., *Clin Cancer Res* 2009;15(22):7092–8, PMID 19861433 (brain-ECF:plasma AUC 0.178); methotrexate — Morse et al., *J Clin Oncol* 1985;3(1):19–24, PMID 3855310 (CSF:plasma 0.013); cisplatin — Nakagawa et al., *Cancer Chemother Pharmacol* 1996;37(4):317–26, PMID 8548876 (peak CSF:plasma total-Pt 0.025–0.029).

**Concordance:** real ordering `methotrexate < cisplatin < temozolomide` vs sim ordering `methotrexate < temozolomide < cisplatin` → Spearman ρ = 0.50, Kendall τ = 0.33 (2/3 pairs concordant).

The ordering is **partially** concordant. The model correctly ranks the least-CNS-penetrant agent (**methotrexate**) lowest, but **inverts** temozolomide vs cisplatin: it over-ranks cisplatin (a small, low-TPSA complex that scores high on physicochemical rules but is in reality hydrophilic / poorly CNS-penetrant) and under-ranks temozolomide (a known CNS-MPO under-prediction). Structure-only rules cannot see the active efflux, ion-trapping and carrier transport that set the true CNS partition.

---

## Overall assessment

PENETRATION (one constant alpha=148.6, fit ONCE on the doxorubicin e-fold-equivalent 58-72 um anchor): out-of-sample 2/2 drugs land in / within-3x of their real cited bands. mitoxantrone L=48 um vs L-equiv 36-108 um (0.67x mid) [PASS]; topotecan L=84 um vs L-equiv 36-108 um (1.16x mid) [PASS]. The single calibrated constant GENERALIZES across the perivascular-anchor drugs at the order-of-magnitude level. PARTITION (bbb_penetrance ORDERING vs real ratio ordering): real methotrexate < cisplatin < temozolomide vs sim methotrexate < temozolomide < cisplatin; Spearman rho=0.50, Kendall tau=0.33 (2/3 pairs concordant). bbb_penetrance ordering is PARTIALLY concordant -- it correctly ranks the least-penetrant agent lowest but INVERTS the pair(s): temozolomide vs cisplatin. This is the known structure-only limitation (CNS-MPO physicochemical rules over-rank cisplatin and under-rank temozolomide; they cannot see active efflux, ion-trapping or the carrier transport that set the true CNS partition).

### Honesty caveats

- **One anchor sets the entire penetration scale.** alpha is fit to a single doxorubicin number; every drug's absolute L rides on that one measurement.
- **The out-of-sample perivascular band is SHARED, not drug-specific.** PMID 23680920 reports the ~25–75 µm half-distance jointly for doxorubicin/mitoxantrone/topotecan at 10 min (both anchors: medium confidence). Landing in the band is an order-of-magnitude agreement, not a fine per-drug discrimination — the real biological fact is that all three have similarly steep tens-of-µm gradients, which the model reproduces without refitting.
- **Mitoxantrone is same-class as the fit drug** (both anthracenedione/anthracycline DNA-intercalators): a magnitude-transfer test within the strongly-bound class. **Topotecan** (camptothecin topo-I inhibitor) is chemically distinct and therefore the more independent penetration test.
- **`bbb_penetrance` is a low-confidence structure-only ordering**, not a partition ratio; only ranks are compared, and even those are only partially concordant (ρ = 0.50). It captures the extreme (methotrexate poorest) but not the cisplatin/temozolomide contrast.
- **No refitting, no synthetic data.** Every anchor is a real, PMID-verified measurement from `pk_anchors.py`; the only fitted quantity is the single alpha, fit once on doxorubicin.
