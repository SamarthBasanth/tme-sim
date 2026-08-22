# Axis inventory — what has been swept, what has not, and which gap breaks a claim next

**This is an inventory, not an experiment.** No sweep was run, no default moved, nothing
under `tmesim/`, `docs/`, `app/`, `landing-3d/` or `data_cache/` was touched. Every number
below is transcribed from a file named beside it. Where a document is ambiguous about
whether something was swept, this file says so rather than deciding.

---

## 0. Why this file exists

Codon publishes an unusually good list of its own limits — IFP 41× low, dx = 12 µm not
converged, R² an interval, the spatial spread 1.08–1.21× not 2.0×, absolute kill not an
output, the full compound ordering not robust. `results/CLAIM_INVENTORY.md` §1d alone
tabulates nine failed absolutes, and its §1e counts *"roughly 18 distinct FAIL rows across
`results/`"*.

On 2026-08-21 an axis that was on none of those lists broke the central ranking claim.
`results/brain_slice_anchor.md` measured that the engine's implied tissue capacity is
**1.0 mL·g⁻¹** against a measured rat brain-slice range of **3.2 – 4184** (n = 33), and
`results/capacity_ranking_robustness.md` §0 then showed that under this repo's own
structural `K_b` rule the **breast top-1 pick moves from `vinorelbine` to `doxorubicin`**,
and that even a *uniform* capacity breaks the published order at **R = 36.3** (GBM 10-drug)
and **R = 5.28** (breast 10-drug).

Capacity was not ruled out and then found wanting. It was **outside the swept set** for as
long as the ranking claim existed, because nobody had picked it. A reader of the published
limit list currently cannot tell the difference between *"we tested this and it held"* and
*"we never looked."* This file is the difference.

**The honest epistemic position, stated once:** on this evidence the remaining claims
should be read as **not yet refuted**, not as *verified*. The next retraction most likely
comes from an axis nobody has named.

---

## 1. How to read this file

Three categories, and one test that decides which one an axis lands in.

* **SWEPT** — varied over a stated range, in a named file, with a recorded outcome.
* **PARTIALLY SWEPT** — swept on one panel, one tissue, one grid or one arm but not
  another; or swept in *magnitude* but never in the *per-compound / per-tissue* direction
  that matters. **This is the category capacity actually lived in**, and it is where the
  next failure most likely lives too.
* **NOT SWEPT** — with the status kept honest and unblurred:
  * **BLOCKED** — considered, registered, and stopped by absent data. Evidence required.
  * **NAMED, NOT PROBED** — the repo identifies it in prose as assumed/uncalibrated, and
    then never varies it. Evidence required.
  * **NEVER CONSIDERED** — no document, docstring or registration in this repo treats it
    as an axis at all. **This status is not softened anywhere below.** Where the honest
    answer is that nobody thought of it, it says so.

### 1.1 The test that matters more than the category

Capacity was predictable in hindsight, and the reason generalises. Sorting axes by "is the
knob per-compound?" is the wrong sort — a *uniform* capacity broke the order too. The
right question is:

> **Does the knob's ACTION factor through a per-compound (or per-tissue) quantity, even
> when the knob itself is a single number?**

Capacity is one number `R`, but `results/capacity_ranking_robustness.md` §2 shows it acts
through `k·T/R`, and the panel's `k` spans **27×** (0.000501 – 0.0134 s⁻¹). That is why
"a constant error cancels out of a ranking" failed on this axis, and §10 of that file says
so outright.

The same test explains the two other places a nominally common-mode knob moved a verdict:
the grid (each drug has its own boundary layer `L`, so refinement acts per-drug — only
docetaxel's verdict is a function of dx) and the normaliser (`tail90` moved docetaxel from
+1.5 % to **+35.2 %** outside its ceiling, `results/tail90_falsification.md` §7.1–7.2,
while leaving eight other verdicts alone).

**Section 7 ranks the unswept axes by this test, and that is the deliverable.**

---

## 2. The pattern the inventory exposes, before the tables

Reading the swept set and the claim set together produces one finding that neither
produces alone.

**Every axis whose action factors through a per-compound quantity has moved a ranking or
verdict claim, essentially every time it has been swept. Every axis whose action is
genuinely common across compounds has not.**

| axis swept | action | outcome |
|---|---|---|
| per-drug IC50 | per-compound | **moved** — the three largest destabilisers of P(cisplatin rank 1), swings 0.161 / 0.157 / 0.108 (`ranking_robustness.md` §5) |
| transvascular delivery law | per-compound (delivery changes up to **9.58×**, 5-FU) | **moved** — GBM top-1 cisplatin → vincristine; P(rank 1) 0.9277 → 0.3965 (`CLAIM_INVENTORY.md` C3, `transvascular_exchange_revalidation.md` §4.1) |
| delivery source geometry | per-compound | **moved** — cisplatin↔temozolomide swap at rank 1 in *both* raw and dose-matched arms; **four** 3-D:2-D verdicts move (`delivery_source_geometry.md` §5.1–5.2) |
| tissue capacity | per-compound *and* uniform, via `k·T/R` | **moved** — breast top-1 flips; full order breaks at R = 36.3 / 5.28 (`capacity_ranking_robustness.md` §0) |
| anchor drug (dox vs mitoxantrone) | per-compound | **moved** — anchor-exempt arm flips top-1 across a 0.03 % margin (`alpha_identifiability.md` §2.1) |
| section × state route | per-tissue | **moved** — 2 distinct orderings; `brain_319` swaps methotrexate/temozolomide (`flagship_markers.md` §4) |
| calibration α | **pure gauge** — proven | held; swing **0.007**, indistinguishable from noise at n = 40 (`ranking_robustness.md` §5; gauge proof `alpha_identifiability.md` §1.2) |
| ECM tortuosity | common | held; swing 0.005 |
| KDE bandwidth | common | held; swing 0.005 |
| endothelial typing threshold | common | held; swing 0.002 |
| oxygen level | common (drug-independent multiplier) | coverage **exactly** invariant; ranking invariant in measurement (`oxygen_level_shift.md` §4) |
| cell-state route | common | ranking invariant across all three routes, on one section (`cell_state_impact.md` §3.3) |
| exposure unit (3.36×) | common | decision probabilities move ≤ **0.004** (`ranking_robustness.md` banner) |

**One per-compound axis was swept and the ranking readout was never taken.**
`logp_provenance.md` §5.2 varies each drug's logP by ±1 log unit and reports that `L` moves
0.75×–1.36× per drug, with the largest composite move (`bbb × L`) at **tamoxifen 6.6×** — a
per-drug spread in delivery larger than several perturbations that *did* flip a top-1 pick.
It reports no re-ranking. So the count is: **seven per-compound axes swept, six measured
against a ranking, six moved it, and the seventh was not measured against one.** There is no
per-compound axis in this corpus that was swept, scored on an ordering, and held.

The published robustness claims — `CLAIM_INVENTORY.md` C1, C2, C4 — are **all** drawn from
the lower block. C2 names its own falsifier correctly (*"a drug-specific transport
parameter"*) and C3 records that the falsifier was found. Capacity is the second instance
of the same falsifier class; it is not a new kind of failure.

**The sharpest single instance.** `results/alpha_identifiability.md` §1 proves that of the
11 parameters in the sink calibration, exactly **1 is identified** (`k_dox`), **9 are
not** (the ratios `b_d/b_dox`, *"zero data constraint"*), and **α itself is a gauge
coordinate** — §1.2 scales every drug's binding by a factor `t`, moves α over four orders
of magnitude (1.4861 → 14860.82), and every downstream number is identical to **2.8e-14
µm**. `results/uncertainty.md` and `results/ranking_robustness.md` then sample **α** with
the widest prior of the four (0.5×–2×).

> The uncertainty quantification swept the one direction that is provably unable to change
> a ranking, and left untouched the nine directions that can — nine numbers that
> `results/binding_sink_refuted.md` §2 shows are *"refuted in direction, not merely in
> magnitude."*

That is not a criticism of the α sweep, which is correct and worth having. It is the
clearest available statement of what "six published limits" is and is not evidence of.

---

## 3. SWEPT

Compressed. Each row: what was varied · range · file §section · outcome. Fuller detail
lives in the cited file; nothing here restates a number that file does not print.

### 3.1 Numerics and discretisation

| axis | varied over | file §section | outcome |
|---|---|---|---|
| Grid `dx`, GBM kill panel | 30/24/18/16/12/10/8/7 µm (dx = 6 refused, CFL 1.129) | `grid_convergence_gbm.md` §0, §2, §3 | **No verdict moved.** 3/3 PASS at every rung; discretisation error at dx = 12 ≤ 1.76 %. Ladder non-monotone; observed orders 0.96 and −3.33, fits flagged DEGENERATE. |
| Grid `dx`, breast kill panel | 8/10/12/16/18/24/30 µm | `grid_convergence_breast.md` §0, §2, §3 | **One verdict moved.** Docetaxel flips on the ladder; Richardson limit 0.1537, power-law 0.1573, both FAIL. 3 of 4 scorable drugs constant. |
| Timestep `dt` | 10 → 5 s at dx = 8 and 12 | `grid_convergence_gbm.md` §7.3; `cfl_resolution.md` §4 | Ratio moves **0.0003 %**. Licenses nothing about absolute AUC or kill. |
| CFL / interstitial speed vs `dx` | dx = 15/12/10/8/7/6; CFL 0.472 → 1.129 | `cfl_resolution.md` §0, §2 | Retracted the claim that dt = 10 s violates CFL at dx ≤ 8. GBM floor dx = 7 (0.955). |
| Boundary-layer resolution `δ/dx` | 54 fitted points on an exact slab, exchange ON and OFF | `boundary_layer_resolution.md` §0, §4.1–4.2 | Boundary layer **refuted** as the cause of verdict flips; `L = √(D_eff/k)` ranks movers at ρ = −0.94, `δ` at +0.47 (wrong sign). |
| Dimensionality 2-D vs 3-D | λ/L ∈ {0.5,1,2,4,8}; r ∈ {25…200 µm}; 18 network configs | `dimensionality.md` §1.3–§4 | 2-D needs no fix; `L` and mean tumour AUC exactly preserved. **Coverage** moves up to 0.205 absolute — the one quantity that materially moves. |
| Sub-grid closure | ON/OFF × 3 oxygen configs × dx ∈ {30,18,12,8,6,5} | `subgrid_closure.md` §5.1 | Works on 1 of 3 configs, **fails with a sign flip** on another, **destructive on the shipped one**. Ships default-off. |
| Grid, Visium reconstruction | dx = 100/50/30/20/12/8 µm | `cross_tissue.md` | dx = 12 adequate for Visium (kill moves 1.3 %). Earlier cross-tissue spreads **retracted** as resolution-confounded. |
| Grid, held-out patient | dx = 12 vs 24 | `heldout_patient_validation.md` §5.1–5.2 | **10/10 verdicts identical.** |
| Grid, immune-sink structure | dx = 12 vs 20 | `immune_sink_impact.md` §5.5 | Structure ratio swing 1.017× — a property of tissue, not mesh. |

### 3.2 Normaliser convention

| axis | varied over | file §section | outcome |
|---|---|---|---|
| Normaliser family | 12 candidates (p50…max, tail90/95/99, mean, rms, exposure-weighted mean) × dx = 1–30 µm | `p99_convention.md` §3.1, §3.3, §3.6 | **p99 is the worst of the twelve** on grid error at dx = 12 (+6.98 % vs tail90 −0.08 %). tail90 recommended; **not adopted**. Verdict set unchanged at dx = 12. |
| Scale invariance | λ ∈ [1e-6, 1e6], 8 factors × 9 drugs | `p99_convention.md` §0.1 | Exact to 6.7e-16 for every degree-1-homogeneous normaliser. The protection is not from the "99". |
| Perfect-normaliser counterfactual | 66 rungs, 9 drug/tissue pairs | `tail90_falsification.md` §0, §4, §5 | **Split.** Confirmed on breast, **REFUTED on GBM**. Docetaxel's FAIL becomes decisive (+1.5 % → **+35.2 %**). |
| The choice of "99" itself | — | `p99_convention.md` §0.3, §2 | *"was never justified… an undocumented constant sets the robustness of every scored verdict in the repo."* |

### 3.3 Oxygen

| axis | varied over | file §section | outcome |
|---|---|---|---|
| `consumption` × `sink_form` | 0.0002 → 0.5 (2,500×, 3.4 decades) × {linear, Michaelis-Menten} | `oxygen_sensitivity_band.md` §2–3 | **The headline spans its entire range**: hypoxic fraction exactly 0.0 on one literature bracket, exactly 1.0 on the other, for both forms. |
| MM half-saturation `K` | 0.010 → 0.500 | `oxygen_sensitivity_band.md` §4 | Hypoxic fraction 1.0 → 0.0 across the literature range. |
| Sink at derived `supply = 194` | consumption ∈ {0.02, 0.5, 2.5, 5, 10, 20, 50} | `oxygen_supply_calibration.md` §5.1 | consumption = 10/20 brackets Lally's measured 5.1 mmHg; `F_state` 1.0000 → 0.3655. |
| Oxygen config × `dx` | dx = 30/18/12/8 × 4 configs | `oxygen_joint_grid_test.md` §1.1–1.4 | **A claim moved with the mesh** — the three-way agreement holds at 18/12/8 and fails at 30. Jointly calibrated config is *more* grid-labile than shipped. |
| Hypoxia threshold, breast | 2.5 / 5.0 / 7.0 / 10.0 / 14.0 mmHg | `oxygen_calibration_breast.md` §4 | Hypoxic fraction 0.0000 % → 95.3734 %. At Vaupel's own 2.5 mmHg the model predicts **exactly zero** hypoxic breast cells. |
| Oxygen level shift | — | `oxygen_level_shift.md` §4 | Coverage **exactly** invariant (`==`, not `isclose`). Relative kill changes span 71.4–77.0 % — near-common, not common. |

### 3.4 Cell state

| axis | varied over | file §section | outcome |
|---|---|---|---|
| State route × mask × dx | 3 routes × 2 masks × dx ∈ {30,20,12}, 3 sections | `cell_state_impact.md` §2.1–2.3, §3.3 | TMZ kill 0.0184–0.0315 (1.71×); ceiling 1.74×. **Cross-drug ranking invariant across routes — on one section.** |
| `F_state` vs dx | 7 → 30 µm | `grid_convergence_gbm.md` §5 | Span 0.0461; hypoxic 0.0000 at every rung. |
| `F_state` mechanism | source rasterisation vs PDE vs sampling | `fstate_mechanism.md` §0, §2–6 | Driver is **vessel-source rasterisation** (cov/var +1.057), not the PDE (−0.056). |
| Classifier aggregation rule | `mean_z`, `max_z`, `marker_vote`, `mean_z_norm`, `matched_null_z` | `classifier_fix.md` §2–3 | Arity bug fixed; endothelial cross-section spread 2.91× → 1.49×. **Gate held fixed at 1.0 for every candidate.** |
| Ranking across sections × routes | 3 sections × 2 routes | `flagship_markers.md` §4 | **Ranking NOT invariant** — 2 orderings; the swapping pair's gap is 1.5σ. |

### 3.5 Calibration, α, transport parameters

| axis | varied over | file §section | outcome |
|---|---|---|---|
| Anchor binding prior `b_dox` | 6 factors spanning 100× | `alpha_identifiability.md` §1.1 | α moves exactly 1/f; `α·b_dox` invariant to 2.1e-16. Only `k_dox` is identified. |
| Global binding gauge `t` | α 1.4861 → 14860.82 (4 orders) | `alpha_identifiability.md` §1.2 | Every downstream number identical to 2.8e-14 µm. **α is a gauge coordinate, not a measurement.** |
| Anchor **drug** | doxorubicin vs mitoxantrone | `alpha_identifiability.md` §2, §2.1 | α 148.61 → 66.29; every `L` moves 1.48–1.50×. Anchor-exempt arm flips top-1 across 0.03 %. |
| `D_tissue`, breast | 3.0e-07 → 1.6e-06 cm²/s (5.3×), `frozen` and `recal` | `dtissue_sensitivity.md` §0, §3–4, §7.3 | 0 flips at dx = 24 in both modes; **1 flip at dx = 12 in `frozen`** (docetaxel FAIL → PASS). §8: *"D_tissue and the perivascular sink strength are not separately identified."* |
| `φ_cell` (hindrance) | 0.40 / 0.50 / 0.60 / 0.72 / 0.7635 / 0.7892 / 0.85 | `cellular_exclusion.md` §4.3 | Fold 2.00× → 9.50×; D_eff 254.7 → 53.6 µm²/s. **Now the dominant uncertainty in the hindrance law.** |
| `φ_cell` scored on Pruijn | 0.40–0.85 | `two_route_transport.md` §4.5 | R² flat at 0.71–0.76. Adopted φ is neither the R² optimum (0.50) nor the bias optimum (0.85). |
| HBD decrement magnitude | 0.000 / 0.579 / 0.941 / 1.000 / **1.340** / 1.412 / 1.774 / 2.173 / 2.245 | `two_route_efflux.md` §3.3 | R² 0.588–0.722 across the whole range. Optimum at 1.000 (0.725); adopted 1.340 gives 0.718. |
| HBD **counting convention** | `NumHDonors` vs `NHOHCount` | `hbd_convention_audit.md` §0, §4, §6 | **α moves 32 %**; **0 of 50 verdicts flip.** |
| `logP` per drug | measured logD endpoints; ±1 log unit for 4 unsourced drugs | `logp_provenance.md` §5.2 | `L` moves 0.75×–1.36× per log unit; largest composite move tamoxifen 6.6×. **No re-ranking reported** — see §2. |

### 3.6 Vascular source and geometry

| axis | varied over | file §section | outcome |
|---|---|---|---|
| Vessel linking distance | 5/7.5/10/15/**20**/25/30 µm (3 sections); 10–40 µm elsewhere | `endothelial_density.md` §4; `vessel_network_impact.md` §1.3; `heldout_patient.md` §2.3 | Endothelial count identical across the whole sweep; **tumour→vessel distance does not move at all**. MVD moves 1.8×. |
| Discrete vs KDE source × dx | dx = 15/12/8/6 | `vessel_promotion.md` §2.1–2.2 | Discrete converged; **KDE fails at the finest step** (6.35 % move, outside the 5 % tolerance). |
| Vessel painting geometry | binary / area-conserving / area_matched × dx = 30…5 | `area_conservation_promotion.md` §1.1–1.2 | `F_state` drift **98.82 % removed**; **zero verdicts moved** on promotion. |
| Off-lattice `dx` | 45/22.5/11.25/7.5/3.75 µm | `vessel_source_conservation.md` §2.1 | Area-conserving 65× better off-lattice. |
| Delivery source geometry | binary / subgrid / subgrid_matched, 9 rows | `delivery_source_geometry.md` §5.1–5.2 | **Four verdicts move.** α is *exactly* degenerate to this axis (`alpha_refit − alpha_frozen = 0.0`). Recommendation: do not promote. |
| Vessel percentile threshold | p90/p95/p98/p99 | `vessel_profile.md` (legacy-metric sensitivity) | Completely flat at all four. |

### 3.7 Pharmacology mechanisms (all swept to a default-OFF verdict)

| axis | varied over | file §section | outcome |
|---|---|---|---|
| Immune (TAM) sink strength | S = 0/1/2/5/10/20 | `immune_sink_impact.md` §4, §4.1 | AUC retained ~0.52× uniformly; kill retained differs by drug (0.47–0.79×). **Cross-drug ranking unchanged.** Ships default-off. |
| Lysosomal `φ_lys` | 0.000 → 5.270 % (14 points) + 5 named configs | `ion_trapping_impact.md` §4.1–4.2 | Monotone, **no interior optimum**. Best is 3.00× at φ = 0, physically unattainable. Refuted; default OFF. |
| Efflux discretionary choices | form × observable × placement (2×2×2) | `efflux_impact.md` §7 | **0 of 8 configurations improve the scorecard.** Default OFF. |
| Intracellular retardation `K_b` | tubulin 8–50 µM × two Kd sources (~30× apart) | `intracellular_retardation.md` §2.5, §4.1–4.2 | Matched-MW rank agreement **77.4 % → 51.6 %**; R² 0.718 → 0.606. Refuted. |

### 3.8 Dose, schedule, ranking, tissue

| axis | varied over | file §section | outcome |
|---|---|---|---|
| Dose multiplier × fractionation | 1/3/10/30× × {1,3,10} doses × both delivery laws | `dose_schedule.md` §4.5, §8 | No verdict moved. Under the corrected law the kill readout is at the resolution limit of 23,812 cells. |
| IC50 prior scenarios | A / B / C, 10,240 re-rankings each | `ranking_robustness.md` §0, §4.1 | P(cisplatin rank 1) 1.000 / **0.931** / 0.587. **Top-1 decision-grade in B; top-2 SET is not** (0.563). |
| Per-input quartile destabilisation | 9 inputs | `ranking_robustness.md` §5 | IC50 terms dominate (0.161 / 0.157 / 0.108); the four transport knobs are 0.002–0.007. |
| 4-knob Monte Carlo | endothelial threshold, KDE bandwidth, α, ECM tortuosity | `uncertainty.md` §1–4 | kill 90 % CI width **129 %**; α dominates variance (S₁ 0.78). `grid_dx_um` deliberately **not** sampled. |
| Tissue / patient | 5 real tumours at dx = 12 | `cross_tissue.md` | Kill 1.5 % → 46.3 % (31.6×). Heterogeneity claim survives; **all earlier numbers retracted**. |
| Second GBM section | s1 vs s2 | `heldout_patient.md` §3.3, §6.1–6.4 | Ordering replicates; absolutes move a lot. |
| Uniform tissue capacity | R = 1 → 4184, 22-point grid | `capacity_ranking_robustness.md` §3, §9 | **Full order breaks at R = 36.3 (GBM-10) and 5.28 (breast-10).** Top-1 holds across the swept grid. |
| Per-compound capacity | 40,000 draws + 4 named assignments | `capacity_ranking_robustness.md` §5, §8 | **Breast top-1 flips.** No panel clears the 0.90/0.10 bar. |

---

## 4. PARTIALLY SWEPT

**This is the category capacity lived in and the one to read most carefully.** An axis
here has a sweep attached to it, which is exactly what makes it easy to mistake for
closed.

| axis | swept where | **not** swept where | file §section |
|---|---|---|---|
| **Grid `dx`** | GBM panel, 8 rungs | dx = 12 was established on **GBM only** and inherited by breast. Breast resolves docetaxel's `L` = 28.2 µm with ~2.4 cells vs GBM's 5.4. *"dx = 12 µm cannot be justified panel-agnostically."* A full **joint (dx, dt)** refinement is *"not made here"*; grid dependence at other reporting concentrations *"not characterised"*; a joint **convention × mesh** sensitivity *"is not measured"*. | `grid_convergence_gbm.md` §6, §8.1, §8.4, §8.5; `grid_dependence_audit.md` §1.4 |
| **Normaliser** | 12 candidates on a slab | `tail90` was **never measured on the tissue ladders** — it was *bounded* by a counterfactual. The 25 stored per-cell exposure vectors are **all at dx = 12**, so the claim that the test costs no new solves is stated as wrong. | `tail90_falsification.md` §0.1, §9.1 |
| **`D_tissue`** | 10 points, 5.3×, at dx = **24** | The published grid is dx = 12, where the reduced 2-point set gives **1 flip**. The §4 flip map is explicitly not a statement about the published verdict. | `dtissue_sensitivity.md` §4 header, §7.2–7.3 |
| **`φ_cell` / `cell_radius_um`** | 0.40–0.85 magnitude sweep | Both are measured on **one Xenium breast section** and applied unchanged to GBM and to a V79 hamster MCL. *"the packing density is imported from the wrong system"*; *"A literature search for an MCL extracellular-volume measurement returned nothing usable."* Never given a per-tissue value. | `cellular_exclusion.md` §3.2; `anchor_registration.md` §5; `two_route.py:123-133` |
| **HBD decrement `1.34`** | magnitude 0.000–2.245; counting convention 2 arms | Never swept **per donor group** — the module's own docstring lists real Abraham `A` values from 0.16 (aliphatic amine) to 0.62 (carboxylic acid), a ~3.9× spread collapsed to one number, then multiplied by a per-compound HBD count of 1–6 across the panel. Neither sweep scored a **ranking**. | `two_route.py:218-250`; `two_route_efflux.md` §3.3; `hbd_convention_audit.md` §6; `descriptors.py:52-61` |
| **`bbb_penetrance` in breast** | replaced by a measured tissue-level PS for **2 drugs**, against Kp anchors | Never applied to the other 8 breast compounds and **never re-ranked**. The measurement says the correct replacement is *drug-independent*: measured breast 5-FU:methotrexate permeability ratio **1.02×** against **28.1×** from `bbb_penetrance`. | `breast_tissue_permeability.md` §1–2 |
| **Immune sink strength** | 6 strengths on the flagship section | Swept on **one section only**, although the immune field is measured on all three — and the other two carry **3.9×** the mean immune area fraction, so the effect would be larger there. | `immune_sink_impact.md` §4 |
| **Vessel linking distance** | 7 values on section 1 | Production 20 µm was fixed on section 1; **section 2's published MVD band is only reached at 25–40 µm**, and it was not re-tuned. | `heldout_patient.md` §2.3 |
| **Classifier gate** | 5 aggregation rules | The gate **stayed fixed at 1.0 for every candidate**. `ranking_robustness.md` §3.1 samples the *endothelial* threshold only; TUMOR / STROMAL / IMMUNE gates are never varied — and TUMOR is the population every kill fraction is averaged over. | `classifier_fix.md` §0, §2; `ranking_robustness.md` §3.1 |
| **Hypoxia banding** | the lower edge `t` swept 2.5–14 mmHg on breast | The **band-width factor 2.0** (`agents.py:83-84`, `QUIESCENT` on `[t, 2t)`) has never been varied independently — sweeping `t` moves both edges together. `kill_calibration.md` §220 records that **90.0 % of tumour cells sit inside that single band** and that `F_state` multiplies directly into the scored ratio. It is not a config field. | `oxygen_calibration_breast.md` §4; `cells/agents.py:82-84`; `kill_calibration.md` §220 |
| **Panel composition** | GBM 5-drug vs 10-drug; breast 10-drug | Never treated as an axis in its own right, yet it decides the answer: the 5-drug GBM sub-panel reproduces the published order under **every** named capacity assignment while the 10-drug panel does not. | `capacity_ranking_robustness.md` §5.2–5.3 |
| **Oxygen config × dx** | dx = 30/18/12/8 | The **scored ratio at dx = 4–6 µm was not measured**; median pO₂ drift below dx = 8 also unmeasured. | `oxygen_joint_grid_test.md` §1.4, §6 |
| **Cross-panel imputation** | 19 targets, LOSO | `proliferation_TOP2A` / `proliferation_composite` **pass on `gse_366` and fail on `brain_319`** — a section-dependent pass. | `panel_gap_imputation_sweep.md` §2, §4.3 |
| **Held-out independence** | a second section exists | Three breaches: the typing rule was selected with s2 in the criterion; the IO marker map was written against s2; **the hypoxia half is not out-of-sample at all, by construction.** | `heldout_patient_validation.md` §2.2–2.4 |

---

## 5. NOT SWEPT

### 5.1 BLOCKED — considered, registered, stopped by absent data

These are real scoping decisions with evidence, and they are labelled as such because the
evidence exists. They are the *smallest* category.

| axis | status | file §section |
|---|---|---|
| **`cell_partition` (K = 1.0)** | Registered as TEST C and **NOT RUNNABLE**: only 4 scoreable brain `Kp,uu,cell` values spanning 54×; runnable at n ≥ 20. *"a hardcoded, load-bearing input to the transcellular route — REMAINS UNCONSTRAINED BY ANY MEASUREMENT."* Explicitly not constrained by TEST B either. | `transport_anchor_prereg.md` §3.2, §3.8; `brain_slice_anchor.md` §8; `capacity_ranking_robustness.md` §10 |
| **Second cell line (TEST A)** | **BLOCKED, registered anyway** — *"Pruijn 2008 does not publish per-compound `Dmcl` for SiHa, HCT8-Ea or HCT8-Ra… There is nothing to digitize."* | `transport_anchor_prereg.md` §5.1 |
| **Tumour `φ_lys`** | *"NOT measured: tumour tissue, on any basis."* All three tumour tissues carry `phi_lys=None`. The whole surface was swept **because** no value is defensible. | `ion_trapping_impact.md` §2, §4.2; `ion_trapping.py:235-271` |
| **Intracellular tubulin concentration** | *"UNVERIFIED, and it is the blocking input"* — no primary mammalian µM value retrievable. | `intracellular_retardation.md` §2.5 |
| **Efflux magnitude in tissue** | *"the literature will not supply the number"* — mechanism works, cannot be applied. | `two_route_efflux.md` §2.4 |
| **Solver output fields** | *"nothing public was found that measures a spatially resolved free-interstitial concentration field in tumour tissue this engine could reconstruct."* | `anchor_registration.md` §7; `external_anchor_survey.md` §5 |
| **Response axis** | *"no open dataset anywhere pairs measured spatial drug distribution with measured treatment response in the same tumour."* | `DATA_HUNT.md`; `VERIFICATION_METHODOLOGY.md` §12d |
| **1-D slab flux geometry** | The engine *"has no 1-D slab configuration"*, so a known **directional disagreement on the engine's own panel drugs** (Tunggal 1999 / Cowan 2001) is unscored. | `anchor_registration.md` §8 row 3 |

### 5.2 NAMED, NOT PROBED — the repo says it is assumed, and then never varies it

This is the largest category, and it is the one most easily mistaken for the one above.
Saying *"this is uncalibrated"* in a docstring is not the same act as sweeping it.

| axis | what the repo says about it | where | ever swept? |
|---|---|---|---|
| **Hill coefficient `hill = 1.2`** | `_DEFAULT_HILL`, applied unconditionally to every compiled drug. *"not a per-drug fit: it is the pipeline's constant, carried through unchanged."* Listed as `assumed` — *"a constant, for every compound."* | `compile_drug.py:67, 367`; `config.py:290-292`; `bring_your_own_compound.md` §2 | **No.** No file in `results/` varies it. |
| **Per-compound plasma PK** (`plasma_peak = 50 µM`, `half_life = 4 h`, `n_doses = 1`) | Every scored panel drives **every compound with the same plasma curve**. The repo names the mismatch itself: *"free-platinum Cmax for cisplatin is roughly an order of magnitude LOWER, so cisplatin's '1×' here is supratherapeutic."* | `config.py:65-68`; `real_drug_panel.py:162`; `breast_drug_panel.py:241`; `ranking_robustness.py:395`; `dose_schedule.py:157-163` | **A common multiplier was**, on 2 drugs (`dose_schedule.md` §2). A **per-compound assignment was never made, and never re-ranked.** |
| **Nine binding ratios `b_d/b_dox`** | *"zero data constraint… not merely unmeasured, they are the ones independently shown to be wrong."* | `alpha_identifiability.md` §1; `binding_sink_refuted.md` §2 | **No.** α — the gauge direction — was swept instead. |
| **`PS_SCALE_S_INV = 1.0`** | *"THIS CONSTANT EXISTS TO MAKE A DIMENSIONAL ERROR VISIBLE RATHER THAN IMPLICIT… an unstated factor of exactly 1 s⁻¹ that nothing derives."* Calibrated value **0.16122 1/s**; residual after rescaling is exactly √8.70 = **2.950×** because *"the engine's cross-drug permeability ORDERING is wrong."* | `treatment/drugs.py:101-134`; `scorecard.py` `ps-scale` row | **No.** |
| **`_LOGP_KP_MAX = 3.0`** (membrane-partition ceiling) | *"the ceiling binds… Cisplatin (22.2 µm/s), imatinib (21.0) and methotrexate (19.6) now sit within 15 % of each other… Physically they differ by ORDERS of magnitude… no longer informative at the top."* Real Caco-2 Papp is **1–4 orders below** what this returns. | `permeability.py:66, 163-174` | **No.** |
| **`CYTOPLASM_D_OVER_WATER = 0.27`** | Measured for **one dye** (BCECF, Kao 1993) and used as the cytoplasmic ratio for every compound. Listed as frozen in every scored arm. | `two_route.py:144`; `transport_anchor_prereg.md` §2.3 | **No.** |
| **`PDConfig.hypoxia_resistance = 0.35` / `quiescent_resistance = 0.5`** | *"still UNVALIDATED PRIORS with no citation."* The literature basis found was **drug-specific**; exactly one per-drug override exists (TMZ 0.80). For `quiescent_resistance`, *"no quantitative fold-change was found anywhere in the literature for any of this project's three drugs."* | `config.py:196-198, 210-239`; `grid_dependence_audit.md` §4.1 | **No.** The state *route* was swept; the multiplier *magnitude* never. |
| **`OxygenConfig.consumption` absolute level** | *"was never calibrated against a measured oxygen tension… Nothing downstream should rest on the absolute hypoxic fraction until `consumption` is anchored to a measurement."* | `MATHEMATICAL_MODEL.md` §4.2; `oxygen.py:183-188` | **Swept in range** (§3.3) but **never anchored** — the sweep spans the entire answer, which is the opposite of a resolution. |
| **Six IFP pressure coefficients** (`_K0`, `_K_ECM_ATTEN`, `_PV`, `_LP_NORMAL`, `_LEAK_GAIN`, `_LAMBDA_LYMPH`) | *"order-of-magnitude biological values, not patient-calibrated"*; **none fitted**; *"nothing constrains the absolute scale."* | `pressure.py:61-99`; `scorecard.py` CONDITIONALS | **No — none of the six.** And IFP is the one PDE anchor with no fitted coefficients, failing 41×. |
| **Per-cell sensitivity σ_c** | *"~99 % of cells sit at σ_c ≈ 0.59 (the **arbitrary** sigmoid(0) = 0.5 midpoint)… Every absolute kill fraction is conditional on this assumption."* Was defect #5 (*"placeholder shipped as model"*). | `MATHEMATICAL_MODEL.md` §6.2; `VERIFICATION_METHODOLOGY.md` §0 | Modes swept (`uniform` vs `determinant`); the **sigmoid's own parameters** were not. |
| **`L_PARACELLULAR_LEAK = 0.238`** | *"fitted to the same contrasts it is scored on… not fully out-of-sample."* | `barrier.py:117`; `scorecard.py` CONDITIONALS | **No.** |
| **Ionisation / extracellular pH** | *"KNOWN MISSING MECHANISM — ionisation and sequestration (pH partitioning)… deliberately not implemented yet."* The membrane leg consumes a **neutral-species logP**; `formal_charge` is computed for every compound and consumed nowhere. | `permeability.py:16-39, 176-188`; `descriptors.py:91`; `transport_anchor_prereg.md` §3.4 E4 | **No.** The engine has no pH field at all. |
| **`binding.py` heuristic midpoints** | Three logistic heuristics (`plasma_free_fraction` k=3.0 s=1.5; `tissue_binding_rate`; `clearance_rate`), each tagged **CONFIDENCE: low**. *"Treat them as priors to be replaced by measured PK values, not as predictions."* | `binding.py:10-27, 85-144` | **No.** |
| **`bbb_penetrance` desirability constants** | Eight hard-coded numbers (sigmoid centres/widths for MW, logP, TPSA, HBD, plus the output logistic) defining "CNS-likeness" for every molecule. Known to mis-rank temozolomide. | `permeability.py:291-333` | **No.** The **aggregation form** was scored (arithmetic → geometric mean, defect #9); the eight constants were not. |
| **`PKA_CLASS`** | One pKa per structural class, 5 entries. *"There is no first-principles pKa predictor in RDKit, so a class table is the honest alternative."* | `ion_trapping.py:507-516` | Module ships off; not swept. |
| **`KILFORD_MG_EQUIV_PER_ML_CELL = 125.0` + Hallifax coefficients** | *"an empirical QSAR fitted to drug-like compounds binding to liver microsomes, not a measurement on any compound in this model's compound set."* **"Both primary papers are paywalled and neither was read."** | `two_route.py:544-618` | **No.** Used as-is in the capacity study's structural arm. |
| **Antibody binding-site-barrier constants** (`_DEFAULT_ANTIGEN_M`, `_DEFAULT_KON`, `_DEFAULT_KOFF`, `_TISSUE_LENGTH_UM = 150`) | Antigen density, k_on and k_off are per-target-per-antibody quantities represented by one number each. | `binding_site_barrier.py:98-116` | **No.** No antibody is on a scored panel. |

### 5.3 NEVER CONSIDERED — no document, docstring or registration treats these as axes

**Not softened.** For each of these I searched `results/`, `docs/`, the module docstrings
and the scorecard, and found no place where the quantity is presented as something that
could have been otherwise. If one of these is in fact registered somewhere I did not look,
that is a correction this file wants.

| axis | why it is an axis | where it lives |
|---|---|---|
| **The QUIESCENT band-width factor `2.0`** | `QUIESCENT` is defined as `[t, 2t)` with `2.0` written inline. It is **not a config field**, so it cannot be varied without editing the module. `kill_calibration.md` §220 records that 90.0 % of tumour cells sit in this one band and that `F_state` multiplies straight into every scored ratio. `cell_state_impact.md` opens by noting ~89 % of the TMZ kill ceiling comes from this band. Nothing anywhere asks whether the proliferating/quiescent boundary sits at exactly twice the hypoxic threshold in any tissue. | `cells/agents.py:82-84` |
| **Interstitial viscosity `1.2e-3 Pa·s` and `T = 310 K`** | One viscosity for every tissue and every compound, feeding Stokes–Einstein and therefore every free diffusivity. Tumour interstitial viscosity is not a universal constant. | `stokes_einstein.py:86-87` |
| **The MW regime threshold, and its inconsistency** | `stokes_einstein.py:75` switches the hydrodynamic-radius law at **MW ≥ 5,000**; `descriptors.py:36` uses **10,000** for the same regime question. Neither is swept, and the two disagree. The prereg separately records the MW ceiling as *"unvalidated"* between the branches, but does not name the inconsistency. | `stokes_einstein.py:75`; `descriptors.py:36` |
| **The independence assumption between the two hindrances** | *"the product is a modelling choice, not a theorem — the two hindrances are assumed statistically independent."* Named as an assumption; never varied, and no alternative composition rule is scored anywhere. | `MATHEMATICAL_MODEL.md` §2.5 |
| **Which compounds are on the panel** | Treated as given, never as a variable — yet `capacity_ranking_robustness.md` §5.2–5.3 shows the 5-drug and 10-drug GBM panels give opposite answers to the same question. Listed here rather than in §4 because no document treats panel composition as an axis at all. | `real_drug_panel.py`; `breast_drug_panel.py` |
| **The exposure statistic `E = AUC/T`** | `pkpd.py:173` maps a whole time course onto its time-average before the Hill curve, so a brief high peak and a long low plateau are indistinguishable by construction. The *window length* is discussed; the *choice of the time-average as the exposure statistic* is not. | `cells/pkpd.py:173-179` |
| **The `1e-300` floors and the zero-clamp** | Numerical guards that could in principle mask a physical negative. `subgrid.py` separately records `beta` reaching **144.0** against a validated max of 5.677 — a clamp doing load-bearing work. Never characterised as a sensitivity. | `pkpd.py:179`; `diffusion.py:455`; `subgrid.py:204-230` |

---

## 6. Sizing the unexplored space

No single ratio is honest here, because axes are not commensurable — `dx` is not one axis
in the sense that `PDConfig.quiescent_resistance` is one axis. What is available is a set
of *closed* counts, each over a set somebody else defined.

* **The pre-registration's own frozen list.** `transport_anchor_prereg.md` §2.3 names **8
  constants held fixed in every scored arm** of the transport tests: `cell_partition`,
  `HBD_LOGP_DECREMENT`, `cell_volume_fraction`, `cell_radius_um`,
  `cytoplasm_D_over_water`, membrane viscosity / thickness / logk-max, `ecm_reference`,
  `_TWO_ROUTE_DEFAULT`. **Three have been swept anywhere** (the HBD decrement, φ_cell, and
  the two-route flag itself). **Five never have.**

* **The IFP pressure model.** **Six** order-of-magnitude coefficients, **zero** swept,
  **zero** fitted (`pressure.py:76-99`). The corresponding anchor is 41× low and is the
  only PDE anchor with no fitted coefficient.

* **The uncertainty quantification.** `uncertainty.md` §1 samples **four** knobs.
  `alpha_identifiability.md` §1.2 proves one of them (α) is a gauge coordinate; §5 of
  `ranking_robustness.md` measures the other three at swings of 0.002–0.005, at or below
  the noise floor at n = 40. **The credible interval on every published kill number is
  built from one gauge direction and three near-null ones.**

* **The scorecard.** The full run scores **33 of 50** rows; **21 of the 33 never invoke
  the transport PDE**; of the rows where a *solved field* meets an *external measurement*
  there are **three, and all three FAIL** (`CLAIM_INVENTORY.md` §1e). `scorecard.py`
  carries **25 `OUT_OF_SCOPE` entries** and a `CONDITIONALS` block naming a further ~14
  unconstrained quantities.

* **The molecular compile.** Every compiled drug receives a per-compound `diffusivity`,
  `binding`, `clearance`, `bbb_penetrance`, `ic50` and `hill`. Of those six, **`hill` is
  a constant**, **`bbb_penetrance` is a brain-specific rule applied to breast**,
  **`binding` is α times a proxy whose nine cross-drug ratios are unconstrained and
  refuted in direction**, and **`clearance` is a logP/MW heuristic tagged CONFIDENCE:
  low**. `ic50` is the one that is measured, and it is the one the ranking study found to
  be the dominant destabiliser.

**The legible summary:** the numerics are thoroughly probed, the reconstruction is well
probed, the oxygen route is probed to the point of showing that it spans its own answer —
and the **per-compound pharmacology layer**, which is what a compound ranking is made of,
is almost entirely unprobed. Of the ~50 swept-axis entries tabulated in §3, the great majority are numerical,
reconstruction-side or common-mode; **seven axes with per-compound action have been swept
anywhere in the corpus, six of them were scored against an ordering, and all six moved it**
(§2).

---

## 7. Ranked: which unswept axis breaks a published claim next

Ranked by the §1.1 test — *does the knob's action factor through a per-compound or
per-tissue quantity that a published claim depends on?* — and by whether the quantity's
real spread is already known to be large.

Cost estimates are stated as estimates. **None of these was run.**

---

### #1 — Per-compound plasma pharmacokinetics
**Shape: exact match to capacity, and cheaper to test.**

* **The constant.** Every scored panel drives every compound with `PKConfig`'s defaults —
  `plasma_peak = 50 µM`, `half_life = 4 h`, `n_doses = 1` (`config.py:65-68`). Verified at
  `real_drug_panel.py:162`, `breast_drug_panel.py:241` and `ranking_robustness.py:395`,
  none of which passes a `pk=` argument; `kill_calibration.py:1036` and
  `kill_calibration_breast.py:273` pass the same nominal values to all drugs.
* **The real spread.** The repo states it itself: *"free-platinum Cmax for cisplatin is
  roughly an order of magnitude LOWER [than temozolomide's], so cisplatin's '1×' here is
  supratherapeutic"* (`dose_schedule.py:157-163`), with reference bands of 30–70 µM for
  TMZ against 2–6 µM for cisplatin.
* **Why it bites.** Under the two shipping defaults — `time_resolved = False`
  (`config.py:87`, plasma collapsed to its window mean) and `transvascular_exchange =
  False` (`config.py:137`, source linear in plasma) — a per-compound PK is **exactly a
  per-compound multiplicative scale on each drug's AUC**. That is algebraically the same
  perturbation class as the delivery-law correction, which flipped the GBM top-1 pick and
  moved P(cisplatin rank 1) from 0.9277 to 0.3965 (`CLAIM_INVENTORY.md` C3).
* **What breaks.** C2 (*"cross-drug ranking is invariant to common-mode perturbations"* —
  this is not common-mode), the top-1 call on both panels, and the dose-limited-vs-
  resistance-limited reading in `dose_schedule.md` §6.3.
* **What does NOT break, and this is the interesting half.** The 3-D:2-D kill verdicts
  (K1/K2/K3) are **immune by construction**: `kill_calibration.decompose` sets each drug's
  plasma so that its *own* p99 exposure equals `c_bath` (`p99_convention.md`, opening),
  so a per-drug PK scale is normalised away per drug. The scored panel has a defence the
  ranking panel does not — and nothing says so on any surface.
* **Why it was never done.** `dose_schedule.md` §2 sweeps a *common* multiplier
  (0.05×–30×) on two drugs and names the per-drug clinical bands in prose. The
  per-compound assignment is one step away and was not taken. Status: **named, not
  probed** — not deferred.
* **Cost.** Very low. `ranking_robustness.md` §2 already establishes that the transport
  solve and the kill mapping *"separate exactly"*. Under the shipping linear law the AUC
  scales exactly with plasma level, so a per-compound PK re-rank needs **no new PDE
  solves** — it is a rescale of stored per-cell AUC plus a re-Hill. *(This linearity is my
  inference from `config.py:137` and `treatment/drugs.py:211-216`, not a measured claim;
  it should be verified algebraically first, exactly as `capacity_ranking_robustness.md`
  §2a did for capacity.)*

---

### #2 — The Hill coefficient `hill = 1.2`
**Shape: a per-compound exponent, held constant, in a regime where it is an exponent.**

* **The constant.** `_DEFAULT_HILL = 1.2` applied unconditionally to every compiled drug
  (`compile_drug.py:67, 367`). `config.py:290-292` states it is *"not a per-drug fit"*.
  `bring_your_own_compound.md` §2 lists it as `assumed` — *"a constant, for every
  compound."*
* **The real spread.** Hill slopes are a per-compound property of a dose–response curve.
  I did **not** find per-drug Hill slopes cached anywhere in this repo, so the magnitude of
  the real spread across *this* panel is not established here — that is the first thing to
  check, and it is cheap.
* **Why it bites.** The panel operates deep in the sub-saturating regime — `config.py:214`
  records median tumour-cell exposure/IC50 = **0.056** for temozolomide — where
  `p_hill ≈ (E/IC50)^h`. `h` is therefore an **exponent on a per-drug ratio that spans
  orders of magnitude across the panel**. Small per-compound differences in `h` produce
  large, order-changing differences in predicted kill.
* **It also moves the bar, not just the prediction.** `grid_convergence_gbm.md` §3 and
  `oxygen_joint_grid_test.md` §4 both state the 3-D:2-D band *"depends only on the
  frozen literature table and on (IC50, hill)"*. So `h` moves the predicted ratio **and**
  the ceiling it is scored against, for all nine 3-D:2-D verdicts. No other axis in this
  document does both.
* **What breaks.** The ranking; potentially K1/K2/K3, which the p99 normaliser does *not*
  protect against this axis.
* **Cost.** Free of solves — the kill mapping is the cheap object.

---

### #3 — The nine unconstrained per-drug binding ratios `b_d/b_dox`
**Shape: exact, and the only axis on this list already proven wrong in direction.**

* **The status.** `alpha_identifiability.md` §1: of 11 parameters, **1 identified, 9 not**;
  the nine come from a logP/HBD sigmoid with *"zero data constraint"*, and
  `binding_sink_refuted.md` §2 shows that sigmoid is *"refuted in direction, not merely in
  magnitude."*
* **Why it was missed.** §1.2 proves α is a **gauge coordinate** — the exact direction the
  uncertainty analysis samples with its widest prior. `ranking_robustness.md` §3.1 even
  explains the consequence correctly: *"because α is shared by every drug in a draw, it
  moves all five kills together, so a wide α prior… costs the ORDER comparatively little."*
  The report identified the structure and did not follow it to the conclusion that the
  *unshared* directions were therefore unsampled.
* **What breaks.** Everything C2 protects, plus the α-conditional PASS on all nine 3-D:2-D
  rows — `anchor_registration.md` notes the kill and transvascular groups *all* cite the
  `alpha-fit` conditional.
* **Cost.** Moderate — the ratios enter the per-drug sink `k`, so this needs real solves.
  Comparable in shape to `capacity_ranking_robustness.md`'s 220 transport solves.

---

### #4 — `bbb_penetrance` as the breast permeability
**Shape: exact, per-compound AND per-tissue, and the discrepancy is already measured.**

* **The mismatch, measured.** `breast_tissue_permeability.md` §1: the measured breast
  5-FU:methotrexate permeability ratio is **1.02×**, against **28.1×** from
  `bbb_penetrance` — *"the measured data says the replacement should be
  drug-independent."* The engine therefore carries a **28× spurious per-drug spread** on
  the delivery side of the breast panel.
* **The repo already agrees.** `config.py:109-111`: *"`PS` is still driven by
  `bbb_penetrance`, a blood-BRAIN-barrier term that is simply wrong for non-brain
  tissue."* `scorecard.py` carries it as the `brain-permeability-in-breast` conditional.
* **Why it is still open.** The substitution was tried for **two** drugs against Kp anchors
  and made them worse (§2) — which is a real result about the *level*. It was never applied
  to the other eight compounds, and the **breast panel was never re-ranked** under a
  drug-independent permeability.
* **What breaks.** The breast top-1 (`vinorelbine`), the breast 3-D:2-D panel, and — since
  removing a 28× per-drug spread is a delivery-side per-compound change — plausibly the
  same class of flip that C3 records for GBM.
* **Cost.** Low: 10 solves on one cached section.

---

### #5 — The membrane-partition ceiling `_LOGP_KP_MAX = 3.0`
**Shape: a global clamp that demonstrably flattens a per-compound axis.**

* **The evidence is already in the module.** `permeability.py:163-174`: *"the ceiling
  binds… Cisplatin (22.2 µm/s), imatinib (21.0) and methotrexate (19.6) now sit within
  15 % of each other… Physically they differ by ORDERS of magnitude… no longer informative
  at the top."* And *"real Caco-2 Papp for these drugs is 0.001–1 µm/s, i.e. 1–4 orders
  below what this returns."*
* **Independent corroboration from a different subsystem.** `treatment/drugs.py:126-129`
  reports that after rescaling `PS_SCALE_S_INV` to its calibrated 0.16122 1/s, the residual
  is exactly **√8.70 = 2.950×** *"because the engine's cross-drug permeability ORDERING is
  wrong."* Two independent places in the codebase say the per-compound permeability order
  is wrong. Neither has been swept.
* **What breaks.** Per-compound delivery ordering on both panels — again the C3 class.
* **Cost.** Low-moderate: one recompile of the panel plus one solve per drug.

---

### #6 — The QUIESCENT band-width factor `2.0`
**Shape: not per-compound — per-tissue, and it hits the verdicts rather than the ranking.**

* **The constant.** `agents.py:83-84` writes `2.0 * t` inline. It is **not a config field**.
* **Why it is load-bearing.** `kill_calibration.md` §220: **90.0 % of tumour cells sit
  inside `[0.2, 0.4)`**, and `F_state` — the fraction those states produce — *"multiplies
  straight into the predicted ratio"*. `cell_state_impact.md` opens by attributing ~89 % of
  the TMZ kill ceiling to this band.
* **Why the existing sweep does not cover it.** `oxygen_calibration_breast.md` §4 sweeps
  the *lower* edge `t` over 2.5–14 mmHg, which moves **both** edges together. The band's
  *width ratio* has never been varied.
* **What breaks.** K1 and K2 — the 3-D:2-D verdicts. Because `F_state` is drug-independent
  it cannot move the ranking (that is C4/C5's territory), but the scored bands are **fixed**,
  so a common multiplier does not cancel out of a verdict the way it cancels out of an
  order. This is the clearest case in the document where "common-mode" is a defence for one
  claim and no defence at all for another.
* **Cost.** Low, but it requires touching `tmesim/` to expose the constant — which is
  itself the finding: an axis that cannot be swept without a code change is an axis nobody
  will sweep.

---

### #7 — `cell_partition` K = 1.0
Ranked below the above only because its status is **already honest and already
registered**: `transport_anchor_prereg.md` §3.8 declares it *"UNCONSTRAINED BY ANY
MEASUREMENT"* and TEST C NOT RUNNABLE at n = 4 spanning 54×. It is per-compound, it is a
load-bearing input to the transcellular route that has shipped by default since `fd2a9e7`,
and `capacity_ranking_robustness.md` §10 confirms that study does not constrain it either.
**Blocked on data, not on thought** — the only axis in the top seven of which that is true.

---

### #8 — Per-compound state-resistance multipliers
`PDConfig.hypoxia_resistance = 0.35` and `quiescent_resistance = 0.5` are uncited global
priors (`config.py:196-198`) documented as **drug-specific in reality**
(`config.py:210-235`), with exactly one per-drug override (TMZ, 0.80).
`CLAIM_INVENTORY.md` C5 correctly names this as C5's falsifier and calls it *"inert today
only while the hypoxic band is empty."*

**That defence has already expired on a second section and nobody noticed.**
`heldout_patient.md` §3.5 records that on section 2, *"34 % of tumour cells fall in the
model's HYPOXIC band versus 0 %"*. On that section the drug-dependent hypoxia multiplier is
live, and it has never been swept there. Cost: low.

---

### #9 — `CYTOPLASM_D_OVER_WATER = 0.27`
Measured for **one dye** (BCECF, Kao 1993) and used as the cytoplasmic diffusivity ratio
for every compound (`two_route.py:144`). Cytoplasmic hindrance is size-, charge- and
binding-dependent, i.e. per-compound. It feeds the transcellular route, which is now the
shipping default. Frozen in every scored arm (`transport_anchor_prereg.md` §2.3), never
varied. Ranked here rather than higher because its effect enters through the same
dimensionless group as `cell_radius_um` and may be partly absorbed by it — that
degeneracy is itself unmeasured.

---

### #10 — Panel composition
`capacity_ranking_robustness.md` §5.2–5.3 is the existence proof: **the same capacity
assignment reproduces the published order exactly on the 5-drug GBM sub-panel and changes
the top-1 pick on the 10-drug breast panel.** No document treats which compounds are on a
panel as an axis, yet it is demonstrably decisive. Cost: free — it is a re-analysis of
sub-panels of results that already exist.

---

### Named but ranked low, with reasons
* **Ionisation / extracellular pH** — genuinely absent (`permeability.py:16-39`), and
  per-compound via pKa. Ranked low **only** because adding it is a modelling exercise, not
  a sweep: *"a pH-partitioning term would have to be designed rather than parameterised."*
  Its absence is a model-form gap, which §12d of the methodology says the whole approach
  does not address.
* **Interstitial viscosity, temperature, the MW regime threshold** — common-mode scales
  that the ranking is structurally insulated against, and whose action does not obviously
  factor through anything per-compound. Low expected yield; listed for completeness and
  because the MW threshold carries a live 5,000-vs-10,000 inconsistency worth fixing on
  correctness grounds alone.
* **The `1e-300` floors and clamps** — no reason to expect a claim depends on them, except
  that `subgrid.py:100-103` records `beta` reaching 144.0 against a validated max of 5.677
  — with 35.77 % of partial cells beyond that maximum — which is a clamp doing real work in
  a module that ships off.

---

## 8. The judgement

**Is the published limit list evidence of thoroughness, or of the axes someone happened to
pick? Both — and the ratio is legible.**

*Evidence of genuine thoroughness, and it is not small.* The numerics are probed past the
point of comfort: eight-rung grid ladders on two panels, an off-lattice conservation
check, a boundary-layer error law fitted to 54 points, a CFL re-measurement that
**retracted a claim in the codebase's own comments**. The oxygen route was swept until it
was shown to span its entire answer. Four pharmacology mechanisms were implemented,
swept, and left default-off because they did not help — including one
(`intracellular_retardation`) that was refuted after being built. Verdicts are
pre-declared, `codon validate` exits 1 on FAIL, and `CLAIM_INVENTORY.md` §5 keeps a
standing list of the project's own corrected beliefs. Very few models are attacked this
hard.

*Evidence of selection, and this is the load-bearing half.* Sort the swept axes by whether
their action factors through a per-compound quantity and the list separates cleanly.
**Seven axes with per-compound action have been swept in the whole corpus; six of them were
scored against an ordering; all six moved it. The seventh (per-drug logP) was swept and the
ordering readout was never taken.** The remaining ~45 swept-axis entries in §3 are
common-mode, numerical or reconstruction-side, and they held — which is what the published
robustness claims C1, C2 and C4 are made of. The uncertainty quantification samples four
knobs, one of which is *provably* a gauge coordinate and three of which measure at or
below the noise floor; meanwhile the nine per-drug ratios that α cannot stand in for are
unconstrained **and already known to be wrong in direction**. That is not a model that was
swept in the directions most likely to break it. It is a model that was swept thoroughly
in the directions that were easiest to parameterise — grids, thresholds, bandwidths,
global scales — and thinly in the direction its central claim actually lives in.

**The ratio, stated plainly: the numerical and reconstruction layers are close to
exhaustively probed; the per-compound pharmacology layer, which is what a compound ranking
is made of, has had roughly six of its axes touched and has failed on every one.**

### What writing this changed about the strongest defensible claim

Before: *"the relative comparison of compounds on one section"* (`CLAIM_INVENTORY.md` §0).

I no longer think that is the strongest defensible claim, and I think it is now slightly
too strong. `capacity_ranking_robustness.md` has already retired the unqualified version,
but the reason it gives is specific to capacity. The inventory says something more general
and less comfortable:

> **The comparative claim is the one claim that sits *entirely* in the layer this project
> has probed least.** Every published invariance supporting it was measured on a
> common-mode axis. Every per-compound axis ever swept has moved it. The comparison is not
> refuted — but "attacked from several directions and held" is not accurate as stated,
> because the directions were mostly the ones it is structurally insulated against.

The claim I would now put forward as strongest is the one
`capacity_ranking_robustness.md` §9 arrived at from a single axis and that this inventory
generalises to all of them:

> **Codon's steady-state and equilibrium-referenced outputs — penetration length,
> interstitial:plasma ratio, converged gradient shape — are exactly invariant to a whole
> class of per-compound parameter error, and that invariance is proven rather than
> measured.** The finite-window kill ranking is not, and the axes on which it has been
> tested are known and enumerable (§3), as are the ones on which it has not (§5, §7).

That is narrower. It is also the only version that survives contact with this file, and it
has one property the current claim lacks: it names the *class* of thing that cannot hurt
it, which means a reader can check whether their axis is in the class.

**And the corollary that matters most for the next surface Codon publishes:** the second
strongest asset here is not any result. It is that this project can produce §5.3 — a list
of things nobody thought of — and print it. A published limit list without an
accompanying *unswept* list invites exactly the misreading that six named failures means
six axes were the ones that mattered. Shipping both is a stronger claim than shipping a
longer list of failures.

---

## 9. What would most change this document

1. **Verify the linearity behind #1.** If AUC is exactly linear in plasma level under the
   shipping law, a per-compound PK re-rank costs no solves and either promotes #1 to a
   retraction or removes it from the list. This is the cheapest item here by a wide margin.
2. **Establish the real per-compound spread of Hill slopes on these two panels.** #2's rank
   rests on the spread being large, which this file does not establish.
3. **Re-rank the breast panel under a drug-independent permeability.** #4 is the only entry
   whose motivating discrepancy is already measured (28.1× vs 1.02×).
4. **Expose the `2.0` band-width factor as a config field.** Not to sweep it — to make it
   sweepable. An axis that requires a code edit to vary will not be varied.
5. **Any correction to §5.3.** The "never considered" list is the most valuable and the
   most fragile part of this document. If an axis there is in fact registered somewhere,
   moving it to §5.1 or §5.2 makes this file more useful, not less.

---

*Compiled 2026-08-21 by reading `results/` (132 markdown documents), `docs/VERIFICATION_
METHODOLOGY.md`, `docs/MATHEMATICAL_MODEL.md`, `tmesim/validation/scorecard.py` and the
engine modules under `tmesim/`. No sweep was run and no default was moved. Test suite not
executed by this pass and not expected to move: 1171 passed / 0 failed.*
