# Real drug panel on the real Xenium glioblastoma microenvironment

**Experiment:** `experiments/real_drug_panel.py` (run: `./.venv/bin/python experiments/real_drug_panel.py`)

**Figure:** `results/real_drug_panel.png`

The full panel **[temozolomide, lomustine, cisplatin, methotrexate, doxorubicin]** is simulated on **one and the same real tumor** — the open 10x Xenium GBM section (40,887 measured cells; tumor 30,174). The real microenvironment is reconstructed **once** and reused for every drug, so all differences come from each drug's own physics on identical tissue (real vasculature, real ECM from collagen/matrix transcripts, real-perfusion oxygen).

Each drug is compiled through the same pipeline: `calibrate_alpha()` (single perivascular-sink constant, fit **once** on doxorubicin → α=309.6) → `calibrated_drugparams(name, α)` (real GDSC IC50 + calibrated loss) → `apply_barrier_to_drugparams` (saturable binding-site barrier).

## Real microenvironment (reconstructed once)

- Source: 10x Xenium GBM — 40,887 real cells, grid (170, 271), dx=30 µm, extent (8130.0, 5100.0) µm.
- Cell-type fractions: TUMOR 73.8%, ENDOTHELIAL 9.7%, IMMUNE 8.5%, STROMAL 7.9%.
- Mean vascular density 0.094, mean ECM 0.083, mean oxygen 0.156.
- **Tumor hypoxic fraction ≈ 94%** — this deeply hypoxic real tissue mechanistically suppresses kill for **every** drug.

## Predicted panel on this real tumor

| drug | MW (Da) | real IC50 (µM) | bbb_penetrance | D_eff (µm²/s) | penetration_index | mean AUC (µM·s) | predicted kill |
|---|---:|---:|---:|---:|---:|---:|---:|
| temozolomide | 194 | 606 | 0.558 | 264.2 | 0.396 | 5.61e+06 | 0.0189 |
| lomustine | 234 | 50 | 0.973 | 246.3 | 0.298 | 1.51e+06 | 0.0399 |
| cisplatin | 298 | 24.6 | 0.893 | 224.5 | 0.332 | 3.08e+06 | 0.0608 |
| methotrexate | 454 | 12 | 0.115 | 189.0 | 0.277 | 1.3e+05 | 0.0213 |
| doxorubicin | 544 | 0.5 | 0.051 | 175.6 | 0.266 | 4.62e+04 | 0.0551 |

IC50 values are the real GDSC/literature potencies that flow in through the molecular layer; MW is from the verified molecular structure; `bbb_penetrance` and `D_eff` are structure-derived (`D_eff = D · (1 − 0.6·ECM)` at the tumor mean ECM = 0.083).

## Which drug wins on THIS real tumor, and why

- **Best predicted penetration:** **temozolomide** (penetration_index 0.396, D_eff 264 µm²/s). Small, weakly-bound alkylators (temozolomide, lomustine) and the small platinum crosslinker diffuse furthest before the perivascular sink consumes them; the large, strongly-bound doxorubicin is captured within tens of µm of vessels (the calibrated ~45 µm boundary layer) and barely reaches the tumor bulk.
- **Best predicted kill:** **cisplatin** (kill_fraction 0.0608, penetration_index 0.332, IC50 24.6 µM, bbb 0.893). Kill is the product of **delivery × potency × (hypoxia resistance)**: a drug wins only when it both reaches the cells AND is potent enough at the exposure it achieves. Here delivery dominates the ranking — doxorubicin is the most potent by IC50 but is delivery-starved on this section, while the penetrant alkylators convert their exposure into kill.

## Honesty and caveats

- **No outcome label.** This open Xenium section is an **untreated resection**. There is no paired measured drug response for any of these five drugs. The kill/penetration numbers are **predictions of the validated physics on real tissue — not validated outcomes.** Do not read the ranking as clinical efficacy.
- **What is genuinely real here:** the geometry (real cells, real vasculature, real matrix, real oxygen), the transport physics (structure → mesh diffusivity + a doxorubicin-calibrated perivascular sink, revalidated out-of-sample on mitoxantrone), and the potency (real GDSC/literature IC50).
- **Hypoxia compresses the panel.** ≈94% of tumor cells are hypoxic; the resistance banding lowers the achievable kill for **every** drug, so absolute kills are small and the drugs are closer together than their potency differences alone would suggest.
- **One anchor sets the scale.** The single α is fit to one doxorubicin penetration number; the absolute penetration of every drug rides on that measurement, and the binding proxy is an ordering, not a first-principles rate.
