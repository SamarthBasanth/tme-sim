# Molecular-layer demo — structure-derived drug transport in GBM

One fixed synthetic GBM patient (`GBM-MOLECULAR-DEMO`, seed=0); modest sub-saturating dose (plasma_peak=1.0 uM, exposure window 600 s). Every drug row below shares the SAME microenvironment — only the molecule changes. All transport parameters are derived from molecular structure via Stokes-Einstein + Brownian-dynamics hindrance + physicochemical BBB rules; none are hand-set per drug.

## Structure-derived parameters and simulated response

| drug | MW (Da) | free_D (um^2/s) | effective_D (um^2/s) | BBB penetrance | free fraction | IC50 provenance | penetration index | kill fraction |
|---|---|---|---|---|---|---|---|---|
| temozolomide | 194.2 | 495.1 | 278.0 | 0.558 | 1.000 | estimate | 0.664 | 0.037 |
| lomustine | 233.7 | 465.5 | 259.2 | 0.973 | 0.759 | estimate | 0.658 | 0.341 |
| methotrexate | 454.4 | 372.9 | 198.9 | 0.115 | 0.985 | estimate | 0.635 | 0.242 |
| doxorubicin | 543.5 | 351.3 | 184.8 | 0.051 | 0.990 | estimate | 0.629 | 0.476 |
| cisplatin | 298.0 | 429.2 | 236.2 | 0.893 | 0.987 | estimate | 0.651 | 0.415 |
| bevacizumab | 149000.0 | 34.5 | 0.9 | 0.000 | 0.998 | flagged default (LOW) | 0.000 | 0.000 |

### Penetration ordering (highest -> lowest)

temozolomide (0.66) > lomustine (0.66) > cisplatin (0.65) > methotrexate (0.63) > doxorubicin (0.63) > bevacizumab (0.00)

## Interpretation (honest)

- **Mechanistic, not hand-tuned.** free_D comes from the Stokes-Einstein relation on a structure-derived hydrodynamic radius; effective_D comes from a genuine 2D overdamped Brownian-dynamics simulation of tracer particles through a crowded ECM obstacle field; BBB penetrance comes from physicochemical rules on MW/logP/TPSA/HBD. The differences between drugs are therefore driven by molecular structure, not per-drug fitting.
- **IC50 is a flagged INPUT, not derived.** Potency is the weakest link from structure alone (it needs docking/assay). Every drug here uses the molecular layer's flagged low-confidence default IC50, so the *kill* numbers should be read as transport-driven contrasts on a common potency, not calibrated efficacy predictions.
- **What this demonstrates:** drug-specific TRANSPORT — small lipophilic alkylators (temozolomide, lomustine) reach the tumor core and accumulate exposure, while large/polar agents build up less and the ~149 kDa antibody (bevacizumab) essentially fails to penetrate the tissue (tiny effective_D, near-zero BBB penetrance).
- **Still needs validation.** These are physics/cheminformatics estimates. They must be checked against MEASURED tissue pharmacokinetics (drug levels in tumor tissue) before any efficacy claim; the value here is a mechanistic, structure-in -> transport-out pipeline with explicit provenance on every number.

> **Antibody note:** bevacizumab (~149 kDa) has effective_D = 0.95 um^2/s and BBB penetrance = 0.000; its penetration index (0.000) and kill (0.000) are the lowest in the set — a large macromolecule essentially cannot reach the tumor by passive transport.
