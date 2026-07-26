# Reality-Check Summary (v2 additions)

Four upgrades turned the tissue simulator into a **multiscale mechanistic engine**
with real potency and a first validation against measured reality.

## 1. Molecular layer — structure → tissue parameters
Each drug's transport parameters are now computed from its actual molecular
structure, not hand-set: RDKit descriptors → Stokes-Einstein free diffusivity →
Brownian-dynamics ECM hindrance → physicochemical BBB/permeability/binding →
`DrugParams`. Six real GBM drugs through the SAME tumor differ mechanistically;
the 149 kDa antibody cannot penetrate. (`results/molecular_*.png`)

## 2. Real potency (measured IC50)
IC50 now resolves **measured GDSC glioma values** where available, else cited
literature, else a flagged default:
- temozolomide **606 µM** (51 glioma cell lines — correctly weak in vitro)
- cisplatin **24.5 µM** (34 lines), methotrexate **12 µM** (49 lines)
- doxorubicin 0.5 µM (literature; not in the glioma screen)
This flips the demo into a realistic tension: doxorubicin is potent but BBB-blocked;
temozolomide reaches tissue but is weak in vitro.

## 3. nm-scale ECM mesh hindrance
Added Ogston-Renkin fiber-matrix theory + nm-scale Brownian dynamics so drug SIZE
gates diffusion (~22× small-molecule vs antibody). Wired into the compile bridge:
the antibody's effective diffusivity now collapses (~1 µm²/s).

## 4. Validation against measured reality — the headline
Checked the sim's predicted doxorubicin tissue penetration against the published
anchor: doxorubicin penetrates only **~40-50 µm** from tumor vessels (Primeau et
al., Clin Cancer Res 2005;11(24):8782-8788, PMID 16361566, confirmed via PubMed).

| metric | sim | published | ratio | band (±3×) |
|---|---|---|---|---|
| analytic penetration length L | ~700 µm | ~45 µm | 15.5× | **FAIL** |
| single-vessel 1/e distance | ~190 µm | ~45 µm | 4.2× | **FAIL** |
| ordering: TMZ reaches further than doxorubicin | yes | yes | — | PASS |
| ordering: small-molecule D ≫ antibody D | ~13× | yes | — | PASS |

**Verdict: directional physics correct, absolute magnitude wrong by ~an order of
magnitude.** The sim over-predicts penetration because the structure-only
binding/clearance priors (~1 h timescale) are ~300× too slow — real perivascular
gradients are driven by rapid high-capacity DNA/tissue binding and cellular uptake.
Diagnosed, not fudged: this is a parameter-calibration gap, not a physics/numerics
error.

## What this means
This is the project validating its own thesis. The engine looks plausible, gets
every *ordering* right, and is numerically sound — yet against **one real measured
number** it is 15× off, and the harness caught it and named the exact parameter to
fix. A trustworthy in-silico microenvironment is not the one that renders nicely;
it's the one that has been checked against reality and told you where it's wrong.
Next calibration target: fit tissue binding/clearance to reproduce the ~45 µm
doxorubicin gradient, then re-validate on a second drug.
