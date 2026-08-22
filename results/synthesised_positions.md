# Does the delivery map survive a section whose cell positions are SYNTHESISED?

**A measurement on two shipping sections. No model change.** Nothing under
`tmesim/` was modified to produce it, both reconstructions are the production
entry point (`reconstruct.real_pipeline.reconstruct_real`) at the same grid, and
the distance field measured is the one the product actually bins by
(`validation.vessel_profile.direct_vessel_distance_um` over
`endothelial_xy(micro.cells)`).

* date: 2026-08-21
* script: `scratchpad/lattice_test.py` (read-only against the engine; not
  committed — the numbers below are the record)
* the disclosure this supports: `app/server/api.py::SPOT_SPATIAL_REFUSAL` and
  the matching capability caveat in `tmesim/data/ingest.py`
* the product surface it governs: the DELIVERY EVIDENCE half of
  `app/web/js/ui/reach.js`

---

## 0. The answer

> ## NO — and the obvious way to check would have said yes.
>
> On a Visium section the pipeline synthesises 12 pseudo-cells per spot at the
> spot centroid with uniform jitter. Every quantity in the delivery half of the
> product is downstream of those positions. **The distance field on such a
> section is a picture of the spot lattice's occupancy, not of the tissue's
> vasculature**, and the discrepancy is **worst exactly where the delivery claim
> lives** — a factor of **5.7x at p99**, not a uniform offset.
>
> **THE OBVIOUS TEST IS A FALSE NEGATIVE — say this first.** The natural check
> is to look for the spot lattice's ~100 um period in the coordinates. It is not
> there. The jitter fully de-quantises them: **99.93 %** of x values are
> distinct on the synthesised section, against **99.91 %** on the measured one.
> A reviewer running that check finds nothing and clears the section.
>
> The check that works is not on the GEOMETRY OF THE COORDINATES. It is on the
> STATISTICS OF THE SPACING.

---

## 1. What was measured

Both sections reconstructed at `dx_um = 30.0`, defaults otherwise
(`tissue_mask_mode='geometric'`, `oxygen_mode='pde'`, `state_mode='oxygen'`).

| | `xenium_319` (positions MEASURED) | `visium_10x_parent_gbm` (positions SYNTHESISED) | ratio |
|---|---:|---:|---:|
| endothelial cells | 4,295 | 4,029 | 0.94x |
| **coordinate quantisation** | | | |
| distinct x / n cells | 0.9991 | 0.9993 | — |
| distinct y / n cells | 0.9995 | 0.9998 | — |
| modal x gap | 1.200 um | 0.684 um | — |
| **nearest-neighbour spacing between endothelial cells** | | | |
| p10 | 11.46 um | **4.90 um** | **0.43x** |
| p50 | 22.95 um | 16.31 um | 0.71x |
| p90 | 72.02 um | 64.13 um | 0.89x |
| coefficient of variation | 0.756 | **0.999** | 1.32x |
| **distance from tissue to nearest endothelial cell** | | | |
| p50 | 44.1 um | 69.7 um | 1.58x |
| p90 | 85.1 um | **273.7 um** | **3.22x** |
| p99 | 129.9 um | **736.6 um** | **5.67x** |
| max | 795.6 um | 1155.7 um | 1.45x |
| tissue fraction beyond 150 um | **0.001** | **0.023** | **23x** |

---

## 2. The mechanism, and why it is not about the lattice

The synthesis fails for a reason that is more general than the grid it happens
to sit on:

> **Real endothelium is a CONNECTED structure.** Vessels are continuous, so the
> distance to the nearest one stays short almost everywhere and the histogram
> decays away by ~150 um.

Twelve pseudo-cells jittered inside a single ~55 um spot disc are a **clump**:
they pack more tightly than real endothelium (p10 nearest-neighbour spacing is
**less than half**, 4.90 um against 11.46 um). And wherever a spot is not called
endothelial, there is a **void the size of the lattice cell**. Clumps plus voids
give a spacing distribution with a coefficient of variation of **0.999** —
indistinguishable from a Poisson process — where measured tissue gives **0.756**,
i.e. materially more regular than random, which is what a connected network
looks like when sampled.

That is why the failure appears in the TAIL rather than in the period. The
positions are not on a grid after jitter; the *connectivity* is what was never
there.

### 2.1 Why the tail is the fatal part

A uniform offset in the distance field would be survivable — the delivery half
of the product reports RATIOS inside one map, and a constant rescaling of the
distance axis would move the falloff curve without inverting any reading. What
is measured is not a rescaling:

* p50 differs by **1.58x**
* p90 by **3.22x**
* p99 by **5.67x**

The distortion grows monotonically into the tail, and **the tail is exactly the
under-reached region the delivery half exists to describe**. The fraction of
tissue more than 150 um from a vessel — the population a customer would read as
"the part the drug cannot reach" — is **23 times larger** on the synthesised
section, and every bit of that excess is an artefact of spot occupancy.

---

## 3. METHOD NOTE — the test that fails, and the test that works

Recorded separately because it is reusable and because the failing test is the
one a reasonable person reaches for first.

**DOES NOT WORK — looking for the lattice in the coordinates.** Periodicity,
modal gaps, hard cutoffs at the spot pitch, quantisation of x or y. Uniform
jitter inside the spot disc destroys all of it. On the numbers above, a
quantisation check would rank the SYNTHESISED section as marginally *less*
quantised than the measured one (0.9993 vs 0.9991 distinct x).

**WORKS — the statistics of the spacing.** Three quantities, in decreasing order
of how cleanly they separated the two sections here:

1. **fraction of tissue beyond a fixed distance** (150 um): 23x separation;
2. **vessel-distance p99**: 5.7x separation;
3. **nearest-neighbour spacing p10 and its coefficient of variation**: 2.3x and
   1.32x — the clumping signature, and the only one of the three that reads the
   endothelial cells directly rather than the field derived from them.

---

## 4. Scope of what this section is, and is not, evidence for

* **N = 2.** One measured section and one synthesised section, one platform
  each. The separation is large on every statistic, but this is not a
  distribution over sections and no threshold is proposed from it.
* **The chain from positions to the product's numbers is read from code, not
  measured here.** `reach.py` reuses `penetration_length_fit(...)["drug"]
  ["perivascular_mean"]` as its reference and bins by the same distance field,
  so the inheritance is structural rather than empirical. What this document
  measures is the INPUT to that chain.
* **It says nothing about whether the Visium reconstruction is useful for
  anything else.** The potency half of the product — the per-cell mechanism
  verdict, the coverage ladder, the exposure shortfall — does not read where
  cells sit and is untouched by this finding.

---

## 5. A LIBRARY FINDING, recorded because it bounds what can be sold

Of the **38** sections the Studio lists as available, **4** carry measured
per-cell positions:

| section | platform | tissue |
|---|---|---|
| `xenium_319` | Xenium | GBM |
| `xenium_io` | Xenium | GBM |
| `gse324268` | Xenium | GBM |
| `xenium_breast` | Xenium | breast |

The other **34** are Visium or GBmap spot sections. Every one of them carries a
blocking `spatial_spread = NO` capability caveat, and after this measurement
that caveat also covers the delivery map.

Consequences, stated plainly:

1. **The delivery half of the product is meaningful on four sections**, three of
   which are the same tumour type.
2. **None of the four is known to be poorly perfused.** On `xenium_319` the
   measured falloff decays only 1.00 -> 0.50 over 145 um at dx = 12 um, and the
   per-cell attribution returns 99.8 % DOSE-limited. A non-degenerate
   DELIVERY-limited verdict has not been demonstrated on any available section:
   the only DELIVERY verdict obtained (bevacizumab) comes from
   `bbb_penetrance = 0`, which is an input rather than a discovery.
3. **The single most valuable acquisition for this product is therefore a
   measured-position section from a poorly perfused tumour** — not another
   compound and not another anchor. It is the only way to demonstrate the
   DELIVERY verdict non-degenerately, and the only way to show the map saying
   something a customer's own imaging does not already tell them.

---

## 6. PROPOSED CAPABILITY GATE — not built, and the reason it is worth building

**This is a proposal, not a result. Nothing below has been implemented or
validated.**

Bring-your-own-tissue currently checks format, layout, gene panel and cell
count. **It does not check whether the endothelial field looks like
vasculature.** A customer can register a correctly-formatted section whose
endothelial calls are sparse, mislabelled, or themselves spot-derived, and
receive a delivery map that is a picture of their annotation rather than of
their vessels — with nothing in the product to catch it.

The statistics in sec. 3 are exactly that check, and they are cheap: they need
the endothelial positions and the tissue mask, both of which ingest already has,
and no simulation at all.

Reference distribution, from the one measured section characterised here
(`xenium_319`, dx = 30 um) — **a single section, so this is a starting point for
a gate and not a specification**:

| statistic | measured (`xenium_319`) | synthesised (`visium_10x_parent_gbm`) |
|---|---:|---:|
| NN spacing p10 | 11.46 um | 4.90 um |
| NN spacing CV | 0.756 | 0.999 |
| vessel-distance p90 | 85.1 um | 273.7 um |
| vessel-distance p99 | 129.9 um | 736.6 um |
| tissue fraction beyond 150 um | 0.001 | 0.023 |

The gate would refuse — or, better, DOWNGRADE the delivery half to unavailable
while leaving the potency half intact, which is the split
`app/web/js/ui/reach.js` already implements — when a section's endothelial field
looks clumped-and-voided rather than connected.

**Before it could ship it needs what this document does not have:** the four
measured sections characterised rather than one, so the reference is a range and
not a point; and a decision, made in advance, about what a section that falls
between the two columns above should be told. A gate calibrated on N = 1 that
refuses a customer's real tissue would be a worse failure than the one it
prevents.

---

## 7. What would falsify this

* A spot section whose vessel-distance p99 and beyond-150 um tissue fraction sit
  inside the measured section's range. That would mean the synthesis can
  reproduce a connected vascular field, and sec. 2's mechanism is wrong.
* A measured-position section with a CV near 1.0 and a long distance tail — i.e.
  genuinely clumped, genuinely poorly perfused real tissue. That would mean the
  proposed gate cannot distinguish "synthesised" from "poorly perfused", which
  is the exact case the product most wants to detect, and the gate as described
  in sec. 6 would have to be abandoned rather than tuned.

The second is the more likely of the two, and it is the reason sec. 6 is a
proposal.
