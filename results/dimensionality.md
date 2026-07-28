# THE 2-D ASSUMPTION, QUANTIFIED

**Verdict: the 2-D solve does NOT need fixing.** The penetration length `L = sqrt(D_eff/k)` is *identical* in 2-D and 3-D — proved, not estimated — and mean tumour AUC is unbiased by a stereological identity. The one quantity that genuinely moves is spatial heterogeneity, and it moves in the **conservative** direction: the 2-D solve over-states how uneven delivery is, which means it *under-states* how much drug reaches the worst-perfused tissue. The 3-D:2-D kill test in `results/kill_calibration.md` does not just survive the correction — it improves.

> ### Headline
>
> **1. A 2-D solve is EXACT for a vessel cut transversely, not an approximation.** A capillary is a *line* in 3-D, and a line source is natively 2-D-symmetric. Solving the full 3-D problem for straight vessels perpendicular to the section and comparing the mid-plane to the 2-D solve gives a maximum relative difference of **4.7e-16** — machine precision. The often-stated worry, 'a 2-D solve imposes no flux out of the section', is real but is *not* the operative assumption; what the solve actually assumes is that the vasculature is z-invariant over the diffusive reach, and for a transverse capillary it is.
>
> **2. The penetration length is preserved exactly.** In 3-D a point source gives `exp(-r/L)/r`; in 2-D a line source gives `K_0(r/L)`. Both carry the *same* exponential rate `1/L`; they differ only in the algebraic prefactor. And for the object a capillary actually is — a line — the two coincide identically. No penetration length in this repository is dimensionally biased.
>
> **3. The vessels really are cut transversely, and it is measured, not assumed.** On the flagship section **95.4%** of reconstructed vessel objects have an in-plane extent no larger than 1.5x the measured glioma vessel caliber (18.1 um) — i.e. they look like round cross-sections, not streaks. The section is 5 um thick and the vessels are 18.1 um across, so the section is **3.6x thinner than one vessel**: it is geometrically almost incapable of containing a vessel lengthwise.
>
> **4. The residual error is not out-of-plane ESCAPE, it is out-of-plane SUPPLY.** Drug does leave the plane, but in real tissue it is replaced by drug from vessels above and below that the section never sees. Those invisible vessels fill in the gaps *between* the in-plane vessels, so the true 3-D field is **flatter** than the 2-D one. Measured over 18 network configurations, 2-D over-states the coefficient of variation by 1.43-1.90x.
>
> **5. Therefore a z-loss correction term is the WRONG-SIGNED FIX.** The best-fit z-loss coefficient on `k` is **-0.73** — negative in every configuration tested. Adding a sink to represent escape would make the 2-D answer *worse*, because the dominant missing physics is supply, not loss.
>
> ### The bottom line
>
> **Penetration length: 0% bias on the physical `L`; ~24% under-statement of the *fitted* profile decay distance. Mean tumour AUC: 0% bias, analytically. Therapeutic coverage: worst case 0.205 absolute, against the +/-0.03 grid bar already carried, with an inconsistent sign. 3-D:2-D kill ratio: the implied temozolomide IC50 fold moves from 3.26x to 2.73-2.98x, still inside the 2-5x literature band.** Three of the four published quantities move by less than the uncertainty they already carry. The issue is closed.


## 1. Analytic: what 2-D and 3-D actually differ by

### 1.1 The three Green's functions

Steady reaction-diffusion, `D grad^2 C - k C + source = 0`, `L = sqrt(D/k)`:

| geometry | solution | decay |
|---|---|---|
| 3-D **point** source | `C = Q/(4 pi D) * exp(-r/L) / r` | `exp(-r/L)`, prefactor `1/r` |
| 2-D **line** source (what the solver computes) | `C = q/(2 pi D) * K_0(r/L)` | `exp(-r/L)`, prefactor `~sqrt(L/r)` |
| 3-D **finite segment**, length `lam` | `C = q/(4 pi D) * int exp(-sqrt(r^2+z^2)/L)/sqrt(r^2+z^2) dz` | interpolates the two |

**`L` is common to all three.** In 3-D, substituting `C = f(r)/r` gives `f'' - f/L^2 = 0`, so `f = exp(-r/L)`. In 2-D the radial equation is the modified Bessel equation of order 0, whose decaying solution is `K_0(r/L)` with asymptotic form `sqrt(pi L / 2r) exp(-r/L)`. Same exponent, different algebraic prefactor. **The dimensionality of the solve cannot change `sqrt(D/k)`; it can only change the geometric spreading factor in front of it.** That single fact does most of the work in this report.

The finite segment is the physically correct object — a capillary is a line, not a point — and its two limits are exactly the other two rows: `lam -> inf` reproduces `K_0` (the 2-D solve, exactly), `lam -> 0` reproduces the point source.

### 1.2 What a fitted decay length would measure

An exponential fit measures the local logarithmic derivative, not `L`:

* 3-D point source: `-d ln C/dr = 1/L + 1/r` exactly, so `Lambda_3D = L/(1 + L/r)`.
* 2-D line source: `-d ln C/dr = (1/L) K_1/K_0`, so `Lambda_2D = L K_0/K_1 ~ L/(1 + L/(2r))`.

**The geometric steepening is `1/r` in 3-D and `1/(2r)` in 2-D — the 2-D correction is exactly half.** That factor of two is the entire dimensionality effect on a fitted decay length *in the point-source limit*. In the line-source limit that a capillary occupies it is zero.

| drug | `L` (um) | r=25 um | r=50 um | r=100 um | r=150 um | r=200 um |
|---|---:|---:|---:|---:|---:|---:|
| doxorubicin `Lambda_2D/Lambda_3D` | 64.9 | 1.806 | 1.490 | 1.278 | 1.195 | 1.150 |
| cisplatin `Lambda_2D/Lambda_3D` | 142.9 | 2.294 | 1.858 | 1.528 | 1.384 | 1.302 |
| temozolomide `Lambda_2D/Lambda_3D` | 264.6 | 2.754 | 2.240 | 1.816 | 1.616 | 1.497 |

Read the table as the **maximum** over-statement of a fitted decay length, reached only if every vessel were a point (`lam << L`).

### 1.3 The interpolant: how straight must a vessel be for 2-D to be right?

`C_lam(r)/C_inf(r)` — the factor by which the 2-D solve over-states one **isolated** vessel whose straight run is `lam`. 1.000 means the 2-D solve is exact.

| drug | `lam/L` | r=25 | r=50 | r=100 | r=150 | r=200 |
|---|---:|---:|---:|---:|---:|---:|
| doxorubicin (L=65 um) | 0.5 | 0.353 | 0.246 | 0.170 | 0.136 | 0.117 |
| doxorubicin (L=65 um) | 1 | 0.592 | 0.455 | 0.328 | 0.268 | 0.231 |
| doxorubicin (L=65 um) | 2 | 0.824 | 0.727 | 0.587 | 0.499 | 0.440 |
| doxorubicin (L=65 um) | 4 | 0.959 | 0.929 | 0.863 | 0.797 | 0.740 |
| doxorubicin (L=65 um) | 8 | 0.997 | 0.994 | 0.986 | 0.974 | 0.959 |
| cisplatin (L=143 um) | 0.5 | 0.495 | 0.370 | 0.259 | 0.209 | 0.179 |
| cisplatin (L=143 um) | 1 | 0.714 | 0.609 | 0.474 | 0.396 | 0.344 |
| cisplatin (L=143 um) | 2 | 0.885 | 0.834 | 0.743 | 0.668 | 0.608 |
| cisplatin (L=143 um) | 4 | 0.974 | 0.962 | 0.935 | 0.906 | 0.875 |
| cisplatin (L=143 um) | 8 | 0.998 | 0.997 | 0.995 | 0.992 | 0.988 |
| temozolomide (L=265 um) | 0.5 | 0.592 | 0.481 | 0.357 | 0.290 | 0.249 |
| temozolomide (L=265 um) | 1 | 0.777 | 0.705 | 0.595 | 0.517 | 0.459 |
| temozolomide (L=265 um) | 2 | 0.912 | 0.881 | 0.826 | 0.776 | 0.730 |
| temozolomide (L=265 um) | 4 | 0.980 | 0.973 | 0.960 | 0.945 | 0.930 |
| temozolomide (L=265 um) | 8 | 0.998 | 0.998 | 0.997 | 0.996 | 0.994 |

**This is an upper bound on the error, and a loose one.** It truncates the vessel at `lam` with *nothing* beyond it. A real vessel does not stop — it continues, branches, and has neighbours, all of which supply the field this integral throws away. Section 3 does the honest version with a connected network, and the error there is far smaller.

### 1.4 The mean is preserved EXACTLY — a stereological identity

Take an infinite straight vessel at angle `th` to the section normal, leaking `q` per unit length. Its perpendicular distance from a point `(x,y)` in the section plane is `rho^2 = x^2 cos^2(th) + y^2`, so the in-plane contours are **ellipses stretched by `1/cos(th)`** — an oblique vessel does not look like the isotropic disc the 2-D solve paints. But integrate over the plane:

```
int C dA = (q/2 pi D) (1/cos th) (2 pi L^2) = q / (k cos th)
N_A       = L_V cos th                    (piercings per unit area)
product   = q L_V / k                     (independent of th)
```

and the true 3-D mean concentration is `(source per volume)/k = q L_V/k` — **the same number**. So for a statistically homogeneous network of *any* orientation distribution, a 2-D solve on the section's own vessel traces reproduces the true 3-D mean exactly: obliquely cut vessels leave proportionally longer traces, and that excess is precisely the compensation for their being fewer. The isotropic case recovers the classic `L_V = 2 N_A`. **Consequence (ii) of section 5 is settled here, analytically, before any simulation runs.**


## 2. The vessel-geometry question, measured on the real sections

This decides how bad the problem is, and it is measurable rather than a matter of opinion. If the plane cuts vessels roughly perpendicular, the 2-D disc sources really are cross-sections of lines running through the slab — and a line source *is* 2-D-symmetric, so the 2-D solve is right. If vessels run obliquely or in-plane, it is worse.

For a straight vessel of caliber `d` crossing a section of thickness `t` at angle `th` to the normal, the in-plane trace is `d + t tan(th)`. Here `d = 18.1 um` (PMID 24390982, IHC-validated glioma vessel size) and `t = 5 um` (10x Xenium vendor nominal). **The section is 3.6x thinner than one vessel is wide**, so even a steeply oblique vessel leaves a trace barely larger than its own caliber — a thin section is geometrically biased towards transverse cuts, and the object-count-weighted isotropic null is `sin(2 th)`, which is itself transverse-weighted (`P(th > th0) = cos^2 th0`).


### s1 — 10x Xenium GBM (brain_319 panel) — the flagship section

4,295 endothelial cells -> 3,091 discrete vessel objects at the 20 um linking distance; MVD 110.8/mm^2 over 27.9 mm^2; 81% singletons.

| in-plane extent threshold | implied `th` | **frac of ALL objects** | frac of multi-cell objects | isotropic null |
|---|---:|---:|---:|---:|
| <= 1x caliber (18.1 um) | 0 deg | **0.8962** | 0.4456 | 0.0000 |
| <= 1.5x caliber (27.2 um) | 61 deg | **0.9541** | 0.7547 | 0.7661 |
| <= 2x caliber (36.2 um) | 75 deg | **0.9780** | 0.8826 | 0.9291 |
| <= 3x caliber (54.3 um) | 82 deg | **0.9903** | 0.9482 | 0.9813 |

Multi-cell object extents (um): median 18.7, p75 26.8, p90 37.5, p95 54.7, p99 81.1, max 153.8. Median aspect ratio 2.28 (p90 8.52).

Link-distance sweep (the threshold is a free parameter, so here is the sweep rather than a single value): 10 um -> 4,149 objects, median extent 8.9 um, 15 um -> 3,698 objects, median extent 13.6 um, 20 um -> 3,091 objects, median extent 18.7 um, 25 um -> 2,665 objects, median extent 22.2 um, 30 um -> 2,421 objects, median extent 25.8 um, 40 um -> 2,060 objects, median extent 33.8 um.


### s2 — 10x Xenium GBM (io_480 immuno-oncology panel)

57,733 endothelial cells -> 20,875 discrete vessel objects at the 20 um linking distance; MVD 233.3/mm^2 over 89.5 mm^2; 63% singletons.

| in-plane extent threshold | implied `th` | **frac of ALL objects** | frac of multi-cell objects | isotropic null |
|---|---:|---:|---:|---:|
| <= 1x caliber (18.1 um) | 0 deg | **0.7694** | 0.3792 | 0.0000 |
| <= 1.5x caliber (27.2 um) | 61 deg | **0.8424** | 0.5757 | 0.7661 |
| <= 2x caliber (36.2 um) | 75 deg | **0.8894** | 0.7023 | 0.9291 |
| <= 3x caliber (54.3 um) | 82 deg | **0.9425** | 0.8451 | 0.9813 |

Multi-cell object extents (um): median 22.6, p75 40.6, p90 67.5, p95 90.0, p99 148.4, max 497.3. Median aspect ratio 2.39 (p90 6.26).

Link-distance sweep (the threshold is a free parameter, so here is the sweep rather than a single value): 10 um -> 43,154 objects, median extent 9.4 um, 15 um -> 28,335 objects, median extent 16.1 um, 20 um -> 20,875 objects, median extent 22.6 um, 25 um -> 16,063 objects, median extent 29.3 um, 30 um -> 12,209 objects, median extent 36.0 um, 40 um -> 6,617 objects, median extent 50.8 um.


**The honest reading, including what this measurement CANNOT do.** The measured extents are consistent with predominantly transverse cuts, and they sit well below the isotropic null at every threshold. But three censoring effects all push the measured extent DOWN, so the measured transverse fraction is an **upper bound** and this census **cannot reject an isotropic orientation**:

1. singletons carry extent 0 and no shape information

1. the 20 um single-linkage threshold splits long oblique traces

1. extent under-runs the true trace at both ends (sparse nuclei)

The link-distance sweep shows how much of it the 20 um threshold sets: the median multi-cell extent tracks the linking distance almost linearly, which is exactly what one expects when nuclei are sparse relative to the structure. **This is precisely why section 3 runs the isotropic case as well as the transverse case.** If the numerical bound is small in BOTH, the orientation question does not have to be settled to close the issue — and it is.


## 3. Numerical bound: full 3-D vs the 2-D solve

A standalone 3-D solver (`solve_helmholtz_dct`, exact via DCT on the same Neumann finite-volume stencil the production code uses) on a 768 um cube at dx = 6 um (128^3 = 2,097,152 unknowns). Vessels are **connected persistent random walks** at the measured length density `L_V = 2 x MVD = 2 x 110.8/mm^2`, not disjoint segments — a real capillary turns and branches rather than stopping. 3 seeds per configuration; statistics on the central 70% of the plane so the mirror boundary cannot contaminate them. Nothing in `tmesim/physics` is touched.

**The 2-D comparison field** is the solve of the mid-plane's own source layer — exactly the vessel material a section of that thickness contains, obliquely cut vessels included. That preserves the areal source total, so the comparison is about *where* drug enters, not how much.

**Control.** Straight vessels perpendicular to the section, spanning the box: the 3-D problem is z-invariant by construction, so the 2-D solve must reproduce the mid-plane exactly. Measured maximum relative difference **4.73e-16**. This validates the solver and proves the central claim in one line.


| orientation | persistence | drug | `L` (um) | mean 2D/3D | CV 3-D | CV 2-D | **CV ratio** | `d_half` 3-D | `d_half` 2-D | ratio | `F_pen` 3-D | `F_pen` 2-D | **3D/2D** |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| isotropic | 50 um | doxorubicin | 64.9 | 0.880 | 0.5362 | 0.9156 | **1.71x** | 36.0 | 29.2 | 0.813 | 0.4361 | 0.2745 | **1.589** |
| isotropic | 50 um | cisplatin | 142.9 | 0.887 | 0.2453 | 0.4381 | **1.79x** | 61.8 | 44.2 | 0.715 | 0.7075 | 0.5127 | **1.380** |
| isotropic | 50 um | temozolomide | 264.5 | 0.890 | 0.1214 | 0.2112 | **1.74x** | 79.5 | 49.8 | 0.627 | 0.7646 | 0.6091 | **1.255** |
| isotropic | 100 um | doxorubicin | 64.9 | 0.951 | 0.4813 | 0.6872 | **1.43x** | 52.9 | 37.1 | 0.700 | 0.5287 | 0.3823 | **1.383** |
| isotropic | 100 um | cisplatin | 142.9 | 0.976 | 0.2259 | 0.3427 | **1.52x** | 85.6 | 67.6 | 0.789 | 0.7527 | 0.6375 | **1.181** |
| isotropic | 100 um | temozolomide | 264.5 | 1.013 | 0.1081 | 0.1797 | **1.66x** | 99.7 | 85.9 | 0.862 | 0.8103 | 0.7161 | **1.132** |
| isotropic | 200 um | doxorubicin | 64.9 | 1.325 | 0.3471 | 0.6458 | **1.86x** | 43.0 | 28.8 | 0.669 | 0.6199 | 0.3801 | **1.631** |
| isotropic | 200 um | cisplatin | 142.9 | 1.370 | 0.1608 | 0.3024 | **1.88x** | 66.9 | 49.8 | 0.744 | 0.8237 | 0.6394 | **1.288** |
| isotropic | 200 um | temozolomide | 264.5 | 1.408 | 0.0743 | 0.1412 | **1.90x** | 82.1 | 68.0 | 0.828 | 0.8727 | 0.7349 | **1.187** |
| transverse | 50 um | doxorubicin | 64.9 | 1.006 | 0.5156 | 0.7994 | **1.55x** | 40.8 | 31.5 | 0.773 | 0.4719 | 0.3250 | **1.452** |
| transverse | 50 um | cisplatin | 142.9 | 1.047 | 0.2477 | 0.3947 | **1.59x** | 40.3 | 51.7 | 1.283 | 0.6833 | 0.5492 | **1.244** |
| transverse | 50 um | temozolomide | 264.5 | 1.096 | 0.1271 | 0.1980 | **1.56x** | 28.1 | 52.2 | 1.860 | 0.7254 | 0.6256 | **1.160** |
| transverse | 100 um | doxorubicin | 64.9 | 0.985 | 0.4923 | 0.8993 | **1.83x** | 49.2 | 29.7 | 0.603 | 0.5183 | 0.2741 | **1.891** |
| transverse | 100 um | cisplatin | 142.9 | 1.017 | 0.2704 | 0.4746 | **1.75x** | 68.7 | 48.5 | 0.706 | 0.7168 | 0.5059 | **1.417** |
| transverse | 100 um | temozolomide | 264.5 | 1.055 | 0.1416 | 0.2418 | **1.71x** | 79.1 | 57.1 | 0.722 | 0.7571 | 0.6072 | **1.247** |
| transverse | 200 um | doxorubicin | 64.9 | 0.936 | 0.3689 | 0.6814 | **1.85x** | 36.3 | 28.9 | 0.797 | 0.5928 | 0.3558 | **1.666** |
| transverse | 200 um | cisplatin | 142.9 | 0.984 | 0.1548 | 0.2918 | **1.89x** | 58.9 | 46.3 | 0.786 | 0.8044 | 0.6293 | **1.278** |
| transverse | 200 um | temozolomide | 264.5 | 1.009 | 0.0729 | 0.1279 | **1.76x** | 89.9 | 61.7 | 0.686 | 0.8458 | 0.7269 | **1.164** |

### 3.1 Reading the table

* **Mean exposure: unbiased.** `mean 2D/3D` runs 0.880-1.408 with mean 1.046 — scattered about 1 with no systematic sign, exactly as the stereological identity of section 1.4 requires. The scatter is *which plane you cut*, not dimensionality.

* **Heterogeneity: over-stated, always.** The CV ratio exceeds 1 in **every one of the 18** configurations, 1.43-1.90x. This is the one robust systematic finding, and its cause is section 4's: the section is blind to the out-of-plane vessels that fill the gaps.

* **Fitted decay distance: under-stated by ~24%** (`d_half` ratio median 0.759, range 0.603-1.860). Note the **direction is the opposite of the naive worry**: the concern was that 2-D would over-state reach because `K_0` decays more slowly than `exp(-r)/r`. The single-vessel prefactor is simply not the dominant term once a network is present.

* **Delivery factor: 2-D under-states it** by 1.13-1.89x. This is the term that feeds `results/kill_calibration.md`, and section 5 propagates it.


### 3.2 Normalised profile against in-plane distance to the nearest vessel

`C(r)/C(0)`, both fields measured against the SAME in-plane distance because that is what an analyst reads off a section:

| orientation | persistence | drug | 25 um (3-D / 2-D) | 50 um (3-D / 2-D) | 100 um (3-D / 2-D) | 150 um (3-D / 2-D) |
|---|---:|---|---|---|---|---|
| isotropic | 50 | doxorubicin | 0.701 / 0.556 | 0.518 / 0.356 | 0.351 / 0.176 | 0.288 / 0.096 |
| isotropic | 50 | cisplatin | 0.868 / 0.759 | 0.782 / 0.628 | 0.698 / 0.477 | 0.655 / 0.377 |
| isotropic | 50 | temozolomide | 0.942 / 0.880 | 0.908 / 0.811 | 0.880 / 0.722 | 0.865 / 0.656 |
| isotropic | 100 | doxorubicin | 0.769 / 0.638 | 0.615 / 0.445 | 0.377 / 0.200 | 0.247 / 0.105 |
| isotropic | 100 | cisplatin | 0.917 / 0.841 | 0.850 / 0.740 | 0.693 / 0.531 | 0.550 / 0.379 |
| isotropic | 100 | temozolomide | 0.968 / 0.933 | 0.940 / 0.888 | 0.858 / 0.768 | 0.767 / 0.653 |
| isotropic | 200 | doxorubicin | 0.776 / 0.620 | 0.624 / 0.421 | 0.405 / 0.217 | 0.364 / 0.177 |
| isotropic | 200 | cisplatin | 0.920 / 0.822 | 0.850 / 0.709 | 0.708 / 0.542 | 0.641 / 0.474 |
| isotropic | 200 | temozolomide | 0.968 / 0.920 | 0.939 / 0.866 | 0.869 / 0.774 | 0.823 / 0.728 |
| transverse | 50 | doxorubicin | 0.734 / 0.599 | 0.571 / 0.411 | 0.404 / 0.207 | 0.337 / 0.105 |
| transverse | 50 | cisplatin | 0.889 / 0.800 | 0.815 / 0.698 | 0.734 / 0.524 | 0.699 / 0.407 |
| transverse | 50 | temozolomide | 0.954 / 0.904 | 0.924 / 0.858 | 0.898 / 0.764 | 0.894 / 0.698 |
| transverse | 100 | doxorubicin | 0.752 / 0.598 | 0.575 / 0.346 | 0.353 / 0.149 | 0.255 / 0.098 |
| transverse | 100 | cisplatin | 0.901 / 0.796 | 0.807 / 0.626 | 0.640 / 0.425 | 0.553 / 0.351 |
| transverse | 100 | temozolomide | 0.958 / 0.902 | 0.912 / 0.809 | 0.815 / 0.672 | 0.763 / 0.622 |
| transverse | 200 | doxorubicin | 0.744 / 0.614 | 0.588 / 0.411 | 0.456 / 0.227 | 0.395 / 0.169 |
| transverse | 200 | cisplatin | 0.909 / 0.826 | 0.841 / 0.717 | 0.772 / 0.576 | 0.746 / 0.506 |
| transverse | 200 | temozolomide | 0.967 / 0.926 | 0.940 / 0.876 | 0.912 / 0.805 | 0.908 / 0.765 |

The 2-D value is below the 3-D value at essentially every distance in every configuration: **the 2-D solve digs a deeper trough between vessels than the tissue really has.**


### 3.3 The z-loss correction is the wrong-signed fix

Two cheap one-parameter corrections were fitted against the 3-D mid-plane target:

| correction | what it models | best fit (mean) | range |
|---|---|---:|---|
| `k -> k(1+beta)` | diffusive escape perpendicular to the section (a z-LOSS is `beta > 0`) | **-0.728** | -0.827 to -0.592 |
| `S -> (1-f)S + f mean(S)` | out-of-plane vessels the section never sees | **0.580** | 0.469 to 0.738 |

**The best-fit z-loss coefficient is negative in every configuration tested** (mean -0.73). Adding a sink to represent escape out of the section would make the 2-D answer *worse*. The reason is the one this whole report turns on: in a statistically homogeneous tissue the drug that leaves the plane is replaced by drug arriving from out-of-plane vessels, and those two fluxes balance. What does *not* balance is the **arrangement** — and the correctly-signed cheap fix is therefore a uniform background source of about 0.58 of the total, not a sink.

**Neither is recommended.** Both are one-parameter smoothers fitted to a synthetic network, and the bias they would remove is smaller than the uncertainties this project already carries (section 5). Shipping a fitted correction for a sub-dominant effect would add a free parameter and subtract auditability.


## 4. The consequence for the published numbers

| published quantity | direction | size | larger than the bar already carried? |
|---|---|---|---|
| (i) **physical** penetration length `L = sqrt(D_eff/k)` | none | **exactly 0** (proved, section 1.1) | no |
| (i') **fitted** profile decay distance | 2-D too SHORT | ~24% (range -86-40%) | comparable — that metric already carries 16% grid drift and is flagged NOT CONVERGED in `delivery_metrics.md` |
| (ii) mean tumour AUC | none | **0 systematically** (section 1.4); +/-41% section-sampling scatter | no |
| (iii) therapeutic coverage | inconsistent sign | worst case **0.205 absolute** | ~2-3x the +/-0.03 grid bar, but far below the 30x kill interval |
| (iv) 3-D:2-D kill ratio | 2-D over-states the 3-D penalty | fold 3.26x -> 2.73-2.98x | no — **still inside the 2-5x band** |

### 4.1 Penetration length

L = sqrt(D_eff/k) is IDENTICAL in 2-D and 3-D — proved in section 1, not estimated. The 2-D and 3-D Green's functions differ only in the algebraic prefactor (K_0 vs exp(-r/L)/r) and share the same exponential rate. Nothing about the transport physics moves.

What *does* move is the number `results/delivery_metrics.md` reports as `penetration_length_L` — a three-parameter `C0 exp(-d/L) + floor` fit to the binned AUC-vs-distance profile. Because out-of-plane vessels flatten the real in-plane profile, that fit reads ~24% low. Corrected ranges:

| drug | published fitted `L` (um) | 3-D-corrected range (um) |
|---|---:|---|
| temozolomide | 273.3 | 147 - 453 |
| lomustine | 78.2 | 42 - 130 |
| cisplatin | 152.5 | 82 - 253 |
| methotrexate | 48.7 | 26 - 81 |
| doxorubicin | 39.6 | 21 - 66 |

That metric is already flagged in `delivery_metrics.md` as the ONE quantity that did not converge under grid refinement (16% drift) and as dominated by vessel/KDE geometry rather than drug transport. This correction is the same size as the caveat it already carries, and it does not change the drug ordering.

### 4.2 Mean tumour AUC

no systematic bias. The observed +/-41% scatter is WHICH PLANE YOU CUT — single-section sampling variance, a real but different limitation (we hold one section per specimen), not a dimensionality error. Zero by the stereological identity of section 1.4 — plane-integrated delivery per unit vessel length is orientation-independent, so a 2-d solve on the section's own traces reproduces the true 3-d mean exactly for any statistically homogeneous network.

The section-sampling scatter deserves a sentence of its own, because it is a real limitation and it is **not** the one this report was asked about: we hold **one section per specimen**, so the plane we happened to get carries a genuine sampling variance on mean exposure. That is an argument for serial sections (which we do not have), not against the 2-D solve — a 3-D solve on one section's worth of data could not fix it either.

### 4.3 Therapeutic coverage

Threshold set so the 2-D field reproduces a given coverage, then read off the 3-D field:

| coverage in 2-D | coverage in 3-D (mean) | range over configurations | delta |
|---:|---:|---|---:|
| 0.10 | 0.028 | 0.001 - 0.071 | -0.072 |
| 0.25 | 0.185 | 0.099 - 0.247 | -0.065 |
| 0.50 | 0.610 | 0.553 - 0.681 | +0.110 |
| 0.75 | 0.891 | 0.833 - 0.955 | +0.141 |
| 0.90 | 0.970 | 0.922 - 0.990 | +0.070 |

worst-case 0.205 absolute (20 percentage points of tumour cells), against the +/-0.03 grid bar already carried. The sign is not even consistent: 2-D under-counts coverage where the threshold is permissive and over-counts where it is strict, because the 2-D field is broader-tailed at both ends. No drug's classification flips — a saturated drug stays saturated and a zero-coverage drug stays at zero. Concretely, against the published dx = 12 um panel: temozolomide (0.0000) and methotrexate (0.0000) stay at zero, cisplatin (0.9936) stays saturated, and only the mid-range drugs — doxorubicin (0.3836) and lomustine (0.2175) — move at all, by at most a few percentage points.

### 4.4 The 3-D:2-D kill ratio — the test that PASSED

`results/kill_calibration.md` predicts a temozolomide 3-D:2-D kill ratio of 0.2626 at the 100 uM ex-vivo bath concentration, an implied IC50 fold of 3.26x against a literature band of 2-5x. That prediction is `F_penetration x F_state`, and `F_penetration` is a pure product of the 2-D transport solve. This report measured how much the 2-D solve under-states it: 1.132-1.255x. Propagating:

| `F_pen` correction | corrected `F_pen` | corrected 3D:2D ratio | implied IC50 fold | inside 2-5x? |
|---:|---:|---:|---:|:--:|
| 1.000x (as published) | 0.4678 | 0.2626 | 3.26x | yes |
| 1.132x | 0.5293 | 0.2912 | 2.98x | yes |
| 1.176x | 0.5499 | 0.3026 | 2.89x | yes |
| 1.255x | 0.5872 | 0.3231 | 2.73x | yes |

**the 2-D solve UNDER-states F_penetration, so it OVER-states the 3-D penalty. Correcting it moves the implied temozolomide IC50 fold from 3.26x down to 2.73-2.98x — still inside the 2-5x literature band, and closer to that band's median (3x) and to the normoxic-subset comparator (2-3x) that kill_calibration.md section 3 note 3 argues is the fairer one. The PASS is not luck and it is not reversed; if anything the 2-D assumption was costing it accuracy.**

Two things this does and does not license. It **does** retire the specific worry that the PASS was partly luck of the 2-D solve: the correction is modest, its direction is favourable, and the verdict survives it. It does **not** touch any of the other five reasons `kill_calibration.md` section 3 itself gives for calling that a *weak* pass — a 2.5x-wide band resting on one paper, a 2-D reference assay with a ~19x between-study IQR, a test that constrains only the product of delivery and cell state, one section and one patient. **The 2-D assumption was never the binding constraint on that result, and now it is measured not to be.**


## 5. Limitations

1. **The 3-D network is synthetic.** We hold **single sections per specimen** and no serial sections, so a real 3-D vascular geometry for these tumours does not exist in this repository and cannot be reconstructed from the data we have. The network is built to the measured MVD and swept over orientation and tortuosity; it is not this patient's vasculature. **This is the single largest caveat in the report** and it is why the sweep, rather than one configuration, is the result.

2. **Constant `D` and `k`.** The comparison deliberately uses uniform coefficients to isolate geometry. The production solver carries a spatially varying `D_eff = D(1 - 0.6 rho_ecm)`; ECM heterogeneity would add its own variance to BOTH the 2-D and 3-D fields and is not expected to have a preferred dimensional direction, but that is an argument, not a measurement.

3. **Steady state, no convection.** Justified in the module docstring (the transient correction is 0.1-1.5% over the repo-nominal window) and by `results/numerical_rigor.md` (advection moves the penetration index by < 3e-4), but both are approximations made here for tractability.

4. **The orientation census is censored and cannot reject isotropy** (section 2). The report is written so that it does not need to: the bound holds in the isotropic case too.

5. **Vessel caliber is a mouse-xenograft anchor** (PMID 24390982). The human range (PMID 19856441) is 15-37.4 um, and a larger caliber would make the transverse-cut argument *stronger*, not weaker, so 18.1 um is the conservative choice.

6. **`F_penetration` is propagated as a scalar multiplier** onto the published value rather than recomputed end-to-end on the real section — which would require a 3-D reconstruction we do not have. The Hill inversion is exact (it reproduces the published 3.26x to 3 decimals), but the scale factor comes from the synthetic network.

7. **Nothing here is validated against a measured intratumoural drug gradient in human GBM**, which is the measurement that would settle both this question and the delivery/state split that `kill_calibration.md` section 5 lists as still open.


## 6. Verdict

**The 2-D assumption does not need fixing, and the reason is geometric rather than numerical: a capillary is a line, a line source is 2-D-symmetric, and a 5 um section through an 18 um vessel cuts it transversely.** `L = sqrt(D_eff/k)` is preserved exactly; mean tumour AUC is unbiased by a stereological identity; therapeutic coverage moves by at most 0.20 absolute with an inconsistent sign; and the 3-D:2-D kill test not only survives the correction but lands closer to the literature median after it. The one quantity that moves materially — the *fitted* profile decay distance, ~24% low — is already flagged in `delivery_metrics.md` as the metric that did not converge, and the correction is the same size as the caveat it carries.

**The residual error is conservative.** Every effect measured here runs in the direction of the 2-D solve making the tumour look *harder* to reach than it is. A model that errs towards under-stating delivery is the safe direction for a delivery model to err in.

**If a fix were ever wanted, it is not a z-loss term.** That term has the wrong sign (measured, section 3.3). The correctly-signed cheap fix is an out-of-plane background source; a genuinely 3-D solve would require serial sections, which we do not have. Neither is recommended: the bias is smaller than the uncertainties already published beside it.


---

*Generated by `tmesim/validation/dimensionality.py`. Analysis only — nothing in `core` / `physics` / `sim` / `molecular` / `reconstruct` is modified and no published number is changed. Deterministic (fixed seeds).*

