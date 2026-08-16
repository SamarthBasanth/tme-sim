# Mathematical model of `tmesim`

**The math the code actually solves.** Every equation below is cross-referenced to
the module and function that implements it. Where a term is a heuristic, a proxy,
or a fitted constant rather than derived physics, it is labelled as such inline —
no equation is stated here that the code does not solve.

**Global unit convention** (enforced throughout the package):

| quantity | unit | note |
|---|---|---|
| length | µm | grid spacing `Domain.dx`; molecular sub-modules work in nm and convert |
| time | s | |
| concentration | µM | drug; oxygen is dimensionless saturation in [0,1] |
| diffusivity | µm²/s | 1 m²/s = 10¹² µm²/s; 1 µm² = 10⁶ nm² |
| first-order rate | 1/s | `clearance`, `binding`, `k_ret` |
| pressure | mmHg | interstitial fluid pressure |
| hydraulic conductivity | µm²/(mmHg·s) | |
| temperature | K | 310 K body default |
| viscosity | Pa·s | |

Grid arrays are indexed `[row = y, col = x]` with shape `(ny, nx)`; cell centres
sit at $x_j = (j + \tfrac12)\,\Delta x$, $y_i = (i + \tfrac12)\,\Delta x$
(`core/domain.py::Domain.coords_um`). The mesh is isotropic ($\Delta y = \Delta x$).
Everything is 2-D by deliberate choice (`core/domain.py` docstring): reconstructing
a faithful 3-D microenvironment from sparse sampling is an unsolved inference
problem, and a whole-slide section is natively 2-D.

---

## 1. Multiscale overview

### 1.1 The two scales and the one-way coupling

`tmesim` is a **two-scale, one-way-coupled** model. The molecular layer is solved
once per drug and produces a small vector of scalars; the tissue layer then solves
PDEs on a fixed grid using those scalars as coefficients. There is no feedback from
tissue to molecule.

```
MOLECULE (structure)                                  TISSUE (grid, PDE)
  descriptors (MW, logP, TPSA, HBD)
    │
    ├─ Stokes–Einstein         → D_free  (µm²/s) ─┐
    ├─ Brownian dynamics φ_ref → D_micro          ├→ DrugParams.diffusivity  D
    ├─ Ogston–Renkin mesh gate → h_mesh  ─────────┘
    ├─ solubility–diffusion    → P (µm/s)   (diagnostic only)
    ├─ CNS-MPO rules           → bbb ∈ [0,1] ────→ DrugParams.bbb_penetrance
    ├─ logP/HBD heuristics     → k_b, k_cl (1/s) ─→ DrugParams.binding, .clearance
    ├─ binding-site barrier    → k_ret   (1/s) ───→ added into DrugParams.binding
    └─ potency lookup          → IC50 (µM) ──────→ DrugParams.ic50
                                                       │
                                                       ▼
                       oxygen PDE → cell states → pressure PDE → velocity
                                                       │
                                     transport PDE (D, v, S, k_cl+k_b) → AUC
                                                       │
                                                Hill kill → kill fraction
```

The **bridge** is `molecular/compile_drug.py::molecule_to_drugparams`, which
returns a `core.DrugParams` — the only object the tissue engine consumes:

$$
\texttt{DrugParams} = \bigl(D_{\mathrm{eff}},\; \mathrm{MW},\; \mathrm{IC}_{50},\; h,\; k_b,\; k_{\mathrm{cl}},\; \beta_{\mathrm{BBB}}\bigr)
$$

Defaults live in `core/config.py::DrugParams`
($D=250\ \mu\mathrm{m^2/s}$, $\mathrm{IC}_{50}=5\ \mu M$, $h=1.5$,
$k_b = 5\times10^{-4}$/s, $k_{\mathrm{cl}} = 2\times10^{-4}$/s, $\beta_{\mathrm{BBB}}=0.3$).

The **pipeline order** in `sim/engine.py::simulate_response` is: oxygen → cell-state
banding → pressure/velocity → delivery source → transport → Hill kill → readouts.
Each stage reads the previous stage's field; no stage is iterated to self-consistency
(oxygen does not respond to kill, pressure does not respond to drug).

### 1.2 Why atomistic MD is not an option (the scale argument)

Atomistic molecular dynamics integrates Newton's equations for every atom. Both
factors — atom count and step count — are out of reach for a tissue.

**Atoms.** Liquid water contains $\approx 3.3\times10^{22}$ molecules/cm³, i.e.
$\approx 1\times10^{11}$ atoms per µm³. One $(20\ \mu m)^3$ voxel — a single cell —
therefore holds

$$
(20\ \mu\mathrm{m})^3 \times 10^{11}\ \mathrm{atoms/\mu m^3}
= 8\times10^{3}\ \mu\mathrm{m^3} \times 10^{11}
\approx 10^{15}\ \text{atoms},
$$

already ~6 orders of magnitude beyond routine all-atom MD (millions of atoms). Even a
small *illustrative* domain of 120×120 at $\Delta x = 15\ \mu m$ — ~1.8 mm across —
is $\sim 10^{20}$–$10^{21}$ atoms.

> **That 120×120 figure is an illustration for this order-of-magnitude argument, not
> the production configuration.** Production runs are on the real measured sections at
> $\Delta x = 12\ \mu m$, which are **20–75× larger in area**:
>
> | section | grid | extent | bounding area |
> |---|---|---|---|
> | `xenium_319` (flagship) | (424, 675) | 8.10 × 5.09 mm | 41.2 mm² |
> | `xenium_gbm_io` (480-panel) | — | 11.5 × 10.9 mm | 124.6 mm² |
> | `gse324268` (366-panel) | — | 10.5 × 14.3 mm | 151.1 mm² |
>
> The `gse324268` row was measured directly from
> `GSM9572843_GBMtumor07_cells.csv.gz` (851,509 cells) and matches
> [`results/classifier_fix.md`](../results/classifier_fix.md) §s3. This note exists
> because the illustrative figure was previously stated in declarative voice and was
> the *only* domain-extent-in-mm number in this document — it misled a careful reader
> into sizing the production domain 4.5–8× too small per edge.

**Steps.** MD is locked to a $\sim1$–2 fs timestep by bond vibrations, while drug
transport and PD play out over $10^3$–$10^4$ s. Reaching a mere 1 ms needs
$10^{-3}/10^{-15} = 10^{12}$ steps; the biologically relevant window needs
$10^{18}$–$10^{19}$.

**The product.** $10^{15}$ atoms × $10^{12}$ steps $\approx 10^{27}$ atom-steps for
1 ms in one voxel. At an optimistic $10^9$ atom-steps/s that is $\sim10^{18}$ s
$\approx 10^{11}$ years. Whole-tissue atomistic MD is not expensive, it is
categorically not computable (see `docs/molecular_layer.md` §1.1).

**What is done instead.** Two honest, computable methods, kept separate:
(i) closed-form physicochemical relations from structure (Stokes–Einstein,
Ogston–Renkin, solubility–diffusion), and (ii) genuine **coarse-grained** particle
simulation at a computable scale — overdamped Langevin tracers through an obstacle
field, which is transport physics without atoms (`molecular/brownian.py`,
`molecular/ecm_mesh.py`). The molecular layer's output is a *coefficient*, not a
trajectory.

---

## 2. Molecular scale

### 2.1 Molecular weight → hydrodynamic radius

`molecular/stokes_einstein.py::hydrodynamic_radius_nm`. Two empirical regimes,
selected by `is_antibody` **or** $\mathrm{MW} \ge 5000$ g/mol:

$$
r_s\,[\mathrm{nm}] =
\begin{cases}
0.066\; \mathrm{MW}^{1/3}, & \text{small molecules (constant-density sphere, } V \propto \mathrm{MW}) \\[4pt]
0.0515\; \mathrm{MW}^{0.392}, & \text{globular proteins / antibodies}
\end{cases}
$$

with MW in g/mol. Checks: MW = 300 → $r_s \approx 0.44$ nm; MW = 150 kDa →
$r_s \approx 5.3$ nm (the known IgG value). **Status: empirical fit**, medium
confidence — the prefactors are literature-calibrated, not derived.

### 2.2 Stokes–Einstein free diffusivity

`molecular/stokes_einstein.py::free_diffusivity`. For a sphere in a continuum
fluid at low Reynolds number,

$$
D_{\mathrm{free}} = \frac{k_B T}{6\pi\eta r_s}
$$

| symbol | meaning | unit | default |
|---|---|---|---|
| $k_B$ | Boltzmann constant | J/K | $1.380649\times10^{-23}$ (exact) |
| $T$ | absolute temperature | K | 310 |
| $\eta$ | dynamic viscosity of interstitial fluid | Pa·s | $1.2\times10^{-3}$ |
| $r_s$ | hydrodynamic radius | m | from §2.1, converted nm → m |
| $D_{\mathrm{free}}$ | dilute-fluid diffusivity | µm²/s | reported after $\times 10^{12}$ |

Dimensional check: $[k_BT] = \mathrm{kg\,m^2\,s^{-2}}$,
$[6\pi\eta r] = \mathrm{(Pa\,s)\,m} = \mathrm{kg\,s^{-1}}$, so
$[D] = \mathrm{m^2/s}$; convert with $1\ \mathrm{m^2/s} = 10^{12}\ \mu\mathrm{m^2/s}$.

**Assumptions/limits:** rigid sphere, no-slip, continuum solvent, infinite dilution.
$D \propto 1/r_s$ and $D \propto 1/\eta$ exactly (asserted in the module self-test).
A ~150 kDa IgG gives $D_{\mathrm{free}} \approx 34\ \mu\mathrm{m^2/s}$ — the module
explicitly **refuses** to shrink the radius to manufacture the ~1–10 µm²/s tissue
value quoted for antibodies; that drop comes from ECM hindrance downstream.

### 2.3 Overdamped Langevin dynamics with excluded volume

`molecular/brownian.py::simulate_effective_diffusivity` (µm scale, ~300 nm
obstacles) and `molecular/ecm_mesh.py::brownian_mesh_hindrance` (nm scale, fibre
cross-sections). Both integrate the **inertia-free (overdamped) Langevin equation**

$$
\mathrm{d}\mathbf{x} = \sqrt{2 D_{\mathrm{free}}\,\mathrm{d}t}\;\boldsymbol{\xi},
\qquad \boldsymbol{\xi}\sim\mathcal{N}(\mathbf{0},\mathbf{I}_2)
$$

i.e. the exact Euler–Maruyama update for free Brownian motion, with per-axis step
standard deviation $\sigma_{\mathrm{step}} = \sqrt{2 D_{\mathrm{free}}\Delta t}$.

**Excluded volume (the only interaction).** Obstacles are impenetrable disks of
radius $R_{\mathrm{obs}}$ ($0.30$ µm in `brownian.py`; the fibre radius
$r_f = 2$ nm in `ecm_mesh.py`) placed by a **Poisson process** (i.i.d. uniform
centres, overlaps allowed). A proposed move whose centre lands within the excluded
distance

$$
r_{\mathrm{excl}} = R_{\mathrm{obs}} + r_s
$$

of any obstacle centre (minimum-image periodic distance) is **rejected** and the
tracer stays put for that step. This is Saxton's rejection scheme (Biophys J 1994);
it preserves the uniform equilibrium measure over the accessible area and yields
proper diffusive motion with $D_{\mathrm{eff}} \le D_{\mathrm{free}}$.

**Boundary conditions:** periodic box of side $L_{\mathrm{box}}$ (20 µm default in
`brownian.py`, 80 nm in `ecm_mesh.py`). The *wrapped* position is used for the
collision test; the *unwrapped* position is accumulated separately so wrapping
cannot corrupt the MSD.

**Initial condition:** tracers placed uniformly in the accessible area by rejection
sampling.

**Timestep choice** (automatic when `dt_s=None`): the rms per-axis step is set to
4% of the obstacle radius so obstacles are well resolved and cannot be tunnelled,

$$
\Delta t = \frac{(0.04\,R_{\mathrm{obs}})^2}{2 D_{\mathrm{free}}}.
$$

**Estimator.** The mean-squared displacement uses **time-origin averaging** over all
stored frames (variance reduction $\sim$ number of origins):

$$
\mathrm{MSD}(\tau) = \bigl\langle \lVert \mathbf{r}(t+\tau) - \mathbf{r}(t) \rVert^2 \bigr\rangle_{t,\,\text{particles}}
\;=\; 4 D_{\mathrm{eff}} \tau \quad \text{(2-D)},
$$

the factor 4 being $2d D$ with $d = 2$. The code fits an ordinary least-squares line
to $\mathrm{MSD}$ vs $\tau$ over the **middle half of the lag range**
($\tau$ index in $[n/4,\,3n/4)$) — skipping the short-time crossover where an
obstructed tracer is still locally free, and the long-lag tail where few origins
remain — and reports

$$
D_{\mathrm{eff}} = \max\!\left(\frac{\text{slope}}{4},\,0\right),
\qquad h_{\mathrm{BD}} = \frac{D_{\mathrm{eff}}}{D_{\mathrm{free}}} \in [0,1].
$$

**Obstacle bookkeeping.** `obstacle_fraction` $\phi$ is the *nominal* area fraction,
so the obstacle count is $N = \mathrm{round}(\phi L_{\mathrm{box}}^2 / \pi R_{\mathrm{obs}}^2)$
and, because placement is Poisson (overlaps allowed), the actually covered fraction
is the Boolean/Swiss-cheese value

$$
\phi_{\mathrm{covered}} = 1 - e^{-\phi},
$$

which coincides with $\phi$ only in the dilute limit. Verified behaviour: $R^2 > 0.99$
for the MSD fit at $\phi = 0$ and $\phi = 0.4$ (diffusive, not sub-diffusive);
$D_{\mathrm{eff}}$ recovers $D_{\mathrm{free}}$ within 5% at $\phi = 0$; hindrance is
monotone decreasing over $\phi \in [0, 0.5]$; the dilute slope
$h \approx 1 - a\phi$ has $a > 0.5$ (the impenetrable effective-medium limit
$(1-\phi)/(1+\phi)$ would give $a = 2$).

**Limits:** 2-D, point-or-sphere tracer, static obstacles, no hydrodynamic
interactions, no electrostatics, no binding. The obstacle geometry is a *stylised
proxy* for ECM, not a reconstruction of it.

### 2.4 Ogston partition and Renkin/Faxén hindrance (the fibre-matrix size gate)

`molecular/ecm_mesh.py`. The ~300 nm obstacles of §2.3 barely differentiate a 0.4 nm
drug from a 5.5 nm antibody — both are effectively points. The real size gate is the
ECM **fibre mesh** (collagen/proteoglycan fibres $r_f \approx 2$ nm at volume
fraction $\phi_f \approx 0.05$ in loose interstitium).

**Ogston (1958) equilibrium partition** — the fraction of void space whose centre a
sphere of radius $r_s$ can occupy without overlapping a fibre of radius $r_f$
(`ogston_partition`):

$$
\Phi(r_s) = \exp\!\left[-\phi_f\left(\frac{r_s + r_f}{r_f}\right)^{\!2}\right]
= \exp\!\left[-\phi_f\left(1 + \frac{r_s}{r_f}\right)^{\!2}\right] \in [0,1]
$$

Monotone decreasing in both $r_s$ and $\phi_f$; $\to 1$ as either $\to 0$.

**Characteristic pore radius** (`mesh_pore_radius_nm`) is defined *self-consistently
from the same law*: the solute radius at which $\Phi = 1/e$,

$$
\phi_f\left(1 + \frac{r_{\mathrm{pore}}}{r_f}\right)^{2} = 1
\;\Longrightarrow\;
r_{\mathrm{pore}} = r_f\left(\frac{1}{\sqrt{\phi_f}} - 1\right)
$$

($r_f = 2$ nm, $\phi_f = 0.05$ → $r_{\mathrm{pore}} \approx 6.9$ nm, mesh diameter
~14 nm), floored at $10^{-6}$ nm so $\lambda$ stays finite.

**Renkin (1954) / Faxén hydrodynamic drag** (`_faxen_drag`), as compiled by Deen
(1987), with $\lambda = r_s / r_{\mathrm{pore}}$:

$$
K_d(\lambda) =
\begin{cases}
1 - 2.104\lambda + 2.09\lambda^3 - 0.95\lambda^5, & 0 \le \lambda < 1 \\
0, & \lambda \ge 1
\end{cases}
\qquad \text{clamped to } [0,1]
$$

Note the code deliberately uses **only the drag polynomial**, not Renkin's own
steric $(1-\lambda)^2$ factor, because the steric exclusion is already supplied by
the Ogston partition — including both would double-count.

**Combined analytic hindrance** (`ogston_renkin_hindrance`):

$$
h_{\mathrm{mesh}}(r_s) = \Phi(r_s)\, K_d(\lambda) \in [0,1]
$$

— fewer accessible sites **and** more near-wall drag. Behaviour asserted in the
self-test: a 5.5 nm IgG in a tight mesh ($\phi_f = 0.30$) has $\Phi < 0.05$; the
small-molecule/antibody hindrance ratio exceeds 5× even in loose matrix, and the
independent nm-scale Brownian check (§2.3 at fibre scale) agrees in trend, sharpening
to $\ge 2\times$ with near-total antibody exclusion ($\phi_f = 0.24$).

### 2.5 Composite effective diffusivity handed to the tissue engine

`molecular/compile_drug.py::build_profile` steps 3 and 3b:

$$
D_{\mathrm{eff}}^{\mathrm{mol}}
= \underbrace{D_{\mathrm{BD}}\!\left(D_{\mathrm{free}},\,\phi_{\mathrm{ref}},\,r_s\right)}_{\text{µm-scale tortuosity (BD)}}
\;\times\;
\underbrace{h_{\mathrm{mesh}}(r_s)}_{\text{nm-scale size gate}},
\qquad
h = \frac{D_{\mathrm{eff}}^{\mathrm{mol}}}{D_{\mathrm{free}}}
$$

with $\phi_{\mathrm{ref}} = 0.4$ by default (the "ECM reference"). The two factors
address different physics (path tortuosity vs. steric/hydrodynamic sieving), which
is why the product is taken rather than either alone. **Honest note:** the product is
a modelling choice, not a theorem — the two hindrances are assumed statistically
independent.

> **Status: PARTIALLY SUPPORTED, and this equation is no longer the whole story.**
> Both factors above depend on the solute only through its hydrodynamic radius,
> $0.066\,\mathrm{MW}^{1/3}$, so on their own they make
> $D_{\mathrm{eff}}^{\mathrm{mol}}$ a function of molecular weight and the ECM
> fraction and nothing else. That model was scored blind against 71 compounds
> (Pruijn 2008, HT29 multicellular layers) and **FALSIFIED**: the measurement spans
> **29.7×** where it spans **2.4×**, at **R² = 0.022**, and at identical MW it is
> flat by construction (**exactly 1.00×**) against a measured **3.7×**.
>
> Since 2026-08-13 the product above is multiplied by a **two-route cell term**
> (§2.5b), which brings $\log P$ and hydrogen-bond donor count into $D_{\mathrm{eff}}$
> through a membrane permeability. On the same 71 compounds, blind and unfitted, the
> composed law reaches **R² = 0.718**, Spearman **+0.821**, **99 %** within 3× of
> measurement, and orders **24 of 31** exact-mass pairs correctly against the old
> model's **0 of 31**. It nonetheless reproduces **8.3×** of the measured **29.7×**
> spread — about a **quarter** — so it *under-separates* a similar-mass series
> systematically, and the whole test is **one cell line**.
> → [`../results/two_route_promotion.md`](../results/two_route_promotion.md),
> [`../results/diffusivity_falsification.md`](../results/diffusivity_falsification.md)

### 2.6 Solubility–diffusion membrane permeability

`molecular/permeability.py::membrane_permeability`. Overton/Meyer picture: the solute
partitions into the bilayer, diffuses across the hydrophobic core, and partitions out.

$$
P = \frac{K_p\, D_m}{d},
\qquad
K_p = 10^{\,\mathrm{clamp}(\log P,\,-3,\,5)},
\qquad
D_m = \frac{k_B T}{6\pi \eta_m r_s}
$$

| symbol | meaning | unit | value |
|---|---|---|---|
| $P$ | passive permeability | µm/s | returned ($\times 10^6$ from m/s) |
| $K_p$ | membrane/water partition coefficient | — | $10^{\log P}$, clamped to $[10^{-3}, 10^{5}]$ |
| $D_m$ | solute diffusivity inside the bilayer core | m²/s | Stokes–Einstein with membrane viscosity |
| $\eta_m$ | bilayer-core microviscosity | Pa·s | $0.10$ (≫ water's $10^{-3}$) |
| $d$ | hydrophobic bilayer thickness | m | $4.0$ nm |
| $r_s$ | solute radius | m | `radius_hydrodynamic_nm`, else from MW at $\rho = 1$ g/cm³ |

Units: $[D_m/d] = \mathrm{m/s}$, $K_p$ dimensionless. Antibodies (`is_antibody`)
return exactly $0$ — macromolecules do not cross by passive diffusion.

**Status: low–medium confidence.** This is *passive* physics only. Active transport,
carrier-mediated uptake, P-gp efflux and metabolism all move real drugs off this
trend line, and the $\log P$ clamp is a numerical guard, not physics.

### 2.7 BBB penetrance (CNS-MPO-style rules — heuristic)

`molecular/permeability.py::bbb_penetrance`, returning $\beta_{\mathrm{BBB}} \in [0,1]$.
Four smooth desirabilities, their **geometric** mean, then a logistic squash.
(This section previously printed an *arithmetic* mean, $\tfrac14\sum d_i$. That
was the pre-fix form and is retracted — the code has computed a geometric mean
since the veto fix; see §8.3. A geometric mean is a **veto**: one disqualifying
descriptor collapses the score, which an average cannot express.)

$$
d_{\mathrm{MW}} = \frac{1}{1+e^{(\mathrm{MW}-450)/60}},\quad
d_{\log P} = e^{-(\log P - 2.5)^2 / (2\cdot 1.7^2)},
$$
$$
d_{\mathrm{TPSA}} = \frac{1}{1+e^{(\mathrm{TPSA}-90)/12}},\quad
d_{\mathrm{HBD}} = \frac{1}{1+e^{(\mathrm{HBD}-3)/0.8}},
$$
$$
\bar d = \Bigl(d_{\mathrm{MW}}\,d_{\log P}\,d_{\mathrm{TPSA}}\,d_{\mathrm{HBD}}\Bigr)^{1/4},
\qquad
\beta_{\mathrm{BBB}} = \frac{1}{1+e^{-8(\bar d - 0.5)}}
$$

Hard gates: $\beta_{\mathrm{BBB}} = 0$ for antibodies and for $\mathrm{MW}\ge1000$ g/mol.

**Status: heuristic, low confidence, with a documented failure.** Temozolomide
(MW ≈ 194, $\log P \approx -1.5$, TPSA ≈ 106) is clinically BBB-penetrant but is
*under-predicted* by any purely physicochemical rule set. The module does **not**
special-case it. `results/pk_multidrug.md` reports the measured consequence: BBB
ordering vs. real CNS-partition ordering gives Spearman $\rho = 0.50$, Kendall
$\tau = 0.33$ (2/3 pairs concordant) — partially concordant, inverting
temozolomide vs. cisplatin.

### 2.8 Plasma binding, tissue retention and clearance (heuristic priors)

`molecular/binding.py`. All three are **low-confidence structure-only priors** whose
*direction* and *order of magnitude* are defensible; their absolute values are not.

Free fraction (albumin-binding proxy, `plasma_free_fraction`):

$$
f_u = \frac{1}{1 + 10^{(\log P - 3.0)/1.5}} \;\in [0.01,\,1.0]\ \text{(clamped)}
$$

Tissue binding rate (`tissue_binding_rate`) and clearance (`clearance_rate`) map a
dimensionless score through a logistic and then **geometrically** onto a documented
band ($\text{lo}\cdot(\text{hi}/\text{lo})^{s}$ — equal score changes give equal
*factor* changes, the natural scale for first-order rate constants):

$$
s_b = \sigma\!\left(\frac{\log P - 2.5}{2.5} + \frac{\mathrm{HBD}-2}{4}\right),
\qquad
k_b = 10^{-5}\left(\frac{5\times10^{-3}}{10^{-5}}\right)^{s_b}\ \mathrm{s^{-1}}
$$

$$
s_c = \sigma\!\left(0.4\,\frac{\log P - 2}{3} - 0.4\,\frac{\mathrm{MW}-300}{700}\right),
\qquad
k_{\mathrm{cl}} = 10^{-5}\left(\frac{10^{-3}}{10^{-5}}\right)^{s_c}\ \mathrm{s^{-1}}
$$

with $\sigma$ the logistic sigmoid. Bands: $k_b \in [10^{-5}, 5\times10^{-3}]$ 1/s,
$k_{\mathrm{cl}} \in [10^{-5}, 10^{-3}]$ 1/s, $f_u \in [0.01, 1]$. Monotonicity in
$\log P$ (and in HBD for $k_b$) is asserted in the self-test. These feed
`DrugParams.binding` and `DrugParams.clearance` — and $k_b$ is later reused as the
*uptake proxy* the calibration constant scales (§8).

---

## 3. Tissue transport: convection–diffusion–reaction

`physics/diffusion.py::simulate_drug_transport`. This is the numerical heart of the
package.

### 3.1 Governing equation

For the drug concentration field $C(x,y,t)$ [µM]:

$$
\frac{\partial C}{\partial t}
= \nabla\!\cdot\!\bigl(D_{\mathrm{eff}}(x,y)\,\nabla C\bigr)
\;-\; \nabla\!\cdot\!\bigl(\mathbf{v}\,C\bigr)
\;+\; S(x,y)
\;-\; \bigl(k_{\mathrm{cl}} + k_b\bigr)\,C
$$

| symbol | meaning | unit |
|---|---|---|
| $C$ | drug concentration | µM |
| $D_{\mathrm{eff}}$ | spatially varying effective diffusivity | µm²/s |
| $\mathbf{v} = (v_x, v_y)$ | interstitial (Darcy) fluid velocity, from §5 | µm/s |
| $S$ | vascular delivery rate | µM/s |
| $k_{\mathrm{cl}}$ | first-order clearance/metabolism | 1/s |
| $k_b$ | first-order tissue binding/retention | 1/s |

The advection term is in **conservative (divergence) form** $\nabla\cdot(\mathbf{v}C)$,
not $\mathbf{v}\cdot\nabla C$; this is what makes the discrete scheme conservative
even where $\nabla\cdot\mathbf{v} \ne 0$.

### 3.2 ECM tortuosity coefficient (heuristic)

`_effective_diffusivity`:

$$
D_{\mathrm{eff}}(x,y) = \max\Bigl(D\,\bigl(1 - \alpha_{\mathrm{ecm}}\,\rho_{\mathrm{ecm}}(x,y)\bigr),\; \max(10^{-6},\,10^{-3}D)\Bigr),
\qquad \alpha_{\mathrm{ecm}} = 0.6
$$

with $\rho_{\mathrm{ecm}} \in [0,1]$ the ECM density field (clipped) and $D$ the
drug's compiled diffusivity. The floor keeps the operator well-posed in dense-ECM
regions where the linear expression could reach zero.

**Status: heuristic.** $\alpha_{\mathrm{ecm}} = 0.6$ is a linear attenuation
coefficient chosen for the right qualitative behaviour (higher ECM → lower $D$ →
steeper gradients → lower penetration index, asserted in
`_test_penetration_monotonicity`). It is **not** derived from the Ogston–Renkin
theory of §2.4. **Double-counting caveat:** the compiled $D$ already contains a
mesh hindrance evaluated at a *fixed reference* $\phi_{\mathrm{ref}} = 0.4$, so the
tissue-scale factor $(1-0.6\rho_{\mathrm{ecm}})$ should be read as a *relative
spatial modulation about that reference*, not as an independent physical hindrance.

### 3.3 Finite-volume diffusion operator with harmonic-mean faces

`_build_diffusion_operator` assembles $\mathbf{L} \in \mathbb{R}^{N\times N}$,
$N = n_y n_x$, in **finite-volume** form. For cell $P$ with neighbours $\mathcal{N}(P)$:

$$
(\mathbf{L}C)_P = \frac{1}{\Delta x^{2}} \sum_{f \in \mathcal{N}(P)} D_f \,\bigl(C_{N_f} - C_P\bigr)
$$

The **face diffusivity is the harmonic mean** of the two adjacent cell values:

$$
D_f = \frac{2\,D_P\,D_N}{D_P + D_N}
$$

The harmonic mean is the correct face value for a *series* resistance: two cells of
thickness $\Delta x/2$ each contribute resistance $\propto 1/D$, so the conductance
of the pair is the harmonic mean. An arithmetic mean would let a high-$D$ cell
short-circuit a low-$D$ neighbour and would over-predict flux across sharp ECM
interfaces.

**Structure of $\mathbf{L}$:** symmetric (each face contributes $+D_f/\Delta x^2$ to
both off-diagonals $(P,N)$ and $(N,P)$, and $-D_f/\Delta x^2$ to both diagonals), and
**every row sums to exactly zero**. Rows of non-tissue cells are identically zero,
so those cells are decoupled and held at $C=0$.

### 3.4 Boundary conditions: no-flux everywhere

$$
\bigl(D_{\mathrm{eff}}\nabla C\bigr)\cdot \hat{\mathbf{n}} = 0
\quad \text{on } \partial\Omega \ \text{ and on every tissue/non-tissue interface}
$$

Implementation: the face coefficient is set to zero unless **both** adjacent cells
are tissue,

$$
D_f \leftarrow D_f \cdot \mathbb{1}\bigl[\texttt{tissue}_P \wedge \texttt{tissue}_N\bigr],
$$

and faces to the domain exterior simply do not exist (they are never assembled).
Zeroing the face coefficient *is* the discrete no-flux condition — no ghost cells,
no special-casing. The same masking is applied to the advective fluxes and to the
source ($S \leftarrow S\cdot\mathbb{1}[\texttt{tissue}]$) and reaction rate
($k \leftarrow 0$ outside tissue).

**Initial condition:** $C(x,y,0) = 0$ everywhere (`domain.zeros()`) — no drug in
tissue before dosing.

### 3.5 Upwind advection and the CFL condition

`_advective_divergence` computes a first-order **donor-cell (upwind)** approximation
of $\nabla\cdot(\mathbf{v}C)$. On the $x$-face between $(i,j)$ and $(i,j+1)$:

$$
v_f = \tfrac12\bigl(v_{x,i,j} + v_{x,i,j+1}\bigr),
\qquad
C_{\mathrm{up}} = \begin{cases} C_{i,j}, & v_f > 0\\ C_{i,j+1}, & v_f \le 0\end{cases},
\qquad
F_f = \frac{v_f\, C_{\mathrm{up}}}{\Delta x}
$$

and $F_f$ is **added** to cell $(i,j)$ and **subtracted** from cell $(i,j+1)$ — the
same face flux with opposite signs, which is precisely discrete antisymmetry. The
$y$-faces are treated identically. Masked faces contribute $F_f = 0$.

Upwinding (rather than central differencing) guarantees a positive, monotone scheme
at any cell Péclet number; its cost is $O(\Delta x)$ numerical diffusion.

Because advection is the **only explicit** piece, it carries a CFL restriction. The
code computes and *enforces* it, raising `ValueError` otherwise:

$$
\mathrm{CFL} = \frac{\Delta t\,\max_{i,j}\bigl(|v_x| + |v_y|\bigr)}{\Delta x} \le 1
$$

Interstitial speeds are $\lesssim 1$ µm/s, so with $\Delta x \sim 10$–30 µm and
$\Delta t \sim 10$ s this is satisfied with wide margin. (`results/numerical_rigor.md`
records that advection changes the penetration index by $<3\times10^{-4}$ at
realistic speeds, which is why the numerical audits run with convection off.)

### 3.6 The IMEX operator split

Diffusion **and** the linear reaction are treated implicitly (backward Euler);
advection and the source are explicit:

$$
\bigl(\mathbf{I} - \Delta t\,\mathbf{L} + \Delta t\,\mathbf{K}\bigr)\,C^{n+1}
\;=\;
C^{n} + \Delta t\bigl(S - \mathrm{Adv}(C^{n})\bigr)
$$

with $\mathbf{K} = \mathrm{diag}(k_{\mathrm{cl}} + k_b)$ (zeroed outside tissue) and
$\mathrm{Adv}$ the upwind divergence of §3.5.

**Why the diffusion part is unconditionally stable.** $-\mathbf{L}$ is symmetric
positive semi-definite (it is a weighted graph Laplacian with non-negative edge
weights $D_f/\Delta x^2$), and $\mathbf{K}\succeq 0$. Hence for any $\Delta t > 0$

$$
\mathbf{A} = \mathbf{I} - \Delta t\,\mathbf{L} + \Delta t\,\mathbf{K}
$$

is symmetric positive definite and weakly diagonally dominant with positive
diagonal, so $\mathbf{A}^{-1}$ exists and is non-negative (M-matrix). Every
eigenvalue of the update operator $\mathbf{A}^{-1}$ lies in $(0,1]$, so the scheme
is **unconditionally stable and monotone**: there is *no* $\Delta t \le \Delta x^2/(4D)$
restriction, which is exactly what makes hour-long windows on a 120×120 grid cheap.

**Cost.** $\mathbf{A}$ is time-independent, so it is LU-factorised **once**
(`scipy.sparse.linalg.factorized`) and the factors are reused for every step;
`_test_runtime_120` asserts a 120×120 grid with convection on and 1080 steps
completes in < 10 s.

After each solve the code applies `np.maximum(c_new, 0.0)` to remove
round-off/upwind undershoot and re-zeros non-tissue cells.

### 3.7 Exact discrete mass conservation

Take $k_{\mathrm{cl}} = k_b = 0$, $S = 0$, no advection. Because every row of
$\mathbf{L}$ sums to zero, $\mathbf{1}^{\mathsf{T}}\mathbf{L} = \mathbf{0}$, hence
$\mathbf{1}^{\mathsf{T}}\mathbf{A} = \mathbf{1}^{\mathsf{T}}$ and therefore
$\mathbf{1}^{\mathsf{T}}\mathbf{A}^{-1} = \mathbf{1}^{\mathsf{T}}$. From
$\mathbf{A}C^{n+1} = C^{n}$:

$$
\mathbf{1}^{\mathsf{T}} C^{n+1} = \mathbf{1}^{\mathsf{T}}\mathbf{A}^{-1} C^{n} = \mathbf{1}^{\mathsf{T}} C^{n}
$$

so the total tissue mass is invariant **to the linear-solver tolerance**, not merely
to truncation error. The upwind advection uses antisymmetric face fluxes that vanish
on all masked/boundary faces, so it too conserves the tissue total. With a source and
no loss, the final mass must equal exactly $n_{\mathrm{steps}}\,\Delta t \sum S$.
Measured: $\sim10^{-13}$–$10^{-15}$ relative error on the real irregular mask (§10).

### 3.8 Exposure (AUC) and the penetration index

The exposure integral is accumulated with the **trapezoidal rule** each step:

$$
\mathrm{AUC}(x,y) = \int_0^{T} C(x,y,t)\,\mathrm{d}t
\;\approx\; \sum_{n} \tfrac12\bigl(C^{n} + C^{n+1}\bigr)\Delta t
\qquad [\mu M\cdot s]
$$

The **penetration index** is a dimensionless flatness statistic over tissue cells:

$$
\Pi = \mathrm{clip}\!\left(\frac{\overline{\mathrm{AUC}}}{\mathrm{AUC}_{p99}},\,0,\,1\right),
\qquad
\mathrm{AUC}_{p99} = \text{99th percentile of AUC over tissue}
$$

$\Pi \to 1$ for a well-mixed (perfectly penetrating) field where mean = peak;
$\Pi \to 0$ when drug is confined to thin perivascular cuffs. The 99th percentile
(rather than the max) is a robustness choice against single-cell outliers.

**Status: a diagnostic, not a physical length.** $\Pi$ has no units and no direct
experimental counterpart; the physically comparable quantity is the penetration
length $L$ of §7. $\Pi$ is also **resolution-sensitive** — see §10.2.

### 3.9 The delivery source: plasma-averaged (default) or time-resolved

`treatment/drugs.py::plasma_profile` gives the systemic 1-compartment multi-dose
profile,

$$
C_{\mathrm{plasma}}(t) = C_{\mathrm{peak}} \sum_{k=0}^{n_d-1} H\!\left(t - k\,\Delta t_{\mathrm{dose}}\right)
\exp\!\left[-\lambda\,\bigl(t - k\,\Delta t_{\mathrm{dose}}\bigr)\right],
\qquad \lambda = \frac{\ln 2}{t_{1/2}}
$$

($H$ = Heaviside; doses superpose linearly, so overlapping doses accumulate), and
`delivery_source` gives the linear vascular-leak source

$$
S(x,y) = \rho_{\mathrm{vasc}}(x,y)\;\beta_{\mathrm{BBB}}\;C_{\mathrm{plasma}}
\qquad [\mu M/s]
$$

**Default — plasma-averaging (`PDEConfig.time_resolved = False`).** The transport
solver applies a *static* $S$ every step, i.e. $C_{\mathrm{plasma}}$ is replaced by
its time-average over the window,

$$
\bar C_{\mathrm{plasma}} = \frac{1}{T}\int_0^T C_{\mathrm{plasma}}(t)\,\mathrm{d}t .
$$

This preserves the *total* delivered drug (mean level × window) and the relative
ranking of drugs/patients, at the cost of the true within-window transient. Its
sharp edge: **any two dosing schedules with the same cumulative plasma exposure
become bit-identical inputs**, so fractionation and infusion duration are
unresolvable *by construction* rather than as a finding. Historical results in
`results/` predating the option below were produced in this mode and any
"schedule has no effect" statement in them is an artefact of it.

**Time-resolved (`PDEConfig.time_resolved = True`).** The source carries its true
time dependence, evaluated on the solver's own step grid at the old time level of
each step (matching the explicit convention already used for advection):

$$
S(x,y,t_n) = \rho_{\mathrm{vasc}}(x,y)\;\beta_{\mathrm{BBB}}\;C_{\mathrm{plasma}}(t_n),
\qquad t_n = n\,\Delta t,\quad n = 0,\dots,N-1
$$

Implemented as the `source_profile` argument of
`physics/diffusion.py::simulate_drug_transport`. The implicit operator
$A = I - \Delta t\,L + \Delta t\,K$ contains no source term, so the LU
factorisation is still computed once and reused — a time-varying source costs one
extra scalar multiply per step. Because the profile is sampled on the step grid,
$\frac{1}{N}\sum_n C_{\mathrm{plasma}}(t_n) = \bar C_{\mathrm{plasma}}$ exactly, so
the two modes deliver **identical total drug** and differ only in timing.

**What timing can and cannot change.** Integrating the transport PDE over $[0,T]$
with $C(0)=0$, and writing $A = \int_0^T C\,\mathrm{d}t$ and
$M = k - \nabla\!\cdot\!(D\nabla\,\cdot) + \nabla\!\cdot\!(v\,\cdot)$:

$$
A = M^{-1}\!\left(S_{\mathrm{total}} - C(T)\right),
\qquad S_{\mathrm{total}} = \int_0^T S\,\mathrm{d}t
$$

$S_{\mathrm{total}}$ is identical across schedules at matched cumulative exposure,
so the **entire** schedule effect is the terminal residual $C(T)$ — drug delivered
but still in the tissue when the window closes. Three consequences:

1. Earlier delivery realises more exposure (a bolus has washed out by $T$; an
   infusion still running has not), so exposure decreases monotonically with the
   number of fractions and with plasma half-life.
2. The relative effect is bounded by $\tau_{\mathrm{tissue}}/T = 1/(kT)$, and the
   constant-infusion limit saturates that bound.
3. As $T \to \infty$ (window past full washout) the schedule effect vanishes and
   the two modes coincide; they also coincide exactly for a constant plasma level.

Note that the naive transport-limited expectation — a slow infusion penetrating
deeper than a bolus — is *false* for this linear operator: the un-counted term
$M^{-1}C(T)$ is a doubly-smoothed, flatter field, so removing it leaves a
**steeper** profile and the bolus wins on relative reach as well as on total
exposure. A saturable perivascular sink (not in this operator) would be needed to
reverse that.

Note also that the source is a one-way leak, not a Starling gradient
$(C_{\mathrm{plasma}} - C)$ — the tissue never back-fluxes into the vessel. This is
unchanged by time resolution.

---

## 4. Oxygen: steady reaction–diffusion and hypoxia banding

### 4.1 Governing equation

`physics/oxygen.py::solve_oxygen`. Oxygen equilibrates on a timescale of seconds,
far faster than drug transport or cell kill, so it is solved at **steady state**:

$$
D_{O_2}\,\nabla^2 O
\;+\; s\,\rho_{\mathrm{vasc}}(x,y)\,\bigl(1 - O\bigr)
\;-\; c\,\rho_{\mathrm{tum}}(x,y)\,O
\;=\; 0
$$

| symbol | meaning | unit | default |
|---|---|---|---|
| $O$ | relative oxygen saturation | dimensionless, [0,1] | — |
| $D_{O_2}$ | oxygen diffusivity | µm²/s | 2000 |
| $s$ | vascular supply rate constant | 1/s | 0.05 |
| $c$ | per-tumor-cell consumption rate | 1/s per cell | 0.02 |
| $\rho_{\mathrm{vasc}}$ | vascular density field, clipped to [0,1] | — | from reconstruction |
| $\rho_{\mathrm{tum}}$ | tumor-cell **count** binned per grid cell | cells | `_tumor_density` |

The source is written $s\rho_{\mathrm{vasc}}(1-O)$ rather than $s\rho_{\mathrm{vasc}}$
so that supply switches off as tissue saturates ($O\to1$) — this is what bounds the
solution in $[0,1]$ and makes the operator positive definite.

### 4.2 Discretisation and the linear solve

`_build_laplacian` builds the 5-point finite-difference Laplacian with **no-flux
(Neumann)** boundaries imposed by ghost cells equal to the interior cell — i.e. a
missing neighbour simply drops out of the stencil, so the centre coefficient is
$-(\text{number of existing neighbours})/\Delta x^2$:

$$
(\mathbf{L}O)_{i,j} = \frac{1}{\Delta x^{2}}\sum_{(p,q)\in\mathcal{N}(i,j)}\bigl(O_{p,q} - O_{i,j}\bigr)
$$

Rearranging the PDE gives the sparse system

$$
\underbrace{\Bigl[-D_{O_2}\mathbf{L} + \mathrm{diag}\bigl(s\rho_{\mathrm{vasc}} + c\,\rho_{\mathrm{tum}}\bigr)\Bigr]}_{\mathbf{A}}\,O \;=\; \underbrace{s\,\rho_{\mathrm{vasc}}}_{\mathbf{b}}
$$

$-\mathbf{L}$ is symmetric positive semi-definite and the diagonal reaction term is
$\ge 0$ and strictly positive wherever there is vasculature or tumor, so $\mathbf{A}$
is **symmetric positive definite** and a single sparse direct solve (`spsolve`) gives
the exact steady state — no iteration, no CFL, fully deterministic (asserted:
solving twice returns a bit-identical field). The result is clipped to $[0,1]$. A
guard returns zeros if both source and reaction vanish everywhere (singular system).

**Note (differs from §3):** this Laplacian uses the full rectangle with Neumann edges
and does **not** apply the tissue mask; masking is handled by the fact that
$\rho_{\mathrm{vasc}}$ and $\rho_{\mathrm{tum}}$ vanish outside tissue.

**Resolution dependence — FIXED.** $\rho_{\mathrm{tum}}$ used to be a *raw count*
per grid cell (`np.add.at` on binned cell indices) rather than a density per unit
area. Counts scale as $\Delta x^{2}$, so the sink strength — and therefore the
hypoxic fraction and the downstream resistance modulation — was **grid-spacing
dependent**. On the real Xenium reconstruction the hypoxic tumor fraction ran
$0.937 / 0.364 / 0.134$ at $\Delta x = 30 / 20 / 12\ \mu m$; on a controlled
synthetic tissue (identical physical extent and cells) it ran $0.998 \to 0.000$
over $\Delta x = 30 \to 10\ \mu m$.

`_tumor_density` now returns an **areal density**,

$$
\rho_{\mathrm{tum}}(x,y) = \frac{N_{\mathrm{tum}}(x,y)}{\Delta x^{2}}\; A_{\mathrm{ref}},
\qquad A_{\mathrm{ref}} = 144\ \mu\mathrm{m}^2 = (12\ \mu m)^2,
$$

i.e. tumor cells per fixed reference tile — dimensionless, of order 1, and free of
the $\Delta x^{2}$ factor. $A_{\mathrm{ref}}$ is the converged production grid
spacing squared (so the correction is a pure re-parameterisation at
$\Delta x = 12\ \mu m$ and changes nothing there) and is also close to a single
glioma-cell footprint ($\sim 100\ \mu\mathrm{m}^2$). On the synthetic tissue the
mean oxygen is now invariant to $10^{-4}$ across a 4× range of $\Delta x$
(asserted in `_test_...` / `tests/test_physics.py::test_oxygen_sink_is_grid_independent`).

**The *tissue mask* was also grid-dependent — now FIXED (geometric mask).**
`reconstruct/real_xenium.py` used to define `Domain.tissue_mask` as "grid cells
containing at least one segmented cell", which covers $50.6\% / 32.7\% / 13.9\%$ of
the rectangle at $\Delta x = 30/20/12\ \mu m$. `real_fields._kde_on_grid` multiplies
the vascular field by that mask, so the total O₂ *supply* shrank as the grid was
refined. The two errors previously **partially cancelled**; with the sink fixed but
the mask still grid-dependent, mean O₂ at the (fixed) real cell positions ran
$0.522 / 0.426 / 0.251$ across those three grids.

`reconstruct/real_xenium.py::tissue_mask_from_cells` now builds the mask
**geometrically** and it is the production default: a grid cell is tissue when its
centre lies within $R = 30\ \mu m$ (`TISSUE_DISC_RADIUS_UM`) of at least one measured
centroid, i.e. the mask is the union of fixed *physical* discs and is a property of
the measured cell cloud rather than of the mesh. $R$ is a declared constant justified
against this section's own geometry — nearest-neighbour centroid spacing is median
$16.7\ \mu m$ (p10 10.4, p90 24.2, p99 34.2), so $R \approx 1.8\times$ the median and
just above the p90 bridges $>90\%$ of real inter-cell gaps, giving a *connected*
domain (an isolated tissue pixel cannot exchange drug or O₂ with anything) while
extending only ~1 cell diameter past the outermost cell. Measured result:

| $\Delta x$ (µm) | coverage (legacy / geometric) | mean O₂ at cells (legacy / geometric) |
|---:|:--|:--|
| 30 | 0.506 / 0.672 | 0.522 / 0.570 |
| 20 | 0.327 / 0.674 | 0.426 / 0.565 |
| 12 | 0.139 / 0.677 | 0.251 / 0.558 |
| **spread** | **3.64× / 1.007×** | **2.08× / 1.02×** |

The legacy rule remains reachable as `mode='occupancy'` for comparison only, and
`tests/test_pipeline.py::test_geometric_tissue_mask_coverage_is_dx_independent`
pins both behaviours.

**Honest caveat on the absolute level — now the dominant uncertainty.**
`OxygenConfig.consumption` was never calibrated against a measured oxygen tension,
so $c\,A_{\mathrm{ref}}$ is a single un-calibrated lumped constant and the *absolute*
hypoxic fraction is very sensitive to it. With the mask fixed, mean O₂ on the real
section rises to $\approx 0.56$ and **the predicted hypoxic tumor fraction is
$0.0$ at every $\Delta x$** — which, for a tumor type partly defined by
pseudopalisading necrosis, is no more credible than the $0.94$ the coarse grid used
to report. Both the old ladder ($0.94 \to 0.07$ over $\Delta x = 30 \to 8\ \mu m$)
and the new $0.0$ are un-calibrated model output. What the two fixes establish is
that the answer no longer depends on the numerical grid — **not** that its level is
calibrated. Nothing downstream should rest on the absolute hypoxic fraction until
`consumption` is anchored to a measurement.

### 4.3 Oxygen → cell-state banding

`cells/agents.py::update_cell_states`. Each cell samples $O$ at its grid index
(`Domain.to_index`), and with $t = $ `OxygenConfig.hypoxia_threshold` (default 0.2):

$$
\text{state} =
\begin{cases}
\textsf{HYPOXIC}, & O < t \\
\textsf{QUIESCENT}, & t \le O < 2t \\
\textsf{VIABLE}, & O \ge 2t
\end{cases}
$$

`DEAD` cells are never relabelled (no resurrection). The update is fully vectorised
and deterministic — no randomness. `physics/oxygen.py::label_hypoxia` is a reduced
variant that only sets the HYPOXIC label.

**Status:** the threshold $t=0.2$ and the "$[t, 2t)$ quiescent shoulder" are
**declared modelling choices**, documented as such in the module — an at-risk
perfusion band just above the hypoxia line, not a measured physiological cutoff.
`state_fractions` reports the per-state fractions of tumor cells (summing to 1).

---

## 5. Interstitial fluid pressure and the Jain effect

`physics/pressure.py::solve_ifp`. Poroelastic-lite, steady state.

### 5.1 Governing equation

Combining **Darcy's law** for interstitial flow with distributed Starling
source/sink terms and linearising about the interstitial pressure gives a
screened-Poisson (Helmholtz-type) equation for $p$ [mmHg]:

$$
-\nabla\!\cdot\!\bigl(K\,\nabla p\bigr) \;+\; (a + d)\,p \;=\; a\,p_v
$$

with the Darcy velocity

$$
\mathbf{v} = -K\,\nabla p \qquad [\mu\mathrm{m/s}]
$$

| symbol | meaning | unit | model |
|---|---|---|---|
| $K$ | interstitial hydraulic conductivity | µm²/(mmHg·s) | $K_0(1 - \kappa_{\mathrm{ecm}}\rho_{\mathrm{ecm}})$, $K_0 = 1.2$, $\kappa_{\mathrm{ecm}} = 0.6$ |
| $a$ | vessel-wall filtration conductance | 1/(mmHg·s) | $\rho_{\mathrm{vasc}}\,L_p^{\mathrm{eff}}$ |
| $L_p^{\mathrm{eff}}$ | effective wall conductivity | 1/(mmHg·s) | $L_p^{0}(1 + G\,\varphi)$, $L_p^{0} = 10^{-3}$, $G = 12$ |
| $d$ | lymphatic drainage conductance | 1/(mmHg·s) | $\lambda_L(1 - \varphi)$, $\lambda_L = 3\times10^{-2}$ |
| $p_v$ | effective microvascular driving pressure | mmHg | 25 |
| $\varphi$ | smoothed tumor occupancy, $[0,1]$ | — | `_tumor_field` |

$\varphi$ is built by binning tumor-cell counts, Gaussian-smoothing with
$\sigma = 5$ grid cells, and normalising by its maximum (so, unlike
$\rho_{\mathrm{tum}}$ in §4, it *is* resolution-robust).

The two tumor pathologies enter through $\varphi$: **leaky vessels** raise $a$ by up
to $1{+}G = 13\times$, and **absent lymphatics** drive $d \to 0$.

### 5.2 Boundary conditions

**Dirichlet $p = 0$ mmHg** on the domain border **and** on every non-tissue node
(normal-tissue reference). Interior closure is provided by the screening term
$(a+d)$; where $a = d = 0$ the equation degenerates to Laplace and no-flux is
implicit in the assembly (a face to a Dirichlet neighbour acts as a $p=0$ ghost).
After the solve $p$ is clipped at $0$ (justified by the maximum principle: with
$b \ge 0$ and an M-matrix $\mathbf{A}$, $p \ge 0$ up to round-off).

### 5.3 Discretisation

5-point finite-volume stencil with **harmonic-mean face conductances**
$K_f = 2K_PK_N/(K_P + K_N + 10^{-30})$ divided by $\Delta x^2$ (same series-resistance
argument as §3.3). For each interior node the face conductance is added to the
diagonal *always* — a Dirichlet neighbour contributes as a zero-valued ghost — while
the off-diagonal $-K_f/\Delta x^2$ is added **only** when the neighbour is itself an
unknown. Dirichlet rows are identity rows with $b = 0$. The reaction term
$(a+d)$ is added to interior diagonals; $b = a\,p_v$ on interior nodes. Solved with
`scipy.sparse.linalg.spsolve`.

The velocity is then evaluated with **central differences** (`np.gradient`, spacing
$\Delta x$) on the pressure field:

$$
v_x = -K\,\frac{\partial p}{\partial x},\qquad v_y = -K\,\frac{\partial p}{\partial y}.
$$

### 5.4 Why elevated flat core IFP ⇒ zero core flow and outward margin flow

Away from strong gradients, the equation is dominated by its algebraic (reaction)
balance, whose local equilibrium is

$$
p_\infty(x,y) = \frac{a}{a + d}\,p_v .
$$

* **Tumor core:** $\varphi \to 1$, so $d = \lambda_L(1-\varphi) \to 0$ while $a$ is
  amplified by the leak gain. Hence $p_\infty \to p_v$ — a **high plateau**. Because
  $a/(a+d)$ is spatially near-constant across the core, $\nabla p \approx 0$ and
  therefore $\mathbf{v} = -K\nabla p \approx \mathbf{0}$: **there is essentially no
  convective transport inside the tumor.**
* **Normal tissue:** $\varphi \to 0$, so $d \gg a$ and $p_\infty \to 0$; the Dirichlet
  border pins the far field at 0 mmHg.
* **Margin:** the transition of $\varphi$ from 1 to 0 happens over a thin rim, so the
  entire $p_v \to 0$ drop is compressed there. The diffusive term $-\nabla\cdot(K\nabla p)$
  smooths this over a screening length $\ell = \sqrt{K/(a+d)}$, but the gradient
  remains steep, giving a large $\mathbf{v} = -K\nabla p$ directed **radially
  outward** (down the pressure gradient).

That outward margin flow is the **Jain effect** (Jain 1987; Baxter & Jain 1989): it
physically opposes inward transport of a blood-borne drug, washing molecules back out
of the rim and starving the core. In §3 it enters through the $-\nabla\cdot(\mathbf{v}C)$
term. The module's self-test asserts all three signatures: core IFP > margin IFP,
mean radial margin velocity $> 0$ (outward), and core speed < margin speed.

**Status: order-of-magnitude parameters, not patient-calibrated.** $K_0$, $p_v$,
$L_p^0$, $G$, $\lambda_L$ are literature-scale values chosen so the *qualitative*
Jain effect emerges robustly; the module says so explicitly. Practically, at these
parameter values interstitial speeds are $< 0.7$ µm/s and change the penetration
index by $< 3\times10^{-4}$ (§10), so in this model IFP is a mechanistically correct
but numerically minor contributor.

---

## 6. Pharmacodynamics: exposure → kill

`cells/pkpd.py::apply_drug_kill`.

### 6.1 Exposure reduction

Each cell samples the AUC field at its own grid index. The AUC (µM·s) is reduced to a
**time-averaged effective concentration** over the exposure window:

$$
E = \frac{\mathrm{AUC}}{T} \quad [\mu M],
\qquad T = \texttt{transport.times[-1]}\ [\mathrm{s}]
$$

(with $T \to 1$ as a guard if the window is empty). Working in a *concentration*
rather than an exposure is what allows the drug's IC₅₀ to be used directly as the
half-effect point: since $E$ and IC₅₀ are both µM,

$$
E_{50} \equiv \mathrm{IC}_{50}
$$

with **no unit conversion** — this is exactly how IC₅₀ is defined pharmacologically.
The per-cell AUC is also stored on `cells.exposure`.

### 6.2 Hill dose–response and modulation

$$
p_{\mathrm{Hill}} = R_{\max}\,\frac{E^{h}}{E_{50}^{h} + E^{h}}
$$

| symbol | meaning | unit | source |
|---|---|---|---|
| $R_{\max}$ | asymptotic kill ceiling as $E\to\infty$ | — | `PDConfig.max_kill_rate` (default 1.0) |
| $h$ | Hill coefficient (curve steepness) | — | `DrugParams.hill` |
| $E_{50}$ | half-effect concentration | µM | `DrugParams.ic50` |

Multiplicative modulation by intrinsic and microenvironmental resistance:

$$
p_{\mathrm{kill}} = \mathrm{clip}\Bigl(p_{\mathrm{Hill}}\;\cdot\;\underbrace{\sigma_c}_{\text{sensitivity}}\;\cdot\;R_{\mathrm{state}},\;0,\;1\Bigr)
$$

$$
R_{\mathrm{state}} =
\begin{cases}
R_{\mathrm{hyp}} = 0.35, & \text{cell is HYPOXIC} \\
R_{\mathrm{qui}} = 0.50, & \text{cell is QUIESCENT} \\
1, & \text{otherwise}
\end{cases}
$$

$\sigma_c \in [0,1]$ is the per-cell intrinsic sensitivity (MGMT status, resistance
clones, …); $\sigma_c = 0$ means the cell is never killed. Both resistance factors
are $<1$, so hypoxic/quiescent cells survive more often at identical exposure — the
classic drug resistance of poorly perfused, non-cycling regions. Already-DEAD cells
are forced to $p_{\mathrm{kill}} = 0$.

**Status: $\sigma_c$ is now an explicit, deterministic, drug-specific assumption —
never random.** It used to be `rng.uniform(0,1)` in the demo paths and a min–max
scaled TUMOR-marker (EGFR/SOX2/OLIG2/MKI67) z-score in the real Xenium
reconstruction — a lineage signal with no mechanistic link to chemosensitivity that
nonetheless scaled every absolute kill fraction. Current behaviour
(`reconstruct/real_xenium.py::per_cell_sensitivity`, table
`SENSITIVITY_DETERMINANTS`):

| drug | determinant on the 319-gene panel | status |
|---|---|---|
| doxorubicin, etoposide, mitoxantrone | **TOP2A** | the drug's actual **molecular target** (TOP2 poisons; PMID 19377506, 18574145) |
| cisplatin, lomustine, carmustine, methotrexate, vincristine | **MKI67** | **mechanistic proxy** — proliferation-gated cytotoxicity, *not* a measured sensitivity |
| temozolomide, procarbazine | *none* (determinant is **MGMT**, absent from the panel) | **unmodelled**: $\sigma_c \equiv 1$ |
| bevacizumab | *none* (no direct tumour-cell cytotoxicity) | $\sigma_c \equiv 1$ |

Where a determinant exists, the per-cell log1p z-score $z_c$ is squashed into a
bounded band,

$$
\sigma_c = \sigma_{\min} + (1-\sigma_{\min})\,\frac{1}{1+e^{-\kappa z_c}},
\qquad \sigma_{\min}=0.2,\ \kappa=1,
$$

so transcript dropout can never fabricate a perfectly resistant ($\sigma_c=0$) cell.
The reconstruction's *default* is uniform $\sigma_c = 1$ (drug-agnostic geometry, no
invented heterogeneity); the drug-specific field is applied at simulation time by
`reconstruct/real_pipeline.py::apply_drug_sensitivity`, leaving `core`, `sim/engine`
and `physics` untouched.

**Honest limitation.** On this *targeted* panel the determinant transcripts are very
sparse — TOP2A is detected in 1.08 % of cells, MKI67 in 0.35 % — so the derived field
is near-binary: ~99 % of cells sit at $\sigma_c \approx 0.59$ (the arbitrary
$\mathrm{sigmoid}(0)\!=\!0.5$ midpoint of the mapping) with a thin tail at 1.0. It
therefore acts mostly as a near-uniform *scale* on kill plus a small genuinely-high
tail, not as a rich heterogeneity map. **Every absolute kill fraction is conditional
on this assumption**, and none of it is calibrated against a measured drug response.
The resistance multipliers $R_{\mathrm{hyp}}, R_{\mathrm{qui}}$ remain declared
constants, not fitted.

### 6.3 The Bernoulli kill draw

For each cell, one uniform draw $u_c \sim \mathcal{U}(0,1)$ from the supplied
`np.random.Generator`:

$$
\text{cell } c \text{ dies} \iff \text{alive}_c \;\wedge\; u_c < p_{\mathrm{kill},c}
$$

Given a seeded generator the outcome is bit-for-bit deterministic. The population
readouts follow:

$$
\text{viable tumor fraction} = \frac{\#\{\text{alive tumor}\}}{\#\{\text{tumor}\}},
\qquad
\text{kill fraction} = 1 - \text{viable tumor fraction} \;=\; \text{predicted response}
$$

(`sim/engine.py`). Because the draw is a Bernoulli sample rather than the expected
value, the kill fraction carries sampling noise $\sim\sqrt{\bar p(1-\bar p)/n}$; the
numerical audits therefore average it over seeds $\{0,1,2,3\}$ while the deterministic
penetration index is used as the primary signal.

The engine's diagnostic viable-fraction *trajectory*
(`sim/engine.py::_viable_history`) reuses the identical Hill mapping at scaled
exposures $f\cdot\mathrm{AUC}$, $f\in\mathrm{linspace}(0,1,6)$, thresholded against a
**single shared** per-cell draw $u_c$, which makes the curve monotone non-increasing
by construction (since $p_{\mathrm{kill}}$ is increasing in $E$) rather than a set of
independent noisy snapshots.

---

## 7. Penetration length, the binding-site barrier, and anchor conversion

### 7.1 Derivation of $L = \sqrt{D_{\mathrm{eff}}/k}$

`validation/pk_reality.py::penetration_length_um`. Consider tissue occupying $x>0$
beside a vessel wall at $x=0$, at steady state, with no convection and a first-order
loss $k = k_{\mathrm{cl}} + k_b$. The 1-D reaction–diffusion balance is

$$
D_{\mathrm{eff}}\,\frac{\mathrm{d}^2 C}{\mathrm{d}x^{2}} - k\,C = 0 .
$$

The general solution is $C(x) = A e^{-x/L} + B e^{+x/L}$ with
$L = \sqrt{D_{\mathrm{eff}}/k}$. Boundedness as $x\to\infty$ forces $B = 0$; the wall
condition $C(0) = C_0$ gives $A = C_0$. Hence

$$
\boxed{\;C(x) = C_0\,e^{-x/L}, \qquad L = \sqrt{\frac{D_{\mathrm{eff}}}{k_{\mathrm{cl}} + k_b}}\;}
\qquad [\mu\mathrm{m}]
$$

$L$ is the **1/e (e-folding) distance**: the depth at which the concentration falls
to $C_0/e \approx 0.368\,C_0$. It is the *principled, grid-independent* penetration
metric — unlike the penetration index $\Pi$ of §3.8, it is a length with a direct
experimental counterpart.

The implementation mirrors the PDE's tortuosity exactly so the analytic number
describes what the solver experiences:

$$
D_{\mathrm{eff}} = \max\bigl(D\,(1 - \alpha_{\mathrm{ecm}}\,\rho_{\mathrm{ecm}}),\,10^{-6}\bigr),
\qquad \alpha_{\mathrm{ecm}} = 0.6,\ \ \rho_{\mathrm{ecm}} = 0.4 \ \text{(default)}
$$

and returns $L = \infty$ when $k = 0$ (a drug that is never consumed penetrates
everywhere).

**Assumptions:** semi-infinite 1-D slab, steady state, constant coefficients, no
advection, linear (unsaturable) first-order loss, single planar source. Real tissue
has curved vessels at finite spacing, so $L$ should be read as the *boundary-layer
scale*, not an exact profile.

### 7.2 The binding-site barrier and the Thiele modulus

`molecular/binding_site_barrier.py`. A low diffusivity alone only *slows* an
antibody — given enough time it would still fill the tissue. That is not what is
observed: anti-tumor antibodies pile up in a thin perivascular cuff and leave an
antibody-free core (Fujimori/Weinstein 1990–92). The mechanism is **saturable,
high-affinity antigen binding**, which is a *sink*, not a *drag*.

**Pseudo-first-order antigen capture:**

$$
k_{\mathrm{fwd}} = k_{\mathrm{on}}\,[\mathrm{Ag}]\quad[1/\mathrm{s}],
\qquad
\text{irrev} = \frac{k_{\mathrm{fwd}}}{k_{\mathrm{fwd}} + k_{\mathrm{off}}}\in[0,1],
\qquad
k_{\mathrm{ret}} = k_{\mathrm{fwd}}\cdot\text{irrev}
$$

| symbol | meaning | unit | default (Thurber/Wittrup lineage) |
|---|---|---|---|
| $[\mathrm{Ag}]$ | local perivascular antigen concentration | mol/L | $3.0\times10^{-7}$ (~0.3 µM) |
| $k_{\mathrm{on}}$ | association rate | 1/(M·s) | $1.0\times10^{5}$ |
| $k_{\mathrm{off}}$ | dissociation rate | 1/s | $1.0\times10^{-4}$ ($K_d = k_{\mathrm{off}}/k_{\mathrm{on}} = 1$ nM) |
| $k_{\mathrm{base}}$ | generic barrier-free loss | 1/s | $5\times10^{-4}$ |
| $R$ | inter-vessel half-distance (Thiele length scale) | µm | 150 |

The `irrev` factor is the **affinity gate**: a high-affinity binder (small
$k_{\mathrm{off}}$) stays bound once captured, so binding is a near-perfect sink
(irrev → 1); a low-affinity binder releases and "hops" forward (irrev → 0). This is
why *raising* $k_{\mathrm{off}}$ is the textbook way to relieve the barrier — and it
is the mechanistic reason high-affinity binders penetrate *worse*.

**Resulting penetration length and Thiele modulus:**

$$
L = \sqrt{\frac{D}{k_{\mathrm{ret}} + k_{\mathrm{base}}}},
\qquad
\phi = \frac{R}{L} = R\sqrt{\frac{k_{\mathrm{ret}} + k_{\mathrm{base}}}{D}}
$$

$\phi$ is the classic Thiele modulus / Damköhler-style criterion (Thurber, Schmidt &
Wittrup 2007–08): $\phi \gg 1$ means binding consumes the antibody long before it
crosses the inter-vessel distance → **barrier-limited, perivascular-only**
distribution; $\phi \lesssim 1$ means the drug saturates the tissue. The code flags

$$
\texttt{is\_barrier\_limited} = \bigl(k_{\mathrm{ret}} > k_{\mathrm{base}}\bigr) \wedge \bigl(\phi > 1\bigr).
$$

**Who experiences it.** Antibodies always; small molecules only if the caller
explicitly passes `high_affinity=True`. For a plain cytotoxic small molecule
$k_{\mathrm{ret}} = 0$ exactly, so its $L$ stays diffusion/clearance-limited
(hundreds of µm) and it is returned unchanged. `apply_barrier_to_drugparams` folds
$k_{\mathrm{ret}}$ into the engine's existing sink via
`binding ← binding + k_ret` — no engine modification is needed, because the PDE
already integrates $k = k_{\mathrm{cl}} + k_b$.

**Status: lumped effective-medium model.** This is a single perivascular sink term,
**not** a spatially resolved moving-front / shrinking-core solve of saturable binding.
It captures the correct mechanism, the correct scaling
($L \sim \sqrt{D/k_{\mathrm{ret}}}$, affinity- and antigen-density-dependent) and the
correct rank (antibody ≪ small molecule); it does not capture the time-dependent
saturation front. The default $[\mathrm{Ag}]$, $k_{\mathrm{on}}$, $k_{\mathrm{off}}$
are literature-anchored order-of-magnitude values, not per-target measurements.

### 7.3 Half-distance ↔ 1/e distance conversion

`validation/pk_anchors.py::half_distance_to_efold_length`. Published perivascular
immunofluorescence studies report $d_{1/2}$ — the distance at which drug intensity
falls to **half**. The model's $L$ is the **1/e** distance. For the exponential
profile $C(x) = C_0 e^{-x/L}$:

$$
\tfrac12 = e^{-d_{1/2}/L}
\;\Longrightarrow\;
d_{1/2} = L\ln 2 \approx 0.693\,L
\;\Longrightarrow\;
\boxed{\,L = \frac{d_{1/2}}{\ln 2} \approx 1.443\,d_{1/2}\,}
$$

Failing to apply this conversion introduces a systematic 44% bias. Concretely, the
Primeau doxorubicin anchor $d_{1/2} = 40$–50 µm corresponds to an L-equivalent band
of 58–72 µm. **Both conventions appear in the results tree and give different fitted
constants** — see §8.3.

---

## 8. Calibration: one mechanistic constant

`validation/calibration.py`.

### 8.1 The diagnosis

The structure-only priors of §2.8 give $k = k_{\mathrm{cl}} + k_b \approx 3\times10^{-4}$/s
(a ~1-hour timescale), hence $L \approx 700$ µm — a ~15× over-prediction against the
published ~40–50 µm doxorubicin half-distance. Reproducing a ~45 µm boundary layer
requires

$$
k = \frac{D_{\mathrm{eff}}}{L^{2}} \approx 7\times10^{-2}\ \mathrm{s^{-1}}
$$

(a ~15-second timescale), i.e. **~300× faster**. The missing physics is not slow
metabolic clearance; it is **rapid cellular uptake plus high-capacity intracellular /
DNA binding**, a near-irreversible perivascular sink.

### 8.2 The single-constant model and its closed-form solve

Keep the engine's first-order sink but split it into a slow passive part and a fast,
molecularly-scaled uptake part:

$$
k_{\mathrm{eff}}(\text{drug}) \;=\; k_{\mathrm{passive}} \;+\; \alpha\; u(\text{drug})
$$

| symbol | meaning | unit | implementation |
|---|---|---|---|
| $k_{\mathrm{passive}}$ | slow metabolic washout, **not** amplified | 1/s | `k_passive()` = `drug.clearance` (~$10^{-4}$/s) |
| $u$ | uptake proxy: structure-based tissue/DNA binding-capacity rate | 1/s | `uptake_proxy()` = `drug.binding` (§2.8) |
| $\alpha$ | the **single** free scaling constant | dimensionless | `calibrate_alpha()` |

Because $L = \sqrt{D_{\mathrm{eff}}/k_{\mathrm{eff}}}$ and $k_{\mathrm{eff}}$ is
**linear in $\alpha$**, the calibration is one closed-form inversion — no optimiser,
no iteration:

$$
k_{\mathrm{eff}}^{*} = \frac{D_{\mathrm{eff}}}{L_{*}^{2}}
\qquad\Longrightarrow\qquad
\boxed{\;\alpha = \frac{k_{\mathrm{eff}}^{*} - k_{\mathrm{passive}}}{u}
= \frac{D_{\mathrm{eff}}/L_{*}^{2} - k_{\mathrm{cl}}}{k_b}\;}
$$

evaluated on freshly-compiled (un-calibrated) doxorubicin `DrugParams` at
$\rho_{\mathrm{ecm}} = 0.4$, with $D_{\mathrm{eff}}$ computed by the *same*
tortuosity expression the PDE uses (`_d_eff_tortuous`) so the fit targets exactly the
metric it is compared against. The code raises if $u \le 0$ or if the required
$k_{\mathrm{eff}}^{*}$ falls below the passive floor (target unreachable).

`calibrated_drugparams` then writes the result back into a normal `core.DrugParams`:

$$
\texttt{binding} \leftarrow \alpha\,u, \qquad \texttt{clearance} \leftarrow k_{\mathrm{passive}},
\qquad\text{so}\qquad \texttt{clearance} + \texttt{binding} = k_{\mathrm{eff}} .
$$

Diffusivity, MW, IC₅₀, Hill and BBB are untouched.

### 8.3 Fitted values and the out-of-sample test

| anchor convention | target $L_*$ | fitted $\alpha$ | status |
|---|---:|---:|---|
| raw half-distance treated as $L$ | 45 µm | 309.6 | **retired** (dimensionally wrong) |
| e-fold-converted, $L = d_{1/2}/\ln 2$ | 64.9 µm | **148.6** | **the one convention in the tree** |

**RESOLVED.** Both values used to be live for the same constant — `calibrate_alpha`
defaulted to the raw 45 µm half-distance (giving 309.6) while `pk_multidrug.py`
passed the converted target (giving 148.6); they differ by exactly $1/\ln^2 2 = 2.081$.
`calibrate_alpha` now takes its target as an **e-folding length**
(`target_efold_um`, default `pk_anchors.half_distance_to_efold_length(45) = 64.92 µm`),
so every caller in the tree gets **α = 148.6**, and `revalidate` compares
out-of-sample predictions to the *converted* anchor band (36–108 µm) rather than the
raw half-distances. Guarded by
`tests/test_reality_layer.py::test_calibrate_alpha_target_is_an_efold_length`. Any
figure or results file still quoting α = 309.6 predates this correction.

With $\alpha$ **frozen**, every other drug's $L$ is a pure prediction:

* **Out-of-sample (α = 148.6, current):** mitoxantrone $L = 48$ µm and topotecan
  $L = 84$ µm vs. the L-equivalent band 36–108 µm → 2/2 pass. Temozolomide
  penetrates far further than doxorubicin — correct direction for the opposite
  (weakly-bound) regime.
* *(Historical, α = 309.6: mitoxantrone $L = 33.4$ µm compared against the raw
  25–75 µm half-distance band. Both sides of that comparison used the wrong
  convention, which is why it also "passed".)*

**Honest caveats** (carried verbatim from the module):
1. **One anchor sets the absolute scale** — every drug's absolute penetration rides
   on a single doxorubicin measurement whose own uncertainty is plausibly ±50–100%.
2. **Mitoxantrone is chemically related to doxorubicin** (both anthracenedione
   intercalators). The measurement is independent, but this is a *magnitude-transfer*
   test within the strongly-bound class, not a maximally independent one. The
   temozolomide contrast probes the opposite regime, and only qualitatively.
3. **The proxy is an ordering, not a rate.** $u = k_b$ is a low-confidence logP/HBD
   heuristic; $\alpha$ corrects its *scale* but inherits its *ranking* errors. A drug
   the heuristic mis-ranks (active transport, efflux, ion-trapping) is mis-predicted.
4. **Nothing was refit** to mitoxantrone, topotecan or temozolomide.

---

## 9. Reconstruction mathematics

The tissue fields the PDEs consume are not painted on — they are estimated from
measured spatial transcriptomics. Three estimators are involved.

### 9.1 Kernel density estimation of the tissue fields

`reconstruct/real_fields.py::_kde_on_grid`. Given $n$ measured cell positions
$(x_c, y_c)$ in µm with optional weights $w_c$, the field is a Gaussian KDE evaluated
on the grid, computed as **bin-then-smooth** (exact for a grid-supported estimate and
$O(N)$ rather than $O(n\,N)$):

$$
n_{ij} = \sum_{c} w_c\,\mathbb{1}\bigl[(\mathrm{row}_c,\mathrm{col}_c) = (i,j)\bigr],
\qquad
\tilde f = G_\sigma * n,
\qquad
\sigma = \frac{b}{\Delta x}\ \text{[grid cells]}
$$

$$
f_{ij} = \mathrm{clip}\!\left(\frac{\tilde f_{ij}\,\mathbb{1}[\texttt{tissue}_{ij}]}{\max_{kl}\bigl(\tilde f_{kl}\,\mathbb{1}[\texttt{tissue}_{kl}]\bigr)},\,0,\,1\right)
$$

with $G_\sigma$ an isotropic Gaussian (`scipy.ndimage.gaussian_filter`, `mode="nearest"`
— i.e. edge-replicating, so the estimate is not artificially damped at the domain
boundary), and bandwidth $b = 40$ µm by default (endothelium marks vessel walls a few
cell diameters apart; 40 µm smooths discrete cell hits into a continuous field
without erasing real heterogeneity). Fully deterministic — no RNG.

**The max-normalisation is a modelling choice with consequences:** it makes each field
a *relative* density in $[0,1]$, so "vascular density 1.0" means "the densest vessel
region in this section", not an absolute perfused-volume fraction. Absolute rates
($s$, $c$, $L_p^0$, $\beta_{\mathrm{BBB}}$) are therefore only interpretable relative
to that per-section normalisation.

**Which fields are real, which are proxies (honest accounting):**

| field | estimator | status |
|---|---|---|
| $\rho_{\mathrm{vasc}}$ | KDE of measured `ENDOTHELIAL` (PECAM1⁺/FLT1⁺) cell positions | **REAL** — vessel locations are measured cells, no geometry invented |
| $\rho_{\mathrm{ecm}}$ | KDE of measured `STROMAL` (CAF-like) cell positions | **PROXY** — ECM is acellular and cannot be binned; matrix-depositing cell density stands in for matrix. `ecm_density_from_positions(weights=…)` exists precisely so a caller can pass ECM-transcript abundance (COL1A2, COL12A1, POSTN, VCAN, DCN, FBLN1, LAMA2 are on the panel) for a higher-fidelity field |
| $O$ | `physics.oxygen.solve_oxygen` on the above | **SOLVED** — real-perfusion-driven, not painted |

### 9.2 NNLS spot deconvolution

`reconstruct/deconvolve.py`. Each Visium spot's expression is modelled as a
non-negative mixture of reference cell-type signatures.

**Signature construction** (`build_signatures`): the Darmanis GBM atlas (GSE84465)
labels are collapsed onto the four coarse `CellType` classes
(Neoplastic→TUMOR; Immune→IMMUNE; Vascular→ENDOTHELIAL;
Astrocyte/Oligodendrocyte/OPC/Neuron→STROMAL). Restricted to the $G$ genes shared
between the Visium panel and the reference, each single cell $m$ is library-size
normalised over those genes,

$$
\tilde{r}_{gm} = \frac{r_{gm}}{\sum_{g'} r_{g'm}},
$$

and averaged within class $k$:

$$
\mathbf{S}_{gk} = \frac{1}{|\mathcal{C}_k|}\sum_{m\in\mathcal{C}_k} \tilde{r}_{gm}
\qquad \mathbf{S}\in\mathbb{R}^{G\times K}_{\ge0}
$$

Normalising *before* averaging puts all columns on a comparable scale, which is what
lets the fitted weights be read as fractions. Genes absent from the reference get an
all-zero row (they carry no discriminative information).

**Per-spot solve** (`deconvolve_spot`): the spot vector is library-size normalised to
$\sum_g b_g = 1$ (scale-invariant for the final fractions, but it keeps the problem
well-conditioned), then

$$
\hat{\mathbf{w}} = \arg\min_{\mathbf{w}\ \ge\ \mathbf{0}} \;\bigl\lVert \mathbf{S}\mathbf{w} - \mathbf{b} \bigr\rVert_2^{2},
\qquad
\mathbf{f} = \frac{\hat{\mathbf{w}}}{\sum_k \hat w_k}
$$

solved with `scipy.optimize.nnls` (Lawson–Hanson active-set). The non-negativity
constraint is the entire regularisation: it is what makes the solution
interpretable as an unmixing rather than an arbitrary least-squares fit, and — since
the signature columns are non-negative and generally not collinear — it typically
yields a sparse, well-determined $\hat{\mathbf{w}}$ without any explicit penalty.
A degenerate all-zero solution (e.g. an empty spot) falls back to the uniform vector
$\mathbf{f} = \mathbf{1}/K$.

**Assumptions:** linear additive mixing of library-size-normalised profiles;
platform-invariant signatures (scRNA-seq reference vs. Visium target); no
per-cell-type capture-efficiency correction; the sum-to-one renormalisation absorbs
any global scale mismatch. $\mathbf{f}$ is a *fraction of normalised expression*,
which equals a cell-count fraction only if cell types have equal total RNA content.

### 9.3 Correlation-based niche-region assignment

`reconstruct/region_label.py`. Each location (a Visium spot or an aggregated Xenium
neighbourhood) is assigned to one of four real GBM anatomic zones — `margin`, `core`,
`perivascular`, `necrotic` — by correlating its expression against **measured Ivy GAP
per-structure signatures** (laser-microdissected FPKM, Puchalski et al. 2018).

**Panel selection** (`region_signatures`): genes with peak FPKM $\ge$ `min_fpkm` (5)
across regions are ranked by their **log-space spread across the four region means**,

$$
\mathrm{spread}(g) = \max_r \log(1 + \mathrm{FPKM}_{gr}) - \min_r \log(1 + \mathrm{FPKM}_{gr}),
$$

and the top $n_{\text{genes}}$ (default 300) are taken, unioned with the canonical
literature markers in `REGION_MARKERS`. High spread = region-discriminative.

**Scoring** (`label_regions`, `_spearman_matrix`): both sides are restricted to the
shared genes and $\log(1+\cdot)$-transformed, then **Spearman (rank) correlation** is
computed as Pearson correlation of the within-row ranks:

$$
\mathbf{X}^{r} = \mathrm{rank}(\mathbf{X}),\quad
\mathbf{S}^{r} = \mathrm{rank}(\mathbf{S}),\quad
\rho_{\ell r} = \frac{\bigl(\mathbf{X}^{r}_{\ell} - \overline{\mathbf{X}^{r}_{\ell}}\bigr)\cdot\bigl(\mathbf{S}^{r}_{r} - \overline{\mathbf{S}^{r}_{r}}\bigr)}
{\bigl\lVert \mathbf{X}^{r}_{\ell} - \overline{\mathbf{X}^{r}_{\ell}} \bigr\rVert_2\;\bigl\lVert \mathbf{S}^{r}_{r} - \overline{\mathbf{S}^{r}_{r}} \bigr\rVert_2}
$$

$$
\text{region}(\ell) = \arg\max_{r}\ \rho_{\ell r}
$$

**Why rank correlation:** the signature is bulk FPKM and the local assay is UMI /
targeted counts. Rank correlation is invariant to any *monotone* transformation, so it
is insensitive to the cross-platform scale difference — which a Pearson correlation on
raw values would not be. Locations with no variance over the shared panel yield NaN
and are labelled `"unlabeled"`. At least 3 shared genes are required.

**Rasterisation** (`region_field`): per-location labels become smooth spatial
fractions by binning, Gaussian-smoothing (bandwidth 100 µm ≈ the Visium spot pitch),
and dividing by the smoothed count of **all labelled** locations:

$$
F_r(i,j) = \mathbb{1}[\texttt{tissue}_{ij}]\cdot
\frac{\bigl(G_\sigma * n_r\bigr)_{ij}}{\bigl(G_\sigma * n_{\mathrm{labelled}}\bigr)_{ij} + \varepsilon},
\qquad \varepsilon = 10^{-12}
$$

so $\sum_r F_r \approx 1$ wherever labelled locations have support and $0$ elsewhere.

**Status: hard assignment.** `argmax` discards the runner-up correlation, so a
location that is nearly equally consistent with `core` and `perivascular` is recorded
as fully one of them; the smoothing in `region_field` partially recovers a soft
composition, but the underlying assignment is a winner-take-all decision, and its
gate (e.g. the endothelial typing threshold) is the single most influential knob in
the sensitivity study (§10.4).

---

## 10. Numerical verification summary

Full results: **`results/numerical_rigor.md`** (generated by
`tmesim/validation/numerical_rigor.py`). Geometry: the real 40,887-cell 10x Xenium
GBM section reconstructed by `reconstruct.real_pipeline.reconstruct_real`; drug: the
α-calibrated doxorubicin; advection off (it moves the penetration index by
$<3\times10^{-4}$), so these audits isolate the implicit diffusion–reaction scheme.

### 10.1 Mass conservation

Real irregular tissue mask, $(170,271)$ at $\Delta x = 30$ µm, 23,302 tissue cells,
$k_{\mathrm{cl}} = k_b = 0$:

| test | relative mass error | min value |
|---|---:|---:|
| (A) seeded blob, pure diffusion, 1080 steps | $9.92\times10^{-14}$ | $0$ |
| (B) public solver, injected source | $2.96\times10^{-15}$ | $0$ |

This is the **linear-solver / round-off floor**, confirming the zero-row-sum argument
of §3.7 holds on the measured mask and not merely on a rectangle: the no-flux
treatment at the tissue interface leaks nothing, and the field never goes negative.

### 10.2 Grid convergence — and the honest resolution requirement

| $\Delta x$ (µm) | grid | penetration index | Δrel | kill fraction | Δrel |
|---:|:--:|---:|---:|---:|---:|
| 60 | 86×136 | 0.2850 | — | 0.0588 | — |
| 40 | 128×204 | 0.2859 | 0.0030 | 0.0581 | 0.0107 |
| 30 | 170×271 | 0.2666 | 0.0676 | 0.0596 | 0.0248 |
| 20 | 255×406 | 0.2291 | 0.141 | 0.0734 | 0.232 |
| 15 | 339×541 | 0.1991 | 0.131 | 0.0764 | 0.0411 |
| 12 | 424×675 | 0.1934 | 0.0287 | 0.0757 | 0.0093 |
| 10 | 508×810 | 0.1869 | 0.0335 | 0.0760 | 0.0035 |

* **`kill_fraction` plateaus** at $\approx 0.076$ on this ladder, moving only 0.0035
  relative between the two finest grids. **This bullet used to read "`kill_fraction`
  is grid-converged"; that wording is withdrawn** — see §10.5 of this document, and
  `results/grid_convergence_gbm.md` §0, where the GBM ladder is **non-monotone**
  (observed orders 0.96 and −3.33, both power-law fits flagged DEGENERATE). A plateau
  on one ladder is not convergence. What the table above does support is that the
  prediction is not a *coarse*-mesh artifact.
* **`penetration_index` converges more slowly and is over-estimated on coarse grids.**
  The reason is physical: the calibrated doxorubicin penetration length is
  $L \approx 45$ µm, so any grid with $\Delta x \gtrsim 45$ µm **cannot resolve the
  perivascular boundary layer** and reports an artificially flat (high) index. The
  requirement is roughly $\Delta x \lesssim L/4$.

  $$\textbf{Quote the penetration index only at } \Delta x \lesssim 12\text{–}15\ \mu\mathrm{m}\ (\Pi \approx 0.187).$$

  The 30 µm value used elsewhere for speed sits ~15–25% high and must be read as an
  upper bound.

This is an **end-to-end resolution study**, not a textbook mesh refinement: $\Delta x$
sets both the PDE mesh *and* the re-rasterised reconstruction (tissue mask + KDE
fields), so the Richardson order estimate

$$
f(h) = f_\infty + C h^{p},\qquad
\frac{f(h_1)-f(h_2)}{f(h_2)-f(h_3)} = \frac{h_1^{p}-h_2^{p}}{h_2^{p}-h_3^{p}}
$$

(solved for $p$ by Brent's method in `_observed_order`) returns **n/a** for both
metrics — the re-rasterisation wobble breaks the single-power ansatz. Convergence is
therefore read off the *successive relative changes*, and the report says so rather
than quoting a fabricated order.

### 10.3 Time-step convergence

Fixed $\Delta x = 30$ µm, fixed source, only $\Delta t$ varied:

| $\Delta t$ (s) | steps | penetration index | rel. err vs finest | mean AUC (µM·s) | AUC rel. err |
|---:|---:|---:|---:|---:|---:|
| 80 | 90 | 0.26660 | $2.0\times10^{-6}$ | 491.9 | $5.4\times10^{-3}$ |
| 40 | 180 | 0.26660 | $9.7\times10^{-7}$ | 493.2 | $2.6\times10^{-3}$ |
| 20 | 360 | 0.26660 | $4.5\times10^{-7}$ | 493.9 | $1.2\times10^{-3}$ |
| 10 | 720 | 0.26660 | $1.9\times10^{-7}$ | 494.3 | $5.2\times10^{-4}$ |
| 5 | 1440 | 0.26660 | $6.5\times10^{-8}$ | 494.4 | $1.7\times10^{-4}$ |
| 2.5 | 2880 | 0.26660 | $0$ | 494.5 | $0$ |

No blow-up at any $\Delta t$ across a 32× range — exactly as the unconditional
stability argument of §3.6 predicts. The AUC error halves as $\Delta t$ halves,
confirming the expected **first-order (backward-Euler) temporal accuracy**.
$\Delta t \approx 10$ s is safely inside the converged regime.

### 10.4 Local parameter sensitivity (elasticities)

Each knob is perturbed by $\pm20\%$ about its base value; the reported quantity is the
**log–log central-difference elasticity** (`_elasticity`):

$$
E = \frac{\mathrm{d}\ln(\text{output})}{\mathrm{d}\ln(\text{input})}
\;\approx\; \frac{\ln y_{\mathrm{hi}} - \ln y_{\mathrm{lo}}}{\ln x_{\mathrm{hi}} - \ln x_{\mathrm{lo}}}
$$

(falling back to a normalised finite difference $\frac{\Delta y/y_0}{\Delta x/x_0}$ if
any log is undefined). $|E| > 1$ amplifies, $|E| < 1$ damps; $E$ is dimensionless and
therefore comparable across knobs with different units.

| knob | $-$pert → $+$pert | $\lvert E_{\mathrm{pen}}\rvert$ | $\lvert E_{\mathrm{kill}}\rvert$ |
|---|:--:|---:|---:|
| `endothelial_typing_threshold` | 0.8→1.2 | **0.520** | 0.120 |
| `grid_dx_um` | 24→36 | 0.235 | 0.354 |
| `kde_bandwidth_um` | 32→48 | 0.198 | **0.444** |
| `calibration_alpha` | 248→372 | 0.110 | 0.203 |

* Penetration is most sensitive to the **endothelial typing threshold** — that gate
  decides which cells become `ENDOTHELIAL` and hence where the delivery source (and
  the O₂ source) sits, reshaping the whole perivascular geometry.
* Kill is most sensitive to the **KDE bandwidth** and **grid $\Delta x$** — both
  control how concentrated vs. smeared the delivered dose is, which is what the Hill
  response integrates.
* **α has the smallest *local* elasticity but the largest *prior* uncertainty.** Over
  ±20% it barely moves the output, but it is pinned to a single anchor whose own
  uncertainty is plausibly ±50–100% (and whose value depends on the $\ln 2$ convention
  of §7.3), and it enters $L = \sqrt{D_{\mathrm{eff}}/k}$ with $k \propto \alpha$. It
  is a low-sensitivity, high-uncertainty knob — the two must not be conflated.

### 10.5 Bottom line

The numerics are trustworthy on the real geometry: mass-conservative to round-off,
unconditionally stable and dt-converged. **The grid-convergence claim that used to
sit in this paragraph is WITHDRAWN** — there is no $\Delta x$ at which this pipeline
has been shown to be grid-converged; the GBM ladder is non-monotone and both
power-law fits are flagged DEGENERATE (`results/grid_convergence_gbm.md`), and
`penetration_index` is retired regardless (§ `results/delivery_metrics.md`).

**The limiting uncertainty is not the solver** — it is, in rough order: the
**partially-supported diffusivity law of §2.5** (still the largest known structural
error, and the only one that changes what a *ranking* means: it now recovers about a
quarter of the measured compound-to-compound spread, having previously recovered
essentially none of it), the reconstruction gates
(endothelial typing threshold, vessel bandwidth), the single-anchor calibration α,
the heuristic terms flagged throughout (§2.7, §2.8, §3.2, §6.2), the unresolved
mesh term, and the placeholder per-cell sensitivity. Error bars built from the
earlier version of this list — which omitted §2.1–§2.5 entirely — **cannot contain
the truth for a similar-mass compound series**. That used to be because the
transport spread there was zero by construction; it is now because the spread is
present but systematically about **4× too small**.

---

## Appendix: where each equation lives

| § | equation | module · function |
|---|---|---|
| 2.1 | MW → $r_s$ | `molecular/stokes_einstein.py::hydrodynamic_radius_nm` |
| 2.2 | $D = k_BT/6\pi\eta r$ | `molecular/stokes_einstein.py::free_diffusivity` |
| 2.3 | overdamped Langevin, MSD $=4Dt$ | `molecular/brownian.py::simulate_effective_diffusivity`; `molecular/ecm_mesh.py::brownian_mesh_hindrance` |
| 2.4 | Ogston partition; Renkin/Faxén $K_d$ | `molecular/ecm_mesh.py::ogston_partition`, `_faxen_drag`, `ogston_renkin_hindrance`, `mesh_pore_radius_nm` |
| 2.5 | $D_{\mathrm{eff}} = D_{\mathrm{BD}}\cdot h_{\mathrm{mesh}}$ | `molecular/compile_drug.py::build_profile` |
| 2.6 | $P = K_pD_m/d$ | `molecular/permeability.py::membrane_permeability` |
| 2.7 | CNS-MPO BBB score | `molecular/permeability.py::bbb_penetrance` |
| 2.8 | $f_u$, $k_b$, $k_{\mathrm{cl}}$ heuristics | `molecular/binding.py` |
| 3.1–3.8 | transport PDE, FV operator, IMEX, AUC, $\Pi$ | `physics/diffusion.py::simulate_drug_transport`, `_build_diffusion_operator`, `_advective_divergence`, `_effective_diffusivity` |
| 3.9 | plasma profile, delivery source | `treatment/drugs.py::plasma_profile`, `delivery_source`; `sim/engine.py` |
| 4.1–4.2 | steady oxygen SPD solve | `physics/oxygen.py::solve_oxygen`, `_build_laplacian` |
| 4.3 | O₂ → state banding | `cells/agents.py::update_cell_states`; `physics/oxygen.py::label_hypoxia` |
| 5 | IFP screened-Poisson, Darcy $\mathbf{v}$ | `physics/pressure.py::solve_ifp`, `_tumor_field` |
| 6 | $E = \mathrm{AUC}/T$, Hill kill, Bernoulli | `cells/pkpd.py::apply_drug_kill` |
| 7.1 | $L = \sqrt{D_{\mathrm{eff}}/k}$ | `validation/pk_reality.py::penetration_length_um` |
| 7.2 | $k_{\mathrm{ret}}$, Thiele $\phi = R/L$ | `molecular/binding_site_barrier.py::binding_site_barrier_factor`, `apply_barrier_to_drugparams` |
| 7.3 | $L = d_{1/2}/\ln 2$ | `validation/pk_anchors.py::half_distance_to_efold_length` |
| 8 | $k_{\mathrm{eff}} = k_{\mathrm{passive}} + \alpha u$; closed-form α | `validation/calibration.py::k_eff`, `calibrate_alpha`, `calibrated_drugparams`, `revalidate` |
| 9.1 | Gaussian KDE fields | `reconstruct/real_fields.py::_kde_on_grid` |
| 9.2 | NNLS deconvolution | `reconstruct/deconvolve.py::build_signatures`, `deconvolve_spot` |
| 9.3 | Spearman region assignment | `reconstruct/region_label.py::region_signatures`, `_spearman_matrix`, `label_regions`, `region_field` |
| 10 | conservation, convergence, elasticity | `validation/numerical_rigor.py::_observed_order`, `_elasticity`; `results/numerical_rigor.md` |
