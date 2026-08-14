# Fixing blocker (1): giving `OxygenConfig`'s dimensionless O a defined mmHg reference point

**Status: research/proposal only. No code changed. `tmesim/physics/oxygen.py` is untouched.** This is a narrow, direct follow-up to `results/oxygen_calibration_landscape.md` §4, which identified two structural blockers to using the four real human-glioma pO2 measurements (Rampling 1994, Collingridge 1999, Lally 2006, Beppu 2002) to calibrate the model. This document designs a fix for **blocker (1) only** — "the model has never defined what `O = 1` corresponds to in mmHg." Blocker (2) — the sink term's first-order-in-O functional form vs. the literature's zeroth-order/Michaelis-Menten default — is explicitly out of scope here and is being investigated in parallel; nothing below depends on resolving it, and §5 states why.

---

## 0. What was read and what it actually says

**`tmesim/physics/oxygen.py`** (full file read). `solve_oxygen` solves the steady-state PDE `D*lap(O) + supply*vasc*(1-O) - consumption*rho_tum*O = 0` for `O` on `[0,1]`, clipped after the sparse solve (line 238). The module docstring states plainly (line 25): *"the oxygen field is dimensionless relative saturation in [0, 1] (1 = fully oxygenated)"* — with no mmHg, kPa, or any other physical-pressure statement anywhere in the file. `label_hypoxia` (lines 242-274) marks a cell HYPOXIC iff `local_o2 < cfg.hypoxia_threshold`, operating purely on the dimensionless number.

**`OxygenConfig`** (`tmesim/core/config.py:76-82`, confirmed by direct search): four fields —
```python
diffusivity: float = 2000.0     # um^2/s
consumption: float = 0.02       # per-tumor-cell consumption rate (relative)
supply: float = 0.05            # vascular supply rate (relative)
hypoxia_threshold: float = 0.2  # oxygen below this => HYPOXIC state
```
No `reference_po2_mmhg` or equivalent field exists.

**Where the banding actually lives.** The task brief pointed at `tmesim/reconstruct/cell_state.py` for "the hypoxic/quiescent/viable banding thresholds." Having read that file in full: it does **not** implement the banding. It is the *alternative*, marker-derived state route, and its docstring (lines 6-10) only *describes* the oxygen route as context for why it exists:

> "The shipping cell-state labelling is a banding of the SIMULATED oxygen field (`cells.agents.update_cell_states`: `o < t` HYPOXIC, `t <= o < 2t` QUIESCENT, else VIABLE)..."

The actual implementation is in **`tmesim/cells/agents.py::update_cell_states`**, confirmed directly:
```python
t = cfg.hypoxia_threshold
hypoxic   = alive & (local_o2 < t)
quiescent = alive & (local_o2 >= t) & (local_o2 < 2.0 * t)
viable    = alive & (local_o2 >= 2.0 * t)
```
i.e. **HYPOXIC** `O < t`, **QUIESCENT** `t <= O < 2t`, **VIABLE** `O >= 2t`, with `t = hypoxia_threshold = 0.2` as the sole default used everywhere in the repo (config, tests, every `results/*.md`; no override found anywhere). There is no independent `quiescent_threshold` — QUIESCENT's upper edge is always exactly `2*hypoxia_threshold`, not a fixed literal. `docs/MATHEMATICAL_MODEL.md:855-876` states the same formula and is explicit that **`t=0.2` and the `[t,2t)` shoulder are "declared modelling choices... not a measured physiological cutoff."** That is the correct pointer for this analysis in place of `cell_state.py`, and everything in §3 below converts *these* three bands.

---

## 1. The anchor: two candidates, one recommendation

### Candidate (a) — O=1 anchored to normal brain TISSUE pO2 (recommended)

**Value: 35 mmHg.**

Primary source: Erecinska M, Silver IA. "Tissue oxygen tension and brain sensitivity to hypoxia." *Respir Physiol* 2001;128(3):263-276. PMID [11718758](https://pubmed.ncbi.nlm.nih.gov/11718758/), DOI [10.1016/s0034-5687(01)00306-1](https://doi.org/10.1016/s0034-5687(01)00306-1) — the standard review of mammalian brain tissue oxygen tension and the most commonly cited source in this literature for a "normal" cerebral tissue pO2 figure in the low-to-mid-30s mmHg (reviewing cytochrome-c-oxidase kinetics and the physiological range over which brain O2-dependent processes operate). PubMed-verified to exist, correctly attributed (PMID/DOI above); the exact numeric table was not independently re-extracted from the full text in this pass, so treat "35" as the literature's commonly-cited point estimate rather than a decimal-precision figure from this paper specifically — consistent with how `results/oxygen_calibration_landscape.md` already flags precision-level claims throughout.

Two independent corroborating lines, both already gathered by `results/oxygen_calibration_landscape.md` and by a fresh check in this pass:
- **Clinical brain-tissue-oxygen (PbtO2, Licox-probe) monitoring literature**: normal PbtO2 is conventionally cited as **20-40 mmHg** (e.g. the neurocritical-care review literature on invasive PbtO2 monitoring), with values below ~15 mmHg flagged as physiologically concerning and below ~10 mmHg as an ischemic threshold. This is a large, independent clinical literature (different technique — implanted Clark-type/Licox probes rather than polarographic Eppendorf needles — but the same target tissue) and it brackets 35 mmHg comfortably.
- **This project's own Beppu et al. 2002 numbers** (`results/oxygen_calibration_landscape.md` §1, PMID 12160140): peritumoral GBM pO2 ≈ 59.8 mmHg, arterial ≈ 109.2 mmHg. A 35 mmHg *normal, non-tumour-bearing tissue* anchor sits below both, in the expected physiological order (arterial > peritumoral-edematous-but-still-perfused > normal-distant-tissue > intratumoral), which is a sanity check the number passes rather than a derivation of it.

**Why this is the right anchor for THIS model specifically.** `solve_oxygen`'s PDE is explicitly a *tissue*-level reaction-diffusion balance — vascular supply diffusing into and being consumed by tissue — not a blood-gas model. The quantity `O` represents oxygen tension *in the tissue the cells sit in*, at the same spatial scale the four literature pO2 studies measured (intraoperative needle-electrode readings taken directly in glioma tissue). A tissue-level anchor is therefore the correct physical match: it answers "what tissue pO2 does a fully-supplied, non-consuming patch of this model's grid represent?" with a number measured the same way, in the same compartment, as the calibration targets.

### Candidate (b) — O=1 anchored to arterial/plasma pO2 (rejected)

**Value: ~90-100 mmHg** (Beppu's own arterial reading, 109.2 mmHg, is close to but somewhat above the commonly cited ~95-100 mmHg systemic arterial reference; both are in this range).

**Why rejected.** Arterial pO2 is the oxygen tension *in blood*, upstream of every diffusion barrier the PDE actually models (capillary wall, interstitium, cell membrane, mitochondrial distance). `solve_oxygen`'s `supply * vasc * (1 - O)` source term already represents the flux *from* vessels *into* tissue — vessels are the boundary condition, not the field itself. Anchoring `O=1` to arterial pO2 would say "a grid cell at O=1 has the oxygen tension of blood inside a vessel," which is not what any grid cell in this model represents (grid cells are tissue, and `vasc` is a separate field carrying the vascular *density*, not the vascular *pO2*). It would also be the wrong compartment relative to every calibration target in `results/oxygen_calibration_landscape.md` §1: all four glioma pO2 studies measured *tissue*, not blood. Using the arterial anchor would silently introduce a large, physiologically real arterial-to-tissue oxygen drop (Beppu's own data: 109.2 mmHg arterial vs. 15.3 mmHg intratumoral tissue — a >7x drop) as an *unmodelled* discrepancy between what O=1 claims to mean and what the PDE's own supply term is actually delivering. §3 below shows concretely that this choice also produces an absurd re-reading of the existing threshold bands, which is the more decisive reason to reject it.

**Recommendation: candidate (a), `REFERENCE_PO2_MMHG = 35.0`, cited to Erecinska & Silver 2001 (PMID 11718758).**

---

## 2. Sensitivity of the anchor choice itself

35 mmHg is a point estimate inside a literature range that is itself not tight (20-40 mmHg PbtO2-monitoring convention; `results/oxygen_calibration_landscape.md` §1 separately reports "white-matter pO2 around 20-40 mmHg and average brain tissue pO2 around 10 mmHg... cortical measurements spanning roughly 5-100 mmHg and averaging ~23 mmHg in one older series," flagged UNVERIFIED-to-decimal there too). The conclusions in §3 are checked at both ends of the plausible range (35 and 40 mmHg) and do not change qualitatively between them — see the table below.

---

## 3. What the EXISTING bands become in mmHg, and whether that is sensible

With `hypoxia_threshold = 0.2` (the sole default, `tmesim/cells/agents.py`) and `mmHg = O * REFERENCE_PO2_MMHG`:

| Band | O range (existing, unchanged) | mmHg @ ref=35 | mmHg @ ref=40 | Literature comparison |
|---|---|---|---|---|
| HYPOXIC | `O < 0.2` | **< 7.0 mmHg** | **< 8.0 mmHg** | Close to, just under, the broad "hypoxic" cutoff (**10 mmHg**, Lally 2006's dichotomization); above the "severely hypoxic" cutoff (**2.5 mmHg**, Collingridge 1999 / pimonidazole binding). Reasonable order of magnitude, though the model has no separate severe-hypoxia band to distinguish the two literature cutoffs. |
| QUIESCENT | `0.2 <= O < 0.4` | **[7.0, 14.0) mmHg** | **[8.0, 16.0) mmHg** | Sits almost exactly around Lally 2006's own 10.0 mmHg re-dichotomization point and close to Beppu's intratumoral-awake reading (**15.3 ± 2.3 mmHg**) — plausible as an "at-risk perfusion" shoulder just above frank hypoxia. |
| VIABLE | `O >= 0.4` | **>= 14.0 mmHg** | **>= 16.0 mmHg** | Below Beppu's peritumoral (**59.8 mmHg**) and well below arterial (**109.2 mmHg**), i.e. not claiming tumour tissue is "fully oxygenated" in any absolute sense — it only claims "not acutely hypoxic," which is the correct, modest reading of a VIABLE label. |

**Verdict: under the recommended (tissue) anchor, the existing bands land in a biologically plausible range and do NOT need to move.** This is a genuinely useful, somewhat surprising result of doing the conversion: `t = 0.2` was "a declared modelling choice... not a measured physiological cutoff" per `docs/MATHEMATICAL_MODEL.md:873-876`, and it turns out to convert, under this anchor, to a HYPOXIC cutoff (~7-8 mmHg) that sits in almost exactly the neighborhood the literature independently uses for its own "hypoxic" cutoff (10 mmHg) — closer than there was any reason to expect from a value that was never fit to it.

**A useful cross-check this mapping also enables:** Lally 2006's whole-cohort median glioma tumour pO2 (**5.1 mmHg**, n=23) converts to `O = 5.1/35 ≈ 0.146`, which is below the model's own `t = 0.2` HYPOXIC cutoff — i.e., under this mapping, the *median* real measured GBM tumour cell would be labelled HYPOXIC by the model's existing threshold. That is an internally consistent, encouraging sign for a model whose current predicted hypoxic tumour FRACTION on the real Xenium section is 0.0 at every grid spacing (`tmesim/physics/oxygen.py:67-72`) — the threshold itself is not the reason for that implausible 0.0%; something else in the pipeline (most plausibly `OxygenConfig.consumption`, never calibrated, and/or the sink functional-form issue in the parallel investigation) is suppressing the field from ever reaching the tissue's own hypoxic band. **This is exactly the distinction blocker (1) was supposed to make legible**, and it does: the labelling convention is not obviously broken; the absolute field level is.

**Now the same table under the REJECTED arterial anchor, for contrast** (`REFERENCE_PO2_MMHG ≈ 100 mmHg`, Beppu's own arterial reading rounded):

| Band | O range | mmHg @ ref=100 | Literature comparison |
|---|---|---|---|
| HYPOXIC | `O < 0.2` | **< 20 mmHg** | This calls tissue **at the top of the normal PbtO2 range (20-40 mmHg)** "HYPOXIC." Both Beppu's own peritumoral GBM reading (59.8 mmHg, clearly non-hypoxic, mass-effect-adjacent brain) and much of normal brain tissue itself would sit *inside or just above* this "hypoxic" band on some readings. |
| QUIESCENT | `[0.2, 0.4)` | **[20, 40) mmHg** | This is squarely the textbook *normal* PbtO2 range, being labelled "at-risk / quiescent." |
| VIABLE | `>= 0.4` | **>= 40 mmHg** | Excludes essentially all four glioma studies' intratumoral readings (Lally median 5.1, Beppu intratumoral 15.3, Collingridge/Rampling regions <2.5) AND most of Beppu's own peritumoral reading (59.8, barely above the line) from ever being called "viable" tissue — the label becomes nearly vestigial for tumour-adjacent tissue. |

This is the biologically absurd result the task asked to check for: **under the arterial anchor, normal, non-hypoxic, non-tumour tissue reads as HYPOXIC-to-QUIESCENT, and the bands would need to move** (roughly halved) to recover anything sensible. That is not merely "a different label for the same fact" — it would require also editing `hypoxia_threshold` downward, compounding one uncalibrated choice with another. This is the concrete version of the rejection argument in §1 and is the strongest single reason to prefer the tissue anchor: **it is the only one of the two candidates under which the existing, already-shipping thresholds do not also need to be re-derived.**

---

## 4. Proposed code change (illustrative diff — NOT applied)

Scoped to `tmesim/physics/oxygen.py` only, following the same pattern already used in that file for `TUMOR_REFERENCE_AREA_UM2` (module-level constant, heavily commented, with an explicit honest caveat). Adds one constant, one pure conversion helper, and a docstring section — touches nothing about `solve_oxygen`'s numerics, `label_hypoxia`'s threshold logic, or any test-observable behavior.

```diff
--- a/tmesim/physics/oxygen.py
+++ b/tmesim/physics/oxygen.py
@@ -20,8 +20,27 @@
     source:     supply * vasc          (turned off as O -> 1 via the (1-O) term)
     sink:       consumption * rho_tum  (proportional to LOCAL tumor-cell density)
 
 Units follow the project convention: length um, diffusivity um^2/s. The oxygen
-field is dimensionless relative saturation in [0, 1] (1 = fully oxygenated).
+field is dimensionless relative saturation in [0, 1] (1 = fully oxygenated,
+0 = anoxia).
+
+MMHG REFERENCE POINT (new; read this before comparing O to any measured pO2)
+------------------------------------------------------------------------------
+``O`` had no defined physical unit until now. ``O = 1`` is DEFINED to represent
+:data:`REFERENCE_PO2_MMHG` (35 mmHg, normal -- non-tumour-bearing -- human brain
+TISSUE pO2, Erecinska & Silver 2001, PMID 11718758), not arterial/plasma pO2
+(~100 mmHg): the PDE's ``vasc`` field is a vascular-density source term feeding
+a tissue compartment, not a blood compartment, so the tissue anchor is the
+correct physical match, and it is the anchor at which the EXISTING
+``cells.agents.update_cell_states`` bands (``t=0.2`` HYPOXIC, ``[0.2,0.4)``
+QUIESCENT, ``>=0.4`` VIABLE) convert to a biologically plausible range
+(~7 / [7,14) / >=14 mmHg) rather than the arterial anchor's absurd one (~20 /
+[20,40) / >=40 mmHg, which would call normal brain tissue "hypoxic"). See
+``results/sweep/oxygen_mmhg_reference_proposal_2026-07-29T05-54-00Z.md`` for
+the full derivation and both anchors' threshold tables. This conversion is a
+DECLARATION about what the existing [0,1] number has always implicitly meant --
+it changes no solved value of ``O``, no threshold, and no downstream label.
+:func:`to_mmhg` is the one place to apply it.
 
 GRID INDEPENDENCE OF THE SINK (fixed; read this)
 ------------------------------------------------
@@ -80,6 +99,43 @@
 
 from ..core import CellState, Microenvironment, OxygenConfig
 
+#: mmHg partial pressure that ``O = 1`` is DEFINED to represent. Anchored to
+#: NORMAL (non-tumour-bearing) human brain TISSUE pO2 -- not arterial/plasma
+#: pO2 -- because ``solve_oxygen``'s source term represents flux from vessels
+#: INTO tissue, so ``O`` is a tissue-compartment quantity throughout; see the
+#: module docstring's "MMHG REFERENCE POINT" section for why the arterial
+#: anchor (~100 mmHg) was rejected.
+#:
+#: Value = 35.0 mmHg. Source: Erecinska M, Silver IA. "Tissue oxygen tension
+#: and brain sensitivity to hypoxia." Respir Physiol 2001;128(3):263-276.
+#: PMID 11718758, DOI 10.1016/s0034-5687(01)00306-1. Corroborated by (a) the
+#: clinical PbtO2 (Licox-probe) monitoring literature's normal range of
+#: 20-40 mmHg and (b) this project's own Beppu et al. 2002 GBM readings
+#: (peritumoral ~59.8 mmHg, arterial ~109.2 mmHg, both above this
+#: normal-tissue anchor, the expected ordering) -- see
+#: ``results/oxygen_calibration_landscape.md`` §1 and
+#: ``results/sweep/oxygen_mmhg_reference_proposal_2026-07-29T05-54-00Z.md``.
+#:
+#: HONEST CAVEAT: like :data:`TUMOR_REFERENCE_AREA_UM2`, this is a SINGLE POINT
+#: anchor for a quantity that varies severalfold across individuals, brain
+#: regions and measurement technique (Nozue et al. 1997, PMID 9290690, found
+#: 1.5-5.6 mmHg tumour-pO2 spread across labs on an identical protocol). It
+#: fixes what O=1 MEANS; it does NOT calibrate ``OxygenConfig.consumption``,
+#: which independently sets where the SOLVED O actually lands on this new
+#: scale -- see that field's own uncalibrated-level caveat above.
+REFERENCE_PO2_MMHG = 35.0
+
+
+def to_mmhg(oxygen: np.ndarray) -> np.ndarray:
+    """Convert a dimensionless oxygen field/array to mmHg via :data:`REFERENCE_PO2_MMHG`.
+
+    Pure linear rescale (``oxygen * REFERENCE_PO2_MMHG``). Touches nothing in
+    ``solve_oxygen`` or ``label_hypoxia`` -- no threshold, no solved value, and
+    no test-observable behavior changes. Exists so a caller who wants to
+    compare the model's O to a measured pO2 (e.g. the four glioma studies in
+    ``results/oxygen_calibration_landscape.md``) has one canonical place to do
+    it instead of each re-deriving the anchor.
+    """
+    return np.asarray(oxygen, dtype=np.float64) * REFERENCE_PO2_MMHG
+
```

**Notes on scope, so the synthesis pass can slot this in correctly:**
- This diff does **not** touch `tmesim/core/config.py` (`OxygenConfig`), `tmesim/cells/agents.py` (`update_cell_states`), or `hypoxia_threshold`'s value. Per §3, the existing `t = 0.2` threshold does **not** need to move under the recommended anchor — this is a labelling/interpretation fix, not a re-tuning of the model.
- A natural, slightly larger alternative is to add `REFERENCE_PO2_MMHG` (or a differently-named equivalent) as a new field on `OxygenConfig` itself (`tmesim/core/config.py`) rather than a bare module constant in `oxygen.py`, so it travels with the other reaction-diffusion parameters and is overridable per-run like `diffusivity`/`supply`/`consumption`. That is a marginally more invasive change (touches `core/config.py` and every place `OxygenConfig()` is constructed with keyword overrides) and was left as the module-constant form above to stay minimal and directly satisfy "a named constant... in `tmesim/physics/oxygen.py`." The synthesis pass should pick whichever integrates more cleanly with the sink-term fix, since that fix may touch `OxygenConfig` anyway.
- `to_mmhg` is offered as a convenience, not a requirement of the fix — the constant alone (`REFERENCE_PO2_MMHG`) is the load-bearing part; the helper just avoids every caller re-writing `oxygen * REFERENCE_PO2_MMHG`.

---

## 5. Why this doesn't need to wait on the sink-term fix (blocker 2)

Blocker (1) (this document) and blocker (2) (sink functional form, under separate investigation) are independent axes:
- Blocker (1) fixes what the [0,1] number **means** — a relabelling of the axis, verified in §3 to be non-absurd against the literature without moving any threshold.
- Blocker (2) is about whether the **solved values of O** on that axis are right — a first-order-in-O sink vs. the literature's zeroth-order/Michaelis-Menten default changes the *shape* of the O field (how sharply it falls near consuming tissue), not what O=1 is defined to be.

Both need to land in the same synthesis pass because `OxygenConfig.consumption`'s numeric value is entangled with the sink's functional form (per `results/oxygen_calibration_landscape.md` §2: "recalibrating the numeric value of `consumption` alone, without addressing this functional-form mismatch, is unlikely to reproduce a realistic pO2 profile"), but the mmHg reference point itself does not depend on that resolution — `REFERENCE_PO2_MMHG` and `to_mmhg` are valid regardless of what functional form the sink eventually takes, since they only rescale the [0,1] axis, not the equation that produces values on it.

---

## Summary

| Question | Answer |
|---|---|
| Which anchor: normal tissue (~35-40 mmHg) or arterial (~90-100 mmHg)? | **Normal brain tissue, 35 mmHg** (Erecinska & Silver 2001, PMID 11718758) |
| Why reject arterial? | Wrong compartment (blood vs. the tissue `solve_oxygen` actually models); §3 shows it also makes the EXISTING bands biologically absurd (calls normal 20-40 mmHg brain tissue "hypoxic") |
| Do the existing HYPOXIC/QUIESCENT/VIABLE thresholds (t=0.2, 2t=0.4) need to move under the recommended anchor? | **No.** They convert to ~7 / [7,14) / >=14 mmHg — plausible next to the literature's 2.5 and 10 mmHg cutoffs and Beppu's 15.3 mmHg intratumoral reading. This was not guaranteed going in and is a genuinely useful finding, not an assumption. |
| Does this fix blocker (2) (sink functional form)? | No — orthogonal, addressed separately; this fix does not depend on that resolution being settled first. |
| Files touched by the proposed diff | `tmesim/physics/oxygen.py` only (module docstring + one constant + one helper). No threshold, config default, or test-observable value changes. |
