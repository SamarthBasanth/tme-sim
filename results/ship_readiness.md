# Ship readiness — reading every claim surface together

**Date: 2026-08-13.** A final consistency pass over every user-facing claim
surface before the narrow claim — *relative comparison of compounds on one
tissue section* — goes to a first real research user.

**The question this file answers, up front: is this ready to put in front of a
first user?**

> **YES, with one condition: the fixes recorded in §3 must ship with it.**
>
> Before this pass the answer was **no**, for one reason. **Not a single surface
> a user could reach told them that the transport model underlying their ranking
> is falsified for similar-mass compounds** — except the bring-your-own-compound
> path, which most users will never touch. A researcher ranking a congeneric
> series against the shipped GBM or breast library would have received a
> confident, well-caveated ranking whose delivery column carried no
> compound-specific information at all, and nothing on the README, the landing
> page, the exported report, the Studio run view or the CSV would have said so.
> That is now fixed on every surface this pass owns.
>
> The remaining defects (§4) are wrong *documentation*, not wrong *product*. None
> of them causes a user to over-read a result; the surfaces that carry a number
> to a user are now consistent with what is true.

---

## 1. What was checked, and how

Everything below was verified **at source** — by running the shipped tool or
reading the shipped code — not by trusting a summary. Three of the briefing
notes I was given turned out to be wrong or stale, and are corrected in §2.4.

| Surface | How it was checked |
|---|---|
| `vivum validate --quick` | Run live. `PASS 7 · FAIL 6 · NOT TESTABLE 4 · DEFERRED 8` of 25 rows, `13 of the 13 scored rows never invoke the transport PDE`, **exit 1**. |
| `vivum validate --only breast-kill-3d2d` | Run live. `PASS 3 · FAIL 1 · NOT TESTABLE 6 · DEFERRED 0` of 10 rows, exit 1. |
| `vivum compare-drugs --csv` | Run live at `dx = 30` and `dx = 80`, GBM, temozolomide + cisplatin. CSV + `.provenance.json` sidecar inspected in full. |
| Studio report generator | `build_limits()` and `render_html()` driven directly from a stored report; limits enumerated and the rendered document read as text. |
| Studio API | Live against the running backend: `/api/sections`, `/api/reports`, `/api/reports/{id}/html`. |
| Studio UI caveat text | `provenance.js`, `inspect.js`, `spec.js`, `compound.js`, `projects.js`, `reports.js`, `protocol.js` read for every disclosure string. |
| `README.md`, `docs/`, `app/README.md` | Read in full and diffed against the above. |
| `landing-3d/lib/copy.ts` | Read; edited; `npx tsc --noEmit` clean; screenshotted and read at 1512×950 and 390×844, with and without `prefers-reduced-motion`. |
| Python suite | `pytest -q` → **950 passed**, 8 warnings, 691 s. |

---

## 2. Contradictions found *between* surfaces

Piecemeal correction by many agents produced exactly the failure mode expected:
individual surfaces are careful, and they disagree with each other.

### 2.1 The one that mattered — nothing told a user the transport model is falsified

`results/diffusivity_falsification.md` establishes that the shipped `D_eff` is a
function of **molecular weight and the ECM fraction only** — `logP`, HBD, HBA,
TPSA and formal charge are computed by `descriptors.compute_descriptors` and
discarded before transport. Against 71 compounds (Pruijn 2008, HT29
multicellular layers, one method, one cell line) the measurement spans **29.7×**
where the engine spans **2.4×**, at **R² = 0.028**. At *identical* molecular
weight the engine is flat **by construction — exactly 1.00×** — where the
measurement spans **3.7×** (four compounds at 263 Da) and **3.8×** (four at 305 Da).

Grepping every claim surface for any mention of this returned:

| Surface | Disclosed? |
|---|---|
| `README.md` | **No** — zero mentions |
| `landing-3d/lib/copy.ts` | **No** — and §03 actively implied the opposite ("descriptors to Stokes-Einstein diffusivity") |
| Exported Studio report (HTML) | **No** — zero mentions across the whole rendered document |
| Studio `STANDING_LIMITS` (shown under every run) | **No** |
| `compare-drugs` CSV export | **No** — and it exports a `diffusivity_um2_s` column |
| `docs/` (every file) | **No** — `diffusivity_falsification.md` and `two_route_transport.md` were cited **zero times** anywhere in `docs/` |
| Studio bring-your-own-compound path | **YES** — `app/server/api.py:681-691`, `app/web/js/ui/spec.js:385-390`, `app/web/js/ui/compound.js:134-136` |

So the disclosure existed, and it was gated on `is_user_compound`. A user ranking
the **curated library** drugs — the documented core use case, the thing
`compare-drugs` is for — saw none of it. The briefing's prediction was correct.

**Worse than silence, four surfaces asserted the opposite**, certifying the
falsified chain as sound:

* `docs/molecular_layer.md` graded `effective_diffusivity` **"First-principles
  physics (simulated)"** — the document's highest grade — and summarised that it
  is one of "**two** quantities [that] are real physics". It is the only quantity
  in that table scored against a large external dataset, and it failed.
* `docs/molecular_layer.md` called the parameters "**mechanistically derived**
  from the drug's structure … a genuine upgrade from 'assumed' to 'computed from
  the molecule.'" It is computed from the molecule's *mass*.
* `docs/molecular_layer.md` described validating against measured tissue PK as
  "**the necessary next step**". That step has been taken and the layer lost.
* `docs/QUICKSTART.md` attributed the comparison table's differences to
  "**structure-derived diffusivity**" — directly above a `D_eff` column spanning
  1.5× across three drugs — and called the output a prediction of "**validated
  transport** + potency physics".

`docs/MATHEMATICAL_MODEL.md` §10.5 additionally instructed readers to build error
bars from a list that **excluded §2.1–§2.5**. Followed literally on a similar-mass
series, that produces an interval that cannot contain the truth, because the
transport spread there is zero by construction.

### 2.2 The breast tally — README frozen mid-reconciliation

`results/zero_width_bands.md` landed **2026-08-12** and settled the question. The
README still said the reconciliation was pending:

> *"Reconciliation is in progress and will be recorded in
> `results/zero_width_bands.md`; until it lands, read both sources and quote
> neither count as settled."*

It had landed. Confirmed by running the shipped tool: `PASS 3 · FAIL 1 · NOT
TESTABLE 6`. Two consequences propagated:

* README §"What has failed" said the full tally "is in flux". The corrected
  full-run tally is **PASS 13 · FAIL 18 · NOT TESTABLE 19** of 50 rows.
* README's improvement list still had *"(0) reconcile `vivum validate` with
  `kill_calibration_breast.md`"* as the top priority. Done.

This is the **understatement** failure mode: a settled, defensible result being
withheld because the page had not caught up.

### 2.3 `docs/QUICKSTART.md` disagreed with `README.md` about the same command

Both describe `vivum validate --quick`. README was right; QUICKSTART was wrong on
every figure, and wrong in the **optimistic** direction on FAIL:

| | QUICKSTART (before) | Live tool |
|---|---|---|
| `--quick` | PASS 7 / FAIL 4 / NOT TESTABLE 3 / DEFERRED 6, "13 registered anchors" | PASS 7 / **FAIL 6** / **NT 4** / **DEFERRED 8**, of 25 rows, **16** anchor groups |
| full run | PASS 13 / FAIL 15 / NT 13, 41 rows | PASS 13 / **FAIL 18** / **NT 19**, **50** rows |
| rows that never touch the solver | "11 rows" | **13 of 13** |
| breast kill panel | "three of the breast 3-D:2-D kill ratios" (fail) | **one** (docetaxel) |

QUICKSTART also claims at the top that *"Every command and every code snippet on
this page was actually run against this repository before being written down;
none of it is aspirational."* That claim was false as written, which is worse
than no claim — it discourages the reader from re-running.

### 2.4 Three briefing notes that were wrong at source

Reported per the instruction to verify rather than trust:

1. **"R² = 0.022"** — both numbers are real and they are different subsets.
   `diffusivity_falsification.md` reports **R² = 0.028** for the engine vs
   measurement over all **71** compounds. `two_route_transport.md` §0 reports
   **0.022** for the shipping engine over the **68** compounds carrying a
   logP₇.₄ and therefore scoreable. I used 0.028/71 on every surface, since that
   is the falsification file's own headline and the larger set.
2. **"2.4× modelled"** — `diffusivity_falsification.md` says **2.4×** over 71;
   `two_route_transport.md` says **2.32×** over 68. Same resolution.
3. **"a 9.58× delivery error moved the scored ratio at most 1.45×"** — verified
   correct, but the README's own qualification is load-bearing and was left
   standing: the residual is never zero (13.8 % at best, ~48 % for five
   low-permeability drugs), and the source file's summary is *"'insulated' is not
   'immune'."*

### 2.5 Two false alarms, chased and cleared

Recorded because both looked like ship-blockers:

* **The served report appeared to contain withdrawn claims.** Fetching
  `/api/reports/rep_00006/html` from the running backend returned a document
  asserting *"The oxygen PDE is FALSIFIED as value-add"* (retracted) and
  *"kill_fraction is grid-converged for dx ≤ 20 µm"* (withdrawn, both halves).
  **Cause: the running server process predates the fix.** `Report.full_view()`
  calls `build_limits()` **live** at call time, so the on-disk code is correct;
  recomputing from current source yields the corrected limits. `first_user_readiness.md`'s
  finding is genuinely fixed and has **not** regressed. *Operational note: restart
  the Studio backend before a demo — a long-lived process serves the code it was
  started with.*
* **`bbb_penetrance` appeared to move with grid spacing** (0.2423 vs 0.1254 for
  temozolomide across two CSVs). Re-running `dx = 30` and `dx = 80` on current
  code returns **0.12537158358322562** at both. The 0.2423 came from a stale
  scratchpad artefact written by another agent under different code. No defect.

---

## 3. What was fixed

All on surfaces this pass owns. `pytest -q` → **950 passed**. `npx tsc --noEmit`
clean.

### 3.1 The falsified transport model now travels with every ranking

| Surface | Change |
|---|---|
| `README.md` | New **"Read this before you rank a compound series"** box directly under the one-sentence claim, stating the 29.7× / 2.4× / R² 0.028 / 1.00×-at-matched-MW result and naming similar-mass analogues as *the* case the product gets wrong. New row in **"What has failed"**. Falsification named in **"The honest summary"**. |
| `app/server/report.py` | New limit `diffusivity_mw_only_falsified`, severity **critical**, positioned as limit **2c** — immediately after the delivery-law limit, with a comment explaining why it must sit with the ranking limits and not in the numerics block at the bottom. Verified: renders as one of **12** limits in the block headed *"read these before any number on this page"*. |
| `app/web/js/ui/provenance.js` | New entry in `STANDING_LIMITS`, which renders under **every** run's caveats regardless of configuration. Phrased with **"MUST NOT"** deliberately, so `CV_HARD` in `inspect.js` counts it as a **hard stop** in the caveats tally. |
| `landing-3d/lib/copy.ts` | Guardrail **1b** added to the header block (the best-calibrated text in the repo). §03 corrected: "descriptors to Stokes-Einstein diffusivity" → "molecular **weight** to a Stokes-Einstein radius", plus "weight is the only property reaching transport, so similar-mass analogues cannot be separated". §05 ("Stated limits") carries the full disclosure, since that is its job and it has the room. See §3.4 for a layout regression this caused and how it was resolved. |
| `docs/QUICKSTART.md` | Falsification box placed directly under the `compare-drugs` output table, pointing at the 1.5× `D_eff` spread visible in the transcript above it. |
| `docs/molecular_layer.md` | §3 box; epistemic grade changed **"First-principles physics (simulated)" → "FALSIFIED against measurement"**; "real physics" summary rewritten; "computed from the molecule" corrected; "necessary next step" replaced with what the step actually returned. |
| `docs/MATHEMATICAL_MODEL.md` | §2.5 **Status: FALSIFIED** box; §10.5 error-budget list corrected to lead with it. |
| `docs/architecture.md` | New entry in §9 "Honest limitations". |
| `docs/DATA_ACCESS_PLAN.md` | Added to the **"What was falsified"** row of the outward-facing summary table. |

In every case the two-route fix is named alongside — **R² 0.725, Spearman +0.842,
blind and unfitted, 24 of 31 exact-MW pairs against 0 of 31** — and in every case
it is stated to be **merged, wired into nothing, and default-OFF**, so it cannot
be read as something the product does today.

### 3.2 Stale pessimism corrected (the README-as-apology failure mode)

Every one verified against source before editing:

| Claim | Reality | Fixed in |
|---|---|---|
| *"No active transport or efflux. There is no P-glycoprotein (P-gp) or other efflux/influx transporter model."* | `tmesim/molecular/efflux.py` and `ion_trapping.py` both exist, are scored, and are opt-in default-OFF | `docs/molecular_layer.md`, `docs/architecture.md` (×2) |
| *"there is no efflux or carrier transport in the model"* — offered as the **causal explanation** for `bbb_penetrance` Spearman 0.50 | That explanation is unavailable; the gap is that *this score* does not read the efflux term | `docs/architecture.md` |
| discrete vessel network *"Built and self-tested but **DORMANT** — nothing imports it"*, *"Wiring it in is the highest-value structural fix outstanding"* | `_DEFAULT_VASCULAR_MODE = "discrete"` (`real_fields.py:179`) — it is the **shipped default** | `docs/architecture.md` (×2), `docs/REPRODUCE.md` (×2) |
| `ranking_robustness` / `vessel_network_impact` *"**INCOMPLETE** — no `__main__`"* | Both have `__main__`; results files and figures exist | `docs/architecture.md` (×3), `docs/REPRODUCE.md` (×2) |
| *"**Until it lands**, report cisplatin (0.241) vs doxorubicin (0.231) as a **TIE**"* | It landed and the measurement **contradicted** that prior (P(cisplatin rank 1) = 0.9277). A live instruction to misreport a result. | `docs/REPRODUCE.md` |
| molecular layer *"specified by the contract"*, only 2 of 8 modules implemented | All implemented; the layer is ~22 modules, 13 of them undocumented in §2 | `docs/molecular_layer.md` |
| breast tally *"under reconciliation — do not quote either figure"* | Settled **3 PASS / 1 FAIL / 6 NOT TESTABLE** | `README.md`, `docs/QUICKSTART.md` |

### 3.3 Overstatements and stale numbers corrected

* `docs/QUICKSTART.md`: all four tallies regenerated from live output; *"validated
  transport + potency physics"* removed; *"structure-derived diffusivity"* removed.
* `docs/MATHEMATICAL_MODEL.md` §2.7: published the **retracted arithmetic mean**
  `¼Σdᵢ` as the shipped equation. Corrected to the geometric mean
  `(∏dᵢ)^(1/4)` that `permeability.py:325` actually computes — in the one document
  that promises *"no equation is stated here that the code does not solve."*
* `docs/MATHEMATICAL_MODEL.md` §10.5: asserted `kill_fraction` grid-converged.
  Withdrawn.
* `docs/REPRODUCE.md` "correction 5": restated **two refuted** diagnoses (the
  26–28 µm degeneracy and the ~12 % flat floor) as current. Marked superseded,
  with what survives.
* `docs/REAL_BUILD_PLAN.md`: *"already validated on real penetration data"* removed.
* `docs/real_data_status.md`: whole file is a pre-correction snapshot quoting a
  retracted `bbb=0.558`, a **retired** `penetration_index`, a bare absolute kill
  fraction, and a legacy pipeline. Given a prominent **ARCHIVED RECORD — DO NOT
  QUOTE** header enumerating each defect.

### 3.4 A layout regression I caused on the landing page, and caught

Worth recording because it is invisible to `tsc`, to the console and to a
horizontal-scroll check, and because commit `d43cfa4` ("Landing: stop the readout
panels covering the finding they belong to") fixed this exact class already.

The §03 copy block is **vertically centred**, so making it longer grows it
**upward** as well as downward. My first draft added ~220 characters, which
raised the block's top from y = 450 to y = **392** and drove it into the
"Interstitial pressure · Solved · relative, uncalibrated" readout chip
(y = 397–414) — a **17 px overlap**, visible in the render as the chip's label
being cut mid-word by the body text drawn over it.

Established by measurement, not by eye: injecting the *original* string at
runtime and re-measuring gave a clear 36 px gap, so the collision was mine. The
copy was tightened until the desktop block returned to **exactly** the original
geometry (top 450, bottom 653) with **zero** collisions. The disclosure that was
cut from §03 is carried in full by §05, which has the room.

Final state, measured at peak opacity in all four required configurations
(1512×950 and 390×844, each with and without `prefers-reduced-motion`), checking
the edited block against every visible text node on the page:

| | §03 | §05 |
|---|---|---|
| desktop | 0 collisions | 0 collisions |
| mobile | **2 pre-existing** (see below) | 0 collisions |
| reduced-motion desktop | 0 collisions | 0 collisions |
| reduced-motion mobile | 0 collisions | 0 collisions |

`hScroll = false` everywhere, **zero console errors** everywhere,
`npx tsc --noEmit` clean.

**The two mobile §03 collisions are pre-existing and are NOT mine.** Same runtime
A/B: with the original copy the overlaps are **17 px and 15 px**; with the current
copy they are **17 px and 15 px** — identical. At 390 px wide the readout chip
sits inside the copy column's vertical span regardless of the text. Reported, not
fixed: the layout owner should decide, and it is the same defect class `d43cfa4`
addressed on desktop.

**No existing honest disclosure was weakened.** The delivery-law caveat, the
`CAVEATS` hard stops, the provenance chips, the NOT TESTABLE reasons, the
"insulated is not immune" qualification and the `NOT_FOR_CLINICAL_USE` banners are
all untouched except to add to them.

---

## 4. What remains

### 4.1 On surfaces this pass does not own — reported, not fixed

1. **The `compare-drugs` CSV export does not carry the falsification.**
   `tmesim/cli.py` — off-limits this pass. The CSV is otherwise the best-disclosed
   artefact in the product: a `#` provenance block carrying section, dx, α and its
   single-anchor provenance, the full delivery-law caveat, `NOT FOR CLINICAL USE`,
   and a machine-readable sidecar. But it **exports a `diffusivity_um2_s` column**
   — the falsified quantity — with no note, and its `what_is_defensible` line reads
   *"The ORDER, under the delivery law named above"* without the similar-mass
   exception. **This is the highest-priority remaining fix**: it is the export most
   likely to end up in a slide. Suggested: add a `diffusivity_law` key to
   `_comparison_provenance()` and extend `what_is_defensible`.
2. **`results/breast_kill_revalidation.md:61`** still carries the superseded
   *"5 PASS + docetaxel FAIL, 4 untestable"*. `results/` files other than this one
   were out of scope. Both superseded breast counts now survive only there.
3. **`app/web/js/ui/spec.js:679`** — "Copy JSON to clipboard" copies the RunSpec
   with **no preamble and no caveat**. Every other export path in the product
   attaches provenance. Minor, but it is a gap in an otherwise complete pattern.
4. **Landing page, mobile §03**: pre-existing 17 px / 15 px overlap between the
   copy column and the interstitial-pressure readout chip at 390×844 (§3.4).
   Layout, not copy — measured to be independent of the text.

### 4.2 Lower-priority documentation debt

* `docs/QUICKSTART.md` and `docs/REPRODUCE.md` route five *"what is the current
  state, honestly?"* pointers at `results/SESSION_AUDIT.md`, which has known-stale
  sections. `results/CLAIM_INVENTORY.md` is the better target.
* Library drug count disagrees across three places: 10 (`REPRODUCE.md:147`),
  11 (`QUICKSTART.md:125`), 10 (`REPRODUCE.md:409`).
* `docs/molecular_layer.md` §2 documents 9 of ~22 modules (now flagged as partial
  rather than enumerated).
* `docs/VERIFICATION_METHODOLOGY.md` "four external validations run" predates the
  diffusivity falsification, the two-route work, the breast panel and the
  zero-width-band fix. Its §8 (structural vs tunable failure) is the natural home
  for the MW-only falsification — arguably the cleanest structural negative the
  project has — and omitting it weakens the document.

### 4.3 Product facts a first user will hit, unchanged by this pass

Not defects; the honest state of the instrument.

* `vivum validate` **exits 1**. It is meant to. Anyone wiring it into CI should
  expect a red build, and `--quick` is a *weaker* gate than the full run, not a
  faster one.
* **`--dx 12` must be passed explicitly.** `compare-drugs` still defaults to
  `--dx 30`, and dx = 30 values are superseded — a coarse grid cannot resolve the
  perivascular boundary layer. dx = 12 is not converged either. This is stated in
  the README but it is a trap for a user who does not read it, and the default is
  the wrong-by-the-project's-own-account value.
* The **top-1 pick is not decision-grade** — the reportable output for the top of
  the panel is a **SET, not a winner**.
* **Restart the Studio backend before a demo** (§2.5).

---

## 5. Verdict

**Ship it — with §3 included.**

The narrow claim the project chose to defend is real, and it is now stated the
same way on every surface a user can reach: the README, the landing page, the
Studio's standing limits, the exported report's limits block, and the docs. The
honesty machinery is genuinely the best part of this product — `vivum validate`
fails loudly and exits non-zero, the report refuses to rank in `determinant` mode
and carries its limits as §1 rather than a footnote, the caveats panel is
non-collapsible by design, and the CSV cannot be separated from its provenance.

What was missing was not honesty in general but honesty about **the specific case
a first user is most likely to run**: a series of close analogues. A competent
reviewer would have found that in an afternoon — the falsification is published in
this very repository, and no surface pointed at it. That gap is closed.

The residual risk is concentrated in one place and it is named: **the CSV export
still ships the falsified `diffusivity_um2_s` column with no caveat** (§4.1). If a
single further change is made before a user touches this, it should be that one.
