# Landing page: the six beats judged as a sequence, and the overlap fixed at its cause

**Scope:** `landing-3d/**` only. `npx tsc --noEmit` clean, zero console errors in
all four configurations, no horizontal scroll, no viewport overflow.

---

## 1. The mobile overlap was not a mobile overlap

`ship_readiness.md` §3.4 reported a 17 px / 15 px collision at 390×844 between
the copy column and the "Interstitial pressure" readout, called it pre-existing,
and left it for the layout owner. Measuring it across the whole scroll window
rather than at settled positions changed the diagnosis three times.

### 1.1 What the sweep found that a centre screenshot cannot

The probe steps 61 scroll positions per configuration and, at each one, tests
**every visible leaf text node against every other** — not one nominated string
against the page. Visibility is cumulative opacity > 0.12, so elements mid-fade
still count.

| | before | after |
|---|---:|---:|
| desktop 1512×950, distinct colliding pairs | **2** | **0** |
| mobile 390×844, distinct colliding pairs | **6** | **0** |
| horizontal scroll, any position | false | false |
| console errors | none | none |

Three findings the previous pass could not have seen:

1. **The collision was on desktop too**, not only mobile — 162 × 17 px at
   p = 0.35. It was missed because the earlier check ran at each beat's
   *settled peak opacity*, and the overlap only exists while beat 02's callout
   is still lit and beat 03's headline is rising. It is a **transient**.
2. **It hit the `<h2>`, not the `<p>`.** The earlier probe compared the
   nominated body paragraph against other nodes; the element actually being
   crossed was beat 03's headline.
3. **The label's screen position depends on scroll *history*, not just scroll
   position.** The camera and tumour damp per frame, so approaching p = 0.333
   by a jump puts the label at x = [241, 418] and approaching it gradually puts
   it at x = [120, 298]. Any measurement taken after a jump is measuring a
   transient. The final sweep approaches every position gradually.

### 1.2 The cause, and why the previous fix could not have worked

The colliding element is not a floating card. It is a drei `<Html>` label
projected from a **world-space anchor** inside the R3F tree, so its screen
position sweeps as the camera moves — ~130 px on desktop, ~190 px on mobile,
across its visible window. The landed fix (`padding-bottom: 17rem` on the two
beats that carry readouts) reserves the **bottom band**; this label lives
**mid-viewport**. Right defect class, wrong axis.

The real defect was in the timeline, not the layout:

```
fieldLabelOpacity = trapezoid(p, 1/6 + 0.02, 2/6 - 0.05, 2/6 + 0.05, 3/6 - 0.02)
```

That holds beat 02's callout at **full opacity from p = 0.383 to 0.480** — the
middle of beat **03**. Beat 02's own copy is gone by p = 0.333. The label
reached full strength exactly as the text it belonged to disappeared, and was
still lit over the next beat's headline. **The callout outlived its beat by
half a beat**, and that is what put it in the copy column.

### 1.3 What was changed

1. **`lib/timeline.ts` — the callout now dies with its own beat.** The window
   moved inside beat 02 and is fully out by p = 1/3, before beat 03's copy
   begins to rise. This alone removes the desktop collision and the beat-03
   half of the mobile one, and it is the editorially correct timing
   independent of the collision.
2. **`Annotation.tsx` — a `reticleOnly` mode, used at ≤ 820 px.** At 390 px
   the copy column is full-bleed and beat 02 now carries a four-row ledger, so
   there is no band left for a floating callout at any offset. Rather than
   hand-tune a 3D offset against a moving projection, the label and its leader
   line are dropped where there is no room and the **reticle stays** — it is
   the part that still works without space. The words are not lost: the beat's
   in-flow ledger names the quantity *and* its three registered verdicts,
   where the chip carried two lines.
3. **Uniform header clearance at ≤ 820 px.** `sectionVisual` translates a beat
   by `(1 - opacity) × 26 px` while it fades, so beats whose kicker sat at
   y ≈ 64 slid **under** the fixed VIVUM mark mid-transition. Clearance is now
   set to survive the transient in both directions.

### 1.4 The failure mode nobody was checking for: silent clipping

`.scroll-pin` is `overflow: hidden` and `.narrative-section` centres its
content vertically. A block that grows therefore grows toward **both** edges,
and anything past the viewport is **cut, not scrolled**. Adding ledgers pushed
beat 05 to 977 px and beat 02 to 803 px against an 844 px viewport — i.e. the
last registered limits would have been silently truncated.

Both ledger beats are now **top-aligned**, which makes the overflow direction
single and predictable and means copy length can no longer move the top edge at
all — closing the exact regression class recorded in `ship_readiness.md` §3.4,
where a longer string grew a centred block upward into a readout.

Measured at peak opacity, with ≥ 26 px of transition headroom required at both
edges:

| beat | 1512×950 | 390×844 |
|---|---:|---:|
| 02 | 835 / 950 | 705 / 844 |
| 05 | 864 / 950 | **806 / 844** |
| others | all clear | all clear |

---

## 2. The six beats judged as a sequence

Read end to end, the sequence was: **a real thing → a feature list → a real
argument → a real finding → a wall of defensive prose → an ask.** Beats 01, 03
and 04 each do something no competitor page can do. Beats 02 and 05 were where
the sequence lost its nerve, in opposite directions — 02 by saying too little,
05 by saying everything at once.

| beat | verdict | action |
|---|---|---|
| 01 Real tissue | **Strong.** The hero is the thesis: 40,887 measured cells, with a legend turning them from confetti into data. | Untouched. |
| 02 Transport physics | **Weakest.** Headline `"Then run physics on it."` was a connective phrase, not a claim; the body was a four-item feature list with **no number in it**; its most interesting content — that the field is genuinely solved but its magnitude fails its own anchor — was sitting in a **code comment**. | Rebuilt. |
| 03 Delivery or potency | **Strong argument, false sentence.** See §3. | Claim corrected. |
| 04 Spatial attribution | **Strong.** Carries both nulls and the 82 % edge confound in its own prose. | Untouched. |
| 05 Stated limits | **Weak by structure, not by content.** One 1,180-character paragraph holding six distinct registered limits. The fifth and sixth were unreachable. A page whose argument is "we publish what breaks" was burying the breaks. | Restructured. |
| 06 CTA | **Plain.** Ended on the ask, with no connection to the five beats above it. | Sharpened. |

### 2.1 The one risk: beat 05 is now a spec sheet of the product's own limits

The single most opinionated call in this pass. Beat 05 stops being prose that
defends the product and becomes a **ledger of what it cannot do** — six rows,
each with a named quantity, a registered verdict set in mono, and the
measurement that produced it. A landing page whose penultimate beat is a table
of its own failures is genuinely unusual, and it is the honest expression of
this product's whole thesis.

**Nothing was cut to do it.** Every clause of the old paragraph is present as a
row with its number intact: the delivery-law exception, R² = 0.718 blind on 71
compounds, ~24 % of the spread, one cell line, the predecessor's 0.022, the
1.08–1.21× against ~1.7, the 82 % edge confound, spot-platform
non-reporting, absolute kill fraction, patient response. The restructure made
each of them *more* prominent, not less.

Deliberate restraint inside the risk: **the state column carries no colour.**
Red/amber/green would render a registered limit as an exception to an otherwise
working product. Here the limits are the product surface, so they are set in the
same ink as everything else and separated by weight and tracking alone. No new
hue enters the palette.

### 2.2 Beat 02 now states which half of its own physics it can stand behind

New headline: **"The shape is solved. The scale is not."** — a claim, and one
that sets up beat 05.

The same ledger device carries four rows from `results/ifp_validation.md`:
field shape `SOLVED`; absolute scale `FAIL · 41× LOW` against 15 ± 9 mmHg
(PMID 1423283); grid resolution `UNDER-RESOLVED`; glioblastoma `NOT TESTABLE`,
because the published intracranial value is statistically indistinguishable
from normal brain. The grid row is quoted **because** `ifp_validation.md`
registers it specifically so the other two cannot be quoted without it.

Using one device for both beats is the point: 02 and 05 now read as the same
instrument appearing twice, which is the rhythm the two floating panels on 03
and 04 already establish.

**In flow, not floating.** All three floating elements on this page have
collided with narrative text at some scroll position. A disclosure that can be
covered by another element is not a disclosure, so the rows that carry
registered limits sit in normal flow inside the beat that owns them.

### 2.3 Beat 06 closes the loop

It now hands back "a ranking, a per-cell mechanism attribution of what killed
each compound, and the list above attached to both — the same list, on your
tissue, with your numbers in it." It ends on the constraint rather than only on
the ask, which is the whole page's argument in one sentence.

---

## 3. The transport claim: the page said the opposite of the truth

Beat 03 shipped this sentence:

> Weight is the only property reaching transport, so similar-mass analogues
> **cannot be separated**.

True of the model the page shipped when it was written; **false as of
`results/two_route_promotion.md`**. The page's disclosure of its own
falsification had itself gone stale.

It now says structure reaches transport — **and says the narrow version**:

> two routes through the cell layer, paracellular and transcellular, fed by
> lipophilicity and hydrogen-bond donors as well as mass … Structure reaches
> transport now — scored blind against 71 measured compounds — but it recovers
> only about a quarter of their measured spread, so close analogues are
> **under-separated rather than resolved**.

The limit travels in the **same sentence** as the claim, because that is the
clause a reader skips if it is moved to a footnote. Beat 05's ledger carries
the full version: R² = 0.718 blind on 71 measured compounds, ~24 % of the
spread, **from one cell line**, under-separating a congeneric series ~4-fold,
predecessor at 0.022, refutation and replacement both published.

Not written anywhere: "validated", "accurate", "gets transport right", or any
absolute diffusivity figure. Verified mechanically — with comments stripped,
the rendered copy contains none of `digital twin`, `machine learning`, ` AI `,
`personalis/z`, `predict patient`, `validated`, `accurate`, or `monoton`
(monotonicity across quintiles is not cited as corroboration), while all of
`0.718`, `71 measured compounds`, `one cell line`, `0.022`, `1.08`, `1.21`,
`2.0×`, `~1.7`, `82%`, `41×` and `1423283` are present.

---

## 4. Also fixed, also measured

**The reduced-motion fallback is the only experience with flowing text**, so it
is the only one where prose actually travels under the fixed header. Swept at
390×844 it produced 26 distinct overlaps of the VIVUM mark with running text,
left half-legible by a gradient scrim. The fallback now gets a solid ground and
a hairline — a deliberate bar, which is what it has to be once it is covering
prose. All remaining flagged pairs in reduced motion are **VIVUM against
flowing text (0 non-header pairs)**, i.e. correct opaque occlusion, not
collision. The cinematic experience pins its sections, so nothing ever passes
under the mark there and its soft scrim is untouched.

---

## 5. Verified

`npx tsc --noEmit` clean. Zero console errors in all four configurations.
`hScroll = false` at every swept position. No visible element extends past the
viewport at either width.

Screenshots captured at **deviceScaleFactor 2** and read, not just taken:
`1512×950` beats 02/03/05/06, `390×844` beats 02/03/05/06, and the
reduced-motion fallback at both widths on the two ledger beats. Every target is
approached by **gradual scroll**, never a jump, because the damped scene lands
on a different position otherwise (§1.1).

### Not fixed, reported

- `AttributionOverlay` still carries run `r_00067`'s shares under the
  **shipped** delivery law with the corrected law in its footnote. That is
  honest and unchanged, but the open decision it references ("which law ships
  is an open decision") is now one more thing that has moved underneath this
  page than when it was written. Worth re-checking against the engine.
- The Next.js dev-mode indicator appears bottom-left in all dev screenshots. It
  is not a page element and does not ship.
