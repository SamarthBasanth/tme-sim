# Landing page: Vivum -> Codon rename, typeface replacement, closing CTA

Scope: `landing-3d/` only. Backend/CLI rename is a separate agent's work and was
not touched here.

Date: 2026-08-14

---

## 1. Rename

Every occurrence of `Vivum` / `vivum` inside `landing-3d/` is now `Codon` /
`CODON`, matching the capitalisation convention at each site.

| File | Site |
| --- | --- |
| `app/page.tsx` | header wordmark `VIVUM` -> `CODON` |
| `app/layout.tsx` | `<title>` and meta description |
| `lib/copy.ts` | guardrail docstring, beat 01 body, beat 06 body |
| `app/globals.css` | file header comment, mobile-padding comment |
| `components/canvas/tumorGeometry.ts` | docstring |
| `components/canvas/shaders/viridis.ts` | docstring |
| `components/scroll/ScrollExperience.tsx` | mail subject line |
| `components/scroll/ReducedMotionExperience.tsx` | mail subject line |

Verified: `grep -rniI "vivum" landing-3d/` (excluding `node_modules`, `.next`,
`package-lock.json`, `tsconfig.tsbuildinfo`) returns **nothing**.

### The CLI was deliberately not named

The brief flagged that some copy invokes the CLI as `vivum validate`. Checked:
**no command string appears anywhere in `landing-3d/`** — not in copy, not in
the components, not in the metadata. So there was nothing to rename and nothing
to invent.

**Corrected 2026-08-14 (later):** the rename has since landed.
`pyproject.toml` `[project.scripts]` now registers **both** `codon` and `vivum`
against `tmesim.cli:main`, with `--help` printing `codon`; the `vivum` alias is
kept so the dated `results/*.md` commands still reproduce. So `codon` is now
safe to assert.

The page still names no command, now by choice rather than by caution. Beat 06
is about handing over the Studio and bringing a section; a CLI invocation would
change its register and lengthen a beat that item 1 below already lengthened. If
a future line would be strengthened by naming it, `codon` is the name.

### Not renamed (third-party strings)

Left alone on purpose: `10x Genomics`, `Xenium`, `GDSC2`, `viridis` /
`matplotlib`, `CC BY 4.0`, PMIDs, `Visium`. None contain the product name; all
refer to third-party things.

---

## 2. Typeface

### What was replaced and why

IBM Plex Sans / Sans Condensed / Mono. The objection was not that Plex is a bad
face for this subject — it was commissioned for technical and scientific
documentation, which is exactly the register — but that it is the family
Anthropic's Claude interface is set in, so the page read as a skin of a chat
product rather than as an instrument.

Excluded by the brief and honoured: Inter, Space Grotesk, Geist. Those are the
current house style of generated design, and swapping into one would trade one
generic look for another.

### What replaced it

| Role | Face | Loaded as |
| --- | --- | --- |
| Display (voice) | **Archivo**, variable, `wdth` axis loaded | `weight: "variable"`, `axes: ["wdth"]` |
| Prose | **Atkinson Hyperlegible Next**, variable | `weight: "variable"`, normal + italic |
| Mono (all data) | **Atkinson Hyperlegible Mono**, variable | `weight: "variable"` |

All three self-hosted via `next/font/google`, which fetches at build time and
serves the woff2 from this origin.

### Why these, in terms of the subject

The page's thesis is not "the simulator is good". It is that the simulator's
**limits are published, legible and checkable** — beat 05 is a spec sheet of the
product's own failures. So the typographic problem is not "look scientific". It
is: *a reader must be able to take a number away without misreading it, and must
actually finish a caveat.* `lib/copy.ts` states this outright — "a limit a reader
does not finish is not a disclosure". That converts legibility from taste into a
product requirement.

**Atkinson Hyperlegible** was commissioned by the Braille Institute and drawn for
low-vision readers, and its explicit design goal is *character disambiguation*:
`0`/`O`, `1`/`l`/`I`, `6`/`8`, `5`/`S` are separated by structural difference — a
slashed zero, a tailed `l` — rather than by a hairline that dies at 10px on a
dark ground. Every quantity on this page is meant to be transcribed and checked
against a `results/` file. A face engineered so a digit cannot be confused for a
letter is an instrument-grade choice, not a stylistic one. This is the rare case
where the type specimen and the product argument are the same sentence.

Taking **both** text roles from one 2025-matched superfamily is load-bearing, not
tidiness. `.ledger__row` sets a sans label and a mono verdict on one shared
baseline, and both instrument panels pair a sans `<dt>` against a mono `<dd>` the
same way. A mono drawn against its own sans keeps x-height and baseline in
agreement, so those pairs sit level without per-element nudging.

**Archivo** (Omnibus-Type) is a revival of late-19th-century American grotesques
drawn for *highlighting* text in newspapers and magazines — headline lettering
whose job is to report a finding, not to sell one. That is the register these
headlines are in: "The shape is solved. The scale is not." is a result, not a
promise.

The `wdth` axis is requested for a structural reason. Hierarchy on this page is
carried by width contrast against the body face, not by size alone, and
`.scroll-pin` clips at `overflow: hidden` — so a display face that sets one line
wider than its predecessor does not degrade, it silently truncates a disclosure.
A continuous width axis lets the display be tuned to the column (headings at
`font-stretch: 80%`, ledger labels at 92%, wordmark at 88%) instead of hoping a
fixed condensed cut fits.

### Tabular numerals

Satisfied structurally: Atkinson Hyperlegible Mono is monospaced, so every glyph
carries the same advance width and columns of figures align by construction. The
existing `font-variant-numeric: tabular-nums` declarations were kept as explicit
statements of intent rather than removed as redundant.

### The three roles were re-tuned, not re-pointed

Every metric was re-set for the new faces rather than inherited:

* Headings: `font-stretch: 80%`, weight 650 (was condensed 600 — at 80% the
  stems thin and 600 read grey), sizes down ~7% because Archivo's x-height is
  larger for the same em, leading opened `1.02 -> 1.06` / `1.06 -> 1.10`,
  tracking `-0.005em -> -0.012em`.
* Prose: `clamp(0.98,1.4vw,1.12) -> clamp(0.92,1.28vw,1.04)`, leading
  `1.62 -> 1.66`. Atkinson's very large x-height set optically larger and
  measurably wider at the old size, adding lines to the two beats with no room
  for them.
* Mono: every uppercase label lost roughly a third of its tracking
  (`0.2em -> 0.13em` kickers, `0.1em -> 0.06em` panel headers,
  `0.12em -> 0.07em` ledger caption, `0.11em -> 0.055em` verdicts,
  `0.34em -> 0.26em` scroll cue). The new mono is wide by design — open
  apertures are how it keeps `0` apart from `O` — so it arrives pre-tracked. Left
  at the old values the panel headers overran a 19rem card.
* Reduced-motion headings now use the display face. They never set a
  `font-family`, so they had been silently inheriting the body sans and the
  fallback was shipping a two-role page while the cinematic one shipped three.

### The mono/prose rule survives — audited, not assumed

`scratchpad/monorule.mjs` walks every visible leaf text node across the whole
scroll and classifies it. Result over 90 distinct strings: **48 mono, 25 sans,
17 display, 0 violations.**

* Every measured quantity is mono: `40,887`, `8.09 x 5.07 mm`, `58.2%`,
  `FAIL · 41x LOW`, `1.08-1.21x, NOT 2.0x`, `0.172`, `0.825`, `Run r_00067`.
* Every prose sentence is sans.
* The display face carries **zero digits** — voice only, never data.

(The audit script flags two h2 headlines as "prose not in sans"; those are the
display face doing its job and are the heuristic's false positives, confirmed by
reading them.)

### Proof the fonts actually load

This matters because the *original* stack failed exactly here: a webfont was
declared and never loaded, so every heading silently resolved to Helvetica Neue
for the life of the page. A font that fails to load does not error — it renders.

`scratchpad/fontcheck.mjs`:

* All three families report `status: "loaded"` in `document.fonts`
  (Archivo `100 900 / 62% 125%`, Atkinson Next normal + italic, Atkinson Mono).
* **Zero requests to `fonts.googleapis.com` or `fonts.gstatic.com`** — the
  `@font-face` `src` values are same-origin `../media/*.woff2`. Nothing for the
  strict CSP to block.
* Width axis is live: the same string at `font-stretch: 80%` measures 968.9px
  against 1168.2px at 100%.
* Computed styles on real elements resolve to the intended family in all three
  roles.

---

## 3. Closing CTA now points at the software

Beat 06's primary action was `mailto:samarthbasanth@gmail.com`. It is now
`STUDIO_URL` from `lib/env.ts` — the same target the hero's "Open the Studio"
button uses, so the first and last frames offer the same door.

A page whose whole argument is "run it and check the numbers yourself" was
terminating in a contact form, which puts a human gate in front of the one thing
a visitor can verify unaided.

Email survives as a deliberately quieter secondary route (`or email us about
your tissue`), set as a text link rather than a second pill so it cannot compete
with the primary action. Changed in **both** experiences — `ScrollExperience`
and `ReducedMotionExperience` — so the two do not drift.

The pair is stacked rather than set side by side: measured on the rendered page,
a row centres the pill and link as a *group*, putting the pill's own centre
~120px left of the axis the headline and paragraph are centred on.

---

## 4. Verification

### Collisions — gradual approach, both widths

Jump-to-position measurement lies here: the camera damps per frame, so arriving
at a scroll position by jumping puts a projected label up to ~120px from where
arriving gradually puts it. `scratchpad/gsweep.mjs` walks the page in 8
`requestAnimationFrame` micro-steps between stops and settles 260ms before
probing, at 81 scroll positions.

| Config | Before | After |
| --- | --- | --- |
| 1512x950 | 0 collisions | **0 collisions** |
| 390x844 | 0 collisions | **0 collisions** |

No horizontal scroll at either width (`scrollWidth == clientWidth`), no viewport
overflow, no console errors.

The sweep was run three times — after the typeface swap, after the closing CTA
was restacked, and after the §5b copy corrections lengthened beat 06 — because
each of those changed layout. All three runs returned **zero at both widths**,
with no hScroll and no console errors. The mono/prose audit was also re-run
after the copy change: still 0 violations.

### Reduced motion, both widths

`scratchpad/rmsweep.mjs`, 61 settled positions each.

| Config | Flowing-content collisions | Header ground |
| --- | --- | --- |
| 1512x950 reduce | **0** | `rgb(7,8,10)`, alpha **1.0** |
| 390x844 reduce | **0** | `rgb(7,8,10)`, alpha **1.0** |

A naive geometric probe reports 14 (desktop) / 35 (mobile) overlaps in the
fallback, all of them the fixed `CODON` wordmark against prose scrolling
beneath it. That is the one experience where content is in normal document flow,
and `globals.css` already handles it on purpose: under
`prefers-reduced-motion: reduce` the header gets a **solid** ground plus a
hairline instead of a gradient scrim, precisely so text passing under it is
occluded rather than left half-legible. The CSS records the same measurement
("the mark overlapped running text at 34 distinct positions" at 390x844), so
this is pre-existing and by design, not introduced by the typeface change. The
sweep now asserts the opacity rather than assuming it — if that ground ever goes
translucent the run fails.

### Clipping

`.scroll-pin` clips silently, so every beat's body box was measured against the
viewport at its own peak opacity, arriving gradually:

| Beat | 1512x950 | 390x844 |
| --- | --- | --- |
| 01 | 170-796 | 210-723 |
| 02 | 120-783 | 92-662 |
| 03 | 293-673 | 143-518 |
| 04 | 307-660 | 143-518 |
| 05 | 120-825 | 92-766 |
| 06 | 308-658 | 289-644 |

Beat 06's figures are post-copy-correction: §5b(a) lengthened that body by a
sentence, growing the box ~27px on desktop and ~13px on mobile. Both still fit
with >180px of headroom. Measured twice by two independent scripts (the full
peak-finding pass and a targeted mobile probe) which agreed exactly on
`[289,644]`.

All fit. Beat 05 at 390px lands at 766px against the 844px viewport — inside the
~805px budget `globals.css` records as the limit for surviving its own fade
transition.

### Typecheck

`npx tsc --noEmit` -> exit 0, clean.

---

## 5. No disclosure was weakened

Nothing in `lib/copy.ts` changed except three occurrences of the product name.
Specifically unchanged:

* All four rows of beat 02's IFP ledger, including `FAIL · 41x LOW` and the
  `UNDER-RESOLVED` row that exists so the other two are never quoted alone, and
  the `NOT TESTABLE` GBM row.
* All six rows of beat 05, including `ABSOLUTE KILL FRACTION / NOT DEFENSIBLE`
  and `PATIENT RESPONSE / OUT OF SCOPE`.
* Beat 03's "recovers only about a quarter of their measured spread" clause, in
  the same sentence as the claim it limits.
* Beat 04's `1.08-1.21x` honest null and the 82% tissue-edge confound.
* Both instrument panels' nulls (`1.0x` well-mixed, `~1.7x` clustering-matched)
  and footnotes.

**Beat 05's state column still carries no colour.** It is set in the same
`--ink` as everything else, separated from the label by mono, uppercase,
tracking and a hairline rule alone. No hue was added to the palette. Verified by
reading the rendered screenshots at both widths, not by reading the CSS: a
red/amber/green vocabulary would render a registered limit as an exception to an
otherwise working product, and on this page the limits are the product surface.

No claim of clinical response prediction, diagnostic use, AI/ML, "digital twin",
or personalised medicine was introduced. No absolute kill fraction appears.

## 5b. Governance copy corrections (2026-08-14, second pass)

Three corrections from a governance audit. Applied; one wording was changed from
what was proposed, for cause.

### (a) The data-locality promise was false — beat 06

Was: *"Codon runs inside your own infrastructure; proprietary tissue and
compound data never cross your firewall."*

Two defects. It described the **customer's perimeter**, which is not ours to
promise, and it was concretely false while the API shipped wildcard CORS with no
default auth — any page the operator browsed could read the compound library, so
the data did have a path out.

Now: *"Codon runs where your tissue already is: the solve is a single process on
your own machine, and it sends no tissue or compound data anywhere.
De-identified tissue only."*

**I did not use the proposed "makes no outbound network requests", and this is
the one place I pushed back.** That half is true of `app/server/` but FALSE of
the package as a whole, and falsifiable by grep in about thirty seconds:

* `tmesim/data/acquire.py` — `import requests`, hits
  `https://api.gdc.cancer.gov/{cases,files,data}`, the Allen Institute API and
  the Sanger GDSC release bucket.
* `tmesim/data/spatial_ingest.py` — `requests.get(url, stream=True)` against
  `cf.10xgenomics.com`.
* `tmesim/data/real/mgmt_spatial.py` — `urllib.request.urlopen` against NCBI GEO
  and its FTP mirror.

Those are deliberate, operator-invoked downloads of **public reference
datasets** that carry nothing of the customer's. But an HTTP GET is an outbound
request, so the sentence would have been literally untrue — and shipping a
second false absolute to replace the first one would have repeated the exact
defect the audit was raised to fix. On a product whose entire pitch is "check our
numbers", a marketing claim that fails a grep is the worst available failure
mode.

So the claim is about what the software **sends**, not about how many sockets it
opens. That is the property that actually matters to a customer, it is narrower,
and it survives the grep. Verified:

* `app/server/` imports no HTTP client, no socket library, and no analytics or
  crash-reporting SDK (`sentry`/`posthog`/`mixpanel`/`segment`/`gtag`: zero
  hits across `app/` and `tmesim/`).
* The Studio front end's only external URLs are an SVG namespace and an OFL
  licence URL — both inert strings, neither fetched.
* This landing page self-hosts all three fonts and issues **zero** third-party
  requests at runtime (measured in §2). The typeface decision is load-bearing
  for the privacy claim, not just for the look.

Also honoured: **"upload" appears nowhere.** There is no tissue upload —
sections are registered by directory path, so the exposed surface is the path,
id and label, and the label becomes the `<h1>` of every exported report
(`docs/DEPLOYMENT.md` §0).

### (b) The guardrail block had no data-provenance clause

That is the hole (a) went through: the block forbade clinical, diagnostic,
treatment-planning, AI and digital-twin claims, and said nothing about data. A
`DATA PROVENANCE` section was added at the head of `lib/copy.ts` stating the
de-identified-tissue scope, banning the word "upload", and carrying the rule
that **any data-locality or privacy claim must be traceable to a verified
property of the software, never to the customer's network**.

It records the failed claim verbatim so it cannot drift back, notes that the
failure was caught by auditing the software rather than by reading the sentence,
and — because a guardrail that only forbids invites a worse replacement — states
what *is* sayable and explicitly blocks the "no outbound network requests"
upgrade, naming the three modules that falsify it.

### (c) "validated" in the meta description

`app/layout.tsx` claimed the simulation was *"validated against real measured
tissue benchmarks"* — the last surviving unqualified "validated" on any
marketing surface, on the highest-visibility one there is, violating guardrail
1b (`Do NOT write "validated"`) in the same repository that defines it, while
the README records that `codon validate` reports more failures than passes.

Now: *"scored against real measured tissue benchmarks with the failures
published alongside the passes."* Not a softened boast — publishing the ones you
lose is the part competitors do not copy.

### Re-verified after the copy edits

Item (a) made beat 06's body materially longer, so the checks were re-run rather
than assumed. `npx tsc --noEmit` exit 0. Results in §4.

## 6. Sibling work preserved

The visual-craft pass already on disk was left intact and re-verified working:
aspect framing, beat-02 retiming, beat-06 scrim, the `:focus-visible` ring, the
rail active tick, and the scroll cue (`components/ui/ScrollCue.tsx` plus its
`globals.css` block). The scroll cue's mono label was re-tracked for the new
mono like every other mono label; nothing was reverted.

## 7. Screenshots read

`scratchpad/codon_{desk,mob}_beat0{1..6}.png` and
`scratchpad/codon2_{desk,mob}_beat06.png`. Read, not merely captured.

Everything left unstaged.
