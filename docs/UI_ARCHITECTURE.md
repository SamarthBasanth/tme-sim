# Vivum Studio — information architecture

Status: normative for `app/web/`. Written before the consolidation pass, and
kept as the rule the front end is measured against afterwards.

This document exists because the front end was built by several independent
agents over one long session. Each of them shipped something that works. What
none of them could do is agree with the others about where a thing goes, how
tall a header is, what a button is called, or which surface a secondary task
opens on. The result reads as six products sharing a window.

Nothing here is about removing features. Every panel, control and shortcut that
existed before this document stays reachable. What changes is that each one
acquires a **predictable place**, and that place is derived from one sentence.

---

## 1. The sentence

> **You reconstruct a specimen, you run a therapy on it, and you find out why
> the therapy did what it did.**

- The **primary object** is a *reconstructed specimen*: a tissue section plus
  the RunSpec that turned it into a simulable grid. There is exactly one of
  these in view at a time.
- The **primary verb** is *run a therapy and read the mechanism*. Its output is
  a `RunResult`: metrics, per-cell mechanism attribution, provenance, caveats.

"Run a therapy" is a simulation, and the sentence is about layout, not about
what the product is licensed to do. The specimen is a de-identified tissue
section; the run is a model prediction that has been scored against no patient
outcome. Nothing this UI displays is for clinical use — that constraint is
carried in the product by `provenance.js`'s standing limits and the reports'
NOT FOR CLINICAL USE banner, and it is repeated here so a reader who arrives at
the sentence first does not read it as a clinical workflow. See
[`DEPLOYMENT.md` §0](DEPLOYMENT.md) for what a field tester may load.

Everything else in the application — run history, projects, the protocol
editor, therapy reports, the mechanism legend, the render diagnostics — is
**secondary**. Secondary does not mean unimportant. It means it is defined by
its relationship to the sentence, and it is placed accordingly.

Two classifying questions decide where anything lives:

1. **Does it change what will be simulated?** → it belongs to SETUP (left).
2. **Does it explain or measure what was simulated?** → it belongs to READ
   (right).

Anything that answers neither is either context (top bar), instrument health
(status bar), an archive of other runs (dock), or a full-attention task
(workspace overlay).

---

## 2. The zones

Seven zones. Fixed. Nothing moves between runs, between themes, or between
backend states. A panel that has no data shows an empty state **in its own
place**; it never collapses and lets its neighbours slide.

```
┌──────────────────────────────────────────────────────────────────────────┐
│ CONTEXT  identity · section · run · cells ·············· backend · [WS]   │  --h-bar
├────────────┬──────────────────────────────────────────┬──────────────────┤
│            │ SPECIMEN toolbar                         │                  │  --h-head
│   SETUP    ├──────────────────────────────────────────┤      READ        │
│   (left    │                                          │   (right rail)   │
│    rail)   │            the specimen                  │                  │
│            │            + stage overlays              │                  │
│            ├──────────────────────────────────────────┤                  │
│            │ RUNS  (dock)                             │                  │
├────────────┴──────────────────────────────────────────┴──────────────────┤
│ INSTRUMENT  mode · fps · frame · drawn · lod · vbo · index · api · msg    │  --h-status
└──────────────────────────────────────────────────────────────────────────┘
      WORKSPACE overlays (protocol, projects, report viewer, keys) cover
      the region between CONTEXT and INSTRUMENT, and only that region.
```

### 2.1 CONTEXT — the top bar

Answers *what am I looking at, and is the instrument connected*. Read-only
except for the workspace launchers.

| Group | Content |
|---|---|
| identity | `Vivum Studio`, version |
| specimen | active section id + its one-line description |
| run | active run id + status chip |
| scale | cell count |
| — spacer — | |
| link | backend LED + state text |
| launchers | `projects` · `protocol` · `keys` · `theme` |

The launcher group is the **only** place a workspace can be opened from, and it
is owned by one module. Before this pass, three different modules injected
buttons into it independently and one of them did so from a dead code path.

### 2.2 SETUP — the left rail

Answers *what will I run*. Everything that mutates the RunSpec, and nothing
that does not.

Two groups, with the same sticky headers and navigator as the READ rail:

| Group | Sections | Question |
|---|---|---|
| **DEFINE** | Specimen, therapy, model · Source | what am I about to run |
| **COMMIT** | Active job | what is it doing |

Within DEFINE the form keeps its own internal rules — `Pharmacokinetics` and
`Integration` are `hr-lbl` dividers inside one dense grid rather than separate
collapsible sections, because splitting a twelve-field form into four
`<details>` costs four header strips and buys nothing: the whole form is
visible at once at the default rail width.

Pinned to the bottom of the rail, **outside the scroll**: the commit bar —
`Submit run`, `Validate`, and the validation badge. It is pinned because it is
the one control whose reachability must not depend on scroll position.

`Active job` — progress, stage, elapsed, log tail — is the consequence of that
bar, not a fifth thing to configure, which is why it is a group of its own
rather than a fourth item under DEFINE.

### 2.3 SPECIMEN — the centre stage

Answers *what does it look like*. The toolbar carries only things that act on
the **viewport itself**: pointer tool, zoom/fit, live cursor read-out, colour
source, unpin.

Stage overlays occupy four fixed corners plus one band, and are never allowed
to overlap each other:

| Slot | Overlay | Shown when |
|---|---|---|
| top-left | hover / pinned cell read-out | always |
| top-right | *(free)* | — |
| bottom-left | scale bar | geometry loaded |
| bottom-right | field legend | a field overlay is selected |
| bottom band, above the corners | time-lapse scrubber | a time-lapse is loaded |
| full stage, under everything | annotation layer (profile / region / distance) | always |
| full stage, over everything | empty state | no geometry |

The empty state is the only overlay allowed to cover the others, because when
it is showing there is nothing underneath it to cover.

The scrubber's clearance above the two bottom corners is **measured**, not a
constant: `ui/overlays.js` publishes `--ov-bottom` from the taller of the scale
bar and the field legend. The constant it replaced went stale the moment the
legend grew a line.

**The time-lapse does not open on its first frame.** Frame 1 is `t = 0`, which
is empty by construction — no drug has entered the tissue yet — so a stage sold
as *see where the drug goes* opened on a picture with no drug in it. It now
opens on the frame carrying the most drug (argmax of the mean over the sampled
grid; see `openingFrameIndex` in `js/ui/mechanism.js`), which is the last frame
for a static plasma source and the peak-exposure frame for a clearing one. The
rule has one obligation attached, and it is not negotiable: **the picture must
never be mistaken for `t = 0`.** Three independent disclosures carry the clock —
the scrubber read-out promotes `t =` to accent at `--fs-md`; the scrubber's
metadata line states the default in words and names the control that returns to
`t = 0`; and the field legend, which is the key to those pixels, prints
`t = … · frame n/N` under its title. The status bar prints it a fourth time.

**The mechanism legend used to hold the top-right slot, and no longer does.**
It was the largest overlay at 300 px wide and ~380 px tall, and this document
booked that as a permanent cost against the specimen, offering "a collapse
control on the card, not a mode coupling" as the escape hatch if it ever became
a complaint. It became one, and both of those options were answers to the wrong
question: the problem was never the card's width, it was that a card sat **on
top of** the product's primary output — occluding roughly the top-right quarter
of the tissue, in every colour mode, including the two it is not the key for.

The rejection of a mode coupling still stands ("tying it to `mech` mode would
hide a result behind a display setting"), which is exactly why the fix is not
to hide the panel but to **move** it: it is now a section of the READ rail,
under RESULT, immediately below `Metrics` (§2.4). Nothing about it is
conditional and nothing about it is collapsed by default. See `js/ui/mechanism.js`
for the full argument; the short version is that the breakdown is read once per
run and then used episodically, that once-per-run reference numbers are what
the READ rail is *for*, and that the breakdown belongs beside the `Kill fraction`
metric it decomposes and reconciles with bit-exactly.

The interactive half of the panel — rows that toggle a category's visibility in
the viewer — moves with it, and has precedent: `type` colouring has no stage key
either, and its per-channel toggles live in the rail under VIEW > Channels.
Selecting the `mech` colour source scrolls the rail to the section, because in
that mode the breakdown is also the key to what is on screen.

### 2.4 READ — the right rail

Answers *what does the run say, and how am I looking at it*. This rail carried
ten flat, unlabelled sections and 2350 px of content in a 1000 px column. It
now carries the same ten sections, plus one that was missing and one that came
in off the stage (`Mechanism`, §2.3), organised into
**three named groups** with a **sticky group navigator** at the top of the rail.

| Group | Sections | Question |
|---|---|---|
| **RESULT** | Metrics · Mechanism · Caveats · Provenance · Reconstruction report | what did the run find, and how much of it is measured |
| **VIEW** | Channels · Field overlay · Window / level · Render | how is the picture being drawn |
| **MEASURE** | Line profile · Regions · Distance | what did I measure on it |

RESULT is first because it is the payoff of the primary verb. VIEW is second
because it is used continuously. MEASURE is last because it is episodic. The
navigator makes the order irrelevant for reaching a section, so the order is
free to express priority.

Rules for this rail:

- Every section is a `<details>` and every one is **open by default**, with a
  single exception: `Render`, which holds point-size, alpha, LOD and the
  benchmark. Those are properties of this machine's renderer, not of any run or
  any specimen, so they are the one thing in the rail that is not about the
  primary object. Collapse is otherwise a user action, never a default, and
  never used to hide primary data.
- `Caveats` is never truncated and never scrolled inside itself.
- Every tool with a stage gesture has a section here. Before this pass,
  `distance` had a toolbar button and a keyboard shortcut but no panel, no
  read-out, and no way to clear it — the only one of the four tools without a
  home.

### 2.5 RUNS — the dock

Answers *what else have I run*. Three panels, three columns, one topic: **runs
that are not the active run.**

| Panel | Content |
|---|---|
| Run history | every run, with A/B comparison selectors |
| Comparison | the A-vs-B spec diff, metric deltas, and caveat differences |
| Therapy reports | multi-arm report jobs and the generated documents |

These three are complements, not alternatives — the comparison panel is driven
by the radio columns in the history table beside it — so they stay side by side
rather than becoming tabs. The dock is resizable; the columns are not
reflowed by run state.

### 2.6 INSTRUMENT — the status bar

Answers *is the renderer healthy*. Render statistics and API latency only. It
is the single owner of that information: a floating "render statistics" card
that duplicated it exactly was removed in an earlier pass and must not return.

### 2.7 WORKSPACE — the overlay class

A workspace is a task that deserves the whole window: writing a protocol,
organising projects, reading a generated report, reading the key map. Before
this pass there were **four different implementations** of this idea — a
`position: fixed` full-bleed div, two `<dialog showModal>` elements with three
different backdrop colours and three different header treatments, and a third
dialog with yet another.

There is now **one** surface contract, `.sheet`:

- covers the region between CONTEXT and INSTRUMENT, and only that region, so
  the operator never loses the specimen identity or the backend state;
- one header height, one border, one backdrop, one elevation, all from tokens;
- `Esc` closes it, and the launcher that opened it shows `aria-pressed`;
- it never becomes the only way to reach a feature.

Members: protocol workspace, projects, report viewer, keyboard map.

Only the first two have a top-bar launcher. The report viewer opens by
selecting a report in the dock and the keyboard map by `?` or the `keys`
button — each is reached from the thing it is about, which is the rule.

---

## 3. The visual system

`css/tokens.css` is the **single source**. Every other stylesheet consumes
tokens and declares none. Before this pass, `protocol.css` and `mechanism.css`
each opened their own `:root` block and shipped their own light-mode variants
with a *different selector form* from the one in `tokens.css`.

### 3.1 Type

Six sizes, no others. Swiss neo-grotesque for prose, monospace with
`tabular-nums` for every numeral, identifier, unit and code path.

`--fs-micro` 9.5 · `--fs-xs` 10.5 · `--fs-sm` 11.5 · `--fs-md` 12.5 ·
`--fs-lg` 15 · `--fs-xl` 19

### 3.2 Space

One scale. `--sp-1` 2 · `--sp-2` 3 · `--sp-3` 4 · `--sp-4` 6 · `--sp-5` 8 ·
`--sp-6` 10 · `--sp-7` 13 · `--sp-8` 18.

The audit found gaps and paddings of 3, 4, 5, 6, 7, 8, 9, 10 and 13 px used
interchangeably across the eight stylesheets. Values off the scale are
normalised to the nearest step.

### 3.3 Height

The audit found **eight distinct header-strip heights** — 19, 20, 20, 21, 23,
23, 26, 26 px — across five stylesheets, all meaning "the title bar of a
thing". There is now one:

| Token | Value | Used by |
|---|---|---|
| `--h-bar` | 30 | the top bar, and only the top bar |
| `--h-head` | 24 | **every** secondary header: panel head, section summary, stage toolbar, sheet head, editor head, tool bar |
| `--h-status` | 21 | the status strip |
| `--h-row` | 22 | a dense table row |
| `--h-th` | 20 | a table header cell |
| `--h-ctl` | 22 | input, select, default button |
| `--h-ctl-sm` | 18 | small button |
| `--h-seg` | 20 | a segmented-control button |
| `--h-chip` | 13 | a chip |

### 3.4 Colour roles

Four families, and they do not mix:

- **neutral** — `--bg-void|0..4`, `--line-0..2`, `--fg-0..3`. Cool-biased
  (hue 210–214, low chroma) so it recedes behind warm tissue imagery.
- **accent** — one desaturated instrument cyan. Interactive affordance,
  selection, focus. Nothing else. Never decoration.
- **semantic** — `--ok` / `--warn` / `--crit` / `--info`, each with a `-bg`.
  State only. Deliberately far in hue from the accent so "selected" can never
  be read as "healthy".
- **categorical** — `--ch-*` (cell type) and `--mech-*` (mechanism
  attribution). Identity, not state. Validated for CVD separation. These now
  live in `tokens.css` with the rest, not in the stylesheet of whichever module
  happened to need them first.

Plus `--syn-*` for the TPL syntax palette — one extra hue (indigo) for numbers
and units, chosen far from both the accent and all three semantic states.

Light mode is **not** an inversion. Panels go to paper white but the viewport
surround stays a deliberate mid-grey, because a white surround around rendered
tissue destroys apparent contrast.

### 3.5 Elevation

`--elev-1` for a panel that raises off the shell. `--elev-2` for something that
genuinely floats: a sheet, a toast. Nothing else casts a shadow. `--backdrop`
is one value, used by every modal surface — the audit found three different
hand-rolled backdrop colours in three stylesheets.

### 3.6 Radius

`--radius: 2px`. That is the ceiling, not a suggestion.

---

## 4. Accessibility and robustness

These are part of "clean", not a separate checklist.

1. **Focus is visible on everything focusable.** `:focus-visible` draws a 2 px
   `--focus` outline. No rule may hide content in order to draw it — the
   canvas carries an explicit exemption from the global promotion rule for
   exactly this reason (promoting `#gl` to `z-index: 5` on focus once made
   every stage overlay disappear behind it when the operator clicked the
   viewport).
2. **`prefers-reduced-motion` must not strobe.** Blanket-compressing every
   animation to 1 ms turns a 1.1 s pulse into a strobe, which is worse than the
   motion it suppresses. Anything that loops opacity or position is stopped
   outright and keeps a static form that still reads: the LED stays lit, the
   indeterminate bar stays a stub.
3. **The body never scrolls horizontally.** Wide content — tables, sweep grids,
   the toolbar — scrolls inside its own container.
4. **Every pane owns its overflow.** `html, body { overflow: hidden }` is the
   contract; a pane that forgets `min-height: 0` pushes the shell instead.
5. **Colour is never the only channel.** Every mechanism legend row ships a
   name and a count; every provenance chip ships its word.
6. **No inline `style=` in markup.** A strict CSP with `style-src 'self'`
   blocks inline style attributes; utility classes exist so `index.html` needs
   none. Styles set from JS via `element.style` are not governed by CSP and are
   used freely.

---

## 5. Naming

- A control is named for **what it does**, not for the module that built it or
  the endpoint it calls. Endpoint names belong in `title`, not in the label.
- Two controls that do the same thing get the **same** label. Two controls that
  do different things get **different** labels. The audit found two buttons
  both labelled `refresh` sitting 600 px apart in the same dock, hitting
  different endpoints, and one action (`download` / `download html`) with two
  labels.
- A toggle states the state it produces, not an ambiguous noun.
- A label must not describe behaviour the code does not have.

---

## 6. What this document does not license

- Removing a feature to make the picture simpler.
- Hiding a control in a menu that used to be visible.
- Adding a dependency, a build step, or a CDN reference.
- Editing `tmesim/` or `app/server/` to make the front end easier. If the
  client and the contract disagree, the disagreement is reported, not
  silently resolved.
