# The first user journey — 2026-09-20, 03:41 → 04:20 CDT

**Who wrote this.** An agent with no prior knowledge of the workspace build,
driving a real Chromium (148.0.7778.96, Playwright) against a real
`app.server` process, as a newcomer. Nothing here is read off the code alone:
every claim below is something that happened on screen, with the step list that
produced it.

**Servers used.** Four processes, all `./.venv/bin/python -m app.server`:

| Port | Started | Purpose |
|---|---|---|
| 8771 | 03:41:15 | the journey server — steps 1–7 |
| 8772 | 03:52 | `TMESIM_MAX_GRID_CELLS=100000` — the small-host refusal |
| 8775 | 04:03 | `TMESIM_STATE_DIR=<tmp>` — the empty first-install workspace |
| 8776 | 04:14:50 | current code, for the re-verification pass |

**THE CODE MOVED UNDER ME.** Roughly ten agents were writing these files while
I read them. `app/web/js/home/projectGrid.js`, `app/web/js/util/dom.js`,
`app/web/js/ui/reach.js` and `app/server/workspace.py` all changed mtime during
the run. Every finding below therefore carries the time I observed it and
whether it still reproduces on the 04:14:50 process. Two findings were fixed by
other agents while I was writing; they are kept, marked FIXED, because the gap
that let them ship is not.

**One deliberate deviation.** Defect 1 below made the home page render nothing
at all. To see steps 2–7 I served `projectGrid.js` with that one selector
repaired **in the browser only**, via a Playwright route interceptor. No file in
the repo was modified by me. Everything after step 1 was walked on that
override until 04:14:50, when the real fix landed and I re-walked the entry
path unpatched.

---

## Part 1 — the walkthrough

### 1. Landing on `/`

**What actually happened first: nothing worked.**

`http://127.0.0.1:8771/` at 03:44. The page painted a card reading
**"Reaching the server · 0.0s elapsed"** with a live indeterminate progress bar
and the sentence *"Asking this server for its projects and its most recent runs
— one call to GET /api/workspace."*

I waited 30 seconds. Then 30 more. It never changed. `#home-body` had zero
children; `aria-busy` had been set back to `false`; the bar kept animating.
`/api/workspace` had in fact answered in under 50 ms — I confirmed it with
`curl` in another window, 200 OK, one project and nine runs.

This is the exact thing contract §6 forbids: *"A refusal says what to do next.
Never a bare dash, never a spinner that is really an error."* The front door of
the product was a spinner that was really an error, and it would have stayed
that way forever.

Cause and fix are under **Defect 1**.

**Once the page renders (from 04:14:50 onwards, unpatched, on current code):**
time from navigation to a populated project grid was **0.23 s**.

What I understood in the first ten seconds, honestly:

* Within ~2 s: this is a list of projects on a server, and there is a Studio
  somewhere else. The H1 "Workspace" plus one subtitle sentence do that work
  well.
* Within ~6 s: I knew what a project *is* — because the ordering-basis
  paragraph under the controls explains `last opened` versus `last modified` in
  plain words, and the card spells out "2 specimens / 9 runs attached / not
  opened yet".
* Within ~10 s: I knew what to press. But only just. **The one thing a new
  arrival should press — `New project` — is 81 × 22 px of 11.5 px type, in the
  far top-right of the controls bar.** Body copy on this page is 11.5 px
  against a 16 px root; the H1 is 19 px. That is the Studio's instrument
  density applied to a page that is a document, not an instrument. At 200 %
  zoom it reads beautifully (see step 10). At 100 % on a 1440 px display it is
  small.

The one thing I did **not** understand in ten seconds: what the product
actually does. "The simulator itself is the Studio" tells me where, not what.
On an empty server that gap is filled superbly (step 11) — but the moment one
project exists, that orientation disappears entirely and a first-time visitor
to a *colleague's* server gets none of it.

### 2. Creating a project

**2 clicks minimum. 0 questions I could not answer. 0 mandatory fields.**

Press `New project` (click 1) → dialog opens with focus already in the name
field, pre-filled `Study 2026-09-20`, labelled *"Suggested — type over it.
Names do not have to be unique; the id is the identity."* Notes are explicitly
optional and explicitly not read by the model. "First specimen" is pre-checked
with a section already selected. Press `Create project` (click 2).

The footer sentence before I committed read: *"Creates a project holding 1
specimen — xenium_319. No run is started and nothing is simulated: the server
assigns an id (p_…) and you choose the therapy in the Studio."* That is the
best single sentence on the whole surface. It told me exactly what pressing the
button would and would not cost.

Two things were wrong here:

* The section picker offers all 38 sections as `id — label` and **drops every
  qualifier the Studio shows on the same list**. In the Studio the same entries
  read `GBM-Space AT10-BRA-5-FO-1_1 [ingested — not validated] [4 refused]`.
  On the home page they read `gbmspace_at10_bra_5_fo_1_1 — GBM-Space
  AT10-BRA-5-FO-1_1`. A newcomer picking from this list has no way to know that
  34 of the 38 are not validated and refuse most compounds. (**Defect 8**)
* After the dialog closes, focus is dumped on `<body>`. (**Defect 7**)

### 3. Opening it in the Studio

**1 click.** The card's `open` is `<a href="/studio?project=p_00003">`.

The project context follows, and it is unambiguous. The Studio's top bar reads
`← workspace | PROJECT p_00003 First-user walkthrough | SECTION xenium_319 | …`
and `POST /api/projects/p_00003/open` fired — I confirmed `last_opened_at` went
from `null` to `1789893884.96` in the workspace payload.

Is it obvious which project a run will land in? **Yes, and the honest answer is
"none of them until you say so."** The first-run tour card states it in as many
words: *"You came in from a project, and the top bar names it. A run is not
filed under it automatically: when one finishes you are offered a specimen to
attach it to, and nothing is attached until you pick one."* That matches
behaviour exactly (step 5).

The tour opens automatically on the first Studio visit and covers the whole
screen. I checked: it fires **once** — visit 2 in the same browser context had
no tour. So a first-timer pays one extra click (`Skip tour`).

**Cold landing to a running simulation: 5 clicks.**
`New project` → `Create project` → `open` → `Skip tour` → `Submit run`.
A returning user: 4.

### 4. Setting up and submitting a run

**Controls I had to touch to submit: zero.** The shipped RunSpec is complete
and valid — xenium_319, temozolomide, plasma peak 50 uM, 21,600 s, dx = 12 um,
9 further fields all at their defaults and labelled `9 FIELDS, ALL DEFAULT`.

The left rail carries **79 visible controls**, but they are not 79 decisions.
The structure is three numbered questions (`1 · WHICH SPECIMEN?`,
`2 · WHICH THERAPY?`, `3 · AT WHAT RESOLUTION?`) exposing **5 primary fields**:
section, drug, plasma peak, total time, grid dx. The remaining ~74 are the
bring-your-own-tissue panel, the bring-your-own-compound panel, the JSON
editor, and one `Remove ingested section` button per ingested section (31 of
them) — all inside collapsed disclosure panels. That is a reasonable structure;
the count is only alarming if you count it the way a machine does.

I touched exactly one field (`GRID DX` 12 → 30) and pressed `Submit run`.

#### `dx = 12` on a small host — is the refusal legible BEFORE submitting?

**Yes. Emphatically yes. This is the best thing in the product.**

On port 8772 (`TMESIM_MAX_GRID_CELLS=100000`), with the shipped default
`dx = 12` still in the box and before I touched anything:

```
GRID BUDGET                                   WILL BE REFUSED
(8,088 / 12) x (5,065 / 12) um
284,498 grid cells
Ceiling 100,000 grid cells (configured on this server through
TMESIM_MAX_GRID_CELLS), reached on this section at dx = 20.24 um.
This spec is at 284.5 % of it.

THE SERVER WILL REFUSE THIS SPEC
dx_um = 12 um on this section gives a grid of ~284,498 voxels, above the
100,000 ceiling this server is configured with (TMESIM_MAX_GRID_CELLS) —
one float64 field alone would be 2.3 MB, and a solve allocates several plus
the cell population. Use a coarser dx_um.

Reproduced from the refusal in app/server/api.py so that what you read
before submitting is what you would read after. The server writes the real
one; if the two ever differ, the server's is the true one.
```

The footer read `1 ERROR`, the `1 error from the validator — press to go to the
first one` button was live, and `#btn-submit-run` was `disabled` — Playwright
sat on it for 30 s and never got a click through. The refusal names the number,
the ceiling, where the ceiling came from, the break-even dx (20.24 um), and
what to do. Nothing is left to be discovered by submitting.

Two small things: the shipped default (`dx = 12`) is invalid out of the box on
such a host, and nothing offers "use the coarsest dx this server allows" in one
press. And `#btn-submit-run` uses `disabled` with no `aria-disabled`, so it
leaves the tab order entirely — a keyboard user tabs past the primary action
and never learns it exists. The adjacent `1 ERROR` button is focusable and does
explain, so this is a wrinkle, not a hole.

### 5. Watching it finish, and reading the result

Run `r_00100` (temozolomide, xenium_319, dx = 30) finished. The READ rail's
first heading is `THE FINDING — what this run says, in cells`:

> Of **23,812** tumour cells measured in this section, **21,142** survived this
> run. Every cell below is one of those survivors, labelled by what the model
> says stood in the way.
>
> **21,111 · 99.9 %** had drug arrive, and it stayed below the half-effect
> concentration
> **31 · 0.1 %** never had drug arrive in a pharmacologically meaningful amount

**Can I tell what the answer is?** Yes. It is stated in counts with the
denominator beside them, and the headline is survival, not kill. The next block
is labelled `This is a POTENCY result, not a delivery result`, and the one
after it is `PART OF THIS VERDICT IS AN INPUT: plasma_peak / IC50 = 50 uM /
606.0 uM = 0.083`.

**Can I tell what it does not claim?** Yes, and this is where the product is at
its most unusual. Without scrolling past the first screen of the rail I was
told, unprompted:

* the ratio the verdict turns on was fixed before any tissue was read;
* the shipped delivery law `S = PS*Cp` has no return path and is **unbounded
  above**, and this run's median tumour cell sits at 3.70× the plasma peak —
  *"That is the bound broken, on screen, by this run"*;
* `THIS VERDICT REVERSES FOR TEMOZOLOMIDE — THE COMPOUND IN THIS RUN`, with the
  census reference and the four compounds on this section whose verdicts hold;
* `WHAT WOULD REFUTE THIS`, naming imaging mass spectrometry or microdialysis
  on a matched block and the concentration that would falsify each half.

I went looking for the prohibited framings on this surface and found none. No
cross-compound ordering; no absolute kill fraction as a headline (the left
rail's first warning is *"ABSOLUTE KILL SHOULD NOT BE QUOTED FROM SECTION
'xenium_319' UNDER ANY state_mode"*); no clinical or treatment-planning
framing; no "digital twin"; no AI/ML claim; no loose "validated"; and every
caveat I measured was set **larger** than the number it qualifies (11.5 px prose
against 10.5 px values on the project card — measured with `getComputedStyle`
on every leaf node in `#projects`).

One thing did go wrong at this moment — twice, `renderLadder is not defined` —
see **Defect 10**.

### 6. Back to `/` — is my work there?

**No. And the page promises it will be.**

I came back the way a user would, by the top bar's `← workspace`. The card for
`First-user walkthrough` read:

```
NO COMPOSITION YET
no run is attached to any of this project's 1 specimen(s), so nothing has
been reconstructed yet.
RUNS            0 runs attached
```

Which is **true** — I had not accepted the attach offer. But four lines down,
the `RECENT RUNS · across every project, newest first` strip listed
`r_00008 … r_00016`, all *"finished 1 month ago"*, and did not contain
`r_00100`, the newest run on this server by 54 days, finished four minutes
earlier from inside that very project.

The strip shows only runs already attached to a project. The server had **107
runs in `app/runs/`**; the workspace header said **"9 runs"**. Both numbers are
defensible in isolation. Together, on one screen, they are two different truths
about the same machine, and the one the user just made is in neither.

And `recentRuns.js`'s own empty-state copy makes the promise explicitly:
*"Open a project and submit a run in the Studio; it will appear here with the
project and specimen it belongs to."* I did exactly that. It did not appear.
(**Defect 3**)

Once I accepted the attach, everything was correct and immediate: the card
gained `MEASURED COMPOSITION xenium_319 · 40,887 cells` with the stacked bar,
`2 runs attached`, `cisplatin  doxorubicin` under *"2 compounds run, listed
alphabetically. This is a set, not an order: nothing here ranks them."* The
composition caveat — *"This is a chart of that measurement. It is not a picture
of the tissue and nothing about its shape corresponds to where cells are"* —
is exactly contract §4 and is set larger than the percentages it qualifies.

One inaccuracy in the other direction: `LAST OPENED not opened yet` renders
correctly for a project that predates the field, and `document updated 1 month
ago (2026-07-28 10:56:48)` always carries the absolute date beside the relative
one. Contract §3 is honoured.

### 7. Deleting something — does it tell the truth about what it destroys?

**Yes, completely, and I verified the claim on disk.**

I duplicated `Two-patient GBM therapy panel` first, so that I would be deleting
something with real runs under it. The duplicate dialog said:

> **RUNS ARE NEITHER RE-EXECUTED NOR DUPLICATED.** A run lives in `app/runs/`
> and a project holds only its id, so after this both projects point at the
> SAME runs. … Two projects agreeing about a run is one run counted twice, not
> two runs that agree.

I checked: `p_00004`'s specimens carried `['r_00008' … 'r_00016']` — the same
ids, not copies. Correct.

Then the delete dialog for `p_00004`:

> **DELETED** · The project document — its name, its notes, and its id p_00004,
> which is not reused. · Its 2 specimens — the bindings between this project and
> those sections. · The 9 run references those specimens hold — the
> associations, not the runs.
>
> **RUNS ARE NOT DELETED.** The 9 runs this project references stay in
> `app/runs/`, which is the single store for run results. A project holds run
> IDS, never results, so deleting one detaches its runs and never deletes them.
> They stay in the Studio's run history, and in any other project that
> references them.

Measured before and after the delete: **107 run directories in `app/runs/`
before, 107 after.** `p_00002`'s specimens still carried all nine ids.
`GET /api/projects/p_00004` → 404. Contract §2 is kept and the UI says so at
the point of deletion, exactly as §2 requires.

Two notes. The delete needs **two presses of the same button inside the
modal**: `Delete project` relabels to `delete p_00004?` and only the second
press issues the `DELETE`. I recorded "delete is broken" before I worked that
out. The arming is sensible; it is announced by nothing — there is no
`aria-live` region in the dialog, so a screen-reader user gets silence.
(**Defect 9**)

And with the duplicate present, the workspace header read **"3 projects on this
server, holding 5 specimens and 18 runs"** — 18 for 9 distinct runs. The
duplicate dialog names that exact failure mode one screen earlier.
(**Defect 5**)

---

## Part 2 — breaking it deliberately

### 8. Reload mid-run — **passes**

Submitted `r_00102` (doxorubicin, dx = 12). At t+7.6 s the server reported
`status: running, stage: compile_drug, progress: 0.5`. I reloaded at t+13.1 s,
with the server reporting `stage: simulate, progress: 0.55`.

After the reload, with no errors:

```
ACTIVE JOB /status poll | r_00102 | RUNNING | Cancel
STAGE simulate 55 % · 12s
04:10:20 queued: {"drug": "doxorubicin", "dx_um": 12, … }
04:10:20 running
```

and the top bar read `RUN r_00102 RUNNING 12s`. The panel's own promise —
*"The run executes in the Vivum server process, not in this tab. Closing the
tab or reloading does not stop it, and this panel picks up a run that is still
in flight when the page loads"* — is exact. When it finished, the attach offer
appeared normally.

### 9. A project id that does not exist — **passes**

`/studio?project=p_99999`. Top bar: `PROJECT p_99999 — unavailable`. A toast:

> **Project p_99999 could not be opened.** The Studio has started WITHOUT a
> project context: the server has no project with that id. It may have been
> deleted, or the link may be mistyped. Everything else works exactly as it does
> with no project at all — runs still execute and are still stored in
> `app/runs/` — but a finished run will not be offered for filing, because
> there is no project to file it under. The address bar still carries
> `?project=p_99999`, so a reload retries it; the "← workspace" link in the top
> bar leaves it behind.

No page errors. `GET /api/projects/p_99999` and `POST …/open` both return a
clean `404 {"detail":{"message":"no such project 'p_99999'"}}`. This is the
model for how every refusal in the product should read.

### 10. `localStorage` unavailable — **passes**

Both `window.localStorage` and `window.sessionStorage` redefined to throw
`SecurityError` before any script ran. `/` rendered the full project grid; 0
page errors. `/studio` rendered both rails, 72,155 characters of body text, 0
page errors. Every access is already inside `try { … } catch { /* private mode
*/ }`.

### 11. Keyboard only — **passes on `/`, one real hole in the Studio**

Tabbing from the top of `/` with nothing else touched: **30 stops, then back to
the document.** Every stop was inside the viewport, every stop had a visible
`outline: solid 2px`, and the order is the reading order — skip link, nav,
search, sort, `New project`, then per card `name → open → rename → duplicate →
delete`, then the nine recent-run links, then the three footer links. No trap,
no ghost stop, nothing reachable only by mouse.

Dialogs: the create dialog opens with focus in the name field, cycles its 6
controls, and `Escape` closes it **and returns focus to `#pgrid-create`**. The
delete dialog opens with focus on `Cancel` — the safe one — and cycles 3.

The holes: focus after a successful create (**Defect 7**), no announcement of
the delete arming (**Defect 9**), and the disabled submit button leaving the tab
order (step 4).

### 12. 200 % zoom — **passes, and is arguably better than 100 %**

Emulated as a 720 × 450 CSS viewport at DPR 2, which is what a 1440 × 900
window at 200 % browser zoom gives.

`document.scrollWidth == clientWidth == 720`. **No horizontal page scroll.**
The only element past the right edge is the `.skip` link at x = −9999, which is
intentional. The grid reflowed from 3 columns to 2, the shell scrolls in y
(1352 px of content in a 450 px viewport, `overflow-y: auto` on `#home-shell`),
and every caveat paragraph became comfortable to read.

### 13. 360 × 740 — **home passes; the Studio loses its way back**

Home: `scrollWidth == clientWidth == 360`, no horizontal scroll, single-column
grid, the top bar wraps to two rows and all four of its controls are
hit-testable. The composition bar, its legend and its caveat all reflow
correctly. This is a genuinely good small-screen page.

Studio: it stacks into one very long column (the two rails measure 1,438 px and
16,705 px tall) and the `Submit run` footer stays pinned and in view. But 995
elements sit past the right edge with no horizontal page scroll to reach them,
and **the top bar's `← workspace` link and the `PROJECT p_00003 First-user
walkthrough` chip are painted over and unreachable.** (**Defect 4**)

### 14. Stub mode — **passes, and is exemplary**

With every `/api/**` request aborted, the nav chip turned to
`STUB DATA · backend unreachable` and the body read:

> **The workspace call did not come back.** … It then asked `js/mock.js` for the
> workspace instead, and `js/mock.js` carries no `/api/workspace` route, so
> there is no fabricated project list to show either. **That is the better of
> the two outcomes: nothing on this page has been invented to fill the gap.**
> … **PROJECTS** — No project list, because this page could not reach the
> server. This is not "there are no projects". It is "nothing answered".

Three actions offered: `Try again`, `Open the Studio`, `Check server health`.
Contract §6's stub-labelling clause is not merely satisfied; it is the best
error state I saw.

### 15. The empty first install — **passes**

A server with `TMESIM_STATE_DIR` pointed at an empty directory: 0 projects, 0
runs, 0 specimens. The page renders `GETTING STARTED`, three numbered steps
with one button each, and:

* *"A project is where that work is kept. It owns nothing but ids … That is why
  deleting a project never deletes its runs, and why two specimens can point at
  the same run."* — §2 taught at the point of first contact.
* *"Tissue is registered by **filesystem path** on the machine running this
  server. This page never moves image data and never holds any. **7** tissue
  sections are already registered on this server."* — no "upload", and the
  number carries its subject.
* *"**14** compounds are loaded on this server."*
* *"You do not need a project to run the simulator. The Studio works on its own
  exactly as it always has."* — §8's compatibility bar, said out loud.

I could not find the word "upload" anywhere on the home surface, and the
`grep` for it over `home.html`, `js/home/` and `home.css` returns nothing.

### 16. The contract's own checks

* §7 token check (`grep var(--…)` against `tokens.css`): clean, no missing
  tokens, and no colour literal anywhere in `home.css`.
* Dark theme (`prefers-color-scheme: dark`): `data-theme="dark"`,
  `body` background `rgb(17, 20, 24)`, 0 errors.
* Prohibition grep over `home.html` + `js/home/` + `home.css`: the only hits for
  `ranking` and `kill fraction` are the comments in `projectCard.js:12-20` and
  `projectGrid.js:25` that *state the prohibition*.
* The cross-specimen comparison (Studio → `projects`) is headed
  `SAME COMPOUND · NOT COMPARABLE` and states *"No ordering of compounds is
  computed anywhere in this document."* The withdrawn claim is not on this
  surface.
* `tests/test_home_page.py`: **28 passed in 5.04 s** — while Defect 1 was live
  on disk, and after it was fixed. It is green either way; see below.

---

## Part 3 — the defects

### Defect 1 — the home page rendered nothing at all, forever · SEVERITY: CRITICAL · **FIXED 04:05:54 by another agent**

`app/web/js/home/projectGrid.js:156` (as of 03:30:42) built the header button
with

```js
h('button.btn.primary#pgrid-create', { … }, 'New project')
```

`h()` at `app/web/js/util/dom.js:9` matches
`/^([a-zA-Z0-9-]+)?(#[^.\s]+)?((?:\.[^.#\s]+)*)$/` — **the id must precede the
classes** — so that selector does not match and line 10 throws
`bad selector: button.btn.primary#pgrid-create`.

The throw propagates out of `createProjectGrid()`, out of `wire()`, out of
`boot()` in `app/web/js/home/main.js:521` — *after* `stopLoading()` has already
run at line 518 and *before* `render(res)` at line 522. Net effect: the loading
card stays, frozen, with a live progress bar, forever.

**Root cause, and the reason this is worth keeping in the record:** the
docstring of `h()` itself, at `app/web/js/util/dom.js:4`, documented the form it
rejects — `h('div.cls#id', {attrs}, ...children)`. The agent that wrote
`projectGrid.js` copied the documented example. The helper's own documentation
was the defect.

**Reproduction (on `projectGrid.js` as of 03:30:42):**
1. `./.venv/bin/python -m app.server --port 8771`
2. Open `http://127.0.0.1:8771/` in Chromium.
3. Wait as long as you like. `#home-body` stays empty; the console carries one
   `bad selector` pageerror.

**Current state.** Fixed at 04:05:54 — `h('button#pgrid-create.btn.primary', …)`
— and `dom.js:3-16` now carries a block comment naming the trap and the
incident. I re-walked `/` unpatched on the 04:14:50 process: grid renders in
0.23 s, 0 page errors.

**The gap that is NOT fixed.** `tests/test_home_page.py` passed 28/28 through
the entire window this was live. Nothing in it executes the page: it asserts
`GET /` is 200, that the six module files exist, that each parses under
`node --input-type=module --check`, the CSS tokens, and the prohibition greps.
A *syntactically valid* module that throws at call time is invisible to all of
it. `tests/test_first_session.py` already has the Playwright harness that would
catch this in four lines — load `/`, wait for `#projects`, assert
`page.on('pageerror')` collected nothing, assert `#home-body` has children. The
acceptance gate for the workspace has no browser test at all, and this is
precisely the class of failure a parallel build produces.

### Defect 2 — the recent-runs strip cannot show a run you just made · SEVERITY: HIGH · reproduces on current code

`GET /api/workspace` builds `recent_runs` only from run ids held by project
specimens (`app/server/workspace.py`, the `cards` list assembled over
`store.list()`'s specimens). Runs not attached to any project do not appear.

On the journey server: **107 run directories in `app/runs/`**, `totals.runs: 9`
(now 12, after I attached three). The Studio's own run history in the same
product reads `RUN HISTORY 97 runs`.

**Reproduction:**
1. Land on `/`, note the `RECENT RUNS · across every project, newest first`
   strip.
2. Create a project, open it in the Studio, submit a run, wait for `done`.
3. Press `← workspace`. The run is not in the strip. The strip's newest entry
   on this server was 54 days older than the run I had just finished.

This is *consistent with contract §1* as written ("across ALL projects"). It is
not consistent with what the page says, or with what a user means by "my recent
runs". The seam is that `app/runs/` is the store and a project is an index over
it, and the home page shows only the index.

Related, same root: `totals.runs` and the header sentence "N projects on this
server, holding M specimens and K runs" read as server-wide facts ("on this
server") while being project-scoped facts. The code comment at
`app/server/workspace.py:447-451` knows this and says the server total lives at
`GET /api/health`. The page does not.

### Defect 3 — the strip promises something the product does not do · SEVERITY: HIGH · reproduces on current code

`app/web/js/home/recentRuns.js:249-252`, the empty-state copy:

> No runs yet, in any project. **Open a project and submit a run in the Studio;
> it will appear here with the project and specimen it belongs to.**

I did exactly that — opened `p_00003` from the home page, submitted `r_00100`
in the Studio, waited for `done`, returned. It did not appear here. It appears
only after the separate attach offer is accepted, which the copy does not
mention. A new user who dismisses the offer (its own button says `not now`) is
told, by the page, that their run will show up; it never will.

This is the same seam as Defect 2, but it is worse, because it is a written
promise rather than an omission.

### Defect 4 — the Studio's way back to the workspace is unreachable below ~505 px · SEVERITY: MEDIUM · reproduces on current code

`app/web/css/layout.css:130-137` pins `.topbar-group.last` with
`position: sticky; right: 0; z-index: 1; background: var(--bg-2)`, over a bar
that is `overflow-x: auto` (line 78). The sticky group is 335 px wide. The
`← workspace` link occupies content-x 128–204 and the `PROJECT p_00003` chip
sits just after it.

At viewport width W the sticky group covers x = W−335 … W **at every scroll
position**, because it is pinned to the scrollport. For W ≤ ~505 that range
includes x = 128–204, and `scrollLeft` cannot go below 0, so those items can
never be brought out from under it.

**Measured with `document.elementFromPoint()` at the centre of the back link:**

| viewport | topmost element at the back link | clickable |
|---|---|---|
| 360 px | the tour popover, then `.topbar-group.last` | **no** |
| 420 px | — | **no** |
| 460 px | — | **no** |
| 480 px | a button belonging to `.topbar-group.last` | **no** |
| 500 px | `.topbar-group.last` | **no** |
| 510 px | the link itself | yes |
| 720 / 900 / 1100 / 1600 px | the link itself | yes |

At 360 × 740 the rendered bar reads `Vi | account | projects | protocol | tour |
keys | theme: light` — the brand clipped to two letters, and no `← workspace`,
no `PROJECT`, no `SECTION`, no `RUN`, no `CELLS`.

Contract §8.2 requires the Studio to "show the project name in the top bar with
a link back to `/`" when `?project=` is present. Below ~505 px it does neither.
The home page's own nav (`js/home/nav.js`) does **not** have this problem — it
wraps to two rows and all four controls stay hit-testable at 360 px.

### Defect 5 — `totals.runs` double-counts a run shared by two projects · SEVERITY: MEDIUM · reproduces on current code

`app/server/workspace.py:344` — `n_runs_total += int(row.get("n_runs") or 0)`,
a plain sum over projects. A run referenced by two projects is counted twice.

**Reproduction:** on a server with one project holding 9 runs, press
`duplicate` on its card and confirm. The header line changes from
`2 projects on this server, holding 3 specimens and 9 runs.`
to
`3 projects on this server, holding 5 specimens and 18 runs.`
There are still 9 distinct runs, and `app/runs/` is unchanged (measured: 107
directories before and after).

The duplicate dialog, one screen earlier, warns in bold: *"Two projects
agreeing about a run is one run counted twice, not two runs that agree."* The
header then does it.

### Defect 6 — `"1 specimen(s)"` on a card that is otherwise meticulous about plurals · SEVERITY: LOW · reproduces on current code

`app/server/projects.py:958-960`:

```python
return _no_thumbnail(
    f"no run is attached to any of this project's {len(specimens)} "
    "specimen(s), so nothing has been reconstructed yet.")
```

**Reproduction:** create a project with one specimen and no runs; the card's
`NO COMPOSITION YET` placeholder reads *"no run is attached to any of this
project's 1 specimen(s)"*. Everything else on the same card gets this right —
`1 specimen`, `0 runs attached`, `1 section touched`, `2 projects`. There is a
second instance at `app/server/projects.py:1689`.

### Defect 7 — focus is lost after creating a project · SEVERITY: LOW · reproduces on current code

**Reproduction:** tab to `New project`, press Enter, press Enter again on
`Create project`. `document.activeElement` is `<body>`. A keyboard user is
returned to the top of the document with no announcement that a project was
made.

Cause: `main.js`'s `refresh()` re-renders the whole grid, destroying the
`#pgrid-create` button that native `<dialog>` would otherwise restore focus to.
The `Escape` path is correct — it returns focus to `#pgrid-create` — so only
the success path is affected, which is the one that matters.

### Defect 8 — the create dialog hides the section qualifiers the Studio shows · SEVERITY: LOW-MEDIUM · reproduces on current code

`app/web/js/home/createProject.js:554-556` renders each option as
`` `${x.id} — ${x.label}` ``. `GET /api/sections` also carries `caveats`,
`n_caveats_blocking`, `spatial_spread_valid` and `reason`, which the Studio's
own picker uses to append `[ingested — not validated]` and `[N refused]`.

**Reproduction:** open `New project`, expand `FIRST SPECIMEN`. 38 sections, all
presented identically. Then open the Studio and look at the same list: 34 of
them carry `[ingested — not validated]` and `[4 refused]`. The one place a
first-time user makes this choice is the one place the qualifiers are missing.

### Defect 9 — the delete arming is announced to nobody · SEVERITY: LOW · reproduces on current code

The modal's `Delete project` button requires two presses: the first relabels it
to `delete p_00004?`, the second issues `DELETE /api/projects/p_00004`. I
confirmed there is no `aria-live` region anywhere in the open dialog
(`[...dialog[open].querySelectorAll('[aria-live]')]` → `[]`), so for a screen
reader the first press produces silence and the second destroys the project.

### Defect 10 — `renderLadder is not defined` · SEVERITY: unknown · TRANSIENT, observed during a live edit

At ~03:49, immediately after `r_00100` finished, the page threw
`renderLadder is not defined` **twice**. The visible consequences: the READ
rail header read `showing no run` while displaying a full result for `r_00100`;
the top bar's `CELLS` read `—` instead of `40,887`; and the stage showed
`NOTHING ON THE STAGE YET · To put something here: choose a section and a
compound in the left rail and press Submit run` — to a user who had just done
precisely that.

`app/web/js/ui/reach.js` has mtime 03:51:35 — sixteen seconds after I captured
the error, and it was rewritten again at 04:06:01. The function is defined at
line 1364 and called at line 1807, both at module scope; `node --input-type=
module --check` parses it clean. **It does not reproduce on current code** — I
reloaded and re-selected `r_00100`, 0 page errors, header `showing r_00100`,
`CELLS 40,887`, stage rendering the field with its time-lapse scrubber.

I record it because it is what a real user would have seen in that window, and
because its blast radius is instructive: one `ReferenceError` in the result
renderer made the product claim nothing had been run.

### Defect 11 — one press of `attach run` files a run you did not choose · SEVERITY: LOW · reproduces on current code

In the Studio's `projects` modal, each specimen row carries a `<select>` of
candidate runs and a bare `<button class="btn sm">attach run</button>`. The
button commits the select's current value with no confirmation. With the Studio
showing `r_00102` in its top bar and the select left at its default, one press
issued `POST /api/projects/p_00003/specimens/s_00001/runs` with body
`{"run_id":"r_00100"}` — a different run from the one on screen.

This is a picker-plus-commit pair behaving as designed, and it is reversible:
each attached run carries an `×` to detach. But the commit is one press, the
default is a real run, and nothing says "this is about to file r_00100".

---

## Part 4 — confusing rather than broken

These cost me more time than several of the defects above.

1. **Two numbers for one machine.** The home page says `9 runs`; the Studio in
   the next tab says `RUN HISTORY 97 runs`. Neither is wrong. Together they are
   the single most disorienting thing in the product, and I had to read
   `workspace.py` to work out which question each was answering.

2. **The delete needs two presses of the same button.** I recorded "delete does
   not work" and only found the arming by dumping the button's label between
   clicks. Inside a modal that is *already* a confirmation, a second press with
   no instruction reads as a dead button.

3. **`SETUP — what will run` versus `READ — showing r_00100`.** When I selected
   an older run, the left rail kept `dx = 12` while the run I was reading used
   `dx = 30`. That is correct and deliberate — the two headers say so, and the
   COMPARISON dock spelled out `SPEC DIFFERENCE — 1 OF 13 FIELDS · dx_um 12.0
   um | 30.0 um`. It still took me a minute to trust that the rail was not
   simply wrong.

4. **The run-history table loses its identity column.** After the dock scrolled
   horizontally to `AUC / PEN / RUNTIME / SUBMITTED`, the `RUN` column was gone
   and five rows of `4.226e6 / 278.5 / 41s` were indistinguishable from one
   another.

5. **A newcomer landing on a populated server gets no orientation.** The
   excellent `GETTING STARTED` block — what Vivum does, what a project is, that
   deleting never deletes runs — renders only when `projects` is empty. Land on
   a colleague's server with one project and none of it exists. The Studio has a
   `tour` button in its top bar; the home page has no equivalent.

6. **Body text is 11.5 px against a 16 px root.** The workspace home is a
   document and reads like an instrument panel. At 200 % zoom it is excellent;
   at 100 % the primary call to action is 81 × 22 px.

7. **A first-run tour that covers the whole Studio.** Well written, honest,
   once-only — but it lands on top of the screen you were sent there to use,
   and `Skip tour` is bottom-right, the last place a first-timer looks.

8. **`Kill fraction` is the default metric of the cross-specimen comparison.**
   The comparison itself is correctly fenced (`SAME COMPOUND · NOT COMPARABLE`,
   "No ordering of compounds is computed anywhere in this document"), and the
   left rail warns that absolute kill must not be quoted from this section under
   any state route. Having it pre-selected in a dropdown is not a violation. It
   is the one place where the default and the warning point in opposite
   directions.

---

## Part 5 — the three things I would fix first

**1. Put a browser in the acceptance gate.**

Not because of Defect 1 specifically — that is already fixed — but because the
gate was green through every minute of it, and would be green through the next
one. `tests/test_home_page.py` cannot distinguish "the workspace home renders
the truth the server sent" from "the workspace home renders nothing and looks
like it is still loading", and those are the two states the §5 file-layout test
was written to tell apart. `tests/test_first_session.py` already carries the
Playwright harness. The test is four assertions: navigate to `/`, wait for
`#projects`, assert `page.on('pageerror')` collected nothing, assert
`#home-body` has children. Add the same for `/studio`.

This is first because it is the only item on this list that changes what the
*next* parallel build can ship.

**2. Close the run/project seam — Defects 2, 3 and 5, which are one defect.**

A user's own work must be visible on the page that claims to hold their work.
Three pieces, in order of value:

* Fix the written promise at `app/web/js/home/recentRuns.js:249-252` first —
  it is one sentence and it is currently false.
* Make `totals.runs` count distinct run ids rather than summing per-project
  `n_runs` (`app/server/workspace.py:344`). The duplicate dialog already
  explains why.
* Then decide what `RECENT RUNS · across every project` should mean. Either it
  covers `app/runs/` and unattached runs carry an honest "not filed under any
  project" marker with the attach affordance right there, or the strip is
  renamed to what it is and the page states the server's own run total beside
  it. Do not invent a placeholder project id to paper over it — an intermediate
  version of `workspace.py` I saw at 03:54 appended unattached runs as
  `project_id: "p_00000", project_name: "unknown project"`, which would have put
  two fabricated ids on the home page. That version was reverted; it should stay
  reverted.

**3. Fix the Studio top bar below ~505 px — Defect 4.**

`app/web/css/layout.css:130-137`. A sticky opaque group that permanently covers
the left end of its own scroll container is a bug in the pattern, not in the
width. The identity of what you are looking at (`PROJECT …`) and the way out
(`← workspace`) are the two things that must survive every viewport; the
protocol/tour/keys/theme launcher is the thing that can afford to scroll. At
360 × 740 the product currently has no visible way back to the workspace, which
is the one navigation the whole contract exists to create.

---

## Appendix — state I left behind

Created through the product, on the journey server's store:

* **`p_00003` "First-user walkthrough"** — 1 specimen (`s_00001`, xenium_319),
  with `r_00100`, `r_00101`, `r_00102` attached. Two of those three attachments
  were made by the one-press `attach run` button while I was probing Defect 11
  and are not runs I would have chosen. Deleting the project detaches them and
  destroys nothing else; I verified that behaviour directly.
* **`r_00100`–`r_00102`** in `app/runs/` — three real runs on xenium_319
  (temozolomide dx = 30, cisplatin dx = 30, doxorubicin dx = 12). `app/runs/`
  went from 104 to 107 directories.
* `p_00004`, the duplicate of `p_00002`, was created and deleted; `p_00002` is
  byte-for-byte as I found it, all 9 run references intact.

I modified no file in the repository. This document is the only thing I wrote.
