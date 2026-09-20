# Workspace contract — the landing page, and what the Studio becomes

**Frozen 2026-09-20.** Agents build against this. If you believe a clause is
wrong, REPORT it rather than diverging — a second interpretation of a shared
contract is how five parallel agents produce four incompatible halves.

---

## 0. The problem this solves

Today `/` serves the Studio: a three-rail workstation that assumes you already
know which section you want, which compound, and what a RunSpec is. There is no
place to *keep* work. `app/server/projects.py` has held a real project store
for weeks — 11 routes, specimens, runs, cross-specimen comparison — and the only
way to reach it is a modal button in a toolbar inside the Studio you must first
know how to drive.

So the product opens on its hardest screen, and its organising concept is
hidden behind it.

**After this contract:** `/` is a workspace home. You land on your projects,
open one, and the Studio opens *inside* that project with its runs already
attached. The simulation stays the centre; it stops being the front door.

---

## 1. Routes (server)

### Pages

| Path | Serves | Notes |
|---|---|---|
| `/` | `app/web/home.html` | the workspace home. **Ungated**, like every static asset — it must load to prompt for a key. |
| `/studio` | `app/web/index.html` | the Studio, unchanged in function |
| `/studio?project=p_00002` | same | Studio scoped to a project |
| `/ingest` | unchanged | bring-your-own-tissue page |

`index.html` KEEPS working at its old URL. Anything bookmarked stays alive.

### `GET /api/workspace` — the home page's single boot call

One call, because the home page needs four things and four round-trips on a
cold Render instance is four chances to look broken.

```
{
  "projects": [ProjectSummary],       // newest-updated first
  "recent_runs": [RunCard],           // across ALL projects, newest first, <= 12
  "sections_available": int,          // how many sections this server can run
  "drugs_available": int,
  "totals": {"projects": int, "runs": int, "specimens": int},
  "server": {"contract": "v1", "stub": false}
}
```

`ProjectSummary` is exactly what `GET /api/projects` already returns, plus:

```
  "last_opened_at": float | null,     // NEW, see §3
  "thumbnail": {...} | null           // NEW, see §4
```

`RunCard` is a *flattened* view for the home page, and it is deliberately not
a `RunResult`:

```
{"run_id", "project_id", "project_name", "specimen_id", "specimen_label",
 "section", "drug", "status", "created_at", "finished_at", "runtime_s"}
```

### New project routes

| Method | Path | Body → Returns |
|---|---|---|
| `POST` | `/api/projects/{id}/open` | `{}` → `{ok, last_opened_at}` — stamps recency. Idempotent. |
| `POST` | `/api/projects/{id}/duplicate` | `{name?}` → `Project` (201). Copies specimens and their run ids. **Runs are NOT re-executed and NOT copied** — a run belongs to the store, and two specimens may reference the same run id. |

Everything else already exists. Do not add a route you have not first checked
for in `app/server/api.py`.

---

## 2. What a project is, restated

A project **owns nothing but ids**. Specimens hold `run_id`s into `app/runs/`,
which stays the single store for run results and binary payloads. Attaching,
detaching and deleting a specimen therefore never touches a run. That is
already true in `projects.py` and this contract does not change it.

The consequence, which the UI must not hide: **deleting a project does not
delete its runs.** The home page says so at the point of deletion, not in a
help file.

---

## 3. `last_opened_at`

A separate field from `updated_at`, and the distinction is load-bearing:

* `updated_at` — the project DOCUMENT changed (a specimen was added, a run
  attached, the notes edited).
* `last_opened_at` — a human looked at it.

Sorting a workspace by `updated_at` puts a project at the top because a
background job attached a run to it. Sorting by `last_opened_at` puts it at the
top because you were working on it. The home page defaults to
**`last_opened_at`, falling back to `updated_at` when it is null** (every
project that predates this field).

`null` is rendered as "not opened yet", never as an old date and never as
"now".

---

## 4. `thumbnail`

A project's visual identity on a card, and it must never be a *fabricated*
picture of tissue.

```
{"kind": "composition", "section": str, "n_cells": int,
 "type_fractions": {label: float}}   // from a real reconstruction
```
or
```
{"kind": "none", "reason": str}      // no run has finished yet
```

The card draws `type_fractions` as a small stacked bar. That is a chart of
measured composition. It is NOT a rendering of the tissue, and the card never
implies it is. `{"kind": "none"}` renders as a neutral placeholder with its
reason on hover — not a grey rectangle that looks like a failed image.

---

## 5. Front-end file layout

```
app/web/home.html              the landing page (NEW)
app/web/css/home.css           its stylesheet (NEW)
app/web/js/home/
    main.js                    boot: fetch /api/workspace, render, wire
    projectCard.js             one card, incl. the composition bar
    projectGrid.js             the grid, search, sort, empty state
    createProject.js           the create/rename/duplicate/delete dialogs
    recentRuns.js              the cross-project recent-runs strip
    nav.js                     the top bar shared with the Studio
app/web/js/util/               SHARED with the Studio — h(), fmt, etc.
app/web/js/api.js              SHARED — add workspace calls here, do not fork
```

**`api.js` is shared and must not be duplicated.** The home page imports the
same module the Studio does; that is what keeps the API key, the stub fallback
and the error shapes identical on both pages.

---

## 6. Design rules the home page inherits

These are not new. They are the Studio's existing rules, and the home page is
the same product.

* **A caveat is never set smaller than the number it qualifies.**
* **Every number carries its unit and its denominator.** "9 runs" not "9".
* **A refusal says what to do next.** Never a bare dash, never a spinner that
  is really an error.
* **No cross-compound ordering, in any form.** Withdrawn 2026-08-21. A project
  card lists the compounds a project has run — as a SET, alphabetical, never
  ordered by any metric, and never with a "best" badge.
* **No absolute kill fraction as a headline.** Not on a card, not in a stat.
* **Tissue registers BY FILESYSTEM PATH.** The word "upload" never appears
  about tissue.
* **Stub mode is labelled.** If `/api/workspace` falls back to mock data the
  home page says so, loudly, exactly as the Studio does.

---

## 7. Theming, and the one hard requirement

`app/web/css/tokens.css` is the palette. **Use its tokens; define no new
colour literal.** The home page must work in both themes and must respect the
`data-theme` attribute the Studio's toggle already sets — the two pages share
`localStorage`, so a user who picked light in the Studio lands on a light home
page.

Verify with the token check that already exists:
```
grep -oE "var\(--[a-z0-9-]+\)" app/web/css/home.css | sort -u \
  | sed 's/var(\(.*\))/\1/' | while read t; do
      grep -q -- "$t:" app/web/css/tokens.css || echo "MISSING $t"; done
```

---

## 8. How the two pages hand off

The Studio reads `?project=<id>` at boot. When present:

1. it calls `POST /api/projects/{id}/open`,
2. it shows the project name in the top bar with a link back to `/`,
3. a finished run is offered for attachment to a specimen in THAT project
   — offered, not attached silently.

With no `?project=`, the Studio behaves exactly as it does today. **A user who
never visits the home page notices no change.** That is the compatibility bar
for this whole contract.

---

## 9. What this contract does NOT authorise

* Rewriting the Studio's three-rail layout. Improve it; do not replace it.
* A new API-key mechanism. There is one shared key; the home page uses the
  same gate and the same `apikey.js`.
* Sign-up, user accounts, or per-user projects. `GET /api/account` already
  states the truth: one shared key, no identities, no attribution. The home
  page must not imply otherwise — no avatars, no "your" projects, no owner
  column.
* Deleting or rewriting anything in `results/` or `archive/`.
