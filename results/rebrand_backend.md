# Rebrand: Vivum → Codon (engine, server, docs, CLI)

**What this file records.** The product is now called **Codon**. This is the
backend half of that rename — the Python engine, the Studio server, the docs,
the README and the CLI. The landing page (`landing-3d/`) and the Studio front
end (`app/web/`) were renamed separately by other workstreams and are not
described here.

It also records the things that were **deliberately not renamed**, because on a
rename the omissions are the part that needs justifying, and an undocumented
omission reads later as an oversight rather than a decision.

---

## 1. What is now Codon

| Surface | Change |
|---|---|
| CLI | `prog` is `codon`; every `--help` string, every printed hint (`codon ingest <dir> --register`, `codon studio`, `codon compare-drugs …`), and the scorecard banner (`CODON EXTERNAL-ANCHOR SCORECARD`) |
| Console script | `codon` is the primary entry point; `vivum` kept as a working alias (§2) |
| Studio server | FastAPI `title` (`Codon Studio API`, shown in `/docs`), the browser `<title>`, the backend-up placeholder page, and the startup banner `codon studio` prints |
| `README.md` | Title and all prose and commands |
| `docs/` | `QUICKSTART.md`, `DEPLOYMENT.md`, `UI_ARCHITECTURE.md`, `molecular_layer.md`, `real_data_status.md` — prose, commands, and the transcripts of CLI output quoted in them |
| Docker | `Dockerfile` header and the documented image tag (`codon-studio`), and the CI image/container names in `.github/workflows/docker.yml` |
| Tests | Assertions on user-facing strings updated to the new string, not deleted (`tests/test_cli.py`, `tests/test_ingest.py`, `tests/test_validate_cli.py`) |

The quoted `usage:` block in `docs/QUICKSTART.md` still lines up: `codon` and
`vivum` are both five characters, so argparse's continuation indent is unchanged.

## 2. The CLI: `codon` is the name, `vivum` still works

`pyproject.toml` registers **both** names against the same entry point:

```toml
[project.scripts]
codon = "tmesim.cli:main"
vivum = "tmesim.cli:main"
```

**Why the alias.** The `results/` directory holds ~200 dated files, and the
reproduce instructions in them are written as `vivum …`. Dropping the old entry
point would have turned every one of those into a command that no longer runs —
which is a worse outcome than one extra line in `[project.scripts]`. The record
stays executable.

The alias is the *same function*, not a shim: there is no deprecation warning
and no behavioural difference. `--help` prints `codon` under either name, so the
alias tells you the new name the first time you use it. `tests/test_cli.py`
pins both entries so neither can be dropped by accident.

**When to remove it.** When nothing in `results/` needs to be runnable as
written — not before.

## 3. What was deliberately NOT renamed

**The Python package is still `tmesim`.** `import tmesim`, the distribution
name, and every module path are untouched. Renaming it would rewrite the imports
in ~1000 tests, the reproduce line of every results file, and the frozen
reference JSONs, in exchange for changing a name a user of the product never
types. `docs/QUICKSTART.md` and `README.md` say "Codon (package name `tmesim`)"
where the distinction matters. This was the right call and it is not close.

**`results/*.md` is untouched.** Those are dated findings. Rewriting the product
name inside them would misrepresent what the project was called when each was
written, and this project's whole discipline is that the record says what was
true when it was measured. Nothing in them is broken by leaving them: the
commands they quote still run, because of the alias in §2. **No file in
`results/` was judged to be current product documentation rather than a dated
finding, so the stated exception was not used.** The two closest calls, and why
they were left alone:

* `results/CLAIM_INVENTORY.md` — reads as a living claims register, but it is
  stamped `Compiled 2026-08-12` and its own precedence rule is *check the date
  first*. It participates in the dated record rather than standing outside it.
* `results/bring_your_own_tissue.md` — reads partly as a how-to, but §1 onward
  is a narrative of work done. Its `vivum ingest` commands still run.

This file is the pointer that resolves the ambiguity: **"Vivum" in any
`results/` file is this product, under its previous name.**

**Internal identifiers are unchanged**, deliberately, because renaming them
breaks running deployments for no user-visible gain: the `TMESIM_STUDIO_API_KEY`
environment variable, the `tmesim-warm-drugs` thread name, and
`app/server/.runs` / `_ingested_sections.json` on disk.

**Not mine, still says Vivum** (owned elsewhere, listed so the next person can
finish the sweep): `app/web/js/ui/welcome.js`, `scripts/build_studio_export.py`,
`app/README.md`, `API_CONTRACT.md`, and `experiments/kill_calibration_breast_report.py`
(which emits a `vivum validate` reference into the report it generates).

## 4. Bring your own tissue is now reachable without `curl`

`results/first_user_readiness.md` §"Can a user bring their own tissue today?"
found the ingest path real, generic and good — and exposed **nowhere** a user
could click:

> The gap is exposure, not capability. `sections/scan` appears **nowhere** in
> `app/web/` — the Studio's section picker is a read-only list from
> `GET /api/sections`.

**This was closable from the server side, and is closed.** `GET /ingest` serves
a self-contained page — its own CSS and script, no external host, no dependency
on `app/web/` having been built — with a field for the sample directory,
optional section id and label, an optional `X-API-Key`, a *measure detection
rates* checkbox, and two buttons: **scan (read-only)** and **scan + register**.

It is a thin surface over the existing endpoint, on purpose:

* It makes **one** `fetch`, to `POST /api/sections/scan` — the same call `curl`
  would make. No parallel implementation, no new capability, no new trust
  boundary. The `X-API-Key` gate applies to it exactly as before.
* It renders the returned `capability_text` **verbatim** — the exact report
  `codon ingest` prints — rather than re-deriving a summary. This surface
  therefore cannot disagree with the CLI about what a sample can support, which
  matters because that report is where a user is told, for example, that
  absolute kill is not quotable from their panel.
* An unusable directory is answered with its reason, never a 500 (already true
  of the endpoint; pinned by a test here).
* A registered section is marked *inspectable only* on the page when the
  platform has no reconstruction route, rather than appearing runnable.

Discoverability, which was the actual failure: `codon studio` now prints the
page's URL at startup (`your tissue  http://…/ingest`), and the backend-up
placeholder page links it.

Registered sections appear in `GET /api/sections` immediately, so a Studio
reload picks them up in the existing picker.

**This is a stopgap, not the destination.** The right end state is the Studio's
own section picker learning to ingest — which is front-end work in `app/web/`
and was not done here. What that needs is specified in §5.

## 5. What the Studio front end still needs (not done — `app/web/` is owned elsewhere)

The backend needs no further change for any of this; every endpoint below
already exists and is already tested.

1. **An "Add your own tissue" affordance in the section picker.** A path input
   plus a button, calling `POST /api/sections/scan` with
   `{path, register: true, section_id?, label?}`. On `ok: true`, re-fetch
   `GET /api/sections` and select the new id.
2. **Show the capability report before the run, verbatim.** The response carries
   `capability_text`, preformatted. It must be displayed, not summarised — it is
   the only place a user is told what their sample *cannot* support.
3. **A read-only scan first.** Same call with `register: false` (add
   `measure_detection_rates: true` for measured firing rates, which costs a full
   ingest). Lets a user inspect before committing.
4. **Render the failure path as an answer.** `ok: false` / `error` is a reason,
   not an exception — show the text. Same for `detection.issues` on a scan.
5. **Honour `runnable: false`.** A registered `mex` section is inspectable but
   has no reconstruction route; the picker must not offer it as runnable.
6. **Removal.** `DELETE /api/sections/{id}` — ingested sections only; built-ins
   are refused server-side.

Until that exists, `GET /ingest` is the working path for a field tester.

## 6. Physics is untouched, and what the suite says

This was a strings-only change: no numerical code was modified.

`calibrate_alpha()` returns **22.223477711172563** — exact, verified after all
edits, unchanged from before them.

**Tests: 963 passed, 10 failed.** All 10 failures are the concurrent
area-conserving vessel-painting promotion (`tmesim/physics/subgrid.py`,
`tmesim/reconstruct/vessel_network.py`), not this rename: a full run taken
*before* any edit here produced the same 10. They are `test_vessel_source` (3),
`test_vessel_network` (3), the three bitwise stored-reference tests in
`test_time_resolved_pk`, and `test_app_server::TestValidate::
test_cfl_is_flagged_once_a_speed_has_been_measured`. The root failure,
`test_binary_default_is_untouched`, is a guard that exists to catch exactly that
default flipping — it is doing its job — and is being re-pinned to the new
contract by the workstream that owns it.

Six tests were added here (957 passing before, 963 after): five for the `/ingest`
page and the endpoint behind it, and one pinning **both** console-script entries
so the `vivum` alias cannot be dropped silently.
