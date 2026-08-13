# Bring your own tissue

*What a researcher must do to run Vivum on their own 10x export, what the validator
catches, and what it refuses.*

Every section in `data_cache/` was prepared by this project. A simulator a
researcher cannot point at their own sample is a demo, not a product. This
documents the path from "a fresh 10x export sits on my disk" to "I can run it."

---

## 1. The honest starting position: most of this already existed

The first job was to establish what was already there rather than build a second
path beside it. The finding: **the ingestion library was already complete and
good.** `tmesim/data/ingest.py` (1,812 lines) already did platform detection,
filename resolution, validation with named bounds, and a capability report;
`app/server/ingest_registry.py` already turned a directory into a section record
that `cache.all_sections()` appends after the seven built-ins; `/api/sections/scan`
already exposed both. None of that was rebuilt.

What was actually missing was narrower, and is what this change adds:

| Gap | Status before | Status now |
| --- | --- | --- |
| A command a user can run | **None.** `grep -c ingest tmesim/cli.py` → 0. The only entry points were the Python API and an HTTP POST to a running Studio. | `vivum ingest <dir>` |
| The spot-data spatial refusal at ingest time | **Absent from the capability report.** `validation/immune_sink_impact.py` refuses to compute the statistic on spot data, but nothing told a user their Visium sample had that limitation *before* a run. | `spatial_spread` capability, refused on every spot section |
| Caveats visible where a section is chosen | Ingested sections were already selectable, but `detect_issues`, `runnable` and the stored `capability` were rendered **nowhere** in the web UI. `getCapability` in `app/web/js/api.js` was a dead export with zero callers. | `caveats` on `GET /api/sections`, shown in the picker |
| Tests against an export this repo has never seen | Every fixture in `tests/test_ingest.py` symlinked real cached bytes, so all of them **skip** without `data_cache/`. | 20 tests built from synthetic bytes |

---

## 2. What a user must do

**Before:** write Python against `tmesim.data.ingest`, or start the Studio and
`POST /api/sections/scan` with a JSON body. Neither is documented in
`docs/QUICKSTART.md`. There was no way to find out what a sample could support
without writing code.

**Now:**

```bash
# Read-only: what is this directory, and what can it NOT support?
vivum ingest /data/PatientX_Xenium

# Same, but load the matrix and measure real detection rates
vivum ingest /data/PatientX_Xenium --deep

# Make it selectable in the Studio and the other subcommands
vivum ingest /data/PatientX_Xenium --register --id patient_x

vivum ingest --list
vivum ingest --remove patient_x
```

The default is read-only and registers nothing. `--json` emits pure JSON on
stdout (no progress line) so it can gate a data-prep script. Exit codes: `0` ok,
`1` refused, `2` usage error.

Accepted layouts, all detected on **structure**, not on this repo's filenames:

- **Xenium** — `cells.csv.gz` + `cell_feature_matrix/` (or the `.tar.gz`), `gene_panel.json` optional.
- **Visium** — a Space Ranger output: `spatial/` + a MEX matrix directory.
- **Bare MEX** — `matrix.mtx.gz` + `features.tsv.gz` + `barcodes.tsv.gz` beside a coordinates table. Ingestable and inspectable, but **registered as INSPECT-ONLY**: `reconstruct_real` has no route for it, and the registry says so rather than failing at run time.

---

## 3. What the validator catches

Every refusal names the file, column or bound that failed, with the measured
value beside it. None of these is new logic invented here except where noted;
what is new is that a user can reach them.

| Condition | What the user is told |
| --- | --- |
| Path is not a directory | Names the layouts it expects, per platform. |
| No recognised layout | Lists all three accepted layouts and their required files. |
| `cell_feature_matrix/` present but short a file | **New.** Names the *missing file* (`missing matrix.mtx[.gz]`) instead of the previous "no cell-feature matrix", which was misleading when the directory was plainly there. This is the commonest real failure — an interrupted download. |
| Matrix still in a tarball | Says it will be extracted in place on first load, and roughly how long. |
| Coordinates in pixels or nanometres | Refused with the measured extent against the 10 cm bound: *"extent 238808.3 x 238331.1 um exceeds 100000 um (10 cm), which is larger than any histology section."* |
| Coordinates possibly in millimetres | Warned, with the 5,100–24,000 um span of this repo's sections as the comparison. |
| Fewer than `MIN_CELLS` (50) | Refused with the count and the floor. |
| Matrix not aligned to coordinates | Refused: *"expression has N rows but there are M cells … every cell would carry another cell's transcripts."* |
| All-zero or negative counts | Refused (raw counts are required; the pipeline log1p-transforms). |
| Incomplete coordinate join (MEX) | Refused rather than silently truncated. |
| Duplicate gene symbols | Warned, stating that marker lookup takes the first column. |

Registration is **refused as a whole, never degraded**: a sample that fails
validation registers nothing. A built-in section id can never be shadowed, so
nothing an ingested section does can change what a published number was measured
from.

---

## 4. What the ingested section IS — and what it cannot support

`vivum ingest` reports cell count, panel size, physical extent, coordinate units,
cell density and median counts per cell, then the capability report: which cell
types have **no marker** on the panel, and what that costs.

The point is that a user's panel will differ from ours, and the consequences are
stated **up front** rather than discovered in a caveat afterwards:

- **No endothelial marker** → `vascular_field: NO`. `vascular_density` would be identically zero, the oxygen PDE would have no source. *This sample cannot support any oxygen or delivery result at all.*
- **No hypoxia marker** → `hypoxia_measured: NO`, and consequently `absolute_kill: NO`. This is the caveat the shipped GBM 319-gene panel already carries — it measures no hypoxia — and it now applies automatically to a user's panel by the same rule rather than being a hardcoded note about ours.
- **No marker map fits** → stated with the best coverage found, rather than silently falling back to a map that does not describe the tissue.
- **No ECM transcript** → the barrier term becomes "a proxy for a proxy", and says so.
- **No MGMT** → TMZ resistance is unmeasurable; every TMZ kill is a no-resistance upper bound.

Without `--deep` only panel **presence** is read, and every rate-dependent verdict
is reported `UNMEASURED` rather than guessed. That distinction is preserved
deliberately: a gene being on a panel does not mean it fires often enough for a
per-cell call.

---

## 5. What it refuses, and why — the spot-data spatial statistic

This is the hard line.

`reconstruct/real_visium.py::_place_cells` does **not** measure where the cells
inside a Visium spot sit. It synthesises `_CELLS_PER_SPOT = 12` pseudo-cells per
spot at the spot centroid with uniform jitter inside a spot-sized disc. Any
statistic computed over the per-cell spatial distribution of such a section
therefore measures **that jitter**, not the tissue.

`results/immune_sink_impact.md` §5.6 already declared the per-cell attenuation
spread (p90/p10) and the immune-neighbourhood structure ratio **INVALID on
Visium/spot data**, and `validation/immune_sink_impact.py` already refuses to
compute them there. What was missing is that nothing said so at ingest time — a
user pointing Vivum at their own Visium export had no way to learn this before a
run.

So `capability_report` now carries a `spatial_spread` capability, keyed on the
**platform** and not on the panel:

- `platform == "visium"` → `available: False`, evidence `platform`. The reason states the 12 synthesised pseudo-cells, cites §5.6, and — importantly — states what *is* still defensible: **the cross-drug ranking on that section is unaffected.**
- single-cell platforms → `available: True`, with the standing caveat that the statistic must be read against a spatially-aware null, not against a well-mixed 1.0 (§5.1).

A test asserts the refusal text quotes the pseudo-cell count that
`real_visium._CELLS_PER_SPOT` actually uses, so the report cannot go stale if that
constant changes.

The refusal reaches the Studio through `GET /api/sections`, which now carries
`spatial_spread_valid` and a `caveats` list for **every** section, built-in or
ingested. All 36 spot sections currently registered (4 built-in Visium + 32
ingested GBM-Space) report `spatial_spread_valid: false`. The Studio does not
compute or display a spatial spread anywhere today, so there was no existing
number to withdraw; this ensures it cannot start reporting one on spot data
without tripping over the flag.

---

## 6. Can the Studio select an ingested section?

**Yes, and it already could** — `cache.all_sections()` appends ingested sections
after the built-ins, so they flow through `GET /api/sections` into the picker in
`app/web/js/ui/spec.js`. That was pre-existing and is unchanged.

What is new is that the picker no longer shows only `id`/`label`/`available`:

- each option carries a `title` tooltip with the panel, the unavailability reason, and every caveat;
- an option with blocking caveats is suffixed `[N refused]`;
- selecting a section paints a caveat line naming each refused readout beneath the picker.

Two limitations, stated rather than papered over:

1. **There is still no UI to register a directory.** Registration is CLI or `POST /api/sections/scan`. Adding a directory picker to the web UI would mean a path-entry surface in a panel owned by another agent this session, so it was left out; the CLI covers it.
2. **The Studio picks up sections registered by another process only on restart.** `ingest_registry.load_store()` reads once per process by design. The CLI prints this after a successful `--register`.

---

## 7. Tests

`tests/test_ingest.py` gains 20 tests built from a **wholly synthetic** export —
random counts, random coordinates, written by `_write_synthetic_xenium` — rather
than symlinked cache bytes. They run in a checkout with no `data_cache/` at all,
which the pre-existing fixtures could not.

They invent bytes but never a result: nothing asserts a scientific number, because
measuring biology from random counts would be meaningless. What they assert is
what the ingester **says** about an unfamiliar directory — which platform it
detected, which named file or column it reports missing, which capabilities it
refuses, and that a registered synthetic section becomes selectable while a
built-in one stays untouched.

Two of them failed on first run and caught real defects, both fixed:

1. `--json` printed a `[ingest] scanning …` progress line before the JSON, making the output unparseable for the scripts that flag exists for.
2. A `cell_feature_matrix/` directory missing one file reported "no cell-feature matrix", naming the directory rather than the absent file.
