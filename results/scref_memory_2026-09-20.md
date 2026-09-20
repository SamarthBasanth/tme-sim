# The scRNA reference, the grid ceiling that cannot catch it, and what 512 MB actually buys

**Date:** 2026-09-20
**Scope:** `tmesim/data/real/scref_gbm.py::load_reference` and the Visium load
path that calls it (`tmesim/reconstruct/deconvolve.py`).
**Status of the headline claim:** the reported defect is REAL and the fix is
PARTIAL — the loader is 3.2x cheaper and bit-exact, and it does **not** make a
Visium section runnable on a 512 MB instance. A refusal was added because the
measurement says a fix cannot get there.

---

## 0. Method, and what the numbers are numbers OF

Peak memory is `resource.getrusage(RUSAGE_SELF).ru_maxrss`.

> **`ru_maxrss` units are platform-dependent.** On **darwin/macOS it is
> BYTES**; on Linux it is KiB. Every figure below was taken on darwin and
> divided by 1024² once, not twice. `scref_gbm._maxrss_bytes()` carries the
> branch, and `tests/test_scref_memory.py` asserts both halves of it.

`ru_maxrss` is a high-water mark for the whole process and never decreases, so
"peak" here means `ru_maxrss_after_the_call − RSS_immediately_before_it`.
Current RSS came from `ps -o rss= -p <pid>` (KiB on macOS). `psutil` is not
installed in this venv.

`tracemalloc` was used once, to check whether it even sees these allocations:
a 47.7 MiB numpy array registered as 47.7 MiB traced, so numpy IS visible to
it on this build. It was then **abandoned for the headline numbers**, because
enabling it inflates the very quantity being measured (with tracemalloc on,
the same load reported a 1,899 MiB peak against 1,682 MiB with it off). Where
a tracemalloc figure appears below it is labelled as such.

Host: darwin 24.6.0, 16 GiB, CPython 3.9, pandas 2.3.3, numpy 2.0.2, scipy
1.13.1. Every run is a **cold subprocess** unless it says otherwise. Scripts
are in the session scratchpad; each is ~20 lines and reproduced by the recipe
in each section.

---

## 1. The reported figure, re-measured

`results/performance_profile_2026-09-20.md` §6 #2 reports **+1,922.6 MiB** in
the pipeline and **+1,612.5 MiB** for `load_reference()` alone.

Re-measured here, `load_reference()` alone, five cold runs:

| run | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| peak over baseline (MiB) | 1586.3 | 1612.8 | 1850.2 | 1960.6 | 1599.5 |
| wall (s) | 2.04 | 2.06 | 2.07 | 2.08 | 2.07 |

**The reported figure reproduces.** The standalone 1,612.5 MiB sits inside the
measured spread and matches run 2 to 0.02 %. The spread itself is the
allocator not returning freed arenas, which is also visible directly: calling
`load_reference()` three times in one process leaves RSS at **377 → 555 → 852
MiB** after each result is deleted and garbage-collected. Repeat loads do not
merely re-pay the transient, they **ratchet the floor**.

The object produced:

| | |
|---|---|
| file on disk | 20.5 MB gzip |
| frame | 23,460 genes x 3,589 cells, **dense int64** |
| `memory_usage(deep=True)` | **643.8 MiB** |
| nonzero fraction | 0.1100 |
| max count | 77,122 (so int32 is exact) |

---

## 2. Why `TMESIM_MAX_GRID_CELLS` structurally cannot catch this

`load_reference()` takes no grid argument and reads no grid state. Its cost is
a function of the reference panel alone, so it is **identical at every
`dx_um`**. `_resolve_max_grid_cells` counts voxels. A visitor who picks a
Visium section therefore defeats the voxel ceiling without trying, and no
value of that ceiling changes the outcome.

This is not an argument, it is the shape of the call: the allocation happens
inside `deconvolve_sample`, before a `Domain` exists.

---

## 3. What a Visium section ACTUALLY costs — the finding that changes the answer

Peak of one cold `deconvolve_sample(sample_id)` per process, **before** any
change (3 cold runs shown where taken):

| sample | spots | peak over baseline (MiB) |
|---|---|---|
| `GBM3` | 1,775 | 3415.6 / 3449.8 / 3264.0 |
| `GBM1` | 4,764 | 3329.4 / 2866.3 / 3206.3 |

And `reconstruct_from_visium(sample, dx_um=30)`, one cold run each:

| section | voxels (`nx*ny`) | peak over baseline (MiB) |
|---|---|---|
| `visium_gbm3` | 18,300 | 3879.7 |
| `visium_gbm1` | 34,571 | 4131.4 |
| `visium_gbm2` | 32,760 | 4034.2 |
| `visium_10x_parent_gbm` | 30,996 | 4258.2 |

The staged breakdown for `visium_gbm3` shows where it goes:

```
imports                        109.0 MiB resident
load_sample(visium)            511.8        (+402.8  -- also grid-independent)
load_reference(scref)         1672.1        (peak 2314.5)
reconstruct_from_visium        ...          (peak 3154.2)
```

**`load_reference` is under half of it.** `load_sample` alone costs +400 to
+620 MiB, and `build_signatures` densifies a ~19,000 x 3,589 submatrix to
float64 **twice** (`ref_found`, then `ref_norm`), which is ~1.1 GB on its own.

> **Therefore: no fix to `load_reference` can make a Visium section fit in
> 512 MB.** The floor without it is still ~2.5 GB. This is the single most
> important result in this file and it is the reason the work below ends in a
> refusal rather than only an optimisation.

*Discrepancy noted, not hidden:* `results/performance_profile_2026-09-20.md`
§5.5 records 2,216–3,004 MiB peak for these sections and voxel counts of
1,665–3,190 at the same `dx = 30`. My voxel counts are `nx*ny` from
`reconstruct_from_visium`; that table's are ~10x smaller, so the two are
counting different things (almost certainly tissue-masked vs. full domain),
and the memory figures likewise measure different call boundaries. **I did not
reproduce that table and make no claim about it**; the numbers above stand on
their own method, which is stated in §0.

---

## 4. What was changed

### 4.1 `load_reference()` — bounded parse, int32, selections applied during the fill

The old body parsed the whole matrix to int64 and then made **two more
whole-matrix copies** (one to drop 5 HTSeq summary rows, one to reorder
columns). Measured: the bare `read_csv` is +1,444 to +1,506 MiB, and the
function around it reached +1,586 to +1,961 MiB.

Three value-preserving changes: **int32** (max count 77,122, and the only
consumers convert to float64; integer reductions still accumulate in int64 —
verified: `expr.sum().dtype` is `int64`, grand total 1,735,935,786);
**chunked** into a preallocated array (`_PARSE_CHUNK_ROWS = 2000`); and the
row/column **selections applied during the fill** so the two extra copies
never happen.

Chunk size was measured, not guessed — parse in isolation:

| chunk rows | 500 | 2,000 | 8,000 |
|---|---|---|---|
| peak (MiB) | 479.5 | **451.8** | 585.9 |
| wall (s) | 5.51 | 3.42 | 3.11 |

**After**, five cold runs of the whole function: **609.8 / 576.1 / 593.2 /
592.5 / 590.4 MiB**, resident frame **322.6 MiB**.

| | before | after | ratio |
|---|---|---|---|
| peak (max of 5 cold runs) | 1,960.6 MiB | **610.0 MiB** | **3.2x less** |
| peak (min of 5) | 1,586.3 MiB | 576.1 MiB | 2.8x less |
| resident frame | 643.8 MiB | **322.6 MiB** | 2.0x less |
| wall | 2.07 s | **3.85 s** | **1.9x SLOWER** |

**Bit-exactness is verified, not asserted.** Against the original
`pd.read_csv(path, sep=r"\s+", index_col=0)` on the real file: gene index
`array_equal` True, columns identical, `np.array_equal(new.astype(int64),
old)` **True**. End to end, `deconvolve_sample` output for `GBM3` and
`10x_parent_gbm` is **bit-identical** (`coords` and `fractions`, max|Δ| = 0.0,
same `type_names`) to the pre-change code.

### 4.2 `reference_gene_symbols()` — the load that never needed to happen

`deconvolve_sample` parsed the entire 23,460 x 3,589 matrix to build
`set(expr_df.index)` — a set of strings. The new accessor reads only the
symbol column:

| | full parse | symbol column |
|---|---|---|
| peak | +1,444 MiB | **+30.0 MiB** |
| wall | 1.94 s | 0.72 s |

Verified `array_equal` to `load_reference()[0].index` **and** to the returned
`gene_names`, order included.

### 4.3 The cache is on the SMALL object, deliberately

`results/performance_profile_2026-09-20.md` §6 #2 (D2-c) already called an
`lru_cache` on `load_reference()` "the wrong fix on the machine that matters".
The measurement here agrees and sharpens it: caching the frame pins **322.6
MiB** resident (not the 643.8 MiB it would have been, but still ~63 % of a
512 MB container) for the life of the process.

**Mutation audit, because a cache is only safe if nothing writes:** the only
callers are `deconvolve.build_signatures` (`.loc`/`[]` — copies),
`deconvolve_sample`, `validation/reconstruction_fidelity` (`expr_df[train_cells]`
— a copy) and `viz/explainer` (shapes only). **Nothing mutates the returned
objects.** Caching the frame would nonetheless be wrong for the memory reason
above, so it is not cached, and `load_reference`'s docstring says why.

What IS cached is `build_signatures`' result: a `(G_shared, 4)` float64 matrix
of a few hundred kB, keyed on `(reference_fingerprint(), tuple(shared_genes))`.

* The fingerprint is `(path, size, mtime_ns)` of both files, so a replaced
  reference invalidates rather than serves stale.
* `shared_genes` **must** be in the key: `_normalize_columns_to_unit_sum`
  divides each cell by its total over the shared genes, so the signature for
  gene set A is not a row-slice of the signature for gene set B. Per-section
  is the largest correct unit.
* The result is **copied on read**, so a caller mutating it cannot poison
  another's.
* Bounded at 8 entries.

Measured effect — second `deconvolve_sample` of the same section in one
process:

| section | 1st run | 2nd run | speedup | extra peak |
|---|---|---|---|---|
| `GBM1` | 6.96 s | 2.47 s | 2.8x | **none** |
| `GBM2` | 6.29 s | 2.00 s | 3.1x | **none** |
| `GBM3` | 5.88 s | 1.60 s | 3.7x | **none** |
| `10x_parent_gbm` | 6.99 s | 2.18 s | 3.2x | **none** |

"None" is literal: `ru_maxrss` after run 2 equals `ru_maxrss` after run 1 to
the decimal, so the second run allocated nothing beyond run 1's high-water.

### 4.4 Two refusals, both quoting measured numbers

Both read the limit from `app/server/account.py::memory_limit()` (cgroup v2 →
v1 → `TMESIM_HOST_MEMORY_MB`) via a lazy, guarded import — there is **no
second copy of the detection logic**. Both refuse only against a **declared**
limit; undeclared memory raises nothing. Both name
`TMESIM_SCREF_ALLOW_OVERSIZE=1` as the escape hatch.

`load_reference` → `ReferenceTooLargeForHost`, threshold 610 MiB (the max of
the five measured runs, not the mean — a refusal threshold set at the average
is wrong half the time in the direction that kills the container).

`deconvolve_sample` → `DeconvolutionTooLargeForHost`, thresholds from the
measured per-section table (max of ≥3 cold runs each):

| sample | measured peak used as the threshold (MiB) |
|---|---|
| `GBM1` | 3454.8 |
| `GBM2` | 3151.8 |
| `GBM3` | 3156.1 |
| `10x_parent_gbm` | 3146.1 |
| *anything else* | 2924.0 — the cheapest Visium deconvolution ever measured here |

This check runs **before `load_sample`**, which alone is +400 to +620 MiB; a
test booby-traps `load_sample` so that reaching it fails rather than merely
costing memory.

Verified live at `TMESIM_HOST_MEMORY_MB=512`: refuses. At 8192: runs. With the
override: runs. With no declared limit: runs.

---

## 5. Net effect on the live 512 MB instance

| | before | after |
|---|---|---|
| Visium section requested | loads ~3–4 GB, **container OOM-killed** | **refused in milliseconds**, with the measured requirement, the detected limit and its source |
| Grid ceiling's role | none — it counts voxels | unchanged, and the refusal says so explicitly |
| Xenium sections | unaffected | unaffected (they never load this reference) |

The four Visium sections are **not fixed** on 512 MB and cannot be by this
work. They now fail loudly and legibly instead of taking the service down for
everyone.

### 4.5 ONE INTEGRATION STEP IS NOT DONE, and it needs a file this pass did not own

The refusal fires and the container is saved — verified through the server's
own entry point, `reconstruct_real(source='visium', sample='GBM3', dx_um=30)`
at `TMESIM_HOST_MEMORY_MB=512`: it raises in milliseconds with the process at
**124.7 MiB** instead of climbing to ~3.2 GB.

But the visitor does not see the refusal's words. `app/server/jobs.py` renders
a failed run's `job.error` as `str(exc)` only when `is_run_refusal(exc)` is
true, which matches class names against `RUN_REFUSALS`:

```
DeconvolutionTooLargeForHost     is_run_refusal=False
ReferenceTooLargeForHost         is_run_refusal=False
RUN_REFUSALS = ['CompoundError', 'CompoundRefused', 'IngestError',
                'ProjectError', 'TplError']
```

So today a visitor on the 512 MB instance gets the **redacted sentence plus a
reference id**, and the measured explanation goes to stderr and
`GET /api/runs/{id}/traceback`. The operator can still read it; the user
cannot.

Closing that is two small edits in files outside this pass's ownership, and
they are **reported rather than made**:

1. `app/server/jobs.py` — add `"ReferenceTooLargeForHost"` and
   `"DeconvolutionTooLargeForHost"` to `RUN_REFUSALS`. Both messages are
   authored strings over integers and an env-var name; neither interpolates a
   filesystem path or a third-party traceback, which is the property that list
   exists to protect.
2. `app/server/api.py` — mirror the check in the submit-time validator, the
   way the comment above `RUN_REFUSALS` says every deliberate refusal is
   mirrored, so a Visium request on a small host is a 400 at `POST /api/runs`
   rather than a job that starts and fails. `deconvolve.check_visium_memory_budget(sample_id)`
   is the whole call.

Neither is required to stop the OOM. Both are required for the refusal to be
readable by the person it is addressed to.

*(Unrelated, found while running the suite and reported for whoever owns it:
`tests/test_response_redaction.py::test_no_undeclared_exception_text_on_the_http_surface`
fails on `app/server/cache.py:295 state_dir_problem`. That is not from this
work — the finding names only that one site, in a file this pass did not
touch.)*

---

## 6. Claims elsewhere that these measurements bear on

* **`results/performance_profile_2026-09-20.md` §6 #2 — CONFIRMED.** Its
  1,612.5 MiB standalone figure reproduces, its "not cached" and "2.5x parse
  transient" claims hold, and its judgement that `lru_cache` on
  `load_reference` is the wrong fix is correct for the reason it gives.
* **Same file, §5.5 — NOT REPRODUCED, and not contradicted either.** Its
  per-section peaks and voxel counts are measured at a different call boundary
  from mine (§3). Flagged so nobody averages the two tables together.
* **Same file, §6 #2's implied remedy.** It frames the reference as "the
  largest single allocation in this product" and reads as though fixing it
  addresses the Visium OOM. It does not: §3 here shows the reference is under
  half the cost and the remainder still exceeds 512 MB by ~5x. That file is a
  dated published record and is **not edited**; this is the correction.
* **`docs/real_data/scref_gbm.md`** states the full atlas "is cached and
  loaded — no subsample". Still true. It says nothing about dtype; the frame
  is now int32 rather than int64, with values unchanged.

---

## 7. Not fixed, and why

* **`build_signatures`' two float64 densifications (~1.1 GB).** This is now
  the dominant term. Both are avoidable in principle (the matrix is 89 %
  zeros; a CSR path would be ~106 MiB) but neither was attempted here: the
  change is inside the deconvolution arithmetic, a sparse mean is not
  trivially bit-identical to a dense one in floating point, and shipping a
  half-verified numerics change to a live scientific product is worse than
  shipping a refusal. **This is the next piece of work and it is the one that
  could actually make Visium fit.**
* **`load_sample`'s +400 to +620 MiB.** Outside this brief's file ownership
  (`tmesim/data/real/visium.py`), and it does not change the 512 MB verdict on
  its own.
* **`validation/reconstruction_fidelity.py`** still calls `load_reference()`
  to build a symbol set, the way `deconvolve_sample` used to. It is a
  validation module, not on the live request path, so it was left alone rather
  than changed without a measurement of its own.
* **The 1.9x slowdown of `load_reference`** (2.07 s → 3.85 s) was accepted,
  not eliminated. The signature cache more than repays it from the second run
  onward, but the first run of a cold process is genuinely slower.
* **Linux figures.** Everything here is darwin. The array sizes are
  platform-independent; the allocator behaviour (the ratchet in §1, the spread
  in §4.1) is not. The refusal thresholds are consequently conservative on a
  platform whose allocator returns memory more readily.
