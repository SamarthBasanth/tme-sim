# The three copies in the Visium loader, and the measured floor that says 512 MB is not coming

**Date:** 2026-10-10
**Scope:** `tmesim/data/real/visium.py::load_from_dirs` (and the two helpers it
grew), plus the refusal thresholds in
`tmesim/reconstruct/deconvolve.py::check_visium_memory_budget`.
**Status of the headline claim:** the reported defect is REAL, the fix is DONE
and it is **bit-identical** — `deconvolve_sample`'s `coords` and `fractions`
are `np.array_equal` to the pre-change values on all four sections, max|Δ|
exactly 0.0 — and `load_sample`'s MEASURED cold peak falls **1.95x to 2.15x**.
It does **not** make a Visium section runnable on a 512 MB instance, the
refusal stays, and this file closes that question with a measurement instead
of an estimate: **512 MB is not reachable, for any of the four sections, even
by a rewrite that does not exist yet.**

This is the piece of work that `results/scref_signatures_2026-09-20.md` §7
named and declined:

> **`load_sample`'s +633.7 MiB, now the largest single term.** In
> `tmesim/data/real/visium.py`, outside this pass's ownership. [...] with the
> caveat that removing all of it still leaves ~600 MiB, which is still over
> 512 MB. **Nothing measured here suggests a Visium section can be made to fit
> in 512 MB; the honest target for that file is a smaller overrun, not a fix.**

Four things that file did not know are established below: the +633.7 MiB was
**three full-matrix copies and a mis-attributed transpose**, not a parse; **two
of the three were free to remove** (one selected every column, the other needed
no data movement at all); the term that the Visium loader was hiding is
`build_signatures`, whose staged cost **went UP** when the Visium loader got
cheaper; and its "~600 MiB" floor is right in the middle of what a streamed
pipeline actually measures — **531.5 to 707.2 MiB**, all four of them over 512.

---

## 0. Method, and what the numbers are numbers OF

Peak memory is `resource.getrusage(RUSAGE_SELF).ru_maxrss`.

> **`ru_maxrss` units are platform-dependent.** On **darwin/macOS it is
> BYTES**; on Linux it is KiB. Every figure here was taken on darwin and
> divided by 1024² **once**, not twice. Current RSS is `ps -o rss= -p <pid>`,
> which on darwin reports **KiB**, so that one is multiplied by 1024 first.
> `psutil` is not installed in this venv.

"Peak" means `ru_maxrss_after − ru_maxrss_before` around the call being
measured. `ru_maxrss` is a process high-water mark and never decreases, so this
is a *rise*, and ordering within a process matters — every headline figure is
therefore a **cold subprocess**, one measurement per process.

Host: darwin 24.6.0, 16 GiB, CPython 3.9, pandas 2.3.3, numpy 2.0.2, scipy
1.13.1 (whose `mmread` is the `fast_matrix_market` backend). Scripts are in the
session scratchpad.

**How before and after are both measured in one session.** The venv carries an
*editable* install of `tmesim` whose finder points at the main checkout
`/Users/samarthbasanth/tme-sim`. A `python script.py` run does NOT put the cwd
on `sys.path`, so a script imports that checkout; inserting the worktree root
on `sys.path` imports this branch instead. The main checkout is byte-identical
to this worktree's HEAD for all three owned files (`diff`-verified:
`visium.py`, `deconvolve.py`, `scref_gbm.py`), and the worktree's
`data_cache/visium_gbm` is a symlink to the main one, so BEFORE and AFTER are
the same interpreter, the same libraries and the same input bytes with one file
different.

> **A trap worth recording, because it produced a clean and completely false
> "no change" result.** The first round of "after" measurements reported
> 641.4 MiB against a 642.0 MiB before — and the first round of bit-exactness
> fingerprints came back "BIT-IDENTICAL" on all four sections. Both were
> measuring the *unmodified* main checkout, because the harness scripts were
> run as scripts and never had the worktree on `sys.path`. The giveaway was an
> `AttributeError` on a helper that demonstrably existed in the file. A
> before/after harness that silently measures "before" twice reports a perfect
> null and a perfect bit-match, which is exactly what a real no-op looks like.

*The run-to-run spread is large and is not hidden.* It is the allocator not
returning freed arenas, the same effect both previous files documented. Every
threshold below is the **maximum** observed and every ratio is computed
**within** this harness; ratios across harnesses are not meaningful and are not
taken.

---

## 1. The reported figure, re-measured

`results/scref_signatures_2026-09-20.md` §4.1 reports **+633.7 MiB** for
`load_sample` on GBM1. Three cold runs here, same staging (peak rise of
`load_sample` over the post-import baseline):

| | run 1 | run 2 | run 3 |
|---|---|---|---|
| GBM1 | 644.2 | 634.3 | 641.5 |

**The reported figure reproduces** — 633.7 sits inside the spread and matches
run 2 to 0.1 %. All four sections, before, maximum of three cold runs each:

| section | spots | nnz | peak (MiB) | the CSR it returns |
|---|---|---|---|---|
| `GBM1` | 4,764 | 15,232,681 | **644.2** | 174.3 MiB |
| `GBM2` | 2,891 | 9,208,983 | 411.7 | 105.4 MiB |
| `GBM3` | 1,775 | 9,089,279 | 410.8 | 104.0 MiB |
| `10x_parent_gbm` | 3,468 | 14,884,287 | 625.4 | 170.4 MiB |

Every section has 36,601 gene columns. The cost tracks `nnz`, not spot count —
GBM3 has 1,775 spots and costs the same as GBM2's 2,891, because they hold
nearly the same number of stored values.

### 1.1 Where it went, statement by statement

Staged inside `load_from_dirs` on GBM1, `ru_maxrss` marked between statements:

```
imports                           96.2 MiB
_read_positions/_features/_barcodes +5.3
_read_matrix (mmread + .tocsr)   +455.5
mat[:, keep_idx]                   +0.3     <-- a whole-matrix copy, invisible
  .T                               +0.0         (it landed in arenas mmread
  .tocsr()                       +174.5         had just freed)
expr.sum(axis=1)                   +0.0
                                 = 630.4
```

Split again, `mmread` from `.tocsr()`:

| | peak rise | object |
|---|---|---|
| `mmread(gzip handle)` → COO | **+286.4 MiB** | 232.4 MiB (`int64` data, `int32` row, `int32` col) |
| `.tocsr()` | +116.0 | 174.5 MiB |

**`mmread` itself is near-optimal**: 62 MiB of internal buffers over the 232.4
MiB it returns. Everything else was the three conversions that followed it.

Note the second line of the staged table: `mat[:, keep_idx]` is a full CSR
column gather and it shows up as **+0.3 MiB**, because `ru_maxrss` is a
high-water mark and the copy fitted inside arenas the COO had just released. A
copy that costs nothing on this measurement still costs on a container that is
already near its limit, and it is the first of the two free removals.

---

## 2. Two facts about the data, both measured, both free

### 2.1 The column gather selects every column

`barcodes.tsv` in a **filtered** Space Ranger matrix already lists only
in-tissue spots, so `keep_idx` — "the barcodes present in both the matrix and
the in-tissue position list" — is `arange(n_spots)`:

| section | barcodes in matrix | in-tissue positions | kept | is the identity |
|---|---|---|---|---|
| `10x_parent_gbm` | 3,468 | 3,468 | 3,468 | **yes** |
| `GBM1` | 4,764 | 4,764 | 4,764 | **yes** |
| `GBM2` | 2,891 | 2,891 | 2,891 | **yes** |
| `GBM3` | 1,775 | 1,775 | 1,775 | **yes** |

So `mat[:, keep_idx]` copied the entire matrix to select 100 % of its columns.
This is the same class of finding as the 321 MiB pair in
`results/scref_signatures_2026-09-20.md` §1.1, and it was worth looking for
again for the same reason: removing a copy that selects everything cannot
change a value, because it changes no value.

It is **not** deleted, only skipped when it is the identity. A *raw* Space
Ranger matrix carries ~1.5x more barcodes than there are in-tissue spots, and
that selection has to keep working for the first person who ingests one — see
§5.4 for the test that covers the branch no cached section exercises.

### 2.2 The transpose needs no data movement

A Matrix Market file from Space Ranger lists its entries **column-major**
(spot-major), with row indices ascending inside each column. In the
**transposed** orientation — which is the `(n_spots, n_genes)` the pipeline
wants — that is exactly CSR canonical order. MEASURED on GBM1's 15,232,681
entries:

| | |
|---|---|
| `coo.col` non-decreasing | **True** |
| `coo.row` strictly increasing within each column | **True** |

Therefore the gene index array **is** the CSR's `indices`, the value array
**is** its `data`, and `indptr` is a `searchsorted` over the spot index array.
Nothing needs to be copied, sorted or summed.

That is a property of these files and not a guarantee of the format, so the
code **checks** it on every load and falls through to the ordinary conversion
when it fails. The check costs **20 to 39 ms** on 9.1 M to 15.2 M entries.

---

## 3. What was changed

`tmesim/data/real/visium.py` only. Three new private helpers, one deleted:

* **`_read_matrix_coo`** replaces `_read_matrix`. It returns `mmread`'s COO
  *unconverted* (and converts a dense-`array` Matrix Market file to COO, which
  the old reader would have raised `AttributeError` on). The deleted function
  had no other caller.
* **`_is_csr_canonical_order(rows, cols, chunk=2**20)`** — rows
  non-decreasing AND cols strictly increasing within a row, compared in chunks
  so the check never allocates anything proportional to `nnz`. `np.diff(rows)`
  would cost 58 MiB on GBM1 and **`np.bincount` silently casts its input to
  `intp`**, which is 116 MiB for an `int32` index array; both were measured and
  both are why this is a loop rather than a one-liner.
* **`_spots_by_genes_csr(coo, spots_axis, keep_idx, n_genes)`** — the
  assembly. On the identity selection it adopts `coo.data` and the gene index
  array outright (`np.shares_memory` is **True** for both, and a test pins
  that); on a non-identity selection it renumbers the kept spots into the
  output block and filters the three flat arrays, which costs the size of the
  output rather than a sparse gather.
* The transpose guard became an **axis number** (`spots_axis`) instead of a
  physical `mat.T.tocsr()`. Nothing between the guard and the assembly reads
  the matrix, so there was nothing for the transpose to serve.
* `del mat` after the assembly. On the identity path the CSR holds two of the
  COO's three arrays, so this frees exactly the third — the spot index,
  58.1 MiB on GBM1.

`indptr` comes from
`np.searchsorted(rows, np.arange(n_keep + 1, dtype=rows.dtype), side="left")`.
The needle shares `rows`' dtype deliberately: a mismatched needle makes
`searchsorted` promote — and therefore copy — the big array.

### 3.1 The one-row trap does not apply here, and why that is worth saying

`results/scref_signatures_2026-09-20.md` §3.3 found that `block[:, idx]` picks
its memory order from the SHAPE, so a chunk of exactly one row reduces
pairwise where every other chunk size reduces sequentially — 45.7 % of
signature entries and 69.9 % of one section's fractions moved, inside any
tolerance anyone would have set. That trap is about a **floating-point
reduction**. Nothing in this path is one: the values here are integers that get
**moved**, never summed, and the only arithmetic is `searchsorted` on indices.
What can move here is an index, and an index that moves silently relabels a
spot or a gene — so every comparison below is `np.array_equal` on the raw
`data` / `indices` / `indptr` arrays, and §5.4's plant shows exactly that
failure mode when the renumbering is wrong.

---

## 4. After

`load_sample`, three cold runs each, interleaved with the before runs:

| section | before (max of 3) | after (max of 3) | ratio | after, all three runs |
|---|---|---|---|---|
| `GBM1` | 644.2 | **299.3** | **2.15x** | 283.1 / 292.5 / 299.3 |
| `GBM2` | 411.7 | **211.4** | 1.95x | 201.0 / 211.4 / 200.7 |
| `GBM3` | 410.8 | **201.8** | 2.04x | 200.8 / 200.1 / 201.8 |
| `10x_parent_gbm` | 625.4 | **292.2** | 2.14x | 289.7 / 292.2 / 292.0 |

Staged inside the new `load_from_dirs`, GBM1:

```
imports                           96.8 MiB
tsv/csv readers                    +6.6
_read_matrix_coo (mmread)        +294.5     <-- all of it, and it is mmread's
keep_idx + barcodes                +0.0
_spots_by_genes_csr                +1.8         shares data: True
del coo                            +0.0         shares indices: True
expr.sum(axis=1)                   +0.1
                                 = 296.4
```

The remaining peak **is `mmread`**. The residual object is 174.3 MiB and the
peak is 296.4, and the 122 MiB between them is the COO's spot index array
(58.1 MiB, discarded immediately) plus `mmread`'s own 62 MiB of buffers. §8
says what it would take to get at that and why it was not taken.

**Wall time improved**, which was not the goal and is recorded because the two
previous passes both had to accept a slowdown:

| | before | after |
|---|---|---|
| `load_sample("GBM1")` | 1.012 / 0.924 s | **0.661 / 0.798 s** |
| `load_sample("GBM3")` | 0.611 / 0.575 s | **0.479 / 0.455 s** |

---

## 5. Bit-exactness: what was checked, and what it cost to be sure

### 5.1 `load_sample`'s whole return value, byte for byte

Every array in the returned dict was fingerprinted with sha256 **of its raw
bytes**, alongside its dtype, the CSR's `has_sorted_indices` /
`has_canonical_format` flags, the gene / barcode / ensembl lists, `umpp`,
`spot_diameter_fullres`, the marker report and the full key set. Before vs
after, cold process each:

| section | `diff` of the fingerprint |
|---|---|
| `GBM1` | **identical** |
| `GBM2` | **identical** |
| `GBM3` | **identical** |
| `10x_parent_gbm` | **identical** |

Hashing the bytes rather than comparing values is deliberate: it makes a dtype
change — `int64` → `int32` data would be the obvious next saving — fail the
comparison instead of passing it.

### 5.2 `deconvolve_sample`, end to end

| section | spots | `coords` | `fractions` | `type_names` | max|Δ| |
|---|---|---|---|---|---|
| `GBM1` | 4,764 | equal | **equal** | equal | 0.0 |
| `GBM2` | 2,891 | equal | **equal** | equal | 0.0 |
| `GBM3` | 1,775 | equal | **equal** | equal | 0.0 |
| `10x_parent_gbm` | 3,468 | equal | **equal** | equal | 0.0 |

Not `allclose`. `np.array_equal`, max absolute deviation exactly 0.0.

### 5.3 That comparison is not vacuous, and here is the proof

A plausible "defensive" defect — `np.sort(gene_of)` on the gene index before
handing it to the constructor, on the grounds that a sorted index array is what
CSR wants — produces a matrix scipy accepts without complaint (a globally
sorted index is per-row sorted too) with every value attached to the wrong
gene. Run end to end on GBM3 against the saved pre-change output:

| | |
|---|---|
| `fractions` equal | **False** |
| entries differing | **6,634 / 7,100 (93.4 %)** |
| max absolute deviation | 1.0 |
| `coords` equal | True |

So the end-to-end comparison in §5.2 can see a wrong answer, and the
`coords`/`fractions` split localises it.

### 5.4 The branch no cached section exercises

Every cached section selects all of its spots in order (§2.1), so the
non-identity selection is dead code on this machine's data and live on anyone's
**raw** Space Ranger matrix. It is covered by a toy-matrix test against the
gather it replaced, over four selections (a subset, a different subset, an
out-of-order pair, a single spot), and the plant that proves the test is the
most instructive one in the file: keeping the FILE's spot numbering instead of
renumbering into the kept block gives

```
got  [[1 0 2 0]        want [[1 0 2 0]
      [0 0 0 0]              [4 0 0 5]
      [4 0 0 5]]             [0 0 7 8]]
```

— spot 4's entries off the end of the block, spot 2's one row late. A **silent
relabelling of which spot is which**, which is the worst shape this defect
could take and the reason the comparisons are on indices and not on sums.

---

## 6. Memory, after, at the pipeline level — and the term that went UP

Cold `deconvolve_sample`, one per process. Three runs before; eight to twelve
after. Maximum and minimum of each, same harness throughout:

| section | before | after | ratio (max/max) | after runs |
|---|---|---|---|---|
| `GBM1` | 1234.4 (min 1182.4) | **895.6** (min 835.0) | 1.38x | 8 |
| `GBM2` | 896.3 (min 680.2) | **782.3** (min 754.0) | 1.15x | 8 |
| `GBM3` | 905.1 (min 893.1) | **787.4** (min 565.8) | 1.15x | 12 |
| `10x_parent_gbm` | 1195.6 (min 1184.5) | **896.7** (min 678.9) | 1.33x | 12 |

A 2.0x–2.15x saving in `load_sample` buys only 1.15x–1.38x at the pipeline
level, and the reason is the most counter-intuitive measurement in this file.
Staged, cold, GBM1, before and after:

```
BEFORE                                   AFTER
imports              109.5 MiB           imports              109.3 MiB
load_sample          746.7  (+637.1)     load_sample          398.3  (+289.0)
reference_gene_syms  770.2  ( +23.5)     reference_gene_syms  425.3  ( +27.0)
build_signatures     864.5  ( +94.3)     build_signatures     952.8  (+527.5)
CSR subset + NNLS    864.5  (  +0.0)     CSR subset + NNLS    952.8  (  +0.0)
                   = 754.8                                  = 843.5
```

**`build_signatures`' staged contribution went from +94.3 MiB to +527.5 MiB
without a line of it changing.** It was never cheap; it was landing in the
~450 MiB of freed-but-retained allocator arenas the old Visium loader left
behind. Taking those away makes its real cost visible, and the real cost is
now the dominant term: `load_reference()` assembling the whole 23,460 x 3,589
`int32` reference frame (321 MiB resident) plus its chunked-parse transient.

The dominant term has therefore moved OUT of the Visium loader and into
`tmesim/data/real/scref_gbm.py::load_reference`. Which is the natural next
piece of work by the same argument that made this one — except that §7 says it
would not be enough.

---

## 7. Can 512 MB be reached? NO, and this is now measured

`results/scref_signatures_2026-09-20.md` §7 projected that removing all of
`load_sample` "still leaves ~600 MiB". That projection is replaced here by a
measurement of the actual floor.

**The floor run** executes the whole computation in the cheapest shape it has:

* `load_sample` as shipped (§3), then the shared-gene column subset built, then
  **the full spot x gene matrix dropped** — the smallest live set the pipeline
  can have at the moment the reference must be read;
* the reference **streamed** — two chunked passes straight off the gzip, column
  sums then per-group means, with the 321 MiB frame **never assembled**. This
  is not a proposed implementation and no code was written for it in `tmesim`;
  it is the arithmetic `_signatures_streamed` already does, with the one
  allocation it cannot currently avoid removed.

Nothing in the run holds an object the computation does not need. Absolute peak
RSS, one cold process each:

| section | shared-gene CSR | RSS at the tightest point | **absolute peak** | over 512 MiB by |
|---|---|---|---|---|
| `GBM3` | 92.4 MiB | 392.0 | **531.5** | 1.04x |
| `GBM2` | 94.3 MiB | 389.4 | **542.1** | 1.06x |
| `GBM1` | 155.4 MiB | 530.4 | **686.6** | 1.34x |
| `10x_parent_gbm` | 152.8 MiB | 522.1 | **707.2** | 1.38x |

**All four are over 512 MiB.** The cheapest section misses by 4 %, the dearest
by 38 %. So no rewrite of `load_reference` or `build_signatures` gets a Visium
section into a 512 MB container either, and the refusal is not a placeholder
waiting for one more optimisation — it is the answer.

**What the gap is made of, honestly.** At GBM1's tightest point the *live*
objects are ~278 MiB (155.4 MiB sub-matrix + 27 MiB gene symbols + ~95 MiB of
interpreter and libraries) against an RSS of 530.4. The 250 MiB difference is
allocator arenas that were freed and not returned, which `ru_maxrss` counts and
a cgroup counts. A sum of live objects would say ~380 MiB and fit; the process
says 530 and does not. On a platform whose allocator releases more readily
(glibc with `malloc_trim`, say) the floor would be lower than the numbers above
— they are darwin figures and §9 says so — but they are the numbers for the
machine these measurements were taken on, and nothing here was measured on the
deployment platform.

---

## 8. Not fixed, and why

* **`mmread`'s own 294.5 MiB peak**, now the whole of `load_sample`. 174.3 MiB
  of it is the answer; 58.1 MiB is the COO's spot index array, which exists
  only to be `searchsorted` into a 19 kB `indptr`; 62 MiB is `mmread`'s
  internal buffering. A hand-written streaming Matrix Market parser could
  write the final CSR arrays directly in file order — the file is already in
  that order (§2.2) — for an estimated ~190 to 200 MiB peak, roughly 100 MiB
  better. **Not attempted.** It does not change the §7 verdict (the floor is
  531 to 707 MiB and this would move it by ~100), and it means owning an
  integer parser whose every value must match `mmread`'s on a file nobody will
  re-check. Trading a library's parser for a local one to shave 15 % off a
  figure that is already over the limit is the wrong trade.
* **`int32` counts.** MEASURED maxima: GBM1 950, 10x_parent_gbm 1,118, GBM3
  1,435, GBM2 2,372 — six orders of magnitude of headroom under `int32`, and
  `scipy.sparse`'s row sum upcasts to `int64` regardless of the stored width,
  so `total_counts` would not narrow. It would save 34.7 MiB (GBM3) to 58.1 MiB
  (GBM1) of resident. **Not taken**: the dtype is part of what ~30 call sites
  read out of `sample["expression"]`, bit-exactness in this pass was defined to
  include it (§5.1 hashes bytes precisely so a dtype change fails), and it buys
  less than the §7 gap.
* **`load_reference()`'s 321 MiB frame / ~527 MiB transient**, now the dominant
  term. A streamed `build_signatures` is measured in §7 as the floor run and
  it does not reach 512 MB, so it is an optimisation and not a fix; it was not
  attempted here. If it is attempted, `results/scref_signatures_2026-09-20.md`
  §3.3's one-row trap applies to it directly and its §2.1 sweep (chunk sizes
  250 through 8,000 and the dense path, all bit-identical) is the evidence that
  the chunk boundary is free to move.
* **`app/server/jobs.py`'s `RUN_REFUSALS`** still does not name
  `DeconvolutionTooLargeForHost` or `ReferenceTooLargeForHost`, so a visitor on
  a small host still gets the redacted sentence while the measured explanation
  goes to stderr and `GET /api/runs/{id}/traceback`. Open since
  `results/scref_memory_2026-09-20.md` §4.5, outside this pass's ownership,
  untouched. The thresholds that refusal quotes are now lower and its wording
  is now correct; neither makes it more visible.
* **`tmesim/validation/reconstruction_fidelity.py`** still calls
  `load_reference()` to build a symbol set. Not on the live request path, not
  measured here, flagged by both previous passes and still true.
* **Linux figures.** Everything here is darwin. The array sizes are
  platform-independent; the allocator retention that §7's floor is half made of
  is not.

---

## 9. The refusal thresholds, re-measured, and verified in BOTH directions

The old table was measured against the old loader and is now 15 % to 37 % too
high. A threshold above the real cost refuses hosts that can serve the section,
so it had to move:

| section | 2026-09-20 | **2026-10-10** | runs behind the new figure |
|---|---|---|---|
| `GBM1` | 1227.9 | **895.6** | 8 |
| `GBM2` | 939.0 | **782.3** | 8 |
| `GBM3` | 922.2 | **787.4** | 12 |
| `10x_parent_gbm` | 1215.7 | **896.7** | 12 |
| `MEASURED_VISIUM_MIN_MIB` | 869.9 | **565.8** | — |

Each section's figure is the **maximum** observed, because a refusal threshold
set at the average is wrong half the time in the direction that kills the
container. `MEASURED_VISIUM_MIN_MIB` is the one constant where the **minimum**
is right — it is the threshold for a section with no measurement of its own,
and its documented meaning is "the figure below which nothing measured here has
ever fitted", so it must never refuse a configuration that was observed to
work. 565.8 is a single low GBM3 run against a mode of 750 to 790 for that same
section; it is kept because the definition demands it, and it is recorded here
as the outlier it is.

The refusal message's explanation was also corrected. It said the cost was
"dominated by loading the Visium section itself — its spot x gene matrix and
image metadata"; that was true when it was written and §6 says it is not true
now. It names the reference panel instead.

**Verified live, through `deconvolve_sample` and through the server's own entry
point `reconstruct_real(source='visium', ...)`:**

| host limit | section | outcome | process at |
|---|---|---|---|
| 512 MiB | `GBM1` | **REFUSED** | 124.6 MiB |
| 512 MiB | `GBM2` | **REFUSED** | 125.4 MiB |
| 512 MiB | `GBM3` | **REFUSED** | 126.0 MiB |
| 512 MiB | `10x_parent_gbm` | **REFUSED** | 124.9 MiB |
| 512 MiB | `some_ingested_section` (no measurement) | **REFUSED** | 125.9 MiB |
| 512 MiB | `GBM3` via `reconstruct_real` | **REFUSED** | 126.7 MiB |
| 8192 MiB | `GBM3` | **ran** | peak 846.4 MiB |

and the half that matters for a lowered threshold — a host the OLD table
refused and the NEW one serves, with the actual run measured to check that the
new threshold is not merely smaller but **right**:

| host limit | section | thresholds | outcome | actual peak |
|---|---|---|---|---|
| 1100 MiB | `GBM1` | 2026-09-20 (1227.9) | **REFUSED** | — |
| 1100 MiB | `GBM1` | 2026-10-10 (895.6) | **ran** | 980.6 MiB absolute, 871.0 rise |

That is 332 MiB of refusal that was pointed at hosts which can now serve a
section, and the run it was blocking fits inside the limit with 119 MiB to
spare.

*(Both rows are run in the AFTER tree with only the table swapped. The
`before` tree cannot be used for this comparison: the limit is detected via a
lazy `from app.server.account import memory_limit`, and `app` is only
importable when the worktree root is on `sys.path`, so a "before" script
reports "not declared" and never refuses anything. A refusal test that silently
cannot detect a limit passes by never raising, which is the same shape of
mistake as §0's.)*

---

## 10. Plant evidence

Every test in `tests/test_visium_matrix_assembly.py` was confirmed RED by
planting the specific defect it targets, confirming the plant **parses**
(`py_compile`) and is reached, running the test, then reverting and confirming
the file is byte-identical to the pre-plant copy. **One of the nine tests was
vacuous on the first attempt and the plant is what caught it.**

| plant | test | observed |
|---|---|---|
| `side="left"` → `side="right"` in the `indptr` `searchsorted` | assembled == converted, real section | **RED** `ValueError: index pointer should start with 0` |
| `_is_csr_canonical_order` → `return False` | the real sections take the assembled path | **RED** `AssertionError: 10x_parent_gbm is not in canonical order` |
| `data = coo.data` → `coo.data.copy()` | the assembly copies neither big array | **RED** `AssertionError: data was copied` |
| the within-row column comparison deleted | only sorted-unique pairs are canonical | **RED** `AssertionError: a column that goes backwards inside a row is not canonical` |
| windowed scan that never compares across a seam | a violation that straddles a chunk | **RED** `assert not True ... chunk=2` (only after the test was fixed — see below) |
| `if _is_csr_canonical_order(...)` → `if True` | a non-canonical matrix falls back and sums duplicates | **RED** `[[5, 11, 9], [0, 0, 0]]` against `[[5, 7, 0], [0, 4, 9]]` |
| `pos[keep_idx] = arange(n_keep)` → `= keep_idx` | a non-identity selection matches the gather | **RED** spot 4 dropped, spot 2 one row late (§5.4) |
| `spot_of = coo.col if spots_axis == 1 else coo.row` → `= coo.col` | a transposed file is read through the other axis | **RED** `ValueError: axis 0 index 2 exceeds matrix dimension 2` |
| `shape=(n_keep, n_genes)` → `(n_genes, n_keep)` | `load_sample` returns a canonical int CSR | **RED** `ValueError: index pointer size (1776) should be (36602)` |
| `np.sort(gene_of)` as the CSR `indices` | the end-to-end output comparison (§5.3) | **RED** 93.4 % of GBM3's fractions moved |

### 10.1 The vacuous check this pass wrote, and then caught

Recorded rather than quietly fixed, because the point of the exercise is that
they are easy to write.

**`test_canonical_order_check_sees_a_violation_that_straddles_a_chunk`** used
`rows = [0, 1, 2, 1]` with `chunk=2`. The violating pair is (index 2, index 3)
— which falls **inside** a window under both the correct indexing and the
planted one, so the seam plant left the test **green** while the whole file
stayed green around it. A seam test has to put the violation on the seam: the
fixture is now `rows = [0, 2, 1, 3]`, whose violation is at (index 1, index 2),
and the repeated-pair case was moved the same way. This is the third time in
three passes that a windowed or chunked check needed its fixture aimed at the
boundary rather than at the middle.

Two of the nine plants were also informative in the other direction: scipy's
own `check_format` rejects a shifted `indptr` and a transposed shape before any
comparison runs. That is a genuine guard and it is not a substitute for the
comparison — it catches damage that changes the array *sizes*, and says nothing
about damage that changes which gene a count belongs to, which is what §5.3's
plant does.

### 10.2 Exit codes

```
tests/test_visium_matrix_assembly.py        9 passed                      (exit 0)
tests/test_signature_streaming.py          15 passed                      (exit 0)
tests/test_scref_memory.py                 15 passed                      (exit 0)
tests/test_ingest.py + test_gbmspace_adapter.py + the three above:
                            3 failed, 85 passed, 41 skipped               (exit 1)
tests/  (whole suite, after)  84 failed, 2210 passed, 153 skipped, 173 errors (exit 1)
tests/  (whole suite, before) 84 failed, 2201 passed, 153 skipped, 173 errors (exit 1)
```

**The two whole-suite failure sets are identical.** Both runs' `FAILED` and
`ERROR` lines were collected and sorted: 257 lines each, `diff` empty. Every
one of the 84 failures and 173 errors is therefore pre-existing and unrelated
— they are the `app/` server fixtures, the Xenium exports this worktree's
`data_cache` does not carry, and the rest of the standing redness of this
branch. The 9-test difference in the pass count is exactly
`tests/test_visium_matrix_assembly.py`.

That comparison was done properly rather than asserted: the two changed files
were replaced with their `HEAD` contents, the whole suite was re-run, and the
changed files were restored and `cmp`-checked against the pre-swap copies. The
three `test_ingest.py` failures named above are the same three the previous
pass recorded — Xenium section registration reading `n_cells` as `None`.

---

## 11. Claims elsewhere that these measurements bear on

* **`results/scref_signatures_2026-09-20.md` §4.1 and §7 — CONFIRMED in
  magnitude, CORRECTED in cause and in consequence.** Its +633.7 MiB for
  `load_sample` reproduces (644.2 MiB here, three cold runs). Its attribution
  is where it goes wrong by implication rather than by statement: the staged
  line reads as a loader that is expensive to *parse*, and §1.1 above shows
  `mmread` is near-optimal and the cost was three conversions after it, two of
  them removable without touching an arithmetic operation. Its judgement that
  "the honest target for that file is a smaller overrun, not a fix" is
  **correct and is now quantified** — §7's floor is 531.5 to 707.2 MiB and its
  "~600 MiB" estimate sits inside that range.
* **Same file's `MEASURED_VISIUM_PEAK_MIB` table (922.2 to 1227.9 MiB) —
  SUPERSEDED, not contradicted.** It measured the code as it then was; §6
  measures the same sections after this change. Both are honest. Do not average
  them and do not take a ratio across them; every ratio in this file is
  within-harness. That file is a dated published record and is **not edited**;
  this is the successor.
* **Same file's §4.1 claim that `build_signatures` "is now a transient
  dominated by `load_reference()` itself" — CONFIRMED and sharpened.** Its
  staged figure of +507.8 MiB for `build_signatures` matches the +527.5 MiB
  measured here; what it could not see is that the same call cost only +94.3
  MiB of *high-water* while the old Visium loader was leaving arenas behind for
  it. The term did not change. What changed is that it is now paid for.
* **`results/scref_memory_2026-09-20.md` §3's "`load_sample` alone costs +400
  to +620 MiB" — CONFIRMED** (410.8 to 644.2 here, before).
* **Same file's §4.5 — STILL OPEN and untouched.** See §8.
* **`docs/real_data/visium_gbm.md`** states that `expression` is "a scipy CSR
  of shape `(n_spots, n_genes)`, integer UMI". Still true, and §5.1 pins every
  part of it.
