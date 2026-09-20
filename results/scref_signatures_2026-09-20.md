# The four densifications, the bit-exactness that survived them, and the one row that does not

**Date:** 2026-09-20
**Scope:** `tmesim/reconstruct/deconvolve.py` — `build_signatures` and the
spot loop in `deconvolve_sample`.
**Status of the headline claim:** the reported defect is REAL, the fix is DONE
and it is **bit-identical**: all four Visium sections produce `coords` and
`fractions` that are `np.array_equal` to the pre-change output, and the
measured cold peak falls **2.2x to 3.3x**. It does **not** make a Visium
section runnable on a 512 MB instance, and the refusal stays.

This is the piece of work that `results/scref_memory_2026-09-20.md` §7 named
and declined:

> `build_signatures` densifies a ~19,000x3,589 submatrix to float64 TWICE
> (~1.1 GB). It is now the dominant term and the only thing that could make
> Visium fit. A sparse mean is not trivially bit-identical in floating point,
> and I would not ship that half-verified. **This is the next piece of work.**

Four things that file did not know are established below: the term was **not**
~1.1 GB but ~1.8 GB; there were **four** densifications, not two, and two of
them were in a different function; **two of the four were removable without
touching the arithmetic at all**; and the mean can be computed exactly, but
the exactness has a **condition that a plausible implementation violates**.

---

## 0. Method, and what the numbers are numbers OF

Peak memory is `resource.getrusage(RUSAGE_SELF).ru_maxrss`.

> **`ru_maxrss` units are platform-dependent.** On **darwin/macOS it is
> BYTES**; on Linux it is KiB. Every figure here was taken on darwin and
> divided by 1024² **once**, not twice. Current RSS is `ps -o rss= -p <pid>`,
> which on darwin reports **KiB**, so that one is multiplied by 1024 first.
> `psutil` is not installed in this venv.

"Peak" means `ru_maxrss_after − ru_maxrss_before` around the call being
measured. `ru_maxrss` is a process high-water mark and never decreases, so
this is a *rise*, and ordering within a process matters — every headline
figure is therefore a **cold subprocess**, one measurement per process.

Host: darwin 24.6.0, 16 GiB, CPython 3.9, pandas 2.3.3, numpy 2.0.2, scipy
1.13.1. Scripts are in the session scratchpad.

**Before/after are both measured here**, with the same harness, in the same
session, on the same host. The "before" state is the code as of the start of
this pass — i.e. the code `results/scref_memory_2026-09-20.md` left behind —
loaded into the live module by path so that every caller sees it.

*The run-to-run spread is large and is not hidden.* Three cold "before" runs
of `GBM2` gave 3,099.8 / 1,703.3 / 2,457.7 MiB. That is the allocator not
returning freed arenas, the same effect that file documented in its §1. Every
threshold below is therefore the **maximum** observed, and every ratio is
computed **within** this harness. Ratios across harnesses are not meaningful
and are not taken.

---

## 1. What was actually there, measured step by step

`build_signatures` on the real panel (17,828 shared genes — identical for all
four built-in sections — against the 23,460 x 3,589 int32 reference), with RSS
and `ru_maxrss` marked between statements:

```
load_reference                 rss  678.3 MiB   maxrss  678.3
  expr_df.loc[~dup]            rss  999.6       maxrss  999.6   (+321.3)
  ref_sub[cells_kept]          rss 1320.9       maxrss 1320.9   (+321.3)
  .loc[found].to_numpy(f8)     rss 1809.1       maxrss 2053.2   (+732.3)
  _normalize_columns_to_unit   rss 1703.5       maxrss 2503.5   (+450.3)
  per-type means               rss 1793.2       maxrss 2503.5   ( +0.0)
```

**Peak rise over the post-`load_reference` state: 1,825.2 MiB**, not the
~1.1 GB the previous file estimated. The estimate counted `ref_found` and
`ref_norm` (488.2 MiB each) and stopped there.

### 1.1 Two of the copies select everything

The two 321.3 MiB copies at the top are whole-matrix int32 duplicates of the
reference, and on this data **both are no-ops**:

| | measured |
|---|---|
| `ref_index.duplicated().sum()` | **0** — nothing to de-duplicate |
| cells with a mappable label | **3,589 of 3,589** — nothing to drop |

So `expr_df.loc[~dup]` and `ref_sub[cells_kept]` copied 642.6 MiB to select
100 % of the rows and 100 % of the columns. Removing them cannot change a
value, because they change no value. This is the "a copy that exists only to
change dtype or layout is free to remove" case, and it was worth looking for.

### 1.2 The normalization held three full-size arrays, not two

```python
out = np.zeros_like(mat, dtype=np.float64)      # 488.2 MiB
out[:, nz] = mat[:, nz] / colsums[nz]           # + the mask copy, + the quotient
```

`mat[:, nz]` is a boolean-fancy-index copy (488.2 MiB) and the quotient is
another (488.2 MiB), both alive at once, on top of `mat` and `out`. Rewritten
as `np.divide(mat, colsums, out=out, where=nz)`, two of the four disappear and
every element is still the same single IEEE-754 double division. Verified
`np.array_equal` against the slice form on real data, zero-count columns and
input-not-mutated included.

### 1.3 The larger pair was in a different function

`results/scref_memory_2026-09-20.md` did not mention this one:

```python
expr_sub = expr[:, cols].toarray().astype(np.float64)
```

`expr` is an **int64** CSR. `toarray()` materialises the whole spot x gene
matrix as int64 and `.astype` materialises it again as float64 — **648.0 MiB
each for GBM1** (4,764 spots x 17,828 genes), 1,296 MiB together, immediately
before a loop that reads one row at a time. The CSR holds 13,576,167 nonzeros
in the selected submatrix: 155.4 MiB of data+indices, 91 % of the dense bytes
being zeros.

---

## 2. What was changed

### 2.1 `build_signatures` — streamed in gene-row chunks, two passes

The row/column selections are resolved to **positions** once
(`_positional_selection`), and the matrix is then walked in chunks of
`_SIGNATURE_CHUNK_ROWS` rows via `DataFrame.take` — chosen over
`expr_df.to_numpy()[rows]` because `to_numpy()` is a free view *only* on a
single-block frame and silently copies the whole matrix otherwise.

Pass 1 accumulates the per-cell column sums in **int64**. Pass 2 re-reads each
chunk, normalizes it in place, and reduces it straight into the `(G, 4)`
result. The 488 MiB `ref_found`/`ref_norm` pair is never allocated.

Chunk size was measured, not guessed — cold peak of `build_signatures` alone,
four runs per value, maximum shown (each includes `load_reference`'s own
~610 MiB floor):

| chunk rows | 250 | 500 | 1,000 | 2,000 | 4,000 | 8,000 | dense |
|---|---|---|---|---|---|---|---|
| peak (MiB) | 579.0 | 620.1 | **601.9** | 684.3 | 786.1 | 1,064.9 | 1,612.2 |
| wall (s) | 4.05 | 4.13 | 4.23 | 4.50 | 3.98 | 4.26 | 4.20 |

**Every value in that sweep, including the dense fallback, produced a
bit-identical signature matrix.** 1,000 is shipped.

### 2.2 `deconvolve_sample` — one spot at a time, into a reused buffer

`iter_spot_vectors(csr, n_genes)` yields each CSR row scattered into a single
reused float64 buffer, which `deconvolve_spot` consumes immediately. Both
648 MiB arrays are gone; the NNLS loop's measured contribution is now
**+0.0 MiB**.

`sum_duplicates()` is called first, and it is load-bearing rather than
defensive: `toarray()` **sums** repeated `(row, col)` entries in a
non-canonical CSR and a scatter assignment keeps only the last. A test pins
that, and the first version of that test was worthless — see §5.

### 2.3 The dense path is kept, and is still the definition

`_signatures_dense` is the original whole-matrix computation, arithmetic
verbatim. `_signatures_streamed` returns `None` and hands over to it whenever
the exactness argument does not apply. All four conditions are checked per
call, not assumed:

| condition | why | verified |
|---|---|---|
| integer dtype | the exactness proof is about integers | float frame -> dense, same answer |
| no negative entry | partial sums are otherwise not bounded by the total | negative entry -> dense, same answer |
| every column sum < 2**53 | otherwise float64 cannot hold it exactly | measured max 1,445,849; headroom 6.2e9x |
| unique column labels | `ref_sub[cells_kept]` would otherwise fan out | duplicate columns -> dense |

---

## 3. Bit-exactness: what was checked, and what it cost to be sure

### 3.1 The result

`deconvolve_sample(sample_id)` before vs after, cold process each, compared
with `np.array_equal` on both returned arrays:

| section | spots | `coords` | `fractions` | `type_names` |
|---|---|---|---|---|
| `GBM1` | 4,764 | equal | **equal** | equal |
| `GBM2` | 2,891 | equal | **equal** | equal |
| `GBM3` | 1,775 | equal | **equal** | equal |
| `10x_parent_gbm` | 3,468 | equal | **equal** | equal |

Not `allclose`. `array_equal`, max|Δ| exactly 0.0, re-run against the final
shipped state of the file.

The `reference=` caller shape used by
`tmesim/validation/reconstruction_fidelity.py` — a shuffled 70 % subset of
cells — was checked separately: streamed equals dense, and the result differs
from the full-reference signature, so the check is not comparing a thing to
itself.

### 3.2 The two reorderings, and which one is actually proved

**Column sums, chunked int64 instead of one float64 pass.** This one is
proved. The counts are non-negative integers, so every partial sum in any
order is an integer; an integer below 2**53 is exactly representable as a
float64; a sum whose every partial result is exact carries no rounding error;
and a computation with no rounding error cannot depend on its order. Measured
on the real file: `np.array_equal(chunked_int64.astype(f8), dense_f8_sum)` is
True, largest column sum 1,445,849, headroom to 2**53 a factor of 6.2e9.

*A plant confirmed the proof rather than the implementation.* Accumulating the
column sums in chunked **float64** instead — the change expected to go red —
left the headline test **green**. It is exact for exactly the same reason. The
int64 accumulator is kept because it makes the premise checkable, not because
float64 was wrong. That negative control is recorded in the test docstring.

**The per-type mean, chunked over gene rows.** The argument is that the
reduction runs along the cell axis while the chunking splits the gene axis, so
no individual reduction is reordered. That argument is sound **and it is not
sufficient**, which is the sharpest thing this pass found.

### 3.3 The one row that breaks it

`block[:, idx]` is an advanced-index copy, and **numpy 2.0.2 chooses its
memory order from the shape**:

| rows in the slice | `C_CONTIGUOUS` of `mat[0:M][:, idx]` | `.mean(axis=1)` reduces |
|---|---|---|
| 1 | **True** | pairwise, over a contiguous run |
| 2, 3, 8, 64 | **False** | sequentially, along a strided axis |

Same values, same operand order, different summation tree. So a chunk of
exactly one row silently computes a different number from every other chunk
size — and, incidentally, the *more* accurate one, which is no comfort when
the requirement is to reproduce the published output.

MEASURED, forcing one-row chunks over the full 17,828-gene panel:

| | vs the shipped chunking |
|---|---|
| signature entries differing | **32,568 / 71,312 (45.7 %)** |
| max absolute deviation | 2.429e-17 |
| max relative deviation | 3.597e-15 (7 ulp) |
| **GBM3 spot fractions differing** | **4,966 / 7,100 (69.9 %)** |
| max absolute deviation (fractions) | 2.331e-15 |
| max relative deviation (fractions) | 5.966e-13 |

**It propagates.** A deviation in the 16th significant figure of the
signatures becomes a deviation in 70 % of the published fractions. No
tolerance-based check would ever have seen it, which is exactly the failure
mode the brief for this work was written against.

`_chunk_bounds` closes it structurally: the floor is **two** rows, not one,
and a lone trailing row (`n % step == 1`) is **absorbed** into the previous
chunk rather than left alone. That second case is the quiet one — with the
shipped 1,000 and the shipped 17,828-gene panel the remainder is 828, so a
panel of 17,001 genes would have been the first anyone heard of it. `n == 1`
is the single legitimate one-row chunk, because the dense path is then a
one-row matrix too and makes the same choice.

### 3.4 What this does and does not entitle anyone to claim

Bit-exactness here is **verified against the dense definition on the installed
numpy**. It is not proved for all numpy versions, and §3.3 is a concrete
demonstration that the property is contingent on a library heuristic rather
than on the arithmetic. The standing guard is
`tests/test_signature_streaming.py::test_streamed_signatures_are_bit_identical_to_the_dense_definition`,
which recomputes both paths on the real reference in the same process and
compares bits. If numpy's heuristic ever changes, that test goes red instead
of the results drifting quietly — which is the defect this project keeps
finding in itself.

---

## 4. Memory, after

Cold `deconvolve_sample`, one per process. Three runs before, five or six
after; maximum of each, same harness throughout.

| section | before | after | ratio |
|---|---|---|---|
| `GBM1` | 2,716.8 | **1,227.9** | 2.2x |
| `GBM2` | 3,099.8 | **939.0** | 3.3x |
| `GBM3` | 2,567.3 | **922.2** | 2.8x |
| `10x_parent_gbm` | 2,923.5 | **1,215.7** | 2.4x |

Wall time is slightly **better**, not worse: GBM1 6.48 s -> 5.97 s, GBM3
5.55 s -> 5.00 s.

### 4.1 Where it goes now, and who owns it

Staged, cold, `GBM1`:

```
imports                     109.7 MiB peak
load_sample(visium)         743.4      (+633.7)   <-- now the largest term
reference_gene_symbols      772.7      ( +29.3)
build_signatures           1280.8      (+507.8 transient, 184 resident)
CSR column subset + NNLS   1280.8      (  +0.0)
```

**The dominant term has moved to `load_sample`**, in
`tmesim/data/real/visium.py`, which this pass does not own. `build_signatures`
is now a transient dominated by `load_reference()` itself, and the spot loop
costs nothing measurable.

### 4.2 The 512 MB verdict is unchanged, and the refusal thresholds were lowered

A Visium section still needs **0.9 to 1.2 GB**. The gap to a 512 MB container
is now **1.7x to 2.4x** instead of ~6x, but it is a gap, and
`MEASURED_VISIUM_PEAK_MIB` was re-measured down to match rather than left
stale — a refusal threshold 2.5x above the real cost would refuse hosts that
can now serve a section. `MEASURED_VISIUM_MIN_MIB` is 869.9, the cheapest run
ever observed here (GBM2), and is deliberately the minimum, not the mean.

The refusal message was also corrected: it used to tell the reader the cost
was "dominated by densifying the scRNA reference to float64 in
build_signatures", which was true when it was written and is not true now.

---

## 5. Plant evidence

Every test in `tests/test_signature_streaming.py` was confirmed RED by
planting the specific defect it targets, confirming the plant parses and
imports, running the test, and reverting. **Two of the tests were vacuous on
the first attempt and the plants are what caught them.**

| plant | test | result |
|---|---|---|
| `block *= 1.0/cs` instead of `np.divide(..., where=)` | streamed == dense | **RED** `array_equal` False |
| chunked **float64** column sums instead of int64 | streamed == dense | **GREEN — negative control**, see §3.2 |
| `_chunk_bounds` floor lowered 2 -> 1 | chunking is a memory knob | **RED** chunk 1 diverges |
| orphan-row absorption deleted | chunk bounds invariant | **RED** `lone row at n=3 step=0: [(0,2),(2,3)]` |
| orphan-row absorption deleted | chunking is a memory knob | **RED** (only after the test was strengthened — see below) |
| `_EXACT_INT_LIMIT` 2**53 -> 2**20 | int64 column sums | **RED** `(1048576.0 / 487894.0) > 1000000.0` |
| integer-dtype guard removed | premises / `float_dtype` | **RED** `streamed path ran on an invalid premise` |
| non-negative guard removed | premises / `negative_value` | **RED** same assertion |
| `cols.is_unique` guard removed | premises / `dup_columns` | **RED** `pandas InvalidIndexError` out of `get_indexer` |
| `_positional_selection` always `(None, None)` | valid toy takes streamed path | **RED** `dense fallback taken` |
| last-occurrence-wins for duplicate symbols | first occurrence wins | **RED** all 40 rows moved, max rel **2.6** |
| `where=nz` dropped from `np.divide` | normalize matches slice form | **RED** `invalid value encountered in divide`, NaN |
| `sum_duplicates()` removed | spot scatter matches `toarray` | **RED** row 1 `[9,0,4,0]` vs `[9,0,11,0]` |
| buffer reset deleted | buffer cleared between rows | **RED** spot 1 `[0,5,7,0]` |
| peak table -> 100.0 MiB | still refuses on 512 MB | **RED** `DID NOT RAISE` |
| `MEASURED_VISIUM_MIN_MIB` -> 100.0 | still refuses on 512 MB | **RED** `DID NOT RAISE` |
| `need <= have` early exit removed | a host with room is not refused | **RED** raised on an 8 GB host |

### 5.1 The two vacuous checks this pass wrote, and then caught

Recorded rather than quietly fixed, because the point of the exercise is that
they are easy to write:

1. **`test_spot_scatter_matches_toarray_including_duplicate_entries`** built
   its matrix with `coo_matrix(...).tocsr()`. That conversion **canonicalises
   on the way**, so the duplicate entry never survived to the call and
   deleting `sum_duplicates()` left the test green. It now builds the CSR from
   raw `indptr/indices/data` and asserts `not m.has_canonical_format` before
   calling.

2. **`test_the_measured_table_still_refuses_every_section_on_512_mb`** was
   first planted with thresholds of 400.0 MiB and stayed green — because
   400 MiB still exceeds the headroom a 512 MB host has once the test process
   is ~156 MiB resident, so the refusal fired anyway. A threshold plant has to
   clear the **headroom**, not the limit. Re-planted at 100.0.

3. A third near-miss: deleting the orphan-row absorption left the chunk-size
   test green on the 1,500-gene panel, because the orphaned row there happens
   to be one of the ~18 % whose value does not move. A three-gene case where
   the row is known to move was added.

### 5.2 Exit codes

```
tests/test_signature_streaming.py      15 passed        (exit 0)
tests/test_scref_memory.py             15 passed        (exit 0)
tests/test_gbmspace_adapter.py         included below   (exit 0 for this file)
tests/test_ingest.py                    3 failed, 45 passed, 18 skipped (exit 1)
```

The three `test_ingest.py` failures are **pre-existing and unrelated**:
confirmed identical against the unmodified `deconvolve.py` by swapping the
file back and re-running. They concern Xenium section registration
(`n_cells == 40887` reading `None`) and this worktree's `data_cache` does not
carry the Xenium exports.

---

## 6. Claims elsewhere that these measurements bear on

* **`results/scref_memory_2026-09-20.md` §7 — CONFIRMED in direction,
  CORRECTED in magnitude and in scope.** Its call ("this is the next piece of
  work and it is the one that could actually make Visium fit") was right that
  the term was dominant and right to refuse to ship it half-verified. Its
  figure of ~1.1 GB undercounts the real 1,825 MiB, it names two
  densifications where there are four, and the two it names are the smaller
  pair. Its prediction that this "could make Visium fit" in 512 MB is **not
  borne out**: the peak fell 2.2x-3.3x and is still 1.7x-2.4x over the limit.
  That file is a dated published record and is **not edited**; this is the
  correction.
* **Same file's `MEASURED_VISIUM_PEAK_MIB` table (3,146-3,455 MiB).** My own
  "before" measurements of the same code are 13-21 % lower (2,567-3,100 MiB).
  Both are honest; the run-to-run spread within a single section here reaches
  1,400 MiB. **Do not average the two tables and do not take a ratio across
  them.** Every ratio in this file is within-harness.
* **Same file's §4.5 — STILL OPEN and untouched.** `app/server/jobs.py`'s
  `RUN_REFUSALS` still does not name `DeconvolutionTooLargeForHost` or
  `ReferenceTooLargeForHost`, so the visitor still sees the redacted sentence.
  This pass does not own `app/` and did not change it. The thresholds those
  refusals quote are now lower, which makes the message more accurate but no
  more visible.

---

## 7. Not fixed, and why

* **`load_sample`'s +633.7 MiB, now the largest single term.** In
  `tmesim/data/real/visium.py`, outside this pass's ownership. It is the next
  piece of work by the same argument that made `build_signatures` this one —
  with the caveat that removing all of it still leaves ~600 MiB, which is
  still over 512 MB. **Nothing measured here suggests a Visium section can be
  made to fit in 512 MB; the honest target for that file is a smaller
  overrun, not a fix.**
* **`load_reference()`'s ~610 MiB peak / 322 MiB resident.** Untouched; it is
  the floor under `build_signatures` and was the subject of the previous pass.
* **`tmesim/validation/reconstruction_fidelity.py`** still calls
  `load_reference()` to build a symbol set, as `deconvolve_sample` used to.
  Not on the live request path, and not measured here.
* **The dense fallback is never exercised by the shipped configuration.** It
  runs only under the four premise failures, which tests construct
  deliberately. It is kept because it is the definition of the answer, and it
  costs 1,612 MiB if it ever runs.
* **Linux figures.** Everything here is darwin. Array sizes are
  platform-independent; allocator behaviour is not. The numpy memory-order
  heuristic in §3.3 is a numpy property rather than a platform one, but it was
  observed only on numpy 2.0.2 and is not assumed to hold elsewhere — which is
  why the guard is a test and not a comment.
