# GBM-Space ingestion, Batch C — 0 of 4 claimed sections ingested, honest null

**Outcome up front**: Batch C claimed 4 sections and attempted to ingest them twice.
Both attempts failed identically before any data reached disk. No new sections were
added to `data_cache/gbmspace_visium/` by this batch. Claims were released; nothing
was left claimed-but-orphaned. This report is a diagnosis of a real, reproducible
extraction failure, not a completed ingestion — reported per the task's honesty
requirement rather than silently omitted or papered over.

---

## 1. Coordination steps followed

1. `data_cache/gbmspace_visium/` listed first: 12 sections already present (5 patients
   — AT10, AT11, AT14, AT15, AT6), matching `results/gbmspace_adapter.md`'s count at
   the start of this pass.
2. `data_cache/.gbmspace_claims/` checked (created — did not yet exist). 4 claims
   already present from sibling agents at first check
   (`AT12-BRA-4-FO-4_1`, `AT5-BRA-5-FO-1_2`, `AT7-BRA-3-FO-1_1`, `AT9-BRA-4-FO-1_1`).
3. Full atlas listed via a streamed `curl -s <url> | tar tzvf -` (zero bytes written
   to disk). This listing itself took several minutes and ended after enumerating
   **65 lines / 30 unique `.h5ad` members** (of the atlas's 97) — the connection
   ended before reaching the rest of the archive; this is normal for a sequential
   tar stream with no central directory, not an error, and was sufent to pick from.
4. Waited ~22s, re-checked claims (unchanged) and re-checked disk (17 GiB free,
   comfortably clear of the 10 GiB floor).
5. From the 30 candidates, removed the 5 already ingested and the 2 already claimed
   that appeared in this list, skipped the next 8 smallest, and claimed 4:
   **`AT5-BRA-5-FO-3_1`, `AT4-BRA5-FO-3_1`, `AT13-BRA-3-FO-1_2`, `AT5-BRA-5-FO-2_1`**
   (sizes 147/148/153/154 MB — patients AT5 and AT13, AT13 being a patient not yet
   represented anywhere in this repo). Claim markers written before any download,
   per protocol. (One process error along the way: the first claim attempt used an
   unquoted shell variable under zsh, which does not word-split unquoted expansions
   the way bash does — this produced one bogus combined-name claim file instead of
   4 separate ones. Caught immediately by listing the claims directory back, deleted,
   and redone correctly with an explicit `for c in a b c d` loop. Flagging this here
   since it's a real shell-portability trap this task's protocol should account for
   if any future batch scripts this the same way.)

## 2. The extraction failure — what actually happened

Ingestion follows the documented pattern exactly: one streamed pass,
`curl -s <atlas-url> | tar xzf - -C <tmpdir> anndata/<name1>.h5ad anndata/<name2>.h5ad ...`,
then `gbmspace.register()` per extracted file, then delete the raw `.h5ad`
immediately — the same method `results/gbmspace_adapter.md` §3/§9 describes.

**Both attempts (a fresh `curl` connection each time) failed identically**:

```
tar: anndata/AT5-BRA-5-FO-3_1: Not found in archive
tar: anndata/AT4-BRA5-FO-3_1: Not found in archive
tar: anndata/AT13-BRA-3-FO-1_2: Not found in archive
tar: anndata/AT5-BRA-5-FO-2_1: Not found in archive
tar: Error exit delayed from previous errors.
```

This is not a case of picking nonexistent names: all 4 exact member paths and their
byte sizes had just been independently confirmed present in the same archive via the
step-3 `tar tzvf -` listing (`anndata/AT5-BRA-5-FO-2_1.h5ad` = 153,717,609 bytes,
`anndata/AT13-BRA-3-FO-1_2.h5ad` = 153,032,020 bytes, etc.) — the "not found" is
bsdtar (macOS's default `tar`) reaching EOF of a stream that ended before it read as
far as these members, then reporting them as absent rather than reporting a truncated
read. Confirmed not a hung/stuck process either time: `ps` showed the `curl` and
`tar` processes accumulating real CPU time (17-18% sustained, ~30s CPU over ~4.5 min
wall) right up to exit, not blocked at 0% CPU.

**This does not look like a server-wide problem or a genuine absence of the data**:
during the same window, sibling batches (confirmed via `ps aux` showing their exact
`curl | tar xzf -C data_cache/_batchB_tmp/_batchD_tmp ...` commands running
concurrently, and confirmed via `data_cache/gbmspace_visium/` gaining 7 new
directories — `AT11-BRA-6-FO-4_2`, `AT12-BRA-4-FO-4_1`, `AT5-BRA-5-FO-1_2`,
`AT6-BRA-3-FO-1_2`... — over the course of this pass) were successfully extracting
their own named members via the identical pattern against the identical URL. The
most likely explanation, stated as a hypothesis rather than a proven cause: multiple
concurrent full-archive streams from this one machine against the same Sanger COG
endpoint appear to make some connections truncate before reaching their target
members' offset in the (non-indexed, sequentially-read) tar stream, while other
concurrent connections succeed — non-deterministic under load, not a fixed failure.
Both of my attempts happened to lose the race; I did not retry a third time with the
same or different targets because two identical failures plus visible heavy
concurrent contention (2-4 other `curl` processes to the same host at any given
check) was enough signal that a third attempt was unlikely to behave differently,
and the task instructed reporting an honest partial result over an indefinite retry
loop.

**This is reported as a real, precise wrinkle, not silently worked around**: no
change was made to `tmesim/data/adapters/gbmspace.py` or any other adapter/core file
to "fix" this — there was nothing in the adapter to fix; the failure is entirely in
the shell-level extraction step that happens *before* any adapter code runs on the
data.

## 3. What was NOT done, honestly

- **0 of the 4 claimed sections were ingested.** No `register()` call was ever
  reached for any of them (extraction failed before that point every time).
- **No reconstruction / hypoxic-fraction run was performed** — nothing new exists to
  reconstruct. The task's request for a reconstruction check on "one of your newly
  ingested sections" does not apply here since none were ingested.
- Claims for all 4 sections were **released** (`rm`'d from
  `data_cache/.gbmspace_claims/`) rather than left in place, so no section is
  permanently blocked from a future, less-contended attempt by any batch.
- Scratchpad temp extraction directory removed; no partial `.h5ad` bytes were left
  anywhere (extraction failed before writing any target file — `find
  .../gbmspace_batch_c_raw -type f` returned empty both times).

## 4. Checks that did pass

- Smoke import: `./.venv/bin/python -c "import tmesim.data.adapters.gbmspace"` →
  `smoke import OK`.
- No file under `tmesim/core`, `tmesim/physics`, `tmesim/sim`, `tmesim/molecular`, or
  `tmesim/data/adapters/gbmspace.py` itself was touched by this batch.
- No `git add`/`commit`/`push` performed.
- Full `pytest` suite was **not** run, per instruction (parallel-CPU-contention
  concern).

## 5. Disk state

| Point | Free space |
|---|---:|
| Start of this pass | 21 GiB |
| After claiming (pre-download) | 17-19 GiB (fluctuating with sibling batches' concurrent activity) |
| End of this pass | 19 GiB |

Never approached the 10 GiB floor. The fluctuation across the pass reflects other
batches' concurrent downloads/materializations/cleanups, not this batch's own usage
— this batch's peak disk footprint was ~0 bytes retained (every extraction failed
before any bytes were written to `data_cache/` or kept on disk).

## 6. Section count after this pass

`data_cache/gbmspace_visium/` holds **19 sections** at the end of this pass (12 at
the start of this pass + 7 added by sibling batches during the same window) — Batch
C itself contributed 0 of those 19.
