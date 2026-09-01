# Stars density subset: data race on `cell->hydro.sort` → out-of-bounds crash

**Status:** root cause identified via instrumented run AND fix validated by A/B test
(FLAMINGO, grammar cluster, single node). See §6.
**Symptom:** segfault in `runner_dosub_pair_subset_stars_density` /
`runner_dopair_subset_stars_density` (via `runner_do_stars_ghost`), reading a gas
particle at an out-of-range index taken from a gas cell's hydro sort array.
**Fix:** hold `cj->hydro.extra_sort_lock` across the sorted-subset read in
`DOPAIR1_SUBSET_BRANCH_STARS` (patch `0002-stars-only-hold-extra_sort_lock.patch`).

---

## 1. Crash signature

The star density loop walks a **gas** cell `cj`'s sorted particle list (near→far) and
fetches each gas particle by the stored sort index:

```c
// runner_doiact_functions_stars.h  (DOPAIR1_SUBSET_STARS, ~line 823/846/909)
const struct sort_entry *restrict sort_j = cell_get_hydro_sorts(cj, sid);
...
struct part *restrict pj = &parts_j[sort_j[pjd].i];   // sort_j[pjd].i is garbage → OOB
if (part_is_inhibited(pj, e)) continue;               // <-- segfault here
```

`sort_j[pjd].i` is expected to be in `[0, count_j)`. In the crash it is ~5400 while
`count_j` is only a few hundred → out-of-bounds read into `parts_j`.

Only **stars** crash, and only when the subset density is driven from the **stars
ghost** (`runner_do_stars_ghost` → `runner_dosub_pair_subset_stars_density`), which
re-runs a subset of the density interaction **directly, outside the task graph**.
Equivalent hydro pair code never crashes.

---

## 2. Instrumented run — guard and data

A guard was inserted right after `pj = &parts_j[sort_j[pjd].i];` in
`DOPAIR1_SUBSET_STARS`, printing the raw index and the cell's sort state whenever
`sort_j[pjd].i >= count_j`:

```
grep STARDENS-BAD output-457134.txt
[0000][00484.1] STARDENS-BAD id=32030723 pjd=0 count_j=762 sort_i=5393 sid=4 sorted=10f0 sort_alloc=18f0 dx_old=0.0 dmin=5.625e-01
[0000][00484.1] STARDENS-BAD id=32030544 pjd=0 count_j=762 sort_i=5393 sid=4 sorted=10f0 sort_alloc=18f0 dx_old=0.0 dmin=5.625e-01
[0000][00639.6] STARDENS-BAD id=32280375 pjd=0 count_j=308 sort_i=5392 sid=4 sorted=1050 sort_alloc=1250 dx_old=0.0 dmin=5.625e-01
...
```

### Decoding the flags (sid = 4)

| case | `sorted`  | bits set        | `sort_alloc` | bits set             | extra alloc bit |
|------|-----------|-----------------|--------------|----------------------|-----------------|
| 484s | `0x10f0`  | {4,5,6,7,12}    | `0x18f0`     | {4,5,6,7,**11**,12}  | **11**          |
| 639s | `0x1050`  | {4,6,12}        | `0x1250`     | {4,6,**9**,12}       | **9**           |

Key observations:

1. **Gate is legitimately satisfied.** Bit `sid=4` is set in *both* `sorted` and
   `sort_allocated`; `dx_old = 0 ≤ space_maxreldx·dmin`. The branch check that decides
   "use the sorted interaction" is correct at the time it runs.
2. **Address is computed correctly.** `cell_get_hydro_sorts` picks the block with
   `j = popcount(sort_allocated & ((1<<sid)-1))`. For `sid=4` the low bits are 0 in both
   cases → `j = 0`. So the star reads `sort[0]`, the *first* block, *first* entry — the
   right place. It is **not** an offset/flag-inconsistency bug.
3. **The buffer *content* is stale/foreign.** `sort[0].i ≈ 5400` is stable across
   different cells, different time steps and different `count_j` (762 vs 308). A stable
   ~5400 is not per-cell random garbage — it is the content of *some other* cell whose
   `count ≈ 5400`, i.e. a **recycled allocation**.
4. **`sort_allocated` has exactly one extra bit over `sorted`** (11, then 9): a direction
   that is *allocated but not yet sorted*. This is the fingerprint of
   `cell_malloc_hydro_sorts` having **just grown** `cj`'s sort meta-array to add a newly
   requested direction.

This rules out the two "static inconsistency" hypotheses (sorted-set-but-not-allocated;
wrong `j` offset; `dx_old` staleness) and points squarely at a **concurrent
reallocation** of `cj->hydro.sort`.

---

## 3. Root cause — lock-free read races the sort-array realloc/free

### The reader (stars subset) drops the lock before reading

`runner_doiact_functions_stars.h`, `DOPAIR1_SUBSET_BRANCH_STARS` (~line 1198):

```c
lock_lock(&cj->hydro.extra_sort_lock);
const int is_sorted = (cj->hydro.sorted & (1<<sid)) &&
                      (cj->hydro.dx_max_sort_old <= space_maxreldx * cj->dmin);
lock_unlock(&cj->hydro.extra_sort_lock);          // <-- unlock
...
if (force_naive || !is_sorted)
  DOPAIR1_SUBSET_STARS_NAIVE(...);
else
  DOPAIR1_SUBSET_STARS(...);                       // <-- reads cj->hydro.sort with NO lock
```

The lock only protects the *decision*; the actual dereference of `cj->hydro.sort`
happens after the unlock.

### The writer frees and swaps the buffer under the lock

`runner_do_hydro_sort` (`runner_sort.c:220, 279`) holds `extra_sort_lock` and calls
`cell_malloc_hydro_sorts`. When a **new** direction is requested, that function
(`cell.h:1400-1424`):

```c
new_array = swift_malloc(..., num_arrays_wanted*(count+1));
memcpy(new_array, c->hydro.sort, ...);        // copy existing directions
swift_free("hydro.sort", c->hydro.sort);      // <-- FREE the old buffer
c->hydro.sort = new_array;                     // <-- swap the pointer
```

### The race

```
star thread                          sort thread (on same cj)
-----------                          ------------------------
lock; read sorted(sid=4)=set; unlock
                                     lock (extra_sort_lock)
                                     cell_malloc_hydro_sorts(cj, newdir 11/9):
                                       swift_free(cj->hydro.sort)   <-- old buffer freed
                                       cj->hydro.sort = new_array
                                     ... fill/sort new dir; unlock
cell_get_hydro_sorts(cj,4)
  reads cj->hydro.sort[0].i  ------> freed-then-recycled memory (a ~5400-part cell)
  => sort_i = 5393 >> count_j = 762
  => parts_j[5393] OOB => segfault
```

The freed block is handed back to the allocator and reused by another cell of the same
size class (count ≈ 5400), which is why the bad index is a stable ~5400 rather than
random.

### Why only stars, and only from the ghost

Normal hydro pair tasks are protected by the **task-graph dependency**: the sort task is
a predecessor of the pair task, so a cell is never being sorted while its pair
interaction runs — the same "unlock then read" pattern in the hydro path can never race.

The stars ghost (`runner_do_stars_ghost`) re-runs the subset density **directly, bypassing
the task graph** (`runner_dosub_pair_subset_stars_density`), so there is **no sort
dependency guarding `cj`**. Another runner thread is free to sort/realloc `cj`'s hydro
sort concurrently. This is why the bug is stars-only and reproduces even on a single node.

### Role of particle splitting

`SPH:particle_splitting` is **not** the direct cause. Splitting forces a rebuild → many
cells acquire new neighbours → many sort directions are requested for the first time →
`cell_malloc_hydro_sorts` takes the free+swap branch far more often → the race fires much
more frequently. Disabling splitting makes the run survive longer but does not remove the
underlying data race.

---

## 4. Suggested fixes (in order of increasing scope)

1. **Extend the critical section (minimal) — IMPLEMENTED & VALIDATED (see §6).** In
   `DOPAIR1_SUBSET_BRANCH_STARS`, hold `cj->hydro.extra_sort_lock` across the
   `DOPAIR1_SUBSET_STARS` call, not just the gate check, so the sort array cannot be
   freed/swapped while it is being read. The naive branch (which does not touch the sort
   array) releases the lock early, so the extra hold is confined to the sorted read.
   Downside: slightly longer lock hold on the sorted path, minor concurrency loss in
   clustered regions only. Patch: `0002-stars-only-hold-extra_sort_lock.patch`.
2. **Don't free in place.** Make `cell_malloc_hydro_sorts` retire the old buffer via
   deferred/epoch reclamation instead of an immediate `swift_free`, so an in-flight
   reader holding the old pointer stays valid. Larger change.
3. **Fix it at the graph level (cleanest).** Give the stars-ghost subset density a proper
   dependency on `cj`'s hydro sort so it never runs concurrently with a sort of `cj`,
   matching the guarantee the normal hydro/stars pair tasks already have.

The same latent unlock-then-read pattern also exists in the **hydro** subset branch
(`DOPAIR1_SUBSET_BRANCH` in `runner_doiact_functions_hydro.h`). It does not crash in
practice because the hydro subset only ever requests already-allocated sort directions
(no grow → no `swift_free`), but for correctness completeness fix (1) should be applied
there too before an upstream MR.

---

## 5. Reproduction / environment

- Model: FLAMINGO, SPHENIX hydro, `--with-subgrid=FLAMINGO`, single node.
- Cluster: grammar (Mellanox IB, OpenHPC gnu14/openmpi5). See
  `[[grammar-cluster-swift-build]]` / `[[grammar-cluster-swift-launch]]`.
- Guard: temporary `STARDENS-BAD` print inserted after
  `pj = &parts_j[sort_j[pjd].i];` in both occurrences in
  `src/runner_doiact_functions_stars.h`, firing when `sort_j[pjd].i >= count_j`,
  dumping `id, pjd, count_j, sort_i, sid, sorted, sort_alloc, dx_old, dmin`. The guard
  prints and `continue`s (non-fatal), so it doubles as a detector: with the bug present it
  logs a bad index instead of crashing; with the fix present it should log nothing.
- Log: `output-457134.txt` (unpatched, bug present).

---

## 6. Validation (A/B test)

Fix (1) was applied (`0002-stars-only-hold-extra_sort_lock.patch`) with the `STARDENS-BAD`
detector guard *left in place*, and the identical FLAMINGO reproduction was rerun with the
only changed variable being the patch (same IC, particle splitting still enabled).

| build | run | result |
|-------|-----|--------|
| unpatched (bug present) | `output-457134.txt` | `STARDENS-BAD` fires reliably and repeatedly at wall-clock **~484 s** and **~639 s**, then crashes |
| **patched (fix 1)** | `output-457155.txt` | ran to wall-clock **≥ 64 307 s (~17.8 h)**, `grep -c STARDENS-BAD` = **0**, no crash |

The patched run passed the previous reliable crash points by a factor of ~100× with **zero**
`STARDENS-BAD` occurrences, while every other variable was held constant. This confirms the
data race on `cj->hydro.sort` is the cause and that holding `extra_sort_lock` across the
sorted-subset read removes it.

**Remaining work:** remove the temporary `STARDENS-BAD` guard and prepare the upstream
SWIFTSIM MR using this evidence. See `UPSTREAM_MR_SORT_UAF.md` for packaging.

---

## 7. The deeper shared root cause — writer-side fix (also validated)

The `extra_sort_lock` reader fix above stops the **stars** reader from racing the free, but
the **free itself** (`swift_free` of the old buffer in `cell_malloc_hydro_sorts`,
`cell.h`) still happens — so any *other* lock-free reader can still hit a use-after-free.
Exactly that turned up next: the EAGLE time-step **limiter** (`runner_dopair1_branch_limiter`)
is a second unprotected reader of the same buffer and crashed with the identical fingerprint
(`sort_idx ≈ 5400`). Details and validation in `LIMITER_MPI_SORT_CRASH_ISSUE.md`.

The unifying fix is **writer-side**: don't free the grown-over sort buffer in place; retire
it and free it at the next rebuild barrier (`0003-writer-side-defer-hydro-sort-free.patch`).
Validated on EAGLE-XL 4-node (`output-468708.txt`): **all** `LIMITER-BAD` and `STARDENS-BAD`
= 0, no crash, ran to z≈2.07 / ~17.8 h. One change covers stars, limiter, and every other
lock-free sort reader.

Relationship of the two fixes (both kept for now):
- **writer-side (`0003`)** removes the use-after-free for *all* readers — the primary fix.
- **stars `extra_sort_lock` read (`0002`)** additionally covers the residual, very narrow
  *grow half-update* window (a reader seeing the new `sort_allocated` bit but the old
  `sort` pointer between the two stores in `cell_malloc_hydro_sorts`). Writer-side alone
  does not close this window; it is nanosecond-wide and usually reads wrong data rather than
  faulting. See `UPSTREAM_MR_SORT_UAF.md` for whether to ship one or both upstream.

> Note: the earlier "MPI foreign cell / rendezvous" hypothesis for the limiter crash
> (`UCX_RENDEZVOUS_LIMITER_ISSUE.md`) was **disproved** — the limiter guard showed all bad
> cells are local. The limiter crash is this same local sort use-after-free.
