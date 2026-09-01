# Time-step limiter: use-after-free on `cell->hydro.sort` → out-of-bounds crash

**Status:** root cause identified AND fix validated. Same free+swap race as the stars bug in
`STARS_SORT_RACE_ISSUE.md`, with the **limiter** as a second, unprotected reader. The
writer-side deferred-reclamation fix (§5) removed **all** `LIMITER-BAD` occurrences,
including the "moderate" burst — so that burst was the same free-race, not a separate bug.
**Symptom:** segfault in `runner_dopair1_branch_limiter` → `DOPAIR1` (via
`runner_dosub_pair1_limiter`), dereferencing a gas particle at an out-of-range index taken
from a **local** cell's hydro sort array.
**Model:** EAGLE-XL, SPHENIX hydro, 4-node MPI. (Splitting OFF — irrelevant to the cause.)

Supersedes the earlier "foreign cell / MPI rendezvous" hypothesis: the instrumentation
below shows the bad cells are **local**, so `UCX_RENDEZVOUS_LIMITER_ISSUE.md` is **not** the
explanation for this crash.

---

## 1. Crash signature

```
[grammar005:...] Caught signal 11 (Segmentation fault: address not mapped ...)
 3  runner_dopair1_branch_limiter()
 4  runner_dosub_pair1_limiter.isra.0()
 6  runner_main()
```

The limiter pair (`DOPAIR1`, `src/runner_doiact_functions_limiter.h`) walks a cell's sorted
list and fetches each particle by `parts[sort[...].i]`. It passes the
`error("Interacting unsorted cells")` gate (so `hydro.sorted` is set) and then faults
dereferencing an out-of-range `sort[...].i`.

---

## 2. Instrumentation and data

A non-fatal `LIMITER-BAD` guard (plain `message()` + `continue`, **not** under `#ifdef
SWIFT_DEBUG_CHECKS`, so it works in the production build) was inserted before all four
sort-index dereferences in the limiter `DOPAIR1` (via `insert_limiter_guard.py`). It fires
when `sort[...].i < 0 || >= count` and dumps side, index, count, `sort_idx`, `sid`,
`ci_foreign`, `cj_foreign`, `sorted`, `sort_alloc`, `dx_old`, `dmin`. Being non-fatal it
also keeps the run alive (turns the crash into logged skips).

Result on a long 4-node EAGLE-XL run (`output-457813.txt`): **132 `LIMITER-BAD`, 0
crashes.** Key aggregates:

- `cj_foreign` — **132 / 132 = 0 (all local)**. Not a foreign/MPI-buffer problem.
- Ranks/steps — spread across **all 4 ranks** and many different steps → **systematic**, not
  a one-off.
- `sid` — 12 (×106), 3 (×14), 2/10 (×4), 7/9 (×2). Multiple directions.
- `count_j` — 92 (×106), plus 9,10,38,42,46,51,192 → multiple distinct cells.
- `sort_idx` — **two regimes:**
  - **moderate:** the `count_j=92`, `sid=12` burst (one step, one cell) with `sort_idx`
    ≈ 198–228 (~2.5× count), `sorted == sort_alloc == 0x1202`, `dx_old = 0`.
  - **wild:** `sort_idx = 5398, 5450, 1936877413 (≈1.9e9)`.

The **wild values match the stars fingerprint** (stars crashed with `sort_idx ≈
5263/5392/5393`): they are reads of a **freed-and-recycled** sort buffer belonging to some
other, larger cell. `1.9e9` is clearly reinterpreted non-sort memory.

---

## 3. Root cause — the limiter is an unprotected reader of the same freed sort buffer

### The writer that frees the buffer

`cell_malloc_hydro_sorts` (`src/cell.h` ~1400-1424), when a cell's sort array is **grown**
to add a newly-requested direction, does:

```c
new_array = swift_malloc(...);
memcpy(new_array, c->hydro.sort, ...);
swift_free("hydro.sort", c->hydro.sort);   // cell.h:1423 — frees the old buffer
c->hydro.sort = new_array;
```

These on-the-fly grows are triggered from the **ghost** phase (stars ghost / hydro subset
ghost requesting a direction not yet allocated).

### The limiter reads sorts with no protection

Task-graph audit (`src/engine_maketasks.c`):

- `hydro.sorts` unlocks only `t_rho`, `t_gradient`, `t_rt_*` (lines 818-820) — **not** the
  limiter.
- The limiter self/pair interaction tasks (created ~2667 / 2976) depend only on
  `hydro.super->hydro.drift` and `super->timestep`, and unlock `kick1` /
  `timestep_limiter` (lines 2910-2913, 3275-3278, 3448-3457). **No edge from
  `hydro.sorts`.**

So the limiter reads `cell->hydro.sort` **lock-free and without any sort dependency**,
relying on the density-phase sorts persisting untouched.

### The race

Different cells sit at different pipeline stages simultaneously. While cell X is still in
its density/stars **ghost** (triggering an on-the-fly grow → `swift_free` + pointer swap of
a shared gas cell G's sort), cell Y can already be at its **limiter**, walking G's sort. The
limiter reads the **freed/recycled** buffer → out-of-range `sort_idx` (the ~5400 / 1.9e9
values) → segfault.

This is the **same free+swap use-after-free as the stars bug** (`STARS_SORT_RACE_ISSUE.md`,
same `swift_free` at `cell.h:1423`). The stars fix (hold `extra_sort_lock` across the stars
read) only protects the **stars** reader; the free still happens, so the **limiter** — a
different unprotected reader — still gets caught.

### Why FLAMINGO never showed a limiter crash

FLAMINGO **does** run the limiter (`timestep_limiter=4096` in its own task counts). It never
surfaced because: (1) the stars race crashed FLAMINGO first at ~484 s — the stars ghost
reads the sort *in the same phase* the grow/free happens, so it catches the free far more
often than the later-running limiter; (2) after the stars fix, the 17 h FLAMINGO run had
**no** `LIMITER-BAD` guard, and a bad index does not always land in unmapped memory (can
read wrong data silently); (3) the limiter race is rarer and 4-node MPI + EAGLE-XL make it
more likely to trigger.

### The moderate burst was the same free-race (resolved)

The `count_j=92` / `sid=12` burst (`sorted == sort_alloc`, `sort_idx` ~2.5× count) was
initially suspected to be a separate count/sort bookkeeping bug. The §5 validation settles
it: with the writer-side fix applied, the moderate burst **also disappeared** (total
`LIMITER-BAD` = 0), so those moderate indices were likewise reads of a freed-then-recycled
buffer that happened to be of a middling size — not a distinct bug. One fix covers both
regimes.

---

## 4. Fixes — writer-side vs reader-side

There is **one writer** (`cell_malloc_hydro_sorts`'s free) and **many readers** of
`hydro.sort` (density / gradient / force / stars / limiter pairs and their subset variants).

- **Writer-side (preferred): defer reclamation to a barrier.** Do NOT `swift_free` the old
  buffer in place; push it onto a per-cell/space "retired sorts" list and free the whole
  list at the next **rebuild** (`space_rebuild` / `cell_free_hydro_sorts`), which is a hard
  barrier where no interaction task is reading. This is **deterministic and provably safe**
  (the free point has no concurrent readers) — it is NOT a tuned time delay — and memory is
  bounded (cleared every rebuild). One change fixes *all* readers at once.
- **Reader-side: lock each reader.** Extend `extra_sort_lock` across each reader's sort read
  (as already done for stars). Provably correct per reader, but: (a) **completeness burden**
  — every reader must be found and locked or it still crashes; (b) the **limiter reads two
  cells** (`ci` and `cj`), so it needs both `extra_sort_lock`s held, which requires
  **address-ordered acquisition** to avoid deadlock, plus extra contention.

Recommendation: writer-side deferred-to-rebuild reclamation as the unifying real fix;
reader-side locking only if an upstream reviewer prefers a surgical limiter-only change.

---

## 5. Fix and validation (writer-side deferred reclamation) — CONFIRMED

The writer-side fix was implemented properly (not the leak-forever test form): retire the
old buffer to a list and free it at the next rebuild barrier. Patch
`0003-writer-side-defer-hydro-sort-free.patch`:

- `cell.h` — `cell_malloc_hydro_sorts` calls `cell_retire_hydro_sort(c->hydro.sort)` instead
  of `swift_free` when growing; declares the two helpers.
- `cell.c` — a mutex-protected retired-buffer list; `cell_retire_hydro_sort()` (append,
  thread-safe, called under each cell's `extra_sort_lock`) and
  `cell_free_retired_hydro_sorts()` (batch free).
- `space_rebuild.c` — calls `cell_free_retired_hydro_sorts()` at the top of `space_rebuild`,
  which runs on the main thread after `engine_launch`'s `wait_barrier` (no task is reading
  sorts) → the free point provably has no concurrent readers.

The `swift_free("hydro.sort", ...)` in `cell_free_hydro_sorts` (normal rebuild free) is left
untouched; no double-free (retired = old buffers, that path frees the current one).

**Validation run `output-468708.txt` (EAGLE-XL, 4-node, guards left in place):**

| metric | before (unpatched) | after (writer-side fix) |
|--------|--------------------|-------------------------|
| `Caught signal` | crash | **0** |
| `LIMITER-BAD` wild (`>1000`, incl. 5398/5450/1.9e9) | present from ~12 749 s | **0** |
| `LIMITER-BAD` moderate (`count=92`, ~228) | present | **0** |
| `STARDENS-BAD` | (n/a) | **0** |
| progress | crashed at ~12 749 s | **step 191 515, z≈2.07 (a≈0.326), ~64 048 s (~17.8 h), still running** |

All `LIMITER-BAD` (both regimes) and all `STARDENS-BAD` are zero, no crash, ran ~5× past the
original crash point. This confirms the free+swap use-after-free is the cause of the limiter
crash and that deferring the free to the rebuild barrier removes it — with the same fix also
covering the stars reader. The moderate burst vanishing shows it was the same free-race.

Cost: negligible. `cell_retire_hydro_sort` is on the cold "grow a new sort direction" path
(rare; skipped once a cell's directions are all allocated), and the retired-buffer memory is
bounded and cleared every rebuild.

---

## 6. Environment

- Build: production (no `--enable-debugging-checks`), `--with-subgrid=EAGLE-XL
  --with-hydro=sphenix --with-kernel=wendland-C2 --enable-mpi`. Restart-dump compatibility
  requires matching the crashed run's struct-affecting flags.
- Cluster: grammar (Mellanox IB, OpenHPC gnu14/openmpi5); see `[[grammar-cluster-swift-build]]`
  / `[[grammar-cluster-swift-launch]]`.
- Guard inserted by `insert_limiter_guard.py`; log `output-457813.txt`.
- The build also carries the stars `extra_sort_lock` fix + `STARDENS-BAD` guard (harmless to
  the limiter path).
- `LEAKPOS`/`LEAKSYNC` (commit `bb1c6216a`) are under `#ifdef SWIFT_DEBUG_CHECKS` → silent in
  this production build; the position/sync-leak angle is unrelated to this use-after-free.

---

## 7. Next steps

- [x] Writer-side deferred reclamation implemented and validated (§5) — free-race resolved
  for both limiter and stars; moderate burst was the same race.
- [ ] Remove the temporary `LIMITER-BAD` / `STARDENS-BAD` guards for a clean production build.
- [ ] Prepare the upstream SWIFTSIM MR. See `UPSTREAM_MR_SORT_UAF.md` for the packaging
  decision (writer-side fix alone vs. writer-side + keep the stars `extra_sort_lock` read for
  the residual narrow grow half-update window).
