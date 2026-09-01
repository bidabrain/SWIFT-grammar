# Upstream MR: use-after-free on `cell->hydro.sort`

Packaging notes for submitting the fix to SWIFTSIM (gitlab.cosma.dur.ac.uk/swift/swiftsim).

## One-paragraph summary (for the MR description)

`cell_malloc_hydro_sorts()` (`src/cell.h`) grows a cell's hydro sort array by allocating a
new buffer, copying the existing directions, and **`swift_free`-ing the old buffer in place**.
Several readers walk `cell->hydro.sort` **without holding `extra_sort_lock` and without a
scheduler dependency on the sort task**: the time-step limiter pair
(`runner_dopair1_branch_limiter`, no `hydro.sorts`→limiter edge) and the stars ghost-driven
subset density (`DOPAIR1_SUBSET_BRANCH_STARS`, which unlocks before reading). When a runner
thread grows a cell's sort array while another thread is mid-walk, the reader dereferences a
freed-then-recycled buffer, reads an out-of-range `sort[].i`, and segfaults (observed indices
~5400 and ~1.9e9). Fix: retire the grown-over buffer and free it at the next rebuild barrier
instead of in place.

## Evidence (already gathered)

- Stars: `STARS_SORT_RACE_ISSUE.md` — reader-side `extra_sort_lock` fix, A/B validated
  (FLAMINGO, `output-457134.txt` bug → `output-457155.txt` clean, ~17.8 h, 0 occurrences).
- Limiter: `LIMITER_MPI_SORT_CRASH_ISSUE.md` — same free-race, all bad cells local (MPI/
  foreign hypothesis disproved); writer-side fix validated (EAGLE-XL 4-node,
  `output-468708.txt`, all `LIMITER-BAD`+`STARDENS-BAD` = 0, z≈2.07 / ~17.8 h, no crash).

## Two mechanisms, two candidate fixes

1. **Use-after-free (primary, causes the crashes).** Reader walks a buffer the writer
   `swift_free`d. Fixed writer-side by deferring the free.
2. **Grow half-update (secondary, narrow).** `cell_malloc_hydro_sorts` updates
   `sort_allocated |= (1<<j)` and later sets `c->hydro.sort = new_array` — not atomic w.r.t. a
   lock-free reader that computes its block offset from `sort_allocated` and then indexes
   `sort`. A reader landing between the two stores can pair a new offset with the old pointer.
   Nanosecond-wide; usually reads wrong data rather than faulting. Only the reader-side lock
   (holding `extra_sort_lock` across the read, as the stars fix does) closes this.

## Packaging options

- **Option A — writer-side only (`0003`).** Minimal, covers every reader's use-after-free
  with one change. Leaves mechanism (2) open (narrow, non-fatal). Recommended as the core MR.
- **Option B — writer-side + keep stars reader lock (`0002` + `0003`).** Also closes (2) for
  the stars reader. `0002` is tiny and already validated. Slight extra lock hold on the stars
  sorted-subset path only.
- **Option C — writer-side + make grow atomic for all readers.** Fully closes (2) everywhere
  (e.g. seqlock/version on the (pointer, sort_allocated) pair, or publish the new pointer
  before the allocated bits with matching acquire/release). Largest change; probably more
  than reviewers want unless (2) is shown to matter.

**Recommendation:** submit **Option A** as the fix, and mention (2) + the `0002` reader-lock
in the MR description as a follow-up the maintainers can take or leave. Let them decide
whether the narrow half-update window justifies B/C.

## Patch contents

### `0003-writer-side-defer-hydro-sort-free.patch` (the fix)
- `src/cell.h`: in `cell_malloc_hydro_sorts` grow branch, `swift_free(...)` →
  `cell_retire_hydro_sort(c->hydro.sort)`; declare `cell_retire_hydro_sort` /
  `cell_free_retired_hydro_sorts`.
- `src/cell.c`: mutex-protected retired-buffer list + the two functions.
- `src/space_rebuild.c`: `cell_free_retired_hydro_sorts()` at the top of `space_rebuild`.

### `0002-stars-only-hold-extra_sort_lock.patch` (optional, Option B)
- `src/runner_doiact_functions_stars.h`: hold `extra_sort_lock` across the sorted-subset read
  in `DOPAIR1_SUBSET_BRANCH_STARS` (naive branch unlocks early).

## Full patch inventory (repo root)

| patch | purpose | ship upstream? |
|-------|---------|----------------|
| `0003-writer-side-defer-hydro-sort-free.patch` | **core fix** — retire grown-over `hydro.sort`, free at rebuild | **yes (Option A)** |
| `0002-stars-only-hold-extra_sort_lock.patch` | reader-side lock for the residual grow half-update (stars) | optional (Option B) |
| `MR-sort-uaf-combined.patch` | 0002 + 0003 combined | if shipping Option B |
| `0004-add-debug-guards-STARDENS-LIMITER.patch` | ADD both debug guards (`STARDENS-BAD` ×2, `LIMITER-BAD` ×4), plain `message()`+`continue`, no `#ifdef` | **no** (debug only) |
| `0005-remove-debug-guards-STARDENS-LIMITER.patch` | REMOVE both guards (round-trip verified with 0004) | **no** (cleanup helper) |

Guard patches apply with `patch -p1` (works on the non-git cluster tree); `0004` baselines on
pristine `stars.h`/`limiter.h` and does not overlap 0002/0003, so order is free. Round-trip
verified: apply 0004 → STARDENS=2/LIMITER=4; apply 0005 → 0/0.

## Before submitting — cleanup checklist

- [ ] Remove temporary guards: apply `0005-remove-debug-guards-STARDENS-LIMITER.patch`
  (`STARDENS-BAD` in stars.h, `LIMITER-BAD` in limiter.h), and revert the `LEAKPOS`/`LEAKSYNC`
  diagnostics (commit `bb1c6216a`). Re-add later for debugging with `0004-...patch`.
- [ ] Rebase onto current SWIFTSIM `master` and re-run `git apply --check`.
- [ ] One clean commit (or two: fix + optional reader-lock), no debug artefacts.
- [ ] Consider a note to maintainers: the same `cell_malloc_hydro_sorts` free pattern exists
  for **stars sorts** (`cell_malloc_stars_sorts`) — check whether stars-sort readers have the
  same exposure and apply the same retire pattern if so.
- [ ] Reproducer/how-validated paragraph from the two issue docs.

## Open follow-ups to mention in the MR

- Stars-sort buffer (`c->stars.sort`) uses the same grow/free pattern
  (`cell_malloc_stars_sorts`) — audit its readers for the same race.
- Whether to add a `hydro.sorts → timestep_limiter` scheduler dependency as defence-in-depth
  (would also remove the limiter's reliance on sorts persisting untouched).
