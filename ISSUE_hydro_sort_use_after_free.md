Possible use-after-free on cell->hydro.sort (crashes in the time-step limiter and in the stars density ghost)

Hi, and thank you for SWIFT!

I have been hitting reproducible segfaults in two different runs and, after some instrumentation, I suspect they may share a common cause: a lock-free read of a cell's hydro sort array (cell->hydro.sort) that races with the array being freed/reallocated in cell_malloc_hydro_sorts(). I am not certain this is the real root cause - I may well be misreading the code - but I wanted to report it with the data I collected and the two small patches I tried, in case it is useful or in case someone can tell me what I got wrong.

I have attached the two patches; below I just describe what they do.


ENVIRONMENT

SWIFT version/commit: c524857e5
Build: production (no --enable-debugging-checks), --enable-mpi --enable-fof, --with-hydro=sphenix, --with-kernel=wendland-C2.
Compiler/stack: GCC 14.2, OpenMPI 5.0.7, HDF5 1.14.6, UCX 1.18 (mt), on an InfiniBand cluster.

Two models:
  - EAGLE (--eagle, --with-subgrid=EAGLE-XL): 4 nodes, 1 MPI rank/node, 32 threads/rank.
  - FLAMINGO (--with-subgrid=FLAMINGO): single node, 32 threads.

Launch (EAGLE): mpirun -np 4 --map-by ppr:1:node:PE=32 ... ./swift_mpi --threads=32 --cosmology --eagle --verbose=1 eagle.yml, with OMP_NUM_THREADS=32, UCX_TLS=rc,sm,self, UCX_RNDV_THRESH=inf.


SYMPTOM 1 - EAGLE, crash in the time-step limiter (multi-node)

Shortly before the crash there are some limiter warnings, then a segfault inside runner_dopair1_branch_limiter:

  [0002] [09624.9] runner_do_limiter: WARNING: Not limiting particle with id 31141401 because it needs to be synced.
  ...
  [0000] [12231.3] runner_do_limiter: WARNING: Not limiting particle with id 33104139 because it needs to be synced.
  [grammar005:...] Caught signal 11 (Segmentation fault: address not mapped to object at address 0x1550b88dac52)
  ==== backtrace ====
   3  runner_dopair1_branch_limiter()
   4  runner_dosub_pair1_limiter.isra.0()
   6  runner_main()
   7  start_thread()


SYMPTOM 2 - FLAMINGO, crash in the stars density (single node)

Here a few "h ~ h_max" warnings precede a segfault reached through the stars ghost (runner_do_stars_ghost -> runner_dosub_pair_subset_stars_density):

  [0000] [34406.7] hydro_prepare_gradient: WARNING: h ~ h_max for particle with ID 24261446 (h: 0.499973)
  [grammar073:...] Caught signal 11 (Segmentation fault ...)
  ==== backtrace ====
   (frames in the stars density interaction, entered from the stars ghost)

In both cases the crash is a bad memory access while walking a cell's sorted particle list and fetching a particle by its stored sort index, e.g. parts_j[sort_j[pjd].i].


WHAT I TRIED (INSTRUMENTATION)

Because these are production builds I added two small non-fatal guards (plain message() + continue, not under #ifdef SWIFT_DEBUG_CHECKS, so they print and skip instead of crashing):

  - STARDENS-BAD in runner_doiact_functions_stars.h, right after pj = &parts_j[sort_j[pjd].i], checking sort_j[pjd].i against count_j.
  - LIMITER-BAD in runner_doiact_functions_limiter.h, before the four parts[sort[...].i] dereferences in the limiter DOPAIR1, checking the index against the cell's count and also printing whether each cell is foreign (nodeID != e->nodeID).

Both dump sort_idx, count, sid, sorted, sort_allocated, dx_max_sort_old, dmin.

FLAMINGO (STARDENS-BAD): the stored sort index is wildly out of range, e.g. sort_idx ~ 5263 / 5392 / 5393 while count_j is only a few hundred, and the sorted flag for that direction is set (so the branch believes the cell is sorted).

EAGLE (LIMITER-BAD): 132 occurrences, and notably:
  - cj_foreign = 0 for all of them - the offending cells are local, not foreign (so this does not look like an MPI/foreign-buffer problem, which is what I first suspected).
  - spread across all 4 ranks and many steps (i.e. systematic, not a one-off).
  - sort_idx came in two flavours: some "wild" values (5398, 5450, even 1936877413) and some "moderate" ones (~228 while count_j = 92), with sorted == sort_allocated and dx_max_sort_old = 0.

The wild EAGLE values look very similar to the FLAMINGO ones (~5400), which made me wonder whether both crashes are reading a sort buffer that has been freed and recycled.


MY (TENTATIVE) INTERPRETATION

Looking at cell_malloc_hydro_sorts() in cell.h, when a cell's sort array is grown to add a newly-requested direction it allocates a new buffer, copies the existing directions, and then frees the old buffer in place:

  /* Swap the pointers */
  swift_free("hydro.sort", c->hydro.sort);
  c->hydro.sort = new_array;

Meanwhile it looks (to me) like a couple of readers walk cell->hydro.sort without holding extra_sort_lock and without a scheduler dependency on the sort task:

  - the time-step limiter pair interaction - I could not find a hydro.sorts -> limiter dependency edge, so it seems to rely on the sorts staying valid; and
  - the stars ghost subset density (DOPAIR1_SUBSET_BRANCH_STARS), which takes extra_sort_lock only to check the sorted flag and then unlocks before reading the array.

If that reading is correct, then a runner thread growing a cell's sort array (which happens from the ghost phase) could free the buffer while another thread is still walking it, giving an out-of-range sort[].i and the segfaults above. This would also explain the near-constant "wild" indices (a freed block reused by a larger cell) and why disabling particle splitting only reduced the frequency (fewer new directions requested) rather than removing the crash.

I want to stress this is a hypothesis - I may be missing a dependency or an invariant that normally prevents this.


TWO THINGS I TRIED, AND THE RESULTS

I attached two patches:

1. 0002-stars-only-hold-extra_sort_lock.patch - reader side: in DOPAIR1_SUBSET_BRANCH_STARS, keep extra_sort_lock held across the sorted-subset read instead of unlocking right after the flag check (the naive branch, which does not touch the sort array, still unlocks early).

2. 0003-writer-side-defer-hydro-sort-free.patch - writer side: instead of swift_free-ing the old sort buffer in place in cell_malloc_hydro_sorts(), retire it to a small list and free that list at the next rebuild (from space_rebuild()), where - as I understand it - no task is reading sort arrays. The intent is that any in-flight reader keeps a valid pointer, and memory is bounded (the list is cleared every rebuild).

Results (guards left in place so I could tell whether the bad reads still happened):

  - FLAMINGO + patch 0002: previously STARDENS-BAD fired reliably at ~484 s and ~639 s and then crashed; with the patch it ran to ~64,000 s (~17.8 h) with zero STARDENS-BAD and no crash.

  - EAGLE (4 nodes) + patch 0003: previously it crashed at ~12,700 s; with the patch it ran to redshift z ~ 2.07 (~64,000 s, ~step 191,500) with zero LIMITER-BAD (both the wild and the moderate values disappeared) and zero STARDENS-BAD, and no crash.

So in my testing the writer-side change (0003) alone removed both the limiter and the stars bad reads, which is what made me think the two crashes might share the same underlying free.


QUESTIONS

  - Is my reading correct that the limiter pair and the stars ghost subset can read cell->hydro.sort concurrently with cell_malloc_hydro_sorts() freeing it? If not, what invariant am I missing that should prevent it?
  - If it is a real race, would you prefer the writer-side approach (defer the free) or adding the missing reader-side protection/dependency? I am happy to adapt the patches to whatever fits the codebase best.
  - I noticed cell_malloc_stars_sorts() uses the same grow-and-free pattern for stars.sort; I did not touch it, but it may be worth checking whether its readers have the same exposure.

Thanks very much for taking a look, and apologies if I have misunderstood something.
