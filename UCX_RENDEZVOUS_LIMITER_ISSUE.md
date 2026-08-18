# SWIFT hang / crash under MPI: limiter buffers vs. UCX rendezvous

**Status:** analysis / not yet root-fixed
**Affected code:** current `master` (byte-identical to upstream `SWIFTSIM/SWIFT` @ `9f4ca58a6`); reproduced on stock SWIFT as well as the local `fix-limiter-sync-mpi` branch.
**Trigger:** InfiniBand fabric via UCX (OpenMPI 4.1.5 `pml=ucx`) with rendezvous protocol enabled (the default for large messages).

---

## 1. Symptoms

Two distinct failure modes, same underlying cause:

### A. With TCP BTL (`OMPI_MCA_pml=ob1`, `OMPI_MCA_btl=self,tcp`)
- Run proceeds for a long time (hours), then **hangs**.
- `top` shows a single `swift_mpi` process at ~100 % on **one core**; the other 31 OpenMPI threads sleep; `load average ≈ 1.0`.
- The busy thread is spinning in `opal_progress` (MPI polling) waiting on a message that never completes — a classic TCP-BTL + `MPI_THREAD_MULTIPLE` deadlock.

### B. With UCX over InfiniBand (`pml=ucx`)
Two crash signatures observed in the same run:

**Node grammar001 — IB remote access error (sender-side buffer torn down mid-transfer):**
```
Remote access error on mlx5_0:1/IB (synd 0x13 vend 0x88)
RC QP ... RDMA_READ  [rva 0x14b9... rkey 0x2c2100] [va 0x154f... len 240 lkey 0x27f600]
  uct_ib_mlx5_completion_with_err
  ucp_worker_progress
  opal_progress
  MPI_Test
  task_lock            <-- SWIFT scheduler polling the MPI request
  queue_gettask
  runner_main
```

**Node grammar020 — segfault reading foreign particle data in the limiter:**
```
Caught signal 11 (Segmentation fault: address not mapped)
  runner_dopair1_branch_limiter
  runner_dosub_pair1_limiter
  runner_main
```

The RDMA_READ error means a **receiver was still reading a remote sender's buffer whose memory registration had already been invalidated** (freed / unmapped). The segfault means a **limiter pair task dereferenced a foreign `parts` array that was no longer mapped**. Both are the same class of use-after-free / stale-registration, exposed only under RDMA rendezvous.

---

## 2. Environment

- Cluster: `grammar*`, Mellanox InfiniBand (`mlx5_0:1`, port Active / LinkUp).
- OpenMPI `4.1.5` (OpenHPC), built `--with-ucx=/opt/ohpc/.../ucx-ohpc/1.14.0`, `--without-verbs` (verbs reached through UCX). `ompi_info` shows `pml: ucx` and `osc: ucx`.
- System UCX `1.18.0` (`--enable-mt`, `--with-mlx5`, `--with-verbs`, rcache/ucm hooks enabled).
- Job layout: `-N 4`, `--ntasks-per-node=1`, `--cpus-per-task=32` → **only 4 MPI ranks, 32 OpenMP threads each**. Few ranks ⇒ large per-rank foreign data and one heavily-loaded progress path.
- The original `sbatch` did `module unload ucx` and forced `pml=ob1 / btl=self,tcp`. That unload removes the UCX that OpenMPI 4.1.5 was built against, which is *why* the run fell back to TCP in the first place.

---

## 3. How SWIFT manages MPI transfer buffers

Three mechanisms, all verified in the source:

1. **Zero-copy sends from the live particle array.** For `xv`/`rho` the send buffer *is* the live cell array:
   `src/scheduler.c:1098` → `buff = t->ci->hydro.parts;` then `MPI_Isend(buff, count, part_mpi_type, ...)`. No copy is made on the sender side.

2. **Completion is detected by lazy `MPI_Test` polling from worker threads.** A send/recv task is only considered "runnable/complete" once its request tests done:
   `src/task.c:734` → `MPI_Test(&t->req, &res, &stat); ... return res;`
   (This is exactly the `task_lock → MPI_Test → opal_progress` frame in the grammar001 backtrace.) MPI requests are *posted* at enqueue time in `scheduler.c` (`MPI_Irecv` @ `1042`, `MPI_Isend` @ `1153`).

3. **Memory is only moved / freed at rebuild.** Foreign arrays are released by `space_free_foreign_parts` (`src/space.c:142`), called from `engine_rebuild` (`src/engine.c:1370`, `1479`). Rebuild is always preceded by an `engine_launch` barrier that drains the task queue.

### The limiter's four-stage exchange (pack → send / recv → unpack)

The time-step limiter exchanges per-particle timebins through a **heap buffer shared between two tasks** via `task_get_unique_dependent`:

- **Send side:** a `pack` task allocates + fills `t->buff`, then hands it to the send task:
  `src/runner_main.c` (pack case) → `runner_do_pack_limiter(r, ci, &t->buff, 1); task_get_unique_dependent(t)->buff = t->buff;`
  The send task frees it after completion: `runner_main.c:444` → `free(t->buff)` (limiter subtype).
- **Recv side:** the recv task `posix_memalign`s the buffer and shares it with its unpack dependent:
  `src/scheduler.c:989-996` → `t->buff = buff; task_get_unique_dependent(t)->buff = buff;`
  The recv itself does nothing (`runner_main.c`: *"Nothing to do here. Unpacking done in a separate task"*); the **unpack task frees it after use**:
  `runner_do_unpack_limiter` → `cell_unpack_timebin(c, buffer); free(buffer);`
- **Dependency wiring:** `engine_maketasks.c:724-730` creates `recv_limiter → unpack_limiter` (`scheduler_addunlock`), and `engine_maketasks.c:254-259` creates `pack_limiter → send_limiter`.
- `task_get_unique_dependent` (`src/task.c:1101`) asserts `nr_unlock_tasks == 1` under `SWIFT_DEBUG_CHECKS`, so the sharing is strictly 1:1.

**Conclusion of the ordering audit:** every transfer buffer is freed only by a task that runs *after* its MPI request has tested complete (send frees after `MPI_Test`; unpack frees after the `recv→unpack` unlock; foreign parts freed only at rebuild, behind the launch barrier). **There is no plain "free-before-MPI-complete" bug in the task graph** — which is consistent with the failure being intermittent rather than deterministic.

---

## 4. Why eager works and rendezvous crashes

The correctness of SWIFT's model implicitly assumes **eager / copy** transport semantics.

- **Eager (and TCP BTL):** `MPI_Isend` copies the payload out of the user buffer into internal/bounce buffers immediately. The user buffer is decoupled within microseconds; the user buffer is **never registered for remote RDMA read**. Any later modification, free, or reuse of that memory is harmless. → SWIFT's buffer handling is safe. (TCP still deadlocks under `MPI_THREAD_MULTIPLE`, giving failure mode A.)

- **Rendezvous (UCX default for large messages):** `MPI_Isend` does **not** copy. It pins/registers the user buffer, gets an `rkey`, sends a small RTS, and the **remote peer RDMA-reads the sender's live buffer directly**. The request completes only after that remote read. The buffer must therefore stay allocated, unmoved, and its registration valid across a long, remote-scheduling-dependent window.

Two structural weak points make this fail:

### Weak point 1 — the limiter/sync **sub-cycle re-uses `engine_launch`**
`src/engine.c:2976` (verbatim comment): *"engine_launch is re-used for the limiter and sync"*. So within one top-level step, after the main force loop, additional `engine_launch` rounds run the limiter and sync sub-graphs. Each round does another **malloc → register → transfer → free** cycle for the limiter timebin buffers, multiplying the registration churn per step. This is precisely the workload UCX's registration cache handles worst.

### Weak point 2 — per-step malloc/free of **registered** buffers vs. UCX rcache + ucm hooks
UCX caches memory registrations keyed by `[address, length]` (rcache) and installs `malloc`/`munmap` hooks (ucm) to invalidate the cache when the application frees memory.

- `free()` of a buffer **below glibc's mmap threshold** returns it to the heap arena **without `munmap`**, so the **ucm hook does not fire** and the **stale registration stays in rcache**.
- The next step's `malloc` of the same size returns the **same virtual address**; SWIFT fills it with new data; UCX reuses the cached registration. Locally this is memory-correct (same physical pages), but it couples the *previous* step's `rkey` to the *current* buffer.
- When rcache is later evicted (memory pressure) or a neighboring large-buffer `munmap` triggers the ucm hook and invalidates a registration **while a remote RDMA_READ using that `rkey` is still in flight**, the remote read lands on memory whose registration no longer maps → **`Remote access error synd 0x13`** (grammar001). Adjacent reuse/misalignment surfaces as the **limiter segfault** (grammar020).

Under `MPI_THREAD_MULTIPLE`, UCX progress runs on arbitrary worker threads (`task_lock → MPI_Test`) while other threads malloc/free, so the hooks and rcache are stressed **concurrently** — widening the race.

**Bottom line:** the bug is *not* a single line that frees a foreign buffer too early. It is SWIFT's pattern of **frequently allocating and freeing already-registered communication buffers every step** (amplified by the limiter/sync sub-cycle) colliding with UCX's registration-cache invalidation race. Eager transport hides it entirely because it never exposes the user buffer to remote RDMA.

---

## 5. How to confirm (turn hypothesis into proof)

1. **SWIFT's own MPI-use reporting.** `scheduler.c:1050 / 1165` call `mpiuse_log_allocation`, which records the post/complete lifecycle of every request. Reconfigure with `--enable-mpiuse-reports`, rebuild, and inspect the report for the step before the crash: a limiter/gpart send or recv that is **posted but never completed** confirms a dangling request + reused buffer.

2. **A/B on the transport, one variable at a time:**
   - `UCX_RNDV_THRESH=inf` (force eager for all sizes). If the crash disappears ⇒ the rendezvous RDMA path is the cause.
   - `UCX_IB_RCACHE=n` (keep rendezvous, disable the registration cache). If that also fixes it ⇒ the culprit is specifically **stale rcache registrations**, not a SWIFT dependency-ordering bug.

---

## 6. Mitigation and recommendations

**Immediate (production workaround) — restore copy semantics that SWIFT relies on:**
```bash
# do NOT `module unload ucx`; use the UCX OpenMPI was built against
export OMPI_MCA_pml=ucx
export OMPI_MCA_osc=ucx
export UCX_TLS=rc,sm,self
export UCX_NET_DEVICES=mlx5_0:1
export UCX_RNDV_THRESH=inf     # force eager; behaves like TCP's copy, no thread deadlock
# export UCX_IB_RCACHE=n       # add if crashes persist
```
This is strictly better than TCP: it avoids both the rendezvous use-after-free **and** the TCP-BTL thread deadlock.

**Secondary — reduce pressure:** run more MPI ranks (e.g. 1 per NUMA socket: `--ntasks-per-node=2 --cpus-per-task=16`) so each rank's progress path and message sizes are smaller.

**Root fix (upstream):** because the code is byte-identical to upstream SWIFT, this is a general SWIFT-on-UCX-rendezvous problem worth reporting to SWIFTSIM with the mpiuse report attached. The durable fix is to stop churning registered memory every step — allocate the limiter/gpart/… transfer buffers **once from a reusable pool** (stable addresses, registered once) instead of `malloc`/`free` per step — or otherwise make buffer registration lifetime explicit and rendezvous-safe, rather than relying on `UCX_RNDV_THRESH`.

---

## 7. Key source references

| What | Location |
|---|---|
| Zero-copy `xv` send from live array | `src/scheduler.c:1098` |
| `MPI_Irecv` / `MPI_Isend` posted at enqueue | `src/scheduler.c:1042`, `1153` |
| `MPI_Test`-gated task completion | `src/task.c:734` |
| Limiter recv buffer alloc + share with unpack | `src/scheduler.c:989-996` |
| Send/pack buffer free | `src/runner_main.c:444` (send), `runner_do_unpack_limiter` (unpack) |
| `pack→send`, `recv→unpack` dependencies | `src/engine_maketasks.c:254-259`, `724-730` |
| `task_get_unique_dependent` (1:1 assert) | `src/task.c:1101` |
| Foreign parts free | `src/space.c:142`; called from `src/engine.c:1370`, `1479` |
| **limiter/sync re-use engine_launch** | `src/engine.c:2976` |
| MPI-use lifecycle logging hook | `src/scheduler.c:1050`, `1165` (`--enable-mpiuse-reports`) |
