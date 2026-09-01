#!/usr/bin/env python3
"""Generate guarded copies of stars.h and limiter.h in /tmp/guardgen for diffing.
Reproduces the two debug guards exactly as deployed:
  - STARDENS-BAD: inserted AFTER  `struct part *restrict pj = &parts_j[sort_j[pjd].i];`
  - LIMITER-BAD : inserted BEFORE  the 4 sort-index derefs in the limiter DOPAIR1
Both are plain message()+continue, no #ifdef, so they work in a production build."""

D = "/tmp/guardgen"

# ---- STARDENS-BAD (stars.h), inserted AFTER the pj line ----
STARS_KEY = "struct part *restrict pj = &parts_j[sort_j[pjd].i];"
STARS_GUARD = [
    "if (cj->hydro.parts == NULL || sort_j[pjd].i < 0 || sort_j[pjd].i >= count_j) {",
    '  message("STARDENS-BAD id=%lld pjd=%d count_j=%d sort_i=%d sid=%d sorted=%x sort_alloc=%x dx_old=%e dmin=%e", (long long)spi->id, pjd, count_j, (int)sort_j[pjd].i, sid, (unsigned int)cj->hydro.sorted, (unsigned int)cj->hydro.sort_allocated, (double)cj->hydro.dx_max_sort_old, (double)cj->dmin);',
    "  continue;",
    "}",
]


def gen_stars(src, dst):
    out, n = [], 0
    for line in open(src):
        out.append(line)
        if STARS_KEY in line:
            ind = line[:len(line) - len(line.lstrip())]
            n += 1
            for g in STARS_GUARD:
                out.append(ind + g + "\n")
    open(dst, "w").write("".join(out))
    print("stars guard sites:", n, "(expect 2)")


# ---- LIMITER-BAD (limiter.h), inserted BEFORE the deref line ----
LIM_TPL = [
    "if (@parts@ == NULL || @srt@[@idx@].i < 0 ||",
    "    @srt@[@idx@].i >= @cnt@) {",
    '  message("LIMITER-BAD side=@s@ @idx@=%d @cnt@=%d "',
    '    "sort_idx=%d sid=%d ci_foreign=%d cj_foreign=%d "',
    '    "sorted=%x sort_alloc=%x dx_old=%e dmin=%e",',
    "    @idx@, @cnt@, (int)@srt@[@idx@].i, sid,",
    "    (ci->nodeID != e->nodeID),",
    "    (cj->nodeID != e->nodeID),",
    "    (unsigned int)@cell@->hydro.sorted,",
    "    (unsigned int)@cell@->hydro.sort_allocated,",
    "    (double)@cell@->hydro.dx_max_sort_old,",
    "    (double)@cell@->dmin);",
    "  continue;",
    "}",
]


def lim_build(s, idx, cnt, srt, parts, cell):
    r = []
    for l in LIM_TPL:
        for a, b in (("@s@", s), ("@idx@", idx), ("@cnt@", cnt),
                     ("@srt@", srt), ("@parts@", parts), ("@cell@", cell)):
            l = l.replace(a, b)
        r.append(l)
    return r


LIM_GI = lim_build("i", "pid", "count_i", "sort_i", "parts_i", "ci")
LIM_GJ = lim_build("j", "pjd", "count_j", "sort_j", "parts_j", "cj")


def gen_limiter(src, dst):
    out, n = [], 0
    for line in open(src):
        g = None
        if "&parts_i[sort_i[pid].i]" in line:
            g = LIM_GI
        elif "&parts_j[sort_j[pjd].i]" in line:
            g = LIM_GJ
        if g:
            ind = line[:len(line) - len(line.lstrip())]
            for gl in g:
                out.append(ind + gl + "\n")
            n += 1
        out.append(line)
    open(dst, "w").write("".join(out))
    print("limiter guard sites:", n, "(expect 4)")


gen_stars(D + "/stars.pristine", D + "/stars.guarded")
gen_limiter(D + "/limiter.pristine", D + "/limiter.guarded")
