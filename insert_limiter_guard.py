#!/usr/bin/env python3
"""Insert non-fatal LIMITER-BAD bounds guards before the 4 sort-index
dereferences in the time-step limiter DOPAIR1. Plain message()+continue,
NOT under #ifdef SWIFT_DEBUG_CHECKS, so it fires in a production build too.
Keyed on substrings, so it is robust to line-number differences."""

f = 'src/runner_doiact_functions_limiter.h'

TPL = [
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


def build(s, idx, cnt, srt, parts, cell):
    r = []
    for l in TPL:
        l = l.replace("@s@", s)
        l = l.replace("@idx@", idx)
        l = l.replace("@cnt@", cnt)
        l = l.replace("@srt@", srt)
        l = l.replace("@parts@", parts)
        l = l.replace("@cell@", cell)
        r.append(l)
    return r


GI = build("i", "pid", "count_i", "sort_i", "parts_i", "ci")
GJ = build("j", "pjd", "count_j", "sort_j", "parts_j", "cj")

out, n = [], 0
for line in open(f):
    g = None
    if "&parts_i[sort_i[pid].i]" in line:
        g = GI
    elif "&parts_j[sort_j[pjd].i]" in line:
        g = GJ
    if g:
        ind = line[:len(line) - len(line.lstrip())]
        for gl in g:
            out.append(ind + gl + "\n")
        n += 1
    out.append(line)

open(f, "w").write("".join(out))
print("patched", n, "sites expect 4")
