#!/usr/bin/env python3
"""Whole-library control, step 1: reads that could hold an insertion shape.

usage: 47_shape_reads.py ALL_READS.paf OUT.names [MIN_READ_GAP=300] [MAX_REF_GAP=8000]

ALL_READS.paf is a chain-level alignment (minimap2 without -c) of every read of a library, or of
a random sample of it, to the assembly plus the organelle genomes. A read is kept when two of its
alignments lie on the same sequence, on the same strand and in the order they have on the read, no
more than MAX_REF_GAP apart on the sequence but with at least MIN_READ_GAP bases of the read
between them, or when one alignment holds at least MIN_READ_GAP more bases of the read than of the
sequence. Whatever the extra DNA is (plastid, mitochondrial, slug DNA from elsewhere, or nothing
recognisable), the read goes on to the base-level alignment of 33_pass2.sh and the shapes of
46_candidates.py, so that every kind of insert is found by the same steps.
"""
import sys, collections
paf, out = sys.argv[1], sys.argv[2]
MINQ = int(sys.argv[3]) if len(sys.argv) > 3 else 300
MAXR = int(sys.argv[4]) if len(sys.argv) > 4 else 8000
def flush(name, segs, fo):
    if name is None: return 0
    segs.sort()
    for qs, qe, strand, t, ts, te in segs:
        if (qe - qs) - (te - ts) >= MINQ: fo.write(name + "\n"); return 1
    for i in range(len(segs)):
        qs1, qe1, s1, t1, ts1, te1 = segs[i]
        for j in range(i + 1, len(segs)):
            qs2, qe2, s2, t2, ts2, te2 = segs[j]
            if t1 != t2 or s1 != s2 or qs2 - qe1 < MINQ: continue
            g = (ts2 - te1) if s1 == "+" else (ts1 - te2)
            if -300 <= g <= MAXR: fo.write(name + "\n"); return 1
    return 0
by = collections.defaultdict(list); n = 0
for line in open(paf):
    c = line.split("\t", 9)
    by[c[0]].append((int(c[2]), int(c[3]), c[4], c[5], int(c[7]), int(c[8])))
with open(out, "w") as fo:
    for name, segs in by.items(): n += flush(name, segs, fo)
print(f"{len(by)} reads with an alignment; {n} could hold an insertion shape")
