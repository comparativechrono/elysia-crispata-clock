#!/usr/bin/env python3
"""E. crispata: reads that hold a short, near-exact piece of organelle DNA away from both ends.

usage: 54_internal_pieces.py PASS1_DIR OUTDIR
PASS1_DIR holds RUN.bait.paf of scripts/10_bait_screen.sh: every read of the library mapped to
the bait panel. A piece of a few hundred bases inside a read does not break the alignment of the
read to the slug sequence, so the second pass of that pipeline (scripts/14, 15) passes over it.
Chosen here: reads with a match to a plastid or mitochondrial sequence of 60 to 5,000 bases at 95%
identity or better (minimap2's estimate) that leaves 200 bases or more of the read on both sides.
Matches to "hot" stretches of the panel are left out: a 300-base stretch of a panel sequence to
which more than six such reads match, or which lies next to one, is a slug repeat that resembles
it (most are transfer RNA genes carried by mobile elements; one such stretch collects 115,000 reads).
Writes OUTDIR/targeted.names and OUTDIR/targeted.p1.paf (the matches of the chosen reads, for
46_candidates.py). The reads themselves are fetched with 52_fetch_reads.py.
"""
import sys, os, glob, collections
p1, out = sys.argv[1:3]
ORG = ("cp", "cpx", "slugmt", "algmt")
sel = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(p1, "*.bait.paf"))):
    for line in open(f):
        c = line.rstrip("\n").split("\t")
        if c[5].split("|")[0] not in ORG: continue
        qlen, qs, qe = int(c[1]), int(c[2]), int(c[3])
        dv = next((float(x[5:]) for x in c[12:] if x.startswith("dv:f:")), 1.0)
        if qs >= 200 and qlen - qe >= 200 and 60 <= qe - qs <= 5000 and dv <= 0.05: sel[c[0]].append((c[5], int(c[7])))
bins = collections.Counter((t, ts // 300) for v in sel.values() for t, ts in v)
keep = sorted(r for r, v in sel.items() if any(bins[(t, ts // 300)] <= 6 and bins[(t, ts // 300 - 1)] <= 6 and bins[(t, ts // 300 + 1)] <= 6 for t, ts in v))
os.makedirs(out, exist_ok=True)
open(os.path.join(out, "targeted.names"), "w").write("\n".join(keep) + "\n")
ks = set(keep); n = 0
with open(os.path.join(out, "targeted.p1.paf"), "w") as fo:
    for f in sorted(glob.glob(os.path.join(p1, "*.bait.paf"))):
        for line in open(f):
            c = line.split("\t", 6)
            if c[0] in ks and c[5].split("|")[0] in ("cp", "slugmt", "algmt"): fo.write(line); n += 1
print(f"{len(sel)} reads with a near-exact organelle piece away from both ends; {len(keep)} outside the hot stretches; {n} matches written")
