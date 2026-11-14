#!/usr/bin/env python3
"""Optional thinning of a pass-1 read set before pass 2 (long reads).

usage: 31b_slim_rdna.py RUN.p1.paf RUN.hits.fa OUT.fa [sample=500]
Libraries rich in ribosomal DNA can leave tens of thousands of reads whose only matches are to
ribosomal DNA. Such a read is kept when one of its matches is to an algal record over >= 1 kb at
<= 2% estimated divergence (long enough to include the variable spacers, so specific to an alga);
of the rest, a random sample is kept. Reads with any other kind of match are all kept.
"""
import sys, collections, random
paf, fa, out = sys.argv[1:4]; sample = int(sys.argv[4]) if len(sys.argv) > 4 else 500
random.seed(7)
only = {}; specific = set()
for line in open(paf):
    c = line.rstrip("\n").split("\t"); cls = c[5].split("|")[0]
    r = c[0]; only[r] = only.get(r, True) and cls in ("algrDNA", "slugrDNA")
    if cls == "algrDNA" and int(c[10]) >= 1000:
        dv = next((float(x[5:]) for x in c[12:] if x.startswith("dv:f:")), 1.0)
        if dv <= 0.02: specific.add(r)
cand = [r for r, o in only.items() if o and r not in specific]
drop = set(cand) - set(random.sample(cand, min(sample, len(cand))))
n = k = 0; keep = True
with open(out, "w") as fo:
    for line in open(fa):
        if line[0] == ">":
            n += 1; keep = line[1:].split()[0] not in drop; k += keep
        if keep: fo.write(line)
print(f"{n} reads in; {len(cand)} match only ribosomal DNA without a specific algal match, {len(drop)} of them dropped; {len(specific)} have a specific algal rDNA match; {k} reads out")
