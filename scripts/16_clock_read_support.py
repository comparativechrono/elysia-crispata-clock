#!/usr/bin/env python3
"""HiFi read support for clock-gene exons, from the pass-1 bait screen.

usage: 16_clock_read_support.py catalogue.tsv clock_exon_baits.tsv out.json pass1/*.bait.paf
For each gene, the exon segments (exon +-150 bp, taken from the Purdue-assembly gene models)
are checked for xgElyCris1 HiFi reads aligning over >= 150 bp at <= 5 % divergence.
"""
import sys, json, collections, statistics

cat, baits, out = sys.argv[1:4]
pafs = sys.argv[4:]
seg_of = collections.defaultdict(list)          # accession -> [segment names]
masked = {}
for line in open(baits):
    name, q, chrom, s, e, strand = line.rstrip("\n").split("\t")
    seg_of[q.split("_")[1]].append(name)
hits = collections.Counter(); reads = collections.defaultdict(set)
for p in pafs:
    for line in open(p):
        c = line.split("\t")
        if not c[5].startswith("clk|"): continue
        if int(c[10]) < 150: continue
        dv = 0.0
        for t in c[12:]:
            if t.startswith("dv:f:"): dv = float(t[5:]); break
        if dv > 0.05: continue
        hits[c[5]] += 1; reads[c[5]].add(c[0])
res = {}
for line in open(cat):
    if line.startswith("#") or line.startswith("symbol\t"): continue
    f = line.rstrip("\n").split("\t")
    accs = [a for a in f[5].split(",") if not a.startswith("ENSUFQG")]
    segs = sorted({s for a in accs for s in seg_of.get(a, [])})
    if not segs: continue
    depth = [len(reads[s]) for s in segs]
    allreads = set().union(*[reads[s] for s in segs]) if segs else set()
    res[f[0]] = dict(segments=len(segs), covered=sum(1 for d in depth if d >= 3), median_reads=statistics.median(depth), min_reads=min(depth),
                     total_reads=len(allreads), per_segment=[dict(segment=s.split("|")[2], reads=d) for s, d in zip(segs, depth)])
json.dump(res, open(out, "w"), indent=1)
alld = [d["median_reads"] for d in res.values()]
print(f"{len(res)} genes; median of per-gene median exon depth = {statistics.median(alld)}")
for k, d in res.items():
    print(f"{k:13s} exon segments {d['covered']:2d}/{d['segments']:2d} covered (>=3 reads); median {d['median_reads']:5.1f} reads, min {d['min_reads']}")
