#!/usr/bin/env python3
"""Collect homologues for each clock-gene family from the DIAMOND seed search.

usage: 03_family_sets.py seeds.fa seeds_vs_all.tsv allprot.fa allprot.tsv outdir [K]
For every seed the top K hits per species (e-value <= 1e-8) are kept; isoforms
of one gene are collapsed to the longest protein.  One FASTA per family is
written (seeds included) for alignment and tree building.
"""
import sys, collections, os

seeds_fa, hits_tsv, all_fa, all_tsv, outdir = sys.argv[1:6]
K = int(sys.argv[6]) if len(sys.argv) > 6 else 4
os.makedirs(outdir, exist_ok=True)

def read_fa(p):
    d, name = {}, None
    for line in open(p):
        if line[0] == ">":
            name = line[1:].split()[0]; d[name] = []
        else:
            d[name].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}

seeds = read_fa(seeds_fa)
fam_of = {s: s.split("|")[2] for s in seeds}
ann = {}
for line in open(all_tsv):
    f = line.rstrip("\n").split("\t")
    ann[f[0]] = f
hits = collections.defaultdict(lambda: collections.defaultdict(list))  # seed -> species -> hits
for line in open(hits_tsv):
    q, s, pid, alen, qlen, slen, qs, qe, ss, se, ev, bs = line.rstrip().split("\t")
    if float(ev) > 1e-8: continue
    hits[q][s.split("|")[0]].append((float(bs), s))
keep = collections.defaultdict(set)
for q, bysp in hits.items():
    for sp, hh in bysp.items():
        hh.sort(reverse=True)
        seen_gene, n = set(), 0
        for bs, s in hh:
            gene = ann[s][3] or s
            if (sp, gene) in seen_gene: continue
            seen_gene.add((sp, gene)); n += 1
            keep[fam_of[q]].add(s)
            if n >= K: break
allp = read_fa(all_fa)
for fam, ids in keep.items():
    # collapse isoforms: (species, gene) -> longest
    best = {}
    for s in ids:
        sp = s.split("|")[0]; gene = ann[s][3] or s
        key = (sp, gene)
        if key not in best or len(allp[s]) > len(allp[best[key]]): best[key] = s
    # make sure the longest isoform overall of each kept gene is used
    with open(os.path.join(outdir, fam + ".fa"), "w") as fo:
        for s, seq in seeds.items():
            if fam_of[s] == fam: fo.write(f">SEED_{s.split('|')[0]}\n{seq}\n")
        for key, s in sorted(best.items()):
            fo.write(f">{s.replace('|','_')}\n{allp[s]}\n")
    print(fam, len(best), "homologues;", sum(1 for k in best if k[0] == "ELYCR"), "ELYCR,", sum(1 for k in best if k[0] == "ELYCRP"), "ELYCRP")
