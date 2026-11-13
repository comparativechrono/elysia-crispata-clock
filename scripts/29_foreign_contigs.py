#!/usr/bin/env python3
"""Summarise the checks of script 29 into foreign_contigs.json.
usage: 29_foreign_contigs.py workdir asm_name feature_table.gz alien_hits.tsv other=fasta ...
"""
import sys, os, json, gzip, collections
wd, asm, ft, hits = sys.argv[1:5]
others = [a.split("=", 1)[0] for a in sys.argv[5:]]
R = json.load(open(f"{wd}/protein_screen_raw.json"))
seqs = {}; name = None
for l in open(f"{wd}/{asm}_foreign.fa"):
    if l[0] == ">": name = l[1:].split()[0]; seqs[name] = []
    else: seqs[name].append(l.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
genes = collections.defaultdict(list)
with gzip.open(ft, "rt") as fh:
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if f[0] == "CDS" and f[6] in seqs: genes[f[6]].append(f[16])
tags = {t for v in genes.values() for t in v}
with_hit = set()
for l in open(hits):
    q = l.split("\t", 1)[0].split("|")
    if len(q) > 2 and q[2] in tags: with_hit.add(q[2])
def paf(path, minlen=500):
    d = collections.defaultdict(list)
    if os.path.exists(path):
        for l in open(path):
            f = l.split("\t")
            if len(f) >= 12 and int(f[10]) >= minlen: d[f[0]].append(f)
    return d
aln = {o: paf(f"{wd}/{asm}_foreign_vs_{o}.paf") for o in others}
rd = collections.defaultdict(list)
rp = f"{wd}/reads_vs_{asm}_foreign.paf"
if os.path.exists(rp):
    for l in open(rp):
        f = l.split("\t")
        if len(f) >= 12: rd[f[5]].append((int(f[7]), int(f[8])))
out = []
for c, s in sorted(seqs.items()):
    loci = [l for l in R["loci"] if l["assembly"] == asm and l["contig"] == c]
    row = dict(assembly=asm, contig=c, length=len(s), gc=round((s.count("G") + s.count("C")) / max(1, len(s) - s.count("N")), 4),
               foreign_loci=[dict(query=l["query"], start=l["start"], end=l["end"], identity=l["identity"], cls=l["class"], frameshifts=l["frameshifts"]) for l in loci if l["class"] in ("alga-like", "microbe-like")],
               animal_loci=sum(1 for l in loci if l["class"] == "animal gene"),
               annotated_genes=len(genes.get(c, [])), annotated_genes_with_protein_match=sum(1 for t in genes.get(c, []) if t in with_hit))
    for o in others:
        row[f"aligned_to_{o}"] = [dict(start=int(f[2]), end=int(f[3]), target=f[5], tstart=int(f[7]), tend=int(f[8]), identity=round(int(f[9]) / int(f[10]), 3)) for f in aln[o].get(c, [])]
    iv = sorted(rd.get(c, [])); covered = 0; last = -1
    for a, b in iv:
        a = max(a, last)
        if b > a: covered += b - a; last = b
    row["read_alignments"] = len(iv); row["bp_covered_by_reads"] = covered
    row["read_cover_span"] = [min(a for a, b in iv), max(b for a, b in iv)] if iv else None
    out.append(row)
    print(f"{c:20s} {len(s):7d} bp GC {100 * row['gc']:.1f}% | foreign loci {len(row['foreign_loci']):2d}, animal loci {row['animal_loci']} | annotated genes {row['annotated_genes']}, with a protein match {row['annotated_genes_with_protein_match']} | "
          + " ".join(f"{o}: {sum(a['end'] - a['start'] for a in row[f'aligned_to_{o}'])} bp" for o in others) + f" | reads: {len(iv)} alignments covering {covered} bp {row['read_cover_span'] or ''}")
json.dump(out, open(f"{wd}/{asm}_foreign_contigs.json", "w"), indent=1)
