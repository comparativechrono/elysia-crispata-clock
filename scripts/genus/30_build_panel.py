#!/usr/bin/env python3
"""Build the reference panel used to screen the reads of every species.

usage: 30_build_panel.py crispata_panel.masked.fa crispata_panel.tsv dasyclad_cp.fa cp_extra.fa extra_rdna.fa elysia_mito.fa outdir

Classes (the part of each sequence name before "|"):
  cp       plastid genomes and plastid gene records: the green-algal plastomes of the E. crispata panel,
           Dasycladales plastid records, Xanthophyceae and Cladophorales plastid records, and a few
           distant plastomes (land plant, diatom, red alga, cyanobacterium) to catch any plastid
  algmt    algal mitochondrial genomes
  algrDNA  nuclear ribosomal DNA of algae (a sign of algal nuclei in the sample)
  slugmt   mitochondrial genomes of Elysia species
  slugrDNA slug ribosomal DNA
The bacterial genomes and clock-gene baits of the E. crispata panel are left out.
"""
import sys, os, re
panel, panel_tsv, dasy, cpx, rdna, mito, out = sys.argv[1:8]
os.makedirs(out, exist_ok=True)

def fasta(p):
    n, d, buf = None, "", []
    for line in open(p):
        if line[0] == ">":
            if n is not None: yield n, d, "".join(buf)
            f = line[1:].rstrip().split(" ", 1); n = f[0]; d = f[1] if len(f) > 1 else ""; buf = []
        else: buf.append(line.strip())
    if n is not None: yield n, d, "".join(buf)

desc = {}
for l in list(open(panel_tsv))[1:]:
    f = l.rstrip("\n").split("\t"); desc[f[0]] = f[4]
seen = set(); rows = []
fo = open(os.path.join(out, "panel.fa"), "w")
def put(cls, acc, d, s, src):
    if acc in seen or len(s) < 200: return
    seen.add(acc); name = f"{cls}|{acc}"
    fo.write(f">{name}\n{s}\n"); rows.append((name, cls, acc, len(s), d, src))
KEEP = {"cp": "cp", "cpx": "cp", "algmt": "algmt", "algrDNA": "algrDNA", "slugmt": "slugmt", "slugrDNA": "slugrDNA"}
for n, d, s in fasta(panel):
    c, acc = n.split("|", 1)
    if c in KEEP: put(KEEP[c], acc, desc.get(n, d), s, "E. crispata panel")
for n, d, s in fasta(dasy): put("cp", n, d, s, "Dasycladales plastid records")
for n, d, s in fasta(cpx):
    if d.startswith("TSA:") or "mitochondri" in d.lower(): continue
    put("cp", n, d, s, "Xanthophyceae and Cladophorales plastid records")
for n, d, s in fasta(rdna):
    if re.search(r"chloroplast|plastid|mitochondri|16S|23S", d, re.I): continue
    put("algrDNA", n, d, s, "algal nuclear rDNA")
for n, d, s in fasta(mito): put("slugmt", n, d, s, "Elysia mitochondrial genomes")
fo.close()
with open(os.path.join(out, "panel.tsv"), "w") as ft:
    ft.write("id\tclass\taccession\tlength\tdescription\tsource\n")
    for r in rows: ft.write("\t".join(str(x) for x in r) + "\n")
import collections
c = collections.Counter(r[1] for r in rows); b = collections.Counter()
for r in rows: b[r[1]] += r[3]
for k in c: print(f"{k:9s} {c[k]:6d} sequences {b[k] / 1e6:8.2f} Mb")
