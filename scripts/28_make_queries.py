#!/usr/bin/env python3
"""Build the query set for the protein screen (script 28).
usage: 28_make_queries.py out.fa TAG=fasta[.gz] ...
e.g.   28_make_queries.py queries.fa UDOFL=klepto/cp_proteins.fa CHLRE=prot/chlre_photosynthesis.fa.gz \\
           ARATH=prot/arath_photosynthesis.fa.gz BRYO=prot/bryopsidales_uniprot.fa.gz VAULI=prot/vaucheria_all.fa.gz
Headers become TAG|accession gene description; proteins shorter than 40 aa are dropped.
"""
import sys, gzip
def fa(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for l in fh:
            if l[0] == ">":
                if name: yield name, "".join(seq)
                name, seq = l[1:].rstrip(), []
            else: seq.append(l.strip())
    if name: yield name, "".join(seq)
with open(sys.argv[1], "w") as out:
    for a in sys.argv[2:]:
        tag, path = a.split("=", 1); k = 0
        for h, s in fa(path):
            s = s.replace("*", "")
            if len(s) < 40: continue
            first = h.split()[0]; acc = first.split("|")[1] if "|" in first else first
            desc = h.split(" ", 1)[1].split(" OS=")[0] if " " in h else ""
            gn = h.split(" GN=")[1].split()[0] if " GN=" in h else ""
            out.write(f">{tag}|{acc} {gn + ' ' if gn else ''}{desc}\n{s}\n"); k += 1
        print(tag, k)
