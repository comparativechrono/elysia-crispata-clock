#!/usr/bin/env python3
"""Build one tagged protein FASTA (+ annotation TSV) from the Ensembl E. crispata
peptides and the UniProt comparison proteomes.

usage: 02_build_protein_db.py ensembl_pep.fa.gz refprot_dir out_prefix
Headers become  >TAG|ID|GENE  where TAG is the species tag (ELYCR for the Sanger
assembly's Ensembl proteins, ELYCRP for the Purdue-assembly proteome).
"""
import sys, gzip, glob, os, re

ens, refdir, out = sys.argv[1:4]

def fasta(path):
    name, seq = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name: yield name, "".join(seq)
                name, seq = line[1:].rstrip(), []
            else:
                seq.append(line.strip())
    if name: yield name, "".join(seq)

n = {}
with open(out + ".fa", "w") as fo, open(out + ".tsv", "w") as ft:
    ft.write("id\tspecies\taccession\tgene\tlength\tdescription\n")
    for h, s in fasta(ens):
        pid = h.split()[0].split(".")[0]
        gene = re.search(r"gene:(\S+)", h).group(1).split(".")[0]
        loc = h.split()[2]
        i = f"ELYCR|{pid}|{gene}"
        fo.write(f">{i}\n{s}\n"); ft.write(f"{i}\tELYCR\t{pid}\t{gene}\t{len(s)}\t{loc}\n")
        n["ELYCR"] = n.get("ELYCR", 0) + 1
    for path in sorted(glob.glob(os.path.join(refdir, "*.fa.gz"))):
        tag = os.path.basename(path)[:-6]
        tag = {"ELYCR_PURDUE": "ELYCRP"}.get(tag, tag)
        for h, s in fasta(path):
            acc = h.split("|")[1]
            m = re.search(r" GN=(\S+)", h)
            gene = m.group(1) if m else ""
            desc = re.sub(r" OS=.*", "", h.split(" ", 1)[1]) if " " in h else ""
            i = f"{tag}|{acc}|{gene}"
            fo.write(f">{i}\n{s}\n"); ft.write(f"{i}\t{tag}\t{acc}\t{gene}\t{len(s)}\t{desc}\n")
            n[tag] = n.get(tag, 0) + 1
for k, v in n.items():
    print(k, v)
