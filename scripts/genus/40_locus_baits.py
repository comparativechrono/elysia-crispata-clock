#!/usr/bin/env python3
"""Baits for testing, with reads, the places where an assembly holds organelle-like DNA.

usage: 40_locus_baits.py SCREEN.json GENOME.fa.gz OUT.fa [flank=3000] [contig:start-end:kind ...]
One bait per locus of 35_assembly_screen.py (plastid-like loci of >= 90% identity and the
mitochondrial NUMT-like loci, which serve as a positive control): the locus with FLANK bases
either side. Bases that the submitters soft-masked as repeats (lower case) are replaced by N so
that the bait collects reads of this locus only. Names: bait|contig:start-end|kind|locus_start-locus_end
Further places can be given by hand (0-based start, kind plastid or numt), for stretches that the
DNA screen reports below 90% identity or that were found in another way.
"""
import sys, json, subprocess, os
screen, genome, out = sys.argv[1:4]; flank = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
extra = sys.argv[5:]
ST = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/samtools"
s = json.load(open(screen)); loci = []
for L in s["loci"]:
    if L["cls"] == "plastid" and L["identity"] >= 0.90: loci.append((L["contig"], L["start"], L["end"], "plastid", L["contig_length"]))
for N in s["numt_like"]: loci.append((N["contig"], N["start"], N["end"], "numt", None))
for x in extra:
    c, r, kind = x.rsplit(":", 2); a, b = (int(v) for v in r.split("-"))
    if not any(l[0] == c and l[1] == a and l[2] == b for l in loci): loci.append((c, a, b, kind, None))
with open(out, "w") as fo, open(out + ".tsv", "w") as ft:
    ft.write("bait\tcontig\twindow_start\twindow_end\tkind\tlocus_start\tlocus_end\n")
    for c, a, b, kind, L in loci:
        ws = max(0, a - flank); we = b + flank
        fa = subprocess.run([ST, "faidx", genome, f"{c}:{ws + 1}-{we}"], capture_output=True, text=True, check=True).stdout
        seq = "".join(fa.split("\n")[1:]); we = ws + len(seq)
        masked = "".join(ch if ch.isupper() else "N" for ch in seq)
        name = f"bait|{c}:{ws}-{we}|{kind}|{a}-{b}"
        fo.write(f">{name}\n{masked}\n"); ft.write(f"{name}\t{c}\t{ws}\t{we}\t{kind}\t{a}\t{b}\n")
        print(name, len(seq), "bp,", sum(1 for ch in masked if ch == "N"), "masked")
