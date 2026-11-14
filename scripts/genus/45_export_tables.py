#!/usr/bin/env python3
"""Copy the small intermediate tables of the read and gene-transfer screens into the repository.

usage: 45_export_tables.py REPO WORK
results/genus/
  panel.tsv.gz                      every sequence of the reference panel (accession, class, length, title)
  barcode_records.tsv.gz            the GenBank rbcL and tufA records the barcodes were compared with
  <library>.junctions.tsv           long reads: every join between two kinds of DNA inside one read
  <library>.plastid_contigs.fa/.tsv short reads: plastid contigs assembled from the reads (also in docs/downloads)
  <library>.barcodes.fa             each distinct rbcL and tufA sequence found (also in docs/downloads)
  protein_screen_loci.json.gz       every place matched by one of the 471 proteins, with its classification
  foreign_sequences.json.gz         every coding stretch on the sequences that carry non-animal matches
  elysia-chlorotica.eggs.locus372.json   egg reads at the 86-base plastid-like element of RQTK01000372.1
"""
import sys, os, gzip, shutil, glob, csv
repo, work = sys.argv[1:3]
RG = os.path.join(repo, "results", "genus"); DL = os.path.join(repo, "docs", "downloads")
os.makedirs(RG, exist_ok=True); os.makedirs(DL, exist_ok=True)
def gz(src, dst):
    if os.path.exists(src):
        with open(src, "rb") as fi, gzip.GzipFile(dst, "wb", mtime=0) as fo: shutil.copyfileobj(fi, fo)
        return True
gz(os.path.join(work, "panel", "panel.tsv"), os.path.join(RG, "panel.tsv.gz"))
gz(os.path.join(work, "protscreen", "protein_screen_raw.json"), os.path.join(RG, "protein_screen_loci.json.gz"))
gz(os.path.join(work, "protscreen", "foreign_sequences.json"), os.path.join(RG, "foreign_sequences.json.gz"))
bc = os.environ.get("BARCODES", "/home/claude/data/genus_refs/barcodes/barcodes.fa")
if os.path.exists(bc):
    with gzip.GzipFile(os.path.join(RG, "barcode_records.tsv.gz"), "wb", mtime=0) as fo:
        fo.write(b"record\ttitle\n")
        for line in open(bc):
            if line[0] == ">":
                f = line[1:].rstrip().split(" ", 1); fo.write((f[0] + "\t" + (f[1] if len(f) > 1 else "") + "\n").encode())
n = 0
for d in sorted(glob.glob(os.path.join(work, "reads", "elysia-*"))):
    if not os.path.isdir(d): continue
    lib = os.path.basename(d); slug = lib.split(".")[0]
    for src, dst in (("junctions.tsv", f"{lib}.junctions.tsv"), ("kcp.contigs.fa", f"{lib}.plastid_contigs.fa"), ("kcp.contigs.tsv", f"{lib}.plastid_contigs.tsv"),
                     (os.path.join("locus372", "locus_test.json"), f"{lib}.locus372.json")):
        p = os.path.join(d, src)
        if os.path.exists(p) and os.path.getsize(p):
            shutil.copy(p, os.path.join(RG, dst)); n += 1
            if src == "kcp.contigs.fa": shutil.copy(p, os.path.join(DL, lib.replace(".", "_") + "_plastid_contigs.fa"))
# barcode sequences, with headers that say what they are (from the collected results of 39_collect_species.py)
import json, re
for kf in sorted(glob.glob(os.path.join(RG, "*.kleptoplast.json"))):
    K = json.load(open(kf)); slug = K["slug"]
    for L in K["libraries"]:
        if not L.get("barcodes"): continue
        lib = slug if L["library"] == "main" else f'{slug}.{L["library"]}'
        out = os.path.join(RG, lib + ".barcodes.fa")
        with open(out, "w") as fo:
            for b in L["barcodes"]:
                src = f'reads={b["reads"]}' if b.get("reads") else f'from_contig_of={b.get("contig_length")}bp'
                fo.write(f'>{lib}|{b["barcode"]}|{b["haplotype"].split("|")[-1]} {src} length={b["length"]} longest_good_match="{b["top_hit"]}" identity={b["identity"]}%\n')
                fo.write(re.sub(r"(.{70})", r"\1\n", b["sequence"]).rstrip("\n") + "\n")
        shutil.copy(out, os.path.join(DL, lib.replace(".", "_") + "_barcodes.fa")); n += 1
print(n, "tables copied to", RG)
