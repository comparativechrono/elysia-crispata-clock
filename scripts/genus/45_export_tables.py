#!/usr/bin/env python3
"""Copy the small intermediate tables of the read and gene-transfer screens into the repository.

usage: 45_export_tables.py REPO WORK [READS_ROOT=WORK/reads]
results/genus/
  panel.tsv.gz                      every sequence of the reference panel (accession, class, length, title)
  barcode_records.tsv.gz            the GenBank rbcL and tufA records the barcodes were compared with
  <library>.junctions.tsv.gz        long reads: every join between two kinds of DNA inside one read
  <library>.junction_pairs.tsv.gz   short reads: every pair and every mate that joins organelle DNA to slug DNA
  <library>.insertion_shapes.tsv    long reads: every read with organelle DNA inside slug DNA, or the mirror image
                                    (the three tables are also in docs/downloads)
  <library>.plastid_contigs.fa/.tsv short reads: plastid contigs assembled from the reads (also in docs/downloads)
  <library>.barcodes.fa             each distinct rbcL and tufA sequence found (also in docs/downloads)
  protein_screen_loci.json.gz       every place matched by one of the 471 proteins, with its classification
  foreign_sequences.json.gz         every coding stretch on the sequences that carry non-animal matches
"""
import sys, os, gzip, shutil, glob, csv
repo, work = sys.argv[1:3]
reads_root = sys.argv[3] if len(sys.argv) > 3 else os.path.join(work, "reads")
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
for d in sorted(glob.glob(os.path.join(reads_root, "elysia-*"))):
    if not os.path.isdir(d): continue
    lib = os.path.basename(d); slug = lib.split(".")[0]
    if slug == "elysia-crispata": continue                      # its tables are made by scripts/26_collect_results.py
    for old_name in (f"{lib}.junctions.tsv", f"{lib}.locus372.json"):          # tables of an earlier layout
        if os.path.exists(os.path.join(RG, old_name)): os.remove(os.path.join(RG, old_name))
    for src, dst in (("kcp.contigs.fa", f"{lib}.plastid_contigs.fa"), ("kcp.contigs.tsv", f"{lib}.plastid_contigs.tsv")):
        p = os.path.join(d, src)
        if os.path.exists(os.path.join(RG, dst)) and os.path.exists(p) and sum(1 for _ in open(p)) < 2: os.remove(os.path.join(RG, dst))          # a table with nothing but its header
        if os.path.exists(p) and os.path.getsize(p) and sum(1 for _ in open(p)) >= 2:
            shutil.copy(p, os.path.join(RG, dst)); n += 1
            if src == "kcp.contigs.fa": shutil.copy(p, os.path.join(DL, lib.replace(".", "_") + "_plastid_contigs.fa"))
    tag = lib.replace(".", "_")
    if os.path.exists(os.path.join(d, "analysis.json")) and gz(os.path.join(d, "junctions.tsv"), os.path.join(RG, f"{lib}.junctions.tsv.gz")):
        shutil.copy(os.path.join(RG, f"{lib}.junctions.tsv.gz"), os.path.join(DL, tag + "_junctions.tsv.gz")); n += 1
    p = os.path.join(d, "junctions.tsv.gz")                     # short reads: written compressed by 37_pair_analysis.py
    if os.path.exists(p):
        import gzip as _g
        with _g.open(p, "rb") as fi, gzip.GzipFile(os.path.join(RG, f"{lib}.junction_pairs.tsv.gz"), "wb", mtime=0) as fo: shutil.copyfileobj(fi, fo)
        shutil.copy(os.path.join(RG, f"{lib}.junction_pairs.tsv.gz"), os.path.join(DL, tag + "_junction_pairs.tsv.gz")); n += 1
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
# reads with the shape of an insertion, from the collected results
for tf in sorted(glob.glob(os.path.join(RG, "*.transfer.json"))):
    T = json.load(open(tf)); slug = T["slug"]
    for R in T.get("reads", []):
        S = R.get("shapes")
        if not S or not S.get("candidates"): continue
        lib = slug if R["library"] == "main" else f'{slug}.{R["library"]}'
        cols = ["read", "read_length", "shape", "insert", "how", "flank_sequence", "left_from", "left_to", "right_from", "right_to", "flank_strand", "flank_mapq", "flank_identity", "reference_gap", "read_gap",
                "piece_sequence", "piece_from", "piece_to", "piece_bp", "piece_identity", "place_position", "other_reads_crossing", "second_pass", "assembly_locus"]
        out = os.path.join(RG, lib + ".insertion_shapes.tsv")
        with open(out, "w") as fo:
            fo.write("\t".join(cols) + "\n")
            for c in S["candidates"]:
                c = dict(c, place_position=(c.get("place") or {}).get("position"))
                fo.write("\t".join("" if c.get(k) is None else str(c[k]) for k in cols) + "\n")
        shutil.copy(out, os.path.join(DL, lib.replace(".", "_") + "_insertion_shapes.tsv")); n += 1
print(n, "tables copied to", RG)
