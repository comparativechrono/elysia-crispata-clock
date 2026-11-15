#!/usr/bin/env python3
"""Reads that span a barcode gene but were set aside when only a sample of the plastid reads was carried on.

usage: 36b_barcode_reads_set_aside.py WORKDIR PANEL.fa FRAMES.fa OUT.fa

Where the first pass keeps one in K of the reads of plastid DNA (31_select_hits.py, thin.tsv), few
of the kept reads span rbcL or tufA. The barcode genes are read from all reads that span them:
FRAMES.fa (WORKDIR/barcodes/frames.fa of a first run of 36_barcodes.py: the record each barcode is
framed by) is placed on the plastid sequences of the panel, the first-pass matches (RUN.p1.paf)
give the reads that cover such a place (all but 15% at either end), and those that are not among the kept
reads are fetched again by name (52_fetch_reads.py). OUT.fa is then given to 36_barcodes.py
together with the kept reads.
"""
import sys, os, glob, subprocess
W, panel, frames, out = sys.argv[1:5]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
HERE = os.path.dirname(os.path.abspath(__file__))
places = []          # (panel sequence, start, end) of each barcode gene
p = subprocess.run([BIN + "/minimap2", "-x", "map-ont", "-k13", "-w5", "-c", "--secondary=yes", "-N", "50", "-t", "2", panel, frames], capture_output=True, text=True, check=True).stdout
for l in p.splitlines():
    c = l.split("\t")
    if int(c[10]) >= 0.75 * int(c[1]): places.append((c[5], int(c[7]), int(c[8])))
kept = set()
for f in glob.glob(os.path.join(W, "*.keep.names")): kept |= {l.strip() for l in open(f)}
want = set()
for f in glob.glob(os.path.join(W, "*.p1.paf")):
    for l in open(f):
        c = l.split("\t", 9)
        for t, a, b in places:
            m = 0.15 * (b - a)          # a first-pass match need not reach the very ends of the gene on a distant reference
            if c[5] == t and int(c[7]) <= a + m and int(c[8]) >= b - m: want.add(c[0])
want -= kept
names = out + ".names"
open(names, "w").write("".join(n + "\n" for n in sorted(want)))
print(f"{len(places)} places of barcode genes on the panel; {len(want)} reads that span one were set aside")
if want: sys.exit(subprocess.run([sys.executable, os.path.join(HERE, "52_fetch_reads.py"), names, out]).returncode)
open(out, "w").close()
