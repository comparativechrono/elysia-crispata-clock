#!/usr/bin/env python3
"""Nanopore adapter sequence inside the sequences of an assembly.

usage: 56_adapter_in_assembly.py GENOME.fa[.gz] OUT.tsv

In a Nanopore ligation library every molecule carries the sequencing adapter at both ends. When
two molecules pass through a pore one after the other and are not separated by the software, the
read holds the first molecule, the adapter and the second molecule. An assembler that follows
such a read builds the adapter, and with it a junction between two unrelated molecules, into a
contig. The script lists every place at which an assembly sequence holds the adapter of the
ligation kits (SQK-LSK108/109/110):
    top strand, at the start of a read      ...AATGTACTTCGTTCAGTTACGTATTGCT
    bottom strand, at the end of a read     GCAATACGTAACTGAACGAAGT
A place is reported when the last 24 bases of the top strand are found with no more than 4
differences, or the bottom strand with no more than 3 (substitutions, insertions or deletions;
either orientation). Places within 150 bases of each other are reported once.
Columns: sequence, length, position, distance to the nearer end, differences (top, bottom),
where (end = within 150 bases of an end of the sequence, as left by an untrimmed read; inside),
and the 80 bases around the place.
"""
import sys, gzip, re
fa, out = sys.argv[1:3]
TOP = "AATGTACTTCGTTCAGTTACGTATTGCT"
BOT = "GCAATACGTAACTGAACGAAGT"
COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(COMP)[::-1]
def best(pat, text):
    """smallest number of differences between pat and any stretch of text"""
    prev = [0] * (len(text) + 1)
    for i in range(1, len(pat) + 1):
        cur = [i] + [0] * len(text)
        for j in range(1, len(text) + 1):
            cur[j] = min(prev[j - 1] + (pat[i - 1] != text[j - 1]), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return min(prev)
# exact 12- and 13-base anchors find the places; the whole adapter is then compared with the surrounding sequence
ANCHORS = ["GTTACGTATTGCT", "TACTTCGTTCAG", "GCAATACGTAAC"]
rx = re.compile("|".join(sorted(set(ANCHORS + [rc(a) for a in ANCHORS]))))
T24 = TOP[-24:]
def scan(name, seq):
    S = seq.upper(); hits = []
    for m in rx.finditer(S):
        w = S[max(0, m.start() - 40):m.end() + 40]
        dt = min(best(T24, w), best(rc(T24), w)); db = min(best(BOT, w), best(rc(BOT), w))
        if dt <= 4 or db <= 3: hits.append([name, len(S), m.start(), dt, db, w])
    merged = []
    for h in hits:
        if merged and h[2] - merged[-1][-1] <= 150:
            merged[-1][3] = min(merged[-1][3], h[3]); merged[-1][4] = min(merged[-1][4], h[4]); merged[-1][-1] = h[2]
        else: merged.append(h + [h[2]])
    return [h[:6] for h in merged]
rows = []; name = None; buf = []
with (gzip.open if fa.endswith(".gz") else open)(fa, "rt") as f:
    for line in f:
        if line[0] == ">":
            if name: rows += scan(name, "".join(buf))
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip())
    if name: rows += scan(name, "".join(buf))
with open(out, "w") as fo:
    fo.write("sequence\tlength\tposition\tdistance_to_end\tdifferences_top\tdifferences_bottom\twhere\tcontext\n")
    for n, L, p, dt, db, w in rows:
        d = min(p, L - p)
        fo.write(f"{n}\t{L}\t{p}\t{d}\t{dt}\t{db}\t{'end' if d < 150 else 'inside'}\t{w}\n")
print(f"{len(rows)} places with adapter sequence on {len({r[0] for r in rows})} sequences; {sum(1 for r in rows if min(r[2], r[1] - r[2]) >= 150)} of them more than 150 bases from an end")
