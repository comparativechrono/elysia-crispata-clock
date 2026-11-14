#!/usr/bin/env bash
# Assemble the plastid reads of a short-read species into contigs (MEGAHIT), keep the contigs
# that are plastid DNA, and name them from rbcL and tufA.
#   pairs used: every pair with a mate that matches a plastid sequence of the panel (>= 80 bases at
#   >= 80% identity) in pass 2
# usage: 38_plastid_contigs.sh WORKDIR PANEL.fa barcodes.fa
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
MEGAHIT=${MEGAHIT:-/home/claude/tools/mamba/envs/asm/bin/megahit}
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(realpath "$1"); PANEL=$(realpath "$2"); BC=$(realpath "$3")
cd "$W"
awk -F'\t' '$6 ~ /^(cp|asmcp)\|/ && $11>=80 && $10/$11>=0.80 {n=$1; sub(/\.[12]$/,"",n); print n}' pass2.paf | sort -u > kcp.spots
awk '{print $1".1"; print $1".2"}' kcp.spots > kcp.names
seqkit grep -j 2 -w 0 -f kcp.names all_hits.fa > kcp.pairs.fa
rm -rf megahit; "$MEGAHIT" --12 kcp.pairs.fa -o megahit --min-contig-len 500 -t 2 -m 0.3 > megahit.log 2>&1
# keep contigs that are plastid: >= 50% of the contig aligns to a panel plastid sequence (any identity minimap2 asm20 reports)
minimap2 -x asm20 -c --secondary=no -t 2 "$PANEL" megahit/final.contigs.fa 2>/dev/null | cut -f1-12 > kcp.vs_panel.paf
python3 - <<'PY'
import collections, re
cov = collections.defaultdict(list); L = {}; best = {}
for l in open("kcp.vs_panel.paf"):
    c = l.split("\t")
    if not c[5].startswith("cp|"): continue
    cov[c[0]].append((int(c[2]), int(c[3]))); L[c[0]] = int(c[1])
    if c[0] not in best or int(c[9]) > best[c[0]][0]: best[c[0]] = (int(c[9]), c[5], int(c[9]) / int(c[10]))
def union(iv):
    iv = sorted(iv); t = 0; cur = None
    for a, b in iv:
        if cur is None: cur = [a, b]
        elif a <= cur[1]: cur[1] = max(cur[1], b)
        else: t += cur[1] - cur[0]; cur = [a, b]
    return t + (cur[1] - cur[0] if cur else 0)
keep = {c for c in cov if union(cov[c]) >= 0.5 * L[c]}
n = 0; tot = 0
with open("kcp.contigs.fa", "w") as fo, open("kcp.contigs.tsv", "w") as ft:
    ft.write("contig\tlength\tdepth\tbest_panel_match\tidentity\n"); name = None
    for l in open("megahit/final.contigs.fa"):
        if l[0] == ">":
            name = l[1:].split()[0]; m = re.search(r"multi=([\d.]+)", l); ln = re.search(r"len=(\d+)", l)
            if name in keep:
                n += 1; tot += int(ln.group(1)); fo.write(f">kcp|{name}\n"); ft.write(f"kcp|{name}\t{ln.group(1)}\t{m.group(1)}\t{best[name][1]}\t{best[name][2]:.4f}\n")
        elif name in keep: fo.write(l)
print(f"plastid contigs: {n}, {tot} bp")
PY
samtools faidx kcp.contigs.fa
python3 "$HERE/36_barcodes.py" contigs kcp.contigs.fa "$BC" barcodes
