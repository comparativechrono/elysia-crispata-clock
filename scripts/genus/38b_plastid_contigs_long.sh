#!/usr/bin/env bash
# Plastid contigs for a library of accurate long reads whose kleptoplasts have no close reference
# sequence. Reads of such plastid DNA match the panel at 80 to 90% identity, below the 97% at
# which a segment of an accurate read is counted, so the search for junctions would not see them.
#   1. the reads that match a plastid sequence of the panel at 80% identity or more over at least
#      500 bases (pass2.paf) are taken, one in THIN of them (default 30: a few thousand reads are
#      enough, and the time needed grows steeply with the coverage), and assembled with hifiasm
#   2. a contig is kept as plastid DNA when panel plastid sequences cover at least 20% of it or at
#      least 5 kb (minimap2 -x asm20), and the assembly of the slug covers less than 10% of it at
#      95% identity or more
#   3. the kept contigs are written as kcp|NAME|ACCESSION|IDENTITY, where ACCESSION is the panel
#      record that the contig matches over the greatest length and IDENTITY its identity to it:
#      34_read_analysis.py counts the reads on a contig for that reference, at that identity
# The contigs are then given to 33_pass2.sh as a further reference (sixth argument), and the
# library is analysed again.
# usage: 38b_plastid_contigs_long.sh WORKDIR PANEL.fa GENOME.fa.gz
# output: kcp.contigs.fa, kcp.contigs.tsv (contig, length, depth in the whole library, best panel match, identity,
#         share of the contig that matches panel plastid sequences, share that matches the assembly)
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
HIFIASM=${HIFIASM:-/home/claude/tools/mamba/envs/lr/bin/hifiasm}
W=$(realpath "$1"); PANEL=$(realpath "$2"); G=$(realpath "$3"); export THIN=${THIN:-30}
cd "$W"
if [ ! -s hifiasm/kcp.bp.p_ctg.gfa ]; then
  awk -F'\t' '$6 ~ /^(cp|asmcp)\|/ && $11>=500 && $10/$11>=0.80 {print $1}' pass2.paf | sort -u > kcp.names
  seqkit grep -j 2 -w 0 -f kcp.names all_hits.fa | awk -v k="$THIN" 'NR % 2 == 1 {keep = (((NR + 1) / 2) % k == 0)} keep' > kcp.sub.fa
  rm -rf hifiasm; mkdir hifiasm
  "$HIFIASM" -o hifiasm/kcp -t 2 -l0 -f0 kcp.sub.fa > hifiasm/hifiasm.log 2>&1       # -f0: no k-mer filter (a small genome); -l0: nothing purged as a duplicate
fi
awk '/^S/{rd = 0; for (i = 4; i <= NF; i++) if ($i ~ /^rd:i:/) rd = substr($i, 6); print ">" $2 " depth=" rd "\n" $3}' hifiasm/kcp.bp.p_ctg.gfa > kcp.assembled.fa
minimap2 -x asm20 -c --secondary=no -t 2 "$PANEL" kcp.assembled.fa 2>/dev/null | cut -f1-12 > kcp.vs_panel.paf
flock /tmp/bigmem.lock minimap2 -x asm5 -c --secondary=no -t 2 -I 500M "$G" kcp.assembled.fa 2>/dev/null | cut -f1-12 > kcp.vs_genome.paf
python3 - <<'PY'
import collections, os, re
THIN = int(os.environ["THIN"])
def union(iv):
    iv = sorted(iv); tot = 0; ce = -1
    for s, e in iv:
        if s > ce: tot += e - s; ce = e
        elif e > ce: tot += e - ce; ce = e
    return tot
seqs = {}; depth = {}; name = None
for l in open("kcp.assembled.fa"):
    if l[0] == ">": name = l[1:].split()[0]; depth[name] = int(re.search(r"depth=(\d+)", l).group(1)); seqs[name] = ""
    else: seqs[name] += l.strip()
cp = collections.defaultdict(list); best = collections.defaultdict(collections.Counter); ident = collections.defaultdict(collections.Counter)
for l in open("kcp.vs_panel.paf"):
    c = l.split("\t")
    if not c[5].startswith("cp|") or int(c[10]) < 300: continue
    cp[c[0]].append((int(c[2]), int(c[3]))); best[c[0]][c[5]] += int(c[10]); ident[c[0]][c[5]] += int(c[9])
slug = collections.defaultdict(list)
for l in open("kcp.vs_genome.paf"):
    c = l.split("\t")
    if int(c[10]) >= 300 and int(c[9]) / int(c[10]) >= 0.95: slug[c[0]].append((int(c[2]), int(c[3])))
kept = []
with open("kcp.contigs.fa", "w") as fo, open("kcp.contigs.tsv", "w") as ft, open("kcp.rejected.tsv", "w") as fr:
    ft.write("contig\tlength\tdepth\tbest_panel_match\tidentity\tplastid_share\tassembly_share\n")
    fr.write("contig\tlength\tdepth\tbest_panel_match\tidentity\tplastid_share\tassembly_share\n")
    for n in sorted(seqs, key=lambda n: -len(seqs[n])):
        L = len(seqs[n]); pc = union(cp[n]) / L; sc = union(slug[n]) / L
        ref, bp = best[n].most_common(1)[0] if best[n] else (None, 0)
        idn = ident[n][ref] / bp if ref else 0.0
        ok = ref is not None and (pc >= 0.20 or union(cp[n]) >= 5000) and sc < 0.10
        row = f"{L}\t{depth[n] * THIN}\t{ref or ''}\t{idn:.4f}\t{pc:.3f}\t{sc:.3f}\n"
        if ok:
            full = f"kcp|{n}|{ref.split('|', 1)[1]}|{idn:.4f}"
            fo.write(f">{full}\n{seqs[n]}\n"); ft.write(f"{full}\t{row}"); kept.append(L)
        else: fr.write(f"{n}\t{row}")
print(f"{len(seqs)} contigs assembled; {len(kept)} kept as plastid DNA, {sum(kept)} bp")
PY
samtools faidx kcp.contigs.fa
