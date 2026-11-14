#!/usr/bin/env bash
# miniprot on an assembly that is too large to index at once in the memory at hand: the sequences
# are divided into N parts of about equal size, each part is searched, and the outputs are joined.
# No gene model spans two sequences, so the result is that of one search, except that the cut-off
# for secondary alignments (-p, --outs) is taken relative to the best alignment within a part;
# the steps that follow test every alignment on its own.
# usage: 10b_miniprot_in_parts.sh GENOME.fa.gz QUERIES.fa OUT.gff N [miniprot options]
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
G=$1; Q=$2; OUT=$3; N=$4; shift 4
[ -s "$OUT" ] && { echo "$OUT exists"; exit 0; }
T=$(mktemp -d "${TMPDIR:-/tmp}/mpparts.XXXXXX"); trap 'rm -rf "$T"' EXIT
awk -v n="$N" -v t="$T" 'BEGIN{FS="\t"} {len[NR]=$2; name[NR]=$1; tot+=$2} END{acc=0; for(i=1;i<=NR;i++){p=int(n*acc/tot); if(p>n-1)p=n-1; print name[i] > (t "/part" p ".txt"); acc+=len[i]}}' "$G.fai"
: > "$OUT.tmp"
for i in $(seq 0 $((N - 1))); do
  [ -s "$T/part$i.txt" ] || continue
  samtools faidx -r "$T/part$i.txt" "$G" | bgzip -@2 -c > "$T/part$i.fa.gz"
  flock /tmp/bigmem.lock miniprot -t2 --gff --trans "$@" "$T/part$i.fa.gz" "$Q" >> "$OUT.tmp" 2> "$OUT.part$i.log"
  rm -f "$T/part$i.fa.gz"
done
mv "$OUT.tmp" "$OUT"
echo "$(basename "$OUT"): $(grep -c $'\tmRNA\t' "$OUT") alignments in $N parts"
