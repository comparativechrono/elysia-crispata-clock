#!/usr/bin/env bash
# Assembly-level screen for organelle-derived DNA: every E. crispata assembly is cut into
# 10 kb windows (1 kb overlap) and mapped against the bait panel (algal plastomes and
# mitogenomes, algal rDNA, bacterial MAGs, slug mitochondrion).
# usage: 22_assembly_organelle_screen.sh panel.mmi outdir name=fasta[.gz] ...
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
MMI=$1; OUT=$2; shift 2
mkdir -p "$OUT"
for kv in "$@"; do
  n=${kv%%=*}; f=${kv#*=}
  [ -s "$OUT/$n.paf" ] && continue
  seqkit sliding -j 2 -W 10000 -s 9000 -g "$f" 2>/dev/null | sed 's/ .*//' \
    | minimap2 -x map-ont -k15 -w10 -c --secondary=no -t 2 "$MMI" - 2> "$OUT/$n.log" \
    | awk 'BEGIN{FS=OFS="\t"} {o=$1; for(i=2;i<=NF;i++) if ($i !~ /^c[gs]:Z:/) o=o OFS $i; print o}' > "$OUT/$n.paf.tmp"
  mv "$OUT/$n.paf.tmp" "$OUT/$n.paf"
  echo "$n: $(wc -l < "$OUT/$n.paf") alignments"
done
echo ASM_SCREEN_DONE
