#!/usr/bin/env bash
# Pass 1 of the read-level screen: map every HiFi read against the "bait" panel of
# non-slug sequences (algal plastomes and mitogenomes, algal nuclear rDNA, the four
# bacterial MAGs, the slug mitochondrion and slug rDNA).  Approximate mapping only
# (no base-level alignment); hit reads are re-aligned properly in pass 2.
# usage: 10_bait_screen.sh panel.mmi reads_dir out_dir RUN [RUN ...]
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
MMI=$1; READS=$2; OUT=$3; shift 3
mkdir -p "$OUT"
for r in "$@"; do
  [ -s "$OUT/$r.bait.paf.done" ] && continue
  echo "[$(date +%T)] $r start"
  zstd -dc "$READS/$r.fa.zst" | minimap2 -x map-ont -k15 -w10 --secondary=no -t2 "$MMI" - 2>"$OUT/$r.bait.log" > "$OUT/$r.bait.paf"
  echo ok > "$OUT/$r.bait.paf.done"
  echo "[$(date +%T)] $r done: $(wc -l < "$OUT/$r.bait.paf") alignments"
done
echo BAIT_SCREEN_DONE
