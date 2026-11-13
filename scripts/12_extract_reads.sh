#!/usr/bin/env bash
# Pull the reads selected by 11_classify_bait_hits.py out of the compressed read sets.
# usage: 12_extract_reads.sh reads_dir sel_dir
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
READS=$1; SEL=$2
cat "$SEL"/cp_like.names "$SEL"/control.names "$SEL"/clock.names | sort -u > "$SEL/all.names"
: > "$SEL/selected.fa"
for z in "$READS"/*.fa.zst; do
  zstd -dc "$z" | seqkit grep -j 2 -w 0 -f "$SEL/all.names" >> "$SEL/selected.fa"
done
for k in cp_like control clock; do
  seqkit grep -j 2 -w 0 -f "$SEL/$k.names" "$SEL/selected.fa" > "$SEL/$k.fa"
  echo "$k: $(grep -c '>' "$SEL/$k.fa") reads"
done
rm -f "$SEL/selected.fa"
