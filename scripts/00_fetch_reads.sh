#!/usr/bin/env bash
# Download PacBio HiFi runs from the NCBI SRA open-data bucket and store them as
# zstd-compressed FASTA (qualities dropped; read names shortened to RUN.N).
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
OUT=${1:?output dir}; shift
mkdir -p "$OUT/tmp"; cd "$OUT"
for r in "$@"; do
  if [ -s "$r.fa.zst" ] && [ -s "$r.done" ]; then echo "$r already done"; continue; fi
  if [ ! -s "$r.sra" ]; then
    echo "[$(date +%T)] downloading $r"
    curl -sS --fail --retry 5 --retry-delay 5 -C - -o "$r.sra.part" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r"
    mv "$r.sra.part" "$r.sra"
  fi
  echo "[$(date +%T)] converting $r"
  fasterq-dump --fasta-unsorted --stdout --threads 2 --seq-defline '>$ac.$si' -t "$OUT/tmp" "$r.sra" 2>"$r.dump.err" \
    | tee >(awk 'NR%2==0{n++; b+=length($0)} END{print n"\t"b}' > "$r.count") \
    | zstd -3 -T2 -q -o "$r.fa.zst"
  wait
  echo "[$(date +%T)] $r reads/bases: $(cat $r.count)  size: $(du -h $r.fa.zst | cut -f1)"
  rm -f "$r.sra"; touch "$r.done"; echo ok > "$r.done"
done
echo "[$(date +%T)] ALL DONE"
