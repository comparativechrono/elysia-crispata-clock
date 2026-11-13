#!/usr/bin/env bash
# Quantify the head RNA-seq library (ERR12342469, same specimen as the genome) against the
# coding sequences of the Purdue E. crispata annotation with salmon.
# Reads are dumped from the SRA archive as two gzipped FASTA files (qualities are not needed),
# the archive is removed to save disk, and salmon is run on the pair.
# usage: 20_rnaseq_quant.sh run.sra salmon_index outdir
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
SRA=$1; IDX=$2; OUT=$3
mkdir -p "$OUT"; cd "$OUT"
if [ ! -s r2.fa.gz ] || [ -s "$SRA" ]; then
  fastq-dump --fasta 0 --split-spot --skip-technical -Z "$SRA" 2> dump.err \
    | awk -v o1="pigz -1 -p 1 > r1.fa.gz" -v o2="pigz -1 -p 1 > r2.fa.gz" '{ if (int((NR-1)/2) % 2 == 0) print | o1; else print | o2 }'
  n1=$(zcat r1.fa.gz | grep -c '>'); n2=$(zcat r2.fa.gz | grep -c '>')
  echo "dumped $n1 / $n2 mates"
  [ "$n1" = "$n2" ] && rm -f "$SRA"
fi
salmon quant -i "$IDX" -l A -1 r1.fa.gz -2 r2.fa.gz -p 2 --validateMappings -o salmon > salmon.log 2>&1
grep -E "Mapping rate|most likely library type" salmon/logs/salmon_quant.log || true
rm -f r1.fa.gz r2.fa.gz
echo RNASEQ_DONE
