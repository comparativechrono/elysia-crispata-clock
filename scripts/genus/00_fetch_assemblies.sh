#!/usr/bin/env bash
# Fetch the genome assembly (and GenBank annotation, where there is one) of every species in
# config/species.tsv from the NCBI genomes FTP site. Assemblies are recompressed with bgzip so
# that samtools faidx, minimap2 and miniprot can all read the same file.
# usage: 00_fetch_assemblies.sh config/species.tsv outdir
set -uo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
CFG=$1; OUT=$2; mkdir -p "$OUT"
tail -n +2 "$CFG" | while IFS=$'\t' read -r slug species tag acc name level ann rest; do
  d="$OUT/$slug"; mkdir -p "$d"
  p=${acc:4:3}/${acc:7:3}/${acc:10:3}; base="${acc}_${name// /_}"
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/$p/$base"
  if [ ! -s "$d/genome.fa.gz.fai" ]; then
    curl -sS --retry 5 --retry-delay 5 -C - -o "$d/genomic.fna.gz" "$url/${base}_genomic.fna.gz" || { echo "FAILED $slug genome"; continue; }
    zcat "$d/genomic.fna.gz" | bgzip -@2 -c > "$d/genome.fa.gz" && samtools faidx "$d/genome.fa.gz" && rm -f "$d/genomic.fna.gz"
  fi
  curl -sS --retry 3 -o "$d/assembly_report.txt" "$url/${base}_assembly_report.txt" || true
  if [ "$ann" = genbank ]; then
    for f in genomic.gff.gz protein.faa.gz; do [ -s "$d/$f" ] || curl -sS --retry 3 -o "$d/$f" "$url/${base}_$f" || true; done
  fi
  echo "$slug: $(wc -l < "$d/genome.fa.gz.fai") sequences, $(awk '{s+=$2} END{printf "%.1f", s/1e6}' "$d/genome.fa.gz.fai") Mb"
done
echo FETCH_DONE
