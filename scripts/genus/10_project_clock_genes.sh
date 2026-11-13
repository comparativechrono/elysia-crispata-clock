#!/usr/bin/env bash
# Align the E. crispata clock-gene proteins (and the seeds of the genes not found in E. crispata)
# to the genome of every other species with miniprot. Introns of up to 500 kb are allowed, and
# secondary alignments are kept down to 15 % of the best chain, so that paralogues, second
# copies and small pieces of a gene that the assembly has split are all seen.
# usage: 10_project_clock_genes.sh config/species.tsv asm_dir queries.fa outdir
set -uo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
CFG=$1; ASM=$2; Q=$3; OUT=$4; mkdir -p "$OUT"
tail -n +2 "$CFG" | while IFS=$'\t' read -r slug rest; do
  [ -s "$OUT/$slug.mp.done" ] && continue
  miniprot -t2 --gff --trans -G 500000 -p 0.15 --outs 0.15 -N 50 "$ASM/$slug/genome.fa.gz" "$Q" > "$OUT/$slug.mp.gff" 2> "$OUT/$slug.mp.log" && touch "$OUT/$slug.mp.done"
  echo "$slug: $(grep -c $'\tmRNA\t' "$OUT/$slug.mp.gff") alignments"
done
echo PROJECT_DONE
