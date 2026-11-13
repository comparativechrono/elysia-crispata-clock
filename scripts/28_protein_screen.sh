#!/usr/bin/env bash
# Search the assemblies directly for algal genes, independent of any gene annotation.
#   queries   plastid-encoded proteins of Udotea flabellum (NC_039528.1, the closest annotated
#             relative of the Penicillus kleptoplasts), photosynthesis proteins of Chlamydomonas and
#             Arabidopsis (nuclear- and plastid-encoded), UniProt proteins of Bryopsidales, and all
#             Vaucheria litorea proteins (the alga of the original Elysia chlorotica gene-transfer
#             claims); assembled by 28_make_queries.py
#   step 1    miniprot: spliced protein-to-genome alignment of every query to every assembly.
#             Secondary hits are kept down to 10 % of the best score (-p 0.1 --outs=0.1), so
#             that a weaker or partial copy is not hidden behind a better copy elsewhere
#   step 2    overlapping hits are merged into loci; the translation of the best hit per locus is
#             searched against the same animal / algal / bacterial / fungal protein database as
#             script 21. A locus whose best match is an animal protein is the slug's own homologue
#             (for example mitochondrial ATP synthase found with chloroplast atpB).
# usage: 28_protein_screen.sh workdir alien_db.dmnd embedded_plastid.json queries.fa name=fasta[.gz] ...
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
HERE=$(cd "$(dirname "$0")" && pwd)
WD=$1; DB=$2; EMB=$3; Q=$4; shift 4
mkdir -p "$WD"
for a in "$@"; do
  n=${a%%=*}; fa=${a#*=}
  [ -s "$WD/$n.done" ] && continue
  miniprot -t2 --gff --trans -p 0.1 --outs=0.1 -N 200 "$fa" "$Q" > "$WD/$n.gff" 2> "$WD/$n.log"
  touch "$WD/$n.done"; echo "$n: $(grep -c $'\tmRNA\t' "$WD/$n.gff" || true) alignments"
done
python3 "$HERE/28_protein_screen.py" "$WD" "$DB" "$EMB" "$Q" "$@"
