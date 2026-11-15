#!/usr/bin/env bash
# Search all eleven assemblies for nuclear genes of plants and algae that regulate the chloroplast or
# signal between chloroplast and nucleus (queries from 70_regulator_queries.py), and classify every
# place found, as the protein screen does for photosynthesis genes (scripts/28_protein_screen.py):
# the translated locus is compared with animal, sacoglossan, algal and plant, cyanobacterial,
# bacterial and fungal proteins (the database of scripts/21_alien_index.py) and called
#   animal gene   the slug's own homologue (a family that animals have as well)
#   alga-like     best match an algal or plant protein, better than any animal protein
#   microbe-like  best match a bacterial, archaeal or fungal protein
#   plastid fragment / plastid-only contig / bacterial contig   on DNA that the DNA screen already
#                 identified as kleptoplast or bacterial
# What else lies on every sequence with an alga-like or microbe-like locus is then listed
# (43_foreign_sequences.py).
#
# usage: 71_regulator_screen.sh WORKDIR QUERIES.fa ALIEN_DB.dmnd
#   assemblies: $DATA/asm/<slug>/genome.fa.gz (bgzip + faidx; the alternate haplotype of E. crispata,
#   elysia-crispata-alt, recompressed from GCA_963854135.1 as fetched by the earlier pipeline); the DNA
#   screens of $WORK/asm_screen (and of the earlier E. crispata pipeline, $HGT/asm_screen, for that haplotype)
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
W=$1; Q=$2; DB=$3
DATA=${DATA:-/home/claude/data}; WORK=${WORK:-/home/claude/work/genus}; HGT=${HGT:-/home/claude/work/hgt}
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
mkdir -p "$W"
REF="elysia-crispata elysia-crispata-purdue elysia-crispata-tampa elysia-chlorotica elysia-marginata elysia-timida elysia-leucolegnote elysia-atroviridis elysia-papillosa"
genome(){ echo "$DATA/asm/$1/genome.fa.gz"; }

# what the DNA screen found: kleptoplast stretches inside sequences, and sequences that are plastid DNA throughout
python3 "$HERE/42_protein_screen_inputs.py" "$WORK/asm_screen" "$W" $REF elysia-papillosa-alt
python3 - "$W" "$HGT/asm_screen" <<'EOF'
import json, sys, os
w, h = sys.argv[1:3]
for f in ("embedded_plastid.json", "assembly_screen.json"):
    d = json.load(open(os.path.join(w, f))); old = json.load(open(os.path.join(h, f)))
    d["elysia-crispata-alt"] = old["alt"] if f == "embedded_plastid.json" else dict(whole_list=old["alt"].get("whole_list", []))
    json.dump(d, open(os.path.join(w, f), "w"))
print("inputs for", len(d), "assemblies")
EOF

args=(); gargs=()
for s in $REF elysia-crispata-alt elysia-papillosa-alt; do
  args+=("$s=x"); gargs+=("$s=$(genome "$s")")
  [ -e "$W/$s.done" ] && continue
  if [ "$s" = elysia-papillosa-alt ]; then
    # 1.2 Gb: more than miniprot can index at once in 6 GB of memory
    bash "$HERE/10b_miniprot_in_parts.sh" "$(genome "$s")" "$Q" "$W/$s.gff" 3 -p 0.1 --outs=0.1 -N 200
  else
    miniprot -t2 --gff --trans -p 0.1 --outs=0.1 -N 200 "$(genome "$s")" "$Q" > "$W/$s.gff" 2> "$W/$s.log"
  fi
  touch "$W/$s.done"; echo "[$(date +%T)] $s: $(grep -c $'\tmRNA\t' "$W/$s.gff" || true) alignments"
done
python3 "$REPO/scripts/28_protein_screen.py" "$W" "$DB" "$W/embedded_plastid.json" "$Q" "${args[@]}"
python3 "$HERE/43_foreign_sequences.py" "$W" "$DB" "$W/foreign_sequences.json" "${gargs[@]}"
