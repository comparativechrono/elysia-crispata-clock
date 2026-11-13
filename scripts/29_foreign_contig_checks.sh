#!/usr/bin/env bash
# Follow-up for contigs that carry alga-like or microbe-like genes outside the known plastid
# fragments (found by script 28). For the contigs of one assembly:
#   - base composition, and the gene models annotated on them (do any match any protein?)
#   - do they align to the other assemblies of the species?            (minimap2 asm20)
#   - do the HiFi reads of the xgElyCris1 animal contain them?         (minimap2 map-hifi)
# A contig that is part of the slug genome should be found in other animals; DNA of an alga
# the sequenced animal had eaten should not.
# usage: 29_foreign_contig_checks.sh workdir asm_name asm.fa feature_table.gz alien_hits.tsv reads_dir "RUN ..." other=fasta ...
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
HERE=$(cd "$(dirname "$0")" && pwd)
WD=$1; ASM=$2; FA=$3; FT=$4; HITS=$5; READS=$6; RUNLIST=$7; shift 7
cd "$WD"
python3 - "$ASM" <<'PY'
import json, sys
R = json.load(open("protein_screen_raw.json"))
c = sorted({l["contig"] for l in R["loci"] if l["assembly"] == sys.argv[1] and l["class"] in ("alga-like", "microbe-like")})
open(f"{sys.argv[1]}_foreign.names", "w").write("\n".join(c) + "\n")
PY
seqkit grep -w 0 -f "${ASM}_foreign.names" "$FA" > "${ASM}_foreign.fa"
for a in "$@"; do
  n=${a%%=*}; fa=${a#*=}
  minimap2 -x asm20 -t 2 --secondary=no "$fa" "${ASM}_foreign.fa" 2> /dev/null | cut -f1-12 > "${ASM}_foreign_vs_$n.paf"
done
# read scan: only contigs not scanned before (the scan is the slow step)
touch "reads_vs_${ASM}_foreign.scanned" "reads_vs_${ASM}_foreign.paf"
grep -v -x -F -f "reads_vs_${ASM}_foreign.scanned" "${ASM}_foreign.names" > "${ASM}_foreign.new.names" || true
if [ -s "${ASM}_foreign.new.names" ]; then
  seqkit grep -w 0 -f "${ASM}_foreign.new.names" "$FA" > "${ASM}_foreign.new.fa"
  for r in $RUNLIST; do
    zstd -dc "$READS/$r.fa.zst" | minimap2 -x map-hifi --secondary=no -t 2 "${ASM}_foreign.new.fa" - 2>> "reads_vs_${ASM}_foreign.log" | awk '$11>=500' | cut -f1-12 >> "reads_vs_${ASM}_foreign.paf"
  done
  cat "${ASM}_foreign.new.names" >> "reads_vs_${ASM}_foreign.scanned"
fi
python3 "$HERE/29_foreign_contigs.py" "$WD" "$ASM" "$FT" "$HITS" "$@"
