#!/usr/bin/env bash
# Align each family (MAFFT), trim (trimAl -gappyout) and infer a tree.
# usage: 04_family_trees.sh famdir treedir [fasttree|iqtree]
set -euo pipefail
FAM=$1; OUT=$2; MODE=${3:-fasttree}
mkdir -p "$OUT"
merge(){ # name src1 src2 ...
  local name=$1; shift
  cat "$@" | awk '/^>/{keep=!seen[$1]++} keep' > "$FAM/$name.fa"
}
[ -s "$FAM/CLOCK.fa" ] && merge BHLHPAS "$FAM/CLOCK.fa" "$FAM/BMAL.fa"
[ -s "$FAM/REVERB.fa" ] && merge NR1 "$FAM/REVERB.fa" "$FAM/ROR.fa"
[ -s "$FAM/CK1DE.fa" ] && merge CK1 "$FAM/CK1DE.fa" "$FAM/CK1A.fa"
[ -s "$FAM/VRI.fa" ] && merge BZIP "$FAM/VRI.fa" "$FAM/PARBZIP.fa"
[ -s "$FAM/TIM.fa" ] && merge TIMFAM "$FAM/TIM.fa" "$FAM/TIMEOUT.fa"
for f in "$FAM"/*.fa; do
  n=$(basename "$f" .fa)
  case $n in CLOCK|BMAL|REVERB|ROR|CK1DE|CK1A|VRI|PARBZIP|TIM|TIMEOUT) continue;; esac
  [ "$(grep -c '>' "$f")" -lt 4 ] && continue
  [ -s "$OUT/$n.nwk" ] && continue
  mafft --auto --thread 2 --quiet "$f" > "$OUT/$n.aln"
  trimal -in "$OUT/$n.aln" -out "$OUT/$n.trim.aln" -gappyout 2>/dev/null || cp "$OUT/$n.aln" "$OUT/$n.trim.aln"
  if [ "$MODE" = iqtree ]; then
    iqtree -s "$OUT/$n.trim.aln" -m MFP -mset LG,WAG,JTT -B 1000 -T 2 --prefix "$OUT/$n" -quiet -redo && cp "$OUT/$n.treefile" "$OUT/$n.nwk"
  else
    FastTree -lg -gamma -quiet "$OUT/$n.trim.aln" > "$OUT/$n.nwk" 2>/dev/null
  fi
  echo "$n $(grep -c '>' "$f") seqs, $(awk '/^>/{if(s){print length(s);exit}; next}{s=s $0}' "$OUT/$n.trim.aln") trimmed cols"
done
