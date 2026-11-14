#!/usr/bin/env bash
# Pass 2 of the read screen: base-level alignment of the reads kept by pass 1 against the
# species' own assembly plus the reference panel.
#   - organelle sequences that are part of the assembly (35_assembly_screen.py, "whole") are
#     renamed CLASS|name so that they count as organelle, not slug
#   - plastid-like stretches inside slug sequences (mask.bed) are hidden, so that a plastid read
#     is not explained away by a plastid fragment that the assembler left in a slug sequence
#   - of the alignment strings only those with an insertion of 30 bases or more are kept (cg tag):
#     46_candidates.py looks inside them for organelle DNA too short to break an alignment
# The index is built in parts (PASS2_INDEX_BASES, default 500M) to keep memory low.
# The first pass must be finished: the script stops if a run of WORKDIR is still being screened.
# usage: 33_pass2.sh WORKDIR GENOME.fa.gz PANEL.fa SCREEN.json hifi|clr|short [EXTRA_REF.fa]
#        SAMPLE_ONLY=1   align only the random sample (RUN.sample.fa -> sample.paf)
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
W=$(realpath "$1"); G=$(realpath "$2"); PANEL=$(realpath "$3"); SCREEN=$(realpath "$4"); KIND=$5; EXTRA=${6:+$(realpath "$6")}
cd "$W"
if ls -d ./*.parts > /dev/null 2>&1 || ls ./*.hits.fa.tmp > /dev/null 2>&1; then echo "a run in $W is still being screened" >&2; exit 1; fi
cat ./*.hits.fa > all_hits.fa.tmp && mv all_hits.fa.tmp all_hits.fa
python3 - "$SCREEN" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))
CL = {"plastid": "asmcp", "algal mitochondrial": "asmalgmt", "slug mitochondrial": "slugmt"}
with open("relabel.tsv", "w") as fo:
    for w in s["whole"]: fo.write(f"{w['contig']}\t{CL[w['cls']]}|{w['contig']}\n")
with open("mask.bed", "w") as fo:
    for L in s["loci"]:
        if L["cls"] == "plastid" and L["identity"] >= 0.90: fo.write(f"{L['contig']}\t{max(0, L['start'] - 20)}\t{L['end'] + 20}\n")
PY
zcat "$G" | sed 's/ .*//' | awk 'BEGIN{while((getline l < "relabel.tsv")>0){split(l,a,"\t"); r[a[1]]=a[2]}} /^>/{n=substr($0,2); if(n in r) print ">" r[n]; else print; next} {print}' > ref.genome.fa
if [ -s mask.bed ]; then bedtools maskfasta -fi ref.genome.fa -bed mask.bed -fo ref.masked.fa && mv ref.masked.fa ref.genome.fa; fi
cat ref.genome.fa "$PANEL" ${EXTRA:+"$EXTRA"} > ref.fa; rm -f ref.genome.fa ref.genome.fa.fai
strip(){ awk 'BEGIN{FS=OFS="\t"} {o=$1; for(i=2;i<=NF;i++){ if ($i ~ /^cs:Z:/) continue; if ($i ~ /^cg:Z:/ && $i !~ /([3-9][0-9]|[0-9][0-9][0-9]+)I/) continue; o=o OFS $i } print o}'; }
case "$KIND" in
  hifi)  PRESET=map-hifi ;;
  clr)   PRESET=map-pb ;;
  short) PRESET=sr ;;
  *) echo "unknown kind $KIND" >&2; exit 1 ;;
esac
map(){ minimap2 -x $PRESET -c --secondary=no -t 2 -I ${PASS2_INDEX_BASES:-500M} --split-prefix p2tmp -K 200M ref.fa "$1" 2> "$2" | strip; }
if [ "${SAMPLE_ONLY:-0}" != 1 ]; then
  map all_hits.fa pass2.log > pass2.paf.tmp || { echo "pass 2 failed (see pass2.log)" >&2; rm -f ref.fa pass2.paf.tmp p2tmp.*; exit 1; }
  mv pass2.paf.tmp pass2.paf
  echo "pass 2: $(grep -c '>' all_hits.fa) reads, $(wc -l < pass2.paf) alignments"
fi
# the random sample of the whole library kept by the first pass (31c_screen_run.py --sample), aligned
# in the same way: what the library is made of, and how often it joins two pieces of slug DNA
# Long reads: no more than SAMPLE_MAX reads of it (default 12,000) are aligned, every k-th one; k is written
# to sample.thin so that 46_candidates.py knows what share of the library the alignments stand for.
if ls ./*.sample.fa > /dev/null 2>&1; then
  k=1
  if [ "$KIND" != short ]; then n=$(cat ./*.sample.fa | grep -c '>'); k=$(( (n + ${SAMPLE_MAX:-12000} - 1) / ${SAMPLE_MAX:-12000} )); fi
  cat ./*.sample.fa | awk -v k="$k" 'NR % 2 == 1 {keep = (((NR + 1) / 2) % k == 0)} keep' > all_sample.fa
  echo "$k" > sample.thin
  map all_sample.fa sample.log > sample.paf.tmp || { echo "alignment of the sample failed (see sample.log)" >&2; rm -f ref.fa sample.paf.tmp p2tmp.*; exit 1; }
  mv sample.paf.tmp sample.paf
  echo "sample: $(grep -c '>' all_sample.fa) reads, $(wc -l < sample.paf) alignments"
fi
rm -f ref.fa
