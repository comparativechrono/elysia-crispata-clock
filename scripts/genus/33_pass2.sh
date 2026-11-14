#!/usr/bin/env bash
# Pass 2 of the read screen: base-level alignment of the reads kept by pass 1 against the
# species' own assembly plus the reference panel.
#   - organelle sequences that are part of the assembly (35_assembly_screen.py, "whole") are
#     renamed CLASS|name so that they count as organelle, not slug
#   - plastid-like stretches inside slug sequences (mask.bed) are hidden, so that a plastid read
#     is not explained away by a plastid fragment that the assembler left in a slug sequence
# The index is built in parts (PASS2_INDEX_BASES, default 500M) to keep memory low.
# usage: 33_pass2.sh WORKDIR GENOME.fa.gz PANEL.fa SCREEN.json hifi|clr|short [EXTRA_REF.fa]
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
W=$(realpath "$1"); G=$(realpath "$2"); PANEL=$(realpath "$3"); SCREEN=$(realpath "$4"); KIND=$5; EXTRA=${6:+$(realpath "$6")}
cd "$W"
cat ./*.hits.fa > all_hits.fa
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
strip(){ awk 'BEGIN{FS=OFS="\t"} {o=$1; for(i=2;i<=NF;i++) if ($i !~ /^c[gs]:Z:/) o=o OFS $i; print o}'; }
case "$KIND" in
  hifi)  minimap2 -x map-hifi -c --secondary=no -t 2 -I ${PASS2_INDEX_BASES:-500M} --split-prefix p2tmp -K 200M ref.fa all_hits.fa 2> pass2.log | strip > pass2.paf.tmp ;;
  clr)   minimap2 -x map-pb   -c --secondary=no -t 2 -I ${PASS2_INDEX_BASES:-500M} --split-prefix p2tmp -K 200M ref.fa all_hits.fa 2> pass2.log | strip > pass2.paf.tmp ;;
  short) minimap2 -x sr       -c --secondary=no -t 2 -I ${PASS2_INDEX_BASES:-500M} --split-prefix p2tmp -K 200M ref.fa all_hits.fa 2> pass2.log | strip > pass2.paf.tmp ;;
esac
mv pass2.paf.tmp pass2.paf; rm -f ref.fa
echo "pass 2: $(grep -c '>' all_hits.fa) reads, $(wc -l < pass2.paf) alignments"
