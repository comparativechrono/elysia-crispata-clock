#!/usr/bin/env bash
# Pass 2 of the read-level screen: base-level alignment of the reads that the bait panel
# did not fully explain (>= 500 bp unaligned, or hits to more than one class).
#   (a) against the Sanger primary + alternate assemblies + the unmasked bait panel
#   (b) reads still carrying >= 500 bp of unexplained sequence go against the Purdue and
#       Tampa assemblies, to recognise slug sequence missing from the Sanger assemblies
# usage: 15_pass2_align.sh sel_dir panel.fa sanger.fa alt.fna.gz purdue.fa tampa.fna.gz
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
SEL=$1; PANEL=$2; SANGER=$3; ALT=$4; PURDUE=$5; TAMPA=$6
cd "$SEL"
[ -s pass2_reads.fa ] || seqkit grep -j 2 -w 0 -f pass2.names selected_reads.fa > pass2_reads.fa
strip(){ awk 'BEGIN{FS=OFS="\t"} {o=$1; for(i=2;i<=NF;i++) if ($i !~ /^c[gs]:Z:/) o=o OFS $i; print o}'; }
if [ ! -s pass2.paf ]; then
  (cat "$SANGER"; zcat "$ALT" | sed 's/ .*//'; cat "$PANEL") > ref_a.fa
  minimap2 -x map-hifi -c --secondary=no -t 2 -I 8G -K 100M ref_a.fa pass2_reads.fa 2> pass2.log | strip > pass2.paf.tmp && mv pass2.paf.tmp pass2.paf
  rm -f ref_a.fa
fi
if [ ! -s pass2_other.paf ]; then
  python3 - <<'PY'
import collections
cov = collections.defaultdict(list); qlen = {}
for line in open("pass2.paf"):
    c = line.split("\t")
    if int(c[10]) < 200: continue
    cov[c[0]].append((int(c[2]), int(c[3]))); qlen[c[0]] = int(c[1])
names = set(l.strip() for l in open("pass2.names"))
res = []
for r in names:
    iv = sorted(cov.get(r, [])); tot = 0; cur = None
    for a, b in iv:
        if cur is None: cur = [a, b]
        elif a <= cur[1]: cur[1] = max(cur[1], b)
        else: tot += cur[1] - cur[0]; cur = [a, b]
    if cur: tot += cur[1] - cur[0]
    if r not in qlen or qlen[r] - tot >= 500: res.append(r)
open("pass2_residual.names", "w").write("\n".join(res) + "\n")
print(len(res), "reads with >= 500 bp still unexplained")
PY
  seqkit grep -j 2 -w 0 -f pass2_residual.names pass2_reads.fa > pass2_residual.fa
  (sed 's/^>/>purdue_/' "$PURDUE"; zcat "$TAMPA" | sed 's/ .*//; s/^>/>tampa_/') > ref_b.fa
  minimap2 -x map-hifi -c --secondary=no -t 2 -I 8G -K 100M ref_b.fa pass2_residual.fa 2> pass2_other.log | strip > pass2_other.paf.tmp && mv pass2_other.paf.tmp pass2_other.paf
  rm -f ref_b.fa
fi
wc -l pass2.paf pass2_other.paf
echo PASS2_DONE
