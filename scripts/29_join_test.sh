#!/usr/bin/env bash
# Is a join inside a contig real? Test it with reads from the animal the contig was assembled from.
# For each join (contig, position, which side is unique sequence):
#   step 1  map every read to a 4 kb bait on the unique side of the join (fast, no base alignment);
#           keep reads whose alignment comes within 1 kb of the join
#   step 2  align those reads base by base to a window of 8 kb either side of the join and count
#             through  one alignment runs >= 300 bp past the join on both sides (supports the join)
#             diverge  the alignment on the unique side stops within 150 bp of the join while
#                      >= 200 bp of the read continues: the molecule goes somewhere else there
#           and records how far past the join each read's alignment reaches (reach_past_join)
# usage: 29_join_test.sh outdir contigs.fa "contig:pos:left|right ..." reads_dir "RUN ..." map-ont|map-hifi
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
OUT=$1; FA=$2; JOINS=$3; READS=$4; RUNLIST=$5; PRESET=${6:-map-ont}
mkdir -p "$OUT"; samtools faidx "$FA"
: > "$OUT/baits.fa"; : > "$OUT/windows.fa"; : > "$OUT/joins.tsv"
for j in $JOINS; do
  IFS=: read -r c pos side <<< "$j"
  len=$(awk -v c="$c" '$1==c{print $2}' "$FA.fai")
  ws=$(( pos > 8000 ? pos - 8000 : 0 )); we=$(( pos + 8000 < len ? pos + 8000 : len ))
  if [ "$side" = right ]; then bs=$(( pos + 200 )); be=$(( pos + 4200 < len ? pos + 4200 : len )); else be=$(( pos - 200 )); bs=$(( pos > 4200 ? pos - 4200 : 0 )); fi
  samtools faidx "$FA" "$c:$((bs + 1))-$be" | sed "1s/.*/>$c/" >> "$OUT/baits.fa"
  samtools faidx "$FA" "$c:$((ws + 1))-$we" | sed "1s/.*/>$c/" >> "$OUT/windows.fa"
  printf '%s\t%s\t%s\t%s\t%s\n' "$c" "$pos" "$side" "$ws" "$(( pos - ws ))" >> "$OUT/joins.tsv"
done
# Repeats in a bait would pull in much of the library. The first 20,000 reads of the first run are
# mapped to the baits and 100-bp bins hit by >= 3 of them are masked (single-copy DNA expects < 1).
set -- $RUNLIST; FIRST=$1
if [ ! -s "$OUT/baits.masked.fa" ]; then
  zstd -dc "$READS/$FIRST.fa.zst" 2>/dev/null | head -40000 > "$OUT/calib.fa" || true
  minimap2 -x "$PRESET" --secondary=no -t 2 "$OUT/baits.fa" "$OUT/calib.fa" 2> /dev/null | cut -f6,8,9 > "$OUT/calib.tsv"
  python3 - "$OUT" <<'PY'
import sys, collections
out = sys.argv[1]; hits = collections.defaultdict(collections.Counter)
for l in open(f"{out}/calib.tsv"):
    c, a, b = l.split("\t")
    for i in range(int(a) // 100, int(b) // 100 + 1): hits[c][i] += 1
name = None
with open(f"{out}/baits.masked.fa", "w") as fo:
    for l in open(f"{out}/baits.fa"):
        if l[0] == ">": name = l[1:].strip(); fo.write(l); continue
        s = list(l.strip()); n = 0
        for i in range(0, len(s), 100):
            if hits[name][i // 100] >= 3: s[i:i + 100] = "N" * len(s[i:i + 100]); n += 1
        fo.write("".join(s) + "\n"); print(f"bait {name}: {n} of {len(s) // 100 + 1} bins masked as repeats")
PY
  rm -f "$OUT/calib.fa"
fi
for r in $RUNLIST; do
  [ -e "$OUT/$r.near.done" ] && continue
  zstd -dc "$READS/$r.fa.zst" | minimap2 -x "$PRESET" --secondary=no -t 2 "$OUT/baits.masked.fa" - 2> "$OUT/$r.step1.log" \
    | awk -v s="$JOINS" 'BEGIN{FS=OFS="\t"; n=split(s,a," "); for(i=1;i<=n;i++){split(a[i],b,":"); side[b[1]]=b[3]}} $11>=500 && ((side[$6]=="right" && $8<=1000) || (side[$6]=="left" && $9>=$7-1000)) {print $1}' \
    | sort -u > "$OUT/$r.near.names"
  zstd -dc "$READS/$r.fa.zst" | seqkit grep -j 2 -w 0 -f "$OUT/$r.near.names" > "$OUT/$r.near.fa"
  touch "$OUT/$r.near.done"; echo "$r: $(wc -l < "$OUT/$r.near.names") reads at a join"
done
cat "$OUT"/*.near.fa > "$OUT/near.fa"
minimap2 -x "$PRESET" -c --secondary=no -t 2 "$OUT/windows.fa" "$OUT/near.fa" 2> "$OUT/step2.log" | cut -f1-12 > "$OUT/joins.paf"
python3 - "$OUT" <<'PY'
import sys, json, collections
out = sys.argv[1]
J = {}
for l in open(f"{out}/joins.tsv"):
    c, pos, side, ws, b = l.rstrip("\n").split("\t"); J[c] = dict(contig=c, position=int(pos), unique_side=side, b=int(b))
A = collections.defaultdict(list)
for l in open(f"{out}/joins.paf"):
    f = l.split("\t")
    if int(f[10]) < 500 or int(f[9]) / int(f[10]) < 0.80: continue
    A[f[5]].append(dict(read=f[0], qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]), strand=f[4], ts=int(f[7]), te=int(f[8])))
res = []
for c, j in J.items():
    b = j["b"]; right = j["unique_side"] == "right"
    through, diverge, ends, reach = set(), set(), [], []
    for a in A.get(c, []):
        uniq_bp = (a["te"] - max(a["ts"], b)) if right else (min(a["te"], b) - a["ts"])
        if uniq_bp >= 1000: reach.append((b - a["ts"]) if right else (a["te"] - b))
        if a["ts"] <= b - 300 and a["te"] >= b + 300: through.add(a["read"]); continue
        if uniq_bp < 1000: continue
        near = a["ts"] if right else a["te"]
        over = (a["qs"] if a["strand"] == "+" else a["qlen"] - a["qe"]) if right else ((a["qlen"] - a["qe"]) if a["strand"] == "+" else a["qs"])
        if over >= 200 and abs(near - b) <= 2000: ends.append(near - b)
        if over >= 200 and abs(near - b) <= 150: diverge.add(a["read"])
    diverge -= through
    hist = collections.Counter(round(e / 50) * 50 for e in ends)
    res.append(dict(**j, through_reads=len(through), diverging_reads=len(diverge), reads_on_unique_side=len(reach), reach_past_join=sorted(reach, reverse=True),
                    clip_positions_relative_to_join=dict(sorted(hist.items()))))
    print(f"{c}:{j['position']} unique side {j['unique_side']:5s} | reads {len(reach):3d} | through {len(through):3d} | diverge within 150 bp {len(diverge):3d} | reach past join: {sorted(reach, reverse=True)[:25]} | clipped ends by offset: {dict(sorted(hist.items()))}")
json.dump(res, open(f"{out}/join_test.json", "w"), indent=1)
PY
