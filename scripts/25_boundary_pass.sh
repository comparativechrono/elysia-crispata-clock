#!/usr/bin/env bash
# Direct test of each assembled plastid/slug boundary with ALL reads of the given runs.
#   windows      8 kb of slug flank + the plastid segment as assembled, one per fragment
#   calibration  the slug flanks are full of repeats, which would pull in most of the library.
#                The first 20,000 reads of the first run are mapped to the windows; 100-bp bins
#                of slug flank hit by >= 4 of them (expected for single-copy DNA: ~0.2) are
#                repeats and are masked, as is the plastid segment itself -> windows.step1.fa
#   step 1       approximate mapping (no base-level alignment) of every read to the masked
#                windows; keep reads whose alignment comes within 1.5 kb of a boundary
#   step 2       base-level alignment of only those reads to the unmasked windows
# usage: 25_boundary_pass.sh embedded_plastid.json outdir reads_dir "RUN [RUN ...]" name=fasta ...
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
JSON=$1; OUT=$2; READS=$3; RUNLIST=$4; shift 4
mkdir -p "$OUT"
if [ ! -s "$OUT/windows.fa" ]; then
python3 - "$JSON" "$OUT" "$@" <<'PY'
import sys, json, gzip
js, out = sys.argv[1:3]
fas = dict(a.split("=", 1) for a in sys.argv[3:])
D = json.load(open(js))
with open(f"{out}/windows.fa", "w") as fo, open(f"{out}/windows.tsv", "w") as ft:
    ft.write("window\tassembly\tcontig\tseg_start\tseg_end\tslug_side\twin_start\twin_end\tboundary_in_window\n")
    for asm, fa in fas.items():
        rows = D.get(asm, [])
        want = {r["contig"] for r in rows}
        seqs = {}; name = None; buf = []
        op = gzip.open if fa.endswith(".gz") else open
        with op(fa, "rt") as fh:
            for line in fh:
                if line[0] == ">":
                    if name in want: seqs[name] = "".join(buf)
                    name = line[1:].split()[0]; buf = []
                elif name in want: buf.append(line.strip())
        if name in want: seqs[name] = "".join(buf)
        for r in rows:
            s = seqs[r["contig"]]
            left_open = r["left"]["nearest_break"] in ("gap", "sequence end") and (r["left"]["distance"] or 0) < 50
            right_open = r["right"]["nearest_break"] in ("gap", "sequence end") and (r["right"]["distance"] or 0) < 50
            if left_open == right_open: continue
            if left_open:   # plastid | slug : window = segment + 8 kb to the right
                a, b = r["start"], min(len(s), r["end"] + 8000); side = "right"; bnd = r["end"] - a
            else:           # slug | plastid : window = 8 kb to the left + segment
                a, b = max(0, r["start"] - 8000), r["end"]; side = "left"; bnd = r["start"] - a
            w = s[a:b].upper()
            if side == "right" and "N" in w: b = a + w.index("N"); w = s[a:b].upper()
            if side == "left" and "N" in w: cut = w.rindex("N") + 1; a += cut; bnd -= cut; w = s[a:b].upper()
            name_ = f"{asm}|{r['contig']}|{r['start']}-{r['end']}"
            fo.write(f">{name_}\n{w}\n")
            ft.write(f"{name_}\t{asm}\t{r['contig']}\t{r['start']}\t{r['end']}\t{side}\t{a}\t{b}\t{bnd}\n")
PY
fi
echo "$(grep -c '>' "$OUT/windows.fa") windows"
set -- $RUNLIST; FIRST=$1
if [ ! -s "$OUT/windows.step1.fa" ]; then
  zstd -dc "$READS/$FIRST.fa.zst" 2>/dev/null | head -40000 > "$OUT/calib.fa" || true
  minimap2 -x map-hifi --secondary=no -t 2 "$OUT/windows.fa" "$OUT/calib.fa" 2> /dev/null > "$OUT/calib.paf"
  python3 - "$OUT" <<'PY'
import sys, collections
out = sys.argv[1]
W = {}
rows = [l.rstrip("\n").split("\t") for l in open(f"{out}/windows.tsv")]
for f in rows[1:]: W[f[0]] = (f[5], int(f[8]))
hits = collections.defaultdict(collections.Counter)
for l in open(f"{out}/calib.paf"):
    f = l.split("\t")
    for b in range(int(f[7]) // 100, int(f[8]) // 100 + 1): hits[f[5]][b] += 1
name = None; seqs = {}
for l in open(f"{out}/windows.fa"):
    if l[0] == ">": name = l[1:].strip()
    else: seqs[name] = l.strip()
uniq = {}
with open(f"{out}/windows.step1.fa", "w") as fo:
    for w, s in seqs.items():
        side, b = W[w]; s = list(s); n = 0
        for i in range(0, len(s), 100):
            j = min(len(s), i + 100)
            plastid = (i >= b) if side == "left" else (j <= b)
            if plastid or hits[w][i // 100] >= 4: s[i:j] = "N" * (j - i)
            elif abs((i + j) // 2 - b) <= 2000: n += j - i
        uniq[w] = n
        fo.write(f">{w}\n{''.join(s)}\n")
with open(f"{out}/windows.tsv", "w") as ft:
    ft.write("\t".join(rows[0][:9] + ["unique_flank_bp_within_2kb"]) + "\n")
    for f in rows[1:]: ft.write("\t".join(f[:9] + [str(uniq[f[0]])]) + "\n")
print("windows with >= 500 bp of single-copy flank within 2 kb of the boundary:", sum(1 for v in uniq.values() if v >= 500), "of", len(uniq))
PY
fi
for r in $RUNLIST; do
  [ -s "$OUT/$r.near.fa" ] && continue
  zstd -dc "$READS/$r.fa.zst" | minimap2 -x map-hifi --secondary=no -t 2 "$OUT/windows.step1.fa" - 2> "$OUT/$r.step1.log" \
    | awk 'BEGIN{FS=OFS="\t"} NR==FNR {if (FNR>1) b[$1]=$9; next} ($6 in b) && $11>=300 && $8 <= b[$6]+1500 && $9 >= b[$6]-1500 {print $1}' "$OUT/windows.tsv" - \
    | sort -u > "$OUT/$r.near.names"
  zstd -dc "$READS/$r.fa.zst" | seqkit grep -j 2 -w 0 -f "$OUT/$r.near.names" > "$OUT/$r.near.fa"
  echo "$r: $(wc -l < "$OUT/$r.near.names") reads near a boundary"
done
cat "$OUT"/*.near.fa > "$OUT/near.fa"
minimap2 -x map-hifi -c --secondary=no -t 2 "$OUT/windows.fa" "$OUT/near.fa" 2> "$OUT/step2.log" \
  | awk 'BEGIN{FS=OFS="\t"} $11>=500 {o=$1; for(i=2;i<=NF;i++) if ($i !~ /^c[gs]:Z:/) o=o OFS $i; print o}' > "$OUT/boundary.paf"
echo "$(wc -l < "$OUT/boundary.paf") base-level alignments"
echo BOUNDARY_PASS_DONE
