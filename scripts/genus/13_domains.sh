#!/usr/bin/env bash
# Pfam domains of every species' clock-protein models (hmmsearch, gathering thresholds).
# usage: 13_domains.sh results/genus Pfam-A.hmm.gz workdir
set -euo pipefail
export PATH=${BIO_BIN:-/home/claude/tools/mamba/envs/bio/bin}:$PATH
R=$1; PFAM=$2; W=$3; mkdir -p "$W"
python3 - "$R" "$W/models.faa" <<'PY'
import sys, json, glob
n = 0
with open(sys.argv[2], "w") as fo:
    for f in sorted(glob.glob(sys.argv[1] + "/elysia-*.clock.json")):
        c = json.load(open(f))
        for g in c["genes"]:
            if g.get("protein"): fo.write(f">{c['slug']}|{g['symbol']}\n{g['protein']}\n"); n += 1
print(n, "protein models")
PY
[ -s "$W/Pfam-A.hmm" ] || zcat "$PFAM" > "$W/Pfam-A.hmm"
hmmsearch --cut_ga --cpu 2 --domtblout "$W/models.domtbl" -o /dev/null "$W/Pfam-A.hmm" "$W/models.faa"
python3 - "$R" "$W/models.domtbl" <<'PY'
import sys, json, glob, collections
dom = collections.defaultdict(list)
for l in open(sys.argv[2]):
    if l[0] == "#": continue
    f = l.split()
    dom[f[0]].append(dict(name=f[3], acc=f[4].split(".")[0], ievalue=float(f[12]), score=float(f[13]), start=int(f[19]), end=int(f[20])))
for f in sorted(glob.glob(sys.argv[1] + "/elysia-*.clock.json")):
    c = json.load(open(f)); n = 0
    for g in c["genes"]:
        keep = []
        for d in sorted(dom.get(f"{c['slug']}|{g['symbol']}", []), key=lambda d: -d["score"]):
            if all(min(d["end"], k["end"]) - max(d["start"], k["start"]) < 0.5 * (d["end"] - d["start"] + 1) for k in keep): keep.append(d)
        g["domains"] = sorted(keep, key=lambda d: d["start"]); n += len(keep)
    json.dump(c, open(f, "w"), indent=1); print(c["slug"], n, "domains")
PY
rm -f "$W/Pfam-A.hmm"; echo DOMAINS_DONE
