#!/usr/bin/env python3
"""Short-read version of 34_read_analysis.py: what the library contains, and read pairs that join
slug DNA to organelle DNA.

usage: 37_pair_analysis.py WORKDIR panel.tsv OUT.json [SCREEN.json]
WORKDIR holds pass2.paf (each mate aligned on its own, minimap2 -x sr -c) and RUN.count files.
A mate is  cp / slugmt / algmt / algrDNA  when >= 80 of its bases align to such a sequence at
>= 97% identity, and  slug  when >= 80 bases align to the assembly at >= 97% with mapping quality
>= 20. A pair is mixed when its two mates are of different kinds; a mate is split when it carries
>= 40 bases of each kind itself. Mixed pairs are grouped by where the slug mate lies (1 kb bins);
a place supported by >= 3 pairs with different positions is listed as recurrent.
"""
import sys, os, json, glob, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sam_segments import paf_segments, best_tiling, relabel_organelle_loci
W, panel_tsv, out = sys.argv[1:4]
SCREEN = sys.argv[4] if len(sys.argv) > 4 else None
desc = {}
for l in list(open(panel_tsv))[1:]:
    f = l.rstrip("\n").split("\t"); desc[f[0]] = f[4]
def taxon(t):
    if t.startswith("asmcp|"): return "plastid sequence in the assembly"
    if t.startswith("kcp|"): return "plastid contigs assembled from these reads"
    d = desc.get(t, "")
    w = d.replace("UNVERIFIED: ", "").replace("UNVERIFIED_ORG: ", "").replace("UNVERIFIED_ASMBLY: ", "").split()
    return " ".join(w[:2]) if len(w) >= 2 else (d or t)
def norm(c): return {"asmcp": "cp", "kcp": "cp", "asmalgmt": "algmt"}.get(c, c)
runs = []
for c in sorted(glob.glob(os.path.join(W, "*.count"))):
    n, b = (int(x) for x in open(c).read().split()); run = os.path.basename(c)[:-6]
    sel = json.load(open(os.path.join(W, run + ".select.json"))) if os.path.exists(os.path.join(W, run + ".select.json")) else {}
    runs.append(dict(run=run, reads=n, bases=b, hit_reads=sel.get("hit_reads"), hot_windows=sel.get("hot_windows", [])[:10]))
tot_reads = sum(r["reads"] for r in runs); tot_bases = sum(r["bases"] for r in runs)
S = paf_segments(os.path.join(W, "pass2.paf"), minlen=40)
relabel_organelle_loci(S, SCREEN)
mate = {}; cls_bp = collections.Counter(); don = collections.defaultdict(lambda: dict(mates=0, bp=0, ids=[], acc=collections.Counter()))
by_taxon = {"algrDNA": collections.defaultdict(list), "algmt": collections.defaultdict(list)}
split = []
for r, ss in S.items():
    for s in ss: s["cls"] = norm(s["cls"])
    conf = best_tiling([s for s in ss if s["aln"] >= 40 and s["ident"] >= 0.97 and (s["cls"] != "slug" or s["mapq"] >= 20)], max_overlap=10)
    big = [s for s in conf if s["qe"] - s["qs"] >= 80]
    if big:
        b = max(big, key=lambda s: s["aln"]); mate[r] = b
        if b["cls"] != "slug": cls_bp[b["cls"]] += b["qe"] - b["qs"]
        if b["cls"] in by_taxon: by_taxon[b["cls"]][taxon(b["t"])].append(b["ident"])
    kinds = {s["cls"] for s in conf}
    if "slug" in kinds and len(kinds) > 1:
        a = [s for s in conf if s["cls"] == "slug"][0]; o = [s for s in conf if s["cls"] != "slug"][0]
        split.append(dict(read=r, other=o["cls"], slug=a["t"], slug_pos=(a["ts"] + a["te"]) // 2, other_t=o["t"], other_pos=(o["ts"] + o["te"]) // 2, slug_len=a["qe"] - a["qs"], other_len=o["qe"] - o["qs"]))
    loose = [s for s in ss if s["cls"] == "cp" and s["aln"] >= 80 and s["ident"] >= 0.80]
    if loose:
        b = max(loose, key=lambda s: s["aln"] * s["ident"]); d = don[taxon(b["t"])]
        d["mates"] += 1; d["bp"] += b["qe"] - b["qs"]; d["ids"].append(b["ident"]); d["acc"][b["t"]] += 1
cp_loose = sum(d["bp"] for d in don.values())
donors = []
for tx, d in sorted(don.items(), key=lambda kv: -kv[1]["bp"])[:15]:
    ids = sorted(d["ids"])
    donors.append(dict(name=tx, reads=d["mates"], bp=d["bp"], fraction=round(d["bp"] / max(1, cp_loose), 5), median_identity=round(ids[len(ids) // 2], 4),
                       identity_q10=round(ids[len(ids) // 10], 4), identity_q90=round(ids[(9 * len(ids)) // 10], 4), accession=d["acc"].most_common(1)[0][0].split("|", 1)[1]))
spots = collections.defaultdict(dict)
for r, b in mate.items():
    sp, i = r.rsplit(".", 1); spots[sp][i] = b
have = collections.Counter(); mixed = collections.defaultdict(list)
for sp, m in spots.items():
    kinds = {b["cls"] for b in m.values()}
    for k in kinds:
        if k != "slug": have[k] += 1
    if len(m) == 2 and "slug" in kinds and len(kinds) == 2:
        a = next(b for b in m.values() if b["cls"] == "slug"); o = next(b for b in m.values() if b["cls"] != "slug")
        mixed[o["cls"]].append(dict(spot=sp, slug=a["t"], slug_pos=(a["ts"] + a["te"]) // 2, other_t=o["t"], other_pos=(o["ts"] + o["te"]) // 2))
rates = {k: dict(pairs_with_mate=have[k], pairs_also_with_slug=len(mixed[k]), rate=round(len(mixed[k]) / have[k], 5) if have[k] else None) for k in have}
recurrent = []
for k, rows in mixed.items():
    bins = collections.defaultdict(list)
    for x in rows: bins[(x["slug"], x["slug_pos"] // 1000)].append(x)
    for (c, b), xs in bins.items():
        distinct = {(x["slug_pos"] // 20, x["other_pos"] // 20) for x in xs}
        if len(distinct) >= 3:
            recurrent.append(dict(kind=k, slug=c, start=b * 1000, pairs=len(xs), distinct_pairs=len(distinct), other_targets=sorted({x["other_t"] for x in xs})[:4],
                                  other_positions=sorted(x["other_pos"] for x in xs)[:20]))
other_taxa = {k: [dict(name=t, mates=len(v), mean_identity=round(sum(v) / len(v), 4)) for t, v in sorted(d.items(), key=lambda kv: -len(kv[1]))[:8]] for k, d in by_taxon.items()}
res = dict(kind="short", runs=runs, reads=tot_reads, bases=tot_bases, mates_examined=len(S), other_taxa=other_taxa,
           classes={c: dict(bp=cls_bp[c], fraction_of_bases=round(cls_bp[c] / max(1, tot_bases), 6), pairs=have[c]) for c in cls_bp},
           plastid_any_identity=dict(bp=cp_loose, fraction_of_bases=round(cp_loose / max(1, tot_bases), 6)), donors=donors,
           chimera=dict(rates=rates, recurrent=sorted(recurrent, key=lambda x: -x["distinct_pairs"])[:60], split_mates=dict(collections.Counter(x["other"] for x in split)),
                        split_examples=split[:40]))
json.dump(res, open(out, "w"), indent=1)
print(f"{tot_reads} mates, {tot_bases / 1e9:.2f} Gb screened; {len(S)} mates examined")
for c in sorted(cls_bp, key=lambda c: -cls_bp[c]): print(f"  {c:9s} bp {cls_bp[c]:12d} ({100 * cls_bp[c] / max(1, tot_bases):.4f}% of bases)  pairs {have[c]}")
for d in donors[:10]: print(f"    {d['name']:42s} {d['reads']:8d} mates {d['bp']:11d} bp {100 * d['fraction']:6.2f}%  identity median {d['median_identity']:.3f} (10-90%: {d['identity_q10']:.3f}-{d['identity_q90']:.3f}) {d['accession']}")
print("  mixed pairs:", {k: (v["pairs_also_with_slug"], v["pairs_with_mate"], v["rate"]) for k, v in rates.items()}, "| split mates:", dict(collections.Counter(x["other"] for x in split)), "| recurrent places:", collections.Counter(x["kind"] for x in recurrent))
for x in sorted(recurrent, key=lambda x: -x["distinct_pairs"])[:10]: print("   ", x["kind"], x["slug"], x["start"], "pairs", x["pairs"], "distinct", x["distinct_pairs"], x["other_targets"][:2])
