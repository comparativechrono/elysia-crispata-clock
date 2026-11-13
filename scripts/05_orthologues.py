#!/usr/bin/env python3
"""Reciprocal-best-hit orthologues of each catalogue gene in the comparison proteomes.

usage: 05_orthologues.py catalogue.tsv allprot.fa allprot.tsv workdir
Writes workdir/orthologues.tsv (symbol, species, protein id, bitscore, reciprocal flag)
and workdir/anchors.fa (all E. crispata anchor proteins, all isoforms).
"""
import sys, os, subprocess, collections

cat, allfa, alltsv, wd = sys.argv[1:5]
os.makedirs(wd, exist_ok=True)
D = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/diamond"
SPECIES = ["ELYCH", "PLAOC", "ELYMA", "APLCA", "BIOGL", "POMCA", "LOTGI", "MAGGI", "MIZYE", "CAPTE", "DROME", "HUMAN"]

ann = {}
for line in open(alltsv):
    f = line.rstrip("\n").split("\t"); ann[f[0]] = f
seqs, name = {}, None
for line in open(allfa):
    if line[0] == ">": name = line[1:].strip(); seqs[name] = []
    else: seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}

entries = []
for line in open(cat):
    if line.startswith("#") or line.startswith("symbol\t"): continue
    f = line.rstrip("\n").split("\t")
    entries.append(dict(symbol=f[0], anchors=f[5].split(",")))
anchor_ids = collections.defaultdict(list)   # symbol -> protein ids in allprot
owner = {}
for e in entries:
    for a in e["anchors"]:
        for pid, f in ann.items():
            if (a.startswith("ENSUFQG") and f[1] == "ELYCR" and f[3] == a) or (f[1] == "ELYCRP" and f[2] == a):
                anchor_ids[e["symbol"]].append(pid); owner[pid] = e["symbol"]
with open(f"{wd}/anchors.fa", "w") as fo:
    for pid in owner: fo.write(f">{pid}\n{seqs[pid]}\n")
# forward search
subprocess.run([D, "blastp", "-q", f"{wd}/anchors.fa", "-d", os.path.splitext(allfa)[0], "-o", f"{wd}/anchors_vs_all.tsv",
                "--very-sensitive", "-e", "1e-5", "-k", "400", "--threads", "2", "--quiet",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"], check=True)
best = collections.defaultdict(dict)   # symbol -> species -> (bits, id)
for line in open(f"{wd}/anchors_vs_all.tsv"):
    q, s, pid, alen, qlen, slen, ev, bs = line.rstrip().split("\t")
    sp = s.split("|")[0]
    if sp not in SPECIES: continue
    sym = owner[q]
    if sp not in best[sym] or float(bs) > best[sym][sp][0]:
        best[sym][sp] = (float(bs), s, float(pid))
# reciprocal search against all E. crispata proteins
with open(f"{wd}/ecr_all.fa", "w") as fo:
    for pid, s in seqs.items():
        if pid.startswith(("ELYCR|", "ELYCRP|")): fo.write(f">{pid}\n{s}\n")
subprocess.run([D, "makedb", "--in", f"{wd}/ecr_all.fa", "-d", f"{wd}/ecr_all", "--threads", "2", "--quiet"], check=True)
hits = sorted({v[1] for d in best.values() for v in d.values()})
with open(f"{wd}/fwd_hits.fa", "w") as fo:
    for h in hits: fo.write(f">{h}\n{seqs[h]}\n")
subprocess.run([D, "blastp", "-q", f"{wd}/fwd_hits.fa", "-d", f"{wd}/ecr_all", "-o", f"{wd}/hits_vs_ecr.tsv",
                "--very-sensitive", "-e", "1e-5", "-k", "10", "--threads", "2", "--quiet",
                "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore"], check=True)
back = collections.defaultdict(list)
for line in open(f"{wd}/hits_vs_ecr.tsv"):
    q, s, pid, alen, ev, bs = line.rstrip().split("\t")
    back[q].append((float(bs), s))
with open(f"{wd}/orthologues.tsv", "w") as fo:
    fo.write("symbol\tspecies\tprotein\tgene\tlength\tbitscore\tpident\treciprocal\tdescription\n")
    for e in entries:
        sym = e["symbol"]
        for sp in SPECIES:
            if sp not in best[sym]: continue
            bs, s, pid = best[sym][sp]
            bb = sorted(back.get(s, []), reverse=True)
            # reciprocal if the top E. crispata hit (or one within 5% of its score) is an anchor of this gene
            rec = any(owner.get(x) == sym for b, x in bb if b >= 0.95 * bb[0][0]) if bb else False
            fo.write(f"{sym}\t{sp}\t{s}\t{ann[s][3]}\t{ann[s][4]}\t{bs:.0f}\t{pid:.1f}\t{'yes' if rec else 'no'}\t{ann[s][5]}\n")
print("anchors:", len(owner), "orthologue rows written")
