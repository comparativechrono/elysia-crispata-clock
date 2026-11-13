#!/usr/bin/env python3
"""Choose a representative protein model for every catalogue gene.

Candidates per gene
  ensembl:<protein>        Ensembl annotation of xgElyCris1.1 (every isoform)
  purdue:<protein>         annotation of the Purdue assembly GCA_033675545.1
  proj:<assembly>:<query>  translation of the miniprot alignment of an orthologue
                           (E. chlorotica, P. ocellatus, E. marginata, Aplysia, Biomphalaria)
                           or of an E. crispata anchor onto the sanger / purdue / tampa assembly
Annotated proteins that do not re-align cleanly to their own assembly (identity < 98 %,
frameshift, in-frame stop or < 95 % coverage) are set aside.  Each candidate is scored by its mean DIAMOND bit score against the reciprocal-best-hit
orthologues of four outgroup molluscs that are never used as projection sources
(Pomacea, Lottia, Magallana, Mizuhopecten).  An existing annotation (Ensembl when it scores within 2 % of Purdue) is kept
unless the best projection beats it by at least 5 % and 30 bits; then the projection is
used (Sanger assembly preferred when within 2 % of the top projection).

usage: 06_choose_models.py catalogue.tsv workdir allprot.fa
"""
import sys, os, subprocess, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from miniprot_parse import parse

cat, wd, allfa = sys.argv[1:4]
D = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/diamond"
SCORE_SP = ["POMCA", "LOTGI", "MAGGI", "MIZYE"]

seqs, name = {}, None
for line in open(allfa):
    if line[0] == ">": name = line[1:].strip(); seqs[name] = []
    else: seqs[name].append(line.strip())
seqs = {k: "".join(v) for k, v in seqs.items()}

entries = []
for line in open(cat):
    if line.startswith("#") or line.startswith("symbol\t"): continue
    f = line.rstrip("\n").split("\t"); entries.append((f[0], f[5].split(",")))
orth = collections.defaultdict(dict)
for line in list(open(f"{wd}/orthologues.tsv"))[1:]:
    f = line.rstrip("\n").split("\t")
    orth[f[0]][f[1]] = (f[2], f[7] == "yes")
anchors = collections.defaultdict(list)
for pid in seqs:
    if pid.startswith(("ELYCR|", "ELYCRP|")):
        for sym, anc in entries:
            sp, acc, gene = pid.split("|")
            if (sp == "ELYCR" and gene in anc) or (sp == "ELYCRP" and acc in anc):
                anchors[sym].append(pid)
MP = {a: parse(f"{wd}/mp_{a}.gff") for a in ("sanger", "purdue", "tampa", "alt")}

cands = collections.defaultdict(dict)       # symbol -> candidate id -> (sequence, info)
qc = {}
for sym, anc in entries:
    for pid in anchors[sym]:
        kind = "ensembl" if pid.startswith("ELYCR|") else "purdue"
        # QC: an annotated protein must re-align cleanly to the assembly it was predicted on
        hh = MP["sanger" if kind == "ensembl" else "purdue"].get(pid.replace("|", "_"), [])
        ok = bool(hh) and hh[0]["ident"] >= 0.98 and hh[0]["fs"] == 0 and hh[0]["stop"] == 0 and hh[0]["cov"] >= 0.95
        qc[(sym, f"{kind}:{pid.split('|')[1]}")] = ok
        cands[sym][f"{kind}:{pid.split('|')[1]}"] = (seqs[pid], dict(kind=kind))
    queries = [p.replace("|", "_") for p in anchors[sym]]
    queries += [p.replace("|", "_") for sp, (p, rec) in orth[sym].items()
                if rec and sp in ("ELYCH", "PLAOC", "ELYMA", "APLCA", "BIOGL")]
    for asm in ("sanger", "purdue", "tampa"):
        for q in queries:
            hh = MP[asm].get(q, [])
            if not hh: continue
            h = hh[0]
            t = (h["trans"] or "").rstrip("*")
            if not t or "*" in t or h["cov"] < 0.5 or h["fs"] > 0: continue
            cands[sym][f"proj:{asm}:{q}"] = (t, dict(kind="proj_" + asm, chrom=h["chrom"], start=h["start"],
                                                    end=h["end"], strand=h["strand"], ident=h["ident"], cov=h["cov"]))
# score
with open(f"{wd}/cands.fa", "w") as fo:
    for sym in cands:
        for cid, (s, info) in cands[sym].items(): fo.write(f">{sym}@@{cid}\n{s}\n")
with open(f"{wd}/score_refs.fa", "w") as fo:
    for sym in orth:
        for sp in SCORE_SP:
            if sp in orth[sym] and orth[sym][sp][1]:
                fo.write(f">{sym}@@{sp}\n{seqs[orth[sym][sp][0]]}\n")
subprocess.run([D, "makedb", "--in", f"{wd}/score_refs.fa", "-d", f"{wd}/score_refs", "--quiet"], check=True)
subprocess.run([D, "blastp", "-q", f"{wd}/cands.fa", "-d", f"{wd}/score_refs", "-o", f"{wd}/cands_vs_refs.tsv",
                "--very-sensitive", "-e", "1e-3", "-k", "500", "--threads", "2", "--quiet", "--masking", "0",
                "--outfmt", "6", "qseqid", "sseqid", "bitscore"], check=True)
sc = collections.defaultdict(dict)
for line in open(f"{wd}/cands_vs_refs.tsv"):
    q, s, bs = line.rstrip().split("\t")
    qs, cid = q.split("@@", 1); ss, sp = s.split("@@")
    if qs != ss: continue
    sc[(qs, cid)][sp] = max(sc[(qs, cid)].get(sp, 0), float(bs))
PREF = {"ensembl": 0, "proj_sanger": 1, "purdue": 2, "proj_purdue": 3, "proj_tampa": 4}
with open(f"{wd}/model_choice.tsv", "w") as fo, open(f"{wd}/representatives.fa", "w") as fr:
    fo.write("symbol\tcandidate\tkind\tlength\tscore\tchosen\tqc\n")
    for sym, anc in entries:
        nsp = [sp for sp in SCORE_SP if sp in orth[sym] and orth[sym][sp][1]]
        rows = []
        for cid, (s, info) in cands[sym].items():
            score = sum(sc[(sym, cid)].get(sp, 0) for sp in nsp) / max(1, len(nsp)) if nsp else float(len(s))
            rows.append((score, cid, info["kind"], len(s)))
        if not rows: print("!! no candidates for", sym); continue
        # existing annotations are preferred unless a projection is clearly more complete
        annot = [r for r in rows if r[2] in ("ensembl", "purdue") and qc.get((sym, r[1]), True)]
        proj = [r for r in rows if r[2].startswith("proj_")]
        best_annot = None
        if annot:
            ens = max([r for r in annot if r[2] == "ensembl"], default=None)
            pur = max([r for r in annot if r[2] == "purdue"], default=None)
            best_annot = ens if ens and (not pur or ens[0] >= 0.98 * pur[0]) else pur
        best_proj = None
        if proj:
            topp = max(r[0] for r in proj)
            okp = sorted([r for r in proj if r[0] >= 0.98 * topp], key=lambda r: (PREF[r[2]], -r[0], -r[3]))
            best_proj = okp[0]
        # a projection replaces an annotation only when it is clearly better: >= 5 % and >= 30 bits
        if best_annot and (not best_proj or best_proj[0] - best_annot[0] < max(0.05 * best_proj[0], 30)):
            best = best_annot
        else:
            best = best_proj
        for r in sorted(rows, reverse=True):
            fo.write(f"{sym}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[0]:.1f}\t{'yes' if r[1] == best[1] else ''}\t{'' if qc.get((sym, r[1]), True) else 'fails re-alignment QC'}\n")
        fr.write(f">{sym} {best[1]} len={best[3]} score={best[0]:.1f}\n{cands[sym][best[1]][0]}\n")
        alt = "; ".join(f"{r[1].split(':')[0]}{':'+r[1].split(':')[1] if r[1].startswith('proj') else ''}:{r[1].split(':')[-1][-22:]}({r[3]}aa,{r[0]:.0f})" for r in sorted(rows, reverse=True)[:4])
        print(f"{sym:12s} -> {best[1][:58]:58s} {best[3]:5d} aa  score {best[0]:7.1f} (n ref sp {len(nsp)}) | top: {alt}")
