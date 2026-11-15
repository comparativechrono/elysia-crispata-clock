#!/usr/bin/env python3
"""Genes of non-animal origin in the annotated proteomes of four sacoglossans.

usage: 72_alien_proteomes.py ALIEN_DB.dmnd ALGAE_DIR WORKDIR OUT.json TAG=proteome.fa[.gz][:genomic.gff.gz] ...
   e.g. ELYCR=Elysia_crispata-...-pep.fa.gz ELYCH=protein.faa.gz:genomic.gff.gz ELYMA=... PLAOC=PLAOC.fa.gz

Every protein (for Ensembl annotations the longest isoform of each gene) is compared with DIAMOND
(--mid-sensitive, E <= 1e-5, 25 targets; block size 0.5 to fit in 7 GB of memory) with the database of scripts/21_alien_index.py:
13 animals other than sea slugs (MET), three sacoglossan proteomes (SAC), 12 algal and plant
proteomes (ALG), a cyanobacterium (CYA), and the bacteria, archaea (BAC) and fungi (FUN) of
Swiss-Prot. A protein is listed when its best match outside the animals is at least 80 bits and
better than its best match among the 13 animals; it is a strong candidate when it is at least
1.2 times as good, or when no animal matches at all. Unlike scripts/21_alien_index.py, matches to
the other sacoglossans are not counted as evidence of animal origin: a gene that an ancestor of
these slugs acquired from an alga would be shared by them, and a sacoglossan match would hide it.
Whether the other sacoglossans have the gene (best SAC match, the species' own proteome excluded)
is reported instead. Locations come from the Ensembl headers or from the GenBank annotation (GFF).
This is the first pass: 75_alien_candidates.py compares the proteins listed here again without the
cap of 25 matches, places them on their assemblies and judges the DNA that carries them.
"""
import collections
import gzip
import json
import os
import re
import subprocess
import sys

db, algdir, wd, out = sys.argv[1:5]
specs = [a.split("=", 1) for a in sys.argv[5:]]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
LINEAGE = {"BRYSP": "green alga (Bryopsidales)", "OSTQU": "green alga (Bryopsidales)", "CHLRE": "green alga", "VOLCA": "green alga",
           "CHLVA": "green alga", "OSTTA": "green alga", "MICPC": "green alga", "ARATH": "land plant", "ECTSI": "stramenopile",
           "PHATC": "stramenopile", "CYAME": "red alga", "CHOCR": "red alga", "SYNY3": "cyanobacterium"}
os.makedirs(wd, exist_ok=True)


def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].rstrip(), []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


qfa = os.path.join(wd, "proteomes.fa")
meta, loc = {}, {}
if not os.path.exists(qfa):
    with open(qfa, "w") as fo:
        for tag, spec in specs:
            prot, gff = (spec.split(":", 1) + [None])[:2]
            if tag == "ELYCR":                            # Ensembl: longest isoform of each gene
                best = {}
                for h, s in fasta(prot):
                    g = re.search(r"gene:(\S+)", h).group(1).split(".")[0]
                    if g not in best or len(s) > len(best[g][1]):
                        best[g] = (h, s)
                for g, (h, s) in best.items():
                    fo.write(f">ELYCR|{h.split()[0]}|{g} {h.split()[2]}\n{s}\n")
            else:
                for h, s in fasta(prot):
                    first = h.split()[0]
                    acc = first.split("|")[1] if first.count("|") >= 2 else first
                    desc = h.split(" ", 1)[1].split(" OS=")[0].split(" [")[0] if " " in h else ""
                    fo.write(f">{tag}|{acc} {desc}\n{s.replace('*', '')}\n")
for h, s in fasta(qfa):
    qid = h.split()[0]
    meta[qid] = dict(length=len(s), description=h.split(" ", 1)[1] if " " in h else "")
    if qid.startswith("ELYCR|"):
        c = meta[qid]["description"].split(":")
        loc[qid] = (c[1], int(c[2]), int(c[3]))
# locations of GenBank proteins from their CDS
for tag, spec in specs:
    if ":" not in spec:
        continue
    gff = spec.split(":", 1)[1]
    span = {}
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r"protein_id=([^;\n]+)", f[8])
            if not m:
                continue
            k = f"{tag}|{m.group(1)}"
            a, b = int(f[3]), int(f[4])
            if k in span:
                span[k] = (f[0], min(span[k][1], a), max(span[k][2], b))
            else:
                span[k] = (f[0], a, b)
    loc.update(span)
hits = os.path.join(wd, "proteomes_vs_db.tsv")
if not os.path.exists(hits):
    subprocess.run([BIN + "/diamond", "blastp", "-q", qfa, "-d", db, "-o", hits + ".tmp", "--mid-sensitive", "-e", "1e-5", "-k", "25", "--threads", "2",
                    "--quiet", "--block-size", "0.5", "--index-chunks", "2", "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"],
                   check=True)
    os.rename(hits + ".tmp", hits)
best = collections.defaultdict(dict)
sac_other = collections.defaultdict(lambda: (0.0, ""))
with open(hits) as fh:
    for line in fh:
        q, s, pid, ln, ql, sl, ev, bs = line.rstrip("\n").split("\t")
        g = s.split("|")[0]
        bs = float(bs)
        if g == "SAC":
            if s.split("|")[1] == q.split("|")[0]:
                continue                                      # the species' own proteome
            if bs > sac_other[q][0]:
                sac_other[q] = (bs, s)
            continue
        if g not in best[q] or bs > best[q][g][0]:
            best[q][g] = (bs, s, float(pid), int(ln))
emb = json.load(open(os.path.join(os.path.dirname(os.path.abspath(out)), "embedded_plastid.json"))) if os.path.exists(
    os.path.join(os.path.dirname(os.path.abspath(out)), "embedded_plastid.json")) else {}
SLUG = {"ELYCR": "elysia-crispata", "ELYCH": "elysia-chlorotica", "ELYMA": "elysia-marginata"}
cands = []
summary = collections.defaultdict(lambda: collections.Counter())
for q, d in best.items():
    tag = q.split("|")[0]
    summary[tag]["with_hits"] += 1
    met = d.get("MET", (0,))[0]
    alien = max(((d[g][0], g) for g in ("ALG", "CYA", "BAC", "FUN") if g in d), default=(0, None))
    if alien[0] < 80 or alien[0] <= met:
        continue
    g = alien[1]
    top = d[g]
    proteome = top[1].split("|")[1] if g in ("ALG", "CYA") else ""
    c = dict(query=q, species=tag, length=meta[q]["length"], description=meta[q]["description"], group=g,
             lineage=LINEAGE.get(proteome, {"BAC": "bacterium or archaeon", "FUN": "fungus"}.get(g, "")), top_hit=top[1], top_bits=top[0],
             top_identity=top[2], top_aligned=top[3], best_animal_bits=met, best_animal_hit=d.get("MET", (0, ""))[1],
             best_alga_bits=d.get("ALG", (0,))[0], best_cyanobacterium_bits=d.get("CYA", (0,))[0], best_bacterium_bits=d.get("BAC", (0,))[0],
             best_fungus_bits=d.get("FUN", (0,))[0], other_sacoglossan_bits=sac_other[q][0], other_sacoglossan_hit=sac_other[q][1],
             strong=bool(met == 0 or alien[0] >= 1.2 * met))
    if q in loc:
        ctg, a, b = loc[q]
        c["location"] = f"{ctg}:{a}-{b}"
        frs = emb.get(SLUG.get(tag, ""), [])
        c["on_plastid_fragment"] = next((f"{r['contig']}:{r['start']}-{r['end']}" for r in frs if r["contig"] == ctg and a < r["end"] and b > r["start"]), None)
    cands.append(c)
    summary[tag]["listed"] += 1
    summary[tag]["strong"] += c["strong"]
    summary[tag]["alga_like" if g in ("ALG", "CYA") else "microbe_like"] += 1
for tag, spec in specs:
    summary[tag]["proteins"] = sum(1 for q in meta if q.startswith(tag + "|"))
# descriptions of the algal top hits
want = collections.defaultdict(list)
for c in cands:
    if c["group"] in ("ALG", "CYA"):
        want[c["top_hit"].split("|")[-1]].append(c)
for p in sorted(os.listdir(algdir)):
    if not p.endswith(".fa.gz") or p.startswith("uniprot_sprot"):
        continue
    with gzip.open(os.path.join(algdir, p), "rt") as fh:
        for line in fh:
            if line[0] != ">":
                continue
            first = line[1:].split()[0]
            acc = first.split("|")[1] if first.count("|") >= 2 else first
            if acc in want:
                dsc = line[1:].rstrip().split(" ", 1)[1].split(" OS=")[0].split(" [")[0] if " " in line else ""
                for c in want[acc]:
                    c["top_hit_description"] = dsc
cands.sort(key=lambda c: (c["species"], not c["strong"], -c["top_bits"]))
json.dump(dict(summary={k: dict(v) for k, v in summary.items()}, candidates=cands), open(out, "w"), indent=1)
print(json.dumps({k: dict(v) for k, v in summary.items()}))
