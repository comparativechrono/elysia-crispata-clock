#!/usr/bin/env python3
"""Query set for the search for nuclear genes of plants and algae that regulate the chloroplast.

usage: 70_regulator_queries.py config/plastid_regulators.tsv ALGAE_DIR WORKDIR OUT.fa OUT.tsv

The Arabidopsis proteins listed in the configuration (UniProt entry names, from the reviewed
Arabidopsis proteome ALGAE_DIR/ARATH.fa.gz) are joined by their homologues in eleven algal proteomes
of ALGAE_DIR (green: Bryopsis sp. and Ostreobium quekettii, the two siphonous Bryopsidales, and
Chlamydomonas, Volvox, Chlorella, Ostreococcus and Micromonas; stramenopiles: Ectocarpus and
Phaeodactylum, relatives of Vaucheria; red: Cyanidioschyzon and Chondrus), because the slugs feed
on algae, not on land plants. For each Arabidopsis protein and each alga, the best DIAMOND match
(E <= 1e-10, covering at least 40% of the shorter protein) is taken when its own best match among
all Arabidopsis proteins is the query or a member of the same family in the configuration
(reciprocal best hit at the level of the family: SIG1 to SIG6 count as one family, and so on).
Aureochromes, the blue-light receptors of stramenopiles (LOV domain and bZIP), have no Arabidopsis
counterpart and are added from Ectocarpus (AUR1) and its homologues in the two stramenopiles
(E <= 1e-25 over at least 90 residues, the length of the LOV and bZIP region).
Writes the queries (header TAG|name|accession) and a table of them.
"""
import collections
import gzip
import os
import re
import subprocess
import sys

cfg, algdir, wd, out_fa, out_tsv = sys.argv[1:6]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
D = os.path.join(BIN, "diamond")
ALGAE = ["BRYSP", "OSTQU", "CHLRE", "VOLCA", "CHLVA", "OSTTA", "MICPC", "ECTSI", "PHATC", "CYAME", "CHOCR"]
LINEAGE = {"BRYSP": "green alga (Bryopsidales)", "OSTQU": "green alga (Bryopsidales)", "CHLRE": "green alga", "VOLCA": "green alga",
           "CHLVA": "green alga", "OSTTA": "green alga", "MICPC": "green alga", "ECTSI": "stramenopile", "PHATC": "stramenopile",
           "CYAME": "red alga", "CHOCR": "red alga", "ARATH": "land plant"}
# families: Arabidopsis paralogues that count as one in the reciprocal test
FAMILY = {}
for group in (["SIG1", "SIG2", "SIG3", "SIG4", "SIG5", "SIG6"], ["FLN1", "FLN2"], ["FSD2", "FSD3"], ["RPOTm", "RPOTmp", "RPOTp"],
              ["CHLI1", "CHLI2"], ["FC1/GUN6", "FC2"], ["EX1", "EX2"], ["GLK1", "GLK2"], ["PORA", "PORB"], ["PHYA", "PHYB"],
              ["PHOT1", "PHOT2"], ["CRY1", "CRY2"], ["ZTL", "LKP2", "FKF1"], ["STN7", "STN8"], ["PGRL1A", "PGRL1B"], ["CP12-1", "CP12-2"],
              ["TRXf1", "TRXm1", "TRXm2", "TRXx", "TRXy1", "TRXz"], ["ALB3", "ALB4"], ["TOC33", "TOC34"],
              ["TOC159", "TOC132", "TOC120", "TOC90"], ["TOC75-III", "OEP80"], ["FTSZ1", "FTSZ2-1", "FTSZ2-2"],
              ["CCA1", "LHY", "RVE8"], ["TOC1", "PRR3", "PRR5", "PRR7", "PRR9"], ["LNK1", "LNK2"], ["ARC6", "PARC6"], ["PDV1", "PDV2"]):
    for n in group:
        FAMILY[n] = group[0]


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


def acc_of(header):
    first = header.split()[0]
    return first.split("|")[1] if first.count("|") >= 2 else first


os.makedirs(wd, exist_ok=True)
rows = [l.rstrip("\n").split("\t") for l in open(cfg) if l.strip() and not l.startswith("#")]
head, rows = rows[0], rows[1:]
conf = {r[0]: dict(zip(head, r)) for r in rows}
# Arabidopsis
arath = {}
for h, s in fasta(os.path.join(algdir, "ARATH.fa.gz")):
    arath[h.split()[0].split("|")[2]] = (acc_of(h), h, s)
missing = [e for e in conf if e not in arath]
assert not missing, missing
entry_of_acc = {v[0]: k for k, v in arath.items()}
name_of = {e: conf[e]["name"] for e in conf}
fam_of = {e: FAMILY.get(conf[e]["name"], conf[e]["name"]) for e in conf}
qa = os.path.join(wd, "arath_queries.fa")
with open(qa, "w") as fo:
    for e in conf:
        fo.write(f">{e}\n{arath[e][2]}\n")
# databases: the algae, tagged, and the whole Arabidopsis proteome
alg_fa = os.path.join(wd, "algae.fa")
seqs, desc = {}, {}
with open(alg_fa, "w") as fo:
    for t in ALGAE:
        for h, s in fasta(os.path.join(algdir, t + ".fa.gz")):
            k = f"{t}|{acc_of(h)}"
            seqs[k] = s.replace("*", "")
            desc[k] = h.split(" ", 1)[1].split(" OS=")[0] if " " in h else ""
            fo.write(f">{k}\n{seqs[k]}\n")
ara_fa = os.path.join(wd, "arath_all.fa")
with open(ara_fa, "w") as fo:
    for e, (a, h, s) in arath.items():
        fo.write(f">{e}\n{s}\n")
for fa in (alg_fa, ara_fa):
    if not os.path.exists(fa[:-3] + ".dmnd"):
        subprocess.run([D, "makedb", "--in", fa, "-d", fa[:-3], "--quiet"], check=True)


def search(q, db, out, k=25):
    if not os.path.exists(out):
        subprocess.run([D, "blastp", "-q", q, "-d", db, "-o", out, "--more-sensitive", "-e", "1e-5", "-k", str(k), "--threads", "2", "--quiet",
                        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue",
                        "bitscore"], check=True)
    return [l.rstrip("\n").split("\t") for l in open(out)]


fwd = search(qa, alg_fa[:-3], os.path.join(wd, "arath_vs_algae.tsv"))
best = {}
for q, s, pid, ln, ql, sl, qs, qe, ss, se, ev, bs in fwd:
    t = s.split("|")[0]
    cov = (int(qe) - int(qs) + 1) / min(int(ql), int(sl))
    if float(ev) > 1e-10 or cov < 0.40:
        continue
    key = (q, t)
    if key not in best or float(bs) > best[key][1]:
        best[key] = (s, float(bs), float(pid))
cand = sorted({v[0] for v in best.values()})
cq = os.path.join(wd, "algal_candidates.fa")
with open(cq, "w") as fo:
    for c in cand:
        fo.write(f">{c}\n{seqs[c]}\n")
rev = search(cq, ara_fa[:-3], os.path.join(wd, "algae_vs_arath.tsv"), k=5)
back = {}
for q, s, *rest in rev:
    bs = float(rest[-1])
    if q not in back or bs > back[q][1]:
        back[q] = (s, bs)
# keep reciprocal hits at the level of the family
kept = collections.defaultdict(list)
for (q, t), (s, bs, pid) in best.items():
    b = back.get(s)
    if b and b[0] in conf and fam_of[b[0]] == fam_of[q]:
        kept[q].append((t, s, bs, pid))
# aureochromes: Ectocarpus AUR1 and its homologues in the two stramenopiles (LOV and bZIP)
AUR = "ECTSI|D7FMH8"
aq = os.path.join(wd, "aureochrome.fa")
open(aq, "w").write(f">{AUR}\n{seqs[AUR]}\n")
aur = [AUR]
for q, s, pid, ln, ql, sl, qs, qe, ss, se, ev, bs in search(aq, alg_fa[:-3], os.path.join(wd, "aureochrome_vs_algae.tsv"), k=50):
    # the LOV-bZIP region is about a third of AUR1, so homologues are taken by the length of that match
    if s.split("|")[0] in ("ECTSI", "PHATC") and float(ev) <= 1e-25 and int(ln) >= 90 and s not in aur:
        aur.append(s)
# write; an algal protein that is the best match of several Arabidopsis paralogues (one sigma factor
# of an alga for SIG1 to SIG6) is written once, with all their names in the table
recs, by_seq = [], {}


def add(qid, seq, row):
    k = by_seq.get(seq)
    if k is not None:
        names = recs[k][2]["name"].split(",")
        if row["name"] not in names:
            recs[k][2]["name"] = ",".join(names + [row["name"]])
        return
    by_seq[seq] = len(recs)
    recs.append((qid, seq, row))


for e in conf:
    c = conf[e]
    a, h, s_ = arath[e]
    add(f"ARATH|{c['name']}|{a}", s_, dict(name=c["name"], category=c["category"], animal_family=c["animal_family"], source="ARATH",
                                          accession=a, identity="", description=c["role"]))
    for t, sid, bs, pid in sorted(kept.get(e, [])):
        acc = sid.split("|", 1)[1]
        add(f"{t}|{c['name']}|{acc}", seqs[sid], dict(name=c["name"], category=c["category"], animal_family=c["animal_family"], source=t,
                                                       accession=acc, identity=f"{pid:.1f}", description=desc[sid]))
for sid in aur:
    t, acc = sid.split("|", 1)
    add(f"{t}|AUREO|{acc}", seqs[sid], dict(name="AUREO", category="light signalling", animal_family="no", source=t, accession=acc,
                                            identity="", description="aureochrome (LOV-bZIP blue-light receptor of stramenopiles): " + desc[sid]))
n = collections.Counter()
with open(out_fa, "w") as fo, open(out_tsv, "w") as ft:
    ft.write("query\tname\tcategory\tanimal_family\tsource\tlineage\taccession\tlength\tidentity_to_arabidopsis\tdescription\n")
    for qid, seq, r in recs:
        fo.write(f">{qid} {r['name']} {r['category']}\n{seq}\n")
        ft.write(f"{qid}\t{r['name']}\t{r['category']}\t{r['animal_family']}\t{r['source']}\t{LINEAGE[r['source']]}\t{r['accession']}\t{len(seq)}\t"
                 f"{r['identity']}\t{r['description']}\n")
        n[r["source"]] += 1
print("queries:", sum(n.values()), dict(n))
by_name = collections.Counter()
for e in conf:
    by_name[conf[e]["name"]] = len(kept.get(e, []))
print("Arabidopsis proteins without an algal homologue:", [k for k, v in by_name.items() if v == 0])
