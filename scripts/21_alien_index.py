#!/usr/bin/env python3
"""Screen both E. crispata annotations for proteins that look algal rather than animal.

usage: 21_alien_index.py allprot.fa allprot.tsv algae_dir prot_dir workdir
Every annotated protein (Ensembl annotation of xgElyCris1.1, longest isoform per gene, and the
Purdue-assembly annotation) is searched with DIAMOND against a database of
  MET  13 non-sacoglossan animal proteomes         SAC  3 other sacoglossan proteomes
  ALG  12 algal / plant proteomes (2 Bryopsidales) CYA  Synechocystis
  BAC  Swiss-Prot bacteria + archaea               FUN  Swiss-Prot fungi
A protein is listed when a non-animal protein is its best hit overall (bit score >= 80, better
than every animal including the other sacoglossans). It is a strong candidate when that hit also
beats the best non-sacoglossan animal hit by >= 20 % ("alien index" style test on bit scores).
Each listed protein is then placed: gene models of xgElyCris1.1 are checked against the plastid
fragments found in that assembly (embedded_plastid.json, optional 6th argument).
"""
import sys, os, gzip, glob, subprocess, collections, json
allfa, alltsv, algdir, protdir, wd = sys.argv[1:6]
embedded = sys.argv[6] if len(sys.argv) > 6 else None
os.makedirs(wd, exist_ok=True)
D = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin") + "/diamond"
MET = {"HUMAN", "DROME", "DANPL", "TRICA", "MAGGI", "APLCA", "LOTGI", "BIOGL", "POMCA", "MIZYE", "CAPTE", "NEMVE", "BRAFL"}
SAC = {"ELYCH", "ELYMA", "PLAOC"}

def fa_iter(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq = None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                if name: yield name, "".join(seq)
                name, seq = line[1:].rstrip(), []
            else: seq.append(line.strip())
    if name: yield name, "".join(seq)

db = f"{wd}/alien_db.fa"; q = f"{wd}/alien_queries.fa"
if not os.path.exists(f"{wd}/alien_db.dmnd"):
    ann = {}
    for line in open(alltsv):
        f = line.rstrip("\n").split("\t"); ann[f[0]] = f
    best_iso = {}
    with open(db, "w") as fd:
        for h, s in fa_iter(allfa):
            tag = h.split("|")[0]
            if tag in MET: fd.write(f">MET|{h}\n{s}\n")
            elif tag in SAC: fd.write(f">SAC|{h}\n{s}\n")
            elif tag == "ELYCR":
                g = h.split("|")[2]
                if g not in best_iso or len(s) > len(best_iso[g][1]): best_iso[g] = (h, s)
        for p in sorted(glob.glob(f"{algdir}/*.fa.gz")):
            tag = os.path.basename(p)[:-6]
            if tag.startswith("uniprot_sprot"): continue
            grp = "CYA" if tag == "SYNY3" else "ALG"
            for h, s in fa_iter(p):
                acc = h.split("|")[1] if "|" in h else h.split()[0]
                fd.write(f">{grp}|{tag}|{acc}\n{s}\n")
        for name, grp in (("sprot_bacteria", "BAC"), ("sprot_archaea", "BAC"), ("sprot_fungi", "FUN")):
            for h, s in fa_iter(f"{protdir}/{name}.fa.gz"):
                fd.write(f">{grp}|{h.split('|')[1]}|{h.split('|')[2].split()[0]}\n{s}\n")
    with open(q, "w") as fq:
        for g, (h, s) in best_iso.items(): fq.write(f">{h}\n{s}\n")
        for h, s in fa_iter(allfa):
            if h.startswith("ELYCRP|"): fq.write(f">{h}\n{s}\n")
    subprocess.run([D, "makedb", "--in", db, "-d", f"{wd}/alien_db", "--threads", "2", "--quiet"], check=True)
    os.remove(db)
out = f"{wd}/alien_hits.tsv"
if not os.path.exists(out):
    subprocess.run([D, "blastp", "-q", q, "-d", f"{wd}/alien_db", "-o", out + ".tmp", "--sensitive", "-e", "1e-5", "-k", "50", "--threads", "2", "--quiet",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "evalue", "bitscore"], check=True)
    os.rename(out + ".tmp", out)
best = collections.defaultdict(dict)
for line in open(out):
    qid, s, pid, ln, ql, sl, ev, bs = line.rstrip().split("\t")
    grp = s.split("|")[0]; bs = float(bs)
    if grp not in best[qid] or bs > best[qid][grp][0]: best[qid][grp] = (bs, s, float(pid), int(ln), int(ql))
loc = {}
for line in open(alltsv):
    f = line.rstrip("\n").split("\t")
    if f[0].startswith("ELYCR|") and f[5].startswith("xgElyCris1.1:"):
        c = f[5].split(":"); loc[f[0]] = (c[1], int(c[2]), int(c[3]))
frags = [(r["contig"], r["start"], r["end"]) for r in json.load(open(embedded)).get("sanger", [])] if embedded else []
n_queries = sum(1 for line in open(q) if line[0] == ">")
cands = []
for qid, d in best.items():
    met = d.get("MET", (0,))[0]; sac = d.get("SAC", (0,))[0]
    alien = max(((d[g][0], g) for g in ("ALG", "CYA", "BAC", "FUN") if g in d), default=(0, None))
    if alien[0] >= 80 and alien[0] > max(met, sac):
        g = alien[1]
        c = dict(query=qid, source="Ensembl xgElyCris1.1" if qid.startswith("ELYCR|") else "Purdue", length=d[g][4], top_group=g, top_hit=d[g][1], top_bits=d[g][0],
                 top_pident=d[g][2], top_alnlen=d[g][3], best_animal_bits=met, best_animal_hit=d.get("MET", (0, ""))[1],
                 best_sacoglossan_bits=sac, best_alga_bits=d.get("ALG", (0,))[0], best_alga_hit=d.get("ALG", (0, ""))[1],
                 best_bacteria_bits=d.get("BAC", (0,))[0], best_fungi_bits=d.get("FUN", (0,))[0], best_cyano_bits=d.get("CYA", (0,))[0],
                 strong=bool(alien[0] >= 1.2 * met))
        if qid in loc:
            ctg, a, b = loc[qid]; c["location"] = f"{ctg}:{a}-{b}"
            c["on_plastid_fragment"] = next((f"{fc}:{fa}-{fb}" for fc, fa, fb in frags if fc == ctg and a < fb and b > fa), None)
        cands.append(c)
# descriptions of the top hits, read from the proteome FASTA headers
want = collections.defaultdict(list)
for c in cands:
    if c["top_group"] in ("ALG", "CYA"): want[c["top_hit"].split("|")[-1]].append(c)
for p in sorted(glob.glob(f"{algdir}/*.fa.gz")):
    if os.path.basename(p).startswith("uniprot_sprot"): continue
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line[0] != ">": continue
            acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
            if acc in want:
                desc = line[1:].rstrip().split(" ", 1)[1] if " " in line else ""
                for c in want[acc]: c["top_hit_description"] = desc.split(" OS=")[0].split(" [")[0]
for c in cands:
    if c.get("on_plastid_fragment"): c["verdict"] = "chloroplast gene on one of the plastid fragments that the assembly joined to slug DNA (artefact, see above)"
    elif not c["strong"]: c["verdict"] = "conserved protein: animal homologues score within 20 % of the algal hit"
    else: c["verdict"] = "unexplained"
cands.sort(key=lambda c: (not c["strong"], -c["top_bits"]))
summary = dict(n_queries=n_queries, n_with_hits=len(best), n_listed=len(cands), n_strong=sum(c["strong"] for c in cands),
               n_on_plastid_fragment=sum(1 for c in cands if c.get("on_plastid_fragment")), n_unexplained=sum(1 for c in cands if c["verdict"] == "unexplained"),
               by_source=dict(collections.Counter(c["source"] for c in cands)))
json.dump(dict(summary=summary, candidates=cands), open(f"{wd}/alien_candidates_annotated.json", "w"), indent=1)
print(summary)
for c in cands:
    print(f"{c['query']:46s} len={c['length']:5d} {c['top_group']} {c['top_hit'][:30]:30s} bits={c['top_bits']:.0f} id={c['top_pident']:.0f}% animal={c['best_animal_bits']:.0f} sacogl={c['best_sacoglossan_bits']:.0f} "
          f"| {c.get('location', '-')} | {c.get('top_hit_description', '')[:50]} | {c['verdict'][:40]}")
