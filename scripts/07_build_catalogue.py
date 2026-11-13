#!/usr/bin/env python3
"""Assemble the clock-gene catalogue (JSON + TSV + FASTA) from all evidence.

usage: 07_build_catalogue.py config_dir work_dir data_dir results_dir
Evidence that is not available yet (read support, RNA-seq, tree placement) is simply left out;
re-run after those steps to fill it in.
"""
import sys, os, json, gzip, re, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from miniprot_parse import parse
from pyfaidx import Fasta
from Bio.Seq import Seq

cfg, wd, data, res = sys.argv[1:5]
os.makedirs(f"{res}/sequences", exist_ok=True)
ASM = ["sanger", "alt", "purdue", "tampa"]
ASM_LABEL = {"sanger": "xgElyCris1.1 primary (GCA_963854125.1)", "alt": "xgElyCris1.1 alternate haplotype (GCA_963854135.1)",
             "purdue": "Purdue ASM3367554v1 (GCA_033675545.1)", "tampa": "Tampa Elysia_crispata_v1 (GCA_984581095.1)"}

def read_fa(path, opener=open):
    d, name = {}, None
    with opener(path, "rt") as fh:
        for line in fh:
            if line[0] == ">":
                name = line[1:].rstrip(); d[name] = []
            else:
                d[name].append(line.strip())
    return {k: "".join(v) for k, v in d.items()}

# ---------- inputs
entries = []
for line in open(f"{cfg}/clock_catalogue.tsv"):
    if line.startswith("#") or line.startswith("symbol\t"): continue
    f = line.rstrip("\n").split("\t")
    entries.append(dict(symbol=f[0], name=f[1], category=f[2], tree=f[3], seeds=f[4].split(","), anchors=f[5].split(","), note=f[6] if len(f) > 6 else ""))
reps = {}
for h, s in read_fa(f"{wd}/representatives.fa").items():
    p = h.split(); reps[p[0]] = dict(source=p[1], seq=s)
choice = collections.defaultdict(list)
for line in list(open(f"{wd}/model_choice.tsv"))[1:]:
    f = line.rstrip("\n").split("\t")
    choice[f[0]].append(dict(candidate=f[1], kind=f[2], length=int(f[3]), score=float(f[4]), chosen=f[5] == "yes", qc=f[6] if len(f) > 6 else ""))
orth = collections.defaultdict(list)
for line in list(open(f"{wd}/orthologues.tsv"))[1:]:
    f = line.rstrip("\n").split("\t")
    orth[f[0]].append(dict(species=f[1], protein=f[2].split("|")[1], gene=f[3], length=int(f[4]), bitscore=float(f[5]), pident=float(f[6]), reciprocal=f[7] == "yes", description=f[8]))
REP = {a: parse(f"{wd}/rep_{a}.gff") for a in ASM}
MPQ = {a: parse(f"{wd}/mp_{a}.gff") for a in ("sanger", "purdue", "tampa")}
genome = {"sanger": Fasta(f"{data}/genome/sanger.fa"), "purdue": Fasta(f"{data}/other_asm/purdue.fa")}

# Ensembl annotation: protein -> transcript/gene, CDS exons, gene spans
ens_pep = {}
for h in read_fa(sorted(p for p in map(lambda x: f"{data}/ensembl/" + x, os.listdir(f"{data}/ensembl")) if p.endswith("pep.fa.gz"))[0], gzip.open):
    pid = h.split()[0].split(".")[0]
    ens_pep[pid] = dict(transcript=re.search(r"transcript:(\S+)", h).group(1).split(".")[0], gene=re.search(r"gene:(\S+)", h).group(1).split(".")[0],
                        symbol=(re.search(r"gene_symbol:(\S+)", h) or [None, ""])[1])
ens_cds = {h.split()[0].split(".")[0]: s for h, s in read_fa([f"{data}/ensembl/" + x for x in os.listdir(f"{data}/ensembl") if x.endswith("cds.fa.gz")][0], gzip.open).items()}
tx_cds = collections.defaultdict(list); tx_info = {}; gene_span = []
with gzip.open([f"{data}/ensembl/" + x for x in os.listdir(f"{data}/ensembl") if x.endswith("genes.gff3.gz")][0], "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        g = line.rstrip("\n").split("\t")
        if g[2] == "CDS":
            tx_cds[re.search(r"Parent=transcript:([^;]+)", g[8]).group(1)].append((int(g[3]), int(g[4])))
        elif g[2] == "mRNA":
            tx_info[re.search(r"ID=transcript:([^;]+)", g[8]).group(1)] = (g[0], int(g[3]), int(g[4]), g[6])
        elif g[2] == "gene":
            gene_span.append((g[0], int(g[3]), int(g[4]), g[6], re.search(r"ID=gene:([^;]+)", g[8]).group(1)))
# Purdue CDS by locus tag
pur_cds = collections.defaultdict(list)
if os.path.exists(f"{wd}/../../rna/purdue_cds.fa"):
    for h, s in read_fa(f"{wd}/../../rna/purdue_cds.fa").items():
        pur_cds[h.split("|")[0]].append(s)
allann = {}
for line in open(f"{wd}/../allprot.tsv"):
    f = line.rstrip("\n").split("\t")
    if f[1] in ("ELYCR", "ELYCRP"): allann[f[2]] = f
# Pfam
dom = collections.defaultdict(list)
if os.path.exists(f"{wd}/reps.pfam.domtbl"):
    for line in open(f"{wd}/reps.pfam.domtbl"):
        if line.startswith("#"): continue
        c = line.split()
        dom[c[0]].append(dict(name=c[3], acc=c[4].split(".")[0], ievalue=float(c[12]), score=float(c[13]), start=int(c[19]), end=int(c[20])))
optional = {}
for key in ("reads_support", "expression", "tree_placement"):
    p = f"{wd}/{key}.json"
    optional[key] = json.load(open(p)) if os.path.exists(p) else {}

def nonoverlap(dd):
    out = []
    for d in sorted(dd, key=lambda d: -d["score"]):
        if all(d["end"] < o["start"] or d["start"] > o["end"] or (min(d["end"], o["end"]) - max(d["start"], o["start"])) < 0.3 * (d["end"] - d["start"]) for o in out):
            out.append(d)
    return sorted(out, key=lambda d: d["start"])

def union_cov(iv, qlen):
    iv = sorted(iv); tot = 0; cur = None
    for a, b in iv:
        if cur is None: cur = [a, b]
        elif a <= cur[1]: cur[1] = max(cur[1], b)
        else: tot += cur[1] - cur[0]; cur = [a, b]
    if cur: tot += cur[1] - cur[0]
    return tot / qlen

def extract_cds(fa, chrom, blocks, strand):
    seq = "".join(str(fa[chrom][s - 1:e]) for s, e in sorted(blocks))
    return str(Seq(seq).reverse_complement()) if strand == "-" else seq

genes = []
for e in entries:
    sym = e["symbol"]
    if sym not in reps:
        continue
    rep = reps[sym]; kind, rid = rep["source"].split(":", 1)[0], rep["source"].split(":")[-1]
    G = dict(symbol=sym, name=e["name"], category=e["category"], note=e["note"], tree=e["tree"], seeds=e["seeds"],
             protein=rep["seq"], length=len(rep["seq"]), model_source=rep["source"], model_kind=kind,
             candidates=choice[sym], orthologues=orth[sym], domains=nonoverlap(dom.get(sym, [])))
    # anchors -> ids
    G["ensembl_genes"] = [a for a in e["anchors"] if a.startswith("ENSUFQG")]
    G["purdue_proteins"] = [dict(accession=a, locus_tag=allann[a][3] if a in allann else "") for a in e["anchors"] if not a.startswith("ENSUFQG")]
    # CDS + description of model
    cds = ""
    if kind == "ensembl":
        cds = ens_cds.get(ens_pep[rid]["transcript"], "")
        G["model_description"] = f"Ensembl protein {rid} (transcript {ens_pep[rid]['transcript']}, gene {ens_pep[rid]['gene']}) annotated on xgElyCris1.1"
    elif kind == "purdue":
        lt = allann[rid][3]
        for s in pur_cds.get(lt, []):
            if str(Seq(s[:len(s) // 3 * 3]).translate()).rstrip("*") == rep["seq"]: cds = s
        G["model_description"] = f"Purdue-assembly protein {rid} (locus {lt}, GCA_033675545.1)"
    else:
        asm, q = rep["source"].split(":")[1], rep["source"].split(":")[2]
        h = MPQ[asm][q][0]
        if asm in genome:
            cds = extract_cds(genome[asm], h["chrom"], h["cds"], h["strand"])
            tr = str(Seq(cds[:len(cds) // 3 * 3]).translate()).rstrip("*")
            if tr != rep["seq"]: cds = ""   # frame bookkeeping differs; leave CDS empty rather than give a wrong one
        src = q.split("_")[0]
        spname = {"BIOGL": "Biomphalaria glabrata", "APLCA": "Aplysia californica", "ELYCH": "Elysia chlorotica", "PLAOC": "Plakobranchus ocellatus", "ELYMA": "Elysia marginata", "ELYCRP": "E. crispata (Purdue)", "ELYCR": "E. crispata (Ensembl)"}.get(src, src)
        G["model_description"] = (f"Translation of the {ASM_LABEL[asm].split(' (')[0]} genome at {h['chrom']}:{h['start']}-{h['end']} ({h['strand']}), "
                                  f"exon structure guided by the {spname} orthologue {q.split('_')[1]} (miniprot); no existing annotation covered the full gene")
    G["cds"] = cds
    # presence in each assembly
    G["assemblies"] = {}
    for a in ASM:
        hh = [h for h in REP[a].get(sym, []) if h["ident"] >= 0.95]
        allh = REP[a].get(sym, [])
        if kind == "ensembl" and a == "sanger":
            chrom, s, en, strand = tx_info[ens_pep[rid]["transcript"]]
            ex = sorted(tx_cds[ens_pep[rid]["transcript"]])
            G["assemblies"][a] = dict(status="complete", coverage=1.0, pieces=[dict(chrom=chrom, start=ex[0][0], end=ex[-1][1], strand=strand, qstart=0, qend=len(rep["seq"]), ident=1.0, fs=0, stop=0, exons=ex)])
            continue
        if not hh:
            bo = allh[0] if allh else None
            st = "absent"
            if bo and bo["ident"] >= 0.85 and bo["cov"] >= 0.5:
                st = "divergent"      # re-labelled "absent" below if that locus belongs to another catalogue gene (a paralogue)
            G["assemblies"][a] = dict(status=st, coverage=0.0, pieces=[],
                                      best_other=(dict(chrom=bo["chrom"], start=bo["start"], end=bo["end"], strand=bo["strand"], ident=round(bo["ident"], 4), cov=round(bo["cov"], 3), fs=bo["fs"], stop=bo["stop"]) if bo else None))
            continue
        qlen = hh[0]["qlen"]
        cov = union_cov([(h["qstart"], h["qend"]) for h in hh], qlen)
        best = max(hh, key=lambda h: h["qend"] - h["qstart"])
        if best["cov"] >= 0.97:
            status = "complete" if best["fs"] == 0 and best["stop"] == 0 else "complete_discrepancy"
            pieces = [best]
        else:
            # keep a minimal set of pieces that add coverage
            pieces, covered = [], []
            for h in sorted(hh, key=lambda h: -(h["qend"] - h["qstart"])):
                before = union_cov(covered, qlen) if covered else 0
                after = union_cov(covered + [(h["qstart"], h["qend"])], qlen)
                if after - before >= 0.03: pieces.append(h); covered.append((h["qstart"], h["qend"]))
            status = "fragmented" if cov >= 0.9 and len(pieces) > 1 else "partial"
        G["assemblies"][a] = dict(status=status, coverage=round(cov, 3),
                                  pieces=[dict(chrom=h["chrom"], start=h["start"], end=h["end"], strand=h["strand"], qstart=h["qstart"], qend=h["qend"],
                                               ident=round(h["ident"], 4), fs=h["fs"], stop=h["stop"], exons=sorted(h["cds"])) for h in pieces])
    # Ensembl genes overlapping the Sanger pieces
    ov = []
    for p in G["assemblies"]["sanger"]["pieces"]:
        for (c, s, en, st, gid) in gene_span:
            if c == p["chrom"] and st == p["strand"] and en >= p["start"] and s <= p["end"] and gid not in ov: ov.append(gid)
    G["ensembl_overlap"] = ov
    for key in optional:
        if sym in optional[key]: G[key] = optional[key][sym]
    genes.append(G)

# a "divergent" locus that is really the locus of another catalogue gene is a paralogue, not a copy of this gene
for G in genes:
    for a in ASM:
        A = G["assemblies"][a]
        if A["status"] != "divergent": continue
        bo = A["best_other"]
        for H in genes:
            if H is G: continue
            for p in H["assemblies"][a]["pieces"]:
                if p["chrom"] == bo["chrom"] and p["end"] >= bo["start"] and p["start"] <= bo["end"]:
                    A["status"] = "absent"; A["best_other"]["paralogue_of"] = H["symbol"]
notfound = []
for line in open(f"{cfg}/clock_not_found.tsv"):
    if line.startswith("#") or line.startswith("symbol\t"): continue
    f = line.rstrip("\n").split("\t"); notfound.append(dict(symbol=f[0], name=f[1], seeds=f[2], interpretation=f[3]))

json.dump(dict(genes=genes, not_found=notfound, assemblies=ASM_LABEL), open(f"{res}/clock_genes.json", "w"), indent=1)
with open(f"{res}/clock_genes.tsv", "w") as fo:
    fo.write("symbol\tname\tcategory\tprotein_length\tmodel_source\tensembl_gene(s)\tpurdue_protein(s)\tsanger_primary_status\tsanger_coverage\tsanger_locus\talt_haplotype_status\tpurdue_status\ttampa_status\tpfam_domains\thuman_orthologue\tdrosophila_orthologue\n")
    for G in genes:
        A = G["assemblies"]; p = A["sanger"]["pieces"]
        locus = "; ".join(f"{x['chrom']}:{x['start']}-{x['end']}({x['strand']})" for x in p) or "-"
        o = {x["species"]: x for x in G["orthologues"]}
        hs = o.get("HUMAN"); dm = o.get("DROME")
        fo.write("\t".join([G["symbol"], G["name"], G["category"], str(G["length"]), G["model_source"], ",".join(G["ensembl_overlap"] or G["ensembl_genes"]) or "-",
                            ",".join(x["accession"] for x in G["purdue_proteins"]) or "-", A["sanger"]["status"], f"{A['sanger']['coverage']:.2f}", locus,
                            A["alt"]["status"], A["purdue"]["status"], A["tampa"]["status"], ",".join(d["name"] for d in G["domains"]) or "-",
                            (hs["gene"] + ("" if hs["reciprocal"] else " (not reciprocal)")) if hs else "-",
                            (dm["gene"] + ("" if dm["reciprocal"] else " (not reciprocal)")) if dm else "-"]) + "\n")
with open(f"{res}/sequences/clock_proteins.faa", "w") as fp, open(f"{res}/sequences/clock_cds.fna", "w") as fc:
    for G in genes:
        fp.write(f">Ecr_{G['symbol']} {G['name']} | {G['model_source']} | {G['length']} aa\n{G['protein']}\n")
        if G["cds"]: fc.write(f">Ecr_{G['symbol']} {G['name']} | {G['model_source']} | {len(G['cds'])} nt\n{G['cds']}\n")
print(len(genes), "genes;", sum(1 for G in genes if G["cds"]), "with CDS")
st = collections.Counter(G["assemblies"]["sanger"]["status"] for G in genes)
print("Sanger primary status:", dict(st))
