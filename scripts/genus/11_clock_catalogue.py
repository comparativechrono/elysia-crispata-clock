#!/usr/bin/env python3
"""Build the clock-gene catalogue of one species.

usage: 11_clock_catalogue.py slug config/species.tsv asm_dir clock_dir refs_dir out.json

Evidence used
  projection   the E. crispata clock proteins aligned to the genome (10_project_clock_genes.sh)
  annotation   the species' own annotated proteins, where GenBank has them (protein.faa.gz and
               genomic.gff.gz in the assembly directory)

For every candidate (a projected locus or an annotated protein):
  1. reciprocal test: searched against all E. crispata proteins plus the 47 catalogue proteins
     (refs_dir/ecr_refs.dmnd). It is kept for gene G only if G is the best-scoring catalogue
     protein and no other E. crispata protein that is not G scores clearly better (> 8 %).
     This separates paralogues (Cry1 / Cry2 / photolyases, Ck1 forms, opsins ...).
  2. scored by a global alignment to the E. crispata protein of G: identical residues, and the
     fraction of the E. crispata protein that is covered.
The model for a gene is the best-scoring candidate; an annotated protein is preferred when it
is within 2 % of the best projection, because it rests on that species' own evidence.

A projected model may be built from several accepted loci that cover different parts of the
protein. If they lie in order on one sequence and strand they are one gene with a very long
intron; otherwise the assembly has the gene in pieces.

status (provisional; 12_compare_species.py re-scores completeness against the genus consensus)
  complete              >= 90 % of the E. crispata protein is covered, no frameshift or stop
  complete_discrepancy  >= 90 % covered, but the genome needs a frameshift or stop codon
  fragmented            >= 90 % covered only by pieces that are out of order or on different sequences
  partial               less than 90 % covered
  absent                no locus passes the reciprocal test
"""
import sys, os, re, json, gzip, subprocess, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from miniprot_parse import parse
from Bio import Align
from Bio.Align import substitution_matrices

slug, cfg, asm_dir, clock_dir, refs, out_json = sys.argv[1:7]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
row = next(l.rstrip("\n").split("\t") for l in open(cfg) if l.startswith(slug + "\t"))
species, tag, acc, asm_name, level = row[1:6]
A = os.path.join(asm_dir, slug); wd = os.path.join(clock_dir, slug + ".work"); os.makedirs(wd, exist_ok=True)

def fasta(path):
    op = gzip.open if path.endswith(".gz") else open
    name, seq, d = None, [], {}
    with op(path, "rt") as fh:
        for l in fh:
            if l[0] == ">":
                if name: d[name] = "".join(seq)
                name, seq = l[1:].split()[0], []
            else: seq.append(l.strip())
    if name: d[name] = "".join(seq)
    return d
def clean(p): return re.sub(r"[^A-Za-z]", "", p or "").upper()
def diamond(q, db, o, k=10):
    subprocess.run([BIN + "/diamond", "blastp", "-q", q, "-d", db, "-o", o, "--sensitive", "-e", "1e-5", "-k", str(k), "--threads", "2", "--quiet",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "qlen", "slen", "bitscore"], check=True)
    return [l.rstrip("\n").split("\t") for l in open(o)]

reps = fasta(os.path.join(refs, "reps.faa"))
equiv = dict(l.rstrip("\n").split("\t")[:2] for l in open(os.path.join(refs, "rep_equiv.tsv")))
aligner = Align.PairwiseAligner(); aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -10; aligner.extend_gap_score = -0.5; aligner.mode = "global"
aligner.end_open_gap_score = 0; aligner.end_extend_gap_score = 0            # end gaps are free
def compare(model, ref):
    """identical residues, fraction of ref covered, identity over aligned columns"""
    if not model: return 0, 0.0, 0.0
    al = aligner.align(ref, model)[0]; ident = cov = 0
    for (r0, r1), (m0, m1) in zip(*al.aligned):
        cov += r1 - r0; ident += sum(1 for a, b in zip(ref[r0:r1], model[m0:m1]) if a == b)
    return ident, cov / len(ref), (ident / cov if cov else 0.0)

# ---- chromosomes
chrom_of = {}; chroms = []
rep_path = os.path.join(A, "assembly_report.txt")
fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(os.path.join(A, "genome.fa.gz.fai"))}
if os.path.exists(rep_path):
    for l in open(rep_path, errors="replace"):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[1] == "assembled-molecule" and f[3] == "Chromosome":
            chrom_of[f[4]] = f[2]; chroms.append(dict(name=f[2], accession=f[4], length=fai.get(f[4], int(f[8]))))

# ---- candidates from the projection
hits = parse(os.path.join(clock_dir, slug + ".mp.gff"))
cand = {}                                   # id -> dict
for q, hs in hits.items():
    for i, h in enumerate(hs):
        pep = clean(h["trans"])
        if len(pep) >= 30: cand[f"{q}__p{i}"] = dict(sym=q, kind="projection", pep=pep, hit=h)

# ---- candidates from the species' own annotation
prot_path = os.path.join(A, "protein.faa.gz"); gff_path = os.path.join(A, "genomic.gff.gz")
annotated = os.path.exists(prot_path) and os.path.getsize(prot_path) > 1000
if annotated:
    prot = fasta(prot_path)
    db = os.path.join(wd, "proteome")
    if not os.path.exists(db + ".dmnd"):
        subprocess.run([BIN + "/diamond", "makedb", "--in", prot_path, "-d", db, "--quiet"], check=True)
    for q, s, pid, ln, ql, sl, bs in diamond(os.path.join(refs, "reps.faa"), db, os.path.join(wd, "reps_vs_proteome.tsv"), k=8):
        cand.setdefault(f"{q}__a__{s}", dict(sym=q, kind="annotation", pep=prot[s].replace("*", ""), protein_id=s))

with open(os.path.join(wd, "cands.faa"), "w") as fo:
    for k, c in cand.items():
        if not c["sym"].startswith("NF_"): fo.write(f">{k}\n{c['pep']}\n")
best_rep = {}; best_other = {}; rep_bits = collections.defaultdict(dict)
for q, s, pid, ln, ql, sl, bs in diamond(os.path.join(wd, "cands.faa"), os.path.join(refs, "ecr_refs"), os.path.join(wd, "cands_vs_refs.tsv"), k=25):
    bs = float(bs)
    if s.startswith("REP|"):
        sym = s[4:]; rep_bits[q][sym] = max(bs, rep_bits[q].get(sym, 0))
        if q not in best_rep or bs > best_rep[q][1]: best_rep[q] = (sym, bs)
    else:
        lab = equiv.get(s)
        if q not in best_other or bs > best_other[q][1]: best_other[q] = (s, bs, lab)
for k, c in cand.items():
    if c["sym"].startswith("NF_"): continue
    Q = c["sym"]; br = best_rep.get(k); bo = best_other.get(k); own = rep_bits[k].get(Q, 0)
    ok = bool(br) and br[0] == Q
    if ok and bo:
        if bo[2] not in (None, Q): ok = bo[1] <= own                       # a different catalogue gene's protein scores higher
        elif bo[2] is None and bo[1] > 1.08 * own: ok = False              # an unrelated E. crispata protein fits clearly better
    c["accepted"] = ok; c["nearest"] = (br[0] if br else None); c["nearest_other"] = bo[2] if bo else None
    if ok: c["n_ident"], c["ref_cov"], c["ref_ident"] = compare(c["pep"], reps[Q])

# ---- GFF of the annotation: exon blocks per protein
gff = collections.defaultdict(list)
want = {c["protein_id"] for c in cand.values() if c["kind"] == "annotation" and c.get("accepted")}
if annotated and want and os.path.exists(gff_path):
    with gzip.open(gff_path, "rt") as fh:
        for l in fh:
            if l[0] == "#": continue
            f = l.split("\t")
            if len(f) > 8 and f[2] == "CDS":
                m = re.search(r"protein_id=([^;\s]+)", f[8])
                if m and m.group(1) in want: gff[m.group(1)].append((f[0], int(f[3]), int(f[4]), f[6]))

def comp(s): return s[::-1].translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))
def cds_seq(chrom, exons, strand):
    regs = [f"{chrom}:{a}-{b}" for a, b in sorted(exons)]
    r = subprocess.run([BIN + "/samtools", "faidx", os.path.join(A, "genome.fa.gz")] + regs, capture_output=True, text=True)
    parts = ["".join(b.split("\n")[1:]) for b in r.stdout.split(">")[1:]]
    s = "".join(parts).upper()
    return comp(s) if strand == "-" else s
def cname(c): return chrom_of.get(c, c)

def piece(h, collinear=None):
    d = dict(chrom=cname(h["chrom"]), accession=h["chrom"], start=h["start"], end=h["end"], strand=h["strand"], qstart=h["qstart"], qend=h["qend"],
             ident=h["ident"], fs=h["fs"], stop=h["stop"], exons=[list(x) for x in sorted(h["cds"])])
    return d
def overlaps(h, g): return h["chrom"] == g["chrom"] and h["start"] < g["end"] and h["end"] > g["start"]

genes = []
for Q in reps:
    L = len(reps[Q])
    acc_p = sorted((c for c in cand.values() if c["sym"] == Q and c["kind"] == "projection" and c.get("accepted")), key=lambda c: -c["n_ident"])
    acc_a = sorted((c for c in cand.values() if c["sym"] == Q and c["kind"] == "annotation" and c.get("accepted")), key=lambda c: -c["n_ident"])
    g = dict(symbol=Q, ref_length=L, candidates_projection=len(acc_p), candidates_annotation=len(acc_a))
    if not acc_p and not acc_a:
        rej = [c for c in cand.values() if c["sym"] == Q and c["kind"] == "projection"]
        near = max(rej, key=lambda c: c["hit"]["score"]) if rej else None
        g.update(status="absent", ref_coverage=0.0, pieces=[], protein="", length=0, model_kind=None, extra_loci=[],
                 nearest=dict(chrom=cname(near["hit"]["chrom"]), start=near["hit"]["start"], end=near["hit"]["end"], coverage=round(near["hit"]["cov"], 3),
                              identity=near["hit"]["ident"], looks_like=near.get("nearest_other") or near.get("nearest")) if near else None)
        genes.append(g); continue
    # ---- projection model: the best locus plus any other accepted loci that add new parts of the protein
    proj = None; extra = []
    if acc_p:
        hs = [acc_p[0]["hit"]]; covered = [(hs[0]["qstart"], hs[0]["qend"])]
        for c in acc_p[1:]:
            h = c["hit"]
            if any(overlaps(h, x) for x in hs): continue
            ov = sum(max(0, min(h["qend"], b) - max(h["qstart"], a)) for a, b in covered); span = h["qend"] - h["qstart"]
            if span - ov >= 0.05 * L and ov < 0.3 * span: hs.append(h); covered.append((h["qstart"], h["qend"]))
            elif h["cov"] >= 0.5: extra.append(dict(chrom=cname(h["chrom"]), start=h["start"], end=h["end"], coverage=round(h["cov"], 3), identity=h["ident"], fs=h["fs"], stop=h["stop"], protein=clean(h["trans"])))
        order = sorted(hs, key=lambda h: h["qstart"])
        collinear = len(hs) == 1 or (all(h["chrom"] == hs[0]["chrom"] and h["strand"] == hs[0]["strand"] for h in hs) and
                                     all((order[i]["start"] < order[i + 1]["start"]) == (hs[0]["strand"] == "+") for i in range(len(order) - 1)))
        pep = "".join(clean(h["trans"]) for h in order)
        n_ident, ref_cov, ref_ident = compare(pep, reps[Q]) if len(hs) > 1 else (acc_p[0]["n_ident"], acc_p[0]["ref_cov"], acc_p[0]["ref_ident"])
        proj = dict(hits=order, collinear=collinear, pep=pep, n_ident=n_ident, ref_cov=ref_cov, ref_ident=ref_ident, fs=sum(h["fs"] for h in hs), stop=sum(h["stop"] for h in hs))
    ba = next((c for c in acc_a if gff.get(c["protein_id"])), None)
    if ba and (not proj or ba["n_ident"] >= 0.98 * proj["n_ident"]):
        ex = gff[ba["protein_id"]]; chrom, strand = ex[0][0], ex[0][3]; exons = sorted((a, b) for c, a, b, s_ in ex if c == chrom)
        pieces = [dict(chrom=cname(chrom), accession=chrom, start=exons[0][0], end=exons[-1][1], strand=strand, qstart=0, qend=len(ba["pep"]), ident=1.0, fs=0, stop=0,
                       exons=[list(x) for x in exons])]
        fs = stop = 0; cds = cds_seq(chrom, exons, strand); kind = "annotation"; pep = ba["pep"]; cov = ba["ref_cov"]; rid = ba["ref_ident"]; nid = ba["n_ident"]; collinear = True
        if proj:
            for h in proj["hits"]:
                if not overlaps(h, dict(chrom=chrom, start=exons[0][0], end=exons[-1][1])) and h["cov"] >= 0.5:
                    extra.append(dict(chrom=cname(h["chrom"]), start=h["start"], end=h["end"], coverage=round(h["cov"], 3), identity=h["ident"], fs=h["fs"], stop=h["stop"], protein=clean(h["trans"])))
    else:
        pieces = [piece(h) for h in proj["hits"]]; fs, stop = proj["fs"], proj["stop"]; kind = "projection"; pep = proj["pep"]
        cov, rid, nid, collinear = proj["ref_cov"], proj["ref_ident"], proj["n_ident"], proj["collinear"]
        cds = "".join(cds_seq(h["chrom"], h["cds"], h["strand"]) for h in proj["hits"])
    if cov >= 0.9 and collinear: status = "complete_discrepancy" if (fs or stop) else "complete"
    elif cov >= 0.9: status = "fragmented"
    else: status = "partial"
    g.update(status=status, ref_coverage=round(cov, 4), ref_identity=round(rid, 4), identical_residues=nid, protein=pep, length=len(pep), cds=cds, model_kind=kind,
             protein_id=ba["protein_id"] if kind == "annotation" else None,
             model_source=(f"GenBank protein {ba['protein_id']}" if kind == "annotation" else f"E. crispata {Q} protein aligned to {asm_name}"),
             frameshifts=fs, stops=stop, pieces=pieces, collinear=collinear, long_intron=bool(kind == "projection" and collinear and len(pieces) > 1), extra_loci=extra,
             projection_ref_coverage=round(proj["ref_cov"], 4) if proj else None, annotation_ref_coverage=round(acc_a[0]["ref_cov"], 4) if acc_a else None,
             annotation_protein=acc_a[0].get("protein_id") if acc_a else None)
    genes.append(g)

# ---- genes not found in E. crispata: does this species have them?
nf = []
seedacc = {}
for l in open(os.path.join(refs, "nf_seeds.tsv")):
    lab, accn, group = l.rstrip("\n").split("\t"); seedacc[lab] = (accn, group)
groups = collections.defaultdict(set)
for lab, (accn, grp) in seedacc.items(): groups[grp].add(accn)
nfc = {k: c for k, c in cand.items() if c["sym"].startswith("NF_") and c["hit"]["cov"] >= 0.25}
found = collections.defaultdict(list)
if nfc:
    with open(os.path.join(wd, "nf.faa"), "w") as fo:
        for k, c in nfc.items(): fo.write(f">{k}\n{c['pep']}\n")
    top = {}
    for q, s, pid, ln, ql, sl, bs in diamond(os.path.join(wd, "nf.faa"), os.path.join(refs, "hsdm"), os.path.join(wd, "nf_vs_hsdm.tsv"), k=1):
        top.setdefault(q, s)
    for k, c in nfc.items():
        lab = c["sym"][3:]; grp = seedacc[lab][1]; s = top.get(k)
        if s and s.split("|")[1] in groups[grp]:
            h = c["hit"]; found[grp].append(dict(seed=lab, chrom=cname(h["chrom"]), start=h["start"], end=h["end"], coverage=round(h["cov"], 3), identity=h["ident"], best_hit=s))
for grp in sorted(groups): nf.append(dict(symbol=grp, found=bool(found[grp]), loci=found[grp]))

res = dict(slug=slug, species=species, tag=tag, annotated=annotated,
           assembly=dict(accession=acc, name=asm_name, level=level, span=sum(fai.values()), sequences=len(fai), chromosomes=chroms),
           genes=genes, not_found=nf)
json.dump(res, open(out_json, "w"), indent=1)
c = collections.Counter(g["status"] for g in genes)
print(f"{slug}: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())) + f" | from annotation {sum(1 for g in genes if g.get('model_kind') == 'annotation')}, projection {sum(1 for g in genes if g.get('model_kind') == 'projection')}"
      + f" | not complete: {', '.join(g['symbol'] + ':' + g['status'] + ('(%d%%)' % round(100 * g['ref_coverage']) if g['status'] != 'absent' else '') for g in genes if g['status'] != 'complete')}"
      + f" | NF found: {[n['symbol'] for n in nf if n['found']]}")
