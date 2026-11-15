#!/usr/bin/env python3
"""What the search for chloroplast regulators found in each assembly, and what kind of DNA carries it.

usage: 73_regulator_summary.py config/plastid_regulators.tsv QUERIES.tsv WORKDIR ALIEN_DB.dmnd OUT.json OUT_MATRIX.tsv
   (also writes OUT_loci.json.gz: every place found, with its code)

WORKDIR is the directory of 71_regulator_screen.sh (protein_screen_raw.json, foreign_sequences.json,
loci.faa, loci_vs_db.tsv). Every place where a regulator aligned is given one code:
  A  alga-like on slug DNA          (best match an algal or plant protein, better than any animal one)
  B  microbe-like on slug DNA       (best match a bacterial, archaeal or fungal protein)
  b  alga- or microbe-like, on a foreign sequence of the assembly (bacterium, phage or other
     contaminant)
  j  alga- or microbe-like, in a foreign stretch at the end of a slug sequence (a join that reads
     would have to confirm)
  P  on kleptoplast DNA (a plastid fragment inside a sequence, or a sequence that is plastid DNA)
  M  a slug gene whose translation matches an animal protein best: a relative of the regulator
     where animals have the family, otherwise a gene that shares a domain with it (a kinase, a
     GTPase) or that it aligned to by chance
  S  a slug gene that matches only other sea slugs, and no algal or bacterial protein
  ?  an alignment that the translated locus does not confirm (it matches no protein at all)
  -  nothing found
Other sea slugs are not counted as animal evidence here: a gene that an ancestor acquired from an
alga would be shared by them, so a locus counts as alga- or microbe-like whenever its best algal,
plant or microbial match is at least 80 bits and better than its best match among the 13 animals
that are not sea slugs.

What kind of DNA carries an alga- or microbe-like locus is judged from the whole sequence by
foreign_dna.py: GC content, the share of it in long open reading frames against random windows of the
same assembly, and how many of its coding stretches match microbes and algae rather than animals or
sea slugs (43_foreign_sequences.py). The locus is then on a foreign sequence (bacterium, phage or
other contaminant), in a foreign stretch at the end of a slug sequence, or on slug DNA; only the last would be a candidate
for a gene transfer into the slug genome.

The alignments that the translated locus did not confirm (no protein at all at E <= 1e-3), when
they are long (at least 120 aligned residues) and not of low complexity (entropy >= 3.3 bits), are
searched again with DIAMOND --ultra-sensitive at E <= 1 against the whole database.
"""
import collections
import csv
import gzip
import json
import os
import subprocess
import sys

cfg_path, q_path, wd, db, out_json, out_tsv = sys.argv[1:7]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
ASM = ["elysia-crispata", "elysia-crispata-alt", "elysia-crispata-purdue", "elysia-crispata-tampa", "elysia-chlorotica",
       "elysia-marginata", "elysia-timida", "elysia-leucolegnote", "elysia-atroviridis", "elysia-papillosa", "elysia-papillosa-alt"]
ORDER = "ABbjPMS?"

rows = [l.rstrip("\n").split("\t") for l in open(cfg_path) if l.strip() and not l.startswith("#")]
cfg = [dict(zip(rows[0], r)) for r in rows[1:]]
reg = {c["name"]: c for c in cfg}
reg["AUREO"] = dict(entry="(Ectocarpus AUR1)", name="AUREO", category="light signalling", animal_family="no",
                    role="aureochrome, the blue-light receptor of stramenopiles")
queries = {r["query"]: r for r in csv.DictReader(open(q_path), delimiter="\t")}
R = json.load(open(os.path.join(wd, "protein_screen_raw.json")))
F = json.load(open(os.path.join(wd, "foreign_sequences.json")))
ctx = {(a, s["contig"]): s for a, L in F["sequences"].items() for s in L}
loci = R["loci"]

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import foreign_dna as fd                                         # noqa: E402

# ---------------------------------------------------------------- code of each locus
def base_code(l):
    met = l["met_bits"] or 0
    alien = max(l["alga_bits"] or 0, l["microbe_bits"] or 0)
    cls = l["class"]
    if cls.startswith("plastid"):
        return "P"
    if cls == "bacterial contig":
        return "b"
    if alien >= 80 and alien > met:
        return "A" if (l["alga_bits"] or 0) >= (l["microbe_bits"] or 0) else "B"
    if cls in ("alga-like", "microbe-like"):                 # below 80 bits: too weak to call either way
        return "?"
    if cls == "animal gene":
        return "M" if met > 0 else "S"
    return "?"


foreign = collections.defaultdict(list)
for l in loci:
    l["code"] = base_code(l)
    if l["code"] in "AB":
        foreign[(l["assembly"], l["contig"])].append(l)

bg = {a: fd.background(a) for a in sorted({a for a, _ in foreign})}
verdicts = {}
for (asm, ctg), ls in sorted(foreign.items()):
    s = ctx[(asm, ctg)]
    st = s["stretches"]
    region = None if s["examined_length"] >= s["contig_length"] else (s["examined_start"], s["examined_start"] + s["examined_length"])
    j = fd.judge(asm, ctg, s, [(l["start"] + l["end"]) // 2 for l in ls], bg[asm], region, F["genome_gc"].get(asm))
    verdict, segment, whole, code = j["verdict"], j["foreign_stretch"], j["orf_share"], j["genetic_code"]
    hits = [x["best_hit"].split("|")[-1] for x in s["list"] if x["cls"] in ("microbe", "alga") and x["best_hit"].split("|")[0] in ("BAC", "ALG", "CYA", "FUN")]
    verdicts[(asm, ctg)] = dict(
        assembly=asm, sequence=ctg, length=s["contig_length"], gc=round(s["gc"], 3), genome_gc=F["genome_gc"].get(asm),
        orf_share=whole, genetic_code=code, background_orf_share={k: v for k, v in bg[asm].items() if k != "list"}, stretches=st, stretch_kinds=j["stretch_kinds"],
        bacterial_identity_median=s["bacterial_identity_median"],
        nearest_bacteria=s["nearest_bacteria"], other_genes=hits[:12], verdict=verdict, foreign_stretch=segment,
        loci=[dict(query=l["query"], regulators=queries[l["query"]]["name"], start=l["start"], end=l["end"], identity=l["identity"],
                   exons=l["exons"], code=l["code"], alga_bits=l["alga_bits"], microbe_bits=l["microbe_bits"], animal_bits=l["met_bits"],
                   sea_slug_bits=l["other_sacoglossan_bits"], best_alga_hit=l["best_alga_hit"]) for l in ls])
    new = {"foreign sequence": "b", "foreign stretch at the end of a slug sequence": "j"}.get(verdict)
    if new:                                                    # weaker alga- or microbe-like loci on the same sequence too
        for l in loci:
            if l["assembly"] == asm and l["contig"] == ctg and (l["code"] in "AB" or l["class"] in ("alga-like", "microbe-like")):
                if new == "b" or (segment and segment["start"] <= l["start"] < segment["end"]):
                    l["code"] = new

# ---------------------------------------------------------------- alignments that nothing confirms
weak = [l for l in loci if l["class"] == "no database match" and l["peptide_entropy"] >= 3.3 and l["aligned_aa"] >= 120]
rd = os.path.join(wd, "recheck")
os.makedirs(rd, exist_ok=True)
fa, tsv = os.path.join(rd, "unconfirmed.faa"), os.path.join(rd, "unconfirmed_vs_db.tsv")
ids = {l["id"] for l in weak}
if not os.path.exists(tsv):
    with open(fa, "w") as fo:
        keep = False
        for line in open(os.path.join(wd, "loci.faa")):
            if line.startswith(">"):
                keep = line[1:].split()[0] in ids
            if keep:
                fo.write(line)
    subprocess.run([BIN + "/diamond", "blastp", "-q", fa, "-d", db, "--ultra-sensitive", "-e", "1", "-k", "10", "-b", "0.4", "-c", "1",
                    "-p", "2", "--quiet", "-o", tsv + ".tmp", "-f", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore"], check=True)
    os.rename(tsv + ".tmp", tsv)
rehit = collections.defaultdict(list)
for line in open(tsv):
    q, s, pid, ln, ev, bs = line.rstrip("\n").split("\t")
    rehit[q].append(dict(hit=s, identity=float(pid), length=int(ln), evalue=float(ev), bits=float(bs)))
recheck = dict(n=len(weak), with_any_match=sum(1 for l in weak if rehit[l["id"]]),
               significant=sum(1 for l in weak if any(h["evalue"] <= 1e-5 for h in rehit[l["id"]])),
               significant_outside_sea_slugs=sum(1 for l in weak if any(h["evalue"] <= 1e-5 and not h["hit"].startswith("SAC|") for h in rehit[l["id"]])),
               loci=[dict(assembly=l["assembly"], sequence=l["contig"], start=l["start"], end=l["end"], regulators=queries[l["query"]]["name"],
                          identity=l["identity"], aligned=l["aligned_aa"], exons=l["exons"],
                          best_match=(sorted(rehit[l["id"]], key=lambda h: -h["bits"]) or [None])[0]) for l in weak])

# ---------------------------------------------------------------- matrix
cells = collections.defaultdict(lambda: collections.defaultdict(set))
for l in loci:
    for n in queries[l["query"]]["name"].split(","):
        cells[n][l["assembly"]].add(l["code"])
names = [c["name"] for c in cfg] + ["AUREO"]
with open(out_tsv, "w") as fo:
    fo.write("# codes: A alga-like on slug DNA; B microbe-like on slug DNA; b on a foreign (bacterial, phage) sequence; j in a foreign stretch at the end\n"
             "# of a slug sequence; P on kleptoplast DNA; M a slug gene matching an animal protein best; S a slug gene matching only other\n"
             "# sea slugs; ? alignment not confirmed by the translated locus; - nothing found (scripts/genus/73_regulator_summary.py)\n")
    fo.write("\t".join(["regulator", "entry", "category", "animal_family"] + ASM) + "\n")
    for n in names:
        r = reg[n]
        fo.write("\t".join([n, r["entry"], r["category"], r["animal_family"]] +
                           ["".join(c for c in ORDER if c in cells[n][a]) or "-" for a in ASM]) + "\n")

per_reg = {}
for n in names:
    found = {a: "".join(c for c in ORDER if c in cells[n][a]) or "-" for a in ASM}
    per_reg[n] = dict(category=reg[n]["category"], animal_family=reg[n]["animal_family"], role=reg[n].get("role", ""),
                      queries=sum(1 for q in queries.values() if n in q["name"].split(",")), assemblies=found,
                      on_slug_dna_alga_or_microbe_like=sum(1 for a in ASM if set(found[a]) & set("AB")))
counts = collections.Counter(l["code"] for l in loci)
summary = dict(
    queries=len(queries), regulators=len(names), assemblies=ASM, loci=len(loci), codes=dict(counts),
    alga_or_microbe_like_on_slug_dna=counts.get("A", 0) + counts.get("B", 0),
    regulators_with_alga_or_microbe_like_copy_on_slug_dna=sorted(n for n in names if per_reg[n]["on_slug_dna_alga_or_microbe_like"]),
    sequences_with_alga_or_microbe_like_loci=collections.Counter(v["verdict"] for v in verdicts.values()),
    unconfirmed_long_alignments={k: recheck[k] for k in ("n", "with_any_match", "significant", "significant_outside_sea_slugs")},
    no_database_match=json.load(open(os.path.join(wd, "protein_screen.json")))["no_match"])
json.dump(dict(summary=summary, regulators=per_reg, foreign_sequences=list(verdicts.values()), unconfirmed=recheck,
               background_orf_share=bg), open(out_json, "w"), indent=1, default=str)
# every place found, with its code (gzip without a time stamp, so that a rerun gives the same file)
KEEP = ("assembly", "contig", "contig_length", "start", "end", "strand", "query", "identity", "query_cover", "aligned_aa", "exons",
        "frameshifts", "stops", "peptide_entropy", "met_bits", "other_sacoglossan_bits", "alga_bits", "microbe_bits", "best_alga_hit",
        "class", "code")
with open(out_json[:-5] + "_loci.json.gz", "wb") as raw_fh, gzip.GzipFile(fileobj=raw_fh, mode="wb", mtime=0, filename="") as gz:
    gz.write(json.dumps([dict({k: l.get(k) for k in KEEP}, regulators=queries[l["query"]]["name"]) for l in loci], indent=0).encode())
print(json.dumps(summary, indent=1, default=str))
