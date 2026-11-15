#!/usr/bin/env python3
"""Foreign DNA in every sequence of the eleven assemblies: one inventory for all the searches.

usage: 76_foreign_census.py WORKDIR ALIEN_DB.dmnd PANEL.fa PANEL.tsv PANEL_EXTRA.tsv RESULTS_DIR OUT.json[.gz] ASSEMBLY ...

The searches with plastid proteins (28_protein_screen.py with 43_foreign_sequences.py), with
chloroplast regulators (73_regulator_summary.py) and with whole gene sets (75_alien_candidates.py)
see a sequence of a bacterium or other organism only when it carries one of the genes they look
for. Here every sequence of each assembly is examined, by the rules of foreign_dna.py:

1. Reading frames (orfscan.py): the share of each sequence in open reading frames of at least 300
   bases (both strands; standard code, or the Mycoplasma code where the standard code gives less
   than 60% and the Mycoplasma code at least 85%; tandem repeats and other simple sequence left out),
   and the 5-kb windows at either end that are coding (at least 60%, or 85% under the Mycoplasma
   code).
2. Candidates: every sequence of which at least 60% lies in such frames, or whose G+C content
   differs from that of the assembly by 0.08 or more (all its bases; 43_foreign_sequences.py, whose
   value 73 and 75 use, takes the first 20 Mb), unless it is a tandem array (fewer than half
   of its 12-base words distinct: arrays of any period lack stop codons in some frame and are not
   genes); and the run of coding windows at an end of a longer sequence whose remainder (at least
   20 kb) has a median share at or below the 90th percentile of 40 random 20-kb windows of the
   assembly (foreign_dna.background), unless the run is itself a tandem array.
3. Translation of every candidate in all six frames (DIAMOND blastx against the same database and
   with the same settings as 43_foreign_sequences.py), and classification of its coding stretches by
   their best match in the same way: another sea slug (the species' own gene set not counted), an
   animal that is not a sea slug, or an alga, plant or microbe (foreign).
4. Verdict on a whole sequence: foreign_dna.judge. Verdict on the coding run at an end: a foreign
   stretch at the end of a slug sequence when it carries at least three foreign coding stretches and
   at least twice as many foreign ones as ones matched best by sea slugs and other animals together.
5. Organelle DNA: a foreign sequence or end of which at least half aligns (minimap2, asm20) to plastid
   or algal mitochondrial genomes of the panel, or to the Acetabularia chloroplast contigs of the
   Purdue assembly (PANEL_EXTRA.tsv), is kleptoplast or algal organelle DNA (as is the algal
   mitochondrial contig of E. papillosa, too divergent to align), reported apart from the DNA of
   bacteria and other organisms; foreign sequences whose coding stretches match algae (or
   cyanobacteria) more often than other microbes are flagged (mostly_algal_matches).
6. In the four assemblies of E. crispata, whose DNA-level search had the genomes of the four bacteria
   assembled from the reference individual in its panel, the sequences that search classed as
   bacterial (RESULTS_DIR/intermediate/crispata_whole_contigs.tsv) are counted too, whatever their
   genes: they are short, and most carry no gene that matches the database.
7. For each foreign sequence, which of the earlier searches found it (RESULTS_DIR: hgt.json,
   genus/foreign_sequences.json.gz, genus/plastid_regulators.json, genus/alien_genes.json and the
   DNA-level list above).

Intermediate files go to WORKDIR (<assembly>.orf.tsv, .distinct.tsv, .ends.tsv, .bg.json, .cand.fa,
.blastx.tsv, .foreign.fa, .panel.paf); they are reused when present. Besides OUT.json[.gz], a table
of every piece of foreign DNA is written next to it (OUT.tsv).
"""
import collections
import csv
import gzip
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import foreign_dna as fd  # noqa: E402
import orfscan  # noqa: E402

BIN = fd.BIN
W = fd.W
OWN_PROTEOME = {"elysia-chlorotica": "ELYCH", "elysia-marginata": "ELYMA"}
# organelle DNA reported by the earlier searches that is too divergent from the panel to align here:
# a fragment of an algal mitochondrial genome in the reference assembly of E. papillosa (README)
KNOWN_ORGANELLE = {("elysia-papillosa", "CDRNNR020002991.1"): "algal mitochondrial DNA"}
GROUP = {"MET": "met", "SAC": "slug", "ALG": "alga", "CYA": "alga", "BAC": "microbe", "FUN": "microbe"}

wd, db, panel_fa, panel_tsv, extra_tsv, resdir, out = sys.argv[1:8]
assemblies = sys.argv[8:]
os.makedirs(wd, exist_ok=True)

genus_name = {}
sl = os.environ.get("SPECLIST", "/home/claude/data/uniprot/speclist.txt")
if os.path.exists(sl):
    for line in open(sl, encoding="latin-1"):
        m = re.match(r"^(\w+)\s+[ABEVO]\s+\d+: N=(\S+)", line)
        if m:
            genus_name[m.group(1)] = m.group(2)
panel_class = {}
for row in open(panel_tsv):
    f = row.rstrip("\n").split("\t")
    if f[0] != "id":
        panel_class[f[0]] = f[1]
# the chloroplast contigs of the Purdue assembly (config/panel_extra.tsv) join the panel, because no
# complete plastid genome of Acetabularia exists
organelle_ref = os.path.join(wd, "panel_with_extra.fa")
if not os.path.exists(organelle_ref):
    with open(organelle_ref + ".tmp", "w") as fo:
        for line in open(panel_fa):
            fo.write(line)
        for row in csv.DictReader(open(extra_tsv), delimiter="\t"):
            seq = fd.fetch(row["assembly"], f'{row["sequence"]}:{int(row["start"]) + 1}-{row["end"]}')
            fo.write(f'>{row["id"]}\n{seq}\n')
    os.replace(organelle_ref + ".tmp", organelle_ref)
for row in csv.DictReader(open(extra_tsv), delimiter="\t"):
    panel_class[row["id"]] = row["class"]


def best(r):
    s, m4 = float(r["orf_standard"]), float(r["orf_mycoplasma"])
    return (m4, "mycoplasma (TGA = Trp)") if (s < 0.6 and m4 >= 0.85) else (s, "standard")


def stretches_of(hits, offset, own_tag):
    """Coding stretches from blastx hits, merged and classified exactly as in 43_foreign_sequences.py."""
    stretches = []
    for bs, a, b, s, pid in sorted(hits, reverse=True):
        g, tag = s.split("|")[0], s.split("|")[1]
        key = "own" if (g == "SAC" and tag == own_tag) else GROUP[g]
        home = None
        for st in stretches:
            ov = min(b, st["end"]) - max(a, st["start"])
            if ov >= 0.5 * (b - a) or ov >= 0.5 * (st["end"] - st["start"]):
                home = st
                break
        if home is None:
            if any(min(b, st["end"]) - max(a, st["start"]) > 0.2 * (b - a) for st in stretches):
                continue
            home = dict(start=a, end=b, best={}, top=s, top_identity=pid)
            stretches.append(home)
        if bs > home["best"].get(key, 0):
            home["best"][key] = bs
    rows = []
    for st in sorted(stretches, key=lambda x: x["start"]):
        d = st["best"]
        met, slug, other = d.get("met", 0), d.get("slug", 0), max(d.get("alga", 0), d.get("microbe", 0))
        if not (met or slug or other):
            cls = "own_only"
        elif (met or slug) and (met >= other or (slug > other and met >= 0.9 * other)):
            cls = "animal" if met else "slug_only"
        elif d.get("alga", 0) >= d.get("microbe", 0):
            cls = "alga"
        else:
            cls = "microbe"
        rows.append(dict(start=st["start"] + offset, end=st["end"] + offset, cls=cls, best_hit=st["top"],
                         best_identity=st["top_identity"], bits={k: round(v) for k, v in d.items()}))
    bac = [r for r in rows if r["cls"] == "microbe" and r["best_hit"].startswith("BAC|")]
    code = lambda r: r["best_hit"].split("|")[2].split("_")[-1]
    gen = collections.Counter(genus_name.get(code(r), code(r)) for r in bac)
    ids = sorted(r["best_identity"] for r in bac)
    return dict(stretches=dict(collections.Counter(r["cls"] for r in rows)), coding_bp=sum(r["end"] - r["start"] for r in rows),
                nearest_bacteria=gen.most_common(5), bacterial_identity_median=(ids[len(ids) // 2] if ids else None), list=rows)


def earlier_searches():
    """Sequences that the earlier searches judged or listed as foreign, by assembly."""
    found = collections.defaultdict(lambda: collections.defaultdict(set))
    hgt = json.load(open(os.path.join(resdir, "hgt.json")))
    name = {"purdue": "elysia-crispata-purdue", "tampa": "elysia-crispata-tampa", "sanger": "elysia-crispata"}
    for r in hgt["protein_screen"]["foreign_contigs"]:
        found[name.get(r["assembly"], r["assembly"])][r["contig"]].add("plastid-protein search")
    F = json.load(gzip.open(os.path.join(resdir, "genus", "foreign_sequences.json.gz")))
    for asm, rows in F["sequences"].items():
        for r in rows:
            found[asm][r["contig"]].add("plastid-protein search")
    R = json.load(open(os.path.join(resdir, "genus", "plastid_regulators.json")))
    for r in R["foreign_sequences"]:
        if r["verdict"] != "slug DNA":
            found[r["assembly"]][r["sequence"]].add("regulator search")
    A = json.load(open(os.path.join(resdir, "genus", "alien_genes.json")))
    for c in A["candidates"]:
        if c.get("verdict", "").startswith("foreign"):
            found[c["assembly"]][c["sequence"]].add("gene-set search")
    for asm, rows in dna_bacteria.items():
        for n in rows:
            found[asm][n].add("DNA-level search")
    return found


# The DNA-level search of the four E. crispata assemblies (scripts 21 to 23) had the genomes of the four
# bacteria assembled from the reference individual in its panel; the sequences that it classed as bacterial
# (at least half of their length aligned to those genomes) are bacterial here too, whatever their genes.
dna_bacteria = collections.defaultdict(dict)
for r in csv.DictReader(open(os.path.join(resdir, "intermediate", "crispata_whole_contigs.tsv")), delimiter="\t"):
    if r["class"] == "bac":
        dna_bacteria[r["assembly"]][r["sequence"]] = r


def fetch_many(asm, regions, path):
    """Write the regions (name, contig, start, end) of an assembly to a FASTA file."""
    with open(path, "w") as fo:
        for name, contig, a, b in regions:
            seq = fd.fetch(asm, f"{contig}:{a + 1}-{b}").upper()
            fo.write(f">{name}\n{seq}\n")


def blastx(fa, tsv):
    if os.path.exists(tsv) and os.path.getmtime(tsv) > os.path.getmtime(fa):
        return
    if os.path.getsize(fa) == 0:
        open(tsv, "w").close()
        return
    subprocess.run([BIN + "/diamond", "blastx", "-q", fa, "-d", db, "-o", tsv + ".tmp", "--sensitive", "-e", "1e-5", "-F", "15", "--range-culling",
                    "-k", "25", "--threads", os.environ.get("THREADS", "2"), "--quiet", "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                    "qstart", "qend", "evalue", "bitscore"], check=True)
    os.replace(tsv + ".tmp", tsv)


def organelle_cover(fa, paf):
    """Share of each query covered by alignments to plastid or algal mitochondrial genomes of the panel."""
    if not (os.path.exists(paf) and os.path.getmtime(paf) > os.path.getmtime(fa)):
        if os.path.getsize(fa) == 0:
            open(paf, "w").close()
        else:
            with open(paf + ".tmp", "w") as fo:
                subprocess.run([BIN + "/minimap2", "-x", "asm20", "-t", os.environ.get("THREADS", "2"), "--secondary=no", organelle_ref, fa],
                               stdout=fo, stderr=subprocess.DEVNULL, check=True)
            os.replace(paf + ".tmp", paf)
    iv = collections.defaultdict(lambda: collections.defaultdict(list))
    qlen = {}
    for line in open(paf):
        f = line.split("\t")
        cls = panel_class.get(f[5])
        if cls in ("cp", "algmt"):
            iv[f[0]][cls].append((int(f[2]), int(f[3])))
            qlen[f[0]] = int(f[1])
    cover = {}
    for q, d in iv.items():
        res = {}
        for cls, xs in d.items():
            tot, end = 0, -1
            for a, b in sorted(xs):
                if b > end:
                    tot += b - max(a, end)
                    end = b
            res[cls] = round(tot / qlen[q], 3)
        cover[q] = res
    return cover


found = earlier_searches()
result = dict(rules=dict(orf_share=0.6, mycoplasma_code=0.85, gc_offset=0.08, distinct_12mers=0.5, end_rest_bp=20000,
                         end_min_foreign=3, end_foreign_ratio=2.0, organelle_cover=0.5), assemblies={})
for asm in assemblies:
    genome = fd.genome(asm)
    orf = os.path.join(wd, f"{asm}.orf.tsv")
    if not os.path.exists(orf):
        orfscan.scan(genome, orf + ".tmp")
        os.replace(orf + ".tmp", orf)
    rows = {r["sequence"]: r for r in csv.DictReader(open(orf), delimiter="\t")}
    acgt = sum(int(r["acgt"]) for r in rows.values())
    genome_gc = round(sum(float(r["gc"]) * int(r["acgt"]) for r in rows.values()) / acgt, 4)
    bgf = os.path.join(wd, f"{asm}.bg.json")
    if not os.path.exists(bgf):
        json.dump(fd.background(asm), open(bgf, "w"))
    bg = json.load(open(bgf))
    # 2. candidates
    whole = [n for n, r in rows.items() if best(r)[0] >= 0.6 or abs(float(r["gc"]) - genome_gc) >= 0.08]
    dtsv = os.path.join(wd, f"{asm}.distinct.tsv")
    if not os.path.exists(dtsv):
        want, dist = set(whole), {}
        for name, s in orfscan.records(genome):
            if name in want:
                dist[name] = orfscan.distinct(s)
        with open(dtsv + ".tmp", "w") as fo:
            for n in whole:
                fo.write(f"{n}\t{dist[n]:.4f}\n")
        os.replace(dtsv + ".tmp", dtsv)
    dist = {l.split("\t")[0]: float(l.split("\t")[1]) for l in open(dtsv)}
    tandem = [n for n in whole if dist[n] < 0.5]
    whole = [n for n in whole if dist[n] >= 0.5]
    ends = []
    for n, r in rows.items():
        L, k = int(r["length"]), int(r["windows"])
        if n in whole or k == 0:
            continue
        for side in ("head", "tail"):
            run = int(r[f"{side}_run"])
            rest = r[f"{side}_rest_median"]
            if run == 0 or run == k or rest == "" or (k - run) * W < 20000 or float(rest) > bg["p90"]:
                continue
            a, b = (0, min(L, run * W)) if side == "head" else ((k - run) * W, L)
            ends.append((f"{n}|{a}|{side}", n, a, b))
    etsv = os.path.join(wd, f"{asm}.ends.tsv")
    if not os.path.exists(etsv):
        with open(etsv + ".tmp", "w") as fo:
            for name, n, a, b in ends:
                fo.write(f"{name}\t{orfscan.distinct(fd.fetch(asm, f'{n}:{a + 1}-{b}').encode()):.4f}\n")
        os.replace(etsv + ".tmp", etsv)
    edist = {l.split("\t")[0]: float(l.split("\t")[1]) for l in open(etsv)}
    tandem_ends = [e for e in ends if edist[e[0]] < 0.5]
    ends = [e for e in ends if edist[e[0]] >= 0.5]
    fa = os.path.join(wd, f"{asm}.cand.fa")
    if not os.path.exists(fa):
        fetch_many(asm, [(f"{n}|0|whole", n, 0, int(rows[n]["length"])) for n in whole] + ends, fa + ".tmp")
        os.replace(fa + ".tmp", fa)
    tsv = os.path.join(wd, f"{asm}.blastx.tsv")
    blastx(fa, tsv)
    hits = collections.defaultdict(list)
    for line in open(tsv):
        q, s, pid, ln, qs, qe, ev, bs = line.rstrip().split("\t")
        a, b = sorted((int(qs), int(qe)))
        hits[q].append((float(bs), a, b, s, float(pid)))
    own = OWN_PROTEOME.get(asm)
    # 4. verdicts
    seqs, end_rows = [], []
    for n in whole:
        r = rows[n]
        rec = stretches_of(hits.get(f"{n}|0|whole", []), 0, own)
        rec["gc"] = float(r["gc"])
        v = fd.judge(asm, n, rec, [], bg, None, genome_gc)
        seqs.append(dict(sequence=n, length=int(r["length"]), gc=round(float(r["gc"]), 4), orf_share=v["orf_share"], genetic_code=v["genetic_code"],
                         distinct_12mers=round(dist[n], 3), stretches=rec["stretches"], stretch_kinds=v["stretch_kinds"], coding_bp=rec["coding_bp"],
                         nearest_bacteria=rec["nearest_bacteria"], bacterial_identity_median=rec["bacterial_identity_median"],
                         verdict=v["verdict"], top_hits=collections.Counter(x["best_hit"].split("|")[-1] for x in rec["list"] if x["cls"] in ("alga", "microbe")).most_common(8)))
    for name, n, a, b in ends:
        rec = stretches_of(hits.get(name, []), a, own)
        k = fd.kinds(rec)
        foreign = k["foreign"] >= 3 and k["foreign"] >= 2 * (k["sea_slug"] + k["animal"])
        end_rows.append(dict(sequence=n, length=int(rows[n]["length"]), start=a, end=b, side=name.split("|")[2], stretches=rec["stretches"], stretch_kinds=k,
                             nearest_bacteria=rec["nearest_bacteria"], bacterial_identity_median=rec["bacterial_identity_median"],
                             verdict="foreign stretch at the end of a slug sequence" if foreign else "slug DNA"))
    # 5. organelle DNA among the foreign sequences and the foreign stretches at sequence ends
    ffa = os.path.join(wd, f"{asm}.foreign.fa")
    fl = [s for s in seqs if s["verdict"] == "foreign sequence"]
    fe = [e for e in end_rows if e["verdict"] != "slug DNA"]
    if not (os.path.exists(ffa) and os.path.getmtime(ffa) > os.path.getmtime(tsv)):
        fetch_many(asm, [(s["sequence"], s["sequence"], 0, s["length"]) for s in fl] +
                   [(f'{e["sequence"]}|{e["start"]}|{e["side"]}', e["sequence"], e["start"], e["end"]) for e in fe], ffa)
    cover = organelle_cover(ffa, os.path.join(wd, f"{asm}.panel.paf"))
    def kind(c):
        if max(c.values(), default=0) >= 0.5:
            return "plastid DNA" if c.get("cp", 0) >= c.get("algmt", 0) else "algal mitochondrial DNA"
        return "bacterium or other organism"
    for s in seqs:
        c = cover.get(s["sequence"], {})
        s["organelle_cover"] = c
        s["found_by"] = sorted(found[asm].get(s["sequence"], ()))
        if s["verdict"] == "foreign sequence":
            s["kind"] = KNOWN_ORGANELLE.get((asm, s["sequence"])) or kind(c)
            # more coding stretches matched best by algae (or cyanobacteria) than by other microbes
            s["mostly_algal_matches"] = s["stretches"].get("alga", 0) > s["stretches"].get("microbe", 0)
    for e in end_rows:
        e["found_by"] = sorted(found[asm].get(e["sequence"], ()))
        if e["verdict"] != "slug DNA":
            e["organelle_cover"] = cover.get(f'{e["sequence"]}|{e["start"]}|{e["side"]}', {})
            e["kind"] = kind(e["organelle_cover"])
    # 5b. sequences that the DNA-level search classed as bacterial (E. crispata assemblies only)
    for n, w in dna_bacteria.get(asm, {}).items():
        x = next((y for y in seqs if y["sequence"] == n), None)
        if x is None:
            r = rows[n]
            sh, code = best(r)
            x = dict(sequence=n, length=int(r["length"]), gc=round(float(r["gc"]), 4), orf_share=round(sh, 3), genetic_code=code, distinct_12mers=None,
                     stretches={}, stretch_kinds=dict(sea_slug=0, animal=0, foreign=0, unmatched=0), coding_bp=0, nearest_bacteria=[],
                     bacterial_identity_median=None, verdict="not a candidate", top_hits=[], organelle_cover={}, found_by=[])
            seqs.append(x)
        x["dna_match"] = w["best_match"]
        if x["verdict"] != "foreign sequence":
            x["verdict_by_genes"] = x["verdict"]
            x["verdict"] = "foreign sequence"
            x["kind"] = "bacterium or other organism"
            x["mostly_algal_matches"] = False
        x["found_by"] = sorted(set(x["found_by"]) | found[asm].get(n, set()))
    # 6. earlier finds that the census does not count as foreign
    census_foreign = {s["sequence"] for s in seqs if s["verdict"] == "foreign sequence"} | {e["sequence"] for e in end_rows if e["verdict"] != "slug DNA"}
    missed = []
    for n, how in sorted(found[asm].items()):
        if n not in census_foreign:
            r = rows.get(n)
            missed.append(dict(sequence=n, found_by=sorted(how), length=int(r["length"]) if r else None,
                               census=next((s["verdict"] for s in seqs if s["sequence"] == n), "not a candidate"),
                               orf_share=round(best(r)[0], 3) if r else None, gc=float(r["gc"]) if r else None))
    other = [s for s in seqs if s.get("kind") == "bacterium or other organism"]
    org = [s for s in seqs if s.get("kind") in ("plastid DNA", "algal mitochondrial DNA")]
    alg = [s for s in other if s["mostly_algal_matches"]]
    und = [s for s in seqs if s["verdict"].startswith("undetermined")]
    fend = [e for e in end_rows if e.get("kind") == "bacterium or other organism"]
    oend = [e for e in end_rows if e.get("kind") in ("plastid DNA", "algal mitochondrial DNA")]
    result["assemblies"][asm] = dict(
        sequences=len(rows), bases=sum(int(r["length"]) for r in rows.values()), genome_gc=genome_gc,
        background_orf_share={k: v for k, v in bg.items() if k != "list"},
        candidates=dict(whole=len(whole), whole_bp=sum(int(rows[n]["length"]) for n in whole), tandem_arrays=len(tandem),
                        tandem_array_bp=sum(int(rows[n]["length"]) for n in tandem), ends=len(ends), ends_bp=sum(b - a for _, _, a, b in ends),
                        tandem_array_ends=len(tandem_ends)),
        verdicts=dict(collections.Counter(s["verdict"] for s in seqs)),
        foreign=dict(sequences=len(other), bp=sum(s["length"] for s in other), of_them_mostly_algal_matches=len(alg),
                     of_them_mostly_algal_matches_bp=sum(s["length"] for s in alg), organelle_sequences=len(org), organelle_bp=sum(s["length"] for s in org),
                     undetermined=len(und), undetermined_bp=sum(s["length"] for s in und),
                     ends=len(fend), ends_bp=sum(e["end"] - e["start"] for e in fend), organelle_ends=len(oend), organelle_ends_bp=sum(e["end"] - e["start"] for e in oend)),
        foreign_sequences=sorted(other, key=lambda s: -s["length"]), organelle_sequences=sorted(org, key=lambda s: -s["length"]),
        foreign_ends=fend, organelle_ends=oend, end_candidates_slug=len(end_rows) - len(fend) - len(oend),
        other_candidates=[s for s in seqs if s["verdict"] != "foreign sequence"],
        earlier_finds_not_foreign_here=missed)
    a = result["assemblies"][asm]
    print(asm, "genome GC", genome_gc, "| candidates", a["candidates"], "| verdicts", a["verdicts"], "| foreign", a["foreign"], "| earlier finds not foreign here", len(missed), flush=True)
if out.endswith(".gz"):
    with gzip.GzipFile(out, "wb", mtime=0) as fo:
        fo.write(json.dumps(result, indent=1).encode())
else:
    json.dump(result, open(out, "w"), indent=1)
# one row for every piece of foreign DNA, for masking or a quick look
with open(re.sub(r"\.json(\.gz)?$", ".tsv", out), "w") as fo:
    fo.write("assembly\tsequence\tstart\tend\tsequence_length\tkind\tgc\torf_share\tgenetic_code\tforeign_stretches\tsea_slug_stretches\t"
             "animal_stretches\tnearest_bacteria\tbacterial_identity_median\tfound_by_earlier_searches\n")
    for asm, a in result["assemblies"].items():
        rows = [(s, 0, s["length"], s["kind"]) for s in a["foreign_sequences"] + a["organelle_sequences"]]
        rows += [(e, e["start"], e["end"], "bacterium or other organism, at the end of a slug sequence" if e["kind"] == "bacterium or other organism" else e["kind"] + ", at the end of a slug sequence")
                 for e in a["foreign_ends"] + a["organelle_ends"]]
        for r, st, en, kd in rows:
            k = r["stretch_kinds"]
            fo.write(f'{asm}\t{r["sequence"]}\t{st}\t{en}\t{r["length"]}\t{kd}\t{r.get("gc", "")}\t{r.get("orf_share", "")}\t{r.get("genetic_code", "")}\t'
                     f'{k["foreign"]}\t{k["sea_slug"]}\t{k["animal"]}\t{", ".join(f"{g} {n}" for g, n in r["nearest_bacteria"][:3])}\t'
                     f'{r["bacterial_identity_median"] if r["bacterial_identity_median"] is not None else ""}\t{"; ".join(r["found_by"])}\n')
