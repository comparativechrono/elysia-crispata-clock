#!/usr/bin/env python3
"""Clock-gene status in the further assemblies of a species (config/species_assemblies.tsv).

usage: 14_other_assemblies.py config/species_assemblies.tsv CLOCK_DIR results/genus

For E. crispata the status of every gene in its four assemblies comes from its own pipeline
(results/clock_genes.json). For any other species with a second assembly, the catalogue of that
assembly is first built like that of the reference assembly (10_project_clock_genes.sh and
11_clock_catalogue.py under the name <slug>-<key>; CLOCK_DIR/<slug>-<key>.clock.json). Each of
its proteins is then added to the alignment of the gene (results/genus/alignments/, mafft --add
--keeplength, which leaves the alignment and everything computed from it unchanged), and the
gene is scored like the others in 12_compare_species.py: the share of the conserved alignment
positions that the model covers, and from it the status.
Output: results/genus/<slug>.<key>.clock.json
"""
import sys, os, csv, json, math, subprocess, tempfile, collections
cfg, clock_dir, out = sys.argv[1:4]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
def fasta(path):
    d, k = {}, None
    for l in open(path):
        if l[0] == ">": k = l[1:].strip(); d[k] = []
        elif k: d[k].append(l.strip())
    return {k: "".join(v).upper() for k, v in d.items()}
for row in csv.DictReader(open(cfg), delimiter="\t"):
    slug, key = row["slug"], row["key"]
    src = os.path.join(clock_dir, f"{slug}-{key}.clock.json")
    if slug == "elysia-crispata" or not os.path.exists(src): continue
    c = json.load(open(src)); genes = []
    for g in c["genes"]:
        o = dict(symbol=g["symbol"], status=g["status"], length=g.get("length", 0), frameshifts=g.get("frameshifts", 0), stops=g.get("stops", 0),
                 collinear=g.get("collinear", True), pieces=g.get("pieces", []), consensus_coverage=0.0, identity_to_crispata=None,
                 protein=g.get("protein", ""), cds=g.get("cds", ""), model_kind=g.get("model_kind"), model_source=g.get("model_source", ""), long_intron=g.get("long_intron", False))
        aln_path = os.path.join(out, "alignments", g["symbol"] + ".aln.fa")
        if g["status"] != "absent" and g.get("protein") and os.path.exists(aln_path):
            aln = fasta(aln_path); n = len(aln); L = len(next(iter(aln.values())))
            core = [i for i in range(L) if sum(1 for v in aln.values() if v[i] != "-") >= math.ceil(0.5 * n)]
            with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fo: fo.write(f">NEW\n{g['protein']}\n"); new = fo.name
            r = subprocess.run([BIN + "/mafft", "--quiet", "--add", new, "--keeplength", aln_path], capture_output=True, text=True, check=True)
            os.remove(new)
            seq = "".join(x.strip() for x in r.stdout.split(">NEW\n", 1)[1].split("\n")).upper()
            assert len(seq) == L, (g["symbol"], len(seq), L)
            cg = sum(1 for i in core if seq[i] != "-") / len(core) if core else 0.0
            ref = aln.get("ELYCR")
            both = [(x, y) for x, y in zip(seq, ref) if x != "-" and y != "-"] if ref else []
            o["consensus_coverage"] = round(cg, 4)
            o["identity_to_crispata"] = round(sum(1 for x, y in both if x == y) / len(both), 4) if len(both) >= 30 else None
            if cg >= 0.9 and o["collinear"]: o["status"] = "complete_discrepancy" if (o["frameshifts"] or o["stops"]) else "complete"
            elif cg >= 0.9: o["status"] = "fragmented"
            else: o["status"] = "partial"
        genes.append(o)
    res = dict(slug=slug, key=key, assembly=dict(accession=row["accession"], name=row["asm_name"], level=row["level"], span=c["assembly"]["span"], sequences=c["assembly"]["sequences"]),
               genes=genes, summary=dict(collections.Counter(g["status"] for g in genes)))
    json.dump(res, open(os.path.join(out, f"{slug}.{key}.clock.json"), "w"), indent=1)
    print(slug, key, res["summary"])
