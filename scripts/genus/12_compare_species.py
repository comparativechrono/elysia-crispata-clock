#!/usr/bin/env python3
"""Compare the clock proteins across the Elysia species and settle each gene's status.

usage: 12_compare_species.py config/species.tsv results/clock_genes.json clock_dir allprot.fa outdir

For every catalogue gene the models of all species are aligned (MAFFT) together with the
orthologues of three other gastropods already identified for the E. crispata catalogue
(Plakobranchus ocellatus, Aplysia californica, Biomphalaria glabrata; reciprocal best hits).

completeness   Projected models can miss ends that are poorly conserved, and the E. crispata
               model is not a gold standard. So completeness is measured against the consensus:
               the alignment columns in which at least half of all sequences have a residue.
               A model is complete when it fills >= 90 % of those columns.
identity       pairwise, over the columns where both proteins have a residue.
species tree   genes complete in every Elysia species and present in the outgroups are
               concatenated (consensus columns only) and analysed with IQ-TREE (LG+F+G4,
               1,000 ultrafast bootstraps), rooted on Aplysia.
Writes outdir/<slug>.clock.json (final), alignments/<gene>.aln.fa, genus.json, clock_matrix.tsv,
species_tree.nwk.
"""
import sys, os, json, subprocess, collections, math
cfg, crisp_json, clock_dir, allprot, out = sys.argv[1:6]
BIN = os.environ.get("BIO_BIN", "/home/claude/tools/mamba/envs/bio/bin")
os.makedirs(os.path.join(out, "alignments"), exist_ok=True)
SP = [l.rstrip("\n").split("\t") for l in list(open(cfg))[1:]]
slugs = [r[0] for r in SP]; TAG = {r[0]: r[2] for r in SP}; NAME = {r[0]: r[1] for r in SP}
OUTG = ["PLAOC", "APLCA", "BIOGL"]
crisp = json.load(open(crisp_json)); CG = {g["symbol"]: g for g in crisp["genes"]}; symbols = [g["symbol"] for g in crisp["genes"]]
cat = {}
for s in slugs:
    p = os.path.join(clock_dir, s + ".clock.json")
    if s != "elysia-crispata" and os.path.exists(p): cat[s] = json.load(open(p))
# outgroup orthologue sequences
want = {}
for g in crisp["genes"]:
    for o in g["orthologues"]:
        if o["species"] in OUTG and o.get("reciprocal") in (True, "yes"): want[(o["species"], o["protein"])] = g["symbol"]
oseq = {}; name = None
with open(allprot) as fh:
    for l in fh:
        if l[0] == ">":
            p = l[1:].split()[0].split("|"); name = (p[0], p[1]) if (p[0], p[1]) in want else None
            if name: oseq[name] = []
        elif name: oseq[name].append(l.strip())
oseq = {k: "".join(v) for k, v in oseq.items()}

def mafft(seqs, path):
    with open(path + ".in", "w") as fo:
        for k, v in seqs.items(): fo.write(f">{k}\n{v}\n")
    r = subprocess.run([BIN + "/mafft", "--auto", "--quiet", "--thread", "2", path + ".in"], capture_output=True, text=True, check=True)
    aln = {}; k = None
    for l in r.stdout.split("\n"):
        if l.startswith(">"): k = l[1:].strip(); aln[k] = []
        elif k: aln[k].append(l.strip())
    os.remove(path + ".in")
    aln = {k: "".join(v).upper() for k, v in aln.items()}
    with open(path, "w") as fo:
        for k, v in aln.items(): fo.write(f">{k}\n{v}\n")
    return aln

G = {}; concat = collections.defaultdict(list); used = []
for sym in symbols:
    seqs = {"ELYCR": CG[sym]["protein"]}
    for s, c in cat.items():
        g = next(x for x in c["genes"] if x["symbol"] == sym)
        if g["status"] != "absent" and g["protein"]: seqs[TAG[s]] = g["protein"]
    for (t, pid), sy in want.items():
        if sy == sym and (t, pid) in oseq and t not in seqs: seqs[t] = oseq[(t, pid)]
    aln = mafft(seqs, os.path.join(out, "alignments", sym + ".aln.fa")) if len(seqs) > 1 else dict(seqs)
    n = len(aln); Lc = len(next(iter(aln.values())))
    core = [i for i in range(Lc) if sum(1 for v in aln.values() if v[i] != "-") >= math.ceil(0.5 * n)]
    comp = {k: (sum(1 for i in core if v[i] != "-") / len(core) if core else 0.0) for k, v in aln.items()}
    ident = {}
    keys = list(aln)
    for a in keys:
        for b in keys:
            if a >= b: continue
            both = [(x, y) for x, y in zip(aln[a], aln[b]) if x != "-" and y != "-"]
            ident[a + "|" + b] = round(sum(1 for x, y in both if x == y) / len(both), 4) if len(both) >= 30 else None
    ely = [TAG[s] for s in slugs if TAG[s] in aln and comp[TAG[s]] >= 0.9]
    pid = [ident["|".join(sorted((a, b)))] for i, a in enumerate(ely) for b in ely[i + 1:]]
    pid = [x for x in pid if x is not None]
    G[sym] = dict(symbol=sym, n_seqs=n, columns=Lc, core_columns=len(core), completeness={k: round(v, 4) for k, v in comp.items()}, identity=ident,
                  mean_identity_elysia=round(sum(pid) / len(pid), 4) if pid else None, min_identity_elysia=min(pid) if pid else None, n_elysia_compared=len(ely),
                  alignment=aln)
    if all(TAG[s] in aln and comp[TAG[s]] >= 0.9 for s in slugs) and "APLCA" in aln and "PLAOC" in aln:
        cols = [i for i in core if sum(1 for v in aln.values() if v[i] == "-") <= 0.2 * n]
        for t in [TAG[s] for s in slugs] + OUTG: concat[t].append("".join(aln[t][i] for i in cols) if t in aln else "-" * len(cols))
        used.append(sym)

# ---- final status per species and gene
matrix = []
for s, c in cat.items():
    for g in c["genes"]:
        cg = G[g["symbol"]]["completeness"].get(TAG[s], 0.0); g["consensus_coverage"] = round(cg, 4)
        g["identity_to_crispata"] = G[g["symbol"]]["identity"].get("|".join(sorted((TAG[s], "ELYCR"))))
        if g["status"] != "absent":
            if cg >= 0.9 and g.get("collinear", True): g["status"] = "complete_discrepancy" if (g["frameshifts"] or g["stops"]) else "complete"
            elif cg >= 0.9: g["status"] = "fragmented"
            else: g["status"] = "partial"
    c["summary"] = dict(collections.Counter(g["status"] for g in c["genes"]))
    json.dump(c, open(os.path.join(out, s + ".clock.json"), "w"), indent=1)
with open(os.path.join(out, "clock_matrix.tsv"), "w") as fo:
    fo.write("gene\t" + "\t".join(f"{NAME[s]} status\t{NAME[s]} consensus coverage\t{NAME[s]} identity to E. crispata\t{NAME[s]} location" for s in slugs if s in cat) + "\n")
    for sym in symbols:
        row = [sym]
        for s in slugs:
            if s not in cat: continue
            g = next(x for x in cat[s]["genes"] if x["symbol"] == sym); p = g["pieces"][0] if g["pieces"] else None
            row += [g["status"], str(g.get("consensus_coverage", "")), str(g.get("identity_to_crispata") or ""), f'{p["chrom"]}:{p["start"]}-{p["end"]}' if p else ""]
        fo.write("\t".join(row) + "\n")

# ---- species tree
tree = None
if len(used) >= 5:
    ph = os.path.join(out, "concat.fa")
    with open(ph, "w") as fo:
        for t, parts in concat.items(): fo.write(f">{t}\n{''.join(parts)}\n")
    subprocess.run([BIN + "/iqtree", "-s", ph, "-m", "LG+F+G4", "-B", "1000", "-T", "2", "-o", "APLCA", "--prefix", os.path.join(out, "species_tree"), "-redo", "-quiet"], check=False)
    tf = os.path.join(out, "species_tree.treefile")
    if os.path.exists(tf):
        tree = open(tf).read().strip(); open(os.path.join(out, "species_tree.nwk"), "w").write(tree + "\n")
    for ext in (".ckp.gz", ".bionj", ".mldist", ".splits.nex", ".contree", ".log", ".model.gz", ".uniqueseq.phy"):
        f = os.path.join(out, "species_tree" + ext)
        if os.path.exists(f): os.remove(f)
json.dump(dict(species=[dict(slug=s, name=NAME[s], tag=TAG[s]) for s in slugs], outgroups=OUTG, genes=G, tree=tree, tree_genes=used,
               tree_columns=len("".join(concat[TAG[slugs[0]]])) if used else 0), open(os.path.join(out, "genus.json"), "w"))
print("genes aligned:", len(G), "| in species tree:", len(used), "genes,", len("".join(concat[TAG[slugs[0]]])) if used else 0, "columns")
for s, c in cat.items():
    print(f"{s:22s}", c["summary"], "| not complete:", ", ".join(f"{g['symbol']}:{g['status']}({round(100 * g.get('consensus_coverage', 0))}%)" for g in c["genes"] if g["status"] != "complete"))
print("tree:", tree)
