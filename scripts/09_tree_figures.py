#!/usr/bin/env python3
"""Draw, for every catalogue gene, the part of its family tree that shows the orthology call.

usage: 09_tree_figures.py catalogue.tsv allprot.tsv seeds.fa.tsv trees_dir out_dir work_dir
Writes out_dir/<symbol>.svg (styled with CSS variables so it follows the page theme) and
work_dir/tree_placement.json (support for the clade that groups the gene with its reference seeds).
"""
import sys, os, json, re
from Bio import Phylo

cat, alltsv, seedtsv, tdir, odir, wd = sys.argv[1:7]
os.makedirs(odir, exist_ok=True)
SPECIES = {"HUMAN": "Homo sapiens", "DROME": "Drosophila melanogaster", "DANPL": "Danaus plexippus", "TRICA": "Tribolium castaneum",
           "MAGGI": "Magallana gigas", "APLCA": "Aplysia californica", "LOTGI": "Lottia gigantea", "BIOGL": "Biomphalaria glabrata",
           "POMCA": "Pomacea canaliculata", "MIZYE": "Mizuhopecten yessoensis", "ELYCH": "Elysia chlorotica", "ELYMA": "Elysia marginata",
           "PLAOC": "Plakobranchus ocellatus", "ELYCRP": "Elysia crispata (Purdue)", "ELYCR": "Elysia crispata (xgElyCris1.1)",
           "CAPTE": "Capitella teleta", "NEMVE": "Nematostella vectensis", "BRAFL": "Branchiostoma floridae"}
ann = {}
for line in open(alltsv):
    f = line.rstrip("\n").split("\t"); ann[f[0].replace("|", "_")] = f
seed_org = {}
for line in list(open(seedtsv))[1:]:
    f = line.rstrip("\n").split("\t"); seed_org[f[2]] = (f[3], f[6])
SEEDSP = {"Hs": "Homo sapiens", "Dm": "Drosophila melanogaster", "Mm": "Mus musculus", "Dp": "Danaus plexippus", "Dr": "Danio rerio", "At": "Arabidopsis thaliana",
          "Cg": "Magallana gigas", "Ac": "Aplysia californica", "Em": "Elysia marginata", "Tp": "Todarodes pacificus", "My": "Mizuhopecten yessoensis",
          "Av": "Ambigolimax valentianus", "La": "Leptochiton asellus"}

def label(name):
    if name.startswith("SEED_"):
        lab = name[5:]; sp = SEEDSP.get(lab[:2], ""); gene = lab[2:]
        return f"{sp} {gene}", "seed"
    f = ann.get(name)
    tag = name.split("_")[0]
    if f is None: return name, "other"
    gene = f[3]; acc = f[2]
    if tag == "ELYCR": return f"Elysia crispata {gene} (Ensembl)", "ecr"
    if tag == "ELYCRP": return f"Elysia crispata {gene} (Purdue)", "ecr"
    g = gene if gene and not gene.startswith(("LOC", "KP79", "C0Q70", "LOTGI", "EGW08", "PoB_", "ElyMa", "CAPTE", "v1g", "BRAFL", "TcasGA2", "KGM")) else acc
    return f"{SPECIES.get(tag, tag)} {g}", "other"

def esc(s): return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

entries = []
for line in open(cat):
    if line.startswith("#") or line.startswith("symbol\t"): continue
    f = line.rstrip("\n").split("\t"); entries.append(dict(symbol=f[0], tree=f[3], seeds=f[4].split(","), anchors=f[5].split(",")))
trees = {}; placement = {}
for e in entries:
    tp = os.path.join(tdir, e["tree"] + ".nwk")
    if not os.path.exists(tp): continue
    if e["tree"] not in trees:
        t = Phylo.read(tp, "newick"); t.root_at_midpoint(); t.ladderize()
        trees[e["tree"]] = t
    t = trees[e["tree"]]
    parents = {}
    for cl in t.find_clades(order="level"):
        for ch in cl.clades: parents[ch] = cl
    tips = t.get_terminals()
    mine = [x for x in tips if any((a.startswith("ENSUFQG") and x.name.startswith("ELYCR_") and x.name.endswith(a)) or
                                   (not a.startswith("ENSUFQG") and x.name.startswith("ELYCRP_" + a + "_")) for a in e["anchors"])]
    seeds = [x for x in tips if x.name.startswith("SEED_") and x.name[5:] in e["seeds"]]
    if not mine: continue
    # smallest clade containing the gene's tips and >= 1 seed of its orthogroup
    node = t.common_ancestor(mine) if len(mine) > 1 else mine[0]
    def n_seeds(n): return [x for x in n.get_terminals() if x in seeds]
    while not n_seeds(node) and node in parents: node = parents[node]
    grp = node
    other_seeds = sorted(x.name[5:] for x in grp.get_terminals() if x.name.startswith("SEED_") and x.name[5:] not in e["seeds"])
    # the orthogroup clade: smallest clade holding the gene and ALL seeds of its orthogroup; used when it
    # contains no seed of another orthogroup, otherwise the nearest-seed clade above is reported
    use = grp
    if seeds:
        fullc = t.common_ancestor(mine + seeds)
        foreign = [x for x in fullc.get_terminals() if x.name.startswith("SEED_") and x.name[5:] not in e["seeds"]]
        if not foreign: use = fullc
    other_seeds = sorted(x.name[5:] for x in use.get_terminals() if x.name.startswith("SEED_") and x.name[5:] not in e["seeds"])
    placement[e["symbol"]] = dict(tree=e["tree"], clade_size=len(use.get_terminals()), support=(use.confidence if use.confidence is not None else None),
                                  seeds_in_clade=sorted(x.name[5:] for x in n_seeds(use)), foreign_seeds_in_clade=other_seeds,
                                  all_seeds=len(n_seeds(use)) == len(seeds))
    # show: the clade containing all seeds if not too large, else grp; plus sister lineage collapsed
    full = t.common_ancestor(mine + seeds) if seeds else grp
    show = full if len(full.get_terminals()) <= 45 else grp
    while len(show.get_terminals()) < 6 and show in parents: show = parents[show]
    # layout
    leaves = show.get_terminals()
    ROW, PADT, LABX = 17, 14, 300
    depth = {show: 0.0}
    for cl in show.find_clades(order="preorder"):
        for ch in cl.clades: depth[ch] = depth[cl] + min(ch.branch_length or 0.0, 1.5)
    maxd = max(depth[l] for l in leaves) or 1.0
    sister = None
    if show in parents:
        sib = [c for c in parents[show].clades if c is not show]
        if sib: sister = sib[0]
    nrows = len(leaves) + (1 if sister is not None else 0)
    H = PADT * 2 + ROW * nrows
    W = 760
    X0 = 26
    def x(n): return X0 + depth[n] / maxd * (LABX - X0 - 10)
    ypos = {l: PADT + ROW * i + ROW / 2 for i, l in enumerate(leaves)}
    def y(n):
        if n in ypos: return ypos[n]
        ypos[n] = sum(y(c) for c in n.clades) / len(n.clades); return ypos[n]
    y(show)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="Family tree placing {esc(e["symbol"])}" class="treefig">']
    out.append('<g fill="none" stroke="var(--tree-line,#8a958f)" stroke-width="1.2">')
    for cl in show.find_clades(order="preorder"):
        if cl.clades:
            ys = [y(c) for c in cl.clades]
            out.append(f'<path d="M{x(cl):.1f} {min(ys):.1f}V{max(ys):.1f}"/>')
            for c in cl.clades: out.append(f'<path d="M{x(cl):.1f} {y(c):.1f}H{x(c):.1f}"/>')
    if sister is not None:
        ys = PADT + ROW * len(leaves) + ROW / 2
        out.append(f'<path d="M{X0 - 14:.1f} {y(show):.1f}V{ys:.1f}M{X0 - 14:.1f} {y(show):.1f}H{X0:.1f}M{X0 - 14:.1f} {ys:.1f}H{X0 + 30:.1f}" stroke-dasharray="3 3"/>')
    out.append('</g>')
    for cl in show.find_clades(order="preorder"):
        if cl.clades and cl.confidence is not None and cl is not show:
            v = cl.confidence
            if v >= (70 if v > 1 else 0.7):
                out.append(f'<text x="{x(cl) - 3:.1f}" y="{y(cl) - 3:.1f}" text-anchor="end" class="sup">{v:.0f}</text>' if v > 1 else
                           f'<text x="{x(cl) - 3:.1f}" y="{y(cl) - 3:.1f}" text-anchor="end" class="sup">{v:.2f}</text>')
    for l in leaves:
        lab, kind = label(l.name)
        mine_flag = l in mine
        cls = "tip-me" if mine_flag else ("tip-seed" if kind == "seed" else ("tip-ecr" if kind == "ecr" else "tip"))
        out.append(f'<path d="M{x(l) + 3:.1f} {y(l):.1f}H{LABX - 6}" stroke="var(--tree-guide,#d5ded8)" stroke-width="1" stroke-dasharray="1 3"/>')
        if mine_flag: out.append(f'<circle cx="{x(l):.1f}" cy="{y(l):.1f}" r="4" class="dot-me"/>')
        elif kind == "seed": out.append(f'<rect x="{x(l) - 3.5:.1f}" y="{y(l) - 3.5:.1f}" width="7" height="7" class="dot-seed"/>')
        out.append(f'<text x="{LABX}" y="{y(l) + 4:.1f}" class="{cls}">{esc(lab)}</text>')
    if sister is not None:
        n = len(sister.get_terminals()); ss = sorted(x_.name[5:] for x_ in sister.get_terminals() if x_.name.startswith("SEED_"))[:4]
        ys = PADT + ROW * len(leaves) + ROW / 2
        txt = f"sister lineage: {n} sequence{'s' if n != 1 else ''}" + (f" (includes {', '.join(ss)})" if ss else "")
        out.append(f'<text x="{X0 + 36}" y="{ys + 4:.1f}" class="tip sister">{esc(txt)}</text>')
    out.append('</svg>')
    open(os.path.join(odir, e["symbol"] + ".svg"), "w").write("\n".join(out))
json.dump(placement, open(os.path.join(wd, "tree_placement.json"), "w"), indent=1)
for k, v in placement.items():
    print(f"{k:13s} {v['tree']:8s} clade n={v['clade_size']:3d} support={v['support']} seeds={','.join(v['seeds_in_clade'])}{' | also: ' + ','.join(v['foreign_seeds_in_clade']) if v['foreign_seeds_in_clade'] else ''}")
