#!/usr/bin/env python3
"""Summarise where E. crispata proteins fall in a family tree.

usage: tree_report.py tree.nwk [--ascii]
For each ELYCR/ELYCRP tip, walks towards the (midpoint) root and reports the
smallest clade that also contains a labelled reference seed, with its support,
size and the seeds it contains.
"""
import sys
from Bio import Phylo

tree = Phylo.read(sys.argv[1], "newick")
tree.root_at_midpoint()
if "--ascii" in sys.argv:
    for t in tree.get_terminals():
        t.name = t.name[:38]
    Phylo.draw_ascii(tree, column_width=150)
    sys.exit()
parents = {}
for cl in tree.find_clades(order="level"):
    for ch in cl.clades:
        parents[ch] = cl
tips = {cl: [t.name for t in cl.get_terminals()] for cl in tree.find_clades()}
for t in tree.get_terminals():
    if not (t.name.startswith("ELYCR_") or t.name.startswith("ELYCRP_")):
        continue
    node, out = t, []
    while node in parents and len(out) < 2:
        node = parents[node]
        seeds = sorted(n[5:] for n in tips[node] if n.startswith("SEED_"))
        if seeds and (not out or seeds != out[-1][2]):
            sup = node.confidence if node.confidence is not None else float("nan")
            out.append((len(tips[node]), sup, seeds))
    d = "  ||  ".join(f"clade n={n} sup={s:.2f}: {','.join(sd[:8])}{'...' if len(sd)>8 else ''}" for n, s, sd in out)
    print(f"{t.name:40s} {d}")
