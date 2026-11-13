"""Helpers to parse miniprot --gff output."""
import collections

def parse(path):
    """Return {query: [hit dicts sorted by rank]}; each hit has exon (CDS) blocks."""
    hits = collections.defaultdict(list)
    paf = None; cur = None; sta = None
    for line in open(path):
        if line.startswith("##PAF"):
            paf = line.rstrip("\n").split("\t")
        elif line.startswith("##STA"):
            sta = line.rstrip("\n").split("\t")[1]
            continue
        elif line.startswith("#"):
            continue
        else:
            g = line.rstrip("\n").split("\t")
            if len(g) < 9: continue
            a = dict(x.split("=", 1) for x in g[8].split(";") if "=" in x)
            if g[2] == "mRNA":
                tags = {t[:2]: t[5:] for t in paf[13:] if len(t) > 5 and t[2] == ":"}
                ql, qs, qe = int(paf[2]), int(paf[3]), int(paf[4])
                cur = dict(query=paf[1], qlen=ql, qstart=qs, qend=qe, cov=(qe - qs) / ql,
                           ident=float(a.get("Identity", 0)), rank=int(a.get("Rank", 1)),
                           chrom=g[0], chrom_len=int(paf[7]), start=int(g[3]), end=int(g[4]), strand=g[6],
                           fs=int(tags.get("fs", 0)), stop=int(tags.get("st", 0)),
                           score=int(tags.get("ms", 0)), cds=[], id=a.get("ID"), trans=sta)
                hits[cur["query"]].append(cur)
            elif g[2] == "CDS" and cur is not None:
                cur["cds"].append((int(g[3]), int(g[4])))
    for q in hits:
        hits[q].sort(key=lambda h: h["rank"])
    return hits
