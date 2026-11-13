#!/usr/bin/env python3
"""Fetch reference clock-protein seed sequences from UniProt.

Reads config/clock_seeds.tsv (family, category, label, taxid, gene) and writes
one FASTA with headers  >LABEL|ACCESSION|family  to the output path.
Reviewed (Swiss-Prot) entries are preferred; TrEMBL is used only as fallback.
"""
import sys, json, ssl, time, urllib.parse, urllib.request, os

cfg, out = sys.argv[1], sys.argv[2]
ca = os.environ.get("SSL_CERT_FILE") or os.environ.get("REQUESTS_CA_BUNDLE")
ctx = ssl.create_default_context(cafile=ca) if ca else ssl.create_default_context()

def query(q):
    url = ("https://rest.uniprot.org/uniprotkb/search?query=" + urllib.parse.quote(q)
           + "&fields=accession,reviewed,gene_names,protein_name,length,sequence,organism_name&format=json&size=5")
    for attempt in range(4):
        try:
            return json.load(urllib.request.urlopen(url, timeout=60, context=ctx)).get("results", [])
        except Exception as e:  # network hiccup: retry
            time.sleep(3 * (attempt + 1))
    return []

rows = [l.rstrip("\n").split("\t") for l in open(cfg) if l.strip() and not l.startswith("#")]
meta = []
with open(out, "w") as fo:
    for fam, cat, label, taxid, gene, *note in rows:
        if gene.startswith("acc:"):
            res = query(f"(accession:{gene[4:]})")
        else:
            res = query(f"(gene_exact:{gene}) AND (organism_id:{taxid}) AND (reviewed:true)")
        if not res:
            res = [] if gene.startswith("acc:") else query(f"(gene_exact:{gene}) AND (organism_id:{taxid})")
        if not res:
            print(f"!! no UniProt entry for {label} ({gene}, {taxid})", file=sys.stderr)
            continue
        res.sort(key=lambda r: (r["entryType"].startswith("UniProtKB reviewed"), r["sequence"]["length"]), reverse=True)
        r = res[0]
        acc = r["primaryAccession"]
        pd = r.get("proteinDescription", {})
        name = (pd.get("recommendedName") or (pd.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value", "")
        fo.write(f">{label}|{acc}|{fam} {name} OS={r['organism']['scientificName']}\n{r['sequence']['value']}\n")
        meta.append((fam, cat, label, acc, r["entryType"], str(r["sequence"]["length"]), name))
        print("\t".join(meta[-1]))
        time.sleep(0.25)
with open(out + ".tsv", "w") as fo:
    fo.write("family\tcategory\tlabel\taccession\tentry_type\tlength\tprotein_name\n")
    for m in meta:
        fo.write("\t".join(m) + "\n")
