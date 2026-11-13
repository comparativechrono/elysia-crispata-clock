#!/usr/bin/env python3
"""Download all nuccore records matching an Entrez query as FASTA.

usage: ncbi_fetch.py "<entrez term>" out.fa [db] [rettype]
"""
import sys, json, ssl, time, os, urllib.parse, urllib.request

term, out = sys.argv[1], sys.argv[2]
db = sys.argv[3] if len(sys.argv) > 3 else "nuccore"
rettype = sys.argv[4] if len(sys.argv) > 4 else "fasta"
ca = os.environ.get("SSL_CERT_FILE")
ctx = ssl.create_default_context(cafile=ca) if ca else ssl.create_default_context()
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get(url, data=None):
    for i in range(6):
        try:
            return urllib.request.urlopen(url, data=data, timeout=180, context=ctx).read()
        except Exception as e:
            time.sleep(3 + 3 * i)
    raise RuntimeError("failed: " + url[:120])

d = json.loads(get(E + "esearch.fcgi?db=%s&retmode=json&usehistory=y&retmax=0&term=%s" % (db, urllib.parse.quote(term))))
r = d["esearchresult"]; n = int(r["count"])
with open(out, "wb") as fo:
    for start in range(0, n, 200):
        url = E + "efetch.fcgi?db=%s&query_key=%s&WebEnv=%s&rettype=%s&retmode=text&retstart=%d&retmax=200" % (db, r["querykey"], r["webenv"], rettype, start)
        fo.write(get(url)); time.sleep(0.4)
print(n, "records ->", out)
