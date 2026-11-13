#!/usr/bin/env python3
"""Find an openly licensed photograph of each species on iNaturalist.

usage: 01_fetch_photos.py config/species.tsv workdir [n_candidates]
Research-grade observations whose photos carry a CC0, CC BY or CC BY-SA licence are preferred,
most-favourited first. Only if a species has none are CC BY-NC / CC BY-NC-SA photos taken
(usable on a non-commercial site with credit; flagged as such). The
candidates (medium size, 500 px) are saved as workdir/<slug>_<n>.jpg with candidates.json,
which records licence, photographer and the observation each photo comes from.
"""
import sys, os, json, time, subprocess, urllib.parse
cfg, wd = sys.argv[1], sys.argv[2]; N = int(sys.argv[3]) if len(sys.argv) > 3 else 4
API = "https://api.inaturalist.org/v1/"
def get(url):
    for i in range(4):
        r = subprocess.run(["curl", "-sS", "-m", "40", "-H", "Accept: application/json", url], capture_output=True, text=True)
        try: return json.loads(r.stdout)
        except Exception: time.sleep(2 + 2 * i)
    return {}
out = {}
for line in list(open(cfg))[1:]:
    slug, species = line.rstrip("\n").split("\t")[:2]
    t = get(API + "taxa?rank=species&per_page=5&q=" + urllib.parse.quote(species)).get("results", [])
    t = [x for x in t if x.get("name") == species] or t
    if not t: print(species, "taxon not found"); continue
    tid = t[0]["id"]; time.sleep(1.1)
    tier = "cc0,cc-by,cc-by-sa"
    obs = get(API + f"observations?taxon_id={tid}&photo_license={tier}&quality_grade=research&photos=true&order_by=votes&order=desc&per_page=30").get("results", [])
    if not obs:
        tier = "cc-by-nc,cc-by-nc-sa"; time.sleep(1.1)
        obs = get(API + f"observations?taxon_id={tid}&photo_license={tier}&quality_grade=research&photos=true&order_by=votes&order=desc&per_page=30").get("results", [])
    cands = []
    for o in obs:
        for p in o.get("photos", [])[:1]:
            if p.get("license_code") not in tier.split(","): continue
            dims = p.get("original_dimensions") or {}
            cands.append(dict(photo_id=p["id"], license=p["license_code"], attribution=p.get("attribution", ""), observer=(o.get("user") or {}).get("login", ""),
                              observer_name=(o.get("user") or {}).get("name") or "", observation="https://www.inaturalist.org/observations/%s" % o["id"],
                              place=o.get("place_guess") or "", date=o.get("observed_on") or "", faves=o.get("faves_count", 0),
                              width=dims.get("width"), height=dims.get("height"), url=p["url"].replace("square", "medium")))
    cands = [c for c in cands if (c["width"] or 1) >= (c["height"] or 1)][:N] or cands[:N]      # landscape first
    for i, c in enumerate(cands):
        f = os.path.join(wd, f"{slug}_{i}.jpg")
        subprocess.run(["curl", "-sS", "-m", "60", "-o", f, c["url"]]); c["file"] = os.path.basename(f)
    out[slug] = dict(species=species, taxon_id=tid, n_open_licence=len(obs), candidates=cands)
    print(f"{species}: taxon {tid}, {len(obs)} openly licensed research-grade observations; kept {len(cands)}: " + ", ".join(f"{c['license']} by {c['observer']} ({c['faves']} faves)" for c in cands))
    time.sleep(1.1)
json.dump(out, open(os.path.join(wd, "candidates.json"), "w"), indent=1)
