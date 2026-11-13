#!/usr/bin/env python3
"""Save the chosen photograph of each species for the site, with its credit.

usage: 02_choose_photos.py candidates.json photo_workdir docs/assets/photos slug=index ...
The choice among the candidates of 01_fetch_photos.py is made by eye (a clear view of the
animal). The large version is fetched and reduced to 800 px across; credits.json records the
licence, photographer and source observation, which the site prints with every photo.
"""
import sys, os, json, subprocess
from PIL import Image
cands = json.load(open(sys.argv[1])); wd, out = sys.argv[2], sys.argv[3]
LIC = {"cc0": ("CC0 1.0", "https://creativecommons.org/publicdomain/zero/1.0/"), "cc-by": ("CC BY 4.0", "https://creativecommons.org/licenses/by/4.0/"),
       "cc-by-sa": ("CC BY-SA 4.0", "https://creativecommons.org/licenses/by-sa/4.0/"), "cc-by-nc": ("CC BY-NC 4.0", "https://creativecommons.org/licenses/by-nc/4.0/"),
       "cc-by-nc-sa": ("CC BY-NC-SA 4.0", "https://creativecommons.org/licenses/by-nc-sa/4.0/")}
credits = {}
for a in sys.argv[4:]:
    slug, i = a.split("="); c = cands[slug]["candidates"][int(i)]
    tmp = os.path.join(wd, slug + "_large.jpg")
    subprocess.run(["curl", "-sS", "-m", "90", "-o", tmp, c["url"].replace("medium", "large")], check=True)
    im = Image.open(tmp).convert("RGB"); im.thumbnail((800, 800)); im.save(os.path.join(out, slug + ".jpg"), quality=82, optimize=True)
    credits[slug] = dict(file=slug + ".jpg", species=cands[slug]["species"], photographer=c["observer_name"] or c["observer"], inaturalist_user=c["observer"],
                         licence=LIC[c["license"]][0], licence_url=LIC[c["license"]][1], non_commercial="nc" in c["license"], source=c["observation"],
                         place=c["place"], date=c["date"], attribution=c["attribution"], width=im.size[0], height=im.size[1])
    print(slug, im.size, credits[slug]["licence"], credits[slug]["photographer"], "|", c["place"][:40], c["date"])
json.dump(credits, open(os.path.join(out, "credits.json"), "w"), indent=1)
