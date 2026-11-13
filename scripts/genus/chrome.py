"""Shared page frame for every page of the site (header, species tabs, footer)."""
import html

LABEL = {"Ck1de": "Ck1δ/ε", "Ck1a": "Ck1α", "Ck2a": "Ck2α", "Ck2b": "Ck2β", "Pp1a": "Pp1α", "Pp1b": "Pp1β", "Pp2aC": "Pp2a-C",
         "Pp2aB55": "Pp2a-B55", "Pp2aB56e": "Pp2a-B56ε", "Pp2aB56d": "Pp2a-B56δ", "Phr6-4": "Phr 6-4", "PhrCPD": "Phr CPD",
         "CryDASH1": "Cry-DASH 1", "CryDASH2": "Cry-DASH 2", "Ropsin1": "r-opsin 1", "Ropsin2": "r-opsin 2", "RopsinL1": "r-opsin-like 1",
         "RopsinL2": "r-opsin-like 2", "GoOpsin1": "Go-opsin 1", "GoOpsin2": "Go-opsin 2"}
CATS = {"core": "Core loop", "photolyase": "Photolyase family", "loop2": "Second loop", "kinase": "Kinase", "phosphatase": "Phosphatase",
        "degradation": "Degradation", "modulator": "Modulator", "output": "Output", "light": "Opsin"}
STATUS_CLASS = {"complete": "ok", "complete_discrepancy": "warn", "fragmented": "part", "partial": "part", "divergent": "miss", "absent": "miss"}
STATUS_WORD = {"complete": "intact", "complete_discrepancy": "whole, sequence differs", "fragmented": "in pieces", "partial": "partial", "divergent": "divergent copy only", "absent": "missing"}

SITE_NAV = [("index.html", "Species"), ("compare.html", "Compare"), ("methods.html", "Methods and data")]
SPECIES_TABS = [("index.html", "Clock genes"), ("kleptoplasts.html", "Whose chloroplasts"), ("gene-transfer.html", "Gene transfer"), ("methods.html", "Methods for this species")]

def page(title, body, root="", description="", site_current=None, species=None, tabs=None, tab_current=None, scripts="", today="", sources=""):
    """root: '' for pages at the top of the site, '../' for pages in a species folder.
    species: scientific name shown in the species bar; tabs: [(href, label)] within that folder."""
    e = html.escape
    nav = "".join(f'<a href="{root}{h}"{" aria-current=\"page\"" if h == site_current else ""}>{t}</a>' for h, t in SITE_NAV)
    sub = ""
    if species:
        links = "".join(f'<a href="{h}"{" aria-current=\"page\"" if h == tab_current else ""}>{t}</a>' for h, t in (tabs or []))
        sub = f'<div class="subhead"><div class="wrap"><a class="sp" href="index.html"><i>{e(species)}</i></a><nav aria-label="Pages for this species">{links}</nav></div></div>'
    desc = description or "Circadian clock genes, kleptoplast identity and tests for algal gene transfer in the sea slugs of the genus Elysia, from public genome data."
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{e(title)}</title>
<meta name="description" content="{e(desc)}">
<link rel="preload" href="{root}assets/fonts/BricolageGrotesque-300-800.woff2" as="font" type="font/woff2" crossorigin>
<link rel="stylesheet" href="{root}assets/site.css">
</head>
<body>
<header class="masthead"><div class="wrap">
<a class="brand" href="{root}index.html"><span class="ld" aria-hidden="true"></span><span><em>Elysia</em> clocks</span></a>
<nav aria-label="Site">{nav}</nav>
</div></header>
{sub}
<main>
{body}
</main>
<footer class="site"><div class="wrap">
{f"<p>{sources}</p>" if sources else ""}
<p>Built {today} from public genome data. The analysis is independent of the groups that produced the genomes and has not been peer reviewed. Gene roles are inferred from homology with flies and mammals; nothing here measures a rhythm.</p>
<p>Scripts, tables and sequences are in the repository that hosts this site. Photographs are credited on the <a href="{root}methods.html#photos">methods page</a>. Typefaces: Bricolage Grotesque and IBM Plex Mono (SIL Open Font Licence).</p>
</div></footer>
{scripts}
</body>
</html>
"""

def redirect(target, title="Moved"):
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<meta http-equiv="refresh" content="0; url={target}"><link rel="canonical" href="{target}">
<script>location.replace("{target}" + location.search + location.hash);</script></head>
<body><p>This page has moved to <a href="{target}">{target}</a>.</p></body></html>
"""
