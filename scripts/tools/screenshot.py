import asyncio, sys
from playwright.async_api import async_playwright
url, out = sys.argv[1], sys.argv[2]
w = int(sys.argv[3]) if len(sys.argv) > 3 else 1280
h = int(sys.argv[4]) if len(sys.argv) > 4 else 900
scroll = int(sys.argv[5]) if len(sys.argv) > 5 else 0
action = sys.argv[6] if len(sys.argv) > 6 else ""
scheme = sys.argv[7] if len(sys.argv) > 7 else "light"
async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch()
        pg = await b.new_page(viewport={'width': w, 'height': h}, color_scheme=scheme)
        errs = []
        pg.on("console", lambda m: errs.append(m.text) if m.type == "error" else None)
        pg.on("pageerror", lambda e: errs.append(str(e)))
        await pg.goto(url); await pg.wait_for_timeout(700)
        if action: await pg.evaluate(action); await pg.wait_for_timeout(500)
        if scroll: await pg.evaluate(f"window.scrollTo(0,{scroll})"); await pg.wait_for_timeout(300)
        await pg.screenshot(path=out)
        if errs: print("PAGE ERRORS:", errs[:5])
        await b.close()
asyncio.run(main())
