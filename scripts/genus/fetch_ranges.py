#!/usr/bin/env python3
"""Download one large file as many short HTTP range requests run side by side.

usage: fetch_ranges.py URL OUT [chunk_MB=128] [parallel=6]
Long single connections are slow through some proxies while short ones are fast; each chunk
is written at its own offset in OUT, and a failed chunk is retried. Resumes if OUT.done lists
finished chunks.
"""
import sys, os, ssl, time, threading, urllib.request, queue

url, out = sys.argv[1], sys.argv[2]
chunk = int(sys.argv[3]) * 1_000_000 if len(sys.argv) > 3 else 128_000_000
par = int(sys.argv[4]) if len(sys.argv) > 4 else 6
ca = os.environ.get("SSL_CERT_FILE")
ctx = ssl.create_default_context(cafile=ca) if ca else ssl.create_default_context()

def head():
    for i in range(6):
        try:
            r = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60, context=ctx)
            return int(r.headers["Content-Length"])
        except Exception as e:
            time.sleep(2 + 2 * i)
    raise SystemExit("cannot reach " + url)

size = head()
done_file = out + ".done"
done = set(int(x) for x in open(done_file).read().split()) if os.path.exists(done_file) and os.path.exists(out) else set()
fd = os.open(out, os.O_RDWR | os.O_CREAT)
os.ftruncate(fd, size)
q = queue.Queue()
for s in range(0, size, chunk):
    if s not in done: q.put(s)
lock = threading.Lock(); got = [sum(min(chunk, size - s) for s in done)]; t0 = time.time(); failed = []

def work():
    while True:
        try: s = q.get_nowait()
        except queue.Empty: return
        e = min(s + chunk, size) - 1
        for attempt in range(8):
            try:
                req = urllib.request.Request(url, headers={"Range": f"bytes={s}-{e}"})
                with urllib.request.urlopen(req, timeout=120, context=ctx) as r:
                    pos = s
                    while True:
                        b = r.read(1 << 20)
                        if not b: break
                        os.pwrite(fd, b, pos); pos += len(b)
                if pos != e + 1: raise IOError(f"short read {pos - s} of {e + 1 - s}")
                with lock:
                    got[0] += e + 1 - s
                    with open(done_file, "a") as f: f.write(f"{s}\n")
                break
            except Exception as ex:
                time.sleep(2 + 2 * attempt)
        else:
            with lock: failed.append(s)

threads = [threading.Thread(target=work) for _ in range(par)]
for t in threads: t.start()
for t in threads: t.join()
os.close(fd)
if failed: raise SystemExit(f"{len(failed)} chunks failed; rerun to resume")
os.remove(done_file)
print(f"{out}: {size / 1e9:.2f} GB in {time.time() - t0:.0f} s ({size / 1e6 / max(1, time.time() - t0):.0f} MB/s)")
