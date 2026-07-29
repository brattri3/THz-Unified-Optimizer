import json, re, sys, time, urllib.error, urllib.parse, urllib.request
from pathlib import Path

UA = {"User-Agent": "thz-metasearch/1.0 (mailto:brattri3@gmail.com)"}
REPO = Path(r"C:\Users\pop\Antigravity IDE\THz-Unified-Optimizer")
FILES = [REPO / "research/papers/INDEX.md",
         REPO / "research/papers/litrev/RETROSPECTIVE_REVIEW.md",
         REPO / "research/papers/litrev/LITERATURE_REVIEW.md"]

DOI_RE = re.compile(r"10\.\d{4,9}/[^\s\)\]\|,;\"'<>`]+", re.I)

found = {}
for f in FILES:
    txt = f.read_text(encoding="utf-8", errors="replace")
    for m in DOI_RE.finditer(txt):
        d = m.group(0).rstrip(".").rstrip(")").lower()
        found.setdefault(d, set()).add(f.name)

print(f"уникальных DOI: {len(found)}\n")
bad = []
for d in sorted(found):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(d, safe="/().-")
    try:
        req = urllib.request.Request(url + "?mailto=brattri3@gmail.com", headers=UA)
        with urllib.request.urlopen(req, timeout=45) as r:
            m = json.loads(r.read().decode("utf-8", "replace"))["message"]
        title = (m.get("title") or [""])[0][:70]
        print(f"OK    {d}\n        {title}")
    except urllib.error.HTTPError as e:
        print(f"DEAD  {d}   [HTTP {e.code}]   в: {', '.join(sorted(found[d]))}")
        bad.append(d)
    except Exception as e:
        print(f"ERR   {d}   {e}")
        bad.append(d)
    time.sleep(0.2)

print(f"\nБИТЫХ: {len(bad)}")
for d in bad:
    print("  " + d)
