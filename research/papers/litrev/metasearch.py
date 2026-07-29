#!/usr/bin/env python3
"""
Free multi-source academic metasearch + legal open-access PDF fetch.

Sources — all free, no mandatory key:
  * OpenAlex   (breadth, graph, OA status)     — email in URL (polite pool)
  * Crossref   (authoritative metadata + refs) — mailto in URL
  * Unpaywall  (best legal OA PDF per DOI)      — email in URL
  * arXiv      (physics preprints + direct PDF) — no key
Optional (set env var to enable, obtained via manual free registration):
  * Semantic Scholar  — S2_API_KEY   (isInfluential / TLDR / recommendations)

Only legal open-access / preprint PDFs are downloaded. No paywalled or
pirated sources (no Sci-Hub, no scrapers, no Google Scholar).

Usage:
  python metasearch.py search  "terahertz wire grid polarizer" --limit 15 --oa-only
  python metasearch.py resolve 10.1109/TAP.2004.838786
  python metasearch.py download "terahertz wire grid polarizer" --limit 10 --out pdfs
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import xml.etree.ElementTree as ET
from typing import Any, Dict, List, Optional

try:
    import requests
except ModuleNotFoundError:  # pragma: no cover
    # Фолбэк на стандартную библиотеку: `requests` не объявлен в requirements.txt
    # и отсутствует и в .venv, и в системном Python (проверено 2026-07-29, Сессия L).
    # Ставить пакет в ОБЩИЙ venv ради инструмента одной зоны нельзя — это тихая
    # мутация окружения других сессий. Поэтому здесь минимальная совместимая
    # обёртка над urllib: ровно та часть API requests, которой пользуется файл
    # (get(url, params, headers, timeout) -> .status_code/.json()/.text/.content,
    # плюс исключение RequestException).
    import urllib.error
    import urllib.parse
    import urllib.request
    import types

    class _Response:
        def __init__(self, status_code: int, content: bytes, encoding: str = "utf-8"):
            self.status_code = status_code
            self.content = content
            self._encoding = encoding

        @property
        def text(self) -> str:
            return self.content.decode(self._encoding, errors="replace")

        def json(self):
            return json.loads(self.text)

    class _RequestException(Exception):
        pass

    def _requests_get(url, params=None, headers=None, timeout=None):
        if params:
            clean = {k: v for k, v in params.items() if v is not None}
            url = url + ("&" if "?" in url else "?") + urllib.parse.urlencode(clean)
        req = urllib.request.Request(url, headers=headers or {})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                charset = resp.headers.get_content_charset() or "utf-8"
                return _Response(resp.status, resp.read(), charset)
        except urllib.error.HTTPError as e:
            # HTTP-ошибку возвращаем как ответ: вызывающий код сам разбирает
            # коды 429/5xx и делает бэкофф.
            return _Response(e.code, e.read() if e.fp else b"")
        except (urllib.error.URLError, OSError) as e:
            raise _RequestException(str(e)) from e

    requests = types.SimpleNamespace(get=_requests_get,
                                     RequestException=_RequestException)

EMAIL = os.getenv("CONTACT_EMAIL", "brattri3@gmail.com")
S2_API_KEY = os.getenv("S2_API_KEY", "").strip()
TIMEOUT = 40
UA = {"User-Agent": f"thz-metasearch/1.0 (mailto:{EMAIL})"}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _get(url: str, params: Optional[Dict] = None, headers: Optional[Dict] = None,
         retries: int = 4) -> Optional[requests.Response]:
    h = dict(UA)
    if headers:
        h.update(headers)
    for i in range(retries):
        try:
            r = requests.get(url, params=params, headers=h, timeout=TIMEOUT)
            if r.status_code == 200:
                return r
            if r.status_code in (429, 500, 502, 503):
                time.sleep(1.5 * (2 ** i))
                continue
            return None
        except requests.RequestException:
            time.sleep(1.0 * (2 ** i))
    return None


def _norm_doi(doi: Optional[str]) -> Optional[str]:
    if not doi:
        return None
    doi = doi.strip().lower()
    doi = re.sub(r"^https?://(dx\.)?doi\.org/", "", doi)
    return doi or None


def _norm_title(t: Optional[str]) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (t or "").lower()).strip()


def _rec(title=None, year=None, doi=None, authors=None, venue=None,
         cited_by=None, source=None, oa_status=None, pdf_url=None,
         is_influential=None, extra=None) -> Dict[str, Any]:
    return {
        "title": title, "year": year, "doi": _norm_doi(doi),
        "authors": authors or [], "venue": venue, "cited_by": cited_by,
        "sources": [source] if source else [], "oa_status": oa_status,
        "pdf_url": pdf_url, "is_influential": is_influential, "extra": extra or {},
    }


# --------------------------------------------------------------------------- #
# source adapters
# --------------------------------------------------------------------------- #
def search_openalex(query: str, limit: int = 20, oa_only: bool = False) -> List[Dict]:
    flt = f"title.search:{query}"
    if oa_only:
        flt += ",is_oa:true"
    sel = "id,doi,display_name,publication_year,authorships,cited_by_count,open_access,best_oa_location,primary_location"
    r = _get("https://api.openalex.org/works",
             {"filter": flt, "sort": "cited_by_count:desc",
              "per-page": min(limit, 200), "select": sel, "mailto": EMAIL})
    out: List[Dict] = []
    if not r:
        return out
    for w in r.json().get("results", []):
        pdf = (w.get("best_oa_location") or {}).get("pdf_url")
        out.append(_rec(
            title=w.get("display_name"), year=w.get("publication_year"),
            doi=w.get("doi"),
            authors=[a["author"]["display_name"] for a in w.get("authorships", [])[:8]],
            venue=((w.get("primary_location") or {}).get("source") or {}).get("display_name"),
            cited_by=w.get("cited_by_count"), source="openalex",
            oa_status=(w.get("open_access") or {}).get("oa_status"), pdf_url=pdf,
        ))
    return out


def search_crossref(query: str, limit: int = 20) -> List[Dict]:
    r = _get("https://api.crossref.org/works",
             {"query": query, "rows": min(limit, 100), "mailto": EMAIL,
              "select": "DOI,title,issued,author,container-title,is-referenced-by-count"})
    out: List[Dict] = []
    if not r:
        return out
    for it in r.json().get("message", {}).get("items", []):
        title = (it.get("title") or [None])[0]
        year = None
        try:
            year = it["issued"]["date-parts"][0][0]
        except Exception:
            pass
        authors = [" ".join(filter(None, [a.get("given"), a.get("family")]))
                   for a in it.get("author", [])[:8]]
        out.append(_rec(
            title=title, year=year, doi=it.get("DOI"), authors=authors,
            venue=(it.get("container-title") or [None])[0],
            cited_by=it.get("is-referenced-by-count"), source="crossref"))
    return out


def search_s2(query: str, limit: int = 20) -> List[Dict]:
    if not S2_API_KEY:
        return []
    h = {"x-api-key": S2_API_KEY}
    r = _get("https://api.semanticscholar.org/graph/v1/paper/search",
             {"query": query, "limit": min(limit, 100),
              "fields": "title,year,externalIds,authors,venue,citationCount,openAccessPdf"},
             headers=h)
    out: List[Dict] = []
    if not r:
        return out
    for p in r.json().get("data", []):
        doi = (p.get("externalIds") or {}).get("DOI")
        pdf = (p.get("openAccessPdf") or {}).get("url")
        out.append(_rec(
            title=p.get("title"), year=p.get("year"), doi=doi,
            authors=[a["name"] for a in p.get("authors", [])[:8]],
            venue=p.get("venue"), cited_by=p.get("citationCount"),
            source="s2", pdf_url=pdf, extra={"s2_id": p.get("paperId")}))
    return out


def arxiv_pdf_for_title(title: str) -> Optional[str]:
    """Best-effort: find an arXiv preprint whose title matches, return its PDF URL."""
    if not title:
        return None
    r = _get("http://export.arxiv.org/api/query",
             {"search_query": f'ti:"{title}"', "max_results": 3})
    if not r:
        return None
    ns = {"a": "http://www.w3.org/2005/Atom"}
    want = _norm_title(title)
    try:
        root = ET.fromstring(r.text)
    except ET.ParseError:
        return None
    for e in root.findall("a:entry", ns):
        et = _norm_title((e.findtext("a:title", default="", namespaces=ns)))
        # accept if titles substantially overlap
        if want and (want in et or et in want or _overlap(want, et) > 0.8):
            for link in e.findall("a:link", ns):
                if link.get("title") == "pdf":
                    return link.get("href")
    return None


def _overlap(a: str, b: str) -> float:
    sa, sb = set(a.split()), set(b.split())
    if not sa or not sb:
        return 0.0
    return len(sa & sb) / max(len(sa), len(sb))


def unpaywall_pdf(doi: str) -> Optional[str]:
    doi = _norm_doi(doi)
    if not doi:
        return None
    r = _get(f"https://api.unpaywall.org/v2/{doi}", {"email": EMAIL})
    if not r:
        return None
    j = r.json()
    loc = j.get("best_oa_location") or {}
    return loc.get("url_for_pdf") or loc.get("url")


# --------------------------------------------------------------------------- #
# merge + resolve
# --------------------------------------------------------------------------- #
def merge(*result_lists: List[Dict]) -> List[Dict]:
    by_key: Dict[str, Dict] = {}
    for lst in result_lists:
        for rec in lst:
            key = rec["doi"] or ("t:" + _norm_title(rec["title"]))
            if not key or key == "t:":
                continue
            if key in by_key:
                cur = by_key[key]
                for s in rec["sources"]:
                    if s not in cur["sources"]:
                        cur["sources"].append(s)
                for f in ("year", "venue", "cited_by", "oa_status", "pdf_url",
                          "is_influential", "doi"):
                    if not cur.get(f) and rec.get(f):
                        cur[f] = rec[f]
                if len(rec.get("authors") or []) > len(cur.get("authors") or []):
                    cur["authors"] = rec["authors"]
            else:
                by_key[key] = dict(rec)
    return sorted(by_key.values(),
                  key=lambda r: (r.get("cited_by") or 0), reverse=True)


def resolve_pdf(rec: Dict) -> Optional[str]:
    """Legal PDF only: OpenAlex OA -> Unpaywall -> arXiv preprint."""
    if rec.get("pdf_url"):
        return rec["pdf_url"]
    if rec.get("doi"):
        u = unpaywall_pdf(rec["doi"])
        if u:
            return u
    return arxiv_pdf_for_title(rec.get("title"))


def download(url: str, path: str) -> bool:
    r = _get(url)
    if not r or not r.content[:4] == b"%PDF":
        # some hosts need a direct stream without our polite headers refusal
        try:
            r = requests.get(url, headers=UA, timeout=TIMEOUT)
        except requests.RequestException:
            return False
    if r and r.content[:4] == b"%PDF":
        with open(path, "wb") as f:
            f.write(r.content)
        return True
    return False


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def _print(records: List[Dict]) -> None:
    for r in records:
        src = "+".join(r["sources"])
        oa = r.get("oa_status") or "-"
        print(f"[{r.get('year')}] c={r.get('cited_by') or 0:<5} {oa:<7} {src:<18} "
              f"{(r.get('title') or '')[:90]}")
        if r.get("doi"):
            print(f"        doi:{r['doi']}")


def cmd_search(a) -> None:
    res = merge(search_openalex(a.query, a.limit, a.oa_only),
                search_crossref(a.query, a.limit),
                search_s2(a.query, a.limit))
    _print(res)
    if a.json:
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=False, indent=2)
        print(f"\nsaved: {a.json} ({len(res)} records)")


def cmd_resolve(a) -> None:
    rec = _rec(doi=a.doi, title=a.title)
    # enrich from openalex if only doi given
    if a.doi and not a.title:
        r = _get(f"https://api.openalex.org/works/doi:{_norm_doi(a.doi)}",
                 {"mailto": EMAIL})
        if r:
            w = r.json()
            rec["title"] = w.get("display_name")
            rec["pdf_url"] = (w.get("best_oa_location") or {}).get("pdf_url")
            rec["oa_status"] = (w.get("open_access") or {}).get("oa_status")
    pdf = resolve_pdf(rec)
    print(json.dumps({"doi": rec["doi"], "title": rec["title"],
                      "oa_status": rec["oa_status"], "pdf": pdf},
                     ensure_ascii=False, indent=2))


def cmd_download(a) -> None:
    os.makedirs(a.out, exist_ok=True)
    res = merge(search_openalex(a.query, a.limit, oa_only=True),
                search_crossref(a.query, a.limit))
    ok = 0
    for r in res:
        pdf = resolve_pdf(r)
        if not pdf:
            print(f"  no legal PDF: {(r.get('title') or '')[:70]}")
            continue
        name = (r.get("doi") or _norm_title(r.get("title"))[:60]).replace("/", "_") + ".pdf"
        path = os.path.join(a.out, name)
        if download(pdf, path):
            ok += 1
            print(f"  OK  {name}")
        else:
            print(f"  fail {name}  ({pdf})")
        time.sleep(0.3)
    print(f"\ndownloaded {ok}/{len(res)} into {a.out}")


def main() -> None:
    p = argparse.ArgumentParser(description="Free academic metasearch (legal OA only)")
    sub = p.add_subparsers(required=True)

    s = sub.add_parser("search")
    s.add_argument("query")
    s.add_argument("--limit", type=int, default=20)
    s.add_argument("--oa-only", action="store_true")
    s.add_argument("--json")
    s.set_defaults(func=cmd_search)

    r = sub.add_parser("resolve")
    r.add_argument("doi")
    r.add_argument("--title")
    r.set_defaults(func=cmd_resolve)

    d = sub.add_parser("download")
    d.add_argument("query")
    d.add_argument("--limit", type=int, default=15)
    d.add_argument("--out", default="pdfs")
    d.set_defaults(func=cmd_download)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
