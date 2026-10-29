"""CPIC's API, fetched once and cached on disk.

CPIC curated content is CC0 1.0 — free of restriction, attribution
requested rather than required. It is the only external data this package
uses, and that is deliberate: one source, one licence, no ShareAlike and
no sale restriction. See ``pipelines/reference/licensing.md``.

The cache is a plain JSON file per endpoint. It exists so the audit runs
offline and so a report can be reproduced months later against the exact
tables it was written from: every cache file records the URL it came from
and the UTC time it was fetched, and the report carries that stamp.

    python -I -m panelaudit.cpic --refresh DPYD TPMT NUDT15
"""

from __future__ import annotations

import datetime as _dt
import json
import pathlib
import urllib.error
import urllib.request
from typing import Any

API = "https://api.cpicpgx.org/v1"
CACHE = pathlib.Path(__file__).resolve().parent / "cache"

#: Attribution CPIC asks for. Carried into every report.
CITATION = (
    "CPIC® (Clinical Pharmacogenetics Implementation Consortium), "
    "www.clinpgx.org. Curated content is CC0 1.0. Content is subject to "
    "update; the version and access date for this run are recorded above."
)


class CpicUnavailable(RuntimeError):
    """No cache and no network. Says which, because the fix differs."""


def _cache_path(endpoint: str, gene: str) -> pathlib.Path:
    return CACHE / f"{gene}.{endpoint}.json"


def _fetch(endpoint: str, gene: str) -> list[dict[str, Any]]:
    url = f"{API}/{endpoint}?genesymbol=eq.{gene}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as fh:  # noqa: S310 — pinned host
        return json.load(fh)


def load(endpoint: str, gene: str, refresh: bool = False) -> tuple[list[dict], dict]:
    """Return ``(rows, provenance)``, from cache unless ``refresh``.

    Provenance is not decoration. A coverage report that cannot name the
    table version it was computed against is not reproducible, and
    reproducibility is most of what a customer is paying for.
    """
    path = _cache_path(endpoint, gene)
    if not refresh and path.exists():
        blob = json.loads(path.read_text())
        return blob["rows"], blob["provenance"]

    try:
        rows = _fetch(endpoint, gene)
    except (urllib.error.URLError, TimeoutError) as exc:
        if path.exists():
            blob = json.loads(path.read_text())
            return blob["rows"], blob["provenance"]
        raise CpicUnavailable(
            f"{gene}/{endpoint}: no cached copy and CPIC's API is unreachable "
            f"({exc}). Run with --refresh while online to populate the cache."
        ) from exc

    provenance = {
        "source": f"{API}/{endpoint}?genesymbol=eq.{gene}",
        "fetched_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "rows": len(rows),
        "version": _version_of(rows),
        "licence": "CC0 1.0",
    }
    CACHE.mkdir(exist_ok=True)
    path.write_text(json.dumps({"provenance": provenance, "rows": rows}, indent=1))
    return rows, provenance


def _version_of(rows: list[dict]) -> str | None:
    """CPIC stamps a version on each row; one table should carry one."""
    versions = {str(r["version"]) for r in rows if r.get("version") is not None}
    if not versions:
        return None
    return sorted(versions)[-1] if len(versions) > 1 else next(iter(versions))


def alleles(gene: str, refresh: bool = False) -> tuple[list[dict], dict]:
    return load("allele", gene, refresh)


def diplotypes(gene: str, refresh: bool = False) -> tuple[list[dict], dict]:
    return load("diplotype", gene, refresh)


def _main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description="Populate the CPIC cache.")
    ap.add_argument("genes", nargs="+")
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args()
    for gene in args.genes:
        for endpoint in ("allele", "diplotype"):
            rows, prov = load(endpoint, gene, refresh=args.refresh)
            print(f"{gene:<8} {endpoint:<10} {len(rows):>5} rows  "
                  f"version {prov.get('version')}  {prov['fetched_utc']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
