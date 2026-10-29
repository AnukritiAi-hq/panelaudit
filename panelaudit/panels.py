"""Where a panel comes from: a customer's file, or a named preset.

The customer's file is the real input to the business. A molecular lab
has a list of the positions its assay interrogates; everything else in
this package is downstream of reading that list honestly.

Accepted:

* **CSV/TSV** with a column named rsid, rsids, dbsnp, snp, variant or
  marker — case and spacing ignored. Other columns are kept out of the
  way, not parsed.
* **A plain list**, one rsID per line, ``#`` for comments.

Anything that is not a well-formed rsID is reported rather than dropped.
A panel file silently losing three positions is exactly the class of
defect this package exists to find in other people's pipelines.
"""

from __future__ import annotations

import csv
import pathlib
import re
from dataclasses import dataclass

RSID = re.compile(r"^rs\d+$")

_COLUMN_NAMES = ("rsid", "rsids", "dbsnp", "dbsnpid", "snp", "variant", "marker")


@dataclass(frozen=True)
class PanelFile:
    name: str
    loci: frozenset[str]
    #: Values that were not well-formed rsIDs, with the line they came from.
    rejected: tuple[tuple[int, str], ...]
    source: str

    @property
    def clean(self) -> bool:
        return not self.rejected


def _column(fieldnames: list[str] | None) -> str | None:
    for field in fieldnames or []:
        if field and field.strip().lower().replace(" ", "_") in _COLUMN_NAMES:
            return field
    return None


def read_panel(path: str | pathlib.Path, name: str | None = None) -> PanelFile:
    path = pathlib.Path(path)
    text = path.read_text()
    loci: set[str] = set()
    rejected: list[tuple[int, str]] = []

    sample = text[:4096]
    delimiter = "\t" if sample.count("\t") > sample.count(",") else ","
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    column = _column(reader.fieldnames)

    if column is not None:
        for lineno, row in enumerate(reader, start=2):
            raw = (row.get(column) or "").strip()
            if not raw:
                continue
            (loci.add(raw) if RSID.match(raw) else rejected.append((lineno, raw)))
    else:
        # No recognised header: treat it as one rsID per line.
        for lineno, line in enumerate(text.splitlines(), start=1):
            raw = line.split("#", 1)[0].strip()
            if not raw:
                continue
            (loci.add(raw) if RSID.match(raw) else rejected.append((lineno, raw)))

    return PanelFile(
        name=name or path.stem,
        loci=frozenset(loci),
        rejected=tuple(rejected),
        source=str(path),
    )


#: Named presets, vendored as data in ``presets/panels.json``. Every entry
#: is a published specification and records its source; none is original to
#: this package.
_PRESETS = pathlib.Path(__file__).resolve().parent / "presets" / "panels.json"


def presets() -> dict[str, tuple[str, tuple[str, ...]]]:
    """The named panels, read from the bundled specification file."""
    import json  # noqa: PLC0415

    raw = json.loads(_PRESETS.read_text())["panels"]
    return {k: (v["label"], tuple(v["loci"])) for k, v in raw.items()}


def preset_sources() -> dict[str, str]:
    """Where each preset's locus list is published. Keyed as ``presets()``."""
    import json  # noqa: PLC0415

    raw = json.loads(_PRESETS.read_text())["panels"]
    return {k: v.get("source", "") for k, v in raw.items()}


def resolve(spec: str) -> tuple[str, frozenset[str], PanelFile | None]:
    """A preset name or a path. Returns ``(label, loci, file_or_None)``."""
    named = presets()
    if spec in named:
        label, loci = named[spec]
        return f"{spec} — {label}", frozenset(loci), None
    path = pathlib.Path(spec)
    if not path.exists():
        raise FileNotFoundError(
            f"{spec!r} is neither a panel file nor a preset. "
            f"Presets: {', '.join(sorted(named))}"
        )
    panel = read_panel(path)
    return f"{panel.name} (from {panel.source})", panel.loci, panel
