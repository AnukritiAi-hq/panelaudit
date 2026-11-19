"""gnomAD v4 allele frequencies, as a second opinion on CPIC's.

**Why this exists.** Every rate this package prints is a frequency-weighted
count, and until now the only weights were CPIC's published ones. Checked
against gnomAD v4 on 2026-11-18, those disagree in both directions and
sometimes by a lot: CPIC overstates ``*2A`` — the most clinically used DPYD
allele — by 4x in Central/South Asian, 5x in European and 12x in African
American, and records ``c.868A>G`` as absent in African ancestry where gnomAD
measures 0.23%.

The tell is in CPIC's own numbers. Four DPYD alleles carry an identical
African-American frequency of 0.0031083, which is 1/321.7 — exactly what one
observed allele in about 161 genotyped people gives. That is a singleton
floor from a small cohort, not a measurement.

Neither source is "the truth". CPIC curates from the pharmacogenomic
literature and is the basis of the guideline; gnomAD is a much larger
convenience sample with its own ascertainment. The useful thing is not to
pick a winner but to show a laboratory where the two disagree, because a
panel decision that flips between them is a decision resting on the choice
of reference rather than on the assay.

**What this module will not do.** It maps a CPIC population to a gnomAD
genetic-ancestry group only where the correspondence is defensible, and
refuses the rest by name. gnomAD has no Oceanian group and its ``afr``
pools African American with continental African, so CPIC's "Sub-Saharan
African" has no counterpart: answering it with ``afr`` would quietly assert
the two populations are one.
"""

from __future__ import annotations

import json
import pathlib
import urllib.error
import urllib.request
from typing import Any

from .cpic import CACHE, _UA

API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r4"

#: CPIC population -> the gnomAD genetic-ancestry group(s) pooled for it.
#: "European" pools non-Finnish and Finnish because CPIC's label is the
#: broader one; the two are reported separately by :func:`subgroups` so a
#: reader can see what the pooling hid.
POPULATIONS: dict[str, tuple[str, ...]] = {
    "African American/Afro-Caribbean": ("afr",),
    "Central/South Asian": ("sas",),
    "East Asian": ("eas",),
    "European": ("nfe", "fin"),
    "Latino": ("amr",),
}

#: CPIC populations with no defensible gnomAD counterpart, and why.
REFUSED: dict[str, str] = {
    "Sub-Saharan African": (
        "gnomAD's afr group pools African American with continental African, "
        "and is already used for CPIC's African American/Afro-Caribbean. "
        "Answering this population with the same group would assert the two "
        "are one."
    ),
    "Oceanian": "gnomAD v4 publishes no Oceanian genetic-ancestry group.",
    "Near Eastern": (
        "gnomAD's mid group is Middle Eastern and much smaller than the "
        "others; CPIC's Near Eastern is not the same construct and the "
        "sample is too thin to stand in for it."
    ),
}


def _cache_path(gene: str) -> pathlib.Path:
    return CACHE / f"{gene}.gnomad_r4.json"


def _fetch(gene: str) -> list[dict[str, Any]]:
    query = (
        "query($g:String!){gene(gene_symbol:$g,reference_genome:GRCh38)"
        "{variants(dataset:" + DATASET + "){rsids "
        "exome{populations{id ac an}} genome{populations{id ac an}}}}}"
    )
    body = json.dumps({"query": query, "variables": {"g": gene}}).encode()
    req = urllib.request.Request(
        API,
        data=body,
        headers={"Content-Type": "application/json", "User-Agent": _UA},
    )
    with urllib.request.urlopen(req, timeout=180) as fh:  # noqa: S310 — pinned host
        payload = json.load(fh)
    if payload.get("errors"):
        raise GnomadUnavailable(str(payload["errors"])[:200])
    return payload["data"]["gene"]["variants"]


class GnomadUnavailable(RuntimeError):
    """No cached copy and the API could not be reached."""


def counts(gene: str, refresh: bool = False) -> tuple[dict[str, dict[str, list[int]]], dict]:
    """``rsid -> gnomAD group -> [ac, an]``, exome and genome summed.

    Summing the two call sets double-counts the minority of samples
    sequenced both ways. gnomAD's own joint frequency would be better; this
    is close enough to answer "do the two sources disagree by a factor that
    changes a decision", which is all it is used for, and it is stated here
    so nobody mistakes it for a published figure.
    """
    path = _cache_path(gene)
    if path.exists() and not refresh:
        raw = json.loads(path.read_text())
        return raw["counts"], raw["provenance"]

    try:
        variants = _fetch(gene)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
        raise GnomadUnavailable(
            f"{gene}: no cached gnomAD copy and the API is unreachable ({exc}). "
            f"Run with --refresh when it is back."
        ) from exc

    out: dict[str, dict[str, list[int]]] = {}
    for variant in variants:
        for rsid in variant.get("rsids") or ():
            per_group = out.setdefault(rsid, {})
            for call_set in ("exome", "genome"):
                block = variant.get(call_set)
                if not block:
                    continue
                for pop in block.get("populations") or ():
                    pid = pop["id"]
                    # skip the sex splits and the sub-population breakdowns,
                    # which would otherwise be counted twice
                    if "_" in pid or pid in ("XX", "XY"):
                        continue
                    slot = per_group.setdefault(pid, [0, 0])
                    slot[0] += pop["ac"]
                    slot[1] += pop["an"]

    provenance = {
        "source": f"{API} · {DATASET} · gene={gene} · GRCh38",
        "variants": len(variants),
        "note": "exome and genome allele counts summed; see counts() docstring",
    }
    path.write_text(json.dumps({"counts": out, "provenance": provenance}, indent=1))
    return out, provenance


def frequencies(gene: str, refresh: bool = False) -> tuple[dict[str, dict[str, float]], dict]:
    """``rsid -> CPIC population name -> allele frequency``.

    Only the populations in :data:`POPULATIONS` appear. A population with no
    observed alleles and a non-zero denominator is a measured zero and is
    reported as 0.0; a population gnomAD did not call at all is absent, not
    zero.
    """
    raw, provenance = counts(gene, refresh=refresh)
    out: dict[str, dict[str, float]] = {}
    for rsid, groups in raw.items():
        per_pop: dict[str, float] = {}
        for cpic_pop, gnomad_groups in POPULATIONS.items():
            ac = an = 0
            for gid in gnomad_groups:
                slot = groups.get(gid)
                if slot:
                    ac += slot[0]
                    an += slot[1]
            if an:
                per_pop[cpic_pop] = ac / an
        if per_pop:
            out[rsid] = per_pop
    return out, provenance


def subgroups(gene: str, rsid: str, refresh: bool = False) -> dict[str, tuple[int, int]]:
    """Raw ``group -> (ac, an)`` for one variant, so pooling can be inspected."""
    raw, _ = counts(gene, refresh=refresh)
    return {k: (v[0], v[1]) for k, v in sorted(raw.get(rsid, {}).items())}
