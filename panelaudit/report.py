"""Rendering an audit: JSON for a machine, text for a reader.

Every report carries its provenance — the CPIC table versions, when they
were fetched, the allele table, the panel and the population. A coverage
report that cannot name what it was computed against is not reproducible,
and reproducibility is most of what a customer is buying.
"""

from __future__ import annotations

import datetime as _dt
import json
from typing import Any

from . import cpic
from .lattice import Case, Summary, UNDERCALL, WRONG
from .model import Gene

DISCLAIMER = (
    "Investigational. Describes an assay, not a person. Computes no dose, "
    "names no patient, and asserts nothing clinical that CPIC does not "
    "publish. Frequencies are Hardy-Weinberg estimates over published "
    "allele frequencies, suitable for ranking panels and not for "
    "predicting an individual."
)


def _case_json(c: Case) -> dict[str, Any]:
    return {
        "carried": c.label,
        "reported_as": c.called_label,
        "true_phenotype": c.true_phenotype,
        "reported_phenotype": c.called_phenotype,
        "verdict": c.verdict,
        "direction": c.direction,
        "frequency": c.frequency,
        "per_100k": None if c.frequency is None else round(c.frequency * 100_000, 3),
        "phase_ambiguous": c.phase_ambiguous,
    }


def as_json(
    gene: Gene,
    panel_name: str,
    panel: frozenset[str],
    population: str | None,
    cases: list[Case],
    summary: Summary,
    additions: list | None = None,
    worst: int = 25,
) -> dict[str, Any]:
    dangerous = sorted(
        (c for c in cases if c.verdict == WRONG and c.direction == UNDERCALL),
        key=lambda c: (-(c.frequency or 0.0), c.true_phenotype != "Poor Metabolizer", c.label),
    )
    return {
        "report": "panel failure lattice",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "disclaimer": DISCLAIMER,
        "gene": gene.symbol,
        "panel": {"name": panel_name, "loci": sorted(panel), "count": len(panel)},
        "population": population,
        "summary": {
            "diplotypes": summary.total,
            "correct": summary.correct,
            "wrong": summary.wrong,
            "no_call": summary.no_call,
            "undercall": summary.undercall,
            "overcall": summary.overcall,
            "incomparable": summary.incomparable,
            "phase_ambiguous": summary.phase_ambiguous,
            "wrong_per_100k": None if summary.wrong_per_100k is None
            else round(summary.wrong_per_100k, 2),
            "missed_poor_metabolizers_per_100k": None
            if summary.missed_poor_per_100k is None
            else round(summary.missed_poor_per_100k, 3),
            "unweighted_cases": summary.unweighted_cases,
        },
        "alleles_without_published_frequency": list(gene.unweighted),
        "reference_frequency": None if population is None else {
            "value": gene.reference_frequency(population)[0],
            "obtained": gene.reference_frequency(population)[1],
        },
        "worst_cases": [_case_json(c) for c in dangerous[:worst]],
        "recommended_additions": [
            {
                "add_loci": list(a.rsids),
                "alleles_gained": list(a.alleles_gained),
                "recovers_per_100k": round(a.recovered_per_100k, 2),
            }
            for a in (additions or [])
            if a.recovered_per_100k > 0
        ],
        "provenance": {
            "cpic": gene.provenance,
            "citation": cpic.CITATION,
            "allele_tables": "bundled per-gene tables, regenerated verbatim from CPIC's API (CC0 1.0)",
        },
    }


def as_text(payload: dict[str, Any], worst: int = 10) -> str:
    s = payload["summary"]
    out: list[str] = []
    add = out.append

    add(f"{payload['gene']} panel audit — {payload['panel']['name']} "
        f"({payload['panel']['count']} loci)")
    if payload["population"]:
        add(f"population: {payload['population']}")
    add("")
    add(f"  {s['correct']} of {s['diplotypes']} diplotypes reported correctly")
    add(f"  {s['wrong']} wrong, {s['no_call']} no-call")
    add(f"  of the wrong: {s['undercall']} undercall, {s['overcall']} overcall, "
        f"{s['incomparable']} incomparable")
    if s["wrong_per_100k"] is not None:
        add("")
        add(f"  {s['wrong_per_100k']:,.1f} per 100,000 people receive a wrong phenotype")
        add(f"  {s['missed_poor_metabolizers_per_100k']:,.2f} per 100,000 are Poor "
            f"Metabolizers reported as something milder")
    ref = payload.get("reference_frequency")
    if ref and ref["value"] is not None and not ref["obtained"].startswith("published"):
        add(f"  reference allele frequency {ref['value']:.4f} — {ref['obtained']}")
    if s["unweighted_cases"]:
        add(f"  ({s['unweighted_cases']} cases carry no published frequency and are "
            f"excluded from the rates above, never counted as zero)")

    worst_cases = payload["worst_cases"][:worst]
    if worst_cases:
        add("")
        add("  Worst cases, by carrier frequency:")
        for c in worst_cases:
            rate = "" if c["per_100k"] is None else f"{c['per_100k']:>8.2f}/100k  "
            add(f"    {rate}{c['carried']:<26} reported {c['reported_as']:<24} "
                f"{c['true_phenotype']} -> {c['reported_phenotype']}")

    adds = payload["recommended_additions"]
    if adds:
        add("")
        add("  Add one locus, ranked by carrier burden recovered:")
        for a in adds:
            gained = ", ".join(a["alleles_gained"]) or "(partial coverage)"
            add(f"    {a['add_loci'][0]:<14} recovers {a['recovers_per_100k']:>8.1f}/100k"
                f"   gains {gained}")

    add("")
    add(f"  {payload['disclaimer']}")
    return "\n".join(out)


def dumps(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2)
