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
        "frequency_source": gene.provenance.get("frequency_source", "cpic"),
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
        add(f"weights:    {payload.get('frequency_source', 'cpic')} allele frequencies")
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
        add("  Add one locus, ranked by fall in modeled wrong-phenotype rate:")
        for a in adds:
            gained = ", ".join(a["alleles_gained"]) or "(partial coverage)"
            add(f"    {a['add_loci'][0]:<14} recovers {a['recovers_per_100k']:>8.1f}/100k"
                f"   gains {gained}")

    add("")
    add(f"  {payload['disclaimer']}")
    return "\n".join(out)


def dumps(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2)


def frequency_comparison(gene_symbol: str, refresh: bool = False) -> str:
    """Where CPIC and gnomAD disagree about this gene's allele frequencies.

    Neither column is the truth. The point is the third one: a panel
    decision that flips depending on which source weighted it is a decision
    resting on the choice of reference rather than on the assay.
    """
    from .gnomad import POPULATIONS, REFUSED, GnomadUnavailable  # noqa: PLC0415
    from .model import load_gene  # noqa: PLC0415

    cpic_gene = load_gene(gene_symbol, refresh=refresh)
    try:
        gnomad_gene = load_gene(gene_symbol, refresh=refresh, frequencies="gnomad")
    except GnomadUnavailable as exc:
        return f"  {exc}"

    out: list[str] = []
    add = out.append
    add(f"{cpic_gene.symbol} — CPIC vs gnomAD v4 allele frequencies")
    add("")
    add("  Rates in every other report are weighted by one of these. Where they")
    add("  disagree by more than a factor of two, the rate is a property of the")
    add("  source as much as of the panel.")
    add("")

    rows: list[tuple[float, str]] = []
    for name in sorted(cpic_gene.alleles):
        allele = cpic_gene.alleles[name]
        if len(allele.rsids) != 1:
            continue
        other = gnomad_gene.alleles.get(name)
        for pop in POPULATIONS:
            c = allele.frequencies.get(pop)
            g = (other.frequencies.get(pop) if other else None)
            if c is None or g is None:
                continue
            if c == 0 and g == 0:
                continue
            if c > 0 and g > 0:
                ratio = max(c / g, g / c)
                verdict = f"{ratio:4.0f}x  {'CPIC over' if c > g else 'CPIC under'}"
            else:
                ratio = float("inf")
                verdict = "CPIC reports zero" if c == 0 else "absent from gnomAD"
            if ratio < 2:
                continue
            # Rank by how much the disagreement could move a rate, not by the
            # ratio. A tenfold disagreement between 0.00002% and 0.0002% changes
            # nothing; CPIC calling a 0.23% allele absent changes a report.
            rows.append((abs(c - g), f"  {name:<16s} {pop:<32s}"
                                     f" {c * 100:9.5f}%  {g * 100:9.5f}%   {verdict}"))

    if not rows:
        add("  No disagreement beyond a factor of two.")
    else:
        add(f"  {'allele':<16s} {'population':<32s} {'CPIC':>10s}  {'gnomAD':>9s}"
            f"   ranked by how far apart they are")
        for _, line in sorted(rows, key=lambda r: -r[0]):
            add(line)

    add("")
    add("  Populations gnomAD is not asked for:")
    for pop, why in REFUSED.items():
        if pop in cpic_gene.populations:
            add(f"    {pop} — {why}")
    multi = sorted(n for n, a in cpic_gene.alleles.items() if len(a.rsids) > 1)
    if multi:
        add("")
        add(f"  Alleles with no gnomAD weight, defined by more than one position:"
            f" {', '.join(multi)}.")
        add("  gnomAD is keyed by position, so a haplotype has no single frequency")
        add("  there. Under --frequencies gnomad these carry no weight at all.")
    add("")
    add(f"  {gnomad_gene.provenance['frequencies']['source']}")
    return "\n".join(out)


def burden_table(gene, panel_name: str, panel: frozenset[str]) -> str:
    """Actionable-burden coverage for one panel, across every population.

    Printed as a table because a single population's figure reads as a
    fact about the panel, and the spread across populations is what shows
    which population the panel was built from.
    """
    from .burden import across_populations

    rows = across_populations(gene, panel)
    out = [
        f"Actionable burden reached by {panel_name} ({len(panel)} loci)",
        "",
        "  CPIC allele frequency carrying a decreased-function, no-function or",
        "  possible-* call, and the share of it this panel's positions detect.",
        "",
    ]
    width = max(len(p) for p in rows) if rows else 10
    for population, b in sorted(rows.items()):
        if b.percent is None:
            out.append(f"  {population:<{width}}   no published frequencies")
            continue
        out.append(
            f"  {population:<{width}}   {b.percent:5.1f}%   "
            f"({b.covered * 100:.3f}% of {b.total * 100:.3f}% allele frequency)"
        )

    unrecognised = sorted({
        f"{fn!r} ({name})"
        for b in rows.values() for name, fn in b.unrecognised
    })
    if unrecognised:
        out.append("")
        out.append(
            "  CPIC function term(s) this build does not classify, excluded "
            "from both figures:"
        )
        for row in unrecognised:
            out.append(f"      {row}")

    # Ranked by how much actionable frequency goes unseen, not by the
    # percentage. Picking the lowest percentage put East Asian at the top
    # of this list — 0% of a 0.226% burden — and hid Sub-Saharan African,
    # where a single allele this panel does not carry is 2.59% of a 2.79%
    # burden. The percentage says how well a panel does; the absolute says
    # how much is at stake.
    worst = max(
        (b for b in rows.values() if b.percent is not None and b.missed),
        key=lambda b: b.total - b.covered,
        default=None,
    )
    if worst is not None:
        unmet = (worst.total - worst.covered) * 100
        out.append("")
        out.append(
            f"  Largest unmet burden — {worst.population}, "
            f"{unmet:.3f}% of allele frequency unseen:"
        )
        for name, value in worst.missed[:5]:
            out.append(f"      {value * 100:6.3f}%  {name}")
        if worst.unpriced:
            out.append("")
            out.append(
                f"  {len(worst.unpriced)} actionable allele(s) have no published "
                f"frequency in that population and are excluded from both figures, "
                f"never counted as zero:"
            )
            out.append("      " + ", ".join(worst.unpriced[:8])
                       + (" ..." if len(worst.unpriced) > 8 else ""))
    out.append("")
    out.append("  " + DISCLAIMER)
    return "\n".join(out)
