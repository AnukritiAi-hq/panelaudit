"""panelaudit — what does this panel report for someone who carries X?

    python -m panelaudit --gene DPYD --panel cpic4 --population "Central/South Asian"
    python -m panelaudit --gene DPYD --panel my_assay.csv --json report.json
    python -m panelaudit --all-panels --gene TPMT
    python -m panelaudit --self-check
"""

from __future__ import annotations

import argparse
import pathlib
import sys

from . import report as _report
from .lattice import audit, rank_additions, summarise
from .model import EXCLUDED, SUPPORTED, load_gene
from .panels import presets, resolve

DEFAULT_POPULATION = "Central/South Asian"


def _run(gene_symbol: str, spec: str, population: str | None, refresh: bool,
         worst: int, suggest: bool) -> dict:
    gene = load_gene(gene_symbol, refresh=refresh)
    label, loci, panel_file = resolve(spec)

    if population and population not in gene.populations:
        raise SystemExit(
            f"{gene.symbol}: CPIC publishes no frequencies for {population!r}. "
            f"Available: {', '.join(gene.populations)}"
        )

    # A panel for another gene produces a full page of confident findings
    # that mean nothing. Refuse rather than render it.
    gene_loci: set[str] = set()
    for allele in gene.alleles.values():
        gene_loci |= allele.rsids
    if loci and not (loci & gene_loci):
        raise SystemExit(
            f"none of the {len(loci)} loci in {spec!r} are {gene.symbol} "
            f"allele-defining positions. Check --gene."
        )

    cases = audit(gene, loci, population)
    summary = summarise(cases, weighted=bool(population))
    additions = rank_additions(gene, loci, population) if (suggest and population) else []
    payload = _report.as_json(gene, label, loci, population, cases, summary, additions, worst)

    if panel_file is not None and panel_file.rejected:
        payload["panel"]["rejected_entries"] = [
            {"line": n, "value": v} for n, v in panel_file.rejected
        ]
    return payload


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        prog="panelaudit", description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__,
    )
    ap.add_argument("--gene", default="DPYD", help=f"one of: {', '.join(SUPPORTED)}")
    ap.add_argument("--genes", action="store_true", help="list supported and excluded genes")
    ap.add_argument("--panel", default="cpic4", help="preset name, or a path to a panel CSV")
    ap.add_argument("--population", default=DEFAULT_POPULATION,
                    help="CPIC population for frequency weighting; 'none' to disable")
    ap.add_argument("--all-panels", action="store_true", help="every preset, summary only")
    ap.add_argument("--json", metavar="PATH", help="write the full report as JSON")
    ap.add_argument("--worst", type=int, default=10, help="worst cases to show")
    ap.add_argument("--no-suggest", action="store_true", help="skip the addition ranking")
    ap.add_argument("--refresh", action="store_true", help="re-fetch CPIC tables")
    ap.add_argument("--self-check", action="store_true", help="run assertions and exit")
    args = ap.parse_args(argv)

    if args.genes:
        print("audited:")
        for g in SUPPORTED:
            print(f"  {g}")
        print("\nnot audited, on purpose:")
        for g, why in EXCLUDED.items():
            print(f"  {g:<9} {why}")
        return 0

    if args.self_check:
        from .selfcheck import run  # noqa: PLC0415

        return run()

    population = None if args.population.lower() == "none" else args.population

    if args.all_panels:
        gene = load_gene(args.gene, refresh=args.refresh)
        print(f"{gene.symbol} — every preset"
              + (f", weighted for {population}" if population else ""))
        print()
        for key in presets():
            cases = audit(gene, frozenset(presets()[key][1]), population)
            s = summarise(cases, weighted=bool(population))
            rate = "" if s.wrong_per_100k is None else f"  {s.wrong_per_100k:>9,.1f}/100k wrong"
            missed = ("" if s.missed_poor_per_100k is None
                      else f"  {s.missed_poor_per_100k:>7.2f}/100k missed Poor")
            print(f"  {key:<16} {len(presets()[key][1]):>3} loci   "
                  f"correct {s.correct:>3}/{s.total}   wrong {s.wrong:>3}   "
                  f"no-call {s.no_call}{rate}{missed}")
        return 0

    payload = _run(args.gene, args.panel, population, args.refresh,
                   args.worst, not args.no_suggest)
    print(_report.as_text(payload, worst=args.worst))

    if args.json:
        pathlib.Path(args.json).write_text(_report.dumps(payload))
        print(f"\n  JSON written to {args.json}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
