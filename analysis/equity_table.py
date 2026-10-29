"""Regenerate the equity table and the CYP3A5 worked case.

The European-optimal panel of a fixed size is built greedily, then measured
in every other population and compared against the panel that population
would have built for itself. The difference is reported as the cost of
being handed a panel designed elsewhere.

    python -m analysis.equity_table            # from pipelines/

Deterministic. Reads only the cached CPIC tables; no network unless the
cache is cold.
"""
from __future__ import annotations

from panelaudit.model import SUPPORTED, load_gene
from panelaudit.lattice import burden, greedy_panel

SIZE = 5
DESIGN = "European"


def populations(gene) -> list[str]:
    seen: set[str] = set()
    for allele in gene.alleles.values():
        seen.update(allele.frequencies)
    return sorted(seen)


def rows():
    for symbol in SUPPORTED:
        gene = load_gene(symbol)
        pops = populations(gene)
        if DESIGN not in pops:
            continue
        design_panel, _ = greedy_panel(gene, DESIGN, SIZE)
        for pop in pops:
            if pop == DESIGN:
                continue
            inherited = burden(gene, design_panel, pop)
            own_panel, _ = greedy_panel(gene, pop, SIZE)
            own = burden(gene, own_panel, pop)
            yield symbol, pop, inherited, own, inherited - own


def main() -> None:
    print(f"Cost of a {DESIGN}-optimal {SIZE}-locus panel, wrong calls per 100,000\n")
    print(f"{'gene':9} {'population':34} {DESIGN:>12} {'own':>10} {'cost':>12}")
    table = sorted(rows(), key=lambda r: -r[4])
    for symbol, pop, inherited, own, cost in table:
        if cost < 1.0:
            continue
        print(f"{symbol:9} {pop:34} {inherited:12,.0f} {own:10,.0f} {cost:12,.0f}")

    gene = load_gene("CYP3A5")
    two = frozenset({"rs776746", "rs10264272"})
    three = two | {"rs41303343"}
    print("\nCYP3A5: two loci (*3, *6) against three (adding *7)\n")
    print(f"{'population':34} {'two':>10} {'three':>10}")
    for pop in populations(gene):
        print(f"{pop:34} {burden(gene, two, pop):10,.0f} {burden(gene, three, pop):10,.0f}")


if __name__ == "__main__":
    main()
