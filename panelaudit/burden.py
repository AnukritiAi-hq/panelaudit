"""How much of a population's actionable burden a panel's positions reach.

A DIFFERENT QUESTION FROM THE REST OF THIS PACKAGE, ON PURPOSE.

The lattice answers "what does this panel report for someone who carries
X", diplotype by diplotype, and ends in a wrong-phenotype rate per 100,000
people. That is the number to quote when the subject is harm.

This module answers the question an external quality scheme asks instead:
of all the allele frequency in a population that CPIC says to act on, what
share sits at a position this panel actually types? It is a property of the
panel and the population, it needs no diplotype model and no
Hardy-Weinberg step, and it is the figure that makes two panels comparable
in one line.

It exists because that figure had been quoted in correspondence while this
package could not produce it. A number that cannot be re-run from the tool
it is attributed to is an assertion, whatever its arithmetic.

WHAT "ACTIONABLE" MEANS HERE, AND WHY IT IS A CHOICE.

An allele counts when CPIC's clinical function is neither normal nor an
admission of not knowing. In CPIC's own vocabulary that leaves decreased
function, no function, and the possible-* variants of each. Normal-function
alleles are excluded because no guideline acts on them. Unknown and
uncertain are excluded because counting them would make a panel look worse
for failing to type something nobody can act on.

That is a modelling decision and reasonable people draw it elsewhere, so
``ACTIONABLE_PREFIXES`` is one tuple in one place, every report states the
denominator it used, and the alleles CPIC gives no frequency for are
reported by name rather than silently treated as zero.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .lattice import project, signature_index
from .model import Gene

#: CPIC clinical-function values that a guideline acts on, lowercased
#: prefix match. "Possible decreased function" and "possible no function"
#: are included: CPIC assigns them an activity value and its diplotype
#: tables give them a phenotype, so a panel that misses one misses an
#: actionable result.
ACTIONABLE_PREFIXES = ("decreased function", "no function", "possible")

#: Excluded, and named so the exclusion is visible rather than implied.
NOT_ACTIONABLE_PREFIXES = ("normal function", "unknown", "uncertain", "n/a")


def is_actionable(function: str) -> bool:
    value = (function or "").strip().lower()
    if value.startswith(NOT_ACTIONABLE_PREFIXES):
        return False
    return value.startswith(ACTIONABLE_PREFIXES)


@dataclass(frozen=True)
class Burden:
    """One panel, one population, one coverage figure and its workings."""

    population: str
    #: Summed allele frequency of actionable alleles the panel can type.
    covered: float
    #: Summed allele frequency of every actionable allele CPIC prices.
    total: float
    #: Actionable alleles with a published frequency that the panel misses,
    #: largest first — the list that says *what* to add, not just how much.
    missed: tuple[tuple[str, float], ...] = field(default_factory=tuple)
    #: Actionable alleles CPIC gives no frequency for in this population.
    #: Excluded from both figures, never counted as zero.
    unpriced: tuple[str, ...] = field(default_factory=tuple)

    @property
    def fraction(self) -> float | None:
        """Covered share of priced actionable burden, or None if unpriced."""
        return (self.covered / self.total) if self.total else None

    @property
    def percent(self) -> float | None:
        f = self.fraction
        return None if f is None else 100.0 * f


def actionable_burden(gene: Gene, panel: frozenset[str], population: str) -> Burden:
    """Share of `population`'s actionable allele frequency this panel types.

    AN ALLELE COUNTS AS TYPED WHEN THE LATTICE SAYS THE PANEL CALLS IT.

    `lattice.project` is reused here rather than reimplemented, and the
    first version of this function did reimplement it: it required the
    panel to carry every position defining an allele. That is wrong for
    exactly the case this package cares most about. The CPIC four-variant
    panel reaches HapB3 through the tag SNP alone, so the strict rule
    scored it as a miss and put European coverage at 33.6% instead of
    98.7% — a figure that would have disagreed with this package's own
    wrong-phenotype report about the same panel.

    Partial visibility is therefore credited exactly as the lattice
    credits it, including the tag. Whether a tag *should* count is a real
    clinical question and this package answers it in one place, not two.

    AND THE TEST IS ON FUNCTION, NOT ON THE NAME. CPIC lists HapB3 three
    times — as HapB3 over both positions, and again as each position
    alone — so a carrier seen through the tag projects to the allele named
    ``c.1236G>A`` rather than to ``HapB3``. Demanding the same name back
    scored that as a miss. What this figure is about is whether the panel
    detects actionable variation at all, so an allele counts as covered
    when what the panel reports is itself actionable. A carrier reported
    as reference is the miss, and it is the only miss that matters: it is
    the one that reaches a prescriber as a normal result.
    """
    index = signature_index(gene)
    covered = total = 0.0
    missed: list[tuple[str, float]] = []
    unpriced: list[str] = []

    for name, allele in gene.alleles.items():
        if name == gene.reference or not is_actionable(allele.function):
            continue
        value = allele.frequencies.get(population)
        if value is None:
            unpriced.append(name)
            continue
        total += value
        reported = project(gene, name, panel, index)
        seen = (
            reported is not None
            and reported != gene.reference
            and is_actionable(gene.alleles[reported].function)
        )
        if seen:
            covered += value
        elif value > 0:
            missed.append((name, value))

    missed.sort(key=lambda row: (-row[1], row[0]))
    return Burden(
        population=population,
        covered=covered,
        total=total,
        missed=tuple(missed),
        unpriced=tuple(sorted(unpriced)),
    )


def across_populations(gene: Gene, panel: frozenset[str]) -> dict[str, Burden]:
    """The same panel measured against every population CPIC prices.

    The comparison is the point. A single population's coverage reads as a
    fact about the panel; the spread across populations is what shows the
    panel was built from one of them.
    """
    return {pop: actionable_burden(gene, panel, pop) for pop in gene.populations}
