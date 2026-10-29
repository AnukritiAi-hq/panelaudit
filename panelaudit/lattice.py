"""The failure lattice: what a panel reports for someone who carries X.

A coverage report lists the positions a panel interrogates. That is not
the question a prescriber has. The question is: *if a patient carries this
genotype, what does this panel say?* There are three answers:

    CORRECT   the panel resolves it and the phenotype is right
    NO_CALL   the panel cannot resolve it and says so. Honest failure.
    WRONG     the panel resolves it to a different, confident phenotype

WRONG is the finding. It needs no patient data: it is combinatorics over
allele definitions, weighted by published allele frequencies, so the whole
lattice is computable before anyone has been sequenced.

Weighting uses Hardy-Weinberg over CPIC's own per-population allele
frequencies. That is an approximation and it is stated as one: real
populations are structured, and consanguinity raises homozygote counts
above HW in some of the groups this is most relevant to. It is the right
order of magnitude for ranking panels against each other, which is what it
is used for, and it is never presented as a prediction for an individual.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

from .model import Gene, rank

CORRECT = "CORRECT"
NO_CALL = "NO_CALL"
WRONG = "WRONG"

#: The panel reports more enzyme activity than the patient has. The patient
#: looks safer than they are and the error favours giving more drug.
UNDERCALL = "UNDERCALL"
#: The panel reports less. Conservative, and defensible.
OVERCALL = "OVERCALL"
#: One side of the comparison is off the activity scale, so there is no
#: direction to state. Reported as its own class rather than guessed.
INCOMPARABLE = "INCOMPARABLE"


@dataclass(frozen=True)
class Case:
    true_diplotype: tuple[str, str]
    called_diplotype: tuple[str, str] | None
    true_phenotype: str
    called_phenotype: str | None
    verdict: str
    direction: str | None
    #: Population frequency of this diplotype under Hardy-Weinberg, or None
    #: when CPIC publishes no frequency for one of the alleles.
    frequency: float | None
    #: Both alleles carry a variant the panel can see, so an unphased assay
    #: cannot tell this compound heterozygote from two variants in cis.
    phase_ambiguous: bool

    @property
    def label(self) -> str:
        return "/".join(self.true_diplotype)

    @property
    def called_label(self) -> str:
        return "/".join(self.called_diplotype) if self.called_diplotype else "—"


def signature_index(gene: Gene) -> dict[frozenset[str], str]:
    """``defining rsid set -> allele name``, for resolving what a panel sees."""
    index: dict[frozenset[str], str] = {}
    for name, allele in gene.alleles.items():
        index.setdefault(allele.rsids, name)
    return index


def project(gene: Gene, allele_name: str, panel: frozenset[str],
            index: dict[frozenset[str], str]) -> str | None:
    """What the panel reports for a carrier of ``allele_name``.

    None means the panel saw variation it cannot assign to a defined
    allele — a genuine no-call.

    The important case is the silent one: an allele whose defining
    positions are all absent is indistinguishable from reference. The
    panel reports no uncertainty, because from its own point of view there
    is none.
    """
    visible = gene.alleles[allele_name].rsids & panel
    if not visible:
        return gene.reference
    return index.get(visible)


def _allele_frequency(gene: Gene, name: str, population: str) -> float | None:
    if name == gene.reference:
        return gene.reference_frequency(population)[0]
    return gene.alleles[name].frequencies.get(population)


def diplotype_frequency(gene: Gene, a: str, b: str, population: str) -> float | None:
    """Hardy-Weinberg. None when either allele has no usable frequency.

    The reference allele is handled by ``Gene.reference_frequency``, which
    derives a value when CPIC publishes none and says so. Without that,
    every heterozygote-with-reference goes unweighted — and for DPYD that
    is most of the real population.
    """
    fa = _allele_frequency(gene, a, population)
    fb = _allele_frequency(gene, b, population)
    if fa is None or fb is None:
        return None
    return fa * fb if a == b else 2.0 * fa * fb


def audit(gene: Gene, panel: frozenset[str], population: str | None = None) -> list[Case]:
    """Classify every diplotype the allele table can express against one panel."""
    index = signature_index(gene)
    cases: list[Case] = []

    # Unordered pairs with replacement: a diplotype is a multiset, so
    # *2A/Reference and Reference/*2A are the same person.
    for a, b in itertools.combinations_with_replacement(sorted(gene.alleles), 2):
        true_pheno = gene.phenotype_of(a, b)
        if true_pheno is None:
            continue  # CPIC publishes no phenotype for this function pair
        freq = diplotype_frequency(gene, a, b, population) if population else None

        pa = project(gene, a, panel, index)
        pb = project(gene, b, panel, index)
        if pa is None or pb is None:
            cases.append(Case((a, b), None, true_pheno, None, NO_CALL, None, freq, False))
            continue

        ambiguous = (
            pa != gene.reference and pb != gene.reference and a != b
        )
        called_pheno = gene.phenotype_of(pa, pb)
        if called_pheno is None:
            cases.append(Case((a, b), (pa, pb), true_pheno, None, NO_CALL, None, freq, ambiguous))
            continue

        if called_pheno == true_pheno:
            cases.append(
                Case((a, b), (pa, pb), true_pheno, called_pheno, CORRECT, None, freq, ambiguous)
            )
            continue

        true_rank, called_rank = rank(true_pheno), rank(called_pheno)
        if true_rank is None or called_rank is None or true_rank[0] != called_rank[0]:
            # Off-scale on either side, or two different scales. There is
            # no direction to state and inventing one would be worse than
            # saying so.
            direction = INCOMPARABLE
        else:
            direction = UNDERCALL if called_rank[1] > true_rank[1] else OVERCALL
        cases.append(
            Case((a, b), (pa, pb), true_pheno, called_pheno, WRONG, direction, freq, ambiguous)
        )
    return cases


@dataclass(frozen=True)
class Summary:
    total: int
    correct: int
    wrong: int
    no_call: int
    undercall: int
    overcall: int
    incomparable: int
    phase_ambiguous: int
    #: Per 100,000 people of the population, under Hardy-Weinberg.
    wrong_per_100k: float | None
    #: Wrong, in the dangerous direction, where the truth is Poor Metabolizer.
    missed_poor_per_100k: float | None
    #: Diplotypes that could not be weighted because CPIC gives no frequency.
    unweighted_cases: int

    @property
    def weighted(self) -> bool:
        return self.wrong_per_100k is not None


def summarise(cases: list[Case], weighted: bool = True) -> Summary:
    counts = {CORRECT: 0, WRONG: 0, NO_CALL: 0}
    directions = {UNDERCALL: 0, OVERCALL: 0, INCOMPARABLE: 0}
    wrong_burden = 0.0
    missed_poor = 0.0
    unweighted = 0
    ambiguous = 0

    for c in cases:
        counts[c.verdict] += 1
        if c.direction:
            directions[c.direction] += 1
        if c.phase_ambiguous:
            ambiguous += 1
        if c.verdict != CORRECT and c.frequency is None:
            unweighted += 1
        if c.verdict == WRONG and c.frequency is not None:
            wrong_burden += c.frequency
            if c.direction == UNDERCALL and c.true_phenotype == "Poor Metabolizer":
                missed_poor += c.frequency

    return Summary(
        total=len(cases),
        correct=counts[CORRECT],
        wrong=counts[WRONG],
        no_call=counts[NO_CALL],
        undercall=directions[UNDERCALL],
        overcall=directions[OVERCALL],
        incomparable=directions[INCOMPARABLE],
        phase_ambiguous=ambiguous,
        wrong_per_100k=wrong_burden * 100_000 if weighted else None,
        missed_poor_per_100k=missed_poor * 100_000 if weighted else None,
        unweighted_cases=unweighted,
    )


def candidate_loci(gene: Gene, panel: frozenset[str]) -> set[str]:
    """Every defining position the panel does not already interrogate."""
    everything: set[str] = set()
    for allele in gene.alleles.values():
        everything |= allele.rsids
    return everything - panel


@dataclass(frozen=True)
class Addition:
    rsids: tuple[str, ...]
    alleles_gained: tuple[str, ...]
    wrong_before: float
    wrong_after: float

    @property
    def recovered_per_100k(self) -> float:
        return self.wrong_before - self.wrong_after


def rank_additions(gene: Gene, panel: frozenset[str], population: str) -> list[Addition]:
    """Which single locus, added next, removes the most wrong calls.

    This is the actionable half. A coverage gap tells a lab what it cannot
    see; this tells it what to buy. Each candidate is scored by the carrier
    burden it moves out of WRONG, so the ranking is by people rather than
    by positions.

    Greedy and one-at-a-time by design. Loci interact — adding the HapB3
    causal variant changes what the tag SNP is worth — so re-run after each
    addition rather than trusting a static top-five.
    """
    base = summarise(audit(gene, panel, population)).wrong_per_100k or 0.0
    out: list[Addition] = []
    for rsid in sorted(candidate_loci(gene, panel)):
        widened = panel | {rsid}
        after = summarise(audit(gene, widened, population)).wrong_per_100k or 0.0
        gained = tuple(
            sorted(
                name for name, a in gene.alleles.items()
                if rsid in a.rsids and a.rsids <= widened
            )
        )
        out.append(Addition((rsid,), gained, base, after))
    out.sort(key=lambda a: (-a.recovered_per_100k, a.rsids))
    return out


def greedy_panel(
    gene: Gene, population: str, size: int, seed: frozenset[str] | None = None
) -> tuple[frozenset[str], list[float]]:
    """Build the ``size``-locus panel that minimises wrong calls here.

    Greedy: take the locus that removes the most carrier burden, then
    re-rank and repeat. Greedy is not guaranteed optimal — loci interact,
    and a pair can be worth more than either alone — but re-ranking after
    every pick captures the common interactions, and the result is a
    concrete panel rather than an unreachable bound.

    Returns the panel and the remaining wrong-per-100k after each pick, so
    a caller can see where adding loci stops paying.
    """
    panel = frozenset(seed or ())
    trace: list[float] = []
    for _ in range(size):
        options = rank_additions(gene, panel, population)
        if not options or options[0].recovered_per_100k <= 0:
            break
        panel = panel | set(options[0].rsids)
        trace.append(options[0].wrong_after)
    return panel, trace


def burden(gene: Gene, panel: frozenset[str], population: str) -> float:
    """Wrong calls per 100,000, for comparing panels across populations."""
    return summarise(audit(gene, panel, population)).wrong_per_100k or 0.0
