"""Gene definitions, assembled from the bundled allele tables and CPIC's.

Two sources, each used for what it is authoritative about:

* The bundled per-gene TSV in ``data/`` says **which alleles can be
  called and from which positions**. It is regenerated verbatim from
  CPIC's API, so it is CPIC's content in a calling table's shape. The
  CPIC version every result was computed against is carried in the
  report's provenance block, not in this filename.
* CPIC's ``allele`` and ``diplotype`` endpoints supply **population
  frequencies and the phenotype of every function pair**.

Phenotype is keyed on the **unordered pair of allele functions**, taken
verbatim from CPIC's diplotype rows, and not on an activity score.
Measured 2026-10-29: only DPYD publishes activity values at all, TPMT and
NUDT15 return ``None`` for every allele, and TPMT's diplotype rows carry
``totalactivityscore: "n/a"``. An activity-score model is therefore a
DPYD model wearing a general name.

The per-gene mapping is not a formality either. ``Decreased/Normal`` is
**Intermediate Metabolizer in DPYD and Normal Metabolizer in TPMT**. One
global table would be correct for exactly one gene.

Unknown frequency is carried as ``None`` and never as zero. An allele CPIC
gives no frequency for is not an allele nobody carries, and that
difference is the thing this package exists to expose in other tools.
"""

from __future__ import annotations

import csv
import pathlib
import re
from dataclasses import dataclass, field

from . import cpic

#: Genes this package audits, and the ones it refuses to. Each exclusion
#: is a modelling limit rather than missing work, and saying so is the
#: point: a tool that quietly produced a CYP2D6 report from rsIDs alone
#: would be wrong in the way this whole package exists to expose.
SUPPORTED = ("DPYD", "TPMT", "NUDT15", "CYP2B6", "CYP2C19", "CYP2C9", "CYP3A5")

EXCLUDED: dict[str, str] = {
    "CYP2D6": (
        "phenotype depends on copy number as well as allele function — CPIC "
        "maps one function pair to several phenotypes, and gene deletions and "
        "duplications cannot be represented by a set of rsIDs at all. An rsID "
        "panel model would be confidently wrong."
    ),
    "G6PD": (
        "X-linked. Hemizygous males have no second allele, so the diplotype "
        "model here does not apply, and CPIC publishes no reference allele row "
        "for it."
    ),
    "SLCO1B1": (
        "the bundled table for this gene is a variant table rather than an "
        "allele table, so there is no allele-to-position mapping to audit."
    ),
    "NAT2": (
        "the bundled allele table records function as \"Slow Acetylator\", a "
        "phenotype term, where CPIC's diplotype rows use \"Decreased function\". "
        "The two tables are not speaking the same language, so no phenotype "
        "resolves. This is a defect in the NAT2 allele table, not a limit of "
        "the model — and NAT2 is CPIC Level C with no guideline anyway."
    ),
    "VKORC1": "CPIC publishes no diplotype table.",
    "CYP1A2": "CPIC publishes no allele or diplotype table.",
    "CYP3A4": "CPIC publishes no allele or diplotype table.",
}


class GeneNotSupported(ValueError):
    """Refused on purpose, with the reason."""


#: Per-gene allele definition tables, vendored with this package. Each is
#: regenerated verbatim from CPIC's API and carries CPIC's version in its
#: filename; CPIC content is CC0 1.0. See LICENSE-DATA.
ENGINE_DATA = pathlib.Path(__file__).resolve().parent / "data"

#: CPIC writes the no-variant allele "Reference" for DPYD and "*1" for the
#: star-allele genes. Both mean: no variant was found.
_REFERENCE_NAMES = ("Reference", "*1")

#: ``c.496A>G (M166V)`` -> ``M166V``.
_PARENTHESISED = re.compile(r"\(([^)]+)\)\s*$")


def _fn_key(function: str) -> str:
    """Canonical key for an allele function.

    CPIC's own tables disagree with themselves about case. Measured
    2026-10-29: the bundled CYP2B6 table carries both ``Normal Function``
    and ``Normal function``, and the diplotype rows use their own casing
    again. Matching on the literal string silently returns no phenotype
    for almost every pair, which reads as a clean gene rather than a
    broken lookup — this package had exactly that bug until CYP2B6,
    CYP3A5 and NAT2 all reported suspiciously perfect scores.

    Whitespace is collapsed too, for the same reason. Collisions are
    checked where the mapping is built, so folding can never merge two
    functions CPIC means to keep apart.
    """
    return " ".join(function.split()).casefold()

#: CPIC does not speak one phenotype language. Measured across its own
#: diplotype tables on 2026-10-29: 19 distinct terms over four unrelated
#: scales. A gene answers in exactly one of them, and comparing a term on
#: one scale to a term on another is meaningless, so the rank carries its
#: scale and callers must check both.
#:
#: Within a scale the ordering is by enzyme activity, least to most. It is
#: a vocabulary ordering, not a clinical judgement: it only says which
#: direction an error runs in.
PHENOTYPE_SCALES: dict[str, dict[str, int]] = {
    # CYP2B6, CYP2C19, CYP2C9, CYP2D6, CYP3A5, DPYD, NAT2, NUDT15, TPMT
    "metabolizer": {
        "Poor Metabolizer": 0,
        "Likely Poor Metabolizer": 0,
        "Intermediate Metabolizer": 1,
        "Possible Intermediate Metabolizer": 1,
        "Likely Intermediate Metabolizer": 1,
        "Normal Metabolizer": 2,
        "Rapid Metabolizer": 3,
        "Ultrarapid Metabolizer": 4,
    },
    # SLCO1B1 answers about transporter function, not metabolism.
    "function": {
        "Poor Function": 0,
        "Decreased Function": 1,
        "Possible Decreased Function": 1,
        "Normal Function": 2,
        "Increased Function": 3,
    },
    # G6PD answers in deficiency terms.
    "deficiency": {
        "Deficient with CNSHA": 0,
        "Deficient": 1,
        "Normal": 2,
    },
}

#: Terms that carry no position on any scale. ``Indeterminate`` means the
#: genotype did not resolve; ``Variable`` is CPIC's term for an allele
#: whose activity is not consistent. Neither is a point on a line, and
#: ordering them would invent a severity CPIC does not state.
OFF_SCALE = ("Indeterminate", "Variable")

#: Flat lookup, built once. A term appearing on two scales would be a data
#: problem rather than a modelling choice, so it is caught here.
PHENOTYPE_RANK: dict[str, tuple[str, int]] = {}
for _scale, _terms in PHENOTYPE_SCALES.items():
    for _term, _value in _terms.items():
        if _term in PHENOTYPE_RANK:  # pragma: no cover - guards the table above
            raise RuntimeError(f"{_term!r} appears on two scales")
        PHENOTYPE_RANK[_term] = (_scale, _value)


class UnknownPhenotype(ValueError):
    """CPIC emitted a term this package has no ordering for.

    Raised rather than guessed. A new CPIC term silently ranked as Normal
    would turn a real finding into a clean report.
    """


@dataclass(frozen=True)
class Allele:
    name: str
    rsids: frozenset[str]
    function: str
    #: population -> allele frequency. Empty when CPIC publishes none.
    frequencies: dict[str, float] = field(default_factory=dict, compare=False)

    @property
    def has_frequency(self) -> bool:
        return bool(self.frequencies)


@dataclass(frozen=True)
class Gene:
    symbol: str
    alleles: dict[str, Allele]
    reference: str
    #: frozenset({function_a, function_b}) -> phenotype, verbatim from CPIC
    phenotype_by_functions: dict[frozenset[str], str]
    populations: tuple[str, ...]
    provenance: dict

    def phenotype(self, function_a: str, function_b: str) -> str | None:
        """CPIC's phenotype for this pair of allele functions."""
        return self.phenotype_by_functions.get(
            frozenset((_fn_key(function_a), _fn_key(function_b)))
        )

    def phenotype_of(self, allele_a: str, allele_b: str) -> str | None:
        return self.phenotype(self.alleles[allele_a].function, self.alleles[allele_b].function)

    def reference_frequency(self, population: str) -> tuple[float | None, str]:
        """Frequency of the no-variant allele, and how it was obtained.

        CPIC publishes this for TPMT and NUDT15 and **not for DPYD**, in
        any population (measured 2026-10-29). Without it every
        heterozygote-with-reference diplotype is unweighted, and those are
        the common ones, so the headline rate would silently describe only
        compound heterozygotes.

        So it is derived as ``1 - sum(published variant frequencies)`` when
        CPIC does not publish it. The derivation absorbs any allele CPIC
        gives no frequency for into the reference, which understates those
        alleles. That is a real error, it is bounded by how many such
        alleles there are, and the caller is told both facts rather than
        handed a bare number.
        """
        published = self.alleles[self.reference].frequencies.get(population)
        if published is not None:
            return published, "published by CPIC"

        total = 0.0
        known = 0
        for name, allele in self.alleles.items():
            if name == self.reference:
                continue
            value = allele.frequencies.get(population)
            if value is not None:
                total += value
                known += 1
        if known == 0:
            return None, "no published frequencies for this population"
        absorbed = len(self.unweighted)
        note = (
            f"derived as 1 - sum of {known} published variant frequencies; "
            f"CPIC publishes none for this gene"
        )
        if absorbed:
            note += f"; {absorbed} allele(s) without a published frequency are absorbed into it"
        return max(0.0, 1.0 - total), note

    @property
    def unweighted(self) -> tuple[str, ...]:
        """Callable alleles CPIC gives no frequency for. Reported, never zeroed."""
        return tuple(
            sorted(
                n for n, a in self.alleles.items()
                if n != self.reference and not a.has_frequency
            )
        )


def rank(phenotype: str) -> tuple[str, int] | None:
    """``(scale, position)``, or None when the term is off-scale.

    Two phenotypes are only comparable when their scales match. A caller
    that ignores the scale will happily report that Deficient is milder
    than Ultrarapid Metabolizer.
    """
    if phenotype in OFF_SCALE:
        return None
    try:
        return PHENOTYPE_RANK[phenotype]
    except KeyError as exc:
        raise UnknownPhenotype(
            f"CPIC emitted {phenotype!r}, which has no ordering here. Add it to "
            "PHENOTYPE_SCALES or OFF_SCALE deliberately; do not let it default."
        ) from exc


def _allele_rows(symbol: str) -> list[dict[str, str]]:
    folder = ENGINE_DATA / symbol.lower()
    try:
        path = next(folder.glob(f"{symbol}_alleles_*.tsv"))
    except StopIteration as exc:  # pragma: no cover - configuration error
        raise FileNotFoundError(f"no allele table for {symbol} under {folder}") from exc
    lines = [ln for ln in path.read_text().splitlines() if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def _phenotype_by_functions(rows: list[dict], symbol: str) -> dict[frozenset[str], str]:
    seen: dict[frozenset[str], set[str]] = {}
    spellings: dict[str, set[str]] = {}
    for r in rows:
        f1, f2 = (r.get("function1") or "").strip(), (r.get("function2") or "").strip()
        result = (r.get("generesult") or "").strip()
        if f1 and f2 and result:
            for raw in (f1, f2):
                spellings.setdefault(_fn_key(raw), set()).add(raw)
            seen.setdefault(frozenset((_fn_key(f1), _fn_key(f2))), set()).add(result)

    # Folding case must never merge two functions CPIC distinguishes. If a
    # key ever covers spellings that differ by more than case or spacing,
    # that is a data change and the audit should stop, not average over it.
    for key, variants in spellings.items():
        if len({" ".join(v.split()).casefold() for v in variants}) > 1:  # pragma: no cover
            raise ValueError(f"{symbol}: function key {key!r} covers {sorted(variants)}")
    ambiguous = {tuple(sorted(k)): sorted(v) for k, v in seen.items() if len(v) > 1}
    if ambiguous:
        raise ValueError(
            f"{symbol}: CPIC gives one function pair several phenotypes, so it "
            f"cannot be resolved deterministically: {ambiguous}"
        )
    return {k: next(iter(v)) for k, v in seen.items()}


def _name_index(cpic_rows: list[dict]) -> dict[str, dict]:
    """Every name the allele table might use, mapped onto CPIC's row."""
    index: dict[str, dict] = {}
    for row in cpic_rows:
        index[row["name"]] = row
        short = _PARENTHESISED.search(row["name"])
        if short:
            index.setdefault(short.group(1), row)
    return index


def load_gene(symbol: str, refresh: bool = False) -> Gene:
    symbol = symbol.upper()
    if symbol in EXCLUDED:
        raise GeneNotSupported(f"{symbol} is not audited here: {EXCLUDED[symbol]}")
    allele_rows, allele_prov = cpic.alleles(symbol, refresh)
    diplo_rows, diplo_prov = cpic.diplotypes(symbol, refresh)

    phenotypes = _phenotype_by_functions(diplo_rows, symbol)
    index = _name_index(allele_rows)

    reference_row = next((index[n] for n in _REFERENCE_NAMES if n in index), None)
    if reference_row is None:
        raise ValueError(f"{symbol}: CPIC publishes no reference allele row")

    populations: set[str] = set()
    for row in allele_rows:
        populations.update((row.get("frequency") or {}).keys())

    reference_name = reference_row["name"]
    alleles: dict[str, Allele] = {
        reference_name: Allele(
            name=reference_name,
            rsids=frozenset(),
            function=reference_row["clinicalfunctionalstatus"].strip(),
            frequencies=dict(reference_row.get("frequency") or {}),
        )
    }

    for row in _allele_rows(symbol):
        name = row["allele"].strip()
        cpic_row = index.get(name)
        alleles[name] = Allele(
            name=name,
            # One allele may be defined by several positions; HapB3 is two.
            rsids=frozenset(r.strip() for r in row["rsid"].split(";") if r.strip()),
            function=row["function"].strip(),
            frequencies=dict((cpic_row or {}).get("frequency") or {}),
        )

    # Every allele function must exist in CPIC's own diplotype vocabulary
    # for this gene. Without this check a mismatch yields no phenotype for
    # almost every pair, and a gene with no resolvable diplotypes reads as
    # a gene with no failures — a clean report produced by a broken join.
    # This is what NAT2 does: the engine table stores "Slow Acetylator", a
    # phenotype term, in its function column, where CPIC writes
    # "Decreased function".
    known = {k for pair in phenotypes for k in pair}
    unmapped = sorted(
        {a.function for a in alleles.values() if _fn_key(a.function) not in known}
    )
    if unmapped:
        raise GeneNotSupported(
            f"{symbol}: allele function(s) {unmapped} appear nowhere in CPIC's "
            f"diplotype vocabulary for this gene ({sorted(known)}). The allele "
            "table and CPIC's are not speaking the same language, so no "
            "phenotype can be resolved."
        )

    gene = Gene(
        symbol=symbol,
        alleles=alleles,
        reference=reference_name,
        phenotype_by_functions=phenotypes,
        populations=tuple(sorted(populations)),
        provenance={"allele": allele_prov, "diplotype": diplo_prov},
    )

    # Fail at load, not halfway through a customer's report.
    for phenotype in set(phenotypes.values()):
        rank(phenotype)
    return gene
