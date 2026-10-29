"""The smallest set of assertions that fails if the audit breaks.

Run with ``python -m panelaudit --self-check``. No framework, no fixtures:
this is a check, not a test suite. Each assertion is here because getting
it wrong would produce a report that looks fine and is not.
"""

from __future__ import annotations

import tempfile
import pathlib

from .lattice import (
    CORRECT, NO_CALL, UNDERCALL, WRONG,
    audit, diplotype_frequency, project, rank_additions, signature_index, summarise,
)
from .model import (
    EXCLUDED, SUPPORTED, GeneNotSupported, UnknownPhenotype, load_gene, rank,
)
from .panels import presets, read_panel

SAS = "Central/South Asian"


def run() -> int:
    # --- the per-gene phenotype model -----------------------------------
    dpyd, tpmt, nudt15 = load_gene("DPYD"), load_gene("TPMT"), load_gene("NUDT15")

    # The reason a global phenotype table would be wrong: the same pair of
    # functions means different things in different genes.
    assert dpyd.phenotype("Decreased function", "Normal function") == "Intermediate Metabolizer"
    assert tpmt.phenotype("Decreased function", "Normal function") == "Normal Metabolizer"

    # CPIC names the no-variant allele differently per gene.
    assert dpyd.reference == "Reference"
    assert tpmt.reference == "*1" and nudt15.reference == "*1"

    # An unrecognised phenotype must raise, never default to healthy.
    assert rank("Indeterminate") is None
    assert rank("Variable") is None
    try:
        rank("Brisk Metabolizer")
    except UnknownPhenotype:
        pass
    else:  # pragma: no cover
        raise AssertionError("an unknown phenotype was silently ranked")

    # CPIC speaks several phenotype languages. Terms carry their scale, and
    # two terms on different scales are not comparable — without this the
    # audit would cheerfully rank Deficient against Ultrarapid Metabolizer.
    assert rank("Poor Metabolizer")[0] == "metabolizer"
    assert rank("Poor Function")[0] == "function"
    assert rank("Deficient")[0] == "deficiency"
    assert rank("Ultrarapid Metabolizer")[1] > rank("Normal Metabolizer")[1]
    assert rank("Poor Metabolizer")[1] < rank("Intermediate Metabolizer")[1]

    # Every gene we claim to support must actually load, and every gene we
    # exclude must refuse with a reason rather than fail obscurely.
    for symbol in SUPPORTED:
        loaded = load_gene(symbol)
        assert loaded.alleles and loaded.phenotype_by_functions, symbol
    for symbol in EXCLUDED:
        try:
            load_gene(symbol)
        except GeneNotSupported as exc:
            assert EXCLUDED[symbol] in str(exc), symbol
        else:  # pragma: no cover
            raise AssertionError(f"{symbol} should be refused")

    # --- projection ------------------------------------------------------
    index = signature_index(dpyd)
    cpic4 = frozenset(presets()["cpic4"][1])

    # An allele the panel does not type is reported as reference, silently.
    assert project(dpyd, "c.2279C>T", cpic4, index) == "Reference"
    # One it does type is reported as itself.
    assert project(dpyd, "*2A", cpic4, index) == "*2A"
    # HapB3 typed by its tag alone resolves to the tag's own allele, which
    # CPIC gives the same decreased function — so this is NOT a wrong call,
    # and a guard that flagged it would be crying wolf.
    assert project(dpyd, "HapB3", cpic4, index) == "c.1236G>A"

    # --- the lattice -----------------------------------------------------
    cases = {c.label: c for c in audit(dpyd, cpic4, SAS)}

    headline = cases["*2A/c.2279C>T"]
    assert headline.verdict == WRONG
    assert headline.direction == UNDERCALL
    assert headline.true_phenotype == "Poor Metabolizer"
    assert headline.called_phenotype == "Intermediate Metabolizer"
    assert headline.called_label == "*2A/Reference"

    # A panel typing every defined allele must get everything right, or the
    # model itself is broken.
    every_locus = frozenset().union(*(a.rsids for a in dpyd.alleles.values()))
    assert all(c.verdict == CORRECT for c in audit(dpyd, every_locus, SAS))

    # The "full" preset is NOT that panel, and the lattice is how you find
    # out: 32 loci against 35 defined alleles. *3 and *7 are absent, both No
    # function, so 71 diplotypes come back a phenotype too mild. This is the
    # recorded 32-vs-35 drift, rediscovered from the other direction. When
    # the preset gains the two loci this fails, and that is the signal to
    # delete this check.
    full_preset = frozenset(presets()["full"][1])
    untyped = {n for n, a in dpyd.alleles.items() if a.rsids and not (a.rsids & full_preset)}
    assert untyped == {"*3", "*7"}, untyped
    wrong = [c for c in audit(dpyd, full_preset, SAS) if c.verdict != CORRECT]
    assert len(wrong) == 71, len(wrong)
    assert all(c.direction == UNDERCALL for c in wrong)

    # --- frequency -------------------------------------------------------
    # Hardy-Weinberg: p^2 for a homozygote, 2pq for a heterozygote.
    p = dpyd.alleles["c.2279C>T"].frequencies[SAS]
    assert abs(diplotype_frequency(dpyd, "c.2279C>T", "c.2279C>T", SAS) - p * p) < 1e-12
    # CPIC publishes no reference frequency for DPYD in any population, so
    # it is derived and labelled. Without this every heterozygote would go
    # unweighted, and those are the common ones.
    ref, how = dpyd.reference_frequency(SAS)
    assert ref is not None and how.startswith("derived"), (ref, how)
    assert "absorbed" in how, how
    assert 0.4 < ref < 0.7, ref
    assert abs(
        diplotype_frequency(dpyd, "Reference", "c.2279C>T", SAS) - 2 * p * ref
    ) < 1e-12

    # TPMT does publish it, and must then be used rather than derived.
    tpmt_ref, tpmt_how = tpmt.reference_frequency(SAS)
    assert tpmt_how == "published by CPIC", tpmt_how
    assert tpmt_ref and tpmt_ref > 0.9, tpmt_ref

    # An allele with no published frequency yields None, never 0.0.
    assert "M166V" in dpyd.unweighted
    assert diplotype_frequency(dpyd, "M166V", "M166V", SAS) is None

    # CPIC's own number for the variant that drives the Indian gap, kept as
    # a canary: if a CPIC refresh moves it, the reports move with it.
    assert abs(p - 0.00603) < 5e-5, p

    summary = summarise(audit(dpyd, cpic4, SAS))
    assert summary.no_call == 0, "a DPYD panel never says 'I don't know'"
    assert summary.overcall == 0, "every miss runs in the dangerous direction"
    assert summary.wrong_per_100k and summary.wrong_per_100k > 0
    assert summary.unweighted_cases > 0, "unweighted cases must be reported, not hidden"

    # --- the optimiser ---------------------------------------------------
    # The whole India thesis, as a number: one locus recovers most of it.
    sas_best = rank_additions(dpyd, cpic4, SAS)[0]
    assert sas_best.rsids == ("rs112766203",), sas_best
    assert sas_best.recovered_per_100k > 500

    # And the ranking is population-specific, which is the point of ranking
    # by carriers rather than by positions.
    eu_best = rank_additions(dpyd, cpic4, "European")[0]
    assert eu_best.rsids != sas_best.rsids, (eu_best, sas_best)

    # --- panel intake ----------------------------------------------------
    with tempfile.TemporaryDirectory() as tmp:
        good = pathlib.Path(tmp) / "assay.csv"
        good.write_text("Position,rsID,Notes\n1,rs3918290,\n2,rs67376798,\n3,BAD-ID,\n")
        panel = read_panel(good)
        assert panel.loci == {"rs3918290", "rs67376798"}, panel.loci
        # A malformed entry is reported, never silently dropped.
        assert panel.rejected == ((4, "BAD-ID"),), panel.rejected
        assert not panel.clean

        plain = pathlib.Path(tmp) / "list.txt"
        plain.write_text("# my panel\nrs3918290\nrs55886062\n\n")
        assert read_panel(plain).loci == {"rs3918290", "rs55886062"}

    # --- the other genes run at all --------------------------------------
    for gene in (tpmt, nudt15):
        star3 = "*3A" if gene.symbol == "TPMT" else "*3"
        loci = frozenset(gene.alleles[star3].rsids)
        result = summarise(audit(gene, loci, gene.populations[0]))
        assert result.total > 0 and result.correct > 0, gene.symbol
        # These genes have Uncertain/Unknown function alleles, so unlike
        # DPYD they can genuinely produce an indeterminate result.
        assert any(c.verdict == NO_CALL or c.direction for c in audit(gene, loci))

    print(f"self-check OK — {len(SUPPORTED)} genes, {len(EXCLUDED)} refused, "
          "scales, projection, lattice, frequency, optimiser, intake")
    return 0
