# panelaudit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23768368.svg)](https://doi.org/10.5281/zenodo.23768368)

**What does this panel report for someone who carries X?**

A coverage report lists the positions an assay interrogates. That is not
the question a prescriber has. This answers the one they do have, for
every genotype at once, **before any patient has been sequenced**.

```bash
python -m panelaudit --gene DPYD --panel cpic4 --population "Central/South Asian"
python -m panelaudit --gene DPYD --panel my_assay.csv --json report.json
python -m panelaudit --all-panels --gene DPYD
python -m panelaudit --self-check
```

## What it does that other tools do not

PharmCAT and every other caller is **sample-first**: a VCF goes in, one
patient's diplotype comes out. They cannot answer a question about an
assay, because they need a sample before they can say anything.

`panelaudit` is **panel-first**. It takes the list of positions an assay
types and classifies every expressible diplotype:

| | |
|---|---|
| `CORRECT` | the panel resolves it and the phenotype is right |
| `NO_CALL` | the panel cannot resolve it and says so. Honest failure. |
| `WRONG` | the panel resolves it to a different, confident phenotype |

Then it weights those by published allele frequency, per population, and
**ranks which single locus to add next by the carrier burden it recovers**.
That last step is the actionable one: not "your panel misses c.2279C>T"
but "adding rs112766203 fixes 1,176 wrong calls per 100,000 people here".

## Genes

Audited: **DPYD, TPMT, NUDT15, CYP2B6, CYP2C19, CYP2C9, CYP3A5.**

Refused, with reasons — `python -m panelaudit --genes`. CYP2D6 because
phenotype depends on copy number, which no set of rsIDs can express. G6PD
because it is X-linked. NAT2 because the engine's table records a
phenotype term in its function column. SLCO1B1, VKORC1, CYP1A2, CYP3A4
because the tables needed do not exist.

Refusing by name beats producing a confidently wrong report, which is the
failure this whole package exists to expose.

## Install

Python 3.10+. **No third-party packages**, nothing to install.

```bash
git clone https://github.com/AnukritiAi-hq/panelaudit.git
cd panelaudit
python -m panelaudit --self-check
```

No patient data, no reference genome, no VCF, no network. A pinned CPIC
snapshot ships in `panelaudit/cache/`; `--refresh` re-fetches it from CPIC
if you want current tables instead of the ones the paper used.

## Data and licence

Code is Apache-2.0. Data is CC0 1.0. See `LICENSE` and `LICENSE-DATA`.

CPIC only, fetched from `api.cpicpgx.org/v1`. One source, one licence, no
ShareAlike, no sale restriction. No PharmVar, no DPWG, no PharmGKB clinical
annotations, no Aldy.

Allele definitions in `panelaudit/data/` are regenerated verbatim from the
same CPIC API. Every report prints the CPIC table version, endpoint and
fetch date it was computed against. Panel
presets in `panelaudit/presets/panels.json` are published specifications —
CPIC's four-variant DPYD set and the AMP/ACMG 2024 joint consensus Tier 1
set (PMID 39032821) — each recording its source. Nothing clinical is
re-typed by hand anywhere in this package.

## Design rules

**Nothing clinical is hardcoded.** Phenotype is keyed on the unordered
pair of allele functions, taken verbatim from CPIC's diplotype rows, and
asserted unambiguous at load. A CPIC update that splits a mapping fails
loudly instead of silently preserving last year's model.

**Phenotype is per gene, never global.** `Decreased/Normal` is
*Intermediate Metabolizer* in DPYD and *Normal Metabolizer* in TPMT. One
global table is correct for exactly one gene.

**An unknown term raises.** A new CPIC phenotype silently ranked as Normal
would turn a real finding into a clean report.

**Unknown frequency is `None`, never `0.0`,** and the count of unweighted
cases is printed with every rate. An allele CPIC gives no frequency for is
not an allele nobody carries.

**Derived numbers say they were derived.** CPIC publishes no reference
frequency for DPYD in any population, so it is computed as
`1 - sum(published variant frequencies)` and the report states that,
including how many alleles the derivation absorbs.

**Hardy-Weinberg is labelled as an approximation.** Real populations are
structured and consanguinity raises homozygote counts above HW in some of
the groups this matters most for. It is right for ranking panels against
each other, which is what it is used for, and it is never presented as a
prediction for an individual.

## Modules

| | |
|---|---|
| `cpic.py` | fetch and cache CPIC tables; the only network code |
| `model.py` | gene definitions, phenotype mapping, frequencies |
| `lattice.py` | classification, weighting, the addition optimiser |
| `burden.py` | share of a population's actionable allele frequency a panel reaches |
| `panels.py` | customer panel files and named presets |
| `report.py` | JSON and text rendering, provenance-stamped |
| `cli.py` | entry point |
| `selfcheck.py` | assertions — `python -m panelaudit --self-check` |
| `data/` | per-gene allele tables, CPIC-derived, CC0 |
| `presets/` | named panel specifications, with sources |

## Status

Investigational. Describes an assay, not a person. Computes no dose, names
no patient, asserts nothing clinical CPIC does not publish. Carries no
clinician sign-off.

## Licence stability

This package stays under Apache-2.0 with its data under CC0 1.0, and will
not be relicensed to more restrictive terms.

That is a commitment, not a default. The accompanying paper argues that
deployed genotyping panels misclassify phenotypes, and its author has a
commercial interest in that conclusion. The competing-interests
declaration rests on the fact that this software is free to obtain,
inspect, modify and use — including commercially — so that no reader has
to take the result on trust. Narrowing the licence later would retro-
actively undercut a published declaration.

A separate engine, `anukriti-pgx-core`, is licensed differently. It is not
used here and is not needed to reproduce any published result.

## Reproducing the published tables

```bash
python -m panelaudit --self-check          # every invariant
python -m analysis.equity_table            # cost-of-transfer tables
python -m panelaudit --gene DPYD --all-panels
python -m panelaudit --burden --panel cpic4   # coverage by population
```

### Actionable-burden coverage

`--burden` answers the question an external quality scheme asks: of all the
allele frequency CPIC says to act on, what share sits at a position this
panel actually types? It is a property of the panel and the population, and
unlike the wrong-phenotype rate it needs no diplotype model.

```
  Central/South Asian                78.1%   (2.538% of 3.250% allele frequency)
  European                           98.7%   (3.595% of 3.641% allele frequency)
```

The whole South Asian gap is one allele. CPIC prices `c.2279C>T` at 0.603%
there and at **0.0** in Europeans, which is why a panel assembled from
European data omits it. Both figures are pinned in `--self-check`, because
they have been quoted outside this repository.

The starker case is African. For Sub-Saharan Africans the same panel
reaches **0.0%** of the actionable frequency CPIC publishes, because
`c.557A>G` alone is 2.59% of a 2.79% burden and `cpic4` does not carry it —
even though the FDA capecitabine label names that variant. The gaps section
ranks populations by frequency unseen rather than by percentage, so this
case leads rather than a population with a small burden the panel also
misses.

An allele counts toward the burden when CPIC's clinical function is
decreased, none, or a possible-\* form of either. Normal, **increased** and
uncertain are excluded, and `burden.py` argues the increased-function case
rather than assuming it: CPIC does act on increased function, but a panel
missing one is a different failure with a different consequence, and DPYD
has no such allele. Classification is total — a CPIC term this build does
not recognise is reported, not rounded into whichever branch an
if-statement ended on. Alleles CPIC publishes no frequency for are listed
by name and excluded from both halves of the fraction rather than counted
as zero.

## Citation

If this software contributed to your work, cite the paper and the archived
release:

```
R B A. panelaudit: panel-first pharmacogenomic audit. Zenodo; 2026.
doi:10.5281/zenodo.23768368
```

The DOI above is the **concept DOI** and always resolves to the latest
version — use it unless you need to pin. The exact snapshot behind the
accompanying preprint is v1.0.0, `doi:10.5281/zenodo.23768369`.

```
[PAPER CITATION — TO COMPLETE ON PREPRINT DEPOSIT]
```

And cite the upstream data as CPIC asks:

> CPIC (Clinical Pharmacogenetics Implementation Consortium),
> www.clinpgx.org.

## Who made this

Built at [Anukriti AI](mailto:abhi@anukritiai.com), which reviews
pharmacogenomic assays and reporting workflows for laboratories and CRO
central-lab teams.

This package describes an assay. It does not describe a person. If you need
the other direction — a real genome called, the uncertain parts refused by
name, and a report that discloses which positions were observed rather than
assumed — that is the service, not this repository. Get in touch.

## Which frequencies weight the rates

Every rate here is a frequency-weighted count, so it inherits whatever
allele frequencies it was weighted by. Two sources are available:

```
python -m panelaudit --population "European"                        # CPIC (default)
python -m panelaudit --population "European" --frequencies gnomad   # gnomAD v4
python -m panelaudit --compare-frequencies                          # where they differ
```

They differ more than you would hope. For DPYD, CPIC overstates `*2A` — the
most clinically used allele in the gene — by 4x in Central/South Asian, 4x in
European and 12x in African American, and records `c.868A>G` as absent in
African ancestry where gnomAD measures 0.23%. Four DPYD alleles share an
identical African-American frequency of 0.0031083, which is 1/321.7: one
observed allele in about 161 genotyped people, a singleton floor rather than
a measurement.

Neither column is the truth. CPIC curates from the pharmacogenomic literature
and is what the guideline rests on; gnomAD is a far larger convenience sample
with its own ascertainment. The useful output is the disagreement, because a
panel decision that flips between the two is resting on the choice of
reference rather than on the assay.

Under `--frequencies gnomad` some weights are deliberately absent:

- **Populations gnomAD cannot answer** are dropped, not approximated. CPIC's
  Sub-Saharan African has no counterpart, because gnomAD's `afr` pools
  African American with continental African and is already in use for CPIC's
  African American/Afro-Caribbean.
- **Haplotypes** have no single position, so HapB3 carries no gnomAD weight.
- **Components of a haplotype** carry none either. `c.1236G>A` is the HapB3
  tag and does have a position-level frequency in gnomAD; weighting it
  independently would count HapB3 carriers twice and treat a tag as causal.

Alleles left unweighted are reported, never silently zeroed.
