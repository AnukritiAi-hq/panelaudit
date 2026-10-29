# panelaudit

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
git clone https://github.com/<org>/panelaudit.git
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

## Reproducing the published tables

```bash
python -m panelaudit --self-check          # every invariant
python -m analysis.equity_table            # cost-of-transfer tables
python -m panelaudit --gene DPYD --all-panels
```

## Citation

If this software contributed to your work, cite the paper and the archived
release:

```
[PAPER CITATION — TO COMPLETE ON PREPRINT DEPOSIT]
Anukriti AI. panelaudit. Zenodo. [DOI — TO COMPLETE ON FIRST RELEASE]
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
