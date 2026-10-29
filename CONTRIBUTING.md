# Contributing

Issues and pull requests are welcome. One thing makes this repository
unusual, and it is worth knowing before you open a pull request.

## This code backs a paper

Release `v1.0.0` is archived at
[10.5281/zenodo.23768369](https://doi.org/10.5281/zenodo.23768369) and is
the snapshot a published manuscript's numbers were computed from. Anyone
should be able to fetch that archive, run two commands, and get the
figures printed in the paper.

That survives changes to this repository, because the archive is
immutable. What does not survive is the claim that the *current* code
still produces those numbers.

So:

- **A change that alters any published number is not a bug fix.** It may
  well be correct, and it is still a change in result. Say so in the pull
  request, show the before and after, and expect it to land under a new
  version rather than quietly.
- CI asserts several headline rows from the paper. If your change moves
  them, CI fails on purpose. Do not update the assertions to match
  without saying what moved and why.
- The pinned CPIC snapshot in `panelaudit/cache/` is deliberate. Refresh
  it in its own commit, never alongside a logic change, so drift in
  CPIC's tables is never confused with a change in this code.

## What this package will not take

- **A new upstream data source.** CPIC only, CC0 1.0. One source and one
  licence is what makes the licensing position defensible, and CI
  enforces it. PharmVar is CC BY-SA; PharmGKB's clinical annotations
  carry a sale restriction. Neither can be added.
- **A third-party dependency.** The package is standard library only and
  CI installs nothing, so the claim stays testable.
- **A hardcoded clinical value.** Phenotype comes from CPIC's own
  diplotype rows, keyed on the unordered pair of allele functions. A
  value typed in by hand is a value that silently preserves last year's
  model.
- **Anything that turns a refusal into a guess.** CYP2D6 is refused
  because phenotype depends on copy number; G6PD because it is X-linked.
  Those are not gaps waiting to be filled with an approximation.

## Before you open a pull request

```bash
python -m panelaudit --self-check
python -m analysis.equity_table
python -m panelaudit --gene DPYD --all-panels
```

No install step. Python 3.10+.

## Reporting a problem with a result

A wrong number matters more here than a crash. If a verdict looks wrong
for a genotype you know, open an issue with the gene, the panel loci and
the diplotype, and say what you expected. A reproducible disagreement
with CPIC is the most useful thing you can send.
