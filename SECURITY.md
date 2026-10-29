# Security

## Scope

This package reads published reference tables and prints a report. It
handles no patient data, no genome, no credentials and no network input
beyond CPIC's public API, and it runs nothing it downloads. The realistic
security surface is small.

What would still be worth reporting:

- A way to make the package execute or write something unexpected from a
  crafted panel file, since a panel file is the one input that comes from
  outside.
- Anything in this repository that looks like a credential, a token, or
  data that should not be public.
- A supply-chain problem in how a release is built or archived.

## A wrong answer is not a security issue, and matters more

If the tool reports a verdict that disagrees with CPIC, that is a
correctness bug, and it is the failure mode this package exists to avoid
in other people's pipelines. Open a normal issue for it, publicly, with
the gene, panel loci and diplotype. It will be taken seriously.

## How to report

Email abhi@anukritiai.com. Please allow a reasonable period to respond
before disclosing publicly. There is no bounty.

## What this software is not

Investigational. It describes an assay, not a person. It computes no dose
and asserts nothing clinical that CPIC does not publish, and it carries no
clinician sign-off. Do not use it to make a prescribing decision.
