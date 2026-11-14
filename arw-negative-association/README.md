# Stationary ARW negative dependence

The [selected technical proofs](report.md) supply the longer calculations
identified in the 20-page reader. They retain exact hypotheses and equation
labels; independent human verification remains outstanding.

From the repository root, `python -B run_checks.py arw-negative-association`
rechecks normalization and positivity of the 44 saved laws, all 7,656 disjoint
event-gap polynomials, and a finite complete-graph source calculation.
`verify_saved_laws.py` does not recompute the absorption solutions.

For the independent reconstruction, run `python -B arw-negative-association/four_site.py`.
For the three-site symbolic certificates and counterexamples, run
`python -B arw-negative-association/proof_checks.py`.
These longer calculations are separate from the default suite. Additional
program options are documented by `--help` and the program docstrings.
Some large optional algebra/complete-graph modes require `python-flint` and
substantial memory; no claim is made that every optional mode was rerun for
this export. The general NA conjecture remains open.
