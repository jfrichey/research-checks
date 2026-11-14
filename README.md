# Research checks

Selected mathematical checks and finite Lean proofs accompanying
[Jacob Richey's research](https://jfrichey.github.io/ai.html).
The code was developed with ChatGPT; source attribution is retained.

| Package | What the download supports |
| --- | --- |
| [ARW monotonicity](arw-monotonicity/) | Theorem 8.1: six exact polynomial identities and 2,232 nonnegative Bernstein coefficients |
| [ARW negative association](arw-negative-association/) | Selected technical proofs, exact saved laws, and reproducible finite checks |
| [Word potentials](word-potentials/) | Symbolic identities and floating-point diagnostics, distinguished in the code |
| [Higher-dimensional SFTs](higher-dimensional-sfts/) | Explicit swap rules and finite agreement-class calculations |
| [Counting clusters](counting-clusters/) | Exact small-strip checks and a separate program for the width-17 enclosure |
| [Combinatorial processes](combinatorial-processes/) | Selected exact, symbolic and numerical calculations in the illustrated notes |
| [Independent coupons](independent-coupons/) | Truncated numerical sums for the stretched-exponential example |
| [Sandpile](sandpile/) | An auxiliary three-site rational obstruction |
| [Balloons](balloons/) | Exact interval computations and a [Lean package](balloons/lean/) for the finite process |

Each package is limited to the stated results. Passing a numerical or symbolic
check does not formally verify the surrounding analytic paper. The balloon
entrance-law and asymptotic noise/history theorems are not Lean-verified.

## Run the default checks

With Python 3.12 and a C++17 compiler (`g++` for the small cluster check):

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-checks.txt
python -B run_checks.py --list
python -B run_checks.py
```

On Windows, activate `.venv\Scripts\activate` and provide a compatible C++
compiler, or run the Python-only projects individually. For example:
`python -B run_checks.py arw-monotonicity`.
Do not use Python's `-O` option: these mathematical programs use assertions.
The runner checks the release's file hashes before executing its explicit list.

The longer enclosure, interval and Lean computations are documented in their
package READMEs. They are separate from the default Python suite. Some optional
large symbolic modes of the ARW programs require `python-flint`; those modes
are not part of the default suite.

## Papers and versions

`checks.json` records each paper's exact PDF hash, versioned download link,
verification commands and scope. `files.sha256.json` identifies the exported
files. The website links to fixed code revisions. A later working-file change
is published only through a new reviewed export.

This repository is a selected public release. It contains no private repository
history, correspondence, agent instructions, continuation notes or build caches.
Readers can browse and download it. Changes to the maintained repository require
write permission from its owner. No additional license grant is included.
