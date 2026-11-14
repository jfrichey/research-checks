# Some combinatorial processes

[Proof and verification map](CLAIMS.md).

Run `python -B run_checks.py combinatorial-processes` from the repository root.
The three suites mix exact finite enumeration, symbolic identities, numerical
integration and asymptotic diagnostics. Saved data and figure/source hashes are
included because the extension suite verifies them. The separate arrow simulator
is outside this export. Passing these checks is not a proof of all asymptotic
statements in the 49-page notes. Figure generation uses the optional dependency
list in `code/requirements.txt`; it is not part of the verification command.
