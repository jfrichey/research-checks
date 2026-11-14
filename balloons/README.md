# Balloon-process companions

[Proof and verification map](CLAIMS.md).

The default suite runs `certify_cutoff_budget.py`, an exact rational evaluation
of the early-label error bound. It does not estimate the limiting noise probability.

The separate interval computations can be reproduced from this directory:

```sh
python3 -B certify_noise_constants.py
python3 -B certify_matching_spectrum.py
```

They use outward integer intervals and explicit quadrature/tail error bounds.
The first evaluates the paper's stated rigorous lower/upper bounds; the second
evaluates the matching-spectrum coefficient pi_2. Neither certifies the proposed
value near 0.658. Saved JSON certificates identify their parameters and scope.

The [Lean development](lean/) proves named finite-process statements through
`Stable.lean`. Its README describes the build and its exact limits. The later
unfinished formal extensions are not included in this release.
