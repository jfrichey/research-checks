# Cluster-counting calculations

[Proof and verification map](CLAIMS.md).

The default suite checks the small-strip and all-length formula certificates,
finite Euler identities, and a width-four C++ residual calculation.

To reconstruct the paper's displayed width-17 enclosure, run from the repository root:

```sh
python -B counting-clusters/code/verify_enclosure.py
```

This separately compiles the C++17 solver into a temporary directory and computes
four residual certificates. It needs several GB of memory and many minutes.
The default suite does not claim to repeat that larger calculation. Its recorded
claim is in `data/enclosure.json`; no extrapolated bound is substituted.
