# Supporting arguments and checks

| Claim in the reader | Companion | Verification scope |
| --- | --- | --- |
| Entropy comparisons, row/layer certificates, and explicit swap conjugacies | [Technical arguments](technical-details.md), using the named theorem/heading in the reader | Ordinary working proofs and exact finite certificates; the general equal-overlap conjugacy question remains open. |

The reader’s swap example and agreement-class corollary are checked by `swap_conjugacy_checks.py` in the default suite. Optional entropy/certificate checks use `minimum_certificates.py`, `alphabet_sweep_checks.py`, `layer_strictness_checks.py`, `ternary_layer_checks.py` and `row_entropy_checks.py`. Run `python3 prepare_certificates.py` first to restore the four losslessly compressed saved JSON inputs, whose byte lengths and SHA-256 hashes are checked against `certificate-inputs.json`. This editorial release has not rerun every large enumeration.

The pinned paper version and PDF hash are recorded in [checks.json](../checks.json).
