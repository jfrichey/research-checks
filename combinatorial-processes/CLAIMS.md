# Supporting arguments and checks

| Claim in the reader | Companion | Verification scope |
| --- | --- | --- |
| Selected identities, finite censuses and diagnostic plots in the illustrated notes | [Result-to-file map](README.md), `code/check_formulas.py`, `code/check_illustrated.py`, `code/check_extensions.py` | Exact combinatorial and symbolic assertions are distinguished from floating-point integrals and asymptotic diagnostics. |

The paper keeps its mathematical arguments. The Egyptian-fraction finite search gives exact minima within its stated complete search; finite data alone do not prove asymptotic claims.

The pinned paper version and PDF hash are recorded in [checks.json](../checks.json).

## Functions by topic

| Topic in the notes | Exact location |
| --- | --- |
| Sock paths and maxima | `code/check_formulas.py:socks`; `code/check_extensions.py:sock_paths` |
| Dishes and permutations | `code/check_illustrated.py:dish_balance, permutation_masses` |
| Coupon and seating laws | `code/check_formulas.py:coupon_finite, permutations_and_seats` |
| Adsorption and its finite distribution | `code/check_extensions.py:adsorption` |
| Arrows, pursuit and circle covering | `code/check_formulas.py:arrows, circle`; `code/check_illustrated.py:arrow_drift, pursuit, circle_targets` |
| Uniform renewal and rapid records | `code/check_extensions.py:uniform_renewal, rapid_records` |
| Egyptian fractions and minimum lengths | `code/check_extensions.py:egyptian_checks`; `code/egyptian_fractions.py:shortest` |
| Saved data and figure provenance | `code/check_extensions.py:manifests` |

Each function's assertions state its finite domain or tolerance. A check of a stored hash establishes file identity, not the mathematical validity of its contents.
