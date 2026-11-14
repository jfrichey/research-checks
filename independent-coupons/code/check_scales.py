"""Diagnostic sums for p_v proportional to exp(-sqrt(v)); standard library only.

Run from any directory: python3 path/to/check_scales.py
These truncated sums corroborate scales; they are not certified error bounds.
"""

from math import exp, expm1, fsum, log, log1p, sqrt


def main():
    raw = [exp(-sqrt(v)) for v in range(12001)]
    normalizer = 1 / fsum(raw)
    probabilities = [normalizer * value for value in raw]
    print(f"alpha=1; C={normalizer:.16g}")
    print("n,R*n/t^2,late_pairs*n/t^2,independent_variance/t^3")

    for n in (10**3, 10**6, 10**9):
        t = log(n * normalizer)
        active = [p for p in probabilities if 1e-10 < n * p < 60]
        future_error = 0.0
        late_pairs = 0.0
        for i, p in enumerate(active):
            for q in active[i + 1 :]:
                combined = p + q
                future_error += (
                    2 * p * q / combined * exp(n * log1p(-combined))
                )
                late_pairs += (
                    p
                    * q
                    / (combined - p * q)
                    * exp(n * (log1p(-p) + log1p(-q)))
                )
        variance = fsum(
            v * exp(n * log1p(-p)) * (-expm1(n * log1p(-p)))
            for v, p in enumerate(probabilities)
        )
        print(
            f"{n},{future_error*n/t**2:.10g},"
            f"{late_pairs*n/t**2:.10g},{variance/t**3:.10g}",
            flush=True,
        )


if __name__ == "__main__":
    main()
