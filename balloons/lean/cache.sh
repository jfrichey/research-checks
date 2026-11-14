#!/usr/bin/env bash
set -euo pipefail
balloons_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [ "$#" -eq 0 ]; then
  set -- Mathlib.Data.Finset.Card Mathlib.Data.Finset.SymmDiff Mathlib.Tactic.Tauto Mathlib.Data.List.Sort Mathlib.Topology.MetricSpace.Defs Mathlib.MeasureTheory.Constructions.BorelSpace.Metric Mathlib.MeasureTheory.Constructions.BorelSpace.Order Mathlib.MeasureTheory.Measure.Prod Mathlib.MeasureTheory.Integral.IntervalIntegral.Periodic Mathlib.MeasureTheory.Constructions.Pi Mathlib.Probability.Independence.Basic Mathlib.MeasureTheory.Group.AddCircle
fi
"$balloons_root/lake.sh" build cache
export GIT_CEILING_DIRECTORIES="$balloons_root/.lake/packages"
"$balloons_root/lake.sh" env bash -c '
cd "$1"
shift
exec .lake/build/bin/cache --repo=leanprover-community/mathlib4 get "$@"
' _ "$balloons_root/.lake/packages/mathlib" "$@"
