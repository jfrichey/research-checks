# Lean verification of the balloon process

The project uses Lean 4.35.0 and mathlib v4.35.0, pinned to mathlib commit
`56a005855ea17e0b0958f73ec801ae140b0a1679`. It formalizes the finite event scan, its metric realization, and the first
exact probabilistic identities.

## Reproduce

On Linux x86-64, with Python 3.12, curl, tar, and zstd:

```sh
python3 bootstrap.py
./cache.sh
python3 check.py
```

The bootstrap downloads archives with recorded SHA256 checksums. Dependencies
and the compiler stay in `.lake/` and `.toolchain/`; archives stay in
`/tmp/balloons-lean-bootstrap`. No Git mutation commands are needed. The
wrapper `lake.sh` selects this project's compiler without changing the user's
shell or home directory. Additional mathlib modules can be fetched with
`./cache.sh Mathlib.Module.Name`. If the default cache endpoint stalls,
mathlib provides the official Azure fallback:
`MATHLIB_CACHE_DEBUG_USE_LEGACY=1 ./cache.sh`.

`check.py` builds the complete library and prints the transitive axiom list
for every public project theorem. It rejects proof placeholders, added axiom
declarations, and native decision shortcuts in the project sources; the
reported axioms must be among `propext`, `Classical.choice`, and `Quot.sound`.
The generated report is `.lake/axiom-audit.txt`.

## Checked scope

This release selects the completed development through `Stable.lean`: 25 modules
and 193 public theorems. Later unfinished extensions are excluded. The build and
axiom audit passed again for this exported selection on November 14, 2026.
The [verification record](verification.json) and [per-theorem axiom report](axiom-audit.txt)
identify that check.

The checked finite foundation includes:

- `Finite.lean`: event scan; exact zero-or-two discrepancy decrease;
  contraction and parity; absorbing agreement; exactly one discrepancy for
  an insertion; the `2*k` replacement bound.
- `Terminal.lean`: a complete schedule leaves one survivor for odd input
  size and none for even input size.
- `Schedule.lean`: an explicit complete pair enumeration, its sorted
  schedules, and contraction between arbitrary prefixes.
- `Relabel.lean`: the scan and discrepancy counts commute with injective
  relabeling.

These are proofs for arbitrary finite configurations and event schedules,
not enumerations of small examples. Schedule completeness and input parity
are explicit hypotheses where needed. The enumeration uses both orientations
of each pair; an already scanned pair cannot remove vertices a second time.

The following additional layers are checked:

- `History.lean` and `HistoryDistance.lean`: disjoint matching edges, exact
  vertex accounting, and `H = N - 1 - 2*J` for equal odd complete inputs.
- `Time.lean`: sorted priority cutoffs, exact event membership, and contraction
  between cutoffs.
- `Origin.lean`: pointwise discrepancy transport, its composition and
  relabeling, and the survivor after insertion into an even background.
- `Restriction.lean`: events outside an initial ambient set can be discarded
  without changing its configuration or history.
- `Metric.lean`: actual distance ordering, separation of survivors, zero-cutoff
  initialization, terminal completion, and edit bounds at metric cutoffs.
- `Measurable.lean`: measurable labeled configurations, histories, and origin
  transport, including measurable random cutoffs. The proof uses the finite
  table of priority comparisons. Measurable locations in a second-countable
  metric Borel space supply the required distance priorities.

- `Basins.lean`: independent diffuse origins agree with mass equal to the sum
  of squared old-vertex basin masses, including the actual even-background
  insertion survivors after removing a null coincidence event.
- `Insertion.lean`: an explicit distance-sorted insertion map with a new
  label and fixed old locations. Its spatial map is measurable, its output is
  the survivor for even backgrounds, and the basin identity applies without
  an assumed measurability premise. Joint measurability and integration over
  an independent random background give the expected basin-square formula.
- `Circle.lean`: normalized uniform circle sampling, diffuseness, completion
  by separation `1/2`, and the exact finite spatial-noise identity
  `p_N = 1 - E(sum basinMass^2)` for `N = 2*m+1`.
- `Sampling.lean`: product laws of distinct iid coordinates and almost-sure
  injectivity of a finite diffuse sample.
- `Ties.lean`: random-sphere avoidance, circle spheres of measure zero, and
  almost-sure absence of distance ties between distinct unordered pairs.

- `Invariance.lean`: sorted schedules representing the same events have the
  same configurations and spatial histories at every cutoff when ties only
  reverse one unordered pair.
- `Enumeration.lean` and `SortedRelabel.lean`: complete pair enumeration has
  no repetitions; sorting the actual process commutes with injective
  relabeling under the proved no-tie condition.

- `Winner.lean` and `UniformWinner.lean`: the unique survivor label is
  measurable and equivariant, and each label has survival probability `1/N`
  for odd iid uniform circle input.
- `InsertionLaw.lean`: the augmented input has the iid product law, its
  insertion label equals the unique winner, and the refreshed label survives
  with probability `1/N`.

- `LabelNoise.lean`: for every odd circle sample size, including one point,
  `p_N^lab ≤ p_N^pos ≤ p_N^lab + 1/N` after one resampling.
- `LabelNoiseExact.lean`: measurable available mass `Theta`, independent-origin
  square identity, and the exact correction
  `p_N^pos = p_N^lab + E[Theta^2]`, hence
  `p_N^lab = 1 - E[sum_old basinMass^2] - E[Theta^2]`.

- `Stable.lean`: the greedy history has no blocking pair; distinct unordered
  pair distances give strictly shorter blocking edges and uniqueness of the
  stable perfect matching on every even finite metric input.

The circle definition uses a fixed label enumeration, with its irrelevance
proved under the almost-sure no-tie condition. The marked entrance limit,
survivor-noise limit, sharp history asymptotics, and fluctuation laws
are not yet formally verified. The [reader article](https://jfrichey.github.io/pagedocs/ai/balloons-sketch.pdf)
explains the informal dependency chain. A successful build certifies only the
statements in these modules, with the hypotheses present in their Lean types.
