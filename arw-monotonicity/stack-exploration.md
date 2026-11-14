# Stack exploration and coupled history switches

Updated 2026-10-10 UTC. The finite all-active single-target occupation
conjecture remains **open**, with the kernel, accessible-sink and
`0 < p_z < 1` assumptions in [the reader’s model and hypotheses](https://raw.githubusercontent.com/jfrichey/jfrichey.github.io/companions-2026-11-14/pagedocs/ai/arw-monotonicity.pdf).

**Local results proved below:** fixed record permutations repair BS/R's
inventory examples; exchanging stopped explorer words gives an adaptive,
invertible repair, including a two-common-particle family where an explorer
wakes a common sleeper before stopping, and a loop-erased second switch repairs
one resulting lower wake. Fixed root-word reversal and a recoverable unused
cycle extend this to a larger two-sleep-endpoint family.
[Complete avalanche exchange](avalanche-exchange.md) handles arbitrary
common-sleeper contacts on its matched-background domain and recovers both
adaptively selected terminal target indices; reserved unused sleeps permit
some unmatched backgrounds. A complete drain gives a common product space
for arrow books and labelled random-walk paths, with a full-book lift of
causal contact exchanges under joint success conditioning. Separate success
conditioning fails. Sharing a common complete FIFO path also fails, even
with one common particle. A component involution cancels disconnected
histories. The general sign and global repair selection remain open.

## 1. Fixed record permutations and their weights

Use R §1's buffer representation. The m-th arrival at z has a record
`R_z(m) = (B_z(m), D_z(m))`, with independent Bernoulli `B_z(m)` of parameter
`p_z` and jump destination `D_z(m)` of law `K(z,.)`. An arrival first emits
the previous resident's waiting arrow, if there is a resident; it then
stays when `B=1`, or emits its own arrow when `B=0`. Pending deliveries
may be postponed as in the saved Abelian buffer construction. This is
not permission to ignore a physically present particle in a raw sleep test.

**Lemma.** A fixed finite permutation of record coordinates between two
copies, preserving the site of every record, preserves the joint product
law, including unused records. Its inverse is the inverse permutation.

**Proof.** Each coordinate at z has weight

`w_z(b,d) = p_z^b (1-p_z)^(1-b) K(z,d)`.

On every finite cylinder containing the permuted coordinates, the product
of these factors is unchanged. The remaining coordinates are untouched.
Cylinder identities establish preservation of the full product measure.
For disjoint transpositions the map is an involution. This argument makes
no symmetry assumption on K or equality assumption on the different p_z.

This also identifies the sleep-odds accounting. If a switch changes sleep
bits `(1,0)` at sites `(z,w)` into `(0,1)`, its **coin-weight** ratio is
`ell_w/ell_z`, where `ell_z=p_z/(1-p_z)`. Jump factors require separate
accounting. In the same-site transpositions below the full ratio is one.

The word *fixed* matters. If addresses are chosen from completed histories,
the resulting piecewise permutation needs disjoint image cylinders and
an inverse that recovers the selected addresses. The lemma does not
establish either property for a history-dependent choice.

## 2. An explicit inverse for the saved inventory examples

Write `-` and `+` for the two copies, take target 0, a common active at 0,
and an extra particle at 1 in the upper copy. First consider an initially
sleeping extra, with waiting arrow `D_1^+(0)`. In the lower copy the first
record at 0 sleeps; in the upper copy it jumps to 1. Let `c=1` denote the
upper arrival index of that first contact.

Apply exactly these two transpositions of **whole potential-arrow records**:

`R_0^-(1) <-> R_0^+(1)`,

`R_1^-(1) <-> R_1^+(c)`.

The lower record at 1 was previously unused. Keep every other record at
its original address, whether it remains used or becomes unused. The
inverse is the same pair of transpositions; no supplementary random mark
or rejection choice is needed.

For the all-active version, process the upper extra at 1 first and restrict
to its initial sleep branch. Its creator record is `R_1^+(1)`, unchanged
by the switch, and now `c=2`. Its sleep factor `p_1` remains in the cylinder
weight. Thus these examples retain an actual creation law from all-active
input, rather than an arbitrary sleeping background.

**BS §5.** The contact record has bit 0 and arrow to the sink, and the
extra resident's arrow also goes to the sink. Before the switch the target
bits are `(1,0)`. Afterwards the upper common particle sleeps immediately
at 0, leaving the extra at 1 untouched. The lower common particle jumps
to 1 and then exits, so the target bits are `(0,1)`. The extra resident's
formerly used sink arrow becomes unused at its unchanged address.

Summing the independent lower padding, the sleeping-extra cylinder weight is

`p_0 (1-p_0) (1-p_1) K(0,1) K(1,sink)^2`.

The image has exactly this weight: one factor `K(1,sink)` now specifies
an unused arrow. The all-active version has the additional factor `p_1`.

**R §1.** On the path `sink—0—1—2—sink`, take the resident arrow at 1 to
be `1->2`, the contact record at 1 to be `(1,2)`, and the visited record
at 2 to be `(0,sink)`. Initially the upper explorer sleeps at 1 while
the resident travels through 2 to the sink; target 0 is empty. After the
same two transpositions, the lower explorer sleeps at 1, the upper common
particle sleeps at 0, and the upper extra remains at 1. The entire record
at 2 is now unused, still stored at its original upper coordinate.

Thus variable amounts of consumed data becoming unused cause no problem
for a specified cylinder permutation. This verifies the proposed operation
on the two obstructions; it does not specify a rule for choosing those
coordinates in an arbitrary history. The repairs themselves were described
in BS/R; the explicit whole-record involution supplies their decoding.

**Common-sleeper extension.** Before these root arrivals, allow any number
of other common particles to stabilize separately in the two copies without
visiting `{0,1}`. Their creation records remain fixed. In the R pattern the
upper explorer sleeps at its first contact at 1; the extra resident then
starts an arbitrary finite avalanche. Restrict only to its final target bit
being zero. This avalanche may wake common sleepers and return to 0 and 1.
After the same fixed transpositions, neither copy starts that avalanche:
the lower root sleeps at 1, and the upper root sleeps at 0 with the extra
untouched. The respective common backgrounds remain stable. Thus the target
bits become `(0,1)`. All later upper records become unused at unchanged
addresses. One fixed involution applies to the entire family, so there is
no length-dependent selection or overlap problem.

For example, with common input `{0,2}` and extra 1, the upper extra wakes
the common sleeper at 2, which sends a packet back to 0; another packet
returns through 1 and 2. The verifier gives original terminal supports
`({0,2},empty)` and switched supports `({1,2},{0,1,2})`. The earlier
restriction excluding returns and other common sleepers was unnecessary
after the explorer's stopping point.

### A recoverable switch of variable-length explorer words

Here is a different involution whose addresses are selected during exploration.
For two independent fresh record fields, start an explorer at the same site a
in each. Follow its own arrow when `B=0`, stopping at its first `B=1` or sink
exit. Retain the whole record, including the unused arrow at a sleep. At z,
write each field as `W_z^- T_z^-` or `W_z^+ T_z^+`, where W is the finite
word read by this explorer and T is the remaining suffix. Define

`W_z^- T_z^-, W_z^+ T_z^+  ->  W_z^+ T_z^-, W_z^- T_z^+`.

**Lemma.** This transformation preserves the joint product law and is an
involution, including with returns and unequal word lengths at the same site.

**Proof.** Following an explorer in the transformed minus field reproduces
the old plus path, and conversely, so both cuts are recovered without guessing.
The remaining tails retain their copy. The two word-cylinder weights have
the same product after exchange: every complete record retains its site.
The stopping rule reads only used prefixes, leaving independent fresh tails.
The countable partition into finite stopped paths therefore gives a
measure-preserving involution on complete fields. Finiteness follows, for
example, from `min_z p_z>0`. This proof also holds conditionally on previously
completed creation histories, using suffixes after their last read records.

**All-active repair family.** Write the common binary input as `C+delta_a`
and add an upper particle at y outside that support; fix target `v!=y`.
First stabilize C in both copies, leaving sleeping supports `B^-`, `B^+`.
Require their union to avoid `{a,y,v}`. In the upper copy deliver y next
and retain its immediate-sleep branch, including its creator record and
factor `p_y`. These finite histories fix the two suffix offsets; they are
not changed or renormalized by the map.

Now explore the common root from a in each copy. Hold any resident's emitted
packet until this explorer stops; this is a legal pending-delivery schedule.
Suppose the lower explorer avoids `B^- union B^+ union {y}` and sleeps at v;
the upper explorer avoids `B^- union B^+` and stops at a different site or
the sink. It may visit y and release the extra resident. Complete the upper
avalanche arbitrarily, and restrict to final target bit zero.

Exchange the stopped explorer words above, below the fixed creator prefixes.
The new lower explorer follows the old upper path through vacant buffers,
ending away from v. The new upper explorer follows the old lower path,
avoids y and both backgrounds, and sleeps at v. No released packet remains
in either image. Hence `(1,0)` becomes `(0,1)`. The later original avalanche,
including any common-sleeper wakes and returns, is stored in its shifted,
unused upper tails. The same stopped-walk rule recovers the inverse. Since
the full involution preserves measure, restricting it to these bad histories
gives an injection of exactly equal mass; no favorable-path conditioning is
being assigned its unconditioned probability.

The restriction concerns contacts **before** the explorers stop. A prior
common-sleeper contact releases a common resident in one or both images;
those packets can later destroy the transplanted target sleep. Their creator
histories must then participate in a further switch. Another uncovered case
has both initial explorer paths sleep at v before the upper avalanche empties
it. Neither case is settled by this local lemma.

### Creator words and a loop-erased continuation

The canonical proof is now in [creator-exchange.md](creator-exchange.md).
With common input `{v,a}` and extra y, exchange the lower creator's stopped
word from v with the upper explorer's stopped word after its first contact
at v. Earlier upper records are retained, so the contact and both cuts are
recovered even with repeated visits. This handles a common sleeper woken
before the explorer stops, under explicit residual-avoidance conditions.

A loop-erased switch and a fixed reversal of the lower root-word order each
repair families of resulting lower wakes. The latter allows arbitrary
returns of a root to its eventual sleep site when the other word avoids it.
For mutual endpoint visits, a further map moves a closed jump word into the
unused books; its first-return decoder recovers every record. This extends
to a full cycle insertion/removal involution Omega.

One fixed composition, root-order reversal after Omega after the creator
switch, handles an exiting second root and a two-sleeper family: both roots
first visit their distinct sleep endpoints at their terminal arrivals, and
the second root visits the first endpoint. Nine full pair witnesses check
the maps and reverse-order inverses; the three cycle witnesses include
repeated interior visits. General endpoint returns, upper residual returns,
and overlaps with other repair families remain open. The detailed domain
and cycle decoder are in CX §§5–6; individual failed switches are no longer
being treated as obstructions to all the available moves.

CX §7 also recovers an upper **last**-target-arrival cut when the lower
creator is a single immediate sleep. The upper word may wake other sleepers:
AX §1 compares the resulting forced resident releases. The upper last index
remains unchanged and supplies the
inverse, allowing earlier upper returns. Three full pairs verify this
extension where the first-contact switch fails. More general lower-only
repair is impossible without the upper conditions: CX gives an exact
three-site source family larger than all possible empty-target images.
In a separate marked variant, a preceding upper sleep recovers the cut
even when a returning lower creator increases the final target index.
The lower word's contacts with old sleepers must be contained in the
upper word's contacts; avoiding all such sleepers is sufficient.
Three further pairs verify the marker and inverse, including an off-site
creator excursion. The two cut rules still need disjoint image accounting
before they can be combined over their full domains.

### What a last-arrival selection does and does not preserve

There is a useful one-sided deterministic fact. Suppose a completed upper
history visits v exactly m times and ends with `B_v(m)=0`. Replacing its last
record `(0,d)` by `(1,d')` leaves its number of visits to v equal to m and
makes v occupied, for any allowed d'. Up to that last arrival the executions
agree. Afterwards the new execution suppresses one emitted packet. With v
held, monotonicity of the outside buffer output prefixes bounds its remaining
arrivals by those of the old execution. The latter never revisits v. Hence
neither does the new one, and its new potential arrow d' remains unused.

Consequently, for a **fixed** lower target index n, exchanging a lower
terminal sleep at n with the upper terminal jump at its adaptively discovered
index m has a recoverable inverse. The image still has upper terminal index
m, while n is fixed in advance. It preserves the cylinder weight, but the
lower copy may return to v and sleep again, giving `(1,1)` instead of `(0,1)`.

Making n the lower history's own last-arrival index loses injectivity.
An exact witness on `sink—0—1—2—sink`, with common input `{0,2}`, extra 1,
and target 0, has two equal-weight bad sources, with terminal-index pairs
`(1,1)` and `(2,1)`, mapping to the same good field. Its target records are

| Field | Lower `R_0(1), R_0(2)` | Upper `R_0(1)` |
|---|---|---|
| Good image | `(0,1), (0,sink)` | `(1,1)` |
| First bad source | `(1,1), (0,sink)` | `(0,1)` |
| Second bad source | `(0,1), (1,1)` | `(0,sink)` |

All other records agree and are specified in `check_terminal_overlap` in
the verifier. These are complete padded cylinders, with no truncation or
simulation inference. Thus last-arrival transpositions cannot simply be
combined over n; the returning segment needs additional reversible handling.
This is a failure of that selection rule, not of occupation monotonicity.

[AX §4](avalanche-exchange.md) supplies the additional reversible handling
on its matched-background domain. After the terminal transposition it
exchanges complete auxiliary avalanches, including initial resident arrows.
The lower image then makes no further target visits, so both original
indices are recovered from the image. A counterfactual lower sleep recovers
the initial background and offsets. The two sources above now have distinct
images, and both full inverses have been checked.

## 3. Separate success conditioning forbids some cross-copy swaps

ST §6 exchanges colours **within one copy** conditioned on at least one
successful extended jump path. Applying the same exchange between two
separately conditioned copies is false.

**Exact obstruction.** Let y jump to target v or the sink with probabilities
1/2 each, and let v jump to the sink. Sleep probabilities may be any values
strictly between zero and one. All particles in this example start at y,
so the two selected labels are co-located before any target visit.

With one particle in the lower copy and two in the upper, the lower colour
is necessarily successful. The upper colour pairs `(1,0),(0,1),(1,1)` each
have probability 1/3 under its own success conditioning. Swapping the lower
colour/path with the first upper colour/path sends the case `(0,1)` outside
the conditional support: the lower copy then has no successful path.

With two lower particles and three upper particles, the corresponding
escape probability is **4/21**. The 3 lower nonempty masks and 7 upper
nonempty masks give 21 equiprobable pairs. Three fail because the lower
selected label carries its sole success and the upper selected label is
unsuccessful; one fails with the two copies reversed.

This refutes only the cross-copy extension under **separate** success
conditioning. The original within-copy ST lemma remains valid, and this
is not an occupation counterexample.

## 4. Joint success conditioning permits cross-copy contact exchange

**Lemma.** Give all labelled particles in both copies independent extended
K-walks, continued to the sink regardless of sleep. Let h(x) be the probability
that such a walk from x visits v. Let E be the event that at least one walk
in the combined collection visits v. At a joint colour-blind causal stopping
time before either copy visits v, the unused paths and success colours of
two labels at the same site can be exchanged, including across copies,
under the conditional law given E.

**Proof.** For current positions `x_i`, put

`F_pair = 1 - product_i(1-h(x_i))`,

where the product runs over both copies. Conditional on the revealed past
and E, a nonempty success set J has probability

`[product_(i in J) h(x_i) product_(i notin J)(1-h(x_i))] / F_pair`.

The likelihood-telescoping argument of ST §5 applies to all labels together:
sleep factors and colour-blind scheduler factors cancel. At a contact the
two h-values agree, this mask law is exchangeable in the two labels, and
their independent conditioned unused walks start at the same vertex.
Exchanging colour and whole unused walk together therefore preserves
the joint conditional law. A copy may lose its last success colour; E is
still satisfied by the other copy. Zero-probability choices are excluded.

If a target visit has already been revealed, E is certified and remaining
walk tails have their ordinary independent laws. Their *future* success
marks are independent Bernoulli `h(x_i)`, without the nonempty restriction.
Equivalently retain a hit flag a and use
`F_a = 1-(1-a) product_i(1-h(x_i))`. One must not keep applying the
pre-hit nonempty-mask formula after a hit.

For the original two independent all-active systems, let `I_-, I_+` be
their final target indicators and let `F_initial=P(E)`. Since either
nonzero indicator implies E,

`q_v(eta+delta_x)-q_v(eta) = F_initial E[I_+ - I_- | E]`.

If `F_initial=0`, both indicators vanish. If an initial particle is at v,
`F_initial=1`, so the conditioning is trivial. This identity supplies no
sign inequality; its use is to choose a probability space closed under
cross-copy path exchanges.

## 5. A common path/stack space obtained by draining all particles

The colour construction and the record constructions can be put on one
explicit probability space. Give every initial active particle a distinct
label, and fix a causal buffer scheduler which never uses an unspent
potential arrow to choose its next operation. After ordinary stabilization,
force every remaining resident to emit its waiting arrow, then ignore all
further sleep bits and route every particle to the original sink. This
**drain is an auxiliary continuation**, not ordinary ARW stabilization.

At any site, the emitted arrows are always an initial segment of its D-book.
After n arrivals the ordinary buffer has emitted `D_1,...,D_(n-B_n)`;
forcing its final resident emits the one remaining `D_n`. In the drain,
every later arrival emits its next arrow. The combined ordinary and drain
execution is therefore a legal routing of all initial particles to the
sink in the arrow-only stack network. Abelianity of that routing network
makes its total consumed prefixes independent of the B-books and scheduler.
Equivalently, they are the prefixes obtained by ignoring sleep from the
start. Sink accessibility gives almost-sure finite drains.

**Product-space lemma.** The complete labelled paths produced by this
execution are independent K-walks from their prescribed initial positions
to the sink. They are independent of the entire B-field. Conditional on
these paths and B, the unused D-tails remain independent original K-books.

**Proof.** Conditional on any fixed B-field, each emitted arrow is the next
unrevealed D-coordinate at its site. Choosing the moving label from the
revealed past does not change its K-row distribution. Thus the successive
steps of every labelled walk have the product walk law, independent of B.
The finite stopping-line property gives independent unused tails. More
explicitly, reconstruct the D-field by running the same ordinary and drain
scheduler on the prescribed labelled paths and B: assign each emitted path
step to its buffer record, then append its site's saved unused D-tail.
This recovers every coordinate. The finite-history density factors into
the product of all walk-step K-factors, the B-field factors, and the
independent residual-tail factors. No arrow or sleep bit is discarded.

Consequently Section 4's joint-success event E is also an event of the
raw D-fields: at least one of the two complete drains visits v. Either
nonzero final target indicator certifies E. Its probability is exactly
`1-product_i(1-h(x_i))`, over all initial labels in both copies.

**Lifting a causal path exchange.** Explore the ordinary labelled histories
until a recoverable causal contact between two labels at the same site.
Exchange their complete, untraversed walk tails, keeping their travelled
prefixes, B-fields and per-copy unused D-tails. A sleeping label's tail
starts with its still-unspent resident arrow. Reconstruct both D-fields
by the preceding algorithm. The two future walks start at the same site,
so their product K-law is unchanged. Their union of past and future visits
is unchanged, hence E is invariant. The actually used history before the
contact is unchanged, so the contact is recovered and the same exchange
gives the full inverse. Resident arrows may change before their records'
old arrival indices, but they were unspent and could not affect that prefix.

This proves a product-law-preserving lift, also under joint success
conditioning. B is fixed; arrows are exchanged within each site, and
unused arrow tails may shift by unequal lengths. One must not instead
shift the B-books along with these tails. The separate per-copy success
conditions still fail, as in Section 3.

A complete measure-preserving bad-to-good record injection also preserves
E on its domain automatically, since both endpoints have a target sleeper.
Thus such a map and the lifted contact moves have compatible probability
spaces and normalizations. This does not prove that an arbitrary adaptive
composition repairs its target bits, has disjoint images, or terminates.

Run `python3 -B verify_drained_paths.py`.
On 256 seeded four-site cylinders, 92 have the specified contact and 87
change their arrow fields; 31 contacts occur after both selected labels
have jumped. Every decomposition/reconstruction and contact inverse is
exact, including unused tails, all B-bits and site-specific arrow inventories.
Complementing every B-bit leaves the drain offsets unchanged. For target
3, absent initially, the 251 joint successes and five failures are preserved.
These audit the representation and local exchange, not occupation positivity.

**Fixing the common full paths is too restrictive.** With all p=1/2,
prescribe the common path `0->2->0->sink` and extra path `1->0->sink`.
Under the same FIFO convention, averaging the sleep bits gives target-0
occupation `5/8` below and `9/16` above. These are positive-probability
paths for `K(0,2)=K(0,sink)=1/2`, `K(1,0)=K(2,0)=1`.
The `frozen_paths` exact recursion in the same verifier checks these values.
This first calculation fixes the whole joint path bundle and keeps product
coin marginals under that conditioning.

The obstruction persists when the extra's complete walk is averaged out.
For the same kernel and FIFO convention, condition only on the common path
`0->2->0->sink`, an event of probability 1/4. The lower occupation remains
`5/8`, while the upper occupation is `269/448`, smaller by `11/448`.
The verifier solves the finite conditional Markov chain exactly, without
truncating the extra's walk. Without this conditioning the occupations are
`4/7` and `41/70`, so the actual increment is positive.
Thus sharing the entire common path identically cannot give an ordered
coupling **in this FIFO path representation**, even with one common
particle and arbitrary correlations between the two coin fields. This
does not exclude other particle-labelling/scheduling conventions or
contradict R's unconditional theorem. The new representation must allow
exchange of the common path geometry rather than fixing it in advance.

## 6. Cancelling components disconnected from the extra source

For two independent completed all-active histories from eta and
`eta+delta_x`, form an undirected graph on the finite vertices V. Include
the edge `{z,w}` for every **emitted** jump `z->w` in either history with
`w in V`. Exclude sink edges; two paths exiting at the sink do not interact.
Also exclude unspent resident arrows. Let Q be the connected component
of the target v, and let `E={x in Q}`.

**Proposition.** On the complement of E, exchanging the two **whole books**
at every site in Q is a measure-preserving involution. It exchanges both
final configurations inside Q and preserves them outside Q. In particular
it pairs every `(1,0)` target outcome there with a `(0,1)` outcome.

**Proof.** There is no emitted jump between Q and its complement in either
original history. Since x is outside Q, their initial particle inputs
restricted to Q are identical. After swapping books, execute each copy's
Q-input using the opposite original Q-history, and retain its old outside
history. These are legal completions: all old interior emissions stay
inside or exit to the sink, and all old exterior emissions stay outside or
exit. Thus their union graph is unchanged and Q is recovered. Applying
the operation again restores every record, including unused records.

For each fixed vertex set Q, whole-book exchange preserves the product
measure. The event that the selected component is exactly Q is invariant
under that exchange, by the preceding legal-execution argument. The finite
partition by Q therefore proves measure preservation for the adaptive map.
No conditional-freshness assumption or choice of favorable future records
is involved. Abelianity makes the emitted record sets and component
independent of the canonical execution order.

Writing `I^-, I^+` for the final target indicators, this proves the exact
connected-history reduction

`q_v(eta+delta_x)-q_v(eta) = E[(I^+-I^-) 1_E]`.

It imposes no sign on the remaining expectation. Unlike the static return
components in LD, Q is built from the two actual completed histories and
can be proper even on an irreducible graph. Unlike a proposed causal
collision test, the inverse uses an invariant completed-history component.
This does not contradict F7's latent-dependence obstruction.

The three-site forced-return cylinder in [AX §7](avalanche-exchange.md)
has `Q={0,1}` and extra x=2 outside Q; this involution repairs it. It can
also repair upper histories that never visit v. The seeded 1,024-cylinder
audit has 542 fields outside E, including 66 bad pairs, of which 27 have
no upper target visit. Every full inverse, component, interior support
exchange and exterior support preservation is checked. These counts are
implementation audits, not bounds on the general connected sector.

The AX maps can change the emitted union graph. Their images must not be
combined with this involution without an overlap argument.

## 7. What a visited-set exploration loses

On a directed acyclic jump graph, a stronger occupied-set coupling is
available. Process sites in topological order. At z let N_z be its total
number of initial and incoming active particles. If N_z is positive,
send `N_z-C_z` independent K-arrows onwards, retaining a sleeper exactly
when `C_z=1`, where `C_z` is a single Bernoulli(p_z) coin. If N_z=0,
send nothing and retain nothing. This is the correct local law: conditional
on N_z, which is determined upstream, the final arrival's sleep bit is
an independent Bernoulli(p_z), and the emitted arrows have the independent
K-row law. Intermediate reset bits do not change the total output word.

Share C_z and the outgoing arrow book between nested all-active inputs.
Both `C_z 1{N_z>0}` and `N_z-C_z 1{N_z>0}` are increasing in N_z.
Induction through the topological order therefore gives both nested
incoming counts and nested final occupied sets. Self-loops can first be
removed: with `t=K(z,z)<1`, the effective sleep probability is
`p_z/[p_z+(1-p_z)(1-t)]`, and an effective jump has law
`K(z,w)/(1-t)` for `w != z`, including the sink. Self-loop instructions
and failed surplus sleeps are null transitions. This covers graphs with
no directed cycle on two or more vertices, under the same all-active
assumption. It does not handle feedback between distinct sites.

An exact two-site example shows why the general problem cannot be reduced
to one completed visited set. Take `sink--0--1--sink`, nearest-neighbor
jumps with equal probabilities, `p_0=p_1=1/2`, and target 0. Suppose a
function phi on subsets represented occupation as

`q_0(A)=E_A[phi(visited set)]`

for every binary all-active input A. The visited set contains A. For a
singleton input it stays that singleton with probability 3/4 and is the
whole graph with probability 1/4. The empty and full inputs have fixed
visited sets. This triangular system has the unique solution

`phi(empty)=0, phi({0})=34/65, phi({1})=-2/195, phi({0,1})=22/39`.

The exact occupation probabilities used here are `8/15,2/15,22/39` for
inputs `{0},{1},{0,1}`. Thus phi is neither nonnegative nor increasing.
Visit monotonicity alone cannot prove the target inequality through this
representation, even though it does prove it on acyclic graphs. The check
`visited_set_obstruction` in `verify_response_renewal.py` verifies the
unique solution against the saved exact raw-ARW law. Section 6 retains two
complete histories and their emitted edges, so it is not ruled out by this
one-history visited-set obstruction.

## 8. Remaining composition problem and next test

Sections 1 and 4 originally used different coordinates. Section 5 now
provides a common augmented space, a full-book lift of causal path exchanges,
and joint-success normalization for complete bad-to-good record injections.
A proposed intermediate record move still needs to preserve E or account
explicitly for changing it. The remaining composition problem is to select
exchanges with a sign argument, recoverable cuts, disjoint images and a
termination or zero-residual proof.

The current [complete-drain inventory problem](drain-inventory.md) asks for
an injection after fixing the combined drain edges and all drain coin
totals. Its finite reduction and last-exit forest enumeration are proved;
an IDLA/closed-walk expansion proves the coefficient inequality when every
off-target coin total is zero or maximal. A stronger one-family pool also
permits positive formulas for one off-target sleeping site and a directed
three-cycle with a permanent extra. General cyclic mixed counts remain
open. Exchanging only opposite terminal
target bits cancels a subset within each pool when deletion makes no
target return. Whole-avalanche exchanges need not preserve these pools.
The exact five-word example in DI shows that a general pool cancellation
must also be allowed to change its last-exit forest. Exact integer
certificates now prove occupation monotonicity for general three-site
kernels with a permanent extra, with finite auxiliary-walk padding. On
directed three-cycles the padding has been removed: positive coefficient
arrays and their finite differences prove all inventory inequalities,
including arbitrary multiplicities and loop counts. Positive count-block
renewals are also proved for one vacancy or two particles on cycles of
every length. The next step is the four-cycle numerator cancellation;
all 28 nontrivial binary increment formulas are generated and audited,
but their signs remain open. Fixing common paths or fixed forests loses
the sign.

A contact traced backward from the final target sleeper is generally not
a forward stopping time. Section 2 now resolves first-contact address
selection for an exchange involving the creator's own stopped word, even
with repeated target visits. [The creator report](creator-exchange.md) composes
it with root-order reversal and a first-return cycle involution, retaining
the jump records which become unused. Continue with more general endpoint
returns in that construction. AX's complete-avalanche construction handles
arbitrary returns when the lower final sleepers outside the target occur
in the repaired upper support. Continue by accounting for missing lower
sleepers and upper histories that never visit the target, retaining their
creation laws and checking image overlaps before combining partial maps.
Continue toward the main conjecture under Jacob's 2026-10-10 instruction;
keep scalar matching and N1 parked absent a changed mechanism.

## Verification and sources

Run `python3 -B verify_stack_exploration.py`.
It checks 32 padded cylinders covering both saved examples and both
initial-extra conventions, their inverse maps, final target bits and exact
weights. Heterogeneous sleep probabilities are `(2/5,3/7,5/11)` where
applicable. The eight cylinders in each family sum respectively to
`3/175`, `9/1225`, `27/7700`, and `81/53900`, with equal image mass.
It also checks the `1/3` and `4/21` colour obstructions and joint-mask
invariance, including unequal h-values away from the exchanged labels.
Three additional creator excursions of lengths 1, 2 and 6 check that the
same permutation leaves all of their records unused in the image.
The common-wake witness checks returns to both switch sites. Six stopped-word
checks use heterogeneous probabilities, unequal numbers of target self-loops,
an additional return through site 3, two padding choices, and a later common
wake. They verify exact record weights, recovered cuts, full-book inversion,
and terminal supports. These finite checks support the general local proofs;
they do not test an unrestricted global injection.
The terminal-swap witness checks two distinct bad cylinders mapping to one
good cylinder, all three having exactly the same product weight.
Three creator-contact checks verify the new early-wake repair, its recovered
causal prefixes, full-book inverse, rational weights, and agreement between
the staged and FIFO executions. Six further pairs check the composed
creator/loop-erased and root-order repairs. Three full pairs check cycle
removal, recovery from unused books, and the fixed combined map. Each of
the two Fomin selectors and the cycle involution is checked on 6,400 short
path pairs, with full inverse and record-inventory checks. Small witnesses
distinguish the scopes of the different moves. See the creator report for hypotheses.
These are exact rational checks of the displayed constructions; the
general local assertions have the proofs above.
AX also checks two complete-avalanche witnesses, repairs both sources of
the terminal-overlap example, and audits 512 seeded finite cylinders:
59 matched bad pairs and 33 with unvisited unmatched sleepers repair with
full inverses; 28 bad pairs remain outside the extended domain.
These counts describe the audit's domain, not a monotonicity test.

Sources: [BS §§2–5](evidence/notes/buffer-switching.txt),
[R §1](evidence/notes/reconstruction.txt),
[ST §§5–6](evidence/notes/success-transform.txt). Imported evidence is
unchanged. This is an AI-written research note, without a novelty or
independent human-review claim.
