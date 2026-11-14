# Complete avalanche exchange with a matched sleeping background

Updated 2026-10-10 UTC. **Proved local injection, valid for arbitrarily many
common particles on its stated domain. The main conjecture remains open.**
Use the finite graph, accessible sink and `0<p_z<1` assumptions from
[the reader’s model and hypotheses](https://raw.githubusercontent.com/jfrichey/jfrichey.github.io/companions-2026-11-14/pagedocs/ai/arw-monotonicity.pdf). Records are `R_z(i)=(B_z(i),D_z(i))`
as in [SX](stack-exploration.md). A sleeping resident's potential arrow
remains attached to the record that created it until a later arrival wakes
it. Pending deliveries may be held.

## 1. A comparison for released resident arrows

Start with a sleeping support B, fixed waiting arrows `(d_b:b in B)`, fixed
future arrival books, and a list Q of pending packets. Compare it with a
system where selected residents in C subset B have already been removed
and their packets `d_b` emitted but held. The second system may also have
one additional sleeping resident w, with its own waiting arrow.

**Lemma.** Every vertex visited while processing Q in the original system
is visited in a legal completion of the second system. The initial waiting
arrows in common are the same; this is not a comparison with arbitrary
resampled resident arrows.

**Proof.** Follow the original deliveries in the second system. At the
first arrival to b in C, release its held packet `d_b` exactly when the
original resident would emit it. Both arrivals then use the same next
record and their buffers agree. At a visit to the extra resident w, hold
its emitted packet as an additional packet. If w also belongs to C, release
the held `d_w` to match the original resident, while holding the extra
resident's own arrow. At the end, release all remaining held packets and
finish. This is a legal schedule for the second system. Therefore its
complete visited set contains that of the original. The same construction
handles any number of extra initial sleepers at distinct available sites.

**Completed-word consequence.** In CX's upper comparison, let L end asleep
at target v and W end away from v. The two explorers start at v in the same
sleeping background B; other packets are held. If

`B intersect range(L) subset B intersect range(W)`,

then, after the exchanged words, the original upper system is obtained
from the image's outside state by additional resident releases and possibly
the extra terminal sleeper of W. Their remaining book words agree. Thus
original residual avoidance of v implies image residual avoidance of v.
In particular, when L is an immediate sleep at v, **W may wake arbitrary
other sleepers**. Avoiding those sleepers was sufficient but unnecessary.
This comparison does not permit new L-only contacts with B.

## 2. What deleting the final target jump leaves available

Take a completed upper history ending empty at v, with last target arrival
index m and record `(0,d)`. Replace that record by `(1,e)`, for any allowed
potential arrow e. The prefix before that arrival is unchanged. Afterwards
one packet at d has been suppressed. Freeze v and compare outside processes
with the same books: adding a pending packet can only increase their used
output prefixes and visited set. The original process never returns to v,
so neither does the new one. Hence its final support sigma contains v and
its last target index is still m. This is SX's terminal-record lemma.

There is a useful additional conclusion. After the new process has fully
stabilized, deliver the suppressed packet d into sigma using its unused
books and the actual resident waiting arrows. This avalanche **avoids v**.
Indeed, with v frozen, Abelianity allows one to stabilize all the other
packets first and d last. The combined outside execution reproduces the
original outside process, which never returned to v after arrival m.
The retained sleeper at v consequently remains untouched. A sink d gives
the empty avalanche.

If `B subset sigma\{v}`, the avalanche from active d among only the
sleepers B, with their same arrows and the same fresh books, also avoids v.
Delay the additional residents in sigma\B as in Section 1. This is a
pathwise visited-set assertion; no occupation attractiveness is used.

## 3. Exchanging full avalanche prefixes

Consider two prefixes which leave the same sleeping support B and the
same active input d. Their resident arrows may differ. Run each avalanche
to complete stabilization, in a fixed legal order. At z write the unused
book at the starting prefix as `U_z T_z`, where U consists of all records
read by this complete avalanche.

Exchange `U_z^-` and `U_z^+`, retaining each copy's T, and also exchange
the initial resident arrows at every b in B. Each copied avalanche now
has exactly the old opposite initial arrows and opposite read records.
It reproduces that avalanche and recovers every cut. Applying the same
operation again restores both books and the resident arrows.

Thus this is a measure-preserving involution on the specified data. Products
of record probabilities agree at each site; the exchanged initial arrows
have the same K-row law. Full stabilization is a finite stopping line.
The proof does not require avoiding contacts, loops, losses, or further
wakes: those are all part of the exchanged U words. Terminal potential
arrows inside U are retained as complete records.

## 4. Repair with both terminal indices selected adaptively

Take any finite all-active lower input eta and upper input `eta+delta_x`.
Consider a bad pair: the lower final support is `B union {v}`, while the
upper target is empty. Require that the upper history visited v; this is
automatic if eta initially has a particle at v. Let n,m be the respective
last target arrival indices. Their terminal records are `(1,e)` and `(0,d)`.
The lower creator and all common histories may have arbitrary returns.

1. Exchange those two complete target records. Let sigma be the fully
   stabilized upper image. Section 2 gives `v in sigma` and preserves m.
   This terminal transposition alone need not be injective when both
   indices are selected adaptively; the next steps are essential.
2. Restrict to the **matched-background event** `B subset sigma`.
   In the lower execution, follow its original history while holding the
   newly emitted packet d. There are no later original target arrivals,
   so this leaves exactly B and the held packet. Its read offsets are
   those of the original complete lower history.
3. In the upper unused books after sigma, explore an auxiliary avalanche
   from d among **only B**, using the resident arrows of sigma at B.
   Section 2 says this auxiliary avalanche avoids v.
4. Apply Section 3's complete-prefix exchange to that auxiliary avalanche
   and the lower avalanche from d among B, including their initial arrows.

The lower now executes the upper auxiliary avalanche and ends empty at v.
In particular, it makes **no target arrivals after n**. The actual upper
history and sigma are unchanged: step 4 alters only unused book suffixes
and arrows belonging to final sleepers, which its completed execution
never used. Thus the pair becomes `(0,1)` and its final target arrival
counts are precisely the original n,m. Abelianity makes these counts
recoverable in any fixed canonical execution order.

**Inverse and background recovery.** Run both image systems and read n,m
from their last target arrivals. Upper record m is still `(1,e)` and lower
record n is still `(0,d)`. Temporarily replace the latter by the former and
stabilize the lower input again. This **counterfactual sleep** reproduces
the original lower creation history, with final support `B union {v}`.
It recovers B, all its suffix offsets, and its resident-arrow addresses.
The arrows at B may have been exchanged, but they are unused terminal
arrows in this counterfactual history, so changing them cannot alter it.
Its used record prefixes otherwise survived step 4 unchanged.

Recreate the actual upper prefix ending at sigma. With n,m fixed, run the
two auxiliary avalanches from d on B using the current words and resident
arrows, and exchange those words and arrows back. Finally exchange target
records n,m back. Retain the decoded n during this inversion; recomputing
it after restoring the old lower avalanche could give a different index.
This recovers every original record and shifted unused tail.

Corresponding finite word cylinders have the same product weight, including
the resident arrows. The explicit inverse makes their images disjoint.
Thus the **complete composite** is a measure-preserving injection on its
stated domain, even though its first transposition alone can have overlapping
images. No external mark recording n, m, B, or the read lengths is needed.

The inclusion test is a domain restriction, not an assumed stochastic
domination or a renormalized sampling step. The result does not claim that
every bad pair satisfies it, nor disjointness from the images of other maps.

## 5. Unvisited unmatched sleepers may be retained

The domain of Section 4 extends without a different map-selection rule.
After its terminal transposition, let `C=B intersect sigma` and
`D=B\sigma`. Explore the upper auxiliary avalanche from d among only C.
It still avoids v by Section 2. Require additionally that it never visits D.
Exchange its complete words and C-resident arrows with a lower auxiliary
avalanche from d among only C, ignoring the lower sleepers in D during
this auxiliary exploration. This is Section 3's involution with background C.

The actual lower image starts with the full B, but its transplanted
avalanche avoids D. Those sleepers therefore remain untouched, and the
actual avalanche agrees with the copied C-avalanche. Hence it avoids v,
and both target indices are preserved. The inverse counterfactual sleep
recovers B and the actual upper history recovers sigma, so C,D and both
auxiliary initial states are recoverable. Apply the same complete-prefix
exchange and undo the terminal transposition. The original lower auxiliary
avalanche is allowed to visit D: it was an auxiliary C-process, not the
actual lower process. This distinction is essential for its stopping rule.

This is a single extension of the Section 4 injection, agreeing with it
when D is empty. Thus their images do not need a separate overlap argument.
It permits unmatched lower sleepers whenever the upper C-avalanche does
not reach them; it does not assume support containment in general.

## 6. Reserved sleep records provide missing resident arrows

A second extension uses the full background B in both auxiliary avalanches.
At each missing site b in `D=B\sigma`, find the first **unused sleep record**
in the upper book after the completed sigma history. Reserve that record,
skip the entire intervening upper prefix, and use its potential arrow as
the auxiliary sleeper's waiting arrow at b. Use sigma's actual resident
arrows on `B intersect sigma`. All auxiliary upper book words start after
the reserved prefixes. There is almost surely a first sleep, since `p_b>0`.
These reservations are data for a comparison process, not actual sleepers
added to the completed upper execution.

Explore this upper B-avalanche from d. **Restrict to its avoiding v.**
Exchange its full words with the lower B-avalanche's full words and exchange
the B-resident arrows. Both arrow coordinates at every b are sleep records,
so this also preserves the combined complete-record inventory. The upper
completed history is unchanged. The actual lower image follows the copied
B-avalanche and avoids v, preserving the target counts n,m.

The Section 4 counterfactual inverse recovers B and sigma. Every reserved
upper prefix has unchanged sleep bits, and its last sleep record has only
changed its unused potential arrow. Hence searching for the first unused
sleep recovers each reserved address. The two auxiliary avalanches and
their initial arrows can then be exchanged back, followed by the terminal
target records. This recovers the full books and proves injectivity and
equal cylinder weights on the stated domain.

When B is contained in sigma this is exactly Section 4's map. Its domain
also includes Section 5's avoidance domain, but their maps need not agree
when D is nonempty: Section 5 explores the lower auxiliary process on C,
whereas this construction uses B. Use this reserved-record rule uniformly
on its domain; do not combine its images with the Section 5 images without
checking overlaps. Trying successive reservations until one succeeds would
be a different, unproved selection rule.

## 7. Why reservations cannot enforce target avoidance in general

There is an exact obstruction to keeping B fixed and merely searching for
other reserved records. Use vertices `v=0, b=1, x=2`, with
`K(0,1)=K(0,s)=1/2`, `K(1,0)=1`, and `K(2,0)=K(2,s)=1/2`.
All vertices can reach the target and the sink. Set every sleep probability
to 1/2. The common roots are `{0,1}` and the upper extra root is 2.
In the lower history both roots sleep immediately, with waiting arrows
`0->s` and `1->0`. Choose upper records

`R_0(1)=(0,s), R_0(2)=(0,1),`
`R_1(1)=(0,0), R_1(2)=(1,0), R_2(1)=(0,s)`.

The original supports are `({0,1},{1})`. The terminal target indices are
`n=1,m=2`; exchanging those records leaves upper support `sigma={0}`.
Thus `B={1}` and the suppressed packet d is at 1. The upper avalanche on
sigma avoids v by sleeping at 1. But an auxiliary avalanche from d among
sleeping B must wake its resident at 1, whose only allowed arrow is to v.
It necessarily visits v, for **every** choice of reserved records and
future instructions. The check `forced_auxiliary_return` verifies this
positive-weight finite cylinder and its first forced emission.

This rules out making Section 6 universal by choosing more favorable
auxiliary arrows. It does not rule out changing B through creator-history
exchanges, permitting lower target returns with a different recoverable
cut, or a different initial coupling. The original conjecture is unaffected.

## 8. Exact checks and the remaining question

Run `python3 -B verify_stack_exploration.py`.
The `matched_avalanches` checks use a six-site complete jump kernel including
the sink, and `p=(2/5,3/7,5/11,7/13,11/17,13/19)`. They exchange full
avalanche blocks and one or two initial resident arrows. Their supports are

`({0,2},{2,4}) -> ({0,2},{0,2,4}) -> ({2},{0,2,4})`,

`({0,2,4},{4,5}) -> ({0,4},{0,2,4,5}) -> ({4},{0,2,4,5})`.

The first record exchange alone leaves a harmful lower wake. The complete
avalanche exchange removes it. The checks recover all offsets and records,
preserve exact weights, and verify target avoidance of both the selected-B
upper auxiliary avalanche and the full-sigma avalanche. These finite checks
audit the general local proof, not the unrestricted conjecture.

The `repaired_terminal_overlap` check applies the general composite to SX's
two distinct bad sources with cuts `(1,1)` and `(2,1)`. Their original
terminal transpositions coincide, but the new images have supports
`({1,2},{0,2})` and `({2},{0,2})`, and retain those two different cut pairs.
Both complete inverses and exact weights are verified. This resolves that
specific selection obstruction on the matched-background domain.

The seeded `general_terminal_repair_samples` audit checks 512 finite
four-site cylinders, using seed 20261010, complete jump rows including the
sink, and `p=(2/5,3/7,5/11,7/13)`. Of 120 bad pairs, 59 satisfy the matching
condition and all repair. Section 5 repairs another 33, leaving 28 outside
the extended domain. Across the 92 repairs, background sizes 0,1,2 occur
21,49,22 times. In 41 repairs n exceeds one; in 49 repairs m exceeds one.
Every image retains n,m, every inverse restores the full
field, and every site-specific record inventory agrees. Finite sink padding
makes each cylinder complete; no long history is discarded. These are
implementation checks, not statistical evidence for a global sign.

The `reserved_terminal_repair_samples` audit uses the same 512 sampled
prefixes, with an additional unused sleep and sink padding at each site.
It repairs 104 of the 120 bad pairs; 45 repaired pairs have missing
background sleepers, and 12 actually wake a reserved auxiliary sleeper.
All full inverses and site-specific complete-record inventories agree.
The 16 excluded upper auxiliary avalanches visit v. This audits Section 6's
single specified reservation rule, not a search over successful choices.

**Next:** account for missing lower sleepers actually contacted by the
upper auxiliary avalanche, upper histories that never visit v, and global
image accounting.
Adding missing sleepers invalidates the visited-set comparison; Section 6
therefore tests its auxiliary target avoidance explicitly. Section 7 shows
why changing only the reservations cannot cover every pair. Retain missing
sleepers' creator histories and the actual upper creation law. A decoder
allowing lower target returns is another possibility. The completed-word
comparison also broadens CX's
first/last-contact families, but their image accounting remains separate.

This is an AI-written research note without a novelty or independent
user-review claim. Imported evidence is unchanged.
