# Creator paths and loop-erased switching

Updated 2026-10-10 UTC. **Proved local history injections; the main occupation
conjecture remains open.** This report gives the current canonical proof of
the creator/explorer construction summarized in [stack-exploration.md](stack-exploration.md).
It uses the finite graph, accessible sink, and `0<p_z<1` assumptions of
[the reader’s model and hypotheses](https://raw.githubusercontent.com/jfrichey/jfrichey.github.io/companions-2026-11-14/pagedocs/ai/arw-monotonicity.pdf), with two common all-active particles.

## 1. The stopped words and their law

Use potential-arrow records `R_z(i)=(B_z(i),D_z(i))`: the incoming particle
stays if B is one; otherwise it follows D. An arrival also releases any old
resident's waiting arrow. Hold those released packets while following the
incoming explorer to its first sleep or sink exit. The explorer's own word
therefore has the ordinary one-particle killed-walk law, even when it wakes
residents. Every visited record includes its potential arrow, used or unused.

At two causal stopping prefixes, suppose the fresh explorers start at the
same vertex. At each site write their remaining books as `W^- T^-` and
`W^+ T^+`, where W is the word read up to first sleep/loss. Exchanging the
W words and retaining each T is a measure-preserving involution. The new
walks reproduce the opposite old walks, so they recover both cuts. Products
of record weights agree site by site, and the stopping prefixes leave fresh
independent tails. This includes unequal local lengths and repeated visits.
See SX §2 for the complete word-cylinder argument.

## 2. A creator exchanged at the first target contact

Take distinct sites v,a,y, common input `{v,a}`, extra upper particle y,
and target v. The two original book fields are independent. Use this legal
pending-delivery schedule:

1. Run the root from v in each copy to its first sleep/loss. Call the lower
   word L and the upper word U, and retain the branch where both sleep at v.
2. In the upper copy deliver y and retain its immediate-sleep branch. Keep
   its complete creator record and its factor `p_y`.
3. Start the upper root from a, holding every resident arrow it releases.
   Retain the branch reaching v before its first sleep/loss. Stop just
   before this first arrival at v is processed. Call the subsequent word
   from v through first sleep/loss W.

The upper prefix in step 3 is a causal stopping prefix and does not depend
on the lower word L. At the next arrival the explorer wakes the common
sleeper created by U. The extra may already have been woken along the
unchanged path from a to this contact.

Let **Phi** exchange L and W as in Section 1. It fixes all earlier upper
records. Consequently U, the extra's creator, and the explorer's first
contact are recovered by the same forward rule in the image. The opposite
new words recover the two variable-length cuts. Thus Phi is an exact
involution, conditionally on the retained upper prefix. Repeated visits to
v inside L or W create no inverse ambiguity.

For the repair statements below impose these additional conditions:

- L and W avoid y; L sleeps at v, while W sleeps at a site w different from
  v and y, or exits at the sink.
- The lower root from a, run after L, avoids v. Denote its ordinary stopped
  word by R. It is an ordinary word because it never wakes the original
  lower sleeper at v.
- The original upper residual avalanche, obtained by delivering the held
  packets after W stops, never visits v.

These events retain their original probabilities. They are not replaced by
renormalized walk laws. The original target bits are `(1,0)`.

**Upper-image lemma.** Under these conditions, Phi makes the upper target
occupied, regardless of the lower continuation R.

**Proof.** The new upper explorer follows L and sleeps at v. Both executions
release U's waiting arrow exactly once at the contact. The extra's arrow is
released in both or neither, according to the fixed pre-contact path, since
neither L nor W subsequently visits y. Hence the same held-packet list remains.

After the exchanged word, the upper outside suffixes are identical as words,
although their absolute indices may shift. The original outside state has
one additional sleeper w, with its original waiting arrow; the image omits
that sleeper. If W exited, the outside states already agree. Keep v frozen
and compare the two residual executions. The delayed-release construction
of [R §2](evidence/notes/reconstruction.txt) includes the image's visited set
in the original one: follow the image packets first in the original system,
holding the extra resident's arrow when w is first visited, then release it
afterwards. This is deterministic for the specified books and waiting arrow.
Since the original residual avoids v, the image never wakes its new sleeper
there. This proves the lemma.

If R also avoids w, Phi already gives a good pair `(0,1)`: the new lower
creator follows W, and the remaining root never wakes it or visits v. The
next section removes avoidance of w on the branch where R exits at the sink.

## 3. A second exchange handles a sink-ending lower root

Use the following classical two-walk switch. Loop-erase the first path and
choose its earliest vertex met by the second path. Cut the first path at
the incoming edge retained in the loop erasure, and the second at its first
visit to that vertex; exchange the suffixes. This is a weight-preserving
involution. The incoming-edge cut is essential: it is generally not the
last visit to the chosen vertex. We use Fomin's proof of Theorem 6.1 in
[Loop-erased walks and total positivity, §6](https://arxiv.org/pdf/math/0004083).

To apply it to stopped record words, give each sleep record `(1,d)` at z
its own absorbing endpoint labelled `(sleep,z,d)`, with edge weight `p_z K(z,d)`.
A jump record has edge weight `(1-p_z)K(z,d)`. A sleep's unused arrow is
therefore retained. Select intersections only at physical vertices. Restrict
the switch to pairs with exactly one physical-sink endpoint and one sleep
endpoint; this restriction is invariant under the switch. On other pairs,
or in the absence of an intersection, use the identity map.

**Book lift.** Read the two lower roots, first v and then a, holding old
resident emissions during the second word. Their stopped words are independent
killed walks. Apply the preceding switch to those words. At each site replace
the concatenation of the two old words by the concatenation of the two new
words, retaining the full old tail. The total word length at that site is
unchanged. The new first and second explorers recover the new words, so the
path involution gives a full-book involution, called **Chi**. Products of
record weights agree and the unused tails remain fresh. Resident arrows
which change from used to unused stay in their complete records.

**Proposition.** Under Section 2's conditions, if R exits at the sink, the
composition `Psi=Chi after Phi` sends every bad pair to a good pair with
exactly the same measure. No condition that R avoid w is needed.

**Proof.** After Phi, the two lower stopped words are W and R. If W also
exits, both roots disappear and Chi is the identity. Otherwise W sleeps
at w and R exits. If R misses the loop erasure of W, it misses w, so the
lower creator remains at w and the other root exits. If they intersect,
Chi exchanges their endpoints: the creator from v now exits, and the root
from a sleeps at w. The creator leaves no resident for the second root to
wake. In either case the complete lower endpoint is `{w}`, or empty in the
two-exit case, and v is empty.

The upper image remains occupied by the upper-image lemma, since Chi changes
only lower records. Thus the target bits are `(0,1)`. Each component is a
measure-preserving involution on its stated full-book domain, so their
composition is injective, with inverse **Chi first, then Phi**. In particular,
one must not run the two operations in their forward order when decoding.

This proves an injection on a complete, explicitly specified family of
all-active histories, including a second wake caused by the first exchange.
It does not sum together images of the different repair families.

## 4. Two sleeping endpoints and the remaining return obstruction

Extend Chi's invariant domain to include two sleep endpoints at distinct
physical sites, as well as one sleep and one sink endpoint. Keep it the
identity outside that domain. The same Fomin involution and book lift apply;
call this version Chi2. Endpoints are exchanged, so the domain is invariant.

**Proposition.** Keep Section 2's assumptions. Suppose W sleeps at w and R
sleeps at r different from w, W avoids r, and R first visits r at its last
arrival. Then `Chi2 after Phi` sends the bad pair to a good pair with the
same measure, even if R visits w and Phi alone creates a harmful wake.

**Proof.** If R misses the loop erasure of W, it misses w and no second
switch is needed. Otherwise let c be Fomin's selected intersection. Since
W avoids r, c is not r. The second-path cut precedes its first visit to r.
After exchange, the first lower root sleeps at r. The second lower path
is R's prefix up to c followed by W's suffix from c; both pieces avoid r.
It therefore reaches w and sleeps without waking the first root. The lower
final support is `{r,w}`, avoiding v. The upper-image lemma still applies.
The full inverse is Chi2 first, then Phi, with unchanged product weights.

The exact sufficient condition is that the switched second path avoid the
switched first endpoint. The stated first-visit and avoidance hypotheses
give a simple checkable version; no general endpoint-avoidance assertion
is being inferred from Fomin's theorem.

**Exact limit of this extension.** Consider two lower words with site
sequences `W: 0,5,3` and `R: 2,4,5,3,4`, ending asleep at 3 and 4.
Both terminal sleepers have waiting arrows to target 0. Fomin selects
c=5. The new first word ends asleep at 4, but the new second word
`2,4,5,3` visits 4 before sleeping at 3. It releases the waiting arrow to 0,
where the next record sleeps. The final lower supports before and after
Chi2 are `{0,4}` and `{0,3}`: the target stays occupied. This exact
full-book involution preserves weights and has an inverse, but is not a
target repair on this family. Erasing the loop through 4 would require a
separate recoverable move and a disjoint accounting of the discarded records.

## 5. Exchanging the fixed order of the two roots

There is a simpler bijection for a larger part of the two-sleeper case.
Read ordinary stopped words W from v and R from a, with resident emissions
held. At each site the book starts with `W_z R_z T_z`. Rewrite it as
`R_z W_z T_z`, and call this map S. In the image, explore a first and v
second: their words are exactly R and W. This is the inverse rule, which
rewrites them back in v,a order. Local total lengths do not change, and
products of record weights agree. Stopped-word cylinders therefore prove
that S is a measure-preserving bijection of full books.

**Proposition.** Under Section 2's conditions, `S after Phi` repairs the
bad pair whenever R exits, or R sleeps at r and W avoids r. No restriction
on earlier visits of R to r is needed.

**Proof.** Execute the image's lower roots in a,v order. An exiting R
leaves no resident. A sleeping R leaves just r, which W avoids. Thus W
finishes without waking a sleeper, and both lower endpoints avoid v.
The upper-image lemma applies as before. The inverse is S's inverse first,
then Phi. This concerns one fixed reversal of root order, not an adaptive
choice among orders based on previewed outcomes.

In particular S repairs Section 4's Fomin obstruction. Conversely Fomin
can repair some mutual endpoint visits which S cannot: take
`W: 0,4,0,5,3` and `R: 2,5,3,4`, asleep at 3 and 4. The visit to 4 is
erased from W before the selected intersection 5; the switched second path
avoids 4. Its final support is `{3,4}`. Reversing the original words' order
instead leaves `{0,3}` when the waiting arrow at 4 returns to target 0.
The more general sufficient Fomin condition is that its new second path
avoid its new first endpoint; the first-visit conditions in Section 4 are
only one way to ensure this.

## 6. Retaining a dependency cycle in the unused books

The following move handles a mutual endpoint dependence and explicitly
recovers the records that become unused. Here a word includes every
visited `(site,bit,arrow)` record, including its terminal sleep record.

Suppose the lower words W,R end asleep at distinct w,r. Require that R
first visits r at its terminal arrival and visits w earlier. Cut W at its
**last** visit to r, and require that its next visit to w is its terminal
arrival. Cut R at its first visit to w. Decompose

`W = A U e_w`, `R = B V e_r`,

where A runs from v up to but not including that last arrival at r, B
runs from a up to but not including its first arrival at w, `e_w,e_r` are
the terminal sleep records, and U,V are jump-only words from r to w and
from w to r. Thus `C=UV` is a closed jump word starting at r, with its
first return to r at the end and its first visit to w at the U,V boundary.
It may loop at other sites and revisit w during V. A may have visited r
or w earlier; B avoids both r and w.

**Cycle map.** Replace the two root words by

`W'=A e_r`, `R'=B e_w`,

and place C at the front of their unused books. Precisely, at every site
rewrite the old concatenation `W_z R_z T_z` as

`W'_z R'_z C_z T_z`.

All records retain their site and the total local length is unchanged.
Both terminal sleep records retain their unused potential arrows.

**Proposition.** This is a weight-preserving injection with an explicit
inverse. Under Section 2's upper-image assumptions its composition after
Phi repairs the bad pair.

**Proof.** The image's first root reads W' and sleeps at r; its second
reads R' and sleeps at w without visiting r. There are no residual packets,
so the lower support is `{r,w}` and avoids v. The upper is unchanged.

For the inverse read those two new stopped words. Their endpoints identify
r and w; their last records identify `e_r,e_w`. Starting at r in the unused
books, follow jump records until the first return to r, stopping before
reading another record there. This uniquely recovers C. Split it at its
first arrival at w to recover U,V. Reinsert the two middle words and put
the terminal sleep records back as in the displayed decomposition. This
recovers the entire original book, including its unchanged tail T.

Both sides are partitions into finite path cylinders, and corresponding
cylinders have exactly the same product of record probabilities. The new
unused cycle remains specified in the image cylinder; integrating it out
would invalidate this weight comparison. This proves the injection and
its law statement without selecting a favorable unused path for free.

**A full involution.** The removal domain D and its image I are disjoint:
in D the second word visits the first word's endpoint; in I it avoids that
endpoint. I is recognized by two distinct sleeping endpoints r,w, a second
word which avoids r and first visits w at its terminal arrival, and an
unused jump word from r that visits w and returns to r before sleeping or
exiting. Stop this auxiliary read at its first return. These conditions
give exactly the inverse above. Define Omega to remove the cycle on D,
insert it on I, and fix all other books. It is a measure-preserving
involution. Testing I terminates almost surely because all p are positive.

**One fixed composition.** `S after Omega after Phi` repairs both the
exiting-R family and the following two-sleeper family: W and R each first
visit their own distinct sleep endpoints at their terminal arrivals, and
R visits w. If W avoids r, Omega fixes the books and S repairs them. If W
visits r, Omega removes the cycle. Its new words avoid each other's sleep
endpoints, since A now avoids w and B avoids r. Reversing their order with
S is therefore harmless. The inverse is S-inverse, then Omega, then Phi.
This combines these branches under one invertible rule; it does not claim
disjointness from images of the earlier, different repair maps.

The simplest mutual cycle `W: 0,4,3`, `R: 2,3,4` defeats both Fomin alone
and root reversal alone when both waiting arrows point to 0. The cycle map
repairs it by retaining `4->3->4` in the unused books. More general repeated
visits to the two endpoints require a further rule: the first-return decoder
does not authorize deleting an arbitrary collection of cycles.

### A limit on repairs using only the two lower words

One cannot extend these maps to **every** pair W,R for which W starts at v
and ends away from v, while R starts at a and avoids v. On vertices 0,1,2,
take target v=0, a=2, jumps `0->1`, `2->1`, and `1->0` or sink equally,
with `p=(1/2,99/100,1/2)`. For independent ordinary stopped words W,R,

`P(W ends away from 0, R avoids 0) = (1-200/399)(1-1/400) = 199/400`.

But the entire probability that ARW from active `{0,2}` ends with target 0
empty is only `99301/238602`. The proposed source family exceeds the total
possible image mass by `3880699/47720400`. Thus no measure-preserving
injection acting only on the lower books can repair that whole family.
The exact raw absorption audit is `creator_capacity_obstruction` in
`verify_response_renewal.py`. This is an obstruction to dropping the upper
conditions, not to the main conjecture or the proved local families. In the
actual paired problem, upper residual avoidance also restricts W and its
creator's pending arrow; that retained joint law must participate in any
extension covering all remaining words.

## 7. A final-contact cut permits earlier upper returns

There is a partial replacement for Section 2's upper residual-avoidance
condition. Require the lower root from v to sleep immediately, so its
creator word L is the single, fixed record `R_v^-(1)`. The remaining lower
word R avoids v. Run the upper all-active system in a fixed canonical
pending-delivery order, for instance FIFO, and assume its final target is
empty. Let m be its last arrival index at v and retain the full prefix
before that arrival. This cut is selected from a completed history; it is
**not** being called a causal stopping time.

At that arrival follow the incoming explorer to its first sleep or exit,
holding all other packets. Call its word W. By the Abelian equality of
arrival counts, W cannot return to v: the prefix already used m-1 target
records. Let B be the sleeping support at the cut. W may wake arbitrary
residents of `B\{v}`. Earlier upper excursions and returns to v are allowed.
Exchange L and W below the respective prefixes, retaining both tails,
as in the stopped-word construction. Call this partial map Phi-last.

**Proposition.** The upper image has final target occupation one and still
has last target arrival index m. Phi-last is injective and preserves the
measure of its domain. It can be followed by any of the proved lower
repairs whose word hypotheses hold, including S when R exits.

**Proof.** At v both executions release the same old resident arrow, if
one is present. After its stop the original outside state differs from
the image by the forced releases of the residents visited by W, and its
extra sleeper w if W slept. The remaining outside book words agree.
The original residual cannot visit v, since its m arrivals were already
used. The released-resident comparison in [AX §1](avalanche-exchange.md)
therefore shows that the image residual also avoids v. Its new sleep at v
survives, and its
target arrival count is exactly m.

For the inverse, first undo the lower repair. Run the upper image to find
its last target index m and replay its canonical prefix to that arrival.
All records of this prefix were retained, so its offsets are recovered.
The word there is L; the lower first-root word is W. Exchange them back.
This recovers every record and every shifted unused tail. Partitioning by
the finite original and recovered histories gives matching word cylinders
with identical product weights. Thus the noncausal choice of m is supported
by an explicit inverse, not by an assumption of fresh conditionally selected
records. The lower index stays fixed at one, so SX's overlap from selecting
both terminal indices does not apply.

The immediate lower sleep makes every off-target old-resident contact
belong only to W. The released-resident comparison is why these contacts
are allowed. An arbitrary returning L needs an additional condition.

**A preceding sleep can mark a returning creator's cut.** There is a
separate injective variant when m=1, or the original upper target record
m-1 is a sleep. Allow L to be any stopped word from v ending asleep at v,
and require `B intersect range(L) subset B intersect range(W)`. Let n be L's number of
target arrivals. Exchange L and W at the original last-arrival cut as
above. Each word releases the old target resident exactly once; any old
resident contacted only by W is an additional forced release in the
original outside state. AX §1 gives the same residual comparison, including
W's possible extra terminal sleeper. Avoiding all other sleepers is a
sufficient special case of the contact inclusion.

The image's last target index is now `m+n-1`, so merely locating its last
arrival would give the wrong cut. However, L has only jump bits before its
terminal sleep. Thus the image's most recent target sleep **before** that
terminal one remains exactly record m-1, or the initial boundary zero when
m=1. Read this marker and add one to recover m, replay the canonical prefix,
and exchange the recovered words back. This gives the full inverse and the
same word-cylinder weight proof. Earlier visits to v inside L, including
excursions through other vacant sites, cause no ambiguity.

The marked variant and Phi-last are separate partial injections. Their
maps agree on the overlap where L is immediate and the marker is present,
but that alone does not make the images of their full domains disjoint.
This does not license merging them by an unchecked adaptive selection.

## 8. Exact checks and the next extension

Run `python3 -B verify_stack_exploration.py`.

The creator-contact checks allow repeated target visits and verify Phi's
causal cut, product weight, complete-book inverse, and agreement between
the staged execution and FIFO stabilization. The loop-erased selector is
checked on 6,400 short path pairs with each invariant domain: Chi makes
2,708 nontrivial switches; Chi2 makes 3,529. Both recover all input pairs
and preserve their full record multisets.

Three compound witnesses use common `{0,2}`, extra 1, target 0, sleep
probabilities `(2/5,3/7,5/11,7/13)`, and jump rows uniform on
`{0,1,2,3,sink}`, `{0,2}`, `{1,3,sink}`, and `{0,3,sink}`, respectively.
The lower creator, upper creator, upper explorer, and lower returning root
have self-loop counts `(0,0,0,0)`, `(2,1,3,2)`, or `(5,3,1,4)`. In each:

`original supports ({0},{3}) -> Phi ({0},{0}) -> Chi ({3},{0})`.

Phi alone creates a harmful lower wake: the transplanted sleeper at 3 sends
a packet to 0, where it sleeps. Chi makes that waiting arrow unused again.
All three fields have equal exact product weight; decoding in reverse order
recovers every padded record. Both execution schedules give the stated
terminal supports. The path checks audit the implementation; the propositions
above, together with the cited path-switching theorem, supply the local proofs.

Three further compound witnesses add site 4, replace the lower root's sink
exit by `3->4` followed by immediate sleep, and use `p_4=11/17`. Row 3 is
uniform on `{0,3,4,sink}` and row 4 jumps to the sink. Their supports are
`({0,4},{3}) -> ({0,4},{0}) -> ({3,4},{0})`; full inverses, exact weights,
and agreement of both execution schedules are checked. A six-site complete
jump kernel, with sink and heterogeneous p, checks the obstruction in
Section 4 without truncation.

The fixed root-order exchange repairs all six exiting/two-sleeper compound
pairs and the earlier Fomin obstruction, with full inverse and weight
checks. Two mutual-dependence examples distinguish the Fomin and root-order
maps. Three additional full pairs verify the cycle map and the fixed
composition S-Omega-Phi. They have supports
`({0,4},{1,3}) -> ({0,4},{0,1}) -> ({3,4},{0,1})`.
Their retained cycles allow 0,2,5 interior loop repetitions, including an
earlier visit to r in A. Both the first-return decoder and all untouched
records are checked. These use the six-site complete jump kernel with
sink and `p=(2/5,3/7,5/11,7/13,11/17,13/19)`.
Omega is also checked on 6,400 short word pairs, with 143 cycle removals:
all are reinserted exactly, and the site-specific record inventory agrees.

Three final-contact pairs have last upper target indices 3,5,8. Their
supports evolve as
`({0},{1,3,4}) -> Phi-last ({0},{0,1,3}) -> S ({4},{0,1,3})`.
The earlier first-contact exchange instead gives `({3},{1,4})`, leaving
the upper target empty because its old creator's arrow returns. The checks
recover the final cut and every book, compare exact weights, and use a
five-site complete jump kernel with sink and the first five heterogeneous
sleep probabilities listed above.

Three marked-cut pairs have `(old last index, L's earlier target visits)`
equal to `(3,2)`, `(5,1)`, `(8,4)`. Their new last indices are 5,6,12; the
decoder recovers the old cuts 3,5,8. One L makes the excursion `0->4->0`.
All full fields and exact weights are recovered, and S makes the resulting
lower target empty. This audits the marker formula with a changed index,
rather than presuming that last-arrival counts are preserved in this variant.

Two completed-word checks allow off-target wakes by W and, in one case,
by a returning L with the same old-resident contact. They check the full
inverse and exact weights under the contact-inclusion condition.

[AX](avalanche-exchange.md) gives a different complete-avalanche repair
with arbitrary lower creator returns, adaptively chosen terminal indices,
and any number of common particles when its background condition holds.
It also resolves SX's terminal-swap overlap on that domain. The preceding
marker is not required by AX's separate construction.

**Open:** general repeated/equal endpoints in the two-word construction,
L-only old-resident contacts in this upper comparison, unmatched backgrounds
in AX, and interactions among different repair families. Selecting
successive repairs needs a disjoint global ledger and a termination or
zero-residual proof. Equal weights of individual moves do not supply either.
The general occupation inequality remains unproved.

This is an AI-written research report, without a novelty or independent
human-review claim. The mathematical evidence retains its original status.
