# Higher-dimensional SFTs: technical arguments

Mathematical selection supporting the reader, with original working-proof and finite-certificate qualifications. These are ordinary proofs and exact finite structural certificates; independent human and priority review remain outstanding. The general equal-overlap conjugacy problem remains open.

## Setup and candidate hierarchy

Fix an alphabet A of size q >= 2, a finite support S in Z^d, and uniform product
letters. For a finite set V of occurrence anchors put
N_w = sum_{z in V} 1{x on S+z equals w}. The physical finite region is S+V;
extra unused sites multiply partition functions by a common factor. G(w) records
nonzero overlapping compatible displacements, excluding disjoint displacements.

The working results now settle both entropy extremes on EVERY two-dimensional
rectangle. Constant patterns are the unique maximizers. Whenever an empty-G
pattern exists, that is exactly the minimizing agreement class. The only
exceptions are binary 2x2, 2x3, 3x3 and their transposes: there the minimizers
have exactly one corner-agreement pair, G={+-(h-1,w-1)} or
G={+-(h-1,1-w)}. Every other agreement class has entropy strictly between the
minimum and maximum. The main remaining question concerns comparisons between
intermediate classes. The exceptional minimum proofs have exact certificates;
they do not assume the unproved agreement-inclusion comparison at q=2.

All 3,953,744 proper agreement inclusions on binary rectangles with both side
lengths at least two and area at most 25 now have all-anchor, nonnegatively
weighted count convex order, over every uniform q>=2. These are working proofs
with finite structural certificates, not extrapolations from small boxes.
Line sweeps give full field order for 3,953,656 of them; Reimer flips close
the remaining 88 count comparisons. At q=2 those 88 field comparisons remain
open; the NEW alphabet-cover certificates settle full field order at q>=3.
The flag-isolation theorem below now
gives STRICT pressure at every finite beta!=0 for ALL 3,953,744, over q>=2.
Strict entropy gaps between intermediate binary classes remain the central
unresolved endpoint problem.
The spare-symbol amplification below now gives strict entropy and strict
negative-coupling pressure for ALL these comparisons over every q>=3. The
remaining intermediate binary entropy issue is at q=2.

The scope now also extends beyond binary patterns on seven small rectangles:
2x2,2x3,2x4,3x3,2x5,2x6,3x4. Every agreement class on these supports has a
ternary representative, and 2,103 exact line-sweep certificates cover all
8,603 proper inclusions between the classes. Thus ALL patterns on these
supports have finite field order and strict finite nonzero pressure under G
inclusion, for every q>=2. The new shared-letter and axis-sweep arguments
below now give strict entropy for every q>=3 throughout area 25. A universal
ternary reduction is false: an explicit 3x7 agreement class requires four
symbols. The alphabet-spectrum section gives its short K4 obstruction.
NEW cover-merger lemma converts the complete alphabet spectra into a broader
theorem: EVERY pattern comparison on EVERY 2D rectangle of area at most 25
has field order and STRICT ENTROPY under proper G inclusion for q>=5.
The new shared-layer theorem below improves the strict-entropy threshold
to q>=4 throughout area 25, closing all 74 four-symbol critical covers.
NEW axis-sweep hard-endpoint amplification and shared-letter certificates
now close q=3 as well. Thus EVERY proper agreement inclusion on EVERY
rectangle through area 25 has STRICT entropy order for EVERY q>=3.
The remaining intermediate entropy issue in this range is binary q=2.
These conclusions use the analytic letter-merging proof, without assuming
the general finite-order conjecture.
The new critical-cover audit closes the finite-order side more broadly:
EVERY pattern on EVERY 2D rectangle of area at most 25 has all-anchor field
order under G inclusion at EVERY q>=3. There are 30,788 certified ternary
covers and 74 four-symbol covers not handled by the cover-merger lemma.
Together with the binary count certificates, this gives weighted count order,
weak entropy, and strict pressure at every finite beta!=0 for ALL patterns
and EVERY q>=2 throughout area 25. The cover reduction localizes small-alphabet hard-endpoint cases to upper
classes requiring all q symbols; the new arguments settle these at q>=3
through area 25.

The main conjecture remains strict entropy increase under strict agreement-set
inclusion for rectangular supports. The published cube version is CMRW,
Conjecture 6.3 (https://arxiv.org/html/2409.09024v1#S6.SS5).
The September working manuscript proves all positive-coupling pressure comparisons
and the strict comparison near zero on both sides; it leaves negative coupling
away from zero and the strict hard endpoint open.

A useful stronger candidate is:

**Finite convex-order conjecture (unproved).** For a fixed support S and every
finite anchor set V, G(u) subset G(v) implies N_u <=_cx N_v. Here convex order
means E f(N_u) <= E f(N_v) for every convex f, with equal means
|V| q^{-|S|}. Rectangular supports are the central target; arbitrary finite
supports are a deliberately stronger version still undergoing tests.

This implies finite-volume partition ordering for every real beta, by
f(n)=exp(beta n), and avoidance-count ordering by f(n)=(1-n)_+ on integer n.
It implies non-strict pressure and entropy ordering in the thermodynamic limit.
The flag-isolation theorem upgrades its all-anchor avoidance consequence to
strict pressure at every finite beta<0 under strict G inclusion. A strict
hard-endpoint entropy gap still needs a separate argument. For a finite maximum count m,
this conjecture is equivalent to nonnegativity of every coefficient of
[Z_v(t)-Z_u(t)]/(t-1)^2, where Z_w(t)=sum_x t^{N_w(x)}. Indeed the coefficient
at t^j is the unnormalised gap of E(N-j-1)_+; the zero total mass and mean
identities cancel the remainder. This is much stronger than a pointwise
partition inequality on t >= 0.

For the pressure/entropy programme, a weaker finite target is sufficient:
**all-anchor avoidance order**, namely Pr(N_u=0)<=Pr(N_v=0) for every finite
V under G inclusion. Retaining each anchor independently with probability
1-t shows that this is equivalent to finite partition order for every V
and every t in [0,1]; the reverse direction is the finite limit t down to
zero. Under strict inclusion the flag theorem then gives strict limiting
pressure for every finite beta<0. Thus full convex order is a stronger
research target, not a prerequisite for the negative-pressure programme.
Agreement inclusion already increases every all-ones mixed moment. The
open finite avoidance target concerns the all-zeros probabilities instead.
This separates the automatic positive-activity comparison from the missing
finite negative-activity input without assuming a total order on incomparable G.

The most useful formulation for further work is the following weighted version,
with its cylinder extension stated separately so that a failure of that
extension would not automatically disprove the translated-pattern assertion.

**Weighted compatibility conjecture (unproved).** For finite supports
S_1,...,S_n in a finite physical coordinate set, let A_i and B_i prescribe one
q-letter word on S_i, under the same uniform product measure. Suppose
A_i intersect A_j nonempty implies B_i intersect B_j nonempty for every pair
i,j. For every a_1,...,a_n>=0 and every convex phi:R->R,

E phi(sum_i a_i 1_{A_i}) <= E phi(sum_i a_i 1_{B_i}).

The means agree automatically, since each marginal is q^{-|S_i|}. Translated
patterns on a common rectangle are the primary specialization. Arbitrary
supports and arbitrary cylinder families are stronger conjectural extensions.
The interval theorem proves this for arbitrary intervals; letter mergers,
the chordal criterion, and binary Reimer flips provide other general sufficient
conditions. The area-25 certificates prove the rectangular specialization for
every binary pattern pair in that range, all finite anchor sets and q>=2.
The later alphabet-cover audit extends this to arbitrary patterns in that
range, and gives the stronger field comparison whenever q>=3.

**Strict rectangular endpoint conjecture (unproved).** If S is a finite
rectangle and G(u) is a proper subset of G(v), then h(X_u)<h(X_v), for the
same uniform q-symbol ambient full shift. This is a separate claim even if
the weighted conjecture holds. The primary programme is to establish these
two statements; full field order is a stronger auxiliary question, and the
scalar K-order and noncrossing candidates below are more speculative.

The local maximal-pair identity below now separates two roles. Convexity
provides an exact finite avoidance defect of one configuration. A separate
amplification mechanism must make its contribution proportional to volume.
The spare-symbol theorem supplies that mechanism over q>=r+1 once the finite
all-anchor comparison over r symbols is available. The binary packing target
below isolates what is still needed without a spare symbol.

**Field strengthening (unproved).** Under the same agreement inclusion, does
the vector of occurrence indicators increase in supermodular order? Here a
set function f is supermodular when f(A union B)+f(A intersect B)>=f(A)+f(B).
This would imply convex order for every nonnegatively weighted occurrence
count, not just the unweighted count. A deliberately broader question replaces
translated patterns by arbitrary prescribed cylinder events on fixed supports,
and replaces agreement inclusion by inclusion of their pair compatibility
graphs. Both versions have passed finite tests. The chordal criterion below
proves some comparisons; the new cycle-completion and deletion argument also
settles every family of at most six cylinders,
and a separate component argument settles any number of binary cylinders on
at most three physical coordinates. Exhaustive structural certificates also
settle at most eight binary cylinders on four physical coordinates. These do
not settle the general version.

Two additional, less established candidates now have a clear scope:

**Rectangular K-order candidate.** For the uniform q-symbol full shift and
patterns on the same finite rectangle, K_u(q)<K_v(q) implies h(X_u)<h(X_v).
This would strengthen the published inclusion conjecture. It has substantial
small-strip support but no proof here, and it is false for nonrectangular S
by the exact certificate below. K ties do not determine entropy.

**Rectangular negative-half-line noncrossing candidate.** If two such pressures
have distinct analytic germs at zero, their difference has the same strict sign
throughout beta<0 that it has for sufficiently small negative beta. In terms
of the first unequal derivative of order r, that sign is the sign of
(-1)^r [P_v^(r)(0)-P_u^(r)(0)]. This is a more speculative candidate. It alone
would imply only the weak endpoint entropy order; strictness at the endpoint
would need another argument. It is false for general finite supports and the
analogous positive-half-line claim is false even for 3x3 rectangles.

The core programme should not depend on either of these two stronger candidates
until further adversarial tests or a structural argument support them. Agreement
inclusion and finite convex order are the more stable primary targets.

**Reduction of the universal cylinder conjectures to fair bits (working proof).**
For each of the weighted-count and full-field conjectures, its version for
arbitrary cylinders on independent fair bits implies its version for arbitrary
finite uniform alphabets. The alphabet sizes may even vary by physical
coordinate. The reverse implication is immediate. This concerns the UNIVERSAL
cylinder conjectures; it does not reduce the single translated-pattern problem
to a single binary translated pattern, and it does not preserve event counts.

At physical coordinate x with alphabet size q_x, choose
M_x=floor(2^k/q_x), for k large enough that M_x>=1. Encode every letter a and
residue r in {1,...,M_x} by a distinct k-bit word e_x(a,r). Replace each event i
by a group of binary cylinders indexed by residue choices (r_x)_{x in S_i}.
Its old and new members request respectively e_x(u_i(x),r_x) and
e_x(v_i(x),r_x) on the SAME support S_i times {1,...,k}. Compatible old lifted
events have compatible original requests and equal residues on intersections;
their new lifted events are therefore compatible too.

At most one member of each group can occur. The probability that a compatible
set J of original groups all occur is
product_{x in union_{i in J} S_i} (M_x/2^k), tending to
product_x q_x^(-1). Incompatible sets have probability zero in both the lifted
and original models. Thus the grouped indicator laws converge to the desired
original laws. For weighted counts, give every member of group i weight a_i
and take the limit of the binary inequality.

For field order one needs the following extension, since taking the OR within
groups is not automatically a supermodular map. Given a supermodular f on
the original indices, let pi(F) be the set of groups hit by a lifted index set
F, and c_i(F) its number of members from group i. Choose
C>=max(0,max_{i notin A}[f(A union {i})-f(A)]) and define

g(F)=f(pi(F))+C sum_i binom(c_i(F),2).

For two added indices in different groups, the elementary mixed difference
is an elementary difference of f or zero. In the same group, it is C minus
a marginal increment of f if that group was empty, and C otherwise. Hence g
is supermodular. On every realized lifted field c_i is zero or one, so g=f
of the grouped field. Apply binary field order to g and then pass to the
finite-law limit. QED. In particular, any genuine finite-alphabet cylinder
counterexample would produce a binary counterexample after a sufficiently
large finite lift, although the number of events and bits may grow sharply.
`private_atom_checks.py` verifies 4,976 elementary inequalities for the
exclusive-group extension, including nonmonotone tests, and 284 equalities on
its realized domain. Literal lifts of the ten-event ternary family below use
10/34 binary events on 6/9 bits; all 1,152 configurations and every mixed
moment satisfy the stated encoding identity. These check the construction,
not the unproved universal binary comparison used as its hypothesis.

**Alphabet-extension identity (proved).** For patterns using only r specified
symbols, let Z_w^{(q)}(V;t) count q-letter configurations on Lambda, with
occurrences at V and weight t^{N_w}. For every integer q>=r,
Z_w^{(q)}(V;t)=sum_{E subset Lambda}(q-r)^{|E|} Z_w^{(r)}(V_E;t),
where V_E={z in V:(S+z) intersect E is empty}; the latter count is on Lambda\E.
Condition on exactly which sites carry symbols outside the specified r-letter
alphabet. This proves the identity, and also a polynomial identity in q.
Consequently finite convex order at alphabet size r, for all finite anchor sets,
lifts to every larger alphabet: condition on E and use the common conditional
mean |V_E|r^{-s}. This is a structural reason to seek the all-anchor-set
finite-volume conjecture. The large-q theorem below cannot be run in reverse.
For binary patterns the same conjecture is equivalent to nonnegative coefficients
of [Z_v(q,t)-Z_u(q,t)]/(t-1)^2 in the two variables q-2 and t, uniformly over
finite regions/anchor sets. The constant coefficient in q-2 gives the converse.
Exact polynomial checks for binary 2x3, 2x4, 3x3 on 3x3 anchors and binary 3x3,
3x4 on 4x4 anchors found no negative coefficient. These finite checks do not
supply the missing general combinatorial argument.

## Replacement conjugacies in every dimension

This section answers Jacob's November 8 question about extending the 1D swap
maps and organizing self-overlaps without a longest-border tree. The theorem
is an all-dimensional deduction from the mechanism of
[CMRW, Proposition 4.11](https://arxiv.org/html/2409.09024v1#S4.SS1),
whose stated setting is one-dimensional. No novelty or external audit is
claimed. The other results in this note compare counting laws; equal counting
laws alone do not supply a conjugacy.

Fix a finite nonempty support S in Z^d and a finite alphabet A. Write X_w for
the subshift of A^{Z^d} avoiding w in A^S. For a,b in A^S define the oriented
cross-agreement set

    C(a,b)={t in (S-S)\{0}: a(s)=b(s-t) for all s in S intersect (S+t)}.

Thus C(a,a)=G(a), and C(b,a)=-C(a,b). Neither the zero displacement nor
disjoint translates is included.

**Simultaneous-swap theorem (complete working proof).** Suppose u!=v and

    C(u,u)=C(u,v)=C(v,u)=C(v,v)=G.

Simultaneously exchanging every occurrence of u and v defines an involutive
automorphism F of the FULL shift A^{Z^d}. Its coding neighborhood is contained
in S-S, and its restriction is a conjugacy X_u -> X_v. The conclusion holds
for arbitrary finite S, with no rectangularity, convexity or cluster-size
assumption.

Proof. Put D={s in S:u(s)!=v(s)}. If t is in G and s in S intersect (S+t),
the four compatibility conditions imply

    u(s)=v(s)=u(s-t)=v(s-t).

Consequently a changed position belongs to no overlap between two distinct
occurrences of either pattern. All replacements are therefore consistent,
and every old occurrence becomes the opposite pattern at the same anchor.

There are no new occurrence anchors. Indeed, suppose y=F(x) has a u/v
occurrence at b that x did not have. Some site z in b+S must have changed.
That change came from an old occurrence at a!=b, which remains an occurrence
of the opposite type in y. The two y-occurrences coexist, so b-a belongs to
their cross-agreement set, which is G. Their intersection cannot contain
any changed position of the occurrence at a, contradicting the choice of z.
Thus the occurrence anchors are exactly preserved and their types exchanged.
Applying F twice restores x. The output at z depends only on old occurrences
covering z; their supports lie in z+(S-S). This gives a finite local rule
commuting with translations. Finally, x has no u occurrence exactly when
F(x) has no v occurrence. This proves the stated conjugacy. QED.

The argument covers infinite overlapping clusters and clusters with cycles:
the local rule never searches or orders a cluster. Pairwise agreement of
fixed symbol prescriptions already implies joint consistency; no additional
triple-consistency condition is needed here.

A nontrivial binary example on S={0,1} x {0,1,2} is

    u=000     v=010
      011       011.

All four cross-agreement sets equal {+-(1,-2)}. The local rule flips x[r,c]
exactly when

    x[r,c-1]=x[r,c+1]=x[r+1,c-1]=0,
    x[r+1,c]=x[r+1,c+1]=1.

The center bit is free, so these markers pick out both patterns. This is
not a global permutation of the alphabet. The exact checker exhausts all
2^15 local input contexts needed to compute F^2(x)[0,0] and the two output
occurrence predicates. All 32,768 contexts verify the involution and exact
exchange of occurrence types. This finite audit checks the displayed local
rule; the general theorem uses the proof above.

**Why equal G alone does not justify direct replacement.** Consider

    u=00     v=10
      01       00.

Both have G={+-(1,-1)}. On the constant-one background insert the patch

    100
    000
    000.

This configuration avoids u and has exactly one v occurrence, at (0,0).
Replacing that v by u produces the patch

    000
    010
    000,

which has a new v occurrence at (1,1). Hence this direct replacement does
not map X_u into X_v. This is a failure of the proposed map, not a proof
of nonconjugacy. The distinction is the same one made for 1D replacement in
CMRW Remark 5.7; its marked-occurrence counting theorem does not assert a
conjugacy.

**Finite graph certificates.** For fixed S,A and G, take all literal
patterns with that self-agreement set. Join two patterns when the theorem's
four-set condition holds, or when a global alphabet permutation takes one
to the other. Each connected component gives conjugate avoidance shifts.
Disconnected components do not imply nonconjugacy by other block codes.

| Shape | q | Realized G classes | Classes disconnected under these moves |
|---|---:|---:|---:|
| 2x2 | 2 | 6 | 2 |
| 2x3 | 2 | 18 | 10 |
| 2x4 | 2 | 44 | 28 |
| 3x3 | 2 | 63 | 45 |
| 2x2 | 3 | 7 | 0 |

In particular, the theorem and exhaustive 81-pattern graph certificate prove
that ANY two ternary 2x2 patterns with equal G define conjugate avoidance
shifts over the same ternary alphabet. For binary 2x2 there are no direct
swap edges between distinct patterns. Every site is a corner. If u and v
differ at a corner, the one-site overlap with the opposite corner would
need a compatibility truth value independent of the two pattern choices;
the two binary choices at the changed corner make this impossible. The
pair in the preceding example also cannot be joined by alphabet permutations
and these swaps on the same support. This does not exclude larger block codes.

**What replaces the longest-border tree?** In one dimension the proper
border words of w are precisely its longest proper border and that word's
proper borders. This makes recursion through one parent possible. In higher
dimensions S intersect (S+t) can be incomparable subshapes, such as a
horizontal strip and a vertical strip. A natural object to retain is the
family

    (t, S intersect (S+t), w restricted to S intersect (S+t)),  t in G(w),

organized by inclusion of overlap domains, keeping translation labels and
actual symbols. It can have several maximal overlaps. Its description as
a partially ordered family does not itself provide a recursive conjugacy
theorem or a renewal state. The longest-border renewal proof additionally
uses the linear order and independent unexplored future in 1D.

For a single swap, every coordinate in the union of the self-overlap domains
must be unchanged, and the cross-agreement condition prevents other windows
from becoming new occurrences. For a broader classification, the next task
is to find moves between the components of the finite graph above, perhaps
using larger markers and recursive overlap restrictions. Any proposed
cluster decoder must have a uniformly bounded coding neighborhood; an
unbounded replacement search does not define a sliding block conjugacy.
Equal G implies equal finite occurrence laws by the earlier arguments,
but a universal equal-G conjugacy claim remains an open target here.

Code and reproducible evidence: swap_conjugacy_checks.py and
swap-conjugacy-evidence.json. The checker uses the Python standard library,
exhausts the local contexts and the specified pattern graphs, verifies
overlap protection at every swap edge, and checks the finite patch in its
constant exterior. Run with --rebuild to regenerate the evidence, or without
arguments for exact replay. Both completed successfully on November 8.

## A proved comparison: merging letters

**Proposition (working proof).** Suppose v=f composed with u for a map f:A->A,
with the reference distribution remaining uniform on the SAME q-letter alphabet.
Then N_u <=_cx N_v for every finite support and finite anchor set. More generally,
the claim holds for any fixed family of prescribed cylinder events on arbitrary
finite supports, if their requested symbols at each physical coordinate are
merged by a coordinate-dependent map. Nonnegative weights on the events are
also allowed.

Proof. Process the physical coordinates one at a time. Conditional on all other
letters, let C be the contribution of events not involving the current site,
and let m_a be the contribution of otherwise satisfied events requesting symbol
a at this site. Before a merge the conditional convex expectation is
q^{-1} sum_a f0(C+m_a), where f0 is the convex test function. Merging two symbol
classes replaces two terms f0(C+r)+f0(C+s) by
f0(C+r+s)+f0(C), which is at least as large by convexity for r,s >= 0.
Disjoint mergers and a permutation of the remaining labels implement any map
of symbols. Iterate over coordinates. Each individual cylinder still has the
same probability q^{-|support|}, so means agree. Intermediate event families
need not be translates of one pattern; the two endpoints are. QED.

In fact the same proof orders the whole indicator vector in supermodular
order. Conditional on the other sites, the active event-index sets requesting
different symbols are disjoint. Merging two such sets replaces their two
indicator vectors by their union and intersection, precisely the defining
supermodular inequality. Thus every supermodular function of finitely many
occurrence indicators has larger expectation after merging. Weighted convex
count comparisons are a consequence. This stronger conclusion also transfers
through the equal-G and sufficient-alphabet reductions below; its extension
to arbitrary agreement inclusions is not proved.

**One-coordinate switching lemma (working proof).** More generally, take two
families of prescribed cylinder events on the same supports which differ only
in their requests at ONE physical coordinate. If every formerly compatible
pair remains compatible, the new indicator vector dominates the old one in
supermodular order. Conditional on all other coordinates, consider the active
events in one old requested-symbol class. Any two are globally compatible
before the switch, so their new requests cannot differ. Thus the switch
coarsens the partition of the active events into symbol classes, even if it
does not coarsen the full request partition. Apply the preceding conditional
supermodular inequality. QED.

This gives a broader possible proof route than uniform letter merging.
Replacing all coordinates by their target requests need not give a monotone
path: intermediate families may gain compatibilities which are later lost.
The endpoint inclusion hypothesis alone does not justify that replacement
argument. The switching lemma therefore remains a sufficient local move,
not a proof of the general conjecture.

**Directional switching theorem (working proof).** Fix a translation-invariant
total order < on Z^d. For an overlapping displacement d, let
M_w(d)={a in S intersect (S+d): w(a)!=w(a-d)} be its conflict sites. Suppose
G(u) subset G(v) and, for each overlapping d not in G(v),

min M_v(d) <= max M_u(d)

in this order. Then the entire finite occurrence field for u is below that for
v in supermodular order, for EVERY finite anchor set and every uniform alphabet
size q containing the patterns. In particular all nonnegatively weighted
counts increase in convex order, and h(X_u)<=h(X_v). Under strict G inclusion,
P_u(beta)<P_v(beta) for every finite beta!=0.

Proof. Sweep the physical coordinates in the chosen order, replacing all
incident requests from u by those from v at the current site. A pair with
displacement outside G(v) always has either an unprocessed old conflict or a
processed new conflict: the displayed inequality excludes a gap between the
last old conflict and first new conflict. Pairs in G(v) can only lose conflicts.
Thus the occurrence compatibility graph increases at every single-site step.
Apply the one-coordinate switching lemma. This works with any finite anchor
set; no box boundary assumption enters the finite-field comparison.

For strict limiting pressure, choose d in G(v)\G(u) and its last old conflict
y. Just before updating a translated y, prescribe the two witness copies on
their union except y, using v at earlier sites and u at later sites. These
requests agree everywhere except y. The two active symbol classes at y are
distinct before its switch and joined afterwards. The conditional partition
ratio therefore has the same lower improvement as in the strict merging proof:
delta=(t-1)^2/[q max(1,t^s)], where s=|S| and t=exp(beta).
The finite-energy lower bound for the witness remains
eta=q^{-(2s-1)} exp(-|beta|s(2s-1)). In free boxes, density-one sites admit the
witness, giving P_v-P_u>=log(1+eta delta)>0. This does not yield a strict hard
endpoint. QED.

The simpler sufficient condition M_u(d) intersect M_v(d) nonempty for every
d outside G(v) allows ANY physical-site processing order. The directional
condition is substantially broader. Equal-G replacement of either pattern is
free, and chains of certified comparisons preserve all the conclusions above.
Thus exact certificates for a finite collection of pattern classes can have
all-volume consequences; they are not finite-box numerical evidence.

**Ordered-interval field theorem (working proof).** Let I_i=[a_i,b_i] be finite
integer intervals indexed so that both a_i and b_i are nondecreasing; equal
intervals are allowed. For two families of prescribed cylinders on these fixed
supports under uniform q-letter product measure, compatibility inclusion implies
supermodular order of the entire indicator field, for ANY number of events and
every q>=2. Equal-length word windows are a special case. The next theorem removes the
endpoint-ordering restriction, allowing arbitrary nested intervals.

Proof. Write H_w for either compatibility graph, p_i=q^{-|I_i|}, and let T_w
be the strictly upper-triangular matrix with

(T_w)_{ij}=1{ij in H_w} q^{-|I_j\I_i|}, for i<j.

For selected indices i_1<...<i_k, compatibility of consecutive selected pairs
implies compatibility of the whole family: any physical coordinate belongs to
a consecutive block of selected intervals. Also the part of a new interval
already covered by earlier selected intervals equals its intersection with
the last selected one. Consequently its joint moment is
p_{i_1} product_{r=1}^{k-1}(T_w)_{i_r,i_{r+1}}.

Put lambda_i=t_i-1 and Lambda=diag(lambda_i). The multivariate probability
generating polynomial is therefore

Z_w(t)=1+p Lambda (I-T_w Lambda)^{-1} 1.

The inverse is a finite polynomial, since T_w Lambda is strictly triangular.
The resolvent identity gives the EXACT difference

Z_v-Z_u = sum_{i<j} [(T_v)_{ij}-(T_u)_{ij}]
                    (t_i-1)(t_j-1) F_i^v(t_{<i}) R_j^u(t_{>j}),

where F_i^v=E_v[1_{A_i} product_{k<i} t_k^{1_{A_k}}] and
R_j^u=E_u[product_{k>j} t_k^{1_{A_k}} | A_j]. These identifications follow by
expanding the left and right resolvents into selected increasing chains.
Both F and R have nonnegative coefficients, and their variables are disjoint.
Compatibility inclusion makes every matrix difference nonnegative. A monomial
t_B(t_i-1)(t_j-1), with i,j outside B, is exactly the signed measure for an
elementary supermodular difference at B. Hence the displayed identity expresses
the field-law difference as a nonnegative combination of these differences.
QED. Setting every t_i=t also gives a positive-coefficient formula for
(Z_v(t)-Z_u(t))/(t-1)^2, directly proving finite convex order of the count.

**General interval-field theorem (working proof).** Compatibility inclusion
orders the indicator field for ANY finite family of prescribed cylinders on
fixed integer intervals, allowing arbitrary nesting, at every uniform q>=2.
Order the intervals I_i=[a_i,b_i] by nondecreasing a_i, breaking ties arbitrarily.
For i<j define

D_ij={k<i:b_k<b_i} union {i<k<j:b_k<=b_i},
k(i,j)=j if b_j>b_i, and k(i,j)=i otherwise,
F_ij^v=E_v[1_{A_i} product_{k in D_ij} t_k^{1_{A_k}}],
R_ij^u=E_u[product_{k>j} t_k^{1_{A_k}} | A_{k(i,j)}].

Then the exact polynomial identity is

Z_v-Z_u=sum_{i<j} (1{ij in H_v}-1{ij in H_u}) q^{-|I_j\I_i|}
                    (t_i-1)(t_j-1) F_ij^v R_ij^u.

Proof. Expand Z in lambda_i=t_i-1 and process potential selected events by
increasing index. For a nonempty selected prefix, retain as its state the FIRST
selected interval with largest right endpoint. If the current state is i,
all previous requests visible to the next interval I_j lie inside I_i: their
starts precede a_j and their right endpoints are at most b_i. Thus adding j
requires only compatibility with i and contributes
lambda_j q^{-|I_j\I_i|}. Its new state is k(i,j). Skipping j has weight one;
selecting the first event j has weight lambda_j q^{-|I_j|}. These transitions
compute the full mixed-moment expansion exactly, even though lambda may be
negative.

Telescope the product of these finite transition matrices, using v-transitions
before the changed step j and u-transitions afterwards. The prefix weight of
state i is lambda_i F_ij^v: apart from i, exactly the indices in D_ij may have
been selected without changing its first-largest-endpoint status. The suffix
weight starting from state k(i,j) is R_ij^u, by the same conditional moment
calculation. Transition differences occur only at newly compatible pairs,
which gives the identity. The variables of F and R are disjoint and exclude
i,j; both are probability-generating polynomials with nonnegative coefficients.
Hence the difference is a nonnegative sum of elementary supermodular
measures. QED. No positivity of an interpolated compatibility kernel is used.

The checker independently verifies this more general identity on 42,518 exact
comparisons at q=2,3 over 22 interval systems, including strict nesting. The
initial endpoint-ordered resolvent proof remains a useful special case.

**Row/column switching corollary (working proof).** For rectangular patterns,
replace the directional criterion by the weak coordinate projection
min pi(M_v(d))<=max pi(M_u(d)), where pi is one of r,-r,c,-c. Then every
finite occurrence field increases in supermodular order. Process whole rows
or columns in that order. Compatibility is preserved at every line update.
Conditional on the other physical lines, the active events on the current
line are prescribed cylinders on equal-length intervals, so the theorem
applies. This proves non-strict finite partition, pressure and entropy order.
Strict limiting pressure requires an additional amplification argument, supplied
next; it does not follow just by counting one improvement per growing row.

**Separated-interval amplification (working proof).** Suppose all intervals
have length L and at most H intervals start at any given integer. Let
Q_w(b)=E_w[t^{sum_{a_i<=b} 1_{A_i}}], for 0<t<infinity, t!=1, and put
R(b)=Q_v(b)/Q_u(b). If the new compatible pairs (i_k,j_k), ordered from left
to right, satisfy a_{j_k}<=a_{i_(k+1)}-L, then

R(b)>=(1+c)^{#{k:a_{j_k}<=b}},
c=(t-1)^2 q^{-2L} min(1,t)^{2HL}/max(1,t)^{3HL}>0.

Proof. Apply the resolvent difference to the field restricted to starts <=b.
For any new pair i<j, their intervals overlap: disjoint pairs were already
compatible. Thus 0<=a_j-a_i<L. Writing m=min(1,t), M=max(1,t), the prefix
factor obeys F_i^v>=q^{-L}m^{HL}Q_v(a_i-L): the far-left events are independent
of A_i, and at most HL other preceding events remain. Similarly
R_j^u>=m^{HL}Q_u(right), where right comprises starts >=a_j+L and <=b.
These events are independent of A_j. The transition gain is at least q^{-L}.
Finally

Q_u(b)<=M^{3HL}Q_u(a_i-L)Q_u(right),

because these two far event families have disjoint physical supports and at
most 3HL intervening events were removed. Dividing the positive resolvent sum
by Q_u(b) gives the recurrence

R(b)>=1+c sum_{new pairs i<j:a_j<=b} R(a_i-L).

Keep only the chosen separated pairs. For the kth pair, all earlier chosen
pairs have ended by a_i-L. Induction therefore gives R(a_i-L)>=(1+c)^{k-1};
summing the geometric progression proves the bound. No interpolated transfer
matrix is assumed to define a probability law. QED.

**Strict line-sweep pressure theorem (working proof).** Under the row/column
switching criterion above, strict G inclusion gives P_v(beta)>P_u(beta) for
EVERY finite beta!=0. For a row sweep of an H by L rectangle, s=HL, one explicit
bound is

P_v-P_u >= [eta/(2L)] log(1+c)>0,
eta=q^{-2s} exp(-2|beta|s^2),
c=(exp(beta)-1)^2 q^{-2L} exp(-3|beta|s).

Proof. Choose a new displacement and the last old-conflict row of its two
copies. At that row update, all conflicts outside the row are already absent;
prescribing the common witness on the other rows activates both row cylinders.
The active cylinders have length L and multiplicity at most H at each start.
Their compatibility increases, so separated-interval amplification applies.
Place translated witnesses at horizontal spacing 2L; their interval pairs
satisfy its separation condition. Conditional on the other physical rows, let
K be the number of active witnesses. The conditional partition ratio is at
least (1+c)^K; the displayed c is a smaller uniform version of the lemma's.

Under the current old Gibbs weight, each prescribed outside-row witness fixes
at most 2s letters and changes at most 2s^2 occurrence indicators. Its marginal
probability is therefore at least eta, uniformly in the volume and sweep stage.
Averaging the conditional ratio and using Jensen gives a row partition ratio
at least (1+c)^{eta times number of witness slots}. In a free box, all but a
boundary fraction of rows admit volume/(2L)+o(volume) witness slots in total.
Multiply the row comparisons and pass to pressure. The same proof works for
columns and reversed directions. The constants vanish as beta tends to minus
infinity, so strict entropy does not follow. QED.

The general interval theorem also permits OBLIQUE line sweeps. For a primitive
integer normal (a,b), let pi(r,c)=ar+bc. If the same projected conflict criterion
holds, process the lattice lines pi=constant. A rectangular support intersects
each such line in an integer interval along the primitive direction (-b,a);
its lengths may vary. Conditioning off the line and applying the general
interval theorem proves field order. More generally this works for any support
whose sections along these lines are intervals, including lattice-convex S.

Strict amplification extends too. If interval lengths are at most L and at
most H events start at each line coordinate, then
F_ij^v>=q^{-L}min(1,t)^{2HL}Q_v(a_i-L), since D_ij includes every far-left
event and at most 2HL additional ones. The suffix bound still costs at most
HL events. The same prefix recurrence therefore holds with
c=(t-1)^2q^{-2L}exp(-3|beta|HL). For a pattern of s sites, H<=s suffices.
Choose witness slots every 2L primitive steps on each line; their lattice
density is 1/(2L). The preceding finite-energy/Jensen argument gives strict
pressure at all finite beta!=0 under strict G inclusion. No new hard-endpoint
strictness follows.

`interval_field_checks.py` checks the resolvent identity with exact rational
supermodular certificates on all covers of three binary interval systems,
including 29,319 comparisons at q=2,3,5. It separately checks unequal/repeated
intervals by literal configurations and tests the quantitative prefix recurrence.
The reproducible record is `interval-field-evidence.json`.

`switching_certificates.py` exhausts binary pattern classes and searches orders
(a*r+b*c,r,c). It saves a transitive reduction of the successful comparisons
in `switching-evidence.json`; default verification independently computes the
literal conflict-site sets for each saved pattern pair and checks the complete
class-poset coverage. Directional sweeps have the following coverage; row/column
sweeps now cover EVERY remaining pair in this table:

| Support | Proper G inclusions | Certified by directional chains |
|---|---:|---:|
| 2x2 | 7 | 7 |
| 2x3 | 53 | 53 |
| 2x4 | 264 | 264 |
| 3x3 | 465 | 465 |
| 2x5 | 891 | 889 |
| 3x4 | 3,357 | 3,357 |
| 2x6 | 2,438 | 2,438 |
| 3x5 | 16,411 | 16,400 |
| 4x4 | 30,672 | 30,668 |
| 2x7 | 6,112 | 6,107 |
| 2x8 | 13,856 | 13,855 |
| 3x6 | 60,725 | 60,721 |
| 4x5 | 189,989 | 189,863 |

Thus the field/convex-order conjecture is proved in working form for EVERY
binary-pattern inclusion on EVERY listed shape, over every uniform q>=2
background. Directional chains and the new separated-interval amplification
also prove STRICT pressure comparison at EVERY finite nonzero coupling for
all 325,240 pairs on the thirteen shapes, including those requiring line sweeps.
All these certificates give the previously missing WEAK entropy comparisons
between intermediate classes. They do not settle strict entropy gaps. The separate
ternary-merger certificates still give strict entropy for binary patterns at
q>=3 on all the shapes covered by those certificates.

**Axis-sweep cover proposal (REFUTED beyond the thirteen certified shapes).** Let G<H be a
cover in the inclusion poset of agreement sets realizable by binary patterns
on a fixed rectangle. There are binary representatives u,v with G(u)=G,
G(v)=H and an axis projection pi in {r,-r,c,-c} such that
min pi(M_v(d))<=max pi(M_u(d)) for every overlapping d outside H.
Equivalently, every binary agreement inclusion admits a chain of axis sweeps,
with equal-G replacements between steps. The ordered-interval theorem and
strict amplification would then prove the full finite-field and finite-coupling
pressure conjectures for binary rectangles, at every uniform q>=2, and weak
entropy inclusion. Strict entropy needs its own argument.

NEW exact enumeration finds four covers on 4x6 and fifteen on 5x5 that have
NO axis sweep between any binary representatives. Thus this proposal is false;
the obstruction is to the proof method, not the field or entropy inequalities.
The negative cases are independently audited by enumerating binary labels of
the forced components and checking every possible intermediate agreement set.
The general interval theorem supplies oblique sweeps for two of the 4x6
covers and four of the 5x5 covers. The remaining two/eleven admit NO linear
line projection between any representatives. Primitive normals perpendicular
to support-site differences exhaust the critical angular directions: the
admissibility conditions are closed and constant within the intervening angular
chambers, so any nonempty feasible chamber has a feasible boundary normal.
Thus the cover proposal remains false even with arbitrary oblique lines.

`line_sweep_checks.cpp` exhausts the larger binary agreement posets, while
`line_sweep_checks.py` independently checks each positive generator using literal
site conflicts and each negative cover by forced-component enumeration. Their
`line-sweep-evidence.json` settles every binary rectangle of area <=23 and also
2x12 and 3x8, up to transpose. The 4x6 and 5x5 certificates cover respectively
946,002/946,006 and 1,499,541/1,499,625 inclusions. Together with the earlier
shapes this is 3,953,656/3,953,744 proper inclusions on every 2D rectangle of
area <=25. The remaining 88 comparisons reduce to thirteen cover relations.
These are all-volume certificates, not finite-anchor numerical tests.

A stronger obstruction is explicit on 4x6. Let
G={+-(1,-5),+-(2,-2),+-(2,4),+-(3,-5),+-(3,1)} and
H=G union {+-(2,-5)}. Only one displacement pair is added, so no intermediate agreement set exists.
Modulo complement its old representatives are
u_0=(110101;101110;010111;111010),
u_1=(101000;000101;100010;010100),
and its new representatives are
v_0=(100000;000100;000010;010000),
v_1=(111101;101111;110111;111110).

Consider ANY assignment of totally ordered update times T(z), replacing all
requests incident to z directly from u_i to v_j while preserving every current
compatible pair. A non-target displacement with singleton old conflict a and
new conflict b forces T(z+b-a)<=T(z) for all z.
For the crossed pairs (i,j)=(0,1),(1,0), these displacements include respectively

(0,2), (2,0), (-1,-1), or their negatives.

Their positive relation (0,2)+(2,0)+2(-1,-1)=0 forces equality along each
step; they generate the even-parity sublattice. Thus each checkerboard class
must be updated simultaneously. For (i,j)=(0,0), singleton constraints instead
force equality along rows and T(r-1)<=T(r). The displacement (1,-4) has old
conflicts in rows {1,2} and new conflicts in row {3}, forcing T(r+3)<=T(r+2),
hence equality between consecutive rows. All sites then have one update time.
The pair (1,1) is the reversed argument. These implications are checked exactly
in `direct_block_obstruction()`.

Consequently NO uniformly bounded-block sweep, straight-line sweep, or sweep
by coordinatewise monotone curves can compare ANY representatives of these two classes
by direct endpoint replacement. It needs an infinite two-dimensional update
block, or a genuinely different intermediate request family. Since G<H is a
cover, chains through binary intermediate agreement classes cannot repair this
particular obstruction. This still does not refute field or entropy ordering.

Every cover on all thirteen listed shapes has such a certificate; oblique
site orders are unnecessary for this coverage. `--axes` reproduces the search,
and the default checker verifies all saved axis generators and their closure.
The stronger ONE-step statement for arbitrary inclusions is FALSE. On 4x5,
exactly ten inclusions lack a direct axis comparison between ANY representatives,
but all have chains. An explicit example, with both signs understood, is

G={+-(1,3),+-(2,2)},
J=G union {+-(0,4)},
H=J union {+-(3,-2)}.

A certified chain uses
u=(10010;00010;00100;11000),
w=(10001;00010;00100;00000),
v=(11001;10011;00110;00100).
The first sweep is in increasing r, the second in decreasing r. There is no
direct axis sweep from the class G to the class H: modulo complement, G has
only two representatives (bit codes 102665,456551), and H only two
(144179,209723), and all four representative pairs fail all four axis orders.
The checker exhaustively verifies these failures. This is a proof-route
counterexample, not a counterexample to field or entropy ordering.

Equal-G replacement is also essential: the specific 3x3 patterns
(101;010;100) and (010;100;000) have strict agreement inclusion but fail all
four axis criteria. Alternative representatives do give the comparison.

The 2x5 obstruction to this proof route is instructive. For
u=(01101;10010), v=(00100;00000), G(u) is a proper subset of G(v), but
two non-target displacements have singleton old/new conflict sites. Requiring
them to remain incompatible during a direct u-to-v physical-site sweep forces,
on every horizontal row, update-time inequalities

T(k)<=T(k+2), and T(k)<=T(k-1), for all integer k.

These force all T(k) to be equal: T(k)<=T(k+2)<=T(k+1)<=T(k).
Thus a direct sweep preserving compatibility would have to update an entire
horizontal row simultaneously. Ordering uniformly bounded physical blocks
cannot fix this obstruction on the infinite lattice. This does not exclude
intermediate request families other than direct u/v choices, or a row-wise
comparison theorem, and is NOT a counterexample to any ordering conjecture.

**Binary block-flip comparison (working application of Reimer).** Start with
binary prescribed cylinders on arbitrary fixed supports under fair independent
bits. Choose a set K of physical coordinates and a set B of event indices.
Complement precisely the requests in the incidence rectangle B times K.
If every formerly compatible pair remains compatible, every nonnegatively
weighted occurrence count increases in convex order. This statement concerns
weighted counts; full field supermodular order is not claimed here.

Proof. Condition on coordinates outside K. The active events are unchanged. Remove events with empty K-support; they
contribute a fixed count c. Write X for the weighted count of
active events outside B, and Y for the weighted count in B, as functions on
the remaining uniform bit cube. The new pair is (X(omega),Y(bar omega)), where
bar denotes bitwise complementation on K. Any old-compatible cross pair
(one index in B, one outside) must have disjoint K-supports: otherwise the
flip would destroy its compatibility. Thus whenever X>=x and Y>=y, selected
occurrence witnesses for these thresholds have disjoint K-supports. For
A={X>=x}, C={Y>=y}, we have A intersect C=A square C, with square denoting
disjoint occurrence. Reimer's inequality gives

P(X>=x,Y>=y)=P(A square C)<=P(A intersect bar C)
                         =P(X>=x,Y(bar omega)>=y).

The two coordinate marginals are unchanged. On their finite ordered value
sets, the expectation of a supermodular function of TWO real variables is its
marginal terms plus a nonnegative linear combination of these joint tails.
Apply this to f(x,y)=phi(c+x+y) for convex phi, where c is the fixed count
outside K; then average the conditioning. QED.

The external input is Reimer, *Proof of the van den Berg--Kesten Conjecture*,
CPC 9 (2000), 27--32, https://doi.org/10.1017/S0963548399004113, Theorem 1.2.
Its precise bit-complement form is also explicitly stated as Proposition 1.3
of van den Berg--Jonasson, https://arxiv.org/abs/1105.3862. Here bar C is the
image under bit complementation, NOT the event-theoretic complement. The
ordinary product bound P(A square C)<=P(A)P(C) alone would not prove the result.

The distinction from field order is substantive even when ALL nonnegative
weights are allowed. Here is an abstract four-indicator example (not cylinder
laws). Encode indicator sets by their bit masks. Let mu be uniform on
{3,5,9,10,12}, and nu uniform on {1,6,8,11,13}. The marginals on blocks
{0,1} and {2,3} are identical, and nu(U times V)>=mu(U times V) for every
pair of increasing subsets U,V of their two-bit Boolean lattices. Consequently
every nonnegatively weighted total is ordered in convex order, by the same
bivariate-tail argument used above. Nevertheless the set function with values
at masks 0,...,15

(0,0,0,1,0,1,0,2,0,1,1,2,1,2,2,4)

is supermodular and has E_nu f-E_mu f=-1/5. The checker verifies all 36
upper-set rectangle inequalities and every elementary supermodular inequality.
Thus neither weighted count order nor this stronger two-block upper-set order
can alone justify the field strengthening. This example does NOT refute the
cylinder or translated-pattern conjectures.

**Binary-character comparison (working corollary).** Let S be ANY finite
support in Z^d and chi:Z^d -> Z/2 a nonzero homomorphism. Suppose, up to a
global symbol complementation, v(s)=u(s)+chi(s) modulo two and G(u) subset
G(v). Then every nonnegatively weighted finite occurrence count for u is
below that for v in convex order, for EVERY uniform q>=2. Hence all finite
partition functions, limiting pressure, and avoidance entropy are ordered
non-strictly. Strict limiting gaps and full field order are not supplied by
the bivariate argument alone. The separate flag-isolation theorem below now
gives strict negative pressure under proper G inclusion, for any finite S
at the same alphabet q. The argument itself does not give full field order;
the later alphabet-cover certificates settle it at q>=3 for the area-25
patterns, leaving the binary q=2 cases separate.

For q=2, gauge the physical configuration by chi(x). An old occurrence at z
then requests v when chi(z)=0 and the complement of v when chi(z)=1. An
old compatible overlapping displacement d cannot have chi(d)=1, because
then v would conflict everywhere on its overlap, contradicting G inclusion.
Thus every old simultaneous pair from opposite anchor classes has disjoint
physical supports. Apply the block-flip theorem with K the whole finite
physical region and B the odd anchor class. The result is exactly the all-v
family. The alphabet-extension identity lifts every finite-anchor comparison
from q=2 to q>=2. QED.

This supplies comparisons beyond the obstruction above: u_0 xor v_0 is
exactly the checkerboard character. It resolves both remaining 4x6 cover
relations and eight of the eleven remaining 5x5 covers at the weighted-count
level. The following periodic construction closes the last three covers and
their 60 dependent comparisons. Combining all these generators with the line
certificates gives weighted finite convex order and WEAK entropy inclusion for
EVERY binary rectangle of area <=25. The full-field coverage remains
946,002/946,006 on 4x6 and 1,499,541/1,499,625 on 5x5; field order is not
inferred from Reimer. The separate flag-isolation theorem below DOES upgrade
all these count comparisons to strict negative-coupling pressure at q>=2.

**Two periodic flips for the residual 5x5 covers (exact structural certificates).**
Index sites (r,c) of the 5x5 square by 5r+c, with (0,0) first, and encode a
binary pattern w by sum w(r,c)2^(5r+c). The three representative pairs are

(u,v)=(10036534,9220), (10159414,132100), (10167606,140292).

For all three, u xor v is 3x3-periodic, equal to one exactly at residues
(0,1),(1,0). Write pi(r,c)=3(r mod 3)+(c mod 3). A periodic cylinder family is
specified by nine nine-bit row masks F_0,...,F_8: the event anchored at z
requests u(s) xor bit_{pi(z+s)} F_{pi(z)} at physical coordinate z+s.
Start with

F=(0,5,344,297,450,258,144,264,294).

For each of the three u, its compatibility with a translated event equals
G(u) at EVERY anchor phase. Now apply the following two incidence-rectangle
flips, with bit masks denoting sets of phases:

| Step | Physical-coordinate phases K | Event-anchor phases B |
|---|---|---|
| 1 | {3,5,6,8} (mask 360) | {2,3,4,8} (mask 284) |
| 2 | {0,1,3,4} (mask 27) | {1,2,3,6} (mask 78) |

Every current compatible pair remains compatible at each step. The final row
masks are (0,30,43,90,170,258,139,264,78), and their compatibility graph is
G(v) at every anchor phase. The checks enumerate all nine anchor phases and
all overlapping displacements using the literal site requests; disjoint pairs
are automatically compatible. Thus they apply to ANY finite anchor set.
On fixed supports, all joint moments of prescribed cylinders are determined by
pair compatibility and union volume. The source and target laws therefore
equal the homogeneous u and v field laws respectively. Apply the block-flip
theorem twice, then alphabet extension. This proves the asserted weighted
count comparisons. The intermediate families need not be translates of one
pattern. Full field order is not inferred; strict negative pressure follows
from the separate flag-isolation theorem below.

`reimer_checks.py` verifies these three paths, the ten character generators,
and their transitive coverage against `line-sweep-evidence.json`. Additional
exact checks comprise 328 bivariate tail comparisons on small pattern systems
at q=2,3, and 968 compatibility-preserving six-cylinder block flips, including
464 strict inclusions. The reproducible record is `reimer-evidence.json`.

**Reimer-only cover route is false for arbitrary cylinders (exact obstruction).**
On five fair bits, take six supports, in index order,

S_0={0,1,3}, S_1={0,3,4}, S_2={0,1,2},
S_3={0,1,4}, S_4={1,2,3}, S_5={2,3,4}.

Let H_u have edges {05,12}, and H_v have edges {05,12,23}. This is a cover
because it adds a single edge. Modulo independent physical-coordinate symbol
complements, the complete lists of request representatives are

u: (0,1,1,2,12,20), (0,1,3,2,12,20), (0,1,3,2,8,20);
v: (0,9,1,17,2,4).

Here entry i is the bit word requested by event i, restricted to S_i. The
checker exhausts all 2^(18-5)=8192 coordinate-normalized request matrices.
For every old/new representative pair, the local request cut changes at both
physical coordinates 0 and 3. Both coordinates must therefore be in any Reimer
flip K. At coordinate 0, the required event bipartition says B_0=B_1; at
coordinate 3 it says B_0!=B_1. This contradiction is invariant under coordinate
symbol complements. Thus no incidence-rectangle flip connects ANY pair of
representatives. Since this is a cover, a monotone chain of such flips and
equal-law replacements cannot connect the two graphs either.

The field inequality itself is particularly simple: the unnormalized law
difference is 2(delta_empty+delta_{2,3}-delta_{2}-delta_{3}). It is one positive
elementary supermodular difference. This shows why a general proof must allow
comparisons beyond Reimer flips; it does not refute any count or field-order
conjecture. The chordal criterion below covers this example.

**Private-atom transport lemma (working proof).** On a finite uniform sample
space Omega, let B_x be a fixed background indicator set at outcome x. Add
pairwise disjoint sets E_x of private event indices, disjoint from all
background indices, so the old field at x is B_x union E_x. Move E_x to an
outcome f(x), where B_x subset B_{f(x)}. Then the whole field increases in
supermodular order.

Proof. The functional graph of f consists of cycles with trees pointing into
them. On a cycle, all B_x are equal, so rotating its private sets to their
destinations merely permutes equally likely field values. Then process tree
vertices in increasing distance from their cycle. At a source x, the private
set E_x is still alone there. Its destination contains B_{f(x)} and possibly
other private sets. Moving E_x replaces these two indicator sets by their
intersection B_x and their union, since B_x subset B_{f(x)} and private sets
are disjoint. Each move is one supermodular inequality. QED. Empty E_x cause
no change. This lemma is for a fixed background; it does not justify changing
the background without a separate comparison.

**Private-atom extension lemma (working proof).** The fixed-background
restriction can be relaxed as follows. On the SAME finite uniform space Omega,
let B_x and C_x be two background indicator sets, with B below C in
supermodular order. Let E_x be pairwise disjoint sets of private indicator
indices, all disjoint from the background indices. Choose a destination d(x)
for each nonempty E_x, such that B_x subset C_{d(x)}. Then the fields

x -> B_x union E_x,
y -> C_y union (union_{x:d(x)=y} E_x)

are in the same supermodular order. Unlike the preceding transport lemma,
the background is allowed to change, but its own comparison is an assumption.

Proof. For a supermodular f, write M_f(E,B)=f(B union E)-f(B), with E and B
disjoint. This increment increases with B. At each destination y,
supermodularity and disjointness of the private sets give

f(C_y union union_{d(x)=y} E_x)-f(C_y)
  >= sum_{d(x)=y} M_f(E_x,C_y)
  >= sum_{d(x)=y} M_f(E_x,B_x).

Sum over y and add sum_y f(C_y)>=sum_x f(B_x), using the background order
for the restriction of f to background coordinates. The result is exactly
the claimed inequality. Uniform weights and the common sample space are
essential to this formulation. QED.

**Full-support deletion principle (working corollary).** For prescribed
cylinders on a fixed finite physical region, remove every event prescribing
ALL its physical coordinates. If compatibility inclusion holds for the whole
family, and the remaining background fields are in supermodular order, then
the whole fields are ordered too. Group the removed old events by their
unique physical configuration x. Identical old prescriptions have identical
new prescriptions, since compatibility is preserved; hence each group has
one destination d(x). Every old background event true at x is compatible
with that full-support event, so its new event is true at d(x). This is exactly
B_x subset C_{d(x)}. Apply the extension lemma.

The same principle holds for all nonnegative weighted-count convex comparisons:
use f(F)=phi(sum_{i in F} a_i) in the preceding proof; only the corresponding
weighted-count background comparison is needed. Thus a minimal counterexample
to either universal cylinder conjecture has no full-support event.

**Two-coordinate and laminar field theorem (working proof, every q>=2).**
Compatibility inclusion orders the whole field for any cylinder family on at
most two physical coordinates, with ANY common uniform alphabet size. Remove
the two-coordinate events. The background consists of one-coordinate events;
at each coordinate, old requested-symbol classes can only merge. Their field
order follows from coordinate-dependent letter merging. Restore all removed
events by the deletion principle. Empty-support indicators are constant.

More generally, the same conclusion holds whenever the supports are laminar:
any two are disjoint or one contains the other. Induct on the support hierarchy.
Distinct maximal supports are disjoint, so their fields are independent in
both models. Within each maximal support, delete its full-support events,
apply induction to the proper subsets, and restore them by the deletion
principle. Supermodular order is preserved under independent products, by
conditioning on all but one factor and replacing factors successively. QED.
The earlier binary proof below remains a special-case argument. The general
two-coordinate theorem also extends its switching corollary to every q>=2.
These two-coordinate and laminar conclusions also follow from the general
interval theorem after ordering the physical coordinates appropriately. The
new ingredient here is the deletion principle for an arbitrary background
comparison, including backgrounds with no such interval representation.
`private_atom_checks.py` verifies 729 ternary two-coordinate comparisons and
36 decorated nonchordal comparisons with 50 indicators, by exact sparse
elementary supermodular decompositions (15,169 terms in total).

**Two-coordinate binary field theorem (working proof).** For any number of
prescribed binary cylinder events depending on at most two physical
coordinates, compatibility inclusion implies supermodular order of the whole
field.

Proof. Regard events prescribing just one coordinate as the background and
events prescribing both as private atoms of the four equally likely outcomes.
Events with empty support are constant. Identical old prescriptions on the
same support remain identical in the target, by compatibility inclusion, so
the private sets E_{xy} each have a well-defined target outcome.

Permute target symbols separately at each coordinate. The single-coordinate
background then agrees with the old one, except possibly at a coordinate
where both old symbols are requested and their two classes are merged to 0.
Call such a coordinate collapsed. If neither is collapsed, compatibility
inclusion says precisely that each private atom moves to an outcome satisfying
all its old background events. Apply the transport lemma. If both are
collapsed, every private atom has target 00, and ordinary coordinate mergers
give the whole target field. If just the first is collapsed, first move each
E_{xy} to (x,y'), where y' is its target second-coordinate request. This
preserves its first-coordinate background and gains only background events
at the second coordinate, so the transport lemma applies. Now merge both
first-coordinate symbols to 0. Every private atom's target first coordinate
is forced to be 0 by its old compatibility with a first-coordinate event.
The resulting field is the target field. The other case is symmetric. QED.

**Two-coordinate switching corollary (working proof).** For binary cylinders
on arbitrary supports, changing requests at at most TWO physical coordinates
increases the whole field in supermodular order whenever no formerly compatible
pair becomes incompatible. Condition on all other coordinates. The events
active there form cylinder families on at most two fair bits; a locally
compatible active pair was globally compatible, so its compatibility persists.
Apply the preceding theorem to each conditional family, including its constant
indicators, then average. This extends the one-coordinate switching lemma.
It does not assert that arbitrary endpoint inclusions admit such a path.

**Corollary.** The constant pattern maximizes all finite-volume convex
occurrence statistics, all finite-volume partition functions at real coupling,
and avoidance entropy among patterns on the same support. This is a non-strict
entropy assertion. Equal agreement sets still imply equal laws by the earlier
mixed-moment identity. Alphabet permutations merely relabel the pattern.

**Proposition (working proof, strict finite coupling).** In the preceding
letter-merging setting, if G(u) is a proper subset of G(v), then
P_u(beta) < P_v(beta) for every finite beta != 0. This holds for any finite S.

Proof. Write s=|S|, t=exp(beta), and choose d in G(v)\G(u). In the two copies
on U=S union (S+d), let y0 be the lexicographically last site where the original
u-prescriptions conflict. All conflicts become agreements after applying f.
Process physical sites in lexicographic order, changing every incident request
from a to f(a). At an interior site y, use the translated witness pair with
last conflict y. At all its other sites the current requests agree: earlier
conflicts have been merged and no later conflicts exist. Fixing the at most
k=2s-1 letters in U+y-y0 except y therefore activates two distinct current
symbol classes at y whose images under f coincide.

For fixed outside letters, put m_a equal to the number of otherwise satisfied
incident requests for a; at most s requests are incident. Their old contribution
is sum_a t^{m_a}; after merging it is sum_b t^{sum_{a:f(a)=b}m_a}. Every merger
is nondecreasing. On the witness assignment at least one merger contributes
(t^r-1)(t^s1-1) >= (t-1)^2, with positive integers r,s1. The old sum is at most
q max(1,t^s). Consequently the conditional ratio is at least 1+delta there,
where delta=(t-1)^2/[q max(1,t^s)].

Under the current inhomogeneous finite-volume Gibbs weight, any prescribed k'
letters have conditional probability at least q^{-k'} exp(-|beta| s k'), since
changing them affects at most s k' indicators. This remains true after summing
out y. Thus the witness event has probability at least
eta=q^{-(2s-1)} exp(-|beta|s(2s-1)), uniformly over the stage and the volume.
Every interior compression multiplies the partition function by at least
1+eta delta; every other compression is nondecreasing. There are volume-o(volume)
interior witness sites. Hence P_v-P_u >= log(1+eta delta)>0. QED.

The bound is intentionally crude. In particular it degenerates at the hard
endpoint; it is not a proof of strict entropy. The assumption v=f composed with
u is essential to this proof and is stronger than agreement inclusion.

**Theorem (working proof, strict hard endpoint for letter merging).** Let S be a
finite lattice-convex set in Z^d, meaning S=conv(S) intersect Z^d, let
v=f composed with u, and suppose G(u) is a proper subset of G(v).
Then h(X_u)<h(X_v). Consequently the constant patterns are the unique entropy
maximizers, up to alphabet permutation, among all patterns on a fixed rectangle.

Proof. Pick d0 in G(v)\G(u), put U=S union (S+d0), and choose any conflicting
site y in S intersect (S+d0) for the two original u-copies. Two geometric facts
are used:
(i) every translate S+z contained in U contains y; and
(ii) every translate meeting U is contained in B=U+(S-S).
For (i), containment S+z subset S union (S+d0) implies
ell(z)<=max(0,ell(d0)) for every linear functional ell, by comparing support
functions and cancelling max_{s in S} ell(s). Hence z=t d0 for some 0<=t<=1.
For y in S intersect (S+d0), y-z=(1-t)y+t(y-d0) lies in conv(S) and is a
lattice point, so lies in S. Fact (ii) follows by subtracting two sites of S.
The segment restriction z in [0,d0] holds for EVERY finite S; lattice convexity
is used only in the following interpolation step.

Let A0 be the symbols used by u. A new agreement entails that f identifies two
members of A0, so |f(A0)|<q. Choose a symbol star outside f(A0). Pack translates
of B into large free boxes on an L-spaced lattice, with L larger than every
coordinate span of B. Process one block at a time: merge the requests at every
site of B except y, and then merge at y. Eventually process the remaining sites.
Every individual step increases or preserves the number of avoiding
configurations by the one-site convex argument, using (1-N)_+.

Consider the last step in one chosen block. Fix an assignment outside B which
admits at least one avoiding extension before this step. All requests at
B\{y} have already been mapped by f, while requests at y are original.
Assign the letters on U\{y} to the consistent v-prescriptions from the two
witness copies. Assign every site of B\U the symbol star. Any occurrence
meeting this collar is impossible, since its request there lies in f(A0).
Every translate meeting U but not contained in U meets the collar, by (ii).
Every translate contained in U contains y, by (i). Occurrences disjoint from B
are absent because the fixed outside assignment was extendible.
Thus every constraint not involving y is satisfied by this filling. Among
constraints involving y, the two witness copies request distinct original
symbols a,b with f(a)=f(b). If R is the set of active original requests at y,
the number of avoiding choices at y changes from q-|R| to q-|f(R)|, an increase
of at least one. For every OTHER filling of B\{y}, the change is nonnegative.
Hence, for this fixed exterior, the number of avoiding extensions increases by
at least one. There were at most q^{|B|} old extensions, so their number increases
by a factor at least 1+q^{-|B|}. Exteriors with no old extension contribute
nonnegatively. Summing over exteriors proves the same partition-count ratio.

There are asymptotically |Lambda|/L^d disjoint blocks. The intermediate count
need not be translation invariant, but its initial and final counts are those
of u and v. Therefore
h(X_v)-h(X_u) >= L^{-d} log(1+q^{-|B|}) > 0.
For a nonconstant u, mapping every used symbol to one symbol strictly enlarges
G: any two differently labeled sites give a displacement missing from G(u),
whereas all nonzero overlapping displacements belong to the constant pattern.
This proves uniqueness of the maximum on lattice-convex supports. QED.

No spare symbol in the original pattern is needed: the collar is processed
before the comparison site, which creates an unused requested symbol there.
This is why the hard-endpoint proof works despite failure of the finite-energy
bound as beta tends to minus infinity. For an arbitrary S, the same argument works if a newly compatible pair has
a conflicting site belonging to every translate contained in its union. That
condition is automatic for a primitive lattice displacement (the segment has
no intermediate lattice points), but not for every nonconvex support. For
example S={0,1,3,4}, d0=2 has S+1 contained in S union (S+2), while their
intersection {3} is missing from S+1.

## Rectangular backgrounds and certified entropy bounds

**Fixed-background insertion theorem (working proof).** Let
S=product_{i=1}^d {0,...,l_i-1}, where d>=2 and every l_i>=2, and let w be
one forbidden q-letter pattern, q>=2. There is a fixed background b, either
constant or of period two in one coordinate, into which EVERY locally
w-avoiding rectangular block can be inserted directly, without changing
any letters of the block. All symbols outside that rectangle are those of b.

Construct b so that no translate of ANY face of w matches the corresponding
restriction of b. A constant symbol c does this unless some face of w is
entirely c. Distinct-colored monochromatic faces cannot meet, so they must be
opposite faces in the same coordinate direction. Consequently at most two
colors occur as monochromatic-face colors. If any alphabet symbol is absent
from that list, use it for the constant background.

The remaining case has q=2 and opposite faces normal to some direction i
monochromatic in different symbols. Choose j!=i and let b(x)=x_j mod 2.
On the two i-faces, b varies along direction j, so cannot equal their constant
requests. Every other face contains both i-endpoints, whose requested symbols
are different, whereas b is independent of x_i. Those faces cannot match
either. This verifies the construction at both translation phases.

Now insert a locally avoiding block p on a rectangle R. Any translate of S
not contained in R has at least one entire face outside R, where the letters
are those of b; the preceding face property prevents an occurrence there.
Translates contained in R are excluded by local admissibility. QED.

More generally, any collection of such rectangular blocks can be inserted into
the SAME b provided no translate of S meets two of their rectangles. In
particular the one-pattern SFT is block gluing at gap max_i(l_i-1), and its
locally admissible rectangular language equals its globally admissible one.
Its periodic points are dense: repeat the inserted rectangle on a sufficiently
spaced grid whose periods preserve b, choosing even periods when necessary.
These assertions concern a chosen background, not filling relative to an
arbitrary fixed exterior. The frozen-corner obstruction is consistent with
them. Boyle--Pavlov--Schraudner, Appendix A, Example A.7, already gives the
one-corner 2x2 instance as a block-gluing shift that is not corner gluing:
https://cs.du.edu/~rpavlov/BPSsofic.pdf. Their Appendix C also shows why block
gluing in general must not be promoted to entropy minimality.

**Certified finite-box bounds (working corollary).** Let A_w(n_1,...,n_d)
count locally avoiding words on the n_1 by ... by n_d box. With natural logs,

log A_w(n) / product_i(n_i+l_i-1)
    <= h(X_w) <= log A_w(n) / product_i n_i.

For the lower bound, independently choose any avoiding word in each box of
a grid with periods n_i+l_i-1, using the fixed background elsewhere. A copy
of S cannot meet two boxes, so all these choices are admissible. Their
asymptotic density gives the lower bound. The upper bound is ordinary box
subadditivity. Thus exact finite counts give explicit convergent two-sided
entropy bounds; for cubic boxes the gap is at most
log(q)[1-product_i n/(n+l_i-1)], tending to zero at rate O(1/n).

This supplies controlled hard-constraint bounds and a common background for
separated blocks. It does NOT turn a fixed finite-volume comparison gap into
a strict limiting gap: the gluing cost lies on the box boundaries. The missing
strict-order mechanism must still give an advantage proportional to volume.
`rectangular_filling_checks.py` audits 101,886 binary/ternary patterns on
two- and three-dimensional boxes, 498,582 exhaustive binary block insertions,
and 3,000 seeded insertions of two separated blocks at q=2,3,4. Every tested
face and literal boundary check passes; the general conclusion uses the proof.

**Universal rectangular collars for q>=3 (working proof).** For a single
forbidden h by l rectangle w, h,l>=2, every inner rectangle R admits a
fixed collar of width g=4 max(h,l) with the following stronger property.
Write B=R+[-g,g]^2. Every translate of the forbidden support that meets B
but is not contained in R has a mismatch at a pinned site of B minus R.
Thus this one collar works with EVERY locally avoiding filling of R and
EVERY locally avoiding exterior of B. In particular X_w has the uniform
filling property (UFP). The construction uses no unused alphabet symbol;
w may use all q letters.

Choose c absent from the monochromatic-face colors as in the background
theorem; q>=3 guarantees it. Initially fill the collar with c. Treat each
outer corner separately, reflecting both coordinates so it is the top-left
corner. In these reflected coordinates, an occurrence crossing both adjacent
outer sides has its bottom-right subrectangle inside B. If the bottom-right
letter of w is not c, leave this corner unchanged.

Otherwise let v be the right column and r the bottom row of the reflected
pattern, both ending in c. Choose a!=c such that v is not a^k c^(h-k) for
any 1<=k<h. At most one non-c letter is excluded by these equations, so such
an a exists. Choose b different from both r[0] and r[l-1]=c. At this corner
insert the h by 3l rectangle whose every row is

E = a^l b^l c^l,

then reflect it back. The row r occurs nowhere in E or its continuation by
c's: it is not constant c (by the choice of c), a constant a or b run has
the wrong final letter, an a-to-b crossing ends in b!=c, and a b-to-c
crossing starts in b!=r[0]. A length-l word cannot cross both transitions.
Consequently this patch itself avoids w.

The four corner patches and R are separated so that no forbidden window
meets two of these rectangles; g=4 max(h,l) is more than sufficient. The
fixed-background insertion argument therefore excludes every window
contained in B except those wholly inside R. For a window crossing only
the top side, its bottom face is a length-l word in the top boundary strip.
It is either constant c or lies in one corner sequence E followed by c's,
possibly reflected. The row argument excludes it. For a window crossing
only the left side, its right face lies in the vertical boundary strip:
near a modified top corner the strip is a^h followed by c's. Its length-h
words are constant a, constant c, or a^k c^(h-k); the choice of a excludes
all of them. Reflection treats the lower corner and the other two sides.
Finally, a window crossing both top and left sides intersects B in a
rectangle of height <h and width <l. At an unmodified corner the requested
bottom-right letter differs from c; at a modified corner that intersection
is all a, whereas the requested bottom-right letter is c. Either way it is
blocked. No window can cross opposite sides of B. This proves the claimed
universal mismatch property. QED.

**Quantitative entropy minimality (working corollary).** If p is any locally
w-avoiding rectangular word on R, let b=|R+[-g,g]^2|. Then

h(X_w)-h(X_w forbidding p) >= -log(1-q^(-b))/b > 0.

Indeed, under the uniform measure on avoiding words in a large free box,
condition on the exterior of any B-cell contained in that box. The fixed
filling consisting of p and its collar is always admissible. There are at
most q^b conditional fillings, so the conditional probability of p at the
cell's designated inner position is at least q^(-b). Tile by disjoint B-cells
and reveal the designated p-absence events sequentially. Earlier events are
measurable outside the current cell; averaging the same conditional bound
gives probability at most (1-q^(-b))^M of absence in all M cells. Global
p-avoidance is a subset of this event. Taking logarithms per volume proves
the bound. It follows that every proper subshift of X_w has smaller entropy,
and every measure of maximal entropy has full support. No uniqueness of the
maximal-entropy measure is asserted.

The general implication UFP => entropy minimality is already recorded in
[Boyle--Pavlov--Schraudner, Appendix B](https://cs.du.edu/~rpavlov/BPSsofic.pdf),
which attributes it to Quas--Sahin, *Entropy gaps and locally maximal entropy
in Z^d subshifts* (2003). The argument above gives a direct explicit bound
for this construction. We do not identify UFP with strong irreducibility.
This is a result about one SFT and its proper subsystems; X_u and X_v are
usually not nested, so it does not by itself strictify their entropy order.
It also does not remove the binary frozen-exterior obstruction.

`uniform_filling_checks.py` checks every ternary 2x2, 2x3 and 3x2 pattern,
every four-letter 2x2 pattern (1,795 in total), and 1,000 seeded patterns
with q=3,...,6 and side lengths 2,...,6. All 3,416,263 translated-window
checks find a mismatch at a pinned collar site. These are universal checks
over the free interior and exterior, rather than sampled exterior tests.

## Strict negative pressure from all-anchor avoidance order

**Flag-isolation theorem (working proof).** Let S be any finite support of
size s, and suppose G(u) is properly contained in G(v). Assume the avoidance
probabilities are ordered for EVERY finite anchor set, over the same uniform
q-symbol background. Then P_u(beta)<P_v(beta) for every finite beta<0,
without changing the alphabet or requiring convexity of S.

Choose any d in G(v) minus G(u), put U=S union (S+d), n=|U|,
J=U-S={a-b:a in U,b in S}, and j=|J|. Let

p0=1-2q^(-s)>0,   rho=1+q^(-n)/p0.

For t=exp(beta) in (0,1), a finite anchor set V, and M pairwise disjoint
translates z+U whose endpoint anchors z,z+d belong to V, one has

log Z_v(V;t)-log Z_u(V;t)
  >= M p0(1-t)^2 t^(j-2) log(rho).

Here Z is unnormalized and both sides use the same physical region; unused
sites multiply both partition functions equally. For a rectangle with side
lengths l_i, translate these unions on a grid with cell volume
b0=product_i(l_i+|d_i|). This yields the explicit limiting bound

P_v(beta)-P_u(beta) >= [p0/b0](1-t)^2 t^(j-2) log(rho)>0.

For an arbitrary finite S, use any rectangular cell containing U. The flag
neighborhoods of different cells may overlap; no independence of good cells
is assumed.

Proof. Write F_w(T)=Pr(no w occurrence at the anchors T), and independently
retain each anchor in V with probability 1-t. For its retained set T,

Q_w(V;t):=E[t^(N_w)] = E_T F_w(T).

Thus all-anchor avoidance order already gives every finite negative-coupling
partition comparison. Tilt the flag law by F_u(T), with normalizer Q_u(V;t),
and call the resulting measure nu_u. In cell z+U put J_z=V intersect (z+J).
Call the cell good if T intersect J_z={z,z+d}. Every other retained support
then misses z+U. The avoidance factors on a good cell are p0 for u and
p0+q^(-n) for v. Different good unions are disjoint, so these factors
multiply. The all-anchor hypothesis applies to the residual retained family.
Consequently

F_v(T)/F_u(T) >= rho^(number of good cells),
Q_v/Q_u >= E_(nu_u) rho^(number of good cells).

The probability of goodness has an exact expression:

nu_u(cell z good)
 = p0(1-t)^2 t^(|J_z|-2) Q_u(V minus J_z;t)/Q_u(V;t)
 >= p0(1-t)^2 t^(j-2).

Indeed, after prescribing these flags, the endpoint events use only z+U,
whereas every remaining retained event has disjoint physical support.
Summing its other flags gives Q_u(V minus J_z;t). The final inequality uses
Q_u(V;t)<=Q_u(V minus J_z;t) and |J_z|<=j. Jensen's inequality now gives
the stated finite bound. A growing box contains a positive density of
disjoint witness unions, proving strict pressure. QED.

The conclusion deliberately stops short of the hard endpoint: the factor
t^(j-2) usually vanishes there. Unlike an extra alphabet symbol, these flags
can suppress surrounding constraints only when t>0. The theorem requires
all-anchor weak avoidance, not the full convex- or field-order conjecture.
Its local ingredients are checked by `thinning_strictness_checks.py`, including
two disjoint witness unions with OVERLAPPING flag neighborhoods and a holed
support whose unsuppressed local gain vanishes at t=0.
The exact rational audit covers 525 binary rectangular inclusions plus the
holed example, at q=2,3 and t=1/4,1/2,3/4: 60,864 conditional factorization
checks, 6,312 thinning identities and 6,312 good-flag probability identities.

**Consequence for the area-25 certificates.** All 3,953,744 binary rectangular
agreement inclusions through area 25 now have strict limiting pressure at
EVERY finite beta!=0, over every q>=2. The positive half-line was already
proved in the shared manuscript; the flag argument supplies the negative
half-line from the certified all-anchor avoidance comparison, including the
88 cases handled only by Reimer flips at q=2. Their full field-order question
remains open at q=2; the later alphabet-cover argument settles it at q>=3.
Strict intermediate entropy at q=2 is still open.

**Conditional avoidance ratios need not be ordered (exact obstruction).**
One cannot replace the all-anchor hypothesis by a claim that every one-anchor
conditional avoidance ratio increases. This fails even for the letter merger
u=(00;11), v=(00;00) on a 2x2 rectangle. Take anchors
i=(0,-1), j=(0,0), k=(1,1). The old compatible pairs are ij and ik;
jk is newly compatible. Their pair-union sizes are respectively 6,8,7,
and the triple union has size 9. Writing a=q^(-4), inclusion-exclusion gives

F_u(ijk)=1-3a+q^(-6)+q^(-8),  F_u(jk)=1-2a,
F_v(ijk)=F_u(ijk)+q^(-7)-q^(-9),  F_v(jk)=1-2a+q^(-7).

Hence the total avoidance probability strictly increases, but

F_v(ijk)/F_v(jk)-F_u(ijk)/F_u(jk)
 = -q^(-9)(1-q^(-2))(1-q^(-4)) /
   [(1-2q^(-4))(1-2q^(-4)+q^(-7))] < 0.

At q=2 the ratios are 429/452 and 213/224, with difference -45/25312.
Thus the conditional-ratio induction used at the overlap-free minimum does
not extend verbatim to general inclusion. The flag proof uses only the
unconditional all-anchor comparison and avoids this false strengthening.
The checker verifies the formulas at q=2,3,5 and performs an independent
literal binary enumeration.

## Strict entropy after adding a spare symbol

**Local maximal-new-pair identity (working proof).** Let S be ANY finite
support and G(u) properly contained in G(v). Choose a newly
compatible displacement d maximizing r=|S intersect (S+d)| among
G(v) minus G(u). Put U=S union (S+d), and include EVERY anchor z for which
S+z is contained in U. In the unique filling of U realizing v at 0 and d,
let C be the set of additional occurrence anchors and k=|C|. Over any uniform
alphabet containing the letters of u and v, the unnormalized count generating
functions satisfy

Z_v(U;t)-Z_u(U;t)=t^k(t-1)^2.

More precisely, the unnormalized field-law difference is the single elementary
supermodular vector delta_C+delta_(C union {0,d})-delta_(C union {0})-
delta_(C union {d}). These identities do not assume the conjectured global
ordering. If S is lattice convex, then k=0, so U has exactly one more
v-avoiding filling than u-avoiding filling.

To prove it, containment S+z subset U implies z lies on the segment [0,d]:
take convex hulls and subtract the support function of conv(S). The union
of any two contained copies, other than the endpoint pair, is a PROPER subset
of U: its convex hull has a strictly smaller extreme in the direction d or
-d. Such a pair has overlap strictly greater than r. Thus maximality makes
{0,d} the only newly compatible anchor pair inside U, even without convexity
of S. All newly compatible cliques consist of this pair plus a subset of C;
the unique endpoint-pair filling realizes all of C simultaneously. Their
unnormalized mixed moments all increase by one. Hence the generating-function
difference is (t-1)^2 sum_{J subset C}(t-1)^|J|=t^k(t-1)^2. Mobius inversion
also gives the asserted field-law identity.

If S is lattice convex there is no common compatible third anchor z.
Otherwise its pairs with 0 and
d would both be old-compatible. Every site in S intersect (S+d) lies in S+z:
if x and x-d belong to S, then x-z belongs to S by lattice convexity. The
third copy would identify the conflicting old requests of the endpoint pair,
a contradiction. Thus C is empty and k=0. QED.

One cannot drop convexity by choosing a different maximizer, or the farthest
maximizer. On S={0,1,2,3,5,6,7} in Z, take u=0000111 in support order and
v=0000000. Then G(u)={+-1}, while v has every overlapping displacement.
Among the new positive displacements, d=2 is the UNIQUE maximal-overlap
choice, with overlap four. Its union U={0,...,9} contains the three copies
at 0,1,2. The constant endpoint-pair filling forces the middle copy, so
Z_v(U;t)-Z_u(U;t)=t(t-1)^2. Both binary avoidance counts are exactly 1004.
This refutes the attempted pair-selection repair even when S has only one
hole. It is not an entropy-order counterexample.

It also pinpoints why the large-alphabet leading-term theorem below needs
convexity. In this example the q^(-10) contribution of the new d=2 pair is
cancelled by its forced triple. Enumerating the normalized compatible anchor
sets with union size <14 gives

h_q(X_v)-h_q(X_u)=2q^(-11)-q^(-12)-q^(-13)+O(q^(-14)).

The first nonzero entropy gain is later than the overlap count alone predicts.

In fact the cancellation exponent can be arbitrarily large, even with a
unique maximizer. For k>=1 put n=2k+2 and
S={0,...,n-1} union {n+k,...,2n+k-1}; label its two runs by 0 and 1 for u,
and label all of S by 0 for v. The old positive agreements are exactly
1,...,k. For k+1<=d<=n the overlap is 2n-k-d, while for d>=n it is at most
n. Hence the unique maximal new overlap is at d=k+1, of size n+1=2k+3.
Its union is the full interval {0,...,6k+4}, containing exactly the copies
at anchors 0,...,k+1. The new endpoint pair forces all k middle occurrences,
and the local partition difference is exactly t^k(t-1)^2. Thus the approach
to zero at the hard endpoint is not uniformly controlled by strict finite-
coupling positivity on arbitrary supports.

**Spare-symbol strictification (working proof).** Suppose u and v use r
symbols and have strict agreement inclusion on a lattice-convex S. Assume
their finite avoidance counts are ordered for EVERY finite anchor set in the
r-symbol background. Then h_q(X_u)<h_q(X_v) for every q>=r+1. Thus a spare
symbol upgrades the all-anchor weak comparison to a strict hard-endpoint gap.

More quantitatively, choose U as above and a rectangular box B containing
U+(S-S); write b=|B|. If the finite partition comparison holds at a fixed
t=exp(beta) in [0,1] for every anchor set, then

P_v^(q)(beta)-P_u^(q)(beta)
  >= (q^(-b)/b) log(1+r^(-|U|)(1-t)^2)
  >= q^(-b) (1-t)^2 / [b(r^(|U|)+1)].

At t=0 this is an entropy gap, with beta interpreted as the hard endpoint.
For a rectangular S with side lengths l_i, one may take B with side lengths
3l_i-2+|d_i|. No independence of the successful boxes is required.

For an arbitrary finite S, the same proof gives the analogous bound with
(1-t)^2 replaced by t^k(1-t)^2, using the general local identity. Hence the
spare-symbol argument still proves STRICT pressure comparison at every
fixed 0<t<1 without convexity. Convexity removes t^k and makes the bound
uniformly positive at the hard endpoint. More generally, any strict local
partition comparison on a finite physical region, with all its contained
anchors, can be amplified this way under the all-anchor weak-comparison
hypothesis; the pair identity is a concrete source of that local defect.

The argument also gives a FINITE packing bound. For a finite anchor set V,
let M_d(V) be the maximum number of pairwise disjoint translates z+B such
that z,z+d belong to V. Padding the physical region by unused sites is harmless.
With rho=1+r^(-|U|)(1-t)^2, one has

log Z_v^(q)(V;t)-log Z_u^(q)(V;t) >= q^(-b) log(rho) M_d(V).

Inside a good cell the local identity still holds if only a subset of the
contained anchors is in V, provided it includes the endpoint pair. The same
factorization and conditional-probability estimate therefore prove this
version too. Tiling a growing box recovers the displayed pressure bound.

Here is the amplification argument. Partition a large free box into disjoint
translates of B. Condition on the set E of sites carrying one of the q-r
extra letters. In the remaining sites, the active forbidden windows are
exactly those avoiding E, so the assumed r-symbol comparison still applies.
Call a cell good if its translate of U contains no extra letter and every
site of B minus U carries an extra letter. Any forbidden window meeting U
is contained in B; if not contained in U, it hits an extra letter. Hence each
good U factors out of BOTH conditional partition functions. The local identity
and Z_u^(r)(U;t)<=r^(|U|) give a ratio at least
rho=1+r^(-|U|)(1-t)^2 per good cell. The residual ratio is at least one by
the all-anchor hypothesis. Thus the total ratio is at least E_u[rho^M], where
M is the number of good cells under the q-symbol u Gibbs measure at t.

For every exterior configuration, there is a fixed filling of B that uses
only the original symbols on U, avoids u there, and uses one fixed extra
symbol on B minus U. It creates NO occurrence meeting B. Such an original-
symbol filling of U exists, for example by choosing a suitable constant word.
Its conditional weight is one after removing the common exterior weight;
every other filling has weight at most one for t<=1. Its conditional
probability is therefore at least q^(-b), including at t=0. Consequently
E_u M is at least q^(-b) times the number of cells. Jensen gives
E_u rho^M >= rho^(E_u M). Dividing logarithms by volume gives the bound,
since the cell density tends to 1/b. Finally log(1+x)>=x/(1+x) yields the
displayed rational lower bound. QED.

**New rectangular consequence.** All 3,953,744 proper binary agreement
inclusions on rectangles of area at most 25 now have STRICT entropy inclusion
over EVERY q>=3, and strict pressure at every beta<0 with the displayed
uniform hard-endpoint bound. This follows from the already-certified binary
all-anchor weighted convex comparisons; it requires no new enumeration of
those pairs. Strict intermediate entropy at q=2 remains open. More generally,
any future all-anchor binary proof immediately yields strict entropy at q=3.

Convexity of the support is material to the zero-background case k=0. On the holed
one-dimensional support S={0,1,2,4,5}, let u=00011 and v=00000. The maximal
new overlap at d=2 has size two, but U={0,...,7} contains anchors 0,1,2 and
Z_v(U;t)-Z_u(U;t)=t(t-1)^2. Both avoidance counts are 236 at q=2: the new
endpoint pair forces the middle occurrence. This refutes the k=0 conclusion
without its lattice-convex hypothesis, not the entropy-inclusion conjecture.

`spare_symbol_checks.py` verifies the local identity for 38,349 binary
agreement inclusions on eight rectangles, a triangle and a convex unequal-row
shape, at q=2,3,5. There are 262 additional literal binary field enumerations.
It also verifies both explicit holed-support obstructions, 1,596 full-field
identities on seeded irregular supports, and the arbitrary-cancellation
family for k=1,...,6, all at q=2,3,5. The seven-site obstruction additionally
has a literal field check and exact entropy coefficients through degree 13.
Saved evidence:
`spare-symbol-evidence.json`. These checks audit the local identity; the
limiting gap uses the separate uniform amplification proof above.

**Binary amplification target (unproved).** Suppose binary rectangular
patterns have strict G inclusion and their avoidance counts are ordered for
every finite anchor set. Fix a maximal-new-overlap displacement d and B as
above. Is there an epsilon(u,v)>0 such that, for EVERY finite V,

log Z_v^(2)(V;0)-log Z_u^(2)(V;0) >= epsilon(u,v) M_d(V)?

This is a concrete sufficient mechanism for the strict rectangular entropy
conjecture, stronger than a strict inequality at each fixed finite volume.
Its spare-symbol counterpart is proved above; its binary version is open.
The holed-support family shows that the maximal-pair argument needs its
geometric hypothesis, while the frozen-exterior example prevents a blanket
binary finite-energy argument. Neither obstruction refutes this rectangular
packing target.

**Strict finite avoidance candidate (unproved).** For rectangular S and
G(u) subset G(v), the finite avoidance counts at V are equal if and only if
V contains no pair with displacement in G(v) minus G(u). The reverse
implication is automatic from equality of every mixed moment on V; the
forward implication is a proposed strengthening of finite avoidance order.
It is false on holed supports by the examples above. The maintained
`finite_avoidance_checks.py` checks every anchor subset of a 4x4 grid for
binary 2x2, 2x3, 2x4 and 3x3 patterns, and every subset of a 3x3 grid for
3x4 patterns. All 48,578,033 inclusion/anchor-subset cases containing a new
pair have STRICT binary avoidance gain; cases seeing no new pair have equal
counts. Integer subset transforms are independently checked against clique
enumeration on selected anchor sets. These remain finite tests and do not
prove either the candidate or the uniform packing target.

Under the all-anchor weak comparison over r symbols, the finite equality
criterion DOES follow over q>=r+1 on lattice-convex S. For a given V seeing
a new pair, choose a newly compatible pair of maximal overlap AMONG pairs
in V. In its union U, the same local proof applies to the contained anchors
that belong to V: every other contained pair has larger overlap and thus
is old-compatible if new-compatible, while convexity excludes a common third
neighbor. The local hard gain is one, and a single spare-symbol cell makes
the full finite gain positive. This proves that part of the candidate for
the certified area-25 binary comparisons over every q>=3, leaving their
q=2 finite equality criterion as a further unproved strengthening.

## Strict entropy with a shared layer

**Product-alphabet theorem (working proof).** Let S be an h by w rectangle,
h,w>=2, and use the uniform product alphabet of size q=rk, r,k>=2. Write
u=(a,m), v=(b,m), where a,b use r symbols and the marker m uses k symbols.
Assume G(a) is contained in G(b), every finite-anchor avoidance count for a
is at most that for b in the uniform r-symbol background, and
G(u)=G(a) intersect G(m) is properly contained in G(v)=G(b) intersect G(m).
Then h_q(u)<h_q(v). The proof gives a positive pressure gap uniformly for
0<=t=exp(beta)<1, including the hard endpoint.

Choose d in G(v) minus G(u) with maximal full-S overlap. Put
U=S union (S+d), n=|U|, and R0=bounding rectangle of U+(S-S). Its side lengths
are 3h-2+|d_1| and 3w-2+|d_2|. Let g=4 max(h,w), B=R0+[-g,g]^2 and b0=|B|.
For any finite anchor set V, let M be the number of pairwise disjoint translates
z+B chosen so that z,z+d belong to V. Extra physical sites can be added without
changing partition-function ratios. Then

    log Z_v(V;t)-log Z_u(V;t)
        >= M q^(-b0) log(1+r^(-n)(1-t)^2),       0<=t<=1.

Tiling growing free boxes by B gives

    P_v(log t)-P_u(log t)
        >= q^(-b0)/b0 * log(1+r^(-n)(1-t)^2),

with entropy at t=0. All constants are independent of the finite volume and
of its exterior. No independence between favorable cells is assumed.

**Marker isolation.** For any rectangular marker m over at least two symbols
and any d in G(m), there is a marker configuration whose occurrences at 0,d
exist and whose every m-occurrence meeting U is contained in U. Reflect
coordinates to put d=(a0,b0)>=0. Fill U with its compatible pair of m copies.
In its bounding rectangle R, fill the upper-right gap with a color different
from m's top-right corner and the lower-left gap with a color different from
its bottom-left corner. Every S-window contained in R but not U has the
corresponding corner in one of these gaps. If a0=0 or b0=0, there are no gaps.
Outside R use the face-avoiding constant or period-two background constructed
above. Its face argument excludes every window not contained in R, regardless
of the interior filling. Thus all remaining occurrences are contained in U.
Restriction to R0 fixes every marker window meeting U, independently of the
exterior. Extra occurrences contained in U are allowed.

**Local identity.** Let C be the active marker anchors contained in U, possibly
restricted to V, but including 0,d. Among pairs in C, the endpoints are the
only pair newly compatible for the first-layer patterns. Indeed, any other
pair of S-windows has union properly contained in U and hence larger overlap.
Its simultaneous marker occurrences force its displacement into G(m), so a
new first-layer agreement would contradict the choice of d. No third active
anchor can be compatible with both endpoints under b. Such an anchor lies on
[0,d]; its pairs with the endpoints have larger overlap and therefore were
already compatible under a. Rectangular convexity puts the endpoint
intersection inside the third window, forcing the old conflicting endpoint
requests to agree. This is a contradiction.

Consequently the only changed mixed moment is that of the endpoint pair.
The pair specifies every first-layer site of U, so its unnormalized moment
increases by exactly one. The local count polynomials therefore satisfy

    Z_b(C;t)-Z_a(C;t)=(t-1)^2.

This remains true for any subset C retaining the endpoints. Since
0<Z_a(C;t)<=r^n for 0<=t<=1, its ratio is at least
rho=1+r^(-n)(1-t)^2. The positive lower bound follows by a constant filling
that avoids a, choosing a different symbol when a is constant.

**Conditioning and amplification.** In a finite physical box Omega containing
S+V and the selected B-cells, condition on the whole marker field eta. Active
anchors are its m-matches, and the first letters are independent and uniform
on r symbols. Call a selected cell good when its marker on R0 matches the
isolating patch. For every good cell, all active first-layer windows meeting
its U are contained in that U. Thus its first-layer partition function factors
out of both systems. The remaining anchors avoid all good U's; their partition
functions are ordered by the all-anchor hypothesis. That hypothesis at t=0
implies it at every t in [0,1] by independent anchor thinning. If J(eta) counts
good cells and f_a(eta),f_b(eta) are the conditional normalized partition
functions, then pointwise

    f_b(eta)>=rho^J(eta) f_a(eta).

Integrating and dividing by Z_u gives
Z_v/Z_u>=E_{mu_u}[rho^J], where mu_u is the OLD composite Gibbs measure,
including its uniform avoidance measure at t=0. This is an expectation under
the tilted marker marginal, not under the unweighted product marker law.

For each selected B-cell, fix its isolating marker on R0 and put a constant
first symbol there that avoids a. This inner joint filling avoids u.
The universal q>=3 collar for u extends it across B and blocks EVERY
u-occurrence meeting B, for every exterior. Conditional on an exterior of
positive mu_u probability, cancel the common weight from occurrences disjoint
from B. The prescribed filling has remaining weight one; every filling has
weight at most one when t<=1, and there are q^b0 possible fillings. Hence
Pr_mu_u(cell good | exterior B)>=q^(-b0), including t=0. This also holds for
arbitrary finite anchor subsets, since the collar blocks more occurrences
than required. Thus E_mu_u J>=Mq^(-b0), and Jensen proves the finite bound.
The B-cells can tile an expanding rectangle with density 1/b0; their U's are
well inside their respective cells, and only a boundary proportion lacks the
endpoint anchors. This proves the limiting bound. QED.

**Exact finite evidence and consequence.** `layer_strictness_checks.py`
checked marker isolation on all binary patterns of 2x2,2x3,2x4,3x3 and all
ternary patterns of 2x2,2x3: 2,765 compatible pairs and 98,107 translated
windows. It checked 28,032 local first-layer identities. All 74 critical q=4
covers through area 25 have a shared binary marker and ordered binary
first-layer classes: 8 on 3x7, 10 on 3x8, 38 on 4x6, 18 on 5x5. Their
literal composite agreement sets, local identities and universal-collar
mismatch certificates were checked. The analytic argument above, rather
than the finite checks alone, supplies the limiting gap.

The factorizations encode the four target symbols as two bits. Any ordered
two of the three balanced cuts (low bit, high bit, xor) encode the target
injectively. Select a binary source class A contained in the new first-layer
class and satisfying A intersect G(marker)=G(source). All 74 searches
succeeded, including the two genuinely four-symbol-to-four-symbol covers
on 5x5. The existing binary area-25 theorem supplies the first-layer
all-anchor comparisons. All other q=4 covers are cover mergers. Therefore
EVERY pattern on EVERY 2D rectangle of area at most 25 has STRICT entropy
order under proper G inclusion at EVERY q>=4. The next axis-sweep argument
also settles q=3; the binary q=2 endpoint remains open.

**Shared-letter line theorem (working proof).** Let q>=3 and let u,v be
patterns on the same 2D rectangle with strict G inclusion. Suppose there is
a symbol c whose position set is the same in u and v. If a primitive linear
projection pi satisfies

    min pi(M_v(d)) <= max pi(M_u(d))

for every overlapping d outside G(v), then h_q(u)<h_q(v). More generally the
same conclusion holds whenever their finite conditional avoidance counts
are ordered for every marker configuration and every anchor subset under
the common marker m=1{pattern=c}. The displayed line condition is a sufficient
structural certificate for that conditional comparison. The symbol c can be
absent from both patterns; proper G inclusion excludes T={m=0} being empty.

To prove the conditional line comparison, fix any physical marker field eta.
At its zero sites the letters are independent uniform on the r=q-1 other
symbols; at its one sites the letter is c. Retain only anchors whose full
marker window equals m. Process physical lines in pi-order, replacing all
requests on each line simultaneously. The projected-conflict condition
preserves compatibility at every step, just as in the original line proof.
Condition additionally on all zero-site letters off the current line. The
active events on the line prescribe letters on

    (S+z) intersect line intersect {eta=0}.

Each full-S line section is an interval. Deleting marker-one sites from the
COMMON physical line and compressing its order makes each nonempty support
an interval again, although T itself may have holes. Empty supports give
identical deterministic indicators and can be set aside. The general interval
field theorem therefore orders each conditional update. Integrating and
composing updates proves the entire field order for every fixed eta and
anchor set. No comparison for arbitrary translates of the holed support T
is assumed.

Now choose a maximal NEW FULL-S displacement d and use the binary marker
isolation construction. In a good U-cell the free sites are exactly
F=U intersect {eta=0}=(T union (T+d)). The local argument above still applies:
compatibility of active anchors is their full composite compatibility, and
the third-window contradiction uses convexity of S, not convexity of T.
Only the endpoint-pair mixed moment changes, and that pair specifies all
free letters in F. Its unnormalized gain is one. Hence the conditional
local ratio is at least 1+r^(-|F|)(1-t)^2, and in particular at least
1+r^(-|U|)(1-t)^2. Residual conditional comparison was just established.

For the uniform good-cell bound, fill zero sites of the isolating marker
patch with a constant non-c symbol avoiding the restriction u|T; choose a
different symbol if that restriction is constant. The resulting inner patch
avoids u. The q>=3 universal u collar extends it across B against every
exterior, again giving probability at least q^(-b0). The previous finite
packing and Jensen argument applies unchanged, even though the unweighted
marker is Bernoulli(1/q) rather than an independent equally sized product
factor. Thus the same bounds hold with r=q-1 (and the sharper |F| in place
of |U| when desired). QED.

At q=3 this permits changing two letters while preserving the third letter's
position set, after independent equal-G relabeling of the endpoint patterns.
Exactly 30,756 of the 30,788 critical ternary covers through area 25 have
such representatives; 32 have none, as independently exhaustively checked. This theorem
does not supply a binary q=2 endpoint mechanism: its nontrivial fiber needs
at least two symbols in addition to the marker symbol.

## Strict entropy along axis line sweeps at q>=3

**Theorem (working proof).** Let u,v be patterns on an h by w rectangle,
h,w>=2, in the uniform q-symbol background, q>=3. Suppose G(u) is properly
contained in G(v), and the row-sweep criterion
min r(M_v(d))<=max r(M_u(d)) holds for every d outside G(v). Then
h_q(u)<h_q(v). Reflections and transposition give the same conclusion for
any axis direction. This upgrades the previously proved strict finite-pressure
comparison to a strict hard endpoint without a spare or shared letter.

Choose a maximal-new-overlap displacement d and orient its endpoints so that
d=(a,b), a,b>=0, reflecting columns if necessary. Let U=S union (S+d),
R=[0,h+a-1] times [0,w+b-1], and let y be the last old conflict row of the
endpoint pair. Then a<=y<h. Put g=4 max(h,w), B=R+[-g,g]^2,
W=w+b+2g, H=h+a+2g, and L=HW. One explicit bound, with natural logs, is

    h_q(v)-h_q(u)
       >= -(1/W) log(1-q^(-L) q^(-2W)/(1+q^(-2W))) > 0.

For 0<=t<=1 the same pressure bound holds with q^(-2W) replaced by
c_t=q^(-2W)(1-t)^2; t=0 denotes entropy. The very small constant is a
uniformity certificate, not an estimate of the numerical gap.

**Uniform collars at a moving interface.** At sweep stage k, an event at
anchor z requests v(s) at physical rows z_1+s_1<=k, and u(s) at rows >k.
This is an inhomogeneous cylinder family, not a single homogeneous SFT.
For a pattern w let C_w be the colors absent from all its monochromatic
faces. At most two face colors occur, so C_u and C_v are nonempty.
Choose (c_v,c_u) in C_v times C_u avoiding the two pairs

    (v(0,0),u(h-1,0)), (v(0,w-1),u(h-1,w-1)).

Such a choice exists. If both sets have size at least two, at most two of
at least four choices are excluded. If one set is a singleton, q=3 and its
pattern has two oppositely placed monochromatic faces in different colors.
Every corner lies on one of those faces, so none of its corners is the
remaining color; neither forbidden pair can then eliminate a choice.

The background c_v above and including row k, c_u below row k, avoids every
entire face of every translated hybrid event. Horizontal faces see a color
in the corresponding C_w. A vertical face lying in one phase is excluded
in the same way; a face crossing the interface would need the excluded pair
at its top and bottom endpoints. Thus the rectangular insertion argument
works at every stage, uniformly in k.

For an inner rectangle containing the interface, fill its collar with this
background. At the two top outer corners use the homogeneous v corner
patches from the universal-collar construction; at the bottom corners use
the u patches. Width g separates all corner patches from each other, the
inner rectangle, and the interface. The original corner and one-side proofs
apply away from the interface; vertical side crossings near the interface
are blocked by the new face-avoiding background. Hence every hybrid event
meeting B but not contained in its inner rectangle has a pinned mismatch.
This holds for arbitrary inner and outer letters. Collars for stages k-1
and k agree everywhere except row k, using the SAME chosen color pair.

**A local row witness.** Translate so that the current row is y. Immediately
before its update rows <y request v and rows >=y request u; afterwards rows
<=y request v. The endpoint pair becomes compatible exactly at this step.
Fill U with this post-update pair. Among S-windows contained in U there is
no compatible third occurrence: every other pair has greater overlap and
therefore has the same compatibility in u and v, while full-S convexity
excludes a common neighbor of the endpoint pair. The stage compatibility
graph lies between the endpoint graphs, so the same exclusion holds at
this hybrid stage.

Fill the top-right gap R minus U with a color different from v's top-right
corner, and the bottom-left gap with a color different from u's bottom-left
corner. Both gaps are OFF row y: the first lies above a<=y, the second below
h>y. Every window contained in R but not U has a mismatching corner there.
Add the post-update universal collar. The resulting full B filling has
exactly the two designated post-update occurrences among all events meeting
B, with arbitrary exterior. Let eta be its restriction off the current row.
A cell is good when its off-row letters equal this eta.

For a good cell, the active cylinders on the current row have length w.
Among those contained in U, only endpoint compatibility changes at this
step, and no active third cylinder is compatible with both endpoints after
it. The endpoints together specify all w+b row sites of U. Thus their local
unnormalized count-polynomial difference is exactly (t-1)^2. In particular
there is at least one post-update row filling avoiding all these local
cylinders, even if no pre-update avoiding row filling exists. Inserting such
a row into the patch and its post-update collar produces a filling of B
with NO post-update occurrence meeting B. Therefore, under the post-update
Gibbs measure at any 0<=t<=1,

    Pr(cell good | exterior B)>=q^(-L).

This uses the full finite-volume conditional law. Cancel the exterior-only
weight; the chosen filling has weight one, every other filling has weight
at most one, and there are at most q^L fillings. At t=0 the same argument
holds for every exterior of positive avoidance probability.

**Conditional row amplification.** Fix all letters off the current row.
Events not touching it are identical before/after; if any occurs at t=0,
this outside configuration has zero post-update mass and can be omitted.
The remaining active events are equal-length interval cylinders. Order them
by (horizontal start, vertical anchor), and let i,j be a good cell's endpoint
indices. Every other window contained in U lies between its endpoints in
this order. Write Q_-(s),Q_+(s) for the pre/post-update row partition functions
restricted to events with horizontal start <=s.

The ordered-interval resolvent identity is a sum of nonnegative terms. Its
term at the new pair i,j is

    q^(-|I_j minus I_i|)(t-1)^2 F_i^+ R_j^-,

where F_i^+ includes occurrence i and excludes/reweights preceding events,
and R_j^- is the suffix partition function conditioned on occurrence j.
If the good B-cell occupies columns l,...,r, r-l+1=W, prescribe its complete
post-update pair row. It realizes i and blocks every preceding event meeting
B. Hence

    F_i^+ >= q^(-W) Q_+(l-w).

For the suffix, use the pre-update collar row and put the old request of
endpoint j on its interval. The off-row pinning is unchanged. There is no
contained-U event after j; all contained-R events outside U are blocked by
the fixed gaps; every crossing event is blocked by the pre-update collar.
The remaining row letters inside U may be arbitrary. Thus

    R_j^- >= q^(-(W-w)) Q_-(starts > r).

These prescriptions need not avoid the events between i and j, since those
are absent from the respective factors. Also |I_j minus I_i|<=w. Restricting
the whole row to any prefix s and deleting its middle events gives

    Q_-(s) <= Q_-(l-w) Q_-(r < starts <=s).

Therefore, when Q_-(s)>0, its ratio D(s)=Q_+(s)/Q_-(s) satisfies

    D(s)>=1+c_t sum_{good pairs j with start(j)<=s} D(l_j-w),
    c_t=q^(-2W)(1-t)^2.

Tile horizontally by adjacent disjoint B-cells. Each earlier endpoint j
satisfies start(j)<=l_next-w, since its whole interval lies in its own cell.
The geometric-prefix induction gives Q_+/Q_->=(1+c_t)^J, where J is the
number of good cells. If Q_-=0 the equivalent inverse inequality is still
valid whenever Q_+>0. No positive old conditional avoidance assumption was
made.

Average the inverse inequality under the POST-update Gibbs measure mu_+:

    Z_-/Z_+ <= E_mu_+[(1+c_t)^(-J)]
             = E_mu_+ product_cells (1-delta_t 1{good}),
    delta_t=c_t/(1+c_t).

Each good event depends only on its B-cell, and the cells are disjoint.
Successively condition on the exterior of a cell. The uniform lower bound
q^(-L) gives an upper factor 1-q^(-L)delta_t each time. Hence
Z_-/Z_+ <= (1-q^(-L)delta_t)^M for M horizontal cells. This multiplicative
upper bound uses conditional probabilities, not independence or Jensen in
the wrong direction.

Finally multiply the successive physical-row comparisons in a large free
box. All but a bounded number of rows admit the witness, and horizontal
cell density tends to 1/W. Dividing log partition functions by volume proves
the displayed pressure and entropy bounds. Transposition and reflections
preserve the statement. QED.

**Audit and rectangular consequence.** `row_entropy_checks.py` checks the
inhomogeneous collar for all 6,561 ordered pairs of ternary 2x2 patterns,
including pairs not satisfying agreement inclusion, and 250 seeded pairs
on larger rectangles at q=3,...,5. The 32 critical ternary covers lacking
any shared-letter representatives all have axis certificates. For each,
the checker verifies the adjacent collars agree off the updated row,
the local row-law difference is (t-1)^2, the new pair extension has exactly
its two designated occurrences, the old suffix row blocks every later event,
and a post-update avoiding filling works against every exterior. These are
literal pinned-mismatch tests and exact integer counts, not strip numerics.

`ternary_layer_checks.py` independently replays all 30,756 shared-letter
line certificates and exhausts all endpoint component colorings for the
32 failures. Together with the new axis theorem and cover mergers when the
upper class has a binary representative, this proves STRICT entropy order
for ALL patterns on EVERY 2D rectangle of area at most 25 at q=3. The
shared-layer q=4 argument and cover-merger bounds handle all larger q.
Thus the entire area-25 strict-entropy problem is now reduced to q=2.
The unbounded-size rectangular conjecture remains open.

**A precise failed structural conjecture.** Not every ternary cover can be
represented by patterns sharing one letter-position class. On 3x8 let the
forced-component array be

    A B C B A B D B
    B C B A A B C B
    B D B A B C B A.

For G with positive directions (1,-7),(1,4),(2,-7),(2,-5),(2,3),(2,7), these
are exactly its four forced components. Missing directions require every
pair of component colors to differ except B,C; with three symbols the ONLY
possible equality partition is {A},{B,C},{D}. Adding direction (2,-3)
produces H, forcing A=D and requiring B,C distinct. Its only ternary equality
partition is {A,D},{B},{C}. The inclusion is a cover because only one
direction was added. The two partitions have no common block, even after
alphabet permutations. This is a shared-letter proof-route obstruction;
the axis theorem above proves its strict entropy comparison at q>=3.

## A chordal sufficient condition for the entire field

**Proposition (working proof).** Consider two families of prescribed cylinder
events on the same finite supports S_1,...,S_n under uniform independent
q-letter coordinates. Let X,Y be their indicator vectors and H_X,H_Y their
pair compatibility graphs, including edges between disjoint supports. Suppose
H_X is a subgraph of H_Y. If H_X[O] is chordal for every clique O of H_Y, then
X<=_sm Y. The same sufficient condition applies to any Bernoulli vectors with
the mixed-moment identity used in the proof.

Proof. Identify a binary vector with its set of ones. For a set function f,
write fhat(J)=sum_{A subset J}(-1)^{|J\A|} f(A). Simultaneous cylinder
occurrence is equivalent to pair compatibility, and a compatible family J
has probability q^{-|union_{j in J} S_j|}. Thus

E f(X) = E [ sum_{J clique in H_X[Y]} fhat(J) ].

For a chordal graph H on vertex set O, put T_H f equal to the clique sum in
brackets, INCLUDING the empty clique. We claim T_H f<=f(O) for supermodular f.
Use a perfect elimination ordering: choose a vertex v with clique neighborhood
B in H, then eliminate it. The cliques containing v contribute exactly
f(B union {v})-f(B), by Mobius inversion. By induction and supermodularity,

T_H f <= f(O\{v})+f(B union {v})-f(B) <= f(O).

The empty graph starts with T_H f=f(empty). Since Y is always a clique of
H_Y, averaging proves the assertion. QED. The standard equivalence between
chordality and a perfect elimination ordering is the only graph-theoretic
input; see Berry--Simonet, https://arxiv.org/abs/1610.09623. The displayed
induction proves the particular supermodular inequality needed here.

**Theorem (working proof).** Compatibility inclusion orders the indicator
field of ANY family of at most FIVE prescribed cylinders, on arbitrary fixed
supports and for every q>=2.

First, if H_Y is complete, its requests agree at every shared coordinate.
Merge all old requested symbols at each coordinate to its common target
symbol. This proves the assertion by coordinate-dependent merging, for
arbitrary n in the complete-target case. If n<=4 and H_Y is not complete,
every H_Y clique has at most three vertices, and the chordal criterion applies.

For n=5 and noncomplete H_Y, every clique of H_Y has size at most four. A
failure of the chordal criterion is therefore an induced old four-cycle O
which is a new clique. We show how to fill one of its diagonals p,r by
one-coordinate switches, preserving old compatibilities and gaining only
target compatibilities. There is at most one event w outside O. Put
A={i in O: iw is an edge of H_Y}, with A empty if w is absent.

At a coordinate where p,r request different symbols, neither of the other
two O-events can be incident: it is old-compatible with both p and r. Thus
only p,r,w can be incident there. The only delicate case is when exactly one
of p,r belongs to A, say p in A and r not in A. Choose a coordinate y where
the target requests of r,w conflict. No event of A is incident at y, because
all O-events are target-compatible. The O-events incident at y lie in
O\{p}, contain r, and induce a connected subgraph of the old four-cycle.
Their old requests are therefore one common symbol. If necessary, change
w's request at y to a different symbol. This changes NO global compatibility:
all potentially affected pairs join w to O\A and were already incompatible.
It supplies a conflict between w,r at a coordinate not incident to p.

Now process every coordinate where p,r still disagree, making their requests
equal there. If w is absent, choose either of their symbols. If w is incident
and at least one of p,r is in A, choose w's current requested symbol. If w is
incident and neither is in A, choose a symbol different from w's request.
Each switch preserves all current compatible pairs and introduces only target
edges. In the one-in-A case the prepared conflict at y prevents the unwanted
edge wr; it is untouched because p is not incident there. In the neither-in-A
case existing local conflicts with w persist. The other cycle diagonal is
untouched throughout. Hence these switches strictly enlarge the compatibility
graph inside H_Y, and increase the field in supermodular order by the switching
lemma. Repeat whenever the chordal criterion fails. The graph strictly gains
edges each time, so the procedure terminates; apply the chordal criterion to
the remaining comparison. QED.

The construction needs at most one event outside the four-cycle. With several
outside events their requests at a shared coordinate can conflict, so the
same symbol choice is not justified. No general-n conclusion follows.
`cylinder_order_checks.py --switches` tests the construction on 450 seeded
support systems at q=2,3,4, with 5,000 sampled request families per system.
All 669 comparisons needing a cycle-filling move passed, including all 786
individual coordinate changes. This checks the construction's implementation;
the all-support/all-q conclusion comes from the written proof.

**Deleting an event with an unchanged old-neighbor graph (working proof).**
Let G subset H be the old and new compatibility graphs. Suppose no new edge
has both endpoints in the OLD neighbor set N_G(i) of an event i. If the
comparison holds after deleting i, it holds for the whole field, for any q>=2.

Write X,Y for the old and new fields with coordinate i removed, and p for the
common probability of event i. Conditional on the old event i, all coordinates
outside N_G(i) vanish. Conditional on the new event i, the restricted field on
N_G(i) has exactly the same law: its compatibility graph is unchanged and
all conditional mixed moments are either zero or q^(-|union S_j minus S_i|).
The new conditional field can therefore be coupled above the old one by
keeping these coordinates identical and adding any newly possible neighbors.
For a supermodular f, the increment f(F union {i})-f(F) increases with F.
Use this conditional stochastic comparison in

E f(X union {i:I_i=1}) = E f(X)
  + p E[f(X union {i})-f(X) | I_i=1],

and use the assumed background supermodular order for the first term. QED.
Full-support deletion is a special case: the old neighbors of a full-support
event form a clique. The earlier private-atom proof also gives the weighted-
count version without assuming background field order.

**Completing a clique against an existing clique (working proof).** Partition
the event indices into O and T. Suppose G subset H, H[O] is complete and G[T]
is complete. Then compatibility inclusion increases the whole field for every
q>=2, by a sequence of one-coordinate switches and an equal-law replacement.
Either part may be empty.

At a physical coordinate, all incident T-events have the same CURRENT request,
since they form an old clique; their requests will be left fixed. In the target,
all incident O-events have one common request, and all incident T-events have
one common request. For every missing target edge ot with o in O,t in T,
choose a target conflict coordinate x. At x, EVERY incident O/T pair is a
target nonedge. Merge all incident O-requests to one symbol different from the
common current T-request. This preserves all current compatible pairs and
creates only edges within O, which are target edges. It also records a conflict
at x for every incident cross pair. Process these finitely many witnesses.

Now process all physical coordinates. Merge the incident O-requests to the
common T-request if some incident cross pair is an H-edge; otherwise choose
a different symbol when T-events are incident, or any symbol when none are.
Previously recorded conflicts remain: their coordinates have no incident
target cross edge. Every current cross compatibility is preserved, and every
new cross compatibility is a target edge. At the end, O is a clique and the
cross edges are exactly H's cross edges; T was a clique throughout. Thus the
final graph is H, and its field law is the target law. Apply the one-coordinate
switching lemma at every step. QED. In particular this handles every target
graph that becomes complete after deleting one vertex.

**Six-cylinder field theorem (working proof).** Compatibility inclusion orders
the whole field of ANY family of at most SIX prescribed cylinders, on arbitrary
fixed supports and every common uniform alphabet q>=2.

Use the five-cylinder theorem as the base case. If the chordal criterion holds,
there is nothing to prove. Otherwise there is an induced old cycle O of length
at least four whose vertices form a target clique. Complete target graphs are
already handled by letter merging. For a noncomplete target on at most six
vertices, O has length four or five. At length five, the remaining part T has
at most one vertex, so the preceding clique-completion theorem applies.

At length four, T has two vertices a,w. If aw is an old edge, the same theorem
applies. Suppose it is not. If no new edge lies inside N_G(a), delete a using
the preceding deletion lemma and the five-cylinder theorem. Otherwise
N_G(a) lies in O and contains the endpoints p,r of an old cycle diagonal.
Thus a is old-compatible with both p and r. We can fill this diagonal by the
five-cylinder construction, treating w as its only active outside event:

- At a coordinate where p,r disagree, the other two cycle events and a are
  absent, because they are all old-compatible with both p and r.
- If w has exactly one target neighbor among p,r, call it p, and choose a
  target conflict coordinate y for w,r. The event p is absent at y. All
  incident cycle events have the common old request of r, since the induced
  cycle with p removed is a path centered at r. If a is incident, it too has
  r's request. When necessary change w's request at y to a different symbol.
  This creates no non-target compatibility: all affected incident cycle
  events are target nonneighbors of w, and the old nonedge aw stays a nonedge.
- At every remaining conflict coordinate of p,r, merge their requests to
  w's current request if w is incident and has a target neighbor among them;
  if w is incident with neither target neighbor, choose a different symbol;
  otherwise choose either p/r request. Only p,r,w can be incident, and the
  prepared conflict at y is untouched because p is absent there.

Each step preserves every current compatible pair and gains only target
edges; in particular it gains pr. Repeat this argument whenever the chordal
criterion fails. Each repetition strictly enlarges the graph inside H, so the
procedure terminates, or one of the deletion/completion cases finishes it.
This proves the comparison. QED. No statement for seven arbitrary cylinders
is inferred from this argument.
The exact structural audit in `six_cylinder_checks.py` checks all 19,926
three-state exterior-edge assignments around a filled four/five-cycle, plus
22,485 seeded comparisons on arbitrary supports at q=2,3,4. All 124,497 literal
coordinate moves preserve current compatibility and stay inside the target
graph, including 2,678 protected diagonal fills. These checks support the
construction; its all-support/all-q scope comes from the proof.

**Adding an edge with a clique outside its common neighbors (working proof).**
Let pr be a new target edge, put A=N_G(p) intersect N_G(r), and let
C=V minus (A union {p,r}). If G[C] is complete, one-coordinate switches can
add pr while preserving every current edge and adding only target edges, for
any number of events and every q>=2.

Keep all C-requests fixed. For every target nonedge between {p,r} and C,
choose a target conflict coordinate x. All incident C-events have one current
request b. All incident events in {p,r} union A have one TARGET request,
different from the target C-request. Change their current requests at x to
any symbol other than b. This preserves their current compatibilities: if an
A-event is incident, it already shares the request of each incident endpoint;
if the endpoints disagree, no A-event can be incident. No A/A local agreement
is lost or gained. Every incident cross pair to C is a target nonedge, and
the change records a conflict there. These prepared coordinates have at most
one endpoint incident or have the two endpoint requests equal afterward.

Now merge p,r at each remaining conflict coordinate. No A-event is incident.
If C is incident and some endpoint/C pair is a target edge, use the common
current C-request; if C is incident and no such target edge exists, use a
different symbol; with C absent use either endpoint request. All prepared
conflicts survive, since those coordinates are no longer endpoint conflicts.
Thus every step preserves compatibility inside the target, and pr is gained.
QED. This generalizes the protected diagonal move; it does not assert that
every comparison has an eligible new edge.
The extension to `six_cylinder_checks.py` checks 3,289 eligible pairs in
seeded families of four through ten events at q=2,3,4, with all 7,703 literal
coordinate changes preserving the required graph inclusions.

**Localization of a field difference (working proof).** For an inclusion
H_X subset H_Y, let A contain both endpoints of each new edge and every common
H_Y-neighbor of those endpoints. Then the difference of the two field laws is
supported on indicator vectors which vanish outside A, and equals the
difference of their A-marginals there. Indeed the mixed moments differ only on
new cliques, and every new clique is contained in A. Mobius inversion gives
the asserted support statement. Hence supermodular comparison reduces to the
two A-marginals. In particular |A|<=6 now proves the full-field comparison, for
arbitrarily many cylinders and every uniform q>=2, by the six-cylinder theorem.

**Alphabet bound for arbitrary cylinders (working proof).** At physical
coordinate x, restrict the old compatibility graph G to the incident event
indices, and let c_x count its connected components. If q>=max_x c_x, the
whole field comparison follows from coordinate-dependent letter merging.
Indeed assign distinct requested symbols to these components at each x.
This preserves exactly G: old edges remain equal at every shared coordinate,
and an old nonedge has an original conflict coordinate whose two events
cannot belong to one component. The target requests are constant on every
component, because G subset H, so the canonical requests merge to the target
at each coordinate. Equal compatibility gives equal old field laws. QED.
In particular q>=n settles every n-cylinder comparison on arbitrary supports.

**Reduction to exceptional physical coordinates (working proof).** Let
K={x:c_x>q}. Give distinct canonical symbols to the old compatibility components
at every coordinate outside K, and retain the original requests on K. This
preserves the entire old compatibility graph by the preceding argument. Outside
K the target is now a coordinate-dependent letter merger. First replace all
requests on K by their target values, leaving the canonical exterior fixed.
Conditional on that exterior, the active cylinders have compatibility inclusion
on K: a locally compatible active pair was globally old-compatible. Thus any
valid field comparison theorem for cylinders on K proves this first step.
Then perform the coordinate mergers outside K. In particular |K|<=2 always
proves the whole field comparison, for every q>=2.

**Three-coordinate binary field theorem (working proof).** Compatibility
inclusion orders ANY number of prescribed cylinders on at most THREE fair
independent bits. Consequently the same holds for binary prescriptions on
three coordinates under every uniform q>=2, by alphabet extension.

First remove full-support events by the private-atom deletion principle;
it suffices to prove the comparison for proper supports. Duplicate old events
remain duplicate in the target and can be grouped into one indicator, since
duplicating a coordinate preserves supermodularity of the test function.
Empty-support events are constant. The remaining distinct event types are the
six one-coordinate prescriptions and twelve two-coordinate prescriptions.

For a physical coordinate x and requested bit b, all incident events requesting
b form a connected compatibility graph if a singleton event occurs there, or
if both pair supports {x,y} and {x,z} occur there. Otherwise this bit class
can have more than one component only when it consists of both prescriptions
(x=b,y=0) and (x=b,y=1) on ONE pair support {x,y}, with no such-bit event on
{x,z}. Call this a private row of that pair, just within this proof.

Suppose c_x>2. Some bit class is disconnected, so choose its private row on
{x,y}. This row supplies both y-bit values. If also c_y>2, a disconnected
y-bit class must therefore be private to the SAME pair {x,y}; it cannot be
private to {y,z}, which would exclude the already present {x,y} event of that
bit. Say these private bits are b at x and d at y. Every event on {x,z} then
requests x=1-b, and every event on {y,z} requests y=1-d. For each fixed z-bit,
there is at most one distinct event on each of these pair supports, and any
such events, together with the possible singleton at z, are mutually compatible.
Hence c_z<=2. Three exceptional coordinates are impossible: choosing a private
row at one of them would give the preceding contradiction at the third.
Thus |K|<=2, and the exceptional-coordinate reduction finishes the proof for
proper supports. Restore the full-support events. QED.

Conditioning outside a three-coordinate block gives a THREE-coordinate switching
lemma for binary cylinder families on arbitrary supports. Its conclusion also
lifts to every uniform q>=2 when all prescriptions use the same two symbols.
Combining this with the exceptional-coordinate reduction, the binary field
comparison holds whenever at most three physical coordinates have c_x>2.
Thus a minimal counterexample on four bits must have c_x>2 at ALL four
coordinates, after deleting full-support events and grouping duplicates.
No corresponding theorem for arbitrary ternary prescriptions on three physical
coordinates is claimed.
The exact component audit in `private_atom_checks.py` exhausts all 262,144
subsets of the 18 proper binary cylinder types. The numbers with zero, one,
two, or three exceptional coordinates are respectively 180,856, 70,416,
10,872, and zero. A further 3,000 seeded comparisons verify the literal
canonical replacement, exceptional-coordinate switch, and final mergers.

**Eight cylinders on four bits (working proof with exhaustive structural
certificates).** Compatibility inclusion orders every family of at most EIGHT
binary prescribed cylinders on at most four physical coordinates in
supermodular order. Alphabet extension gives the same conclusion for these
binary prescriptions over every uniform q>=2. This does not settle seven
events on arbitrary supports.

Induct on the number n of events, starting with the six-event theorem. For
n=7 and then n=8, remove full-support events by private-atom extension, group
duplicate old events, and discard constants. A simplicial vertex of the old
compatibility graph can be deleted, since its old neighbor set contains no
new edge. If any of the four coordinates has c_x<=2, the preceding
three-exceptional-coordinate corollary applies. Thus only n-element subsets
of the 64 distinct nonconstant proper binary cylinder types remain, with no
simplicial old vertex and c_x>2 at every coordinate.

The exhaustive source enumeration and symmetry reduction give:

| n | All source families | No simplicial vertex | All four exceptional | Symmetry orbits | Interval orbits | Other target graphs |
|---|---:|---:|---:|---:|---:|---:|
| 7 | 621,216,192 | 9,085,584 | 2,328 | 11 | 6 | 768 |
| 8 | 4,426,165,368 | 84,199,082 | 93,651 | 298 | 97 | 48,425 |

Symmetries are coordinate permutations and independent bit complements.
An interval orbit has a coordinate ordering in which every physical support
is an interval, so the general interval theorem applies. For each other
source, every target request is constant on each old component at each
coordinate. Enumerating the 2^(c_x-1) binary assignments modulo complementation
at coordinate x therefore gives ALL target graphs.

Process these realized target graphs H by increasing number of edges. A
comparison G to H is settled directly if it satisfies the chordal criterion,
permits deletion to n-1 events, or has complete target (coordinate mergers).
Otherwise the checker finds an already-settled realized graph K with
G properly contained in K properly contained in H, and verifies K to H by
one of these criteria or the three-exceptional-coordinate theorem. Equal
compatibility gives equal laws on the fixed supports, so the intermediate
representative need not equal the original endpoint words. Every graph is
reached. For n=7 the counts are 5 equal, 607 chordal, 126 deletion, 5 complete
and 25 chain comparisons. For n=8 they are 201, 28,330, 17,155, 201 and 2,538.
The n=8 deletion step uses the just-established n=7 result on four bits.

The proof uses no LP certificates. Thirty earlier direct rational comparisons
(492 elementary terms) are retained as independent checks. The Python replay
verifies every retained source property, symmetry orbit, target graph and
structural step; `--enumerate` also rebuilds and runs both complete C++ source
enumerations. Both have passed. Code/data: `binary_four_coordinate_checks.cpp`,
`binary_four_coordinate_checks.py`, `binary-four-coordinate-evidence.json`.

**Ternary stress test, exact finite result.** The binary component bound fails
for ternary proper cylinders. A source family on coordinates x,y,z is

xy=00; xz=10,11; yz=01,02,10,12; x=0,2; y=2.

These ten events have c_x=c_y=c_z=4. Enumerate all partitions of the four old
components at each coordinate into at most three requested-symbol classes:
there are 14^3 assignments up to independent target-symbol permutations and
2,386 distinct target compatibility graphs. EVERY one has field supermodular
order above this source. The chordal criterion handles 1,009; deletion followed
by the component/interval/six-event criteria handles 1,333; the remaining 44
have exact rational certificates totaling 2,106 elementary supermodular terms.
`ternary_three_coordinate_checks.py` replays all target classes and certificates
using only the standard library. The exact laws count the 27 equally likely
physical outcomes. Alphabet extension lifts each of these comparisons to
every uniform q>=3 with the same prescriptions. This is one fixed source
family, not an arbitrary-ternary three-coordinate theorem.

A minimal counterexample to the UNIVERSAL field conjecture must therefore
have n>=7 events, at least three physical coordinates (at least FOUR for binary
prescriptions, by the new theorem), q<max_x c_x<=n, and
no full-support event. Every vertex's OLD neighbor set must contain a new
edge, by the deletion principle. Its localization set must contain every
vertex, and its target graph is noncomplete. These are necessary restrictions
on a counterexample, not a conjecture that every surviving graph is realizable
or fails the comparison. They suggest testing seven-event families with
genuinely overlapping proper supports before increasing rectangle sizes again.

The chordal criterion is not universal, even on covers of the realizable
compatibility poset. Here is an explicit binary example on physical coordinates
0,...,4. Rows list each event support and the corresponding requests in that
support's displayed order:

| i | S_i | old | new |
|---|---|---|---|
| 0 | 0,2,4 | 000 | 000 |
| 1 | 2,3,4 | 100 | 001 |
| 2 | 0,2,3 | 000 | 001 |
| 3 | 1,3,4 | 010 | 000 |
| 4 | 0,1,4 | 011 | 010 |
| 5 | 0,1,3 | 011 | 011 |
| 6 | 0,3,4 | 011 | 010 |

The new compatible edges are exactly 04,06,25,26. The old graph induced by
{0,2,4,5} is a chordless four-cycle, and this set is a clique in the new graph.
Exact enumeration of all 27,209 realizable graphs on these supports finds no
strict intermediate graph. Nevertheless this comparison DOES satisfy
supermodular order. Thus failure of this structural criterion is not a
counterexample to the ordering conjecture.

`cylinder_order_checks.py` exhausts the compatibility classes for two fixed
support systems on five independent bits. For six triple-support cylinders,
all 14,045 covers of the 2,826 classes satisfy the chordal criterion; hence all
193,172 strict inclusions increase the whole field in supermodular order.
For the seven-cylinder supports above there are 4,940,162 strict inclusions
and 167,038 covers; exactly twelve covers fail that criterion. Each has an
additional exact certificate expressing 32 times the
expectation difference as a nonnegative rational combination of elementary
supermodular differences. All other covers satisfy the chordal criterion.
Each of the twelve exceptional covers also admits a two-coordinate binary
switch between equal-law representatives, giving a structural explanation
through the new switching corollary. The exact supermodular identities remain
independent checks of these finite comparisons.
Alternatively, each of those twelve covers has a localization set of size five,
so the localization lemma gives their comparison for every uniform q>=2.
Two-coordinate switch paths are not universal either. The same certificate
file stores another seven-cylinder cover on five fair bits for which exhaustive
local-cut enumeration finds no two-coordinate switch between ANY representatives
of its two classes. Its new edges are 34 and 56, and its localization set is
{3,4,5,6}. Its field comparison is therefore a four-cylinder comparison; an
independent certificate is the sum of four elementary differences
Delta_34({5})+Delta_34({6})+Delta_56({3,4})+Delta_56(empty), where
Delta_ij(B)=f(B)+f(B union {i,j})-f(B union {i})-f(B union {j}).
The cover property rules out chains of two-coordinate switches too. This is
another obstruction to a proposed proof route, while the inequality survives.
The certificate file `cylinder-order-evidence.json` stores the full support
data and the twelve exceptional certificates. Default verification needs only
the standard library; optional `--solve` uses numerical linear programming
to discover coefficients, which are then checked as exact rational identities.
These are finite-field results on these specific support systems.

Realizability cannot be discarded in an attempted general proof. For example,
take four supports consisting of the three-element subsets of four physical
coordinates, and let Y indicate their all-zero cylinders. Under fair bits,
Y is all four events with probability 1/16, each singleton with probability
1/16, and empty with probability 11/16. Form an abstract law X by putting
probability 1/16 on each edge of a four-cycle and 12/16 on the empty set.
Its mixed moments are exactly those of Y with all noncliques of that cycle
deleted. Both laws are nonnegative and have the same marginals, but

E(|X|-1)_+=4/16 > 3/16=E(|Y|-1)_+.

The old four-cycle cannot be realized by cylinder prescriptions on those
supports: at each physical coordinate the three incident event indices induce
a connected path, so their old compatible pairs force all three requests to
agree. All pairs would then be compatible. This excludes an argument based
only on moment deletion and positivity of the two endpoint laws; it does not
refute the cylinder-field or pattern conjecture.

## Agreement inclusion when there are enough alphabet symbols

**Theorem (working proof).** Form a graph on S by joining distinct a,b whenever
a-b belongs to G(u), and let c(u) be its number of connected components. If
q>=c(u) and G(u) subset G(v), then N_u<=_cx N_v on every finite anchor set.
If the inclusion is strict, P_u(beta)<P_v(beta) at every finite beta!=0.
If S is lattice-convex, the strict inclusion also gives h(X_u)<h(X_v).
In particular, strict agreement-inclusion entropy ordering holds on every
finite lattice-convex support whenever q>=|S|.

Proof. The original u is constant on every component: each graph edge is an
actual equality of two sites. Label the components by distinct alphabet symbols
to obtain u*. This is possible when q>=c(u). Every old agreement d remains an
agreement of u*, since all pairs of sites compared at displacement d are graph
edges. A displacement missing from G(u) compares some pair with different
u-symbols; those sites must lie in different components and still differ in u*.
Thus G(u*)=G(u), so their occurrence fields have identical laws under the
uniform q-symbol background.

Because G(u) subset G(v), v is also constant on every component. Hence v=f
composed with u* for a map of alphabet symbols. Apply the letter-merging convex
comparison, its strict finite-coupling version, and (for lattice-convex supports) its strict
entropy version. QED.

This argument uses the original q-symbol reference distribution throughout;
it is not an enlargement of the alphabet followed by an unjustified transfer
back to smaller q. If q<c(u), the canonical refinement is not a pattern over
that alphabet. This is a precise obstruction to this proof, not a counterexample
to the small-alphabet conjecture. The original letter-merging theorem can still
apply in individual pairs even when q<c(u).

**Corollary (working proof, paired labels).** If the number of distinct pairs
(u(s),v(s)), s in S, is at most q, agreement inclusion again reduces to merging:
label each occurring pair distinctly to obtain w. Then G(w)=G(u) intersect
G(v)=G(u), while v is a letter merger of w. All conclusions of the theorem
follow. Thus for patterns using at most r and t symbols respectively, q>=rt
suffices independently of |S|. In particular, for BINARY patterns the entire
inclusion comparison, including strict entropy on lattice-convex supports,
holds over every uniform background with q>=4. The remaining background sizes
for this restricted class are q=2 and q=3; this does not settle patterns using
all q symbols.

The agreement sets themselves form a partially ordered collection of forced
partitions. For a realized G, let Pi(G) be the components of the graph above.
A displacement belongs to G exactly when all the site pairs it compares lie
in common blocks of Pi(G): one implication uses graph edges, the other uses
constancy of the original pattern on each component. Consequently G subset H
if and only if Pi(G) refines Pi(H). Strict inclusion strictly decreases the
number of components. This identifies the inclusion poset without imposing
an entropy order on incomparable elements. The alphabet restriction determines
which forced partitions can be realized: merging component labels is allowed
only when it creates no additional full-displacement agreement.

The core remaining inclusion problem is therefore the low-alphabet case,
especially binary rectangles. The lower extreme with G empty has c(u)=|S|,
but the separate Shearer argument below covers its non-strict comparison even
when q is much smaller.

### Alphabet-dependent agreement posets

For fixed S, let Gamma_q(S) be the collection of agreement sets realizable
over at most q symbols, and let a_S(G) be the smallest such q for a given G.
Under a fixed uniform q-symbol background, the finite field law and entropy
depend only on G. The main entropy conjecture therefore says that h_q is
strictly increasing on the inclusion poset Gamma_q(S). As q grows, additional
agreement classes can appear. Equal-G replacement is legitimate at the same
ambient q only when the replacement uses at most q symbols.

There is an exact finite description of realizability. Given a set D of
displacements, join all site pairs separated by an element of D, and let
Pi(D) be the resulting component partition. Define cl(D) to contain every
displacement whose compared pairs all lie inside components of Pi(D).
Then the agreement sets over arbitrary alphabets are exactly the fixed points
of cl: distinct labels on the components realize cl(D), while any word with
agreements D is constant on those components. Starting from cl(empty) and
repeatedly adding one displacement and closing enumerates every fixed point.
No bound on the original alphabet is imposed by this enumeration.

For a fixed point G, color the components of Pi(G). For every d outside G,
at least one compared pair of DISTINCT components must receive distinct
colors. These clauses give an exact finite coloring problem for a_S(G).
Equivalently, every arbitrary-alphabet G is an intersection of binary
agreement sets: take the binary indicator of each letter class of a realizing
word. This does not mean the intersection itself is binary-realizable.

**Cover-merger lemma (working proof).** Fix q. Suppose G is covered by H in
Gamma_q(S), and H has a representative v using fewer than q symbols. Then
there is a q-symbol representative u of G such that v is a letter merger of
u. More generally the same holds for a cover in the full alphabet-independent
poset whenever q>=a_S(H)+1.

Proof. The partition Pi(G) properly refines Pi(H). Choose a Pi(G) component
C lying strictly inside one Pi(H) component. The representative v is constant
on C. Replace that letter on C by one new symbol, obtaining u. Every G
agreement survives because C is a whole Pi(G) component; no new agreement
outside H is possible because merging the new symbol back gives v. Hence
G subset G(u) subset H. Some H edge crosses from C to the rest of its Pi(H)
component; at its displacement u disagrees, so G(u) is properly contained
in H. It lies in Gamma_q(S), and the cover assumption forces G(u)=G. QED.

**Consequence (working theorem).** Put a(S)=max_G a_S(G), over all realizable
agreement classes on S. For EVERY uniform background q>=a(S)+1, agreement
inclusion implies finite field supermodular order. Under proper inclusion it
gives strict pressure at every finite beta!=0, and, when S is lattice-convex,
strict entropy. Take a saturated chain of covers in the full poset and use
the cover-merger lemma at every step, with equal-G replacement between steps.
Each comparison uses the SAME q-symbol background. More locally it suffices
that a_S(K)<q for every upper class K on one saturated chain from G to H.

This also directs the small-alphabet programme: a comparison failure must
already occur on a cover of Gamma_q(S) whose UPPER class requires all q
symbols. Every other cover is settled by a literal one-letter merger.
The conclusion applies separately to avoidance, convex order, field order,
finite pressure, and, on convex supports, entropy, by their transitivity.
It neither assumes nor proves a universal bound on a(S).

**A universal quantitative alphabet bound (working proof).** Let
D=(|S-S|-1)/2 be the number of nonzero overlapping displacement pairs.
For every realizable G,

binom(a_S(G),2) <= D-|G_+|,

where G_+ chooses one of each pair d,-d. Indeed, in a representative using
the minimum number r=a_S(G) of letters, merging ANY pair of letters must
create a new agreement. Its old mismatches involve precisely that unordered
letter pair. A displacement cannot be newly enabled by two different pairs
of letter mergers. Assign one such missing displacement pair to each of the
binom(r,2) mergers, giving the injection and the inequality. The same argument
applies to any representative admitting no agreement-preserving letter merger.

Consequently the whole inclusion comparison is proved whenever

binom(q,2) >= D-|G(u)_+|.

To see this, take a saturated chain above G(u) in Gamma_q(S). Every upper
class K in the chain has at most D-|G(u)_+|-1 missing displacement pairs.
The preceding bound rules out a_S(K)=q. Every cover is therefore a merger
by the lemma. This gives finite FIELD order and strict finite nonzero pressure;
on lattice-convex supports it gives STRICT ENTROPY under proper inclusion.
No assumed avoidance inequality or cluster expansion enters this argument.

In particular, binom(q,2)>=D settles ALL comparisons on S. For an h by w
rectangle, D=2hw-h-w, so the explicit sufficient alphabet size is of order
2 sqrt(hw), improving the canonical-component bound q>=hw for large boxes.
For a d-dimensional box, use D=(product_i(2l_i-1)-1)/2. This is an elementary
shape-dependent bound; it does not prove a uniform four-symbol realization
bound in two dimensions.
The alphabet bound is sharp for general finite supports. Take n points
S={1,2,4,...,2^(n-1)} in Z. Their positive differences are all distinct, so
every nonzero overlapping displacement compares just one site pair. Realizing
G empty forces all n site letters to differ, and a(S)=n with
D=binom(n,2). Thus any size-independent alphabet-reduction conjecture must
use support geometry; it cannot extend to arbitrary supports even in Z.

**Exact all-alphabet small-rectangle comparison.** The closure enumeration and
component-coloring search give the following complete posets. Every class has
a representative on at most three symbols; every cover has a ternary
line-sweep certificate satisfying the literal projected-conflict criterion.

| Support | All agreement classes | Proper inclusions | Certified covers |
|---|---:|---:|---:|
| 2x2 | 7 | 13 | 9 |
| 2x3 | 19 | 71 | 35 |
| 2x4 | 48 | 302 | 111 |
| 3x3 | 72 | 634 | 179 |
| 2x5 | 100 | 946 | 266 |
| 2x6 | 206 | 2,687 | 619 |
| 3x4 | 275 | 3,950 | 884 |

Consequently every G inclusion between arbitrary patterns on these supports
has all-anchor field supermodular order at every q>=3, and strict pressure
at every finite beta!=0 under proper inclusion. At q=2 the previously
certified binary comparisons apply. Equal-G replacement followed by the
spare-symbol theorem gives STRICT entropy for EVERY q>=4. The new
shared-letter and axis-sweep theorems additionally settle q=3.

**Four symbols are necessary already in two dimensions (exact obstruction).**
On a 3x7 rectangle put

    A B A C A D A
    C A D A D A C
    A D A C A B A

with A,B,C,D distinct. Its positive-half agreement set is

G_+={(0,6),(1,-3),(2,-6),(2,-4),(2,4),(2,6)}.

The equalities at these six displacements force exactly the four displayed
components. Any other representative with this G must therefore have the
same displayed form, allowing initially coincident labels. The following
missing displacements force their pairwise inequalities:

| Missing displacement | Necessary inequality |
|---|---|
| (1,-6) | A != C |
| (1,-5) | C != D |
| (1,3) | B != D |
| (1,5) | B != C |
| (2,-5) | A != D |
| (2,5) | A != B |

For each row, all nontrivial compared pairs are precisely that component
pair; all remaining comparisons are equalities already forced by G. Thus
three colors cannot realize G, and four suffice. This refutes a universal
ternary representation lemma for 2D rectangles. It is not an entropy-order
counterexample. The exhaustive spectra show this is the FIRST possible area
among rectangles with both sides at least two: every class through area 20
has a ternary representative. On rectangles through area 25, four suffice;
the numbers of classes requiring four symbols are 2 on 3x7, 2 on 3x8,
8 on 4x6, and 4 on 5x5, with none on the other shapes. No universal
four-symbol bound is claimed.

In dimension three, the cube pattern in lexicographic site order
(A,B,C,D,D,C,B,A) already requires four symbols. Its agreements are exactly
the four body-diagonal pairs. They force equality on each antipodal pair;
the six face-diagonal directions impose all six inequalities between these
four components. The full cube spectrum has 106 classes: 1 needs one color,
71 need two, 33 need three, and 1 needs four. On 2x2x3 the corresponding
counts are 1,741,445,8, totaling 1,195. No all-comparison theorem in dimension
three at small q is inferred from these spectra. The cover-merger lemma DOES
give all comparisons, including strict entropy, for both 3D boxes at q>=5.

**Five symbols can be necessary on a higher-dimensional box (exact finite
certificate).** On S={0,1}^5 let G consist exactly of the body-diagonal
displacements, with all five coordinates nonzero. They force equality on
each antipodal vertex pair, leaving 16 components. The pattern, in binary
lexicographic site order,

    0112123012344023 3204432103212110

has exactly G and uses five symbols. Four are impossible. Missing directions
with exactly one zero coordinate require a proper coloring of the graph on
the 16 antipodal pairs, adjacent when cube vertices are neighbors. Label its
vertices x=0,...,15; its neighbors are x xor 1, x xor 2, x xor 4, x xor 8,
and x xor 15. Exhaustive enumeration gives 392 nonempty independent sets,
with size histogram 16,80,160,120,16 for sizes 1,...,5. Partitioning the
vertices into at most four independent sets gives exactly 1,190 colorings
modulo color permutation. Of these, 1,110 have a bichromatic square face,
which gives an unwanted agreement. Each of the other 80 has exactly three
unwanted displacement pairs with three zero coordinates. Thus none realizes
G. The independent-set partition audit and the literal five-color word are
in `alphabet_sweep_checks.py`. This refutes a dimension-independent
four-symbol reduction for boxes, while leaving the two-dimensional question
open. It asserts a_S(G)=5 for this class, not a(S)=5 for the entire cube.

**New uniform-alphabet consequences of the exact spectra.** The computed
bound a(S)<=4 on EVERY 2D rectangle of area at most 25 proves ALL pattern
inclusions there for every q>=5, including strict entropy and field order.
For all rectangles through area 20, and additionally 2x11 and 2x12, a(S)<=3
improves this to q>=4. This extends the earlier binary-only area-25 result to
arbitrary patterns, including patterns initially using all q letters. Equal-G
replacement is what allows those patterns to use the small representatives.
At q=4, the exceptional four-symbol upper classes remain possible obstacles
on 3x7,3x8,4x6,5x5; at q=3 the analogous upper classes are genuinely ternary.

**All-pattern finite-order theorem through area 25 (working proof with exact
certificates).** For every 2D rectangle S with both side lengths at least two
and area at most 25, every q>=3, and every G(u) subset G(v), the occurrence
field of u lies below that of v in supermodular order on EVERY finite anchor
set. Under proper inclusion pressure is strictly ordered for every finite
beta!=0. At q=2, weighted count order and the same strict pressure conclusion
hold; the 88 previously identified field comparisons remain open there.

Proof. For q=3, enumerate Gamma_3(S). Every cover whose upper class uses at
most two symbols is a merger by the analytic lemma. The remaining 30,788
covers across the 21 shapes all have exact ternary line-sweep certificates.
For q=4, every class is realizable, and only upper classes needing four
symbols require separate treatment. There are 74 incoming covers: 8 on
3x7, 10 on 3x8, 38 on 4x6 and 18 on 5x5. All have exact four-symbol
line-sweep certificates. For q>=5 all covers are mergers. Transitivity proves
the all-anchor field assertion; flag isolation and the known positive-
coupling comparison give strict finite pressure. The q=2 count assertion is
the earlier binary line/Reimer theorem. QED.

The certificate search needs only the saved minimum-alphabet representatives,
with a change of one forced component to an existing letter in the few
cases where those representatives fail. This is a successful finite search,
not a universal line-sweep theorem. Literal agreement/mismatch replay and
exhaustive intermediate-submask checks audit all 30,862 exceptional covers.
All other covers are handled analytically, so a separate enormous list of
their generators is unnecessary.

The new shared-layer, shared-letter and axis-sweep theorems now give strict
entropy for ALL patterns through area 25 at EVERY q>=3. Cover mergers handle
upper classes with an unused ambient symbol; the 74 q=4 critical covers have
binary-layer factorizations, and the 30,788 q=3 critical covers are covered
by 30,756 shared-letter certificates and 32 axis witnesses. A possible entropy
counterexample in this finite size range is therefore restricted to q=2.
This does not settle the general unbounded-size finite-order or entropy
conjectures.

`alphabet_sweep_checks.py` reconstructs every agreement closure system and
minimum alphabet through rectangular area 25 and on these two 3D supports.
It independently checks the binary layers against the earlier exhaustive
enumerations, and checks all equality partitions on supports of size at most
nine. Every saved line certificate is replayed through the original literal
agreement routine and projected mismatch sets; all covers are accounted for.
Both four-color examples also have direct K4 lower-bound checks. Evidence:
`alphabet-sweep-evidence.json`. These are exact structural certificates with
all-volume consequences from the working interval/strictification proofs.
The cover-merger construction has an additional literal audit on all 2,103
small-rectangle covers: 192/6,377/610 valid component splits use respectively
2/3/4 symbols. The general cover argument is analytic, not inferred from that
finite audit.
The alphabet-bound audit additionally checks every pair merger of every saved
minimum-alphabet representative and the disjointness of the newly enabled
displacement sets for different letter pairs.

### Three-symbol merger chains: exact certificates and an obstruction

At q=3, form a directed graph whose vertices are agreement classes and whose
edges come from merging two letters of some ternary representative. Equal-G
replacement is free. Every directed path gives convex order for ALL finite
anchor sets; on lattice-convex supports it gives strict entropy under strict
inclusion. This conclusion uses the merging proofs, rather than any passage
from a single finite-box check to infinite volume.

`merger_reachability.py --large` exhausts ternary alphabet partitions of the
support and verifies every stored edge against the original agreement routine.
`merger-evidence.json` records the results. Every binary inclusion pair is
reachable on 2x2, 2x3, 2x4, 2x5, 2x6, 3x3, 3x4, 3x5 and 4x4 supports (up to
transpose). The last two involve respectively 16,411 and 30,672 binary inclusion
pairs. Thus the working comparison theorems hold for binary patterns on those
shapes at q=3 as well as q>=4. On 2x2 and 2x3 every ternary inclusion pair is
reachable too. These are exact, finite pattern-class certificates with
all-volume consequences, not a proof for arbitrary rectangles.

The possible structural lemma is: **for binary patterns on a rectangle,
agreement inclusion is generated by equal-G relabelings and mergers using
at most three letters.** It would settle the q=3 binary-pattern case. Three
letters cannot connect all inclusions between genuinely ternary patterns,
since a proper merger has an image using at most two symbols.

Even lattice convexity is insufficient for this proposed structural lemma.
Take S={(0,0),(0,1),(1,0),(1,1),(1,2),(1,3)} and patterns
u=(10;1100), v=(01;1010), where the unequal row lengths specify S.
Then G(u)={+-(1,1)} and G(v)={+-(0,2),+-(1,-1),+-(1,1),+-(1,3)}.
Every representative of G(u) has the form (AB;CABD), with B!=C, C!=D,
A!=D; these inequalities are also sufficient. Its four forced components
are A,B,C,D. The target v joins B with C and A with D. To refine that target
while retaining G(u) requires four distinct labels: both target fibers must
split. Thus there is no direct ternary merger to its class. The only other
strict intermediate agreement classes inside this interval add (1,-1) alone
or (1,3) alone, and each requires three symbols. They cannot be the endpoint
of a proper ternary merger. Hence no directed merger path exists. This is a
proof-method obstruction, NOT a counterexample to entropy or convex order.
The same obstruction occurs on the five-point one-dimensional support
{0,1,2,5,7}, for u=(1,1,0,0,0), v=(0,1,1,1,0).

## The lower extreme when an overlap-free pattern exists

**Proposition (Shearer bound applied to occurrences).** If there is u on S with
G(u) empty, then for every other w on S, every finite anchor set V, and every
real beta, Z_u(beta) <= Z_w(beta). The same order holds for avoidance counts
and hence for entropy. This does not assert strictness or full convex order.

Here is a self-contained special-case proof. Make a graph H on V by joining
anchors whose supports overlap. For every B subset V put
Q(B)=sum_{I independent in H[B]}(-p)^{|I|}, p=q^{-s}. For u, pairwise
exclusivity on edges and independence on nonedges give
Q(B)=Pr(N_{u,B}=0)>0. Positivity follows also from a constant configuration
avoiding u. For any other w let F(B)=Pr(N_{w,B}=0).
Inductively prove F(B)/F(B\{i}) >= Q(B)/Q(B\{i}). Indeed, independence
outside the closed neighborhood N[i] gives
F(B) >= F(B\{i}) - p F(B\N[i]). The induction applied while adding the
neighbors of i gives F(B\N[i])/F(B\{i}) <= Q(B\N[i])/Q(B\{i}).
Use the independent-set recursion Q(B)=Q(B\{i})-p Q(B\N[i]).
All denominators are positive (constant configurations avoid w as well).
Multiplying the ratio inequalities gives F(V)>=Q(V).

For beta<=0, independently retain each occurrence with probability
lambda=1-exp(beta). The probability no retained occurrence remains is
E exp(beta N). The same proof applies with p replaced by lambda p; its Q(B)
is the corresponding strictly positive expectation for u. For beta>=0 use
the marked-occurrence expansion and G(u) subset G(w). QED.

This is a direct application of Shearer's extremal avoidance mechanism, not a
claim to a new general probability theorem. A primary reference spelling out
the extremal measure and bound is Achlioptas--Zampetakis (ICALP 2021), Section
1.3, https://doi.org/10.4230/LIPIcs.ICALP.2021.8. The overlap-free realizability
hypothesis is necessary to identify this bound with a pattern in the class.

**Proposition (working proof, strict lower extreme).** If G(u) is empty and
G(w) is nonempty, then P_u(beta)<P_w(beta) at every finite beta!=0. Moreover,
if some alphabet symbol is absent from w, then h(X_u)<h(X_w). The finite-coupling
statement holds for arbitrary finite supports; so does the spare-symbol entropy
statement. This particular hard-endpoint argument needs a spare symbol; the maximal-pair
argument below removes that restriction for every rectangle and for arbitrary
supports of size at least 19.

For clarity, refine the preceding induction. For retained events A_i of marginal
p=lambda q^{-s}, write F(B) for avoidance and Q(B) for the overlap-free comparison.
There is an exact correction term
E_i(B)=Pr(A_i, no event in B\N[i], at least one event in B intersect N(i)),
so F(B)=F(B\{i})-p F(B\N[i])+E_i(B).
The proved ratio comparison then gives
F(B)/F(B\{i}) >= Q(B)/Q(B\{i}) + E_i(B)/F(B\{i}).

Choose a compatible displacement d for w, with d lexicographically positive,
and j=i+d in B. Let U be the union of the two supports and J all occurrence
anchors whose supports meet U. Modify letters on U to realize the two w-copies,
retain i,j, and suppress all other retention flags in J. Any exterior admitting
an assignment avoiding B\{i} also admits this assignment contributing to E_i.
The product probability of the prescribed local assignment is at least
c=q^{-2s} lambda^2 (1-lambda)^{2s^2}>0 for 0<lambda<1,
since |U|<=2s and |J|<=2s^2. Conditional old avoidance probability is at most
one, so E_i(B)>=c F(B\{i}). Since Q(B)/Q(B\{i})<=1, this yields a
multiplicative ratio improvement by at least 1+c. Eliminate anchors in
lexicographic order; the fixed forward neighbor j is still present at all but
a boundary fraction of steps. Thus P_w-P_u>=log(1+c)>0 at every beta<0.
The positive-coupling strict comparison is already in the shared manuscript.

At the hard endpoint, suppose star is absent from w. Prescribe the compatible
pair on U and fill the collar [U+(S-S)]\U with star. Any other copy lying
inside U must intersect the first support: a disjoint copy would have to fit
inside the second support minus its nonempty overlap with the first, a set
of cardinality smaller than s. Copies crossing the collar are impossible.
Hence this prescription again contributes to E_i(B) for every extendible
exterior. It has probability q^{-|U+(S-S)|}, giving the same uniform ratio
argument and a strictly positive entropy gap. QED.

The retention flags are auxiliary independent Bernoulli variables, used only
to represent exp(beta N) as a no-retained-occurrence probability. They are not
new lattice symbols. At lambda=1 the suppression step disappears, which is
the exact obstruction to using this proof without a spare symbol.

## A quantitative K-order theorem for weak overlaps

**Theorem (working proof).** This result allows arbitrary finite supports and
includes binary patterns. Write s=|S|, D=|S-S|-1 and p=q^{-s}. Let
L_w=|union_{d in G(w)} [S intersect (S+d)]| count the sites that can be covered
by another compatible copy. For a pair u,v put
L=max(L_u,L_v), g=max(1,|G(u)|,|G(v)|), and

R = min{q^{s-L}/(2e^2 g), q^s/[4e(D+1)]}.

If R>=2, then, uniformly for -1<=lambda<=0,

P_w(log(1+lambda))-log(q)
= p lambda + [(K_w-D-1)/2] p^2 lambda^2 + E_w(lambda),

with |E_w(lambda)| <= 4ep R^{-2}|lambda|^3. At lambda=-1 the left-hand
side means h(X_w)-log(q). Consequently, if

K_v-K_u > 16e q^s/R^2,

then P_u(beta)<P_v(beta) for every beta<0 and h(X_u)<h(X_v). The entropy
gap is at least p^2(K_v-K_u)/2 - 8ep/R^2. No spare symbol is required.

Proof. Expand E(1+lambda)^N in marked occurrences and group each marked set
into its connected components in the support-overlap graph. A polymer C is
a nonempty compatible connected set of anchors, with activity
z(C)=lambda^{|C|}q^{-|S+C|}; two polymers are incompatible when their supports
meet. Every newly inserted compatible copy overlaps the union of previous
copies only at its L_w specified sites. Thus
|S+C| >= s+(|C|-1)(s-L).
There are at most (eg)^{k-1} connected compatible k-vertex anchor sets
containing a fixed anchor: use the graph with edge displacements in G(w),
and count rooted embedded trees. Consequently

sum_{C contains x} |z(C)|e^{|C|}
<= e|lambda|p/[1-e^2 g |lambda|q^{-s+L}].

For |lambda|<=R the denominator is at least 1/2. A polymer with k anchors
has at most k(D+1) anchor-neighborhood sites where an incompatible polymer
can start. The Kotecky--Preiss sum is therefore at most
2k(D+1)eRp <= k/2.
Apply the standard polymer convergence criterion with a(C)=|C|. It gives an analytic
free-box limiting function F_w(lambda), defined using the partition polynomials
and the logarithm branch at zero, on this disk. For real lambda>-1 it equals
P_w(log(1+lambda))-log(q), and |F_w(lambda)|<=2eRp throughout the disk. The single-polymer
sum per unit volume is bounded by the displayed anchored sum; the usual
pinned-cluster bound supplies the factor e^{|C|}. Uniform convergence
justifies the limiting analytic expansion, including lambda=-1.

The first coefficient is p. The second is
p^2(K_w-D-1)/2: subtract p^2/2 for a repeated anchor, subtract Dp^2/2 for
overlapping pairs, and add p^2 K_w/2 for compatible overlapping pairs.
Cauchy's bound on the remaining coefficients gives, for |lambda|<=1,
|E_w(lambda)| <= 2eRp (|lambda|/R)^3/(1-1/R)
<=4ep R^{-2}|lambda|^3.
Subtract the two expansions and use |lambda|^3<=lambda^2. QED.

The convergence theorem and pinned-cluster bound are standard inputs, not new
claims: Kotecky--Preiss, CMP 103 (1986), 491--498,
https://doi.org/10.1007/BF01211762; see also Fernandez--Procacci,
https://arxiv.org/abs/math-ph/0605041, equations (2.12)--(2.15), checked in the
primary paper. The constants are deliberately loose.

**Concrete binary consequence.** Put
B_w=(1/2)sum_{d in G(w)}|S intersect (S+d)|. On EVERY two-dimensional
rectangle of area s>=24, B_u,B_v<=2 and K_u<K_v imply the strict comparisons
above. Indeed L<=4, g<=4, D+1<=4s, and the comparison threshold is at most
max{1024e^5 2^{-s+8},4096e^3 s^2 2^{-s}}. At s=24 these are e^5/64 and
9e^3/64, both less than 4 (use e<11/4), and both decrease thereafter. A
positive K difference is at least 4. This estimate also implies R>=2.
In particular, it gives strict comparison against the overlap-free minimum
for such w with G(w) nonempty, even when w uses both symbols. Strongly
periodic patterns and K ties remain outside this argument.

## Strict entropy minimality: a maximal pair and the conditional local lemma

**Separate lower-extreme question, unproved.** Let D be a finite dependency
graph and p=(p_i) a marginal vector for which the Shearer law mu exists as a
probability measure. In occurrence-indicator convention its mixed moments are
mu(X_i=1 for i in I)=product_{i in I} p_i for D-independent I, and zero
otherwise. Suppose another Bernoulli field Y has the same marginals and each
Y_i is independent of the ENTIRE field outside its closed D-neighborhood.
Must mu <=_sm Law(Y)? The weaker question asks only for convex order of every
nonnegative weighted total. An affirmative answer would settle the remaining
finite convex-order question at the overlap-free extreme, because its field
is the Shearer law for the physical-overlap dependency graph.

**Elementary count case (working observation).** If the independence number
of D is at most two, the UNWEIGHTED count comparison already follows from
the usual Shearer avoidance bound, for arbitrary feasible marginals p_i.
Indeed N_mu is supported on {0,1,2}. At threshold one,
E(N-1)_+=E N-1+P(N=0), so Shearer's avoidance inequality and equal means
give the needed hinge comparison. At every integer threshold at least two,
the source hinge is zero; the remaining thresholds follow from equal means
and piecewise linear interpolation. Thus N_mu<=_cx N_Y. The same applies
to any induced subsystem of independence number at most two. This argument
does not establish arbitrary weighted-count or field comparisons.

The general field and weighted-count claims remain auxiliary conjectures,
beyond the known avoidance comparison. They use the weak single-event/nonneighbor independence
condition; independence of two separated multi-event families is stronger.
See Temmel's definitions (4)--(6), https://arxiv.org/abs/1105.1683, after
reversing his zero/one convention. The supermodular local lemma of
Achlioptas--Zampetakis, Theorem 5, https://arxiv.org/abs/2111.08837, bounds a
supermodular FUNCTION satisfying a multiplicative factorization condition.
It does not directly compare expectations of arbitrary supermodular TESTS
under two field laws. That distinction is essential here.

A concrete sufficient-and-necessary finite dual target is: for each
supermodular f, find functions h_i on subsets of V minus N_D[i] such that,
for every A subset V,

f(A) >= E_mu f + sum_i (1_{i in A}-p_i) h_i(A minus N_D[i]).

Equality must hold on every positive-mass atom of mu. Averaging cancels each
summand by the stated independence. Conversely finite linear-programming
duality gives such h_i if mu minimizes the expectation under these marginal
and independence constraints. This reformulation is exact, but construction
of the h_i for general D and f is OPEN. Preliminary floating-point searches
have found no violation; they are exploratory, not exact universal certificates.

**Partial Shearer field theorem (working proof).** The preceding field-minimum
assertion holds when the COMPLEMENT of D is chordal. Its class of valid
dependency graphs is also closed under disjoint unions and graph joins.
Consequently it holds for every graph constructed from single vertices by
these two operations, and for unions/joins of graphs with chordal complements.

For the chordal-complement case, weak independence forces the mixed moment
product_{i in I} p_i on every D-independent set I. Thus E_mu f equals the
expectation under Y of the clique-Mobius sum T_{complement D} f on its active
set. The chordal inequality proved above bounds that sum by f of the active
set, pointwise. This proves the assertion.

For a disjoint union D=D_1 disjoint-union D_2, condition on the whole second
field. Each event in the first part still has marginal p_i and is independent
of its first-part nonneighbors, because originally it was independent of
their JOINT field with the whole second part. Apply the assumed first-part
comparison to each conditional supermodular test. It replaces that part by
an independent Shearer field. Apply the second-part comparison afterward.
The product of the two Shearer laws is the Shearer law on the union.

For a join D=D_1 join D_2, every cross pair is a dependency edge, so the
Shearer field has events active in at most one part. Supermodularity gives
f(A union B)>=f(A)+f(B)-f(empty) for the disjoint parts. Average under Y and
apply the assumed comparison to both part-marginals. The resulting bound is
exactly E_mu f, because the displayed inequality is equality on the support
of mu. Singleton graphs start the induction. QED. This uses only the stated
weak independence conditions and Shearer feasibility, not cylinder geometry.
The connected five-vertex path and five-cycle lie outside these arguments.

**Five-vertex Shearer theorem (working proof with exact polynomial
certificates).** On EVERY dependency graph with at most five vertices, at
EVERY feasible common marginal p, the Shearer field is minimal in
supermodular order among ALL fields satisfying the stated weak independence.
Unequal marginals remain outside the polynomial certificates.

The analytic reductions leave just three graph types: P_5, C_5, and a triangle
with a pendant two-edge path, whose edges we label 01,12,23,24,34. To see
this, suppose D and its complement are connected and the complement is
nonchordal. An induced five-cycle uses every vertex. Otherwise the complement
contains an induced four-cycle, with at most one vertex left. If that vertex
has zero or four neighbors on the cycle, the graph or its complement is
disconnected. With two opposite neighbors the complement D is disconnected,
as it is with three neighbors. The remaining one-neighbor and two-adjacent-
neighbor cases give the triangle-with-path and P_5. The four-vertex exception
is the disjoint union of two edges. This proves the reduction; the checker
also audits all 1,099 labeled graphs on one through five vertices.

Here is the reduction used by `shearer_order_checks.py`. Write delta=nu-mu.
For every D-nonindependent active set A, put r_A=delta(A)=nu(A)>=0, since
mu(A)=0. All independent-set mixed moments of delta vanish. Therefore its
remaining atom weights are uniquely recovered, in decreasing cardinality, by

delta(I)=-sum_{A strictly containing I} delta(A), for independent I.

The weak independence conditions for vertex i and a nonindependent set
T contained in V minus N_D[i] are exactly

sum_{A containing T} (4*1_{i in A}-1) r_A=0.

Here 4 is replaced by q=1/p for the general parameter. Conditions for
independent T already follow from the fixed independent-set
moments. Thus every admissible law difference comes from the nonnegative
kernel cone Br=0. There are 21 variables and 5 equations for C_5, and 19
variables and 8 equations for P_5. The independent atom nonnegativity
conditions restrict this cone further and need not be used for this sufficient
comparison.

The checker enumerates its rays using exact rational arithmetic: start with
the coordinate rays and intersect successively with the equation hyperplanes.
Keep zero rays and combine every pair of opposite-sign rays; retain a candidate
with support S precisely when the restricted equation matrix has rank |S|-1.
This is complete because a new ray of a hyperplane section comes from a
two-dimensional face of the preceding cone, and this nonnegative kernel cone
has an extreme ray on S exactly when its restricted kernel is one-dimensional.
All candidates are normalized to primitive integer vectors.

Every resulting full delta has a saved nonnegative rational decomposition
into elementary supermodular differences. C_5 has 502 extreme rays and 5,991
certificate terms; P_5 has 32 rays and 306 terms. The standard-library replay
regenerates all rays and verifies all 6,297 terms exactly against
`shearer-order-evidence.json`. Hence every conic combination is positive on
every supermodular test. The Shearer probabilities themselves are also
checked exactly.

The parameter extension is also exact. Write t=q-3. At every stage of the
cone construction, the checker supplies a polynomial vector on each ray
support. All its nonzero coordinate polynomials have nonnegative rational
coefficients in t, and its next-equation evaluation is identically zero or
has coefficients of one strict sign. Thus all ray supports and all cut signs
remain the same for every t>0. Completeness does not assume constant matrix
rank: a cone generated by rays, cut by one homogeneous equation, is generated
by the zero rays and pairwise positive/negative intersections. For a
nonnegative kernel cone the extreme supports are precisely the inclusion-
minimal supports among those generators. The checker verifies this support
recurrence explicitly at every stage.

Graph symmetries reduce the final polynomial rays to 65 orbits for C_5 and
20 for P_5. The triangle-with-path has 35 rays in 27 orbits and is handled by
the same construction. Reconstruct delta coefficient by coefficient; every coefficient
has a nonnegative rational elementary supermodular decomposition. There are
197/52/59 such coefficient certificates, totaling 1,699/370/345 terms. Hence each
ray is positive on every supermodular test for all q>3. This covers every
feasible common marginal: at p=1/3 the formal all-zero masses are -1/9 for
C_5 and the triangle-with-path, and -1/27 for P_5, while a feasible Shearer law with larger p could be
independently thinned to give a feasible law at 1/3. Thus feasibility forces
p<1/3; p=0 is trivial. The replay uses only rational polynomial arithmetic
and the standard library. `--solve-symbolic` regenerates the polynomial
certificates using SymPy and SciPy. This settles every graph through five
vertices after the analytic reductions, not arbitrary dependency graphs.

**Theorem (working proof; small rectangular cases use exact rational
certificates).** If an overlap-free pattern exists on a two-dimensional
rectangle S, then its agreement class is exactly the class of entropy
minimizers, for every q>=2. More generally, this strict minimum statement
holds on ANY finite support of size at least 19 whenever an overlap-free
pattern exists there.

The key sufficient condition for an arbitrary support is the following. Put
s=|S|, D=|S-S|-1, and let r be the MAXIMUM overlap size among the compatible
displacements of a pattern w with G(w) nonempty. If

q^{-max(r,s-r)} <= 1/[e(D+1)],

then h(X_w)>h(X_u) for an overlap-free u, with gap at least
log(1+e^{-1}q^{-(2s-r)}).

Proof of the sufficient condition. We use the elementary conditional
local-lemma estimate: if events have a dependency graph of maximum degree D
and probabilities at most x(1-x)^D, every conditional occurrence probability,
given avoidance of any other subset, is at most x. Induct on the size of that
subset: split it into neighbors and nonneighbors, use independence for the
latter, and bound the denominator by (1-x) for each neighbor using the
induction hypothesis. This also proves positivity of every avoidance
probability. Set x=1/(D+1), so x(1-x)^D>=1/[e(D+1)].

For the Shearer defect E_i(B) above, choose j=i+d in B with d forward and
having maximal overlap r. Condition on realizing both w-copies, fixing the
letters on U=(S+i) union (S+j), with probability q^{-(2s-r)}.
Let C=B\N[i]. Its supports miss S+i, so each meets U in at most s-r sites
and has at least r unfixed sites. If such a copy meets S+j, its conditional
probability is zero unless it is compatible with that copy. In the latter
case its overlap with S+j is at most r BY MAXIMALITY, so it also has at
least s-r unfixed sites. Thus every conditional occurrence probability at
C is at most q^{-max(r,s-r)}. The conditioned letters outside U remain a
product measure and the dependency degree is still at most D.

Let J be the anchors in C whose supports meet U; they overlap S+j, so
|J|<=D. Add their avoidance conditions last. The conditional local lemma gives
Pr(no C | copies at i,j)>=(1-x)^D Pr(no C\J).
The last probability is unconditional because C\J is disjoint from U.
Since F_w(B\{i})<=Pr(no C\J),
E_i(B)>=q^{-(2s-r)}(1-x)^D F_w(B\{i}).
Insert this uniform defect into the Shearer ratio inequality and eliminate
anchors lexicographically. All but a boundary fraction have their forward
neighbor j still available. The gap is at least
log(1+q^{-(2s-r)}(1-x)^D), proving the sufficient condition.

Since max(r,s-r)>=ceil(s/2) and D+1<=s(s-1)+1, the sufficient condition holds
for s=19 and every s>=21: use e<11/4, check s=19,21,22, and then increment
by two. At s=20 the same estimate works unless r=10. In that case the r
ordered pairs realizing d, and their reversals, give
D+1<=s(s-1)+1-2(r-1)=363; now 2^{10}>(11/4)363. Thus every support with
s>=19 is covered. For a 2D rectangle, D+1<4s, and
2^{ceil(s/2)}>=11s for every integer s>=16 (check 16,17 and again increment
by two). This proves the rectangular statement for area at least 16.

The remaining nondegenerate rectangles have area below 16. Binary 2x2, 2x3,
3x3 have no empty-G pattern; they are not assertions of this theorem.
For all other cases, the maximal-overlap sufficient condition leaves only
the displacements below, up to reflections and transpose. The q=3 rows also
cover every larger q, by monotonicity of the probability bounds.

| Shape | q | Remaining representative pair displacements d |
|---|---:|---|
| 2x4 | 2 | (1,0), (1,1), (0,2) |
| 2x5 | 2 | (0,3), (1,1), (1,0), (0,2) |
| 2x6 | 2 | (0,3), (1,0) |
| 3x4 | 2 | (0,2), (1,1) |
| 2x2 | 3 | (0,1) |
| 2x3 | 3 | (1,0) |

Here use the inhomogeneous local lemma directly on the infinite set C of
anchors whose supports miss the first fixed copy. For k in C set
p'_k=q^{-| (S+k)\(S+d) |}, an upper bound on its conditional occurrence
probability, and set x_k=2p'_k. Only the finite set J of anchors meeting
S+d differs from the constant value 2q^{-s}. The required inequalities are

product_{l in C: (S+l) meets (S+k), l!=k}(1-x_l) >= 1/2.

`minimum_certificates.py` checks these inequalities with exact fractions on
J and its neighbor set. Outside that finite set the uniform lower bound
(1-2q^{-s})^D suffices. It also checks x_k<1 and the complete small-shape
case coverage. All THIRTEEN certificates pass; the smallest neighbor-product
bound is greater than 0.58. The exact fractions and coverage counts are in
`minimum-certificates.json`; displayed decimals are not used in the tests.
The inhomogeneous conditional local lemma follows from the same induction,
with a factor 1-x_l for each neighbor. The defect estimate now has the
positive factor product_{k in J}(1-x_k), uniformly in B. This finishes every
small rectangular case. Degenerate one-dimensional rectangles are already
covered by the strict full-shift word comparison in the shared manuscript.
QED.

This establishes the lower EXTREME, not the full agreement-inclusion entropy
conjecture. The weak-overlap K theorem above also compares nonminimal patterns
and remains a separate result. The general strict minimum on supports of fewer
than 19 sites with holes is not settled by this argument.

## When the lower extreme is realizable in two dimensions

**Proposition (working proof).** A binary m by n rectangle with m,n>=2 admits
a pattern with empty G if and only if max(m,n)>=4. Every such rectangle admits
an overlap-free pattern over three symbols. One-dimensional rectangles admit
an overlap-free binary word at every length.

For the binary construction, transpose if necessary so n>=4. Set
A=0^{n-2}11 and B=01 0^{n-3}1, and use m-1 rows A followed by one row B.
A has no proper border. No prefix of A equals a suffix of B of the same positive
length, and no suffix of A equals a prefix of B (including full length).
For lengths <=n-2, the former fails at the last symbol; at n-1 it fails at the
first symbol; full words differ. For the latter, length 1 fails, length 2
compares 11 with 01, lengths 3 through n-1 fail at the last symbol, and the full
words differ. Horizontal overlaps therefore fail by A, and every nonzero
vertical displacement has a row comparing a translate of B with A, so fails.
For three symbols, take a zero top row and fill all remaining sites with 1
except the bottom right site, labeled 2. A vertical displacement compares a
zero top-row symbol with a nonzero symbol; a purely horizontal displacement
conflicts at the unique 2.

For nonexistence it suffices, up to transpose, to inspect 2x2, 2x3, 3x3.
For 2x2, write rows ab and cd. Avoiding both diagonal corner agreements forces
d=1-a and c=1-b. If a=b, the columns agree; otherwise the rows agree.
For 2x3, write rows abc and def. Corner conflicts force f=1-a,d=1-c.
If a=c, the displacement by two columns agrees. If a!=c, then d=a,f=c.
Avoiding vertical agreement requires b!=e, whereas avoiding both diagonal
two-site overlaps forces e=b (use b=a or b=c), a contradiction.
For 3x3 write rows abc, def, ghi. Corner conflicts give i=1-a,g=1-c.
If a!=c, the same top/bottom-row argument forces b=h and hence agreement at
vertical displacement two. Otherwise complement if needed so a=c=0,g=i=1.
Avoiding horizontal displacement two requires d!=f; whichever of d,f is 0,
one of the displacements (1,2),(1,-2) agrees. The middle symbol e is irrelevant.
This proves the claimed classification. QED.

Thus the overlap-free theorem identifies the minimum outside binary 2x2,
2x3, 3x3 and their transposes. In those shapes several incomparable minimal
agreement classes exist, so inclusion alone does not identify the global
minimum. The following arguments settle all three exceptions.

## The exceptional binary 3x3 minimum

**Theorem (working proof, with exact rational constants).** Among binary 3x3
patterns, the entropy minimizers are exactly those whose agreement set is
{+-(2,2)} or {+-(2,-2)}. These two classes are related by reflection. They
also strictly minimize limiting pressure at every beta<0. For example
(110;110;000) represents the second class. The other binary 3x3 patterns
have entropy larger by more than 0.0000033076. This is an analytic planar
comparison, not a strip extrapolation or an appeal to the open inclusion
conjecture.

Put p=2^{-9}, theta=1-exp(beta), and let
Q_theta(V)=sum_{I independent in the support-overlap graph}(-theta p)^{|I|}.
Write P_0(theta)=log(2)+lim |V|^{-1}log Q_theta(V) along free boxes.
This formal Shearer comparison is positive and its pressure exists by the
polymer bound below. It is NOT the pressure of an overlap-free binary 3x3
pattern, since such a pattern does not exist. The Shearer ratio induction and
its quantitative defect apply whenever all these Q's are positive; their
earlier proof used realizability only to ensure this positivity.

First let w have a single corner-agreement pair. Every compatible connected
anchor set is a consecutive chain in that corner direction. A chain of k
copies has union size 8k+1, since adjacent copies share one site and
nonadjacent copies are disjoint. Its polymer activity is
z_k=(-theta)^k p rho^{k-1}, where rho=1/256. Exactly k such chains contain
any given anchor. With a=1/16, the Kotecky--Preiss sum for a k-anchor polymer
is bounded by

25k p exp(a)/(1-rho exp(a))^2 <= (3000/57121) k < ak.

Here exp(a)<=16/15, and 25 is the closed overlap-neighborhood size. This
also proves convergence and positivity for the monomer gas Q_theta and for
every interpolation multiplying all nonmonomer activities by t in [0,1].
For a fixed polymer C, differentiating the finite partition function with
respect to its activity gives the ratio R_C of the partition function
restricted to polymers compatible with C to the original partition function.
Both are positive. The pinned-cluster bound gives
|log R_C|<=a|C|: bound the clusters removed by the restriction by the sum of
the pinned bounds over polymers incompatible with C. Thus R_C<=exp(a|C|).
Integrate the interpolation derivative from t=0 to t=1. Odd-chain terms
are negative, so discarding them and bounding the even terms gives

P_w(beta)-P_0(theta)
<= theta^2 p rho exp(2a)/(1-rho^2 exp(2a))
<= theta^2 U,                    U=1/115198.

The same inequality holds at theta=1 for entropy. All steps are uniform in
free boxes, with their usual boundary corrections.

For any other pattern let r be its maximal compatible overlap size. The
possible values on 3x3 are 1,2,3,4,6. The conditional local-lemma/Shearer
argument above gives

P_w(beta)-P_0(theta) >= log(1+c theta^2)
                       >= theta^2 c/(1+c),

with the following explicit positive constants c. Retention only decreases
the conditional event probabilities, so their hard-endpoint bounds work for
every theta in [0,1].

* If r=2, use the uniform conditional local lemma with x=1/100:
  x(1-x)^24>=1/128. Then c=2^{-16}(99/100)^24.
* If r=3,4,6, the inhomogeneous certificates for displacements (0,2),
  (1,1),(0,1), respectively, give
  c=2^{-(18-r)} product_{k in J}(1-x_k), with x_k=2p'_k as above.
  Their c values exceed 0.00001801, 0.00005253, 0.00022561, respectively;
  the actual comparisons use exact fractions.
* If r=1 and w is not in a single-corner class, both corner pairs are
  compatible; empty G is impossible. At an interior anchor i, take the two
  forward neighbors j_+=(2,2)+i and j_-=(2,-2)+i. Their supports are disjoint.
  After fixing the pair i,j_+, every remaining event in C=B\N[i], and also
  the event at j_-, has conditional probability at most 1/256. The latter
  overlaps the fixed union at just one site. Use x=1/200, for which
  x(1-x)^24>=1/256, and require avoidance of C and of j_-. At most 25 of
  these events meet the fixed union. Repeat with + and - exchanged. The
  resulting two witness events are disjoint, so their defects add:
  c=2*2^{-17}(199/200)^25.

In every case c/(1+c)>U. `minimum_certificates.py` checks these rational
inequalities, the chain convergence bound, and the three additional local
certificates. The smallest lower bound on c/(1+c)-U is greater than
0.0000033076, attained by the r=2 estimate. Consequently every other pattern
has a pressure gap at least theta^2 times this constant, and an entropy gap
at least this constant. Equality of the two single-corner pressures follows
from reflection and the equal-G occurrence-law identity. QED.

**Theorem (working proof, complete finite-class certificate).** For binary
2x3 patterns, the classes G={+-(1,2)} and G={+-(1,-2)} are exactly the
entropy minimizers. They also strictly minimize limiting pressure at every
beta<0. Every other class has a pressure gap greater than
0.00004430 (1-exp(beta))^2 and an entropy gap greater than 0.00004430.
The two classes are reflections; (110;000) represents the second one.

Here p=1/64, rho=1/32. The single-corner polymers again form chains, now
with activity (-theta)^k p rho^{k-1}. Their closed overlap neighborhoods
have size 15+9(k-1). Use a(C)=alpha |N[C]| with alpha=1/40. The required
pinned sum per possible anchor is

p exp(15 alpha)/(1-rho exp(9 alpha))^2 < 0.024624 < alpha.

The same interpolation proof therefore gives the upper bound theta^2 U_23
above the formal Shearer pressure, where

U_23 = p rho exp(24 alpha)/(1-rho^2 exp(18 alpha)) < 0.000891072.

The code uses a rational upper bound for each exponential, obtained from its
Taylor polynomial through degree ten and a geometric tail. Neither this comparison
nor the convergence test uses floating-point arithmetic.

Two refinements of the previous defect bound complete all other classes.
Both are useful independently of this small-shape calculation.

First, a conditional local-lemma certificate after fixing ONE copy gives a
uniform lower bound b on its conditional occurrence probability given avoidance
of any other finite set. To see this, condition on that copy, apply the
certificate to the neighboring events, and use independence for all remaining
events. If x_k are its certificate values, one can take
b=p product_k(1-x_k). For retained events the corresponding bound is theta b,
since scaling every x_k by theta preserves the local-lemma inequalities.
If no such certificate is used, put b=0.

Eliminate anchors in lexicographic order. At an interior step, B contains
exactly seven overlap neighbors of the current anchor i. With C=B\N[i],
the conditional rate bound gives

F(B\{i}) <= (1-theta b)^7 Pr(no event in C).

Second, choose an ordered list of compatible forward neighbors j_1,...,j_m.
For each l prescribe the pair i,j_l, require avoidance of C, and also require
that j_1,...,j_{l-1} do not occur. These witness events are disjoint. After
fixing the pair, an event at k has conditional probability zero if its requests
conflict with either fixed copy; otherwise its probability before retention
is 2^{-| (S+k)\((S+i) union (S+j_l)) |}. The certificate domain is the
lexicographically forward anchors disjoint from S+i, plus the earlier j's.
Only finitely many probabilities differ from p. Exact local-lemma checks on
this set and its neighbors certify values x_k; outside it the constant value
Mp is checked, where M is the multiplier chosen by the certificate.

Let J_l be the certificate events meeting the fixed union, let a_l be the
probability of the fixed unretained pair, and put

c_l(theta) = a_l product_{k in J_l}(1-theta x_k)/(1-theta b)^7.

The two preceding inequalities give E_i(B)/F(B\{i}) >= theta^2 sum_l c_l(theta).
All nonzero x_k in these certificates satisfy x_k>=b. Thus c_l(theta) is
bounded below on [0,1] by the smaller of its two endpoint values. For b=0
this is monotonicity. For b>0, the sign of its logarithmic derivative is the
sign of
7-sum_k (x_k/b)(1-theta b)/(1-theta x_k), a decreasing function of theta.
Any interior extremum is therefore a maximum. Consequently the uniform
defect constant is

c = sum_l min{a_l, a_l product_{k in J_l}(1-x_k)/(1-b)^7}.

As before the pressure gap above the formal Shearer comparison is at least
theta^2 c/(1+c). The script exhausts all 64 binary words, checks their 18
agreement classes against the independent agreement routine, excludes the
two single-corner classes, and supplies a list of at most three neighbors
and exact certificates for EACH of the other 16 classes. Every resulting
c/(1+c) exceeds the rational upper bound on U_23. The least certified gap
is greater than 0.00004430. The witnesses, conditional-rate certificates,
pair certificates, and exact fractions are in `minimum-certificates.json`.
This proves the stated planar comparison without assuming agreement inclusion
at q=2. QED.

## The binary 2x2 minimum: explicit counting bounds

**Theorem (working proof, exact certificate).** The binary 2x2 entropy
minimizers are exactly G={+-(1,1)} and G={+-(1,-1)}. Up to lattice symmetry
and alphabet permutation there are just four classes to compare: a corner
pattern (10;00), an axis pattern (11;00), a checkerboard (10;01), and the
constant pattern. The constant is already strictly maximal. The exact bounds
in `binary_two_by_two_certificate.py` are

| Class | Certified planar entropy bound |
|---|---:|
| Single corner | < 0.623432202 |
| Axis | > 0.623749256 |
| Checkerboard | > 0.624881660 |

Thus the corner minimum is strict. This also settles every strict
agreement-inclusion entropy comparison on binary 2x2. This argument does NOT
order the axis and checkerboard classes and does not compare all finite
negative couplings.

For the upper bound, remove constraints crossing the boundaries of vertical
strips of width W=16. This gives h_corner <= W^{-1} log rho(T_W), where
T_W is the free-boundary row transfer. For a binary row a let A(a) be its
starts of 10 and B(a) its starts of 00. Then
T_W(a,b)=1{A(a) intersect B(b) is empty}. Its rows depend only on A(a).
Group equal A masks and evaluate transfer actions by exact subset sums over
B masks. Ten normalized INTEGER iterations produce a positive vector v with

T_W v <= (87985395/4096) v.

The code checks this inequality exactly at all 1,597 matching-mask states.
It also checks the subset-sum transfer against direct row-pair enumeration
for widths 2 through 6. Rational logarithm bounds give the stated upper bound;
no numerical Perron approximation or planar extrapolation is used.

For both lower bounds, assign independent fair bits to the even parity sites.
Each plaquette then either imposes no condition on its two odd sites or
forbids one assignment on a diagonal edge between them. Build the resulting
graph on the odd sites. Count allowed assignments independently on its small
path components, and fill the remaining components by a fixed allowed
assignment. Summing over the even bits and applying Jensen gives

h >= (1/2)log(2) + expected log(number of chosen odd assignments) per site.

Only components of a fixed bounded size are used, so the omitted boundary
contribution is o(volume). The construction is a direct lower bound on
numbers of avoiding configurations; it does not require a claim about the
measure of maximal entropy.

For the checkerboard (10;01), an edge in one diagonal direction is activated
by two even zeros and forbids two odd ones; in the other direction it is
activated by two even ones and forbids two odd zeros. Edges of different
types cannot meet. Each component is therefore a straight path of a single
type. Along either diagonal, successive activation events use disjoint even
bits, hence are iid with probability 1/4. The probability an odd vertex is
isolated is 1/8. A run of k active edges with inactive ends has probability
(3/4)^2 4^{-k}, and its k+1 vertices have F_{k+3} valid colorings, where
F_0=0,F_1=1 are the Fibonacci numbers. There are two directions and half as
many odd sites as lattice sites. Keeping runs through length m gives

h_checker >= (9/16)log(2) + (9/16)sum_{k=1}^m 4^{-k} log F_{k+3}.

The displayed lower bound uses only m=4. Larger components have a constant
valid filling, chosen according to their type.

For the axis pattern (11;00), an activated edge requires that its two even
sites, ordered from top to bottom, have values 1,0; its odd endpoints must
then satisfy top<=bottom. At any odd vertex, the two left edges cannot both
be active, nor can the two right edges. Components are paths moving
monotonically horizontally. An odd vertex is isolated with probability 3/8.
For a specified first edge to be active and its initial vertex an endpoint
has probability 1/16. On exploring the path, a straight continuation fixes
two fresh even bits (probability 1/4), a turn fixes one (probability 1/2),
and termination has probability 1/4. Thus a specified k-edge path with t
turns has probability (1/64) 4^{-(k-1)} 2^t of being a complete component.

For its vertical step sequence epsilon in {+,-}^k, the number Z(epsilon)
of odd colorings is computed from (a,b)=(1,1) by replacing (a,b) with
(a,a+b) on a downward step and with (a+b,b) on an upward step; the final
count is a+b. Fix the first step downward, using symmetry for the other
three initial directions. Counting components at their two endpoints, and
including the density 1/2 of odd sites, yields the finite lower bound

h_axis >= (11/16)log(2)
          + (1/64)sum_{k=1}^m 4^{-(k-1)}
                    sum_{epsilon_1=+} 2^{t(epsilon)} log Z(epsilon).

The certificate takes m=19. The sum of integer sequence weights at length k
is checked to be 3^{k-1}. Larger components can be filled with zeros, which
satisfy every inequality. All local degree claims are also checked against
the sixteen assignments of the four neighboring even bits.

The logarithm bounds use
log((1+t)/(1-t))=2 sum_{j>=0} t^{2j+1}/(2j+1), after normalization to
0<=t<=1/3. Fourteen terms and an explicit geometric remainder enclose each
logarithm; bounds are rounded outward to multiples of 2^{-40}. The resulting
axis-minus-corner and checkerboard-minus-corner gaps exceed 0.000317055 and
0.001449458, respectively. All arithmetic entering these comparisons is
integer or rational. The reproducible parameters and bounds are in
`binary-two-by-two-certificate.json`. QED.

## Large-alphabet entropy coefficients

**Proposition (working proof).** For any fixed finite nonempty S and a fixed
pattern w, h_q(X_w)-log(q) has a convergent expansion in z=1/q near z=0.
For s<=k<2s, its coefficient is

c_k(w) = sum_{[J]: |S+J|=k, J compatible} (-1)^{|J|},

where [J] ranges over translation classes of nonempty finite anchor sets.
In particular c_s=-1. This is an entropy expansion in alphabet size, distinct
from the pressure expansion in coupling at a fixed q.

Use the preceding polymer representation at lambda=-1. The elementary
sumset inequality |S+C|>=s+|C|-1 gives the anchored absolute-activity bound
e|z|^s/[1-e^2D|z|]. The Kotecky--Preiss argument applies, for example, when
|z|<1/[4e^2(D+1)], proving analytic convergence uniformly in free boxes.
Every nonempty marked-set term has degree at least s. Products of two such
terms in the formal logarithm therefore have degree at least 2s, so lower
degrees equal those of the marked partition polynomial itself. A union with
size below 2s cannot contain two disjoint copies, hence all its anchor pairs
overlap and its translation class has bounded diameter. The number of its
translates in a box divided by box volume tends to one. This proves the formula.
The sumset inequality follows by ordering S and C with a generic linear
functional: the sums of successive S-points with the first C-point, followed
by the last S-point with the remaining C-points, are s+|C|-1 distinct sums.

**Triangular dependence on displacement SETS (working proof).** Suppose S is
lattice-convex. For 1<=r<s, the coefficient c_{2s-r}(w)-A_r(w), where
A_r(w)=(1/2)|{d in G(w): |S intersect (S+d)|=r}|, depends only on S and
the compatible displacements whose overlaps have size STRICTLY greater than r.
This statement retains the displacements themselves, not just their sizes.

Indeed a marked set J with |S+J|=2s-r has every pair overlap at least r.
If every pair overlap is greater than r, compatibility is already determined
by those larger overlaps. Otherwise one pair has overlap exactly r and its
two supports cover S+J. Every other anchor lies strictly between the two
endpoints on their joining segment, by the support-function argument used in
the compression proof. Every pair except the two endpoints then has overlap
greater than r: on each line parallel to the primitive displacement, a
lattice-convex S is an integer interval, and its intersection with its shift
strictly decreases with the shift length while nonempty. If |J|>=3, choose
an intermediate anchor. Its support contains the intersection of the endpoint
supports, again by lattice convexity. Compatibility with this intermediate
copy therefore forces compatibility of the endpoints. Thus compatibility of
every term with at least three anchors is determined by the larger-overlap
displacements. The two-anchor terms contribute exactly A_r. QED.

Consequently, for fixed patterns u,v on a lattice-convex S, let r be the
largest overlap size at which their displacement sets differ. If
a=A_r(v)-A_r(u) is nonzero, then

h_q(X_v)-h_q(X_u) = a q^{-(2s-r)} + O(q^{-(2s-r+1)}) as q tends to infinity.

In particular, under strict G inclusion the leading coefficient is the
number of newly compatible displacement pairs with maximal overlap. This
quantifies the already proved large-alphabet comparison. When displacements
first differ but their counts tie, the theorem does not determine the later
coefficients or the entropy order.

For a rectangle, every class in this coefficient formula has a unique
coordinatewise-minimum-zero representative lying in an anchor box of the same
side lengths as S. `high_alphabet_checks.py --large` enumerates these sets
exactly, recording the required G mask, union size and sign. Its output is
`high-alphabet-evidence.json`.

The checks suggest a sharper structural question, using the same A_r.
Is c_{2s-r}(w)-A_r(w) determined by A_{r+1}(w),...,A_{s-1}(w) for rectangular
S? Such a triangular relation would imply that the first differing overlap
volume determines the eventual entropy order as q tends to infinity. It would
not by itself prove ordering by K at any specified small alphabet size.
This relation is NOT proved here. It passed exact checks for all binary classes
on 2x4, 3x3, 3x4, 4x4, 2x8, 3x6, 4x5; also 2x2x2, 2x2x3, 2x2x4, 2x3x3,
and 2x2x2x2 boxes. Further checks use canonical component labels generated by
at most two prescribed periods on 5x5, 6x6, 5x8 and 7x7 (respectively 591,
1,415, 1,772 and 2,825 agreement classes). They are bounded structural evidence;
these generated classes do not exhaust all patterns on the larger supports.

## Exact tests and their limitations

`ordering_checks.py` uses integer marked-occurrence/binomial moments and exact
binomial inversion to compute the occurrence-count law. Its small-case check
independently enumerates every letter configuration on a 3x3 box for all binary
2x2 agreement classes and on a 3x4 box for all binary 2x3 classes.
The commands have deterministic output; the following checks passed:

- Binary 2x2 patterns, 4x4 anchors: 6 agreement classes, 7 proper-inclusion pairs.
- Binary 2x3 patterns, 3x3 anchors: 18 classes, 53 pairs.
- Binary 2x4 patterns, 4x4 anchors: 44 classes, 264 pairs.
- Binary 3x3 patterns, 4x4 anchors: 63 classes, 465 pairs.
- Binary 3x4 patterns, 4x4 anchors: 247 classes, 3,357 pairs.
- Binary 4x4 patterns, 4x4 anchors: 1,190 classes, 30,672 pairs.
- Ternary 2x2 patterns, 3x3 anchors: 7 classes, 13 pairs.
- Ternary 2x3 patterns, 3x3 anchors: 19 classes, 71 pairs.
- 500 seeded irregular-support/anchor tests: 73,139 proper-inclusion comparisons,
  with support a subset of a 3x3 box and anchors a subset of a 4x4 box, seed
  20261020. No negative convex hinge gap was found.

These are bounded checks, not proofs of the conjecture. Trivial agreement sets
are not realizable for every shape/alphabet: no binary 2x2, 2x3, or 3x3 pattern
has G empty. Binary 2x4 does, e.g. (0011;0101). Thus any claim that the minimum
is always attained by an overlap-free pattern needs a realizability hypothesis.

`strip_scan.py` computes floating-point Perron estimates on horizontally periodic
strips for binary patterns of height two. No interval arithmetic or rigorous
planar truncation bound is supplied. Small widths can produce false impressions:
for 2x5 patterns, width 8 gives reversals of K-order which disappear at widths
10--12. These are genuine strip phenomena, not counterexamples for planar entropy.
All binary 2x3, 2x4, and 2x5 agreement classes at widths 8, 8, and 10 respectively
showed no strip entropy reversal under agreement inclusion or under strictly
ordered K. Ties of K can and do split. The consolidated `numerical-evidence.json` preserves these strip tables and
cubic comparisons, with their numerical status. Larger checks now include all 193 binary 2x6
classes at width 14, 362 binary 2x7 classes at width 16, 247 binary 3x4 classes
at width 10, and 680 binary 3x5 classes at width 12. No strict K-order or
agreement-inclusion reversal was detected at the hard strip endpoint. The
height-two quotient groups rows with identical two matching masks; the
height-three quotient retains pending one-row and two-row matches. Both were
checked against uncompressed transfers on smaller widths. These are exact
quotient constructions, with floating-point spectral evaluation.

For 2x5, 3x4, and 3x5 at those widths, comparison with beta=-1 and beta=-0.2
found no negative-half-line order reversal above the chosen numerical detection
thresholds (hard gap >1e-10 and opposite soft gap >1e-12). Small gaps below those
thresholds are unresolved. Positive-coupling reversals DO occur; the exact
3x3 obstruction below makes the qualitative point without using strip limits.

Exact cubic cumulants are computed by `cumulant_checks.py`, using connected
placements and integer moments. The known 2x3 third-derivative polynomials were
reproduced exactly for q=2,3,5. Comparing cubic tie-breaks with hard-strip values
revealed tiny apparent reversals for 3x5 at width 12 (up to 7.9e-12). The largest
selected cases reverse back to the cubic-predicted order at widths 14 and 16;
these are again boundary-size artifacts, not planar counterexamples. The
remaining sub-threshold pairs have not been certified either way.

## A new obstruction inside the positive half-line

**Proposition (working proof).** For the 3x3 patterns
u=(001;001;001) and v=(010;101;010), pressures cross at a strictly positive
coupling for every uniform alphabet size q>=2 containing symbols 0,1.
Their agreement sets are incomparable. The example has unequal K, unlike the
previous cubic crossing.

For u, G consists of +-(1,0), +-(2,0). For v, G consists of all nonzero
displacements in {-2,-1,0,1,2}^2 with even coordinate sum. Therefore
K_u=2q^6+2q^3 and K_v=4q^4+4q^3+4q. The difference is
2q^6-4q^4-2q^3-4q>0 for q>=2 (it equals 40 for q=2).
For example, q^6-2q^4=q^4(q^2-2)>=2q^4>q^3+2q for q>=2.
The covariance theorem in the recovered manuscript gives P_u>P_v for all
sufficiently small nonzero beta.

Let rho_max(w) denote maximal asymptotic occurrence density. Within each row,
001 occurrences have horizontal gaps at least 3, so rho_max(u)<=1/3. Repeating
001 horizontally and using identical rows attains 1/3. For v, horizontally
adjacent occurrences conflict, so rho_max(v)<=1/2; the infinite checkerboard
attains 1/2. Consequently P_u(beta)<=log(q)+beta/3 for beta>0, while
P_v(beta)>=beta/2. Hence P_v>P_u whenever beta>6 log(q). Continuity of pressure
forces at least one positive crossing. No numerical approximation or exchange
of derivatives with a hard-limit calculation is needed for this conclusion.

This refutes a total pressure order even on the positive half-line and refutes
using strict K-order for all positive couplings. It does not refute negative
coupling or entropy ordering by K. The two distinct controls are local
occurrence clustering near beta=0 and maximal packing density at beta=+infinity.
The known inclusion comparisons remain intact.

## A certified negative crossing when the support has a hole

**Proposition (exact finite-state certificate).** Let S={0,1,2,3,5} in Z, with
uniform binary letters, and assign values u=(0,0,0,0,1), v=(0,1,0,1,1) in the
listed site order. Then G(u)={+-1}, G(v)={+-2,+-4}, so K_u=16<K_v=20.
Nevertheless h(X_u)>h(X_v), and their pressures cross at a negative coupling.
This is also a two-dimensional example on the nonrectangular support
{(0,j): j in S}, with independent rows. It does not address rectangular supports
or agreement inclusion, because these two agreement sets are incomparable.

An occurrence is respectively a match of 0000?1 or 0101?1, where ? denotes the
unobserved site at position 4. Thus hard avoidance forbids {000001,000011} or
{010101,010111}. The prefix automata have seven states each. With t=exp(beta),
their characteristic polynomials are

F_u(x,t)=x^7-2x^6+2(1-t)x-(1-t)^2,
F_v(x,t)=x[x^6-2x^5+(1-t)(x^4-2x^3+x^2-2x+2)].

Pressure is log of the Perron root. All proper-prefix states are reachable
from the empty state before any length-six match, so no inaccessible component
can inflate this growth rate. `sparse_support_certificate.py` derives the
matrices and polynomials, checks weighted counts against direct configuration
enumeration at lengths 6,8,10, and isolates all real roots by exact rational
arithmetic. `sparse-support-certificate.json` records the rational intervals.
The certified root orders are:

| t | Perron root for u (approximate display) | for v | Certified order |
|---|---:|---:|---|
| 0 | 1.94685626827 | 1.94630604273 | h_u > h_v |
| 1/16 | 1.94941360678 | 1.94919703290 | P_u > P_v |
| 1/2 | 1.97061884346 | 1.97106986840 | P_u < P_v |

The displayed decimals are only for reading; the ordering certificate uses
disjoint rational root intervals of width less than 10^{-12}. Continuity gives
a crossing at beta in (-log(16),-log(2)). The hard entropy gap is about
0.000282662526. This rules out both strict-K entropy ordering and negative-half-
line noncrossing for arbitrary finite supports. Some geometric restriction on support is essential. This example does not
decide whether lattice convexity is enough for either proposed strengthening.
It does not refute the general-support finite convex-order conjecture under
agreement inclusion.

## Known obstruction retained

For u=(001;011), v=(000;101), G(u)={+-(1,-1)}, G(v)={+-(0,2)}, and K_u=K_v=2q^2.
The recovered working manuscript proves P_v-P_u=q^{-16} beta^3+O(beta^4).
The pair covariance sums tie, but their labeled pair fields differ. Their G's
are incomparable. No general coupling-independent total order survives this
example. The geometric difference is the count of placements of a third
occurrence disjoint from a compatible pair, not an extra triple-consistency rule.

