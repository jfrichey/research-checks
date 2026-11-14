# A finite decorated-drain inventory problem

Updated 2026-10-10 UTC. **The reduction below is proved; the proposed
inventory inequality is OPEN and stronger than the occupation conjecture.**
The project assumptions are those in [the reader’s model and hypotheses](https://raw.githubusercontent.com/jfrichey/jfrichey.github.io/companions-2026-11-14/pagedocs/ai/arw-monotonicity.pdf).
This is a possible common framework for stack and path exchanges, not a
completed injection or a proof of monotonicity. No novelty or independent
user-review claim is made.

## 1. A parameter-free finite pool

Use the complete drain representation in [SX §5](stack-exploration.md).
There are two independent systems, with labelled all-active roots eta and
`eta+delta_x`. A complete path runs from its specified root to the original
sink. Fix one canonical legal labelled FIFO convention for ordinary
stabilization and the subsequent auxiliary drain.

Attach a Bernoulli sleep coin to every nonsink vertex occurrence of every
complete labelled path. This includes occurrences reached only during the
drain, whose coins ordinary ARW never tests. Conditional on the complete
paths, these coins have the independent site-specific Bernoulli laws.
Indeed, the causal execution reads one new site coin at each arrival;
after stabilization assign the remaining path-occurrence coins in the
canonical drain order. The inverse runs that execution on the decorated
paths and reconstructs each site's used coin prefix and its unused tail.
Every permutation is within a site. This is the coin analogue of SX's
arrow-book/path correspondence; it retains all ignored coins.

Let `J(z,w)` be the combined number of directed steps `z->w` in the two
complete path families, including exits to the sink s. Write
`m_z=sum_w J(z,w)`. A feasible J has a finite decomposition into paths
from the prescribed labelled roots, all ending at s. Let `d(J)` be the
number of such decompositions, with repeated steps to the same destination
indistinguishable. For each site prescribe a total `k_z` of sleep coins
equal to one among its m_z path occurrences.

The resulting pool `Omega(J,k)` contains

`d(J) product_z binomial(m_z,k_z)`

decorated path pairs. It is finite and uniform under the original product
law conditional on J,k: every atom has the same factor

`product_(z,w) K(z,w)^J(z,w) product_z p_z^k_z (1-p_z)^(m_z-k_z)`.

Unrestricted residual books are integrated out using SX's full inverse.
In particular the common paths in the two copies are **not** constrained
to be identical. That restriction already fails in SX's one-common-particle
example. The pool also differs from F3's combined *ordinary consumed*
inventory: it includes complete drains and their ignored sleep coins.

## 2. The proposed inequality and its implication

For an atom omega, run each ordinary system from its labelled roots using
its decorated paths, and write `I^-(omega), I^+(omega)` for its final target
indicators. Define the integer

`c(J,k) = sum_(omega in Omega(J,k)) [I^+(omega)-I^-(omega)]`.

**Proposed sufficient statement:** `c(J,k)>=0` for every feasible J,k and
every binary common input A with an extra site x outside A. Equivalently,
each pool has at least as many good `(0,1)` atoms as bad `(1,0)` atoms.
This assertion is not proved, even for all pools with two common particles.

If proved, conditioning on J,k gives the desired occupation inequality for
binary inputs. Almost-sure finite drains partition the full probability
space into countably many such pools; bounded indicators justify summing.
The saved surplus-IDLA reduction then gives arbitrary finite all-active
inputs. There is no truncation or unproved normalization in this reduction.
The inequality does not assert occupied-set domination or chance-count order.

It is useful to first omit the coin-count conditioning. The conditional
increment at fixed J is the polynomial

`Gamma_p(J) = (1/d(J)) sum_k c(J,k) product_z p_z^k_z (1-p_z)^(m_z-k_z)`.

Coefficientwise positivity is stronger than positivity of this polynomial
at a particular p, and stronger than the original inequality after averaging
over J. A negative coefficient would refute this proposed pool statement,
not necessarily either weaker inequality or the main conjecture.

### A stronger pool using one complete path family

There is also a smaller, **OPEN** sufficient inequality. Take one family of
complete decorated paths with roots `eta+delta_x`, and designate the extra's
label. Run the upper system on all paths, and the lower system on just the
common paths, reusing their path-occurrence coins. Each marginal has the
correct ARW law. Condition on the edge inventory J and coin totals k of this
one family, and call its signed count `c_1(J,k)`. Its root divergence is
`eta+delta_x`, whereas the earlier two-family J has divergence `2 eta+delta_x`.
Sharing a whole fixed common path still fails as in SX; this pool allows
its geometry to change while retaining only the aggregate inventory.

With the walk generating functions h_z and H_eta defined in Section 6,
the two signed generating series are related by

`C_1(X,p)=H_eta h_x [Q_v^(K^X)(eta+delta_x)-Q_v^(K^X)(eta)]`,
`C(X,p)=H_eta C_1(X,p)`.

The unused coins on the extra H_eta paths supply positive binomial degree
elevations. Hence coefficientwise `c_1>=0` would imply the earlier inventory
inequality and the main conjecture. It is a stronger candidate, not a
consequence of the two-family finite checks.

Section 6's endpoint theorem and its no-external-return corollary hold also
for c_1. In the proof of (3)--(5), remove the single extra factor H_eta from
each expression. The remaining factor `h_C h_y` still contains h_v in
every nonzero term, so the same positive homogenization applies. The
residual drain roots in (5) are now C together with y. No new sign assumption
is used in this extension.

These counts need only one ordinary-buffer recursion for each marginal.
For the upper, process all initial labels and close the residual drain by
(1). For the lower, process the common labels and retain the extra as an
unprocessed drain root; the residual inventory's divergence includes that
root automatically. At either ordinary stopping state, multiply the drain
count by the binomial choices of its untested coins. Thus both marginal
counts use the same full J and coin degree. For Section 3's example with
`J(0,s)=2` instead of 3, there are three path families and
`sum_k c_1(J,k)t^k=t_0(1+t_0)^2 t_2^2`.

## 3. A complete small cancellation

Take lower root 0 and upper roots 0,1. Fix

`J(0,2)=J(2,0)=2, J(0,s)=3, J(1,0)=1`,

with all other counts zero. The three labelled walks share two excursions
`0->2->0`, distributed with counts `(a,b,c)` among the lower walk, upper
common walk, and upper extra walk; `a+b+c=2`. There are six decompositions.
The extra walk begins `1->0`, and each walk finally exits from 0.

One decomposition has the common path `0->2->0->s` in both copies and
extra path `1->0->s`. At p=1/2 its conditional target increment is `-1/16`,
as in SX. Nevertheless the **whole** six-decomposition pool has

`Gamma_p(J) = p_0(1-p_0)(1-p_1)p_2^2 / 6 >= 0`.

The nonzero integer coefficients c(J,k) are precisely

| `(k_0,k_1,k_2)` | c(J,k) |
|---|---:|
| `(1,0,2)` | 1 |
| `(2,0,2)` | 3 |
| `(3,0,2)` | 3 |
| `(4,0,2)` | 1 |

Thus the signed count polynomial is `t_0(1+t_0)^3 t_2^2`. This is an exact
finite cancellation, not yet a general cycle-transfer rule.

## 4. Counting by last-exit forests

There is a smaller exact calculation which avoids listing complete path
families. For a nonnegative edge inventory J, put
`r_z=m_z-sum_w J(w,z)`, excluding the sink from the incoming sum. If some
r_z is negative, there is no decomposition into paths starting at r.
Otherwise let U be the vertices with m_z>0, and let L be the matrix on U
with entries `L(z,w)=m_z 1{z=w}-J(z,w)`. Then

`d(J) = det(L) product_(z in U) (m_z-1)! / product_(z,w) J(z,w)!`.       (1)

Empty products and the empty determinant are one. A positive incoming
count at a vertex with no outgoing edge is infeasible. Formula (1) counts
ordered labelled paths for any fixed ordering of roots with multiplicities
r; identical destinations are not distinguished.

**Proof.** Put the J(z,.) destinations in a local word at each z. Route
the labelled roots using these words and a fixed legal arrow-only schedule.
It uses every word exactly when the last destinations form a directed
forest toward the sink. For necessity, a cycle of last departures would
force their departure times to increase strictly around that cycle. For
sufficiency, if routing stopped with unused edges, their residual inventory
would have zero divergence. Its support would contain a cycle of the
remaining last edges, contradicting the forest condition. Routing cannot
require an unavailable edge: at any prefix, the remaining divergence
equals the nonnegative pending-particle vector.

For a chosen last-exit forest T, the number of local words is
`product_z (m_z-1)! / product_(z,w) (J(z,w)-1{T(z)=w})!`.
Factoring out the factorials in (1) leaves the sum of
`product_z J(z,T(z))` over forests T. Determinant expansion gives this sum:
terms containing a directed cycle cancel, leaving precisely the forests
to the sink. Finally, routing and reconstructing the local words are
inverse operations between these words and the labelled path families.
This proves (1).

For the signed count polynomial, recursively expose ordinary buffer
arrivals in the lower copy and then the upper. A state records the remaining
J, each queue, and every resident's waiting destination. At a delivery to z,
choose one distinct available destination d, subtract its edge, and branch
over the two bits. Emitting an old resident's waiting destination needs no
further edge choice. The sleep branch contributes t_z. Destinations carry
weight one, not their remaining multiplicity, because identical edges were
never individually labelled.

When both ordinary queues are empty, the sign delta of that state is known.
Force all waiting resident arrows. The residual inventory J' has divergence
equal to these drain starting positions; this follows from conservation at
every preceding arrival. Its exact contribution is

`delta d(J') product_z (1+t_z)^(sum_w J'(z,w))`.                       (2)

Thus (1) closes the recursion, retaining all untested coins. A fixed ordinary
prefix and an arbitrary complete labelled drain uniquely reconstruct the
original decorated paths, so (2) neither discards nor repeats histories.
This proves the counting method, not positivity of its signed sum.

Equivalently, route and drain the entire lower copy before starting the
upper. Concatenating the two consumed arrow words at each site gives one
pooled word whose last arrows form a forest. The first copy's split lengths
are determined by arrow-only routing and are independent of every sleep
bit. Conversely, a pooled word satisfying the forest condition, together
with its m_z coins at each site, recovers both complete copies. This gives
a second finite description of Omega(J,k): admissible local arrow words
and sitewise coin words with their prescribed totals. Any proof may change
the legal schedule, provided it retains this full bijection.

**The forest cannot be frozen.** Take lower root 0, upper roots 0,1, and
`J(0,s)=2`, with one edge each `0->1, 1->s, 1->2, 2->0`. There are five
admissible pooled words. Grouping their signed count polynomials by last
destinations `(T(0),T(1),T(2))` gives

| Last destinations | Word count | Signed count polynomial |
|---|---:|---|
| `(s,s,0)` | 2 | `t_0(1+t_0)(t_1^2-1)` |
| `(s,2,0)` | 2 | `t_0(1+t_0)(1+t_1)` |
| `(1,s,0)` | 1 | `0` |

The full polynomial is `t_0(1+t_0)t_1(1+t_1)`, coefficientwise nonnegative.
The first forest has negative coefficients at `(k_0,k_1,k_2)=(1,0,0)` and
`(2,0,0)`. Its conditional increment is
`p_0(1-p_0)(2p_1-1)(1-p_2)/2`, negative when p_1<1/2. In contrast, the
full pool has increment `p_0(1-p_0)p_1(1-p_2)/5`. These identities are
verified by enumerating the five words and all 64 assignments of their
six coins. Only the refinement by one fixed last-exit forest is refuted;
the proposed inequality for complete J,k remains open. An injection must
permit the last-exit forest to change.

## 5. A cancellation which stays inside every pool

Suppose both ordinary histories visit v and their final target bits differ.
Let n,m be their last target arrival indices. Exchange **only** their bits
`B_v^-(n),B_v^+(m)`, keeping all potential arrows fixed. Accept this exchange
exactly when the new ordinary histories have the same two terminal indices;
otherwise leave the pair unchanged.

**Proposition.** This is a measure-preserving involution within Omega(J,k).
Every accepted exchange reverses the target sign.

**Proof.** At unchanged positive terminal indices the exchanged bits are
the new final occupation indicators. The same indices are therefore
recovered from the image, and exchanging them again restores both complete
books. Acceptance is invariant, so rejected points cannot collide with an
accepted image. All D-books are fixed; hence both complete drain lengths
and their entire edge inventories are fixed. The two changed coordinates
lie inside those drain prefixes, since both arose at ordinary arrivals.
One target one-bit and one target zero-bit exchange places, preserving k
and the joint product weight. Unused suffixes are unchanged.

The jump-to-sleep change always keeps its last index, by the suppression
argument in [RN §6](response-renewal.md). Thus acceptance is equivalent to
the opposite forced-deletion avalanche avoiding v. This is a partial
cancellation only: deletion can return to v, and the empty-target copy may
never have visited v. The complete-avalanche maps in AX may change J and
cannot automatically be used as moves in these fixed pools.

## 6. Positivity at all off-target coin-count endpoints

**Theorem.** Fix any finite all-active input eta, extra site x and target v.
For every feasible J, suppose `k_z` is either 0 or m_z at each z other than
v. Then `c(J,k)>=0`, for every `0<=k_v<=m_v`.

The assumption allows a different endpoint choice at every off-target site.
It concerns finite coin counts, not an extension of the main conjecture's
strict sleep assumptions. In the proof, evaluating a count polynomial at
sleep probabilities zero or one selects exactly these endpoint coefficients.

Let S be the off-target sites with k_z=m_z, put `T=S union {v}` and
`R=V\T`, and write p for the target sleep probability and q=1-p. Sites in
S sleep immediately when alone; sites in R never sleep. An IDLA procedure
first routes eta to the original sink or unoccupied sites of T. Its settled
cluster C is a subset of T. This is a legal initial phase of ARW: hold one
active particle at each reached T-site, move surplus particles there, and
always move particles in R. The retained particles then sleep immediately
at C\{v}; if v is present, activate its retained particle with parameter p.
The unused instructions are fresh. Equivalently this is the saved surplus
IDLA argument with R eliminated by tracing walks to T or the sink.

For any C, let r_C be the probability that a forced departure from v returns
to v before reaching the original sink or a vacant T-site. Its first-return
path has all interior vertices in `R union (C\{v})`. With v initially active
and the other C-sites sleeping, there is only one moving particle until it
stops. An excursion returning to v preserves C\{v}. Therefore the conditional
target payoff after IDLA is

`F_p(C) = 1{v in C} p/(1-q r_C)`.

It is increasing in C. Running the common IDLA particles first and the extra
last couples C with `C union {y}`, where y is the extra's newly settled site;
if the extra exits, leave C unchanged. This already proves occupation
monotonicity for these boundary sleep parameters. The next argument proves
the stronger assertion coefficient by coefficient in J and k_v.

Give each directed edge an indeterminate weight X(z,w), and let h_z(X) be
the sum of the weights of all complete walks from z to the sink; set h_s=1.
For sufficiently small positive X these sums converge and
`K^X(z,w)=X(z,w)h_w/h_z` is a stochastic kernel. A walk's K^X-weight is
its X-weight multiplied by h_endpoint/h_start. Write
`H_eta=product_z h_z^eta_z` and `h_C=product_(z in C) h_z`.
The generating series of the joint signed counts is consequently

`sum_J X^J sum_(k_v) c(J,k) p^k_v q^(m_v-k_v)`
`    = H_eta^2 h_x [Q_v^(K^X)(eta+delta_x)-Q_v^(K^X)(eta)]`.           (3)

Here off-target k is the selected endpoint vector. This identity follows
by weighting every complete decorated path, not by assuming an IDLA
coupling preserves a particular drain inventory.

Let gamma be a complete common-plus-extra IDLA history, with common cluster
C and extra endpoint y. Its X-edge inventory is E_gamma. Telescoping the
h-factors in its K^X-probability changes the prefactor in (3) to
`X^E_gamma H_eta h_C h_y`. A nonzero payoff difference has one of two forms:

- If v is absent from C and y=v, take every finite closed walk W from v
  whose vertices lie in `R union C union {v}`, including the empty walk.
- If v belongs to C and y is a new off-target T-site, take every finite
  closed walk W from v in `R union C union {y}` which visits y.

All other histories contribute zero. Let ell(W) count W's departures from
v, and E_W be its edge inventory. Expanding `p/(1-q r_C)` into successive
first-return excursions gives precisely the positive terms

`X^(E_gamma+E_W) H_eta h_C h_y p q^ell(W)`.                           (4)

In the second case the permitted W are exactly the excursions available
after adding y but unavailable before it; their coefficients are nonnegative.
Closed-walk h-factors cancel, so these are ordinary X-path weights.

For a fixed J and chosen gamma,W, the coefficient of the remaining edges
in `H_eta h_C h_y` is `d(J-E_gamma-E_W)`, by the path count in Section 4.
Its labelled roots are eta together with one at each site of C and one at y.
Flow conservation gives precisely that residual divergence. Whenever the
contribution is nonzero, this collection has a root at v, since either
v is in C or y=v. Thus `b=m_v-ell(W)>=1`. Insert
`1=(p+q)^(b-1)` into (4). The target factor becomes a homogeneous polynomial
of degree m_v with nonnegative coefficients. Uniqueness of the basis
`p^k q^(m_v-k)`, `0<=k<=m_v`, on q=1-p now gives the explicit formula

`c(J,k) = sum_(gamma,W) d(J-E_gamma-E_W)`
`                    * binomial(m_v-ell(W)-1,k_v-1) >= 0`.           (5)

The sum uses the two permitted history/walk classes above and nonnegative
residual inventories; infeasible residuals and out-of-range binomials
contribute zero. It is finite for a fixed J. The convergent generating
identities near X=0 justify coefficient extraction and hence (5).
If m_v=0, both target indicators vanish directly. This completes the proof.

In particular, when S is empty and eta initially contains v, the sum has
no permitted histories and all these coefficients are zero. For the
six-path example in Section 3, the positive corner with S={2} has one
contributing IDLA prefix and one marked return `0->2->0`; formula (5)
gives the coefficients `binomial(3,k_0-1)` directly.

**Remaining gap:** an off-target site with `0<k_z<m_z` has both kinds of
coin. It no longer acts as a permanent IDLA settling site or a permanent
passage. The endpoint theorem does not control these interior Bernstein
coefficients; interpolation of endpoint values would not prove their sign.

### Mixed sites without a return through another vertex

**Corollary.** The inventory inequality also holds if every mixed off-target
site lies on no directed cycle involving another vertex in the support of J.
Self-loops at these sites are permitted. Other sites, including the target,
may belong to arbitrary directed cycles.

First suppose a selected mixed site w has no self-loop. Process all its
upstream particles before processing w. No output from w can affect this
incoming count, since that would give a directed return to w. Conditional
on its N arrivals, its final bit is a fresh Bernoulli(p_w) when N>0; its
output is the first `N-B_w(N)` arrows. All other local bits can be discarded
for this calculation. Thus its entire effect on the target law is affine
in p_w. This holds separately for every such w, even with arbitrary sleep
parameters elsewhere. Work on J's directed support and apply the generating
identity (3): each fixed-J target coefficient has the same affine dependence
before homogenizing to degree m_w. If A and B are its endpoint coefficients
at k_w=m_w and k_w=0, respectively, the interior coefficient is

`A binomial(m_w-1,k_w-1) + B binomial(m_w-1,k_w)`.

This uses actual affinity, not interpolation of an arbitrary polynomial.
Applying it successively reduces all mixed sites to Section 6's endpoints.

For self-loops write `sigma=X(w,w)`, let ell be their count in J, and let
`r=m_w-ell` be the nonself outgoing count. Feasibility and a visit to w
give r>=1. Removing self-loops changes the sleep parameter to
`p_w/(1-(1-p_w)sigma)`. In the edge generating series it replaces every
nonself edge weight out of w by `X(w,z)/(1-sigma)`. This follows either by
summing repeated self-loops or by the Doob normalization used in (3).
Consequently the fixed nonself inventory contributes

`[sigma^ell] (1-sigma)^(-r)`
`       * [p_w A + (1-p_w)(1-sigma) B] / [1-(1-p_w)sigma]`.          (6)

Both parts have positive expansions: their factors in sigma are respectively
`(1-sigma)^(-r)` and `(1-sigma)^(-(r-1))`, followed by the geometric series
with ratio `(1-p_w)sigma`. A term using i self-loops from the first factor
and j from the second has i+j=ell and degree j+1 in p_w,1-p_w.
Multiply it by `(p_w+(1-p_w))^(r+i-1)` to reach degree m_w. The exponent
is nonnegative and every resulting coefficient is nonnegative. For several
selected sites these factors combine independently after all their self-loops
are removed. This proves the corollary.

The unresolved mixed sites therefore require a return through another vertex.
A simple monotone interpolation rule is false even though its counts remain
positive: for common `{0,2}`, extra 1, take
`J(0,2)=J(0,s)=J(2,s)=2` and
`J(2,1)=J(1,s)=J(2,0)=J(1,0)=1`.
At k_0=1,k_1=2, the normalized counts
`c(J,k)/binomial(4,k_2)`, for k_2=0,...,4, are
`10,17/2,8,17/2,10`. The mixed site 2 lies on a directed return cycle,
so this does not contradict the corollary.

## 7. One off-target site with nonzero sleep counts

**Theorem.** Both c_1 and c are nonnegative whenever all off-target sleep
counts are zero except possibly at one site w. Its count may be strictly
between zero and m_w, and w may lie on directed return cycles. The input
eta may have arbitrary finite multiplicities. This does not cover an
additional off-target site with a positive, even maximal, sleep count.

First consider two sites, target 0 and site 1, without self-loops. Write
`p_i` for their sleep probabilities, `q_i=1-p_i`, and
`u=K(0,1), w=K(1,0)`. Let Q_0,Q_1,Q_01 be target occupations from those
binary inputs. Put `D=1-q_0 q_1 u w` and
`D_2=1-(1-p_0 p_1)u w`. Single-particle equations give
`Q_0=p_0/D`, `Q_1=p_0 q_1 w/D`. Expanding the two-site complete-block
inverse from CB gives

`p_0 p_1 = D_2 Q_01`
` - q_1[1-u w-q_0(1-u)w] Q_0`
` - q_0[1-u w-q_1 u(1-w)] Q_1`.

For example, this follows by using
`a_0^2=(1-u)a_0+u a_0 a_1` and
`a_1^2=(1-w)a_1+w a_0 a_1` in that inverse. Substitution and simplification
give

`Q_01-Q_0 = p_0 p_1 q_0 u w(1-q_1 w)/(D D_2)`,
`Q_01-Q_1 = (Q_01-Q_0)+p_0(1-q_1 w)/D`.                            (7)

To obtain coefficient positivity, use edge variables
`alpha=X(0,1), beta=X(1,0), a=X(0,s), b=X(1,s)` and t=alpha beta.
Then `h_0=(a+alpha b)/(1-t)` and `h_1=b+beta h_0`. Under K^X,
`u w=t` and `h_1(1-q_1 w)=b+p_1 beta h_0`. Therefore the two one-family
signed generating functions are

`C_(0;1) = p_0 p_1 q_0 t h_0 (b+p_1 beta h_0)`
`          / [(1-q_0 q_1 t)(1-(1-p_0 p_1)t)]`,
`C_(1;0) = C_(0;1) + p_0 h_0(b+p_1 beta h_0)/(1-q_0 q_1 t)`.         (8)

Expand both denominators geometrically, writing
`1-p_0 p_1=q_0 q_1+p_0 q_1+q_0 p_1`. Each such denominator step consumes
one edge from each site and carries one coin factor from each site. All
coefficients are nonnegative. In every numerator term, the number of
specified coin factors at each site is at most the number of its edge
departures: use the displayed sink/connecting edges and the obligatory
first departure of each h_0 path. Multiply by the required powers of
`p_i+q_i` to fill the unused degrees. Thus (8) proves coefficientwise
positivity for both binary one-common-particle comparisons.

If both sites occur in eta, its IDLA reduction fills both sites and every
surplus particle exits; the increment is zero. If eta consists of n>=1
particles at site i alone, its IDLA cluster remains {i} exactly when all
n-1 surplus walkers exit before reaching the other site j. In the loopless
graph this probability is `K(i,s)^(n-1)`. After multiplying by the one-family
walk factor, the generating function for adding at j is
`X(i,s)^(n-1) C_(i;j)`; adding at i supplies the additional positive factor
`X(i,j)`. The empty-input case has only a positive single-particle payoff.
These observations extend the coefficient result to every two-site input.

**Self-loop insertion preserves coefficient positivity generally.** This
fact does not require the selected vertex to avoid other return cycles.
In a fixed inventory without self-loops at z, let r be its outgoing count
and k its sleep count. Insert a self-loop variable sigma there. The change
of kernel and sleep parameter used in (6) sends its coin monomial to

`p_z^k q_z^(r-k) (1-sigma)^(-k) (1-q_z sigma)^(-r)`.                (9)

The coefficients in sigma are nonnegative. A term using i loops from the
first factor and j from the second has coin degree r+j; multiply it by
`(p_z+q_z)^i` to reach the new degree r+i+j. Thus a nonnegative coefficient
array stays nonnegative after any prescribed self-loop counts are inserted.
This applies to either one-family or two-family generating series. It
completes the two-site result with arbitrary self-loops.

Finally, allow arbitrarily many vertices at which all coins are zero.
Trace each initial particle to its first hit of `{v,w}` or the original
sink before moving the retained particles. Let gamma be these initial
prefixes, let xi be the resulting common multiset on `{v,w}`, and let y
be the extra's endpoint. Trace all subsequent walk steps between these two
sites in the same way. The traced edge weights Y are positive power series
in X. Each traced step has exactly one departure from its starting retained
site and no other retained-site departure before its endpoint.

The initial-prefix K^X-weight is
`X^E_gamma H_xi h_y/(H_eta h_x)`. Consequently

`C_1(X,p) = sum_gamma X^E_gamma C_1^(traced)(Y,p;xi,y)`,              (10)

with zero contribution when y is the sink. Here `H_xi` uses the traced
walk functions, which agree with the original h at v,w. The two-site
theorem applies to each xi, including multiplicities. Positive substitution
of Y preserves the edge coefficients and their retained-site coin degrees;
all other coins are fixed zero. This proves the asserted result for c_1.
Multiplication by H_eta and positive unused-coin padding proves it for c.

The next cyclic case has a mixed off-target site **and** another off-target
site with positive sleep count. Unlike a never-sleep vertex, that latter
site cannot be eliminated by the trace identity (10).

## 8. A permanent extra on a directed three-cycle

Here the additional positive-count site can be included without restricting
the mixed count. Take nonsink edges only `0->1, 1->2, 2->0`, together with
arbitrary sink exits. The common input is `{0,1}`, the extra site is 2, and
all site-2 coins equal one. Then c_1 and c are nonnegative for either target
0 or 1 and every inventory J; the counts at 0 and 1 may both be interior.
Self-loops may subsequently be inserted using (9).

Write a_i for the jump probability from i to its cyclic successor, and
put `t=a_0 a_1 a_2`, `D=1-q_0 q_1 t`,
`E=1-(1-p_0 p_1)t`. In the upper system, site 2 is permanently occupied,
so moving particles trace through it. The weak sites have the two-site
kernel with forward probability a_0 and backward probability a_1 a_2.
Thus (7) gives its two target occupations.

For the lower, process the initial particle at 0 first. If it sleeps,
the other particle can only sleep at 1, exit, or fill the vacant site 2.
If the first particle instead jumps to 1, there are two active particles
there. Their output consists of one forced jump and, with probability q_1,
a second jump. A further weak-site arrival is possible only if both outputs
reach 2, where the second arrival releases one particle toward 0. At that
point 2 is occupied and only one weak particle is moving. This gives

`L_0=p_0+p_0 q_0 q_1 a_0 a_1^2 a_2/D`,
`L_1=p_1+p_1 q_0^2 q_1 a_0^2 a_1^2 a_2/D`.

Subtracting from the upper two-site occupations yields

`Delta_0=p_0 q_0 q_1 t(1-a_1)/D`
`          +p_0 p_1 q_0 t(1-q_1 a_1 a_2)/(D E)`,
`Delta_1=p_1 q_0 q_1 t(1-q_0 a_0 a_1)/D`
`          +p_0 p_1 q_1 t(1-q_0 a_0)/(D E)`.                        (11)

These already prove the occupation inequalities in this family. For the
stronger inventory sign, use edge weights `r_i=X(i,i+1)` cyclically and
`d_i=X(i,s)`, so `h_i=d_i+r_i h_(i+1)` and now t=r_0 r_1 r_2.
Multiplying (11) by h_0 h_1 h_2 and using these identities gives

`C_(target 0)=p_0 q_0 t h_0 h_2`
`    * [q_1 d_1/D + p_1(d_1+r_1 d_2+p_1 r_1 r_2 h_0)/(D E)]`,
`C_(target 1)=p_1 q_1 t h_1 h_2`
`    * [q_0(d_0+r_0 d_1+p_0 r_0 r_1 h_2)/D`
`       +p_0(d_0+p_0 r_0 h_1)/(D E)]`.                             (12)

As in (8), expand D and E geometrically with
`1-p_0 p_1=q_0 q_1+p_0 q_1+q_0 p_1`. Every cycle supplies one departure
from each weak site. In each numerator term, its explicit edges and the
first steps of its h-paths provide at least as many departures at each
weak site as specified coin factors. Positive powers of p_i+q_i supply
the unused coins. All site-2 coins are fixed one. This proves coefficient
positivity of (12), hence the stated one-family inequality; the two-family
inequality follows by multiplying by H_eta.

This proof uses the directed-cycle geometry. An extra return direction can
leave two weak particles moving before the vacant site 2 is first filled,
so the displayed lower formulas cannot be reused without a new argument.

## 9. General three-site kernels: an exact certificate with finite padding

**Proved by an exact symbolic certificate:** with arbitrary accessible-sink
three-site K, common input `{0,1}`, extra site 2, p_2=1 and 0<p_0,p_1<1,
the occupation increment is nonnegative at either weak target. Moreover,
for target 0 a positive inventory expansion exists after adding ten auxiliary
complete walks which neither ordinary ARW execution uses. This does **not**
prove the unpadded c_1 or c inequality for general three-site kernels, and
does not cover an interior sleep count at site 2.

Remove self-loops initially. Write the rows of K, with sink last, as

`(0,a,b,1-a-b)`, `(d,0,e,1-d-e)`, `(g,h,0,1-g-h)`.

Put q_i=1-p_i, `s_0=b g`, `s_1=e h`, `u=a+b h`, `w=d+e g`,
`A_i=1-q_i s_i`, `B_i=1-s_i`, and

`G=A_0 A_1-q_0 q_1 u w`,
`D=1-q_0 q_1 a d`, `E=1-(1-p_0 p_1)a d`,
`M=A_0 A_1 B_0 B_1-(A_0 A_1-p_0 p_1)u w`,
`S=1-q_1 d(a+b)-q_0 q_1 a d(e+p_1 d)`,
`R=q_1[p_0 d b A_1+q_0(a e+b)w]`.

The lower and upper target-0 occupations are exactly

`L=p_0(S G+R D)/(E G D)`,
`U=p_0[A_1 M+p_1 q_0 u w B_0(A_1-q_1 w)]/(G M)`.                  (13)

For the upper, trace through its permanently occupied site 2 and apply
the two-site formula (7), including self-loop removal. For the lower,
stop at the first arrival to the initially vacant site 2. Before that
arrival, the two weak particles see 2 as a sink. Afterwards at most one
weak particle remains, and it sees the upper traced kernel. Their single
particle occupations are `p_0/D,p_0 q_1 d/D` before that event and
`p_0 A_1/G,p_0 q_1 w/G` afterwards. First-step equations for two active
particles at 0, at 1, or one at each, and for one active plus one weak
sleeper, give

`E L = p_0(1-q_1 d)`
` +p_0 q_1 d[b(p_0 A_1/G)+(1-a-b)(p_0/D)]`
` +q_0[(a e+b)(p_0 q_1 w/G)`
`       +(a(1-d-e)+1-a-b)(p_0 q_1 d/D)]`.

This simplifies to (13). In particular these are finite first-step
identities, not fitted formulas from parameter samples.

Define the polynomial P by the exact division

`q_0 P = M(E D A_1-S G-R D)`
`          + E D p_1 q_0 u w B_0(A_1-q_1 w)`.

Then `U-L=p_0 q_0 P/(E D G M)`. All four denominator factors are positive.
For M, explicitly

`M=A_0 A_1 B_0 B_1(1-T)`,
`T=u w[p_0 q_1/(A_0 A_1 B_0)+q_0/(A_0 B_1)]`.

The two terms of T are disjoint return histories to the state with one
active particle at each weak site, after summing self-loops. The event
that both particles sleep before such a return has positive probability,
so T<1. The other factors are single-particle or two-site killed-walk
denominators and are positive under the stated assumptions.

Here is the finite nonnegativity certificate for P. It has 1,603 ordinary
integer monomials. Homogenize its three jump rows to degrees `(4,5,4)`,
using each row's sum as one, and express its sleep variables in the basis
`p_0^i q_0^(2-i) p_1^j q_1^(4-j)`, with i=0,...,2 and j=0,...,4.
The resulting polynomial has **47,808 positive integer coefficients and
23,067 zero coefficients**, with no negative coefficient. Every basis
monomial is nonnegative for every admissible jump row and sleep parameter.
This proves P>=0 and hence the asserted occupation comparison. Relabelling
0 and 1 gives the other target.

The standard-library verifier [verify_three_site_permanent.py](verify_three_site_permanent.py)
constructs P from the displayed expression, performs and checks the exact
division by q_0, and computes every coefficient by integer arithmetic.
Its sorted nonzero coefficient table has SHA256
`a6a5d7fbef290ee1350455e0e8a7bee98b162c3b9240c46882c30d0598eb4a8d`.
The certificate is a uniform polynomial identity and nonnegativity check;
the separate 24 rational raw-ARW comparisons only audit its connection to
the occupation formulas. No independent user-validation claim is made.

**Why ten auxiliary walks suffice for an inventory expansion.** Multiply
the increment by `h_0^4 h_1^5 h_2^4` and substitute K^X. Each homogenized
jump monomial has exactly 4,5,4 departures from the respective sites. Its
negative h-powers cancel against this prefactor, leaving positive walk
tails from its jump destinations. The external p_0 q_0 and certificate's
sleep basis have degrees 4 and 4 at sites 0 and 1, within the available
4 and 5 departures. Supply the remaining coins by positive powers of
p_i+q_i; all site-2 coins are one.

Every denominator also has a positive closed-walk expansion. For example,

`G=1-q_0 b g-q_1 e h-q_0 q_1(a d+a e g+b h d)`.

D and E expand as in Section 7, and the displayed expression for M uses
the positive geometric expansions of A_i, B_i and 1-T. Each loop specifies
at most one coin per weak-site departure; unused departures receive the
same positive padding. Closed-walk h-factors cancel. Hence all inventory
coefficients of this augmented generating series are nonnegative.

The root multiset `(4,5,4)` consists of the three designated ordinary
labels at 0,1,2 and ten unused labels, with counts `(3,4,3)`. Route these
auxiliary labels only in the drain. Their complete paths and all their
ignored coins remain part of the finite pool, with their original product
weights. They do not change either ARW marginal. Formula (9) then permits
arbitrary self-loops as well. This is explicit finite padding, not selection
of favorable unused paths or an unproved normalization.

## 10. All sleep parameters on a directed three-cycle

**Proved by an exact integer certificate:** the occupation conjecture holds
for every finite all-active input on the directed three-cycle, with arbitrary
site-dependent sleep probabilities in (0,1), sink exits and self-loops.
This removes Section 8's permanent-extra assumption for the occupation
comparison. Section 11 proves the stronger unpadded inventory inequality
for this geometry, including all three interior sleep counts.

After removing self-loops, let i jump to i+1 modulo 3 with probability a_i,
or to the sink with probability 1-a_i. Set q_i=1-p_i and t=a_0 a_1 a_2.
It suffices to treat target 0. Write Q_i for its occupation probability
from one active particle at i, B_i for the input {i,i+1}, and U for the
full input. The single-particle equations give

`D_1=1-q_0 q_1 q_2 t`,
`(Q_0,Q_1,Q_2)=p_0(1,q_1 a_1 q_2 a_2,q_2 a_2)/D_1`.

The pair occupations satisfy B=T B+g, where the only nonzero entries are

`T_(i,i+1)=q_i p_(i+1) a_i a_(i+1)`,
`T_(i,i-1)=q_(i+1) a_(i+1)`,
`g_i=p_i p_j 1{0 in {i,j}}+q_j(1-a_j)Q_i`
`       +q_i p_j a_i(1-a_j)Q_j+q_i p_j(1-a_i)1{j=0}`, j=i+1.       (14)

For clarity, this is an identity from the complete-block inverse, not a
conditional scheduling choice. Its pair row initially has the additional
single-particle expression `q_i(1-a_i)[Q_j-q_j a_j Q_(j+1)]`.
The single-particle equation replaces the bracket by `p_j 1{j=0}`,
giving (14). The remaining pair coefficients are exactly the two displayed
entries of T. Every row sum of T is at most 1-p_i p_j<1. Hence its inverse
is positive and its determinant is positive:

`D_2=det(I-T)`
`   =1-t[q_0 q_1 q_2+sum_i p_i q_(i+1)q_(i+2)]`
`       -t^2 p_0 p_1 p_2 q_0 q_1 q_2`.

To specify the full-input formula without a long expanded polynomial, let
R_S(C) be the IDLA final-support law obtained by forcing one jump from each
site of S in the full binary input, then routing all surplus particles.
Put `F_C=sum_(S subset V) (-1)^|S| q_S R_S(C)`. The complete-block inverse
gives `sum_C F_C Q_C=p_0 p_1 p_2`, with Q_empty=0 and Q_V=U.
For nonempty S the only way to retain full occupancy is to traverse the
whole directed cycle, with probability t. Thus

`F_V=D_3=1-(1-p_0 p_1 p_2)t`,
`U=[p_0 p_1 p_2-sum_(C proper subset V) F_C Q_C]/D_3`.              (15)

Equations (14)--(15) define finite polynomials
`N_x=D_1 D_2 D_3[U-B_(x+1)]`, where x is the extra site. The verifier
constructs each F_C by finite IDLA enumeration, uses Cramer's rule in (14),
and forms N_x with integer polynomial arithmetic. Extract p_0 from N_0,
and `t p_0 q_0 p_x` from N_x for x=1,2. In variables
`(p_0,p_1,p_2,a_0,a_1,a_2)`, the residual polynomial has the following
nonnegative unnormalized Bernstein coefficients:

| Extra x | Basis degrees | Positive coefficients | Negative |
|---|---|---:|---:|
| 0 | (2,4,4,3,4,4) | 7,027 | 0 |
| 1 | (2,3,4,2,3,3) | 2,756 | 0 |
| 2 | (1,3,2,2,3,3) | 982 | 0 |

Explicitly, a monomial `c z^i` contributes
`c product_r binomial(d_r-i_r,j_r-i_r)` to the coefficient of
`product_r z_r^j_r(1-z_r)^(d_r-j_r)` whenever j>=i. The integer verifier
checks every coefficient, the extracted factors and the determinant
identities. All basis factors and the extracted factors are nonnegative;
D_1,D_2,D_3 are positive. This proves every full-versus-pair increment.
The saved general one-common-particle theorem handles the other binary
comparisons. Surplus IDLA couples arbitrary inputs to nested binary inputs,
and self-loop removal preserves their occupation laws. These two reductions
complete the stated directed-three-cycle theorem.

There is also a finite padded inventory conclusion. Multiply the increment
by `h_0^4 h_1^4 h_2^4` and substitute K^X. The certificate's total jump
row degrees, including the extracted cycle, are (3,4,4), while its total
sleep degrees are at most (4,4,4). The extra h_0 supplies the potentially
missing target departure. All negative h-powers cancel, leaving positive
walk tails. The reciprocals of D_1,D_2,D_3 have positive closed-cycle
expansions: for D_3 use the sum over all coin assignments other than three
sleeps; the two terms in 1-D_2 specify respectively one or two departures
and at most that many coins at each site. Pad all unused coins positively.
Thus nine auxiliary walks, three from each site in addition to the three
ordinary labels, suffice for a nonnegative inventory expansion for each x.
This argument alone uses padding; the following certificate removes it on
the directed three-cycle. General c_1>=0 remains open.

## 11. Removing padding on the directed three-cycle

**Theorem, proved by a finite exact certificate:** on the directed
three-cycle, c_1(J,k) and c(J,k) are nonnegative for every feasible inventory,
every coin-count vector, every finite all-active input and each target.
Self-loops are allowed. The theorem concerns this graph geometry; it does
not cover general three-site kernels or longer directed cycles.

Use forward weights r_i, sink weights d_i, t=r_0 r_1 r_2, and

`h_i=c_i/(1-t)`, `c_i=d_i+r_i d_(i+1)+r_i r_(i+1)d_(i+2)`.

For a binary input B let `Z_B=H_B Q_B^(K^X)`. Multiplying the
complete-block inverse by H_B cancels the Doob factors on every forced
jump and every surplus-IDLA path. Thus its entries are the same path
expressions with forward weights r_i and sink weights d_i, without
negative powers of h. The single-particle Z-vector is
`p_0 h_0(1,q_1 r_1 q_2 r_2,q_2 r_2)/D_1`.
The pair equations have the matrix I-T of (14), with a replaced by r,
and right-hand side

`p_i p_j h_i h_j 1{0 in {i,j}}+q_j d_j Z_i`
` +q_i p_j r_i d_j Z_j+q_i p_j d_i h_j 1{j=0}`, j=i+1.            (16)

The full equation is (15) with each F_C replaced by the sum of the actual
X-weights of its forced-jump/IDLA histories, and right-hand side
`p_0 p_1 p_2 h_0 h_1 h_2`. The denominators D_1,D_2,D_3 are unchanged,
since all their edge products are closed cycles.

For every nonempty binary common input A and x outside A, set N=|A|+1.
There are six comparisons with N=2 and three with N=3. Cramer's rule in
(16) and the full equation give

`C_(A;x)=Z_(A union {x})-h_x Z_A`
`       = P_(A;x)(r,d,p)/[(1-t)^N D_1 ... D_N]`.                  (17)

The numerator P is a finite integer polynomial, homogeneous of degree N
in the sink weights. Each term has edge divergence 1_A+delta_x.
For a fixed vector e of sink counts, this forces all possible forward
counts to differ by a common cycle count. Consequently write

`P_(A;x)=sum_e d^e r^b P_e(t,p)`,

where b is the componentwise minimum forward exponent in that sector.
There are six choices of e when N=2 and ten when N=3. The verifier
constructs these polynomials directly from (16)--(17), restores the third
sink exponent by homogeneity, and checks every divergence identity.
No fitting from sampled inventories is involved.

Write `P_e(t,p)=sum_(j=0)^L v_j(p)t^j`, with L<=5, and put

`G_ell(p)=sum_(j<=ell) binomial(ell-j+N-1,N-1) v_j(p)`.

These are the coefficients after division by (1-t)^N. Their required
sitewise coin degrees are `m_i(ell)=b_i+e_i+ell`. For ell>=L, G_ell is
a polynomial of degree at most N-1 in ell. Therefore, for u>=0,

`G_(L+u)=sum_(j=0)^(N-1) binomial(u,j) Delta^j G_L`,               (18)

where Delta is the forward difference in ell. The certificate checks the
Bernstein coefficients of each Delta^j G_L at the **fixed** degrees m(L).
All are nonnegative. Elevating by u degrees at each site is positive, so
(18) certifies every subsequent coefficient, without truncating the number
of loops.

For N=2, and for N=3 with x!=0, the same certificate checks every G_ell
for 0<=ell<=L at its own degrees m(ell). For N=3,x=0, a few coefficients
of G_1 are negative. Include the full-input renewal denominator D_3 first:
set H_(-1)=0 and

`H_ell=G_ell+(1-p_0 p_1 p_2)H_(ell-1)`.

All H_ell for 0<=ell<=L have nonnegative Bernstein coefficients. The
tail result (18) then extends H positivity by induction. The multiplier
`1-p_0 p_1 p_2` is the sum of the seven coin assignments other than three
sleeps, so it uses exactly one coin per site and preserves the required
degrees. This is a cancellation across full-input renewals, not positivity
of the uncancelled numerator.

Finally, the remaining reciprocals D_j expand positively. For D_1 the
cycle weight is q_0 q_1 q_2. For D_2 the one-cycle weight is the sum of
assignments with at most one sleep, and the two-cycle weight is
`p_0 p_1 p_2 q_0 q_1 q_2`. For D_3 use the seven assignments just described.
Each term uses exactly as many coins per site as new cycle departures.
Multiplying these positive expansions by the certified G or H coefficients
proves nonnegativity of every coefficient in (17).

The certificate consists of **413 Bernstein arrays across 66 comparison/
exit-pattern sectors: 36,765 positive integer coefficients, 29,256 zeros,
and none negative**. Their canonical table has SHA256
`efa3de81fda2e90791804102d0e474198398cb00610bcf2ce047e89be3a8d11a`.
The script [verify_cycle_inventory.py](verify_cycle_inventory.py) reconstructs
the entire certificate using integer arithmetic. Its independent forest/
buffer enumeration agrees on 198 complete pools, with 23,802 visited states.
The finite-difference argument, together with the coefficient certificate,
proves all loop counts; the 198 checks separately audit the generating
identity.

To extend from binary inputs to arbitrary eta, run the common surplus-IDLA
phase first and the extra particle last. Let gamma be this joint prefix,
C its common settled cluster and y the extra's new site. Telescoping gives

`H_eta h_x Delta Q_eta = sum_gamma X^E_gamma C_(C;y)`.

An extra which exits contributes zero; an empty C gives a nonnegative
single-particle count. The other terms are the nine binary comparisons
above. Every prefix departure has an ignored sleep coin, so multiply by
positive powers of p_i+q_i to restore its prescribed coin degree. This
proves c_1 positivity for all multiplicities. Multiplying by H_eta proves
c positivity, and (9) inserts arbitrary self-loops. This finishes the proof.

## 12. Positive count-block renewals on longer directed cycles

The following identities are **proved for every cycle length n>=2**.
They concern the part of the complete-block inverse which preserves the
number of particles. They do not settle the signed forcing from smaller
particle counts, or the add-one increment on longer cycles.

Let r_i be the forward edge weight, p_i the sleep probability, q_i=1-p_i,
`t=product_i r_i` and `P=product_i p_i`. An initial forced-jump/IDLA history
which retains its particle count has no sink jump. Thus its homogeneous
count block depends only on r,p. Write F for this unnormalized inverse
block: its row for A is the signed sum over subsets S of initially forced
jumps, with weight `(-1)^|S| q_S`. The usual row factor 1/p_A is omitted.

**One vacancy.** Index inputs by their missing site k. If a nonempty S is
forced, let i be its first site in cyclic order after k. All forced jumps
and surplus moves traverse precisely the arc from i to k, leaving vacancy
i. Summing over the other selected sites in that arc gives

`F=I-T`,
`T_(k,i)=q_i product_(z in [i,k)) r_z product_(z in (i,k)) p_z`, i!=k,
`T_(k,k)=0`.                                                     (19)

The intervals run forward around the cycle, without a full revolution.
At r_i<=1 the row sum is at most `1-product_(i!=k) p_i<1`.
In particular F has a positive inverse by the geometric series. Every
factor in T specifies at most one coin per arc departure, so this is also
a positive expansion with the required edge and coin counts.

There is a useful closed determinant. Let E_j be the probability of at
least j zero bits among independent site coins with probabilities p_i
of one. Equivalently E_j is the sum of their degree-one coin monomials
with at least j q-factors. Then

`D_(n-1)=det F=1-sum_(ell=1)^(n-1) P^(ell-1) E_(ell+1) t^ell`.    (20)

To prove this, put `W_(i+1,i)=p_i r_i`, let R have entries
`R_(i+1,i)=r_i`, and set `b_i=t q_i product_(z!=i) p_z`.
Direct telescoping of (19) gives

`(I-W)F=diag(1+b_i)-R`.

Hence `det F=[product_i(1+b_i)-t]/(1-Pt)`. The coefficient of t^ell
in the product is P^(ell-1) times the probability of exactly ell zero
bits. Division by 1-Pt gives (20). This also explains why every term
subtracted in (20) has exactly ell coins and ell departures at each site.

For completeness, writing D for (20), the inverse entries are

`(F^(-1))_(k,i)=q_i product_(z in [i,k)) r_z`
`                 * product_(z outside [i,k]) (1+b_z)/D`, i!=k,
`(F^(-1))_(i,i)=1+q_i sum_(ell=1)^(n-1)`
`                 P^(ell-1) E_ell^(-i) t^ell/D`.                  (21)

Here [i,k] includes both endpoints, and E_ell^(-i) counts zero bits among
the sites other than i. Inverting the diagonal-minus-cycle matrix above
gives these formulas directly. All numerators and the reciprocal of D
have positive edge/coin expansions. Formula (20) specializes to D_2 in
Section 10 when n=3.

**Two particles.** Choose the cut order 0<...<n-1, and index states by
unordered pairs {i,j} with i<j. Define a unit triangular matrix L by

`L_({i,j},{k,j})=product_(z=i)^(k-1) q_z r_z`, i<=k<j,

and zero elsewhere. Then `LF=I-T'`, with the two contributions in row
{i,j}

`q_j r_j` to the pair {i,j+1},
`p_j r_j product_(z=i)^(j-1) q_z r_z` to the pair {j,j+1}.          (22)

Successors are modulo n. If the first pair has a repeated site, replace
`{i,i}` by `{i,i+1}` and multiply its weight by r_i, the forced surplus
jump. Contributions to the same state are added, including self-loops.
Every row sum is at most `q_j+p_j product_(z=i)^(j-1)q_z<1`.
Thus `F^(-1)=(I-T')^(-1)L` has a positive path/coin expansion for arbitrary n.

Here is the cancellation behind (22). In the homogeneous IDLA algebra,
let u_i denote a particle at i and use `u_i^2=r_i u_i u_(i+1)`.
The row F_{k,j} represents
`(u_k-q_k r_k u_(k+1))(u_j-q_j r_j u_(j+1))`.
Its weighted sum over k=i,...,j-1 telescopes to

`[u_i-(product_(z=i)^(j-1)q_z r_z)u_j]`
`             *[u_j-q_j r_j u_(j+1)]`.

Reducing u_j^2 gives exactly (22). This is a positive row elimination
of simultaneous-jump terms, not a new legal ARW scheduling assumption.
For four sites, this and the one-vacancy formula handle every nontrivial
homogeneous count block. Their lower-count forcing still requires a sign
argument; neither formula proves the four-cycle increment by itself.

One explicit four-cycle denominator is useful for the next calculation.
Let Q=product_i q_i, let A be the probability of at least three zero bits,
let `B=Q Pr(at most one zero bit)` and `C=P Q^2`. Direct determinant
expansion of the six-state two-particle block gives

`D_2=1-A t-B t^2+C t^3`.                                         (23)

The exact polynomial identity is checked by the four-cycle generator below.
Although (23) has a positive nonconstant coefficient, its reciprocal has
a positive expansion:

`1/D_2 = [1/(1-A t)]`
`          * [1-B t^2-(AB-C)t^3/(1-A t)]^(-1)`.

Indeed A contains Q and B contains P Q as individual positive coin
monomials, so AB-C is coefficientwise nonnegative. A, B and AB-C have
respectively one, two and three coins per site, matching their cycle
counts. Thus every count-block denominator needed for the four-cycle
calculation has a positive expansion with the correct coin degrees.

The exact rational verifier [verify_cycle_renewal.py](verify_cycle_renewal.py)
independently enumerates the signed forced-jump/IDLA blocks for cycle
lengths 2 through 8. It checks (19)--(22), including 203 explicit inverse
entries and 1,596 two-particle matrix entries. These audit the formulas;
the telescoping arguments prove them for all n.


## Verification scope

The companion retains the named inventory verifiers. Their exact checks have the scope stated in the reader and program headers; this publication does not rerun every optional mode.
