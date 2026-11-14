# Technical proofs for stationary ARW negative dependence

Selected mathematical sections supporting the November 14, 2026 reader article.
These arguments were generated with ChatGPT from research questions posed by
Jacob Richey. References to an authoring or reviewing AI describe AI checks;
independent human verification and priority review remain outstanding.
The general stationary negative-association conjecture remains open.

The excerpts preserve their equation labels. They cover the small- and large-rate
arguments, harmonic and complete-graph identities, full occupation, weak killing,
weak internal motion, and the common-entrance gluing theorem. Unrelated research
branches and operational continuation notes are omitted. Code and certificates
in this directory support the finite calculations; they do not turn the general
analytic proofs into computer-verified theorems.

## Small sleep rates on every fixed finite graph

**Theorem (complete working proof).** Let `K` be a substochastic transition
matrix on a fixed finite vertex set `V`, with an accessible absorbing sink
from every vertex. Fix positive numbers `alpha_v`, and give vertex `v`
sleep rate `lambda_v=t alpha_v`. There is an `epsilon(K,alpha)>0` such that
the stationary stable law `mu_t` is NA for every `0<t<epsilon(K,alpha)`.
In particular this proves the main conjecture for sufficiently small common
sleep rate on each fixed finite undirected sink graph. The threshold is
not uniform in graph size. This scoped small-rate theorem also allows
nonreversible kernels; it does not assert the withdrawn arbitrary-rate
general-kernel conjecture.

Let `G=(I-K)^(-1)` be the killed-walk Green matrix, so `G(x,y)` is the expected
number of visits to `y` before absorption by a walk starting at `x`, counting
time zero. Write `G_A` for its principal submatrix on `A`, with empty
determinant one, and `alpha_A=product_{v in A} alpha_v`. We first prove

\[
 \mu_t(X\supseteq A)
   =t^{|A|}\alpha_A\det G_A+O(t^{|A|+1}).                 \tag{1}
\]

All constants in this section may depend on the fixed graph and `alpha`.

**Rare-sleep bound.** Uniformly over configurations with at most `|V|+1`
particles, the probability of at least `j` successful sleep events during
one stabilization is `O(t^j)`. To see this, use a fixed legal priority order
and skip ineffective sleep instructions at sites with surplus particles.
At rate zero the resulting finite chain stabilizes with uniformly finite
expected number of topplings. One explicit reason is to prescribe, at each
toppling, a jump decreasing distance to the sink: the sum of particle
distances decreases at each such jump, so stabilization has a uniformly
positive chance in a bounded block of topplings. This also applies with
initial sleepers; they are either left asleep or awakened and routed out.

Until the first successful sleep, the process follows this zero-rate chain
with a sleep coin inserted whenever the toppled site has one active
particle. Each coin has probability
`t alpha_v/(1+t alpha_v) <= t max_w alpha_w`. A union bound over the expected
zero-rate number of topplings gives probability at most `C t` for the first
sleep. Restart after each sleep with the same uniform bound. Thus the
probability of at least `j` sleeps is at most `(C t)^j`.

The saturated sampler starts with all particles active. Every final sleeper
requires a successful sleep, so

`mu_t(X=B)=O(t^|B|)`                                              (2)

for every exact sleeper set `B`. In particular, replacing an exact set `A`
by the event `X superset A` changes its probability by `O(t^(|A|+1))`.

**Coefficient induction.** Fix nonempty `A`, choose `v in A`, set
`B=A\{v}`, and put `k=|A|`. Apply stationarity under addition at `v`, then
partition the initial stationary configuration by its exact sleeper set
`C`. The only contribution of order `t^k` can come from `C=B`:

- If `|C|<=k-2`, mass conservation prevents the final state from containing
  `A`, since only one particle is added.
- If `|C|=k-1` and `C!=B`, at least two successful sleeps are necessary.
  If `v` was asleep, addition wakes it, and a missing site of `A` must also
  sleep. If `v` was absent, at least two sites of `A` were absent. By (2)
  and the rare-sleep bound these contributions are `O(t^(k+1))`.
- If `|C|=k`, at least one successful sleep is necessary: either a site of
  `A` was absent or `C=A` and addition wakes `v`. Again the contribution is
  `O(t^(k+1))`. Configurations with more than `k` sleepers already have
  probability `O(t^(k+1))`.

For `C=B`, paths with exactly one successful sleep must avoid `B` and sleep
at `v`. A hit on `B` wakes a sleeper and would require at least two sleeps
to finish with all of `A` occupied. Before the first hit on `B`, only the
added particle moves. Hence, writing `G^B` for the Green matrix of the walk
additionally killed on hitting `B`,

\[
 \Pr_B(X_{\rm final}\supseteq A\mid\hbox{addition at }v)
       =t\alpha_vG^B(v,v)+O(t^2).                       \tag{3}
\]

For the error in (3), insert the no-sleep factors along the walk up to its
hitting time of `B` or the sink. Replacing their product by one costs at
most a constant times `t` times the path length. This hitting time has a
finite second moment, so multiplication by the terminal sleep probability
gives `O(t^2)`. Paths with two or more sleeps have the same error bound by
the rare-sleep lemma.

Equations (2)–(3), starting with `mu_t(empty)=1+O(t)`, give recursively

`c_A=alpha_v G^B(v,v)c_B`,

where `mu_t(X=A)=t^|A| c_A+O(t^(|A|+1))`. The taboo Green identity and Schur
complement formula say

\[
 G^B(v,v)=G(v,v)-G(v,B)G_B^{-1}G(B,v)
         =\frac{\det G_A}{\det G_B}.
\]

All these principal minors are positive: Jacobi's identity expresses
`det G_B` as `det(I-K)_(V\B)/det(I-K)`, and every principal restriction of
`K` has spectral radius less than one. This proves (1), including for
nonreversible kernels.

**Strict determinant inequality.** Suppose first that `K` is irreducible
on `V`. For disjoint nonempty `I,J subset V`,

\[
 \det G_{I\cup J}<\det G_I\det G_J.                    \tag{4}
\]

For reversible `K` this is strict Hadamard–Fischer after diagonal
symmetrization. The following direct proof also treats nonreversible `K`.
For `D subset V` set `h(D)=log det(I-K_D)`, where `K_D` is the principal
restriction. Since its spectral radius is less than one,

\[
 h(D)=-\sum_{m\ge1}\frac{\operatorname{tr}(K_D^m)}m.
\]

Each trace is the sum of nonnegative weights of rooted closed walks of
length `m` lying in `D`. Jacobi's identity gives

\[
 \log\frac{\det G_I\det G_J}{\det G_{I\cup J}}
 =h(V\setminus I)+h(V\setminus J)
   -h(V\setminus(I\cup J))-h(V).
\]

In this expression a closed walk contributes its weight divided by its
length exactly when it visits both `I` and `J`, and contributes zero
otherwise. Irreducibility supplies such a positive-weight closed walk,
proving strict positivity and (4). All series converge absolutely.

**From (1) to full NA.** Let `E,F` be nonconstant increasing events depending
on disjoint nonempty coordinate sets `U,W`. Let `a` be the smallest size
of a subset of `U` on which `E` occurs, and let `mathcal E` be the family
of such subsets of size `a`. Define `b,mathcal F` similarly. Both `a,b`
are positive. Formula (1) and the exact-set version of it give

\[
\begin{split}
 \mu_t(E)&=t^a\sum_{I\in\mathcal E}\alpha_I\det G_I+O(t^{a+1}),\\
 \mu_t(F)&=t^b\sum_{J\in\mathcal F}\alpha_J\det G_J+O(t^{b+1}),\\
 \mu_t(E\cap F)&=t^{a+b}
   \sum_{I\in\mathcal E,J\in\mathcal F}
       \alpha_I\alpha_J\det G_{I\cup J}+O(t^{a+b+1}).
\end{split}
\]

Each summand in the leading coefficient of
`mu_t(E)mu_t(F)-mu_t(E intersection F)` is strictly positive by (4).
There are finitely many such event pairs, so one common positive `epsilon`
makes every inequality hold. Every increasing real function on a finite
cube is a constant plus a nonnegative linear combination of increasing
event indicators. Thus these event inequalities imply NA for all the
functions in the conjecture.

For reducible `K`, use the exact factorization over strongly connected
components proved below. Apply the preceding argument to each component
with exits killed; singleton components are automatically NA. Take the
minimum of their thresholds. Independent products preserve NA, as follows
by conditioning on one component at a time. This completes the theorem.

For distinct vertices the expansion also gives the explicit diagnostic

\[
 \operatorname{Cov}_{\mu_t}(X_i,X_j)
   =-t^2\alpha_i\alpha_jG(i,j)G(j,i)+O(t^3).
\]

The three-site path symbolic calculation checks (1) for every nonempty
subset: singleton Green minors are `3/2,2,3/2`, each two-site minor is `2`,
and the full determinant is `2`. The computation supports the derivation;
the theorem uses the proof above, not extrapolation from finite tests.

**The small-rate argument also proves BK on irreducible kernels.** For
increasing events `E,F`, let `E square F` mean that there are disjoint
occupied subsets witnessing the two events. The BK inequality is
`mu(E square F)<=mu(E)mu(F)`. It implies NA when the events depend on
disjoint coordinate sets; it does not require conditional NA. These
distinctions are illustrated by
[Markström, Definition 1.2 and Theorem 3.4](https://emis.de/ft/41869).

Assume `K` is irreducible and `lambda_v=t alpha_v`, with fixed positive
`alpha_v`. Exclude constant events, which give trivial inequalities.
Let `a,b` be the minimum witness cardinalities for `E,F`, and write
`E_a,F_b` for their witnesses of those sizes. Put
`w(I)=alpha_I det(G_I)`. The product has leading coefficient

`mu_t(E)mu_t(F)=t^(a+b) sum_(I in E_a,J in F_b) w(I)w(J)
                 +O(t^(a+b+1))`.

Every witness for `E square F` has size at least `a+b`. Its witnesses of
exactly that size are the **distinct unions** `C=I union J` of disjoint
`I in E_a,J in F_b`; call this set of unions `U`. Thus

`mu_t(E square F)=t^(a+b) sum_(C in U) w(C)+O(t^(a+b+1))`.

For each `C in U`, choose one disjoint decomposition. Its product term
is strictly greater than `w(C)` by the irreducible disjoint Green-minor
inequality. Different unions select different ordered pairs, and every
unselected term in the first sum is nonnegative. If `U` is empty, the
first sum is still positive. Hence the leading BK gap is strictly
positive in every nonconstant case. There are finitely many pairs of
increasing events, so one sufficiently small interval works for all of
them. This proves BK for the stationary **sleeper indicators** in that
interval. No assertion of BK closure under products or of arbitrary-rate
BK is needed for this irreducible-kernel statement.

**A version allowing zero leading gaps.** More generally, suppose sparse
Boolean variables `Z_t` have cylinder expansions

`Pr(Z_t superset I)=t^|I| c_I(1+t d_I+O(t^2))`, `c_I>0`.

Assume: (i) `c_I c_J>=c_(I union J)` for disjoint nonempty sets;
(ii) whenever equality holds, `d_I+d_J-d_(I union J)>0`;
(iii) equality is closed under unions on either side, provided the two
sides remain disjoint. Then `Z_t` satisfies BK for all sufficiently small
`t>0`.

Here is the additional argument beyond the preceding leading-order proof.
If the leading BK gap for events `E,F` vanishes, there can be no overlapping
pair of minimum witnesses, no repeated union of disjoint pairs, and no
strict leading cylinder inequality. Hence the unions `U,W` of their
minimum witnesses are disjoint, and each pair `I in E_a,J in F_b` has
`c_I c_J=c_(I union J)`. At order `t^(a+b+1)`, these pairs contribute
the strictly positive sum

`sum_(I in E_a,J in F_b) c_I c_J(d_I+d_J-d_(I union J))`.

The remaining contributions cannot make this sum negative. For a set `T`
of size `a+1`, the next event coefficient is `1_E(T)-h_E(T)`, where
`h_E(T)` counts size-a witnesses contained in `T`. A negative coefficient
means `h_E(T)>=2`; then `T` is their union and lies in `U`. In the
disjoint-occurrence event its correction pairs with each fixed `J in F_b`
in exactly the same way: two leading unions whose union has size
`a+b+1` can differ on only one side, since `U,W` are disjoint. The
negative terms thus multiply leading cylinder gaps that vanish by (iii).
The same reasoning applies on the other side.

Every positive next event coefficient is a new minimal witness `T` of
size `a+1`, or a new size-`b+1` witness on the other side. If `T` overlaps
the opposite minimum witness, its product contribution is nonnegative
and creates no disjoint witness. Otherwise their union is a new witness
of size `a+b+1`: it cannot contain a leading disjoint union, since that
would force `T` to contain a minimum witness on its own side. Group equal
new unions and choose one pair for each. Assumption (i) bounds that
union's leading mass by the selected product; any extra pairs contribute
nonnegative terms. This accounts for every next-order term and proves
strict positivity. Finitely many event pairs give the claimed interval.

**Complementation for sparse BK families (complete working proof).**
Let `nu_t` be BK measures on a fixed Boolean lattice for all sufficiently
small `t>0`. Suppose every exact-state mass is positive and

`nu_t(S)/nu_t(T) -> 0` whenever `T` is a proper subset of `S`.

Then the complemented measures are also BK for all sufficiently small
`t>0`. In particular, the cylinder expansions above imply the required
mass ratios, since `nu_t(S)=t^|S| c_S+O(t^(|S|+1))`, with `c_S>0`.
This is a local statement near the empty configuration, not a claim that
arbitrary BK measures are closed under complementation.

First note that if increasing events `A,B` satisfy
`A intersect B=A square B` pointwise, they depend on disjoint coordinate
sets. Indeed, under the uniform product measure, product BK bounds the
intersection probability above by the product of its marginals, whereas
Harris bounds it below. If coordinate `i` is essential to both events,
conditioning on `i` makes the second term in the covariance decomposition
equal to `Delta_i Pr(A) Delta_i Pr(B)/4>0`; the first term is nonnegative
by Harris on the remaining coordinates. Equality is therefore impossible.
The product BK theorem is the classical input, stated as
[van den Berg–Jonasson, Theorem 1.1](https://arxiv.org/pdf/1105.3862).

Now fix increasing events `E,F` for `X=1-Y`, and define increasing events
`G={Y:1-Y notin E}`, `H={Y:1-Y notin F}`. Put

`Q={Y:1-Y in (E intersect F)\(E square F)}`,
`D=(G intersect H)\(G square H)`.

Every state `S in D` strictly contains a state `T in Q`. To prove this,
suppose there is no such `T subset S`. Restrict `E,F` to coordinates in
`S`, fixing all coordinates outside `S` to one. At every restricted
configuration where both events hold, the corresponding vacancy set
is a subset of `S`, so the absence of `Q` gives disjoint global witnesses.
Their intersections with `S` are disjoint restricted witnesses. The
restricted events therefore satisfy intersection equals disjoint
occurrence, and the preceding observation makes their essential supports
disjoint. Their complements under bit complementation, namely the
restrictions of `G,H` with coordinates outside `S` set to zero, also
have disjoint supports. Since `S in G intersect H`, these provide disjoint
witnesses for `G,H` contained in `S`, contrary to `S in D`. Finally,
`S notin Q`, because `G intersect H` means that neither `E` nor `F`
holds. The containment is strict.

Choose one such `T(S)` for each `S in D`. Finiteness and the mass-ratio
assumption give `nu_t(D)=o(nu_t(Q))` whenever `D` is nonempty. If `Q` is
empty, then `D` is empty as well. Writing all probabilities on the right
with respect to `Y`, the desired BK gap for `X` is exactly

`Pr(E)Pr(F)-Pr(E square F)`
` =Pr(G)Pr(H)-Pr(G square H)+Pr(Q)-Pr(D)`.

The first difference is nonnegative by BK for `Y`. The last difference
is eventually positive if `Q` is nonempty and is zero otherwise.
Finitely many event pairs give one common interval. This proves the claim.
Consequently, on every fixed irreducible killed kernel, **both sleeper
and vacancy indicators are BK at sufficiently small sleep rates**.
`bk.py` also checks the strict-containment statement directly on every
four-site event pair; this validates its encoding, not the all-size proof.

**Finite-rate BK diagnostic.** On the five-site tree with edges
`sink--0--1`, `1--2`, `1--3--4`, `bk.py` checks all 7,579 nonconstant
increasing events and all 28,724,410 unordered pairs, including equal and
overlapping events. Every BK inequality passes at `lambda=1,64`, separately
for sleeper and vacancy indicators. Thus the positive conditional
correlation on this graph at `lambda=64` does not obstruct BK there.
These are four finite rational certificates, not an arbitrary-rate proof.
The truth-table operation is checked against direct disjoint-witness
enumeration on all three-site event pairs. It also detects the gap
`-1/2500` in [Markström's Example 2.1](https://emis.de/ft/41869), whose law
is NA but not BK. The law, comparisons, and validation are saved in
`bk_certificates.json`.

The two-arm tree `sink--0--1`, `1--2--3`, `1--4--5` gives a further
**all-rate NA certificate on six non-sink sites**. The purpose of this
case was to test whether two pendant arms could amplify conditional
dependence into a disjoint-event violation. At `lambda=1,64`, all 58,674
nonconstant increasing-event pairs over unordered bipartitions pass.
Every disjoint-support pair occurs in this list after unused coordinates
are assigned to one side. A symbolic check then proves the same result
for every common `lambda>0`: all gap numerators have nonnegative
coefficients in `t=1/lambda`, and none is identically zero.

For this certificate, exchanging the two arms reduces 64 stable masks to
40 orbits. The basis matrix is triangular by particle count; solving its
successive diagonal blocks by fraction-free elimination produces a
common positive-coefficient denominator of degree 112. The largest
block has dimension 12. Lifting the orbit probabilities back to all
64 states verifies `M mu=e_V` symbolically, as well as normalization;
specialization at rates 1 and 64 agrees with the separate rational
solver. Thus no unverified symmetry assumption enters the final law.
For each increasing event on the smaller side of a bipartition, the
checker forms a signed polynomial measure on the other side whose
event sums are precisely the NA numerators. It then checks every upset.
The checker also accepts a product law and detects the exact two-site
insertion counterexample. The law polynomials and complete comparison
counts are in `bk_certificates.json` under `balanced_tree_NA_symbolic`.
This certifies the specified six-site graph, not arbitrary trees.

Two further **all-rate finite certificates** pass. For the original
three-site Rayleigh counterexample (`0--1--2`, sink edges at `0,1`) and
the four-site endpoint path, every BK-gap numerator has nonnegative
coefficients in `t=1/lambda` after clearing a common positive polynomial
denominator. Both sleeper and vacancy orientations are checked, with
171 and 13,861 event pairs respectively. The exact denominators, law
numerators, and complete coefficient-test counts are saved in the same
certificate. No claim for every four-site graph is made here.

**An obstruction to separating witnesses between two copies (proved).**
For increasing events `E,F`, let `D(E,F)` be the set of pairs `(x,y)`
admitting disjoint occupied witnesses `I subset x`, `J subset y` for
the respective events. The strengthened inequality
`mu(E square F)<=(mu tensor mu)(D(E,F))` is false for stationary ARW.
It must not be substituted for ordinary BK.

Take distinct `a,b,c`, with `E={X_a=1 or X_c=1}` and
`F={X_b=1 or X_c=1}`. Write `v_I=Pr(X_i=0 for all i in I)`.
The one-copy disjoint occurrence is the event that at least two of these
three sites are occupied. In `E times F`, the pairs lacking disjoint
witnesses are exactly those with `x_a=0,x_c=1` and `y_b=0,y_c=1`.
Consequently, direct subtraction gives

`(mu tensor mu)(D(E,F))-mu(E square F)`
` =v_ab-v_a v_b+v_a v_bc+v_b v_ac-2v_abc`.

If `v_I=t^|I| c_I+O(t^(|I|+1))`, the leading term is
`t^2(c_ab-c_a c_b)`. Thus strictly negative leading vacancy covariance
defeats this strengthening. On the complete graph with three non-sink
vertices and one sink, common rate `lambda=2`, the masses of each exact
configuration of cardinality `0,1,2,3` are respectively
`(1051,3747,11592,30996)/78064`. The strengthened gap is
`-1803573/6093988096`, whereas the ordinary BK gap for these same events
is `58370337/1523497024>0`. Here `c_a=2/3`, `c_ab=32/81`, so the
strengthened gap is also `-4/(81 lambda^2)+O(lambda^-3)` at large rates.
`bk.py --two-copy-obstruction` checks the rational law by raw stabilization,
both toppling orders, every addition equation, and the independent IDLA
algebra. This refutes neither BK nor NA.

In particular, no probability-preserving switching of two independent
instruction fields can always send a one-copy disjoint occurrence into
`D(E,F)`. Any successful ordinary-BK switch must sometimes allow overlapping
witness supports across its two outputs. Fresh suffixes at an adapted
stopping line do not remove this law-level obstruction. Cutting just
before a terminal sleep also need not be adapted: on one site,
`final odometer minus terminal occupation` is zero exactly when the first
instruction is sleep, a fact not measurable from the empty prefix.

## Large sleep rates: escape probabilities and region volumes

**Theorem (complete working proof).** Suppose the killed kernel `K` is
irreducible on `V`, and, for every two distinct vertices `i,j`, the sink is
reachable from `i` without visiting `j` and from `j` without visiting `i`.
Fix positive `alpha_v`. The stationary law with `lambda_v=alpha_v/t` is NA
for all sufficiently small `t>0`. In particular, every fixed finite simple
undirected **biconnected** graph including its sink has an interval of large
common sleep rates on which NA holds. The threshold depends on the graph;
this does not settle intermediate rates or graphs with sink-separating
articulation vertices. Independent component factorization also applies
when the condition holds inside each component with its exits killed.

Write `Y_v=1-X_v` for vacancy, and retain `G=(I-K)^(-1)`. For nonempty
`A subset V`, let

`e_A(v)=Pr_v(T_sink < T_A^+)`, `v in A`,

where `T_A^+` is the first visit to `A` after at least one jump. The last
visit to `A` before absorption gives `G_A e_A=1`. In particular
`e_{ {v} }(v)=1/G(v,v)`. These escape probabilities are nonnegative,
not all zero on any nonempty `A`, and satisfy
`e_B(v)>=e_A(v)` for `v in B subset A`.

**Expansion for every killed kernel, without the extra escape condition.**
There are positive coefficients `C(A)`, depending only on `K`, such that

\[
 \Pr(Y\supseteq A)
  =\frac{t^{|A|}}{\prod_{v\in A}\alpha_v}C(A)
       +O(t^{|A|+1}).                                      \tag{6}
\]

They are specified by `C(empty)=1` and the triangular recurrence

\[
 \sum_{H\subseteq A}(-1)^{|A\setminus H|}C(H)
          \prod_{v\in A\setminus H}e_A(v)=0,
       \qquad A\ne\varnothing.                            \tag{7}
\]

Here is a justification of the expansion and recurrence. At `t=0` every
lone active particle sleeps; only surplus particles jump. With a fixed
legal toppling priority this finite absorbing chain has uniformly finite
mean absorption time for configurations with at most `|V|` particles.
Prescribe sink-directed jumps: each decreases total particle distance to
the sink, and between jumps there are at most `|V|` lone sleep topplings.
Thus absorption has uniformly positive probability in a bounded block.
The finite absorption equations are analytic in `t` near zero. Before
the first lone jump, the process is the zero-parameter chain with rare
jump coins of probability `t/(alpha_v+t)`. Its finite expected length and
restarting after each such jump bound the probability of `m` lone jumps
by `O(t^m)`. Every final vacancy needs a lone departure, so the probability
of an exact vacancy set `H` is `O(t^|H|)`. Exact sets and containing events
therefore have the same leading coefficient.

First assume direct killing is possible at every vertex. Apply the
harmonic identity (5) with `J=V\A`. For `v in J` its value `a_v` is the
strictly positive probability of reaching the sink before `A`; for
`v in A` its sleeper weight is `-t e_A(v)/alpha_v`. In the resulting zero
of `P`, a term with vacancy set `H` has order at least
`t^(|A|+|H\A|)`. At order `t^|A|` only `H subset A` contributes.
Divide by the positive product of the `a_v` on `J`, and by
`product_(v in A) alpha_v^(-1)`, to obtain (7).

For general `K`, replace it by `(1-epsilon)K`, adding direct killing
everywhere. The same finite absorption matrix is invertible at
`(t,epsilon)=(0,0)`, so its Taylor coefficients in `t` are continuous in
`epsilon`. The Green matrices and escape probabilities are also continuous.
Letting `epsilon` decrease to zero proves (6)–(7) without direct killing.

**Region interpretation and disjoint inequalities.** For each `A` define
the bounded region

\[
 \mathcal R_A=\{u\in[0,\infty)^A:
   \text{for every nonempty }B\subseteq A,
   \text{ some }v\in B\text{ has }u_v<e_B(v)\}.
\]

Then `C(A)` is its `|A|`-dimensional volume. Indeed, every point of
`mathcal R_A` has some coordinate below `e_A(v)`. For any nonempty
`D subset A`, imposing `u_v<e_A(v)` for all `v in D` automatically
satisfies every constraint meeting `D`, since `e_A(v)<=e_B(v)` whenever
`v in B subset A`. The remaining constraints are precisely those for
`mathcal R_(A\D)`. Inclusion–exclusion over these coordinate strips gives
(7). Also the region has positive volume: take all coordinates smaller
than the minimum, over nonempty `B subset A`, of `max_(v in B)e_B(v)>0`.

For disjoint `I,J`, restricting coordinates gives

`mathcal R_(I union J) subset mathcal R_I times mathcal R_J`,

and hence

\[
 C(I\cup J)\le C(I)C(J).                                  \tag{8}
\]

This proves nonnegativity of every leading disjoint vacancy-cylinder gap
on **every** killed kernel. Equality is possible, so (8) alone does not
prove NA at large finite rates on arbitrary graphs.

**Strictness under the theorem's hypotheses.** Write
`c_i=C({i})=1/G(i,i)`. For `A={i,j}`, (7) and `G_A e_A=1` give

\[
 \begin{split}
 c_i c_j-C(\{i,j\})
  &=(c_i-e_A(i))(c_j-e_A(j))\\
  &=\frac{G(i,j)G(j,i)e_A(i)e_A(j)}{G(i,i)G(j,j)}.
 \end{split}                                               \tag{9}
\]

Irreducibility and the pairwise escape assumption make this strictly
positive. To lift strictness from pairs to arbitrary disjoint nonempty
`I,J`, choose `i in I,j in J` and numbers
`e_{ {i,j} }(i)<u_i<c_i`, `e_{ {i,j} }(j)<u_j<c_j`.
Each can be extended to an interior point of its own region by making
all other coordinates sufficiently small and positive. To verify this
extension, consider any nonempty `B subset I`. If some `v in B\{i}`
has `e_B(v)>0`, a sufficiently small `u_v` satisfies the constraint.
If no such vertex exists, then `i in B` and the equation `G_B e_B=1`
in row `i` gives `e_B(i)=c_i>u_i`. The same argument handles `J`.
There are only finitely many constraints, so these choices give an open
box in `mathcal R_I times mathcal R_J`. It is excluded from
`mathcal R_(I union J)` by the constraint for `{i,j}`. Thus (8) is strict.

Finally expand increasing events in `Y` by their smallest triggering
sets, exactly as in the small-rate proof. Each leading gap is a positive
sum of the strict quantities in (8), multiplied by positive products of
`alpha_v^(-1)`. Finitely many event pairs give one common interval of
small `t`. Threshold decomposition gives NA for arbitrary increasing
functions. NA is invariant under complementing **all** coordinates:
negating two decreasing functions leaves their covariance unchanged.
Consequently `X` is NA too, proving the theorem.

**Checks and limitation.** `vacancy.py` verifies (6) against all eight
three-site symbolic stabilization laws and independently integrates the
regions over rational boxes. It also checks the duplicated/shared-gate
motifs described in its certificate. None violates (8), as the proof now
explains. On a path with the sink attached only at one endpoint all the
recorded disjoint leading gaps vanish, so higher orders are essential.
On the previous conditional-correlation example, the coefficients are

`C_0=3/4, C_1=1/2, C_2=3/5, C_01=1/3, C_02=77/180,`
`C_12=3/10, C_012=26/135`.

They satisfy `C_0 C_012-C_01 C_02=1/540>0`, recovering the leading
conditional-correlation failure. Region inclusion proves the **disjoint**
inequality (8); it does not imply log-submodularity for overlapping sets.

## Large sleep rates on endpoint paths

**Theorem (complete working proof).** On every fixed path
`sink--0--...--(n-1)`, the law with `lambda_i=alpha_i/t`, for fixed
positive `alpha_i`, is NA for all sufficiently small `t>0`. The proof also
allows positive edge conductances and self-loops, with killing only at
the left endpoint. It resolves the vanishing leading gaps in this family;
it gives neither an all-rate result nor a volume-uniform threshold.

Let edge `(i-1,i)` have conductance `g_i>0`, interpreting `-1` as the sink,
and let the self-loop weight at `i` be `ell_i>=0`. Put `g_n=0`,
`m_i=g_i+g_(i+1)+ell_i`; each transition probability is its weight divided
by `m_i`. Define

`R_i=sum_(j=0)^i 1/g_j`, `c_i=1/(m_i R_i)`, `w_i=c_i/alpha_i`.

Thus `R_i=i+1` and `m_i` is the degree in the simple-walk case. For
`F_A(t)=Pr(Y superset A)`, the expansion through one extra order is

`F_A(t)=t^|A| (product_(i in A) w_i)[1+t d(A)+O(t^2)]`,

`d(A)=sum_(i in A)[-w_i+R_i sum_(i<v<i^+) w_v/(R_v-R_i)]`,       (10)

where `i^+` is the next selected site, or `n` for the last site, and
`d(empty)=0`. In particular, with

`kappa_(ij)=R_i sum_(v=j)^(n-1) w_v/(R_v-R_i)`, `i<j`,

`d(A)=sum_(i in A)d({i})-sum_(i,j consecutive in A)kappa_(ij)`.
Every `kappa_(ij)` is positive and strictly decreases as `j` increases.

Here is a derivation of (10), including the degenerate harmonic points.
Regularize `K` to `K_epsilon=(1-epsilon)K`. Write

`F_A(t,epsilon)=t^|A|/alpha_A [C_A(epsilon)+t B_A(epsilon)+O(t^2)]`.

The finite absorption argument used for (6) makes these coefficients
analytic near `epsilon=0`. For positive epsilon all vertices can kill
directly. Let `a_A(v)=Pr_v(T_sink<T_A)` off `A`, and
`h_A(v)=(1-a_A(v))/a_A(v)`. All quantities in the next identity use the
regularized kernel. Expand the harmonic zero via
`P(b)=sum_H F_H product_H(1-b_v) product_(V\H)b_v`, with
`b_v=-t e_A(v)/alpha_v` on `A` and `b_v=a_A(v)` elsewhere. At order
`t^(|A|+1)`, after dividing out the outside `a_A` factors and `alpha_A^-1`,
this gives

```
sum_(H subset A) (-1)^|A\H| (product_(A\H) e_A(v))
  [B_H+C_H sum_(v in H)e_A(v)/alpha_v
       +sum_(v outside A)h_A(v) C_(H union {v})/alpha_v] = 0.    (11)
```

At epsilon zero, if `i=min A`, only `e_A(i)=c_i` is nonzero. Recurrence
(7) therefore gives `C_A=product_A c_i`. Primes below mean derivatives in
epsilon at zero. For `i<j`, set

`eta_(ij)=R_i E_j[T_i]/(R_j-R_i)`, and `eta_(i,n)=0`.

If `j` is the next site of `A`, then
`e_A(i)'/c_i=c_i'/c_i-eta_(ij)`. Indeed, the difference between the
singleton escape probability and `e_A(i)` consists of excursions hitting
`j` before returning to `i`, then being killed before hitting `i`.
The former probability at zero is `1/[m_i(R_j-R_i)]`; the latter has
derivative `E_j[T_i]`. Differentiate (7). Terms involving derivatives of
the other escape coordinates cancel because `C_A=product_A c_i`, giving

`C_A'/C_A=sum_(i in A)c_i'/c_i-sum_(i,j consecutive in A)eta_(ij)`. (12)

To take the apparently singular limit in (11), put

`S_A(v)=sum_(H subset A)(-1)^|A\H| C_(H union {v}) product_(A\H)e_A`.

Always `S_A(v)=0` at epsilon zero. For `v<i`, `h_A(v)` stays finite, so
its contribution vanishes. For `v>i`, let `j=i^+`. Formula (12) and the
same cancellation of the other escape derivatives give

`S_A(v)'/(C_A c_v)=eta_(ij)-eta_(i,min(j,v))`.

Also `a_A(v)=epsilon E_v[T_A]+O(epsilon^2)`. If `v>j`, the displayed
derivative is zero. If `i<v<j`, the strong Markov identity

`E_v[T_i]=E_v[T_{ {i,j} }]+[(R_v-R_i)/(R_j-R_i)]E_j[T_i]`

shows that
`[eta_(ij)-eta_(iv)]/E_v[T_A]=-R_i/(R_v-R_i)`.
For the terminal interval `j=n`, use `T_A=T_i` and `eta_(i,n)=0` instead.
Consequently (11) tends to

`B_A-c_i B_(A\{i})+C_A w_i-C_A R_i sum_(i<v<i^+)w_v/(R_v-R_i)=0`.

Divide by `C_A` and iterate over the selected sites. This proves (10).

For disjoint nonempty `I,J`, their cylinder gap has first possible term

`F_I F_J-F_(I union J)`
` =t^(|I|+|J|+1)(product_(I union J)w_i)`
`   *[d(I)+d(J)-d(I union J)]+O(t^(|I|+|J|+2))`.              (13)

The bracket is strictly positive. Color the ordered union by membership
in `I` or `J`. Same-color edges consecutive in the union cancel from its
expression in kappas. Each remaining same-color edge `(i,j)`, consecutive
only in its own color, is charged to the cross-color edge from `i` to its
next union site `v<j`. These charges are distinct, and
`kappa_(iv)>kappa_(ij)`. The last cross-color edge is uncharged and positive.
This proves the strict inequality for all disjoint blocks.

Finally, write any increasing Boolean event as a multilinear polynomial
in vacancy indicators. Its smallest nonzero monomials are precisely its
minimum-size triggering sets, each with coefficient one. Every cylinder
gap has a zero leading term and the strictly positive next term (13).
Thus the first event-gap coefficient is a positive sum over pairs of
minimum-size triggers; all other monomials enter later. Finitely many
event pairs give a common small-t interval, and threshold decomposition
and complementing all coordinates give NA for the sleeper indicators.

`vacancy.py` checks every moment and disjoint-cylinder first correction
for simple endpoint paths with two through six sites, a four-site path
with unequal `alpha`, and a weighted four-site path with self-loops and
unequal `alpha`. The three-site coefficients also agree independently
with the earlier symbolic raw-absorption law. These are checks on the
formula; the argument above supplies its all-size proof.

## Removing the large-rate articulation restriction

**Theorem (complete working proof).** Let the underlying finite sink graph
be undirected, with the sink accessible from every vertex. For fixed
positive `alpha_v`, the stationary law with `lambda_v=alpha_v/t` is NA
for all sufficiently small `t>0`. Positive edge weights and self-loops are
allowed. Thus the main undirected conjecture holds at large sleep rates
on every fixed graph, without a biconnectivity assumption. The threshold
may depend on the graph and `alpha`; intermediate rates remain open.

By component factorization, assume the non-sink graph is connected.
Say `v` dominates `u` if every path from `u` to the sink passes through
`v`. Put `D_v={u!=v: v dominates u}` and
`O_v=V\({v} union D_v)`. Write `C_A=C(A)`, `c_v=C_{ {v} }`, and define
the relative next coefficient by

`F_A(t)=t^|A| C_A/alpha_A [1+t d_A+O(t^2)]`, `d_empty=0`.

**Equality in the leading inequality.** For disjoint nonempty `I,J`,

`C_(I union J)=C_I C_J`

if and only if every pair `i in I,j in J` is comparable by domination.
If a cross-pair is incomparable, both pair escape probabilities in (9)
are positive. The same open-box argument used after (9) makes the
region inclusion strict for `I,J`. Conversely, a vertex dominating all
other selected vertices is the only nonzero escape coordinate, so (7)
gives `C_( {v} union W)=c_v C_W` whenever `W subset D_v`.
Ancestors of any vertex under domination form a chain: their order along
a path to the sink determines which separates the other from the sink.
If all cross-pairs are comparable, the union has a unique minimal vertex
until one block is empty. Otherwise two incomparable minimal vertices
would have the same color and would both dominate any opposite-colored
vertex, contradicting this chain property. Peel off this minimal vertex
and use the displayed factorization repeatedly. This proves equality.

**Correction at a separating vertex.** Fix nonempty `W subset D_v`.
For `u!=v`, define

`H_v(u)=Pr_v(T_u<T_v^+ and T_u<T_sink)`,

and, for `u outside W`, let `H_v(u;W)` be the same probability with the
additional requirement `T_u<T_W`. Put `r_W(u)=C_(W union {u})/C_W`
for `u outside W`. For `u in O_v`, also put

`h_v(u)=Pr_u(T_v<T_sink)/Pr_u(T_sink<T_v)`.

The denominator here is positive. Then

`d_{ {v} }+d_W-d_( {v} union W)=K(v,W)`,                       (14)

where

```
c_v K(v,W)
 = sum_(u in W) c_u H_v(u)/alpha_u
   +sum_(u in D_v\W) [c_u H_v(u)-r_W(u) H_v(u;W)]/alpha_u
   +sum_(u in O_v) h_v(u)[c_v-e_{ {v,u} }(v)]
                         [c_u-r_W(u)]/alpha_u.              (15)
```

All summands are nonnegative: (8) gives `r_W(u)<=c_u`, and
`H_v(u;W)<=H_v(u)`. The first sum is strictly positive. Set `K(v,empty)=0`.
There is also the precise monotonicity needed below: for disjoint
`W,J subset D_v` satisfying `C_(W union J)=C_W C_J`,

`K(v,W union J)>=K(v,W)`.                                    (16)

Indeed, for `u outside W union J`, (8) gives
`r_(W union J)(u)<=r_W(u)`, while the taboo hitting probability decreases.
For `u in J` the subtracted term in the second sum of (15) disappears.
Every other summand is unchanged or increases. This uses the stated
leading-factorization hypothesis; no conditional NA is asserted.

Here is a proof of (14)–(15). Use the regularization and identity (11),
with `B_A=C_A d_A`. At epsilon zero, the only escape coordinate on
`A={v} union W` is `e_A(v)=c_v`. The first two terms of (11) become

`B_A-c_v B_W+C_A c_v/alpha_v`.

For a site `u in O_v`, undirectedness ensures every path from `u` to
`W` passes through `v`. Thus `h_A(u)=h_v(u)`. In the recurrence for
`C_( {v} union W union {u})`, only `v,u` have nonzero escape coordinates,
equal to those for `{v,u}`. Consequently, writing
`q_v=e_{ {v,u} }(v)` and `q_u=e_{ {v,u} }(u)`,

`C_(vWu)=q_v C_(Wu)+q_u(c_v-q_v)C_W`.

Compare its contribution in (11) with the singleton equation at `v`.
Their difference is exactly the `O_v` term of (15).

For the apparently singular terms in `D_v`, let

`eta(v,W)=(c_v'-e_(vW)(v)')/c_v`,

where primes now denote epsilon derivatives at zero. By conditioning at
the first visit to `W`, this is `1/c_v` times the expected duration from
that visit until the next return to `v`, counted only when `W` is visited
before return or absorption. In particular

`eta(v,{u})=H_v(u) E_u[T_v]/c_v`,

`eta(v,W union {u})-eta(v,W)`
` =H_v(u;W) E_u[T_( {v} union W)]/c_v`.                      (17)

The second identity counts just the extra excursion time between the
first visit to `u` and the first subsequent visit to `{v} union W`.
Differentiate (7) for `{v} union Z`, where `Z subset D_v`. All derivatives
of the other escape coordinates cancel using `C_(vZ)=c_v C_Z`, leaving

`C_(vZ)'=e_(vZ)(v)' C_Z+c_v C_Z'`.

It follows that the derivative of `S_A(u)` from the proof of (10) is

`S_A(u)'=c_v C_(W union {u})[eta(v,W)-eta(v,W union {u})]`.

Also `a_A(u)=epsilon E_u[T_A]+O(epsilon^2)`. Divide the last two
expressions and use (17). The normalized contribution to `d_(vW)-d_W`
is `r_W(u) H_v(u;W)/(c_v alpha_u)`; the singleton contribution to `d_v`
is `c_u H_v(u)/(c_v alpha_u)`. Subtracting the block equation from the singleton
equation gives the `D_v\W` terms in (15); sites in `W` contribute only
to the singleton equation and give its first sum. This proves (14).

**All disjoint cylinders.** If their leading gap is zero, the comparability
criterion allows the same minimal-vertex peeling. For `v in I` minimal
in `I union J`, put `I'=I\{v}`. Then

`d_I+d_J-d_(I union J)`
` =d_(I')+d_J-d_(I' union J)+K(v,I' union J)-K(v,I')`.

The last difference is nonnegative by (16), since the remaining cross-pairs
are still comparable. Continue until one block becomes empty; that last
step contributes a strictly positive `K(v,W)`. Thus every zero leading
cylinder gap has a strictly positive coefficient at its next order.

**From cylinders to arbitrary increasing events.** Let the minimum
trigger sizes of two nonconstant increasing events on disjoint blocks
be `a,b`. If the sum of leading cylinder gaps over their minimum-size
triggers is positive, it already settles the event gap. Otherwise every
such pair has zero leading gap, and hence comparable cross-pairs.
At the next order their correction gaps are all strictly positive.
The only potentially negative extra terms come from event monomials of
size `a+1` or `b+1`. A coefficient at size `a+1` equals
`1_E(S)` minus the number of size-a triggers contained in `S`. If it is
negative, `S` is the union of at least two such triggers. Its cross-pairs
with each size-b trigger are therefore still comparable, so its leading
cylinder gap is zero. The same argument applies on the other block.
All other next-order contributions are nonnegative. The event gap is
therefore strictly positive at one of these two orders. Finitely many
event pairs, threshold decomposition, and complementing all coordinates
finish the proof of NA.

The same hypotheses also prove **BK for the vacancy indicators** at all
sufficiently large rates, when the non-sink component is connected.
Apply the sparse BK criterion above with `c_I=C(I)/product_I alpha_v`.
The disjoint inequality gives (i); the cylinder correction just proved
gives (ii); and the cross-pair comparability criterion gives (iii).
The exact vacancy masses have the required sparse orders, so the sparse
complementation lemma also gives **BK for sleeper indicators** in this
large-rate interval. Thus both orientations satisfy BK near both ends
of the rate range on every fixed connected undirected non-sink graph.
The argument does not establish BK at arbitrary rates or unrestricted
closure of BK under complementation or products.

**Checks.** `vacancy.py` compares (15) with the independent IDLA-algebra
series on a five-site rooted fork and a five-site graph with a triangle
below a cut vertex, using unequal `alpha` in the latter. It checks every
dominated subset, (16), every zero leading cylinder gap, and all disjoint
increasing-event gaps through order six. The general result depends on
the proof, not on extrapolating these cases. Independent series also check
all 13,861 vacancy BK event pairs on each of the four-site endpoint path
and the tree `sink--0--1`, `1--2`, `1--3`. Every gap is positive at one
of its first two possible orders, as required by the sparse BK criterion.

## Stationary sampling and the DAG case

**Established input, checked in the source.** Write `S` for stabilization,
`1_V` for one active particle at each non-sink vertex, and
`p_v = lambda_v/(1+lambda_v)`. The stationary law is `Law(S(1_V))`.
It is invariant under addition at every vertex. For independent additions
whose support can reach every vertex, it is the unique stationary law and
does not depend on the drive. This includes uniform additions. See
[Levine–Liang, Theorem 1 and Lemma 6](https://arxiv.org/pdf/2110.14008),
also available in `../papers/2110.14008v3-source.txt`.

For the invariance argument, start with `1_V` and any extra particles.
Move only surplus particles until all leave, keeping one active particle
at every vertex. These moves are legal; ineffective sleep instructions can
be consumed. The accessible sink ensures termination. The unused independent
stack suffixes then give a fresh stabilization of `1_V`. Abelianity gives
the same law as adding the extra particles after first stabilizing `1_V`.
This explains the exact sampler; it is not a new theorem of this project.

**Proved here (no novelty claim): DAG product law.** If the transition graph
on `V` is acyclic, then

`Law(S(1_V)) = product_{v in V} Bernoulli(p_v)`.

Indeed, process vertices in a topological order. When vertex `v` is processed,
it has its original active particle plus all arrivals from predecessors, so
its number of particles is at least one. Consume instructions until all
surplus particles have jumped. With one particle remaining, the next
instruction is a fresh sleep/jump trial: a sleep leaves occupation one,
and a jump leaves occupation zero. There are no self-loops or later returns.
Conditioned on everything exposed at earlier vertices, this bit is
Bernoulli(`p_v`), since the stack at `v` has not been inspected and its suffix
after the surplus jumps is fresh. Iterating these conditional laws gives
independence. The generating polynomial is a product of linear factors
`1-p_v+p_v z_v`, hence is real stable and the law is strongly Rayleigh and NA.

**Proved here: factorization over strongly connected components.** More
generally, let `C_1,...,C_m` be the strongly connected components of the
transition graph, ordered so that edges between components point forward.
For each `C_i`, collapse all jumps out of `C_i` to a sink and let `mu_i`
be its saturated stationary law. Then `mu = product_i mu_i`.
To prove this, stabilize entire components in that order. Conditional on
the exposed earlier components, the input to `C_i` is one active particle
per vertex plus a fixed nonnegative vector of extra arrivals. Its stacks
are unexposed and independent. The saturated sampling identity above says
that its final law is `mu_i`, independently of that vector. Later components
cannot send particles back. Iterated conditional expectations prove the
factorization. The DAG product law is the singleton, no-self-loop case.

These are direct stationary-law arguments, not a verification of the
historical add-a-sink gluing outline, which remains conditional in the
shared compendium.



## Dynamics preservation fails even for two sites

**Exact counterexample, distinct from stationarity.** On the directed chain
`0 -> 1 -> sink`, start with a sleeper at vertex 1 and add an active particle
at 0. The initial stable law is deterministic, hence strongly Rayleigh.
With sleep probabilities `p_0,p_1`, write `q_i=1-p_i`. The final law has
mass `p_0` on `{0,1}`, `q_0 p_1` on `{1}`, and `q_0 q_1` on the empty set.
Consequently

`Cov(X_0,X_1) = p_0 q_0 q_1 > 0`.

The same failure occurs for the initial undirected model: take edges
`{sink,0}`, `{0,1}`, `{1,sink}`, common `lambda=1`, the deterministic sleeper
set `{1}`, and addition at 0. In bitmask order the final masses are
`(25,19,136,210)/390`; covariance is `1333/76050 > 0`.
Even uniform drive fails: on the undirected graph its output masses from
the same initial sleeper set are `(25,19,58,54)/156`, with covariance
`31/3042 > 0`. On the directed chain at `p_0=p_1=1/2`, uniform drive gives
`(3,0,3,2)/8`, with covariance `3/32 > 0`.

Thus a single-site addition/stabilization operator does not preserve NA,
Rayleigh, or strong Rayleigh. This verifies a precise version of the
recollection; it does not identify its historical example. The stationary
two-site law `(7,10,10,12)/39` is strongly Rayleigh. Any preservation proof
needs a narrower class or a different evolution.

Even restricting to exchangeable starting laws on complete graphs does
not repair preservation. On `K_4` including the sink, with `lambda=1`,
start with all three non-sink sites asleep and make one uniform addition.
The rank probabilities are
`(3721/102660,12737/61596,31556/76995,151/435)`, giving pair covariance
`7366296671/853665123600>0`. Thus even this deterministic exchangeable
ULC law loses NA and ULC. `complete_graph.py` verifies the rank law by
raw absorption; its certificate records this separately from stationarity.

There is a stronger obstruction for the all-active stabilization kernel.
On `sink--0--1--2--3`, start with one active particle at each of `0,1`
and no particles at `2,3`. Set `t=1/lambda`. The exact expansion gives

`E[X_2]=3t/4+O(t^2)`, `E[X_3]=7t^2/16+O(t^3)`,
`E[X_2 X_3]=t^2/16+O(t^3)`.

Thus `Cov(X_2,X_3)=t^2/16+O(t^3)>0`. Even a deterministic binary
all-active input can lose NA. The exact `lambda=64` covariance is positive
and its entire stabilization column is checked by raw absorption in
`algebra_certificates.json`. This is a nonsaturated input; the stationary
sampler starts with one active at **every** vertex.

## Two-site stationary law for all parameters

**Proved by the finite equations below.** Allow site-specific sleep
probabilities `p,r` in `(0,1)`, and let `q=1-p`, `s=1-r`. From 0 a jump goes
to 1 with probability `a`, otherwise to the sink; from 1 it goes to 0 with
probability `b`, otherwise to the sink. Assume `0<=a,b<=1` and `ab<1`.
Define `D=1-ab+abpr` and `H=1-abqs`, both strictly positive. The stationary
masses satisfy

```
mu(11) = pr/D,
mu(10) = ps[1-ab+abr(1-b)q]/(DH),
mu(01) = rq[1-ab+abp(1-a)s]/(DH),
mu(00) = 1-mu(11)-mu(10)-mu(01).
```

Here `10` means only vertex 0 is occupied. Their Rayleigh determinant is

```
mu(10)mu(01)-mu(00)mu(11)
 = abprqs [1-a+a(1-b)q] [1-b+b(1-a)s] / (D^2 H^2) >= 0.
```

For verification without trusting symbolic elimination, let `F_0,F_1`
be the generating polynomials for a lone active particle at 0 and 1, with
the other site empty, and let `F` be the law from one active at each site.
First-step analysis gives

```
F_0 = p z_0 + q(a F_1 + 1-a),
F_1 = r z_1 + s(b F_0 + 1-b),
D F = pr z_0 z_1 + psb(1-a) F_0 + ps(1-b) z_0 + q(1-ab) F_1.
```

The last equation follows by toppling vertex 0 first, then vertex 1, and
immediately moving any surplus. The two return terms have total weight
`ab(ps+q)=ab(1-pr)` and are moved to the left. Solving the first two equations
and comparing coefficients proves the displayed masses and determinant.
For a bivariate multiaffine polynomial this determinant criterion is
equivalent to real stability. Thus every such two-site law is strongly
Rayleigh, with no uniform-rate or reversibility assumption.

Fixed self-loops can be erased: replace `lambda_i` by
`lambda_i/(1-K_ii)` and every other jump probability in row `i` by
`K_ij/(1-K_ii)`. A self-loop does not change the configuration, and
conditioning on the next sleep or non-loop jump gives these parameters.
Together with component factorization, the two-site result therefore proves
strong Rayleigh whenever every strongly connected component has at most
two vertices, including self-loops. This is a restricted theorem, not
the general-kernel conjecture.



## All rates on three sites and stationary stronger-property failures

**Exact symbolic proof of NA on at most three non-sink vertices.** For an
irreducible simple undirected graph on three non-sink vertices, the internal
graph is a path or a triangle. Up to reflection, a path has five nonempty
sink-attachment sets (an endpoint; the center; both endpoints; an endpoint
and center; all three). A triangle has three (one, two, or all vertices).
`proof_checks.py` constructs the symbolic absorption law for each of these
eight cases, checks the full absorption equations and normalization, and
certifies all 15 increasing-event gaps as rational functions with
nonnegative numerator and denominator coefficients in `lambda`. Thus the
gaps are nonnegative for every `lambda>0`. Disconnected internal graphs are
covered by component factorization and the two-site theorem. This proves
NA on every simple sink graph with at most three non-sink vertices, at
every common sleep rate. This is a finite symbolic proof, not an empirical
claim for arbitrary graphs. Exact expressions are in `certificates.json`.

Some, but not all, of these laws are strongly Rayleigh. The script tests
each quadratic `a z^2+b z+c` by sign certificates for `a,c,4ac-b^2`.
Failure of that sufficient coefficient test alone is not a counterexample;
the following separate rational witnesses are checked against direct
absorption, every addition kernel, and the opposite toppling order.

**Strong Rayleigh fails on a plain path.** Let the full graph be
`sink--0--1--2` with common `lambda=3`. If `P` is its stationary generating
polynomial, then

`Delta_(0,1) P(z_2=-1/12) = -68136427/519532701427600 < 0`.

The multiaffine stability criterion therefore rules out strong Rayleigh.
The negative field does not disprove Rayleigh, and this law is NA by the
three-site certificate.

**Rayleigh itself fails in the undirected uniform-rate model.** Use internal
edges `{0,1},{1,2}` and sink edges `{0,sink},{1,sink}`, with common
`lambda=64`. At the positive external field `(z_0,z_1,z_2)=(1/100,1,1)`,

```
Delta_(1,2) P
 = -1356631018599034504351686468980948992
   /33379990838030144205088397496430705187032921875 < 0.
```

A more transparent explanation comes from the exact symbolic large-rate
expansion on the same graph:

\[
 \Pr(X_0=0)=\frac{3}{4\lambda}+O(\lambda^{-2}),\qquad
 \operatorname{Cov}(X_1,X_2\mid X_0=0)
       =\frac{4}{1215\lambda^2}+O(\lambda^{-3})>0
\]

for all sufficiently large `lambda`. Thus conditional negative correlation
also fails. The underlying unconditioned law remains NA for every
`lambda>0`. This rules out proving the full conjecture through Rayleigh,
strong Rayleigh, conditional NA, or negative regression in this generality.
It does not refute NA itself.

**Conditional negative correlation fails even on trees.** There is a
general domination obstruction. Work in one connected non-sink component
of an undirected sink graph. Suppose `b` dominates `c`, while `a,b` are
incomparable under domination. For the leading vacancy coefficients write
`C_i=C({i})`, and set `e_a=e_{ab}(a)`, `e_b=e_{ab}(b)`. Then

`Cov(X_b,X_c | X_a=0)`
` =lambda^(-2) e_a(C_b-e_b)(C_a C_c-C_{ac})/C_a^2
   +O(lambda^(-3))`,

and the displayed coefficient is strictly positive. To prove the formula,
domination gives `e_{abc}(c)=0` and leaves the other two escape
probabilities equal to `e_a,e_b`. Also `C_{bc}=C_b C_c`. The vacancy
recurrence therefore gives

`C_{abc}=C_c e_a(C_b-e_b)+C_{ac}e_b`,
`C_{ab}=e_a(C_b-e_b)+C_a e_b`.

Subtract `C_{ab}C_{ac}/C_a` and divide by `C_a`. This is precisely the
conditional covariance coefficient, since covariance is unchanged when
both occupation indicators are complemented. Incomparability and
connectedness give `e_a>0`, `C_b>e_b`, and `C_a C_c>C_{ac}`: `a,c` are
also incomparable, and the strict pair criterion applies. This argument
uses the proved leading vacancy expansion, not a finite-rate extrapolation.

For example, take edges `sink--0--1`, `1--2`, and `1--3--4`, and
`(a,b,c)=(2,3,4)`. Here

`C_2=1/3`, `C_{23}=7/150`, `C_{24}=7/96`, `C_{234}=17/1600`,

so `Cov(X_3,X_4 | X_2=0)=1/(800 lambda^2)+O(lambda^(-3))`.
`vacancy.py` checks these four coefficients by the escape recurrence,
independent region integration, and the IDLA series. At `lambda=64` the
exact conditional covariance is positive, while all 2,280 unconditioned
NA comparisons pass. Thus a tree-specific Rayleigh or conditional-NA
induction is unavailable in this generality. The star reduction below
concerns a narrower family and is unaffected.

**Targeted NA diagnostics.** Three four-site graphs add one pendant vertex
to 0, 1, or 2 in this last graph, keeping its two sink edges and common
`lambda=64`. These test whether extra vertices turn the conditional positive
correlation into an ordinary disjoint-event violation.

A separate parameter diagnostic keeps the same three-site graph and sets
`lambda_1=lambda_2=64` while varying `lambda_0` through `1/64,1,64,4096`.
This tests whether a fluctuating site 0 transmits the conditional correlation
to the unconditioned pair. These site-dependent-rate variants are distinct
from the main common-rate conjecture.

**Outcome of these diagnostics:** all three pendant extensions pass all
174 exact NA comparisons at `lambda=64`. All four rate vectors
`(1/64,64,64)`, `(1,64,64)`, `(64,64,64)`, and `(4096,64,64)` pass all 15
comparisons on the base graph. Thus these specific ways of trying to
promote the conditional violation to ordinary NA do not work.

## Harmonic polynomial identities

The following finite algebra also handles the cases where distinct harmonic
evaluation points coincide. It is used for higher large-rate coefficients
and the four-site classification.

**Exact change of basis.** Index a vector basis `e_B` by subsets of `V`.
Let `A_v` be the IDLA addition matrix: launch a walk at `v`, stop at the
first vacant vertex and occupy it, or at the sink. If `v` is vacant, it
is occupied immediately. These column-stochastic matrices commute by
abelian surplus routing, and

`A_v^2=A_v J_v`, where `J_v=K(v,sink) I+sum_w K(v,w) A_w`.

Let `N` have column `B` equal to the ARW stabilization law from one active
particle at each vertex of `B`, and let `T_v` be addition followed by ARW
stabilization on stable states. Moving only surplus particles first is
legal in ARW and leaves fresh stack suffixes, so `N A_v=T_v N`.
For sleep probability `p_v` and `q_v=1-p_v`, put

`B_v=(A_v-q_v J_v)/p_v`,
`M e_B=(product_(v in B) B_v)e_empty`.

Then `N M=I`. Indeed, for an empty site `v` in a stable configuration `B`,
the first toppling after addition gives

`T_v e_B=p_v e_(B union {v})+q_v J_v(T)e_B`.

Thus `(T_v-q_v J_v(T))/p_v` adds a sleeper deterministically whenever `v`
is empty. Intertwining through `N` and adding the distinct sites of `B`
one at a time proves the assertion for every column. In particular `M`
is invertible, and the stationary law is exactly `M^(-1)e_V`.
This identity does not require diagonalizability or distinct eigenvalues.

For `lambda_v=alpha_v/t`, `B_v=A_v+(t/alpha_v)(A_v-J_v)`. Thus
`M(t)=I+sum_(j=1)^n t^j M_j` is polynomial. Writing
`mu(t)=sum_(r>=0) t^r u_r`, its Taylor coefficients are given without
symbolic rational-function expansion by

`u_0=e_V`, `u_r=-sum_(j=1)^min(n,r) M_j u_(r-j)`.

**Exact rate evolution.** Write `x_v=1/lambda_v` and `p_v=1/(1+x_v)`.
Let `D_v e_B=1_(v in B)e_(B\{v})` and
`U_v e_B=1_(v notin B)e_(B union {v})`. Differentiating the column products
defining `M` gives

`partial_(x_v) M=(A_v-J_v) M D_v`.

Using `N=M^-1`, intertwining, and the first-toppling equation on empty
sites yields

`partial_(x_v) N=p_v (R_v-I)N`,
`R_v=I+[J_v(T)-U_v]D_v`.                                    (18)

Here `R_v` is a stochastic kernel: leave an empty `v` unchanged; otherwise
wake its sleeper, force its first jump according to `K`, then stabilize.
It never increases particle number. In particular, along `x_v=t/alpha_v`,

`d mu_t/dt=sum_v (alpha_v+t)^(-1)(R_v(t)-I)mu_t`,
`mu_0=delta_V`.

This supplies an exact all-rate evolution of the stationary family.
`algebra.py` checks the three separate derivatives of the sitewise
three-site law at inverse rates `(1,2,3)` and verifies stochasticity of
the kernels. It does not supply an NA-preservation theorem: the product
input counterexample above already rules out such a theorem for `N` on
arbitrary NA inputs. An argument using (18) must exploit the stationary
family, or a narrower initial class.

`algebra.py` implements these identities and checks the IDLA matrix
relations. Its laws and seventh-order expansions were compared with the
independent symbolic absorption laws on the three-site endpoint path and
the three-site conditional-correlation example. At `lambda=3`, the
four-site endpoint-path law agrees with raw stabilization, both toppling
orders, and every addition-stationarity equation. On endpoint paths with
three through six non-sink vertices, all computed pair gaps have a strictly
positive coefficient at order `t^3`; their order-`t^2` coefficients vanish.
These checks agree with the all-size formula (10).

**Finite symbolic theorem.** Every simple sink graph with at most four
non-sink vertices has an NA stationary law for every common `lambda>0`.
For connected internal graphs, enumeration up to permutations fixing the
sink gives exactly 44 four-site types. `four_site.py` computes each law by
the change of basis above. Write its sixteen masses as `W_B(t)/D(t)`,
where `t=1/lambda` and `D` has nonnegative coefficients and is positive
for `t>0`. For all 174 unordered disjoint nonconstant increasing-event
pairs `E,F`, the polynomial

`(sum_(B in E) W_B)(sum_(B in F) W_B)-D sum_(B in E intersect F) W_B`

has nonnegative coefficients. All 44 cases pass. The complete exact laws
and coverage flag are saved in `four_site_certificates.json`; the script
reconstructs every gap from them. Disconnected internal graphs factor
over components, each already covered by the smaller cases. Increasing
event inequalities imply the full function formulation by threshold
decomposition. This is a finite symbolic proof, not an all-size result.

The selected five-site graph with edges
`sink--0, 0--1, 1--2, 1--3, 2--3, 3--4` also passes all 2,280
increasing-event comparisons by the same coefficientwise argument, for
every common sleep rate. Its exact law and certificate are in
`five_site_certificates.json`. This tests whether the all-rate polynomial
positivity extends through a cyclic block below a cut vertex; it is not a
classification of five-site graphs.

A second five-site test joins two copies of the conditional-correlation
example at their common vertex 0: edges `0--1--2, 0--3--4`, with sink
edges at `0,1,3`. All 2,280 inequalities again have coefficientwise
nonnegative numerators for every common rate. Thus this specific way of
combining the positive conditional correlations does not produce an
unconditional NA counterexample. Its law is also in
`five_site_certificates.json`.

A separate path calculation proves **Rayleigh**, as well as NA, on
`sink--0--...--(n-1)` for `n=3,4,5` and every common sleep rate. After
clearing a common positive polynomial denominator in `t=1/lambda`, all
coefficients of every Rayleigh difference are nonnegative jointly in
`t` and the remaining positive external-field variables. The full NA
event gaps were checked separately; their counts are `15,174,2280`.
The law and both certificates are in `endpoint_path_certificates.json`.
This is an all-rate, bounded-length statement. The three-site negative
real-field witness above still rules out strong Rayleigh, and the
branching-tree obstruction prevents extending this conclusion to all trees.



## Complete-graph identities and deletion sources

**Proved identity.** Fix real numbers `a_v` and put
`j_v=K(v,sink)+sum_w K(v,w)a_w`. Suppose `a_v^2=a_v j_v`, and set
`b_v=(a_v-(1-p_v)j_v)/p_v`. For the stationary generating polynomial `P`,

`P(b)=product_v a_v`.                                           (5)

To check it, assign an active count `m` at `v` weight `a_v^m`, a sleeper
weight `b_v`, and an empty site weight one, and multiply weights across
vertices. This function is harmonic under legal topplings. A lone particle
uses `a_v=p_v b_v+(1-p_v)j_v`; a surplus uses `a_v^2=a_v j_v`. When a jump
wakes a sleeper, use `a_v b_v=a_v^2`, which follows from the same identities.
There are finitely many configurations with at most `|V|` particles, so
the harmonic function is bounded and optional stopping at stabilization
is justified. Start from `1_V` to obtain (5).

One family of solutions is indexed by `J subset V`: set `a=0` on `V\J`
and, on `J`, let `a_v` be the probability of reaching the sink before
`V\J`. Then `a=j` on `J`, so `b=a` there, while `b=-j/lambda_v` outside.
Every proper `J` therefore gives a real zero of `P`; `J=V` gives `P(1)=1`.
Some points can coincide when parts of `J` cannot reach the sink inside `J`;
these identities alone must not be assumed to determine `P` in that case.

For the complete graph on `n` non-sink vertices plus one sink, simple random
walk, and common `lambda`, exchangeability reduces (5) to `n+1` linear
equations for `r_k=Pr(|X|=k)`. If `j=|J|<n`, the test values are

`a_j=1/(n+1-j)`, `b_j=-(n+1)/(n lambda (n+1-j))`,

and the equation is

`sum_k r_k sum_l binom(j,l)binom(n-j,k-l) a_j^l b_j^(k-l)/binom(n,k)=0`.

Together with `sum_k r_k=1`, these determine the rank law. Here is why the
system is invertible, not just a numerical observation. The function
`h_J(B)=product_{v in B} b_v` is a right eigenfunction of each addition
operator, with eigenvalue `a_v`: apply the same harmonic argument to an
initial sleeper set plus one added particle, using `a_v b_v=a_v^2` when
the addition wakes a sleeper. Each `h_J` is nonzero at the empty set.
With direct killing possible at every vertex, `a_v>0` exactly on `J`.
The `2^n` joint eigenvalue vectors are distinct. A generic linear combination
of the commuting addition operators therefore has these `2^n` eigenvectors
with distinct eigenvalues, so they form a basis. Their expectations determine
the stationary law. For an exchangeable signed law, all subsets `J` of the
same size impose the same condition; thus the reduced system is invertible.

**Independent exact verification and remaining lemma.**
[Junge–Kaufman–Meisel, Proposition `prop:stop` and its proof](https://arxiv.org/html/2604.04747v2)
give a second representation in terms of a binomial-refresh chain and a
descending boundary. Their looped model must be converted: to obtain our
simple random walk with sleep rate `lambda`, use refresh parameter

`rho=n lambda/(n+1+n lambda)`.

Let `Y` have birth rate `beta_y=rho(n-y)` and death rate
`delta_y=(1-rho)y`, with initial law `pi_y=binom(n,y)rho^y(1-rho)^(n-y)`.
Independently start a boundary at `n` and let it decrease at rate one.
Stop when the two meet. The terminal value has rank law `r`. The cited
proof routes particles through a holding vertex; each arrival refreshes
one Bernoulli coordinate and each sink event lowers the boundary. No NA
theorem is imported from that paper.

The resulting recursion is particularly small. Let `L_b` be the generator
of `Y` restricted to `0,...,b-1`, killed on hitting `b`. Put
`f_(n+1)(y)=1` and, for `b=n,...,1`, define

`f_b=(I-L_b)^(-1)(f_(b+1) restricted to {0,...,b-1})`.

Reversibility with respect to `pi` and boundary flux give

`r_k/pi_k=f_(k+1)(k)+delta_k f_k(k-1)`,

with the second term zero when `k=0`. Thus the next precise candidate lemma
is log-concavity of this last sequence in `k`. It is equivalent to
ultra-log-concavity (ULC) of `r`, meaning log-concavity of
`r_k/binom(n,k)`.

For completeness, ULC suffices for NA here. Given a fixed pattern on some
coordinates with `j` occupied sites, the remaining configuration has weight
`q_(j+ell)` on each set of size `ell`, where `q_k=r_k/binom(n,k)`.
Log-concavity makes `q_(j+1+ell)/q_(j+ell)` decrease with `ell`, so increasing
the conditioned number of occupied coordinates stochastically decreases
the remaining count. Uniform subsets of increasing sizes have a nested
coupling. Hence the conditional expectation of an increasing function
on a disjoint block decreases with the conditioned count. The conditional
average of the other increasing function increases with that same count;
their scalar covariance is nonpositive. This proves NA.

`complete_graph.py` verifies agreement of harmonic interpolation and the
boundary recursion for `n=2,3,4,5,8,12` at `lambda=1`, and `n=8` at
`lambda=1/64,64`. For `n<=4` it also agrees with direct stabilization and
every addition-stationarity equation. All these exact laws are ULC, hence
NA. Root counts at `lambda=1` are `n` for `n<=5`, six for `n=8`, and ten
for `n=12`; at `n=8,lambda=1/64` there are only two real roots. Thus strong
Rayleigh fails even for these symmetric complete graphs. Exact rank laws,
ULC gaps, and root counts are in `complete_graph_certificates.json`.

**Gap:** log-concavity of the boundary-flux sequence is unproved for general
`n,rho`. It cannot be justified merely by claiming that each resolvent
preserves NA: independent Bernoulli coordinates started empty, refreshed
for one common independent exponential time, have positive covariance
because their shared random time mixes the success probability. The proof
must use the specific incoming densities and consecutive boundary levels.

**Direct pair inequalities from bounded stopping harmonics (proved at a
specified field; the unmodified sign is proved below by deletion).** In the descending-boundary
representation write p=rho, q=1-p, T for the meeting time and S for
the terminal count. Put L=n-S and let N_t count boundary jumps.
Since T is bounded by the nth boundary arrival, stopped generator
identities are justified. The pair-NC slack has the exact forms

`(n-1)(E S)^2/n-E[S(S-1)]`
` =(n-1)(E T)^2/n-2 E integral_0^T N_s ds`
` =(n-1)(E S-np)^2/n+2 E integral_0^T K_2(Y_s) ds`,

where `K_1(y)=y-np` and
`K_2(y)=y(y-1)-2(n-1)p y+n(n-1)p^2`.
Indeed E L=E T, `E[L(L-1)]=2 E integral_0^T N_s ds`, while
the count generator gives `L_Y K_j=-j K_j` for j=1,2. Their initial
binomial means vanish. Neither occupation integrand has the needed sign.

The bounded functions `2^(-b)K_1(y)` and `3^(-b)K_2(y)` are harmonic
for the joint chain: the boundary generator contributes respectively
1 and2, cancelling the count eigenvalues. Optional stopping gives

`E[2^(-S)(S-np)]=0`, `E[3^(-S)K_2(S)]=0`.

For the rank law tilted by z^S let u(z) be its mean. The first identity
gives u(1/2)=np. The second implies, for distinct sites and n>=2,

`Cov_(z=1/3)(X_i,X_j)=-[u(1/3)/n-p]^2<0`.

Strictness follows because the full-support rank law has strictly increasing
u(z), with derivative Var_z(S)/z. This proves an all-size pair inequality
at that uniform positive field; it does not establish it at z=1 or prove NA.

The two terminal identities alone cannot do so. At n=3,p=1/2, the
artificial rank law `(1,8,13,10)/32` satisfies both, but its mean is2
and second factorial moment43/16, exceeding `(2/3)*2^2` by1/48.
Its exchangeable pair covariance is1/288. This is not the stopping law:
it fails the next bounded harmonic identity. The authoring AI checked the
identities and this exact arithmetic independently. Any continuation
must use the actual stopped occupation law or additional structure.

**Maximum and occupation representation (proved; the curvature inequality
follows from the deletion theorem below).**
For the n-bit chain started at i define
`M_i=sup_(t>=0)(Y_t-N_t)` and `u_i=E_i M_i`, with an independent
rate-one Poisson N. The maximum is attained before its nth arrival and
belongs to [i,n]. With stationary initial count pi=Bin(n,p), this maximum
has exactly the terminal-rank law S. Indeed
`f_b(i)=Pr_i(M_i<b)` obeys the rising-boundary recursion. Extending it
by zero at i>=b, reversibility and pi L_Y=0 give

`sum_i pi_i[f_(k+1)(i)-f_k(i)]`
` =pi_k[f_(k+1)(k)+q k f_k(k-1)]`,

the previously verified descending-boundary flux. Also the actual phase
occupation density before its stopping time is

`rho_i=E integral_0^T 1_(Y_t=i)dt=pi_i(n-u_i)`.

This follows by summing the forward stage occupations
`pi_i sum_(b=i+1)^n f_b(i)` and applying the integer tail-sum identity.

Couple all initial active sets A using common future refreshes. If U_t
is the nested set of unrefreshed sites, then
`M(A)=max_t[Y_empty(t)-N_t+|A intersect U_t|]`.
This set function is increasing, integer valued, one-Lipschitz and
supermodular. For the last assertion take maximizing times for A,C,
ordered so their unrefreshed sets are U_j superset U_l. Evaluating the
union at j and the intersection at l improves their sum by
`|(C minus A) intersect (U_j minus U_l)|>=0`. Thus
`r_E(V)=n-M(V^c)` is a matroid rank function for each future environment E.
In particular, uniform-set averaging proves

`0<d_0<=...<=d_(n-1)<1`, `d_i=u_(i+1)-u_i`.

The strict bounds follow from positive-probability environments giving
M(A)=|A| or M(A)=n. If the first event is a boundary jump followed by
n distinct successful refreshes before the next boundary jump, then
M(A)=max(|A|,n-1). This has probability
`n! p^n/(n+1)^(n+1)>0`, so the final second difference is strictly positive.

Substitute rho into the stopped moment identities and integrate by parts
under the binomial law. With mu=E S,

`mu-np=npq E_(Bin(n-1,p)) d_J`,
`E K_2(S)=2n(n-1)p^2 q^2 E_(Bin(n-2,p))(d_(J+1)-d_J)`.

Consequently pair NC is exactly the scalar inequality

`[E_(Bin(n-1,p)) d_J]^2`
` >=2 E_(Bin(n-2,p))(d_(J+1)-d_J)`.

The common factor in its slack is n(n-1)p^2 q^2. In particular
`np<mu<np+npq`, while `E K_2(S)>0` and
`E integral_0^T K_2(Y_t)dt<0`. Hence the mean-square term must
compensate a strictly negative quadratic occupation term.

There is a second exact description of the missing control. For independent
Bernoulli(q) initial vacancies V put L=r_E(V), `m_E=E[L|E]`, m=E L.
Conditional Bernoulli Poincare gives `Var(L|E)<=p m_E`: single-element
rank differences are zero or one, and q times their sum counts expected
coloops, at most expected rank. Put `C_E=p m_E-Var(L|E)>=0`. Then

`m(1-m/n)-Var(L)=m(q-m/n)+E C_E-Var(m_E)`.

Random-matroid structure alone cannot bound the last variance. For example,
at n=2,p=q=1/2 mix the free, rank-one uniform and rank-zero matroids with
probabilities1/2,1/4,1/4. This gives u=(3/4,5/4,2), d=(1/2,3/4) and
`Pr(M=0,1,2)=(1/8,7/16,7/16)`. It satisfies both displayed maximum-moment
identities and the global occupation balance
`Pr(M=i)=pi_i[1+L_Y(n-u)(i)]`, but its pair covariance is7/1024>0.
This is not an Ehrenfest future: it fails the bounded terminal harmonics.
Thus these abstract shape and balance properties do not suffice on their
own. The authoring AI independently checked the maximum proof, all factors and
this exact witness. The deletion argument below uses the actual refresh
and boundary dynamics to prove the required curvature inequality.

**Two initial bits (proved reduction; clock pairing open).** Fix the
other initial bits and start two distinguished sites at zero, with common
refreshes for all four choices of those bits. Conditional on site1
refreshing first, let t1<t2 be their first refresh times. Let A,B,C be
the maxima of Y_00-N on [0,t1), [t1,t2), [t2,infinity), respectively.
Then, pathwise,

`M(a,b)=max(A+a+b,B+b,C)`.

The mixed initial difference is exactly `1_(A>=B,C=A+1)`: if A<B,
the B term dominates the A term even after increasing the first bit;
if A>=B, check the three consecutive maxima max(A+j,C), j=0,1,2.
Average the other initial bits independently with Bernoulli(p), and put

`gamma=Pr(A>=B,C=A+1)`, `U=Pr(A>=B,C<=A)`,
`V=Pr(max(A,B)>=C)`.

The two influences, conditional on this refresh order and averaged over
the other distinguished bit, are a=U+p gamma and b=V+p gamma.
Exchangeability gives `dbar=(U+V)/2+p gamma`, and gamma is the
binomially averaged second difference above. There is no further
order factor. Pair NC is therefore exactly `dbar^2>=2 gamma`.

The sufficient comparison `2 gamma<=UV` is false for actual clocks.
For n=2, write D=1+p+p^2 and H=(1+2p)D. The actual first-record
recursion gives

`U=(1+p)/H`, `V=(1+2p+2p^2)/H`, `gamma=p/H`.

At p=2/3 (ordinary sleep rate3), these are (45,87,18)/133, so
`UV-2 gamma=-873/17689`, whereas
`dbar^2-2 gamma=1296/17689>0`. Retaining the initial-bit contribution
leaves the potentially sufficient comparison `2 gamma<=ab` open.
For n=2 its slack is `(1-p^2)/[(1+2p)^2 D]>0`, with
`a=1/(1+2p)` and `b=(1+2p+3p^2)/H`.

On the gamma event, inserting a boundary point anywhere in (t1,t2)
makes both U and V hold: A is unchanged, B can only decrease and
C drops to A. This is not a product-probability bound. Poisson insertion
weights the input by t2-t1 and the output by its number of inverse
deletions, and both successes occur in one environment. The interval
length can be arbitrarily small even on the actual gamma event.
No measure-controlled construction of two independent environments was
obtained. The reduction and two-site obstruction were independently
checked; this coupling route is parked pending such a construction.

**Ascending records and a weighted Stein identity (proved; pair sign
open).** Put beta_i=p(n-i), delta_i=qi. For initial phase i let
`R_i=M_i-i`, and let `L_(j,i)` be the probability that the first new
maximum of Y-N above i occurs at phase j. Necessarily j>i, so L is
strictly lower triangular in column convention and substochastic.
With `G_i(z)=E z^(R_i)`, the strong Markov property gives

`G_i(z)=1-sum_(j>i)L_(j,i)+z sum_(j>i)L_(j,i)G_j(z)`.

Define matrices by `B e_i=beta_i e_(i+1)`,
`D e_i=delta_i e_(i-1)`, C=I+D, and
`A=diag(1+beta_i+delta_i)`. First-step decomposition gives the actual
Riccati identity

`L A=B+L^2 C`.

A birth immediately sets a record; a death or boundary jump drops the
current level by one, requiring two successive records to exceed the
original maximum. In particular
`L_(i+1,i)=beta_i/[1+beta_i+delta_i-delta_i L_(i,i-1)]`,
omitting the death term at i=0. Since Q=B+C-A is the count generator,
`(I-L)[A-(I+L)C]=-Q`. Apply this to pi=Bin(n,p) and use detailed
balance to obtain

`L[(1+beta) odot pi]=delta odot pi`.

Now draw I with law pi and R with conditional law R_I. Weight the
first-record decomposition by `(1+beta_i)pi_i`. For every function h,

`E[(1+np-pI)h(R-1)1_(R>=1)]=q E[I h(R)]`.

Writing g=E R=E S-np and taking h=1,r,r^2 gives

`E[(1+np-pI)1_(R>=1)]=npq`,
`E[IR]=(1+np)g-npq`,
`Cov(I,R^2)=E R^2-2q(1+np)g+npq(1-2p)`.

The common-refresh coupling has `R_i-R_(i+1)` in {0,1}, with strictly
positive mean. Therefore E R_i^2 strictly decreases with i and
`Cov(I,R^2)<0`. Substitution gives another exact form of the pair slack:

`Delta=(n-1)g^2/n-(3+2npq)g+npq(3-2p)-Cov(I,R^2)`.

The negative sign of this covariance alone is insufficient: at the
actual n=2,p=1/2 law, g=9/28 and Cov(I,R^2)=-1/4; the polynomial
part is -367/1568 and the full slack is25/1568. A quantitative
estimate using the actual record dynamics is still needed. Positivity
of L is an additional constraint absent from the earlier n=2 synthetic
matroid witness: its record PGFs would force L_(1,0)=1 and
L_(2,0)=-1/2. Neither this observation nor the weighted identity proves
pair NC. The authoring AI independently checked the Riccati, flow and moment
identities.

The same joint law has a stationary reflected-queue realization. Extend
stationary Y and independent Poisson increments to the whole real line,
and put `W_t=sup_(s<=t)[Y_s-(N_t-N_s)]`. The supremum is finite, lies
between Y_t and n, and is stationary. Time reversal of Y and independent
Poisson increments identifies `(W_t,Y_t)` with `(M_I,I)`. Thus I=Y,
R=W-Y has the finite-state generator

`(i,r)->(i+1,(r-1)_+)` at rate p(n-i),
`(i,r)->(i-1,r+1)` at rate qi,
`(i,r)->(i,(r-1)_+)` at rate1,

on i,r>=0, i+r<=n; omit zero-gap self-events. The chain is irreducible:
deaths and boundary events reach (0,0), while births followed by deaths
reach every state. Its R-level crossing balance gives exactly the
weighted Stein identity above. The terminal rank W increases only at
r=0 and decreases only at r>0.

There is an exact labelled particle realization of this queue. Give each
site colour A (active), B (buffered), or0 (empty), with counts i,r,h.
Each A becomes B at rate q and each B becomes A at rate p. Each empty
site requests activation at rate p: if r=0 it becomes A; otherwise
choose a B uniformly, change the requesting site to A and the donor to0.
At total rate1 remove a uniformly chosen B, if one exists. The count
projection has exactly the preceding generator. Moreover the A-indicators
themselves have the independent refresh generator: each inactive site
becomes active at rate p and each active site becomes inactive at rate q;
a donor was already inactive. This is an intertwining of the entire
active field, not just its count.

The coloured chain is irreducible. Convert all A to B and remove them
to reach empty; from there, create any prescribed occupied set as A
before converting its desired B sites. Its unique stationary law is
exchangeable, and its occupied count i+r has law S. Conditional occupied
sets are uniform, so its stationary occupied-set law is exactly the
complete-graph ARW law. Occupancy alone is not Markov: an all-A occupied
set admits births, whereas an equal occupied set containing B has removal
and relocation. A negative-dependence proof for this coloured process
would therefore transfer, but none is established here.

Two familiar representation shortcuts fail. A product over the three
compartment counts, `mu(i,r,h)=constant*a_i*b_r*c_h`, would imply
`mu(2,1,0)mu(0,2,1)mu(1,0,2)`
` =mu(2,0,1)mu(0,1,2)mu(1,2,0)`.
At the actual n=3,p=1/2 chain, the factors on the left are
(11/120,1/30,11/60), and those on the right are
(17/60,1/20,1/24). The product difference is -13/432000.
Multinomial orbit factors can be absorbed into the three functions and
do not avoid this obstruction. The ten stationary balance equations
and the labelled projection on all27 states were checked exactly.

Nor can an irreducible finite symmetric-exclusion process with independent
site reservoirs realize S as the count on a fixed subset, even with extra
sites. Symmetric exchanges preserve strongly Rayleigh laws by
[Borcea--Branden--Liggett, Proposition5.1](https://arxiv.org/pdf/0707.2340).
Independent site births and deaths preserve stability as well: if
F=A+zB and their time-t flip probabilities are u,v, the new polynomial is
`[(1-u)+uz]F((v+(1-v)z)/((1-u)+uz))`.
Here 1-u-v>0, so the fractional linear map preserves the upper half-plane.
Trotter products and stationary limits give a strongly Rayleigh stationary
law. Its fixed-subset count polynomial must have only real zeros, but the
existing n=8,lambda=1 certificate has a square-free degree-eight rank
polynomial with exactly six real zeros. This rules out that representation,
not NC/NA, arbitrary nonlinear observables, or an arbitrary mixture of
conserved exclusion sectors. The latter already prescribes a rank law.

The labelled process admits an actual one-discrepancy coupling. If two
colour configurations agree except that one has an extra A or B at an
empty site, their occupied sets can remain nested with at most one extra
particle until coalescence. For buffer sets D and D union {z}, choose
the same uniform element of D with probability |D|/(|D|+1); otherwise
choose z in the larger buffer and an independent uniform element of D
in the smaller. If D is empty only the larger selects. This has uniform
marginals. A differing choice moves the discrepancy to the smaller
system's donor. Service can remove it, and a successful refresh at its
site either moves it to a donor or coalesces the systems. An unsuccessful
refresh turns an extra A into B. All other events agree. Consequently,
for increasing occupancy f, the backward equation from empty gives
`d E_empty f(O(t))/dt=np[E_(one uniform A)f(O(t))-E_empty f(O(t))]>=0`.
This is stochastic growth in time, not pathwise growth.

There is also a conditional negative-dependence statement for the actual
empty-start law. Write A_x(t) for activity and O_x(t) for occupancy. For
every t>0 and every tag x,

`Law(O_-x(t) | A_x(t)=1)<=st Law(O_-x(t) | A_x(t)=0)`.   (55act)

It therefore holds in the stationary limit as well. Here is the full
coupling argument. Represent all activity by independent rate-one refresh
clocks and Bernoulli(p) coins, adding independent rate-one service. Let
H mean that x refreshes before t. Conditional on H and either final
activity, the tag clock times and all earlier coins have the same law;
only the last coin is forced to1 or0. Couple their identical prefixes.
After the last tag refresh, keep x active in the plus system and inactive
in the minus system. The possible discrepancy forms are

* x is0 below and A above, with all other colours equal;
* x is B below and A above, with all other colours equal;
* x is0 below and A above, and one other y is B below and0 above.

The uniform buffer-choice coupling preserves these forms while the tag
does not refresh. In the second form, selecting the tag below and a
common donor above creates the third form, or the first if no common
buffer remains. In the third form, selecting the discrepant B below
moves the other-site discrepancy to a common donor, or removes it.
A successful refresh at that other site does the same. All remaining
events agree. Thus the plus occupied set outside x stays contained in
the minus one. The three forms are exactly those produced by the final
coin change, according to the pre-refresh tag colour.

The remaining branch A_x(t)=0,H^c has x permanently empty, and the other
sites run the size-(n-1) process from empty. Couple their active paths
and service clock N with the plus process conditioned on A_x(t)=1.
The independent tag path changes the total active count Y by either0
or1 at each time. Occupied counts from empty satisfy, pathwise,
`W(t)=sup_(0<=s<=t)[Y(s)-(N_t-N_s)]`.
Hence the two total W values differ by0 or1. The plus tag is occupied
at t, so its other-site occupied count is at most the minus count.
The tag conditioning does not alter the latter marginal. Both conditional
laws are exchangeable outside x; coupling their ordered counts with a
common uniform permutation lifts this comparison to occupied sets.
Mixing the H and H^c branches proves (55act).

The same domination holds after prescribing every other terminal active
bit. The frozen-tag coupling keeps their whole active paths identical;
the no-refresh reflection comparison still works with those prescribed
endpoints. The sites prescribed active are occupied in both copies, and
conditional buffer sets are uniform among the prescribed inactive sites.
Thus ordered buffer counts again lift to configuration domination.

**Every vacant block dominates its trace (proved, all sizes and common
rates).** For an ordinary n-site complete graph at common lambda>0 and
a prescribed deleted block D, the traced stationary law nu_D satisfies

`nu_D<=st Law(X_(V minus D) | X_D=0)`.                 (55vac)

This proves the right conditional trace comparison; the left comparison
is obtained below by unconditioned deletion. Here is a positive-generator proof.
Write pi_n(i,r) for the stationary coloured count law,
rho_i=binom(n,i)p^i q^(n-i), theta_i(r)=Pr(R=r|I=i), and
F_i(r)=Pr(R<=r|I=i). The refined active-bit comparison says that R
given I=j+1 is dominated by R given I=j after deleting one uniformly
chosen inactive site. Hence, for1<=r<=n-j,

`F_(j+1)(r-1)>=F_j(r-1)+r theta_j(r)/(n-j)`.

For1<=s<=n-1 and0<=j<s define
`kappa_(s,j)=p(n-s)pi_n(j,s-j)-sum_(a=0)^j pi_n(a,s+1-a)`.
Balance across {I<=j,I+R<=s}, followed by binomial detailed balance,
with r=s-j, rewrites this as
`p(n-j)rho_j[F_(j+1)(r-1)-F_j(r-1)-r theta_j(r)/(n-j)]>=0`.
Set kappa_(s,-1)=kappa_(s,s)=0; the latter equality is occupied-level
crossing balance.

Let k=|D|, m=n-k>=1, h=n-i-r, and H_k=E[(n-I-R)_k], where
(h)_k=h(h-1)...(h-k+1). Conditional on D being empty, the remaining
count law is `nu_k(i,r)=(h)_k pi_n(i,r)/H_k`, i+r<=m.
Freeze D empty and evolve with the size-m generator G_m at the same p.
Its forward source is exactly

`G_m^*nu_k=(1/H_k)sum_(s=1)^m k(n-s-1)_(k-1)`
`        *sum_(j=0)^(s-1) kappa_(s,j)`
`          [e_(j,s-j)-e_(j+1,s-j-1)]`.                 (55src)

For verification, put w(h)=(h)_k. Birth rates drop by pk, while the
service inflow changes weight by
`w(h-1)-w(h)=-k(h-1)_(k-1)`. The component of H_k G_m^*nu_k is
`k(h-1)_(k-1)[p h pi_n(i,r)-p h pi_n(i-1,r+1)-pi_n(i,r+1)]`.
On level s=i+r the bracket is kappa_(s,i)-kappa_(s,i-1).
At r=0 the extra birth inflow cancels because
`w(h+1)p(h+1-k)=w(h)p(h+1)`; at s=m the excluded service source
has w(k-1)=0 and the factor is k!=w(k). The s=0 component vanishes
by balance at (0,0). This checks every boundary in (55src).

At a fixed initial occupied count s, changing one B to A stochastically
increases every later occupied count. Couple the active bits in order
and use the same service clock; the reflection is
`W(t)=max{s-N_t,sup_(0<=u<=t)[Y(u)-(N_t-N_u)]}`.
It preserves this order. Pair (55src) with the semigroup applied to an
increasing rank function. Every coefficient is nonnegative and every
paired difference is nonpositive, so the occupied count from nu_k
decreases stochastically in time. Its limit is the stationary size-m
law at p. Exchangeability outside D lifts the rank comparison to sets.

This smaller law is exactly the physical trace. Tracing k deleted sites
gives holding probability k/[n(m+1)] and jump probability
(n+1)/[n(m+1)] to each other surviving vertex, including the sink.
Erasing holding therefore gives the ordinary m-site graph at
`lambda_tr=lambda n(m+1)/[m(n+1)]`. Its refresh parameter is the same
p=n lambda/(n+1+n lambda), proving (55vac). Deleting all sites is trivial.
The comparison does not use the original unadjusted loopless rate on
the smaller graph. authoring AI and reviewing AI separately checked the region balance,
factorial source, boundary terms and trace conversion.

There is an auxiliary left comparison for every block as well:
`Law(O_(V minus D) | A_D=1)<=st nu_D`.
From empty, share the outside active paths and service clock with a
size-m system. The k independent tagged paths add between0 and k to
the active count, hence between0 and k to the reflected occupied count.
At observation all k tags are occupied; subtracting them gives the
outside count inequality. Conditioning their terminal activity leaves
the smaller marginal unchanged. Exchangeability and stationary limits
give the displayed comparison. A_D is auxiliary activity, not occupancy;
the occupied-block replacement now follows from the threshold theorem below.

**Unconditioned deletion, the trace sandwich and strict pair NC (proved,
two independent analytic audits).** For every ordinary complete graph
with n>=2 non-sink sites and every common positive finite sleep rate,
every pair of distinct occupation indicators has strictly negative
covariance. More generally, for every prescribed nonempty proper D,
write E_D={some site of D is occupied}. Then

`Law(X_outside)<=st nu_D<=st Law(X_outside | X_D=0)`,
`Law(X_outside | E_D)<=st nu_D<=st Law(X_outside | E_D^c)`. (55sand)

Thus `Cov(1_(E_D),f(X_outside))<=0` for every increasing f. For one
deleted site this proves one-coordinate NA and all three comparisons
in (29). For larger blocks E_D is an OR event, not full occupation.
Full NA and full ULC remain open; occupied-block trace comparisons are
strengthened to every upper threshold below.

Keep p fixed and put m=n-1. The unconditioned coloured outside count
law after deleting a specified tag is
`nu(i,r)=[(n-i-r)pi_n(i,r)+(i+1)pi_n(i+1,r)`
`                              +(r+1)pi_n(i,r+1)]/n`.
For1<=s<=m and0<=i<s set
`D_(s,i)=e_(i,s-i)-e_(i,s-i-1)`.
Its exact source under the smaller generator is

`G_m^*nu=(1/n)sum_(s=1)^m sum_(i=0)^(s-1)`
`                 [kappa_(s,i)-kappa_(s,i-1)]D_(s,i)`. (55del)

To verify it, compare the full generator on outside functions with
G_m and use full stationarity. An active tag makes no difference. An
empty tag removes an outside buffer at rate p, contributing the positive
source `p(n-i-r)pi_n(i,r)D_(s,i)/n`. A buffered tag protects the outside
from service with probability1/(r+1), contributing
`-pi_n(i,r+1)D_(s,i)/n`. It can also supply the donor at an empty
outside site's successful refresh. After shifting the active index this
contributes `-p(n-i-r)pi_n(i-1,r+1)D_(s,i)/n`.
The bracket is exactly the difference of kappa values from (55src).
When the outside buffer is zero, both buffered-tag discrepancies vanish:
service affects no outside site, and an activation has the same outside
output whether it creates a particle or imports the tag. Thus no i=s
term exists. Negative-index probabilities vanish. At the top s=m the
full-state inflows remain valid; the bottom-state component is the lower
endpoint of D_(1,0), not zero.

Summation by parts in (55del) gives, at each s,
`sum_(i=0)^(s-2)kappa_(s,i)[D_(s,i)-D_(s,i+1)]`
`                         +kappa_(s,s-1)D_(s,s-1)`.
Every kappa is nonnegative by the actual active-bit coupling above.
Let P_t be the smaller semigroup and f any increasing rank function.
For identical initial active count i and initial occupied counts s,s-1,
use the same active paths Y^i and service N. Put
`A_i(t)=sup_(0<=u<=t)[Y^i(u)-(N_t-N_u)]`.
Reflection gives `W_(i,s)(t)=max{s-N_t,A_i(t)}`, so

`<D_(s,i),P_t f>=E[(f(s-N_t)-f(s-1-N_t))`
`                              *1_(A_i(t)<=s-1-N_t)]`.

Set the integrand to zero when s-1-N_t<0; the event is then impossible.
Nested initial active sets and common refreshes give A_i<=A_(i+1).
The nonnegative f increment depends only on N, hence these pairings
are nonnegative and nonincreasing in i. The adjacent comparison uses
only i<=s-2, where both lower initial states are valid. Therefore
`d<nu,P_t f>/dt=<G_m^*nu,P_t f>>=0` for every t>=0.
The finite irreducible chain converges to its stationary law at the same
p. Exchangeability lifts ordered occupied counts to ordered sets by a
common uniform permutation. This proves the unconditioned one-deletion
comparison. Project and iterate it for any D; trace is associative and
preserves p, with the physical rate conversion already given after
(55src). Combining with (55vac) proves the first line of (55sand).
If a_f is the unconditional outside mean, t_f its traced mean and
b_f^0,b_f^1 its means on E_D^c,E_D, then a_f<=t_f<=b_f^0 and
`a_f=Pr(E_D^c)b_f^0+Pr(E_D)b_f^1` imply b_f^1<=t_f. This proves the
second line and its stated covariance consequence.

Strict pair NC follows from a strictly positive initial mean derivative.
Writing beta_i=p(n-i), delta_i=qi, (55del) gives
`(d E_nu W_t/dt)|_0=(1/n)sum_(s=1)^(n-1)kappa_(s,s-1)`.
The sum equals `sum_i(beta_i-p)pi_n(i,1)-Pr(R>=2)`.
Boundary crossing and occupied-count balance give
`E[delta_I;R=0]=E[beta_I+1;R=1]` and
`E[beta_I;R=0]=Pr(R>0)`. Substituting reduces the sum to
`E[(I-np)1_(R=0)]-p Pr(R=1)`.
In the matroid representation below, C=R+K counts non-coloops and
Bernoulli integration gives
`-Cov(I,f(R))=p E[C(f(R)-f(R-1))]`.
Applying it to1_(R>0) proves the exact identity

`(d E_nu W_t/dt)|_0=(p/n)E[K 1_(R=1)]>0`.            (55strict)

For positivity, the event consisting of one boundary point followed by
n distinct successful refreshes before the next boundary point makes
the vacancy matroid rank min(|V|,1). Independently, exactly two initial
vacancies have positive probability. On their intersection R=K=1.
The mean derivative is continuous and nonnegative at all later times,
so, with mu_n(p)=E S_n, the stationary limit satisfies
`mu_(n-1)(p)>(n-1)mu_n(p)/n`.
The right conditional trace comparison and exchangeability also give
`mu_(n-1)<=E[S_n | X_tag=0]=mu_n-Var(S_n)/(n-mu_n)`.
Combining the two proves

`Var(S_n)<mu_n(1-mu_n/n)`, `Cov(X_i,X_j)<0` for i!=j.  (55NC)

This is an all-size, all-common-rate theorem. The density decrease is
at fixed refresh parameter p, so smaller graphs use the adjusted trace
sleep rate, not the original loopless lambda. The proof uses neither
preservation of a generic cone nor a stronger LC assumption. Both audits
checked the kappa input, source boundaries, reflection, trace conversion
and strictness. No strictness for every increasing outside function is
asserted.

**Arbitrary homogeneous sink ratios (proved extension; authoring AI checked).**
All conclusions (55sand)–(55NC) also hold for the holding kernel
`K_ij=1/(n+a)`, including i=j, with sink probability a/(n+a), every
finite a>0 and every common holding-model sleep rate lambda>0.
Here `p=lambda/(1+lambda)`, and the service clock has rate a.
Indeed the binomial-update sampler of
[Junge–Kaufman–Meisel, Proposition5 and Section6.2](https://arxiv.org/html/2604.04747v2#S6.SS2)
refreshes each site at rate1 and sends particles to the sink at rate
`n kappa/(1-kappa)=a`. Poissonizing the release probabilities at total
rate n+a gives precisely these independent clocks.

The maximum representation remains valid with N of rate a. Its boundary
functions obey `(a I-L_b)f_b=a f_(b+1)`, with f_(n+1)=1 and zero
extension above b-1. Reversibility gives the descending-boundary
occupation density `w_i f_b(i)/a`, w=Bin(n,p). The terminal rank mass
at k is `w_k[f_(k+1)(k)+(qk/a)f_k(k-1)]`, also the kth mass of
sup(Y-N) by summing that recursion. Time reversal therefore identifies
the same reflected queue and exchangeable occupied law, with service a.

The refined active-tag coupling uses the same common service clock at
any rate. Rectangle balance now yields the nonnegative coefficients
`kappa_(s,j)=p(n-s)pi_n(j,s-j)-a sum_(ell=0)^j pi_n(ell,s+1-ell)`.
Both (55src) and (55del) hold verbatim with these coefficients: their
component bracket is
`p h pi_n(i,r)-p h pi_n(i-1,r+1)-a pi_n(i,r+1)`.
The r=0 cancellations and the semigroup reflection signs are unchanged.
For strictness the two balances become
`E[delta_I;R=0]=E[beta_I+a;R=1]` and
`E[beta_I;R=0]=a Pr(R>0)`. Their a terms cancel in the sum of the
last kappa coefficients, leaving exactly `(p/n)E[K;R=1]` in
(55strict). The rank-one matroid cylinder now has probability
`a n! p^n/(n+a)^(n+1)>0`. Thus strictness holds at every positive
finite a; there is no uniform gap claimed at its endpoints.

Tracing k sites, m=n-k, gives K^tr_ij=1/(m+a) and killing a/(m+a),
with lambda unchanged. Hence both p and a remain fixed in the dimension
comparison. Equivalently, the loopless weighted complete graph has
K_ij=1/(n-1+a) for i!=j and killing a/(n-1+a). If its physical sleep
rate is ell, the equivalent holding rate is
`lambda=ell(n+a-1)/(n+a)`; after tracing, its loopless rate is
`ell_tr=ell(n+a-1)(m+a)/[(n+a)(m+a-1)]`.
At a=1 these formulas recover the ordinary graph. The extension asserts
neither unequal sleep rates nor arbitrary kernels, full NA or full ULC.

**Reader consequence (proved).** The singleton-versus-increasing-outside
comparison already gives product exponential-moment bounds for both signs of
the parameter, occupied and vacant cylinder bounds, and Hoeffding concentration
for nonnegative weighted counts. This needs no full NA assumption. The precise
statements and induction are in Corollary 8.2 of the
[reader](https://raw.githubusercontent.com/jfrichey/jfrichey.github.io/companions-2026-11-14/pagedocs/ai/arw-negative-association.pdf), following Proposition 8.1;
the Bernoulli means there are the actual stationary occupation probabilities.

**Every occupied-block threshold lies below its trace (proved,
independently audited).** Fix a proper block D of size k, put m=n-k,
and let K_D count its occupied sites. For every1<=ell<=k, at every
refresh parameter0<p<1 and service rate alpha>0,

`Law(X_outside | K_D>=ell)<=st nu_D`.                 (55thr)

This includes conditioning on all D occupied. Combined with (55vac),
it gives the full occupied-block/trace/vacant-block sandwich. It does
not compare the upper-threshold conditional law with the unconditional
outside law: both lie below the trace. The following exact source also
isolates that remaining NA comparison.

Write N_ell for the unnormalized stationary outside coloured law weighted
by1_(K_D>=ell), Z_ell=Pr(K_D>=ell), and N_0 for the unconditional
outside law. All outside generators use the same p,alpha. For1<=s<=m,
0<=j<s, set
`D_(s,j)=e_(j,s-j)-e_(j,s-j-1)`,
`E_(s,j)=e_(j+1,s-j-1)-e_(j,s-j)`.
For1<=t<=k define the explicitly nonnegative coefficients

`a^(t)_(s,j)=binom(n-s-t,k-t)/binom(n,k)`
` *sum_(J=j)^(j+t-1)binom(J,j)binom(s+t-J-2,s-j-1)`
`                                      *kappa_(s+t-1,J)`.        (55lift)

Here kappa uses service alpha as in the preceding extension. Put
`A^(ell)=sum_(t=ell)^k a^(t)`, `A^(0)=A^(1)`, and A_(s,-1)=0.
The exact source is

`G_m^*N_ell=sum_(s,j)[A^(ell)_(s,j)-A^(ell)_(s,j-1)]D_(s,j)`
`                              +sum_(s,j)a^(ell)_(s,j)E_(s,j)`. (55tsrc)

For ell=0 only the first sum is present, using A^(0).

Here is a verification from actual service flux. At outside counts (i,r),
put h=m-i-r and
`w_t(i,r)=E[1_(outside=(i,r),K_D=t) B_D/R_tot]`,
where B_D counts inside buffers and the ratio is zero when R_tot=0.
Exchangeability gives
`w_t(i,r)=binom(h+k-t,k-t)/binom(n,k)`
` *sum_(u=0)^(t-1)binom(i+u,u)binom(r+t-u-1,t-u-1)`
`                                      *pi_n(i+u,r+t-u)`.
Multiplication by the inside-buffer fraction explains the second binomial
factor; w_t contains no factor alpha. For any outside array v define
`K_(s,j)[v]=p(m-s+1)v(j,s-j-1)-alpha sum_(i=0)^j v(i,s-i)`.
Then

`K_(s,j)[w_t]=a^(t)_(s,j)`.                           (55flux)

To check the positive lift, put T=s+t-1 and
C=binom(n-s-t,k-t)/binom(n,k). The first term has weights
`C p(n-T)H(J)pi_n(J,T-J)`,
`H(J)=binom(J,j)binom(T-J-1,s-j-1)`.
This uses `(m-s+1)binom(n-s-t+1,k-t)`
`=(n-T)binom(n-s-t,k-t)`. The service term has coefficient
`C alpha F(I)` at pi_n(I,T+1-I), where
`F(I)=sum_(u=0)^j binom(I,u)binom(T-I,s-u)`.
Changing one inactive population label to active gives
`F(I)-F(I+1)=H(I)`, F(T)=0. Therefore F(I)=sum_(J>=I)H(J);
interchanging these finite sums proves (55flux). All indices satisfy
T<=n-1, J<T and n-s-t>=k-t. At s=m,t=k the prefactor is
1/binom(n,k), so the top boundary is retained.

For the source calculation put v_ell=sum_(t>=ell)w_t and v_0=v_1.
At outside (i,r) and inside counts (u,b,c), with R=r+b and r>0,
the full-minus-smaller generator on outside functions is
`(alpha b-p c r)D_(i+r,i)/R+p h b D_(i+r+1,i+1)/R`.
When r=0 it is zero: import and creation have the same outside output.
Swapping a designated inside hole and outside B preserves the full
stationary colour weight and gives
`E[c r/R;outside=(i,r),K_D>=ell]=(h+1)v_(ell+1)(i,r-1)`.
Negating the weighted generator discrepancy by stationarity gives
`sum_(i,r>=1)[p(h+1)v_(ell+1)(i,r-1)-alpha v_ell(i,r)]D_(i+r,i)`
`                         -sum_(i,r>=1)p h v_ell(i,r)D_(i+r+1,i+1)`.

Threshold crossings contribute two additional measures. Pair transfers
across the block boundary by swapping their B donor and empty recipient;
their outside poststates give
`sum_(i,r)p h w_ell(i,r)[e_(i+1,r)-e_(i,r)]`.
Inside service and genuine creation instead contribute
`-alpha sum_(i,r)w_ell(i,r)[e_(i+r,0)-e_(i,r)]`.
For the latter identity, on fixed outside occupied level s multiply
`alpha Pr(S_n=s+ell,R_tot>0)=p(n-s-ell+1)pi_n(s+ell-1,0)`
by `binom(k,ell)binom(m,s)ell/[binom(n,s+ell)(s+ell)]`.
This matches the total service flux to the birth flux, whose outside
state is all active. In a threshold crossing the outside **poststate**
must be used. The displayed service/creation term is an opposing term
and cannot be discarded.

Insert v_(ell+1)=v_ell-w_ell. Combining the extra w_ell term with
the transfer crossings leaves
`sum_(i,r)p h w_ell(i,r)[e_(i+1,r)-e_(i,r+1)]`.
Telescoping the service/creation term along a fixed occupied level gives
K[w_ell]E. The remaining terms give the adjacent K[v_ell] differences
on D. Linearity and (55flux) prove (55tsrc), including r=0 and both
occupied-level boundaries. At k=1, a^(1)=kappa/n; at t=1 the coefficient
is binom(n-s-1,k-1)kappa_(s,j)/binom(n,k), recovering both earlier
deletion sources with their correct factorial factors.

For increasing rank f, reflection gives E_(s,j)P_t f>=0 and
D_(s,j)P_t f>=0 decreasing in j. Summation by parts in (55tsrc) then
makes every pairing nonnegative. Thus N_ell/Z_ell increases stochastically
under the smaller semigroup to its stationary law. This is evolution
from the conditional initial law, not a time-dependent conditioning of
the full process. Exchangeability and the trace conversion prove (55thr).

For the actual NA target center the source: eta_ell=N_ell-Z_ell N_0
has zero mass, and its D coefficients are
`C^(ell)_(s,j)=A^(ell)_(s,j)-Z_ell A^(0)_(s,j)`;
the favorable E coefficients remain a^(ell). A sufficient new inequality is

`sum_(t=ell)^k a^(t)_(s,j)>=Z_ell sum_(t=1)^k a^(t)_(s,j)`. (55tail)

It says that the distribution with these kappa weights dominates K_D.
It is automatic for ell=1 and **open** for higher thresholds. If it
holds, the same semigroup sign and exponential convergence give
`Cov(1_(K_D>=ell),f(S_outside))`
`                    =-integral_0^infinity <G_m^*eta_ell,P_t f>dt<=0`.
These threshold comparisons imply full NA: conditional on disjoint block
counts, their occupied subsets are independent and uniform; permutation
averages of increasing functions increase with count. Thresholds of the
whole outside count suffice, since a smaller outside block then has an
increasing hypergeometric upper-tail probability. Condition (55tail) is
sufficient, not necessary; a negative coefficient could be compensated
by the E sector or by time integration.

The complementary lower-threshold source is
`F(sum_(t<ell)a^(t))-sum_(s,j)a^(ell)_(s,j)E_(s,j)`,
where F(A)=sum(A_(s,j)-A_(s,j-1))D_(s,j). Its two sectors oppose
one another. Monotone downward evolution is actually false: the saved
n=3,p=1/2,alpha=1 law, with k=ell=2, has
`Z_2=19/45`, `N_2(empty,A,B)=(7/45,1/5,1/15)`,
`a^(1)=1/180`, `a^(2)=1/90`, `C^(2)=11/2700`.
Direct singleton-generator substitution gives
`G_1^*N_2=(e_A-e_empty)/90`. Conditional on K_D<2, its initial outside
mean is35/52 and its derivative is1/104>0, although its stationary
limit is2/3, smaller by1/156. This refutes that trajectory shortcut,
not a static trace inequality or NA.

The independent audit checked the crossing signs, hypergeometric factors,
boundary terms, semigroup reduction and this exact obstruction.
Two additional exact checks of (55tail), using15 states at n=4,p=1/2
and21 at n=5,p=4/5, service1, found all5 and15 nontrivial centered
coefficients positive. They verify the stationary equations, binomial
active marginal, and the positive lift against direct service flux.
The preserved `trace_source_check.py` and `trace_source_checks.json`
reproduce them in under0.02s and18MiB each. The same two laws also
verify (55res) exactly on their6 and10 reduced states.
These two discriminators are evidence, not an all-size proof or a sweep.

**Integrating the favorable source sector (proved identities; authoring AI
checked and independently audited 2026-11-01, all-rate sign open).** Let G_m be the backward
three-colour queue generator and Q_m={(i,r):i,r>=0,i+r<=m}. On Q_(m-1),
define differences of a function u on Q_m by

`D u(i,r)=u(i,r+1)-u(i,r)`,
`E u(i,r)=u(i+1,r)-u(i,r+1)`.

Let L_m have diagonal `-[p(m-i)+qi+alpha]`, death
`(i,r)->(i-1,r+1)` at qi, and, for r>0, birth to(i+1,r-1)
at p(m-i) and service to(i,r-1) at alpha. At r=0, births and
services kill. Put
`B_m v(i,r)=1{r=0}[alpha v(i,0)+p(m-1-i)v(i+1,0)]`,
omitting zero-rate boundary terms. Direct row subtraction gives

`D G_m=L_m D`,
`E G_m=(G_(m-1)-I)E+B_m D`.                         (55diff)

For D, births and services at zero common buffer coalesce the two
copies, giving killing. For E, the interior birth/death differences
leave rates p(m-1-i),qi and an extra killing p+q=1. At r=0,
the extra birth term is p(m-1-i)D u(i+1,0) and the service term is
alpha D u(i,0). These are exactly B_m D, including m=1.

For an increasing rank function f, set h_f(i,r)=f(i+r+1)-f(i+r).
Initially D f=h_f and E f=0. The two identities and variation of
constants give

`integral_0^infinity D exp(tG_m)f dt=R_m h_f`,
`integral_0^infinity E exp(tG_m)f dt`
`                 =(I-G_(m-1))^-1 B_m R_m h_f`,
`R_m=(-L_m)^-1>=0`.                                 (55diff/int)

The killed discrepancy dies by its mth service event, so its mean
lifetime is at most m/alpha; the discounted ordinary queue resolvent
is finite as well. Thus the integrals and positive inverses are justified.

Fix k=n-m,ell and the centered C_(s,j) above. Regard
`a(i,r)=a^(ell)_(i+r+1,i)` as a column on Q_(m-1), and solve
`z=(I-G_(m-1)^*)^-1 a>=0`. Define

`b_s=alpha z(s-1,0)+p(m-s+1)z(s-2,0)>=0`,           (55bc)

with negative indices zero. These are the only nonzero entries of
B_m^*z, at differences D_(s,s-1). If
`v_(s,j)=(R_m h_f)(j,s-j-1)`, reflection gives
`v_(s,0)>=...>=v_(s,s-1)>=0`. Integrating the centered source and
summing its D differences by parts therefore gives the exact identity

`-Cov(1{K_D>=ell},f(S_out))`
` =sum_s [sum_(j=0)^(s-2) C_(s,j)(v_(s,j)-v_(s,j+1))`
`                         +(C_(s,s-1)+b_s)v_(s,s-1)]`. (55comp)

Consequently `C_(s,j)>=0` for j<s-1 together with
`C_(s,s-1)+b_s>=0` suffices for NA. This is weaker than (55tail).
Allowing compensation between these summands gives the exact criterion:
substitute every threshold increment h_f(i,r)=1{i+r+1=d} and sign
the resulting scalar for all k,ell,d. This reformulation is not itself
a positivity proof; its possible advantage is the explicit positive
smaller-queue contribution, which requires no historywise sign.

For m=2 the only remaining nonconstant outside test beyond OR is AND.
Write a_10,a_20,a_21 for its three a^(ell) coefficients and set
`Delta=(2p+alpha)(1+alpha)-2pq>0`. The relevant potentials solve
`(2p+alpha)v_20-2p v_21=1`,
`-q v_20+(1+alpha)v_21=1`, giving
`v_20-v_21=p/Delta`, `v_21=(1+p+alpha)/Delta`.
The discounted singleton resolvent gives

`2(1+p+alpha)b_2=p(alpha+1)(alpha+2)a_10`
` +alpha p(alpha+3)a_20`
` +alpha[1+3p+alpha(1+p)]a_21`.

Thus, for arbitrary n=k+2 and ell,

`-Cov(1{K_D>=ell},AND_out)=J/(2 Delta)`,
`J=2p C_20+2(1+p+alpha)C_21+p(alpha+1)(alpha+2)a_10`
`       +alpha p(alpha+3)a_20`
`       +alpha[1+3p+alpha(1+p)]a_21`.                (55and)

The all-rate sign J>=0 is open. In the already-saved n4,p1/2,alpha1,
k=ell=2 law, b_2=271/19840 and
`J/(2 Delta)=2671408081/215533633536`
`           =(193481/464256)^2-5/31>0`, agreeing with direct covariance.

**Occupied-particle marks (proved representation, authoring AI checked;
averaged comparison open).** Direct queue balance in(55and) leaves the
same unsigned centered kappa tail, so that scalar-rewriting route is
parked. A distinct representation follows persistent particle identities
in the labelled queue. A buffer donor retains its identity on moving to
an empty requesting site. Genuine creation occurs only when the buffer
is empty. Mark an ordered tuple of d distinct particles, and kill the
process when service removes a mark. It suffices here to take d=1,2.

Let(i,r,a) count unmarked active, unmarked buffered, and marked active
particles; put b=d-a, h=n-d-i-r, R=r+b. The killed backward generator
L_d has deaths(i-1,r+1,a) at qi and(i,r,a-1) at qa. If R>0 it has

`(i,r,a)->(i+1,r-1,a)` at `p r(1+h/R)`,
`(i,r,a)->(i,r,a+1)` at `p b(1+h/R)`,
`(i,r,a)->(i,r-1,a)` at `alpha r/R`,
`kill` at `alpha b/R`.

If R=0, replace those transitions by(i+1,0,d) at ph. Omit
zero-rate transitions. A marked donor moving to a hole leaves its
old site empty, so the unmarked counts do not change in that move.
For falling factorials(t)_j, the tuple counting measure is

`nu_d(i,r,a)=binom(d,a)(i+a)_a(r+d-a)_(d-a)pi_n(i+a,r+d-a)`.

Every existing tuple follows L_d until its death. A birth creates a
new tuple in any of d positions, with its other d-1 positions chosen
among existing particles. Thus the exact stationary source identity is

`-L_d^* nu_d=b_d`,
`b_d(i,0,d)=d p(n-i-d+1)(i+d-1)_(d-1)pi_n(i+d-1,0)`, (55marks)

with zero source elsewhere. A finite positive-probability sequence
buffers and removes every mark; the finite killed chain is transient.
Hence nu_d=(-L_d^*)^-1 b_d. Its mass is E(S_n)_d, and exchangeability
identifies its normalized remaining colour counts with those conditional
on d fixed sites occupied. This is a law-level identification; the
marked particle locations themselves move.

For d=2, AND versus an outside threshold becomes

`<b_2,(-L_2)^-1 1{i+r>=ell}>`
` <=Pr(K_out>=ell)<b_2,(-L_2)^-1 1>`.              (55marks/NA)

This is an exact reduction, not a sign proof. The autonomous active
field used in the unconditioned deletion proof is absent: unmarked
activation has rate p r(1+h/R). The source weights and marked lifetime
must be controlled together. A sufficient alternative is one-coordinate
NA in the occupied-site Palm law, since
`Cov(X_x X_y,f)=Pr(X_x=1)Cov(X_y,f|X_x=1)`
`                 +E[X_y|X_x=1]Cov(X_x,f)`;
the second term is already signed, the first remains open.

The entrance-by-entrance shortcut is false in the same saved
n4,p1/2,alpha1 law. In its18-state d=2 system, the largest all-active
entrance(2,0,2) has mean lifetime200206/67291 and mean full-outside
occupation time130080/67291. Their ratio65040/100103 exceeds the
unconditional two-site full probability193481/464256 by
10827181697/46473418368. This entrance has positive source weight.
The actual source weights at i=0,1,2 are
(3115/12896,1981/2418,147/248), and their combined occupation law
has full probability74880/193481, smaller than the unconditional one
by2671408081/89824715136. Thus averaging is essential even here.
`trace_source_check.py n4-half --particle-palm` checks both d=1,2
source equations, factorial masses, conditional rank marginals and
the backward potentials with exact rational arithmetic. The certificate
is in trace_source_checks.json; this reuses the existing model and
does not add a parameter search or support an all-size sign claim.

The natural lifetime-weighted Markov coupling is also false. Put
h_d=(-L_d)^-1 1 and M_d=E(S_n)_d. One exact sampler starts at e with
probability b_d(e)h_d(e)/M_d, jumps x->y at rate L_d(x,y)h_d(y)/h_d(x),
and stops at rate1/h_d(x). The identity L_d h_d=-1 and resolvent
conjugation give recorded law nu_d/M_d. At birth,
`b_2(i,0,2)=2(i+1)b_1(i+1,0,1)`; the second mark chooses an older
particle, but its shorter lifetime does not automatically sign the
normalized law.

In the same n4 model compare two marked active particles with one mark
and one deleted active site. Let both common sites have colours(B,0).
Their common occupied count can decrease only by service of that B.
The transformed service rates are
`gamma_2=h_2(0,0,2)/h_2(0,1,2)=43333173/43894216` and
`gamma_1=h_1(1,0,1)/h_1(1,1,1)=2400804630/2411640121`, with
`gamma_2-gamma_1=-877418425986147/105857052385440136<0`.
Thus the two-mark chain's derivative on the increasing common-OR test
is strictly larger at an equal configuration. No graphical coupling
preserving the proposed order from every naturally matched state can
work. These states are reachable from positive-source entrances.
The same checker independently verifies all lifetime equations and
four values. A coupling tailored to the averaged entrances, or a
non-Markov construction, is not excluded; none is established.
The natural graphical coupling is parked, not the Palm NA target.

**The actual kappa potential (proved, independently audited).** Put
N=n-2 and define a positive exchangeable coloured measure on N sites by
`Omega(i,r,h)=kappa_(i+r+1,i)`, i+r+h=N. It removes one B and one
hole from the kappa state; removing only the B gives the wrong weights.
If F has m-1 sites and H has k-1 sites, hypergeometric sampling gives
`a_t(s,j)=km/[n(n-1)] Omega(F has counts(j,s-j-1), K_H=t-1)`.
Thus (55tail) asks whether1+K_H under this conditional marked law
dominates the original K_D. The following actual-clock representation
proves positivity but does not yet sign that conditional comparison.

On the N-site common counts define a substochastic generator L. Put
b_i=p(n-1-i). For r>0 its transitions are
`(i,r)->(i+1,r-1)` at b_i, `(i,r)->(i-1,r+1)` at qi,
and `(i,r)->(i,r-1)` at alpha. At r=0 the birth and service events
kill the process, while active death remains. Its diagonal is
`-(b_i+qi+alpha)` in every state. Then exactly

`(I-L^*)Omega(i,r)=p q(h+1)pi_n(i,r+1)`
`                  +p alpha pi_n(i,r+2)`
`                  +p^2(h+1)pi_n(i-1,r+2)=:b(i,r)>0`, (55res)
`Omega=integral_0^infinity exp(-t)exp(tL^*)b dt`.

For the proof use the stationary last-refresh coupling for a tag y,
forcing its final coin to1 above and0 below. Its age is independent
Exp(1); the pre-refresh state has the stationary law by the Poisson
Palm decomposition. All refreshes, including null ones, must be counted.
The three discrepancy forms are the earlier ones: E has tag0 below/A
above and identical outside; C has tag B below/A above and identical
outside; X additionally has one outside B below/0 above. E cannot
return to X while the tag is frozen, and X cannot enter C.

In X remove the tagged hole and the discrepant B below. Common active
failures occur at qi. Successful refreshes at the common inactive sites
and the discrepancy site have total rate p(n-1-i); for r>0 they move
to(i+1,r-1), possibly moving the mark, and for r=0 they coalesce the
outside configurations. Service likewise moves to(i,r-1) or coalesces.
This proves the killed generator L. Its independent age discount is1.

The initial forced refresh enters X with mass
`x_init(i,r)=(h+1)pi_n(i,r+1)/n`: the tag was empty and a donor B
was present. No factor p is included, since the successful coin is
forced. At observation the C mass at outside counts(u,v) is
`(v+1)pi_n(u,v+1)/(n q)`, exactly the minus tag-buffer probability
conditioned on tagged inactivity. Exceptional donor choices create X
by service from C(i,r+1), rate alpha/(r+2), or by an empty-site
refresh from C(i-1,r+1), rate p(h+1)/(r+2). These give its two further
inflows. The Exp(1) resolvent equation for X therefore has right side
`x_init+alpha pi_n(i,r+2)/(nq)+p(h+1)pi_n(i-1,r+2)/(nq)`.
Finally X is the unique form crossing the conditional outside-buffer
CDF threshold r. The active-pattern formula and
`p(n-i)rho_i=npq binom(n-1,i)p^i q^(n-1-i)` identify
`x(i,r)=Omega(i,r)/(npq)`. Multiplication proves (55res).
At invalid indices the masses vanish; births and service still kill
at r=0, and h=0 removes inflows outside the common state space.
The positive discounted resolvent is finite because L is substochastic.
In particular every interior kappa is strictly positive at finite parameters.

Consequently every expectation comparison in (55sand) and (55thr) is
strict for a nonconstant increasing outside function. To see this,
average that function over outside permutations, obtaining an increasing
rank function f with Delta f(s)>0 for some1<=s<=m. The last D term
in the telescoped source pairs at least
`exp(-(m+alpha)t)Delta f(s)>0` with P_t f, by no refresh or service.
All its coefficients are strictly positive, including A^(ell) for the
threshold source. Thus these trajectories strictly increase to the trace.
For the vacant-block source use its negative E_(s,s-1) term instead:
no refresh and exactly one service gives a pairing at least
`alpha t exp(-(m+alpha)t)Delta f(s)>0` for t>0. Its trajectory strictly
decreases to the trace. The mixture identity then proves strict
`Cov(OR_D,f(X_outside))<0` for every nonempty proper D and every
nonconstant increasing f. This also gives a strict pair proof independent
of the quantitative matroid identity (55strict). It does not sign a
general upper-block-threshold covariance with the outside marginal.

The normalized Omega is not simply a smaller stationary queue. At
n=4,p=1/2,alpha=1 its rank masses are (778,3385,3088)/7251, so its
mean is3187/2417, below the actual two-site mean37/28 by193/67676.
Its mean active count is7210/7251, not1. Thus even unconditional
domination of the smaller queue and an unchanged product active field
are false; neither failure refutes (55tail) with its mandatory occupied
site and its different conditioning.

A stronger averaged coupling to the occupied-site Palm law is also
false. Split Omega's n-2 sites into F of size m-1 and H of size k-1,
and choose x in the original deleted block D. The proposal was
`Law_Omega(K_H|F colours)>=st Law_mu(K_(D\{x})|X_x=1)`.
It would imply (55tail), since restoring the occupied x gives the
size-biased K_D law, which dominates unconditioned K_D. But at the same
n4,p1/2,alpha1 point, with k=m=2 and F prescribed active,

`Pr_Omega(H occupied|F=A)=460/721`,
`Pr_mu(X_2=1|X_1=1)=193481/302253`,
`460/721-193481/302253=-66203/31132059<0`.          (55palm)

Indeed the three labelled marked masses for H=0,B,A are respectively
kappa_(2,1)/2, kappa_(3,1)/2, kappa_(3,2). Their ratios give the
displayed value. The actual required shifted comparison still has
positive gap `460/721-Pr_mu(K_D=2)=74057959/334728576`.
This rules out the Palm coupling, not the weaker centered tail or
compensated criteria. It uses the already-saved exact stationary law.

There is a useful quantitative consequence as alpha decreases to zero.
Fix1<=ell<=n-1 and put i=n-ell-1. The nonnegative inflows in (55res)
give
`kappa_(n-ell,n-ell-1)`
` >=[p q ell/(1+p ell+qi+alpha)]pi_n(i,1)`.
Start the stationary chain in the full occupied state(i,ell+1).
No refreshes over a fixed interval T, together with exactly ell service
points, ends at(i,1), since a B remains throughout. Therefore
`pi_n(i,1)>=pi_n(i,ell+1)exp(-(n+alpha)T)(alpha T)^ell/ell!`.
At alpha=0 occupancy cannot decrease and the autonomous active field
hits the all-active state. The unique closed class has full occupancy
and binomial active law. Finite-state continuity gives
`pi_n(i,ell+1)->binom(n,i)p^i q^(ell+1)>0`. Consequently

`kappa_(n-ell,n-ell-1)>=c_ell alpha^ell`, c_ell>0,      (55pow)

for all sufficiently small alpha. This convergence and the lower bound
are uniform for p in any fixed compact subset of(0,1). For example,
at T=1 and alpha<=1 one can take c_ell to be half the limiting state
mass times `p q ell exp(-(n+1))/[(2+p ell+qi)ell!]`.
The separate audit checked the last-refresh law, resolvent normalization,
every boundary and this unconditioned clock lower bound. This supplies
the exact leading order, not an inference from finite-alpha strictness.

There is also an exact local matroid identity. In the simultaneous
maximum representation below let E be the future clock environment,
A the independent Bernoulli(p) initial active set and V its complement.
Write S=M_E(A)=n-rank_E(V), R=|V|-rank_E(V), and let C_E(V) count
non-coloops of the restricted matroid. Its surplus is K_E(V)=C_E(V)-R.
Then for every interior index,

`kappa_(s,j)=p E[K_E(V); |A|=j,S=s]`.                 (55mat)

Indeed activate a uniformly chosen vacancy of a uniform j-set. The
result is a uniform(j+1)-set with the same independent environment,
and removal of v changes R to R-1+1_(v is a coloop). Thus
`F_(j+1)(r-1)-F_j(r-1)` equals
`E[C_E(V)1_(R=r)||A|=j]/(n-j)`.
Substitute this in the rectangle formula for kappa, which subtracts
r Pr(R=r||A|=j)/(n-j). This proves (55mat), and defines nonnegative
environment contributions kappa^E by conditioning its right side on E.
No stationary balance at a fixed environment is being asserted.

These actual environments cannot be signed separately in (55tail).
For n>=3 prescribe the first two superposed events to be services and
the next n to be successful refreshes at distinct sites. This future
clock cylinder has probability alpha^2 n! p^n/(n+alpha)^(n+2)>0.
For every A it gives M_E(A)=max(|A|,n-2), so the associated matroid
is uniform of rank2 and
`kappa^E_(s,j)=2p rho_j 1_(s=n-2)` for0<=j<=n-3.
Take k=ell=2,m=n-2,s=m in the positive lift. Then
`a_1^E(m,j)=2p rho_j/binom(n,2)` and a_2^E(m,j)=0. Its contribution
to the centered coefficient, with the original unconditional Z_2, is

`C_2^E(m,j)=-2p rho_j Z_2/binom(n,2)<0`.              (55bad)

At time zero, pairing the entire centered source with1_(S_out=m)
also gives a negative contribution: active replacements vanish and
buffer differences telescope to C_2^E(m,m-1). Similarly, signing each
positive resolvent inflow atom fails: the all-empty atom remains empty
until killing and contributes negatively to C_ell(1,0) for ell>=2.
These are obstructions to decomposing the proof into separately positive
histories, not counterexamples to the averaged centered inequality or
NA. The bounded all-rate tail attempt is parked pending a mechanism
for compensation between histories; no further parameter tests were run.

A precise logical obstruction separates the simpler OR comparisons from full
NA. The artificial exchangeable four-bit law with rank masses
`(1785,8204,12198,7852,1961)/32000` has full support and satisfies
every OR-block comparison against every increasing outside function,
but `Cov(X_1 X_2,X_3 X_4)=1/40000>0`.
For verification, averaging any increasing function under permutations
of its support gives an increasing function of its rank. It is therefore
enough to check outside thresholds. For support sizes (a,b) and outside
threshold k, the nonredundant OR covariances are

| (a,b,k) | Covariance |
|---|---:|
| (1,1,1) | -1/400 |
| (1,2,1) | -31/8000 |
| (1,2,2) | -9/8000 |
| (1,3,1) | -133/32000 |
| (1,3,2) | -53/16000 |
| (1,3,3) | -1/32000 |
| (2,2,1) | -219/40000 |
| (2,2,2) | -91/40000 |

The remaining cases follow by symmetry. These fractions follow directly
by assigning each k-subset probability equal to the kth rank mass divided
by binom(4,k); a small exact Fraction check confirmed them. Equivalently
the density relative to fair independent bits, with Z_i=2X_i-1, is
`1-(1/100)sum_(i<j)Z_i Z_j+(11/1000)sum_(i<j<k)Z_i Z_j Z_k`
`                         -(7/2000)Z_1 Z_2 Z_3 Z_4`.
This is not an ARW law, and does not refute the conjecture. It shows why
the OR consequence alone cannot close the remaining arbitrary-block gap.

A direct active/buffered split does not prove this theorem. The buffered contribution
can have the opposite sign even in the actual process: at n=2,p=1/2,
the stationary (i,r) masses in order
(0,0),(0,1),(0,2),(1,0),(1,1),(2,0) are
`(3/28,3/28,1/28,5/14,1/7,1/4)`. They give
`Cov(B_1,O_2)=3/3136>0`, `Cov(A_1,O_2)=-1/112`, and
`Cov(O_1,O_2)=-25/3136`. Convergence from empty makes the first
covariance positive at sufficiently large finite times as well. Thus
separate nonpositive active and buffered contributions cannot prove the
target; (55del) supplies a different compensation. The authoring AI checked
every coupling case, the conditioning and these fractions.
Conditioning the whole active pattern does not repair this: given I=0,
the buffer-rank probabilities are (3,3,1)/7, so the two buffer indicators
have covariance3/196>0. In the total-covariance decomposition by active
bits, the conditional-covariance term is12/3136 and the covariance of
conditional means is-37/3136. Thus conditional buffer NA is false even
though their sum is the actual negative occupied covariance.
Nor can the buffered conditional law simply be replaced or bounded above
by the trace. At n=2 the empty-outside probability given a buffered tag is
`2q(1+p)/(2+3p+4p^2)`, whereas trace emptiness is q/(1+p).
Indeed this is pi_2(0,1)/E R, using
`pi_2(0,1)=2p q^2(1+p)/[(1+2p)(1+p+p^2)]` and
`E R=pq(2+3p+4p^2)/[(1+2p)(1+p+p^2)]`.
At p=2/3, ordinary lambda3 and traced lambda4, buffered outside
occupation is21/26, exceeding trace occupation4/5 by1/130. The apparent
equality at p=1/2 is not a general identity. This rules out the simple
buffered-law closure, not the left occupied conditional comparison.

The matroid representation gives a quantitative improvement over the
mere covariance sign. With independent initial vacancies V~Bernoulli(q),
R is the nullity `|V|-r_E(V)`. Let C_E(V) count elements of V which
are not coloops in the restricted matroid. Deleting such an element
reduces nullity by one; deleting a coloop leaves it unchanged. Bernoulli
integration by parts therefore gives

`-Cov(I,R^2)=p E[(2R-1)C_E(V)]`
`             >=p[2 E R^2-g]`,

since C_E(V)>=R and both vanish when R=0. Combining with the Stein
identity yields
`(1+2p)E R^2 <= [2q(1+np)+p]g-npq(1-2p)`.
At n=2,p=1/2,g=9/28 this already gives Delta>=11/1568>0.
It does not establish the sign at all sizes.

Retain the surplus K=C_E(V)-R, and write v=npq. Applying the same
integration by parts to R gives
`p E C_E(V)=v-g`, hence `p E K=v-(1+p)g`.
For n>=2 the positive-probability environment already used above has
rank `min(|V|,1)`. On |V|>=2 it gives K=1. Thus E K>0 and the mean
bound improves strictly to

`np<E S<np+npq/(1+p)` for n>=2.

Put J=E[(R-1)K]>=0. Keeping the surplus in the square identity gives
`(1+2p)E R^2=(3+2v)g-2vq-2pJ`, and consequently

`(1+2p)Delta=[(n-1)/n](1+2p)g^2`
`              -[4(1+p+p^2)+2v]g+2v(2+p)+2pJ`.

An independent direct estimate of J is not obtained by this identity;
(55NC) now signs the complete expression. In size two J=0, because
K>0 requires R=1; the known g=9/28 at p=1/2 then recovers the
exact Delta=25/1568. This is not a new all-size sign theorem.

The surplus has an actual marked-queue description. In the forward
maximum construction let A_t be the refreshed set and let Y_full start
with all bits active, using the same future. Then
`R(V)=max_t[|V intersect A_t|-b_t]`,
`b_t=n-Y_full(t)+N_t`. The non-coloops are exactly V intersect A_tau
for the earliest maximizing prefix: deleting a vacancy reduces every
maximum iff it belongs to this smallest nested prefix. Hence K=b_tau.
Reversing the marked stationary refresh paths shows that, in the queue,
C_E(V) counts current vacancies refreshed during the present positive-R
busy period. Include the initiating death refresh; the preceding last
contact R=0 may be a left limit. Reset the marks whenever R returns to0.
At stationary (I,R,K)=(i,r,k), r>0, put l=n-i-r-k. The transitions are

| Event | Rate | New state before reset |
|---|---:|---|
| Active becomes vacant | qi | (i-1,r+1,k) |
| Marked vacancy becomes active | p(r+k) | (i+1,r-1,k) |
| Unmarked vacancy becomes active | pl | (i+1,r-1,k+1) |
| Unmarked vacancy refreshes vacant | ql | (i,r,k+1) |
| Boundary point | 1 | (i,r-1,k+1) |

If the new r is zero, replace k by zero. At r=0, a birth stays on
that boundary and a death starts (i-1,1,0). All other null refreshes
are irrelevant. These rules preserve i+r+k<=n and depend only on the
three counts, so this is a finite Markov chain. K increases within
each busy period and resets only at its end. No estimate of the
weighted occupation J sufficient for pair NC follows yet. Independent
audit checked the marked reversal, all transition rates and the
endpoint convention.

There is an exact scalar hierarchy. Let Z(z)=E z^R and
`H(z)=E[K z^(R-1);R>=1]`. The weighted Stein identity gives
`(p+qz)E[I z^R]=(1+np)Z(z)-1`; integration by parts gives

`[1+v(1-z)]Z(z)-1=p(1-z)(p+qz)[Z'(z)+H(z)]`.

Equivalently, for `F_l=E binom(R,l)` and
`h_l=E[K binom(R-1,l);R>=1]`, h_(-1)=0,
`(1+pl)F_l=pq(n-l+1)F_(l-1)-p h_(l-1)-pq h_(l-2)`.
H has nonnegative coefficients and degree at most n-2. This positivity
does not supply the required lower bound on h_1=J.

A full joint concavity invariant for this queue is impossible. If
`F(A,B,C)=E A^I B^R C^(n-I-R)` were log-concave, homogeneity would
make F^(1/n) concave. At (1,1,1) this requires
`Cov(counts)-diag(means)+means*means^T/n` to be negative semidefinite.
Its I,I entry is zero because I is exactly Bin(n,p), so its I,R entry
must vanish. But that entry equals `(1+p)g-v=-p E K<0`.
Thus the stationary joint polynomial is not log-concave, and is not
Lorentzian, for every n>=2 and 0<p<1. This does not refute LC or ULC
of the terminal rank marginal. Independent analytic review checked the
authoring AI queue, nullity and concavity arguments; the authoring AI checked the
reviewing AI's strict surplus and factorial-moment consequences.

Even the full scalar hierarchy omits essential dynamics. In an artificial
matroid environment draw Theta~Beta(1,1/p), then independently declare
each ground element a loop with probability Theta and a coloop otherwise.
Sample V~Bernoulli(q) independently. Here R counts selected loops,
C_E(V)=R, and K=0. The conditional record PGFs satisfy
`G_i=1-ell_i+z ell_i G_(i+1)`,
`ell_i=p(n-i)/[1+p(n-i)]`; this follows directly from the beta-binomial
probabilities. The positive first-subdiagonal ladder satisfies the
weighted phase-flow identity, so all scalar Stein and nullity identities
above hold. Yet conditional on Theta, M is Bin(n,p+q Theta), giving
strictly positive exchangeable pair covariance
`p^2q^2/[(1+p)^2(1+2p)]`.
This is not an actual future. It fails the Riccati subdiagonal equation
at i>=1, and the actual quadratic occupation identity: its u_i is
affine, so E[M K_2(I)]=0, while
`E K_2(M)=n(n-1)q^2 E Theta^2>0`.
Thus this obstruction concerns the scalar relaxation, not the conjunction
of all actual constraints or the desired pair inequality.

**Parameter coupling and the published large-size pair inequality
(proved coupling; source-derived asymptotic consequence).** The maximum representation
also permits a quantitative comparison of two refresh parameters p,p'.
Use the same clocks and uniform coins, including the initial coins, and
put delta=|p'-p|. The count of discrepant coordinates at every fixed
refresh epoch has law Bin(n,delta), even conditional on all clock times
and site labels: the last coin used at each site is still an independent
uniform variable. Count all refresh clock rings, including unchanged bits;
conditioning on actual state-change times would not justify this argument.
Both maxima are attained before the nth boundary
arrival. There are n^2 refreshes before that time in expectation, so the
union bound gives, for every integer r>=1,

`Pr(|S_(n,p')-S_(n,p)|>=r)`
` <=min(1,(1+n^2)Pr(Bin(n,delta)>=r))`.

Here the two S variables are coupled through their maximum realizations;
no pathwise assertion about the two descending stopping times is needed.
If n delta<=C, the elementary bound
`Pr(Bin(n,delta)>=r)<=(e C/r)^r` yields
`||S_(n,p')-S_(n,p)||_2=O_C(log(n+1))` by summing the tail.
For example, above max(2eC,6 log_2(n+1)) the tail is at most
`(1+n^2)2^(-r)`.

Fix a loopless common sleep rate lambda>0 and put p*=lambda/(1+lambda).
The corresponding refresh parameter is
`p_n=n lambda/(n+1+n lambda)=n p*/(n+1-p*)`, so
`n|p_n-p*|<=p*(1-p*)<=1/4`.
[Junge--Kaufman--Meisel, Theorem 3 and Section 1.2](https://arxiv.org/html/2604.04747v2)
give, at fixed p* and sink probability1/(n+1),
`Var(S_(n,p*))~p*(1-p*) pi^2 n/(6 log n)` and
`E S_(n,p*)/n->p*`. Their Remark7 explicitly distinguishes the
moving-parameter issue when removing self-loops. The preceding coupling
handles it without assuming a uniform version of their theorem:
`O(log n)=o(sqrt(n/log n))`, and the difference of standard deviations
is at most the L^2 norm of the difference. Thus the ordinary loopless
stationary law has the same variance asymptotic and, for distinct sites,

`Cov(X_i,X_j)~ -p*(1-p*)/n`.

This is an asymptotic consequence of the cited theorem with an elementary
parameter comparison, not a new all-size sign theorem or an NA result.
For every fixed lambda it gives strict pair NC for all sufficiently
large n; neither an effective threshold nor uniformity in lambda follows.
An independent audit checked the conditional coin law, second-moment
bound, loop conversion, source hypotheses and variance constant.



## Auxiliary complete-graph survival inequalities

**The top infinite-horizon survival gap is now proved.** For every
m>=2 and 0<p<1, with boundary rate one,

`a_(m-1)^2>a_(m-2)a_m=a_(m-2)`.                     (55g)

This is an auxiliary survival inequality, distinct from stationary
top ULC: the latter concerns h_k=a_k+qk a_(k-1) and also uses
a_(m-3). The next survival gap is proved below; the lower gaps and
finite-time assertion remain open.

Write F,G for the ordinary transforms of the last two passages,
evaluated at z=1 unless indicated. The queue representation gives
`a_(m-1)=1-F`, `a_(m-2)=1-G+FG'`, so its gap is
`H=G-2F+F^2-FG'`. In the return measure from (55e), now at z=1,
put u=E theta=F and v=E theta^2=FG. The pushforward by theta has
density
`[p^(m+1)/R_m(1)](p+q theta)^(-m-2)`, 0<theta<1.
It is a positive mixture of normalized densities exp(-L theta) on
[0,1]: expand the negative power by its elementary gamma integral
and absorb each component's normalizer into the mixing measure.
All tilts L are finite almost surely.

The passage continued fraction gives
`C=-G'-G^2/(2p)>=0`. Equivalently, the killed reversible resolvent
has `G(z)=2p sum_j w_j/(z+lambda_j)` with w_j>=0, sum w_j=1;
C is 2p times the variance of 1/(1+lambda_j). Direct substitution yields

`uH=K(u,v)+q v^2/(2p)+u^2 C`,
`K(u,v)=v+v^2/2-u^2(2-u)`.

The needed moment lemma holds for every mixture just described:
`K(u,v)>0`. For a single component with tilt L>0,
`u_L=1/L-1/(exp(L)-1)` and
`v_L=2/L^2-(1+2/L)/(exp(L)-1)`.
Set `phi(u)=sqrt(1+4u^2-2u^3)-1`. On 0<=u<=1/2,
`phi''(u)=(4-6u-8u^3+3u^4)/(1+4u^2-2u^3)^(3/2)>0`:
the numerator decreases to 3/16. The scalar inequality
v_L>phi(u_L) therefore passes strictly to arbitrary mixtures by
Jensen. At L=0 it holds since K(1/2,1/3)=1/72.

Here is a fixed positive integral certificate for the scalar inequality.
Define `N(L)=2 L^4(exp(L)-1)^3 K(u_L,v_L)` and let
`N_0=N`, `N_(j+1)=(partial_L-c_j)N_j`, with
`(c_0,...,c_9)=(0,0,0,0,0,1,1,1,1,1)`.
Exact differentiation gives
`(N_0(0),...,N_9(0))=(0,0,0,0,0,0,0,140,700,4088)` and

`N_10=exp(2L)[648L^5+6524L^4+21968L^3+46864L^2`
` +50528L+19040+(15552L+95904)(exp(L)-sum_(j=0)^4 L^j/j!)]>0`.

The exponential remainder is the integral of
`exp(s)(L-s)^4/4!` over 0<=s<=L. Each preceding first-order equation
has a positive integrating factor and nonnegative initial value,
so backward integration proves N(L)>0 for L>0. This completes
the mixture lemma and (55g). At p=1 the gap is 1/36.

An independent certificate uses only three neighboring first-step
equations. For m>=4 put h=m-2, and let J be the preceding passage
transform at one. Then
`G=2p/[1+2p+hq(1-J)]`,
`F=p/[1+p+(h+1)q(1-G)]`, and
`-G'>=G^2[1+hq J^2/(3p)]/(2p)`.
The next preceding transform is at least
`K_0=4p/[1+4p+q(m-4)]`; hence
`ell=3p/[1+3p+q(m-3)(1-K_0)]<=J<=v=3p/(1+3p)`.
Substitute `m=w+4`, `p=t/(1+t)`, `J=(ell+zv)/(1+z)` in the
resulting lower bound for H. Its positive-denominator numerator has
389 positive integer coefficients in t,w,z, of total degree 17.
Their counts by z-degree are 95,98,98,98, and the coefficient of t
at w=z=0 is positive. The cases m=2,3 have direct positive rational
formulas after the same substitution for p. The selected
`complete_graph.py --survival-top-lc` calculation preserves all 389
coefficients and independently checks the queue normalization against
killed resolvents at m=4,p=2/5 and m=6,p=3/7. The moment proof and
the scalar-bound certificate agree without an all-size numerical sweep.

The return representation also retains the full hierarchy. For 0<=s<=m,
`M_s(z)=product_(i=m-s)^(m-1)phi_i(z)=E_(mu_z)theta^s`, where
`mu_z(dtheta)` is proportional to
`(1-theta)^(z-1)(p+q theta)^(-m-z-1)dtheta`.
This follows by renewing the all-ones transition probability from s
specified zeros. Higher integral moments exist but are not additional
physical passages. Put
`tau_p(theta)=log[(p+q theta)/(p(1-theta))]`.
Then mu_z is the exponential tilt exp(-(z-1)tau_p) of mu_1, so
`M_s'=-Cov(theta^s,tau_p)` and
`M_s''=Cov(theta^s,(tau_p-E tau_p)^2)`.

A useful comparison survives every required size bias. Let nu_(s,z)
be mu_z tilted by theta^s, 0<=s<m, and let
`gamma_(s,z,L)(dtheta)` be proportional to
`theta^s(1-theta)^(z-1)exp(-L theta)dtheta`.
There is a unique finite L>=0 with the same mean eta as nu_(s,z):
the reference mean decreases continuously from (s+1)/(s+z+1) to zero,
and the actual decreasing tilt puts eta strictly between these limits.
The log density ratio is
`constant-(m+z+1)log(p+q theta)+L theta`, strictly convex.
Equal mass and mean therefore force the density difference to have
the signs +,-,+ at two interior crossings. Subtracting a convex
function's secant through them proves strict convex-order domination
for strictly convex integrable tests.

In particular, psi(theta)=(theta-eta)tau_p(theta) is strictly convex:

`psi''(theta)=[(1-theta)(2p(1-eta)+eta)+theta(1-eta)]`
`                         /[(1-theta)^2(p+q theta)^2]>0`.

Since eta=phi_(m-s-1)(z), this proves the actual derivative bound

`-phi_(m-s-1)'(z)>Cov_(gamma_(s,z,L))(theta,tau_p(theta))`. (55h)

At z=1 the reference is a truncated gamma density. Its moments use
`Z_s(L)=s! L^(-s-1)[1-exp(-L)sum_(j=0)^s L^j/j!]`, with
Z_s(0)=1/(s+1). Under p=1-c/m, fixed c>0 and fixed s, both actual
and matched reference densities tend to theta^s exp(-c theta), and
tau_p tends to -log(1-theta). Dominated logarithmic moments show that
(55h) is sharp in this limit. The weaker derivative floors below
lose information precisely in this scaling.

For the next gap let F,G,J be the last three passage transforms at one,
`y=-J'`, `zeta=J''/2`, `d=-G'`, and A,B,C the last three nontrivial
survival values. The independent count failures give exactly
`A=1-F`, `B=1-G-Fd`, `C=1-J-yG-ydF-zeta GF`.
If a=M_1,b=M_2,c=M_3 and kappa=(1+a')/a, these simplify to
`B=1+b'-kappa b` and
`C=1-c''/2+kappa c'-(1+kappa b'-b''/2)c/b`.
Thus the first three return moments and their joint derivatives
suffice; replacing these derivatives by unrelated bounds is lossy.

For example set beta=q(m-2), d_0=G^2/(2p), d_1=beta d_0,
`B_0=1-G-Fd_0`. Then
`B^2-AC=c_0+c_1 y+c_2 y^2+AGF e`, where
`c_0=B_0^2-A(1-J)`,
`c_1=-2B_0 Fd_1+A(G+Fd_0)`,
`c_2=F^2d_1^2+AF(G/J+d_1)>0`, and e=zeta-y^2/J>=0.
The floor y>=J^2/(3p) alone is insufficient even at actual F,G,J.
Along p=1-2/m its relaxed constant term tends to
`P(exp(2))/[288(exp(2)-5)(exp(2)-3)^2(exp(2)-1)^2]<0`, with
`P(E)=387E^5-8613E^4+76620E^3-344404E^2+776137E-694687`.
Indeed P' is positive on [7,infinity), P(37/5)<0, and the exponential
series gives 7<exp(2)<37/5. This refutes that relaxation only.
Retaining the preceding passage K, with alpha=q(m-3), gives
`y>=J^2[1+alpha K^2/(4p)]/(3p)` and
`e>=alpha J^2 K^3/(48p^3)`. Bound (55h) supplies a different,
asymptotically sharp lower bound for y.

**The next infinite-horizon survival gap is also proved.** For every
integer m>=3 and 0<p<1, with boundary rate one,

`a_(m-2)^2>a_(m-3)a_(m-1)`.                         (55i)

Use the preceding notation, and call the two fourth-passage floors
y_0 and e_0. Differentiating the actual first-step equations gives
`y=J^2(1-alpha K')/(3p)` and
`e=alpha J^2 K''/(6p)`. In turn
`-K'>=K^2/(4p)` and `K''/2>=K'^2/K>=K^3/(16p^2)`.
These statements use only the preceding transform's decreasing first
derivative and nonnegative second derivative; they remain valid when
K is the initial exponential passage. Hence s=y-y_0>=0 and e-e_0>=0,
and the exact gap decomposes as

`B^2-AC=R_0+R_1 s+c_2 s^2+AGF(e-e_0)`,
`R_0=c_0+c_1 y_0+c_2 y_0^2+AGF e_0`,
`R_1=c_1+2c_2 y_0`.

The linear coefficient is positive analytically. The last-passage
equation gives `F=p/[2+beta-(beta+q)G]`, hence

`c_1=G(1-F-G)+Fd_0(5-F-2qG)+2F^2 beta d_0^2`.

Since F<=1/2, G<=2/3, d_0>=G^2/2 and
5-F-2qG>=19/6>3, this is at least
`G[(1-F)-G(1-3F/2)]>=G/3>0`. Thus R_1>0.
It remains to sign the constant R_0 using the actual range of K.

First write m=w_0+4, p=t/(1+t), and
`K=(ell+u h)/(1+u)`, where
`ell=4t/(1+5t+w_0)`, `h=4t/(1+5t)`.
Substitution in R_0 gives N(t,w_0,u)/D(t,w_0,u).
Both N and D have degree five in u. The denominator is a product
of three positive-coefficient polynomials, with constant terms
2,6,24 and powers 1,2,2. N has 1,285 terms, of which 65 have
negative coefficients: this wider interval alone does not suffice.

For m=w+6 retain two more actual first-step equations. Their lower
values are
`k_6=6t/(1+7t+w)`,
`k_5=5t/[1+6t+(w+1)(1-k_6)]`,
`k_4=4t/[1+5t+(w+2)(1-k_5)]`.
Monotonicity of each first-step fraction gives k_4<=K<h. At m=6
all three lower values are exact, so this includes the endpoint.
Writing `K=(k_4+v h)/(1+v)` gives

`w_0=w+2`, `u=u_0+(1+u_0)v`,
`u_0=alpha_0/beta_0`
` =5t(5t+1)(7t+w+1)/[(5t+w+3)(7t^2+2tw+9t+w^2+3w+2)]`.

Here t>0,w,v>=0 and beta_0>0. The homogenized numerator
`Ntilde=beta_0^5 N(t,w+2,[alpha_0+(alpha_0+beta_0)v]/beta_0)`
has exactly 3,895 positive integer coefficients, individual degrees
(36,30,5) and total degree 41. Its coefficient of t at w=v=0 is
5,079,158,784,000, proving strict positivity throughout the domain.
The denominator remains positive under this substitution. The ordered
coefficient-list SHA256 is
`07bf0de0d0a198ac1110368c95dde732bdf3a083d405a92be77da3a2d29483a4`.
Thus R_0>0 for all m>=6, including m=6 itself.

For m=4, K=4p/(1+4p) and the same numerator is N(t,0,0):
its 21 nonzero coefficients are positive, beginning with 17280t.
For m=5 the initial birth-5p passage is exponential and
`K=2t(6t+1)/(15t^2+6t+1)`. Put
`alpha_5=5t(5t+1)`, `beta_5=(t+1)(5t+2)`.
The numerator `beta_5^5 N(t,1,alpha_5/beta_5)` has 31 positive
coefficients, beginning with 318504960t. The derivative floors need
not be equalities in this case. Finally, the exact m=3 gap is

`t(t+1)^4(424t^6+2064t^5+3636t^4+3265t^3+1626t^2+432t+48)`
` /[36(2t+1)^2(4t+1)(2t^2+2t+1)^2(6t^2+4t+1)^2]>0`.

This covers every integer size. At p=1 the gap is 53/10368.
The selected `complete_graph.py --survival-next-lc` calculation
reconstructs the formulas and preserves the exact coefficient lists;
it also verifies all three survival values against independent killed
resolvents at the base sizes and both a refined-interval endpoint and
an interior point. A separate polynomial-ring implementation reproduced
the entire 3,895-coefficient list, denominators and base cases. The
1,170 positive coefficients of R_1 provide an additional algebraic check,
although the displayed analytic bound already proves its sign.
The third survival gap is addressed below. Finite-time LC and full
stationary ULC remain open.

**Reduction of the third survival gap, with a quantitative obstruction.**
For m>=4 retain the same actual transforms F,G,J,K and put
`A=a_(m-1)`, `B=a_(m-2)`, `C=a_(m-3)`, `D=a_(m-4)`.
Use xi=-K', v=K''/2 and w=-K'''/6, at argument one. The independent
Poisson passage-count queue has first-failure prefix counts 1,1,2,5
at lengths one through four. Summing their disjoint probabilities gives

`D=1-K-xi(1-C)-v[J(1-B)+yGF]-wJGF`.                 (55j)

Set alpha=J^2/(3p), r=q(m-3)alpha, beta=GF r and Q=JGF. The
first-step identities are
`y=alpha+r xi`, `d=d_0+d_1 y`, `zeta=y^2/J+r v`, with the previous
d_0,d_1. If c=q(m-4) and L is the fifth passage transform, then

`xi=K^2(1-cL')/(4p)`,
`v=xi^2/K+epsilon`, `w=xi^3/K^2+2xi epsilon/K+eta`,
`epsilon=cK^2 L''/(8p)>=0`, `eta=cK^2(-L'''/6)/(4p)>=0`.

At m=4 omit L and its zero-rate terms; then
xi=xi_0=K^2/(4p) and epsilon=eta=0 exactly. For formal xi>=xi_0 define

`B(xi)=1-G-F(d_0+d_1 y)`,
`U(xi)=1-J-Gy-Fy(d_0+d_1 y)-GF[y^2/J+r xi^2/K]`,
`D_0(xi)=1-K-xi(1-U)-xi^2[J(1-B)+yGF]/K-xi^3 Q/K^2`,
`H_0(xi)=U(xi)^2-B(xi)D_0(xi)`,
`Gamma=J(1-B)+GF(2y-alpha+2xi J/K)`.

Thus H_0 is an explicit quartic, C=U-beta epsilon, and the actual
third gap has the exact expansion

`C^2-BD=H_0(xi)+epsilon[B Gamma-2beta U]`
`                       +beta^2 epsilon^2+BQ eta`.

The stronger certificate in (55i) justifies lowering both higher
slacks to zero, not merely bounding the preceding actual gap.
With its R_0>0, R_1>0 and c_2>0, all evaluated at the actual zeroth
transforms, direct substitution gives

`B(xi)^2-AU(xi)=R_0+R_1 r(xi-xi_0)+c_2 r^2(xi-xi_0)^2`
`                           +AGF r(xi^2-xi_0^2)/K>0`.

This holds for every formal xi>=xi_0. If B(xi)>0, then
`partial_epsilon(C^2-BD)=B Gamma-2beta C>0`:
use C<=B^2/A, B<=1-G, Fq(m-3)(1-G)<p and 2J<=3A to get
`2beta C<BGJ`, while
`B Gamma-BGJ=BJF d+BGF(2y-alpha+2xi J/K)>0`.
Also `partial_eta(C^2-BD)=BQ>0`. Actual B is positive, hence

`actual third gap >= H_0(xi_actual)`,

with equality at m=4 and strict comparison at m>=5. This reduces
the target but does not yet prove its sign.

The actual derivative has a finite upper bound. Reversibility of the
ordinary chain killed on upward exit writes its passage density as
`4p sum_j c_j exp(-lambda_j s)`, with c_j>=0, sum c_j=1 and
lambda_j>0. The weights omega_j=4p c_j/lambda_j sum to one.
Putting R_j=1/(1+lambda_j) gives
`(K,xi,v,w)=E_omega[(1-R),R(1-R),R^2(1-R),R^3(1-R)]`.
Variance and Cauchy--Schwarz therefore yield

`xi_0<=xi<=K(1-K)`,
`0<=epsilon<=xi[1-xi/(K(1-K))]`, `eta>=epsilon^2/xi`.

An unbounded shift of xi cannot give a positive-coefficient proof:
at each fixed m, as q=1-p decreases to zero, the quartic's leading
coefficient is `-25(m-3)(m+1)q^2/12288+O(q^3)<0`.
Even positivity on the whole displayed bounded interval is false.
Along actual parameters p=1-2/m, put E=exp(2). The zeroth transforms
tend to
`F=(E-3)/(2(E-1))`, `G=(E-5)/(E-3)`,
`J=(3E-19)/(2(E-5))`, `K=6(E-7)/(3E-19)`.
At the relaxed lower endpoint xi_0=K^2/4,

`H_0(xi_0)=N(E)/[82944(E-5)^2(E-3)^2(E-1)^2(3E-19)]<0`,
`N(E)=11298771E^7-479855421E^6+8727946911E^5-88141154049E^4`
` +533853269433E^3-1939918623255E^2+3917702255253E-3393875313259`.

Indeed 369/50<E<37/5 by the exponential series. Substituting
`E=(369/50+(37/5)u)/(1+u)` and clearing `50^7(1+u)^7` gives
eight strictly negative integer coefficients. At the same endpoint
B and U are positive, and 8/11<K<4/5, compatible with the sixth-passage
interval in (55i). Continuity gives the same obstruction for all
sufficiently large finite m. No actual derivative has been substituted
here, so this is not a survival-LC counterexample.

The missing quantitative information is available in principle.
For m>=5, c>0 and `L=[1+4p+c-4p/K]/c`. Its first derivative bound gives
`xi>=xi_5=K^2[1+cL^2/(5p)]/(4p)>xi_0`, excluding the bad endpoint.
The higher bounds
`epsilon>=cK^2 L^3/[4p(5p)^2]` and
`eta>=cK^2 L^4/[4p(5p)^3]` can also be retained using the proved
monotonicity. These floors follow by expanding the fifth first-step
reciprocal at argument 1-h: its denominator is
`(5p/L)[1-Ah-sum_(j>=2) B_j h^j]`, with A>=L/(5p), B_j>=0.
Thus its j-th normalized jet is at least L^(j+1)/(5p)^j.
Write epsilon_5,eta_5 for the two displayed floors and
`P(X)=H(X,epsilon_5,eta_5)`. Then

`actual third gap >= P(xi_actual)`,
`xi_5<=xi_actual<=xi_+=K(1-K)`,
`xi_+-xi_5=cK^2[(1-L)-L^2/(5p)]/(4p)>0`.

The final bracket is at least 1/(1+5p)^2, since
0<L<=5p/(1+5p). At m=5 all three jet floors are exact.

**One simultaneous-limit regime is proved.** Along p=1-2/m, the fifth
transform tends to `L=(15E-109)/[6(E-7)]`, E=exp(2). The first-jet-only
bound H_0(xi_5) is still strictly negative. However, retaining both
higher floors makes P positive on the entire interval [xi_5,xi_+].
More precisely,

`(1+u)^4 P((xi_5+u xi_+)/(1+u))=sum_(j=0)^4 N_j(E)u^j/D(E)`,
`D(E)=51840000(E-7)^2(E-5)^2(E-3)^2(E-1)^2(3E-19)>0`.

Each N_j has degree nine. The exponential series gives
`ell=691283/93555<E<134800213/18243225=h`:
retain terms through degree 12 and bound the subsequent tail by a
geometric series of ratio 1/7. Substitution E=(ell+h z)/(1+z), followed
by positive denominator clearing, makes all ten coefficients of each
N_j strictly positive. The certificate retains all 50 signs and their
original coefficient rows, whose SHA256 is
`b3809f54abee37e020eefc62d8303ed98a99c38c658cc35dddc5ca83d20c824e`.
The same interval proves B(xi_+)>0; B decreases in X. Hence all five
Bernstein coefficients are strictly positive on this interval.
Dominated convergence of the actual return-density moments and
continuity of these rational functions prove the actual third gap
for every sufficiently large m along p=1-2/m. No threshold is claimed.

The remaining all-size candidate is explicit. For m=w+8 and p=t/(1+t),
start k_9=0 and recursively set
`k_r=rt/[1+(r+1)t+(w+8-r)(1-k_(r+1))]`, r=8,7,6,5.
Then `k_5<=L_actual<5t/(1+6t)`. Substituting
`L=[k_5+u*5t/(1+6t)]/(1+u)` into the five Bernstein coefficients of
P gives rational inequalities in t>0,w,u>=0, but the enlarged domain is
too weak. Along t->infinity,w=3t,u=0 its lower endpoint tends to
L=115/153, and the reconstructed transforms tend to
(F,G,J,K)=(1999/7117,2878/5997,879/1439,204/293). Here
xi_5=41792/257547<xi_+=18156/85849, yet
`P(xi_5)=-31670939693550051419413/29841111934881500548931031<0`.
This exact relaxed-endpoint obstruction does not substitute the actual
fifth transform. Stronger value or derivative bounds are needed. The actual smaller bases m=4,5,6,7
are proved: after p=t/(1+t), the exact third-gap numerators have
respectively 19,26,33,40 strictly positive nonzero coefficients, with
positive denominator coefficients. Their four survival values also match
independent boundary resolvents at p=2/5. The certificate retains every
coefficient and its hash. The queue normalization,
jet identities, zero-slack reduction and moment bounds passed independent
audits. The selected `--survival-third-reduction` and
`--survival-third-fifth-floors` certificates retain the obstructions and
the scoped positive result; neither claims all-size third-gap LC.

The derivative floors themselves also need refinement. At the actual
zeroth transforms along p=1-3/m, the old endpoint
H(xi_5,epsilon_5,eta_5) is strictly negative. Thus refining only the
zeroth-value interval cannot make that particular endpoint positive.
For m>=6 put d=q(m-5), M=P_6(1) and
`L=5p/[1+5p+d(1-M)]`. The sixth-passage jet floors give

`s=M^2/(6p)`, `e=M^3/(6p)^2`, `h=M^4/(6p)^3`,
`ell_1=L^2(1+ds)/(5p)`, `e_L=dL^2 e/(5p)`, `h_L=dL^2 h/(5p)`,
`ell_2=ell_1^2/L+e_L`,
`ell_3=ell_1^3/L^2+2ell_1 e_L/L+h_L`.

The positive reciprocal first-step expansion proves that these bound
the actual first three normalized L jets from below. Therefore

`xi_6=K^2(1+c ell_1)/(4p)`,
`epsilon_6=cK^2 ell_2/(4p)`, `eta_6=cK^2 ell_3/(4p)`

are valid floors, exact at m=6. The existing higher-slack monotonicity
gives `actual third gap>=P_6(xi_actual)`, where
`P_6(X)=H(X,epsilon_6,eta_6)`. Moreover

`xi_+-xi_6=cK^2[(1-L)-ell_1]/(4p)>0`,
`(1-L)-ell_1=L[(1-L)+d((1-M)-LM^2/(6p))]/(5p)`.

The sign follows from M<=6p/(1+6p), which gives
M^2/(6p)<=M(1-M), and 0<L<1.

These stronger floors prove a second simultaneous-limit regime.
At p=1-3/m, put E=exp(3). The limiting values are
`F=(E-4)/[3(E-1)]`, `G=(2E-17)/[3(E-4)]`,
`J=2(E-13)/(2E-17)`, `K=(8E-131)/[6(E-13)]`,
`L=8(5E-92)/[3(8E-131)]`, `M=(80E-1553)/[8(5E-92)]`.
They are successive moment ratios of exp(-3 theta) on [0,1].
Let ell=sum_(j=0)^18 3^j/j! and
h=ell+(3^19/19!)/(1-3/20); then ell<E<h by positive series and the
geometric tail bound. For the old floor package the constant Bernstein
numerator has ten strictly negative coefficients after
E=(ell+h z)/(1+z). For the sixth package all five degree-nine
Bernstein numerators have ten strictly positive coefficients.
Their common positive denominator is
`10883911680000(E-13)(E-4)^2(E-1)^2(2E-17)^2(8E-131)^2`.
The original 50 coefficient rows have SHA256
`3e654515dffec88bd65821f82cda72faf7a56dfc046475d8ff222c6f161540dc`;
all transformed signs and clearing factors are in the selected certificate.
Also `B(xi_+)=(23E^2+23E+1889)/[81(E-4)(E-1)]>0`.
Thus P_6 is positive on its whole interval. Dominated convergence of
the actual moment ratios and strictness prove the third gap for all
sufficiently large m along p=1-3/m, with no explicit threshold.

**The third survival gap at every size and rate.** The sixth floors
above, together with the depth-eight value interval, now prove

`a_(m-3)^2>a_(m-4)a_(m-2)` for every m>=4 and 0<p<1.       (55k)

This is the actual infinite-horizon survival diagonal with boundary
rate one. For m>=8 put m=w+8 and p=t/(1+t), so w>=0 and t>0. Define

`b_8=1+9t+w`,
`b_7=(2+8t+w)b_8-8t(w+1)`,
`b_6=(3+7t+w)b_7-7t(w+2)b_8`,
`M_-=6t b_7/b_6`, `M_+=6t/(1+7t)`.

The increasing first-step fractions at passages eight, seven and six,
started with the ninth transform replaced by zero, give
`M_-<=M_actual<M_+`. At m=8 the lower value is exact. Hence every
actual value has the representation

`M=(M_-+u M_+)/(1+u)=N_M/D_M`, u>=0,
`N_M=6t[(1+7t)b_7+u b_6]`, `D_M=(1+u)(1+7t)b_6`.

All displayed polynomial denominators are positive. Recover L,K,J,G,F
by their exact first-step fractions. The previous width identity gives
`xi_6<xi_+=K(1-K)`, and the actual first jet lies in this interval.
At the actual first jet the actual B is positive; the stronger preceding
gap and (55j) therefore justify lowering epsilon and eta to their sixth
floors. Thus the actual gap is at least P_6(xi_actual). No positivity
of B at other formal first jets is needed for this comparison.

Define the unnormalized Bernstein coefficients C_j by

`(1+v)^4 P_6((xi_6+v xi_+)/(1+v))=sum_(j=0)^4 C_j v^j`.

Their positive denominators and numerator signs have a finite exact
certificate. Before substituting M, use affine continuants

`kappa_6=1`, `kappa_5=4+6t+w-(w+3)M`,
`kappa_r=[1+(r+1)t+w+8-r]kappa_(r+1)`
`        -(w+8-r)(r+1)t kappa_(r+2)`, r=4,3,2,1.

Then `P_r(1)=rt kappa_(r+1)/kappa_r` for r=1,...,5. The reduced
denominator of C_0 is `9 kappa_1^2 kappa_2^2 kappa_3^2 kappa_4`;
for C_1,...,C_4 it has the additional factor kappa_5^2. The generic
numerator degrees in M are respectively 7 and 9, equal to these total
denominator powers. Homogenize each numerator with D_M to that degree.
The substituted continuants
`kappa_tilde_r=D_M kappa_r(t,w,N_M/D_M)` have positive nonzero
coefficients and constants 40320,5040,720,120,24, respectively.
The corresponding products are consequently positive denominators;
all powers of D_M cancel exactly.

The five homogenized numerators have respectively
11,592,20,250,20,250,20,250,20,250 strictly positive integer coefficients.
Their degrees in (t,w,u) are (54,45,7) for j=0 and (64,53,9) otherwise.
Put `A_*=1981728150031171584000000`. Their coefficients of t*w^0*u^0
are A_* for j=0 and `576 binom(4,j) A_*` for j=1,...,4.
These anchors give strict positivity even at w=u=0. The ordered
coefficient-list SHA256 values are:

| j | SHA256 |
|---|---|
| 0 | `11469380c0c03823b735c6ab684af669f0b2bffd82553e0165c50a776c6640b7` |
| 1 | `648aadf78a1f13abab77a3f1dc5a4ccfdcc10f68bf037d04a2a3c2f9f1387c4b` |
| 2 | `cacb4a482f80683523c1c46e43d32dc18b51064a8079f1c1be45aa2a3002c854` |
| 3 | `7ccb56fbd46cfb5f13be810106cc85476174d65d08e29e6b36214e361db35152` |
| 4 | `ed04e9c10948e1018a7ac62f1e7c77fb968e3fc426aabdbfe86a2355638243a8` |

Each digest uses ascending rows [t-degree,w-degree,u-degree,integer]
serialized as compact JSON. Positivity of C_0,...,C_4 proves positivity
on the closed interval: finite v covers [xi_6,xi_+), while C_4=P_6(xi_+).
Together with the already proved actual m=4,5,6,7 bases this proves (55k).

The selected `complete_graph.py --survival-third-lc` reconstructs all
92,592 coefficients, the four bases and their independent boundary-resolvent
checks. It preserves compact statistics, hashes and strict anchors rather
than the large raw coefficient lists. It checks the queue quartic against
direct polynomial arithmetic, normalized jets against a truncated continued
fraction, and homogeneous substitution against an independent Horner form.
A separate checker, importing none of this construction, recovered every
coefficient and denominator from ordinary series division, the nine
ungrouped queue-failure terms and Taylor-to-Bernstein conversion.
The parameter coverage, slack comparison and endpoints passed independent
analytic audit. This proves the first three survival gaps at all sizes
and rates; the fourth is proved in (55q) below. Finite-time LC and full
stationary ULC/NA remain open.

**Responses of the two highest passage slacks.** There is an all-depth
comparison behind the preceding reductions. Fix actual Ehrenfest zeroth
transforms f_i=P_i(1), write d_i=-P_i'(1), and put delta_i=(m-i)q.
The return-moment representation gives

`f_i=I_i/I_(i-1)`, `I_j=integral_0^1 theta^j(p+q theta)^(-m-2)dtheta`.

Integrating the derivative of `theta^i/(p+q theta)^(m+1)` proves
`ip I_(i-1)-(m+1-i)q I_i=1`. In particular
`mq f_1=p-1/I_0<p`, hence `delta_i f_1<p`.
This identity uses actual Ehrenfest zeroth values, not arbitrary
birth-and-death transforms.

For 4<=r<=m put S_j=a_(m-j), S_0=1, and parameterize the earliest
passage by nonnegative reciprocal slacks:

`P_r(1-h)=f_r/[1-f_r^(-1) sum_(j>=1) sigma_j h^j]`, `sigma_1=d_r`.

For j>=2 the actual slack is
`sigma_j=delta_r f_r^2 (-1)^j P_(r+1)^(j)(1)/(rp j!)`.
Propagate the earlier passages by their first-step fractions. Define
`Q=product_(i=1)^(r-1)f_i`,
`rho=delta_(r-1)f_(r-1)^2/((r-1)p)`,
`Acal=rho Q/f_(r-1)`, and x=d_r. Then

`-partial_(sigma_(r-1))S_r=Q`,
`partial_(sigma_(r-1))S_(r-1)=partial_(sigma_(r-1))S_(r-2)=0`,
`-partial_(sigma_(r-2))S_(r-1)=Acal`,
`partial_(sigma_(r-2))S_(r-2)=0`,
`-partial_(sigma_(r-2))S_r=Bcal`,
`Bcal=Q[1/f_1+sum_(j=2)^(r-1)d_j/f_j+2x/f_r]+x Acal`.       (55l)

Here the first jets may be formal and nonnegative. To prove the formulas,
enumerate first queue failures: a failure at prefix length j has total
count j-1, and every proper-prefix count sum is at least its length.
A jet of order r-2 can therefore occur only at the first position, or
at the second after a first count one. Its first-position coefficient
is `Q[1/f_1+sum_(j=2)^(r-1)d_j/f_j]`. The order r-1 term contributes
`2xQ/f_r`, since
`partial_(sigma_k)P_r(1-h)=h^k P_r(1-h)^2/f_r^2`.
Propagation to the second position contributes x Acal. The highest
order has only the all-zero remaining suffix, giving Q. This proves
(55l) for every r, without extrapolating finite checks.

The actual zeroth-value identity supplies a strict response gain:

`Bcal/Acal >= (r-1)p/[delta_(r-1)f_(r-1)f_1]`
` > (r-1)/f_(r-1) >= r-1+1/p > r`.

Let T=S_(r-3), U=S_(r-2), V=S_(r-1), W=S_r and H=V^2-UW.
Hold the lower slacks actual. Then T,U are fixed actual values with
0<U<=T, by the coupled Ehrenfest survival order, and

`partial_(sigma_(r-1))H=UQ>0`,
`partial_(sigma_(r-2))H=U Bcal-2 Acal V`.

Both partials are strictly positive at every actual point, since V<=U.
They stay positive along a formal lowering path if `U^2>=TV` holds
throughout: then V<=U and the second partial exceeds
`(r-2)Acal U`. The weaker ceiling V<=rU/2 also suffices. An actual
preceding LC inequality alone does not establish this path condition;
lowering the slack increases V.

**The fourth survival gap reduces to two derivative parameters.** Put
F,G,J,K,L,M=P_1(1),...,P_6(1), and (A,B,C,D,E)=(S_1,...,S_5).
For m>=8 set d=(m-5)q and use the sixth-passage floors

`ell_1=L^2[1+dM^2/(6p)]/(5p)`,
`e_L=dL^2 M^3/[(5p)(6p)^2]`,
`h_L=dL^2 M^4/[(5p)(6p)^3]`,
`j_L=dL^2 M^5/[(5p)(6p)^4]`.

Keep x=-L' and s=sigma_2 actual. Lowering sigma_3 to h_L preserves
the three L-jet floors used in (55k), while the first K jet stays
actual. The strong formal third-gap certificate thus supplies C^2>BD
throughout this path. The highest slack sigma_4 does not affect C,D.
Consequently, writing H_5 for D^2-CE as a function of the L slacks,

`(D^2-CE)_actual >= H_5(x_actual,s_actual,h_L,j_L)`, m>=5.       (55m)

The same reduction is equality at m=5,6, where these high-slack floors
are already actual (omit M and set d=0 at m=5). At m=7 a separate
formal endpoint certificate closes the gap left by the actual base
in (55k). Write D_*=D_actual+Acal(sigma_(3,actual)-h_L). With
p=t/(1+t), the reduced rational function C_actual^2-B_actual D_*
has numerator and denominator of degree 44. Its 44 nonzero numerator
coefficients and all 45 denominator coefficients are positive; the
lowest numerator term is 3175200000t and the denominator constant is
23814000000. The sorted nonzero [degree,integer] rows have SHA256
`5932ae6dd97b373783d3ac4ed9b46b4774309ed631a372e6d44e064e42f90d62`
and `af5d0d999cfee91327d0d9ed09b2fdf88d9962a4cfddb0a9cbfbfd3df424bf4d`.
Independent tridiagonal killed-resolvent determinants, Taylor division
and ungrouped queue failures recover both complete coefficient arrays.
Since D_* is the largest D on the lowering path, (55l) now justifies
that path at m=7 too. The selected `--survival-high-slacks` preserves
this endpoint and the depth-four/five symbolic checks without running
the earlier large certificates. authoring AI reproduced it and merged its
compact packet, preserving seventeen earlier certificate entries.
No fourth-gap theorem follows yet. The two remaining actual parameters satisfy

`ell_1<=x<=L(1-L)`,
`e_L<=s<=s_max(x)=x[1-x/(L(1-L))]`.

There is a compact exact test for monotonicity in s. Put k_1=-K',
k_2=K''/2, y=-J', Q_3=FGJ, Q_4=KQ_3, and

`T_*=J(1-B)+yGF`,
`T_2=K(1-C)+k_1 T_*+k_2 Q_3`, `T_3=K T_*+k_1 Q_3`.

With L_j=(-1)^j L^(j)(1)/j!, grouping first queue counts gives
`E=1-L-x(1-D)-L_2 T_2-L_3 T_3-L_4 Q_4`.
At fixed x,sigma_3,sigma_4, write the values at s=0 with subscript zero.
Let `rho_3=(m-3)qJ^2/(3p)`, `rho_4=(m-4)qK^2/(4p)` and

`alpha=rho_3 rho_4 FG`,
`beta=rho_4[T_*+2Q_3(k_1/K+x/L)]+k_1 alpha`,
`R_*=K alpha+rho_4 Q_3`, `u_0=x^2/L`,
`T_20=K(1-C_0)+k_1 T_*+k_20 Q_3`,
`gamma=x beta+T_20+u_0 R_*+2xT_3/L+3x^2 Q_4/L^2`,
`zeta=R_*+Q_4/L`.

All four response coefficients alpha,beta,gamma,zeta are positive.
The exact identities are
`C=C_0-alpha s`, `D=D_0-beta s`, `E=E_0-gamma s-zeta s^2`, hence

`H_5(s)=D_0^2-C_0E_0+(C_0 gamma+alpha E_0-2D_0 beta)s`
` +(beta^2+C_0 zeta-alpha gamma)s^2-alpha zeta s^3`.       (55n)

Thus H_5' is strictly concave. On a nonempty closed interval [s_-,s_+],
it is nonnegative throughout exactly when both endpoint values are
nonnegative. This criterion applies on [e_L,s_max(x)] when that interval
is nonempty. The actual lowering path admits the stronger response
argument below; the separate lower value still has to be signed.
Equivalently `H_5'=alpha E+C(gamma+2zeta s)-2beta D`.
On a region with B>0, C,E>=0 and C^2>=BD, the response bound
`gamma+2zeta s>=2beta C/B` is sufficient. These profile hypotheses
must hold on the formal interval, not merely at the actual endpoint.
There is no unrestricted-half-line version: H_5' has negative leading
term -3alpha zeta s^2, whereas C^2-BD has positive leading term
alpha^2 s^2. The eventual negative C,D are outside the physical domain.

The all-depth prefix argument, response gain, cubic and endpoint criterion
passed independent analytic audit; exact depth-four/five prefix checks
verified the indices and algebra. The audit identified the m=7 scope
gap, which the separate independently reconstructed endpoint closes.
These are reductions and local monotonicity results,
not full survival LC, stationary ULC or NA.

**Response growth and a one-jet reduction.** For m>=8 the remaining
response inequality has an analytic proof. Keep all zeroth transforms
actual, propagate the lower first jets by first step, and allow formal
x,s>=0. Then

`gamma+2zeta s>2beta`.                              (55o)

Put delta_i=(m-i)q, Q=FGJ and
`theta=rho_4/K`, `eta=rho_3 K/J`, `xi=delta_2 GJ/(2p)`,
`w_0=K/(4p)`. The actual return identities used in (55l) imply

`theta F<K/4`, `theta G<K/2`, `0<eta<K<1`, `0<xi<J<1`,
`G<2/3`, `w_0=1-K-theta(1-L)<1-K`.

Normalize the lower first jets by w=k_1/K, v=y/J, u=(-G')/G.
First step gives `w=w_0+theta x`, `v=J/(3p)+eta w`,
`u=G/(2p)+xi v`. Define

`v_0=J/(3p)+eta w_0`, `h=eta(1+xi)`, `tau=2+eta+h`,
`U_0=1/F+G/(2p)+(1+xi)v_0`,
`B_0=U_0+(2+eta)w_0`, `B_1=theta tau+2/L`.

The response formula (55n) simplifies to
`beta/(rho_4 Q)=B_0+B_1 x` and

`(gamma+2zeta s-2beta)/(rho_4 Q)`
` =a_0+a_1 x+a_2 x^2+2[1+eta+1/(theta L)]s`,

where

`theta a_0=1/(FG)+(v_0+w_0-2theta)U_0`
`             +(1+eta)w_0^2-2theta(2+eta)w_0`,
`a_1=(2+eta)U_0+(4+3eta)w_0+h(v_0+w_0)-2theta tau`
`       +2(U_0+w_0-2theta)/(theta L)`,
`a_2=B_1+theta(1+eta)(1+h)+2(h+2+eta)/L+3/(theta L^2)`.

These are direct rational identities: T_*=Q(1/F+u+v),
k_20=Kw^2+rho_4 x^2/L and J_20=Jv^2+rho_3 k_20 suffice
to derive them. The coefficient a_2 is positive. Since U_0>4theta
and h<2eta<2, the first and negative terms of a_1 together exceed
`2theta(2+eta-h)>2theta`; its other terms are nonnegative.

For a_0 put g=G/(2p)=2G w_0/K and A_0=1/F+g. Expanding in v_0 gives

`theta a_0=(1+xi)v_0^2+[A_0+(1+xi)(w_0-2theta)]v_0+P_0`,
`P_0=(1-2theta G)/(FG)+w_0/F+g w_0+(1+eta)w_0^2`
`       -2theta g-2theta(2+eta)w_0`.

The linear coefficient of v_0 is positive because 1/F>4theta
and xi<1. Also `1-2theta G>1-K>w_0`. Dropping two positive
quadratic terms and then using 1/F>4theta/K gives

`P_0>theta w_0[4(1/G+1-G)/K-4-2eta]>(4/3)theta w_0`.

The last strict inequality follows from the exact positive remainder

`4(1/G+1-G)/K-4-2eta-4/3`
` =4(1-K)(1/G+1-G)/K+2(2-3G)(3+2G)/(3G)+2(1-eta)`.

Thus a_0>4w_0/3, proving (55o), in fact with margin
`gamma+2zeta s-2beta>(4/3)rho_4 FGJ w_0`. This proof requires
actual zeroth values, but no derivative floors or restriction on x,s
beyond nonnegativity.

To apply it, fix actual x, set sigma_3=h_L,sigma_4=j_L, and take
`e_L<=s<=sigma_(2,actual)`. Lowering reciprocal slacks decreases
every propagated passage jet. Queue failure polynomials have positive
coefficients, so C,D,E are at least their positive actual values.
A,B stay actual. The strong second-gap identity, at the unchanged
actual K first jet, gives `B^2>A U(k_1)>=A C`, hence C<B.
The third-gap sixth floors remain valid throughout, giving C^2>BD
and therefore D<C. Equation (55o) now yields

`H_5'(s)=alpha E+2beta(C-D)+C(gamma+2zeta s-2beta)>0`.

Combining with (55m) proves
`(D^2-CE)_actual>=H_5(x_actual,e_L,h_L,j_L)` for m>=8.
The same argument permits any stronger componentwise floors lying
between these sixth floors and the actual slacks. It does not give
positive or ordered profiles at arbitrary formal x or s. The sign of
the remaining one-variable lower bound is a separate problem.

**Seventh-passage floors and a fourth-gap scaling theorem.** The full
sixth-floor lower bound need not be positive. At the already critical
scaling p=1-3/m, the actual zeroth transforms converge to

`f_i=(i/3)[E_3-sum_(k=0)^i 3^k/k!]/[E_3-sum_(k=0)^(i-1)3^k/k!]`,
`E_3=exp(3)`.

Putting x=ell_1 and s=e_L,sigma_3=h_L,sigma_4=j_L at these limiting
zeroths gives a strictly negative fourth gap, approximately
-0.00249344543308. All five formal survivals are positive, and the
second-slack envelope margin is positive. This is a failed lower bound,
not an actual survival counterexample or a claim of full moment realizability.

Retain N=P_7(1) as well as M=P_6(1). The geometric seventh-passage
series is a coefficientwise floor. Feed it into the sixth passage:

`Mhat(h)=6p/[1-h+6p+delta_6(1-N/(1-Nh/(7p)))]`,
`mu_j=[h^j]Mhat(h)`, `delta_i=(m-i)q`.

Its zeroth value is the actual M, and
`M^(j+1)/(6p)^j<=mu_j<=(-1)^j P_6^(j)(1)/j!`.
Thus the stronger fifth-passage floors are

`ell_7=L^2(1+delta_5 mu_1)/(5p)`,
`e_7=delta_5 L^2 mu_2/(5p)`,
`h_7=delta_5 L^2 mu_3/(5p)`,
`j_7=delta_5 L^2 mu_4/(5p)`.

They dominate the sixth floors and remain below the actual slacks.
The audited coordinatewise version of (55o) proves, for m>=8,

`(D^2-CE)_actual>=Phi_m(x_actual)`,
`Phi_m(X)=H_5(X,e_7,h_7,j_7)`,
`ell_7<=x_actual<=L(1-L)`.

For the limiting c=3 zeroths, Phi is a degree-six polynomial in X.
Improving only the first floor still gives a negative lower-endpoint
value, approximately -0.000186411190753. Improving all four gives
Phi(ell_7)>0, approximately 0.000613329485744. More strongly, every
coefficient C_j in

`(1+v)^6 Phi((ell_7+v L(1-L))/(1+v))=sum_(j=0)^6 C_j v^j`

is strictly positive. This is an exact certificate: bracket E_3 by
`ell=sum_(k=0)^18 3^k/k!` and
`h=ell+(3^19/19!)/(1-3/20)`. Substitute E_3=(ell+h z)/(1+z)
in each rational C_j and homogenize its numerator and denominator.
All nonzero coefficients of both are positive. Their degrees are
9,11,13,13,13,13,13. The same transformation certifies the negative
sixth-floor examples, the positive moment-envelope margins and the
strict interval width. A separate checker reconstructed all 19 rational
functions using an absorbing queue, then recovered the 38 sign packets
and their 320 strict coefficient signs by Taylor-to-Bernstein conversion.

This yields an actual finite-size result, not just a formal limit:

`a_(m-4)(m,1-3/m)^2 > a_(m-3)(m,1-3/m)a_(m-5)(m,1-3/m)`
`for all sufficiently large integers m`.                         (55p)

Indeed the return-integral density
`[1-3(1-theta)/m]^(-m-2)` converges dominatedly on [0,1] to
`exp(3(1-theta))`. The zeroths, floor coefficients and interval endpoints
therefore converge. Every coefficient of the degree-six Phi_m and each
of its seven Bernstein coefficients depends continuously on this finite
data. Strict limiting signs give positive finite coefficients for every
sufficiently large m. This signs Phi_m on its entire interval containing
the actual derivative; convergence of that derivative is unnecessary.
The comparison above then proves (55p). No explicit threshold or
all-size statement follows from this limit argument alone. The limiting signs, response comparison
and this continuity step have separate independent audits.

The actual fourth gap is also strictly positive at m=5,6,7,8 for every
0<p<1. With p=t/(1+t), its four reduced positive rational expressions
have numerator degrees 29,38,47,56, respectively, with every coefficient
from t through the leading term positive. The primitive denominators
have positive coefficients and constant one; their extracted scalar
prefactors are positive. These are 170 positive numerator coefficients.
The calculation checks normalized first-step series through order four,
agrees with both grouped and ungrouped five-stage queue failures, and
matches all five survival values from independent boundary resolvents
at p=2/5 in each size. Intermediate sizes and rates are not covered
by combining these bases with the eventual scaling theorem.
The selected `--survival-fourth-reduction` reconstructs the response
identities, the limiting floor obstructions and interval certificate,
and these four actual bases. authoring AI reproduced and merged the packet,
preserving all eighteen earlier certificate entries and thirty-one
earlier mathematical helper bodies.

**The fourth survival gap at every size and rate.** For every integer
m>=5 and 0<p<1, with boundary rate one,

`a_(m-4)^2>a_(m-3)a_(m-5)`.                                  (55q)

Use the seventh floors above. At actual zeroths, (55o) bounds the
actual gap below by the sextic Phi(X), evaluated at x_actual in
`[ell_7,U]`, where U=L(1-L). It remains to sign this entire interval.
For m>=10 put m=w+10 and p=t/(1+t), with w>=0,t>0. Define

`b_10=1+11t+w`,
`b_9=(2+10t+w)b_10-10t(w+1)`,
`b_8=(3+9t+w)b_9-9t(w+2)b_10`,
`b_7=(4+8t+w)b_8-8t(w+3)b_9`.

Iteration of the increasing first-step fractions from P_11(1)=0 gives
`N_lo=7t b_8/b_7<=N=P_7(1)`; this is exact at m=10. Dropping the
positive return term at passage seven gives `N<N_hi=7t/(1+8t)`.
Consequently every actual N is represented by finite u>=0 as

`N=R/D`, `R=7t[(1+8t)b_8+u b_7]`, `D=(1+u)(1+8t)b_7`.

The enlarged real domain needs no interpretation as an actual chain.
Its first-step denominators are positive. Moreover
`n_1=N^2/(7p)<=N(1-N)<1-N`, so
`mu_1=M^2(1+delta_6 n_1)/(6p)<M(1-M)<1-M`.
The exact identity

`U-ell_7=delta_5 L^2[(1-M)-mu_1]/(5p)>0`

shows that the first-jet interval is nondegenerate throughout this domain.

The sextic is specified without interpolation. Set g_(5,k) to the
coefficients through h^4 of
`L/[1-(Xh+e_7 h^2+h_7 h^3+j_7 h^4)/L]`.
For i=4,...,1 retain g_(i,0)=f_i and propagate

`g_(i,k)=f_i/(ip)[g_(i,k-1)+delta_i sum_(a=1)^k g_(i+1,a)g_(i,k-a)]`

through k=i-1. If W_l is the set of nonnegative words of total l-1
whose every proper prefix of length b has sum at least b, then

`S_r=1-sum_(l=1)^r sum_(n in W_l) product_(a=0)^(l-1) g_(r-a,n_a)`,
`Phi(X)=S_4^2-S_3 S_5`.

Define the seven unnormalized Bernstein coefficients by

`(1+v)^6 Phi((ell_7+Uv)/(1+v))=sum_(j=0)^6 C_j v^j`.

Their positivity has the following exact integer certificate. To specify
the reduced fractions, use W=m-9=w+1 and continuants

`kappa_7=1`, `kappa_6=4+7t+W-(W+3)N`,
`kappa_i=[10-i+(i+1)t+W]kappa_(i+1)`
`             -(W+9-i)(i+1)t kappa_(i+2)`, i=5,...,1.

Thus f_i=it kappa_(i+1)/kappa_i for i<=6. Exact cancellation gives

`C_j=P_j/[7^6 kappa_1^2 kappa_2^2 kappa_3^2 kappa_4^2 kappa_5 kappa_6^(s_j)]`,
`(s_0,...,s_6)=(0,2,4,4,4,4,4)`.

Both numerator and denominator have N-degree d_j=9+s_j. There are no
additional t or 1+t factors. Put
`K_i=D kappa_i(t,w+1,R/D)` and
`Ptilde_j=D^(d_j) P_j(t,w+1,R/D)`.
The cleared denominator is the same product of K_i, with scalar 7^6
and no remaining D power. Each K_i has positive nonzero coefficients;
their constant terms for i=1,...,6 are
3628800,362880,40320,5040,720,120. All nonzero coefficients of every
Ptilde_j are positive:

| j | N degree | Positive numerator coefficients | Coefficient SHA256 |
|---|---:|---:|---|
| 0 | 9 | 34,310 | `325ebe4c94f70f792b7e846577f1ea79b78e1e63bc43584991230d36330d86de` |
| 1 | 11 | 53,784 | `30ce74d346ea522af141ef1898931f02c0d636d8d24f22d38bd0bea72e829269` |
| 2 | 13 | 79,422 | `c488b252a3e427a08d5a994ebd134e35e843bf200eaf86fbee11460a2f20ef96` |
| 3 | 13 | 79,422 | `93309dad45b19548dc2ac4d17038b55e48c64a343c46c1e1d84193ea99d1e781` |
| 4 | 13 | 79,422 | `944e69b6ac741d4d453e4821d2c6a04d644236a27a1be09cbc54a1dcea630f02` |
| 5 | 13 | 79,422 | `401cabccaa89700c7a9ffdd4673afa35f3fefd4485f3f64bbf7a68310003a02f` |
| 6 | 13 | 79,422 | `385dc012621c1de507a97081031350ee8b7de92772638cb09d6411612ec3b365` |

Hashes encode ascending nonzero [t-degree,w-degree,u-degree,coefficient]
rows as compact JSON. Strictness on w=u=0 follows from

`[t w^0 u^0]Ptilde_j=A_0 binom(6,j)120^(s_j)`,
`A_0=397153497263171567682510237471237734400000000000>0`.

Thus C_j>0 for all t>0,w,u>=0. Finite v covers [ell_7,U), and
C_6=Phi(U)>0 covers the other endpoint. No positivity or ordering
of the free-X survival profiles is assumed here: those properties
were needed only on the actual-x lowering path in (55o).
Applying the interval certificate at x_actual proves (55q) for m>=10.

The remaining actual bases m=5,...,9 have positive rational gap
certificates after p=t/(1+t). Their primitive numerator and denominator
degrees are 29,38,47,56,65. Every numerator coefficient from t through
its leading term is positive; every denominator coefficient is positive,
with constant one. The extracted positive prefactors are respectively
`1/518400`, `1/33592320000`, `1/6858432000000`,
`1/29872281600000000`, and `1/1422757028044800000`.
These 235 positive numerator coefficients cover all omitted sizes.

An independent checker rebuilt the generic fractions by scalar Taylor
division and unshifted-X jets, enumerated queue failures independently,
and used homogeneous endpoint Horner instead of the shifted Taylor
conversion. Native translation and power-sum homogenization recovered
all 485,204 coefficients, all denominator factors, hashes and strict
anchors. Separate scalar endpoint calculations checked C_0 and C_6.
A second checker used killed determinant continuants and queue absorption
recursion for the five actual bases; every coefficient, prefactor and
hash agreed. Original tridiagonal boundary equations at the prescribed
p=2/5 checked all five survivals in each base. Parameter coverage,
seventh-floor lowering and the boundary cases passed analytic audit.
authoring AI reproduced `--survival-fourth-lc` and preserved all 32 earlier
helper bodies and nineteen earlier certificate entries. The compact
`survival_fourth_lc` packet records the construction and audit evidence.

The depth-ten interval is substantive. Its depth-nine predecessor has
a negative lower-endpoint fourth gap at the c=4 Poisson scaling,
although the actual c=4 zeroths give a positive seventh-floor gap.
For the depth-nine endpoint N=511/652 the exact formal gap is
`-6986790400910203906628511329039262245955/`
`9514702004171936119107440549083140704436224`.
This refutes that looser value bound, not the actual inequality.
Together with (55g,i,k), (55q) gives the first four infinite-horizon
survival gaps. Later gaps, finite-time LC and full stationary ULC/NA
remain open.



## Full occupation and principal Green minors

**Complete working proof, all rates and all finite killed kernels.** Write
`x_v=1/lambda_v`, `G=(I-K)^(-1)`, and give the empty determinant value one.
Then

`mu(X=V)=1/D_V(x)`,
`D_V(x)=sum_(A subset V) (product_(v in A) x_v)/det(G_A)`.      (22)

This is an exact full-occupation formula, not a formula for the other
cylinder probabilities. It permits self-loops, nonreversible walks, and
independent positive site rates.

Here is the IDLA constant-term identity needed in the proof. Let `Q` be
any killed kernel on `m` vertices, set
`ell_i(z)=z_i-sum_j Q_ij z_j`, and define, for homogeneous degree-`m`
polynomials `P`,

`Phi_Q(P)=CT [P(z)/(product_i z_i)`
`                   *product_i (1-(Qz)_i/z_i)^(-1)]`.

The geometric expansions converge absolutely on a torus `|z_i|=r_i`
with `Qr<r`; one may take `r=(I-Q)^(-1)1`. The classical
[MacMahon master theorem, Theorem 1.1 and (1.3) in Konvalinka–Pak](https://www.math.ucla.edu/~pak/papers/qmm7.pdf)
gives `Phi_Q(product_i z_i)=1/det(I-Q)`.
Also `Phi_Q(z_i ell_i(z) H(z))=0` for every homogeneous degree-`m-2`
polynomial `H`: cancel the `i`th geometric factor. Every remaining
Laurent monomial then has strictly positive exponent of `z_i`.

Consequently `det(I-Q) Phi_Q(z^a)`, for nonnegative integer vectors
`a` with sum `m`, satisfies precisely the IDLA surplus-routing equations:
if `a_i>=2`, its value is the `Q_ij`-weighted sum of the values at
`a-e_i+e_j`, and a particle sent to the sink contributes zero. Its value
at `a=(1,...,1)` is one. The finite absorbing routing chain gives uniqueness,
so this value is the probability of filling all `m` sites without loss.
This also identifies the full-state coefficient of any homogeneous
polynomial in the IDLA addition operators. In particular, cancellation
of every geometric factor gives

`[full state] (product_i (A_i-J_i)) e_empty = det(I-Q)`,

where sink terms in `J_i` do not contribute: fewer than `m` additions
cannot fill `m` initially empty sites.

Return to the change of basis `M` above. Expanding its last column in the
variables `x_v`, the coefficient of `product_(v in A) x_v` is obtained by
first filling `V\A`, then applying `product_(v in A)(A_v-J_v)`.
IDLA never vacates a filled site. Its remaining walks are therefore walks
on `A` with trace kernel

`Q^A=K_AA+K_A,A^c (I-K_A^c,A^c)^(-1) K_A^c,A`.

The constant-term identity gives coefficient `det(I-Q^A)`. Schur
complementation gives `(I-Q^A)^(-1)=G_A`, proving
`M_(V,V)=D_V(x)`. Since each column `B` of `M` has support on states
with at most `|B|` particles, its last row has no other nonzero entry.
Thus `mu=M^(-1)e_V` has last entry `1/D_V`, proving (22).

**Reader consequence (proved).** Full occupation is strictly increasing in each
sleep rate and nondecreasing under entrywise increase of a killed kernel.
Corollary 2.4 of the [reader](https://raw.githubusercontent.com/jfrichey/jfrichey.github.io/companions-2026-11-14/pagedocs/ai/arw-negative-association.pdf) proves this
from (22) and the positive closed-walk expansion of log principal Green minors.
This is a full-occupation comparison, not configuration stochastic domination.

For disjoint `I,J`, the Green-minor inequality gives the coefficientwise
inequality `D_(I union J)>=D_I D_J`, where all minors use the same ambient
`G`. Equivalently, full-occupation probabilities on the **trace kernels**
satisfy the corresponding product bound. These are not the original law's
marginals, so this does not prove the NA cylinder inequality. For example,
on `sink--0--1` with both sleep rates one, the trace singleton occupation
probabilities are `2/3`, whereas the original marginals are `22/35,23/35`.

`algebra.py` checks (22) against all 44 saved four-site symbolic laws.
A separate directed three-site kernel with self-loops passes the full
identity in three independent inverse rates; at sleep rates `(2,3,5)`,
raw stabilization, both toppling orders, and every addition-stationarity
equation agree. The certificate is in `algebra_certificates.json`.

**Trace laws at rate endpoints (proved).** Fix a nonempty retained set
`A`, keep its sleep rates positive, and send each rate outside `A` to
either zero or infinity. The limiting outside configuration is the
chosen deterministic zero/one pattern, while the law on `A` is the
stationary law of the traced killed walk with its original sleep rates.
This endpoint statement allows general killed kernels.

To see it, skip ineffective sleep instructions at sites with surplus
particles. A zero-rate removed site never retains a sleeper; an
infinite-rate removed site can first be filled from its initial active
particle and then retains one sleeper, routing every extra particle
onward. In both cases subsequent traffic through the removed set follows
the traced walk. The particles initially at zero-rate removed sites can
first be routed to `A` or the sink, giving extra additions to the saturated
sampler on `A`. Fresh instruction tails make these additions independent
of the remaining traced stabilization, and addition-stationarity removes
their effect on its law. The finite absorbing chain with ineffective
sleeps skipped converges to this endpoint chain; the endpoint chain also
absorbs almost surely. The resolvent therefore gives the asserted limit.

**Exact vertex deletion and insertion (complete working proof).** Fix a
vertex `v`, put `A=V\{v}`, and let `nu` be the stationary law on the
traced kernel. Keep arbitrary positive sleep rates and any accessible
killed kernel. Write `a,b` for the unnormalized outside laws with `X_v=0,1`,
respectively, and put `x=1/lambda_v`. All operators below act on column
vectors indexed by stable configurations on `A`.

Let `T_w` be addition at `w` for the traced ARW system, and define

`H=[K(v,sink) I+sum_(w in A) K(v,w) T_w]/[1-K(v,v)]`.

This is the stochastic operator obtained by starting one active at a
zero-rate removed vertex and tracing its first exit. In particular,
`H nu=nu`. Let `T_w^0` be the substochastic addition kernel on `A` which
keeps only stabilizations with **no particle hitting v**; the original
sink still absorbs normally. Define

`J=K(v,sink) I+sum_(w in A) K(v,w) T_w^0`.

A self-loop at `v` is excluded from `J`, since it returns to `v`.
There are two exact identities:

`a=x J b`, `a+H b=nu`.                                  (30)

For the first, stabilize `A` while holding particles at `v`, and repeat
after every toppling at `v`. Let `w` be the expected occupation measure
of the stages with one active at `v` and `A` stable. Sleeping terminates
with law `b=p_v w`. A jump terminates with `v` empty precisely when no
particle returns there, giving `a=(1-p_v)Jw`. This proves the first
identity. For the second, continue an original stabilization by waking
the final sleeper at `v`, if present, and setting its rate to zero.
Every previous sleep there can instead be regarded as holding an active
particle until the next arrival; outside topplings are unchanged.
This is a valid partial stabilization of the zero-rate system. Fresh
unused instructions then give outside law `a+Hb`, and the endpoint
result identifies it with `nu`.



## Complete-graph rank inequalities

**Equal-row vacancy sectors and strict top ULC (complete working proofs).**
The equal-row calculation extends to every vacancy sector, for any a>0.
For hole sets I,J put s=|I|, t=|J|, m=s-t, delta=|J\I|, and
`P_T(z)=product_(k notin T)(1+x_k z)`. If t>s, the basis entry is zero;
otherwise it is exactly

`M_(V\I,V\J)=(-1)^|I\J| (a+s)a^m x_(I\J)/m!`
`   *integral_0^1 z^(a+s-1)(1-z)^delta(-log z)^m P_(I union J)(z) dz`.

To derive this, select a subset A of common occupied sites in the
constant-term formula (23). These cancel their geometric denominators.
Let R=V\(I union A), sigma=det L_R=rho(a+s+|A|), and
k=|I\J|=m+delta. The residual constant term is the probability that
k identical common-jump particles fill the delta prescribed empty targets
and send the other m particles to the original sink, with other exits
forbidden, divided by sigma. With d targets left, the next successful
exit has weight rho/(sigma+rho d) for each target, and
kappa/(sigma+rho d) for the original sink. Choose the target order
and allocate the sink exits among the delta+1 stages. The probability is

`delta! rho^delta kappa^m product_(d=1)^delta(sigma+rho d)^-1`
`       *h_m(1/sigma,...,1/(sigma+rho delta))`,

where h_m is the complete homogeneous polynomial. Differentiating the
beta integral m times gives the displayed logarithmic integral and its
1/m!. Multiplication by det L_(V\I)=rho(a+s) cancels all remaining
rho factors. The uncancelled selected factors give the displayed sign
and x_(I\J); summing A gives P_(I union J). This includes both the
full normalizer and the preceding one-vacancy equations.

For two holes define `b_ij=mu(V\{i,j})/[x_i x_j mu(V)]` and retain
u_i from (35w). Fix I={i,j}, and set
`t_k(z)=x_k(1-z)/(1+x_k z)` for k outside I, ell=-log z.
The exact normalized pair equation L_2 b=f(u) is

`(L_2 b)_ij=integral z^(a+1)P_I`
`   *[b_ij-sum_(k outside I)t_k(b_ik+b_jk)`
`                         +sum_(k<l outside I)t_k t_l b_kl] dz`,
`f_ij=integral z^(a+1)P_I`
`   *[a ell(u_i+u_j-sum_(k outside I)t_k u_k)-a^2 ell^2/2] dz`.

If H_2 is the corresponding same-rank basis block and
X_2=diag(x_i x_j), then `L_2=(a+2)^-1 X_2^-1 H_2 X_2` and
`L_2^-1=(a+2)X_2^-1 N_2 X_2>=0`. Here N_2=H_2^-1 is the
substochastic stabilization kernel retaining all n-2 initial particles.
Nevertheless L_2 is not a killed-network matrix: for four distinct sites
its disjoint-pair entry is
`x_k x_l integral z^(a+1)(1-z)^2P_{i,j,k,l} dz>0`.
Positive diagonal scaling and elimination of extra states from an
M-matrix cannot produce this sign. This obstructs that specific network
extension, not conditional NC.

Indeed, putting `Z_I=a ell+sum_(k outside I)t_k u_k`, the exact
conditional-NC gap vector is `L_2^-1 R`, where

`R_ij=integral z^(a+1)P_I`
` *[u_i u_j-(u_i+u_j)Z_I+Z_I^2/2-sum_(k outside I)t_k^2 u_k^2/2] dz`.

Entrywise R>=0 would suffice, but its sign at unequal rates remains
open. At the single minimum-size discriminator n=4,a=1,x=(1,2,3,4),
all 66 entries in the pair rows and full/one-hole/pair columns agree
with the independent IDLA basis. The disjoint-pair normalized entry
from {2,3} to {0,1} is 2/5>0; all six actual conditional-NC gaps and
all six R entries at that point are positive.

For homogeneous x_i=x, symmetry closes the scalar comparison and proves
strict top ULC at every n>=2,a>0,x>0. Write u_i=u, b_ij=b,
`P(z)=(1+xz)^(n-2)`, and let eta have density proportional to
`z^(a-1)(1-z)^2P(z)`. The one-hole equation gives

`u=integral z^(a-1)(1-z)(1+xz)P v dz`
`       /integral z^(a-1)(1-z)(1+xz)P dz`,

where v(z)=-z log z/(1-z). Put Y=(1-z)^2P and m=n-2.
With t(z)=x(1-z)/(1+xz), the pair diagonal row sum has integrand
`P[1-2m t+m(m-1)t^2/2]=Y''/2`, while
`P(2-m t)=-Y'/(1-z)`. Two integrations by parts in the row sum
and one in the forcing yield

`L=a(a+1)/2 *integral z^(a-1)Y dz`,
`b=[2u E_eta(a v+zv')-a E_eta(v^2)]/(a+1)`.

The boundary terms vanish for every a>0: Y(1)=Y'(1)=0, v(1)=1,
and at zero z^aY, z^(a+1)Y', and z^a vY tend to zero. Rearranging,

`(a+1)(u^2-b)=a E_eta[(v-u)^2]+u[u-2E_eta(zv')]`.       (35aa)

Both terms are positive. To sign the bracket, for every A>0 the
elementary beta integrals give

`E_(Beta(A,2))v=2E_(Beta(A,3))(zv')=A/(A+1)`.

Expand P=sum c_k z^k with c_k>=0. The eta expression mixes these
increasing values at A=a+k with weights c_k B(A,3). Replacing eta
by eta_1 proportional to z^(a-1)(1-z)P tilts the discrete weights by
`B(A,2)/B(A,3)=(A+2)/2`, an increasing multiplier. Hence
`E_eta_1 v>=2E_eta(zv')`. The measure defining u is eta_1 tilted by
1+xz. Strict increase of v gives `u>E_eta_1 v`, so the bracket in
(35aa) is positive. The square term is also strictly positive because
v is nonconstant under eta's positive interval density. Thus u^2>b.

Writing r_j for stationary rank probabilities, this is exactly

`(r_(n-1)/n)^2>r_n r_(n-2)/binom(n,2)`.

Equivalently, two sites have strict negative covariance conditional on
all other sites being occupied. The result covers homogeneous simple
complete graphs after the holding-rate conversion, and all positive
sink weights in the equal-row family. It proves only the top ULC
inequality, not all ranks, unequal-site-rate conditional NC, or monotonicity
of top curvature with feedback. The sector derivation and homogeneous
comparison passed separate independent audits.

`algebra.py --equal-row-vacancies` preserves the 66-entry pair-sector
check and verifies (35aa) against the already saved n=4, original-rate-one
rank law. At holding inverse x=5/4 it gives
`u=797/1239`, `b=18043/50445`,
`u^2-b=2727514/48612165>0`.
Both beta-tilt gains and both positive terms in (35aa) are checked exactly;
no new homogeneous rank law is computed.

There is an exact homogeneous recurrence at every vacancy level. For
1<=s<=n define
`b_s=mu(V\I)/(x^s mu(V))`, |I|=s, with b_0=1,
`Y_s=(1-z)^s(1+xz)^(n-s)`, and let eta_s have density proportional
to z^(a-1)Y_s. Write D=z partial_z and (a)_t for the rising factorial.
Then

`sum_(t=0)^s (-1)^t binom(s,t)a^(s-t)b_t`
`      *E_(eta_s)[product_(j=0)^(t-1)(D+a+j) v^(s-t)]=0`. (35ac)

The coefficient of b_s is (-1)^s(a)_s. The s=0 row instead gives
the full-state normalization `a mu(V) integral z^(a-1)Y_0 dz=1`.
To verify (35ac), group source hole sets J of size t by delta=|J\I|.
Their polynomial sum in the sector formula is

`(1+xz)^(n-s) sum_delta(-1)^delta binom(s,t-delta)binom(n-s,delta)`
`                  *[x(1-z)/(1+xz)]^delta`
` =(-1)^t(1-z)^(t-s)Y_s^(t)/t!`.

Leibniz's rule gives this equality. Combining the sector signs leaves
the integrals `integral z^(a+t-1)v^(s-t)Y_s^(t) dz` with coefficients
a^(s-t)b_t/[(s-t)!t!]. Integrate by parts t times and use

`partial_z^t[z^(a+t-1)f]=z^(a-1)product_(j=0)^(t-1)(D+a+j)f`.

Every endpoint term vanishes even for noninteger a>0. At zero the
k-th derivative of the first factor is
`O(z^(a+s-1-k)(1+|log z|)^(s-t))`; boundary terms have k<=t-1,
so the exponent is at least a+s-t>0. At one, v extends analytically
and the accompanying derivative of Y_s has a zero of order at least
s-t+1. The final integrands have smallest zero-endpoint exponent
a+s-t-1>-1. This proves (35ac) with all endpoint limits justified.

Every differentiated kernel in (35ac) is itself strictly positive.
For m>=1,
`v(z)=z integral_0^infinity dc/[(1+c)(z+c)]`.
Multiplying and using the positive simplex identity
`product_i(z+c_i)^-1=(m-1)! integral_Delta(z+sum theta_i c_i)^-m dtheta`
gives a positive measure rho_m with
`v(z)^m=integral z^m/(z+c)^m rho_m(dc)`.
The simplex identity follows by writing each reciprocal as a Laplace
integral and substituting total time and simplex proportions. The measure
rho_m can have infinite total mass, but its integral against (z+c)^-m
is finite for every z>0.

Fix m, put beta=a+m-1, and write
`partial_z^t[z^(beta+t)/(z+c)^m]=z^beta P_t(z,c)/(z+c)^(m+t)`.
The product rule gives the exact polynomial recursion

`P_0=1`,
`P_(t+1)=z(z+c)partial_z P_t+[a z+(beta+t+1)c]P_t`.

All coefficients are nonnegative, with positive leading coefficient
(a)_t, since a>0 and beta+t+1=a+m+t>0. Thus P_t(z,c)>0 for
z>0,c>=0. On any compact interval [l,u] inside (0,1), each derivative
through order t is bounded by a constant times (l+c)^-m, which is
rho_m-integrable. Differentiation therefore passes through the mixture,
proving `product_(j<t)(D+a+j)v^m>0`. For m=0 the value is simply
(a)_t. This identity and the recurrence were independently audited.

Positivity of these kernels does not sign the alternating sum or prove
ULC of its solution. The proof of b_2^2>b_1 b_3 below retains
the distinct eta_3 average in the s=3 equation.

**The next ULC inequality at every size (computer-assisted proof).**
For every n>=3,a>0,x>0 in the
homogeneous equal-row family,

`b_2^2>b_1 b_3`, or equivalently
`[r_(n-2)/binom(n,2)]^2>[r_(n-1)/n][r_(n-3)/binom(n,3)]`.

This is the next inequality after (35aa), not all-rank ULC. Equivalently,
two sites have strict negative covariance conditional on any specified
outside configuration with exactly one vacant site and all other outside
sites occupied. We prove a stronger auxiliary mixture inequality below.

For n>=3 the
first three measures can be written
`eta_s proportional to z^(a-1)(1-z)^s(1+xz)^(3-s)P(z) dz`,
where `P=(1+xz)^(n-3)`. The following reduction also makes sense for
any nonzero polynomial `P=sum p_k z^k` with p_k>=0. Put u=b_1,b=b_2,
c=b_3. Then

`u=E_(eta_1)v`,
`b=[2u E_(eta_2)(D+a)v-a E_(eta_2)v^2]/(a+1)`,
`c=[a^2 E_(eta_3)v^3-3a u E_(eta_3)(D+a)v^2`
`                 +3b E_(eta_3)(D+a)(D+a+1)v]/[(a+1)(a+2)]`.

There is an exact rational kernel table for this problem. Set A=a+k and
define

`zeta_s(k)=s!/[A(A+1)...(A+s)]`, `n_1(k)=1/(A+1)^2`,
`m_2(k)=(k+2)/(A+2)^2-k/(A+1)^2`, `h_2(k)=2/(A+2)^3`,
`h_3(k)=6/(A+3)^4`,
`j_3(k)=2[(k+3)/(A+3)^3-k/(A+2)^3]`,
`k_3(k)=(k+2)(k+3)/(A+3)^2-2k(k+2)/(A+2)^2+k(k-1)/(A+1)^2`.

These are the integrals against z^(A-1) of, respectively,
`(1-z)^s`, `(1-z)v`, `(1-z)^2(D+a)v`, `(1-z)^2v^2`,
`(1-z)^3v^3`, `(1-z)^3(D+a)v^2`, and
`(1-z)^3(D+a)(D+a+1)v`. The adjoint of D+a under z^(A-1)dz
is -D-k. Moving derivatives onto the polynomial factors cancels the
denominators of v; then use
`integral z^(A+r-1)(-log z)^j dz=j!/(A+r)^(j+1)`.
At zero each boundary term has a positive power of z times a logarithm;
at one the polynomial factors and their needed first derivatives vanish.
Thus these identities hold for every a>0 and k>=0, including real k
for an auxiliary single power.

For any table entry f, set
`T_r[f]=sum_k p_k sum_(ell=0)^r binom(r,ell)x^ell f(k+ell)`.
The unnormalized forms are

`Z_1=T_2[zeta_1]`, `N=T_2[n_1]`,
`Z_2=T_1[zeta_2]`, `S=T_1[m_2]`, `T=T_1[h_2]`,
`Z_3=T_0[zeta_3]`, `M=T_0[h_3]`, `F=T_0[j_3]`, `G=T_0[k_3]`.

With `B=2NS-aZ_1T`, direct substitution gives

`b^2-u c=Phi(P)/[(a+1)^2(a+2)Z_1^2 Z_2^2 Z_3]`,
`Phi(P)=(a+2)Z_3 B^2-a^2(a+1)N Z_1 Z_2^2 M`
`                    +3a(a+1)N^2 Z_2^2 F-3N B Z_2 G`. (35ad)

The denominator is positive, and Phi is homogeneous of degree five in
the coefficients of P. Two auxiliary families have a complete positive
coefficient proof for every a>0,x>=0:

* For P=z^k, k>=0, the reduced denominator of b^2-u c is
  `(A+1)^4(A+2)^4(A+3)^2(a+1)^2(a+2)(Ax+A+3)^2`
  `*[A^2x^2+2A^2x+A^2+Ax^2+6Ax+5A+6]^2`.
  After A=a+k, every one of the numerator's 1,074 coefficients
  in a,k,x is nonnegative. Its constant term in k,x is
  `a^4(a+1)^4(a+2)^3(a+3)^3(a^3+13a^2+45a+47)>0`.
  Thus b^2>u c, including the actual n=3 case k=0.
* More generally, `P=z^i(c_0+c_1 z^d)`, with i,d>=0 and
  c_0,c_1>=0 not both zero, satisfies the same strict inequality.
  To prove this, let f_k denote the vector of the nine unnormalized
  forms for the single component z^k. Define
  `D(A)=A(A+1)^2(A+2)^3(A+3)^4` and `M_k=D(a+k)f_k`.
  Each entry of M_k is polynomial: every denominator in the table,
  including the shifts in T_r, divides D(a+k). If c_0>0, set
  `w=c_1 D(a+i)/[c_0 D(a+i+d)]`. The original moment vector is
  `c_0[D(a+i)]^-1(M_i+w M_(i+d))`. Since Phi is homogeneous
  of degree five, its sign equals that of the polynomial
  `Phi(M_i+w M_(i+d))`. All 296,701 nonzero coefficients of this
  polynomial in a,i,d,x,w are nonnegative. By degree in w their
  counts are `5044,36965,55850,64350,67632,66860`. Its constant
  term in w is strictly positive by the single-component result;
  c_0=0 is the other single-component case. This also handles d=0.
  Thus every two-power mixture is covered, including the actual
  n=4 factor P=1+xz. No integer restriction on i,d is needed for
  this auxiliary integral assertion.

The selected `--third-sector-lc` certificate constructs and signs the
entire polynomials, retaining compact counts and hashes. It also retains
the independently reconstructed 2,890-coefficient specialization
P=c_0+c_2 z^2, checks the quintic identity, and differentiates the
defining integrals independently. For
a=x=1,P=1+z^2 the exact gap is `7146358579/435025920000>0`.
The monomial rational moments were separately audited directly from their
integrands. At that same point the general two-power construction uses
w=4/1125; all nine rescaled forms and its quintic agree exactly with
the independent integrals. These are coefficient proofs, not parameter
samples.

It remains to handle three powers; a linear support reduction will then
give every finite mixture. Order their exponents as i,i+d,i+d+e,
where i,d,e>=0, and substitute `sum_(j=0)^2 q_j M_(k_j)` into Phi,
with the same componentwise clearing M as above. The weight monomials
using all three indices have exactly the following multiplicities.
Every coefficient in a,i,d,e,x of each resulting polynomial is
nonnegative:

| Multiplicities of q_0,q_1,q_2 | Nonzero coefficients |
|---|---:|
| 3,1,1 | 340009 |
| 1,3,1 | 445471 |
| 1,1,3 | 666730 |
| 2,2,1 | 416123 |
| 2,1,2 | 559671 |
| 1,2,2 | 614743 |

These six exact polynomials exhaust the positive compositions of five
into three parts. To construct them, expand B in (35ad), obtaining
seven products of five linear forms. For each product sum all distinct
assignments of the specified index multiset to its five factors. Repeated
factors still retain their multiplicities. This is ordinary coefficient
extraction, with no parameter sampling. An independent implementation
first uses separate variables A_j and then substitutes
`A_0=a+i`, `A_1=a+i+d`, `A_2=a+i+d+e`; the full 3,1,1 polynomial
has the identical 340,009-term hash. An independent audit verifies the
seven-product expansion, all assignment multiplicities, case coverage,
and positive rescaling. At a=x=1 and exponents 0,1,3, all 27 component
moments and all six coefficients agree with directly differentiated
integrals and direct expansion in the three weight variables. For
example the cleared 3,1,1 coefficient is exactly
`4819528811309696318177280`. These finite audits check the construction;
the all-parameter signs come from the complete coefficient scans above.
The canonical `algebra.py --third-sector-lc-all` reproduces all six full
polynomials, their original hashes, these 27 integrals and six coefficient
values, and the one- and two-power certificates. Its selected run passed
all 3,042,747 three-power coefficients, with no negative coefficient.

Together with the one- and two-power cases, these signs prove
`Phi(sum_j q_j M_(k_j))>0` for every nonzero nonnegative combination
of at most three powers. Indeed at least one q_j is positive, and its
pure fifth-power term is strictly positive. Returning to the original
coefficients uses `q_j=p_j/D(a+k_j)>0` on positive coefficients, so
componentwise denominator clearing preserves all these signs.

Now fix any finite P and its u,b, and normalize its coefficients by
Z_3(P). On the original finite exponent set consider nonnegative
coefficient vectors q satisfying

`Z_3(q)=1`, `N(q)-uZ_1(q)=0`,
`2uS(q)-aT(q)-(a+1)bZ_2(q)=0`.

This nonempty set is compact: each coefficient is bounded by the
reciprocal of its strictly positive Z_3 component. Every feasible q
has Z_1,Z_2>0, hence its first two sector values are exactly u,b.
Its third value is the linear functional
`c(q)=[a^2 M(q)-3a u F(q)+3b G(q)]/[(a+1)(a+2)]`.
Choose a maximizer with minimal support. If more than three entries
were positive, a nonzero variation supported there would annihilate
the three equality rows. Its entries have both signs because the Z_3
row is strictly positive. Both sufficiently small directions are feasible,
so maximality makes the objective constant along this variation. Moving
to an endpoint of the feasible segment deletes an entry, a contradiction.
Thus a maximizer q_* has at most three positive entries; dependence
among the constraint rows only improves this support bound.

Since u>0, the three-power result gives
`b(P)^2-u(P)c(P)>=b(q_*)^2-u(q_*)c(q_*)>0`.
This proves the auxiliary inequality for every nonzero finite positive
combination of nonnegative real powers, and in particular for
`P=(1+xz)^(n-3)` at every n>=3. The support-reduction argument, including
strictness, has passed an independent audit. No uniform lower bound on
the gap, degree-independent compactness, or infinite-mixture claim is used.
Lower-rank ULC and full NA remain open.

The support argument itself extends to every vacancy level. For fixed
s use the common-factor measures proportional to
`z^(a-1)(1-z)^r(1+xz)^(s+1-r)P(z)`, 1<=r<=s+1.
Fix the values b_1,...,b_s of an original finite P and normalize
Z_(s+1)=1. The first s unnormalized equations (35ac) are linear
constraints on its coefficients. Every feasible nonnegative vector has
those same b_r: check successively, since the coefficient
`(-1)^r(a)_r Z_r` of b_r is nonzero. The last equation then makes
b_(s+1), and hence `b_s^2-b_(s-1)b_(s+1)`, an affine functional.
Minimizing the gap on the compact feasible set gives a vector with at
most s+1 powers, by the same nullspace argument. Thus it suffices to
prove the next gap on that many powers. Minimizing the gap directly
avoids any assumption about the sign of the auxiliary b_(s-1).
This reduces the higher-sector problem but does not prove its signs.

There is also an analytic mixture comparison at an auxiliary boundary.
Set x=0, fix `P=sum_(k>=1)p_k z^k` with p_k>=0 and P nonzero, and
define b_s by (35ac) for every s using this common P. This is a formal
extension of the recurrence, without an asserted stationary law. For
each fixed s the limit b_s^* as a decreases to zero exists, and

`b_s^*/b_(s-1)^* = [sum_k p_k B(k,s+1) k/(k+1)]`
`                         /[sum_k p_k B(k,s+1)]`.        (35ae)

Indeed put `R_s(z)=integral_z^1 (1-t)^(s-1)dt/t`. Induction, using
`R_(s+1)=R_s-(1-z)^s/s`, gives
`product_(j=0)^(s-2)(D+j)v=(s-1)!z R_s/(1-z)^s`.
The formula includes s=1 with an empty product. Tonelli then gives
`E_(Beta(k,s+1))[product_(j<s-1)(D+j)v]=(s-1)!k/[s(k+1)]`.
Divide (35ac) by a. The coefficient (a)_s/a tends to (s-1)!;
only t=s-1 and t=s survive. Induction proves (35ae). The absence
of a constant term in P gives an integrable bound at zero, uniform
for a in a bounded interval, for all the finitely many differentiated
kernels at each fixed s; at one they are bounded. This justifies the
moment limits and the induction.

The beta weights from s to s+1 are tilted by
`(s+1)/(k+s+1)`, decreasing in k, while k/(k+1) increases. Hence
`(b_s^*)^2>=b_(s-1)^*b_(s+1)^*`, strictly if two distinct powers have
positive coefficients; a monomial gives `b_s^*=[k/(k+1)]^s`.
For each fixed mixed P and fixed s, strictness persists for sufficiently
small positive a by continuity. No threshold uniform in P or s is
claimed. The actual binomial P has a positive constant term and a
singular normalization as a decreases to zero, so this boundary argument
does not apply to it. The kernel identity and limiting argument passed
an independent audit.



## Weak killing at fixed sleep rates

**Positive loss expansion (complete working proof); general NA remains
open.** Let P be a finite irreducible stochastic kernel with no sink,
and set `K_e=(1-e)P`, so every attempted jump is killed with probability
e. Keep every sleep rate fixed and positive. This is distinct from weak
internal transitions and from large sleep rates. Reversibility is assumed
only for the proposed undirected NA application, not for the expansion.

Closed ARW with at most n particles stabilizes almost surely on these n
vertices. A finite legal sequence routes surplus particles to vacant sites,
retaining one at each visited occupied site, and then sleeps every active
singleton. Irreducibility and positive sleep rates give this sequence
positive probability. Abelian least action transfers termination on that
stack cylinder to any fixed legal priority. Finitely many configurations
give a uniform geometric absorption bound for that priority.

Use row kernels in this paragraph. For m particles let Q_m be the closed
transient matrix, R_m its transitions to stable states, and
`G_m=(I-Q_m)^-1`, `H_m=G_m R_m`. Extend H_m by the identity on stable
inputs. Let D_m delete one active at the chosen toppling vertex, weighted
by its next effective jump probability: `1/(1+lambda_v)` for a singleton,
one for a surplus. Its target retains all other active particles and
sleepers. For a final stable B of size m-k, the coefficient of e^k is

`h_(m,0)=H_m`,
`h_(m,k)=G_m D_m h_(m-1,k-1)` for k>=1.

For stable inputs the positive-loss coefficients are zero. The finite
absorption resolvents are analytic at e=0. Exactly k losses are needed,
and at their first possible order every internal factor 1-e can be
replaced by one. This proves the recursion; its factors are nonnegative.
There is no factorial because the losses are chronologically ordered.
For saturated input s_* define
`c_A=h_(n,|A|)(s_*,V\A)`. Then

`Pr(A all vacant)=e^|A| c_A+O(e^(|A|+1))`, `c_A>0`.

Only exact vacancy set A contributes at this order. Strict positivity
follows by prescribing immediate killing at its sites and immediate sleep
at all remaining sites. This proves the valuation for every A.

**A direct finite algebra for all leading coefficients (proved).** There
is no need to expand a symbolic killed-kernel law. Let M_e be the exact
IDLA change of basis, set `D_e(S,S)=e^(n-|S|)`, and define
`Mhat=lim_(e down to 0) D_e^-1 M_e D_e`. Every entry from a k-site
column to a j-site row is O(e^(k-j)): each IDLA term needs at least k-j
losses. Consequently this limit exists. Its diagonal particle-number
blocks are those of M_0, the inverse of closed all-active stabilization.
The latter stabilizes for every input with at most n particles, so these
blocks are invertible. Therefore Mhat is invertible and

`Mhat c=e_V`, where `c(S)=c_(V\S)`.

For explicit operators, fix a proper occupied set S. If i lies outside S,
let Ahat_i fill i immediately. Otherwise its diagonal entry is
`[(I-P_S)^-1 1]_i`, the mean number of jumps before leaving S, and its
entry at S union {j} is `[(I-P_S)^-1 P_(S,S^c)]_ij`.
These are the limits of the sink-exit probability divided by e and of
the closed exit probabilities, respectively. On proper columns put

`Jhat_i=I+sum_j P_ij Ahat_j`,
`Bhat_i=(1+x_i)Ahat_i-x_i Jhat_i`,
`Mhat e_C=(product_(i in C)Bhat_i)e_empty`, `x_i=1/lambda_i`.

One may set full columns to zero: they are never used before the last
of at most n operations from empty. These truncated operators are not
asserted to commute globally; their admissible products inherit the
order independence of the original operators. This gives an exact
rational computation without the weak-killing parameter.

The first vacancy layer has a positive source formula. Write
`L=I-P`, `x_A=product_(v in A)x_v`, and for i outside R define
`t_R(i)=1+P_(i,R)(I-P_R)^-1 1`, with t_empty(i)=1. Then

`u_i=x_i det L_(V\{i})`
`       *sum_(A subset V\{i}) x_A t_(V\({i} union A))(i)`
`                                     /det L_(V\({i} union A))>0`.

This is `-Mhat_(V\{i},V)`, obtained directly from (24). If N_closed
is all-active stabilization on the (n-1)-particle sector, the singleton
vacancy coefficients satisfy `c_single=N_closed u`. This algebraic
source u is distinct from the chronological first-loss departure law.

There is also a closed formula for the expected total number M of jumps
in closed saturated stabilization, counting self-loop jumps. Put
`kappa=sum_i det L_(V\{i})`. Jacobi's identity and
`det(I-(1-e)P)=e kappa+O(e^2)` give, for every nonempty A,
`1/det G_e[A]=e kappa/det L_(A^c)+O(e^2)`.
Differentiating (22), or summing the positive first-vacancy sources, yields

`M=sum_i c_i=sum_i u_i`
` =kappa sum_(nonempty A subset V) x_A/det L_(A^c)`.

The determinant derivative follows from
`adj(L)=kappa 1 pi^T`, where pi is the invariant probability of P.
Thus M is an explicit multiaffine polynomial in the inverse sleep rates,
even though individual c_i need not be polynomials. The departure balance
below gives expected vertex jump counts M pi_v.

Write `w=delta_(s_*) G_n D_n`, `M=mass(w)`, and, for an (n-1)-particle
post-loss state z, put
`u_i(z)=H_(n-1)(z,V\{i})` and
`v_ij(z)=(G_(n-1)D_(n-1)H_(n-2))(z,V\{i,j})`.
The leading pair inequality is precisely

`c_i=w u_i`, `c_ij=w v_ij`, `w v_ij <= (w u_i)(w u_j)`.

Stabilizing z before applying the next loss operator would discard its
transient activity and is not a valid replacement. There is nevertheless
an exact first-loss marginal: if pi is the invariant probability of P,
the departure-site masses are `M pi_v`. Indeed, expected closed jump
counts J_v satisfy `J=P^T J`, since initial and final particle counts
are both one at every site. Irreducibility gives J=M pi. Marking one jump
as the first loss has exactly this occupation measure. The normalized
departure law is independent of sleep rates and priority; the post-loss
configuration given the departure site is not fixed by it.

The natural pointwise bound
`v_ij(z)<=[c_i u_j(z)+c_j u_i(z)]/2` is false. On the closed path
0--1--2 at rate one, with priority 0,1,2, the state z=(0,2,0) has
first-loss weight one. Here `c_0=c_2=33/5`,
`u_0(z)=u_2(z)=2/5`, and `v_02(z)=8/3`, giving gap -2/75.
The actual averaged pair gap is nevertheless `c_0c_2-c_02=414/25>0`.
Thus an averaged argument must use the special first-loss measure w.

Nor can one fix the closed instruction stacks and pair loss marks with
distinct final holes. A two-site edge supplies an exact obstruction at
every positive finite pair of sleep rates. Use priority 0 before 1 and
prescribe initial closed stacks

`0: J,J,S`, `1: J,S,J,S`,

where J jumps to the other site and S sleeps (ineffectively at a surplus).
This cylinder has positive probability. Name the four jump instructions
in unmarked chronological order A=(0,1), B=(1,1), C=(0,2), D=(1,3).
Mark each jump independently for killing with probability e. With no
marks both particles sleep. Killing any one of A,B,C,D leaves the sole
sleeper at 1; killing A and either B or C empties the system. These runs
terminate within the prescribed stacks. In particular a two-loss history
can create a hole at 1 with no corresponding single-mark witness.

Conditional on this cylinder, the occupied-set law is exactly
`Pr({0,1})=(1-e)^4`, `Pr(empty)=2e^2-e^3`, `Pr({0})=0`, with
remaining mass at {1}. The vacancy pair gap is
`-e^2(2-e)(1-e)^4<0`. Its leading coefficients are 4,0,2, so their
pair gap is -2. This excludes conditional pair-NC on each fixed stack
background and a pairing that keeps that background and assigns one
loss to each distinct final hole. It does not refute the unconditioned
weak-killing inequality or rule out a transformation that exchanges
instructions between backgrounds. Effectiveness of sleep instructions,
which changes after deletion, must be included in such a transformation.
The authoring AI checked all branches directly; no simulation was used.

Uniform fixed-rank input is not preserved by closed restabilization,
even for a symmetric kernel, common rate and a vertex-transitive graph.
On the closed four-cycle, choose a uniform two-subset and activate both
sites. Its terminal opposite-pair probability at every lambda>0 is
`1/3+1/[3(2lambda^2+6lambda+5)]`.

Here is a five-state proof. For opposite-pair absorption probabilities
write a,b,z for two actives adjacent, opposite or co-located, and d,e
for an active and sleeper adjacent or opposite. The legal continuous-time
scheduler gives
`(2+2lambda)a=z+b+2lambda d`,
`(2+2lambda)b=2a+2lambda e`, `z=a`,
`(1+lambda)d=(z+e)/2`, `(1+lambda)e=d+lambda`.
With D=2lambda^2+6lambda+5 the unique bounded solution is
`a=z=2(1+lambda)/D`, `b=2(lambda^2+lambda+1)/D`,
`d=(lambda+2)/D`, `e=2(1+lambda)^2/D`.
The uniform input mixture gives `(2a+b)/3=(D+1)/(3D)`.
At lambda1 each adjacent pair has terminal mass25/156 and each opposite
pair7/39. The separate raw absorption check agrees under both priorities.
Symmetry gives uniformity within distance orbits but not between them.

In column convention for the two-particle stabilization kernel N, this
means `(N1)_T=25/26` on adjacent pairs and14/13 on opposite pairs.
Thus a source bound `s<=C1` cannot imply `Ns<=C1` by positivity alone.
Transitivity still makes the actual first-vacancy coefficients c_i equal,
but does not close this proposed maximum principle for the second layer.
The actual second-layer source may compensate the unequal row sums;
no sign for it is established or refuted by this example.

**The product-gap residual is not entrywise positive (exact obstruction,
2026-11-02).** In the closed IDLA block convention write
`c^(1)=N_1u`, `c^(2)=N_2(Uc^(1)-b)`, and
`M_22=N_2^-1`. For `W_ij=c_i c_j`, the candidate
`R=M_22 W-Uc^(1)+b>=0` would imply pair NC because
`W-c^(2)=N_2 R` and `N_2>=0`. This sufficient condition is false even
on the closed three-site path with common `lambda=1/2`.

Substituting this rate into the already-audited two-rate path certificate
gives `c^(1)=(43/2,17,43/2)` and, in pair order `01,02,12`,
`c^(2)=(281/2,226,281/2)`. The pair-gap vector is
`W-c^(2)=(225,945/4,225)>0`. Here `N_2` is one-particle stabilization:
in occupied-site coordinates it is `p(I-qP^T)^-1`, so its inverse in
the stated pair order is

```
M_22 = [[ 3,-1, 0],
        [-2, 3,-2],
        [ 0,-1, 3]].
```

Consequently `R=(1755/4,-765/4,1755/4)`. This includes the actual b
correction, using `M_22 c^(2)=Uc^(1)-b`; omitting b is not the issue.
The source has a negative component although its positive resolvent
image has all components positive. authoring AI independently substituted the
saved rational formulas and checked the matrix product. No new raw
absorption solve or parameter sweep was used. The entrywise-residual
route is closed; the source-averaged pair inequality remains open.

**A restricted small-killing NA theorem (symbolic proof).** For the closed
three-site path with leaf rates r>0 and center rate h>0, the stationary law
is NA for every sufficiently small positive killing probability, with a
threshold depending on r,h. The exact Green/loss equations give, with
`d=2hr+2h+r`,

`c_0=c_2=(8hr^2+10hr+2h+5r^2+6r+2)/(r^2 d)`,
`c_1=2(2hr^2+5hr+h+r^2+2r+1)/(hr d)`.

All six disjoint nonempty cylinder gaps `c_A c_B-c_(A union B)` have
rational expressions with nonzero coefficientwise positive numerators and
denominators in r,h. The finite symbolic certificate checks every matrix
residual and all these coefficient signs, including the singleton versus
pair inequalities. To pass from cylinders to events, take the smallest
minimal witnesses in two increasing vacancy events on disjoint supports.
Their leading event gap is the sum of the corresponding strict cylinder
gaps: the union uniquely determines each witness because the supports
are disjoint. There are finitely many event pairs, so one small-killing
threshold works for all. Threshold decomposition gives arbitrary increasing
functions, and simultaneous complementation gives sleeper NA. This argument
also applies in any fixed size whenever all these leading cylinder gaps
are strict; their general sign remains open.

At rate one on the closed four-cycle, the raw hierarchy in both opposite
priorities gives `c_i=61/6`, `c_ij=205/4` for neighbors and 57 for opposite
vertices, `c_ijk=1079/6`, and `c_V=440`. These agree with the independent
IDLA calculation. The three-site rate-one expansion also agrees with the
full symbolic killed-kernel law. These checks are exact, not simulations.
All 25 disjoint cylinder gaps on this four-cycle are strict, so the same
argument proves NA for sufficiently small killing at this fixed sleep rate.
`weak_killing.py --path3-symbolic`, `--cycle4-raw`, and `--audit` preserve
the symbolic theorem, independent raw checks, rescaled algebra, and the
first-loss/source distinctions in `weak_killing_certificates.json`.

Reversibility of P does not make either the transient ARW chain or the
all-active restabilization kernel reversible. For the latter, on the same
closed path with rates (1,2,3), the two-particle kernel in mask order
3,5,6 has entries
`[[65/116,17/116,11/116],[33/232,105/232,27/232],`
` [69/232,93/232,183/232]]`.
Its opposite three-cycle products differ by -9/26912, precluding every
reversible measure. Both legal priorities give this matrix. Thus a
reversible quadratic-form argument for the loss hierarchy needs more
than reversibility of the underlying walk.

**An aggregate vacancy bound and the weak-killing loss-PGF theorem
(complete working proof, independently audited).** For every finite
irreducible stochastic P and arbitrary fixed positive finite site sleep
rates, put `M=sum_i c_i`, `C2=sum_(i<j)c_ij`. Then

`M^2-C2>=U2>=sum_(i<j) x_i x_j>0` for n>=2,             (36a)

where U2 is the positive two-deletion source defined below. Consequently,
for `K_e=(1-e)P`, one sufficiently small e interval gives

`Psi_e(u)=Cov(N(Y),u^L)>0` for every0<=u<1,              (36b)

simultaneously for every entrance distribution, where X->Y is a stationary
addition and `L=1+N(X)-N(Y)`. Directed edges, self-loops and unequal rates
are allowed. No transitivity or reversibility is required. This proves
the loss-PGF sign in (67) in a new regime, not individual pair NC or NA.

Here is the aggregate cancellation. Let N_e be all-active stabilization,
with column C starting one active on C, and M_e=N_e^-1. The positive
loss expansion gives
`Nhat=lim D_e^-1 N_e D_e=Mhat^-1`, with D_e as above.
Index its blocks by numbers of holes, so the matrices are block lower
triangular. The diagonal block `Nhat_(r,r)=N_r` is closed all-active
stabilization with n-r particles, and is column stochastic. Put

`u=-Mhat_(1,0)>0`, `b=Mhat_(2,0)>0`, `H=Nhat_(2,1)`.

The full-full block equals1. The (1,0) and (2,0) blocks of
Nhat Mhat=I therefore give
`c^(1)=N_1 u`, `c^(2)=H u-N_2 b`. In particular sum u=M. Column C
of H is its leading one-loss law, whose mass is J(C), the expected number
of jumps in closed stabilization from all-active C. Marking the first
killed jump and summing the subsequent closed absorption law proves this
identity, with self-loop jumps included and no factorial. Thus

`C2=sum_(|C|=n-1) u_C J(C)-U2`, `U2=sum_B b_B`.         (36c)

Strict positivity of b follows explicitly from (24). Write L=I-P,
B=V minus {i,j}, R=B minus A, and
`t_R(v)=1+P_(v,R)(I-P_R)^-1 1`. Then

`b_B=x_i x_j det(L_B)`
`       *sum_(A subset B) x_A t_R(i)t_R(j)/det(L_R)>0`.  (36d)

Indeed `a_R^e(v)=e t_R(v)+O(e^2)`. Every displayed principal set is
proper, so irreducibility makes its determinant positive. Empty
determinants and t_empty equal1. The A=empty term is at least x_i x_j.
The earlier one-deletion formula similarly proves u>0.

Initial-particle monotonicity gives `J(C)<=J(V)=M`. For completeness,
couple the full unskipped instruction stacks and use
`0<sleeping<1<2<...`. Every toppling legal in the smaller all-active
configuration C is legal in the larger V. Applying the same instruction
preserves order: a sleep turns a singleton into a sleeper and leaves
surplus unchanged; a jump preserves order at departure and arrival,
including wakeups. After the smaller stabilizes, the larger completes
a legal stabilization. Abelianity identifies its odometer with the usual
one, whose read prefix is longer at every site. Counts of jump instructions
in these prefixes are ordered as well. Ineffective sleeps must be retained
for this coupling. The finite-state closed absorption proof gives finite
expected jump counts for both inputs. Hence

`M^2-C2=sum_C u_C[M-J(C)]+U2>=U2`,

proving (36a). The source u is algebraic and positive; it is not the
chronological post-loss measure w, whose activity distribution remains
complicated. No uniform-sector preservation was used.

To prove (36b), write H_X=n-N(X), H_Y=n-N(Y). The leading vacancy
expansion gives
`Pr(H_X=1)=eM+O(e^2)`, `Pr(H_X=2)=e^2 C2+O(e^3)`, and
`Pr(H_X>=3)=O(e^3)`. Every nonfull-input addition has at most n particles,
so closed absorption gives no-loss probability1-O(e), uniformly over
inputs and pure entrances. A full input cannot have no loss. Since
H_Y=H_X-1 on L=0,

`Pr(L=0)=eM+O(e^2)`,
`E[H_Y;L=0]=e^2 C2+O(e^3)`,
`Psi_e(0)=e^2(M^2-C2)+O(e^3)`.                         (36e)

Stationarity gives the first-order joint masses of (N(Y),L):
(n,0) has mass eM, (n,1) has mass1-2eM, and (n-1,2) has mass eM,
all modulo O(e^2). Every other mass is O(e^2). Thus
`Psi_e(u)=eM u(1-u)+O(e^2)` coefficientwise in u.
There is no assumption of closed absorption with n+1 particles:
for each pure entrance its stationary full-input joint column is
`nu_e(Y)-sum_(X!=V)nu_e(X)T_e(Y,X)`.
All terms on the right are analytic by closed absorption with at most
n particles. This proves analyticity of the full joint array. Its
linearity in the entrance gives uniform coefficient bounds over all
entrances, even those depending on e.

Put d=M^2-C2>0. Since Psi_e(1)=0, divide by1-u polynomially. The
two expansions imply
`Psi_e(u)/(1-u)=eM u+e^2[d+u a_rho(u)]+e^3 b_(e,rho)(u)`,
with the last functions uniformly bounded for small e, u in[0,1] and
all entrances. For a common sufficiently small e this gives

`Psi_e(u)>=(1-u)[(M/2)e u+(d/2)e^2]>0`.                (36f)

This includes u=0 and arbitrarily small1-u. For n=1 the same expansion
works directly with C2=0 and M>0. The threshold depends on P and the
fixed rates. The block orientations, positive source, unskipped-stack
monotonicity and uniform-u argument passed an independent analytic audit.
Conditional stochastic loss order, higher aggregate coefficient ratios,
individual leading pair inequalities and arbitrary-strength gluing remain
open. This theorem does not settle the general weak-killing NA conjecture.

For a fixed increasing rank function h with
`Delta=h(n)-h(n-1)>0`, the same proof gives strict
`Cov(h(N(Y)),u^L)>0` uniformly in u and entrance, after possibly reducing
the e threshold depending on h. The coefficients in (36e) and the first
order in u are multiplied by Delta: lower rank increments require at
least three initial holes on the no-loss event. In particular this
includes the full-occupation indicator. It does not handle all increasing
h uniformly or those with zero top increment.

**Arbitrary weak killing profiles, including boundary killing (proved,
independently audited).** Fix a nonzero profile a_v>=0 and put
`K_e=(I-e diag(a))P`, `k_e=e a`. Take e small enough that e max a<1.
All preceding identities extend with
`J_a(C)=E sum_v a_v #(closed jumps from v)` and `M_a=J_a(V)`.
The leading loss operator now weights a marked departure from v by a_v.
In the nested source formulas replace t_R(v) by

`t_R^a(v)=a_v+P_(v,R)(I-P_R)^-1 a_R`.                  (36g)

The one-deletion source u_a is strictly positive: its A=empty term
uses R=V minus {i}, and a simple path from i to a vertex with a_v>0,
without first returning to i, makes t_R^a(i)>0. The two-deletion source
b_a is nonnegative but may vanish. In particular t_empty^a(v)=a_v,
and the uniform-killing lower bound U2>=sum x_i x_j cannot be retained.

Strictness instead follows from `J_a(C)<J_a(V)` for every proper C.
Common unskipped stacks order jump counts at each vertex. On the
positive-probability cylinder where every site of C immediately sleeps,
the smaller system stops, while the larger has those sleepers and active
particles outside C. Carry one extra active along a simple positive
directed path to v with a_v>0, and prescribe one jump from v, holding
all other actives. These are legal positive-probability instructions;
arrivals at sleepers create a surplus. The larger system has then used
strictly more a-weighted jumps and subsequently stabilizes almost surely.
Finite expectations and the nonnegative pathwise difference prove the
strict inequality. Consequently the same block identity gives

`d(a)=M_a^2-C2(a)`
` =sum_(|C|=n-1) u_a(C)[M_a-J_a(C)]+U2(a)>0`.          (36h)

The entire no-loss/analytic argument is unchanged. Thus (36b) and the
fixed-h extension hold for each fixed nonzero profile, including killing
at a single vertex. Moreover one e interval works for **all normalized
profiles** a>=0,sum a=1, all entrances and all0<=u<1. Indeed M_a is
continuous linear and strictly positive; d(a) is continuous quadratic
and strictly positive even at the simplex boundary. Their minima on
that compact simplex are positive. Closed absorption with at most n
particles has uniform analytic bounds near e=0 because its limiting
matrix is independent of a. Stationary reconstruction retains that
uniformity for the full-input column. The polynomial estimate therefore
holds simultaneously with M_a,d(a) in (36f). No uniformity over arbitrarily
large unnormalized profiles is asserted. For n=1 use M_a=a_1/lambda_1.

There is a direct conductance formulation. For a connected internal
undirected graph with nonzero degrees d_i and conductances c_ij, put
P_ij=c_ij/d_i. Weaken prescribed sink conductances s_i>=0, not all zero, by e,
keeping sleep rates fixed, and set a_i=s_i/d_i. Then
`K_e(i,j)=P_ij/(1+e a_i)`, `k_e(i)=e a_i/(1+e a_i)`.
Every leading vacancy coefficient agrees with the linear-profile model:
at order e^r, a history with r losses uses only their leading factors
e a_i and the closed zeroth-order internal transitions P. Analytic
absorption and stationary reconstruction apply as before. Hence the
same loss-PGF sign holds for sufficiently weak sink conductances,
uniformly over entrances and u. This is a fixed graph/rate/profile
statement, not a conclusion at the original conductances e=1.

**Higher aggregate vacancies and all rank thresholds through three sites
(proved, independently audited).** In this paragraph killing is uniform,
`K_e=(1-e)P`. Put `C_k=sum_(|A|=k)c_A`, `C_0=1`, `C_(n+1)=0`,
and M=C_1. For an all-active input C, let F_k(C) be the coefficient of
e^k in its probability of at least k losses. Common unskipped stacks,
with the order used in (36a), give nested odometers and hence
`L_C(e)<=L_V(e)` pathwise. This orders losses, not terminal sets.
For proper C choose v outside C and a k-set T containing v. Prescribe
the first instruction to kill on T and to sleep elsewhere. This cylinder
has probability `e^k delta_T`, where
`delta_T=product_(i in T)(1+lambda_i)^(-1)`
`        *product_(j notin T)lambda_j/(1+lambda_j)>0`.
The full input loses k particles, whereas C loses fewer. Analyticity of
the at-most-n-particle absorption systems and pathwise nonnegativity on
all other stacks therefore prove

`C_k-F_k(C)>=delta_T>0`, for proper C and1<=k<=n.       (36i)

Let `b_2=Mhat_(2,0)>0`, `b_3=-Mhat_(3,0)>0`, and `U_3=sum b_3`.
The (3,0) block of Nhat Mhat=I, summed over terminal states, gives

`C_3=sum_(|C|=n-1)u_C F_2(C)-sum_(|B|=n-2)(b_2)_B J(B)+U_3`,
`M C_2-C_3=sum_C u_C[C_2-F_2(C)]+R_3`,
`R_3=sum_B(b_2)_B J(B)-U_3`.                          (36j)

There is no factorial: the relevant columns have chronological two-loss
mass F_2(C), one-loss mass J(B), and closed absorption mass1.
The first term in (36j) is strictly positive. Thus R_3>=0 suffices,
but is not necessary. For r=2,3, formula (24), with I=V minus B and
`tau_R(v)=1+P_(v,R)(I-P_R)^(-1)1`, reads

`(b_r)_B=x_I det L_B sum_(A subset B)`
`           x_A product_(i in I)tau_(B minus A)(i)/det L_(B minus A)`.

When n=3, write {i,j,k}=V. Then U_3=x_1 x_2 x_3 and
`(b_2)_{ {i} }/U_3=lambda_i tau_{ {i} }(j)tau_{ {i} }(k)+1-P_ii`
`                     >=lambda_i+1-P_ii`.
The single-particle jump count g_i=J({i}) solves
`[diag(lambda)+I-P]g=1`, including self-loop jumps as rewards. Hence
`(lambda_i+1-P_ii)g_i=1+sum_(j!=i)P_ij g_j>=1`.
Each of the three terms (b_2)_{ {i} }J({i}) is at least U_3, so

`R_3>=2x_1 x_2 x_3`, `M C_2-C_3>2x_1 x_2 x_3>0`.     (36k)

At every size there is also a high-proportional-sleep result. Fix
alpha_i>0, set lambda_i=alpha_i/s and a_i=1/alpha_i, and first form
the e-leading coefficients at fixed s. For proper B,

`J_s(B)=s sum_(v in B)a_v tau_(B minus {v})(v)+O(s^2)`.

At first order exactly one initial singleton makes an exceptional jump,
of weight s a_v. Keep one resident at each site of B minus {v} and
route the surplus carrier until it exits that set; then all particles
sleep. The expected jump count, including the initial jump, is tau.
For justification use the effective chain with ineffective surplus
sleeps removed. At s=0 surplus carriers perform IDLA to vacancies;
irreducibility and at most n particles imply finite closed absorption.
The finite absorption/reward resolvent is therefore analytic near0.
Two exceptional singleton jumps contribute O(s^2).

For a two-hole set {i,j}, the leading source is
`(b_2)_B=s^2 a_i a_j tau_B(i)tau_B(j)+O(s^3)`.
Put `T_3=sum_(|I|=3)a_I product_(i in I)tau_(V minus I)(i)>0`.
Then U_3=s^3 T_3+O(s^4). Group the s^3 coefficient of sum b_2 J
by triples I and a distinguished v in I. With R=V minus I and
{i,j}=I minus {v}, its term is
`a_I tau_(R union {v})(i)tau_(R union {v})(j)tau_R(v)`.
Increasing the resident set delays the walk's first exit, so each term
is at least `a_I product_(w in I)tau_R(w)`. There are three choices of v.
Consequently

`R_3(s)=s^3 r_0+O(s^4)`, `r_0>=2T_3>0`.              (36l)

This proves C_3<M C_2 for sufficiently high proportional rates at each
fixed size and P. It is not an all-rate result or an e-threshold uniform
in s. The all-size all-rate residual in (36j) remains unsigned.

For the threshold `h_k(N)=1_(N>=n-k+1)`, put `d_k=M C_k-C_(k+1)`.
On no loss, H_Y=H_X-1, giving
`Pr(H_Y>=k)=e^k C_k+O(e^(k+1))`,
`Pr(L=0)=eM+O(e^2)`, and
`Pr(H_Y>=k,L=0)=e^(k+1)C_(k+1)+O(e^(k+2))`.
Thus `Psi_(h_k,e)(0)=e^(k+1)d_k+O(e^(k+2))`.
An input with r>=1 holes ending with j>=k holes loses1+j-r particles;
its joint probability is O(e^(j+1)). The leading terminal k-hole mass
therefore comes from full initial input, with L=k+1. As before,
stationarity reconstructs this column without closed n+1-particle
absorption. Uniformly over entrances, coefficientwise in u,
`Psi_(h_k,e)(u)=e^k C_k[u-u^(k+1)]+O(e^(k+1))`.
Dividing by1-u and separating the next-order constant gives, when d_k>0,

`Psi_(h_k,e)(u)>=(1-u)[(C_k/2)e^k u+(d_k/2)e^(k+1)]>0` (36m)

for one sufficiently small e interval, all entrances and0<=u<1.
Indeed the leading quotient is C_k u(1+...+u^(k-1))>=C_k u;
the remaining polynomials are uniformly bounded as in (36f).
For n<=3 all d_k are positive: use (36a), (36k), and C_(n+1)=0.
Taking the minimum of the finitely many thresholds proves the strict
PGF covariance simultaneously for **every nonconstant increasing rank
function h**, all entrances and all0<=u<1. Every such h is a constant
plus a nonzero nonnegative combination of h_k. This does not prove
conditional stochastic loss order. The higher-layer assertion here is
for uniform killing; the profile theorem (36h) only supplies d_1.

**All aggregate gaps for homogeneous complete kernels (proved,
independently audited).** Let P be invariant under every permutation,
not merely transitive, with common fixed positive sleep rate. For n>=2
write P_ii=d and P_ij=c=(1-d)/(n-1) for i!=j, where0<=d<1.
At each fixed0<e<1, erasing holding in K_e=(1-e)P gives a loopless
complete graph with sink ratio `a_e=e/[(1-e)c]` and sleep rate
`ell_e=lambda/[1-(1-e)d]`. The all-sink-ratio extension of (55sand)
therefore applies pointwise. For P=J/n the direct holding representation
has a_e=n e/(1-e) and unchanged lambda; for loopless P the ratio is
e(n-1)/(1-e) with unchanged physical rate. No uniform limit of the
transformed parameters is required.

For a prescribed set A put H_A=1_(all A vacant)=1-OR_A. For disjoint
A,B, (55sand) gives Cov(H_A,H_B)=Cov(OR_A,OR_B)<=0. Let gamma_k
be the common leading exact-hole coefficient for a prescribed k-set.
The containment probability has the same leading coefficient:
`Pr(all A vacant)=e^k gamma_k+O(e^(k+1))`, since all strictly larger
exact hole sets have higher order. Dividing the cylinder inequality by
e^(r+s), for |A|=r and |B|=s, gives

`gamma_(r+s)<=gamma_r gamma_s`, r+s<=n.

In particular, with C_k=binom(n,k)gamma_k and M=n gamma_1,

`C_(k+1)<=[(n-k)/(n(k+1))] M C_k`,
`d_k>= [k(n+1)/(n(k+1))] M C_k>0`, 1<=k<=n.           (36n)

Use C_(n+1)=0 at k=n. Positivity of every gamma_k follows already
from immediate killing on its set and immediate sleeping elsewhere.
The strict aggregate gap comes from the binomial factor; strictness of
individual leading cylinder gaps is not asserted. For n=1 use d_1=M^2.
Thus (36m), with one minimum over all rank thresholds, proves the
weak-killing PGF comparison for **every increasing rank function at
every fixed complete-graph size**, all entrances and0<=u<1. The uniform
quotient bound used in (67q) follows as well, extending its product-recipient
two-count NA theorem to these first blocks of arbitrary size. The e
threshold depends on fixed P,lambda; no site-level NA, conditional loss
order, arbitrary feedback or separate all-rate source-residual sign is
inferred. An independent audit checked the coefficient conversion and
all gluing quantifiers, including zero entrance weights.

**Full complete-graph NA and both BK orientations at weak killing
(proved; independent audits of the input and consequence).** Fix n>=2,
a fully permutation-invariant irreducible stochastic P, and a common
positive finite physical sleep rate lambda. For every sufficiently small
e>0 the saturated stationary law for K_e=(1-e)P is strictly NA, and both
its vacancy and sleeper indicators satisfy BK. One interval works for
all event pairs; it depends on the fixed size, kernel and rate. This is
not an all-sink-strength NA theorem. The n=1 statements are immediate.

First use the holding representation with refresh parameter p and service
alpha. All bounds below are uniform for p in any fixed compact interval
inside(0,1). Delete one tag and put m=n-1. For a prescribed ell-set B
outside, the exchangeable rank average of OR_B is
`f(s)=1-(m-s)_ell/(m)_ell`, with
`Delta f(n-ell)=1/binom(m,ell)`.
In the telescoped deletion source retain the last D term at s=n-ell.
No outside refresh and no service before t gives its pairing at least
`exp(-(m+alpha)t)/binom(m,ell)`; all other pairings are nonnegative.
Integration and (55pow) yield the quantitative trace gap

`nu_trace(OR_B)-E OR_B`
` >=kappa_(n-ell,n-ell-1)/[n(m+alpha)binom(m,ell)]`
` >=d_(n,ell)alpha^ell`, d_(n,ell)>0.                 (36wc)

Let A be any disjoint nonempty r-set, r+ell<=n. Delete one member of A
first and then the rest. Iterated unconditioned trace order only
increases the OR_B expectation. The right vacant-block trace comparison
therefore gives
`E[OR_B | X_A=0]-E OR_B>=d_(n,ell)alpha^ell`.
Writing H_D for the event that D is entirely vacant, multiply by Pr(H_A):

`Pr(H_A)Pr(H_B)-Pr(H_(A union B))`
`                 >=d_(n,ell)alpha^ell Pr(H_A)`.       (36wd)

For a prescribed k-set, the positive weak-killing expansion in this
scale is
`Pr(exact holes D)=alpha^k gamma_k(p)+O(alpha^(k+1))`.
Containment H_D has the same leading term, since larger exact hole sets
have higher order. All gamma_k are positive: immediate killing on D
and sleeping elsewhere has probability
`[q alpha/(n+alpha)]^k p^(n-k)`, so
`gamma_k(p)>=q^k p^(n-k)/n^k`.
Finite closed absorption gives analytic expansions and uniform remainders
on compact p intervals. Dividing (36wd) by alpha^(r+ell) proves

`gamma_r gamma_ell-gamma_(r+ell)`
`                     >=d_(n,ell)gamma_r>0`.           (36we)

This is strict at the leading order, including r+ell=n; finite-alpha
strictness alone would not suffice. The gaps have a positive common
lower bound on a compact p interval for the finitely many r,ell.

For nonconstant increasing vacancy events E,F on disjoint supports, let
r,ell be their smallest witness sizes, and a_r,b_ell their numbers of
witnesses at those sizes. Their marginal leading coefficients are
a_r gamma_r and b_ell gamma_ell. The joint coefficient is
a_r b_ell gamma_(r+ell): disjoint supports make the union determine each
witness uniquely. Thus every event gap has positive leading coefficient
`a_r b_ell(gamma_r gamma_ell-gamma_(r+ell))`.
There are finitely many pairs. One small interval gives all inequalities,
and nonnegative superlevel decompositions give strict NA for disjoint
nonconstant increasing real functions. Simultaneous complementation
transfers vacancy NA to sleepers.

For vacancy BK allow overlapping supports. At the leading product order
the disjoint-occurrence witnesses are distinct unions of disjoint minimum
witness pairs. Choose one ordered decomposition for each distinct union.
Different unions choose different pairs; by (36we) each selected product
coefficient exceeds its union coefficient. Unselected product terms are
positive. If no minimum witnesses can be disjoint, the joint leading
coefficient is zero; if r+ell>n, disjoint occurrence is impossible.
This proves vacancy BK uniformly over event pairs. For sleeper BK use
the proved sparse-complementation lemma: for exact hole states T properly
contained in S, their mass ratio is
`alpha^(|S|-|T|)gamma_|S|/gamma_|T|[1+O(alpha)]->0`,
uniformly on the same compact p interval. Every exact state has positive
mass. These are the required extra hypotheses; unrestricted BK closure
under complementation is not assumed. Shrink the common interval if needed.

Finally write the original kernel as P_ii=d, P_ij=w=(1-d)/(n-1)>0.
Its exact holding-model conversion is
`alpha_e=e/[(1-e)w]`,
`lambda_H(e)=lambda/[n(1-e)w+e]`,
`p_e=lambda/[lambda+n(1-e)w+e]`.
Hence alpha_e/e->1/w and p_e tends to the interior value
lambda/(lambda+n w). The compact-p uniformity just proved applies along
this curve. Equivalently, its e-scale leading coefficients are
`gamma_k^physical=w^(-k)gamma_k(lambda/(lambda+n w))`.
For P=J/n the holding rate stays lambda; for loopless P it is
lambda(n-1)/(n-e). This distinction prevents silently keeping the wrong
rate fixed during the parameter conversion. The theorem proves neither
unequal-rate or general-kernel weak-killing NA nor all-rate complete-graph
NA; the latter remains the centered comparison (55tail) or a weaker sign
for the full centered source.

**Robustness around the complete kernel (proved; independently audited).**
The weak-killing conclusion holds on a relative parameter neighborhood,
chosen before the killing parameter tends to zero. Fix a fully
permutation-invariant irreducible P0 and a common positive rate lambda0.
There are a neighborhood U of (P0,lambda0 1,1) and e0>0 such that

`K_e(i,j)=(1-e a_i)P(i,j)`,
`(P,lambda,a) in U`, `0<e<e0`                         (36wn)

has strict NA and BK in both orientations. Nearby stochastic P may
be nonreversible, and the rates and positive killing profiles may be
unequal. This is a neighborhood of the fixed reference parameters,
not a statement for arbitrary kernels or profiles.

Choose a compact relative neighborhood with P_ij>=w_*>0 for i!=j,
lambda_*<=lambda_i<=lambda^*, and a_*<=a_i<=a^*. A uniform active-site
scheduler has a uniformly positive bounded-length stabilization prefix
from every state with at most n particles: move each surplus particle
directly to an empty site, then sleep each remaining active singleton.
At most2n steps suffice. Thus the finite closed transient absorption
matrix has a uniformly bounded inverse, including at e=0. Its analytic
extension is relative to the affine row-sum constraints; a zero diagonal
at a loopless reference causes no issue. No n+1-particle closed
absorption is used.

Grade this transient matrix by total particle count. Every downward
transition has a factor e, so block inversion gives the exact analytic
factorization
`Pr(exact holes I)=e^|I|[c_I(P,lambda,a)+e R_I(e,P,lambda,a)]`.
Remainders are uniform on a smaller compact neighborhood. Immediate
killing on I and sleeping elsewhere bounds c_I below by
`product_(i in I)[a_i/(1+lambda_i)]`
` *product_(j notin I)[lambda_j/(1+lambda_j)]>0`.
At the reference tuple c_I is the physical gamma_|I|. The finitely
many strict leading gaps in (36we), by continuity, stay uniformly
positive on a neighborhood chosen independently of e. Minimum-witness
and distinct-union expansions therefore give every NA and vacancy-BK
inequality on one e interval. Exact-state ratios for T properly
contained in S are bounded by C e^(|S|-|T|), uniformly in parameters.
The sparse-complementation proof is thus uniform too and supplies
sleeper BK. This retains the strict leading-order input; continuity
at a single positive e would not give these quantifiers.

One bounded attempt at the remaining source sign gives the legal carrier
prefix bound `J(B)>=tau_(B minus {v})(v)/[lambda_v+ell_(B minus {v})(v)]`
for each v in B, where
`ell_R(v)=1-P_vv-P_(v,R)(I-P_R)^(-1)P_(R,v)`.
Hold the residents of B minus {v} and follow a carrier until it sleeps
at v or exits B. Its departure excursions have mean tau and return
probability1-ell, proving the bound by one-step renewal. These prefixes
reuse the same residents and their costs cannot simply be added.
The source ratio `(b_3)_(B minus {v})/(b_2)_B` is likewise not controlled
by tau-monotonicity alone: grouping the nested source terms introduces
ratios ell_(B minus {v})(v)/ell_R(v). This route is parked pending a
source-weighted jump-cost estimate; no actual counterexample is claimed.

## Weak internal transitions at fixed sleep rates

Here the sleep rates stay fixed and positive, while `K_epsilon=epsilon P`
is the internal jump kernel; the remaining row mass goes to the sink.
Assume `P_ii=0`, `P>=0`, and choose `epsilon` small enough that its rows
are substochastic. Write `p_i=lambda_i/(1+lambda_i)`, `q_i=1-p_i`, and
`pi_p` for the independent Bernoulli law with these parameters. This is
a separate parameter regime from the small/large sleep-rate theorems.

**Analytic support lemma (proved).** Give each allowed directed arc an
independent variable `z_ij`. The stationary law is analytic near zero:
the basis matrix is analytic there and its zero-kernel value is invertible.
If the directed support of a Taylor monomial has no strongly connected
component meeting both supports of `f,g`, its coefficient in
`Cov(f,g)` is zero. Indeed, set every other arc to zero. Component
factorization makes `f` and `g` independent throughout this restricted
parameter family. Likewise, a coefficient of the law itself vanishes
if its arc support contains an arc between distinct strong components:
the component laws do not depend on that outgoing arc's destination.

The two-site formula above therefore determines the entire second-order
law. With `pi^(ij)` obtained from `pi_p` by replacing `p_i,p_j` by their
squares, it gives

`mu_epsilon=pi_p+epsilon^2 sum_(i<j) P_ij P_ji (pi_p-pi^(ij))`
`             +O(epsilon^3)`.                              (37)

For disjoint increasing functions `f,g`, define
`I_i(f)=E_pi_p[f(X with X_i=1)-f(X with X_i=0)]`. Expanding covariance
in (37) gives

`Cov(f,g)=-epsilon^2 sum_(i in supp f,j in supp g)`
`          P_ij P_ji p_i q_i p_j q_j I_i(f) I_j(g)`
`          +O(epsilon^3)`.

Consequently, if `P_ij P_ji>0` for every distinct pair, the stationary
law is **NA for all sufficiently small positive epsilon**, at any fixed
positive site rates. Every nonconstant increasing event has a positive
influence at an essential coordinate; the leading coefficient is then
strictly negative for every disjoint nonconstant event pair. There are
only finitely many such pairs. This proves the assertion, with a
threshold depending on `P` and all rates. On a sparse graph the same
argument handles event pairs whose essential supports have an edge
between them. The deletion argument below handles separated supports.

**Single-destination deletion lemma (proved).** Let `v` have self-loop
probability `r<1`, exactly one possible other internal destination `w`
with probability `z`, and sink probability `1-r-z`. Fix the other rows
and rates, with an accessible sink at `z=0`. Let `nu(z)` be the
stationary law on the trace deleting `v`, and `mu_out(z)` the original
outside marginal. Then `mu_out(z)-nu(z)` is divisible by `z^2` as an
analytic vector near zero. This is a marginal statement, without
conditioning on the deleted coordinate.

To prove it, use the excursion notation (30)–(32), let `x=1/lambda_v`,
and denote traced addition at `w` by `T`. The first forced jump gives
`E=r I+z F`, with `F` analytic: a self-loop immediately returns to the
same one-active-`v` stage, while every other return begins with the
jump to `w`. Also
`H=I+z(T-I)/(1-r)`, `T nu=nu`. Set `c=1+x(1-r)>0`. Then

`b=(c I-xz F)^(-1)nu`
` =c^(-1)nu+(xz/c)(c I-xz F)^(-1)F nu`.

Since `mu_out-nu=(I-H)b`, this yields the exact factorization

`mu_out-nu= -x z^2/[c(1-r)]`
`             *(T-I)(c I-xz F)^(-1)F nu`.                  (38)

All operators and `nu` may depend on `z`; only analyticity and
`T nu=nu` were used. Thus the identity is also valid modulo the square
of an individual formal outgoing arc variable. When contraction creates
parallel arcs their sum need not have square zero; that distinction
prevents an unrestricted repeated-deletion argument.

**Cut-vertex deletion modulo edge squares (proved).** Suppose deleting
`v` splits the outside graph into blocks `B_1,...,B_k`, with no transitions
between blocks. Block `B_i` communicates with `v` only through one vertex
`w_i`: write `a_i=K(w_i,v)`, `b_i=K(v,w_i)`, and `r=K(v,v)<1`.
Fix the other internal transitions and adjust sink probabilities as these
variables change. Assume each killed block at zero ports is absorbing.
In the ring with `a_i^2=b_i^2=0` for every `i`,

`mu_out=nu`,

where `nu` is the stationary law on the trace deleting `v`.

For a complete proof, stabilize the outside during the first forced
departure from `v`, holding all particles returning there. After a
departure into `B_i`, two or more returns have a factor `a_i^2` and
vanish in this quotient. With exactly one return no further toppling
at `v` is needed. Therefore
`E=r I+sum_i a_i b_i D_i`, where `D_i` is the coefficient of a single
return from an addition at `w_i`, evaluated at `a_i=0`. It acts only
on `B_i` and is independent of all port variables. The `D_i` commute.
This expansion also accounts for later surplus routing: it occurs only
after at least two returns from the initially visited block, already
carrying the vanishing factor `a_i^2`.

Put `c=1+x(1-r)` and, for `S subset {1,...,k}`, let
`h_S=product_(i in S) a_i b_i`, `D_S=product_(i in S) D_i`.
The finite resolvent expansion is

`b=sum_S |S|! x^|S| c^(-|S|-1) h_S D_S nu`.

Also `H-I=(1-r)^(-1)sum_j b_j(T_(w_j)-I)`, where the additions are on
the traced system. Its application to the `S=empty` term vanishes by
stationarity. In a nonempty term, indices `j in S` give `b_j^2=0`.
For `j notin S`, multiplication by `h_S` permits setting `a_i=b_i=0`
for all `i in S` inside the remaining factors. Those blocks then
disconnect from the trace, whose law factors as
`(product_(i in S) nu_i) tensor nu_remaining`. The operator `D_S`
changes only the disconnected blocks. Traced addition at `w_j` acts
only on the remaining system and leaves `nu_remaining` invariant.
Thus every term of `(H-I)b` vanishes. Since
`mu_out-nu=(I-H)b`, the lemma follows.

The argument is an identity of analytic coefficients in the port
variables, with the other kernel entries as parameters. It therefore
remains valid after substituting nilpotent products for ports and
analytic nilpotent expressions for other entries, whenever the
constant denominators remain invertible. In particular it can be
applied successively along a bidirected path, including the self-loops
created by tracing. This substitution point matters because successive
trace entries need not be independent variables.

**One directed cycle, arbitrary length (proved).** Suppose the only
internal arcs form one oriented simple cycle `C`. Work in the ring
`Q[z_e:e in C]/(z_e^2:e in C)` and put `c_C=product_(e in C) z_e`.
For every retained nonempty set `A`, successively delete vertices outside
`A`. Each has one other internal destination, so (38) preserves the
retained marginal in this ring. Contracted arc probabilities are products
along disjoint segments, whose squares still vanish. The remaining law
is the traced law. For every nonempty `B`, principal minors give
`1/det(G_B)=1-c_C`; hence the full-occupation identity on that trace gives

`E[product_(i in A) X_i]=p_A+p_A(1-p_A)c_C`,
`p_A=product_(i in A) p_i`.

Moments determine the law, so in the quotient

`mu=pi_p+c_C(pi_p-pi_(p^2))`.                              (39)

For `K_epsilon=epsilon P` on this cycle of length `m`, all lower-order
terms vanish by the support lemma. Thus for disjoint increasing `f,g`,

`Cov(f,g)= -epsilon^m (product_(e in C) P_e)`
`              *(E_pi_p f-E_pi_(p^2) f)`
`              *(E_pi_p g-E_pi_(p^2) g)+O(epsilon^(m+1))`.  (40)

Both expectation differences are strictly positive for nonconstant
increasing functions. The same finite-event argument proves NA for
all sufficiently small positive epsilon on any such fixed cycle, with
arbitrary positive site rates.

**Weak-interaction NA theorem (complete working proof).** On every fixed
finite bidirected support and at every fixed collection of positive
site rates, the stationary law of `K_epsilon=epsilon P` is NA for all
sufficiently small positive `epsilon`. Reversibility of the weights
is unnecessary; the support must satisfy `P_ij>0 iff P_ji>0`.

Let `d` be the undirected distance between the essential supports of
nonconstant increasing functions `f,g` with disjoint supports. If no
component meets both, component factorization gives independence.
Otherwise
every nonzero covariance monomial has degree at least `2d`.
Indeed, a strong component meeting both supports must cross each of the
`d` distance-level cuts in both directions. An arc crosses at most one
such cut. At degree exactly `2d`, there is one upward and one downward
arc at each level, so the support is the union of two oppositely directed
geodesics. Their intersections occur in distance order: they form a chain
of directed cycles sharing at most cut vertices. In particular the
monomial is squarefree, and this minimal support meets each function's
essential support at only one endpoint.

Each such monomial contributes `-p_i q_i p_j q_j` to the endpoint
covariance. To prove this, restrict the kernel to its chain of cycles
and work modulo the squares of these directed-edge variables. Delete
non-cut vertices other than the two endpoints using (38); each has one
other internal destination. Each directed segment contracts to a product
of its arc variables. The result is a bidirected path, possibly just
one edge. Delete its internal vertices using the cut-vertex lemma.
Each new arc retains a factor equal to the product along its oriented
segment, so its probability has square zero. The lemma allows the
self-loops and shared nilpotent parameters produced by previous traces.
Thus the endpoints' joint law equals their two-site traced law in the
original quotient ring.

Let `Q` be that final two-site kernel. The product `Q_ij Q_ji` is exactly
the product of all directed arc variables in the original chain: each
direction must traverse its successive segments, and any additional
closed excursion repeats at least one of these arcs. Both diagonal
entries of `Q` have zero constant term. Remove self-loops by replacing
the two sleep rates by `lambda_i/(1-Q_ii)` and
`lambda_j/(1-Q_jj)`, and the opposite jump probabilities by
`Q_ij/(1-Q_ii)`, `Q_ji/(1-Q_jj)`. The exact two-site covariance formula
then gives coefficient `-p_i q_i p_j q_j`; all nonconstant corrections
to its remaining factors multiply the full arc product and vanish.

The rest of the observed coordinates lie outside this minimal support
and have their independent Bernoulli laws. Their averaging supplies
the influences `I_i(f),I_j(g)`. Summing over the two geodesics gives

`Cov(f,g)=-epsilon^(2d) sum_(i,j:dist(i,j)=d)`
` (P^d)_ij (P^d)_ji p_i q_i p_j q_j I_i(f) I_j(g)`
` +O(epsilon^(2d+1))`,                                    (41)

The coefficient is strictly negative: a nearest pair of essential
coordinates has both positive walk weights and positive influences.
Take the minimum validity threshold over the finitely many disjoint
increasing-event pairs, and then use threshold decomposition for general
increasing functions. This proves NA, with no volume-uniform threshold
and no continuation to `epsilon=1` claimed.

**BK in both orientations near a product law (proved corollary).** Let
`pi` be any full-support Bernoulli product law on a fixed Boolean cube,
and set `m=min_x pi(x)>0`. Every NA probability law `rho` satisfying
`||rho-pi||_TV<m^2/3` is BK, and its coordinatewise complement is BK.
Here total variation is the supremum of event-probability differences.

The input is strengthened product BK: if `D` consists of pairs with
disjoint positive witnesses for increasing events `A,B`, then
`pi(A square B)<=(pi tensor pi)(D)`. See
[Cerf–Tesio, Theorem 1.5 and Section 2.4](https://arxiv.org/html/2607.02272v1).
Their coordinate-splitting proof applies also to unequal Bernoulli
parameters: conditional on the other coordinates, the indicator of `D`
is an increasing two-bit function `phi` that cannot be AND, because
disjoint witnesses use the coordinate on at most one side. Splitting
coordinate `i` changes its expectation by
`p_i(1-p_i)[phi(1,0)+phi(0,1)-phi(0,0)-phi(1,1)]>=0`.
If an essential coordinate belongs to both events, choose minimal
occupied witnesses `I,J` containing it. The pair `(I,J)` is outside `D`,
so the product BK gap is at least `pi(I)pi(J)>=m^2`.

For the corollary, changing from `pi` to `rho` changes a BK gap by at
most `3||rho-pi||_TV`. Thus every event pair with overlapping essential
supports still has a positive gap. For disjoint essential supports,
disjoint occurrence equals intersection and NA gives the required gap.
This covers all pairs, including the immediate constant cases. The
complemented law is NA as well: negate both decreasing functions in
the covariance inequality. Its distance from the complemented product
law is unchanged, as is that product's smallest atom `m`. Apply the same
argument again.

Since the weak-interaction laws above converge to `pi_p`, both sleeper
and vacancy indicators are therefore **BK for all sufficiently small
positive epsilon**, on every fixed bidirected support and on a single
oriented cycle. This is a neighborhood result; it does not assert
unrestricted BK closure under complementation or that NA always implies BK.

**A strict BK core with independent sites (proved; independently audited).**
Let mu be a full-support measure on a finite core C, strictly BK for
every pair of nonconstant increasing core events. Let pi be a
full-support Bernoulli product law on a disjoint finite set W. Every
NA law sufficiently close in total variation to mu tensor pi is BK.
If the complemented core has the same strict BK property, both
orientations hold on one neighborhood. Its radius may depend on mu,pi.

First, adjoining one independent Bernoulli(p) coordinate preserves
ordinary BK for any BK measure. For increasing E,F, write their core
slices as E_0 subset E_1 and F_0 subset F_1. Put
`D_st=E_s square F_t`, `d_st=mu(D_st)`,
`g_st=mu(E_s)mu(F_t)-d_st`,
`U=D_10 union D_01`, `u=mu(U)`, and
`v=mu(D_10 intersect D_01)`.
The disjoint-occurrence slices are exactly D_00 and U: at most one
positive witness uses the new coordinate. Monotonicity gives
`D_00 subset D_10 intersect D_01` and `U subset D_11`. The exact gap is

`G_(mu tensor Ber(p))(E,F)`
` =(1-p)^2 g_00+p^2 g_11+p(1-p)(g_10+g_01)`
`   +p^2(d_11-u)+p(1-p)(v-d_00)>=0`.                  (41bk)

Indeed d_10+d_01=u+v, so the occurrence mass subtracted on the right
is (1-p)d_00+p u. Iterating over W yields
`G_(mu tensor pi)(E,F)>=sum_(x,y)pi(x)pi(y)G_mu(E_x,F_y)`,
where x,y are separately fixed leaf states. If both events depend
essentially on the core, some pair of slices is nonconstant and this
sum is strictly positive.

If at most one depends on the core, say F is leaf-only, then for each
core state c, `(E square F)_c=E_c square F`. A leaf witness for E_c
lifts by adjoining all occupied core sites of c; monotonicity makes
it a global witness disjoint from the leaf witness for F. Thus
`G_(mu tensor pi)(E,F)=sum_c mu(c)G_pi(E_c,F)`.
When E,F share an essential leaf coordinate, some E_c shares it with
F, and full-support product BK is strict for that pair. Hence every
pair with overlapping essential supports has a positive gap. The zero
nonconstant pairs are precisely those with disjoint essential supports
and at most one event depending on the core; these pairs are independent.

Take delta>0 to be the minimum gap over overlapping essential-support
pairs. There are finitely many. Perturbation by total variation less
than delta/3 keeps those gaps positive. For disjoint essential supports,
disjoint occurrence equals intersection, so NA of the perturbed law
supplies BK. Apply the same argument to the complemented reference law
and take the smaller radius when both core orientations are strict BK.
This proves the assertion without conditional NA, arbitrary dependent
BK-product closure, or a witness-preserving coupling. The weak-killing
complete core in (36wc)–(36we) satisfies both strict core hypotheses.

**Exact evidence and a failed stronger claim.** The proposed equality
of all occupation moments with trace full-occupation probabilities modulo
edge squares holds in the saved three-/four-site paths, the four-site
star, and two oriented triangles sharing one vertex (five sites), at the
listed unequal rates. It fails on the complete three-site graph and
the four-site diamond: respectively 30 of 448 and 114 of 3,840 compared
coefficients disagree. On complete3 with rates `(2,3,5)`, the singleton
`{0}` difference has coefficient `-1/12` at
`z_01 z_10 z_12 z_20`. These failures concern a stronger formal moment
identity, not the trace order comparisons (29) or (35).

All computed leading pair coefficients still agree with (41). For
example, on the path with `P_01=P_21=1`, `P_10=P_12=1/2`, and rates
`(2,3,5)`, the endpoint covariance is
`-5 epsilon^4/648+O(epsilon^5)`. The full symbolic-kernel law at
`epsilon=1/2` agrees with raw absorption under both toppling orders and
all addition-stationarity equations. `weak_interaction.py` also checks
the nilpotent IDLA relations, commutation, and law normalization; its
saved certificate lists the arc order, rates, comparisons, and exact
leading coefficients. A separate 128-coefficient check deletes a vertex
between a singleton branch and a two-site branch while keeping fixed
internal jumps and self-loops; only the four port arcs are nilpotent.
It verifies the cut-vertex marginal identity with unscaled outside
dynamics. These checks support the working proofs above, rather than
replacing their all-size arguments.



## A weak bridge between interacting components

**Exact response identity (proved).** Let a finite killed kernel `K`
have an irreducible internal graph, stationary law `nu`, and a distinguished
vertex `u` with positive direct killing. Fix all positive sleep rates.
Let `T` be addition at `u`, and let `D` be its addition kernel with one
absorbing departure from `u` marked, with coefficient one in place of
the killing probability. Equivalently, weight each addition trajectory
by its number of such absorbing departures divided by `K(u,sink)`.
The resulting `D` is nonnegative but need not be stochastic.

Put `g=G_uu`, `G=(I-K)^(-1)`. Stationary expected jump counts solve
`(I-K^T)h=e_u`, by particle balance at each vertex, so
`1^T D nu=g`. Let `U` be the derivative of the stationary law when
`K_uu` is increased by `z` and `K(u,sink)` is decreased by `z`.
Then

`U=g nu-D nu`.                                             (42)

Indeed the derivative of the addition kernel is `(T-I)D`: keep the
single marked departure until the original addition has stabilized,
then either return it as one fresh addition at `u` or leave it absorbed.
Abelianity permits this order, and the unused stacks are independent
even after weighting the exposed trajectory by its marked exits.
Differentiating `T(z)nu(z)=nu(z)` gives
`(T-I)(U+Dnu)=0`. Addition at `u` has a unique stationary law because
`u` reaches every vertex. Thus `U+Dnu` is a scalar multiple of `nu`;
its mass is `g`, proving (42). All differentiations concern finite
absorbing chains and hence analytic functions near the base kernel.

Now take two such independent components `A,B`, with ports `u,v`,
laws `nu_A,nu_B`, Green diagonals `g_A,g_B`, and responses `U_A,U_B`.
Replace an amount `a` of killing at `u` by a jump to `v`, and an amount
`b` of killing at `v` by a jump to `u`. In the quotient `a^2=b^2=0`,
the joined stationary law is

`mu=nu_A tensor nu_B + ab Delta`,
`Delta=g_B U_A tensor nu_B + g_A nu_A tensor U_B`
`                         -U_A tensor U_B`.                (43)

In particular, for functions `f` on `A` and `h` on `B`,

`Cov_mu(f,h)=-ab U_A(f)U_B(h)+O(ab(a+b))`.                  (44)

For a proof of (43), either one-way coupling leaves the product law
unchanged by component factorization. Thus only the mixed coefficient
`Delta` remains. The addition-at-`u` kernel through this order is

`T_A tensor I + a D_A tensor (T_B-I)`
`              + ab (T_A-I)D_A tensor D_B`.

Route the first crossing into `B` only after the original addition in
`A` has stabilized; then route its single return after stabilizing `B`.
Two crossings in either direction vanish in the quotient. The same
unused-stack argument justifies each composition. Stationarity implies

`((T_A-I) tensor I)(Delta+D_A nu_A tensor D_B nu_B)=0`.

The corresponding equation for addition at `v` puts this vector in the
stationary space of both independent operators, hence it equals
`g_A g_B nu_A tensor nu_B` by its mass. Substituting (42) proves (43).
Subtracting the two marginal responses gives (44). Analyticity and the
vanishing of cross covariance whenever `a=0` or `b=0` justify its stated
remainder for small nonnegative `a,b`.

**Open order target:** does `U(f)>=0` hold for every increasing `f`,
equivalently `Dnu <=st g nu`? This is stochastic monotonicity under
replacing killing by a self-loop. It differs from the proved monotonicity
in a site's sleep rate. A positive answer would sign the leading
correlations between arbitrary interacting components; zero responses
would still require higher terms, and (44) alone is not full NA.
Neither this order assertion nor an all-rate gluing theorem is claimed.

**Conditional weak-gluing theorem (proved).** Suppose each component law
has full support and is strictly NA on disjoint nonconstant increasing
events. Choose a nonempty set of vertices with positive direct killing
in each component, and assume each chosen vertex has `U(f)>0` for every
nonconstant increasing `f`. Join chosen ports by the bridge above.
For all sufficiently small positive `a,b`, the joined law is strictly
NA and retains strictly positive self-loop responses at every chosen
vertex of either component.

For NA, at the product law a disjoint event pair has a strictly negative
covariance whenever both events depend essentially on the same component.
Condition on the other component and use its NA for the covariance of
the conditional means. Strictness in the first component occurs for some
conditioning configuration: the two events use disjoint coordinate sets
in the other component, so their nonconstant sections can be selected
simultaneously, and full support assigns positive probability to that
configuration. These strict gaps persist by continuity. The remaining
nonconstant pairs depend on separate components, and (44) gives a strictly
negative leading coefficient. Finitely many pairs give one threshold.

For a chosen vertex `w in A`, an increasing event that depends essentially
on `A` has a strictly positive response already at the product law,
by averaging its sections. If the event depends only on `B`, differentiate
the `B` marginal of (43) with respect to the self-loop at `w`. Its response is

`ab G_A(u,w)G_A(w,u) U_B(f)+O(ab(a+b))`,

strictly positive for small positive `a,b` by irreducibility. Interchange
the components for vertices of `B`. Analyticity permits the differentiated
remainder, and there are finitely many chosen vertices and events.
This proves the asserted closure and allows successive weak joining
along any finite tree of components. The permitted bridge probabilities
depend on the components and rates; they cannot be increased arbitrarily
using this argument. The construction changes kernel weights, so it is
not an all-rate theorem for arbitrary simple unweighted graphs.

Two restricted sign statements do follow from existing formulas. For
every positive collection of rates, `U(1_{X=V})>0`. Indeed,
`partial_(K_uu) log det(G_A)=G_uu-G^killed_A_uu>0` for every nonempty
`A`, where the second Green function stops on hitting `A` and is zero
if `u in A`. Strictness follows from irreducibility: a walk can hit `A`
and subsequently return to `u`. Differentiating the full-occupation
normalizer (22) proves the assertion. Further, along fixed proportional
small rates `lambda_i=t alpha_i`, the complete order `Dnu <=st g nu`
holds for all sufficiently small positive `t`, strictly on nonconstant
increasing events. Differentiate the small-rate cylinder expansion;
each minimum trigger has a strictly positive leading coefficient by
the same Green-minor derivative. Joint analyticity in the rates and
kernel entries controls the differentiated remainder, and the finite-event
argument applies.



## Complete core with product recipients

**Product recipients and uniform entrance faces (proved; authoring AI audited).**
The weak-killing result gives a stronger conclusion for this recipient.
Fix a first block K_e=(1-e)P with irreducible P and positive finite rates,
and a finite recipient W with K_2=0 and positive finite rates. For each
sufficiently small fixed e, one b0>0 works simultaneously for all
forward strengths0<a<=1, reverse strengths0<b<b0, forward entrances
rho on W, and reverse entrances eta into block1:

`Cov(N_1,g(X_W))<=0` for every increasing g.            (67p)

For a first block of at most three sites, N_1 can be replaced by **every
increasing rank function h(N_1)**, with one b interval for all h and g.
Consequently (N_1,N_2) is NA, strictly for nonconstant functions of both
counts. The transfer kernel, in row convention, is
`K=[[ (1-e)P, a e 1 rho^T ],[ b 1 eta^T, 0 ]]`.
These entrances specify transfers, not the external addition drive.
Use the canonical saturated invariant law; a full-support external drive
also makes it the unique stationary law. No uniqueness for an arbitrary
restricted external drive is asserted at entrance faces.

Here is the uniformity argument, including the degenerate faces. Write
`Phi_(h,eta)(u)=Cov(h(N(Y)),u^L)/(1-u)`, polynomially extended to u=1.
For n_1<=3, (36m) gives for the finite threshold basis h_k a uniform
lower bound `Phi_(h_k,eta)(u)>=delta_e>0` for all eta and u in[0,1].
One may take `delta_e=(1/2)min_k d_k e^(k+1)` after choosing e small.
For arbitrary n_1 and h=N_1 alone, (36f) supplies such a bound instead.
At u=1 the value is `-Cov(h(N(Y)),L)`; the formula Var(L)/2 is specific
to linear h and is not used for other thresholds.

For increasing g put `c_i(g)=Cov_(nu_2)(X_i,g)>=0` and let E(g) be
its essential coordinates. Since nu_2 is a full-support product law,
c_i(g)>0 exactly for i in E(g). The preceding product eigenvector
calculation, now with f=h in (66), gives
`partial_b Cov(h(N_1),g)|_0`
` =-a sum_i rho_i c_i(g) Phi_(h,eta)(1-a rho_i)`.
If rho(E(g))=0, all sites on which g depends have no incoming edge.
They are upstream singleton SCCs, so SCC factorization gives exact
independence of X_E from the rest for every a,b. The covariance thus
vanishes on this entrance face, as well as on a=0 and b=0.

For a nonconstant upper-set indicator g with E proper in W, write
rho=t sigma on E and rho=(1-t)tau off E. All three parameter domains
are compact. On0<=a<=1 and0<=b<=b_*<1, the walk has a uniform
positive probability of killing within two jumps: from block1 at least
e(1-b_*), and from W at least1-b_*. The finite raw ARW absorption
resolvent is therefore invertible and smooth throughout the parameter
set, including a=1, t=0 and all entrance faces. Exact vanishing on
the three faces factors the covariance as `C=abt F`, with continuous
F, by integrating its mixed derivative along the three intervals.
At b=0, including a t=0 by continuity,
`F=-sum_(i in E)c_i(g)sigma_i Phi_(h_k,eta)(1-a t sigma_i)`
` <=-delta_e sum_(i in E)c_i(g)sigma_i`.
The last sum is at least min_(i in E)c_i(g)>0. Uniform continuity
then gives
`C<=-(delta_e/2)ab sum_i rho_i c_i(g)` for0<b<b_(k,g).
When E=W omit t and factor only ab; the same positive minimum applies.

For n_1<=3, take the minimum over finitely many thresholds h_k and upper-set
indicators g. Nonnegative threshold and superlevel decompositions extend
the estimate to all increasing observables:

`Cov(h(N_1),g(X_W))`
` <=-(delta_e/2)ab [h(n_1)-h(0)]sum_i rho_i c_i(g)`.    (67q)

It is strict when h is nonconstant and g depends on a positive-entrance
site. Otherwise the covariance is exactly zero by the upstream-SCC
argument. For any nonconstant increasing G(N_2), every recipient site
is essential, so strictness holds for every probability entrance rho.
The threshold b0 may depend on e and the fixed blocks/rates, with no
uniformity as e decreases, volume grows or rates approach their endpoints.
This proves cross-block comparisons and two-count NA, not site-level NA:
an arbitrary first-block event need not be a rank function, and an NA
test may distribute recipient coordinates between both observables.

**Full site NA for a complete core and independent recipients (proved;
authoring AI checked and independent audit passed).** Fix a fully
permutation-invariant irreducible core kernel P on n sites, a common
positive finite core sleep rate, and m recipient sites with arbitrary
fixed positive finite rates. Let eta be any common reverse entrance law
on the core. Write
t_i=a rho_i, where t_i>=0 and sum t_i=a<=1, and use

`K_VV=(1-e)P`, `K_Vi=e t_i 1_V`,
`K_iV=b eta^T`, `K_ij=0` for recipient i,j.            (67core)

For every sufficiently small e>0, one b0(e)>0 gives full site NA for
all0<a<=1, all rho, all eta and all0<b<b0(e). Both BK orientations hold after
possibly shrinking b0. Covariance of disjoint nonconstant increasing
tests is strictly negative exactly when both depend essentially on the
interacting SCC V union {i:t_i>0}; a test using only zero-weight
recipients is independent of every disjoint test. The thresholds depend
on the fixed core, recipient rates and size. No uniform b interval as
e decreases to zero or arbitrary-feedback conclusion is asserted.
The base argument below uses uniform eta for configuration tests. The
audited configuration-PGF extension (A19)–(A21) below removes that
restriction, with the same threshold uniform over the full eta simplex.
Different recipient entrance laws and different return probabilities
remain separate working extensions awaiting completion/review.

Let nu_e be the isolated core law, A_e its column addition matrix at
eta, and N the core-count diagonal matrix. The loss-weighted matrix is
`A_e(z)_(Y,X)=z^(1+|X|-|Y|)A_e(Y,X)`; impossible negative exponents
have zero entries. Thus it is a polynomial, including z=0. Proper
input columns are analytic through e=0 by closed absorption with at
most n particles. The saturated law is analytic and nu_0=e_V. Common
invariance reconstructs the remaining column:

`nu_e(V)A_e(:,V)=nu_e-sum_(X!=V)nu_e(X)A_e(:,X)`.       (67ext)

It too is analytic and tends to e_V. Therefore
`A_0(z)e_V=z e_V`: the limiting full core emits exactly one particle.
This is a terminal-kernel continuation, not closed stabilization of
n+1 particles. These facts hold for any fixed irreducible core and
positive rates, without symmetry.

For completeness the whole coupled terminal law has a uniform analytic
extension, even at a=1 and every entrance face. Temporarily give the
recipients separate returns0<=b_i<=b_*<1. First stabilize the saturated
core and each original recipient particle separately, holding routed
exits in queues. Their laws are nu_e and independent Bernoulli(p_i),
where p_i=mu_i/(1+mu_i). Their exit counts are n-|X| and1-X_i.
Thereafter process a queued particle by its block addition kernel and
route each exit independently; use fresh unused stack suffixes. This
is the legal recapture sampler for e>0. Occupied count plus queue size
is at most N_tot=n+m, so the state space is finite.

At e=0 a proper core input absorbs the new particle without loss; a
full input emits one. Process all core queues first, in at most N_tot
steps. At each subsequent recipient step prescribe terminal bit1 and
kill the at most one emitted particle. With
delta=(min_i p_i)(1-b_*)>0, this has probability at least delta and
decreases queue size. Hence uniformly over all parameters,
`Pr(absorption within2N_tot steps)>=delta^N_tot`.
Continuity preserves half this bound for small e. The finite coarse
resolvent is uniformly invertible and analytic. Raw toppling times
need not remain bounded. In the saturated e=0 continuation the core
stays full; eliminating this router gives the actual limiting leaf kernel

`K^leaf_ij=b_i t_j`.                                  (67router)

Here are the first two responses, retaining the changing means. Put
xi_i=X_i-p_i, q_i=1-p_i, v_i=p_i q_i, and let pi_W be the product
recipient law. At zero returns the complete signed first law is

`partial_(b_i)nu|_0`
` =t_i(I-A_e)N nu_e tensor pi_W`
`  +(I-A_e(1-t_i))nu_e tensor (xi_i pi_W)`.             (67first)

Indeed the marked recipient-i sink measure is
`[1-X_i+t_i(n-N)](nu_e tensor pi_W)`. Replace its marked killing by a
core addition. The xi_i bias survives just when none of the L core
losses chooses i, with probability(1-t_i)^L; every other recipient
remains product. This proves (67first) and
`[b_i]E xi_i=v_i[1-G_e(1-t_i)]`, where
`G_e(z)=1^T A_e(z)nu_e=E z^L`.

For distinct i,j put t=t_i,u=t_j,s=1-t-u,x=1-t,y=1-u and H_s=I-A_e(s).
Let M_f be the core-state signed column E[1_(core=X)f]. The exact
two-target harmonic identity is
`M_(xi_i xi_j)=q_i b_i H_s M_(X_i xi_j)`
` +q_j b_j H_s M_(X_j xi_i)`
` -q_i q_j b_i b_j H_s^2 M_(X_i X_j)`.
To verify it, test against a base core harmonic eigenrow phi A_e=h phi.
With both recipient active weights zero, core active and sleeper
weights scale by s; its stable row is phi D_s, D_s=diag(s^|X|).
The recipient factors, after multiplication by -p_i,-p_j, are
`xi_i-q_i b_i(1-sh)X_i` and `xi_j-q_j b_j(1-sh)X_j`.
Moreover `(phi D_s)H_s=(1-sh)phi D_s`. The base harmonic rows form a
joint addition eigenbasis at e>0; no simple spectrum for A_e is needed.
This proves the vector identity for s>0, and continuity supplies s=0.
It assumes no commutation of differently tilted A_e(z).

Substituting (67first) gives
`[b_i b_j]M_(xi_i xi_j)`
` =v_i v_j(I-A_e(s))[I+A_e(s)-A_e(x)-A_e(y)]nu_e`.
Define the ordered, successive-addition PGF
`G_(2,e)(z,w)=E[z^L_1 w^L_2]=1^T A_e(w)A_e(z)nu_e`.
Taking mass and subtracting the first mean responses yields exactly

`[b_i b_j]Cov(X_i,X_j)=-v_i v_j t u H_e(t,u)`,
`H_e(t,u)=[G_e(1-t)G_e(1-u)+G_(2,e)(s,s)`
` -G_(2,e)(1-t,s)-G_(2,e)(1-u,s)]/(t u)`.             (67two)

These losses are successive, not independent. There is no factor two
in the mixed coefficient. Stationarity gives G_e(1)=1 and
G_(2,e)(1,z)=G_e(z), so the polynomial numerator is exactly divisible
by tu. Its quotient extends analytically to both faces. By (67ext),
G_0(z)=z and G_(2,0)(z,w)=zw, so

`H_0(t,u)=1`, `H_e(t,u)>=1/2` for all small e,          (67H)

uniformly on t,u>=0,t+u<=1. This second-return lemma applies to any
fixed irreducible weakly killed core and any fixed entrance, with
arbitrary positive site rates. It does not alone prove core NA.

It remains to handle all disjoint tests, rather than block counts.
For F(core,S) and leaf-only G(T) on disjoint coordinates, let bar F
average F over its own base leaf bits and c_i(G)=Cov_piW(X_i,G).
Centering (67first) gives
`[b]Cov(F,G)=-sum_(i in T)c_i(G)Cov(bar F(Y),(1-t_i)^L)`.
Uniform eta makes the actual joint (X,Y,L) exchangeable in the core,
so bar F can be replaced by its increasing rank average h_F. If F
depends on the core, Delta_F=h_F(n)-h_F(0)>0. The all-rank quotient
bound (36m,n) supplies a common delta_e>0, and hence

`[b]Cov(F,G)<=-delta_e Delta_F sum_i t_i c_i(G)`.       (67mix)

For two leaf-only tests on disjoint supports S,T, the first derivative
is zero. A second-order monomial must contain one return edge from
each support: without such an edge the corresponding leaves are
independent downstream singleton SCCs. With only returns i,j present,
all other leaves average independently and the two functions become
affine in X_i,X_j. Therefore (67two) gives

`[b^2]Cov(F,G)`
` =-sum_(i in S,j in T)t_i t_j H_e(t_i,t_j)c_i(F)c_j(G)`.(67leaf)

Choose e small enough for strict core NA, (67mix) and (67H), then fix
it. If both tests use the core, their base leaf averages are nonconstant
on disjoint core supports. Their covariance at b=0 is strictly negative,
independent of t, and stays so uniformly for small b. If exactly one
uses the core, factor its covariance as b r J, where r is the total
entrance mass of the leaf-only test's essential coordinates. This
factorization is smooth: covariance vanishes at b=0 and at r=0, where
those recipients are independent upstream singleton SCCs. Equation
(67mix) makes J uniformly negative at b=0. For two leaf-only tests,
the zero first derivative and both upstream-SCC zeros similarly give
b^2 r s J; (67leaf,H) makes this quotient uniformly negative at b=0.
Group entrance masses can be decreased along every integration rectangle
without leaving the simplex, justifying division on its faces.

Positive minimum product influences on essential coordinates and compact
continuity give one b interval for each pair of indicators. There are
finitely many pairs; their common minimum and nonnegative superlevel
decomposition give full NA and the stated strictness. Finally the coarse
resolvent gives uniform convergence to nu_e tensor pi_W as b tends to
zero. The core is strictly BK in both orientations by (36wc)–(36we).
The independently audited neighborhood lemma (41bk) therefore proves
both BK orientations after a further uniform decrease of b0(e).

Stationarity, upward skip-freeness, and even reversibility with nonnegative
spectrum do not imply (67). The abstract column kernel
`T=[[3/4,1/4,0],[1/4,37/50,1/800],[0,1/100,799/800]]`
has invariant law `(1/10,1/10,4/5)`, is reversible and lazy, and its two
nontrivial linear rank responses are 1/40 and 1/1000. Nevertheless
`Psi(u)=(1-u)(108u-43)/2500`, negative at u=1/4. This is not an ARW
kernel counterexample. The unresolved sign in (67) needs more of the
ARW structure than these general count-chain properties.



### Closed additions and configuration loss-PGF

**Proved; independent audit completed 2026-11-01, authoring AI reviewed.** The
closed-addition lemma, uniform configuration-PGF theorem, and arbitrary
common reverse-entrance gluing consequence have passed analytic audits.
The gluing consequence includes both BK orientations by (41bk). Neither
separate recipient entrance laws nor separate return probabilities are
included in this audit.

#### 1. Closed addition identities, including the occupied-tag step

First take the closed holding kernel P=J/n, n>=2, and common sleep rate
lambda=p/q, where 0<p<1 and q=1-p. Let U_m be the uniform law on stable
m-subsets of V. Fix v. Write E_m^0 and E_m^1 for U_m conditioned on v
being absent or present, respectively. Their domains are

    E_m^0: 0<=m<=n-1;       E_m^1: 1<=m<=n.

Let A_v denote a point addition followed by closed stabilization, and
Abar=n^(-1)sum_v A_v. These kernels are used only with total particle
number at most n. Closed absorption is finite: route surplus particles
directly to empty sites until there is no surplus, then sleep all active
singletons. This gives a positive finite legal prefix from each state.
The finite transient system therefore has finite mean absorption time.
The abelian property and fresh unused stack suffixes give commutation of
two addition kernels whenever the initial count is at most n-2.

Conservation and permutation symmetry give

    Abar U_m=U_(m+1),                    0<=m<=n-1.       (A1)

The two first-step identities are

    A_v E_m^0=p E_(m+1)^1+q Abar E_m^0, 0<=m<=n-1,
    A_v E_m^1=A_v Abar E_(m-1)^0,        1<=m<=n-1.      (A2)

For the second identity, adding to the sleeper at v creates two active
particles there. Discard ineffective sleep instructions until the first
jump. Its destination is uniform on V, including v, and the unused stack
suffixes are fresh and independent. After this jump the state consists
of the original m-1 outside sleepers, one active particle at v, and one
active particle at that uniform destination. Stabilizing these two active
additions is exactly A_v Abar E_(m-1)^0. This is a forced first jump from
an overfull site; it has no extra factor q. Its total particle count is
m+1<=n, so commutation invokes no closed n+1-particle stabilization.

Induction now proves

    A_v E_m^0=A_v E_m^1=A_v U_m,          1<=m<=n-1.      (A3)

The induction base is only A_v E_0^0=A_v U_0; E_0^1 is undefined and is
not asserted. Assuming that equality at m-1, (A2), commutation and (A1)
give

    A_v E_m^1=Abar A_v E_(m-1)^0
              =Abar A_v U_(m-1)=A_v U_m.

Since U_m=(m/n)E_m^1+((n-m)/n)E_m^0 and n-m>0, the same equality follows
for E_m^0. At m=n-1 the common output is the deterministic full set.

#### 2. The tag bias and its positive recurrence

For 1<=m<=n-1 put T_m=E_m^1-E_m^0. The signed measures invariant under
permutations fixing v, of mass zero on a fixed nontrivial rank sector,
form a one-dimensional space. Thus

    Abar T_m=r_m T_(m+1),                1<=m<=n-2,
    A_v U_m=U_(m+1)+c_m T_(m+1),        0<=m<=n-2.      (A4)

The coefficient c_m is precisely the excess terminal occupancy
probability at v over (m+1)/n. At m=0, (A2) gives

    c_0=p(n-1)/n>0.                                    (A5)

For 1<=m<=n-2, commutation applied to U_(m-1) gives
c_m=c_(m-1) r_m. Also

    E_m^0=U_m-(m/n)T_m,
    E_(m+1)^1=U_(m+1)+((n-m-1)/n)T_(m+1).

Substitution in the first identity (A2), using (A3), gives

    n c_m=p(n-m-1)-q m r_m.

Consequently

    r_m=p(n-m-1)/(n c_(m-1)+q m)>0,
    c_m=p(n-m-1)c_(m-1)/(n c_(m-1)+q m)>0.             (A6)

This proves the claimed bias at every nonfull output rank. The terminal
case is Abar T_(n-1)=0 and A_v U_(n-1)=U_n; no T_n or c_(n-1) is
needed. In particular there is no strict tag excess at the full output.

#### 3. Exact no-loss column sums

Fix a terminal hole set I of size 1<=r<=n-1, and put m=n-r-1. Define

    S_v(I)=sum_(|J|=r+1) A_v(V\I,V\J)
          =binom(n,r+1) Pr_(U_m,A_v)(Y=V\I).           (A7)

The residual symmetry fixing v makes the output uniform within each of
the two classes according to the presence of v. Equation (A4) therefore
gives the exact formulas

    v in I:
      S_v(I)=(n-r)/(r+1) [1-n c_m/r],
    v not in I:
      S_v(I)=(n-r+n c_m)/(r+1).                        (A8)

For the first formula use binom(n-1,r-1) possible terminal hole sets and
total class probability r/n-c_m. For the second use binom(n-1,r) sets
and class probability (n-r)/n+c_m. These also verify the factorial and
uniform-input normalizations in (A7).

By (A6), c_m>=0, and by nonnegativity of the first class, c_m<=r/n.
Hence

    S_v(I)<=(n-r)/(r+1) if v in I,
    S_v(I)<=n/(r+1)     if v not in I.                 (A9)

In particular S_eta(I)<=n/(r+1) for every entrance distribution eta,
where S_eta=sum_v eta_v S_v. No exchangeability of the eta-addition law
is assumed. At r=n there is no possible (n+1)-hole input, and the
corresponding no-loss sum is zero.

These bounds also hold for every fully permutation-invariant irreducible
closed kernel with common sleep rate. Indeed, write P_ii=d and P_ij=w>0
for i!=j, so d+(n-1)w=1. Erasing null self-jumps identifies its stable
closed addition kernels with those of J/n at sleep rate

    lambda_H=lambda/(n w).                             (A10)

Both models have the same singleton sleep-to-nonself-jump ratio
lambda/((n-1)w), and the same uniform nonself jump destination. This
argument is for closed addition kernels; it does not freeze a transformed
sleep rate in a varying killed family.

#### 4. Weak-killing configuration-event no-loss gap

Fix one such homogeneous closed P, its common positive finite lambda,
and K_e=(1-e)P. Let nu_e be the canonical stationary law. For a prescribed
k-set I, write its exact-hole leading coefficient as gamma_k:

    nu_e(V\I)=e^k gamma_k+O(e^(k+1)),  gamma_k>0,
    M=n gamma_1.

We use the already-proved complete-graph OR-cylinder consequence

    gamma_(r+1)<=gamma_1 gamma_r,             r+1<=n.   (A11)

Only this non-strict leading-cylinder bound is needed here, not the
later strict-cylinder theorem.

Let f be a nonconstant increasing indicator on the full configuration
cube. Thus f(V)=1 and f(empty)=0. Put

    F(I)=1-f(V\I),
    r=min{|I|: F(I)=1},
    A_r=#{I: |I|=r, F(I)=1},       C_f=gamma_r A_r>0.

The failure event F is increasing in the holes. In a stationary
eta-addition X->Y, put L=1+|X|-|Y| and
Psi_(f,e,eta)(u)=Cov(f(Y),u^L), with 0^0=1. Then

    Pr(f(Y)=0)=e^r C_f+O(e^(r+1)),
    Pr(L=0)=e M+O(e^2).

On no loss, the initial hole count is one greater than the final count.
Only the (r+1)-hole input sector can contribute at order e^(r+1) to
Pr(f(Y)=0,L=0). Its leading input coefficient is gamma_(r+1) per set,
and its e=0 transition is the closed addition from sections 1-3. Thus,
for r<n,

    Pr(f(Y)=0,L=0)
      =e^(r+1) gamma_(r+1)
           sum_(|I|=r,F(I)=1) S_eta(I)+O(e^(r+2)).      (A12)

For r=n the probability in (A12) is identically zero. Consequently

    Psi_(f,e,eta)(0)=e^(r+1)d_(f,eta)+O(e^(r+2)),
    d_(f,eta)=M C_f-gamma_(r+1) sum_I S_eta(I),         (A13)

with the second term zero at r=n. Equations (A9) and (A11) prove

    d_(f,eta)>=(r/(r+1)) M C_f=:d_f^*>0                (A14)

uniformly in eta; at r=n the exact value is M C_f. This is the desired
configuration-level no-loss gap. It is not obtained by symmetrizing f.

#### 5. Analyticity, all u, and simultaneous configuration tests

Every nonfull-input addition has at most n particles. Its finite closed
absorption system is invertible at e=0, giving analytic transition
probabilities. For a transition from h>=1 initial holes to j terminal
holes, conservation requires exactly ell=1+j-h kills. When ell>=0 its
probability is divisible by e^ell: grade transient states by particle
number, and each downward step of the block resolvent has one factor e.
For ell<0 the transition is impossible.

The stationary full-input joint column is reconstructed, separately for
each pure entrance, by

    nu_e(V) A_(e,v)(Y,V)
       =nu_e(Y)-sum_(X!=V)nu_e(X)A_(e,v)(Y,X).          (A15)

Thus all stationary joint masses are analytic without assuming closed
n+1-particle absorption. There are finitely many pure entrances, and
mixtures are linear in eta, so all remainder bounds below are uniform
over the entire entrance simplex, including entrances depending on e.

For each terminal r-hole configuration, every nonfull-input joint mass
is O(e^(r+1)): its input factor e^h is multiplied by e^(1+r-h).
The leading e^r gamma_r mass therefore comes from full input, with
L=r+1. Also E[u^L]=u+O(e) coefficientwise, since the limiting stationary
transition is full to full and has loss 1. It follows that

    Psi_(f,e,eta)(u)
       =e^r C_f [u-u^(r+1)]+O_coeff(e^(r+1)).          (A16)

Since 0<=L<=n+1, these are polynomials of uniformly bounded degree.
Divide by 1-u polynomially, using Psi(1)=0, and combine (A13),(A16):

    Psi(u)/(1-u)
      =e^r C_f u(1+u+...+u^(r-1))
         +e^(r+1)[d_(f,eta)+u a_(f,eta)(u)]
         +e^(r+2)b_(e,f,eta)(u).                      (A17)

The a,b polynomials are uniformly bounded over eta and u in [0,1] for
small e. Since d_(f,eta)>=d_f^*>0, a sufficiently small e interval gives

    Psi_(f,e,eta)(u)/(1-u)
       >=(C_f/2)e^r u+(d_f^*/2)e^(r+1)>0,             (A18)

including the continuous quotient at u=1. This proves strict covariance
for every 0<=u<1, including u=0.

There are finitely many nonconstant increasing indicator functions on
the fixed cube. Taking the minimum of their e thresholds and applying
nonnegative superlevel decomposition gives the simultaneous theorem:

For every fixed homogeneous irreducible complete P and common positive
finite sleep rate, there is e_0>0 such that for 0<e<e_0,

    Cov(f(Y),u^L)>0                                    (A19)

for every nonconstant increasing configuration function f, every entrance
eta, and every 0<=u<1. Moreover for each such fixed e there is delta_e>0,
independent of eta,u,f, such that

    Cov(f(Y),u^L)/(1-u)
       >=delta_e [f(V)-f(empty)].                     (A20)

For n=1 the closed addition to the empty state ends full. The r=n case
above gives d=M gamma_1=M^2>0, so the same theorem holds directly.
No conditional stochastic ordering of L given individual configurations
is asserted.

#### 6. Consequence for complete core plus product recipients

The only use of uniform reverse core entrance in the argument around
(67mix) above was converting the configuration function bar F into its
rank average. Replace that step
by (A20) directly. The exact first response remains

    [b]Cov(F,G)
      =-sum_i Cov_(pi_W)(X_i,G)
          Cov_coreadd(bar F(Y),(1-t_i)^L).

It is bounded above by

    -delta_e [bar F(V)-bar F(empty)]
          sum_i t_i Cov_(pi_W)(X_i,G),                 (A21)

uniformly over all common reverse entrance laws eta, including point
entrances and simplex faces.

The other uniformity requirements remain valid with eta variable:

* A_(e,eta)(z) is linear in eta, and its full-state router limit is
  A_(0,eta)(z)e_V=z e_V for every eta. Thus the normalized two-leaf
  coefficient H_(e,eta)(t,u) is continuous on the compact joint parameter
  set, equals 1 at e=0, and has a common positive lower bound for small e.
* The finite coarse-chain absorption bound before (67router) is independent of eta:
  a proper core absorbs an arrival at e=0, while a full core emits exactly
  one particle; the prescribed leaf-success/kill sequence creates no
  returns. Hence its absorption and analytic bounds are uniform also
  over eta, with the same a=1 and forward-face coverage.
* Since the closed core is irreducible, arbitrary nonempty reverse
  entrance support leaves the interacting SCC unchanged. Exact zeros
  on forward entrance faces, and the b t_E or b^2 t_E t_F
  factorizations, are unchanged. Compactness can now include eta.

It follows that the full site-NA gluing theorem holds for arbitrary
common reverse core entrance eta, with one small-e threshold and, for
each fixed allowed e, one b_0(e)>0 uniform in a,rho,eta and all disjoint
increasing tests. This does not claim a b threshold uniform as e->0 or
allow different reverse entrance laws at different recipient sites.
Both BK orientations also follow after shrinking b_0(e): apply (41bk)
to the strict-BK isolated complete core and its product recipient, using
the same uniform convergence over t and eta. This uses the separate
BK-neighborhood lemma, not a consequence of the PGF inequality alone.

#### Boundary and normalization audit

The authoring AI's formulas need only explicit domain conventions, not a
mathematical repair: E_0^1 is absent; T_n is absent; m=n-1 gives a full
deterministic output; r=n has zero no-loss failure term. The occupied-tag
step has no q factor, the column sum has exactly binom(n,r+1), and its
two conditional terminal class sizes are binom(n-1,r-1) and binom(n-1,r).
The general homogeneous-kernel extension uses lambda/(n w) only for
the closed leading addition. All e expansions use the actual physical
uniform-killing parameter and the fixed original sleep rate.

